// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;

using BenchmarkPropagator = Acts::Propagator<Acts::EigenStepper<>>;
using Linearizer = Acts::HelicalTrackLinearizer;

struct AthenaConfiguration {
  std::string era;
  double tracksMaxZinterval;
  double finderMinWeight;
  int finderMaxIterations;
};

AthenaConfiguration athenaConfiguration(const std::string& era) {
  if (era == "run3") {
    return {era, 3_mm, 0.0001, 100};
  }
  if (era == "run4") {
    return {era, 0.5_mm, 0.02, 200};
  }
  throw std::invalid_argument("athena-era must be 'run3' or 'run4'");
}

struct InputData {
  std::vector<BoundTrackParameters> tracks;
  std::vector<InputTrack> inputTracks;
};

struct EventMetadata {
  Vector3 beamPosition = Vector3::Zero();
  SquareMatrix3 beamCovariance = SquareMatrix3::Zero();
  Transform3 perigee = Transform3::Identity();
  bool useConstraint = true;
  std::size_t selectedTracks = 0;
};

EventMetadata readMetadata(const std::string& path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("Cannot read metadata: " + path);
  EventMetadata metadata;
  std::map<std::string, std::vector<double>> fields;
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    std::string key;
    row >> key;
    double value;
    while (row >> value)
      fields[key].push_back(value);
  }
  auto values = [&](const std::string& key, std::size_t size) -> const auto& {
    if (fields[key].size() != size)
      throw std::runtime_error("Invalid metadata: " + key);
    return fields[key];
  };
  if (values("schema", 1)[0] != 1.)
    throw std::runtime_error("Unknown metadata schema");
  metadata.selectedTracks =
      static_cast<std::size_t>(values("selected_tracks", 1)[0]);
  metadata.useConstraint = values("use_beam_constraint", 1)[0] != 0.;
  for (int i = 0; i < 3; ++i) {
    metadata.beamPosition[i] = values("beam_position", 3)[i];
    for (int j = 0; j < 3; ++j)
      metadata.beamCovariance(i, j) = values("beam_covariance", 9)[3 * i + j];
  }
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      metadata.perigee.matrix()(i, j) =
          values("perigee_transform", 16)[4 * i + j];
  }
  return metadata;
}

InputData readTracks(const std::string& path, const Transform3& reference) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("Could not open input track CSV: " + path);
  InputData data;
  auto surface = Surface::makeShared<PerigeeSurface>(reference);
  std::string line;
  std::getline(input, line);
  if (!line.empty() && line.back() == '\r')
    line.pop_back();
  std::istringstream header(line);
  std::map<std::string, std::size_t> columns;
  std::string cell;
  while (std::getline(header, cell, ','))
    columns.emplace(cell, columns.size());
  const bool fullCovariance = columns.contains("trackId");
  constexpr std::array names{"d0", "z0", "phi", "theta", "qop"};
  while (std::getline(input, line)) {
    if (line.empty())
      continue;
    std::istringstream row(line);
    std::vector<double> values;
    while (std::getline(row, cell, ','))
      values.push_back(std::stod(cell));
    if (values.size() != columns.size())
      throw std::runtime_error("Malformed track row in " + path);
    auto value = [&](const std::string& name) {
      return values.at(columns.at(name));
    };
    BoundVector parameters = BoundVector::Zero();
    BoundMatrix covariance = BoundMatrix::Zero();
    if (fullCovariance) {
      // Read the existing ActsExamples TrackParameterData schema without
      // bringing the examples sequencer into the core-only timing harness.
      for (int i = 0; i < 5; ++i) {
        parameters[i] = value(names[i]);
        for (int j = 0; j < 5; ++j) {
          const auto key = i == j ? std::string("var_") + names[i]
                                  : std::string("cov_") + names[i] + names[j];
          covariance(i, j) = value(key);
        }
      }
      covariance(eBoundTime, eBoundTime) = 1.;
    } else {
      constexpr std::array legacyNames{"loc0",  "loc1", "phi",
                                       "theta", "qop",  "time"};
      for (int i = 0; i < 6; ++i) {
        parameters[i] = value(legacyNames[i]);
        const double sigma = value(std::string("sigma_") + legacyNames[i]);
        covariance(i, i) = sigma * sigma;
      }
    }
    if (!parameters.allFinite() || !covariance.allFinite())
      throw std::runtime_error("Nonfinite track input");
    data.tracks.emplace_back(surface, parameters, covariance,
                             ParticleHypothesis::pion());
  }
  data.inputTracks.reserve(data.tracks.size());
  for (const auto& track : data.tracks)
    data.inputTracks.emplace_back(&track);
  return data;
}

AdaptiveMultiVertexFinder makeFinder(
    const std::shared_ptr<ConstantBField>& field, Linearizer& linearizer,
    const ImpactPointEstimator& ipEstimator,
    const AthenaConfiguration& athena) {
  AnnealingUtility::Config annealingConfig{9., {1.}};

  AdaptiveMultiVertexFitter::Config fitterConfig(ipEstimator);
  fitterConfig.annealingTool = AnnealingUtility(annealingConfig);
  fitterConfig.maxIterations = 30;
  fitterConfig.maxDistToLinPoint = 0.5;
  fitterConfig.minWeight = 0.001;
  fitterConfig.maxRelativeShift = 0.01;
  fitterConfig.doSmoothing = true;
  fitterConfig.useTime = false;
  fitterConfig.extractParameters.connect<&InputTrack::extractParameters>();
  fitterConfig.trackLinearizer.connect<&Linearizer::linearizeTrack>(
      &linearizer);

  GaussianTrackDensity::Config densityConfig;
  densityConfig.d0MaxSignificance = 3.5;
  densityConfig.z0MaxSignificance = 12.;
  densityConfig.extractParameters.connect<&InputTrack::extractParameters>();
  auto seedFinder = std::make_shared<TrackDensityVertexFinder>(
      TrackDensityVertexFinder::Config{
          GaussianTrackDensity(std::move(densityConfig))});

  AdaptiveMultiVertexFinder::Config finderConfig(
      AdaptiveMultiVertexFitter(std::move(fitterConfig)), std::move(seedFinder),
      ipEstimator, field);
  finderConfig.initialVariances = Vector4::Constant(1e8);
  finderConfig.tracksMaxZinterval = athena.tracksMaxZinterval;
  finderConfig.tracksMaxSignificance = 5.;
  finderConfig.maxVertexChi2 = 18.42;
  finderConfig.doRealMultiVertex = true;
  finderConfig.useFastCompatibility = true;
  finderConfig.minWeight = athena.finderMinWeight;
  finderConfig.maxIterations = athena.finderMaxIterations;
  finderConfig.addSingleTrackVertices = false;
  finderConfig.doFullSplitting = false;
  finderConfig.maxMergeVertexSignificance = 3.;
  finderConfig.maximumVertexContamination = 0.5;
  finderConfig.useVertexCovForIPEstimation = false;
  finderConfig.useSeedConstraint = false;
  finderConfig.useTime = false;
  finderConfig.doNotBreakWhileSeeding = false;
  finderConfig.extractParameters.connect<&InputTrack::extractParameters>();
  return AdaptiveMultiVertexFinder(std::move(finderConfig));
}

}  // namespace

int runBenchmark(int argc, char* argv[]) {
  std::string inputPath;
  std::string metadataPath;
  std::string inputDumpPath;
  bool seederOnly = false;
  std::string outputPath;
  std::string athenaEra = "run4";
  std::size_t warmup = 1;
  std::size_t repetitions = 7;
  double beamSpotSigmaX = 0.15_mm;
  double beamSpotSigmaY = 0.15_mm;
  double beamSpotSigmaZ = 53_mm;

  po::options_description options("Allowed options");
  options.add_options()("help,h", "show this help")(
      "input", po::value<std::string>(&inputPath)->required(),
      "CSV with six bound parameters followed by their six uncertainties")(
      "metadata", po::value<std::string>(&metadataPath),
      "Athena event metadata")(
      "output-input", po::value<std::string>(&inputDumpPath),
      "Write parsed parameters and covariance for round-trip verification")(
      "seeder-only", po::bool_switch(&seederOnly),
      "Time only Gaussian seeding, independent of field and geometry")(
      "output-vertices", po::value<std::string>(&outputPath),
      "Write the last result at full precision for regression comparisons")(
      "athena-era",
      po::value<std::string>(&athenaEra)->default_value(athenaEra),
      "Athena offline AMVF configuration: run3 or run4")(
      "beamspot-sigma-x",
      po::value<double>(&beamSpotSigmaX)->default_value(beamSpotSigmaX),
      "beam-spot x width in mm")(
      "beamspot-sigma-y",
      po::value<double>(&beamSpotSigmaY)->default_value(beamSpotSigmaY),
      "beam-spot y width in mm")(
      "beamspot-sigma-z",
      po::value<double>(&beamSpotSigmaZ)->default_value(beamSpotSigmaZ),
      "beam-spot z width in mm")(
      "warmup", po::value<std::size_t>(&warmup)->default_value(warmup),
      "warm-up iterations")(
      "repetitions",
      po::value<std::size_t>(&repetitions)->default_value(repetitions),
      "measured iterations");

  try {
    po::variables_map variables;
    po::store(po::parse_command_line(argc, argv, options), variables);
    if (variables.contains("help")) {
      std::cout << options << '\n';
      return 0;
    }
    po::notify(variables);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n\n" << options << '\n';
    return 1;
  }
  if (repetitions == 0) {
    std::cerr << "error: repetitions must be positive\n";
    return 1;
  }

  AthenaConfiguration athena;
  try {
    athena = athenaConfiguration(athenaEra);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
  if (beamSpotSigmaX <= 0. || beamSpotSigmaY <= 0. || beamSpotSigmaZ <= 0.) {
    std::cerr << "error: beam-spot widths must be positive\n";
    return 1;
  }

  EventMetadata metadata;
  if (!metadataPath.empty())
    metadata = readMetadata(metadataPath);
  const InputData data = readTracks(inputPath, metadata.perigee);
  if (!metadataPath.empty() && metadata.selectedTracks != data.tracks.size()) {
    throw std::runtime_error("Track count does not match metadata");
  }
  if (!inputDumpPath.empty()) {
    std::ofstream dump(inputDumpPath);
    if (!dump)
      throw std::runtime_error("Cannot open input dump");
    dump << std::setprecision(17);
    for (const auto& track : data.tracks) {
      for (int i = 0; i < 6; ++i)
        dump << track.parameters()[i] << ' ';
      for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j)
          dump << (*track.covariance())(i, j) << ' ';
      }
      dump << '\n';
    }
  }
  const GeometryContext geometryContext =
      GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext magneticFieldContext;
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  auto propagator =
      std::make_shared<BenchmarkPropagator>(EigenStepper<>(field));
  ImpactPointEstimator::Config ipEstimatorConfig(field, propagator);
  ipEstimatorConfig.maxIterations = 20;
  ipEstimatorConfig.precision = 1e-10;
  ImpactPointEstimator ipEstimator(ipEstimatorConfig);
  Linearizer::Config linearizerConfig;
  linearizerConfig.bField = field;
  linearizerConfig.propagator = propagator;
  Linearizer linearizer(linearizerConfig);
  auto finder = makeFinder(field, linearizer, ipEstimator, athena);
  const Vector3 beamSpotPosition =
      metadataPath.empty() ? Vector3::Zero().eval() : metadata.beamPosition;
  Vertex beamSpotConstraint(beamSpotPosition);
  SquareMatrix3 beamSpotCovariance = SquareMatrix3::Zero();
  beamSpotCovariance.diagonal() << beamSpotSigmaX * beamSpotSigmaX,
      beamSpotSigmaY * beamSpotSigmaY, beamSpotSigmaZ * beamSpotSigmaZ;
  if (!metadataPath.empty())
    beamSpotCovariance = metadata.beamCovariance;
  if (!metadata.useConstraint) {
    beamSpotConstraint.setPosition(Vector3::Zero());
    beamSpotCovariance = SquareMatrix3::Identity() * 1e8;
  }
  beamSpotConstraint.setCovariance(beamSpotCovariance);
  beamSpotSigmaX = std::sqrt(beamSpotCovariance(0, 0));
  beamSpotSigmaY = std::sqrt(beamSpotCovariance(1, 1));
  beamSpotSigmaZ = std::sqrt(beamSpotCovariance(2, 2));
  const VertexingOptions vertexingOptions(geometryContext, magneticFieldContext,
                                          beamSpotConstraint,
                                          metadata.useConstraint);

  GaussianTrackDensity::Config densityConfig;
  densityConfig.extractParameters.connect<&InputTrack::extractParameters>();
  GaussianTrackDensity density(densityConfig);
  std::optional<std::pair<double, double>> seed;
  std::size_t seedingTracks = 0;
  std::vector<double> samplesMs;
  samplesMs.reserve(repetitions);
  std::size_t reconstructedVertices = 0;
  for (std::size_t iteration = 0; iteration < warmup + repetitions;
       ++iteration) {
    if (seederOnly) {
      GaussianTrackDensity::State densityState(data.tracks.size());
      const auto start = std::chrono::steady_clock::now();
      auto result =
          density.globalMaximumWithWidth(densityState, data.inputTracks);
      const auto stop = std::chrono::steady_clock::now();
      if (!result.ok())
        throw std::runtime_error(result.error().message());
      seed = *result;
      seedingTracks = densityState.trackEntries.size();
      if (iteration >= warmup)
        samplesMs.push_back(
            std::chrono::duration<double, std::milli>(stop - start).count());
      continue;
    }
    auto state = finder.makeState(magneticFieldContext);
    const auto start = std::chrono::steady_clock::now();
    auto result = finder.find(data.inputTracks, vertexingOptions, state);
    const auto stop = std::chrono::steady_clock::now();
    if (!result.ok()) {
      std::cerr << "AMVF failed: " << result.error().message() << '\n';
      return 2;
    }
    reconstructedVertices = result->size();
    // Serialize outside the timed region, including association identities and
    // fitted parameters rather than relying only on the vertex count.
    if (!outputPath.empty() && iteration + 1 == warmup + repetitions) {
      std::ofstream output(outputPath);
      if (!output) {
        throw std::runtime_error("Could not open vertex output: " + outputPath);
      }
      output << std::setprecision(17);
      for (const auto& vertex : *result) {
        output << "vertex\n"
               << vertex.fullPosition() << '\n'
               << vertex.fullCovariance() << '\n'
               << vertex.fitQuality().first << ' ' << vertex.fitQuality().second
               << '\n';
        for (const auto& track : vertex.tracks()) {
          const auto* original =
              track.originalParams.as<BoundTrackParameters>();
          output << "track " << (original - data.tracks.data()) << ' '
                 << track.trackWeight << ' ' << track.chi2Track << ' '
                 << track.ndf << ' ' << track.vertexCompatibility << '\n'
                 << track.fittedParams.parameters() << '\n';
          if (track.fittedParams.covariance()) {
            output << *track.fittedParams.covariance() << '\n';
          }
        }
      }
    }
    if (iteration >= warmup) {
      samplesMs.push_back(
          std::chrono::duration<double, std::milli>(stop - start).count());
    }
  }

  std::ranges::sort(samplesMs);
  const double median = samplesMs.size() % 2 == 0
                            ? (samplesMs[samplesMs.size() / 2 - 1] +
                               samplesMs[samplesMs.size() / 2]) /
                                  2.
                            : samplesMs[samplesMs.size() / 2];
  const double mean = std::accumulate(samplesMs.begin(), samplesMs.end(), 0.) /
                      static_cast<double>(samplesMs.size());

  std::cout << std::setprecision(9) << "{\n"
            << "  \"input\": " << std::quoted(inputPath) << ",\n"
            << "  \"athena_era\": " << std::quoted(athena.era) << ",\n"
            << "  \"use_beam_constraint\": "
            << (metadata.useConstraint ? "true" : "false") << ",\n"
            << "  \"tracks_max_z_interval_mm\": "
            << athena.tracksMaxZinterval / 1_mm << ",\n"
            << "  \"finder_min_weight\": " << athena.finderMinWeight << ",\n"
            << "  \"finder_max_iterations\": " << athena.finderMaxIterations
            << ",\n"
            << "  \"beamspot_sigma_mm\": [" << beamSpotSigmaX << ", "
            << beamSpotSigmaY << ", " << beamSpotSigmaZ << "],\n"
            << "  \"tracks\": " << data.tracks.size() << ",\n"
            << "  \"reconstructed_vertices\": " << reconstructedVertices
            << ",\n"
            << "  \"warmup_runs\": " << warmup << ",\n"
            << "  \"measured_runs\": " << repetitions << ",\n"
            << "  \"median_ms\": " << median << ",\n"
            << "  \"mean_ms\": " << mean << ",\n"
            << "  \"min_ms\": " << samplesMs.front() << ",\n"
            << "  \"max_ms\": " << samplesMs.back() << ",\n"
            << "  \"samples_ms\": [";
  for (std::size_t i = 0; i < samplesMs.size(); ++i) {
    std::cout << (i == 0 ? "" : ", ") << samplesMs[i];
  }
  std::cout << "],\n  \"seeder_only\": " << (seederOnly ? "true" : "false")
            << ",\n  \"seed_z_width\": ";
  if (seed)
    std::cout << std::setprecision(17) << '[' << seed->first << ','
              << seed->second << ']';
  else
    std::cout << "null";
  std::cout << ",\n  \"metadata\": " << std::quoted(metadataPath)
            << ",\n  \"field_tesla\": 2,\n  \"seeding_tracks\": "
            << seedingTracks << ",\n  \"beamspot_position_mm\": [";
  for (int i = 0; i < 3; ++i)
    std::cout << (i ? "," : "") << beamSpotConstraint.position()[i];
  std::cout << "],\n  \"beamspot_covariance_mm2\": [";
  for (int i = 0; i < 9; ++i)
    std::cout << (i ? "," : "") << beamSpotCovariance(i / 3, i % 3);
  std::cout << "]\n}\n";
  return 0;
}

int main(int argc, char* argv[]) {
  try {
    return runBenchmark(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
