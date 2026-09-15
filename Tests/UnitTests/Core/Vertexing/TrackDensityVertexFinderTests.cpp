// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <system_error>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

namespace ActsTests {

using Covariance = BoundMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext::dangerouslyDefaultConstruct();
MagneticFieldContext magFieldContext = MagneticFieldContext();

BOOST_AUTO_TEST_SUITE(VertexingSuite)
///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation, i.e.
/// tests if reordering tracks returns the same result
///
BOOST_AUTO_TEST_CASE(track_density_finder_test) {
  // Define some track parameter properties
  Vector3 pos0{0, 0, 0};
  Vector3 pos1a{1.86052_mm, -1.24035_mm, -10_mm};
  Vector3 mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3 pos1b{-1.24035_mm, 1.86052_mm, -3_mm};
  Vector3 mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3 pos1c{1.69457_mm, -0.50837_mm, -7_mm};
  Vector3 mom1c{300_MeV, 1000_MeV, 100_MeV};

  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect<&InputTrack::extractParameters>();
  TrackDensityVertexFinder finder{
      TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)}};
  auto state = finder.makeState(magFieldContext);

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  auto params1a =
      BoundTrackParameters::create(
          geoContext, perigeeSurface, makeVector4(pos1a, 0), mom1a.normalized(),
          1_e / mom1a.norm(), covMat, ParticleHypothesis::pion())
          .value();
  auto params1b =
      BoundTrackParameters::create(
          geoContext, perigeeSurface, makeVector4(pos1b, 0), mom1b.normalized(),
          -1_e / mom1b.norm(), covMat, ParticleHypothesis::pion())
          .value();
  auto params1c =
      BoundTrackParameters::create(
          geoContext, perigeeSurface, makeVector4(pos1c, 0), mom1c.normalized(),
          1_e / mom1c.norm(), covMat, ParticleHypothesis::pion())
          .value();

  // Vectors of track parameters in different orders
  std::vector<InputTrack> vec1 = {InputTrack{&params1a}, InputTrack{&params1b},
                                  InputTrack{&params1c}};
  std::vector<InputTrack> vec2 = {InputTrack{&params1c}, InputTrack{&params1a},
                                  InputTrack{&params1b}};

  auto res1 = finder.find(vec1, vertexingOptions, state);
  auto res2 = finder.find(vec2, vertexingOptions, state);

  if (!res1.ok()) {
    std::cout << res1.error().message() << std::endl;
  }

  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }

  if (res1.ok() && res2.ok()) {
    BOOST_CHECK(!(*res1).empty());
    BOOST_CHECK(!(*res2).empty());
    Vector3 result1 = (*res1).back().position();
    Vector3 result2 = (*res2).back().position();
    BOOST_CHECK_EQUAL(result1, result2);
  }
}

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation
///
BOOST_AUTO_TEST_CASE(track_density_finder_constr_test) {
  // Define some track parameter properties
  Vector3 pos0{0, 0, 0};
  Vector3 pos1a{1.86052_mm, -1.24035_mm, -10_mm};
  Vector3 mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3 pos1b{-1.24035_mm, 1.86052_mm, -3_mm};
  Vector3 mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3 pos1c{1.69457_mm, -0.50837_mm, -7_mm};
  Vector3 mom1c{300_MeV, 1000_MeV, 100_MeV};

  // From Athena VertexSeedFinderTestAlg
  double const expectedZResult = -13.013;

  // Create constraint for seed finding
  Vector3 constraintPos{1.7_mm, 1.3_mm, -6_mm};
  SquareMatrix3 constrCov = SquareMatrix<3>::Identity();

  Vertex constraint(constraintPos);
  constraint.setCovariance(constrCov);

  // Finder options
  VertexingOptions vertexingOptions(geoContext, magFieldContext, constraint);
  GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect<&InputTrack::extractParameters>();
  TrackDensityVertexFinder finder{
      TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)}};
  auto state = finder.makeState(magFieldContext);

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  auto params1a =
      BoundTrackParameters::create(
          geoContext, perigeeSurface, makeVector4(pos1a, 0), mom1a.normalized(),
          1_e / mom1a.norm(), covMat, ParticleHypothesis::pion())
          .value();
  auto params1b =
      BoundTrackParameters::create(
          geoContext, perigeeSurface, makeVector4(pos1b, 0), mom1b.normalized(),
          -1_e / mom1b.norm(), covMat, ParticleHypothesis::pion())
          .value();
  auto params1c =
      BoundTrackParameters::create(
          geoContext, perigeeSurface, makeVector4(pos1c, 0), mom1c.normalized(),
          -1_e / mom1c.norm(), covMat, ParticleHypothesis::pion())
          .value();

  // Vector of track parameters
  std::vector<InputTrack> vec1 = {InputTrack{&params1a}, InputTrack{&params1b},
                                  InputTrack{&params1c}};

  auto res = finder.find(vec1, vertexingOptions, state);

  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3 result = (*res).back().position();

    BOOST_CHECK_EQUAL(result[eX], constraintPos[eX]);
    BOOST_CHECK_EQUAL(result[eY], constraintPos[eY]);
    CHECK_CLOSE_ABS(result[eZ], expectedZResult, 0.001_mm);
  }
}

const double zVertexPos = 12.;
// x position
std::normal_distribution<double> xdist(1_mm, 0.1_mm);
// y position
std::normal_distribution<double> ydist(-0.7_mm, 0.1_mm);
// z1 position
std::normal_distribution<double> z1dist(zVertexPos * 1_mm, 1_mm);
// z2 position
std::normal_distribution<double> z2dist(-3_mm, 0.5_mm);
// Track pT distribution
std::uniform_real_distribution<double> pTDist(0.1_GeV, 100_GeV);
// Track phi distribution
std::uniform_real_distribution<double> phiDist(-std::numbers::pi,
                                               std::numbers::pi);
// Track eta distribution
std::uniform_real_distribution<double> etaDist(-4., 4.);

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation
///
BOOST_AUTO_TEST_CASE(track_density_finder_random_test) {
  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3 pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions vertexingOptions(geoContext, magFieldContext);
  GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect<&InputTrack::extractParameters>();
  TrackDensityVertexFinder finder{
      TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)}};
  auto state = finder.makeState(magFieldContext);

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundTrackParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3 pos(xdist(gen), ydist(gen), 0);

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    double charge = std::copysign(1., etaDist(gen));

    // project the position on the surface
    Vector3 direction = makeDirectionFromPhiEta(phi, eta);
    Intersection3D intersection =
        perigeeSurface->intersect(geoContext, pos, direction).closest();
    pos = intersection.position();

    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    pos[eZ] = ((i % 4) == 0) ? z2dist(gen) : z1dist(gen);

    trackVec.push_back(BoundTrackParameters::create(
                           geoContext, perigeeSurface, makeVector4(pos, 0),
                           direction, charge / pt, covMat,
                           ParticleHypothesis::pion())
                           .value());
  }

  std::vector<InputTrack> inputTracks;
  for (const auto& trk : trackVec) {
    inputTracks.emplace_back(&trk);
  }

  auto res3 = finder.find(inputTracks, vertexingOptions, state);
  if (!res3.ok()) {
    std::cout << res3.error().message() << std::endl;
  }

  if (res3.ok()) {
    BOOST_CHECK(!(*res3).empty());
    Vector3 result = (*res3).back().position();
    CHECK_CLOSE_ABS(result[eZ], zVertexPos, 1_mm);
  }
}

// Dummy user-defined InputTrackStub type
struct InputTrackStub {
  explicit InputTrackStub(const BoundTrackParameters& params)
      : m_parameters(params) {}

  const BoundTrackParameters& parameters() const { return m_parameters; }

  // store e.g. link to original objects here

 private:
  BoundTrackParameters m_parameters;
};

///
/// @brief Unit test for TrackDensityVertexFinder with user-defined input track
/// type with same values as in other tests
///
BOOST_AUTO_TEST_CASE(track_density_finder_usertrack_test) {
  // Define some track parameter properties
  Vector3 pos0{0, 0, 0};
  Vector3 pos1a{1.86052_mm, -1.24035_mm, -10_mm};
  Vector3 mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3 pos1b{-1.24035_mm, 1.86052_mm, -3_mm};
  Vector3 mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3 pos1c{1.69457_mm, -0.50837_mm, -7_mm};
  Vector3 mom1c{300_MeV, 1000_MeV, 100_MeV};

  // From Athena VertexSeedFinderTestAlg
  double const expectedZResult = -13.013;

  // Create constraint for seed finding
  Vector3 constraintPos{1.7_mm, 1.3_mm, -6_mm};
  SquareMatrix3 constrCov = SquareMatrix3::Identity();

  Vertex constraint(constraintPos);
  constraint.setCovariance(constrCov);

  // Finder options
  VertexingOptions vertexingOptions(geoContext, magFieldContext, constraint);

  auto extractParameters = [](const InputTrack& params) {
    return params.as<InputTrackStub>()->parameters();
  };

  GaussianTrackDensity::Config densityCfg;
  densityCfg.extractParameters.connect(extractParameters);
  TrackDensityVertexFinder finder{
      TrackDensityVertexFinder::Config{Acts::GaussianTrackDensity(densityCfg)}};
  auto state = finder.makeState(magFieldContext);

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  InputTrackStub params1a(BoundTrackParameters::create(
                              geoContext, perigeeSurface, makeVector4(pos1a, 0),
                              mom1a, 1_e / mom1a.norm(), covMat,
                              ParticleHypothesis::pion())
                              .value());
  InputTrackStub params1b(BoundTrackParameters::create(
                              geoContext, perigeeSurface, makeVector4(pos1b, 0),
                              mom1b, -1_e / mom1b.norm(), covMat,
                              ParticleHypothesis::pion())
                              .value());
  InputTrackStub params1c(BoundTrackParameters::create(
                              geoContext, perigeeSurface, makeVector4(pos1c, 0),
                              mom1c, -1_e / mom1c.norm(), covMat,
                              ParticleHypothesis::pion())
                              .value());

  // Vector of track parameters
  std::vector<InputTrack> vec1 = {InputTrack{&params1a}, InputTrack{&params1b},
                                  InputTrack{&params1c}};

  auto res = finder.find(vec1, vertexingOptions, state);

  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3 result = (*res).back().position();

    BOOST_CHECK_EQUAL(result[eX], constraintPos[eX]);
    BOOST_CHECK_EQUAL(result[eY], constraintPos[eY]);
    CHECK_CLOSE_ABS(result[eZ], expectedZResult, 0.001_mm);
  }
}

// Independent exhaustive search: evaluate every support interval at each trial
// point. The indexed implementation must preserve both the maximum and width
// exactly, including floating-point summation and tie-breaking order.
std::optional<std::pair<double, double>> exhaustiveMaximum(
    const GaussianTrackDensity::State& state, bool gaussian) {
  double maxZ = 0.;
  double maxDensity = 0.;
  double maxCurvature = 0.;
  for (const auto& trial : state.trackEntries) {
    double z = trial.z;
    for (int step = 0; step < 3; ++step) {
      double density = 0.;
      double slope = 0.;
      double curvature = 0.;
      for (const auto& track : state.trackEntries) {
        if (!(track.lowerBound < z && z < track.upperBound)) {
          continue;
        }
        const double value = std::exp(track.c0 + z * (track.c1 + z * track.c2));
        const double qPrime = track.c1 + 2. * z * track.c2;
        const double derivative = value * qPrime;
        density += value;
        slope += derivative;
        curvature += 2. * track.c2 * value + qPrime * derivative;
      }
      if (curvature >= 0. || density <= 0.) {
        break;
      }
      if (density > maxDensity) {
        maxZ = z;
        maxDensity = density;
        maxCurvature = curvature;
      }
      z += gaussian ? (density * slope) / (slope * slope - density * curvature)
                    : -slope / curvature;
    }
  }
  if (maxCurvature == 0.) {
    return std::nullopt;
  }
  return std::pair{maxZ, std::sqrt(-maxDensity / maxCurvature)};
}

BOOST_AUTO_TEST_CASE(track_density_interval_lookup_test) {
  std::mt19937 generator(981);
  std::uniform_real_distribution<double> position(-100., 100.);
  std::uniform_real_distribution<double> logWidth(-4., 4.);
  for (bool gaussian : {false, true}) {
    GaussianTrackDensity::Config config;
    config.isGaussianShaped = gaussian;
    config.extractParameters.connect<&InputTrack::extractParameters>();
    GaussianTrackDensity density(config);
    for (unsigned int size : {0u, 1u, 2u, 31u, 32u, 200u, 1258u}) {
      GaussianTrackDensity::State state(size);
      for (unsigned int i = 0; i < size; ++i) {
        const double z = position(generator);
        const double sigma = std::exp(logWidth(generator));
        const double c2 = -0.5 / (sigma * sigma);
        state.trackEntries.emplace_back(z, c2 * z * z, -2. * c2 * z, c2,
                                        z - 12. * sigma, z + 12. * sigma);
      }
      for (int repeat = 0; repeat < 3; ++repeat) {
        std::shuffle(state.trackEntries.begin(), state.trackEntries.end(),
                     generator);
        auto expected = exhaustiveMaximum(state, gaussian);
        auto actual = density.globalMaximumWithWidth(state, {});
        BOOST_REQUIRE(actual.ok());
        BOOST_REQUIRE_EQUAL(actual->has_value(), expected.has_value());
        if (expected) {
          BOOST_CHECK_EQUAL(actual->value().first, expected->first);
          BOOST_CHECK_EQUAL(actual->value().second, expected->second);
        }
      }
    }
    // Strict endpoints, equal lower bounds, nested and empty intervals, and
    // non-finite trials/bounds. NaN intervals must not enter a sorted lookup.
    GaussianTrackDensity::State state(80);
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    state.trackEntries = {{-1., 0., 0., -0.5, -1., 1.},
                          {1., 0., 0., -0.5, -1., 1.},
                          {0., 0., 0., -1., -0.1, 0.1},
                          {2., -2., 2., -0.5, -infinity, infinity},
                          {0., 0., 0., -1., 0., 0.},
                          {nan, 0., 0., -1., nan, 1.},
                          {infinity, 0., 0., -1., -1., nan}};
    // Force the indexed path, with support endpoints on either side of bin
    // boundaries and trial positions exactly on the boundaries.
    for (int i = -32; i <= 32; ++i) {
      const double z = i;
      state.trackEntries.emplace_back(z, -0.5 * z * z, z, -0.5,
                                      std::nextafter(z - 1., -infinity),
                                      std::nextafter(z + 1., infinity));
    }
    auto expected = exhaustiveMaximum(state, gaussian);
    auto actual = density.globalMaximumWithWidth(state, {});
    BOOST_REQUIRE(actual.ok());
    BOOST_REQUIRE(expected.has_value());
    BOOST_REQUIRE(actual->has_value());
    BOOST_CHECK_EQUAL(actual->value().first, expected->first);
    BOOST_CHECK_EQUAL(actual->value().second, expected->second);

    // All starting positions lie below the density maximum. Refinements must
    // still evaluate tracks outside the range covered by the lookup.
    state.trackEntries.clear();
    for (int i = 0; i < 40; ++i) {
      state.trackEntries.emplace_back(0.1 * i / 40., -0.02, 0.2, -0.5, -10.,
                                      10.);
    }
    expected = exhaustiveMaximum(state, gaussian);
    actual = density.globalMaximumWithWidth(state, {});
    BOOST_REQUIRE(actual.ok());
    BOOST_REQUIRE(expected.has_value());
    BOOST_REQUIRE(actual->has_value());
    BOOST_CHECK_EQUAL(actual->value().first, expected->first);
    BOOST_CHECK_EQUAL(actual->value().second, expected->second);

    // A degenerate trial range also falls back to the exhaustive scan.
    for (auto& entry : state.trackEntries) {
      entry.z = 0.;
    }
    expected = exhaustiveMaximum(state, gaussian);
    actual = density.globalMaximumWithWidth(state, {});
    BOOST_REQUIRE(actual.ok());
    BOOST_REQUIRE(expected.has_value());
    BOOST_REQUIRE(actual->has_value());
    BOOST_CHECK_EQUAL(actual->value().first, expected->first);
    BOOST_CHECK_EQUAL(actual->value().second, expected->second);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
