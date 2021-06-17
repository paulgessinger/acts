#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsExamples/Printers/HitsPrinter.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

using namespace ActsExamples;
using namespace Acts;

namespace {

ACTS_PYTHON_COMPONENT(ExampleAlgorithms, ctx) {
  auto& [m, mex, prop] = ctx;
  {
    using Config = ActsExamples::FatrasAlgorithm::Config;

    auto alg =
        py::class_<ActsExamples::FatrasAlgorithm, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::FatrasAlgorithm>>(
            mex, "FatrasAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config",
                                   &ActsExamples::FatrasAlgorithm::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputParticlesInitial);
    ACTS_PYTHON_MEMBER(outputParticlesFinal);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(imputParametrisationNuclearInteraction);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(pMin);
    ACTS_PYTHON_MEMBER(emScattering);
    ACTS_PYTHON_MEMBER(emEnergyLossIonisation);
    ACTS_PYTHON_MEMBER(emEnergyLossRadiation);
    ACTS_PYTHON_MEMBER(emPhotonConversion);
    ACTS_PYTHON_MEMBER(generateHitsOnSensitive);
    ACTS_PYTHON_MEMBER(generateHitsOnMaterial);
    ACTS_PYTHON_MEMBER(generateHitsOnPassive);
    ACTS_PYTHON_MEMBER(averageHitsPerParticle);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = ActsExamples::DigitizationConfig;

    py::class_<ActsExamples::DigitizationAlgorithm, ActsExamples::BareAlgorithm,
               std::shared_ptr<ActsExamples::DigitizationAlgorithm>>(
        mex, "DigitizationAlgorithm")
        .def(py::init<Config&, Acts::Logging::Level>(), py::arg("config"),
             py::arg("level"));

    auto c = py::class_<Config>(mex, "DigitizationConfig")
                 .def(py::init<Acts::GeometryHierarchyMap<
                          ActsExamples::DigiComponentsConfig>>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(outputSourceLinks);
    ACTS_PYTHON_MEMBER(outputMeasurements);
    ACTS_PYTHON_MEMBER(outputClusters);
    ACTS_PYTHON_MEMBER(outputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_MEMBER(digitizationConfigs);
    ACTS_PYTHON_STRUCT_END();

    c.def_readonly("isSimpleSmearer", &Config::isSimpleSmearer);
    c.def_readonly("doMerge", &Config::doMerge);
    c.def_readonly("mergeNsigma", &Config::mergeNsigma);
    c.def_readonly("mergeCommonCorner", &Config::mergeCommonCorner);

    py::class_<Acts::GeometryHierarchyMap<ActsExamples::DigiComponentsConfig>>(
        mex, "GeometryHierarchy_DigiComponentsConfig");

    mex.def("readDigiConfigFromJson", ActsExamples::readDigiConfigFromJson);
  }

  {
    using Config = ActsExamples::SpacePointMaker::Config;
    auto alg =
        py::class_<ActsExamples::SpacePointMaker, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::SpacePointMaker>>(
            mex, "SpacePointMaker")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config",
                                   &ActsExamples::SpacePointMaker::config);

    auto c = py::class_<ActsExamples::SpacePointMaker::Config>(alg, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(outputSpacePoints);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(geometrySelection);
    ACTS_PYTHON_STRUCT_END();

    mex.def("readJsonGeometryList", ActsExamples::readJsonGeometryList);
  }

  {
    using Config = ActsExamples::SeedingAlgorithm::Config;

    auto alg =
        py::class_<ActsExamples::SeedingAlgorithm, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::SeedingAlgorithm>>(
            mex, "SeedingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config",
                                   &ActsExamples::SeedingAlgorithm::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSpacePoints);
    ACTS_PYTHON_MEMBER(outputSeeds);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_MEMBER(beamPosX);
    ACTS_PYTHON_MEMBER(beamPosY);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_STRUCT_END();

    c.def_property(
        "deltaR",
        [](Config& cfg) {
          return std::pair{cfg.deltaRMin, cfg.deltaRMax};
        },
        [](Config& cfg, std::pair<double, double> values) {
          cfg.deltaRMin = values.first;
          cfg.deltaRMax = values.second;
        });
    c.def_property(
        "z",
        [](Config& cfg) {
          return std::pair{cfg.zMin, cfg.zMax};
        },
        [](Config& cfg, std::pair<double, double> values) {
          cfg.zMin = values.first;
          cfg.zMax = values.second;
        });

    c.def_property(
        "collisionRegion",
        [](Config& cfg) {
          return std::pair{cfg.collisionRegionMin, cfg.collisionRegionMax};
        },
        [](Config& cfg, std::pair<double, double> values) {
          cfg.collisionRegionMin = values.first;
          cfg.collisionRegionMax = values.second;
        });
  }

  {
    using Alg = ActsExamples::TrackParamsEstimationAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackParamsEstimationAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSeeds);
    ACTS_PYTHON_MEMBER(inputSpacePoints);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(outputTrackParameters);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(bFieldMin);
    ACTS_PYTHON_MEMBER(sigmaLoc0);
    ACTS_PYTHON_MEMBER(sigmaLoc1);
    ACTS_PYTHON_MEMBER(sigmaPhi);
    ACTS_PYTHON_MEMBER(sigmaTheta);
    ACTS_PYTHON_MEMBER(sigmaQOverP);
    ACTS_PYTHON_MEMBER(sigmaT0);
    ACTS_PYTHON_MEMBER(initialVarInflation);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::ParticlesPrinter;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "ParticlesPrinter")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    py::class_<Alg::Config>(alg, "Config")
        .def(py::init<>())
        .def_readwrite("inputParticles", &Alg::Config::inputParticles);
  }

  {
    using Alg = ActsExamples::HitsPrinter;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "HitsPrinter")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputHitIds);
    ACTS_PYTHON_MEMBER(selectIndexStart);
    ACTS_PYTHON_MEMBER(selectIndexLength);
    ACTS_PYTHON_MEMBER(selectVolume);
    ACTS_PYTHON_MEMBER(selectLayer);
    ACTS_PYTHON_MEMBER(selectModule);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TrackParametersPrinter;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackParametersPrinter")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    py::class_<Alg::Config>(alg, "Config")
        .def(py::init<>())
        .def_readwrite("inputTrackParameters",
                       &Alg::Config::inputTrackParameters);
  }

  {
    using Alg = ActsExamples::TruthSeedSelector;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TruthSeedSelector")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_MEMBER(rhoMin);
    ACTS_PYTHON_MEMBER(rhoMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(etaMin);
    ACTS_PYTHON_MEMBER(etaMax);
    ACTS_PYTHON_MEMBER(absEtaMin);
    ACTS_PYTHON_MEMBER(absEtaMax);
    ACTS_PYTHON_MEMBER(ptMin);
    ACTS_PYTHON_MEMBER(ptMax);
    ACTS_PYTHON_MEMBER(keepNeutral);
    ACTS_PYTHON_MEMBER(nHitsMin);
    ACTS_PYTHON_MEMBER(nHitsMax);
    ACTS_PYTHON_STRUCT_END();

    pythonRangeProperty(c, "rho", &Config::rhoMin, &Config::rhoMax);
    pythonRangeProperty(c, "z", &Config::zMin, &Config::zMax);
    pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
    pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
    pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
    pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
    pythonRangeProperty(c, "nHits", &Config::nHitsMin, &Config::nHitsMax);
  }
}

}  // namespace