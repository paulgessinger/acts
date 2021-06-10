#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

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
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
#define _MEMBER(name) PY_MEMBER(c, Config, name)
    _MEMBER(inputParticles);
    _MEMBER(outputParticlesInitial);
    _MEMBER(outputParticlesFinal);
    _MEMBER(outputSimHits);
    _MEMBER(imputParametrisationNuclearInteraction);
    _MEMBER(randomNumbers);
    _MEMBER(trackingGeometry);
    _MEMBER(magneticField);
    _MEMBER(pMin);
    _MEMBER(emScattering);
    _MEMBER(emEnergyLossIonisation);
    _MEMBER(emEnergyLossRadiation);
    _MEMBER(emPhotonConversion);
    _MEMBER(generateHitsOnSensitive);
    _MEMBER(generateHitsOnMaterial);
    _MEMBER(generateHitsOnPassive);
    _MEMBER(averageHitsPerParticle);
#undef _MEMBER
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

    PY_MEMBER(c, Config, inputSimHits);
    PY_MEMBER(c, Config, outputSourceLinks);
    PY_MEMBER(c, Config, outputMeasurements);
    PY_MEMBER(c, Config, outputClusters);
    PY_MEMBER(c, Config, outputMeasurementParticlesMap);
    PY_MEMBER(c, Config, outputMeasurementSimHitsMap);
    PY_MEMBER(c, Config, trackingGeometry);
    PY_MEMBER(c, Config, randomNumbers);
    PY_MEMBER(c, Config, digitizationConfigs);

    // PY_MEMBER(c, Config, isSimpleSmearer);
    // PY_MEMBER(c, Config, doMerge);
    // PY_MEMBER(c, Config, mergeNsigma);
    // PY_MEMBER(c, Config, mergeCommonCorner);
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
    PY_MEMBER(c, Config, inputSourceLinks);
    PY_MEMBER(c, Config, inputMeasurements);
    PY_MEMBER(c, Config, outputSpacePoints);
    PY_MEMBER(c, Config, trackingGeometry);
    PY_MEMBER(c, Config, geometrySelection);

    mex.def("readJsonGeometryList", ActsExamples::readJsonGeometryList);
  }

  {
    using Config = ActsExamples::SeedingAlgorithm::Config;

    auto alg =
        py::class_<ActsExamples::SeedingAlgorithm, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::SeedingAlgorithm>>(
            mex, "SeedingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    PY_MEMBER(c, Config, inputSpacePoints);
    PY_MEMBER(c, Config, outputSeeds);
    PY_MEMBER(c, Config, outputProtoTracks);
    PY_MEMBER(c, Config, rMax);
    PY_MEMBER(c, Config, deltaRMin);
    PY_MEMBER(c, Config, deltaRMax);
    PY_MEMBER(c, Config, collisionRegionMin);
    PY_MEMBER(c, Config, collisionRegionMax);
    PY_MEMBER(c, Config, zMin);
    PY_MEMBER(c, Config, zMax);
    PY_MEMBER(c, Config, maxSeedsPerSpM);
    PY_MEMBER(c, Config, cotThetaMax);
    PY_MEMBER(c, Config, sigmaScattering);
    PY_MEMBER(c, Config, radLengthPerSeed);
    PY_MEMBER(c, Config, minPt);
    PY_MEMBER(c, Config, bFieldInZ);
    PY_MEMBER(c, Config, beamPosX);
    PY_MEMBER(c, Config, beamPosY);
    PY_MEMBER(c, Config, impactMax);

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
}

}  // namespace