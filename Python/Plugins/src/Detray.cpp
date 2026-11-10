// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetrayJsonConverter.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceGridPolicy.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <optional>
#include <string>

#include <detray/detectors/itk_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/io/frontend/detector_reader_config.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/io/json/detector_reader.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

using namespace ActsPlugins;

/// Bind the detector type and its JSON reader for one metadata as
/// DetrayDetector<tag> and read<tag>
template <typename metadata_t>
void bindDetector(py::module_& m, const std::string& tag) {
  using Detector = detray::host::detector<metadata_t>;

  // Registered so shared_ptr<Detector> can be passed to functions in other
  // pybind11 modules (e.g. StraightLinePropagator<tag>)
  py::class_<Detector, std::shared_ptr<Detector>>(
      m, ("DetrayDetector" + tag).c_str())
      .def("volumes", [](Detector& self) { return self.volumes(); })
      .def("surfaces", [](Detector& self) { return self.surfaces(); })
      .def("checkConsistency",
           [](Detector& self) { detray::detail::check_consistency(self); })
      .def("writeToJson",
           [](Detector& self, const typename Detector::name_map& names,
              const std::string& source, const std::string& fname) {
             auto cfg = detray::io::detector_writer_config{}
                            .format(detray::io::format::json)
                            .path(fname)
                            .source(source)
                            .replace_files(true);
             detray::io::write_detector(self, names, cfg);
           });

  m.def(("read" + tag).c_str(),
        [](vecmem::memory_resource& mr, const std::vector<std::string>& files) {
          detray::io::detector_reader_config cfg{};
          cfg.do_check(false).add_files(files);
          return detray::io::read_detector_json<Detector>(mr, cfg);
        },
        "mr"_a, "files"_a);
}

/// Bind the conversion result for one metadata as @p name, and the
/// conversion itself as the converter method @p method
template <typename metadata_t, typename converter_class_t>
void bindConversion(converter_class_t& converter, const std::string& name,
                    const std::string& method) {
  using Geometry = DetrayGeometryConverter::DetrayGeometry<metadata_t>;
  py::class_<Geometry>(converter, name.c_str())
      .def_readonly("detector", &Geometry::detector)
      .def_readonly("names", &Geometry::names)
      .def_readonly("surfaceIds", &Geometry::surfaceIds);

  converter.def(
      method.c_str(),
      [](const DetrayGeometryConverter& self, vecmem::memory_resource& mr,
         const Acts::GeometryContext& gctx,
         std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
        return self.convert<metadata_t>(mr, gctx, std::move(trackingGeometry));
      },
      "mr"_a, "gctx"_a, "trackingGeometry"_a);
}

}  // namespace

PYBIND11_MODULE(ActsPluginsPythonBindingsDetray, detray) {
  using namespace ActsPlugins;

  using DetrayMetaDataODD = detray::odd_metadata<detray::array<float>>;
  using DetrayDetectorODD = detray::host::detector<DetrayMetaDataODD>;

  using DetrayMetaDataITk = detray::itk_metadata<detray::array<float>>;

  // Shared by all detector types; the name predates the ITk type
  py::class_<DetrayDetectorODD::name_map>(detray, "DetrayDetectorODDNameMap");

  bindDetector<DetrayMetaDataODD>(detray, "ODD");
  bindDetector<DetrayMetaDataITk>(detray, "ITk");

  // ── Surface grid policy ────────────────────────────────────────────────
  {
    py::class_<DetraySurfaceGridConfig>(detray, "DetraySurfaceGridConfig")
        .def(py::init<>())
        .def_readwrite("vertexZHalfRange",
                       &DetraySurfaceGridConfig::vertexZHalfRange)
        .def_readwrite("minPt", &DetraySurfaceGridConfig::minPt)
        .def_readwrite("bField", &DetraySurfaceGridConfig::bField)
        .def_readwrite("nVertexZ", &DetraySurfaceGridConfig::nVertexZ)
        .def_readwrite("nCurvature", &DetraySurfaceGridConfig::nCurvature)
        .def_readwrite("maxSampleSpacing",
                       &DetraySurfaceGridConfig::maxSampleSpacing)
        .def_readwrite("numberOfBinsFactor",
                       &DetraySurfaceGridConfig::numberOfBinsFactor)
        .def_readwrite("lookupOnLayer",
                       &DetraySurfaceGridConfig::lookupOnLayer);

    py::class_<DetraySurfaceGridPolicy>(detray, "DetraySurfaceGridPolicy")
        .def_property_readonly("config", &DetraySurfaceGridPolicy::config);

    detray.def(
        "addDetraySurfaceGridPolicy",
        [](Acts::NavigationPolicyFactory& factory,
           const DetraySurfaceGridConfig& config) {
          return std::move(factory).template add<DetraySurfaceGridPolicy>(
              config);
        },
        "factory"_a, "config"_a = DetraySurfaceGridConfig{},
        "Add a DetraySurfaceGridPolicy carrying @p config to the factory");

    detray.def("addDetrayNavigationPolicies", &addDetrayNavigationPolicies,
               "config"_a,
               "Register the Detray navigation policies with a "
               "TrackingGeometryJsonConverter config");
  }

  // ── Payload converter ──────────────────────────────────────────────────
  auto payloadConverter = py::class_<DetrayPayloadConverter,
                                     std::shared_ptr<DetrayPayloadConverter>>(
      detray, "DetrayPayloadConverter");

  auto payloadConfig =
      py::class_<DetrayPayloadConverter::Config>(payloadConverter, "Config")
          .def(py::init<>())
          .def_readwrite("sensitiveStrategy",
                         &DetrayPayloadConverter::Config::sensitiveStrategy)
          .def_readwrite("includeWorld",
                         &DetrayPayloadConverter::Config::includeWorld)
          .def_readwrite(
              "portalSegmentTolerance",
              &DetrayPayloadConverter::Config::portalSegmentTolerance)
          .def_readwrite("convertMaterial",
                         &DetrayPayloadConverter::Config::convertMaterial)
          .def_readwrite("convertSurfaceGrids",
                         &DetrayPayloadConverter::Config::convertSurfaceGrids)
          .def_property(
              "detectorName",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.detectorName;
              },
              [](DetrayPayloadConverter::Config& cfg, const std::string& name) {
                cfg.detectorName = name;
              },
              py::return_value_policy::reference)
          .def_property(
              "beampipeVolume",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.beampipeVolume;
              },
              [](DetrayPayloadConverter::Config& cfg,
                 const Acts::TrackingVolume* volume) {
                cfg.beampipeVolume = volume;
              },
              py::return_value_policy::reference)
          .def_property(
              "replacementWorld",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.replacementWorld;
              },
              [](DetrayPayloadConverter::Config& cfg,
                 const Acts::TrackingVolume* volume) {
                cfg.replacementWorld = volume;
              },
              py::return_value_policy::reference);

  payloadConfig.def_readwrite("surfaceGridFill",
                              &DetrayPayloadConverter::Config::surfaceGridFill);

  py::enum_<DetrayPayloadConverter::Config::SensitiveStrategy>(
      payloadConfig, "SensitiveStrategy")
      .value("Identifier",
             DetrayPayloadConverter::Config::SensitiveStrategy::Identifier)
      .value(
          "DetectorElement",
          DetrayPayloadConverter::Config::SensitiveStrategy::DetectorElement);

  payloadConverter.def(
      py::init([](const DetrayPayloadConverter::Config& config,
                  Acts::Logging::Level level) {
        return std::make_shared<DetrayPayloadConverter>(
            config, Acts::getDefaultLogger("DetrayPayloadConverter", level));
      }),
      "config"_a, "level"_a = Acts::Logging::INFO);

  // ── Geometry converter ─────────────────────────────────────────────────
  auto geometryConverter = py::class_<DetrayGeometryConverter,
                                      std::shared_ptr<DetrayGeometryConverter>>(
      detray, "DetrayGeometryConverter");

  py::class_<DetrayGeometryConverter::Config>(geometryConverter, "Config")
      .def(py::init<>())
      .def_readwrite("payloadConverter",
                     &DetrayGeometryConverter::Config::payloadConverter)
      .def_readwrite("deduplicateMaterial",
                     &DetrayGeometryConverter::Config::deduplicateMaterial);

  py::class_<DetraySurfaceIdMap, std::shared_ptr<DetraySurfaceIdMap>>(
      detray, "DetraySurfaceIdMap")
      // detray identifiers are exposed as their integer value
      .def(
          "toDetray",
          [](const DetraySurfaceIdMap& self, Acts::GeometryIdentifier id)
              -> std::optional<detray::geometry::identifier::value_t> {
            if (auto detrayId = self.toDetray(id)) {
              return detrayId->value();
            }
            return std::nullopt;
          },
          "id"_a)
      .def(
          "fromDetray",
          [](const DetraySurfaceIdMap& self,
             detray::geometry::identifier::value_t detrayId) {
            return self.fromDetray(detray::geometry::identifier{detrayId});
          },
          "detrayId"_a)
      .def("__len__", &DetraySurfaceIdMap::size);

  geometryConverter.def(
      py::init([](DetrayGeometryConverter::Config config,
                  Acts::Logging::Level level) {
        return std::make_shared<DetrayGeometryConverter>(
            std::move(config),
            Acts::getDefaultLogger("DetrayGeometryConverter", level));
      }),
      "config"_a, "level"_a = Acts::Logging::INFO);

  bindConversion<DetrayMetaDataODD>(geometryConverter, "DetrayGeometry",
                                    "convert");
  bindConversion<DetrayMetaDataITk>(geometryConverter, "DetrayGeometryITk",
                                    "convertITk");
}
