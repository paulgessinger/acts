// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Detray/DetrayPropagator.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

// detray includes
#include <detray/core/detector.hpp>
#include <detray/detectors/itk_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/navigation/volume_graph.hpp>
#include <detray/propagator/line_stepper.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/propagator/rk_stepper.hpp>
#include <detray/test/common/bfield.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

#ifdef ACTS_DETRAY_NAVIGATION_VALIDATION
#include <detray/test/cpu/detector_scan.hpp>
#include <detray/test/cpu/navigation_validation.hpp>
#include <detray/test/framework/register_checks.hpp>
#include <detray/test/framework/whiteboard.hpp>
#include <gtest/gtest.h>
#endif

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPlugins;
using namespace ActsPython;

namespace {

using DetrayScalar = float;
using DetrayConstField =
    typename detray::bfield::const_field_t<DetrayScalar>::view_t;

/// Bind StraightLinePropagator<tag> and HelixPropagator<tag> for a detector
template <typename detector_t>
void bindPropagators(py::module_& m, const std::string& tag) {
  using Algebra = typename detector_t::algebra_type;

  m.def(
      ("StraightLinePropagator" + tag).c_str(),
      [](std::shared_ptr<detector_t> detector,
         vecmem::memory_resource& memoryResource, bool sterile,
         Logging::Level logLevel, std::array<unsigned int, 2> searchWindow)
          -> std::shared_ptr<PropagatorInterface> {
        using Stepper = detray::line_stepper<Algebra>;
        auto logger = getDefaultLogger("DetrayPropagator", logLevel);
        if (sterile) {
          return std::make_shared<DetraySterilePropagator<Stepper, detector_t>>(
              detector, memoryResource, std::move(logger), searchWindow);
        }
        return std::make_shared<DetrayPropagator<Stepper, detector_t>>(
            detector, memoryResource, std::move(logger), searchWindow);
      },
      "detrayDetector"_a, "memoryResource"_a, "sterile"_a,
      "logLevel"_a = Logging::INFO,
      "searchWindow"_a = std::array<unsigned int, 2>{0u, 0u});

  m.def(("HelixPropagator" + tag).c_str(),
        [](std::shared_ptr<detector_t> detector,
           vecmem::memory_resource& memoryResource, bool sterile,
           std::array<double, 3> bField, Logging::Level logLevel,
           std::array<unsigned int, 2> searchWindow)
            -> std::shared_ptr<PropagatorInterface> {
          using Stepper = detray::rk_stepper<DetrayConstField, Algebra>;
          // ACTS and detray share native field units; the constant backend's
          // view holds the field value itself
          const detray::dvector3D<Algebra> b{
              static_cast<DetrayScalar>(bField[0]),
              static_cast<DetrayScalar>(bField[1]),
              static_cast<DetrayScalar>(bField[2])};
          const auto field = detray::create_const_field<DetrayScalar>(b);
          DetrayConstField view(field);
          auto logger = getDefaultLogger("DetrayPropagator", logLevel);
          if (sterile) {
            return std::make_shared<
                DetraySterilePropagator<Stepper, detector_t, DetrayConstField>>(
                detector, memoryResource, std::move(logger), searchWindow,
                view);
          }
          return std::make_shared<
              DetrayPropagator<Stepper, detector_t, DetrayConstField>>(
              detector, memoryResource, std::move(logger), searchWindow, view);
        },
        "detrayDetector"_a, "memoryResource"_a, "sterile"_a, "bField"_a,
        "logLevel"_a = Logging::INFO,
        "searchWindow"_a = std::array<unsigned int, 2>{0u, 0u});
}

#ifdef ACTS_DETRAY_NAVIGATION_VALIDATION

/// Options of the detray navigation validation, see runNavigationValidation
struct NavigationValidationConfig {
  bool rays = true;
  bool helices = true;
  std::size_t nTracks = 10000;
  std::uint64_t seed = 42;
  std::array<double, 2> etaRange = {-4., 4.};
  std::array<double, 2> pTRange = {1. * UnitConstants::GeV,
                                   1. * UnitConstants::GeV};
  std::array<double, 3> origin = {0., 0., 0.};
  /// Gaussian vertex smearing, no smearing if all zero
  std::array<double, 3> originStddev = {0., 0., 0.};
  std::array<double, 3> bField = {0., 0., 2. * UnitConstants::T};
  std::array<unsigned int, 2> searchWindow = {0u, 0u};
  std::string dataDir = ".";
  bool sensitivesOnly = false;
  /// Write an SVG of the detector for every failed track; for large
  /// detectors these are tens of MB each
  bool displayOnFailure = false;
};

/// Scan the detector with rays and/or helices by intersecting every surface,
/// then propagate the same tracks with the detray navigator and compare.
/// Runs through gtest, so it can only be called once per process.
template <typename detector_t>
int runNavigationValidation(const detector_t& det,
                            const detray::name_map& names,
                            const NavigationValidationConfig& cfg) {
  static bool alreadyRun = false;
  if (alreadyRun) {
    throw std::runtime_error(
        "The navigation validation can only run once per process");
  }
  alreadyRun = true;

  using Scalar = detray::dscalar<typename detector_t::algebra_type>;
  using Vector3 = detray::dvector3D<typename detector_t::algebra_type>;

  const Vector3 bField{static_cast<Scalar>(cfg.bField[0]),
                       static_cast<Scalar>(cfg.bField[1]),
                       static_cast<Scalar>(cfg.bField[2])};

  typename detray::test::ray_scan<detector_t>::config rayScanCfg{};
  typename detray::test::helix_scan<detector_t>::config helScanCfg{};
  typename detray::test::straight_line_navigation<detector_t>::config
      rayNavCfg{};
  typename detray::test::helix_navigation<detector_t>::config helNavCfg{};

  auto& trackGen = helScanCfg.track_generator();
  trackGen.n_tracks(cfg.nTracks)
      .seed(cfg.seed)
      .eta_range(static_cast<Scalar>(cfg.etaRange[0]),
                 static_cast<Scalar>(cfg.etaRange[1]))
      .pT_range(static_cast<Scalar>(cfg.pTRange[0] / UnitConstants::GeV),
                static_cast<Scalar>(cfg.pTRange[1] / UnitConstants::GeV))
      .origin(static_cast<Scalar>(cfg.origin[0]),
              static_cast<Scalar>(cfg.origin[1]),
              static_cast<Scalar>(cfg.origin[2]))
      .randomize_charge(true);
  if (cfg.originStddev != std::array<double, 3>{0., 0., 0.}) {
    trackGen.origin_stddev(static_cast<Scalar>(cfg.originStddev[0]),
                           static_cast<Scalar>(cfg.originStddev[1]),
                           static_cast<Scalar>(cfg.originStddev[2]));
  }
  rayScanCfg.track_generator() = trackGen;

  const std::string detName = det.name(names);
  const std::string prefix = cfg.dataDir + "/" + detName;

  rayScanCfg.display_on_failure(cfg.displayOnFailure);
  helScanCfg.display_on_failure(cfg.displayOnFailure);
  rayNavCfg.display_on_failure(cfg.displayOnFailure);
  helNavCfg.display_on_failure(cfg.displayOnFailure);

  rayScanCfg.name(detName + "_ray_scan");
  rayScanCfg.intersection_file(prefix + "_ray_scan_intersections");
  rayScanCfg.track_param_file(prefix + "_ray_scan_track_parameters");

  helScanCfg.name(detName + "_helix_scan");
  helScanCfg.B_vector(bField);
  // Let the Newton solver choose its tolerance from the approximation error
  helScanCfg.mask_tolerance(detray::detail::invalid_value<Scalar>());
  helScanCfg.intersection_file(prefix + "_helix_scan_intersections");
  helScanCfg.track_param_file(prefix + "_helix_scan_track_parameters");

  rayNavCfg.name(detName + "_straight_line_navigation");
  rayNavCfg.fail_on_diff(false);
  rayNavCfg.collect_sensitives_only(cfg.sensitivesOnly);
  rayNavCfg.n_tracks(trackGen.n_tracks());
  rayNavCfg.propagation().navigation.search_window = cfg.searchWindow;
  // Same mask tolerance as the scan, so both find the same intersections
  rayNavCfg.propagation().navigation.intersection.min_mask_tolerance =
      static_cast<float>(rayScanCfg.mask_tolerance());
  rayNavCfg.propagation().navigation.intersection.max_mask_tolerance =
      static_cast<float>(rayScanCfg.mask_tolerance());
  rayNavCfg.intersection_file(rayScanCfg.intersection_file());
  rayNavCfg.track_param_file(rayScanCfg.track_param_file());

  helNavCfg.name(detName + "_helix_navigation");
  helNavCfg.fail_on_diff(false);
  helNavCfg.collect_sensitives_only(cfg.sensitivesOnly);
  helNavCfg.n_tracks(trackGen.n_tracks());
  helNavCfg.p_range(trackGen.mom_range());
  helNavCfg.B_vector(bField);
  helNavCfg.propagation().navigation.search_window = cfg.searchWindow;
  helNavCfg.intersection_file(helScanCfg.intersection_file());
  helNavCfg.track_param_file(helScanCfg.track_param_file());

  typename detector_t::geometry_context gctx{};
  auto whiteBoard = std::make_shared<detray::test::whiteboard>();

  ::testing::InitGoogleTest();
  if (cfg.rays) {
    detray::test::register_checks<detray::test::ray_scan>(
        det, names, rayScanCfg, gctx, whiteBoard);
    detray::test::register_checks<detray::test::straight_line_navigation>(
        det, names, rayNavCfg, gctx, whiteBoard);
  }
  if (cfg.helices) {
    detray::test::register_checks<detray::test::helix_scan>(
        det, names, helScanCfg, gctx, whiteBoard);
    detray::test::register_checks<detray::test::helix_navigation>(
        det, names, helNavCfg, gctx, whiteBoard);
  }
  return RUN_ALL_TESTS();
}

#endif

}  // namespace

PYBIND11_MODULE(ActsExamplesPythonBindingsDetray, detray) {
  using DetrayDetectorODD =
      detray::host::detector<detray::odd_metadata<detray::array<float>>>;
  using DetrayDetectorITk =
      detray::host::detector<detray::itk_metadata<detray::array<float>>>;

  bindPropagators<DetrayDetectorODD>(detray, "ODD");
  bindPropagators<DetrayDetectorITk>(detray, "ITk");

#ifdef ACTS_DETRAY_NAVIGATION_VALIDATION
  using Cfg = NavigationValidationConfig;
  py::class_<Cfg>(detray, "NavigationValidationConfig")
      .def(py::init<>())
      .def_readwrite("rays", &Cfg::rays)
      .def_readwrite("helices", &Cfg::helices)
      .def_readwrite("nTracks", &Cfg::nTracks)
      .def_readwrite("seed", &Cfg::seed)
      .def_readwrite("etaRange", &Cfg::etaRange)
      .def_readwrite("pTRange", &Cfg::pTRange)
      .def_readwrite("origin", &Cfg::origin)
      .def_readwrite("originStddev", &Cfg::originStddev)
      .def_readwrite("bField", &Cfg::bField)
      .def_readwrite("searchWindow", &Cfg::searchWindow)
      .def_readwrite("dataDir", &Cfg::dataDir)
      .def_readwrite("sensitivesOnly", &Cfg::sensitivesOnly)
      .def_readwrite("displayOnFailure", &Cfg::displayOnFailure);

  detray.def("runNavigationValidationODD",
             &runNavigationValidation<DetrayDetectorODD>, "detector"_a,
             "names"_a, "config"_a);
  detray.def("runNavigationValidationITk",
             &runNavigationValidation<DetrayDetectorITk>, "detector"_a,
             "names"_a, "config"_a);
#endif
}
