#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
// using namespace ActsExamples;

template <typename stepper_t>
void addStepper(const std::string& prefix, py::module_& m, py::module_& prop) {
  auto stepper = py::class_<stepper_t>(m, (prefix + "Stepper").c_str());

  if constexpr (std::is_same_v<stepper_t, Acts::StraightLineStepper>) {
    stepper.def(py::init<>());
  } else {
    stepper.def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());
  }

  using propagator_t = Acts::Propagator<stepper_t, Acts::Navigator>;
  py::class_<propagator_t>(prop, (prefix + "Propagator").c_str())
      .def(py::init<stepper_t, Acts::Navigator>());

  using Algorithm = typename ActsExamples::PropagationAlgorithm<propagator_t>;
  using Config = typename Algorithm::Config;
  auto alg = py::class_<Algorithm, ActsExamples::BareAlgorithm,
                        std::shared_ptr<Algorithm>>(
                 prop, (prefix + "PropagationAlgorithm").c_str())
                 .def(py::init<const Config&, Acts::Logging::Level>());

  py::class_<Config>(alg, "Config")
      .def(py::init<propagator_t>())
      .def_readwrite("randomNumberSvc", &Config::randomNumberSvc)
      .def_readwrite("propagator", &Config::propagator)
      .def_readwrite("ntests", &Config::ntests);
}

void addExamplesAlgorithms(py::module_& m) {
  py::class_<Acts::Navigator, std::shared_ptr<Acts::Navigator>>(m, "Navigator")
      .def(py::init<std::shared_ptr<const Acts::TrackingGeometry>>())
      .def_readwrite("resolveMaterial", &Acts::Navigator::resolveMaterial)
      .def_readwrite("resolvePassive", &Acts::Navigator::resolvePassive)
      .def_readwrite("resolveSensitive", &Acts::Navigator::resolveSensitive);

  auto iAlgorithm =
      py::class_<ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::IAlgorithm>>(m, "IAlgorithm");

  auto bareAlgorithm =
      py::class_<ActsExamples::BareAlgorithm, ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::BareAlgorithm>>(m,
                                                               "BareAlgorithm");

  auto prop = m.def_submodule("propagator");

  addStepper<Acts::EigenStepper<>>("Eigen", m, prop);
  addStepper<Acts::AtlasStepper>("Atlas", m, prop);
  addStepper<Acts::StraightLineStepper>("StraightLine", m, prop);
}