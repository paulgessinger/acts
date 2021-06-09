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
}

void addPropagation(py::module_& m, py::module_& prop) {
  py::class_<Acts::Navigator, std::shared_ptr<Acts::Navigator>>(m, "Navigator")
      .def(py::init<std::shared_ptr<const Acts::TrackingGeometry>>())
      .def_readwrite("resolveMaterial", &Acts::Navigator::resolveMaterial)
      .def_readwrite("resolvePassive", &Acts::Navigator::resolvePassive)
      .def_readwrite("resolveSensitive", &Acts::Navigator::resolveSensitive);

  addStepper<Acts::EigenStepper<>>("Eigen", m, prop);
  addStepper<Acts::AtlasStepper>("Atlas", m, prop);
  addStepper<Acts::StraightLineStepper>("StraightLine", m, prop);
}