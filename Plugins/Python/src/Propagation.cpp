#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

namespace {
template <typename stepper_t>
void addStepper(const std::string& prefix, py::module_& m, py::module_& prop,
                py::module_& mex) {
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
                 mex, (prefix + "PropagationAlgorithm").c_str())
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("propagator"), py::arg("level"));

  auto c = py::class_<Config>(alg, "Config").def(py::init<propagator_t>());
#define _MEMBER(name) PY_MEMBER(c, Config, name)
  _MEMBER(propagator);
  _MEMBER(randomNumberSvc);
  _MEMBER(mode);
  _MEMBER(sterileLogger);
  _MEMBER(debugOutput);
  _MEMBER(energyLoss);
  _MEMBER(multipleScattering);
  _MEMBER(recordMaterialInteractions);
  _MEMBER(ntests);
  _MEMBER(d0Sigma);
  _MEMBER(z0Sigma);
  _MEMBER(phiSigma);
  _MEMBER(thetaSigma);
  _MEMBER(qpSigma);
  _MEMBER(tSigma);
  _MEMBER(phiRange);
  _MEMBER(etaRange);
  _MEMBER(ptRange);
  _MEMBER(ptLoopers);
  _MEMBER(maxStepSize);
  _MEMBER(propagationStepCollection);
  _MEMBER(propagationMaterialCollection);
  _MEMBER(covarianceTransport);
  _MEMBER(covariances);
  _MEMBER(correlations);
#undef _MEMBER
}

ACTS_PYTHON_COMPONENT(Propagation, ctx) {
  auto& [m, mex, prop] = ctx;
  py::class_<Acts::Navigator, std::shared_ptr<Acts::Navigator>>(m, "Navigator")
      .def(py::init<std::shared_ptr<const Acts::TrackingGeometry>>())
      .def_readwrite("resolveMaterial", &Acts::Navigator::resolveMaterial)
      .def_readwrite("resolvePassive", &Acts::Navigator::resolvePassive)
      .def_readwrite("resolveSensitive", &Acts::Navigator::resolveSensitive);

  addStepper<Acts::EigenStepper<>>("Eigen", m, prop, mex);
  addStepper<Acts::AtlasStepper>("Atlas", m, prop, mex);
  addStepper<Acts::StraightLineStepper>("StraightLine", m, prop, mex);
}

}  // namespace
