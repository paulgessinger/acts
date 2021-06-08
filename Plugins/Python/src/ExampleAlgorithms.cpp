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
#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

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