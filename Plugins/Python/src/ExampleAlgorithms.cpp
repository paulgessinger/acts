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
void addStepper(const std::string& prefix, py::module_& prop) {
  using propagator_t = Acts::Propagator<stepper_t, Acts::Navigator>;
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

void addExamplesAlgorithms(py::module_& mex, py::module_& prop) {
  auto iAlgorithm =
      py::class_<ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::IAlgorithm>>(mex, "IAlgorithm");

  auto bareAlgorithm =
      py::class_<ActsExamples::BareAlgorithm, ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::BareAlgorithm>>(mex,
                                                               "BareAlgorithm");

  addStepper<Acts::EigenStepper<>>("Eigen", prop);
  addStepper<Acts::AtlasStepper>("Atlas", prop);
  addStepper<Acts::StraightLineStepper>("StraightLine", prop);
}