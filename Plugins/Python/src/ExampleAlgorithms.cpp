#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
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
}

}  // namespace