#include "ActsExamples/Geant4/GeantinoRecording.hpp"
#include "ActsExamples/Geant4/PrimaryGeneratorAction.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace {

ACTS_PYTHON_COMPONENT(Geant4, ctx) {
  auto& [m, mex, prop] = ctx;

  {
    using Alg = GeantinoRecording;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "GeantinoRecording")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(outputMaterialTracks);
    ACTS_PYTHON_MEMBER(gdmlInputPath);
    ACTS_PYTHON_MEMBER(tracksPerEvent);
    ACTS_PYTHON_MEMBER(generationConfig);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto cls =
        py::class_<PrimaryGeneratorAction>(mex, "PrimaryGeneratorAction");
    auto c = py::class_<PrimaryGeneratorAction::Config>(cls, "Config")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, PrimaryGeneratorAction::Config);
    ACTS_PYTHON_MEMBER(particleName);
    ACTS_PYTHON_MEMBER(energy);
    ACTS_PYTHON_MEMBER(randomSeed1);
    ACTS_PYTHON_MEMBER(randomSeed2);
    ACTS_PYTHON_MEMBER(vertexPosX);
    ACTS_PYTHON_MEMBER(vertexPosY);
    ACTS_PYTHON_MEMBER(vertexPosZ);
    ACTS_PYTHON_MEMBER(vertexSigmaX);
    ACTS_PYTHON_MEMBER(vertexSigmaY);
    ACTS_PYTHON_MEMBER(vertexSigmaZ);
    ACTS_PYTHON_MEMBER(phiRange);
    ACTS_PYTHON_MEMBER(etaRange);
    ACTS_PYTHON_MEMBER(samplingVariable);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace