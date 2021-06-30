#include "Acts/Plugins/Python/ActsModule.hpp"
#include "Acts/Utilities/PolymorphicValue.hpp"
#include "ActsExamples/Geant4/GeantinoRecording.hpp"
#include "ActsExamples/Geant4/PrimaryGeneratorAction.hpp"

#include <memory>

#include <G4VUserDetectorConstruction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

PYBIND11_DECLARE_HOLDER_TYPE(T, Acts::PolymorphicValue<T>)

void addGeant4DD4hep(py::module_& m);

namespace {

ACTS_PYTHON_COMPONENT(Geant4, ctx) {
  auto& [m, mex, prop] = ctx;

  auto g4 = mex.def_submodule("geant4");

  {
    using Alg = GeantinoRecording;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            g4, "GeantinoRecording")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(outputMaterialTracks);
    ACTS_PYTHON_MEMBER(detectorConstruction);
    ACTS_PYTHON_MEMBER(tracksPerEvent);
    ACTS_PYTHON_MEMBER(generationConfig);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto cls = py::class_<PrimaryGeneratorAction>(g4, "PrimaryGeneratorAction");
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

  py::class_<G4VUserDetectorConstruction,
             Acts::PolymorphicValue<G4VUserDetectorConstruction>>(
      g4, "G4VUserDetectorConstruction");

  addGeant4DD4hep(g4);
  py::module::import("acts._adapter").attr("_patch_config")(g4);
}

}  // namespace