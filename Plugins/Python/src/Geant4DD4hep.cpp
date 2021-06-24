#include "Acts/Utilities/PolymorphicValue.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/Geant4DD4hep/DD4hepDetectorConstruction.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <G4VUserDetectorConstruction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

PYBIND11_DECLARE_HOLDER_TYPE(T, Acts::PolymorphicValue<T>)

namespace {

ACTS_PYTHON_COMPONENT(Geant4, ctx) {
  auto& [m, mex, prop] = ctx;

  auto cc = py::class_<DD4hepDetectorConstruction, G4VUserDetectorConstruction,
                       Acts::PolymorphicValue<DD4hepDetectorConstruction>>(
      mex, "DD4hepDetectorConstruction");

  // add a special constructor so we don't have to expose the internals of
  // DD4hepGeometryService
  cc.def(py::init([](DD4hep::DD4hepGeometryService& geometrySvc) {
    return DD4hepDetectorConstruction{*geometrySvc.lcdd()};
  }));
}

}  // namespace