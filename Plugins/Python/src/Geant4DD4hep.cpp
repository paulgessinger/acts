#include "Acts/Plugins/Python/ActsModule.hpp"
#include "Acts/Utilities/PolymorphicValue.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/Geant4DD4hep/DD4hepDetectorConstruction.hpp"

#include <memory>

#include <G4VUserDetectorConstruction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

PYBIND11_DECLARE_HOLDER_TYPE(T, Acts::PolymorphicValue<T>)

void addGeant4DD4hep(py::module_& m) {
  auto cc = py::class_<DD4hepDetectorConstruction, G4VUserDetectorConstruction,
                       Acts::PolymorphicValue<DD4hepDetectorConstruction>>(
      m, "DD4hepDetectorConstruction");

  // add a special constructor so we don't have to expose the internals of
  // DD4hepGeometryService
  cc.def(py::init([](DD4hep::DD4hepGeometryService& geometrySvc) {
    return DD4hepDetectorConstruction{*geometrySvc.lcdd()};
  }));
}
