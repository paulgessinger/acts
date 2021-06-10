#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace {
ACTS_PYTHON_COMPONENT(DD4hep, ctx) {
  auto& [m, mex, prop] = ctx;
  {
    using Config = ActsExamples::DD4hep::DD4hepGeometryService::Config;

    auto gd =
        py::class_<DD4hepDetector, std::shared_ptr<DD4hepDetector>>(
            mex, "DD4hepDetector")
            .def(py::init<>())
            .def("finalize",
                 py::overload_cast<
                     Config, std::shared_ptr<const Acts::IMaterialDecorator>>(
                     &DD4hepDetector::finalize));

    py::class_<Config>(gd, "Config")
        .def(py::init<>())
        .def_readwrite("logLevel", &Config::logLevel)
        .def_readwrite("xmlFileNames", &Config::xmlFileNames)
        .def_readwrite("name", &Config::name)
        .def_readwrite("bTypePhi", &Config::bTypePhi)
        .def_readwrite("bTypeR", &Config::bTypeR)
        .def_readwrite("bTypeZ", &Config::bTypeZ)
        .def_readwrite("envelopeR", &Config::envelopeR)
        .def_readwrite("envelopeZ", &Config::envelopeZ)
        .def_readwrite("defaultLayerThickness", &Config::defaultLayerThickness);
  }
}
}  // namespace