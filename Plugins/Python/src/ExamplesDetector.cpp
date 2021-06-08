#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

void addExamplesDetector(py::module_& m) {
  {
    py::class_<IContextDecorator, std::shared_ptr<IContextDecorator>>(
        m, "IContextDecorator");
  }

  {
    using Config = GenericDetector::Config;

    auto gd = py::class_<GenericDetector, std::shared_ptr<GenericDetector>>(
                  m, "GenericDetector")
                  .def(py::init<>())
                  .def("finalize",
                       py::overload_cast<
                           const Config&,
                           std::shared_ptr<const Acts::IMaterialDecorator>>(
                           &GenericDetector::finalize));

    py::class_<Config>(gd, "Config")
        .def(py::init<>())
        .def_readwrite("buildLevel", &Config::buildLevel)
        .def_readwrite("surfaceLogLevel", &Config::surfaceLogLevel)
        .def_readwrite("layerLogLevel", &Config::layerLogLevel)
        .def_readwrite("volumeLogLevel", &Config::volumeLogLevel)
        .def_readwrite("buildProto", &Config::buildProto);
  }

  {
    using Config = ActsExamples::DD4hep::DD4hepGeometryService::Config;

    auto gd = py::class_<DD4hepDetector, std::shared_ptr<DD4hepDetector>>(
                  m, "DD4hepDetector")
                  .def(py::init<>())
                  .def("finalize",
                       py::overload_cast<
                           const Config&,
                           std::shared_ptr<const Acts::IMaterialDecorator>>(
                           &GenericDetector::finalize));

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