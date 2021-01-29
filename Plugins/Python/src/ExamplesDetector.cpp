#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

void addExamplesDetector(py::module_& m) {
  using Config = GenericDetector::Config;

  auto gd =
      py::class_<GenericDetector, std::shared_ptr<GenericDetector>>(
          m, "GenericDetector")
          .def(py::init<>())
          .def("finalize", py::overload_cast<
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