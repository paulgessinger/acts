#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

PYBIND11_MODULE(ActsPythonBindingsDD4hep, m) {
  {
    using Config = ActsExamples::DD4hep::DD4hepGeometryService::Config;
    auto s = py::class_<DD4hep::DD4hepGeometryService,
                        std::shared_ptr<DD4hep::DD4hepGeometryService>>(
                 m, "DD4hepGeometryService")
                 .def(py::init<const Config&>());

    auto c = py::class_<Config>(s, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_MEMBER(xmlFileNames);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(bTypePhi);
    ACTS_PYTHON_MEMBER(bTypeR);
    ACTS_PYTHON_MEMBER(bTypeZ);
    ACTS_PYTHON_MEMBER(envelopeR);
    ACTS_PYTHON_MEMBER(envelopeZ);
    ACTS_PYTHON_MEMBER(defaultLayerThickness);
    ACTS_PYTHON_STRUCT_END();

    py::module::import("acts._adapter").attr("_patch_config_constructor")(c);
  }

  {
    auto gd = py::class_<DD4hepDetector, std::shared_ptr<DD4hepDetector>>(
                  m, "DD4hepDetector")
                  .def(py::init<>())
                  .def("finalize",
                       py::overload_cast<
                           DD4hep::DD4hepGeometryService::Config,
                           std::shared_ptr<const Acts::IMaterialDecorator>>(
                           &DD4hepDetector::finalize));
  }

  // auto parent = py::module::import("acts._adapter");
  // parent.attr("_patch_config")(m);
}