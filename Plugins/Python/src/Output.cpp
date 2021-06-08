#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

void addOutput(py::module_& m) {
  {
    using Writer = ActsExamples::ObjPropagationStepsWriter<Acts::detail::Step>;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 m, "ObjPropagationStepsWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level") = Acts::Logging::INFO);

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("collection", &Writer::Config::collection)
        .def_readwrite("outputDir", &Writer::Config::outputDir)
        .def_readwrite("outputScalor", &Writer::Config::outputScalor)
        .def_readwrite("outputPrecision", &Writer::Config::outputPrecision);
  }

  {
    using Writer = ActsExamples::RootPropagationStepsWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 m, "RootPropagationStepsWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level") = Acts::Logging::INFO);

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("collection", &Writer::Config::collection)
        .def_readwrite("filePath", &Writer::Config::filePath)
        .def_readwrite("fileMode", &Writer::Config::fileMode)
        .def_readwrite("treeName", &Writer::Config::treeName);
  }
}