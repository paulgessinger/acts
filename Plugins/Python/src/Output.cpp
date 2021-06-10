#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace {
ACTS_PYTHON_COMPONENT(Output, ctx) {
  auto& [m, mex, prop] = ctx;
  {
    using Writer = ActsExamples::ObjPropagationStepsWriter<Acts::detail::Step>;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "ObjPropagationStepsWriter")
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
                 mex, "RootPropagationStepsWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level") = Acts::Logging::INFO);

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("collection", &Writer::Config::collection)
        .def_readwrite("filePath", &Writer::Config::filePath)
        .def_readwrite("fileMode", &Writer::Config::fileMode)
        .def_readwrite("treeName", &Writer::Config::treeName);
  }

  {
    using Writer = ActsExamples::CsvParticleWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvParticleWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level") = Acts::Logging::INFO);

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("inputParticles", &Writer::Config::inputParticles)
        .def_readwrite("outputDir", &Writer::Config::outputDir)
        .def_readwrite("outputStem", &Writer::Config::outputStem)
        .def_readwrite("outputPrecision", &Writer::Config::outputPrecision);
  }

  {
    using Writer = ActsExamples::RootParticleWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootParticleWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level") = Acts::Logging::INFO);

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("inputParticles", &Writer::Config::inputParticles)
        .def_readwrite("filePath", &Writer::Config::filePath)
        .def_readwrite("fileMode", &Writer::Config::fileMode)
        .def_readwrite("treeName", &Writer::Config::treeName);
  }
}
}  // namespace