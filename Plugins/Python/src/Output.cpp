#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

namespace {
ACTS_PYTHON_COMPONENT(Output, ctx) {
  auto& [m, mex, prop] = ctx;
  {
    using Writer = ActsExamples::ObjPropagationStepsWriter<Acts::detail::Step>;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "ObjPropagationStepsWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

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
                      py::arg("cfg"), py::arg("level"));

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
                      py::arg("cfg"), py::arg("level"));

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
                      py::arg("cfg"), py::arg("level"));

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("inputParticles", &Writer::Config::inputParticles)
        .def_readwrite("filePath", &Writer::Config::filePath)
        .def_readwrite("fileMode", &Writer::Config::fileMode)
        .def_readwrite("treeName", &Writer::Config::treeName);
  }

  {
    using Writer = ActsExamples::TrackFinderPerformanceWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "TrackFinderPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("inputProtoTracks", &Writer::Config::inputProtoTracks)
        .def_readwrite("inputMeasurementParticlesMap",
                       &Writer::Config::inputMeasurementParticlesMap)
        .def_readwrite("inputParticles", &Writer::Config::inputParticles)
        .def_readwrite("outputDir", &Writer::Config::outputDir)
        .def_readwrite("outputFilename", &Writer::Config::outputFilename);
  }

  {
    using Writer = ActsExamples::SeedingPerformanceWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "SeedingPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    py::class_<Writer::Config>(w, "Config")
        .def(py::init<>())
        .def_readwrite("inputProtoTracks", &Writer::Config::inputProtoTracks)
        .def_readwrite("inputMeasurementParticlesMap",
                       &Writer::Config::inputMeasurementParticlesMap)
        .def_readwrite("inputParticles", &Writer::Config::inputParticles)
        .def_readwrite("outputDir", &Writer::Config::outputDir)
        .def_readwrite("outputFilename", &Writer::Config::outputFilename)
        .def_readwrite("effPlotToolConfig", &Writer::Config::effPlotToolConfig)
        .def_readwrite("duplicationPlotToolConfig",
                       &Writer::Config::duplicationPlotToolConfig);
  }

  {
    using Writer = ActsExamples::RootTrackParameterWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootTrackParameterWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputFilename);
    ACTS_PYTHON_MEMBER(outputTreename);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootMaterialTrackWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootMaterialTrackWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(collection);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(recalculateTotals);
    ACTS_PYTHON_MEMBER(prePostStep);
    ACTS_PYTHON_MEMBER(storeSurface);
    ACTS_PYTHON_MEMBER(storeVolume);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace