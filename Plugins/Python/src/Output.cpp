#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootBFieldWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {
void addOutput(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");
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
        .def_readwrite("filePath", &Writer::Config::filePath);
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

  {
    using Writer = ActsExamples::RootBFieldWriter;
    auto w =
        py::class_<Writer>(mex, "RootBFieldWriter")
            .def_static(
                "run",
                [](const Writer::Config& config, Acts::Logging::Level level) {
                  Writer::run(config, Acts::getDefaultLogger("RootBFieldWriter",
                                                             level));
                },
                py::arg("config"), py::arg("level"));

    py::enum_<Writer::GridType>(w, "GridType")
        .value("rz", Writer::GridType::rz)
        .value("xyz", Writer::GridType::xyz);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileName);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(bField);
    ACTS_PYTHON_MEMBER(gridType);
    ACTS_PYTHON_MEMBER(rBounds);
    ACTS_PYTHON_MEMBER(zBounds);
    ACTS_PYTHON_MEMBER(rBins);
    ACTS_PYTHON_MEMBER(zBins);
    ACTS_PYTHON_MEMBER(phiBins);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootMeasurementWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootMeasurementWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    c.def("addBoundIndicesFromDigiConfig",
          [](Writer::Config& self, const DigitizationConfig& digiCfg) {
            self.boundIndices =
                Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(
                    digiCfg.getBoundIndices());
          });

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(boundIndices);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CsvMeasurementWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvMeasurementWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python