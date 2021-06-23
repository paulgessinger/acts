#include "ActsModule.hpp"

#include "Acts/ActsVersion.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

constexpr std::tuple<unsigned int, unsigned int, unsigned int> version{
    Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

using namespace ActsExamples;

namespace {
class PyIAlgorithm : public IAlgorithm {
 public:
  using IAlgorithm::IAlgorithm;

  std::string name() const override {
    py::gil_scoped_acquire acquire{};
    PYBIND11_OVERRIDE_PURE(std::string, IAlgorithm, name);
  }

  ProcessCode execute(const AlgorithmContext& ctx) const override {
    py::gil_scoped_acquire acquire{};
    PYBIND11_OVERRIDE_PURE(ProcessCode, IAlgorithm, execute, ctx);
  }
};
}  // namespace

PYBIND11_MODULE(_acts, m) {
  m.doc() = "Acts";

  m.attr("__version__") =
      std::tuple{Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

  {
    py::class_<Acts::Logger, std::unique_ptr<Acts::Logger>>(m, "Logger");
    m.def("getDefaultLogger", &Acts::getDefaultLogger);
  }

  auto mex = m.def_submodule("_examples");

  py::class_<ActsExamples::IWriter, std::shared_ptr<ActsExamples::IWriter>>(
      mex, "IWriter");

  py::class_<ActsExamples::IReader, std::shared_ptr<ActsExamples::IReader>>(
      mex, "IReader");

  // py::enum_<ProcessCode>(mex, "ProcessCode")
  //     .value("SUCCESS", ProcessCode::SUCCESS)
  //     .value("ABORT", ProcessCode::ABORT)
  //     .value("END", ProcessCode::END);

  // py::class_<AlgorithmContext>(mex, "AlgorithmContext")
  //     .def_readonly("algorithmNumber", &AlgorithmContext::algorithmNumber)
  //     .def_readonly("eventNumber", &AlgorithmContext::eventNumber)
  //     // .def_readonly("eventStore", &AlgorithmContext::eventStore)
  //     // .def_readonly("magFieldContext", &AlgorithmContext::magFieldContext)
  //     // .def_readonly("calibContext", &AlgorithmContext::calibContext)
  //     ;

  auto iAlgorithm =
      py::class_<ActsExamples::IAlgorithm, PyIAlgorithm,
                 std::shared_ptr<ActsExamples::IAlgorithm>>(mex, "IAlgorithm")
          .def(py::init<>())
      // .def("execute", &IAlgorithm::execute)
      // .def("name", &ActsExamples::IAlgorithm::name)
      ;

  auto bareAlgorithm =
      py::class_<ActsExamples::BareAlgorithm, ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::BareAlgorithm>>(mex,
                                                               "BareAlgorithm");

  using ActsExamples::Sequencer;
  using Config = Sequencer::Config;
  auto sequencer =
      py::class_<Sequencer>(mex, "Sequencer")
          .def(py::init([](Config cfg) {
            cfg.iterationCallback = []() {
              if (PyErr_CheckSignals() != 0) {
                throw py::error_already_set{};
              }
            };
            return Sequencer{cfg};
          }))
          .def("run", &Sequencer::run)
          .def("addContextDecorator", &Sequencer::addContextDecorator)
          .def("addAlgorithm", &Sequencer::addAlgorithm)
          .def("addAlgorithms",
               [](Sequencer& s, py::args algs) {
                 for (auto& alg : algs) {
                   s.addAlgorithm(
                       alg.cast<std::shared_ptr<ActsExamples::IAlgorithm>>());
                 }
               })
          .def("addReader", &Sequencer::addReader)
          .def("addWriter", &Sequencer::addWriter);

  py::class_<Config>(sequencer, "Config")
      .def(py::init<>())
      .def_readwrite("skip", &Config::skip)
      .def_readwrite("events", &Config::events)
      .def_readwrite("logLevel", &Config::logLevel)
      .def_readwrite("numThreads", &Config::numThreads)
      .def_readwrite("outputDir", &Config::outputDir);

  using ActsExamples::RandomNumbers;
  auto randomNumbers =
      py::class_<RandomNumbers, std::shared_ptr<RandomNumbers>>(mex,
                                                                "RandomNumbers")
          .def(py::init<const RandomNumbers::Config&>());

  py::class_<ActsExamples::RandomEngine>(mex, "RandomEngine").def(py::init<>());

  py::class_<RandomNumbers::Config>(randomNumbers, "Config")
      .def(py::init<>())
      .def_readwrite("seed", &RandomNumbers::Config::seed);

  auto prop = m.def_submodule("_propagator");

  Acts::PythonContext ctx{m, mex, prop};

  Acts::ComponentRegistry::instance().apply(ctx);
}