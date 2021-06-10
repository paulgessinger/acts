#include "ActsModule.hpp"

#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_acts, m) {
  m.doc() = "Acts";

  auto mex = m.def_submodule("_examples");

  py::class_<ActsExamples::IWriter, std::shared_ptr<ActsExamples::IWriter>>(
      mex, "IWriter");

  py::class_<ActsExamples::IReader, std::shared_ptr<ActsExamples::IReader>>(
      mex, "IReader");

  auto iAlgorithm =
      py::class_<ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::IAlgorithm>>(mex, "IAlgorithm");

  auto bareAlgorithm =
      py::class_<ActsExamples::BareAlgorithm, ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::BareAlgorithm>>(mex,
                                                               "BareAlgorithm");

  using ActsExamples::Sequencer;
  using Config = Sequencer::Config;
  auto sequencer =
      py::class_<Sequencer>(mex, "Sequencer")
          .def(py::init<const Config&>())
          .def("run", &Sequencer::run)
          .def("addContextDecorator", &Sequencer::addContextDecorator)
          .def("addAlgorithm", &Sequencer::addAlgorithm)
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