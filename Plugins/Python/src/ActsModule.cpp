#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void addExamplesDetector(py::module_& m);
void addMaterial(py::module_& m);
void addGeometry(py::module_& m);

PYBIND11_MODULE(_acts, m) {
  m.doc() = "Acts";

  auto logging = m.def_submodule("logging", "");
  py::enum_<Acts::Logging::Level>(logging, "Level")
      .value("VERBOSE", Acts::Logging::VERBOSE)
      .value("DEBUG", Acts::Logging::DEBUG)
      .value("INFO", Acts::Logging::INFO)
      .value("WARNING", Acts::Logging::WARNING)
      .value("ERROR", Acts::Logging::ERROR)
      .value("FATAL", Acts::Logging::FATAL)
      .export_values();

  using ActsExamples::Sequencer;
  using Config = Sequencer::Config;
  auto sequencer =
      py::class_<Sequencer>(m, "Sequencer")
          .def(py::init<const Config&>())
          .def("run", &Sequencer::run)
          .def("addContextDecorator", &Sequencer::addContextDecorator);

  py::class_<Config>(sequencer, "Config")
      .def(py::init<>())
      .def_readwrite("skip", &Config::skip)
      .def_readwrite("events", &Config::events)
      .def_readwrite("logLevel", &Config::logLevel)
      .def_readwrite("numThreads", &Config::numThreads)
      .def_readwrite("outputDir", &Config::outputDir);

  using ActsExamples::RandomNumbers;
  auto randomNumbers = py::class_<RandomNumbers>(m, "RandomNumbers");

  py::class_<RandomNumbers::Config>(randomNumbers, "Config")
      .def(py::init<>())
      .def_readwrite("seed", &RandomNumbers::Config::seed);

  randomNumbers.def(py::init<const RandomNumbers::Config&>());

  addExamplesDetector(m);
  addMaterial(m);
  addGeometry(m);
}