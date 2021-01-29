#include "ActsExamples/Framework/Sequencer.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) {
  return i + j + 1;
}

PYBIND11_MODULE(_acts, m) {
  m.doc() = "Acts";
  m.def("add", &add, "Add stuff");

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
  auto sequencer = py::class_<Sequencer>(m, "Sequencer");

  using Config = Sequencer::Config;
  py::class_<Config>(sequencer, "Config")
      .def(py::init<>())
      .def_readwrite("skip", &Config::skip)
      .def_readwrite("events", &Config::events)
      .def_readwrite("logLevel", &Config::logLevel)
      .def_readwrite("numThreads", &Config::numThreads)
      .def_readwrite("outputDir", &Config::outputDir);

  sequencer.def(py::init<const Config&>());

  
}