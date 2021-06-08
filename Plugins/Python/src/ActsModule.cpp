#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void addExamplesDetector(py::module_& m);
void addMaterial(py::module_& m);
void addGeometry(py::module_& m);
void addExamplesAlgorithms(py::module_& m);
void addMagneticField(py::module_& m);
void addOutput(py::module_& m);

void addUnits(py::module_& m) {
  auto u = m.def_submodule("UnitConstants");
  using namespace Acts::UnitConstants;

#define UNIT(x) u.attr(#x) = x;

  UNIT(fm)
  UNIT(pm)
  UNIT(um)
  UNIT(nm)
  UNIT(mm)
  UNIT(cm)
  UNIT(m)
  UNIT(km)
  UNIT(mm2)
  UNIT(cm2)
  UNIT(m2)
  UNIT(mm3)
  UNIT(cm3)
  UNIT(m3)
  UNIT(s)
  UNIT(fs)
  UNIT(ps)
  UNIT(ns)
  UNIT(us)
  UNIT(ms)
  UNIT(min)
  UNIT(h)
  UNIT(mrad)
  UNIT(rad)
  UNIT(degree)
  UNIT(eV)
  UNIT(keV)
  UNIT(MeV)
  UNIT(GeV)
  UNIT(TeV)
  UNIT(J)
  UNIT(u)
  UNIT(g)
  UNIT(kg)
  UNIT(e)
  UNIT(C)
  UNIT(T)
  UNIT(Gauss)
  UNIT(kGauss)
  UNIT(mol)

#undef UNIT
}

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

  py::class_<ActsExamples::IWriter, std::shared_ptr<ActsExamples::IWriter>>(
      m, "IWriter");

  using ActsExamples::Sequencer;
  using Config = Sequencer::Config;
  auto sequencer =
      py::class_<Sequencer>(m, "Sequencer")
          .def(py::init<const Config&>())
          .def("run", &Sequencer::run)
          .def("addContextDecorator", &Sequencer::addContextDecorator)
          .def("addAlgorithm", &Sequencer::addAlgorithm)
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
      py::class_<RandomNumbers, std::shared_ptr<RandomNumbers>>(
          m, "RandomNumbers");

  py::class_<RandomNumbers::Config>(randomNumbers, "Config")
      .def(py::init<>())
      .def_readwrite("seed", &RandomNumbers::Config::seed);

  randomNumbers.def(py::init<const RandomNumbers::Config&>());

  addUnits(m);
  addGeometry(m);
  addMagneticField(m);
  addExamplesDetector(m);
  addMaterial(m);
  addExamplesAlgorithms(m);
  addOutput(m);
}