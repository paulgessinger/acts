#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void addExamplesDetector(py::module_& m);
void addMaterial(py::module_& m, py::module_& mex);
void addGeometry(py::module_& m);
void addExamplesAlgorithms(py::module_& m, py::module_& prop);
void addPropagation(py::module_& m, py::module_& prop, py::module_& mex);
void addMagneticField(py::module_& m);
void addOutput(py::module_& m);
void addGenerators(py::module_& mex);

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

void addAlgebra(py::module_& m) {
  py::class_<Acts::Vector3>(m, "Vector3")
      .def(py::init<double, double, double>())
      .def(py::init([](std::array<double, 3> a) {
        Acts::Vector3 v;
        v << a[0], a[1], a[2];
        return v;
      }));

  py::class_<Acts::Vector4>(m, "Vector4")
      .def(py::init<double, double, double, double>())
      .def(py::init([](std::array<double, 4> a) {
        Acts::Vector4 v;
        v << a[0], a[1], a[2], a[3];
        return v;
      }));
}

PYBIND11_MODULE(_acts, m) {
  m.doc() = "Acts";

  auto mex = m.def_submodule("_examples");

  auto logging = m.def_submodule("logging", "");
  py::enum_<Acts::Logging::Level>(logging, "Level")
      .value("VERBOSE", Acts::Logging::VERBOSE)
      .value("DEBUG", Acts::Logging::DEBUG)
      .value("INFO", Acts::Logging::INFO)
      .value("WARNING", Acts::Logging::WARNING)
      .value("ERROR", Acts::Logging::ERROR)
      .value("FATAL", Acts::Logging::FATAL)
      .export_values();

  py::enum_<Acts::PdgParticle>(m, "PdgParticle")
      .value("eInvalid", Acts::PdgParticle::eInvalid)
      .value("eElectron", Acts::PdgParticle::eElectron)
      .value("eAntiElectron", Acts::PdgParticle::eAntiElectron)
      .value("ePositron", Acts::PdgParticle::ePositron)
      .value("eMuon", Acts::PdgParticle::eMuon)
      .value("eAntiMuon", Acts::PdgParticle::eAntiMuon)
      .value("eTau", Acts::PdgParticle::eTau)
      .value("eAntiTau", Acts::PdgParticle::eAntiTau)
      .value("eGamma", Acts::PdgParticle::eGamma)
      .value("ePionZero", Acts::PdgParticle::ePionZero)
      .value("ePionPlus", Acts::PdgParticle::ePionPlus)
      .value("ePionMinus", Acts::PdgParticle::ePionMinus)
      .value("eNeutron", Acts::PdgParticle::eNeutron)
      .value("eAntiNeutron", Acts::PdgParticle::eAntiNeutron)
      .value("eProton", Acts::PdgParticle::eProton)
      .value("eAntiProton", Acts::PdgParticle::eAntiProton)
      .export_values();

  py::class_<ActsExamples::IWriter, std::shared_ptr<ActsExamples::IWriter>>(
      mex, "IWriter");

  py::class_<ActsExamples::IReader, std::shared_ptr<ActsExamples::IReader>>(
      mex, "IReader");

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

  addUnits(m);
  addAlgebra(m);
  addGeometry(m);
  addMagneticField(m);
  addExamplesDetector(mex);
  addMaterial(m, mex);
  addExamplesAlgorithms(mex, prop);
  addPropagation(m, prop, mex);
  addOutput(mex);
  addGenerators(mex);
}