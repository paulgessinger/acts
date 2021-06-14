#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Generators/Pythia8ProcessGenerator.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace {
ACTS_PYTHON_COMPONENT(Pythia8, ctx) {
  auto& [m, mex, prop] = ctx;

  using Gen = ActsExamples::Pythia8Generator;
  auto gen = py::class_<Gen, ActsExamples::EventGenerator::ParticlesGenerator,
                        std::shared_ptr<Gen>>(mex, "Pythia8Generator")
                 .def(py::init<const Gen::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

  py::class_<Gen::Config>(gen, "Config")
      .def(py::init<>())
      .def_readwrite("pdgBeam0", &Gen::Config::pdgBeam0)
      .def_readwrite("pdgBeam1", &Gen::Config::pdgBeam1)
      .def_readwrite("cmsEnergy", &Gen::Config::cmsEnergy)
      .def_readwrite("settings", &Gen::Config::settings);
}
}  // namespace