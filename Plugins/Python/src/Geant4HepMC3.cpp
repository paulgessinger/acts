#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Geant4HepMC/EventRecording.hpp"

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {
void addGeant4HepMC3(Context& ctx) {
  auto m = ctx.get("geant4");

  auto h3 = m.def_submodule("hepmc3");

  {
    using Alg = EventRecording;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            h3, "EventRecording")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputHepMcTracks);
    ACTS_PYTHON_MEMBER(detectorConstructionFactory);
    ACTS_PYTHON_MEMBER(seed1);
    ACTS_PYTHON_MEMBER(seed2);
    ACTS_PYTHON_MEMBER(processesCombine);
    ACTS_PYTHON_MEMBER(processSelect);
    ACTS_PYTHON_MEMBER(processesReject);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python