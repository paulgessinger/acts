#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/HepMC/HepMCProcessExtractor.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {
void addHepMC3(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  auto hepmc3 = mex.def_submodule("_hepmc3");

  {
    using Alg = ActsExamples::HepMCProcessExtractor;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            hepmc3, "HepMCProcessExtractor")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(inputEvents);
    ACTS_PYTHON_MEMBER(outputSimulationProcesses);
    ACTS_PYTHON_MEMBER(extractionProcess);
    ACTS_PYTHON_MEMBER(absPdgMin);
    ACTS_PYTHON_MEMBER(absPdgMax);
    ACTS_PYTHON_MEMBER(pMin);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python