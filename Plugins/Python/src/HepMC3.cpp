#include "Acts/Plugins/Python/Utilities.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {
void addHepMC3(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  // {
  //   using Alg = ActsExamples::TrackParametersPrinter;

  //   auto alg =
  //       py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
  //           mex, "TrackParametersPrinter")
  //           .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
  //                py::arg("config"), py::arg("level"));

  //   py::class_<Alg::Config>(alg, "Config")
  //       .def(py::init<>())
  //       .def_readwrite("inputTrackParameters",
  //                      &Alg::Config::inputTrackParameters);
  // }
}
}  // namespace Acts::Python