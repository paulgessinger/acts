#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
// using namespace ActsExamples;

using EigenPropagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;

void addExamplesAlgorithms(py::module_& m) {
  py::class_<EigenPropagator, std::shared_ptr<EigenPropagator>>(
      m, "EigenPropagator")
      //   .def(py::init<Acts::EigenStepper<>, Acts::Navigator>());
      .def(py::init([](std::shared_ptr<Acts::EigenStepper<>> stepper,
                       std::shared_ptr<Acts::Navigator> navigator) {
        return std::make_shared<EigenPropagator>(*stepper, *navigator);
      }));

  py::class_<Acts::EigenStepper<>, std::shared_ptr<Acts::EigenStepper<>>>(
      m, "EigenStepper")
      .def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());

  py::class_<Acts::Navigator, std::shared_ptr<Acts::Navigator>>(m, "Navigator")
      .def(py::init<std::shared_ptr<const Acts::TrackingGeometry>>())
      .def_readwrite("resolveMaterial", &Acts::Navigator::resolveMaterial)
      .def_readwrite("resolvePassive", &Acts::Navigator::resolvePassive)
      .def_readwrite("resolveSensitive", &Acts::Navigator::resolveSensitive);

  auto iAlgorithm =
      py::class_<ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::IAlgorithm>>(m, "IAlgorithm");

  auto bareAlgorithm =
      py::class_<ActsExamples::BareAlgorithm, ActsExamples::IAlgorithm,
                 std::shared_ptr<ActsExamples::BareAlgorithm>>(m,
                                                               "BareAlgorithm");

  {
    using Algorithm = ActsExamples::PropagationAlgorithm<EigenPropagator>;
    using Config = Algorithm::Config;
    auto alg =
        py::class_<Algorithm, ActsExamples::BareAlgorithm,
                   std::shared_ptr<Algorithm>>(m, "EigenPropagationAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>());

    py::class_<Config>(alg, "Config")
        .def(py::init<EigenPropagator>())
        .def(py::init<>())
        .def_readwrite("randomNumberSvc", &Config::randomNumberSvc)
        .def_readwrite("propagator", &Config::propagator)
        .def_readwrite("ntests", &Config::ntests);
  }
}