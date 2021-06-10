#include "ActsExamples/MagneticField/MagneticField.hpp"

#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

ACTS_PYTHON_COMPONENT(MagneticField, ctx) {
  auto& [m, mex, prop] = ctx;
  py::class_<Acts::MagneticFieldProvider,
             std::shared_ptr<Acts::MagneticFieldProvider>>(
      m, "MagneticFieldProvider");

  m.def("solenoidFieldMap", &Acts::solenoidFieldMap, py::arg("rlim"),
        py::arg("zlim"), py::arg("nbins"), py::arg("field"));

  py::class_<Acts::ConstantBField, Acts::MagneticFieldProvider,
             std::shared_ptr<Acts::ConstantBField>>(m, "ConstantBField")
      .def(py::init<Acts::Vector3>());

  py::class_<ActsExamples::detail::InterpolatedMagneticField2,
             Acts::MagneticFieldProvider,
             std::shared_ptr<ActsExamples::detail::InterpolatedMagneticField2>>(
      m, "InterpolatedMagneticField2");

  py::class_<ActsExamples::detail::InterpolatedMagneticField3,
             Acts::MagneticFieldProvider,
             std::shared_ptr<ActsExamples::detail::InterpolatedMagneticField3>>(
      m, "InterpolatedMagneticField3");

  py::class_<Acts::NullBField, Acts::MagneticFieldProvider,
             std::shared_ptr<Acts::NullBField>>(m, "NullBField")
      .def(py::init<>());

  {
    using Config = Acts::SolenoidBField::Config;

    auto sol =
        py::class_<Acts::SolenoidBField, Acts::MagneticFieldProvider,
                   std::shared_ptr<Acts::SolenoidBField>>(m, "SolenoidBField")
            .def(py::init<Config>())
            .def(py::init([](double radius, double length, size_t nCoils,
                             double bMagCenter) {
                   return Acts::SolenoidBField{
                       Config{radius, length, nCoils, bMagCenter}};
                 }),
                 py::arg("radius"), py::arg("length"), py::arg("nCoils"),
                 py::arg("bMagCenter"));

    py::class_<Config>(sol, "Config")
        .def(py::init<>())
        .def_readwrite("radius", &Config::radius)
        .def_readwrite("length", &Config::length)
        .def_readwrite("nCoils", &Config::nCoils)
        .def_readwrite("bMagCenter", &Config::bMagCenter);
  }
}

}  // namespace