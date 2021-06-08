#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

void addMagneticField(py::module_& m) {
  py::class_<Acts::MagneticFieldProvider,
             std::shared_ptr<Acts::MagneticFieldProvider>>(
      m, "MagneticFieldProvider");

  py::class_<Acts::Vector3>(m, "Vector3")
      .def(py::init<double, double, double>());

  py::class_<Acts::ConstantBField, Acts::MagneticFieldProvider,
             std::shared_ptr<Acts::ConstantBField>>(m, "ConstantBField")
      .def(py::init<Acts::Vector3>());
}