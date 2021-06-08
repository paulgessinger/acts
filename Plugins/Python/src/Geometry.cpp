#include "Acts/Geometry/TrackingGeometry.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

// #define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

void addGeometry(py::module_& m) {
  {
    py::class_<Acts::TrackingGeometry, std::shared_ptr<Acts::TrackingGeometry>>(
        m, "TrackingGeometry");
  }
}