#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

// #define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

namespace {

ACTS_PYTHON_COMPONENT(Geometry, ctx) {
  auto& [m, mex, prop] = ctx;
  {
    py::class_<Acts::TrackingGeometry, std::shared_ptr<Acts::TrackingGeometry>>(
        m, "TrackingGeometry");
  }
}

}  // namespace