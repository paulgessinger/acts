#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {
void addGeometry(Context& ctx) {
  auto m = ctx.get("main");

  py::class_<Acts::Surface, std::shared_ptr<Acts::Surface>>(m, "Surfaces");

  {
    py::class_<Acts::TrackingGeometry, std::shared_ptr<Acts::TrackingGeometry>>(
        m, "TrackingGeometry")
        .def("visitSurfaces", &Acts::TrackingGeometry::visitSurfaces);
  }
}

}  // namespace Acts::Python