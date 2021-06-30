#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {
void addGeometry(Context& ctx) {
  auto m = ctx.get("main");

  {
    py::class_<Acts::TrackingGeometry, std::shared_ptr<Acts::TrackingGeometry>>(
        m, "TrackingGeometry");
  }
}

}  // namespace Acts::Python