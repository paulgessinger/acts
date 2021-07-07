#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {
void addJson(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<JsonMaterialDecorator, Acts::IMaterialDecorator,
               std::shared_ptr<JsonMaterialDecorator>>(m,
                                                       "JsonMaterialDecorator")
        .def(py::init<const MaterialMapJsonConverter::Config&,
                      const std::string&, Acts::Logging::Level, bool, bool>(),
             py::arg("rConfig"), py::arg("jFileName"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true);
  }

  {
    auto cls =
        py::class_<MaterialMapJsonConverter>(m, "MaterialMapJsonConverter")
            .def(py::init<const MaterialMapJsonConverter::Config&,
                          Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<MaterialMapJsonConverter::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, MaterialMapJsonConverter::Config);
    ACTS_PYTHON_MEMBER(context);
    ACTS_PYTHON_MEMBER(processSensitives);
    ACTS_PYTHON_MEMBER(processApproaches);
    ACTS_PYTHON_MEMBER(processRepresenting);
    ACTS_PYTHON_MEMBER(processBoundaries);
    ACTS_PYTHON_MEMBER(processVolumes);
    ACTS_PYTHON_MEMBER(processDenseVolumes);
    ACTS_PYTHON_MEMBER(processNonMaterial);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python