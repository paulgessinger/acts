#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {
void addInput(Context& ctx) {
  auto mex = ctx.get("examples");

  {
    using Reader = ActsExamples::RootParticleReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "RootParticleReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(particleCollection);
    ACTS_PYTHON_MEMBER(vertexPrimaryCollection);
    ACTS_PYTHON_MEMBER(vertexSecondaryCollection);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(inputFile);
    ACTS_PYTHON_MEMBER(inputDir);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python