#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsModule.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

namespace {
ACTS_PYTHON_COMPONENT(Input, ctx) {
  auto& [m, mex, prop] = ctx;

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
    PY_MEMBER(c, Config, particleCollection);
    PY_MEMBER(c, Config, vertexPrimaryCollection);
    PY_MEMBER(c, Config, vertexSecondaryCollection);
    PY_MEMBER(c, Config, treeName);
    PY_MEMBER(c, Config, inputFile);
    PY_MEMBER(c, Config, inputDir);
  }
}
}  // namespace