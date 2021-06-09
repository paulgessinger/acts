#include "Acts/Material/IMaterialDecorator.hpp"
#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

#define PY_MEMBER(obj, t, name) obj.def_readwrite(#name, &t::name)

void addMaterial(py::module_& m, py::module_& mex) {
  {
    py::class_<Acts::IMaterialDecorator,
               std::shared_ptr<Acts::IMaterialDecorator>>(m,
                                                          "IMaterialDecorator");
  }

  {
    auto rmd =
        py::class_<RootMaterialDecorator,
                   std::shared_ptr<RootMaterialDecorator>,
                   Acts::IMaterialDecorator>(mex, "RootMaterialDecorator")
            .def(py::init<RootMaterialDecorator::Config>());

    using Config = RootMaterialDecorator::Config;
    auto c =
        py::class_<Config>(rmd, "Config")
            .def(py::init<const std::string&, Acts::Logging::Level>(),
                 "lname"_a = "MaterialReader", "lvl"_a = Acts::Logging::INFO);

    PY_MEMBER(c, Config, folderNameBase);
    PY_MEMBER(c, Config, voltag);
    PY_MEMBER(c, Config, boutag);
    PY_MEMBER(c, Config, laytag);
    PY_MEMBER(c, Config, apptag);
    PY_MEMBER(c, Config, sentag);
    PY_MEMBER(c, Config, ntag);
    PY_MEMBER(c, Config, vtag);
    PY_MEMBER(c, Config, otag);
    PY_MEMBER(c, Config, mintag);
    PY_MEMBER(c, Config, maxtag);
    PY_MEMBER(c, Config, ttag);
    PY_MEMBER(c, Config, x0tag);
    PY_MEMBER(c, Config, l0tag);
    PY_MEMBER(c, Config, atag);
    PY_MEMBER(c, Config, ztag);
    PY_MEMBER(c, Config, rhotag);
    PY_MEMBER(c, Config, fileName);
    PY_MEMBER(c, Config, logger);
    PY_MEMBER(c, Config, name);
  }
}