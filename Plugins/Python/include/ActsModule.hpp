#pragma once

#include <functional>
#include <memory>
#include <vector>

namespace pybind11 {
class module_;
}

namespace Acts {

struct PythonContext {
  pybind11::module_& main;
  pybind11::module_& examples;
  pybind11::module_& propagator;
};

class ComponentRegistry {
 public:
  static ComponentRegistry& instance() {
    static ComponentRegistry reg;
    return reg;
  }

  void add(std::function<void(PythonContext&)> comp) {
    m_components.push_back(std::move(comp));
  }

  void apply(PythonContext& ctx) {
    for (auto& comp : m_components) {
      comp(ctx);
    }
  }

 private:
  std::vector<std::function<void(PythonContext&)>> m_components;
};

namespace detail {

struct ModuleRegisterer {
  ModuleRegisterer(std::function<void(PythonContext&)> comp) {
    auto& reg = ComponentRegistry::instance();
    reg.add(comp);
  }
};
}  // namespace detail

#define ACTS_PYTHON_COMPONENT(name, ctxvar)                                  \
  void _acts_init_##name(Acts::PythonContext& ctxvar);                       \
  static Acts::detail::ModuleRegisterer _acts_reg##name{&_acts_init_##name}; \
  void _acts_init_##name(Acts::PythonContext& ctxvar)

}  // namespace Acts