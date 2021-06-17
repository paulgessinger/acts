#pragma once

#include <string>
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
  using Function = void (*)(PythonContext&);

  static ComponentRegistry& instance() {
    static ComponentRegistry reg;
    return reg;
  }

  void add(Function comp) { m_components.push_back(std::move(comp)); }

  void apply(PythonContext& ctx) {
    for (auto& comp : m_components) {
      comp(ctx);
    }
  }

 private:
  std::vector<Function> m_components;
};

namespace detail {

struct ModuleRegisterer {
  ModuleRegisterer(ComponentRegistry::Function comp) {
    auto& reg = ComponentRegistry::instance();
    reg.add(comp);
  }
};
}  // namespace detail

#define ACTS_PYTHON_COMPONENT(name, ctxvar)                                   \
  void _acts_init_##name(Acts::PythonContext& ctxvar);                        \
  static Acts::detail::ModuleRegisterer _acts_reg_##name{&_acts_init_##name}; \
  void _acts_init_##name(Acts::PythonContext& ctxvar)

#define ACTS_PYTHON_MEMBER(name) \
  _binding_instance.def_readwrite(#name, &_struct_type::name)

#define ACTS_PYTHON_STRUCT_BEGIN(obj, cls) \
  {                                        \
    auto& _binding_instance = obj;         \
    using _struct_type = cls;              \
    do {                                   \
    } while (0)

#define ACTS_PYTHON_STRUCT_END() \
  }                              \
  do {                           \
  } while (0)

template <typename T, typename Ur, typename Ut>
void pythonRangeProperty(T& obj, const std::string& name, Ur Ut::*begin,
                         Ur Ut::*end) {  // Ut::_;
  obj.def_property(
      name.c_str(),
      [=](Ut& self) {
        return std::pair{self.*begin, self.*end};
      },
      [=](Ut& self, std::pair<Ur, Ur> p) {
        self.*begin = p.first;
        self.*end = p.second;
      });
}

}  // namespace Acts