// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicies.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TypeTag.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

#include <boost/core/demangle.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {

struct AnyNavigationPolicyFactory : public Acts::NavigationPolicyFactory {
  virtual std::unique_ptr<AnyNavigationPolicyFactory> add(
      TypeTag<TryAllPortalNavigationPolicy> /*type*/) = 0;

  virtual std::unique_ptr<AnyNavigationPolicyFactory> add(
      TypeTag<TryAllSurfaceNavigationPolicy> /*type*/) = 0;

  virtual std::unique_ptr<AnyNavigationPolicyFactory> add(
      TypeTag<SurfaceArrayNavigationPolicy> /*type*/,
      SurfaceArrayNavigationPolicy::Config config) = 0;
};

template <typename Factory = detail::NavigationPolicyFactoryImpl<>,
          typename... Policies>
struct NavigationPolicyFactoryT : public AnyNavigationPolicyFactory {
  NavigationPolicyFactoryT(Factory impl)
    requires(sizeof...(Policies) > 0)
      : m_impl(std::move(impl)) {}

  NavigationPolicyFactoryT()
    requires(sizeof...(Policies) == 0)
      : m_impl{} {}

  std::unique_ptr<AnyNavigationPolicyFactory> add(
      TypeTag<TryAllPortalNavigationPolicy> /*type*/) override {
    return add<TryAllPortalNavigationPolicy>();
  }

  std::unique_ptr<AnyNavigationPolicyFactory> add(
      TypeTag<TryAllSurfaceNavigationPolicy> /*type*/) override {
    return add<TryAllSurfaceNavigationPolicy>();
  }

  std::unique_ptr<AnyNavigationPolicyFactory> add(
      TypeTag<SurfaceArrayNavigationPolicy> /*type*/,
      SurfaceArrayNavigationPolicy::Config config) override {
    return add<SurfaceArrayNavigationPolicy>(std::move(config));
  }

  std::unique_ptr<INavigationPolicy> build(
      const GeometryContext& gctx, const TrackingVolume& volume,
      const Logger& logger) const override {
    if constexpr (sizeof...(Policies) > 0) {
      return m_impl.build(gctx, volume, logger);
    } else {
      throw std::runtime_error("No policies added to the factory");
    }
  }

 private:
  template <typename T, typename... Args>
  std::unique_ptr<AnyNavigationPolicyFactory> add(Args&&... args) {
    if constexpr (!((std::is_same_v<T, Policies> || ...))) {
      auto impl = m_impl.template add<T>(std::forward<Args>(args)...);
      return std::make_unique<
          NavigationPolicyFactoryT<decltype(impl), Policies..., T>>(
          std::move(impl));
    } else {
      throw std::invalid_argument("Policy already added to the factory");
    }
  }

  Factory m_impl;
};

class NavigationPolicyFactory : public Acts::NavigationPolicyFactory {
 public:
  // This overload is for all the navigation policies that don't have extra
  // arguments
  NavigationPolicyFactory& addNoArguments(const py::object& cls) {
    auto m = py::module_::import("acts");
    if (py::object o = m.attr("TryAllPortalNavigationPolicy"); cls.is(o)) {
      m_impl = m_impl->add(Type<TryAllPortalNavigationPolicy>);
    } else if (o = m.attr("TryAllSurfaceNavigationPolicy"); cls.is(o)) {
      m_impl = m_impl->add(Type<TryAllSurfaceNavigationPolicy>);
    } else {
    }
    return *this;
  }

  NavigationPolicyFactory& addSurfaceArray(
      const py::object& /*cls*/,
      const SurfaceArrayNavigationPolicy::Config& config) {
    m_impl = m_impl->add(Type<SurfaceArrayNavigationPolicy>, config);
    return *this;
  }

  std::unique_ptr<INavigationPolicy> build(
      const GeometryContext& gctx, const TrackingVolume& volume,
      const Logger& logger) const override {
    return m_impl->build(gctx, volume, logger);
  }

 private:
  std::unique_ptr<AnyNavigationPolicyFactory> m_impl =
      std::make_unique<NavigationPolicyFactoryT<>>();
};

void addNavigation(Context& ctx) {
  auto m = ctx.get("main");

  py::class_<Acts::NavigationPolicyFactory,
             std::shared_ptr<Acts::NavigationPolicyFactory>>(
      m, "_NavigationPolicyFactory");

  py::class_<TryAllPortalNavigationPolicy>(m, "TryAllPortalNavigationPolicy");
  py::class_<TryAllSurfaceNavigationPolicy>(m, "TryAllSurfaceNavigationPolicy");

  py::class_<NavigationPolicyFactory, Acts::NavigationPolicyFactory,
             std::shared_ptr<NavigationPolicyFactory>>(
      m, "NavigationPolicyFactory")
      .def(py::init<>())
      .def("add", &NavigationPolicyFactory::addNoArguments)
      .def("add", &NavigationPolicyFactory::addSurfaceArray)
      .def("_buildTest", [](NavigationPolicyFactory& self) {
        auto vol1 = std::make_shared<TrackingVolume>(
            Transform3::Identity(),
            std::make_shared<CylinderVolumeBounds>(30, 40, 100));
        vol1->setVolumeName("TestVolume");

        std::unique_ptr<INavigationPolicy> result =
            self.build(GeometryContext{}, *vol1,
                       *getDefaultLogger("Test", Logging::VERBOSE));
      });

  {
    auto saPolicy = py::class_<SurfaceArrayNavigationPolicy>(
        m, "SurfaceArrayNavigationPolicy");

    using LayerType = SurfaceArrayNavigationPolicy::LayerType;
    py::enum_<LayerType>(saPolicy, "LayerType")
        .value("Cylinder", LayerType::Cylinder)
        .value("Disc", LayerType::Disc)
        .value("Plane", LayerType::Plane);

    using Config = SurfaceArrayNavigationPolicy::Config;
    auto c = py::class_<Config>(saPolicy, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(surfaceArrayConfig);
    ACTS_PYTHON_MEMBER(layerType);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python