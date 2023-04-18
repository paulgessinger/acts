// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
namespace Acts {

template <typename surface_t>
class ConstSurfacePtrT;

template <typename surface_t>
class SurfacePtrT {
 public:
  using element_type = surface_t;

  SurfacePtrT(std::shared_ptr<surface_t> surface)
      : m_shared{std::move(surface)} {}
  SurfacePtrT(surface_t* surface) : m_shared{surface} {}
  SurfacePtrT() = default;
  SurfacePtrT(const SurfacePtrT&) = default;
  SurfacePtrT(SurfacePtrT&&) = default;

  SurfacePtrT& operator=(const SurfacePtrT& other) = default;
  SurfacePtrT& operator=(SurfacePtrT&& other) = default;

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  SurfacePtrT& operator=(const SurfacePtrT<T>& other) {
    m_shared = other.m_shared;
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  SurfacePtrT& operator=(SurfacePtrT<T>&& other) {
    m_shared = std::move(other.m_shared);
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  SurfacePtrT(const SurfacePtrT<T>& other) : m_shared{other.m_shared} {}

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  SurfacePtrT(SurfacePtrT<T>&& other) : m_shared{std::move(other.m_shared)} {}

  surface_t& operator*() { return *m_shared; }
  surface_t* operator->() { return m_shared.get(); }
  operator bool() const { return !!m_shared; }
  surface_t* get() { return m_shared.get(); }

  const surface_t& operator*() const { return *m_shared; }
  const surface_t* operator->() const { return m_shared.get(); }
  const surface_t* get() const { return m_shared.get(); }

  template <typename T>
  friend class ConstSurfacePtrT;
  template <typename T>
  friend class SurfacePtrT;

 private:
  std::shared_ptr<surface_t> m_shared;
};

template <typename surface_t>
class ConstSurfacePtrT {
 public:
  using element_type = surface_t;

  ConstSurfacePtrT(surface_t* surface) : m_shared{surface} {}
  ConstSurfacePtrT(std::shared_ptr<surface_t> surface)
      : m_shared{std::move(surface)} {}
  ConstSurfacePtrT(const surface_t* surface) : m_shared{surface} {}
  ConstSurfacePtrT(std::shared_ptr<const surface_t> surface)
      : m_shared{std::move(surface)} {}
  ConstSurfacePtrT() = default;
  ConstSurfacePtrT(const ConstSurfacePtrT&) = default;
  ConstSurfacePtrT(ConstSurfacePtrT&&) = default;

  ConstSurfacePtrT(const SurfacePtrT<surface_t>& other)
      : m_shared{other.m_shared} {}
  ConstSurfacePtrT(SurfacePtrT<surface_t>&& other)
      : m_shared{std::move(other.m_shared)} {}

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT& operator=(const ConstSurfacePtrT<T>& other) {
    m_shared = other.m_shared;
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT& operator=(ConstSurfacePtrT<T>&& other) {
    m_shared = std::move(other.m_shared);
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT& operator=(const SurfacePtrT<T>& other) {
    m_shared = other.m_shared;
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT& operator=(SurfacePtrT<T>&& other) {
    m_shared = std::move(other.m_shared);
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT(const ConstSurfacePtrT<T>& other)
      : m_shared{other.m_shared} {}

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT(ConstSurfacePtrT<T>&& other)
      : m_shared{std::move(other.m_shared)} {}

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT(const SurfacePtrT<T>& other) : m_shared{other.m_shared} {}

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<surface_t, T>>>
  ConstSurfacePtrT(SurfacePtrT<T>&& other)
      : m_shared{std::move(other.m_shared)} {}

  ConstSurfacePtrT& operator=(const ConstSurfacePtrT& other) = default;
  ConstSurfacePtrT& operator=(ConstSurfacePtrT&& other) = default;

  const surface_t& operator*() const { return *m_shared; }
  const surface_t* operator->() const { return m_shared.get(); }
  operator bool() const { return !!m_shared; }
  const surface_t* get() const { return m_shared.get(); }

 private:
  std::shared_ptr<const surface_t> m_shared;
};

class Surface;

using SurfacePtr = SurfacePtrT<Surface>;
using ConstSurfacePtr = ConstSurfacePtrT<Surface>;

}  // namespace Acts
