// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <memory>

namespace Acts {

template <typename element_t>
class ConstAsymHandle;

namespace detail {

struct AsymHandleControlBlock {
  bool sharedOwnership = true;
  std::atomic<size_t> refCount;
};

}  // namespace detail

template <typename element_t>
class AsymHandle {
 public:
  using element_type = element_t;

  // AsymHandle(std::shared_ptr<element_t> surface)
  // : m_shared{std::move(surface)} {}

  // AsymHandle(std::nullptr_t) {}

  AsymHandle(std::unique_ptr<element_t> up) noexcept {
    m_cb = new detail::AsymHandleControlBlock{true, 1};
    m_ptr = up.release();
  }

 public:
  AsymHandle() = default;
  AsymHandle(const AsymHandle& other) {
    m_ptr = other.m_ptr;
    m_cb = other.m_cb;
    increment();
  }

  AsymHandle(AsymHandle&& other) {
    m_ptr = other.m_ptr;
    m_cb = other.m_cb;
    // don't need to increment, ref count stays the same
  }

  AsymHandle& operator=(const AsymHandle& other) = default;
  AsymHandle& operator=(AsymHandle&& other) = default;

  friend bool operator==(const AsymHandle<element_t>& lhs, std::nullptr_t) {
    return !lhs;
  }

  friend bool operator!=(const AsymHandle<element_t>& lhs, std::nullptr_t) {
    return lhs;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
  AsymHandle& operator=(const AsymHandle<T>& other) {
    m_shared = other.m_shared;
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
  AsymHandle& operator=(AsymHandle<T>&& other) {
    m_shared = std::move(other.m_shared);
    return *this;
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
  AsymHandle(const AsymHandle<T>& other) : m_shared{other.m_shared} {}

  template <typename T,
            typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
  AsymHandle(AsymHandle<T>&& other) : m_shared{std::move(other.m_shared)} {}

 private:
  void increment() {}
  void decrement() {}

 public:
  element_t& operator*() { return *m_shared; }
  element_t* operator->() { return m_shared.get(); }
  operator bool() const { return !!m_shared; }
  element_t* get() { return m_shared.get(); }

  const element_t& operator*() const { return *m_shared; }
  const element_t* operator->() const { return m_shared.get(); }
  const element_t* get() const { return m_shared.get(); }

  size_t refCount() const { return m_shared.use_count(); }

  template <typename T>
  friend class ConstAsymHandle;
  template <typename T>
  friend class AsymHandle;

 private:
  detail::AsymHandleControlBlock* m_cb{nullptr};
  element_t* m_ptr{nullptr};
  // std::shared_ptr<element_t> m_shared;
};

template <typename T, typename... Args>
AsymHandle<T> makeAsymHandle(Args&&... args) {
  auto up = std::make_unique<T>(std::forward<Args>(args)...);
  return AsymHandle<T>{up.release()};
}

// template <typename element_t>
// class ConstAsymHandle {
// public:
// using element_type = element_t;

// ConstAsymHandle(element_t* surface) : m_shared{surface} {}
// ConstAsymHandle(std::shared_ptr<element_t> surface)
// : m_shared{std::move(surface)} {}

// ConstAsymHandle(std::nullptr_t) {}

// ConstAsymHandle(const element_t* surface) : m_shared{surface} {}
// ConstAsymHandle(std::shared_ptr<const element_t> surface)
// : m_shared{std::move(surface)} {}
// ConstAsymHandle() = default;
// ConstAsymHandle(const ConstAsymHandle&) = default;
// ConstAsymHandle(ConstAsymHandle&&) = default;

// ConstAsymHandle(const AsymHandle<element_t>& other)
// : m_shared{other.m_shared} {}
// ConstAsymHandle(AsymHandle<element_t>&& other)
// : m_shared{std::move(other.m_shared)} {}

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle& operator=(const ConstAsymHandle<T>& other) {
// m_shared = other.m_shared;
// return *this;
// }

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle& operator=(ConstAsymHandle<T>&& other) {
// m_shared = std::move(other.m_shared);
// return *this;
// }

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle& operator=(const AsymHandle<T>& other) {
// m_shared = other.m_shared;
// return *this;
// }

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle& operator=(AsymHandle<T>&& other) {
// m_shared = std::move(other.m_shared);
// return *this;
// }

// friend bool operator==(const ConstAsymHandle<element_t>& lhs,
// std::nullptr_t) {
// return !lhs;
// }

// friend bool operator!=(const ConstAsymHandle<element_t>& lhs,
// std::nullptr_t) {
// return lhs;
// }

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle(const ConstAsymHandle<T>& other) : m_shared{other.m_shared}
// {}

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle(ConstAsymHandle<T>&& other)
// : m_shared{std::move(other.m_shared)} {}

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle(const AsymHandle<T>& other) : m_shared{other.m_shared} {}

// template <typename T,
// typename = std::enable_if_t<std::is_base_of_v<element_t, T>>>
// ConstAsymHandle(AsymHandle<T>&& other)
// : m_shared{std::move(other.m_shared)} {}

// ConstAsymHandle& operator=(const ConstAsymHandle& other) = default;
// ConstAsymHandle& operator=(ConstAsymHandle&& other) = default;

// const element_t& operator*() const { return *m_shared; }
// const element_t* operator->() const { return m_shared.get(); }
// operator bool() const { return !!m_shared; }
// const element_t* get() const { return m_shared.get(); }

// size_t refCount() const { return m_shared.use_count(); }

// template <typename T>
// friend class ConstAsymHandle;

// private:
// std::shared_ptr<const element_t> m_shared;
// };

// template <typename T>
// bool operator==(const ConstAsymHandle<T>& lhs, std::nullptr_t) {
// return !lhs;
// }

}  // namespace Acts
