// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>

namespace Acts {

template <typename>
class Delegate;

/// Delegate type that allows type erasure of a callable without allocation and
/// with a single level of indirection.
/// This type can support:
/// - a free function pointer
/// - a pointer to a member function alongside an instance pointer
/// @note @c Delegate does not assume ownership of the instance.
///          You need to ensure that the lifetime of the callable
///          instance is longer than that of the @c Delegate.
/// @note Currently @c Delegate only supports callables that are ``const``
/// @tparam R Return type of the function signature
/// @tparam Args Types of the arguments of the function signatures
///
template <typename R, typename... Args>
class Delegate<R(Args...)> {
  /// Alias of the return type
  using return_type = R;
  /// Alias to the function pointer type this class will store
  using function_type = return_type (*)(const void*, Args...);

 public:
  Delegate() = default;

  /// Connect a free function pointer.
  /// @note The function pointer must be ``constexpr`` for @c Delegate to accept it
  /// @tparam Callable The compile-time free function pointer
  template <auto Callable>
  void connect() {
    m_function = [](const void* /*payload*/, Args... args) -> return_type {
      return std::invoke(Callable, std::forward<Args>(args)...);
    };
  }

  /// Connect a member function to be called on an instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate does not assume owner ship over @p instance. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <auto Callable, typename Type>
  void connect(const Type* instance) {
    m_payload = instance;
    m_function = [](const void* payload, Args... args) -> return_type {
      const auto* concretePayload = static_cast<const Type*>(payload);
      return std::invoke(Callable, concretePayload,
                         std::forward<Args>(args)...);
    };
  }

  /// The call operator that exposes the functionality of the @c Delegate type.
  /// @param args The arguments to call the contained function with
  /// @return Return value of the contained function
  return_type operator()(Args... args) const {
    assert(m_function != nullptr && "Delegate is not connected");
    return std::invoke(m_function, m_payload, std::forward<Args>(args)...);
  }

 private:
  /// Stores the instance pointer
  const void* m_payload{nullptr};
  /// Stores the function pointer wrapping the compile time function pointer given in @c connect().
  function_type m_function{nullptr};
};
}  // namespace Acts
