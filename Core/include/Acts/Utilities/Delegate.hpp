// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename>
struct Delegate;

template <typename Ret, typename... Args>
struct Delegate<Ret(Args...)> {
  using function_type = Ret (*)(const void*, Args...);
  using type = Ret(Args...);
  using return_type = Ret;

 public:
  Delegate() = default;

  template <auto Callable>
  void connect() {
    m_function = [](const void*, Args... args) -> Ret {
      return std::invoke(Callable, std::forward<Args>(args)...);
    };
  }

  template <auto Callable, typename Type>
  void connect(Type* _type) {
    m_payload = _type;
    m_function = [](const void* payload, Args... args) -> Ret {
      const auto* __type = static_cast<const Type*>(payload);
      return std::invoke(Callable, __type, std::forward<Args>(args)...);
    };
  }

  Ret operator()(Args... args) const {
    assert(m_function != nullptr && "Delegate is not connected");
    return std::invoke(m_function, m_payload, std::forward<Args>(args)...);
  }

 private:
  void* m_payload{nullptr};
  function_type m_function{nullptr};
};
}  // namespace Acts