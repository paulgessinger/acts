// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <functional>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Acts {

template <typename identifier_t, typename value_t>
class IdentifyableContainer {
 public:
  using value_store_t = std::vector<value_t>;

  using Iterator = typename value_store_t::iterator;

  auto rangesForIdentifier(const identifier_t& identifier) const {
    return m_ranges.at(identifier);
  }

  IdentifyableContainer(value_store_t values,
                        std::function<identifier_t(value_t)> mapper)
      : m_valueStore{std::move(values)} {
    // std::sort(m_valueStore.begin(), m_valueStore.end());

    std::optional<identifier_t> prev = nullptr;
    Iterator groupStart = m_valueStore.begin();
    // Iterator groupEnd;
    for (Iterator it = m_valueStore.begin(); it != m_valueStore.end(); ++it) {
      identifier_t id = mapper(*it);
      // const auto& value = *it;
      if (prev && (*prev) != id) {
        // new range
        m_ranges[id] = {groupStart, it};
        groupStart = it;
      }
    }
    // close the last range
    m_ranges[mapper(m_valueStore.back())] = {groupStart, m_valueStore.end()};
  }

 private:
  std::unordered_map<identifier_t, std::pair<Iterator, Iterator>> m_ranges;
  std::vector<value_t> m_valueStore;
};

}  // namespace Acts
