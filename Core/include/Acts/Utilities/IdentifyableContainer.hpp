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
#include <iostream>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

template <typename T>
struct RangeTag {};

template <typename value_t>
struct RangeProxy {
  struct iterator {
    uint8_t* ptr;
    size_t size;

    value_t& operator*() { return *reinterpret_cast<value_t*>(ptr); }

    void operator++() { ptr += size; }

    bool operator==(const iterator& other) { return ptr == other.ptr; }
    bool operator!=(const iterator& other) { return !(*this == other); }
  };

  iterator start;
  iterator stop;

  iterator begin() { return start; }
  iterator end() { return stop; }
};

template <typename identifier_t, typename value_t, size_t inline_size = 10>
class IdentifyableContainerAccess {
 public:
  virtual boost::container::small_vector<RangeProxy<value_t>, inline_size>
  rangesForIdentifier(const identifier_t& identifier,
                      RangeTag<value_t>) const = 0;

  virtual ~IdentifyableContainerAccess() = 0;
};

template <typename identifier_t, typename value_t, size_t inline_size>
inline IdentifyableContainerAccess<
    identifier_t, value_t, inline_size>::~IdentifyableContainerAccess() =
    default;

template <typename identifier_t, typename value_t, size_t inline_size = 10>
class IdentifyableContainer
    : public IdentifyableContainerAccess<identifier_t, value_t, inline_size> {
 public:
  using value_store_t = std::vector<value_t>;

  using Iterator = typename value_store_t::iterator;
  using Range = std::pair<Iterator, Iterator>;

  auto rangesForIdentifierDirect(const identifier_t& identifier) const {
    return m_ranges.at(identifier);
  }

  boost::container::small_vector<RangeProxy<value_t>, inline_size>
  rangesForIdentifier(const identifier_t& identifier,
                      RangeTag<value_t>) const override {
    boost::container::small_vector<RangeProxy<value_t>, inline_size> ranges;

    for (auto [r_start, r_end] : rangesForIdentifierDirect(identifier)) {
      typename RangeProxy<value_t>::iterator start{
          reinterpret_cast<uint8_t*>(&*r_start), sizeof(value_t)};
      typename RangeProxy<value_t>::iterator end{
          reinterpret_cast<uint8_t*>(&*r_end), sizeof(value_t)};

      ranges.push_back(RangeProxy<value_t>{start, end});
    }

    return ranges;
  }

  IdentifyableContainer(value_store_t values,
                        std::function<identifier_t(value_t)> mapper)
      : m_valueStore{std::move(values)} {
    // std::sort(m_valueStore.begin(), m_valueStore.end());

    std::optional<identifier_t> prev = std::nullopt;
    Iterator groupStart = m_valueStore.begin();
    // Iterator groupEnd;
    for (Iterator it = m_valueStore.begin(); it != m_valueStore.end(); ++it) {
      identifier_t id = mapper(*it);
      // const auto& value = *it;
      if (prev && (*prev) != id) {
        // new range
        std::cout << "add range: " << *prev << ": "
                  << std::distance(m_valueStore.begin(), groupStart) << " -> "
                  << std::distance(m_valueStore.begin(), it) << std::endl;
        m_ranges[*prev].emplace_back(groupStart, it);
        groupStart = it;
      }
      prev = id;
    }
    // close the last range
    std::cout << "add range: " << *prev << ": "
              << std::distance(m_valueStore.begin(), groupStart) << " -> "
              << std::distance(m_valueStore.begin(), m_valueStore.end())
              << std::endl;
    m_ranges[*prev].emplace_back(groupStart, m_valueStore.end());
  }

 private:
  std::unordered_map<identifier_t,
                     boost::container::small_vector<Range, inline_size>>
      m_ranges;
  std::vector<value_t> m_valueStore;
};

}  // namespace Acts
