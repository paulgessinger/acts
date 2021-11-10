// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/Utilities/IdentifyableContainer.hpp"

#include <numeric>
#include <optional>
#include <random>
#include <tuple>

using namespace Acts;

namespace bd = boost::unit_test::data;

BOOST_AUTO_TEST_SUITE(DelegateTests)

class PseudoSourceLink {
 public:
  PseudoSourceLink(int id, int value_) : m_id{id}, m_value{value_} {}

  int identifier() const { return m_id; }

  int value() const { return m_value; }

 private:
  int m_id;
  int m_value;
};

BOOST_AUTO_TEST_CASE(AccessPattern) {
  using container_t = IdentifyableContainer<int, PseudoSourceLink>;

  std::vector<PseudoSourceLink> values = {{1, 1}, {1, 2},  {1, 3},  {1, 4},
                                          {2, 5}, {2, 6},  {2, 7},  {1, 8},
                                          {1, 9}, {3, 10}, {3, 11}, {2, 12}};

  container_t container{
      values, [](const PseudoSourceLink& sl) { return sl.identifier(); }};

  for (int i = 1; i < 4; i++) {
    std::cout << "lookup: " << i << std::endl;
    auto ranges = container.rangesForIdentifier(i);
    for (auto& [it, end] : ranges) {
      for (; it != end; ++it) {
        const auto& sl = *it;
        std::cout << sl.identifier() << " " << sl.value() << std::endl;
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
