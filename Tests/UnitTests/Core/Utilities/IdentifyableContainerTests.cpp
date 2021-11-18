// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
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

// BOOST_AUTO_TEST_CASE(AccessPattern) {
// using container_t = IdentifyableContainer<int, PseudoSourceLink>;

// std::vector<PseudoSourceLink> values = {{1, 1}, {1, 2},  {1, 3},  {1, 4},
// {2, 5}, {2, 6},  {2, 7},  {1, 8},
// {1, 9}, {3, 10}, {3, 11}, {2, 12}};

// container_t container{
// values, [](const PseudoSourceLink& sl) { return sl.identifier(); }};

// for (int i = 1; i < 4; i++) {
// std::cout << "lookup: " << i << std::endl;
// auto ranges = container.rangesForIdentifier(i);
// for (auto& [it, end] : ranges) {
// for (; it != end; ++it) {
// const auto& sl = *it;
// std::cout << sl.identifier() << " " << sl.value() << std::endl;
// }
// }
// }
// }

BOOST_AUTO_TEST_CASE(Identifier) {
  using module_identifier_t = int;

  struct RDO {
    module_identifier_t module_id;
    int rdo_id;
  };

  std::vector<RDO> rdos;
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 6; j++) {
      rdos.push_back(RDO{i, j});  // RDO number j on module number i
    }
  }

  using rdo_container_t = IdentifyableContainer<module_identifier_t, RDO>;

  rdo_container_t rdoContainer{rdos,
                               [](const RDO& sl) { return sl.module_id; }};

  struct Cluster {
    rdo_container_t::Range rdoRange;
    Vector2 localPosition;
    SymMatrix2 localCovariance;
    module_identifier_t module_id;
  };

  std::vector<Cluster> clusters;

  // loop over modules
  for (int i = 0; i < 10; i++) {
    std::cout << "lookup rdos for module: " << i << std::endl;
    auto ranges = rdoContainer.rangesForIdentifier(i);

    for (auto& [start, end] : ranges) {
      auto it = start;
      for (; it != end; ++it) {
        const auto& rdo = *it;
        std::cout << " - " << rdo.module_id << " " << rdo.rdo_id << std::endl;
      }
      clusters.push_back(Cluster{
          {start, std::next(start, 3)}, Vector2{1, 2}, SymMatrix2{}, i});
      clusters.push_back(
          Cluster{{std::next(start, 3), end}, Vector2{1, 2}, SymMatrix2{}, i});
    }
    // for (auto& [it, end] : ranges) {
    // clusters.push_back(Cluster{{it, end}, Vector2{1, 2}, SymMatrix2{}, i});
    // clusters.push_back(
    // Cluster{{it, std::next(it, 3)}, Vector2{1, 2}, SymMatrix2{}, i});
    // clusters.push_back(
    // Cluster{{std::next(it, 3), end}, Vector2{1, 2}, SymMatrix2{}, i});
    // }
  }

  using cluster_container_t =
      IdentifyableContainer<module_identifier_t, Cluster>;

  cluster_container_t clusterConainer{
      clusters, [](const Cluster& c) { return c.module_id; }};

  // loop over modules
  for (int i = 0; i < 10; i++) {
    std::cout << "lookup clusters for module: " << i << std::endl;
    auto ranges = clusterConainer.rangesForIdentifier(i);
    int c = 1;
    for (auto& [it, end] : ranges) {
      for (; it != end; ++it) {
        auto& cluster = *it;
        std::cout << "- cluster " << c << std::endl;
        auto [rdo_it, rdo_end] = cluster.rdoRange;
        for (; rdo_it != rdo_end; ++rdo_it) {
          std::cout << "  - rdo: " << rdo_it->rdo_id << std::endl;
        }
        c++;
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
