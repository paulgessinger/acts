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

#include <boost/container/small_vector.hpp>

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

using module_identifier_t = int;
struct RDO {
  module_identifier_t module_id;
  int rdo_id;
};
using rdo_container_t = IdentifyableContainer<module_identifier_t, RDO>;
struct PRD {
  boost::container::small_vector<RDO, 10> rdos;
  Vector2 localPosition;
  SymMatrix2 localCovariance;
  module_identifier_t module_id;
};

using prd_container_t = IdentifyableContainer<module_identifier_t, PRD>;

struct PixelCluster : public PRD {
  double totalCharge;
  boost::container::small_vector<int, 10> tot;

  PixelCluster(boost::container::small_vector<RDO, 10> rdos_,
               Vector2 localPosition_, SymMatrix2 localCovariance_,
               module_identifier_t module_id_, double totalCharge_,
               boost::container::small_vector<int, 10> tot_)
      : PRD{rdos_, localPosition_, localCovariance_, module_id_},
        totalCharge{totalCharge_},
        tot{tot_} {}
};

struct PixelClusterContainer
    : public IdentifyableContainer<module_identifier_t, PixelCluster>,
      public IdentifyableContainerAccess<module_identifier_t, PRD> {
  using IdentifyableContainer<module_identifier_t,
                              PixelCluster>::IdentifyableContainer;

  boost::container::small_vector<RangeProxy<PRD>, 10> rangesForIdentifier(
      const module_identifier_t& identifier, RangeTag<PRD>) const override {
    boost::container::small_vector<RangeProxy<PRD>, 10> ranges;

    for (auto [r_start, r_end] : rangesForIdentifierDirect(identifier)) {
      typename RangeProxy<PRD>::iterator start{
          reinterpret_cast<uint8_t*>(&*r_start), sizeof(PixelCluster)};
      typename RangeProxy<PRD>::iterator end{
          reinterpret_cast<uint8_t*>(&*r_end), sizeof(PixelCluster)};

      ranges.push_back(RangeProxy<PRD>{start, end});
    }

    return ranges;
  }
};

void processPRDs(
    const IdentifyableContainerAccess<module_identifier_t, PRD>& prds) {
  for (auto [start, stop] : prds.rangesForIdentifier(2, RangeTag<PRD>{})) {
    for (auto it = start; it != stop; ++it) {
      std::cout << "prd: " << (*it).module_id << " rdos:";
      for (const auto& rdo : (*it).rdos) {
        std::cout << " " << rdo.rdo_id;
      }
      std::cout << std::endl;
    }
  }
}

void processPixelClusters(
    const IdentifyableContainerAccess<module_identifier_t, PixelCluster>&
        clusters) {
  for (auto [start, stop] :
       clusters.rangesForIdentifier(2, RangeTag<PixelCluster>{})) {
    for (auto it = start; it != stop; ++it) {
      std::cout << "prd: " << (*it).module_id << " rdos:";
      for (const auto& rdo : (*it).rdos) {
        std::cout << " " << rdo.rdo_id;
      }
      std::cout << std::endl;
      std::cout << "totalCharge: " << (*it).totalCharge << std::endl;
    }
  }
}

BOOST_AUTO_TEST_CASE(Identifier) {
  std::vector<RDO> rdos;
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 6; j++) {
      rdos.push_back(RDO{i, j});  // RDO number j on module number i
    }
  }

  rdo_container_t rdoContainer{rdos,
                               [](const RDO& sl) { return sl.module_id; }};

  std::vector<PRD> prds;

  // loop over modules
  for (int i = 0; i < 10; i++) {
    // std::cout << "lookup rdos for module: " << i << std::endl;
    auto ranges = rdoContainer.rangesForIdentifierDirect(i);

    for (auto& [start, end] : ranges) {
      auto it = start;
      for (; it != end; ++it) {
        // const auto& rdo = *it;
        // std::cout << " - " << rdo.module_id << " " << rdo.rdo_id <<
        // std::endl;
      }
      prds.push_back(
          PRD{{start, std::next(start, 3)}, Vector2{1, 2}, SymMatrix2{}, i});
      prds.push_back(
          PRD{{std::next(start, 3), end}, Vector2{1, 2}, SymMatrix2{}, i});
    }
  }

  prd_container_t prdContainer{prds, [](const PRD& c) { return c.module_id; }};

  // loop over modules
  // for (int i = 0; i < 10; i++) {
  // std::cout << "lookup clusters for module: " << i << std::endl;
  // auto ranges = prdContainer.rangesForIdentifier(i);
  // int c = 1;
  // for (auto& [it, end] : ranges) {
  // for (; it != end; ++it) {
  // auto& cluster = *it;
  // std::cout << "- cluster " << c << std::endl;
  // for (const auto rdo : cluster.rdos) {
  // std::cout << "  - rdo: " << rdo.rdo_id << std::endl;
  // }
  // c++;
  // }
  // }
  // }

  //
  std::vector<PixelCluster> clusters;

  // loop over modules
  for (int i = 0; i < 10; i++) {
    // std::cout << "lookup rdos for module: " << i << std::endl;
    auto ranges = rdoContainer.rangesForIdentifierDirect(i);

    for (auto& [start, end] : ranges) {
      auto it = start;
      // for (; it != end; ++it) {
      // const auto& rdo = *it;
      // std::cout << " - " << rdo.module_id << " " << rdo.rdo_id <<
      // std::endl;
      // }
      clusters.push_back(PixelCluster{{start, std::next(start, 3)},
                                      Vector2{1, 2},
                                      SymMatrix2{},
                                      i,
                                      5,
                                      {3, 4}});
      clusters.push_back(PixelCluster{{std::next(start, 3), end},
                                      Vector2{1, 2},
                                      SymMatrix2{},
                                      i,
                                      9,
                                      {3, 1}});
    }
  }

  PixelClusterContainer pixClusterContainer{
      clusters, [](const PixelCluster& c) { return c.module_id; }};

  processPRDs(prdContainer);
  processPRDs(pixClusterContainer);

  processPixelClusters(pixClusterContainer);
}

BOOST_AUTO_TEST_SUITE_END()
