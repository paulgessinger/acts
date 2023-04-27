// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/Podio/PodioTrackStateContainer.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Tests/CommonHelpers/MultiTrajectoryTestsCommon.hpp"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::Test;
namespace bd = boost::unit_test::data;

std::default_random_engine rng(31415);

class NullHelper : public PodioUtil::ConversionHelper {
 public:
  std::optional<PodioUtil::Identifier> surfaceToIdentifier(
      const Surface&) const override {
    return {};
  }
  const Surface* identifierToSurface(PodioUtil::Identifier) const override {
    return nullptr;
  }
};

struct MapHelper : public PodioUtil::ConversionHelper {
  std::optional<PodioUtil::Identifier> surfaceToIdentifier(
      const Surface& surface) const override {
    for (auto&& [id, srf] : surfaces) {
      if (srf == &surface) {
        return id;
      }
    }
    return {};
  }
  const Surface* identifierToSurface(PodioUtil::Identifier id) const override {
    auto it = surfaces.find(id);
    if (it == surfaces.end()) {
      return nullptr;
    }

    return it->second;
  }
  std::unordered_map<PodioUtil::Identifier, const Surface*> surfaces;
};

struct Factory {
  using trajectory_t = MutablePodioTrackStateContainer;
  using const_trajectory_t = ConstPodioTrackStateContainer;

  ActsPodioEdm::TrackStateCollection m_collection;
  ActsPodioEdm::BoundParametersCollection m_params;
  NullHelper m_helper;

  MutablePodioTrackStateContainer create() {
    return {m_helper, m_collection, m_params};
  }
};

using CommonTests = MultiTrajectoryTestsCommon<Factory>;

}  // namespace

BOOST_AUTO_TEST_SUITE(PodioTrackStateContainerTest)

BOOST_AUTO_TEST_CASE(Build) {
  CommonTests ct;
  ct.testBuild();
}

// BOOST_AUTO_TEST_CASE(ConstCorrectness) {
// CommonTests ct;
// ct.testConstCorrectness();
// }

BOOST_AUTO_TEST_CASE(Clear) {
  CommonTests ct;
  ct.testClear();
}

BOOST_AUTO_TEST_CASE(ApplyWithAbort) {
  CommonTests ct;
  ct.testApplyWithAbort();
}

BOOST_AUTO_TEST_CASE(AddTrackStateWithBitMask) {
  CommonTests ct;
  ct.testAddTrackStateWithBitMask();
}

// assert expected "cross-talk" between trackstate proxies
BOOST_AUTO_TEST_CASE(TrackStateProxyCrossTalk) {
  CommonTests ct;
  ct.testTrackStateProxyCrossTalk(rng);
}

#if 0

BOOST_AUTO_TEST_CASE(TrackStateReassignment) {
  CommonTests ct;
  ct.testTrackStateReassignment(rng);
}

BOOST_DATA_TEST_CASE(TrackStateProxyStorage, bd::make({1u, 2u}),
                     nMeasurements) {
  CommonTests ct;
  ct.testTrackStateProxyStorage(rng, nMeasurements);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyAllocations) {
  CommonTests ct;
  ct.testTrackStateProxyAllocations(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyGetMask) {
  CommonTests ct;
  ct.testTrackStateProxyGetMask();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopy) {
  CommonTests ct;
  ct.testTrackStateProxyCopy(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopyDiffMTJ) {
  CommonTests ct;
  ct.testTrackStateProxyCopyDiffMTJ();
}

BOOST_AUTO_TEST_CASE(ProxyAssignment) {
  CommonTests ct;
  ct.testProxyAssignment();
}

BOOST_AUTO_TEST_CASE(CopyFromConst) {
  CommonTests ct;
  ct.testCopyFromConst();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyShare) {
  CommonTests ct;
  ct.testTrackStateProxyShare(rng);
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumns) {
  CommonTests ct;
  ct.testMultiTrajectoryExtraColumns();
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumnsRuntime) {
  CommonTests ct;
  ct.testMultiTrajectoryExtraColumnsRuntime();
}

#endif

BOOST_AUTO_TEST_SUITE_END()
