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
#include "Acts/Tests/CommonHelpers/MultiTrajectoryTestsCommon.hpp"

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::Test;
namespace bd = boost::unit_test::data;

std::default_random_engine rng(31415);

using CommonTests = MultiTrajectoryTestsCommon<MutablePodioTrackStateContainer,
                                               ConstPodioTrackStateContainer>;

}  // namespace

BOOST_AUTO_TEST_SUITE(PodioTrackStateContainerTest)

BOOST_AUTO_TEST_CASE(Build) {
  CommonTests::testBuild();
}

// BOOST_AUTO_TEST_CASE(ConstCorrectness) {
// CommonTests::testConstCorrectness();
// }

BOOST_AUTO_TEST_CASE(Clear) {
  CommonTests::testClear();
}

BOOST_AUTO_TEST_CASE(ApplyWithAbort) {
  CommonTests::testApplyWithAbort();
}

BOOST_AUTO_TEST_CASE(AddTrackStateWithBitMask) {
  CommonTests::testAddTrackStateWithBitMask();
}

// assert expected "cross-talk" between trackstate proxies
BOOST_AUTO_TEST_CASE(TrackStateProxyCrossTalk) {
  CommonTests::testTrackStateProxyCrossTalk(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateReassignment) {
  CommonTests::testTrackStateReassignment(rng);
}

// BOOST_DATA_TEST_CASE(TrackStateProxyStorage, bd::make({1u, 2u}),
// nMeasurements) {
// CommonTests::testTrackStateProxyStorage(rng, nMeasurements);
// }

BOOST_AUTO_TEST_CASE(TrackStateProxyAllocations) {
  CommonTests::testTrackStateProxyAllocations(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyGetMask) {
  CommonTests::testTrackStateProxyGetMask();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopy) {
  CommonTests::testTrackStateProxyCopy(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopyDiffMTJ) {
  CommonTests::testTrackStateProxyCopyDiffMTJ();
}

BOOST_AUTO_TEST_CASE(ProxyAssignment) {
  CommonTests::testProxyAssignment();
}

BOOST_AUTO_TEST_CASE(CopyFromConst) {
  CommonTests::testCopyFromConst();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyShare) {
  CommonTests::testTrackStateProxyShare(rng);
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumns) {
  CommonTests::testMultiTrajectoryExtraColumns();
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumnsRuntime) {
  CommonTests::testMultiTrajectoryExtraColumnsRuntime();
}

BOOST_AUTO_TEST_SUITE_END()
