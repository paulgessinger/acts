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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Podio/PodioTrackContainer.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsPodioEdm/Surface.h"
#include <ActsPodioEdm/TrackCollection.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <random>

using namespace Acts;
using namespace Acts::UnitLiterals;
BOOST_AUTO_TEST_SUITE(PodioTrackConversion)

BOOST_AUTO_TEST_CASE(ConvertSurface) {
  auto rBounds = std::make_shared<RectangleBounds>(15, 20);

  auto trf = Transform3::Identity();
  trf.translation().setRandom();

  auto free = Acts::Surface::makeShared<PlaneSurface>(trf, rBounds);

  Acts::GeometryContext gctx;
  auto surface = convertSurfaceToPodio(gctx, *free);

  auto free2 = convertSurfaceFromPodio(surface);

  BOOST_REQUIRE(free2);
  BOOST_CHECK_EQUAL(free->type(), free2->type());
  BOOST_CHECK_EQUAL(free->bounds().type(), free2->bounds().type());
  BOOST_CHECK_EQUAL(free->center(gctx), free2->center(gctx));

  const auto* rBounds2 = dynamic_cast<const RectangleBounds*>(&free2->bounds());
  BOOST_REQUIRE_NE(rBounds2, nullptr);

  BOOST_CHECK_EQUAL(rBounds2->halfLengthX(), rBounds->halfLengthX());
  BOOST_CHECK_EQUAL(rBounds2->halfLengthY(), rBounds->halfLengthY());

  // this could probably use some more complete checks
}

BOOST_AUTO_TEST_CASE(ConvertTrack) {
  ActsPodioEdm::TrackCollection tracks;

  Acts::VectorMultiTrajectory mtj{};
  Acts::PodioTrackContainer ptc{tracks};

  Acts::TrackContainer tc{ptc, mtj};

  BOOST_CHECK_EQUAL(tc.size(), 0);

  auto t = tc.getTrack(tc.addTrack());
  BOOST_CHECK_EQUAL(t.tipIndex(), MultiTrajectoryTraits::kInvalid);
  t.tipIndex() = 5;
  BOOST_CHECK_EQUAL(t.tipIndex(), 5);

  BOOST_CHECK_EQUAL(tc.size(), 1);

  auto pTrack = tracks.at(0);
  BOOST_CHECK_EQUAL(pTrack.data().tipIndex, 5);

  t.parameters() << 1, 2, 3, 4, 5, 6;
  Eigen::Map<BoundVector> pars{pTrack.data().parameters.data()};
  BOOST_CHECK_EQUAL(pars, (BoundVector{1, 2, 3, 4, 5, 6}));

  auto ref = BoundMatrix::Random().eval();
  t.covariance() = ref;

  Eigen::Map<const BoundMatrix> cov{pTrack.data().covariance.data()};
  BOOST_CHECK_EQUAL(ref, cov);

  t.nMeasurements() = 17;
  BOOST_CHECK_EQUAL(pTrack.data().nMeasurements, 17);

  t.nHoles() = 34;
  BOOST_CHECK_EQUAL(pTrack.data().nHoles, 34);

  t.chi2() = 882.3f;
  BOOST_CHECK_EQUAL(pTrack.data().chi2, 882.3f);

  t.nDoF() = 9;
  BOOST_CHECK_EQUAL(pTrack.data().ndf, 9);

  t.nOutliers() = 77;
  BOOST_CHECK_EQUAL(pTrack.data().nOutliers, 77);

  t.nSharedHits() = 99;
  BOOST_CHECK_EQUAL(pTrack.data().nSharedHits, 99);

  auto rBounds = std::make_shared<RectangleBounds>(15, 20);
  auto trf = Transform3::Identity();
  trf.translation().setRandom();
  auto free = Acts::Surface::makeShared<PlaneSurface>(trf, rBounds);

  Acts::GeometryContext gctx;
  t.setReferenceSurface(free);
  const auto& free2 = t.referenceSurface();
  BOOST_CHECK_EQUAL(free->center(gctx), free2.center(gctx));

  const auto* rBounds2 = dynamic_cast<const RectangleBounds*>(&free2.bounds());
  BOOST_REQUIRE_NE(rBounds2, nullptr);

  BOOST_CHECK_EQUAL(rBounds2->halfLengthX(), rBounds->halfLengthX());
  BOOST_CHECK_EQUAL(rBounds2->halfLengthY(), rBounds->halfLengthY());

  // std::cout << track.getL0() << std::endl;
  // std::cout << track.getT() << std::endl;

  // auto refSurface = Surface::makeShared<PerigeeSurface>(Vector3{50, 30, 20});

  // BoundVector par;
  // par << 1_mm, 5_mm, 0, M_PI_2, -1 / 1_GeV,
  // 5_ns;  // -> perpendicular to perigee and pointing right, should be PCA

  // BoundMatrix cov;
  // cov.setIdentity();
  // cov(5, 5) = 25_ns;

  // SingleBoundTrackParameters<SinglyCharged> boundPar{refSurface, par, cov};

  // double Bz = 2_T;

  // Acts::GeometryContext gctx;

  // EDM4hepUtil::detail::Parameters converted =
  // EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  // BOOST_CHECK(converted.covariance.has_value());
  // BOOST_CHECK(converted.surface);

  // // input is already on perigee, should not be modified
  // BOOST_CHECK_EQUAL(par.template head<2>(),
  // converted.values.template head<2>());
  // BOOST_CHECK_EQUAL(
  // (converted.covariance.value().template topLeftCorner<4, 4>()),
  // ActsSymMatrix<4>::Identity());
  // BOOST_CHECK(converted.covariance.value()(4, 4) > 0);
  // BOOST_CHECK_EQUAL(converted.covariance.value()(5, 5), 25_ns);

  // // convert back for roundtrip test

  // SingleBoundTrackParameters<SinglyCharged> roundtripPar =
  // EDM4hepUtil::detail::convertTrackParametersFromEdm4hep(Bz, converted);

  // BOOST_CHECK(roundtripPar.parameters().isApprox(boundPar.parameters()));
  // BOOST_CHECK(roundtripPar.covariance().value().isApprox(
  // boundPar.covariance().value()));
}

BOOST_AUTO_TEST_SUITE_END()
