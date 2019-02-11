// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE GenericCuboidVolumeBounds Test

#include <boost/test/included/unit_test.hpp>

#include <chrono>
#include <iostream>
#include <memory>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Visualization.hpp"
#include "Acts/Volumes/GenericCuboidVolumeBounds.hpp"

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_SUITE(Volumes)

  BOOST_AUTO_TEST_CASE(construction_test)
  {

    std::array<Vector3D, 8> vertices;
    vertices = {{{0, 0, 0},
                 {2, 0, 0},
                 {2, 1, 0},
                 {0, 1, 0},
                 {0, 0, 1},
                 {2, 0, 1},
                 {2, 1, 1},
                 {0, 1, 1}}};
    GenericCuboidVolumeBounds cubo(vertices);

    BOOST_CHECK(cubo.inside({0.5, 0.5, 0.5}));
    BOOST_CHECK(cubo.inside({1.5, 0.5, 0.5}));
    BOOST_CHECK(!cubo.inside({2.5, 0.5, 0.5}));
    BOOST_CHECK(!cubo.inside({0.5, 1.5, 0.5}));
    BOOST_CHECK(!cubo.inside({0.5, 0.5, 1.5}));
    BOOST_CHECK(!cubo.inside({-0.5, 0.5, 0.5}));
  }

  BOOST_AUTO_TEST_CASE(ply_test)
  {
    std::array<Vector3D, 8> vertices;
    vertices = {{{0, 0, 0},
                 {2, 0, 0},
                 {2, 1, 0},
                 {0, 1, 0},
                 {0, 0, 1},
                 {2, 0, 1},
                 {2, 1, 1},
                 {0, 1, 1}}};
    GenericCuboidVolumeBounds cubo(vertices);
    ply_helper<double>        ply;
    cubo.draw(ply);

    std::ofstream os("cuboid.ply");
    os << ply << std::flush;
    os.close();
  }

  BOOST_AUTO_TEST_CASE(bounding_box_creation)
  {
    float tol = 1e-6;
    std::array<Vector3D, 8> vertices;
    vertices = {{{-1, -1, -2},
                 {2, 0, 0},
                 {2, 1, 0},
                 {0, 1, 0},
                 {0, 0, 1},
                 {2, 0, 1},
                 {2, 1, 1},
                 {0, 1, 1}}};

    GenericCuboidVolumeBounds gcvb(vertices);
    auto                      bb = gcvb.boundingBox();

    Transform3D rot;
    rot = AngleAxis3D(M_PI / 2., Vector3D::UnitX());

    BOOST_CHECK_EQUAL(bb.entity(), nullptr);
    BOOST_CHECK_EQUAL(bb.max(), Vector3F(2, 1, 1));
    BOOST_CHECK_EQUAL(bb.min(), Vector3F(-1, -1, -2));

    bb = gcvb.boundingBox(&rot);

    BOOST_CHECK_EQUAL(bb.entity(), nullptr);
    BOOST_CHECK_EQUAL(bb.max(), Vector3F(2, 2, 1));
    BOOST_CHECK_EQUAL(bb.min(), Vector3F(-1, -1, -1));

    rot = AngleAxis3D(M_PI / 2., Vector3D::UnitZ());
    bb  = gcvb.boundingBox(&rot);
    BOOST_CHECK_EQUAL(bb.entity(), nullptr);
    BOOST_CHECK_EQUAL(bb.max(), Vector3F(1, 2, 1));
    BOOST_CHECK_EQUAL(bb.min(), Vector3F(-1, -1, -2));

    rot = AngleAxis3D(0.542, Vector3D::UnitZ())
        * AngleAxis3D(1.424, Vector3D::UnitX())
        * AngleAxis3D(0.9941, Vector3D::UnitY());

    bb = gcvb.boundingBox(&rot);

    ply_helper<float>  plyf;
    ply_helper<double> plyd;
    gcvb.draw(plyd, &rot);
    bb.draw(plyf);

    std::ofstream os("test_gcvb.ply");
    os << plyd;
    os.close();
    os = std::ofstream("test_bb.ply");
    os << plyf;
    os.close();

    BOOST_CHECK_EQUAL(bb.entity(), nullptr);
    BOOST_CHECK_EQUAL(bb.max(), Vector3F(1, 2, 1));
    BOOST_CHECK_EQUAL(bb.min(), Vector3F(-1, -1, -2));
  }

  BOOST_AUTO_TEST_SUITE_END()
}
}
