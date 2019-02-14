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
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

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

    BOOST_CHECK(!cubo.inside({2.2, 1, 1}, 0.1));
    BOOST_CHECK(cubo.inside({2.2, 1, 1}, 0.21));
    BOOST_CHECK(cubo.inside({2.2, 1, 1}, 0.3));

  }

  BOOST_AUTO_TEST_CASE(decomposeToSurfaces)
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

    ply_helper<double> ply;
    cubo.draw(ply);

    std::ofstream os("deco_direct.ply");
    os << ply;
    os.close();
    ply.clear();

    auto surfaces = cubo.decomposeToSurfaces(nullptr);
    for(const auto& srf : surfaces) {
      auto pbounds = dynamic_cast<const PlanarBounds*>(&srf->bounds());
      std::vector<Vector3D> vertices;
      for (const auto& vtx : pbounds->vertices()) {
        Vector3D glob;
        srf->localToGlobal(vtx, {}, glob);
        vertices.push_back(glob);
      }
      ply.face(vertices);
    }

    os = std::ofstream("deco_decomp.ply");
    os << ply;
    os.close();

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
    float tol = 1e-4;
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
        * AngleAxis3D(M_PI / 5., Vector3D(1, 3, 6).normalized());

    bb = gcvb.boundingBox(&rot);
    BOOST_CHECK_EQUAL(bb.entity(), nullptr);
    CHECK_CLOSE_ABS(bb.max(), Vector3F(1.09416, 2.35122, 1.11988), tol);
    CHECK_CLOSE_ABS(bb.min(), Vector3F(-0.871397, -1.61236, -1.84328), tol);
  }

  BOOST_AUTO_TEST_SUITE_END()
}
}
