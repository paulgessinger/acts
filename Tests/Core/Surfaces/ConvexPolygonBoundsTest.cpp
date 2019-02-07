// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE ConvexPolygonBounds Test

#include <boost/test/included/unit_test.hpp>

#include <chrono>
#include <iostream>
#include <memory>

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

using vec2 = Acts::Vector2D;
template <int N>
using poly = Acts::ConvexPolygonBounds<N>;

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_SUITE(Surfaces)

  BOOST_AUTO_TEST_CASE(construction_test)
  {

    std::vector<vec2> vertices;

    // triangle
    vertices = {{0, 0}, {1, 0}, {0.5, 1}};
    poly<3> triangle(vertices);
    std::cout << triangle << std::endl;
    std::cout << triangle.boundingBox() << std::endl;

    // rectangular poly
  }
  BOOST_AUTO_TEST_SUITE_END()
}
}
