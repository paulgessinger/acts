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
#include "Acts/Volumes/GenericCuboidVolumeBounds.hpp"

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_SUITE(Volumes)

    BOOST_AUTO_TEST_CASE(construction_test) {

      std::array<Vector3D, 8> vertices;
      vertices = {{
        {0, 0, 0},
        {2, 0, 0},
        {2, 1, 0},
        {0, 1, 0},
        {0, 0, 1},
        {2, 0, 1},
        {2, 1, 1},
        {0, 1, 1}
      }};
      GenericCuboidVolumeBounds cubo(vertices);

      BOOST_CHECK(cubo.inside({0.5, 0.5, 0.5}));

    }
  
  BOOST_AUTO_TEST_SUITE_END()
}
}
