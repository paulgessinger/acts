// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AsymHandle.hpp"

#include <utility>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(AsymHandleTests)

struct A {
  std::function<void()> onDestroy;

  A(std::function<void()> fn) : onDestroy{fn} {}

  ~A() { onDestroy(); }
};

BOOST_AUTO_TEST_CASE(RefCount) {
  bool destructorCalled = false;
  {
    std::cout << __LINE__ << std::endl;
    auto h1 = makeAsymHandle<A>([&]() { destructorCalled = true; });
    BOOST_CHECK_EQUAL(h1.refCount(), 1);
    BOOST_CHECK(!destructorCalled);
    {
      std::cout << __LINE__ << std::endl;
      auto h2 = h1;
      BOOST_CHECK_EQUAL(h1.refCount(), 2);
      BOOST_CHECK(!destructorCalled);

      {
        std::cout << __LINE__ << std::endl;
        auto h3 = h2;
        BOOST_CHECK_EQUAL(h1.refCount(), 3);
        BOOST_CHECK(!destructorCalled);

        {
          std::cout << __LINE__ << std::endl;
          auto h4 = std::move(h3);
          BOOST_CHECK(!h3);
          BOOST_CHECK_EQUAL(h1.refCount(), 3);
          BOOST_CHECK(!destructorCalled);
          std::cout << __LINE__ << std::endl;
        }
        BOOST_CHECK_EQUAL(h1.refCount(), 2);
        std::cout << __LINE__ << std::endl;
      }
      std::cout << __LINE__ << std::endl;
      BOOST_CHECK_EQUAL(h1.refCount(), 2);
      BOOST_CHECK(!destructorCalled);
      std::cout << __LINE__ << std::endl;
    }

    BOOST_CHECK_EQUAL(h1.refCount(), 1);
    BOOST_CHECK(!destructorCalled);
  }
  BOOST_CHECK(destructorCalled);
}

BOOST_AUTO_TEST_SUITE_END()
