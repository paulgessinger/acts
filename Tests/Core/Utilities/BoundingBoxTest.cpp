// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE BoundingBoxTests

#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <map>
#include <memory>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>
#include <set>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

struct Object {};

using ObjectBBox = Acts::AABB<Object>;

BOOST_AUTO_TEST_CASE(intersect_points)
{
  using vertex_type = ObjectBBox::vertex_type;

  Object o;
  ObjectBBox bb(o, {0, 0, 0}, {1, 1, 1});
  vertex_type p; 

  p = {0.5, 0.5, 0.5};
  BOOST_TEST(bb.intersect(p));
  p = {0.25, 0.25, 0.25};
  BOOST_TEST(bb.intersect(p));
  p = {0.75, 0.75, 0.75};
  BOOST_TEST(bb.intersect(p));

  // lower bound is inclusive
  p = {0, 0, 0};
  BOOST_TEST(bb.intersect(p));
  // upper bound is exclusive
  p = {1.0, 1.0, 1.0};
  BOOST_TEST(!bb.intersect(p));

  // some outsides
  p = {2, 0, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, 2, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, 0, 2};
  BOOST_TEST(!bb.intersect(p));
  p = {2, 2, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {2, 0, 2};
  BOOST_TEST(!bb.intersect(p));
  p = {2, 2, 2};
  BOOST_TEST(!bb.intersect(p));
  
  p = {-1, 0, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, -1, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, 0, -1};
  BOOST_TEST(!bb.intersect(p));
  p = {-1, -1, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {-1, 0, -1};
  BOOST_TEST(!bb.intersect(p));
  p = {-1, -1, -1};
  BOOST_TEST(!bb.intersect(p));

}

BOOST_AUTO_TEST_CASE(intersect_rays)
{
  BOOST_TEST_CONTEXT("2D") 
  {
    using Box = AxisAlignedBoundingBox<Object, float, 2>;

    Object o;
    Box bb(o, {-1, -1}, {1, 1});

    // ray in positive x direction

    Ray<float, 2> ray({-2, 0}, {1, 0});
    BOOST_TEST(bb.intersect(ray));

    ray = {{-2, 2}, {1, 0}};
    BOOST_TEST(!bb.intersect(ray));

    ray = {{-2, -2}, {1, 0}};
    BOOST_TEST(!bb.intersect(ray));

    // upper bound is exclusive
    ray = {{-2, 1}, {1, 0}};
    BOOST_TEST(!bb.intersect(ray));

    // lower bound is inclusive
    ray = {{-2, -1}, {1, 0}};
    BOOST_TEST(bb.intersect(ray));

    // ray faces away from box
    ray = {{2, 0}, {1, 0}};
    BOOST_TEST(!bb.intersect(ray));


    // ray in negative x direction 

    ray = {{2, 0}, {-1, 0}};
    BOOST_TEST(bb.intersect(ray));

    ray = {{2, 2}, {-1, 0}};
    BOOST_TEST(!bb.intersect(ray));

    ray = {{2, -2}, {-1, 0}};
    BOOST_TEST(!bb.intersect(ray));

    // upper bound is exclusive
    ray = {{2, 1}, {-1, 0}};
    BOOST_TEST(!bb.intersect(ray));

    // lower bound is inclusive
    ray = {{2, -1}, {-1, 0}};
    BOOST_TEST(bb.intersect(ray));


    // ray in positive y direction 

    ray = {{0, -2}, {0, 1}};
    BOOST_TEST(bb.intersect(ray));

    ray = {{2, -2}, {0, 1}};
    BOOST_TEST(!bb.intersect(ray));

    ray = {{-2, -2}, {0, 1}};
    BOOST_TEST(!bb.intersect(ray));

    // upper bound is exclusive
    ray = {{1, -2}, {0, 1}};
    BOOST_TEST(!bb.intersect(ray));

    // lower bound is not inclusive, 
    // due to Eigen's NaN handling.
    ray = {{-1, -2}, {0, 1}};
    BOOST_TEST(!bb.intersect(ray));

    // other direction
    ray = {{0, -2}, {0, -1}};
    BOOST_TEST(!bb.intersect(ray));


    // ray in positive y direction 

    ray = {{0, 2}, {0, -1}};
    BOOST_TEST(bb.intersect(ray));

    ray = {{2, 2}, {0, -1}};
    BOOST_TEST(!bb.intersect(ray));

    ray = {{-2, 2}, {0, -1}};
    BOOST_TEST(!bb.intersect(ray));

    // upper bound is exclusive
    ray = {{1, 2}, {0, -1}};
    BOOST_TEST(!bb.intersect(ray));

    // lower bound is not inclusive, 
    // due to Eigen's NaN handling.
    ray = {{-1, 2}, {0, -1}};
    BOOST_TEST(!bb.intersect(ray));
    
    // other direction
    ray = {{0, 2}, {0, 1}};
    BOOST_TEST(!bb.intersect(ray));


    // some off axis rays

    ray = {{-2, 0}, {0.5, 0.25}};
    BOOST_TEST(bb.intersect(ray));
    
    ray = {{-2, 0}, {0.5, 0.4}};
    BOOST_TEST(bb.intersect(ray));

    ray = {{-2, 0}, {0.5, 0.6}};
    BOOST_TEST(!bb.intersect(ray));
    
    ray = {{-2, 0}, {0.5, 0.1}};
    BOOST_TEST(bb.intersect(ray));
    
    ray = {{-2, 0}, {0.5, -0.4}};
    BOOST_TEST(bb.intersect(ray));
    
    ray = {{-2, 0}, {0.5, -0.6}};
    BOOST_TEST(!bb.intersect(ray));
    
    ray = {{-2, 0}, {0.1, 0.5}};
    BOOST_TEST(!bb.intersect(ray));

    // starting point inside
    ray = {{0, 0,}, {-1, 0}};
    BOOST_TEST(bb.intersect(ray));
    ray = {{0, 0,}, {1, 0}};
    BOOST_TEST(bb.intersect(ray));
    ray = {{0, 0,}, {0, -1}};
    BOOST_TEST(bb.intersect(ray));
    ray = {{0, 0,}, {0, 1}};
    BOOST_TEST(bb.intersect(ray));
  }

  BOOST_TEST_CONTEXT("3D")
  {
    using vertex_type3 = ObjectBBox::vertex_type;
    Object o;

    // lets make sure it also works in 3d
    ObjectBBox bb3(o, {-1, -1, -1}, {1, 1, 1});
    Ray<float, 3> ray3({0, 0, -2}, {0, 0, 1});
    BOOST_TEST(bb3.intersect(ray3));

    // facing away from box
    ray3 = {{0, 0, -2}, {0, 0, -1}};
    BOOST_TEST(!bb3.intersect(ray3));

    ray3 = {{0, 2, -2}, {0, 0, 1}};
    BOOST_TEST(!bb3.intersect(ray3));
    
    ray3 = {{0, -2, -2}, {0, 0, 1}};
    BOOST_TEST(!bb3.intersect(ray3));

    // right on slab
    ray3 = {{0, 1, -2}, {0, 0, 1}};
    BOOST_TEST(!bb3.intersect(ray3));

    // right on slab
    ray3 = {{0, -1, -2}, {0, 0, 1}};
    BOOST_TEST(bb3.intersect(ray3));

    // right on slab
    ray3 = {{-1, 0, -2 }, {0, 0, 1}};
    BOOST_TEST(!bb3.intersect(ray3));
    
    // right on slab
    ray3 = {{1, 0, -2 }, {0, 0, 1}};
    BOOST_TEST(!bb3.intersect(ray3));

    ray3 = {{-0.95, 0, -2 }, {0, 0, 1}};
    BOOST_TEST(bb3.intersect(ray3));
    
    // some off-axis rays
    ObjectBBox::vertex_type p(0, 0, -2);

    ray3 = {p, vertex_type3(1, 1, 1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {p, vertex_type3(-1, 1, 1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {p, vertex_type3(-1, -1, 1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {p, vertex_type3(1, -1, 1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {p, vertex_type3(1.1, 0, -1)-p};
    BOOST_TEST(!bb3.intersect(ray3));
    
    ray3 = {p, vertex_type3(-1.1, 0, -1)-p};
    BOOST_TEST(!bb3.intersect(ray3));

    ray3 = {p, vertex_type3(0, 1.1, -1)-p};
    BOOST_TEST(!bb3.intersect(ray3));

    ray3 = {p, vertex_type3(0, -1.1, -1)-p};
    BOOST_TEST(!bb3.intersect(ray3));

    ray3 = {p, vertex_type3(0.9, 0, -1)-p};
    BOOST_TEST(bb3.intersect(ray3));
    
    ray3 = {p, vertex_type3(-0.9, 0, -1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {p, vertex_type3(0, 0.9, -1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {p, vertex_type3(0, -0.9, -1)-p};
    BOOST_TEST(bb3.intersect(ray3));

    ray3 = {{0, 0, 0}, {1, 0, 0}};
    BOOST_TEST(bb3.intersect(ray3));
    ray3 = {{0, 0, 0}, {0, 1, 0}};
    BOOST_TEST(bb3.intersect(ray3));
    ray3 = {{0, 0, 0}, {0, 0, 1}};
    BOOST_TEST(bb3.intersect(ray3));
    
    ray3 = {{0, 0, 0}, {-1, 0, 0}};
    BOOST_TEST(bb3.intersect(ray3));
    ray3 = {{0, 0, 0}, {0, -1, 0}};
    BOOST_TEST(bb3.intersect(ray3));
    ray3 = {{0, 0, 0}, {0, 0, -1}};
    BOOST_TEST(bb3.intersect(ray3));
  }


}


BOOST_AUTO_TEST_CASE(frustum_intersect)
{
  BOOST_TEST_CONTEXT("2D")
  {

    using Frustum2 = Frustum<float, 2, 2>;
    Frustum2 fr({2, 2}, {1, 1}, M_PI/4.);

    std::ofstream os("frust2d.svg");
    fr.svg(os);
    os.close();

  }
  
  BOOST_TEST_CONTEXT("3D")
  {

    using Frustum3 = Frustum<float, 3, 3>;
    Frustum3 fr({2, 0, 0}, {1, 0, 1}, {0, 1, 0}, 3*M_PI/4.);
    
    std::ofstream os("frust3d.obj");
    size_t n_vtx = 1;
    fr.obj(os, n_vtx);
    os.close();



  }
}


}
}
