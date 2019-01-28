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

    // BEGIN VISUAL PARAMETER TEST

    std::ofstream os("frust2d.svg");
    float w = 1000;
    os << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    os << "<svg width=\"" << w << "\" height=\"" << w << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";

    size_t n=10;
    float min = -20, max = 20;
    float step = (max-min)/float(n);

    for(size_t i=0;i<=n;i++) {
      for(size_t j=0;j<=n;j++) {
        ActsVectorF<2> dir = {1, 0};
        ActsVectorF<2> origin = {min + step*i, min + step*j};
        origin.x() *= 1.10; // visual
        Eigen::Rotation2D<float> rot(2*M_PI/float(n) * i);
        float angle = 0.5*M_PI/n * j;
        Frustum2 fr(origin, rot * dir, angle);
        fr.svg(os, w, w, 2);
      }
    }

    os << "</svg>";
    os.close();

    // END VISUAL PARAMETER TEST

    os = std::ofstream("frust2d_test.svg");
    w = 1000;
    float unit = 20;
    os << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    os << "<svg width=\"" << w << "\" height=\"" << w << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";

    using Box = AxisAlignedBoundingBox<Object, float, 2>;
    Object o;
    Box::Size size(ActsVectorF<2>(2, 2));

    n = 10;
    float minx = -20;
    float miny = -20;
    float maxx = 20;
    float maxy = 20;
    float stepx = (maxx-minx)/float(n);
    float stepy = (maxy-miny)/float(n);

    std::vector<Box> boxes;
    boxes.reserve((n+1)*(n+1));
    for(size_t i=0;i<=n;i++) {
      for(size_t j=0;j<=n;j++) {
        boxes.emplace_back(o, ActsVectorF<2>{minx + i*stepx, miny + j*stepy}, size);
        std::stringstream ss;
        ss << boxes.size()-1;
        boxes.back().svg(os, w, w, unit, ss.str());
      }
    }

    Frustum2 fr({0, 0}, {1, 0}, M_PI/2.);
    fr.svg(os, w, w, 6, unit);
    
    os << "</svg>";

    os.close();
  }
  
  BOOST_TEST_CONTEXT("3D - 3 Sides")
  {

    using Frustum3 = Frustum<float, 3, 3>;
    size_t n_vtx = 1;
    //Frustum3 fr({2, 0, 0}, {1, 0, 1}, {0, 1, 0}, 1*M_PI/4.);
    auto make = [&](double angle, ActsVectorF<3> origin, std::ofstream& os) {
      float far = 1;
      Frustum3 fr(origin, {0, 0, 1}, angle);
      fr.obj(os, n_vtx, far);
      fr = Frustum3(origin, {0, 0, -1}, angle);
      fr.obj(os, n_vtx, far);
      fr = Frustum3(origin, {1, 0, 0}, angle);
      fr.obj(os, n_vtx, far);
      fr = Frustum3(origin, {-1, 0, 0}, angle);
      fr.obj(os, n_vtx, far);
      
      fr = Frustum3(origin, {0, 1, 0}, angle);
      fr.obj(os, n_vtx, far);
      fr = Frustum3(origin, {0, -1, 0}, angle);
      fr.obj(os, n_vtx, far);
    };

    std::ofstream os("frust3d_dir.obj");
    size_t s = 5;
    double min = -10, max = 10;
    double step = (max-min)/double(s);
    for(size_t i=0;i<=s;i++) {
      for(size_t j=0;j<=s;j++) {
        for(size_t k=0;k<=s;k++) {
          ActsVectorF<3> origin(min + i*step, min+j*step, min+k*step);
          //std::cout << origin.transpose() << std::endl;
          make(M_PI/4., origin, os);
        }
      }
    }
    os.close();

    os = std::ofstream("frust3D_angle.obj");
    n_vtx = 1;
    size_t n = 10;
    Eigen::Affine3f rot;
    for(size_t i=0;i<=n;i++) {
      ActsVectorF<3> origin(i*4, 0, 0);
      rot = Eigen::AngleAxisf(M_PI/float(n)*i, ActsVectorF<3>::UnitY());
      float angle = (M_PI/2.)/float(n) * (1+i);
      ActsVectorF<3> dir(1, 0, 0);
      Frustum3 fr(origin, rot*dir, angle);
      fr.obj(os, n_vtx, 2);

    }

    os.close();

  }
  
  BOOST_TEST_CONTEXT("3D - 4 Sides")
  {
    using Frustum34 = Frustum<float, 3, 4>;
    size_t n_vtx = 1;
    
    std::ofstream os("frust3d-4s_dir.obj");

    size_t s = 5;
    double min = -10, max = 10;
    double step = (max-min)/double(s);
    double angle = M_PI/4.;
    for(size_t i=0;i<=s;i++) {
      for(size_t j=0;j<=s;j++) {
        for(size_t k=0;k<=s;k++) {
          ActsVectorF<3> origin(min + i*step, min+j*step, min+k*step);
          ActsVectorF<3> dir(1, 0, 0);

          Eigen::Affine3f rot;
          rot = Eigen::AngleAxisf(M_PI/float(s)*i, ActsVectorF<3>::UnitX())
              * Eigen::AngleAxisf(M_PI/float(s)*j, ActsVectorF<3>::UnitY())
              * Eigen::AngleAxisf(M_PI/float(s)*k, ActsVectorF<3>::UnitZ());

          Frustum34 fr(origin, rot*dir, angle);
          fr.obj(os, n_vtx, 1);
        }
      }
    }
    os.close();
    os = std::ofstream("frust3d-4s_angle.obj");

    n_vtx = 1;
    size_t n = 10;
    for(size_t i=0;i<=n;i++) {
      ActsVectorF<3> origin(i*4, 0, 0);
      Eigen::Affine3f rot;
      rot = Eigen::AngleAxisf(M_PI/float(n)*i, ActsVectorF<3>::UnitY());
      angle = (M_PI/2.)/float(n) * (1+i);
      ActsVectorF<3> dir(1, 0, 0);
      Frustum34 fr(origin, rot*dir, angle);
      fr.obj(os, n_vtx, 2);

    }

    os.close();
  }

}


}
}
