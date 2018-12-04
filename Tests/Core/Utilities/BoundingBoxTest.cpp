#include <iostream>
#include <map>
#include <memory>
#include <fstream>
#include <random>
#include <chrono>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"


namespace Acts {
namespace Test {

struct Object {};

using ObjectBBox = Acts::AABB<Object>;
using vertex_type = ObjectBBox::vertex_type;
using vertex_array_type = ObjectBBox::vertex_array_type;
using value_type = ObjectBBox::value_type;


//void test_general() 
//{
  //std::cout << "hello world" << std::endl;
  //Object o;

  //vertex_type vmin(0, 0, 0);
  //vertex_type vmax(10, 10, 10);

  //ObjectBBox bbox(o, vmin, vmax);
  //bbox.name = "bbox";
  //std::cout << bbox << std::endl;


  //ObjectBBox bb1(o, {0, 0, 0}, {1, 1, 1});
  //ObjectBBox bb2(o, {0, 0, 1}, {1, 1, 2});

  //bb1.name = "bb1";
  //bb2.name = "bb2";
  
  //ObjectBBox bb4(o, {-3, -1, 0}, {2, 0, 1});
  //bb4.name = "bb4";
  //ObjectBBox bbcomp3({&bb2,&bb4});
  //bbcomp3.name = "bbcomp3";

  //ObjectBBox bbcomp({&bb1, &bbcomp3});
  //bbcomp.name = "bbcomp";


  //std::cout << bbcomp << std::endl;

  //ObjectBBox bb3(o, {0, 2, 0}, {2, 2, 1});
  //bb3.name = "bb3";
  //ObjectBBox bbcomp2({&bbcomp, &bb3});
  //bbcomp2.name = "bbcomp2";
  //std::cout << bbcomp2 << std::endl;

  ////auto nodes = bbcomp2.flatten();
  ////for(const auto& node : nodes) {
    ////std::cout << node << std::endl;
  ////}
  
  //std::cout << "----" << std::endl;

  //std::map<const ObjectBBox*, bool> hitmap;
  //hitmap[&bb1] = true;
  //hitmap[&bb2] = true;
  //hitmap[&bb3] = true;
  //hitmap[&bb4] = true;
  //hitmap[&bbcomp] = true;
  //hitmap[&bbcomp2] = true;
  //hitmap[&bbcomp3] = true;

  //const auto* node = &bbcomp2;
  //do {
    //std::cout << *node << std::endl;
    //if(hitmap[node] == true) {
      //if(node->hasEntity()) {
        //// do entity intersect check or collect entity
        //node = node->getSkip();
      //}
      //else {
        //node = node->getLeftChild();
      //}
    //}
    //else {
      //node = node->getSkip();
    //}
  //}
  //while(node != nullptr);


//}

void
make_grid(std::vector<std::unique_ptr<ObjectBBox>>& boxes,
          vertex_type             vmin,
          vertex_type             vmax,
          vertex_array_type       n,
          vertex_type             hw)
{
  Object      o;
  vertex_type step = (vmax - vmin).array() / (n - 1);
  // vertex_type hw = step / 2.;

  for (size_t i = 0; i < n[0]; i++) {
    for (size_t j = 0; j < n[1]; j++) {
      for (size_t k = 0; k < n[2]; k++) {
        vertex_type v(step.x() * i, step.y() * j, step.z() * k);

        vertex_type bmin = vmin + v - hw;
        vertex_type bmax = vmin + v + hw;
        //std::cout << bmin.x() << ", " << bmin.y() << ", " << bmin.z() <<
        //std::endl; std::cout << bmax.x() << ", " << bmax.y() << ", " <<
        //bmax.z() << std::endl; std::cout << std::endl;

        boxes.emplace_back(std::make_unique<ObjectBBox>(o, bmin, bmax));
      }
    }
  }
}

std::array<std::vector<ObjectBBox*>, 8>
make_octants(const std::vector<ObjectBBox*>& boxes) {

  std::array<std::vector<ObjectBBox*>, 8> octants;
  // calc center of boxes
  //ObjectBBox wrap = ObjectBBox::wrap(boxes);
  ObjectBBox::vertex_type vmin, vmax;
  std::tie(vmin, vmax) = ObjectBBox::wrap(boxes);
  ObjectBBox::vertex_type glob_ctr = (vmin + vmax)/2.;

  //std::cout << vmin << std::endl << vmax << std::endl;

  for(auto* box : boxes) {
    vertex_type ctr = box->center() - glob_ctr;
    if(ctr.x() < 0 && ctr.y() < 0 && ctr.z() < 0) {octants[0].push_back(box);}
    if(ctr.x() > 0 && ctr.y() < 0 && ctr.z() < 0) {octants[1].push_back(box);}
    if(ctr.x() < 0 && ctr.y() > 0 && ctr.z() < 0) {octants[2].push_back(box);}
    if(ctr.x() > 0 && ctr.y() > 0 && ctr.z() < 0) {octants[3].push_back(box);}

    if(ctr.x() < 0 && ctr.y() < 0 && ctr.z() > 0) {octants[4].push_back(box);}
    if(ctr.x() > 0 && ctr.y() < 0 && ctr.z() > 0) {octants[5].push_back(box);}
    if(ctr.x() < 0 && ctr.y() > 0 && ctr.z() > 0) {octants[6].push_back(box);}
    if(ctr.x() > 0 && ctr.y() > 0 && ctr.z() > 0) {octants[7].push_back(box);}
  }

  return octants;
}


//std::tuple<std::vector<ObjectBBox>, std::vector<ObjectBBox>, ObjectBBox*>
std::unique_ptr<ObjectBBox>
make_geometry(std::vector<std::unique_ptr<ObjectBBox>>& boxes, std::vector<std::unique_ptr<ObjectBBox>>& octant_boxes) {

  //std::vector<ObjectBBox> boxes;
  make_grid(boxes, {-10, -10, -10}, {10, 10, 10}, {8, 8, 8}, {1, 1, 1});
  //std::cout << __LINE__ << " boxes.size() = " << boxes.size() << std::endl;

  std::vector<ObjectBBox*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  std::transform(boxes.begin(),
                 boxes.end(),
                 std::back_inserter(box_ptrs),
                 [](auto& box) { return box.get(); });
  
  std::ofstream prims("prims.obj");
  size_t n_vtx = 1;
  
  for(auto& box : boxes) {
    box->obj(prims, n_vtx);
  }
  
  prims.close();
  size_t n_vtx_oct1 = 1;
  size_t n_vtx_oct2 = 1;
  std::ofstream oct1("oct1.obj");
  std::ofstream oct2("oct2.obj");


  vertex_array_type envelope = vertex_array_type::Constant(0.1);

  std::array<std::vector<ObjectBBox*>, 8> octants;
  //std::cout << "box_ptrs.size() = " << box_ptrs.size() << std::endl;
  octants = make_octants(box_ptrs);
  //std::vector<ObjectBBox> octant_boxes;
  std::vector<ObjectBBox*> octant_boxes_ptr;
  octant_boxes.reserve(8*8);

  for(const auto& octant : octants) {
    //std::cout << octant.size() << std::endl;

    std::array<std::vector<ObjectBBox*>, 8> sub_octants;
    sub_octants = make_octants(octant);

    std::vector<ObjectBBox*> sub_octant_boxes;

    //std::cout << " - ";
    for(auto& sub_octant: sub_octants) {
      //std::cout << " " << sub_octant.size();;
      octant_boxes.emplace_back(std::make_unique<ObjectBBox>(sub_octant, envelope));
      sub_octant_boxes.push_back(octant_boxes.back().get());
      sub_octant_boxes.back()->obj(oct2, n_vtx_oct2);
    }
    //std::cout << std::endl;

    octant_boxes.emplace_back(std::make_unique<ObjectBBox>(sub_octant_boxes, envelope));
    octant_boxes_ptr.push_back(octant_boxes.back().get());
    octant_boxes.back()->obj(oct1, n_vtx_oct1);
  }


  auto world = std::make_unique<ObjectBBox>(octant_boxes_ptr);
  const ObjectBBox* node = world.get();
  size_t objects_found = 0;
  do {
    //std::cout << *node << std::endl;
    if(node->hasEntity()) {
      // do entity intersect check or collect entity
      node = node->getSkip();
      objects_found++;
    }
    else {
      node = node->getLeftChild();
    }
  }
  while(node != nullptr);


  assert(boxes.size() == objects_found);
  oct1.close();

  return world;
}

void test_points() {

  std::vector<std::unique_ptr<ObjectBBox>> boxes;
  std::vector<std::unique_ptr<ObjectBBox>> octants;
  std::unique_ptr<ObjectBBox> world = make_geometry(boxes, octants);

  std::mt19937 rng(42);
  std::uniform_real_distribution<value_type> dist;

  std::vector<vertex_type> points;
  
  for(const auto& box : boxes) {
    points.emplace_back(box->center());
  }

  std::shuffle(points.begin(), points.end(), rng);
  

  auto bench = [&boxes, &world](const std::vector<vertex_type>& points) -> std::tuple<size_t, size_t, double, double> {
    size_t n_intersects = 0;

    auto dumb_search = [&](const vertex_type& v) -> const ObjectBBox* {
      for(const auto& bb : boxes) {
        n_intersects++;
        if(bb->intersect(v)) {
          return bb.get();
        }
      }
      return nullptr;
    };

    auto bvh_search = [&](const vertex_type& v) -> const ObjectBBox* {
      const ObjectBBox* node = world.get();
      do {
        n_intersects++;
        if(node->intersect(v)) {
          
          if(node->hasEntity()) {
            // found primitive
            //hit = node;
            //break;
            return node;
          }
          else {
            // go over children
            node = node->getLeftChild();
          }
        }
        else {
          node = node->getSkip();
        }
      }
      while(node != nullptr);
      return nullptr;
    };

    std::vector<const ObjectBBox*> dumb_hits;
    dumb_hits.reserve(points.size());

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // dumb search algorithm
    for(const auto& v : points) {
      dumb_hits.push_back(dumb_search(v));
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double diff_dumb = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000.;
    std::cout << "Time difference = " << diff_dumb << "s, " << n_intersects << " intersects" <<std::endl;

    size_t n_int_dumb = n_intersects;

    std::vector<const ObjectBBox*> bvh_hits;
    bvh_hits.reserve(points.size());
    n_intersects = 0;

    begin = std::chrono::steady_clock::now();

    for(const auto& v : points) {
      bvh_hits.push_back(bvh_search(v));
    }

    end = std::chrono::steady_clock::now();
    double diff_bvh = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000.;
    std::cout << "Time difference = " << diff_bvh << "s, " << n_intersects << " intersects" <<std::endl;
    
    size_t n_int_bvh = n_intersects;

    for(size_t i=0;i<points.size();i++) {
      auto dumb = dumb_hits.at(i);
      auto bvh = bvh_hits.at(i);
      //std::cout << "dumb: " << dumb << ", bvh: " << bvh << std::endl;
      assert(dumb == bvh);
    }

    return {n_int_dumb, n_int_bvh, diff_dumb, diff_bvh};

  };

  bench(points);

  //std::vector<size_t> sizes = {10, 100, 1000, 10000, 100000};
  size_t powers = 6;
  std::vector<size_t> sizes;
  for(size_t i=0;i<powers;i++) {
    size_t n = 1;
    for(size_t j=0;j<=i;j++) {
      n *= 10;
    }
    sizes.push_back(n);
    std::cout << n << std::endl;
  }

  std::ofstream bch("times.csv");
  bch << "n,n_int_dumb,n_int_bvh,t_dumb,t_bvh" << std::endl;

  for(size_t n : sizes) {
    //size_t n = 1e5;
    points.resize(0);
    points.reserve(n);
    for(size_t i=0;i<n;i++) {
      points.emplace_back(dist(rng), dist(rng), dist(rng));
    }
    size_t n_int_dumb, n_int_bvh;
    double diff_dumb, diff_bvh;
    std::tie(n_int_dumb, n_int_bvh, diff_dumb, diff_bvh) = bench(points);
    bch << n << "," << n_int_dumb << "," << n_int_bvh << "," << diff_dumb << "," << diff_bvh << std::endl;
  }
  bch.close();


  
}

void test_intersect_points()
{
  Object o;
  ObjectBBox bb(o, {0, 0, 0}, {1, 1, 1});
  vertex_type p; 

  p = {0.5, 0.5, 0.5};
  assert(bb.intersect(p));
  p = {0.25, 0.25, 0.25};
  assert(bb.intersect(p));
  p = {0.75, 0.75, 0.75};
  assert(bb.intersect(p));

  // lower bound is inclusive
  p = {0, 0, 0};
  assert(bb.intersect(p));
  // upper bound is exclusive
  p = {1.0, 1.0, 1.0};
  assert(!bb.intersect(p));

  // some outsides
  p = {2, 0, 0};
  assert(!bb.intersect(p));
  p = {0, 2, 0};
  assert(!bb.intersect(p));
  p = {0, 0, 2};
  assert(!bb.intersect(p));
  p = {2, 2, 0};
  assert(!bb.intersect(p));
  p = {2, 0, 2};
  assert(!bb.intersect(p));
  p = {2, 2, 2};
  assert(!bb.intersect(p));
  
  p = {-1, 0, 0};
  assert(!bb.intersect(p));
  p = {0, -1, 0};
  assert(!bb.intersect(p));
  p = {0, 0, -1};
  assert(!bb.intersect(p));
  p = {-1, -1, 0};
  assert(!bb.intersect(p));
  p = {-1, 0, -1};
  assert(!bb.intersect(p));
  p = {-1, -1, -1};
  assert(!bb.intersect(p));
}

void test_intersect_rays()
{

  using Box = AxisAlignedBoundingBox<Object, float, 2>;
  Object o;
  Box bb(o, {-1, -1}, {1, 1});

  // ray in positive x direction

  Ray<float, 2> ray({-2, 0}, {1, 0});
  assert(bb.intersect(ray));

  ray = {{-2, 2}, {1, 0}};
  assert(!bb.intersect(ray));

  ray = {{-2, -2}, {1, 0}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{-2, 1}, {1, 0}};
  assert(!bb.intersect(ray));

  // lower bound is inclusive
  ray = {{-2, -1}, {1, 0}};
  assert(bb.intersect(ray));

  // ray faces away from box
  ray = {{2, 0}, {1, 0}};
  assert(!bb.intersect(ray));


  // ray in negative x direction 

  ray = {{2, 0}, {-1, 0}};
  assert(bb.intersect(ray));

  ray = {{2, 2}, {-1, 0}};
  assert(!bb.intersect(ray));

  ray = {{2, -2}, {-1, 0}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{2, 1}, {-1, 0}};
  assert(!bb.intersect(ray));

  // lower bound is inclusive
  ray = {{2, -1}, {-1, 0}};
  assert(bb.intersect(ray));


  // ray in positive y direction 

  ray = {{0, -2}, {0, 1}};
  assert(bb.intersect(ray));

  ray = {{2, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  ray = {{-2, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{1, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  // lower bound is not inclusive, 
  // due to Eigen's NaN handling.
  ray = {{-1, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  // other direction
  ray = {{0, -2}, {0, -1}};
  assert(!bb.intersect(ray));


  // ray in positive y direction 

  ray = {{0, 2}, {0, -1}};
  assert(bb.intersect(ray));

  ray = {{2, 2}, {0, -1}};
  assert(!bb.intersect(ray));

  ray = {{-2, 2}, {0, -1}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{1, 2}, {0, -1}};
  assert(!bb.intersect(ray));

  // lower bound is not inclusive, 
  // due to Eigen's NaN handling.
  ray = {{-1, 2}, {0, -1}};
  assert(!bb.intersect(ray));
  
  // other direction
  ray = {{0, 2}, {0, 1}};
  assert(!bb.intersect(ray));


  // some off axis rays

  ray = {{-2, 0}, {0.5, 0.25}};
  assert(bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, 0.4}};
  assert(bb.intersect(ray));

  ray = {{-2, 0}, {0.5, 0.6}};
  assert(!bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, 0.1}};
  assert(bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, -0.4}};
  assert(bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, -0.6}};
  assert(!bb.intersect(ray));
  
  ray = {{-2, 0}, {0.1, 0.5}};
  assert(!bb.intersect(ray));


  // lets make sure it also works in 3d
  ObjectBBox bb3(o, {-1, -1, -1}, {1, 1, 1});
  Ray<float, 3> ray3({0, 0, -2}, {0, 0, 1});
  assert(bb3.intersect(ray3));

  ray3 = {{0, 0, -2}, {0, 0, -1}};
  assert(!bb3.intersect(ray3));
  
  ray3 = {{0, 2, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));
  
  ray3 = {{0, -2, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  ray3 = {{0, 1, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  ray3 = {{0, -1, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  ray3 = {{-1, 0, -2 }, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  ray3 = {{-0.95, 0, -2 }, {0, 0, 1}};
  assert(bb3.intersect(ray3));
  
  // some off-axis rays
  ObjectBBox::vertex_type p(0, 0, -2);

  ray3 = {p, vertex_type(1, 1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type(-1, 1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type(-1, -1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type(1, -1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type(1.1, 0, -1)-p};
  assert(!bb3.intersect(ray3));
  
  ray3 = {p, vertex_type(-1.1, 0, -1)-p};
  assert(!bb3.intersect(ray3));

  ray3 = {p, vertex_type(0, 1.1, -1)-p};
  assert(!bb3.intersect(ray3));

  ray3 = {p, vertex_type(0, -1.1, -1)-p};
  assert(!bb3.intersect(ray3));

  ray3 = {p, vertex_type(0.9, 0, -1)-p};
  assert(bb3.intersect(ray3));
  
  ray3 = {p, vertex_type(-0.9, 0, -1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type(0, 0.9, -1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type(0, -0.9, -1)-p};
  assert(bb3.intersect(ray3));

}


}
}


int main() {

  Acts::Test::test_intersect_points();
  Acts::Test::test_intersect_rays();

  Acts::Test::test_points();
  //Acts::Test::test_rays();
  

  return 0;
}
