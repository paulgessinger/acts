#include <iostream>
#include <map>
#include <memory>
#include <fstream>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

struct Object {};

using ObjectBBox = Acts::AABB<Object>;
using vertex_t = ObjectBBox::vertex_t;
using vertex_array_t = ObjectBBox::vertex_array_t;


//void test_general() 
//{
  //std::cout << "hello world" << std::endl;
  //Object o;

  //vertex_t vmin(0, 0, 0);
  //vertex_t vmax(10, 10, 10);

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

std::pair<std::vector<ObjectBBox>, ObjectBBox*>
make_grid(vertex_t vmin, vertex_t vmax, Eigen::Array3f n, vertex_t hw) {
  Object o;
  vertex_t step = (vmax-vmin).array() / (n-1);
  //vertex_t hw = step / 2.;

  std::vector<ObjectBBox> boxes;


  for(size_t i=0;i<n[0];i++) {
    for(size_t j=0;j<n[1];j++) {
      for(size_t k=0;k<n[2];k++) {
        vertex_t v(step.x()*i, step.y()*j, step.z()*k);

        vertex_t bmin = vmin + v - hw;
        vertex_t bmax = vmin + v + hw;
        //std::cout << bmin.x() << ", " << bmin.y() << ", " << bmin.z() << std::endl;
        //std::cout << bmax.x() << ", " << bmax.y() << ", " << bmax.z() << std::endl;
        //std::cout << std::endl;

        boxes.emplace_back(o, bmin, bmax);
      }
    }
  }


  return {std::move(boxes), nullptr};
}

//std::vector<ObjectBBox*> make_octant(std::vector<ObjectBBox>& boxes) {

  //// calc center of boxes
  //ObjectBBox wrap = ObjectBBox::wrap(boxes);

//}

void test_points() {

  std::vector<ObjectBBox> boxes;
  ObjectBBox* top;
  std::tie(boxes, top) = make_grid({-10, -10, -10}, {10, 10, 10}, {8, 8, 8}, {1, 1, 1});
  
  std::ofstream prims("prims.obj");
  size_t n_vtx = 1;

  using octant_t = std::vector<ObjectBBox*>;
  std::array<octant_t, 8> octants;

  std::cout << "boxes.size() = " << boxes.size() << std::endl;

  for(auto& box : boxes) {
    vertex_t ctr = box.center();
    box.obj(prims, n_vtx);
    //std::cout << ctr << "\n" << std::endl;
    if(ctr.x() < 0 && ctr.y() < 0 && ctr.z() < 0) {octants[0].push_back(&box);}
    if(ctr.x() > 0 && ctr.y() < 0 && ctr.z() < 0) {octants[1].push_back(&box);}
    if(ctr.x() < 0 && ctr.y() > 0 && ctr.z() < 0) {octants[2].push_back(&box);}
    if(ctr.x() > 0 && ctr.y() > 0 && ctr.z() < 0) {octants[3].push_back(&box);}

    if(ctr.x() < 0 && ctr.y() < 0 && ctr.z() > 0) {octants[4].push_back(&box);}
    if(ctr.x() > 0 && ctr.y() < 0 && ctr.z() > 0) {octants[5].push_back(&box);}
    if(ctr.x() < 0 && ctr.y() > 0 && ctr.z() > 0) {octants[6].push_back(&box);}
    if(ctr.x() > 0 && ctr.y() > 0 && ctr.z() > 0) {octants[7].push_back(&box);}
  }

  vertex_array_t envelope = vertex_array_t::Constant(0.1);

  prims.close();
  n_vtx = 1;
  std::ofstream oct1("oct1.obj");

  std::vector<std::unique_ptr<ObjectBBox>> octant_boxes;
  std::vector<ObjectBBox*> octant_boxes_ptr;
  for(auto& octant : octants) {
    std::cout << octant.size() << std::endl;
    //assert(octant.size() == 8);
    auto octbox = std::make_unique<ObjectBBox>(ObjectBBox(octant, envelope));
    octbox->obj(oct1, n_vtx);
    octant_boxes_ptr.push_back(octbox.get());
    octant_boxes.push_back(std::move(octbox));
  }

  ObjectBBox world(octant_boxes_ptr);
  const ObjectBBox* node = &world;
  do {
    std::cout << *node << std::endl;
    if(node->hasEntity()) {
      // do entity intersect check or collect entity
      node = node->getSkip();
    }
    else {
      node = node->getLeftChild();
    }
  }
  while(node != nullptr);

  //std::cout << world << std::endl;

  oct1.close();

  
  //Object o;
  //ObjectBBox box(o, {0, 0, 0}, {10, 10, 10});
  //std::ofstream file("grid.obj");
  //size_t n=1;
  //box.obj(file, n);

}

}
}

int main() {

  //Acts::Test::test_general();
  Acts::Test::test_points();

  return 0;
}
