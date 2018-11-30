#include <iostream>
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

struct Object {};

using ObjectBBox = Acts::AABB<Object>;
using vertex_t = ObjectBBox::vertex_t;


void test() 
{
  std::cout << "hello world" << std::endl;
  Object o;

  vertex_t vmin(0, 0, 0);
  vertex_t vmax(10, 10, 10);

  ObjectBBox bbox(o, vmin, vmax);
  std::cout << bbox << std::endl;


  ObjectBBox bb1(o, {0, 0, 0}, {1, 1, 1});
  ObjectBBox bb2(o, {0, 0, 1}, {1, 1, 2});

  ObjectBBox bbcomp({bb1, bb2});
  std::cout << bbcomp << std::endl;

  ObjectBBox bb3(o, {0, 2, 0}, {2, 2, 1});
  ObjectBBox bbcomp2({bbcomp, bb3});

  std::cout << bbcomp2 << std::endl;


  auto nodes = bbcomp2.flatten();
  for(const auto& node : nodes) {
    std::cout << node << std::endl;
  }


}

}
}

int main() {

  Acts::Test::test();

  return 0;
}
