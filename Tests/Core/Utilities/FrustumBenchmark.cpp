#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  struct Object
  {
  };

  using real_t = float;
  using vec_t = ActsVector<real_t, 3>;
  using Box = Acts::AxisAlignedBoundingBox<Object, real_t, 3>;
  using clock = std::chrono::steady_clock;

  template <typename T>
  double to_s(const T& delta) {
    return std::chrono::duration_cast<std::chrono::microseconds>(delta).count();
  }

  template <size_t sides>
  std::tuple<size_t, double>
  bench_frustum(size_t octree_depth)
  {
    std::cout << "Bench frustum: sides: " << sides << " octd: " << octree_depth << std::endl;
    using Frustum = Frustum<real_t, 3, sides>;

    size_t n    = 19;
    real_t  min  = -100;
    real_t  max  = 100;
    real_t  step = (max - min) / real_t(n);
    
    ply_helper<real_t> ply;

    Object o;
    Box::Size size(ActsVectorF<3>(2, 2, 2));

    std::vector<std::unique_ptr<Box>> boxes;
    boxes.reserve((n+1)*(n+1)*(n+1));
    std::vector<Box*> prims;
    prims.reserve((n+1)*(n+1)*(n+1));

    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          ActsVectorF<3> pos(min + i * step, min + j * step, min + k * step);
          //boxes.emplace_back(o, pos, size);
          boxes.push_back(std::make_unique<Box>(o, pos, size));
          //boxes.back()->draw(ply);
          prims.push_back(boxes.back().get());
        }
      }
    }

    //std::vector<std::unique_ptr<Box>> box_store;
    //std::vector<Box*> prims;
    //prims.reserve(boxes.size());
    //for(Box& bb : boxes) {
      //prims.push_back(&bb);
    //}

    Box* top_node = make_octree(boxes, prims, octree_depth);
    std::mt19937 rng(42);
    std::uniform_real_distribution<real_t> pos_dist(min*1.5, max*1.5);
    std::uniform_real_distribution<real_t> dir_dist(-1, 1);
    std::uniform_real_distribution<real_t> angle_dist(0.1*M_PI, M_PI/2.*0.9);

    size_t n_frust = 1e6;


    std::vector<Frustum> frustums;
    frustums.reserve(n_frust);

    for(size_t i=0;i<n_frust;i++) {
      vec_t pos(pos_dist(rng), pos_dist(rng), pos_dist(rng));
      vec_t dir(dir_dist(rng), dir_dist(rng), dir_dist(rng));
      dir.normalize();
      real_t angle = angle_dist(rng);

      //std::cout << "pos: " << pos.transpose() << std::endl;
      //std::cout << "dir: " << dir.transpose() << std::endl;
      //std::cout << "angle: " << angle << std::endl;

      frustums.emplace_back(pos, dir, angle);
      if(i<50) {
        frustums.back().draw(ply);
      }
    }

    std::ofstream os("prims.ply");

    for(Box* box : prims) {
      box->draw(ply);
    }

    os << ply << std::flush;
    os.close();

    ply.clear();
    for(auto& box : boxes) {
      box->draw(ply);
    }

    os = std::ofstream("octree.ply");
    os << ply << std::flush;
    os.close();


    clock::time_point begin = clock::now();

    //size_t n_ix = frustums.size() * boxes.size();
    size_t n_ix = 0;

    //for(const Frustum& fr : frustums) {
      //for(const Box* bb : prims) {
        //volatile bool ix = bb->intersect(fr);
      //}
    //}

    for(const Frustum& fr : frustums) {
      std::vector<const Box*> result;
      const Box* lnode = top_node;
      do {
        n_ix++;
        if (lnode->intersect(fr)) {

          if (lnode->hasEntity()) {
            // found primitive
            result.push_back(lnode);
            //std::cout << "intersect prim: " << lnode << " -> " << lnode->getSkip() << std::endl;
            lnode = lnode->getSkip();
          } else {
            // go over children
            //std::cout << "intersect bvh: " << lnode << " -> " << lnode->getLeftChild() << std::endl;
            lnode = lnode->getLeftChild();
          }
        } else {
          //std::cout << "no intersect: " << lnode << " -> " << lnode->getSkip() << std::endl;
          lnode = lnode->getSkip();
        }

        //if(n_intersects > 80) throw std::runtime_error("");
      } while (lnode != nullptr);
    }

    clock::time_point end = clock::now();
    double diff = to_s(end - begin);
    std::cout << n_ix << " intersects, " << diff/1000000. << "s" << std::endl;
    double t_per_ix = diff / n_ix;
    std::cout << "=> " << t_per_ix << "us per intersect" << std::endl;
    double ix_per_sec = n_ix / (diff/1000000.);
    std::cout << "=> " << ix_per_sec << " intersects per second" << std::endl;


    return {n_ix, diff};
  }

}  // namespace Test
}  // namespace Acts

int
main()
{
  using res_t = std::tuple<size_t, double>;

  std::vector<std::pair<size_t, res_t>> results;

  // Tested with sides = 10, NaN / inf behavior!
  
  //results = {
    //{3, Acts::Test::bench_frustum<3>()},
    //{8, Acts::Test::bench_frustum<8>()},
    //{9, Acts::Test::bench_frustum<9>()},
    //{10, Acts::Test::bench_frustum<10>()},
    //{11, Acts::Test::bench_frustum<11>()},
    //{16, Acts::Test::bench_frustum<16>()},
    //{19, Acts::Test::bench_frustum<19>()},
    //{20, Acts::Test::bench_frustum<20>()},
    //{21, Acts::Test::bench_frustum<21>()},
    //{100, Acts::Test::bench_frustum<100>()}
  //};

  //std::ofstream os("sides.csv");
  //os << "sides,n_ix,diff\n";
  //for(const auto& pr : results) {
    //size_t sides = pr.first;
    //size_t n_ix;
    //double diff;
    //std::tie(n_ix, diff) = pr.second;
    //os << sides << "," << n_ix << "," << diff << "\n";
    //sides++;
  //}
  //os.close();
  
  std::ofstream os("octree_fr.csv");
  os << "octd,nix,time\n";
  std::vector<size_t> octree_depths = {0, 1, 2, 3, 4, 5};
  for(size_t octree_depth : octree_depths) {
    size_t n_ix;
    double diff;
    std::tie(n_ix, diff) = Acts::Test::bench_frustum<9>(octree_depth);
    os << octree_depth << "," << n_ix << "," << diff << "\n";
  }
  os.close();

  return 0;
}
