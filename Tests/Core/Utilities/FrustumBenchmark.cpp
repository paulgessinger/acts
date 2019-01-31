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

#define ASSERT(condition) \
  if(!(condition)) {\
    std::cerr << "Assertion failed in line #" << __LINE__ << std::endl; \
    abort(); \
  }

namespace Acts {
namespace Test {

  struct Object
  {
    using id_t = std::tuple<size_t, size_t, size_t>;
    Object(id_t id_) : id(id_) {}
    id_t id;
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
  std::tuple<size_t, double, std::vector<std::set<Object::id_t>>>
  bench_frustum(int octree_depth)
  {
    std::cout << "Bench frustum: sides: " << sides << " octd: " << octree_depth << std::endl;
    using Frustum = Frustum<real_t, 3, sides>;

    size_t n    = 19;
    real_t  min  = -100;
    real_t  max  = 100;
    real_t  step = (max - min) / real_t(n);
    
    ply_helper<real_t> ply;

    Box::Size size(ActsVectorF<3>(2, 2, 2));

    std::vector<std::unique_ptr<Box>> boxes;
    boxes.reserve((n+1)*(n+1)*(n+1));
    std::vector<Box*> prims;
    prims.reserve((n+1)*(n+1)*(n+1));
    std::vector<std::unique_ptr<Object>> objects;

    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          ActsVectorF<3> pos(min + i * step, min + j * step, min + k * step);
          //boxes.emplace_back(o, pos, size);
          //Object o{{i, j, k}};
          //Object o;
          //o.id = {i, j, k};
          objects.push_back(std::make_unique<Object>(std::make_tuple(i ,j, k)));
          boxes.push_back(std::make_unique<Box>(*objects.back(), pos, size));
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

    Box* top_node = nullptr;

    if(octree_depth >= 0) {
      top_node = make_octree(boxes, prims, octree_depth);
    }
    else {
      // no octree, just link all boxes sequentially
      for(size_t i=1;i<prims.size();i++) {
        prims.at(i-1)->setSkip(prims.at(i));
      }
      top_node = prims.front();
    }



    std::mt19937 rng(42);
    std::uniform_real_distribution<real_t> pos_dist(min*1.5, max*1.5);
    std::uniform_real_distribution<real_t> dir_dist(-1, 1);
    std::uniform_real_distribution<real_t> angle_dist(0.1*M_PI, M_PI/2.*0.9);

    size_t n_frust = 1e4;

    std::cout << "Generating " << n_frust << " frustums" << std::endl;


    std::vector<Frustum, Eigen::aligned_allocator<Frustum>> frustums;
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


    std::cout << "begin intersect timing" << std::endl;
    clock::time_point begin = clock::now();

    //size_t n_ix = frustums.size() * boxes.size();
    size_t n_ix = 0;

    //for(const Frustum& fr : frustums) {
      //for(const Box* bb : prims) {
        //volatile bool ix = bb->intersect(fr);
      //}
    //}

    std::vector<std::set<Object::id_t>> result;

    for(const Frustum& fr : frustums) {
      const Box* lnode = top_node;
      std::set<Object::id_t> hits;
      do {
        n_ix++;
        if (lnode->intersect(fr)) {

          if (lnode->hasEntity()) {
            // found primitive
            hits.insert(lnode->entity()->id);
            auto id = lnode->entity()->id;
            //std::cout << "intersect prim: idx: (" << std::get<0>(id) << ", " << std::get<1>(id) << ", " << std::get<2>(id) << ")" << std::endl;
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

      result.push_back(std::move(hits));
    }

    clock::time_point end = clock::now();
    double diff = to_s(end - begin);
    std::cout << n_ix << " intersects, " << diff/1000000. << "s" << std::endl;
    double t_per_ix = diff / n_ix;
    std::cout << "=> " << t_per_ix << "us per intersect" << std::endl;
    double ix_per_sec = n_ix / (diff/1000000.);
    std::cout << "=> " << ix_per_sec << " intersects per second" << std::endl;


    return {n_ix, diff, result};
  }

}  // namespace Test
}  // namespace Acts

int
main()
{
  using oid_t = std::tuple<size_t, size_t, size_t>;
  //using res_t = std::tuple<size_t, double, id_t>;


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

  constexpr size_t n_sides = 9;
  
  std::ofstream os("octree_fr.csv");
  os << "octd,nix,time\n";
  std::vector<int> octree_depths = {0, 1, 2, 3, 4, 5};
  size_t n_ix;
  double diff;
  std::vector<std::set<oid_t>> hits;
  
  // get ref
  std::tie(n_ix, diff, hits) = Acts::Test::bench_frustum<n_sides>(-1);

  std::vector<std::set<oid_t>> ref = hits;
  //std::cout << "n hits ref: " << ref.size() << std::endl;

  auto print = [](const auto& tup) {
    std::cout << "(" << std::get<0>(tup) << ", " << std::get<1>(tup) << ", " << std::get<2>(tup) << ")";
  };

  for(int octree_depth : octree_depths) {
    std::tie(n_ix, diff, hits) = Acts::Test::bench_frustum<n_sides>(octree_depth);

    std::cout << "n hits: " << hits.size() << std::endl;
    //ASSERT(ref.size() == hits.size());

    for(size_t i=0;i<ref.size();i++) {
      std::set<oid_t>& ref_set = ref[i];
      std::set<oid_t>& act_set = hits[i];
      if(ref_set != act_set) {
        std::cout << "mismatch at idx " << i << ":" << std::endl;
        //print(ref[i]);
        //std::cout << " ";
        //print(hits[i]);
        //std::cout << std::endl;
        
        std::set<oid_t> diff;
        std::set_difference(ref_set.begin(), ref_set.end(), act_set.begin(), act_set.end(),
            std::inserter(diff, diff.end()));

        std::cout << "in ref but not in act" << std::endl;
        for(oid_t id : diff) {
          print(id);
          std::cout << std::endl;
        }

        abort();
      }

    }

    os << octree_depth << "," << n_ix << "," << diff << "\n";
  }
  os.close();

  return 0;
}
