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

namespace {
struct Object
{
  using id_t = size_t;
  Object(id_t id_) : id(id_) {}
  id_t id;
};

using real_t = float;
using vec_t  = Acts::ActsVector<real_t, 3>;
using Box    = Acts::AxisAlignedBoundingBox<Object, real_t, 3>;
using clock  = std::chrono::steady_clock;
}  // namespace

#define ASSERT(condition)                                                      \
  if (!(condition)) {                                                          \
    std::cerr << "Assertion failed in line #" << __LINE__ << std::endl;        \
    abort();                                                                   \
  }

namespace Acts {
namespace Test {

  template <typename T>
  double
  to_s(const T& delta)
  {
    return std::chrono::duration_cast<std::chrono::microseconds>(delta).count();
  }

  template <typename T, typename coll_t>
  std::tuple<size_t, double, std::vector<std::set<Object::id_t>>>
  bench_frustum(int                                                octree_depth,
                std::function<std::vector<std::unique_ptr<Box>>()> boxFactory,
                const coll_t& objects)
  {

    static_assert(std::is_same<T, typename coll_t::value_type>::value, "Invalid type combo");
    std::cout << "Bench octd: " << octree_depth << std::endl;

    ply_helper<real_t> ply;

    std::vector<Box*> prims;
    std::vector<std::unique_ptr<Box>> boxes;
    boxes = boxFactory();

    for(auto& box : boxes) {
      prims.push_back(box.get());
    }

    Box* top_node = nullptr;

    if (octree_depth >= 0) {
      top_node = make_octree(boxes, prims, octree_depth);
    } else {
      // no octree, just link all boxes sequentially
      for (size_t i = 1; i < prims.size(); i++) {
        prims.at(i - 1)->setSkip(prims.at(i));
      }
      top_node = prims.front();
    }

    std::ofstream os("prims.ply");

    //for (Box* box : prims) {
       //box->draw(ply);
    //}

    //// os << ply << std::flush;
    //os.close();

    ply.clear();
    for (auto& box : boxes) { box->draw(ply); }

    std::stringstream fn;
    fn << "octree_d" << octree_depth << ".ply";
    os = std::ofstream(fn.str());
    os << ply << std::flush;
    os.close();

    ply.clear();
    // print objects
    for(size_t o=0;o<objects.size();o++) {
      objects[o].draw(ply, 3000);
      if(o>50) break;
    }
    fn = std::stringstream();
    fn << "octree_d" << octree_depth << "_objects.ply";
    os = std::ofstream(fn.str());
    os << ply << std::flush;
    os.close();



    std::cout << "begin intersect timing" << std::endl;
    clock::time_point begin = clock::now();

    // size_t n_ix = frustums.size() * boxes.size();
    size_t n_ix = 0;

    // for(const Frustum& fr : frustums) {
    // for(const Box* bb : prims) {
    // volatile bool ix = bb->intersect(fr);
    //}
    //}



    std::vector<std::set<Object::id_t>> result;

    for (const T& obj : objects) {
      const Box*             lnode = top_node;
      std::set<Object::id_t> hits;
      do {
        n_ix++;
        if (lnode->intersect(obj)) {

          if (lnode->hasEntity()) {
            // found primitive
            // hits.insert(lnode->entity()->id);
            auto id = lnode->entity()->id;
            // std::cout << "intersect prim: idx: (" << std::get<0>(id) << ", "
            // << std::get<1>(id) << ", " << std::get<2>(id) << ")" << std::endl;
            // std::cout << "intersect prim: " << lnode << " -> " <<
            // lnode->getSkip() << std::endl;
            lnode = lnode->getSkip();
          } else {
            // go over children
            // std::cout << "intersect bvh: " << lnode << " -> " <<
            // lnode->getLeftChild() << std::endl;
            lnode = lnode->getLeftChild();
          }
        } else {
          // std::cout << "no intersect: " << lnode << " -> " <<
          // lnode->getSkip() << std::endl;
          lnode = lnode->getSkip();
        }

        // if(n_intersects > 80) throw std::runtime_error("");
      } while (lnode != nullptr);

      result.push_back(std::move(hits));
    }

    clock::time_point end  = clock::now();
    double            diff = to_s(end - begin);
    std::cout << n_ix << " intersects, " << diff / 1000000. << "s" << std::endl;
    double t_per_ix = diff / n_ix;
    std::cout << "=> " << t_per_ix << "us per intersect" << std::endl;
    double ix_per_sec = n_ix / (diff / 1000000.);
    std::cout << "=> " << ix_per_sec << " intersects per second" << std::endl;

    return {n_ix, diff, result};
  }

}  // namespace Test
}  // namespace Acts

int
main()
{
  using oid_t = size_t;
  // using res_t = std::tuple<size_t, double, id_t>;

  // std::vector<std::pair<size_t, res_t>> results;

  // Tested with sides = 10, NaN / inf behavior!

  // results = {
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

  // std::ofstream os("sides.csv");
  // os << "sides,n_ix,diff\n";
  // for(const auto& pr : results) {
  // size_t sides = pr.first;
  // size_t n_ix;
  // double diff;
  // std::tie(n_ix, diff) = pr.second;
  // os << sides << "," << n_ix << "," << diff << "\n";
  // sides++;
  //}
  // os.close();

  constexpr size_t n_sides = 9;
  using Frustum            = Acts::Frustum<real_t, 3, n_sides>;
  using Ray                = Acts::Ray<real_t, 3>;

  std::vector<int>             octree_depths = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  size_t                       n_ix;
  double                       diff;
  std::vector<std::set<oid_t>> hits;

  size_t n_tests = 1e4;

  std::mt19937 rng(42);

  std::uniform_real_distribution<real_t> dir_dist(-1, 1);
  //std::uniform_real_distribution<real_t> angle_dist(0.1 * M_PI,
                                                    //M_PI / 4. * 0.9);

  auto frustFactory = [&](real_t min, real_t max, real_t angle) {
    std::uniform_real_distribution<real_t> pos_dist(min * 1.5, max * 1.5);
    vec_t pos(pos_dist(rng), pos_dist(rng), pos_dist(rng));
    vec_t dir(dir_dist(rng), dir_dist(rng), dir_dist(rng));
    dir.normalize();
    //real_t angle = angle_dist(rng);

    // std::cout << "pos: " << pos.transpose() << std::endl;
    // std::cout << "dir: " << dir.transpose() << std::endl;
    // std::cout << "angle: " << angle << std::endl;

    return Frustum(pos, dir, angle);
  };

  auto rayFactory = [&](real_t min, real_t max) {
    std::uniform_real_distribution<real_t> pos_dist(min * 1.5, max * 1.5);
    vec_t pos(pos_dist(rng), pos_dist(rng), pos_dist(rng));
    vec_t dir(dir_dist(rng), dir_dist(rng), dir_dist(rng));
    dir.normalize();
    return Ray(pos, dir);
  };

  size_t n    = 59;
  real_t min  = -200;
  real_t max  = 200;
  real_t step = (max - min) / real_t(n);

  std::vector<std::unique_ptr<Object>> entities;
  auto gridBoxFactory = [&]() {
    Box::Size size(Acts::ActsVectorF<3>(2, 2, 2));

    std::vector<std::unique_ptr<Box>> boxes;
    boxes.reserve((n + 1) * (n + 1) * (n + 1));

    std::cout << "generating: " << (n + 1) * (n + 1) * (n + 1)
              << " bounding boxes" << std::endl;

    size_t idx = 0;
    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          Acts::ActsVectorF<3> pos(min + i * step, min + j * step, min + k * step);
          // boxes.emplace_back(o, pos, size);
          // Object o{{i, j, k}};
          // Object o;
          // o.id = {i, j, k};
          entities.push_back(
              std::make_unique<Object>(idx));
          idx++;
          boxes.push_back(std::make_unique<Box>(*entities.back(), pos, size));
          // boxes.back()->draw(ply);
        }
      }
    }

    return boxes;
  };

  auto atlasCaloFactory = [&]() {
    min = 1e10;
    max = -1e10;

    float l = 100;
    Box::Size size(Acts::ActsVectorF<3>(l, l, l));
    std::vector<std::unique_ptr<Box>> boxes;

    std::ifstream is("../output_geo.csv");
    std::string line("");
    size_t idx = 0;
    while(std::getline(is, line)) {
      std::istringstream iss(line);
      char del;
      float x, y, z;
      iss >> x;
      iss >> del;
      iss >> y;
      iss >> del;
      iss >> z;

      //std::cout << line << std::endl;
      //std::cout << x << ";" << y << ";" << z << std::endl;
          
      Acts::ActsVectorF<3> pos(x, y, z);
      entities.push_back(
          std::make_unique<Object>(idx));
      boxes.push_back(std::make_unique<Box>(*entities.back(), pos, size));
      Box& b = *boxes.back();
      min = std::min(b.min().x(), min);
      min = std::min(b.min().y(), min);
      min = std::min(b.min().z(), min);
      max = std::max(b.max().x(), max);
      max = std::max(b.max().y(), max);
      max = std::max(b.max().z(), max);

      idx++;
    }

    return boxes;
  };


  //auto boxFactory = gridBoxFactory
  auto boxFactory = atlasCaloFactory;

   //auto objectFactory = rayFactory;
   //using object_t = Ray;
  auto objectFactory = frustFactory;
  using object_t     = Frustum;
  
  std::cout << "Generating " << n_tests << " test objects" << std::endl;

  std::vector<object_t, Eigen::aligned_allocator<object_t>> objects;
  objects.reserve(n_tests);

  real_t angle = M_PI/2.;
  
  for (size_t i = 0; i < n_tests; i++) {

    // objects.emplace_back(pos, dir, angle);
    objects.push_back(objectFactory(min, max, angle));
    if (i < 50) {
      // frustums.back().draw(ply);
    }
  }

  // get ref
  std::tie(n_ix, diff, hits)
      = Acts::Test::bench_frustum<object_t>(-1, boxFactory, objects);
  std::vector<std::set<oid_t>> ref = hits;
  
  std::ofstream os("octree_fr.csv");
  os << "octd,ntests,angle,nix,time\n";
  os << -1 << "," << n_tests << "," << angle << "," << n_ix << "," << diff << "\n";

  auto print = [](const auto& tup) {
    std::cout << "(" << std::get<0>(tup) << ", " << std::get<1>(tup) << ", "
              << std::get<2>(tup) << ")";
  };

  for (int octree_depth : octree_depths) {
    std::tie(n_ix, diff, hits) = Acts::Test::bench_frustum<object_t>(
        octree_depth, boxFactory, objects);

    for (size_t i = 0; i < ref.size(); i++) {
      std::set<oid_t>& ref_set = ref[i];
      std::set<oid_t>& act_set = hits[i];
      if (ref_set != act_set) {
        std::cout << "mismatch at idx " << i << ":" << std::endl;

        std::set<oid_t> diff;
        std::set_difference(ref_set.begin(),
                            ref_set.end(),
                            act_set.begin(),
                            act_set.end(),
                            std::inserter(diff, diff.end()));

        std::cout << "in ref but not in act" << std::endl;
        for (oid_t id : diff) {
          std::cout << id;
          std::cout << std::endl;
        }

        abort();
      }
    }

    os << octree_depth << "," << n_tests << "," << angle << "," << n_ix << "," << diff << "\n";
  }
  os.close();

  return 0;
}
