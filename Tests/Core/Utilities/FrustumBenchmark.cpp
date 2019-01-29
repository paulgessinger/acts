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
  void
  bench_frustum()
  {
    using Frustum = Frustum<real_t, 3, sides>;

    size_t n    = 10;
    real_t  min  = -33;
    real_t  max  = 33;
    real_t  step = (max - min) / real_t(n);
    
    ply_helper<real_t> ply;

    Object o;
    Box::Size size(ActsVectorF<3>(2, 2, 2));

    std::vector<Box> boxes;
    boxes.reserve((n+1)*(n+1)*(n+1));

    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          ActsVectorF<3> pos(min + i * step, min + j * step, min + k * step);
          boxes.emplace_back(o, pos, size);
          boxes.back().draw(ply);
        }
      }
    }

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

    std::ofstream os("test.ply");
    os << ply << std::flush;
    os.close();


    clock::time_point begin = clock::now();

    size_t n_ix = 0;

    for(const Frustum& fr : frustums) {
      for(const Box& bb : boxes) {
        bb.intersect(fr);
        n_ix++;
      }
    }

    clock::time_point end = clock::now();
    double diff = to_s(end - begin);
    std::cout << n_ix << " intersects, " << diff/1000000. << "s" << std::endl;
    double t_per_ix = diff / n_ix;
    std::cout << "=> " << t_per_ix << "us per intersect" << std::endl;
    double ix_per_sec = n_ix / (diff/1000000.);
    std::cout << "=> " << ix_per_sec << " intersects per second" << std::endl;


  }

}  // namespace Test
}  // namespace Acts

int
main()
{

  // Tested with sides = 10, NaN / inf behavior!
  Acts::Test::bench_frustum<4>();

  return 0;
}
