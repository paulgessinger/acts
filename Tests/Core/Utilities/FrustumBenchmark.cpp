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
  id_t        id;
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
                const coll_t&                                      objects)
  {

    static_assert(std::is_same<T, typename coll_t::value_type>::value,
                  "Invalid type combo");
    std::cout << "Bench octd: " << octree_depth << std::endl;

    ply_helper<real_t> ply;

    std::vector<Box*>                 prims;
    std::vector<std::unique_ptr<Box>> boxes;
    boxes = boxFactory();

    for (auto& box : boxes) {
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

    // for (Box* box : prims) {
    // box->draw(ply);
    //}

    //// os << ply << std::flush;
    // os.close();

    ply.clear();
    for (auto& box : boxes) {
      box->draw(ply);
    }

    std::stringstream fn;
    fn << "octree_d" << octree_depth << ".ply";
    os = std::ofstream(fn.str());
    os << ply << std::flush;
    os.close();

    ply.clear();
    // print objects
    for (size_t o = 0; o < objects.size(); o++) {
      objects[o].draw(ply, 3000);
      if (o > 50) break;
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
            // << std::get<1>(id) << ", " << std::get<2>(id) << ")" <<
            // std::endl;
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

  std::vector<int> octree_depths = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  size_t           n_ix;
  double           diff;
  std::vector<std::set<oid_t>> hits;

  size_t n_tests = 1e4;

  std::mt19937 rng(42);

  std::uniform_real_distribution<real_t> dir_dist(-1, 1);
  // std::uniform_real_distribution<real_t> angle_dist(0.1 * M_PI,
  // M_PI / 4. * 0.9);

  auto frustFactory = [&](real_t min, real_t max, real_t angle) {
    std::uniform_real_distribution<real_t> pos_dist(min * 1.5, max * 1.5);
    vec_t pos(pos_dist(rng), pos_dist(rng), pos_dist(rng));
    vec_t dir(dir_dist(rng), dir_dist(rng), dir_dist(rng));
    dir.normalize();
    // real_t angle = angle_dist(rng);

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
  auto                                 gridBoxFactory = [&]() {
    Box::Size size(Acts::ActsVectorF<3>(2, 2, 2));

    std::vector<std::unique_ptr<Box>> boxes;
    boxes.reserve((n + 1) * (n + 1) * (n + 1));

    std::cout << "generating: " << (n + 1) * (n + 1) * (n + 1)
              << " bounding boxes" << std::endl;

    size_t idx = 0;
    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          Acts::ActsVectorF<3> pos(
              min + i * step, min + j * step, min + k * step);
          // boxes.emplace_back(o, pos, size);
          // Object o{{i, j, k}};
          // Object o;
          // o.id = {i, j, k};
          entities.push_back(std::make_unique<Object>(idx));
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

    float                             l = 100;
    Box::Size                         size(Acts::ActsVectorF<3>(l, l, l));
    std::vector<std::unique_ptr<Box>> boxes;

    std::ifstream is("../output_geo.csv");
    std::string   line("");
    size_t        idx = 0;
    //,x,y,z,r,phi_raw,eta_raw,dphi,deta,dr,dx,dy,dz,calosample

    size_t row;
    float  x, y, z, r, phi_raw, eta_raw, dphi, deta, dr, dx, dy, dz;
    size_t calosample;
    char   del;
    float  scale;

    // strip header row
    std::getline(is, line);

    auto eta_to_theta
        = [](double eta) { return 2 * std::atan(std::exp(-eta)); };

    auto build_endcap = [&eta_to_theta](auto&  ply,
                                        double r,
                                        double z,
                                        double dz,
                                        double eta,
                                        double deta,
                                        double phi,
                                        double dphi) {
      std::cout << "build endcap" << std::endl;

      double eta_max   = eta + deta * 0.5;
      double eta_min   = eta - deta * 0.5;
      double theta_max = eta_to_theta(eta_max);
      double theta_min = eta_to_theta(eta_min);
      double phi_max   = phi + dphi * 0.5;
      double phi_min   = phi - dphi * 0.5;
      double z_min     = z - dz;
      double z_max     = z + dz;

      double         r_min, r_max;
      Acts::Vector3D p1, p2, p3, p4, p5, p6, p7, p8;

      // inner face
      r_min = std::tan(theta_min) * z_min;
      r_max = std::tan(theta_max) * z_min;

      p1 << r_min * std::cos(phi_min), r_min * std::sin(phi_min), z_min;
      p2 << r_min * std::cos(phi_max), r_min * std::sin(phi_max), z_min;
      p3 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_min;
      p4 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_min;

      ply.vertex(p1);
      ply.vertex(p2);
      ply.vertex(p3);
      ply.vertex(p4);
      ply.face(std::vector<Acts::Vector3D>({p1, p2, p3, p4}));

      // outer face
      r_min = std::tan(theta_min) * z_max;
      r_max = std::tan(theta_max) * z_max;

      p5 << r_min * std::cos(phi_min), r_min * std::sin(phi_min), z_max;
      p6 << r_min * std::cos(phi_max), r_min * std::sin(phi_max), z_max;
      p7 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_max;
      p8 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_max;

      ply.vertex(p5);
      ply.vertex(p6);
      ply.vertex(p7);
      ply.vertex(p8);
      ply.face(std::vector<Acts::Vector3D>({p5, p6, p7, p8}));

      // top face
      ply.face(std::vector<Acts::Vector3D>({p3, p4, p8, p7}));

      // bottom face
      ply.face(std::vector<Acts::Vector3D>({p1, p2, p6, p5}));

      // left face
      ply.face(std::vector<Acts::Vector3D>({p1, p4, p8, p5}));

      // right face
      ply.face(std::vector<Acts::Vector3D>({p2, p3, p7, p6}));

    };

    auto build_barrel = [&eta_to_theta](auto&  ply,
                                        double r,
                                        double dr,
                                        double z,
                                        double eta,
                                        double deta,
                                        double phi,
                                        double dphi) {
      std::cout << "build barrel" << std::endl;
      double eta_max   = eta + deta * 0.5;
      double eta_min   = eta - deta * 0.5;
      double theta_max = eta_to_theta(eta_max);
      double theta_min = eta_to_theta(eta_min);
      double phi_max   = phi + dphi * 0.5;
      double phi_min   = phi - dphi * 0.5;

      double r_min = r - dr;
      double r_max = r + dr;

      double         z_min, z_max;
      Acts::Vector3D p1, p2, p3, p4, p5, p6, p7, p8;

      // inner face
      z_min = r_min / std::tan(theta_min);
      z_max = r_min / std::tan(theta_max);

      p1 << r_min * std::cos(phi_min), r_min * std::sin(phi_min), z_min;
      p2 << r_min * std::cos(phi_min), r_min * std::sin(phi_min), z_max;
      p3 << r_min * std::cos(phi_max), r_min * std::sin(phi_max), z_max;
      p4 << r_min * std::cos(phi_max), r_min * std::sin(phi_max), z_min;

      ply.vertex(p1);
      ply.vertex(p2);
      ply.vertex(p3);
      ply.vertex(p4);
      ply.face(std::vector<Acts::Vector3D>({p1, p2, p3, p4}));

      // outer face
      z_min = r_max / std::tan(theta_min);
      z_max = r_max / std::tan(theta_max);

      p5 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_min;
      p6 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_max;
      p7 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_max;
      p8 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_min;

      ply.face(std::vector<Acts::Vector3D>({p5, p6, p7, p8}));

      // top face
      ply.face(std::vector<Acts::Vector3D>({p3, p4, p8, p7}));

      // bottom face
      ply.face(std::vector<Acts::Vector3D>({p1, p2, p6, p5}));

      // left face
      ply.face(std::vector<Acts::Vector3D>({p1, p4, p8, p5}));

      // right face
      ply.face(std::vector<Acts::Vector3D>({p2, p3, p7, p6}));

    };

    auto build_box = [](auto&  ply,
                        double x,
                        double dx,
                        double y,
                        double dy,
                        double z,
                        double dz) {
      std::cout << "build box" << std::endl;

      double x_min, x_max, y_min, y_max, z_min, z_max;
      x_min = x - dx;
      x_max = x + dx;
      y_min = y - dy;
      y_max = y + dy;
      z_min = z - dz;
      z_max = z + dz;

      Acts::Vector3D p1, p2, p3, p4, p5, p6, p7, p8;

      // inner face
      p1 << x_min, y_min, z_min;
      p2 << x_min, y_max, z_min;
      p3 << x_max, y_max, z_min;
      p4 << x_max, y_min, z_min;

      ply.face(std::vector<Acts::Vector3D>({p1, p2, p3, p4}));

      // outer face
      p5 << x_min, y_min, z_max;
      p6 << x_min, y_max, z_max;
      p7 << x_max, y_max, z_max;
      p8 << x_max, y_min, z_max;

      ply.face(std::vector<Acts::Vector3D>({p5, p6, p7, p8}));

      // top face
      ply.face(std::vector<Acts::Vector3D>({p2, p3, p7, p6}));

      // bottom face
      ply.face(std::vector<Acts::Vector3D>({p1, p4, p8, p5}));

      // left face
      ply.face(std::vector<Acts::Vector3D>({p1, p2, p6, p5}));

      // right face
      ply.face(std::vector<Acts::Vector3D>({p3, p4, p8, p7}));

    };

    Acts::ply_helper<double> ply_lar;
    Acts::ply_helper<double> ply_tile;
    Acts::ply_helper<double> ply_fcal;

    while (std::getline(is, line)) {
      std::istringstream iss(line);
      iss >> row;
      iss >> del;
      iss >> x;
      iss >> del;
      iss >> y;
      iss >> del;
      iss >> z;
      iss >> del;
      iss >> r;
      iss >> del;
      iss >> phi_raw;
      iss >> del;
      iss >> eta_raw;
      iss >> del;
      iss >> dphi;
      iss >> del;
      iss >> deta;
      iss >> del;
      iss >> dr;
      iss >> del;
      iss >> dx;
      iss >> del;
      iss >> dy;
      iss >> del;
      iss >> dz;
      iss >> del;
      iss >> calosample;

      // ec = [4, 5, 6, 7, 8, 9, 10, 11, 17]
      // brl = [0, 1, 2, 3, 12, 13, 14, 15, 16, 18, 19, 20]

      scale = 1.;
      if (calosample >= 12 && calosample <= 20) {
        scale = 0.5;
      }

      Acts::ply_helper<double>* ply;
      if (calosample <= 11) {
        ply = &ply_lar;
      } else if (calosample <= 20) {
        ply = &ply_tile;
      } else {
        ply = &ply_fcal;
      }

      switch (calosample) {
      case 4:
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
      case 10:
      case 17:
        dz *= scale;
        build_endcap(*ply, r, z, dz, eta_raw, deta, phi_raw, dphi);
        break;
      case 0:
      case 1:
      case 2:
      case 3:
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 18:
      case 19:
      case 20:
        dr *= scale;
        build_barrel(*ply, r, dr, z, eta_raw, deta, phi_raw, dphi);
        break;
      case 21:
      case 22:
      case 23:
        scale = 0.5;
        dx *= scale;
        dy *= scale;
        // dz *= scale;
        build_box(*ply, x, dx, y, dy, z, dz);
        break;
      default:
        std::stringstream ss;
        ss << "Unkown calo sample " << calosample;
        std::runtime_error(ss.str());
      }

      // std::cout << x << " " << y << " " << z << " " << r << " " << phi_raw <<
      // " ";
      // std::cout << eta_raw << " " << dphi << " " << deta  << " " << dr << "
      // ";
      // std::cout << dx << " " << dy << " " << dz << " " << calosample <<
      // std::endl;

      // std::cout << line << std::endl;
      // std::cout << x << ";" << y << ";" << z << std::endl;

      // Acts::ActsVectorF<3> pos(x, y, z);
      // entities.push_back(
      // std::make_unique<Object>(idx));
      // boxes.push_back(std::make_unique<Box>(*entities.back(), pos, size));
      // Box& b = *boxes.back();
      // min = std::min(b.min().x(), min);
      // min = std::min(b.min().y(), min);
      // min = std::min(b.min().z(), min);
      // max = std::max(b.max().x(), max);
      // max = std::max(b.max().y(), max);
      // max = std::max(b.max().z(), max);

      idx++;
    }

  endloop:

    std::ofstream os("lar.ply");
    os << ply_lar << std::flush;
    os.close();

    os = std::ofstream("tile.ply");
    os << ply_tile << std::flush;
    os.close();

    os = std::ofstream("fcal.ply");
    os << ply_fcal << std::flush;
    os.close();

    return boxes;
  };

  atlasCaloFactory();

  /*
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
    os << -1 << "," << n_tests << "," << angle << "," << n_ix << "," << diff <<
    "\n";

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

      os << octree_depth << "," << n_tests << "," << angle << "," << n_ix << ","
    << diff << "\n";
    }
    os.close();
    */

  return 0;
}
