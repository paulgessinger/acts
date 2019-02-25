// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/AbstractVolume.hpp"
#include "Acts/Volumes/GenericCuboidVolumeBounds.hpp"
#include "Acts/Volumes/Volume.hpp"

namespace {
struct Object
{
  using id_t = size_t;
  Object(id_t                 id_) : id(id_) {}
  id_t                        id;
  const Acts::AbstractVolume* volume;
};

using oid_t = size_t;

using real_t = double;
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
    std::cout << "Build geometry" << std::endl;
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
            // auto id = lnode->entity()->id;
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

double
eta_to_theta(double eta)
{
  return 2 * std::atan(std::exp(-eta));
}

Acts::AbstractVolume
build_endcap(

    double z,
    double dz,
    double eta,
    double deta,
    double phi,
    double dphi)
{
  // std::cout << "build endcap" << std::endl;

  double eta_max   = eta + deta * 0.5;
  double eta_min   = eta - deta * 0.5;
  double theta_max = eta_to_theta(eta_max);
  double theta     = eta_to_theta(eta);
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

  // outer face
  r_min = std::tan(theta_min) * z_max;
  r_max = std::tan(theta_max) * z_max;

  p5 << r_min * std::cos(phi_min), r_min * std::sin(phi_min), z_max;
  p6 << r_min * std::cos(phi_max), r_min * std::sin(phi_max), z_max;
  p7 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_max;
  p8 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_max;

  double         r_mid = std::tan(theta) * z_min;
  Acts::Vector3D center;
  center.x() = r_mid * std::cos(phi);
  center.y() = r_mid * std::sin(phi);
  center.z() = z;

  Acts::Transform3D glob2vol = Acts::Transform3D::Identity();
  glob2vol *= Acts::AngleAxis3D(-phi, Acts::Vector3D::UnitZ());
  glob2vol *= Acts::AngleAxis3D(
      -theta, Acts::Vector3D::UnitZ().cross(center).normalized());
  glob2vol
      *= Acts::Translation3D(-(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / 8.);

  p1 = glob2vol * p1;
  p2 = glob2vol * p2;
  p3 = glob2vol * p3;
  p4 = glob2vol * p4;
  p5 = glob2vol * p5;
  p6 = glob2vol * p6;
  p7 = glob2vol * p7;
  p8 = glob2vol * p8;

  auto globalToLocal = std::make_shared<Acts::Transform3D>(glob2vol.inverse());

  auto cubo = std::make_shared<Acts::GenericCuboidVolumeBounds>(
      std::array<Acts::Vector3D, 8>({{p1, p2, p3, p4, p5, p6, p7, p8}}));
  Acts::AbstractVolume vol(std::move(globalToLocal), std::move(cubo));

  return vol;
}

Acts::AbstractVolume
build_barrel(double r,
             double dr,
             double eta,
             double deta,
             double phi,
             double dphi)
{
  // std::cout << "build barrel" << std::endl;
  double eta_max   = eta + deta * 0.5;
  double eta_min   = eta - deta * 0.5;
  double theta     = eta_to_theta(eta);
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

  // outer face
  z_min = r_max / std::tan(theta_min);
  z_max = r_max / std::tan(theta_max);

  p5 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_min;
  p6 << r_max * std::cos(phi_min), r_max * std::sin(phi_min), z_max;
  p7 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_max;
  p8 << r_max * std::cos(phi_max), r_max * std::sin(phi_max), z_min;

  Acts::Vector3D center;
  center.x() = r * std::cos(phi);
  center.y() = r * std::sin(phi);
  center.z() = r / std::tan(theta);

  Acts::Transform3D glob2vol = Acts::Transform3D::Identity();
  glob2vol *= Acts::AngleAxis3D(-phi, Acts::Vector3D::UnitZ());
  glob2vol *= Acts::AngleAxis3D(
      -theta, Acts::Vector3D::UnitZ().cross(center).normalized());
  glob2vol
      *= Acts::Translation3D(-(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / 8.);

  p1 = glob2vol * p1;
  p2 = glob2vol * p2;
  p3 = glob2vol * p3;
  p4 = glob2vol * p4;
  p5 = glob2vol * p5;
  p6 = glob2vol * p6;
  p7 = glob2vol * p7;
  p8 = glob2vol * p8;

  auto globalToLocal = std::make_shared<Acts::Transform3D>(glob2vol.inverse());

  auto cubo = std::make_shared<Acts::GenericCuboidVolumeBounds>(
      std::array<Acts::Vector3D, 8>({{p1, p2, p3, p4, p5, p6, p7, p8}}));

  Acts::AbstractVolume vol(std::move(globalToLocal), std::move(cubo));

  return vol;
}

Acts::AbstractVolume
build_box(double x, double dx, double y, double dy, double z, double dz)
{
  // std::cout << "build box" << std::endl;

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

  // outer face
  p5 << x_min, y_min, z_max;
  p6 << x_min, y_max, z_max;
  p7 << x_max, y_max, z_max;
  p8 << x_max, y_min, z_max;

  Acts::Transform3D glob2vol = Acts::Transform3D::Identity();
  glob2vol
      *= Acts::Translation3D(-(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / 8.);

  p1 = glob2vol * p1;
  p2 = glob2vol * p2;
  p3 = glob2vol * p3;
  p4 = glob2vol * p4;
  p5 = glob2vol * p5;
  p6 = glob2vol * p6;
  p7 = glob2vol * p7;
  p8 = glob2vol * p8;

  auto globalToLocal = std::make_shared<Acts::Transform3D>(glob2vol.inverse());

  auto cubo = std::make_shared<Acts::GenericCuboidVolumeBounds>(
      std::array<Acts::Vector3D, 8>({{p1, p2, p3, p4, p5, p6, p7, p8}}));
  Acts::AbstractVolume vol(std::move(globalToLocal), std::move(cubo));

  return vol;
};

template <typename object_t>
void
do_octree_scan(size_t                                             n_tests,
               std::function<std::vector<std::unique_ptr<Box>>()> boxFactory,
               std::function<object_t()>                          objectFactory)
{
  std::vector<int> octree_depths = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  size_t           n_ix;
  double           diff;
  std::vector<std::set<oid_t>> hits;
  std::cout << "Generating " << n_tests << " test objects" << std::endl;

  std::vector<object_t, Eigen::aligned_allocator<object_t>> objects;
  objects.reserve(n_tests);

  for (size_t i = 0; i < n_tests; i++) {

    // objects.emplace_back(pos, dir, angle);
    objects.push_back(objectFactory());
    if (i < 50) {
      // frustums.back().draw(ply);
    }
  }

  // get ref
  std::tie(n_ix, diff, hits)
      = Acts::Test::bench_frustum<object_t>(-1, boxFactory, objects);
  std::vector<std::set<oid_t>> ref = hits;

  std::ofstream os("octree_fr.csv");
  os << "octd,ntests,nix,time\n";
  os << -1 << "," << n_tests << "," << n_ix << "," << diff << "\n";

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

        std::set<oid_t> diff_;
        std::set_difference(ref_set.begin(),
                            ref_set.end(),
                            act_set.begin(),
                            act_set.end(),
                            std::inserter(diff_, diff_.end()));

        std::cout << "in ref but not in act" << std::endl;
        for (oid_t id : diff_) {
          std::cout << id;
          std::cout << std::endl;
        }

        abort();
      }
    }

    os << octree_depth << "," << n_tests << "," << n_ix << "," << diff << "\n";
  }
  os.close();
}

int
main()
{
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

  std::mt19937 rng(42);

  // std::uniform_real_distribution<real_t> angle_dist(0.1 * M_PI,
  // M_PI / 4. * 0.9);

  size_t n    = 59;
  real_t min  = -10;
  real_t max  = 10;
  real_t step = (max - min) / real_t(n);

  real_t angle      = M_PI / 2.;
  auto frustFactory = [&](auto& pos_dist, auto& dir_dist) {
    vec_t pos(pos_dist(rng), pos_dist(rng), pos_dist(rng));
    vec_t dir(dir_dist(rng), dir_dist(rng), dir_dist(rng));
    dir.normalize();
    // real_t angle = angle_dist(rng);

    // std::cout << "pos: " << pos.transpose() << std::endl;
    // std::cout << "dir: " << dir.transpose() << std::endl;
    // std::cout << "angle: " << angle << std::endl;

    return Frustum(pos, dir, angle);
  };
  (void)frustFactory;

  auto rayFactory = [&](auto& pos_dist, auto& dir_dist) {
    vec_t pos(pos_dist(rng), pos_dist(rng), pos_dist(rng));
    vec_t dir(dir_dist(rng), dir_dist(rng), dir_dist(rng));
    dir.normalize();
    return Ray(pos, dir);
  };
  (void)rayFactory;

  std::vector<std::unique_ptr<Object>> entities;
  auto                                 gridBoxFactory = [&]() {
    Box::Size size(Acts::ActsVectorD<3>(2, 2, 2));

    std::vector<std::unique_ptr<Box>> boxes;
    boxes.reserve((n + 1) * (n + 1) * (n + 1));

    std::cout << "generating: " << (n + 1) * (n + 1) * (n + 1)
              << " bounding boxes" << std::endl;

    size_t idx = 0;
    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          Acts::ActsVectorD<3> pos(
              min + i * step, min + j * step, min + k * step);
          // boxes.emplace_back(o, pos, size);
          // Object o{{i, j, k}};
          // Object o;
          // o.id = {i, j, k};
          entities.push_back(std::make_unique<Object>(idx));
          idx++;
          boxes.push_back(
              std::make_unique<Box>(entities.back().get(), pos, size));
          // boxes.back()->draw(ply);
        }
      }
    }

    return boxes;
  };
  (void)gridBoxFactory;

  std::vector<std::unique_ptr<Acts::AbstractVolume>> cells;
  auto                                               atlasCaloFactory = [&]() {
    min = 1e10;
    max = -1e10;

    float                             l = 100;
    Box::Size                         size(Acts::ActsVectorD<3>(l, l, l));
    std::vector<std::unique_ptr<Box>> boxes;

    std::ifstream is("../output_geo.csv");
    std::string   line("");
    size_t        idx = 0;

    Acts::ply_helper<double> ply_lar;
    Acts::ply_helper<double> ply_tile;
    Acts::ply_helper<double> ply_fcal;

    //,x,y,z,r,phi_raw,eta_raw,dphi,deta,dr,dx,dy,dz,calosample
    size_t row;
    float  x, y, z, r, phi_raw, eta_raw, dphi, deta, dr, dx, dy, dz;
    size_t calosample;
    char   del;
    float  scale;

    // strip header row
    std::getline(is, line);

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
        cells.push_back(std::make_unique<Acts::AbstractVolume>(
            build_endcap(z, dz, eta_raw, deta, phi_raw, dphi)));
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
        cells.push_back(std::make_unique<Acts::AbstractVolume>(
            build_barrel(r, dr, eta_raw, deta, phi_raw, dphi)));
        break;
      case 21:
      case 22:
      case 23:
        scale = 0.5;
        dx *= scale;
        dy *= scale;
        // dz *= scale;
        cells.push_back(std::make_unique<Acts::AbstractVolume>(
            build_box(x, dx, y, dy, z, dz)));
        break;
      default:
        std::stringstream ss;
        ss << "Unkown calo sample " << calosample;
        std::runtime_error(ss.str());
      }

      // cells.back().boundingBox({0.1, 0.1, 0.1}).draw(*ply);
      // need to copy box again, because we use our custom object here
      // @TODO: change this

      auto vbox = cells.back()->boundingBox({0.1, 0.1, 0.1});

      entities.push_back(std::make_unique<Object>(idx));
      entities.back()->volume = cells.back().get();
      boxes.push_back(
          std::make_unique<Box>(entities.back().get(), vbox.min(), vbox.max()));
      // boxes.push_back(std::make_unique<Box>(cells.back().boundingBox({0.1,
      // 0.1, 0.1})));

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

  // auto boxFactory = gridBoxFactory
  auto boxFactory = atlasCaloFactory;
  // auto objectFactory = rayFactory;
  // using object_t     = Ray;
  // auto objectFactory = frustFactory;
  // using object_t     = Frustum;

  // size_t n_tests = 1e4;
  // do_octree_scan<object_t>(n_tests, boxFactory, objectFactory);

  std::uniform_real_distribution<real_t> dir_dist(-1, 1);
  std::uniform_real_distribution<real_t> pos_dist(0., 0.);

  auto              boxes = boxFactory();
  std::vector<Box*> prims;
  for (auto& box : boxes) {
    prims.push_back(box.get());
  }
  Box* top_node = nullptr;
  top_node      = make_octree(boxes, prims, 6);

  Acts::ply_helper<float> ply;

  auto intersections = [&](const auto& obj, const Box* top) {
    const Box*              lnode = top;
    std::vector<const Box*> hits;
    do {
      if (lnode->intersect(obj)) {

        if (lnode->hasEntity()) {
          // found primitive
          hits.push_back(lnode);
          lnode = lnode->getSkip();
        } else {
          // go over children
          lnode = lnode->getLeftChild();
        }
      } else {
        lnode = lnode->getSkip();
      }
    } while (lnode != nullptr);
    return hits;
  };

  // for (const auto& box : hits) { box->draw(ply); }

  // std::ofstream os("boxes.ply");
  // os << ply;
  // os.close();
  // ply.clear();

  size_t nrays = 1e6;

  std::ofstream os("ray_efficiency.csv");
  os << "n,eta,phi,boxes,obb,cells\n";

  std::uniform_real_distribution<float> eta_dist(-5, 5);
  std::uniform_real_distribution<float> phi_dist(-M_PI, M_PI);
  for (size_t i = 0; i < nrays; i++) {
    float eta = eta_dist(rng);
    float phi = phi_dist(rng);

    float          theta = 2 * std::atan(std::exp(-eta));
    Acts::Vector3D dir;
    dir << std::cos(phi), std::sin(phi), 1. / std::tan(theta);
    dir.normalize();

    // std::cout << "eta phi: " << eta << ", " << phi << std::endl;
    // std::cout << dir.transpose() << std::endl;
    // std::cout << Acts::VectorHelpers::phi(dir) << std::endl;
    // std::cout << Acts::VectorHelpers::eta(dir) << std::endl;

    Ray ray({0, 0, 0}, dir);

    auto hits = intersections(ray, top_node);

    size_t boxes_hit = 0;
    size_t cells_hit = 0;
    size_t obb_hit   = 0;
    for (const auto& box : hits) {
      const Acts::AbstractVolume* vol;
      vol = box->entity()->volume;
      // auto vbo = dynamic_cast<const Acts::GenericCuboidVolumeBounds*>(
      //&vol->volumeBounds());
      // std::cout << *vol << std::endl;

      boxes_hit++;

      auto obb = vol->orientedBoundingBox();
      // do we hit the obb?
      if (obb.intersect(ray.transformed(vol->transform().inverse()))) {
        obb_hit++;
      }

      // check if we actually hit it
      auto surfaces = vol->boundarySurfaces();
      for (const auto& boundarySurface : surfaces) {
        auto& srf = boundarySurface->surfaceRepresentation();
        auto ix = srf.intersectionEstimate(ray.origin().template cast<double>(),
                                           ray.dir().template cast<double>(),
                                           Acts::forward,
                                           true);

        if (ix) {
          // std::cout << "valid surface intersect" << std::endl;
          cells_hit++;
          // vbo->draw(ply, &vol->transform());
          break;  // no need to test other boundary surfaces
        }
      }
    }

    os << i << "," << eta << "," << phi << "," << boxes_hit << "," << obb_hit
       << "," << cells_hit << "\n";

    // std::cout << "boxes hit: " << boxes_hit << " cells_hit: " << cells_hit
    //<< std::endl;
    // os = std::ofstream("cells.ply");
    // os << ply;
    // os.close();

    // ply.clear();
    // ray.draw(ply, 10000);
    // os = std::ofstream("rays.ply");
    // os << ply;
    // os.close();
  }
  os.close();

  return 0;
}
