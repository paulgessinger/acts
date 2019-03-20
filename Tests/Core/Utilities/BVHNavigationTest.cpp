// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Utilities/Visualization.hpp"
#include "Acts/Volumes/AbstractVolume.hpp"
#include "Acts/Volumes/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"
#include "Acts/Volumes/GenericCuboidVolumeBounds.hpp"
#include "Acts/Volumes/Volume.hpp"

#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Tools/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"

#include <thread>
#include "tbb/tbb.h"

using namespace Acts::VectorHelpers;

namespace {
using clock = std::chrono::steady_clock;
}

double
eta_to_theta(double eta)
{
  return 2 * std::atan(std::exp(-eta));
}

Acts::AbstractVolume
build_endcap(double z,
             double dz,
             double eta,
             double deta,
             double phi,
             double dphi)
{
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
}

std::vector<std::unique_ptr<Acts::AbstractVolume>>
atlasCaloFactory(std::string input_file)
{
  std::ifstream is(input_file);
  std::string   line("");
  size_t        idx = 0;

  // Acts::ply_helper<double> ply_lar;
  // Acts::ply_helper<double> ply_tile;
  // Acts::ply_helper<double> ply_fcal;

  size_t row;
  float  x, y, z, r, phi_raw, eta_raw, dphi, deta, dr, dx, dy, dz;
  size_t calosample;
  char   del;
  float  scale;

  // storage of cells we will produce
  std::vector<std::unique_ptr<Acts::AbstractVolume>> cells;
  cells.reserve(180000);  // about 180k

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

    scale = 1.;
    if (calosample >= 12 && calosample <= 20) { scale = 0.5; }

    // Acts::ply_helper<double>* ply;
    // if (calosample <= 11) {
    // ply = &ply_lar;
    //} else if (calosample <= 20) {
    // ply = &ply_tile;
    //} else {
    // ply = &ply_fcal;
    //}

    switch (calosample) {
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
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
      scale = 1.;
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
    // auto cvb = dynamic_cast<const
    // Acts::GenericCuboidVolumeBounds*>(&cells.back()->volumeBounds());
    // cvb->draw(*ply, cells.back()->transform());
    idx++;
  }

  // std::ofstream os("lar.ply");
  // os << ply_lar << std::flush;
  // os.close();

  // os = std::ofstream("tile.ply");
  // os << ply_tile << std::flush;
  // os.close();

  // os = std::ofstream("fcal.ply");
  // os << ply_fcal << std::flush;
  // os.close();

  return cells;
}

int
main(int argc, char* argv[])
{

  size_t n_rays = 1e6;

  if (argc < 2) {
    std::cerr << "Provide filename for geo building" << std::endl;
    return 1;
  }

  std::string filename = argv[1];
  if (argc > 2) { n_rays = std::stoi(argv[2]); }

  std::cout << "Reading from: " << filename << std::endl;
  std::cout << "Build calo geometry..." << std::flush;
  std::vector<std::unique_ptr<Acts::AbstractVolume>> cells;
  cells = atlasCaloFactory(filename);
  std::cout << " => done!" << std::endl;

  // draw calo geometry
  // Acts::ply_helper<double> ply;
  // for (const auto& cell : cells) {
  // auto cvb = dynamic_cast<const
  // Acts::GenericCuboidVolumeBounds*>(&cell->volumeBounds()); cvb->draw(ply,
  // cell->transform());
  //}
  // std::ofstream os("calo.ply");
  // os << ply;
  // os.close();
  using Box = Acts::Volume::BoundingBox;
  using Ray = Acts::Ray<double, 3>;

  {
    std::mt19937                           rng(42);
    std::uniform_real_distribution<double> eta_dist(-5, 5);
    std::uniform_real_distribution<double> phi_dist(-M_PI, M_PI);
    std::uniform_real_distribution<double> z_dist(-100, 100);
    std::vector<Acts::Vector3D>            dirs;

    std::ofstream csv("../ray_purity.csv");
    csv << "n,eta,phi,aabb,obb,cells\n";

    Acts::ply_helper<double> ply_cells;
    Acts::ply_helper<double> ply_cells_all;
    Acts::ply_helper<double> ply_aabb;
    Acts::ply_helper<double> ply_obb;
    Acts::ply_helper<double> ply_ray;

    tbb::concurrent_vector<
        std::tuple<size_t, double, double, size_t, size_t, size_t>>
        rows;
    rows.reserve(n_rays);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, n_rays), [&](const auto& r) {
          thread_local std::mt19937 engine(
              std::hash<std::thread::id>{}(std::this_thread::get_id()));
          for (size_t i = r.begin(); i != r.end(); ++i) {
            size_t         n_cells = 0, n_aabb = 0, n_obb = 0;
            double         eta   = eta_dist(engine);
            double         phi   = phi_dist(engine);
            double         theta = 2 * std::atan(std::exp(-eta));
            Acts::Vector3D dir;
            dir << std::cos(phi), std::sin(phi), 1. / std::tan(theta);
            dir.normalize();

            Ray ray({0, 0, 0}, dir);
            // ray.draw(ply_ray, 10000);

            for (const auto& cell : cells) {

              // do we hit the bb?
              auto bb = cell->boundingBox();
              if (bb.intersect(ray)) {
                n_aabb++;
                auto vb = dynamic_cast<const Acts::GenericCuboidVolumeBounds*>(
                    &cell->volumeBounds());
                auto surfaces = vb->decomposeToSurfaces(&cell->transform());
                bool is_hit   = false;
                for (const auto& srf : surfaces) {
                  Acts::NavigationOptions<Acts::Surface> no(Acts::forward,
                                                            true);
                  if (srf->surfaceIntersectionEstimate(
                          ray.origin(), ray.dir(), no)) {
                    is_hit = true;
                    // vb->draw(ply_cells, cell->transform());
                    n_cells++;
                    break;
                  }
                }

                // std::cout << "hit" << std::endl;
                // if (!is_hit) { vb->draw(ply_cells_all, cell->transform()); }
                // bb.draw(ply_aabb);

                // draw obb
                auto obb = cell->orientedBoundingBox();
                if (obb.intersect(ray.transformed(cell->itransform()))) {
                  // obb.draw(ply_obb, {120, 120, 120}, cell->transform());
                  n_obb++;
                }
              }
            }

            // csv << i << "," << eta << "," << phi << "," << n_aabb << ","
            //<< n_obb << "," << n_cells << "\n";
            rows.push_back({i, eta, phi, n_aabb, n_obb, n_cells});
          }
        });

    for (const auto& row : rows) {
      auto [i, eta, phi, n_aabb, n_obb, n_cells] = row;
      csv << i << "," << eta << "," << phi << "," << n_aabb << "," << n_obb
          << "," << n_cells << "\n";
    }

    csv.close();

    std::ofstream eff("../eff_cells.ply");
    eff << ply_cells;
    eff = std::ofstream("../eff_cells_all.ply");
    eff << ply_cells_all;
    eff = std::ofstream("../eff_aabb.ply");
    eff << ply_aabb;
    eff = std::ofstream("../eff_obb.ply");
    eff << ply_obb;

    eff = std::ofstream("../eff_ray.ply");
    eff << ply_ray;
  }

  // auto intersections = [](const auto& obj, const Box* top) {
  // const Box*              lnode = top;
  // std::vector<const Box*> hits;
  // do {
  // if (lnode->intersect(obj)) {

  // if (lnode->hasEntity()) {
  //// found primitive
  //// check obb to limit false positivies
  // auto obb = lnode->entity()->orientedBoundingBox();
  // if (obb.intersect(obj.transformed(
  // lnode->entity()->transform().inverse().cast<float>()))) {
  // hits.push_back(lnode);
  //}
  //// we skip in any case, whether we actually hit the OBB or not
  // lnode = lnode->getSkip();
  //} else {
  //// go over children
  // lnode = lnode->getLeftChild();
  //}
  //} else {
  // lnode = lnode->getSkip();
  //}
  //} while (lnode != nullptr);
  // return hits;
  //};

  return 0;

  // create BVH for the calo geo
  std::vector<std::unique_ptr<Box>> boxStore;
  std::vector<Box*>                 prims;
  for (const auto& cell : cells) {
    boxStore.push_back(
        std::make_unique<Box>(cell->boundingBox({0.1, 0.1, 0.1})));
    prims.push_back(boxStore.back().get());
  }

  // CYLINDER VOLUME BOUNDS

  // double zmin = std::numeric_limits<double>::max();
  // double zmax = std::numeric_limits<double>::lowest();
  // double rmin = std::numeric_limits<double>::max();
  // double rmax = -1;
  // for (const auto& box : boxStore) {
  // Acts::Vector3D vmin = box->min().cast<double>();
  // Acts::Vector3D vmax = box->max().cast<double>();

  // zmin = std::min(zmin, vmin.z());
  // zmin = std::min(zmin, vmax.z());
  // zmax = std::max(zmax, vmin.z());
  // zmax = std::max(zmax, vmax.z());

  // rmin = std::min(rmin, Acts::VectorHelpers::perp(vmin));
  // rmin = std::min(rmin, Acts::VectorHelpers::perp(vmax));
  // rmax = std::max(rmax, Acts::VectorHelpers::perp(vmin));
  // rmax = std::max(rmax, Acts::VectorHelpers::perp(vmax));
  //}
  //// the cylinder volume bounds for the TV need to wrap around all the
  //// bounding box assuming this is symmetric right now

  // double halez = (zmax - zmin) / 2.;
  // rmin         = 0;
  // auto volBds  = std::make_shared<Acts::CylinderVolumeBounds>(
  // rmin, rmax + 10, halez + 10);
  // volBds->dump(std::cout);

  // CUTOUT CYLINDER VOLUME BOUNDS
  double rmin_at_center = std::numeric_limits<double>::max();
  double rmin_at_choke  = std::numeric_limits<double>::max();
  double rmax           = std::numeric_limits<double>::lowest();
  double zmin           = std::numeric_limits<double>::max();
  double zmax           = std::numeric_limits<double>::lowest();
  for (const auto& box : boxStore) {
    Acts::Vector3D vmin = box->min().cast<double>();
    Acts::Vector3D vmax = box->max().cast<double>();

    double vrmin = perp(vmin);
    double vrmax = perp(vmax);

    rmin_at_choke = std::min(rmin_at_choke, std::min(vrmin, vrmax));

    rmax = std::max(rmax, std::max(vrmin, vrmax));
    zmin = std::min(zmin, std::min(vmin.z(), vmax.z()));
    zmax = std::max(zmax, std::max(vmin.z(), vmax.z()));

    if (std::abs(vmin.z()) < 100) {
      rmin_at_center = std::min(vrmin, rmin_at_center);
    }
    if (std::abs(vmax.z()) < 100) {
      rmin_at_center = std::min(vrmax, rmin_at_center);
    }
  }

  double cutout_zmin_abs = std::numeric_limits<double>::max();
  // double cutout_zmax_abs = std::numeric_limits<double>::lowest();
  for (const auto& box : boxStore) {
    Acts::Vector3D vmin  = box->min().cast<double>();
    Acts::Vector3D vmax  = box->max().cast<double>();
    double         vrmin = perp(vmin);
    double         vrmax = perp(vmax);

    if (vrmin < rmin_at_center * 0.9) {
      cutout_zmin_abs = std::min(cutout_zmin_abs, std::abs(vmin.z()));
      // cutout_zmax_abs = std::max(cutout_zmax_abs, std::abs(vmin.z()));
    }
    if (vrmax < rmin_at_center * 0.9) {
      cutout_zmin_abs = std::min(cutout_zmin_abs, std::abs(vmax.z()));
      // cutout_zmax_abs = std::max(cutout_zmax_abs, std::abs(vmax.z()));
    }
  }

  double dz1 = (zmax - zmin) / 2.;
  // double dz2 = (cutout_zmax_abs - cutout_zmin_abs) / 2.;
  double dz2 = cutout_zmin_abs;

  std::cout << "rmin_at_center: " << rmin_at_center
            << " rmin at choke: " << rmin_at_choke;
  std::cout << " rmax: " << rmax << " zmin: " << zmin << " zmax: " << zmax;
  std::cout << " coutout_zmin_abs: " << cutout_zmin_abs << std::endl;

  std::shared_ptr<Acts::CutoutCylinderVolumeBounds> volBds = nullptr;
  volBds = std::make_shared<Acts::CutoutCylinderVolumeBounds>(
      rmin_at_choke, rmin_at_center, rmax, dz1, dz2);
  std::cout << *volBds << std::endl;

  for (const auto& srf : volBds->decomposeToSurfaces(nullptr)) {
    std::cout << "NORMAL: " << srf->normal().transpose() << " ";
    std::cout << "CTR: " << srf->center().transpose() << std::endl;
  }

  // return 0;

  // this cylinder vol bounds should fit right inside the cutout
  auto innerCylCtrBds
      = std::make_shared<Acts::CylinderVolumeBounds>(0, rmin_at_center, dz2);

  auto innerCyl = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(Acts::Transform3D::Identity()),
      innerCylCtrBds,
      nullptr,
      "ID");

  auto chokeCylBds = std::make_shared<Acts::CylinderVolumeBounds>(
      0, rmin_at_choke, (dz1 - dz2) / 2.);

  double chokePosZ = dz2 + (dz1 - dz2) / 2.;  // mid point of choke cylinder

  auto posChokeVol = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(
          Acts::Translation3D(Acts::Vector3D::UnitZ() * chokePosZ)),
      chokeCylBds,
      nullptr,
      "calo_choke_pos");

  auto negChokeVol = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(
          Acts::Translation3D(Acts::Vector3D::UnitZ() * -1 * chokePosZ)),
      chokeCylBds,
      nullptr,
      "calo_choke_neg");

  Acts::ply_helper<double> ply;
  std::ofstream            os("cyl.ply");
  innerCylCtrBds->draw(ply);
  chokeCylBds->draw(ply, posChokeVol->transform());
  chokeCylBds->draw(ply, negChokeVol->transform());
  os << ply;
  os.close();

  Box* top;
  top = Acts::make_octree(boxStore, prims, 1, 0.1);
  // tracking volume needs to store the box store, and the top vol
  auto tvTrf
      = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());

  // repack cells
  std::vector<std::unique_ptr<const Acts::Volume>> cells_vol;
  cells_vol.reserve(cells.size());
  for (auto& cell : cells) {
    cells_vol.push_back(std::unique_ptr<const Acts::Volume>(cell.release()));
  }
  std::shared_ptr<Acts::TrackingVolume> tv
      = Acts::TrackingVolume::create(std::move(tvTrf),
                                     volBds,
                                     std::move(boxStore),
                                     std::move(cells_vol),
                                     top,
                                     nullptr,  // no material
                                     "calo");

  tv->glueTrackingVolume(Acts::tubeInnerCover, innerCyl, Acts::tubeOuterCover);
  tv->glueTrackingVolume(Acts::index7, posChokeVol, Acts::tubeOuterCover);
  tv->glueTrackingVolume(Acts::index6, negChokeVol, Acts::tubeOuterCover);

  // create binutility for ID to Calo
  std::vector<float> rBoundaries   = {0, float(rmin_at_choke), float(rmax)};
  auto               binUtilityPos = std::make_unique<const Acts::BinUtility>(
      rBoundaries, Acts::open, Acts::binR);
  auto binUtilityNeg = std::make_unique<const Acts::BinUtility>(
      rBoundaries, Acts::open, Acts::binR);

  Acts::Vector3D caloBinPos(
      volBds->rMin() + (volBds->rMed() - volBds->rMin()) / 2., 0, 0);

  std::vector<Acts::TrackingVolumeOrderPosition> tVolOrderedPos;
  tVolOrderedPos.push_back(std::make_pair(tv, caloBinPos));
  tVolOrderedPos.push_back(std::make_pair(
      posChokeVol, Acts::Vector3D(rmin_at_choke / 2., 0, -chokePosZ)));

  std::vector<Acts::TrackingVolumeOrderPosition> tVolOrderedNeg;
  tVolOrderedNeg.push_back(std::make_pair(tv, caloBinPos));
  tVolOrderedNeg.push_back(std::make_pair(
      negChokeVol, Acts::Vector3D(rmin_at_choke / 2., 0, -chokePosZ)));

  auto tVolArrPos
      = std::make_shared<const Acts::BinnedArrayXD<Acts::TrackingVolumePtr>>(
          tVolOrderedPos, std::move(binUtilityPos));
  auto tVolArrNeg
      = std::make_shared<const Acts::BinnedArrayXD<Acts::TrackingVolumePtr>>(
          tVolOrderedNeg, std::move(binUtilityNeg));

  using bnd_srf_t = Acts::BoundarySurfaceT<Acts::TrackingVolume>;
  // now attach those to the ID volume
  std::const_pointer_cast<bnd_srf_t>(
      innerCyl->boundarySurfaces().at(Acts::positiveFaceXY))
      ->attachVolumeArray(tVolArrPos, Acts::outsideVolume);
  std::const_pointer_cast<bnd_srf_t>(
      innerCyl->boundarySurfaces().at(Acts::negativeFaceXY))
      ->attachVolumeArray(tVolArrNeg, Acts::outsideVolume);

  // the calo volume needs to point in
  std::const_pointer_cast<bnd_srf_t>(tv->boundarySurfaces().at(Acts::index5))
      ->attachVolume(innerCyl, Acts::outsideVolume);
  std::const_pointer_cast<bnd_srf_t>(tv->boundarySurfaces().at(Acts::index6))
      ->attachVolume(innerCyl, Acts::outsideVolume);

  std::const_pointer_cast<bnd_srf_t>(
      posChokeVol->boundarySurfaces().at(Acts::negativeFaceXY))
      ->attachVolume(innerCyl, Acts::outsideVolume);
  std::const_pointer_cast<bnd_srf_t>(
      negChokeVol->boundarySurfaces().at(Acts::positiveFaceXY))
      ->attachVolume(innerCyl, Acts::outsideVolume);

  // now we need to create a new container volume that contains everything

  // create three pseudo container, that contain the neg and pos choke, and the
  // center
  auto posContainer = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(posChokeVol->transform()),
      std::make_shared<Acts::CylinderVolumeBounds>(
          0, volBds->rMax(), (volBds->dZ1() - volBds->dZ2()) / 2.),
      tVolArrPos);
  auto negContainer = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(negChokeVol->transform()),
      std::make_shared<Acts::CylinderVolumeBounds>(
          0, volBds->rMax(), (volBds->dZ1() - volBds->dZ2()) / 2.),
      tVolArrNeg);

  std::vector<Acts::TrackingVolumeOrderPosition> tVolOrderedCtr;
  tVolOrderedCtr.push_back(
      std::make_pair(innerCyl, Acts::Vector3D(volBds->rMed() / 2., 0, 0)));
  tVolOrderedCtr.push_back(std::make_pair(
      tv,
      Acts::Vector3D(
          volBds->rMed() + (volBds->rMax() - volBds->rMed()) / 2., 0, 0)));

  std::vector<float> ctrBoundariesR
      = {0, float(volBds->rMed()), float(volBds->rMax())};
  auto binUtilityCtr = std::make_unique<const Acts::BinUtility>(
      ctrBoundariesR, Acts::open, Acts::binR);

  auto tVolArrCtr
      = std::make_shared<const Acts::BinnedArrayXD<Acts::TrackingVolumePtr>>(
          tVolOrderedCtr, std::move(binUtilityCtr));

  auto ctrContainer = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(Acts::Transform3D::Identity()),
      std::make_shared<Acts::CylinderVolumeBounds>(
          0, volBds->rMax(), volBds->dZ2()),
      tVolArrCtr);

  // and now combine those together into another one
  Acts::TrackingVolumeArrayCreator tvac;

  auto mainContainer = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(Acts::Transform3D::Identity()),
      std::make_shared<Acts::CylinderVolumeBounds>(
          0, volBds->rMax(), volBds->dZ1()),
      tvac.trackingVolumeArray({negContainer, ctrContainer, posContainer},
                               Acts::binZ));

  auto tg = std::make_shared<Acts::TrackingGeometry>(mainContainer);

  // draw it!
  // std::vector<std::shared_ptr<const Acts::Surface>> volBSrf
  //= tv->volumeBounds().decomposeToSurfaces(&tv->transform());

  // for (const auto& srf : volBSrf) {

  // const Acts::CylinderSurface* cyl
  //= dynamic_cast<const Acts::CylinderSurface*>(srf.get());
  // const Acts::DiscSurface* disc
  //= dynamic_cast<const Acts::DiscSurface*>(srf.get());
  // Acts::PolyhedronRepresentation poly({}, {});
  // if (cyl != nullptr) {
  // poly = cyl->polyhedronRepresentation(50);
  //} else if (disc != nullptr) {
  // poly = disc->polyhedronRepresentation(50);
  //} else {
  // throw std::runtime_error("Incompatible surface");
  //}

  // poly.draw(ply);
  //}

  // outer
  ply.clear();
  os = std::ofstream("ccyl.ply");
  volBds->draw(ply, tv->transform());
  os << ply;
  os.close();

  // return 0;

  // Acts::NavigationOptions<Acts::Surface> opt(Acts::forward, true);
  // auto surfaces = tv->compatibleSurfacesFromHierarchy(
  //{0, 0, 0}, Acts::Vector3D(0, 1500, 7000).normalized(), opt);

  // std::cout << "FOUND " << surfaces.size() << " SURFACES" << std::endl;
  // for(const auto& srf : surfaces) {
  // std::cout << *srf.object << std::endl;
  //}

  Acts::Vector3D origin(0, 0, 0);

  std::cout << "Testing " << n_rays << " rays..." << std::endl;

  ply.clear();

  std::mt19937                           rng(42);
  std::uniform_real_distribution<double> eta_dist(-5, 5);
  std::uniform_real_distribution<double> phi_dist(-M_PI, M_PI);
  std::uniform_real_distribution<double> z_dist(-100, 100);
  std::vector<Acts::Vector3D>            dirs;
  dirs.reserve(n_rays);
  for (size_t i = 0; i < n_rays; i++) {
    double eta = eta_dist(rng);
    // double         eta   = 2.0;
    double         phi   = phi_dist(rng);
    double         theta = 2 * std::atan(std::exp(-eta));
    Acts::Vector3D dir;
    dir << std::cos(phi), std::sin(phi), 1. / std::tan(theta);
    dir.normalize();
    dirs.push_back(std::move(dir));
  }

  using SteppingLogger    = Acts::detail::SteppingLogger;
  using BField_type       = Acts::ConstantBField;
  using EigenStepper_type = Acts::EigenStepper<BField_type>;
  using EigenPropagatorType
      = Acts::Propagator<EigenStepper_type, Acts::Navigator>;

  BField_type         bField(0, 0, 0);
  EigenStepper_type   stepper(bField);
  Acts::Navigator     navigator(tg);
  EigenPropagatorType propagator(std::move(stepper), navigator);

  using DebugOutput     = Acts::detail::DebugOutputActor;
  using ActionList      = Acts::ActionList<SteppingLogger, DebugOutput>;
  using AbortConditions = Acts::AbortList<>;

  // setup propagation options

  // std::ofstream os("compatibleSurfacesFromHierarchy_ray.ply");
  bool debug = true;

  std::ofstream step_os("propsteps.csv");
  step_os << "step_x,step_y,step_z,step_r,vol_id,bnd_id,lay_id,app_id,sen_id\n";

  clock::time_point start = clock::now();
  for (size_t i = 0; i < n_rays; i++) {
    const auto& dir = dirs[i];
    double      mom = 50 * Acts::units::_GeV;

    // Acts::Vector3D zshift(0, 0, z_dist(rng));
    Acts::Vector3D zshift(0, 0, 0);

    Acts::CurvilinearParameters startPar(
        nullptr, origin + zshift, dir * mom, +1);

    Acts::PropagatorOptions<ActionList, AbortConditions> options;
    options.debug     = debug;
    options.pathLimit = 20 * Acts::units::_m;

    const auto& result = propagator.propagate(startPar, options);

    const auto debugString
        = result.template get<DebugOutput::result_type>().debugString;
    // if (debug) { std::cout << debugString << std::endl; }
    // ply.line(origin, (origin + dir * 10000).eval());
    // os << ply;
    // ply.clear();

    // Acts::NavigationOptions<Acts::Surface> opt(Acts::forward, true);
    // auto sfis = tv->compatibleSurfacesFromHierarchy(origin, dir, opt);
    // for (const auto& sfi : sfis) {
    // const auto* pb
    //= dynamic_cast<const Acts::PlanarBounds*>(&sfi.object->bounds());
    // std::vector<Acts::Vector3D> vvtx;
    // for (const auto& vtx : pb->vertices()) {
    // Acts::Vector3D glob;
    // sfi.object->localToGlobal(vtx, {}, glob);
    // vvtx.push_back(glob);
    //}
    // ply.face(vvtx);
    // std::cout << "SRFIX: at:" << sfi.intersection.pathLength;
    // std::cout << " with: " << *sfi.object << std::endl;
    //}

    // get the steps!
    auto steppingResults
        = result.template get<SteppingLogger::result_type>().steps;
    using ag = Acts::GeometryID;

    auto last = steppingResults.back();
    if (last.position.z() < -7000 && perp(last.position) < 10000) {
      // print and end
      if (last.surface != nullptr) { std::cout << *last.surface << std::endl; }
      std::cout << last.position.transpose() << std::endl;
      std::cout << "DEBUG:" << std::endl;
      std::cout << debugString << std::endl;
      // break;
    }

    for (const auto& step : steppingResults) {
      geo_id_value volumeID    = 0;
      geo_id_value boundaryID  = 0;
      geo_id_value layerID     = 0;
      geo_id_value approachID  = 0;
      geo_id_value sensitiveID = 0;
      // get the identification from the surface first
      if (step.surface) {
        auto geoID  = step.surface->geoID();
        sensitiveID = geoID.value(ag::sensitive_mask);
        approachID  = geoID.value(ag::approach_mask);
        layerID     = geoID.value(ag::layer_mask);
        boundaryID  = geoID.value(ag::boundary_mask);
        volumeID    = geoID.value(ag::volume_mask);
      }
      // a current volume overwrites the surface tagged one
      if (step.volume) {
        volumeID = step.volume->geoID().value(ag::volume_mask);
      }
      // now fill
      step_os << step.position.x() << ",";
      step_os << step.position.y() << ",";
      step_os << step.position.z() << ",";
      step_os << Acts::VectorHelpers::perp(step.position) << ",";

      step_os << volumeID << "," << boundaryID << "," << layerID << ","
              << approachID << "," << sensitiveID << "\n";
    }
  }

  // os = std::ofstream("compatibleSurfacesFromHierarchy_surf.ply");
  // os << ply;
  // os.close();

  clock::time_point end = clock::now();
  double            duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();
  std::cout << "=> done in " << duration * 1e-6 << "s !" << std::endl;

  double t_per_ray = duration / n_rays;
  std::cout << " => time per ray: " << t_per_ray << "us" << std::endl;

  // auto hits = intersections(ray, top);

  // Acts::ply_helper<float> ply;
  // for (const auto& hit : hits) {
  // auto vol = hit->entity();
  // auto cvb = dynamic_cast<const Acts::GenericCuboidVolumeBounds*>(
  //&vol->volumeBounds());
  // cvb->draw(ply, vol->transform());
  //}

  // std::ofstream os("hits.ply");
  // os << ply;
  // os.close();

  // ply.clear();
  // ray.draw(ply, 10000);
  // os = std::ofstream("rays.ply");
  // os << ply;
  // os.close();
}

