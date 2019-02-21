#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Utilities/Visualization.hpp"
#include "Acts/Volumes/AbstractVolume.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"
#include "Acts/Volumes/GenericCuboidVolumeBounds.hpp"
#include "Acts/Volumes/Volume.hpp"

#include "Acts/Detector/TrackingVolume.hpp"

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
      std::array<Acts::Vector3D, 8>({p1, p2, p3, p4, p5, p6, p7, p8}));
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
      std::array<Acts::Vector3D, 8>({p1, p2, p3, p4, p5, p6, p7, p8}));

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
      std::array<Acts::Vector3D, 8>({p1, p2, p3, p4, p5, p6, p7, p8}));
  Acts::AbstractVolume vol(std::move(globalToLocal), std::move(cubo));

  return vol;
};

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
    if (calosample >= 12 && calosample <= 20) {
      scale = 0.5;
    }

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
};

int
main()
{
  std::cout << "Build calo geometry..." << std::flush;
  std::vector<std::unique_ptr<Acts::AbstractVolume>> cells;
  cells = atlasCaloFactory("../output_geo.csv");
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

  using Box = Acts::AxisAlignedBoundingBox<Acts::Volume, float, 3>;
  using Ray = Acts::Ray<float, 3>;

  auto intersections = [](const auto& obj, const Box* top) {
    const Box*              lnode = top;
    std::vector<const Box*> hits;
    do {
      if (lnode->intersect(obj)) {

        if (lnode->hasEntity()) {
          // found primitive
          // check obb to limit false positivies
          auto obb = lnode->entity()->orientedBoundingBox();
          if (obb.intersect(obj.transformed(
                  lnode->entity()->transform().inverse().cast<float>()))) {
            hits.push_back(lnode);
          }
          // we skip in any case, whether we actually hit the OBB or not
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

  // create BVH for the calo geo
  std::vector<std::unique_ptr<Box>> boxStore;
  std::vector<Box*>                 prims;

  double zmin = std::numeric_limits<double>::max();
  double zmax = std::numeric_limits<double>::lowest();
  double rmin = std::numeric_limits<double>::max();
  double rmax = -1;
  for (const auto& cell : cells) {
    boxStore.push_back(
        std::make_unique<Box>(cell->boundingBox({0.1, 0.1, 0.1})));
    prims.push_back(boxStore.back().get());

    Acts::Vector3D vmin = boxStore.back()->min().cast<double>();
    Acts::Vector3D vmax = boxStore.back()->max().cast<double>();

    zmin = std::min(zmin, vmin.z());
    zmin = std::min(zmin, vmax.z());
    zmax = std::max(zmax, vmin.z());
    zmax = std::max(zmax, vmax.z());

    rmin = std::min(rmin, Acts::VectorHelpers::perp(vmin));
    rmin = std::min(rmin, Acts::VectorHelpers::perp(vmax));
    rmax = std::max(rmax, Acts::VectorHelpers::perp(vmin));
    rmax = std::max(rmax, Acts::VectorHelpers::perp(vmax));
  }

  Box* top;
  top = Acts::make_octree(boxStore, prims, 6, 0.1);
  // tracking volume needs to store the box store, and the top vol
  auto tvTrf
      = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());

  // the cylinder volume bounds for the TV need to wrap around all the bounding
  // boxe assuming this is symmetric right now
  double halez   = (zmin + zmax) / 2.;
  auto cylVolBds = std::make_shared<Acts::CylinderVolumeBounds>(0, rmax, halez);

  std::shared_ptr<Acts::TrackingVolume> tv
      = Acts::TrackingVolume::create(std::move(tvTrf),
                                     cylVolBds,
                                     std::move(boxStore),
                                     top,
                                     nullptr,  // no material
                                     "calo");

  Ray ray({0, 0, 0}, {1, 1, 1});
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
