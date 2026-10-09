// Gen3 (blueprint) build of the sPHENIX TGeo geometry, mirroring what the
// Gen1 TGeoDetector + tgeo-sphenix-mms json produce: four barrel subsystems
// (MVTX, INTT "Silicon", TPC, MICROMEGAS), stacked in R, each a stack of
// cylindrical layers built from sensor clusters.
//
// Usage: sphenix_gen3 <geometry.root> <out-prefix> [--zgaps] [loglevel]

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/BlueprintBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <ranges>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

using namespace Acts;
using namespace Acts::UnitLiterals;
using TGeoBuilder = ActsPlugins::BlueprintBuilder;
using Element = ActsPlugins::TGeoBlueprintBuilderBackend::Element;

namespace {

// Same sensitive volume names as the Gen1 json config
const std::regex kSensitive{
    "MVTXSensor.*|siactive.*|tpc_gas_measurement.*|micromegas_measurement.*"};

struct SubsystemSpec {
  std::string name;    // container name (= Gen1 volume name)
  std::regex sensors;  // sensor volume name pattern
  std::string axes;    // sensor axes (Gen1 "geo-tgeo-sensitive-axes")
  ActsPlugins::TGeoAxes axesDef;  // same, as compile-time checked TGeoAxes
  double rTolerance;              // merge tolerance for the radial clustering
  std::size_t matPhiBins;
  std::size_t matZBins;
};

/// Radial extent (min, max) of the sensor *surface* (2D plane spanned by
/// the first two axes letters), same as what ProtoLayer sees in Gen1.
std::pair<double, double> surfaceRExtent(const TGeoBuilder& builder,
                                         const Element& el,
                                         const std::string& axes) {
  const auto& backend = builder.backend();
  const auto* box =
      dynamic_cast<const TGeoBBox*>(backend.nodeOf(el).GetVolume()->GetShape());
  const auto m = backend.transformOf(el);
  auto half = [&](char a) {
    switch (std::toupper(a)) {
      case 'X':
        return box->GetDX();
      case 'Y':
        return box->GetDY();
      default:
        return box->GetDZ();
    }
  };
  auto idx = [](char a) { return std::toupper(a) - 'X'; };
  double rmin = 1e12, rmax = -1e12;
  for (int s0 : {-1, 1}) {
    for (int s1 : {-1, 1}) {
      double local[3] = {0, 0, 0};
      local[idx(axes[0])] = s0 * half(axes[0]);
      local[idx(axes[1])] = s1 * half(axes[1]);
      double global[3];
      m.LocalToMaster(local, global);
      const double r = std::hypot(global[0], global[1]) * 1_cm;
      rmin = std::min(rmin, r);
      rmax = std::max(rmax, r);
    }
  }
  return {rmin, rmax};
}

/// Transitive radial clustering of sensors by their surface extent: any two
/// sensors whose (tolerance-padded) r ranges overlap end up in the same
/// layer. This is what Gen1's ProtoLayerHelper is meant to do (and does after
/// the fix to merge bridged clusters). Returns sensors sorted by cluster and a
/// path -> layer key map.
std::pair<std::vector<Element>, std::map<std::string, std::string>> clusterInR(
    const TGeoBuilder& builder, std::vector<Element> sensors,
    const SubsystemSpec& spec, const Logger& log) {
  auto logger = [&]() -> const Logger& { return log; };
  struct Item {
    Element el;
    double rmin, rmax;
  };
  std::vector<Item> items;
  for (auto& el : sensors) {
    auto [rmin, rmax] = surfaceRExtent(builder, el, spec.axes);
    items.push_back({el, rmin, rmax});
  }
  std::ranges::sort(items, {}, &Item::rmin);

  std::vector<Element> sorted;
  std::map<std::string, std::string> keys;
  int layer = -1;
  double hi = -1e12, lo = 0;
  std::size_t count = 0;
  auto report = [&]() {
    if (layer >= 0) {
      ACTS_INFO(
          std::format("  {} layer {}: {} sensors, surface r [{:.2f}, "
                      "{:.2f}] mm",
                      spec.name, layer, count, lo, hi));
    }
  };
  for (auto& it : items) {
    if (it.rmin - spec.rTolerance > hi) {
      hi -= spec.rTolerance;
      report();
      ++layer;
      lo = it.rmin;
      count = 0;
    }
    hi = std::max(hi, it.rmax + spec.rTolerance);
    ++count;
    keys[builder.backend().pathOf(it.el)] = std::format("L{}", layer);
    sorted.push_back(it.el);
  }
  hi -= spec.rTolerance;
  report();
  return {std::move(sorted), std::move(keys)};
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: " << argv[0]
              << " <geometry.root> <out-prefix> [--zgaps] [loglevel]\n";
    return 1;
  }
  const std::string outPrefix = argv[2];
  bool zGaps = false;
  int logLevel = Logging::INFO;
  for (int i = 3; i < argc; ++i) {
    if (std::string_view{argv[i]} == "--zgaps") {
      zGaps = true;
    } else {
      logLevel = std::stoi(argv[i]);
    }
  }
  auto loggerPtr = getDefaultLogger("sPHENIXGen3", Logging::Level(logLevel));
  const Logger& LOGGER = *loggerPtr;
  auto logger = [&]() -> const Logger& { return LOGGER; };

  TGeoManager::Import(argv[1]);
  if (gGeoManager == nullptr) {
    return 1;
  }

  ActsPlugins::TGeoBlueprintBuilderBackend::Config cfg;
  cfg.root = gGeoManager->GetTopNode();
  cfg.lengthScale = 1_cm;  // Gen1 "geo-tgeo-unit-scalor": 10
  cfg.sensitivePredicate = [](const Element& el) {
    return std::regex_match(el.context->node->GetVolume()->GetName(),
                            kSensitive);
  };
  // The TPC "measurement" volumes are full gas slabs (~5.6 mm thick), and the
  // sensitive thickness enters the ProtoLayer extent, so neighbouring TPC
  // layers overlap in r (Gen1 tolerated that: its approach surfaces overlap).
  // Gen3 can't stack overlapping volumes, so treat them as thin planes.
  cfg.elementFactory =
      [](const ActsPlugins::TGeoDetectorElement::Identifier& id,
         const TGeoNode& node, const TGeoMatrix& matrix,
         ActsPlugins::TGeoAxes axes, double scale,
         std::shared_ptr<const ISurfaceMaterial> material) {
        auto de =
            ActsPlugins::TGeoBlueprintBuilderBackend::defaultElementFactory(
                id, node, matrix, axes, scale, std::move(material));
        if (std::string_view{node.GetVolume()->GetName()}.starts_with(
                "tpc_gas_measurement")) {
          de->surface().assignThickness(0.);
        }
        return de;
      };
  TGeoBuilder builder{cfg, LOGGER.clone("TGeoBlpBld")};

  const std::vector<SubsystemSpec> subsystems = {
      {"MVTX", std::regex{"MVTXSensor"}, "XZY", "XZY", 1_mm, 36, 50},
      {"Silicon", std::regex{"siactive_volume_.*"}, "YZX", "YZX", 1_mm, 36, 50},
      {"TPC", std::regex{"tpc_gas_measurement_.*"}, "YZX", "YZX", 1_mm, 36, 50},
      {"MICROMEGAS", std::regex{"micromegas_measurement_.*"}, "YZX", "YZX",
       1_mm, 36, 50},
  };

  Blueprint::Config bpCfg;
  bpCfg.envelope = ExtentEnvelope{{.z = {20_mm, 20_mm}, .r = {0_mm, 20_mm}}};
  Blueprint root{bpCfg};

  auto& world = root.addCylinderContainer("sPHENIX", AxisDirection::AxisR);
  world.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  world.setResizeStrategies(VolumeResizeStrategy::Gap,
                            VolumeResizeStrategy::Gap);

  // Inner volume down to r = 0 so the origin is inside the world (Gen1's MVTX
  // volume extends to r = 0). No beampipe is built, as in the Gen1 config
  // ("geo-tgeo-build-beampipe": false); the real Be pipe (r = 20.8 mm) would
  // overlap the off-axis MVTX surface extent anyway. Its z length is
  // synchronized by the R-stack.
  world.addChild(
      std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 15_mm, 100_mm),
          "BeamPipeRegion")));

  const ExtentEnvelope layerEnvelope{{.z = {1_mm, 1_mm}, .r = {1_mm, 1_mm}}};

  for (const auto& spec : subsystems) {
    auto sensors = builder.findDetElementByNamePattern(
        builder.backend().world(), spec.sensors);
    ACTS_INFO(spec.name << ": " << sensors.size() << " sensors");
    auto [sorted, keyMap] =
        clusterInR(builder, std::move(sensors), spec, LOGGER);
    const int nLayers = static_cast<int>(
        std::set<std::string>(std::from_range, keyMap | std::views::values)
            .size());

    auto container =
        builder.layersFromSensors()
            .barrel()
            .setSensorAxes(spec.axesDef)
            .setSensors(std::move(sorted))
            .groupBy([&builder, keys = keyMap](const Element& el) {
              return keys.at(builder.backend().pathOf(el));
            })
            .setContainerName(spec.name)
            .setEnvelope(layerEnvelope)
            .setAttachmentStrategy(VolumeAttachmentStrategy::Gap)
            .onLayer([&spec, &layerEnvelope, zGaps, nLayers](
                         const std::optional<Element>&,
                         std::shared_ptr<LayerBlueprintNode> layer)
                         -> std::shared_ptr<BlueprintNode> {
              // keep layers on the beam axis (sensor CoG is off-axis)
              layer->setUseCenterOfGravity(false, false, true);
              SurfaceArrayNavigationPolicy::Config navCfg;
              navCfg.layerType =
                  SurfaceArrayNavigationPolicy::LayerType::Cylinder;
              navCfg.bins = {0, 0};  // auto: one bin per module
              navCfg.envelope = layerEnvelope;
              layer->setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SurfaceArrayNavigationPolicy>(navCfg)
                      .asUniquePtr());

              // Gen1 approach surfaces <-> Gen3 inner/outer layer portals.
              // With --zgaps the subsystem's outermost faces are merged with
              // the z-gap volumes, and material can't be designated on a
              // partial merged face: those go on the subsystem level instead.
              using enum CylinderVolumeBounds::Face;
              const int idx = std::stoi(
                  layer->name().substr(layer->name().rfind("|L") + 2));
              std::vector<CylinderVolumeBounds::Face> faces;
              if (!zGaps || idx > 0) {
                faces.push_back(InnerCylinder);
              }
              if (!zGaps || idx + 1 < nLayers) {
                faces.push_back(OuterCylinder);
              }
              if (faces.empty()) {
                return layer;
              }
              auto mat = std::make_shared<MaterialDesignatorBlueprintNode>(
                  layer->name() + "_Mat");
              for (auto face : faces) {
                mat->configureFace(
                    face,
                    AxisSpec::DeferredEquidistant(spec.matPhiBins,
                                                  AxisDirection::AxisRPhi),
                    AxisSpec::DeferredEquidistant(spec.matZBins,
                                                  AxisDirection::AxisZ));
              }
              mat->addChild(std::move(layer));
              return mat;
            })
            .build();
    if (!zGaps) {
      // Default: the R-stack synchronizes z by *expanding* every layer to the
      // full world length, so all layer portals span the full z range.
      world.addChild(std::move(container));
      continue;
    }
    // R-stacks always synchronize z by expanding their children. Wrap each
    // subsystem in a z-container that absorbs the z-resize with gap volumes
    // instead, giving Gen1's {fGap | Barrel | sGap} structure.
    auto zWrap = std::make_shared<CylinderContainerBlueprintNode>(
        spec.name + "_Z", AxisDirection::AxisZ);
    zWrap->setResizeStrategies(VolumeResizeStrategy::Gap,
                               VolumeResizeStrategy::Gap);
    zWrap->addChild(std::move(container));
    // Material for the subsystem's innermost/outermost faces, designated on
    // the merged (full-z) portal.
    auto zMat =
        std::make_shared<MaterialDesignatorBlueprintNode>(spec.name + "_Mat");
    for (auto face : {CylinderVolumeBounds::Face::InnerCylinder,
                      CylinderVolumeBounds::Face::OuterCylinder}) {
      zMat->configureFace(
          face,
          AxisSpec::DeferredEquidistant(spec.matPhiBins,
                                        AxisDirection::AxisRPhi),
          AxisSpec::DeferredEquidistant(spec.matZBins, AxisDirection::AxisZ));
    }
    zMat->addChild(std::move(zWrap));
    world.addChild(std::move(zMat));
  }

  {
    std::ofstream dot{outPrefix + "_blueprint.dot"};
    root.graphviz(dot);
  }

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto tg = root.construct(BlueprintOptions{}, gctx, LOGGER);

  // ---- dump: volumes, sensitive counts, material-designated surfaces ----
  std::ofstream vols{outPrefix + "_volumes.csv"};
  vols << "vol,name,rmin,rmax,hz,cz,nsens\n";
  tg->apply([&](const TrackingVolume& v) {
    const auto& b = v.volumeBounds().values();
    std::size_t nsens = 0;
    for (const auto& s : v.surfaces()) {
      nsens += s.isSensitive() ? 1 : 0;
    }
    vols << std::format("{},{},{:.2f},{:.2f},{:.2f},{:.2f},{}\n",
                        v.geometryId().volume(), v.volumeName(), b[0], b[1],
                        b[2], v.localToGlobalTransform(gctx).translation().z(),
                        nsens);
  });

  std::ofstream mats{outPrefix + "_material_surfaces.csv"};
  mats << "geoid,vol,type,r0,r1,z0,z1\n";
  tg->visitSurfaces(
      [&](const Surface* s) {
        if (s->surfaceMaterial() == nullptr) {
          return;
        }
        const auto& t = s->localToGlobalTransform(gctx);
        const double cz = t.translation().z();
        std::string type;
        double r0 = 0, r1 = 0, z0 = cz, z1 = cz;
        if (const auto* cb =
                dynamic_cast<const CylinderBounds*>(&s->bounds())) {
          type = "cyl";
          r0 = r1 = cb->get(CylinderBounds::eR);
          z0 = cz - cb->get(CylinderBounds::eHalfLengthZ);
          z1 = cz + cb->get(CylinderBounds::eHalfLengthZ);
        } else if (const auto* rb =
                       dynamic_cast<const RadialBounds*>(&s->bounds())) {
          type = "disc";
          r0 = rb->rMin();
          r1 = rb->rMax();
        } else {
          type = "other";
        }
        std::ostringstream gid;
        gid << s->geometryId();
        mats << std::format("\"{}\",{},{},{:.2f},{:.2f},{:.2f},{:.2f}\n",
                            gid.str(), s->geometryId().volume(), type, r0, r1,
                            z0, z1);
      },
      false);

  // ---- straight-line navigation check: sensitive hits per subsystem ----
  {
    std::map<unsigned, std::string> volSubsystem;
    tg->apply([&](const TrackingVolume& v) {
      std::string n = v.volumeName();
      volSubsystem[v.geometryId().volume()] =
          n.substr(0, n.find_first_of("|:_"));
    });
    std::shared_ptr<const TrackingGeometry> tgShared = std::move(tg);
    Navigator::Config navCfg{tgShared};
    Navigator nav{navCfg, LOGGER.clone("Nav", Logging::WARNING)};
    Propagator prop{StraightLineStepper{}, std::move(nav)};
    using Actors = ActorList<SurfaceCollector<>, EndOfWorldReached>;
    auto mctx = MagneticFieldContext{};
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> uEta{-1.2, 1.2}, uPhi{-M_PI, M_PI};
    std::ofstream hits{outPrefix + "_hits.csv"};
    hits << "eta,phi,ok,MVTX,Silicon,TPC,MICROMEGAS\n";
    std::size_t nFail = 0;
    for (int i = 0; i < 5000; ++i) {
      const double eta = uEta(rng), phi = uPhi(rng);
      const double theta = 2 * std::atan(std::exp(-eta));
      auto start = BoundTrackParameters::createCurvilinear(
          Vector4::Zero(), phi, theta, 1. / 10_GeV, std::nullopt,
          ParticleHypothesis::pion());
      decltype(prop)::Options<Actors> opts{gctx, mctx};
      opts.pathLimit = 5_m;
      auto res = prop.propagate(start, opts);
      std::map<std::string, int> n;
      if (res.ok()) {
        for (const auto& h :
             res->get<SurfaceCollector<>::result_type>().collected) {
          ++n[volSubsystem[h.surface->geometryId().volume()]];
        }
      } else {
        ++nFail;
      }
      hits << std::format("{:.4f},{:.4f},{},{},{},{},{}\n", eta, phi,
                          res.ok() ? 1 : 0, n["MVTX"], n["Silicon"], n["TPC"],
                          n["MICROMEGAS"]);
    }
    ACTS_INFO("Propagation: " << nFail << " / 5000 failed");
  }

  ACTS_INFO("Done, wrote " << outPrefix << "_{blueprint.dot,volumes.csv,"
                           << "material_surfaces.csv}");
  return 0;
}
