// Gen3 (blueprint) build of the sPHENIX TGeo geometry, mirroring what the
// Gen1 TGeoDetector + tgeo-sphenix-mms json produce: four barrel subsystems
// (MVTX, INTT "Silicon", TPC, MICROMEGAS), stacked in R.
//
// MVTX, INTT and MICROMEGAS are single volumes navigated with the TryAll
// policy. Their physical layers overlap in r (MVTX is off the beam axis, INTT
// ladders are staggered), so they can't be separate layer volumes. Instead,
// each physical layer gets a passive cylinder "material carrier" surface in
// the volume (centred on the subsystem's own axis, extent from its sensors),
// like Gen1 approach surfaces but decoupled from the volume structure.
//
// The TPC is a stack of layer volumes with surface-array navigation, and
// material is designated on the layer portals.
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
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsPlugins/Root/BlueprintBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
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

const ExtentEnvelope kLayerEnvelope{{.z = {1_mm, 1_mm}, .r = {1_mm, 1_mm}}};
constexpr std::size_t kMatPhiBins = 36;
constexpr std::size_t kMatZBins = 50;

AxisSpec phiBins() {
  return AxisSpec::DeferredEquidistant(kMatPhiBins, AxisDirection::AxisRPhi);
}
AxisSpec zBins() {
  return AxisSpec::DeferredEquidistant(kMatZBins, AxisDirection::AxisZ);
}

/// Global corners (mm) of the sensor *surface*: the 2D plane spanned by the
/// first two axes letters, which is what the TGeo surface converter builds.
std::array<Vector3, 4> surfaceCorners(const TGeoBuilder& builder,
                                      const Element& el,
                                      const std::string& axes) {
  const auto& backend = builder.backend();
  const auto* box =
      dynamic_cast<const TGeoBBox*>(backend.nodeOf(el).GetVolume()->GetShape());
  const auto m = backend.transformOf(el);
  const double half[3] = {box->GetDX(), box->GetDY(), box->GetDZ()};
  auto idx = [](char a) { return std::toupper(a) - 'X'; };
  std::array<Vector3, 4> corners;
  std::size_t i = 0;
  for (int s0 : {-1, 1}) {
    for (int s1 : {-1, 1}) {
      double local[3] = {0, 0, 0};
      local[idx(axes[0])] = s0 * half[idx(axes[0])];
      local[idx(axes[1])] = s1 * half[idx(axes[1])];
      double global[3];
      m.LocalToMaster(local, global);
      corners[i++] = Vector3{global[0], global[1], global[2]} * 1_cm;
    }
  }
  return corners;
}

Vector3 centerOf(const TGeoBuilder& builder, const Element& el) {
  const double* t = builder.backend().transformOf(el).GetTranslation();
  return Vector3{t[0], t[1], t[2]} * 1_cm;
}

/// Transitive radial clustering of sensors by their surface extent: any two
/// sensors whose (tolerance-padded) r ranges overlap end up in the same
/// layer. This is what Gen1's ProtoLayerHelper is meant to do (and does after
/// the fix to merge bridged clusters). Returns sensors sorted by cluster and a
/// path -> layer key map.
std::pair<std::vector<Element>, std::map<std::string, std::string>>
clusterInR(const TGeoBuilder& builder, std::vector<Element> sensors,
           const std::string& axes, double tolerance) {
  struct Item {
    Element el;
    double rmin, rmax;
  };
  std::vector<Item> items;
  for (auto& el : sensors) {
    Item it{el, 1e12, -1e12};
    for (const auto& c : surfaceCorners(builder, el, axes)) {
      it.rmin = std::min(it.rmin, c.head<2>().norm());
      it.rmax = std::max(it.rmax, c.head<2>().norm());
    }
    items.push_back(it);
  }
  std::ranges::sort(items, {}, &Item::rmin);

  std::vector<Element> sorted;
  std::map<std::string, std::string> keys;
  int layer = -1;
  double hi = -1e12;
  for (auto& it : items) {
    if (it.rmin - tolerance > hi) {
      ++layer;
    }
    hi = std::max(hi, it.rmax + tolerance);
    keys[builder.backend().pathOf(it.el)] = std::format("L{}", layer);
    sorted.push_back(it.el);
  }
  return {std::move(sorted), std::move(keys)};
}

struct Carrier {
  std::string subsystem;
  std::string name;
  Vector2 axis;  // (x, y) of the cylinder axis
  double r = 0;
  double z0 = 1e12, z1 = -1e12;
  std::size_t nSensors = 0;
};

/// One passive cylinder per physical layer: radius = mean sensor-centre
/// distance from @p axis, z = sensor surface extent.
std::vector<Carrier> makeCarriers(
    const TGeoBuilder& builder, const std::string& subsystem,
    const std::vector<Element>& sensors, const std::string& axes,
    const Vector2& axis,
    const std::function<std::string(const Element&, double)>& physLayer) {
  std::map<std::string, Carrier> carriers;
  for (const auto& el : sensors) {
    const double r = (centerOf(builder, el).head<2>() - axis).norm();
    const std::string key = physLayer(el, r);
    auto& c = carriers
                  .try_emplace(key, Carrier{.subsystem = subsystem,
                                            .name = key,
                                            .axis = axis})
                  .first->second;
    c.r += r;
    ++c.nSensors;
    for (const auto& corner : surfaceCorners(builder, el, axes)) {
      c.z0 = std::min(c.z0, corner.z());
      c.z1 = std::max(c.z1, corner.z());
    }
  }
  std::vector<Carrier> result;
  for (auto& [key, c] : carriers) {
    c.r /= c.nSensors;
    result.push_back(c);
  }
  std::ranges::sort(result, {}, &Carrier::r);
  return result;
}

std::shared_ptr<Surface> carrierSurface(const Carrier& c) {
  Transform3 t = Transform3::Identity();
  t.translation() = Vector3{c.axis.x(), c.axis.y(), 0.5 * (c.z0 + c.z1)};
  auto srf = Surface::makeShared<CylinderSurface>(
      t, std::make_shared<CylinderBounds>(c.r, 0.5 * (c.z1 - c.z0) + 1_mm));
  srf->assignSurfaceMaterial(std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D{{phiBins(), zBins()}}, MappingType::Default,
      c.subsystem + "_" + c.name));
  return srf;
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
  const auto& backend = builder.backend();
  auto find = [&](const char* pattern) {
    return builder.findDetElementByNamePattern(backend.world(),
                                               std::regex{pattern});
  };

  Blueprint::Config bpCfg;
  bpCfg.envelope = ExtentEnvelope{{.z = {20_mm, 20_mm}, .r = {0_mm, 20_mm}}};
  Blueprint root{bpCfg};

  auto& top = root.addCylinderContainer("sPHENIX", AxisDirection::AxisR);
  top.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  top.setResizeStrategies(VolumeResizeStrategy::Gap, VolumeResizeStrategy::Gap);

  // Inner volume down to r = 0 so the origin is inside the world (Gen1's MVTX
  // volume extends to r = 0). No beampipe is built, as in the Gen1 config
  // ("geo-tgeo-build-beampipe": false); the real Be pipe (r = 20.8 mm) would
  // overlap the off-axis MVTX surface extent anyway. Its z length is
  // synchronized by the R-stack.
  top.addChild(std::make_shared<StaticBlueprintNode>(
      std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 15_mm, 100_mm),
          "BeamPipeRegion")));

  // Adds a subsystem node to the top-level R-stack. With --zgaps, it is
  // wrapped in a z-container that absorbs the z-synchronization of the
  // R-stack with gap volumes (Gen1's {fGap | Barrel | sGap}); otherwise the
  // R-stack expands every volume to the full world length.
  auto addSubsystem = [&](const std::string& name,
                          std::shared_ptr<BlueprintNode> node) {
    if (!zGaps) {
      top.addChild(std::move(node));
      return;
    }
    auto zWrap = std::make_shared<CylinderContainerBlueprintNode>(
        name + "_Z", AxisDirection::AxisZ);
    zWrap->setResizeStrategies(VolumeResizeStrategy::Gap,
                               VolumeResizeStrategy::Gap);
    zWrap->addChild(std::move(node));
    top.addChild(std::move(zWrap));
  };

  // A single TryAll-navigated volume holding all sensors of a subsystem plus
  // one passive material carrier per physical layer.
  auto addTryAllSubsystem =
      [&](const std::string& name, const std::vector<Element>& sensors,
          ActsPlugins::TGeoAxes axesDef, const std::string& axes,
          const Vector2& axis,
          const std::function<std::string(const Element&, double)>&
              physLayer) {
        auto layer = builder.layerFromSensors()
                         .barrel()
                         .setSensorAxes(axesDef)
                         .setSensors(sensors)
                         .setLayerName(name)
                         .setEnvelope(kLayerEnvelope)
                         .build();
        auto& layerNode = dynamic_cast<LayerBlueprintNode&>(*layer);
        auto surfaces = layerNode.surfaces();
        auto carriers =
            makeCarriers(builder, name, sensors, axes, axis, physLayer);
        for (const auto& c : carriers) {
          ACTS_INFO(std::format(
              "  {} carrier {}: {} sensors, r = {:.2f} mm around ({:.2f}, "
              "{:.2f}), z [{:.1f}, {:.1f}] mm",
              name, c.name, c.nSensors, c.r, c.axis.x(), c.axis.y(), c.z0,
              c.z1));
          surfaces.push_back(carrierSurface(c));
        }
        layerNode.setSurfaces(std::move(surfaces));
        // The sensor CoG is off-axis for the MVTX: keep the volume coaxial.
        layerNode.setUseCenterOfGravity(false, false, true);
        layerNode.setNavigationPolicyFactory(
            NavigationPolicyFactory{}
                .add<TryAllNavigationPolicy>(TryAllNavigationPolicy::Config{
                    .portals = true, .sensitives = true, .passives = true})
                .asUniquePtr());
        ACTS_INFO(name << ": " << sensors.size() << " sensors, "
                       << carriers.size() << " carriers (TryAll)");
        addSubsystem(name, layer);
      };

  // ---- MVTX: one volume, 3 carriers around the MVTX axis ----
  {
    const auto wrapper = builder.findDetElementByName("log_MVTX_Wrapper");
    const double* t = backend.transformOf(*wrapper).GetTranslation();
    const Vector2 mvtxAxis = Vector2{t[0], t[1]} * 1_cm;
    // Physical layers are well separated in radius around the MVTX axis:
    // cut between them.
    addTryAllSubsystem("MVTX", find("MVTXSensor"), "XZY", "XZY", mvtxAxis,
                       [](const Element&, double r) {
                         return r < 28_mm ? "L0" : (r < 36_mm ? "L1" : "L2");
                       });
  }

  // ---- INTT: one volume, 4 carriers (one per ladder layer) ----
  // Ladders are named ladder_<layer>_<type>_<phi>_<side>; the two staggered
  // sub-layers of a barrel overlap radially, which is why they can't be
  // separate volumes.
  addTryAllSubsystem(
      "Silicon", find("siactive_volume_.*"), "YZX", "YZX", Vector2::Zero(),
      [&](const Element& el, double) {
        const std::string ladder = backend.nodeOf(backend.parent(el)).GetName();
        return "L" + ladder.substr(std::string{"ladder_"}.size(), 1);
      });

  // ---- TPC: 48 layer volumes, surface arrays, material on layer portals ----
  {
    auto [sorted, keys] =
        clusterInR(builder, find("tpc_gas_measurement_.*"), "YZX", 1_mm);
    const int nLayers = static_cast<int>(
        std::set<std::string>(std::from_range, keys | std::views::values)
            .size());
    auto container =
        builder.layersFromSensors()
            .barrel()
            .setSensorAxes("YZX")
            .setSensors(std::move(sorted))
            .groupBy([&builder, keys = keys](const Element& el) {
              return keys.at(builder.backend().pathOf(el));
            })
            .setContainerName("TPC")
            .setEnvelope(kLayerEnvelope)
            .setAttachmentStrategy(VolumeAttachmentStrategy::Gap)
            .onLayer([&](const std::optional<Element>&,
                         std::shared_ptr<LayerBlueprintNode> layer)
                         -> std::shared_ptr<BlueprintNode> {
              SurfaceArrayNavigationPolicy::Config navCfg;
              navCfg.layerType =
                  SurfaceArrayNavigationPolicy::LayerType::Cylinder;
              navCfg.bins = {0, 0};  // auto: one bin per module
              navCfg.envelope = kLayerEnvelope;
              layer->setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SurfaceArrayNavigationPolicy>(navCfg)
                      .asUniquePtr());

              // Gen1 approach surfaces <-> inner/outer layer portals. With
              // --zgaps the stack's outermost faces are merged with the
              // z-gap volumes, where material can't be designated on a
              // partial face: those are designated on the subsystem instead.
              using enum CylinderVolumeBounds::Face;
              const int idx = std::stoi(
                  layer->name().substr(layer->name().rfind("|L") + 2));
              auto mat = std::make_shared<MaterialDesignatorBlueprintNode>(
                  layer->name() + "_Mat");
              if (!zGaps || idx > 0) {
                mat->configureFace(InnerCylinder, phiBins(), zBins());
              }
              if (!zGaps || idx + 1 < nLayers) {
                mat->configureFace(OuterCylinder, phiBins(), zBins());
              }
              mat->addChild(std::move(layer));
              return mat;
            })
            .build();
    ACTS_INFO("TPC: " << nLayers << " layer volumes (surface arrays)");
    if (zGaps) {
      auto mat = std::make_shared<MaterialDesignatorBlueprintNode>("TPC_Mat");
      for (auto face : {CylinderVolumeBounds::Face::InnerCylinder,
                        CylinderVolumeBounds::Face::OuterCylinder}) {
        mat->configureFace(face, phiBins(), zBins());
      }
      auto zWrap = std::make_shared<CylinderContainerBlueprintNode>(
          "TPC_Z", AxisDirection::AxisZ);
      zWrap->setResizeStrategies(VolumeResizeStrategy::Gap,
                                 VolumeResizeStrategy::Gap);
      zWrap->addChild(std::move(container));
      mat->addChild(std::move(zWrap));
      top.addChild(std::move(mat));
    } else {
      top.addChild(std::move(container));
    }
  }

  // ---- MICROMEGAS: one volume, 2 carriers (phi tiles inner, z tiles outer)
  addTryAllSubsystem(
      "MICROMEGAS", find("micromegas_measurement_.*"), "YZX", "YZX",
      Vector2::Zero(), [&](const Element& el, double) {
        return backend.nameOf(el).find("inner") != std::string::npos ? "phi"
                                                                     : "z";
      });

  {
    std::ofstream dot{outPrefix + "_blueprint.dot"};
    root.graphviz(dot);
  }

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  std::shared_ptr<const TrackingGeometry> tg =
      root.construct(BlueprintOptions{}, gctx, LOGGER);

  // ---- dumps for plotting ----
  std::map<unsigned, std::string> volSubsystem;
  {
    std::ofstream vols{outPrefix + "_volumes.csv"};
    vols << "vol,name,rmin,rmax,hz,cx,cy,cz,nsens,npassive\n";
    tg->apply([&](const TrackingVolume& v) {
      const auto& b = v.volumeBounds().values();
      std::size_t nsens = 0, npass = 0;
      for (const auto& s : v.surfaces()) {
        (s.isSensitive() ? nsens : npass) += 1;
      }
      const Vector3 c = v.localToGlobalTransform(gctx).translation();
      vols << std::format(
          "{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{},{}\n",
          v.geometryId().volume(), v.volumeName(), b[0], b[1], b[2], c.x(),
          c.y(), c.z(), nsens, npass);
      const std::string n = v.volumeName();
      volSubsystem[v.geometryId().volume()] =
          n.substr(0, n.find_first_of("|:_"));
    });
  }
  {
    std::ofstream sens{outPrefix + "_sensors.csv"};
    sens << "subsystem,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3\n";
    std::ofstream mats{outPrefix + "_material_surfaces.csv"};
    mats << "geoid,subsystem,kind,cx,cy,r,z0,z1\n";
    std::set<const Surface*> seen;
    tg->visitSurfaces(
        [&](const Surface* s) {
          if (!seen.insert(s).second) {
            return;
          }
          const auto& sub = volSubsystem[s->geometryId().volume()];
          if (s->isSensitive()) {
            const auto vtx = s->polyhedronRepresentation(gctx, 1).vertices;
            sens << sub;
            for (std::size_t i = 0; i < 4; ++i) {
              sens << std::format(",{:.2f},{:.2f},{:.2f}", vtx[i].x(),
                                  vtx[i].y(), vtx[i].z());
            }
            sens << "\n";
            return;
          }
          const auto* cb = dynamic_cast<const CylinderBounds*>(&s->bounds());
          if (s->surfaceMaterial() == nullptr || cb == nullptr) {
            return;
          }
          const Vector3 c = s->localToGlobalTransform(gctx).translation();
          const double hz = cb->get(CylinderBounds::eHalfLengthZ);
          std::ostringstream gid;
          gid << s->geometryId();
          mats << std::format(
              "\"{}\",{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}\n", gid.str(),
              sub, s->geometryId().boundary() != 0 ? "portal" : "carrier",
              c.x(), c.y(), cb->get(CylinderBounds::eR), c.z() - hz,
              c.z() + hz);
        },
        false);
  }

  // ---- straight-line navigation check: hits + material crossings ----
  {
    Navigator::Config navCfg{tg};
    Navigator nav{navCfg, LOGGER.clone("Nav", Logging::WARNING)};
    Propagator prop{StraightLineStepper{}, std::move(nav)};
    using Actors = ActorList<SurfaceCollector<>, EndOfWorldReached>;
    auto mctx = MagneticFieldContext{};
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> uEta{-1.2, 1.2}, uPhi{-M_PI, M_PI};
    const std::vector<std::string> subs = {"MVTX", "Silicon", "TPC",
                                           "MICROMEGAS"};
    std::ofstream hits{outPrefix + "_hits.csv"};
    hits << "eta,phi,ok";
    for (const auto& s : subs) {
      hits << "," << s << "," << s << "_mat";
    }
    hits << "\n";
    std::size_t nFail = 0;
    for (int i = 0; i < 5000; ++i) {
      const double eta = uEta(rng), phi = uPhi(rng);
      const double theta = 2 * std::atan(std::exp(-eta));
      auto start = BoundTrackParameters::createCurvilinear(
          Vector4::Zero(), phi, theta, 1. / 10_GeV, std::nullopt,
          ParticleHypothesis::pion());
      decltype(prop)::Options<Actors> opts{gctx, mctx};
      opts.pathLimit = 5_m;
      opts.actorList.get<SurfaceCollector<>>().selector.selectMaterial = true;
      auto res = prop.propagate(start, opts);
      std::map<std::string, int> nSens, nMat;
      if (res.ok()) {
        for (const auto& h :
             res->get<SurfaceCollector<>::result_type>().collected) {
          const auto& sub = volSubsystem[h.surface->geometryId().volume()];
          if (h.surface->isSensitive()) {
            ++nSens[sub];
          } else if (h.surface->surfaceMaterial() != nullptr) {
            ++nMat[sub];
          }
        }
      } else {
        ++nFail;
      }
      hits << std::format("{:.4f},{:.4f},{}", eta, phi, res.ok() ? 1 : 0);
      for (const auto& s : subs) {
        hits << "," << nSens[s] << "," << nMat[s];
      }
      hits << "\n";
    }
    ACTS_INFO("Propagation: " << nFail << " / 5000 failed");
  }

  ACTS_INFO("Done, wrote " << outPrefix << "_*");
  return 0;
}
