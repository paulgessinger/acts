// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Portal-crossing validation and timing harness for Gen3 navigation.
//
// Each navigated propagation is checked against a fine-step RK4 replay of the
// trajectory it actually followed (step by step, from the recorded stepper
// states), so the reference is independent of the navigator, of the field
// model and of trajectory changes applied by actors.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {

// ---------------------------------------------------------------- fields

/// Curl- and divergence-free solenoid-like field with axial fall-off.
class VaryingField final : public MagneticFieldProvider {
 public:
  VaryingField(double b0, double zScale)
      : m_b0(b0), m_a(1. / (2. * zScale * zScale)) {}
  Cache makeCache(const MagneticFieldContext& /*mctx*/) const override {
    return Cache{std::in_place_type<int>, 0};
  }
  Result<Vector3> getField(const Vector3& p, Cache& /*cache*/) const override {
    return Result<Vector3>::success(field(p));
  }
  Vector3 field(const Vector3& p) const {
    const double r2 = p.x() * p.x() + p.y() * p.y();
    const double bz = m_b0 * (1. - m_a * p.z() * p.z() + 0.5 * m_a * r2);
    const double brOverR = m_b0 * m_a * p.z();
    return {brOverR * p.x(), brOverR * p.y(), bz};
  }

 private:
  double m_b0;
  double m_a;
};

struct FieldModel {
  std::shared_ptr<const MagneticFieldProvider> provider;
  std::function<Vector3(const Vector3&)> eval;
};

FieldModel makeField(const std::string& name) {
  if (name == "const") {
    Vector3 b{0, 0, 2_T};
    return {std::make_shared<ConstantBField>(b),
            [b](const Vector3&) { return b; }};
  }
  if (name == "zero") {
    Vector3 b{0, 0, 0};
    return {std::make_shared<ConstantBField>(b),
            [b](const Vector3&) { return b; }};
  }
  if (name == "tilted") {
    Vector3 b{0.25_T, -0.15_T, 2_T};
    return {std::make_shared<ConstantBField>(b),
            [b](const Vector3&) { return b; }};
  }
  if (name == "varying") {
    auto f = std::make_shared<VaryingField>(2_T, 3000_mm);
    return {f, [f](const Vector3& p) { return f->field(p); }};
  }
  throw std::invalid_argument("unknown field " + name);
}

// ---------------------------------------------------------------- geometry

struct VolInfo {
  const TrackingVolume* vol{};
  std::string name;
  double rMin{}, rMax{}, zMin{}, zMax{};
  int depth{};
};

struct Geo {
  std::shared_ptr<const TrackingGeometry> tg;
  std::vector<VolInfo> vols;
  std::unordered_map<const TrackingVolume*, int> index;
  std::vector<double> radii;
  std::vector<double> planes;
  double worldR{}, worldZ{};

  int idx(const TrackingVolume* v) const {
    if (v == nullptr) {
      return -1;
    }
    return index.at(v);
  }

  /// Deepest volume that contains the point, -1 outside of the world.
  int classify(const Vector3& p) const {
    const double r = std::hypot(p.x(), p.y());
    int best = -1;
    int bestDepth = -1;
    for (std::size_t i = 0; i < vols.size(); ++i) {
      const VolInfo& v = vols[i];
      if (r >= v.rMin && r <= v.rMax && p.z() >= v.zMin && p.z() <= v.zMax &&
          v.depth > bestDepth) {
        best = static_cast<int>(i);
        bestDepth = v.depth;
      }
    }
    return best;
  }
};

void dedupe(std::vector<double>& v) {
  std::ranges::sort(v);
  std::vector<double> out;
  for (double x : v) {
    if (out.empty() || x - out.back() > 1e-7) {
      out.push_back(x);
    }
  }
  v = std::move(out);
}

Geo loadGeometry(const GeometryContext& gctx, const std::string& path) {
  Geo geo;
  TrackingGeometryJsonConverter converter(
      TrackingGeometryJsonConverter::Config::defaultConfig(),
      getDefaultLogger("Converter", Logging::WARNING));
  geo.tg = converter.fromFile(gctx, path);
  geo.tg->visitVolumes([&](const TrackingVolume* vol) {
    const auto* cyl =
        dynamic_cast<const CylinderVolumeBounds*>(&vol->volumeBounds());
    if (cyl == nullptr) {
      throw std::runtime_error("non-cylinder volume " + vol->volumeName());
    }
    const Transform3& trf = vol->localToGlobalTransform(gctx);
    if (!trf.linear().isApprox(SquareMatrix3::Identity()) ||
        std::hypot(trf.translation().x(), trf.translation().y()) > 1e-6) {
      throw std::runtime_error("volume not coaxial " + vol->volumeName());
    }
    VolInfo info;
    info.vol = vol;
    info.name = vol->volumeName();
    info.rMin = cyl->get(CylinderVolumeBounds::eMinR);
    info.rMax = cyl->get(CylinderVolumeBounds::eMaxR);
    const double hl = cyl->get(CylinderVolumeBounds::eHalfLengthZ);
    info.zMin = trf.translation().z() - hl;
    info.zMax = trf.translation().z() + hl;
    geo.index[vol] = static_cast<int>(geo.vols.size());
    geo.vols.push_back(info);
  });
  // Nesting depth from containment; Gen3 volumes do not link their mother.
  for (VolInfo& v : geo.vols) {
    for (const VolInfo& o : geo.vols) {
      if (&o != &v && o.rMin <= v.rMin && o.rMax >= v.rMax &&
          o.zMin <= v.zMin && o.zMax >= v.zMax &&
          (o.rMin < v.rMin || o.rMax > v.rMax || o.zMin < v.zMin ||
           o.zMax > v.zMax)) {
        ++v.depth;
      }
    }
  }
  for (const VolInfo& v : geo.vols) {
    if (v.rMin > 0) {
      geo.radii.push_back(v.rMin);
    }
    geo.radii.push_back(v.rMax);
    geo.planes.push_back(v.zMin);
    geo.planes.push_back(v.zMax);
    if (v.depth == 0) {
      geo.worldR = v.rMax;
      geo.worldZ = v.zMax;
    }
  }
  dedupe(geo.radii);
  dedupe(geo.planes);
  return geo;
}

// ---------------------------------------------------------------- tracks

struct Start {
  Vector3 pos;
  Vector3 dir;
  double qOverP{};
  double pathLimit{};
  Direction direction = Direction::Forward();
};

/// RK4 on dr/ds = t, dt/ds = q/p * t x B.
struct Rk4 {
  const std::function<Vector3(const Vector3&)>* field{};
  double qOverP{};
  void step(Vector3& p, Vector3& t, double h) const {
    auto f = [&](const Vector3& pos, const Vector3& dir) {
      return Vector3(qOverP * dir.cross((*field)(pos)));
    };
    const Vector3 k1 = f(p, t);
    const Vector3 k2 =
        f(p + 0.5 * h * t + 0.125 * h * h * k1, t + 0.5 * h * k1);
    const Vector3 k3 =
        f(p + 0.5 * h * t + 0.125 * h * h * k1, t + 0.5 * h * k2);
    const Vector3 k4 = f(p + h * t + 0.5 * h * h * k3, t + h * k3);
    p += h * t + h * h / 6. * (k1 + k2 + k3);
    t += h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
    t.normalize();
  }
};

std::vector<Start> makeStarts(const std::string& category, std::size_t n,
                              std::uint64_t seed, const Geo& geo,
                              const FieldModel& field) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uni(0., 1.);
  std::normal_distribution<double> gauss(0., 1.);
  auto logUniform = [&](double lo, double hi) {
    return lo * std::pow(hi / lo, uni(rng));
  };
  auto direction = [&](double etaMax) {
    const double eta = (2. * uni(rng) - 1.) * etaMax;
    const double phi = (2. * uni(rng) - 1.) * std::numbers::pi;
    const double theta = 2. * std::atan(std::exp(-eta));
    return Vector3(std::sin(theta) * std::cos(phi),
                   std::sin(theta) * std::sin(phi), std::cos(theta));
  };
  auto qop = [&](const Vector3& dir, double pt) {
    const double p = pt / std::hypot(dir.x(), dir.y());
    const double q = uni(rng) < 0.5 ? -1. : 1.;
    return q / p;
  };

  std::vector<Start> out;
  out.reserve(n);
  const bool backward = category == "backward";
  for (std::size_t i = 0; i < n; ++i) {
    Start s;
    if (category == "prompt" || backward) {
      s.pos = Vector3(0.01 * gauss(rng), 0.01 * gauss(rng), 50. * gauss(rng));
      s.dir = direction(4.);
      s.qOverP = qop(s.dir, logUniform(0.5_GeV, 50_GeV));
      s.pathLimit = 4_m;
    } else if (category == "displaced") {
      const double r = 400. * uni(rng);
      const double phi = 2. * std::numbers::pi * uni(rng);
      s.pos = Vector3(r * std::cos(phi), r * std::sin(phi),
                      (2. * uni(rng) - 1.) * 1000.);
      s.dir = direction(3.);
      s.qOverP = qop(s.dir, logUniform(0.5_GeV, 20_GeV));
      s.pathLimit = 4_m;
    } else if (category == "stress") {
      const bool displaced = uni(rng) < 0.5;
      const double r = displaced ? 300. * uni(rng) : 0.01 * uni(rng);
      const double phi = 2. * std::numbers::pi * uni(rng);
      s.pos =
          Vector3(r * std::cos(phi), r * std::sin(phi),
                  displaced ? (2. * uni(rng) - 1.) * 800. : 50. * gauss(rng));
      s.dir = direction(2.5);
      s.qOverP = qop(s.dir, logUniform(0.1_GeV, 0.6_GeV));
      s.pathLimit = 6_m;
    } else {
      throw std::invalid_argument("unknown category " + category);
    }
    if (backward) {
      // Fly out with the reference integrator, then come back.
      Rk4 rk{&field.eval, s.qOverP};
      Vector3 p = s.pos;
      Vector3 t = s.dir;
      const double target = 300. + 1700. * uni(rng);
      double len = 0.;
      while (len < target) {
        Vector3 pn = p;
        Vector3 tn = t;
        rk.step(pn, tn, 0.5);
        if (std::hypot(pn.x(), pn.y()) > geo.worldR - 5. ||
            std::abs(pn.z()) > geo.worldZ - 5.) {
          break;
        }
        p = pn;
        t = tn;
        len += 0.5;
      }
      s.pos = p;
      s.dir = t;
      s.pathLimit = len + 200.;
      s.direction = Direction::Backward();
    }
    out.push_back(s);
  }
  return out;
}

// ---------------------------------------------------------------- recording

struct Sample {
  double s{};   // signed stepper path length on arrival
  Vector3 pos;  // arrival state
  Vector3 dir;
  Vector3 depPos;  // departure state (after a kick)
  Vector3 depDir;
  double depQop{};
  int vol{-1};
  std::uint64_t surface{};
  bool portal{};
  bool sensitive{};
  bool kicked{};
};

struct KickConfig {
  bool enabled = false;
  double angle = 2e-3;
  double eloss = 3e-3;
  double shift = 0.02;
};

struct Recorder {
  struct this_result {
    std::vector<Sample> samples;
    std::uint64_t kickState = 0;
  };
  using result_type = this_result;

  const Geo* geo{};
  KickConfig kick{};
  std::uint64_t seed{};

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, result_type& result,
                   const Logger& /*logger*/) const {
    if (state.stage == PropagatorStage::postPropagation) {
      return Result<void>::success();
    }
    Sample sample;
    sample.s = stepper.pathLength(state.stepping);
    sample.pos = stepper.position(state.stepping);
    sample.dir = stepper.direction(state.stepping);
    sample.vol = geo->idx(navigator.currentVolume(state.navigation));
    const Surface* surface = navigator.currentSurface(state.navigation);
    if (surface != nullptr) {
      sample.surface = surface->geometryId().value();
      sample.portal = surface->geometryId().boundary() != 0;
      sample.sensitive = surface->geometryId().sensitive() != 0;
    }
    sample.depPos = sample.pos;
    sample.depDir = sample.dir;
    sample.depQop = stepper.qOverP(state.stepping);

    if (kick.enabled && sample.sensitive &&
        state.stage == PropagatorStage::postStep) {
      // Deterministic pseudo-random kick: scattering, energy loss and an
      // in-plane position update, as material and a filter step would do.
      std::mt19937_64 rng(seed * 7919 + (++result.kickState));
      std::uniform_real_distribution<double> uni(-1., 1.);
      const Vector3 normal =
          surface->normal(state.options.geoContext, sample.pos, sample.dir);
      Vector3 axis = Vector3(uni(rng), uni(rng), uni(rng));
      axis -= axis.dot(sample.dir) * sample.dir;
      axis.normalize();
      const double angle = kick.angle * uni(rng);
      sample.depDir = (AngleAxis3(angle, axis) * sample.dir).normalized();
      Vector3 shift = kick.shift * Vector3(uni(rng), uni(rng), uni(rng));
      shift -= shift.dot(normal) * normal;
      sample.depPos = sample.pos + shift;
      sample.depQop *= 1. + kick.eloss * std::abs(uni(rng));
      stepper.update(state.stepping, sample.depPos, sample.depDir,
                     sample.depQop, stepper.time(state.stepping));
      sample.kicked = true;
    }
    result.samples.push_back(sample);
    return Result<void>::success();
  }
};

struct Transition {
  double u{};
  int from{-1};
  int to{-1};
  bool portal{};
};

// ---------------------------------------------------------------- reference

struct Hermite {
  Vector3 p0, m0, p1, m1;
  Vector3 at(double t) const {
    const double t2 = t * t;
    const double t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * p0 + (t3 - 2 * t2 + t) * m0 +
           (-2 * t3 + 3 * t2) * p1 + (t3 - t2) * m1;
  }
};

template <typename F>
double bisect(F&& f, double a, double b) {
  double fa = f(a);
  for (int i = 0; i < 60; ++i) {
    const double m = 0.5 * (a + b);
    const double fm = f(m);
    if ((fa < 0) == (fm < 0)) {
      a = m;
      fa = fm;
    } else {
      b = m;
    }
  }
  return 0.5 * (a + b);
}

struct Reference {
  std::vector<Transition> transitions;
  double maxReplayDeviation = 0.;
  bool kickRetraced = false;
};

Reference replay(const Geo& geo, const FieldModel& field,
                 const std::vector<Sample>& samples, double sign,
                 double fineStep) {
  Reference ref;
  if (samples.empty()) {
    return ref;
  }
  int current = geo.classify(samples.front().depPos);
  double uMax = sign * samples.front().s;
  double uLastKick = -std::numeric_limits<double>::max();

  auto roots = [&](const Hermite& h, std::vector<double>& out) {
    auto rAt = [&](double t) {
      const Vector3 p = h.at(t);
      return std::hypot(p.x(), p.y());
    };
    auto zAt = [&](double t) { return h.at(t).z(); };
    // The radius can turn within the interval, split it at the turning point
    std::array<double, 3> knots{0., 1., 1.};
    std::size_t nKnots = 2;
    auto radial = [&](double t) {
      const Vector3 p = h.at(t);
      const Vector3 d =
          h.at(std::min(1., t + 1e-6)) - h.at(std::max(0., t - 1e-6));
      return p.x() * d.x() + p.y() * d.y();
    };
    if ((radial(0.) < 0) != (radial(1.) < 0)) {
      knots = {0., bisect(radial, 0., 1.), 1.};
      nKnots = 3;
    }
    for (std::size_t k = 0; k + 1 < nKnots; ++k) {
      const double ta = knots[k];
      const double tb = knots[k + 1];
      const double ra = rAt(ta);
      const double rb = rAt(tb);
      auto lo = std::ranges::lower_bound(geo.radii, std::min(ra, rb));
      auto hi = std::ranges::upper_bound(geo.radii, std::max(ra, rb));
      for (auto it = lo; it != hi; ++it) {
        const double R = *it;
        if (ra != R && (ra - R) * (rb - R) <= 0) {
          out.push_back(bisect([&](double t) { return rAt(t) - R; }, ta, tb));
        }
      }
    }
    const double za = zAt(0.);
    const double zb = zAt(1.);
    auto lo = std::ranges::lower_bound(geo.planes, std::min(za, zb));
    auto hi = std::ranges::upper_bound(geo.planes, std::max(za, zb));
    for (auto it = lo; it != hi; ++it) {
      const double Z = *it;
      if (za != Z && (za - Z) * (zb - Z) <= 0) {
        out.push_back(bisect([&](double t) { return zAt(t) - Z; }, 0., 1.));
      }
    }
    std::ranges::sort(out);
  };

  // Collect the covered trajectory and all boundary-surface crossings first.
  // The owning volume is then classified in the middle between consecutive
  // crossings, well away from the boundaries themselves.
  struct Segment {
    double u0{}, u1{};
    Hermite h;
  };
  std::vector<Segment> segments;
  std::vector<double> crossings;
  const double uStart = uMax;
  static_cast<void>(current);

  std::vector<double> ts;
  for (std::size_t i = 0; i + 1 < samples.size(); ++i) {
    const Sample& a = samples[i];
    const Sample& b = samples[i + 1];
    if (a.kicked) {
      uLastKick = sign * a.s;
    }
    const double ds = b.s - a.s;
    if (ds == 0.) {
      continue;
    }
    if (sign * b.s < uLastKick - 1e-9) {
      ref.kickRetraced = true;
    }
    Rk4 rk{&field.eval, a.depQop};
    Vector3 p = a.depPos;
    Vector3 t = a.depDir;
    const int n =
        std::max(1, static_cast<int>(std::ceil(std::abs(ds) / fineStep)));
    const double h = ds / n;
    for (int j = 0; j < n; ++j) {
      const Vector3 p0 = p;
      const Vector3 t0 = t;
      rk.step(p, t, h);
      if (j + 1 == n) {
        // Close the step on the recorded arrival state, so the reference
        // stays continuous across steps.
        ref.maxReplayDeviation =
            std::max(ref.maxReplayDeviation, (p - b.pos).norm());
        p = b.pos;
        t = b.dir;
      }
      const double u0 = sign * (a.s + j * h);
      const double u1 = j + 1 == n ? sign * b.s : sign * (a.s + (j + 1) * h);
      if (u1 <= uMax) {
        continue;  // retracing an already covered part
      }
      Hermite herm{p0, h * t0, p, h * t};
      ts.clear();
      roots(herm, ts);
      if (std::getenv("NAVCHECK_DEBUG") != nullptr && std::abs(ds) < 1e-2) {
        std::cout.precision(12);
        std::cout << "tiny u0=" << u0 << " u1=" << u1 << " uMax=" << uMax
                  << " r0=" << std::hypot(p0.x(), p0.y())
                  << " r1=" << std::hypot(p.x(), p.y())
                  << " roots=" << ts.size() << "\n";
      }
      for (double tk : ts) {
        const double uk = u0 + tk * (u1 - u0);
        // A portal arrival can sit within an ulp of the boundary, which puts
        // the crossing at the very start of the following step.
        if (uk >= uMax) {
          crossings.push_back(uk);
        }
      }
      segments.push_back({u0, u1, herm});
      uMax = u1;
    }
  }
  if (segments.empty()) {
    return ref;
  }
  auto positionAt = [&](double u) {
    auto it = std::ranges::partition_point(
        segments, [&](const Segment& seg) { return seg.u1 < u; });
    if (it == segments.end()) {
      --it;
    }
    return it->h.at(std::clamp((u - it->u0) / (it->u1 - it->u0), 0., 1.));
  };
  std::ranges::sort(crossings);
  double uPrev = uStart;
  int volume = geo.classify(
      positionAt(crossings.empty() ? 0.5 * (uStart + uMax)
                                   : 0.5 * (uStart + crossings.front())));
  for (std::size_t k = 0; k < crossings.size(); ++k) {
    const double uNext = k + 1 < crossings.size() ? crossings[k + 1] : uMax;
    const int next = geo.classify(positionAt(0.5 * (crossings[k] + uNext)));
    if (std::getenv("NAVCHECK_DEBUG") != nullptr) {
      const Vector3 q = positionAt(0.5 * (crossings[k] + uNext));
      std::cout.precision(10);
      std::cout << "crossing u=" << crossings[k]
                << " mid r=" << std::hypot(q.x(), q.y()) << " z=" << q.z()
                << " -> "
                << (next >= 0 ? geo.vols[static_cast<std::size_t>(next)].name
                              : "none")
                << " depth="
                << (next >= 0 ? geo.vols[static_cast<std::size_t>(next)].depth
                              : -1)
                << "\n";
    }
    if (next != volume) {
      ref.transitions.push_back({crossings[k], volume, next, true});
      volume = next;
    }
    uPrev = crossings[k];
  }
  static_cast<void>(uPrev);
  return ref;
}

std::vector<Transition> navigated(const std::vector<Sample>& samples,
                                  double sign) {
  std::vector<Transition> out;
  for (std::size_t i = 1; i < samples.size(); ++i) {
    if (samples[i].vol != samples[i - 1].vol) {
      out.push_back({sign * samples[i].s, samples[i - 1].vol, samples[i].vol,
                     samples[i].portal});
    }
  }
  return out;
}

// ---------------------------------------------------------------- main

struct Args {
  std::map<std::string, std::string> kv;
  std::string get(const std::string& k, const std::string& def) const {
    auto it = kv.find(k);
    return it == kv.end() ? def : it->second;
  }
  double num(const std::string& k, double def) const {
    auto it = kv.find(k);
    return it == kv.end() ? def : std::stod(it->second);
  }
};

struct Summary {
  std::size_t tracks = 0, failed = 0, steps = 0, renavigations = 0;
  std::size_t expected = 0, matched = 0, late = 0, missing = 0, spurious = 0;
  std::size_t worldExpected = 0, worldMatched = 0;
  std::size_t affectedTracks = 0, hits = 0, kickRetraced = 0, backSteps = 0;
  double maxDeviation = 0., maxBackStep = 0.;
};

template <typename stepper_t>
int run(const Args& args, stepper_t stepper, const Geo& geo,
        const FieldModel& field, const GeometryContext& gctx,
        const MagneticFieldContext& mctx) {
  Navigator::Config navCfg;
  navCfg.trackingGeometry = geo.tg;
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator(
      navCfg,
      getDefaultLogger("Navigator", Logging::Level(static_cast<int>(args.num(
                                        "navlog", Logging::WARNING)))));
  Propagator<stepper_t, Navigator> propagator(
      std::move(stepper), std::move(navigator),
      getDefaultLogger("Propagator", Logging::Level(static_cast<int>(args.num(
                                         "proplog", Logging::WARNING)))));

  const std::string category = args.get("category", "prompt");
  const std::size_t n = static_cast<std::size_t>(args.num("n", 1000));
  const std::uint64_t seed = static_cast<std::uint64_t>(args.num("seed", 1));
  const std::vector<Start> starts = makeStarts(category, n, seed, geo, field);
  const double maxStep = args.num("maxstep", 100.);
  const bool withCov = args.num("cov", 0) != 0;
  const long only = static_cast<long>(args.num("only", -1));

  auto makeParams = [&](const Start& s) {
    std::optional<BoundMatrix> cov;
    if (withCov) {
      cov = BoundMatrix::Identity() * 1e-4;
    }
    Vector4 pos4;
    pos4 << s.pos, 0.;
    return BoundTrackParameters::createCurvilinear(pos4, s.dir, s.qOverP, cov,
                                                   ParticleHypothesis::pion());
  };
  auto configure = [&](auto& options, const Start& s) {
    options.direction = s.direction;
    options.pathLimit = s.pathLimit;
    options.stepping.maxStepSize = maxStep;
    options.maxSteps = 100000;
#ifdef NAVCHECK_HAVE_BACKTRACK
    options.navigation.portalBacktracking = args.num("fix", 1) != 0;
#endif
#ifdef NAVCHECK_HAVE_GUARD
    options.navigation.experimentalCylindricalPortalGuard =
        args.num("fix", 1) != 0;
#endif
  };

  // ------------------------------------------------------------ timing
  if (args.get("mode", "check") == "bench") {
    using Actors = ActorList<EndOfWorldReached>;
    using Options =
        typename Propagator<stepper_t, Navigator>::template Options<Actors>;
    std::vector<BoundTrackParameters> params;
    for (const Start& s : starts) {
      params.push_back(makeParams(s));
    }
    const int batches = static_cast<int>(args.num("batches", 10));
    for (int b = 0; b < batches + 1; ++b) {
      std::size_t steps = 0;
      double path = 0.;
      const auto t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; i < starts.size(); ++i) {
        Options options(gctx, mctx);
        configure(options, starts[i]);
        auto res = propagator.propagate(params[i], options, false);
        if (res.ok()) {
          steps += res.value().steps;
          path += res.value().pathLength;
        }
      }
      const auto t1 = std::chrono::steady_clock::now();
      if (b == 0) {
        continue;  // warm-up
      }
      std::cout << std::chrono::duration<double, std::micro>(t1 - t0).count() /
                       static_cast<double>(starts.size())
                << " " << steps << " " << path << "\n";
    }
    return 0;
  }

  // ------------------------------------------------------------ checking
  using Actors = ActorList<Recorder, EndOfWorldReached>;
  using Options =
      typename Propagator<stepper_t, Navigator>::template Options<Actors>;

  const double fineStep = args.num("finestep", 0.25);
  const double tol = args.num("tol", 1e-3);
  KickConfig kick;
  kick.enabled = args.num("kick", 0) != 0;

  std::ofstream hitsOut;
  if (auto p = args.get("hits", ""); !p.empty()) {
    hitsOut.open(p);
  }
  std::ofstream detailOut;
  if (auto p = args.get("detail", ""); !p.empty()) {
    detailOut.open(p);
    detailOut.precision(10);
  }

  Summary sum;
  for (std::size_t i = 0; i < starts.size(); ++i) {
    if (only >= 0 && static_cast<long>(i) != only) {
      continue;
    }
    const Start& start = starts[i];
    Options options(gctx, mctx);
    configure(options, start);
    auto& recorder = options.actorList.template get<Recorder>();
    recorder.geo = &geo;
    recorder.kick = kick;
    recorder.seed = seed * 1000003 + i;

    ++sum.tracks;
    auto state = propagator.makeState(options);
    const BoundTrackParameters startParams = makeParams(start);
    auto init = propagator.initialize(state, startParams, nullptr);
    bool ok = init.ok();
    if (ok) {
      ok = propagator.propagate(state).ok();
    }
    if (!ok) {
      ++sum.failed;
      if (detailOut.is_open()) {
        detailOut << i << " failed steps=" << state.steps << "\n";
      }
    }
    const auto& samples = state.template get<Recorder::result_type>().samples;
    sum.steps += state.steps;
    sum.renavigations += state.navigation.statistics.nRenavigations;

    const double sign = start.direction == Direction::Forward() ? 1. : -1.;
    const Reference ref = replay(geo, field, samples, sign, fineStep);
    const std::vector<Transition> nav = navigated(samples, sign);
    sum.maxDeviation = std::max(sum.maxDeviation, ref.maxReplayDeviation);
    sum.kickRetraced += ref.kickRetraced ? 1 : 0;
    for (std::size_t k = 1; k < samples.size(); ++k) {
      const double du = sign * (samples[k].s - samples[k - 1].s);
      if (du < -tol) {
        ++sum.backSteps;
        sum.maxBackStep = std::max(sum.maxBackStep, -du);
      }
    }

    // The navigator does not get to act on the final step, unless it ended on
    // a portal
    const double uEnd = samples.empty() ? 0.
                        : (samples.size() > 1 && !samples.back().portal)
                            ? sign * samples[samples.size() - 2].s
                            : sign * samples.back().s;
    std::vector<bool> used(nav.size(), false);
    bool affected = false;
    for (const Transition& exp : ref.transitions) {
      if (exp.u > uEnd - tol && exp.to != -1) {
        continue;  // beyond where the propagation stopped
      }
      const bool world = exp.to == -1;
      (world ? sum.worldExpected : sum.expected) += 1;
      int best = -1;
      for (std::size_t k = 0; k < nav.size(); ++k) {
        if (used[k] || nav[k].from != exp.from || nav[k].to != exp.to) {
          continue;
        }
        if (best < 0 ||
            std::abs(nav[k].u - exp.u) <
                std::abs(nav[static_cast<std::size_t>(best)].u - exp.u)) {
          best = static_cast<int>(k);
        }
      }
      std::string verdict;
      if (best >= 0 &&
          std::abs(nav[static_cast<std::size_t>(best)].u - exp.u) <= tol &&
          nav[static_cast<std::size_t>(best)].portal) {
        used[static_cast<std::size_t>(best)] = true;
        (world ? sum.worldMatched : sum.matched) += 1;
        continue;
      }
      if (world) {
        continue;
      }
      affected = true;
      if (best >= 0) {
        used[static_cast<std::size_t>(best)] = true;
        ++sum.late;
        verdict = "late";
      } else {
        ++sum.missing;
        verdict = "missing";
      }
      if (detailOut.is_open()) {
        const double delta =
            best >= 0 ? nav[static_cast<std::size_t>(best)].u - exp.u : 0.;
        detailOut << i << " " << verdict << " u=" << exp.u << " "
                  << geo.vols[static_cast<std::size_t>(exp.from)].name << " -> "
                  << (exp.to >= 0
                          ? geo.vols[static_cast<std::size_t>(exp.to)].name
                          : "none")
                  << " delta=" << delta << " pt="
                  << std::hypot(start.dir.x(), start.dir.y()) /
                         std::abs(start.qOverP)
                  << "\n";
      }
    }
    for (std::size_t k = 0; k < nav.size(); ++k) {
      if (!used[k] && nav[k].to != -1 && nav[k].u <= uEnd - tol) {
        ++sum.spurious;
        affected = true;
        if (detailOut.is_open()) {
          detailOut
              << i << " spurious u=" << nav[k].u << " "
              << (nav[k].from >= 0
                      ? geo.vols[static_cast<std::size_t>(nav[k].from)].name
                      : "none")
              << " -> " << geo.vols[static_cast<std::size_t>(nav[k].to)].name
              << " portal=" << nav[k].portal << "\n";
        }
      }
    }
    sum.affectedTracks += affected ? 1 : 0;

    for (const Sample& s : samples) {
      if (s.sensitive) {
        ++sum.hits;
        if (hitsOut.is_open()) {
          hitsOut << i << " " << s.surface << "\n";
        }
      }
    }
    if (only >= 0) {
      std::cout.precision(10);
      const std::size_t first = samples.size() > 400 ? samples.size() - 60 : 0;
      for (std::size_t k = first; k < samples.size(); ++k) {
        const Sample& s = samples[k];
        std::cout << "s=" << s.s << " r=" << std::hypot(s.pos.x(), s.pos.y())
                  << " z=" << s.pos.z() << " vol="
                  << (s.vol >= 0
                          ? geo.vols[static_cast<std::size_t>(s.vol)].name
                          : "none")
                  << (s.portal ? " PORTAL" : "") << (s.sensitive ? " SENS" : "")
                  << "\n";
      }
      for (const Transition& t : ref.transitions) {
        std::cout << "ref u=" << t.u << " "
                  << (t.from >= 0
                          ? geo.vols[static_cast<std::size_t>(t.from)].name
                          : "none")
                  << " -> "
                  << (t.to >= 0 ? geo.vols[static_cast<std::size_t>(t.to)].name
                                : "none")
                  << "\n";
      }
    }
  }

  std::cout << "{\"category\":\"" << category << "\",\"field\":\""
            << args.get("field", "const") << "\",\"stepper\":\""
            << args.get("stepper", "eigen") << "\",\"kick\":" << kick.enabled
            << ",\"fix\":" << args.num("fix", 1) << ",\"tracks\":" << sum.tracks
            << ",\"failed\":" << sum.failed << ",\"steps\":" << sum.steps
            << ",\"renavigations\":" << sum.renavigations
            << ",\"expected\":" << sum.expected
            << ",\"matched\":" << sum.matched << ",\"late\":" << sum.late
            << ",\"missing\":" << sum.missing
            << ",\"spurious\":" << sum.spurious
            << ",\"affectedTracks\":" << sum.affectedTracks
            << ",\"worldExpected\":" << sum.worldExpected
            << ",\"worldMatched\":" << sum.worldMatched
            << ",\"hits\":" << sum.hits << ",\"backSteps\":" << sum.backSteps
            << ",\"maxBackStep\":" << sum.maxBackStep
            << ",\"kickRetraced\":" << sum.kickRetraced
            << ",\"maxReplayDeviation\":" << sum.maxDeviation << "}\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const auto eq = a.find('=');
    if (eq == std::string::npos) {
      std::cerr << "expected key=value, got " << a << "\n";
      return 1;
    }
    args.kv[a.substr(0, eq)] = a.substr(eq + 1);
  }
  const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext mctx;
  const Geo geo = loadGeometry(gctx, args.get("geometry", ""));
  const FieldModel field = makeField(args.get("field", "const"));

  if (args.get("stepper", "eigen") == "sympy") {
    return run(args, SympyStepper(field.provider), geo, field, gctx, mctx);
  }
  return run(args, EigenStepper<>(field.provider), geo, field, gctx, mctx);
}
