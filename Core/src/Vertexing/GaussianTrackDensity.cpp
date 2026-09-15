// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/GaussianTrackDensity.hpp"

#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory_resource>
#include <numbers>

namespace Acts {

// Keep the per-candidate calculation visible to the compiler in the hot loop.
inline void
Acts::GaussianTrackDensity::GaussianTrackDensityStore::addTrackToDensity(
    const TrackEntry& entry) {
  // Take track only if it's within bounds
  if (entry.lowerBound < m_z && m_z < entry.upperBound) {
    double delta = std::exp(entry.c0 + m_z * (entry.c1 + m_z * entry.c2));
    double qPrime = entry.c1 + 2. * m_z * entry.c2;
    double deltaPrime = delta * qPrime;
    m_density += delta;
    m_firstDerivative += deltaPrime;
    m_secondDerivative += 2. * entry.c2 * delta + qPrime * deltaPrime;
  }
}

// Cache complete sums, never subtract contributions from a previous sum.
struct GaussianTrackDensity::EvaluationCache {
  using Value = std::tuple<double, double, double>;
  std::pmr::unsynchronized_pool_resource pool;
  std::vector<TrackEntry> previous;
  std::pmr::map<double, Value> values;
  std::size_t capacity = 0;

  EvaluationCache() : values(&pool) {}
  EvaluationCache(const EvaluationCache& other)
      : previous(other.previous),
        values(other.values, &pool),
        capacity(other.capacity) {}

  static bool same(const TrackEntry& a, const TrackEntry& b) {
    const auto bits = [](double x) { return std::bit_cast<std::uint64_t>(x); };
    return bits(a.z) == bits(b.z) && bits(a.c0) == bits(b.c0) &&
           bits(a.c1) == bits(b.c1) && bits(a.c2) == bits(b.c2) &&
           bits(a.lowerBound) == bits(b.lowerBound) &&
           bits(a.upperBound) == bits(b.upperBound);
  }

  void prepare(const std::vector<TrackEntry>& entries) {
    // Only a subsequence preserves the accumulation order. Any addition,
    // coefficient change or reordering conservatively resets the cache.
    std::size_t next = 0;
    for (const auto& old : previous) {
      if (next < entries.size() && same(old, entries[next])) {
        ++next;
      } else if (old.lowerBound < old.upperBound) {
        values.erase(values.upper_bound(old.lowerBound),
                     values.lower_bound(old.upperBound));
      }
    }
    capacity = 8 * entries.size();
    if (next != entries.size() || values.size() > capacity) {
      values.clear();
    }
    previous = entries;
  }
};

// Index support intervals, not density values: every query still evaluates the
// original Gaussian and its derivatives with the original strict bounds and
// summation order. A fixed maximum number of bins bounds storage by O(nTracks).
struct GaussianTrackDensity::DensityIndex {
  explicit DensityIndex(const State& state) {
    if (state.trackEntries.size() < 32) {
      return;
    }
    minZ = std::numeric_limits<double>::infinity();
    maxZ = -std::numeric_limits<double>::infinity();
    for (const auto& entry : state.trackEntries) {
      if (std::isfinite(entry.z)) {
        minZ = std::min(minZ, entry.z);
        maxZ = std::max(maxZ, entry.z);
      }
    }
    width = maxZ - minZ;
    if (!(width > 0.) || !std::isfinite(width)) {
      return;
    }
    bins.resize(64);
    for (std::size_t i = 0; i < state.trackEntries.size(); ++i) {
      const auto& entry = state.trackEntries[i];
      if (!(entry.lowerBound < entry.upperBound) || entry.upperBound <= minZ ||
          entry.lowerBound >= maxZ) {
        continue;
      }
      const auto first = bin(entry.lowerBound);
      const auto last = bin(entry.upperBound);
      for (auto j = first; j <= last; ++j) {
        bins[j].push_back(i);
      }
    }
  }

  std::size_t bin(double z) const {
    // Clamping also handles infinite support bounds. Use the same monotonic
    // mapping for construction and queries so rounding at bin edges cannot
    // omit a contributing interval.
    const double fraction = std::clamp((z - minZ) / width, 0., 1.);
    return std::min(static_cast<std::size_t>(fraction * bins.size()),
                    bins.size() - 1);
  }

  const std::vector<std::size_t>* query(double z) const {
    // Refinement steps can leave the range of the initial trial positions.
    // Evaluate those uncommon queries with the exhaustive scan.
    if (bins.empty() || !std::isfinite(z) || z < minZ || z > maxZ) {
      return nullptr;
    }
    return &bins[bin(z)];
  }

  double minZ = 0.;
  double maxZ = 0.;
  double width = 0.;
  std::vector<std::vector<std::size_t>> bins;
};

Result<std::optional<std::pair<double, double>>>
Acts::GaussianTrackDensity::globalMaximumWithWidth(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto result = addTracks(state, trackList);
  if (!result.ok()) {
    return result.error();
  }

  if (state.cacheQueries) {
    if (!state.evaluationCache) {
      state.evaluationCache = std::make_shared<EvaluationCache>();
    } else if (!state.evaluationCache.unique()) {
      // Copying a State must not share subsequent cache mutations.
      state.evaluationCache =
          std::make_shared<EvaluationCache>(*state.evaluationCache);
    }
    state.evaluationCache->prepare(state.trackEntries);
  }
  DensityIndex index(state);

  double maxPosition = 0.;
  double maxDensity = 0.;
  double maxSecondDerivative = 0.;

  for (const auto& track : state.trackEntries) {
    double trialZ = track.z;

    auto [density, firstDerivative, secondDerivative] =
        trackDensityAndDerivatives(state, index, trialZ);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);

    trialZ += stepSize(density, firstDerivative, secondDerivative);
    std::tie(density, firstDerivative, secondDerivative) =
        trackDensityAndDerivatives(state, index, trialZ);

    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
    trialZ += stepSize(density, firstDerivative, secondDerivative);
    std::tie(density, firstDerivative, secondDerivative) =
        trackDensityAndDerivatives(state, index, trialZ);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    std::tie(maxPosition, maxDensity, maxSecondDerivative) =
        updateMaximum(trialZ, density, secondDerivative, maxPosition,
                      maxDensity, maxSecondDerivative);
  }

  if (maxSecondDerivative == 0.) {
    return std::nullopt;
  }

  return std::pair{maxPosition, std::sqrt(-(maxDensity / maxSecondDerivative))};
}

Result<std::optional<double>> Acts::GaussianTrackDensity::globalMaximum(
    State& state, const std::vector<InputTrack>& trackList) const {
  auto maxRes = globalMaximumWithWidth(state, trackList);
  if (!maxRes.ok()) {
    return maxRes.error();
  }
  const auto& maxOpt = *maxRes;
  if (!maxOpt.has_value()) {
    return std::nullopt;
  }
  return maxOpt->first;
}

Result<void> Acts::GaussianTrackDensity::addTracks(
    State& state, const std::vector<InputTrack>& trackList) const {
  for (auto trk : trackList) {
    const BoundTrackParameters& boundParams = m_cfg.extractParameters(trk);
    // Get required track parameters
    const double d0 = boundParams.parameters()[BoundIndices::eBoundLoc0];
    const double z0 = boundParams.parameters()[BoundIndices::eBoundLoc1];
    // Get track covariance
    if (!boundParams.covariance().has_value()) {
      return VertexingError::NoCovariance;
    }
    const auto perigeeCov = *(boundParams.covariance());
    const double covDD =
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc0);
    const double covZZ =
        perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1);
    const double covDZ =
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1);
    const double covDeterminant = (perigeeCov.block<2, 2>(0, 0)).determinant();

    // Do track selection based on track cov matrix and m_cfg.d0SignificanceCut
    if ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) ||
        (covZZ <= 0) || (covDeterminant <= 0)) {
      continue;
    }

    // Calculate track density quantities
    double constantTerm =
        -(d0 * d0 * covZZ + z0 * z0 * covDD + 2. * d0 * z0 * covDZ) /
        (2. * covDeterminant);
    const double linearTerm =
        (d0 * covDZ + z0 * covDD) /
        covDeterminant;  // minus signs and factors of 2 cancel...
    const double quadraticTerm = -covDD / (2. * covDeterminant);
    double discriminant =
        linearTerm * linearTerm -
        4. * quadraticTerm * (constantTerm + 2. * m_cfg.z0SignificanceCut);
    if (discriminant < 0) {
      continue;
    }

    // Add the track to the current maps in the state
    discriminant = std::sqrt(discriminant);
    const double zMax = (-linearTerm - discriminant) / (2. * quadraticTerm);
    const double zMin = (-linearTerm + discriminant) / (2. * quadraticTerm);
    constantTerm -= std::log(2. * std::numbers::pi * std::sqrt(covDeterminant));

    state.trackEntries.emplace_back(z0, constantTerm, linearTerm, quadraticTerm,
                                    zMin, zMax);
  }
  return Result<void>::success();
}

std::tuple<double, double, double>
Acts::GaussianTrackDensity::trackDensityAndDerivatives(const State& state,
                                                       DensityIndex& index,
                                                       double z) const {
  auto* cache = state.evaluationCache.get();
  // A map treats +0 and -0 as equal; keep their arithmetic independent.
  const bool cacheable = cache != nullptr && std::isfinite(z) && z != 0.;
  if (cacheable) {
    if (auto it = cache->values.find(z); it != cache->values.end()) {
      return it->second;
    }
  }
  GaussianTrackDensityStore densityResult(z);
  if (const auto* candidates = index.query(z)) {
    for (const auto i : *candidates) {
      densityResult.addTrackToDensity(state.trackEntries[i]);
    }
  } else {
    for (const auto& entry : state.trackEntries) {
      densityResult.addTrackToDensity(entry);
    }
  }
  auto result = densityResult.densityAndDerivatives();
  if (cacheable && cache->values.size() < cache->capacity) {
    cache->values.emplace(z, result);
  }
  return result;
}

std::tuple<double, double, double> Acts::GaussianTrackDensity::updateMaximum(
    double newZ, double newValue, double newSecondDerivative, double maxZ,
    double maxValue, double maxSecondDerivative) const {
  if (newValue > maxValue) {
    maxZ = newZ;
    maxValue = newValue;
    maxSecondDerivative = newSecondDerivative;
  }
  return {maxZ, maxValue, maxSecondDerivative};
}

double Acts::GaussianTrackDensity::stepSize(double y, double dy,
                                            double ddy) const {
  return (m_cfg.isGaussianShaped ? (y * dy) / (dy * dy - y * ddy) : -dy / ddy);
}

}  // namespace Acts
