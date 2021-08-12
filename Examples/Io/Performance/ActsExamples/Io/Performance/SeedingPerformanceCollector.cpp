// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "SeedingPerformanceCollector.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <stdexcept>
#include <unordered_map>

#include <TFile.h>

namespace {
using SimParticleContainer = ActsExamples::SimParticleContainer;
using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;
using ProtoTrackContainer = ActsExamples::ProtoTrackContainer;
}  // namespace

ActsExamples::SeedingPerformanceCollector::SeedingPerformanceCollector(
    ActsExamples::SeedingPerformanceCollector::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputProtoTracks, "SeedingPerformanceCollector", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
}

ActsExamples::SeedingPerformanceCollector::~SeedingPerformanceCollector() {}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceCollector::endRun() {
  float eff = float(m_nTotalMatchedParticles) / m_nTotalParticles;
  float fakeRate = float(m_nTotalSeeds - m_nTotalMatchedSeeds) / m_nTotalSeeds;
  float duplicationRate =
      float(m_nTotalDuplicatedParticles) / m_nTotalMatchedParticles;
  float aveNDuplicatedSeeds =
      float(m_nTotalMatchedSeeds - m_nTotalMatchedParticles) /
      m_nTotalMatchedParticles;
  ACTS_DEBUG("nTotalSeeds               = " << m_nTotalSeeds);
  ACTS_DEBUG("nTotalMatchedSeeds        = " << m_nTotalMatchedSeeds);
  ACTS_DEBUG("nTotalParticles           = " << m_nTotalParticles);
  ACTS_DEBUG("nTotalMatchedParticles    = " << m_nTotalMatchedParticles);
  ACTS_DEBUG("nTotalDuplicatedParticles = " << m_nTotalDuplicatedParticles);

  ACTS_INFO("Efficiency (nMatchedParticles / nAllParticles) = " << eff);
  ACTS_INFO("Fake rate (nUnMatchedSeeds / nAllSeeds) = " << fakeRate);
  ACTS_INFO(
      "Duplication rate (nDuplicatedMatchedParticles / nMatchedParticles) = "
      << duplicationRate);
  ACTS_INFO(
      "Average number of duplicated seeds ((nMatchedSeeds - nMatchedParticles) "
      "/ nMatchedParticles) = "
      << aveNDuplicatedSeeds);

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceCollector::writeT(
    const AlgorithmContext& ctx, const ProtoTrackContainer& tracks) {
  // Read truth information collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  size_t nSeeds = tracks.size();
  size_t nMatchedSeeds = 0;
  // Map from particles to how many times they were successfully found by a seed
  std::unordered_map<ActsFatras::Barcode, std::size_t> truthCount;

  for (size_t itrack = 0; itrack < tracks.size(); ++itrack) {
    const auto& track = tracks[itrack];
    std::vector<ParticleHitCount> particleHitCounts;
    identifyContributingParticles(hitParticlesMap, track, particleHitCounts);
    // All hits matched to the same particle
    if (particleHitCounts.size() == 1) {
      auto prt = particleHitCounts.at(0);
      auto it = truthCount.try_emplace(prt.particleId, 0u).first;
      it->second += 1;
      nMatchedSeeds++;
    }
  }

  int nMatchedParticles = 0;
  int nDuplicatedParticles = 0;
  // Fill the effeciency and fake rate plots
  for (const auto& particle : particles) {
    const auto it1 = truthCount.find(particle.particleId());
    int nMatchedSeedsForParticle = 0;
    if (it1 != truthCount.end()) {
      nMatchedParticles++;
      nMatchedSeedsForParticle = it1->second;
      if (nMatchedSeedsForParticle > 1) {
        nDuplicatedParticles++;
      }
    }
  }
  ACTS_DEBUG("Number of seeds: " << nSeeds);
  m_nTotalSeeds += nSeeds;
  m_nTotalMatchedSeeds += nMatchedSeeds;
  m_nTotalParticles += particles.size();
  m_nTotalMatchedParticles += nMatchedParticles;
  m_nTotalDuplicatedParticles += nDuplicatedParticles;

  return ProcessCode::SUCCESS;
}

size_t ActsExamples::SeedingPerformanceCollector::nTotalSeeds() const {
  return m_nTotalSeeds;
}
size_t ActsExamples::SeedingPerformanceCollector::nTotalMatchedSeeds() const {
  return m_nTotalMatchedSeeds;
}
size_t ActsExamples::SeedingPerformanceCollector::nTotalParticles() const {
  return m_nTotalParticles;
}
size_t ActsExamples::SeedingPerformanceCollector::nTotalMatchedParticles()
    const {
  return m_nTotalMatchedParticles;
}
size_t ActsExamples::SeedingPerformanceCollector::nTotalDuplicatedParticles()
    const {
  return m_nTotalDuplicatedParticles;
}