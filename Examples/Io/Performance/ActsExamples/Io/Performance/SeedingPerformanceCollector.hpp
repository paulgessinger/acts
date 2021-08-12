// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"

#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {
class SeedingPerformanceCollector final : public WriterT<ProtoTrackContainer> {
 public:
  struct Config {
    /// Input reconstructed proto tracks collection.
    std::string inputProtoTracks;
    /// Input hit to particles map.
    std::string inputMeasurementParticlesMap;
    /// Input truth particles collection.
    std::string inputParticles;
  };

  /// Construct from configuration and log level.
  /// @param config The configuration
  /// @param level
  SeedingPerformanceCollector(Config config, Acts::Logging::Level level);

  ~SeedingPerformanceCollector() override;

  /// Finalize plots.
  ProcessCode endRun() final override;

  size_t nTotalSeeds() const;
  size_t nTotalMatchedSeeds() const;
  size_t nTotalParticles() const;
  size_t nTotalMatchedParticles() const;
  size_t nTotalDuplicatedParticles() const;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ProtoTrackContainer& tracks) final override;

  Config m_cfg;

  std::atomic<size_t> m_nTotalSeeds = 0;
  std::atomic<size_t> m_nTotalMatchedSeeds = 0;
  std::atomic<size_t> m_nTotalParticles = 0;
  std::atomic<size_t> m_nTotalMatchedParticles = 0;
  std::atomic<size_t> m_nTotalDuplicatedParticles = 0;
};

}  // namespace ActsExamples
