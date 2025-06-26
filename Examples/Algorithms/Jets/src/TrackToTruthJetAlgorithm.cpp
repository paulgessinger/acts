// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TrackToTruthJetAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ostream>
#include <stdexcept>

namespace ActsExamples {

TrackToTruthJetAlgorithm::TrackToTruthJetAlgorithm(const Config& cfg,
                                                   Acts::Logging::Level lvl)
    : IAlgorithm("TrackToTruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Input tracks are not configured");
  }
  if (m_cfg.inputJets.empty()) {
    throw std::invalid_argument("Input tracks are not configured");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputJets.initialize(m_cfg.inputJets);
  m_outputTrackJets.initialize(m_cfg.outputTrackJets);
}

ProcessCode ActsExamples::TrackToTruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  ACTS_DEBUG("Executing track to truth jet matching algorithm");

  const auto& tracks = m_inputTracks(ctx);
  const auto& truthJets = m_inputJets(ctx);
  // Take a copy that we will modify
  TrackJetContainer jets = truthJets;

  ACTS_DEBUG("TrackToTruthJetAlg - Number of tracks: " << tracks.size());
  ACTS_DEBUG("TrackToTruthJetAlg - Number of truth jets: " << truthJets.size());
  ACTS_DEBUG("TrackToTruthJetAlg - Number of jets: " << jets.size());

  for (const auto& track : tracks) {
    double minDeltaR = m_cfg.maxDeltaR;

    ACTS_VERBOSE("Track index: " << track.index()
                                 << ", momentum: " << track.momentum().x()
                                 << ", " << track.momentum().y() << ", "
                                 << track.momentum().z());

    TrackJet* matchedJet = nullptr;

    // Loop over the jets to find the closest one
    std::size_t i = 0;
    for (auto& jet : jets) {
      Acts::Vector3 jetMom = jet.getFourMomentum().head<3>();
      double deltaR = Acts::VectorHelpers::deltaR(jetMom, track.momentum());

      if (deltaR < minDeltaR) {
        minDeltaR = deltaR;
        matchedJet = &jet;
      }

      ACTS_DEBUG("Track " << track.index() << " delta R to jet " << i << ": "
                          << deltaR);

      i++;
    }  // loop over jets

    if (matchedJet != nullptr) {
      ACTS_VERBOSE("Adding track " << track.index() << " to jet "
                                   << matchedJet->getFourMomentum());
      matchedJet->addTrack(track.index());
    }

  }  // loop over tracks

  // Write the matched jets to the output
  m_outputTrackJets(ctx, std::move(jets));

  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TrackToTruthJetAlgorithm::finalize() {
  ACTS_INFO("Finalizing track to truth jet matching algorithm");
  return ProcessCode::SUCCESS;
}

};  // namespace ActsExamples
