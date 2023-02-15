// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackSelector.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

ActsExamples::TrackSelector::TrackSelector(const Config& config,
                                           Acts::Logging::Level level)
    : BareAlgorithm("TrackSelector", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Input track collection is empty");
  }

  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Output track collection is empty");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackSelector::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // helper functions to select tracks
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidTrack = [&](const auto& trk) {
    const auto theta = trk.theta();
    const auto eta = -std::log(std::tan(theta / 2));
    // define charge selection
    const bool validNeutral = (trk.charge() == 0) and not m_cfg.removeNeutral;
    const bool validCharged = (trk.charge() != 0) and not m_cfg.removeCharged;
    const bool validCharge = validNeutral or validCharged;
    return validCharge and
           within(trk.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
           within(eta, m_cfg.etaMin, m_cfg.etaMax) and
           within(trk.phi(), m_cfg.phiMin, m_cfg.phiMax) and
           within(trk.loc0(), m_cfg.loc0Min, m_cfg.loc0Max) and
           within(trk.loc1(), m_cfg.loc1Min, m_cfg.loc1Max) and
           within(trk.time(), m_cfg.timeMin, m_cfg.timeMax);
  };

  ACTS_VERBOSE("Reading tracks from: " << m_cfg.inputTracks);

  const auto& inputTracks =
      ctx.eventStore.get<ConstTrackContainer>(m_cfg.inputTracks);

  std::shared_ptr<Acts::ConstVectorMultiTrajectory> trackStateContainer =
      inputTracks.trackStateContainerHolder();

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>(
      inputTracks.container());  // mutable copy of the immutable
                                 // source track container
  auto tempTrackStateContainer =
      std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer filteredTracks{trackContainer, tempTrackStateContainer};

  ACTS_VERBOSE(
      "Track container size before filtering: " << filteredTracks.size());

  for (Acts::MultiTrajectoryTraits::IndexType itrack = 0;
       itrack < filteredTracks.size();) {
    auto track = filteredTracks.getTrack(itrack);
    if (isValidTrack(track)) {
      // Track is valid, go to next index
      ACTS_VERBOSE(" - Keeping track #" << itrack);
      itrack++;
      continue;
    }
    ACTS_VERBOSE(" - Removing track #" << itrack);
    filteredTracks.removeTrack(itrack);
    // Do not increment track index
  }

  ACTS_VERBOSE(
      "Track container size after filtering: " << filteredTracks.size());

  // @TODO: Make track container with original MTJ + write to whiteboard

  // for (const auto& srcTrack : inputTracks) {
  // auto dstTrack = tempTracks.getTrack(tempTracks.addTrack());
  // }

  // TrajectoriesContainer outputTrajectories;
  // outputTrajectories.reserve(inputTrajectories.size());

  // std::size_t inputCount = 0;
  // std::size_t outputCount = 0;
  // for (const auto& trajectories : inputTrajectories) {
  // std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;
  // Trajectories::IndexedParameters parameters;

  // for (auto tip : trajectories.tips()) {
  // if (!trajectories.hasTrackParameters(tip)) {
  // continue;
  // }
  // ++inputCount;
  // if (!isValidTrack(trajectories.trackParameters(tip))) {
  // continue;
  // }
  // tips.push_back(tip);
  // parameters.emplace(tip, trajectories.trackParameters(tip));
  // ++outputCount;
  // }

  // outputTrajectories.emplace_back(trajectories.multiTrajectory(), tips,
  // parameters);
  // }

  // ACTS_DEBUG("event " << ctx.eventNumber << " selected " << outputCount
  // << " from " << inputCount << " tracks in trajectories");

  // ctx.eventStore.add(m_cfg.outputTrajectories,
  // std::move(outputTrajectories));
  // }

  return ProcessCode::SUCCESS;
}
