// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/TracksToTrajectories.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

namespace ActsExamples {

ProcessCode TracksToTrajectories::execute(const AlgorithmContext& ctx) const {
  const auto& tracks = ctx.eventStore.get<TrackContainer>(m_cfg.inputTracks);

  // Prepare the output data with MultiTrajectory
  TrajectoriesContainer trajectories;
  // trajectories.reserve(initialParameters.size());

  Trajectories::IndexedParameters parameters;
  parameters.reserve(tracks.size());
  std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;
  tips.reserve(tracks.size());

  for (const auto& track : tracks) {
    tips.push_back(track.tipIndex());
    parameters.emplace(
        std::pair{track.tipIndex(),
                  TrackParameters{track.referenceSurface().getSharedPtr(),
                                  track.parameters(), track.covariance()}});
  }

  trajectories.emplace_back(tracks.trackStateContainer(), std::move(tips),
                            std::move(parameters));

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
