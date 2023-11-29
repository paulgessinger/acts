// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Portal.hpp"

namespace Acts::Experimental {

void updateCandidates(const GeometryContext& gctx, NavigationState& nState) {
  const auto& position = nState.position;
  const auto& direction = nState.direction;
  auto& nCandidates = nState.surfaceCandidates;

  for (auto& c : nCandidates) {
    // Get the surface representation: either native surfcae of portal
    const Surface& sRep =
        c.surface != nullptr ? *c.surface : c.portal->surface();

    // Get the intersection @todo make a templated intersector
    // TODO surface tolerance
    auto sIntersection = sRep.intersect(gctx, position, direction,
                                        c.boundaryCheck, s_onSurfaceTolerance);
    c.objectIntersection = sIntersection[c.objectIntersection.index()];
  }
}

}  // namespace Acts::Experimental
