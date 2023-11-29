// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <boost/core/span.hpp>

namespace Acts::Experimental {

class DetectorVolume;
struct NavigationState;

class IDetectorVolumeUpdater {
 public:
  virtual ~IDetectorVolumeUpdater() = default;
  virtual boost::span<const DetectorVolume* const> volumes() const {
    return {};
  };
};

/// Declare a Detector Volume finding or switching delegate
///
/// @param gctx is the current geometry context
/// @param nState [in, out] is the navigation state to be updated
///
/// @return the new DetectorVolume into which one changes at this switch
using DetectorVolumeUpdater =
    OwningDelegate<void(const GeometryContext& gctx, NavigationState& nState),
                   IDetectorVolumeUpdater>;

}  // namespace Acts::Experimental
