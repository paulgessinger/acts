// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <memory>

namespace ActsPodioEdm {
class Surface;
}

namespace Acts {

std::shared_ptr<const Surface> convertSurfaceFromPodio(
    const ActsPodioEdm::Surface& surface);

ActsPodioEdm::Surface convertSurfaceToPodio(const Acts::GeometryContext& gctx,
                                            const Acts::Surface& surface);

}  // namespace Acts
