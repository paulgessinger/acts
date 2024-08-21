// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace Acts {

class RegularSurface;
class GeometryContext;
class TrackingVolume;
class CylinderSurface;
class PlaneSurface;
class DiscSurface;

class PortalLinkBase;

class Portal {
 public:
  Portal(std::shared_ptr<RegularSurface> surface);

  static std::shared_ptr<Portal> fuse(const std::shared_ptr<Portal>& aPortal,
                                      const std::shared_ptr<Portal>& bPortal);

  static std::shared_ptr<Portal> mergeAdjacent(
      const std::shared_ptr<Portal>& aPortal,
      const std::shared_ptr<Portal>& bPortal);

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector3& position,
                                      const Vector3& direction) const;

 private:
  // @TODO: Potentially short circuit the virtual call
  // using VolumeResolver = Delegate<const TrackingVolume*(const Vector3&
  // position)>;

  std::shared_ptr<RegularSurface> m_surface;

  std::unique_ptr<PortalLinkBase> m_alongNormal;
  std::unique_ptr<PortalLinkBase> m_oppositeNormal;
};

template <typename S>
concept PortalSurfaceConcept =
    std::is_same_v<S, CylinderSurface> || std::is_same_v<S, DiscSurface> ||
    std::is_same_v<S, PlaneSurface>;

class PortalLinkBase {
 public:
  PortalLinkBase(std::shared_ptr<RegularSurface> surface)
      : m_surface(std::move(surface)) {}

  virtual ~PortalLinkBase() = default;

  // @TODO: Does this need boundary tolerance?
  virtual const TrackingVolume* resolveVolume(
      const GeometryContext& gctx, const Vector2& position) const = 0;

  static std::unique_ptr<PortalLinkBase> merge(
      const std::shared_ptr<PortalLinkBase>& a,
      const std::shared_ptr<PortalLinkBase>& b, BinningValue direction,
      const Logger& logger = getDummyLogger());

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link);

  const RegularSurface& surface() const { return *m_surface; }

 protected:
  static void checkMergePreconditions(const PortalLinkBase& a,
                                      const PortalLinkBase& b,
                                      BinningValue direction);

  std::shared_ptr<RegularSurface> m_surface;
};

}  // namespace Acts