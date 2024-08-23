// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>

namespace Acts {

Portal::Portal(const GeometryContext& gctx, Direction direction,
               std::unique_ptr<PortalLinkBase> link) {
  setLink(gctx, direction, std::move(link));
  assert(m_surface != nullptr);
}

Portal::Portal(const GeometryContext& gctx, Direction direction,
               std::shared_ptr<RegularSurface> surface, TrackingVolume& volume)
    : Portal(gctx, direction,
             std::make_unique<TrivialPortalLink>(std::move(surface), volume)) {}

Portal::Portal(const GeometryContext& gctx,
               std::unique_ptr<PortalLinkBase> alongNormal,
               std::unique_ptr<PortalLinkBase> oppositeNormal) {
  if (alongNormal == nullptr && oppositeNormal == nullptr) {
    throw std::invalid_argument("At least one link must be provided");
  }

  if (alongNormal != nullptr) {
    setLink(gctx, Direction::AlongNormal, std::move(alongNormal));
  }
  if (oppositeNormal != nullptr) {
    setLink(gctx, Direction::OppositeNormal, std::move(oppositeNormal));
  }
}

Portal::Portal(const GeometryContext& gctx, Config&& config) {
  if (!config.alongNormal.m_surface && !config.oppositeNormal.m_surface) {
    throw std::invalid_argument("At least one link must be provided");
  }

  if (config.alongNormal.m_surface) {
    setLink(gctx, Direction::AlongNormal,
            std::make_unique<TrivialPortalLink>(
                std::move(config.alongNormal.m_surface),
                *config.alongNormal.m_volume));
  }
  if (config.oppositeNormal.m_surface) {
    setLink(gctx, Direction::OppositeNormal,
            std::make_unique<TrivialPortalLink>(
                std::move(config.oppositeNormal.m_surface),
                *config.oppositeNormal.m_volume));
  }
}

void Portal::setLink(const GeometryContext& gctx, Direction direction,
                     std::unique_ptr<PortalLinkBase> link) {
  assert(link != nullptr);

  auto& target =
      direction == Direction::AlongNormal ? m_alongNormal : m_oppositeNormal;
  auto& other =
      direction == Direction::AlongNormal ? m_oppositeNormal : m_alongNormal;

  // check if surfaces are identical
  if (m_surface != nullptr &&
      !isSameSurface(gctx, link->surface(), *m_surface)) {
    throw PortalFusingException();
  }

  // check if they both have material but are not the same surface
  if (m_surface != nullptr && (m_surface.get() != &link->surface()) &&
      link->surface().surfaceMaterial() != nullptr &&
      m_surface->surfaceMaterial() != nullptr) {
    throw PortalFusingException();
  }

  target = std::move(link);

  if (other == nullptr) {
    // We don't have an existing surface, take the one we just got
    m_surface = target->surfacePtr();
    return;
  }

  if (target->surface().surfaceMaterial() != nullptr) {
    // new link has material: assign that to existing link
    m_surface = target->surfacePtr();
    other->setSurface(m_surface);
  } else {
    // none have material, or the existing surface had material: assign the
    // existing surface by convention
    target->setSurface(m_surface);
  }
}

void Portal::setLink(const GeometryContext& gctx, Direction direction,
                     std::shared_ptr<RegularSurface> surface,
                     TrackingVolume& volume) {
  setLink(gctx, direction,
          std::make_unique<TrivialPortalLink>(std::move(surface), volume));
}

const PortalLinkBase* Portal::getLink(Direction direction) const {
  if (direction == Direction::AlongNormal) {
    return m_alongNormal.get();
  } else {
    return m_oppositeNormal.get();
  }
}

Result<const TrackingVolume*> Portal::resolveVolume(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(m_surface != nullptr);
  const Vector3 normal = m_surface->normal(gctx, position);
  Direction side = Direction::fromScalar(normal.dot(direction));

  const PortalLinkBase* link = side == Direction::AlongNormal
                                   ? m_alongNormal.get()
                                   : m_oppositeNormal.get();

  if (link == nullptr) {
    // no link is attached in this direction => this is the end of the world as
    // we know it. (i feel fine)
    return nullptr;
  } else {
    auto res = link->resolveVolume(gctx, position);
    if (!res.ok()) {
      return res.error();
    }
    return *res;
  }
}

bool Portal::isValid() const {
  return m_alongNormal != nullptr || m_oppositeNormal != nullptr;
}

const RegularSurface& Portal::surface() const {
  assert(m_surface != nullptr);
  return *m_surface;
}

Portal Portal::merge(const GeometryContext& gctx, Portal& aPortal,
                     Portal& bPortal, BinningValue direction,
                     const Logger& logger) {
  ACTS_DEBUG("Merging to portals along " << direction);

  if (&aPortal == &bPortal) {
    ACTS_ERROR("Cannot merge a portal with itself");
    throw PortalMergingException{};
  }

  if (aPortal.m_surface->surfaceMaterial() != nullptr ||
      bPortal.m_surface->surfaceMaterial() != nullptr) {
    ACTS_ERROR("Cannot merge portals with material");
    throw PortalMergingException{};
  }

  std::unique_ptr<PortalLinkBase> mergedAlongNormal = nullptr;
  std::unique_ptr<PortalLinkBase> mergedOppositeNormal = nullptr;

  bool aHasAlongNormal = aPortal.m_alongNormal != nullptr;
  bool aHasOppositeNormal = aPortal.m_oppositeNormal != nullptr;
  bool bHasAlongNormal = bPortal.m_alongNormal != nullptr;
  bool bHasOppositeNormal = bPortal.m_oppositeNormal != nullptr;

  if (aHasAlongNormal != bHasAlongNormal ||
      aHasOppositeNormal != bHasOppositeNormal) {
    ACTS_ERROR("Portals do not have the same links attached");
    throw PortalMergingException();
  }

  if (aPortal.m_alongNormal != nullptr) {
    if (bPortal.m_alongNormal == nullptr) {
      ACTS_ERROR(
          "Portal A has link along normal, while b does not. This is not "
          "supported");
      throw PortalMergingException();
    }

    ACTS_VERBOSE("Portals have links along normal, merging");
    mergedAlongNormal = PortalLinkBase::merge(std::move(aPortal.m_alongNormal),
                                              std::move(bPortal.m_alongNormal),
                                              direction, logger);
  }

  if (aPortal.m_oppositeNormal != nullptr) {
    if (bPortal.m_oppositeNormal == nullptr) {
      ACTS_ERROR(
          "Portal A has link opposite normal, while b does not. This is not "
          "supported");
      throw PortalMergingException();
    }

    ACTS_VERBOSE("Portals have links opposite normal, merging");
    mergedOppositeNormal = PortalLinkBase::merge(
        std::move(aPortal.m_oppositeNormal),
        std::move(bPortal.m_oppositeNormal), direction, logger);
  }

  aPortal.m_surface.reset();
  bPortal.m_surface.reset();
  return Portal{gctx, std::move(mergedAlongNormal),
                std::move(mergedOppositeNormal)};
}

Portal Portal::fuse(const GeometryContext& gctx, Portal& aPortal,
                    Portal& bPortal, const Logger& logger) {
  ACTS_DEBUG("Fusing two portals");
  if (&aPortal == &bPortal) {
    ACTS_ERROR("Cannot merge a portal with itself");
    throw PortalMergingException{};
  }

  bool aHasAlongNormal = aPortal.m_alongNormal != nullptr;
  bool aHasOppositeNormal = aPortal.m_oppositeNormal != nullptr;
  bool bHasAlongNormal = bPortal.m_alongNormal != nullptr;
  bool bHasOppositeNormal = bPortal.m_oppositeNormal != nullptr;

  if (aPortal.m_surface == nullptr || bPortal.m_surface == nullptr) {
    ACTS_ERROR("Portals have no surface");
    throw PortalFusingException();
  }

  if (aPortal.m_surface->associatedDetectorElement() != nullptr ||
      bPortal.m_surface->associatedDetectorElement() != nullptr) {
    ACTS_ERROR("Cannot fuse portals with detector elements");
    throw PortalFusingException();
  }

  if (!isSameSurface(gctx, *aPortal.m_surface, *bPortal.m_surface)) {
    ACTS_ERROR("Portals have different surfaces");
    throw PortalFusingException();
  }

  if (aPortal.m_surface->surfaceMaterial() != nullptr &&
      bPortal.m_surface->surfaceMaterial() != nullptr) {
    ACTS_ERROR("Cannot fuse portals if both have material");
    throw PortalFusingException();
  }

  if (aHasAlongNormal == bHasAlongNormal ||
      aHasOppositeNormal == bHasOppositeNormal) {
    ACTS_ERROR("Portals have the same links attached");
    throw PortalFusingException();
  }

  aPortal.m_surface.reset();
  bPortal.m_surface.reset();
  if (aHasAlongNormal) {
    ACTS_VERBOSE("Taking along normal from lhs, opposite normal from rhs");
    return Portal{gctx, std::move(aPortal.m_alongNormal),
                  std::move(bPortal.m_oppositeNormal)};
  } else {
    ACTS_VERBOSE("Taking along normal from rhs, opposite normal from lhs");
    return Portal{gctx, std::move(bPortal.m_alongNormal),
                  std::move(aPortal.m_oppositeNormal)};
  }
}

bool Portal::isSameSurface(const GeometryContext& gctx, const Surface& a,
                           const Surface& b) {
  if (&a == &b) {
    return true;
  }

  if (a.type() != b.type()) {
    return false;
  }

  if (a.bounds() != b.bounds()) {
    return false;
  }

  if (!a.transform(gctx).isApprox(b.transform(gctx), 1e-9)) {
    return false;
  }

  return true;
};

}  // namespace Acts
