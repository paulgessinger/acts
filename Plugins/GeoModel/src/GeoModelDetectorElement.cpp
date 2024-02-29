// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"

#include "Acts/Surfaces/Surface.hpp"

#include <utility>

Acts::GeoModelDetectorElement::GeoModelDetectorElement(
    std::shared_ptr<Acts::Surface> surface, const GeoVPhysVol& geoPhysVol,
    const Acts::Transform3& toGlobal, Acts::ActsScalar thickness)
    : m_surface(std::move(surface)),
      m_geoPhysVol(&geoPhysVol),
      m_toGlobal(toGlobal),
      m_thickness(thickness) {}

const Acts::Transform3& Acts::GeoModelDetectorElement::transform(
    const GeometryContext& /*gctx*/) const {
  return m_toGlobal;
}

const Acts::Surface& Acts::GeoModelDetectorElement::surface() const {
  return *m_surface;
}

Acts::Surface& Acts::GeoModelDetectorElement::surface() {
  return *m_surface;
}

Acts::ActsScalar Acts::GeoModelDetectorElement::thickness() const {
  return m_thickness;
}

const GeoVPhysVol& Acts::GeoModelDetectorElement::geoVPhysicalVolume() const {
  return *m_geoPhysVol;
}
