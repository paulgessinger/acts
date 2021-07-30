// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::TrapezoidBounds::~TrapezoidBounds() = default;

Acts::SurfaceBounds::BoundsType Acts::TrapezoidBounds::type() const {
  return SurfaceBounds::eTrapezoid;
}

bool Acts::TrapezoidBounds::inside(const Acts::Vector2& lposition,
                                   const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, m_vertices);
}

void Acts::TrapezoidBounds::vertices(std::vector<Acts::Vector2>& result,
                                     unsigned int /*lseg*/) const {
  double minhx = get(TrapezoidBounds::eHalfLengthXnegY);
  double maxhx = get(TrapezoidBounds::eHalfLengthXposY);
  double hy = get(TrapezoidBounds::eHalfLengthY);
  result = {{-minhx, -hy}, {minhx, -hy}, {maxhx, hy}, {-maxhx, hy}};
}

const Acts::RectangleBounds& Acts::TrapezoidBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& Acts::TrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (halfXnegY, halfXposY, halfY) = "
     << "(" << get(eHalfLengthXnegY) << ", " << get(eHalfLengthXposY) << ", "
     << get(eHalfLengthY) << ")";
  sl << std::setprecision(-1);
  return sl;
}

void Acts::TrapezoidBounds::fillVertices() {
  std::vector<Vector2> vtxs;
  vertices(vtxs);
  for (size_t i = 0; i < 4; i++) {
    m_vertices[i] = vtxs[i];
  }
}