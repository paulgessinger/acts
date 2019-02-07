// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

template <int N>
class ConvexPolygonBounds : public PlanarBounds
{
  static constexpr size_t num_vertices = N;
  using vertex_array                   = std::array<Vector2D, num_vertices>;

public:
  static_assert(N >= 3, "ConvexPolygonBounds needs at least 3 sides.");

  /// Trapezoid bounds default constructor is forbidden
  ConvexPolygonBounds() = delete;

  ConvexPolygonBounds(const std::vector<Vector2D>& vertices);

  ConvexPolygonBounds(const vertex_array& vertices);

  ~ConvexPolygonBounds() override = default;

  ConvexPolygonBounds<N>*
  clone() const final;

  BoundsType
  type() const final;

  bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const final;

  double
  distanceToBoundary(const Vector2D& lpos) const final;

  /// Return the vertices - or, the points of the extremas
  std::vector<Vector2D>
  vertices() const final;

  std::vector<TDD_real_t>
  valueStore() const final;

  // Bounding box representation
  const RectangleBounds&
  boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream&
  dump(std::ostream& sl) const final;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  virtual variant_data
  toVariantData() const final;

private:
  vertex_array    m_vertices;
  RectangleBounds m_boundingBox;

  template <typename coll_t>
  static RectangleBounds
  makeBoundingBox(const coll_t& vertices);
};

}  // namespace

#include "Acts/Surfaces/ConvexPolygonBounds.ipp"
