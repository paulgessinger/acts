// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <ostream>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

namespace Acts {

class GenericCuboidVolumeBounds : public VolumeBounds
{

public:
  GenericCuboidVolumeBounds() = delete;

  /// Constructor from a set of vertices
  /// The ordering is considered to be:
  /// - the first 4 vertices are the "top" face
  /// - the second 4 vertices are the "bottom" face
  /// - both faces are given in counter clock wise order
  GenericCuboidVolumeBounds(const std::array<Acts::Vector3D, 8>& vertices);

  virtual ~GenericCuboidVolumeBounds() = default;

  ///  clone() method to make deep copy in Volume copy constructor and for
  /// assigment operator  of the Surface class.
  virtual VolumeBounds*
  clone() const;

  /// Checking if position given in volume frame is inside
  ///
  /// @param gpos is the global position to be checked
  /// @param tol is the tolerance applied for the inside check
  ///
  /// @return boolean indicating if the position is inside
  virtual bool
  inside(const Vector3D& gpos, double tol = 0.) const;

  /// Method to decompose the Bounds into Surfaces
  /// the Volume can turn them into BoundarySurfaces
  ///
  /// @param transform is the 3D transform to be applied to the boundary
  /// surfaces to position them in 3D space
  /// @note this is factory method
  ///
  /// @return a vector of surfaces bounding this volume
  virtual std::vector<std::shared_ptr<const Surface>>
  decomposeToSurfaces(const Transform3D* transform) const;

  AABB3F<Volume>
  boundingBox(const Transform3D* trf      = nullptr,
              const Vector3F&    envelope = {0, 0, 0}) const final;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the output stream to be dumped into
  virtual std::ostream&
  dump(std::ostream& sl) const;

  template <typename helper_t>
  void
  draw(helper_t&          helper,
       const Transform3D& transform = Transform3D::Identity()) const;

private:
  std::array<Vector3D, 8> m_vertices;
  // which vertices go with which normal (this is an implementation detail)
  // static constexpr std::array<size_t, 6> s_planeVertexIndices = {0, 4, 0, 1,
  // 2, 1};
  std::array<Vector3D, 6> m_normals;
};
}  // namespace Acts

template <typename helper_t>
void
Acts::GenericCuboidVolumeBounds::draw(helper_t&          helper,
                                      const Transform3D& transform) const
{
  auto draw_face
      = [&](const auto& a, const auto& b, const auto& c, const auto& d) {
          helper.face(std::vector<Vector3D>(
              {transform * a, transform * b, transform * c, transform * d}));
        };

  draw_face(m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
  draw_face(m_vertices[4], m_vertices[5], m_vertices[6], m_vertices[7]);
  draw_face(m_vertices[0], m_vertices[3], m_vertices[7], m_vertices[4]);
  draw_face(m_vertices[1], m_vertices[2], m_vertices[6], m_vertices[5]);
  draw_face(m_vertices[2], m_vertices[3], m_vertices[7], m_vertices[6]);
  draw_face(m_vertices[1], m_vertices[0], m_vertices[4], m_vertices[5]);
}
