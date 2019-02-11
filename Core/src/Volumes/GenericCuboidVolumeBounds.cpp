// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Volumes/GenericCuboidVolumeBounds.hpp"

#include <array>
#include <ostream>

#include "Acts/Utilities/Definitions.hpp"

Acts::GenericCuboidVolumeBounds::GenericCuboidVolumeBounds(
    const std::array<Acts::Vector3D, 8>& vertices)
  : m_vertices(vertices)
{
  // calculate approximate center of gravity first, so we can make sure
  // the normals point inwards
  Vector3D cog(0, 0, 0);

  for (size_t i = 0; i < 8; i++) {
    cog += m_vertices[i];
  }

  cog *= 0.125;  // 1/8.

  //std::cout << "cog: " << cog.transpose() << std::endl;

  size_t idx = 0;

  auto handle_face
      = [&](const auto& a, const auto& b, const auto& c, const auto& d) {
          // we assume a b c d to be counter clockwise
          const Vector3D ab = b - a, ac = c - a;
          Vector3D       normal = ab.cross(ac).normalized();

          //std::cout << "normal: " << normal.transpose() << std::endl;
          //std::cout << (cog - d).dot(normal) << std::endl;
          if ((cog - d).dot(normal) < 0) {
            // normal points outwards, flip normal
            normal *= -1.;
          }

          // get rid of -0 values if present
          normal += Vector3D::Zero();

          m_normals[idx] = normal;
          idx++;

          //std::cout << "a: " << a.transpose() << std::endl;
          //std::cout << "b: " << b.transpose() << std::endl;
          //std::cout << "c: " << c.transpose() << std::endl;
          //std::cout << "d: " << d.transpose() << std::endl;
          //std::cout << "normal: " << normal.transpose() << std::endl;
          //std::cout << std::endl;
        };

  // handle faces
  handle_face(m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
  handle_face(m_vertices[4], m_vertices[5], m_vertices[6], m_vertices[7]);
  handle_face(m_vertices[0], m_vertices[3], m_vertices[7], m_vertices[4]);
  handle_face(m_vertices[1], m_vertices[2], m_vertices[6], m_vertices[5]);
  handle_face(m_vertices[2], m_vertices[3], m_vertices[7], m_vertices[6]);
  handle_face(m_vertices[1], m_vertices[0], m_vertices[4], m_vertices[5]);
}

Acts::VolumeBounds*
Acts::GenericCuboidVolumeBounds::clone() const
{
  return new GenericCuboidVolumeBounds(*this);
}

bool
Acts::GenericCuboidVolumeBounds::inside(const Acts::Vector3D& gpos,
                                        double                tol) const
{
  constexpr std::array<size_t, 6> vtxs = {0, 4, 0, 1, 2, 1};
  // needs to be on same side, get ref
  bool ref = (gpos-m_vertices[vtxs[0]]).dot(m_normals[0]) > 0;
  for(size_t i=1;i<6;i++) {
    if(((gpos-m_vertices[vtxs[i]]).dot(m_normals[i]) > 0) != ref) {
      return false;
    }
  }
  return true;
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::GenericCuboidVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<const Acts::Transform3D> transform) const
{
  return {};
}

std::ostream&
Acts::GenericCuboidVolumeBounds::dump(std::ostream& sl) const
{
  sl << "Acts::GenericCuboidVolumeBounds: vertices (x, y, z) =\n";
  for (size_t i = 0; i < 8; i++) {
    if (i > 0) {
      sl << ",\n";
    }
    sl << "[" << m_vertices[i].transpose() << "]";
  }
  return sl;
}
