// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Definitions.hpp"

template <typename value_t, size_t DIM>
Acts::Ray<value_t, DIM>::Ray(const vertex_type& origin, const vertex_type& dir)
  : m_origin(origin), m_dir(dir.normalized()), m_idir(1 / m_dir.array())
{
}
template <typename value_t, size_t DIM>
std::ostream&
Acts::Ray<value_t, DIM>::dump(std::ostream& os) const
{
  os << "Ray(";
  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_origin[i];
  }
  os << " -> ";
  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_dir[i];
  }
  os << ")";

  // os << m_idir;

  return os;
}

template <typename value_t, size_t  DIM>
template <typename helper_t, size_t D, std::enable_if_t<D == 3, int>>
void
Acts::Ray<value_t, DIM>::draw(helper_t& helper, value_type far_distance) const
{
  static_assert(DIM == 3, "OBJ is only supported in 3D");
  static_assert(std::is_same<typename helper_t::value_type, value_type>::value,
                "not the same value type");

  helper.line(m_origin, m_origin + m_dir * far_distance);
}

template <typename U, size_t V>
std::ostream&
operator<<(std::ostream& os, const Acts::Ray<U, V>& ray)
{
  ray.dump(os);
  return os;
}
