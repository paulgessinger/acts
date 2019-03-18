// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

#include <ostream>

namespace Acts {

template <typename value_t, size_t DIM, size_t SIDES>
class Frustum
{
  using translation_t = Eigen::Translation<value_t, DIM>;

public:
  using transform_type = Eigen::Transform<value_t, DIM, Eigen::Affine>;
  using value_type        = value_t;
  using vertex_type       = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;

  static const size_t dim   = DIM;
  static const size_t sides = SIDES;

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  Frustum(const vertex_type& origin,
          const vertex_type& dir,
          value_type         opening_angle);

  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  Frustum(const vertex_type& origin,
          const vertex_type& dir,
          value_type         opening_angle);

  template <typename helper_t,
            size_t D = DIM,
            std::enable_if_t<D == 3, int> = 0>
  void
  draw(helper_t& helper, value_type far_distance = 10) const;

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  std::ostream&
  svg(std::ostream& os,
      value_type    w,
      value_type    h,
      value_type    far_distance = 1,
      value_type    unit         = 20.) const;

  const vertex_type&
  origin() const
  {
    return m_origin;
  }

  const vertex_type&
  dir() const
  {
    return m_normals[0];
  }

  const std::array<vertex_type, SIDES + 1>&
  normals() const
  {
    return m_normals;
  }

private:
  vertex_type m_origin;
  // need one more for direction we're facing
  std::array<vertex_type, SIDES + 1> m_normals;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace Acts

#include "Acts/Utilities/Frustum.ipp"
