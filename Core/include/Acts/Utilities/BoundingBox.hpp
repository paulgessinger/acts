// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Utilities/Visualization.hpp"

namespace Acts {

template <typename entity_t, typename value_t, size_t DIM>
class AxisAlignedBoundingBox
{
private:
  using self_t = AxisAlignedBoundingBox<entity_t, value_t, DIM>;

  // strong type helper, not public
  template <typename T, typename P>
  class NamedType
  {
  public:
    explicit NamedType(const T& value) : m_value(value) {}
    explicit NamedType(T&& value) : m_value(std::move(value)) {}
    T&
    get()
    {
      return m_value;
    }
    const T&
    get() const
    {
      return m_value;
    }

  private:
    T m_value;
  };

  struct SizeParameter
  {
  };

public:
  using vertex_type       = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;
  using entity_type       = entity_t;
  using value_type        = value_t;

  /// Strong type to select the correct constructor
  using Size = NamedType<vertex_type, struct SizeParameter>;

  static const size_t dim = DIM;

  // if we construct this with an entity, the entity can not be null
  AxisAlignedBoundingBox(const entity_t&    entity,
                         const vertex_type& vmin,
                         const vertex_type& vmax);

  // if we construct this with an entity, the entity can not be null
  AxisAlignedBoundingBox(const entity_t&    entity,
                         const vertex_type& center,
                         const Size&        size);

  AxisAlignedBoundingBox(const std::vector<self_t*>& boxes,
                         vertex_array_type           envelope
                         = vertex_array_type::Zero());

  static std::pair<vertex_type, vertex_type>
  wrap(const std::vector<const self_t*>& boxes,
       vertex_array_type                 envelope = vertex_array_type::Zero());

  static std::pair<vertex_type, vertex_type>
  wrap(const std::vector<self_t*>& boxes,
       vertex_array_type           envelope = vertex_array_type::Zero());

  static std::pair<vertex_type, vertex_type>
  wrap(const std::vector<self_t>& boxes,
       vertex_array_type          envelope = vertex_array_type::Zero());

  bool
  intersect(const vertex_type& point) const;

  /// @brief Implements the slab method for Ray/AABB intersections.
  ///
  /// See https://tavianator.com/fast-branchless-raybounding-box-intersections/,
  /// https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/,
  /// https://medium.com/@bromanz/another-view-on-the-classic-ray-aabb-intersection-algorithm-for-bvh-traversal-41125138b525
  ///
  /// @note This implementation may treat parallel rays on any of the slabs
  ///       as **outside** due to how @c NaNs are handled by Eigen.
  ///       See http://eigen.tuxfamily.org/bz/show_bug.cgi?id=564
  bool
  intersect(const Ray<value_type, DIM>& ray) const;

  template <size_t sides>
  bool
  intersect(const Frustum<value_type, DIM, sides>& fr) const;

  void
  setSkip(self_t* skip);

  const self_t*
  getLeftChild() const;

  const self_t*
  getSkip() const;

  bool
  hasEntity() const;

  const entity_t*
  entity() const;

  const vertex_type&
  center() const;

  const vertex_type&
  min() const;

  const vertex_type&
  max() const;

  std::ostream&
  dump(std::ostream& os) const;

  template <typename helper_t,
            size_t D = DIM,
            std::enable_if_t<D == 3, int> = 0>
  void
  draw(helper_t& helper, std::array<int, 3> color = {120, 120, 120}) const;

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  std::ofstream&
  svg(std::ofstream& os,
      value_type     w,
      value_type     h,
      value_type     unit      = 10,
      std::string    label     = "",
      std::string    fillcolor = "grey") const;

private:
  const entity_t*   m_entity;
  vertex_type       m_vmin;
  vertex_type       m_vmax;
  vertex_type       m_center;
  vertex_array_type m_width;
  vertex_array_type m_iwidth;

  self_t* m_left_child{nullptr};
  self_t* m_right_child{nullptr};
  self_t* m_skip{nullptr};

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename box_t>
box_t*
make_octree(std::vector<std::unique_ptr<box_t>>& store,
            const std::vector<box_t*>&           prims,
            size_t                               max_depth = 1,
            typename box_t::value_type           envelope1 = 0);

template <typename T, typename U, size_t V>
std::ostream&
operator<<(std::ostream& os, const AxisAlignedBoundingBox<T, U, V>& box);

// default type and dim
template <typename T>
using AABB = AxisAlignedBoundingBox<T, float, 3>;

}  // namespace Acts

#include "Acts/Utilities/BoundingBox.ipp"
