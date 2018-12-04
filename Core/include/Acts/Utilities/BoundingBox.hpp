// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <limits>
#include <tuple>
#include <algorithm>

namespace Acts {

template <typename value_t, size_t DIM>
class Ray 
{
public:
  using value_type = value_t;
  using vertex_type = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;

  Ray(const vertex_type& origin, const vertex_type& dir)
    : m_origin(origin),
      m_dir(dir.normalized()),
      m_idir(1/m_dir.array())
  {}

  const vertex_type& origin() const { return m_origin; }
  const vertex_type& dir() const { return m_dir; }
  const vertex_array_type& idir() const { return m_idir; }

  std::ostream& dump(std::ostream& os) const
  {
    os << "Ray(";
    for(size_t i=0;i<DIM;i++) {
      if(i>0) {
        os << ", ";
      }
      os << m_origin[i];
    }
    os << " -> ";
    for(size_t i=0;i<DIM;i++) {
      if(i>0) {
        os << ", ";
      }
      os << m_dir[i];
    }
    os << ")";

    os << m_idir;

    return os;
  }

private:
  vertex_type m_origin;
  vertex_type m_dir;
  vertex_array_type m_idir;
};

template <typename U, size_t V>
std::ostream& operator<<(std::ostream& os, const Ray<U, V>& ray) 
{
  ray.dump(os);
  return os;
}

using Ray3F = Ray<float, 3>;

template <typename entity_t, typename value_t, size_t DIM>
class AxisAlignedBoundingBox
{
private:
  using self_t = AxisAlignedBoundingBox<entity_t, value_t, DIM>;
public:
  using vertex_type = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;
  using entity_type = entity_t;
  using value_type = value_t;

  // if we construct this with an entity, the entity can not be null
  AxisAlignedBoundingBox(const entity_t& entity, const vertex_type& vmin, const vertex_type& vmax)
    : m_entity(&entity),
      m_vertices({vmin, vmax}),
      m_center((vmin + vmax)/2.),
      m_width(vmax - vmin),
      m_iwidth(1/m_width)
  {
  }

  AxisAlignedBoundingBox(const std::vector<self_t*>& boxes,
                         vertex_array_type envelope = vertex_array_type::Zero())
    : m_entity(nullptr)
  {
    assert(boxes.size() > 1);


    for(size_t i=0;i<boxes.size();i++) {
      if(i<boxes.size()-1) {
        // set next on i to i+1
        boxes[i]->setSkip(boxes[i+1]);
      }
      else {
        // make sure last is set to nullptr, this marks end
        //boxes[i]->m_next = nullptr;
        boxes[i]->setSkip(nullptr);
      }
    }

    m_left_child = boxes.front();
    m_right_child = boxes.back();
    m_skip = nullptr;

    std::tie(m_vertices[0], m_vertices[1]) = wrap(boxes, envelope);

    m_center = (m_vertices[0] + m_vertices[1])/2.;
    m_width = m_vertices[1] - m_vertices[0];
    m_iwidth = 1/m_width;
  }

  static 
  std::pair<vertex_type, vertex_type>
  wrap(const std::vector<const self_t*>& boxes,
              vertex_array_type envelope = vertex_array_type::Zero()) 
  {
    assert(boxes.size() > 1);
    // figure out extent of boxes
    // use array for Eigen coefficient wise min/max
    vertex_array_type vmax(
        vertex_array_type::Constant(std::numeric_limits<value_type>::lowest()));
    vertex_array_type vmin(
        vertex_array_type::Constant(std::numeric_limits<value_type>::max()));

    for(size_t i=0;i<boxes.size();i++) {
      vmin = vmin.min(boxes[i]->m_vertices[0].array());
      vmax = vmax.max(boxes[i]->m_vertices[1].array());
    }

    vmax += envelope;
    vmin -= envelope;

    return {vmin, vmax};
  }
  
  static 
  std::pair<vertex_type, vertex_type>
  wrap(const std::vector<self_t*>& boxes,
       vertex_array_type envelope = vertex_array_type::Zero())
  {
    assert(boxes.size() > 1);
    std::vector<const self_t*> box_ptrs;
    box_ptrs.reserve(boxes.size());
    std::transform(boxes.begin(), boxes.end(), std::back_inserter(box_ptrs),
        [](const auto* box) { return box; });
    return wrap(box_ptrs, envelope);
  }
  
  static 
  std::pair<vertex_type, vertex_type>
  wrap(const std::vector<self_t>& boxes,
              vertex_array_type envelope = vertex_array_type::Zero())
  {
    assert(boxes.size() > 1);
    std::vector<const self_t*> box_ptrs;
    box_ptrs.reserve(boxes.size());
    std::transform(boxes.begin(), boxes.end(), std::back_inserter(box_ptrs),
        [](auto& box) { return &box; });
    return wrap(box_ptrs, envelope);
  }

  bool
  intersect(const vertex_type& point) const
  {
    vertex_array_type t = (point - m_vertices[0]).array() * m_iwidth;
    return t.minCoeff() >= 0 && t.maxCoeff() < 1;
  }

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
  intersect(const Ray<value_type, DIM>& ray) const
  {
    //std::cout << "--- " << ray << std::endl;

    const vertex_type& origin = ray.origin();
    const vertex_array_type& idir = ray.idir();
    // this is NaN origin is on box boundary and ray is parallel to
    // that boundary, since 0*inf = NaN. 
    vertex_array_type t0s = (m_vertices[0] - origin).array() * idir;
    vertex_array_type t1s = (m_vertices[1] - origin).array() * idir;
    
    //std::cout << "t0s:\n" << t0s << "\nt1s\n" << t1s << "\n\n";

    // this is non-compliant with IEEE-754-2008, NaN gets propagated through
    // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=564
    // this means that rays parallel to boundaries might not be considered 
    // to intersect.
    vertex_array_type tsmaller = t0s.min(t1s);
    vertex_array_type tbigger = t0s.max(t1s);

    //std::cout << "tsmaller:\n" << tsmaller << "\ntbigger\n" << tbigger << "\n\n";

    value_type tmin = tsmaller.maxCoeff();
    value_type tmax = tbigger.minCoeff();

    //std::cout << "tmin:" << tmin << "\ntmax:" << tmax << std::endl;
    
    return tmin < tmax && tmin > 0.0;
  }

private:

  void setSkip(self_t* skip) {
    // set next on this
    m_skip = skip;
    // find last child and set its next
    if(m_right_child != nullptr) {
      m_right_child->setSkip(skip);
    }
  }

public:

  const self_t* getLeftChild() const
  {
    return m_left_child;
  }

  const self_t* getSkip() const 
  {
    return m_skip;
  }

  bool hasEntity() const 
  {
    return m_entity != nullptr;
  }

  const entity_t* entity() const
  {
    assert(m_entity != nullptr);
    return m_entity;
  }

  const vertex_type& center() const
  {
    return m_center;
    //return (m_vertices[0] + m_vertices[1]) /2.;
  }

  const vertex_type& min() const
  {
    return m_vertices[0];
  }

  const vertex_type& max() const
  {
    return m_vertices[1];
  }

  std::ostream& dump(std::ostream& os) const
  {
    os << "AxisAlignedBoundingBox(center=(";

    for(size_t i=0;i<DIM;i++) {
      if (i>0) {
        os << ", ";
      }
      os << m_center[i];
    }


    os << ") vmin=(";
    for(size_t i=0;i<DIM;i++) {
      if (i>0) {
        os << ", ";
      }
      os << m_vertices[0][i];
    }

    os << ") vmax=(";

    for(size_t i=0;i<DIM;i++) {
      if (i>0) {
        os << ", ";
      }
      os << m_vertices[1][i];
    }

    os << "))";

    return os;
  }

  void obj(std::ostream& os, size_t& vtx_offset) const 
  {
    assert(DIM == 3);
    using face_t = std::array<vertex_type, 4>;
    const vertex_type& vmin = m_vertices[0];
    const vertex_type& vmax = m_vertices[1];

    face_t min_x = {
      vertex_type(vmin.x(), vmin.y(), vmin.z()),
      vertex_type(vmin.x(), vmax.y(), vmin.z()),
      vertex_type(vmin.x(), vmax.y(), vmax.z()),
      vertex_type(vmin.x(), vmin.y(), vmax.z())
    };
    
    face_t max_x = {
      vertex_type(vmax.x(), vmin.y(), vmin.z()),
      vertex_type(vmax.x(), vmax.y(), vmin.z()),
      vertex_type(vmax.x(), vmax.y(), vmax.z()),
      vertex_type(vmax.x(), vmin.y(), vmax.z())
    };
    
    face_t min_y = {
      vertex_type(vmin.x(), vmin.y(), vmin.z()),
      vertex_type(vmax.x(), vmin.y(), vmin.z()),
      vertex_type(vmax.x(), vmin.y(), vmax.z()),
      vertex_type(vmin.x(), vmin.y(), vmax.z())
    };

    face_t max_y = {
      vertex_type(vmin.x(), vmax.y(), vmin.z()),
      vertex_type(vmax.x(), vmax.y(), vmin.z()),
      vertex_type(vmax.x(), vmax.y(), vmax.z()),
      vertex_type(vmin.x(), vmax.y(), vmax.z())
    };

    face_t min_z = {
      vertex_type(vmin.x(), vmin.y(), vmin.z()),
      vertex_type(vmax.x(), vmin.y(), vmin.z()),
      vertex_type(vmax.x(), vmax.y(), vmin.z()),
      vertex_type(vmin.x(), vmax.y(), vmin.z())
    };

    face_t max_z = {
      vertex_type(vmin.x(), vmin.y(), vmax.z()),
      vertex_type(vmax.x(), vmin.y(), vmax.z()),
      vertex_type(vmax.x(), vmax.y(), vmax.z()),
      vertex_type(vmin.x(), vmax.y(), vmax.z())
    };

    auto write = [&](const face_t& face) {
      for(const vertex_type& v : face) {
        os << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
      }
      os << "f " << vtx_offset << " " << vtx_offset+1 << " " << vtx_offset+2 << " " 
                 << vtx_offset+3 << std::endl;
      os << "l " << vtx_offset+0 << " " << vtx_offset+1 << std::endl;
      os << "l " << vtx_offset+1 << " " << vtx_offset+2 << std::endl;
      os << "l " << vtx_offset+2 << " " << vtx_offset+3 << std::endl;
      os << "l " << vtx_offset+3 << " " << vtx_offset+0 << std::endl;
      vtx_offset += 4;
    };

    write(min_x);
    write(max_x);
    write(min_y);
    write(max_y);
    write(min_z);
    write(max_z);

  }

private:
  const entity_t* m_entity;
  std::array<vertex_type, 2> m_vertices;
  vertex_type m_center;
  vertex_array_type m_width;
  vertex_array_type m_iwidth;

  self_t* m_left_child{nullptr};
  self_t* m_right_child{nullptr};
  self_t* m_skip{nullptr};
};

template <typename T, typename U, size_t V>
std::ostream& operator<<(std::ostream& os, const AxisAlignedBoundingBox<T, U, V>& box) 
{
  box.dump(os);
  return os;
}

// default type and dim
template <typename T>
using AABB = AxisAlignedBoundingBox<T, float, 3>;


}
