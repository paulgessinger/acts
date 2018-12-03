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

namespace Acts {

template <typename box_t>
class Node : box_t
{

};


template <typename entity_t, size_t DIM, typename value_t>
class AxisAlignedBoundaryBox
{
private:
  using self_t = AxisAlignedBoundaryBox<entity_t, DIM, value_t>;
public:
  using vertex_t = ActsVector<value_t, DIM>;
  using vertex_array_t = Eigen::Array<value_t, DIM, 1>;

  // if we construct this with an entity, the entity can not be null
  AxisAlignedBoundaryBox(const entity_t& entity, const vertex_t& vmin, const vertex_t& vmax)
    : m_entity(&entity),
      m_vertices({vmin, vmax}),
      m_center((vmin + vmax)/2.)
  {
  }

  AxisAlignedBoundaryBox(const std::vector<self_t*>& boxes,
                         vertex_array_t envelope = vertex_array_t::Zero())
    : m_entity(nullptr)
  {
    assert(boxes.size() > 1);

    // figure out extent of boxes
    // use array for Eigen coefficient wise min/max
    using array_t = Eigen::Array<value_t, DIM, 1>;
    array_t vmax(array_t::Constant(-std::numeric_limits<value_t>::max()));
    array_t vmin(array_t::Constant(std::numeric_limits<value_t>::max()));

    for(size_t i=0;i<boxes.size();i++) {
      vmin = vmin.min(boxes[i]->m_vertices[0].array());
      vmax = vmax.max(boxes[i]->m_vertices[1].array());

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

    m_vertices[0] = std::move(vmin);
    m_vertices[1] = std::move(vmax);

    m_vertices[0] = m_vertices[0].array() - envelope;
    m_vertices[1] = m_vertices[1].array() + envelope;

    m_center = (m_vertices[0] + m_vertices[1])/2.;
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

  const vertex_t& center() const
  {
    return m_center;
    //return (m_vertices[0] + m_vertices[1]) /2.;
  }

  const vertex_t& min() const
  {
    return m_vertices[0];
  }

  const vertex_t& max() const
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
    using face_t = std::array<vertex_t, 4>;
    const vertex_t& vmin = m_vertices[0];
    const vertex_t& vmax = m_vertices[1];

    face_t min_x = {
      vertex_t(vmin.x(), vmin.y(), vmin.z()),
      vertex_t(vmin.x(), vmax.y(), vmin.z()),
      vertex_t(vmin.x(), vmax.y(), vmax.z()),
      vertex_t(vmin.x(), vmin.y(), vmax.z())
    };
    
    face_t max_x = {
      vertex_t(vmax.x(), vmin.y(), vmin.z()),
      vertex_t(vmax.x(), vmax.y(), vmin.z()),
      vertex_t(vmax.x(), vmax.y(), vmax.z()),
      vertex_t(vmax.x(), vmin.y(), vmax.z())
    };
    
    face_t min_y = {
      vertex_t(vmin.x(), vmin.y(), vmin.z()),
      vertex_t(vmax.x(), vmin.y(), vmin.z()),
      vertex_t(vmax.x(), vmin.y(), vmax.z()),
      vertex_t(vmin.x(), vmin.y(), vmax.z())
    };

    face_t max_y = {
      vertex_t(vmin.x(), vmax.y(), vmin.z()),
      vertex_t(vmax.x(), vmax.y(), vmin.z()),
      vertex_t(vmax.x(), vmax.y(), vmax.z()),
      vertex_t(vmin.x(), vmax.y(), vmax.z())
    };

    face_t min_z = {
      vertex_t(vmin.x(), vmin.y(), vmin.z()),
      vertex_t(vmax.x(), vmin.y(), vmin.z()),
      vertex_t(vmax.x(), vmax.y(), vmin.z()),
      vertex_t(vmin.x(), vmax.y(), vmin.z())
    };

    face_t max_z = {
      vertex_t(vmin.x(), vmin.y(), vmax.z()),
      vertex_t(vmax.x(), vmin.y(), vmax.z()),
      vertex_t(vmax.x(), vmax.y(), vmax.z()),
      vertex_t(vmin.x(), vmax.y(), vmax.z())
    };

    auto write = [&](const face_t& face) {
      for(const vertex_t& v : face) {
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
  std::array<vertex_t, 2> m_vertices;
  vertex_t m_center;

  self_t* m_left_child{nullptr};
  self_t* m_right_child{nullptr};
  self_t* m_skip{nullptr};
};

template <typename T, size_t U, typename V>
std::ostream& operator<<(std::ostream& os, const AxisAlignedBoundaryBox<T, U, V>& box) 
{
  box.dump(os);
  return os;
}

// default type and dim
template <typename T>
using AABB = AxisAlignedBoundaryBox<T, 3, float>;

}
