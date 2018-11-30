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

namespace Acts {

template <typename T>
struct node {
  node(T* self_, T* next_, T* skip_) : self(self_), next(next_), skip(skip_) {}
  T* self;
  T* next;
  T* skip;
};

  template <typename T>
std::ostream& operator<<(std::ostream& os, const node<T>& n) 
{
  os << "node(" << n.self << ", " << n.next << ", " << n.skip << ")";
  return os;
}

template <typename entity_t, size_t DIM, typename value_t>
class AxisAlignedBoundaryBox
{
private:
  using self_t = AxisAlignedBoundaryBox<entity_t, DIM, value_t>;
public:
  using vertex_t = ActsVector<value_t, DIM>;

  // if we construct this with an entity, the entity can not be null
  AxisAlignedBoundaryBox(const entity_t& entity, const vertex_t& vmin, const vertex_t& vmax)
    : m_entity(&entity),
      m_vertices({vmin, vmax})
  {
    // no children, need no capacity
    m_children.shrink_to_fit();
  }

  AxisAlignedBoundaryBox(std::vector<self_t> boxes) 
    : m_entity(nullptr),
      m_children(std::move(boxes))
  {
  
    assert(m_children.size() > 1);
    m_children.shrink_to_fit(); // we'll never need more

    // figure out extent of boxes
    //vertex_t vmax(vertex_t::Constant(std::numeric_limits<value_t>::min()));
    //vertex_t vmin(vertex_t::Constant(std::numeric_limits<value_t>::max()));
    //auto vmax = vertex_t(vertex_t::Constant(std::numeric_limits<value_t>::min())).array();
    //auto vmin = vertex_t(vertex_t::Constant(std::numeric_limits<value_t>::max())).array();

    // use array for Eigen coefficient wise min/max
    using array_t = Eigen::Array<value_t, DIM, 1>;
    array_t vmax(array_t::Constant(std::numeric_limits<value_t>::min()));
    array_t vmin(array_t::Constant(std::numeric_limits<value_t>::max()));

    for(const self_t& child : m_children) {
      //vmin = child.m_vertices[0].array().min(vmin.array());
      //vmax = child.m_vertices[1].array().max(vmax.array());
      vmin = vmin.min(child.m_vertices[0].array());
      vmax = vmax.max(child.m_vertices[1].array());
    }

    m_vertices[0] = std::move(vmin);
    m_vertices[1] = std::move(vmax);
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

  std::ostream& dump(std::ostream& os) const
  {
    os << "AxisAlignedBoundingBox(vmin=(";
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

  std::vector<node<const self_t>>
  flatten() const {
    std::vector<node<const self_t>> result;
    flatten_impl(result, nullptr);
    return result;
  }




private:
  void flatten_impl(std::vector<node<const self_t>>& list, const self_t* skip) const
  {
    for(size_t i=0;i<m_children.size()-1;i++) {
      list.emplace_back(&m_children[i], &m_children[i+1], skip);
      m_children[i].flatten_impl(list, &m_children[i+1]);
    }
    list.emplace_back(&m_children.back(), skip, skip);
  }

  const entity_t* m_entity;
  std::array<vertex_t, 2> m_vertices;

  // children are owned via vector
  std::vector<self_t> m_children;
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

