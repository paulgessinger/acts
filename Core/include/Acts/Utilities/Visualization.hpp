// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

template <typename T>
struct ply_helper
{
  using value_type  = T;
  using vertex_type = ActsVector<value_type, 3>;

  using face_type  = std::vector<size_t>;
  using color_type = std::array<int, 3>;

  void
  vertex(const vertex_type& vtx, color_type color = {120, 120, 120})
  {
    m_vertices.emplace_back(vtx, color);
  }

  template <typename coll_t>
  void
  face(const coll_t& vtxs, color_type color = {120, 120, 120})
  {

    static_assert(std::is_same<typename coll_t::value_type, vertex_type>::value,
                  "not a collection of vertex_type");

    face_type idxs;
    idxs.reserve(vtxs.size());
    for (const vertex_type& vtx : vtxs) {
      vertex(vtx, color);
      idxs.push_back(m_vertices.size() - 1);
    }
    m_faces.push_back(std::move(idxs));
  }

  void
  line(const vertex_type& a,
       const vertex_type& b,
       color_type         color = {120, 120, 120})
  {
    vertex(a, color);
    size_t idx_a = m_vertices.size() - 1;
    vertex(b, color);
    size_t idx_b = m_vertices.size() - 1;
    m_edges.emplace_back(std::make_pair(std::make_pair(idx_a, idx_b), color));
  }

  void
  write(std::ostream& os) const
  {
    os << "ply\n";
    os << "format ascii 1.0\n";
    os << "element vertex " << m_vertices.size() << "\n";
    os << "property float x\n";
    os << "property float y\n";
    os << "property float z\n";
    os << "property uchar red\n";
    os << "property uchar green\n";
    os << "property uchar blue\n";
    os << "element face " << m_faces.size() << "\n";
    os << "property list uchar int vertex_index\n";
    os << "element edge " << m_edges.size() << "\n";
    os << "property int vertex1\n";
    os << "property int vertex2\n";
    os << "property uchar red\n";
    os << "property uchar green\n";
    os << "property uchar blue\n";
    os << "end_header\n";

    for (const std::pair<vertex_type, color_type>& vtx : m_vertices) {
      os << vtx.first.x() << " " << vtx.first.y() << " " << vtx.first.z()
         << " ";
      os << vtx.second[0] << " " << vtx.second[1] << " " << vtx.second[2]
         << "\n";
    }

    for (const face_type& fc : m_faces) {
      os << fc.size();
      for (size_t i = 0; i < fc.size(); i++) {
        os << " " << fc[i];
      }
      os << "\n";
    }

    for (const std::pair<std::pair<size_t, size_t>, color_type>& edge :
         m_edges) {
      std::pair<size_t, size_t> idxs = edge.first;
      os << idxs.first << " " << idxs.second << " ";
      os << edge.second[0] << " " << edge.second[1] << " " << edge.second[2]
         << "\n";
    }
  }

  void
  clear()
  {
    m_vertices.clear();
    m_faces.clear();
    m_edges.clear();
  }

private:
  std::vector<std::pair<vertex_type, color_type>> m_vertices;
  std::vector<face_type> m_faces;
  std::vector<std::pair<std::pair<size_t, size_t>, color_type>> m_edges;
};

template <typename T>
std::ostream&
operator<<(std::ostream& os, const ply_helper<T>& ply)
{
  ply.write(os);
  return os;
}

template <typename T>
struct obj_helper
{
  using value_type  = T;
  using vertex_type = ActsVector<value_type, 3>;

  using face_type = std::vector<size_t>;
  // unsupported
  using color_type = std::array<int, 3>;

  void
  vertex(const vertex_type& vtx, color_type color = {120, 120, 120})
  {
    (void)color;
    m_vertices.push_back(vtx);
  }

  void
  line(const vertex_type& /*a*/,
       const vertex_type& /*b*/,
       color_type         /*color*/)
  {
  // not implemented
  }

  template <typename coll_t>
  void
  face(const coll_t& vtxs, color_type color = {120, 120, 120})
  {
    (void)color;
    static_assert(std::is_same<typename coll_t::value_type, vertex_type>::value,
                  "not a collection of vertex_type");

    face_type idxs;
    idxs.reserve(vtxs.size());
    for (const vertex_type& vtx : vtxs) {
      vertex(vtx);
      idxs.push_back(m_vertices.size() - 1);
    }
    m_faces.push_back(std::move(idxs));
  }

  void
  write(std::ostream& os) const
  {
    for (const vertex_type& vtx : m_vertices) {
      os << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
    }

    for (const face_type& fc : m_faces) {
      os << "f";
      for (size_t i = 0; i < fc.size(); i++) {
        os << " " << fc[i] + 1;
      }
      os << "\n";
    }
  }

  void
  clear()
  {
    m_vertices.clear();
    m_faces.clear();
  }

private:
  std::vector<vertex_type> m_vertices;
  std::vector<face_type>   m_faces;
};
}
