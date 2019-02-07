// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/VariantData.hpp"

template <int N>
Acts::ConvexPolygonBounds<N>::ConvexPolygonBounds(
    const std::vector<Acts::Vector2D>& vertices)
  : m_vertices(), m_boundingBox(makeBoundingBox(vertices))
{
  assert(vertices.size() == N);
  for (size_t i = 0; i < N; i++) { m_vertices[i] = vertices[i]; }
  assert(convex());
}

template <int N>
Acts::ConvexPolygonBounds<N>::ConvexPolygonBounds(const vertex_array& vertices)
  : m_vertices(vertices), m_boundingBox(makeBoundingBox(vertices))
{
  assert(convex());
}

template <int N>
template <typename coll_t>
Acts::RectangleBounds
Acts::ConvexPolygonBounds<N>::makeBoundingBox(const coll_t& vertices)
{
  assert(vertices.size() == N);
  Vector2D vmax, vmin;
  vmax = vertices[0];
  vmin = vertices[0];

  for (size_t i = 1; i < N; i++) {
    vmax = vmax.cwiseMax(vertices[i]);
    vmin = vmin.cwiseMin(vertices[i]);
  }

  return {vmin, vmax};
}

template <int N>
Acts::ConvexPolygonBounds<N>*
Acts::ConvexPolygonBounds<N>::clone() const
{
  return new ConvexPolygonBounds<N>(*this);
}

template <int N>
Acts::SurfaceBounds::BoundsType
Acts::ConvexPolygonBounds<N>::type() const
{
  return SurfaceBounds::ConvexPolygon;
}

template <int N>
bool
Acts::ConvexPolygonBounds<N>::inside(const Acts::Vector2D&      lpos,
                                     const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(lpos, m_vertices);
}

template <int N>
double
Acts::ConvexPolygonBounds<N>::distanceToBoundary(
    const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(lpos, m_vertices);
}

template <int N>
std::vector<Acts::Vector2D>
Acts::ConvexPolygonBounds<N>::vertices() const
{
  return {m_vertices.begin(), m_vertices.end()};
}

template <int N>
std::vector<TDD_real_t>
Acts::ConvexPolygonBounds<N>::valueStore() const
{
  std::vector<TDD_real_t> values;
  for (const auto& vtx : m_vertices) {
    values.push_back(vtx.x());
    values.push_back(vtx.y());
  }
  return values;
}

template <int N>
Acts::variant_data
Acts::ConvexPolygonBounds<N>::toVariantData() const
{
  using namespace std::string_literals;
  variant_map payload;
  payload["sides"] = N;

  variant_vector vertices;
  for (const auto& vtx : m_vertices) { vertices.push_back(to_variant(vtx)); }

  payload["vertices"] = vertices;

  variant_map data;
  data["type"]    = "ConvexPolygonBounds";
  data["payload"] = payload;

  return data;
}

template <int N>
const Acts::RectangleBounds&
Acts::ConvexPolygonBounds<N>::boundingBox() const
{
  return m_boundingBox;
}

template <int N>
bool
Acts::ConvexPolygonBounds<N>::convex() const
{
  for (size_t i = 0; i < N; i++) {
    size_t          j = (i + 1) % N;
    const Vector2D& a = m_vertices[i];
    const Vector2D& b = m_vertices[j];

    const Vector2D ab     = b - a;
    const Vector2D normal = Vector2D(ab.y(), -ab.x()).normalized();

    bool first = true;
    bool ref;
    // loop over all other vertices
    for (size_t k = 0; k < N; k++) {
      if (k == i || k == j) { continue; }

      const Vector2D& c   = m_vertices[k];
      double          dot = normal.dot(c - a);

      if (first) {
        ref   = std::signbit(dot);
        first = false;
        continue;
      }

      if (std::signbit(dot) != ref) { return false; }
    }
  }
  return true;
}

template <int N>
std::ostream&
Acts::ConvexPolygonBounds<N>::dump(std::ostream& sl) const
{
  sl << "Acts::ConvexPolygonBounds<" << num_vertices << ">: vertices: [x, y]\n";
  for (size_t i = 0; i < m_vertices.size(); i++) {
    const auto& vtx = m_vertices[i];
    if (i > 0) {
      sl << ",";
      sl << "\n";
    }
    sl << "[" << vtx.x() << ", " << vtx.y() << "]";
  }
  return sl;
}
