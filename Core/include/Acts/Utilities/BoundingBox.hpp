// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <limits>
#include <tuple>
#include <algorithm>

namespace Acts {

template <typename value_t, size_t DIM, size_t SIDES>
class Frustum
{
  //static_assert(DIM >= SIDES, "Cannot have DIM<=SIDES");
  using transform_t = Eigen::Transform<value_t, DIM, Eigen::Affine>;
  using translation_t = Eigen::Translation<value_t, DIM>;
public:
  using value_type = value_t;
  using vertex_type = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;

  static const size_t dim = DIM;
  static const size_t sides = SIDES;

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  Frustum(const vertex_type& origin, const vertex_type& dir, value_type opening_angle)
    : m_origin(origin)
  {
    using rotation_t = Eigen::Rotation2D<value_type>;

    static_assert(SIDES==2, "2D frustum can only have 2 sides");
    assert(opening_angle < M_PI);

    translation_t translation(origin);
    value_type angle = VectorHelpers::phi(dir);
    Eigen::Rotation2D<value_type> rot(angle);
    m_transform = transform_t(translation * rot);

    value_type normal_angle = 0.5*M_PI - 0.5*opening_angle;
    vertex_type normal1 = rotation_t(normal_angle) * vertex_type::UnitX();
    vertex_type normal2 = rotation_t(-normal_angle) * vertex_type::UnitX();

    m_normals = {
      rot * vertex_type::UnitX(),
      rot * normal1,
      rot * normal2
    };

  }
  
  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  Frustum(const vertex_type& origin, const vertex_type& dir, value_type opening_angle)
    : m_origin(origin)
  {
    static_assert(SIDES>2, "3D frustum must have 3 or more sides");
    assert(opening_angle < M_PI);
    using angle_axis_t = Eigen::AngleAxis<value_type>;

    const vertex_type ldir = vertex_type::UnitZ();
    const vertex_type lup = vertex_type::UnitX();

    m_transform = (
        Eigen::Quaternion<value_type>().setFromTwoVectors(ldir, dir)
    );

    m_normals[0] = ldir;

    const value_type phi_sep = 2*M_PI / sides;
    transform_t rot;
    rot = angle_axis_t(phi_sep, ldir);

    value_type half_opening_angle = opening_angle/2.;
    auto calculate_normal = [&ldir, &half_opening_angle](const vertex_type& out) 
      -> vertex_type
    {
      const vertex_type tilt_axis = -1*out.cross(ldir);
      return (-1 * (angle_axis_t(half_opening_angle, tilt_axis) * out)).normalized();
    };

    vertex_type current_outward = lup;
    m_normals[1] = calculate_normal(current_outward);

    for(size_t i=1;i<sides;i++) {
      current_outward = rot * current_outward;
      m_normals[i+1] = calculate_normal(current_outward);
    }

    for(auto& normal : m_normals) {
      normal = m_transform * normal;
    }

  }

  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  std::ofstream& obj(std::ofstream& os, size_t& n_vtx, value_type far_distance = 10) const
  {
    //using angle_axis_t = Eigen::AngleAxis<value_type>;
    static_assert(DIM == 3, "OBJ is only supported in 3D");

    //auto draw_line = [&os, &n_vtx](const vertex_type& start,
                                   //const vertex_type& end)
    //{
    
      //os << "v " << start.x() << " " << start.y() << " " << start.z() << std::endl;
      //os << "v " << end.x() << " " << end.y() << " " << end.z() << std::endl;

      //os << "l " << n_vtx << " " << n_vtx+1 << std::endl;
      //n_vtx += 2;
    //};

    auto draw_plane = [&os, &n_vtx](const std::vector<vertex_type>& vtxs)
    {
      for(const auto& v : vtxs) {
        os << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
      }

      os << "f";
      for(size_t i=0;i<vtxs.size();i++) {
        os << " " << n_vtx+i;
      }
      os << std::endl;
      n_vtx += vtxs.size();
    };
    

    //translation_t trans(m_origin);
    //const auto& normal = m_normals.at(0);
    //transform_t l2g(trans);
    //l2g = l2g * Eigen::Quaternion<value_type>().setFromTwoVectors(vertex_type::UnitZ(), normal);

    //std::vector<vertex_type> vtxs = {
      //{10, 10, 0},
      //{10, -10, 0},
      //{-10, -10, 0},
      //{-10, 10, 0},
    //};

    //for(auto& v : vtxs) {
      //v = l2g * (v * 2.);
    //}
   
      
    // iterate around normals, calculate cross with "far" plane
    // to get intersection lines.
    // Work in O = (0, 0) and shift draw vertices at the end
    vertex_type far_normal = m_normals[0]; // far has same normal as pseudo-near
    vertex_type far_center = m_normals[0] * far_distance; // center defined as 10 "forward"
    std::array<std::pair<vertex_type, vertex_type>, SIDES> planeFarIXs;

    auto ixPlanePlane = [](const auto& n1, const auto& p1, const auto& n2, const auto& p2) 
      -> std::pair<vertex_type, vertex_type> {
      const vertex_type m = n1.cross(n2).normalized();
      const double j = ( n2.dot( p2-p1 ) ) / ( n2.dot( n1.cross(m) )  );
      const vertex_type q = p1 + j * n1.cross(m);
      return {m, q};
    };

    auto ixLineLine = [](const auto& p1, const auto& d1, const auto& p2, const auto& d2)
      -> vertex_type {
       return p1 + ( ( (p2-p1).cross(d2) ).norm() / (d1.cross(d2)).norm() ) * d1;
    };

    // skip i=0 <=> pseudo-near
    for(size_t i=1;i<m_normals.size();i++) {
      //vertex_type cross = m_normals[i].cross(far_normal).normalize();
      const auto ixLine = ixPlanePlane(far_normal, far_center, m_normals[i], vertex_type::Zero());
      //std::cout << ixLine.first.transpose() << " " << ixLine.second.transpose() << std::endl;
      //std::cout << "idx: " << i << std::endl;
      planeFarIXs.at(i-1) = ixLine; 

      //draw_line(m_origin+ixLine.second-ixLine.first*15, m_origin+ixLine.second+ixLine.first*15);
    }

    std::array<vertex_type, SIDES> points;

    for(size_t i=0;i<planeFarIXs.size();i++) {
      size_t j = (i+1) % planeFarIXs.size();
      //std::cout << i << " " << j << std::endl;
      //const vertex_type ix = planeFarIXs.at(i).intersection(planeFarIXs.at(j));
      const auto& l1 = planeFarIXs.at(i);
      const auto& l2 = planeFarIXs.at(j);
      const vertex_type ix = m_origin + ixLineLine(l1.second, l1.first, l2.second, l2.first);
      //std::cout << "llix: " << ix.transpose() << std::endl;

      points.at(i) = ix;

      //os << "v " << ix.x() << " " << ix.y() << " " << ix.z() << std::endl;
      //n_vtx++;
    }

    for(size_t i=0;i<points.size();i++) {
      size_t j = (i+1) % points.size();
      draw_plane({m_origin, points.at(i), points.at(j)});
    }




    return os;
  }

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  std::ofstream& svg(std::ofstream& os, value_type w, value_type h, value_type far_distance = 1) const
  {
    static_assert(DIM == 2, "SVG is only supported in 2D");
    const value_type unit = 20;

    vertex_type mid(w/2., h/2.);

    //using transform_t = Eigen::Transform<value_type, 2, Eigen::Affine>;
    transform_t trf = transform_t::Identity();
    //Eigen::Matrix<value_type, 2, 2> flip;
    //flip << 1,  0,
            //0, -1;
    trf.translate(mid);
    trf = trf * Eigen::Scaling(vertex_type(1, -1));
    trf.scale(unit);

    //trf =  * Eigen::Translation<value_type, 2>(mid);

    std::array<std::string, 3> colors({"orange", "blue", "red"});

    auto draw_line = [&](const vertex_type& left_, 
                           const vertex_type& right_,
                           std::string color, 
                           size_t width) {

      vertex_type left = trf*left_;
      vertex_type right = trf*right_;
      os << "<line ";

      os << "x1=\"" << left.x() << "\" ";
      os << "y1=\"" << left.y() << "\" "; 
      os << "x2=\"" << right.x() << "\" ";
      os << "y2=\"" << right.y() << "\" ";

      os <<" stroke=\"" << color << "\" stroke-width=\"" << width << "\"/>\n";

    };

    auto draw_point = [&](const vertex_type& p_, std::string color, size_t r) {
      vertex_type p = trf*p_;
      os << "<circle ";
      os << "cx=\"" << p.x() << "\" cy=\"" << p.y() << "\" r=\"" << r << "\"";
      os << " fill=\"" << color << "\"";
      os << "/>\n";
    };
    
    using vec3 = ActsVector<value_type, 3>;
    auto ixLineLine = [](const vertex_type& p1_2, const vertex_type& d1_2, const vertex_type& p2_2, const vertex_type& d2_2)
      -> vertex_type {
      
      const vec3 p1(p1_2.x(), p1_2.y(), 0);
      const vec3 p2(p2_2.x(), p2_2.y(), 0);
      const vec3 d1(d1_2.x(), d1_2.y(), 0);
      const vec3 d2(d2_2.x(), d2_2.y(), 0);

      vec3 num = (p2-p1).cross(d2);
      vec3 den = d1.cross(d2);

      value_type f = 1.;
      value_type dot = num.normalized().dot(den.normalized());
      if(std::abs(dot)-1 < 1e-9 && dot < 0) {
        f = -1.;
      }

      const vec3 p = p1 + f*( num.norm() / den.norm() ) * d1;
      assert(std::abs(p.z()) < 1e-9);
      return {p.x(), p.y()};
    };

    const vertex_type far_dir = {m_normals[0].y(), -m_normals[0].x()};
    const vertex_type far_point = m_normals[0]*far_distance;

    //draw_line(m_origin+far_point - 5*far_dir, m_origin+far_point + 5*far_dir, "black", 2);

    std::array<vertex_type, 2> points;

    for(size_t i=1;i<m_normals.size();i++) {
      vertex_type plane_dir(m_normals[i].y(), -m_normals[i].x());
      
      const vertex_type ix = ixLineLine(far_point, far_dir, {0, 0}, plane_dir);
      draw_point(m_origin + ix, "black", 3);
      draw_line(m_origin, m_origin + ix, "black", 2);
      points.at(i-1) = ix;
    }
  
    draw_line(m_origin + points.at(0), m_origin + points.at(1), "black", 2);

    //for(size_t i=1;i<points.size();i++) {
      //size_t j = (i+1)%points.size();
    //}
    

    draw_line(m_origin, m_origin+m_normals[0]*2, colors[0], 3);
    
    draw_point({0, 0}, "yellow", 5);
    draw_point(m_origin, "green", 5);


    return os;

  }

private:
  vertex_type m_origin;
  vertex_type m_dir;
  transform_t m_transform;
  // need one more for direction we're facing
  std::array<vertex_type, SIDES+1> m_normals;

};




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

    //os << m_idir;

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
  
  // strong type helper, not public
  template <typename T, typename P>
  class NamedType
  {
    public:
      explicit NamedType(const T& value) : m_value(value) {}
      explicit NamedType(T&& value) : m_value(std::move(value)) {}
      T& get() { return m_value; }
      const T& get() const { return m_value; }
    private:
      T m_value;
  };

  struct SizeParameter{};


public:
  using vertex_type = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;
  using entity_type = entity_t;
  using value_type = value_t;
  
  using Size = NamedType<vertex_type, struct SizeParameter>;

  static const size_t dim = DIM;

  // if we construct this with an entity, the entity can not be null
  AxisAlignedBoundingBox(const entity_t& entity, const vertex_type& vmin, const vertex_type& vmax)
    : m_entity(&entity),
      m_vmin(vmin),
      m_vmax(vmax),
      m_center((vmin + vmax)/2.),
      m_width(vmax - vmin),
      m_iwidth(1/m_width)
  {
  }

  AxisAlignedBoundingBox(const entity_t&    entity,
                         const vertex_type& center,
                         const Size& size)
    : m_entity(&entity)
    , m_vmin(center - size.get() * 0.5)
    , m_vmax(center + size.get() * 0.5)
    , m_center(center)
    , m_width(size.get())
    , m_iwidth(1 / m_width)
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

    std::tie(m_vmin, m_vmax) = wrap(boxes, envelope);

    m_center = (m_vmin + m_vmax)/2.;
    m_width = m_vmax - m_vmin;
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
      vmin = vmin.min(boxes[i]->min().array());
      vmax = vmax.max(boxes[i]->max().array());
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
    vertex_array_type t = (point - m_vmin).array() * m_iwidth;
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

    const vertex_type& origin = ray.origin();
    const vertex_array_type& idir = ray.idir();
    // this is NaN origin is on box boundary and ray is parallel to
    // that boundary, since 0*inf = NaN. 
    vertex_array_type t0s = (m_vmin - origin).array() * idir;
    vertex_array_type t1s = (m_vmax - origin).array() * idir;
    

    // this is non-compliant with IEEE-754-2008, NaN gets propagated through
    // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=564
    // this means that rays parallel to boundaries might not be considered 
    // to intersect.
    vertex_array_type tsmaller = t0s.min(t1s);
    vertex_array_type tbigger = t0s.max(t1s);


    value_type tmin = tsmaller.maxCoeff();
    value_type tmax = tbigger.minCoeff();

    //std::cout << "--- " << ray << " -> " << this << std::endl;
    //std::cout << "t0s:\n" << t0s << "\nt1s\n" << t1s << "\n\n";
    //std::cout << "tsmaller:\n" << tsmaller << "\ntbigger\n" << tbigger << "\n\n";
    //std::cout << "tmin:" << tmin << "\ntmax:" << tmax << std::endl;
    
    return tmin < tmax && tmax > 0.0;// ((tmin > 0.0 && tmax > 0.0) || (tmin < 0.0 && tmax > 0.0));
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
  }

  const vertex_type& min() const
  {
    return m_vmin;
  }

  const vertex_type& max() const
  {
    return m_vmax;
  }

  std::ostream& dump(std::ostream& os) const
  {
    os << "AABB(ctr=(";

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
      os << m_vmin[i];
    }

    os << ") vmax=(";

    for(size_t i=0;i<DIM;i++) {
      if (i>0) {
        os << ", ";
      }
      os << m_vmax[i];
    }

    os << "))";

    return os;
  }

  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  void obj(std::ostream& os, size_t& vtx_offset) const 
  {
    static_assert(DIM == 3, "OBJ output only supported in 3D");
    using face_t = std::array<vertex_type, 4>;
    const vertex_type& vmin = m_vmin;
    const vertex_type& vmax = m_vmax;

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

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  std::ofstream& svg(std::ofstream& os, value_type w, value_type h, value_type unit = 10) const
  {
    static_assert(DIM == 2, "SVG is only supported in 2D");

    vertex_type mid(w/2., h/2.);

    transform_t trf = transform_t::Identity();
    trf.translate(mid);
    trf = trf * Eigen::Scaling(vertex_type(1, -1));
    trf.scale(unit);


    auto draw_line = [&](const vertex_type& left_, 
                           const vertex_type& right_,
                           std::string color, 
                           size_t width) {

      vertex_type left = trf*left_;
      vertex_type right = trf*right_;
      os << "<line ";

      os << "x1=\"" << left.x() << "\" ";
      os << "y1=\"" << left.y() << "\" "; 
      os << "x2=\"" << right.x() << "\" ";
      os << "y2=\"" << right.y() << "\" ";

      os <<" stroke=\"" << color << "\" stroke-width=\"" << width << "\"/>\n";

    };

    auto draw_point = [&](const vertex_type& p_, std::string color, size_t r) {
      vertex_type p = trf*p_;
      os << "<circle ";
      os << "cx=\"" << p.x() << "\" cy=\"" << p.y() << "\" r=\"" << r << "\"";
      os << " fill=\"" << color << "\"";
      os << "/>\n";
    };


  }

private:
  const entity_t* m_entity;
  vertex_type m_vmin;
  vertex_type m_vmax;
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
