// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::AxisAlignedBoundingBox(
    const entity_t*    entity,
    const vertex_type& vmin,
    const vertex_type& vmax)
  : m_entity(entity)
  , m_vmin(vmin)
  , m_vmax(vmax)
  , m_center((vmin + vmax) / 2.)
  , m_width(vmax - vmin)
  , m_iwidth(1 / m_width)
{
}

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::AxisAlignedBoundingBox(
    const entity_t*    entity,
    const vertex_type& center,
    const Size&        size)
  : m_entity(entity)
  , m_vmin(center - size.get() * 0.5)
  , m_vmax(center + size.get() * 0.5)
  , m_center(center)
  , m_width(size.get())
  , m_iwidth(1 / m_width)
{
}

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::AxisAlignedBoundingBox(
    const std::vector<self_t*>& boxes,
    vertex_array_type           envelope)
  : m_entity(nullptr)
{
  assert(boxes.size() > 1);

  for (size_t i = 0; i < boxes.size(); i++) {
    if (i < boxes.size() - 1) {
      // set next on i to i+1
      boxes[i]->setSkip(boxes[i + 1]);
    } else {
      // make sure last is set to nullptr, this marks end
      // boxes[i]->m_next = nullptr;
      boxes[i]->setSkip(nullptr);
    }
  }

  m_left_child  = boxes.front();
  m_right_child = boxes.back();
  m_skip        = nullptr;

  std::tie(m_vmin, m_vmax) = wrap(boxes, envelope);

  m_center = (m_vmin + m_vmax) / 2.;
  m_width  = m_vmax - m_vmin;
  m_iwidth = 1 / m_width;
}

template <typename entity_t, typename value_t, size_t DIM>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::wrap(
    const std::vector<const self_t*>& boxes,
    vertex_array_type                 envelope)
{
  assert(boxes.size() > 1);
  // figure out extent of boxes
  // use array for Eigen coefficient wise min/max
  vertex_array_type vmax(
      vertex_array_type::Constant(std::numeric_limits<value_type>::lowest()));
  vertex_array_type vmin(
      vertex_array_type::Constant(std::numeric_limits<value_type>::max()));

  for (size_t i = 0; i < boxes.size(); i++) {
    vmin = vmin.min(boxes[i]->min().array());
    vmax = vmax.max(boxes[i]->max().array());
  }

  vmax += envelope;
  vmin -= envelope;

  return {vmin, vmax};
}

template <typename entity_t, typename value_t, size_t DIM>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::wrap(
    const std::vector<self_t*>& boxes,
    vertex_array_type           envelope)
{
  assert(boxes.size() > 1);
  std::vector<const self_t*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  std::transform(boxes.begin(),
                 boxes.end(),
                 std::back_inserter(box_ptrs),
                 [](const auto* box) { return box; });
  return wrap(box_ptrs, envelope);
}

template <typename entity_t, typename value_t, size_t DIM>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::wrap(
    const std::vector<self_t>& boxes,
    vertex_array_type          envelope)
{
  assert(boxes.size() > 1);
  std::vector<const self_t*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  std::transform(boxes.begin(),
                 boxes.end(),
                 std::back_inserter(box_ptrs),
                 [](auto& box) { return &box; });
  return wrap(box_ptrs, envelope);
}

template <typename entity_t, typename value_t, size_t DIM>
bool
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::intersect(
    const vertex_type& point) const
{
  vertex_array_type t = (point - m_vmin).array() * m_iwidth;
  return t.minCoeff() >= 0 && t.maxCoeff() < 1;
}

template <typename entity_t, typename value_t, size_t DIM>
bool
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::intersect(
    const Ray<value_type, DIM>& ray) const
{

  const vertex_type&       origin = ray.origin();
  const vertex_array_type& idir   = ray.idir();
  // this is NaN origin is on box boundary and ray is parallel to
  // that boundary, since 0*inf = NaN.
  vertex_array_type t0s = (m_vmin - origin).array() * idir;
  vertex_array_type t1s = (m_vmax - origin).array() * idir;

  // this is non-compliant with IEEE-754-2008, NaN gets propagated through
  // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=564
  // this means that rays parallel to boundaries might not be considered
  // to intersect.
  vertex_array_type tsmaller = t0s.min(t1s);
  vertex_array_type tbigger  = t0s.max(t1s);

  value_type tmin = tsmaller.maxCoeff();
  value_type tmax = tbigger.minCoeff();

  // std::cout << "--- " << ray << " -> " << this << std::endl;
  // std::cout << "t0s:\n" << t0s << "\nt1s\n" << t1s << "\n\n";
  // std::cout << "tsmaller:\n" << tsmaller << "\ntbigger\n" << tbigger <<
  // "\n\n";
  // std::cout << "tmin:" << tmin << "\ntmax:" << tmax << std::endl;

  return tmin < tmax
      && tmax
      > 0.0;  // ((tmin > 0.0 && tmax > 0.0) || (tmin < 0.0 && tmax > 0.0));
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t sides>
bool
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::intersect(
    const Frustum<value_type, DIM, sides>& fr) const
{
  // std::cout << __FUNCTION__ << std::endl;
  // std::cout << "sides: " << sides << std::endl;

  const auto&             normals = fr.normals();
  const vertex_array_type vmin    = m_vmin - fr.origin();
  const vertex_array_type vmax    = m_vmax - fr.origin();

  // std::cout << "vmin: " << vmin.transpose() << "\nvmax: " <<
  // vmax.transpose()
  //<< std::endl;

  vertex_type p_vtx;
  // for loop, we could eliminate this, probably,
  // but sides+1 is known at compile time, so the compiler
  // will most likely unroll the loop
  for (size_t i = 0; i < sides + 1; i++) {
    const vertex_type& normal = normals[i];

    // for (size_t j=0;j<DIM;j++) {
    // p_vtx[j] = normal[j] < 0 ? vmin[j] : vmax[j];
    // std::cout << p_vtx[j] << std::endl;
    //}

    p_vtx = (normal.array() < 0).template cast<value_t>() * vmin
        + (normal.array() >= 0).template cast<value_t>() * vmax;
    // std::cout << p_vtx.transpose() << std::endl;

    if (p_vtx.dot(normal) < 0) {
      // p vertex is outside on this plane, box must be outside
      return false;
    }
  }

  // If we get here, no p-vertex was outside, so box intersects or is
  // contained. We don't care, so report 'intsersect'
  return true;
}

template <typename entity_t, typename value_t, size_t DIM>
void
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::setSkip(self_t* skip)
{
  // set next on this
  m_skip = skip;
  // find last child and set its skip
  if (m_right_child != nullptr) {
    m_right_child->setSkip(skip);
  }
}

template <typename entity_t, typename value_t, size_t DIM>
const Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>*
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::getLeftChild() const
{
  return m_left_child;
}

template <typename entity_t, typename value_t, size_t DIM>
const Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>*
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::getSkip() const
{
  return m_skip;
}

template <typename entity_t, typename value_t, size_t DIM>
bool
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::hasEntity() const
{
  return m_entity != nullptr;
}

template <typename entity_t, typename value_t, size_t DIM>
const entity_t*
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::entity() const
{
  return m_entity;
}

template <typename entity_t, typename value_t, size_t DIM>
void
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::setEntity(
    const entity_t* entity)
{
  m_entity = entity;
}

template <typename entity_t, typename value_t, size_t DIM>
const typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::
    vertex_type&
    Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::center() const
{
  return m_center;
}

template <typename entity_t, typename value_t, size_t DIM>
const typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::
    vertex_type&
    Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::min() const
{
  return m_vmin;
}

template <typename entity_t, typename value_t, size_t DIM>
const typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::
    vertex_type&
    Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::max() const
{
  return m_vmax;
}

template <typename entity_t, typename value_t, size_t DIM>
std::ostream&
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::dump(
    std::ostream& os) const
{
  os << "AABB(ctr=(";

  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_center[i];
  }

  os << ") vmin=(";
  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_vmin[i];
  }

  os << ") vmax=(";

  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_vmax[i];
  }

  os << "))";

  return os;
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 3, int>>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transformVertices(
    const transform_type& trf) const
{
  // we need to enumerate all the vertices, transform,
  // and then recalculate min and max

  std::array<vertex_type, 8> vertices({{
      {m_vmin.x(), m_vmin.y(), m_vmin.z()},
      {m_vmin.x(), m_vmax.y(), m_vmin.z()},
      {m_vmax.x(), m_vmax.y(), m_vmin.z()},
      {m_vmax.x(), m_vmin.y(), m_vmin.z()},
      {m_vmin.x(), m_vmin.y(), m_vmax.z()},
      {m_vmin.x(), m_vmax.y(), m_vmax.z()},
      {m_vmax.x(), m_vmax.y(), m_vmax.z()},
      {m_vmax.x(), m_vmin.y(), m_vmax.z()},
  }});

  vertex_type vmin = trf * vertices[0];
  vertex_type vmax = trf * vertices[0];

  for (size_t i = 1; i < 8; i++) {
    const vertex_type vtx = trf * vertices[i];
    vmin                  = vmin.cwiseMin(vtx);
    vmax                  = vmax.cwiseMax(vtx);
  }

  return {vmin, vmax};
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 2, int>>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transformVertices(
    const transform_type& trf) const
{
  // we need to enumerate all the vertices, transform,
  // and then recalculate min and max

  std::array<vertex_type, 4> vertices({{{m_vmin.x(), m_vmin.y()},
                                        {m_vmin.x(), m_vmax.y()},
                                        {m_vmax.x(), m_vmax.y()},
                                        {m_vmax.x(), m_vmin.y()}}});

  vertex_type vmin = trf * vertices[0];
  vertex_type vmax = trf * vertices[0];

  for (size_t i = 1; i < 4; i++) {
    const vertex_type vtx = trf * vertices[i];
    vmin                  = vmin.cwiseMin(vtx);
    vmax                  = vmax.cwiseMax(vtx);
  }

  return {vmin, vmax};
}

template <typename entity_t, typename value_t, size_t DIM>
void
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transform(
    const transform_type& trf)
{
  std::tie(m_vmin, m_vmax) = transformVertices(trf);
}

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transformed(
    const transform_type& trf) const
{
  vertex_type vmin, vmax;
  std::tie(vmin, vmax) = transformVertices(trf);
  return self_t(m_entity, vmin, vmax);
}

template <typename entity_t, typename value_t, size_t DIM>
template <typename helper_t, size_t D, std::enable_if_t<D == 3, int>>
void
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::draw(
    helper_t& helper,
    std::array<int, 3> color,
    const transform_type& trf) const
{
  static_assert(DIM == 3, "PLY output only supported in 3D");

  const vertex_type& vmin = m_vmin;
  const vertex_type& vmax = m_vmax;

  auto write = [&](const vertex_type& a,
                   const vertex_type& b,
                   const vertex_type& c,
                   const vertex_type& d) {
    helper.face(std::vector<vertex_type>({trf * a, trf * b, trf * c, trf * d}),
                color);
  };

  write({vmin.x(), vmin.y(), vmin.z()},
        {vmin.x(), vmax.y(), vmin.z()},
        {vmin.x(), vmax.y(), vmax.z()},
        {vmin.x(), vmin.y(), vmax.z()});

  write({vmax.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmax.z()},
        {vmax.x(), vmin.y(), vmax.z()});

  write({vmin.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmin.y(), vmax.z()},
        {vmin.x(), vmin.y(), vmax.z()});

  write({vmin.x(), vmax.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmax.z()},
        {vmin.x(), vmax.y(), vmax.z()});

  write({vmin.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmin.z()},
        {vmin.x(), vmax.y(), vmin.z()});

  write({vmin.x(), vmin.y(), vmax.z()},
        {vmax.x(), vmin.y(), vmax.z()},
        {vmax.x(), vmax.y(), vmax.z()},
        {vmin.x(), vmax.y(), vmax.z()});
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 2, int>>
std::ostream&
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::svg(
    std::ostream& os,
    value_type    w,
    value_type    h,
    value_type    unit,
    std::string   label,
    std::string   fillcolor) const
{
  static_assert(DIM == 2, "SVG is only supported in 2D");

  vertex_type mid(w / 2., h / 2.);

  using transform_t = Eigen::Transform<value_t, DIM, Eigen::Affine>;

  transform_t trf = transform_t::Identity();
  trf.translate(mid);
  trf = trf * Eigen::Scaling(vertex_type(1, -1));
  trf.scale(unit);

  // auto draw_line = [&](const vertex_type& left_,
  // const vertex_type& right_,
  // std::string color,
  // size_t width) {

  // vertex_type left = trf*left_;
  // vertex_type right = trf*right_;
  // os << "<line ";

  // os << "x1=\"" << left.x() << "\" ";
  // os << "y1=\"" << left.y() << "\" ";
  // os << "x2=\"" << right.x() << "\" ";
  // os << "y2=\"" << right.y() << "\" ";

  // os <<" stroke=\"" << color << "\" stroke-width=\"" << width << "\"/>\n";

  //};

  auto draw_point = [&](const vertex_type& p_, std::string color, size_t r) {
    vertex_type p = trf * p_;
    os << "<circle ";
    os << "cx=\"" << p.x() << "\" cy=\"" << p.y() << "\" r=\"" << r << "\"";
    os << " fill=\"" << color << "\"";
    os << "/>\n";
  };

  auto draw_rect = [&](
      const vertex_type& center_, const vertex_type& size_, std::string color) {
    vertex_type size   = size_ * unit;
    vertex_type center = trf * center_ - size * 0.5;

    os << "<rect ";
    os << "x=\"" << center.x() << "\" y=\"" << center.y() << "\" ";
    os << "width=\"" << size.x() << "\" height=\"" << size.y() << "\"";
    os << " fill=\"" << color << "\"";
    os << "/>\n";
  };

  auto draw_text = [&](const vertex_type& center_,
                       std::string        text,
                       std::string        color,
                       size_t             size) {
    vertex_type center = trf * center_;
    os << "<text dominant-baseline=\"middle\" text-anchor=\"middle\" ";
    os << "fill=\"" << color << "\" font-size=\"" << size << "\" ";
    os << "x=\"" << center.x() << "\" y=\"" << center.y() << "\">";
    os << text << "</text>\n";
  };

  // std::array<vertex_type, 4> points = {
  // m_vmin,bb
  //{m_vmin.x(), m_vmax.y()},
  // m_vmax,
  //{m_vmax.x(), m_vmin.y()}
  //};

  draw_rect(m_center, m_width, fillcolor);
  draw_point(m_vmin, "black", 2);
  draw_point(m_vmax, "black", 2);
  draw_text(m_center, label, "white", 10);

  // for(size_t i=0;i<points.size();i++) {
  // size_t j = (i+1)%points.size();
  // draw_point(points.at(i), "black", 3);
  ////draw_line(points.at(i), points.at(j), "black", 3);
  //}

  return os;
}

template <typename box_t>
box_t*
octree_inner(std::vector<std::unique_ptr<box_t>>& store,
             size_t                               max_depth,
             typename box_t::vertex_array_type    envelope,
             const std::vector<box_t*>&           lprims,
             size_t                               depth)
{

  using vertex_type = typename box_t::vertex_type;
  // using vertex_array_type = typename box_t::vertex_array_type;

  assert(lprims.size() > 0);
  if (lprims.size() == 1) {
    // just return
    return lprims.front();
  }

  if (depth >= max_depth) {
    // just wrap them all up
    auto bb = std::make_unique<box_t>(lprims, envelope);
    store.push_back(std::move(bb));
    return store.back().get();
  }

  std::array<std::vector<box_t*>, 8> octants;
  // calc center of boxes
  vertex_type vmin, vmax;
  std::tie(vmin, vmax) = box_t::wrap(lprims);
  vertex_type glob_ctr = (vmin + vmax) / 2.;

  for (auto* box : lprims) {
    vertex_type ctr = box->center() - glob_ctr;
    if (ctr.x() < 0 && ctr.y() < 0 && ctr.z() < 0) {
      octants[0].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() < 0 && ctr.z() < 0) {
      octants[1].push_back(box);
      continue;
    }
    if (ctr.x() < 0 && ctr.y() > 0 && ctr.z() < 0) {
      octants[2].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() > 0 && ctr.z() < 0) {
      octants[3].push_back(box);
      continue;
    }

    if (ctr.x() < 0 && ctr.y() < 0 && ctr.z() > 0) {
      octants[4].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() < 0 && ctr.z() > 0) {
      octants[5].push_back(box);
      continue;
    }
    if (ctr.x() < 0 && ctr.y() > 0 && ctr.z() > 0) {
      octants[6].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() > 0 && ctr.z() > 0) {
      octants[7].push_back(box);
      continue;
    }

    // not in any quadrant (numerics probably)
    octants[0].push_back(box);
  }

  std::vector<box_t*> sub_octs;
  for (const auto& sub_prims : octants) {
    if (sub_prims.size() <= 8) {
      if (sub_prims.size() < 1) {
        // done
      } else if (sub_prims.size() == 1) {
        sub_octs.push_back(sub_prims.front());
      } else {
        store.push_back(std::make_unique<box_t>(sub_prims, envelope));
        sub_octs.push_back(store.back().get());
      }
    } else {
      // recurse
      sub_octs.push_back(
          octree_inner(store, max_depth, envelope, sub_prims, depth + 1));
    }
  }

  if (sub_octs.size() == 1) {
    return sub_octs.front();
  }

  // std::cout << "sub_octs.size() = " << sub_octs.size() << std::endl;
  auto bb = std::make_unique<box_t>(sub_octs, envelope);
  store.push_back(std::move(bb));
  return store.back().get();
}

template <typename box_t>
box_t*
Acts::make_octree(std::vector<std::unique_ptr<box_t>>& store,
                  const std::vector<box_t*>&           prims,
                  size_t                               max_depth,
                  typename box_t::value_type           envelope1)
{
  static_assert(box_t::dim == 3, "Octree can only be created in 3D");

  // using vertex_type       = typename box_t::vertex_type;
  using vertex_array_type = typename box_t::vertex_array_type;

  vertex_array_type envelope(vertex_array_type::Constant(envelope1));

  // std::function<box_t*(const std::vector<box_t*>&, size_t)> oct;

  /*oct = [&store, &max_depth, &oct, &envelope](const std::vector<box_t*>&
  lprims,
                                              size_t depth) -> box_t* {


  };*/

  box_t* top = octree_inner(store, max_depth, envelope, prims, 0);
  return top;
}

template <typename T, typename U, size_t V>
std::ostream&
operator<<(std::ostream& os, const Acts::AxisAlignedBoundingBox<T, U, V>& box)
{
  box.dump(os);
  return os;
}
