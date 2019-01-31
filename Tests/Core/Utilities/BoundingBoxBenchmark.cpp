#include <iostream>
#include <map>
#include <memory>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>
#include <set>

#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"


namespace Acts {
namespace Test {

struct Object {};

using ObjectBBox = Acts::AABB<Object>;

void
make_grid(std::vector<std::unique_ptr<ObjectBBox>>& boxes,
          ObjectBBox::vertex_type             vmin,
          ObjectBBox::vertex_type             vmax,
          ObjectBBox::vertex_array_type       n,
          ObjectBBox::vertex_type             hw)
{
  using vertex_type = ObjectBBox::vertex_type;
  Object      o;
  vertex_type step = (vmax - vmin).array() / (n - 1);
  // vertex_type hw = step / 2.;

  for (size_t i = 0; i < n[0]; i++) {
    for (size_t j = 0; j < n[1]; j++) {
      for (size_t k = 0; k < n[2]; k++) {
        vertex_type v(step.x() * i, step.y() * j, step.z() * k);

        vertex_type bmin = vmin + v - hw;
        vertex_type bmax = vmin + v + hw;
        //std::cout << bmin.x() << ", " << bmin.y() << ", " << bmin.z() <<
        //std::endl; std::cout << bmax.x() << ", " << bmax.y() << ", " <<
        //bmax.z() << std::endl; std::cout << std::endl;

        boxes.emplace_back(std::make_unique<ObjectBBox>(o, bmin, bmax));
      }
    }
  }
}


void test_points() {

  std::vector<std::unique_ptr<ObjectBBox>> boxstore;
  make_grid(boxstore, {-10, -10, -10}, {10, 10, 10}, {8, 8, 8}, {1, 1, 1});
  std::vector<ObjectBBox*> prim_boxes;
  prim_boxes.reserve(boxstore.size());
  std::transform(boxstore.begin(), boxstore.end(), std::back_inserter(prim_boxes),
      [](const auto& box) { return box.get(); });
  
  ObjectBBox* top;
  top = make_octree(boxstore, prim_boxes, 0);

  std::ofstream obj("octree.obj");
  size_t n_vtx = 1;
  
  const ObjectBBox* node = top;
  do {
    //std::cout << *node << std::endl;
    node->obj(obj, n_vtx);
    if(node->hasEntity()) {
      node = node->getSkip();
    }
    else {
      node = node->getLeftChild();
    }
  } while(node != nullptr);

  obj.close();
    
  auto dumb_search = [](const ObjectBBox::vertex_type& v, 
                        const std::vector<ObjectBBox*>& prims,
                        size_t& n_intersects,
                        const ObjectBBox*) -> const ObjectBBox*
  {
    for (const auto& bb : prims) {
      n_intersects++;
      if (bb->intersect(v)) { return bb; }
    }
    return nullptr;
  };

  auto bvh_search = [](const ObjectBBox::vertex_type& v, 
                        const std::vector<ObjectBBox*>& /*prims*/,
                        size_t& n_intersects,
                        const ObjectBBox* world) -> const ObjectBBox*
  {
    const ObjectBBox* lnode = world;
    do {
      n_intersects++;
      if (lnode->intersect(v)) {

        if (lnode->hasEntity()) {
          // found primitive
          return lnode;
        } else {
          // go over children
          lnode = lnode->getLeftChild();
        }
      } else {
        lnode = lnode->getSkip();
      }
    } while (lnode != nullptr);
    return nullptr;
  };

  using searcher = std::function<const ObjectBBox*(const ObjectBBox::vertex_type&, 
                                             const std::vector<ObjectBBox*>&,
                                             size_t&,
                                             const ObjectBBox*)>;
  
  auto bench = [](const std::vector<ObjectBBox::vertex_type>& points,
                  const std::vector<ObjectBBox*>& prims,
                  const ObjectBBox*               world,
                  searcher                        search)
      -> std::tuple<size_t, double, std::vector<const ObjectBBox*>> 
  {

    size_t n_intersects = 0;

    std::vector<const ObjectBBox*> hits;
    hits.reserve(points.size());

    std::chrono::steady_clock::time_point begin
        = std::chrono::steady_clock::now();

    for (const auto& v : points) { hits.push_back(search(v, prims, n_intersects, world)); }

    std::chrono::steady_clock::time_point end
        = std::chrono::steady_clock::now();
    double diff
        = std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
              .count()
        / 1000000.;



    return {n_intersects, diff, hits};
  };

  auto compare_hits = [](const std::vector<ObjectBBox::vertex_type>& points,
                         const std::vector<const ObjectBBox*>& a, 
                         const std::vector<const ObjectBBox*>& b)
  {
    if (a.size() != b.size()) {
      throw std::runtime_error("compare: nonequal size");
    }

    for (size_t i = 0; i < a.size(); i++) {
      auto hit_a = a.at(i);
      auto hit_b  = b.at(i);
      if(hit_a != hit_b) {
        const auto& p = points.at(i);
        std::cout << "p: " << p.x() << " " << p.y() << " " << p.z() << " -> ";
        std::cout << "a: " << hit_a << ", b: " << hit_b << std::endl;

        if(hit_a != nullptr) {
          std::cout << *hit_a << std::endl;
        }

        throw std::runtime_error("Unequal hits pointer");
      }
    }
  };

    
  ObjectBBox::vertex_array_type width(1, 1, 1);
  ObjectBBox::vertex_array_type space(0.2, 0.2, 0.2);

  std::mt19937 rng(42);

  std::vector<size_t> npoints = {/*10, 
                                 100, 
                                 1000, 
                                 10000, 
                                 100000,*/ 
                                 1000000
                                 /*10000000*/};
  size_t n_boxes_1d_max = 20;

  std::vector<size_t> octree_depths = {0, 1, 2, 3, 4, 5};

  std::ofstream csv_centers("centers.csv");
  csv_centers << "n_boxes,octree_depth,n_int_dumb,n_int_bvh,diff_dumb,diff_bvh" << std::endl;
  csv_centers << std::setprecision(6);
  
  std::ofstream csv_rand("random.csv");
  csv_rand << "n_points,n_boxes,octree_depth,n_int_dumb,n_int_bvh,diff_dumb,diff_bvh" << std::endl;
  csv_rand << std::setprecision(6);

  for(size_t n_boxes_1d=4;n_boxes_1d<n_boxes_1d_max;n_boxes_1d++) {
    if(n_boxes_1d%2 != 0) {
      continue;
    }

    for(size_t octree_depth : octree_depths) {
      // produce primitives for problem size
      std::vector<std::unique_ptr<ObjectBBox>> lboxstore;
      std::cout << "n_box_1D: " << n_boxes_1d << " oct_d: " << octree_depth << std::endl;

      ObjectBBox::vertex_array_type n(n_boxes_1d, n_boxes_1d, n_boxes_1d);
      ObjectBBox::vertex_type vmin = -1*n*(width+space);
      ObjectBBox::vertex_type vmax = +1*n*(width+space);

      //std::cout << vmin << "\n" << vmax << std::endl;
      make_grid(lboxstore, vmin, vmax, n, width);

      std::vector<ObjectBBox*> lprim_boxes;
      lprim_boxes.reserve(lboxstore.size());
      std::transform(lboxstore.begin(),lboxstore.end(), std::back_inserter(lprim_boxes),
          [](const auto& box) { return box.get(); });

      ObjectBBox* ltop;
      ltop = make_octree(lboxstore, lprim_boxes, octree_depth, 0.1);
      
      std::stringstream ss;
      ss << "obj/" << n_boxes_1d << "_" << octree_depth << ".obj";
      size_t ln_vtx = 1;
      std::ofstream lobj(ss.str());
      for(const auto& box : lboxstore) {
        box->obj(lobj, ln_vtx);
      }

      std::vector<ObjectBBox::vertex_type> centers;
      centers.reserve(lprim_boxes.size());
      for(const auto* prim: lprim_boxes) {
        centers.push_back(prim->center());
      }
      std::shuffle(centers.begin(), centers.end(), rng);

      size_t n_int_dumb, n_int_bvh;
      double diff_dumb, diff_bvh;
      std::vector<const ObjectBBox*> hits_dumb, hits_bvh;

      std::tie(n_int_dumb, diff_dumb, hits_dumb) = bench(centers, lprim_boxes, ltop, dumb_search);
      std::tie(n_int_bvh, diff_bvh, hits_bvh) = bench(centers, lprim_boxes, ltop, bvh_search);
      compare_hits(centers, hits_dumb, hits_bvh);

      csv_centers << n_boxes_1d << "," << octree_depth << "," << n_int_dumb << "," << n_int_bvh
                  << "," << diff_dumb << "," << diff_bvh << std::endl;;

      std::cout << "o centers: dumb/bvh: " << diff_dumb << "s " << diff_bvh << "s, " << n_int_dumb << ", " << n_int_bvh;
      std::cout << std::endl;
  
      std::uniform_real_distribution<ObjectBBox::value_type> dist(ltop->min().minCoeff()*0.9, ltop->max().maxCoeff()*1/0.9);

      for(size_t size : npoints) {
        std::vector<ObjectBBox::vertex_type> rand_points;
        rand_points.reserve(size);
        for(size_t i=0;i<size;i++) {
          ObjectBBox::vertex_type v(dist(rng), dist(rng), dist(rng));
          rand_points.push_back(std::move(v));
        }

        hits_dumb.resize(0);
        hits_bvh.resize(0);

        std::tie(n_int_dumb, diff_dumb, hits_dumb) = bench(rand_points, lprim_boxes, ltop, dumb_search);
        std::tie(n_int_bvh, diff_bvh, hits_bvh) = bench(rand_points, lprim_boxes, ltop, bvh_search);
        compare_hits(rand_points, hits_dumb, hits_bvh);
        
        csv_rand << size << "," << n_boxes_1d << "," << octree_depth << "," << n_int_dumb << "," << n_int_bvh
                    << "," << diff_dumb << "," << diff_bvh << std::endl;;
      
        std::cout << "o rand: size: " << size << " dumb/bvh: " << diff_dumb << "s " << diff_bvh
          << "s, " << n_int_dumb << ", " << n_int_bvh;
        std::cout << std::endl;
      

      }
      
    }

  }

  csv_centers.close();
  csv_rand.close();

  
}

void test_rays() 
{
  using Ray = Ray3F;
  using vertex_type = ObjectBBox::vertex_type;
  

  auto dumb_search = [](const Ray& ray, 
                        const std::vector<ObjectBBox*>& prims,
                        size_t& n_intersects,
                        const ObjectBBox*) -> std::vector<const ObjectBBox*>
  {
    std::vector<const ObjectBBox*> result;
    for (const auto& bb : prims) {
      n_intersects++;
      if (bb->intersect(ray)) { 
        result.push_back(bb);
      }
    }
    return result;
  };

  auto bvh_search = [](const Ray& ray, 
                        const std::vector<ObjectBBox*>& /*prims*/,
                        size_t& n_intersects,
                        const ObjectBBox* world) -> std::vector<const ObjectBBox*>
  {
    //std::cout << "go for " << ray << std::endl;
    std::vector<const ObjectBBox*> result;
    const ObjectBBox* lnode = world;
    do {
      n_intersects++;
      if (lnode->intersect(ray)) {

        if (lnode->hasEntity()) {
          // found primitive
          result.push_back(lnode);
          //std::cout << "intersect prim: " << lnode << " -> " << lnode->getSkip() << std::endl;
          lnode = lnode->getSkip();
        } else {
          // go over children
          //std::cout << "intersect bvh: " << lnode << " -> " << lnode->getLeftChild() << std::endl;
          lnode = lnode->getLeftChild();
        }
      } else {
        //std::cout << "no intersect: " << lnode << " -> " << lnode->getSkip() << std::endl;
        lnode = lnode->getSkip();
      }

      //if(n_intersects > 80) throw std::runtime_error("");
    } while (lnode != nullptr);
    //std::cout << "done" << std::endl;
    return result;
  };

  using searcher = std::function<std::vector<const ObjectBBox*>(const Ray&, 
                                             const std::vector<ObjectBBox*>&,
                                             size_t&,
                                             const ObjectBBox*)>;
  
  auto bench = [](const std::vector<Ray>& rays,
                  const std::vector<ObjectBBox*>& prims,
                  const ObjectBBox*               world,
                  searcher                        search)
      -> std::tuple<size_t, double, std::vector<std::vector<const ObjectBBox*>>> 
  {

    size_t n_intersects = 0;

    std::vector<std::vector<const ObjectBBox*>> hits;
    hits.reserve(rays.size());

    std::chrono::steady_clock::time_point begin
        = std::chrono::steady_clock::now();

    for(size_t i=0;i<rays.size();i++) {
      const auto& r = rays[i];
      hits.push_back(search(r, prims, n_intersects, world));
    }

    std::chrono::steady_clock::time_point end
        = std::chrono::steady_clock::now();
    double diff
        = std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
              .count()
        / 1000000.;



    return {n_intersects, diff, hits};
  };

  auto compare_hits = [](const std::vector<Ray>& rays,
                         const std::vector<std::vector<const ObjectBBox*>>& a, 
                         const std::vector<std::vector<const ObjectBBox*>>& b)
  {
    if (a.size() != b.size()) {
      throw std::runtime_error("compare: nonequal size");
    }

    for (size_t i = 0; i < a.size(); i++) {
      const auto& ray = rays.at(i);
      auto hits_a = a.at(i);
      auto hits_b  = b.at(i);

      if(hits_a.size() != hits_b.size()) {
        std::cout << ray << std::endl;
        throw std::runtime_error("Ray hit different number of primitives");
      }

      std::set<const ObjectBBox*> hits_set_a(hits_a.begin(), hits_a.end());
      std::set<const ObjectBBox*> hits_set_b(hits_b.begin(), hits_b.end());

      std::set<const ObjectBBox*> hits_diff;
      std::set_difference(hits_set_a.begin(), hits_set_a.end(), hits_set_b.begin(), hits_set_b.end(),
          std::inserter(hits_diff, hits_diff.begin()));

      if(hits_diff.size() != 0) {
        std::cout << ray << std::endl;
        for(const auto* hit : hits_diff) {
          std::cout << *hit << std::endl;
        }
        throw std::runtime_error("Ray has different set of hits");
      }

    }
  };
    
  ObjectBBox::vertex_array_type width(1, 1, 1);
  ObjectBBox::vertex_array_type space(0.2, 0.2, 0.2);

  std::mt19937 rng(43);

  //std::vector<size_t> nrays = {[>10, 
                                 //100, 
                                 //1000, 
                                 //10000, 
                                 //100000,*/ 
                                 //1000000
                                 //[>10000000<]};
  size_t nrays = 1e6;
  size_t n_boxes_1d_min = 22;
  size_t n_boxes_1d_max = 23;

  std::vector<size_t> octree_depths = {0, 1, 2, 3, 4, 5};

  std::ofstream csv_rand("rays_random.csv");
  csv_rand << "n_rays,n_boxes,octree_depth,n_int_dumb,n_int_bvh,diff_dumb,diff_bvh" << std::endl;
  csv_rand << std::setprecision(6);

  for(size_t n_boxes_1d=n_boxes_1d_min;n_boxes_1d<n_boxes_1d_max;n_boxes_1d++) {
    if(n_boxes_1d%2 != 0) {
      continue;
    }

    for(size_t octree_depth : octree_depths) {
      // produce primitives for problem size
      std::vector<std::unique_ptr<ObjectBBox>> lboxstore;
      std::cout << "n_box_1D: " << n_boxes_1d << " oct_d: " << octree_depth << std::flush;

      ObjectBBox::vertex_array_type n(n_boxes_1d, n_boxes_1d, n_boxes_1d);
      ObjectBBox::vertex_type vmin = -1*n*(width+space);
      ObjectBBox::vertex_type vmax = +1*n*(width+space);

      //std::cout << vmin << "\n" << vmax << std::endl;
      make_grid(lboxstore, vmin, vmax, n, width);

      std::vector<ObjectBBox*> lprim_boxes;
      lprim_boxes.reserve(lboxstore.size());
      //std::cout << " prims: " << std::endl;
      std::transform(lboxstore.begin(),lboxstore.end(), std::back_inserter(lprim_boxes),
          [](const auto& box) { 
            //std::cout << box.get() << std::endl;
            return box.get(); 
          });

      //std::cout << " -- " << std::endl;

      ObjectBBox* ltop;
      ltop = make_octree(lboxstore, lprim_boxes, octree_depth, 0.1);
  
      std::uniform_real_distribution<ObjectBBox::value_type> dist(ltop->min().minCoeff()*0.9, ltop->max().maxCoeff()*1/0.9);
      std::uniform_real_distribution<ObjectBBox::value_type> dirDist(-1, 1);

      std::vector<Ray> rays;
      rays.reserve(nrays);

      for(size_t i=0;i<nrays;i++) {
        vertex_type pos(dist(rng), dist(rng), dist(rng));
        vertex_type dir(dirDist(rng), dirDist(rng), dirDist(rng));
        dir.normalize();
        Ray ray(pos, dir);

        //std::cout << ray << std::endl;

        rays.push_back(ray);
      }

      size_t n_int_dumb, n_int_bvh;
      double diff_dumb, diff_bvh;
      std::vector<std::vector<const ObjectBBox*>> hits_dumb, hits_bvh;

      //std::tie(n_int_dumb, diff_dumb, hits_dumb) = bench(rays, lprim_boxes, ltop, dumb_search);
      std::tie(n_int_bvh, diff_bvh, hits_bvh) = bench(rays, lprim_boxes, ltop, bvh_search);
      //compare_hits(rays, hits_dumb, hits_bvh);
        
      csv_rand << nrays << "," << n_boxes_1d << "," << octree_depth << "," << n_int_dumb << "," << n_int_bvh
                        << "," << diff_dumb << "," << diff_bvh << std::endl;;
      std::cout << " rand: nrays: " << nrays << " dumb/bvh: " << diff_dumb << "s " << diff_bvh
        << "s, " << n_int_dumb << ", " << n_int_bvh;
      std::cout << " -> int/sec:" << (n_int_dumb / float(diff_dumb));
      std::cout << std::endl;

    }

  }

  csv_rand.close();
}

void test_intersect_points()
{
  using vertex_type = ObjectBBox::vertex_type;

  Object o;
  ObjectBBox bb(o, {0, 0, 0}, {1, 1, 1});
  vertex_type p; 

  p = {0.5, 0.5, 0.5};
  assert(bb.intersect(p));
  p = {0.25, 0.25, 0.25};
  assert(bb.intersect(p));
  p = {0.75, 0.75, 0.75};
  assert(bb.intersect(p));

  // lower bound is inclusive
  p = {0, 0, 0};
  assert(bb.intersect(p));
  // upper bound is exclusive
  p = {1.0, 1.0, 1.0};
  assert(!bb.intersect(p));

  // some outsides
  p = {2, 0, 0};
  assert(!bb.intersect(p));
  p = {0, 2, 0};
  assert(!bb.intersect(p));
  p = {0, 0, 2};
  assert(!bb.intersect(p));
  p = {2, 2, 0};
  assert(!bb.intersect(p));
  p = {2, 0, 2};
  assert(!bb.intersect(p));
  p = {2, 2, 2};
  assert(!bb.intersect(p));
  
  p = {-1, 0, 0};
  assert(!bb.intersect(p));
  p = {0, -1, 0};
  assert(!bb.intersect(p));
  p = {0, 0, -1};
  assert(!bb.intersect(p));
  p = {-1, -1, 0};
  assert(!bb.intersect(p));
  p = {-1, 0, -1};
  assert(!bb.intersect(p));
  p = {-1, -1, -1};
  assert(!bb.intersect(p));

}

void test_intersect_rays()
{
  using Box = AxisAlignedBoundingBox<Object, float, 2>;

  Object o;
  Box bb(o, {-1, -1}, {1, 1});

  // ray in positive x direction

  Ray<float, 2> ray({-2, 0}, {1, 0});
  assert(bb.intersect(ray));

  ray = {{-2, 2}, {1, 0}};
  assert(!bb.intersect(ray));

  ray = {{-2, -2}, {1, 0}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{-2, 1}, {1, 0}};
  assert(!bb.intersect(ray));

  // lower bound is inclusive
  ray = {{-2, -1}, {1, 0}};
  assert(bb.intersect(ray));

  // ray faces away from box
  ray = {{2, 0}, {1, 0}};
  assert(!bb.intersect(ray));


  // ray in negative x direction 

  ray = {{2, 0}, {-1, 0}};
  assert(bb.intersect(ray));

  ray = {{2, 2}, {-1, 0}};
  assert(!bb.intersect(ray));

  ray = {{2, -2}, {-1, 0}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{2, 1}, {-1, 0}};
  assert(!bb.intersect(ray));

  // lower bound is inclusive
  ray = {{2, -1}, {-1, 0}};
  assert(bb.intersect(ray));


  // ray in positive y direction 

  ray = {{0, -2}, {0, 1}};
  assert(bb.intersect(ray));

  ray = {{2, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  ray = {{-2, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{1, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  // lower bound is not inclusive, 
  // due to Eigen's NaN handling.
  ray = {{-1, -2}, {0, 1}};
  assert(!bb.intersect(ray));

  // other direction
  ray = {{0, -2}, {0, -1}};
  assert(!bb.intersect(ray));


  // ray in positive y direction 

  ray = {{0, 2}, {0, -1}};
  assert(bb.intersect(ray));

  ray = {{2, 2}, {0, -1}};
  assert(!bb.intersect(ray));

  ray = {{-2, 2}, {0, -1}};
  assert(!bb.intersect(ray));

  // upper bound is exclusive
  ray = {{1, 2}, {0, -1}};
  assert(!bb.intersect(ray));

  // lower bound is not inclusive, 
  // due to Eigen's NaN handling.
  ray = {{-1, 2}, {0, -1}};
  assert(!bb.intersect(ray));
  
  // other direction
  ray = {{0, 2}, {0, 1}};
  assert(!bb.intersect(ray));


  // some off axis rays

  ray = {{-2, 0}, {0.5, 0.25}};
  assert(bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, 0.4}};
  assert(bb.intersect(ray));

  ray = {{-2, 0}, {0.5, 0.6}};
  assert(!bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, 0.1}};
  assert(bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, -0.4}};
  assert(bb.intersect(ray));
  
  ray = {{-2, 0}, {0.5, -0.6}};
  assert(!bb.intersect(ray));
  
  ray = {{-2, 0}, {0.1, 0.5}};
  assert(!bb.intersect(ray));

  // starting point inside
  ray = {{0, 0,}, {-1, 0}};
  assert(bb.intersect(ray));
  ray = {{0, 0,}, {1, 0}};
  assert(bb.intersect(ray));
  ray = {{0, 0,}, {0, -1}};
  assert(bb.intersect(ray));
  ray = {{0, 0,}, {0, 1}};
  assert(bb.intersect(ray));



  using vertex_type3 = ObjectBBox::vertex_type;

  // lets make sure it also works in 3d
  ObjectBBox bb3(o, {-1, -1, -1}, {1, 1, 1});
  Ray<float, 3> ray3({0, 0, -2}, {0, 0, 1});
  assert(bb3.intersect(ray3));

  // facing away from box
  ray3 = {{0, 0, -2}, {0, 0, -1}};
  assert(!bb3.intersect(ray3));

  ray3 = {{0, 2, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));
  
  ray3 = {{0, -2, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  // right on slab
  ray3 = {{0, 1, -2}, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  // right on slab
  ray3 = {{0, -1, -2}, {0, 0, 1}};
  assert(bb3.intersect(ray3));

  // right on slab
  ray3 = {{-1, 0, -2 }, {0, 0, 1}};
  assert(!bb3.intersect(ray3));
  
  // right on slab
  ray3 = {{1, 0, -2 }, {0, 0, 1}};
  assert(!bb3.intersect(ray3));

  ray3 = {{-0.95, 0, -2 }, {0, 0, 1}};
  assert(bb3.intersect(ray3));
  
  // some off-axis rays
  ObjectBBox::vertex_type p(0, 0, -2);

  ray3 = {p, vertex_type3(1, 1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type3(-1, 1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type3(-1, -1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type3(1, -1, 1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type3(1.1, 0, -1)-p};
  assert(!bb3.intersect(ray3));
  
  ray3 = {p, vertex_type3(-1.1, 0, -1)-p};
  assert(!bb3.intersect(ray3));

  ray3 = {p, vertex_type3(0, 1.1, -1)-p};
  assert(!bb3.intersect(ray3));

  ray3 = {p, vertex_type3(0, -1.1, -1)-p};
  assert(!bb3.intersect(ray3));

  ray3 = {p, vertex_type3(0.9, 0, -1)-p};
  assert(bb3.intersect(ray3));
  
  ray3 = {p, vertex_type3(-0.9, 0, -1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type3(0, 0.9, -1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {p, vertex_type3(0, -0.9, -1)-p};
  assert(bb3.intersect(ray3));

  ray3 = {{0, 0, 0}, {1, 0, 0}};
  assert(bb3.intersect(ray3));
  ray3 = {{0, 0, 0}, {0, 1, 0}};
  assert(bb3.intersect(ray3));
  ray3 = {{0, 0, 0}, {0, 0, 1}};
  assert(bb3.intersect(ray3));
  
  ray3 = {{0, 0, 0}, {-1, 0, 0}};
  assert(bb3.intersect(ray3));
  ray3 = {{0, 0, 0}, {0, -1, 0}};
  assert(bb3.intersect(ray3));
  ray3 = {{0, 0, 0}, {0, 0, -1}};
  assert(bb3.intersect(ray3));

}


}
}


int main() {

  Acts::Test::test_intersect_points();
  Acts::Test::test_intersect_rays();

  //Acts::Test::test_points();
  Acts::Test::test_rays();
  

  return 0;
}
