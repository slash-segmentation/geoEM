#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are basic geometry utility functions, most are simple one-liners that get used in many
// places.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "GeometryTypes.hpp"

#include <CGAL/Handle_hash_function.h>

#include <list>
#include <unordered_map>
#include <unordered_set>


// Basic Utility Macros
#ifndef ARRAYSIZE
#define ARRAYSIZE(a) sizeof(a)/sizeof(a[0])
#endif
#define CONCATENATE_DETAIL(x, y) x##y
#define CONCATENATE(x, y) CONCATENATE_DETAIL(x, y)
#define MAKE_UNIQUE(x) CONCATENATE(x, __LINE__)


// Easy sqrt (uses sqrt if available for the given numeric type, otherwise uses the double approximation of sqrt)
template <class FT> inline double dbl_sqrt(const FT x) { return CGAL::sqrt(CGAL::to_double(x)); } // always a double-precision approximation of sqrt
template <class FT> inline FT _ft_sqrt(const FT x, CGAL::Integral_domain_without_division_tag) { return Kernel::FT(dbl_sqrt(x)); }
template <class FT> inline FT _ft_sqrt(const FT x, CGAL::Field_with_sqrt_tag) { return CGAL::sqrt(x); }
template <class FT> inline FT ft_sqrt(const FT x) { return _ft_sqrt(x, typename CGAL::Algebraic_structure_traits<FT>::Algebraic_category()); }


// Normalize vectors and directions
template <class Kernel> inline CGAL::Vector_2<Kernel> normalized(const CGAL::Vector_2<Kernel> &v) { typename Kernel::FT m = v.squared_length(); return m == 0 ? v : v / ft_sqrt(m); }
template <class Kernel> inline CGAL::Vector_3<Kernel> normalized(const CGAL::Vector_3<Kernel> &v) { typename Kernel::FT m = v.squared_length(); return m == 0 ? v : v / ft_sqrt(m); }


// Get the distance (not-squared) between two points (approximating to double always)
template <class Kernel> inline double distance(const CGAL::Point_2<Kernel>& p, const CGAL::Point_2<Kernel>& q) { return dbl_sqrt(CGAL::squared_distance(p, q)); }
template <class Kernel> inline double distance(const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q) { return dbl_sqrt(CGAL::squared_distance(p, q)); }


// Project 3D shapes into 2D using a plane
template <class Kernel> inline CGAL::Point_2<Kernel>    to_2d(const CGAL::Point_3<Kernel>& p,    const CGAL::Plane_3<Kernel>& h) { return h.to_2d(p); }
template <class Kernel> inline CGAL::Segment_2<Kernel>  to_2d(const CGAL::Segment_3<Kernel>& s,  const CGAL::Plane_3<Kernel>& h) { return CGAL::Segment_2<Kernel>(h.to_2d(s.source()), h.to_2d(s.target())); }
template <class Kernel> inline CGAL::Triangle_2<Kernel> to_2d(const CGAL::Triangle_3<Kernel>& t, const CGAL::Plane_3<Kernel>& h) { return CGAL::Triangle_2<Kernel>(h.to_2d(t[0]), h.to_2d(t[1]), h.to_2d(t[2])); }


// Removing (nearly) collinear points
// In all of these, cyclic means the first and last points are connected (and should be checked as well)
// Nearly collinear points are determined if the squared area of the triangle made by 3 consecutive points is less than the threshold
void remove_collinear_points(std::list<Point2>* pts, bool cyclic = true);
void remove_collinear_points(std::list<Point3>* pts, bool cyclic = true);
void remove_nearly_collinear_points(std::list<Point2>* pts, Kernel::FT threshold = Kernel::FT(0.5), bool cyclic = true);
void remove_nearly_collinear_points(std::list<Point3>* pts, Kernel::FT threshold = Kernel::FT(0.5), bool cyclic = true);


// Handle lookup collections and geometric hash functions for unordered_map/unordered_set
// Use handle_map and handle_set like "typedef handle_map<Handle_type, Value> my_handle_map;"
template <typename H, typename T> using handle_map = std::unordered_map<H, T, CGAL::Handle_hash_function>;
template <typename H>             using handle_set = std::unordered_set<H,    CGAL::Handle_hash_function>;
#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
namespace boost
#else
namespace CGAL
#endif
{
inline size_t hash_value(const EPEC_Kernel::FT x) { std::size_t h = 0; boost::hash_combine(h, CGAL::to_double(x)); return h; }
inline size_t hash_value(const Point2& p) { std::size_t h = 0; boost::hash_combine(h, p.x()); boost::hash_combine(h, p.y()); return h; }
inline size_t hash_value(const Point3& p) { std::size_t h = 0; boost::hash_combine(h, p.x()); boost::hash_combine(h, p.y()); boost::hash_combine(h, p.z()); return h; }
}


// Common Types and Sets for Skeleton and Polyhedron parts
typedef Skeleton3::vertex_descriptor S3VertexDesc;
typedef Skeleton3::edge_descriptor   S3EdgeDesc;
typedef Polyhedron3::Vertex_const_handle    P3CVertex;
typedef Polyhedron3::Halfedge_const_handle  P3CHalfedge;
typedef Polyhedron3::Facet_const_handle     P3CFacet;
typedef Polyhedron3::Vertex_handle    P3Vertex;
typedef Polyhedron3::Halfedge_handle  P3Halfedge;
typedef Polyhedron3::Facet_handle     P3Facet;
typedef handle_set<P3CVertex>   P3CVertexSet;
typedef handle_set<P3Vertex>    P3VertexSet;
typedef handle_set<P3CHalfedge> P3CHalfedgeSet;
typedef handle_set<P3Halfedge>  P3HalfedgeSet;
typedef handle_set<P3CFacet>    P3CFacetSet;
typedef handle_set<P3Facet>     P3FacetSet;


// Bbox Utilities
inline bool is_inside(const Bbox2& out, const Bbox2& in) { return out.xmin() <= in.xmin() && out.ymin() <= in.ymin() && out.xmax() >= in.xmax() && out.ymax() >= in.ymax(); }
inline bool is_inside(const Bbox3& out, const Bbox3& in) { return out.xmin() <= in.xmin() && out.ymin() <= in.ymin() && out.zmin() <= in.zmin() && out.xmax() >= in.xmax() && out.ymax() >= in.ymax() && out.zmax() >= in.zmax(); }
inline bool is_inside(const Bbox2& out, const Point2& p) { return p.x() >= out.xmin() && p.x() <= out.xmax() && p.y() >= out.ymin() && p.y() <= out.ymax(); }
inline bool is_inside(const Bbox3& out, const Point3& p) { return p.x() >= out.xmin() && p.x() <= out.xmax() && p.y() >= out.ymin() && p.y() <= out.ymax() && p.z() >= out.zmin() && p.z() <= out.zmax(); }


// Graph utilities
template <class Graph, class vd = typename Graph::vertex_descriptor, class ed = typename Graph::edge_descriptor>
inline vd opposite(const Graph& g, const ed e, const vd v)
{
    const vd& u = source(e, g);
    return (u != v) ? u : target(e, g);
}
template <class Graph, class vd = typename Graph::vertex_descriptor>
inline vd next_vertex(const Graph& g, const vd v)
{
    // Assumes that the degree of v is 1
    auto e = *out_edges(v, g).first; // the first outgoing edge
    return opposite(g, e, v);
}

template <class Graph, class vd = typename Graph::vertex_descriptor>
inline vd next_vertex(const Graph& g, const vd u, const vd v)
{
    // Gets the vertex after v coming from u
    // Assumes that the degree of v is 2
    BOOST_FOREACH(auto e, out_edges(v, g))
    {
        auto w = opposite(g, e, v);
        if (w != u) { return w; }
    }
    return u; // fallback...
}
