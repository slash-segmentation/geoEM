#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are basic geometry utility functions, most are simple one-liners that get used in many
// places.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "GeometryTypes.hpp"

#include <CGAL/Random.h>
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


// Fast min and max of 2 or 3 numbers
template <class FT> inline FT min2(const FT x, const FT y) { return x < y ? x : y; }
template <class FT> inline FT max2(const FT x, const FT y) { return x > y ? x : y; }
template <class FT> inline FT min3(const FT x, const FT y, const FT z) { return x < y ? (x < z ? x : z) : (y < z ? y : z); }
template <class FT> inline FT max3(const FT x, const FT y, const FT z) { return x > y ? (x > z ? x : z) : (y > z ? y : z); }


// Inter-convert vectors and directions (vector->direction is easy [just constructor] but provided for completion)
template <class Kernel> inline CGAL::Direction_2<Kernel> to_direction(const CGAL::Vector_2<Kernel> &v) { return CGAL::Direction_2<Kernel>(v); }
template <class Kernel> inline CGAL::Direction_3<Kernel> to_direction(const CGAL::Vector_3<Kernel> &v) { return CGAL::Direction_3<Kernel>(v); }
template <class Kernel> inline CGAL::Vector_2<Kernel> to_vector(const CGAL::Direction_2<Kernel> &d) { return CGAL::Vector_2<Kernel>(d.dx(),d.dy(),d.dz()); }
template <class Kernel> inline CGAL::Vector_3<Kernel> to_vector(const CGAL::Direction_3<Kernel> &d) { return CGAL::Vector_3<Kernel>(d.dx(),d.dy(),d.dz()); }


// Normalize vectors and directions
template <class Kernel> inline CGAL::Vector_2<Kernel> normalized(const CGAL::Vector_2<Kernel> &v) { typename Kernel::FT m = v.squared_length(); return m == 0 ? v : v / ft_sqrt(m); }
template <class Kernel> inline CGAL::Vector_3<Kernel> normalized(const CGAL::Vector_3<Kernel> &v) { typename Kernel::FT m = v.squared_length(); return m == 0 ? v : v / ft_sqrt(m); }
template <class Kernel> inline CGAL::Direction_2<Kernel> normalized(const CGAL::Direction_2<Kernel> &d) { return CGAL::Direction_2<Kernel>(normalized(to_vector(d))); }
template <class Kernel> inline CGAL::Direction_3<Kernel> normalized(const CGAL::Direction_3<Kernel> &d) { return CGAL::Direction_3<Kernel>(normalized(to_vector(d))); }


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


// Get vertex and facet indices from Polyhedron3 and Polygon3 data structures.
// When marked DEPRECATED they should only be used in debugging since they are really, really slow (using list iterators).
//template <class I> DEPRECATED(inline typename std::iterator_traits<I>::difference_type _dist(I a, I b, std::input_iterator_tag)) { return std::distance(a, b); }
template <class I> inline typename std::iterator_traits<I>::difference_type _dist(I a, I b, std::input_iterator_tag) { return std::distance(a, b); }
template <class I> inline typename std::iterator_traits<I>::difference_type _dist(I a, I b, std::random_access_iterator_tag) { return b - a; }
template <class DS> inline size_t vertex_number  (const typename DS::Vertex_const_handle   v, const DS* ds) { return _dist(ds->vertices_begin(),  v, typename std::iterator_traits<typename DS::Vertex_const_handle>  ::iterator_category()); }
template <class DS> inline size_t facet_number   (const typename DS::Facet_const_handle    f, const DS* ds) { return _dist(ds->facets_begin(),    f, typename std::iterator_traits<typename DS::Facet_const_handle>   ::iterator_category()); }
template <class DS> inline size_t halfedge_number(const typename DS::Halfedge_const_handle e, const DS* ds) { return _dist(ds->halfedges_begin(), e, typename std::iterator_traits<typename DS::Halfedge_const_handle>::iterator_category()); }
// If you need to use these functions for list iterators for all elements, the below class is much better.
// It is constructed using the beginning iterator and the length and defines the [] operator to take an iterator and give an index.
// For random-access iterators it is equivalent to the above functions (so does not do lookups).
namespace irl_detail
{
    template <class I, class IT>
    struct Internal
    {
        handle_map<I, size_t> x;
        Internal(I begin, const size_t len) : x(size_t(len*1.2)) { this->x.reserve(len); for (size_t i = 0; i < len; ++begin, ++i) { this->x[begin] = i; } }
        size_t get(const I iter) const { return this->x.at(iter); }
    };
    template <class I>
    struct Internal<I, std::random_access_iterator_tag>
    {
        const I begin;
        Internal(I begin, const size_t len) : begin(begin) { }
        size_t get(const I iter) const { return iter - this->begin; }
    };
};
template <class I>
class IteratorReverseLookup
{
    irl_detail::Internal<I, typename std::iterator_traits<I>::iterator_category> x;
public:
    IteratorReverseLookup(I begin, const size_t len) : x(begin, len) { }
    size_t operator[](const I iter) const { return this->x.get(iter); }
};


// Get a string version of a point / vector identical to how the original CurveSk program did it (for debugging)
inline std::string tostr(const Point3& p)  { std::stringstream s; s << p; return "( " + s.str() + " )"; }
inline std::string tostr(const Vector3& v) { std::stringstream s; s << v; return "( " + s.str() + " )"; }


// Cartesian converters
//extern CGAL::Cartesian_converter<Kernel, EPEC_Kernel> main_to_exact;
//extern CGAL::Cartesian_converter<Kernel, EPIC_Kernel> main_to_inexact;
//extern CGAL::Cartesian_converter<EPEC_Kernel, Kernel> exact_to_main;
//extern CGAL::Cartesian_converter<EPIC_Kernel, Kernel> inexact_to_main;
//extern CGAL::Cartesian_converter<EPEC_Kernel, EPIC_Kernel> exact_to_inexact;
//extern CGAL::Cartesian_converter<EPIC_Kernel, EPEC_Kernel> inexact_to_exact;


// Bbox3 creation from points (significantly faster than adding the bboxes of each point)
// There are specializations for 1, 2, and 3 points, and a general function that takes an iterator of points
inline Bbox3 bbox3(const Point3& p) { return p.bbox(); }
#define RETURN_BBOX3 return Bbox3(min_x, min_y, min_z, max_x, max_y, max_z)
#define PT3_TO_DBL(P) P##x = CGAL::to_double(P.x()), P##y = CGAL::to_double(P.y()), P##z = CGAL::to_double(P.z())
#define MIN_MAX_2(d) if (p##d > q##d) { min_##d = q##d; max_##d = p##d; } else { min_##d = p##d; max_##d = q##d; }
inline Bbox3 bbox3(const Point3& p, const Point3& q) { double PT3_TO_DBL(p), PT3_TO_DBL(q), min_x, max_x, min_y, max_y, min_z, max_z; MIN_MAX_2(x); MIN_MAX_2(y); MIN_MAX_2(z); RETURN_BBOX3; }
#undef MIN_MAX_2
#define MIN_MAX_3(d) if (p##d < q##d) { if (p##d < r##d) { min_##d = p##d; max_##d = q##d < r##d ? r##d : q##d; } else { min_##d = r##d; max_##d = q##d; } } else if (p##d < r##d) { min_##d = q##d; max_##d = r##d; } else { min_##d = q##d < r##d ? q##d : r##d; max_##d = p##d; }
inline Bbox3 bbox3(const Point3& p, const Point3& q, const Point3& r) { double PT3_TO_DBL(p), PT3_TO_DBL(q), PT3_TO_DBL(r), min_x, max_x, min_y, max_y, min_z, max_z; MIN_MAX_3(x); MIN_MAX_3(y); MIN_MAX_3(z); RETURN_BBOX3; }
#undef MIN_MAX_3
#undef PT3_TO_DBL
#define PT3_TO_DBL(d) min_##d = CGAL::to_double(i->d()), max_##d = min_##d
#define MIN_MAX_2(d) double d = CGAL::to_double(i->d()); if (d < min_##d) { min_##d = d; } else if (d > max_##d) { max_##d = d; }
template <class InputIterator> inline Bbox3 bbox3_from_points(InputIterator i, InputIterator end) { if (i == end) { return Bbox3(); } double PT3_TO_DBL(x), PT3_TO_DBL(y), PT3_TO_DBL(z); while (++i != end) { MIN_MAX_2(x); MIN_MAX_2(y); MIN_MAX_2(z); } RETURN_BBOX3; }
#undef MIN_MAX_2
#undef PT3_TO_DBL
#undef RETURN_BBOX3


// Bbox Utilities
inline bool is_inside(const Bbox2& out, const Bbox2& in) { return out.xmin() <= in.xmin() && out.ymin() <= in.ymin() && out.xmax() >= in.xmax() && out.ymax() >= in.ymax(); }
inline bool is_inside(const Bbox3& out, const Bbox3& in) { return out.xmin() <= in.xmin() && out.ymin() <= in.ymin() && out.zmin() <= in.zmin() && out.xmax() >= in.xmax() && out.ymax() >= in.ymax() && out.zmax() >= in.zmax(); }
inline bool is_inside(const Bbox2& out, const Point2& p) { return p.x() >= out.xmin() && p.x() <= out.xmax() && p.y() >= out.ymin() && p.y() <= out.ymax(); }
inline bool is_inside(const Bbox3& out, const Point3& p) { return p.x() >= out.xmin() && p.x() <= out.xmax() && p.y() >= out.ymin() && p.y() <= out.ymax() && p.z() >= out.zmin() && p.z() <= out.zmax(); }
inline double bbox_max_length(const Bbox3& bb) { return max3(bb.max(0)-bb.min(0), bb.max(1)-bb.min(1), bb.max(2)-bb.min(2)); }


// Random constructs
extern CGAL::Random Rand;
inline Point3 random_point(const Bbox3& bbox) { return Point3(Rand.uniform_real(bbox.xmin(),bbox.xmax()),Rand.uniform_real(bbox.ymin(),bbox.ymax()),Rand.uniform_real(bbox.zmin(),bbox.zmax())); }
inline Vector3 random_vector() { return Vector3(Rand.uniform_real(0.0,1.0),Rand.uniform_real(0.0,1.0),Rand.uniform_real(0.0,1.0)); }
inline Vector3 random_nonnull_vector() { Vector3 v = random_vector(); while (v == CGAL::NULL_VECTOR) { v = random_vector(); } return v; }


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
