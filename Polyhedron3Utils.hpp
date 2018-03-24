#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include <algorithm>

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are utilities for Polyhedron3 objects.
///////////////////////////////////////////////////////////////////////////////////////////////////

// Polyhedron part converters
Polygon2 facet_to_polygon2(Polyhedron3::Facet_const_handle f);
inline Plane3 facet_to_plane3(Polyhedron3::Facet_const_handle f) { const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next(); return Plane3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()); }
inline Triangle3 facet_to_triangle3(Polyhedron3::Facet_const_handle f) { const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next(); return Triangle3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()); }
inline Bbox3 facet_to_bbox3(Polyhedron3::Facet_const_handle f) { const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next(); return bbox3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()); }
inline Segment3 halfedge_to_segment3(Polyhedron3::Halfedge_const_handle he) { return Segment3(he->prev()->vertex()->point(), he->vertex()->point()); }


// Polyhedron part normals
inline Direction3 normal(const Polyhedron3::Facet_const_handle &f)
{
    // assumes triangle face
    Polyhedron3::Halfedge_around_facet_const_circulator he = f->facet_begin();
    Vector3 n = CGAL::normal(he->prev()->vertex()->point(), he->vertex()->point(), he->next()->vertex()->point());
    return Direction3(normalized(n));
}
inline Direction3 normal(const Polyhedron3::Vertex_const_handle &v)
{
    Vector3 n = CGAL::NULL_VECTOR;
    Polyhedron3::Halfedge_around_vertex_const_circulator he = v->vertex_begin(), end(he);
    CGAL_For_all(he, end) { if (!he->is_border()) { n = n + to_vector(normal(he->facet())); } }
    return Direction3(normalized(n));
}


// Calculate polyhedron facet planes
template <class Facet>
inline typename Facet::Plane_3 __facet_to_plane3(Facet& f)
{
    const Polyhedron3::Halfedge_const_handle &a = f.facet_begin(), &b = a->next(), &c = b->next(); 
    return typename Facet::Plane_3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point());
}
inline void calculate_facet_planes(Polyhedron3* P)
{
    std::transform(P->facets_begin(), P->facets_end(), P->planes_begin(), __facet_to_plane3<Polyhedron3::Facet>);
}


// Polyhedron looping around a vertex/facet
#define FOR_EDGES_AROUND_VERTEX(V, e) \
    auto MAKE_UNIQUE(_he) = V->vertex_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto e = MAKE_UNIQUE(_he); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), e = MAKE_UNIQUE(_he))
#define FOR_FACETS_AROUND_VERTEX(V, f) \
    auto MAKE_UNIQUE(_he) = V->vertex_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto f = MAKE_UNIQUE(_he)->facet(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), f = MAKE_UNIQUE(_he)->facet())
#define FOR_VERTICES_AROUND_VERTEX(V, v) \
    auto MAKE_UNIQUE(_he) = V->vertex_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto v = MAKE_UNIQUE(_he)->prev()->vertex(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), v = MAKE_UNIQUE(_he)->prev()->vertex())
#define FOR_EDGES_AROUND_FACET(F, e) \
    auto MAKE_UNIQUE(_he) = F->facet_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto e = MAKE_UNIQUE(_he); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), e = MAKE_UNIQUE(_he))
#define FOR_FACETS_AROUND_FACET(F, f) \
    auto MAKE_UNIQUE(_he) = F->facet_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto f = MAKE_UNIQUE(_he)->opposite()->facet(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), f = MAKE_UNIQUE(_he)->opposite()->facet())
#define FOR_VERTICES_AROUND_FACET(F, v) \
    auto MAKE_UNIQUE(_he) = F->facet_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto v = MAKE_UNIQUE(_he)->vertex(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), v = MAKE_UNIQUE(_he)->vertex())


// Other utilities
inline Kernel::FT avg_edge_length(const Polyhedron3 *P)
{
    Kernel::FT length = 0;
    for (auto e = P->edges_begin(), end = P->edges_end(); e != end; ++e)
    {
        length += distance(e->vertex()->point(), e->opposite()->vertex()->point());
    }
    return length / Kernel::FT(P->size_of_halfedges() / 2);
}
Kernel::FT surface_area_covered(const Polyhedron3* P, P3CVertexSet verts);
bool is_not_degenerate(const Polyhedron3* P);
void triangulate_holes(Polyhedron3* P);
