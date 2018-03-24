#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>
#include <iterator>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
namespace PMP = CGAL::Polygon_mesh_processing;

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are utilities for Polyhedron3 objects. Many of the utilities are complex and take multiple
// lines.
///////////////////////////////////////////////////////////////////////////////////////////////////

Polygon2 facet_to_polygon2(Polyhedron3::Facet_const_handle f)
{
    // Assumes the facet is planar
    const Plane3 h = facet_to_plane3(f);
    std::vector<Point2> pts;
    FOR_VERTICES_AROUND_FACET(f, v) { pts.push_back(h.to_2d(v->point())); }
    return Polygon2(pts.begin(), pts.end());
}

FilteredPolyhedron3* filter_faces(Polyhedron3* P, const P3FacetSet& facets)
{
    return new FilteredPolyhedron3(*P, facets);
}
FilteredPolyhedron3* filter_vertices(Polyhedron3* P, const P3VertexSet& verts)
{
    P3FacetSet facets, processed;
    facets.reserve(verts.size()*2);
    processed.reserve(verts.size()*2);
    for (P3Vertex v : verts)
    {
        FOR_FACETS_AROUND_VERTEX(v, f)
        {
            if (processed.count(f)) { continue; } // facet already processed
            processed.insert(f); // now mark it as processed
            bool include = true;
            FOR_VERTICES_AROUND_FACET(f, vf)
            {
                if (!verts.count(vf)) { include = false; break; } // do not include facet - not all vertices are present
            }
            // Include it
            if (include) { facets.insert(f); }
        }
    }
    return filter_faces(P, facets);
}

// Checks uniqueness of vertex points and if each face is a simple polygon
bool is_not_degenerate(const Polyhedron3* P)
{
    // Check uniqueness of vertices
    std::unordered_set<Point3, boost::hash<Point3>> points;
    for (Polyhedron3::Point_const_iterator p = P->points_begin(), end = P->points_end(); p != end; ++p)
    {
        if (!points.insert(*p).second) { return false; }
    }

    // Check geometry of facets
    for (Polyhedron3::Facet_const_iterator f = P->facets_begin(), end = P->facets_end(); f != end; ++f)
    {
        const auto &a = f->halfedge(), &b = a->next(), &c = b->next();
        if (c->next() == a) // f->is_triangle()
        {
            // This is significantly more accurate and faster then the other method
            if (CGAL::collinear(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()))
            {
                return false;
            }
        }
        else
        {
            const Polygon2& p = facet_to_polygon2(f);
            if (!p.is_simple() || p.area() == 0)
            {
                return false;
            }
        }
    }

    return true;
}


// For triangulating holes (don't need the output)
struct null_output_iterator : std::iterator<std::output_iterator_tag, null_output_iterator>
{
    template<typename T> void operator=(T const&) { }
    null_output_iterator & operator++() { return *this; }
    null_output_iterator operator++(int) { return *this; }
    null_output_iterator & operator*() { return *this; }
};

// Triangulate all holes in the given polyhedron
void triangulate_holes(Polyhedron3* P)
{
    for (auto he = P->halfedges_begin(), end = P->halfedges_end(); he != end; ++he)
    {
        if (he->is_border())
        {
            PMP::triangulate_hole(*P, he, null_output_iterator());
        }
    }
    assert(P->is_closed());
    CGAL::set_halfedgeds_items_id(*P);
}
