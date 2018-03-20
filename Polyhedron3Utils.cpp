#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>

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

// Calculate the volume of a polyhedron using Gauss's theorem / divergence theorem
// Polyhedron must be pure triangular
Kernel::FT volume(const Polyhedron3* P)
{
    // To increase numerical stability make the 4th point in every tetrahedron the centriod of the
    // polyhedron instead of some random point like the origin.
    Kernel::FT x = 0, y = 0, z = 0;
    for (auto v = P->vertices_begin(), v_end = P->vertices_end(); v != v_end; ++v)
    {
        x += v->point().x(); y += v->point().y(); z += v->point().z();
    }
    Point3 center(x / Kernel::FT(P->size_of_vertices()),
                  y / Kernel::FT(P->size_of_vertices()),
                  z / Kernel::FT(P->size_of_vertices()));
    Kernel::FT volume = 0;
    for (auto f = P->facets_begin(), f_end = P->facets_end(); f != f_end; ++f)
    {
        const auto &a = f->halfedge(), &b = a->next(), &c = b->next();
        volume += CGAL::volume(center, a->vertex()->point(), b->vertex()->point(), c->vertex()->point());
        assert(a == c->next()); // must be a triangular mesh
    }
    return volume;
}

// Calculate the surface area of a single facet
inline Kernel::FT surface_area(P3CFacet f)
{
    const auto &a = f->halfedge(), &b = a->next(), &c = b->next();
    if (c->next() == a) // f->is_triangle()
    {
        return ft_sqrt(CGAL::squared_area(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()));
    }
    else
    {
        return CGAL::abs(facet_to_polygon2(f).area());
    }
}

// Calculate the surface area of a polyhedron by summing up the areas of each facet
Kernel::FT surface_area(const Polyhedron3* P)
{
    Kernel::FT sa = 0;
    for (auto f = P->facets_begin(), f_end = P->facets_end(); f != f_end; ++f)
    {
        sa += surface_area(f);
    }
    return sa;
}

// Calculate the surface area of a polyhedron covered by the given vertices (where every facet that
// has all vertices in verts in counted).
Kernel::FT surface_area_covered(const Polyhedron3* P, P3CVertexSet verts)
{
    Kernel::FT sa = 0;
    P3CFacetSet processed;
    for (P3CVertex v : verts)
    {
        FOR_FACETS_AROUND_VERTEX(v, f)
        {
            if (processed.count(f)) { continue; } // facet already processed
            processed.insert(f); // now mark it as processed
            bool count = true;
            FOR_VERTICES_AROUND_FACET(f, vf)
            {
                if (!verts.count(vf)) { count = false; break; } // do not count facet - not all vertices are present
            }
            // Count it
            if (count) { sa += surface_area(f); }
        }
    }
    return sa;
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
