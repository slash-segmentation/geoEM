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

// Calculate the surface area of a polyhedron by summing up the areas of each facet
Kernel::FT surface_area(const Polyhedron3* P)
{
    Kernel::FT sa = 0;
    for (auto f = P->facets_begin(), f_end = P->facets_end(); f != f_end; ++f)
    {
        //if (f->is_triangle())
        //{
        //    const auto &a = f->halfedge(), &b = a->next(), &c = b->next();
        //    sa += CGAL::squared_area(a->vertex()->point(), b->vertex()->point(), c->vertex()->point());
        //}
        //else
        //{
            sa += facet_to_polygon2(f).area();
        //}
    }
    return CGAL::abs(sa);
}

// Checks uniqueness of vertex points and if each face is a simple polygon
bool is_not_degenerate(const Polyhedron3* P)
{
    // Check uniqueness of vertices
    std::unordered_set<Point3, boost::hash<Point3>> points;
    for (Polyhedron3::Point_const_iterator p = P->points_begin(), end = P->points_end(); p != end; ++p) { if (!points.insert(*p).second) { return false; } }

    // Check geometry of facets
    for (Polyhedron3::Facet_const_iterator f = P->facets_begin(), end = P->facets_end(); f != end; ++f)
    {
        if (f->is_triangle())
        {
            // This is significantly more accurate and faster then the other method
            const auto &a = f->halfedge(), &b = a->next(), &c = b->next();
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
