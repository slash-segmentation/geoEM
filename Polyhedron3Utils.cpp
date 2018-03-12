#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are utilities for Polyhedron3 objects. Many of the utilities are complex and take multiple
// lines.
///////////////////////////////////////////////////////////////////////////////////////////////////

Polygon2 facet_to_polygon2(Polyhedron3::Facet_const_handle f)
{
    const Polyhedron3::Halfedge_const_handle &a = f->halfedge(), &b = a->next(), &c = b->next();
    const Plane3 h = Plane3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point());
    std::vector<Point2> pts;
    FOR_VERTICES_AROUND_FACET(f, v) { pts.push_back(h.to_2d(v->point())); }
    return Polygon2(pts.begin(), pts.end());
}

// Calculate the volume of a polyhedron using Gauss's theorem / divergence theorem
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
        auto he = f->facet_begin();
        volume += CGAL::volume((he++)->vertex()->point(),
                               (he++)->vertex()->point(),
                               (he++)->vertex()->point(), center);
        assert(he == f->facet_begin()); // must be a triangular mesh
    }
    return volume;
}
bool is_not_degenerate(const Polyhedron3* P)
{
    // Check uniqueness of vertices
    std::unordered_set<Point3, boost::hash<Point3>> points;
    for (Polyhedron3::Point_const_iterator p = P->points_begin(), end = P->points_end(); p != end; ++p) { if (!points.insert(*p).second) { return false; } }

    // Check geometry of facets
    for (Polyhedron3::Facet_const_iterator f = P->facets_begin(), end = P->facets_end(); f != end; ++f)
    {
        const Polygon2& p = facet_to_polygon2(f);
        if (!p.is_simple() || p.area() == 0) { return false; }
    }

    return true;
}

// Point-in-Polyhedron Algorithm
// Based on Shulin et al. Point-in-polyhera test with direct handling of degeneracies. Geo-spatial Information Science. 2001.
struct pip_calculator
{
    typedef Polyhedron3::Facet_const_handle value_type;
    typedef const value_type& const_reference;
    const Point3 &p, &q;
    const Plane3 h;
    bool inside;
    pip_calculator(const Point3& p, const Point3& q, const Point3& vf) : p(p), q(q), h(p, q, vf), inside(false) { }
    void push_back(const Polyhedron3::Facet_const_handle& f)
    {
        const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next();
        const Point3 &A = b->vertex()->point(), &B = c->vertex()->point(), &C = a->vertex()->point();
        const bool AA = h.has_on_negative_side(A), AB = h.has_on_negative_side(B), AC = h.has_on_negative_side(C);
        inside = inside != ((((AA != AB) && (AA ? Plane3(q,B,A).has_on_negative_side(p) : Plane3(q,A,B).has_on_negative_side(p)))  !=
                             ((AB != AC) && (AB ? Plane3(q,C,B).has_on_negative_side(p) : Plane3(q,B,C).has_on_negative_side(p)))) !=
                             ((AC != AA) && (AC ? Plane3(q,A,C).has_on_negative_side(p) : Plane3(q,C,A).has_on_negative_side(p))));
    }
};
bool point_in_polyhedron(const Point3& p, const FacetTree& tree)
{
    const static Vector3 dir = random_nonnull_vector();
    const Point3 q = p + dir;
    static Point3 vf(1.0,1.0,1.0);
    while (CGAL::collinear(p, q, vf)) { vf = random_point(tree.bbox()); } // arbitrary point not collinear with p and q to define h
    pip_calculator pipc(p, q, vf);
    tree.all_intersected_primitives(Ray3(p, q), std::back_inserter(pipc));
    return pipc.inside;
}
