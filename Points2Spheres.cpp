#include "Points2Spheres.hpp"

#include "GeometryUtils.hpp"

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <math.h>

typedef std::unordered_map<Point3, size_t, boost::hash<Point3>> pt2ind;

size_t ipow(int base, int exp)
{
    size_t x = 1;
    for (int i = 0; i < exp; ++i) { x *= base; }
    return x;
}

template <class HDS>
class Build_sphere : public CGAL::Modifier_base<HDS>
{
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    const Point3& pt;
    const double r;
    const size_t n;
public:
    Build_sphere(const Point3& pt, const double r, const size_t n) : pt(pt), r(r), n(n) {}
    void operator()(HDS& hds)
    {
        Builder B(hds, true);
        B.begin_surface(2+2*ipow(4, n), ipow(4, n+1), 6*ipow(4, n));

        // Build a refined tetrahedron
        Point3 a(0, 0, -1), b(0, 0.9428, 1.0/3), c(-0.8165, -0.4714, 1.0/3), d(0.8165, -0.4714, 1.0/3);
        pt2ind map;
        _us_recurse(a, b, c, B, map, n);
        _us_recurse(d, c, b, B, map, n);
        _us_recurse(a, d, b, B, map, n);
        _us_recurse(a, c, d, B, map, n);
        
        // Make that tetrahedron a unit sphere (the noramlization) and then scale and translate it appropiately
        for (size_t i = 0; i < map.size(); ++i)
        {
            B.vertex(i)->point() = pt + r*normalized(B.vertex(i)->point() - CGAL::ORIGIN);
        }

        B.end_surface();
    }
    static void _us_recurse(const Point3& a, const Point3& b, const Point3& c,
                            Builder& B, pt2ind& map, int n)
    {
        if (n == 0)
        {
            // Base case: add a triangle
            size_t va = vert_ind(a, B, map);
            size_t vb = vert_ind(b, B, map);
            size_t vc = vert_ind(c, B, map);
            B.begin_facet();
            B.add_vertex_to_facet(va);
            B.add_vertex_to_facet(vb);
            B.add_vertex_to_facet(vc);
            B.end_facet();
        }
        else
        {
            // Bisect each side of the triangle
            Point3 ab((a.x() + b.x())/2, (a.y() + b.y())/2, (a.z() + b.z())/2);
            Point3 ac((a.x() + c.x())/2, (a.y() + c.y())/2, (a.z() + c.z())/2);
            Point3 bc((b.x() + c.x())/2, (b.y() + c.y())/2, (b.z() + c.z())/2);
            // Recurse into the four sub-triangles
            _us_recurse(a, ab, ac, B, map, n-1);
            _us_recurse(ab, b, bc, B, map, n-1);
            _us_recurse(bc, c, ac, B, map, n-1);
            _us_recurse(ab, bc, ac, B, map, n-1);
        }
    }
    static size_t vert_ind(const Point3& v, Builder& B, pt2ind& map)
    {
        pt2ind::iterator itr = map.find(v);
        if (itr != map.end()) { return itr->second; }
        size_t idx = map.size();
        B.add_vertex(v);
        map.insert(itr, std::make_pair(v, idx));
        return idx;
    }
};

Polyhedron3* point2sphere(const Point3& pt, double r, size_t n)
{
    Build_sphere<Polyhedron3::HalfedgeDS> bs(pt, r, n);
    Polyhedron3* P = new Polyhedron3();
    P->delegate(bs);
    return P;
}

std::vector<Polyhedron3*> points2spheres(const std::vector<Point3>& pts, double r, size_t n)
{
    std::vector<Polyhedron3*> spheres;
    spheres.reserve(pts.size());
    for (std::vector<Point3>::const_iterator itr = pts.begin(), end = pts.end(); itr != end; ++itr)
    {
        spheres.push_back(point2sphere(*itr, r, n));
    }
    return spheres;
}
