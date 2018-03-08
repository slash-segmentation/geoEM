#include "Segments2Cylinders.hpp"

#include "GeometryUtils.hpp"

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <math.h>


inline static size_t n_cyl_verts(const size_t n) { return 2*n; }
inline static size_t n_cyl_faces(const size_t n) { return 2*n; }
inline static size_t n_cyl_end_faces(const size_t n) { return n-2+(n%2); }

template <class HDS>
static void add_cyl_vertices(CGAL::Polyhedron_incremental_builder_3<HDS>& B, const Segment3& seg, const double r, const size_t n)
{
    const Point3& p1 = seg.source(), p2 = seg.target();
    const Plane3 h(p1, p2-p1);
    const Vector3 a = r * normalized(h.base1()), b = r * normalized(h.base2());
    // Top circle vertices
    for (size_t i = 0; i < n; ++i)
    {
        double theta = i * 2 * M_PI / n;
        B.add_vertex(p1 + cos(theta) * a + sin(theta) * b);
    }
    // Bottom circle vertices
    for (size_t i = 0; i < n; ++i)
    {
        double theta = i * 2 * M_PI / n;
        B.add_vertex(p2 + cos(theta) * a + sin(theta) * b);
    }
}
template <class HDS>
static void add_cyl_faces(CGAL::Polyhedron_incremental_builder_3<HDS>& B, size_t n, size_t off=0)
{
    off *= n_cyl_verts(n);
    
    // Add triangles around cylinder
    for (size_t i = 0; i < n-1; ++i)
    {
        // clockwise... should it be counter-clockwise?
        B.begin_facet();
        B.add_vertex_to_facet(off+i+0);
        B.add_vertex_to_facet(off+i+1);
        B.add_vertex_to_facet(off+i+n);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(off+i+n);
        B.add_vertex_to_facet(off+i+1);
        B.add_vertex_to_facet(off+i+n+1);
        B.end_facet();
    }
    B.begin_facet();
    B.add_vertex_to_facet(off+n-1);
    B.add_vertex_to_facet(off+0);
    B.add_vertex_to_facet(off+2*n-1);
    B.end_facet();
    B.begin_facet();
    B.add_vertex_to_facet(off+0);
    B.add_vertex_to_facet(off+n);
    B.add_vertex_to_facet(off+2*n-1);
    B.end_facet();
}
template <class HDS>
static void add_cyl_end_face_ccw(CGAL::Polyhedron_incremental_builder_3<HDS>& B, size_t n, size_t off=0)
{
    off *= n_cyl_verts(n);
    const size_t n2 = (n + 1) / 2 - 1;
    
    // Add top of cylinder (counter-clockwise)
    B.begin_facet();
    B.add_vertex_to_facet(off+0);
    B.add_vertex_to_facet(off+n-1);
    B.add_vertex_to_facet(off+1);
    B.end_facet();
    for (size_t i = 1; i < n2; ++i)
    {
        B.begin_facet();
        B.add_vertex_to_facet(off+i+0);
        B.add_vertex_to_facet(off+n-i);
        B.add_vertex_to_facet(off+i+1);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(off+i+1);
        B.add_vertex_to_facet(off+n-i);
        B.add_vertex_to_facet(off+n-i-1);
        B.end_facet();
    }
    if ((n % 2) == 0)
    {
        B.begin_facet();
        B.add_vertex_to_facet(off+n2+0);
        B.add_vertex_to_facet(off+n2+2);
        B.add_vertex_to_facet(off+n2+1);
        B.end_facet();
    }
}
template <class HDS>
static void add_cyl_end_face_cw(CGAL::Polyhedron_incremental_builder_3<HDS>& B, size_t n, size_t off=0)
{
    off *= n_cyl_verts(n);
    const size_t n2 = (n + 1) / 2 - 1;
    
    // Add bottom of cylinder (clockwise)
    B.begin_facet();
    B.add_vertex_to_facet(off+n+0);
    B.add_vertex_to_facet(off+n+1);
    B.add_vertex_to_facet(off+n+n-1);
    B.end_facet();
    for (size_t i = 1; i < n2; ++i)
    {
        B.begin_facet();
        B.add_vertex_to_facet(off+n+i+0);
        B.add_vertex_to_facet(off+n+i+1);
        B.add_vertex_to_facet(off+n+n-i);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(off+n+i+1);
        B.add_vertex_to_facet(off+n+n-i-1);
        B.add_vertex_to_facet(off+n+n-i);
        B.end_facet();
    }
    if ((n % 2) == 0)
    {
        B.begin_facet();
        B.add_vertex_to_facet(off+n+n2+0);
        B.add_vertex_to_facet(off+n+n2+1);
        B.add_vertex_to_facet(off+n+n2+2);
        B.end_facet();
    }
}

template <class HDS>
class Build_cylinder : public CGAL::Modifier_base<HDS>
{
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    const Segment3& seg;
    const double r;
    const size_t n;
public:
    Build_cylinder(const Segment3& seg, const double r, const size_t n) : seg(seg), r(r), n(n) {}
    void operator()(HDS& hds)
    {
        Builder B(hds, true);
        B.begin_surface(n_cyl_verts(n), n_cyl_faces(n) + 2*n_cyl_end_faces(n));
        add_cyl_vertices(B, seg, r, n);
        add_cyl_faces(B, n);
        add_cyl_end_face_ccw(B, n);
        add_cyl_end_face_cw(B, n);
        B.end_surface();
    }
};

template <class HDS>
class Build_cylinder_merged : public CGAL::Modifier_base<HDS>
{
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    const std::vector<Segment3>& segs;
    const double r;
    const size_t n;
public:
    Build_cylinder_merged(const std::vector<Segment3>& segs, const double r, const size_t n) : segs(segs), r(r), n(n) {}
    void operator()(HDS& hds)
    {
        size_t nsegs = segs.size();
        Builder B(hds, true);
        B.begin_surface(nsegs*n_cyl_verts(n), nsegs*n_cyl_faces(n) + 2*n_cyl_end_faces(n));
        for (size_t i = 0; i < nsegs; ++i) { add_cyl_vertices(B, segs[i], r, n); }
        add_cyl_end_face_ccw(B, n, 0);
        for (size_t i = 0; i < nsegs; ++i) { add_cyl_faces(B, n, i); }
        add_cyl_end_face_cw(B, n, nsegs-1);
        B.end_surface();
    }
};

Polyhedron3* segment2cylinder(const Segment3& seg, double r, size_t n)
{
    Build_cylinder<Polyhedron3::HalfedgeDS> bc(seg, r, n);
    Polyhedron3* P = new Polyhedron3();
    P->delegate(bc);
    return P;
}

Polyhedron3* segments2cylinders_merged(const std::vector<Segment3>& segs, double r, size_t n)
{
    Build_cylinder_merged<Polyhedron3::HalfedgeDS> bcm(segs, r, n);
    Polyhedron3* P = new Polyhedron3();
    P->delegate(bcm);
    return P;
}

std::vector<Polyhedron3*> segments2cylinders(const std::vector<Segment3>& segs, double r, size_t n)
{
    std::vector<Polyhedron3*> cyls;
    cyls.reserve(segs.size());
    for (std::vector<Segment3>::const_iterator itr = segs.begin(), end = segs.end(); itr != end; ++itr)
    {
        cyls.push_back(segment2cylinder(*itr, r, n));
    }
    return cyls;
}
