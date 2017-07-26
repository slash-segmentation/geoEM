#include "Segments2Cylinders.hpp"

#include <GeometryUtils.hpp>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <math.h>

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
        B.begin_surface(2*n, 4*n-4+(n%2)*2);
		
		// Add vertices
		const Point3& p1 = seg.source();
		const Point3& p2 = seg.target();
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

		// Add triangles around cylinder
		for (size_t i = 0; i < n-1; ++i)
		{
			// clockwise... should it be counter-clockwise?
			B.begin_facet();
			B.add_vertex_to_facet(i+0);
			B.add_vertex_to_facet(i+1);
			B.add_vertex_to_facet(i+n);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(i+n);
			B.add_vertex_to_facet(i+1);
			B.add_vertex_to_facet(i+n+1);
			B.end_facet();
		}
		B.begin_facet();
		B.add_vertex_to_facet(n-1);
		B.add_vertex_to_facet(0);
		B.add_vertex_to_facet(2*n-1);
		B.end_facet();
		B.begin_facet();
		B.add_vertex_to_facet(0);
		B.add_vertex_to_facet(n);
		B.add_vertex_to_facet(2*n-1);
		B.end_facet();
		
		// Add top of cylinder (counter-clockwise)
		size_t n2 = (n + 1) / 2 - 1;
		B.begin_facet();
		B.add_vertex_to_facet(0);
		B.add_vertex_to_facet(n-1);
		B.add_vertex_to_facet(1);
		B.end_facet();
		for (size_t i = 1; i < n2; ++i)
		{
			B.begin_facet();
			B.add_vertex_to_facet(i+0);
			B.add_vertex_to_facet(n-i);
			B.add_vertex_to_facet(i+1);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(i+1);
			B.add_vertex_to_facet(n-i);
			B.add_vertex_to_facet(n-i-1);
			B.end_facet();
		}
		if ((n % 2) == 0)
		{
			B.begin_facet();
			B.add_vertex_to_facet(n2+0);
			B.add_vertex_to_facet(n2+2);
			B.add_vertex_to_facet(n2+1);
			B.end_facet();
		}		

		// Add bottom of cylinder (clockwise)
		B.begin_facet();
		B.add_vertex_to_facet(n+0);
		B.add_vertex_to_facet(n+1);
		B.add_vertex_to_facet(n+n-1);
		B.end_facet();
		for (size_t i = 1; i < n2; ++i)
		{
			B.begin_facet();
			B.add_vertex_to_facet(n+i+0);
			B.add_vertex_to_facet(n+i+1);
			B.add_vertex_to_facet(n+n-i);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(n+i+1);
			B.add_vertex_to_facet(n+n-i-1);
			B.add_vertex_to_facet(n+n-i);
			B.end_facet();
		}
		if ((n % 2) == 0)
		{
			B.begin_facet();
			B.add_vertex_to_facet(n+n2+0);
			B.add_vertex_to_facet(n+n2+1);
			B.add_vertex_to_facet(n+n2+2);
			B.end_facet();
		}		
		
		// Done
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
