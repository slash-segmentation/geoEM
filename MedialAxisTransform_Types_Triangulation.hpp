#pragma once

#include "GeometryTypes.hpp"
#include "MedialAxisTransform_Types_MAT.hpp"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

#include <array>
#include <boost/logic/tribool.hpp>

//////////////////////////////////////////////////////////////////////////////////////////////
// Internal classes and typedefs for triangulation data during MAT construction.
//////////////////////////////////////////////////////////////////////////////////////////////

template <class Traits, class Vb = CGAL::Triangulation_vertex_base_3<Traits>>
class Triangulation_vertex_with_extras : public Vb
{
	typedef Vb Base;
	typedef Triangulation_vertex_with_extras<Traits> Self;
public:
	typedef typename Vb::Vertex_handle  Vertex_handle;
	typedef typename Vb::Cell_handle    Cell_handle;
	typedef typename Vb::Point          Point;
	template <class TDS2> struct Rebind_TDS { typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2; typedef Triangulation_vertex_with_extras<Traits, Vb2> Other; };
	Triangulation_vertex_with_extras() : Base(), original() { }
	Triangulation_vertex_with_extras(const Point& p) : Base(p), original() { }
	Triangulation_vertex_with_extras(const Point& p, Cell_handle cell) : Base(p, cell), original() { }
	Triangulation_vertex_with_extras(Cell_handle cell) : Base(cell), original() { }
	Polyhedron3::Vertex_handle original; // vertex in the original polyhedral mesh
};

template <class Traits, class Cb = CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<Traits>>
class Triangulation_cell_with_extras : public Cb
{
	// circumcenters are the Voronoi vertices (the point at which a maximum number of Voronoi regions touch)
	typedef Cb Base;
	typedef Triangulation_vertex_with_extras<Traits> Self;
public:
	typedef typename Cb::Vertex_handle  Vertex_handle;
	typedef typename Cb::Cell_handle    Cell_handle;
	template <class TDS2> struct Rebind_TDS { typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2; typedef Triangulation_cell_with_extras<Traits, Cb2> Other; };

	Triangulation_cell_with_extras() : Base(), inside(false), facet_intersects_mesh(), mat_vertex() { this->facet_intersects_mesh.fill(boost::logic::indeterminate); }
	Triangulation_cell_with_extras(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) : Base(v0, v1, v2, v3), inside(false), facet_intersects_mesh(), mat_vertex() { this->facet_intersects_mesh.fill(boost::logic::indeterminate); }
	Triangulation_cell_with_extras(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3) : Base(v0, v1, v2, v3, n0, n1, n2, n3), inside(false), facet_intersects_mesh(), mat_vertex() { this->facet_intersects_mesh.fill(boost::logic::indeterminate); }

	bool inside;
	std::array<boost::logic::tribool, 4> facet_intersects_mesh;
	MAT::Vertex_handle mat_vertex;
};

typedef Triangulation_vertex_with_extras<Kernel> TriVertex;
typedef Triangulation_cell_with_extras<Kernel> TriCell;
extern template class CGAL::Triangulation_data_structure_3<TriVertex, TriCell>;
typedef CGAL::Triangulation_data_structure_3<TriVertex, TriCell> TDS;
extern template class CGAL::Delaunay_triangulation_3<Kernel, TDS>;
typedef CGAL::Delaunay_triangulation_3<Kernel, TDS> Triangulation;

