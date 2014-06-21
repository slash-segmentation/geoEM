#pragma once

#include "GeometryTypes.hpp"
#include "Polygons_3.hpp"
#include "MedialAxisTransform_Types_Geodesic.hpp"
#include "Flags.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: description
///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename P3>
class MAT_Polygons_3_Facet : public Polygons_3_Facet<P3>
{
	friend P3;
	typedef Polygons_3_Facet<P3> base;
public:
	FLAGS_DECL(Flags, StableSkeleton, NotStableSkeleton, Boundary, Eaten, InHeap);
	Flags flags;
private:
//#ifdef DETAILED_SKELETON
	Polyhedron3::Vertex_const_handle _vtchs, _vtche; // vertices of the dual Delaunay
//#endif
	double _weight;
	Vector3 _wgrad;
	//Vector3 _normal;
public:
	MAT_Polygons_3_Facet(const Polyhedron3::Vertex_const_handle V1, const Polyhedron3::Vertex_const_handle V2)
		: Polygons_3_Facet<P3>(), flags(0), _vtchs(V1), _vtche(V2) //, _normal(normalized(V1->point() - V2->point()))
	{
		const GeodesicDistance &gd_12 = V1->extra()->geod.at(V2), &gd_21 = V2->extra()->geod.at(V1);
		this->_weight = gd_12.distance;

		const Vector3 normal = normalized(V1->point() - V2->point());
		const Vector3 &wgrad_a = gd_21.derivative, &wgrad_b = gd_12.derivative;
		const Vector3 wgrad_a_proj = wgrad_a - (wgrad_a * /*this->_*/normal) * /*this->_*/normal; // normalized?
		const Vector3 wgrad_b_proj = wgrad_b - (wgrad_b * /*this->_*/normal) * /*this->_*/normal; // normalized?
		//std::cerr<<"dot(a, b): "<<wgrad_a_proj*wgrad_b_proj<<std::endl;
		this->_wgrad = wgrad_a_proj + wgrad_b_proj;
	}
	inline MAT_Polygons_3_Facet(const Polyhedron3::Vertex_const_handle V1, const Polyhedron3::Vertex_const_handle V2, double weight, const Vector3& wgrad) //, const Vector3& normal
		: Polygons_3_Facet<P3>(), flags(0),
//#ifdef DETAILED_SKELETON
		_vtchs(V1), _vtche(V2),
//#endif
		_weight(weight), _wgrad(wgrad) //, _normal(normal)
	{
	}
//#ifdef DETAILED_SKELETON
	inline Polyhedron3::Vertex_const_handle touching_vertex_start() const { return this->_vtchs; }
	inline Polyhedron3::Vertex_const_handle touching_vertex_end()   const { return this->_vtche; }
//#endif
	inline double weight() const { return this->_weight; }
	inline const Vector3& wgrad() const { return this->_wgrad; }
	//inline const Vector3& normal() const { return this->_normal; }
protected:
	virtual void finish() override;
};
template<typename P3>
void MAT_Polygons_3_Facet<P3>::finish()
{
	// Try to compute flux going out of along each edge. Since the shape of
	// Voronoi facet could be really "bad" we need to be careful here -
	// we assign 0 flux value to those. This is okay since e_flux values
	// only are used to identify stable edges and the algorithm does not
	// require to capture every stable edge.

	// Assume the degree isn't too large because in that case we would have problems with underflow.
	size_t degree = this->degree();
	Kernel::FT degree_ft = Kernel::FT(degree);

	// Get facet "center" (not really centroid [except for triangle]... I don't know which center this is)
	Kernel::FT x = 0.0, y = 0.0, z = 0.0;
	typename base::Edge_around_facet_circulator e = this->edges_circ(), e_end = e;
	CGAL_For_all(e, e_end) { const Point3& p = e->vertex()->point(); x += p.x(); y += p.y(); z += p.z(); }
	const Point3 center(x/degree_ft, y/degree_ft, z/degree_ft);

	// Calculate the flux of each edge
	double ddn = 0.0, area = 0.0;
	const Vector3 wgrad = normalized(this->_wgrad);
	CGAL_For_all(e, e_end)
	{
		const Point3 &p0 = e->vertex()->point(), &p1 = e->next_around_facet()->vertex()->point();
		double l_edge = dbl_sqrt(CGAL::squared_distance(p0, p1));
		const Vector3 nvedge = (p0-p1)/l_edge, wgrad_nproj = wgrad - (wgrad * nvedge) * nvedge;
		double flux = CGAL::sign(wgrad_nproj * (p0 - center)) * dbl_sqrt(wgrad_nproj.squared_length());
		e->_flux = flux;
		ddn += l_edge * e->flux();
		area += dbl_sqrt(CGAL::squared_area(center, p0, p1));
	}

	// For normal Voronoi facet, ddn should be 0; additionally make sure we had a significant area
	if (ddn > 1e-10 || area < 1e-15) // TODO: better limits?
	{
		typename base::Edge_around_facet_circulator e = this->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end) { e->_flux = 0.0; }
	}
}

template<typename P3>
class MAT_Polygons_3_Edge : public Polygons_3_Edge<P3>
{
	friend P3;
	friend typename P3::Facet;
	typedef Polygons_3_Facet<P3> base;

	double _flux;
public:
	FLAGS_DECL(Flags, Skeleton, StableSkeleton, Visited);
	Flags flags;
	MAT_Polygons_3_Edge(typename base::Vertex_handle v, typename base::Facet_handle f) : Polygons_3_Edge<P3>(v, f), flags(0) { }
	inline double flux() const { return this->_flux; }
};

extern template class Polygons_3<Kernel, Polygons_3_Vertex, MAT_Polygons_3_Edge, MAT_Polygons_3_Facet>;
extern template class MAT_Polygons_3_Facet<Polygons_3<Kernel, Polygons_3_Vertex, MAT_Polygons_3_Edge, MAT_Polygons_3_Facet>>;
typedef Polygons_3<Kernel, Polygons_3_Vertex, MAT_Polygons_3_Edge, MAT_Polygons_3_Facet> MAT;
