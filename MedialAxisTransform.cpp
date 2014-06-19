#include "MedialAxisTransform.h"
#include "MedialAxisTransform_Types_Geodesic.h"
#include "MedialAxisTransform_Types_Triangulation.h"

#include <iostream>
#include <vector>

template class Polygons_3<Kernel, Polygons_3_Vertex, MAT_Polygons_3_Edge, MAT_Polygons_3_Facet>;
template class MAT_Polygons_3_Facet<MAT>;
template class CGAL::Triangulation_data_structure_3<TriVertex, TriCell>;
template class CGAL::Delaunay_triangulation_3<Kernel, TDS>;

void init(Triangulation& T, Polyhedron3* mesh);
void cleanup(Polyhedron3* mesh);
void compute_geodesic(Polyhedron3* mesh);
void compute_geodesic_using_cache(Polyhedron3* mesh, const std::string filename, const std::string objname);

//// Debugging functions:
//void dump_geod(Polyhedron3* mesh)
//{
//	for (Polyhedron3::Vertex_const_handle i = mesh->vertices_begin(), end = mesh->vertices_end(); i != end; ++i)
//	{
//		std::cout << tostr(i->point()) << ": " << std::endl;
//		std::vector<Polyhedron3::Vertex_const_handle> keys;
//		keys.reserve(i->extra()->geod.size());
//		for (GeodesicDistances::const_iterator j = i->extra()->geod.begin(), end = i->extra()->geod.end(); j != end; ++j) { keys.push_back(j->first); }
//		std::sort(keys.begin(), keys.end());
//		for (std::vector<Polyhedron3::Vertex_const_handle>::const_iterator j = keys.begin(), end = keys.end(); j != end; ++j)
//		{
//			const GeodesicDistance& gd = i->extra()->geod.at(*j);
//			std::cout << "\t" << tostr((*j)->point()) << ": " << gd.distance << ", " << tostr(gd.derivative) << std::endl;
//		}
//	}
//	std::cout << std::endl;
//}
//static size_t cell_number(Triangulation::Cell_handle c, const Triangulation& T) { return std::distance(T.cells_begin(), c); }


///////////////////////////////////////////////////////////////////////////////////////
// Compute the medial axis transform from the geodesic distances in the triangulation.
///////////////////////////////////////////////////////////////////////////////////////
static void compute_mat(Triangulation& T, MAT* mat)
{
	std::vector<MAT::Vertex_handle> vs;
	for (Triangulation::Finite_edges_iterator e = T.finite_edges_begin(), end = T.finite_edges_end(); e != end; ++e)
	{
		const Triangulation::Cell_handle ch = e->first;
		const Polyhedron3::Vertex_const_handle V1 = ch->vertex(e->second)->original, V2 = ch->vertex(e->third)->original;
		const Triangulation::Cell_circulator begin = T.incident_cells(*e);

		// Check whether this Delaunay edge is on the surface or not (if it is, don't use it)
		// If not an inner cell or if the segment between neighboring circumcenters intersects the mesh then the Delaunay edge is on the surface
		bool x = true;
		size_t n = 0;
		for (Triangulation::Cell_circulator c = begin; x && c->inside && !(bool)c->facet_intersects_mesh[c->index(std::next(c))]; x = ++c != begin, ++n);
		if (x || n < 3) { continue; }

		// Create the facet with the MAT information and add the facet vertices
		vs.clear();
		Triangulation::Cell_circulator c = begin;
		do
		{
			if (c->mat_vertex == MAT::Vertex_handle()) { c->mat_vertex = mat->add_vertex(c->circumcenter()); }
			vs.push_back(c->mat_vertex);
		}
		while (++c != begin);
		mat->add_facet(vs.begin(), vs.end(), MAT::Facet(V1, V2));
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// Construct the medial axis of the polygon mesh, possibly using cached geodesic data.
///////////////////////////////////////////////////////////////////////////////////////
MAT* construct_medial_axis_transform(Polyhedron3* mesh, const std::string filename, const std::string obj_name)
{
	std::cerr << "Initializing medial axis transform calculator... " << std::endl;
	Triangulation T;
	init(T, mesh);
	compute_geodesic_using_cache(mesh, filename, obj_name); // uses information in mesh based on T

	std::cerr << "Computing medial axis transform..." << std::endl;
	MAT* mat = new MAT(T.number_of_finite_cells(), 3*T.number_of_finite_edges(), T.number_of_finite_edges()); // over/under/over-estimates, but the best we can really do
	compute_mat(T, mat); // uses information in T based on mesh
	mat->shrink_to_fit(); // remove wasted space in storage

	// Done with triangulation and mesh, everything we need is in mat
	T.clear();
	cleanup(mesh);
	return mat;
}
///////////////////////////////////////////////////////////////////////////////////////
// Construct the medial axis of the polygon mesh.
///////////////////////////////////////////////////////////////////////////////////////
MAT* construct_medial_axis_transform(Polyhedron3* mesh)
{
	std::cerr << "Initializing medial axis transform calculator... " << std::endl;
	Triangulation T;
	init(T, mesh);
	compute_geodesic(mesh); // uses information in mesh based on T

	std::cerr << "Computing medial axis transform..." << std::endl;
	MAT* mat = new MAT(T.number_of_finite_cells(), 3*T.number_of_finite_edges(), T.number_of_finite_edges()); // over/under/over-estimates, but the best we can really do
	compute_mat(T, mat); // uses information in T based on mesh
	mat->shrink_to_fit(); // remove wasted space in storage

	// Done with triangulation and mesh, everything we need is in mat
	T.clear();
	cleanup(mesh);
	return mat;
}
