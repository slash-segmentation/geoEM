#include "MedialAxisTransform_Types_Geodesic.hpp"
#include "MedialAxisTransform_Types_Triangulation.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////
// Building and annotating the triangulation and mesh.
//////////////////////////////////////////////////////////////////////////////////////////////
struct AddIncidentVertices
{
	Vertex_extra* ve;
	inline AddIncidentVertices(Vertex_extra* ve) : ve(ve) { }
	inline AddIncidentVertices& operator=(const Triangulation::Vertex_handle& v) { this->ve->geod[v->original] = GeodesicDistance(-1); return (*this); }
	inline AddIncidentVertices& operator*()     { return (*this); }
	inline AddIncidentVertices& operator++()    { return (*this); }
	inline AddIncidentVertices& operator++(int) { return (*this); }
};
void init(Triangulation& T, Polyhedron3* mesh)
{
	// Build Delaunay triangulation (the dual of the Voronoi diagram)
	// We record the original polyhedral vertex in each Delaunay vertex
	T.tds().vertices().reserve(mesh->size_of_vertices());
	T.tds().cells().reserve(8*mesh->size_of_vertices()); // crude approximation
	//   The following are simple definitions:
	//     V_infinite = 1
	//     V_all = V_finite + V_infinite
	//     E_all = E_finite + E_infinite
	//     F_all = F_finite + F_infinite
	//     C_all = C_finite + C_infinite
	//     E_infinite = V_convex_hull <= V_finite
	//     F_infinite = E_convex_hull <= E_finite
	//     C_infinite = F_convex_hull <= F_finite
	//   We have V_finite.
	//
	//   Euler's characteristic for the simplexes made: (found through example)
	//      1 = V_finite   - E_finite   + F_finite   - C_finite
	//     -1 = V_infinite - E_infinite + F_infinite - C_infinite
	//      0 = V_all      - E_all      + F_all      - C_all
	//
	//   Euler's characteristic for the convex hull made:
	//      2 = V_convex_hull - E_convex_hull + F_convex_hull
	//
	//   Some other properties found through example: (and make sense with informal proofs)
	//     C_infinite = 2/3*F_infinite      [F_convex_hull = 2/3*E_convex_hull]
	//     C_all = 1/2*F_all
	//     F_infinite = 3*(E_infinite - 2)  [E_convex_hull = 3*(V_convex_hull - 2)]
	//     C_infinite = 2*(E_infinite - 2)  [F_convex_hull = 2*(V_convex_hull - 2)]
	//
	//   Derived from the above:
	//     C_all = E_all - V_all
	//     C_all = E_infinite + E_finite - V_finite - 1
	//     C_all = V_convex_hull + E_finite - V_finite - 1
	//
	//     C_finite = 1/2*F_finite - 1/2*V_convex_hull - 1
	//
	//   We need more relationships for the finite ones.
	for (Polyhedron3::Vertex_iterator v = mesh->vertices_begin(), end = mesh->vertices_end(); v != end; ++v) { T.insert(v->point())->original = v; }
	//std::cout << "Mesh has: " << mesh->size_of_vertices() << " vertices, " << mesh->size_of_halfedges() / 2 << " edges, " << mesh->size_of_facets() << " facets" << std::endl;
	//std::cout << "Triangulation has: [finite] " << T.number_of_vertices()   << " vertices, " << T.number_of_finite_edges() << " edges, " << T.number_of_finite_facets() << " facets, " << T.number_of_finite_cells() << " cells" << std::endl;
	//std::cout << "Triangulation has: [all]    " << T.number_of_vertices()+1 << " vertices, " << T.number_of_edges()        << " edges, " << T.number_of_facets()        << " facets, " << T.number_of_cells()        << " cells" << std::endl;

	// Setup the mesh with extra data
	assert(mesh->extras_are_not_used());
	for (Triangulation::Finite_vertices_iterator vit = T.finite_vertices_begin(), end = T.finite_vertices_end(); vit != end; ++vit)
	{
		T.finite_adjacent_vertices(vit, AddIncidentVertices(vit->original->extra() = new Vertex_extra(vit->original)));
	}
	for (Polyhedron3::Halfedge_iterator e = mesh->halfedges_begin(), end = mesh->halfedges_end(); e != end; ++e) { e->extra() = new Halfedge_extra(e); }

	// Make inside vs. outside cells and facet intersections
#ifdef CGAL_OLD_FACET_TREE
	FacetTree tree(mesh->facets_begin(), mesh->facets_end());
#else
	FacetTree tree(mesh->facets_begin(), mesh->facets_end(), *mesh);
#endif
	for (Triangulation::Finite_cells_iterator c = T.finite_cells_begin(), end = T.finite_cells_end(); c != end; ++c)
	{
		const Point3& p = c->circumcenter();
		c->inside = point_in_polyhedron(p, tree); // for computing the out curve-skeletons: c->inside = !...;
		//if (c->inside) { continue; }
		for (int j = 0; j < 4; ++j)
		{
			const Triangulation::Cell_handle cn = c->neighbor(j);
			if (T.is_infinite(cn) /*|| cn->inside */ || !boost::logic::indeterminate(c->facet_intersects_mesh[j])) { continue; }
			const Point3& q = cn->circumcenter();
			bool intersects = CGAL::to_double(CGAL::squared_distance(p, q)) >= 1e-20 && tree.do_intersect(Segment3(p, q)); // TODO: use another method for checking distance
			c ->facet_intersects_mesh[j           ] = intersects;
			cn->facet_intersects_mesh[cn->index(c)] = intersects;
		}
	}
}
void cleanup(Polyhedron3* mesh)
{
	for (Polyhedron3::Vertex_iterator   v = mesh->vertices_begin(),  end = mesh->vertices_end();  v != end; ++v) { delete v->extra(); v->extra() = nullptr; }
	for (Polyhedron3::Halfedge_iterator e = mesh->halfedges_begin(), end = mesh->halfedges_end(); e != end; ++e) { delete e->extra(); e->extra() = nullptr; }
}
