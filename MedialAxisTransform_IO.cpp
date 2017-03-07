#include "MedialAxisTransform_IO.hpp"

#include <iostream>
#include <fstream>

template <typename Iterator> static Iterator add(Iterator i, ptrdiff_t n) { std::advance(i, n); return i; }

MAT* load_mat(const std::string filename, Polyhedron3* mesh)
{
	std::fstream file(filename.c_str(), std::ios_base::in | std::ios_base::binary);
	CGAL::set_binary_mode(file);

	// Setup base structure
	unsigned int n_verts, n_edges, n_facets;
	file.read(reinterpret_cast<char*>(&n_verts), sizeof(n_verts));
	file.read(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));
	file.read(reinterpret_cast<char*>(&n_facets), sizeof(n_facets));
	MAT* mat = new MAT(n_verts, n_edges, n_facets);

	// Read vertices
	std::vector<MAT::Vertex_handle> vertices;
	vertices.reserve(n_verts);
	for (unsigned int i = 0; i < n_verts; ++i)
	{
		Point3 p;
		file >> p;
		vertices.push_back(mat->add_vertex(p));
	}

	// Read facets and fluxes
	std::vector<MAT::Vertex_handle> f_verts;
	for (unsigned int i = 0; i < n_facets; ++i)
	{
		unsigned int V1, V2, n_f_verts;
		double weight;
		Vector3 wgrad;
		file.read(reinterpret_cast<char*>(&V1), sizeof(V1));
		file.read(reinterpret_cast<char*>(&V2), sizeof(V2));
		file.read(reinterpret_cast<char*>(&weight), sizeof(weight));
		file >> wgrad;

		std::vector<double> fluxes;
		file.read(reinterpret_cast<char*>(&n_f_verts), sizeof(n_f_verts));
		f_verts.clear();
		f_verts.reserve(n_f_verts);
		for (unsigned int j = 0; j < n_f_verts; ++j)
		{
			unsigned int v;
			double f;
			file.read(reinterpret_cast<char*>(&v), sizeof(v));
			file.read(reinterpret_cast<char*>(&f), sizeof(f));
			f_verts.push_back(vertices[v]);
			fluxes.push_back(f);
		}

		MAT::Facet_handle f = mat->add_facet(f_verts.begin(), f_verts.end(), MAT::Facet(add(mesh->vertices_begin(), V1), add(mesh->vertices_begin(), V2), weight, wgrad));
		MAT::Edge_around_facet_const_circulator e = f->edges_circ();
		for (unsigned int j = 0; j < n_f_verts; ++j, ++e) { if (abs(fluxes[j] - e->flux()) > abs(0.0001*fluxes[j])) { std::cerr << "Warning! flux mismatch: " << fluxes[j] << " != " << e->flux() << "     " << abs(fluxes[j] - e->flux()) << " > " << 0.0001*fluxes[j] << std::endl; } }
	}

	return mat;
}
void save_mat(const std::string filename, MAT* mat, Polyhedron3* mesh)
{
	std::fstream file(filename.c_str(), std::ios_base::out | std::ios_base::binary);
	CGAL::set_binary_mode(file);

	// Write base structure
	unsigned int n_verts = (unsigned int)mat->size_of_vertices(), n_edges = (unsigned int)mat->size_of_edges(), n_facets = (unsigned int)mat->size_of_facets();
	file.write(reinterpret_cast<char*>(&n_verts), sizeof(n_verts));
	file.write(reinterpret_cast<char*>(&n_edges), sizeof(n_edges));
	file.write(reinterpret_cast<char*>(&n_facets), sizeof(n_facets));

	// Write vertices
	for (MAT::Vertex_const_iterator i = mat->vertices_begin(), end = mat->vertices_end(); i != end; ++i)
	{
		file << i->point();
	}

	// Write facets and fluxes
	for (MAT::Facet_const_iterator i = mat->facets_begin(), end = mat->facets_end(); i != end; ++i)
	{
		unsigned int V1 = (unsigned int)std::distance(Polyhedron3::Vertex_const_iterator(mesh->vertices_begin()), i->touching_vertex_start()),
			V2 = (unsigned int)std::distance(Polyhedron3::Vertex_const_iterator(mesh->vertices_begin()), i->touching_vertex_end()), n_f_verts = (unsigned int)i->degree();
		double weight = i->weight();
		Vector3 wgrad = i->wgrad();
		file.write(reinterpret_cast<char*>(&V1), sizeof(V1));
		file.write(reinterpret_cast<char*>(&V2), sizeof(V2));
		file.write(reinterpret_cast<char*>(&weight), sizeof(weight));
		file << wgrad;

		file.write(reinterpret_cast<char*>(&n_f_verts), sizeof(n_f_verts));
		MAT::Edge_around_facet_const_circulator j = i->edges_circ(), j_end = j;
		CGAL_For_all(j, j_end)
		{
			unsigned int v = (unsigned int)(j->vertex() - mat->vertices_begin());
			double f = j->flux();
			file.write(reinterpret_cast<char*>(&v), sizeof(v));
			file.write(reinterpret_cast<char*>(&f), sizeof(f));
		}
	}
}
void dump_mat(MAT* mat)
{
	for (MAT::Facet_const_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		std::cout << tostr(f->touching_vertex_start()->point()) << "\t" << tostr(f->touching_vertex_end()->point()) << "\t" << f->weight() << "\t" << tostr(f->wgrad()); //<< "\t" << tostr(f->normal);
		MAT::Facet::Edge_around_facet_const_circulator e = f->edges_circ(), eend = e;
		do
		{
			std::cout << "\t" << tostr(e->vertex()->point()) << "\t" << e->flux();
		} while (++e != eend);
		std::cout << std::endl;
	}
}
