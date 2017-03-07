///////////////////////////////////////////////////////////////////////////////
// Simple program to simplify a mesh, by removing edges between neighboring
// facets that are nearly coplanar. It tries to keep the surface area and
// volume of the polyhedron consistent.
///////////////////////////////////////////////////////////////////////////////

#include "GeometryTypes.hpp"
#include "IO.hpp"

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h> // Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>              // Simplification function
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h> // Stop-condition policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

namespace SMS = CGAL::Surface_mesh_simplification;

template <class ShouldStop>
int simplify(Polyhedron3& P, const ShouldStop& stop)
{
	/*
	TODO: should these be added to the end of the params line?
	#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
	#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
		.get_cost(SMS::Edge_length_cost<Polyhedron3>())
		.get_placement(SMS::Midpoint_placement<Polyhedron3>())
	*/

	auto params = CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, P)).halfedge_index_map(get(CGAL::halfedge_external_index, P));
	return SMS::edge_collapse(P, stop, params);
}

static int simplify(Polyhedron3* P, char* param)
{
	if (strrchr(param, '.') >= 0)
	{
		double p = atof(param);
		if (p >= 1.0 || p <= 0.0) { return -1; }
		return simplify(*P, SMS::Count_ratio_stop_predicate<Polyhedron3>(p));
	}
	else
	{
		long n = atol(param);
		if (n <= 0) { return -1; }
		return simplify(*P, SMS::Count_stop_predicate<Polyhedron3>(n));
	}
}

static void usage(const char* err = NULL)
{
	if (err) { std::cerr << err << std::endl << std::endl; }

	std::cerr << "Usage: mesh-simplifier [in_file] [n|p] [optional:out_file]" << std::endl << std::endl;
	std::cerr << "in_file  file to read mesh from" << std::endl;
	std::cerr << "n        maximum number of edges to keep, integer >0" << std::endl;
	std::cerr << "p        percentage of edges to keep, given as a number in (0.0, 1.0)" << std::endl;
	std::cerr << "out_file file to save to, defaults to the '${in_file}_simplified.off'" << std::endl;
	std::cerr << "" << std::endl;
	std::cerr << mesh_file_usage << std::endl;
	std::cerr << mesh_file_read_usage << std::endl;
	std::cerr << mesh_file_write_usage << std::endl;
}

int main(int argc, char** argv)
{
	if (argc == 1) { usage(); return 1; }
	if (argc < 3 || argc > 4) { usage("ERROR: invalid number of arguments"); return 1; }

	// Read polyhedron
	Polyhedron3* P;
	std::cout << "Reading mesh..." <<std::endl;
	try { P = read_mesh(argv[1]); }
	catch (std::exception& ex)
	{
		std::cerr << ex.what() << std::endl << std::endl;
		usage("ERROR: unable to read mesh from input file");
		return 2;
	}

	// Simplify polyhedron
	std::cout << "Simplifying mesh..." <<std::endl;
	int r = simplify(P, argv[2]);
	if (r == -1) { usage("ERROR: invalid number/percentage of edges to keep"); return 3; }

	// Write polyhedron
	std::cout << "Saving mesh..." <<std::endl;
	char* output = argc == 3 ? strcat(strcpy(new char[strlen(argv[1]) + 16], argv[1]), "_simplified.off") : argv[3];
	try { write_mesh(P, output); }
	catch (std::exception& ex)
	{
		std::cerr << ex.what() << std::endl << std::endl;
		usage("ERROR: unable to write mesh from output file");
		return 4;
	}
	if (argc == 3) { delete[] output; }

	// Output information about simplification process
	size_t original_nedges = P->size_of_halfedges()/2 + r;
	std::cout << r << " out of " << original_nedges << " (" << std::setprecision(3) << r * 100.0 / original_nedges << "%)" << " edges removed" << std::endl;
	return 0;
}
