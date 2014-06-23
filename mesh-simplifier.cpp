///////////////////////////////////////////////////////////////////////////////
// Simple program to simplify a mesh, by removing edges between neighboring
// facets that are nearly coplanar. It tries to keep the surface area and
// volume of the polyhedron consistent.
///////////////////////////////////////////////////////////////////////////////

// TODO: use generalized I/O functions to support non-OFF files

#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h> // Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>              // Simplification function
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h> // Stop-condition policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron3;

namespace SMS = CGAL::Surface_mesh_simplification;

template <class ShouldStop>
int simplify(Polyhedron3& P, const ShouldStop& stop)
{
	return SMS::edge_collapse(P, stop, CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index, P)).edge_index_map(boost::get(CGAL::edge_external_index, P)));
}

static int simplify(Polyhedron3& P, char* param)
{
	if (strrchr(param, '.') >= 0)
	{
		double p = atof(param);
		if (p >= 1.0 || p <= 0.0) { return -1; }
		return simplify(P, SMS::Count_ratio_stop_predicate<Polyhedron3>(p));
	}
	else
	{
		long n = atol(param);
		if (n <= 0) { return -1; }
		return simplify(P, SMS::Count_stop_predicate<Polyhedron3>(n));
	}
}

static void usage(const char* err = NULL)
{
	if (err) { std::cerr << err << std::endl << std::endl; }

	std::cerr << "Usage: mesh-simplifier [in_file] [n|p] [optional:out_file]" << std::endl << std::endl;
	std::cerr << "in_file  an OFF file to read with a single polyhedron" << std::endl;
	std::cerr << "n        maximum number of edges to keep, integer >0" << std::endl;
	std::cerr << "p        percentage of edges to keep, given as a number in (0.0, 1.0)" << std::endl;
	std::cerr << "out_file output OFF file to save to, defaults to the '${in_file}_simplified.off'" << std::endl;
}

int main(int argc, char** argv)
{
	if (argc == 1) { usage(); return 1; }
	if (argc < 3 || argc > 4) { usage("ERROR: invalid number of arguments"); return 1; }

	// Read polyhedron
	std::ifstream in(argv[1]);
	if (!in) { usage("ERROR: unable to read from input file"); return 2; }
	Polyhedron3 P; in >> P;
	if (!in) { in.close(); usage("ERROR: input file could not be parsed"); return 2; }
	in.close();

	// Simplify polyhedron
	int r = simplify(P, argv[2]);
	if (r == -1) { usage("ERROR: invalid number/percentage of edges to keep"); return 3; }

	// Write polyhedron
	char* output = argc == 3 ? strcat(strcpy(new char[strlen(argv[1]) + 16], argv[1]), "_simplified.off") : argv[3];
	std::ofstream out(output);
	if (!in) { usage("ERROR: unable to write to output file"); return 4; }
	out << P;
	out.close();
	if (argc == 3) { delete[] output; }

	// Output information about simplification process
	size_t original_nedges = P.size_of_halfedges()/2 + r;
	std::cout << r << " out of " << original_nedges << " (" << std::setprecision(3) << r * 100.0 / original_nedges << "%)" << " edges removed" << std::endl;
	return 0;
}
