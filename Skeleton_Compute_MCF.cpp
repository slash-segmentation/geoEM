#include "Skeleton.hpp"

#include <CGAL/Subdivision_method_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>

#include <boost/graph/adjacency_list.hpp>

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron3> Skeletonization;

///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton using mean curvature flow.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P, int loop_subdivisions)
{
	assert(loop_subdivisions >= 0);

	// Refine mesh
	Polyhedron3 Pr;
	if (loop_subdivisions > 0) {
		// If we refine the mesh we need to a copy first
		Pr = *P;
		CGAL::Subdivision_method_3::Loop_subdivision(Pr, loop_subdivisions);
		P = &Pr;
	}

	// Calculate skeletonization
	Skeletonization::Skeleton *skel = new Skeletonization::Skeleton();
	CGAL::extract_mean_curvature_flow_skeleton(*P, *skel);
	return skel;
}
