#include "Skeleton.hpp"

#include <CGAL/extract_mean_curvature_flow_skeleton.h>

#include <boost/graph/adjacency_list.hpp>

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron3> Skeletonization;

///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton using mean curvature flow.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P)
{
	Skeletonization::Skeleton *skel = new Skeletonization::Skeleton();
	CGAL::extract_mean_curvature_flow_skeleton(*P, *skel);
	return skel;
}
