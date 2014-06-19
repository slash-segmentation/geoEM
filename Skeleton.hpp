#pragma once

#include "GeometryTypes.hpp"
#include "MedialAxisTransform.hpp"

///////////////////////////////////////////////////////////////////////////////////////
// This file declares many functions for creating and working with skeletons.
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
// Construct the curve-skeleton of a polygon mesh using its medial axis with geodesic.
// The flux and weight threshold parameters control the aggressiveness of the trimming
// process.
//
// Weight threshold must be >= 0 (0 ensures nothing is trimmed due to its weight)
// Flux threshold must be [0,1)
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(MAT* mat, double flux_thd = 0.5, double wt_thd = 0.0);

///////////////////////////////////////////////////////////////////////////////////////
// Construct the "graph" version of a skeleton from the "raw" skeleton. The two
// representations are not linked in any way, so modifying one will not modify the
// other.
///////////////////////////////////////////////////////////////////////////////////////
SkeletonGraph3* construct_skeleton_graph(const Skeleton3* S);


///////////////////////////////////////////////////////////////////////////////////////
// Reduce the number of coordinates in the skeleton by removing all collinear
// coordinates and coordinates where the two neighbors form a triangle with squared
// area less than the threshold. It will not remove branch points or entire branches
// but may simplify a branch down to just its endpoints.
///////////////////////////////////////////////////////////////////////////////////////
void skeleton_reduce(SkeletonGraph3* SG, double threshold = 0.5);
