#pragma once

#include <iostream>

#include "GeometryTypes.hpp"

///////////////////////////////////////////////////////////////////////////////////////
// This file declares many functions for creating and working with skeletons.
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
// Construct the curve-skeleton of a polygon mesh using the mean curvature flow method.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P);


///////////////////////////////////////////////////////////////////////////////////////
// Write a skeleton to a file.
///////////////////////////////////////////////////////////////////////////////////////
void write_skeleton(const Skeleton3* S, std::ofstream& output);
inline void write_skeleton(const Skeleton3* S, const char* filename, int precision=10)
{
    std::ofstream output(filename);
    output << std::setprecision(10);
    write_skeleton(S, output);
    output.close();
}


///////////////////////////////////////////////////////////////////////////////////////
// Construct the "graph" version of a skeleton from the "raw" skeleton. The two
// representations are not linked in any way, so modifying one will not modify the
// other.
///////////////////////////////////////////////////////////////////////////////////////
SkeletonGraph3* construct_skeleton_graph(const Skeleton3* S);


///////////////////////////////////////////////////////////////////////////////////////
// Remove coordinates in the skeleton by removing all collinear coordinates. This will
// not remove branch points or entire branches but may simplify a branch down to just
// its endpoints.
///////////////////////////////////////////////////////////////////////////////////////
void skeleton_remove_collinear(SkeletonGraph3* SG);

///////////////////////////////////////////////////////////////////////////////////////
// Reduce coordinates in the skeleton by removing all collinear coordinates and
// coordinates where the two neighbors form a triangle with squared area less than the
// threshold. If a negative value is provided for the threshold then the threshold is
// automatically determined as the lowest 20% of the squared areas which means that at
// least 20% of the vertices will be removed (but likely more). Default threshold is
// negative. This will not remove branch points or entire branches but may simplify a
// branch down to just its endpoints.
///////////////////////////////////////////////////////////////////////////////////////
void skeleton_reduce(SkeletonGraph3* SG, double threshold = -1.0);
