#include "Skeleton.hpp"

#include <CGAL/extract_mean_curvature_flow_skeleton.h>

///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton using mean curvature flow.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P)
{
    Skeleton3* S = new Skeleton3();
    CGAL::extract_mean_curvature_flow_skeleton(*P, *S);
    return S;
}
