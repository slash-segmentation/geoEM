#pragma once

#include <iostream>

#include "GeometryTypes.hpp"

#include <CGAL/boost/graph/split_graph_into_polylines.h>

///////////////////////////////////////////////////////////////////////////////////////
// This file declares many functions for creating and working with skeletons.
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
// Construct the curve-skeleton of a polygon mesh using the mean curvature flow method.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P, double quality_factor=0.1);


///////////////////////////////////////////////////////////////////////////////////////
// Enumerate the branches of a skeleton. In many cases this allievates the need for
// creating a skeleton graph (except for enumerating branch points themselves and the
// branches from a branch point). It also provides the vertices on the surface of the
// mesh which the skeleton graph does not.
//
// The callback recieves a vector of the skeleton points (in order) and the mesh
// surface vertices for each of those points as a vector of vectors.
///////////////////////////////////////////////////////////////////////////////////////
template <class F>
struct SkeletonVisitor
{
    const Skeleton3& S;
    F& callback;
    std::vector<Point3> points;
    std::vector<std::vector<Polyhedron3::Vertex_handle> > vertices;
    SkeletonVisitor(const Skeleton3& S, F& callback) : S(S), callback(callback) {}
    void start_new_polyline() { points.clear(); vertices.clear(); }
    void add_node(Skeleton3::vertex_descriptor v)
    {
        points.push_back(S[v].point);
        vertices.push_back(S[v].vertices);
    }
    void end_polyline()
    {
        if (points.size() == 0) { return; }
        callback(points, vertices);
    }
};
template <class F>
void skeleton_enum_branches(const Skeleton3* S, F& callback)
{
    SkeletonVisitor<F> sv(*S, callback);
    CGAL::split_graph_into_polylines(*S, sv);
}


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
