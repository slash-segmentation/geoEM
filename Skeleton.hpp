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
// Enumerate the branches of a skeleton. The callback recieves a vector of the skeleton
// points (in order) and the mesh surface vertices for each of those points as a vector
// of vectors.
///////////////////////////////////////////////////////////////////////////////////////
template <class _F>
struct SkeletonVisitor
{
    //typedef typename std::add_lvalue_reference<_F>::type F; // TODO: find a way to add & automatically if possible (change in function as well)
    typedef _F F;
    F callback;
    std::vector<Skeleton3::vertex_descriptor> vertices;
    SkeletonVisitor(F callback) : callback(callback) {}
    void start_new_polyline() { vertices.clear(); }
    void add_node(Skeleton3::vertex_descriptor v) { vertices.push_back(v); }
    void end_polyline() { if (vertices.size() != 0) { callback(vertices); } }
};
template <class F>
void skeleton_enum_branches(const Skeleton3* S, F callback)
{
    SkeletonVisitor<F> sv(callback);
    CGAL::split_graph_into_polylines(*S, sv);
}
template <class F>
void skeleton_enum_branches_as_pts(const Skeleton3* S, F callback)
{
    skeleton_enum_branches(S, [&](std::vector<Skeleton3::vertex_descriptor> vertices)
    {
        std::vector<Point3> pts;
        pts.reserve(vertices.size());
        BOOST_FOREACH(Skeleton3::vertex_descriptor v, vertices) { pts.push_back((*S)[v].point); }
        callback(pts);
    });
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
