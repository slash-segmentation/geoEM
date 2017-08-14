#pragma once

#include "GeometryTypes.hpp"

#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// This file declares many functions for creating and working with segmentations.
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
// Construct the segmentation of a polygon mesh using SDF values. There are several
// parameters that can be set, see:
//   http://doc.cgal.org/latest/Surface_mesh_segmentation/group__PkgSurfaceSegmentation.html#gabc864f396347009726b858434c6d8659
///////////////////////////////////////////////////////////////////////////////////////
std::vector<Polyhedron3*> compute_segmentation(const Polyhedron3* P,
	double cone_angle = 2.0/3.0*CGAL_PI, size_t number_of_rays = 25,
	size_t number_of_clusters = 5, double smoothing_lambda = 0.26);

///////////////////////////////////////////////////////////////////////////////////////
// Construct the segmentation of a polygon mesh using SDF created from a skeleton.
// There are a few parameters to control this process, see:
//    http://doc.cgal.org/latest/Surface_mesh_segmentation/group__PkgSurfaceSegmentation.html#ga8a429857a748922d0e8460619db69764
// The skeleton must come from the curvature flow method and cannot be saved and loaded
// as some metadata on the skeleton will be lost.
///////////////////////////////////////////////////////////////////////////////////////
std::vector<Polyhedron3*> compute_segmentation(/*const*/ Polyhedron3* P, const Skeleton3* S,
	size_t number_of_clusters = 5, double smoothing_lambda = 0.26);
