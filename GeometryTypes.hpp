#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// All of the basic 2D and 3D geometry types we will use and the CGAL Kernel we will be using.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <CGAL/utils.h>

#ifdef _DEBUG
#define CGAL_DONT_USE_LAZY_KERNEL
#define CGAL_NO_STATIC_FILTERS
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC_Kernel;
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel EPEC_Kernel;

// Using the inexact kernel is incredibly faster, but some algorithms may fail with it. Things are
// carefully written everything so far to support the inexact kernel but it wasn't always that way.
#ifdef KERNEL_INEXACT
typedef EPIC_Kernel Kernel;
#else
typedef EPEC_Kernel Kernel;
#endif

#ifdef INSTANTIATE_GEOM_TYPES
#define DECLARE_TEMPLATE(...)        template class __VA_ARGS__;
#else
#define DECLARE_TEMPLATE(...) extern template class __VA_ARGS__;
#endif

// 2D
#include <CGAL/Direction_2.h>
DECLARE_TEMPLATE(CGAL::Direction_2<Kernel>)
typedef CGAL::Direction_2<Kernel>   Direction2;

#include <CGAL/Vector_2.h>
//DECLARE_TEMPLATE(CGAL::Vector_2<Kernel>)
typedef CGAL::Vector_2<Kernel>      Vector2;

#include <CGAL/Point_2.h>
DECLARE_TEMPLATE(CGAL::Point_2<Kernel>)
typedef CGAL::Point_2<Kernel>       Point2;

#include <CGAL/Segment_2.h>
DECLARE_TEMPLATE(CGAL::Segment_2<Kernel>)
typedef CGAL::Segment_2<Kernel>     Segment2;

#include <CGAL/Ray_2.h>
DECLARE_TEMPLATE(CGAL::Ray_2<Kernel>)
typedef CGAL::Ray_2<Kernel>         Ray2;

#include <CGAL/Line_2.h>
DECLARE_TEMPLATE(CGAL::Line_2<Kernel>)
typedef CGAL::Line_2<Kernel>		Line2;

#include <CGAL/Triangle_2.h>
DECLARE_TEMPLATE(CGAL::Triangle_2<Kernel>)
typedef CGAL::Triangle_2<Kernel>    Triangle2;

#include <CGAL/Polygon_2.h>
DECLARE_TEMPLATE(CGAL::Polygon_2<Kernel>)
typedef CGAL::Polygon_2<Kernel>     Polygon2;

#include <CGAL/Bbox_2.h>
typedef CGAL::Bbox_2                Bbox2;


// 3D
#include <CGAL/Direction_3.h>
DECLARE_TEMPLATE(CGAL::Direction_3<Kernel>)
typedef CGAL::Direction_3<Kernel>   Direction3;

#include <CGAL/Vector_3.h>
//DECLARE_TEMPLATE(CGAL::Vector_3<Kernel>)
typedef CGAL::Vector_3<Kernel>      Vector3;

#include <CGAL/Point_3.h>
DECLARE_TEMPLATE(CGAL::Point_3<Kernel>)
typedef CGAL::Point_3<Kernel>       Point3;

#include <CGAL/Segment_3.h>
DECLARE_TEMPLATE(CGAL::Segment_3<Kernel>)
typedef CGAL::Segment_3<Kernel>     Segment3;

#include <CGAL/Ray_3.h>
DECLARE_TEMPLATE(CGAL::Ray_3<Kernel>)
typedef CGAL::Ray_3<Kernel>         Ray3;

#include <CGAL/Line_3.h>
DECLARE_TEMPLATE(CGAL::Line_3<Kernel>)
typedef CGAL::Line_3<Kernel>        Line3;

#include <CGAL/Triangle_3.h>
DECLARE_TEMPLATE(CGAL::Triangle_3<Kernel>)
typedef CGAL::Triangle_3<Kernel>    Triangle3;

#include <CGAL/Plane_3.h>
DECLARE_TEMPLATE(CGAL::Plane_3<Kernel>)
typedef CGAL::Plane_3<Kernel>       Plane3;

#include "Polyhedron_3.hpp"
//DECLARE_TEMPLATE(Polyhedron_3<Kernel>)
typedef Polyhedron_3<Kernel> Polyhedron3;

#include "Polygons_3.hpp"
DECLARE_TEMPLATE(Polygons_3<Kernel>)
typedef Polygons_3<Kernel>			Polygons3;

#include "Skeleton_3.hpp"
DECLARE_TEMPLATE(SKELETON_3(Kernel))
typedef SKELETON_3(Kernel)			Skeleton3;

#include "SkeletonGraph_3.hpp"
DECLARE_TEMPLATE(SkeletonGraph_3<Kernel>)
typedef SkeletonGraph_3<Kernel>		SkeletonGraph3;

#include <CGAL/Bbox_3.h>
typedef CGAL::Bbox_3                Bbox3;


// Other
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
DECLARE_TEMPLATE(CGAL::AABB_tree<CGAL::AABB_traits<Kernel, CGAL::AABB_face_graph_triangle_primitive<const Polyhedron3>>>)
typedef CGAL::AABB_tree<CGAL::AABB_traits<Kernel, CGAL::AABB_face_graph_triangle_primitive<const Polyhedron3>>> FacetTree;
