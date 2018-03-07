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

// Using a vector backing makes everything slightly faster but prevents polyhedrons being used for
// many operations such as surface triangulation and mesh simplification. If vector is not used,
// we at least provide a memory pool for the list to allocate memory from to speed it up.
//#define POLYHEDRON_USE_VECTOR
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#if defined(POLYHEDRON_USE_VECTOR)
#include <CGAL/HalfedgeDS_vector.h>
//DECLARE_TEMPLATE(CGAL::Polyhedron_3<Kernel, ...>)
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_vector> Polyhedron3;
#else
#include <boost/pool/pool_alloc.hpp>
//DECLARE_TEMPLATE(CGAL::Polyhedron_3<Kernel, ...>)
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, boost::fast_pool_allocator<int>> Polyhedron3;
#endif

#include <boost/graph/adjacency_list.hpp>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
//DECLARE_TEMPLATE(CGAL::Mean_curvature_flow_skeletonization<Polyhedron3>)
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron3>::Skeleton Skeleton3;

#include "SkeletonGraph_3.hpp"
DECLARE_TEMPLATE(SkeletonGraph_3<Kernel>)
typedef SkeletonGraph_3<Kernel>		SkeletonGraph3;

#include <CGAL/Bbox_3.h>
typedef CGAL::Bbox_3                Bbox3;


// Other
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
//DECLARE_TEMPLATE(CGAL::AABB_tree<CGAL::AABB_traits<Kernel, CGAL::AABB_face_graph_triangle_primitive<const Polyhedron3>>>)
typedef CGAL::AABB_tree<CGAL::AABB_traits<Kernel, CGAL::AABB_face_graph_triangle_primitive<const Polyhedron3>>> FacetTree;

