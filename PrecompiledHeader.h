#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// General includes and defines. This is expected to be a pre-compiled header to speed up
// compiling. Won't change much. All files do not depend on this header giving certain headers.
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
// Don't need Windows, but it may be included later and we want to make it usable
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#ifndef _CRT_NONSTDC_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#endif
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#undef min
#undef max
#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#define _nextafter(a,b) ::nextafter(a,b) // workaround for a CGAL bug when compiling with GCC on Windows
#endif
#else // #ifdef _WIN32
// *nix-specific includes and macros...
#endif // #ifdef _WIN32

// Include a bunch of CGAL stuff (heavy on the templates...)
#include <CGAL/utils.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Random.h>
#include <CGAL/Handle_hash_function.h>

// Common STL Includes
#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

// Common Boost Includes
#include <boost/config.hpp>
#include <boost/thread.hpp>

