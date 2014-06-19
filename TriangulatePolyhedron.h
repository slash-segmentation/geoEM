#pragma once

#include "GeometryTypes.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Triangulate the surface of a polyhedron, making sure that every facet is a triangle.
// Many tools like the AABB_tree (thus the Intersection tool and point-in-polyhedron function)
// require that every facet is a triangle.
// Assumes (P->is_valid() && is_not_degenerate(P)) is true
///////////////////////////////////////////////////////////////////////////////////////////////////

void triangulate_polyhedron(Polyhedron3* P);
