#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include <utility>
#include <vector>

typedef std::pair<std::vector<S3VertexDesc>, Polyhedron3*> Slice;
typedef std::vector<Slice> Slices;

Slices slice(const int group_sz, const Polyhedron3* P, const Skeleton3* S,
    const P3CVertexSet& facet_verts, bool verbose=false);
