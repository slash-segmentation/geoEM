#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include <limits>
#include <vector>

//////////////////////////////////////////////////////////////////////////////////////////////
// Internal classes and typedefs for polyhedron extra geodesic data.
//////////////////////////////////////////////////////////////////////////////////////////////

#define MAX_DIST std::numeric_limits<double>::max()

struct GeodesicDistance
{
    double distance;
    Vector3 derivative;
    explicit GeodesicDistance(double distance = MAX_DIST, const Vector3& derivative = Vector3()) : distance(distance), derivative(derivative) {}
};
typedef handle_map<Polyhedron3::Vertex_const_handle, GeodesicDistance> GeodesicDistances;
inline double get_dist(const GeodesicDistances& gds, const Polyhedron3::Vertex_const_handle& v)
{
    GeodesicDistances::const_iterator gd = gds.find(v); return (gd == gds.end()) ? MAX_DIST : gd->second.distance;
}

struct VertexData
{
    GeodesicDistances geod;
    size_t n_inc_facets;
    VertexData() : geod(), n_inc_facets() { }
    VertexData(Polyhedron3::Vertex_const_handle v) : geod(), n_inc_facets(num_inc_facets(v)) { }
private:
    static size_t num_inc_facets(Polyhedron3::Vertex_const_handle v) { size_t n = 0; FOR_FACETS_AROUND_VERTEX(v, f) { if (f != Polyhedron3::Face_handle()) { ++n; } } return n; }
};
typedef std::vector<VertexData> VData;

struct HalfedgeData
{
    Kernel::FT squared_length;
    double length;
    Polyhedron3::Halfedge_const_handle eid; // the min2() of e and e->opposite(), basically making the edge unique regardless of half-edge
    HalfedgeData() : squared_length(), length(), eid() { }
    HalfedgeData(Polyhedron3::Halfedge_const_handle e) :
        squared_length(CGAL::squared_distance(e->vertex()->point(), e->prev()->vertex()->point())),
        length(dbl_sqrt(squared_length)), eid(min2(e, e->opposite())) { }
};
typedef std::vector<HalfedgeData> HEData;
