#include "MedialAxisTransform.hpp"
#include "MedialAxisTransform_Types_Geodesic.hpp"
#include "MedialAxisTransform_Types_Triangulation.hpp"

#include <iostream>
#include <vector>

template class Polygons_3<Kernel, Polygons_3_Vertex, MAT_Polygons_3_Edge, MAT_Polygons_3_Facet>;
template class MAT_Polygons_3_Facet<MAT>;
template class CGAL::Triangulation_data_structure_3<TriVertex, TriCell>;
template class CGAL::Delaunay_triangulation_3<Kernel, TDS>;

FLAGS_DEFN(MAT_Polygons_3_Facet<MAT>::Flags, StableSkeleton, NotStableSkeleton, Boundary, Eaten, InHeap);
FLAGS_DEFN(MAT_Polygons_3_Edge<MAT>::Flags, Skeleton, StableSkeleton, Visited);

void mat_init(Triangulation& T, Polyhedron3* mesh, VData& vdata);
void compute_geodesic(Polyhedron3* mesh, VData& vdata);
void compute_geodesic_using_cache(Polyhedron3* mesh, VData& vdata, const std::string filename, const std::string objname);

//// Debugging functions:
//void dump_geod(const Polyhedron3* mesh, const VData& vdata)
//{
//  for (Polyhedron3::Vertex_const_handle i = mesh->vertices_begin(), end = mesh->vertices_end(); i != end; ++i)
//  {
//      std::cout << tostr(i->point()) << ": " << std::endl;
//      std::vector<Polyhedron3::Vertex_const_handle> keys;
//      keys.reserve(vdata[i->id()]->geod.size());
//      for (GeodesicDistances::const_iterator j = vdata[i->id()]->geod.begin(), end = vdata[i->id()]->geod.end(); j != end; ++j) { keys.push_back(j->first); }
//      std::sort(keys.begin(), keys.end());
//      for (std::vector<Polyhedron3::Vertex_const_handle>::const_iterator j = keys.begin(), end = keys.end(); j != end; ++j)
//      {
//          const GeodesicDistance& gd = vdata[i->id()]->geod.at(*j);
//          std::cout << "\t" << tostr((*j)->point()) << ": " << gd.distance << ", " << tostr(gd.derivative) << std::endl;
//      }
//  }
//  std::cout << std::endl;
//}
//static size_t cell_number(Triangulation::Cell_handle c, const Triangulation& T) { return std::distance(T.cells_begin(), c); }


///////////////////////////////////////////////////////////////////////////////////////
// Compute the medial axis transform from the geodesic distances in the triangulation.
///////////////////////////////////////////////////////////////////////////////////////
static void compute_mat(Triangulation& T, MAT* mat, const VData& vdata)
{
    std::vector<MAT::Vertex_handle> vs;
    for (Triangulation::Finite_edges_iterator e = T.finite_edges_begin(), end = T.finite_edges_end(); e != end; ++e)
    {
        const Triangulation::Cell_handle ch = e->first;
        const Polyhedron3::Vertex_const_handle V1 = ch->vertex(e->second)->original, V2 = ch->vertex(e->third)->original;
        const Triangulation::Cell_circulator begin = T.incident_cells(*e);

        // Check whether this Delaunay edge is on the surface or not (if it is, don't use it)
        // If not an inner cell or if the segment between neighboring circumcenters intersects the mesh then the Delaunay edge is on the surface
        bool x = true;
        size_t n = 0;
        for (Triangulation::Cell_circulator c = begin; x && c->inside && !(bool)c->facet_intersects_mesh[c->index(std::next(c))]; x = ++c != begin, ++n);
        if (x || n < 3) { continue; }

        // Create the facet with the MAT information and add the facet vertices
        vs.clear();
        Triangulation::Cell_circulator c = begin;
        do
        {
            if (c->mat_vertex == MAT::Vertex_handle()) { c->mat_vertex = mat->add_vertex(c->circumcenter()); }
            vs.push_back(c->mat_vertex);
        }
        while (++c != begin);
        mat->add_facet(vs.begin(), vs.end(), MAT::Facet(V1, V2, vdata));
    }
}


///////////////////////////////////////////////////////////////////////////////////////
// Construct the medial axis of the polygon mesh, possibly using cached geodesic data.
///////////////////////////////////////////////////////////////////////////////////////
MAT* construct_medial_axis_transform(Polyhedron3* mesh, const std::string filename, const std::string obj_name)
{
    std::cerr << "Initializing medial axis transform calculator... " << std::endl;
    Triangulation T;
    VData vdata;
    mat_init(T, mesh, vdata);
    compute_geodesic_using_cache(mesh, vdata, filename, obj_name); // uses information in mesh based on T

    std::cerr << "Computing medial axis transform..." << std::endl;
    MAT* mat = new MAT(T.number_of_finite_cells(), 3*T.number_of_finite_edges(), T.number_of_finite_edges()); // over/under/over-estimates, but the best we can really do
    compute_mat(T, mat, vdata); // uses information in T based on mesh
    mat->shrink_to_fit(); // remove wasted space in storage

    // Done with triangulation and mesh, everything we need is in mat
    T.clear();
    return mat;
}
///////////////////////////////////////////////////////////////////////////////////////
// Construct the medial axis of the polygon mesh.
///////////////////////////////////////////////////////////////////////////////////////
MAT* construct_medial_axis_transform(Polyhedron3* mesh)
{
    std::cerr << "Initializing medial axis transform calculator... " << std::endl;
    Triangulation T;
    VData vdata;
    mat_init(T, mesh, vdata);
    compute_geodesic(mesh, vdata); // uses information in mesh based on T

    std::cerr << "Computing medial axis transform..." << std::endl;
    MAT* mat = new MAT(T.number_of_finite_cells(), 3*T.number_of_finite_edges(), T.number_of_finite_edges()); // over/under/over-estimates, but the best we can really do
    compute_mat(T, mat, vdata); // uses information in T based on mesh
    mat->shrink_to_fit(); // remove wasted space in storage

    // Done with triangulation and mesh, everything we need is in mat
    T.clear();
    return mat;
}
