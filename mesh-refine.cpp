///////////////////////////////////////////////////////////////////////////////
// Simple program to refine a mesh using various techniques.
///////////////////////////////////////////////////////////////////////////////

#include "GeometryTypes.hpp"
#include "Polyhedron3Utils.hpp"
#include "IO.hpp"
#include "Strings.hpp"

#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/subdivision_method_3.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
//#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
namespace PMP = CGAL::Polygon_mesh_processing;

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h> // Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>              // Simplification function
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h> // Stop-condition policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
namespace SMS = CGAL::Surface_mesh_simplification;

static void usage(const char* err = nullptr, int exit_code=0)
{
    if (err) { std::cerr << err << std::endl << std::endl; }

    std::cerr << "Usage: mesh-refine [in_file] options... [optional:out_file]" << std::endl << std::endl;
    std::cerr << "in_file  file to read mesh from" << std::endl;
    std::cerr << "out_file file to save to, defaults to the '${in_file}_refined.off'" << std::endl;
    std::cerr << "Specify one or more options on how to refine the mesh:" << std::endl;
    std::cerr << "  -f     fill holes and fix self-intersections in the mesh" << std::endl;
    std::cerr << "  -i f   isotropic remeshing, targeting f*(average edge length)" << std::endl;
    std::cerr << "  -I l   isotropic remeshing, targeting l" << std::endl;
    std::cerr << "  -l #   # loop subdivisions" << std::endl;
    std::cerr << "  -e f   edge collapse, keeping f*(number of edges) edges" << std::endl;
    std::cerr << "  -E n   edge collapse, keeping n edges" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << mesh_file_usage << std::endl;
    std::cerr << mesh_file_read_usage << std::endl;
    std::cerr << mesh_file_write_usage << std::endl;
    exit(exit_code);
}

static void show_info(const Polyhedron3* P)
{
    // Mesh details
    std::cout << "  Mesh contains " << P->size_of_facets() << " facets, " << P->size_of_halfedges() << " edges, and " << P->size_of_vertices() << " vertices" << std::endl;
    std::cout << "    Volume: " << PMP::volume(*P) << ", Surface Area: " << PMP::area(*P) << ", Avg Edge Length: " << avg_edge_length(P) << std::endl;

    // Check mesh
    if (!P->is_valid())         { std::cerr << "  Warning: mesh is not valid" << std::endl; return; }
    if (!is_not_degenerate(P))  { std::cerr << "  Warning: mesh is degenerate" << std::endl; return; }
    if (!P->is_closed())
    {
        std::cerr << "  Warning: mesh is not closed" << std::endl;
    }
    else if (!PMP::is_outward_oriented(*P)) { std::cerr << "  Warning: mesh is not outward oriented" << std::endl; }
    if (!P->is_pure_triangle()) { std::cerr << "  Warning: mesh is not pure triangle" << std::endl; }
    if (PMP::does_self_intersect(*P))
    {
        std::cerr << "  Warning: mesh is self-intersecting" << std::endl;
    }
    else if (!PMP::does_bound_a_volume(*P)) { std::cerr << "  Warning: mesh does not bound a volume" << std::endl; }
    //if (!is_single_component(P)) { std::cerr << "Warning: mesh is not single connected component" << std::endl; }
}

static void isotropic_remeshing_fraction(Polyhedron3* P, double fraction)
{
    std::cout << "Isotropic remeshing..." << std::endl;
    if (fraction <= 0) { throw std::invalid_argument("fraction must be >0"); }
    PMP::isotropic_remeshing(faces(*P), CGAL::to_double(avg_edge_length(P)*fraction), *P);
}

static void isotropic_remeshing_length(Polyhedron3* P, double length)
{
    std::cout << "Isotropic remeshing..." << std::endl;
    if (length <= 0) { throw std::invalid_argument("length must be >0"); }
    PMP::isotropic_remeshing(faces(*P), CGAL::to_double(length), *P);
}

static void loop_subdivisions(Polyhedron3* P, long long num)
{
    std::cout << "Loop subdivisions..." << std::endl;
    if (num <= 0) { throw std::invalid_argument("number must be >0"); }
    CGAL::Subdivision_method_3::Loop_subdivision(*P, num);
}

template <class ShouldStop>
static int edge_collapse(Polyhedron3& P, const ShouldStop& stop)
{
    auto params = CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, P)).halfedge_index_map(get(CGAL::halfedge_external_index, P));
    return SMS::edge_collapse(P, stop, params);
}

static void edge_collapse_fraction(Polyhedron3* P, double fraction)
{
    std::cout << "Edge collapsing..." << std::endl;
    if (fraction <= 0 || fraction >= 1) { throw std::invalid_argument("fraction must be between 0 and 1"); }
    edge_collapse(*P, SMS::Count_ratio_stop_predicate<Polyhedron3>(fraction));
}

static void edge_collapse_count(Polyhedron3* P, long long count)
{
    std::cout << "Edge collapsing..." << std::endl;
    if (count <= 0) { throw std::invalid_argument("count must be >0"); }
    edge_collapse(*P, SMS::Count_stop_predicate<Polyhedron3>(count));
}

static void remove_self_intersections(Polyhedron3* P, bool aggressive=false)
{
    CGAL::set_halfedgeds_items_id(*P);
    if (PMP::does_self_intersect(*P))
    {
        std::cout << "Removing self-intersections..." << (aggressive ? " (aggressively)" : "") << std::endl;
        std::vector<std::pair<P3Facet, P3Facet>> out;
        PMP::self_intersections(*P, std::back_inserter(out));
        handle_set<P3Facet> to_remove;
        for (auto& si : out) { to_remove.insert(si.first); to_remove.insert(si.second); }
        if (aggressive)
        {
            for (auto& si : out)
            {
                FOR_FACETS_AROUND_FACET(si.first, f) { to_remove.insert(f); }
                FOR_FACETS_AROUND_FACET(si.second, g) { to_remove.insert(g); }
            }
        }
        for (auto f : to_remove) { P->erase_facet(f->halfedge()); }
        std::cout << "  " << to_remove.size() << " facets causing " << out.size() << " self intersections removed" << std::endl;
    }
}

static void fill_holes(Polyhedron3* P)
{
    CGAL::set_halfedgeds_items_id(*P);
    if (!P->is_closed())
    {
        std::cout << "Filling holes..." << std::endl;
        size_t num_borders = 0, num_holes = 0;
        handle_set<P3Halfedge> processed;
        for (auto he = P->halfedges_begin(), end = P->halfedges_end(); he != end; ++he)
        {
            num_borders += he->is_border();
            if (he->is_border() && !processed.count(he))
            {
                num_holes++;
                P3Halfedge start = he;
                do { processed.insert(he); he = he->next(); } while (start != he);
            }
        }
        std::cout << "  " << num_holes << " holes with " << num_borders << " halfedges needed filling" << std::endl;
        triangulate_holes(P);
    }
}

static void fix_mesh(Polyhedron3* P)
{
    std::cout << "Fixing mesh..." << std::endl;
    for (int i = 0; i < 2; ++i)
    {
        remove_self_intersections(P, i == 1);
        fill_holes(P);
        if (P->is_closed() && !PMP::is_outward_oriented(*P))
        {
            std::cout << "Reversing facet orientations..." << std::endl;
            PMP::reverse_face_orientations(*P);
        }
    }
}

static long long get_int(char* a)
{
    if (!a) { usage("ERROR: missing value for option", 1); }
    return atoll(a);
}

static float get_float(char* a)
{
    if (!a) { usage("ERROR: missing value for option", 1); }
    return atof(a);
}

int main(int argc, char** argv)
{
    if (argc == 1) { usage(); }
    if (argc < 4) { usage("ERROR: invalid number of arguments", 1); }

    // Read polyhedron
    Polyhedron3* P;
    std::cout << "Reading mesh..." << std::endl;
    try { P = read_mesh(argv[1]); }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl << std::endl;
        usage("ERROR: unable to read mesh from input file", 2);
    }
    CGAL::set_halfedgeds_items_id(*P);
    show_info(P);
    
    // Make sure the faces are facing outward
    if (P->is_closed() && !PMP::is_outward_oriented(*P))
    {
        std::cout << "Reversing facet orientations..." << std::endl;
        PMP::reverse_face_orientations(*P);
    }
    
    // Refine the mesh
    int i;
    for (i = 2; i < argc; ++i)
    {
        if      (streq(argv[i], "-i")) { isotropic_remeshing_fraction(P, get_float(argv[++i])); }
        else if (streq(argv[i], "-I")) { isotropic_remeshing_length(P, get_float(argv[++i])); }
        else if (streq(argv[i], "-l")) { loop_subdivisions(P, get_int(argv[++i])); }
        else if (streq(argv[i], "-e")) { edge_collapse_fraction(P, get_float(argv[++i])); }
        else if (streq(argv[i], "-E")) { edge_collapse_count(P, get_int(argv[++i])); }
        else if (streq(argv[i], "-f")) { fix_mesh(P); }
        else { break; }
        CGAL::set_halfedgeds_items_id(*P);
        show_info(P);
    }
    if (i < argc-1) { usage("ERROR: invalid arguments", 1); }
    
    // Write polyhedron
    std::cout << "Saving mesh..." << std::endl;
    char* output = argc == i ? strcat(strcpy(new char[strlen(argv[1]) + 16], argv[1]), "_refined.off") : argv[argc-1];
    try { write_mesh(P, output); }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl << std::endl;
        usage("ERROR: unable to write mesh to output file", 7);
    }
    if (argc == i) { delete[] output; }

    return 0;
}
