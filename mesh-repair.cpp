///////////////////////////////////////////////////////////////////////////////
// Simple program to repair a mesh by ensuring all facets are triangular and
// oriented outwards, non-manifold points are corrected, and holes are filed.
///////////////////////////////////////////////////////////////////////////////

#include "GeometryTypes.hpp"
#include "Polyhedron3Utils.hpp"
#include "IO.hpp"
#include "Strings.hpp"
#include "PolyBuilder.hpp"

#include <string.h>
#include <iostream>
#include <vector>

#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/OFF_reader.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
namespace PMP = CGAL::Polygon_mesh_processing;

static void usage(const char* err = nullptr, int exit_code=0)
{
    if (err) { std::cerr << err << std::endl << std::endl; }

    std::cerr << "Usage: mesh-repair [in_file] [optional:out_file]" << std::endl << std::endl;
    std::cerr << "in_file  file to read mesh from" << std::endl;
    std::cerr << "out_file file to save to, defaults to the '${in_file}_repaired.off'" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << mesh_file_usage << std::endl;
    //std::cerr << mesh_file_read_usage << std::endl;
    std::cerr << "Read Options:" << std::endl;
    std::cerr << "  None supported, file is read at a very basic level" << std::endl;
    std::cerr << "  All objects/groups are merged into one" << std::endl;
    std::cerr << mesh_file_write_usage << std::endl;
    exit(exit_code);
}

static void read_mesh(const char* input_fn, std::vector<Point3>& pts, std::vector<std::vector<size_t>>& polygons)
{
    std::cout << "Reading mesh..." << std::endl;
    const char* ext = strrchr(input_fn, '.');
    bool off_file;
    if (ext == nullptr || !((off_file = strieq(ext, ".off")) || strieq(ext, ".obj")))
    {
        usage("ERROR: unknown file extension", 2);
    }
    std::ifstream input(input_fn);
    if (!input)
    {
        usage("ERROR: unable to open input file", 2);
    }
    if (!(off_file ? CGAL::read_OFF(input, pts, polygons) : CGAL::read_OBJ(input, pts, polygons)))
    {
        usage("ERROR: invalid input file format", 2);
    }
    input.close();
}

struct triangle_output
{
    const int a, b, c;
    triangle_output(int a, int b, int c) : a(a), b(b), c(c) { }
};
struct triangle_output_iterator : std::iterator<std::output_iterator_tag, triangle_output>
{
    std::vector<std::vector<size_t>>& polygons;
    const std::vector<size_t>& verts;
    triangle_output_iterator(std::vector<std::vector<size_t>>& polygons, const std::vector<size_t>& verts) : polygons(polygons), verts(verts) { }
    void operator=(triangle_output const& t) { polygons.push_back(std::vector<size_t>({{verts[t.a], verts[t.b], verts[t.c]}})); }
    triangle_output_iterator& operator++() { return *this; }
    triangle_output_iterator operator++(int) { return *this; }
    triangle_output_iterator& operator*() { return *this; }
};
static void triangulate_polygons(std::vector<Point3>& pts, std::vector<std::vector<size_t>>& polygons)
{
    std::cout << "Triangulating facets..." << std::endl;
    auto get_pt = [&pts] (size_t pi) -> const Point3& { return pts[pi]; };
    std::vector<std::vector<size_t>> polys_to_add;
    size_t n_triangulated = 0;
    size_t n_removed = 0;
    for (auto itr = polygons.begin(); itr != polygons.end(); )
    {
        if (itr->size() < 3) { n_removed++; }
        else if (itr->size() == 3)
        {
            const size_t i = (*itr)[0], j = (*itr)[1], k = (*itr)[2];
            const Point3 &p = pts[i], &q = pts[j], &r = pts[k];
            if (i != j && j != k && k != i && p != q && q != r && r != p) { ++itr; continue; }
            n_removed++;
        }
        else
        {
            // TODO: check for simple polygon?
            PMP::triangulate_hole_polyline(boost::adaptors::transform(*itr, get_pt),
                                           triangle_output_iterator(polys_to_add, *itr));
            n_triangulated++;
        }
        itr = polygons.erase(itr);
    }
    if (polys_to_add.size())
    {
        polygons.insert(polygons.begin(), polys_to_add.begin(), polys_to_add.end());
        std::cout << "  " << n_triangulated << " facets triangulated to " << polys_to_add.size() << " triangular facets" << std::endl;
    }
    if (n_removed) { std::cout << "  " << n_removed << " degenerate facets removed" << std::endl; }
}

static Polyhedron3* create_mesh(std::vector<Point3>& pts, std::vector<std::vector<size_t>>& polygons)
{
    std::cout << "Creating polyhedral mesh..." << std::endl;
    if (!PMP::orient_polygon_soup(pts, polygons))
    {
        std::cerr << "  Warning: unable to consistently orient facets" << std::endl;
    }
    if (!PMP::is_polygon_soup_a_polygon_mesh(polygons))
    {
        usage("ERROR: unable to create polyhedral mesh", 3);
    }
    Polyhedron3* P = new Polyhedron3();
    PMP::polygon_soup_to_polygon_mesh(pts, polygons, *P);
    return P;
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
        std::cout << "  " << num_holes << " holes with " << num_borders << " halfedges need filling" << std::endl;
        triangulate_holes(P);
    }
}

static void reorient(Polyhedron3* P)
{
    if (!PMP::is_outward_oriented(*P))
    {
        std::cout << "Reversing faces..." << std::endl;
        PMP::reverse_face_orientations(*P);
    }
}

int main(int argc, char** argv)
{
    if (argc == 1) { usage(); }
    if (argc > 3) { usage("ERROR: invalid number of arguments", 1); }

    // Read polyhedron
    std::vector<Point3> pts;
    std::vector<std::vector<size_t>> polygons;
    read_mesh(argv[1], pts, polygons);
    std::cout << "  Input file contains " << pts.size() << " points and " << polygons.size() << " polygons" << std::endl;
    
    // Triangulate polygons
    triangulate_polygons(pts, polygons);
    
    // Create polyhedral mesh from polygons and points
    Polyhedron3* P = create_mesh(pts, polygons);
    
    // Remove self intersections
    remove_self_intersections(P);

    // Fill any holes if necessary
    fill_holes(P);

    // Reorient mesh if necessary
    if (P->is_closed()) { reorient(P); }
    else { std::cerr << "Warning: mesh not closed even after filling holes" << std::endl; }

    // Now try to remove self intersections aggressively
    remove_self_intersections(P, true);
    
    // One final hole filling
    fill_holes(P);
    
    // Check mesh
    if (!P->is_valid())         { std::cerr << "Warning: mesh is not valid" << std::endl; }
    if (!is_not_degenerate(P))  { std::cerr << "Warning: mesh is degenerate" << std::endl; }
    if (!P->is_closed())        { std::cerr << "Warning: mesh is not closed" << std::endl; }
    else if (!PMP::is_outward_oriented(*P)) { std::cerr << "Warning: mesh is not outward oriented" << std::endl; }
    if (!P->is_pure_triangle()) { std::cerr << "Warning: mesh is not pure triangle" << std::endl; }
    if (PMP::does_self_intersect(*P))
    {
        std::cerr << "Warning: mesh is self-intersecting" << std::endl;
        std::vector<std::pair<P3Facet, P3Facet>> out;
        PMP::self_intersections(*P, std::back_inserter(out));
        std::cout << "  " << out.size() << " self intersections" << std::endl;
        for (auto& si : out)
        {
            std::cout << "    " << facet_to_triangle3(si.first) << " & " << facet_to_triangle3(si.second) << std::endl;
        }
    }
    //else if (!PMP::does_bound_a_volume(*P)) { std::cerr << "Warning: mesh does not bound a volume" << std::endl; }
    //if (!is_single_component(P)) { std::cerr << "Warning: mesh is not single connected component" << std::endl; }
    
    
    // Write polyhedron
    std::cout << "Saving mesh..." << std::endl;
    char* output = argc == 2 ? strcat(strcpy(new char[strlen(argv[1]) + 16], argv[1]), "_repaired.off") : argv[2];
    try { write_mesh(P, output); }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl << std::endl;
        usage("ERROR: unable to write mesh to output file");
        return 7;
    }
    if (argc == 2) { delete[] output; }

    // Output information about the output mesh
    std::cout << "Output polyhedral mesh contained " << P->size_of_vertices() << " vertices, " << P->size_of_facets() << " facets, and " << P->size_of_halfedges() << " halfedges" << std::endl;
    typedef boost::property_map<Polyhedron3, boost::face_index_t>::const_type P3_facet_index_map_t;
    typedef boost::vector_property_map<size_t, P3_facet_index_map_t> P3_facet_int_map;
    P3_facet_int_map component_id(get(boost::face_index, *P));
    size_t n_components = PMP::connected_components(*P, component_id);
    std::cout << "  Connected components: " << n_components << std::endl;
    for (size_t i = 0; i < n_components; ++i)
    {
        size_t count = 0;
        for (P3Facet f = P->facets_begin(), end = P->facets_end(); f != end; ++f) { count += component_id[f] == i; }
        std::cout << "    " << count << " facets in component #" << i <<  std::endl;
    }
    
    return 0;
}
