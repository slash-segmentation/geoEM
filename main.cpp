#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include "Skeleton.hpp"
#include "Polyhedron3Utils.hpp"
#include "Slice.hpp"

#include "IO.hpp"
#include "IO_OBJ.hpp"
#include "PolyBuilder.hpp"
#include "Segments2Cylinders.hpp"
#include "Points2Spheres.hpp"

#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include <utility>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#ifdef CREATE_GUI
#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt.h>
#endif

#include <CGAL/boost/graph/properties.h>

#include <CGAL/subdivision_method_3.h>

#include <CGAL/polyhedron_cut_plane_3.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
namespace PMP = CGAL::Polygon_mesh_processing;


typedef handle_map<P3CVertex, S3VertexDesc> VertexMap;
typedef std::unordered_map<Point3, S3VertexDesc, boost::hash<Point3>> PointMap;
// TODO: should this use a Segment3 instead of a pair of Point3s?
typedef std::pair<Point3, Point3> EdgeMapKey;
typedef std::unordered_map<EdgeMapKey, P3CHalfedge, boost::hash<EdgeMapKey>> EdgeMap;

// Polyhedron Property Maps
typedef boost::property_map<Polyhedron3, boost::face_index_t>::const_type P3_facet_index_map_t;
typedef boost::vector_property_map<size_t, P3_facet_index_map_t> P3_facet_int_map;

// Skeleton Property Maps
typedef boost::property_map<Skeleton3, boost::vertex_index_t>::const_type S3_vert_index_map_t;
typedef boost::vector_property_map<Kernel::FT, S3_vert_index_map_t> S3_vert_number_map;
typedef boost::vector_property_map<double, S3_vert_index_map_t> S3_vert_double_map;
typedef boost::vector_property_map<bool, S3_vert_index_map_t> S3_vert_bool_map;
typedef boost::vector_property_map<S3VertexDesc, S3_vert_index_map_t> S3_vert_vert_map;

// Self-Intersection Results
typedef std::vector<std::pair<P3CFacet, P3CFacet>> self_intersections;

#include <list>
#include <CGAL/AABB_segment_primitive.h>

typedef std::vector<Segment3> SkelSearchStorage;
typedef CGAL::AABB_tree<CGAL::AABB_traits<Kernel, CGAL::AABB_segment_primitive<Kernel, SkelSearchStorage::iterator>>> SkelSearchTree; // TODO: support skeleton edges directly instead of segments


std::vector<Point3> get_bp_mesh_points(const Polyhedron3* P, const Skeleton3* S)
{
    // Gets all of the mesh vertex points associated with the branch points in a skeleton
    std::vector<Point3> pts;
    BOOST_FOREACH (S3VertexDesc v, vertices(*S))
    {
        if (degree(v, *S) > 2)
        {
            // Branch point
            for (auto mesh_v : (*S)[v].vertices)
            {
                // For my Polyhedron3 can just do mesh_v->point() but this solution is more general
                pts.push_back(get(CGAL::vertex_point, *P, mesh_v));
            }
        }
    }
    return pts;
}

class Polyline3 : public std::vector<Point3>
{
    bool _open;
    typedef std::vector<Point3> Base;
public:
    Polyline3(bool is_open = false) : Base(), _open(is_open) {}
    Polyline3(const Polyline3& p) : Base(p.begin(), p.end()), _open(p._open) {}
    template <class InputIterator> Polyline3(InputIterator first, InputIterator last, bool is_open = false) : Base(first, last), _open(is_open) {}
};

std::vector<Polyhedron3*> get_bp_surface(const Polyhedron3* P, const Skeleton3* S, bool rings=false)
{
    // This creates a mesh that contains all of the facets that contributed to a branch point.
    // This does end up dropping important vertices where there were insufficient vertices to make
    // a facet. It is possible that what is needed is some way to create "degenerate facets" in the
    // builder above that connect neighboring vertices that are otherwise being dropped. That may
    // still leave some vertices that are completely isolated however.
    std::vector<Polyhedron3*> segs;
    BOOST_FOREACH (S3VertexDesc v, vertices(*S))
    {
        if (degree(v, *S) > 2)
        {
            P3CVertexSet verts;

            // Add all vertices that contributed to the branch point itself
            verts.insert((*S)[v].vertices.begin(), (*S)[v].vertices.end());

            // Also add all vertices that contributed to the neighbors of the branch point
            BOOST_FOREACH (S3EdgeDesc e, out_edges(v, *S))
            {
                S3VertexDesc u = opposite(*S, e, v);
                verts.insert((*S)[u].vertices.begin(), (*S)[u].vertices.end());
            }

            // TODO: what may be better is instead of adding the vertices that contributed to the
            // neighboring skeletal points at this point add any facet that has a vertex that contributed
            // to the branch point directly.

            // Build the mesh that represent the entire "branch point area"
            Polyhedron3* seg = extract_mesh_from_vertices(verts);
            if (!rings)
            {
                // Just the meshes of the branch point areas
                // These are not capped, but that could likely be added
                segs.push_back(seg);
            }
            else
            {
                // Build up the borders of the branch point areas
                seg->normalize_border();
                PolyBuilder<Point3, Polyline3> pb;
                for (auto i = seg->border_edges_begin(), end = seg->edges_end(); i != end; ++i)
                {
                    pb.add_seg(i->vertex()->point(), i->opposite()->vertex()->point());
                }
                std::vector<Polyline3> polylines = pb.finish_up();

                for (auto i = polylines.begin(), i_end = polylines.end(); i != i_end; ++i)
                {
                    auto j = i->begin(), j_end = i->end();
                    Point3 pt_end = *j++;
                    std::vector<Segment3> ring;
                    ring.reserve(i->size());
                    while (j != j_end)
                    {
                        Point3 pt = pt_end; pt_end = *j++;
                        ring.push_back(Segment3(pt, pt_end));
                    }
                    segs.push_back(segments2cylinders_merged(ring, 2.5, 8));
                }
            }
        }
    }
    return segs;
}


template <class Graph>
bool is_cyclic(const Graph& g)
{
    // Assumes that the graph is a single connected component
    return num_vertices(g) - 1 != num_edges(g);
}

S3_vert_double_map calc_distances_to_soma(const Polyhedron3* P, const Skeleton3* S)
{
    // Find the endpoint that represents the "soma" (which has largest normalized surface area)
    bool first = true;
    Kernel::FT max = 0;
    S3VertexDesc soma = 0; // 0 is to suppress warnings about possible uninitialized variables
    BOOST_FOREACH (S3VertexDesc v, vertices(*S))
    {
        if (degree(v, *S) == 1)
        {
            // Endpoint found, calculate its normalized surface area
            S3VertexDesc u = next_vertex(*S, v);
            Point3 a = (*S)[v].point, b = (*S)[u].point;
            P3CVertexSet vs((*S)[v].vertices.begin(), (*S)[v].vertices.end());
            vs.insert((*S)[u].vertices.begin(), (*S)[u].vertices.end());
            Polyhedron3* seg = extract_mesh_from_vertices(vs);
            Kernel::FT sa = surface_area(seg)/distance(a, b);
            delete seg;
            // Check if it is a new maximum
            if (first || sa > max) { soma = v; max = sa; first = false; }
        }
    }

    // Calculate the distance for each skeleton point
    S3_vert_double_map distances(get(boost::vertex_index, *S));
    S3_vert_bool_map discovered(get(boost::vertex_index, *S));
    std::vector<std::pair<S3VertexDesc, double>> stack;
    stack.push_back(std::make_pair(soma, 0.0));
    discovered[soma] = true;
    while (!stack.empty())
    {
        S3VertexDesc v = stack.back().first;
        double dist = stack.back().second;
        stack.pop_back();
        Point3 p = (*S)[v].point;
        distances[v] = dist;
        BOOST_FOREACH (S3EdgeDesc e, out_edges(v, *S))
        {
            S3VertexDesc u = opposite(*S, e, v);
            if (discovered[u]) { continue; }
            stack.push_back(std::make_pair(u, dist + distance(p, (*S)[u].point)));
            discovered[u] = true;
        }
    }

    return distances;
}

S3_vert_vert_map calc_next_bps(const Polyhedron3* P, const Skeleton3* S, const S3_vert_double_map& dist_to_soma)
{
    S3_vert_vert_map next_bp(get(boost::vertex_index, *S));
    S3_vert_bool_map valid(get(boost::vertex_index, *S));

    // Recursive function
    std::function<S3VertexDesc(S3VertexDesc)> f = [&] (S3VertexDesc v) -> S3VertexDesc
    {
        if (valid[v]) { return next_bp[v]; }
        int deg = degree(v, *S);
        if ((deg == 1 && dist_to_soma[v] == 0.0) || deg > 2)
        {
            // A branch point or the soma itself
            valid[v] = true;
            next_bp[v] = v;
            return v;
        }

        // Travel towards the soma/branch-point
        S3VertexDesc nxt;
        if (deg == 1) { nxt = next_vertex(*S, v); }
        else
        {
            // Get the two neighboring vertices
            auto e_itr = out_edges(v, *S).first;
            S3VertexDesc v1 = opposite(*S, *e_itr++, v);
            S3VertexDesc v2 = opposite(*S, *e_itr, v);
            // Select the one closer to the soma
            nxt = dist_to_soma[v1] < dist_to_soma[v2] ? v1 : v2;
        }
        // Recursively call the function
        S3VertexDesc u = f(nxt);
        valid[v] = true;
        next_bp[v] = u;
        return u;
    };

    // Apply function to every vertex
    BOOST_FOREACH (S3VertexDesc v, vertices(*S)) { f(v); }
    return next_bp;
}

void write_obj_cmap_segments(std::ofstream& f, std::vector<Polyhedron3*>& segs, std::vector<double>& values,
                             size_t& off, const std::string& cmap = "hot")
{
    double max = *std::max_element(values.begin(), values.end());
    double min = *std::min_element(values.begin(), values.end());
    std::vector<std::string> colors;
    for (auto value : values)
    {
        char num[4];
        sprintf(num, "%d", (int)((value - min) / (max - min) * 255 + 0.5));
        colors.push_back(cmap + num);
    }
    write_obj_file(f, segs, off, false, cmap + ".mtl", colors);
}


int main(int argc, char **argv)
{
    CGAL::set_pretty_mode(std::cout);

    ///////////////////////////////////////////////////////////////////////////
    // Set the file to read
    ///////////////////////////////////////////////////////////////////////////
    bool assume_good = false;
    //std::string filename = "example-data/other/cube-triangles.off";
    //std::string filename = "example-data/other/cube-quads.off";
    //std::string filename = "example-data/other/elephant.off";
    //std::string filename = "example-data/trunc_cone.off";
    //std::string filename = "example-data/small.off";
    //std::string filename = "example-data/big.off";
    std::string filename = "example-data/big-nn.off"; // no nucleus
    //std::string filename = "example-data/big_simplified.off";


    ///////////////////////////////////////////////////////////////////////////
    // Settings
    ///////////////////////////////////////////////////////////////////////////
    double remesh_size = 0.75; // perform isotropic remeshing to make the size of triangles more equal
    int loop_subdivisions = 1; // add extra vertices to the input mesh
    const char* output_obj = "output.obj"; // the output OBJ file
    //bool output_bp_surface = true; // output the branch point surface
    const char* output_skel = "skel.cgal"; // the output CGAL file for the skeleton points


    ///////////////////////////////////////////////////////////////////////////
    // Read in the 3D polyhedral mesh
    ///////////////////////////////////////////////////////////////////////////
    Polyhedron3 *P;
    {
        std::cout << "Reading mesh..." << std::endl;
        boost::timer::auto_cpu_timer t;
        try
        {
            P = read_mesh(filename, assume_good); // TODO: delete P; somewhere...
        }
        catch (std::invalid_argument& err) { std::cerr << err.what() << std::endl; return -1; }
        CGAL::set_halfedgeds_items_id(*P);
        calculate_facet_planes(P);
        if (!assume_good && PMP::connected_components(*P, P3_facet_int_map(get(boost::face_index, *P))) != 1)
        {
            std::cerr << "ERROR: model is not a single connected component" << std::endl;
            return -1;
        }
    }
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // Refine the mesh
    ///////////////////////////////////////////////////////////////////////////
    if (remesh_size || loop_subdivisions > 0)
    {
        std::cout << "Refining mesh..." << std::endl;
        boost::timer::auto_cpu_timer t;
        Kernel::FT avg_edge_len = avg_edge_length(P);
        std::cout << "    Original: " << num_faces(*P) << " faces and an average edge length of " << avg_edge_len << std::endl;

        if (remesh_size)
        {
            PMP::isotropic_remeshing(faces(*P), CGAL::to_double(avg_edge_len*remesh_size), *P, PMP::parameters::number_of_iterations(3));
            std::cout << "    After remeshing: " << num_faces(*P) << " faces and an average edge length of " << avg_edge_length(P) << std::endl;
        }

        if (loop_subdivisions > 0)
        {
            CGAL::Subdivision_method_3::Loop_subdivision(*P, loop_subdivisions);
            std::cout << "    After subdivision: " << num_faces(*P) << " faces and an average edge length of " << avg_edge_length(P) << std::endl;
        }

        if (!P->is_valid())         { std::cerr << "WARNING: mesh is not valid after refining" << std::endl; }
        if (!P->is_closed())        { std::cerr << "WARNING: mesh is not closed after refining" << std::endl; }
        if (!is_not_degenerate(P))  { std::cerr << "WARNING: mesh is degenerate after refining" << std::endl; }
        if (PMP::does_self_intersect(*P))  { std::cerr << "WARNING: mesh is self-intersecting after refining" << std::endl; }
        if (!P->is_pure_triangle()) { std::cerr << "WARNING: mesh is not pure triangle after refining" << std::endl; }
    }
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // Add facet "vertices"
    ///////////////////////////////////////////////////////////////////////////
    // The middle of each facet will have a vertex added to it and that is the
    // vertex we will care about when looking at associations between the
    // skeleton and the mesh.
    /*P3CVertexSet facet_verts;
    {
        std::cout << "Adding facet vertices..." << std::endl;
        boost::timer::auto_cpu_timer t;
        // Get the current set of facets
        std::vector<P3Facet> facets;
        facets.reserve(P->size_of_facets());
        for (auto f = P->facets_begin(), f_end = P->facets_end(); f != f_end; ++f)
        {
            facets.push_back(f);
        }
        // Go through and add split each current facet
        for (auto f : facets)
        {
            auto he = f->facet_begin();
            auto p = (he++)->vertex()->point(), q = (he++)->vertex()->point(), r = (he)->vertex()->point();
            auto v = P->create_center_vertex(f->halfedge())->vertex();
            v->point() = CGAL::centroid(p, q, r);
            facet_verts.insert(v);
            // The first halfedge of each sub-facet points to the "facet vertex"
            FOR_EDGES_AROUND_VERTEX(v, e) { e->facet()->set_halfedge(e); }
        }
        CGAL::set_halfedgeds_items_id(*P); // after adding lots of vertices/facets all of the IDs probably need fixing
        calculate_facet_planes(P);
    }
    std::cout << std::endl;*/

    ///////////////////////////////////////////////////////////////////////////
    // Construct the skeleton
    ///////////////////////////////////////////////////////////////////////////
    Skeleton3* S;
    {
        std::cout << "Constructing skeleton..." << std::endl;
        boost::timer::auto_cpu_timer t;
        S = construct_skeleton(P, 0.5); // TODO: delete S; somewhere...
    }
    std::cout << std::endl;
    // Save the skeleton in an easy-to-use format
    // NOTE: This loses the surface mesh vertices associated with every point
    write_skeleton(S, output_skel, 10);
    // read with: S = read_cg(output_skel);
    if (is_cyclic(*S)) { std::cerr << "ERROR: Skeleton is cyclic!" << std::endl; return -1; }


    ///////////////////////////////////////////////////////////////////////////
    // Map mesh vertices to skeleton vertices
    ///////////////////////////////////////////////////////////////////////////
    // There is a many to one relationship between the vertices of the mesh and
    // the skeleton. We can get all mesh vertices from a skeleton vertex with
    // (*S)[v].vertices. This sets up the reverse direction, going from a mesh
    // vertex to its single associated skeleton vertex.
    // Additionally we set of maps that go from points (not vertices) on the
    // mesh to skeleton vertices and pairs of points on the mesh to halfedges
    // on the mesh. The edges go from source to target.
    /*VertexMap vmap;
    PointMap pmap;
    EdgeMap emap;
    {
        std::cout << "Creating mesh to skeleton mapping..." << std::endl;
        boost::timer::auto_cpu_timer t;
        vmap.reserve(P->size_of_vertices());
        pmap.reserve(P->size_of_vertices());
        emap.reserve(P->size_of_halfedges());
        BOOST_FOREACH (S3VertexDesc v, vertices(*S))
        {
            auto& verts = (*S)[v].vertices;
            if (verts.size() == 0) { std::cerr << "Warning: a skeleton vertex has no mesh vertex associated with it" << std::endl; }
            for (P3CVertex u : verts)
            {
                auto itr = vmap.find(u);
                if (itr != vmap.end()) { std::cerr << "ERROR: Skeleton mapping not many-to-one!" << std::endl; return -1; }
                vmap.insert(itr, std::make_pair(u, v));
                pmap.insert({{u->point(), v}});
            }
        }
        if (vmap.size() != P->size_of_vertices()) { std::cerr << "ERROR: Not all mesh vertices mapped to skeleton!" << std::endl; return -1; }
        for (auto he = P->halfedges_begin(), he_end = P->halfedges_end(); he != he_end; ++he)
        {
            emap.insert({{std::make_pair(he->opposite()->vertex()->point(), he->vertex()->point()), he}});
        }
        if (emap.size() != P->size_of_halfedges()) { std::cerr << "ERROR: Not all mesh edges mapped handles!" << std::endl; return -1; }
    }
    std::cout << std::endl;*/


    ///////////////////////////////////////////////////////////////////////////
    // Create the skeleton search tree
    ///////////////////////////////////////////////////////////////////////////
    /*SkelSearchTree* skel_search;
    SkelSearchStorage skel_search_storage;
    {
        std::cout << "Creating skeleton search tree..." << std::endl;
        boost::timer::auto_cpu_timer t;
        skel_search_storage.reserve(num_edges(*S));
        BOOST_FOREACH (auto e, edges(*S))
        {
            skel_search_storage.push_back(Segment3((*S)[source(e, *S)].point,
                                                   (*S)[target(e, *S)].point));
        }
        skel_search = new SkelSearchTree(skel_search_storage.begin(), skel_search_storage.end());
        skel_search->accelerate_distance_queries();
    }
    std::cout << std::endl;*/


    ///////////////////////////////////////////////////////////////////////////
    // Save the model as an OBJ file
    ///////////////////////////////////////////////////////////////////////////
    std::ofstream f(output_obj);
    f << std::setprecision(10);
    size_t off = 0;

    // The skeleton as a series of black cylinders
    std::vector<Polyhedron3*> skeleton_cyl;
    skeleton_enum_branches_as_pts(S, [&skeleton_cyl] (const std::vector<Point3>& pts) mutable
    {
        std::vector<Segment3> branch;
        branch.reserve(pts.size());
        auto i = pts.begin(), end = pts.end();
        Point3 pt_end = *i++;
        while (i != end)
        {
            Point3 pt = pt_end; pt_end = *i++;
            branch.push_back(Segment3(pt, pt_end));
        }
        skeleton_cyl.push_back(segments2cylinders_merged(branch, 2.5, 8));
    });
    write_obj_file(f, skeleton_cyl, off, false, "colors.mtl", std::vector<std::string>{ "Black" });
    for (auto itr = skeleton_cyl.begin(), end = skeleton_cyl.end(); itr != end; ++itr) { delete *itr; }

    // Calculate the basic values that follow the skeleton
    S3_vert_double_map dists = calc_distances_to_soma(P, S);
    S3_vert_vert_map next_bps = calc_next_bps(P, S, dists);

    // The mesh as a series of colored segments
    auto slices = slice(5, P, S);
    auto capped_slices = cap_slices(slices);

    std::vector<Polyhedron3*> segs;
    std::vector<double> values;
    for (auto& slc : capped_slices)
    {
        // Volume:
        //double val = CGAL::to_double(volume(slc->mesh()));

        // Normalized volume:
        //double val = CGAL::to_double(volume(slc->mesh())/slc->length());

        // Distance to soma:
        double val = 0.0;
        for (auto sv : slc->skeleton_vertices()) { val += dists[sv]; }
        val /= slc->skeleton_vertices().size();

        // Distance to branch point:
        S3VertexDesc near = *slc->skeleton_vertices().begin();
        for (auto sv : slc->skeleton_vertices()) { if (dists[sv] < dists[near]) { near = sv; }  }
        val -= dists[next_bps[near]];

        values.push_back(val);
        segs.push_back(slc->mesh());
    }
    write_obj_cmap_segments(f, segs, values, off);

    // The segments in the mesh in a myriad of colors
    /*std::vector<std::string> colors { // missing Black which is used for skeleton
        "Maroon", "Red", "Orange-Red", "Orange", "Gold", "Yellow", "Yellow-Green", "Lime",
        "Green", "Spring-Green", "Cyan", "Dark-Cyan", "Dark-Turquoise", "Deep-Sky-Blue", "Navy",
        "Medium-Blue", "Blue", "Royal-Blue", "Blue-Violet", "Indigo", "Purple", "Violet",
        "Magenta", "Hot-Pink", "Pink", "Brown", "Sienna", "Dark-Gray", "Gray", "Silver", "White",
    };*/

    // Draw the branch point surface
    /*if (output_bp_surface)
    {
        //std::vector<Polyhedron3*> s = points2spheres(get_bp_mesh_points(P, S), 5.0, 3); // points
        //std::vector<Polyhedron3*> s = get_bp_surface(P, S, false); // surface
        //std::vector<Polyhedron3*> s = get_bp_surface(P, S, true); // rings
        write_obj_file(f, s, off, false, "", std::vector<std::string>{ "Black" });
        for (auto x : s) { delete x; }
    }*/

    // Close the OBJ file
    f.close();


    std::cout << "There are " << slices.size() << " slices" << std::endl;
    std::cout << "The original mesh contains:" << std::endl;
    std::cout << "    " << P->size_of_facets() << " facets" << std::endl;
    std::cout << "    " << P->size_of_halfedges() << " halfedges" << std::endl;
    std::cout << "    " << P->size_of_vertices() << " vertices" << std::endl;
    size_t total_facets = 0;
    size_t total_halfedges = 0;
    size_t total_vertices = 0;
    for (auto& slc : slices)
    {
        total_facets += slc->mesh()->size_of_facets();
        total_halfedges += slc->mesh()->size_of_halfedges();
        total_vertices += slc->mesh()->size_of_vertices();
    }
    std::cout << "The slices contain:" << std::endl;
    std::cout << "    " << total_facets << " facets" << std::endl;
    std::cout << "    " << total_halfedges << " halfedges" << std::endl;
    std::cout << "    " << total_vertices << " vertices" << std::endl;
    std::cout << "The volume of the entire mesh is: " << volume(P) << std::endl;
    Kernel::FT total_volume = 0;
    for (auto& slc : capped_slices)
    {
        total_volume += volume(slc->mesh());
    }
    std::cout << "The volume of the sum of the capped slices is: " << total_volume << std::endl;
    std::cout << "The surface area of the entire mesh is: " << surface_area(P) << std::endl;
    Kernel::FT total_sa = 0;
    for (auto& slc : slices)
    {
        total_sa += surface_area(slc->mesh());
    }
    std::cout << "The surface area of the sum of the slices is: " << total_sa << std::endl;


#ifdef CREATE_GUI
    QApplication app(argc, argv);
    app.setApplicationName("Cross Section");

    // Import resources from libCGALQt.
    // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
    CGAL_QT_INIT_RESOURCES;

    /////////////////////////////////////////////////////////////////////////////
    //// Setup the viewer
    /////////////////////////////////////////////////////////////////////////////
    Viewer viewer;
    viewer.setWindowTitle("Cross Section");
    viewer.set_polyhedron(P);
    viewer.set_skeleton(S);
    //viewer.set_intersection(intersection);
    //viewer.set_point(pt, 100);
    viewer.show();

    /////////////////////////////////////////////////////////////////////////////
    //// Run the visual part of the program
    /////////////////////////////////////////////////////////////////////////////
    QStringList args = app.arguments();
    args.removeAt(0);
    if (!args.empty() && args[0] == "--use-meta")
    {
        viewer.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
        args.removeAt(0);
    }
    return app.exec();
#else
    return 0;
#endif

}
