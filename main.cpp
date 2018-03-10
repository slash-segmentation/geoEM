#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include "IO.hpp"
#include "IO_OBJ.hpp"
#include "PolyBuilder.hpp"
#include "Intersection.hpp"
#include "Skeleton.hpp"
#include "Segments2Cylinders.hpp"
#include "Points2Spheres.hpp"
#include "Polyhedron3Utils.hpp"

#include <boost/timer/timer.hpp>

#include <tuple>
#include <iostream>

#ifdef CREATE_GUI
#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt.h>
#endif

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/subdivision_method_3.h>

typedef Skeleton3::vertex_descriptor S3VertexDesc
typedef handle_map<Polyhedron3::Vertex_const_handle, S3VertexDesc> VertexMap;

std::vector<Point3> get_bp_mesh_points(const Polyhedron3* P, const Skeleton3* S)
{
    // Gets all of the mesh vertex points associated with the branch points in a skeleton
    std::vector<Point3> pts;
    BOOST_FOREACH(S3VertexDesc v, vertices(*S))
    {
        if (degree(v, *S) > 2)
        {
            // Branch point
            BOOST_FOREACH(auto mesh_v, (*S)[v].vertices)
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
    BOOST_FOREACH(S3VertexDesc v, vertices(*S))
    {
        if (degree(v, *S) > 2)
        {
            handle_set<Polyhedron3::Vertex_const_handle> verts;

            // Add all vertices that contributed to the branch point itself
            verts.insert((*S)[v].vertices.begin(), (*S)[v].vertices.end());

            // Also add all vertices that contributed to the neighbors of the branch point
            BOOST_FOREACH(Skeleton3::edge_descriptor e, out_edges(v, *S))
            {
                S3VertexDesc u = opposite(*S, e, v);
                verts.insert((*S)[u].vertices.begin(), (*S)[u].vertices.end());
            }

            // TODO: what may be better is instead of adding the vertices that contributed to the
            // neighboring skeletal points at this point add any facet that has a vertex that contributed
            // to the branch point directly.

            // Build the mesh that represent the entire "branch point area"
            Polyhedron3* seg = extract_mesh(verts);
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

Polyhedron3* extract_mesh_capped(const Point3& a, const Point3& b,
                                 const std::vector<Polyhedron3::Vertex_handle>& va,
                                 const std::vector<Polyhedron3::Vertex_handle>& vb)
{
    // Extract the mesh patch
    handle_set<Polyhedron3::Vertex_const_handle> vs(va.begin(), va.end());
    vs.insert(vb.begin(), vb.end());
    Polyhedron3* seg = extract_mesh(vs);

    // Cap the ends
    std::unordered_set<Point3, boost::hash<Point3> > a_map;
    BOOST_FOREACH(auto pv, va) { a_map.insert(pv->point()); }
    std::unordered_set<Point3, boost::hash<Point3> > b_map;
    BOOST_FOREACH(auto pv, vb) { b_map.insert(pv->point()); }
    seg->normalize_border();
    for (auto i = seg->border_halfedges_begin(), end = seg->halfedges_end(); i != end; ++i)
    {
        if (i->is_border()) // since we are removing borders we need to check this
        {
            seg->fill_hole(i);
            seg->create_center_vertex(i)->vertex()->point() = a_map.count(i->vertex()->point()) ? a : b;
        }
    }

    // Complete
    return seg;
}

std::vector<std::tuple<S3VertexDesc, S3VertexDesc, Polyhedron3*> > break_into_segments(const Polyhedron3* P, const Skeleton3* S)
{
    std::vector<std::tuple<S3VertexDesc, S3VertexDesc, Polyhedron3*> > segs;
    skeleton_enum_branches(S, [&segs, S] (const std::vector<S3VertexDesc>& verts) mutable
    {
        for (size_t i = 0; i < verts.size() - 1;)
        {
            auto& v = (*S)[verts[i]];
            auto& u = (*S)[verts[i+1]];
            //assert(v.vertices.size() != 0 && u.vertices.size() != 0);
            Polyhedron3* seg = extract_mesh_capped(v.point, u.point, v.vertices, u.vertices);
            segs.push_back(std::make_tuple(verts[i], verts[i+1], seg));
        }
    });
    return segs;
}

template <class Graph>
bool is_cyclic(const Graph& g)
{
    // TODO
    return false;
}

typedef boost::property_map<Skeleton3, boost::vertex_index_t>::const_type Skel_vert_index_map_t;
typedef boost::vector_property_map<Kernel::FT, Skel_vert_index_map_t> Skel_vert_number_map;
typedef boost::vector_property_map<double, Skel_vert_index_map_t> Skel_vert_double_map;
typedef boost::vector_property_map<bool, Skel_vert_index_map_t> Skel_vert_bool_map;
typedef boost::vector_property_map<S3VertexDesc, Skel_vert_index_map_t> Skel_vert_vert_map;

Skel_vert_double_map calc_distances_to_soma(const Polyhedron3* P, const Skeleton3* S)
{
    // Find the endpoint that represents the "soma" (which has largest normalized volume)
    bool first = true;
    Kernel::FT max = 0;
    S3VertexDesc soma = 0; // 0 is to supress warnings about possible uninitialized variables
    BOOST_FOREACH(S3VertexDesc v, vertices(*S))
    {
        if (degree(v, *S) == 1)
        {
            // Endpoint found, calculate its noramlized volume
            S3VertexDesc u = next_vertex(*S, v);
            Point3 a = (*S)[v].point, b = (*S)[u].point;
            Polyhedron3* seg = extract_mesh_capped(a, b, (*S)[v].vertices, (*S)[u].vertices);
            Kernel::FT vol = volume(seg)/distance(a, b);
            delete seg;
            // Check if it is a new maximum
            if (first || vol > max) { soma = v; max = vol; first = false; }
        }
    }
    
    // Calcualte the distance for each skeleton point
    Skel_vert_double_map distances(get(boost::vertex_index, *S));
    Skel_vert_bool_map discovered(get(boost::vertex_index, *S));
    std::vector<std::pair<S3VertexDesc, double> > stack;
    stack.push_back(std::make_pair(soma, 0.0));
    discovered[soma] = true;
    while (!stack.empty())
    {
        S3VertexDesc v = stack.back().first;
        double dist = stack.back().second;
        stack.pop_back();
        Point3 p = (*S)[v].point;
        distances[v] = dist;
        BOOST_FOREACH(Skeleton3::edge_descriptor e, out_edges(v, *S))
        {
            S3VertexDesc u = opposite(*S, e, v);
            if (discovered[u]) { continue; }
            stack.push_back(std::make_pair(u, dist + distance(p, (*S)[u].point)));
            discovered[u] = true;
        }
    }
    
    return distances;
}

Skel_vert_vert_map calc_next_bps(const Polyhedron3* P, const Skeleton3* S, const Skel_vert_double_map& dist_to_soma)
{
    Skel_vert_vert_map next_bp(get(boost::vertex_index, *S));
    Skel_vert_bool_map valid(get(boost::vertex_index, *S));
    
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
    BOOST_FOREACH(S3VertexDesc v, vertices(*S)) { f(v); }
    return next_bp;
}

void write_obj_cmap_segments(std::ofstream& f, std::vector<Polyhedron3*>& segs, std::vector<double>& values,
                             size_t& off = 0, const std::string& cmap = "hot")
{
    double max = *std::max_element(values.begin(), values.end());
    double min = *std::min_element(values.begin(), values.end());
    std::vector<std::string> colors;
    BOOST_FOREACH(auto value, values)
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
    bool assume_good = true;
    //std::string filename = "example-data/other/cube-triangles.off";
    //std::string filename = "example-data/other/cube-quads.off";
    //std::string filename = "example-data/other/elephant.off";
    //std::string filename = "example-data/small.off";
    //std::string filename = "example-data/big.off";
    std::string filename = "example-data/big-nn.off"; // no nucleus
    //std::string filename = "example-data/big_simplified.off";


    ///////////////////////////////////////////////////////////////////////////
    // Settings
    ///////////////////////////////////////////////////////////////////////////
    bool remesh = true; // perform isotropic remeshing to make the size of triangles more equal
    int loop_subdivisions = 2; // add extra vertices to the input mesh
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
        P = read_mesh(filename, assume_good); // TODO: delete P; somewhere...
    }
    std::cout << std::endl;

    
    ///////////////////////////////////////////////////////////////////////////
    // Refine the mesh
    ///////////////////////////////////////////////////////////////////////////
    if (remesh || loop_subdivisions > 0)
    {
        std::cout << "Refining mesh..." << std::endl;
        boost::timer::auto_cpu_timer t;
        Kernel::FT avg_edge_len = avg_edge_length(P);
        std::cout << "    Original: " << num_faces(*P) << " faces and an average edge length of " << avg_edge_len << std::endl;
        
        if (remesh)
        {
            namespace PMP = CGAL::Polygon_mesh_processing;
            std::size_t i = 0;
            for (auto itr = P->facets_begin(), end = P->facets_end(); itr != end; ++itr, ++i) { itr->id() = i; }
            PMP::isotropic_remeshing(faces(*P), avg_edge_len, *P, PMP::parameters::number_of_iterations(3));
            std::cout << "    After remeshing: " << num_faces(*P) << " faces and an average edge length of " << avg_edge_length(P) << std::endl;
        }
        
        if (loop_subdivisions > 0)
        {
            CGAL::Subdivision_method_3::Loop_subdivision(*P, loop_subdivisions);
            std::cout << "    After subdivision: " << num_faces(*P) << " faces and an average edge length of " << avg_edge_length(P) << std::endl;
        }
    }
    std::cout << std::endl;

    
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
    VertexMap vmap;
    {
        std::cout << "Creating mesh to skeleton mapping..." << std::endl;
        boost::timer::auto_cpu_timer t;
        vmap.reserve(P->size_of_vertices());
        BOOST_FOREACH(S3VertexDesc v, vertices(*S))
        {
            auto& verts = (*S)[v].vertices;
            if (verts.size() == 0) { std::cerr << "Warning: a skeleton vertex has no mesh vertex associated with it" << std::endl; }
            BOOST_FOREACH(Polyhedron3::Vertex_const_handle u, verts)
            {
                auto itr = vmap.find(u);
                if (itr != vmap.end()) { std::cerr << "ERROR: Skeleton mapping not many-to-one!" << std::endl; return -1; }
                vmap.insert(itr, make_pair(u, v));
            }
        }
        if (vmap.size() != P->size_of_vertices()) {std::cerr << "ERROR: Not all mesh vertices mapped to skeleton!" << std::endl; return -1; }
    }


    ///////////////////////////////////////////////////////////////////////////
    // Save the model as an OBJ file
    ///////////////////////////////////////////////////////////////////////////
    // Model includes:
    //  * The skeleton (solid black)
    //  * The segments (semi-transparent colored [or white])
    //  * The "branch point surface" (optional, not complete yet)
    std::ofstream f(output_obj);
    f << std::setprecision(10);
    size_t off = 0;

    // The skeleton as a series of cylinders
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

    
    auto segments = break_into_segments(P, S);
    Skel_vert_double_map dists = calc_distances_to_soma(P, S);
    Skel_vert_vert_map next_bps = calc_next_bps(P, S, dists);
    
    std::vector<Polyhedron3*> segs;
    std::vector<double> values;
    S3VertexDesc v, u;
    Polyhedron3* seg;
    BOOST_FOREACH(auto x, segments)
    {
        std::tie(v, u, seg) = x;
        // Volume:
        //double val = CGAL::to_double(volume(seg));
        // Normalized volume:
        //double val = CGAL::to_double(volume(seg)/distance((*S)[v].point, (*S)[u].point));
        // Distance to soma:
        //double val = (dists[v] + dists[u]) / 2;
        // Distance to branch point:
        double val = (dists[v] + dists[u]) / 2 - dists[dists[v] < dists[u] ? next_bps[v] : next_bps[u]];
        values.push_back(val);
        segs.push_back(seg);
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
        BOOST_FOREACH(auto x, s) { delete x; }
    }*/

    // Close the OBJ file
    f.close();

    ///////////////////////////////////////////////////////////////////////////
    // Calculate the intersection of planes along the skeleton
    ///////////////////////////////////////////////////////////////////////////
    /*FacetTree ft(P->facets_begin(), P->facets_end(), *P);
    std::ofstream out("output.csv");
    skeleton_enum_branches_as_pts(S, [&ft, &out] (const std::vector<Point3>& pts) mutable
    {
        out << "Branch," << pts.front() << "," << pts.back() << std::endl;
        auto i = pts.begin(), end = pts.end();
        Point3 pt_end = *i++;
        const char* sep = "";
        while (i != end)
        {
            Point3 pt = pt_end; pt_end = *i++;
            Point3 mid = CGAL::midpoint(pt, pt_end);
            Intersection intrstn = Intersection(ft, Plane3(mid, normalized(pt_end - pt))).keep_around_point(mid);
            out << sep << intrstn.area();
            sep = ",";
        }
        out << std::endl;
    });
    out.close();*/


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
