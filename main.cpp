#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include "Skeleton.hpp"
#include "Polyhedron3Utils.hpp"
#include "Slice.hpp"

#include "IO.hpp"
#include "IO_OBJ.hpp"
#include "Segments2Cylinders.hpp"
#include "Points2Spheres.hpp"

#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#ifdef CREATE_GUI
#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt.h>
#endif

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
namespace PMP = CGAL::Polygon_mesh_processing;

// Polyhedron Property Maps
typedef boost::property_map<Polyhedron3, boost::face_index_t>::const_type P3_facet_index_map_t;
typedef boost::vector_property_map<size_t, P3_facet_index_map_t> P3_facet_int_map;

// Skeleton Property Maps
typedef boost::property_map<Skeleton3, boost::vertex_index_t>::const_type S3_vert_index_map_t;
typedef boost::vector_property_map<double, S3_vert_index_map_t> S3_vert_double_map;
typedef boost::vector_property_map<bool, S3_vert_index_map_t> S3_vert_bool_map;
typedef boost::vector_property_map<S3VertexDesc, S3_vert_index_map_t> S3_vert_vert_map;

Slice* find_largest_endpoint_slice(Slices slices, S3VertexDesc& sv)
{
    // Finds the endpoint that has largest normalized volume
    Kernel::FT max = 0;
    Slice* max_slc = nullptr;
    for (Slice* slc : slices)
    {
        if (slc->is_endpoint())
        {
            // Calculate the normalized volume
            Kernel::FT vol = PMP::volume(*slc->mesh()) / slc->length();
            // Check if it is a new maximum
            if (max_slc == nullptr || vol > max) { max_slc = slc; max = vol; sv = slc->endpoint(); }
            break;
        }
    }
    return max_slc;
}

S3_vert_double_map calc_distances(const Skeleton3* S, S3VertexDesc target)
{
    // Calculate the distance for each skeleton point to the given target
    S3_vert_double_map distances(get(boost::vertex_index, *S));
    S3_vert_bool_map discovered(get(boost::vertex_index, *S));
    std::vector<std::pair<S3VertexDesc, double>> stack;
    stack.push_back(std::make_pair(target, 0.0));
    discovered[target] = true;
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

S3_vert_vert_map calc_next_bps(const Skeleton3* S, const S3_vert_double_map& dist_to_soma)
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

template <class Graph>
bool is_cyclic(const Graph& g)
{
    // Assumes that the graph is a single connected component
    return num_vertices(g) - 1 != num_edges(g);
}

int main(int argc, char **argv)
{
    CGAL::set_pretty_mode(std::cout);

    ///////////////////////////////////////////////////////////////////////////
    // Set the file to read
    ///////////////////////////////////////////////////////////////////////////
    bool assume_good = false;
    //std::string filename = "example-data/other/elephant.off";
    //std::string filename = "example-data/trunc_cone.off";
    //std::string filename = "example-data/small.off";
    //std::string filename = "example-data/big.off";
    std::string filename = "example-data/big-nn-refined.off"; // no nucleus
    //std::string filename = "example-data/big_simplified.off";

    bool process_organelles = false;
    std::string filename_organelles = "...";


    ///////////////////////////////////////////////////////////////////////////
    // Settings
    ///////////////////////////////////////////////////////////////////////////
    int slice_sz = 3; // number of skeleton vertices to group together to form a slice
    const char* output_obj = "output.obj"; // the output OBJ file
    const char* output_obj_int = "intersections.obj"; // the output OBJ file for the intersections
    const char* output_skel = "skel.cgal"; // the output CGAL file for the skeleton points


    ///////////////////////////////////////////////////////////////////////////
    // Read in the 3D polyhedral mesh
    ///////////////////////////////////////////////////////////////////////////
    Polyhedron3* P;
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
        if (!assume_good && !is_single_component(P))
        {
            std::cerr << "ERROR: model is not a single connected component" << std::endl;
            return -1;
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
        if (is_cyclic(*S)) { std::cerr << "ERROR: Skeleton is cyclic!" << std::endl; return -1; }
    }
    std::cout << std::endl;
    // Save the skeleton in an easy-to-use format
    // NOTE: This loses the surface mesh vertices associated with every point
    write_skeleton(S, output_skel, 10);
    // read with: S = read_cg(output_skel);


    ///////////////////////////////////////////////////////////////////////////
    // Slice
    ///////////////////////////////////////////////////////////////////////////
    Slices slices;
    {
        std::cout << "Slicing..." << std::endl;
        boost::timer::auto_cpu_timer t;
        slices = slice(slice_sz, P, S); // TODO: for (Slice* slc : slices) { delete slc; } somewhere...
    }
    std::cout << std::endl;


    ///////////////////////////////////////////////////////////////////////////
    // Calculate skeleton vertex information
    ///////////////////////////////////////////////////////////////////////////
    S3VertexDesc soma_sv;
    Slice* soma_slc;
    S3_vert_double_map dists;
    S3_vert_vert_map next_bps;
    {
        std::cout << "Calculating skeleton vertex information..." << std::endl;
        boost::timer::auto_cpu_timer t;
        soma_slc = find_largest_endpoint_slice(slices, soma_sv);
        (void)soma_slc; // suppress unused variable warning
        dists = calc_distances(S, soma_sv);
        next_bps = calc_next_bps(S, dists);
    }
    std::cout << std::endl;


    ///////////////////////////////////////////////////////////////////////////
    // Calculate slice metrics
    ///////////////////////////////////////////////////////////////////////////
    std::vector<double> vols, vols_norm, surf_areas, surf_areas_norm, svrs, dist_avgs, dist_to_bps;
    {
        std::cout << "Calculating slice metrics..." << std::endl;
        boost::timer::auto_cpu_timer t;

        // Pre-allocate space
        vols.reserve(slices.size()); vols_norm.reserve(slices.size());
        surf_areas.reserve(slices.size()); surf_areas_norm.reserve(slices.size());
        svrs.reserve(slices.size());
        dist_avgs.reserve(slices.size()); dist_to_bps.reserve(slices.size());

        for (auto& slc : slices)
        {
            // Mesh metrics
            Kernel::FT vol = PMP::volume(*slc->mesh());
            Kernel::FT sa = PMP::area(*slc->uncapped_mesh());
            double len = slc->length();

            vols.push_back(CGAL::to_double(vol));
            vols_norm.push_back(CGAL::to_double(vol/len));
            surf_areas.push_back(CGAL::to_double(sa));
            surf_areas_norm.push_back(CGAL::to_double(sa/len));
            svrs.push_back(CGAL::to_double(sa/vol));

            // Distance metrics
            double dist = 0.0;
            for (auto sv : slc->skeleton_vertices()) { dist += dists[sv]; }
            dist /= slc->skeleton_vertices().size();
            dist_avgs.push_back(dist);

            if (slc->degree() > 2) { dist = 0; } // branch point
            else { dist -= dists[next_bps[*slc->skeleton_vertices().begin()]]; }
            dist_to_bps.push_back(dist);
        }
    }
    std::cout << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // Save the model as an OBJ file
    ///////////////////////////////////////////////////////////////////////////
    std::ofstream f(output_obj);
    f << std::setprecision(10);
    size_t off = 0;

    // The skeleton as a series of black cylinders
    std::vector<Polyhedron3*> skeleton_cyl;
    skeleton_enum_branches_as_pts(S, [&skeleton_cyl] (const std::vector<Point3>& pts)
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

    // The mesh as a series of colored segments based on their "value"
    std::vector<Polyhedron3*> segs;
    std::vector<double> values = dist_to_bps;
    for (auto& slc : slices)
    {
        //segs.push_back(slc->mesh());
        segs.push_back(slc->uncapped_mesh());
    }
    write_obj_cmap_segments(f, segs, values, off);

    // Close the OBJ file
    f.close();


    ///////////////////////////////////////////////////////////////////////////
    // Show statistics
    ///////////////////////////////////////////////////////////////////////////
    size_t total_facets = 0, total_halfedges = 0, total_vertices = 0;
    Kernel::FT total_volume = 0, total_sa = 0;
    for (auto& slc : slices)
    {
        total_facets += slc->mesh()->size_of_facets();
        total_halfedges += slc->mesh()->size_of_halfedges();
        total_vertices += slc->mesh()->size_of_vertices();
        total_volume += PMP::volume(*slc->mesh());
        total_sa += PMP::area(*slc->uncapped_mesh());
    }
    std::cout << "There are " << slices.size() << " slices" << std::endl;
    std::cout << "The original mesh contains:" << std::endl;
    std::cout << "    " << P->size_of_facets() << " facets" << std::endl;
    std::cout << "    " << P->size_of_halfedges() << " halfedges" << std::endl;
    std::cout << "    " << P->size_of_vertices() << " vertices" << std::endl;
    std::cout << "The slices contain:" << std::endl;
    std::cout << "    " << total_facets << " facets" << std::endl;
    std::cout << "    " << total_halfedges << " halfedges" << std::endl;
    std::cout << "    " << total_vertices << " vertices" << std::endl;
    std::cout << "The volume of the entire mesh is: " << PMP::volume(*P) << std::endl;
    std::cout << "The volume of the sum of the slices is: " << total_volume << std::endl;
    std::cout << "The surface area of the entire mesh is: " << PMP::area(*P) << std::endl;
    std::cout << "The surface area of the sum of the slices is: " << total_sa << std::endl;

    if (process_organelles)
    {
        ///////////////////////////////////////////////////////////////////////////
        // Read in the 3D polyhedral mesh for the organelles
        ///////////////////////////////////////////////////////////////////////////
        Polyhedron3 *Po;
        {
            std::cout << "Reading organelle mesh..." << std::endl;
            boost::timer::auto_cpu_timer t;
            try
            {
                Po = read_mesh(filename_organelles, assume_good); // TODO: delete Po; somewhere...
            }
            catch (std::invalid_argument& err) { std::cerr << err.what() << std::endl; return -1; }
            CGAL::set_halfedgeds_items_id(*Po);
            calculate_facet_planes(Po);
        }
        std::cout << std::endl;


        ///////////////////////////////////////////////////////////////////////////
        // Intersection of meshes
        ///////////////////////////////////////////////////////////////////////////
        // Both meshes must be non-self-intersecting and must bound a volume.
        // The output parameter can be one of the inputs for in-place operation.
        // TODO: Both input meshes will end up "refined" I believe.
        std::vector<Polyhedron3*> intrsctns;
        {
            std::cout << "Computing intersections..." << std::endl;
            boost::timer::auto_cpu_timer t;
            intrsctns.reserve(slices.size());
            for (auto& slc : slices)
            {
                Polyhedron3* intrsctn = new Polyhedron3(*Po);
                Polyhedron3 input(*slc->mesh());
                bool result = PMP::corefine_and_compute_intersection(input, *intrsctn, *intrsctn);
                if (!result)
                {
                    std::cerr << "ERROR: Failed to compute intersection" << std::endl;
                    delete intrsctn;
                    intrsctn = new Polyhedron3();
                }
                intrsctns.push_back(intrsctn); // TODO: for (Polyhedron3* intrsctn : intrsctns) { delete intrsctn; } somewhere...
            }
        }
        std::cout << std::endl;


        ///////////////////////////////////////////////////////////////////////////
        // Calculate the intersection metrics
        ///////////////////////////////////////////////////////////////////////////
        std::vector<double> volume_percs, organelle_svrs;
        {
            std::cout << "Calculating intersection metrics..." << std::endl;
            boost::timer::auto_cpu_timer t;
            volume_percs.reserve(intrsctns.size()); organelle_svrs.reserve(intrsctns.size());
            for (size_t i = 0; i < slices.size(); ++i)
            {
                Slice* slice = slices[i];
                Polyhedron3* intrsctn = intrsctns[i];

                Kernel::FT slc_vol = PMP::volume(*slice->mesh());
                Kernel::FT vol = PMP::volume(*intrsctn);
                Kernel::FT sa = PMP::area(*intrsctn); // TODO: needs to be uncapped!

                volume_percs.push_back(CGAL::to_double(vol/slc_vol));
                organelle_svrs.push_back(CGAL::to_double(sa/vol));
            }
        }
        std::cout << std::endl;

    }

#ifdef CREATE_GUI
    ///////////////////////////////////////////////////////////////////////////
    // Show GUI
    ///////////////////////////////////////////////////////////////////////////
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
