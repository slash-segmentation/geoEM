#include "GeometryTypes.hpp"

#include "IO.hpp"
#include "IO_OBJ.hpp"
#include "Intersection.hpp"
#include "MedialAxisTransform.hpp"
#include "MedialAxisTransform_IO.hpp"
#include "Skeleton.hpp"
#include "Segmentation.hpp"
#include "Segments2Cylinders.hpp"
#include "Points2Spheres.hpp"

#include <boost/timer/timer.hpp>

#include <iostream>

#ifdef CREATE_GUI
#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt.h>
#endif

#include <CGAL/subdivision_method_3.h>

Intersection get_intersection(const FacetTree& ft, Point3 pt, Vector3 n)
{
    // Gets an intersection of an plane (defined by pt and n) using the given FacetTree.
    // Only keeps the parts of the intersections that surround the given point (or subtract
    // from the areas that surround the given point).
    Plane3 h(pt, n);
    Point2 pt2 = to_2d(pt, h);
    Intersection intersection_full(ft, h), intersection(h);

    // During the first pass add all of the positive polygons which contain the point to the list
    // This really should just be a single polygon, but maybe not...
    for (size_t i = 0; i < intersection_full.count(); ++i)
    {
        IntersectionPolygon2& p = intersection_full[i];
        if (p.orientation() == CGAL::COUNTERCLOCKWISE && p.bounded_side(pt2) != CGAL::ON_UNBOUNDED_SIDE)
        {
            intersection.add(p);
        }
    }

    // Now add all polygons (positive or negative) which are contained within any of the already
    // found polygons
    size_t n_outside = intersection.count();
    if (n_outside == 0) { return intersection; }
    for (size_t i = 0; i < intersection_full.count(); ++i)
    {
        IntersectionPolygon2& p = intersection_full[i];
        bool already_added = false;
        for (size_t j = 0; j < n_outside; ++j)
        {
            if (p == intersection[j]) { already_added = true; break; }
        }
        if (already_added) { continue; }
        for (size_t j = 0; j < n_outside; ++j)
        {
            if (p.is_inside(intersection[j]))
            {
                intersection.add(p);
                break;
            }
        }
    }
    return intersection;
}

typedef boost::graph_traits<const Polyhedron3> P3_traits;

std::vector<Point3> compute_rings(const Polyhedron3* P, const Skeleton3* S)
{
    // TODO: not quite to the rings yet, just the points that form the mesh
    // Next steps are:
    //   turn those points into a mesh using the connectivity in the original mesh
    //   get the boundary edges of the mesh
    //      M->normalize_border(); go through M->border_edges_begin() to M->edges_end(); (or halfedges)
    //   sort those edges into the rings

    //std::vector<Segment3> segs;
    std::vector<Point3> pts;
    BOOST_FOREACH(Skeleton3::vertex_descriptor v, vertices(*S))
    {
        if (degree(v, *S) > 2)
        {
            //const Point3& skel_pt = (*S)[v].point;
            BOOST_FOREACH(P3_traits::vertex_descriptor mesh_v, (*S)[v].vertices)
            {
                //segs.push_back(Segment3(mesh_v->point(), skel_pt));
                pts.push_back(mesh_v->point()); // get(CGAL::vertex_point, *P, mesh_v) ?
            }
        }
    }
    //return segs;
    return pts;
}

template <class HDS>
class Build_segment : public CGAL::Modifier_base<HDS>
{
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    typedef handle_map<Polyhedron3::Vertex_const_handle, size_t> Vertex_int_map;
    typedef handle_set<Polyhedron3::Facet_const_handle> Facet_set;
    typedef boost::graph_traits<const Polyhedron3>::vertex_descriptor vertex_desc;
    const Polyhedron3* P;
    const std::vector<vertex_desc> verts;
public:
    Build_segment(const Polyhedron3* P, const std::vector<vertex_desc> verts) : P(P), verts(verts) {}
    void operator()(HDS& hds)
    {
        // Start building the polyhedron
        Builder B(hds, true);
        B.begin_surface(this->verts.size(), 2*this->verts.size()); // number of faces is very approximate

        // Add all of the vertices
        Vertex_int_map lookup;
        BOOST_FOREACH(vertex_desc v, this->verts)
        {
            B.add_vertex(v->point());
            lookup.insert({{v, lookup.size()}});
        }

        // Add faces (all faces that have all vertices in the lookup table)
        Facet_set faces; // only add each face once
        std::vector<size_t> v_inds;
        BOOST_FOREACH(vertex_desc v, this->verts)
        {
            FOR_FACETS_AROUND_VERTEX(v, f)
            {
                if (faces.count(f)) { continue; } // face already processed
                faces.insert(f); // now mark it as processed
                FOR_VERTICES_AROUND_FACET(f, vf)
                {
                    auto itr = lookup.find(vf);
                    if (itr != lookup.end()) { v_inds.push_back(itr->second); }
                    else { v_inds.clear(); break; } // do not add face - not all vertices are present
                }
                if (v_inds.size() != 0)
                {
                    B.begin_facet();
                    BOOST_FOREACH(size_t idx, v_inds) { B.add_vertex_to_facet(idx); }
                    B.end_facet();
                    v_inds.clear();
                }
            }
        }
        B.end_surface();
    }
};

std::vector<Polyhedron3*> compute_bp_surface(const Polyhedron3* P, const Skeleton3* S)
{
    // This creates a mesh that contains all of the facets that contributed to a branch point.
    // This does end up dropping important vertices where there were insufficient vertices to make
    // a facet. It is possible that what is needed is some way to create "degenerate facets" in the
    // builder above that connect neighboring vertices that are otherwise being dropped. That may
    // still leave some vertices that are completely isolated however.
    // 
    // If this was more complete I would also want to perform the following steps:
    //     get the boundary edges of the mesh
    //        M->normalize_border(); go through M->border_edges_begin() to M->edges_end(); (or halfedges?)
    //     sort those edges into the rings

    std::vector<Polyhedron3*> segs;
    BOOST_FOREACH(Skeleton3::vertex_descriptor v, vertices(*S))
    {
        if (degree(v, *S) > 2)
        {
            Build_segment<Polyhedron3::HalfedgeDS> bs(P, (*S)[v].vertices);
            Polyhedron3* seg = new Polyhedron3();
            seg->delegate(bs);
            segs.push_back(seg);

            // TODO: rest of steps
            //seg->normalize_border();
        }
    }
    return segs;
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
    std::string filename = "example-data/big.off";
    //std::string filename = "example-data/big_simplified.off";
    std::string obj_name = "";

    
    ///////////////////////////////////////////////////////////////////////////
    // Settings
    ///////////////////////////////////////////////////////////////////////////
    bool use_mat = false; // use Medial Axis Transform instead of Mean Curvature Flow to calculate skeleton
    int loop_subdivisions = use_mat ? 0 : 1; // add extra vertices to the input mesh
    bool use_sdf = false; // always use SDF instead of MCF skeleton to calculate segments (only used if use_mat is false)
    const char* output_obj = "output.obj"; // the output OBJ file
    bool output_bp_surface = true; // output the branch point surface which is currently not complete but adds lots of rendering time
    const char* output_skel = "skel.cgal"; // the output CGAL file for the skeleton points
    bool reduce_skel_graph = use_mat; // if the skeleton graph should have vertices reduced before being analyzed


    ///////////////////////////////////////////////////////////////////////////
    // Read in the 3D polyhedral mesh
    ///////////////////////////////////////////////////////////////////////////
    std::cout << "Reading mesh..." << std::endl;
    Polyhedron3 *P;
    {
        boost::timer::auto_cpu_timer t;
        P = read_mesh(filename, assume_good); // TODO: delete P; somewhere...
    }
    std::cout << std::endl;

    
    ///////////////////////////////////////////////////////////////////////////
    // Refine the mesh
    ///////////////////////////////////////////////////////////////////////////
    if (loop_subdivisions > 0) {
        CGAL::Subdivision_method_3::Loop_subdivision(*P, loop_subdivisions);
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // Construct the skeleton
    ///////////////////////////////////////////////////////////////////////////
    // read with: S = read_cg("skeleton.cg");
    Skeleton3* S;
    if (use_mat)
    {
        // Construct the medial axis transform from the mesh
        std::cout << "Constructing medial axis transform..." << std::endl;
        MAT* mat;
        {
            boost::timer::auto_cpu_timer t;
            mat = construct_medial_axis_transform(P, filename, obj_name);
            //dump_mat(mat);
            //mat = load_mat("output.mat", P);
        }
        std::cout << std::endl;

        // Construct the skeleton from the medial axis transform
        std::cout << "Constructing skeleton from MAT..." << std::endl;
        {
            boost::timer::auto_cpu_timer t;
            S = construct_skeleton(mat, 0.7, 0.0); // TODO: how to determine parameters?
        }
        std::cout << std::endl;
        delete mat; // done with the medial axis transform
    }
    else
    {
        // Construct the skeleton using mean curvature flow
        std::cout << "Constructing skeleton using mean curvature flow..." << std::endl;
        {
            boost::timer::auto_cpu_timer t;
            S = construct_skeleton(P);
        }
        std::cout << std::endl;
    }
    // Save the skeleton in an easy-to-use format for analysis
    write_skeleton(S, output_skel, 10);

    
    ///////////////////////////////////////////////////////////////////////////
    // Compute Segmention
    ///////////////////////////////////////////////////////////////////////////
    // TODO: how to determine parameters?
    std::vector<Polyhedron3*> segments;
    if (use_mat || use_sdf)
    {
        // Segmention using SDF
        std::cout << "Segmenting..." << std::endl;
        {
            boost::timer::auto_cpu_timer t;
            segments = compute_segmentation(P,
                2.0/3.0*CGAL_PI, 25, // SDF value calculation parameters (2/3*pi, 25)
                15, 0.5);            // clustering parameters            (5, 0.26)
        }
        std::cout << std::endl;
    }
    else
    {
        // Segmention using MCF skeleton
        std::cout << "Segmenting using skeleton..." << std::endl;
        {
            boost::timer::auto_cpu_timer t;
            segments = compute_segmentation(P, S, 15 /*5*/, 0.5 /*0.26*/); 
        }
        std::cout << std::endl;
    }
    

    ///////////////////////////////////////////////////////////////////////////
    // Get the "graph" version of the skeleton - easier to use
    ///////////////////////////////////////////////////////////////////////////
    std::cout << "Constructing skeleton graph..." << std::endl;
    SkeletonGraph3* SG;
    {
        boost::timer::auto_cpu_timer t;
        SG = construct_skeleton_graph(S); // TODO: delete SG; somewhere...

        // Simplify the data by removing some of the coordinates
        //   Always remove collinear coordinates
        //   Optionally remove coordinates that form small triangles
        size_t nverts_before = SG->total_vertices();
        if (reduce_skel_graph) { skeleton_reduce(SG); }
        else { skeleton_remove_collinear(SG); }
        size_t nverts_after = SG->total_vertices();
        std::cerr << "Removed " << nverts_before - nverts_after << " vertices/edges out of " << nverts_before << std::endl;
    }
    std::cout << std::endl;


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
    
    // The skeleton graph as a series of cylinders
    std::vector<Polyhedron3*> skeleton_cyl;
    for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), Bend = SG->branches_end(); B != Bend; ++B)
    {
        SkeletonGraph3::Branch::const_iterator i = B->begin(), end = B->end();
        Point3 pt_end = *i++;
        std::vector<Segment3> branch;
        branch.reserve(B->size());
        while (i != end)
        {
            Point3 pt = pt_end; pt_end = *i++;
            branch.push_back(Segment3(pt, pt_end));
        }
        skeleton_cyl.push_back(segments2cylinders_merged(branch, 2.5, 8));
    }
    write_obj_file(f, skeleton_cyl, off, false, "", std::vector<std::string>{ "Black" });
    for (auto itr = skeleton_cyl.begin(), end = skeleton_cyl.end(); itr != end; ++itr) { delete *itr; }

    // The segments in the mesh in a myriad of colors
    std::vector<std::string> colors { // missing Black which is used for skeleton
        "Maroon", "Red", "Orange-Red", "Orange", "Gold", "Yellow", "Yellow-Green", "Lime",
        "Green", "Spring-Green", "Cyan", "Dark-Cyan", "Dark-Turquoise", "Deep-Sky-Blue", "Navy",
        "Medium-Blue", "Blue", "Royal-Blue", "Blue-Violet", "Indigo", "Purple", "Violet",
        "Magenta", "Hot-Pink", "Pink", "Brown", "Sienna", "Dark-Gray", "Gray", "Silver", "White",
    };
    write_obj_file(f, segments, off, false, "colors.mtl", colors);

    // The "branch point surface" as black spheres at the moment - should be polylines as cylinders
    if (output_bp_surface)
    {
        /*std::vector<Point3> rings = compute_rings(P, S);
        std::vector<Polyhedron3*> rings_spheres = points2spheres(rings, 5.0, 3);
        write_obj_file(f, rings_spheres, off, false, "", std::vector<std::string>{ "Black" });
        for (std::vector<Polyhedron3*>::iterator itr = rings_spheres.begin(), end = rings_spheres.end(); itr != end; ++itr) { delete *itr; }*/

        std::vector<Polyhedron3*> rings = compute_bp_surface(P, S);
        write_obj_file(f, rings, off, false, "", std::vector<std::string>{ "Black" });
        for (auto itr = rings.begin(), end = rings.end(); itr != end; ++itr) { delete *itr; }
    }

    // Close the OBJ file
    f.close();

    delete S; // done with the skeleton

    
    ///////////////////////////////////////////////////////////////////////////
    // Calculate the intersection of planes along the skeleton
    ///////////////////////////////////////////////////////////////////////////
    FacetTree ft(P->facets_begin(), P->facets_end(), *P);
    std::ofstream out("output.csv");
    for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), Bend = SG->branches_end(); B != Bend; ++B)
    {
        const Point3& a = B->source()->coord();
        const Point3& b = B->target()->coord();
        out << "Branch," << a << "," << b << std::endl;
        SkeletonGraph3::Branch::const_iterator i = B->begin(), end = B->end();
        Point3 pt_end = *i++;
        const char* sep = "";
        while (i != end)
        {
            Point3 pt = pt_end; pt_end = *i++;
            Intersection intersection = get_intersection(ft,
                CGAL::midpoint(pt, pt_end), normalized(pt_end - pt));
            out << sep << intersection.area();
            sep = ",";
        }
        out << std::endl;
    }
    out.close();


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
    viewer.set_skeleton(SG);
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
