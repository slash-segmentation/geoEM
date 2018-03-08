#include "Skeleton.hpp"

#include "GeometryUtils.hpp"

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>

#include <algorithm>


///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton using mean curvature flow.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P)
{
    Skeleton3* S = new Skeleton3();
    CGAL::Mean_curvature_flow_skeletonization<Polyhedron3> skeletonizer(*P);

    // Default values:
    // max triangle angle = 110 degrees
    // min edge length = 0.002 * diagonal of bounding box
    // max iterations = 500
    // area variation factor = 0.0001
    // quality speed tradeoff = 0.1             - larger values result in high quality skeletons
    // is medially centered = true
    // medially centered speed tradeoff = 0.2   - larger value is more accurate to medial but not as smooth or fast
    //skeletonizer.set_max_iterations(...);
    //skeletonizer.set_quality_speed_tradeoff(...);

    skeletonizer.contract_until_convergence();
    skeletonizer.convert_to_skeleton(*S);
    return S;
}


///////////////////////////////////////////////////////////////////////////////////////
// Other skeleton utilties.
///////////////////////////////////////////////////////////////////////////////////////
struct Skeleton2Graph
{
    typedef SkeletonGraph3::BranchPoint_handle BPH;
    typedef std::unordered_map<Point3, BPH, boost::hash<Point3>> pt2bp_map;
    const Skeleton3& s;
    SkeletonGraph3* g;
    pt2bp_map pt2bp;
    SkeletonGraph3::Branch branch;
    Skeleton2Graph(const Skeleton3& skel)
        : s(skel), g(new SkeletonGraph3())
    {}
    void start_new_polyline()
    {
        branch.clear();
    }
    void add_node(Skeleton3::vertex_descriptor v)
    {
        branch.push_back(s[v].point);
    }
    void end_polyline()
    {
        if (branch.size() == 0) { return; }

        Point3 a = branch.front(), b = branch.back();
        pt2bp_map::iterator itr;

        // TODO: this creates branch points at the terminals!

        // Find or add the starting branch point
        itr = pt2bp.find(a);
        if (itr == pt2bp.end()) { itr = pt2bp.insert(itr, std::make_pair(a, g->add_branch_point(a))); }
        BPH bp_0 = itr->second;

        // Find or add the ending branch point
        itr = pt2bp.find(b);
        if (itr == pt2bp.end()) { itr = pt2bp.insert(itr, std::make_pair(b, g->add_branch_point(b))); }
        BPH bp_1 = itr->second;

        g->add_branch(bp_0, bp_1, branch);
    }
};

SkeletonGraph3* construct_skeleton_graph(const Skeleton3* S)
{
    Skeleton2Graph s2g(*S);
    CGAL::split_graph_into_polylines(*S, s2g);
    s2g.g->shrink_to_fit();
    return s2g.g;
}

/*
SkeletonGraph3* construct_skeleton_graph(const Skeleton3* S)
{
    typedef std::unordered_map<Skeleton3::vertex_descriptor, SkeletonGraph3::BranchPoint_handle> vert2bp;
    vert2bp branch_points;

    SkeletonGraph3* SG = new SkeletonGraph3();

    // Find and create all branch points (have degree >2)
    Skeleton3::vertex_iterator V, Vend;
    for (boost::tie(V, Vend) = boost::vertices(*S); V != Vend; ++V)
    {
        if (boost::out_degree(*V, *S) > 2)
        {
            branch_points.insert(std::make_pair(*V, SG->add_branch_point((*S)[*V].point)));
        }
    }

    // Search out from every branch point to make branches
    for (vert2bp::iterator bp = branch_points.begin(), end = branch_points.end(); bp != end; ++bp)
    {
        Skeleton3::out_edge_iterator E, Eend;
        for (boost::tie(E, Eend) = boost::out_edges(bp->first, *S); E != Eend; ++E)
        {
            // Travel along edge until we hit a vertex that doesn't have degree 2
            SkeletonGraph3::Branch b;
            Skeleton3::vertex_descriptor v = bp->first;
            Skeleton3::edge_descriptor e = *E;
            for(;;)
            {
                b.push_back((*S)[v].point);
                v = (v == boost::source(e, *S)) ? boost::target(e, *S) : boost::source(e, *S);
                if (boost::out_degree(v, *S) != 2) { break; }
                Skeleton3::out_edge_iterator EX, EXend;
                boost::tie(EX, EXend) = boost::out_edges(v, *S);
                e = *((e == *EX) ? ++EX : EX);
            }
            b.push_back((*S)[v].point);
            SG->add_branch(bp->second, boost::out_degree(v, *S) == 1 ? SkeletonGraph3::BranchPoint_handle() : branch_points[v], b);
        }
    }

    SG->shrink_to_fit();
    return SG;
}
*/


struct Output_polylines {
    const Skeleton3& skeleton;
    std::ofstream& out;
    int polyline_size;
    std::stringstream sstr;
    Output_polylines(const Skeleton3& skeleton, std::ofstream& out)
        : skeleton(skeleton), out(out)
    {}
    void start_new_polyline()
    {
        polyline_size = 0;
        sstr.str("");
        sstr.clear();
    }
    void add_node(Skeleton3::vertex_descriptor v)
    {
        ++polyline_size;
        sstr << " " << skeleton[v].point;
    }
    void end_polyline() { out << polyline_size << sstr.str() << std::endl; }
};
void write_skeleton(const Skeleton3* S, std::ofstream& output)
{
    Output_polylines outputter(*S, output);
    CGAL::split_graph_into_polylines(*S, outputter);
}

void skeleton_remove_collinear(SkeletonGraph3* SG)
{
    for (SkeletonGraph3::Branch_iterator b = SG->branches_begin(), end = SG->branches_end(); b != end; ++b)
    {
        remove_collinear_points(&*b, false);
    }
}

void skeleton_reduce(SkeletonGraph3* SG, double threshold)
{
    if (threshold < 0) {
        std::vector<double> areas;
        // Calculate mean and save all areas
        for (SkeletonGraph3::Branch_iterator b = SG->branches_begin(), end = SG->branches_end(); b != end; ++b)
        {
            SkeletonGraph3::Branch_handle pts = b;
            if (pts->size() < 3) { continue; }
            // <p,q,r> starts out with the first three items of the list
            SkeletonGraph3::Branch::iterator r = pts->begin(), p = r++, q = r++;
            for (; r != pts->end(); ++r)
            {
                areas.push_back(CGAL::squared_area(*p, *q, *r));
            }
        }
        // Get the value at 20% through the areas
        size_t i = (size_t)(areas.size() * 0.2);
        std::nth_element(areas.begin(), areas.begin() + i, areas.end());
        threshold = areas[i];
    }

    // Perform reduction
    for (SkeletonGraph3::Branch_iterator b = SG->branches_begin(), end = SG->branches_end(); b != end; ++b)
    {
        remove_collinear_points(&*b, false);
        remove_nearly_collinear_points(&*b, threshold, false);
    }
}
