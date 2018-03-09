#include "Skeleton.hpp"

#include "GeometryUtils.hpp"

#include <CGAL/Mean_curvature_flow_skeletonization.h>

#include <algorithm>


///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton using mean curvature flow.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(const Polyhedron3* P, double quality_factor)
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
    skeletonizer.set_quality_speed_tradeoff(quality_factor);

    skeletonizer.contract_until_convergence();
    skeletonizer.convert_to_skeleton(*S);
    return S;
}


///////////////////////////////////////////////////////////////////////////////////////
// Other skeleton utilties.
///////////////////////////////////////////////////////////////////////////////////////
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
