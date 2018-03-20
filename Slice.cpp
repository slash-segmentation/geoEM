#include "Slice.hpp"

#include "Skeleton.hpp"
#include "PolyBuilder.hpp"
#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>
#include <iostream>
#include <iterator>

#include <boost/foreach.hpp>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
namespace PMP = CGAL::Polygon_mesh_processing;

// Groups
typedef std::vector<std::vector<S3VertexDesc>> Groups;

template <class Container>
typename Container::value_type pop(Container& c)
{
    typename Container::value_type val = *c.begin();
    c.erase(c.begin());
    return val;
}

///////////////////////////////////////////////////////////////////////////////
// Slice class basic functions and methods
///////////////////////////////////////////////////////////////////////////////
void Slice::init()
{
    // Initializes a new slice. After this all fields should be ready except end_neighbors which
    // is filled in with nullptrs after this and must be corrected with set_neighbors.
    for (auto sv : svs)
    {
        if (boost::degree(sv, *S) > deg) { deg = boost::degree(sv, *S); }
        // If the vertex sv has an outgoing edge to a vertex not in svs it is an endpoint
        BOOST_FOREACH(auto e, out_edges(sv, *S))
        {
            S3VertexDesc op = opposite(*S, e, sv);
            if (!svs.count(op))
            {
                end_verts.push_back(sv);
                Point3 p = (*S)[op].point, q = (*S)[sv].point;
                end_planes.push_back(Plane3(CGAL::midpoint(p, q), q - p));
                end_neighbors.push_back(nullptr);
                break;
            }
        }
    }
    if (deg > 2)
    {
        // Setup aux planes at branch points
        for (auto& h : end_planes)
        {
            std::vector<Plane3> aux;
            for (auto& h2 : end_planes)
            {
                if (h == h2) { continue; }

                Vector3 n = h2.orthogonal_vector() - h.orthogonal_vector();
                auto result = CGAL::intersection(h2, h);
                Line3* l;
                if (result && (l = boost::get<Line3>(&*result)))
                {
                    aux.push_back(Plane3(l->point(0), n));
                }
                // else: planes are parallel and we will ignore them
            }
            aux_planes.push_back(aux);
        }
    }
}
std::vector<Plane3> Slice::all_planes() const
{
    // This include this->end_planes along with auxilary planes from branch points
    if (deg > 2) { return this->end_planes; } // branch points don't try to grab any other planes
    std::vector<Plane3> planes = this->end_planes;

    // Look for branch points by hopping along neighbors
    std::unordered_set<const Slice*> processed, stack;
    stack.insert(this);
    while (!stack.empty())
    {
        const Slice* slc = pop(stack);
        processed.insert(slc);
        for (Slice* neighbor : slc->end_neighbors)
        {
            if (processed.count(neighbor) || stack.count(neighbor) || neighbor->deg == 1) { continue; }
            if (neighbor->deg == 2) { stack.insert(neighbor); continue; }

            // Found a branch point
            // Look for which of its neighbor slc is
            size_t i;
            for (i = 0; i < neighbor->end_neighbors.size(); ++i)
            {
                if (slc == neighbor->end_neighbors[i]) { break; }
            }
            assert(i != neighbor->end_neighbors.size());

            // Add the BP's respective auxilary planes
            planes.insert(planes.end(), neighbor->aux_planes[i].begin(), neighbor->aux_planes[i].end());
        }
    }
    return planes;
}
Kernel::FT Slice::length() const
{
    // Calculates the skeletal length of the slice, including to the midpoints towards neighboring
    // groups.
    Kernel::FT len = 0;
    std::unordered_set<S3VertexDesc> processed;
    for (auto sv : svs)
    {
        const Point3& p = (*S)[sv].point;
        BOOST_FOREACH (auto e, out_edges(sv, *S))
        {
            auto sv2 = opposite(*S, e, sv);
            if (!processed.count(sv2))
            {
                Kernel::FT dist = distance(p, (*S)[sv2].point);
                if (!svs.count(sv2)) { dist /= 2; }  // part of other group, half distance (to midpoint)
                len += dist;
            }
        }
        processed.insert(sv);
    }
    return len;
}
void Slice::set_neighbors(std::vector<Slice*> slices)
{
    // Set all of the neighbors of all of the slices. This uses information from each slice to
    // determine which other slices are neighbors.

    // Create quick lookup of skeleton vertex to slice
    std::unordered_map<S3VertexDesc, Slice*> sv2slc;
    for (Slice* slc : slices) { for (auto sv : slc->svs) { sv2slc.insert({{sv, slc}}); } }
    // Now fill in all neighbors of every endpoint
    for (Slice* slc : slices)
    {
        for (size_t i = 0; i < slc->end_verts.size(); ++i)
        {
            S3VertexDesc sv = slc->end_verts[i];
            BOOST_FOREACH (auto e, out_edges(sv, *slc->S))
            {
                S3VertexDesc op = opposite(*slc->S, e, sv);
                if (!slc->svs.count(op)) { slc->end_neighbors[i] = sv2slc[op]; }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Slice class meshing methods
///////////////////////////////////////////////////////////////////////////////
inline static P3CVertex find_vertex_with_pt(const Polyhedron3* P, const Point3& p)
{
    for (auto v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v)
    {
        if (p == v->point())
        {
            return v;
        }
    }
    throw std::invalid_argument("point not found in polyhedron");
}
inline static bool on_all_pos_sides(const std::vector<Plane3> planes, const Point3& p)
{
    for (auto& h : planes) { if (!h.has_on_positive_side(p)) { return false; } }
    return true;
}
class CutAtPlane : public CGAL::Modifier_base<Polyhedron3::HalfedgeDS>
{
    typedef CGAL::Polyhedron_incremental_builder_3<Polyhedron3::HalfedgeDS> Builder;
    typedef std::unordered_map<Point3, size_t, boost::hash<Point3>> PointMap;
    typedef handle_map<P3CHalfedge, size_t> EdgeMap;

    const Polyhedron3* P;
    const Plane3 h;
    const P3CVertex v;
    Builder* B;
    PointMap lookup;
    EdgeMap edge_lookup;
    PolyBuilder<size_t> hole_edges;

    handle_set<P3CFacet> finished, stack;

    //mutable handle_map<P3CVertex, CGAL::Oriented_side> side_cache;
    inline CGAL::Oriented_side side(P3CVertex v) const
    {
        return this->h.oriented_side(v->point());
        // TODO: Caching the calculation of the side values seems to make the times more
        // consistently in the 8 sec range, without it the time range from 6.5 to 10.5 seconds...
        // I really don't know which way to go for this.
        /*auto itr = side_cache.find(v);
        if (itr != side_cache.end()) { return itr->second; }
        CGAL::Oriented_side side = this->h.oriented_side(v->point());
        side_cache.insert({{v, side}});
        return side;*/
    }
    inline bool on_pos_side(const P3CVertex& v) const { return this->side(v) == CGAL::Oriented_side::ON_POSITIVE_SIDE; }
    inline bool on_neg_side(const P3CVertex& v) const { return this->side(v) == CGAL::Oriented_side::ON_NEGATIVE_SIDE; }
    inline bool on_pos_side(const P3CHalfedge& he) const { return this->side(he->vertex()) == CGAL::Oriented_side::ON_POSITIVE_SIDE; }
    inline bool on_neg_side(const P3CHalfedge& he) const { return this->side(he->vertex()) == CGAL::Oriented_side::ON_NEGATIVE_SIDE; }

    inline size_t n_verts_on_pos_side(const P3CFacet& f) const
    {
        size_t n = 0;
        FOR_EDGES_AROUND_FACET(f, he) { n += this->on_pos_side(he); }
        return n;
    }
    inline size_t n_verts_on_neg_side(const P3CFacet& f) const
    {
        size_t n = 0;
        FOR_EDGES_AROUND_FACET(f, he) { n += this->on_neg_side(he); }
        return n;
    }
    inline P3CHalfedge find_neg_vert(const P3CFacet& f)
    {
        FOR_EDGES_AROUND_FACET(f, he) { if (this->on_neg_side(he)) { return he; } }
        throw std::invalid_argument("negative-side vertex not found");
    }
    inline P3CHalfedge find_pos_vert(const P3CFacet& f)
    {
        FOR_EDGES_AROUND_FACET(f, he) { if (this->on_pos_side(he)) { return he; } }
        throw std::invalid_argument("positive-side vertex not found");
    }
    inline size_t get_pt_id(const Point3& p)
    {
        auto itr = this->lookup.find(p);
        if (itr == this->lookup.end())
        {
            B->add_vertex(p);
            itr = this->lookup.insert(std::make_pair(p, this->lookup.size())).first;
        }
        return itr->second;
    }
    inline size_t get_pt_id(const P3CVertex& v) { return this->get_pt_id(v->point()); }
    inline size_t get_pt_id(const P3CHalfedge& he) { return this->get_pt_id(he->vertex()->point()); }
    inline size_t get_int_pt_id(const P3CHalfedge& he)
    {
        auto itr = this->edge_lookup.find(he);
        if (itr != this->edge_lookup.end()) { return itr->second; }
        for (;;)
        {
            auto intrsctn = CGAL::intersection(h, halfedge_to_segment3(he));
            if (intrsctn)
            {
                const Point3* p = boost::get<Point3>(&*intrsctn);
                if (!p) { throw std::domain_error("intersection was a segment"); }
                size_t id = this->get_pt_id(*p);
                this->edge_lookup.insert({{he, id}});
                this->edge_lookup.insert({{he->opposite(), id}});
                return id;
            }
        }
        throw std::domain_error("no intersection");
    }
    inline void add_op_facet(const P3CHalfedge& he)
    {
        auto e = he->opposite();
        if (!e->is_border() && !this->finished.count(e->facet()))
        {
            this->stack.insert(e->facet());
        }
    }
    inline void add_hole_edge(size_t v, size_t u)
    {
        this->hole_edges.add_seg(v, u);
    }
    inline void create_facet(size_t a, size_t b, size_t c, const Plane3& h)
    {
        B->begin_facet()->plane() = h;
        B->add_vertex_to_facet(a);
        B->add_vertex_to_facet(b);
        B->add_vertex_to_facet(c);
        B->end_facet();
    }
public:
    // Polyhedron must be pure triangle
    // Keep positive side of plane
    // The vertex must be on positive side of plane
    CutAtPlane(const Polyhedron3* P, Plane3 h, P3CVertex v) : P(P), h(h), v(v) {}
    void operator()(Polyhedron3::HalfedgeDS& hds)
    {
        assert(P->is_pure_triangle());
        assert(this->on_pos_side(v));

        // Start building the polyhedron
        Builder B(hds, true);
        this->B = &B;
        size_t facet_count = 0, vertex_count = 0;
        for (auto f = P->facets_begin(), end = P->facets_end(); f != end; ++f)
        {
            size_t n = this->n_verts_on_pos_side(f);
            if (n == 1) { facet_count += 1; vertex_count += 1; }
            else if (n == 2) { facet_count += 2; vertex_count += 2; }
            else if (n == 3) { facet_count += 1; }
        }
        this->edge_lookup.reserve(vertex_count);
        for (auto v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v) { vertex_count += this->on_pos_side(v); }
        B.begin_surface(vertex_count, facet_count, this->P->size_of_halfedges());
        this->lookup.reserve(vertex_count);

        // Start processing based on the given seed vertex
        FOR_EDGES_AROUND_VERTEX(this->v, he) { if (!he->is_border()) { this->stack.insert(he->facet()); } }
        while (!this->stack.empty())
        {
            P3CFacet f = pop(this->stack);
            this->finished.insert(f);
            size_t n_neg = this->n_verts_on_neg_side(f);
            assert(n_neg < 3);
            if (n_neg == 0)
            {
                // Add the entire facet
                const auto &a = f->halfedge(), &b = a->next(), &c = b->next();
                assert(a == c->next());
                create_facet(this->get_pt_id(a), this->get_pt_id(b), this->get_pt_id(c), f->plane());
                this->add_op_facet(a); this->add_op_facet(b); this->add_op_facet(c);
            }
            else if (n_neg == 1)
            {
                // The edge pointing to negative side and two vertices on the non-negative side
                const P3CHalfedge& he = this->find_neg_vert(f);
                size_t a = this->get_pt_id(he->next());
                size_t b = this->get_pt_id(he->prev());
                // Number of positive vertices
                size_t n_pos = this->n_verts_on_pos_side(f);
                if (n_pos == 0) { continue; } // 1 negative, 2 on the plane -> no facet
                else if (n_pos == 1)
                {
                    // 1 negative, 1 on the plane, 1 positive -> single facet
                    const P3CHalfedge& he_cross = this->on_pos_side(he->next()) ? he->next() : he;
                    size_t c = this->get_int_pt_id(he_cross);
                    create_facet(a, b, c, f->plane());
                    // Record the edge as part of the hole
                    this->add_hole_edge(this->get_pt_id(he_cross->next()), c);
                    // Add neighboring facets
                    this->add_op_facet(he_cross);
                }
                else // if (n_pos == 2)
                {
                    // 1 negative, 2 positive -> two facets:
                    //   1 vertex already on the positive side and two on the plane
                    //   2 vertices already on the positive side and one on the plane
                    // The two vertices on the plane
                    size_t c = this->get_int_pt_id(he);
                    size_t d = this->get_int_pt_id(he->next());
                    // Add the facets
                    create_facet(b, c, d, f->plane());
                    create_facet(b, d, a, f->plane());
                    // Record the edge as part of the hole
                    this->add_hole_edge(c, d);
                    // Add neighboring facets
                    this->add_op_facet(he);
                    this->add_op_facet(he->next());
                }
                this->add_op_facet(he->prev());

            }
            else //if (n_neg == 2)
            {
                if (this->n_verts_on_pos_side(f) == 0) { continue; } // 2 negative 1 on the plane -> no facet

                // Add a facet with the vertex that is already on the positive side plus two
                // vertices that lie on the plane itself from the edges.
                const P3CHalfedge& he = this->find_pos_vert(f);
                size_t a = this->get_pt_id(he), b = this->get_int_pt_id(he->next()), c = this->get_int_pt_id(he);
                create_facet(a, b, c, f->plane());
                // Record the edge as part of the hole
                this->add_hole_edge(b, c);
                // Add neighboring facets
                this->add_op_facet(he);
                this->add_op_facet(he->next());
            }
        }

        // Finish up with filling the holes
        for (auto& hole : this->hole_edges.finish_up())
        {
            std::reverse(hole.begin(), hole.end()); // TODO: always necessary? or use B.test_facet?
            std::vector<Point3> points;
            points.reserve(hole.size());
            for (size_t i : hole) { points.push_back(B.vertex(i)->point()); }
            PMP::triangulate_hole_polyline(points, triangle_output_iterator(this, points));
        }

        B.end_surface();
    }
private:
    struct triangle_output
    {
        int a, b, c;
        triangle_output() : a(0), b(0), c(0) { }
        triangle_output(int a, int b, int c) : a(a), b(b), c(c) { }
    };
    struct triangle_output_iterator : std::iterator<std::output_iterator_tag, triangle_output>
    {
        CutAtPlane* cut;
        std::vector<Point3>& pts;
        triangle_output_iterator(CutAtPlane* cut, std::vector<Point3>& points) : cut(cut), pts(points) { }
        void operator=(triangle_output const& t)
        {
            cut->create_facet(cut->get_pt_id(pts[t.a]), cut->get_pt_id(pts[t.b]), cut->get_pt_id(pts[t.c]), cut->h);
        }
        triangle_output_iterator& operator++() { return *this; }
        triangle_output_iterator operator++(int) { return *this; }
        triangle_output_iterator& operator*() { return *this; }
    };
};
void Slice::build_mesh(const Polyhedron3* P)
{
    assert(_mesh->empty());

    // Get all planes that will restrict the mesh
    std::vector<Plane3> planes = all_planes();

    // Find all of the polyhedron vertices that should be in the final result
    P3CVertex seed = P->vertices_end();
    for (S3VertexDesc sv : svs)
    {
        for (P3CVertex v : (*S)[sv].vertices)
        {
            if (on_all_pos_sides(planes, v->point())) { seed = v; break; }
        }
        if (seed != P->vertices_end()) { break; }
    }
    if (seed == P->vertices_end())
    {
        // TODO: this is occuring somewhat frequently (even going on to the error) - why?
        std::cerr << "Warning: no seed found using shortcut, using fallback..." << std::endl;
        for (auto v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v)
        {
            if (on_all_pos_sides(planes, v->point())) { seed = v; break; }
        }
        if (seed == P->vertices_end()) { std::cerr << "ERROR: no seed found to start building mesh from" << std::endl; }
    }

    // Cut up the mesh into a new mesh
    Polyhedron3* mesh = new Polyhedron3();
    bool first = true;
    Point3 p = seed->point();
    for (auto& h : planes)
    {
        if (first)
        {
            // First cut - uses original mesh
            CutAtPlane cut(P, h, seed);
            mesh->delegate(cut);
            first = false;
        }
        else
        {
            // All other cuts - use previous result
            P3CVertex v = find_vertex_with_pt(mesh, p);
            Polyhedron3* temp = new Polyhedron3();
            CutAtPlane cut(mesh, h, v);
            temp->delegate(cut);
            delete mesh;
            mesh = temp;
        }
    }
    CGAL::set_halfedgeds_items_id(*mesh);

    // Check mesh
    if (!mesh->is_valid())         { std::cerr << "Warning: slice is not valid" << std::endl; }
    if (!mesh->is_closed())        { std::cerr << "Warning: slice is not closed" << std::endl; }
    if (!is_not_degenerate(mesh))  { std::cerr << "Warning: slice is degenerate" << std::endl; }
    if (PMP::does_self_intersect(*mesh)) { std::cerr << "Warning: slice is self-intersecting" << std::endl; }
    else if (!PMP::does_bound_a_volume(*mesh)) { std::cerr << "Warning: slice does not bound a volume" << std::endl; } // crashes in many cases...
    if (!mesh->is_pure_triangle()) { std::cerr << "Warning: slice is not pure triangle" << std::endl; }

    // Set the class's mesh
    delete _mesh;
    _mesh = mesh;

    // Create uncapped version of mesh
    std::unordered_set<P3Facet> facets;
    for (auto f = mesh->facets_begin(), end = mesh->facets_end(); f != end; ++f)
    {
        // This would add capping-only facets (which is an interesting view)
        //for (auto& h : planes) { if (f->plane() == h) { facets.insert(f); break; } }
        // This adds original facets only
        bool any_match = false;
        for (auto& h : planes) { if (f->plane() == h) { any_match = true; break; } }
        if (!any_match) { facets.insert(f); }
    }
    CGAL::Face_filtered_graph<Polyhedron3> uncapped_graph(*mesh, facets);
    if (uncapped) { delete uncapped; uncapped = nullptr; }
    CGAL::copy_face_graph(uncapped_graph, *(uncapped = new Polyhedron3()));
}

///////////////////////////////////////////////////////////////////////////////
// Group creation
///////////////////////////////////////////////////////////////////////////////
static void create_groups(const int group_sz, const int bp_group_sz, const Skeleton3* S, Groups& groups)
{
    // Creates groups of skeleton vertices. The branch points in the skeleton will always be in put
    // into groups so that the distance from one endpoint to another is bp_group_sz (or
    // bp_group_sz+1 is bp_group_sz is even). This means that the branch part groups will contain
    // 1 + (bp_group_sz/2*deg) skeleton vertices where deg is the degree of the branch point.
    //
    // Other groups will contain at most group_sz consecutive skeleton vertices. In the case
    // branches cannot be evenly divided by group_sz the groups will be made smaller while trying to
    // keep all of the groups roughly the same size (at most different by 1).
    const int bp_size = bp_group_sz / 2;

    groups.reserve(num_vertices(*S) / group_sz);
    // First add all of the branch points with their neighbors
    BOOST_FOREACH(auto sv, vertices(*S))
    {
        if (degree(sv, *S) > 2)
        {
            std::unordered_set<S3VertexDesc> svs;
            svs.reserve(degree(sv, *S)*bp_size + 1);
            std::unordered_set<S3VertexDesc> next;
            next.insert(sv);
            for (int distance = 1; distance <= bp_size; ++distance)
            {
                std::unordered_set<S3VertexDesc> now(next);
                svs.insert(next.begin(), next.end());
                next.clear();
                for (auto sv : now)
                {
                    BOOST_FOREACH(auto e, out_edges(sv, *S))
                    {
                        auto sv_n = opposite(*S, e, sv);
                        // TODO: how to deal with BPs close together?
                        if (!svs.count(sv_n) && degree(sv_n, *S) <= 2)
                        {
                            next.insert(sv_n);
                        }
                    }
                }
            }
            svs.insert(next.begin(), next.end());
            groups.push_back(std::vector<S3VertexDesc>(svs.begin(), svs.end()));
        }
    }
    // Then add all others
    skeleton_enum_branches(S, [&, S, group_sz] (const std::vector<S3VertexDesc>& verts) mutable
    {
        int start = (degree(verts.front(), *S) != 1)*(bp_size+1);
        int end = (degree(verts.back(), *S) != 1)*(bp_size+1);
        ssize_t total = (ssize_t)verts.size()-end-start;
        if (total <= 0) { return; }
        std::vector<int> sizes(total/group_sz, group_sz);
        // Distribute the remainder
        int rem = total%group_sz;
        if (rem != 0)
        {
            // TODO: redistribute smaller groups to be where point distance is greater?
            sizes.push_back(group_sz);
            for (size_t to_remove = group_sz - rem, i = sizes.size() - 1; to_remove > 0; to_remove -= 1)
            {
                sizes[i] -= 1;
                if (i == 0) { i = sizes.size(); }
                i -= 1;
            }
        }
        // Add the groups of the correct sizes
        auto itr = verts.begin() + start;
        for (size_t sz : sizes)
        {
            groups.push_back(std::vector<S3VertexDesc>(itr, itr + sz));
            itr += sz;
        }
        assert(itr + end == verts.end());
    });
}

///////////////////////////////////////////////////////////////////////////////
// Function that takes a mesh and skeleton and creates the slices
///////////////////////////////////////////////////////////////////////////////
Slices slice(const int group_sz, const int bp_group_sz, const Polyhedron3* P, const Skeleton3* S)
{
    // Create the groups
    Groups groups;
    create_groups(group_sz, bp_group_sz, S, groups);

    // Create each slice
    Slices slices;
    slices.reserve(groups.size());
    for (auto& group : groups) { slices.push_back(new Slice(S, group.begin(), group.end())); }

    // Once all slices are created setup neighbors
    Slice::set_neighbors(slices);

    // Finally the meshes for each of the slices can be built
    size_t i = 0;
    for (Slice* slc : slices)
    {
        std::cout << "#" << i++ << " ";
        slc->build_mesh(P);
    }
    std::cout << std::endl;

    return slices;
}
