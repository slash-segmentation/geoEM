#include "Slice.hpp"

#include "Skeleton.hpp"
#include "PolyBuilder.hpp"
#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>
#include <iostream>
#include <iterator>

#include <boost/foreach.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
namespace PMP = CGAL::Polygon_mesh_processing;

// Incremental builder
typedef CGAL::Polyhedron_incremental_builder_3<Polyhedron3::HalfedgeDS> Builder;

// Polyhedron Property Map for connected components
typedef boost::property_map<Polyhedron3, boost::face_index_t>::const_type P3_facet_index_map_t;
typedef boost::vector_property_map<size_t, P3_facet_index_map_t> P3_facet_int_map;

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
    // This includes this->end_planes along with auxilary planes from branch points
    if (deg > 2) { return this->end_planes; } // branch points don't try to grab any other planes
    std::vector<Plane3> planes = this->end_planes;

    // Look for branch points by hopping along neighbors
    // TODO: there should likely be a limit on how many jumps to look so that hair-pin branches
    // don't get destroyed. However the limit needs to be more than 1, just not sure how much more.
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
double Slice::length() const
{
    // Calculates the skeletal length of the slice, including to the midpoints towards neighboring
    // groups.
    double len = 0;
    std::unordered_set<S3VertexDesc> processed;
    for (auto sv : svs)
    {
        const Point3& p = (*S)[sv].point;
        BOOST_FOREACH (auto e, out_edges(sv, *S))
        {
            auto sv2 = opposite(*S, e, sv);
            if (!processed.count(sv2))
            {
                double dist = distance(p, (*S)[sv2].point);
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
        if (p == v->point()) { return v; }
    }
    throw std::invalid_argument("point not found in polyhedron");
}
inline static bool on_all_pos_sides(const std::vector<Plane3>& planes, const Point3& p)
{
    for (auto& h : planes) { if (!h.has_on_positive_side(p)) { return false; } }
    return true;
}

// These are for triangulating holes during building
struct triangle_output
{
    const int a, b, c;
    triangle_output(int a, int b, int c) : a(a), b(b), c(c) { }
};
struct triangle_output_iterator : std::iterator<std::output_iterator_tag, triangle_output>
{
    Builder& B;
    const Plane3& h;
    const std::vector<size_t>& hole;
    triangle_output_iterator(Builder& B, const Plane3& h, const std::vector<size_t>& hole)
        : B(B), h(h), hole(hole) { }
    void operator=(triangle_output const& t)
    {
        B.begin_facet()->plane() = h;
        B.add_vertex_to_facet(hole[t.a]);
        B.add_vertex_to_facet(hole[t.b]);
        B.add_vertex_to_facet(hole[t.c]);
        B.end_facet();
    }
    triangle_output_iterator& operator++() { return *this; }
    triangle_output_iterator operator++(int) { return *this; }
    triangle_output_iterator& operator*() { return *this; }
};
enum struct Reverse { no, yes, check };
void triangulate_hole(Builder& B, const Plane3& h, std::vector<size_t>& hole, Reverse rev=Reverse::check)
{
    if (rev == Reverse::yes || (rev == Reverse::check && !B.test_facet(hole.begin(), hole.end())))
    {
        std::reverse(hole.begin(), hole.end());
    }
    auto get_pt = [&B] (size_t i) -> const Point3& { return B.vertex(i)->point(); };
    PMP::triangulate_hole_polyline(boost::adaptors::transform(hole, get_pt),
                                   triangle_output_iterator(B, h, hole));
}

// For triangulating holes outside of building (don't need the output)
struct null_output_iterator : std::iterator<std::output_iterator_tag, null_output_iterator>
{
    template<typename T> void operator=(T const&) { }
    null_output_iterator & operator++() { return *this; }
    null_output_iterator operator++(int) { return *this; }
    null_output_iterator & operator*() { return *this; }
};


class CutAtPlane : public CGAL::Modifier_base<Polyhedron3::HalfedgeDS>
{
    const Polyhedron3* P;
    const Plane3 h;
    const P3CVertex v;
    Builder* B;
    size_t lookup_count = 0;
    std::vector<size_t> lookup, edge_lookup; // (size_t)-1 indicates a point/edge that hasn't been added to the mesh yet
    std::vector<char> reached; // not <bool> since that has odd behaviors
    std::vector<P3CFacet> stack;
    PolyBuilder<size_t> hole_edges;

    std::vector<CGAL::Oriented_side> side_cache;
    inline CGAL::Oriented_side side(const P3CVertex& v) const { return this->side_cache[v->id()]; }
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
    inline size_t get_pt_id(const P3CHalfedge& he)
    {
        const P3CVertex& v = he->vertex();
        size_t id = v->id(), i = this->lookup[id];
        if (i == (size_t)-1)
        {
            B->add_vertex(v->point());
            this->lookup[id] = i = this->lookup_count++;
        }
        return i;
    }
    inline size_t get_int_pt_id(const P3CHalfedge& he)
    {
        size_t id = he->id(), i = this->edge_lookup[id];
        if (i == (size_t)-1)
        {
            auto intrsctn = CGAL::intersection(h, halfedge_to_segment3(he));
            if (!intrsctn) { throw std::domain_error("no intersection"); }
            const Point3* p = boost::get<Point3>(&*intrsctn);
            if (!p) { throw std::domain_error("intersection was a segment"); } // TODO: handle this situation?
            B->add_vertex(*p);
            this->edge_lookup[id] = i = this->lookup_count++;
            this->edge_lookup[he->opposite()->id()] =  i;
        }
        return i;
    }
    inline void add_op_facet(const P3CHalfedge& he)
    {
        const P3CFacet& f = he->opposite()->facet();
        size_t id = f->id();
        if (!this->reached[id])
        {
            this->reached[id] = true;
            this->stack.push_back(f);
        }
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
    // Polyhedron must be pure triangle and closed
    // Keep positive side of plane
    // The vertex must be on positive side of plane
    // The resulting mesh will be pure triangle and closed
    CutAtPlane(const Polyhedron3* P, Plane3 h, P3CVertex v) : P(P), h(h), v(v)
    {
        assert(P->is_closed());
        assert(P->is_pure_triangle());
        assert(h.has_on_positive_side(v->point()));
        this->side_cache.resize(this->P->size_of_vertices());
        for (auto v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v)
        {
            this->side_cache[v->id()] = this->h.oriented_side(v->point());
        }
    }
    void operator()(Polyhedron3::HalfedgeDS& hds)
    {
        // Start building the polyhedron
        Builder B(hds, true);
        this->B = &B;
        B.begin_surface(P->size_of_vertices(), P->size_of_facets()*2); // *2 for the caps

        // Start processing based on the given sees, pre-adding it as vertex
        this->lookup.resize(P->size_of_vertices(), (size_t)-1);
        this->edge_lookup.resize(P->size_of_halfedges(), (size_t)-1);
        this->reached.resize(P->size_of_facets(), false);
        B.add_vertex(this->v->point());
        this->lookup[this->v->id()] = this->lookup_count++;
        FOR_FACETS_AROUND_VERTEX(this->v, f)
        {
            this->reached[f->id()] = true;
            this->stack.push_back(f);
        }

        // Process the stack until empty
        while (!this->stack.empty())
        {
            const P3CFacet f = this->stack.back(); this->stack.pop_back();
            size_t n_neg = this->n_verts_on_neg_side(f);
            if (n_neg == 0)
            {
                // Add the entire facet
                const P3CHalfedge &a = f->halfedge(), &b = a->next(), &c = b->next();
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
                    const P3CHalfedge he_cross = this->on_pos_side(he->next()) ? he->next() : he;
                    size_t c = this->get_int_pt_id(he_cross);
                    create_facet(a, b, c, f->plane());
                    // Record the edge as part of the hole
                    this->hole_edges.add_seg(this->get_pt_id(he_cross->next()), c);
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
                    this->hole_edges.add_seg(c, d);
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
                this->hole_edges.add_seg(b, c);
                // Add neighboring facets
                this->add_op_facet(he);
                this->add_op_facet(he->next());
            }
        }

        // Finish up with filling the holes
        for (auto& hole : this->hole_edges.finish_up()) { triangulate_hole(B, h, hole, Reverse::yes); }

        // Done
        B.end_surface();
    }
};
class QuickCutAtPlanes : public CGAL::Modifier_base<Polyhedron3::HalfedgeDS>
{
    const Polyhedron3* P;
    const std::vector<Plane3>& planes;
    const std::vector<P3CVertex>& seeds;
    Builder* B;
    size_t lookup_count = 0;
    std::vector<size_t> lookup; // (size_t)-1 indicates a point that hasn't been added to the mesh yet
    std::vector<char> reached; // not <bool> since that has odd behaviors
    std::vector<P3CFacet> stack;

    inline bool on_pos_side(const P3CHalfedge& he) const
    {
        const Point3& p = he->vertex()->point();
        for (const Plane3& h : this->planes) { if (!h.has_on_positive_side(p)) { return false; } }
        return true;
    }
    inline bool any_vert_on_pos_side(const P3CFacet& f) const
    {
        FOR_EDGES_AROUND_FACET(f, he) { if (this->on_pos_side(he)) { return true; } }
        return false;
    }
    inline size_t get_pt_id(const P3CHalfedge& he)
    {
        const P3CVertex& v = he->vertex();
        size_t id = v->id();
        size_t i = this->lookup[id];
        if (i == (size_t)-1)
        {
            B->add_vertex(v->point());
            this->lookup[id] = i = this->lookup_count++;
        }
        return i;
    }
    inline void add_op_facet(const P3CHalfedge& he)
    {
        const P3CFacet& f = he->opposite()->facet();
        size_t id = f->id();
        if (!this->reached[id])
        {
            this->reached[id] = true;
            this->stack.push_back(f);
        }
    }
    // The quick-and-dirty version of CutAtPlane that can deal with multiple planes but lets facets
    // cross the planes a little bit and does not close the holes.
public:
    // Polyhedron must be pure triangle and closed
    // Keep positive side of all given planes (along with some parts of facets on the negative sides)
    // The seed vertices must be on positive sides of all planes
    // The resulting mesh will be pure triangle and NOT closed (all attempts at fixing it failed so just doing triangulation outside of the cutter)
    QuickCutAtPlanes(const Polyhedron3* P, const std::vector<Plane3>& planes, const std::vector<P3CVertex>& seeds)
        : P(P), planes(planes), seeds(seeds)
    {
        assert(P->is_closed());
        assert(P->is_pure_triangle());
        #ifdef _DEBUG
        for (const P3CVertex& v : seeds) { assert(this->on_pos_side(v->halfedge())); }
        #endif
    }
    void operator()(Polyhedron3::HalfedgeDS& hds)
    {
        // Start building the polyhedron
        Builder B(hds, true);
        this->B = &B;
        B.begin_surface(P->size_of_vertices(), P->size_of_facets(), P->size_of_halfedges()); // very conservative estimate

        // Start processing based on the given seeds, pre-adding them as vertices
        this->lookup.resize(P->size_of_vertices(), (size_t)-1);
        this->reached.resize(P->size_of_facets(), false);
        this->stack.reserve(seeds.size()*2);
        for (P3CVertex v : this->seeds)
        {
            B.add_vertex(v->point());
            this->lookup[v->id()] = this->lookup_count++;
            FOR_FACETS_AROUND_VERTEX(v, f)
            {
                size_t id = f->id();
                if (!this->reached[id])
                {
                    this->reached[id] = true;
                    this->stack.push_back(f);
                }
            }
        }

        // Process the stack until empty
        while (!this->stack.empty())
        {
            const P3CFacet f = this->stack.back(); this->stack.pop_back();
            const P3CHalfedge &a = f->halfedge(), &b = a->next(), &c = b->next();
            const size_t ai = this->get_pt_id(a), bi = this->get_pt_id(b), ci = this->get_pt_id(c);
            B.begin_facet()->plane() = f->plane();
            B.add_vertex_to_facet(ai); B.add_vertex_to_facet(bi); B.add_vertex_to_facet(ci);
            B.end_facet();
            if (this->any_vert_on_pos_side(f))
            {
                this->add_op_facet(a); this->add_op_facet(b); this->add_op_facet(c);
            }
        }

        // Done
        B.end_surface();
    }
};
void Slice::build_mesh(const Polyhedron3* P)
{
    assert(_mesh == nullptr || _mesh->empty());

    // Get all planes that will restrict the mesh
    std::vector<Plane3> planes = all_planes();

    // Find all of the polyhedron vertices that should be in the final result
    std::vector<P3CVertex> seeds;
    seeds.reserve(512);
    for (S3VertexDesc sv : svs)
    {
        for (const P3CVertex& v : (*S)[sv].vertices)
        {
            if (on_all_pos_sides(planes, v->point())) { seeds.push_back(v); }
        }
    }
    if (seeds.size() == 0)
    {
        // TODO: this is occurring somewhat frequently (sometimes there are no points in the entire mesh that match) - why?
        std::cerr << "ERROR: no seed found to start building mesh from" << std::endl;
        if (_mesh) { delete _mesh; _mesh = nullptr; }
        _mesh = new Polyhedron3();
        if (uncapped) { delete uncapped; uncapped = nullptr; }
        uncapped = new Polyhedron3();
        return;
    }

    // Get ready to cut up the mesh
    Polyhedron3* mesh = new Polyhedron3();
    bool first = true;

    // Quickly cut up with all planes simultaneously so that the whole process is much faster.
    // This step is completely optional but does cut the computation time in half.
    {
        QuickCutAtPlanes cut(P, planes, seeds);
        mesh->delegate(cut);
        assert(mesh->is_valid());

        // Quick-Cut doesn't seem amenable to closing the mesh itself so we do it here
        CGAL::set_halfedgeds_items_id(*mesh);
        for (auto he = mesh->halfedges_begin(), end = mesh->halfedges_end(); he != end; ++he)
        {
            if (he->is_border())
            {
                PMP::triangulate_hole(*mesh, he, null_output_iterator());
            }
        }
        assert(mesh->is_closed());

        // If a single connected component is formed than we can continue with the results of the
        // quick cut, otherwise we need to start over from the beginning. At the moment this never
        // seems to be triggered, but is a very fast check so we will leave it in.
        CGAL::set_halfedgeds_items_id(*mesh);
        if (PMP::connected_components(*mesh, P3_facet_int_map(get(boost::face_index, *mesh))) == 1)
        {
            first = false;
        }
        else { mesh->clear(); }
    }

    CGAL::set_halfedgeds_items_id(*mesh);

    // Cut up the mesh into a new mesh
    P3CVertex seed = pop(seeds);
    const Point3& p = seed->point();
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
        assert(mesh->is_valid());
        assert(mesh->is_closed());
        CGAL::set_halfedgeds_items_id(*mesh);
    }

    // Check mesh
    #ifdef _DEBUG // these checks are incredibly unlikely to fail and/or expensive to compute so usually don't do them
    if (!mesh->is_valid())         { std::cerr << "Warning: slice is not valid" << std::endl; }
    if (!mesh->is_closed())        { std::cerr << "Warning: slice is not closed" << std::endl; }
    if (!is_not_degenerate(mesh))  { std::cerr << "Warning: slice is degenerate" << std::endl; }
    if (!mesh->is_pure_triangle()) { std::cerr << "Warning: slice is not pure triangle" << std::endl; }
    #endif
    if (PMP::does_self_intersect(*mesh)) { std::cerr << "Warning: slice is self-intersecting" << std::endl; } // this one is expensive but necessary
    #ifdef _DEBUG
    else if (!PMP::does_bound_a_volume(*mesh)) { std::cerr << "Warning: slice does not bound a volume" << std::endl; }
    #endif

    // Set the class's mesh
    if (_mesh) { delete _mesh; }
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

    // TODO: instead of a fixed size this should be a bit dynamic so that large branch points can
    // get completely covered without forcing smaller ones to get overwhelmed.
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
    skeleton_enum_branches(S, [&, S, group_sz, bp_group_sz] (const std::vector<S3VertexDesc>& verts) mutable
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
