#include "Slice.hpp"

#include "Skeleton.hpp"
#include "PolyBuilder.hpp"
#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <numeric>

#include <math.h>

#include <boost/foreach.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
namespace PMP = CGAL::Polygon_mesh_processing;

// Incremental builder
typedef CGAL::Polyhedron_incremental_builder_3<Polyhedron3::HalfedgeDS> Builder;

// Groups
typedef std::vector<std::unordered_set<S3VertexDesc>> Groups;

template <class Container>
typename Container::value_type pop(Container& c)
{
    typename Container::value_type val = *c.begin();
    c.erase(c.begin());
    return val;
}

// Utilities
inline static bool has_all_on_pos_side(const Plane3& h, const std::unordered_set<S3VertexDesc>& svs, const Skeleton3* S)
{
    for (S3VertexDesc sv : svs)
    {
        for (const P3CVertex& v : (*S)[sv].vertices)
        {
            if (!h.has_on_positive_side(v->point())) { return false; }
        }
    }
    return true;
}
inline static bool has_all_on_pos_sides(const std::vector<Plane3>& planes, const Point3& p)
{
    for (auto& h : planes) { if (!h.has_on_positive_side(p)) { return false; } }
    return true;
}
inline static std::vector<P3CVertex> get_all_on_pos_sides(const std::vector<Plane3>& planes, const std::unordered_set<S3VertexDesc>& svs, const Skeleton3* S)
{
    std::vector<P3CVertex> verts;
    for (S3VertexDesc sv : svs)
    {
        for (const P3CVertex& v : (*S)[sv].vertices)
        {
            if (has_all_on_pos_sides(planes, v->point())) { verts.push_back(v); }
        }
    }
    return verts;
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
    _all_planes.insert(_all_planes.end(), end_planes.begin(), end_planes.end());
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
S3VertexDesc Slice::endpoint() const
{
    for (S3VertexDesc sv : this->svs)
    {
        if (boost::degree(sv, *this->S) == 1) { return sv; }
    }
    return Skeleton3::null_vertex();
}
size_t Slice::is_endpoint() const
{
    for (S3VertexDesc sv : this->svs)
    {
        if (boost::degree(sv, *this->S) == 1) { return true; }
    }
    return false;
}

bool Slice::is_relevant(const Plane3& h, const P3Side& side) const
{
    for (size_t i = 0; i < this->end_planes.size(); ++i)
    {
        S3VertexDesc sv = neighbor_sv(this->end_verts[i]);
        Point3 p = CGAL::midpoint((*S)[sv].point, (*S)[this->end_verts[i]].point);

        // Find the line formed by the intersection of h and g
        auto result = CGAL::intersection(h, this->end_planes[i]);
        //if (!result) { continue; } // plane is parallel to the current plane, ???? TODO
        const Line3* l = boost::get<Line3>(&*result);
        if (!l) { continue; } // identical plane, not relevant

        // If the given point projected onto the line is inside the mesh then it is relevant
        if (side(l->projection(p)) != CGAL::ON_UNBOUNDED_SIDE) { return true; }
    }
    return false;
}

void Slice::set_neighbors(std::vector<Slice*> slices, const Polyhedron3* P)
{
    // Set all of the neighbors of all of the slices. This uses information from each slice to
    // determine which other slices are neighbors. This also calculates the extra planes.

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

    // Calculate the extra planes
    P3Side side(*P);
    for (Slice* slc : slices)
    {
        // For each branch point
        if (slc->deg > 2)
        {
            // Work down the branches copying the relevant auxilary planes
            for (size_t i = 0; i < slc->end_neighbors.size(); ++i)
            {
                std::vector<Plane3> planes = slc->aux_planes[i];
                Slice *last = slc, *nghbr = slc->end_neighbors[i];
                size_t dist = 0;
                while (!planes.empty() && nghbr && nghbr->deg <= 2)
                {
                    if (++dist > 2)
                    {
                        // Remove any planes that are no longer relevant
                        for (auto itr = planes.begin(); itr != planes.end(); )
                        {
                            if (nghbr->is_relevant(*itr, side)) { ++itr; }
                            else { itr = planes.erase(itr); }
                        }
                    }
                    // Add the planes
                    nghbr->_all_planes.insert(nghbr->_all_planes.end(), planes.begin(), planes.end());
                    // Move down to the next slice
                    Slice* next = nghbr->other_neighbor(last);
                    last = nghbr; nghbr = next;
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Slice class meshing methods
///////////////////////////////////////////////////////////////////////////////

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


class CutAtPlane : public CGAL::Modifier_base<Polyhedron3::HalfedgeDS>
{
    const Polyhedron3* P;
    const Plane3 h;
    const std::vector<P3CVertex>& seeds;
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
    CutAtPlane(const Polyhedron3* P, Plane3 h, const std::vector<P3CVertex>& seeds) : P(P), h(h), seeds(seeds)
    {
        assert(P->is_closed());
        assert(P->is_pure_triangle());
        #ifdef _DEBUG
        for (const P3CVertex& v : seeds) { assert(h.has_on_positive_side(v->point())); }
        #endif
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

        // Start processing based on the given seeds, pre-adding them as vertices
        this->lookup.resize(P->size_of_vertices(), (size_t)-1);
        this->edge_lookup.resize(P->size_of_halfedges(), (size_t)-1);
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

inline static std::vector<P3CVertex> find_vertices(const Polyhedron3* P, std::vector<P3CVertex>& old)
{
    std::unordered_map<Point3, P3CVertex, boost::hash<Point3>> map;
    map.reserve(P->size_of_vertices());
    for (P3CVertex v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v) { map.insert({{v->point(), v}}); }
    std::vector<P3CVertex> vs;
    vs.reserve(old.size());
    for (P3CVertex v : old) { vs.push_back(map.at(v->point())); }
    return vs;
}

void Slice::build_mesh(const Polyhedron3* P)
{
    assert(_mesh == nullptr || _mesh->empty());

    // Find all of the vertices mapped to skeletal vertices that should be in the final result
    // We also include skeletal vertices in neighboring slices since there may be many "exchanges"
    // especially near branch points.
    std::unordered_set<S3VertexDesc> all_svs(this->svs);
    // less aggressive: for (S3VertexDesc sv : this->end_verts) { all_svs.insert(this->neighbor_sv(sv)); }
    for (Slice* neighbor : this->end_neighbors) { all_svs.insert(neighbor->svs.begin(), neighbor->svs.end()); }
    std::vector<P3CVertex> seeds = get_all_on_pos_sides(_all_planes, all_svs, this->S);
    if (seeds.size() == 0)
    {
        // TODO: this is can still occur (sometimes there are no points in the entire mesh that match)
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
        QuickCutAtPlanes cut(P, _all_planes, seeds);
        mesh->delegate(cut);
        assert(mesh->is_valid());

        // Quick-Cut doesn't seem amenable to closing the mesh itself so we do it here
        CGAL::set_halfedgeds_items_id(*mesh);
        triangulate_holes(mesh);

        // If a single connected component is formed than we can continue with the results of the
        // quick cut, otherwise we need to start over from the beginning. At the moment this never
        // seems to be triggered, but is a very fast check so we will leave it in.
        if (is_single_component(mesh)) { first = false; } else { mesh->clear(); }
    }

    // Cut up the mesh into a new mesh
    for (auto& h : _all_planes)
    {
        if (first)
        {
            // First cut - uses original mesh (this is only used when quick-cut wasn't used)
            CutAtPlane cut(P, h, seeds);
            mesh->delegate(cut);
            first = false;
        }
        else
        {
            // All other cuts - use previous result
            std::vector<P3CVertex> seeds2 = find_vertices(mesh, seeds);
            Polyhedron3* temp = new Polyhedron3();
            CutAtPlane cut(mesh, h, seeds2);
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
    if (!is_not_degenerate(mesh))  { std::cerr << "Warning: slice is degenerate" << std::endl; }
    if (!mesh->is_closed())        { std::cerr << "Warning: slice is not closed" << std::endl; }
    else if (!PMP::is_outward_oriented(*mesh)) { std::cerr << "Warning: slice is not outward oriented" << std::endl; }
    if (!mesh->is_pure_triangle()) { std::cerr << "Warning: slice is not pure triangle" << std::endl; }
    #endif
    if (PMP::does_self_intersect(*mesh)) { std::cerr << "Warning: slice is self-intersecting" << std::endl; } // this one is expensive but necessary
    #ifdef _DEBUG
    else if (!PMP::does_bound_a_volume(*mesh)) { std::cerr << "Warning: slice does not bound a volume" << std::endl; }
    if (!is_single_component(*mesh)) { std::cerr << "Warning: slice is not single connected component" << std::endl; }
    #endif

    // Set the class's mesh
    if (_mesh) { delete _mesh; }
    _mesh = mesh;

    // Create uncapped version of mesh
    std::unordered_set<P3Facet> facets;
    for (auto f = mesh->facets_begin(), end = mesh->facets_end(); f != end; ++f)
    {
        // This would add capping-only facets (which is an interesting view)
        //for (auto& h : _all_planes) { if (f->plane() == h) { facets.insert(f); break; } }
        // This adds original facets only
        bool any_match = false;
        for (auto& h : _all_planes) { if (f->plane() == h) { any_match = true; break; } }
        if (!any_match) { facets.insert(f); }
    }
    CGAL::Face_filtered_graph<Polyhedron3> uncapped_graph(*mesh, facets);
    if (uncapped) { delete uncapped; uncapped = nullptr; }
    CGAL::copy_face_graph(uncapped_graph, *(uncapped = new Polyhedron3()));
}

///////////////////////////////////////////////////////////////////////////////
// Group creation
///////////////////////////////////////////////////////////////////////////////
// Info about the group formed around a branch point
struct BPInfo
{
    size_t index;
    S3VertexDesc sv;
    std::unordered_map<S3VertexDesc, std::vector<S3VertexDesc>> branches;
};
inline static Kernel::FT min_dist2_to_skel(const S3VertexDesc sv, const Skeleton3* S)
{
    // Calculate the minimum squared distance of every mesh vertex associated with the skeleton
    // vertex to the skeleton vertex.
    Kernel::FT min_dist2 = -1;
    const Point3& p = (*S)[sv].point;
    BOOST_FOREACH (P3CVertex v, (*S)[sv].vertices)
    {
        Kernel::FT dist2 = CGAL::squared_distance(v->point(), p);
        if (min_dist2 == -1 || dist2 < min_dist2) { min_dist2 = dist2; }
    }
    return min_dist2;
}
inline static double sv_length(const S3VertexDesc sv, const Skeleton3* S)
{
    // Calculates the skeletal "length" of a vertex by going to the midpoints of the neighboring
    // vertices.
    double len = 0;
    const Point3& p = (*S)[sv].point;
    BOOST_FOREACH (auto e, out_edges(sv, *S))
    {
        len += distance(p, (*S)[opposite(*S, e, sv)].point);
    }
    return len / 2; // only go to midpoints
}
template <class Callback>
inline static void for_each_perm(size_t n1, size_t n2, int val1, int val2, Callback cb)
{
    // Iterates through each permutation that has n1 or val1 and n2 or val2 in it. For each
    // permutation the callback is called with a vector containing the permutation.
    std::vector<int> cur;
    cur.reserve(n1 + n2);
    std::function<void(size_t,size_t)> recurse = [&, val1, val2] (size_t n1, size_t n2)
    {
        if (n1 == 0)
        {
            // All remaining n2 values are val2
            cur.insert(cur.end(), n2, val2);
            cb(cur);
            cur.erase(cur.end()-n2, cur.end());
        }
        else if (n2 == 0)
        {
            // All remaining n1 values are val1
            cur.insert(cur.end(), n1, val1);
            cb(cur);
            cur.erase(cur.end()-n1, cur.end());
        }
        else
        {
            // Need to recurse with both a val1 and val2 at the front
            cur.push_back(val1);
            recurse(n1 - 1, n2);
            cur.back() = val2;
            recurse(n1, n2 - 1);
            cur.pop_back();
        }
    };

}
inline static double ptp_of_groups(const std::vector<double>& values, const std::vector<int>& sizes)
{
    // Calculate the peak-to-peak distance (max-min) for each of the groups of values given by
    // sizes. E.g. if the first values in sizes is 4, then the first 4 values are grouped into one.
    double min = INFINITY, max = -INFINITY;
    size_t i = 0;
    for (int size : sizes)
    {
        double val = std::accumulate(values.begin()+i, values.begin()+i+size, 0.0);
        if (val < min) { min = val; }
        if (val > max) { max = val; }
        i += size;
    }
    return max - min;
}
static void create_groups(const int group_sz, const Skeleton3* S, Groups& groups)
{
    // Creates groups of skeleton vertices. The branch points in the skeleton will always be in put
    // into groups based on how large the mesh is within their region, with at least degree + 1
    // skeleton vertices but possibly many more. Other groups will contain at most group_sz
    // consecutive skeleton vertices. In the case branches cannot be evenly divided by group_sz the
    // groups will be made smaller while trying to keep all of the groups roughly the same size (at
    // most different by 1).

    groups.reserve(num_vertices(*S) / group_sz);

    // First add all of the branch points with their neighbors
    std::unordered_map<S3VertexDesc, BPInfo> bps;
    BOOST_FOREACH(auto sv, vertices(*S))
    {
        if (degree(sv, *S) > 2)
        {
            BPInfo bp;
            bp.index = groups.size();
            bp.sv = sv;
            std::unordered_set<S3VertexDesc> svs;
            svs.insert(sv);
            BOOST_FOREACH(auto e, out_edges(sv, *S))
            {
                std::vector<S3VertexDesc> branch;
                S3VertexDesc prev = sv, cur = opposite(*S, e, prev);
                const Point3& p = (*S)[sv].point;
                Kernel::FT min_dist2 = min_dist2_to_skel(cur, S);
                do
                {
                    branch.push_back(cur);
                    svs.insert(cur);
                    S3VertexDesc temp = cur; cur = next_vertex(*S, prev, cur); prev = temp;
                }
                while (degree(cur, *S) <= 2 && CGAL::squared_distance((*S)[cur].point, p) <= min_dist2);
                bp.branches.insert({{branch[0], branch}});
            }
            bps.insert({{sv, bp}});
            groups.push_back(svs);
        }
    }
    // Then add all others
    skeleton_enum_branches(S, [&, S, group_sz] (const std::vector<S3VertexDesc>& verts)
    {
        // Is the start/end a branch point and how many skeleton vertices does that BP claim?
        BPInfo* bp1 = degree(verts.front(), *S) > 1 ? &bps[verts.front()] : nullptr;
        BPInfo* bp2 = degree(verts.back(), *S) > 1 ? &bps[verts.back()] : nullptr;
        int start = bp1 ? (bp1->branches.at(verts[1]).size() + 1) : 0;
        int end = bp2 ? (bp2->branches.at(verts[verts.size()-2]).size() + 1) : 0;
        // Total number of skeleton vertices on this branch
        ssize_t total = (ssize_t)verts.size()-end-start;
        if (total <= 0) // no groups to add since branch is too short
        {
            if (total < 0 && bp1 && bp2)
            {
                // We have a branch that is so short that the branch points overlap
                std::unordered_set<S3VertexDesc>& g1 = groups[bp1->index];
                std::unordered_set<S3VertexDesc>& g2 = groups[bp2->index];
                std::vector<S3VertexDesc> overlap;
                overlap.reserve(-2*total);
                for (S3VertexDesc sv : verts)
                {
                    if (g1.count(sv) && g2.count(sv)) { overlap.push_back(sv); }
                }
                // Half to each (middle point to smaller or first if equal)
                size_t n = overlap.size() / 2 + (overlap.size() % 2 == 1 && g2.size() > g1.size());
                for (size_t i = 0; i < n; ++i) { g2.erase(overlap[i]); }
                for (size_t i = n; i < overlap.size(); ++i) { g1.erase(overlap[i]); }
            }
            return;
        }
        // The number of skeleton vertices in each group
        std::vector<int> sizes(total/group_sz, group_sz);
        // Distribute the remainder
        int rem = total%group_sz;
        if (rem != 0)
        {
            // Need one more group for the remainder
            // We also need to distribute the remainder losses around
            sizes.push_back(group_sz);
            size_t n = sizes.size(), i = n - 1;
            for (size_t to_remove = group_sz - rem; to_remove > 0; to_remove -= 1)
            {
                sizes[i] -= 1;
                if (i == 0) { i = n; }
                i -= 1;
            }
            // Redistribute smaller groups to be where length is greater
            if (i != n-1)
            {
                // There are two different group sizes, attempt to reduce dispersion of data
                int sz1 = sizes[i+1], sz2 = sizes[i]; // sz1 = sz2 - 1, 2 <= sz2 <= group_sz
                size_t n1 = n-i-1, n2 = i+1; // n1 + n2 = n
                // Calculate the metric we will be minimizing
                std::vector<double> values;
                values.reserve(total);
                for (size_t i = start; i < verts.size()-end; ++i)
                {
                    //double val = sv_length(verts[i], S)); // "length" of each sv
                    double val = (*S)[verts[i]].vertices.size(); // number of mapped mesh vertices
                    values.push_back(val);
                }
                // Go through each permutation for the sizes and find the one that reduces the
                // peak-to-peak distance of the metric across the branch.
                double min_val = INFINITY;
                for_each_perm(n1, n2, sz1, sz2, [&] (const std::vector<int>& szs)
                {
                    double val = ptp_of_groups(values, szs);
                    if (val < min_val) { min_val = val; sizes = szs; }
                });
            }
        }
        // Add the groups of the correct sizes
        auto itr = verts.begin() + start;
        for (size_t sz : sizes)
        {
            groups.push_back(std::unordered_set<S3VertexDesc>(itr, itr + sz));
            itr += sz;
        }
        assert(itr + end == verts.end());
    });
}

///////////////////////////////////////////////////////////////////////////////
// Function that takes a mesh and skeleton and creates the slices
///////////////////////////////////////////////////////////////////////////////
Slices slice(const int group_sz, const Polyhedron3* P, const Skeleton3* S)
{
    // Create the groups
    Groups groups;
    create_groups(group_sz, S, groups);

    // Create each slice
    Slices slices;
    slices.reserve(groups.size());
    for (auto& group : groups) { slices.push_back(new Slice(S, group.begin(), group.end())); }

    // Once all slices are created setup neighbors
    Slice::set_neighbors(slices, P);

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
