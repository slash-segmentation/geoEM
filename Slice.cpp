#include "Slice.hpp"

#include "Skeleton.hpp"
#include "Polyhedron3Utils.hpp"

#include <vector>
#include <unordered_set>
#include <iostream>
#include <iterator>

#include <boost/foreach.hpp>

#include <CGAL/boost/graph/properties.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
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
struct null_output_iterator : std::iterator<std::output_iterator_tag, null_output_iterator>
{
	// An output iterator that discards all data, useful for PMP::triangulate_hole.
	// From: https://stackoverflow.com/questions/335930/discarding-the-output-of-a-function-that-needs-an-output-iterator
    template<typename T> void operator=(T const&) { }
    null_output_iterator & operator++() { return *this; }
    null_output_iterator operator++(int) { return *this; }
    null_output_iterator & operator*() { return *this; }
};
Slice* Slice::capped() const
{
	// Creates a new slice that has all holes capped using PMP::triangulate_hole()
	Slice* slc = new Slice(this);
	for (auto he = slc->_mesh->halfedges_begin(); he != slc->_mesh->halfedges_end(); ++he)
	{
		if (he->is_border()) { PMP::triangulate_hole(*slc->_mesh, he, null_output_iterator()); }
	}

	// Check mesh
	if (!slc->_mesh->is_valid())         { std::cerr << "Warning: slice is not valid" << std::endl; }
	if (!slc->_mesh->is_closed())        { std::cerr << "Warning: slice is not closed" << std::endl; }
	if (!is_not_degenerate(slc->_mesh))  { std::cerr << "Warning: slice is degenerate" << std::endl; }
	if (PMP::does_self_intersect(*slc->_mesh)) { std::cerr << "Warning: slice is self-intersecting" << std::endl; }
	//if (!PMP::does_bound_a_volume(*slc->_mesh)) { std::cerr << "Warning: slice does not bound a volume" << std::endl; } // crashes in some cases...
	if (!slc->_mesh->is_pure_triangle()) { std::cerr << "Warning: slice is not pure triangle" << std::endl; }

	return slc;
}

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
    handle_set<P3CFacet> finished, stack;

    inline bool on_pos_side(const P3CVertex& v) const { return this->h.has_on_positive_side(v->point()); }
    inline bool on_neg_side(const P3CVertex& v) const { return this->h.has_on_negative_side(v->point()); }
    inline bool on_pos_side(const P3CHalfedge& he) const { return this->h.has_on_positive_side(he->vertex()->point()); }
    inline bool on_neg_side(const P3CHalfedge& he) const { return this->h.has_on_negative_side(he->vertex()->point()); }

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
        for (auto v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v) { vertex_count += this->on_pos_side(v); }
        B.begin_surface(vertex_count, facet_count, this->P->size_of_halfedges());

        // Start processing based on the given seed vertex
        FOR_EDGES_AROUND_VERTEX(this->v, he) { if (!he->is_border()) { this->stack.insert(he->facet()); } }
        std::vector<size_t> v_inds;
        while (!this->stack.empty())
        {
            P3CFacet f = *this->stack.begin();
            this->stack.erase(this->stack.begin());
            this->finished.insert(f);
            size_t n_neg = this->n_verts_on_neg_side(f);
            assert(n_neg < 3);
            if (n_neg == 0)
            {
                // Add the entire facet
                FOR_EDGES_AROUND_FACET(f, he)
                {
                    v_inds.push_back(this->get_pt_id(he));
                    this->add_op_facet(he);
                }
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
                    v_inds.push_back(a); v_inds.push_back(b);
                    const P3CHalfedge& he_cross = this->on_pos_side(he->next()) ? he->next() : he;
                    v_inds.push_back(this->get_int_pt_id(he_cross));
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
                    // Add the first facet
                    v_inds.push_back(b); v_inds.push_back(c); v_inds.push_back(d);
                    B.add_facet(v_inds.begin(), v_inds.end());
                    v_inds.clear();
                    // Add the second facet
                    v_inds.push_back(b); v_inds.push_back(d); v_inds.push_back(a);
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
                // The vertex already on the positive side
                v_inds.push_back(this->get_pt_id(he));
                // The vertices on the plane
                v_inds.push_back(this->get_int_pt_id(he->next()));
                v_inds.push_back(this->get_int_pt_id(he));
                // Add neighboring facets
                this->add_op_facet(he);
                this->add_op_facet(he->next());
            }
            if (v_inds.size() != 0)
            {
                B.add_facet(v_inds.begin(), v_inds.end());
                v_inds.clear();
            }
        }
        B.end_surface();
    }
};
void Slice::build_mesh(const Polyhedron3* P)
{
	assert(_mesh->empty());

	// Get all planes that will restrict the mesh
	std::vector<Plane3> planes = all_planes();

	// Find all of the polyhedron vertices that should be in the final result
	std::unordered_map<Point3, P3CVertex, boost::hash<Point3>> seeds;
	for (S3VertexDesc sv : svs)
    {
        for (P3CVertex v : (*S)[sv].vertices)
        {
			if (on_all_pos_sides(planes, v->point()))
			{
				seeds.insert({{v->point(), v}});
			}
		}
	}

	// Cut up the mesh into a new mesh
	Polyhedron3* mesh = new Polyhedron3();
	bool first = true;
	Point3 p; P3CVertex seed;
	std::tie(p, seed) = pop(seeds);
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

	// Update the seeds, removing any that are in the mesh
	for (auto v = mesh->vertices_begin(), end = mesh->vertices_end(); v != end; ++v)
	{
		seeds.erase(v->point());
	}
	if (!seeds.empty()) { std::cerr << "Error: not all points ended up in final mesh" << std::endl; }

	// Check mesh
	if (!mesh->is_valid())         { std::cerr << "Warning: slice is not valid" << std::endl; }
    if (!is_not_degenerate(mesh))  { std::cerr << "Warning: slice is degenerate" << std::endl; }
    if (PMP::does_self_intersect(*mesh)) { std::cerr << "Warning: slice is self-intersecting" << std::endl; }
    if (!mesh->is_pure_triangle()) { std::cerr << "Warning: slice is not pure triangle" << std::endl; }

	// Set the class's mesh
	delete _mesh;
	_mesh = mesh;
}

///////////////////////////////////////////////////////////////////////////////
// Group creation
///////////////////////////////////////////////////////////////////////////////
static void create_groups(const int group_sz, const Skeleton3* S, Groups& groups)
{
    // Creates groups of skeleton vertices. The branch points in the skeleton will always be in
    // groups of one more than the degree of the branch point, including the branch point itself and
    // a single neighboring point in each direction.
    //
    // Other groups will contain at most group_sz consecutive skeleton vertices. In the case
    // branches cannot be evenly divided by group_sz the groups will be made smaller while trying to
    // keep all of the groups roughly the same size (at most different by 1).

    static const int BP_DISTANCE = 2; // TODO: based on group_sz?

    groups.reserve(num_vertices(*S) / group_sz);
    // First add all of the branch points with their neighbors
    BOOST_FOREACH(auto sv, vertices(*S))
    {
        if (degree(sv, *S) > 2)
        {
            std::unordered_set<S3VertexDesc> svs;
            svs.reserve(degree(sv, *S)*BP_DISTANCE + 1);
            std::unordered_set<S3VertexDesc> next;
            next.insert(sv);
            for (int distance = 1; distance <= BP_DISTANCE; ++distance)
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
        int start = (degree(verts.front(), *S) != 1)*(BP_DISTANCE+1);
        int end = (degree(verts.back(), *S) != 1)*(BP_DISTANCE+1);
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
    Slice::set_neighbors(slices);

	// Finally the meshes for each of the slices can be built
	for (Slice* slc : slices) { slc->build_mesh(P); }

    return slices;
}

///////////////////////////////////////////////////////////////////////////////
// Function that takes slices and creates capped version of all of them
///////////////////////////////////////////////////////////////////////////////
Slices cap_slices(const Slices& slices)
{
    Slices capped;
    capped.reserve(slices.size());
    for (Slice* slc : slices)
    {
        capped.push_back(slc->capped());
    }
    Slice::set_neighbors(capped);
    return capped;
}
