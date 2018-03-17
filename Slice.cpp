#include "Slice.hpp"

#include "Skeleton.hpp"
#include "Polyhedron3Utils.hpp"

#include <utility>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <boost/foreach.hpp>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
namespace PMP = CGAL::Polygon_mesh_processing;

// Polyhedron Property Maps
typedef boost::property_map<Polyhedron3, boost::face_index_t>::const_type P3_facet_index_map_t;
typedef boost::vector_property_map<size_t, P3_facet_index_map_t> P3_facet_int_map;

// Filtered Polyhedron3 Graph
typedef CGAL::Face_filtered_graph<Polyhedron3> Filtered_P3;

// Groups
typedef std::vector<std::vector<S3VertexDesc>> Groups;
typedef std::unordered_map<size_t, P3FacetSet> GroupFacets;

static void create_groups(const int group_sz, const Skeleton3* S, Groups& groups)
{
    // Creates groups of skeleton vertices. The branch points in the skeleton will always be in
    // groups of one more than the degree of the branch point, including the branch point itself and
    // a single neighboring point in each direction.
    //
    // Other groups will contain at most group_sz consecutive skeleton vertices. In the case
    // branches cannot be evenly divided by group_sz the groups will be made smaller while trying to
    // keep all of the groups roughly the same size (at most different by 1).

    groups.reserve(num_vertices(*S) / group_sz);
    // First add all of the branch points with their neighbors
    BOOST_FOREACH(auto sv, vertices(*S))
    {
        if (degree(sv, *S) > 2)
        {
            groups.push_back({{sv}});
            BOOST_FOREACH(auto e, out_edges(sv, *S))
            {
                groups.back().push_back(opposite(*S, e, sv));
            }
            groups.back().shrink_to_fit();
        }
    }
    // Then add all others
    skeleton_enum_branches(S, [&, S, group_sz] (const std::vector<S3VertexDesc>& verts) mutable
    {
        int start = (degree(verts.front(), *S) != 1)*2, end = (degree(verts.back(), *S) != 1)*2;
        size_t total = verts.size()-end-start;
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
        BOOST_FOREACH(size_t sz, sizes)
        {
            groups.push_back(std::vector<S3VertexDesc>(itr, itr + sz));
            itr += sz;
        }
        assert(itr + end == verts.end());
    });
}

static size_t __find_facets_to_skip(/*const*/ Polyhedron3* P, handle_map<P3CVertex, size_t>& v2grp, P3CVertexSet& skip_facets)
{
    // Fills in the skip_facets set with facets that should be skipped during the first pass of
    // determining the segments. Facets that should be skipped are either trivially not a
    // connected component (they have no neighbors of the same group) or because they connect
    // to other facets of the same group at single vertex and not an edge (which is not
    // supported by Face_filtered_graph).
    //
    // The v2grp parameter is the mapping of facet-vertices to group ids. Any facet-vertices to be
    // skipped will be set to (size_t)-1 in this map. Only facet-vertices will be added to the
    // skip_facets set. The return value is the number of values added to the skip_facets set.
    //
    // This internal function may need to be called multiple times to find all facets to skip.
    const size_t SKIPPED = (size_t)-1; // this becomes the group id of any skipped facet
    const size_t initially_skipped = skip_facets.size();
    for (P3Vertex v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v)
    {
        auto itr = v2grp.find(v);
        if (itr != v2grp.end())
        {
            if (itr->second == SKIPPED) { continue; }

            // Check if the facet is isolated
            int any_same = false; // true if any of the neighbors have the same group id as us
            FOR_EDGES_AROUND_VERTEX(v, e)
            {
                if (v2grp.at(e->prev()->opposite()->next()->vertex()) == itr->second) { any_same = true; break; }
            }
            if (!any_same) { skip_facets.insert(v); itr->second = SKIPPED; } // we are isolated
        }
        else
        {
            // Check if vertex has a non-continous set of group facets around it
            std::vector<size_t> around; // contains the group ids as we go around the vertex with repeats removed
            FOR_VERTICES_AROUND_VERTEX(v, u)
            {
                auto itr = v2grp.find(u);
                if (itr != v2grp.end() && (around.size() == 0 || around.back() != itr->second)) { around.push_back(itr->second); }
            }
            // Correct for the fact that it is a cycle around
            if (around.back() == around.front()) { around.erase(around.end()-1); }
            // around is now a list of hopefully unique values - check it
            std::sort(around.begin(), around.end()); // moves any non-adjecent duplicates next to each other
            for (size_t i = 1; i < around.size(); ++i)
            {
                if (around[i] != SKIPPED && around[i] == around[i-1])
                {
                    // Found a bad group id - mark it as skipped
                    FOR_VERTICES_AROUND_VERTEX(v, u)
                    {
                        auto itr = v2grp.find(u);
                        if (itr != v2grp.end() && around[i] == itr->second) { skip_facets.insert(u); itr->second = SKIPPED; }
                    }
                }
            }
        }
    }
    return skip_facets.size() - initially_skipped;
}

static void _find_facets_to_skip(/*const*/ Polyhedron3* P, const Skeleton3* S,
    const P3CVertexSet& facet_verts, const Groups& groups, P3CVertexSet& skip_facets)
{
    // Fills in the skip_facets set with facets that should be skipped during the first pass of
    // determining the segments. Facets that should be skipped are either trivially not a
    // connected component (they have no neighbors of the same group) or because they connect
    // to other facets of the same group at single vertex and not an edge (which is not
    // supported by Face_filtered_graph).

    // Create a map of mesh vertices to group ids
    handle_map<P3CVertex, size_t> v2grp;
    for(size_t gid = 0; gid < groups.size(); ++gid)
    {
        BOOST_FOREACH(S3VertexDesc sv, groups[gid])
        {
            BOOST_FOREACH(P3CVertex v, (*S)[sv].vertices)
            {
                if (facet_verts.count(v)) { v2grp.insert({{v, gid}}); }
            }
        }
    }

    // Call the internal function until no additional facets are found to skip
    while (__find_facets_to_skip(P, v2grp, skip_facets) != 0);
}

static void _primary_facets_for_groups(/*const*/ Polyhedron3* P, const Skeleton3* S,
    const P3CVertexSet& facet_verts, const Groups& groups, const P3CVertexSet& skip_facets, GroupFacets& group_facets)
{
    // Finds the facets on the mesh P that are associated with each of the groups of skeleton
    // vertices. It will skip any facets specified.
    //
    // This is accomplished by creating meshes for each group that includes all of the facets for
    // the mesh then selecting the facets in the largest connected component.

    Filtered_P3 P_filtered(*P, P3FacetSet());
    P3_facet_int_map component_map(get(boost::face_index, *P));
    P3FacetSet facets;
    for (size_t gid = 0; gid < groups.size(); ++gid)
    {
        facets.clear();
        BOOST_FOREACH(S3VertexDesc sv, groups[gid])
        {
            BOOST_FOREACH(P3Vertex v, (*S)[sv].vertices)
            {
                if (facet_verts.count(v) && !skip_facets.count(v))
                {
                    FOR_FACETS_AROUND_VERTEX(v, f) { facets.insert(f); }
                }
            }
        }
        if (facets.size() == 0) { continue; }

        // Perform connected-components analysis on just the select facets
        // Keep only the facets that are part of the largest connected-component
        P_filtered.set_selected_faces(facets);
        if (!P_filtered.is_selection_valid())
        {
            std::cerr << "Error: Selection not valid!" << std::endl;
        }
        size_t n_components = PMP::connected_components(P_filtered, component_map);
        std::vector<size_t> component_sizes(n_components, 0); // TODO: move out of loop
        BOOST_FOREACH (auto f, faces(P_filtered)) { component_sizes[component_map[f]] += 1; }
        size_t largest_component = 0, largest_component_size = component_sizes[0];
        for (size_t i = 1; i < n_components; ++i)
        {
            if (component_sizes[i] > largest_component_size)
            {
                largest_component_size = component_sizes[largest_component = i];
            }
        }
        BOOST_FOREACH (auto f, faces(P_filtered))
        {
            if (component_map[f] != largest_component) { facets.erase(f); }
        }

        // Done with this group, add the claimed facets for this group
        group_facets.insert({{gid, facets}});
    }
}

static void __get_neighbor_claims(P3Vertex v, handle_map<P3Vertex, size_t>& claimed, std::vector<size_t>& out)
{
    // Fills the out vector with all of the group id for the groups that claim the neighboring
    // facets. Unclaimed neighbors are not put in the vector. Groups that claim more than one
    // neighbor will show up multiple times in the vector.
    out.clear();
    //assert(facet_verts.count(v));
    FOR_EDGES_AROUND_VERTEX(v, e)
    {
        auto itr = claimed.find(e->prev()->opposite()->next()->vertex());
        if (itr != claimed.end()) { out.push_back(itr->second); }
        //else { assert(unclaimed.count(e->prev()->opposite()->next()->vertex())); }
    }
}

static bool _expand_facets_for_groups(/*const*/ Polyhedron3* P, P3VertexSet& unclaimed,
    handle_map<P3Vertex, size_t>& claimed, GroupFacets& group_facets, bool force_random=false)
{
    // Expands the facets claimed by the groups of skeletal vertices. This is done by taking any
    // facet which has a majority of neighbors claimed by the same group and assigning it to that
    // group. If no facets can be claimed in that manner, facets are claimed by the a randomly
    // selected group that claims one of their neighbors. The return value is true unless facets
    // had to be claimed randomly.
    //
    // This will have to be called many times until unclaimed.size() is 0.
    //
    // Special care must be taken for new claims that only have once neighbor already claimed by the
    // group as in these situations a group can form that doesn't form a proper set of facets.
    handle_map<P3Vertex, size_t> to_claim, to_claim_1;
    to_claim.reserve(unclaimed.size());
    std::vector<size_t> claims;
    if (!force_random)
    {
        for (auto v = unclaimed.begin(), end = unclaimed.end(); v != end; ++v)
        {
            __get_neighbor_claims(*v, claimed, claims);
            if (claims.size() == 1) { to_claim_1.insert({{*v, claims[0]}}); } // only one neighbor - special
            else if (claims.size() >= 2 && claims[0] == claims[1])
            {
                to_claim.insert({{*v, claims[0]}});
            }
            else if (claims.size() == 3 && (claims[0] == claims[2] || claims[1] == claims[2]))
            {
                to_claim.insert({{*v, claims[2]}});
            }
        }
    }

    bool randomly = to_claim.size() == 0 && to_claim_1.size() == 0;
    if (randomly)
    {
        // At this point randomly assign each unclaimed facet to one of its neighbors that is
        // claimed. If there are any facets with all 3 unclaimed neighbors then skip it and
        // we will get it the next time through.
        //size_t i = 0;
        for (auto f = unclaimed.begin(), end = unclaimed.end(); f != end; ++f)
        {
            __get_neighbor_claims(*f, claimed, claims);
            //if (claims.size() > 0) { to_claim_1.insert({{*f, claims[i++ % claims.size()]}}); }
            if (claims.size() > 0) { to_claim_1.insert({{*f, claims[Rand.get_int(0, claims.size())]}}); }
        }
    }

    // Actually claim the facets
    if (to_claim.size() == 0)
    {
        // Use to_claim_1 but with checking of the results after each addition that the facets are valid
        bool all_thrown_out = true;
        Filtered_P3 P_filtered(*P, handle_set<P3Facet>());
        BOOST_FOREACH(auto c, to_claim_1)
        {
            P3Vertex v = c.first;
            size_t gid = c.second;
            FOR_FACETS_AROUND_VERTEX(v, f) { group_facets[gid].insert(f); }
            P_filtered.set_selected_faces(group_facets[gid]);
            if (!P_filtered.is_selection_valid())
            {
                FOR_FACETS_AROUND_VERTEX(v, f) { group_facets[gid].erase(f); }
            }
            else
            {
                all_thrown_out = false;
                unclaimed.erase(v);
                claimed.insert({{v, gid}});
            }
        }
        if (all_thrown_out && !randomly)
        {
            // Try again with random forced
            return _expand_facets_for_groups(P, unclaimed, claimed, group_facets, true);
        }
    }
    else
    {
        // Use to_claim and quickly add all of the claims
        BOOST_FOREACH(auto c, to_claim)
        {
            P3Vertex v = c.first;
            size_t gid = c.second;
            FOR_FACETS_AROUND_VERTEX(v, f) { group_facets[gid].insert(f); }
            unclaimed.erase(v);
            claimed.insert({{v, gid}});
        }
    }

    return randomly;
}

static void __declaim(P3VertexSet& unclaimed, handle_map<P3Vertex, size_t>& claimed, GroupFacets& group_facets)
{
    // Goes through each unclaimed vertex and "de-claims" all of its neighbors.
    P3VertexSet declaim;
    BOOST_FOREACH (auto v, unclaimed)
    {
        FOR_EDGES_AROUND_VERTEX(v, e) { declaim.insert(e->prev()->opposite()->next()->vertex()); }
    }
    BOOST_FOREACH (auto v, declaim)
    {
        unclaimed.insert(v);
        auto itr = claimed.find(v);
        if (itr != claimed.end())
        {
            FOR_FACETS_AROUND_VERTEX(v, f) { group_facets[itr->second].erase(f); }
            claimed.erase(itr);
        }
    }
}

static void facets_for_groups(/*const*/ Polyhedron3* P, const Skeleton3* S,
    const P3CVertexSet& facet_verts, const Groups& groups, GroupFacets& group_facets, bool verbose=false)
{
    // Assigns every facet in P to a group based on the groups given. This makes sure that the
    // facets assigned to a group for a single, manifold, connected component and that every facet
    // in P is placed into a group.

    P3CVertexSet skip_facets;
    _find_facets_to_skip(P, S, facet_verts, groups, skip_facets);

    // Apply first pass to get the majority of facets claimed by a group
    _primary_facets_for_groups(P, S, facet_verts, groups, skip_facets, group_facets);

    // Setup the unclaimed and claimed set/map
    P3VertexSet unclaimed;
    for (auto v = P->vertices_begin(), end = P->vertices_end(); v != end; ++v)
    {
        // Convert from P3CVertex to P3Vertex...
        if (facet_verts.count(v)) { unclaimed.insert(v); }
    }
    handle_map<P3Vertex, size_t> claimed;
    BOOST_FOREACH(auto& gid_facets, group_facets)
    {
        BOOST_FOREACH(auto f, gid_facets.second)
        {
            P3Vertex v = f->halfedge()->vertex();
            unclaimed.erase(v);
            claimed.insert({{v, gid_facets.first}});
        }
    }

    // Display information
    if (verbose)
    {
        std::cout << unclaimed.size() << " (" << (unclaimed.size() * 100.0 * 3 / P->size_of_facets()) << "%) unclaimed facets after connected components" << std::endl;
    }

    int no_progress_iterations = 0, num_rewinds = 0;
    for (size_t iterations = 1; unclaimed.size(); ++iterations)
    {
        if (no_progress_iterations >= 2)
        {
            // We have a problem that seems unresolvable with the current setup, rewind a little bit
            // We will go through each unclaimed vertex and "de-claim" all of its neighbors
            // If this happens a bunch we get more aggressive with rewinding
            if (num_rewinds > 9)
            {
                std::cerr << "ERROR: Expansion of group facets is stuck with no hope of recovering..." << std::endl;
                throw std::domain_error("expansion stuck");
            }
            std::cerr << "WARNING: Expansion of group facets seems to be stuck so rewinding... (" << unclaimed.size() << " -> ";
            for (int i = 0; i < (num_rewinds/3) + 1; ++i) { __declaim(unclaimed, claimed, group_facets); }
            std::cerr << unclaimed.size() << ")" << std::endl;
            num_rewinds++;
        }

        // Perform the expansion of groups
        size_t before = unclaimed.size();
        bool randomly = _expand_facets_for_groups(P, unclaimed, claimed, group_facets);
        no_progress_iterations = (before == unclaimed.size()) ? no_progress_iterations + 1 : 0;

        // Display information
        if (verbose)
        {
            std::cout << unclaimed.size() << " (" << (unclaimed.size() * 100.0 * 3 / P->size_of_facets()) << "%) unclaimed facets after " << iterations << " iterations expansion";
            if (randomly) { std::cout << " (resorted to random assignment)"; }
            std::cout << std::endl;
        }
    }
}

Slices slice(const int group_sz, const Polyhedron3* P, const Skeleton3* S,
    const P3CVertexSet& facet_verts, bool verbose)
{
	Polyhedron3* P_ = const_cast<Polyhedron3*>(P); // the Face_filtered_graph doesn't cooperate well with the const...
	
    Groups groups;
    create_groups(group_sz, S, groups);

    GroupFacets group_facets;
    facets_for_groups(P_, S, facet_verts, groups, group_facets, verbose);

    // Extract each slice and pair it with its group of skeleton vertices
    Filtered_P3 P_filtered(*P_, handle_set<P3Facet>());
    Slices slices;
    slices.reserve(groups.size());
    BOOST_FOREACH(auto& gid_facets, group_facets)
    {
        if (gid_facets.second.size() == 0) { continue; }
        // Extract the facets into a new mesh
        P_filtered.set_selected_faces(gid_facets.second);
        if (!P_filtered.is_selection_valid()) { std::cerr << "ERROR: An invalid facet arrangement was generated for a group" << std::endl; }
        Polyhedron3* slc = new Polyhedron3();
        CGAL::copy_face_graph(P_filtered, *slc);
        // Add the slice
        slices.push_back(std::make_pair(groups[gid_facets.first], slc));
    }
    slices.shrink_to_fit();
    return slices;
}
