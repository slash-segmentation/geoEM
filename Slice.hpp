#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include <utility>
#include <vector>
#include <unordered_set>

#include <CGAL/Side_of_triangle_mesh.h>
typedef CGAL::Side_of_triangle_mesh<Polyhedron3, Kernel> P3Side;

class Slice
{
    // Represents a slice of the mesh corresponding to several skeleton vertices. The slice()
    // function will create the slices based on a Polyhedron3 and Skeleton3. If creating manually
    // the general order is Slice constructor for each slice, then set_neighbors(), then finally
    // build_mesh() for each slice.

private:
    // TODO: the ends of a slice can only have one outgoing edge each for the moment (i.e. branch
    // points cannot be at the ends and unless the slice is degree 1 it cannot have a single
    // vertex).
    const Skeleton3* S;
    std::unordered_set<S3VertexDesc> svs;
    size_t deg = 0;
    // All of the end_* are in the same order
    std::vector<S3VertexDesc> end_verts;
    std::vector<Plane3> end_planes;
    std::vector<Slice*> end_neighbors;
    // Branch points cause additional planes for each end to be passed on along the branches
    // Same order as end_* if available (i.e a branch point)
    std::vector<std::vector<Plane3>> aux_planes;
    std::vector<Plane3> _all_planes; // planes plus aux_planes passed along branches
    Polyhedron3* _mesh = nullptr;
    Polyhedron3* uncapped = nullptr;

    void init();

    // Used by copy constructor to either take a, copy b, or just give nullptr.
    inline static Polyhedron3* get_mesh(Polyhedron3* a, Polyhedron3* b) { return a ? a : b ? new Polyhedron3(*b) : nullptr; }

    // Used by set_neighbors, gets a neighbor that is not the given neighbor
    inline Slice* other_neighbor(Slice* nghbr)
    {
        for (Slice* slc : this->end_neighbors) { if (slc != nghbr) { return slc; } }
        return nullptr;
    }

    // Used to calculate _all_planes
    bool is_relevant(const Plane3& h, const P3Side& side) const;

    // Gets the skeletal vertex that is adjacent to sv but in one of the neighbors
    inline S3VertexDesc neighbor_sv(S3VertexDesc sv) const
    {
        BOOST_FOREACH (auto e, out_edges(sv, *S))
        {
            S3VertexDesc sv2 = opposite(*S, e, sv);
            if (!this->svs.count(sv2)) { return sv2; }
        }
        throw std::invalid_argument("no neighbor sv");
    }

public:
    // Create a slice from the skeleton vertices given by the iterators. If a mesh or uncapped mesh
    // is provided then it will become owned by this object (i.e. will be tied to the lifetime of
    // this object).
    template <class ForwardIterator>
    inline Slice(const Skeleton3* S, ForwardIterator begin, ForwardIterator end,
                 Polyhedron3* mesh=nullptr, Polyhedron3* uncapped=nullptr)
        : S(S), svs(begin, end), _mesh(mesh), uncapped(uncapped) { init(); }

    // Copy constructors copy the mesh and uncapped unless given a new mesh/uncapped in which case
    // this becomes the owner of those meshes like the regular constructor. Neighbors are not
    // copied, you must still call set_neighbors.
    inline Slice(const Slice* slc, Polyhedron3* mesh=nullptr, Polyhedron3* uncapped=nullptr)
        : S(slc->S), svs(slc->svs), deg(slc->deg),
        end_verts(slc->end_verts), end_planes(slc->end_planes), end_neighbors(slc->end_neighbors.size(), nullptr), aux_planes(slc->aux_planes),
        _mesh(get_mesh(mesh, slc->_mesh)), uncapped(get_mesh(uncapped, slc->uncapped)) { }
    Slice(const Slice&) = delete;
    Slice& operator=(const Slice&) = delete;
    inline ~Slice() { if (_mesh) { delete _mesh; } if (uncapped) { delete uncapped; } }

    // Setup functions
    static void set_neighbors(std::vector<Slice*> slices, const Polyhedron3* P);
    void build_mesh(const Polyhedron3* P);

    // Access to basic properties (a few of these are calculated, but most just return a field)
    double length() const;
    inline size_t degree() const { return this->deg; }
    inline const Skeleton3* skeleton() const { return this->S; }
    inline const std::unordered_set<S3VertexDesc>& skeleton_vertices() const { return this->svs; }
    inline const std::vector<Slice*>& neighbors() const { return this->end_neighbors; }
    inline const std::vector<S3VertexDesc>& ends() const { return this->end_verts; }
    inline const std::vector<Plane3>& planes() const { return this->end_planes; }
    inline const std::vector<std::vector<Plane3>>& bp_aux_planes() const { return this->aux_planes; }
    inline const std::vector<Plane3>& all_planes() const { return this->_all_planes; }; // planes() plus planes for neighboring branch points if they make sense

    // Access to the meshes which are returned as mutable objects
    inline Polyhedron3* mesh() { return this->_mesh; }
    inline const Polyhedron3* mesh() const { return this->_mesh; }
    inline Polyhedron3* uncapped_mesh() { return this->uncapped; }
    inline const Polyhedron3* uncapped_mesh() const { return this->uncapped; }
};

typedef std::vector<Slice*> Slices;

Slices slice(const int group_sz, const int bp_group_sz, const Polyhedron3* P, const Skeleton3* S);
