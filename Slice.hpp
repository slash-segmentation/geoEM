#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include <utility>
#include <vector>
#include <unordered_set>

class Slice
{
private:
    // Note: the ends of a slice must only have one outgoing edge each
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
    Polyhedron3* _mesh;
    Polyhedron3* uncapped = nullptr;

    void init();

public:
    template <class ForwardIterator>
    inline Slice(const Skeleton3* S, ForwardIterator begin, ForwardIterator end, Polyhedron3* mesh=nullptr)
        : S(S), svs(begin, end), _mesh(mesh == nullptr ? new Polyhedron3() : mesh) { init(); }

    // Copy constructors copy the mesh unless given a new mesh
    // Neighbors are not copied, must call set_neighbors
    // TODO: copy uncapped?
    inline Slice(const Slice& slc, Polyhedron3* mesh=nullptr)
        : S(slc.S), svs(slc.svs), deg(slc.deg),
        end_verts(slc.end_verts), end_planes(slc.end_planes), end_neighbors(slc.end_neighbors.size(), nullptr), aux_planes(slc.aux_planes),
        _mesh(mesh == nullptr ? new Polyhedron3(*slc._mesh) : mesh) { }
    inline Slice(const Slice* slc, Polyhedron3* mesh=nullptr)
        : S(slc->S), svs(slc->svs), deg(slc->deg),
        end_verts(slc->end_verts), end_planes(slc->end_planes), end_neighbors(slc->end_neighbors.size(), nullptr), aux_planes(slc->aux_planes),
        _mesh(mesh == nullptr ? new Polyhedron3(*slc->_mesh) : mesh) { }
    Slice& operator=(const Slice&) = delete;
    inline ~Slice() { delete _mesh; if (uncapped != nullptr) { delete uncapped; } }

    void build_mesh(const Polyhedron3* P);
    static void set_neighbors(std::vector<Slice*> slices);

    Kernel::FT length() const;
    inline size_t degree() const { return this->deg; }
    inline const Skeleton3* skeleton() const { return this->S; }

    inline std::unordered_set<S3VertexDesc>& skeleton_vertices() { return this->svs; }
    inline const std::unordered_set<S3VertexDesc>& skeleton_vertices() const { return this->svs; }

    inline std::vector<S3VertexDesc>& ends() { return this->end_verts; }
    inline const std::vector<S3VertexDesc>& ends() const { return this->end_verts; }

    inline std::vector<Plane3>& planes() { return this->end_planes; }
    inline const std::vector<Plane3>& planes() const { return this->end_planes; }

    inline std::vector<std::vector<Plane3>>& bp_aux_planes() { return this->aux_planes; }
    inline const std::vector<std::vector<Plane3>>& bp_aux_planes() const { return this->aux_planes; }

    std::vector<Plane3> all_planes() const; // planes for this slice and some of the aux planes from neighboring branch points

    inline std::vector<Slice*>& neighbors() { return this->end_neighbors; }
    inline const std::vector<Slice*>& neighbors() const { return this->end_neighbors; }

    inline Polyhedron3* mesh() { return this->_mesh; }
    inline const Polyhedron3* mesh() const { return this->_mesh; }

    inline Polyhedron3* uncapped_mesh() { return this->uncapped; }
    inline const Polyhedron3* uncapped_mesh() const { return this->uncapped; }
};

typedef std::vector<Slice*> Slices;

Slices slice(const int group_sz, const int bp_group_sz, const Polyhedron3* P, const Skeleton3* S);
