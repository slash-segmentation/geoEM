#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

#include <utility>
#include <vector>
#include <unordered_set>

class Slice
{
    // Note: the ends of a slice must only have one outgoing edge each
    const Skeleton3* S;
    std::unordered_set<S3VertexDesc> svs;
    size_t deg = 0;
    // All of the end_* are in the same order
    std::vector<S3VertexDesc> end_verts;
    std::vector<Plane3> end_planes;
    std::vector<Slice*> end_neighbors;
    Polyhedron3* _mesh = nullptr;
    
public:
    template <class ForwardIterator>
    inline Slice(const Skeleton3* S, ForwardIterator begin, ForwardIterator end, Polyhedron3* mesh=nullptr)
        : S(S), svs(begin, end), _mesh(mesh == nullptr ? new Polyhedron3() : mesh)
    {
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
    }
    
    // Copy constructors copy the mesh unless given a new mesh
    // Neighbors are not copied, must call set_neighbors
    inline Slice(const Slice& slc, Polyhedron3* mesh=nullptr)
        : S(slc.S), svs(slc.svs), deg(slc.deg), 
        end_verts(slc.end_verts), end_planes(slc.end_planes), end_neighbors(slc.end_neighbors.size(), nullptr),
        _mesh(mesh == nullptr ? new Polyhedron3(*slc._mesh) : mesh) { }
    inline Slice(const Slice* slc, Polyhedron3* mesh=nullptr)
        : S(slc->S), svs(slc->svs), deg(slc->deg),
        end_verts(slc->end_verts), end_planes(slc->end_planes), end_neighbors(slc->end_neighbors.size(), nullptr),
        _mesh(mesh == nullptr ? new Polyhedron3(*slc->_mesh) : mesh) { }
    Slice& operator=(const Slice&) = delete;
    inline ~Slice() { delete _mesh; }
    
    inline size_t degree() const { return this->deg; }
    inline std::unordered_set<S3VertexDesc>& skeleton_vertices() { return this->svs; }
    inline const std::unordered_set<S3VertexDesc>& skeleton_vertices() const { return this->svs; }
    
    inline std::vector<S3VertexDesc>& ends() { return this->end_verts; }
    inline const std::vector<S3VertexDesc>& ends() const { return this->end_verts; }

    inline std::vector<Plane3>& planes() { return this->end_planes; }
    inline const std::vector<Plane3>& planes() const { return this->end_planes; }

    inline std::vector<Slice*>& neighbors() { return this->end_neighbors; }
    inline const std::vector<Slice*>& neighbors() const { return this->end_neighbors; }

    inline Polyhedron3* mesh() { return this->_mesh; }
    inline const Polyhedron3* mesh() const { return this->_mesh; }
    
    /* TODO: enable
    inline Slice* capped() const
    {
        Slice* slc = new Slice(this);
        for (auto he = slc->_mesh->halfedges_begin(); he != slc->_mesh->halfedges_end(); ++he)
        {
            if (he->is_border())
            {
                std::vector<boost::graph_traits<Polyhedron3>::face_descriptor> x;
                CGAL::Polygon_mesh_processing::triangulate_hole(*slc->_mesh, he, std::back_inserter(x));
            }
        }
        return slc;
    }*/
    
    // TODO: add distance measurement (total distance between skeletal points (including midpoints)
    
    inline static void set_neighbors(std::vector<Slice*> slices)
    {
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
};

typedef std::vector<Slice*> Slices;

Slices slice(const int group_sz, const Polyhedron3* P, const Skeleton3* S,
    const P3CVertexSet& facet_verts, bool verbose=false);
