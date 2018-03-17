#pragma once

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are utilities for Polyhedron3 objects.
///////////////////////////////////////////////////////////////////////////////////////////////////

// Polyhedron part converters
Polygon2 facet_to_polygon2(Polyhedron3::Facet_const_handle f);
inline Plane3 facet_to_plane3(Polyhedron3::Facet_const_handle f) { const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next(); return Plane3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()); }
inline Triangle3 facet_to_triangle3(Polyhedron3::Facet_const_handle f) { const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next(); return Triangle3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()); }
inline Bbox3 facet_to_bbox3(Polyhedron3::Facet_const_handle f) { const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next(); return bbox3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point()); }
inline Segment3 halfedge_to_segment3(Polyhedron3::Halfedge_const_handle he) { return Segment3(he->prev()->vertex()->point(), he->vertex()->point()); }


// Polyhedron part normals
inline Direction3 normal(const Polyhedron3::Facet_const_handle &f)
{
    // assumes triangle face
    Polyhedron3::Halfedge_around_facet_const_circulator he = f->facet_begin();
    Vector3 n = CGAL::normal(he->prev()->vertex()->point(), he->vertex()->point(), he->next()->vertex()->point());
    return Direction3(normalized(n));
}
inline Direction3 normal(const Polyhedron3::Vertex_const_handle &v)
{
    Vector3 n = CGAL::NULL_VECTOR;
    Polyhedron3::Halfedge_around_vertex_const_circulator he = v->vertex_begin(), end(he);
    CGAL_For_all(he, end) { if (!he->is_border()) { n = n + to_vector(normal(he->facet())); } }
    return Direction3(normalized(n));
}

inline void fix_ids(Polyhedron3* P)
{
    // Make sure all of the facets and vertices have correct ids
    // This is needed for several PMP functions
    std::size_t i = 0;
    for (auto itr = P->facets_begin(), end = P->facets_end(); itr != end; ++itr, ++i) { itr->id() = i; }
    i = 0;
    for (auto itr = P->vertices_begin(), end = P->vertices_end(); itr != end; ++itr, ++i) { itr->id() = i; }
}


// Polyhedron looping around a vertex/facet
#define FOR_EDGES_AROUND_VERTEX(V, e) \
    auto MAKE_UNIQUE(_he) = V->vertex_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto e = MAKE_UNIQUE(_he); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), e = MAKE_UNIQUE(_he))
#define FOR_FACETS_AROUND_VERTEX(V, f) \
    auto MAKE_UNIQUE(_he) = V->vertex_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto f = MAKE_UNIQUE(_he)->facet(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), f = MAKE_UNIQUE(_he)->facet())
#define FOR_VERTICES_AROUND_VERTEX(V, v) \
    auto MAKE_UNIQUE(_he) = V->vertex_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto v = MAKE_UNIQUE(_he)->prev()->vertex(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), v = MAKE_UNIQUE(_he)->prev()->vertex())
#define FOR_EDGES_AROUND_FACET(F, e) \
    auto MAKE_UNIQUE(_he) = F->facet_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto e = MAKE_UNIQUE(_he); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), e = MAKE_UNIQUE(_he))
#define FOR_FACETS_AROUND_FACET(F, f) \
    auto MAKE_UNIQUE(_he) = F->facet_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto f = MAKE_UNIQUE(_he)->opposite()->facet(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), f = MAKE_UNIQUE(_he)->opposite()->facet())
#define FOR_VERTICES_AROUND_FACET(F, v) \
    auto MAKE_UNIQUE(_he) = F->facet_begin(), MAKE_UNIQUE(_end) = MAKE_UNIQUE(_he); \
    auto v = MAKE_UNIQUE(_he)->vertex(); \
    for (bool _circ_loop_flag = ! ::CGAL::is_empty_range(MAKE_UNIQUE(_he), MAKE_UNIQUE(_end)); _circ_loop_flag; _circ_loop_flag = (++MAKE_UNIQUE(_he) != MAKE_UNIQUE(_end)), v = MAKE_UNIQUE(_he)->vertex())


// Other utilities
inline Kernel::FT avg_edge_length(const Polyhedron3 *P)
{
    Kernel::FT length = 0;
    for (auto e = P->edges_begin(), end = P->edges_end(); e != end; ++e)
    {
        length += distance(e->vertex()->point(), e->opposite()->vertex()->point());
    }
    return length / Kernel::FT(P->size_of_halfedges() / 2);
}
Kernel::FT volume(const Polyhedron3* P);
Kernel::FT surface_area(const Polyhedron3* P);
bool is_not_degenerate(const Polyhedron3* P);


// Mesh Extracting
template <class Polyhedron_3, class VertexCollection>
class MeshExtractorV : public CGAL::Modifier_base<typename Polyhedron3::HalfedgeDS>
{
    typedef typename CGAL::Polyhedron_incremental_builder_3<typename Polyhedron_3::HalfedgeDS> Builder;
    typedef handle_map<typename Polyhedron_3::Vertex_const_handle, size_t> Vertex_int_map;
    typedef handle_set<typename Polyhedron_3::Facet_const_handle> Facet_set;
    typedef typename Polyhedron_3::Vertex_const_handle Vertex;
    const VertexCollection verts;
    const bool only_enclosed;
public:
    MeshExtractorV(const VertexCollection verts, bool only_enclosed=true) :
        verts(verts), only_enclosed(only_enclosed) {}
    void operator()(typename Polyhedron_3::HalfedgeDS& hds)
    {
        // Start building the polyhedron
        Builder B(hds, true);
        B.begin_surface(this->verts.size(), 2*this->verts.size()); // number of faces is very approximate

        // Add all of the vertices
        Vertex_int_map lookup;
        BOOST_FOREACH(Vertex v, this->verts)
        {
            B.add_vertex(v->point());
            lookup.insert({{v, lookup.size()}});
        }

        // Add faces (all faces that have all vertices in the lookup table)
        Facet_set faces; // only add each face once
        std::vector<size_t> v_inds;
        BOOST_FOREACH(Vertex v, this->verts)
        {
            FOR_FACETS_AROUND_VERTEX(v, f)
            {
                if (faces.count(f)) { continue; } // face already processed
                faces.insert(f); // now mark it as processed
                FOR_VERTICES_AROUND_FACET(f, vf)
                {
                    auto itr = lookup.find(vf);
                    if (itr != lookup.end()) { v_inds.push_back(itr->second); }
                    else if (this->only_enclosed) { v_inds.clear(); break; } // do not add face - not all vertices are present
                    else
                    {
                        // Add vertex to polyhedron and then proceed with adding it to the facet
                        B.add_vertex(vf->point());
                        size_t idx = lookup.size();
                        lookup.insert({{vf, idx}});
                        v_inds.push_back(idx);
                    }
                }
                if (v_inds.size() != 0)
                {
                    B.add_facet(v_inds.begin(), v_inds.end());
                    v_inds.clear();
                }
            }
        }
        B.end_surface();
    }
};

template <class VertexCollection>
Polyhedron3* extract_mesh_from_vertices(VertexCollection verts, bool only_enclosed=true)
{
    MeshExtractorV<Polyhedron3, VertexCollection> extractor(verts, only_enclosed);
    Polyhedron3* mesh = new Polyhedron3();
    mesh->delegate(extractor);
    return mesh;
}
