#pragma once

#include <assert.h>
#include <vector>
#include "Iterators.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// A skeleton is simply a graph with vertices and edges in which every vertex has a 3D point.
// This class represents the skeleton such that every edge between any two spatial coordinates is
// explicit. In general it is easier to work with a SkeletonGraph_3 (which can be constructed from
// a Skeleton_3).
///////////////////////////////////////////////////////////////////////////////////////////////////

//#define DETAILED_SKELETON // Defining this makes it so the skeleton data structure has a lot of detailed information on each edge


///////////////////////////////////////////////////////////////////////////////////////////////////
// The vertices in a Skeleton_3, which possibly has many edges incident to it.
// The standard version only stores the point of the vertex and the edges incident to it.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename S3>
class Skeleton_3_Vertex
{
	friend S3;
public:
	typedef typename S3::Point  Point;
	typedef typename S3::Vertex Vertex;
	typedef typename S3::Edge   Edge;
	typedef typename S3::Vertex_handle Vertex_handle;
	typedef typename S3::Edge_handle   Edge_handle;
	typedef typename S3::Vertex_const_handle Vertex_const_handle;
	typedef typename S3::Edge_const_handle   Edge_const_handle;
	typedef std::vector<Edge_handle> Edge_vector;
	typedef handle_iterator<typename Edge_vector::iterator,       Edge_handle>       Edge_iterator;
	typedef handle_iterator<typename Edge_vector::const_iterator, Edge_const_handle> Edge_const_iterator;
private:
	Point p;
	Edge_vector es;
public:
	inline Skeleton_3_Vertex(const Point& p) : p(p), es() { }
	inline const Point& point() const { return this->p; }
	inline size_t degree() const { return this->es.size(); } // number of edges that are incident
	inline         Edge_iterator edges_begin()       { return this->es.begin(); }
	inline         Edge_iterator edges_end()         { return this->es.end();   }
	inline   Edge_const_iterator edges_begin() const { return this->es.begin(); }
	inline   Edge_const_iterator edges_end()   const { return this->es.end();   }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// The edges in a Skeleton_3. The standard version only stores the two end vertices.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename S3>
class Skeleton_3_Edge
{
	friend S3;
public:
	typedef typename S3::Vertex Vertex;
	typedef typename S3::Edge   Edge;
	typedef typename S3::Vertex_handle Vertex_handle;
	typedef typename S3::Edge_handle   Edge_handle;
	typedef typename S3::Vertex_const_handle Vertex_const_handle;
	typedef typename S3::Edge_const_handle   Edge_const_handle;
private:
	Vertex_handle u, v;
public:
	inline Skeleton_3_Edge() : u(), v() { }
	inline Vertex_handle       source()       { return this->u; }
	inline Vertex_const_handle source() const { return this->u; }
	inline Vertex_handle       target()       { return this->v; }
	inline Vertex_const_handle target() const { return this->v; }
	inline size_t index(Vertex_const_handle v) const { if (v == this->u) { return 0; } else { assert(v == this->v); return 1; } }
	inline Vertex_handle       opposite(Vertex_const_handle v)       { if (v == this->u) { return this->v; } else { assert(v == this->v); return this->u; } }
	inline Vertex_const_handle opposite(Vertex_const_handle v) const { if (v == this->u) { return this->v; } else { assert(v == this->v); return this->u; } }
	inline Vertex_handle       vertex(size_t i)       { if (i == 0) { return this->u; } else { assert(i == 1); return this->v; } }
	inline Vertex_const_handle vertex(size_t i) const { if (i == 0) { return this->u; } else { assert(i == 1); return this->v; } }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Represents a 3D 'skeleton' - essentially a graph with 3D points associated with each vertex
// where every edge and vertex is explicitly represented.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, template<typename> class V = Skeleton_3_Vertex, template<typename> class E = Skeleton_3_Edge>
class Skeleton_3
{
public:
	typedef Skeleton_3<K, V, E> Self;
	typedef typename K::Point_3 Point;
	typedef K       Kernel;
	typedef V<Self> Vertex;
	typedef E<Self> Edge;
	
	typedef std::vector<Vertex> Vertex_vector;
	typedef std::vector<Edge>   Edge_vector;

	typedef stable_vector_iterator<      Vertex_vector, typename Vertex_vector::iterator>       Vertex_iterator;
	typedef stable_vector_iterator<const Vertex_vector, typename Vertex_vector::const_iterator> Vertex_const_iterator;
	typedef stable_vector_iterator<      Edge_vector,   typename Edge_vector::iterator>         Edge_iterator;
	typedef stable_vector_iterator<const Edge_vector,   typename Edge_vector::const_iterator>   Edge_const_iterator;

	typedef Vertex_iterator       Vertex_handle;
	typedef Vertex_const_iterator Vertex_const_handle;
	typedef Edge_iterator         Edge_handle;
	typedef Edge_const_iterator   Edge_const_handle;

private:
	Vertex_vector vs;
	Edge_vector   es;

public:
	inline explicit Skeleton_3(size_t v = 0, size_t e = 0) : vs(), es() { this->vs.reserve(v); this->es.reserve(e); }

	inline size_t size_of_vertices()     const { return this->vs.size();     }
	inline size_t size_of_edges()        const { return this->es.size();     }
	inline size_t capacity_of_vertices() const { return this->vs.capacity(); }
	inline size_t capacity_of_edges()    const { return this->es.capacity(); }
	inline bool empty() const { return this->vs.empty(); }
	inline void shrink_to_fit()
	{
		this->vs.shrink_to_fit();
		this->es.shrink_to_fit();
		for (typename Vertex_vector::iterator v = this->vs.begin(), end = this->vs.begin(); v != end; ++v) { v->es.shrink_to_fit(); }
	}
	
	inline Vertex_iterator vertices_begin() { return Vertex_iterator(&this->vs, 0); }
	inline Vertex_iterator vertices_end()   { return Vertex_iterator(&this->vs, this->vs.size()); }
	inline Edge_iterator   edges_begin()    { return Edge_iterator  (&this->es, 0); }
	inline Edge_iterator   edges_end()      { return Edge_iterator  (&this->es, this->es.size()); }
	inline Vertex_const_iterator vertices_begin() const { return Vertex_const_iterator(&this->vs, 0); }
	inline Vertex_const_iterator vertices_end()   const { return Vertex_const_iterator(&this->vs, this->vs.size()); }
	inline Edge_const_iterator   edges_begin()    const { return Edge_const_iterator  (&this->es, 0); }
	inline Edge_const_iterator   edges_end()      const { return Edge_const_iterator  (&this->es, this->es.size()); }

	// Add an isolated vertex at the given point to the collection and return a handle to it.
	inline Vertex_handle add_vertex(const Point& p) { this->vs.push_back(Vertex(p)); return Vertex_handle(&this->vs, this->vs.size() - 1); }
	// Add an isolated vertex to the collection and return a handle to it.
	inline Vertex_handle add_vertex(const Vertex& v) { this->vs.push_back(v); return Vertex_handle(&this->vs, this->vs.size() - 1); }
	// Add an edge to the collection. The edge is defined by the two vertices given. Both vertices
	// are updated to know about the new incident edge. Return a handle to the new edge.
	inline Edge_handle add_edge(Vertex_handle u, Vertex_handle v, const Edge& e = Edge())
	{
		this->es.push_back(e);
		Edge_handle eh = Edge_handle(&this->es, this->es.size() - 1);
		eh->u = u; u->es.push_back(eh);
		eh->v = v; v->es.push_back(eh);
		//this->check();
		return eh;
	}
	//inline void check() const
	//{
	//	for (Vertex_const_iterator v = this->vertices_begin(), end = this->vertices_end(); v != end; ++v)
	//	{
	//		for (Vertex::Edge_const_iterator e = v->edges_begin(), end = v->edges_end(); e != end; ++e)
	//		{
	//			assert((e->source() == v) != (e->target() == v)); // make sure one of the endpoints of each incident edge is the vertex
	//		}
	//	}
	//	for (Edge_const_iterator e = this->edges_begin(), e_end = this->edges_end(); e != e_end; ++e)
	//	{
	//		assert(e->u != Vertex_handle() && e->v != Vertex_handle());
	//		Vertex::Edge_const_iterator i, end;
	//		for (i = e->u->edges_begin(), end = e->u->edges_end(); i != end && (Edge_const_handle)i != (Edge_const_handle)e; ++i); assert(i != end);
	//		for (i = e->v->edges_begin(), end = e->v->edges_end(); i != end && (Edge_const_handle)i != (Edge_const_handle)e; ++i); assert(i != end);
	//	}
	//}
};


#ifdef DETAILED_SKELETON
///////////////////////////////////////////////////////////////////////////////////////
// The edge of a detailed skeleton structure, includes eccentricity, geodesic distance
// angle, stability, and touching vertices.
///////////////////////////////////////////////////////////////////////////////////////
template<typename S3>
class Detailed_Skeleton_3_Edge : public Skeleton_3_Edge<S3>
{
	friend S3;
	double _eccentricity;
	double _geod;
	double _angle;
	bool _stable;
	unsigned int _n_touching_verts;
	Polyhedron3::Vertex_const_handle _touching_verts[3];
	template<class set> static inline bool toggle(set s, typename set::key_type k) { set::iterator i = s.find(k); if (i != s.end()) { s.erase(i); return false; } else { s.insert(k); return true; } } // returns true if inserted, false if removed
public:
	Detailed_Skeleton_3_Edge() : Skeleton_3_Edge<S3>(), _stable(false), _n_touching_verts(0), _eccentricity(0.0), _geod(0.0), _angle(0.0)
	{
	}
	Detailed_Skeleton_3_Edge(MAT::Facet_handle f, MAT::Edge_handle e, double wt_thd) : Skeleton_3_Edge<S3>(), _stable(e->flags & MAT::Edge::Flags::StableSkeleton()), _n_touching_verts(2)
	{
		typedef handle_set<Polyhedron3::Vertex_const_handle>::type p3vertices;

		double geod = f->weight(), angle = e->flux(), circumference;
		Polyhedron3::Vertex_const_handle tvs = f->touching_vertex_start(), tve = f->touching_vertex_end(), tvx;
		this->_touching_verts[0] = tvs;
		this->_touching_verts[1] = tve;
		if (e->is_unique())
		{
			// Compute the eccentricity for the edge if there is no adjacent facet along the edge
			// Means two vertices on the mesh forming the Delaunay triangle for the skeleton edge
			geod *= 2;
			circumference = CGAL_PI * distance(tvs->point(), tve->point());
		}
		else
		{
			p3vertices vt; // max size is 3
			vt.insert(tvs);
			vt.insert(tve);
			MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
			for (; ei != e; ++ei)
			{
				ei->flags |= MAT::Edge::Flags::Visited();
				MAT::Facet_iterator fa = ei->facet();
				geod += fa->weight(); // ? should have three rings here
				if (fa->weight() > wt_thd && ei->flux() < angle) { angle = ei->flux(); }
				Polyhedron3::Vertex_const_handle vs = fa->touching_vertex_start(), ve = fa->touching_vertex_end();
				if (toggle(vt, vs)) { assert(tvx == Polyhedron3::Vertex_handle()); tvx = vs; }
				if (toggle(vt, ve)) { assert(tvx == Polyhedron3::Vertex_handle()); tvx = ve; }
			}
			this->_touching_verts[2] = tvx;
			++this->_n_touching_verts;
			if (vt.size() != 0)
			{
				assert(vt.size() == 2);
				p3vertices::const_iterator vit = vt.begin();
				geod += distance((*vit)->point(), (*++vit)->point());
			}
			circumference = 2 * CGAL_PI * distance(CGAL::circumcenter(tvs->point(), tve->point(), tvx->point()), tvs->point());
		}
		this->_eccentricity = geod / circumference;
		this->_geod = geod;
		this->_angle = (angle > 0 ? angle : 0);
	}
	inline bool   has_details()  const { return this->_n_touching_verts >= 2; }
	inline bool   stable()       const { return this->_stable; }
	inline double eccentricity() const { return this->_eccentricity; }
	inline double geod()         const { return this->_geod; }
	inline double angle()        const { return this->_angle; }
	inline unsigned int size_of_touching_vertices() const { return this->_n_touching_verts; }
	inline Polyhedron3::Vertex_handle touching_vertex(unsigned int i) const { assert(i < this->_n_touching_verts); return this->_touching_verts[i]; }
};
#define SKELETON_3(Kernel) Skeleton_3<Kernel, Skeleton_3_Vertex, Detailed_Skeleton_3_Edge>
#else
#define SKELETON_3(Kernel) Skeleton_3<Kernel>
#endif
