#pragma once

#include <CGAL/circulator_bases.h>
#include <vector>
#include "Iterators.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Similar to Polyhedron_3 except it is really just a collection of planar-polygons that may share
// vertices/edges. The reason why Polyhedron_3 and Linear_cell_complex<3> cannot be used is because
// they require each edge to have 1 or 2 incident facets or facets cannot be incident at a single
// vertex respectively which we cannot assume. By not having those restrictions though things are
// made a little more complicated/slower. The general queries we support are:
//   All edges that end on a vertex (e->next_around_facet()->vertex() == v), indirectly provides:
//     All facets incident to vertex (e->facet())
//     All edges incident to vertex (e && e->next_around_facet())
//   The vertices of an edge (using e->vertex() and e->next_around_facet()->vertex())
//   The facets incident to an edge
//   The edges around a facet
//   The vertices around a facet
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// First we define a circulator we need for Polygons_3.
///////////////////////////////////////////////////////////////////////////////////////////////////
// These unary functors take an edge handle and get the next edge in a circulator (either around facet or identified)
template <typename EH> struct Polygons_3_Next_edge_around_facet { inline static EH next(EH e) { return e->next_around_facet(); } };
template <typename EH> struct Polygons_3_Next_identified_edge   { inline static EH next(EH e) { return e->next_identifed();    } };

// Circulator for edges either around a facet or an identified edge depending on _Next (one of the above unary functors above)
template <typename EH, template<typename> class _Next = Polygons_3_Next_edge_around_facet>
class Polygons_3_Edge_circulator : public CGAL::Circulator_base<CGAL::Forward_circulator_tag, typename EH::value_type, ptrdiff_t, size_t, EH>
{
	template <typename EH2, template<typename> class _Next2> friend class Polygons_3_Edge_circulator;
public:
	typedef EH Edge_handle;
	typedef _Next<Edge_handle> Next;
	typedef Polygons_3_Edge_circulator<EH, _Next> type;
	typedef CGAL::Circulator_base<CGAL::Forward_circulator_tag, typename EH::value_type, ptrdiff_t, size_t, EH> base;
private:
	Edge_handle e;
public:
	inline Polygons_3_Edge_circulator() : e() { } // an empty circulator (C() == nullptr is true)
	inline Polygons_3_Edge_circulator(type const& rhs) : e(rhs.e) { }
	inline type& operator=(const type& rhs) { this->e = rhs.e; return *this; }
	inline explicit Polygons_3_Edge_circulator(Edge_handle e) : e(e) { }
	template <typename EH2> inline Polygons_3_Edge_circulator(const Polygons_3_Edge_circulator<EH2,_Next>& i) : e(i.e) { }
	inline bool operator==(CGAL::Nullptr_t CGAL_assertion_code(rhs)) const { CGAL_assertion(rhs == 0); return this->e == Edge_handle(); }
	inline bool operator!=(CGAL::Nullptr_t CGAL_assertion_code(rhs)) const { CGAL_assertion(rhs == 0); return this->e != Edge_handle(); }
	inline bool operator==(const type& rhs) const { return this->e == rhs.e; }
	inline bool operator!=(const type& rhs) const { return this->e != rhs.e; }
	inline bool operator==(const Edge_handle& rhs) const { return this->e == rhs; }
	inline bool operator!=(const Edge_handle& rhs) const { return this->e != rhs; }
	//friend inline bool operator==(const Edge_handle& lhs, const type& rhs) { return lhs == rhs.e; } // these are covered by casts
	//friend inline bool operator!=(const Edge_handle& lhs, const type& rhs) { return lhs != rhs.e; }
	inline operator typename base::reference() const { return *this->e; }
	inline operator typename base::pointer() const { return this->e; }
	inline typename base::reference operator*() const { return *this->e; }
	inline typename base::pointer operator->() const { return this->e; }
	inline type& operator++() { this->e = Next::next(e); return *this; }
	inline type operator++(int) { type x = *this; this->e = Next::next(e); return x; }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// The vertices in a Polygons_3, which possibly has many facets incident to it.
// The standard version only stores the point of the vertex and the edges that terminate on it.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename P3>
class Polygons_3_Vertex
{
	friend P3;
public:
	typedef typename P3::Point  Point;
	typedef typename P3::Vertex Vertex;
	typedef typename P3::Edge   Edge;
	typedef typename P3::Facet  Facet;
	typedef typename P3::Vertex_handle Vertex_handle;
	typedef typename P3::Edge_handle   Edge_handle;
	typedef typename P3::Facet_handle  Facet_handle;
	typedef typename P3::Vertex_const_handle Vertex_const_handle;
	typedef typename P3::Edge_const_handle   Edge_const_handle;
	typedef typename P3::Facet_const_handle  Facet_const_handle;
	typedef std::vector<Edge_handle> Edge_vector;
	typedef handle_iterator<typename Edge_vector::iterator,       Edge_handle>       Edge_iterator;
	typedef handle_iterator<typename Edge_vector::const_iterator, Edge_const_handle> Edge_const_iterator;
private:
	Point p;
	Edge_vector es;
public:
	inline Polygons_3_Vertex(const Point& p) : p(p), es() { }
	inline const Point& point() const { return this->p; }
	inline size_t degree() const { return this->es.size(); } // number of facets that are incident
	inline         Edge_iterator edges_begin()          { return this->es.begin(); }
	inline         Edge_iterator edges_end()            { return this->es.end();   }
	inline   Edge_const_iterator edges_begin()    const { return this->es.begin(); }
	inline   Edge_const_iterator edges_end()      const { return this->es.end();   }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// The (semi-)edges in a Polygons_3. These edges are unique to a facet (similar to a halfedge in
// Polyhedron_3). In a Polygons_3 there may be many "identified" edges (spatially identical but
// with different incident facets). The edge goes from e->vertex() to e->next_around_facet()->vertex().
// To get the other identified edges use next_identified()/circ_identified().
// The standard version stores the starting vertex, the incident facet, the next edges around the
// facet, and the next identified edge.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename P3>
class Polygons_3_Edge
{
	friend P3;
public:
	typedef typename P3::Vertex Vertex;
	typedef typename P3::Edge   Edge;
	typedef typename P3::Facet  Facet;
	typedef typename P3::Vertex_handle Vertex_handle;
	typedef typename P3::Edge_handle   Edge_handle;
	typedef typename P3::Facet_handle  Facet_handle;
	typedef typename P3::Vertex_const_handle Vertex_const_handle;
	typedef typename P3::Edge_const_handle   Edge_const_handle;
	typedef typename P3::Facet_const_handle  Facet_const_handle;
	typedef typename P3::Edge_around_facet_circulator       Edge_around_facet_circulator;
	typedef typename P3::Edge_around_facet_const_circulator Edge_around_facet_const_circulator;
	typedef typename P3::Edge_identified_circulator         Edge_identified_circulator;
	typedef typename P3::Edge_identified_const_circulator   Edge_identified_const_circulator;
private:
	Vertex_handle v;
	Facet_handle f;
	Edge_handle ef, ei;
protected:
	inline Polygons_3_Edge(Vertex_handle v, Facet_handle f) : v(v), f(f) { }
public:
	inline size_t degree() const { return CGAL::circulator_size(this->circ_identified()); } // number of identified edges
	inline bool is_unique() const { return this->ei == this->ei; } // has no other identified edges (degree() == 0)
	inline Vertex_handle       vertex()                  { return this->v; }
	inline Vertex_const_handle vertex()            const { return this->v; }
	inline Facet_handle        facet()                   { return this->f; }
	inline Facet_const_handle  facet()             const { return this->f; }
	inline Edge_handle         next_around_facet()       { return this->ef; }
	inline Edge_const_handle   next_around_facet() const { return this->ef; }
	inline Edge_handle         next_identifed()          { return this->ei; }
	inline Edge_const_handle   next_identifed()    const { return this->ei; }
	inline Edge_around_facet_circulator       circ_around_facet()       { return Edge_around_facet_circulator      (this->ef); }
	inline Edge_around_facet_const_circulator circ_around_facet() const { return Edge_around_facet_const_circulator(this->ef); }
	inline Edge_identified_circulator         circ_identified()         { return Edge_identified_circulator        (this->ei); }
	inline Edge_identified_const_circulator   circ_identified()   const { return Edge_identified_const_circulator  (this->ei); }
	inline bool is_identified_with(const Edge& e) const        { return (this->vertex() == e.v  && this->ef->v == e.ef->v ) || (this->v == e.ef->v  && this->ef->v == e.v);  }
	inline bool is_identified_with(const Edge_handle& e) const { return (this->vertex() == e->v && this->ef->v == e->ef->v) || (this->v == e->ef->v && this->ef->v == e->v); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// The facets in a Polygons_3. The standard version only stores a single incident edge.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename P3>
class Polygons_3_Facet
{
	friend P3;
public:
	typedef typename P3::Vertex Vertex;
	typedef typename P3::Edge   Edge;
	typedef typename P3::Facet  Facet;
	typedef typename P3::Vertex_handle Vertex_handle;
	typedef typename P3::Edge_handle   Edge_handle;
	typedef typename P3::Facet_handle  Facet_handle;
	typedef typename P3::Vertex_const_handle Vertex_const_handle;
	typedef typename P3::Edge_const_handle   Edge_const_handle;
	typedef typename P3::Facet_const_handle  Facet_const_handle;
	typedef typename P3::Edge_around_facet_circulator         Edge_around_facet_circulator;
	typedef typename P3::Edge_around_facet_const_circulator   Edge_around_facet_const_circulator;
private:
	Edge_handle e;
protected:
	virtual void finish() { /* no-op, can be implemented by derived classes, called once the facet is done being created */ }
public:
	inline Polygons_3_Facet() : e() { }
	inline size_t degree() const { return CGAL::circulator_size(this->edges_circ()); } // number of edges/vertices around facet
	inline   Edge_around_facet_circulator       edges_circ()          { return Edge_around_facet_circulator        (this->e); }
	inline   Edge_around_facet_const_circulator edges_circ()    const { return Edge_around_facet_const_circulator  (this->e); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Represents a collection of planar-polygons in 3D space. Contains all the vertices, edges, and
// facets of this collection.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, template<typename> class V = Polygons_3_Vertex, template<typename> class E = Polygons_3_Edge, template<typename> class F = Polygons_3_Facet>
class Polygons_3
{
public:
	typedef Polygons_3<K, V, E, F> Self;
	typedef typename K::Point_3 Point;
	typedef K       Kernel;
	typedef V<Self> Vertex;
	typedef E<Self> Edge;
	typedef F<Self> Facet;
	
	typedef std::vector<Vertex> Vertex_vector;
	typedef std::vector<Edge>   Edge_vector;
	typedef std::vector<Facet>  Facet_vector;

	typedef stable_vector_iterator<      Vertex_vector, typename Vertex_vector::iterator>       Vertex_iterator;
	typedef stable_vector_iterator<const Vertex_vector, typename Vertex_vector::const_iterator> Vertex_const_iterator;
	typedef stable_vector_iterator<      Edge_vector,   typename Edge_vector::iterator>         Edge_iterator;
	typedef stable_vector_iterator<const Edge_vector,   typename Edge_vector::const_iterator>   Edge_const_iterator;
	typedef stable_vector_iterator<      Facet_vector,  typename Facet_vector::iterator>        Facet_iterator;
	typedef stable_vector_iterator<const Facet_vector,  typename Facet_vector::const_iterator>  Facet_const_iterator;

	typedef Vertex_iterator       Vertex_handle;
	typedef Vertex_const_iterator Vertex_const_handle;
	typedef Edge_iterator         Edge_handle;
	typedef Edge_const_iterator   Edge_const_handle;
	typedef Facet_iterator        Facet_handle;
	typedef Facet_const_iterator  Facet_const_handle;

	typedef Polygons_3_Edge_circulator<Edge_handle      > Edge_around_facet_circulator;
	typedef Polygons_3_Edge_circulator<Edge_const_handle> Edge_around_facet_const_circulator;
	typedef Polygons_3_Edge_circulator<Edge_handle,       Polygons_3_Next_identified_edge> Edge_identified_circulator;
	typedef Polygons_3_Edge_circulator<Edge_const_handle, Polygons_3_Next_identified_edge> Edge_identified_const_circulator;

private:
	Vertex_vector vs;
	Edge_vector   es;
	Facet_vector  fs;

public:
	inline explicit Polygons_3(size_t v = 0, size_t e = 0, size_t f = 0) : vs(), es(), fs() { this->vs.reserve(v); this->es.reserve(e); this->fs.reserve(f); }

	inline size_t size_of_vertices()     const { return this->vs.size();     }
	inline size_t size_of_edges()        const { return this->es.size();     }
	inline size_t size_of_facets()       const { return this->fs.size();     }
	inline size_t capacity_of_vertices() const { return this->vs.capacity(); }
	inline size_t capacity_of_edges()    const { return this->es.capacity(); }
	inline size_t capacity_of_facets()   const { return this->fs.capacity(); }
	inline bool empty() const { return this->vs.empty(); }
	inline void shrink_to_fit()
	{
		this->vs.shrink_to_fit();
		this->es.shrink_to_fit();
		this->fs.shrink_to_fit();
		for (typename Vertex_vector::iterator v = this->vs.begin(), end = this->vs.begin(); v != end; ++v) { v->es.shrink_to_fit(); }
	}

	inline Vertex_iterator vertices_begin() { return Vertex_iterator(&this->vs, 0); }
	inline Vertex_iterator vertices_end()   { return Vertex_iterator(&this->vs, this->vs.size()); }
	inline Edge_iterator   edges_begin()    { return Edge_iterator  (&this->es, 0); }
	inline Edge_iterator   edges_end()      { return Edge_iterator  (&this->es, this->es.size()); }
	inline Facet_iterator  facets_begin()   { return Facet_iterator (&this->fs, 0); }
	inline Facet_iterator  facets_end()     { return Facet_iterator (&this->fs, this->fs.size()); }
	inline Vertex_const_iterator vertices_begin() const { return Vertex_const_iterator(&this->vs, 0); }
	inline Vertex_const_iterator vertices_end()   const { return Vertex_const_iterator(&this->vs, this->vs.size()); }
	inline Edge_const_iterator   edges_begin()    const { return Edge_const_iterator  (&this->es, 0); }
	inline Edge_const_iterator   edges_end()      const { return Edge_const_iterator  (&this->es, this->es.size()); }
	inline Facet_const_iterator  facets_begin()   const { return Facet_const_iterator (&this->fs, 0); }
	inline Facet_const_iterator  facets_end()     const { return Facet_const_iterator (&this->fs, this->fs.size()); }

	// Add an isolated vertex at the given point to the collection and return a handle to it.
	inline Vertex_handle add_vertex(const Point& p) { this->vs.push_back(Vertex(p)); return Vertex_handle(&this->vs, this->vs.size() - 1); }
	// Add an isolated vertex to the collection and return a handle to it.
	inline Vertex_handle add_vertex(const Vertex& v) { this->vs.push_back(v); return Vertex_handle(&this->vs, this->vs.size() - 1); }
private:
	inline Edge_handle add_edge(const Edge& e) { this->es.push_back(e); return Edge_handle(&this->es, this->es.size() - 1); }
public:
	// Add a facet to the collection. The facet is defined by a set of vertices given by the
	// iterator which is an iterator of Vertex_handles. All of the edges of the facet are
	// created and the vertices and edges are updated with the necessary information about
	// incidence and identification.
	template <typename InputIter>
	inline Facet_handle add_facet(InputIter vbegin, InputIter vend, const Facet& f = Facet())
	{
		// Get a handle for the facet
		this->fs.push_back(f);
		Facet_handle fh = Facet_handle(&this->fs, this->fs.size() - 1);

		// Create the cycle of edges around the facet, adding each to the vertices lists
		Vertex_handle vh = *vbegin++;
		Edge_handle first_edge = this->add_edge(Edge(vh, fh)), eh = first_edge;
		fh->e = eh;
		for (; vbegin != vend; ++vbegin)
		{
			Vertex_handle vh = *vbegin;
			vh->es.push_back(eh);
			eh = (eh->ef = this->add_edge(Edge(vh, fh)));
		}
		eh->ef = first_edge;
		vh->es.push_back(eh);

		// Find identified edges
		eh = first_edge;
		do
		{
			eh->ei = eh;
			Vertex_handle vh = eh->vertex();
			for (typename Vertex::Edge_iterator ei = vh->edges_begin(), end = vh->edges_end(); ei != end; ++ei)
			{
				// Insert this edge into the edge circulator of identified edges
				if (eh != ei && ei->is_identified_with(eh)) { eh->ei = ei->ei; ei->ei = eh; break; }
				Edge_handle en = ei->next_around_facet();
				if (eh != en && en->is_identified_with(eh)) { eh->ei = en->ei; en->ei = eh; break; }
			}
		}
		while ((eh = eh->next_around_facet()) != first_edge);

		// Finish up
		fh->finish();
		return fh;
	}
};
