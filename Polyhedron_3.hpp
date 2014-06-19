///////////////////////////////////////////////////////////////////////////////////////////////////
// An extended Polyhedron class (from CGAL) with the following extra features:
//   * facet and vertex normal calculations (and optionally caching)
//   * internal data structures using vectors instead of lists (optional)
//   * arbitrary extra data can be stored in each facet, halfedge, and vertex
//   * the bounding box of the polyhedron can be calculated and is cached
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef POLYHEDRON_3_H
#define POLYHEDRON_3_H

// By defining this normals are cached in the vertices and facets of polyhedron.
// This has the benefit of possibly speeding up and increasing the quality of visual rendering at
// the expense of increased loading times and memory usage. By defining it the times when the
// polyhedron structure can be modified are limited.
//#define POLYHEDRON_CACHED_NORMALS

// This makes everything faster but disables the use of polyhedron surface triangulation for
// arbitrary polygons (only meshes with triangles and quads can be loaded).
#define POLYHEDRON_USE_VECTOR


// Cached normals notes:
// The normals can be set directly or calculated (and cached) upon first use. Calling inside_out()
// flips all normals. In general this modified Polyhedron class does NOT support modifications once
// fully created and the normal() functions of facets and vertices should not be called until it is
// fully created unless has_normal() returns true.

#include <CGAL/Direction_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polyhedron_3.h>
#ifdef POLYHEDRON_USE_VECTOR
#include <CGAL/HalfedgeDS_vector.h>
#endif
#include <CGAL/IO/Polyhedron_iostream.h>

// Extras class (adds generic extras to Face, halfedge, and vertices)
// To use it, locally (to the compilation unit) define a the class Facet_extra, Edge_extra, or Vertex_extra.
template <class Data>
class Extra
{
private:
	Data* _extra;
public:
	inline Extra() : _extra(nullptr) { }
	inline bool has_extra() const { return this->_extra != nullptr; }
	inline Data*& extra() { return this->_extra; }
	inline const Data* const extra() const { return this->_extra; }
};

// Forward declarations of functions in GeometryUtils.h
template <class Kernel> inline CGAL::Vector_3<Kernel> to_vector(const CGAL::Direction_3<Kernel> &d);
template <class Kernel> inline CGAL::Vector_3<Kernel> normalized(const CGAL::Vector_3<Kernel> &v);
template <class Kernel> inline CGAL::Direction_3<Kernel> normalized(const CGAL::Direction_3<Kernel> &d);

// Face class
struct Facet_extra;
template <class Refs, class Traits>
class HalfedgeDS_face : public CGAL::HalfedgeDS_face_base<Refs, CGAL::Tag_true, CGAL::Tag_false>, public Extra<Facet_extra>
{
public:
	typedef typename Traits::Kernel      Kernel;
	typedef typename Kernel::Direction_3 Direction;

private:
	typedef typename CGAL::HalfedgeDS_face_base<Refs, CGAL::Tag_true, CGAL::Tag_false> Base;
	typedef HalfedgeDS_face<Refs, Traits> Self;
	typedef typename CGAL::I_Polyhedron_facet<Self> Facet;
	static Direction compute_normal(const Facet &f)
	{
		typename Facet::Halfedge_around_facet_const_circulator he = f.facet_begin();
		return Direction(normalized(CGAL::normal(he->prev()->vertex()->point(), he->vertex()->point(), he->next()->vertex()->point())));
	}
	
#ifdef POLYHEDRON_CACHED_NORMALS
	Direction _normal;
#endif

public:
#ifdef POLYHEDRON_CACHED_NORMALS
	inline static Direction const& NO_DIRECTION() { static Direction nd(0,0,0); return nd; }
	HalfedgeDS_face() : Base(), Extra(), _normal(NO_DIRECTION()) {}
	HalfedgeDS_face(const Direction& n) : Base(), Extra(), _normal(normalized(d)) {}
	Direction&       normal()       { if (!has_normal()) { this->_normal = compute_normal(*static_cast<Facet*>(this)); } return _normal; }
	const Direction& normal() const { if (!has_normal()) { Self* t = const_cast<Self*>(this); t->_normal = compute_normal(*static_cast<Facet*>(t)); } return _normal; }
	const Direction& calculated_normal() const { return compute_normal(*static_cast<Facet*>(t)); }
	void             set_normal(const Direction& n) { _normal = normalized(n); }
	inline bool      has_normal() const             { return _normal != NO_DIRECTION(); }
#else
	HalfedgeDS_face() : Base(), Extra() {}
	Direction normal() const { return compute_normal(*static_cast<const Facet*>(this)); }
#endif
};

// Halfedge class
struct Halfedge_extra;
template <class Refs, class Traits>
class HalfedgeDS_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs, CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true>, public Extra<Halfedge_extra>
{
public:
	typedef typename Traits::Kernel      Kernel;
	
private:
	typedef typename CGAL::HalfedgeDS_halfedge_base<Refs, CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true> Base;
	typedef HalfedgeDS_halfedge<Refs, Traits> Self;
	typedef typename CGAL::I_Polyhedron_halfedge<Self> Halfedge;

public:
	HalfedgeDS_halfedge() : Base(), Extra() {}
};

// Vertex class
struct Vertex_extra;
template <class Refs, class Traits>
class HalfedgeDS_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, typename Traits::Point_3>, public Extra<Vertex_extra>
{
public:
	typedef typename Traits::Kernel      Kernel;
	typedef typename Traits::Point_3     Point;
	typedef typename Kernel::Direction_3 Direction;

private:
	typedef typename CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> Base;
	typedef HalfedgeDS_vertex<Refs, Traits> Self;
	typedef typename CGAL::I_Polyhedron_vertex<Self> Vertex;
	static Direction compute_normal(Vertex &v)
	{
		typename Kernel::Vector_3 normal = CGAL::NULL_VECTOR;
		typename Vertex::Halfedge_around_vertex_circulator he = v.vertex_begin(), end(he);
		CGAL_For_all(he, end) { if (!he->is_border()) { normal = normal + to_vector(he->facet()->normal()); } }
		return Direction(normalized(normal));
	}
	static Direction compute_normal(const Vertex &v)
	{
		typename Kernel::Vector_3 normal = CGAL::NULL_VECTOR;
		typename Vertex::Halfedge_around_vertex_const_circulator he = v.vertex_begin(), end(he);
		CGAL_For_all(he, end) { if (!he->is_border()) { normal = normal + to_vector(he->facet()->normal()); } }
		return Direction(normalized(normal));
	}


#ifdef POLYHEDRON_CACHED_NORMALS
	Direction _normal;
#endif

public:
#ifdef POLYHEDRON_CACHED_NORMALS
	inline static Direction const& NO_DIRECTION() { static Direction nd(0,0,0); return nd; }
	HalfedgeDS_vertex() : Base(), Extra(), _normal(NO_DIRECTION()) {}
	HalfedgeDS_vertex(const Point& p) : Base(p), Extra(), _normal(NO_DIRECTION()) {}
	HalfedgeDS_vertex(const Point& p, const Direction &n) : Base(p), Extra(), _normal(normalized(n)) {}
	Direction&       normal()		{ if (!has_normal()) { this->_normal = compute_normal(*static_cast<Vertex*>(this)); } return _normal; }
	const Direction& normal() const { if (!has_normal()) { Self* t = const_cast<Self*>(this); t->_normal = compute_normal(*static_cast<Vertex*>(t)); } return _normal; }
	const Direction& calculated_normal() const { return compute_normal(*static_cast<Vertex*>(t)); }
	void             set_normal(const Direction& n) { _normal = normalized(n); }
	inline bool      has_normal() const             { return _normal != NO_DIRECTION(); }
#else
	HalfedgeDS_vertex() : Base(), Extra() {}
	HalfedgeDS_vertex(const Point& p) : Base(p), Extra() {}
	Direction normal() const { return compute_normal(*static_cast<const Vertex*>(this)); }
#endif
};

// Items class
struct Polyhedron_items_3 : public CGAL::Polyhedron_items_3
{
	template <class Refs, class Traits>
	struct Face_wrapper { typedef HalfedgeDS_face<Refs, Traits> Face; };
	
	template <class Refs, class Traits>
	struct Halfedge_wrapper { typedef HalfedgeDS_halfedge<Refs, Traits> Halfedge; };

	template <class Refs, class Traits>
	struct Vertex_wrapper { typedef HalfedgeDS_vertex<Refs, Traits> Vertex; };
};

// Polyhedron class
template<class Traits>
#ifdef POLYHEDRON_USE_VECTOR
class Polyhedron_3 : public CGAL::Polyhedron_3<Traits, Polyhedron_items_3, CGAL::HalfedgeDS_vector>
{
	typedef typename CGAL::Polyhedron_3<Traits, Polyhedron_items_3, CGAL::HalfedgeDS_vector> Base;
#else
class Polyhedron_3 : public CGAL::Polyhedron_3<Traits, Polyhedron_items_3>
{
	typedef typename CGAL::Polyhedron_3<Traits, Polyhedron_items_3> Base;
#endif

	CGAL::Bbox_3 _bbox;

public:
	typedef Polyhedron_3<Traits> Self;

	Polyhedron_3(const Traits& trts = Traits()) : Base(trts), _bbox(0,0,0,0,0,0) {}
	Polyhedron_3(size_t v, size_t h, size_t f, const Traits& traits = Traits()) : Base(v,h,f,traits), _bbox(0,0,0,0,0,0) {}

#ifdef POLYHEDRON_CACHED_NORMALS
	void inside_out()
	{
		Base::inside_out();
		for (Base::Vertex_iterator i = Base::vertices_begin(), end = Base::vertices_end(); i != end; ++i) { if (i->has_normal()) { i->set_normal(-i->normal()); } }
		for (Base::Facet_iterator  i = Base::facets_begin(),   end = Base::facets_end();   i != end; ++i) { if (i->has_normal()) { i->set_normal(-i->normal()); } }
	}
#endif

	const CGAL::Bbox_3& bbox() const
	{
		if (this->_bbox == CGAL::Bbox_3(0,0,0,0,0,0)) { const_cast<Self*>(this)->_bbox = bbox3_from_points(this->points_begin(), this->points_end()); }
		return this->_bbox;
	}

	bool extras_are_not_used() const
	{
		for (typename Base::Vertex_const_iterator   i = Base::vertices_begin(),  end = Base::vertices_end();  i != end; ++i) { if (i->has_extra()) { return false; } }
		for (typename Base::Facet_const_iterator    i = Base::facets_begin(),    end = Base::facets_end();    i != end; ++i) { if (i->has_extra()) { return false; } }
		for (typename Base::Halfedge_const_iterator i = Base::halfedges_begin(), end = Base::halfedges_end(); i != end; ++i) { if (i->has_extra()) { return false; } }
		return true;
	}
};

#endif
