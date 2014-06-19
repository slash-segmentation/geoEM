#pragma once

#include "GeometryTypes.h"
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// Computes the intersection between an infinite plane and a polyhedral mesh. It
// converts the results into polygons with positive or negative area.
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// A polygon representing the intersection of the 3D plane with the polyhedral mesh.
// The polygon itself is 2D (since it contained within a plane). If the polygon is
// 'open' then the first and last points don't have an edge between them. Note that if
// the area is negative if the intersection is a 'hole'.
///////////////////////////////////////////////////////////////////////////////////////
class IntersectionPolygon2 : public Polygon2
{
	bool _open;
	typedef Polygon2 Base;
public:
	IntersectionPolygon2(bool is_open = false, const Traits& p_traits = Traits()) : Base(p_traits), _open(is_open) {}
	IntersectionPolygon2(const IntersectionPolygon2& p) : Base(p), _open(p._open) {}
	template <class InputIterator> IntersectionPolygon2(InputIterator first, InputIterator last, bool is_open = false, Traits p_traits = Traits()) : Base(first, last, p_traits), _open(is_open) {}
	inline bool is_open()   const { return  this->_open; }
	inline bool is_closed() const { return !this->_open; }
};
extern template class std::vector<IntersectionPolygon2>;

///////////////////////////////////////////////////////////////////////////////////////
// A collection of intersection polygons, computed using a plane and a facet tree.
///////////////////////////////////////////////////////////////////////////////////////
class Intersection
{
	Plane3 h;
	std::vector<IntersectionPolygon2> polys;
	
	explicit inline Intersection(const Plane3& h) : h(h) { assert(!this->h.is_degenerate()); }
	//Intersection(const Intersection&);
	//Intersection& operator=(const Intersection&);

public:
	typedef std::vector<IntersectionPolygon2>::iterator iterator;
	typedef std::vector<IntersectionPolygon2>::const_iterator const_iterator;
	inline Intersection() { } // construct an empty and degenerate intersection collection
	Intersection(const FacetTree& t, const Plane3& h);
	inline bool is_degenerate() const { return this->h.is_degenerate(); }
	inline bool is_empty() const { return this->polys.size() == 0; }
	inline bool is_valid() const { return !this->h.is_degenerate() && this->polys.size() != 0; }
	inline const Plane3& plane() const { return this->h; }
	inline       IntersectionPolygon2& operator[](size_t i)       { return this->polys[i]; }
	inline const IntersectionPolygon2& operator[](size_t i) const { return this->polys[i]; }
	inline iterator       begin()       { return  this->polys.begin(); }
	inline const_iterator begin() const { return  this->polys.begin(); }
	inline iterator       end()         { return  this->polys.end(); }
	inline const_iterator end() const   { return  this->polys.end(); }
	inline void add(const IntersectionPolygon2& s) { assert(!s.is_empty()); this->polys.push_back(s); }
	inline void add(IntersectionPolygon2&& s)      { assert(!s.is_empty()); this->polys.push_back(s); }
	inline size_t count() const { return this->polys.size(); } // Number of polygons
	size_t total_count() const; // The total number of vertices across all polygons
	Bbox2 bbox() const; // The 2D bounding box of all polygons
	Kernel::FT area() const; // The sum of all areas giving the overall cross-sectional area (since holes are negative)
};
