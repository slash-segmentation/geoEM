#include "Intersection.h"

#include "GeometryUtils.h"

#include <CGAL/Boolean_set_operations_2.h>

#include <unordered_map>
#include <unordered_set>

template class std::vector<IntersectionPolygon2>;

size_t Intersection::total_count() const
{
	if (this->is_empty()) { return 0; }
	size_t n = polys[0].size();
	for (size_t i = 1; i < this->polys.size(); ++i) { n += polys[i].size(); }
	return n;
}
Bbox2 Intersection::bbox() const
{
	if (this->is_empty()) { return Bbox2(); }
	Bbox2 bb = polys[0].bbox();
	for (size_t i = 1; i < this->polys.size(); ++i) { bb += polys[i].bbox(); }
	return bb;
}
Kernel::FT Intersection::area() const
{
	if (this->is_empty()) { return 0; }
	Kernel::FT a = polys[0].area();
	for (size_t i = 1; i < this->polys.size(); ++i) { a += polys[i].area(); }
	return a;
}

// Object that builds the polygons one intersection at a time
struct PolyBuilder
{
	typedef FacetTree::Primitive_id value_type;
	typedef const value_type& const_reference;

	typedef std::list<Point2> PartialPolygon;
	typedef std::unordered_map<Point2, PartialPolygon*, boost::hash<Point2> > pt2partialpoly;

	const Plane3& h;

	typedef std::pair<const Point2&, const Point2&> simple_seg; // just a pair of points, should be "smaller" point first then "larger" (as defined by <).
	std::unordered_set<simple_seg, boost::hash<simple_seg> > have_processed; // for finding duplicates

	Intersection& ints; // the built-up intersection (completed) polygons
	pt2partialpoly starters, enders; // the partial polygons
	
	std::unordered_map<const Point3, const Point2, boost::hash<Point3> > pts;
	inline const Point2& to_2d(const Point3& p) { return this->pts.insert(std::make_pair(p,::to_2d(p,this->h))).first->second; }

	inline void add_seg(const Point2& p, const Point2& q) { this->add_seg(p < q ? simple_seg(p, q) : simple_seg(q, p)); }
	inline void add_seg(const Segment2& s) { this->add_seg(s.source(), s.target()); }
	void add_seg(const simple_seg& seg)
	{
		// Add a segment to the set of intersections

		if (!have_processed.insert(seg).second) { return; } // already done

		pt2partialpoly::iterator ss = starters.find(seg.first);
		pt2partialpoly::iterator st = starters.find(seg.second);
		pt2partialpoly::iterator es = enders.find(seg.first);
		pt2partialpoly::iterator et = enders.find(seg.second);
		if (ss != starters.end())
		{
			PartialPolygon* p = ss->second;
			starters.erase(ss);
			if (et != enders.end())
			{
				PartialPolygon* p2 = et->second;
				enders.erase(et);
				if (p == p2) { this->add_polygon(p, false); delete p; } // closed the polygon
				else { p2->splice(p2->end(), *p); enders[p2->back()] = p2; } // combined multiple partial polygons
			}
			else if (st != starters.end())
			{
				// combined multiple partial polygons (requires reversing)
				PartialPolygon* p2 = st->second;
				starters.erase(st);
				p->reverse();
				p2->splice(p2->begin(), *p);
				enders.erase(p2->front());
				starters[p2->front()] = p2;
				//enders[p2->back()] = p2; // already valid
			}
			else { p->push_front(seg.second); starters[seg.second] = p; } // add the target to the beginning of the partial polygon
		}
		else if (es != enders.end())
		{
			PartialPolygon* p = es->second;
			enders.erase(es);
			if (st != starters.end())
			{
				PartialPolygon* p2 = st->second;
				starters.erase(st);
				if (p == p2) { this->add_polygon(p, false); delete p; } // closed the polygon
				else { p->splice(p->end(), *p2); enders[p->back()] = p; } // combined multiple partial polygons
			}
			else if (et != enders.end())
			{
				// combined multiple partial polygons (requires reversing)
				PartialPolygon* p2 = et->second;
				enders.erase(et);
				p->reverse();
				p2->splice(p2->end(), *p);
				starters.erase(p2->back());
				//starters[p2->front()] = p2; // already valid
				enders[p2->back()] = p2;
			}
			else { p->push_back(seg.second); enders[seg.second] = p; } // add the target to the end of the partial polygon 
		}
		else if (st != starters.end()) { PartialPolygon* p = st->second; starters.erase(st); p->push_front(seg.first); starters[seg.first] = p; } // add the source to the beginning of the partial polygon
		else if (et != enders.end())   { PartialPolygon* p = et->second;   enders.erase(et); p->push_back (seg.first);   enders[seg.first] = p; } // add the source to the end of the partial polygon
		else
		{
			// neither endpoint currently exists, create a new partial polygon
			PartialPolygon* p = new PartialPolygon();
			p->push_front(seg.first);
			p->push_back(seg.second);
			starters[seg.first] = p;
			enders[seg.second] = p;
		}
	}

	void add_polygon(std::list<Point2>* pts, bool is_open)
	{
		// This is only possible if the original mesh is not a manifold, which we disallow.
		//// Pairwise check for points that are duplicates (and then should be split)
		//for (std::list<Point2>::const_iterator i = pts->begin(), end = pts->end(); i != end; ++i)
		//{
		//	for (std::list<Point2>::const_iterator j = std::next(i); j != end; ++j)
		//	{
		//		if (*i == *j) { /* coincident point found, split by doing [i,j) in another call and removing [i,j) from this one */ }
		//	}
		//}

		// Remove collinear points and make sure it is oriented clockwise (positive area)
		remove_collinear_points(pts);
		this->ints.add((CGAL::orientation_2(pts->begin(),pts->end()) == CGAL::CLOCKWISE) ? IntersectionPolygon2(pts->rbegin(), pts->rend(), is_open) : IntersectionPolygon2(pts->begin(), pts->end(), is_open));
	}

	void finish_up()
	{
		// Add the open contours
		for (pt2partialpoly::iterator s = starters.begin(); s != starters.end(); ++s) { this->add_polygon(s->second, true); delete s->second; }

		// Check to see which polygons are holes (and reverse them so they have negative area)
		size_t n = ints.count();
		std::vector<Bbox2>      boxes; boxes.reserve(n);
		std::vector<Kernel::FT> areas; areas.reserve(n);
		for (Intersection::const_iterator i = ints.begin(), end = ints.end(); i != end; ++i) { boxes.push_back(i->bbox()); areas.push_back(i->area()); }
		for (size_t i = 0; i < n; ++i)
		{
			bool hole = false;
			for (size_t j = 0; j < n; ++j) { if (is_inside(boxes[j], boxes[i]) && areas[j] > areas[i] && ints[j].bounded_side(ints[i][0]) == CGAL::ON_BOUNDED_SIDE) { hole = !hole; } }
			if (hole) { ints[i].reverse_orientation(); }
		}
	}

	const Point2& intersection(const Point3& p, const Point3& q)
	{
		// Calculate when line segment PQ intersections H (member variable)
		// Assumption is that P is on the negative side of the plane while Q is on the positive side of plane H
		Vector3 n = this->h.orthogonal_vector(), l = q-p;
		Kernel::FT d = ((this->h.point()-p)*n) / (l*n); // distance along line defined by p and (p-q)
		return to_2d(p+l*d);
	}
	inline const Point2& intersection(const Point3& p, const Point3& q, CGAL::Oriented_side p_os)
	{
		// Calculate when line segment PQ intersections H (member variable) given the side of the plane H that P is on
		// Assumption is that P and Q are on opposite sides of the plane H and neither on the plane H
		return p_os == CGAL::ON_NEGATIVE_SIDE ? intersection(p, q) : intersection(q, p);
	}
	
	PolyBuilder(Intersection& ints) : h(ints.plane()), ints(ints) { }

	void push_back(const value_type& f)
	{
		// Called with the next intersected primitive, must be called push_back to work with the back_inserter

		const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next();
		const Point3 &A = a->vertex()->point(), &B = b->vertex()->point(), &C = c->vertex()->point();
		const CGAL::Oriented_side Aos = this->h.oriented_side(A), Bos = this->h.oriented_side(B), Cos = this->h.oriented_side(C);

		if (Aos == CGAL::ON_ORIENTED_BOUNDARY)
		{
			if (Bos == CGAL::ON_ORIENTED_BOUNDARY)
			{
				if (Cos == CGAL::ON_ORIENTED_BOUNDARY)
				{
					// TODO: face ABC
				}
				else { this->add_seg(to_2d(A),to_2d(B)); } // edge AB
			}
			else if (Cos == CGAL::ON_ORIENTED_BOUNDARY) { this->add_seg(to_2d(A),to_2d(C)); } // edge AC
			else if (Bos != Cos) { this->add_seg(to_2d(A),intersection(B,C,Bos)); } // segment A to [point between BC]
			// else point A (ignored)
		}
		else if (Bos == CGAL::ON_ORIENTED_BOUNDARY)
		{
			if (Cos == CGAL::ON_ORIENTED_BOUNDARY) { this->add_seg(to_2d(B),to_2d(C)); } // edge BC
			else if (Aos != Cos) { this->add_seg(to_2d(B),intersection(A,C,Aos)); } // segment B to [point between AC]
			// else point B (ignored)
		}
		else if (Cos == CGAL::ON_ORIENTED_BOUNDARY)
		{
			if (Aos != Cos) { this->add_seg(to_2d(C),intersection(A,B,Aos)); } // segment C to [point between AB]
			// else point C (ignored)
		}
		else if (Aos == Bos) { this->add_seg(intersection(A,C,Aos),intersection(B,C,Bos)); } // segment [point between AC] to [point between BC]
		else if (Aos == Cos) { this->add_seg(intersection(A,B,Aos),intersection(B,C,Bos)); } // segment [point between AB] to [point between BC]
		else if (Bos == Cos) { this->add_seg(intersection(A,B,Aos),intersection(A,C,Aos)); } // segment [point between AB] to [point between AC]
		else
		{
			// impossible?
		}
	}

	// Version if use "all_intersections" instead of "all_intersected_primitives" (however, it is ~3x slower and fails with inexact kernels)
	//typedef Facet_Tree::Intersection_and_primitive_id<Plane3>::Type value_type;
	//typedef const value_type& const_reference;
	//void push_back(const value_type& x)
	//{
	//	const Segment3* s = boost::get<Segment3>(&(x.first));
	//	if (s) { this->add_seg(to_2d(*s, this->h)); }

	//	// We ignore faces and points since faces will have all their edges intersected and work that way, and singleton points would give a polygon of area 0 and really complicate things
	//	// const Triangle3 *t = boost::get<Triangle3>(&(x.first));
	//	// const Point3 *p = boost::get<Point3>(&(x.first));
	//}
};

Intersection::Intersection(const FacetTree& tree, const Plane3& h) : h(h)
{
	assert(!h.is_degenerate());
	PolyBuilder pb(*this);
	tree.all_intersected_primitives(h, std::back_inserter(pb));
	pb.finish_up();
}
