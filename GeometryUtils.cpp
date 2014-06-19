#include "GeometryUtils.hpp"

#include <vector>
#include <set>

#include "Progress.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are the few basic geometry utility functions that were not one-liners.
///////////////////////////////////////////////////////////////////////////////////////////////////

CGAL::Random Rand;

//CGAL::Cartesian_converter<Kernel, EPEC_Kernel> main_to_exact;
//CGAL::Cartesian_converter<Kernel, EPIC_Kernel> main_to_inexact;
//CGAL::Cartesian_converter<EPEC_Kernel, Kernel> exact_to_main;
//CGAL::Cartesian_converter<EPIC_Kernel, Kernel> inexact_to_main;
//CGAL::Cartesian_converter<EPEC_Kernel, EPIC_Kernel> exact_to_inexact;
//CGAL::Cartesian_converter<EPIC_Kernel, EPEC_Kernel> inexact_to_exact;

void remove_collinear_points(std::list<Point2>* pts, bool cyclic)
{
	typedef std::list<Point2>::iterator iter;
	if (pts->size() < 3) { return; }
	// The three points we are considering at one point in time are p, q, r (in that order)
	// If collinear we remove q and make: <p,q,r> = <p,r,r+1>
	// If not collinear we advance each: <p,q,r> = <q,r,r+1>
	// If the list is cyclic (e.g. making a closed polygon) then we need to make sure to also consider the triplets that cross the ends (at least two to check - more if one of those is collinear)
	iter r = pts->begin(), p = cyclic ? std::prev(pts->end()) : r++, q = r++; // <p,q,r> starts out with the first three items of the list (if acyclic) or with the last item and the first two items of the list (if cyclic) 
	for (; r != pts->end(); ++r) { if (CGAL::collinear(*p, *q, *r)) { q = r = pts->erase(q); } else { p = q; q = r; } }
	if (cyclic && pts->size() >= 3 && CGAL::collinear(*p, *q, pts->front())) { pts->erase(q); }
}

void remove_collinear_points(std::list<Point3>* pts, bool cyclic)
{
	typedef std::list<Point3>::iterator iter;
	if (pts->size() < 3) { return; }
	// The three points we are considering at one point in time are p, q, r (in that order)
	// If collinear we remove q and make: <p,q,r> = <p,r,r+1>
	// If not collinear we advance each: <p,q,r> = <q,r,r+1>
	// If the list is cyclic (e.g. making a closed polygon) then we need to make sure to also consider the triplets that cross the ends (at least two to check - more if one of those is collinear)
	iter r = pts->begin(), p = cyclic ? std::prev(pts->end()) : r++, q = r++; // <p,q,r> starts out with the first three items of the list (if acyclic) or with the last item and the first two items of the list (if cyclic) 
	for (; r != pts->end(); ++r) { if (CGAL::collinear(*p, *q, *r)) { q = r = pts->erase(q); } else { p = q; q = r; } }
	if (cyclic && pts->size() >= 3 && CGAL::collinear(*p, *q, pts->front())) { pts->erase(q); }
}

void remove_nearly_collinear_points(std::list<Point2>* pts, Kernel::FT threshold, bool cyclic)
{
	typedef std::list<Point2>::iterator iter;
	if (pts->size() < 3) { return; }
	// The three points we are considering at one point in time are p, q, r (in that order)
	// If collinear we remove q and make: <p,q,r> = <p,r,r+1>
	// If not collinear we advance each: <p,q,r> = <q,r,r+1>
	// If the list is cyclic (e.g. making a closed polygon) then we need to make sure to also consider the triplets that cross the ends (at least two to check - more if one of those is collinear)
	iter r = pts->begin(), p = cyclic ? std::prev(pts->end()) : r++, q = r++; // <p,q,r> starts out with the first three items of the list (if acyclic) or with the last item and the first two items of the list (if cyclic) 
	for (; r != pts->end(); ++r) { Kernel::FT a = CGAL::area(*p,*q,*r); if (a*a < threshold) { q = r = pts->erase(q); } else { p = q; q = r; } }
	if (cyclic && pts->size() >= 3) { Kernel::FT a = CGAL::area(*p,*q,pts->front()); if (a*a < threshold) { pts->erase(q); } }
}

void remove_nearly_collinear_points(std::list<Point3>* pts, Kernel::FT threshold, bool cyclic)
{
	typedef std::list<Point3>::iterator iter;
	if (pts->size() < 3) { return; }
	// The three points we are considering at one point in time are p, q, r (in that order)
	// If collinear we remove q and make: <p,q,r> = <p,r,r+1>
	// If not collinear we advance each: <p,q,r> = <q,r,r+1>
	// If the list is cyclic (e.g. making a closed polygon) then we need to make sure to also consider the triplets that cross the ends (at least two to check - more if one of those is collinear)
	iter r = pts->begin(), p = cyclic ? std::prev(pts->end()) : r++, q = r++; // <p,q,r> starts out with the first three items of the list (if acyclic) or with the last item and the first two items of the list (if cyclic) 
	for (; r != pts->end(); ++r) { if (CGAL::squared_area(*p, *q, *r) < threshold) { q = r = pts->erase(q); } else { p = q; q = r; } }
	if (cyclic && pts->size() >= 3 && CGAL::squared_area(*p, *q, pts->front()) < threshold) { pts->erase(q); }
}

Polygon2 facet_to_polygon2(Polyhedron3::Facet_const_handle f)
{
	const Polyhedron3::Halfedge_const_handle &a = f->halfedge(), &b = a->next(), &c = b->next();
	const Plane3 h = Plane3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point());
	std::vector<Point2> pts;
	FOR_VERTICES_AROUND_FACET(f, v) { pts.push_back(h.to_2d(v->point())); }
	return Polygon2(pts.begin(), pts.end());
}

bool is_not_degenerate(const Polyhedron3* P)
{
	// Check uniqueness of vertices
	std::set<Point3> points;
	for (Polyhedron3::Point_const_iterator p = P->points_begin(), end = P->points_end(); p != end; ++p) { if (!points.insert(*p).second) { return false; } }

	// Check geometry of facets
	for (Polyhedron3::Facet_const_iterator f = P->facets_begin(), end = P->facets_end(); f != end; ++f)
	{
		const Polygon2& p = facet_to_polygon2(f);
		if (!p.is_simple() || p.area() == 0) { return false; }
	}

	return true;
}
static bool get_vertex_he(const Polyhedron3::Vertex_const_handle v, Polyhedron3::Halfedge_around_facet_const_circulator& e)
{
	const Polyhedron3::Halfedge_around_facet_const_circulator end = e;
	do { if (e->vertex() == v) { return true; } } while (++e != end);
	return false;
}
static bool edge_intersects_face(const Polyhedron3::Halfedge_const_handle e, const Polyhedron3::Facet_const_handle f)
{
	if (e->opposite()->facet() == f) { return false; }
	const Polyhedron3::Vertex_const_handle vE1 = e->vertex(), vE2 = e->prev()->vertex();
	Polyhedron3::Halfedge_around_facet_const_circulator a = f->facet_begin();
	bool has_vertex = get_vertex_he(vE1, a) || get_vertex_he(vE2, a);
	const Polyhedron3::Halfedge_const_handle b = a->next(), c = b->next();
	const Point3 &A = a->vertex()->point(), &B = b->vertex()->point(), &C = c->vertex()->point(), &E1 = vE1->point(), &E2 = vE2->point();
	if (has_vertex)
	{
		const Point3 &E = (A == E1) ? E2 : E1;
		const Vector3 &CE = C-E, &EA = E-A;
		return CGAL::coplanar(A, B, C, E) && CGAL::cross_product(CE,EA)*CGAL::cross_product(B-E,EA) <= 0.0 && CGAL::cross_product(CE,C-B)*CGAL::cross_product(C-A,C-B) <= 0.0;
	}
	else if (CGAL::orientation(A, B, C, E1) != CGAL::orientation(A, B, C, E2))
	{
		CGAL::Orientation e2 = CGAL::orientation(A, B, E1, E2), e3 = CGAL::orientation(B, C, E1, E2), e1 = CGAL::orientation(C, A, E1, E2);
		return (e1 != CGAL::COPLANAR || e2 != CGAL::COPLANAR || e3 != CGAL::COPLANAR) && ((e1 >= CGAL::COPLANAR && e2 >= CGAL::COPLANAR && e3 >= CGAL::COPLANAR) || (e1 <= CGAL::COPLANAR && e2 <= CGAL::COPLANAR && e3 <= CGAL::COPLANAR));
	}
	return false;
}
bool is_manifold(const Polyhedron3* P)
{
	// Manifold requirements that are given and checked with P->is_valid()
	// * Every edge belongs to two faces
	// * Every vertex is surrounded by one sequence of edges and faces

#if POLYHEDRON_CACHED_NORMALS
	// * All normals face inside or outside, but not both (given, unless custom normals are present)
	Polyhedron3::Facet_const_iterator f = P->facets_begin(), end = P->facets_end();
	CGAL::Oriented_side os;
	for (; f != end; ++f) { if (f->has_normal()) { Plane3 h = facet_to_plane3(f); os = h.oriented_side(Ray3(h.point(), f->normal()).point(1)); ++f; break; } }
	for (; f != end; ++f) { if (f->has_normal()) { Plane3 h = facet_to_plane3(f); if (os != h.oriented_side(Ray3(h.point(), f->normal()).point(1))) { return false; } } }
#endif

	// * Faces only intersect each other in common edges/vertices
	// TODO: make faster
	Progress progress("Checking manifold...", P->size_of_facets());
	progress.start();
	for (Polyhedron3::Face_const_iterator f = P->facets_begin(), end = P->facets_end(); f != end; ++f)
	{
		const Bbox3 bb = facet_to_bbox3(f);
		for (Polyhedron3::Face_const_iterator g = std::next(f); g != end; ++g)
		{
			const Polyhedron3::Halfedge_const_handle a = g->facet_begin(), b = a->next(), c = b->next();
			if (CGAL::do_overlap(bb, bbox3(a->vertex()->point(), b->vertex()->point(), c->vertex()->point())) &&
				(edge_intersects_face(a, f) || edge_intersects_face(b, f) || edge_intersects_face(c, f))) { return false; }
		}
		progress.update();
	}
	progress.done();

	return true;
}

// Point-in-Polyhedron Algorithm
// Based on Shulin et al. Point-in-polyhera test with direct handling of degeneracies. Geo-spatial Information Science. 2001.
struct pip_calculator
{
	typedef Polyhedron3::Facet_const_handle value_type;
	typedef const value_type& const_reference;
	const Point3 &p, &q;
	const Plane3 h;
	bool inside;
	pip_calculator(const Point3& p, const Point3& q, const Point3& vf) : p(p), q(q), h(p, q, vf), inside(false) { }
	void push_back(const Polyhedron3::Facet_const_handle& f)
	{
		const Polyhedron3::Halfedge_const_handle &a = f->facet_begin(), &b = a->next(), &c = b->next();
		const Point3 &A = b->vertex()->point(), &B = c->vertex()->point(), &C = a->vertex()->point();
		const bool AA = h.has_on_negative_side(A), AB = h.has_on_negative_side(B), AC = h.has_on_negative_side(C);
		inside = inside != ((((AA != AB) && (AA ? Plane3(q,B,A).has_on_negative_side(p) : Plane3(q,A,B).has_on_negative_side(p)))  !=
							 ((AB != AC) && (AB ? Plane3(q,C,B).has_on_negative_side(p) : Plane3(q,B,C).has_on_negative_side(p)))) !=
							 ((AC != AA) && (AC ? Plane3(q,A,C).has_on_negative_side(p) : Plane3(q,C,A).has_on_negative_side(p))));
	}
};
bool point_in_polyhedron(const Point3& p, const FacetTree& tree)
{
	const static Vector3 dir = random_nonnull_vector();
	const Point3 q = p + dir;
	static Point3 vf(1.0,1.0,1.0);
	while (CGAL::collinear(p, q, vf)) { vf = random_point(tree.bbox()); } // arbitrary point not collinear with p and q to define h
	pip_calculator pipc(p, q, vf);
	tree.all_intersected_primitives(Ray3(p, q), std::back_inserter(pipc));
	return pipc.inside;
}
//// Original, straight-forward, method.
//// It is faster if you don't handle any degeneracies, but even just detecting degeneracies makes the times nearly the same (and sometimes the new one still beats this one out).
//// Never finished handling all the degenerate cases as the new method is way superior. Also, never abandoned the vector for a custom back_insterer.
//static bool point_in_polyhedron_y(const Point3& p, const Polyhedron3* mesh, const Facet_Tree& tree)
//{
//	const static Vector3 dir = random_nonnull_vector();
//	const Point3 q = p + dir;
//	const Ray3 R(p, q);
//	//return (tree.number_of_intersected_primitives(R) % 2) != 0; // fastest: don't even care about possible degeneracies
//	int n_half_intersections = 0;
//	//int n_intersections = 0;
//	std::vector<Polyhedron3::Facet_const_handle> intersected_facets;
//	tree.all_intersected_primitives(R, std::back_inserter(intersected_facets));
//	for (std::vector<Polyhedron3::Facet_const_handle>::const_iterator fi = intersected_facets.begin(), fend = intersected_facets.end(); fi != fend; ++fi)
//	{
//		const Polyhedron3::Halfedge_const_handle &a = (*fi)->facet_begin(), &b = a->next(), &c = b->next();
//		const Point3 &A = b->vertex()->point(), &B = c->vertex()->point(), &C = a->vertex()->point();
//		//if (CGAL::coplanar(A, B, C, p) || CGAL::coplanar(A, B, p, q) || CGAL::coplanar(B, C, p, q) || CGAL::coplanar(C, A, p, q)) // this condition is to simply fail for any degeneracies (instead of handling the easy edge intersection)
//		if (CGAL::coplanar(A, B, C, p) || R.has_on(A) || R.has_on(B) || R.has_on(C))
//		{
//			std::cerr << "! Warning: Degeneracy during the intersection test (test ray intersected vertex or is coplanar with facet)." << std::endl;
//			throw std::invalid_argument("degeneracy");
//		}
//
//		//// INCOMPLETE: Check if R intersects a vertex, it is an intersection if S is on the same side of all incident vertices
//		//if (R.has_on(A))
//		//{
//		//	Polyhedron3::Halfedge_around_vertex_const_circulator he = b->vertex_begin(), end = he;
//		//	CGAL::Oriented_side os = facet_to_plane3(he++->facet()).oriented_side(p);
//		//	while (++he != end) { if (facet_to_plane3(he->facet()).oriented_side(p) != os) { /* INCOMPLETE: no intersection*/ } }
//		//	/* INCOMPLETE: fractional intersection */
//		//}
//
//		// Check if R intersects an edge and only count 1/2 intersection if the incident facets are on opposite sides (since the incident facet will come up a second time for the second half)
//		const Plane3 ABp(A, B, p); if (ABp.has_on(q)) { if (ABp.oriented_side(C) != ABp.oriented_side(c->opposite()->next()->vertex()->point())) { ++n_half_intersections; } continue; }
//		const Plane3 BCp(B, C, p); if (BCp.has_on(q)) { if (BCp.oriented_side(A) != BCp.oriented_side(a->opposite()->next()->vertex()->point())) { ++n_half_intersections; } continue; }
//		const Plane3 CAp(C, A, p); if (CAp.has_on(q)) { if (CAp.oriented_side(B) != CAp.oriented_side(b->opposite()->next()->vertex()->point())) { ++n_half_intersections; } continue; }
//		n_half_intersections += 2;
//		//++n_intersections;
//	}			//for i
//	return ((n_half_intersections/2)%2) != 0;
//	//return (n_intersections%2) != 0;
//}


#include <iostream>
void print_tri_of_unit_cube(const Triangle3& t)
{
	static std::map<std::string, Triangle3> cube_parts;
	if (cube_parts.size() == 0)
	{
		cube_parts["a1"] = Triangle3(Point3(1,0,0), Point3(0,0,0), Point3(0,1,0));
		cube_parts["a2"] = Triangle3(Point3(0,1,0), Point3(1,1,0), Point3(1,0,0));
		cube_parts["b1"] = Triangle3(Point3(1,1,1), Point3(1,1,0), Point3(0,1,0));
		cube_parts["b2"] = Triangle3(Point3(0,1,0), Point3(0,1,1), Point3(1,1,1));
		cube_parts["c1"] = Triangle3(Point3(1,0,0), Point3(1,1,0), Point3(1,1,1));
		cube_parts["c2"] = Triangle3(Point3(1,1,1), Point3(1,0,1), Point3(1,0,0));
		cube_parts["d1"] = Triangle3(Point3(1,1,1), Point3(0,1,1), Point3(0,0,1));
		cube_parts["d2"] = Triangle3(Point3(0,0,1), Point3(1,0,1), Point3(1,1,1));
		cube_parts["e1"] = Triangle3(Point3(0,0,0), Point3(0,0,1), Point3(0,1,1));
		cube_parts["e2"] = Triangle3(Point3(0,1,1), Point3(0,1,0), Point3(0,0,0));
		cube_parts["f1"] = Triangle3(Point3(1,0,0), Point3(1,0,1), Point3(0,0,1));
		cube_parts["f2"] = Triangle3(Point3(0,0,1), Point3(0,0,0), Point3(1,0,0));
	}

	for (std::map<std::string, Triangle3>::iterator i = cube_parts.begin(), end = cube_parts.end(); i != end; ++i)
	{
		if (i->second == t) { std::cout << i->first << std::endl; return; }
	}
	std::cout << t << std::endl;
}
