#include "GeometryUtils.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// These are the few basic geometry utility functions that were not one-liners.
///////////////////////////////////////////////////////////////////////////////////////////////////

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
    if (cyclic && pts->size() >= 3) { Kernel::FT a = CGAL::area(*p,*q, pts->front()); if (a*a < threshold) { pts->erase(q); } }
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
