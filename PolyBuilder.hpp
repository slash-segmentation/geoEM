#include "GeometryTypes.hpp"

#include "GeometryUtils.hpp"

#include <list>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/circulator.h>

// Simple helper for PolyBuilder to be used as the default Polygon type
template <class Point>
class Polyline3 : public std::vector<Point>
{
    bool _open;
    typedef typename std::vector<Point> Base;
public:
    Polyline3(bool is_open=false) : Base(), _open(is_open) {}
    Polyline3(const Polyline3& p) : Base(p.begin(), p.end()), _open(p._open) {}
    template <class InputIterator> Polyline3(InputIterator first, InputIterator last, bool is_open=false) : Base(first, last), _open(is_open) {}
    bool is_open() const { return this->_open; }
};

template <class Point, class Polygon=Polyline3<Point>>
class PolyBuilder
{
    // The only requirement for Point is that it is comparable.
    // The only requirement for Polygon is that there is a constructor:
    // Polygon(InputIterator begin, InputIterator pts, bool is_open)
    // where InputIterator has a value_type of Point. A default is provided.
protected:
    typedef std::list<Point> PartialPolygon;
    typedef std::unordered_map<Point, PartialPolygon*, boost::hash<Point>> pt2partialpoly;
    pt2partialpoly starters, enders; // the partial polygons

    std::vector<Polygon> polys;

    typedef std::pair<const Point&, const Point&> simple_seg; // just a pair of points, should be "smaller" point first then "larger" (as defined by <).
    std::unordered_set<simple_seg, boost::hash<simple_seg>> have_processed; // for finding duplicates

    void add_seg(const simple_seg& seg)
    {
        // Add a segment to the set of intersections
        if (!have_processed.insert(seg).second) { return; } // already done

        auto ss = starters.find(seg.first), st = starters.find(seg.second);
        auto es = enders.find(seg.first), et = enders.find(seg.second);
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
        else if (st != starters.end())
        {
            // add the source to the beginning of the partial polygon
            PartialPolygon* p = st->second;
            starters.erase(st);
            p->push_front(seg.first);
            starters[seg.first] = p;
        }
        else if (et != enders.end())
        {
            // add the source to the end of the partial polygon
            PartialPolygon* p = et->second;
            enders.erase(et);
            p->push_back(seg.first);
            enders[seg.first] = p;
        }
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
    virtual void add_polygon(std::list<Point>* pts, bool is_open)
    {
        this->remove_collinear_points(pts);
        this->polys.push_back(Polygon(pts->begin(), pts->end(), is_open));
    }
    inline void remove_collinear_points(std::list<Point>* pts)
    {
        // does nothing unless we have Point2 or Point3
        typedef std::integral_constant<bool, std::is_base_of<Point2, Point>::value || std::is_base_of<Point3, Point>::value> is_point23;
        this->remove_collinear_points(pts, is_point23());
    }
    inline void remove_collinear_points(std::list<Point>* pts, const std::false_type &) { }
    void remove_collinear_points(std::list<Point>* pts, const std::true_type &)
    {
        typedef std::list<Point2>::iterator iter;
        if (pts->size() < 3) { return; }
        // The three points we are considering at one point in time are p, q, r (in that order)
        // If collinear we remove q and make: <p,q,r> = <p,r,r+1>
        // If not collinear we advance each: <p,q,r> = <q,r,r+1>
        // We need to make sure to also consider the triplets that cross the ends (at least two to check - more if one of those is collinear)
        iter r = pts->begin(), p = std::prev(pts->end()), q = r++; // <p,q,r> starts out with the last item and the first two items of the list (if cyclic)
        for (; r != pts->end(); ++r) { if (CGAL::collinear(*p, *q, *r)) { q = r = pts->erase(q); } else { p = q; q = r; } }
        if (pts->size() >= 3 && CGAL::collinear(*p, *q, pts->front())) { pts->erase(q); }
    }
    inline void fix_holes() { this->fix_holes(std::is_base_of<Polygon2, Polygon>()); } // does nothing unless we have a Polygon2
    inline void fix_holes(const std::false_type &) { }
    void fix_holes(const std::true_type &)
    {
        // Check to see which polygons are holes (and reverse them so they have negative area)
        size_t n = this->polys.size();
        std::vector<Bbox2>      boxes; boxes.reserve(n);
        std::vector<Kernel::FT> areas; areas.reserve(n);
        for (auto i = this->polys.begin(), end = this->polys.end(); i != end; ++i)
        {
            boxes.push_back(i->bbox()); areas.push_back(i->area());
        }
        for (size_t i = 0; i < n; ++i)
        {
            bool hole = false;
            for (size_t j = 0; j < n; ++j)
            {
                if (is_inside(boxes[j], boxes[i]) && areas[j] > areas[i] &&
                    this->polys[j].bounded_side(this->polys[i][0]) == CGAL::ON_BOUNDED_SIDE)
                {
                    hole = !hole;
                }
            }
            if (hole) { this->polys[i].reverse_orientation(); }
        }
    }

public:
    inline void add_seg(const Point& p, const Point& q) { this->add_seg(p < q ? simple_seg(p, q) : simple_seg(q, p)); }

    std::vector<Polygon>& finish_up()
    {
        // Add the open contours
        for (auto s = starters.begin(); s != starters.end(); ++s)
        {
            this->add_polygon(s->second, true);
            delete s->second;
        }

        // Possibly fix holes
        this->fix_holes();

        // All done
        return this->polys;
    }
};
