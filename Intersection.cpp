#include "Intersection.hpp"

#include "GeometryUtils.hpp"
#include "PolyBuilder.hpp"

#include <unordered_map>
#include <vector>

template class std::vector<IntersectionPolygon2>;

bool IntersectionPolygon2::is_inside(const IntersectionPolygon2& p) const
{
    // Checks if this polygon is inside of another polygon. This assumes that both of the polygons
    // are from the same set of intersections and short-cuts the calculations considerably.
    return ::is_inside(bbox(), p.bbox()) && bounded_side(p[0]) == CGAL::ON_UNBOUNDED_SIDE;
}

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
class IntersectionPolyBuilder : public PolyBuilder<Point2, IntersectionPolygon2>
{
    typedef PolyBuilder<Point2, IntersectionPolygon2> Base;
    typedef typename Base::PartialPolygon PartialPolygon;
    typedef CGAL::Circulator_from_container<PartialPolygon> Circulator;

    const Plane3& h;
    
    std::unordered_map<const Point3, const Point2, boost::hash<Point3>> pts;
    inline const Point2& to_2d(const Point3& p) { return this->pts.insert(std::make_pair(p,::to_2d(p,this->h))).first->second; }

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
    
protected:
    void add_polygon(std::list<Point2>* pts, bool is_open) override
    {
        // This is only possible if the original mesh is not a manifold, which we disallow.
        // Pairwise check for points that are duplicates (and then should be split)
        /*for (auto i = pts->begin(), end = pts->end(); i != end; ++i)
        {
            for (auto j = std::next(i); j != end; ++j)
            {
                if (*i == *j)
                {
                    // coincident point found
                    // could handle by splitting [i,j) in another call and removing [i,j) from this one
                    throw std::runtime_error("coincident point found");
                }
            }
        }*/
    
#ifdef KERNEL_INEXACT
        // When the kernel is inexact we can inadvertently create very small "twists" in the polygon
        // making it not simple and unusable. We correct these by eliminating the twists (only as
        // long as the twist was very small, otherwise we let the problem propogate).
        if (!CGAL::is_simple_2(pts->begin(), pts->end())) { filter_polygon(*pts); }
#endif

        // Remove collinear points and make sure it is oriented clockwise (positive area)
        // We correct for holes later
        remove_collinear_points(pts);
        this->polys.push_back((CGAL::orientation_2(pts->begin(), pts->end()) == CGAL::CLOCKWISE) ?
             IntersectionPolygon2(pts->rbegin(), pts->rend(), is_open) :
             IntersectionPolygon2(pts->begin(), pts->end(), is_open));
    }

#ifdef KERNEL_INEXACT
    #define AREA_THRESHOLD 0.0001 // largest seen has been ~1.5e-5
    #define MAX_WINDOW 16 // largest seen has been 6 so this is plenty large

    static void filter_polygon(PartialPolygon& P)
    {
        // When the kernel is inexact we can inadvertently create very small "twists" in the polygon
        // making it not simple and unusable. We correct these by eliminating the twists (only as
        // long as the twist was very small, otherwise we let the problem propogate).
        //
        // The twists must contain at most MAX_WINDOW vertices and be at most AREA_THRESHOLD proportion
        // of the total area. The MAX_WINDOW value allows the algorithm used here to go from O(n^3) to
        // O(n).

        if (P.size() < 4) { return; }
        typename IntersectionPolygon2::Traits traits; // needed for CGAL::polygon_area_2

        // i goes over all points from begin() to end() - 3
        Circulator i(&P), i_end = i;
        do
        {
            Segment2 a(*i, *std::next(i)); ++i;
            Circulator j = std::next(i), j_end = P.size() < MAX_WINDOW ? std::prev(i) : std::next(j, MAX_WINDOW);
            do
            {
                Segment2 b(*j, *std::next(j)); ++j;
                if (auto intersection = CGAL::intersection(a, b))
                {
                    if (intersection->which() == 0)
                    {
                        // Segments intersect at a point
                        const Point2 x = boost::get<Point2>(*intersection);

                        // Extract the polygon that goes x -> i -> ... -> (j-1)
                        // The remaining polygon will be ... -> (i-1) -> x -> j -> ...
                        auto ix = P.insert(i.current_iterator(), x); // ix is pointing to x, the intersection
                        PartialPolygon Q(Circulator(&P, ix), j);
                        if (__is_before(j, i))
                        {
                            // remove ... -> (j-1) and i -> ... from original polygon
                            P.erase(i.current_iterator(), P.end());
                            P.erase(P.begin(), j.current_iterator());
                        }
                        else
                        {
                            // remove i -> ... -> (j-1) from original polygon
                            P.erase(i.current_iterator(), j.current_iterator());
                        }

                        // Get the areas of the two polygons
                        Kernel::FT area1 = CGAL::abs(CGAL::polygon_area_2(P.begin(), P.end(), traits));
                        Kernel::FT area2 = CGAL::abs(CGAL::polygon_area_2(Q.begin(), Q.end(), traits));

                        // Check the area and select the correct polygon
                        if (area1 < AREA_THRESHOLD*area2)
                        {
                            // Q contains the real polygon, need to replace P with its data
                            P.swap(Q);
                            i_end = Circulator(&P);
                            i = std::next(i_end);
                        }
                        else if (area2 < AREA_THRESHOLD*area1)
                        {
                            // polygon is being kept in P
                            i_end = Circulator(&P); // in case we removed the beginning elements
                            i = Circulator(&P, ix); // make i point to x
                        }
                        else { continue; } // let the problem propogate but finishing processing the polygon juts in case
                        if (P.size() < 4) { return; } // completed

                        // Update a, j, and j_end
                        a = Segment2(*std::prev(i), *i);
                        j = std::next(i); j_end = P.size() < MAX_WINDOW ? std::prev(i) : std::next(j, MAX_WINDOW);
                    }
                    else // if (intersection->which() == 1)
                    {
                        // TODO: Segments overlap, never seen and unsure how to handle, for now we shall skip it
                        continue;
                    }
                }
            } while (j != j_end);
        } while (i != i_end);
    }
    static bool __is_before(Circulator i, Circulator j)
    {
        // Checks if the node pointed to be i comes before j. This assumes that if this is true then i
        // will point to one of the first MAX_WINDOW nodes.
        Circulator x = i.min_circulator();
        if (x != j.min_circulator()) { throw std::runtime_error("i and j must be from the same container"); }
        // For very small collections do it a very easy way
        if (i.container()->size() <= MAX_WINDOW) { return std::distance(x, i) < std::distance(x, j); }
        for (size_t n = 0; n < MAX_WINDOW; ++n, ++x)
        {
            if (x == j) { return false; }
            else if (x == i) { return true; }
        }
        return false;
    }
#endif
    
public:
    typedef FacetTree::Primitive_id value_type;
    typedef const value_type& const_reference;
    
    IntersectionPolyBuilder(const Plane3& h) : h(h) {}
    
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
                    throw std::runtime_error("face intersections not implemented yet");
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
            throw std::runtime_error("this should be impossible");
        }
    }

    // Version if use "all_intersections" instead of "all_intersected_primitives" (however, it is ~3x slower and fails with inexact kernels)
    //typedef Facet_Tree::Intersection_and_primitive_id<Plane3>::Type value_type;
    //typedef const value_type& const_reference;
    //void push_back(const value_type& x)
    //{
    //  const Segment3* s = boost::get<Segment3>(&(x.first));
    //  if (s) { this->add_seg(to_2d(*s, this->h)); }

    //  // We ignore faces and points since faces will have all their edges intersected and work that way, and singleton points would give a polygon of area 0 and really complicate things
    //  // const Triangle3 *t = boost::get<Triangle3>(&(x.first));
    //  // const Point3 *p = boost::get<Point3>(&(x.first));
    //}
    //inline void add_seg(const Segment2& s) { this->add_seg(s.source(), s.target()); }
};


Intersection::Intersection(const FacetTree& tree, const Plane3& h) : h(h)
{
    assert(!h.is_degenerate());
    IntersectionPolyBuilder pb(h);
    tree.all_intersected_primitives(h, std::back_inserter(pb));
    this->polys = pb.finish_up();
}

Intersection Intersection::keep_around_point(const Point3& pt)
{
    // Create a new intersection object that only has the parts of this intersection that
    // surround the given point (or subtract from the areas that surround the given point).

    Intersection intersection(this->h);
    Point2 pt2 = to_2d(pt, this->h);
    
    // During the first pass add all of the positive polygons which contain the point to the list
    // This really should just be a single polygon, but maybe not...
    for (size_t i = 0; i < this->count(); ++i)
    {
        IntersectionPolygon2& p = (*this)[i];
        if (p.orientation() == CGAL::COUNTERCLOCKWISE && p.bounded_side(pt2) != CGAL::ON_UNBOUNDED_SIDE)
        {
            intersection.add(p);
        }
    }

    // Now add all polygons (positive or negative) which are contained within any of the already
    // found polygons
    size_t n_outside = intersection.count();
    if (n_outside == 0) { return intersection; }
    for (size_t i = 0; i < this->count(); ++i)
    {
        IntersectionPolygon2& p = (*this)[i];
        bool already_added = false;
        for (size_t j = 0; j < n_outside; ++j)
        {
            if (p == intersection[j]) { already_added = true; break; }
        }
        if (already_added) { continue; }
        for (size_t j = 0; j < n_outside; ++j)
        {
            if (p.is_inside(intersection[j]))
            {
                intersection.add(p);
                break;
            }
        }
    }
    return intersection;
}
