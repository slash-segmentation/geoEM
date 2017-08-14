#include <GeometryTypes.hpp>

Polyhedron3* segment2cylinder(const Segment3& seg, double r=1, size_t n=32);
Polyhedron3* segments2cylinders_merged(const std::vector<Segment3>& segs, double r=1, size_t n=32);
std::vector<Polyhedron3*> segments2cylinders(const std::vector<Segment3>& segs, double r=1, size_t n=32);
