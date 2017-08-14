#include <GeometryTypes.hpp>

Polyhedron3* point2sphere(const Point3& pt, double r=1, size_t n=3);
std::vector<Polyhedron3*> points2spheres(const std::vector<Point3>& pts, double r=1, size_t n=3);
