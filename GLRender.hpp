#ifndef _GL_RENDER_
#define _GL_RENDER_

#include "GeometryTypes.hpp"
#include "PolygonSorter.hpp"
#include "Intersection.hpp"
#include <QColor>

///////////////////////////////////////////////////////////////////////////////////////////////////
// The functions for efficiently rendering objects to the screen. They work by buffering all
// necessary data to the graphics card memory when constructed and when actually rendering do the
// bare minimum to update the view. This results in very fast rendering even for very large
// objects. Note that the current system disables lighting and depth testing in all cases and only
// properly sorts the polyhedron faces relative to themselves, so rendering artifacts are expected.
///////////////////////////////////////////////////////////////////////////////////////////////////

struct GlPolyhedron
{
    size_t nverts, nedges, nfaces;
    unsigned int bufs[2];
    unsigned int* edges;
    unsigned int* faces;
    PolygonSorter* ps;

    GlPolyhedron(const Polyhedron3* P);
    ~GlPolyhedron();
    void render_edges(const QColor& color);
    void render_faces(const Ray3& view, const QColor& color);
};

struct GlSkeleton
{
    size_t nverts, nedges;
    unsigned int bufs[1];
    unsigned int* edges;

    GlSkeleton(const Skeleton3* S);
    ~GlSkeleton();
    void render(const double weight, const QColor& color);
};

struct GlIntersection
{
    Point3 plane_corners[4];
    size_t nintersections, *npoints;
    unsigned int bufs[1];
    unsigned int** polys;

    GlIntersection(const Intersection& I);
    ~GlIntersection();
    void render_plane(const QColor& color);
    void render_polygons(const QColor& color);
};

struct GlPoint
{
    Point3 pt;
    double radius;
    GlPoint(const Point3& pt, double radius);
    ~GlPoint();
    void render(const QColor& color);
};

void gl_setup();

#endif
