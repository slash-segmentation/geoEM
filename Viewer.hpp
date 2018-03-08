#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// A simple 3D viewer tailored for our data. This class handles mostly just the window and the
// basic setup while the classes in gl_render actually handle the 3D rendering.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "GeometryTypes.hpp"
#include "Intersection.hpp"

#include <QGLViewer/qglviewer.h>

class Viewer : public QGLViewer
{
private:
    Polyhedron3* P;
    void* glP;
    QColor edge_color, fill_color;
    
    SkeletonGraph3* S;
    void* glS;
    double line_weight;
    QColor line_color;

    Intersection intersection;
    void* glI;
    QColor i_color, ip_color;

    Point3 pt;
    double pt_radius = 0;
    void* glPt;
    QColor pt_color;

public:
    static QColor DefaultPolyhedronFillColor;
    static QColor DefaultPolyhedronEdgeColor;
    static QColor DefaultSkeletonColor;
    static QColor DefaultIntersectionColor;
    static QColor DefaultIntersectionPlaneColor;
    static QColor DefaultPointColor;

    Viewer();
    ~Viewer();
    void set_polyhedron(Polyhedron3* p, const QColor& fill_color = DefaultPolyhedronFillColor, const QColor& edge_color = DefaultPolyhedronEdgeColor);
    void set_skeleton(SkeletonGraph3* S, const double line_weight = 5.0, const QColor& line_color = DefaultSkeletonColor);
    void set_intersection(const Intersection& intersection, const QColor& intersection_color = DefaultIntersectionColor, const QColor& plane_color = DefaultIntersectionPlaneColor);
    void set_point(const Point3& pt, const double radius = 1.0, const QColor& color = DefaultPointColor);

protected:
    virtual void draw();
    virtual void init();
    virtual void keyPressEvent(QKeyEvent *e);
};

