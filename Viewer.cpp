#include "Viewer.hpp"

#include "GeometryUtils.hpp"

#include <QKeyEvent>

#include <CGAL/Min_sphere_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_3.h>

#include "GLRender.hpp"

QColor Viewer::DefaultPolyhedronFillColor(0, 128, 0, 64);
QColor Viewer::DefaultPolyhedronEdgeColor(0, 0, 0);
QColor Viewer::DefaultSkeletonColor(255, 0, 255);
QColor Viewer::DefaultIntersectionColor(255, 0, 0);
QColor Viewer::DefaultIntersectionPlaneColor(0, 127, 255, 32);
QColor Viewer::DefaultPointColor(255, 0, 0);

inline qglviewer::Vec pt2vec(const Point3& p) { return qglviewer::Vec(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())); }
inline Point3 vec2pt(const qglviewer::Vec& v) { return Point3(v.x, v.y, v.z); }
inline qglviewer::Vec vec2vec(const Vector3& v) { return qglviewer::Vec(CGAL::to_double(v.x()), CGAL::to_double(v.y()), CGAL::to_double(v.z())); }
inline Vector3 vec2vec(const qglviewer::Vec& v) { return Vector3(v.x, v.y, v.z); }

void Viewer::set_polyhedron(Polyhedron3* P, const QColor& fill_color, const QColor& edge_color)
{
	if (P->empty()) { /* TODO: error! */ return; }
	
	// Cleanup previous polyhedron
	if (this->P)   { delete this->P; this->P = nullptr; }
	if (this->glP) { delete (GlPolyhedron*)this->glP; this->glP = nullptr; }

	// Set the polyhedron object
	this->P = P;
	this->fill_color = fill_color;
	this->edge_color = edge_color;

	// Setup the scene
	CGAL::Min_sphere_d<CGAL::Min_sphere_annulus_d_traits_3<Kernel>> ms(P->points_begin(), P->points_end());
	Point3 center = ms.center();
	double radius = dbl_sqrt(ms.squared_radius());
	this->setSceneCenter(pt2vec(center));
	this->setSceneRadius(radius);

	// Setup the camera
	qglviewer::Camera *cam = this->camera();
	//cam->setType(qglviewer::Camera::ORTHOGRAPHIC);
	cam->lookAt(this->sceneCenter());
	cam->showEntireScene();
	cam->setFOVToFitScene();
}

void Viewer::set_skeleton(SkeletonGraph3* S, const double line_weight, const QColor& line_color)
{
	if (this->glS) { delete (GlSkeleton*)this->glS; this->glS = nullptr; }
	this->S = S;
	this->line_weight = line_weight;
	this->line_color = line_color;
}

void Viewer::set_intersection(const Intersection& intersection, const QColor& intersection_color, const QColor& plane_color)
{
	if (this->glI) { delete (GlIntersection*)this->glI; this->glI = nullptr; }
	this->intersection = intersection;
	this->i_color = intersection_color;
	this->ip_color = plane_color;
}

void Viewer::set_point(const Point3& pt, const double radius, const QColor& color)
{
	if (this->glPt) { delete (GlPoint*)this->glPt; this->glPt = nullptr; }
	this->pt = pt;
	this->pt_radius = radius;
	this->pt_color = color;
}

void Viewer::draw()
{
	// Create GL objects (mostly loading buffers into graphics card)
	if (!this->glP && this->P) { this->glP = new GlPolyhedron(this->P); }
	if (!this->glS && this->S) { this->glS = new GlSkeleton(this->S); }
	if (!this->glI && this->intersection.is_valid()) { this->glI = new GlIntersection(this->intersection); }
	if (!this->glPt && this->pt_radius) { this->glPt = new GlPoint(this->pt, this->pt_radius); }

	const Ray3 view(vec2pt(this->camera()->position()), vec2vec(this->camera()->orientation()*qglviewer::Vec(0,0,-1)));
	
	// Render!
	if (this->glPt)
	{
		((GlPoint*)this->glPt)->render(this->pt_color);
	}
	if (this->glI)
	{
		((GlIntersection*)this->glI)->render_plane(this->ip_color);
		((GlIntersection*)this->glI)->render_polygons(this->i_color);
	}
	if (this->glS)
	{
		((GlSkeleton*)this->glS)->render(this->line_weight, this->line_color);
	}
	if (this->glP)
	{
		((GlPolyhedron*)this->glP)->render_faces(view, this->fill_color);
		((GlPolyhedron*)this->glP)->render_edges(this->edge_color);
	}
}

Viewer::Viewer() : P(nullptr), glP(nullptr), glS(nullptr), glI(nullptr), glPt(nullptr)
{
	setKeyDescription(Qt::Key_I, "Flip the polyhedron inside/out");
}

Viewer::~Viewer()
{
	if (this->P)   { delete this->P; this->P = nullptr; }
	if (this->glP) { delete (GlPolyhedron*)this->glP; this->glP = nullptr; }
	if (this->glI) { delete (GlIntersection*)this->glI; this->glI = nullptr; }
	if (this->glPt) { delete (GlPoint*)this->glPt; this->glPt = nullptr; }
}

void Viewer::init()
{
	gl_setup();
	setBackgroundColor(Qt::white);
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_I)
	{
		inside_out(this->P);
		if (this->glP) { delete (GlPolyhedron*)this->glP; this->glP = nullptr; }
		update(); // Refresh display
	}
	else
	{
		QGLViewer::keyPressEvent(e);
	}
}
