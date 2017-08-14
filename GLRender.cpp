#include "GLRender.hpp"

#include "GeometryUtils.hpp"
#include "PolygonSorter.hpp"

////////// Basics //////////
#if defined(_WIN32)
// Windows OpenGL
#include <GL/gl.h>
#include "glext.h"
#define GL_GET_PROC_ADDR(n) wglGetProcAddress(n) // requires a valid OpenGL context
#elif defined(__APPLE__)
// Apple OpenGL (does not need anything special for extensions)
#include <OpenGL/gl.h>
#else
// Linux OpenGL
#include <GL/glx.h>
#include <GL/glxext.h>
#define GL_GET_PROC_ADDR(n) glXGetProcAddressARB((const GLubyte *)n) // does not require a valid OpenGL context
#endif
static PFNGLGENBUFFERSPROC glGenBuffers = nullptr;
static PFNGLBINDBUFFERPROC glBindBuffer = nullptr;
static PFNGLBUFFERDATAPROC glBufferData = nullptr;
static PFNGLDELETEBUFFERSPROC glDeleteBuffers = nullptr;
static PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray = nullptr;
static PFNGLDISABLEVERTEXATTRIBARRAYPROC glDisableVertexAttribArray = nullptr;
static PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer = nullptr;
static void _gl_init()
{
#ifndef __APPLE__
	if (glGenBuffers == nullptr)
	{
		glGenBuffers = (PFNGLGENBUFFERSPROC)GL_GET_PROC_ADDR("glGenBuffers");
		glBindBuffer = (PFNGLBINDBUFFERPROC)GL_GET_PROC_ADDR("glBindBuffer");
		glBufferData = (PFNGLBUFFERDATAPROC)GL_GET_PROC_ADDR("glBufferData");
		glDeleteBuffers = (PFNGLDELETEBUFFERSPROC)GL_GET_PROC_ADDR("glDeleteBuffers");
		
		glEnableVertexAttribArray = (PFNGLENABLEVERTEXATTRIBARRAYPROC)GL_GET_PROC_ADDR("glEnableVertexAttribArray");
		glDisableVertexAttribArray = (PFNGLDISABLEVERTEXATTRIBARRAYPROC)GL_GET_PROC_ADDR("glDisableVertexAttribArray");
		glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC)GL_GET_PROC_ADDR("glVertexAttribPointer");
	}
#endif
}
void gl_setup()
{
	_gl_init();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
}
inline void _gl_set_color(const QColor &c) { glColor4d(c.redF(), c.greenF(), c.blueF(), c.alphaF()); }
//inline void _gl_set_normal(const Direction3 &n) { glNormal3d(CGAL::to_double(n.dx()), CGAL::to_double(n.dy()), CGAL::to_double(n.dz())); }
inline void _gl_add_point(const Point3 &p) { glVertex3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())); }

////////// GlPolyhedron //////////
GlPolyhedron::GlPolyhedron(const Polyhedron3* P) : nverts(0), nedges(0), nfaces(0), edges(nullptr), faces(nullptr), ps(nullptr)
{
	size_t i;
	this->bufs[0] = 0;
	this->bufs[1] = 0;

	// Get all vertices and their normals
	this->nverts = P->size_of_vertices();
	GLdouble* verts = new GLdouble[this->nverts*3];
	GLdouble* norms = new GLdouble[this->nverts*3];
	i = 0;
	for (Polyhedron3::Vertex_const_iterator V = P->vertices_begin(), end = P->vertices_end(); V != end; ++V)
	{
		const Point3& v = V->point();
		const Direction3& n = normal(V);
		verts[i  ] = CGAL::to_double(v.x());
		norms[i++] = CGAL::to_double(n.dx());
		verts[i  ] = CGAL::to_double(v.y());
		norms[i++] = CGAL::to_double(n.dy());
		verts[i  ] = CGAL::to_double(v.z());
		norms[i++] = CGAL::to_double(n.dz());
	}

	IteratorReverseLookup<Polyhedron3::Vertex_const_handle> vertex_lookup(P->vertices_begin(), P->size_of_vertices());

	// Get all edges
	this->nedges = P->size_of_halfedges();
	this->edges = new unsigned int[this->nedges];
	i = 0;
	for (Polyhedron3::Edge_const_iterator HE = P->edges_begin(), end = P->edges_end(); HE != end; ++HE)
	{
		this->edges[i++] = (unsigned int)vertex_lookup[HE->vertex()];
		this->edges[i++] = (unsigned int)vertex_lookup[HE->opposite()->vertex()];
	}

	// Get all faces
	this->nfaces = P->size_of_facets();
	this->faces = new unsigned int[this->nfaces*3]; // used when getting sorted face order, not below
	unsigned int** faces = new unsigned int*[this->nfaces];
	i = 0;
	for (Polyhedron3::Facet_const_iterator F = P->facets_begin(), end = P->facets_end(); F != end; ++F)
	{
		const Polyhedron3::Halfedge_const_handle &a = F->halfedge(), &b = a->next(), &c = b->next();
		faces[i] = new unsigned int[3];
		faces[i][0] = (unsigned int)vertex_lookup[a->vertex()];
		faces[i][1] = (unsigned int)vertex_lookup[b->vertex()];
		faces[i][2] = (unsigned int)vertex_lookup[c->vertex()];
		++i;
	}

	// Save vertex and normal data to the graphics card
	glGenBuffers(2, this->bufs);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble)*nverts*3, verts, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[1]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble)*nverts*3, norms, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Create the polygon sorter
	this->ps = new PolygonSorter(verts, faces, this->nfaces);

	// Cleanup temporary data
	delete[] verts;
	delete[] norms;
	for (size_t i = 0; i < this->nfaces; ++i) { if (faces[i]) { delete[] faces[i]; faces[i] = nullptr; } }
	delete[] faces;
}
GlPolyhedron::~GlPolyhedron()
{
	glDeleteBuffers(2, this->bufs); this->bufs[0] = 0; this->bufs[1] = 0;
	if (this->edges) { delete[] this->edges; this->edges = nullptr; }
	if (this->faces) { delete[] this->faces; this->faces = nullptr; }
	if (this->ps) { delete this->ps; this->ps = nullptr; }
}
void GlPolyhedron::render_edges(const QColor& color)
{
	// Setup the GL environment
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glLineWidth(0.25f);
	_gl_set_color(color);

	// Draw all edges
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, 0);
	glDrawElements(GL_LINES, (GLsizei)this->nedges, GL_UNSIGNED_INT, this->edges);
	glDisableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//// Indirect/slow/easy way
	//glBegin(GL_LINES);
	//for (Polyhedron3::Edge_const_iterator he = P->edges_begin(), end = P->edges_end(); he != end; ++he)
	//{
	//	_gl_add_point(he->vertex()->point());
	//	_gl_add_point(he->opposite()->vertex()->point());
	//}
	//glEnd();
}
void GlPolyhedron::render_faces(const Ray3& view, const QColor& color)
{
	// Calculate the ordering of the faces
	this->ps->query(view.source(), view.direction(), this->faces);

	// Setup the GL environment
	//glEnable(GL_LIGHTING); // TODO: alpha and lighting
	glDisable(GL_DEPTH_TEST);
	_gl_set_color(color);
	
	// Draw all faces
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[1]);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(1);
	glDrawElements(GL_TRIANGLES, (GLsizei)this->nfaces*3, GL_UNSIGNED_INT, this->faces);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
		
	//// Indirect/slow/easy way
	//glBegin(GL_TRIANGLES);
	//for (; f != end; ++f)
	//{
	//	// revolve around current face to get vertices
	//  VERTICES_AROUND_FACET(f, v)
	//	{
	//		// Gouraud (smooth) shading has 1 normal per vertex
	//		_gl_set_normal(normal(v));
	//		_gl_add_point(v->point());
	//	}
	//}
	//glEnd();
}


////////// GlSkeleton //////////
GlSkeleton::GlSkeleton(const Skeleton3* S) : nverts(boost::num_vertices(*S)), nedges(boost::num_edges(*S)), edges(nullptr)
{
	size_t i;
	this->bufs[0] = 0;
	
	// Get all vertices
	GLdouble* verts = new GLdouble[this->nverts*3];
	i = 0;
	Skeleton3::vertex_iterator V, Vend;
	std::unordered_map<Skeleton3::vertex_descriptor, size_t> vertex_lookup;
	for (boost::tie(V, Vend) = boost::vertices(*S); V != Vend; ++V)
	{
		vertex_lookup[*V] = i;
		const Point3& v = (*S)[*V].point;
		verts[3*i+0] = CGAL::to_double(v.x());
		verts[3*i+1] = CGAL::to_double(v.y());
		verts[3*i+2] = CGAL::to_double(v.z());
		++i;
	}
	
	// Get all edges
	this->edges = new unsigned int[2*this->nedges];
	i = 0;
	Skeleton3::edge_iterator E, Eend;
	for (boost::tie(E, Eend) = boost::edges(*S); E != Eend; ++E)
	{
		this->edges[i++] = (unsigned int)vertex_lookup[boost::source(*E, *S)];
		this->edges[i++] = (unsigned int)vertex_lookup[boost::target(*E, *S)];
	}
	
	// Save vertex data to the graphics card
	glGenBuffers(1, this->bufs);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble)*nverts*3, verts, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	// Cleanup temporary data
	delete[] verts;
}
GlSkeleton::GlSkeleton(const SkeletonGraph3* SG) : nverts(0), nedges(0), edges(nullptr)
{
	// Note: we end up duplicating all branch point points with the current method, but that doesn't really matter

	size_t i, j;
	this->bufs[0] = 0;

	// Get sizes of vertices and edges
	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), end = SG->branches_end(); B != end; ++B)
	{
		size_t s = B->size();
		this->nverts += s;
		this->nedges += s - 1;
	}

	// Get all vertices
	GLdouble* verts = new GLdouble[this->nverts*3];
	i = 0;
	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), end = SG->branches_end(); B != end; ++B)
	{
		for (SkeletonGraph3::Branch::const_iterator C = B->begin(), end = B->end(); C != end; ++C)
		{
			const Point3& v = *C;
			verts[i++] = CGAL::to_double(v.x());
			verts[i++] = CGAL::to_double(v.y());
			verts[i++] = CGAL::to_double(v.z());
		}
	}
	
	// Get all edges
	this->edges = new unsigned int[2*this->nedges];
	i = 0; j = 0;
	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), end = SG->branches_end(); B != end; ++B)
	{
		for (SkeletonGraph3::Branch::const_iterator C = B->begin(), end = std::prev(B->end()); C != end; ++C)
		{
			this->edges[i++] = (unsigned int)j++;
			this->edges[i++] = (unsigned int)j;
		}
		j++;
	}
	
	// Save vertex data to the graphics card
	glGenBuffers(1, this->bufs);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble)*nverts*3, verts, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	// Cleanup temporary data
	delete[] verts;
}
GlSkeleton::~GlSkeleton()
{
	if (this->bufs) { glDeleteBuffers(1, this->bufs); this->bufs[0] = 0; }
	if (this->edges) { delete[] this->edges; this->edges = nullptr; }
}
void GlSkeleton::render(const double weight, const QColor& color)
{
	// Setup the GL environment
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glLineWidth(weight);
	_gl_set_color(color);

	// Draw all edges
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, 0);
	glDrawElements(GL_LINES, (GLsizei)this->nedges*2, GL_UNSIGNED_INT, this->edges);
	glDisableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}


////////// GlIntersection //////////
GlIntersection::GlIntersection(const Intersection& I) : nintersections(0), npoints(nullptr)
{
	this->bufs[0] = 0;

	const Plane3& h = I.plane();

	// Calculate the plane corners
	// TODO: instead do an off-axis best-fit box or a circle.
	CGAL::Bbox_2 bb = I.bbox();
	Point2 center = Point2((bb.xmax()+bb.xmin())/2, (bb.ymax()+bb.ymin())/2);
	double x_radius = 0.75*(bb.xmax() - bb.xmin());
	double y_radius = 0.75*(bb.ymax() - bb.ymin());
	this->plane_corners[0] = h.to_3d(Point2(center.x() + x_radius, center.y() + y_radius));
	this->plane_corners[1] = h.to_3d(Point2(center.x() + x_radius, center.y() - y_radius));
	this->plane_corners[2] = h.to_3d(Point2(center.x() - x_radius, center.y() - y_radius));
	this->plane_corners[3] = h.to_3d(Point2(center.x() - x_radius, center.y() + y_radius));

	// Buffer the intersection lines
	this->nintersections = I.count();
	size_t pts_len = I.total_count()*3;
	GLdouble* pts = new GLdouble[pts_len];
	this->npoints = new size_t[this->nintersections];
	this->polys = new unsigned int*[this->nintersections];
	size_t i = 0, p = 0;
	for (Intersection::const_iterator IP2 = I.begin(), ip2_end = I.end(); IP2 != ip2_end; ++IP2)
	{
		size_t start_p = p / 3;
		for (Polygon2::Vertex_const_iterator V = IP2->vertices_begin(), end = IP2->vertices_end(); V != end; ++V)
		{
			const Point3& v = h.to_3d(*V);
			pts[p++] = CGAL::to_double(v.x());
			pts[p++] = CGAL::to_double(v.y());
			pts[p++] = CGAL::to_double(v.z());
		}
		size_t n = IP2->size()*2, j = 1, end_p = p / 3, last_p;
		if (IP2->is_open()) { n -= 2; end_p -= 1; last_p = end_p;   }
		else                {                     last_p = start_p; }
		this->npoints[i] = n;
		this->polys[i] = new unsigned int[n];
		this->polys[i][0] = (unsigned int)start_p;
		this->polys[i][n-1] = (unsigned int)last_p;
		for (size_t x = start_p + 1; x < end_p; ++x) { this->polys[i][j++] = (unsigned int)x; this->polys[i][j++] = (unsigned int)x; }
		i++;
	}
	
	// Save points data to the graphics card
	glGenBuffers(1, bufs);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble)*pts_len, pts, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Cleanup temporary data
	delete[] pts;
}
GlIntersection::~GlIntersection()
{
	if (this->bufs) { glDeleteBuffers(1, this->bufs); this->bufs[0] = 0; }
	if (this->npoints) { delete[] this->npoints; this->npoints = nullptr; }
	if (this->polys) { for (size_t i = 0; i < this->nintersections; ++i) { if (this->polys[i]) { delete[] this->polys[i]; this->polys[i] = nullptr; } } delete[] this->polys; this->polys = nullptr; }
}
void GlIntersection::render_plane(const QColor& color)
{
	// Setup the GL environment
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	_gl_set_color(color);

	// Draw the plane
	glBegin(GL_QUADS);
	_gl_add_point(this->plane_corners[0]);
	_gl_add_point(this->plane_corners[1]);
	_gl_add_point(this->plane_corners[2]);
	_gl_add_point(this->plane_corners[3]);
	glEnd();
}
void GlIntersection::render_polygons(const QColor& color)
{
	// Setup the GL environment
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glLineWidth(5.0f);
	_gl_set_color(color);
	
	// Draw all edges
	glBindBuffer(GL_ARRAY_BUFFER, this->bufs[0]);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, 0);
	for (size_t i = 0; i < this->nintersections; ++i)
	{
		glDrawElements(GL_LINES, (GLsizei)this->npoints[i], GL_UNSIGNED_INT, this->polys[i]);
	}
	glDisableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//// Indirect/slow/easy way
	//glBegin(GL_LINES);
	//Polygon2::Edge_const_circulator e = this->p.edges_circulator(), end = e;
	//CGAL_For_all(e, end)
	//{
	//	_gl_add_point(this->h.to_3d(e->source()));
	//	_gl_add_point(this->h.to_3d(e->target()));
	//}
	//glEnd();
}



////////// GlPoint //////////
GlPoint::GlPoint(const Point3& pt, double radius) : pt(pt), radius(radius)
{
}
GlPoint::~GlPoint()
{
}
#include <iostream>

void GlPoint::render(const QColor& color)
{
	// Setup the GL environment
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glPointSize(2*this->radius);
	_gl_set_color(color);

	// Draw the point
	glBegin(GL_POINTS);
	_gl_add_point(this->pt);
	glEnd();
}
