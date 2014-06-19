#include "GeometryTypes.hpp"

#include "IO.hpp"
#include "Intersection.hpp"
#include "MedialAxisTransform.hpp"
#include "MedialAxisTransform_IO.hpp"
#include "Skeleton.hpp"

#include <iostream>
#include <string>
#include <algorithm>
#include <ctype.h>

#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt4.h>

std::string str_tolower(std::string s)
{
	std::transform(s.begin(), s.end(), s.begin(), ::tolower);
	return s;
}

int main(int argc, char **argv)
{
	QApplication app(argc, argv);
	//app.setOrganizationDomain("inria.fr");
	//app.setOrganizationName("INRIA");
	app.setApplicationName("Cross Section");

	// Import resources from libCGALQt4.
	// See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
	CGAL_QT4_INIT_RESOURCES;

	CGAL::set_pretty_mode(std::cout);
	
	///////////////////////////////////////////////////////////////////////////
	// Set the file to read
	///////////////////////////////////////////////////////////////////////////
	bool assume_good = true;
	//std::string filename = "example-meshes/cube-triangles.off";
	//std::string filename = "example-meshes/cube-quads.off";
	//std::string filename = "small.off";
	//std::string filename = "big.off";
	std::string filename = "big_simplified.off";
	std::string obj_name = "";

	///////////////////////////////////////////////////////////////////////////
	// Read in the 3D polyhedral mesh
	///////////////////////////////////////////////////////////////////////////
	Polyhedron3 *P;
	size_t dot = filename.find_last_of('.');
	std::string ext = dot == std::string::npos ? "" : str_tolower(filename.substr(dot));
	if (ext == ".obj")
	{
		ObjFile obj(filename.c_str(), false, assume_good);
		obj_name = obj.names().at(0);
		P = obj[obj_name];
	}
	else if (ext == ".off")
	{
		P = read_off_file(filename.c_str(), assume_good);
	}
	else
	{
		throw std::invalid_argument("unknown file format");
	}
	// TODO: delete P; somewhere...

	///////////////////////////////////////////////////////////////////////////
	// Construct the medial axis transform from the mesh
	///////////////////////////////////////////////////////////////////////////
	MAT* mat = construct_medial_axis_transform(P, filename, obj_name);
	//dump_mat(mat);
	//MAT* mat = load_mat("output.mat", P);
	// TODO: delete mat; somewhere...

	///////////////////////////////////////////////////////////////////////////
	// Construct the skeleton from the medial axis transform
	///////////////////////////////////////////////////////////////////////////
	Skeleton3* S = construct_skeleton(mat, 0.7, 0.0);
	//Skeleton3* S = read_cg("skeleton-small.cg");
	// TODO: delete S; somewhere...

	///////////////////////////////////////////////////////////////////////////
	// Get the "graph" version of the skeleton - easier to use
	///////////////////////////////////////////////////////////////////////////
	SkeletonGraph3* SG = construct_skeleton_graph(S); // TODO: delete SG; somewhere...
	skeleton_reduce(SG, 25.0);
	
	///////////////////////////////////////////////////////////////////////////
	// Pick a random plane to use with intersections.
	// It is through a skeleton edge midpoint, perpendicular to that edge.
	///////////////////////////////////////////////////////////////////////////
	Plane3 h;
	srand(time(nullptr));
	do
	{
		SkeletonGraph3::Branch_handle B = SG->branches_begin() + (rand() % SG->size_of_branches());
		SkeletonGraph3::Branch::iterator i = B->begin();
		std::advance(i, rand() % B->size());
		Point3 p = *i++, q = *i;
		Vector3 v = p-q;
		h = Plane3(p+v/2, v);
	}
	while (h.is_degenerate());
	std::cout << "Chosen plane: " << h.point() << " " << h.orthogonal_vector() << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Calculate the intersection of the plane with the polyhedron
	///////////////////////////////////////////////////////////////////////////
#ifdef CGAL_OLD_FACET_TREE
	FacetTree t(P->facets_begin(), P->facets_end());
#else
	FacetTree t(P->facets_begin(), P->facets_end(), *P);
#endif
	Intersection intersection(t, h);
	std::cout << "Intersections: " << intersection.count() << ", Total area: " << intersection.area() << std::endl;
	for (size_t i = 0; i < intersection.count(); ++i)
	{
		std::cout << "  Area: " << intersection[i].area() << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////
	// Setup the viewer
	///////////////////////////////////////////////////////////////////////////
	Viewer viewer;
	viewer.setWindowTitle("Cross Section");
	viewer.set_polyhedron(P);
	viewer.set_skeleton(SG);
	viewer.set_intersection(intersection);
	viewer.show();

	///////////////////////////////////////////////////////////////////////////
	// Run the visual part of the program
	///////////////////////////////////////////////////////////////////////////
	QStringList args = app.arguments();
	args.removeAt(0);
	if (!args.empty() && args[0] == "--use-meta")
	{
		viewer.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
		args.removeAt(0);
	}
	//Q_FOREACH(QString filename, args)
	//	viewer.open(filename);
	return app.exec();
}

