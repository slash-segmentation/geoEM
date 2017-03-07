#include "GeometryTypes.hpp"

#include "IO.hpp"
#include "Intersection.hpp"
#include "MedialAxisTransform.hpp"
#include "MedialAxisTransform_IO.hpp"
#include "Skeleton.hpp"

#include <boost/timer/timer.hpp>

#include <iostream>

#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt.h>


int main(int argc, char **argv)
{
	QApplication app(argc, argv);
	//app.setOrganizationDomain("inria.fr");
	//app.setOrganizationName("INRIA");
	app.setApplicationName("Cross Section");

	// Import resources from libCGALQt4.
	// See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
	CGAL_QT_INIT_RESOURCES;

	CGAL::set_pretty_mode(std::cout);
	
	///////////////////////////////////////////////////////////////////////////
	// Set the file to read
	///////////////////////////////////////////////////////////////////////////
	bool assume_good = true;
	//std::string filename = "example-data/other/cube-triangles.off";
	//std::string filename = "example-data/other/cube-quads.off";
	std::string filename = "example-data/other/elephant.off";
	//std::string filename = "example-data/small.off";
	//std::string filename = "example-data/big.off";
	//std::string filename = "example-data/big_simplified.off";
	std::string obj_name = "";

	///////////////////////////////////////////////////////////////////////////
	// Read in the 3D polyhedral mesh
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Reading mesh..." << std::endl;
	Polyhedron3 *P;
	{
		boost::timer::auto_cpu_timer t;
		P = read_mesh(filename, assume_good); // TODO: delete P; somewhere...
	}
	std::cout << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Construct the medial axis transform from the mesh
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Constructing medial axis transform..." << std::endl;
	MAT* mat;
	{
		boost::timer::auto_cpu_timer t;
		mat = construct_medial_axis_transform(P, filename, obj_name);
		//dump_mat(mat);
		//MAT* mat = load_mat("output.mat", P);
	}
	std::cout << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Construct the skeleton from the medial axis transform
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Constructing skeleton..." << std::endl;
	Skeleton3* S;
	{
		boost::timer::auto_cpu_timer t;
		S = construct_skeleton(mat, 0.7, 0.0); // TODO: how to determine parameters?
		//Skeleton3* S = read_cg("skeleton-small.cg");
	}
	std::cout << std::endl;
	delete mat; // done with the medial axis transform

	///////////////////////////////////////////////////////////////////////////
	// Get the "graph" version of the skeleton - easier to use
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Constructing skeleton graph..." << std::endl;
	SkeletonGraph3* SG;
	{
		boost::timer::auto_cpu_timer t;
		SG = construct_skeleton_graph(S); // TODO: delete SG; somewhere...
		skeleton_reduce(SG);
	}
	std::cout << std::endl;
	delete S; // done with the skeleton
	
	///////////////////////////////////////////////////////////////////////////
	// Pick a random plane to use with intersections.
	// It is through a skeleton edge midpoint, perpendicular to that edge.
	///////////////////////////////////////////////////////////////////////////
	Plane3 h;
	Point3 p;
	srand(time(nullptr));
	do
	{
		SkeletonGraph3::Branch_handle B = SG->branches_begin() + (rand() % SG->size_of_branches());
		SkeletonGraph3::Branch::iterator i = B->begin();
		std::advance(i, rand() % B->size());
		p = *i++;
		Point3 q = *i;
		Vector3 v = p-q;
		h = Plane3(p+v/2, v);
	}
	while (h.is_degenerate());
	std::cout << "Chosen plane: " << h.point() << " " << h.orthogonal_vector() << std::endl;
	std::cout << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Calculate the intersection of the plane with the polyhedron
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Calculating intersection..." << std::endl;
	Intersection intersection;
	{
		boost::timer::auto_cpu_timer t_;
		FacetTree t(P->facets_begin(), P->facets_end(), *P);
		intersection = Intersection(t, h);
		std::cout << "Intersections: " << intersection.count() << ", Total area: " << intersection.area() << std::endl;
		for (size_t i = 0; i < intersection.count(); ++i)
		{
			std::cout << "  Area: " << intersection[i].area() << std::endl;
		}
	}
	std::cout << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Find which intersection is actually around the point of interest
	///////////////////////////////////////////////////////////////////////////
	Point2 p2 = to_2d(p, h);
	std::cout << h << " " << p2 << std::endl;
	Intersection intersection_isolated;
	for (size_t i = 0; i < intersection.count(); ++i)
	{
		if (intersection[i].has_on_bounded_side(p2))
		{
			intersection_isolated.add(intersection[i]);
			std::cout << "The intersection #" << i << " is the correct one" << std::endl;
			std::cout << "It has an area of " << intersection[i].area() << std::endl;
			std::cout << "And the points:";
			for (auto itr = intersection[i].vertices_begin(), end = intersection[i].vertices_end(); itr != end; ++itr)
			{
				std::cout << *itr << ", ";
			}
			std::cout << std::endl;
		}
	}

	//return 0;

	/////////////////////////////////////////////////////////////////////////////
	//// Setup the viewer
	/////////////////////////////////////////////////////////////////////////////
	Viewer viewer;
	viewer.setWindowTitle("Cross Section");
	viewer.set_polyhedron(P);
	viewer.set_skeleton(SG);
	viewer.set_intersection(intersection);
	viewer.show();

	/////////////////////////////////////////////////////////////////////////////
	//// Run the visual part of the program
	/////////////////////////////////////////////////////////////////////////////
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

