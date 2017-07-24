#include "GeometryTypes.hpp"

#include "IO.hpp"
#include "IO_OBJ.hpp"
#include "Intersection.hpp"
#include "MedialAxisTransform.hpp"
#include "MedialAxisTransform_IO.hpp"
#include "Skeleton.hpp"
#include "Segmentation.hpp"

#include <boost/timer/timer.hpp>

#include <iostream>

#ifdef CREATE_GUI
#include "Viewer.hpp"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/auto_link/Qt.h>
#endif

#include <CGAL/boost/graph/split_graph_into_polylines.h>

struct Display_polylines {
	const Skeleton3& skeleton;
	std::ofstream& out;
	int polyline_size;
	std::stringstream sstr;
	Display_polylines(const Skeleton3& skeleton, std::ofstream& out)
		: skeleton(skeleton), out(out)
	{}
	void start_new_polyline() {
		polyline_size = 0;
		sstr.str("");
		sstr.clear();
	}
	void add_node(Skeleton3::vertex_descriptor v) {
		++polyline_size;
		sstr << " " << skeleton[v].point;
	}
	void end_polyline()
	{
		out << polyline_size << sstr.str() << "\n";
	}
};


Intersection get_intersection(const FacetTree& ft, Point3 pt, Vector3 n)
{
	Plane3 h(pt, n);
	Point2 pt2 = to_2d(pt, h);
	Intersection intersection_full(ft, h), intersection(h);
	for (size_t i = 0; i < intersection_full.count(); ++i)
	{
		if (intersection_full[i].has_on_bounded_side(pt2))
		{
			intersection.add(intersection_full[i]);
		}
	}
	return intersection;
}


int main(int argc, char **argv)
{
#ifdef CREATE_GUI
	QApplication app(argc, argv);
	//app.setOrganizationDomain("inria.fr");
	//app.setOrganizationName("INRIA");
	app.setApplicationName("Cross Section");

	// Import resources from libCGALQt4.
	// See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
	CGAL_QT_INIT_RESOURCES;
#endif

	CGAL::set_pretty_mode(std::cout);
	
	///////////////////////////////////////////////////////////////////////////
	// Set the file to read
	///////////////////////////////////////////////////////////////////////////
	bool assume_good = true;
	//std::string filename = "example-data/other/cube-triangles.off";
	//std::string filename = "example-data/other/cube-quads.off";
	//std::string filename = "example-data/other/elephant.off";
	//std::string filename = "example-data/small.off";
	//std::string filename = "example-data/big.off";
	std::string filename = "example-data/big_simplified.off";
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

	bool use_mat = false;

	Skeleton3* S;
	std::vector<Polyhedron3*> segments;
	if (use_mat)
	{
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
		std::cout << "Constructing skeleton from MAT..." << std::endl;
		{
			boost::timer::auto_cpu_timer t;
			S = construct_skeleton(mat, 0.7, 0.0); // TODO: how to determine parameters?
			//Skeleton3* S = read_cg("skeleton-small.cg");
		}
		std::cout << std::endl;
		delete mat; // done with the medial axis transform
		
		///////////////////////////////////////////////////////////////////////////
		// Segmention
		///////////////////////////////////////////////////////////////////////////
		std::cout << "Segmenting..." << std::endl;
		{
			boost::timer::auto_cpu_timer t;
			segments = compute_segmentation(P, 2.0/3.0*CGAL_PI, 25, 5, 0.26);
		}
		std::cout << std::endl;
	}
	else
	{
		///////////////////////////////////////////////////////////////////////////
		// Construct the skeleton using mean curvature flow
		///////////////////////////////////////////////////////////////////////////
		std::cout << "Constructing skeleton using mean curvature flow..." << std::endl;
		{
			boost::timer::auto_cpu_timer t;
			S = construct_skeleton(P, 1); // TODO: how to determine parameter?
			//Skeleton3* S = read_cg("skeleton-small.cg");
		}
		std::cout << std::endl;

		///////////////////////////////////////////////////////////////////////////
		// Segmention
		///////////////////////////////////////////////////////////////////////////
		std::cout << "Segmenting using skeleton..." << std::endl;
		{
			boost::timer::auto_cpu_timer t;
			segments = compute_segmentation(P, S, 5, 0.26);
		}
		std::cout << std::endl;
	}
	
	std::ofstream output("skel.cgal");
	Display_polylines display(*S, output);
	CGAL::split_graph_into_polylines(*S, display);
	output.close();


	///////////////////////////////////////////////////////////////////////////
	// Get the "graph" version of the skeleton - easier to use
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Constructing skeleton graph..." << std::endl;
	SkeletonGraph3* SG;
	{
		boost::timer::auto_cpu_timer t;
		SG = construct_skeleton_graph(S); // TODO: delete SG; somewhere...
		//size_t nverts_before = SG->total_vertices();
		//skeleton_reduce(SG);
		//size_t nverts_after = SG->total_vertices();
		//std::cerr << "Removed " << nverts_before - nverts_after << " vertices/edges out of " << nverts_before << std::endl;
	}
	std::cout << std::endl;
	delete S; // done with the skeleton
	

	
	FacetTree ft(P->facets_begin(), P->facets_end(), *P);
	double total_length = 0.0;
	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), Bend = SG->branches_end(); B != Bend; ++B)
	{
		total_length += B->length();
	}
	size_t increments = 4096;
	double inc_length = total_length / increments;
	std::cout << total_length << ", " << increments << ", " << inc_length << std::endl;
	std::ofstream out("output.csv");
	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), Bend = SG->branches_end(); B != Bend; ++B)
	{
		const Point3& a = B->source()->coord();
		const Point3& b = B->target()->coord();
		out << "Branch," << a << "," << b << std::endl;
		SkeletonGraph3::Branch::const_iterator i = B->begin(), end = B->end();
		Point3 pt_end = *i++;
		double dist_carryover = 0.0;
		while (i != end)
		{
			Point3 pt = pt_end;
			pt_end = *i++;
			Vector3 n = normalized(pt_end - pt);
			double len_line = distance(pt, pt_end), dist_line = dist_carryover;
			const char* sep = "";
			while (dist_line < len_line)
			{
				Intersection intersection = get_intersection(ft, pt, n);
				out << sep << intersection.area();
				sep = ",";
				pt = pt + inc_length * n;
				dist_line += inc_length;
			}
			dist_carryover = dist_line - len_line;
		}
		out << std::endl;
	}
	out.close();

	write_obj_file("output.obj", segments, false);


	///////////////////////////////////////////////////////////////////////////
	// Calculate the intersection of a random plane which is perpendicular to
	// a point along the skeleton with the polyhedron
	///////////////////////////////////////////////////////////////////////////
	std::cout << "Calculating intersection..." << std::endl;
	Point3 pt;
	Intersection intersection;
	{
		boost::timer::auto_cpu_timer t_;

		srand(time(nullptr));
		SkeletonGraph3::Branch_handle B = SG->branches_begin() + (rand() % SG->size_of_branches());
		SkeletonGraph3::Branch::iterator itr = B->begin();
		std::advance(itr, rand() % (B->size() - 1));
		Point3 p = *itr++, q = *itr;
		Vector3 n = q - p;
		pt = p + n / 2;

		FacetTree ft(P->facets_begin(), P->facets_end(), *P);
		intersection = get_intersection(ft, pt, n);

		std::cout << "Intersections: " << intersection.count() << ", Total area: " << intersection.area() << std::endl;
		for (size_t i = 0; i < intersection.count(); ++i)
		{
			std::cout << "  Area: " << intersection[i].area() << std::endl;
			std::cout << "  Polygon:";
			for (auto itr = intersection[i].vertices_begin(), end = intersection[i].vertices_end(); itr != end; ++itr)
			{
				std::cout << *itr << ", ";
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;

#ifdef CREATE_GUI
	/////////////////////////////////////////////////////////////////////////////
	//// Setup the viewer
	/////////////////////////////////////////////////////////////////////////////
	Viewer viewer;
	viewer.setWindowTitle("Cross Section");
	viewer.set_polyhedron(P);
	viewer.set_skeleton(SG);
	viewer.set_intersection(intersection);
	viewer.set_point(pt, 100);
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
#else
	return 0;
#endif
}

