#include "PolygonSorter.hpp"

#include <CGAL/Cartesian.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Incremental_neighbor_search.h>
//#include <CGAL/Orthogonal_incremental_neighbor_search.h>

#include "GeometryUtils.hpp"

////static Kernel::FT distance(const Point3& p, const Plane3& h)
////{
////	//Kernel::FT top = Vector3(p.x(), p.y(), p.z())*h.orthogonal_vector();
////	//return top*top/(h.orthogonal_vector()*h.orthogonal_vector())
////	Kernel::FT top = h.a()*p.x()+h.b()*p.y()+h.c()*p.z();
////	return top*top/(h.a()*h.a()+h.b()*h.b()+h.c()*h.c());
////}
//static Kernel::FT rel_distance(const Point3& p, const Plane3& h) { return h.a()*p.x()+h.b()*p.y()+h.c()*p.z(); } // a distance that is comparable for anything using the same plane

///// Base Types for Searching (including a custom distance class) /////
typedef CGAL::Cartesian<double> DoubleKernel; // we are using a simpler, faster, kernel since nothing here really matters
typedef DoubleKernel::Point_3 Point;
typedef DoubleKernel::Plane_3 Plane;
typedef DoubleKernel::Vector_3 Vector;
typedef boost::tuple<Point, unsigned int, unsigned int, unsigned int> PointAndFaceVertices;
class Plane_Euclidean_Distance
{
	// Distance from point to plane is (a*x+b*y+c*z+d)/sqrt(a*a+b*b+c*c)
	// The transformed distances will be (a*x+b*y+c*z)
public:
	typedef DoubleKernel::FT      FT;
	typedef CGAL::Dimension_tag<3> D;
	typedef PointAndFaceVertices  Point_d;
	typedef DoubleKernel::Plane_3 Query_item;
	typedef DoubleKernel::Point_3 Point3;
private:
	const Query_item h;
public:
	Plane_Euclidean_Distance(const Query_item& h):h(h) {} // the distance calculator is fixed to the plane that will be used to query
	
	inline FT transformed_distance(double x, double y, double z) const { return h.a()*x+h.b()*y+h.c()*z; }
	inline FT transformed_distance(const Query_item& h, const Point_d& pf) const { assert(h==this->h); const Point3& p = pf.get_head(); return transformed_distance(p.x(), p.y(), p.z()); }
	inline FT min_distance_to_rectangle(const Query_item& h, const CGAL::Kd_tree_rectangle<FT, D>& r) const
	{
		assert(h==this->h);
		// TODO: more efficient?
		const FT dist1 = this->transformed_distance(r.min_coord(0), r.min_coord(1), r.min_coord(2));
		const FT dist2 = this->transformed_distance(r.min_coord(0), r.min_coord(1), r.max_coord(2));
		const FT dist3 = this->transformed_distance(r.min_coord(0), r.max_coord(1), r.min_coord(2));
		const FT dist4 = this->transformed_distance(r.min_coord(0), r.max_coord(1), r.max_coord(2));
		const FT dist5 = this->transformed_distance(r.max_coord(0), r.min_coord(1), r.min_coord(2));
		const FT dist6 = this->transformed_distance(r.max_coord(0), r.min_coord(1), r.max_coord(2));
		const FT dist7 = this->transformed_distance(r.max_coord(0), r.max_coord(1), r.min_coord(2));
		const FT dist8 = this->transformed_distance(r.max_coord(0), r.max_coord(1), r.max_coord(2));
		FT min_dist = dist1;
		if (dist2 < min_dist) { min_dist = dist2; }
		if (dist3 < min_dist) { min_dist = dist3; }
		if (dist4 < min_dist) { min_dist = dist4; }
		if (dist5 < min_dist) { min_dist = dist5; }
		if (dist6 < min_dist) { min_dist = dist6; }
		if (dist7 < min_dist) { min_dist = dist7; }
		if (dist8 < min_dist) { min_dist = dist8; }
		return min_dist;
		
		//FT distance = FT(0);

		//if      (q[0] < r.min_coord(0))	distance += (r.min_coord(0)-q[0])*(r.min_coord(0)-q[0]);
		//else if (q[0] > r.max_coord(0))	distance += (q[0]-r.max_coord(0))*(q[0]-r.max_coord(0));

		//if      (q[1] < r.min_coord(1))	distance += (r.min_coord(1)-q[1])*(r.min_coord(1)-q[1]);
		//else if (q[1] > r.max_coord(1))	distance += (q[1]-r.max_coord(1))*(q[1]-r.max_coord(1));

		//if      (q[2] < r.min_coord(2))	distance += (r.min_coord(2)-q[2])*(r.min_coord(2)-q[2]);
		//else if (q[2] > r.max_coord(2))	distance += (q[2]-r.max_coord(2))*(q[2]-r.max_coord(2));
	}
	inline FT max_distance_to_rectangle(const Query_item& h, const CGAL::Kd_tree_rectangle<FT, D>& r) const
	{
		assert(h==this->h);
		// TODO: more efficient?
		const FT dist1 = this->transformed_distance(r.min_coord(0), r.min_coord(1), r.min_coord(2));
		const FT dist2 = this->transformed_distance(r.min_coord(0), r.min_coord(1), r.max_coord(2));
		const FT dist3 = this->transformed_distance(r.min_coord(0), r.max_coord(1), r.min_coord(2));
		const FT dist4 = this->transformed_distance(r.min_coord(0), r.max_coord(1), r.max_coord(2));
		const FT dist5 = this->transformed_distance(r.max_coord(0), r.min_coord(1), r.min_coord(2));
		const FT dist6 = this->transformed_distance(r.max_coord(0), r.min_coord(1), r.max_coord(2));
		const FT dist7 = this->transformed_distance(r.max_coord(0), r.max_coord(1), r.min_coord(2));
		const FT dist8 = this->transformed_distance(r.max_coord(0), r.max_coord(1), r.max_coord(2));
		FT max_dist = dist1;
		if (dist2 > max_dist) { max_dist = dist2; }
		if (dist3 > max_dist) { max_dist = dist3; }
		if (dist4 > max_dist) { max_dist = dist4; }
		if (dist5 > max_dist) { max_dist = dist5; }
		if (dist6 > max_dist) { max_dist = dist6; }
		if (dist7 > max_dist) { max_dist = dist7; }
		if (dist8 > max_dist) { max_dist = dist8; }
		return max_dist;

		//FT distance = FT(0);

		//if (q[0] <= (r.min_coord(0)+r.max_coord(0))/FT(2.0)) distance += (r.max_coord(0)-q[0])*(r.max_coord(0)-q[0]);
		//else                                                 distance += (q[0]-r.min_coord(0))*(q[0]-r.min_coord(0));

		//if (q[1] <= (r.min_coord(1)+r.max_coord(1))/FT(2.0)) distance += (r.max_coord(1)-q[1])*(r.max_coord(1)-q[1]);
		//else                                                 distance += (q[1]-r.min_coord(1))*(q[1]-r.min_coord(1));

		//if (q[2] <= (r.min_coord(2)+r.max_coord(2))/FT(2.0)) distance += (r.max_coord(2)-q[2])*(r.max_coord(2)-q[2]);
		//else                                                 distance += (q[2]-r.min_coord(2))*(q[2]-r.min_coord(2));
	}
	//inline FT new_distance(FT dist, FT old_off, FT new_off, int /* cutting_dimension */) const // TODO: needed for Orthogonal NN searches (which should be must faster)
	//{
	//	FT new_dist = dist + (new_off*new_off - old_off*old_off);
	//	return new_dist;
	//}
	inline FT transformed_distance(FT d) const { return d*ft_sqrt(h.a()*h.a()+h.b()*h.b()+h.c()*h.c()) - h.d(); }
	inline FT inverse_of_transformed_distance(FT d) const { return (d+h.d())/ft_sqrt(h.a()*h.a()+h.b()*h.b()+h.c()*h.c()); }
};
typedef CGAL::Search_traits_adapter<PointAndFaceVertices, CGAL::Nth_of_tuple_property_map<0, PointAndFaceVertices>, CGAL::Search_traits_3<DoubleKernel>> NNSearchTraits;
typedef CGAL::/*Orthogonal_*/Incremental_neighbor_search<NNSearchTraits, Plane_Euclidean_Distance> NNIncrementalSearch;


///// Iterator Transformer /////
struct CalculateTriangleCentriod
{
	const double* vs;
	inline CalculateTriangleCentriod(const double* vs) : vs(vs) { }
	inline const Point pt(unsigned int v) const { v*= 3; return Point(vs[v], vs[v+1], vs[v+2]); }
	inline const PointAndFaceVertices operator() (PolygonSorter::ConstFace f) const { return boost::make_tuple(CGAL::centroid(pt(f[0]), pt(f[1]), pt(f[2])), f[0], f[1], f[2]); }
};
typedef boost::transform_iterator<CalculateTriangleCentriod, PolygonSorter::ConstFaces, PointAndFaceVertices, PointAndFaceVertices> CalculateTriangleCentriodIter;


///// PolygonSorter /////
struct PolygonSorter::internal_data
{
	NNIncrementalSearch::Tree* tree;
	internal_data(NNIncrementalSearch::Tree* t = nullptr) : tree(t) {}
	~internal_data() { if (this->tree) { delete this->tree; this->tree = nullptr; } }
};
PolygonSorter::PolygonSorter(const double* vertices, ConstFaces faces, size_t faces_len) : data(new PolygonSorter::internal_data())
{
	CalculateTriangleCentriod ctc = CalculateTriangleCentriod(vertices);
	this->data->tree = new NNIncrementalSearch::Tree(CalculateTriangleCentriodIter(faces, ctc), CalculateTriangleCentriodIter(faces + faces_len, ctc));
}
PolygonSorter::~PolygonSorter() { if (this->data) { delete this->data; this->data = nullptr; } }
void PolygonSorter::query(const Point3& viewing_source, const Direction3& viewing_dir, unsigned int* buf)
{
	Plane h(Point(CGAL::to_double(viewing_source.x()), CGAL::to_double(viewing_source.y()), CGAL::to_double(viewing_source.z())), Vector(CGAL::to_double(viewing_dir.dx()), CGAL::to_double(viewing_dir.dy()), CGAL::to_double(viewing_dir.dz())));
	NNIncrementalSearch search(*this->data->tree, h, 0.1, false, Plane_Euclidean_Distance(h));
	size_t i = 0;
	for (NNIncrementalSearch::iterator si = search.begin(), end = search.end(); si != end; ++si)
	{
		buf[i++] = boost::get<1>(si->first);
		buf[i++] = boost::get<2>(si->first);
		buf[i++] = boost::get<3>(si->first);
	}
}
