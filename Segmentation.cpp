#include "Segmentation.hpp"

#include "GeometryUtils.hpp"

#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

typedef Polyhedron3::HalfedgeDS P3_HDS;
typedef Polyhedron3::Vertex_const_handle P3_Vertex;
typedef Polyhedron3::Facet_const_handle P3_Facet;

typedef std::vector<P3_Vertex> Vertices;
typedef std::vector<P3_Facet> Facets;

typedef handle_map<P3_Vertex, double> Vertex_double_map;
typedef handle_map<P3_Vertex, size_t> Vertex_int_map;
typedef handle_map<P3_Facet, double> Facet_double_map;
typedef handle_map<P3_Facet, size_t> Facet_int_map;

typedef boost::associative_property_map<Facet_double_map> Facet_double_prop_map;
typedef boost::associative_property_map<Facet_int_map> Facet_int_prop_map;

// Builds a Polyhedron3 from a list of vertices and facets and, a reverse lookup map of vertices.
class Builder : public CGAL::Modifier_base<P3_HDS>
{
    typedef typename P3_HDS::Vertex::Point Point;
	const Vertices& verts;
	Vertex_int_map& vert_lookup;
	const Facets& faces;
public:
    Builder(const Vertices& v, Vertex_int_map& vl, const Facets& f) :
		verts(v), vert_lookup(vl), faces(f) { }
    void operator()(P3_HDS& hds)
	{
        CGAL::Polyhedron_incremental_builder_3<P3_HDS> B(hds, true);
        B.begin_surface(verts.size(), faces.size());
		for (Vertices::const_iterator itr = verts.begin(), end = verts.end(); itr != end; ++itr)
		{
			B.add_vertex((*itr)->point());
		}
		for (Facets::const_iterator F = faces.begin(), Fend = faces.end(); F != Fend; ++F)
		{
			B.begin_facet();
			Polyhedron3::Halfedge_around_facet_const_circulator H = (*F)->facet_begin(), Hstart = H;
			do
			{
				B.add_vertex_to_facet(vert_lookup[H->vertex()]);
			} while (++H != Hstart);
			B.end_facet();
		}
        B.end_surface();
    }
};

std::vector<Polyhedron3*> get_meshes(const Polyhedron3* P, const Facet_int_prop_map segs, const size_t n)
{
	// Get a list of faces/verts for each segmented portion of the polyhedron
	std::vector<Facets> faces(n);
	std::vector<Vertices> verts(n);
	std::vector<Vertex_int_map> vert_lookup(n);
    for (P3_Facet F = P->facets_begin(), end = P->facets_end(); F != end; ++F)
	{
		size_t i = segs[F];
		faces[i].push_back(F);

		Vertices& vs = verts[i];
		Vertex_int_map& vls = vert_lookup[i];
		Polyhedron3::Halfedge_around_facet_const_circulator H = F->facet_begin(), Hstart = H;
		do
		{
			if (vls.insert(std::make_pair(H->vertex(), vs.size())).second)
			{
				vs.push_back(H->vertex());
			}
		} while (++H != Hstart);
	}
	
	// Now we can build the polyhedrons themselves
	std::vector<Polyhedron3*> p3s(n);
    for (size_t i = 0; i < n; ++i)
	{
		p3s[i] = new Polyhedron3();
		Builder builder(verts[i], vert_lookup[i], faces[i]);
		p3s[i]->delegate(builder);
	}
	return p3s;
}

std::vector<Polyhedron3*> compute_segmentation(const Polyhedron3* P, double cone_angle, size_t number_of_rays, size_t number_of_clusters, double smoothing_lambda)
{
	// Property-map for segment-ids
	Facet_int_map seg_map;
	Facet_int_prop_map seg_prop_map(seg_map);

	// Calculate SDF values and segment the mesh
	size_t num_segs = CGAL::segmentation_via_sdf_values(*P, seg_prop_map, cone_angle, number_of_rays, number_of_clusters, smoothing_lambda);

	// Break the original mesh into the segments
	return get_meshes(P, seg_prop_map, num_segs);
}

std::vector<Polyhedron3*> compute_segmentation(const Polyhedron3* P, const Skeleton3* S, size_t number_of_clusters, double smoothing_lambda)
{
	// For each input vertex compute its distance to the skeleton
	Vertex_double_map distances(P->size_of_vertices());
	Skeleton3::vertex_iterator V, Vend;
	for (boost::tie(V, Vend) = boost::vertices(*S); V != Vend; ++V)
	{
		const Point3& skel_pt = (*S)[*V].point;
		const std::vector<Polyhedron3::Vertex_handle>& verts = (*S)[*V].vertices;
		for (std::vector<Polyhedron3::Vertex_handle>::const_iterator itr = verts.begin(), end = verts.end(); itr != end; ++itr)
		{
			const Point3& mesh_pt = (*itr)->point();
			distances[*itr] = distance(skel_pt, mesh_pt);
		}
	}

	// Compute SDF values with skeleton
	Facet_double_map sdf_map;
	Facet_double_prop_map sdf_prop_map(sdf_map);
	for (P3_Facet F = P->facets_begin(), Fend = P->facets_end(); F != Fend; ++F)
	{
		double dist = 0;
		Polyhedron3::Halfedge_around_facet_const_circulator H = F->facet_begin(), Hstart = H;
		do { dist += distances[H->vertex()]; } while (++H != Hstart);
		sdf_prop_map[F] = dist / 3.0;
	}

	// Post-process the SDF values
	CGAL::sdf_values_postprocessing(*P, sdf_prop_map);

	// Segment the mesh
	Facet_int_map seg_map;
	Facet_int_prop_map seg_prop_map(seg_map);
	size_t num_segs = CGAL::segmentation_from_sdf_values(*P, sdf_prop_map, seg_prop_map, number_of_clusters, smoothing_lambda);
	
	// Break the original mesh into the segments
	return get_meshes(P, seg_prop_map, num_segs);
}
