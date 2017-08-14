#include "Skeleton.hpp"

#include <queue>
#include <vector>
#include <utility>
#include "Heap.hpp"
#include "GeometryUtils.hpp"

///////////////////////////////////////////////////////////////////////////////////////
// Computes the "stable" skeleton from the medial axis transformation which is the set
// of the most obvious edges that will be part of the skeleton. They are marked such
// that there is at most one stable skeleton edge per facet.
///////////////////////////////////////////////////////////////////////////////////////
void compute_stable_skeleton(MAT* mat, double e_flux_thd, double wt_thd)
{
	static const MAT::Facet::Flags NSS_SS = MAT::Facet::Flags::NotStableSkeleton | MAT::Facet::Flags::StableSkeleton;

	// Mark facets that cannot have have any stable edges
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		if (f->weight() < wt_thd || (f->flags & MAT::Facet::Flags::NotStableSkeleton)) { continue; }
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end)
		{
			if (e->flux() > 0) { continue; }
			MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
			for (; ei != e && (ei->facet()->weight() < wt_thd || ei->flux() > 0); ++ei);
			if (ei != e)
			{
				f->flags           |= MAT::Facet::Flags::NotStableSkeleton;
				ei->facet()->flags |= MAT::Facet::Flags::NotStableSkeleton;
				break;
			}
		}
	}
	
	// Mark the stable edges
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		if (f->weight() < wt_thd || (f->flags & NSS_SS)) { continue; }
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end)
		{
			if (e->flux() <= e_flux_thd || e->next_identifed() == e) { continue; }
			MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
			for (; ei != e && (ei->facet()->weight() < wt_thd || (!(ei->facet()->flags & NSS_SS) && ei->flux() > e_flux_thd)); ++ei);
			
			// If we get here then along e there is not already a stable facet
			if (ei == e)
			{
				MAT::Edge::Edge_identified_circulator ei = e->circ_identified(), ei_end = ei;
				CGAL_For_all(ei, ei_end)
				{
					MAT::Facet_handle fa = ei->facet();
					if (fa->weight() < wt_thd) { continue; }
					fa->flags |= MAT::Facet::Flags::StableSkeleton;
					ei->flags |= MAT::Edge::Flags::StableSkeleton;
				}
				break;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// Computes the skeleton by retraction (as opposed to eating or "plain"). This works
// by continually removing boundary facets based on some criteria.
// Calls compute_stable_skeleton after reseting flags.
///////////////////////////////////////////////////////////////////////////////////////
void compute_skeleton(MAT* mat, double e_flux_thd, double wt_thd)
{
	// Initialize the flags
	static const MAT::Facet::Flags inv_mask = ~(MAT::Facet::Flags::Boundary | MAT::Facet::Flags::StableSkeleton | MAT::Facet::Flags::Eaten | MAT::Facet::Flags::NotStableSkeleton);
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		f->flags &= inv_mask;
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end) { e->flags = MAT::Edge::Flags::None; }
	}
	
	//-----------------------------------------------------
	//Compute the stable part of the skeleton
	compute_stable_skeleton(mat, e_flux_thd, wt_thd);
	//-----------------------------------------------------

	//-------------------------------------------------------------
	// Generate the boundary information. Ignore facet with weights less than "wt_thd".
	// Set the boundary flag and add the boundary facet into the queue.
	heap<double, MAT::Facet_handle> facet_heap;
	facet_heap.reserve(mat->size_of_facets());
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		if (f->weight() < wt_thd) { continue; }
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end)
		{
			MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
			for (; ei != e && ei->facet()->weight() <= wt_thd; ++ei);
			if (ei == e) { f->flags |= MAT::Facet::Flags::Boundary; break; }
		}
		if (f->flags & MAT::Facet::Flags::Boundary)
		{
			facet_heap.push_raw(f, f->weight());
			f->flags |= MAT::Facet::Flags::InHeap;
		}
	}
	facet_heap.heapify();
	//-------------------------------------------------------------

	//---------------------------------------------------------------
	// Start to retract.
	typedef std::vector<MAT::Edge_handle> Edges;
	std::vector<Edges> boundaries;
	//size_t fcount = 0;
	while (!facet_heap.empty())
	{
		MAT::Facet_handle f = facet_heap.top(); facet_heap.pop();
		//++fcount;
		f->flags |= MAT::Facet::Flags::Eaten;

		// Calculate the boundary information
		size_t degree = f->degree();
		if (boundaries.size() <= degree) { boundaries.resize(degree + 1); }
		boundaries[0].clear();
		std::vector<Edges>::iterator boundaries_end = boundaries.begin();
		{
			bool connect_last = false, connect_first = false;
			MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
			CGAL_For_all(e, e_end)
			{
				MAT::Vertex::Edge_iterator ev = e->vertex()->edges_begin(), ev_end = e->vertex()->edges_end();
				for (MAT::Facet_handle fa; ev != ev_end && ((fa = ev->facet()) == f || fa->weight() < wt_thd || ((fa->flags & MAT::Facet::Flags::Eaten) && !((ev->flags | ev->next_around_facet()->flags) & MAT::Edge::Flags::Skeleton))); ++ev);
				bool connect = ev != ev_end; // true if this edge is connected, false if it is a boundary
				if (e == e_end) { connect_first = connect; } else if (connect_last && !connect) { (++boundaries_end)->clear(); }
				connect_last = connect;

				MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
				for (MAT::Facet_handle fa; ei != e && ((fa = ei->facet())->weight() < wt_thd || (fa->flags & MAT::Facet::Flags::Eaten)); ++ei);
				connect = ei != e;
				if (!connect) { if (connect_last) { (++boundaries_end)->clear(); } boundaries_end->push_back(e); }
				connect_last = connect;
			}
			if (!connect_first && connect_last) { (++boundaries_end)->clear(); }
		}

		// If there is one facet component
		if (boundaries_end == boundaries.begin()) { std::cerr << facet_number(f, mat) << " is a one facet component." << std::endl; continue; }

		// Connect the possibly disconnected boundaries
		std::reverse(boundaries[0].begin(), boundaries[0].end());
		boundaries[0].insert(boundaries[0].end(), boundaries_end->rbegin(), boundaries_end->rend());
		boundaries_end->clear();

		// If there is only one connected boundary component, set the Skeleton flag only if there is a stable skeleton edge.
		if (boundaries_end == std::next(boundaries.begin()))
		{
			size_t i = 0, count = boundaries[0].size();
			for (; i != count && !(boundaries[0][i]->flags & MAT::Edge::Flags::StableSkeleton); ++i);
			if (i != count)
			{
				// i is the position of the stable skeleton edge, we set the Skeleton flag for the smaller of the ranges [0,i] or [i,count)
				size_t first = i, second = i;
				if (2*i <= count) { first = 0; } else { second = count - 1; }
				for (size_t i = first; i <= second; ++i)
				{
					MAT::Edge_handle e = boundaries[0][i];
					e->flags |= MAT::Edge::Flags::Skeleton;
					MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
					for (; ei != e; ++ei) { if (ei->facet()->weight() >= wt_thd) { ei->flags |= MAT::Edge::Flags::Skeleton; } }
				}
			}
		}
		else
		{
			// Select a connected boundary component which is not on the stable skeleton
			std::vector<Edges>::const_iterator notske;
			double min_inflow = DBL_MAX;
			for (std::vector<Edges>::iterator l = boundaries.begin(); l != boundaries_end; ++l)
			{
				double inflow = -DBL_MAX;
				for (Edges::iterator e = l->begin(), end = l->end(); e != end; ++e)
				{
					if ((*e)->flags & MAT::Edge::Flags::StableSkeleton) { inflow = DBL_MAX; break; } // always keep the component with stable one
					else if (inflow < (*e)->flux()) { inflow = (*e)->flux(); }
				}
				if (inflow < min_inflow) { notske = l; min_inflow = inflow; }
			}
			for (std::vector<Edges>::iterator l = boundaries.begin(); l != boundaries_end; ++l)
			{
				if (l == notske) { continue; }
				for (Edges::iterator e = l->begin(), end = l->end(); e != end; ++e)
				{
					(*e)->flags |= MAT::Edge::Flags::Skeleton;
					MAT::Edge::Edge_identified_circulator ei = (*e)->circ_identified();
					for (; ei != (*e); ++ei) { if (ei->facet()->weight() >= wt_thd) { ei->flags |= MAT::Edge::Flags::Skeleton; } }
				}
			}
		}

		// Delete the possible non-stable skeleton
		std::queue<MAT::Vertex_handle> queue;
		{
			MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
			CGAL_For_all(e, e_end) { queue.push(e->vertex()); }
		}
		while (!queue.empty())
		{
			MAT::Facet_handle fr;
			MAT::Edge_handle er;
			MAT::Vertex_handle vr, v = queue.front(); queue.pop();
			MAT::Vertex::Edge_iterator ev = v->edges_begin(), end = v->edges_end();
			for (; ev != end; ++ev)
			{
				MAT::Facet_handle fa = ev->facet();
				if (fa->weight() < wt_thd) { continue; }
				if (!(fa->flags & MAT::Facet::Flags::Eaten)) { break; }
		
				MAT::Edge_handle en = ev->next_around_facet();
				if (en->flags & MAT::Edge::Flags::StableSkeleton) { break; }
				else if (en->flags & MAT::Edge::Flags::Skeleton)
				{
					MAT::Vertex_handle vv = en->next_around_facet()->vertex();
					if (vr == MAT::Vertex_handle() || vr == vv) { vr = vv; fr = fa; er = en; } else { break; }
				}

				if (ev->flags & MAT::Edge::Flags::StableSkeleton) { break; }
				else if (ev->flags & MAT::Edge::Flags::Skeleton)
				{
					MAT::Vertex_handle vv = ev->vertex();
					if (vr == MAT::Vertex_handle() || vr == vv) { vr = vv; fr = fa; er = ev; } else { break; }
				}
			}
	  
			if (ev == end && vr != MAT::Vertex_handle())
			{
				er->flags &= ~MAT::Edge::Flags::Skeleton;
				MAT::Edge::Edge_identified_circulator ei = er->circ_identified();
				for (; ei != er; ++ei) { if (ei->facet()->weight() >= wt_thd) { ei->flags &= ~MAT::Edge::Flags::Skeleton; } }
				queue.push(vr);
			}
		}
		//----------------------------------------------------------------------------

		// Propagate
		{
			MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
			CGAL_For_all(e, e_end)
			{
				MAT::Facet_handle fr;
				MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
				for (; ei != e; ++ei)
				{
					MAT::Facet_handle fa = ei->facet();
					if (fa->weight() < wt_thd || fa->flags & MAT::Facet::Flags::Eaten) { continue; }
					if (fr != MAT::Facet_handle()) { fr = MAT::Facet_handle(); break; }
					fr = fa;
				}
				if (fr != MAT::Facet_handle() && !(fr->flags & MAT::Facet::Flags::InHeap))
				{
					facet_heap.push(fr, fr->weight());
					fr->flags |= MAT::Facet::Flags::InHeap;
				}
			}
		}
	}

	// Unset flag
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f) { f->flags &= ~MAT::Facet::Flags::InHeap; }
}

///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton data structure using the Skeleton flags in the MAT data.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton_internal(MAT* mat, double wt_thd)
{
	typedef handle_map<MAT::Vertex_handle, Skeleton3::vertex_descriptor> mat2skel;
	mat2skel verts_map;

	// Set visited flag on non-skeleton edges
	size_t n = 0;
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end) { if (!(e->flags & MAT::Edge::Flags::Skeleton)) { e->flags |= MAT::Edge::Flags::Visited; } else { ++n; } }
	}

	// Go through each edge and add it to the output (if it is a skeleton edge)	
	Skeleton3* S = new Skeleton3();
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		if (f->weight() < wt_thd) { continue; }
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end)
		{
			if (e->flags & MAT::Edge::Flags::Visited) { continue; }

			// Get and create endpoints of edge in the skeleton
			const MAT::Vertex_handle v0 = e->vertex(), v1 = e->next_around_facet()->vertex();
			Skeleton3::vertex_descriptor vs, ve;
			mat2skel::iterator i;
			i = verts_map.find(v0);
			if (i != verts_map.end()) { vs = i->second; }
			else
			{
				Skeleton3::vertex_property_type v = { v0->point() };
				verts_map.insert(std::make_pair(v0, vs = boost::add_vertex(v, *S)));
			}
			i = verts_map.find(v1);
			if (i != verts_map.end()) { ve = i->second; }
			else
			{
				Skeleton3::vertex_property_type v = { v1->point() };
				verts_map.insert(std::make_pair(v1, ve = boost::add_vertex(v, *S)));
			}

			// Mark the edge as visited and add it to the skeleton
			e->flags |= MAT::Edge::Flags::Visited;
			MAT::Edge::Edge_identified_circulator ei = e->circ_identified();
			for (; ei != e; ++ei) { ei->flags |= MAT::Edge::Flags::Visited; }
			boost::add_edge(vs, ve, *S);
		}
	}

	// Clear visited flag
	for (MAT::Facet_iterator f = mat->facets_begin(), end = mat->facets_end(); f != end; ++f)
	{
		MAT::Facet::Edge_around_facet_circulator e = f->edges_circ(), e_end = e;
		CGAL_For_all(e, e_end) { e->flags &= ~MAT::Edge::Flags::Visited; }
	}

	return S;
}

///////////////////////////////////////////////////////////////////////////////////////
// Constructs the skeleton data structure using the Skeleton flags in the MAT data.
///////////////////////////////////////////////////////////////////////////////////////
Skeleton3* construct_skeleton(MAT* mat, double flux_thd, double wt_thd)
{
	assert(flux_thd >= 0.0 && flux_thd < 1.0);
	assert(wt_thd >= 0.0);
	compute_skeleton(mat, flux_thd, wt_thd);
	return construct_skeleton_internal(mat, wt_thd);
}
