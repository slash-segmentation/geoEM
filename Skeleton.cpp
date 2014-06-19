#include "Skeleton.hpp"

#include "GeometryUtils.hpp"

SkeletonGraph3* construct_skeleton_graph(const Skeleton3* S)
{
	typedef handle_map<Skeleton3::Vertex_const_handle, SkeletonGraph3::BranchPoint_handle> vert2bp;
	vert2bp branch_points;

	SkeletonGraph3* SG = new SkeletonGraph3();

	// Find and create all branch points (have degree >2)
	for (Skeleton3::Vertex_const_iterator v = S->vertices_begin(), end = S->vertices_end(); v != end; ++v)
	{
		if (v->degree() > 2) { branch_points.insert(std::make_pair(v, SG->add_branch_point(v->point()))); }
	}
	
	// Search out from every branch point to make branches
	for (vert2bp::iterator bp = branch_points.begin(), end = branch_points.end(); bp != end; ++bp)
	{
		for (Skeleton3::Vertex::Edge_const_iterator i = bp->first->edges_begin(), end = bp->first->edges_end(); i != end; ++i)
		{
			// Travel along edge until we hit a vertex that doesn't have degree 2
			SkeletonGraph3::Branch b;
			Skeleton3::Vertex_const_handle v = bp->first;
			Skeleton3::Edge_const_handle e = i;
			for(;;)
			{
				b.push_back(v->point());
				v = e->opposite(v);
				if (v->degree() != 2) { break; }
				Skeleton3::Vertex::Edge_const_iterator ex = v->edges_begin();
				e = (e == ex) ? std::next(ex) : ex;
			}
			b.push_back(v->point());
			SG->add_branch(bp->second, v->degree() == 1 ? SkeletonGraph3::BranchPoint_handle() : branch_points[v], b);
		}
	}

	SG->shrink_to_fit();
	return SG;
}

void skeleton_reduce(SkeletonGraph3* SG, double threshold)
{
	size_t nverts_before = 0, nedges_before = 0, nverts_after = 0, nedges_after = 0;
	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), end = SG->branches_end(); B != end; ++B) { size_t s = B->size(); nverts_before += s; nedges_before += s-1; }

	for (SkeletonGraph3::Branch_iterator b = SG->branches_begin(), end = SG->branches_end(); b != end; ++b)
	{
		remove_collinear_points(&*b, false);
		remove_nearly_collinear_points(&*b, threshold, false);
	}

	for (SkeletonGraph3::Branch_const_iterator B = SG->branches_begin(), end = SG->branches_end(); B != end; ++B) { size_t s = B->size(); nverts_after += s; nedges_after += s-1; }
	std::cerr << "Removed " << nverts_before - nverts_after << " vertices/edges out of " << nverts_before << std::endl;
}
