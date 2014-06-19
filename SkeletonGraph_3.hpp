#pragma once

#include "GeometryTypes.hpp"
#include "Iterators.hpp"
#include <list>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////////////////
// A skeleton is simply a graph with vertices and edges in which every vertex has a 3D point.
// This class represents the skeleton such that all sets of sequential vertices of degree 2
// (together with the neighboring non-degree-2 vertices) are a single branch. Vertices with degree
// greater than 2 are made into branch points. Vertices of degree 1 are implicitly represented.
//
// Note: to reduce confusion between branch points and points on branches, spatial points are
// called coordinates.
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
// A branch point in a SkeletonGraph_3. It will have at least 3 incident branches and a coordinate.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename SG3>
class SkeletonGraph_3_BranchPoint
{
	friend SG3;
public:
	typedef typename SG3::Coordinate  Coordinate; // not called Point here due to possible confusion of Branch::Point with BranchPoint
	typedef typename SG3::Branch      Branch;
	typedef typename SG3::BranchPoint BranchPoint;
	typedef typename SG3::Branch_handle      Branch_handle;
	typedef typename SG3::BranchPoint_handle BranchPoint_handle;
	typedef typename SG3::Branch_const_handle      Branch_const_handle;
	typedef typename SG3::BranchPoint_const_handle BranchPoint_const_handle;
	typedef std::list<Branch_handle> Branch_list;
	typedef handle_iterator<typename Branch_list::iterator,      Branch_handle>       Branch_iterator;
	typedef handle_iterator<typename Branch_list::const_iterator,Branch_const_handle> Branch_const_iterator;
private:
	Coordinate c;
	Branch_list bs;
public:
	inline SkeletonGraph_3_BranchPoint(const Coordinate& c) : c(c), bs() { }
	inline bool is_valid() const
	{
		for (Branch_const_iterator b = this->bs.begin(), end = this->bs.end(); b != end; ++b)
		{
			bool is_source = *b->source() == this, is_target = *b->target() == this;
			if (is_source == is_target || (is_source && b->front() != this->c) || (is_target && b->back() != this->c)) { return false; }
		}
		return true;
	}
	inline size_t degree() const { return this->bs.size(); } // number of branches
	inline const Coordinate& coord() const { return this->c; }
	inline       Branch_iterator branches_begin()       { return this->bs.begin(); }
	inline Branch_const_iterator branches_begin() const { return this->bs.begin(); }
	inline       Branch_iterator branches_end()         { return this->bs.end();   }
	inline Branch_const_iterator branches_end()   const { return this->bs.end();   }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// A branch in a SkeletonGraph_3. It will have at least 2 coordinates. It may have source and
// target branch points. It can be used anywhere a list of 3D coordinates/points can be used. The
// endpoint coordinate values should always be kept consistent with the source and target branch
// points.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename SG3>
class SkeletonGraph_3_Branch : public std::list<typename SG3::Coordinate> // a branch is just a list of Point3s, but you really shouldn't modify the first and last points once the source and target branch points are added
{
	// TODO: overload some of the functions from list to make sure they keep the data behaved - like reverse
	friend SG3;
public:
	typedef typename SG3::Coordinate  Coordinate; // not called Point here due to possible confusion of Branch::Point with BranchPoint
	typedef typename SG3::Branch      Branch;
	typedef typename SG3::BranchPoint BranchPoint;
	typedef typename SG3::Branch_handle      Branch_handle;
	typedef typename SG3::BranchPoint_handle BranchPoint_handle;
	typedef typename SG3::Branch_const_handle      Branch_const_handle;
	typedef typename SG3::BranchPoint_const_handle BranchPoint_const_handle;
private:
	typedef std::list<Coordinate> Base;
	BranchPoint_handle bps, bpt;
public:
	inline SkeletonGraph_3_Branch() : Base(), bps(), bpt() { }
	inline bool is_valid() const { return (this->bps == BranchPoint_handle() || this->bps.coord() == this->cs.front()) && (this->bpt == BranchPoint_handle() || this->bps.coord() == this->cs.back()); }
	inline bool has_source() const { return this->bps != BranchPoint_handle(); }
	inline bool has_target() const { return this->bpt != BranchPoint_handle(); }
	inline BranchPoint_handle       source()       { return this->bps; }
	inline BranchPoint_handle       target()       { return this->bpt; }
	inline BranchPoint_const_handle source() const { return this->bps; }
	inline BranchPoint_const_handle target() const { return this->bpt; }
	inline size_t index(BranchPoint_const_handle bp) const { if (bp == this->bps) { return 0; } else { assert(bp == this->bpt); return 1; } }
	inline BranchPoint_handle       opposite(BranchPoint_const_handle bp)       { if (bp == this->bps) { return this->bpt; } else { assert(bp == this->bpt); return this->bps; } }
	inline BranchPoint_const_handle opposite(BranchPoint_const_handle bp) const { if (bp == this->bps) { return this->bpt; } else { assert(bp == this->bpt); return this->bps; } }
	inline BranchPoint_handle       branch_point(size_t i)       { if (i == 0) { return this->bps; } else { assert(i == 1); return this->bpt; } }
	inline BranchPoint_const_handle branch_point(size_t i) const { if (i == 0) { return this->bps; } else { assert(i == 1); return this->bpt; } }
	//inline Base::iterator       begin_from(BranchPoint_const_handle bp)       { if (bp == this->bps) { return this->begin(); } else { assert(bp == this->bpt); return this->rbegin(); } }
	//inline Base::const_iterator begin_from(BranchPoint_const_handle bp) const { if (bp == this->bps) { return this->begin(); } else { assert(bp == this->bpt); return this->rbegin(); } }
	//inline Base::iterator       end_from  (BranchPoint_const_handle bp)       { if (bp == this->bps) { return this->end();   } else { assert(bp == this->bpt); return this->rend();   } }
	//inline Base::const_iterator end_from  (BranchPoint_const_handle bp) const { if (bp == this->bps) { return this->end();   } else { assert(bp == this->bpt); return this->rend();   } }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Represents a 3D 'skeleton' - essentially a graph with 3D points associated with each vertex
// but only the overall structure is explicitly represented.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, template<typename> class BP = SkeletonGraph_3_BranchPoint, template<typename> class B = SkeletonGraph_3_Branch>
class SkeletonGraph_3
{
public:
	typedef SkeletonGraph_3<K, BP, B> Self;
	typedef typename K::Point_3 Coordinate; // not called Point here due to possible confusion of Branch::Point with BranchPoint
	typedef K        Kernel;
	typedef BP<Self> BranchPoint;
	typedef B<Self>  Branch;
	
	typedef std::vector<BranchPoint> BranchPoint_vector;
	typedef std::vector<Branch>      Branch_vector;

	typedef stable_vector_iterator<      BranchPoint_vector, typename BranchPoint_vector::iterator>       BranchPoint_iterator;
	typedef stable_vector_iterator<const BranchPoint_vector, typename BranchPoint_vector::const_iterator> BranchPoint_const_iterator;
	typedef stable_vector_iterator<      Branch_vector,      typename Branch_vector::iterator>            Branch_iterator;
	typedef stable_vector_iterator<const Branch_vector,      typename Branch_vector::const_iterator>      Branch_const_iterator;

	typedef BranchPoint_iterator       BranchPoint_handle;
	typedef BranchPoint_const_iterator BranchPoint_const_handle;
	typedef Branch_iterator            Branch_handle;
	typedef Branch_const_iterator      Branch_const_handle;

private:
	BranchPoint_vector bps;
	Branch_vector      bs;

public:
	inline explicit SkeletonGraph_3(size_t bp = 0, size_t b = 0) : bps(), bs() { this->bps.reserve(bp); this->bs.reserve(b); }

	inline size_t size_of_branch_points()     const { return this->bps.size();     }
	inline size_t size_of_branches()          const { return this->bs.size();      }
	inline size_t capacity_of_branch_points() const { return this->bps.capacity(); }
	inline size_t capacity_of_branches()      const { return this->bs.capacity();  }
	inline bool empty() const { return this->bps.empty(); }
	inline void shrink_to_fit() { this->bps.shrink_to_fit(); this->bs.shrink_to_fit(); }
	
	inline BranchPoint_iterator branch_points_begin() { return BranchPoint_iterator(&this->bps, 0); }
	inline BranchPoint_iterator branch_points_end()   { return BranchPoint_iterator(&this->bps, this->bps.size()); }
	inline Branch_iterator      branches_begin()      { return Branch_iterator     (&this->bs, 0); }
	inline Branch_iterator      branches_end()        { return Branch_iterator     (&this->bs, this->bs.size()); }
	inline BranchPoint_const_iterator branch_points_begin() const { return BranchPoint_const_iterator(&this->bps, 0); }
	inline BranchPoint_const_iterator branch_points_end()   const { return BranchPoint_const_iterator(&this->bps, this->bps.size()); }
	inline Branch_const_iterator      branches_begin()      const { return Branch_const_iterator     (&this->bs, 0); }
	inline Branch_const_iterator      branches_end()        const { return Branch_const_iterator     (&this->bs, this->bs.size()); }

	// Add an isolated branch point at the given coordinate to the collection and return a handle to it.
	inline BranchPoint_handle add_branch_point(const Coordinate& c) { this->bps.push_back(BranchPoint(c)); return BranchPoint_handle(&this->bps, this->bps.size() - 1); }
	// Add an isolated branch point to the collection and return a handle to it.
	inline BranchPoint_handle add_branch_point(const BranchPoint& bp) { this->bps.push_back(bp); return BranchPoint_handle(&this->bps, this->bps.size() - 1); }
	// Add a branch to the collection. The branch is defined by the two branch points given. Both
	// branch points are are updated to know about the new incident branch and the branch is
	// updated to make sure it has the branch points at the end of its coordinate list. Return a
	// handle to the new branch. The branch points can be left out to create a "leaf" branch or
	// an isolated branch.
	inline Branch_handle add_branch(BranchPoint_handle bp0 = BranchPoint_handle(), BranchPoint_handle bp1 = BranchPoint_handle(), const Branch& b = Branch())
	{
		this->bs.push_back(b);
		Branch_handle bh = Branch_handle(&this->bs, this->bs.size() - 1);
		bool have_bp0 = bp0 != BranchPoint_handle(), have_bp1 = bp1 != BranchPoint_handle();
		if      (have_bp0 && bp0->coord() == bh->front()) { bh->bps = bp0; bh->bpt = bp1; if (have_bp1 && bp1->coord() != bh->back ()) { bh->push_back (bp1->coord()); } }
		else if (have_bp0 && bp0->coord() == bh->back ()) { bh->bps = bp1; bh->bpt = bp0; if (have_bp1 && bp1->coord() != bh->front()) { bh->push_front(bp1->coord()); } }
		else if (have_bp1 && bp1->coord() == bh->front()) { bh->bps = bp1; bh->bpt = bp0; if (have_bp0) { bh->push_back (bp0->coord()); } }
		else if (have_bp1 && bp1->coord() == bh->back ()) { bh->bps = bp0; bh->bpt = bp1; if (have_bp0) { bh->push_front(bp0->coord()); } }
		else                                              { bh->bps = bp0; bh->bpt = bp1; if (have_bp0) { bh->push_front(bp0->coord()); } if (have_bp1) { bh->push_back(bp1->coord()); } }
		if (have_bp0) { bp0->bs.push_back(bh); }
		if (have_bp1) { bp1->bs.push_back(bh); }
		return bh;
	}
};
