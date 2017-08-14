#include "TriangulatePolyhedron.hpp"

#ifdef POLYHEDRON_USE_VECTOR
// When the Polyhedron uses a vector we cannot change the polyhedron so we just fail.

void triangulate_polyhedron(Polyhedron3* P) { throw std::invalid_argument("Cannot triangulate polyhedrons with vector-based polyhedrons"); }

#else

// Based on the undocumented CGAL/triangulate_polyhedron.h by Laurent Rineau

#include "GeometryUtils.hpp"

#include <CGAL/Modifier_base.h>
#include <CGAL/HalfedgeDS_decorator.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>

#include <queue>
#include <vector>
#include <utility> // make_pair

class Triangulate_modifier : public CGAL::Modifier_base<Polyhedron3::HalfedgeDS> 
{
    typedef Polyhedron3::HalfedgeDS      HDS;

    typedef Polyhedron3::Halfedge        Halfedge;
    typedef Polyhedron3::Halfedge_handle Halfedge_handle;
    typedef Polyhedron3::Facet           Facet;
    typedef Polyhedron3::Facet_iterator  Facet_iterator;
    typedef Polyhedron3::Facet_handle    Facet_handle;

    struct Face_info { Polyhedron3::Halfedge_handle e[3]; bool is_external; };

    typedef CGAL::Triangulation_2_projection_traits_3<Kernel>            P_traits;
    typedef CGAL::Triangulation_data_structure_2<CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle, P_traits>, CGAL::Constrained_triangulation_face_base_2<P_traits, CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits>>> TDS;
    typedef CGAL::Constrained_triangulation_plus_2<CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, CGAL::No_intersection_tag>> CDT;

public:
    Triangulate_modifier() { }
    void operator()(HDS& hds)
    {
        // One needs to store facet handles into a vector, because the list of
        // facets of the polyhedron will be modified during the loop, and
        // that invalidates the range [facets_begin(), facets_end()[.
        std::vector<Facet_handle> facets;
        facets.reserve(hds.size_of_faces());
        for (Facet_iterator fit = hds.faces_begin(), end = hds.faces_end(); fit != end; ++fit) { facets.push_back(fit); }

        // Iterates on the vector of facet handles
        for (std::vector<Facet_handle>::iterator fit_it = facets.begin(), end = facets.end(); fit_it != end; ++fit_it)
        {
            Facet_handle fit = *fit_it;
            CDT cdt(P_traits(to_vector(normal(fit))));

            Facet::Halfedge_around_facet_circulator he = fit->facet_begin(), he_end(he);
            CDT::Vertex_handle previous, first;
            do
            {
                CDT::Vertex_handle vh = cdt.insert(he->vertex()->point());
                if (first == 0) { first = vh; }
                vh->info() = he;
                if (previous != 0 && previous != vh) { cdt.insert_constraint(previous, vh); }
                previous = vh;
            }
            while(++he != he_end);
            cdt.insert_constraint(previous, first);

            // sets mark is_external
            for (CDT::All_faces_iterator fit = cdt.all_faces_begin(), end = cdt.all_faces_end(); fit != end; ++fit) { fit->info().is_external = false; }
            std::queue<CDT::Face_handle> face_queue;
            face_queue.push(cdt.infinite_vertex()->face());
            while (!face_queue.empty())
            {
                CDT::Face_handle fh = face_queue.front(); face_queue.pop();
                if (fh->info().is_external) { continue; }
                fh->info().is_external = true;
                for (int i = 0; i <3; ++i) { if (!cdt.is_constrained(std::make_pair(fh, i))) { face_queue.push(fh->neighbor(i)); } }
            }

            // then modify the polyhedron
            CGAL::HalfedgeDS_decorator<HDS> decorator(hds);
            decorator.make_hole(fit->halfedge());
            for(CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(), end = cdt.finite_edges_end(); eit != end; ++eit)
            {
                CDT::Face_handle fh = eit->first, opposite_fh = fh->neighbor(eit->second);
                const int index = eit->second, opposite_index = opposite_fh->index(fh);
                const CDT::Vertex_handle va = fh->vertex(cdt. cw(index)), vb = fh->vertex(cdt.ccw(index));

                if (!(fh->info().is_external && opposite_fh->info().is_external) && !cdt.is_constrained(*eit))
                {
                    // strictly internal edge
                    Halfedge_handle h = hds.edges_push_back(Halfedge(), Halfedge());
                    fh->info().e[index] = h;
                    opposite_fh->info().e[opposite_index] = h->opposite();
                    decorator.set_vertex(h, va->info()->vertex());
                    decorator.set_vertex(h->opposite(), vb->info()->vertex());
                }
                if( cdt.is_constrained(*eit) )
                {
                    if (!fh->info().is_external) { fh->info().e[index] = va->info(); }
                    if (!opposite_fh->info().is_external) { opposite_fh->info().e[opposite_index] = vb->info(); }
                }
            }
            for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(), end = cdt.finite_faces_end(); fit != end; ++fit)
            {
                if (!fit->info().is_external)
                {
                    Halfedge_handle h0 = fit->info().e[0], h1 = fit->info().e[1], h2 = fit->info().e[2];
                    CGAL_assertion( h0 != Halfedge_handle() );
                    CGAL_assertion( h1 != Halfedge_handle() );
                    CGAL_assertion( h2 != Halfedge_handle() );
                    typedef Halfedge::Base HBase;
                    h0->HBase::set_next(h1); decorator.set_prev(h1, h0);
                    h1->HBase::set_next(h2); decorator.set_prev(h2, h1);
                    h2->HBase::set_next(h0); decorator.set_prev(h0, h2);
                    decorator.fill_hole(h0);
                }
            }
        }
    }
};

void triangulate_polyhedron(Polyhedron3* P)
{
    Triangulate_modifier modifier;
    P->delegate(modifier);
}

#endif
