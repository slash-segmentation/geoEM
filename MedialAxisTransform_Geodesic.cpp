#include "MedialAxisTransform_Types_Geodesic.hpp"

#include "Heap.hpp"
#include "Progress.hpp"

#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// Compute the geodesic distance and derivative on mesh between two endpoints of all
// Delaunay edges. The results are stored in the mesh itself and uses data about the
// triangulation stored in the mesh.
// This uses the "flatexact_fast" method (instead of flatexact or Dijkstra methods).
///////////////////////////////////////////////////////////////////////////////////////
void compute_geodesic(Polyhedron3* mesh, VData& vdata)
{
    Progress progress("Computing geodesic distances...", mesh->size_of_vertices());
    progress.start();

    // Basic collections that are reused frequently
    typedef heap< double, Polyhedron3::Halfedge_const_handle, std::less<double>, handle_map> EdgeHeap;
    EdgeHeap edge_heap; edge_heap.reserve(1024);
    GeodesicDistances geod((size_t)(1.1*mesh->size_of_vertices())); // making the max load factor ~0.91
    typedef handle_map<Polyhedron3::Vertex_const_handle, handle_set<Polyhedron3::Facet_const_handle>> LECF;
    LECF le_cf(128); // the checking flag of the linking edges
    handle_set<Polyhedron3::Vertex_const_handle> dv_df(128); // the done flag of the vertices having the Delaunay edge connected with
    HEData edata(mesh->size_of_halfedges());
    for (Polyhedron3::Halfedge_iterator e = mesh->halfedges_begin(), end = mesh->halfedges_end(); e != end; ++e)
    {
        edata[e->id()] = HalfedgeData(e);
    }

    for (Polyhedron3::Vertex_iterator v = mesh->vertices_begin(), vend = mesh->vertices_end(); v != vend; ++v)
    {
        //std::cout << std::endl << std::endl << "source point: " << vertex_number(v, mesh) << std::endl;
        GeodesicDistances& v_geod = vdata[v->id()].geod;

        // Init flags
        size_t n_inc_tri_verts = v_geod.size();
        le_cf.clear(); le_cf.rehash((size_t)(1.1*n_inc_tri_verts));
        dv_df.clear(); dv_df.rehash((size_t)(1.1*n_inc_tri_verts));
        for (GeodesicDistances::const_iterator i = v_geod.begin(), end = v_geod.end(); i != end; ++i)
        {
            le_cf.insert(std::make_pair(i->first, handle_set<Polyhedron3::Facet_const_handle>(
                (size_t)(1.1*vdata[i->first->id()].n_inc_facets))));
        }
        
        // Compute distance/direction to all adjacent vertices and get the starting edges (they are the opposite edges of the vertex we are circulating, not including border edges)
        geod.clear();
        geod[v] = GeodesicDistance(0, Vector3());
        edge_heap.clear();
        //std::cout << "vid: " << vertex_number(v, mesh) << " #adj: " << v->degree() << std::endl;
        FOR_EDGES_AROUND_VERTEX(v, e)
        {
            Polyhedron3::Halfedge_const_handle ep = e->prev(); // the edge opposite v
            Polyhedron3::Vertex_const_handle va = ep->vertex();
            double d = edata[e->id()].length;
            geod[va] = GeodesicDistance(d, (va->point()-v->point()) / d);
            //std::cout << "the new geod for " << vertex_number(va, mesh) << " : " << d <<std::endl;
            if (v_geod.count(va)) { dv_df.insert(va); }

            //std::cout << "eid: " << edge_number(e, mesh) << " fid_inc:  " << facet_number(e->facet(), mesh) << std::endl;
            if (!ep->is_border_edge()) { edge_heap.push_raw(edata[ep->id()].eid, d); }
        }
        edge_heap.heapify();

        // Terminate the while loop if exhaust the edge_heap or we compute the geodesic distance to all adjacent vertices
        while (!edge_heap.empty() && dv_df.size() != n_inc_tri_verts)
        {
            // Get the closest edge along with its vertices
            Polyhedron3::Halfedge_const_handle e = edge_heap.top(); edge_heap.pop();
            const Polyhedron3::Vertex_const_handle v1 = e->prev()->vertex(), v2 = e->vertex();

            // If either vertex is the v itself, there is no need to compute
            if (v1 == v || v2 == v) { continue; }

            //const double L = distance(v2->point(), v1->point());
            const double L = edata[e->id()].length;
            const Kernel::FT L_2 = edata[e->id()].squared_length;
            //std::cout << std::endl << std::endl;
            //std::cout << edge_number(e, mesh) << " length " << L << std::endl;
            double C1 = get_dist(geod, v1);
            double C2 = get_dist(geod, v2);
            //CGAL_assertion(C1 + C2 >= L);
            CGAL_assertion(C1 >= 0 && C2 >= 0 && L >= 0);
            if      (L + C1 < C2) { C2 = L + C1; }
            else if (L + C2 < C1) { C1 = L + C2; }
            else if (C1 + C2 < L)
            {
                std::cerr << " v: " << v->point() << " e: " << halfedge_to_segment3(e) << std::endl;
                std::cerr << " L: " << L << " C1: " << C1 << " C2: " << C2 << " C1+C2: " << C1 + C2 << std::endl;
                C1 = L - C2;
            }
            const Kernel::FT C1_2 = C1 * C1, C2_2 = C2 * C2;
            const Kernel::FT ncos_a1 = (C1_2 + L_2 - C2_2), dcos_a1_2 = 4 * L_2 * C1_2;
            const Kernel::FT ncos_a2 = (C2_2 + L_2 - C1_2), dcos_a2_2 = 4 * L_2 * C2_2;
            
            Polyhedron3::Halfedge_const_handle e1 = e->prev(), e2 = e->next();
            for (int s = 0; s < 2; ++s)
            {
                const Polyhedron3::Facet_const_handle f = e->facet();
                const Polyhedron3::Vertex_const_handle v0 = e->next()->vertex();

                //std::cout << std::endl;
                //std::cout << "fid: " << facet_number(f, mesh) << " eid: " << edge_number(e, mesh) << " index: " << "?" << std::endl;
                
                //std::cout << "vid0: " << vertex_number(v0, mesh) << " " << v0->point() << std::endl;
                //std::cout << "vid1: " << vertex_number(v1, mesh) << " " << v1->point() << std::endl;
                //std::cout << "vid2: " << vertex_number(v2, mesh) << " " << v2->point() << std::endl;

                //const Kernel::FT d1_2 = CGAL::squared_distance(v0->point(), v1->point());
                //const Kernel::FT d2_2 = CGAL::squared_distance(v0->point(), v2->point());
                const Kernel::FT d1_2 = edata[e1->id()].squared_length;
                const Kernel::FT d2_2 = edata[e2->id()].squared_length;
                const Kernel::FT ncos_b1 = (d1_2 + L_2 - d2_2), dcos_b1_2 = 4 * L_2 * d1_2;
                const Kernel::FT ncos_b2 = (d2_2 + L_2 - d1_2), dcos_b2_2 = 4 * L_2 * d2_2;
                
                //std::cout << edge_number(e->prev(), mesh) << " length" << ft_sqrt(d1_2) << std::endl;
                //std::cout << edge_number(e->next(), mesh) << " length" << ft_sqrt(d2_2) << std::endl;


                double dist, v0_dist = get_dist(geod, v0);
                Vector3 ddist;
                if (ft_sqrt(CGAL::abs((dcos_a1_2 - ncos_a1 * ncos_a1) / (dcos_b1_2 - ncos_b1 * ncos_b1))) * ncos_b1 < -ncos_a1)
                {
                    const double d1 = dbl_sqrt(d1_2);
                    dist = C1 + d1;
                    if (dist < v0_dist) { ddist = (v0->point() - v1->point()) / d1; }
                    //std::cout << "vid1: " << "dist: " << dist << std::endl;
                }
                else if (ft_sqrt(CGAL::abs((dcos_a2_2 - ncos_a2 * ncos_a2) / (dcos_b2_2 - ncos_b2 * ncos_b2))) * ncos_b2 < -ncos_a2)
                {
                    const double d2 = dbl_sqrt(d2_2);
                    dist = C2 + d2;
                    if (dist < v0_dist) { ddist = (v0->point() - v2->point()) / d2; }
                    //std::cout << "vid2: " << "dist: " << dist << std::endl;
                }
                else if (L < 1e-5)
                {
                    const double d1 = dbl_sqrt(d1_2);
                    dist = C1 + d1;
                    ddist = (v0->point() - v1->point()) / d1;
                    std::cerr << "L: " << L << std::endl;
                }
                else
                {
                    dist = dbl_sqrt(CGAL::abs(C1_2 + d1_2 - (ncos_a1 * ncos_b1 - ft_sqrt(CGAL::abs((dcos_a1_2 - ncos_a1 * ncos_a1) * (dcos_b1_2 - ncos_b1 * ncos_b1)))) / (2 * L_2)));
                    if (dist < v0_dist)
                    {
                        Vector3 nn = normalized(CGAL::cross_product(v0->point() - v1->point(), v2->point() - v1->point()));
                        Vector3 xx = (v2->point() - v1->point()) / L;
                        Vector3 yy = CGAL::cross_product(nn, xx);
                        Vector3 cc1 = (xx * ncos_a1 + yy * ft_sqrt(CGAL::abs(dcos_a1_2 - ncos_a1 * ncos_a1))) * (1 / (2 * L));
                        ddist = normalized(v0->point() - (v1->point() + cc1));
                    }
                    //std::cout << "vid12: " << "dist: " << dist << std::endl;
                }

                // If the new geod_dist is less than the current one
                if (dist + 1e-10 < v0_dist)
                {
                    // Assign the new geod dist to vertex
                    geod[v0] = GeodesicDistance(dist, ddist);
                    
                    //std::cout << "the new geod for " << vertex_number(v0, mesh) << " : " << dist << " from edge: " << vertex_number(v1, mesh) << "--" << vertex_number(v2, mesh) << std::endl;

                    // Add the incident edges of v0 to heap if necessary
                    FOR_EDGES_AROUND_VERTEX(v0, e_circ)
                    {
                        //std::cout << "try to propagate to the adjacent facet: " << facet_number(e_circ->facet(), mesh) << std::endl;
                        double dist_circ = min2(dist, get_dist(geod, e_circ->prev()->vertex())); // endpoint of e_circ that is not v0
                        Polyhedron3::Halfedge_const_handle _e_circ = edata[e_circ->id()].eid;
                        EdgeHeap::handle h = edge_heap.find(_e_circ);
                        if (h != nullptr)
                        {
                            if (edge_heap.decrease(h, dist_circ))
                            {
                                //std::cout << "modify the edge: " << edge_number(e_circ, mesh) << " dist: " << dist_circ << std::endl;
                            }
                        }
                        else
                        {
                            edge_heap.push(_e_circ, dist_circ);
                            //std::cout << "insert the edge: " << edge_number(e_circ, mesh) << " dist: " << dist_circ << std::endl;
                        }
                    }
                }

                LECF::iterator i = le_cf.find(v0);
                if (i != le_cf.end() && i->second.insert(f).second && vdata[i->first->id()].n_inc_facets == i->second.size()) { dv_df.insert(v0); }
                
                // Switch to the incident facet
                e = e->opposite();
                if (e->is_border()) { break; }
                e1 = e->next(); e2 = e->prev();
            }
        }

        // Assign geod
        for (GeodesicDistances::iterator i = v_geod.begin(), end = v_geod.end(); i != end; ++i)
        {
            i->second = geod[i->first];
            //double d = (i->second = geod[i->first]).distance;
            //if (max_geod < d) { max_geod = d; }
            //if (min_geod > d) { min_geod = d; }
        }

        //std::cerr << std::endl << vertex_number(v, mesh) << ": #v " << ", d "; 
        //for(GeodesicDistances::iterator i = v_geod.begin(), end = v_geod.end(); i != end; ++i) { std::cerr << vertex_number(i->first, mesh) << " # " << i->second.distance << ",  "; }
        //std::cerr<<std::endl;

        progress.update();
    }
    progress.done();
}
// Same as compute_geodesic but attempts to cache the data to a file (or read it from that file)
void compute_geodesic_using_cache(Polyhedron3* mesh, VData& vdata, const std::string filename, const std::string objname)
{
    std::string gdc_name = filename + (objname == "" ? "" : ("_" + objname)) + "_geodesic_cache.dat";
    struct stat org_stat, gdc_stat;
    int org_exists = stat(filename.c_str(), &org_stat);
    int gdc_exists = stat(gdc_name.c_str(), &gdc_stat);
    if (org_exists == 0 && gdc_exists == 0 && org_stat.st_mtime < gdc_stat.st_mtime)
    {
        std::fstream f(gdc_name, std::ios_base::in | std::ios_base::binary);
        if (f.bad()) { f.close(); goto CALCULATE; }
        else
        {
            CGAL::set_binary_mode(f);
            for (Polyhedron3::Vertex_handle i = mesh->vertices_begin(), end = mesh->vertices_end(); i != end; ++i)
            {
                GeodesicDistances& gds = vdata[i->id()].geod;
                std::vector<Polyhedron3::Vertex_const_handle> keys;
                keys.reserve(gds.size());
                for (GeodesicDistances::const_iterator j = gds.begin(), end = gds.end(); j != end; ++j) { keys.push_back(j->first); }
                std::sort(keys.begin(), keys.end());
                for (std::vector<Polyhedron3::Vertex_const_handle>::iterator j = keys.begin(), end = keys.end(); j != end; ++j)
                {
                    GeodesicDistance& gd = gds[*j];
                    f.read(reinterpret_cast<char*>(&gd.distance), sizeof(gd.distance));
                    f >> gd.derivative;
                    //if (gd.distance > max_geod) { max_geod = gd.distance; }
                    //if (gd.distance < min_geod) { min_geod = gd.distance; }
                }
            }
        }
        f.close();
    }
    else
    {
CALCULATE:
        compute_geodesic(mesh, vdata);

        std::fstream f(gdc_name, std::ios_base::out | std::ios_base::binary);
        if (!f.bad())
        {
            CGAL::set_binary_mode(f);
            for (Polyhedron3::Vertex_const_handle i = mesh->vertices_begin(), end = mesh->vertices_end(); i != end; ++i)
            {
                const GeodesicDistances& gds = vdata[i->id()].geod;
                std::vector<Polyhedron3::Vertex_const_handle> keys;
                keys.reserve(gds.size());
                for (GeodesicDistances::const_iterator j = gds.begin(), end = gds.end(); j != end; ++j) { keys.push_back(j->first); }
                std::sort(keys.begin(), keys.end());
                for (std::vector<Polyhedron3::Vertex_const_handle>::const_iterator j = keys.begin(), end = keys.end(); j != end; ++j)
                {
                    const GeodesicDistance& gd = gds.at(*j);
                    f.write(reinterpret_cast<const char*>(&gd.distance), sizeof(gd.distance));
                    f << gd.derivative;
                }
            }
        }
        f.close();
    }
}
