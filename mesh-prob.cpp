///////////////////////////////////////////////////////////////////////////////
// Meshes a probability map.
///////////////////////////////////////////////////////////////////////////////

#ifndef BOOST_PARAMETER_MAX_ARITY
#define BOOST_PARAMETER_MAX_ARITY 12
#endif

#include "GeometryTypes.hpp"
#include "IO.hpp"
#include "Strings.hpp"
#include "Image3.hpp"

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Domain
#include <CGAL/Gray_image_mesh_domain_3.h>
typedef CGAL::Gray_image_mesh_domain_3<Image3, Kernel, double> Mesh_domain;

// Triangulation
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/make_mesh_3.h>
#ifdef MULTITHREADED
#include "tbb/task_scheduler_init.h"
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3T3;

// Criteria
#include <CGAL/Mesh_criteria_3.h>
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
using namespace CGAL::parameters;

// Conversion of a triangulation to a polyhedron mesh
#include "facets_in_complex_3_to_triangle_mesh.h" // normally part of CGAL, but only in newest versions


static void usage(const char* err = nullptr, int exit_code=0)
{
    if (err) { std::cerr << err << std::endl << std::endl; }

    std::cerr << "Usage: mesh-data image-stack [options...] out-model" << std::endl << std::endl;
    std::cerr << "image-stack image stack to read, either an MRC file or set of PNG files" << std::endl;
    std::cerr << "            surrounded by [ and ] or another format supported by CGAL::ImageIO" << std::endl;
    std::cerr << "out-mesh    file to save mesh to" << std::endl;
    std::cerr << "The available options are:" << std::endl;
    std::cerr << "  -s x,y,z  the size of the voxels of the image data, some formats like MRC have" << std::endl;
    std::cerr << "            this information embedded but other will likely need this set" << std::endl;
    std::cerr << "  -t thresh the threshold to use (from 0.0 to 1.0), default is 0.5" << std::endl;
    std::cerr << "  -o value  the value outside the image to assume, default is 0.0" << std::endl;
#ifdef MULTITHREADED
    std::cerr << "  -n nt     the number of threads to use, default is all available (" << tbb::task_scheduler_init::default_num_threads() << ")" << std::endl;
#endif
    std::cerr << "" << std::endl;
    std::cerr << mesh_file_usage << std::endl;
    std::cerr << mesh_file_write_usage << std::endl;
    exit(exit_code);
}

int main(int argc, char** argv)
{
    if (argc == 1) { usage(); }
    if (argc < 3) { usage("ERROR: invalid number of arguments", 1); }
    int argi = 1;

    registerMRCFormat();
    registerPNGFormat();

    
    // Read image stack
    std::cout << "Reading image stack..." << std::endl;
    Image3 im;
    if (streq(argv[argi], "["))
    {
        // Set of files
        ++argi; // skip [
        std::vector<Image3> images;
        for (; argi < argc && !streq(argv[argi], "]"); ++argi)
        {
            images.emplace_back();
            if (!images.back().read(argv[argi])) { usage("failed to read input image", 2); }
        }
        if (argi == argc) { usage("ERROR: no matching ] for stack of images", 2); }
        if (!stack_images(images, im)) { usage("ERROR: failed to stack input images", 2); }
    }
    else if (!im.read(argv[argi])) { usage("ERROR: failed to read input image", 2); } // single file
    if (im.image()->vdim != 1) { usage("ERROR: input image stack is not grayscale", 2); }
    if (++argi == argc) { usage("ERROR: no output file given", 5); }


    // Convert image stack to doubles
    std::cout << "Converting image stack..." << std::endl;
    Image3 im_dbl;
    try { convert_image<double>(im, im_dbl); }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl << std::endl;
        usage("ERROR: unable to convert image stack to doubles", 2);
    }
    
    
    // Handle arguments
    double threshold = 0.5;
    double value_outside = 0;
#ifdef MULTITHREADED
    int nthreads = tbb::task_scheduler_init::default_num_threads();
#endif
    while (argi < argc)
    {
        // Spacing
        if (streq(argv[argi], "-s"))
        {
            if (++argi == argc) { usage("ERROR: -s requires argument", 3); }
            double vx, vy, vz;
            size_t n;
            if (sscanf(argv[argi], "%lf,%lf,%lf%zn", &vx, &vy, &vz, &n) != 3 || argv[argi][n] ||
                vx == 0 || vy == 0 || vz == 0) { usage("ERROR: -s has invalid argument ", 3); }
            im_dbl.image()->vx = vx; im_dbl.image()->vy = vy; im_dbl.image()->vz = vz;
            ++argi;
        }
        // Threshold
        else if (streq(argv[argi], "-t"))
        {
            if (++argi == argc) { usage("ERROR: -t requires argument", 3); }
            threshold = atof(argv[argi]);
            if (threshold <= 0.0 || threshold > 1.0) { usage("ERROR: -t has invalid argument ", 3); }
            ++argi;
        }
        // Value Outside
        else if (streq(argv[argi], "-o"))
        {
            if (++argi == argc) { usage("ERROR: -o requires argument", 3); }
            value_outside = atof(argv[argi]);
            if (value_outside <= 0 || value_outside > 1.0) { usage("ERROR: -o has invalid argument ", 3); }
            ++argi;
        }
#ifdef MULTITHREADED
        // Number of threads
        else if (streq(argv[argi], "-n"))
        {
            if (++argi == argc) { usage("ERROR: -n requires argument", 3); }
            nthreads = atoi(argv[argi]);
            if (nthreads <= 0 || nthreads > tbb::task_scheduler_init::default_num_threads()) { usage("ERROR: -n has invalid argument", 3); }
            ++argi;
        }
#endif
        // No more options
        else if (argi != argc-1) { usage("ERROR: multiple output files given", 5); }
        else { break; }
    }
    if (argi == argc) { usage("ERROR: no output file given", 5); }

    
    // Perform meshing
#ifdef MULTITHREADED
    tbb::task_scheduler_init init(nthreads);
#endif
    std::cout << "Performing meshing..." << std::endl;
    Mesh_domain domain(im_dbl, threshold, value_outside);
    Mesh_criteria criteria(facet_angle=30, facet_size=6, facet_distance=2,
                           cell_radius_edge_ratio=3, cell_size=8);
    C3T3 c3t3 = CGAL::make_mesh_3<C3T3>(domain, criteria);

    
    // Convert to surface mesh
    std::cout << "Converting volume mesh to surface mesh..." << std::endl;
    Polyhedron3 *P = new Polyhedron3();
    facets_in_complex_3_to_triangle_mesh(c3t3, *P);
    
    
    // Write polyhedron
    std::cout << "Saving mesh..." << std::endl;
    char* output = argv[argc-1];
    try { write_mesh(P, output); }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl << std::endl;
        usage("ERROR: unable to write mesh to output file", 5);
    }

    return 0;
}
