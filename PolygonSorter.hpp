#pragma once

#include "GeometryTypes.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Sorts a collection of polygons by their distance from a plane (defined by a point and orthogonal
// direction). This is designed to be used by the GL rendering utilities which need to render
// things from bottom-to-top.
///////////////////////////////////////////////////////////////////////////////////////////////////

class PolygonSorter
{
    struct internal_data;
    internal_data* data;

    PolygonSorter(const PolygonSorter&);
    PolygonSorter& operator= (PolygonSorter const&);

public:
    typedef const unsigned int* ConstFace;
    typedef const ConstFace* ConstFaces;

    PolygonSorter(const double* vertices, ConstFaces faces, size_t faces_len);
    // vertices is a x0, y0, z0, x1, y1, z1, ... array
    // faces is an array of 3-element indices into vertices and is of length faces_len
    // vertices must have 3 times the max values in faces
    
    ~PolygonSorter();

    void query(const Point3& viewing_source, const Direction3& viewing_dir, unsigned int* buf);
    // the faces vertex indices are written out to buffer, which must be 3 times the faces_len given to the constructor
};
