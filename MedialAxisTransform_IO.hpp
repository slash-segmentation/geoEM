#pragma once

#include <string>
#include "GeometryTypes.hpp"
#include "MedialAxisTransform_Types_MAT.hpp"

// Load a MAT from a binary file
MAT* load_mat(const std::string filename, Polyhedron3* mesh);

// Save a MAT from a binary file
void save_mat(const std::string filename, MAT* mat, Polyhedron3* mesh);

// Dump a MAT to stdout for debugging
void dump_mat(MAT* mat);

