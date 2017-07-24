#pragma once

#include "GeometryTypes.hpp"
#include <string>
#include <iostream>
#include <unordered_map>

///////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for reading and writing various geometric formats.
//
// Supported polyhedral mesh formats:
//   OBJ - common in the 3D modeling world, each file can contain multiple, possibly named, meshes
//   OFF - common in the computational geometry field, they specify a single 3D mesh each
// In both formats normals, colors, textures, etc are all ignored. Multi-object OFF files are not
// supported. Reading OBJ files support normals only if compiled with polyhedron support for custom
// normals.
//
// OFF File Reference: http://www.geomview.org/docs/html/OFF.html
//
// Reading functions:
//   These return a new object allocated with "new" and must be deleted with "delete" when no loner
//   needed. They take an optional "assume_good" boolean, defaulting to false. When false, meshes
//   are fully checked to make sure they are suitable for being processed (which can take a long
//   time for large meshes). In general, as long as you know the mesh is good, you will want to
//   give true to greatly speed this up. You can always call check_mesh on your own.
//
// Writing functions:
//   ... TODO
//   TODO: OBJ file writing is currently unsupported
//
// Format-specific options are given as a string which is a comma-separated list of name=value
// items (or name items for "boolean" options). The names are case-insensitive. If creating
// file_options objects yourself, make sure to use lowercase names.
//
// Read Options:
//   OBJ - ObjectName  - load the given object by name
//         GroupName   - load the given group by name
//         ObjectIndex - load the given object by index
//         GroupIndex  - load the given group by index
//         Default is to load the first group, GroupIndex=0
//   OFF - (no options)
// 
// Writing Options:
//   OBJ - ObjectName, GroupName, ObjectIndex, GroupIndex as per reading options
//         Append      - if given and file exists, data is appended to existing file
//         Default is to create a new file with a single, unnamed, group with the mesh
//         If Append is given, one of ObjectName or GroupName must be given as well
//   OFF - Binary      - if given the file in written in the binary OFF format
//
///////////////////////////////////////////////////////////////////////////////////////////////////

extern const std::string mesh_file_usage;
extern const std::string mesh_file_read_usage;
extern const std::string mesh_file_write_usage;

void check_mesh(Polyhedron3* P);

typedef std::unordered_map<std::string, std::string> file_options;
enum class file_type { UNKNOWN, OBJ, OFF };

// Parses a filename and options string into a filename and options. A colon (:) separates the
// parts. The : cannot be the first character (or second on Windows). This means that on Windows
// single-character filenames followed by options must be specified like ./x:options. The options
// are parsed as per parse_file_options.
file_options parse_filename_and_options(const std::string& filename_and_options, std::string& filename);
// Parses options which are a comma-separated (,) list of name=value pairs. The name and values
// cannot contain commas and the name cannot contain =.
file_options parse_file_options(const std::string& options);
// Gets the file type of a file simply based on file extension.
file_type get_file_type(const std::string& filename);

// Read a filename with options in the string as a Polyhedron3
Polyhedron3* read_mesh(const std::string& filename_and_options, bool assume_good = false);
// Read a filename with a separate options string as a Polyhedron3
Polyhedron3* read_mesh(const std::string& filename, const std::string& options, bool assume_good = false);
// Read a filename with a separate options object as a Polyhedron3
Polyhedron3* read_mesh(const std::string& filename, const file_options& options, bool assume_good = false);
// Read a stream object as a Polyhedron3
Polyhedron3* read_mesh(std::istream& in, file_type type, const file_options& options = file_options(), bool assume_good = false);

void write_mesh(const Polyhedron3* P, const std::string& filename_and_options);
void write_mesh(const Polyhedron3* P, const std::string& filename, const std::string& options);
void write_mesh(const Polyhedron3* P, const std::string& filename, const file_options& options);
void write_mesh(const Polyhedron3* P, std::ostream& out, file_type type, const file_options& options = file_options());


///////////////////////////////////////////////////////////////////////////////////////////////////
// CG Files - don't know much about this, but they have a single skeleton each and are similar to
// OBJ files except use e for edge instead of f for face
///////////////////////////////////////////////////////////////////////////////////////////////////
Skeleton3* read_cg(const std::string& filename);
Skeleton3* read_cg(std::istream& in);
void write_cg(const Skeleton3* S, const std::string& filename);
void write_cg(const Skeleton3* S, std::ostream& out);
