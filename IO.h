#pragma once

#include "GeometryTypes.h"
#include <iostream>

///////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for reading and writing various geometric formats.
// Many of the reading functions take an optional "assume_good" boolean, defaulting to false. When
// false, large meshes will take a very long time to read because they are fully checked to make
// sure they are suitable for being processed. In general, as long as you know the mesh is good,
// you will want to give true to greatly speed this up.
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// OFF Files - common in the computational geometry field, they specify a single 3D mesh each.
// Note that only the basic GeomView format is supported, but both ASCII and binary versions can
// be read (only ASCII version is written). Normals, colors, textures, etc are all ignored. Lists
// are not supported.
///////////////////////////////////////////////////////////////////////////////////////////////////
Polyhedron3* read_off_file(const char* filename, bool assume_good = false);
Polyhedron3* read_off(std::istream &in, bool assume_good = false);
void write_off_file(const char* filename, Polyhedron3* P, bool verbose=true);
void write_off(std::ostream &out, Polyhedron3* P, bool verbose=true);

///////////////////////////////////////////////////////////////////////////////////////////////////
// CG Files - don't know much about this, but they have a single skeleton each.
///////////////////////////////////////////////////////////////////////////////////////////////////
Skeleton3* read_cg(const char* filename);

///////////////////////////////////////////////////////////////////////////////////////////////////
// OBJ Files - common in the 3D modeling world, each file can contain multiple meshes which you
// can look up by index or name. Currently writing multiple meshes to a single OBJ file is not
// directly supported (but could be done with the write_obj() functions).
// Note that normals, colors, textures, etc are all ignored (except if caching normals, then they
// are used).
///////////////////////////////////////////////////////////////////////////////////////////////////
class ObjFile
{
	bool _as_objects, _assume_good;

	struct internal_data;
	internal_data* data;

public:
	ObjFile(const char *filename, bool as_objects = false, bool assume_good = false);
	~ObjFile();

	inline const bool assume_good() const { return this->_assume_good; }
	inline void set_assume_good(bool x) { this->_assume_good = x; }

	inline const bool as_objects() const { return this->_as_objects; }
	inline void set_as_objects(bool x) { this->_as_objects = x; }

	const std::vector<std::string> names() const;
	size_t count() const;

	Polyhedron3* operator[](const std::string&) const; // need to "delete" return value
	Polyhedron3* operator[](size_t) const;
};
//void write_obj_file(const char* filename, const PolyhedronCollection& P, bool as_objects = false);
void write_obj_file(const char* filename, Polyhedron3* P);
void write_obj(std::ostream &out, Polyhedron3* P);
void write_obj(std::ostream &out, Polyhedron3* P, size_t& off);
