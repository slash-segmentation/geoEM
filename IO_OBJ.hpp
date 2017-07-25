#pragma once

#include "GeometryTypes.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>

class ObjFileDetail;

///////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for reading and writing OBJ Files.
// Common in the 3D modeling world, each file can contain multiple meshes which you can look up by
// index or name. Note that normals, colors, textures, etc are all ignored (except if caching
// normals, then they are used).
///////////////////////////////////////////////////////////////////////////////////////////////////
class ObjFile
{
	bool _as_objects;
	ObjFileDetail* detail;

public:
	ObjFile(const char *filename, bool as_objects = false);
	ObjFile(std::istream& in, bool as_objects = false);
	~ObjFile();

	inline const bool as_objects() const { return this->_as_objects; }
	inline void set_as_objects(bool x) { this->_as_objects = x; }

	const std::vector<std::string> names() const;
	size_t size() const;

	Polyhedron3* operator[](const std::string&) const; // need to "delete" return value
	Polyhedron3* operator[](size_t) const;
};


void write_obj_file(const char* filename, const std::map<std::string, Polyhedron3*>& P, bool as_objects = false,
					const std::string& matfile = "", const std::vector<std::string>& mtls = std::vector<std::string>());
void write_obj_file(const char* filename, const std::vector<Polyhedron3*>& P, bool as_objects = false,
					const std::string& matfile = "", const std::vector<std::string>& mtls = std::vector<std::string>());
void write_obj_file(const char* filename, const Polyhedron3* P);
void write_obj(std::ostream &out, const Polyhedron3* P);
void write_obj(std::ostream &out, const Polyhedron3* P, size_t& off);
