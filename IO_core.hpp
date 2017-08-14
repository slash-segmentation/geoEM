#pragma once

#include "GeometryTypes.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Vertex reader class used by both ObjReader and CgReader
class VertexReader
{
protected:
	typedef Point3     Point;
	typedef Direction3 Direction;

	static char* read_command(std::istream &f, std::string& buf, char **params);
	void v(const char *params);

	std::vector<Point> vertices;

	virtual void parse_cmd(const char* cmd, char* params, void* extra) = 0;
	void init(std::istream &in, void* extra);
};
