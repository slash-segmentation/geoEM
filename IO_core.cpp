#include "IO.hpp"
#include "GeometryUtils.hpp"

#include "Strings.hpp"

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <exception>
#include <stdexcept>

// Vertex reader class used by both ObjReader and CgReader
class VertexReader
{
protected:
	typedef Point3     Point;
	typedef Direction3 Direction;
	
	static char *read_command(char *buf, int buf_len, FILE *f, char **params)
	{
		if (feof(f)) { return nullptr; }

		// Read an entire line (possibly)
		if (fgets(buf, buf_len, f))
		{
			// Test to see if we got an entire line
			char *comment = strchr(buf, '#');
			size_t len = strlen(buf);
			if ((len == (size_t)(buf_len - 1)) && !feof(f) && buf[len - 1] != '\n' && buf[len - 1] != '\r')
			{
				if (comment)
				{
					// Just read and toss characters until we get to a newline
					char xbuf[256];
					while (!feof(f) && fgets(xbuf, 256, f) && strlen(buf) == 255 && !feof(f) && buf[255] != '\n' && buf[255] != '\r') { }
					if (ferror(f)) { fclose(f); throw std::invalid_argument("Error: problem reading file\n"/*, ferror(f)*/); }
				}
				else
				{
					// We need to increase buffer size (although a line of length 1024+ without using a comment would be ridiculous)
					throw std::invalid_argument("Error: line too long (max line length excluding comments is 1024)\n");
					//buf = (char*)realloc(buf, buf_len <<= 1);
					//if (buf == nullptr) { fclose(f); throw std::invalid_argument("Error: problem reading file\n", ferror(f)); }
					//*_buf = buf;
					//*_buf_len = buf_len;
					//if (fgets(buf + len, buf_len - len, f))
					//{
					//	goto PARSE_LINE;
					//}
					//else if (ferror(f)) { fclose(f); throw std::invalid_argument("Error: problem reading file\n", ferror(f)); }
				}
			}

			// Strip the comment and trim the line
			if (comment) { *comment = '\0'; }
			char *line = trim(buf);

			// Get the command and parameters
			char *space = strchr(line, " \t\v");
			if (space) { *space = '\0'; *params = ltrim(space + 1); }
			else { *params = line + strlen(line); } // no params
			return line; // the command
		}
		else if (ferror(f)) { fclose(f); throw std::invalid_argument("Error: problem reading file\n"/*, ferror(f)*/); }

		return nullptr;
	}

	void v(const char *params)
	{
		double x, y, z, w = 1.0;
		int c1, c2, n = sscanf(params, "%lf %lf %lf%n %lf%n", &x, &y, &z, &c1, &w, &c2);
		if (n == EOF || n < 3 || (n == 3 && params[c1]) || (n == 4 && params[c2])) { throw std::invalid_argument("Error: invalid file format"); }
		this->vertices.push_back(Point(x, y, z)); // NOTE: w is dropped since it is only used in free-form (which is not supported at the moment)
	}
#ifdef POLYHEDRON_CACHED_NORMALS
	void vn(const char *params)
	{
		double x, y, z;
		int c, n = sscanf(params, "%lf %lf %lf%n", &x, &y, &z, &c);
		if (n != 3 || params[c]) { throw std::invalid_argument("Error: invalid file format"); }
		this->normals.push_back(Direction(x, y, z));
	}
#else
	inline void vn(const char*) { }
#endif
	
	FILE *file;

	std::vector<Point> vertices;
#ifdef POLYHEDRON_CACHED_NORMALS
	std::vector<Direction> normals;
#endif

	virtual void parse_cmd(const char* cmd, char* params, void* extra) = 0;
	void init(void* extra)
	{
		char *cmd, *params, buf[1024];

		// Get all vertices and validate most of the file structure
		while ((cmd = read_command(buf, ARRAYSIZE(buf), this->file, &params)) != nullptr)
		{
			if (cmd[0] == 0) { continue; }
			else if (streq_case_insensitive(cmd, "v")) { this->v(params); }
			else if (streq_case_insensitive(cmd, "vn")) { this->vn(params); } // no-op if not caching normals
			else { this->parse_cmd(cmd, params, extra); }
		}

		// Cleanup
		this->vertices.shrink_to_fit();
#ifdef POLYHEDRON_CACHED_NORMALS
		this->normals.shrink_to_fit();
#endif
	}

public:
	VertexReader(const char* filename) : file(nullptr)
	{
		FILE *f = fopen(filename, "rb");
		if (f == nullptr || ferror(f)) { throw std::invalid_argument("Error: failed to open file"); }
		this->file = f;
	}
	~VertexReader() { if (this->file) { fclose(this->file); this->file = nullptr; } }
};

#ifdef _MSC_VER
#pragma region Reading CG Files
#endif
class CgReader : public VertexReader
{
private:
	void e(char *s)
	{
		// First make sure all vertices have been added to the skeleton
		for (size_t i = this->vert_handles.size(), len = this->vertices.size(); i < len; ++i)
		{
			this->vert_handles.push_back(this->skel->add_vertex(this->vertices[i]));
		}

		ssize_t v;
		Skeleton3::Vertex_handle vs[2];
		int i = 0;
		while (i < 2 && str_read_int(s, &s, &v) && (!*s || isspace(*s)))
		{
			if (v == 0 || (size_t)CGAL::abs(v) > this->vert_handles.size()) { break; }
			vs[i++] = this->vert_handles[(v > 0) ? v - 1 : v + this->vert_handles.size()];
			s = ltrim(s);
			if (!*s)
			{
				if (i != 2) { break; }
				this->skel->add_edge(vs[0], vs[1]);
				return;
			}
		}
		throw std::invalid_argument("Error: invalid file format");
	}
	virtual void parse_cmd(const char* cmd, char* params, void*)
	{
		switch (tolower(cmd[0]))
		{
		case 'e': e(params); break; // v and vn are already done, just need to get the edges
		default: throw std::invalid_argument("Error: CG file contains unimplemented commands");
		}
	}

	Skeleton3* skel;
	std::vector<Skeleton3::Vertex_handle> vert_handles;

public:
	CgReader(const char* filename) : VertexReader(filename), skel(new Skeleton3()) { this->init(nullptr); }

	Skeleton3* skeleton() { return this->skel; }
	const Skeleton3* skeleton() const { return this->skel; }
};

Skeleton3* read_cg(const char* filename)
{
	CgReader cg(filename);
	return cg.skeleton();
}

#ifdef _MSC_VER
#pragma endregion
#endif
