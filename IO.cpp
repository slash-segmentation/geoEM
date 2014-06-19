#include "IO.hpp"
#include "TriangulatePolyhedron.hpp"
#include "GeometryUtils.hpp"

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <exception>
#include <stdexcept>

#ifdef _MSC_VER
#pragma region String Utilities
#endif
inline char *trim(char *s) { while (isspace(*s)) { ++s; } size_t l = strlen(s); while (l && isspace(s[l - 1])) { s[--l] = 0; } return s; }
inline char *ltrim(char *s) { while (isspace(*s)) { ++s; } return s; }
inline const char *ltrim(const char *s) { while (isspace(*s)) { ++s; } return s; }
inline char *strchr(char *s, const char *cs) { while (*s) { const char *c = cs; while (*c) { if (*s == *c) { return s; } ++c; } ++s; } return nullptr; }
inline const char *strchr(const char *s, const char *cs) { while (*s) { const char *c = cs; while (*c) { if (*s == *c) { return s; } ++c; } ++s; } return nullptr; }
inline bool equal_case_insensitive(const char *a, const char *b)
{
	while (*a && *b && tolower(*a) == tolower(*b)) { ++a; ++b; }
	return *a == 0 && *b == 0;
}
inline bool str2int(const char *s, char **endptr, ssize_t *x)
{
	bool negative = *s == '-';
	if (negative || *s == '+') { ++s; }
	if (!isdigit(*s)) { return false; }
	ssize_t _x = 0;
	do
	{
		int64_t temp = ((int64_t)_x) * 10 + (*s++ - '0');
		if (temp > INT_MAX) { return false; } // overflowed
		_x = (int)temp;
	}
	while (isdigit(*s));
	*x = (negative ? -_x : _x);
	*endptr = (char*)s;
	return true;
}
inline bool strskipint(const char *s, char **endptr)
{
	if (*s == '-' || *s == '+') { ++s; }
	if (!isdigit(*s)) { return false; }
	while (isdigit(*++s));
	*endptr = (char*)s;
	return true;
}
#ifdef _MSC_VER
#pragma endregion
#endif

static void check(Polyhedron3 *P)
{
	// TODO: add check for isolated vertices
	// expensive operations that won't be necessary, usually
	if (!P->is_closed())       { throw std::invalid_argument("Error: non-closed polyhedron."); }
	if (!is_not_degenerate(P)) { throw std::invalid_argument("Error: degenerate polyhedron."); }
	if (!is_manifold(P))       { throw std::invalid_argument("Error: non-manifold polyhedron."); }
	if (!P->is_pure_triangle())
	{
		std::cerr << "Polyhedron was not pure triangle, will take some time to triangulate." << std::endl;
		triangulate_polyhedron(P);
	}
}

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
			else if (equal_case_insensitive(cmd, "v")) { this->v(params); }
			else if (equal_case_insensitive(cmd, "vn")) { this->vn(params); } // no-op if not caching normals
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

// OFF File Reference:
//   http://www.geomview.org/docs/html/OFF.html
#ifdef _MSC_VER
#pragma region OFF Files
#endif
Polyhedron3* read_off_file(const char* filename, bool assume_good)
{
	std::ifstream f(filename);
	Polyhedron3* P = read_off(f, assume_good);
	return P;
}
Polyhedron3* read_off(std::istream &in, bool assume_good)
{
	Polyhedron3* P = new Polyhedron3();
	in >> *P;
	if (!assume_good) { check(P); }
	return P;
}
void write_off_file(const char* filename, Polyhedron3* P, bool verbose)
{
	std::ofstream f(filename);
	write_off(f, P, verbose);
	f.close();
}
void write_off(std::ostream &out, Polyhedron3* P, bool verbose)
{
	CGAL::set_ascii_mode(out);
	if (verbose) { CGAL::set_pretty_mode(out); }
	out << *P;
}
#ifdef _MSC_VER
#pragma endregion
#endif

// OBJ File References:
//   http://en.wikipedia.org/wiki/Wavefront_.obj_file
//   http://www.martinreddy.net/gfx/3d/OBJ.spec

#ifdef _MSC_VER
#pragma region Writing OBJ Files
#endif
//void write_obj_file(const char* filename, const PolyhedronCollection& P, bool as_objects)
//{
//	size_t off = 0;
//	std::ofstream f(filename);
//	// TODO: handle no-name and "default" polyhedrons
//	for (PolyhedronCollection::const_iterator i = P.begin(), end = P.end(); i != end; ++i)
//	{
//		f << (as_objects ? 'o' : 'g') << ' ' << i->first << std::endl;
//		write_obj(f, i->second, off);
//	}
//	f.close();
//}
void write_obj_file(const char* filename, Polyhedron3* P)
{
	size_t off = 0;
	std::ofstream f(filename);
	write_obj(f, P, off);
	f.close();
}
void write_obj(std::ostream &out, Polyhedron3* P)
{
	size_t off = 0;
	write_obj(out, P, off);
}
void write_obj(std::ostream &out, Polyhedron3* P, size_t& off)
{
	CGAL::set_ascii_mode(out);
	for (Polyhedron3::Point_iterator i = P->points_begin(), end = P->points_end(); i != end; ++i) { out << "v " << *i << std::endl; }
	for (Polyhedron3::Facet_iterator i = P->facets_begin(), end = P->facets_end(); i != end; ++i)
	{
		out << 'f';
		FOR_VERTICES_AROUND_FACET(i, v) { out << ' ' << (vertex_number(v, P) + 1 + off); }
		out << std::endl;
	}
	off += P->size_of_vertices();
}
#ifdef _MSC_VER
#pragma endregion
#endif

#ifdef _MSC_VER
#pragma region Reading OBJ Files
#endif
class ObjReader : public VertexReader, public CGAL::Modifier_base<Polyhedron3::HalfedgeDS>
{
private:
	typedef Polyhedron3::HalfedgeDS HDS;
	typedef HDS::Vertex   Vertex;
	typedef CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
	
	typedef std::map<std::string, size_t> name2count;
	typedef name2count::iterator name2count_iter;
	
	inline static name2count_iter get_entry(const char *name, name2count &map)
	{
		name2count_iter i = map.find(name);
		return (i == map.end()) ? map.insert(i, name2count::value_type(name, 0)) : i;
	}
	std::vector<name2count_iter> g(char *params)
	{
		std::vector<name2count_iter> grps;
		char *name = strtok(params, " \t\v");
		while (name != nullptr)
		{
			grps.push_back(get_entry(name, this->grps));
			name = strtok(nullptr, " \t\v");
		}
		if (grps.size() == 0) { grps.push_back(get_entry("default", this->grps)); }
		return grps;
	}
	name2count_iter o(const char *params) { return get_entry(params, this->objs); }
	
#ifdef POLYHEDRON_CACHED_NORMALS
	void f_v(char *s, bool has_vt, Builder *b, size_t num_vertices)
	{
		b->begin_facet();
		ssize_t v;
		size_t count = 0;
		while (str2int(s, &s, &v) && // read "v"
				(!has_vt || *s == '/' && strskipint(s+1, &s)) && // skip "vt"
				(!*s || isspace(*s)))
		{
			// Add vertex to facet
			if (v == 0 || (size_t)CGAL::abs(v) > num_vertices) { break; }
			b->add_vertex_to_facet((v > 0) ? v - 1 : v + num_vertices);
			++count;
			s = ltrim(s);
			if (!*s) { if (count < 3) { break; } b->end_facet(); return; } // Finished
		}
		b->end_facet();
		throw std::invalid_argument("Error: invalid file format");
	}
	void f_v_vn(char *s, bool has_vt, Builder *b, size_t num_vertices, size_t num_normals, size_t& total_vertices)
	{
		b->begin_facet();
		Builder::Vertex_handle vv;
		ssize_t v, n;
		size_t count = 0;
		while (str2int(s, &s, &v) && // read "v"
				*s++ == '/' && (!has_vt || strskipint(s, &s)) && // skip "vt"
				*s == '/' && str2int(s+1, &s, &n) && (!*s || isspace(*s))) // read "vn"
		{
			// Add vertex with normal
			if (v == 0 || (size_t)CGAL::abs(v) > num_vertices || n == 0 || (size_t)CGAL::abs(n) > num_normals) { break; }
			v = (v > 0) ? v - 1 : v + num_vertices;
			n = (n > 0) ? n - 1 : n + num_normals;
			vv = b->vertex(v);
			if (vv->has_normal() && vv->normal() != this->normals[n])
			{
				vv = b->add_vertex(vv->point());
				v = total_vertices++;
			}
			vv->set_normal(this->normals[n]);
			b->add_vertex_to_facet(v);
			++count;
			s = ltrim(s);
			if (!*s) { if (count < 3) { break; } b->end_facet(); return; } // Finished
		}
		b->end_facet();
		throw std::invalid_argument("Error: invalid file format");
	}
	void f(char *params, Builder *b, size_t num_vertices, size_t num_normals, size_t& total_vertices)
	{
		// At least 3 points
		// Each point is of one of the following forms, and in a face all points have the same form
		//   v, v/vt, v//vn, v/vt/vn

		// Call the proper function for the given form by parsing the first vertex
		char *s = params;
		if (strskipint(s, &s))
		{
			if (isspace(*s)) { f_v(params, false, b, num_vertices); return; } // v
			else if (*s++ == '/')
			{
				if (*s == '/' && strskipint(s+1, &s)) { f_v_vn(params, false, b, num_vertices, num_normals, total_vertices); return; } // v//vn
				else if (strskipint(s, &s))
				{
					if (isspace(*s)) { f_v(params, true, b, num_vertices); return; } // v/vt
					else if (*s == '/' && strskipint(s+1, &s) && isspace(*s)) { f_v_vn(params, true, b, num_vertices, num_normals, total_vertices); return; } // v/vt/vn
				}
			}
		}
		throw std::invalid_argument("Error: invalid file format");
	}
#else
	void f_v(char *s, bool has_vt, bool has_vn, Builder *b, size_t num_vertices)
	{
		b->begin_facet();
		ssize_t v;
		size_t count = 0;
		while (str2int(s, &s, &v) && // read "v"
				((!has_vt && !has_vn) || *s++ == '/') &&
				(!has_vt || strskipint(s, &s)) && // skip "vt"
				(!has_vn || (*s == '/' && strskipint(s+1, &s))) && // skip "vn"
				(!*s || isspace(*s)))
		{
			// Add vertex to facet
			if (v == 0 || (size_t)CGAL::abs(v) > num_vertices) { break; }
			b->add_vertex_to_facet((v > 0) ? v - 1 : v + num_vertices);
			++count;
			s = ltrim(s);
			if (!*s) { if (count < 3) { break; } b->end_facet(); return; } // Finished
		}
		b->end_facet();
		throw std::invalid_argument("Error: invalid file format");
	}
	void f(char *params, Builder *b, size_t num_vertices)
	{
		// At least 3 points
		// Each point is of one of the following forms, and in a face all points have the same form
		//   v, v/vt, v//vn, v/vt/vn

		// Call the proper function for the given form by parsing the first vertex
		char *s = params;
		if (strskipint(s, &s))
		{
			if (isspace(*s)) { f_v(params, false, false, b, num_vertices); return; } // v
			else if (*s++ == '/')
			{
				if (*s == '/' && strskipint(s+1, &s)) { f_v(params, false, true, b, num_vertices); return; } // v//vn
				else if (strskipint(s, &s))
				{
					if (isspace(*s)) { f_v(params, true, false, b, num_vertices); return; } // v/vt
					else if (*s == '/' && strskipint(s+1, &s) && isspace(*s)) { f_v(params, true, true, b, num_vertices); return; } // v/vt/vn
				}
			}
		}
		throw std::invalid_argument("Error: invalid file format");
	}
#endif

	inline static bool o_matches(const char *params, const char* obj_name) { return strcmp(params, obj_name) == 0; }
	inline static bool g_matches(char *params, const char* grp_name) { char *name = strtok(params, " "); while (name != nullptr) { if (strcmp(name, grp_name) == 0) { return true; } name = strtok(nullptr, " "); } return false; }
	
	name2count grps, objs;

	std::string name;
	bool as_object;

	struct Current
	{
		name2count_iter obj;
		std::vector<name2count_iter> grps;
		Current(name2count_iter o, name2count_iter g) : obj(o), grps(1, g) {}
	};

	virtual void parse_cmd(const char* cmd, char* params, void* extra)
	{
		Current* cur = (Current*)extra;
		switch (tolower(cmd[0]))
		{
		// v and vn are already done, just need to count the number of faces in each group/object and get all vertices and validate ignored commands
		case 'o':
			if (cmd[1] == '\0') { cur->obj = o(params); }
			else { throw std::invalid_argument("Error: OBJ file contains unimplemented commands"); }
			break;
		case 'g':
			if (cmd[1] == '\0') { cur->grps = g(params); }
			else { throw std::invalid_argument("Error: OBJ file contains unimplemented commands"); }
			break;
		case 'f':
			if (cmd[1] == '\0' || equal_case_insensitive(cmd+1, "o")) // TODO: validation
			{
				++cur->obj->second;
				for (size_t i = 0, len = cur->grps.size(); i < len; ++i)
					++cur->grps[i]->second;
			}
			else { throw std::invalid_argument("Error: OBJ file contains unimplemented commands"); }
			break;
		default:
			if (!equal_case_insensitive(cmd, "vt") && !equal_case_insensitive(cmd, "s") && !equal_case_insensitive(cmd, "p") && !equal_case_insensitive(cmd, "l") &&
				!equal_case_insensitive(cmd, "mtllib") && !equal_case_insensitive(cmd, "usemtl") && !equal_case_insensitive(cmd, "maplib") && !equal_case_insensitive(cmd, "usemap") && 
				!equal_case_insensitive(cmd, "bevel") && !equal_case_insensitive(cmd, "c_interp") && !equal_case_insensitive(cmd, "d_interp") && !equal_case_insensitive(cmd, "log") && !equal_case_insensitive(cmd, "shadow_obj") && !equal_case_insensitive(cmd, "trace_obj"))
			{
				throw std::invalid_argument("Error: OBJ file contains unimplemented commands");
			}
			// else they are purposely ignored commands (materials, rendering tips, etc)
		}
	}

public:
	ObjReader(const char* filename) : VertexReader(filename), name("default"), as_object(false)
	{
		// Initialize
		Current cur(get_entry("", this->objs), get_entry("default", this->grps));
		this->init(&cur);

		// Cleanup empty objects and groups
		for (name2count_iter i = this->grps.begin(); i != this->grps.end(); ++i) { if (i->second == 0) { i = this->grps.erase(i); } }
		for (name2count_iter i = this->objs.begin(); i != this->objs.end(); ++i) { if (i->second == 0) { i = this->objs.erase(i); } }
	}
	inline const size_t count_objects() const { return this->objs.size(); }
	inline const size_t count_groups() const { return this->grps.size(); }
	const std::vector<std::string> object_names() const
	{
		std::vector<std::string> names;
		names.reserve(this->objs.size());
		for (name2count::const_iterator i = this->objs.begin(); i != this->objs.end(); ++i) { names.push_back(i->first); }
		return names;
	}
	const std::vector<std::string> group_names() 
	{
		std::vector<std::string> names;
		names.reserve(this->grps.size());
		for (name2count::const_iterator i = this->grps.begin(); i != this->grps.end(); ++i) { names.push_back(i->first); }
		return names;
	}
	const std::string get_obj_name(size_t i) const
	{
		size_t x = 0;
		for (name2count::const_iterator it = this->objs.begin(); it != this->objs.end(); ++it) { if (x++ == i) { return it->first; } }
		throw std::invalid_argument("Error: index out of range for objects");
	}
	const std::string get_grp_name(size_t i) const
	{
		size_t x = 0;
		for (name2count::const_iterator it = this->grps.begin(); it != this->grps.end(); ++it) { if (x++ == i) { return it->first; } }
		throw std::invalid_argument("Error: index out of range for groups");
	}
	void set_group(std::string name, bool as_object)
	{
		this->as_object = as_object;
		//if (name) { name = ltrim(name); } // rtrim doesn't work on const strings
		this->name = name = (name.empty() || name == "default") ? (as_object ? "" : "default") : name;
		if (as_object) { if (!this->objs.count(name)) { throw std::invalid_argument("Error: name doesn't exist as object"); } }
		else           { if (!this->grps.count(name)) { throw std::invalid_argument("Error: name doesn't exist as group");  } }
	}
	void operator()(HDS& hds)
	{
		// Builds the polyhedron

		// Setup basic variables
		Builder b(hds, true);
		char *cmd, *params, buf[1024];
		size_t num_vertices = 0, num_faces = 0, total_vertices = this->vertices.size();
#ifdef POLYHEDRON_CACHED_NORMALS
		size_t num_normals = 0;
#endif
		rewind(this->file);

		// Setup variables related to which group and object are to be converted
		char g_or_o;
		size_t total_faces;
		bool in_group;
		if (this->as_object) { total_faces = this->objs[this->name]; g_or_o = 'o'; in_group = this->name == ""       ; }
		else                 { total_faces = this->grps[this->name]; g_or_o = 'g'; in_group = this->name == "default"; }
		
		// Create the surface
		b.begin_surface(total_vertices, total_faces);
		for (size_t i = 0, len = this->vertices.size(); i < len; ++i)
			b.add_vertex(this->vertices[i]);
		while ((cmd = read_command(buf, ARRAYSIZE(buf), this->file, &params)) != nullptr)
		{
			char c = (char)tolower(cmd[0]);
			if (c == g_or_o)
			{
				in_group = this->as_object ? o_matches(params, this->name.c_str()) : g_matches(params, this->name.c_str());
			}
			else switch (c)
			{
			case 'v':
				// need to count these for negative/relative indexing
				if (cmd[1] == '\0') { ++num_vertices; }
#ifdef POLYHEDRON_CACHED_NORMALS
				else if (equal_case_insensitive(cmd+1, "n")) { ++num_normals; }
#endif
				break;
			case 'f':
#ifdef POLYHEDRON_CACHED_NORMALS
				if (cmd[1] == '\0' || equal_case_insensitive(cmd+1, "o")) { if (in_group) { this->f(params, &b, num_vertices, num_normals, total_vertices); if (total_faces == ++num_faces) { goto AFTER_LOOP; } } }
#else
				if (cmd[1] == '\0' || equal_case_insensitive(cmd+1, "o")) { if (in_group) { this->f(params, &b, num_vertices); if (total_faces == ++num_faces) { goto AFTER_LOOP; } } }
#endif
				break;
			}
		}

AFTER_LOOP:
		b.end_surface();

		// Remove all the vertices we added that weren't part of faces
		// This is going to be either 0 for single group files or massive for multiple-group files
		assert(b.remove_unconnected_vertices());
	}
};

struct ObjFile::internal_data
{
	ObjReader* converter;
	internal_data(ObjReader* c = nullptr) : converter(c) {}
	~internal_data() { if (this->converter) { delete this->converter; this->converter = nullptr; } }
};
ObjFile::ObjFile(const char *filename, bool as_objects, bool assume_good) : _as_objects(as_objects), _assume_good(assume_good), data(new ObjFile::internal_data(new ObjReader(filename))) { }
ObjFile::~ObjFile() { if (this->data) { delete this->data; this->data = nullptr; } }
const std::vector<std::string> ObjFile::names() const { return this->_as_objects ? this->data->converter->object_names() : this->data->converter->group_names(); }
size_t ObjFile::count() const { return this->_as_objects ? this->data->converter->count_objects() : this->data->converter->count_groups(); }

Polyhedron3* ObjFile::operator[](const std::string& name) const
{
	this->data->converter->set_group(name, this->_as_objects);
	Polyhedron3 *P = new Polyhedron3();
	P->delegate(*this->data->converter);
	if (!this->_assume_good) { check(P); }
	return P;
}
Polyhedron3* ObjFile::operator[](size_t i) const { return this->operator[](this->_as_objects ? this->data->converter->get_obj_name(i) : this->data->converter->get_grp_name(i)); }
#ifdef _MSC_VER
#pragma endregion
#endif

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
		while (i < 2 && str2int(s, &s, &v) && (!*s || isspace(*s)))
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

