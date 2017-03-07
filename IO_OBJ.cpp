#include "IO_OBJ.hpp"
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
    
    static char* read_command(std::istream &f, std::string& buf, char **params)
    {
        if (f.eof()) { return nullptr; }

        // Read an entire line
        std::getline(f, buf);
        size_t comment = buf.find_first_of('#');
        if (comment != std::string::npos) { buf.erase(comment); }
        char *cmd = trim(const_cast<char*>(buf.c_str()));

        // Get the command and parameters
        char *space = strchr(cmd, " \t\v");
        if (space) { *space = '\0'; *params = ltrim(space + 1); }
        else { *params = cmd + buf.length(); } // no params
        return cmd; // the command
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
    
    std::vector<Point> vertices;
#ifdef POLYHEDRON_CACHED_NORMALS
    std::vector<Direction> normals;
#endif

    virtual void parse_cmd(const char* cmd, char* params, void* extra) = 0;
    void init(std::istream &in, void* extra)
    {
        if (in.bad()) { throw std::invalid_argument("Error: bad file"); }

        std::string buf;
        buf.reserve(1024);
        char *cmd, *params;

        // Get all vertices and validate most of the file structure
        while ((cmd = read_command(in, buf, &params)) != nullptr)
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
};

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
}
#ifdef _MSC_VER
#pragma endregion
#endif

#ifdef _MSC_VER
#pragma region Reading OBJ Files
#endif
class ObjFileDetail : public VertexReader, public CGAL::Modifier_base<Polyhedron3::HalfedgeDS>
{
private:
    typedef Polyhedron3::HalfedgeDS HDS;
    typedef HDS::Vertex   Vertex;
    typedef CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    
    typedef std::vector<size_t> facet;
    struct raw_mesh
    {
        std::string name;
        std::vector<size_t> facets;
    };
    typedef std::list<raw_mesh> raw_meshes;
    typedef std::map<std::string, raw_meshes::iterator> name_lookup;
    
    inline static raw_meshes::iterator get_entry(const char *name, name_lookup &lookup, raw_meshes &meshes)
    {
        name_lookup::const_iterator i = lookup.find(name);
        if (i != lookup.end()) { return i->second; }
        raw_meshes::iterator x = meshes.insert(meshes.end(), raw_mesh());
        lookup.insert(i, name_lookup::value_type(name, x));
        return x;
    }
    std::vector<raw_meshes::iterator> g(char *params)
    {
        std::vector<raw_meshes::iterator> grps;
        char *name = strtok(params, " \t\v");
        while (name != nullptr)
        {
            grps.push_back(get_entry(name, this->grp_lookup, this->grps));
            name = strtok(nullptr, " \t\v");
        }
        if (grps.size() == 0) { grps.push_back(get_entry("default", this->grp_lookup, this->grps)); }
        return grps;
    }
    raw_meshes::iterator o(const char *params) { return get_entry(params, this->obj_lookup, this->objs); }
    
#ifdef POLYHEDRON_CACHED_NORMALS
    facet f_v(char *s, bool has_vt)
    {
        facet f;
        ssize_t v;
        while (str_read_int(s, &s, &v) && // read "v"
                (!has_vt || *s == '/' && str_skip_int(s+1, &s)) && // skip "vt"
                (!*s || isspace(*s)))
        {
            // Add vertex to facet
            if (v == 0 || (size_t)CGAL::abs(v) > this->vertices.size()) { break; }
            f.push_back((v > 0) ? v - 1 : v + this->vertices.size());
            s = ltrim(s);
            if (!*s) { if (f.size() < 3) { break; } return f; } // Finished
        }
        throw std::invalid_argument("Error: invalid OBJ file format");
    }
    facet f_v_vn(char *s, bool has_vt)
    {
        // TODO: update for new reading system
        b->begin_facet();
        Builder::Vertex_handle vv;
        ssize_t v, n;
        size_t count = 0;
        while (str_read_int(s, &s, &v) && // read "v"
                *s++ == '/' && (!has_vt || str_skip_int(s, &s)) && // skip "vt"
                *s == '/' && str_read_int(s+1, &s, &n) && (!*s || isspace(*s))) // read "vn"
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
        throw std::invalid_argument("Error: invalid OBJ file format");
    }
    facet f(char *params)
    {
        // At least 3 points
        // Each point is of one of the following forms, and in a face all points have the same form
        //   v, v/vt, v//vn, v/vt/vn

        // Call the proper function for the given form by parsing the first vertex
        char *s = params;
        if (str_skip_int(s, &s))
        {
            if (isspace(*s)) { return f_v(params, false); } // v
            else if (*s++ == '/')
            {
                if (*s == '/' && str_skip_int(s+1, &s)) { return f_v_vn(params, false); } // v//vn
                else if (str_skip_int(s, &s))
                {
                    if (isspace(*s)) { return f_v(params, true); } // v/vt
                    else if (*s == '/' && str_skip_int(s+1, &s) && isspace(*s)) { return f_v_vn(params, true); } // v/vt/vn
                }
            }
        }
        throw std::invalid_argument("Error: invalid OBJ file format");
    }
#else
    facet f_v(char *s, bool has_vt, bool has_vn)
    {
        facet f;
        ssize_t v;
        while (str_read_int(s, &s, &v) && // read "v"
                ((!has_vt && !has_vn) || *s++ == '/') &&
                (!has_vt || str_skip_int(s, &s)) && // skip "vt"
                (!has_vn || (*s == '/' && str_skip_int(s+1, &s))) && // skip "vn"
                (!*s || isspace(*s)))
        {
            // Add vertex to facet
            if (v == 0 || (size_t)CGAL::abs(v) > this->vertices.size()) { break; }
            f.push_back((v > 0) ? v - 1 : v + this->vertices.size());
            s = ltrim(s);
            if (!*s) { if (f.size() < 3) { break; } return f; } // Finished
        }
        throw std::invalid_argument("Error: invalid OBJ file format");
    }
    facet f(char *params)
    {
        // At least 3 points
        // Each point is of one of the following forms, and in a face all points have the same form
        //   v, v/vt, v//vn, v/vt/vn

        // Call the proper function for the given form by parsing the first vertex
        char *s = params;
        if (str_skip_int(s, &s))
        {
            if (isspace(*s)) { return f_v(params, false, false); } // v
            else if (*s++ == '/')
            {
                if (*s == '/' && str_skip_int(s+1, &s)) { return f_v(params, false, true); } // v//vn
                else if (str_skip_int(s, &s))
                {
                    if (isspace(*s)) { return f_v(params, true, false); } // v/vt
                    else if (*s == '/' && str_skip_int(s+1, &s) && isspace(*s)) { return f_v(params, true, true); } // v/vt/vn
                }
            }
        }
        throw std::invalid_argument("Error: invalid OBJ file format");
    }
#endif

    inline static bool o_matches(const char *params, const char* obj_name) { return strcmp(params, obj_name) == 0; }
    inline static bool g_matches(char *params, const char* grp_name) { char *name = strtok(params, " "); while (name != nullptr) { if (strcmp(name, grp_name) == 0) { return true; } name = strtok(nullptr, " "); } return false; }
    
    std::vector<facet> facets;
    raw_meshes objs;
    raw_meshes grps;
    name_lookup obj_lookup;
    name_lookup grp_lookup;
    raw_meshes::iterator current;

    struct Current
    {
        raw_meshes::iterator obj;
        std::vector<raw_meshes::iterator> grps;
        Current(raw_meshes::iterator o, raw_meshes::iterator g) : obj(o), grps(1, g) {}
    };

    virtual void parse_cmd(const char* cmd, char* params, void* extra)
    {
        Current* cur = (Current*)extra;
        switch (tolower(cmd[0]))
        {
        // v and vn are already done, need to read the faces into each group/object and validate ignored commands
        case 'o':
            if (cmd[1] == '\0') { cur->obj = o(params); }
            else { throw std::invalid_argument("Error: OBJ file contains unimplemented commands"); }
            break;
        case 'g':
            if (cmd[1] == '\0') { cur->grps = g(params); }
            else { throw std::invalid_argument("Error: OBJ file contains unimplemented commands"); }
            break;
        case 'f':
            if (cmd[1] == '\0' || streq_case_insensitive(cmd+1, "o")) // TODO: validation
            {
                facet f = this->f(params);
                size_t fi = this->facets.size();
                this->facets.push_back(f);
                cur->obj->facets.push_back(fi);
                for (size_t i = 0, len = cur->grps.size(); i < len; ++i)
                    cur->grps[i]->facets.push_back(fi);
            }
            else { throw std::invalid_argument("Error: OBJ file contains unimplemented commands"); }
            break;
        default:
            if (!streq_case_insensitive(cmd, "vt") && !streq_case_insensitive(cmd, "s") && !streq_case_insensitive(cmd, "p") && !streq_case_insensitive(cmd, "l") &&
                !streq_case_insensitive(cmd, "mtllib") && !streq_case_insensitive(cmd, "usemtl") && !streq_case_insensitive(cmd, "maplib") && !streq_case_insensitive(cmd, "usemap") && 
                !streq_case_insensitive(cmd, "bevel") && !streq_case_insensitive(cmd, "c_interp") && !streq_case_insensitive(cmd, "d_interp") && !streq_case_insensitive(cmd, "log") && !streq_case_insensitive(cmd, "shadow_obj") && !streq_case_insensitive(cmd, "trace_obj"))
            {
                throw std::invalid_argument("Error: OBJ file contains unimplemented commands");
            }
            // else they are purposely ignored commands (materials, rendering tips, etc)
        }
    }

public:
    ObjFileDetail(std::istream& in)
    {
        // Initialize
        Current cur(get_entry("", this->obj_lookup, this->objs), get_entry("default", this->grp_lookup, this->grps));
        this->init(in, &cur);

        // Remove default names if they are empty
        if (this->objs.front().facets.size() == 0) { this->objs.pop_front(); this->obj_lookup.erase(""); }
        if (this->grps.front().facets.size() == 0) { this->grps.pop_front(); this->grp_lookup.erase("default"); }

        // Remove excess allocated data
        for (auto i = this->objs.begin(), end = this->objs.end(); i != end; ++i) { i->facets.shrink_to_fit(); }
        for (auto i = this->grps.begin(), end = this->grps.end(); i != end; ++i) { i->facets.shrink_to_fit(); }
        this->facets.shrink_to_fit();
    }
    inline const size_t count_objects() const { return this->objs.size(); }
    inline const size_t count_groups()  const { return this->grps.size(); }
    const std::vector<std::string> object_names() const
    {
        std::vector<std::string> names;
        names.reserve(this->objs.size());
        for (auto i = this->obj_lookup.begin(), end = this->obj_lookup.end(); i != end; ++i) { names.push_back(i->first); }
        return names;
    }
    const std::vector<std::string> group_names() const
    {
        std::vector<std::string> names;
        names.reserve(this->grps.size());
        for (auto i = this->grp_lookup.begin(), end = this->grp_lookup.end(); i != end; ++i) { names.push_back(i->first); }
        return names;
    }
    void set_current(size_t i, bool as_object)
    {
        raw_meshes& rm = as_object ? this->objs : this->grps;
        if (i >= rm.size()) { std::invalid_argument("Error: index out of range in OBJ file"); }
        raw_meshes::iterator x = rm.begin();
        std::advance(x, i);
        this->current = x;
    }
    void set_current(const std::string& name, bool as_object)
    {
        name_lookup& nl = as_object ? this->obj_lookup : this->grp_lookup;
        name_lookup::iterator x = nl.find((name.empty() || name == "default") ? (as_object ? "" : "default") : name);
        if (x == nl.end()) { throw std::invalid_argument("Error: name doesn't exist in OBJ file"); }
        this->current = x->second;
    }
    void operator()(HDS& hds)
    {
        // Builds the polyhedron

        // Setup basic variables
        Builder b(hds, true);
        size_t total_faces = this->current->facets.size(), total_vertices = 0;

        // Calculate the vertices and their new indices
        for (std::vector<size_t>::const_iterator i = this->current->facets.begin(), end = this->current->facets.end(); i != end; ++i)
        {
            total_vertices += this->facets[*i].size();
        }
        
        // Create the surface
        // TODO: need to calculate the actual set of vertices we want to keep in the mesh and the new vertex offsets
        b.begin_surface(total_vertices, total_faces);
        for (size_t i = 0, len = this->vertices.size(); i < len; ++i)
        {
            b.add_vertex(this->vertices[i]);
        }
        for (std::vector<size_t>::const_iterator i = this->current->facets.begin(), end = this->current->facets.end(); i != end; ++i)
        {
            b.begin_facet();
            facet& f = this->facets[*i];
            for (size_t j = 0; j < f.size(); ++j) { b.add_vertex_to_facet(f[j]); }
            b.end_facet();
        }
        b.end_surface();

#ifndef POLYHEDRON_USE_VECTOR
        // Remove all the vertices we added that weren't part of faces
        // This is going to be either 0 for single group files or massive for multiple-group files
        assert(b.remove_unconnected_vertices());
#endif
    }
};

ObjFile::ObjFile(const char *filename, bool as_objects) : _as_objects(as_objects), detail(new ObjFileDetail(std::ifstream(filename, std::ifstream::binary))) { }
ObjFile::ObjFile(std::istream& in, bool as_objects) : _as_objects(as_objects), detail(new ObjFileDetail(in)) { }
ObjFile::~ObjFile() { if (this->detail) { delete this->detail; this->detail = nullptr; } }
const std::vector<std::string> ObjFile::names() const { return this->_as_objects ? this->detail->object_names() : this->detail->group_names(); }
size_t ObjFile::size() const { return this->_as_objects ? this->detail->count_objects() : this->detail->count_groups(); }
Polyhedron3* ObjFile::operator[](const std::string& name) const
{
    this->detail->set_current(name, this->_as_objects);
    Polyhedron3 *P = new Polyhedron3();
    P->delegate(*this->detail);
    return P;
}
Polyhedron3* ObjFile::operator[](size_t i) const
{
    this->detail->set_current(i, this->_as_objects);
    Polyhedron3 *P = new Polyhedron3();
    P->delegate(*this->detail);
    return P;
}
#ifdef _MSC_VER
#pragma endregion
#endif