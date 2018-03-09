#include "IO.hpp"

#include "IO_OBJ.hpp"

#include "GeometryTypes.hpp"
#include "GeometryUtils.hpp"
#include "Polyhedron3Utils.hpp"
#include "TriangulatePolyhedron.hpp"

#include "Strings.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const std::string mesh_file_usage = 
    "Supported polyhedral mesh formats determined from file extension:\n"
    "  OBJ - each file can contain multiple, possibly named, meshes\n"
    "  OFF - specify a single 3D mesh each\n"
    "In both formats normals, colors, textures, etc are all ignored. Multi-object OFF\n"
    "files are not supported."
    "\n\n"
    "Filenames can also include options for reading/writing the files. They are given\n"
    "after a colon (:) after the filename as a comma-separated list of name=value\n"
    "items (or just names for boolean options). The names are case-insensitive.\n";

const std::string mesh_file_read_usage = 
    "Read Options:\n"
    "  OBJ - ObjectName  - load the given object by name\n"
    "        GroupName   - load the given group by name\n"
    "        ObjectIndex - load the given object by index\n"
    "        GroupIndex  - load the given group by index\n"
    "        Default is to load the first group, GroupIndex=0\n"
    "  OFF - (no options)\n";

const std::string mesh_file_write_usage = 
    "Writing Options:\n"
    "  OBJ - ObjectName, GroupName, ObjectIndex, GroupIndex as per reading options\n"
    "        Append      - if given data is appended to existing file\n"
    "        Default is to create a new file with a single unnamed group\n"
    "        If Append is given, ObjectName/GroupName must be given as well\n"
    "  OFF - Binary      - if given the file in written in the binary OFF format\n";


void check_mesh(Polyhedron3* P)
{
    // expensive operations that won't be necessary, usually
    // TODO: add check for isolated vertices
    if (!P->is_closed())       { throw std::invalid_argument("Error: non-closed polyhedron."); }
    if (!is_not_degenerate(P)) { throw std::invalid_argument("Error: degenerate polyhedron."); }
    if (!is_manifold(P))       { throw std::invalid_argument("Error: non-manifold polyhedron."); }
    if (!P->is_pure_triangle())
    {
#ifdef POLYHEDRON_USE_VECTOR
        throw std::invalid_argument("Error: Polyhedron has non-triangle facets and unable to triangulate");
#else
        // TODO: add check for planarity of faces?
        std::cerr << "Polyhedron has non-triangle facets, triangulating..." << std::endl;
        triangulate_polyhedron(P);
#endif
    }
}

#ifdef _MSC_VER
#pragma region General filename and options parsing
#endif
file_options parse_filename_and_options(const std::string& fao, std::string& filename)
{
    const static int start = 
#ifdef _WIN32
        2; // avoid drive paths like C:
#else
        1;
#endif
    size_t c = fao.find_first_of(':', start);
    filename = fao.substr(0, c);
    return (c == std::string::npos) ? file_options() : parse_file_options(fao.substr(c+1));
}
file_options parse_file_options(const std::string& options)
{
    std::stringstream opts_raw(options);
    file_options opts;
    std::string o;
    while (std::getline(opts_raw, o, ','))
    {
        size_t eq = o.find_first_of('=');
        opts[tolower(trim(o.substr(0, eq)))] = (eq == std::string::npos) ? "" : trim(o.substr(eq+1));
    }
    return opts;
}
file_type get_file_type(const std::string& filename)
{
    size_t dot = filename.find_last_of('.');
    if (dot == std::string::npos) { return file_type::UNKNOWN; }
    std::string ext = tolower(filename.substr(dot+1));
    if      (ext == "obj") { return file_type::OBJ; }
    else if (ext == "off") { return file_type::OFF; }
    return file_type::UNKNOWN;
}
inline bool get_bool_opt(const std::string& value)
{
    return !((value.length() == 1 && (value[0] == '0' || value[0] == 'f' || value[0] == 'F')) ||
             (value.length() == 5 && tolower(value) == "false"));
}
#ifdef _MSC_VER
#pragma endregion
#endif

#ifdef _MSC_VER
#pragma region Reading
#endif

class P3Reader
{
public:
    virtual Polyhedron3* operator()(std::istream &in) = 0;
    virtual ~P3Reader() { }
};
class P3ReaderFunc : public P3Reader
{
public:
    typedef Polyhedron3* (*Func)(std::istream &in);
    P3ReaderFunc(Func f) : f(f) { }
    virtual Polyhedron3* operator()(std::istream &in) { return this->f(in); }
private:
    Func f;
};

class P3Reader_OBJ : public P3Reader
{
public:
    P3Reader_OBJ(bool as_objects, size_t i) : from_name(false), as_objects(as_objects), i(i), name("") { }
    P3Reader_OBJ(bool as_objects, const std::string& name) : from_name(true), as_objects(as_objects), i(0), name(name) { }
    virtual Polyhedron3* operator()(std::istream &in) 
    {
        return this->from_name ? ObjFile(in, this->as_objects)[this->name] : ObjFile(in, this->as_objects)[this->i];
    }
private:
    const bool from_name, as_objects;
    const size_t i;
    const std::string name;
};
static P3Reader* get_obj_reader(const file_options& options)
{
    if (options.size() > 1) { throw std::invalid_argument("Error: invalid file options given for reading mesh"); }
    if (options.size() == 0) { return new P3Reader_OBJ(false, 0); }
    auto& entry = *options.begin();
    const std::string& name = entry.second;
    if      (entry.first == "objectname")  { return new P3Reader_OBJ(true,  name); }
    else if (entry.first == "groupname")   { return new P3Reader_OBJ(false, name); }
    else if (entry.first == "objectindex") { size_t i = read_int(name); return new P3Reader_OBJ(true,  i); }
    else if (entry.first == "groupindex")  { size_t i = read_int(name); return new P3Reader_OBJ(false, i); }
    else { throw std::invalid_argument("Error: invalid file options given for reading mesh"); }
}

class P3Reader_OFF : public P3Reader
{
public:
    virtual Polyhedron3* operator()(std::istream &in)
    { 
        Polyhedron3* P = new Polyhedron3();
        in >> *P;
        return P;
    }
};
static P3Reader* get_off_reader(const file_options& options)
{
    if (options.size()) { throw std::invalid_argument("Error: invalid file options given for reading mesh"); }
    return new P3Reader_OFF();
}


Polyhedron3* read_mesh(const std::string& filename_and_options, bool assume_good)
{
    std::string filename;
    file_options opts = parse_filename_and_options(filename_and_options, filename);
    return read_mesh(filename, opts, assume_good);
}
Polyhedron3* read_mesh(const std::string& filename, const std::string& options, bool assume_good)
{
    return read_mesh(filename, parse_file_options(options), assume_good);
}
Polyhedron3* read_mesh(const std::string& filename, const file_options& options, bool assume_good)
{
    std::ifstream f(filename.c_str());
    if (f.bad()) { throw std::invalid_argument("Error: unable to open file for reading"); }
    return read_mesh(f, get_file_type(filename), options, assume_good);
}
Polyhedron3* read_mesh(std::istream &in, file_type type, const file_options& options, bool assume_good)
{
    if (in.bad()) { throw std::invalid_argument("Error: file handle is bad"); }
    std::unique_ptr<P3Reader> reader;
    switch (type)
    {
    case file_type::OBJ: reader.reset(get_obj_reader(options)); break;
    case file_type::OFF: reader.reset(get_off_reader(options)); break;
    default: throw std::invalid_argument("Error: invalid file type given the reading mesh");
    }
    Polyhedron3* P  = reader->operator()(in);
    if (!assume_good) { check_mesh(P); }
    return P;
}
#ifdef _MSC_VER
#pragma endregion
#endif

#ifdef _MSC_VER
#pragma region Writing
#endif

class P3Writer
{
public:
    virtual void operator()(const Polyhedron3* P, std::ostream &out) = 0;
    virtual ~P3Writer() { }
};
class P3WriterFunc : public P3Writer
{
public:
    typedef void (*Func)(const  Polyhedron3* P, std::ostream &out);
    P3WriterFunc(Func f) : f(f) { }
    virtual void operator()(const Polyhedron3* P, std::ostream &out) { return this->f(P, out); }
private:
    Func f;
};

class P3Writer_OBJ : public P3Writer
{
public:
    virtual void operator()(const Polyhedron3* P, std::ostream &out)
    {
        IteratorReverseLookup<Polyhedron3::Vertex_const_handle> vertex_lookup(P->vertices_begin(), P->size_of_vertices());
        CGAL::set_ascii_mode(out);
        for (Polyhedron3::Point_const_iterator i = P->points_begin(), end = P->points_end(); i != end; ++i) { out << "v " << *i << std::endl; }
        for (Polyhedron3::Facet_const_iterator i = P->facets_begin(), end = P->facets_end(); i != end; ++i)
        {
            out << 'f';
            FOR_VERTICES_AROUND_FACET(i, v) { out << ' ' << (vertex_lookup[v] + 1 /*+ off*/); }
            out << std::endl;
        }
        //off += P->size_of_vertices();
    }
};
static P3Writer* get_obj_writer(const file_options& options)
{
    return NULL;
}


class P3Writer_OFF : public P3Writer
{
    bool binary;
public:
    P3Writer_OFF(bool binary) : binary(binary) { }
    virtual void operator()(const Polyhedron3* P, std::ostream &out)
    {
        if (this->binary) { CGAL::set_binary_mode(out); }
        else
        {
            CGAL::set_ascii_mode(out);
            //if (this->verbose) { CGAL::set_pretty_mode(out); }
        }
        out << *P;
    }
};
static P3Writer* get_off_writer(const file_options& options)
{
    if (options.size() > 1) { throw std::invalid_argument("Error: invalid file options given for writing mesh"); }
    if (options.size() == 0) { return new P3Writer_OFF(false); }
    auto& entry = *options.begin();
    if (entry.first == "binary") { return new P3Writer_OFF(get_bool_opt(entry.second)); }
    else { throw std::invalid_argument("Error: invalid file options given for writing mesh"); }
}

void write_mesh(const Polyhedron3* P, const std::string& filename_and_options)
{
    std::string filename;
    file_options opts = parse_filename_and_options(filename_and_options, filename);
    write_mesh(P, filename, opts);
}
void write_mesh(const Polyhedron3* P, const std::string& filename, const std::string& options)
{
    write_mesh(P, filename, parse_file_options(options));
}
void write_mesh(const Polyhedron3* P, const std::string& filename, const file_options& options)
{
    std::ofstream f(filename.c_str());
    if (f.bad()) { throw std::invalid_argument("Error: unable to open file for writing"); }
    write_mesh(P, f, get_file_type(filename), options);
}
void write_mesh(const Polyhedron3* P, std::ostream &out, file_type type, const file_options& options)
{
    if (out.bad()) { throw std::invalid_argument("Error: file handle is bad"); }
    std::unique_ptr<P3Writer> writer;
    switch (type)
    {
    case file_type::OBJ: writer.reset(get_obj_writer(options)); break;
    case file_type::OFF: writer.reset(get_off_writer(options)); break;
    default: throw std::invalid_argument("Error: invalid file type given the writing mesh");
    }
    writer->operator()(P, out);
}

#ifdef _MSC_VER
#pragma endregion
#endif
