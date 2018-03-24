#include "IO.hpp"
#include "IO_core.hpp"
#include "GeometryUtils.hpp"

#include "Strings.hpp"

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <map>
#include <exception>
#include <stdexcept>

#ifdef _MSC_VER
#pragma region General VertexReader
#endif
// Vertex reader class used by both ObjReader and CgReader
char* VertexReader::read_command(std::istream &f, std::string& buf, char **params)
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

void VertexReader::v(const char *params)
{
    double x, y, z, w = 1.0;
    int c1, c2, n = sscanf(params, "%lf %lf %lf%n %lf%n", &x, &y, &z, &c1, &w, &c2);
    if (n == EOF || n < 3 || (n == 3 && params[c1]) || (n == 4 && params[c2])) { throw std::invalid_argument("Error: invalid file format"); }
    this->vertices.push_back(Point(x, y, z)); // NOTE: w is dropped since it is only used in free-form (which is not supported at the moment)
}

void VertexReader::init(std::istream &in, void* extra)
{
    if (in.bad()) { throw std::invalid_argument("Error: bad file"); }

    std::string buf;
    buf.reserve(1024);
    char *cmd, *params;

    // Get all vertices and validate most of the file structure
    while ((cmd = read_command(in, buf, &params)) != nullptr)
    {
        if (cmd[0] == 0) { continue; }
        else if (strieq(cmd, "v")) { this->v(params); }
        else { this->parse_cmd(cmd, params, extra); }
    }

    // Cleanup
    this->vertices.shrink_to_fit();
}
#ifdef _MSC_VER
#pragma endregion
#endif

#ifdef _MSC_VER
#pragma region Reading/Writing CG Files
#endif
class CgReader : public VertexReader
{
private:
    void e(char *s)
    {
        // First make sure all vertices have been added to the skeleton
        for (size_t i = this->vert_handles.size(), len = this->vertices.size(); i < len; ++i)
        {
            Skeleton3::vertex_property_type v = { this->vertices[i] };
            this->vert_handles.push_back(boost::add_vertex(v, *this->skel));
        }

        ssize_t v;
        Skeleton3::vertex_descriptor vs[2];
        int i = 0;
        while (i < 2 && str_read_int(s, &s, &v) && (!*s || isspace(*s)))
        {
            if (v == 0 || (size_t)CGAL::abs(v) > this->vert_handles.size()) { break; }
            vs[i++] = this->vert_handles[(v > 0) ? v - 1 : v + this->vert_handles.size()];
            s = ltrim(s);
            if (!*s)
            {
                if (i != 2) { break; }
                boost::add_edge(vs[0], vs[1], *this->skel);
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
    std::vector<Skeleton3::vertex_descriptor> vert_handles;

public:
    CgReader(std::istream& in) : skel(new Skeleton3())
    {
        this->init(in, nullptr);
    }

    Skeleton3* skeleton() { return this->skel; }
    const Skeleton3* skeleton() const { return this->skel; }
};

Skeleton3* read_cg(const std::string& filename)
{
    std::ifstream in(filename);
    return read_cg(in);
}
Skeleton3* read_cg(std::istream& in)
{
    CgReader cg(in);
    return cg.skeleton();
}

void write_cg(const Skeleton3* S, const std::string& filename)
{
    std::ofstream out(filename);
    write_cg(S, out);
}
void write_cg(const Skeleton3* S, std::ostream& out)
{
    out << "# Skeleton written by slash-segmentation's geoEM" << std::endl << std::endl;

    // Write all of the vertices
    int i = 0;
    Skeleton3::vertex_iterator V, Vend;
    std::unordered_map<Skeleton3::vertex_descriptor, size_t> vertex_lookup;
    for (boost::tie(V, Vend) = boost::vertices(*S); V != Vend; ++V)
    {
        vertex_lookup[*V] = i++;
        Point3 v = (*S)[*V].point;
        out << "v " << v << std::endl; // v.x() << " " << v.y() << " " << v.z()
    }

    // Write all of the edges
    Skeleton3::edge_iterator E, Eend;
    for (boost::tie(E, Eend) = boost::edges(*S); E != Eend; ++E)
    {
        size_t src = (unsigned int)vertex_lookup[boost::source(*E, *S)];
        size_t trgt = (unsigned int)vertex_lookup[boost::target(*E, *S)];
        out << "e " << src << " " << trgt << std::endl;
    }
}

#ifdef _MSC_VER
#pragma endregion
#endif
