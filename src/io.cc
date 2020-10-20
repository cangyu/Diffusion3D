#include "../inc/element.h"
#include "../inc/io.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

/**
 * Load computation mesh in custom format,
 * which is converted from FLUENT "msh" file.
 * @param fin Input stream of the mesh file.
 */
void read_mesh(std::istream &fin)
{
    /// Update counting of geom elements.
    size_t NumOfPatch, NumOfNode, NumOfFace, NumOfCell;
    fin >> NumOfNode >> NumOfFace >> NumOfCell >> NumOfPatch;

    /// Allocate memory for geom entities and related physical variables.
    node.resize(NumOfNode);
    for (size_t i = 1; i <= NumOfNode; ++i)
        node.at(i - 1) = new Node();

    face.resize(NumOfFace);

    cell.resize(NumOfCell);
    for (size_t i = 1; i <= NumOfCell; ++i)
        cell.at(i - 1) = new Cell();

    patch.resize(NumOfPatch);

    /// Update nodal information.
    for (size_t i = 1; i <= NumOfNode; ++i)
    {
        auto n_dst = node.at(i - 1);

        // 1-based global index
        n_dst->index = i;

        // Boundary flag
        int flag;
        fin >> flag;
        if (flag == 1)
            n_dst->at_boundary = true;
        else if (flag == 0)
            n_dst->at_boundary = false;
        else
            throw invalid_boundary_flag("node", i, flag);

        // 3D location
        fin >> n_dst->coordinate.x() >> n_dst->coordinate.y() >> n_dst->coordinate.z();

        // Adjacent nodes
        size_t n_adj_node;
        fin >> n_adj_node;
        for (size_t j = 0; j < n_adj_node; ++j)
        {
            size_t tmp;
            fin >> tmp;
        }

        // Dependent faces
        size_t n_dep_face;
        fin >> n_dep_face;
        for (size_t j = 0; j < n_dep_face; ++j)
        {
            size_t tmp;
            fin >> tmp;
        }

        // Dependent cells
        size_t n_dep_cell;
        fin >> n_dep_cell;
        n_dst->cell_dependency.resize(n_dep_cell);
        for (size_t j = 0; j < n_dep_cell; ++j)
        {
            size_t tmp;
            fin >> tmp;
            n_dst->cell_dependency.at(j) = cell.at(tmp - 1);
        }
    }

    /// Update face information.
    for (size_t i = 1; i <= NumOfFace; ++i)
    {
        // Boundary flag
        int flag;
        fin >> flag;
        if (flag == 1)
            face.at(i - 1) = new BoundaryFace();
        else if (flag == 0)
            face.at(i - 1) = new InternalFace();
        else
            throw invalid_boundary_flag("face", i, flag);

        auto f_dst = face.at(i - 1);

        // 1-based global index
        f_dst->index = i;

        // Shape
        int shape;
        fin >> shape;
        if (shape != 3 && shape != 4)
            throw unsupported_shape("face", i, shape);

        // Centroid
        fin >> f_dst->centroid.x() >> f_dst->centroid.y() >> f_dst->centroid.z();

        // Area
        fin >> f_dst->area;

        // Included nodes
        f_dst->vertex.resize(shape);
        for (int j = 0; j < shape; ++j)
        {
            size_t tmp;
            fin >> tmp;
            f_dst->vertex.at(j) = node.at(tmp - 1);
        }

        // Adjacent cells
        size_t c0, c1;
        fin >> c0 >> c1;
        if (c0 == 0)
            f_dst->cell_dependency[0] = nullptr;
        else
            f_dst->cell_dependency[0] = cell.at(c0 - 1);

        if (c1 == 0)
            f_dst->cell_dependency[1] = nullptr;
        else
            f_dst->cell_dependency[1] = cell.at(c1 - 1);

        // Unit normal vector
        fin >> f_dst->n[0].x() >> f_dst->n[0].y() >> f_dst->n[0].z();
        fin >> f_dst->n[1].x() >> f_dst->n[1].y() >> f_dst->n[1].z();
    }

    /// Update cell information.
    for (size_t i = 1; i <= NumOfCell; ++i)
    {
        auto c_dst = cell.at(i - 1);

        // 1-based global index
        c_dst->index = i;

        // Shape
        int shape;
        fin >> shape;
        int N1, N2;
        if (shape == 2) // tet
        {
            N1 = 4;
            N2 = 4;
        }
        else if (shape == 4) // hex
        {
            N1 = 8;
            N2 = 6;
        }
        else if (shape == 5) // pyr
        {
            N1 = 5;
            N2 = 5;
        }
        else if (shape == 6) // wedge
        {
            N1 = 6;
            N2 = 5;
        }
        else
            throw unsupported_shape("cell", i, shape);

        // Centroid
        fin >> c_dst->centroid.x() >> c_dst->centroid.y() >> c_dst->centroid.z();

        // Volume
        fin >> c_dst->volume;

        // Included nodes
        c_dst->vertex.resize(N1);
        for (int j = 0; j < N1; ++j)
        {
            size_t tmp;
            fin >> tmp;
            c_dst->vertex.at(j) = node.at(tmp - 1);
        }

        // Included faces
        c_dst->surface.resize(N2);
        for (int j = 0; j < N2; ++j)
        {
            size_t tmp;
            fin >> tmp;
            c_dst->surface.at(j) = face.at(tmp - 1);
        }

        // Adjacent cells
        c_dst->cell_adjacency.resize(N2);
        for (int j = 0; j < N2; ++j)
        {
            size_t tmp;
            fin >> tmp;
            if (tmp == 0)
                c_dst->cell_adjacency.at(j) = nullptr;
            else
                c_dst->cell_adjacency.at(j) = cell.at(tmp - 1);
        }

        // Surface unit normal vectors
        c_dst->n.resize(N2);
        for (int j = 0; j < N2; ++j)
        {
            auto &loc_n = c_dst->n.at(j);
            fin >> loc_n.x() >> loc_n.y() >> loc_n.z();
        }
    }

    /// Update boundary patch information.
    for (size_t i = 1; i <= NumOfPatch; ++i)
    {
        auto p_dst = patch.at(i - 1) = new Patch();
        fin >> p_dst->name;

        size_t n_face, n_node;
        fin >> n_face >> n_node;
        p_dst->surface.resize(n_face);
        p_dst->vertex.resize(n_node);
        for (size_t j = 0; j < n_face; ++j)
        {
            size_t tmp;
            fin >> tmp;
            auto ptr = dynamic_cast<BoundaryFace *>(face.at(tmp - 1));
            if (ptr == nullptr)
                throw wrong_face(tmp);

            p_dst->surface.at(j) = ptr;
            ptr->set_parent(p_dst);
        }
        for (size_t j = 0; j < n_node; ++j)
        {
            size_t tmp;
            fin >> tmp;
            p_dst->vertex.at(j) = node.at(tmp - 1);
        }
    }
}

void write_data(std::ostream &out, size_t iter, FLM_SCALAR t)
{
    static const char SEP = ' ';

    out << iter << SEP << t << std::endl;
    out << node.size() << SEP << face.size() << SEP << cell.size() << std::endl;

    for (auto e : node)
    {
        out << e->T << std::endl;
    }

    for (auto e : face)
    {
        out << e->T << std::endl;
    }

    for (auto e : cell)
    {
        out << e->T << std::endl;
    }
}

void read_data(std::istream &in, size_t &iter, FLM_SCALAR &t)
{
    in >> iter >> t;

    size_t n_node, n_face, n_cell;
    in >> n_node >> n_face >> n_cell;

    if (n_node != node.size() || n_face != face.size() || n_cell != cell.size())
        throw inconsistent_mesh();

    for (auto e : node)
    {
        in >> e->T;
    }

    for (auto e : face)
    {
        in >> e->T;
    }

    for (auto e : cell)
    {
        in >> e->T;
    }
}
