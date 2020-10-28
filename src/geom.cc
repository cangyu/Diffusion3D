#include "../inc/element.h"
#include "../inc/noc.h"
#include "../inc/geom.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

void calc_geom_val()
{
    /// Cell-to-Node interpolation coefficients.
    for (auto n_dst : node)
    {
        FLM_SCALAR s = 0.0;
        for (size_t j = 0; j < n_dst->cell_weighting.size(); ++j)
        {
            auto curAdjCell = n_dst->cell_dependency.at(j);
            const FLM_SCALAR weighting = 1.0 / (n_dst->coordinate - curAdjCell->centroid).norm();
            n_dst->cell_weighting.at(j) = weighting;
            s += weighting;
        }
        for (auto &val : n_dst->cell_weighting)
            val /= s;
    }

    /// Cell centroid to face centroid vectors and ratios.
    for (auto f_dst : face)
    {
        /// Displacement vector
        if (f_dst->c0)
            f_dst->r0 = f_dst->centroid - f_dst->c0->centroid;
        else
            f_dst->r0.setZero();

        if (f_dst->c1)
            f_dst->r1 = f_dst->centroid - f_dst->c1->centroid;
        else
            f_dst->r1.setZero();

        /// Displacement ratio
        if (f_dst->at_boundary())
        {
            if (f_dst->c0)
            {
                f_dst->ksi0 = 1.0;
                f_dst->ksi1 = 0.0;
                f_dst->g0 = 1.0;
                f_dst->g1 = 0.0;
            }
            else
            {
                f_dst->ksi0 = 0.0;
                f_dst->ksi1 = 1.0;
                f_dst->g0 = 0.0;
                f_dst->g1 = 1.0;
            }
        }
        else
        {
            const FLM_SCALAR l0 = 1.0 / f_dst->r0.norm();
            const FLM_SCALAR l1 = 1.0 / f_dst->r1.norm();
            const FLM_SCALAR w = l0 + l1;
            f_dst->ksi0 = l0 / w;
            f_dst->ksi1 = l1 / w;
            f_dst->g0 = 1.0 / f_dst->c0->volume;
            f_dst->g1 = 1.0 / f_dst->c1->volume;
        }
    }

    /// Vectors used for Non-Orthogonal correction within each cell.
    for (auto cur_cell : cell)
    {
        const size_t Nf = cur_cell->surface.size();
        if (cur_cell->adjCell.size() != Nf)
            throw std::runtime_error("Inconsistency detected!");

        // Allocate storage
        cur_cell->S_E.resize(Nf);
        cur_cell->S_T.resize(Nf);
        cur_cell->d.resize(Nf);
        cur_cell->d_mod.resize(Nf);

        // Vector d, E, T
        for (size_t j = 0; j < Nf; ++j)
        {
            auto cur_face = cur_cell->surface.at(j);
            auto cur_adj_cell = cur_cell->adjCell.at(j);

            // Displacement vector
            auto &cur_d = cur_cell->d.at(j);
            if (cur_adj_cell == nullptr)
                cur_d = cur_face->centroid - cur_cell->centroid;
            else
                cur_d = cur_adj_cell->centroid - cur_cell->centroid;

            // Non-Orthogonal correction
            noc_decompose(cur_d, cur_cell->S.at(j), cur_cell->S_E.at(j), cur_cell->S_T.at(j));
        }

        // ||d||
        for (size_t j = 0; j < Nf; ++j)
            cur_cell->d_mod.at(j) = cur_cell->d.at(j).norm();
    }
}
