#include <iostream>
#include <iomanip>
#include <cmath>
#include "../inc/element.h"
#include "../inc/noc.h"
#include "../inc/geom.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

/**
 * Cell-to-Node interpolation coefficients.
 */
static void helper1()
{
    for (auto n_dst : node)
    {
        const size_t N = n_dst->cell_dependency.size();

        /// Allocate storage
        n_dst->cell_weighting1.resize(N);
        n_dst->cell_weighting2.resize(N);
        n_dst->cell_weighting3.resize(N);
        n_dst->cell_weighting4.resize(N);

        /// Weighting1: 1/||r||
        /// Weighting2: 1/||r||^2
        FLM_SCALAR s1 = 0.0;
        FLM_SCALAR s2 = 0.0;
        for (size_t j = 0; j < N; ++j)
        {
            auto curAdjCell = n_dst->cell_dependency.at(j);
            const FLM_SCALAR w = 1.0 / (n_dst->coordinate - curAdjCell->centroid).norm();
            n_dst->cell_weighting1.at(j) = w;
            n_dst->cell_weighting2.at(j) = w * w;
            s1 += w;
            s2 += n_dst->cell_weighting2.at(j);
        }
        for (auto &val : n_dst->cell_weighting1)
            val /= s1;
        for (auto &val : n_dst->cell_weighting2)
            val /= s2;

        /// Weighting3: 1/V
        FLM_SCALAR s3 = 0.0;
        for (size_t j = 0; j < N; ++j)
        {
            auto curAdjCell = n_dst->cell_dependency.at(j);
            const FLM_SCALAR w = 1.0 / curAdjCell->volume;
            n_dst->cell_weighting3.at(j) = w;
            s3 += w;
        }
        for (auto &val : n_dst->cell_weighting3)
            val /= s3;

        /// Weighting4: Linear preserving
        /// TODO
    }
}

/**
 * Cell centroid to face centroid vectors and ratios.
 */
static void helper2()
{
    for (auto f_dst : face)
    {
        /// Displacement vector
        if (f_dst->c0 == nullptr)
            f_dst->r0.setZero();
        else
            f_dst->r0 = f_dst->centroid - f_dst->c0->centroid;

        if (f_dst->c1 == nullptr)
            f_dst->r1.setZero();
        else
            f_dst->r1 = f_dst->centroid - f_dst->c1->centroid;

        /// Ratio
        if (f_dst->at_boundary())
        {
            if (f_dst->c0 == nullptr)
            {
                f_dst->cell_weighting1 = {0.0, 1.0};
                f_dst->cell_weighting2 = {0.0, 1.0};
                f_dst->cell_weighting3 = {0.0, 1.0};
            }
            else
            {
                f_dst->cell_weighting1 = {1.0, 0.0};
                f_dst->cell_weighting2 = {1.0, 0.0};
                f_dst->cell_weighting3 = {1.0, 0.0};
            }
        }
        else
        {
            /// Weighting1: 1/||r||
            const FLM_SCALAR rl0 = 1.0 / f_dst->r0.norm();
            const FLM_SCALAR rl1 = 1.0 / f_dst->r1.norm();
            const FLM_SCALAR s1 = rl0 + rl1;
            f_dst->cell_weighting1 = {rl0 / s1, rl1 / s1};

            /// Weighting2: 1/||r||^2
            const FLM_SCALAR rll0 = 1.0 / std::pow(f_dst->r0.norm(), 2);
            const FLM_SCALAR rll1 = 1.0 / std::pow(f_dst->r1.norm(), 2);
            const FLM_SCALAR s2 = rll0 + rll1;
            f_dst->cell_weighting2 = {rll0 / s2, rll1 / s2};

            /// Weighting3: 1/V
            const FLM_SCALAR rv0 = 1.0 / f_dst->c0->volume;
            const FLM_SCALAR rv1 = 1.0 / f_dst->c1->volume;
            const FLM_SCALAR s3 = rv0 + rv1;
            f_dst->cell_weighting3 = {rv0 / s3, rv1 / s3};
        }
    }
}

/**
 * Displacement vectors within each cell.
 */
static void helper3()
{
    for (auto c : cell)
    {
        const size_t Nf = c->surface.size();
        if (c->cell_adjacency.size() != Nf)
            throw std::runtime_error("Inconsistency detected!");

        /// Allocate storage
        c->d.resize(Nf);

        /// Vector d
        for (size_t j = 0; j < Nf; ++j)
        {
            auto cur_face = c->surface.at(j);
            auto cur_adj_cell = c->cell_adjacency.at(j);

            if (cur_adj_cell == nullptr)
                c->d.at(j) = cur_face->centroid - c->centroid;
            else
                c->d.at(j) = cur_adj_cell->centroid - c->centroid;
        }
    }
}

/**
 * Decomposition for Non-Orthogonal correction.
 */
static void helper4()
{
    for (auto c : cell)
    {
        const size_t Nf = c->surface.size();

        /// Allocate storage
        c->S_E.resize(Nf);
        c->S_T.resize(Nf);

        /// Vector S_E, S_T
        for (size_t j = 0; j < Nf; ++j)
            noc_decompose(c->d.at(j), c->S.at(j), c->S_E.at(j), c->S_T.at(j));
    }
}

void calculate_geometric_value()
{
    helper1();
    helper2();
    helper3();
    helper4();
}

static FLM_SCALAR to_degree(const FLM_SCALAR &x)
{
    static const FLM_SCALAR PI = 3.14159265;
    return x * 180.0 / PI;
}

void check_skewness()
{
    std::vector<size_t> stat(91, 0);

    for (auto f : face)
    {
        FLM_SCALAR ct;
        if (f->at_boundary())
        {
            if (f->c0 == nullptr)
                ct = f->r1.dot(f->n10) / f->r1.norm();
            else
                ct = f->r0.dot(f->n01) / f->r0.norm();
        }
        else
        {
            const FLM_VECTOR d01 = f->r0 - f->r1;
            ct = d01.dot(f->n01) / d01.norm();
        }
        f->alpha = 1.0 / ct;
        const FLM_SCALAR ang = to_degree(std::acos(ct));
        const auto tag = std::lround(ang + 0.5);
        ++stat[tag];
    }

    const auto N = face.size();

    std::cout << "==============================" << std::endl;
    std::cout << "| theta |  count  | ratio(%) |" << std::endl;
    std::cout << "------------------------------" << std::endl;
    for (size_t i = 0; i < stat.size(); ++i)
    {
        if (stat[i] == 0)
            continue;

        const FLM_SCALAR ratio = 100.0 * stat[i] / N;
        std::cout << "|" << std::setw(7) << i;
        std::cout << "|" << std::setw(9) << stat[i];
        std::cout << "|" << std::setw(10) << std::fixed << std::setprecision(4) << ratio;
        std::cout << "|" << std::endl;
    }
    std::cout << "==============================" << std::endl;
}
