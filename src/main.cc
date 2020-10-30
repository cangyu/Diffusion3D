#include <iostream>
#include <fstream>
#include <cstdlib>
#include <regex>
#include <filesystem>
#include "../inc/element.h"
#include "../inc/io.h"
#include "../inc/geom.h"
#include "../inc/bc.h"
#include "../inc/lsq.h"
#include "../inc/ic.h"
#include "../inc/spatial.h"
#include "../inc/temporal.h"
#include "../inc/diagnose.h"
#include "../inc/misc.h"

std::vector<Patch *> patch;
std::vector<Node *> node;
std::vector<Face *> face;
std::vector<Cell *> cell;

/// Iteration timing and counting
static size_t OUTPUT_GAP = 100;
static size_t MAX_ITER = 1000000;
static FLM_SCALAR MAX_TIME = 10000.0; /// s
static size_t iter = 0;
static FLM_SCALAR t = 0.0; /// s
FLM_SCALAR dt = 1e-4; /// s

static void banner()
{
    std::cout << "================================================================================" << std::endl;
    std::cout << "                                  Diffusion3D                                   " << std::endl;
    std::cout << "            Solve 3D Poisson equation using FVM on unstructured mesh.           " << std::endl;
    std::cout << "                        (2nd-Order, Sequential, Steady)                         " << std::endl;
    std::cout << "================================================================================" << std::endl;
}

int main(int argc, char *argv[])
{
    std::string MESH_PATH, DATA_PATH, RUN_TAG;
    std::string OUTPUT_PREFIX = "ITER";
    bool resume_mode = false;
    clock_t tick_begin, tick_end;

    banner();

    /// Parse parameters
    int cnt = 1;
    while (cnt < argc)
    {
        if (!std::strcmp(argv[cnt], "--mesh"))
        {
            MESH_PATH = argv[cnt + 1];
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--data"))
        {
            /// Will initialize from certain data file
            DATA_PATH = argv[cnt + 1];
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--tag"))
        {
            RUN_TAG = argv[cnt + 1];
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--iteration"))
        {
            char *pEnd;
            MAX_ITER = std::strtol(argv[cnt + 1], &pEnd, 10);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--time-span"))
        {
            char *pEnd;
            MAX_TIME = std::strtod(argv[cnt + 1], &pEnd); /// s
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--time-step"))
        {
            char *pEnd;
            dt = std::strtod(argv[cnt + 1], &pEnd); /// s
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--write-interval"))
        {
            char *pEnd;
            OUTPUT_GAP = std::strtol(argv[cnt + 1], &pEnd, 10);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--resume-from"))
        {
            RUN_TAG = argv[cnt + 1];
            resume_mode = true;
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--help"))
        {
            /// TODO
            return 0;
        }
        else if (!std::strcmp(argv[cnt], "--version"))
        {
            std::cout << "V2.0.0" << std::endl;
            return 0;
        }
        else if (!std::strcmp(argv[cnt], "--output-prefix"))
        {
            OUTPUT_PREFIX = argv[cnt + 1];
            cnt += 2;
        }
        else
            throw std::invalid_argument("Unrecognized option: \"" + std::string(argv[cnt]) + "\".");
    }

    std::cout << "\nOutput directory set to: ";
    {
        if (RUN_TAG.empty())
            runtime_str(RUN_TAG);

        if (!std::filesystem::exists(RUN_TAG) && !std::filesystem::create_directory(RUN_TAG))
            throw failed_to_create_folder(RUN_TAG);
    }
    std::cout << "\"" << RUN_TAG << "\"" << std::endl;

    if (resume_mode)
    {
        size_t latest = 0;
        const std::regex data_file("([a-zA-Z]*)(\\d+)\\.dat");
        for (auto &it : std::filesystem::directory_iterator(RUN_TAG))
        {
            std::string s = it.path().filename().string();
            std::smatch sm;
            if (std::regex_match(s, sm, data_file))
            {
                std::string prefix = sm[1];
                if (prefix == OUTPUT_PREFIX)
                {
                    std::string idx = sm[2];
                    char *pEnd;
                    size_t idx10 = std::strtol(idx.c_str(), &pEnd, 10);
                    if (idx10 > latest)
                        latest = idx10;
                }
            }
        }
        auto output_dir = std::filesystem::path(RUN_TAG);
        std::string data_name = OUTPUT_PREFIX + std::to_string(latest) + ".dat";
        auto p_data = output_dir.append(data_name);
        DATA_PATH = p_data.string();
    }

    /// Report
    std::cout << "\ndt=" << dt << "s" << std::endl;
    std::cout << "\nMax iterations: " << MAX_ITER << std::endl;
    std::cout << "\nMax run time: " << MAX_TIME << "s" << std::endl;
    std::cout << "\nRecord solution every " << OUTPUT_GAP << " iteration" << std::endl;

    /// Init
    std::cout << "\nLoading mesh from \"" << MESH_PATH << "\" ... ";
    {
        std::ifstream in(MESH_PATH);
        if (in.fail())
            throw failed_to_open_file(MESH_PATH);
        tick_begin = clock();
        read_mesh(in);
        tick_end = clock();
        in.close();
    }
    std::cout << duration(tick_begin, tick_end) << "s" << std::endl;

    std::cout << "\nPreparing geometric quantities ... ";
    {
        tick_begin = clock();
        calculate_geometric_value();
        tick_end = clock();
    }
    std::cout << duration(tick_begin, tick_end) << "s" << std::endl;

    std::cout << "\nCalculating skewness factor on each face ... " << std::endl;
    check_skewness();

    std::cout << "\nSetting B.C. for each patch ... ";
    {
        set_bc_desc();
        set_bc_val();
    }
    std::cout << "Done!" << std::endl;

    std::cout << "\nPreparing Least-Square coefficients ... ";
    {
        tick_begin = clock();
        //calc_lsq_coefficient_matrix();
        tick_end = clock();
    }
    std::cout << duration(tick_begin, tick_end) << "s" << std::endl;

    std::cout << "\nPreparing Poisson equation coefficients ... ";
    {
        tick_begin = clock();
        /// TODO
        tick_end = clock();
    }
    std::cout << duration(tick_begin, tick_end) << "s" << std::endl;

    if (DATA_PATH.empty())
    {
        std::cout << "\nSetting I.C. ... ";
        {
            tick_begin = clock();
            zero_init();
            interpolate_nodal_value();
            tick_end = clock();
        }
        std::cout << duration(tick_begin, tick_end) << "s" << std::endl;

        std::cout << "\nWriting initial output ... ";
        {
            std::filesystem::path p_output(RUN_TAG);
            p_output.append(OUTPUT_PREFIX + "0.txt");
            std::ofstream dts(p_output);
            if (dts.fail())
                throw failed_to_open_file(p_output.filename());
            tick_begin = clock();
            write_data(dts, 0, 0.0);
            tick_end = clock();
            dts.close();
        }
        std::cout << duration(tick_begin, tick_end) << "s" << std::endl;
    }
    else
    {
        std::cout << "\nSetting I.C. from \"" + DATA_PATH + "\" ... ";
        std::ifstream dts(DATA_PATH);
        if (dts.fail())
            throw failed_to_open_file(DATA_PATH);
        tick_begin = clock();
        read_data(dts, iter, t);
        tick_end = clock();
        dts.close();
        std::cout << duration(tick_begin, tick_end) << "s" << std::endl;
    }

    /// Solve
    std::cout << "\nStarting calculation ... " << std::endl;
    while (iter <= MAX_ITER && t <= MAX_TIME)
    {
        ++iter;
        t += dt;

        /// Time-Stepping
        std::cout << "\nIter" << iter << ": " << "t=" << t << "s, dt=" << dt << "s" << std::endl;
        {
            tick_begin = clock();
            //ForwardEuler(dt);
            tick_end = clock();
        }
        std::cout << duration(tick_begin, tick_end) << "s CPU time" << std::endl;

        /// Check
        bool diverge_flag = false;
        diagnose(diverge_flag);
        if (diverge_flag)
        {
            /// TODO
        }

        /// Output
        if (!(iter % OUTPUT_GAP))
        {
            const std::string fn = OUTPUT_PREFIX + std::to_string(iter) + ".txt";
            std::filesystem::path p_output(RUN_TAG);
            p_output.append(fn);
            std::ofstream dts(p_output);
            if (dts.fail())
                throw failed_to_open_file(p_output.filename());
            write_data(dts, iter, t);
            dts.close();
        }
    }

    /// Finalize
    std::cout << "\nReleasing Memory ... " << std::endl;
    {
        for (auto e : node)
            delete e;
        for (auto e : face)
            delete e;
        for (auto e : cell)
            delete e;
        for (auto e : patch)
            delete e;
    }
    std::cout << "\nFinished!" << std::endl;

    return 0;
}
