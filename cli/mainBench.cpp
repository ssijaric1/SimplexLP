//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  mainBench.cpp - console benchmark (target: lpBench).
//
//  Generates a sweep of synthetic standard-form LP instances of increasing
//  size, solves each with the revised simplex and with the interior-point
//  method, cross-checks optimality (KKT residuals + objective gap), and
//  writes results.csv plus the .mps files needed for the external COIN-OR
//  CLP comparison (scripts/run_clp.sh consumes exactly these files).
//
//  Usage:
//      lpBench [mStart mEnd steps nOverM nnzPerCol seed outDir]
//  Defaults:
//      lpBench 20 200 6 2.5 4 42 ./SimplexLP_bench
#include <mu/Application.h>

#include "../src/core/LPGenerator.h"
#include "../src/core/RevisedSimplex.h"
#include "../src/core/InteriorPoint.h"
#include "../src/core/MPS.h"
#include "../src/core/Verify.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <cmath>

int main(int argc, const char* argv[])
{
    mu::Application app(argc, argv);

    int      mStart = 20, mEnd = 200, steps = 6;
    double   nOverM = 2.5;
    int      nnzPerCol = 4;
    unsigned seed = 42;
    std::string outDir = "./SimplexLP_bench";

    if (argc > 1) mStart    = std::atoi(argv[1]);
    if (argc > 2) mEnd      = std::atoi(argv[2]);
    if (argc > 3) steps     = std::atoi(argv[3]);
    if (argc > 4) nOverM    = std::atof(argv[4]);
    if (argc > 5) nnzPerCol = std::atoi(argv[5]);
    if (argc > 6) seed      = (unsigned)std::atoi(argv[6]);
    if (argc > 7) outDir    = argv[7];

    if (mStart < 2) mStart = 2;
    if (mEnd < mStart) mEnd = mStart;
    if (steps < 1) steps = 1;
    if (nOverM < 1.1) nOverM = 1.1;
    if (nnzPerCol < 1) nnzPerCol = 1;

    std::error_code ec;
    std::filesystem::create_directories(outDir, ec);
    if (ec)
    {
        std::cerr << "cannot create output dir " << outDir << ": "
                  << ec.message() << "\n";
        return 1;
    }

    const std::string csvPath = outDir + "/results.csv";
    std::ofstream csv(csvPath);
    csv << "m,n,spx_status,spx_ms,spx_iters,spx_factorizations,spx_obj,"
           "ipm_status,ipm_ms,ipm_iters,ipm_obj,obj_gap,spx_kkt_ok,ipm_kkt_ok\n";

    std::cout << "SimplexLP benchmark  m=" << mStart << ".." << mEnd
              << "  steps=" << steps << "  n/m=" << nOverM
              << "  nnz/col=" << nnzPerCol << "  seed=" << seed << "\n"
              << "output: " << outDir << "\n\n";

    char buf[320];
    for (int s = 0; s < steps; ++s)
    {
        const double t = (steps == 1) ? 0.0 : (double)s / (double)(steps - 1);
        const int m = (int)std::lround(mStart + t * (mEnd - mStart));
        const int n = std::max(m + 1, (int)std::lround(m * nOverM));

        lp::GeneratorParams gp;
        gp.m = m;
        gp.n = n;
        gp.nnzPerCol = nnzPerCol;
        gp.seed = seed + (unsigned)s * 1000u;
        lp::StandardLP prob = lp::generateInstance(gp);

        std::snprintf(buf, sizeof(buf), "%s/bench_m%04d_n%05d.mps",
                      outDir.c_str(), m, n);
        std::string err;
        if (!lp::writeMPS(prob, buf, &err))
            std::cerr << "MPS export failed: " << err << "\n";

        // ---- revised simplex ----
        lp::RevisedSimplex spx;
        lp::SimplexOptions so;
        lp::LPResult rs = spx.solve(prob, so);
        lp::KKTReport ks = lp::verifyKKT(prob, rs.x, rs.y);

        // ---- interior point ----
        lp::InteriorPoint ipm;
        lp::IPMOptions io;
        lp::LPResult ri = ipm.solve(prob, io);
        lp::KKTReport ki = lp::verifyKKT(prob, ri.x, ri.y);

        const double gap = std::fabs(rs.objective - ri.objective) /
                           (1.0 + std::fabs(rs.objective));

        std::snprintf(buf, sizeof(buf),
            "m=%-4d n=%-5d | simplex %-9s %8.2f ms %5d it (%d LU) | "
            "ipm %-9s %8.2f ms %3d it | gap %.2e",
            m, n,
            lp::toString(rs.status), rs.stats.msTotal, rs.stats.iterations,
            rs.stats.factorizations,
            lp::toString(ri.status), ri.stats.msTotal, ri.stats.iterations,
            gap);
        std::cout << buf << "\n";

        csv << m << ',' << n << ','
            << lp::toString(rs.status) << ',' << rs.stats.msTotal << ','
            << rs.stats.iterations << ',' << rs.stats.factorizations << ','
            << rs.objective << ','
            << lp::toString(ri.status) << ',' << ri.stats.msTotal << ','
            << ri.stats.iterations << ',' << ri.objective << ','
            << gap << ','
            << (ks.ok(1e-6) ? 1 : 0) << ',' << (ki.ok(1e-6) ? 1 : 0) << '\n';
    }

    std::cout << "\nwrote " << csvPath << "\n"
              << "compare with CLP:  scripts/run_clp.sh " << outDir << "\n";
    return 0;
}
