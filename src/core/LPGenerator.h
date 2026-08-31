//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  LPGenerator.h
//  Synthetic benchmark instances for the Simplex-vs-IPM comparison.
//
//  The generator builds standard-form problems  min c^T x, Ax = b, x >= 0
//  that are feasible AND bounded *by construction*:
//
//      1. sparse A (m x n, n > m) with a guaranteed structural rank
//         (one round-robin entry per row first, then random fill),
//      2. a strictly positive primal point  x* > 0   ->  b := A x*,
//      3. a dual pair (y*, s* > 0)               ->  c := A^T y* + s*.
//
//  x* is primal feasible and (y*, s*) is dual feasible, so by strong
//  duality an optimal solution exists. Both solvers therefore always have
//  a well-defined target, and any disagreement between them on the optimal
//  value indicates a bug - which is exactly what the test suite checks.
#pragma once

#include "LPProblem.h"
#include <random>

namespace lp
{

struct GeneratorParams
{
    int      m         = 50;     // rows (equalities)
    int      n         = 100;    // columns, must be > m
    int      nnzPerCol = 5;      // target nonzeros per column (>= 1)
    unsigned seed      = 42;
};

inline StandardLP generateInstance(const GeneratorParams& gp)
{
    std::mt19937 rng(gp.seed);
    std::uniform_real_distribution<double> uVal(-1.0, 1.0);
    std::uniform_real_distribution<double> uPos(0.1, 10.0);
    std::uniform_int_distribution<int>     uRow(0, gp.m - 1);

    StandardLP P;
    P.name        = "GEN_M" + std::to_string(gp.m) + "_N" + std::to_string(gp.n);
    P.nStructural = gp.n;
    P.A.resize(gp.m, gp.n);
    P.b.assign(gp.m, 0.0);
    P.c.assign(gp.n, 0.0);
    P.slackOfRow.assign(gp.m, -1); // equality form: no crash slacks
    P.colNames.resize(gp.n);
    P.rowNames.resize(gp.m);

    for (int j = 0; j < gp.n; ++j)
        P.colNames[j] = "X" + std::to_string(j + 1);
    for (int i = 0; i < gp.m; ++i)
        P.rowNames[i] = "R" + std::to_string(i + 1);

    // ---- sparse A ----------------------------------------------------------
    // Column j gets one deterministic entry on row (j mod m) - this covers
    // every row at least ceil(n/m) times and keeps A full row rank with very
    // high probability - plus (nnzPerCol-1) random extra rows.
    std::vector<char> used(gp.m, 0);
    for (int j = 0; j < gp.n; ++j)
    {
        std::fill(used.begin(), used.end(), 0);

        const int anchor = j % gp.m;
        double v = uVal(rng);
        if (std::fabs(v) < 0.1) v = (v < 0 ? -0.1 : 0.1) ;
        P.A.add(anchor, j, v);
        used[anchor] = 1;

        for (int t = 1; t < gp.nnzPerCol; ++t)
        {
            int i = uRow(rng);
            if (used[i])
                continue; // slightly fewer nnz is fine; no duplicates
            used[i] = 1;
            double w = uVal(rng);
            if (std::fabs(w) < 0.05) w = (w < 0 ? -0.05 : 0.05);
            P.A.add(i, j, w);
        }
    }

    // ---- primal certificate: x* > 0,  b = A x* -----------------------------
    std::vector<double> xStar(gp.n);
    for (int j = 0; j < gp.n; ++j)
        xStar[j] = uPos(rng);
    P.A.mulVec(xStar.data(), P.b.data());

    // ---- dual certificate: c = A^T y* + s*,  s* > 0 -------------------------
    std::vector<double> yStar(gp.m), sStar(gp.n);
    for (int i = 0; i < gp.m; ++i)
        yStar[i] = uVal(rng);
    for (int j = 0; j < gp.n; ++j)
        sStar[j] = std::uniform_real_distribution<double>(0.05, 2.0)(rng);

    P.A.mulTransVec(yStar.data(), P.c.data());
    for (int j = 0; j < gp.n; ++j)
        P.c[j] += sStar[j];

    return P;
}

} // namespace lp
