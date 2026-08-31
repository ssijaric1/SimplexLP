//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  Verify.h
//  Independent optimality check used by the test suite and shown in the
//  application log. A claimed solution (x, y) of  min c^T x, Ax=b, x>=0
//  is optimal iff the KKT conditions hold:
//
//      primal feasibility   A x = b,  x >= 0
//      dual   feasibility   s := c - A^T y >= 0
//      complementarity      x_j * s_j = 0   for all j
//
//  equivalently (strong duality)  c^T x = b^T y. All measures are reported
//  relative so they are meaningful across problem scales. This check is
//  solver-agnostic: it validates simplex and IPM results the same way and
//  is the ground truth when the two methods are compared.
#pragma once

#include "LPProblem.h"

namespace lp
{

struct KKTReport
{
    double primalResidual = 0.0; // ||Ax-b||_inf / (1+||b||_inf)
    double primalNegativity = 0.0; // max(0, -min_j x_j)
    double dualNegativity   = 0.0; // max(0, -min_j s_j)
    double complementarity  = 0.0; // max_j |x_j s_j| / (1+|c^T x|)
    double dualityGap       = 0.0; // |c^T x - b^T y| / (1+|c^T x|)
    double primalObjective  = 0.0;
    double dualObjective    = 0.0;

    bool ok(double tol = 1e-6) const
    {
        return primalResidual   < tol &&
               primalNegativity < tol &&
               dualNegativity   < tol &&
               dualityGap       < tol;
    }
};

inline KKTReport verifyKKT(const StandardLP& P,
                           const cnt::SafeFullVector<double>& x,
                           const cnt::SafeFullVector<double>& y)
{
    KKTReport rep;
    const int m = P.rows();
    const int n = P.colsTotal();

    std::vector<double> Ax(m, 0.0);
    std::vector<double> xs((size_t)n);
    for (int j = 0; j < n; ++j)
        xs[j] = x[(size_t)j];
    P.A.mulVec(xs.data(), Ax.data());

    double normB = 0.0, rPrim = 0.0;
    for (int i = 0; i < m; ++i)
    {
        normB = std::max(normB, std::fabs(P.b[i]));
        rPrim = std::max(rPrim, std::fabs(Ax[i] - P.b[i]));
    }
    rep.primalResidual = rPrim / (1.0 + normB);

    double minX = 0.0;
    for (int j = 0; j < n; ++j)
        minX = std::min(minX, x[(size_t)j]);
    rep.primalNegativity = std::max(0.0, -minX);

    // s = c - A^T y
    std::vector<double> yv((size_t)m);
    for (int i = 0; i < m; ++i)
        yv[i] = y[(size_t)i];
    std::vector<double> s(n, 0.0);
    P.A.mulTransVec(yv.data(), s.data());

    double cx = 0.0, by = 0.0, minS = 0.0, comp = 0.0;
    for (int j = 0; j < n; ++j)
    {
        s[j] = P.c[j] - s[j];
        cx  += P.c[j] * x[(size_t)j];
        minS = std::min(minS, s[j]);
        comp = std::max(comp, std::fabs(s[j] * x[(size_t)j]));
    }
    for (int i = 0; i < m; ++i)
        by += P.b[i] * y[(size_t)i];

    rep.primalObjective = cx;
    rep.dualObjective   = by;
    rep.dualNegativity  = std::max(0.0, -minS);
    rep.complementarity = comp / (1.0 + std::fabs(cx));
    rep.dualityGap      = std::fabs(cx - by) / (1.0 + std::fabs(cx));
    return rep;
}

} // namespace lp
