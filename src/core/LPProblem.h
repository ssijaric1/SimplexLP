//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  LPProblem.h
//  Linear program in standard form:
//
//        minimize    c^T x
//        subject to  A x = b,   x >= 0,        A in R^{m x n} (sparse)
//
//  plus a converter from the user-friendly inequality form
//
//        min/max  c^T x   s.t.   a_i^T x {<=,>=,=} b_i ,  x >= 0
//
//  which appends slack/surplus columns. Columns created as slacks are
//  remembered (slackOfRow) so the simplex can start from the all-slack
//  basis ("crash basis") whenever it is feasible and skip phase I.
#pragma once

#include "LPCommon.h"

namespace lp
{

enum class Relation : unsigned char { LessEq = 0, GreaterEq, Equal };

inline const char* toString(Relation r)
{
    switch (r)
    {
        case Relation::LessEq:    return "<=";
        case Relation::GreaterEq: return ">=";
        case Relation::Equal:     return "=";
    }
    return "?";
}

// One row of the user-entered model (kept for visualization, MPS export
// and for building the standard form).
struct IneqConstraint
{
    std::vector<double> a;  // dense coefficients over the structural variables
    Relation rel = Relation::LessEq;
    double   b   = 0.0;
};

// User-level model with named structural variables.
struct UserModel
{
    bool maximize = false;
    std::vector<double>       c;       // objective over structural variables
    std::vector<IneqConstraint> rows;
    std::vector<std::string>  varNames;
    std::string               name = "USERLP";

    int nVars() const { return (int)c.size(); }
};

// Standard-form LP (what both solvers consume).
class StandardLP
{
public:
    SparseColMatrix A;                  // m x n
    std::vector<double> b;              // m
    std::vector<double> c;              // n
    std::vector<std::string> colNames;  // n (for MPS / logs)
    std::vector<std::string> rowNames;  // m
    std::string name = "LP";

    // Crash-basis support: slackOfRow[i] = column index of the slack of row i
    // with coefficient +1, or -1 when the row has no such slack.
    std::vector<int> slackOfRow;

    // Bookkeeping for reporting in the original (user) space.
    int nStructural = 0;     // first nStructural columns are the user variables
    double objSign   = 1.0;  // +1 if user minimized, -1 if user maximized
    double objShift  = 0.0;  // constant term (unused for now, kept for MPS RHS of obj)

    int rows() const { return A.nRows; }
    int colsTotal() const { return A.nCols; }

    bool valid() const
    {
        return A.nRows > 0 && A.nCols > 0 &&
               (int)b.size() == A.nRows && (int)c.size() == A.nCols;
    }

    // User-space objective value from a standard-form x.
    double userObjective(const double* x) const
    {
        double v = 0.0;
        for (int j = 0; j < (int)c.size(); ++j)
            v += c[j] * x[j];
        return objSign * v + objShift;
    }
};

// ---------------------------------------------------------------------------
// Conversion: user inequality model -> standard form (min, equalities, x>=0)
//  * maximize c^T x   becomes  minimize (-c)^T x   (objSign restores value)
//  * a^T x <= b       becomes  a^T x + s = b,  s >= 0   (slack, crash basis)
//  * a^T x >= b       becomes  a^T x - s = b,  s >= 0   (surplus)
//  * a^T x  = b       stays as it is
// ---------------------------------------------------------------------------
inline StandardLP toStandardForm(const UserModel& um)
{
    StandardLP lp;
    lp.name        = um.name;
    lp.nStructural = um.nVars();
    lp.objSign     = um.maximize ? -1.0 : 1.0;

    const int m = (int)um.rows.size();
    const int n0 = um.nVars();

    // Count extra columns (one slack/surplus per inequality row).
    int nExtra = 0;
    for (const auto& r : um.rows)
        if (r.rel != Relation::Equal)
            ++nExtra;

    const int n = n0 + nExtra;
    lp.A.resize(m, n);
    lp.b.assign(m, 0.0);
    lp.c.assign(n, 0.0);
    lp.slackOfRow.assign(m, -1);
    lp.colNames.resize(n);
    lp.rowNames.resize(m);

    for (int j = 0; j < n0; ++j)
    {
        lp.c[j] = lp.objSign * um.c[j]; // min form
        lp.colNames[j] = (j < (int)um.varNames.size() && !um.varNames[j].empty())
                       ? um.varNames[j]
                       : ("X" + std::to_string(j + 1));
    }

    int extra = n0;
    for (int i = 0; i < m; ++i)
    {
        const IneqConstraint& r = um.rows[i];
        lp.rowNames[i] = "R" + std::to_string(i + 1);
        for (int j = 0; j < n0 && j < (int)r.a.size(); ++j)
            lp.A.add(i, j, r.a[j]);
        lp.b[i] = r.b;

        if (r.rel == Relation::LessEq)
        {
            lp.A.add(i, extra, 1.0);
            lp.colNames[extra] = "S" + std::to_string(i + 1);
            lp.slackOfRow[i] = extra;       // +1 slack -> usable for crash basis
            ++extra;
        }
        else if (r.rel == Relation::GreaterEq)
        {
            lp.A.add(i, extra, -1.0);
            lp.colNames[extra] = "S" + std::to_string(i + 1);
            // surplus has coefficient -1: NOT a crash-basis slack
            ++extra;
        }
    }
    return lp;
}

} // namespace lp
