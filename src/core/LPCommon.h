//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  LPCommon.h
//  Shared basic types for the LP core: solver status, options, statistics,
//  iteration records (for visualization) and a small column-major sparse
//  matrix container that feeds natID's sparse::ISolver with triples.
//
//  The heavy numerical work (LU / LLT factorization and triangular solves)
//  is delegated to natID (sparse::createDblSolver). This file holds only
//  problem data and bookkeeping.
#pragma once

#include <cnt/SafeFullVector.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <chrono>

namespace lp
{

// ---------------------------------------------------------------- status --
enum class Status : unsigned char
{
    NotSolved = 0,
    Optimal,
    Infeasible,
    Unbounded,
    IterationLimit,
    NumericalError
};

inline const char* toString(Status s)
{
    switch (s)
    {
        case Status::NotSolved:      return "NotSolved";
        case Status::Optimal:        return "Optimal";
        case Status::Infeasible:     return "Infeasible";
        case Status::Unbounded:      return "Unbounded";
        case Status::IterationLimit: return "IterationLimit";
        case Status::NumericalError: return "NumericalError";
    }
    return "?";
}

// ------------------------------------------------------------- stopwatch --
// Thin wrapper around std::chrono (the same clock mu::Timer uses on POSIX).
class StopWatch
{
    std::chrono::high_resolution_clock::time_point _t0;
public:
    StopWatch() { restart(); }
    void restart() { _t0 = std::chrono::high_resolution_clock::now(); }
    double ms() const
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(t1 - _t0).count();
    }
};

// ------------------------------------------------- sparse column storage --
// Column-major triplet storage: exactly what sparse::ISolver::addTriple
// wants, and what revised simplex needs (fast access to columns a_j).
struct SparseEntry
{
    int    row;
    double val;
};

class SparseColMatrix
{
public:
    int nRows = 0;
    int nCols = 0;
    std::vector<std::vector<SparseEntry>> cols; // cols[j] -> entries of column j

    SparseColMatrix() = default;
    SparseColMatrix(int m, int n) { resize(m, n); }

    void resize(int m, int n)
    {
        nRows = m;
        nCols = n;
        cols.assign((size_t)n, {});
    }

    void clear() { resize(0, 0); }

    inline void add(int i, int j, double v)
    {
        if (v != 0.0)
            cols[(size_t)j].push_back({i, v});
    }

    inline const std::vector<SparseEntry>& column(int j) const { return cols[(size_t)j]; }

    int nnz() const
    {
        size_t s = 0;
        for (const auto& c : cols)
            s += c.size();
        return (int)s;
    }

    // y += alpha * A(:,j)
    inline void axpyColumn(int j, double alpha, double* y) const
    {
        for (const auto& e : cols[(size_t)j])
            y[e.row] += alpha * e.val;
    }

    // dot( A(:,j), y )
    inline double dotColumn(int j, const double* y) const
    {
        double s = 0.0;
        for (const auto& e : cols[(size_t)j])
            s += e.val * y[e.row];
        return s;
    }

    // dense y = A * x  (y has nRows entries, x has nCols entries)
    void mulVec(const double* x, double* y) const
    {
        for (int i = 0; i < nRows; ++i)
            y[i] = 0.0;
        for (int j = 0; j < nCols; ++j)
        {
            const double xj = x[j];
            if (xj != 0.0)
                axpyColumn(j, xj, y);
        }
    }

    // dense y = A^T * x  (y has nCols entries, x has nRows entries)
    void mulTransVec(const double* x, double* y) const
    {
        for (int j = 0; j < nCols; ++j)
            y[j] = dotColumn(j, x);
    }
};

// --------------------------------------------------------------- options --
struct SimplexOptions
{
    double tolOpt        = 1e-9;   // reduced-cost optimality tolerance
    double tolPivot      = 1e-9;   // smallest acceptable pivot |w_r|
    double tolFeas       = 1e-7;   // phase-I objective threshold for feasibility
    int    maxIterations = 50000;
    int    blandAfter    = 60;     // consecutive degenerate pivots before Bland's rule
    bool   recordPath    = false;  // store x at every iteration (small demos only)
    bool   verbose       = false;
};

struct IPMOptions
{
    double tolMu         = 1e-9;   // complementarity target (relative)
    double tolResidual   = 1e-8;   // relative primal/dual residual target
    int    maxIterations = 200;
    double stepFraction  = 0.99;   // fraction-to-boundary
    double regularize    = 1e-10;  // diagonal regularization added to A D^2 A^T
    bool   useCholesky   = true;   // LLT on normal equations; falls back to LU
    bool   recordPath    = false;
    bool   verbose       = false;
};

// ------------------------------------------------------------ statistics --
struct SolveStats
{
    int    iterations     = 0;  // total pivot / IPM iterations
    int    phase1Iters    = 0;  // simplex only
    int    factorizations = 0;
    double msTotal        = 0.0;
    double msFactorize    = 0.0;
    double msSolves       = 0.0; // triangular solves (FTRAN/BTRAN or IPM solves)
};

// One snapshot per iteration; used by the 2D/3D pivot visualization.
struct IterSnapshot
{
    int    iter     = 0;
    int    phase    = 2;    // 1 / 2 = simplex phases, 0 = interior point
    int    entering = -1;   // column index entering the basis (-1 for IPM)
    int    leaving  = -1;   // column index leaving the basis  (-1 for IPM)
    double objective = 0.0;
    std::vector<double> x;  // full primal vector (structural + slacks)
};

// ----------------------------------------------------------------- result --
struct LPResult
{
    Status status = Status::NotSolved;
    double objective = 0.0;

    cnt::SafeFullVector<double> x;  // primal solution (size = problem nCols)
    cnt::SafeFullVector<double> y;  // dual solution / simplex multipliers (size = nRows)

    SolveStats stats;
    std::vector<IterSnapshot> path; // filled only when options.recordPath
    std::string message;            // human-readable detail (errors etc.)
};

} // namespace lp
