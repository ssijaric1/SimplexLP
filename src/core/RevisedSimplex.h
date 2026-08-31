//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  RevisedSimplex.h
//  Two-phase revised simplex method for   min c^T x,  A x = b,  x >= 0.
//
//  The "revised" formulation never forms the full tableau. Per iteration it
//  keeps only the basis matrix B = A(:, basis) and computes
//
//      x_B = B^{-1} b                  (basic solution,        FTRAN)
//      y   = B^{-T} c_B                (simplex multipliers,   BTRAN)
//      d_j = c_j - y^T a_j             (reduced costs, sparse dot products)
//      w   = B^{-1} a_q                (entering column,       FTRAN)
//
//  All B-solves go through natID's sparse LU (see BasisFactor.h), exactly as
//  stated in the project proposal ("the basis will be updated by solving the
//  linear system B d = a_j ... using natID's sparse LU solver").
//
//  Pivoting rules:
//    * entering: Dantzig (most negative reduced cost); after `blandAfter`
//      consecutive degenerate steps the code switches to Bland's rule,
//      which provably prevents cycling.
//    * leaving: minimum ratio test  theta = min x_B(i)/w_i over w_i > tol,
//      ties broken by the largest |w_i| (numerical stability).
//
//  Phase I uses as few artificial variables as possible: rows whose slack
//  has coefficient +1 (after sign-normalizing b >= 0) start with the slack
//  in the basis; only the remaining rows receive an artificial with cost 1.
//  If every row is covered by a slack, phase I is skipped entirely
//  ("crash basis" - the common case for <= models with b >= 0).
#pragma once

#include "LPProblem.h"
#include "BasisFactor.h"

namespace lp
{

class RevisedSimplex
{
public:
    LPResult solve(const StandardLP& P, const SimplexOptions& opt = {})
    {
        LPResult res;
        if (!P.valid())
        {
            res.status  = Status::NumericalError;
            res.message = "invalid problem dimensions";
            return res;
        }

        StopWatch total;

        _opt = opt;
        _m = P.rows();
        _n = P.colsTotal();

        // ---- 1. row-sign normalization: bw >= 0 --------------------------
        // Working system (S A) x = S b with S = diag(rowSign); duals map back
        // through y_orig = S y_w at the end.
        _rowSign.assign(_m, 1.0);
        _bw.assign(_m, 0.0);
        for (int i = 0; i < _m; ++i)
        {
            _rowSign[i] = (P.b[i] < 0.0) ? -1.0 : 1.0;
            _bw[i] = _rowSign[i] * P.b[i];
        }

        _Aw.resize(_m, _n); // filled below; artificial columns appended after
        for (int j = 0; j < _n; ++j)
            for (const SparseEntry& e : P.A.column(j))
                _Aw.add(e.row, j, _rowSign[e.row] * e.val);

        // ---- 2. initial basis: slacks where possible, artificials else ---
        _basis.assign(_m, -1);
        std::vector<int> artRow; // rows that need an artificial column
        for (int i = 0; i < _m; ++i)
        {
            const int sj = (i < (int)P.slackOfRow.size()) ? P.slackOfRow[i] : -1;
            if (sj >= 0 && _rowSign[i] > 0.0)
                _basis[i] = sj;      // slack column is +e_i, x_slack = bw_i >= 0
            else
                artRow.push_back(i);
        }

        const int nArt = (int)artRow.size();
        _nTotal = _n + nArt;
        for (int t = 0; t < nArt; ++t)
        {
            _Aw.cols.push_back({});
            _Aw.nCols++;
            _Aw.add(artRow[t], _n + t, 1.0);
            _basis[artRow[t]] = _n + t;
        }

        _inBasis.assign(_nTotal, 0);
        for (int i = 0; i < _m; ++i)
            _inBasis[_basis[i]] = 1;

        _xB.assign(_m, 0.0);
        _y.assign(_m, 0.0);
        _w.assign(_m, 0.0);
        _colDense.assign(_m, 0.0);
        _cB.assign(_m, 0.0);

        // ---- 3. phase I ---------------------------------------------------
        res.stats.phase1Iters = 0;
        if (nArt > 0)
        {
            std::vector<double> c1(_nTotal, 0.0);
            for (int t = 0; t < nArt; ++t)
                c1[_n + t] = 1.0;

            Status st = runSimplexLoop(c1, /*phase=*/1, /*banArtificials=*/false, res);
            if (st != Status::Optimal)
            {
                // Phase I objective is bounded below by 0 - "unbounded" here
                // can only mean numerical trouble.
                res.status = (st == Status::Unbounded) ? Status::NumericalError : st;
                res.stats.msTotal = total.ms();
                finalizeSolution(P, res);
                return res;
            }

            const double phase1Obj = basicObjective(c1);
            if (phase1Obj > _opt.tolFeas)
            {
                res.status  = Status::Infeasible;
                res.message = "phase I optimum = " + std::to_string(phase1Obj);
                res.stats.msTotal = total.ms();
                finalizeSolution(P, res);
                return res;
            }

            if (!driveOutArtificials(res))
            {
                res.status = Status::NumericalError;
                res.stats.msTotal = total.ms();
                finalizeSolution(P, res);
                return res;
            }
        }

        // ---- 4. phase II --------------------------------------------------
        std::vector<double> c2(_nTotal, 0.0);
        for (int j = 0; j < _n; ++j)
            c2[j] = P.c[j];

        Status st = runSimplexLoop(c2, /*phase=*/2, /*banArtificials=*/true, res);
        res.status = st;

        res.stats.msTotal = total.ms();
        finalizeSolution(P, res);
        return res;
    }

private:
    // ---- workspace --------------------------------------------------------
    SimplexOptions _opt;
    int _m = 0, _n = 0, _nTotal = 0;

    SparseColMatrix     _Aw;       // sign-normalized A, artificials appended
    std::vector<double> _bw;       // sign-normalized b (>= 0)
    std::vector<double> _rowSign;
    std::vector<int>    _basis;    // basis[i] = column in row position i
    std::vector<char>   _inBasis;

    std::vector<double> _xB, _y, _w, _colDense, _cB;
    BasisFactor _factor;
    int _iterGlobal = 0;

    // ------------------------------------------------------------------------
    double basicObjective(const std::vector<double>& c) const
    {
        double obj = 0.0;
        for (int i = 0; i < _m; ++i)
            obj += c[_basis[i]] * _xB[i];
        return obj;
    }

    void gatherColumnDense(int j)
    {
        std::fill(_colDense.begin(), _colDense.end(), 0.0);
        for (const SparseEntry& e : _Aw.column(j))
            _colDense[e.row] = e.val;
    }

    bool refactorizeAndSolve(LPResult& res, const std::vector<double>& c)
    {
        StopWatch sw;
        if (!_factor.factorize(_Aw, _basis))
        {
            res.message = _factor.lastError();
            return false;
        }
        res.stats.msFactorize += sw.ms();
        res.stats.factorizations++;

        for (int i = 0; i < _m; ++i)
            _cB[i] = c[_basis[i]];

        sw.restart();
        const bool ok1 = _factor.ftran(_bw.data(), _xB.data());
        const bool ok2 = _factor.btran(_cB.data(), _y.data());
        res.stats.msSolves += sw.ms();

        if (!ok1 || !ok2)
        {
            res.message = "triangular solve failed";
            return false;
        }
        return true;
    }

    void recordSnapshot(LPResult& res, int phase, int entering, int leaving,
                        const std::vector<double>& c)
    {
        if (!_opt.recordPath)
            return;
        IterSnapshot s;
        s.iter      = _iterGlobal;
        s.phase     = phase;
        s.entering  = entering;
        s.leaving   = leaving;
        s.objective = basicObjective(c);
        s.x.assign(_n, 0.0);
        for (int i = 0; i < _m; ++i)
            if (_basis[i] < _n)
                s.x[_basis[i]] = _xB[i];
        res.path.push_back(std::move(s));
    }

    // The main pivoting loop, shared by both phases.
    Status runSimplexLoop(const std::vector<double>& c, int phase,
                          bool banArtificials, LPResult& res)
    {
        int degenerateRun = 0;

        if (!refactorizeAndSolve(res, c))
            return Status::NumericalError;

        // snapshot of the starting vertex
        recordSnapshot(res, phase, -1, -1, c);

        while (true)
        {
            if (_iterGlobal >= _opt.maxIterations)
                return Status::IterationLimit;

            // ---- pricing: choose entering column -------------------------
            const bool useBland = (degenerateRun >= _opt.blandAfter);
            int    q     = -1;
            double bestD = -_opt.tolOpt;

            for (int j = 0; j < _nTotal; ++j)
            {
                if (_inBasis[j])
                    continue;
                if (banArtificials && j >= _n)
                    continue;

                const double dj = c[j] - _Aw.dotColumn(j, _y.data());

                if (useBland)
                {
                    if (dj < -_opt.tolOpt) { q = j; break; }       // Bland: first index
                }
                else if (dj < bestD)
                {
                    bestD = dj; q = j;                              // Dantzig: most negative
                }
            }

            if (q < 0)
                return Status::Optimal; // no improving column

            // ---- FTRAN: w = B^{-1} a_q ------------------------------------
            gatherColumnDense(q);
            {
                StopWatch sw;
                if (!_factor.ftran(_colDense.data(), _w.data()))
                {
                    res.message = "FTRAN failed";
                    return Status::NumericalError;
                }
                res.stats.msSolves += sw.ms();
            }

            // ---- ratio test ------------------------------------------------
            int    r     = -1;
            double theta = 0.0;
            for (int i = 0; i < _m; ++i)
            {
                if (_w[i] > _opt.tolPivot)
                {
                    const double ratio = _xB[i] / _w[i];
                    if (r < 0 || ratio < theta - 1e-12 ||
                        (ratio < theta + 1e-12 && std::fabs(_w[i]) > std::fabs(_w[r])))
                    {
                        theta = ratio;
                        r = i;
                    }
                }
            }

            if (r < 0)
                return Status::Unbounded; // entering direction never blocked

            degenerateRun = (theta <= _opt.tolPivot) ? degenerateRun + 1 : 0;

            // ---- pivot ------------------------------------------------------
            const int leaving = _basis[r];
            _inBasis[leaving] = 0;
            _inBasis[q]       = 1;
            _basis[r]         = q;

            ++_iterGlobal;
            res.stats.iterations++;
            if (phase == 1)
                res.stats.phase1Iters++;

            // ---- refactorize new basis & resolve ---------------------------
            // One sparse LU per pivot (plus one for B^T). Product-form /
            // Forrest-Tomlin updates are noted as future work in the report.
            if (!refactorizeAndSolve(res, c))
                return Status::NumericalError;

            recordSnapshot(res, phase, q, leaving, c);

            if (_opt.verbose)
            {
                std::fprintf(stderr, "[simplex] it=%d ph=%d in=%d out=%d theta=%.3g obj=%.8g\n",
                             _iterGlobal, phase, q, leaving, theta, basicObjective(c));
            }
        }
    }

    // After a successful phase I, pivot remaining artificials out of the
    // basis whenever possible. A row whose artificial cannot be replaced is
    // linearly dependent (redundant); its artificial stays basic at value 0
    // and is simply never priced in phase II.
    bool driveOutArtificials(LPResult& res)
    {
        for (int r = 0; r < _m; ++r)
        {
            if (_basis[r] < _n)
                continue; // not an artificial

            // z = B^{-T} e_r  =>  (B^{-1} a_j)_r = z^T a_j  for every column j
            std::fill(_colDense.begin(), _colDense.end(), 0.0);
            _colDense[r] = 1.0;
            std::vector<double> z(_m, 0.0);
            if (!_factor.btran(_colDense.data(), z.data()))
            {
                res.message = "BTRAN failed while driving out artificials";
                return false;
            }

            int replacement = -1;
            for (int j = 0; j < _n; ++j)
            {
                if (_inBasis[j])
                    continue;
                if (std::fabs(_Aw.dotColumn(j, z.data())) > _opt.tolPivot)
                {
                    replacement = j;
                    break;
                }
            }

            if (replacement >= 0)
            {
                const int art = _basis[r];
                _inBasis[art]         = 0;
                _inBasis[replacement] = 1;
                _basis[r]             = replacement; // degenerate pivot, theta = 0

                if (!refactorizeAndSolve(res, std::vector<double>(_nTotal, 0.0)))
                    return false;
            }
            // else: redundant row, artificial remains basic at 0 - safe.
        }
        return true;
    }

    void finalizeSolution(const StandardLP& P, LPResult& res)
    {
        res.x.reserve((size_t)_n);
        res.y.reserve((size_t)_m);
        for (int j = 0; j < _n; ++j)
            res.x[j] = 0.0;
        for (int i = 0; i < _m; ++i)
        {
            if (_basis[i] < _n)
                res.x[_basis[i]] = _xB[i];
            res.y[i] = _rowSign[i] * _y[i]; // map duals back to original rows
        }

        double obj = 0.0;
        for (int j = 0; j < _n; ++j)
            obj += P.c[j] * res.x[j];
        res.objective = obj; // standard-form (min) objective
    }
};

} // namespace lp
