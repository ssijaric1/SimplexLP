//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  InteriorPoint.h
//  Primal-dual interior-point (barrier) method with Mehrotra's
//  predictor-corrector for   min c^T x,  A x = b,  x >= 0.
//
//  Per iteration the KKT Newton system is reduced to the normal equations
//
//      M dy = rhs,     M = A D A^T,   D = diag(x_j / s_j)  (SPD)
//
//  M is factorized once per iteration with natID's sparse solver in
//  symmetric positive-definite mode (Symmetry::SymmetricPosDef - only the
//  lower triangle is inserted, exactly like the SDK's own symmetric-matrix
//  tests - with SolverType::LU as the factorization engine). Predictor and
//  corrector right-hand sides reuse the same factorization through
//  solveExt(). If the SPD factorization fails numerically, the solver
//  transparently falls back to plain LU on the (mirrored) symmetric matrix.
//
//  This is the in-framework "interior point" counterpart for the
//  Simplex-vs-IPM comparison required by the project statement. The method
//  assumes the instance is feasible and bounded (guaranteed by construction
//  for the synthetic benchmark instances; see LPGenerator.h). Rigorous
//  infeasibility certificates would require a homogeneous self-dual
//  embedding, which is noted as future work in the report.
#pragma once

#include "LPProblem.h"
#include <sparse/ISolver.h>
#include <unordered_map>
#include <algorithm>

namespace lp
{

// Factorization helper for M(D) = A D A^T + reg*I.
// The sparsity pattern of M is computed once (symbolic step) and only the
// numerical values are re-accumulated every IPM iteration.
class NormalEqFactor
{
    const SparseColMatrix* _A = nullptr;
    int _m = 0;
    double _reg = 0.0;

    // pattern of the lower triangle (i >= k)
    std::vector<std::pair<int,int>> _pattern;
    std::unordered_map<long long, int> _index; // (i,k) -> position in _pattern
    std::vector<double> _values;

    sparse::DblSolverReleaser _solver;
    bool _useCholesky = true;
    std::string _lastError;

    static long long key(int i, int k) { return (long long)i * 1000000007LL + k; }

public:
    NormalEqFactor() : _solver(nullptr) {}

    const std::string& lastError() const { return _lastError; }
    bool usingCholesky() const { return _useCholesky; }

    void prepare(const SparseColMatrix& A, double regularize, bool useCholesky)
    {
        _A = &A;
        _m = A.nRows;
        _reg = regularize;
        _useCholesky = useCholesky;

        _pattern.clear();
        _index.clear();

        // diagonal first - guarantees every (i,i) exists for regularization
        for (int i = 0; i < _m; ++i)
        {
            _index.emplace(key(i, i), (int)_pattern.size());
            _pattern.emplace_back(i, i);
        }

        // symbolic A D A^T: every column contributes a small clique
        for (int j = 0; j < A.nCols; ++j)
        {
            const auto& col = A.column(j);
            const int   nc  = (int)col.size();
            for (int p = 0; p < nc; ++p)
            {
                for (int q = 0; q <= p; ++q)
                {
                    int i = col[p].row, k = col[q].row;
                    if (i < k) std::swap(i, k);
                    const long long kk = key(i, k);
                    if (_index.find(kk) == _index.end())
                    {
                        _index.emplace(kk, (int)_pattern.size());
                        _pattern.emplace_back(i, k);
                    }
                }
            }
        }
        _values.assign(_pattern.size(), 0.0);
    }

    // Accumulate values of A D A^T for D = diag(d) and factorize.
    bool factorize(const std::vector<double>& d)
    {
        _lastError.clear();
        std::fill(_values.begin(), _values.end(), 0.0);

        const SparseColMatrix& A = *_A;
        for (int j = 0; j < A.nCols; ++j)
        {
            const double dj = d[j];
            if (dj == 0.0)
                continue;
            const auto& col = A.column(j);
            const int   nc  = (int)col.size();
            for (int p = 0; p < nc; ++p)
            {
                const double vp = dj * col[p].val;
                for (int q = 0; q <= p; ++q)
                {
                    int i = col[p].row, k = col[q].row;
                    double contrib = vp * col[q].val;
                    if (i < k) std::swap(i, k);
                    _values[_index[key(i, k)]] += contrib;
                }
            }
        }
        // diagonal regularization (positions 0.._m-1 by construction)
        for (int i = 0; i < _m; ++i)
            _values[i] += _reg;

        if (_useCholesky)
        {
            if (load(true) && _solver->factorize())
                return true;
            // graceful fallback: LU on the mirrored matrix
            _useCholesky = false;
        }
        if (!load(false))
            return false;
        if (!_solver->factorize())
        {
            const char* e = _solver->getLastError();
            _lastError = std::string("normal-equation factorization failed: ") + (e ? e : "?");
            return false;
        }
        return true;
    }

    bool solve(const double* rhs, double* x)
    {
        return _solver->solveExt(rhs, x);
    }

private:
    bool load(bool cholesky)
    {
        const int nnz = (int)_pattern.size();
        const int reserveNZ = cholesky ? (nnz + _m + 8) : (2 * nnz + _m + 8);
        const double zero = 0.0;

        if (cholesky)
            // SymmetricPosDef tells the solver to exploit SPD structure
            // (insert only the lower triangle); the factorization engine is
            // LU. This build asserts on SolverType::LLT, while
            // SymmetricPosDef + LU is the pattern the SDK's own SPD tests use.
            _solver = sparse::createDblSolver(_m, reserveNZ,
                                              sparse::Symmetry::SymmetricPosDef,
                                              sparse::SolverType::LU,
                                              sparse::Pivoting::No);
        else
            _solver = sparse::createDblSolver(_m, reserveNZ,
                                              sparse::Symmetry::NonSymmetric,
                                              sparse::SolverType::LU,
                                              sparse::Pivoting::MarkowitzSinglePass);
        if (!_solver)
        {
            _lastError = "createDblSolver returned null";
            return false;
        }

        _solver->populateDiagonals(zero);
        for (size_t t = 0; t < _pattern.size(); ++t)
        {
            const int i = _pattern[t].first;
            const int k = _pattern[t].second;
            const double v = _values[t];
            _solver->addTriple(i, k, v);             // lower triangle
            if (!cholesky && i != k)
                _solver->addTriple(k, i, v);         // mirror for LU
        }
        return true;
    }
};

// ---------------------------------------------------------------------------
class InteriorPoint
{
public:
    LPResult solve(const StandardLP& P, const IPMOptions& opt = {})
    {
        LPResult res;
        if (!P.valid())
        {
            res.status  = Status::NumericalError;
            res.message = "invalid problem dimensions";
            return res;
        }

        StopWatch total;
        _recordPath = opt.recordPath;

        const int m = P.rows();
        const int n = P.colsTotal();
        const SparseColMatrix& A = P.A;
        const std::vector<double>& b = P.b;
        const std::vector<double>& c = P.c;

        const double normB = 1.0 + infNorm(b);
        const double normC = 1.0 + infNorm(c);

        NormalEqFactor neq;
        neq.prepare(A, opt.regularize, opt.useCholesky);

        std::vector<double> x(n), s(n), y(m);
        std::vector<double> rp(m), rd(n), rxs(n);
        std::vector<double> dxAff(n), dsAff(n), dyAff(m);
        std::vector<double> dx(n), ds(n), dy(m);
        std::vector<double> d(n), rhs(m), tmpN(n), tmpM(m);

        // ---- Mehrotra starting point -------------------------------------
        // x~ = A^T (A A^T)^{-1} b ,  y~ = (A A^T)^{-1} A c ,  s~ = c - A^T y~
        {
            std::vector<double> ones(n, 1.0);
            StopWatch sw;
            if (!neq.factorize(ones)) // M = A A^T (+reg)
            {
                res.status  = Status::NumericalError;
                res.message = neq.lastError();
                res.stats.msTotal = total.ms();
                return res;
            }
            res.stats.msFactorize += sw.ms();
            res.stats.factorizations++;

            sw.restart();
            neq.solve(b.data(), tmpM.data());          // (AA^T)^{-1} b
            A.mulTransVec(tmpM.data(), x.data());      // x~ = A^T(..)

            A.mulVec(c.data(), rhs.data());            // A c
            neq.solve(rhs.data(), y.data());           // y~
            res.stats.msSolves += sw.ms();

            A.mulTransVec(y.data(), tmpN.data());
            for (int j = 0; j < n; ++j)
                s[j] = c[j] - tmpN[j];                 // s~

            shiftToPositive(x, s);
        }

        recordSnapshot(res, P, x, 0);

        // ---- main loop -----------------------------------------------------
        int it = 0;
        for (; it < opt.maxIterations; ++it)
        {
            // residuals
            A.mulVec(x.data(), tmpM.data());
            for (int i = 0; i < m; ++i)
                rp[i] = b[i] - tmpM[i];
            A.mulTransVec(y.data(), tmpN.data());
            for (int j = 0; j < n; ++j)
                rd[j] = c[j] - tmpN[j] - s[j];

            double mu = 0.0;
            for (int j = 0; j < n; ++j)
                mu += x[j] * s[j];
            mu /= n;

            const double objP = dot(c, x);
            const double relP = infNorm(rp) / normB;
            const double relD = infNorm(rd) / normC;
            const double relG = mu / (1.0 + std::fabs(objP));

            if (opt.verbose)
                std::fprintf(stderr, "[ipm] it=%d mu=%.3e rp=%.3e rd=%.3e obj=%.8g %s\n",
                             it, mu, relP, relD, objP,
                             neq.usingCholesky() ? "SPD-LU" : "LU");

            if (relP < opt.tolResidual && relD < opt.tolResidual && relG < opt.tolMu)
            {
                res.status = Status::Optimal;
                break;
            }

            // ---- factorize M = A D A^T,  D = diag(x/s) ---------------------
            for (int j = 0; j < n; ++j)
                d[j] = x[j] / s[j];

            {
                StopWatch sw;
                if (!neq.factorize(d))
                {
                    res.status  = Status::NumericalError;
                    res.message = neq.lastError();
                    break;
                }
                res.stats.msFactorize += sw.ms();
                res.stats.factorizations++;
            }

            // ---- predictor (affine) ----------------------------------------
            for (int j = 0; j < n; ++j)
                rxs[j] = x[j] * s[j];
            if (!newtonStep(A, neq, d, rp, rd, rxs, x, s, dxAff, dyAff, dsAff, res))
            {
                res.status = Status::NumericalError;
                break;
            }

            const double aPaff = maxStep(x, dxAff);
            const double aDaff = maxStep(s, dsAff);

            double muAff = 0.0;
            for (int j = 0; j < n; ++j)
                muAff += (x[j] + aPaff * dxAff[j]) * (s[j] + aDaff * dsAff[j]);
            muAff /= n;

            const double ratio = muAff / std::max(mu, 1e-300);
            const double sigma = ratio * ratio * ratio; // Mehrotra heuristic

            // ---- corrector --------------------------------------------------
            for (int j = 0; j < n; ++j)
                rxs[j] = x[j] * s[j] - sigma * mu + dxAff[j] * dsAff[j];
            if (!newtonStep(A, neq, d, rp, rd, rxs, x, s, dx, dy, ds, res))
            {
                res.status = Status::NumericalError;
                break;
            }

            const double aP = std::min(1.0, opt.stepFraction * maxStep(x, dx));
            const double aD = std::min(1.0, opt.stepFraction * maxStep(s, ds));

            if (aP < 1e-12 && aD < 1e-12)
            {
                res.status  = Status::NumericalError;
                res.message = "step length collapsed (mu = " + std::to_string(mu) + ")";
                break;
            }

            for (int j = 0; j < n; ++j) x[j] += aP * dx[j];
            for (int i = 0; i < m; ++i) y[i] += aD * dy[i];
            for (int j = 0; j < n; ++j) s[j] += aD * ds[j];

            res.stats.iterations++;
            recordSnapshot(res, P, x, it + 1);
        }

        if (it >= opt.maxIterations && res.status == Status::NotSolved)
            res.status = Status::IterationLimit;

        // ---- pack result ----------------------------------------------------
        res.x.reserve((size_t)n);
        res.y.reserve((size_t)m);
        for (int j = 0; j < n; ++j) res.x[j] = x[j];
        for (int i = 0; i < m; ++i) res.y[i] = y[i];
        res.objective = dot(c, x);
        res.stats.msTotal = total.ms();
        return res;
    }

private:
    static double infNorm(const std::vector<double>& v)
    {
        double s = 0.0;
        for (double e : v)
            s = std::max(s, std::fabs(e));
        return s;
    }

    static double dot(const std::vector<double>& a, const std::vector<double>& v)
    {
        double s = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            s += a[i] * v[i];
        return s;
    }

    // max alpha in (0,1] with v + alpha*dv >= 0
    static double maxStep(const std::vector<double>& v, const std::vector<double>& dv)
    {
        double a = 1.0;
        for (size_t j = 0; j < v.size(); ++j)
            if (dv[j] < 0.0)
                a = std::min(a, -v[j] / dv[j]);
        return a;
    }

    // Mehrotra's positivity shift for the starting point.
    static void shiftToPositive(std::vector<double>& x, std::vector<double>& s)
    {
        const int n = (int)x.size();
        double minX = 0.0, minS = 0.0;
        for (int j = 0; j < n; ++j)
        {
            minX = std::min(minX, x[j]);
            minS = std::min(minS, s[j]);
        }
        const double dx = std::max(-1.5 * minX, 0.0) + 1e-2;
        const double ds = std::max(-1.5 * minS, 0.0) + 1e-2;
        for (int j = 0; j < n; ++j) { x[j] += dx; s[j] += ds; }

        double xs = 0.0, sumX = 0.0, sumS = 0.0;
        for (int j = 0; j < n; ++j)
        {
            xs   += x[j] * s[j];
            sumX += x[j];
            sumS += s[j];
        }
        if (xs <= 0.0) xs = (double)n;
        const double dx2 = 0.5 * xs / std::max(sumS, 1e-300);
        const double ds2 = 0.5 * xs / std::max(sumX, 1e-300);
        for (int j = 0; j < n; ++j) { x[j] += dx2; s[j] += ds2; }
    }

    // Solve one Newton system given the current factorization of A D A^T:
    //   M dy  = rp + A D (rd + X^{-1} rxs)
    //   ds    = rd - A^T dy
    //   dx    = -D ds - S^{-1} rxs
    bool newtonStep(const SparseColMatrix& A, NormalEqFactor& neq,
                    const std::vector<double>& d,
                    const std::vector<double>& rp, const std::vector<double>& rd,
                    const std::vector<double>& rxs,
                    const std::vector<double>& x, const std::vector<double>& s,
                    std::vector<double>& dx, std::vector<double>& dy,
                    std::vector<double>& ds, LPResult& res)
    {
        const int m = A.nRows, n = A.nCols;
        std::vector<double> t(n), rhs(m);

        for (int j = 0; j < n; ++j)
            t[j] = d[j] * (rd[j] + rxs[j] / x[j]);
        A.mulVec(t.data(), rhs.data());
        for (int i = 0; i < m; ++i)
            rhs[i] += rp[i];

        StopWatch sw;
        if (!neq.solve(rhs.data(), dy.data()))
        {
            res.message = "normal-equation solve failed";
            return false;
        }
        res.stats.msSolves += sw.ms();

        A.mulTransVec(dy.data(), t.data());
        for (int j = 0; j < n; ++j)
        {
            ds[j] = rd[j] - t[j];
            dx[j] = -d[j] * ds[j] - rxs[j] / s[j];
        }
        return true;
    }

    void recordSnapshot(LPResult& res, const StandardLP& P,
                        const std::vector<double>& x, int it)
    {
        if (!_recordPath)
            return;
        IterSnapshot snap;
        snap.iter  = it;
        snap.phase = 0;
        snap.x     = x;
        snap.objective = 0.0;
        for (size_t j = 0; j < P.c.size(); ++j)
            snap.objective += P.c[j] * x[j];
        res.path.push_back(std::move(snap));
    }

    bool _recordPath = false;
};

} // namespace lp
