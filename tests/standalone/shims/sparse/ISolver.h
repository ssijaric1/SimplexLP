// SimplexLP standalone test shim - mirrors the sparse::ISolver interface of
// natID (same method names, same factory signature, same accumulate
// semantics for addTriple after populateDiagonals) so the LP core compiles
// unchanged. The factorization behind it is a straightforward dense LU with
// partial pivoting - perfectly adequate for unit-test problem sizes.
//
// The real application links natID's Matrix library instead of this file.
#pragma once
#include <cnt/SafeFullVector.h>
#include <mem/PointerReleaser.h>
#include <vector>
#include <cmath>
#include <complex>
#include <cstring>
#include <ostream>

namespace td { typedef int INT4; }

#define MATRIX_API

namespace sparse
{

enum class SolverType : unsigned char { LU = 0, LLT, LDLT, Ctrl, Pardiso };
enum class Pivoting   : unsigned char { No = 0, DiagonalSinglePass, DiagonalMultiPass,
                                        MarkowitzSinglePass, MarkowitzMultiPass,
                                        AlterMatrixIfIndefinite };
enum class Ordering   : unsigned char { Own = 0, OwnRadial, ExternSym, ExternAtPlusA,
                                        ExternAtMulA, ExternAMulAt };
enum class Symmetry   : unsigned char { NonSymmetric = 0, SymmetricPosDef, SymmetricGeneral,
                                        SymmetricIndef, HermitianPosDef, HermitianIndef };
enum class Format     : unsigned char { Matlab = 0, Fortran };

template <typename TVAL, typename INDEX>
class ISolver
{
public:
    virtual void addTriple(INDEX i, INDEX j, const TVAL& val) = 0;
    virtual void setRHS(INDEX i, const TVAL& val) = 0;
    virtual const TVAL& x(INDEX i) const = 0;
    virtual void populateDiagonals(const TVAL& val) = 0;
    virtual void clearRHS() = 0;
    virtual INDEX size() const = 0;
    virtual bool factorize() = 0;
    virtual const char* getLastError() const = 0;
    virtual bool solve() = 0;
    virtual bool solveExt(const TVAL* rhs, TVAL* x) = 0;
    virtual void release() = 0;
    virtual int getNoOfNonZero() const = 0;
    virtual ~ISolver() = default;
};

typedef ISolver<double, td::INT4>               DblSolver;
typedef ISolver<std::complex<double>, td::INT4> CmplxSolver;
typedef mem::PointerReleaser<DblSolver>   DblSolverReleaser;
typedef mem::PointerReleaser<CmplxSolver> CmplxSolverReleaser;

// ------------------------------------------------------- test-only backend --
namespace shim
{

class DenseLUSolver : public DblSolver
{
    int _n = 0;
    std::vector<double> _A;    // accumulated input, row-major n x n
    std::vector<double> _LU;   // factors
    std::vector<int>    _perm;
    std::vector<double> _rhs, _x;
    bool _factorized = false;
    const char* _err = nullptr;

    // counts symmetry mode so a lower-triangle-only insert (Cholesky style)
    // still factorizes correctly through the dense LU
    bool _symmetricMirror = false;

public:
    DenseLUSolver(int n, bool symmetricMirror)
    : _n(n)
    , _A((size_t)n * n, 0.0)
    , _perm(n)
    , _rhs(n, 0.0)
    , _x(n, 0.0)
    , _symmetricMirror(symmetricMirror)
    {
    }

    void addTriple(td::INT4 i, td::INT4 j, const double& val) override
    {
        _A[(size_t)i * _n + j] += val; // accumulate, like the real solver
        if (_symmetricMirror && i != j)
            _A[(size_t)j * _n + i] += val;
        _factorized = false;
    }

    void populateDiagonals(const double& val) override
    {
        for (int i = 0; i < _n; ++i)
            _A[(size_t)i * _n + i] += val;
    }

    void setRHS(td::INT4 i, const double& val) override { _rhs[i] = val; }
    const double& x(td::INT4 i) const override { return _x[i]; }
    void clearRHS() override { std::fill(_rhs.begin(), _rhs.end(), 0.0); }
    td::INT4 size() const override { return _n; }

    int getNoOfNonZero() const override
    {
        int nz = 0;
        for (double v : _A)
            if (v != 0.0)
                ++nz;
        return nz;
    }

    bool factorize() override
    {
        _LU = _A;
        for (int i = 0; i < _n; ++i)
            _perm[i] = i;

        for (int k = 0; k < _n; ++k)
        {
            int    piv  = k;
            double pval = std::fabs(_LU[(size_t)k * _n + k]);
            for (int i = k + 1; i < _n; ++i)
            {
                const double v = std::fabs(_LU[(size_t)i * _n + k]);
                if (v > pval) { pval = v; piv = i; }
            }
            if (pval < 1e-13)
            {
                _err = "singular matrix";
                return false;
            }
            if (piv != k)
            {
                for (int j = 0; j < _n; ++j)
                    std::swap(_LU[(size_t)k * _n + j], _LU[(size_t)piv * _n + j]);
                std::swap(_perm[k], _perm[piv]);
            }
            const double dkk = _LU[(size_t)k * _n + k];
            for (int i = k + 1; i < _n; ++i)
            {
                const double f = _LU[(size_t)i * _n + k] / dkk;
                _LU[(size_t)i * _n + k] = f;
                if (f != 0.0)
                    for (int j = k + 1; j < _n; ++j)
                        _LU[(size_t)i * _n + j] -= f * _LU[(size_t)k * _n + j];
            }
        }
        _factorized = true;
        _err = nullptr;
        return true;
    }

    const char* getLastError() const override { return _err; }

    bool solveExt(const double* rhs, double* x) override
    {
        if (!_factorized)
            return false;
        std::vector<double> ywork(_n);
        for (int i = 0; i < _n; ++i)
            ywork[i] = rhs[_perm[i]];
        for (int i = 1; i < _n; ++i)
        {
            double s = ywork[i];
            for (int j = 0; j < i; ++j)
                s -= _LU[(size_t)i * _n + j] * ywork[j];
            ywork[i] = s;
        }
        for (int i = _n - 1; i >= 0; --i)
        {
            double s = ywork[i];
            for (int j = i + 1; j < _n; ++j)
                s -= _LU[(size_t)i * _n + j] * x[j];
            x[i] = s / _LU[(size_t)i * _n + i];
        }
        return true;
    }

    bool solve() override { return solveExt(_rhs.data(), _x.data()); }
    void release() override { delete this; }
};

} // namespace shim

inline DblSolver* createDblSolver(int NMAT, int /*NZ*/,
                                  Symmetry sym = Symmetry::NonSymmetric,
                                  SolverType /*solverType*/ = SolverType::LU,
                                  Pivoting /*pivoting*/ = Pivoting::DiagonalSinglePass,
                                  Ordering /*ordering*/ = Ordering::Own)
{
    // a SymmetricPosDef insert provides only one triangle; mirror it so the
    // dense LU sees the full matrix
    const bool mirror = (sym == Symmetry::SymmetricPosDef ||
                         sym == Symmetry::SymmetricGeneral ||
                         sym == Symmetry::SymmetricIndef);
    return new shim::DenseLUSolver(NMAT, mirror);
}

} // namespace sparse
