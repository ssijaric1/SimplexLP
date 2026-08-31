//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  BasisFactor.h
//  Wraps natID's sparse LU (sparse::createDblSolver, SolverType::LU) around
//  the simplex basis matrix B. The revised simplex needs two kinds of
//  triangular solves per iteration:
//
//      FTRAN:  B  w = a_q      (entering column, ratio test)
//      FTRAN:  B  x_B = b      (current basic solution)
//      BTRAN:  B^T y = c_B     (simplex multipliers / pricing)
//
//  sparse::ISolver factors one matrix and offers solveExt() for arbitrary
//  right-hand sides, but has no transpose solve, so a second solver instance
//  is loaded with B^T. Both share the same sparsity work per iteration and
//  keep everything inside natID, as required by the project proposal.
//
//  Solvers are created fresh at every (re)factorization - the ISolver
//  interface has no "reset values" call, and this mirrors how the SDK's own
//  MatrixTests use the factory (create -> populate -> factorize -> solve).
#pragma once

#include "LPCommon.h"
#include <sparse/ISolver.h>

namespace lp
{

class BasisFactor
{
    sparse::DblSolverReleaser _sB;   // factor of B
    sparse::DblSolverReleaser _sBT;  // factor of B^T
    int _m = 0;
    std::string _lastError;

public:
    BasisFactor() : _sB(nullptr), _sBT(nullptr) {}

    const std::string& lastError() const { return _lastError; }

    // Load and factorize B = A(:, basis) and its transpose.
    // 'A' is the (row-sign normalized) working matrix.
    bool factorize(const SparseColMatrix& A, const std::vector<int>& basis)
    {
        _m = (int)basis.size();
        _lastError.clear();

        int nnzB = 0;
        for (int k = 0; k < _m; ++k)
            nnzB += (int)A.column(basis[k]).size();

        // +_m: populateDiagonals() seeds the diagonal pattern as in the SDK
        // examples; a small headroom keeps the allocator comfortable.
        const int reserveNZ = nnzB + _m + 8;
        const double zero = 0.0;

        _sB = sparse::createDblSolver(_m, reserveNZ,
                                      sparse::Symmetry::NonSymmetric,
                                      sparse::SolverType::LU,
                                      sparse::Pivoting::MarkowitzSinglePass);
        _sBT = sparse::createDblSolver(_m, reserveNZ,
                                       sparse::Symmetry::NonSymmetric,
                                       sparse::SolverType::LU,
                                       sparse::Pivoting::MarkowitzSinglePass);
        if (!_sB || !_sBT)
        {
            _lastError = "createDblSolver returned null";
            return false;
        }

        _sB->populateDiagonals(zero);
        _sBT->populateDiagonals(zero);

        for (int k = 0; k < _m; ++k)
        {
            for (const SparseEntry& e : A.column(basis[k]))
            {
                _sB->addTriple(e.row, k, e.val);   // B(i,k)
                _sBT->addTriple(k, e.row, e.val);  // B^T(k,i)
            }
        }

        if (!_sB->factorize())
        {
            const char* err = _sB->getLastError();
            _lastError = std::string("LU(B) failed: ") + (err ? err : "?");
            return false;
        }
        if (!_sBT->factorize())
        {
            const char* err = _sBT->getLastError();
            _lastError = std::string("LU(B^T) failed: ") + (err ? err : "?");
            return false;
        }
        return true;
    }

    // w = B^{-1} rhs
    bool ftran(const double* rhs, double* w) { return _sB->solveExt(rhs, w); }

    // y = B^{-T} rhs
    bool btran(const double* rhs, double* y) { return _sBT->solveExt(rhs, y); }
};

} // namespace lp
