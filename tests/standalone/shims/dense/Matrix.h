// SimplexLP standalone test shim - mirrors the slice of natID's
// dense::Matrix API that Polytope.h calls: construction, element access via
// getManipulator()/getFirstColumnManipulator(), and the in-place solve()
// where the right-hand-side matrix is overwritten by the solution.
#pragma once
#include <vector>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace dense
{

template <typename T> class Matrix;

template <typename T>
class MatrixIOProxy
{
    Matrix<T>* _m;
public:
    explicit MatrixIOProxy(Matrix<T>* m) : _m(m) {}
    inline T& operator()(unsigned i, unsigned j);
};

template <typename T>
class FirstColumnIOProxy
{
    Matrix<T>* _m;
public:
    explicit FirstColumnIOProxy(Matrix<T>* m) : _m(m) {}
    inline T& operator()(unsigned i);
};

template <typename T>
class Matrix
{
    unsigned _rows = 0, _cols = 0;
    std::vector<T> _d; // row-major

public:
    Matrix() = default;
    Matrix(unsigned nRows, unsigned nCols, void* /*pMgr*/ = nullptr, bool initZero = false)
    : _rows(nRows), _cols(nCols), _d((size_t)nRows * nCols, T{})
    {
        (void)initZero; // shim always zero-initializes
    }

    unsigned getNoOfRows() const { return _rows; }
    unsigned getNoOfCols() const { return _cols; }

    void zeros() { std::fill(_d.begin(), _d.end(), T{}); }

    T&       at(unsigned i, unsigned j)       { assert(i < _rows && j < _cols); return _d[(size_t)i * _cols + j]; }
    const T& at(unsigned i, unsigned j) const { assert(i < _rows && j < _cols); return _d[(size_t)i * _cols + j]; }

    MatrixIOProxy<T>      getManipulator()            { return MatrixIOProxy<T>(this); }
    FirstColumnIOProxy<T> getFirstColumnManipulator() { return FirstColumnIOProxy<T>(this); }

    // In-place solve: on input BX holds B of A x = B, on output it holds X.
    // Gaussian elimination with partial pivoting on a copy of *this.
    bool solve(Matrix<T>& BX)
    {
        if (_rows != _cols || BX.getNoOfRows() != _rows)
            return false;
        const unsigned n = _rows, k = BX.getNoOfCols();
        std::vector<T> A = _d;

        for (unsigned c = 0; c < n; ++c)
        {
            unsigned piv = c;
            double   pv  = std::fabs((double)A[(size_t)c * n + c]);
            for (unsigned i = c + 1; i < n; ++i)
            {
                const double v = std::fabs((double)A[(size_t)i * n + c]);
                if (v > pv) { pv = v; piv = i; }
            }
            if (pv < 1e-13)
                return false;
            if (piv != c)
            {
                for (unsigned j = 0; j < n; ++j)
                    std::swap(A[(size_t)c * n + j], A[(size_t)piv * n + j]);
                for (unsigned j = 0; j < k; ++j)
                    std::swap(BX.at(c, j), BX.at(piv, j));
            }
            const T d = A[(size_t)c * n + c];
            for (unsigned i = c + 1; i < n; ++i)
            {
                const T f = A[(size_t)i * n + c] / d;
                if (f != T{})
                {
                    for (unsigned j = c; j < n; ++j)
                        A[(size_t)i * n + j] -= f * A[(size_t)c * n + j];
                    for (unsigned j = 0; j < k; ++j)
                        BX.at(i, j) -= f * BX.at(c, j);
                }
            }
        }
        for (int i = (int)n - 1; i >= 0; --i)
        {
            for (unsigned j = 0; j < k; ++j)
            {
                T s = BX.at((unsigned)i, j);
                for (unsigned t = (unsigned)i + 1; t < n; ++t)
                    s -= A[(size_t)i * n + t] * BX.at(t, j);
                BX.at((unsigned)i, j) = s / A[(size_t)i * n + i];
            }
        }
        return true;
    }
};

template <typename T>
inline T& MatrixIOProxy<T>::operator()(unsigned i, unsigned j) { return _m->at(i, j); }

template <typename T>
inline T& FirstColumnIOProxy<T>::operator()(unsigned i) { return _m->at(i, 0); }

typedef Matrix<double> DblMatrix;

} // namespace dense
