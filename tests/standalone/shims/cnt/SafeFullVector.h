// SimplexLP standalone test shim - mirrors the small API surface of
// natID's cnt::SafeFullVector that the LP core uses. Semantics copied from
// the real header: reserve(size) (re)allocates and SETS the size; element
// data starts uninitialized; zeros() clears; clean() frees.
//
// This file is used ONLY by tests/standalone - the application builds
// against the real natID SDK headers.
#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>

namespace cnt
{

template <typename T>
class SafeFullVector
{
    T*     _data = nullptr;
    size_t _size = 0;

public:
    SafeFullVector() = default;
    SafeFullVector(const SafeFullVector& o) { *this = o; }
    SafeFullVector& operator=(const SafeFullVector& o)
    {
        if (this == &o)
            return *this;
        reserve(o._size);
        for (size_t i = 0; i < _size; ++i)
            _data[i] = o._data[i];
        return *this;
    }
    ~SafeFullVector() { clean(); }

    inline void clean()
    {
        delete[] _data;
        _data = nullptr;
        _size = 0;
    }

    inline void reserve(size_t size)
    {
        if (_size == size)
            return;
        delete[] _data;
        _data = (size > 0) ? new T[size] : nullptr;
        _size = size;
    }

    inline size_t size() const { return _size; }

    inline const T& operator[](size_t pos) const { assert(pos < _size); return _data[pos]; }
    inline T&       operator[](size_t pos)       { assert(pos < _size); return _data[pos]; }
    inline const T& operator()(size_t pos) const { assert(pos < _size); return _data[pos]; }
    inline T&       operator()(size_t pos)       { assert(pos < _size); return _data[pos]; }

    inline void zeros() { std::fill(_data, _data + _size, T{}); }

    T* begin() { return _data; }
    T* end()   { return _data + _size; }
};

} // namespace cnt
