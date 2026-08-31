// SimplexLP standalone test shim - mirrors natID's mem::PointerReleaser.
#pragma once

namespace mem
{

template <typename T>
class PointerReleaser
{
protected:
    T* _ptr;

public:
    PointerReleaser(T* ptr) : _ptr(ptr) {}
    PointerReleaser() : _ptr(nullptr) {}

    ~PointerReleaser()
    {
        if (_ptr)
            _ptr->release();
    }

    inline void operator=(T* ptr)
    {
        if (_ptr)
            _ptr->release();
        _ptr = ptr;
    }

    inline T*       operator->()       { return _ptr; }
    inline const T* operator->() const { return _ptr; }
    inline bool     operator!()  const { return _ptr == nullptr; }
    inline T&       ref()              { return *_ptr; }
    inline T*       ptr()              { return _ptr; }
};

} // namespace mem
