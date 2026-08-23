#pragma once

#include <cstddef>
#include <cstdlib>

#include "memory/alloc.hpp"

template<typename T>
class MallocRAII {
public:
    MallocRAII() : m_buffer(nullptr)
    {
    }

    explicit MallocRAII(size_t n) : m_buffer(static_cast<T*>(malloc(n * sizeof(T))))
    {
        if (m_buffer == nullptr && n != 0) {
            out_of_memory("a buffer", n * sizeof(T));
        }
    }

    ~MallocRAII()
    {
        free(m_buffer);
    }

    MallocRAII(const MallocRAII&) = delete;
    MallocRAII& operator=(const MallocRAII&) = delete;
    MallocRAII(MallocRAII&& rhs) noexcept : m_buffer(rhs.m_buffer)
    {
        rhs.m_buffer = nullptr;
    }
    MallocRAII& operator=(MallocRAII&& rhs) noexcept
    {
        if (this != &rhs) {
            free(m_buffer);
            m_buffer = rhs.m_buffer;
            rhs.m_buffer = nullptr;
        }
        return *this;
    }

    T* operator->() const
    {
        return m_buffer;
    }

    T& operator[](size_t i) const
    {
        return m_buffer[i];
    }

    T* get() const
    {
        return m_buffer;
    }

    void reset()
    {
        free(m_buffer);
        m_buffer = nullptr;
    }

    void reset(T* buffer)
    {
        free(m_buffer);
        m_buffer = buffer;
    }


private:
    T* m_buffer;
};
