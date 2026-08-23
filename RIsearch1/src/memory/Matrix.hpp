#pragma once

#include <cstdint>
#include <cstdlib>
#include <utility>

#include "memory/alloc.hpp"

/**
 * A rows x cols matrix, allocated as a single block: the row pointers first,
 * the values straight after them.
 *
 * One allocation rather than one per row, and the values contiguous, which is
 * what lets a fill walk a row and a vector kernel load one.
 *
 * Move-only. A copy would have to allocate a second block, and every matrix
 * here is scratch that one owner fills and reads.
 */
template<typename T>
class Matrix {
public:
    Matrix(std::uint32_t rows, std::uint32_t cols)
    {
        const auto pointer_bytes = static_cast<std::size_t>(rows) * sizeof(T*);
        const auto data_bytes = static_cast<std::size_t>(rows) * cols * sizeof(T);

        m_data = static_cast<T**>(std::malloc(pointer_bytes + data_bytes));
        if (m_data == nullptr) {
            out_of_memory("an alignment matrix", pointer_bytes + data_bytes);
        }

        T* values = reinterpret_cast<T*>(m_data + rows);
        for (auto i = 0u; i < rows; i++) {
            m_data[i] = values + static_cast<std::size_t>(i) * cols;
        }
    }

    ~Matrix()
    {
        std::free(m_data);
    }

    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;

    Matrix(Matrix&& other) noexcept : m_data(std::exchange(other.m_data, nullptr))
    {
    }

    Matrix& operator=(Matrix&& other) noexcept
    {
        if (this != &other) {
            std::free(m_data);
            m_data = std::exchange(other.m_data, nullptr);
        }
        return *this;
    }

    /* One row, so that m[j][i] reads as it does on a plain array. */
    T* operator[](std::size_t row) const
    {
        return m_data[row];
    }

    /* The row pointers, for the fills and the backtrack, which take T**. */
    T** get() const
    {
        return m_data;
    }

private:
    T** m_data = nullptr;
};
