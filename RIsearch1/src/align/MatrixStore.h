#pragma once

#include <cstdint>

#include "align/int16_safety.h"
#include "memory/Matrix.hpp"
#include "nucleotide.h" /* NEGINF */


/* Row 0 and column 0 hold the same values whatever the window, so they are set
   once for the whole matrix and no fill writes them. */
static void ris_fill_bounds(std::int32_t** M, std::int32_t** Ix, std::int32_t** Iy, int rows,
                            int cols)
{
    M[0][0] = Ix[0][0] = Iy[0][0] = 0;

    // Target position 0, every query position
    for (auto i = 1; i < cols; i++) {
        Iy[0][i] = M[0][i] = NEGINF; /* not possible before beginning of target seq */
        Ix[0][i] = 0;
    }

    // Query position 0, every target position
    for (auto j = 1; j < rows; j++) {
        Ix[j][0] = M[j][0] = NEGINF;
        Iy[j][0] = 0;
    }
}


class MatrixStore {
public:
    explicit MatrixStore(int traceback_len)
        : m_M(traceback_len + 1, traceback_len + 1), m_Ix(traceback_len + 1, traceback_len + 1),
          m_Iy(traceback_len + 1, traceback_len + 1)
    {
        ris_fill_bounds(m_M.get(), m_Ix.get(), m_Iy.get(), traceback_len + 1, traceback_len + 1);
    }

    MatrixStore(const MatrixStore&) = delete;
    MatrixStore& operator=(const MatrixStore&) = delete;

    std::int32_t** M() const
    {
        return m_M.get();
    }
    std::int32_t** Ix() const
    {
        return m_Ix.get();
    }
    std::int32_t** Iy() const
    {
        return m_Iy.get();
    }

private:
    Matrix<std::int32_t> m_M;
    Matrix<std::int32_t> m_Ix;
    Matrix<std::int32_t> m_Iy;
};
