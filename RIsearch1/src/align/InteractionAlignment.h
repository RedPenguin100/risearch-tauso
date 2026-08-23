#pragma once

#include <cstdint>

#include "memory/MallocRAII.hpp"

/* How much room the three alignment strings get.
 *
 * An alignment can come out longer than the window it was recovered from: every
 * bulge adds a position that only one of the two sequences has. Half again the
 * traceback length is the allowance.
 */
inline int alignment_capacity(std::uint32_t tblen)
{
    return static_cast<int>(1.5 * tblen);
}

/* The longest alignment a backtrack may write. It writes one character at l and
 * the terminator at l + 1 afterwards, so it has to stop two short of the end.
 */
inline int max_alignment_length(std::uint32_t tblen)
{
    return alignment_capacity(tblen) - 2;
}

struct IA {
    int qbeg, qend, tbeg, tend;
    int max;

    MallocRAII<char> ali_seq1, ali_seq2, ali_ia;

    explicit IA(int capacity)
        : qbeg(0), qend(0), tbeg(0), tend(0), max(0), ali_seq1(capacity), ali_seq2(capacity),
          ali_ia(capacity)
    {
    }

    [[nodiscard]] int nucleotide_count() const
    {
        return qend - qbeg + 1 + tend - tbeg + 1;
    }
};
