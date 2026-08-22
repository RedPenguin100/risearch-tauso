#pragma once

#include <climits>
#include <cstdint>
#include <cstdlib>
#include <utility>

#include "RunningMax.h"
#include "int16_safety.h"
#include "nucleotide.h"
#include "operations.h"
#include "optimization/QueryProfile.h"

/* Scoring the target: the whole j loop over target positions, run once per
 * alignment.
 *
 * For every target position it records the best score of an alignment ending
 * there and the query position it ended at -- run_scores and run_positions below -- and keeps the
 * single best across the whole target. It reconstructs nothing; that is RIs.
 *
 * score_target_scalar is the reference and the definition of what the output
 * must be. score_target_avx2 computes the same recurrence eight query positions
 * at a time and has to agree with it bit for bit, so the two are meant to be
 * read side by side.
 *
 * BOTH KEEP TWO ROWS AS A PAIR OF POINTERS, swapped at the end of each target
 * position. Row j belongs in row j % 2, so j = 2 writes row 0 and the swap says
 * the same thing, while keeping the row addresses in registers rather than
 * reloading them from the array every row.
 */

#include "avx2/primitives.h"


/* What a row leaves behind once its columns are done: the run entries at this
 * target position, and the query's best so far.
 *
 * The column the row max was reached at is recovered by rescanning the row
 * rather than tracked alongside the max, because it is read only where the row
 * clears the threshold or improves on the query's best -- which almost no row
 * does, while every row would have paid for tracking it.
 */
template<typename int_type>
__attribute__((always_inline)) static inline void
record_row(const int_type* m_row, const typename QueryProfile<int_type>::RowView& T, unsigned m,
           int row_max, unsigned j, int threshold, int* run_scores, int* run_positions,
           RunningMax& running_max)
{
    auto row_pos = 1u;
    if (row_max > threshold || row_max > running_max.score) {
        for (auto i = 1u; i <= m; ++i) {
            if (m_row[i] + T.close[i] == row_max) {
                row_pos = i;
                break;
            }
        }
    }

    run_scores[j - 1] = row_max;
    run_positions[j - 1] = static_cast<int>(row_pos);
    running_max.set_if_better(row_max, static_cast<int>(row_pos), static_cast<int>(j));
}


/* Inlined into the caller on purpose: there the rows and the profile's tables
   are visibly separate allocations, and out of line the compiler has to assume a
   store through one could land in the other and re-load everything per row. */
template<typename int_type>
__attribute__((always_inline)) static inline void
score_target_scalar(const unsigned char* target_sequence, const QueryProfile<int_type>& profile,
                    int_type* const* M, int_type* const* Ix, int_type* const* Iy, int* run_scores,
                    int* run_positions, std::size_t n, int threshold, RunningMax& running_max)
{
    const auto m = profile.query_length();


    int_type* m_cur = M[0];
    int_type* m_last = M[1];
    int_type* ix_cur = Ix[0];
    int_type* ix_last = Ix[1];
    int_type* iy_cur = Iy[0];
    int_type* iy_last = Iy[1];

    for (auto j = 2u; j <= n; j++) {
        const auto target_current = target_sequence[n - j];
        const auto target_prev = target_sequence[n - j + 1];

        const auto context = QueryProfile<int_type>::context(target_prev, target_current);
        const auto T = profile.row(context);
        const auto iy_ext = T.iy_extend;


        /* Column 1 is the query's first nt, nothing can precede it */
        m_cur[1] = MAX(0, T.m_open[1]);

        // Track only the best value, the position is recovered later (OPTIMIZATION)
        auto row_max = m_cur[1] + T.close[1];

        // Ix bulges a query nt, impossible in 1st nucleotide
        ix_cur[1] = neg_inf<int_type>();

        // Iy bulges a target nt, j>=2 so bulge possible.
        // Only M can win row max, so we don't take this as a candidate
        iy_cur[1] = MAX(
            // We are now opening the bulge
            m_last[1] + T.iy_from_m[1],
            // We are extending a bulge
            iy_last[1] + iy_ext);

        /* finished init of i=1 col */


        for (auto i = 2u; i <= m; i++) {
            // Assign M with a 4-way max
            m_cur[i] = max4(
                /* coming from a match */
                m_last[i - 1] != 0 ? m_last[i - 1] + T.m_from_m[i] : -1,
                /* coming from gap in target */
                ix_last[i - 1] + T.m_from_ix[i],
                /* coming from gap in query */
                iy_last[i - 1] + T.m_from_iy[i],
                /* start fresh */
                T.m_open[i]);

            // Set max now, position is recovered later
            row_max = MAX(row_max, m_cur[i] + T.close[i]);

            // Iy: target nt against a gap
            iy_cur[i] = MAX(
                // pair at previous row, now bulge
                m_last[i] + T.iy_from_m[i],
                // already bulging, add one more
                iy_last[i] + iy_ext);
        }

        // Ix reads the row it writes, so it gets a pass of its own.
        for (auto i = 2u; i <= m; ++i) {
            // Ix: query nt against a gap
            ix_cur[i] = MAX(
                // pair at i - 1, now bulge
                m_cur[i - 1] + T.ix_from_m[i],
                // already bulging, add one more
                ix_cur[i - 1] + T.ix_extend[i]);
        }

        record_row<int_type>(m_cur, T, m, row_max, j, threshold, run_scores, run_positions,
                             running_max);

        /* The row just written becomes the row read. */
        std::swap(m_cur, m_last);
        std::swap(ix_cur, ix_last);
        std::swap(iy_cur, iy_last);
    } /*next row j */
}


#if RISEARCH1_HAS_AVX2

/* M and Iy for the eight columns starting at i, and their contribution to the
   row max. Returns M, which the Ix scan reads. */
template<typename int_type>
__attribute__((target("avx2"), always_inline)) static inline __m256i
main_dp_loop_avx2(unsigned i, const typename QueryProfile<int_type>::RowView& T, int_type* m_cur,
                  int_type* iy_cur, const int_type* m_last, const int_type* ix_last,
                  const int_type* iy_last, __m256i v_iy_ext, __m256i* v_row_max)
{
    /* For each column this block writes, its diagonal predecessor: one target
       position back and one query position back. Adjacent because the columns
       are adjacent; the diagonal is only the -1. */
    const __m256i m_diag = v_vec_load<int_type>(m_last + i - 1);


    // Assign M with a 4-way max
    const __m256i m_new = v_max4<int_type>(
        /* coming from a match */
        v_add_unless_zero_or_neg1<int_type>(m_diag, v_vec_load<int_type>(T.m_from_m + i)),
        /* coming from gap in target */
        v_add<int_type>(v_vec_load<int_type>(ix_last + i - 1),
                        v_vec_load<int_type>(T.m_from_ix + i)),
        /* coming from gap in query */
        v_add<int_type>(v_vec_load<int_type>(iy_last + i - 1),
                        v_vec_load<int_type>(T.m_from_iy + i)),
        /* start fresh */
        v_vec_load<int_type>(T.m_open + i));
    v_vec_store<int_type>(m_cur + i, m_new);

    // Set max now, position is recovered later
    *v_row_max =
        v_max<int_type>(*v_row_max, v_add<int_type>(m_new, v_vec_load<int_type>(T.close + i)));

    // Iy: target nt against a gap. Predecessors are vertical -- previous row,
    // same column -- so no -1 on the address, unlike M's diagonal above.
    v_vec_store<int_type>(iy_cur + i,
                          v_max<int_type>(
                              /* pair at previous row, now bulge */
                              v_add<int_type>(v_vec_load<int_type>(m_last + i),
                                              v_vec_load<int_type>(T.iy_from_m + i)),
                              /* already bulging, add one more */
                              v_add<int_type>(v_vec_load<int_type>(iy_last + i), v_iy_ext)));

    return m_new;
}

/* The Ix scan for a block's eight columns, starting at column i.
 *
 * m_left holds M at the column just left of each of the eight, which is how the
 * loop hands the scan the M it has in a register instead of putting the row
 * back through memory. The value returned is the carry the next block starts
 * from: lane 7 of the running max, which is its maximum over the whole block,
 * broadcast so that the max above needs no shuffling. */
template<typename int_type>
__attribute__((target("avx2"), always_inline)) static inline __m256i
ix_dp_loop_avx2(int_type* ix_out, const int_type* ix_from_m_scan, const int_type* ix_prefix,
                __m256i m_left, __m256i ix_carry)
{
    // Ix: query nt against a gap
    const __m256i candidates = v_add<int_type>(m_left, v_vec_load<int_type>(ix_from_m_scan));
    const __m256i best = v_max<int_type>(v_prefix_max<int_type>(candidates), ix_carry);
    /* Adding ix_prefix back turns the carried quantity into the real Ix, which
       is what the next row and the traceback read. */
    v_vec_store<int_type>(ix_out, v_add<int_type>(best, v_vec_load<int_type>(ix_prefix)));
    return v_broadcast_last<int_type>(best);
}

/* The same recurrence, eight query positions at a time. Differences from the
 * scalar version above, and nothing else:
 *
 *  - M and Iy read only the previous row, so eight of their columns are
 *    independent and a block computes them outright. The blocks cover columns
 *    2..m, the last backed up to end at m and recomputing the columns it
 *    overlaps -- to the same values, since it reads only the previous row --
 *    which is cheaper than a scalar epilogue. That needs m >= 9, which
 *    score_target_is_vectorized has already checked;
 *  - Ix reads the row it writes, so it is rewritten as a plain running max over
 *    candidates with ix_prefix taken out (see QueryProfile). A running max
 *    cannot revisit a column already folded into its carry, so only whole
 *    blocks go through it and the last few columns finish serially;
 *  - the row max is accumulated in eight lanes and reduced once per row.
 *
 * Column 1 stays scalar throughout: the neighbouring column every block reads
 * does not exist there.
 *
 * THE SCAN TAKES M FROM THE REGISTER THE BLOCK BUILT IT IN, one column to the
 * right, rather than loading the row it has just been stored to. A load that
 * reads back what the same row has just written cannot be served from the store
 * that wrote it, and the scan's candidates are the first thing on the row's
 * serial chain, so that one load sits in front of everything else. Removing it
 * is worth around 30% of this kernel at m = 12 or 16, where a row is one or two
 * blocks wide. At m >= 80 a row is ten blocks and there is enough independent
 * work to cover the load anyway, and the shift that replaces it costs about as
 * much as it saves.
 */
template<typename int_type>
__attribute__((target("avx2"))) static void
score_target_avx2(const unsigned char* target_sequence, const QueryProfile<int_type>& profile,
                  int_type* const* M, int_type* const* Ix, int_type* const* Iy, int* run_scores,
                  int* run_positions, std::size_t n, int threshold, RunningMax& running_max)
{
    const auto m = profile.query_length();

    int_type* m_cur = M[0];
    int_type* m_last = M[1];
    int_type* ix_cur = Ix[0];
    int_type* ix_last = Ix[1];
    int_type* iy_cur = Iy[0];
    int_type* iy_last = Iy[1];

    for (auto j = 2u; j <= n; j++) {
        const auto target_current = target_sequence[n - j];
        const auto target_prev = target_sequence[n - j + 1];

        const auto context = QueryProfile<int_type>::context(target_prev, target_current);
        const auto T = profile.row(context);
        const auto iy_ext = T.iy_extend;

        /* Column 1 is the query's first nt, nothing can precede it */
        m_cur[1] = MAX(0, T.m_open[1]);
        auto row_max = m_cur[1] + T.close[1];
        ix_cur[1] = neg_inf<int_type>();
        iy_cur[1] = MAX(m_last[1] + T.iy_from_m[1], iy_last[1] + iy_ext);
        /* finished init of i=1 col */

        const __m256i v_iy_ext = v_int_to_avx2<int_type>(iy_ext);
        __m256i v_row_max = v_int_to_avx2<int_type>(row_max);

        /* M at the column just left of the scan's first block, and the carry it
           enters with -- what Ix[1] holds, so the first block sees exactly what
           the serial recurrence would have carried into it. */
        __m256i m_left = v_int_to_avx2<int_type>(m_cur[1]);
        __m256i v_ix_carry = v_int_to_avx2<int_type>(neg_inf<int_type>());

        /* Begin main DP */
        constexpr auto kBlock = v_lanes<int_type>();
        auto start = 2u;
        for (; start + kBlock - 1 <= m; start += kBlock) {
            const __m256i m_block = main_dp_loop_avx2<int_type>(
                start, T, m_cur, iy_cur, m_last, ix_last, iy_last, v_iy_ext, &v_row_max);
            // The separate ix loop in a separate function
            v_ix_carry = ix_dp_loop_avx2<int_type>(
                ix_cur + start, T.ix_from_m_scan + start, T.ix_prefix + start,
                v_shifted_left_one<int_type>(m_left, m_block), v_ix_carry);
            m_left = m_block;
        }

        // Fewer than a block's columns need special case for main DP
        if (start <= m) {
            main_dp_loop_avx2<int_type>(m - (kBlock - 1), T, m_cur, iy_cur, m_last, ix_last,
                                        iy_last, v_iy_ext, &v_row_max);
        }
        row_max = v_hmax<int_type>(v_row_max);

        // Fewer than a block's columns need special case for main IX DP
        // Ix can't re-use the function like main_dp because an overlap behaves differently.
        for (auto i = start; i <= m; ++i) {
            ix_cur[i] = MAX(m_cur[i - 1] + T.ix_from_m[i], ix_cur[i - 1] + T.ix_extend[i]);
        }
        /* End main DP */


        record_row<int_type>(m_cur, T, m, row_max, j, threshold, run_scores, run_positions,
                             running_max);

        /* The row just written becomes the row read. */
        std::swap(m_cur, m_last);
        std::swap(ix_cur, ix_last);
        std::swap(iy_cur, iy_last);
    }
}

#endif /* RISEARCH1_HAS_AVX2 */


/* The entry point the alignment calls. Inlined for the same reason
   score_target_scalar is. */
template<typename int_type>
__attribute__((always_inline)) static inline void
score_target(const unsigned char* target_sequence, const QueryProfile<int_type>& profile,
             int_type* const* M, int_type* const* Ix, int_type* const* Iy, int* run_scores,
             int* run_positions, std::size_t n, int threshold, RunningMax& running_max)
{
#if RISEARCH1_HAS_AVX2
    /* A block covers the query positions from 2 upwards, so a query shorter
       than one block has none to fill and the scalar sweep is the cheaper way. */
    if (profile.query_length() <= v_lanes<int_type>()) {
        score_target_scalar<int_type>(target_sequence, profile, M, Ix, Iy, run_scores,
                                      run_positions, n, threshold, running_max);
        return;
    }
    // Only use AVX2 if the CPU supports
    if (CPU_HAS_AVX2) {
        score_target_avx2<int_type>(target_sequence, profile, M, Ix, Iy, run_scores, run_positions,
                                    n, threshold, running_max);
        return;
    }
#endif
    score_target_scalar<int_type>(target_sequence, profile, M, Ix, Iy, run_scores, run_positions, n,
                                  threshold, running_max);
}
