#pragma once

#include <climits>
#include <cstdint>
#include <utility>

#include "RunningMax.h"
#include "align/avx2/primitives.h"
#include "align/int16_safety.h"
#include "optimization/BatchedQueryProfile.h"
#include "optimization/QueryProfile.h"

/* The sweep of ScoreTarget.h with one query per lane instead of one query spread
 * across the lanes.
 *
 * The kernels there put sixteen query positions of a single query in a register,
 * so Ix reads the row it writes and has to be carried as a prefix max, and M and
 * Iy need blocks whose geometry fits the query's length. At the lengths this
 * program is used at most of those lanes hold nothing.
 *
 * Sixteen different queries fill every lane whatever their length, and they are
 * independent: nothing here crosses lanes, so there is no prefix max, no block
 * geometry and no length-dependent dispatch. Every statement below is one of
 * score_target_scalar's with each scalar replaced by sixteen.
 *
 * IX IS NEVER STORED. Ix[j] is read in one place only, by the row below it and
 * one column to the left, and its own recurrence runs along the row. With
 * ix_prefix taken out it is a plain running max whose inputs at column i are
 * M[j][i - 1] and a term -- both of which the row below is already reading -- so
 * that row runs the max itself, one column behind, and carries it in a register.
 * What ix_prefix owes is folded into m_from_ix, which is where the profile keeps
 * it.
 */

#if RISEARCH1_HAS_AVX2

/* RunningMax for sixteen queries at once. The row number needs more than a short,
   so it is held as two halves of eight ints. */
struct BatchedRunningMax {
    __m256i score;
    __m256i pos_i;
    __m256i pos_j_lo;
    __m256i pos_j_hi;
};

/* A row's scores are its columns' M plus their close; the largest of them, read
   back out of the row the sweep has just written. Max is associative, so taking
   the columns in any order gives the same answer as accumulating them. */
__attribute__((target("avx2"), always_inline)) static inline __m256i
row_max_exact(const std::int16_t* m_cur, const BatchedQueryProfile::RowView& T, unsigned m)
{
    constexpr auto queries = BatchedQueryProfile::kLanes;
    __m256i row_max =
        v_add<std::int16_t>(v_vec_load(m_cur + queries), v_vec_load(T.column(1).close));
    for (auto i = 2u; i <= m; i++) {
        row_max = v_max<std::int16_t>(row_max, v_add<std::int16_t>(v_vec_load(m_cur + i * queries),
                                                                   v_vec_load(T.column(i).close)));
    }
    return row_max;
}

/* The column each query's row max was reached at, first one wins. Walking the
   columns backwards makes "first" fall out of the order the blends happen in: a
   lower column overwrites what a higher one selected, and the row max is a max
   over exactly these candidates, so every query matches somewhere. */
__attribute__((target("avx2"), always_inline)) static inline __m256i
row_max_position(const std::int16_t* m_cur, const BatchedQueryProfile::RowView& T, unsigned m,
                 __m256i row_max)
{
    constexpr auto queries = BatchedQueryProfile::kLanes;
    __m256i pos = v_int_to_avx2<std::int16_t>(1);
    for (auto i = m; i >= 1; i--) {
        const __m256i candidate =
            v_add<std::int16_t>(v_vec_load(m_cur + i * queries), v_vec_load(T.column(i).close));
        pos = v_select(pos, v_int_to_avx2<std::int16_t>(i),
                       v_equals<std::int16_t>(candidate, row_max));
    }
    return pos;
}

/* M and Iy for one column, sixteen queries at a time, and the column's
   contribution to the row max. Mirrors main_dp_loop_avx2, which does the same
   for eight columns of one query. */
template<bool kDefer>
__attribute__((target("avx2"), always_inline)) static inline void
main_dp_column_batched(const BatchedQueryProfile::ColumnTerms& t, __m256i m_last_prev,
                       __m256i iy_last_prev, __m256i m_last_here, __m256i iy_last_here,
                       __m256i ix_scan, __m256i iy_ext, std::int16_t* m_out, std::int16_t* iy_out,
                       __m256i* row_max)
{
    const __m256i m_new = v_max4<std::int16_t>(
        /* coming from a match, and its M[lastRow][i-1] != 0 test */
        v_add_unless_zero_or_neg1<std::int16_t>(m_last_prev, v_vec_load(t.m_from_m)),
        /* coming from gap in target, as the row above's scan left it */
        v_add<std::int16_t>(ix_scan, v_vec_load(t.m_from_ix)),
        /* coming from gap in query */
        v_add<std::int16_t>(iy_last_prev, v_vec_load(t.m_from_iy)),
        /* start fresh */
        v_vec_load(t.m_open));

    v_vec_store(m_out, m_new);

    // set max now, position is recovered later.
    *row_max = v_max<std::int16_t>(
        *row_max, kDefer ? m_new : v_add<std::int16_t>(m_new, v_vec_load(t.close)));

    // Iy: target nt against a gap. Predecessors are vertical -- previous row,
    // same column -- so no -1 on the address, unlike M's diagonal above.
    v_vec_store(iy_out, v_max<std::int16_t>(
                            /* pair at previous row, now bulge */
                            v_add<std::int16_t>(m_last_here, v_vec_load(t.iy_from_m)),
                            /* already bulging, add one more */
                            v_add<std::int16_t>(iy_last_here, iy_ext)));
}

/* M and Iy hold two target positions each, as they do in the single query sweep.
   Ix is not among them: see IX IS NEVER STORED above.

   WHAT ROW 1 HANDS OVER. The scan a row runs on behalf of the row above it is

       U[j][i] = max(U[j][i - 1], M[j][i - 1] + ix_from_m_scan[i]),

   and the terms it reads come from the profile, at the target nucleotide of the
   row the bulge is in. Row 1 is the exception: it is not written by this kernel
   at all. The caller computes it in int32 under the recurrence that governs a
   first row -- where a predecessor of exactly 0 is refused and the whole thing
   is floored at 0 -- which is not the recurrence this scan states. So the caller
   converts what it computed into terms of its own choosing, such that running
   the scan above over them reproduces U[1][i] exactly. That is what arrives here
   as ix_from_m_scan_row1, and it is read only when j == 2.

   run_scores and run_positions take target position j at [(j - 1) * queries]. */
/* kDefer: a row's close is added once at the end as a bound rather than per
   column. Not available where the reporting reads a non-clearing row's score,
   which is what a vicinity window does. */
template<bool kDefer>
__attribute__((target("avx2"))) static void
score_target_batched(const unsigned char* target_sequence, const BatchedQueryProfile& profile,
                     std::int16_t* const* M, std::int16_t* const* Iy,
                     const std::int16_t* ix_from_m_scan_row1, std::int16_t* run_scores,
                     std::int16_t* run_positions, std::uint16_t* row_clears, std::size_t n,
                     int threshold, BatchedRunningMax& running_max)
{
    const auto m = profile.m();
    constexpr auto queries = BatchedQueryProfile::kLanes;
    constexpr auto solo_group = BatchedQueryProfile::kSoloGroup;

    std::int16_t* m_cur = M[0];
    std::int16_t* m_last = M[1];
    std::int16_t* iy_cur = Iy[0];
    std::int16_t* iy_last = Iy[1];

    const __m256i zero = v_zero_to_avx2();
    const __m256i one = v_int_to_avx2<std::int16_t>(1);

    /* A row's position only has to be right where it is read, which is where the
       row max clears the threshold or improves on the query's best so far -- the
       same test the single query kernels defer their position scan behind. */
    const __m256i v_threshold = v_int_to_avx2<std::int16_t>(static_cast<std::int16_t>(
        threshold > SHRT_MAX ? SHRT_MAX : (threshold < SHRT_MIN ? SHRT_MIN : threshold)));

    /* A row is wanted where it clears the threshold or improves on the query's
       best, which is one comparison against the smaller of the two. The best only
       rises, and only where a row was wanted, so this rises with it. */
    __m256i gate = v_min<std::int16_t>(running_max.score, v_threshold);

    for (auto j = 2u; j <= n; j++) {
        const auto target_current = target_sequence[n - j];
        const auto target_prev = target_sequence[n - j + 1];

        const auto context = QueryProfile<std::int32_t>::context(target_prev, target_current);
        const auto T = profile.row(context, target_current);
        const __m256i iy_ext = v_int_to_avx2<std::int16_t>(T.iy_extend);

        /* The run the row above's scan reads. Row 2's is the caller's; every
           later row's sits at the target nucleotide of the row above. */
        const auto scan_stride = j == 2 ? queries : solo_group;
        /* Already at column 1's term, so that the step at the top of the column
           loop leaves it on column i. */
        const std::int16_t* scan =
            (j == 2 ? ix_from_m_scan_row1 : profile.scan_terms(target_prev)) + scan_stride;

        const auto t1 = T.column(1);

        /* Column 1 is the query's first nt, nothing can precede it. */
        const __m256i m_col1 = v_max<std::int16_t>(zero, v_vec_load(t1.m_open));
        v_vec_store(m_cur + queries, m_col1);

        /* Deferring, this is a row max over M alone; the close each column owes
           it is put back once, as a bound, when the row is finished. */
        __m256i row_max = kDefer ? m_col1 : v_add<std::int16_t>(m_col1, v_vec_load(t1.close));

        __m256i m_last_prev = v_vec_load(m_last + queries);
        __m256i iy_last_prev = v_vec_load(iy_last + queries);

        /* Iy bulges a target nt, j >= 2 so a bulge is possible. */
        v_vec_store(iy_cur + queries,
                    v_max<std::int16_t>(v_add<std::int16_t>(m_last_prev, v_vec_load(t1.iy_from_m)),
                                        v_add<std::int16_t>(iy_last_prev, iy_ext)));

        /* Ix[j-1][1] is unreachable, and carried that is what the row above's
           scan enters its first column with. */
        __m256i ix_scan = v_int_to_avx2<std::int16_t>(NEG_INF_SHORT);

        const std::int16_t* pair_at = T.pair + 2 * BatchedQueryProfile::kPairGroup;
        const std::int16_t* solo_at = T.solo + 2 * BatchedQueryProfile::kSoloGroup;

        for (auto i = 2u; i < m; i++, pair_at += BatchedQueryProfile::kPairGroup,
                  solo_at += BatchedQueryProfile::kSoloGroup) {
            const auto off = i * queries;
            scan += scan_stride;

            const __m256i m_last_here = v_vec_load(m_last + off);
            const __m256i iy_last_here = v_vec_load(iy_last + off);

            main_dp_column_batched<kDefer>(BatchedQueryProfile::RowView::column_at(pair_at, solo_at),
                                           m_last_prev, iy_last_prev, m_last_here,
                                           iy_last_here, ix_scan, iy_ext, m_cur + off, iy_cur + off,
                                           &row_max);

            /* One step of the row above's Ix scan, a column behind this row so
               that column i is ready for column i + 1, where it is read. */
            ix_scan =
                v_max<std::int16_t>(ix_scan, v_add<std::int16_t>(m_last_prev, v_vec_load(scan)));

            m_last_prev = m_last_here;
            iy_last_prev = iy_last_here;
        }

        /* THE LAST COLUMN NEEDS ONLY M. Iy[j][m] is read by exactly one thing,
           M's target-bulge arm at column m + 1, and by its own recurrence, which
           the same thing reads one row further down; there being no column m + 1,
           the whole chain leaves nothing behind. The Ix scan's step at column m
           is the value column m + 1 would enter with, and goes for the same
           reason. */
        if (m >= 2) {
            const auto t = T.column(m);
            const __m256i m_new = v_max4<std::int16_t>(
                v_add_unless_zero_or_neg1<std::int16_t>(m_last_prev, v_vec_load(t.m_from_m)),
                v_add<std::int16_t>(ix_scan, v_vec_load(t.m_from_ix)),
                v_add<std::int16_t>(iy_last_prev, v_vec_load(t.m_from_iy)), v_vec_load(t.m_open));

            v_vec_store(m_cur + m * queries, m_new);

            row_max = v_max<std::int16_t>(
                row_max, kDefer ? m_new : v_add<std::int16_t>(m_new, v_vec_load(t.close)));
        }

        /* Deferring, what the row loop has is the largest M of the row; adding
           the largest close any of its columns could have owed bounds the row's
           score from above, since a saturating add is monotone in both arms. A
           row whose bound clears neither the threshold nor the query's best is
           read nowhere, so the bound stands in for its score. */
        const __m256i bound =
            kDefer ? v_add<std::int16_t>(row_max, v_vec_load(profile.close_max(target_current)))
                   : row_max;

        if (!kDefer || v_any(v_greater_than<std::int16_t>(bound, gate))) {
            const __m256i exact = kDefer ? row_max_exact(m_cur, T, m) : row_max;
            const __m256i improved = v_greater_than<std::int16_t>(exact, running_max.score);
            const __m256i wanted = v_greater_than<std::int16_t>(exact, gate);
            const __m256i row_pos = v_any(wanted) ? row_max_position(m_cur, T, m, exact) : one;

            v_vec_store(run_scores + (j - 1) * queries, exact);
            v_vec_store(run_positions + (j - 1) * queries, row_pos);
            row_clears[j - 1] = static_cast<std::uint16_t>(
                v_lane_bits16(v_greater_than<std::int16_t>(exact, v_threshold)));

            running_max.score = v_max<std::int16_t>(running_max.score, exact);
            gate = v_min<std::int16_t>(running_max.score, v_threshold);
            running_max.pos_i = v_select(running_max.pos_i, row_pos, improved);
            /* The row number does not fit a short, so the sixteen query mask is
               widened to two eight wide ones to select it. */
            const __m256i v_j = v_int_to_avx2<std::int32_t>(static_cast<std::int32_t>(j));
            running_max.pos_j_lo = v_select(running_max.pos_j_lo, v_j, v_widen_low(improved));
            running_max.pos_j_hi = v_select(running_max.pos_j_hi, v_j, v_widen_high(improved));
        } else {
            /* The position is read only where a row reported, and such a row went
               the other way, so this one's is left as it lies. */
            v_vec_store(run_scores + (j - 1) * queries, bound);
            /* The gate is at most the threshold and this row did not clear the
               gate, so no lane of it clears the threshold either. */
            row_clears[j - 1] = 0;
        }

        /* The row just written becomes the row read. */
        std::swap(m_cur, m_last);
        std::swap(iy_cur, iy_last);
    }
}

#endif /* RISEARCH1_HAS_AVX2 */
