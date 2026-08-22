#pragma once

#include <cstdint>
#include <cstring>

#include "InteractionAlignment.h"
#include "RunningMax.h"
#include "avx2/primitives.h"
#include "cli/cli.h"
#include "int16_safety.h"
#include "nucleotide.h" /* GAP, NEGINF */
#include "operations.h"
#include "optimization/QueryProfile.h"
#include "string_util.h"


// For sequences that iterate target[n - 1 - i]..
// we can point start to target[n - 1] and iterate regularly
// also type safety helps with confusion.
struct ReversedSequence {
    const unsigned char* start;

    unsigned char operator[](int i) const
    {
        return start[-i];
    }
};

enum class TraceState { TRACE_M = 0, TRACE_IX = 1, TRACE_IY = 2, TRACE_DONE = 3 };

/* '|' for a Watson-Crick pair, '.' for a G:U wobble, ' ' for a mismatch.
   The indices are chosen so a pair sums to 3 and a wobble to 5. */
static char pair_symbol(unsigned char q, unsigned char t)
{
    if (q + t == 3) {
        return '|';
    }
    if (q + t == 5) {
        return '.';
    }
    return ' ';
}

static void emit_pair(IA* hit, int l, const unsigned char* query, ReversedSequence target, int i,
                      int j)
{
    const auto qn = query[i];
    const auto tn = target[j];
    hit->ali_seq1[l] = index2nt(qn);
    hit->ali_ia[l] = pair_symbol(qn, tn);
    hit->ali_seq2[l] = index2nt(tn);
}

static void emit_query_bulge(IA* hit, int l, const unsigned char* q, int i)
{
    hit->ali_seq1[l] = index2nt(q[i]);
    hit->ali_ia[l] = ' ';
    hit->ali_seq2[l] = '-';
}

static void emit_target_bulge(IA* hit, int l, ReversedSequence t, int j)
{
    hit->ali_seq1[l] = '-';
    hit->ali_ia[l] = ' ';
    hit->ali_seq2[l] = index2nt(t[j]);
}


/**
 * The old version was transposed. When transposing we unlocked performance, but changed
 * the order of the "best hits", so here we transpose it back without losing major performance.
 */
static RunningMax transpose_best_cell(ReversedSequence target_seq, int m, int n, std::int32_t** M,
                                      const QueryProfile<std::int32_t>& profile, int q_offset,
                                      const std::int32_t* best, const RunningVectorMax& first_row)
{
    RunningVectorMax inverted_best{};
    inverted_best.set(best[1], 1);
    for (auto i = 2; i <= m; i++) {
        inverted_best.set_if_better(best[i], i);
    }

    /* Target position 1 came first, so it keeps ties. */
    RunningMax running_max{first_row.score, first_row.pos_i, 1};

    // if inverted_best.score not better, return the running max of the first row
    if (!running_max.test(inverted_best.score)) {
        return running_max;
    }

    const auto qp = q_offset + inverted_best.pos_i;
    for (auto j = 2; j <= n; j++) {
        const auto ctx = QueryProfile<std::int32_t>::context(target_seq[j - 2], target_seq[j - 1]);
        if (M[j][inverted_best.pos_i] + profile.row(ctx).close[qp] == inverted_best.score) {
            running_max.set(inverted_best.score, inverted_best.pos_i, j);
            break;
        }
    }
    return running_max;
}

/* Target position 1 of the window fill, which the row loop cannot state.
 *
 * The recurrence a row runs reads the target position above it, and there is
 * none here, so this row follows different rules: M can only open on the pair,
 * Iy needs a previous target position and is unreachable, and Ix is the only
 * state left with a recurrence -- along the query. The terms come from the
 * (GAP, target_seq[0]) context, which is how a row with no previous target
 * nucleotide is spelled.
 *
 * best starts empty here rather than at the fill's first row, because best
 * only ever takes a max and every later row folds into it.
 *
 * first_row takes the row's best score and the column it was reached at; this
 * row has one target position, so the position is a column and nothing more.
 */
static void fill_first_row(ReversedSequence target_seq, int m, std::int32_t** M, std::int32_t** Ix,
                           std::int32_t** Iy, const QueryProfile<std::int32_t>& profile,
                           int q_offset, std::int32_t* best, RunningVectorMax& first_row)
{
    const auto T = profile.row(QueryProfile<std::int32_t>::context(GAP, target_seq[0]));
    const std::int32_t* const ix_from_m_1 = T.ix_from_m;
    const std::int32_t* const ix_ext = profile.ix_extend();

    M[1][1] = T.m_open[q_offset + 1];
    best[1] = NEGINF;
    first_row.set(M[1][1] + T.close[q_offset + 1], 1);

    /* The (1,1) cell can be in neither bulge state. */
    Ix[1][1] = Iy[1][1] = NEGINF;

    for (auto i = 2; i <= m; i++) {
        const auto qp = q_offset + i;
        M[1][i] = T.m_open[qp];
        best[i] = NEGINF;
        first_row.set_if_better(M[1][i] + T.close[qp], i);

        Ix[1][i] = max3(M[1][i - 1] != 0 ? M[1][i - 1] + ix_from_m_1[qp] : -1,
                        Ix[1][i - 1] != 0 ? Ix[1][i - 1] + ix_ext[qp] : -1, 0);

        Iy[1][i] = NEGINF;
    }
}

/* Fills the three matrices and, per query column, the best M + close seen in any
 * target position. RIs backtracks through what this leaves behind, so every cell
 * is kept rather than two rows as the sweep keeps.
 */
static void ris_fill_scalar(ReversedSequence target_seq, int m, int n, std::int32_t** M,
                            std::int32_t** Ix, std::int32_t** Iy,
                            const QueryProfile<std::int32_t>& profile, int q_offset,
                            std::int32_t* best, RunningVectorMax& first_row)
{
    const std::int32_t* const ix_ext = profile.ix_extend();

    fill_first_row(target_seq, m, M, Ix, Iy, profile, q_offset, best, first_row);

    const auto qp_first = q_offset + 1;

    for (auto j = 2; j <= n; j++) {
        const auto t =
            profile.row(QueryProfile<std::int32_t>::context(target_seq[j - 2], target_seq[j - 1]));
        const auto iy_ext = t.iy_extend;

        M[j][1] = t.m_open[qp_first];
        best[1] = MAX(best[1], M[j][1] + t.close[qp_first]);

        Ix[j][1] = NEGINF;

        Iy[j][1] = max3(M[j - 1][1] != 0 ? M[j - 1][1] + t.iy_from_m[qp_first] : -1,
                        Iy[j - 1][1] != 0 ? Iy[j - 1][1] + iy_ext : -1, 0);

        for (auto i = 2; i <= m; i++) {
            const auto qp = q_offset + i;

            M[j][i] = max4(
                // continue from a pair
                M[j - 1][i - 1] != 0 ? M[j - 1][i - 1] + t.m_from_m[qp] : -1,
                // close a bulge in query
                Ix[j - 1][i - 1] != 0 ? Ix[j - 1][i - 1] + t.m_from_ix[qp] : -1,
                // close a bulge in target
                Iy[j - 1][i - 1] != 0 ? Iy[j - 1][i - 1] + t.m_from_iy[qp] : -1,
                // start fresh
                t.m_open[qp]);

            best[i] = MAX(best[i], M[j][i] + t.close[qp]);


            Ix[j][i] = max3(
                // Open a bulge in query
                M[j][i - 1] != 0 ? M[j][i - 1] + t.ix_from_m[qp] : -1,
                // Continue a bulge in query
                Ix[j][i - 1] != 0 ? Ix[j][i - 1] + ix_ext[qp] : -1, 0);

            Iy[j][i] = max3(
                // Open a bulge in target
                M[j - 1][i] != 0 ? M[j - 1][i] + t.iy_from_m[qp] : -1,
                // Continue a bulge in target
                Iy[j - 1][i] != 0 ? Iy[j - 1][i] + iy_ext : -1, 0);
        }
    }
}


#if RISEARCH1_HAS_AVX2

/* The same fill, eight query positions at a time. It is meant to be read beside
 * ris_fill_scalar: every statement below is one of that function's with each
 * scalar replaced by eight, and the differences are only these.
 *
 *  - Target position 1 and query position 1 stay scalar. The first has no
 *    previous target nt for its terms, the second no column to its left for a
 *    block to read.
 *
 *  - M and Iy read only the previous target position, so eight of their columns
 *    are independent and a block computes them outright. Blocks cover columns
 *    2..m with the last backed up to end at m; it recomputes the columns it
 *    overlaps to the same values, and best only ever takes a max, so overlapping
 *    changes nothing.
 *
 *  - Ix and Iy drop their tests for a predecessor of exactly 0. Both end in a max
 *    with 0, so the -1 a failed test yields and the bulge it avoids are both <= 0
 *    and neither can win. That holds only while no bulge term is positive, which
 *    is what has_positive_gap() rules out and what the dispatch below requires. M
 *    keeps its tests: its fourth candidate is an opening score rather than 0, so
 *    -1 can still win there.
 *
 *  - Ix reads the column to its left in the row it writes, so it is rewritten as
 *    a running max over candidates with ix_prefix taken out -- the same rewrite
 *    the sweep uses. A running max cannot revisit a column already folded into
 *    its carry, so only whole blocks go through it and the last columns finish
 *    serially.
 */
__attribute__((target("avx2"))) static void
ris_fill_avx2(ReversedSequence target_seq, int m, int n, std::int32_t** M, std::int32_t** Ix,
              std::int32_t** Iy, const QueryProfile<std::int32_t>& profile, int q_offset,
              std::int32_t* best, RunningVectorMax& first_row)
{
    constexpr auto kBlock = static_cast<int>(v_lanes<std::int32_t>());

    fill_first_row(target_seq, m, M, Ix, Iy, profile, q_offset, best, first_row);

    /* The single block a window shorter than one block runs writes columns 2
       through kBlock + 1, and its running max reads best where it is about to
       write, so those columns start from the same place the others do. */
    for (auto i = m + 1; i <= kBlock + 1; i++) {
        best[i] = NEGINF;
    }

    const auto qp_first = q_offset + 1;
    const __m256i zero = v_zero_to_avx2();

    for (auto j = 2; j <= n; j++) {
        const auto t =
            profile.row(QueryProfile<std::int32_t>::context(target_seq[j - 2], target_seq[j - 1]));
        const auto iy_ext = t.iy_extend;

        const std::int32_t* const m_prev = M[j - 1];
        const std::int32_t* const ix_prev = Ix[j - 1];
        const std::int32_t* const iy_prev = Iy[j - 1];
        std::int32_t* const m_cur = M[j];
        std::int32_t* const ix_cur = Ix[j];
        std::int32_t* const iy_cur = Iy[j];

        m_cur[1] = t.m_open[qp_first];
        best[1] = MAX(best[1], m_cur[1] + t.close[qp_first]);

        ix_cur[1] = NEGINF;

        iy_cur[1] = max3(m_prev[1] != 0 ? m_prev[1] + t.iy_from_m[qp_first] : -1,
                         iy_prev[1] != 0 ? iy_prev[1] + iy_ext : -1, 0);

        /* Offset once, so a block's terms are indexed by the same i as its rows. */
        const std::int32_t* const from_m = t.m_from_m + q_offset;
        const std::int32_t* const from_ix = t.m_from_ix + q_offset;
        const std::int32_t* const from_iy = t.m_from_iy + q_offset;
        const std::int32_t* const open = t.m_open + q_offset;
        const std::int32_t* const close = t.close + q_offset;
        const std::int32_t* const iy_from_m = t.iy_from_m + q_offset;
        const std::int32_t* const ix_scan = t.ix_from_m_scan + q_offset;
        const std::int32_t* const ix_pref = t.ix_prefix + q_offset;

        const __m256i v_iy_ext = v_int_to_avx2<std::int32_t>(iy_ext);

        for (auto start = 2; start <= m; start += kBlock) {
            /* The clamp is what makes the last block end exactly at m, and what
               keeps a window shorter than one block on the single block that
               starts at column 2. That block reaches past m, into the slack the
               profile keeps; nothing reads those columns, since a column is read
               only by the row below it and the column to its right. */
            const auto i = MIN(start, MAX(2, m - (kBlock - 1)));

            /* For each column the block writes, its diagonal predecessor: one
               target position back and one query position back. */
            const __m256i m_diag = v_vec_load<std::int32_t>(m_prev + i - 1);
            const __m256i ix_diag = v_vec_load<std::int32_t>(ix_prev + i - 1);
            const __m256i iy_diag = v_vec_load<std::int32_t>(iy_prev + i - 1);

            const __m256i m_new = v_max4<std::int32_t>(
                // continue from a pair
                v_add_unless_zero_or_neg1<std::int32_t>(m_diag,
                                                        v_vec_load<std::int32_t>(from_m + i)),
                // close a bulge in query
                v_add_unless_zero_or_neg1<std::int32_t>(ix_diag,
                                                        v_vec_load<std::int32_t>(from_ix + i)),
                // close a bulge in target
                v_add_unless_zero_or_neg1<std::int32_t>(iy_diag,
                                                        v_vec_load<std::int32_t>(from_iy + i)),
                // start fresh
                v_vec_load<std::int32_t>(open + i));
            v_vec_store<std::int32_t>(m_cur + i, m_new);

            v_vec_store<std::int32_t>(
                best + i, v_max<std::int32_t>(
                              v_vec_load<std::int32_t>(best + i),
                              v_add<std::int32_t>(m_new, v_vec_load<std::int32_t>(close + i))));

            /* Iy's predecessors are vertical: previous target position, same
               column, so no -1 on the address. */
            v_vec_store<std::int32_t>(
                iy_cur + i,
                v_max3<std::int32_t>(
                    // Open a bulge in target
                    v_add<std::int32_t>(v_vec_load<std::int32_t>(m_prev + i),
                                        v_vec_load<std::int32_t>(iy_from_m + i)),
                    // Continue a bulge in target
                    v_add<std::int32_t>(v_vec_load<std::int32_t>(iy_prev + i), v_iy_ext), zero));
        }

        /* Ix, as a running max over the same candidates with ix_prefix taken
           out. The 0 arm becomes -ix_prefix in that space; the carry starts
           below every candidate, and Ix[j][1] is unreachable so it cannot win. */
        __m256i carry = v_int_to_avx2<std::int32_t>(NEGINF);
        auto i = 2;
        for (; i + kBlock <= m + 1; i += kBlock) {
            const __m256i prefix = v_vec_load<std::int32_t>(ix_pref + i);
            const __m256i candidates = v_max<std::int32_t>( // Open a bulge in query
                v_add<std::int32_t>(v_vec_load<std::int32_t>(m_cur + i - 1),
                                    v_vec_load<std::int32_t>(ix_scan + i)),
                // or hold nothing yet, which is this column's 0
                v_sub<std::int32_t>(zero, prefix));
            const __m256i winner =
                v_max<std::int32_t>(v_prefix_max<std::int32_t>(candidates), carry);
            /* Putting ix_prefix back turns the carried quantity into the real Ix,
               which is what the next column and the backtrack read. */
            v_vec_store<std::int32_t>(ix_cur + i, v_add<std::int32_t>(winner, prefix));
            carry = v_broadcast_last<std::int32_t>(winner);
        }
        /* The columns left over go through one more block, backed up to end at m.
           A running max cannot revisit a column already folded into its carry, so
           the carry is rebuilt from the column just below the block: taking
           ix_prefix back out of the Ix stored there is what the scan would have
           been carrying. */
        if (i <= m) {
            const auto back = MAX(2, m - (kBlock - 1));
            const __m256i prefix = v_vec_load<std::int32_t>(ix_pref + back);
            const __m256i candidates =
                v_max<std::int32_t>(v_add<std::int32_t>(v_vec_load<std::int32_t>(m_cur + back - 1),
                                                        v_vec_load<std::int32_t>(ix_scan + back)),
                                    v_sub<std::int32_t>(zero, prefix));
            const __m256i winner = v_max<std::int32_t>(
                v_prefix_max<std::int32_t>(candidates),
                v_int_to_avx2<std::int32_t>(ix_cur[back - 1] - ix_pref[back - 1]));
            v_vec_store<std::int32_t>(ix_cur + back, v_add<std::int32_t>(winner, prefix));
        }
    }
}

#endif

/* The vector fill needs eight columns to the right of column 1 for a block, and
   it may only drop Ix's and Iy's zero tests where no bulge term is positive. */
/* A window shorter than one block still runs one, which writes past m into
   columns the matrices hold but no window of this size fills. The matrices are
   square in the traceback length, and both m and n are bounded by it, so the
   longer of the two says whether a whole block fits. */
static bool ris_fill_is_vectorized(int m, int n, const QueryProfile<std::int32_t>& profile)
{
#if RISEARCH1_HAS_AVX2
    return m >= 2 && MAX(m, n) >= static_cast<int>(v_lanes<std::int32_t>()) + 1 &&
           !profile.has_positive_gap() && CPU_HAS_AVX2;
#else
    (void)m;
    (void)n;
    (void)profile;
    return false;
#endif
}

static void ris_fill(ReversedSequence target_seq, int m, int n, std::int32_t** M, std::int32_t** Ix,
                     std::int32_t** Iy, const QueryProfile<std::int32_t>& profile, int q_offset,
                     std::int32_t* best, RunningVectorMax& first_row)
{
#if RISEARCH1_HAS_AVX2
    if (ris_fill_is_vectorized(m, n, profile)) {
        ris_fill_avx2(target_seq, m, n, M, Ix, Iy, profile, q_offset, best, first_row);
        return;
    }
#endif
    ris_fill_scalar(target_seq, m, n, M, Ix, Iy, profile, q_offset, best, first_row);
}

static void
RIs(const unsigned char* query_seq, /* query sequence - numeric representation */
    ReversedSequence target_seq,    /* target sequence - reversed */
    int m,                          /* query seq length */
    int n,                          /* target seq length */
    IA* hit,                        /* pointer to struct, fill results */
    const config_st& config, std::int32_t** M, std::int32_t** Ix, std::int32_t** Iy,
    const QueryProfile<std::int32_t>& profile, int q_offset,
    std::int32_t* best // We use `best` parameter to preserve output order of old C version
)
{
    /* The matrices are indexed [target][query] so that a row fixes the target
       context and the terms it reads are runs over consecutive query positions. */
    /* The backtrack re-derives which arm produced each cell, so it reads the
       same extension terms the fill did. */
    const std::int32_t* const ix_ext = profile.ix_extend();

    RunningVectorMax first_row{};
    ris_fill(target_seq, m, n, M, Ix, Iy, profile, q_offset, best, first_row);
    const auto running_max =
        transpose_best_cell(target_seq, m, n, M, profile, q_offset, best, first_row);


    /*backtrack*/
    const auto capacity = max_alignment_length(config.tblen);

    auto i = running_max.pos_i;
    auto j = running_max.pos_j;
    auto k = TraceState::TRACE_M;
    auto l = 0; /* alilen so far -used in backtrack */

    while ((i > 0 && j > 0) && (M[j][i] > 0 || Ix[j][i] > 0 || Iy[j][i] > 0)) {
        if (l > capacity) {
            printf("Interaction longer than max, so the following is only the end of the full "
                   "alignment:\n");
            /*alt: stop here / reallocate? prevent creation of longer alignments in first place? */
            break;
        }


        /*l++ in end of while instead of every sub? */

        const auto qp = q_offset + i;


        switch (k) {
        case TraceState::TRACE_M: {
            if (const auto open_score =
                    profile.row(QueryProfile<std::int32_t>::context(GAP, target_seq[j - 1]))
                        .m_open[qp];
                M[j][i] == open_score) {
                k = TraceState::TRACE_DONE;
            } else {
                const auto& t = profile.row(
                    QueryProfile<std::int32_t>::context(target_seq[j - 2], target_seq[j - 1]));

                if (M[j][i] == M[j - 1][i - 1] + t.m_from_m[qp]) {
                    k = TraceState::TRACE_M;
                } else if (M[j][i] == Ix[j - 1][i - 1] + t.m_from_ix[qp]) {
                    k = TraceState::TRACE_IX;
                } else if (M[j][i] == Iy[j - 1][i - 1] + t.m_from_iy[qp]) {
                    k = TraceState::TRACE_IY;
                } else {
                    printf("unexpected value in k=0.\n");
                }
            }
            emit_pair(hit, l++, query_seq, target_seq, --i, --j);
            break;
        }
        case TraceState::TRACE_IX: {
            // in this case, t_prev doesn't matter so we put it as GAP
            // avoids bugs in the j == 1 case.
            /* seq1(query) paired to a gap (in target) */
            const auto gap_row =
                profile.row(QueryProfile<std::int32_t>::context(GAP, target_seq[j - 1]));
            if (Ix[j][i] == M[j][i - 1] + gap_row.ix_from_m[qp]) {
                k = TraceState::TRACE_M; /* open a new gap coming from match */
            } else if (Ix[j][i] == Ix[j][i - 1] + ix_ext[qp]) {
                k = TraceState::TRACE_IX; /* extend existing gap */
            } else if (Ix[j][i] ==
                       profile.row(QueryProfile<std::int32_t>::context(GAP, GAP)).m_open[qp]) {
                fprintf(stderr, "\nErr: This alignment starts in a gap - not even an option!?\n");
                k = TraceState::TRACE_DONE; /* start new alignment with gap; not possible, prevented
                                               by scoring... */
            } else {
                printf("unexpected case in k=1 : %d\n", Ix[j][i]);
            }
            emit_query_bulge(hit, l++, query_seq, --i);
            break;
        }
        case TraceState::TRACE_IY: {
            const auto context =
                QueryProfile<std::int32_t>::context(target_seq[j - 2], target_seq[j - 1]);

            /* seq2(target) paired to a gap (in query) */
            if (Iy[j][i] == M[j - 1][i] + profile.row(context).iy_from_m[qp]) {
                k = TraceState::TRACE_M; /* open a new gap coming from match */
            } else if (Iy[j][i] == Iy[j - 1][i] + profile.row(context).iy_extend) {
                k = TraceState::TRACE_IY; /* extend existing gap */
            } else if (Iy[j][i] ==
                       profile.row(QueryProfile<std::int32_t>::context(GAP, target_seq[j - 1]))
                           .iy_extend) {
                k = TraceState::TRACE_DONE;
                fprintf(stderr, "\nErr: This alignment starts in a gap - not even an option!?\n");
            } else {
                printf("unexpected case in k=2 : %d\n", Iy[j][i]);
            }
            emit_target_bulge(hit, l++, target_seq, --j);
            break;
        }
        default: {
            fprintf(stderr, "\nThis should really NEVER happen!\n");
            break;
        }
        }

        if (TraceState::TRACE_DONE == k) {
            break;
        }
    }
    hit->ali_seq1[l] = '\0';
    hit->ali_ia[l] = '\0';
    hit->ali_seq2[l] = '\0';

    /* reverse sequences in the end*/

    if (l > 0) { // fixes a pre-existing bug, access to invalid memory [-1]
        reverse_inplace(hit->ali_seq1.get(), l - 1);
        reverse_inplace(hit->ali_ia.get(), l - 1);
        reverse_inplace(hit->ali_seq2.get(), l - 1);
    }

    hit->qbeg = i + 1;
    hit->qend = running_max.pos_i;
    hit->tbeg = n + 1 - running_max.pos_j;
    hit->tend = n - j;
    hit->max = running_max.score;
}
