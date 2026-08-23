#pragma once

#include <cstdint>

#include "RunningMax.h"
#include "align/int16_safety.h"
#include "dsm.h"
#include "nucleotide.h" /* GAP */
#include "operations.h"

/* Target position 1 of a sweep, which no sweep kernel computes.
 *
 * The recurrence a sweep states reads the target position before the current
 * one, and at the first there is none, so this row follows different rules: M
 * can only open on the pair, Iy needs a previous target position and is
 * unreachable, and Ix is the only state left with a recurrence -- along the
 * query rather than down the target.
 *
 * M, Ix and Iy are one row each, indexed by query position, column 0 included.
 * t_last is the target's 3' nucleotide, which is where a sweep begins, since it
 * walks the target backwards. Returned is the row's best score and the query
 * column it was reached at; the row has one target position, so the position is
 * a column and nothing more.
 *
 * Both sweeps start from this. The single-query one writes it into the first of
 * its two rows; the batched one runs it per lane in int32 and converts what it
 * gets into the terms its kernel enters with.
 *
 * The window fill has a first row of its own, in traceback.h: same recurrence,
 * but it reads profile terms at a query offset and seeds the per-column best
 * that the backtrack later reads.
 */
template<typename int_type>
static RunningVectorMax score_first_row(int_type* M, int_type* Ix, int_type* Iy,
                                        const unsigned char* query, std::uint32_t m,
                                        unsigned char t_last, const Dsm& dsm)
{
    /* Column 0: an alignment must open on a base pair, never on a bulge. */
    Iy[0] = 0;
    Ix[0] = M[0] = neg_inf<int_type>();

    M[1] = dsm[GAP][query[0]][GAP][t_last];
    /* The (1,1) cell can be in neither bulge state. */
    Ix[1] = Iy[1] = neg_inf<int_type>();

    RunningVectorMax row_max{};
    row_max.set(M[1] + dsm[query[0]][GAP][t_last][GAP], 1);

    for (auto i = 2u; i <= m; i++) {
        const auto q_prev = query[i - 2];
        const auto q_cur = query[i - 1];

        M[i] = dsm[GAP][q_cur][GAP][t_last];
        row_max.set_if_better(M[i] + dsm[q_cur][GAP][t_last][GAP], static_cast<int>(i));

        /* One option is missing on purpose: starting a new alignment already in
           a gap, which is (-, Xi; -, -). */
        Ix[i] = max3(0,
                     /* pair at i - 1, now a gap -- add (Xi-1, Xi; Y1, -) */
                     M[i - 1] != 0 ? M[i - 1] + dsm[q_prev][q_cur][t_last][GAP] : -1,
                     /* already bulging, add one more -- add (Xi-1, Xi; -, -) */
                     Ix[i - 1] != 0 ? Ix[i - 1] + dsm[q_prev][q_cur][GAP][GAP] : -1);

        /* There is no previous target position, so there can be no bulge. */
        Iy[i] = neg_inf<int_type>();
    }
    return row_max;
}
