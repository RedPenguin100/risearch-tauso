#pragma once

#include <climits>
#include <cstdint>
#include <cstdio>

#include "HitReporter.h"
#include "RunningMax.h"
#include "align/ScoreTarget.h"
#include "align/first_row.h"
#include "align/int16_safety.h"
#include "align/optimization/QueryProfile.h"
#include "cli/cli.h"
#include "energy.hpp"
#include "memory/ByteBuffer.hpp"
#include "memory/MallocRAII.hpp"
#include "nucleotide.h"
#include "operations.h"

/* Compiled once for AVX2 and once for the baseline, with the CPU checked at load
   time. The sweep and the window fill this calls carry target("avx2"), and gcc
   will not inline a function that has the attribute into a caller that does not,
   so without the clone both stay out-of-line calls and nothing around them is
   optimised across the boundary. Giving them a caller that already has AVX2 is
   what lets them inline, and the baseline clone keeps the binary running where
   AVX2 is missing. */
template<typename int_type>
__attribute__((target_clones("avx2", "default"))) static void
RIs_linSpace(const ByteBuffer& query_sequence_ix,  // query sequence numerical representation
             const ByteBuffer& target_sequence_ix, // target sequence numerical representation
             Dsm& dsm,                             /* scoring matrix -- TODO variable length!? */
             int threshold,                        /* give out hits higher than that */
             const char* qname,                    /* query name */
             const char* tname,                    /* target name */
             const config_st& config)
{
    const auto m = query_sequence_ix.size();
    const auto n = target_sequence_ix.size();
    const auto* target_sequence = target_sequence_ix.unsigned_data();
    const auto* query_sequence = query_sequence_ix.unsigned_data();

    RunningMax running_max{};

    MallocRAII<int> hits_score(n);
    MallocRAII<int> hits_pos(n);

    const auto alignment_capacity = static_cast<int>(1.5 * config.tblen);
    IA maxHit(alignment_capacity);

    /* matrices for alignment scores ending in different states */
    // Since we only need 2 rows we can optimize the memory layout.
    MallocRAII<int_type> dp_rows(6 * (m + 1));
    int_type* const M[2] = {dp_rows.get() + 0 * (m + 1), dp_rows.get() + 1 * (m + 1)};
    int_type* const Ix[2] = {dp_rows.get() + 2 * (m + 1), dp_rows.get() + 3 * (m + 1)};
    int_type* const Iy[2] = {dp_rows.get() + 4 * (m + 1), dp_rows.get() + 5 * (m + 1)};


    /* The sweep runs at whichever width the query allows; the traceback always
       runs at int32, because a window is only as long as the hit it re-aligns and
       most of them are too short for a sixteen wide block. */
    const QueryProfile<int_type> profile(query_sequence, m, dsm, has_positive_gap(dsm));

    M[0][0] = Ix[0][0] = Iy[0][0] = 0;

    /**
     * Row j=0 -- the empty prefix before the first target nt processed.
     * (target_sequence[n-1], the 3' end, since the target is walked backwards).
     *
     * M and Iy need a target nt to exist, so neither state is reachable here.
     *
     * Ix is 0 because an alignemnt must open on a base pair, and never on a bulge.
     *
     * This scheme can't score dangling ends.
     *
     */
    for (auto i = 1u; i <= m; ++i) {
        Iy[0][i] = M[0][i] = neg_inf<int_type>();
        Ix[0][i] = 0;
    }


    const RunningVectorMax running_row_max = score_first_row<int_type>(
        M[1], Ix[1], Iy[1], query_sequence, m, target_sequence[n - 1], dsm);
    running_max.set(running_row_max.score, running_row_max.pos_i, 1);


    /**
     * The first row can only hold 1nt, so NN terms only now begin to apply.
     */

    // run_scores[j - 1] and run_positions[j - 1] hold the best alignment ending at pos j
    int* const run_scores = hits_score.get();
    int* const run_positions = hits_pos.get();

    run_scores[0] = running_row_max.score;
    run_positions[0] = running_row_max.pos_i;


    score_target<int_type>(target_sequence, profile, M, Ix, Iy, run_scores, run_positions, n,
                           threshold, running_max);


    const QueryProfile<std::int32_t> wide_profile(query_sequence, m, dsm, has_positive_gap(dsm));
    HitReporter reporter(query_sequence, target_sequence, n, wide_profile, config, qname, tname);

    reporter.report_sweep(run_scores, run_positions, threshold, running_max);
}
