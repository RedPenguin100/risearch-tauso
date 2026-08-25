#pragma once

#include <cstdint>

#include "align/force_start.h"
#include "align/int16_safety.h"
#include "align/linspace.h"
#include "cli/cli.h"
#include "memory/ByteBuffer.hpp"

static void run_alignment(const ByteBuffer& query_seq, const ByteBuffer& target_seq, Dsm& dsm,
                          const char* nameQ, const char* nameT, const config_st& config)
{
    if (uses_force_start(config)) {
        RIs_force_start_end_init(config.force_start_val, config.pos_weights, query_seq, target_seq,
                                 dsm, config.mat_name);
    } else {
        /* Sixteen columns to a register rather than eight, where the scores this
           query and this matrix can reach leave room in a short and there is at
           least one full block of them. Without AVX2 there are no blocks at all
           and the scalar sweep is the whole implementation: it gains nothing from
           a narrower score and pays to widen every one it reads. */
        const auto m = static_cast<std::uint32_t>(query_seq.size());
        if (CPU_HAS_AVX2 && m > v_lanes<std::int16_t>() &&
            fits_int16(dsm, query_seq.unsigned_data(), m)) {
            RIs_linSpace<std::int16_t>(query_seq, target_seq, dsm, config.min_score, nameQ, nameT,
                                       config);
        } else {
            RIs_linSpace<std::int32_t>(query_seq, target_seq, dsm, config.min_score, nameQ, nameT,
                                       config);
        }
    }
}
