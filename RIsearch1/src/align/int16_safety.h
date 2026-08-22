#pragma once

#include <climits>
#include <cstdint>

#include "dsm.h"
#include "nucleotide.h" /* GAP */

/**
 * The original software was written with int32 as the scores. However, in practice, it
 * rarely reaches 10k, let alone scores of 30k+. So we can use int16 lanes most of the time.
 * However, for large queries, that would not be true. So we have safety functions that
 * will check the maximal possible score for query length, and will dispatch long queries
 * to wider lanes.
 */


static constexpr short NEG_INF_SHORT = SHRT_MIN;

/**
 * What an unreachable state holds, at whichever width the rows are kept. NEGINF
 * is INT_MIN / 2 and truncates to 0 in a short -- which is a real score meaning
 * no alignment -- so a row written at one width must never take the other's.
 */
template<typename int_type>
constexpr int_type neg_inf();

template<>
constexpr std::int32_t neg_inf<std::int32_t>()
{
    return NEGINF;
}

template<>
constexpr std::int16_t neg_inf<std::int16_t>()
{
    return NEG_INF_SHORT;
}

/**
 * Does extending a target bulge pay? Iy recurses down the target -- Iy[j][i] reads
 * Iy[j-1][i] -- so such a run is as long as the target and the score has no bound
 * at all. Reachable: a negative -d turns the extension into a bonus.
 */
inline bool is_pos_target_bulge(const Dsm& dsm)
{
    for (auto t_prev = 0u; t_prev < DSM_SIDE; t_prev++) {
        for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
            if (dsm[GAP][GAP][t_prev][t_cur] > 0) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Does extending a query bulge pay? Ix recurses along the query -- Ix[j][i] reads
 * Ix[j][i-1] -- so a run is at most m long and the score stays bounded either way.
 * int16_bound still needs this to be false, because it accumulates P as a sum of
 * terms it takes to be at most 0.
 */
inline bool is_pos_query_bulge(const Dsm& dsm)
{
    for (auto q_prev = 0u; q_prev < DSM_SIDE; q_prev++) {
        for (auto q_cur = 0u; q_cur < DSM_SIDE; q_cur++) {
            if (dsm[q_prev][q_cur][GAP][GAP] > 0) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Given dsm matrix and query seq/lenm, what is the possible maximal score?
 * Only meaningful where neither bulge extension pays -- that is what makes those
 * two terms at most 0, which the derivation rests on.
 */
inline long long int16_bound(const Dsm& dsm, const unsigned char* query_sequence, std::uint32_t m)
{
    const short* const flat = &dsm[0][0][0][0];
    long long s = 0; /* largest positive entry              */
    long long l = 0; /* magnitude of largest negative entry */
    for (auto k = 0u; k < DSM_SIDE * DSM_SIDE * DSM_SIDE * DSM_SIDE; k++) {
        const long long v = flat[k];
        if (v > s) {
            s = v;
        }
        if (-v > l) {
            l = -v;
        }
    }

    /* P is the same sum ix_prefix accumulates. */
    long long prefix = 0;
    for (auto i = 2u; i <= m; i++) {
        prefix -= dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][GAP];
    }

    long long bound = 2 * static_cast<long long>(m) * s + s + prefix;
    if (3 * l > bound) {
        bound = 3 * l;
    }
    if (2 * l + s > bound) {
        bound = 2 * l + s;
    }
    return bound;
}

/**
 * Given dsm matrix, query sequence and its length, returns if the maximal score fits in
 * int16 wide lanes. We choose 30k as an upper boundary out of abundance of safety.
 */
inline bool fits_int16(const Dsm& dsm, const unsigned char* query_sequence, std::uint32_t m)
{
    return !is_pos_target_bulge(dsm) && !is_pos_query_bulge(dsm) &&
           int16_bound(dsm, query_sequence, m) <= 30000;
}
