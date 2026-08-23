#pragma once

#include <cstdint>
#include <cstdlib>

#include "align/int16_safety.h" /* NEG_INF_SHORT */
#include "dsm.h"
#include "memory/GrowableBuffer.hpp"
#include "nucleotide.h" /* GAP */
#include "optimization/QueryProfile.h"

/**
 * The terms sixteen queries need, laid out so that one load reaches the same
 * term at the same query position in all sixteen.
 *
 * The single query profile stores a term as a run over query positions. Here a
 * term is a run over lanes, sixteen wide, and those runs sit one after another
 * inside a group. A column of the sweep therefore reads one group and has every
 * lane's value for every term it needs.
 *
 * TWO TABLES, BECAUSE ONLY THREE OF THE SEVEN TERMS READ BOTH TARGET
 * NUCLEOTIDES. m_from_m, m_from_iy and iy_from_m are indexed by the dinucleotide
 * (t_prev, t_cur); the other four have t_prev nowhere in them and are indexed by
 * t_cur alone, which is six rows of the table rather than thirty-six.
 *
 * ix_prefix has no target dependence and gets its own run over query positions;
 * iy_extend has no query dependence and is one per dinucleotide.
 */
class BatchedQueryProfile {
public:
    static constexpr unsigned kLanes = 16;
    static constexpr unsigned kContexts = DSM_SIDE * DSM_SIDE;

    /* Term order inside a dinucleotide group. */
    enum : unsigned {
        kMFromM = 0, /* dsm[q_prev][q_cur][t_prev][t_cur] -- extend a pair        */
        kMFromIy,    /* dsm[GAP][q_cur][t_prev][t_cur]    -- close a target bulge */
        kIyFromM,    /* dsm[q_cur][GAP][t_prev][t_cur]    -- open a target bulge  */
        kPairTerms
    };

    /* And inside a t_cur group. */
    enum : unsigned {
        kMFromIx = 0, /* dsm[q_prev][q_cur][GAP][t_cur] plus ix_prefix[i - 1],
                         which is what the Ix it is added to is carried without  */
        kMOpen,       /* dsm[GAP][q_cur][GAP][t_cur]       -- open on this pair   */
        kClose,       /* dsm[q_cur][GAP][t_cur][GAP]       -- terminate after it  */
        kIxFromMScan, /* ix_from_m with ix_prefix taken out                       */
        kSoloTerms
    };

    static constexpr unsigned kPairGroup = kPairTerms * kLanes;
    static constexpr unsigned kSoloGroup = kSoloTerms * kLanes;

    /* What a lane carrying no query, or a column past the end of a shorter one,
       holds in every term. Each of M's arms then saturates to the magic value
       and the row max candidate with it, so a padded column can never win a row
       max nor an argmax. Queries of different lengths cost nothing in the sweep
       because of this: it reads every lane to the batch's longest m and the
       short ones lose on their own. */
    static constexpr std::int16_t kDead = NEG_INF_SHORT;

    /* What a term is allowed to reach, which is the cap fits_int16 puts on every
       value the sweep computes. The folded terms are sums of two quantities that
       are each inside it, so their own width is checked rather than argued. */
    static constexpr int kTermLimit = 30000;

    /**
     * queries[l] / lengths[l] for l < count; the remaining lanes are dead. Every
     * column of the table is filled, so a lane's own length is the only thing
     * that distinguishes it from the batch's.
     *
     * False when a folded term does not fit a short, which leaves the batch with
     * no sweep and its queries to be aligned one at a time.
     */
    bool build(const unsigned char* const* queries, const std::uint32_t* lengths, unsigned count,
               Dsm& dsm, std::uint32_t m)
    {
        m_m = m;
        m_stride = m + 1;

        m_pair.reserve(kContexts * m_stride * kPairGroup);
        m_solo.reserve(DSM_SIDE * m_stride * kSoloGroup);
        m_ix_prefix.reserve(m_stride * kLanes);

        std::int16_t* const ix_prefix = m_ix_prefix.get();

        /* The running total of a query bulge over a gap on both sides, which the
           Ix scan takes out of its candidates. No target dependence. */
        for (auto lane = 0u; lane < kLanes; lane++) {
            ix_prefix[0 * kLanes + lane] = 0;
            ix_prefix[1 * kLanes + lane] = 0;
        }
        for (auto i = 2u; i <= m; i++) {
            for (auto lane = 0u; lane < kLanes; lane++) {
                std::int16_t prefix = ix_prefix[(i - 1) * kLanes + lane];
                if (lane < count && i <= lengths[lane]) {
                    prefix = static_cast<std::int16_t>(
                        prefix + dsm[queries[lane][i - 2]][queries[lane][i - 1]][GAP][GAP]);
                }
                ix_prefix[i * kLanes + lane] = prefix;
            }
        }

        for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
            std::int16_t* const base = m_solo.get() + t_cur * m_stride * kSoloGroup;
            m_solo_base[t_cur] = base;
            /* Column 0 exists only so that a column index can be used directly;
               nothing reads it. */
            for (auto k = 0u; k < kSoloGroup; k++) {
                base[k] = kDead;
            }
            for (auto lane = 0u; lane < kLanes; lane++) {
                m_close_max[t_cur * kLanes + lane] = kDead;
            }
            for (auto i = 1u; i <= m; i++) {
                std::int16_t* const group = base + i * kSoloGroup;
                for (auto lane = 0u; lane < kLanes; lane++) {
                    if (lane >= count || i > lengths[lane]) {
                        for (auto term = 0u; term < kSoloTerms; term++) {
                            group[term * kLanes + lane] = kDead;
                        }
                        continue;
                    }
                    const unsigned char* const q = queries[lane];
                    const auto q_cur = q[i - 1];
                    /* Column 1 has no predecessor, so the q_prev terms are never
                       read there; GAP is a placeholder. */
                    const auto q_prev = i >= 2 ? q[i - 2] : static_cast<unsigned char>(GAP);

                    const int from_ix =
                        dsm[q_prev][q_cur][GAP][t_cur] + ix_prefix[(i - 1) * kLanes + lane];
                    if (from_ix < -kTermLimit || from_ix > kTermLimit) {
                        return false;
                    }
                    group[kMFromIx * kLanes + lane] = static_cast<std::int16_t>(from_ix);
                    group[kMOpen * kLanes + lane] = dsm[GAP][q_cur][GAP][t_cur];
                    const std::int16_t close = dsm[q_cur][GAP][t_cur][GAP];
                    group[kClose * kLanes + lane] = close;
                    if (close > m_close_max[t_cur * kLanes + lane]) {
                        m_close_max[t_cur * kLanes + lane] = close;
                    }
                    group[kIxFromMScan * kLanes + lane] = static_cast<std::int16_t>(
                        dsm[q_prev][q_cur][t_cur][GAP] - ix_prefix[i * kLanes + lane]);
                }
            }
        }

        for (auto t_prev = 0u; t_prev < DSM_SIDE; t_prev++) {
            for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
                const auto ctx = QueryProfile<std::int32_t>::context(t_prev, t_cur);
                m_iy_extend[ctx] = dsm[GAP][GAP][t_prev][t_cur];

                std::int16_t* const base = m_pair.get() + ctx * m_stride * kPairGroup;
                m_pair_base[ctx] = base;
                for (auto k = 0u; k < kPairGroup; k++) {
                    base[k] = kDead;
                }

                for (auto i = 1u; i <= m; i++) {
                    std::int16_t* const group = base + i * kPairGroup;
                    for (auto lane = 0u; lane < kLanes; lane++) {
                        if (lane >= count || i > lengths[lane]) {
                            for (auto term = 0u; term < kPairTerms; term++) {
                                group[term * kLanes + lane] = kDead;
                            }
                            continue;
                        }
                        const unsigned char* const q = queries[lane];
                        const auto q_cur = q[i - 1];
                        const auto q_prev = i >= 2 ? q[i - 2] : static_cast<unsigned char>(GAP);

                        group[kMFromM * kLanes + lane] = dsm[q_prev][q_cur][t_prev][t_cur];
                        group[kMFromIy * kLanes + lane] = dsm[GAP][q_cur][t_prev][t_cur];
                        group[kIyFromM * kLanes + lane] = dsm[q_cur][GAP][t_prev][t_cur];
                    }
                }
            }
        }
        return true;
    }

    /* One column's terms, sixteen lanes of each, named as the single query
       profile names them so a sweep written against either reads the same. */
    struct ColumnTerms {
        const std::int16_t* m_from_m;
        const std::int16_t* m_from_ix;
        const std::int16_t* m_from_iy;
        const std::int16_t* m_open;
        const std::int16_t* close;
        const std::int16_t* iy_from_m;
    };

    /* The two tables a row reads, and the term that has no query dependence. */
    struct RowView {
        const std::int16_t* pair;
        const std::int16_t* solo;
        std::int16_t iy_extend;

        ColumnTerms column(unsigned i) const
        {
            const std::int16_t* const p = pair + i * kPairGroup;
            const std::int16_t* const s = solo + i * kSoloGroup;
            return {p + kMFromM * kLanes, s + kMFromIx * kLanes, p + kMFromIy * kLanes,
                    s + kMOpen * kLanes,  s + kClose * kLanes,   p + kIyFromM * kLanes};
        }
    };

    RowView row(unsigned ctx, unsigned t_cur) const
    {
        return {m_pair_base[ctx], m_solo_base[t_cur], m_iy_extend[ctx]};
    }

    /* The largest close a query can be terminated with, over its whole run of
       columns. A row's scores are its columns' M plus their close, so this bounds
       a row max from the row's largest M alone. */
    const std::int16_t* close_max(unsigned t_cur) const
    {
        return m_close_max + t_cur * kLanes;
    }

    /* ix_from_m_scan belongs to the row the bulge is in, which is the row above
       the one running the scan, so it is read at that row's target nucleotide. */
    const std::int16_t* scan_terms(unsigned t_cur) const
    {
        return m_solo_base[t_cur] + kIxFromMScan * kLanes;
    }

    std::uint32_t m() const
    {
        return m_m;
    }

    /* First group of a dinucleotide, and of a t_cur: the sweep adds i times the
       matching group size to reach column i. */
    const std::int16_t* pair_base(unsigned ctx) const
    {
        return m_pair_base[ctx];
    }

    const std::int16_t* solo_base(unsigned t_cur) const
    {
        return m_solo_base[t_cur];
    }

    const std::int16_t* ix_prefix() const
    {
        return m_ix_prefix.get();
    }

    std::int16_t iy_extend(unsigned ctx) const
    {
        return m_iy_extend[ctx];
    }

private:
    std::uint32_t m_m = 0;
    std::uint32_t m_stride = 0;

    GrowableBuffer<std::int16_t> m_pair;
    GrowableBuffer<std::int16_t> m_solo;
    GrowableBuffer<std::int16_t> m_ix_prefix;

    const std::int16_t* m_pair_base[kContexts]{};
    const std::int16_t* m_solo_base[DSM_SIDE]{};
    std::int16_t m_iy_extend[kContexts]{};
    std::int16_t m_close_max[DSM_SIDE * kLanes]{};
};
