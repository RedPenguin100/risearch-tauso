#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "dsm.h"
#include "optimization/QueryProfile.h"

/**
 * The int32 profiles the reporting needs, kept between targets.
 *
 * A profile is resolved from the query and the scoring matrix alone -- no term
 * in it reads the target -- so one built for a query stays good for every target
 * that query is ever swept against. Building it per query/target pair was most
 * of the run wherever targets are many and short: the cost is fixed per pair
 * while the sweep it pays for grows with the target.
 *
 * Keyed by the query's record number rather than by a lane, because a run with
 * more than one batch puts a different query in the same lane.
 *
 * BOUNDED. A profile is about 8 * DSM_SIDE^2 * (m + 17) * 4 bytes -- some 43 kB
 * for a 20 nt query -- so a file of many thousands of queries would otherwise
 * grow without limit. Past the budget the profile is built into a scratch slot
 * and thrown away, which is what every caller did before this existed.
 */
class QueryProfileCache {
public:
    /* Enough for a few thousand ASO-length queries; past it, correctness is
       unchanged and only the saving stops. */
    static constexpr std::size_t kBudgetBytes = 256u * 1024u * 1024u;

    explicit QueryProfileCache(const Dsm& dsm) : m_has_positive_gap(has_positive_gap(dsm))
    {
    }

    /* dsm decayed, not by reference: see the note on QueryProfile's constructor. */
    const QueryProfile<std::int32_t>& get(int query_id, const unsigned char* query,
                                          std::uint32_t m, Dsm dsm)
    {
        const auto bytes = profile_bytes(m);

        if (query_id > 0 && m_bytes + bytes <= kBudgetBytes) {
            const auto slot = static_cast<std::size_t>(query_id);
            if (m_cached.size() <= slot) {
                m_cached.resize(slot + 1);
            }
            if (!m_cached[slot]) {
                m_cached[slot] = std::make_unique<QueryProfile<std::int32_t>>(
                    query, m, dsm, m_has_positive_gap);
                m_bytes += bytes;
            }
            return *m_cached[slot];
        }

        m_scratch =
            std::make_unique<QueryProfile<std::int32_t>>(query, m, dsm, m_has_positive_gap);
        return *m_scratch;
    }

private:
    /* What QueryProfile allocates: eight runs over every target context, and two
       that have no target dependence. */
    static std::size_t profile_bytes(std::uint32_t m)
    {
        const std::size_t stride = m + 1 + 16;
        return (8 * DSM_SIDE * DSM_SIDE + 2) * stride * sizeof(std::int32_t);
    }

    bool m_has_positive_gap;
    std::vector<std::unique_ptr<QueryProfile<std::int32_t>>> m_cached;
    std::unique_ptr<QueryProfile<std::int32_t>> m_scratch;
    std::size_t m_bytes = 0;
};
