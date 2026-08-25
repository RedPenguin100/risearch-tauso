#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "dsm.h"
#include "optimization/BatchedQueryProfile.h"

/**
 * The batched profiles, kept between targets.
 *
 * BatchedQueryProfile::build resolves every term of sixteen queries against
 * every target context. Nothing in it reads the target, so a profile built for
 * one batch stays good for every target that batch is swept against -- and it is
 * about 23% of the work of the sweep it sets up, so rebuilding it per target was
 * most of what remained around the kernel.
 *
 * Keyed by the record number of the batch's first query, which is what
 * identifies a batch: a run with more than sixteen queries cycles through
 * several batches per target, and they must not share an entry.
 *
 * Whether a batch can be swept at all is decided from the queries and the matrix
 * too, so the answer is cached beside the profile rather than re-derived per
 * target.
 *
 * BOUNDED, like the single-query cache, though a batch profile is about 5.5 kB
 * per query against that one's 42 kB. Past the budget a batch builds into a
 * scratch slot, as it did before this existed.
 */
class BatchedProfileCache {
public:
    static constexpr std::size_t kBudgetBytes = 256u * 1024u * 1024u;

    /* What a batch resolved to. profile is null where the batch cannot be swept,
       which is a decision this remembers rather than repeats. */
    struct Entry {
        BatchedQueryProfile* profile = nullptr;
        std::uint32_t m = 0; /* the batch's longest query */
        bool decided = false;
    };

    Entry& entry(int first_query_id)
    {
        if (first_query_id > 0 && static_cast<std::size_t>(first_query_id) < kMaxKey) {
            const auto slot = static_cast<std::size_t>(first_query_id);
            if (m_slots.size() <= slot) {
                m_slots.resize(slot + 1);
            }
            return m_slots[slot];
        }
        m_scratch = Entry{};
        return m_scratch;
    }

    /* Storage for a profile this cache will keep, or null when the budget is
       spent and the caller should use scratch(). */
    BatchedQueryProfile* make_kept(std::uint32_t m)
    {
        if (m_bytes + profile_bytes(m) > kBudgetBytes) {
            return nullptr;
        }
        m_kept.push_back(std::make_unique<BatchedQueryProfile>());
        m_bytes += profile_bytes(m);
        return m_kept.back().get();
    }

    /* The profile a batch past the budget builds into, reused each time. */
    BatchedQueryProfile& scratch()
    {
        if (!m_scratch_profile) {
            m_scratch_profile = std::make_unique<BatchedQueryProfile>();
        }
        return *m_scratch_profile;
    }

private:
    /* A record number past this is not worth a slot in the index. */
    static constexpr std::size_t kMaxKey = 1u << 22;

    static std::size_t profile_bytes(std::uint32_t m)
    {
        const std::size_t stride = m + 1;
        const std::size_t pair = DSM_SIDE * DSM_SIDE * stride * BatchedQueryProfile::kPairGroup;
        const std::size_t solo = DSM_SIDE * stride * BatchedQueryProfile::kSoloGroup;
        return (pair + solo + stride * BatchedQueryProfile::kLanes) * sizeof(std::int16_t);
    }

    std::vector<Entry> m_slots;
    std::vector<std::unique_ptr<BatchedQueryProfile>> m_kept;
    std::unique_ptr<BatchedQueryProfile> m_scratch_profile;
    Entry m_scratch;
    std::size_t m_bytes = 0;
};
