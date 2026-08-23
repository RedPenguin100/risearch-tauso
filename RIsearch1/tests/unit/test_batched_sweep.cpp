// The batched sweep against the sweep it has to agree with, query by query and
// target position by target position.
//
// End-to-end byte-identity says only that something differs. This says which
// query and which row, which is the difference between reading a diff of hits
// and reading a diff of the recurrence.

#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

#include "align/QueryBatch.h"
#include "align/ScoreTarget.h"
#include "align/optimization/QueryProfile.h"
#include "dsm.h"
#include "memory/ByteBuffer.hpp"
#include "nucleotide.h"

#if RISEARCH1_HAS_AVX2

namespace {

using Profile = QueryProfile<std::int32_t>;

std::vector<unsigned char> sequence(std::uint32_t len, unsigned seed)
{
    std::vector<unsigned char> s(len);
    auto x = seed * 2654435761u + 1u;
    for (auto i = 0u; i < len; i++) {
        x = x * 1103515245u + 12345u;
        s[i] = static_cast<unsigned char>((x >> 16) & 3u);
    }
    return s;
}

void to_buffer(const std::vector<unsigned char>& s, ByteBuffer& out)
{
    out.clear();
    out.append(reinterpret_cast<const char*>(s.data()), s.size());
}

/* The ordinary sweep for one query, set up as linspace sets it up. */
struct SingleSweep {
    std::vector<int> hs, hp;

    SingleSweep(const std::vector<unsigned char>& query, const std::vector<unsigned char>& target,
                Dsm& dsm, int threshold)
        : hs(target.size()), hp(target.size())
    {
        const auto m = static_cast<std::uint32_t>(query.size());
        const auto n = static_cast<int>(target.size());
        const auto* q = query.data();
        const auto* t = target.data();

        const Profile profile(q, m, dsm, has_positive_gap(dsm));
        std::vector<int> rows(6 * (m + 1));
        int* const M[2] = {rows.data() + 0 * (m + 1), rows.data() + 1 * (m + 1)};
        int* const Ix[2] = {rows.data() + 2 * (m + 1), rows.data() + 3 * (m + 1)};
        int* const Iy[2] = {rows.data() + 4 * (m + 1), rows.data() + 5 * (m + 1)};

        M[0][0] = Ix[0][0] = Iy[0][0] = 0;
        for (auto i = 1u; i <= m; ++i) {
            Iy[0][i] = M[0][i] = NEGINF;
            Ix[0][i] = 0;
        }
        Iy[1][0] = 0;
        Ix[1][0] = M[1][0] = NEGINF;

        const auto t_last = t[n - 1];
        M[1][1] = dsm[GAP][q[0]][GAP][t_last];
        RunningVectorMax first_row{};
        first_row.set(M[1][1] + dsm[q[0]][GAP][t_last][GAP], 1);
        Ix[1][1] = Iy[1][1] = NEGINF;

        for (auto i = 2u; i <= m; i++) {
            const auto q_prev = q[i - 2];
            const auto q_cur = q[i - 1];
            M[1][i] = dsm[GAP][q_cur][GAP][t_last];
            first_row.set_if_better(M[1][i] + dsm[q_cur][GAP][t_last][GAP], static_cast<int>(i));
            Ix[1][i] =
                max3(0, M[1][i - 1] != 0 ? M[1][i - 1] + dsm[q_prev][q_cur][t_last][GAP] : -1,
                     Ix[1][i - 1] != 0 ? Ix[1][i - 1] + dsm[q_prev][q_cur][GAP][GAP] : -1);
            Iy[1][i] = NEGINF;
        }

        RunningMax running_max{};
        running_max.set(first_row.score, first_row.pos_i, 1);
        hs[0] = first_row.score;
        hp[0] = first_row.pos_i;

        score_target<std::int32_t>(t, profile, M, Ix, Iy, hs.data(), hp.data(), n, threshold,
                                   running_max);
    }
};

} // namespace

TEST(BatchedSweep, EveryQueryMatchesTheSweepItReplaces)
{
    Dsm dsm;
    getMat("su95_noGU", &dsm[0][0][0][0], 0, 0);

    // Sixteen queries, deliberately of different lengths so the dead columns of
    // the shorter ones are exercised beside the longest.
    const std::uint32_t lengths[16] = {20, 20, 12, 31, 17, 20, 44, 20,
                                       9,  25, 20, 20, 33, 20, 20, 18};
    std::vector<std::vector<unsigned char>> queries;
    for (auto k = 0u; k < 16; k++) {
        queries.push_back(sequence(lengths[k], k + 1));
    }
    const auto target = sequence(4000, 99);

    ByteBuffer target_buf, query_buf;
    to_buffer(target, target_buf);

    QueryBatch batch;
    for (auto k = 0u; k < 16; k++) {
        to_buffer(queries[k], query_buf);
        batch.add(query_buf, "q", static_cast<int>(k), lengths[k]);
    }

    constexpr int threshold = 900;
    ASSERT_TRUE(batch.sweep(target_buf, dsm, threshold)) << "the batch refused to sweep";
    ASSERT_EQ(batch.swept_length(), static_cast<int>(target.size()));

    for (auto k = 0u; k < 16; k++) {
        const SingleSweep one(queries[k], target, dsm, threshold);
        const auto* got_hs = batch.scores(k);
        const auto* got_hp = batch.positions(k);

        for (auto j = 0; j < batch.swept_length(); j++) {
            ASSERT_EQ(got_hs[j], one.hs[j]) << "score differs: query " << k << " (m=" << lengths[k]
                                            << ") at target position " << j + 1;
            /* A position is only read where its score is worth reporting. */
            if (one.hs[j] > threshold) {
                ASSERT_EQ(got_hp[j], one.hp[j])
                    << "position differs: query " << k << " at target position " << j + 1;
            }
        }
    }
}

TEST(BatchedSweep, RefusesABatchTooSmallToPayForItsLanes)
{
    Dsm dsm;
    getMat("su95_noGU", &dsm[0][0][0][0], 0, 0);

    const auto q = sequence(20, 5);
    const auto target = sequence(2000, 11);
    ByteBuffer target_buf, query_buf;
    to_buffer(target, target_buf);
    to_buffer(q, query_buf);

    QueryBatch batch;
    for (auto k = 0u; k < 4; k++) {
        batch.add(query_buf, "q", static_cast<int>(k), 20);
    }
    EXPECT_FALSE(batch.sweep(target_buf, dsm, 900));
}

#endif
