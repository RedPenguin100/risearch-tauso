// The lane-interleaved profile against the single-query one it must agree with.
// Every term, every context, every query position, lane for lane -- plus what a
// lane holds where it has no query, or where its query has already ended.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

#include "align/optimization/BatchedQueryProfile.h"
#include "align/optimization/QueryProfile.h"
#include "dsm.h"
#include "nucleotide.h"

namespace {

using Profile = QueryProfile<std::int32_t>;
constexpr unsigned kLanes = BatchedQueryProfile::kLanes;

std::vector<unsigned char> sequence(std::uint32_t m, unsigned seed)
{
    std::vector<unsigned char> q(m);
    for (auto i = 0u; i < m; i++) {
        q[i] = static_cast<unsigned char>((seed * 7 + i * 3 + i / 4) % 4);
    }
    return q;
}

void matrix(Dsm& out)
{
    getMat("su95_noGU", &out[0][0][0][0], 0, 0);
}

const std::int16_t* solo(const BatchedQueryProfile& b, unsigned t_cur, unsigned i, unsigned term)
{
    return b.solo_base(t_cur) + i * BatchedQueryProfile::kSoloGroup + term * kLanes;
}

const std::int16_t* pair(const BatchedQueryProfile& b, unsigned ctx, unsigned i, unsigned term)
{
    return b.pair_base(ctx) + i * BatchedQueryProfile::kPairGroup + term * kLanes;
}

} // namespace

TEST(BatchedQueryProfile, HoldsWhatSixteenSingleProfilesHold)
{
    Dsm dsm;
    matrix(dsm);

    std::vector<std::vector<unsigned char>> qs;
    std::vector<const unsigned char*> ptrs;
    std::vector<std::uint32_t> lens;
    for (auto lane = 0u; lane < kLanes; lane++) {
        qs.push_back(sequence(20, lane));
    }
    for (auto lane = 0u; lane < kLanes; lane++) {
        ptrs.push_back(qs[lane].data());
        lens.push_back(20);
    }

    BatchedQueryProfile batch;
    ASSERT_TRUE(batch.build(ptrs.data(), lens.data(), kLanes, dsm, 20));

    for (auto lane = 0u; lane < kLanes; lane++) {
        const Profile one(qs[lane].data(), 20, dsm, has_positive_gap(dsm));
        for (auto t_prev = 0u; t_prev < DSM_SIDE; t_prev++) {
            for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
                const auto ctx = Profile::context(t_prev, t_cur);
                const auto row = one.row(ctx);
                EXPECT_EQ(batch.iy_extend(ctx), row.iy_extend) << "ctx " << ctx;
                for (auto i = 1u; i <= 20; i++) {
                    EXPECT_EQ(pair(batch, ctx, i, BatchedQueryProfile::kMFromM)[lane],
                              row.m_from_m[i])
                        << "lane " << lane << " ctx " << ctx << " i " << i;
                    EXPECT_EQ(pair(batch, ctx, i, BatchedQueryProfile::kMFromIy)[lane],
                              row.m_from_iy[i]);
                    EXPECT_EQ(pair(batch, ctx, i, BatchedQueryProfile::kIyFromM)[lane],
                              row.iy_from_m[i]);
                }
            }
        }
        for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
            const auto row = one.row(Profile::context(0, t_cur));
            for (auto i = 1u; i <= 20; i++) {
                EXPECT_EQ(solo(batch, t_cur, i, BatchedQueryProfile::kMOpen)[lane], row.m_open[i]);
                EXPECT_EQ(solo(batch, t_cur, i, BatchedQueryProfile::kClose)[lane], row.close[i]);
            }
        }
        /* The running total of ix_extend, which is what the Ix scan takes out of
           its candidates and the only form of it either profile keeps. From 2: a
           query bulge needs two positions, so neither writes the first two. */
        const auto ix_row = one.row(Profile::context(0, 0));
        for (auto i = 2u; i <= 20; i++) {
            EXPECT_EQ(batch.ix_prefix()[i * kLanes + lane], ix_row.ix_prefix[i])
                << "lane " << lane << " i " << i;
        }
    }
}

TEST(BatchedQueryProfile, ALaneWithNoQueryIsDeadInEveryTerm)
{
    Dsm dsm;
    matrix(dsm);

    const auto q = sequence(20, 1);
    const unsigned char* ptrs[kLanes] = {};
    std::uint32_t lens[kLanes] = {};
    for (auto lane = 0u; lane < kLanes; lane++) {
        ptrs[lane] = q.data();
        lens[lane] = 20;
    }

    BatchedQueryProfile batch;
    ASSERT_TRUE(batch.build(ptrs, lens, 5, dsm, 20)); // only five lanes carry a query

    for (auto lane = 5u; lane < kLanes; lane++) {
        for (auto i = 1u; i <= 20; i++) {
            EXPECT_EQ(pair(batch, 0, i, BatchedQueryProfile::kMFromM)[lane],
                      BatchedQueryProfile::kDead)
                << "lane " << lane << " i " << i;
            EXPECT_EQ(solo(batch, 0, i, BatchedQueryProfile::kMOpen)[lane],
                      BatchedQueryProfile::kDead);
        }
    }
}

TEST(BatchedQueryProfile, AShortQueryIsDeadPastItsOwnEnd)
{
    Dsm dsm;
    matrix(dsm);

    // Lane 0 runs to 20, lane 1 stops at 12; the batch is built to 20.
    const auto q_long = sequence(20, 3);
    const auto q_short = sequence(12, 4);
    const unsigned char* ptrs[2] = {q_long.data(), q_short.data()};
    const std::uint32_t lens[2] = {20, 12};

    BatchedQueryProfile batch;
    ASSERT_TRUE(batch.build(ptrs, lens, 2, dsm, 20));

    const Profile one(q_short.data(), 12, dsm, has_positive_gap(dsm));
    const auto row = one.row(Profile::context(0, 1));

    for (auto i = 1u; i <= 12; i++) {
        EXPECT_EQ(pair(batch, Profile::context(0, 1), i, BatchedQueryProfile::kMFromM)[1],
                  row.m_from_m[i])
            << "i " << i;
    }
    for (auto i = 13u; i <= 20; i++) {
        EXPECT_EQ(pair(batch, Profile::context(0, 1), i, BatchedQueryProfile::kMFromM)[1],
                  BatchedQueryProfile::kDead)
            << "i " << i;
        EXPECT_EQ(solo(batch, 1, i, BatchedQueryProfile::kClose)[1], BatchedQueryProfile::kDead);
    }
    // The longer lane is untouched by its neighbour ending early.
    const Profile full(q_long.data(), 20, dsm, has_positive_gap(dsm));
    const auto full_row = full.row(Profile::context(0, 1));
    for (auto i = 1u; i <= 20; i++) {
        EXPECT_EQ(pair(batch, Profile::context(0, 1), i, BatchedQueryProfile::kMFromM)[0],
                  full_row.m_from_m[i])
            << "i " << i;
    }
}
