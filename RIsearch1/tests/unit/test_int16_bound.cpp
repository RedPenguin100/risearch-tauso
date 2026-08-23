// Tests for the bound the int16 kernel needs. A wrong bound is a wrong answer
// and not a crash, so these pin both the arithmetic and the rejections.

#include <gtest/gtest.h>

#include <vector>

#include "align/avx2/primitives.h"
#include "align/int16_safety.h"
#include "dsm.h"
#include "nucleotide.h"

namespace {

std::vector<unsigned char> alternating(unsigned char a, unsigned char b, std::uint32_t m)
{
    std::vector<unsigned char> q(m);
    for (auto i = 0u; i < m; i++) {
        q[i] = (i % 2) ? b : a;
    }
    return q;
}

void matrix(const char* name, Dsm& out)
{
    getMat(name, &out[0][0][0][0], 0, 0);
}

} // namespace

TEST(Int16Bound, HoldsForTheLengthsAsosAreWritten)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);

    for (auto m = 15u; m <= 25u; m++) {
        const auto q = alternating(0, 1, m);
        EXPECT_TRUE(fits_int16(dsm, q.data(), m)) << "m = " << m;
    }
}

TEST(Int16Bound, GrowsWithTheQueryUntilItStopsFitting)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);

    const auto short_q = alternating(0, 1, 20);
    const auto long_q = alternating(0, 1, 60);

    EXPECT_LT(int16_bound(dsm, short_q.data(), 20), int16_bound(dsm, long_q.data(), 60));
    EXPECT_TRUE(fits_int16(dsm, short_q.data(), 20));
    EXPECT_FALSE(fits_int16(dsm, long_q.data(), 60));
}

TEST(Int16Bound, LeavesRoomBelowWhatAShortHolds)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);

    const auto q = alternating(0, 1, 20);
    const auto bound = int16_bound(dsm, q.data(), 20);

    EXPECT_GT(bound, 0);
    EXPECT_LE(bound, 30000);
    // The cap is what keeps a saturating add off every real value.
    EXPECT_LT(bound, 32767);
}

TEST(BulgeExtension, HoldsForEveryMatrixRisearchShips)
{
    for (const char* name : {"su95", "su95_noGU", "t04", "t99", "slh04_noGU"}) {
        Dsm dsm;
        matrix(name, dsm);
        const auto q = alternating(0, 1, 20);
        EXPECT_TRUE(!is_pos_target_bulge(dsm) && !is_pos_query_bulge(dsm)) << name;
    }
}

TEST(BulgeExtension, FailsWhenExtendingATargetBulgePays)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);
    const auto q = alternating(0, 1, 20);
    ASSERT_TRUE(!is_pos_target_bulge(dsm) && !is_pos_query_bulge(dsm));

    dsm[GAP][GAP][0][1] = 1;

    EXPECT_FALSE(!is_pos_target_bulge(dsm) && !is_pos_query_bulge(dsm));
    EXPECT_FALSE(fits_int16(dsm, q.data(), 20));
}

TEST(BulgeExtension, FailsWhenExtendingTheQuerysOwnBulgePays)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);
    const auto q = alternating(0, 1, 20);
    ASSERT_TRUE(!is_pos_target_bulge(dsm) && !is_pos_query_bulge(dsm));

    dsm[q[0]][q[1]][GAP][GAP] = 1;

    EXPECT_FALSE(!is_pos_target_bulge(dsm) && !is_pos_query_bulge(dsm));
    EXPECT_FALSE(fits_int16(dsm, q.data(), 20));
}

TEST(BulgeExtension, IsAPropertyOfTheMatrixAndNoQuery)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);
    const auto q = alternating(0, 1, 20);

    /* A dinucleotide this query cannot produce still rejects the matrix, which is
       what lets the check run once before any sequence is read. */
    dsm[1][1][GAP][GAP] = 1;

    EXPECT_TRUE(is_pos_query_bulge(dsm));
    EXPECT_FALSE(fits_int16(dsm, q.data(), 20));
}

TEST(Int16Bound, HoldsForEveryMatrixRisearchShips)
{
    for (const char* name : {"su95", "su95_noGU", "t04", "t99", "slh04_noGU"}) {
        Dsm dsm;
        matrix(name, dsm);
        const auto q = alternating(0, 1, 20);
        EXPECT_TRUE(fits_int16(dsm, q.data(), 20)) << name;
    }
}

// The conditions dispatch.h uses to pick the narrow sweep. Nothing else fails if
// these stop holding -- the run stays correct and quietly loses the width.

TEST(Int16Dispatch, AnAsoOfTheUsualLengthTakesTheNarrowSweep)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);
    const auto q = alternating(0, 1, 20);

    EXPECT_GT(20u, v_lanes<std::int16_t>()) << "needs one full block past column 1";
    EXPECT_TRUE(fits_int16(dsm, q.data(), 20));
}

TEST(Int16Dispatch, TheTwoWidthsAreWhatTheKernelAssumes)
{
    // The block loops and their clamps are written in terms of these.
    EXPECT_EQ(v_lanes<std::int32_t>(), 8u);
    EXPECT_EQ(v_lanes<std::int16_t>(), 16u);
}

TEST(Int16Dispatch, ALongQueryFallsBackRatherThanOverflowing)
{
    Dsm dsm;
    matrix("su95_noGU", dsm);
    const auto q = alternating(0, 1, 60);

    EXPECT_GT(int16_bound(dsm, q.data(), 60), 30000);
    EXPECT_FALSE(fits_int16(dsm, q.data(), 60));
}
