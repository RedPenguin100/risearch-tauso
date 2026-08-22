// Tests for the max helpers in operations.h.

#include <gtest/gtest.h>

#include <limits>

#include "operations.h"

TEST(Operations, Max3PicksTheLargestFromAnyPosition)
{
    EXPECT_EQ(max3(3, 1, 2), 3);
    EXPECT_EQ(max3(1, 3, 2), 3);
    EXPECT_EQ(max3(1, 2, 3), 3);
}

TEST(Operations, Max4PicksTheLargestFromAnyPosition)
{
    EXPECT_EQ(max4(4, 1, 2, 3), 4);
    EXPECT_EQ(max4(1, 4, 2, 3), 4);
    EXPECT_EQ(max4(1, 2, 4, 3), 4);
    EXPECT_EQ(max4(1, 2, 3, 4), 4);
}

TEST(Operations, MaxHandlesNegativesAndTies)
{
    // The recursion compares deeply negative sentinels, so negatives matter.
    EXPECT_EQ(max3(-3, -1, -2), -1);
    EXPECT_EQ(max4(-4, -4, -4, -4), -4);
}

TEST(Operations, MaxSurvivesIntExtremes)
{
    constexpr int kMin = std::numeric_limits<int>::min();
    constexpr int kMax = std::numeric_limits<int>::max();
    EXPECT_EQ(max4(kMin, kMin, kMax, kMin), kMax);
    EXPECT_EQ(max3(kMin, kMin, kMin), kMin);
}

TEST(Operations, FloatVariantsMatchTheIntegerOnes)
{
    EXPECT_FLOAT_EQ(max3f(1.5f, -2.5f, 0.5f), 1.5f);
    EXPECT_FLOAT_EQ(max4f(-1.5f, -2.5f, -0.5f, -3.5f), -0.5f);
}
