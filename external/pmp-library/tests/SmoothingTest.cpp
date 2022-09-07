// Copyright 2017-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "gtest/gtest.h"

#include "pmp/algorithms/Smoothing.h"
#include "pmp/algorithms/DifferentialGeometry.h"
#include "Helpers.h"

using namespace pmp;

TEST(SmoothingTest, implicit_smoothing)
{
    auto mesh = open_cone();
    auto area_before = surface_area(mesh);
    Smoothing ss(mesh);
    ss.implicit_smoothing(0.01, false, false);
    auto area_after = surface_area(mesh);
    EXPECT_LT(area_after, area_before);
}

TEST(SmoothingTest, implicit_smoothing_with_rescale)
{
    auto mesh = open_cone();
    auto area_before = surface_area(mesh);
    Smoothing(mesh).implicit_smoothing(0.01, false, true);
    auto area_after = surface_area(mesh);
    EXPECT_FLOAT_EQ(area_after, area_before);
}

TEST(SmoothingTest, explicit_smoothing)
{
    auto mesh = open_cone();
    auto area_before = surface_area(mesh);
    Smoothing ss(mesh);
    ss.explicit_smoothing(10, true);
    ss.explicit_smoothing(10, false);
    auto area_after = surface_area(mesh);
    EXPECT_LT(area_after, area_before);
}
