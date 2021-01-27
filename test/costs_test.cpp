
#include "../costs.hpp"
#include "CppUTest/TestHarness.h"

using namespace marzone;

TEST_GROUP(CostsTestsGroup)
{
};

TEST(CostsTestsGroup, defaultCost_test)
{
    sfname fnames = {}; // default, no cost file
    Costs c(fnames);

    CHECK_EQUAL(1, c.costCount);
    CHECK_EQUAL(0, c.GetCostIndex("cost")); // default cost name is "cost"
    CHECK(c.Contains("cost"));
    CHECK(!c.Contains("cost1"));
}

TEST(CostsTestsGroup, customCost_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);

    CHECK_EQUAL(3, c.costCount);
    CHECK(c.Contains("area"));
    CHECK(c.Contains("salmonfishing"));
    CHECK(c.Contains("squidfishing"));
}