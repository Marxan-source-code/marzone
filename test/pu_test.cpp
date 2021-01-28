#include "../pu.hpp"
#include "CppUTest/TestHarness.h"

using namespace marzone;

TEST_GROUP(PuTestsGroup)
{
};

// test ensures pu file gets correctly parsed.
// no pulock or puzone supplied
TEST(PuTestsGroup, ReadPuData_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);
    
    Pu pu(fnames, c, 0);

    // Check costs are set correctly. Note that salmonfishing is missing so should be set to 1 by default.
    CHECK_EQUAL(5, pu.puno);
    CHECK(pu.LookupIndex(2) != -1);
    CHECK(pu.LookupIndex(5) != -1);
    CHECK(pu.LookupIndex(10) != -1);
    CHECK(pu.LookupIndex(11) != -1);
    CHECK(pu.LookupIndex(8) != -1);
    CHECK_EQUAL(-1, pu.LookupIndex(1)); // non existent id.

    // Costs sum should just be summation of all the costs found in the file.
    CHECK_EQUAL(7, pu.puList[0].cost);
    CHECK_EQUAL(67, pu.puList[1].cost);
    CHECK_EQUAL(667, pu.puList[2].cost);
    CHECK_EQUAL(3, pu.puList[3].cost);
    CHECK_EQUAL(3, pu.puList[4].cost);

    // Cost breakdown should be c1, 0, c2. since one of the costs is missing.
    vector<double> c1 {5, 1, 1};
    vector<double> c2 {55, 1, 11};
    vector<double> c3 {555, 1, 111};
    vector<double> c4 {1, 1, 1};
    CHECK(c1 == pu.puList[0].costBreakdown);
    CHECK(c2 == pu.puList[1].costBreakdown);
    CHECK(c3 == pu.puList[2].costBreakdown);
    CHECK(c4 == pu.puList[3].costBreakdown);
    CHECK(c4 == pu.puList[4].costBreakdown);
}

// Test pulock and puzone
TEST(PuTestsGroup, PuLockPuZone_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.pulockname = "data/pulock_test1.dat";
    fnames.puzonename = "data/puzone_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);
    
    Pu pu(fnames, c, 0);
    
    // Ensure locked pus get returned as such.
    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);

    CHECK_EQUAL(2, pu.puLockCount);
    CHECK_EQUAL(1, pu.GetPuLock(ind1));
    CHECK_EQUAL(2, pu.GetPuLock(ind2));
    CHECK_EQUAL(-1, pu.GetPuLock(ind3));

    // ensure locked indices are returned
    vector<int> lockedInd = pu.GetPuLockedIndices();
    CHECK_EQUAL(2, lockedInd.size());
    CHECK(lockedInd[0] == ind1 || lockedInd[0] == ind2);
    CHECK(lockedInd[1] == ind1 || lockedInd[1] == ind2);

    // check puzone correct entered for each pu.
    CHECK_EQUAL(1, pu.puList[ind1].numZones);
    CHECK_EQUAL(1, pu.puList[ind2].numZones);
    CHECK_EQUAL(2, pu.puList[ind3].numZones);
    CHECK_EQUAL(0, pu.puList[ind4].numZones); // 0 means no zone limitation on pu.
    CHECK_EQUAL(0, pu.puList[ind5].numZones);
}

// Test connections TODO

// Test puvspr TODO