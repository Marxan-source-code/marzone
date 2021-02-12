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
    LoggerMock logger;
    Costs c(fnames, logger);    
    // arbitrary zones map
    map<int, ZoneName> zoneNames;
    Pu pu(fnames, c, 0, zoneNames, logger);

    // Check costs are set correctly. Note that salmonfishing is missing so should be set to 1 by default.
    CHECK_EQUAL(5, pu.puno);

    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);

    CHECK(ind1 != -1);
    CHECK(ind2 != -1);
    CHECK(ind3 != -1);
    CHECK(ind4 != -1);
    CHECK(ind5 != -1);
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

    // No connections entered, ensure all connections are 0
    CHECK(!pu.connectionsEntered);
    CHECK_EQUAL(0, pu.ConnectionCost1(ind1));
    CHECK_EQUAL(0, pu.ConnectionCost1(ind2));
    CHECK_EQUAL(0, pu.ConnectionCost1(ind3));
    CHECK_EQUAL(0, pu.ConnectionCost1(ind4)); 
    CHECK_EQUAL(0, pu.ConnectionCost1(ind5)); 
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
    LoggerMock logger;
    Costs c(fnames, logger);

    // build zones map
    map<int, ZoneName> zoneNames;
    zoneNames[1] = {"available", 0};
    zoneNames[2] = {"available", 1};
    zoneNames[3] = {"available", 2};

    Pu pu(fnames, c, 0, zoneNames, logger);
    
    // Ensure locked pus get returned as such.
    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);

    // check fPuLock is set to the zone INDEX (not id)
    CHECK_EQUAL(2, pu.puLockCount);
    CHECK_EQUAL(0, pu.puList[ind1].iPULock);
    CHECK_EQUAL(1, pu.puList[ind2].iPULock);
    CHECK_EQUAL(-1, pu.puList[ind3].iPULock);
    CHECK_EQUAL(1, pu.puList[ind1].fPULock);
    CHECK_EQUAL(1, pu.puList[ind2].fPULock);
    CHECK_EQUAL(0, pu.puList[ind3].fPULock);

    // check puzone correct entered for each pu.
    CHECK_EQUAL(1, pu.puList[ind1].numZones);
    CHECK_EQUAL(1, pu.puList[ind2].numZones);
    CHECK_EQUAL(2, pu.puList[ind3].numZones);
    CHECK_EQUAL(0, pu.puList[ind4].numZones); // 0 means no zone limitation on pu.
    CHECK_EQUAL(0, pu.puList[ind5].numZones);

    // test valid indices
    CHECK_EQUAL(3, pu.validPuIndices.size());
    CHECK_EQUAL(ind3, pu.validPuIndices[0]);
    CHECK_EQUAL(ind4, pu.validPuIndices[1]);
    CHECK_EQUAL(ind5, pu.validPuIndices[2]);
}

// Test connections and connectionCost1
TEST(PuTestsGroup, Connections_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.connectionname = "data/connections_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    
    // arbitrary zones map
    map<int, ZoneName> zoneNames;
    Pu pu(fnames, c, 0, zoneNames, logger);

    // check connectionsEntered flag
    CHECK(pu.connectionsEntered);

    CHECK_EQUAL(5, pu.connections.size());

    // Check neighbours entered correctly.
    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);

    CHECK_EQUAL(2, pu.connections[ind1].nbrno);
    CHECK_EQUAL(2, pu.connections[ind2].nbrno);
    CHECK_EQUAL(2, pu.connections[ind3].nbrno);
    CHECK_EQUAL(0, pu.connections[ind4].nbrno);
    CHECK_EQUAL(0, pu.connections[ind5].nbrno);

    // Test cost
    CHECK_EQUAL(2.0, pu.ConnectionCost1(ind4)); // pu11 has no neighbours, only itself.
    CHECK_EQUAL(100.7, pu.ConnectionCost1(ind1));
    CHECK_EQUAL(101.4, pu.ConnectionCost1(ind2));
    CHECK_EQUAL(2.3, pu.ConnectionCost1(ind3));
    CHECK_EQUAL(0.0, pu.ConnectionCost1(ind5)); // pu8 has no cost
}

// Test puvspr 
TEST(PuTestsGroup, Puvspr_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

    // arbitrary zones map
    map<int, ZoneName> zoneNames;
    Pu pu(fnames, c, 0, zoneNames, logger);

    //load puvspr file
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");

    // Check counts - it should ignore the last 2 rows which contains a puid/spid not in the original set.
    CHECK_EQUAL(6, pu.puvspr.size());

    // Check amounts for all species
    vector<double> totals = pu.TotalSpeciesAmount(spec.spno); // spno = 3
    CHECK_EQUAL(3, totals.size());

    int spind1 = spec.LookupIndex(1);
    int spind2 = spec.LookupIndex(2);
    int spind3 = spec.LookupIndex(3);

    // check species totals add up.
    CHECK_EQUAL(70.5, totals[spind1]);
    CHECK_EQUAL(133.8, totals[spind2]);
    CHECK_EQUAL(212, totals[spind3]);

    // Check occurence totals
    vector<int> occTotals = pu.TotalOccurrenceAmount(spec.spno);
    CHECK_EQUAL(3, occTotals.size());
    CHECK_EQUAL(2, occTotals[spind1]);
    CHECK_EQUAL(2, occTotals[spind2]);
    CHECK_EQUAL(2, occTotals[spind3]);
}

// Test puvspr get functions
// uses the same data files as the function above.
TEST(PuTestsGroup, Puvspr_RtnAmountSpecAtPu_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

    // arbitrary zones map
    map<int, ZoneName> zoneNames;
    Pu pu(fnames, c, 0, zoneNames, logger);

    //load puvspr file
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");

    // Get all the lookup indices
    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);
    int spind1 = spec.LookupIndex(1);
    int spind2 = spec.LookupIndex(2);
    int spind3 = spec.LookupIndex(3);

    // Check richnesses
    CHECK_EQUAL(1, pu.puList[ind1].richness);
    CHECK_EQUAL(1, pu.puList[ind2].richness);
    CHECK_EQUAL(2, pu.puList[ind3].richness);
    CHECK_EQUAL(1, pu.puList[ind4].richness);
    CHECK_EQUAL(1, pu.puList[ind5].richness);

    // Check indices
    CHECK(pu.RtnIndexSpecAtPu(ind1, spind1) != -1);
    CHECK_EQUAL(-1, pu.RtnIndexSpecAtPu(ind2, spind1));
    CHECK(pu.RtnIndexSpecAtPu(ind3, spind1) != -1);
    CHECK_EQUAL(-1, pu.RtnIndexSpecAtPu(ind3, spind2));
    CHECK(pu.RtnIndexSpecAtPu(ind5, spind2) != -1);
    CHECK(pu.RtnIndexSpecAtPu(ind4, spind2) != -1);
    CHECK(pu.RtnIndexSpecAtPu(ind2, spind3) != -1);
    CHECK(pu.RtnIndexSpecAtPu(ind3, spind3) != -1);
    CHECK_EQUAL(-1, pu.RtnIndexSpecAtPu(ind1, spind3));

    // Check most assignments including some zeroes
    CHECK_EQUAL(20.0, pu.RtnAmountSpecAtPu(ind1, spind1));
    CHECK_EQUAL(0, pu.RtnAmountSpecAtPu(ind2, spind1));
    CHECK_EQUAL(50.5, pu.RtnAmountSpecAtPu(ind3, spind1));
    CHECK_EQUAL(0, pu.RtnAmountSpecAtPu(ind3, spind2));
    CHECK_EQUAL(15.8, pu.RtnAmountSpecAtPu(ind5, spind2));
    CHECK_EQUAL(118.0, pu.RtnAmountSpecAtPu(ind4, spind2));
    CHECK_EQUAL(200.0, pu.RtnAmountSpecAtPu(ind2, spind3));
    CHECK_EQUAL(12.0, pu.RtnAmountSpecAtPu(ind3, spind3));
    CHECK_EQUAL(0, pu.RtnAmountSpecAtPu(ind1, spind3));
}

// Test puvspr RtnAmountAllSpecAtPu functions
// uses the same data files as the function above.
TEST(PuTestsGroup, RtnAmountAllSpecAtPu_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

    // arbitrary zones map
    map<int, ZoneName> zoneNames;
    Pu pu(fnames, c, 0, zoneNames, logger);

    //load puvspr file
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");

    // Get all the lookup indices
    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);

    // check totals for each spindex 
    vector<double> expected1 {20, 0, 0};
    vector<double> expected2 {0, 0, 200.0};
    vector<double> expected3 {50.5, 0, 12};
    vector<double> expected4 {0, 118, 0};
    vector<double> expected5 {0, 15.8, 0};

    CHECK(expected1 == pu.RtnAmountAllSpecAtPu(ind1, spec.spno));
    CHECK(expected2 == pu.RtnAmountAllSpecAtPu(ind2, spec.spno));
    CHECK(expected3 == pu.RtnAmountAllSpecAtPu(ind3, spec.spno));
    CHECK(expected4 == pu.RtnAmountAllSpecAtPu(ind4, spec.spno));
    CHECK(expected5 == pu.RtnAmountAllSpecAtPu(ind5, spec.spno));
}

// tests getPuAmountsSorted which is an important function for penalty calculations.
TEST(PuTestsGroup, getPuAmountsSorted_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

    // arbitrary zones map
    map<int, ZoneName> zoneNames;
    Pu pu(fnames, c, 0, zoneNames, logger);

    //load puvspr file
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");

    // each vect in penaltySorted should be in amount/cost ratio.
    vector<vector<lockedPenaltyTerm>> lockedPenalties;
    vector<vector<penaltyTerm>> penaltySorted = pu.getPuAmountsSorted(spec.spno, lockedPenalties);

    CHECK_EQUAL(3, penaltySorted.size());

    int spind1 = spec.LookupIndex(1);
    int spind2 = spec.LookupIndex(2);
    int spind3 = spec.LookupIndex(3);

    // ensure number of terms is correct
    CHECK_EQUAL(2, penaltySorted[spind1].size());
    CHECK_EQUAL(2, penaltySorted[spind2].size());
    CHECK_EQUAL(2, penaltySorted[spind3].size());

    // ensure ordering of each vector<penaltyTerm> is correct.
    CHECK_EQUAL(20.0/7, penaltySorted[spind1][0].amount/penaltySorted[spind1][0].cost);
    CHECK_EQUAL(50.5/667, penaltySorted[spind1][1].amount/penaltySorted[spind1][1].cost);
    CHECK_EQUAL(118.0/3, penaltySorted[spind2][0].amount/penaltySorted[spind2][0].cost);
    CHECK_EQUAL(15.8/3, penaltySorted[spind2][1].amount/penaltySorted[spind2][1].cost);
    CHECK_EQUAL(200.0/67, penaltySorted[spind3][0].amount/penaltySorted[spind3][0].cost);
    CHECK_EQUAL(12.0/667, penaltySorted[spind3][1].amount/penaltySorted[spind3][1].cost);

    // No locked penalties since file is not supplied
    CHECK_EQUAL(3, lockedPenalties.size());
    CHECK_EQUAL(0, lockedPenalties[0].size());
    CHECK_EQUAL(0, lockedPenalties[1].size());
    CHECK_EQUAL(0, lockedPenalties[2].size());
}

// testing valud zone return
TEST(PuTestsGroup, RtnValidZoneForPu_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.pulockname = "data/pulock_test1.dat";
    fnames.puzonename = "data/puzone_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);

    // build zones map
    map<int, ZoneName> zoneNames;
    zoneNames[1] = {"available", 0};
    zoneNames[2] = {"available", 1};
    zoneNames[3] = {"available", 2};

    Pu pu(fnames, c, 0, zoneNames, logger);

    int ind1 = pu.LookupIndex(2);
    int ind2 = pu.LookupIndex(5);
    int ind3 = pu.LookupIndex(10);
    int ind4 = pu.LookupIndex(11);
    int ind5 = pu.LookupIndex(8);

    uniform_int_distribution<int> randomDist(0, 99);
    mt19937 rng(1); // arbitrary seed

    // Ensure this function returns correct zones.
    // pu2 and 5 are locked, and 10 is zoned to 0 or 2. 11 and 8 can be any zone.
    CHECK_EQUAL(0, pu.RtnValidZoneForPu(ind1, 0, randomDist, rng, 3));
    CHECK_EQUAL(1, pu.RtnValidZoneForPu(ind2, 0, randomDist, rng, 3));
    CHECK_EQUAL(2, pu.RtnValidZoneForPu(ind3, 0, randomDist, rng, 3)); // zoneind2 only because starts from 0.

    // test 100 times. 
    int randZone = 0;
    for (int i = 0; i < 100; i++) {
        randZone = pu.RtnValidZoneForPu(ind4, 0, randomDist, rng, 3);
        CHECK(randZone == 1 || randZone == 2);

        randZone = pu.RtnValidZoneForPu(ind5, 0, randomDist, rng, 3);
        CHECK(randZone == 1 || randZone == 2);
    }
}