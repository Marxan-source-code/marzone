#include "../reserve.hpp"
#include "CppUTest/TestHarness.h"

using namespace marzone;

TEST_GROUP(ReserveTestsGroup)
{
};

// ensure loadZones correctly parses.
TEST(ReserveTestsGroup, Reserve_InitializeSolution_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);
    Species spec(fnames);

    Reserve r(spec, 3, 1); // 3 zones
    CHECK_EQUAL(3, r.speciesAmounts.size());

    r.InitializeSolution(5);
    CHECK_EQUAL(5, r.solution.size());

    r.InitializeSolution(2);
    CHECK_EQUAL(2, r.solution.size());
}


// tests basic ComputeSpeciesAmounts with the _test1 suite.
TEST(ReserveTestsGroup, Reserve_ComputeSpeciesAmounts_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);
    Species spec(fnames);
    Zones zones(fnames, c);
    LoggerMock logger;
    Pu pu(fnames, c, 0, zones.zoneNames, logger);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    zones.BuildZoneContributions(spec, pu);

    Reserve r(spec, 3, 1); // 3 zones
    // set a solution
    r.InitializeSolution(pu.puno);

    // compute species amounts
    r.ComputeSpeciesAmounts(pu, spec, zones);

    int spind1 = spec.LookupIndex(1);
    int spind2 = spec.LookupIndex(2);
    int spind3 = spec.LookupIndex(3);

    CHECK(spind1 != -1 && spind2 != -1 && spind3 != -1);

    CHECK_EQUAL(70.5, r.speciesAmounts[spind1].amount); //20.0*1 + 50.5*1
    CHECK_EQUAL(15.8*0.5 + 118*0.5, r.speciesAmounts[spind2].amount); //15.8*0.5 + 118*0.5
    CHECK_EQUAL(212, r.speciesAmounts[spind3].amount); //200.0*1 + 12*1

    // Check occs
    CHECK_EQUAL(2, r.speciesAmounts[spind1].occurrence);
    CHECK_EQUAL(2, r.speciesAmounts[spind2].occurrence);
    CHECK_EQUAL(2, r.speciesAmounts[spind3].occurrence);

    r.solution[0] = 1; // flip zone to 1 indirectly
    r.ComputeSpeciesAmounts(pu, spec, zones); // recompute.

    CHECK_EQUAL(50.5, r.speciesAmounts[spind1].amount); //20.0*0 + 50.5*1
    CHECK_EQUAL(15.8*0.5 + 118*0.5, r.speciesAmounts[spind2].amount); // same as before
    CHECK_EQUAL(212, r.speciesAmounts[spind3].amount); // same as before.
}

// tests basic CheckChangeValue with the _test1 suite.
TEST(ReserveTestsGroup, Reserve_CheckChangeValue_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);
    Species spec(fnames);
    Zones zones(fnames, c);
    LoggerMock logger;
    Pu pu(fnames, c, 0, zones.zoneNames, logger);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    zones.BuildZoneContributions(spec, pu);

    Reserve r(spec, 3, 1); // 3 zones
    // set a solution
    r.InitializeSolution(pu.puno);
    r.ComputeSpeciesAmounts(pu, spec, zones);

    int spind1 = spec.LookupIndex(1);
    int spind2 = spec.LookupIndex(2);
    int spind3 = spec.LookupIndex(3);

    schange change1 = r.InitializeChange(spec, zones);
    r.CheckChangeValue(change1, 0, 0, 1, pu, zones, spec, 0); //puindex 0, preZone 0, postZone 1.

    // Check species amounts not changed
    CHECK_EQUAL(70.5, r.speciesAmounts[spind1].amount); //20.0*1 + 50.5*1
    CHECK_EQUAL(66.9, r.speciesAmounts[spind2].amount); //15.8*0.5 + 118*0.5
    CHECK_EQUAL(212, r.speciesAmounts[spind3].amount); //200.0*1 + 12*1

    // check change value is 20 for species 1
    CHECK_EQUAL(1, change1.specListChangeTarget.size());
    CHECK_EQUAL(spind1, change1.specListChangeTarget[0].first);
    CHECK_EQUAL(-20, change1.specListChangeTarget[0].second);

    // check change occurrence
    CHECK_EQUAL(1, change1.specListChangeOcc.size());
    CHECK_EQUAL(-1, change1.specListChangeOcc[0]);

    // Check others are 0 because zoneTargets don't exist
    CHECK_EQUAL(0, change1.zoneTargetChange.size());
    CHECK_EQUAL(0, change1.zoneOccChange.size());
    CHECK_EQUAL(0, change1.speciesClumpChange.size());

    // Check change overall values
    CHECK_EQUAL(20, change1.shortfall);
    CHECK_EQUAL(0, change1.penalty); // no penalties set for this example.
    CHECK_EQUAL(0, change1.cost);
    CHECK_EQUAL(0, change1.connection);
    CHECK_EQUAL(0, change1.total); // no penalty, cost or connection changes here.

    // Apply the change.
    r.ApplyChange(0, 1, change1, pu, zones, spec);

    CHECK_EQUAL(1, r.solution[0]);
    CHECK_EQUAL(50.5, r.speciesAmounts[spind1].amount); //20.0*0 + 50.5*1
    CHECK_EQUAL(15.8*0.5 + 118*0.5, r.speciesAmounts[spind2].amount); //15.8*0.5 + 118*0.5
    CHECK_EQUAL(212, r.speciesAmounts[spind3].amount); //200.0*1 + 12*1
}

// tests basic EvaluateObjectiveValue with the _test1 suite.
TEST(ReserveTestsGroup, Reserve_EvaluateObjectiveValue_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.inputdir = "";
    Costs c(fnames);
    Species spec(fnames);
    Zones zones(fnames, c);
    LoggerMock logger;
    Pu pu(fnames, c, 0, zones.zoneNames, logger);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    zones.BuildZoneContributions(spec, pu);

    Reserve r(spec, 3, 1); // 3 zones
    // set a solution to 0 0 0 0 0 
    r.InitializeSolution(pu.puno);

    // Evaluate objective value. and check shortfall, connections, cost
    r.EvaluateObjectiveValue(pu, spec, zones);

    // We know given this configuration, the species amounts are 70.5, 66.9, 212
    // we have not called SetProportionTargets so the regular (non proportion) targets are 100, 110 and 120
    double expected = (100-70.5)+(110-66.9) + (3-2); // 3-2 = shortfall occurrence for species 3
    CHECK_EQUAL(expected, r.objective.shortfall);

    // Since all pu are being used, and there's no zoneCost supplied, it should be total pu cost in this case.
    CHECK_EQUAL(6+66+666+2+2, r.objective.cost);

    // Connection cost - 0 in this case since no connections entered
    CHECK_EQUAL(0, r.objective.connection);
}