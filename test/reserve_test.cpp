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
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

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
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);
    Zones zones(fnames, c, logger);
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
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);
    Zones zones(fnames, c, logger);
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
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);
    Zones zones(fnames, c, logger);
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
    // note that the 5 missing costs get defaulted to 1, hence +5.
    CHECK_EQUAL(6+66+666+2+2+5, r.objective.cost);

    // Connection cost - 0 in this case since no connections entered
    CHECK_EQUAL(0, r.objective.connection);
}

// Test shortfall and penalty calculation when there is zone targets. 
TEST(ReserveTestsGroup, Reserve_EvaluateObjectiveValue_ZoneTargets_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.zonetargetname = "data/zonetarget_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

    // set dummy penalties
    vector<double> penalties {1000.0, 1000.0, 1000.0};
    spec.SetPenalties(penalties);

    Zones zones(fnames, c, logger);
    Pu pu(fnames, c, 0, zones.zoneNames, logger);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    zones.BuildZoneContributions(spec, pu);
    zones.BuildZoneTarget(spec, pu, fnames);

    Reserve r(spec, 3, 1); // 3 zones
    // set a solution to 0 0 0 0 0 
    r.InitializeSolution(pu.puno);

    // Evaluate objective value. and check shortfall, connections, cost
    r.EvaluateObjectiveValue(pu, spec, zones);

    double expected = (100-70.5)+(110-66.9) + (3-2); // regular shortfall
    expected += 100+ 500+ (1000-70.5); // shortfall with zone targets. Zone targ only supplied for species 1, but for all zones.
    CHECK_EQUAL(expected, r.objective.shortfall);
    CHECK(r.objective.penalty > 0);

    // Test change with pu2 to zone 1.
    schange change1 = r.InitializeChange(spec, zones);
    r.CheckChangeValue(change1, 0, 0, 1, pu, zones, spec, 0); //puindex 0, preZone 0, postZone 1.

    // ensure targets adjusted. Zone1 loses 20 and zone2 gains 20
    CHECK_EQUAL(2, change1.zoneTargetChange.size());
    CHECK_EQUAL(2, change1.zoneOccChange.size());
    CHECK_EQUAL(-20, change1.zoneTargetChange[0].second);
    CHECK_EQUAL(20, change1.zoneTargetChange[1].second);
    CHECK_EQUAL(-1, change1.zoneOccChange[0]);
    CHECK_EQUAL(1, change1.zoneOccChange[1]);

    // overall change in shortfall is positive (i.e. shortfall increased), since zone2 has a zonecontrib of 0
    // we should also ensure symmetry of changing a planning unit back and forth.
    double pre_change = change1.shortfall, pre_penalty = change1.penalty;
    CHECK(pre_change > 0);
    CHECK(pre_penalty > 0);

    r.ApplyChange(0, 1, change1, pu, zones, spec);
    r.CheckChangeValue(change1, 0, 1, 0, pu, zones, spec, 0); //puindex 0, preZone 1, postZone 0. Opposite of before.

    // ensure symmetry of values
    CHECK(change1.shortfall == -pre_change);
    CHECK(change1.penalty == -pre_penalty);
}

// This is a generic test that ensures flipping pu back and forth gives symmetrical results.
// This test does not check for specific values
TEST(ReserveTestsGroup, Reserve_Symmetry_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.zonetargetname = "data/zonetarget_test1.dat";
    fnames.inputdir = "";
    LoggerMock logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);

    // set dummy penalties
    vector<double> penalties {1000.0, 1000.0, 1000.0};
    spec.SetPenalties(penalties);

    Zones zones(fnames, c, logger);
    
    Pu pu(fnames, c, 0, zones.zoneNames, logger);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    zones.BuildZoneContributions(spec, pu);
    zones.BuildZoneTarget(spec, pu, fnames);

    Reserve r(spec, 3, 1); // 3 zones
    // set a solution to 0 0 0 0 0 
    r.InitializeSolution(pu.puno);

    // Evaluate objective value. and check shortfall, connections, cost
    r.EvaluateObjectiveValue(pu, spec, zones);

    // Set up objects needed for symmetry testing.
    uniform_int_distribution<int> randomDist(0, zones.zoneCount);
    mt19937 rngEngine(1); // arbitrary seed.
    int preZone, postZone;
    schange preChange = r.InitializeChange(spec, zones), postChange = r.InitializeChange(spec, zones);
    for (int i = 0; i < pu.puno; i++) {
        preZone = r.solution[i];
        postZone = pu.RtnValidZoneForPu(i, preZone, randomDist, rngEngine, 3);
        r.CheckChangeValue(preChange, i, preZone, postZone, pu, zones, spec, 0);
        r.ApplyChange(i, postZone, preChange, pu, zones, spec);
        r.CheckChangeValue(postChange, i, postZone, preZone, pu, zones, spec, 0);

        // Ensure symmetry in changes
        CHECK_EQUAL(preChange.total, -postChange.total);
        CHECK_EQUAL(preChange.shortfall, -postChange.shortfall);
        CHECK_EQUAL(preChange.penalty, -postChange.penalty);
        CHECK_EQUAL(preChange.cost, -postChange.cost);
        CHECK_EQUAL(preChange.connection, -postChange.connection);
    }
}