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
    Pu pu(fnames, c, 0);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    Zones zones(fnames, c);
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
    Pu pu(fnames, c, 0);
    pu.LoadSparseMatrix(spec, "data/puvspr_test1.dat");
    Zones zones(fnames, c);
    zones.BuildZoneContributions(spec, pu);

    Reserve r(spec, 3, 1); // 3 zones
    // set a solution
    r.InitializeSolution(pu.puno);
    r.ComputeSpeciesAmounts(pu, spec, zones);

    schange change1 = r.CheckChangeValue(0, 0, 1, pu, zones, spec, 1); //preZone 0, postZone 1.

    // check change value is 20 for species 1
    int spind1 = spec.LookupIndex(1);

    CHECK_EQUAL(1, change1.specListChangeTarget.size());
    CHECK_EQUAL(spind1, change1.specListChangeTarget[0].first);
    CHECK_EQUAL(-20, change1.specListChangeTarget[0].second);
}