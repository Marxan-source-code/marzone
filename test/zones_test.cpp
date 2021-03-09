#include "../zones.hpp"
#include "CppUTest/TestHarness.h"

using namespace marzone;

TEST_GROUP(ZonesTestsGroup)
{
};

// ensure loadZones correctly parses.
TEST(ZonesTestsGroup, LoadZones_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.inputdir = "";
    Logger logger;
    Costs c(fnames, logger);
    Zones z(fnames, c, logger);

    CHECK_EQUAL(3, z.zoneCount);

    // check indices assigned properly
    CHECK_EQUAL(3, z.zoneNames.size());
    CHECK_EQUAL("available", z.zoneNames[1].name);
    CHECK_EQUAL("partial", z.zoneNames[2].name);
    CHECK_EQUAL("reserve", z.zoneNames[3].name);
    CHECK_EQUAL(0, z.zoneNames[1].index);
    CHECK_EQUAL(1, z.zoneNames[2].index);
    CHECK_EQUAL(2, z.zoneNames[3].index);

    // Check backward index map
    CHECK_EQUAL(3, z.zoneNameIndexed.size());
    CHECK_EQUAL(1, z.zoneNameIndexed[0].id);
    CHECK_EQUAL(2, z.zoneNameIndexed[1].id);
    CHECK_EQUAL(3, z.zoneNameIndexed[2].id);
    CHECK_EQUAL("available", z.zoneNameIndexed[0].name);
    CHECK_EQUAL("partial", z.zoneNameIndexed[1].name);
    CHECK_EQUAL("reserve", z.zoneNameIndexed[2].name);
}

// Tests default connection and zone cost values.
TEST(ZonesTestsGroup, DefaultCosts_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.inputdir = "";
    Logger logger;
    Costs c(fnames, logger);
    Zones z(fnames, c, logger);

    // Ensure connectionCosts and zoneCosts are non-empty and set to 1 if file NOT supplied.
    CHECK_EQUAL(c.costCount*z.zoneCount, z.zoneCost.size());
    CHECK_EQUAL(z.zoneCount*z.zoneCount, z.zoneConnectionCost.size());

    for (int i= 0; i < z.zoneCost.size(); i++)
    {
        CHECK_EQUAL(1, z.zoneCost[i]);
    }

    for (int i= 0; i < z.zoneConnectionCost.size(); i++)
    {
        CHECK_EQUAL(1, z.zoneConnectionCost[i]);
    }
}

// Tests that zoneContrib are assigned to 1 by default
TEST(ZonesTestsGroup, DefaultZoneContrib_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.inputdir = "";
    Logger logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);
    Zones z(fnames, c, logger);
    Pu pu(fnames, c, 0, z.zoneNames, logger);

    // Build zone contribs
    z.BuildZoneContributions(spec, pu);

    // Check zone contribs are 1
    // Available zone 0 should be set to 0 (zoneid1)
    CHECK_EQUAL(0, z.GetZoneContrib(0,0));
    CHECK_EQUAL(1, z.GetZoneContrib(0,1));
    CHECK_EQUAL(1, z.GetZoneContrib(0,2));
    CHECK_EQUAL(0, z.GetZoneContrib(1,0));
    CHECK_EQUAL(1, z.GetZoneContrib(1,1));
    CHECK_EQUAL(1, z.GetZoneContrib(1,2));
    CHECK_EQUAL(0, z.GetZoneContrib(2,0));
    CHECK_EQUAL(1, z.GetZoneContrib(2,1));
    CHECK_EQUAL(1, z.GetZoneContrib(2,2));

    // Check pu version too for a sample of pu
    CHECK_EQUAL(0, z.GetZoneContrib(0, pu.puno, 0, 0));
    CHECK_EQUAL(1, z.GetZoneContrib(1, pu.puno, 2, 2));
    CHECK_EQUAL(1, z.GetZoneContrib(2, pu.puno, 1, 1));
}

// ensure zone contrib1 correctly parses.
TEST(ZonesTestsGroup, LoadZoneContrib_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.inputdir = "";
    Logger logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);
    Zones z(fnames, c, logger);
    Pu pu(fnames, c, 0, z.zoneNames, logger);

    // Build zone contribs
    z.BuildZoneContributions(spec, pu);

    // species index
    int ind1 = spec.LookupIndex(1);
    int ind2 = spec.LookupIndex(2);

    // Here the file exists, so missing contribs are 0 instead of 1. 
    // Final row of file should be ignored
    // species 2 contribs are all set to 0.5
    CHECK_EQUAL(1, z.GetZoneContrib(ind1,0));
    CHECK_EQUAL(0, z.GetZoneContrib(ind1,1));
    CHECK_EQUAL(1, z.GetZoneContrib(ind1,2));
    CHECK_EQUAL(0.5, z.GetZoneContrib(ind2,0));
    CHECK_EQUAL(0.5, z.GetZoneContrib(ind2,1));
    CHECK_EQUAL(0.5, z.GetZoneContrib(ind2,2));
}

TEST(ZonesTestsGroup, LoadZoneTargets_fileparsing_test)
{
    sfname fnames = {};
    fnames.costsname = "data/costs_test1.dat";
    fnames.specname = "data/species_test1.dat";
    fnames.puname = "data/pu_test1.dat";
    fnames.zonesname = "data/zones_test1.dat";
    fnames.zonecontribname = "data/zonecontrib_test1.dat";
    fnames.zonetargetname = "data/zonetarget_test1.dat";
    fnames.inputdir = "";
    Logger logger;
    Costs c(fnames, logger);
    Species spec(fnames, logger);
    Zones z(fnames, c, logger);
    Pu pu(fnames, c, 0, z.zoneNames, logger);

    z.BuildZoneTarget(spec, pu, fnames, logger);

    // ensure 4th and 5th row skipped because invalid species id 
    CHECK_EQUAL(3, z.zoneTarget.size());
    CHECK_EQUAL(3, z.zoneTarget[0].size());
    CHECK_EQUAL(3, z.zoneTarget[1].size());
    CHECK_EQUAL(3, z.zoneTarget[2].size());

    // ensure targets loaded correctly
    CHECK_EQUAL(1000, z.zoneTarget[0][0].target);
    CHECK_EQUAL(500, z.zoneTarget[0][1].target);
    CHECK_EQUAL(100, z.zoneTarget[0][2].target);
    CHECK_EQUAL(0, z.zoneTarget[1][0].target);
    CHECK_EQUAL(0, z.zoneTarget[1][1].target);
    CHECK_EQUAL(0, z.zoneTarget[1][2].target);
    CHECK_EQUAL(0, z.zoneTarget[2][0].target);
    CHECK_EQUAL(0, z.zoneTarget[2][1].target);
    CHECK_EQUAL(0, z.zoneTarget[2][2].target);
}
