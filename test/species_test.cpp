
#include "../species.hpp"
#include "CppUTest/TestHarness.h"

using namespace marzone;

TEST_GROUP(SpeciesTestsGroup)
{
};

// test ensures species file gets correctly parsed.
TEST(SpeciesTestsGroup, ReadSpeciesData_fileparsing_test)
{
    sfname fnames = {}; // default, no cost file
    fnames.inputdir = "";
    fnames.specname = "data/species_test1.dat";

    Species sp(fnames);
    CHECK_EQUAL(3, sp.spno);

    // ensure ids were transformed into indices
    int ind1 = sp.LookupIndex(1);
    int ind2 = sp.LookupIndex(2);
    int ind3 = sp.LookupIndex(3);

    CHECK(ind1 != -1);
    CHECK(ind2 != -1);
    CHECK(ind3 != -1);
    CHECK_EQUAL(-1, sp.LookupIndex(4));

    // There was a prop measure so the flag should be set
    CHECK(sp.fSpecPROPLoaded);

    // Check regular targets set.
    CHECK_EQUAL(100, sp.specList[ind1].target);
    CHECK_EQUAL(110, sp.specList[ind2].target);
    CHECK_EQUAL(120, sp.specList[ind3].target);
    CHECK_EQUAL(1, sp.specList[ind1].targetocc);
    CHECK_EQUAL(2, sp.specList[ind2].targetocc);
    CHECK_EQUAL(3, sp.specList[ind3].targetocc);

    // Check target2 set to 0.
    CHECK_EQUAL(0, sp.specList[ind1].target2);
    CHECK_EQUAL(0, sp.specList[ind2].target2);
    CHECK_EQUAL(0, sp.specList[ind3].target2);

    vector<double> sums = {90, 100, 10};
    // Compute prop targets
    sp.SetSpeciesProportionTarget(sums);

    CHECK_EQUAL(27, sp.specList[ind1].target);
    CHECK_EQUAL(90, sp.specList[ind2].target);
    CHECK_EQUAL(5, sp.specList[ind3].target);

    // Check that fpf correctly read
    CHECK_EQUAL(1, sp.specList[ind1].spf);
    CHECK_EQUAL(100, sp.specList[ind2].spf);
    CHECK_EQUAL(1, sp.specList[ind3].spf);
}

// Tests how negative values are treated.
TEST(SpeciesTestsGroup, SetSpeciesDefaults_test)
{
    sfname fnames = {}; // default, no cost file
    fnames.inputdir = "";
    fnames.specname = "data/species_test2.dat";

    Species sp(fnames);

    // ensure ids were transformed into indices
    int ind1 = sp.LookupIndex(2);
    int ind2 = sp.LookupIndex(7);
    int ind3 = sp.LookupIndex(5);
    CHECK(ind1 != -1);
    CHECK(ind2 != -1);
    CHECK(ind3 != -1);

    // Check targets were read correctly. 
    // Check regular targets set.
    CHECK_EQUAL(-100, sp.specList[ind1].target);
    CHECK_EQUAL(110, sp.specList[ind2].target);
    CHECK_EQUAL(120, sp.specList[ind3].target);
    CHECK_EQUAL(1, sp.specList[ind1].targetocc);
    CHECK_EQUAL(2, sp.specList[ind2].targetocc);
    CHECK_EQUAL(-3, sp.specList[ind3].targetocc);
    CHECK_EQUAL(-1, sp.specList[ind2].spf);

    // Clean negatives
    sp.SetSpeciesDefaults();
    CHECK_EQUAL(0, sp.specList[ind1].target);
    CHECK_EQUAL(0, sp.specList[ind3].targetocc);
    CHECK_EQUAL(1, sp.specList[ind2].spf);
}

// Tests that setpenalties applies penalties. 
TEST(SpeciesTestsGroup, SetPenalties_test)
{
    sfname fnames = {}; // default, no cost file
    fnames.inputdir = "";
    fnames.specname = "data/species_test1.dat";

    Species sp(fnames);
    CHECK_EQUAL(3, sp.spno);

    // ensure ids were transformed into indices
    int ind1 = sp.LookupIndex(1);
    int ind2 = sp.LookupIndex(2);
    int ind3 = sp.LookupIndex(3);

    // Check penalty is default 0
    CHECK_EQUAL(0, sp.specList[ind1].penalty);
    CHECK_EQUAL(0, sp.specList[ind2].penalty);
    CHECK_EQUAL(0, sp.specList[ind3].penalty);

    vector<double> pens = {9.1, 1.01, 10000};
    // Compute prop targets
    sp.SetPenalties(pens);

    CHECK_EQUAL(9.1, sp.specList[0].penalty);
    CHECK_EQUAL(1.01, sp.specList[1].penalty);
    CHECK_EQUAL(10000, sp.specList[2].penalty);
}

// ensure that both the spf and fpf headers are parsed correctly. 
TEST(SpeciesTestsGroup, TestSpfFpfHeader_test)
{
    sfname fnames = {}; // default, no cost file
    fnames.inputdir = "";
    fnames.specname = "data/species_test3_spf.dat";

    Species sp(fnames);
    CHECK_EQUAL(2, sp.spno);

        // ensure ids were transformed into indices
    int ind1 = sp.LookupIndex(1);
    int ind2 = sp.LookupIndex(2);

    // Check spf assigned even with "spf" header
    CHECK_EQUAL(101.1, sp.specList[ind1].spf);
    CHECK_EQUAL(99, sp.specList[ind2].spf);
}