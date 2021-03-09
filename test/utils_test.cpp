
#include "../util.hpp"
#include "CppUTest/TestHarness.h"

using namespace marzone;

TEST_GROUP(UtilTestsGroup)
{
};

// Test cleanDirectoryString which adds a '/' to end of a given path.
TEST(UtilTestsGroup, cleanDirectoryString_test1)
{
    CHECK_EQUAL("out/", cleanDirectoryString("out"));
    CHECK_EQUAL("/", cleanDirectoryString(""));
    CHECK_EQUAL("/", cleanDirectoryString("/"));
    CHECK_EQUAL("\\", cleanDirectoryString("\\"));
    CHECK_EQUAL("out space/", cleanDirectoryString("out space"));
    CHECK_EQUAL("out space/", cleanDirectoryString("out space/"));
    CHECK_EQUAL("out space\\", cleanDirectoryString("out space\\"));
}

// test regular space is stripped
TEST(UtilTestsGroup, trim_test1)
{
    std::string original = " original s ";
    trim(original);
    CHECK_EQUAL("original s", original);
}

// Test tab is stripped
TEST(UtilTestsGroup, trim_test2)
{
    std::string original = " original s\t";
    trim(original);
    CHECK_EQUAL("original s", original);

    original = "\toriginal\t\t";
    trim(original);
    CHECK_EQUAL("original", original);
}

// Testing readInputOption
TEST(UtilTestsGroup, readInputOption_test1)
{
    stringstream w, e; 
    bool present;
    vector<string> lines;

    // Construct file lines
    lines.push_back("VAR1 1.23345");
    lines.push_back("VAR2 \t100"); // test uneven spacing
    lines.push_back("VAR3 filename with spaces");
    lines.push_back("VAR4  "); // empty
    lines.push_back("VAR5 \"filename with spaces\"");
    lines.push_back("VAR1 2.23345"); // repeat variable

    // read non existing var
    double var0;
    readInputOption(lines, "VAR0", var0, false, present, w, e);
    CHECK_EQUAL(false, present);
    CHECK(e.str().empty()); // no error message

    // read non existing crit var
    readInputOption(lines, "VAR0", var0, true, present, w, e);
    CHECK_EQUAL(false, present);
    CHECK(!e.str().empty()); // error message exists

    // read var 1
    double var1;
    readInputOption(lines, "VAR1", var1, true, present, w, e);
    CHECK_EQUAL(true, present);
    CHECK_EQUAL(2.23345, var1); // check that it's the last variable
    CHECK(!w.str().empty()); // check there's a warning message

    // read var2 - pretty standard
    int var2;
    readInputOption(lines, "VAR2", var2, true, present, w, e);
    CHECK_EQUAL(100, var2);

    // read var 3 and ensure full text is read.
    string var3;
    readInputOption(lines, "VAR3", var3, false, present, w, e);
    CHECK_EQUAL("filename with spaces", var3);

    // read var 5 and ensure full text is read.
    string var5;
    readInputOption(lines, "VAR5", var5, false, present, w, e);
    CHECK_EQUAL("\"filename with spaces\"", var5);

    // clear debug buffers
    e.clear();
    w.clear(); 

    // read var4 - the empty var
    string var4;
    readInputOption(lines, "VAR4", var4, false, present, w, e);
    CHECK(var4.empty());
    CHECK(!w.str().empty());

    // try to read invalid type
    double var6;
    readInputOption(lines, "VAR3", var6, false, present, w, e);
    CHECK(!e.str().empty());
}

TEST(UtilTestsGroup, getFileHeaders_test1) {
    char tempHeader[] = "h1,h2,h3";
    char invalidHeader[] = "h1,h1";
    char tempHeader2[] = "h1\t\th2";
    char tempHeader3[] = "h1, h2"; // test extra chars get trimmed.
    char tempHeader4[] = "id, prop, fpf, name, target, targetocc";
    stringstream errorBuf;

    vector<string> t1 = getFileHeaders(tempHeader, "", errorBuf);
    CHECK_EQUAL(3, t1.size());
    
    getFileHeaders(invalidHeader, "", errorBuf);
    CHECK(!errorBuf.str().empty());

    vector<string> t2 = getFileHeaders(tempHeader2, "", errorBuf);
    CHECK_EQUAL(2, t2.size());

    vector<string> t3 = getFileHeaders(tempHeader3, "", errorBuf);
    CHECK_EQUAL(2, t3.size());
    CHECK(t3[0] == "h1" || t3[0] == "h2");
    CHECK(t3[1] == "h1" || t3[1] == "h2");

    vector<string> t4 = getFileHeaders(tempHeader4, "", errorBuf);
    CHECK_EQUAL(6, t4.size());
    CHECK("targetocc" == t4[5]);
}

TEST(UtilTestsGroup, Range_test1) {
    vector<int> v1 = Range(0,1);
    CHECK_EQUAL(1, v1.size());
    CHECK_EQUAL(0, v1[0]);

    v1 = Range(0,0);
    CHECK_EQUAL(0, v1.size());

    v1 = Range(5,7);
    CHECK_EQUAL(2, v1.size());
    CHECK_EQUAL(5, v1[0]);
}

TEST(UtilTestsGroup, intToPaddedString_test)
{
    CHECK_EQUAL("01", intToPaddedString(1,2));
    CHECK_EQUAL("0100", intToPaddedString(100,4));
    CHECK_EQUAL("100", intToPaddedString(100,2));
    CHECK_EQUAL("1", intToPaddedString(1,1));
    CHECK_EQUAL("0000056", intToPaddedString(56,7));
    CHECK_EQUAL("56", intToPaddedString(56,0));
    CHECK_EQUAL("1", intToPaddedString(1,0));
    CHECK_EQUAL("0", intToPaddedString(0,0));
}