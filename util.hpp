#pragma once

// Handles utility operations like parsing arguments, strings etc.
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include "common.hpp"
#include <iostream>
#include <iomanip>

namespace marzone {

using namespace std;

inline
// parses command line options specifically for marxan
void HandleOptions(int argc, char *argv[], string& sInputFileName,int& marxanIsSecondary)
{
    if (argc>4) // if more than three commandline argument then exit
    { 
        string programName(argv[0]);
        throw invalid_argument(programName + " usage: " + programName + " -[o] -[c] [input file name]\n"); 
    }

    for (int i=1;i<argc;i++)
    {
        if (argv[i][0] == '/' || argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
            // No longer supporting oldstyle in new marzone
            //case '0':
            //case 'o':
            //case 'O':
            //    oldStyle = 1;
            //    break;
            case 'C':
            case 'c':
            case 'S':
            case 's':
                marxanIsSecondary = 1;
                break;
            default:
                fprintf(stderr,"unknown option %s\n",argv[i]);
            }
        }
        else
            sInputFileName = argv[i]; // If not a -option then must be input.dat name
    }  // Deal with all arguments
}

inline
void SetDefaultOptions(srunoptions &runoptions, sanneal &anneal, sfname &fnames) {
    /* Setup all of the default parameter variables */
    runoptions.prop = 0;
    anneal.type = 1;
    anneal.iterations = 0;
    anneal.Tinit = 1;
    anneal.Tcool = 0;
    anneal.Titns = 1;
    runoptions.iseed = -1;
    runoptions.costthresh = 0;
    runoptions.tpf1 = 0;
    runoptions.tpf2 = 0;
    runoptions.repeats = 0;
    runoptions.PopulationAnnealingOn = false;
    fnames.saverun = 0;
    fnames.savebest = 0;
    fnames.savesum = 0;
    fnames.savesen = 0;
    fnames.savespecies = 0;
    fnames.savesumsoln = 0;
    fnames.savepenalty = 0;
    fnames.savetotalareas = 0;
    fnames.savesolutionsmatrix = 0;
    fnames.solutionsmatrixheaders = 1;
    fnames.saveannealingtrace = 0;
    fnames.annealingtracerows = 0;
    fnames.suppressannealzones = 0;
    fnames.itimptracerows = 0;
    fnames.saveitimptrace = 0;
    fnames.savespeciesdata = 0;
    fnames.savezoneconnectivitysum = 0;
    fnames.savelog = 0;
    fnames.savename = "temp";
    runoptions.misslevel = 1;
    runoptions.heurotype = 1;
    runoptions.clumptype = 0;
    runoptions.verbose = 1;
    runoptions.blm = 1;
}

inline
/* Open file and store all lines in vector */
vector<string> GetFileLines(string& filename, stringstream& error) {
    vector<string> fileLines;
    ifstream fp(filename.c_str());
    // Check if object is valid
    if(!fp)
    {
        error << "File " << filename << " not found\nAborting Program.\n\n";
        return vector<string>(); // return empty.
    }

    string line;
    while (getline(fp, line))
    {
        // Store non empty lines
        if(line.size() > 0)
            fileLines.push_back(line);
    }
    fp.close();

    return fileLines;
}

// String trimming functions credited to https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// trim from start (in place)
inline
void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline
void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// Adds '/' to end of dir string if not existent
inline
std::string cleanDirectoryString(std::string dirName) {
    if (dirName.back() != '\\' && dirName.back() != '/') {
        return dirName + "/";
    }

    return dirName;
}

// given the parsed value of an option and its given value storage, parse it in.
// For string types, the overload below is given.
inline
void readInputOptionValue(stringstream& parsed, string& value) {
    value = parsed.str();
}

template<class T> inline
void readInputOptionValue(stringstream& parsed, T& value) {
    parsed >> value;
}

// read value for a single parameter specified in file lines "infile"
template<class T> inline
void readInputOption(vector<string>& infile, string varname, T& value, bool crit, bool& present, stringstream& warningBuf, stringstream& errorBuf)
// Given a varname, scans the infile to see if this field was configured with primitive type T. 
// Here the lines of the file are supplied as a vector of strings.
// If configured, returns the value in "value". 
// crit = whether variable is mandatory
// present = stores whether variable was found.
{
    int foundit = 0;
    string modifiedVar = varname;

    for (string line : infile)
    {   /* loop through file looking for varname */
        // if line contains varname
        if (line.find(modifiedVar)!=string::npos) {
            // try to parse the value next to the variable by erasing first modifiedVar characters of line
            string varValue = line.erase(0, modifiedVar.size());
            trim(varValue);

            if (varValue.empty()) { // variable defined but no value - keep looping.
                warningBuf << "WARNING: found bad/empty value for variable " << varname << ". Value ignored\n";
                continue;
            }

            // read string representation to specified type T
            stringstream linestream(varValue);
            readInputOptionValue(linestream, value);

            if (linestream.fail()) {
                // some error occurred when reading in value 
                errorBuf << "Invalid parameter type for given param " << varname << ". Expecting type: " << typeid(value).name() << "\n";
                continue;
            }

            foundit++;
        }
    }

    if (!foundit)
        if (crit)
           errorBuf << "Unable to find a valid value for " << varname << " in input file.\n";

    if (foundit > 1)
        warningBuf << "WARNING variable: " << varname << " appears more than once in the input file. Final value taken\n";

    present = foundit > 0 ? true : false;
} // readInputOption

inline
ifstream openFile(string filename) {
    ifstream fp;
    fp.open(filename);
    if (!fp.is_open())
        throw invalid_argument("Cannot find or open file " + filename + ". Aborting program.");
    return fp;
}

inline 
string getFileSuffix(int mode) {
    if (mode == 3) {
        return ".csv";
    }
    else if (mode == 2) {
        return ".txt";
    }
    else {
        return ".dat";
    }
}

inline 
vector<int> Range(int start, int end) {
    vector<int> r;
    for (int i = start; i < end; i++)
        r.push_back(i);
    return r;
}

inline
std::vector<std::string> get_tokens(const std::string& str)
{
    static const std::string delimeters(" ,;:^*\"/\t\'\\\n");
    std::vector<std::string> tokens;
    std::string word;
    for (char ch : str)
    {
        if (delimeters.find_first_of(ch) == std::string::npos)
            word.push_back(ch);
        else
        {
            if (!word.empty())
            {
                tokens.push_back(word);
                word.clear();
            }
        }
    }
    if (!word.empty())
        tokens.push_back(word);
    return tokens;
}

inline
std::stringstream stream_line(const std::string& str)
{
    static const std::string delimeters(" ,;:^*\"/\t\'\\\n");
    std::stringstream ss;
    for (char ch : str)
    {
        if (delimeters.find_first_of(ch) == std::string::npos)
            ss << ch;
        else
        {
            ss << ' ';
        }
    }
    return ss;
}

inline
bool is_like_numerical_data(const std::string& str)
{
    //check if there are chars other then delimeters and ones used to represent numbers
    size_t letter_pos = str.find_first_not_of(" ,;:^*\"/\t\'\\\n\r.+-Ee0123456789");
    if (letter_pos != std::string::npos)
        return false;
    return true;
}

inline
// Toks the given header line and returns an ordered list of header names. 
// We could use set instead of vector, but the size of these sets is very small (i.e under 10) so vector should be equal perf.
vector<string> getFileHeaders(const string& header, const string& filename, stringstream& errorBuf) {
    vector<string> headers;
    vector<string> tokens = get_tokens(header);
    for (string cleaned : tokens)
    {
        trim(cleaned);
        if (find(headers.begin(), headers.end(), cleaned) == headers.end()) {
            // Add to list if not already present.
            headers.push_back(cleaned);
        }
        else {
            errorBuf << "Header name " << cleaned << " has been defined twice in data file " << filename << ".\n";
        }
    }

    return headers;
}

// converts number to a string padded with leading zeros
// does nothing if stringLength is less than the digits in number.
inline 
std::string intToPaddedString(int number, int stringLength)
{
    std::ostringstream ss;
    ss << std::setw(stringLength) << std::setfill('0') << number;
    return ss.str();
}

} // namespace marzone