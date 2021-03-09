#pragma once

// Contains the spec class
// Handles all operations relating to species lookups

#include <fstream>
#include <map>
#include <string>
#include "common.hpp"
#include "util.hpp"
#include "logger.hpp"

namespace marzone {
using namespace std;

class Species {
    public:
    Species(sfname& fnames, Logger& logger): spno(0), gspno(0), aggexist(false), sepexist(false), fSpecPROPLoaded(false) {
        ReadSpeciesData(fnames.inputdir + fnames.specname, specList, spno, true, logger);

        if (!fnames.blockdefname.empty())
        {
            // read the general species file too
            ReadSpeciesData(fnames.inputdir + fnames.specname, genSpecList, gspno, false, logger);
        }
    }

    // given a species id, returns index of species. Returns -1 if not originally defined.
    // This operation is O(log) by number of species.
    int LookupIndex(int spid) {
        auto it = lookup.find(spid);
        if (it != lookup.end()) {
            return it->second;
        }

        return -1; // species not found.
    }

    // Prereq - spindex must be in range of 0 and < spno, otherwise behaviour undefined.
    void SetRichness(int spindex, int richness) {
        specList[spindex].richness = richness;
    }

    // Prereq - spindex must be in range of 0 and < spno, otherwise behaviour undefined.
    void SetOffset(int spindex, int offset) {
        specList[spindex].offset = offset;
    }

    // Sets a species target if proportion was set. 
    // No action if proportion target not used
    void SetSpeciesProportionTarget(vector<double>& speciesSums) {
        if (fSpecPROPLoaded) {
            for (int i = 0; i < speciesSums.size(); i++) {
                specList[i].target = specList[i].prop > 0 ? speciesSums[i]*specList[i].prop : 0.0;
            }
        }
    }

    void SetSpeciesBlockDefinitions(vector<double>& speciesSums) {
        int igsp, isp, ipu;

        for (igsp = 0; igsp < gspno; igsp++)
        {
            if (genSpecList[igsp].prop > 0) // deal with percentage in a different way
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].target < 0)
                    {
                        specList[isp].target = speciesSums[isp]* genSpecList[igsp].prop;
                    } // Setting target with percentage
            if (genSpecList[igsp].target > 0)
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].target < 0)
                        specList[isp].target = genSpecList[igsp].target;
            if (genSpecList[igsp].target2 > 0)
                for (isp = 0; isp < spno; isp++)
                {
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].target2 < 0)
                    {
                        specList[isp].target2 = genSpecList[igsp].target2;
                    }
                }
            if (genSpecList[igsp].targetocc > 0)
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].targetocc < 0)
                        specList[isp].targetocc = genSpecList[igsp].targetocc;
            if (genSpecList[igsp].spf > 0)
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].spf < 0)
                        specList[isp].spf = genSpecList[igsp].spf;
            // Percentage is not dealt with here yet. To do this identify
            // target species then determine their total abundance then set target
            // according to percentage
        }
    }

    // Scans through the species list and sets all defaults
    // * * *  Set Defaults * * *
    // If '-1' values haven't been set yet then this one will do it
    // Also refreshes aggexist and sepexist
    void SetSpeciesDefaults() {
        for (sspecies& spec : specList)
        {
            if (spec.target < 0)
                spec.target = 0;
            if (spec.target2 < 0)
                spec.target2 = 0;
            if (spec.targetocc < 0)
                spec.targetocc = 0;
            if (spec.spf < 0)
                spec.spf = 1;

            // agg exist
            if (spec.target2)
                aggexist = true;
        }
    }

    void SetPenalties(vector<double>& specPenalties) {
        for (int i = 0; i < spno; i++) {
            specList[i].penalty = specPenalties[i];
        }
    }

    void SetTotalAreas(vector<double>& specTotals) {
        for (int i = 0; i < spno; i++) {
            specList[i].totalarea = specTotals[i];
        }
    }

    void LoadCustomPenalties(string filename, Logger& logger) {
        ifstream fp = openFile(filename);
        string sLine;

        int i=0, spid;

        bool file_is_empty = true;
        for (int line_num = 1; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (line_num == 1)
            {
                if (is_like_numerical_data(sLine))
                    logger.ShowWarningMessage("File " + filename + " has no header in the first line.\n");
                else
                    continue;//skip header
            }
            if (sLine.empty())
                continue;

            stringstream ss = stream_line(sLine);
            ss >> spid >> specList[i].penalty;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");
            i++;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;

    }

    // Writes species penalties into filename
    void WritePenalties(string filename, int delimType) {
        ofstream myfile;
        myfile.open(filename);
        string delim = delimType == 3 ? "," : "    ";
        myfile << "spid" << delim << "penalty\n"; // write header

        // Ouput the Summary Statistics
        for (int i=0;i<spno;i++)
            myfile << specList[i].name << delim << specList[i].penalty << "\n";

        myfile.close();
    }

    void WriteSpeciesData(string filename, vector<reservespecies>& reserveStat) {
        ofstream myfile;
        myfile.open(filename);

        myfile << "i,name,type,sname,target,prop,targetocc,spf,penalty,amount,occurrence,clumps,target2,richness,offset\n";
        for (int i=0;i<spno;i++) {
            myfile << i << "," << specList[i].name << "," << specList[i].type << ","
            << specList[i].sname << "," << specList[i].target << "," << specList[i].prop << "," 
            << specList[i].targetocc << "," << specList[i].spf << "," << specList[i].penalty << "," << reserveStat[i].amount << "," << reserveStat[i].occurrence << "," 
            << reserveStat[i].clumps << "," << specList[i].target2 << "," << specList[i].richness << "," << specList[i].offset << "\n";
        }

        myfile.close();
    }

    void WriteTotalAreasAndOccs(string filename, vector<int>& TotalOccurrences, 
    vector<int>& TO_2, vector<int>& TO_3, vector<double>& TA_2, vector<double>& TA_3) {
        ofstream myfile;
        myfile.open(filename);

        myfile << "spname,spindex,totalarea,reservedarea,excludedarea,targetarea,totalocc,reservedocc,excludedocc,targetocc\n";
        for (int i = 0; i < spno; i++)
            myfile << specList[i].name << "," << i << "," << specList[i].totalarea << "," << TA_2[i] << "," << TA_3[i] << "," 
                << specList[i].target << "," << TotalOccurrences[i] <<"," << TO_2[i] << "," << TO_3[i] << "," << specList[i].targetocc << "\n";
        
        myfile.close();
    }

    unsigned spno;
    unsigned gspno;
    bool aggexist;
    bool sepexist;
    bool fSpecPROPLoaded;
    vector<sspecies> specList;
    map<int, int> lookup; // original specId -> index

    /* General Species Data */
    /* This function reads in a fixed file named file the general species info. It requires that
        species are typed and changes the value of various characteristics of species
        of that type. */
    vector<sgenspec> genSpecList; // contains species with only the general features.

    private:
    // precondition - T must be of type sgenspec or one of its inheritors
    template<typename T>
    void ReadSpeciesData(string filename, vector<T>& typeList, unsigned& count, bool populateLookup, Logger& logger) {
        ifstream fp = openFile(filename);
        string sLine, unusedHeader;
        stringstream errorBuf;

        int numvars = 10;
        /* Scan header */
        getline(fp, sLine);
        vector<string> headerNames = getFileHeaders(sLine, filename, errorBuf);

        if (!errorBuf.str().empty())
            logger.ShowErrorMessage(errorBuf.str());

        // check for sepnum and sepdistance
        if (find(headerNames.begin(), headerNames.end(), "sepnum") != headerNames.end()
            || find(headerNames.begin(), headerNames.end(), "sepdistance") != headerNames.end())
        {
            logger.ShowWarningMessage("Warning: sepdistance, sepnum features no longer supported and will be ignored in v4. Please use a previous version of Marxan with Zones.");
        }

        /* While there are still lines left populate specList */
        bool file_is_empty = true;
        for (int line_num = 2; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (sLine.empty())
                continue;

            T spectemp = {};
            /** Clear important species stats **/
            spectemp.target = -1;
            spectemp.type = -1;
            spectemp.spf = -1;
            spectemp.target2 = 0;
            spectemp.targetocc = -1;

            stringstream ss = stream_line(sLine);

            for (string temp : headerNames)
            { /* Tok out the next Variable */

                if (temp.compare("id") == 0)
                {
                    ss >> spectemp.name;
                }
                else if (temp.compare("name") == 0)
                {
                    string speciesname, word;
                    ss >> word;
                    speciesname += word;
                    //read the rest words of multiple word names
                    while(!ss.eof())
                    {
                        char delim;
                        ss.get(delim); //do not need to put delimeter back into stream
                        if (ss.eof())
                            break;
                        char first_letter = ss.peek();
                        if (isalpha(first_letter) || first_letter == '(' || first_letter == ')')
                        {
                            ss >> word;
                            speciesname += " ";
                            speciesname += word;
                        }
                        else
                            break;
                    }
                    spectemp.sname = speciesname;
                } 
                else if (temp.compare("type") == 0)
                {
                    ss >> spectemp.type;
                }
                else if (temp.compare("target") == 0)
                {
                    ss >> spectemp.target;
                }
                else if (temp.compare("prop") == 0)
                {
                    ss >> spectemp.prop;
                    if (spectemp.prop > 0)
                        fSpecPROPLoaded = true;
                }
                else if (temp.compare("fpf") == 0 || temp.compare("spf") == 0) // fpf, spf same header meaning.
                {
                    ss >> spectemp.spf;
                }
                else if (temp.compare("target2") == 0)
                {
                    ss >> spectemp.target2;
                }
                else if (temp.compare("targetocc") == 0)
                {
                    ss >> spectemp.targetocc;
                }
                else {
                    // un-enforced header
                    ss >> unusedHeader;
                }
            } /* looking for ivar different input variables */
            
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            typeList.push_back(spectemp);

            if (populateLookup)
                lookup[spectemp.name] = count;

            count++;
        } /* Scanning through each line of file */

        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;

    }
};

} // namespace marzone