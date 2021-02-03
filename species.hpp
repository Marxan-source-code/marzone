#pragma once

// Contains the spec class
// Handles all operations relating to species lookups

#include <fstream>
#include <map>
#include <string>
#include "common.hpp"
#include "util.hpp"

namespace marzone {
using namespace std;

class Species {
    public:
    Species(sfname& fnames): spno(0), gspno(0), aggexist(false), sepexist(false), fSpecPROPLoaded(false) {
        ReadSpeciesData(fnames.inputdir + fnames.specname, specList, spno, true);

        if (!fnames.blockdefname.empty())
        {
            // read the general species file too
            ReadSpeciesData(fnames.inputdir + fnames.specname, genSpecList, gspno, false);
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
            if (genSpecList[igsp].sepnum > 0)
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].sepnum < 0)
                        specList[isp].sepnum = genSpecList[igsp].sepnum;
            if (genSpecList[igsp].spf > 0)
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].spf < 0)
                        specList[isp].spf = genSpecList[igsp].spf;
            if (genSpecList[igsp].sepdistance > 0)
                for (isp = 0; isp < spno; isp++)
                    if (specList[isp].type == genSpecList[igsp].type && specList[isp].sepdistance < 0)
                        specList[isp].sepdistance = genSpecList[igsp].sepdistance;
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
            if (spec.sepnum < 0)
                spec.sepnum = 0;
            if (spec.sepdistance < 0)
                spec.sepdistance = 0;
            if (spec.spf < 0)
                spec.spf = 1;

            // agg and sep exist
            if (spec.sepdistance)
                sepexist = true;
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

    void LoadCustomPenalties(string filename) {
        FILE *fp = openFile(filename);
        char sLine[500], *sVarVal;

        /* Scan header */
        fgets(sLine,499,fp);

        int i=0, _spid;
        while (fgets(sLine,499,fp))
        {
            sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal,"%d",&_spid);
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal,"%lf",&specList[i].penalty);

            i++;
        }

        fclose(fp);
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

        myfile << "i,name,type,sname,target,prop,targetocc,spf,penalty,amount,occurrence,sepdistance,sepnum,separation,clumps,target2,richness,offset\n";
        for (int i=0;i<spno;i++) {
            myfile << i << "," << specList[i].name << "," << specList[i].type << ","
            << specList[i].sname << "," << specList[i].target << "," << specList[i].prop << "," 
            << specList[i].targetocc << "," << specList[i].spf << "," << specList[i].penalty << "," << reserveStat[i].amount << "," << reserveStat[i].occurrence << "," 
            << specList[i].sepdistance << "," << specList[i].sepnum << "," << reserveStat[i].separation << "," << reserveStat[i].clumps << "," << specList[i].target2 << "," << specList[i].richness << "," << specList[i].offset << "\n";
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
    void ReadSpeciesData(string filename, vector<T>& typeList, unsigned& count, bool populateLookup) {
        FILE *fp = openFile(filename);
        char sLine[1000];

        int numvars = 10,namespecial = 0;
        char *sVarName,*sVarVal;

        /* Scan header */
        fgets(sLine,999,fp);
        vector<string> headerNames = getFileHeaders(sLine, filename);

        /* While there are still lines left populate specList */
        while (fgets(sLine,999,fp))
        {
            T spectemp = {};
            /** Clear important species stats **/
            spectemp.target = -1;
            spectemp.type = -1;
            spectemp.spf = -1;
            spectemp.target2 = 0;
            spectemp.targetocc = -1;
            spectemp.sepdistance = -1;
            spectemp.sepnum = -1;

            for (string temp : headerNames)
            { /* Tok out the next Variable */

                if (namespecial)
                { /* used for special name handling function */
                    namespecial = 0;
                }
                else
                {
                    if (temp == headerNames.front())
                    {
                        sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
                    }
                    else {
                        sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
                    }
                }

                if (temp.compare("id") == 0)
                {
                    sscanf(sVarVal, "%d", &spectemp.name);
                }
                else if (temp.compare("name") == 0)
                {
                    /* Cpy first part of name into this */
                    string speciesname(sVarVal);
                    /* get next part of name */
                    do
                    {
                        sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
                        if (!sVarVal)
                            namespecial = 2;
                        else
                        {
                            if (isalpha(sVarVal[0]))
                            {
                                speciesname += " ";
                                speciesname += string(sVarVal);
                            }
                            else
                                namespecial = 1;
                        }
                    } while (!namespecial);
                    if (namespecial == 2)
                        namespecial = 0; /* Handles end of line case */
                    /* namespecial == 1 means not at end of line and next variable should be processed*/

                    spectemp.sname = speciesname;
                } /* Special name handling routine */
                else if (temp.compare("type") == 0)
                {
                    sscanf(sVarVal, "%d", &spectemp.type);
                }
                else if (temp.compare("target") == 0)
                {
                    sscanf(sVarVal, "%lf", &spectemp.target);
                }
                else if (temp.compare("prop") == 0)
                {
                    sscanf(sVarVal, "%lf", &spectemp.prop);
                    if (spectemp.prop > 0)
                        fSpecPROPLoaded = true;
                }
                else if (temp.compare("fpf") == 0 || temp.compare("spf")) // fpf, spf same header meaning.
                {
                    sscanf(sVarVal, "%lf", &spectemp.spf);
                }
                else if (temp.compare("sepdistance") == 0)
                {
                    sscanf(sVarVal, "%lf", &spectemp.sepdistance);
                }
                else if (temp.compare("sepnum") == 0)
                {
                    sscanf(sVarVal, "%d", &spectemp.sepnum);
                }
                else if (temp.compare("target2") == 0)
                {
                    sscanf(sVarVal, "%lf", &spectemp.target2);
                }
                else if (temp.compare("targetocc") == 0)
                {
                    sscanf(sVarVal, "%d", &spectemp.targetocc);
                }
            } /* looking for ivar different input variables */
            
            typeList.push_back(spectemp);

            if (populateLookup)
                lookup[spectemp.name] = count;

            count++;
        } /* Scanning through each line of file */

        fclose(fp);
    }
};

} // namespace marzone