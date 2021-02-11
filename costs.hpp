#pragma once

#include <map>

#include "common.hpp"
#include "util.hpp"
#include "logger.hpp"
/*
 Short class that handles cost information.
*/
namespace marzone {

extern Logger logger;

typedef struct costField
{
    int id; // original cost id
    unsigned index; // zero-indexed cost id. Costs are indexed based on the order they appear.
} costField;

class Costs {
    public:
    Costs(sfname& fnames) : costCount(0) {
        if (!fnames.costsname.empty())
        {
            LoadCostNames(fnames.inputdir + fnames.costsname);
        }
        else
        {
            DefaultCostNames();
        }
    }

    // Prerequisite - costName must be defined in original cost file
    unsigned GetCostIndex(string costName) {
        return costNames[costName].index;
    }

    unsigned GetCostIndex(int costId) {
        return costIds[costId];
    }

    // Checks if costName is defined
    bool Contains(string costName) {
        if (costNames.find(costName) == costNames.end()) {
            return false;
        }
        return true;
    }

    // Debugging functions (for dumping)
    void DumpCostNames(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "costid,costname\n";
        for (auto& [name, term]: costNames)
        {
            myfile << term.id << "," << name << "\n";
        }
        myfile.close();
    }

    unsigned costCount;

    private:
    void LoadCostNames(string filename)
    {
        ifstream fp = openFile(filename);
        string sLine;
        int id;
        string name;

        // create the CostNames array
        // load the data to an array
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
            ss >> id;
            string costName;
            ss >> costName;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");
            trim(costName);
            costNames[costName] = {id, costCount};
            costIds[id] = costCount;
            costCount++;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;
    }

    void DefaultCostNames()
    {
        costCount= 1;
        costIds[1] = 0u;
        costNames["cost"] = {1, 0u};
    }

    map<string, costField> costNames;
    map<int, unsigned> costIds; // cost id to index.
};

} // namespace marzone