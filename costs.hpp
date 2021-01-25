#pragma once
#include "common.hpp"
#include "util.hpp"

/*
 Short class that handles cost information.
*/
namespace marzone {

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
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int id;
        string name;

        // ignore header
        fgets(sLine,999,fp);

        // create the CostNames array
        // load the data to an array
        while (fgets(sLine,999,fp))
        {
            sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &id);

            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            costNames[string(sVarVal)] = {id, costCount};
            costCount++;
        }
        fclose(fp);
    }

    void DefaultCostNames()
    {
        costCount= 1;
        costNames["cost"] = {1, 0u};
    }

    map<string, costField> costNames;
};

} // namespace marzone