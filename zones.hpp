#pragma once
/*
  Class that handles operations relating to zones. Includes parsing.
*/

#include <stdexcept>
#include <map>
#include "common.hpp"
#include "pu.hpp"
#include "species.hpp"
#include "util.hpp"

namespace marzone {

typedef struct ZoneName {
    string name;
    uint64_t index;
} ZoneName;

typedef struct ZoneIndex : ZoneName {
    int id;
} ZoneIndex;

typedef struct zonecoststruct
{
    int zoneid;
    int costid;
    double fraction;
} zonecoststruct;

typedef struct relconnectioncoststruct
{
    int zoneid1;
    int zoneid2;
    double fraction;
} relconnectioncoststruct;

typedef struct zonetargetstructtemp
{
    int zoneid;
    int speciesid;
    double target;
    int targettype;
} zonetargetstructtemp;

typedef struct zonetarget
{
    double target;
    int occurrence;
} zonetarget;

class Zones {
    public:
    Zones(sfname& fnames, Costs& costs) :
        zoneCount(0), zoneContribCount(0), zoneContrib2Count(0), zoneContrib3Count(0), availableZoneInput(false)
    {
        // Parse all zone files 
        if (!fnames.zonesname.empty())
        {
            LoadZones(fnames.inputdir + fnames.zonesname);
        }
        else
        {
            DefaultZones();
            availableZoneInput = true; // double check true/false and use cases?
        }

        // Zone conntrib files. Only one should be supplied.
        if (!fnames.zonecontribname.empty())
        {
            LoadZoneContrib(fnames.inputdir + fnames.zonecontribname);
        }
        if (!fnames.zonecontrib2name.empty())
        {
            LoadZoneContrib2(fnames.inputdir + fnames.zonecontrib2name);
        }
        if (!fnames.zonecontrib3name.empty())
        {
            LoadZoneContrib3(fnames.inputdir + fnames.zonecontrib3name);
        }

        // Build zone cost with inputted information
        PopulateZoneCosts(fnames, costs.costCount);
        PopulateConnectionCosts(fnames);
    }

    // Constructs the zone contributions map depending on which type of zone contribution was used. 
    void BuildZoneContributions(Species& spec, Pu& pu) {
        if (!zoneContribCount && !zoneContrib2Count && !zoneContrib3Count) {
            DefaultZoneContributions(spec);
        }
        else {
            PopulateZoneContributions(spec, pu);
        }
    }

    void BuildZoneTarget(Species& spec, Pu& pu, sfname& fnames) {
        int i, j, iSpeciesIndex;

        // Zone target files
        if (!fnames.zonetargetname.empty())
        {
            LoadZoneTarget(fnames.inputdir + fnames.zonetargetname, 1, zoneTargetsTemp);
        }
        if (!fnames.zonetarget2name.empty())
        {
            LoadZoneTarget(fnames.inputdir + fnames.zonetarget2name, 2, zoneTargetsTemp);
        }

        if (zoneTargetsTemp.empty())
            return; // do nothing if target file not supplied.

        // init arrays of species area and occurrence totals
        vector<double> SpecArea = pu.TotalSpeciesAmount(spec.spno);
        vector<int> SpecOcc = pu.TotalOccurrenceAmount(spec.spno);

        // create and initialise _ZoneTarget
        zoneTarget.resize(spec.spno);
        for (auto& zone: zoneTarget) {
            zone.assign(zoneCount, {0,0});
        }

        // populate _ZoneTarget from ZoneTarget
        vector<int> speciesInd;
        for (i = 0; i < zoneTargetsTemp.size(); i++)
        {
            // note that speciesid could be -1, which means it applies to all species.
            if (zoneTargetsTemp[i].speciesid != -1)
            {
                iSpeciesIndex = spec.LookupIndex(zoneTargetsTemp[i].speciesid);
                if (iSpeciesIndex == -1) // invalid species - skip
                    continue;
                speciesInd.clear();
                speciesInd.push_back(iSpeciesIndex);
            }
            else {
                speciesInd = Range(0, spec.spno);
            }

            for (int spindex: speciesInd) 
            {
                // .zoneid .speciesid .target
                if (zoneTargetsTemp[i].targettype == 0) // area target as hectare
                    zoneTarget[spindex][zoneTargetsTemp[i].zoneid - 1].target = zoneTargetsTemp[i].target;
                if (zoneTargetsTemp[i].targettype == 1) // area target as proportion
                    zoneTarget[spindex][zoneTargetsTemp[i].zoneid - 1].target = zoneTargetsTemp[i].target * SpecArea[spindex];
                if (zoneTargetsTemp[i].targettype == 2) // occurrence target as occurrences
                    zoneTarget[spindex][zoneTargetsTemp[i].zoneid - 1].occurrence = ceil(zoneTargetsTemp[i].target);
                if (zoneTargetsTemp[i].targettype == 3) // occurrence target as proportion
                    zoneTarget[spindex][zoneTargetsTemp[i].zoneid - 1].occurrence = ceil(zoneTargetsTemp[i].target * SpecOcc[spindex]);
            }
        }
    }

    // given a species index and zoneid, return the contrib fraction
    double GetZoneContrib(int spindex, int zoneid) {
        // only spno*zoneCount used
        return zoneContribValues[(spindex*zoneCount)+(zoneid-1)];
    }

    // Get zone contrib but for a specific pu.
    double GetZoneContrib(int puindex, int puno, int spindex, int zoneid) {
        if (!zoneContrib3Count) {
            // only spno*zoneCount used
            return GetZoneContrib(spindex, zoneid);
        }
        else {
            return zoneContribValues[(spindex*puno*zoneCount)+(puindex*zoneCount)+(zoneid-1)];
        }
    }

    // returns an aggregate vector by species containing total zone targets and occurences.
    vector<double> AggregateTargetAreaBySpecies(int spno) {
        // If zone target is not specified, return a vector of zeroes
        if (zoneTarget.size() == 0)
            return vector<double>(spno, 0.0);
        
        vector<double> toReturn(spno, 0.0);

        for (int i = 0; i < spno; i++) {
            for (zonetarget& zoneTerm: zoneTarget[i]) {
                toReturn[i] += zoneTerm.target;
            }
        }

        return toReturn;
    }

    vector<int> AggregateTargetOccurrenceBySpecies(int spno) {
        // If zone target is not specified, return a vector of zeroes
        if (zoneTarget.size() == 0)
            return vector<int>(spno, 0);
        
        vector<int> toReturn(spno, 0);

        for (int i = 0; i < spno; i++) {
            for (zonetarget& zoneTerm: zoneTarget[i]) {
                toReturn[i] += zoneTerm.occurrence;
            }
        }

        return toReturn;
    }

    // Given a pu cost vector, returns the total cost for that vector.
    double AggregateTotalCostByPuAndZone(int zoneIndex, vector<double>& puCost) {
        double rCost = 0;

        for (int i=0;i<puCost.size();i++)
        {
            rCost += puCost[i] * zoneCost[(i*zoneCount)+(zoneIndex)];
        }

        return rCost;
    }

    // Returns the connection cost of a puindex, given the zones of the other pu.
    // imode = check specifics of this param. For now I am assuming 1, 0 or -1
    double ConnectionCost2Linear(Pu& pu, int puindex, int imode, vector<int> &solution)
    {
        if (pu.connectionsEntered) {
            double fcost, rResult, rZoneConnectionCost;
            int iCurrentZone = solution[puindex];

            fcost = pu.connections[puindex].fixedcost * imode;
            for (sneighbour &p : pu.connections[puindex].first)
            {
                if (p.nbr > puindex)
                {
                    rZoneConnectionCost = GetZoneConnectionCost(iCurrentZone, solution[p.nbr]);
                    fcost += imode * p.cost * rZoneConnectionCost;
                }
            }

            return fcost;
        }

        return 0;
    }

    // Connection cost but we can specify which zone to calculate the current pu for
    double ConnectionCost2(Pu& pu, int puindex, int imode, vector<int>& solution, int curZone) {
        if (pu.connectionsEntered)
        {
            double fcost, rZoneConnectionCost;
            // Initial fixed cost
            fcost = pu.connections[puindex].fixedcost * imode;

            // Asymmetric connectivity not supported in marzone, so we can ignore it.
            // We can add it back in the future if needed.
            for (sneighbour &p : pu.connections[puindex].first)
            {
                rZoneConnectionCost = GetZoneConnectionCost(curZone, solution[p.nbr]);
                fcost += imode * p.cost * rZoneConnectionCost;
            }

            return fcost;
        }
        return 0;
    }

    vector<vector<double>> InitializeZoneMatrix() {
        // Returns an empty matrix of size zoneCount x zoneCount
        vector<vector<double>> matrix(zoneCount);
        for (auto& row: matrix) {
            row.resize(zoneCount, 0);
        }
        return matrix;
    }

    /*
        Getter functions
    */
    string IndexToName(int index) {
        return zoneNameIndexed[index].name;
    }

    // Formats all zone names into a string
    string ZoneNameHeaders(int imode, string othertext) {
        string d = imode > 1 ? "," : "    ";
        stringstream text;

        for (int i = 0; i < zoneCount; i++)
        {
            text << d;
            if (imode == 2)
                text << "\"";
            text << zoneNameIndexed[i].name << othertext;
            if (imode == 2)
                text << "\"";
        }
        return text.str();
    }

    void WriteZoneTargetHeaders(string filename, int imode) {
        ofstream myfile;
        myfile.open(filename);

        if (imode > 1)
        {
            myfile << "\"Feature\",\"Feature Name\",\"Target\"";
            myfile << ",\"Total Amount\",\"Contributing Amount Held\",\"Occurrence Target \",\"Occurrences Held\",\"Target Met\"";
            for (int i=0;i<zoneCount;i++)
            {
                string zoneName = zoneNameIndexed[i].name;
                myfile << ",\"Target " << zoneName << "\",\"Amount Held " << zoneName << "\",\"Contributing Amount Held " << zoneName 
                << "\",\"Occurrence Target " << zoneName << "\"" << ",\"Occurrences Held " << zoneName << "\",\"Target Met " << zoneName << "\"";
            }
            myfile << ",MPM\n";
        }
        else {
            myfile << "Feature\tFeature Name\tTarget";
            myfile << "\tAmount Held\tContributing Amount Held\tOccurrence Target \tOccurrences Held\tTarget Met";
            for (int i=0;i<zoneCount;i++)
            {
                string zoneName = zoneNameIndexed[i].name;
                myfile << "\tTarget " << zoneName << "\tAmount Held " << zoneName << "\tContributing Amount Held " << zoneName 
                << "\tOccurrence Target " << zoneName << "\tOccurrences Held " << zoneName << "\tTarget Met " << zoneName;
            }
            myfile << "\tMPM\n";
        }

        myfile.close();
    }

    /*
      Debugging functions (for dumping)
    */
    void DumpZoneNames(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,zonename\n";
        for (int i = 0; i < zoneCount; i++)
        {
            myfile << zoneNameIndexed[i].id << "," << zoneNameIndexed[i].name << "\n";
        }
        myfile.close();
    }

    void DumpZoneContrib(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,speciesid,fraction\n";
        for (zonecontribstruct& z: zoneContrib)
        {
            myfile << z.zoneid << "," << z.speciesid << "," << z.fraction << "\n";
        }
        myfile.close();
    }

    void DumpZoneContrib2(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,fraction\n";
        for (zonecontrib2struct& z: zoneContrib2)
        {
            myfile << z.zoneid  << "," << z.fraction << "\n";
        }
        myfile.close();
    }

    void DumpZoneContrib3(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,puid,speciesid,fraction\n";
        for (zonecontrib3struct& z: zoneContrib3)
        {
            myfile << z.zoneid<<","<<z.puid<<","<< z.speciesid<<","<< z.fraction << "\n";
        }
        myfile.close();
    }

    // dumps zone target with a species id
    void DumpZoneTarget(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,speciesid,target,targettype\n";
        for (zonetargetstructtemp& z: zoneTargetsTemp) {
            if (z.speciesid)
                myfile<<z.zoneid<<","<<z.speciesid<<","<<z.target<<","<<z.targettype<<"\n";
        }
    }

    // dumps zone target WITHOUT a species id (zonetarget2)
    void DumpZoneTarget2(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,target,targettype\n";
        for (zonetargetstructtemp& z: zoneTargetsTemp) {
            if (z.speciesid == 0)
            myfile<<z.zoneid<<","<<z.target<<","<<z.targettype<<"\n";
        }
    }

    void DumpZoneCost(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid,costid,fraction\n";
        for (zonecoststruct& z: zoneCostFileLoad)
        {
            myfile <<z.zoneid<<","<<z.costid<<","<<z.fraction<<"\n";
        }
    }

    void DumpRelConnectionCost(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneid1,zoneid2,fraction\n";
        for (auto& r: zoneRelConnectionCost)
        {
            myfile <<r.zoneid1<<","<<r.zoneid2<<","<< r.fraction<<"\n";
        }
        myfile.close();
    }

    bool availableZoneInput;
    unsigned zoneCount; // number of available zones in the system.
    unsigned zoneContribCount;
    unsigned zoneContrib2Count;
    unsigned zoneContrib3Count;
    unsigned zoneTarget2Count;
    vector<ZoneIndex> zoneNameIndexed; // backward map of index to zoneid/zonename
    map<int, ZoneName> zoneNames; // zoneid to zonename/index mapping

    // Zone contribution maps. TODO - merge into 1
    vector<zonecontribstruct> zoneContrib;
    vector<zonecontrib2struct> zoneContrib2;
    vector<zonecontrib3struct> zoneContrib3;
    vector<double> zoneContribValues; // size depends on whether zonecontrib3 is used. Usually it is spno*zoneCount. Species, then zone. 

    // Zone target maps
    vector<zonecoststruct> zoneCostFileLoad; // temp zonecost from file
    vector<zonetargetstructtemp> zoneTargetsTemp; // needed for dump debug.
    vector<vector<zonetarget>> zoneTarget; //species -> zone target matrix, ordering is [species][zoneid].

    // Zone cost
    vector<relconnectioncoststruct> zoneRelConnectionCost; // temp raw figures.
    vector<double> zoneCost; // flattened cost -> zone -> fraction map.
    vector<double> zoneConnectionCost; // zone to zone cost matrix

    private:
    // connection cost between two zone indices
    double GetZoneConnectionCost(int zoneindex1, int zoneindex2) {
        return zoneConnectionCost[zoneindex1*zoneCount + zoneindex2];
    }

    void PopulateConnectionCosts(sfname& fnames) {
        if (!fnames.relconnectioncostname.empty())
        {
            zoneRelConnectionCost = LoadRelConnectionCost(fnames.inputdir + fnames.relconnectioncostname);
        }

        zoneConnectionCost.assign(zoneCount*zoneCount, 0); // default zone boundary cost is 0

        for (relconnectioncoststruct& term: zoneRelConnectionCost) {
            zoneConnectionCost[((term.zoneid1-1)*zoneCount)+(term.zoneid2-1)] = term.fraction;
            zoneConnectionCost[((term.zoneid2-1)*zoneCount)+(term.zoneid1-1)] = term.fraction;
        }
    }

    void PopulateZoneCosts(sfname& fnames, int costCount) {
        // read zone cost files (if any)
        if (!fnames.zonecostname.empty())
        {
            zoneCostFileLoad = LoadZoneCost(fnames.inputdir + fnames.zonecostname);
        }
        else
        {
            zoneCostFileLoad = DefaultZoneCost();
        }

        zoneCost.assign(costCount*zoneCount, 0.0);

        for (zonecoststruct& cost: zoneCostFileLoad) {
            zoneCost[((cost.costid-1)*zoneCount)+(cost.zoneid-1)] = cost.fraction;
        }
    }

    void PopulateZoneContributions(Species& spec, Pu& pu) {
        // Here we populate differently depending on which contrib was used. 
        // Prioritize contrib1 >> contrib2 >> contrib3
        int puindex, spindex;

        uint64_t arraySize = 0;
        if (zoneContribCount || zoneContrib2Count) {
            arraySize = spec.spno * zoneCount;
        }
        else {
            arraySize = pu.puno * spec.spno * zoneCount; // note this could explode in size.
        }

        zoneContribValues.assign(arraySize, 0); // Default to 0 contrib if not specified.

        if (zoneContribCount) {
            for (zonecontribstruct& contrib: zoneContrib) {
                spindex = spec.LookupIndex(contrib.speciesid); 
                if (spindex != -1) // we ignore undefined species. TODO - should we check for valid zone too?
                    zoneContribValues[(spindex*zoneCount)+(contrib.zoneid-1)] = contrib.fraction;
            }
        }
        else if (zoneContrib2Count) {
            // Contribs in this file apply to all species
            for (zonecontrib2struct& contrib : zoneContrib2) {
                for (int i = 0; i < spec.spno; i++) {
                    zoneContribValues[(i*zoneCount)+(contrib.zoneid-1)] = contrib.fraction;
                }
            }
        }
        else {
            for (zonecontrib3struct& contrib : zoneContrib3)
            {
                puindex = pu.LookupIndex(contrib.puid);
                spindex = spec.LookupIndex(contrib.speciesid); 
                if (spindex != -1 && puindex != -1)
                    zoneContribValues[(spindex*pu.puno*zoneCount)+(puindex*zoneCount)+(contrib.zoneid-1)] = contrib.fraction;
            }
        }
    }

    void DefaultZoneContributions(Species& spec) {
        // neither zonecontrib.dat or zonecontrib2.dat exist so we are using defaults of 1 for each zone and species
        zoneContribValues.resize(spec.spno * zoneCount);
        for (int j=0;j<spec.spno;j++)
        {
            for (int i=0; i<zoneCount; i++)
            {
                if (availableZoneInput == (i+1))
                {
                    zoneContribValues[(j*zoneCount)+i] = 0;
                } else {
                    zoneContribValues[(j*zoneCount)+i] = 1;
                }
            }
        }
    }

    vector<relconnectioncoststruct> LoadRelConnectionCost(string filename) {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int zoneid1, zoneid2;
        double fraction;
        vector<relconnectioncoststruct> zoneRelConnectionCost;

        // count the number of records
        fgets(sLine,999,fp);

        // create the RelConnectionCost array
        while (fgets(sLine,999,fp))
        {
            // read the integer zoneid1 from this line
            sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid1);

            // read the integer zoneid2 from this line
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid2);

            // read the double fraction from this line
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &fraction);

            zoneRelConnectionCost.push_back({zoneid1, zoneid2, fraction});
        }
        fclose(fp);

        return zoneRelConnectionCost;
    }

    vector<zonecoststruct> LoadZoneCost(string filename)
    {
        vector<zonecoststruct> tempZoneCost;
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int zoneid, costid; 
        double fraction;

        fgets(sLine,999,fp);

        // load the data to an array
        while (fgets(sLine,999,fp))
        {
            // read the integer zoneid from this line
            sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // read the integer costid from this line
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &costid);

            // read the double fraction from this line
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &fraction);

            tempZoneCost.push_back({zoneid, costid, fraction});
        }
        fclose(fp);

        return tempZoneCost;
    }

    void LoadZoneContrib(string filename) 
    {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int zoneid, speciesid;
        double fraction;

        // skip header
        fgets(sLine, 999, fp);

        // create the ZoneContrib array
        while (fgets(sLine, 999, fp))
        {
            // read the integer zoneid from this line
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // read the integer speciesid from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &speciesid);

            // read the double fraction from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &fraction);

            zoneContrib.push_back({zoneid, speciesid, fraction});
            zoneContribCount++;
        }
        fclose(fp);
    }

    // This file should just be num zones sized. 
    void LoadZoneContrib2(string filename)
    {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int zoneid;
        double fraction;

        // skip header
        fgets(sLine, 999, fp);

        // create the ZoneContrib array
        // load the data to an array
        while (fgets(sLine, 999, fp))
        {
            // read the integer zoneid from this line
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // read the double fraction from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &fraction);

            zoneContrib2.push_back({zoneid, fraction});
            zoneContrib2Count++;
        }
        fclose(fp);
    }

    void LoadZoneContrib3(string filename)
    {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int zoneid, puid, speciesid;
        double fraction;

        // skip header
        fgets(sLine, 999, fp);

        // load the data to an array
        while (fgets(sLine, 999, fp))
        {
            // read the integer zoneid from this line
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // read the integer puid from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &puid);

            // read the integer speciesid from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &speciesid);

            // read the double fraction from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &fraction);

            zoneContrib3.push_back({zoneid, puid, speciesid, fraction});
            zoneContrib3Count++;
        }
        fclose(fp);
    }

    // type = determines which kind of target is being loaded. targetMode=1 means regular target. Else target2.
    void LoadZoneTarget(string filename, int targetMode, vector<zonetargetstructtemp>& zoneTargets)
    {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int zoneid, speciesid = -1, targettype;
        double target;

        // skip header
        fgets(sLine,999,fp);
        while (fgets(sLine,999,fp))
        {
            // read the integer zoneid from this line
            sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // read the integer speciesid from this line if type1
            if (targetMode == 1) {
                sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
                sscanf(sVarVal, "%d", &speciesid);
            }

            // read the double fraction from this line
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &target);

            // read the integer targettype from this line if it exists
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            if (sVarVal == NULL)
                targettype = 0;
            else
                sscanf(sVarVal, "%d", &targettype);

            zoneTargets.push_back({zoneid, speciesid, target, targettype});
        }
        fclose(fp);
    }

    void LoadZones(string filename) {
        FILE *fp = openFile(filename);
        char sLine[5000], *sVarVal;
        int tempId;

        // skip header
        fgets(sLine,4999,fp);

        // load the data to a zones map
        while (fgets(sLine,4999,fp))
        {
            // read the integer id from this line
            sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &tempId);

            // read the string name from this line
            sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            zoneNames[tempId] = {string(sVarVal), zoneCount};

            // Construct index object
            ZoneIndex z;
            z.id = tempId;
            z.name = string(sVarVal);
            zoneNameIndexed.push_back(z);
            
            zoneCount++;
        }
        fclose(fp);
    }

    // Non supplied zones - create the default zones id to name map
    void DefaultZones() {
        zoneCount = 2;
        zoneNames[1] = {"available", 0};
        zoneNames[2] = {"reserved", 1};

        zoneNameIndexed.resize(2);
        zoneNameIndexed[0].id = 1;
        zoneNameIndexed[0].name = "available";

        zoneNameIndexed[1].id = 2;
        zoneNameIndexed[1].name = "reserved";
    }

    // Non supplied zone costs - create defaults
    vector<zonecoststruct> DefaultZoneCost()
    {
        vector<zonecoststruct> tempZoneCost;
        // create the ZoneCost array
        zonecoststruct z = {2, 1, 1.0}; // zoneid, costid, fraction.
        tempZoneCost.push_back(z);
        return tempZoneCost;
    }
};

} // namespace marzone