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
#include "logger.hpp"

namespace marzone {
    
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
    Zones(sfname& fnames, Costs& costs, Logger& logger) :
        zoneCount(0), zoneContribCount(0), zoneContrib2Count(0), zoneContrib3Count(0),
        zoneContribSupplied(false), availableZoneCost(false)
    {
        // Parse all zone files 
        if (!fnames.zonesname.empty())
        {
            LoadZones(fnames.inputdir + fnames.zonesname, logger);
        }
        else
        {
            DefaultZones();
        }

        // Zone conntrib files. Only one should be supplied.
        if (!fnames.zonecontribname.empty())
        {
            zoneContribSupplied = true;
            LoadZoneContrib(fnames.inputdir + fnames.zonecontribname, logger);
        }
        if (!fnames.zonecontrib2name.empty())
        {
            zoneContribSupplied = true;
            LoadZoneContrib2(fnames.inputdir + fnames.zonecontrib2name, logger);
        }
        if (!fnames.zonecontrib3name.empty())
        {
            zoneContribSupplied = true;
            LoadZoneContrib3(fnames.inputdir + fnames.zonecontrib3name, logger);
        }

        // Build zone cost with inputted information
        PopulateZoneCosts(fnames, costs, logger);
        PopulateConnectionCosts(fnames, logger);
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

    void BuildZoneTarget(Species& spec, Pu& pu, sfname& fnames, Logger& logger) {
        int i, j, iSpeciesIndex;

        // Zone target files
        if (!fnames.zonetargetname.empty())
        {
            LoadZoneTarget(fnames.inputdir + fnames.zonetargetname, 1, zoneTargetsTemp, logger);
        }
        if (!fnames.zonetarget2name.empty())
        {
            LoadZoneTarget(fnames.inputdir + fnames.zonetarget2name, 2, zoneTargetsTemp, logger);
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
        int zoneind;
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

            zoneind = LookupIndex(zoneTargetsTemp[i].zoneid);
            if (zoneind == -1)
                continue; // ignore invalid zone ids

            for (int spindex: speciesInd) 
            {
                // .zoneid .speciesid .target
                if (zoneTargetsTemp[i].targettype == 0) // area target as hectare
                    zoneTarget[spindex][zoneind].target = zoneTargetsTemp[i].target;
                if (zoneTargetsTemp[i].targettype == 1) // area target as proportion
                    zoneTarget[spindex][zoneind].target = zoneTargetsTemp[i].target * SpecArea[spindex];
                if (zoneTargetsTemp[i].targettype == 2) // occurrence target as occurrences
                    zoneTarget[spindex][zoneind].occurrence = ceil(zoneTargetsTemp[i].target);
                if (zoneTargetsTemp[i].targettype == 3) // occurrence target as proportion
                    zoneTarget[spindex][zoneind].occurrence = ceil(zoneTargetsTemp[i].target * SpecOcc[spindex]);
            }
        }
    }

    // given a species index and zoneid, return the contrib fraction
    // given a zone index.
    // If zone contrib file was not supplied, all contribs are 1. 
    double GetZoneContrib(int spindex, int zoneind) {
        // only spno*zoneCount used
        return zoneContribValues[(spindex*zoneCount)+zoneind];
    }

    // Get zone contrib but for a specific pu.
    // If zone contrib file was not supplied, all contribs are 1. 
    double GetZoneContrib(int puindex, int puno, int spindex, int zoneind) {
        if (!zoneContrib3Count) {
            // only spno*zoneCount used
            return GetZoneContrib(spindex, zoneind);
        }
        else {
            return zoneContribValues[(spindex*puno*zoneCount)+(puindex*zoneCount)+zoneind];
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
    double ConnectionCost2Linear(Pu& pu, int puindex, int imode, vector<int> &solution, double blm)
    {
        if (pu.connectionsEntered && blm != 0) {
            double fcost, rZoneConnectionCost;
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

            return fcost*blm;
        }

        return 0;
    }

    // Connection cost but we can specify which zone to calculate the current pu for
    double ConnectionCost2(Pu& pu, int puindex, int imode, vector<int>& solution, int curZone, double blm) {
        if (pu.connectionsEntered && blm != 0)
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

            return fcost*blm;
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

    int IndexToId(int index) {
        return zoneNameIndexed[index].id;
    }

    // Given a zoneid, returns index of zone or -1 if invalid.
    int LookupIndex(int zoneid) {
        auto it = zoneNames.find(zoneid);
        if (it != zoneNames.end()) {
            return it->second.index;
        }

        return -1; // zoneid not found.
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

    // usually named debug_Zone<number>_Contrib.csv or debug_ZoneContrib.csv
    void DumpZoneContribFinalValues(string filename, Species& spec) {
        ofstream myfile;

        if (zoneContrib3Count) {
            // TODO - write debug file for each zone
            // each file should contain the matrix pu x spec contribs.

        }
        else {
            // write contrib values normally.
            myfile.open(filename);
            myfile << "spname,spindex\n";
            for (int i=0;i<zoneCount;i++)
                myfile << ",contrib" << IndexToId(i);
            myfile << "\n";

            for (int j=0; j < spec.spno; j++) {
                myfile << spec.specList[j].name << "," << j;
                for (int i =0; i<zoneCount; i++)
                    myfile << "," << GetZoneContrib(j, i);
                myfile << "\n";
            }

            myfile.close();
        }
    }

    void DumpZoneCostFinalValues(string filename, Costs& costs) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "costindex";
        for (int i=0; i < zoneCount; i++)
            myfile << "," << IndexToId(i);
        myfile << "\n";

        for (int j = 0; j < costs.costCount; j++)  
        {
            myfile << "," << j;
            for (int i=0; i< zoneCount; i++)
                myfile << "," << zoneCost[j*zoneCount + i];
            myfile << "\n";
        }
        myfile.close();
    }

    void DumpRelConnectionCostFinalValues(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "zoneindex";

        for(int j=0; j < zoneCount; j++)
            myfile << "," << IndexToId(j);
        myfile << "\n";

        for (int j=0; j<zoneCount; j++) {
            myfile << "," << IndexToId(j);
            for (int i=0; i<zoneCount;i++) {
                myfile << "," << GetZoneConnectionCost(j, i);
            }
            myfile <<"\n";
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

    void DumpZoneTargetFinalValues(string filename, Species& spec) 
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "spname,spindex";
        
        int id;
        for (int i=0; i< zoneCount; i++) {
            id = IndexToId(i);
            myfile << ",zone" << id << "target,zone" << id << "occurrence";
        }
        myfile << "\n";

        for (int j=0; j<spec.spno; j++) {
            myfile << spec.specList[j].name << "," << j;
            for (int i = 0; i < zoneCount; i++) {
                myfile << "," << zoneTarget[j][i].target << "," << zoneTarget[j][i].occurrence;
            }
            myfile << "\n";
        }

        myfile.close();
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

    unsigned zoneCount; // number of available zones in the system.
    vector<ZoneIndex> zoneNameIndexed; // backward map of index to zoneid/zonename
    map<int, ZoneName> zoneNames; // zoneid to zonename/index mapping
    vector<double> zoneContribValues; // size depends on whether zonecontrib3 is used. Usually it is spno*zoneCount. Species, then zone. 

    // Zone target maps
    vector<vector<zonetarget>> zoneTarget; //species -> zone target matrix, ordering is [species][zoneid].

    // Zone cost
    vector<double> zoneCost; // flattened cost -> zone -> fraction map.
    bool availableZoneCost;
    vector<double> zoneConnectionCost; // zone to zone cost matrix

    private:
    vector<zonetargetstructtemp> zoneTargetsTemp; // needed for dump debug.
    bool zoneContribSupplied; // whether zone contribution modifiers are supplied or not.
    vector<zonecontribstruct> zoneContrib;
    vector<zonecontrib2struct> zoneContrib2;
    vector<zonecontrib3struct> zoneContrib3;
    unsigned zoneContribCount;
    unsigned zoneContrib2Count;
    unsigned zoneContrib3Count;
    vector<zonecoststruct> zoneCostFileLoad; // temp zonecost from file
    vector<relconnectioncoststruct> zoneRelConnectionCost; // temp raw figures.

    // connection cost between two zone indices
    double GetZoneConnectionCost(int zoneindex1, int zoneindex2) {
        return zoneConnectionCost[zoneindex1*zoneCount + zoneindex2];
    }

    void PopulateConnectionCosts(sfname& fnames, Logger& logger) {
        if (!fnames.relconnectioncostname.empty())
        {
            zoneRelConnectionCost = LoadRelConnectionCost(fnames.inputdir + fnames.relconnectioncostname, logger);
        }
        else {
            zoneConnectionCost.assign(zoneCount*zoneCount, 1); // default zone boundary cost is 1 (to not affect pu costs)
            return;
        }

        zoneConnectionCost.assign(zoneCount*zoneCount, 0);

        int zoneind1, zoneind2;
        for (relconnectioncoststruct& term: zoneRelConnectionCost) {
            zoneind1 = zoneNames[term.zoneid1].index;
            zoneind2 = zoneNames[term.zoneid2].index;
            zoneConnectionCost[zoneind1*zoneCount+zoneind2] = term.fraction;
            zoneConnectionCost[zoneind2*zoneCount+zoneind1] = term.fraction;
        }
    }

    void PopulateZoneCosts(sfname& fnames, Costs& costs, Logger& logger) {
        // read zone cost files (if any)
        if (!fnames.zonecostname.empty())
        {
            zoneCostFileLoad = LoadZoneCost(fnames.inputdir + fnames.zonecostname, logger);
            availableZoneCost = true;
        }
        else {
            zoneCost.assign(costs.costCount*zoneCount, 1.0);
            return;
        }

        zoneCost.assign(costs.costCount*zoneCount, 0.0);
        int zoneind1;
        for (zonecoststruct& cost: zoneCostFileLoad) {
            zoneind1 = zoneNames[cost.zoneid].index;
            zoneCost[costs.GetCostIndex(cost.costid)*zoneCount+zoneind1] = cost.fraction;
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
        int zoneind;
        if (zoneContribCount) {
            for (zonecontribstruct& contrib: zoneContrib) {
                zoneind = LookupIndex(contrib.zoneid);
                spindex = spec.LookupIndex(contrib.speciesid); 
                if (spindex != -1 && zoneind != -1) 
                    zoneContribValues[(spindex*zoneCount)+zoneind] = contrib.fraction;
            }
        }
        else if (zoneContrib2Count) {
            // Contribs in this file apply to all species
            for (zonecontrib2struct& contrib : zoneContrib2) {
                zoneind = LookupIndex(contrib.zoneid);
                if (zoneind != -1)
                    for (int i = 0; i < spec.spno; i++) {
                        zoneContribValues[(i*zoneCount)+zoneind] = contrib.fraction;
                    }
            }
        }
        else {
            for (zonecontrib3struct& contrib : zoneContrib3)
            {
                puindex = pu.LookupIndex(contrib.puid);
                spindex = spec.LookupIndex(contrib.speciesid); 
                zoneind = LookupIndex(contrib.zoneid);
                if (spindex != -1 && puindex != -1 && zoneind != -1)
                    zoneContribValues[(spindex*pu.puno*zoneCount)+(puindex*zoneCount)+zoneind] = contrib.fraction;
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
                if (i == 0) // the "available" zone.
                {
                    zoneContribValues[(j*zoneCount)+i] = 0;
                } else {
                    zoneContribValues[(j*zoneCount)+i] = 1;
                }
            }
        }
    }

    vector<relconnectioncoststruct> LoadRelConnectionCost(string filename, Logger& logger) {
        ifstream fp = openFile(filename);
        string sLine;
        int zoneid1, zoneid2;
        double fraction;
        vector<relconnectioncoststruct> zoneRelConnectionCost;

        // create the RelConnectionCost array
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

            // read the integer zoneid1 from this line
            stringstream ss = stream_line(sLine);
            ss >> zoneid1;
            // read the integer zoneid2 from this line
            ss >> zoneid2;
            // read the double fraction from this line
            ss >> fraction;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            zoneRelConnectionCost.push_back({zoneid1, zoneid2, fraction});
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");

        return zoneRelConnectionCost;
    }

    vector<zonecoststruct> LoadZoneCost(string filename, Logger& logger)
    {
        vector<zonecoststruct> tempZoneCost;
        ifstream fp = openFile(filename);
        string sLine;
        int zoneid, costid;
        double fraction;

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
            // read the integer zoneid from this line
            stringstream ss = stream_line(sLine);
            ss >> zoneid;
            // read the integer costid from this line
            ss >> costid;
            // read the double fraction from this line
            ss >> fraction;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");
            tempZoneCost.push_back({zoneid, costid, fraction});
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");

        return tempZoneCost;
    }

    void LoadZoneContrib(string filename, Logger& logger) 
    {
        ifstream fp = openFile(filename);
        string sLine;
        int zoneid, speciesid;
        double fraction;

        // create the ZoneContrib array
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
            // read the integer zoneid from this line
            stringstream ss = stream_line(sLine);
            ss >> zoneid;
            // read the integer speciesid from this line
            ss >> speciesid;
            // read the double fraction from this line
            ss >> fraction;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            zoneContrib.push_back({zoneid, speciesid, fraction});
            zoneContribCount++;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");
    }

    // This file should just be num zones sized. 
    void LoadZoneContrib2(string filename, Logger& logger)
    {
        ifstream fp = openFile(filename);
        string sLine;
        int zoneid;
        double fraction;

        // create the ZoneContrib array
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
            // read the integer zoneid from this line
            stringstream ss = stream_line(sLine);
            ss >> zoneid;

            // read the double fraction from this line
            ss >> fraction;

            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            zoneContrib2.push_back({zoneid, fraction});
            zoneContrib2Count++;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");
    }

    void LoadZoneContrib3(string filename, Logger& logger)
    {
        ifstream fp = openFile(filename);
        string sLine;
        int zoneid, puid, speciesid;
        double fraction;

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
            // read the integer zoneid from this line
            stringstream ss = stream_line(sLine);
            ss >> zoneid;

            // read the integer puid from this line
            ss >> puid;

            // read the integer speciesid from this line
            ss >> speciesid;

            // read the double fraction from this line
            ss >> fraction;

            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");


            zoneContrib3.push_back({zoneid, puid, speciesid, fraction});
            zoneContrib3Count++;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");
    }

    // type = determines which kind of target is being loaded. targetMode=1 means regular target. Else target2.
    void LoadZoneTarget(string filename, int targetMode, vector<zonetargetstructtemp>& zoneTargets, Logger& logger)
    {
        ifstream fp = openFile(filename);
        string sLine;
        int zoneid, speciesid = -1, targettype;
        double target;

        // skip header
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
            // read the integer zoneid from this line
            stringstream ss = stream_line(sLine);
            ss >> zoneid;
            // read the integer speciesid from this line if type1
            if (targetMode == 1) {
                ss >> speciesid;
            }

            // read the double fraction from this line
            ss >> target;

            // read the integer targettype from this line if it exists
            if (!ss)
                targettype = 0;
            else
                ss >> targettype;
           
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");
            zoneTargets.push_back({zoneid, speciesid, target, targettype});
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;
    }

    void LoadZones(string filename, Logger& logger) {
        ifstream fp = openFile(filename);
        string sLine;
        int tempId;

        // load the data to a zones map
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
            // read the integer id from this line
            stringstream ss = stream_line(sLine);
            ss >> tempId;

            // read the string name from this line
            string sval;
            ss >> sval;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");
            trim(sval);
            zoneNames[tempId] = {sval, zoneCount};

            // Construct index object
            ZoneIndex z;
            z.id = tempId;
            z.name = sval;
            zoneNameIndexed.push_back(z);
            
            zoneCount++;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");
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
};

} // namespace marzone