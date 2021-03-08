#pragma once

// Contains the pu class
// Handles all operations relating to pu, puzone, pulock etc. and pu connections.
// Also includes the sparse matrix. 

#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <fstream>
#include "costs.hpp"
#include "common.hpp"
#include "species.hpp"
#include "util.hpp"
#include "logger.hpp"

namespace marzone {
using namespace std;

class Pu {
    public:
    Pu(sfname& fnames, Costs& costs, int asymmetric, map<int, ZoneName>& zoneLookup, Logger& logger) : puno(0), puLockCount(0), puZoneCount(0) {
        ReadPUData(fnames.inputdir + fnames.puname, costs, logger);

        if (!fnames.pulockname.empty())
        {
            LoadPuLock(fnames.inputdir + fnames.pulockname, logger);
        }

       // Persist pulock and puzone counts
       PersistPuLock(zoneLookup);

        // Read and persist puzone.
        if (!fnames.puzonename.empty())
        {
            puZoneTemp = LoadPuZone(fnames.inputdir + fnames.puzonename, logger);
            PersistPuZone(puZoneTemp, zoneLookup, logger);
        }

       // Read and set connections
        if (!fnames.connectionname.empty())
        {
            connectionsEntered = true;
            ReadConnections(fnames.inputdir + fnames.connectionname, asymmetric, logger);
        }
        else {
            connectionsEntered = false;
        }

        // construct valid pu list
        ConstructValidPuList();

        asymmetric = asymmetric;
    }

    // given a planning unit id, returns index of pu. Returns -1 if not originally defined.
    // This operation is O(log) by number of pu.
    int LookupIndex(int puid) {
        auto it = lookup.find(puid);
        if (it != lookup.end()) {
            return it->second;
        }

        return -1; // pu not found.
    }

    void LoadSparseMatrix(Species& spec, string filename, Logger& logger) {
        vector<map<int,spu>> SMTemp(puno); // temporarily storing in this structure prevents the need for ordering.
        string sLine;
        int _puid, _spid, puindex, spindex;
        double amount;

        ifstream fp;
        fp.open(filename);
        if (!fp.is_open())
            logger.ShowErrorMessage("Cannot open file "+ filename + ".\n");

        bool file_is_empty = true;
        for (int line_num = 1; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (line_num == 1)
            {
                if (is_like_numerical_data(sLine))
                    logger.ShowWarningMessage("File " + filename + " has no header in the first line.\n");
                else
                    continue;
            }
            if (sLine.empty())
                continue;

            stringstream ss = stream_line(sLine);
            ss >> _spid >> _puid >> amount;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) +".\n");


            spindex = spec.LookupIndex(_spid);
            if (spindex == -1)
            {
                logger.ShowWarningMessage(filename + " found undefined spid " + to_string(_spid) + ". Ignoring this species.\n");
                continue;
            }

            puindex = LookupIndex(_puid);
            if (puindex == -1)
            {
                logger.ShowWarningMessage(filename + " found undefined puid " + to_string(_puid) + ". Ignoring this pu.\n");
                continue;
            }

            /* increment richness for planning unit containing this feature */
            puList[puindex].richness += 1;

            spu temp;
            temp.amount = amount;
            temp.clump = 0;
            SMTemp[puindex][spindex] = temp;
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;

        // load all into real spu vector
        int j = 0;
        for (int i = 0; i < puno; i++) {
            puList[i].offset = j;
            for (auto &[spindex, val] : SMTemp[i]) {
                spu temp = val;
                temp.spindex = spindex;
                puvspr.push_back(temp);
                j++;
            }
        }

        density = puvspr.size()*100.0/(puno*spec.spno);
    }

    // Returns a vector containing the summation of all species amounts in all the pu
    vector<double> TotalSpeciesAmount(int spno) {
        vector<double> specSums(spno, 0.0);
        for (spu& term: puvspr) {
            if (term.spindex >= 0 && term.spindex < spno) { // ensure species within given bound
                specSums[term.spindex] += term.amount;
            }
        }

        return specSums;
    }

    vector<int> TotalOccurrenceAmount(int spno) {
        vector<int> specOcc(spno, 0);
        for (spu& term: puvspr) {
            if (term.spindex >= 0 && term.spindex < spno) { // ensure species within given bound
                specOcc[term.spindex] ++;
            }
        }

        return specOcc;
    }

    // Total species amount vector based on status and pu for pu that do not have zone limitations
    // puStatusZone is a given reserve configuration.
    vector<double> TotalSpeciesAmountByAvailable(int spno, vector<int> puStatusZone) {
        vector<double> specAmount(spno, 0);
        int ism, isp;

        for (int ipu = 0; ipu < puno; ipu++)
            // only count species amounts for pu that meet a certain status
            if ((puList[ipu].richness) && (puStatusZone[ipu] >= 0) 
            && (puList[ipu].status) < 2 && (puList[ipu].fPULock != 1) && (puList[ipu].numZones == 0))
                for (int i = 0; i < puList[ipu].richness; i++)
                {
                    ism = puList[ipu].richness + i;
                    isp = puvspr[ism].spindex;
                    specAmount[isp] += puvspr[ism].amount;
                }
        
        return specAmount;
    }

    // Aggregates areas and occurrences only pu.status == 2 and pu.status == 3
    void TotalSpeciesAmountByStatus(vector<int>& TO_2, vector<int>& TO_3, vector<double>& TA_2, vector<double>& TA_3) {
        int isp,ism;

        for (int ipu=0;ipu<puno;ipu++)
        {
            if (puList[ipu].status == 2 || puList[ipu].status == 3) {
                for (int i=0;i<puList[ipu].richness;i++)
                {
                    ism = puList[ipu].offset + i;
                    isp = puvspr[ism].spindex;

                    if (puList[ipu].status == 2)
                    {
                        TO_2[isp]++;
                        TA_2[isp] += puvspr[ism].amount;
                    }

                    if (puList[ipu].status == 3)
                    {
                        TO_3[isp]++;
                        TA_3[isp] += puvspr[ism].amount;
                    }
                }
            }
        }
    }

    // * * * * * Connection Cost Type 1 * * * * *
    // Total cost of all connections for PU independant of neighbour status * * * * *
    double ConnectionCost1(int ipu)
    {
        if (connectionsEntered)
        {
            double fcost = connections[ipu].fixedcost;

            for (sneighbour &p : connections[ipu].first)
            {
                if (asymmetric)
                {
                    if (p.connectionorigon)
                        fcost += p.cost;
                }
                else
                    fcost += p.cost;
            }
            return (fcost);
        }

        return 0; // no connections entered. 
    } // * * * * Connection Cost Type 1 * * * *
    
    // returns the amount of a species at a planning unit, if the species doesn't occur here, returns 0
    // precondition - puindex and iSpexIndex must be within spno and puno.
    double RtnAmountSpecAtPu(int puindex, int iSpecIndex)
    {
        int smIndex = RtnIndexSpecAtPu(puindex, iSpecIndex);
        if (smIndex != -1) {
            return puvspr[smIndex].amount;
        }

        return 0;
    }

    // returns index of a species at puvspr
    int RtnIndexSpecAtPu(int puindex, int iSpecIndex)
    {
        if (puList[puindex].richness > 0)
        {
            auto start_it = puvspr.begin() + puList[puindex].offset;
            auto end_it = start_it + puList[puindex].richness;
            auto spindex_cmp = [](const spu &lhs, int rhs) -> bool { return lhs.spindex < rhs; };
            auto elem_it = std::lower_bound(start_it, end_it, iSpecIndex, spindex_cmp);
            if (elem_it != end_it && elem_it->spindex == iSpecIndex)
            {
                return elem_it - puvspr.begin();
            }
        }
        return -1;
    }

    // Returns all the spec info for a pu
    vector<double> RtnAmountAllSpecAtPu(int puindex, int spno)
    {
        vector<double> amounts(spno, 0.0);
        spu term;

        if (puList[puindex].richness > 0)
            for (int i=0;i<puList[puindex].richness;i++)
            {   
                term = puvspr[puList[puindex].offset + i];
                amounts[term.spindex] += term.amount;
            }

        return amounts;
    }

    // Given a pu and current zone, returns a different but valid zone for the pu.
    // Used for randomly selecting the next zone.
    int RtnValidZoneForPu(int puindex, int iPreviousZone, uniform_int_distribution<int>& randomDist, mt19937 &rngEngine, int zoneCount) {
        int chosenZoneInd,iZone;

        if (puList[puindex].numZones > 1)
        {
            // pick a random available zone for this pu that is different to the current zone
            chosenZoneInd = randomDist(rngEngine) % puList[puindex].numZones;
            iZone = puZone[puindex][chosenZoneInd];

            // if zone is already chosen, just increment it
            if (iZone == iPreviousZone)
            {
                chosenZoneInd = (chosenZoneInd + 1) % puList[puindex].numZones;
                iZone = puZone[puindex][chosenZoneInd];
            }
        }
        else if (puList[puindex].numZones == 0)
        {
            // pu can be in  any zone.
            iZone = randomDist(rngEngine) % zoneCount;
            if (iZone == iPreviousZone)
            {
                iZone = (iZone + 1) % zoneCount;
            }                          
        }
        else {
            iZone = puList[puindex].iPULock; // assume pulocked
        }

        return iZone;
    }

    void WriteSolutionsMatrixHeader(string filename, int delimType, int iIncludeHeaders) {
        ofstream myfile;
        myfile.open(filename);

        if (iIncludeHeaders == 1) 
        {
            string delim = delimType == 3 ? "," : "    ";
            myfile << "SolutionsMatrix"; // write header first time

            for (int i=0;i<puno;i++)
                myfile << delim << puList[i].id;

            myfile << "\n";

        }

        myfile.close();
    }

    vector<double>& GetCostBreakdown(int puindex) {
        return puList[puindex].costBreakdown;
    }

    // Gets the species -> pu map sorted in amount/cost order. 
    // It also populates lockedPu with the species terms for a locked pu.
    vector<vector<penaltyTerm>> getPuAmountsSorted(int spno, vector<vector<lockedPenaltyTerm>>& lockedPu) {
        vector<vector<penaltyTerm>> penaltyTerms(spno);
        lockedPu.resize(spno);
        int ipu;

        for (int i = 0; i < puno; i++) {
            ipu = puList[i].offset;

            if (puList[i].fPULock)
            {
                // append amount to lockedPu
                for (int j = 0; j < puList[i].richness; j++) {
                    // append {zoneid, puindex, amount}
                    lockedPu[puvspr[ipu+j].spindex].push_back({puList[i].iPULock, i, puvspr[ipu+j].amount});
                }
            }
            else {
                for (int j = 0; j < puList[i].richness; j++) {
                    penaltyTerms[puvspr[ipu+j].spindex].push_back({puvspr[ipu+j].amount, puList[i].cost});
                }
            }
        }

        // Sort them
        for (vector<penaltyTerm>& v: penaltyTerms) {
            sort(v.begin(), v.end(), penaltyTermSortOperator);
        }
        
        return penaltyTerms;
    }

    // Debugging functions (dumps original files)
    void DumpPuLock(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "puid,zoneid\n";
        for (auto& [puid, zoneid]: puLock)
        {
            myfile << puid << "," << zoneid << "\n";
        }
        myfile.close();
    }

    void DumpPuZone(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "puid,zoneid\n";
        for (puzonestruct& term: puZoneTemp)
        {
            myfile << term.puid << "," << term.zoneid << "\n";
        }
        myfile.close();
    }

    void DumpPuLockZoneData(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "id,fPULock,iPULock,numZones\n";
        for (int i = 0; i < puno; i++)
        {
            myfile << puList[i].id << "," << puList[i].fPULock << "," << puList[i].iPULock
                << "," << puList[i].numZones << "\n";
        }
        myfile.close();
    }

    void DumpCostValues(string filename) {
        if (puList.size() == 0)
            return;

        int costCount = puList[0].costBreakdown.size();
        ofstream myfile;
        myfile.open(filename);
        myfile << "puindex\n";
        for (int j=0; j<costCount; j++) {
            myfile << "," << j;
        }
        myfile << "\n";

        for (int i =0; i < puno; i++) {
            myfile << "," << puList[i].id;
            for (int j=0; j < costCount; j++) {
                myfile << "," << puList[i].costBreakdown[j];
            }
            myfile << "\n";
        }
        myfile.close();
    }

    int puno;
    int puLockCount;
    int puZoneCount;
    double density; // matrix density.
    int asymmetric; // whether asymmetric connectivity is on.
    bool connectionsEntered;

    // List of pu indices that can be changed (i.e. not locked to one zone or status > 1)
    vector<int> validPuIndices;

    map<int, int> puLock; // puid -> zoneid
    vector<vector<int>> puZone; 
    // puindex -> list of available zones. If pu has empty puzone, it is allowed in ALL zones (unless pulock).
    // contains the zone index (not zone id.)

    vector<spustuff> puList;
    map<int, int> lookup; // original puid -> pu index
    vector<sconnections> connections; 

    //matrix
    vector<spu> puvspr;
    vector<spusporder> puvspr_sporder;

    private:
    vector<puzonestruct> puZoneTemp;

    static bool penaltyTermSortOperator(penaltyTerm& t1, penaltyTerm& t2) {
        if (t1.cost == 0 && t2.cost == 0) {
            return t1.amount > t2.amount; // return higher amount
        }
        else if (t1.cost == 0) {
            return true; // t1 has lower cost and therefore higher priority
        }
        else if (t2.cost == 0) {
            return false; 
        }
        else {
            return t1.amount/t1.cost > t2.amount/t2.cost; //prioritise higher ratio 
        }
    }

    // constructs a vector containing only the valid pus that can be assigned zones.
    void ConstructValidPuList() {
        validPuIndices.reserve(puno-puLock.size()); // pre-allocate
        int i = 0;
        for (spustuff& pu: puList) {
            if (!((pu.status > 1) || (pu.fPULock == 1))) {
                validPuIndices.push_back(i);
            }
            i++;
        }
    }

    void PersistPuLock(map<int, ZoneName>& zoneLookup) {
        int ind, zoneind;
        for (auto& [puid, zoneid]: puLock) {
            ind = lookup[puid];
            
            auto zoneit = zoneLookup.find(zoneid);
            if (zoneit == zoneLookup.end())
                continue;

            puList[ind].fPULock = 1;
            puList[ind].iPULock = zoneit->second.index; // store the index of the zone, not the zone itself.
            puList[ind].numZones = 1;
        }
    }

    void PersistPuZone(vector<puzonestruct>& puZoneTemp, map<int, ZoneName>& zoneLookup, Logger& logger) {
        // resize puzone to puno size
        int zoneind;
        puZone.resize(puno);
        for (puzonestruct& term: puZoneTemp) {
            int puindex = lookup[term.puid];
            auto zoneit = zoneLookup.find(term.zoneid);
            if (zoneit == zoneLookup.end())
                continue; // ignore invalid zoneid

            puZone[puindex].push_back(zoneit->second.index);
        }

        // Check puZone to make sure there are no single zoned pu. 
        // we could in the future convert these to lock.
        int i = 0;
        for (vector<int>& zoneMap: puZone) {
            // exlcude already pulocked indices
            if (puList[i].numZones != 1)
            {
                puList[i].numZones = zoneMap.size();
                if (zoneMap.size() == 1)
                {
                    logger.ShowErrorMessage("Planning units found in puzone locked to a single zone. Do not use this file to lock planning units! Planning unit " + to_string(puList[i].id) + " is locked to a single zone " + to_string(zoneMap[0]) + "\n");
                }
            }
            i++;
        }
    }

    void ReadConnections(string filename, int asymmetric, Logger& logger) {
        ifstream fp;
        fp.open(filename);
        if (!fp.is_open())
            logger.ShowErrorMessage("Cannot open file " + filename + ".\n");
        connections.resize(puno);
        string sLine;
        int id1, id2, tempid1, tempid2, bad; 
        double fcost;

        bool file_is_empty = true;
        for (int line_num = 1; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (line_num == 1)
            {
                if (is_like_numerical_data(sLine))
                    logger.ShowWarningMessage("File " + filename + " has no header in the first line.\n");
                else
                    continue;
            }
            if (sLine.empty())
                continue;

            stringstream ss = stream_line(sLine);
            ss >> tempid1 >> tempid2 >> fcost;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            id1 = LookupIndex(tempid1);
            id2 = LookupIndex(tempid2);

            if (id1 == -1 || id2 == -1) // one of the puid is not found - ignore.
            {
                logger.ShowErrorMessage("A connection is out of range " + to_string(fcost) + " " + 
                    to_string(id1) + " " + to_string(id2) + "\n");
                continue;
            }

            if (id1 == id2)
            {
                if (asymmetric)
                    connections[id1].fixedcost = 0;
                else
                    connections[id1].fixedcost += fcost;

                continue;
            } /* My code for an irremovable connection */

            if (id1 >= 0 && id1 < puno)
            { /* Is I a sensible number ?*/
                bad = 0;

                if (!asymmetric)
                    for (sneighbour& p: connections[id1].first) {
                        if (p.nbr == id2)
                            bad = 1;
                    }

                if (bad)
                {
                    logger.ShowWarningMessage("Double connection definition " + to_string(puList[id1].id) + " " + to_string(puList[id2].id) + "\n");
                }
                else
                {
                    connections[id1].nbrno++;
                    sneighbour p = {
                        id2, // nbr
                        fcost, // cost
                        1 // connectionorigon
                    };

                    connections[id1].first.push_back(p);
                }
            }

            if (id2 >= 0 && id2 < puno)
            { /* Is I a sensible number ?*/
                bad = 0;

                if (!asymmetric)
                    for (sneighbour& p: connections[id2].first) {
                        if (p.nbr == id1)
                            bad = 1;
                    }

                if (!bad) 
                {
                    connections[id2].nbrno++;
                    sneighbour p = {
                        id1, // nbr
                        fcost, // cost
                        1 // connectionorigon
                    };

                    if (asymmetric)
                        p.connectionorigon = 0;

                    connections[id2].first.push_back(p);
                }
            }
        }
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;
    }

    /* Read Planning Unit Data */
    /* The pu.dat for marzone contains only the puid, and any cost columns */
    void ReadPUData(string filename, Costs& costs, Logger& logger)
    {
        ifstream fp;
        fp.open(filename);
        if (!fp.is_open())
            logger.ShowErrorMessage("Cannot open file " + filename + ".\n");
        string sLine, unusedHeader;
        stringstream errorBuf;

        getline(fp, sLine);
        vector<string> headerNames = getFileHeaders(sLine, filename, errorBuf);

        if (!errorBuf.str().empty())
            logger.ShowErrorMessage(errorBuf.str());

        /*Feed information into temporary link list */
        bool file_is_empty = true;
        for (int line_num = 2; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (sLine.empty())
                continue;
            spustuff putemp = {}; // init all to 0
            putemp.id = -1;
            putemp.cost = 0;
            putemp.xloc = -1;
            putemp.yloc = -1;
            putemp.iPULock = -1;
            putemp.costBreakdown.resize(costs.costCount, 1.0); // default cost to 1 if not supplied.

            stringstream ss = stream_line(sLine);

            for (string temp: headerNames)
            {
                if (temp.compare("id") == 0)
                {
                    ss >> putemp.id;
                }
                else if (costs.Contains(temp)) {
                    // Cost is defined
                    ss >> putemp.costBreakdown[costs.GetCostIndex(temp)];
                }
                else {
                    // un-enforced header
                    ss >> unusedHeader;
                }

            } /* looking for ivar different input variables */
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            if (putemp.id == -1)
                logger.ShowErrorMessage("ERROR: Unable to parse planning unit id for line "+ to_string(puList.size()+1) + "\n");

            // Aggregate pu costs 
            for (double c: putemp.costBreakdown) {
                putemp.cost += c;
            }

            puList.push_back(putemp);
            lookup[putemp.id] = puno;
            puno++;

        } /* while still lines in data file */

        fp.close();
        puno = puList.size();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;
    } // ReadPUData

    void LoadPuLock(string filename, Logger& logger)
    {
        ifstream fp;
        fp.open(filename);
        if (!fp.is_open())
            logger.ShowErrorMessage("Cannot open file " + filename + ".\n");

        string sLine;
        int puid, zoneid;

        bool file_is_empty = true;
        for (int line_num = 1; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (line_num == 1)
            {
                if (is_like_numerical_data(sLine))
                    logger.ShowWarningMessage("File " + filename + " has no header in the first line.\n");
                else
                    continue;
            }
            if (sLine.empty())
                continue;

            stringstream ss = stream_line(sLine);
            ss >> puid >> zoneid;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");

            // persist only if pu is already defined
            if (lookup.find(puid) != lookup.end())
            {
                puLock[puid] = zoneid;
            }
            else {
                logger.ShowWarningMessage(filename + ": puid found that was not defined in pu.dat. Puid will be ignored: " + to_string(puid) + "\n");
            }
        }
        puLockCount = puLock.size();
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;
    }

    vector<puzonestruct> LoadPuZone(string filename, Logger& logger)
    {
        ifstream fp;
        fp.open(filename);
        if (!fp.is_open())
            logger.ShowErrorMessage("Cannot open file " + filename + ".\n");

        string sLine;
        int puid, zoneid;
        vector<puzonestruct> puZoneTemp;

        bool file_is_empty = true;
        for (int line_num = 1; getline(fp, sLine); line_num++)
        {
            file_is_empty = false;
            if (line_num == 1)
            {
                if (is_like_numerical_data(sLine))
                    logger.ShowWarningMessage("File " + filename + " has no header in the first line.\n");
                else
                    continue;
            }
            if (sLine.empty())
                continue;
            stringstream ss = stream_line(sLine);
            ss >> puid >> zoneid;
            if (ss.fail())
                logger.ShowErrorMessage("File " + filename + " has incorrect values at line " + to_string(line_num) + ".\n");
            // persist only if pu is already defined
            if (lookup.find(puid) != lookup.end())
                puZoneTemp.push_back({puid, zoneid});
            else
                logger.ShowWarningMessage(filename + ": puid found that was not defined in pu.dat. Puid will be ignored: " + to_string(puid) + "\n");
        }
        puZoneCount = puZoneTemp.size();
        fp.close();
        if (file_is_empty)
            logger.ShowErrorMessage("File " + filename + " cannot be read or is empty.\n");;
        return puZoneTemp;
    }
};
    
}