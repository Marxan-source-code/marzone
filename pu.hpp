#pragma once

// Contains the pu class
// Handles all operations relating to pu, puzone, pulock etc. and pu connections.
// Also includes the sparse matrix. 

#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include "costs.hpp"
#include "common.hpp"
#include "species.hpp"
#include "util.hpp"

namespace marzone {
using namespace std;

class Pu {
    public:
    Pu(sfname& fnames, Costs& costs, int asymmetric) : puno(0), puLockCount(0), puZoneCount(0) {
        ReadPUData(fnames.inputdir + fnames.puname, costs);

        if (!fnames.pulockname.empty())
        {
            LoadPuLock(fnames.inputdir + fnames.pulockname);
        }

        // TODO 
        /*
        if (iVerbosity > 3)
            DumpCostValues(iCostCount,puno,CostValues,fnames);
        */

       // Persist pulock and puzone counts
       PersistPuLock();

        // Read and persist puzone.
        if (!fnames.puzonename.empty())
        {
            vector<puzonestruct> puZoneTemp = LoadPuZone(fnames.inputdir + fnames.puzonename);
            PersistPuZone(puZoneTemp);
        }

       // Read and set connections
        if (!fnames.connectionname.empty())
        {
            connectionsEntered = true;
            ReadConnections(fnames.inputdir + fnames.connectionname, asymmetric);
        }
        else {
            connectionsEntered = false;
        }

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

    void LoadSparseMatrix(Species& spec, string filename) {
        FILE *fp = openFile(filename);
        vector<map<int,spu>> SMTemp(puno); // temporarily storing in this structure prevents the need for ordering.
        char sLine[1000];
        char *sVarVal;
        int _puid, _spid, puindex, spindex;
        double amount;

        fgets(sLine,999,fp); // skip header
        while (fgets(sLine, 999, fp))
        {
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &_spid);
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &_puid);
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &amount);

            spindex = spec.LookupIndex(_spid);
            if (spindex == -1)
                continue;

            puindex = LookupIndex(_puid);
            if (puindex == -1)
                continue;

            /* increment richness for planning unit containing this feature */
            puList[puindex].richness += 1;

            spu temp;
            temp.amount = amount;
            temp.clump = 0;
            SMTemp[puindex][spindex] = temp;
        }

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

    // TODO - delete as no longer needed.
    void LoadSparseMatrix_sporder(Species& spec, string filename) {
        FILE *fp = openFile(filename);
        vector<map<int,spusporder>> SMTemp(spec.spno); // temporarily storing in this structure prevents the need for ordering.
        char sLine[1000];
        char *sVarVal;
        int _puid, _spid, puindex, spindex;
        double amount;

        fgets(sLine,999,fp); // skip header
        while (fgets(sLine, 999, fp))
        {
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &_spid);
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &_puid);
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &amount);

            puindex = LookupIndex(_puid);
            spindex = spec.LookupIndex(_spid);

            if (puindex == -1 || spindex == -1)
                continue;

            SMTemp[spindex][puindex].amount = amount;
        }

        // load all into real spu vector and set species richness. 
        int j = 0;
        for (int i = 0; i < spec.spno; i++) {
            spec.SetOffset(i, j); // set offset j to species i
            spec.SetRichness(i, SMTemp[i].size());
            for (auto &[puindex, val] : SMTemp[i]) {
                spusporder temp = val;
                temp.puindex = puindex;
                puvspr_sporder.push_back(temp);
                j++;
            }
        }
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
            iZone = puZone[puindex][chosenZoneInd] - 1;

            // if zone is already chosen, just increment it
            if (iZone == iPreviousZone)
            {
                chosenZoneInd = (chosenZoneInd + 1) % puList[puindex].numZones;
                iZone = puZone[puindex][chosenZoneInd] - 1;
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

    // Given a pu index, returns -1 if pu is not locked, or zoneid if pu is locked.
    int GetPuLock(int puindex) {
        if (puList[puindex].fPULock) {
            return puList[puindex].iPULock;
        }
        return -1;
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
        for (int i = 0; i < puno; i++)
        {
            auto& zoneList = puZone[i];
            for (auto& zoneid: zoneList)
                myfile << puList[i].id << "," << zoneid << "\n";
        }
        myfile.close();
    }

    int puno;
    int puLockCount;
    int puZoneCount;
    double density; // matrix density.
    int asymmetric; // whether asymmetric connectivity is on.
    bool connectionsEntered;

    // TODO - make private.
    map<int, int> puLock; // puid -> zoneid
    vector<vector<int>> puZone; // puindex -> list of available zones. If pu has empty puzone, it is allowed in ALL zones (unless pulock).
    vector<spustuff> puList;
    map<int, int> lookup; // original puid -> pu index
    vector<sconnections> connections; 

    //matrix
    vector<spu> puvspr;
    vector<spusporder> puvspr_sporder;

    private:

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

    void PersistPuLock() {
        for (auto& [puid, zoneid]: puLock) {
            int ind = lookup[puid];
            puList[ind].fPULock = 1;
            puList[ind].iPULock = zoneid;
            puList[ind].numZones = 1;
        }
    }

    void PersistPuZone(vector<puzonestruct>& puZoneTemp) {
        // resize puzone to puno size
        puZone.resize(puno);
        for (puzonestruct& term: puZoneTemp) {
            int puindex = lookup[term.puid];
            puZone[puindex].push_back(term.zoneid); // for the public map, we use actual zone number
        }

        // Check puZone to make sure there are no single zoned pu. TODO - If so, convert to lock.
        int i = 0;
        for (vector<int>& zoneMap: puZone) {
            // exlcude already pulocked indices
            if (puList[i].numZones != 1)
            {
                puList[i].numZones = zoneMap.size();
                if (zoneMap.size() == 1)
                {
                    throw new invalid_argument("Planning units found in puzone locked to a single zone. Do not use this file to lock planning units!");
                    // TODO - add this when logging library is complete.
                    // ShowErrorMessage("Error: planning unit %i is locked to a single zone in %s.\nDo not use this file to lock planning units to a single zone; use pulock.dat for this purpose.\nAborting Program.\n",pu[i].id,fnames.puzonename);
                }
            }
            i++;
        }
    }

    void ReadConnections(string filename, int asymmetric) {
        FILE *fp = openFile(filename);
        connections.resize(puno);
        char sLine[1000];
        char *sVarVal;
        int id1, id2, tempid1, tempid2, bad; 
        double fcost;

        fgets(sLine,999,fp); /* Skip header line */

        while (fgets(sLine, 999, fp))
        {
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &tempid1);
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &tempid2);
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%lf", &fcost);

            id1 = LookupIndex(tempid1);
            id2 = LookupIndex(tempid2);

            if (id1 == -1 || id2 == -1) // one of the puid is not found - ignore.
            {
                //    ShowErrorMessage("A connection is out of range %f %i %i \n", fcost, id1, id2);
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
                    cout << "Double connection definition " << puList[id1].id << " " << puList[id2].id << "\n";
                    // TODO - re-enable.
                    //ShowDetProg("Double connection definition %i %i \n", pu[id1].id, pu[id2].id);
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
        fclose(fp);
    }

    /* Read Planning Unit Data */
    /* The pu.dat for marzone contains only the puid, and any cost columns */
    void ReadPUData(string filename, Costs& costs)
    {
        FILE *fp = openFile(filename);
        char sLine[1000];
        char *sVarVal;

        /* Scan header. We expect id header, and number of cost headers */
        fgets(sLine, 999, fp);
        vector<string> headerNames = getFileHeaders(sLine, filename);

        /* While there are still lines left feed information into temporary link list */
        while (fgets(sLine, 999, fp))
        {
            spustuff putemp = {}; // init all to 0
            putemp.id = -1;
            putemp.cost = 0;
            putemp.xloc = -1;
            putemp.yloc = -1;
            putemp.costBreakdown.resize(costs.costCount, 1.0); // default cost to 1 if not supplied.

            for (string temp: headerNames)
            {
                if (temp == headerNames.front())
                    sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
                else
                    sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");

                if (temp.compare("id") == 0)
                {
                    sscanf(sVarVal, "%d", &putemp.id);
                }

                // Cost is defined
                if (costs.Contains(temp)) {
                    sscanf(sVarVal, "%lf", &putemp.costBreakdown[costs.GetCostIndex(temp)]);
                }

            } /* looking for ivar different input variables */

            if (putemp.id == -1)
                throw invalid_argument("ERROR: Unable to parse planning unit id for line "+ to_string(puList.size()+1) + "\n");

            // Aggregate pu costs 
            for (double c: putemp.costBreakdown) {
                putemp.cost += c;
            }

            puList.push_back(putemp);
            lookup[putemp.id] = puno;
            puno++;

        } /* while still lines in data file */

        fclose(fp);
        puno = puList.size();
    } // ReadPUData

    void LoadPuLock(string filename)
    {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int puid, zoneid;

        // count the number of records
        fgets(sLine, 999, fp);
        // create the PuLock map
        while (fgets(sLine, 999, fp))
        {
            // read the integer puid from this line
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &puid);

            // read the integer zoneid from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // persist only if pu is already defined
            if (lookup.find(puid) != lookup.end())
            {
                puLock[puid] = zoneid;
            }
        }
        puLockCount = puLock.size();
        fclose(fp);
    }

    vector<puzonestruct> LoadPuZone(string filename)
    {
        FILE *fp = openFile(filename);
        char sLine[1000], *sVarVal;
        int puid, zoneid;
        vector<puzonestruct> puZoneTemp;

        // skip header
        fgets(sLine, 999, fp);

        // create the PuLock array
        while (fgets(sLine, 999, fp))
        {
            // read the integer puid from this line
            sVarVal = strtok(sLine, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &puid);

            // read the integer zoneid from this line
            sVarVal = strtok(NULL, " ,;:^*\"/|\t\'\\\n");
            sscanf(sVarVal, "%d", &zoneid);

            // persist only if pu is already defined
            if (lookup.find(puid) != lookup.end())
                puZoneTemp.push_back({puid, zoneid});
        }
        puZoneCount = puZoneTemp.size();
        fclose(fp);

        return puZoneTemp;
    }

};
    
}