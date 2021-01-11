#pragma once

// Stores information relating to the reserve configuration.
#include <limits>
#include <random>

#include "marzone.hpp"
#include "pu.hpp"
#include "zones.hpp"

namespace marzone {
using namespace std;

class Reserve {
    public:
    Reserve(Species& spec, Zones& zones) {
        InitializeZoneSpec(spec.spno, zones.zoneCount);
        speciesAmounts.resize(spec.spno, {}); // init amounts struct
    }

    // Sets all amounts and occurences to 0
    void InitializeZoneSpec(uint64_t spno, uint64_t zoneCount) {
        zoneSpec.resize(spno*zoneCount);

        for (int j=0; j<spno; j++)
        {
            for (int i=0; i<zoneCount; i++)
            {
                zoneSpec[(j*zoneCount)+i].amount = 0;
                zoneSpec[(j*zoneCount)+i].occurrence = 0;
            }
        }
    }

    void InitializeSolution(uint64_t puno) {
        solution.assign(puno, 0);
    }

    void RandomiseSolution(Pu& pu, mt19937& rngEngine) {
        uniform_int_distribution<int> random_dist(0, numeric_limits<int>::max());

        for (int i = 0; i < pu.puno; i++)
        {
            spustuff& puTerm = pu.puList[i];
            if (puTerm.fPULock == 1) {
                solution[i] = puTerm.iPULock;
                continue; // pu lock takes precedence above all.
            }

            solution[i] = 0; // default setting is 0

            if (puTerm.status > 0) // starting pu state detected
                solution[i] = puTerm.status;

            if (puTerm.numZones > 0) { // enforce puzone
                int zoneInd = random_dist(rngEngine)%puTerm.numZones;
                solution[i] = pu.puZone[i][zoneInd]-1; //TODO - change this depending on later design, but for now it is the zero-indexed zone.
            }
        }
    }
    
    void WriteSolution(string filename, Pu& pu, int imode) {
        ofstream myfile;
        myfile.open(filename);

        if (imode > 1)
            myfile << "planning_unit,zone\n";
        
        for (int i= 0; i < pu.puno; i++)
        {
            if (imode > 1)
                myfile << pu.puList[i].id << "," << solution[i]+1 << "\n";
        }

        myfile.close();
    }

    void WriteZoneConnectivitySum(string filename, Pu& pu, Zones& zones, int imode) {
        vector<vector<double>> ZCS = zones.InitializeZoneMatrix();
        //string delim = imode > 1 ? "," : "    ";

        ComputeZoneConnectivitySum(pu, ZCS); // populate ZCS

        ofstream myfile;
        myfile.open(filename);

        if (imode > 1)
            myfile << "\"Zone_Connectivity_Sum\"";
        for (int i=0;i<zones.zoneCount;i++)
            myfile << ",\"" << zones.IndexToName(i) << "\"";
        myfile << "\n";

        // write a data row for each zone
        for (int i=0;i<zones.zoneCount;i++)
        {
            myfile << "\"" << zones.IndexToName(i) << "\"";
            for (int j=0;j<zones.zoneCount;j++)
                myfile << "," << ZCS[i][j];
            myfile << "\n";
        }

        myfile.close();
    }
    
    // Given zones and pu, and using current solution.
    void ComputeZoneConnectivitySum(Pu& pu, vector<vector<double>>& ZCS) {
        for (int i=0;i<pu.puno;i++)
        {
            // fixed cost for ipu is between this zone and itself
            ZCS[solution[i]][solution[i]] += pu.connections[i].fixedcost;

            // traverse connections for this ipu
            for (sneighbour& p: pu.connections[i].first)
            {
                if (p.nbr > i) // avoid double counting connnections
                {
                    if (solution[i] != solution[p.nbr]) // ignore internal connections within a zone
                    {
                        // connections are symmetric
                        ZCS[solution[i]][solution[p.nbr]] += p.cost;
                        ZCS[solution[p.nbr]][solution[i]] += p.cost;
                    }
                }
            }
        }
    }

    // Writes amounts related to species amounts in this particular reserve to a file.
    void WriteSpeciesAmounts(string filename, Species& spec, Zones& zones, int imode, double misslevel) {
        ofstream myfile;
        myfile.open(filename);
        double rMPM, rTestMPM, rTarget, rAmountHeld;
        int rOccurrenceTarget, rOccurrencesHeld;
        string temp, d = imode > 1 ? "," : "\t";

        for (int isp=0;isp<spec.spno;isp++)
        {
            rMPM = 1;

            // Write species statis
            myfile << spec.specList[isp].name << d << spec.specList[isp].sname << d << spec.specList[isp].target;

            if (imode > 1)
                myfile << d << spec.specList[isp].totalarea;
            else
                myfile << d << speciesAmounts[isp].amount;

            myfile << speciesAmounts[isp].amount << d << spec.specList[isp].targetocc << d << speciesAmounts[isp].occurrence;
            
            temp = ReturnStringIfTargetMet(spec.specList[isp].target, speciesAmounts[isp].amount, 
                spec.specList[isp].targetocc, speciesAmounts[isp].occurrence, rMPM, misslevel);

            if (imode > 1)
            {
                if (spec.specList[isp].sepnum)
                {
                    if (speciesAmounts[isp].separation/spec.specList[isp].sepnum < misslevel)
                        temp = "no";
                }
            }
            myfile << d << temp;

            int iZoneArrayIndex;
            for (int i=0;i<zones.zoneCount;i++)
            {
                iZoneArrayIndex = (isp * zones.zoneCount) + i;
                rTarget = zones.zoneTarget[isp][i].target;
                rAmountHeld = zoneSpec[iZoneArrayIndex].amount;
                rOccurrenceTarget = zones.zoneTarget[isp][i].occurrence;
                rOccurrencesHeld = zoneSpec[iZoneArrayIndex].occurrence;

                myfile << d << rTarget << d << rAmountHeld << d << zones.zoneContribValues[iZoneArrayIndex] << d << rOccurrenceTarget << d << rOccurrencesHeld;

                temp = ReturnStringIfTargetMet(rTarget, rAmountHeld, rOccurrenceTarget, rOccurrencesHeld, rMPM, misslevel);
                myfile << d << temp;
            }

            myfile << d << rMPM << "\n";
        }

        myfile.close();
    }

    vector<zonespecstruct> zoneSpec;
    
    // Species data and amounts
    vector<reservespecies> speciesAmounts;

    // Pu configuration (assignment)
    vector<int> solution; // contains the zone index in the corresponding to a supplied Pu 

    private:
    string ReturnStringIfTargetMet(double targetArea, double area, int targetOcc, int occ, double& rMPM, double& misslevel) {
        string temp = "yes";
        double rTestMPM;

        if (targetArea)
        {
            if (area/ targetArea < misslevel)
                temp = "no";

            rTestMPM = area/ targetArea;
            if (rTestMPM < rMPM)
                rMPM = rTestMPM;
        }
        if (targetOcc)
        {
            if (occ / targetOcc < misslevel)
                temp = "no";

            rTestMPM = occ / targetOcc;
            if (rTestMPM < rMPM)
                rMPM = rTestMPM;
        }
        return temp;
    }

};

} // namespace marzone