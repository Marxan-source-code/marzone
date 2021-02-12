#pragma once

#include <chrono>
#include <limits>
#include <random>

#include "../marzone.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"

namespace marzone {
using namespace std;

/* * * * ***** Main Iterative Improvement Engine * * * * * * * * ****/
class IterativeImprovement {
    public:
    IterativeImprovement(mt19937& rngEngine, sfname& fnames, int iterativeImprovementMode) 
    : rngEngine(rngEngine), iterativeImprovementMode(iterativeImprovementMode), optimise(optimise)
    {
        savename = fnames.savename;
        saveItimpTrace = fnames.saveitimptrace;
        itimpTraceRows = fnames.itimptracerows;
    }

    void Run(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh) {
        // run iterative improvement
        RunOptimisedItImp(r, spec, pu, zones, tpf1, tpf2, costthresh);

        if (iterativeImprovementMode == 3) // run again if mode = 3
            RunOptimisedItImp(r, spec, pu, zones, tpf1, tpf2, costthresh);
    }

    vector<double> rare;

    private:
    void RunOptimisedItImp(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh) {
        int puvalid = 0, ipu = 0, imode, ichoice, iZone, iPreviousR, ichanges = 0;
        struct iimp *iimparray;
        double debugfloat;
        ofstream ttfp, zonefp;
        string writename;

        int iSamplesForEachPu = (zones.zoneCount-1)*2; // allow sampling to each zone and back to available for each non available zone
        vector<int> validPuIndicies; // list of valid pu indicies that we want to modify

        // counting pu's we need to test
        for (int i = 0; i < pu.puno; i++)
        {
            if ((r.solution[ipu] >= 0) && (pu.puList[ipu].status < 2) && 
            (pu.puList[ipu].fPULock != 1) && (pu.puList[ipu].numZones != 1))
            {
                puvalid += iSamplesForEachPu;
                for (int j= 0; j < iSamplesForEachPu; j++)
                    validPuIndicies.push_back(i);
            }
        }

        schange change = r.InitializeChange(spec, zones);
        InitSaveItImpTrace(ttfp, zonefp, pu, r, validPuIndicies.size());

        if (puvalid > 0)
        {
            // Shuffle the validPuIndices array
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            shuffle(validPuIndicies.begin(), validPuIndicies.end(), default_random_engine(seed));
            uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());;

            for (int& ichoice: validPuIndicies) {
                // Pick a valid zone for the current pu
                iPreviousR = r.solution[ichoice];
                iZone = pu.RtnValidZoneForPu(ichoice, iPreviousR, randomDist, rngEngine, zones.zoneCount);

                // Check change value
                r.CheckChangeValue(change, ichoice, iPreviousR, iZone, pu, zones, spec, costthresh, tpf1, tpf2);
                if (change.total < 0) {
                    ichanges++;
                    r.ApplyChange(ichoice, iZone, change, pu, zones, spec);
                }

                if (saveItimpTrace)
                    WriteItImpTraceRow(ttfp, zonefp, r, ichoice);
            }
        }

        CloseTraceFiles(ttfp, zonefp);
    }

    void InitSaveItImpTrace(ofstream& ttfp, ofstream& zonefp, Pu& pu, Reserve& r, int puvalid) {
        string tempname;
        string d = saveItimpTrace > 1 ? "," : " ";
        if (saveItimpTrace)
        {
            tempname = savename + "_itimp_objective" + to_string(r.id) + getFileSuffix(saveItimpTrace);
            ttfp.open(tempname);
            ttfp << "improvement" << d << "total" << d << "cost" << d << "connection" << d << "penalty\n";

            tempname = savename + "_itimp_zones" + to_string(r.id) + getFileSuffix(saveItimpTrace);
            zonefp.open(tempname);
            zonefp << "configuration";

            for (int i = 0; i < pu.puno; i++)
            {
                zonefp << d << pu.puList[i].id;
            }
            zonefp << "\n";

            for (int i = 0; i < pu.puno; i++)
            {
                zonefp << d << r.solution[i];
            }
            zonefp << "\n";

            rowCounter = 0;
            if (itimpTraceRows == 0)
                rowLimit = 0;
            else
                rowLimit = floor(puvalid / itimpTraceRows);
        }
    }

    void CloseTraceFiles(ofstream& ttfp, ofstream& zonefp) {
        if (saveItimpTrace)
        {
            ttfp.close();
            zonefp.close();
        }
    }
    
    void WriteItImpTraceRow(ofstream& ttfp, ofstream& zonefp, Reserve& r, int ipu) {
        rowCounter++;
        string d = saveItimpTrace > 1 ? "," : " ";
        if (rowCounter > rowLimit)
            rowCounter = 1;

        if (rowCounter == 1)
        {
            zonefp << ipu;
            // i,costthresh,cost,connection,penalty
            ttfp << ipu << d << r.objective.total << d << r.objective.cost << d << r.objective.connection << d << r.objective.penalty << "\n";
            for (int j = 0; j <r.solution.size(); j++)
                zonefp << d << r.solution[j];
            zonefp << "\n";
        }
    }

    mt19937 &rngEngine;
    int iterativeImprovementMode;
    int saveItimpTrace;
    int optimise;
    unsigned rowCounter;
    unsigned itimpTraceRows;
    unsigned rowLimit;

    string savename;
};

} // namespace marzone