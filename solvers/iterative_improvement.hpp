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
        // TODO - init trace files based on fnames.
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
        // FILE *ttfp, *zonefp; - move to other func
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
                //schange change = r.CheckChangeValue(ichoice, iPreviousR, iZone, pu, zones, spec, costthresh, tpf1, tpf2);
                r.CheckChangeValue(change, ichoice, iPreviousR, iZone, pu, zones, spec, costthresh, tpf1, tpf2);
                if (change.total < 0) {
                    ichanges++;
                    r.ApplyChange(ichoice, iZone, change, pu, zones, spec);
                }

                if (!itImpTraceFileName.empty())
                    WriteItImpTraceRow();
            }
        }
    }

    // TODO
    void InitSaveItImpTrace() {
        /*
            if (fnames.saveitimptrace)
    {
        if (fnames.saveitimptrace==3)
        {
            sprintf(tempname2,"%s_itimp_objective%05i.csv",savename,irun%10000);
        } else {
            if (fnames.saveitimptrace==2)
                sprintf(tempname2,"%s_itimp_objective%05i.txt",savename,irun%10000);
            else
                sprintf(tempname2,"%s_itimp_objective%05i.dat",savename,irun%10000);
        }

        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcpy(writename,tempname2);
        if ((ttfp = fopen(writename,"w"))==NULL)
            ShowErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
        if (fnames.saveitimptrace > 1)
            fprintf(ttfp,"improvement,total,cost,connection,penalty\n");
        else
            fprintf(ttfp,"improvement total cost connection penalty\n");

        if (fnames.saveitimptrace==3)
            sprintf(tempname2,"%s_itimp_zones%05i.csv",savename,irun%10000);
        else
        if (fnames.saveitimptrace==2)
            sprintf(tempname2,"%s_itimp_zones%05i.txt",savename,irun%10000);
        else
            sprintf(tempname2,"%s_itimp_zones%05i.dat",savename,irun%10000);

        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcat(writename,tempname2);
        if ((zonefp = fopen(writename,"w"))==NULL)
            ShowErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
        fprintf(zonefp,"configuration");
        if (fnames.saveitimptrace > 1)
        {
            for (i = 0;i<puno;i++)
                fprintf(zonefp,",%i",pu[i].id);
            fprintf(zonefp,"\n0");

            for (i = 0;i<puno;i++)
                fprintf(zonefp,",%i",R[i]);
        } else {
            for (i = 0;i<puno;i++)
                fprintf(zonefp," %i",pu[i].id);
            fprintf(zonefp,"\n0");

            for (i = 0;i<puno;i++)
                fprintf(zonefp," %i",R[i]);
        }
        fprintf(zonefp,"\n");

        iRowCounter = 0;
        if (fnames.itimptracerows == 0)
            iRowLimit = 0;
        else
            iRowLimit = floor(puvalid / fnames.itimptracerows);
    }

        */
       
    }
    
    void WriteItImpTraceRow() {
        /*
                        if (fnames.saveitimptrace)
            {
                iRowCounter++;
                if (iRowCounter > iRowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    fprintf(zonefp,"%i",i);

                    if (fnames.saveitimptrace > 1)
                    {
                        fprintf(ttfp,"%i,%f,%f,%f,%f\n",
                                    i,reserve->total,
                                    reserve->cost,reserve->connection,reserve->penalty); // i,costthresh,cost,connection,penalty

                        for (j = 0;j<puno;j++)
                            fprintf(zonefp,",%i",R[j]);
                    } else {
                        fprintf(ttfp,"%i %f %f %f %f\n",
                                     i,reserve->total,reserve->cost,reserve->connection,reserve->penalty);

                        for (j = 0;j<puno;j++)
                            fprintf(zonefp," %i",R[j]);
                    }

                    fprintf(zonefp,"\n");
                }
            }
        }/* no untested PUs left */
        
    }

    ostringstream debugBuffer; // TODO - print buffer once logging library
    mt19937 &rngEngine;
    string itImpTraceFileName;
    int iterativeImprovementMode;
    int optimise;

};

} // namespace marzone