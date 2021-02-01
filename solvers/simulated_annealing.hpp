#pragma once

#include <limits>
#include <random>

#include "../common.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"

// Debug settings
#undef TRACE_ZONE_TARGETS

namespace marzone {
using namespace std;

class SimulatedAnnealing {
    public:
    SimulatedAnnealing(int annealingOn, sanneal& anneal, mt19937& rngEngine, int annealTrace, int id) 
    : rngEngine(rngEngine), annealTraceType(annealTrace), id(id)
    {
        if (annealingOn)
        {
            settings = anneal; // copy construct
        }
    }

    void Initialize(Species& spec, Pu& pu, Zones& zones, int clumptype) {
        if (settings.type >= 2)
        {
            if (settings.type == 2)
            {
                ConnollyInit(spec, pu, zones, clumptype);
            }
            else if (settings.type == 3)
            {
                AdaptiveInit(spec, pu, zones, clumptype);
            }
        }

        settings.temp = settings.Tinit;
    }
    
    // TODO - complete.
    void DumpBuffer() {
    }

    void RunAnneal(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh) {
        uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());
        uniform_int_distribution<int> randomPuDist(0, pu.puno-1);
        uniform_real_distribution<double> float_range(0.0, 1.0);
        int ipu, iZone, itemp, iPreviousZone, iGoodChange;
        uint64_t ichanges;

        for (int itime = 1; itime <= settings.iterations; itime++)
        {
            // toggle a random planning unit between reserved and available
            pair<int, int> next = GetPuAndZone(r, pu, randomPuDist, randomDist, zones.zoneCount);
            ipu = next.first;
            iZone = next.second;

            itemp = 1;
            iPreviousZone = r.solution[ipu];

#ifdef TRACE_ZONE_TARGETS
            debugbuffer << "annealing time " << itime << " of " << settings.iterations << "\n";
#endif

            schange change = r.CheckChangeValue(ipu, iPreviousZone, iZone, pu, zones, spec, costthresh, 
            tpf1, tpf2, (double)itime / (double)settings.iterations);

#ifdef TRACE_ZONE_TARGETS
            debugbuffer << "annealing after CheckChange\n";
#endif

            /* Need to calculate Appropriate temperature in GoodChange or another function */
            /* Upgrade temperature */
            if (itime % settings.Tlen == 0)
            {
                if (settings.type == 3)
                    AdaptiveDec();
                else
                    settings.temp = settings.temp * settings.Tcool;
                
                /* TODO enable after logging lib
                ShowDetProg("time %i temp %f Complete %i%% currval %.4f\n",
                            itime, anneal.temp, (int)itime * 100 / anneal.iterations, reserve->total);
                */
            } /* reduce temperature */

            /* TODO - enable 
            if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
            {
                if (repeats > 1)
                    sprintf(tempname1, "_r%05i", irun);
                else
                    tempname1[0] = 0;
                if (fnames.savesnapchanges == 3)
                    sprintf(tempname2, "%s_snap%st%05i.csv", savename, tempname1, ++snapcount % 10000);
                else if (fnames.savesnapchanges == 2)
                    sprintf(tempname2, "%s_snap%st%05i.txt", savename, tempname1, ++snapcount % 10000);
                else
                    sprintf(tempname2, "%s_snap%st%05i.dat", savename, tempname1, ++snapcount % 10000);

                OutputSolution(puno, R, pu, tempname2, fnames.savesnapsteps);
            } /* Save snapshot every savesnapfreq timesteps */
            
            if (GoodChange(change, float_range) == 1)
            {
                iGoodChange = 1;
                ++ichanges;

                r.ApplyChange(ipu, iZone, change, pu, zones, spec);

                /* TODO - re-enable
                if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
                {
                    if (repeats > 1)
                        sprintf(tempname1, "_r%05i", irun);
                    else
                        tempname1[0] = 0;
                    if (fnames.savesnapchanges == 3)
                        sprintf(tempname2, "%s_snap%sc%05i.csv", savename, tempname1, ++snapcount % 10000);
                    else if (fnames.savesnapchanges == 2)
                        sprintf(tempname2, "%s_snap%sc%05i.txt", savename, tempname1, ++snapcount % 10000);
                    else
                        sprintf(tempname2, "%s_snap%sc%05i.dat", savename, tempname1, ++snapcount % 10000);

                    OutputSolution(puno, R, pu, tempname2, fnames.savesnapchanges);
                } /* Save snapshot every savesnapfreq changes */

            } /* Good change has been made */
            else
                iGoodChange = 0;

#ifdef TRACE_ZONE_TARGETS
            debugBuffer << "annealing after DoChange\n";
#endif

            if (settings.type == 3)
            {
                settings.sum += r.objective.total;
                settings.sum2 += r.objective.total * r.objective.total;
            } /* Keep track of scores for averaging stuff */
/*
#ifdef DEBUGTRACEFILE
            if (verbose > 4)
                fprintf(fp, "%i,%i,%i,%i,%i,%i,%i,%f,%f,%f,%f,%f\n", itime, ipu, pu[ipu].id, iPreviousR, itemp, iZone, iGoodChange, change->total, change->cost, change->connection, change->penalty, anneal.temp);
#endif

            if (fnames.saveannealingtrace)
            {
                iRowCounter++;
                if (iRowCounter > iRowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    if (fnames.suppressannealzones == 0)
                        fprintf(zonefp, "%i", itime);

                    if (fnames.saveannealingtrace > 1)
                    {
                        fprintf(ttfp, "%i,%f,%i,%f,%f,%f,%f,%f,%i,%f\n", itime, costthresh, iGoodChange, reserve->total, reserve->cost, reserve->connection, reserve->penalty, reserve->shortfall, ipu, anneal.temp); // itime,costthresh,cost,connection,penalty

#ifdef DEBUG_PEW_CHANGE_PEN
                        AppendDebugTraceFile("iteration,threshold,dochange,total,cost,connection,penalty,puindex\n");
                        sprintf(debugbuffer, "%i,%f,%i,%f,%f,%f,%f,%i\n", itime, costthresh, iGoodChange, reserve->total, reserve->cost, reserve->connection, reserve->penalty, ipu);
                        AppendDebugTraceFile(debugbuffer);
#endif

                        if (fnames.suppressannealzones == 0)
                            for (i = 0; i < puno; i++)
                                fprintf(zonefp, ",%i", R[i]);
                    }
                    else
                    {
                        fprintf(ttfp, "%i %f %i %f %f %f %f %f %i %f\n", itime, costthresh, iGoodChange, reserve->total, reserve->cost, reserve->connection, reserve->penalty, reserve->shortfall, ipu, anneal.temp);

                        if (fnames.suppressannealzones == 0)
                            for (i = 0; i < puno; i++)
                                fprintf(zonefp, " %i", R[i]);
                    }

                    if (fnames.suppressannealzones == 0)
                        fprintf(zonefp, "\n");
                }
            }
            */

        } /* Run Through Annealing */
    }

    // Randomly chooses the next pu and zone to flip
    // Should only return puindex that is not locked, and a zone that is different from current zone.
    pair<int, int> GetPuAndZone(Reserve& r, Pu& pu, uniform_int_distribution<int>& randomPuDist, uniform_int_distribution<int>& randomDist, int zoneCount) {
        // pick pu at random
        int ipu,iZone, iPreviousR; 

        do
        {
            ipu = randomPuDist(rngEngine);
        } while ((pu.puList[ipu].status > 1) || (pu.puList[ipu].fPULock == 1)); // ignore locked/ non avail pu

        iPreviousR = r.solution[ipu];
        iZone = pu.RtnValidZoneForPu(ipu, iPreviousR, randomDist, rngEngine, zoneCount);

        return pair<int,int>(ipu, iZone);
    }

    sanneal settings;

    private:
    // * * * * Good Change * * * *
    int GoodChange(schange& change, uniform_real_distribution<double>& float_range)
    {
        // avoid math operations if possible
        if (change.total < 0) {
            return 1;
        }

        return (exp(-change.total / settings.temp) > float_range(rngEngine)) ? 1 : 0;

    } // Good Change

    /**** Function to decrement T and decide if it is time to stop *****/
    /* Adaptive Decrement. Sets the new temperature based on old values */
    void AdaptiveDec() {
        double omega = 0.7; /* Control parameter */
        double sigmanew, sigmamod;
        double lambda = 0.7; /* control parameter*/

        sigmanew = (settings.sum2 - pow((settings.sum / settings.Tlen), 2)) / (settings.Tlen - 1);
        sigmamod = (1 - omega) * sigmanew + omega * settings.sigma * (settings.temp / settings.tempold);
        settings.tempold = settings.temp;
        settings.temp = exp(-lambda * settings.temp / sigmamod);
        settings.sigma = sigmamod;
        settings.sum = 0;
        settings.sum2 = 0;
    }

    void ConnollyInit(Species& spec, Pu& pu, Zones& zones, int clumptype) {
        double localdelta = numeric_limits<double>::epsilon();
        uniform_int_distribution<int> random_dist(0, numeric_limits<int>::max());
        uniform_int_distribution<int> random_pu_dist(0, pu.puno-1);

        // Set reserve to a random and evaluate 
        Reserve r(spec, zones.zoneCount, clumptype);
        r.InitializeSolution(pu.puno);
        r.RandomiseSolution(pu, rngEngine, zones.zoneCount);
        r.EvaluateObjectiveValue(pu, spec, zones);

        int ipu, iPreviousR, chosenZoneInd, chosenZone, imode = 1;
        double deltamin = 0,deltamax = 0;
        for (int i=1;i<= settings.iterations/100; i++)
        {
            // pick pu at random
            do {
                ipu = random_pu_dist(rngEngine);
            } while ((pu.puList[ipu].status > 1) || (pu.puList[ipu].fPULock == 1)); // ignore locked/ non avail pu

            if (pu.puList[ipu].numZones > 1)
            {
                // pick a random available zone for this pu that is different to the current zone
                chosenZoneInd = random_dist(rngEngine) % pu.puList[ipu].numZones;
                chosenZone = pu.puZone[ipu][chosenZoneInd]-1;

                // if zone is already chosen, just increment it
                if (chosenZone == iPreviousR)
                {
                    chosenZoneInd = (chosenZoneInd + 1)%pu.puList[ipu].numZones;
                    chosenZone = pu.puZone[ipu][chosenZoneInd]-1;
                }
            }
            else if (pu.puList[ipu].numZones == 0) {
                // pu can be in  any zone.
                chosenZone = random_dist(rngEngine) % zones.zoneCount;
                if (chosenZone == iPreviousR)
                {
                    chosenZone = (chosenZone + 1)%zones.zoneCount;
                }
            }

            // Check change in zone on penalty
            schange change = r.CheckChangeValue(ipu, r.solution[ipu], chosenZone, pu, zones, spec, 0);

            // apply change
            r.ApplyChange(ipu, chosenZone, change, pu, zones, spec);

            if (change.total > deltamax)
                deltamax = change.total;
            if (change.total >localdelta && (deltamin < localdelta || change.total < deltamin))
                deltamin = change.total;
        }

        // Set snneal settings
        settings.Tinit = deltamax;
        if (deltamax)
        {
            if (settings.Titns)
            {
                settings.Tcool = exp(log(deltamin / settings.Tinit) / (double)settings.Titns);
            }
            else
            {
                settings.Tcool = 1;
            }
        }

        /*
        #ifdef DEBUGTRACEFILE
        sprintf(debugbuffer,"Tinit %g Titns %i Tcool %g\n",settings.Tinit,settings.Titns,settings.Tcool);
        AppendDebugTraceFile(debugbuffer);
        AppendDebugTraceFile("ConnollyInit end\n");
        #endif
        */
    }

    /* * * * Adaptive Annealing 2 * * * * * * * * *****/
    /**** Initial Trial Runs. Run for some small time to establish sigma. ****/
    void AdaptiveInit(Species& spec, Pu& pu, Zones& zones, int clumptype) {
        int i, isamples = 1000; /* Hardwired number of samples to take */
        double sum = 0, sum2 = 0;
        double sigma;
        Reserve r(spec, zones.zoneCount, clumptype);
        double c = 10; /* An initial temperature acceptance number */

        r.InitializeSolution(pu.puno);
        for (i = 0; i < isamples; i++)
        { /* Generate Random Reserve */
            r.RandomiseSolution(pu, rngEngine, zones.zoneCount);
            /* Score Random reserve */
            r.EvaluateObjectiveValue(pu, spec, zones);

            /* Add Score to Sum */
            sum += r.objective.total;
            sum2 += r.objective.total * r.objective.total;
        } /* Sample space iterations/100 times */

        sigma = sqrt(sum2 - pow(sum / isamples, 2)) / (isamples - 1);

        settings.Tinit = c * sigma;
        settings.sigma = sigma;
        settings.temp = settings.Tinit;
        settings.tempold = settings.temp;
        settings.sum = 0;
        settings.sum2 = 0;
    }

    ostringstream debugBuffer; // TODO - print buffer once logging library
    mt19937 &rngEngine;
    int id;

    // anneal trace settings
    int annealTraceType;
    string annealTraceFileName;
    ostringstream annealTraceBuffer;
};

} // namespace marzone