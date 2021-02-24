#pragma once

#include <chrono>
#include <limits>
#include <random>

#include "../common.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"
#include "../logger.hpp"

// Debug settings
#undef TRACE_ZONE_TARGETS

namespace marzone {
using namespace std;

class SimulatedAnnealing {
    public:
    SimulatedAnnealing(sfname& fnames, LoggerBase& logger, int annealingOn, sanneal& anneal, mt19937& rngEngine, int annealTrace, int id) 
    : rngEngine(rngEngine), annealTraceType(annealTrace), id(id)
    {
        logger = logger;
        if (annealingOn)
        {
            settings = anneal; // copy construct
        }
        else {
            settings = {}; // empty/0
        }

        saveSnapSteps = fnames.savesnapsteps;
        saveSnapChanges = fnames.savesnapchanges;
        saveSnapFrequency = fnames.savesnapfrequency;
        suppressAnnealZones = fnames.suppressannealzones;
        savename = fnames.savename;
        annealingTraceRows = fnames.annealingtracerows;

        if (annealingTraceRows) {
            rowLimit = settings.iterations/annealingTraceRows;
        }
    }

    void Initialize(Species& spec, Pu& pu, Zones& zones, int clumptype, double blm) {
        if (settings.type >= 2)
        {
            if (settings.type == 2)
            {
                ConnollyInit(spec, pu, zones, clumptype, blm);
            }
            else if (settings.type == 3)
            {
                AdaptiveInit(spec, pu, zones, clumptype, blm);
            }
        }

        settings.temp = settings.Tinit;
    }

    void RunAnneal(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh, double blm) {
        uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());
        uniform_int_distribution<int> randomPuDist(0, pu.validPuIndices.size()-1);
        uniform_real_distribution<double> float_range(0.0, 1.0);
        int ipu, iZone, itemp, iPreviousZone, iGoodChange, iRowCounter = 0;
        uint64_t ichanges = 0, snapcount = 0;
        string tempname, paddedRun = intToPaddedString(r.id, 5);
        ofstream zonefp, tracefp;
        ostringstream debugBuffer, annealTraceBuffer; // debug buffer and progress buffer

        // Initialize objects and files (if any)
        schange change = r.InitializeChange(spec, zones);

        if (annealTraceType)
        {
            tempname = savename + "_anneal_objective" + paddedRun + getFileSuffix(annealTraceType);
            tracefp.open(tempname);
            WriteAnnealTraceHeader(tracefp, r, costthresh);

            if (suppressAnnealZones == 0)
            {
                tempname = savename + "_anneal_zones" + paddedRun + getFileSuffix(annealTraceType);
                zonefp.open(tempname);
                WriteZoneTraceHeader(zonefp, r, pu);
            }
        }

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
            r.CheckChangeValue(change, ipu, iPreviousZone, iZone, pu, zones, spec, costthresh, blm, 
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

                if (logger.GetVerbosity() > 3)
                {
                    annealTraceBuffer << "run number: " << r.id << " time " << itime << " temp " << settings.temp << " Complete "
                        << (int)itime * 100 / settings.iterations << " currval " << r.objective.total << "\n";
                    logger.ShowWarningMessage(annealTraceBuffer.str());
                    annealTraceBuffer.clear();
                }

            } /* reduce temperature */

            if (saveSnapSteps && !(itime % saveSnapFrequency))
            {
                tempname = savename + "_snap_r" + paddedRun
                    + "t" + intToPaddedString(++snapcount, 5) + getFileSuffix(saveSnapChanges); 

                r.WriteSolution(tempname, pu, zones, saveSnapSteps);
            } /* Save snapshot every savesnapfreq timesteps */
            
            if (GoodChange(change, float_range) == 1)
            {
                iGoodChange = 1;
                ++ichanges;

                r.ApplyChange(ipu, iZone, change, pu, zones, spec);

                if (saveSnapChanges && !(ichanges % saveSnapFrequency))
                {
                    tempname = savename + "_snap_r" + paddedRun + "c" 
                        + intToPaddedString(++snapcount, 5) + getFileSuffix(saveSnapChanges);

                    r.WriteSolution(tempname, pu, zones, saveSnapChanges);
                } /* Save snapshot every savesnapfreq changes */

            } /* Good change has been made */
            else
                iGoodChange = 0;

#ifdef TRACE_ZONE_TARGETS
            debugBuffer << "annealing after DoChange\n";
            logger.AppendDebugTraceFile(debugBuffer.str());
            debugBuffer.clear();
#endif

            if (settings.type == 3)
            {
                settings.sum += r.objective.total;
                settings.sum2 += r.objective.total * r.objective.total;
            } /* Keep track of scores for averaging stuff */

            if (logger.GetVerbosity() > 4)
            {
                debugBuffer << "run: " << r.id << "," << itime << "," << ipu << ","
                    << pu.puList[ipu].id << "," << iPreviousZone << "," << itemp << ","
                    << iZone << "," << iGoodChange << "," << change.total << "," << change.cost << ","
                    << change.connection << "," << change.penalty << "," << settings.temp << "\n";
                logger.AppendDebugTraceFile(debugBuffer.str());
                debugBuffer.clear();
            }

            if (annealTraceType)
            {
                iRowCounter++;
                if (iRowCounter > rowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    string d = annealTraceType > 1 ? "," : " ";

                    tracefp << itime << d << costthresh << d << iGoodChange << d
                            << r.objective.total << d << r.objective.cost << d << r.objective.connection << d
                            << r.objective.penalty << d << r.objective.shortfall << d << ipu << d << settings.temp << "\n"; // itime,costthresh,cost,connection,penalty
                                                                                                                                  /*
#ifdef DEBUG_PEW_CHANGE_PEN
                        AppendDebugTraceFile("iteration,threshold,dochange,total,cost,connection,penalty,puindex\n");
                        sprintf(debugbuffer, "%i,%f,%i,%f,%f,%f,%f,%i\n", itime, costthresh, iGoodChange, reserve->total, reserve->cost, reserve->connection, reserve->penalty, ipu);
                        AppendDebugTraceFile(debugbuffer);
#endif
*/

                    if (suppressAnnealZones == 0)
                    {
                        zonefp << itime;
                        for (int i = 0; i < pu.puno; i++)
                            zonefp << d << zones.IndexToId(r.solution[i]);

                        zonefp <<  "\n";
                    }
                }
            }
        } /* Run Through Annealing */

        if (annealTraceType)
        {
            tracefp.close();
            if (suppressAnnealZones == 0)
            {
                zonefp.close();
            }
        }
    }

    // Randomly chooses the next pu and zone to flip
    // Should only return puindex that is not locked, and a zone that is different from current zone.
    pair<int, int> GetPuAndZone(Reserve& r, Pu& pu, uniform_int_distribution<int>& randomPuDist, uniform_int_distribution<int>& randomDist, int zoneCount) {
        // pick pu at random
        int ipu,iZone, iPreviousR; 

        // pick an index from the list of valid pu.
        ipu = pu.validPuIndices[randomPuDist(rngEngine)];

        iPreviousR = r.solution[ipu];
        iZone = pu.RtnValidZoneForPu(ipu, iPreviousR, randomDist, rngEngine, zoneCount);

        return pair<int,int>(ipu, iZone);
    }

    sanneal settings;

    private:
    void WriteAnnealTraceHeader(ofstream& fp, Reserve& r, double costthresh) {
        string d = annealTraceType > 1 ? "," : " ";
        fp << "iteration" << d << "threshold"<< d << "dochange" << d << "total" << d << "cost" << d << "connection" 
           << d << "penalty" << d << "shortfall" << d << "puindex" << d << "anneal.temp" << "\n";

        // write iteration 0
        fp << 0 << d << costthresh << d << false << d << r.objective.total << d << r.objective.cost << d << r.objective.connection 
            << d << r.objective.penalty << d << r.objective.shortfall << d << -1 << d << settings.temp << "\n";
    }   

    void WriteZoneTraceHeader(ofstream& fp, Reserve& r, Pu& pu) {
        string d = annealTraceType > 1 ? "," : " ";
        fp << "configuration";
        for (int i = 0; i < pu.puno; i++) {
            fp << d << pu.puList[i].id;
        }
        fp << "\n0";
        for (int i = 0; i < pu.puno; i++) {
            fp << d << r.solution[i];
        }
        fp << "\n";
    }

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

    void ConnollyInit(Species& spec, Pu& pu, Zones& zones, int clumptype, double blm) {
        double localdelta = numeric_limits<double>::epsilon();
        uniform_int_distribution<int> random_dist(0, numeric_limits<int>::max());
        uniform_int_distribution<int> random_pu_dist(0, pu.validPuIndices.size()-1);

        // Set reserve to a random and evaluate 
        Reserve r(spec, zones.zoneCount, clumptype);
        r.InitializeSolution(pu.puno);
        r.RandomiseSolution(pu, rngEngine, zones.zoneCount);
        r.EvaluateObjectiveValue(pu, spec, zones, blm);

        schange change = r.InitializeChange(spec, zones);

        int ipu, iPreviousR, chosenZoneInd, chosenZone, imode = 1;
        double deltamin = 0,deltamax = 0;
        for (int i=1;i<= settings.iterations/100; i++)
        {
            // pick pu at random
            pair<int, int> next = GetPuAndZone(r, pu, random_pu_dist, random_dist, zones.zoneCount);
            ipu = next.first;
            chosenZone = next.second;

            // Check change in zone on penalty
            r.CheckChangeValue(change, ipu, r.solution[ipu], chosenZone, pu, zones, spec, 0, blm);

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
    void AdaptiveInit(Species& spec, Pu& pu, Zones& zones, int clumptype, double blm) {
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
            r.EvaluateObjectiveValue(pu, spec, zones, blm);

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

    mt19937 &rngEngine;
    int id;

    // anneal trace settings
    int annealTraceType;
    string savename;
    LoggerBase logger;

    private:
    int saveSnapSteps;
    int saveSnapChanges;
    int saveSnapFrequency; // frequency in terms of iterations of how often to save snaps.
    int suppressAnnealZones;
    int annealingTraceRows;

    // constants for trace storing
    int rowLimit = 0;
};

} // namespace marzone