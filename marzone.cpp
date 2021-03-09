

// C++ file for Marxan with Zones

#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <omp.h>
#include <chrono>

#define MARZONE
#undef MEMDEBUG

#define DEBUGTRACEFILE
#undef ANNEALINGDETAILDEBUGTRACE
#undef EXTRADEBUGTRACE
#undef ANNEALING_TEST
#undef DEBUG_CM  // this definition for debugging run time error in CM>0
#undef DEBUGFPERROR
#undef DEBUGCHECKCHANGE
#undef DEBUGCALCPENALTIES
#undef DEBUGNEWPENALTY
#undef DEBUGSHOWTIMEPASSED
#undef DEBUG_ZONECONNECTIONCALCS
#undef DEBUG_CONNECTIONCOST2

#undef DEBUG_COUNT_MISSING
#undef DEBUG_CONNECTIONCOST2_LINEAR

#undef DEBUGINITIALISERESERVE
#undef DEBUG_PuNotInAllowedZone
#undef DEBUG_CALC_PENALTIES
#undef DEBUG_CHANGE_PEN

#undef DEBUG_PEW_CHANGE_PEN
#undef ASYMCON
#undef PENX_MOD

#undef DEBUG_ZONATION_COST
#undef DEBUG_CONNOLLYINIT

#define DEBUG_PENALTY_NEGATIVE

#include "analysis.hpp"
#include "logger.hpp"
#include "marzone.hpp"

// Solver files
#include "solvers/simulated_annealing.hpp"
#include "solvers/population_annealing.hpp"
#include "solvers/heuristic.hpp"
#include "solvers/iterative_improvement.hpp"

namespace marzone {

// Version specific constants
string sVersionString = "Marxan with Zones v 4.0.5";
string sMarxanWebSite = "https://marxansolutions.org/";
string sDebugTraceFileName = "DebugTraceFile_Marxan_with_Zones.txt";

int iMemoryUsed=0;
int fSpecPROPLoaded = 0;
// version 2.0 introduces these features;
//   enhanced flexibility in objectives
//   probabilistic treatment of threats (1D prob)
//   probabilistic treatment of species (2D prob)
//   asymmetric connectivity
//   multiple connectivity files
int iVerbosity;
int iOptimisationIterativeImprovement = 1, iPuZoneCount = 0;
int iOptimisationCalcPenalties = 1;
int asymmetricconnectivity = 0;
int iZoneContrib3On = 0;
int iZonationCost = 0;
mt19937 rngEngine;
Logger logger;

chrono::high_resolution_clock::time_point startTime;

string StartMessage()
{
    stringstream myfile;
    myfile << "        " << sVersionString << " \n\n  Spatial Prioritization via Zoning and Annealing\n\n";
    myfile << "   Marxan with Zones coded by Matthew Watts\n";
    myfile << "   Written by Ian Ball, Hugh Possingham and Matthew Watts\n\n";
    myfile << "   Based on Marxan coded by Ian Ball, modified by Matthew Watts\n";
    myfile << "   Written by Ian Ball and Hugh Possingham\n\n";
    myfile << "   Marxan website\n\n";
    myfile << sMarxanWebSite << "\n\n";

    return myfile.str();
}

int MarZone(string sInputFileName, int marxanIsSecondary)
{
    // Set start time
    startTime = chrono::high_resolution_clock::now();

    srunoptions runoptions = {};
    sfname fnames = {};
    sanneal anneal = {};

    int iMessageCounter = 0, iZoneSumSolnIndex;
    int seedinit, aggexist=0,sepexist=0;
    int itemp, ipu;
    string tempname;
    double rBestScore;

    ShowStartupScreen();

    SetOptions(sInputFileName, runoptions, anneal, fnames);
    logger.verbosity = runoptions.verbose;
    iVerbosity = runoptions.verbose;
    SetRunOptions(runoptions);

    #ifdef DEBUGTRACEFILE
    logger.StartDebugTraceFile(sDebugTraceFileName);
    logger.AppendDebugTraceFile(sVersionString + " begin execution\n\nLoadOptions\n");
    if (iVerbosity > 3)
       DumpFileNames(fnames, logger);
    #endif

    #ifdef DEBUGCHECKCHANGE
    StartDebugFile("debug_MarZone_CheckChange.csv","ipu,puid,R,total,cost,connection,penalty,threshpen\n",fnames);
    #endif

    #ifdef DEBUG_CHANGE_PEN
    StartDebugFile("debug_MarZone_ChangePen.csv","ipu,puid,isp,spid,cost,newamount,famount\n",fnames);
    #endif

    #ifdef DEBUG_PEW_CHANGE_PEN
    StartDebugFile("debug_MarZone_PewChangePen.csv","iteration,ipu,isp,puid,spid,Zone,newZone,ZoneTarget,PUAmount,Amount,newAmount,Shortfall,newShortfall,rSF,rNSF,iCSF,iNSF,zone\n",fnames);
    #endif

    #ifdef DEBUGCALCPENALTIES
    StartDebugFile("debug_MarZone_CalcPenalties.csv","\n",fnames);
    #endif

    delta = 10e-12; 
    rngEngine = mt19937(runoptions.iseed); // Initialize rng engine

    #ifdef DEBUG
    SaveSeed(runoptions.iseed);
    #endif

    logger.AppendDebugTraceFile("RandSeed iseed " + to_string(runoptions.iseed) + "\n");

    //  ****     Data File Entry    * * * * * * *
    logger.ShowGenProg("\nEntering in the data files \n");
    logger.AppendDebugTraceFile("before Loading Cost Names\n");
    Costs costs(fnames, logger);

    logger.AppendDebugTraceFile("before Loading Zone Files\n");
    Zones zones(fnames, costs, logger); // load all zone files

    // read in the MarZone files
    logger.AppendDebugTraceFile("before Loading Pu Files (Pu lock, pu zone etc.)\n");
    Pu pu(fnames, costs, asymmetricconnectivity, zones.zoneNames, logger);
    logger.ShowDetProg("    Reading in the Planning Unit names \n");
    logger.ShowGenProg("   There are " + to_string(pu.puno) + " Planning units.\n  " + to_string(pu.puno) + " Planning Unit names read in \n");
    logger.ShowDetProg("    Reading in the species file \n");

    Species spec(fnames, logger);
    logger.AppendDebugTraceFile("After Loading species files\n");
    logger.ShowGenProg("  " + to_string(spec.spno) + " species read in \n");

    #ifdef DEBUGTRACEFILE
    if (iVerbosity > 3)
    {
        pu.DumpCostValues(fnames.outputdir + "debugCostValues.csv");
        zones.DumpZoneNames(fnames.outputdir + "debugZoneNames.csv");
        costs.DumpCostNames(fnames.outputdir + "debugCostNames.csv");
        zones.DumpZoneContrib(fnames.outputdir + "debugZoneContrib.csv");
        zones.DumpZoneContrib2(fnames.outputdir + "debugZoneContrib2.csv");
        zones.DumpZoneContrib3(fnames.outputdir + "debugZoneContrib3.csv");
        zones.DumpZoneTarget(fnames.outputdir + "debugZoneTarget.csv");
        zones.DumpZoneTarget2(fnames.outputdir + "debugZoneTarget2.csv");
        zones.DumpZoneCost(fnames.outputdir + "debugZoneCost.csv");

        if (!fnames.pulockname.empty())
            pu.DumpPuLock(fnames.outputdir + "debugPuLock.csv"); 
        if (!fnames.puzonename.empty())
            pu.DumpPuZone(fnames.outputdir + "debugPuZone.csv");
        if (!fnames.relconnectioncostname.empty())
            zones.DumpRelConnectionCost(fnames.outputdir + "debugZoneConnectionCost.csv");
    }
    #endif

    // Build zone contributions. 
    zones.BuildZoneContributions(spec, pu);

    if (fnames.zonetargetname.empty() && fnames.zonetarget2name.empty())
    {
        logger.ShowGenProg("Warning: No targets specified for zones.\n");
        logger.AppendDebugTraceFile("Warning: No targets specified for zones.\n");
    }

    #ifdef DEBUGTRACEFILE
    if (iVerbosity > 3)
    {
        zones.DumpZoneContribFinalValues(fnames.outputdir + "debug_ZoneContrib.csv", spec);
        zones.DumpZoneCostFinalValues(fnames.outputdir + "debug_ZoneCost.csv", costs);
        zones.DumpRelConnectionCostFinalValues(fnames.outputdir + "debug_ZoneConnectionCost.csv");
        pu.DumpPuLockZoneData(fnames.outputdir + "debugPuLockZone.csv");
    }
    #endif

    // Init analysis object
    Analysis analysis;

    if (fnames.savesumsoln)
    {
        logger.AppendDebugTraceFile("before InitSumSoln\n");
        analysis.initSumSolution(pu.puno, zones.zoneCount);
        logger.AppendDebugTraceFile("after InitSumSoln\n");
    }

    logger.ShowGenProg("  " + to_string(pu.connections.size()) + " connections entered \n");
    logger.ShowDetProg("    Reading in the Planning Unit versus Species File \n");
    logger.AppendDebugTraceFile("before LoadSparseMatrix\n");
    pu.LoadSparseMatrix(spec, fnames.inputdir + fnames.puvsprname, logger);
    logger.ShowGenProg(to_string(pu.puvspr.size()) + " conservation values counted, " + 
        to_string(pu.puno*spec.spno) + " big matrix size, " + to_string(pu.density) + "% density of matrix \n");
    logger.AppendDebugTraceFile("after LoadSparseMatrix\n");

    #ifdef DEBUGTRACEFILE
    logger.AppendDebugTraceFile("before CalcTotalAreas\n");
    CalcTotalAreas(pu,spec);
    logger.AppendDebugTraceFile("after CalcTotalAreas\n");
    #endif

    
    if (fnames.savetotalareas)
    {
        tempname = fnames.savename + "_totalareas" + getFileSuffix(fnames.savetotalareas);
        CalcTotalAreas(pu, spec, tempname, true);
    }

    // finalise zone and non-zone targets now that matrix has been loaded
    if (spec.fSpecPROPLoaded)
    {
        logger.AppendDebugTraceFile("before ApplySpecProp\n");

        // species have prop value specified
        ApplySpecProp(spec, pu);

        logger.AppendDebugTraceFile("after ApplySpecProp\n");
    }

    logger.AppendDebugTraceFile("before Build_ZoneTarget\n");
    zones.BuildZoneTarget(spec, pu, fnames, logger);
    logger.AppendDebugTraceFile("after Build_ZoneTarget\n");

    if (iVerbosity > 3)
        zones.DumpZoneTargetFinalValues(fnames.outputdir + "debug_ZoneTarget.csv", spec);

    // Read and process species block definitions
    logger.AppendDebugTraceFile("before process block definitions\n");

    if (!fnames.blockdefname.empty())
    {
        logger.ShowDetProg("    Setting Block Definitions \n");
        vector<double> totalSpecAmount = pu.TotalSpeciesAmount(spec.spno);
        spec.SetSpeciesBlockDefinitions(totalSpecAmount);
    }

    spec.SetSpeciesDefaults();

    logger.AppendDebugTraceFile("after process block definitions\n");
    logger.ShowGenProgInfo("Checking to see if there are aggregating or separating species.\n");

    if (fnames.savesen)
    {
        logger.AppendDebugTraceFile("before OutputScenario\n");

        tempname = fnames.savename + "_sen.dat";
        OutputScenario(pu.puno,spec.spno,zones.zoneCount, costs.costCount, logger, anneal, runoptions, tempname);

        logger.AppendDebugTraceFile("after OutputScenario\n");
    }

    if (runoptions.verbose > 1)
        logger.ShowTimePassed(startTime);

    logger.AppendDebugTraceFile("before initializing reserve\n");
    Reserve reserve(spec, zones.zoneCount, runoptions.clumptype);

   // Init reserve object
    reserve.InitializeSolution(pu.puno);
    logger.AppendDebugTraceFile("after initializing reserve\n");

    logger.AppendDebugTraceFile("before InitialiseReserve\n");
    reserve.RandomiseSolution(pu, rngEngine, zones.zoneCount);
    logger.AppendDebugTraceFile("after InitialiseReserve\n");

    // * * *  Pre-processing    * * * * ***
    logger.ShowGenProg("\nPre-processing Section. \n");
    logger.ShowGenProgInfo("    Calculating all the penalties \n");
  

    if (fnames.penaltyname.empty())
    {
        if (fnames.matrixspordername.empty())
        {
            logger.AppendDebugTraceFile("before CalcPenalties\n");

            // we don't have sporder matrix available, so use slow CalcPenalties method
            itemp = CalcPenalties(pu, spec, zones, reserve, runoptions.clumptype);

            logger.AppendDebugTraceFile("after CalcPenalties\n");
        }
        else
        {
            logger.AppendDebugTraceFile("before CalcPenalties\n");

            itemp = CalcPenalties(pu, spec, zones, reserve, runoptions.clumptype);

            logger.AppendDebugTraceFile("after CalcPenalties\n");
        }

        if (itemp > 0)
            logger.ShowProg(to_string(itemp) + " species cannot meet target%c.\n");
    }
    else
    {
        logger.AppendDebugTraceFile("before LoadPenalties\n");
        spec.LoadCustomPenalties(fnames.inputdir + fnames.penaltyname, logger);
        logger.AppendDebugTraceFile("after LoadPenalties\n");
    }

    logger.ShowTimePassed(startTime);

    if (runoptions.AnnealingOn)
    {
       logger.ShowGenProgInfo("    Calculating temperatures.\n");
       if (!anneal.Titns)
          logger.ShowErrorMessage("Initial Temperature is set to zero. Fatal Error \n");

       anneal.Tlen = anneal.iterations/anneal.Titns;
       logger.ShowGenProgInfo("  Temperature length " + to_string(anneal.Tlen) + " \n");
       logger.ShowGenProgInfo("  iterations " + to_string(anneal.iterations) + ", repeats " + to_string(runoptions.repeats) +" \n");
    } // Annealing Preprocessing. Should be moved to SetAnnealingOptions

    if (fnames.savepenalty)
    {
        string savePenaltyName = fnames.savename + "_penalty" + getFileSuffix(fnames.savepenalty);
        spec.WritePenalties(savePenaltyName, fnames.savepenalty);
    }

    if (fnames.savespeciesdata)
    {
       string saveSpeciesName = fnames.savename + "_spec.csv";
       spec.WriteSpeciesData(saveSpeciesName, reserve.speciesAmounts);
    }

    if (fnames.savesolutionsmatrix)
    {
        string saveSolutionMatrixName = fnames.savename + "_solutionsmatrix" + getFileSuffix(fnames.savesolutionsmatrix);
        pu.WriteSolutionsMatrixHeader(saveSolutionMatrixName,fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);

       for (int i=0;i<zones.zoneCount;i++)
       {
           string saveSolutionMatrixNameByZone = fnames.savename + "_solutionsmatrix_zone" + to_string(zones.IndexToId(i)) + getFileSuffix(fnames.savesolutionsmatrix);
           // init solutions matrix for each zone separately
           pu.WriteSolutionsMatrixHeader(saveSolutionMatrixNameByZone,fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);
       }
    }

    //   The larger repetition loop
    // for each repeat run
    Reserve bestR;
    vector<string> summaries(runoptions.repeats);
    bestR.objective.total = numeric_limits<double>::max();

    // lock for bestR
    omp_lock_t bestR_write_lock;
    omp_init_lock(&bestR_write_lock);

    //create seeds for local rng engines
    vector<unsigned int> seeds(runoptions.repeats);
    for (int run_id = 1; run_id <= runoptions.repeats; run_id++)
        seeds[run_id - 1] = run_id*runoptions.iseed;

    bool quitting_loop = false;
    int maxThreads = omp_get_max_threads();
    logger.ShowGenProg("Running " + to_string(runoptions.repeats) + " runs multithreaded over number of threads: " + to_string(maxThreads) + "\n");
    logger.ShowGenProg("Runs will show as they complete, and may not be in sequential order.\n");
    #pragma omp parallel for schedule(dynamic)
    for (int irun = 1;irun <= runoptions.repeats;irun++)
    {
        if(quitting_loop)
            continue; //skipping iterations. It is not allowed to break or throw out omp for loop.

        mt19937 rngThread(seeds[irun-1]);
        stringstream debugbuffer; // buffer to print at the end
        stringstream progbuffer; // buffer to print at the end
        try
        {
            // Create new reserve object and init
            Reserve reserveThread(spec, zones.zoneCount, runoptions.clumptype, irun);
            reserveThread.InitializeSolution(pu.puno);
            reserveThread.RandomiseSolution(pu, rngThread, zones.zoneCount);

            debugbuffer << "annealing start run " << irun << "\n";
            if (runoptions.verbose > 1)
                progbuffer << "\nRun: " << irun << " ";
          
            SimulatedAnnealing sa(fnames, runoptions.AnnealingOn, anneal,
                rngEngine, fnames.saveannealingtrace, irun);
            if (runoptions.AnnealingOn && !runoptions.PopulationAnnealingOn)
            {
                debugbuffer << "before Annealling Init run " << irun << "\n";

                // init sa parameters if setting is appropriate
                sa.Initialize(spec, pu, zones, runoptions.clumptype, runoptions.blm);
                debugbuffer << "after Annealling Init run " << irun << "\n";

                if (runoptions.verbose > 1)
                    progbuffer << "  Using Calculated Tinit = " << sa.settings.Tinit << "Tcool = " << sa.settings.Tcool << "\n";
            } // Annealing Setup only for sa


            if (runoptions.verbose > 1)
                progbuffer << "  creating the initial reserve \n";

            debugbuffer << "before ZonationCost run " << irun << "\n";
            reserveThread.EvaluateObjectiveValue(pu, spec, zones, runoptions.blm);
            debugbuffer << "after ZonationCost run " << irun << "\n";

            if (runoptions.verbose > 1)
            {
                progbuffer << "\n  Init:";
                PrintResVal(reserveThread, spec, zones, runoptions.misslevel, progbuffer);
            }

            if (runoptions.verbose > 5)
            {
                logger.ShowTimePassed(startTime);
            }

            // * * * * * * * * * * * * * * * * * * * ***
            // * * *  main annealing algorithm * * * * *
            // * * * * * * * * * * * * * * * * * * * ***
            
            if (runoptions.AnnealingOn && !runoptions.PopulationAnnealingOn)
            {
                debugbuffer << "before Annealing run " << irun << "\n";
                if (runoptions.verbose > 1)
                    progbuffer << "  Main Annealing Section.\n";

                sa.RunAnneal(reserveThread, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh, runoptions.blm, logger);

                if (runoptions.verbose > 1)
                {
                    progbuffer << "  ThermalAnnealing:";
                    PrintResVal(reserveThread, spec, zones, runoptions.misslevel, progbuffer);
                }

                debugbuffer << "after Annealing run " << irun << "\n";
            } 
            else if (runoptions.PopulationAnnealingOn)
            {
                // run population annealing instead of regular thermal annealing
                debugbuffer << "before population annealing run " << irun << "\n";
                progbuffer << "  Main Population Annealing Section.\n";

                PopulationAnnealing popAnneal(anneal, rngEngine, irun, fnames);
                popAnneal.Run(reserveThread, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh, runoptions.blm, logger);
                if (runoptions.verbose > 1)
                {
                    progbuffer << "  PopAnnealing:";
                    PrintResVal(reserveThread, spec, zones, runoptions.misslevel, progbuffer);
                }
            }
            
            if (runoptions.HeuristicOn)
            {
                debugbuffer << "before Heuristics run " << irun << "\n";
                Heuristic heur = Heuristic(rngThread, runoptions.heurotype);
                heur.RunHeuristic(reserveThread, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh, runoptions.blm);

                if (runoptions.verbose > 1 && (runoptions.runopts == 2 || runoptions.runopts == 5))
                {
                    progbuffer << "\n  Heuristic:";
                    PrintResVal(reserveThread, spec, zones, runoptions.misslevel, progbuffer);
                }

                debugbuffer << "after Heuristics run " << irun << "\n";
            } // Activate Greedy

            if (runoptions.ItImpOn)
            {
                debugbuffer << "before IterativeImprovementOptimise run " << irun << "\n";
                IterativeImprovement itImp = IterativeImprovement(rngThread, fnames, runoptions.itimptype);
                itImp.Run(reserveThread, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh, runoptions.blm);
                debugbuffer << "after IterativeImprovementOptimise run " << irun << "\n";

                if (runoptions.verbose > 1)
                {
                    progbuffer << "  Iterative Improvement:";
                    PrintResVal(reserveThread, spec, zones, runoptions.misslevel, progbuffer);
                }
            } // Activate Iterative Improvement

            debugbuffer << "before file output run " << irun << "\n";

            // Start writing output
            string tempname2;
            string paddedRun = intToPaddedString(irun, 5);
            if (fnames.saverun)
            {
                tempname2 = fnames.savename + "_r" + paddedRun + getFileSuffix(fnames.saverun);
                reserveThread.WriteSolution(tempname2, pu, zones, fnames.saverun);
            }

            debugbuffer << "WriteSolution ran " << irun << "\n";

            if (fnames.savespecies && fnames.saverun)
            {
                // output species distribution for a run (specific reserve)
                tempname2 = fnames.savename + "_mv" + paddedRun + getFileSuffix(fnames.savespecies);
                OutputFeatures(tempname2, zones, reserveThread, spec, fnames.savespecies, runoptions.misslevel);
            }

            if (fnames.savesum)
            {
                summaries[irun - 1] = OutputSummaryString(pu, spec, zones, reserveThread, runoptions.misslevel, fnames.savesum, runoptions.blm);
            }

#ifdef DEBUGFPERROR
            debugbuffer << "OutputSummary ran\n";
#endif

            // Saving the best from all the runs
            if (fnames.savebest)
            {
                // re-evaluate entire system in case of any floating point/evaluation errors.
                reserveThread.EvaluateObjectiveValue(pu, spec, zones, runoptions.blm);
                if (reserveThread.objective.total < bestR.objective.total)
                {
                    omp_set_lock(&bestR_write_lock);
                    if (reserveThread.objective.total < bestR.objective.total)
                    {
                        bestR = reserveThread; // deep copy

                        if (runoptions.verbose > 1)
                        {
                            progbuffer << "  Best:";
                            PrintResVal(bestR, spec, zones, runoptions.misslevel, progbuffer);
                        }
                    }
                    omp_unset_lock(&bestR_write_lock);
                }
            }

            if (fnames.savesumsoln) // Add current run to my summed solution
                analysis.ApplyReserveToSumSoln(reserveThread);

            if (fnames.savesolutionsmatrix)
            {
                string solutionsMatrixName = fnames.savename + "_solutionsmatrix" + getFileSuffix(fnames.savesolutionsmatrix);
                reserveThread.AppendSolutionsMatrix(solutionsMatrixName, zones.zoneCount, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);

                for (int i = 0; i < zones.zoneCount; i++)
                {
                    string solutionsMatrixZoneName = fnames.savename + "_solutionsmatrix_zone" + to_string(zones.IndexToId(i)) + getFileSuffix(fnames.savesolutionsmatrix);
                    // append solutions matrix for each zone separately
                    reserveThread.AppendSolutionsMatrixZone(solutionsMatrixZoneName, i, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);
                }
            }

            if (fnames.savezoneconnectivitysum)
            {
                string zoneConnectivityName = fnames.savename + "_zoneconnectivitysum" + paddedRun + getFileSuffix(fnames.savezoneconnectivitysum);
                reserveThread.WriteZoneConnectivitySum(zoneConnectivityName, pu, zones, fnames.savezoneconnectivitysum);
            }

            debugbuffer << "after file output run " << irun << "\n";
            debugbuffer << "annealing end run " << irun << "\n";

            if (marxanIsSecondary == 1)
                WriteSecondarySyncFileRun(irun);
        }
        catch (const exception &e)
        {
            // No matter how thread runs, record the message and exception if any.
            string msg = "Exception on run " + to_string(irun) + " Message: " + e.what() + "\n";
            progbuffer << msg;
            debugbuffer << msg;

            //cannot throw or break out of omp loop
            quitting_loop = true;
            continue;
        }

        // print the logs/debugs
        logger.ShowProg(progbuffer.str());
        logger.AppendDebugTraceFile(debugbuffer.str());

        if (runoptions.verbose > 1)
            logger.ShowTimePassed(startTime);
    } // ** the repeats  **

    if (quitting_loop)
    {
        logger.ShowErrorMessage("\nRuns were aborted due to error.\n"); 
    }

    // Output summary
    if (fnames.savesum)
    {
        string tempname2 = fnames.savename + "_sum" + getFileSuffix(fnames.savesum);
        OutputSummary(pu, zones, summaries, tempname2, fnames.savesum);
    }

    logger.AppendDebugTraceFile("before final file output\n");

    if (fnames.savebest)
    {
        string saveBestName = fnames.savename + "_best" + getFileSuffix(fnames.savebest);
        bestR.WriteSolution(saveBestName, pu, zones, fnames.savebest);

        logger.AppendDebugTraceFile("Best solution is run " + to_string(bestR.id) + "\n");
        logger.ShowProg("\nBest solution is run " + to_string(bestR.id) + "\n");

        // Print the best solution
        stringstream bestbuf; 
        PrintResVal(bestR, spec, zones, runoptions.misslevel, bestbuf);
        logger.ShowWarningMessage(bestbuf.str());

        if (fnames.savezoneconnectivitysum)
        {
            string saveZoneConnectivityName = fnames.savename + "_zoneconnectivitysumbest" + getFileSuffix(fnames.savezoneconnectivitysum);
            bestR.WriteZoneConnectivitySum(saveZoneConnectivityName, pu, zones, fnames.savezoneconnectivitysum);
        }

        if (fnames.savespecies)
        {
            string saveBestSpeciesName = fnames.savename + "_mvbest" + getFileSuffix(fnames.savespecies);
            OutputFeatures(saveBestSpeciesName,zones, bestR, spec,fnames.savespecies, runoptions.misslevel);
        }
    }

    if (fnames.savesumsoln)
    {
        string saveSumSolnName = fnames.savename + "_ssoln" + getFileSuffix(fnames.savesumsoln);
        analysis.WriteAllSumSoln(saveSumSolnName, pu, zones, fnames.savesumsoln);
    }

    ShowShutdownScreen();

    logger.CloseLogFile();

    #ifdef DEBUGTRACEFILE
    logger.AppendDebugTraceFile("end final file output\n");
    logger.AppendDebugTraceFile("\nMarxan with Zones end execution\n");
    #endif

    return 0;

} // MarZone

// *** Set run options. Takes an integer runopts value and returns flags ***
void SetRunOptions(srunoptions& runoptions)
{
    if (runoptions.runopts < 0)
        return; // runopts < 0 indicates that these are set in some other way
    switch (runoptions.runopts)
    {
    case 0:
        runoptions.AnnealingOn = 1;
        runoptions.HeuristicOn = 1;
        runoptions.ItImpOn = 0;
        break;
    case 1:
        runoptions.AnnealingOn = 1;
        runoptions.HeuristicOn = 0;
        runoptions.ItImpOn = 1;
        break;
    case 2:
        runoptions.AnnealingOn = 1;
        runoptions.HeuristicOn = 1;
        runoptions.ItImpOn = 1;
        break;
    case 3:
        runoptions.AnnealingOn = 0;
        runoptions.HeuristicOn = 1;
        runoptions.ItImpOn = 0;
        break;
    case 4:
        runoptions.AnnealingOn = 0;
        runoptions.HeuristicOn = 0;
        runoptions.ItImpOn = 1;
        break;
    case 5:
        runoptions.AnnealingOn = 0;
        runoptions.HeuristicOn = 1;
        runoptions.ItImpOn = 1;
        break;
    case 6:
        runoptions.AnnealingOn = 1;
        runoptions.HeuristicOn = 0;
        runoptions.ItImpOn = 0;
        break;
    default:
        runoptions.AnnealingOn = 0;
        runoptions.HeuristicOn = 0;
        runoptions.ItImpOn = 0;
        break;
    }

} // Set Run Options

// * * * * Calculate Initial Penalties * * * *
// This routine calculates the initial penalties or the penalty if you had no representation
int CalcPenalties(Pu& pu, Species& spec, Zones& zones, Reserve& r, int clumptype) {
    int badspecies = 0, goodspecies = 0, itargetocc;
    double rZoneSumTarg, iZoneSumOcc, penalty, ftarget, rAmount, ftemp;

    vector<double> specTargetZones = zones.AggregateTargetAreaBySpecies(spec.spno);
    vector<int> specOccurrenceZones = zones.AggregateTargetOccurrenceBySpecies(spec.spno);

    vector<vector<lockedPenaltyTerm>> lockedSpecAmounts;
    vector<vector<penaltyTerm>> specPuAmounts = pu.getPuAmountsSorted(spec.spno, lockedSpecAmounts); // list of pus that contribute to a species.

    for (int i = 0; i < spec.spno; i++)
    {
        rZoneSumTarg = specTargetZones[i];
        iZoneSumOcc = specOccurrenceZones[i];

        sspecies &specTerm = spec.specList[i];

        if (specTerm.target > rZoneSumTarg)
            rZoneSumTarg = specTerm.target;
        if (specTerm.targetocc > iZoneSumOcc)
            iZoneSumOcc = specTerm.targetocc;

        if (specTerm.target2)
        {
            int j = r.ComputePenaltyType4(specTerm, specPuAmounts[i], i, rZoneSumTarg, specTerm.target2, iZoneSumOcc); 
            badspecies += (j > 0);
            goodspecies += (j < 0);

            continue;
        } // Species has aggregation requirements

        itargetocc = 0, ftarget = 0.0, penalty = 0.0;
        int lockedZone = -1;
        double lockedContrib = 0.0;
        // For this species, sum up all locked occurrences and areas.
        for (lockedPenaltyTerm& term: lockedSpecAmounts[i]) {
            lockedZone = term.lockedZoneId;
            //lockedContrib = zones.GetZoneContrib(i, lockedZone);

            //if (lockedContrib) {
                // Do this calculation by taking into account zonecontrib of the locked zone. If positive contrib, we count it.
            //    rAmount = term.amount*lockedContrib;
            rAmount = term.amount;
                if (rAmount > 0) {
                    ftarget += rAmount;
                    itargetocc++;
                    penalty += rtnMaxNonAvailableCost(term.puindex, pu, zones);
                }
            //}
        }
        spec.specList[i].penalty = penalty;

        // Already adequately represented on type 2 planning unit
        if (ftarget >= rZoneSumTarg && itargetocc >= iZoneSumOcc)
        {
            goodspecies++;
            logger.ShowGenProgInfo("Species " + to_string(spec.specList[i].name) +"(" +
                spec.specList[i].sname + ") has already met target " + to_string(rZoneSumTarg) +"\n");

            continue;
        }

        // cycle through all pus that contain this spec, and keep adding until spec target is met
        bool targetMet = false;
        for (penaltyTerm& p: specPuAmounts[i]) {
            if (p.amount) {
                ftarget += p.amount;
                itargetocc++;
                spec.specList[i].penalty += p.cost;
            }

            // Check if targets met
            if (ftarget >= rZoneSumTarg && itargetocc >= iZoneSumOcc) {
                targetMet = true;
                break;
            }
        }

        // If target not met with available pu, scale the penalty. 
        if (!targetMet) {
            logger.ShowGenProgInfo("Species " + to_string(spec.specList[i].name) + "(" +
                spec.specList[i].sname + ") cannot reach target " + to_string(rZoneSumTarg) + " there is only " + to_string(ftarget) + " available.\n");

            if (ftarget == 0)
                ftarget = delta; // Protect against divide by zero
            ftemp = 0;
            if (ftarget < rZoneSumTarg)
                ftemp = rZoneSumTarg / ftarget;
            if (itargetocc < iZoneSumOcc && itargetocc) // If ! itargetocc then also !ftarget
                ftemp += (double)iZoneSumOcc / (double)itargetocc;
            spec.specList[i].penalty = spec.specList[i].penalty * ftemp; // Scale it up
            // This value will be ~ 1/delta when there are no occ's of target species in system
            badspecies++;
        }
    }

    if (goodspecies)
        logger.ShowGenProg(to_string(goodspecies) + " species are already adequately represented.\n");

    logger.AppendDebugTraceFile("CalcPenalties end\n");

    return(badspecies);
}

double rtnMaxNonAvailableCost(int ipu, Pu& pu, Zones& zones)
{
    double fcost = 0, rMaxCost = 0;

    if (zones.availableZoneCost) // only needed if there is zone cost specified.
        for (int iZone = 0; iZone < zones.zoneCount; iZone++)
        {
            fcost = ReturnPuZoneCost(ipu, iZone, pu, zones);
            if (fcost > rMaxCost)
                rMaxCost = fcost;
        }
    else
        rMaxCost = pu.puList[ipu].cost;

    rMaxCost += pu.ConnectionCost1(ipu);
    return(rMaxCost);
} // Cost of Planning Unit

// return cost of planning unit in given zone
// parameter iZone is zero base
double ReturnPuZoneCost(int ipu,int iZone, Pu& pu, Zones& zones)
{
    return zones.AggregateTotalCostByPuAndZone(iZone, pu.puList[ipu].costBreakdown);
}

// * * * * * * * * * * * * * * * * * * * * * * * * *****
// * * * * *** Post Processing * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * *****

// * * * * Reporting Value of a Reserve * * * *
void PrintResVal(Reserve& reserve, Species& spec, Zones& zones, double misslevel, stringstream& buffer)
{
    int iMissing;
    double rMPM;

    iMissing = reserve.CountMissing(spec, zones, misslevel, rMPM);
    string sPuZones = reserve.CountPuZones(zones);

    // Attach message to buffer
    buffer << "Value " << reserve.objective.total << " "
        << "Cost " <<  reserve.objective.cost << " " << sPuZones << " "
        << "Connection " << reserve.objective.connection << " "
        << "Missing " << iMissing << " "
        << "Shortfall " << reserve.objective.shortfall << " "
        << "Penalty " << reserve.objective.penalty << " "
        << "MPM " << rMPM << "\n";
} /* * * * Print Reserve Value * * * * */

/*
    Output functions
*/

/* * * * ***** Output Solutions * * * * * * * */
/** imode = 1   Output Summary Stats only ******/
/** imode = 2    Output Everything * * * * *****/
void OutputSummary(Pu& pu, Zones& zones, vector<string>& summaries, string filename, int imode) {
    ofstream fp;  /* Imode = 1, REST output, Imode = 2, Arcview output */
    fp.open(filename);

    if (!fp.is_open())
        logger.ShowErrorMessage("Cannot save output to " + filename + " \n");

    string sZoneNames = zones.ZoneNameHeaders(imode, " PuCount");
    string sZoneCostNames = zones.ZoneNameHeaders(imode, " Cost");

    if (imode > 1)
    {
        fp << "\"Run Number\",\"Score\",\"Cost\",\"Planning Units\"" << sZoneNames << sZoneCostNames;
        fp << ",\"Connection Strength\",\"Penalty\",\"Shortfall\",\"Missing_Values\",\"MPM\"\n";
    }
    else
    {
        fp << "Run no.    Score      Cost   Planning Units  " << sZoneNames << sZoneCostNames;
        fp << "  Connection_Strength   Penalty  Shortfall Missing_Values MPM\n";
    }

    // write all summaries in run order
    for (string& s: summaries) {
        fp << s;
    }

    fp.close();
}

// formats a summary string for a particular run.
string OutputSummaryString(Pu& pu, Species& spec, Zones& zones, Reserve& r, double misslevel, int imode, double blm)
{
    int ino=0,isp;
    double connectiontemp = 0,rMPM;
    stringstream s;
    string d = imode > 1 ? "," : "    ";

    string sZoneNames,sZonePuCount,sZoneCostNames,sZoneCost;
    r.CountPuZones2(zones, imode, sZonePuCount);
    r.CostPuZones(pu, zones, sZoneCost, imode);

    /*** Ouput the Summary Statistics *****/
    for (int i=0;i<pu.puno;i++)
        if (r.solution[i] < 0)
            ino ++;

    isp = r.CountMissing(spec, zones, misslevel, rMPM);
    for (int i=0;i<pu.puno;i++)
    {
        connectiontemp += zones.ConnectionCost2Linear(pu, i, 1, r.solution, blm);
    } /* Find True (non modified) connection */

    s << r.id << d << r.objective.total << d << r.objective.cost << d << ino << sZonePuCount 
      << sZoneCost << d << connectiontemp << d << r.objective.penalty << d << r.objective.shortfall << d << isp << d << rMPM << "\n";

    return s.str();
} // OutputSummary

/* * * * * * * * * * * * * * * * * * * * ****/
/* * * *   Main functions * * * * * * * * ***/
/* * * * * * * * * * * * * * * * * * * * ****/

/* The following function shows the startup screen information.
   The program title and authors */

void ShowStartupScreen(void)
{
    logger.ShowProg(StartMessage());
}  /* Show Startup Screen */


/* Show ShutDown Screen displays all the end of program information. It only displays when
    the iVerbosity has been set to 1 or higher */
void ShowShutdownScreen(void)
{
        logger.ShowProg("\n");
        logger.ShowTimePassed(startTime);
        logger.ShowProg("\n              The End. \n");
}

void SaveSeed(int iseed)
{
     ofstream fp;
     fp.open("debug.out");
     fp << "Debugging Output! \n";
     fp << "iseed is " << iseed << " \n";
     fp.close();
}

/* ShowPauseExit delivers a message prior to exiting */
void ShowPauseExit(void)
{
     logger.ShowProg("Press return to exit.\n");
     std::cin.get();
}  /* Show Pause Exit  */

void WriteSecondarySyncFileRun(int iSyncRun)
{
     FILE* fsync;
     char sSyncFileName[80];

     sprintf(sSyncFileName,"sync%i",iSyncRun);

     fsync = fopen(sSyncFileName,"w");
     fprintf(fsync,"%s",sSyncFileName);
     fclose(fsync);
}

/* SecondaryExit does not deliver a message prior to exiting, but creates a file so C-Plan knows marxan has exited */
void SecondaryExit(void)
{
    WriteSecondarySyncFileRun(0);
}

void ShowPauseProg(void)
{
     logger.ShowProg("Press return to continue.\n");
     getchar();
} /** Pause **/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ***/
/*        Set Options    */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ***/
void SetOptions(string &sInputFileName, srunoptions &runoptions, sanneal &anneal, sfname &fnames)
{
    double version = 0.1;
    bool present;
    stringstream errorMessage;
    stringstream warningMessage;
    string stemp;

    logger.verbosity = 1; /* This enables local warning messages */
    SetDefaultOptions(runoptions, anneal, fnames);

    /* Open file and then feed in each variable type */
    vector<string> lines = GetFileLines(sInputFileName, errorMessage);
    if (!errorMessage.str().empty())
        logger.ShowErrorMessage(errorMessage.str());

    readInputOption(lines,"VERSION",version, false, present, warningMessage, errorMessage);
    readInputOption(lines,"PROP",runoptions.prop, false, present, warningMessage, errorMessage);
    readInputOption(lines,"RANDSEED",runoptions.iseed, false, present, warningMessage, errorMessage); /* The random seed. -1 to set by clock */
    if (!present || runoptions.iseed == -1) { //if seed not present or -1, set as time based.
        runoptions.iseed = (long int)time(NULL);
    }
    readInputOption(lines,"BLM",runoptions.blm, false, present, warningMessage, errorMessage);

    /* Annealing Controls */
    readInputOption(lines, "NUMITNS", anneal.iterations, false, present, warningMessage, errorMessage);
    readInputOption(lines, "STARTTEMP", anneal.Tinit, false, present, warningMessage, errorMessage);
    readInputOption(lines, "COOLFAC", anneal.Tcool, false, present, warningMessage, errorMessage);
    readInputOption(lines, "NUMTEMP", anneal.Titns, false, present, warningMessage, errorMessage);

    anneal.type = 1;
    if (anneal.iterations < 1)
        anneal.type = 0;
    if (anneal.Tinit < 0)
        anneal.type = (int)(-anneal.Tinit) + 1; /* type is negative of Tinit */

    /* Various controls */
    readInputOption(lines, "NUMREPS", runoptions.repeats, false, present, warningMessage, errorMessage);
    readInputOption(lines, "COSTTHRESH", runoptions.costthresh, false, present, warningMessage, errorMessage);
    readInputOption(lines, "THRESHPEN1", runoptions.tpf1, false, present, warningMessage, errorMessage);
    readInputOption(lines, "THRESHPEN2", runoptions.tpf2, false, present, warningMessage, errorMessage);

    /* SaveFiles */
    readInputOption(lines, "SCENNAME", fnames.savename, false, present, warningMessage, errorMessage);
    /* SaveFiles New Method */
    readInputOption(lines, "SAVERUN", fnames.saverun, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVEBEST", fnames.savebest, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESUMMARY", fnames.savesum, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESCEN", fnames.savesen, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVETARGMET", fnames.savespecies, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESUMSOLN", fnames.savesumsoln, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESPECIESDATA", fnames.savespeciesdata, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVEPENALTY", fnames.savepenalty, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVETOTALAREAS", fnames.savetotalareas, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESOLUTIONSMATRIX", fnames.savesolutionsmatrix, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SOLUTIONSMATRIXHEADERS", fnames.solutionsmatrixheaders, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVELOG", fnames.savelog, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESNAPSTEPS", fnames.savesnapsteps, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESNAPCHANGES", fnames.savesnapchanges, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVESNAPFREQUENCY", fnames.savesnapfrequency, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVEANNEALINGTRACE", fnames.saveannealingtrace, false, present, warningMessage, errorMessage);
    readInputOption(lines, "ANNEALINGTRACEROWS", fnames.annealingtracerows, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SUPPRESSANNEALZONES", fnames.suppressannealzones, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVEITIMPTRACE", fnames.saveitimptrace, false, present, warningMessage, errorMessage);
    readInputOption(lines, "SAVEZONECONNECTIVITYSUM", fnames.savezoneconnectivitysum, false, present, warningMessage, errorMessage);
    readInputOption(lines, "ITIMPTRACEROWS", fnames.itimptracerows, false, present, warningMessage, errorMessage);
    if (!fnames.savesnapfrequency)
        fnames.savesnapfrequency = 1;

    /* Filenames */
    readInputOption(lines, "INPUTDIR", fnames.inputdir, true, present, warningMessage, errorMessage);
    fnames.inputdir = cleanDirectoryString(fnames.inputdir);

    readInputOption(lines, "OUTPUTDIR", fnames.outputdir, true, present, warningMessage, errorMessage);
    fnames.outputdir = cleanDirectoryString(fnames.outputdir);

    if (!fnames.outputdir.empty())
    {
        fnames.savename = fnames.outputdir + fnames.savename;
    }

    if (fnames.savelog)
    {
        logger.SetLogFile(fnames.savename + "_log.dat");
        logger.AppendLogFile(StartMessage());
    }

    readInputOption(lines, "PUNAME", fnames.puname, true, present, warningMessage, errorMessage);

    stemp = "spec.dat";
    readInputOption(lines, "SPECNAME", stemp, false, present, warningMessage, errorMessage);
    if (!present) // SPECNAME not present
    {
        readInputOption(lines, "FEATNAME", stemp, false, present, warningMessage, errorMessage);
    }
    fnames.specname = stemp;

    readInputOption(lines, "PENALTYNAME", fnames.penaltyname, false, present, warningMessage, errorMessage);

    fnames.puvsprname = "puvspr2.dat";
    readInputOption(lines, "PUVSPRNAME", fnames.puvsprname, false, present, warningMessage, errorMessage);
    readInputOption(lines, "MATRIXSPORDERNAME", fnames.matrixspordername, false, present, warningMessage, errorMessage);
    if (present)
        logger.ShowWarningMessage("input.dat option: MATRIXSPORDERNAME no longer needed and will be ignored. Please refer to the version 4 feature changelog.md");

    readInputOption(lines, "BOUNDNAME", fnames.connectionname, false, present, warningMessage, errorMessage);
    if (!present)
        readInputOption(lines, "CONNECTIONNAME", fnames.connectionname, false, present, warningMessage, errorMessage);

    readInputOption(lines, "BLOCKDEFNAME", fnames.blockdefname, false, present, warningMessage, errorMessage);
    readInputOption(lines, "ZONESNAME", fnames.zonesname, false, present, warningMessage, errorMessage);
    readInputOption(lines, "COSTSNAME", fnames.costsname, false, present, warningMessage, errorMessage);

    // Read zone contrib settings. Only one of these files should be defined.
    int numZoneContribsDefined = 0;
    readInputOption(lines, "ZONECONTRIBNAME", fnames.zonecontribname, false, present, warningMessage, errorMessage);
    if (present) numZoneContribsDefined++;
    readInputOption(lines, "ZONECONTRIB2NAME", fnames.zonecontrib2name, false, present, warningMessage, errorMessage);
    if (present) numZoneContribsDefined++;
    readInputOption(lines, "ZONECONTRIB3NAME", fnames.zonecontrib3name, false, present, warningMessage, errorMessage);
    if (present) numZoneContribsDefined++;
    
    if (numZoneContribsDefined > 1 )
        logger.ShowErrorMessage("Multiple Zone Contribution files defined. Please only define one ZONECONTRIBNAME, ZONECONTRIB2NAME or ZONECONTRIB3NAME");
    
    readInputOption(lines, "ZONETARGETNAME", fnames.zonetargetname, false, present, warningMessage, errorMessage);
    readInputOption(lines, "ZONETARGET2NAME", fnames.zonetarget2name, false, present, warningMessage, errorMessage);   
    readInputOption(lines, "ZONECOSTNAME", fnames.zonecostname, false, present, warningMessage, errorMessage);
    readInputOption(lines, "PULOCKNAME", fnames.pulockname, false, present, warningMessage, errorMessage);    
    readInputOption(lines, "PUZONENAME", fnames.puzonename, false, present, warningMessage, errorMessage);

    readInputOption(lines, "ZONEBOUNDCOSTNAME", fnames.relconnectioncostname, false, present, warningMessage, errorMessage);
    if (!present)
        readInputOption(lines, "ZONECONNECTIONCOSTNAME", fnames.relconnectioncostname, false, present, warningMessage, errorMessage);

    /* various other controls */
    readInputOption(lines, "RUNMODE", runoptions.runopts, true, present, warningMessage, errorMessage);
    readInputOption(lines, "MISSLEVEL", runoptions.misslevel, false, present, warningMessage, errorMessage);
    readInputOption(lines, "HEURTYPE", runoptions.heurotype, false, present, warningMessage, errorMessage);
    readInputOption(lines, "CLUMPTYPE", runoptions.clumptype, false, present, warningMessage, errorMessage);
    readInputOption(lines, "ITIMPTYPE", runoptions.itimptype, false, present, warningMessage, errorMessage);
    readInputOption(lines, "VERBOSITY", runoptions.verbose, false, present, warningMessage, errorMessage);
    readInputOption(lines, "POPULATIONANNEALINGON", runoptions.PopulationAnnealingOn, false, present, warningMessage, errorMessage);

    
    // Check and print any warning/error messages
    logger.ShowWarningMessage(warningMessage.str());
    if (!errorMessage.str().empty())
        logger.ShowErrorMessage(errorMessage.str());

} /***** Set Options 2* * * */

// use the prop value from the conservation feature file to set a proportion target for species
void ApplySpecProp(Species& spec, Pu& pu)
{
    vector<double> speciesSums = pu.TotalSpeciesAmount(spec.spno);
    spec.SetSpeciesProportionTarget(speciesSums);
}

void CalcTotalAreas(Pu& pu, Species& spec, string filename /*= "MarZoneTotalAreas.csv"*/, bool save /*= false*/)
{
    int ipu, i, ism, isp;
    vector<int> TotalOccurrences, TO_2, TO_3;
    vector<double> TA_2, TA_3, TotalAreas = pu.TotalSpeciesAmount(spec.spno);

    // Set total areas in species.
    spec.SetTotalAreas(TotalAreas);

    if (iVerbosity > 3 || save) 
    {
        TotalOccurrences = pu.TotalOccurrenceAmount(spec.spno);
        TO_2.resize(spec.spno, 0);
        TO_3.resize(spec.spno, 0);
        TA_2.resize(spec.spno, 0.0);
        TA_3.resize(spec.spno, 0.0);

        pu.TotalSpeciesAmountByStatus(TO_2, TO_3, TA_2, TA_3);

        // Write calculated arrays
        spec.WriteTotalAreasAndOccs(filename, TotalOccurrences, TO_2, TO_3, TA_2, TA_3);
    }
} // CalcTotalAreas

void StartDebugFile(string sFileName,string sHeader, sfname& fnames)
{
    string writename = fnames.outputdir + sFileName;
    ofstream myfile;
    myfile.open(writename);
    myfile << sHeader;
    myfile.close();
}

void AppendDebugFile(string sFileName,string& sLine, sfname& fnames)
{
    string writename = fnames.outputdir + sFileName;
    ofstream myfile;
    myfile.open(writename, ofstream::app);
    myfile << sLine;
    myfile.close();
}

void DumpFileNames(sfname& fnames, Logger& logger)
{
    FILE *fp;
    string writename;

    writename = fnames.outputdir + "debugFileNames.csv";
    fp = fopen(writename.c_str(),"w");
    if (fp==NULL)
        logger.ShowErrorMessage("cannot create DumpFileNames file " + writename + "\n");

    fprintf(fp,"input name,file name\n");

    fprintf(fp,"zonesname,%s\n",fnames.zonesname.c_str());
    fprintf(fp,"costsname,%s\n",fnames.costsname.c_str());
    fprintf(fp,"zonecontribname,%s\n",fnames.zonecontribname.c_str());
    fprintf(fp,"zonecontrib2name,%s\n",fnames.zonecontrib2name.c_str());
    fprintf(fp,"zonecontrib3name,%s\n",fnames.zonecontrib3name.c_str());
    fprintf(fp,"zonetargetname,%s\n",fnames.zonetargetname.c_str());
    fprintf(fp,"zonetarget2name,%s\n",fnames.zonetarget2name.c_str());
    fprintf(fp,"zonecostname,%s\n",fnames.zonecostname.c_str());
    fprintf(fp,"pulockname,%s\n",fnames.pulockname.c_str());
    fprintf(fp,"puzonename,%s\n",fnames.puzonename.c_str());
    fprintf(fp,"zoneconnectioncostname,%s\n",fnames.relconnectioncostname.c_str());

    fclose(fp);
}

void OutputScenario(int puno,int spno, int zoneCount, int costCount, Logger& logger,
                    sanneal& anneal, srunoptions& runoptions,
                    string filename)
{
    ofstream fp;
    string temp;
    fp.open(filename);
    if (!fp.is_open())
    {
        logger.ShowErrorMessage("Error: Cannot save to log file " + filename + "\n");
    } /* open failed */

    fp << "Number of Planning Units " << puno << "\n";
    fp << "Number of Conservation Values " << spno << "\n";
    fp << "Number of Zones " << zoneCount << "\n";
    fp << "Number of Costs " << costCount << "\n";
    fp << "Starting proportion " << runoptions.prop << "\n";
    switch (runoptions.clumptype)
    {
        case 0: temp = "Clumping - default step function\n";break;
        case 1: temp = "Clumping - two level step function.\n";break;
        case 2: temp = "Clumping - rising benefit function\n";break;
    }
    fp << temp;

    /* Use character array here and set up the name of the algorithm used */
    switch (runoptions.runopts)
    {
        case 0: temp="Annealing and Heuristic";break;
        case 1: temp="Annealing and Iterative Improvement";break;
        case 2: temp="Annealing and Both";break;
        case 3: temp="Heuristic only";break;
        case 4: temp="Iterative Improvement only";break;
        case 5: temp="Heuristic and Iterative Improvement";
    }
    fp << "Algorithm Used :" << temp << "\n";
    if (runoptions.runopts == 0 || runoptions.runopts == 3 || runoptions.runopts == 5)
    {
        switch (runoptions.heurotype)
        {
            case 0: temp="Richness";break;
            case 1: temp="Greedy";break;
            case 2: temp="Maximum Rarity";break;
            case 3: temp="Best Rarity";break;
            case 4: temp="Average Rarity";break;
            case 5: temp="Summation Rarity";break;
            case 6: temp="Product Irreplaceability";break;
            case 7: temp="Summation Irreplaceability";break;
            default: temp="Unkown Heuristic Type";
        }
        fp << "Heuristic type : " << temp << "\n";
    }
    else
        fp << "No Heuristic used \n";

    if (runoptions.runopts <= 2)
    {
        fp << "Number of iterations " << anneal.iterations << "\n";
        if (anneal.Tinit >= 0)
        {
            fp << "Initial temperature " << anneal.Tinit << "\n";
            fp << "Cooling factor " << anneal.Tcool << "\n";
        }
        else
        {
            fp << "Initial temperature set adaptively" << "\n";
            fp << "Cooling factor set adaptively" << "\n";
        }
        fp << "Number of temperature decreases " << anneal.Titns << "\n\n";
    }
    else
    {
        fp << "Number of iterations N/A\nInitial temperature N/A\nCooling Factor N/A\n";
        fp << "Number of temperature decreases N/A\n\n";
    }

    if (runoptions.costthresh)
    {
        fp << "Cost Threshold Enabled: " << runoptions.costthresh << "\n";
        fp << "Threshold penalty factor A " << runoptions.tpf1 << "\n";
        fp << "Threshold penalty factor B " << runoptions.tpf2 << "\n\n";
    }
    else
    {
        fp << "Cost Threshold Disabled\nThreshold penalty factor A N/A\n";
        fp << "Threshold penalty factor B N/A\n\n";
    }

    fp << "Random Seed " << runoptions.iseed << "\n";
    fp << "Number of runs " << runoptions.repeats << "\n";
    fp.close();
}  /*** OutputScenario ****/

// Output Feature report (missing values report)
void OutputFeatures(std::string filename, marzone::Zones& zones, marzone::Reserve& reserve, marzone::Species& spec, int imode, double misslevel)
{
    zones.WriteZoneTargetHeaders(filename, imode);
    reserve.WriteSpeciesAmounts(filename,spec, zones, imode, misslevel);
} // OutputFeatures

} // namespace marzone

int main(int argc,char *argv[])
{
    std::string sInputFileName = "input.dat"; // If no arguments then assume the default file name of 'input.dat'
    int marxanIsSecondary = 0;

    if (argc > 1)
        marzone::HandleOptions(argc,argv,sInputFileName,marxanIsSecondary);  // handle the program options

    // Calls the main annealing unit
    try {
        if (marzone::MarZone(sInputFileName, marxanIsSecondary))
        {
            if (marxanIsSecondary == 1)
                marzone::SecondaryExit();
            return 1;
        } // Abnormal Exit
        if (marxanIsSecondary == 1)
            marzone::SecondaryExit();
    } 
    catch (const std::exception& e)
    {
        marzone::logger.ShowErrorMessage("Error occurred in execution: " + std::string(e.what()));
        exit(EXIT_FAILURE);
    }

    return 0;
}

