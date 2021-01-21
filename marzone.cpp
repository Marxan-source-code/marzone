

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
#include "marzone.hpp"
#include "input.hpp"
#include "output.hpp"
#include "util.hpp"
#include "costs.hpp"
#include "pu.hpp"
#include "reserve.hpp"
#include "species.hpp"
#include "zones.hpp"

// Solver filers
#include "solvers/simulated_annealing.hpp"
#include "solvers/heuristic.hpp"
#include "solvers/iterative_improvement.hpp"

namespace marzone {

// Version specific constants
string sVersionString = "Marxan with Zones v 4.0 alpha";
string sIanBallEmail = "ian.ball@aad.gov.au";
string sHughPossinghamEmail = "hugh.possingham@tnc.org";
string sMattWattsEmail = "matt.watts@une.edu.au";
string sMarxanWebSite = "http://marxan.net";
string sDebugTraceFileName = "DebugTraceFile_Marxan_with_Zones.txt";

jmp_buf jmpbuf;
int iMemoryUsed=0;
int iCostCount = 0;
int iZoneCount = 0;
int fSpecPROPLoaded = 0;
long int RandSeed1;
// version 2.0 introduces these features;
//   enhanced flexibility in objectives
//   probabilistic treatment of threats (1D prob)
//   probabilistic treatment of species (2D prob)
//   asymmetric connectivity
//   multiple connectivity files
int iVerbosity;
int savelog;
char* savelogname;
FILE* fsavelog;
int iOptimisationIterativeImprovement = 1, iPuZoneCount = 0;
int iOptimisationCalcPenalties = 1;
double rClocksPerSec;
int asymmetricconnectivity = 0;
int iZoneContrib3On = 0;
int iZonationCost = 0;
mt19937 rngEngine;

int MarZone(string sInputFileName, int marxanIsSecondary)
{
    srunoptions runoptions = {};
    sfname fnames = {};
    sanneal anneal = {};

    int iMessageCounter = 0, iZoneSumSolnIndex, puno,spno,gspno;
    int seedinit, aggexist=0,sepexist=0;
    int *R, *R_CalcPenalties; //,*bestrun,
    int itemp, ipu, iBestRun = 1;
    int i;
    string tempname;
    double rBestScore;
    #ifdef DEBUGTRACEFILE
    string debugbuffer;
    #endif

    rClocksPerSec = CLOCKS_PER_SEC;

    // Handle Error driven termination
    if (setjmp(jmpbuf))
        return 1;

    ShowStartupScreen();

    SetOptions(sInputFileName, runoptions, anneal, fnames);
    SetVerbosity(runoptions.verbose);
    SetRunOptions(runoptions);

    #ifdef DEBUGTRACEFILE
    StartDebugTraceFile();
    debugbuffer = sVersionString + " begin execution\n\nLoadOptions\n";
    AppendDebugTraceFile(debugbuffer);
    if (iVerbosity > 3)
       DumpFileNames(fnames);
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

    if (fnames.savelog)
    {
        SetLogFile(fnames.savelog, fnames.savename + "_log.dat");
    }

    delta = numeric_limits<double>::epsilon(); 
    rngEngine = mt19937(runoptions.iseed); // Initialize rng engine

    #ifdef DEBUG
    SaveSeed(runoptions.iseed);
    #endif

    #ifdef DEBUGTRACEFILE
    AppendDebugTraceFile("RandSeed iseed " + to_string(runoptions.iseed) + "\n");
    #endif

    //  ****     Data File Entry    * * * * * * *
    ShowGenProg("\nEntering in the data files \n");
    AppendDebugTraceFile("before Loading Cost Names\n");
    Costs costs(fnames);

    AppendDebugTraceFile("before Loading Zone Files\n");
    Zones zones(fnames, costs); // load all zone files

    // read in the MarZone files
    AppendDebugTraceFile("before Loading Pu Files (Pu lock, pu zone etc.)\n");
    Pu pu(fnames, costs, asymmetricconnectivity);
    ShowDetProg("    Reading in the Planning Unit names \n");
    ShowGenProg("   There are %i Planning units.\n  %i Planning Unit names read in \n",pu.puno, pu.puno);
    ShowDetProg("    Reading in the species file \n");

    Species spec(fnames);
    AppendDebugTraceFile("After Loading species files\n");
    ShowGenProg("  %i species read in \n",spec.spno);

    /* TODO - update
    #ifdef DEBUGTRACEFILE
    if (iVerbosity > 3)
    {
        DumpZoneNames(iZoneCount,Zones,fnames);
        DumpCostNames(iCostCount,CostNames,fnames);
        DumpZoneContrib(iZoneContribCount,ZoneContrib,fnames);
        DumpZoneContrib2(iZoneContrib2Count,ZoneContrib2,fnames);
        DumpZoneContrib3(iZoneContrib3Count,ZoneContrib3,fnames);
        DumpZoneTarget(iZoneTargetCount,ZoneTarget,fnames);
        DumpZoneTarget2(iZoneTarget2Count,ZoneTarget2,fnames);
        DumpZoneCost(iZoneCostCount,ZoneCost,fnames);
        if (strcmp("NULL",fnames.pulockname) != 0)
            DumpPuLock(iPuLockCount,PuLock,fnames);
        if (strcmp("NULL",fnames.puzonename) != 0)
            DumpPuZone(iPuZoneCount,PuZone,fnames);
        if (strcmp("NULL",fnames.relconnectioncostname) != 0)
            DumpRelConnectionCost(iRelConnectionCostCount,RelConnectionCost,fnames);
    }
    #endif
    */

    // Build zone contributions. 
    zones.BuildZoneContributions(spec, pu);

    if (fnames.zonetargetname.empty() && fnames.zonetarget2name.empty())
    {
        ShowGenProg("Warning: No targets specified for zones.\n");
        AppendDebugTraceFile("Warning: No targets specified for zones.\n");
    }

    AppendDebugTraceFile("before initializing reserve\n");
    Reserve reserve(spec, zones, runoptions.clumptype);
    AppendDebugTraceFile("after initializing reserve\n");


    // parse PuLock and PuZone
    /* TODO - update after refactoring.
    #ifdef DEBUGTRACEFILE
    if (iVerbosity > 3)
    {
        Dump_ZoneContrib(puno,spno,spec,iZoneCount,_ZoneContrib,fnames);
        Dump_ZoneCost(iCostCount,iZoneCount,_ZoneCost,fnames);
        Dump_RelConnectionCost(iZoneCount,_RelConnectionCost,fnames);
        DumpZoneSpec(iMessageCounter,spno,iZoneCount,ZoneSpec,spec,fnames);
        DumpPuLockZone(puno,pu);
    }
    #endif
    */

   // Init reserve object
    Reserve bestyet(spec, zones, runoptions.clumptype);
    reserve.InitializeSolution(pu.puno);

    // Init analysis object
    Analysis analysis;

    if (fnames.savesumsoln)
    {
        AppendDebugTraceFile("before InitSumSoln\n");
        analysis.initSumSolution(puno, iZoneCount);
        AppendDebugTraceFile("after InitSumSoln\n");
    }

    /* TODO - re-enable
    if (asymmetricconnectivity)
    {
        AppendDebugTraceFile("Asymmetric connectivity is on.\n");
        DumpAsymmetricConnectionFile(puno,connections,pu,fnames);

        ShowGenProg("  Asymmetric connectivity is on.\n");
    }
    */ 

    ShowGenProg("  %i connections entered \n", pu.connections.size());
    ShowDetProg("    Reading in the Planning Unit versus Species File \n");
    AppendDebugTraceFile("before LoadSparseMatrix\n");
    pu.LoadSparseMatrix(spec, fnames.inputdir + fnames.puvsprname);
    ShowGenProg("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                pu.puvspr.size(),pu.puno*spec.spno,pu.density);
    AppendDebugTraceFile("after LoadSparseMatrix\n");

    /* Removing the need for this function.
    if (!fnames.matrixspordername.empty())
    {
        AppendDebugTraceFile("before LoadSparseMatrix_sporder\n");
        pu.LoadSparseMatrix_sporder(spec, fnames.inputdir + fnames.matrixspordername);
        AppendDebugTraceFile("after LoadSparseMatrix_sporder\n");
        ShowGenProg("after LoadSparseMatrix_sporder\n");
    }
    */

    #ifdef DEBUGTRACEFILE
    AppendDebugTraceFile("before CalcTotalAreas\n");
    CalcTotalAreas(pu,spec);
    AppendDebugTraceFile("after CalcTotalAreas\n");
    #endif

    /*
    if (fnames.savetotalareas)
    {
        tempname = fnames.savename + "_totalareas" + getFileSuffix(fnames.savetotalareas);
        OutputTotalAreas(puno,spno,pu,spec,SM,tempname,fnames.savepenalty);
    }
    */ 

    // finalise zone and non-zone targets now that matrix has been loaded
    if (spec.fSpecPROPLoaded)
    {
        AppendDebugTraceFile("before ApplySpecProp\n");

        // species have prop value specified
        ApplySpecProp(spec, pu);

        AppendDebugTraceFile("after ApplySpecProp\n");
    }

    AppendDebugTraceFile("before Build_ZoneTarget\n");
    zones.BuildZoneTarget(spec, pu, fnames);
    AppendDebugTraceFile("after Build_ZoneTarget\n");

    /* TODO enable
    if (iVerbosity > 3)
        Dump_ZoneTarget(spno,iZoneCount,_ZoneTarget,fnames);
    #endif
    */

    // Read and process species block definitions
    AppendDebugTraceFile("before process block definitions\n");

    if (!fnames.blockdefname.empty())
    {
        ShowDetProg("    Setting Block Definitions \n");
        spec.SetSpeciesBlockDefinitions(pu.TotalSpeciesAmount(spec.spno));
    }

    spec.SetSpeciesDefaults();

    AppendDebugTraceFile("after process block definitions\n");
    ShowGenProgInfo("Checking to see if there are aggregating or separating species.\n");
    aggexist = spec.aggexist;
    sepexist = spec.sepexist;

    if (fnames.savesen)
    {
        /* TODO - enable
        #ifdef DEBUGTRACEFILE
        AppendDebugTraceFile("before OutputScenario\n");
        #endif

        sprintf(tempname2,"%s_sen.dat",savename);
        OutputScenario(puno,spno,prop,anneal,seedinit,repeats,clumptype,
                      runopts,heurotype,costthresh,tpf1,tpf2,tempname2);

        #ifdef DEBUGTRACEFILE
        AppendDebugTraceFile("after OutputScenario\n");
        #endif
        */
    }

    if (runoptions.verbose > 1)
        ShowTimePassed();


    AppendDebugTraceFile("before InitialiseReserve\n");
    reserve.RandomiseSolution(pu, rngEngine);
    AppendDebugTraceFile("after InitialiseReserve\n");

    // * * *  Pre-processing    * * * * ***
    ShowGenProg("\nPre-processing Section. \n");
    ShowGenProgInfo("    Calculating all the penalties \n");

    if (fnames.penaltyname.empty())
    {
        if (fnames.matrixspordername.empty())
        {
            AppendDebugTraceFile("before CalcPenalties\n");

            // we don't have sporder matrix available, so use slow CalcPenalties method
            itemp = CalcPenalties(pu, spec, zones, reserve, runoptions.clumptype);

            AppendDebugTraceFile("after CalcPenalties\n");
        }
        else
        {
            /* Removing optimise since we no longer need it.
            // we have sporder matrix available, so use optimised CalcPenalties method
            if (iOptimisationCalcPenalties == 1)
            {
                AppendDebugTraceFile("before CalcPenaltiesOptimise\n");

                itemp = CalcPenaltiesOptimise(pu, spec, zones, reserve);

                AppendDebugTraceFile("after CalcPenaltiesOptimise\n");
            }
            else
            {*/
                AppendDebugTraceFile("before CalcPenalties\n");

                itemp = CalcPenalties(pu, spec, zones, reserve, runoptions.clumptype);

                AppendDebugTraceFile("after CalcPenalties\n");
            //}
        }

        if (itemp > 0)
            ShowProg("%d species cannot meet target%c.\n", itemp, itemp == 1 ? ' ' : 's');
    }
    else
    {
        AppendDebugTraceFile("before LoadPenalties\n");
        spec.LoadCustomPenalties(fnames.inputdir + fnames.penaltyname);
        AppendDebugTraceFile("after LoadPenalties\n");
    }

    if (runoptions.AnnealingOn)
    {
       ShowGenProgInfo("    Calculating temperatures.\n");
       if (!anneal.Titns)
          ShowErrorMessage("Initial Temperature is set to zero. Fatal Error \n");

       anneal.Tlen = anneal.iterations/anneal.Titns;
       ShowGenProgInfo("  Temperature length %d \n",anneal.Tlen);
       ShowGenProgInfo("  iterations %i, repeats %i \n",anneal.iterations, runoptions.repeats);
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

       for (int i=1;i<=zones.zoneCount;i++)
       {
           string saveSolutionMatrixNameByZone = fnames.savename + "_solutionsmatrix" + to_string(i) + getFileSuffix(fnames.savesolutionsmatrix);
           // init solutions matrix for each zone separately
           pu.WriteSolutionsMatrixHeader(saveSolutionMatrixNameByZone,fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);
       }
    }

    //   The larger repetition loop
    // TODO - parallelize
    // for each repeat run
    Reserve bestR;
    bestR.objective.total = numeric_limits<double>::max();
    int maxThreads = omp_get_max_threads();
    ShowGenProg("Running multithreaded over number of threads: %d\n", maxThreads);
    for (int irun = 1;irun <= runoptions.repeats;irun++)
    {
        stringstream debugbuffer; 

        // Create new reserve object and init
        Reserve reserveThread(spec, zones, irun);
        reserveThread.InitializeSolution(pu.puno);
        reserveThread.RandomiseSolution(pu, rngEngine);

        debugbuffer << "annealing start run " << irun << "\n";

        ShowGenProg("\n");
        ShowProg("Run %i ",irun);

        SimulatedAnnealing sa(runoptions.AnnealingOn, anneal, rngEngine, fnames.saveannealingtrace);
        if (runoptions.AnnealingOn)
        {
            debugbuffer << "before Annealling Init run "<< irun << "\n";

            // init sa parameters if setting is appropriate
            sa.Initialize(spec, pu, zones, runoptions.clumptype);
            debugbuffer << "after Annealling Init run "<< irun << "\n";

            ShowGenProg("  Using Calculated Tinit = %.4f Tcool = %.8f \n", anneal.Tinit,anneal.Tcool);
           anneal.temp = anneal.Tinit;
        }  // Annealing Settup

        ShowGenProg("  creating the initial reserve \n");

        debugbuffer << "before ZonationCost run " << irun << "%i\n";
        reserveThread.EvaluateObjectiveValue(pu, spec, zones);
        debugbuffer << "after ZonationCost run " << irun << "\n";

        if (runoptions.verbose > 1)
        {
           ShowGenProg("\n  Init:");
           PrintResVal(reserve, spec, zones, runoptions.misslevel);
        }

        if (runoptions.verbose > 5)
        {
           ShowTimePassed();
        }

        // * * * * * * * * * * * * * * * * * * * ***
        // * * *  main annealing algorithm * * * * *
        // * * * * * * * * * * * * * * * * * * * ***

        if (runoptions.AnnealingOn)
        {
            debugbuffer << "before Annealing run " << irun << "\n";
            ShowGenProgInfo("  Main Annealing Section.\n");

            sa.RunAnneal(reserve, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh);

            debugbuffer << "after Annealing run " << irun << "\n";
        } // End of Annealing On

        if (runoptions.HeuristicOn)
        {
           debugbuffer << "before Heuristics run " << irun << "\n";
           Heuristic heur = Heuristic(rngEngine, runoptions.heurotype);
           heur.RunHeuristic(reserve, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh);

           if (runoptions.verbose > 1 && (runoptions.runopts == 2 || runoptions.runopts == 5))
           {
              ShowGenProg("\n  Heuristic:");
              PrintResVal(reserve,spec,zones,runoptions.misslevel);
           }

           debugbuffer << "after Heuristics run " << irun << "\n";
        }    // Activate Greedy

        if (runoptions.ItImpOn)
        {
            debugbuffer << "before IterativeImprovementOptimise run " << irun << "\n";
            IterativeImprovement itImp = IterativeImprovement(rngEngine, fnames, runoptions.itimptype);
            itImp.Run(reserve, spec, pu, zones, runoptions.tpf1, runoptions.tpf2, runoptions.costthresh);
            debugbuffer << "after IterativeImprovementOptimise run " << irun << "\n";

            /* TODO - reenable.
           if (aggexist)
              ClearClumps(spno,spec,pu,SM);
            */

           if (runoptions.verbose > 1)
           {
              ShowGenProg("  Iterative Improvement:");
              PrintResVal(reserve,spec,zones,runoptions.misslevel);
           }
        } // Activate Iterative Improvement

        debugbuffer << "before file output run " << irun << "\n";

        // Start writing output
        string tempname2;
        if (fnames.saverun)
        {
            tempname2 = fnames.savename + "_r" + to_string(irun) + getFileSuffix(fnames.saverun);
            reserve.WriteSolution(tempname2, pu, fnames.saverun);
        }

        debugbuffer << "OutputSolution ran\n";

        if (fnames.savespecies && fnames.saverun)
        {   
            // output species distribution for a run (specific reserve)
            tempname2 = fnames.savename + "_mv" + to_string(irun) + getFileSuffix(fnames.savespecies);
            OutputFeatures(tempname2, zones, reserve, spec,fnames.savespecies, runoptions.misslevel);
        }

        if (fnames.savesum)
        {
            tempname2 = fnames.savename + "_sum" + to_string(irun) + getFileSuffix(fnames.savesum);

            /* TODO
           if (irun == 1)
              OutputSummary(puno,spno,R,spec,reserve,irun,tempname2,misslevel,fnames.savesum);
           else
               OutputSummary(puno,spno,R,spec,reserve,irun,tempname2,misslevel,fnames.savesum);
            */ 
        }

        #ifdef DEBUGFPERROR
        debugbuffer << "OutputSummary ran\n";
        #endif

        // Saving the best from all the runs
        if (fnames.savebest)
        {
            if (reserveThread.objective.total < bestR.objective.total)
            {
                bestR = reserveThread; // deep copy. TODO LOCK

                if (runoptions.verbose > 1)
                {
                    ShowGenProg("  Best:");
                    PrintResVal(bestR, spec, zones, runoptions.misslevel);
                }
            }
        }

        if (fnames.savesumsoln) // Add current run to my summed solution
            analysis.ApplyReserveToSumSoln(reserveThread);

        if (fnames.savesolutionsmatrix)
        {
            string solutionsMatrixName = fnames.savename + "_solutionsmatrix" + getFileSuffix(fnames.savesolutionsmatrix);
            reserveThread.AppendSolutionsMatrix(solutionsMatrixName, zones.zoneCount, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);

           for (i=1;i<=zones.zoneCount;i++)
           {
               string solutionsMatrixZoneName = fnames.savename + "_solutionsmatrix_zone" + to_string(i) + getFileSuffix(fnames.savesolutionsmatrix);
               // append solutions matrix for each zone separately
               reserveThread.AppendSolutionsMatrixZone(solutionsMatrixZoneName, i-1, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);
           }
        }

        if (fnames.savezoneconnectivitysum)
        {
            string zoneConnectivityName = fnames.savename + "_zoneconnectivitysum" + to_string(irun) + getFileSuffix(fnames.savezoneconnectivitysum);
            reserveThread.WriteZoneConnectivitySum(zoneConnectivityName, pu, zones, fnames.savezoneconnectivitysum);
        }

        //if (aggexist)
        //   ClearClumps(spno,spec,pu,SM);

        debugbuffer << "after file output run " << irun << "\n";
        debugbuffer << "annealing end run " << irun << "\n";

        if (marxanIsSecondary == 1)
           WriteSecondarySyncFileRun(irun);

        if (runoptions.verbose > 1)
           ShowTimePassed();

    } // ** the repeats  **

    AppendDebugTraceFile("before final file output\n");

    if (fnames.savebest)
    {
        string saveBestName = fnames.savename + "_best" + getFileSuffix(fnames.savebest);
        reserve.WriteSolution(saveBestName, pu, fnames.savebest);

        AppendDebugTraceFile("Best solution is run " + to_string(iBestRun) + "\n");
        ShowGenProg("\nBest solution is run %i\n",iBestRun);

        if (fnames.savezoneconnectivitysum)
        {
            string saveZoneConnectivityName = fnames.savename + "_zoneconnectivitysumbest" + getFileSuffix(fnames.savezoneconnectivitysum);
            reserve.WriteZoneConnectivitySum(saveZoneConnectivityName, pu, zones, fnames.savezoneconnectivitysum);
        }

        if (fnames.savespecies)
        {
            string saveBestSpeciesName = fnames.savename + "_mvbest" + getFileSuffix(fnames.savespecies);
            OutputFeatures(saveBestSpeciesName,zones, reserve, spec,fnames.savespecies, runoptions.misslevel);
        }
    }

    if (fnames.savesumsoln)
    {
        string saveSumSolnName = fnames.savename + "_ssoln" + getFileSuffix(fnames.savesumsoln);
        analysis.WriteAllSumSoln(saveSumSolnName, pu, zones, fnames.savesumsoln);
    }

    ShowShutdownScreen();

    if (fnames.savelog)
        SetLogFile(0,NULL);  /* tidy up files */

    #ifdef DEBUGTRACEFILE
    AppendDebugTraceFile("end final file output\n");
    AppendDebugTraceFile("\nMarxan with Zones end execution\n");
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

    vector<double>& specTargetZones = zones.AggregateTargetAreaBySpecies();
    vector<int>& specOccurrenceZones = zones.AggregateTargetOccurrenceBySpecies();

    vector<vector<penaltyTerm>>& specPuAmounts = pu.getPuAmountsSorted(spec.spno, true); // list of pus that contribute to a species.

    for (int i = 0; i < spec.spno; i++)
    {
        rZoneSumTarg = specTargetZones[i];
        iZoneSumOcc = specOccurrenceZones[i];

        sspecies &specTerm = spec.specList[i];

        if (specTerm.target > rZoneSumTarg)
            rZoneSumTarg = specTerm.target;
        if (specTerm.targetocc > iZoneSumOcc)
            iZoneSumOcc = specTerm.targetocc;

        if (specTerm.target2 || specTerm.sepnum)
        {
            int j = r.ComputePenaltyType4(specTerm, specPuAmounts[i], i, rZoneSumTarg, specTerm.target2, iZoneSumOcc, clumptype); // TODO
            badspecies += (j > 0);
            goodspecies += (j < 0);

            continue;
        } // Species has aggregation requirements

        itargetocc = 0, ftarget = 0.0, penalty = 0.0;
        int lockedZone = -1;
        double lockedContrib = 0.0;
        // For this species, sum up all locked occurrences and areas. TODO - move this logic somewhere general.
        for (int j : pu.GetPuLockedIndices())
        {
            // Get zoneid of locked pu and zone contrib of this zone for this species
            lockedZone = pu.GetPuLock(j);
            lockedContrib = zones.GetZoneContrib(i, lockedZone);

            if (lockedContrib) {
                // Do this calculation by taking into account zonecontrib of the locked zone. If positive contrib, we count it.
                rAmount = pu.RtnAmountSpecAtPu(j,i)*lockedContrib;
                if (rAmount > 0) {
                    ftarget += rAmount;
                    itargetocc++;
                }

                penalty += rtnMaxNonAvailableCost(j, pu, zones); // TODO - figure out what this is intended to do. Probably should be the cost of the zone if it's counted.
            }

        } // reset PUtemp and also target
        spec.specList[i].penalty = penalty;

        // Already adequately represented on type 2 planning unit
        if (ftarget >= rZoneSumTarg && itargetocc >= iZoneSumOcc)
        {
            goodspecies++;
            ShowGenProgInfo("Species %i (%s) has already met target %.2f\n",
                spec.specList[i].name,spec.specList[i].sname,rZoneSumTarg);

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
            ShowGenProgInfo("Species %d (%s) cannot reach target %.2f there is only %.2f available.\n",
                            spec.specList[i].name, spec.specList[i].sname, rZoneSumTarg, ftarget);

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
         ShowGenProg("%i species are already adequately represented.\n",goodspecies);

    AppendDebugTraceFile("CalcPenalties end\n");

    return(badspecies);
}

double rtnMaxNonAvailableCost(int ipu, Pu& pu, Zones& zones)
{
    double fcost = 0, rMaxCost = 0;

    for (int iZone = 0; iZone < zones.zoneCount; iZone++)
    {
        fcost = ReturnPuZoneCost(ipu, iZone, pu, zones);
        if (fcost > rMaxCost)
            rMaxCost = fcost;
    }

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
// TODO - move to Reserve object
void PrintResVal(Reserve& reserve, Species& spec, Zones& zones, double misslevel)
{
    int i, iMissing;
    double shortfall, rMPM;

    iMissing = reserve.CountMissing(spec, zones, misslevel, rMPM);
    string sPuZones = reserve.CountPuZones(zones);

    ShowProg("Value %.1f Cost %.1f %s Connection %.1f Missing %i Shortfall %.2f Penalty %.1f MPM %.1f\n",
             reserve.objective.total, reserve.objective.cost, sPuZones.c_str(), reserve.objective.connection, iMissing, shortfall, reserve.objective.penalty, rMPM);

} /* * * * Print Reserve Value * * * * */

/* * * * * * * * * * * * * * * * * * * * ****/
/* * * *   Main functions * * * * * * * * ***/
/* * * * * * * * * * * * * * * * * * * * ****/

/* The following function shows the startup screen information.
   The program title and authors */

void ShowStartupScreen(void)
{
    printf("        %s \n\n   Marine Reserve Design with Zoning and Annealing\n\n", sVersionString.c_str());
    printf("   Marxan with Zones coded by Matthew Watts\n");
    printf("   Written by Ian Ball, Hugh Possingham and Matthew Watts\n\n");
    printf("   Based on Marxan coded by Ian Ball, modified by Matthew Watts\n");
    printf("   Written by Ian Ball and Hugh Possingham\n\n");
    printf("%s\n%s\n%s\n\n", sIanBallEmail.c_str(), sHughPossinghamEmail.c_str(), sMattWattsEmail.c_str());
    printf("   Marxan website\n\n");
    printf("%s\n\n", sMarxanWebSite.c_str());

}  /* Show Startup Screen */


/* Show ShutDown Screen displays all the end of program information. It only displays when
    the iVerbosity has been set to 1 or higher */
void ShowShutdownScreen(void)
{
     if (iVerbosity > 0)
     {
        printf("\n");
        ShowTimePassed();
        printf("\n              The End \n");
        if (savelog)
           fprintf(fsavelog,"\n              The End \n");
     }
}

/* ShowPauseExit delivers a message prior to exiting */
void ShowPauseExit(void)
{
     printf("Press return to exit.\n");
     getchar();
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
     WriteSecondarySyncFileRun();
}

void ShowPauseProg(void)
{
     printf("Press return to continue.\n");
     getchar();
} /** Pause **/

/* Set Verbosity sets the module variable iVerbosity to decide how to handle different
    user specified verbosity options */

void SetVerbosity(int verb)
{
     iVerbosity = verb;

} /* Set Verbosity */

/* * * *  ShowTimePassed displays the time passed so far * * * * */

void ShowTimePassed(void)
{
     int itemp, iClock;
     double rTemp;
     #ifdef DEBUGSHOWTIMEPASSED
     char debugbuffer[1000];
     #endif

     #ifdef DEBUGSHOWTIMEPASSED
     AppendDebugTraceFile("ShowTimePassed start\n");
     #endif

     iClock = clock();

     #ifdef DEBUGSHOWTIMEPASSED
     sprintf(debugbuffer,"ShowTimePassed iClock %i\n",iClock);
     AppendDebugTraceFile(debugbuffer);
     sprintf(debugbuffer,"ShowTimePassed rClocksPerSec %g\n",rClocksPerSec);
     AppendDebugTraceFile(debugbuffer);
     #endif

     rTemp = iClock/rClocksPerSec;

     #ifdef DEBUGSHOWTIMEPASSED
     sprintf(debugbuffer,"ShowTimePassed rTemp %g\n",rTemp);
     AppendDebugTraceFile(debugbuffer);
     #endif

     itemp = floor(rTemp);

     #ifdef DEBUGSHOWTIMEPASSED
     sprintf(debugbuffer,"ShowTimePassed itemp %i\n",itemp);
     AppendDebugTraceFile(debugbuffer);
     #endif

     printf("Time passed so far is ");
     if (itemp >= 60*60)
        printf(" %i hour%c,%i min%c and %i secs \n",
               itemp/3600,((itemp/3600==1)?' ':'s'),
               (itemp/60)%60,((itemp/60==1)?' ':'s'),itemp%60);
     else
     {
         if (itemp >=60 )
            printf(" %i min%c and %i secs \n",itemp/60,((itemp/60==1)?' ':'s'),itemp%60);
         else
             printf("%i secs \n",itemp);
     }

     if (savelog)
     {
        fprintf(fsavelog,"Time passed so far is ");
        if (itemp >= 60*60)
           fprintf(fsavelog," %i hour%c,%i min%c and %i secs \n",
                   itemp/3600,((itemp/3600==1)?' ':'s'),
                   (itemp/60)%60,((itemp/60==1)?' ':'s'),itemp%60);
        else
        {
            if (itemp >=60 )
               fprintf(fsavelog," %i min%c and %i secs \n",itemp/60,((itemp/60==1)?' ':'s'),itemp%60);
            else
                fprintf(fsavelog,"%i secs \n",itemp);
        }
     }

     #ifdef DEBUGSHOWTIMEPASSED
     AppendDebugTraceFile("ShowTimePassed end\n");
     #endif
} /* Show Time Passed */

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

    SetVerbosity(1); /* This enables local warning messages */
    SetDefaultOptions(runoptions, anneal, fnames);

    /* Open file and then feed in each variable type */
    vector<string> lines = GetFileLines(sInputFileName, errorMessage);
    if (!errorMessage.str().empty())
        ShowErrorMessage(errorMessage.str());

    readInputOption(lines,"VERSION",version, false, present, warningMessage, errorMessage);
    readInputOption(lines,"PROP",runoptions.prop, false, present, warningMessage, errorMessage);
    readInputOption(lines,"RANDSEED",runoptions.iseed, false, present, warningMessage, errorMessage); /* The random seed. -1 to set by clock */
    if (!present) { //if seed not present, set as time based.
        runoptions.iseed = (long int)time(NULL);
    }

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
        ShowErrorMessage("Multiple Zone Contribution files defined. Please only define one ZONECONTRIBNAME, ZONECONTRIB2NAME or ZONECONTRIB3NAME");
    
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
    if (!fnames.outputdir.empty())
    {
        fnames.savename = fnames.outputdir + fnames.savename;
    }
    
    // Check and print any warning/error messages
    ShowWarningMessage(warningMessage.str());
    if (!errorMessage.str().empty())
        ShowErrorMessage(errorMessage.str());

} /***** Set Options 2* * * */

void CountPuZones2(char *sNames,char *sCounts,int imode,int puno,int R[])
{
     int i,*ZoneCount;
     char sCount[100];
     char sDelimiter[100];

     if (imode > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"   ");

     ZoneCount = (int *) calloc(iZoneCount,sizeof(int));
     for (i=0;i<iZoneCount;i++)
         ZoneCount[i] = 0;

     for (i=0;i<puno;i++)
         if (R[i] > 0)
            ZoneCount[R[i]-1] += + 1;

     strcpy(sNames,"");
     strcpy(sCounts,"");
     for (i=0;i<iZoneCount;i++)
     {
         strcat(sNames,sDelimiter);
         if (imode == 2)
            strcat(sNames,"\"");
         strcat(sNames,Zones[i].name);
         strcat(sNames," PuCount");
         if (imode == 2)
            strcat(sNames,"\"");

         strcat(sCounts,sDelimiter);
         sprintf(sCount,"%i",ZoneCount[i]);
         strcat(sCounts,sCount);
     }

     free(ZoneCount);
}

void CostPuZones(char *sNames,char *sCounts,int imode,int puno,int R[])
{
     int i;
     double *ZoneCosts;
     char sCount[1000];
     char sDelimiter[20];
     double rZoneCost;
     #ifdef DEBUGTRACEFILE
     char debugbuffer[1000];
     #endif

     if (imode > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"   ");

     ZoneCosts = (double *) calloc(iZoneCount,sizeof(double));
     for (i=0;i<iZoneCount;i++)
         ZoneCosts[i] = 0;

     for (i=0;i<puno;i++)
     {
         rZoneCost = ReturnPuZoneCost(i,R[i]);

         #ifdef DEBUGTRACEFILE
         //sprintf(debugbuffer,"CostPuZones ipu %i R %i zonecost %f\n",i,R[i],rZoneCost);
         //AppendDebugTraceFile(debugbuffer);
         #endif

         ZoneCosts[R[i]-1] += rZoneCost;
     }

     strcpy(sNames,"");
     strcpy(sCounts,"");
     for (i=0;i<iZoneCount;i++)
     {
         #ifdef DEBUGTRACEFILE
         sprintf(debugbuffer,"CostPuZones zone %i zonecost %f\n",i,ZoneCosts[i]);
         AppendDebugTraceFile(debugbuffer);
         #endif

         strcat(sNames,sDelimiter);
         if (imode == 2)
            strcat(sNames,"\"");
         strcat(sNames,Zones[i].name);
         strcat(sNames," Cost");
         if (imode == 2)
            strcat(sNames,"\"");

         strcat(sCounts,sDelimiter);
         sprintf(sCount,"%f",ZoneCosts[i]);
         strcat(sCounts,sCount);
     }

     free(ZoneCosts);
}

void InitRandSeed(int iSeed)
{
    if (iSeed>0)
        RandSeed1 = iSeed;
    else
        RandSeed1 = (long int)time(NULL);
    if (RandSeed1 > 0)
        RandSeed1 = -RandSeed1;
}

// penalty associated with separation
/*
double SepPenalty(int ival)
{
       // here ival = 1, 2 or 3. being number of separate locations for speceis

       switch (ival)
       {
              case 1: return(0.5);
              case 2: return(0.2);
              case 3: return (0.0);
       }

       return(0); // This line should never be reached

} // SepPenalty
*/

// * * * **** Sep Penalty 2 * * * * * * * *
// This returns the penalty for not meeting separation requirments. Feed in sepnum and current
//    separation and returns a value from 0 to 1 which is an artificial shortfall.
/*
double SepPenalty2(int ival,int itarget)
{
    double fval;

    if (!itarget)
        return (0); // no penalty if no separation requirement
    fval = (double) ival / (double) itarget;
    if (!ival)
        fval = 1.0 /(double) itarget;

    return(1/(7*fval+0.2)-(1/7.2)); // Gives a nice hyperbole with fval = 1 return 0 and
                                    // fval = 0 or 0.1 returning almost 1
} // SepPenalty2
*/

/* TODO - for Sepnum/Sepdistance
int CheckDistance(int i, int j,struct spustuff pu[],double squaretarget)
{
    // compare x1*x2+y1*y2 with squaretarget
    if ((pu[i].xloc-pu[j].xloc)*(pu[i].xloc-pu[j].xloc) + (pu[i].yloc-pu[j].yloc)* (pu[i].yloc-pu[j].yloc) >= squaretarget)
        return(1);
    else
        return(0);
} // Is Distant returns true if PU's are big enough distance apart
*/

/*
int CountSeparation(int isp,struct sclumps *newno,
                    struct spustuff pu[],struct spu SM[],sspecies spec[],int imode)
{
    // imode 0 = count separation on current
    // imode 1 = count separation if ipu were included
    // imode -1 = count separation if ipu were excluded
    // The following assumes imode = 0 for starters

    struct slink{int id; struct slink *next;} *first = NULL, *second = NULL,*ptemp,*ptest;
    struct sclumps *pclump;
    struct sclumppu *ppu;
    int sepcount = 1,test;
    double targetdist;
    targetdist = spec[isp].sepdistance * spec[isp].sepdistance;

    if (targetdist == 0)
        return(3); // Shortcut if sep not apply to this species
                   // This assumes that 3 is highest sep count

    // Set up the first list
    if (imode == 1)
    {
        if (ValidPU(newno->clumpid,isp,newno,spec,pu,SM,imode))
        {
            ptemp = (struct slink *) malloc(sizeof(struct slink));
            ptemp->id = newno->clumpid;
            ptemp->next = first;
            first = ptemp;
        }
    }
    for (pclump = spec[isp].head;pclump;pclump = pclump->next)
    {
        for (ppu = pclump->head;ppu;ppu= ppu->next)
        {
            if (ValidPU(ppu->puid,isp,newno,spec,pu,SM,imode))
            {
                ptemp = (struct slink *) malloc(sizeof(struct slink));
                ptemp->id = ppu->puid;
                ptemp->next = first;
                first = ptemp;
            }  // Add all valid species bearing PU's to list
        }
    }
    // need to worry about added pu which is not on spec[isp].head list

    // cycle through this list
    while (first)
    {
        test = first->id;
        ptemp = first;
        first = first->next;
        free(ptemp);
        DebugFree(sizeof(struct slink));

        for (ptemp = first;ptemp;ptemp = ptemp->next)
        {
            if (CheckDistance(ptemp->id,test,pu,targetdist))
            {
                for (ptest=second;ptest;ptest = ptest->next)
                {
                    if (CheckDistance(ptemp->id,ptest->id,pu,targetdist))
                    {
                        // Clear all lists
                        while (first)
                        {
                            ptemp = first;
                            first = ptemp->next;
                            free(ptemp);
                            DebugFree(sizeof(struct slink));
                        }
                        while (second)
                        {
                            ptemp = second;
                            second = ptemp->next;
                            free(ptemp);
                            DebugFree(sizeof(struct slink));
                        }
                        return(3);
                    } // I have succeeded in finding what I'm looking for
                }

                 ptest = (struct slink *) malloc(sizeof(struct slink));
                 ptest->id = ptemp->id;
                 ptest->next = second;
                  second = ptest;
                 sepcount = 2; // there is a separation of at least2.
                               // This should be used to cut down calls to this function
            } // I am sufficient distance from test location
        }

        while (second)
        {
            ptemp = second;
            second = ptemp->next;
            free(ptemp);
            DebugFree(sizeof(struct slink));
        } // clear second between tests
    } // finished scanning through list. first is neccessarily empty now

    while (second)
    {
        ptemp = second;
        second = ptemp->next;
        free(ptemp);
        DebugFree(sizeof(struct slink));
    }

    return(sepcount);
} // CountSeparation
*/

// use the prop value from the conservation feature file to set a proportion target for species
void ApplySpecProp(Species& spec, Pu& pu)
{
    vector<double> speciesSums = pu.TotalSpeciesAmount(spec.spno);
    spec.SetSpeciesProportionTarget(speciesSums);
}

void CalcTotalAreas(Pu& pu, Species& spec)
{
    int ipu, i, ism, isp;
    vector<int> TotalOccurrences, TO_2, TO_3;
    vector<double> TA_2, TA_3, TotalAreas = pu.TotalSpeciesAmount(spec.spno);

    if (iVerbosity > 3) 
    {
        TotalOccurrences = pu.TotalOccurrenceAmount(spec.spno);
        TO_2.resize(spec.spno, 0);
        TO_3.resize(spec.spno, 0);
        TA_2.resize(spec.spno, 0.0);
        TA_3.resize(spec.spno, 0.0);

        pu.TotalSpeciesAmountByStatus(TO_2, TO_3, TA_2, TA_3);

        // Write calculated arrays
        spec.WriteTotalAreasAndOccs("MarZoneTotalAreas.csv", TotalOccurrences, TO_2, TO_3, TA_2, TA_3);
    }

    // Set total areas in species.
    spec.SetTotalAreas(TotalAreas);

} // CalcTotalAreas

} // namespace marzone

int main(int argc,char *argv[])
{
    std::string sInputFileName = "input.dat"; // If no arguments then assume the default file name of 'input.dat'
    int marxanIsSecondary = 0;

    if (argc > 1)
        marzone::HandleOptions(argc,argv,sInputFileName,marxanIsSecondary);  // handle the program options

    // Calls the main annealing unit
    if (marzone::MarZone(sInputFileName, marxanIsSecondary))
    {
        if (marxanIsSecondary == 1)
          marzone::SecondaryExit();
        return 1;
    }  // Abnormal Exit
    if (marxanIsSecondary == 1)
        marzone::SecondaryExit();

    return 0;
}

