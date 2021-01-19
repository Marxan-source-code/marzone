

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

    // TODO move into struct
    int iSparseMatrixFileLength = 0, iSparseMatrixFileLength_sporder = 0;
    int iZoneContribCount = 0, iZoneTargetCount = 0, iZoneCostCount = 0, iPuLockCount = 0;
    int iZoneContrib2Count = 0, iZoneTarget2Count = 0, iRelConnectionCostCount = 0;
    int iZoneContrib3Count = 0;
    int iMessageCounter = 0, iZoneSumSolnIndex, puno,spno,gspno;
    int seedinit, aggexist=0,sepexist=0;
    int *R, *R_CalcPenalties; //,*bestrun,
    int *sumsoln, *ZoneSumSoln, itemp, ipu, iBestRun = 1;
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
    debugbuffer = "RandSeed iseed " + to_string(runoptions.iseed) + "\n";
    AppendDebugTraceFile(debugbuffer);
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
    Reserve reserve(spec, zones);
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
    Reserve bestyet(spec, zones);
    reserve.InitializeSolution(pu.puno);

    // Init analysis object
    Analysis analysis;

    if (fnames.savesumsoln)
    { // TODO - move to analysis.hpp
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
            itemp = CalcPenalties(pu, spec, zones, reserve);

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

                itemp = CalcPenalties(pu, spec, zones, reserve);

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
    int maxThreads = omp_get_max_threads();
    ShowGenProg("Running multithreaded over number of threads: %d\n", maxThreads);
    for (int irun = 1;irun <= runoptions.repeats;irun++)
    {
        stringstream debugbuffer; 

        // Create new reserve object and init
        Reserve reserveThread(spec, zones);
        reserveThread.InitializeSolution(pu.puno);
        reserveThread.RandomiseSolution(pu, rngEngine);

        debugbuffer << "annealing start run " << irun << "\n";
        //AppendDebugTraceFile(debugbuffer);

        ShowGenProg("\n");
        ShowProg("Run %i ",irun);

        SimulatedAnnealing sa(runoptions.AnnealingOn, anneal, rngEngine, fnames.saveannealingtrace);
        if (runoptions.AnnealingOn)
        {
            debugbuffer << "before Annealling Init run "<< irun << "\n";

            // init sa parameters if setting is appropriate
            sa.Initialize(spec, pu, zones);

            debugbuffer << "after Annealling Init run "<< irun << "\n";

            ShowGenProg("  Using Calculated Tinit = %.4f Tcool = %.8f \n", anneal.Tinit,anneal.Tcool);
           anneal.temp = anneal.Tinit;
        }  // Annealing Settup

        ShowGenProg("  creating the initial reserve \n");

        /* TODO reenable
        if (aggexist)
           ClearClumps(spno,spec,pu,SM);
        */

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
           ZonationCost(irun,puno,spno,R,pu,connections,SM,spec,aggexist,&change,clumptype,prop,0);

           if (irun == 1)
           {
              rBestScore = change.total;
              memcpy(bestyet,R,sizeof(int)* puno);
              iBestRun = irun;

              if (verbose >1)
              {
                 ShowGenProg("  Best:");
                 PrintResVal(puno,spno,bestyet,change,spec,misslevel);
              }
           }
           else
           {
               if (change.total <= rBestScore)
               {
                  rBestScore = change.total;
                  memcpy(bestyet,R,sizeof(int)* puno);
                  iBestRun = irun;

                  if (verbose >1)
                  {
                     ShowGenProg("  Best:");
                     PrintResVal(puno,spno,bestyet,change,spec,misslevel);
                  }
               }
           }
        }
        #ifdef DEBUGFPERROR
        sprintf(debugbuffer,"ZonationCost ran\n");
        AppendDebugTraceFile(debugbuffer);
        #endif

        if (fnames.savesumsoln) // Add current run to my summed solution
           for (ipu=0;ipu<puno;ipu++)
               if (R[ipu] > 0)
               {
                  sumsoln[ipu] += 1;
                  iZoneSumSolnIndex = (puno * (R[ipu]-1)) + ipu;
                  ZoneSumSoln[iZoneSumSolnIndex] += 1;
               }

        if (fnames.savesolutionsmatrix)
        {
           if (fnames.savesolutionsmatrix==3)
              sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
           else
           if (fnames.savesolutionsmatrix==2)
              sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
           else
               sprintf(tempname2,"%s_solutionsmatrix.dat",savename);

           AppendSolutionsMatrix(irun,puno,R,tempname2,fnames.savesolutionsmatrix,fnames.solutionsmatrixheaders);

           for (i=1;i<=iZoneCount;i++)
           {
               // append solutions matrix for each zone separately
               if (fnames.savesolutionsmatrix==3)
                  sprintf(tempname2,"%s_solutionsmatrix_zone%i.csv",savename,i);
               else
               if (fnames.savesolutionsmatrix==2)
                  sprintf(tempname2,"%s_solutionsmatrix_zone%i.txt",savename,i);
               else
                   sprintf(tempname2,"%s_solutionsmatrix_zone%i.dat",savename,i);

               AppendSolutionsMatrixZone(i,irun,puno,R,tempname2,fnames.savesolutionsmatrix,fnames.solutionsmatrixheaders);
           }
        }

        if (fnames.savezoneconnectivitysum)
        {
           if (fnames.savezoneconnectivitysum == 3)
              sprintf(tempname2,"%s_zoneconnectivitysum%05i.csv",savename,irun%10000);
           else
           if (fnames.savezoneconnectivitysum == 2)
              sprintf(tempname2,"%s_zoneconnectivitysum%05i.txt",savename,irun%10000);
           else
               sprintf(tempname2,"%s_zoneconnectivitysum%05i.dat",savename,irun%10000);

           OutputZoneConnectivitySum(puno,R,tempname2,fnames.savezoneconnectivitysum);
        }

        if (aggexist)
           ClearClumps(spno,spec,pu,SM);

        #ifdef DEBUGTRACEFILE
        sprintf(debugbuffer,"after file output run %i\n",irun);
        AppendDebugTraceFile(debugbuffer);
        sprintf(debugbuffer,"annealing end run %i\n",irun);
        AppendDebugTraceFile(debugbuffer);
        #endif

        if (marxanisslave == 1)
           WriteSlaveSyncFileRun(irun);

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
    sprintf(debugbuffer,"end final file output\n");
    AppendDebugTraceFile(debugbuffer);
    AppendDebugTraceFile("\nMarxan with Zones end execution\n");
    #endif

    return 0;

} // MarZone

// test the function rtnAmountSpecAtPu
// TODO - move to unit testing
/*
void TestrtnAmountSpecAtPu(int puno, int spno, struct spustuff pu[], struct sspecies spec[], struct spu SM[],
                           struct sfname fnames)
{
     FILE *fp;
     int i,j;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugTestrtnAmountSpecAtPu.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugTestrtnAmountSpecAtPu.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create TestrtnAmountSpecAtPu file %s\n",writename);
     free(writename);
     fputs("puindex,specindex,puname,specname,amount\n",fp);
     for (i=0;i<puno;i++)
         for (j=0;j<spno;j++)
             fprintf(fp,"%d,%d,%d,%d,%g\n",i,j,pu[i].id,spec[j].name,rtnAmountSpecAtPu(pu,SM,i,j));

     fclose(fp);
}
*/

// returns the 0-base index of a species at a planning unit, if the species doesn't occur here, returns -1
int rtnIdxSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;
    //int iTop, iBottom, iCentre, iCount;

    if (PU[iPUIndex].richness > 0)
       //if (PU[iPUIndex].richness < 100)
       {
          for (i=0;i<PU[iPUIndex].richness;i++)
              if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
                 return (PU[iPUIndex].offset + i);
       }

    return -1;
}

// returns the clump number of a species at a planning unit, if the species doesn't occur here, returns 0
int rtnClumpSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;

    if (PU[iPUIndex].richness > 0)
       for (i=0;i<PU[iPUIndex].richness;i++)
           if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
              return SM[PU[iPUIndex].offset + i].clump;

    return 0;
}

// sets the clump number of a species at a planning unit
void setClumpSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex, int iSetClump)
{
     int i;

     if (PU[iPUIndex].richness > 0)
        for (i=0;i<PU[iPUIndex].richness;i++)
            if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
               SM[PU[iPUIndex].offset + i].clump = iSetClump;
}

// returns the amount of a species at a planning unit, if the species doesn't occur here, returns 0
double rtnAmountSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
       int i;

       if (PU[iPUIndex].richness > 0)
          for (i=0;i<PU[iPUIndex].richness;i++)
              if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
                 return SM[PU[iPUIndex].offset + i].amount;

       return 0;
}

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
int CalcPenalties(Pu& pu, Species& spec, Zones& zones, Reserve& r) {
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
            int j = CalcPenaltyType4(i, puno, SM, connections, spec, pu, clumptype, R); // TODO
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

    if (spec.aggexist)
        ClearClumps(spno,spec,pu,SM);

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

/* * * * * * * * * * * * * * * * * * * * *****/
/* * * * *  Clump Utilities * * * * * * * * **/
/* * * * * * * * * * * * * * * * * * * * ****/

// Clear a single Clump
void ClearClump(int isp,struct sclumps *target,struct spustuff pu[],
                struct spu SM[])
{
     struct sclumppu *ppu;

      // Remove all links from this clump
     while (target->head)
     {
           ppu = target->head;
           if (rtnClumpSpecAtPu(pu,SM,ppu->puid,isp) == target->clumpid) // in case pu is in new clump
              setClumpSpecAtPu(pu,SM,ppu->puid,isp,0);
           target->head = ppu->next;
           free(ppu);
           DebugFree(sizeof(struct sclumppu));
     } // Remove all links from this clump

} // Clear Clump (single clump removal)

/**** Function does this cut clump? * * * */
/**** Returns the value of the fragmented clumps if the given PU were removed ***/
/**** If imode == 1 then it will also do a separation count ***/

int ClumpCut(int isp,struct spustuff pu[],
             struct sspecies spec[],struct sclumps *clump,
             struct sclumppu *clumppu,struct sconnections connections[],struct spu SM[],
             double *totalamount,int *totalocc,
             int *iseparation, int imode,int clumptype)
{
    int ineighbour = 0,iclumps = 0;
    struct slink{int id; struct slink *next;} *head = NULL, *newhead,*thead, *clumplist, *clumpcurr;
    struct sneighbour *pnbr;
    struct sclumps *spclump = NULL, *newpclump;
    struct sclumppu *pclumppu;
    double clumpamount, rAmount;
    int iocc;

    *totalamount = 0;
    *totalocc = 0;

    /* Set up spclump for counting Separation */
    if (imode)
    {
        newpclump = (struct sclumps *) malloc(sizeof(struct sclumps));
        newpclump->clumpid = clumppu->puid;
        newpclump->amount = 0;
        newpclump->next = spclump;
        spclump = newpclump;
    }

    /** Generate list of all neighbours and count them **/
    /*First check if there are no neighbours then exit. **/
      /* return null for no clump cut done and need to do separation count */
    if (connections[clumppu->puid].nbrno == 0)
    {
        if (imode)
        {
            *iseparation = CountSeparation(isp,spclump,pu,SM,spec,0);
            free(spclump);
            DebugFree(sizeof(struct sclumps));
        }
        return(0);
    }

    for (pnbr = connections[clumppu->puid].first; pnbr;pnbr = pnbr->next)
    {
        if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == clump->clumpid)
        {
            ineighbour++;
            newhead = (struct slink *) malloc(sizeof(struct slink));
            newhead->id = pnbr->nbr;
            newhead->next = head;
            head = newhead;
        } /** If neighbour is part of the same clump **/
    } /** For cycling through all neighbours **/

    if (ineighbour <= 1) { /* One or fewer neighbours */
        if (imode)
        { /* separation distance called */
            for(pclumppu=clump->head;pclumppu;pclumppu=pclumppu->next)
             if (pclumppu != clumppu)
             {
                newpclump = (struct sclumps *) malloc(sizeof(struct sclumps));
                newpclump->clumpid = pclumppu->puid;
                newpclump->amount = clump->amount - rtnAmountSpecAtPu(pu,SM,clumppu->puid,isp);
                newpclump->next = spclump;
                spclump = newpclump;

             } /* found someone in the clump who is not being removed */

              *iseparation = CountSeparation(isp,spclump,pu,SM,spec,0);
        }
             else (*iseparation = spec[isp].sepnum);
        if (head)
        {
            free(head);
            DebugFree(sizeof(struct slink));
        }
        while (spclump)
        {
            newpclump = spclump;
            spclump = spclump->next;
            free(newpclump);
            DebugFree(sizeof(struct sclumps));
        }  /* clearing up spclump */
        rAmount = rtnAmountSpecAtPu(pu,SM,clumppu->puid,isp);
        *totalamount = PartialPen4(isp,clump->amount-rAmount,spec,clumptype);
        *totalocc = (clump->occs - (rAmount > 0))*(*totalamount > 0); /* count only if still valid size */
    return(0);
    } /** Only one neighbour **/

    /** More than one neighbour. Can they form their own clump? **/
    /* Put first neighbour at head of new list */
    while (head)
    {
    clumpamount = 0;
    iclumps++;
    clumplist = (struct slink *) malloc(sizeof(struct slink));
    clumplist->next = NULL;
    clumplist->id = head->id;
    clumpcurr = clumplist;
    newhead = head;
    head = head->next;
    free(newhead);  /* move first site from head to clumplist */
    DebugFree(sizeof(struct slink));
    ineighbour--;
    do
    {
      for (pnbr = connections[clumpcurr->id].first;pnbr;pnbr = pnbr->next)
      {
          if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == clump->clumpid && pnbr->nbr != clumppu->puid) /* if neighbour in clump but not cut out one */
          {
             for (newhead = clumplist;newhead && newhead->id != pnbr->nbr;newhead= newhead->next)
                 ; /* Cycle through clumplist looking to see if this fellow is already in it */
             if (!newhead)
             {
                newhead = (struct slink *) malloc(sizeof(struct slink));
                newhead->id = pnbr->nbr;
                newhead->next = clumpcurr->next;
                clumpcurr->next = newhead;  /* put this item in my clumplist */
                /* go through neighbour list and see if this one is there */
                for (newhead=head;newhead && newhead->id != pnbr->nbr;newhead = newhead->next)
                    ; /* find this item on the neighbour list */
                if (newhead && newhead->id == pnbr->nbr)
                {
                   ineighbour--;
                   if (newhead == head)
                      head = newhead->next;
                   else
                   {
                       for (thead=head;thead->next != newhead; thead = thead->next)
                           ; /* find link before the one to be removed */
                       thead->next = newhead->next;
                   } /* remove link that is not head */
                   free(newhead);
                   DebugFree(sizeof(struct slink));
                } /* A new neighbour is taken into account*/
             } /* Adding a novel neighbour to list */
          } /* found a neighbour in clump which isn't the one being cut */
      } /* cycling through every neighbour on this clump */

      /* point to next one on list but keep clump head where it is */
      clumpcurr = clumpcurr->next;

    } while (clumpcurr); /* if youv'e run out of new list then...*/

    iocc = 0;
    for(newhead=clumplist;newhead;newhead=newhead->next) {
        rAmount = rtnAmountSpecAtPu(pu,SM,newhead->id,isp);
        clumpamount += rAmount; /* find total amount */
        iocc += (rAmount > 0);
    }
    *totalamount += PartialPen4(isp,clumpamount,spec,clumptype);
    if (PartialPen4(isp,clumpamount,spec,clumptype))
        *totalocc += iocc;

    if (imode)
        for(newhead=clumplist;newhead;newhead=newhead->next) {
            newpclump = (struct sclumps *)malloc(sizeof(struct sclumps));
            newpclump->clumpid = newhead->id;
            newpclump->amount = clumpamount;
            newpclump->next = spclump;
            spclump = newpclump;
    } /* stick this clump into my clump list for separation purposes */

    /* clean up all lists */
    while (clumplist) {
        clumpcurr = clumplist;
        clumplist = clumplist->next;
        free(clumpcurr);
        DebugFree(sizeof(struct slink));
    } /** clean up clumplist **/
    } /*** Continue clump formation whilst there are members in the list*/

    if (imode)
    {
        *iseparation = CountSeparation(isp,spclump,pu,SM,spec,0);
      while (spclump)
      {
        newpclump = spclump;
        spclump = spclump ->next;
        free(newpclump);
        DebugFree(sizeof(struct sclumps));
      } /* clean up separation clump list */
    }
    else
        *iseparation = spec[isp].sepnum;

    while (head) {
        newhead = head;
        head = head->next;
        free(newhead);
        DebugFree(sizeof(struct slink));
    } /** clean up neighbour list **/
    return(iclumps);
} /*** Function Clump Cut.. Do I cut ? ****/


/* * * * * * * * Clear Clumps * * * * * * * */
/*** This is for clean up purposes * * * * */
void ClearClumps(int spno,struct sspecies spec[],struct spustuff pu[],
                 struct spu SM[])
{
     int i;
     struct sclumps *pclump;

     for (i=0;i<spno;i++)
     {
         while (spec[i].head)
         {
               ClearClump(i,spec[i].head,pu,SM);
               pclump = spec[i].head;
               spec[i].head = spec[i].head->next;
               free(pclump);

         }  // Remove each clump
         spec[i].clumps = 0;

     } /** Clear clump for each species **/
} /* * * * Clear\n Clumps ******/


/***** Add New Clump ******/
struct sclumps *AddNewClump(int isp,int ipu,struct sspecies spec[],struct spustuff pu[],struct spu SM[])
{
       int iclumpno = 0;
       struct sclumps *pclump,*pnewclump;
       struct sclumppu *pnewclumppu;
       double rAmount;

       /** find good clump number **/
       pclump = spec[isp].head;
       if (!pclump)
          iclumpno = 1;
       while (!iclumpno)
       {
             if (!pclump->next)
             {
                iclumpno = pclump->clumpid+1;
                break;
             } /* I've found the end of the list */
             if (pclump->next->clumpid - pclump->clumpid > 1)
             {
                iclumpno = pclump->clumpid+1;
                continue;
             } /* Looking for good number */
             pclump = pclump->next;
       } /*  Find first available clump number */

       setClumpSpecAtPu(pu,SM,ipu,isp,iclumpno);
       pnewclump = (struct sclumps *) malloc(sizeof(struct sclumps));
       pnewclump->clumpid = iclumpno;
       if (spec[isp].head)
       {
          pnewclump->next = pclump->next;
          pclump->next = pnewclump;
       } /* Stick clump into correct location */
       else
       {
           spec[isp].head = pnewclump;
           pnewclump->next = NULL;
       } /* First clump on the block */
       /** Add first clumppu to this new clump **/
       pnewclumppu = (struct sclumppu *) malloc(sizeof(struct sclumppu));
       pnewclumppu->puid = ipu;
       pnewclumppu->next = NULL;
       pnewclump->head = pnewclumppu;
       rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
       pnewclump->amount = rAmount;
       pnewclump->occs = (rAmount > 0);

       spec[isp].clumps++;

       return(pnewclump);

}  /*(* * * * *** Add New Clump * * * * ******/

/* * * * * * * * * ADD NEW PU * * * * * * * */
/* * * * **** Add New Planning Unit for a given Species * * * */
/* * * * * * * * * * * * * * * * * * * * ****/
void AddNewPU(int ipu,int isp,struct sconnections connections[],struct sspecies spec[],struct spustuff pu[],
              struct spu SM[],int clumptype)
{
     int ineighbours = 0;
     int iclumpno, iClump;
     struct sneighbour *pnbr;
     struct sclumps *pclump, *pnewclump, *ptempclump;
     struct sclumppu *pnewclumppu;
     double ftemp, rAmount;

     pnbr = connections[ipu].first;
     while (pnbr)
     {
           // Check all the neighbours to see if any are already in clumps
           iClump = rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp);
           if (iClump > 0)
           {  // Neighbour that is part of clump
              ineighbours++;
              if (ineighbours == 1)
              {  // Join to the first clump that is also a neighbour
                 iclumpno = iClump;
                 for (pclump = spec[isp].head; pclump->clumpid != iclumpno;pclump = pclump->next)
                     ;
                 pnewclumppu = (struct sclumppu *) malloc(sizeof(struct sclumppu));
                 pnewclumppu->puid = ipu;
                 pnewclumppu->next = pclump->head;
                 setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,iclumpno);
                 pclump->head = pnewclumppu;

                 /* Remove old value for this clump */
                 ftemp = PartialPen4(isp,pclump->amount,spec,clumptype);
                 spec[isp].amount -= ftemp;
                 spec[isp].occurrence -= pclump->occs *(ftemp > 0);
                 rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
                 pclump->occs += (rAmount > 0);
                 pclump->amount += rAmount;
              }  // Adding the pu to the clump
              else
              {   // pclump points to the good clump
                  if (pclump->clumpid != iClump)
                  {  // Check if this is a different clump
                     // Join this new clump to the old one
                     for (pnewclump= spec[isp].head; pnewclump->clumpid != iClump;pnewclump = pnewclump->next)
                         ;  // point pnewclump to the joining clump
                     // Run through joining clump and tell all pu's their new number
                     for (pnewclumppu = pnewclump->head;pnewclumppu->next;pnewclumppu=pnewclumppu->next)
                         setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,pclump->clumpid);
                     setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,pclump->clumpid);
                     /** cut out this clump and join it to pclump **/
                     pnewclumppu->next = pclump->head;
                     pclump->head = pnewclump->head;
                     pclump->amount += pnewclump->amount;
                     pclump->occs += pnewclump->occs;
                     ftemp = PartialPen4(isp,pnewclump->amount,spec,clumptype);
                     spec[isp].amount -= ftemp;
                     spec[isp].occurrence -= pnewclump->occs * (ftemp > 0);

                     /** Remove clump head and free memory **/
                     if (pnewclump == spec[isp].head)
                        spec[isp].head = pnewclump->next;
                     else
                     {
                         for (ptempclump = spec[isp].head;ptempclump->next != pnewclump;ptempclump = ptempclump->next)
                             ; /** Find clump just before redundant clump **/
                         ptempclump->next = pnewclump->next;
                     }

                     free(pnewclump);
                     DebugFree(sizeof(struct sclumps));

                  } /** Join the two clumps together **/
              } /** Found another neighbour **/
           }
           pnbr = pnbr->next;
     } /** cycling through all the neighbours **/

     /*** Adding a New clump ***/
     if (!ineighbours)
     {
        AddNewClump(isp,ipu,spec,pu,SM);
        ftemp = PartialPen4(isp,rAmount,spec,clumptype);
        spec[isp].amount += ftemp;
        spec[isp].occurrence += (ftemp>0);
     } /** Adding a new clump **/

     /*** Correcting Amount if new clump not added ***/
     if (ineighbours)
     {
        ftemp = PartialPen4(isp,pclump->amount,spec,clumptype);
        spec[isp].amount += ftemp;
        spec[isp].occurrence += pclump->occs * (ftemp > 0);
     }
} /* * * * * * * * Add New Pu * * * * * * * * */

/* * * * ****** REM PU * * * * * * * * * * * * * * * * * * * */
/* * * * *** Remove a planning unit. Note it is similar to CutClump but actually does action */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void RemPu(int ipu, int isp,struct sconnections connections[], struct sspecies spec[],struct spustuff pu[],
           struct spu SM[],int clumptype)
{
    int ineighbours = 0;
    struct slink{int id;struct slink *next;} *head = NULL, *newhead, *thead;
    struct sclumps *oldclump,*pclump;
    struct sclumppu *cppu,*ppu, *clumpcurr, *tppu;
    struct sneighbour *pnbr;
    double oldamount,newamount = 0.0, rAmount;
    int newoccs;

    for (oldclump = spec[isp].head;oldclump && oldclump->clumpid != rtnClumpSpecAtPu(pu,SM,ipu,isp); oldclump= oldclump->next)
    ; /* Find the correct clump to remove */
    if (!oldclump)
       ShowErrorMessage("Serious error in Remove Type 4 species routine\n");

    for (cppu = oldclump->head;cppu->puid != ipu; cppu = cppu->next)
        ; /* Locate the correct clumppu */
    setClumpSpecAtPu(pu,SM,cppu->puid,isp,0);

    oldamount = PartialPen4(isp,oldclump->amount,spec,clumptype);
    spec[isp].amount -= oldamount;
    spec[isp].occurrence -= oldclump->occs * (oldamount > 0);

    for (pnbr = connections[cppu->puid].first;pnbr;pnbr = pnbr->next)
        if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == oldclump->clumpid)
        {
            ineighbours++;
            newhead = (struct slink *)malloc(sizeof(struct slink));
            newhead->id = pnbr->nbr;
            newhead->next = head;
            head = newhead;
        } /* Building the neighbour list */

    if (ineighbours <= 1)
    {
        rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
        oldclump->amount -= rAmount;
        oldclump->occs -= (rAmount > 0);
        newamount = PartialPen4(isp,oldclump->amount,spec,clumptype);
        newoccs = oldclump->occs * (newamount > 0);
        /* remove clumppu */
        if (cppu == oldclump->head)
        {
            oldclump->head = cppu->next;
        }
        else
        {
            for (ppu= oldclump->head;ppu->next != cppu; ppu = ppu->next)
                ; /* find preceding clumppu;*/
            ppu->next = cppu->next;
        }
        free(cppu);
        DebugFree(sizeof(struct sclumppu));
        if (ineighbours < 1)
        {
            if (oldclump == spec[isp].head)
                spec[isp].head = oldclump->next;
            else {
                for (pclump = spec[isp].head;pclump->next != oldclump;pclump = pclump->next)
                ;
                pclump->next = oldclump->next;
            } /* find preceeding clump */
            free(oldclump);
            DebugFree(sizeof(struct sclumps));
            spec[isp].clumps--;
        } /* Removing redundant clump */
        spec[isp].amount += newamount;
        spec[isp].occurrence += newoccs;
        if (head)
        {
            free(head);
            DebugFree(sizeof(struct slink));
        } /* Only need to free head if ineighbours ==1. Then only 1 item in list */
        return;
    } /* I'm not cutting a clump */

    /* Else create new clumps */
    while (head){
        /* take first element as seed for new clump number */
        pclump = AddNewClump(isp,head->id,spec,pu,SM);
        clumpcurr = pclump->head;
        do {
            for (pnbr=connections[clumpcurr->puid].first;pnbr;pnbr=pnbr->next) {
              if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == oldclump->clumpid)
              {
                  if (oldclump->head->puid == pnbr->nbr)
                  {
                     ppu = oldclump->head;
                     oldclump->head = ppu->next;
                  } /* cut out old clump of head */
                  else
                    {
                        for (tppu= oldclump->head;tppu->next->puid != pnbr->nbr;tppu= tppu->next)
                        ; /* find preceeding pu in clump */
                        ppu = tppu->next;
                        tppu->next = ppu->next;
                    } /* cut from middle of old clump */
                     ppu->next = clumpcurr->next;
                     clumpcurr->next = ppu;
                     setClumpSpecAtPu(pu,SM,ppu->puid,isp,pclump->clumpid);
                     rAmount = rtnAmountSpecAtPu(pu,SM,ppu->puid,isp);
                     pclump->amount += rAmount;
                     pclump->occs += (rAmount>0);
                     /* Check if it is on neighbours list and if so then remove it from that list*/
                    if (head->id == ppu->puid)
                    {
                        newhead = head;
                        head = newhead->next;
                        free(newhead);
                        DebugFree(sizeof(struct slink));
                    }
                    else {
                        for (newhead= head;newhead->next && newhead->next->id != ppu->puid;newhead= newhead->next)
                            ; /* check if next one on list is same person */
                        if (newhead->next && newhead->next->id == ppu->puid)
                        {
                            thead = newhead->next;
                            newhead->next = thead->next;
                            free(thead);
                            DebugFree(sizeof(struct slink));
                        } /* cut out none head element */
                    }
                } /* This one is worth removing */
            } /* Cycling through all neighbours */
            clumpcurr = clumpcurr->next;
        } while (clumpcurr); /* Continue until you've added every conceivable neighbour */
        spec[isp].amount += PartialPen4(isp,pclump->amount,spec,clumptype);
        spec[isp].occurrence += pclump->occs * (PartialPen4(isp,pclump->amount,spec,clumptype)>0);
        newhead = head;
        head = newhead->next;
        free(newhead);
        DebugFree(sizeof(struct slink));
    } /** Account for every neighbour in my list **/

    /* Every neighbour in local list has been used and every clump formed*/
    /* Remove old clump */
    /* Worry about change in amount and hence score */
    if (oldclump == spec[isp].head)
    {
        spec[isp].head = oldclump->next;
    }
    else {
        for(pclump=spec[isp].head;pclump->next != oldclump;pclump=pclump->next)
        ; /* find neighbouring clump */
        pclump->next = oldclump->next;
    } /* removing old clump */
    ClearClump(isp,oldclump,pu,SM);
    free(oldclump);
    DebugFree(sizeof(struct sclumps));

} /* Remove a Planning Unit ***/

/* * * * * * * * * * * * * * * * * * * * ****/
/* * * *   Main functions * * * * * * * * ***/
/* * * * * * * * * * * * * * * * * * * * ****/

/* * * * * * * * Setting Clumps for Species Aggregation Rule* * * * **/
void SetSpeciesClumps(int puno,int R[],struct sspecies spec[],struct spustuff pu[],
                      struct spu SM[],struct sconnections connections[],int clumptype)
{
     int i, ipu, isp, ism;

     for (ipu=0;ipu<puno;ipu++)
         if (pu[ipu].richness)
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                if (spec[isp].target2)
                {
                   spec[isp].clumps = 0;
                   //if ((R[ipu]!=iAvailableEquivalentZone) && (R[ipu] != 0) && (SM[ism].amount > 0) && (SM[ism].clump  == 0))
                   if ((SM[ism].amount > 0) && (SM[ism].clump  == 0))
                   {
                      AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype);
                   }// Add a New planning unit
                }// For each type 4 species
            }

} /* * * *  Set Species clumps * * * * *****/

/* * * * **** Species Amounts Type 4 * * * * ******/
/** Assumes Set Species Clumps has been called **/
void SpeciesAmounts4(int isp,struct sspecies spec[],int clumptype)
{
     double ftemp;
     struct sclumps *pclump;

     for (pclump = spec[isp].head;pclump;pclump= pclump->next)
     {
         ftemp = PartialPen4(isp,pclump->amount,spec,clumptype);
         spec[isp].amount += ftemp;
         spec[isp].occurrence += pclump->occs*(ftemp>0);
     }

} /*** Species Amounts 4 **/

/*** Remove Clump Check ***/
/** returns 0 if any member of clump is non-removable, Ie status == 2 **/
int RemClumpCheck(struct sclumps *pclump,struct spustuff pu[])
{
    struct sclumppu *pcpu;

    for (pcpu = pclump->head;pcpu;pcpu = pcpu->next)
        if (pu[pcpu->puid].status == 2)
            return(0);

    return(1);
}

/* * * * * Set Penalties for a given Type 4 Species ***/
/* Returns 1 if the species is a 'bad species' and -1 if it is a 'good species' */
/* Also sticks the penalty into spec[isp].penalty */
int CalcPenaltyType4(int isp,int puno, struct spu SM[],struct sconnections connections[],
                     struct sspecies spec[],struct spustuff pu[],int clumptype,int inputR[])
{
    int i,j,ipu,iputotal = 0,iZone;
    int /*ineighbours = 0,*/iclumpno,badspecies = 0;
    int *R;
    double totalamount,dummy = 0;
    int idummy;
    double cost = 0.0, connection = 0.0, rAmount, rZoneConnectionCost;
    struct slink {int id; struct slink *next;} *plist,*plisthead = NULL,*pdiscard;
    struct sneighbour *pnbr;
    struct sclumps *pclump, *pnewclump;
    struct sclumppu *pnewclumppu, *pcpu;


    R = (int *) calloc(puno,sizeof(int)); /* needed for separation */
    for (i=0;i<puno;i++)
        R[i] = inputR[i];
    /*memcpy(R,pustat,sizeof(struct spustuff)*puno);*/

    /*** Step 1. Make a link list of all the possible PUs to be included ****/
    /*** This might change if I change the species v site into link lists ****/
    plisthead = NULL;
    for (i=0;i<puno;i++)
      if (rtnAmountSpecAtPu(pu,SM,i,isp) > 0)
      {
         if (R[i] == 0)
            continue; // not allowed to consider excluded planning unit
         //if (R[i] != iAvailableEquivalentZone)
         {  /* add to clumps and remove from list */
            AddNewPU(i,isp,connections,spec,pu,SM,clumptype);
            continue;
         } /* checking if PU forced into reserve */
         iputotal++;
         plist = (struct slink *) malloc(sizeof(struct slink));
         plist->id = i;
         plist->next = plisthead;  /* Insert on list */
         plisthead = plist;  /* point head to new number */
      } /** Made link list of all sites with this species **/

    /* Check first to see if I've already satisfied targets for this species */
    SpeciesAmounts4(isp,spec,clumptype);
    if (spec[isp].sepnum>0)
       spec[isp].separation = CountSeparation2(isp,0,0,puno,R,pu,SM,spec,0);
    if ((spec[isp].amount >= spec[isp].target) && (spec[isp].occurrence >= spec[isp].targetocc) && (spec[isp].separation >= spec[isp].sepnum))
    {
       spec[isp].amount = 0;
       spec[isp].occurrence = 0;
       spec[isp].separation = 0;
       /** Clean out all the clump numbers for this species.*/
       while (spec[isp].head)
       {
             ClearClump(isp,spec[isp].head,pu,SM);
             pclump = spec[isp].head;
             spec[isp].head = spec[isp].head->next;
             free(pclump);
             DebugFree(sizeof(struct sclumps));
             spec[isp].clumps = 0;
       }  /** Remove each clump ***/
       free(R); /* dummy array for separation */
       DebugFree(puno * sizeof(int));
       return(-1);
    }  /* return when all targets already met. */

    if (iputotal)
       do
       { /*** take all pu's at random until satisfied or I've run out **/
         /* Pluck a PU out at random */
         ipu = RandNum(iputotal);
         plist = plisthead;
         for (;ipu>0;ipu--)
         {
               plist = plist->next;
         }
         iputotal--;

         /** Add this PU to our system **/
         //do
           iZone = RandNum(iZoneCount) + 1;

         //while (iZone != iAvailableEquivalentZone);

         R[plist->id] = iZone;
         AddNewPU(plist->id,isp,connections,spec,pu,SM,clumptype);

         /** Remove the chosen site from my site list **/
         if (plisthead == plist)
         {
            plisthead = plist->next;
         } /* special case for head of list */
         else
         {
             for (pdiscard = plisthead; pdiscard->next != plist; pdiscard = pdiscard->next)
             {
             }; /*** Find link before plist ***/
             pdiscard->next = plist->next;
         } /* remove plist from the list */
         free(plist);
         DebugFree(sizeof(struct slink));

         /*** Check to see if I should continue by calculating current holdings **/
         SpeciesAmounts4(isp,spec,clumptype);
         if (spec[isp].sepnum>0)
              spec[isp].separation = CountSeparation2(isp,0,0,puno,R,pu,SM,spec,0);

       } while ((spec[isp].amount < spec[isp].target ||
                    spec[isp].separation < spec[isp].sepnum ||
                 spec[isp].occurrence < spec[isp].targetocc)
                && iputotal >= 1  );

    if (spec[isp].amount < spec[isp].target || spec[isp].occurrence < spec[isp].targetocc)
    {
       badspecies = 1;
       ShowGenProg("Species %i cannot be fully represented!\n",spec[isp].name);
    } /*** Record the fact that the species is unrepresentable ***/
    if (spec[isp].separation < spec[isp].sepnum && spec[isp].amount >= spec[isp].target && spec[isp].occurrence >= spec[isp].targetocc)
    {
       badspecies = 1;
       ShowGenProg("Species %i can only get %i separate valid clumps where %i are wanted!\n",
                   spec[isp].name,spec[isp].separation,spec[isp].sepnum);
    } /*** Record the fact that the species is unrepresentable ***/


    /* Search through the clumps looking for any which can be removed */
    /* But only do this if occurrence target met. Otherwise every single pu is neccessary*/
    if (spec[isp].occurrence >= spec[isp].targetocc)
      {
       pclump = spec[isp].head;
       while (pclump)
       {
             i = 0; /* if i becomes and stays 1 then this clump is removable */
             if (RemClumpCheck(pclump,pu))
                i = 1;
             if (i)
             {
                if (spec[isp].occurrence - pclump->occs >= spec[isp].targetocc)
                   i = 1;  /* if pclump-amount < target2 is caught in next step */
                else
                    i = 0;
             } /* Check is occurrence decrease ok? */
             if (i)
             {
                if ((spec[isp].amount - pclump->amount >= spec[isp].target) || (pclump->amount < spec[isp].target2))
                   i = 1;
                else
                    i = 0;
             } /* Check is amount decrease OK? */
             if (i && spec[isp].sepnum)
             {
                j = CountSeparation2(isp,0,pclump,puno,R,pu,SM,spec,-1);
                if ((j < spec[isp].separation) && (j < spec[isp].sepnum))
                   i = 0;
                else
                    i = 1;
                if (!spec[isp].target2)
                   i = 0; /* cannot elegantly remove clumps if species is listed as non-clumping */
             }
             if (i)   /* This is a clump which can be safely removed */
             {  /* cut clump if uneccessary or it is too small */
                if (spec[isp].head == pclump)
                {
                   spec[isp].head = pclump->next;
                }
                else
                {
                    for (pnewclump = spec[isp].head;pnewclump->next != pclump;pnewclump = pnewclump->next)
                        ; /** find clump before pclump **/
                    pnewclump->next = pclump->next;
                }
                while (pclump->head)
                {
                      pnewclumppu = pclump->head;
                      pclump->head = pnewclumppu->next;
                      setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,0);
                      free(pnewclumppu);
                      DebugFree(sizeof(struct sclumppu));
                }
                totalamount -= pclump->amount;
                /* cut out clump and progress pclump*/
                pnewclump = pclump;
                pclump = pclump->next;
                free(pnewclump);
                DebugFree(sizeof(struct sclumps));
                spec[isp].clumps--;
             } /** removing unneccessary pclump **/
             else
                 pclump = pclump->next;

       }
    } /*** Remove unneccesary clumps and links****/


    /** Test all PU's to see if any one of them are superfluous **/
    /* But only do this if occurrence target met. Otherwise every single pu is neccessary*/
    if (spec[isp].occurrence >= spec[isp].targetocc)
    {
       pclump = spec[isp].head;
       while (pclump)
       {
             pcpu = pclump->head;
             while (pcpu)
             {     /** Test to see if this pcpu is necessary **/
                   i = 0;
                   //if (R[pcpu->puid] != 2)
                   if ((R[ipu] > 0) && (pu[ipu].status < 2) && (pu[ipu].fPULock != 1) && (pu[ipu].fPUZone != 1))
                      i = 1;
                   if (i)
                   {
                      rAmount = rtnAmountSpecAtPu(pu,SM,pcpu->puid,isp);
                      if ((pclump->amount - rAmount > spec[isp].target2) && (spec[isp].amount - rAmount > spec[isp].target))
                         i = 1;
                      else
                          i = 0;
                   }  /* doesn't drop amount below min clump size or target */
                   if (i)
                   {
                      if (spec[isp].occurrence > spec[isp].targetocc)
                         i = 1;
                      else
                          i = 0;
                   } /* Does it drop occurrences too low? */
                   if (i)
                   {
                      pnewclump = (struct sclumps *)malloc(sizeof(struct sclumps));
                      pnewclump->clumpid = pcpu->puid;  /* sclump used to store clumpPU info */
                      pnewclump->amount = 0;
                      pnewclump->next = NULL;
                      j = CountSeparation2(isp,pcpu->puid,pnewclump,puno,R,pu,SM,spec,-1);
                      free(pnewclump);
                      if ((j < spec[isp].separation) && (j < spec[isp].sepnum))
                         i = 0;
                      else
                          i = 1;
                   } /* How does count sep fare? */
                   if (i)
                   {
                      if (ClumpCut(isp,pu,spec,pclump,pcpu,connections,SM,&dummy,&idummy,&j,0,clumptype))
                         i = 0;
                      else
                          i = 1;
                   } /* Does it cut the clump? these are not allowed to remove */
                   /* Theoretically they could possible be removed */
                   if (i)  /* Is this removable? */
                   {  /* remove pcpu */
                      setClumpSpecAtPu(pu,SM,pcpu->puid,isp,0);
                      totalamount -= rAmount;
                      pclump->amount -= rAmount;
                      if (pcpu == pclump->head)
                      {
                         pclump->head = pcpu->next;
                         free(pcpu);
                         DebugFree(sizeof(struct sclumppu));
                         pcpu = pclump->head;
                      } /* removing first clump */
                      else
                      {
                          for (pnewclumppu = pclump->head;pnewclumppu->next != pcpu;pnewclumppu = pnewclumppu->next)
                              ; /* find previous pcpu */
                          pnewclumppu->next = pcpu->next;
                          free(pcpu);
                          DebugFree(sizeof(struct sclumppu));
                          pcpu = pnewclumppu->next;
                      } /* removing pcpu when it is not the head */
                   }  /** remove unneccessary clumppu **/
                   else
                       pcpu = pcpu->next; /* moving pointer when it is not removable */
             } /* Checking each pcpu in clump */
             pclump = pclump->next;
       }
    } /** Cycle over each pclump **/

    while (plisthead)
    {
          plist = plisthead;
          plisthead = plisthead->next;
          free(plist);
          DebugFree(sizeof(struct slink));
    } /* Cleaing link list */


    /*** Now count the cost of this particular reserve ****/
    /*** For each clump figure out connection cost ***/
    pclump = spec[isp].head;
    while (pclump)
    {
          iclumpno = pclump->clumpid;
          pcpu = pclump->head;
          while (pcpu)
          {
                if (pu[pcpu->puid].status != 2)
                {
                   cost += pu[pcpu->puid].cost;
                   connection += connections[pcpu->puid].fixedcost;
                } /* only count fixed costs if PU not forced into reserve */
                if (connections[pcpu->puid].nbrno)
                {
                   pnbr = connections[pcpu->puid].first;
                   while (pnbr)
                   {
                         rZoneConnectionCost = _RelConnectionCost[((R[pcpu->puid]-1) * iZoneCount) + (R[pnbr->nbr] - 1)];

                         if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) != iclumpno)
                            connection += (pnbr->cost * rZoneConnectionCost);
                         pnbr = pnbr->next;
                   } /** Counting each individual connection **/
                } /** Counting connection strength if neccessary **/
                pcpu = pcpu->next;
          } /** Checking each PU in clump **/
          pclump = pclump->next;
    } /*** Count cost for each clump ***/

    /* Finally. Calculate penalty from all of this.*/
    spec[isp].penalty = cost + connection;

    /* Consider case where targets cannot be met */
    totalamount = 0;
    if (spec[isp].amount < spec[isp].target)
       totalamount = spec[isp].target / spec[isp].amount;
    if (spec[isp].occurrence < spec[isp].targetocc)
       totalamount += (double) spec[isp].targetocc/(double) spec[isp].occurrence;
    if (totalamount)
       spec[isp].penalty *= totalamount;  /* Scale it up */

    if (spec[isp].sepdistance)
       spec[isp].separation = 1;
    spec[isp].amount = 0; /* because my routines add it in */
    spec[isp].occurrence = 0;
    /** Clean out all the clump numbers for this species.*/
    while (spec[isp].head)
    {
          ClearClump(isp,spec[isp].head,pu,SM);
          pclump = spec[isp].head;
          spec[isp].head = spec[isp].head->next;
          free(pclump);
          DebugFree(sizeof(struct sclumps));
          spec[isp].clumps = 0;
    }  /** Remove each clump ***/

    free(R); /* dummy array for separation */
    DebugFree(puno * sizeof(int));
    return(badspecies);

} /*** Calculate Penalty for a Type 4 Species ***/

/**** Partial Penalty for type 4 species ***/
double PartialPen4(int isp, double amount,struct sspecies spec[],int clumptype)
{
       if (amount >= spec[isp].target2)
          return (amount);    /* I'm not a partial penalty */
       else
           switch (clumptype)
           {
                  case 0: return(0.0); /* default step function */
                  case 1: return(amount/ 2.0); /* nicer step function */
                  case 2: if (spec[isp].target2)
                             return (amount/spec[isp].target2 * amount);
                  default: return(0.0);
           }
}  /* Partial Penalty for type 4 species */

/*** Value for Adding a Planning Unit ****/
double ValueAdd(int isp,int ipu,int puno, int R[],struct sconnections connections[],struct spustuff pu[],
                struct spu SM[],struct sspecies spec[],int clumptype)
{
       int iclumpid,iseparation,oldoccs = 0,occs, iClump, i, iArrayIndex;
       struct sneighbour *pnbr;
       struct slink {int clumpid;double amount;
                     int occs; struct slink *next;} *head = NULL,*plink;
       struct sclumps *pclump,*sepclump=NULL,*psclump;
       struct sclumppu *ppu;
       double amount,oldamount = 0.0,shortfall,zshortfall;

       /* Count neighbours */
       if (connections[ipu].nbrno > 0)
       {
          pnbr = connections[ipu].first;
          while (pnbr)
          {
                iClump = rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp);
                if (iClump)
                {
                   iclumpid = 1;
                   /* Is nbr on my list ?*/
                   for (plink = head;plink;plink=plink->next)
                       if (plink->clumpid == iClump)
                          iclumpid = 0;

                   if (iclumpid)
                   {
                      //ineighbours++;
                      plink = (struct slink *) malloc(sizeof(struct slink));
                      plink->clumpid = iClump;
                      /* find amount for this clump */
                      for (pclump = spec[isp].head;plink->clumpid != pclump->clumpid;pclump = pclump->next)
                          ; /* find the right clump */
                      plink->amount = pclump->amount;
                      plink->occs = pclump->occs;
                      plink->next = head;
                      head = plink;
                      if (spec[isp].sepnum)
                      for (ppu = pclump->head;ppu;ppu=ppu->next)
                      {
                          psclump = (struct sclumps *) malloc(sizeof(struct sclumps));
                          psclump->clumpid = ppu->puid;
                          psclump->next = sepclump;
                          sepclump = psclump;  /* glue to sep list. Still need amount */
                      } /* stick whole clump onto separation clump for later */
                   } /* new neighbour found */
                } /* neighbour of clump */
                pnbr = pnbr->next;
          } /** count all neighbours if they have a clump **/
       }  /* If There are neighbours */

       if (spec[isp].sepnum)
       {
          psclump = (struct sclumps *) malloc(sizeof(struct sclumps));
          psclump->clumpid = ipu;
          psclump->next = sepclump;
          sepclump = psclump;
       } /* Add ipu to my sepclump list */

       /* now I know number and names of neighbouring clumps */
       amount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
       occs = (amount > 0);
       for (plink = head;plink;plink = plink->next)
       {
           amount += plink->amount;
           occs += plink->occs;
           oldamount += PartialPen4(isp,plink->amount,spec,clumptype);
           oldoccs += plink->occs * (PartialPen4(isp,plink->amount,spec,clumptype)>0);
       }

       /* set the sepclump amounts to this new amount */
       if (spec[isp].sepnum)
          for (psclump = sepclump;psclump;psclump = psclump->next)
              psclump->amount = amount;

       amount = PartialPen4(isp,amount,spec,clumptype);
       occs = occs * (amount > 0);

       amount = amount - oldamount; /* amount is change in amount for this species */
       occs = occs - oldoccs;

       if (spec[isp].sepnum)
       {
          iseparation = CountSeparation2(isp,0,sepclump,puno,R,pu,SM,spec,1);  /* imode = 1 doesn't do anything*/
          while (sepclump)
          {
                psclump = sepclump;
                sepclump = sepclump->next;
                free(psclump);
                DebugFree(sizeof(struct sclumps));
          }
       } /* clean up sepcount link list */

       while (head)
       {
             plink = head;
             head = head->next;
             free(plink);
             DebugFree(sizeof(struct slink));
       }  /* Clean up link list */

       /* Return the effective amount for this species */
       /* Old amount + change in amount + separation penalty for changed structure */

       amount = spec[isp].amount + amount;
       shortfall = 0;
       if (spec[isp].target)
          shortfall = amount >= spec[isp].target ? 0 : (spec[isp].target - amount)/spec[isp].target;
       if (spec[isp].targetocc)
       {
          occs = occs + spec[isp].occurrence;
          amount = occs >= spec[isp].targetocc ? 0: ((double)spec[isp].targetocc - (double) occs)/(double)spec[isp].targetocc;
          shortfall += amount;
       }
       /*if (spec[isp].target && spec[isp].targetocc)
          shortfall /= 2;*/
       for (i=0;i<iZoneCount;i++)
       {
           iArrayIndex = (isp * iZoneCount) + i;
           zshortfall = 0;
           if (_ZoneTarget[iArrayIndex].target)
              if (ZoneSpec[iArrayIndex].amount < _ZoneTarget[iArrayIndex].target)
                 zshortfall += (_ZoneTarget[iArrayIndex].target - ZoneSpec[iArrayIndex].amount) / _ZoneTarget[iArrayIndex].target;
           if (_ZoneTarget[iArrayIndex].occurrence)
              if (ZoneSpec[iArrayIndex].occurrence < _ZoneTarget[iArrayIndex].occurrence)
                 zshortfall += (_ZoneTarget[iArrayIndex].occurrence - ZoneSpec[iArrayIndex].occurrence) / _ZoneTarget[iArrayIndex].occurrence;
           /*if (_ZoneTarget[iArrayIndex].target && _ZoneTarget[iArrayIndex].occurrence)
              zshortfall /= 2;*/
           shortfall += zshortfall;
       }
       return(shortfall + SepPenalty2(iseparation,spec[isp].sepnum));
} /*** Value for Adding a Planning Unit ****/


/** Value Remove. The amount of species loss for removing a single pu */
double ValueRem(int ipu,int isp,struct sspecies spec[],struct sconnections connections[],
                struct spustuff pu[],struct spu SM[],int clumptype)
{
       double newamount = 0,amount,shortfall=0,zshortfall;
       struct sclumps *pclump;
       struct sclumppu *ppu;
       int iseparation,newocc = 0,i,iArrayIndex;

       /* locate the clump and clumppu of the target site ipu */
       for (pclump = spec[isp].head; pclump && pclump->clumpid != rtnClumpSpecAtPu(pu,SM,ipu,isp); pclump = pclump->next)
           ; /* locate correct clump list */

       for (ppu = pclump->head;ppu->puid != ipu; ppu = ppu->next)
           ; /* locate the correct pclump pu */

       if (spec[isp].sepnum)
          ClumpCut(isp,pu,spec,pclump,ppu,connections,SM,&newamount,&newocc,&iseparation,1,clumptype);
       else
           ClumpCut(isp,pu,spec,pclump,ppu,connections,SM,&newamount,&newocc,&iseparation,0,clumptype);

       if (spec[isp].target)
       {
          amount = spec[isp].amount + newamount -PartialPen4(isp,pclump->amount,spec,clumptype) ;
          shortfall = amount > spec[isp].target ? 0 : (spec[isp].target - amount)/spec[isp].target;
       }  /* there is an abundance amount */

       if (spec[isp].targetocc)
       {  /* Handle the case where there is a targetocc */
          amount = spec[isp].occurrence +newocc - pclump->occs * (PartialPen4(isp,pclump->amount,spec,clumptype)>0);
          if (amount < spec[isp].targetocc)
             shortfall += ((double) spec[isp].targetocc - amount)/(double) spec[isp].targetocc;
          /*if (spec[isp].target)
             shortfall /= 2;*/
       }
       for (i=0;i<iZoneCount;i++)
       {
           iArrayIndex = (isp * iZoneCount) + i;
           zshortfall = 0;
           if (_ZoneTarget[iArrayIndex].target)
              if (ZoneSpec[iArrayIndex].amount < _ZoneTarget[iArrayIndex].target)
                 zshortfall += (_ZoneTarget[iArrayIndex].target - ZoneSpec[iArrayIndex].amount) / _ZoneTarget[iArrayIndex].target;
           if (_ZoneTarget[iArrayIndex].occurrence)
              if (ZoneSpec[iArrayIndex].occurrence < _ZoneTarget[iArrayIndex].occurrence)
                 zshortfall += (_ZoneTarget[iArrayIndex].occurrence - ZoneSpec[iArrayIndex].occurrence) / _ZoneTarget[iArrayIndex].occurrence;
           /*if (_ZoneTarget[iArrayIndex].target && _ZoneTarget[iArrayIndex].occurrence)
              zshortfall /= 2;*/
           shortfall += zshortfall;
       }

       return(shortfall + SepPenalty2(iseparation,spec[isp].sepnum));
} /** Value for removing a planning unit ****/


/* * * * * * * *   NewPenalty4   * * * * * * * * *****/
/* Calculates the new penalty for adding or removing a PU for species which have
    clumping requirements */


double NewPenalty4(int ipu,int isp,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],
                    int R[],struct sconnections connections[],int imode,int clumptype)
{
       double amount;

       if (imode == 1)
       {
          if (spec[isp].penalty == 0)
             return (0);  // Targets have all already been met
          amount = ValueAdd(isp,ipu,puno,R,connections,pu,SM,spec,clumptype);
       }
       else
       {
           // determine change in this amount
           amount = ValueRem(ipu,isp,spec,connections,pu,SM,clumptype);
       } // removing a planning unit
       return(amount);

}  /*** The new penalty for type 4 species ***/

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

/* SlaveExit does not deliver a message prior to exiting, but creates a file so C-Plan knows marxan has exited */
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

// * * * **** Sep Penalty 2 * * * * * * * *
// This returns the penalty for not meeting separation requirments. Feed in sepnum and current
//    separation and returns a value from 0 to 1 which is an artificial shortfall.
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

int ValidPU(int ipu,int isp,struct sclumps *newno,struct sspecies spec[],struct spustuff pu[],
            struct spu SM[],int imode)
{
    // Returns true if ipu is acceptable as a planning unit
    int i = rtnIdxSpecAtPu(pu,SM,ipu,isp);
    struct sclumps *pclump, *ppu;
    if (newno)
    {
        if (imode == -2)
        {
            if (SM[i].clump == newno->clumpid)
                return(0); // This whole clump is to be removed
        }
        for (ppu=newno;ppu;ppu=ppu->next)
        {
            if (ipu == ppu->clumpid)
            {
                if (ppu->amount < spec[isp].target2)
                    return(0);
                else
                    return(1);
            }// debugging braces
        } // ipu is on list of changed pus
    }  // Above used only when newno is not NULL
    // Find clump
    for (pclump = spec[isp].head;pclump && (SM[i].clump != pclump->clumpid);pclump= pclump->next)
        ; // scan through to find clump
    if (pclump)
    {
        if (pclump->amount <spec[isp].target2)
            return(0);
        else
            return(1);
    } else {
        if (SM[i].amount < spec[isp].target2)
            return(0);
        else
            return(1);
    }
} // ValidPU

int CheckDistance(int i, int j,struct spustuff pu[],double squaretarget)
{
    // compare x1*x2+y1*y2 with squaretarget
    if ((pu[i].xloc-pu[j].xloc)*(pu[i].xloc-pu[j].xloc) + (pu[i].yloc-pu[j].yloc)* (pu[i].yloc-pu[j].yloc) >= squaretarget)
        return(1);
    else
        return(0);
} // Is Distant returns true if PU's are big enough distance apart

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


// * * * ******Make List * * * * * * * * *
// This makes a list of all the valid PUs which occur on the reserve and on which species
//    isp is present (or NULL), in the form of a slink link list

struct slink *makelist(int isp,int ipu,int puno,int R[],
                       struct sclumps *newno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
// imode: 0 : as is. -1 ipu being removed, +1 ipu being added
{
    struct sclumps *pclump;
    struct sclumppu *ppu;
    struct slink *ptemp,*head=NULL;
    int i;

    double rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);

    if (spec[isp].target2)
    { /* if target2 */ /* deal with clumping species differently from non-clumping*/
        if ((imode == 1) && (ValidPU(newno->clumpid,isp,newno,spec,pu,SM,imode)))
        {
            ptemp = (struct slink *) malloc(sizeof(struct slink));
            ptemp->id = newno->clumpid;
            ptemp->next = head;
            head = ptemp;
        }
        for (pclump = spec[isp].head;pclump;pclump = pclump->next)
        {
            for (ppu = pclump->head;ppu;ppu= ppu->next)
            {
                if (ValidPU(ppu->puid,isp,newno,spec,pu,SM,imode))
                {
                    ptemp = (struct slink *) malloc(sizeof(struct slink));
                    ptemp->id = ppu->puid;
                    ptemp->next = head;
                    head = ptemp;
                }  /* Add all valid species bearing PU's to list */
            }
        }
    } else {   /* non clumping species */
        if ((imode ==1) && rAmount)
        {
            ptemp = (struct slink *)malloc(sizeof(struct slink));
            ptemp->id = ipu;
            ptemp->next = head;
            head = ptemp;
        } /* deal with imode == 1 case */

        for (i=0;i<puno;i++)
        {
            if (rAmount && !(imode == -1 && ipu == i))
            {
                ptemp = (struct slink *)malloc(sizeof(struct slink));
                ptemp->id = i;
                ptemp->next = head;
                head = ptemp;
            }
        }
    } /* non clumping species */

    return(head);
} // makelist

/* * * * ***** Sep Deal List * * * * * * * * *****/
/* This funciton is called by count separation2. It takes a link list of sites and 'deals' them
    out on to the seplist */
int SepDealList(struct slink *head, typeseplist *Dist,struct spustuff *pu,
                struct sspecies spec[],int first,int sepnum,double targetdist,
                int isp)
/* Currsep is the current separation maximum it is 0 up to sepnum */
/* first is only needed if maximum is at 0, sepnum is the target separation */
{
    int placefound,currtarget,bestsep=0;
    int currsep;
    struct slink *temp;

    while (head)
    {
        placefound = 0;
        currtarget = first;
        currsep = sepnum;
        do
        {
            if (CheckDistance(head->id,currtarget,pu,targetdist))
            { /* Good distance */
                currsep++;
                if (currsep == spec[isp].sepnum-1)
                {
                    while (head)
                    {
                        temp = head->next;
                        head->next = Dist[currsep].head;
                        Dist[currsep].head = head;
                        head = temp;
                    }  /* glue remaining list on to bottom of Dist. ignoring size and tail as useless */
                    return(currsep);
                } /* Just found valid separation */
                if (Dist[currsep].head)
                {
                    currtarget = Dist[currsep].head->id;
                } else {
                    placefound = 1;
                    Dist[currsep].head = head;
                    Dist[currsep].tail = head;
                    Dist[currsep].size++;
                    head = head->next;
                    Dist[currsep].tail->next = NULL;
                } /* I'm at the end of the line */
            } else {
                placefound = 1;
                Dist[currsep].tail->next = head;
                Dist[currsep].tail = head;
                Dist[currsep].size++;
                head = head->next;
                Dist[currsep].tail->next = NULL;
            } // bad distance
        } while (!placefound); // Doing each individual
        if (currsep > bestsep)
            bestsep = currsep;
    }

    return(bestsep);
} // SepDealList

// This is a modified form of count separation where the user can specify any
// maximum separation distance rather than just assuming a sep distance of three
// ipu and newno used when imode <> 0. When counting as if ipu were added or removed
// ipu used for non-clumping and newno for clumping species

int CountSeparation2(int isp,int ipu,struct sclumps *newno,int puno,int R[],
                     struct spustuff pu[],struct spu SM[],sspecies spec[],int imode)
{
    typeseplist *Dist;
    struct slink *head = NULL,*temp;
    int sepcount,bestsep = 0,i,currcol;
    double targetdist;

    targetdist = spec[isp].sepdistance * spec[isp].sepdistance;

    if (targetdist == 0)
        return(spec[isp].sepnum); /*Shortcut if sep not apply to this species */

    /* Set up array for counting separation */
    Dist = (typeseplist *) calloc(spec[isp].sepnum,sizeof(typeseplist));
    /*First scan through sites. Grab first valid and place rest in lists */
    head = makelist(isp,ipu,puno,R,newno,spec,pu,SM,imode);

    if (!head)
    {
        free(Dist);
        return(0);
    } /* There was nothing to put in the list */


    Dist[0].head = head;
    Dist[0].size = 1;
    Dist[0].tail = head;
    head = head->next;
    Dist[0].tail->next = NULL;
    if (!head)
    {
        free(Dist[0].head);
        free(Dist);
        return(1);
    }  /* There was only one item in the list */

    /* Deal out link list */
    sepcount = SepDealList(head,Dist,pu,spec,Dist[0].head->id,0,targetdist,isp);
    if (sepcount >= spec[isp].sepnum-1)
    {
        /* clean up arrays */
        for (i=0;i<spec[isp].sepnum;i++)
        {
            while(Dist[i].head)
            {
                temp = Dist[i].head;
                Dist[i].head = Dist[i].head->next;
                free(temp);
            }
        }
        free(Dist);
        return(spec[isp].sepnum);
    }  /* I'm at maximum separation */
    bestsep = sepcount;


    do
    {  /* The main Loop */
        for (currcol=0;Dist[currcol+1].head && currcol < spec[isp].sepnum-2;currcol++)
            ;
        if (currcol == 0)
        {
            if (Dist[0].size < spec[isp].sepnum)
            {
                while (Dist[0].head)
                {
                    temp = Dist[0].head;
                    Dist[0].head = Dist[0].head->next;
                    free(temp);
                }
                free(Dist);
                return(bestsep + 1);
            } /* cannot increase separation terminate function */
            else
            {
                temp = Dist[0].head;
                Dist[0].head = Dist[0].head->next;
                head = Dist[0].head->next;
                Dist[0].head->next = NULL;
                Dist[0].size = 1;
                Dist[0].tail = Dist[0].head;
                free(temp);
                sepcount = SepDealList(head,Dist,pu,spec,Dist[0].head->id,0,targetdist,isp);
            }
        } /* Deal with first column */
        else
        {
            if (Dist[currcol].size + currcol  < spec[isp].sepnum)
            {
                Dist[currcol-1].tail->next = Dist[currcol].head;
                Dist[currcol-1].tail = Dist[currcol].tail;
                Dist[currcol-1].tail->next = NULL;
                Dist[currcol-1].size += Dist[currcol].size;
                Dist[currcol].head = NULL;
                Dist[currcol].size = 0;
                Dist[currcol].tail = NULL;
                sepcount = 0;
            } /* list is not long enough to increase sepcount */
            else
            {
                Dist[currcol-1].tail->next = Dist[currcol].head;
                Dist[currcol-1].tail = Dist[currcol].head;
                Dist[currcol-1].size++;
                Dist[currcol].head = Dist[currcol].head->next;
                head = Dist[currcol].head->next;
                Dist[currcol].head->next = NULL;
                Dist[currcol-1].tail->next = NULL;
                Dist[currcol].tail = Dist[currcol].head;
                Dist[currcol].size = 1;
                sepcount = SepDealList(head,Dist,pu,spec,Dist[currcol].head->id,currcol,targetdist,isp);
            } /* else this column might be long enough */
        } /* Deal with columns other than the first */
        if (sepcount > bestsep)
            bestsep = sepcount;
    } while (bestsep < spec[isp].sepnum-1); /* Main loop. */

    for (i=0;i<spec[isp].sepnum;i++)
    {
        while (Dist[i].head)
        {
            temp = Dist[i].head;
            Dist[i].head = Dist[i].head->next;
            free(temp);
        }
    }
    free(Dist);
    return(bestsep+1);
} // CountSeparation2

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

