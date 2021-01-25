#pragma once
// headers for Marxan with Zones

#define DebugFree(x)

#include "common.hpp"
#include "util.hpp"
#include "costs.hpp"
#include "pu.hpp"
#include "reserve.hpp"
#include "species.hpp"
#include "zones.hpp"

#include <random>

namespace marzone
{
    using namespace std;

    double delta;

    /* Protytpe function Headers */
    void SetRunOptions(srunoptions &runoptions);
    int CalcPenalties(Pu &pu, Species &spec, Zones &zones, Reserve &r, int clumptype);

    void PrintResVal(Reserve& reserve, Species& spec, Zones& zones, double misslevel,stringstream& buffer);
    void SetOptions(string &sInputFileName, srunoptions &runoptions, sanneal &anneal, sfname &fnames);

    // Marxan with Zones debugging functions
    /*
    void DumpPuZone_Debug(int iPuZoneCount, struct puzonestruct PuZone[], struct sfname fnames, int iMessage);
    void DumpRelConnectionCost(int iRelConnectionCostCount, struct relconnectioncoststruct RelConnectionCost[], struct sfname fnames);
    void DumpCostValues(int iCostCount, int puno, double **CostValues, struct sfname fnames);
    void DumpCostFieldNumber(int iFieldCount, int CostFieldNumber[], char *sFields, struct sfname fnames);
    void Dump_ZoneContrib(int puno, int spno, sspecies spec[], int iZoneCount, double _ZoneContrib[], struct sfname fnames);
    void Dump_ZoneTarget(int spno, int iZoneCount, struct _zonetargetstruct _ZoneTarget[], struct sfname fnames);
    void Dump_ZoneCost(int iCostCount, int iZoneCount, double _ZoneCost[], struct sfname fnames);
    void DumpPuLockZone(int puno, struct spustuff pu[]);
    void Dump_RelConnectionCost(int iZoneCount, double _RelConnectionCost[], struct sfname fnames);
    void DumpZoneSpec(int iMessage, int spno, int iZoneCount, struct zonespecstruct ZoneSpec[], struct sspecies spec[], struct sfname fnames);
    */

    void BuildZoneSpec(int spno, int iZoneCount, struct zonespecstruct *ZoneSpec[]);
    void InitZoneSpec(int spno, int iZoneCount, struct zonespecstruct ZoneSpec[]);
    void InitSumSoln(int puno, int iZoneCount, int sumsoln[], int ZoneSumSoln[]);
    /* FILEIN.H END */
    /* ************************************************************************** */
    /* FILEOUT.H START */

    void SaveSeed(int iseed);
    void OutputSummary(int puno, int spno, int R[], struct sspecies spec[], struct scost reserve,
                       int itn, char savename[], double misslevel, int imode);
    void OutputScenario(int puno,int spno, int zoneCount, int costCount, Logger& logger,
                        sanneal& anneal, srunoptions& runoptions,
                        string filename);
    void OutputSolution(int puno, int R[], struct spustuff pu[], char savename[], int imode);
    void OutputFeatures(std::string filename, marzone::Zones &zones, marzone::Reserve &reserve, marzone::Species &spec, int imode, double misslevel);
    void OutputSumSoln(int puno, int sumsoln[], int ZoneSumSoln[], int R[], struct spustuff pu[], char savename[], int imode);
    void OutputSpeciesData(int spno, struct sspecies spec[], char savename[]);

    /* FILEOUT.H END */
    /* ************************************************************************** */

    /* ************************************************************************** */
    /* RANDOM.H START */

    void InitRandSeed(int iSeed);

    /* RANDOM.H END */
    /* ************************************************************************** */
    /* SCREENOUT.H START */

    /* Prototype Headers */

    void ShowStartupScreen(void);  /* Displays splash screen info */
    void ShowShutdownScreen(void); /* displays program termination info */
    void ShowPauseExit(void);      /* Press any key to terminate */

    /* Re-enable when separation is reimplemented. 
    double SepPenalty(int ival);
    double SepPenalty2(int ival, int itarget);
    int ValidPU(int ipu, int isp, struct sclumps *newno, struct sspecies spec[], struct spustuff pu[],
                struct spu SM[], int imode);
    int CheckDistance(int i, int j, struct spustuff pu[], double squaretarget);
    int CountSeparation(int isp, struct sclumps *newno,
                        struct spustuff pu[], struct spu SM[], sspecies spec[], int imode);
    struct slink *makelist(int isp, int inpu, int puno, int R[], struct sclumps *newno, struct sspecies spec[],
                           struct spustuff pu[], struct spu SM[], int imode);
    int SepDealList(struct slink *head, typeseplist *Dist, struct spustuff *pu,
                    struct sspecies spec[], int first, int sepnum, double targetdist, int isp);
    int CountSeparation2(int isp, int ipu, struct sclumps *newno, int puno, int R[],
                         struct spustuff pu[], struct spu SM[], sspecies spec[], int imode);
    void CheckDist(struct sseplist *Dist, int sepnum);
    */


    double rtnMaxNonAvailableCost(int ipu, Pu& pu, Zones& zones);

    void StartDebugFile(char sFileName[], char sHeader[], struct sfname fnames);
    void AppendDebugFile(char sFileName[], char sLine[], struct sfname fnames);
    void OutputPenalty(int spno, struct sspecies spec[], char savename[], int iOutputType);
    void LoadPenalty(int spno, struct sspecies spec[], struct sfname fnames);
    void InitSolutionsMatrix(int puno, struct spustuff pu[], char savename[], int iOutputType, int iIncludeHeaders);
    void AppendSolutionsMatrix(int iRun, int puno, int R[], char savename[], int iOutputType, int iIncludeHeaders);
    void AppendSolutionsMatrixZone(int iZone, int iRun, int puno, int R[], char savename[], int iOutputType, int iIncludeHeaders);

    void WriteSecondarySyncFileRun(int iSyncRun);
    void OutputZoneConnectivitySum(int puno, int R[], char savename[], int imode);

    void OutputTotalAreas(int puno, int spno, struct spustuff pu[], struct sspecies spec[], struct spu SM[], char savename[], int iOutputType);
    void WriteStopErrorFile(char sMess[]);
    void CountPuZones2(char *sNames, char *sCounts, int imode, int puno, int R[]);
    void CostPuZones(char *sNames, char *sCounts, int imode, int puno, int R[]);

} // namespace marzone
