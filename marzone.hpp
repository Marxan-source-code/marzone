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
    void CalcTotalAreas(Pu& pu, Species& spec, string filename = "MarZoneTotalAreas.csv", bool save = false);
    double ReturnPuZoneCost(int ipu,int iZone, Pu& pu, Zones& zones);

    void ApplySpecProp(Species& spec, Pu& pu);

    void PrintResVal(Reserve& reserve, Species& spec, Zones& zones, double misslevel,stringstream& buffer);
    void SetOptions(string &sInputFileName, srunoptions &runoptions, sanneal &anneal, sfname &fnames);

    /* ************************************************************************** */
    /* FILEOUT.H START */

    void SaveSeed(int iseed);
    void DumpFileNames(sfname& fnames, Logger& logger);
    string OutputSummaryString(Pu& pu, Species& spec, Zones& zones, Reserve& r, double misslevel, int imode, double blm);
    void OutputSummary(Pu& pu, Zones& zones, vector<string>& summaries, string filename, int imode);
    void OutputScenario(int puno,int spno, int zoneCount, int costCount, Logger& logger,
                        sanneal& anneal, srunoptions& runoptions,
                        string filename);
    void OutputFeatures(std::string filename, marzone::Zones &zones, marzone::Reserve &reserve, marzone::Species &spec, int imode, double misslevel);

    /* FILEOUT.H END */
    /* ************************************************************************** */

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

    void StartDebugFile(string sFileName,string sHeader, sfname& fnames);
    void AppendDebugFile(string sFileName,string& sLine, sfname& fnames);

    void WriteSecondarySyncFileRun(int iSyncRun);

} // namespace marzone
