#pragma once
// Common data structures for Marxan with Zones, accessed by most files/classes

#include <string>
#include <vector>

namespace marzone {
    using namespace std;

    // type definitions for Marxan with Zones data structures
    typedef struct stringname
    {
        int id;
        string name;
    } stringname;

    typedef struct zonecontribstruct
    {
        int zoneid;
        int speciesid;
        double fraction;
    } zonecontribstruct;

    typedef struct zonecontrib2struct
    {
        int zoneid;
        double fraction;
    } zonecontrib2struct;

    typedef struct zonecontrib3struct
    {
        int zoneid;
        int puid;
        int speciesid;
        double fraction;
    } zonecontrib3struct;

    typedef struct puzonestruct
    {
        int puid;
        int zoneid;
    } puzonestruct;

    typedef struct zonespecstruct
    {
        double amount;
        int occurrence;
    } typezonespec;

    typedef struct spu
    {
        double amount;
        int clump;
        int spindex;
    } spu;

    typedef struct spusporder
    {
        double amount;
        int puindex;
    } typepusporder;

    typedef struct penaltyTerm {
        double amount; // amount contributed
        double cost; // cost of pu
    } penaltyTerm; // used for penalty calculation sorting. 

    typedef struct lockedPenaltyTerm {
        int lockedZoneId;
        int puindex;
        double amount;
    } lockedPenaltyTerm;

    // type definitions for original Ian Ball Marxan data structures

    typedef struct spustuff
    {
        int id;
        int status;
        double xloc,yloc;
        double cost;
        int richness,offset;
        int fPULock;
        int iPULock;
        int numZones;
        int iPreviousStatus;
        vector<double> costBreakdown;
    } spustuff;

    typedef struct scost
    {
        double total;
        double connection;
        int missing;
        double penalty;
        double cost;
        double threshpen;
        double shortfall;
    } scost;

    /* General Species information structure */
    typedef struct sgenspec
    {
        int name; // id of species
        string sname;
        int type;
        int targetocc;
        double target;
        double target2; /* Only some species need this (clumping species) */
        int sepnum;
        double sepdistance;
        double prop;
        double spf;
    } typegenspec;

    // contains information about species specific to a reserve
    typedef struct reservespecies {
        double amount;
        int occurrence;
        int clumps;
        int separation;
    } reservespecies;

    typedef struct sspecies : sgenspec
    {
        double penalty;
        int richness,offset;
        double totalarea;
    }sspecies;

    /* Connection Structure. Fixed connection number. Should replace with link list! */
    typedef struct sneighbour
    {
        int nbr;
        double cost;
        int connectionorigon;
    } sneighbour;

    typedef struct sconnections
    {
        vector<sneighbour> first;
        double fixedcost;
        int nbrno;
    } sconnections;

    struct sclumps
    {
      int clumpid;
      double amount;
      int occs;
    }; /* Clump nodes for species Clump Structure */


/*** Annealing Control ****/

    typedef struct sanneal
    {
        int Titns;
        int iterations;
        int Tlen;
        double Tinit;    /* Initial Temperature */
        double Tcool;    /* Cooling Factor */
        double temp; /* Current Temperature */
        double tempold;
        int type;    /* Type of annealing. 0 = none, 1 = fixed, 2 = adaptive */
        double sigma; /*Used in adaptive annealing */
        double sum; /* Used in adaptive annealing */
        double sum2; /* used in adaptive annealing */
    } sanneal; /* Annealing Control handler */

    /* Input File Name Structure */
    typedef struct sfname
    {
        string inputdir;
        string outputdir;
        string specname;
        string puname;
        string puvsprname;
        string matrixspordername;
        string connectionname;
        string blockdefname;
        string zonesname;
        string costsname;
        string zonecontribname;
        string zonecontrib2name;
        string zonecontrib3name;
        string zonetargetname;
        string zonetarget2name;
        string zonecostname;
        string pulockname;
        string puzonename;
        string relconnectioncostname;
        string penaltyname;
        string connectivityfilesname;
        string savename; // generic savename prefix
        int savebest;
        int saverun;
        int savesum;
        int savesen;
        int savespecies;
        int savesumsoln;
        int savelog;
        int savesnapsteps;
        int savesnapchanges;
        int savesnapfrequency;
        int savepenalty;
        int savetotalareas;
        int savesolutionsmatrix;
        int solutionsmatrixheaders;
        int saveannealingtrace;
        int annealingtracerows;
        int suppressannealzones;
        int saveitimptrace;
        int itimptracerows;
        int savespeciesdata;
        int savezoneconnectivitysum;
    } sfname;

    // Contains run options about which algorithms to run, as well as other run constants
    typedef struct srunoptions
    {
        // Algorithms to run
        int AnnealingOn;
        int HeuristicOn;
        int ItImpOn;

        // Algorithm settings
        int heurotype;
        int clumptype;
        int itimptype;

        // Run constants
        double prop;
        long int iseed;
        double misslevel;
        double costthresh;
        double tpf1, tpf2;
        int repeats, runopts, verbose;

    } srunoptions;

    typedef struct ZoneName
    {
        string name;
        uint64_t index;
    } ZoneName;

} // namespace marzone