#pragma once

#include <string>

#include "common.hpp"
#include "logger.hpp"
#include "zones.hpp"

namespace marzone {

void DumpCostValues(int iCostCount, int puno,double **CostValues,struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i,j;

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugCostValues.csv") + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debugCostValues.csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create DumpCostValues file %s\n",writename);
    free(writename);

    fprintf(fp,"puindex");
    for (j=0;j<iCostCount;j++)
        fprintf(fp,",%d",j);
    fprintf(fp,"\n");

    for (i=0;i<puno;i++)
    {
        fprintf(fp,"%d,",i);

        for (j=0;j<iCostCount;j++)
        {
            fprintf(fp,"%lf",CostValues[i][j]);
            if (j != (iCostCount-1))
                fprintf(fp,",");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

/* TODO - rewrite
void Dump_ZoneContrib(int puno,int spno,typesp spec[],int iZoneCount,double _ZoneContrib[],struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i,j,k;
    char messagebuffer[1000];

    if (strcmp("NULL",fnames.zonecontrib3name) != 0)
    {
        for (i=0;i<iZoneCount;i++)
        {
            // create one debug file for each zone. Each file is a 2d matrix of puXfeature
            sprintf(messagebuffer,"%i",i+1);
            writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_Zone_Contrib.csv") + 5, sizeof(char));
            strcpy(writename,fnames.outputdir);
            strcat(writename,"debug_Zone");
            strcat(writename,messagebuffer);
            strcat(writename,"_Contrib.csv");
            fp = fopen(writename,"w");
            if (fp==NULL)
                ShowErrorMessage("cannot create Dump_ZoneContrib file %s\n",writename);
            free(writename);

            // header row
            fprintf(fp,"pu/feature");
            for (j=0;j<spno;j++)
                fprintf(fp,",%i",(j+1));
            fprintf(fp,"\n");

            for (k=0;k<puno;k++)
            {
                fprintf(fp,"%i",k+1);

                for (j=0;j<spno;j++)
                {
                    fprintf(fp,",%lf",_ZoneContrib[(j*puno*iZoneCount)+(k*iZoneCount)+i]);

                    //_ZoneContrib index = (zero_base_feature_index * number_of_planning_units * number_of_zones) +
                    //                     (zero_base_planning_unit_index * number_of_zones) +
                    //                     zero_base_zone_index
                }

                fprintf(fp,"\n");
            }

            fclose(fp);
        }
    } else {
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_ZoneContrib.csv") + 2, sizeof(char));
        strcpy(writename,fnames.outputdir);
        strcat(writename,"debug_ZoneContrib.csv");
        fp = fopen(writename,"w");
        if (fp==NULL)
            ShowErrorMessage("cannot create Dump_ZoneContrib file %s\n",writename);
        free(writename);

        // header row
        fprintf(fp,"spname,spindex");
        for (i=0;i<iZoneCount;i++)
            fprintf(fp,",contrib%i",(i+1));
        fprintf(fp,"\n");

        for (j=0;j<spno;j++)
        {
            fprintf(fp,"%i,%i",spec[j].name,j);

            for (i=0;i<iZoneCount;i++)
                fprintf(fp,",%lf",_ZoneContrib[(j*iZoneCount)+i]);

            fprintf(fp,"\n");
        }

        fclose(fp);
    }
}

void Dump_ZoneTarget(int spno,int iZoneCount,struct _zonetargetstruct _ZoneTarget[],struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i,j;
    char debugbuffer[1000];

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_ZoneTarget.csv") + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debug_ZoneTarget.csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create Dump_ZoneTarget file %s\n",writename);
    free(writename);

    // write header row
    fprintf(fp,"spname,spindex,");
    for (i=0;i<iZoneCount;i++)
    {
        fprintf(fp,"zone%dtarget,zone%doccurrence",i+1,i+1);
        if (i != (iZoneCount-1))
            fprintf(fp,",");
    }
    fprintf(fp,"\n");

    for (j=0;j<spno;j++)
    {
        fprintf(fp,"%i,%i,",spec[j].name,j);
        for (i=0;i<iZoneCount;i++)
        {
            fprintf(fp,"%lf,%i",_ZoneTarget[(j*iZoneCount)+i].target,_ZoneTarget[(j*iZoneCount)+i].occurrence);
            if (i != (iZoneCount-1))
                fprintf(fp,",");
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
}

void Dump_ZoneCost(int iCostCount,int iZoneCount,double _ZoneCost[],struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i,j;

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_ZoneCost.csv") + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debug_ZoneCost.csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create Dump_ZoneCost file %s\n",writename);
    free(writename);

    fprintf(fp,"costindex");
    for (i=0;i<iZoneCount;i++)
        fprintf(fp,",%d",i);
    fprintf(fp,"\n");

    for (j=0;j<iCostCount;j++)
    {
        fprintf(fp,"%d",j);

        for (i=0;i<iZoneCount;i++)
        {
            fprintf(fp,",%lf",_ZoneCost[(j*iZoneCount)+i]);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
}

void DumpPuZone_Debug(int iPuZoneCount,struct puzonestruct PuZone[],struct sfname fnames,int iMessage)
{
    FILE *fp;
    char *writename;
    char messagebuffer[1000];
    int i;

    sprintf(messagebuffer,"%i",iMessage);

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugPuZone_.csv") + strlen(messagebuffer) + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debugPuZone_");
    strcat(writename,messagebuffer);
    strcat(writename,".csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create DumpPuZone file %s\n",writename);
    free(writename);

    fprintf(fp,"puid,zoneid\n");

    for (i=0;i<iPuZoneCount;i++)
    {
        fprintf(fp,"%d,%d\n",PuZone[i].puid,PuZone[i].zoneid);
    }

    fclose(fp);
}

void Dump_RelConnectionCost(int iZoneCount,double _RelConnectionCost[],struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i,j;

    #ifdef DEBUGTRACEFILE
    AppendDebugTraceFile("Dump_ZoneConnectionCost start\n");
    #endif

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_ZoneConnectionCost.csv") + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debug_ZoneConnectionCost.csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create Dump_ZoneConnectionCost file %s\n",writename);
    free(writename);

    fprintf(fp,"zoneindex");
    for (j=0;j<iZoneCount;j++)
        fprintf(fp,",%d",j);
    fprintf(fp,"\n");

    for (j=0;j<iZoneCount;j++)
    {
        fprintf(fp,"%d",j);

        for (i=0;i<iZoneCount;i++)
        {
            fprintf(fp,",%lf",_RelConnectionCost[(j*iZoneCount)+i]);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);

    #ifdef DEBUGTRACEFILE
    AppendDebugTraceFile("Dump_ZoneConnectionCost end\n");
    #endif
}

void DumpZoneSpec(int iMessage,int spno,int iZoneCount,struct zonespecstruct ZoneSpec[],struct sspecies spec[],struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i,j;
    char messagebuffer[1000];
    char debugbuffer[1000];

    sprintf(messagebuffer,"%i",iMessage);

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneSpec_.csv") + strlen(messagebuffer) + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debugZoneSpec_");
    strcat(writename,messagebuffer);
    strcat(writename,".csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneSpec file %s\n",writename);
    free(writename);

    // write header row, species identifier is spec.name integer
    fprintf(fp,"spname,spindex,amount,occurrence,");
    for (i=0;i<iZoneCount;i++)
    {
        fprintf(fp,"amount%i,occurrence%i",(i+1),(i+1));
        if (i != (iZoneCount-1))
            fprintf(fp,",");
    }
    fprintf(fp,"\n");

    for (j=0;j<spno;j++)
    {
        fprintf(fp,"%i,%i,%lf,%i,",spec[j].name,j,spec[j].amount,spec[j].occurrence);

        for (i=0;i<iZoneCount;i++)
        {
            fprintf(fp,"%lf,%i",ZoneSpec[(j*iZoneCount)+i].amount,ZoneSpec[(j*iZoneCount)+i].occurrence);
            if (i != (iZoneCount-1))
                fprintf(fp,",");
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
}

*/
#ifdef DEBUGTRACEFILE
void DumpR(int iMessage,char sMessage[],int puno,int reservedarray[],struct spustuff pu[],struct sfname fnames)
{
    FILE *fp;
    char *writename;
    int i;
    char messagebuffer[1000];
    char debugbuffer[1000];

    sprintf(messagebuffer,"%s%i",sMessage,iMessage);

    sprintf(debugbuffer,"DumpR %s %i\n",sMessage,iMessage);
    AppendDebugTraceFile(debugbuffer);

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugR_.csv") + strlen(messagebuffer) + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debugR_");
    strcat(writename,messagebuffer);
    strcat(writename,".csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create DumpR file %s\n",writename);
    free(writename);

    // write header row
    fprintf(fp,"puid,reservedarray\n");

    for (i=0;i<puno;i++)
    {
        fprintf(fp,"%i,%i\n",pu[i].id,reservedarray[i]);
    }

    fclose(fp);

    sprintf(debugbuffer,"DumpR %s %i end\n",sMessage,iMessage);
    AppendDebugTraceFile(debugbuffer);
}
#endif

void DumpFileNames(sfname& fnames, Logger& logger)
{
    FILE *fp;
    string writename;

    writename = fnames.outputdir + "debugFileNames.csv";
    fp = fopen(writename.c_str(),"w");
    if (fp==NULL)
        logger.ShowErrorMessage("cannot create DumpFileNames file " + writename + "\n";

    fprintf(fp,"input name,file name\n");

    fprintf(fp,"zonesname,%s\n",fnames.zonesname);
    fprintf(fp,"costsname,%s\n",fnames.costsname);
    fprintf(fp,"zonecontribname,%s\n",fnames.zonecontribname);
    fprintf(fp,"zonecontrib2name,%s\n",fnames.zonecontrib2name);
    fprintf(fp,"zonecontrib3name,%s\n",fnames.zonecontrib3name);
    fprintf(fp,"zonetargetname,%s\n",fnames.zonetargetname);
    fprintf(fp,"zonetarget2name,%s\n",fnames.zonetarget2name);
    fprintf(fp,"zonecostname,%s\n",fnames.zonecostname);
    fprintf(fp,"pulockname,%s\n",fnames.pulockname);
    fprintf(fp,"puzonename,%s\n",fnames.puzonename);
    fprintf(fp,"zoneconnectioncostname,%s\n",fnames.relconnectioncostname);

    fclose(fp);
}

/*
void DumpPuLockZone(int puno,struct spustuff pu[])
{
    FILE *fp;
    char *writename;
    int i;

    writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugPuLockZone.csv") + 2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debugPuLockZone.csv");
    fp = fopen(writename,"w");
    if (fp==NULL)
        ShowErrorMessage("cannot create DumpPuLockZone file %s\n",writename);
    free(writename);

    fprintf(fp,"id,fPULock,iPULock,fPUZone,iPUZone,iPUZones,iPreviousStatus\n");
    for (i=0;i<puno;i++)
    {
        fprintf(fp,"%d,%d,%d,%d,%d,%d,%d\n",pu[i].id,pu[i].fPULock,pu[i].iPULock,pu[i].fPUZone,pu[i].iPUZone,pu[i].iPUZones,pu[i].iPreviousStatus);
    }

    fclose(fp);
}
*/

void OutputZonationCostDebugTable(int spno,struct sspecies spec[],char savename[])
{
    FILE *fp;
    int i, k, iExistingArrayIndex;
    double rShortfall, rTotalShortfall;

    fp = fopen(savename,"w");

    fprintf(fp,"SPID,target,amount,shortfall");
    for (k=1;k<=iZoneCount;k++)
        fprintf(fp,",Z%itarget,Z%iamount,Z%ishortfall",k,k,k);
    fprintf(fp,",total shortfall\n");

    for (i=spno-1;i>=0;i--)
    {
        rShortfall = 0;
        rTotalShortfall = 0;

        if (spec[i].target > 0)
        {
            if (spec[i].target > spec[i].amount)
            {
                rShortfall = spec[i].target - spec[i].amount;
            }
        }

        fprintf(fp,"%i,%f,%f,%f",spec[i].name,spec[i].target,spec[i].amount,rShortfall);

        rTotalShortfall += rShortfall;
        rShortfall = 0;

        for (k=0;k<iZoneCount;k++)
        {
            iExistingArrayIndex = (i * iZoneCount) + k;

            if (_ZoneTarget[iExistingArrayIndex].target > 0)
            {
                if (_ZoneTarget[iExistingArrayIndex].target > ZoneSpec[iExistingArrayIndex].amount)
                {
                    rShortfall = _ZoneTarget[iExistingArrayIndex].target - ZoneSpec[iExistingArrayIndex].amount;
                }
            }

            fprintf(fp,",%f,%f,%f",_ZoneTarget[iExistingArrayIndex].target,ZoneSpec[iExistingArrayIndex].amount,rShortfall);

            rTotalShortfall += rShortfall;
            rShortfall = 0;
        }

        fprintf(fp,",%f\n",rTotalShortfall);
    }

    fclose(fp);
}

void WriteSecondarySyncFile(void)
{
    FILE* fsync;

    fsync = fopen("sync","w");
    fprintf(fsync,"sync");
    fclose(fsync);
}

/* ShowProg displays fundamental progress information. Basic run summary */

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
    myfile.open(writename);
    myfile << sLine;
    myfile.close();
}

/* * * * ***** Output Solutions * * * * * * * */
/** imode = 1   Output Summary Stats only ******/
/** imode = 2    Output Everything * * * * *****/

void OutputSummary(int puno,int spno,int R[],struct sspecies spec[],struct scost reserve,
                   int itn,char savename[],double misslevel,int imode)
{
    FILE *fp;  /* Imode = 1, REST output, Imode = 2, Arcview output */
    int i,ino=0,isp;
    double shortfall,connectiontemp,rMPM;
    char sZoneNames[1000],sZonePuCount[1000],sZoneCostNames[1000],sZoneCost[1000];

    CountPuZones2(sZoneNames,sZonePuCount,imode,puno,R);
    CostPuZones(sZoneCostNames,sZoneCost,imode,puno,R);

    if (itn==1)
        fp = fopen(savename,"w");
    else
        fp = fopen(savename,"a");
    if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

    /*** Ouput the Summary Statistics *****/
    for (i=0;i<=puno-1;i++)
        if (R[i] != 0)
            ino ++;
    isp = CountMissing(spno,spec,misslevel,&shortfall,&rMPM);
    for (i=0,connectiontemp = 0;i<puno;i++)
    {
        connectiontemp += ConnectionCost2Linear(i,R[i],pu,connections,R,1,0);
    } /* Find True (non modified) connection */

    if (itn==1)
    {
        if (imode > 1)
        {
            fprintf(fp,"\"Run Number\",\"Score\",\"Cost\",\"Planning Units\"%s%s",sZoneNames,sZoneCostNames);
            fprintf(fp,",\"Connection Strength\",\"Penalty\",\"Shortfall\",\"Missing_Values\",\"MPM\"\n");
        }
        else
        {
            fprintf(fp,"Run no.    Score      Cost   Planning Units  %s%s",sZoneNames,sZoneCostNames);
            fprintf(fp,"  Connection_Strength   Penalty  Shortfall Missing_Values MPM\n");
        }
    }
    if (imode > 1)
        fprintf(fp,"%i,%f,%f,%i%s%s,%f,%f,%f,%i,%f\n",
                   itn,reserve.total,reserve.cost,ino,sZonePuCount,sZoneCost,connectiontemp,reserve.penalty,shortfall,isp,rMPM);
    else
        fprintf(fp,"%-4i    %8.2f  %8.2f  %8i%s%s     %8.2f         %8.2f   %8.2f       %i   %8.2f\n",
                   itn,reserve.total,reserve.cost,ino,sZonePuCount,sZoneCost,connectiontemp,reserve.penalty,shortfall,isp,rMPM);
    fclose(fp);
    return;
} // OutputSummary

/* * * * ***** Scenario Output File * * * * * * * */
/*** OutputScenario ****/
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

void DumpAsymmetricConnectionFile(int puno,struct sconnections connections[],struct spustuff pu[],struct sfname fnames)
{
    #ifdef ASYMCON
    int i;
    FILE *fp;
    char *writename;
    struct sneighbour *p;

    writename = (char *) calloc(22 + strlen(fnames.outputdir)+2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debug_asymmetric_connection.csv");
    if ((fp = fopen(writename,"w"))==NULL)
    {
        ShowGenProg("Warning: Cannot create file %s",writename);
        free(writename);
    }
    free(writename);

    fprintf(fp,"idA,idB,connectionorigon\n");
    for (i=0;i<puno;i++)
    {
        for (p = connections[i].first;p;p=p->next)
            fprintf(fp,"%i,%i,%i,%lf\n",pu[i].id,pu[p->nbr].id,p->connectionorigon,p->cost);
    }

    fclose(fp);
    #endif
}

} // namespace marzone