#pragma once

#include <string>

#include "common.hpp"
#include "logger.hpp"
#include "zones.hpp"

namespace marzone {
// This file is deprecated. If a function is commented out here, it has yet to be rewritten.

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
        // already migrated.
    }
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

/*
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
*/

/* * * * ***** Scenario Output File * * * * * * * */
/*** OutputScenario ****/
/*
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
*/

} // namespace marzone