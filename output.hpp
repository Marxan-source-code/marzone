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


} // namespace marzone