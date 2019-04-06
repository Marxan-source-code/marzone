void DumpZoneNames(int iZoneCount,struct stringname Zones[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneNames.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneNames.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneNames file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,zonename\n");

     for (i=0;i<iZoneCount;i++)
     {
         fprintf(fp,"%d,%s\n",Zones[i].id,Zones[i].name);
     }

     fclose(fp);
}

void DumpCostNames(int iCostCount,struct stringname CostNames[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugCostNames.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugCostNames.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpCostNames file %s\n",writename);
     free(writename);

     fprintf(fp,"costid,costname\n");

     for (i=0;i<iCostCount;i++)
     {
         fprintf(fp,"%d,%s\n",CostNames[i].id,CostNames[i].name);
     }

     fclose(fp);
}

void DumpCostFieldNumber(int iFieldCount,int CostFieldNumber[],char *sFields,struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugCostFieldNumber.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugCostFieldNumber.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpCostFieldNumber file %s\n",writename);
     free(writename);

     fprintf(fp,"%s",sFields);
     for (i=0;i<iFieldCount;i++)
     {
         fprintf(fp,"%d",CostFieldNumber[i]);
         if (i != (iFieldCount-1))
            fprintf(fp,",");
     }
     fprintf(fp,"\n");

     fclose(fp);
}

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

void Dump_ZoneContrib(int puno,int spno,typesp spec[],int iZoneCount,double _ZoneContrib[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i,j,k;
     char messagebuffer[1000];

     #ifdef DEBUGTRACEFILE
     //AppendDebugTraceFile("Dump_ZoneContrib start\n");
     #endif

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
     }
     else
     {
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

     #ifdef DEBUGTRACEFILE
     //AppendDebugTraceFile("Dump_ZoneContrib end\n");
     #endif
}

void DumpZoneContrib(int iZoneContribCount,struct zonecontribstruct ZoneContrib[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneContrib.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneContrib.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneContrib file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,speciesid,fraction\n");

     for (i=0;i<iZoneContribCount;i++)
     {
         fprintf(fp,"%d,%d,%lf\n",ZoneContrib[i].zoneid,ZoneContrib[i].speciesid,ZoneContrib[i].fraction);
     }

     fclose(fp);
}

void DumpZoneContrib2(int iZoneContrib2Count,struct zonecontrib2struct ZoneContrib2[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneContrib2.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneContrib2.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneContrib2 file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,fraction\n");

     for (i=0;i<iZoneContrib2Count;i++)
     {
         fprintf(fp,"%d,%lf\n",ZoneContrib2[i].zoneid,ZoneContrib2[i].fraction);
     }

     fclose(fp);
}

void DumpZoneContrib3(int iZoneContrib3Count,struct zonecontrib3struct ZoneContrib3[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneContrib3.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneContrib3.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneContrib3 file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,puid,speciesid,fraction\n");

     for (i=0;i<iZoneContrib3Count;i++)
     {
         fprintf(fp,"%d,%d,%d,%lf\n",ZoneContrib3[i].zoneid,ZoneContrib3[i].puid,ZoneContrib3[i].speciesid,ZoneContrib3[i].fraction);
     }

     fclose(fp);
}

void Dump_ZoneTarget(int spno,int iZoneCount,struct _zonetargetstruct _ZoneTarget[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i,j;
     char debugbuffer[1000];

     #ifdef DEBUGTRACEFILE
     //sprintf(debugbuffer,"Dump_ZoneTarget start iZoneCount %i\n",iZoneCount);
     //AppendDebugTraceFile(debugbuffer);
     #endif

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
             #ifdef DEBUGTRACEFILE
             //sprintf(debugbuffer,"Dump_ZoneTarget j %i i %i index %i\n",j,i,(j*iZoneCount)+i);
             //AppendDebugTraceFile(debugbuffer);
             #endif

             fprintf(fp,"%lf,%i",_ZoneTarget[(j*iZoneCount)+i].target,_ZoneTarget[(j*iZoneCount)+i].occurrence);
             if (i != (iZoneCount-1))
                fprintf(fp,",");
         }
         fprintf(fp,"\n");
     }

     fclose(fp);

     #ifdef DEBUGTRACEFILE
     //AppendDebugTraceFile("Dump_ZoneTarget end\n");
     #endif
}

void DumpZoneTarget(int iZoneTargetCount,struct zonetargetstruct ZoneTarget[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneTarget.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneTarget.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneTarget file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,speciesid,target,targettype\n");

     for (i=0;i<iZoneTargetCount;i++)
     {
         fprintf(fp,"%d,%d,%lf,%d\n",ZoneTarget[i].zoneid,ZoneTarget[i].speciesid,ZoneTarget[i].target,ZoneTarget[i].targettype);
     }

     fclose(fp);
}

void DumpZoneTarget2(int iZoneTarget2Count,struct zonetarget2struct ZoneTarget2[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneTarget2.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneTarget2.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneTarget2 file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,target,targettype\n");

     for (i=0;i<iZoneTarget2Count;i++)
     {
         fprintf(fp,"%d,%lf,%d\n",ZoneTarget2[i].zoneid,ZoneTarget2[i].target,ZoneTarget2[i].targettype);
     }

     fclose(fp);
}

void Dump_ZoneCost(int iCostCount,int iZoneCount,double _ZoneCost[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i,j;

     #ifdef DEBUGTRACEFILE
     //AppendDebugTraceFile("Dump_ZoneCost start\n");
     #endif

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

     #ifdef DEBUGTRACEFILE
     //AppendDebugTraceFile("Dump_ZoneCost end\n");
     #endif
}

void DumpZoneCost(int iZoneCostCount,struct zonecoststruct ZoneCost[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneCost.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneCost.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneCost file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid,costid,fraction\n");

     for (i=0;i<iZoneCostCount;i++)
     {
         fprintf(fp,"%d,%d,%lf\n",ZoneCost[i].zoneid,ZoneCost[i].costid,ZoneCost[i].fraction);
     }

     fclose(fp);
}

void DumpPuLock(int iPuLockCount,struct pulockstruct PuLock[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugPuLock.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugPuLock.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpPuLock file %s\n",writename);
     free(writename);

     fprintf(fp,"puid,zoneid\n");

     for (i=0;i<iPuLockCount;i++)
     {
         fprintf(fp,"%d,%d\n",PuLock[i].puid,PuLock[i].zoneid);
     }

     fclose(fp);
}

void DumpPuZone(int iPuZoneCount,struct puzonestruct PuZone[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugPuZone.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugPuZone.csv");
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

void DumpRelConnectionCost(int iRelConnectionCostCount,struct relconnectioncoststruct RelConnectionCost[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugZoneConnectionCost.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugZoneConnectionCost.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpZoneConnectionCost file %s\n",writename);
     free(writename);

     fprintf(fp,"zoneid1,zoneid2,fraction\n");

     for (i=0;i<iRelConnectionCostCount;i++)
     {
         fprintf(fp,"%d,%d,%lf\n",RelConnectionCost[i].zoneid1,RelConnectionCost[i].zoneid2,RelConnectionCost[i].fraction);
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

     #ifdef DEBUGTRACEFILE
     //sprintf(debugbuffer,"DumpZoneSpec %i start\n",iMessage);
     //AppendDebugTraceFile(debugbuffer);
     #endif

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

     #ifdef DEBUGTRACEFILE
     //sprintf(debugbuffer,"DumpZoneSpec %i end\n",iMessage);
     //AppendDebugTraceFile(debugbuffer);
     #endif
}

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

void DumpFileNames(struct sfname fnames)
{
     FILE *fp;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debugFileNames.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debugFileNames.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create DumpFileNames file %s\n",writename);
     free(writename);

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
            if (spec[i].target > spec[i].amount)
            {
               rShortfall = spec[i].target - spec[i].amount;
            }

         fprintf(fp,"%i,%f,%f,%f",spec[i].name,spec[i].target,spec[i].amount,rShortfall);

         rTotalShortfall += rShortfall;
         rShortfall = 0;

         for (k=0;k<iZoneCount;k++)
         {
             iExistingArrayIndex = (i * iZoneCount) + k;

             if (_ZoneTarget[iExistingArrayIndex].target > 0)
                if (_ZoneTarget[iExistingArrayIndex].target > ZoneSpec[iExistingArrayIndex].amount)
                {
                   rShortfall = _ZoneTarget[iExistingArrayIndex].target - ZoneSpec[iExistingArrayIndex].amount;
                }

             fprintf(fp,",%f,%f,%f",_ZoneTarget[iExistingArrayIndex].target,ZoneSpec[iExistingArrayIndex].amount,rShortfall);

             rTotalShortfall += rShortfall;
             rShortfall = 0;
         }

         fprintf(fp,",%f\n",rTotalShortfall);
     }

     fclose(fp);
}

void DumpBinarySearchArrays(char *sName,struct sfname fnames, int puno, int spno, struct binsearch PULookup[],
                            struct binsearch SPLookup[])
{
     FILE *pufp, *specfp;
     int i;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug") + strlen(sName) + strlen("pu.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debug");
     strcat(writename,sName);
     strcat(writename,"pu.csv");
     pufp = fopen(writename,"w");
     if (pufp==NULL)
        ShowErrorMessage("cannot create BinarySearchArrays pu file %s\n",writename);
     free(writename);
     fputs("name,index\n",pufp);
     for (i=0;i<puno;i++)
     {
         fprintf(pufp,"%d,%d\n",PULookup[i].name,PULookup[i].index);
     }
     fclose(pufp);

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug") + strlen(sName) + strlen("spec.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"debug");
     strcat(writename,sName);
     strcat(writename,"spec.csv");
     specfp = fopen(writename,"w");
     if (specfp==NULL)
        ShowErrorMessage("cannot create BinarySearchArrays spec file %s\n",writename);
     free(writename);
     fputs("name,index\n",specfp);
     for (i=0;i<spno;i++)
     {
         fprintf(specfp,"%d,%d\n",SPLookup[i].name,SPLookup[i].index);
     }
     fclose(specfp);
}

void WriteStopErrorFile(char sMess[])
{
     //
     FILE* fsync;

     fsync = fopen("stop_error.txt","w");
     fprintf(fsync,"%s",sMess);
     fclose(fsync);
}

void WriteSlaveSyncFile(void)
{
     FILE* fsync;

     fsync = fopen("sync","w");
     fprintf(fsync,"sync");
     fclose(fsync);
}

/* ShowProg displays fundamental progress information. Basic run summary */
void ShowProg(char sMess[],...)
{
     va_list args;
     char theMessage[1000];
     
     if (iVerbosity > 0)
     {
        va_start(args,sMess);
        
        //printf("S_1\n");
        
        //vprintf(sMess,args);
        vsprintf(theMessage,sMess,args);
        
        //printf("S_2\n");
        
        printf("%s",theMessage);
                
        //printf("S_3\n");
        
        if (savelog)
        {
           //vfprintf(fsavelog,sMess,args);
           fprintf(fsavelog,"%s",theMessage);
           fflush(fsavelog);
        }
           
        //printf("S_4\n");
        
        va_end(args);
        
        //printf("S_5\n");
     }
} /* Show Progress Message */

#ifdef DEBUGTRACEFILE
void StartDebugTraceFile(void)
{
     FILE* fdebugtrace;


     if (iVerbosity > 2)
     {
        fdebugtrace = fopen(sDebugTraceFileName,"w");
        fflush(fdebugtrace);
        fclose(fdebugtrace);
     }
}
void AppendDebugTraceFile(char sMess[],...)
{
     FILE* fdebugtrace;

     if (iVerbosity > 2)
     {
        fdebugtrace = fopen(sDebugTraceFileName,"a");
        fprintf(fdebugtrace,"%s",sMess);
        fflush(fdebugtrace);
        fclose(fdebugtrace);
     }
}
#endif

void StartDebugFile(char sFileName[],char sHeader[],struct sfname fnames)
{
     FILE* fdebugtrace;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen(sFileName) + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,sFileName);
     fdebugtrace = fopen(writename,"w");
     free(writename);

     fprintf(fdebugtrace,"%s",sHeader);
     fflush(fdebugtrace);
     fclose(fdebugtrace);
}

void AppendDebugFile(char sFileName[],char sLine[],struct sfname fnames)
{
     FILE* fdebugtrace;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen(sFileName) + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,sFileName);
     fdebugtrace = fopen(writename,"a");
     free(writename);

     fprintf(fdebugtrace,"%s",sLine);
     fclose(fdebugtrace);
}

/* ShowGenProg displays a general progress message when verbosity > 1 */

void ShowGenProg(char sMess[],...)
{
     va_list args;

     if (iVerbosity > 1)
     {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
           vfprintf(fsavelog,sMess,args);
        va_end(args);
     }

}  /* Show General Progress Message */


/* ShowGenProgInfo displays a general progress with information
    message when verbosity > 2 */

void ShowGenProgInfo(char sMess[],...)
{
     va_list args;

     if (iVerbosity > 5)
     {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
           vfprintf(fsavelog,sMess,args);
        va_end(args);
     }
} /*Show General Progress Information Message */


/* ShowDetailedProgress shows detailed progress information
    message when verbosity > 3 */

void ShowDetProg(char sMess[],...)
{
     va_list args;

     if (iVerbosity > 5)
     {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
           vfprintf(fsavelog,sMess,args);
        va_end(args);
     }
} /* Show Detailed Progess Message */


/* ShowErrorMessage displays an error message. No matter what verbosity these are
    always displayed. The program is terminated following a ShowPauseExit*/

void ShowErrorMessage(char sMess[],...)
{
     extern jmp_buf jmpbuf;
     va_list args;

     va_start(args,sMess);
     vprintf(sMess,args);
     if (savelog)
        vfprintf(fsavelog,sMess,args);
     va_end(args);
     longjmp(jmpbuf,1);
} /* Show Error Message */

/* ShowWarningMessage displays a warning message no matter what verbosity level */

void ShowWarningMessage(char sMess[],...)
{
     va_list args;

     if (iVerbosity > 0)
     {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
           vfprintf(fsavelog,sMess,args);
        va_end(args);
     }
} /* Show Warning Message */

/* * * *  Set logged file. Also resets log file ****/

void SetLogFile(int my_savelog, char* my_savelogname)
{
     if (savelog)
     {
        fclose(fsavelog);
        free(savelogname);
     } /* close and delete old savelog info */

     savelog = my_savelog;

     if (savelog)
     {
        savelogname = calloc(strlen(my_savelogname)+1,sizeof(char));
        strcpy(savelogname,my_savelogname);
        /* Try to open file and complain if it don't work */
        if ((fsavelog = fopen(savelogname,"w"))==NULL)
        {
           free(savelogname);
           savelog = 0;
           ShowErrorMessage("Error: Cannot save to log file %s \n",savelogname);
        }  /* open failed */

        fprintf(fsavelog,"        %s \n\n   Marine Reserve Design with Zoning and Annealing\n\n",sVersionString);
        fprintf(fsavelog,"   Marxan with Zones coded by Matthew Watts\n");
        fprintf(fsavelog,"   Written by Ian Ball, Hugh Possingham and Matthew Watts\n\n");
        fprintf(fsavelog,"   Based on Marxan coded by Ian Ball, modified by Matthew Watts\n");
        fprintf(fsavelog,"   Written by Ian Ball and Hugh Possingham\n\n");
        fprintf(fsavelog,"%s\n%s\n%s\n\n",sIanBallEmail,sHughPossinghamEmail,sMattWattsEmail);
        fprintf(fsavelog,"   Marxan website\n\n");
        fprintf(fsavelog,"%s\n\n",sMarxanWebSite);

     } /* save log has just been turned on */
}  /* Set Log File */

/* creates an output file showing the planning unit richness and offset */
void DumpPU_richness_offset(int puno, struct spustuff PU[],struct sfname fnames)
{
     FILE *fp;
     int i;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("pu_richness_offset.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"pu_richness_offset.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create PU_richness_offset file %s\n",writename);
     free(writename);

     fputs("puindex,richness,offset\n",fp);
     for (i=0;i<puno;i++)
     {
         fprintf(fp,"%i,%i,%i\n",i,PU[i].richness,PU[i].offset);
     }

     fclose(fp);
}
/* creates an output file showing the planning unit richness and offset */

/* creates an output file from the loaded Sparse Matrix */
void DumpSparseMatrix(int iSMno, struct spu SM[],struct sfname fnames)
{
     FILE *fp;
     int i;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen("sm.csv") + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,"sm.csv");
     fp = fopen(writename,"w");
     if (fp==NULL)
        ShowErrorMessage("cannot create PUvSpecies file %s\n",writename);
     free(writename);

     fputs("amount,clump,spid\n",fp);
     for (i=0;i<iSMno;i++)
     {
         fprintf(fp,"%g,%i,%i\n",SM[i].amount,SM[i].clump,SM[i].spindex);
     }

     fclose(fp);
}
/* creates an output file from the loaded Sparse Matrix */

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

} /** Output Summary ***/

void OutputSpeciesData(int spno,struct sspecies spec[],char savename[])
{
     FILE *fp;
     int i;

     fp = fopen(savename,"w");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     fprintf(fp,"i,name,type,sname,target,prop,targetocc,spf,penalty,amount,occurrence,sepdistance,sepnum,separation,clumps,target2,richness,offset\n");

     for (i=0;i<spno;i++)
        fprintf(fp,"%i,%i,%i,%s,%f,%f,%i,%f,%f,%f,%i,%f,%i,%i,%i,%f,%i,%i"
                  ,i,spec[i].name,spec[i].type,spec[i].sname,spec[i].target,spec[i].prop,spec[i].targetocc
                  ,spec[i].spf,spec[i].penalty,spec[i].amount,spec[i].occurrence,spec[i].sepdistance
                  ,spec[i].sepnum,spec[i].separation,spec[i].clumps,spec[i].target2,spec[i].richness,spec[i].offset);

     fclose(fp);
} // Output Species Data

void OutputPenalty(int spno,struct sspecies spec[],char savename[],int iOutputType)
{
     FILE *fp;  // Imode = 1, REST output, Imode = 2, Arcview output
     int i;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     if (iOutputType == 3)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"    ");

     fprintf(fp,"spid%spenalty\n",sDelimiter);

     // Ouput the Summary Statistics
     for (i=0;i<spno;i++)
         fprintf(fp,"%i%s%g\n",spec[i].name,sDelimiter,spec[i].penalty);

     fclose(fp);
} // Output Penalty

void InitSolutionsMatrix(int puno,struct spustuff pu[],char savename[],int iOutputType,int iIncludeHeaders)
{
     FILE *fp;
     int i;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     if (iIncludeHeaders == 1)
     {
        if (iOutputType == 3)
           strcpy(sDelimiter,",");
        else
            strcpy(sDelimiter,"    ");

        fprintf(fp,"SolutionsMatrix");

        for (i=0;i<puno;i++)
            fprintf(fp,"%s%i",sDelimiter,pu[i].id);

        fprintf(fp,"\n");
     }

     fclose(fp);
}

void AppendSolutionsMatrix(int iRun,int puno,int R[],char savename[],int iOutputType,int iIncludeHeaders)
{
     FILE *fp;
     int i,j;
     char sDelimiter[20];

     fp = fopen(savename,"a");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     if (iOutputType == 3)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"    ");

     for (i=1;i<=iZoneCount;i++)
     {
         if (iIncludeHeaders == 1)
         {
            fprintf(fp,"Z%iS%i%s",i,iRun,sDelimiter);
         }

         for (j=0;j<puno;j++)
         {
             if (j > 0)
                fprintf(fp,"%s",sDelimiter);

             if (R[j] == i)
                fprintf(fp,"1");
             else
                 fprintf(fp,"0");
         }

         fprintf(fp,"\n");
     }

     fclose(fp);
}

void AppendSolutionsMatrixZone(int iZone,int iRun,int puno,int R[],char savename[],int iOutputType,int iIncludeHeaders)
{
     FILE *fp;
     int j;
     char sDelimiter[20];

     fp = fopen(savename,"a");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     if (iOutputType == 3)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"    ");

     if (iIncludeHeaders == 1)
     {
        fprintf(fp,"S%i%s",iRun,sDelimiter);
     }

     for (j=0;j<puno;j++)
     {
         if (j > 0)
            fprintf(fp,"%s",sDelimiter);

         if (R[j] == iZone)
            fprintf(fp,"1");
         else
             fprintf(fp,"0");
     }

     fprintf(fp,"\n");

     fclose(fp);
}

/*** Output A Solution ***/
void OutputSolution(int puno,int R[],struct spustuff pu[],char savename[],int imode)
{
     FILE *fp;  /* Imode = 1, REST output, Imode = 2, Arcview output */
     int i;
     fp = fopen(savename,"w");
     if (!fp)  ShowErrorMessage("Cannot save output to %s \n",savename);

     if (imode > 1)
        fprintf(fp,"planning_unit,zone\n");
     //for (i=0;i<puno;i++)
     for (i=puno-1;i>-1;i--)
         if (R[i] != 0)
         {
            fprintf(fp,"%i",pu[i].id);
            if (imode > 1)
               fprintf(fp,",%i",R[i]);
            fprintf(fp,"\n");
         }
     fclose(fp);
} /* Output Solution */

/* * * * ***** Scenario Output File * * * * * * * */
/*** OutputScenario ****/
void OutputScenario(int puno,int spno,double prop,
                    struct sanneal anneal,int seedinit,int repeats,int clumptype,
                    int runopts,int heurotype,double costthresh, double tpf1, double tpf2,
                    char savename[])
{
     FILE *fp;
     char temp[100];

     fp = fopen(savename,"w");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     fprintf(fp,"Number of Planning Units %i\n",puno);
     fprintf(fp,"Number of Conservation Values %i\n",spno);
     fprintf(fp,"Number of Zones %i\n",iZoneCount);
     fprintf(fp,"Number of Costs %i\n",iCostCount);
     fprintf(fp,"Starting proportion %.2f\n",prop);
     switch (clumptype)
     {
            case 0:strcpy(temp,"Clumping - default step function");break;
            case 1: strcpy(temp,"Clumping - two level step function.");break;
            case 2: strcpy(temp,"Clumping - rising benefit function");break;
     }
     fprintf(fp,"%s\n",temp);

     /* Use character array here and set up the name of the algorithm used */
     switch (runopts)
     {
            case 0: strcpy(temp,"Annealing and Heuristic");break;
            case 1: strcpy(temp,"Annealing and Iterative Improvement");break;
            case 2: strcpy(temp,"Annealing and Both");break;
            case 3: strcpy(temp,"Heuristic only");break;
            case 4: strcpy(temp,"Iterative Improvement only");break;
            case 5: strcpy(temp,"Heuristic and Iterative Improvement");
     }
     fprintf(fp,"Algorithm Used :%s\n",temp);
     if (runopts == 0 || runopts == 3 || runopts == 5)
     {
        switch (heurotype)
        {
               case 0: strcpy(temp,"Richness");break;
               case 1: strcpy(temp,"Greedy");break;
               case 2: strcpy(temp,"Maximum Rarity");break;
               case 3: strcpy(temp,"Best Rarity");break;
               case 4: strcpy(temp,"Average Rarity");break;
               case 5: strcpy(temp,"Summation Rarity");break;
               case 6: strcpy(temp,"Product Irreplaceability");break;
               case 7: strcpy(temp,"Summation Irreplaceability");break;
               default: strcpy(temp,"Unkown Heuristic Type");
        }
        fprintf(fp,"Heuristic type : %s\n",temp);
     }
     else
         fprintf(fp,"No Heuristic used \n");

     if (runopts <= 2)
     {
        fprintf(fp,"Number of iterations %i\n",anneal.iterations);
        if (anneal.Tinit >= 0)
        {
           fprintf(fp,"Initial temperature %.2f\n",anneal.Tinit);
           fprintf(fp,"Cooling factor %.6f\n",anneal.Tcool);
        }
        else
        {
            fprintf(fp,"Initial temperature set adaptively\n");
            fprintf(fp,"Cooling factor set adaptively\n");
        }
        fprintf(fp,"Number of temperature decreases %i\n\n",anneal.Titns);
     }
     else
     {
         fprintf(fp,"Number of iterations N/A\nInitial temperature N/A\nCooling Factor N/A\n");
         fprintf(fp,"Number of temperature decreases N/A\n\n");
     }

     if (costthresh)
     {
        fprintf(fp,"Cost Threshold Enabled: %lf\n",costthresh);
        fprintf(fp,"Threshold penalty factor A %.2f\n",tpf1);
        fprintf(fp,"Threshold penalty factor B %.2f\n\n",tpf2);
     }
     else
     {
         fprintf(fp,"Cost Threshold Disabled\nThreshold penalty factor A N/A\n");
         fprintf(fp,"Threshold penalty factor B N/A\n\n");
     }

     fprintf(fp,"Random Seed %i\n",seedinit);
     fprintf(fp,"Number of runs %i\n",repeats);
     fclose(fp);
}  /*** Output Scenario ****/

// Output Feature report (missing values report)
void OutputFeatures(int spno,struct sspecies spec[],char savename[],int imode,double misslevel)
{
     FILE *fp; // Imode = 1, Tab Delimitted Text output, Imode = 2, Arcview output, Imode = 3, Arcview output with CSV file extension
     int i,isp, iZoneArrayIndex;
     char temp[4];
     double rTarget, rAmountHeld, rOccurrenceTarget, rOccurrencesHeld, rMPM, rTestMPM;

     fp = fopen(savename,"w");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     if (imode > 1)
     {
        fprintf(fp,"\"Feature\",\"Feature Name\",\"Target\"");
        fprintf(fp,",\"Total Amount\",\"Contributing Amount Held\",\"Occurrence Target \",\"Occurrences Held\"");
        fprintf(fp,",\"Target Met\"");
        for (i=0;i<iZoneCount;i++)
        {
            fprintf(fp,",\"Target %s\",\"Amount Held %s\",\"Contributing Amount Held %s\",\"Occurrence Target %s\"",
                         Zones[i].name,Zones[i].name,Zones[i].name,Zones[i].name);
            fprintf(fp,",\"Occurrences Held %s\",\"Target Met %s\"",Zones[i].name,Zones[i].name);
        }
        fprintf(fp,",MPM\n");
     }
     else
     {
         fprintf(fp,"Feature\tFeature Name\tTarget");
         fprintf(fp,"\tAmount Held\tContributing Amount Held\tOccurrence Target \tOccurrences Held");
         fprintf(fp,"\tTarget Met");
         for (i=0;i<iZoneCount;i++)
         {
             fprintf(fp,"\tTarget %s\tAmount Held %s\tContributing Amount Held %s\tOccurrence Target %s",
                        Zones[i].name,Zones[i].name,Zones[i].name,Zones[i].name);
             fprintf(fp,"\tOccurrences Held %s\tTarget Met %s",Zones[i].name,Zones[i].name);
         }
         fprintf(fp,"\tMPM\n");
     }

     for (isp=0;isp<spno;isp++)
     {
         rMPM = 1;

         if (imode > 1)
         {
            fprintf(fp,"%i,%s,%lf",spec[isp].name,
                                   spec[isp].sname,
                                   spec[isp].target);
            fprintf(fp,",%lf,%lf,%i,%i",TotalAreas[isp],
                                        spec[isp].amount,
                                        spec[isp].targetocc,
                                        spec[isp].occurrence);
            strcpy(temp,"");  // use MISSLEVEL when computing target met
            if (spec[isp].target)
            {
               strcpy(temp,"yes");
               if (spec[isp].amount/spec[isp].target < misslevel)
                  strcpy(temp,"no");

               rTestMPM = spec[isp].amount/spec[isp].target;
               if (rTestMPM < rMPM)
                  rMPM = rTestMPM;
            }
            if (spec[isp].targetocc)
            {
               strcpy(temp,"yes");
               if (spec[isp].occurrence/spec[isp].targetocc < misslevel)
                  strcpy(temp,"no");

               rTestMPM = spec[isp].occurrence/spec[isp].targetocc;
               if (rTestMPM < rMPM)
                  rMPM = rTestMPM;
            }
            if (spec[isp].sepnum)
            {
               strcpy(temp,"yes");
               if (spec[isp].separation/spec[isp].sepnum < misslevel)
                  strcpy(temp,"no");
            }
            fprintf(fp,",%s",temp);
            for (i=0;i<iZoneCount;i++)
            {
                iZoneArrayIndex = (isp * iZoneCount) + i;
                rTarget = _ZoneTarget[iZoneArrayIndex].target;
                rAmountHeld = ZoneSpec[iZoneArrayIndex].amount;
                rOccurrenceTarget = _ZoneTarget[iZoneArrayIndex].occurrence;
                rOccurrencesHeld = ZoneSpec[iZoneArrayIndex].occurrence;
                fprintf(fp,",%lf,%lf,%lf,%lf,%lf",rTarget,
                                                  rAmountHeld,
                                                  rAmountHeld * _ZoneContrib[(isp * iZoneCount) + i],
                                                  rOccurrenceTarget,
                                                  rOccurrencesHeld);
                strcpy(temp,"");  // use MISSLEVEL when computing target met
                if (rTarget)
                {
                   strcpy(temp,"yes");  // use MISSLEVEL when computing target met
                   if (rAmountHeld/rTarget < misslevel)
                      strcpy(temp,"no");

                   rTestMPM = rAmountHeld/rTarget;
                   if (rTestMPM < rMPM)
                      rMPM = rTestMPM;
                }
                if (rOccurrenceTarget)
                {
                   strcpy(temp,"yes");  // use MISSLEVEL when computing target met
                   if (rOccurrencesHeld/rOccurrenceTarget < misslevel)
                      strcpy(temp,"no");

                   rTestMPM = rOccurrencesHeld/rOccurrenceTarget;
                   if (rTestMPM < rMPM)
                      rMPM = rTestMPM;
                }
                fprintf(fp,",%s",temp);
            }
            fprintf(fp,",%lf\n",rMPM);
         }
         else
         {
             fprintf(fp,"%i\t%s\t",spec[isp].name,
                                   spec[isp].sname);
             fprintf(fp,"%lf\t%lf\t%lf\t%i\t%i\t",spec[isp].target,
                                                  spec[isp].amount,
                                                  spec[isp].amount,
                                                  spec[isp].targetocc,
                                                  spec[isp].occurrence);
             strcpy(temp,"");  // use MISSLEVEL when computing target met
             if (spec[isp].target)
             {
                strcpy(temp,"yes");
                if (spec[isp].amount/spec[isp].target < misslevel)
                   strcpy(temp,"no");

                rTestMPM = spec[isp].amount/spec[isp].target;
                if (rTestMPM < rMPM)
                   rMPM = rTestMPM;
             }
             if (spec[isp].targetocc)
             {
                strcpy(temp,"yes");
                if (spec[isp].occurrence/spec[isp].targetocc < misslevel)
                   strcpy(temp,"no");

                rTestMPM = spec[isp].occurrence/spec[isp].targetocc;
                if (rTestMPM < rMPM)
                   rMPM = rTestMPM;
             }
             fprintf(fp,"\t%s",temp);
             for (i=0;i<iZoneCount;i++)
             {
                 iZoneArrayIndex = (isp * iZoneCount) + i;
                 rTarget = _ZoneTarget[iZoneArrayIndex].target;
                 rAmountHeld = ZoneSpec[iZoneArrayIndex].amount;
                 rOccurrenceTarget = _ZoneTarget[iZoneArrayIndex].occurrence;
                 rOccurrencesHeld = ZoneSpec[iZoneArrayIndex].occurrence;
                 fprintf(fp,"\t%lf\t%lf\t%lf\t%lf\t%lf",rTarget,
                                                        rAmountHeld,
                                                        rAmountHeld,
                                                        rOccurrenceTarget,
                                                        rOccurrencesHeld);
                 strcpy(temp,"");  // use MISSLEVEL when computing target met
                 if (rTarget)
                 {
                    strcpy(temp,"yes");
                    if (rAmountHeld/rTarget < misslevel)
                       strcpy(temp,"no");

                    rTestMPM = rAmountHeld/rTarget;
                    if (rTestMPM < rMPM)
                       rMPM = rTestMPM;
                 }
                 if (rOccurrenceTarget)
                 {
                    strcpy(temp,"yes");
                    if (rOccurrencesHeld/rOccurrenceTarget < misslevel)
                       strcpy(temp,"no");

                    rTestMPM = rOccurrencesHeld/rOccurrenceTarget;
                    if (rTestMPM < rMPM)
                       rMPM = rTestMPM;
                 }
                 fprintf(fp,"\t%s",temp);
             }
             fprintf(fp,"\t%lf\n",rMPM);
         }
     }

     fclose(fp);
} // Output Feature report (missing values report)

/*** Output A Solution ***/
void OutputSumSoln(int puno,int sumsoln[],int ZoneSumSoln[],int R[],struct spustuff pu[],char savename[],int imode)
{
     FILE *fp;  /* Imode = 1, REST output, Imode = 2, Arcview output */
     int i,j;

     fp = fopen(savename,"w");
     if (!fp)
        ShowErrorMessage("Cannot save output to %s \n",savename);

     if (imode > 1)
     {
        fprintf(fp,"\"planning unit\",\"number\"");
        for (j=0;j<iZoneCount;j++)
            fprintf(fp,",\"%s\"",Zones[j].name);
        fprintf(fp,"\n");
     }

     for (i=0;i<puno;i++)
         if (imode > 1)
         {
            fprintf(fp,"%i,%i",pu[i].id,sumsoln[i]);
            for (j=0;j<iZoneCount;j++)
                fprintf(fp,",%i",ZoneSumSoln[(puno * j) + i]);
            fprintf(fp,"\n");
         }
         else
         {
             fprintf(fp,"%i %i",pu[i].id,sumsoln[i]);
             for (j=0;j<iZoneCount;j++)
                 fprintf(fp," %i",ZoneSumSoln[(puno * j) + i]);
             fprintf(fp,"\n");
         }
     fclose(fp);
} /* Output Solution */

void OutputZoneConnectivitySum(int puno,int R[],char savename[],int imode)
{
     FILE *fp;
     int i,j;
     char sDelimiter[20];
     double **ZCS;

     ZCS = (double **) calloc(iZoneCount,sizeof(double *));
     for (i=0;i<iZoneCount;i++)
         ZCS[i] = (double *) calloc(iZoneCount,sizeof(double));
     for (i=0;i<iZoneCount;i++)
         for (j=0;j<iZoneCount;j++)
             ZCS[i][j] = 0;

     ComputeZoneConnectivitySum(ZCS,puno,R);

     if (imode > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"   ");

     fp = fopen(savename,"w");
     if (!fp)  ShowErrorMessage("Cannot save output to %s \n",savename);

     // write header row
     if (imode > 1)
        fprintf(fp,"\"Zone_Connectivity_Sum\"");
     for (i=0;i<iZoneCount;i++)
        fprintf(fp,",\"%s\"",Zones[i].name);
     fprintf(fp,"\n");

     // write a data row for each zone
     for (i=0;i<iZoneCount;i++)
     {
         fprintf(fp,"\"%s\"",Zones[i].name);
         for (j=0;j<iZoneCount;j++)
             fprintf(fp,",%f",ZCS[i][j]);
         fprintf(fp,"\n");
     }

     fclose(fp);

     for (i=0;i<iZoneCount;i++)
         free(ZCS[i]);
     free(ZCS);
}

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

void OutputTotalAreas(int puno,int spno,struct spustuff pu[],struct sspecies spec[],struct spu SM[],char savename[],int iOutputType)
{
     int ipu, i, ism, isp, *TotalOccurrences, *TO_2, *TO_3;
     double *TotalAreas, *TA_2, *TA_3;
     FILE* TotalAreasFile;
     char sDelimiter[20];

     TotalOccurrences = (int *) calloc(spno,sizeof(int));
     TO_2 = (int *) calloc(spno,sizeof(int));
     TO_3 = (int *) calloc(spno,sizeof(int));
     TotalAreas = (double *) calloc(spno,sizeof(double));
     TA_2 = (double *) calloc(spno,sizeof(double));
     TA_3 = (double *) calloc(spno,sizeof(double));

     for (i=0;i<spno;i++)
     {
         TotalAreas[i] = 0;
         TA_2[i] = 0;
         TA_3[i] = 0;
     }

     for (ipu=0;ipu<puno;ipu++)
         if (pu[ipu].richness)
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;

                TotalOccurrences[isp]++;
                TotalAreas[isp] += SM[ism].amount;

                if (pu[ipu].status == 2)
                {
                   TO_2[isp]++;
                   TA_2[isp] += SM[ism].amount;
                }

                if (pu[ipu].status == 3)
                {
                   TO_3[isp]++;
                   TA_3[isp] += SM[ism].amount;
                }
            }

     if (iOutputType > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"    ");

     TotalAreasFile = fopen(savename,"w");
     fprintf(TotalAreasFile,"spname%stotalarea%sreservedarea%sexcludedarea%stargetarea%stotalocc%sreservedocc%sexcludedocc%stargetocc\n"
                           ,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter);
     for (i=0;i<spno;i++)
         fprintf(TotalAreasFile,"%i%s%g%s%g%s%g%s%g%s%i%s%i%s%i%s%i\n"
                               ,spec[i].name,sDelimiter,TotalAreas[i],sDelimiter,TA_2[i],sDelimiter,TA_3[i],sDelimiter
                               ,spec[i].target,sDelimiter,TotalOccurrences[i],sDelimiter,TO_2[i],sDelimiter,TO_3[i],sDelimiter,spec[i].targetocc);
     fclose(TotalAreasFile);

     free(TotalOccurrences);
     free(TO_2);
     free(TO_3);
     free(TotalAreas);
     free(TA_2);
     free(TA_3);
}
