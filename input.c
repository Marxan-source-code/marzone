void LoadZones(int *iZoneCount,struct stringname *Zones[],struct sfname fnames)
{
    FILE *fp;
    char *readname, *writename, sLine[5000], *sVarVal;
    int i = 0, iLineCount = 0,j;
    char debugbuffer[1000];

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonesname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonesname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open Zones file %s\n",readname);

    // count the number of records
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the zones array
    *iZoneCount = iLineCount;
    *Zones = (struct stringname *) calloc(iLineCount,sizeof(struct stringname));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer id from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*Zones)[i].id);

        // read the string name from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
           
        // remove last character from string if it is carriage return
        if (sVarVal[strlen(sVarVal)-1]==13)
        {
            sVarVal[strlen(sVarVal)-1]=0;
        }
        strcpy((*Zones)[i].name,sVarVal);
        
        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadCostNames(int *iCostCount,struct stringname *CostNames[],struct sfname fnames)
{

    FILE *fp;
    char *readname, *writename, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.costsname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.costsname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open CostNames file %s\n",readname);

    // count the number of records
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the CostNames array
    *iCostCount = iLineCount;
    *CostNames = (struct stringname *) calloc(iLineCount,sizeof(struct stringname));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*CostNames)[i].id);

        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        strcpy((*CostNames)[i].name,sVarVal);
        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadCostValues(int iCostCount,double **CostValues,struct stringname CostNames[],struct sfname fnames,int puno)
{

    FILE *fp;
    char *readname, *writename, sLine[1000], sFields[1000], *sVarVal, sField[1000];
    int i, iLine = 0, /*iLineCount = 0,*/ iFieldCount = 0, *CostFieldNumber, iCostIndex;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.puname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.puname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open Cost Values file %s\n",readname);

    fgets(sLine,999,fp);
    strcpy(sFields,sLine);
    // count the number of fields
    strtok(sFields," ,;:^*\"/|\t\'\\\n");
    iFieldCount++;
    while ((sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n")) != NULL)
        iFieldCount++;
    // determine the field number for each of our costs
    i = 0;
    strcpy(sFields,sLine);
    CostFieldNumber = (int *) calloc(iFieldCount,sizeof(int));
    sVarVal = strtok(sFields," ,;:^*\"/|\t\'\\\n");
    CostFieldNumber[i] = rtnCostIndex(iCostCount,CostNames,sVarVal);
    while ((sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n")) != NULL)
    {
        i++;
        CostFieldNumber[i] = rtnCostIndex(iCostCount,CostNames,sVarVal);
    }
    #ifdef DEBUGTRACEFILE
    if (iVerbosity > 3)
        DumpCostFieldNumber(iFieldCount,CostFieldNumber,sLine,fnames);
    #endif

    // load the cost value data to an array from the pu.dat file
    while (fgets(sLine,999,fp))
    {
        i = 0;
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        if (CostFieldNumber[i] != -1)
            sscanf(sVarVal, "%lf", &CostValues[puno-iLine-1][CostFieldNumber[i]]);
        while ((sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n")) != NULL)
        {
            i++;
            if (CostFieldNumber[i] != -1)
                sscanf(sVarVal,"%lf",&CostValues[puno-iLine-1][CostFieldNumber[i]]);
        }
        iLine++;
    }
    fclose(fp);

    free(readname);
    free(CostFieldNumber);
}

void LoadZoneContrib(int *iZoneContribCount,struct zonecontribstruct *ZoneContrib[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonecontribname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonecontribname);
    fp = fopen(readname,"r");
    if (fp==NULL)
         ShowErrorMessage("cannot open ZoneContrib file %s\n",readname);

    // count the number of record
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the ZoneContrib array
    *iZoneContribCount = iLineCount;
    *ZoneContrib = (struct zonecontribstruct *) calloc(iLineCount,sizeof(struct zonecontribstruct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneContrib)[i].zoneid);

        // read the integer speciesid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneContrib)[i].speciesid);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*ZoneContrib)[i].fraction);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadZoneContrib2(int *iZoneContrib2Count,struct zonecontrib2struct *ZoneContrib2[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonecontrib2name) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonecontrib2name);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open ZoneContrib2 file %s\n",readname);

    // count the number of record
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the ZoneContrib array
    *iZoneContrib2Count = iLineCount;
    *ZoneContrib2 = (struct zonecontrib2struct *) calloc(iLineCount,sizeof(struct zonecontrib2struct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneContrib2)[i].zoneid);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*ZoneContrib2)[i].fraction);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadZoneContrib3(int *iZoneContrib3Count,struct zonecontrib3struct *ZoneContrib3[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonecontrib3name) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonecontrib3name);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open ZoneContrib3 file %s\n",readname);

    // count the number of record
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the ZoneContrib array
    *iZoneContrib3Count = iLineCount;
    *ZoneContrib3 = (struct zonecontrib3struct *) calloc(iLineCount,sizeof(struct zonecontrib3struct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneContrib3)[i].zoneid);

        // read the integer puid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneContrib3)[i].puid);

        // read the integer speciesid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneContrib3)[i].speciesid);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*ZoneContrib3)[i].fraction);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadZoneTarget(int *iZoneTargetCount,struct zonetargetstruct *ZoneTarget[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonetargetname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonetargetname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open ZoneTarget file %s\n",readname);

    // count the number of record
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the ZoneTarget array
    *iZoneTargetCount = iLineCount;
    *ZoneTarget = (struct zonetargetstruct *) calloc(iLineCount,sizeof(struct zonetargetstruct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    // does the table have a 4th column called prop?

    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneTarget)[i].zoneid);

        // read the integer speciesid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneTarget)[i].speciesid);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*ZoneTarget)[i].target);

        // read the integer targettype from this line if it exists
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        if (sVarVal == NULL)
            (*ZoneTarget)[i].targettype = 0;
        else
            sscanf(sVarVal, "%d", &(*ZoneTarget)[i].targettype);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadZoneTarget2(int *iZoneTarget2Count,struct zonetarget2struct *ZoneTarget2[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonetarget2name) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonetarget2name);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open ZoneTarget2 file %s\n",readname);

    // count the number of record
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the ZoneTarget array
    *iZoneTarget2Count = iLineCount;
    *ZoneTarget2 = (struct zonetarget2struct *) calloc(iLineCount,sizeof(struct zonetarget2struct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    // does the table have a 4th column called prop?

    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneTarget2)[i].zoneid);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*ZoneTarget2)[i].target);

        // read the integer targettype from this line if it exists
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        if (sVarVal == NULL)
            (*ZoneTarget2)[i].targettype = 0;
        else
            sscanf(sVarVal, "%i", &(*ZoneTarget2)[i].targettype);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadZoneCost(int *iZoneCostCount,struct zonecoststruct *ZoneCost[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.zonecostname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.zonecostname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open ZoneCost file %s\n",readname);

    // count the number of records
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the ZoneCost array
    *iZoneCostCount = iLineCount;
    *ZoneCost = (struct zonecoststruct *) calloc(iLineCount,sizeof(struct zonecoststruct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneCost)[i].zoneid);

        // read the integer costid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*ZoneCost)[i].costid);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*ZoneCost)[i].fraction);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadPuLock(int *iPuLockCount,struct pulockstruct *PuLock[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.pulockname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.pulockname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open PuLock file %s\n",readname);

    // count the number of records
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the PuLock array
    *iPuLockCount = iLineCount;
    *PuLock = (struct pulockstruct *) calloc(iLineCount,sizeof(struct pulockstruct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer puid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*PuLock)[i].puid);

        // read the integer zoneid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*PuLock)[i].zoneid);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadPuZone(int *iPuZoneCount,struct puzonestruct *PuZone[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.puzonename) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.puzonename);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open PuZone file %s\n",readname);

    // count the number of records
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the PuLock array
    *iPuZoneCount = iLineCount;
    *PuZone = (struct puzonestruct *) calloc(iLineCount,sizeof(struct puzonestruct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer puid from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*PuZone)[i].puid);

        // read the integer zoneid from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*PuZone)[i].zoneid);

        i++;
    }
    fclose(fp);

    free(readname);
}

void LoadRelConnectionCost(int *iRelConnectionCostCount,struct relconnectioncoststruct *RelConnectionCost[],struct sfname fnames)
{

    FILE *fp;
    char *readname, sLine[1000], *sVarVal;
    int i = 0, iLineCount = 0, id;

    readname = (char *) calloc(strlen(fnames.inputdir) + strlen(fnames.relconnectioncostname) + 2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.relconnectioncostname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("cannot open RelConnectionCost file %s\n",readname);

    // count the number of records
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
        iLineCount++;
    fclose(fp);

    // create the RelConnectionCost array
    *iRelConnectionCostCount = iLineCount;
    *RelConnectionCost = (struct relconnectioncoststruct *) calloc(iLineCount,sizeof(struct relconnectioncoststruct));

    // load the data to an array
    fp = fopen(readname,"r");
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        // read the integer zoneid1 from this line
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*RelConnectionCost)[i].zoneid1);

        // read the integer zoneid2 from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%d", &(*RelConnectionCost)[i].zoneid2);

        // read the double fraction from this line
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal, "%lf", &(*RelConnectionCost)[i].fraction);

        i++;
    }
    fclose(fp);

    free(readname);
}



/* * * * ** Read Planning Unit Costs * * * * ******/
int ReadPUCosts(int puno,struct spustuff pu[],struct binsearch PULookup[],int verbose,char indir[])
{
    FILE *fp;
    int i,n;
    double ftemp;
    char readname[1000];

    if (indir[0] != '0')
        strcpy(readname,indir);
    strcat(readname,"cost.dat");

    fp = fopen(readname,"r");
    if (fp==NULL)
    {
        if (verbose > 1)
            ShowWarningMessage("File %s not found. All P.U.s set to cost of 1\n",readname);
        for (i=0;i<puno;i++)
            pu[i].cost = 1;
        return(0);
    } /* no PUcost file */

    /** Clear the cost structure **/
    i = 0;
    while (fscanf(fp,"%d %lf",&n,&ftemp)==2)
    {
        n = FastPUIDtoPUINDEX(puno,n,PULookup);
        if (n<0 || n>=puno)
            ShowErrorMessage("Invalid planning unit number %d \n",n);
        pu[n].cost += ftemp;
        if (ftemp == 0)
            pu[n].cost = delta; /* Don't like zero cost. This is temporary line */
        i++;
    } /** Found another valid looking line **/
    fclose(fp);
    return(i);
} // ReadPUCosts

/****** Read Species Information Data  ****/
int ReadSpeciesData(int *spno,struct sspecies *spec[],int *aggexist, int *sepexist,char indir[])
{
    FILE *fp;
    int i,n=0;
    struct slink{struct sspecies stemp;struct slink *next;} *head = NULL,*newlink;
    struct sspecies spectemp;
    char readname[1000];
    char speciesname[1000];

    if (indir[0] != '0')
        strcpy(readname,indir);
    strcat(readname,"species.dat");

    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("Species file %s has not been found.\nAborting Program.",readname);

    /** Clear important species stats **/
    while (fscanf(fp,"%d %d %lf %lf %lf %lf %s",&spectemp.name,&i,&spectemp.target,&spectemp.spf,&spectemp.target2,&spectemp.sepdistance,speciesname)==7)
    {
        newlink = (struct slink *) malloc(sizeof(struct slink));
        newlink->stemp.name = spectemp.name;
        newlink->stemp.target = spectemp.target;
        newlink->stemp.spf = spectemp.spf;
        newlink->stemp.target2 = spectemp.target2;
        newlink->stemp.sepdistance = spectemp.sepdistance;
        if (spectemp.target2)
            *aggexist = 1;
        if (spectemp.sepdistance)
        {
            *sepexist = 1;
            newlink->stemp.separation = 1;
        }
        newlink->stemp.sname = (char *) calloc(strlen(speciesname)+1,sizeof(char));
        strcpy(newlink->stemp.sname,speciesname);
        n++;
        newlink->next = head;
        head = newlink;
    }
    fclose(fp);

    /* Now do as Name.dat in forming species array */
    *spno = n;
    *spec = (struct sspecies *) calloc(*spno,sizeof(struct sspecies));
    /* put each link into namelist and free it */
    n = 0;
    while (head)
    {
        (* spec)[n].name = head->stemp.name;
        (* spec)[n].target = head->stemp.target;
        (* spec)[n].spf = head->stemp.spf;
        (* spec)[n].target2 = head->stemp.target2;
        (* spec)[n].sepdistance = head->stemp.sepdistance;
        (* spec)[n].separation = head->stemp.separation;
        (* spec)[n].sname = (char *) calloc(strlen(head->stemp.sname)+1,sizeof(char));
        (* spec)[n].richness = 0;
        strcpy((* spec)[n].sname,head->stemp.sname);
        (* spec)[n].targetocc = 0;
        n++;
        newlink = head;
        head = head->next;
        free(newlink->stemp.sname);
        free(newlink);
    }

    return(n);
} // ReadSpeciesData


/***** Planning Unit Information File * * * * ******/
/** Note status = 0     Not in Reserve
        Status = 1        In Reserve
        Status = 2        In Reserve, non-removable
        Status = 3        Not in Reserve, can not be added **/
int ReadPUFile(int puno,struct spustuff pu[],struct binsearch PULookup[],int verbose,char indir[])
{
    FILE *fp;
    int i=0,n,ireserved =0,iproscribed = 0,iinit= 0,idup=0;
    int itemp;
    char readname[1000];

    if (indir[0] != '0')
       strcpy(readname,indir);
    strcat(readname,"pustat.dat");

    fp = fopen(readname,"r");
    if (fp==NULL)
    {
        ShowGenProg("No PU Status file \n");
        return(0);
    }

    while (fscanf(fp,"%d %d",&n,&itemp)==2)
    {
        i++;
        n = FastPUIDtoPUINDEX(puno,n,PULookup);
        if (n<0 || n>=puno)
            ShowErrorMessage("planning unit id out of bounds %d \n",n);
        if (pu[n].status)
            idup++;
        pu[n].status = itemp;
        if (itemp == 1)
            iinit++;
        if (itemp == 2)
            ireserved++;
        if (itemp == 3)
            iproscribed++;
    }
    fclose(fp);

    if (verbose > 1)
    {
        if (iinit || ireserved || iproscribed)
            ShowGenProg("Reserve Status:");
        if (iinit)
            ShowGenProg(" initial reserve %i.",iinit);
        if (ireserved)
            ShowGenProg(" iremovable  %i.",ireserved);
        if (iproscribed)
            ShowGenProg(" Not available %i.",iproscribed);
        if (idup)
            ShowGenProg("There were %i planning units duplicated",idup);
        ShowGenProg("\n");
    }
    return(i);
} // ReadPUFile

int ReadPUXYfile(int puno,struct spustuff pu[],struct binsearch PULookup[],char indir[])
{
    FILE *fp;
    int i = 0,n;
    double x,y;
    char readname[1000];

    if (indir[0] != '0')
        strcpy(readname,indir);
    strcat(readname,"puxy.dat");

    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("PU x-y data file %s not available but required.",readname);


    while (fscanf(fp,"%i%lf%lf",&n,&x,&y)==3)
    {
        n = FastPUIDtoPUINDEX(puno,n,PULookup);
        if (n<0 || n >= puno)
            ShowErrorMessage("planning unit id out of bounds %d \n",n);
        pu[i].xloc = x;
        pu[i].yloc = y;
        i++;
    }  /* Found another valid looking line */

    return(i);
} // ReadPUXYfile

/* Read Planning Unit Data */
/* This file reads in the data relating to each planning unit. Namely, the ID, the cost,
    The status, x and y coordinates, where appropriate. */
int ReadPUData(int *puno,struct spustuff *pu[],int iCostCount,struct stringname CostNames[],struct sfname fnames)
{
    FILE *fp;

    struct pustruct{int id;double cost; int status; double xloc; double yloc;};
    char *readname;
    char sLine[1000];
    char *varlist[5] = {"id","cost","status","xloc","yloc"};
    int numvars = 5,i=0;
    char *sVarName,*sVarVal;
    struct snlink *head = NULL, *temp = NULL;
    struct spustuff putemp;
    struct spulink{struct spustuff stemp;struct spulink *next;} *spuhead = NULL, *newspulink;

    readname = (char *) calloc(strlen(fnames.puname) + strlen(fnames.inputdir)+2, sizeof(char));

    #ifdef MEMDEBUG
    iMemoryUsed += (strlen(fnames.puname) + strlen(fnames.inputdir)+2) * sizeof(char);
    ShowGenProg("memory used %i\n",iMemoryUsed);
    #endif

    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.puname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("Planning Unit file %s has not been found.\nAborting Program.",readname);
    free(readname);

    /* Scan header */
    fgets(sLine,999,fp);

    sVarName = strtok(sLine," ,;:^*\"/|\t\'\\\n");
    head = GetVarNamePU(varlist,numvars,CostNames,iCostCount,sVarName,head,fnames.puname);

    temp = head;
    while ((sVarName = strtok(NULL," ,;:^*\"/|\t\'\\\n")) != NULL)
    {
        temp->next = GetVarNamePU(varlist,numvars,CostNames,iCostCount,sVarName,head,fnames.puname);
        temp = temp->next;
    }  /* tokking out all the variable names from the header line. There are numVars of them*/

    /* While there are still lines left feed information into temporary link list */
    while (fgets(sLine,999,fp))
    {
        i++;
        putemp.id = -1; /* Set defaults for any missing values */
        putemp.cost = 1;
        putemp.status = 0;
        putemp.xloc = -1;
        putemp.yloc = -1;
        for (temp = head;temp;temp = temp->next)
        {
            if (temp == head)
                sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            else
                sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
            if (strcmp("id",temp->name)==0)
            {
                sscanf(sVarVal,"%d",&putemp.id);
            }
            else if (strcmp("status",temp->name)==0)
            {
                sscanf(sVarVal,"%d",&putemp.status);
            }
            else if (strcmp("xloc",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&putemp.xloc);
            }
            else if (strcmp("yloc",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&putemp.yloc);
            }
            else if (strcmp("cost",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&putemp.cost);
            }
        } /* looking for ivar different input variables */

        if (putemp.id == -1)
            ShowErrorMessage("ERROR: Missing planning unit id for line %d. \n",i);
        /* Stick everything from putemp into link list */
        newspulink = (struct spulink *) malloc(sizeof(struct spulink));
        newspulink->stemp.id = putemp.id;
        newspulink->stemp.status = putemp.status;
        newspulink->stemp.cost = putemp.cost;
        newspulink->stemp.xloc = putemp.xloc;
        newspulink->stemp.yloc = putemp.yloc;
        newspulink->next = spuhead;
        spuhead = newspulink;

    } /* while still lines in data file */

    fclose(fp);

    /* Create array to store the information */
    *puno = i;
    *pu = (struct spustuff *) calloc(*puno,sizeof(struct spustuff));

    #ifdef MEMDEBUG
    iMemoryUsed += (*puno) * sizeof(struct spustuff);
    ShowGenProg("memory used %i\n",iMemoryUsed);
    #endif

    for (i=0;i<*puno;i++)
    {
        (* pu)[i].id = spuhead->stemp.id;
        (* pu)[i].cost = spuhead->stemp.cost;
        (* pu)[i].status = spuhead->stemp.status;
        (* pu)[i].xloc = spuhead->stemp.xloc;
        (* pu)[i].yloc = spuhead->stemp.yloc;
        (* pu)[i].richness = 0;
        (* pu)[i].offset = 0;
        (* pu)[i].fPULock = 0;
        (* pu)[i].iPULock = 0;
        (* pu)[i].fPUZone = 0;
        (* pu)[i].iPUZone = 0;
        (* pu)[i].iPUZones = 0;
        (* pu)[i].iPreviousStatus = 0;
        newspulink = spuhead;
        spuhead = spuhead->next;
        free(newspulink);
    }
    return(i);
} // ReadPUData

/****** Read Species Information Data  ****/
int ReadSpeciesData2(int *spno,struct sspecies *spec[],struct sfname fnames)
{
    FILE *fp;
    int n=0;
    struct snlink *snhead= NULL,*temp;
    struct slink{struct sspecies stemp;struct slink *next;} *head = NULL,*newlink;
    struct sspecies spectemp;
    char *readname;
    char speciesname[255];
    char sLine[1000];
    char *varlist[10] = {"id","type","target","spf","target2","sepdistance","sepnum","name","targetocc","prop"};
    int numvars = 10,namespecial = 0;
    char *sVarName,*sVarVal;

    readname = (char *) calloc(strlen(fnames.specname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.specname);
    fp = fopen(readname,"r");
    if (fp==NULL)
       ShowErrorMessage("Species file %s has not been found.\nAborting Program.",readname);
    free(readname);

    /* Scan header */
    fgets(sLine,999,fp);

    sVarName = strtok(sLine," ,;:^*\"/|\t\'\\\n");
    snhead = GetVarName(varlist,numvars,sVarName,snhead,fnames.specname);
    temp = snhead;
    while ((sVarName = strtok(NULL," ,;:^*\"/|\t\'\\\n")) != NULL)
    {
        temp->next = GetVarName(varlist,numvars,sVarName,snhead,fnames.specname);
        temp = temp->next;
    }  /* tokking out all the variable names from the header line. There are numVars of them*/

    /* While there are still lines left feed information into temporary link list */
    while (fgets(sLine,999,fp))
    {
        n++;
        /** Clear important species stats **/
        spectemp.name = -1;
        spectemp.target = -1;
        spectemp.type = -1;
        spectemp.spf = -1;
        spectemp.target2 = -1;
        spectemp.targetocc = -1;
        spectemp.sepdistance = -1;
        spectemp.sepnum = -1;
        spectemp.prop = -1;
        speciesname[0] = '\0';
        for (temp = snhead;temp;temp = temp->next)
        {  /* Tok out the next Variable */

            if (namespecial)
            {  /* used for special name handling function */
                namespecial = 0;
            }
            else {
                if (temp == snhead)
                {
                       sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
                }
                else
                    sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
                }
            if (strcmp("id",temp->name)==0)
            {
                sscanf(sVarVal,"%d",&spectemp.name);
            }
            else if (strcmp("name",temp->name)==0)
            {
                /* Cpy first part of name into this */
                strcpy(speciesname,sVarVal);
                /* get next part of name */
                do
                {
                    sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
                    if (!sVarVal)
                       namespecial = 2;
                    else
                    {
                        if(isalpha(sVarVal[0]))
                        {
                            strcat(speciesname," ");
                            strcat(speciesname,sVarVal);
                        }
                        else
                            namespecial = 1;
                    }
                } while (!namespecial);
                if (namespecial == 2)
                   namespecial = 0; /* Handles end of line case */
                /* namespecial == 1 means not at end of line and next variable should be processed*/
            } /* Special name handling routine */
            else if (strcmp("type",temp->name)==0)
            {
                sscanf(sVarVal,"%d",&spectemp.type);
            }
            else if (strcmp("target",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&spectemp.target);
            }
            else if (strcmp("prop",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&spectemp.prop);
                if (spectemp.prop > 0)
                   fSpecPROPLoaded = 1;
            }
            else if (strcmp("spf",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&spectemp.spf);
            }
            else if (strcmp("sepdistance",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&spectemp.sepdistance);
            }
            else if (strcmp("sepnum",temp->name)==0)
            {
                sscanf(sVarVal,"%d",&spectemp.sepnum);
            }
            else if (strcmp("target2",temp->name)==0)
            {
                sscanf(sVarVal,"%lf",&spectemp.target2);
            }
            else if (strcmp("targetocc",temp->name)==0)
            {
                sscanf(sVarVal,"%d",&spectemp.targetocc);
            }
        } /* looking for ivar different input variables */
        newlink = (struct slink *) malloc(sizeof(struct slink));
        newlink->stemp.name = spectemp.name;
        newlink->stemp.target = spectemp.target;
        newlink->stemp.prop = spectemp.prop;
        newlink->stemp.spf = spectemp.spf;
        newlink->stemp.type = spectemp.type;
        newlink->stemp.targetocc = spectemp.targetocc;
        newlink->stemp.target2 = spectemp.target2;
        newlink->stemp.sepdistance = spectemp.sepdistance;
        newlink->stemp.sepnum = spectemp.sepnum;
        newlink->stemp.sname = (char *) calloc(strlen(speciesname)+1,sizeof(char));
        strcpy(newlink->stemp.sname,speciesname);
        newlink->next = head;
        head = newlink;
    } /* Scanning through each line of file */

    fclose(fp);

    /* Now do as Name.dat in forming species array */
    *spno = n;
    *spec = (struct sspecies *) calloc(*spno,sizeof(struct sspecies));
    /* put each link into namelist and free it */
    n = 0;
    while (head)
    {
        (* spec)[n].name = head->stemp.name;
        (* spec)[n].type = head->stemp.type;
        (* spec)[n].target = head->stemp.target;
        (* spec)[n].prop = head->stemp.prop;
        (* spec)[n].spf = head->stemp.spf;
        (* spec)[n].target2 = head->stemp.target2;
        (* spec)[n].targetocc = head->stemp.targetocc;
        (* spec)[n].sepdistance = head->stemp.sepdistance;
        (* spec)[n].sepnum = head->stemp.sepnum;
        (* spec)[n].richness = 0;
        (* spec)[n].sname = (char *) calloc(strlen(head->stemp.sname)+1,sizeof(char));
        strcpy((* spec)[n].sname,head->stemp.sname);
        n++;
        newlink = head;
        head = head->next;
        free(newlink->stemp.sname);
        free(newlink);
    }

    return(n);
} // ReadSpeciesData2


/* Read General Species Data */
/* This function reads in a fixed file named file the general species info. It requires that
    species are typed and changes the value of various characteristics of species
    of that type. This is done in a separate function */

int ReadGenSpeciesData(int *gspno,struct sgenspec *gspec[],struct sfname fnames)
{
    FILE *fp;
    char *readname;
    char sLine[1000];
    char *varlist[8] = {"type","target","target2","targetocc","sepnum","sepdistance","prop","spf"};
    int numvars = 8,i=0;
    char *sVarName,*sVarVal;
    struct snlink *head = NULL, *temp = NULL;
    struct sgenspec gstemp;
    struct sgslink{struct sgenspec stemp;struct sgslink *next;} *gshead = NULL, *newgslink;


    /* Find and Open File */
    readname = (char *) calloc(strlen(fnames.blockdefname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.blockdefname);
    fp = fopen(readname,"r");
    if (fp==NULL)
    {
        ShowWarningMessage("Warning: Block Definition File %s not found.\n",fnames.blockdefname);
        free(readname);
        return(0);
    }
    free(readname);

    /* Scan header */
    fgets(sLine,999,fp);
    sVarName = strtok(sLine," ,;:^*\"/|\t\'\\\n");
    head = GetVarName(varlist,numvars,sVarName,head,fnames.blockdefname);
    temp = head;
    while ((sVarName = strtok(NULL," ,;:^*\"/|\t\'\\\n")) != NULL)
    {
        temp->next = GetVarName(varlist,numvars,sVarName,head,fnames.blockdefname);
        temp = temp->next;
    }  /* tokking out all the variable names from the header line. There are numVars of them*/
    /* While there are still lines left feed information into temporary link list */
    while (fgets(sLine,999,fp))
    {
        i++;
        gstemp.type = -1; /* Set defaults for any missing values */
        gstemp.targetocc = -1;
        gstemp.target = -1;
        gstemp.target2 = -1;
        gstemp.sepnum = -1;
        gstemp.sepdistance = -1;
        gstemp.prop = -1;
        for (temp = head;temp->next;temp = temp->next)
        {
            if (temp == head)
                sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
            else
                sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");

            if (strcmp("type",temp->name)==0)
            {
                 sscanf(sVarVal,"%d",&gstemp.type);
            } else {
                if (strcmp("targetocc",temp->name)==0)
                {
                    sscanf(sVarVal,"%d",&gstemp.targetocc);
                } else {
                    if (strcmp("target",temp->name)==0)
                    {
                        sscanf(sVarVal,"%lf",&gstemp.target);
                    } else {
                        if (strcmp("target2",temp->name)==0)
                        {
                            sscanf(sVarVal,"%lf",&gstemp.target2);
                        } else {
                            if (strcmp("sepnum",temp->name)==0)
                            {
                                sscanf(sVarVal,"%d",&gstemp.sepnum);
                            } else {
                                if (strcmp("sepdistance",temp->name)==0)
                                {
                                    sscanf(sVarVal,"%lf",&gstemp.sepdistance);
                                } else {
                                    if (strcmp("prop",temp->name)==0)
                                    {
                                        sscanf(sVarVal,"%lf",&gstemp.prop);
                                    } else {
                                        if (strcmp("spf",temp->name)==0)
                                        {
                                            sscanf(sVarVal,"%lf",&gstemp.spf);
                                        } else {
                                            ShowWarningMessage("Cannot find |%s| \n",temp->name);
                                            ShowErrorMessage("Serious error in GenSpecies data reading function.\n");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } /* looking for ivar different input variables */

        if (gstemp.type== -1)
            ShowErrorMessage("ERROR: Missing Gen Species type for line %d. \n",i);
        /* Stick everything from gstemp into link list */
        newgslink = (struct sgslink *) malloc(sizeof(struct sgslink));
        newgslink->stemp.type = gstemp.type;
        newgslink->stemp.targetocc = gstemp.targetocc;
        newgslink->stemp.target = gstemp.target;
        newgslink->stemp.target2 = gstemp.target2;
        newgslink->stemp.sepnum = gstemp.sepnum;
        newgslink->stemp.sepdistance = gstemp.sepdistance;
        newgslink->stemp.prop = gstemp.prop;
        newgslink->next = gshead;
        gshead = newgslink;
    } /* while still lines in data file */

    fclose(fp);

    /* Now do as Name.dat in forming species array */
    *gspno = i;
    *gspec = (struct sgenspec *) calloc(*gspno,sizeof(struct sgenspec));
    /* put each link into namelist and free it */
    i = 0;
    while (gshead)
    {
        (* gspec)[i].type = gshead->stemp.type;
        (* gspec)[i].targetocc = gshead->stemp.targetocc;
        (* gspec)[i].target = gshead->stemp.target;
        (* gspec)[i].target2 = gshead->stemp.target2;
        (* gspec)[i].sepnum = gshead->stemp.sepnum;
        (* gspec)[i].sepdistance = gshead->stemp.sepdistance;
        (* gspec)[i].prop = gshead->stemp.prop;
        i++;
        newgslink = gshead;
        gshead = gshead->next;
        free(newgslink);
    }

    return(i);
} // ReadGenSpeciesData

int ReadConnections(int puno,struct sconnections connections[],int verbose,
                    struct spustuff pu[],struct binsearch PULookup[],struct sfname fnames)
{
    FILE *fp;
    int id1,id2;
    double fcost;
    int icount = 0,idup = 0;
    int bad;
    struct sneighbour *p;
    char *readname;
    int numvars = 3,ivars;
    char *sVarName,*sVarVal;
    char sLine[1000];

    readname = (char *) calloc(strlen(fnames.connectionname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.connectionname);
    fp = fopen(readname,"r");
    if (fp==NULL)
    {
        ShowGenProg("Warning: Connection File %s not found ",fnames.connectionname);
        free(readname);
        return(0);
    }
    free(readname);

    fgets(sLine,999,fp); /* Scan header line */

    while (fgets(sLine,999,fp))
    {
        icount++;

        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&id1);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&id2);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%lf",&fcost);

        id1 = FastPUIDtoPUINDEX(puno,id1,PULookup);
        id2 = FastPUIDtoPUINDEX(puno,id2,PULookup);

        if (id1==id2)
        {
            #ifdef ASYMCON
            if (asymmetricconnectivity)
                connections[id1].fixedcost = 0;
            else
            #endif
                connections[id1].fixedcost += fcost;

            continue;
        } /* My code for an irremovable connection */
        if (id1>=0 && id1<puno)
        {  /* Is I a sensible number ?*/
            p = connections[id1].first;
            bad = 0;
            while (p)
            {
                if (p->nbr == id2)
                    bad = 1;

                p = p->next;
            }

            #ifdef ASYMCON
            if (asymmetricconnectivity)
                bad = 0;
            #endif

            if (bad)
                ShowDetProg("Double connection definition %i %i \n",pu[id1].id,pu[id2].id);
            else
            {
                connections[id1].nbrno++;
                p = (struct sneighbour *) malloc (sizeof(struct sneighbour));

                #ifdef MEMDEBUG
                iMemoryUsed += sizeof(struct sneighbour);
                #endif

                p->cost = fcost;
                p->nbr = id2;
                p->next = connections[id1].first;

                #ifdef ASYMCON
                if (asymmetricconnectivity)
                {
                    p->connectionorigon = 1;
                }
                else
                    p->connectionorigon = 1;
                #endif

                connections[id1].first = p;
            }
        }
        else
            ShowErrorMessage("A connection is out of range %f %i %i \n",fcost,id1,id2);

        if (id2>=0 && id2<puno)
        {  /* Is I a sensible number ?*/
            p = connections[id2].first;
            bad = 0;
            while (p)
            {
                if (p->nbr == id1)
                    bad = 1;
                p = p->next;
            }

            #ifdef ASYMCON
            if (asymmetricconnectivity)
                bad = 0;
            #endif

            if (bad && verbose > 4)
                ShowDetProg("Double connection definition %i %i \n",id1,id2);
            if (bad)
                idup++;
            else
            {
                connections[id2].nbrno++;
                p = (struct sneighbour *) malloc (sizeof(struct sneighbour));
                #ifdef MEMDEBUG
                iMemoryUsed += sizeof(struct sneighbour);
                #endif
                p->cost = fcost;
                p->nbr = id1;
                p->next = connections[id2].first;

                #ifdef ASYMCON
                p->connectionorigon = 1;
                if (asymmetricconnectivity)
                    p->connectionorigon = 0;
                else
                    p->connectionorigon = 1;
                #endif

                connections[id2].first = p;
            }
        }
        else
            ShowErrorMessage("A connection is out of range %f %i %i \n",fcost,id1,id2);
    }

    fclose(fp);

    if (idup)
        ShowGenProg("There were %i duplicate connection definitions.\n",idup);

    return(icount);

} /*** Read Connections ***/


/* * * *  Reading in the Planning Unit versus Species Data Type 2. Relational data *****/
/* * * *  into a Sparse Matrix data structure *****/
void LoadSparseMatrix(int *iSMSize, struct spu *SM[], int puno, int spno, struct spustuff PU[],
                      struct binsearch PULookup[],struct binsearch SPLookup[],
                      struct sfname fnames)
{
    FILE *fp;
    char *readname;
    char sLine[1000];
    char sValVal;
    int ivars, iLastPUID;
    char *sVarName,*sVarVal;
    int i, _spid, _puid, iInternalSMSize = 0, iBigMatrixSize;
    double amount, rDensity, rInternalSMSize, rBigMatrixSize;
    #ifdef DEBUGTRACEFILE
    char debugbuffer[300];
    #endif

    readname = (char *) calloc(strlen(fnames.puvsprname) + strlen(fnames.inputdir)+2, sizeof(char));

    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.puvsprname);
    fp = fopen(readname,"r");
    if (fp==NULL)
        ShowErrorMessage("PU v Species file %s not found\nAborting Program.",readname);
    free(readname);

    /* read through the file first to see how many lines */
    fgets(sLine,999,fp);
    while (fgets(sLine,999,fp))
    {
        iInternalSMSize++;
    }

    rewind(fp);
    fgets(sLine,999,fp);

    *iSMSize = iInternalSMSize;

    /* create the sparse matrix */
    *SM = (struct spu *) calloc(iInternalSMSize,sizeof(struct spu));

    /* planning unit richness and offset are already set to zero */

    iLastPUID = -1;

    /* init with zero values */
    for (i=0;i<iInternalSMSize;i++)
    {
        fgets(sLine,999,fp);

        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&_spid);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&_puid);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%lf",&amount);

        if (_puid < iLastPUID)
        {
            // error condition exists, file is not in ascending order for PUID

            #ifdef DEBUGTRACEFILE
            sprintf(debugbuffer,"Error: PU v Species file %s is not in ascending order for PUID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
            AppendDebugTraceFile(debugbuffer);
            #endif

            ShowErrorMessage("Error: PU v Species file %s is not in ascending order for PUID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
        }

        iLastPUID = _puid;

        _puid = FastPUIDtoPUINDEX(puno,_puid,PULookup);
        _spid = FastSPIDtoSPINDEX(spno,_spid,SPLookup);

        /* increment richness for planning unit containing this feature */
        PU[_puid].richness += 1;
        /* if planning units richness is one, set its offset */
        if (PU[_puid].richness == 1)
            PU[_puid].offset = i;

        (* SM)[i].amount = amount;
        (* SM)[i].clump = 0;
        (* SM)[i].spindex = _spid;
    }

    fclose(fp);


    iBigMatrixSize = puno * spno;
    rInternalSMSize = iInternalSMSize;
    rBigMatrixSize = iBigMatrixSize;
    rDensity = rInternalSMSize / rBigMatrixSize * 100;

    ShowGenProg("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                iInternalSMSize,iBigMatrixSize,rDensity);
} // LoadSparseMatrix

void LoadSparseMatrix_sporder(int *iSMSize, struct spusporder *SM[], int puno, int spno,
                              struct binsearch PULookup[],struct binsearch SPLookup[],// typesp spec[],
                              struct sfname fnames)
{
    FILE *fp;
    char *readname,sLine[500],*sVarName,*sVarVal;
    int i, _spid,spid, _puid, iInternalSMSize = 0, iBigMatrixSize, iLastSPID;
    double amount, rDensity, rInternalSMSize, rBigMatrixSize;
    #ifdef DEBUGTRACEFILE
    char debugbuffer[300];
    #endif

    readname = (char *) calloc(strlen(fnames.matrixspordername) + strlen(fnames.inputdir)+2, sizeof(char));

    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.matrixspordername);
    if((fp = fopen(readname,"r"))==NULL)
        ShowErrorMessage("PU v Species file %s not found\nAborting Program.",readname);
    free(readname);

    /* read through the file first to see how many lines */
    fgets(sLine,500-1,fp);
    while (fgets(sLine,500-1,fp))
        iInternalSMSize++;

    rewind(fp);
    fgets(sLine,500-1,fp);

    *iSMSize = iInternalSMSize;

    /* create the sparse matrix */
    *SM = (struct spusporder *) calloc(iInternalSMSize,sizeof(struct spusporder));

    iLastSPID = -1;
    /* planning unit richness and offset are already set to zero */

    /* init with zero values */
    for (i=0;i<iInternalSMSize;i++)
    {
        fgets(sLine,500-1,fp);

        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&_spid);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&_puid);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%lf",&amount);

        if (_spid < iLastSPID)
        {
            // error condition exists, file is not in ascending order for SPID

            #ifdef DEBUGTRACEFILE
            sprintf(debugbuffer,"Error: PU v Species file %s is not in ascending order for SPID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
            AppendDebugTraceFile(debugbuffer);
            #endif

            ShowErrorMessage("Error: PU v Species file %s is not in ascending order for SPID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
        }

        iLastSPID = _spid;

        _puid = FastPUIDtoPUINDEX(puno,_puid,PULookup);
        spid = FastSPIDtoSPINDEX(spno,_spid,SPLookup);

        /* increment richness for planning unit containing this feature */
        spec[spid].richness += 1;
        /* if planning units richness is one, set its offset */
        if (spec[spid].richness == 1)
            spec[spid].offset = i;

        (* SM)[i].amount = amount;
        (* SM)[i].puindex = _puid;
    }

    fclose(fp);

    iBigMatrixSize = puno * spno;
    rInternalSMSize = iInternalSMSize;
    rBigMatrixSize = iBigMatrixSize;
    rDensity = rInternalSMSize / rBigMatrixSize * 100;

    ShowGenProg("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                iInternalSMSize,iBigMatrixSize,rDensity);
} // LoadSparseMatrix_sporder

void LoadPenalty(int spno,struct sspecies spec[],struct sfname fnames)
{
    FILE *fp;
    int i, _spid;
    char *readname,sLine[500],*sVarName,*sVarVal;

    readname = (char *) calloc(strlen(fnames.penaltyname) + strlen(fnames.inputdir)+2, sizeof(char));

    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.penaltyname);
    if((fp = fopen(readname,"r"))==NULL)
        ShowErrorMessage("penalty file %s not found\nAborting Program.",readname);
    free(readname);

    // read header row
    fgets(sLine,500-1,fp);

    i=0;
    // read data rows
    while (fgets(sLine,500-1,fp))
    {
        sVarVal = strtok(sLine," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%d",&_spid);
        sVarVal = strtok(NULL," ,;:^*\"/|\t\'\\\n");
        sscanf(sVarVal,"%lf",&spec[i].penalty);

        i++;
    }

    fclose(fp);
} // LoadPenalty

