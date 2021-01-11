
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

