/* 
Portions Copyright (c) 2007-2012 Hogue Laboratory, National University of Singapore
Portions Copyright (c) 1997-2007 Hogue Laboratory, University of Toronto
Portions Copyright (c) 1997-2005 Hogue Laboratory, Mount Sinai Hospital, Toronto


All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE COPYRIGHT HOLDERS HAVE NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS OR MODIFICATIONS.

Principal Contact: Christopher W.V. Hogue  chogue@blueprint.org

Contributing Hogue Laboratory Developers and Institutional Collaborators: 
(SLRI and TraDES Directories)
Christopher Hogue (mmdbapi, b-d tree, Cn3D, trades, solvate, valmerge, str2trj, seq2trj). 
Howard Feldman (salign, base foldtraj & maketraj libraries)
Contributions from John J. Salama (VisTraj) and Kevin Snyder (valmerge)
Phillipe Phan (analyzeMovie) Adrian Heilbut, Mark Kotowycz, Van Le (trj visualization),
Michel Dumontier (scoring function libraries, statistical functions) 
Elena Garderman (software maintenance), Mingxi Yao (valtopdb),
Gil Alteroviz (GOR implementation),  Boris Steipe (University of Toronto) (fragment based construction)
Brendan McConkey (University of Waterloo) & Michael Brougham(VSCORE scoring function)
Jonathan Kans (NCBI) (vistraj)

*/




/* includes */
#include "benchtraj.h"

#ifdef OS_UNIX
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif
#ifdef OS_UNIX_QNX
#include <sched.h>
#endif



#define BENCH_TRAJ "bench"

#include <d4all.h>  /* required to match 32 or 64 bit codebase headers and libraries properly */

static Int4 totalstruc=0,experiment=0;
static Int2 generation=0;
static Int4 strucs_per_run=5000;
static Char fnamtrj[PATH_MAX];
Boolean volatile timetoquit=FALSE;
Char handle[HANDLE_LENGTH+1];
#define STRUCS_PER_RUN strucs_per_run
#define ZHANG_WINDOWSIZE 4
#define BRYANT_WINDOWSIZE 4

/* Global Variables */



Boolean CheckForAbort(void)
{

	return FALSE;
}



#ifdef OS_MSWIN
TrajErr GetMSWinVerInfo(CharPtr buf)
{
	OSVERSIONINFOEX osvi;
	BOOL bOsVersionInfoEx;
	Char tmpbuf[PATH_MAX];
	HKEY hKey;
	char szRegbuf[255];
	DWORD dwBufLen;

	StringCpy(buf,"");
	RegOpenKeyEx(HKEY_LOCAL_MACHINE,"SYSTEM\\CurrentControlSet\\Services\\Tcpip\\Parameters",0,KEY_QUERY_VALUE,&hKey);
	dwBufLen=255;
	RegQueryValueEx(hKey,"Hostname",NULL,NULL,(LPBYTE)szRegbuf,&dwBufLen);
	StringCat(buf,szRegbuf);
	dwBufLen=255;
	RegQueryValueEx(hKey,"Domain",NULL,NULL,(LPBYTE)szRegbuf,&dwBufLen);
	if (dwBufLen>0)
		StringCat(buf,".");
	StringCat(buf,szRegbuf);
	if (dwBufLen>0)
		StringCat(buf," ");
	RegCloseKey(hKey);
	ZeroMemory(&osvi,sizeof(OSVERSIONINFOEX));
	osvi.dwOSVersionInfoSize=sizeof(OSVERSIONINFOEX);
	if (!(bOsVersionInfoEx=GetVersionEx((OSVERSIONINFO *) &osvi))) {
		osvi.dwOSVersionInfoSize=sizeof(OSVERSIONINFO);
		if (!GetVersionEx((OSVERSIONINFO *) &osvi))
			return ERR_FAIL;
	}
	switch (osvi.dwPlatformId) {
		case VER_PLATFORM_WIN32_NT:
			if (osvi.dwMajorVersion<=4)
				StringCat(buf,"Microsoft Windows NT ");
			if (osvi.dwMajorVersion==5)
				StringCat(buf,"Microsoft Windows 2000 ");
			/*if (bOsVersionInfoEx) {
				if (osvi.wProductType==VER_NT_WORKSTATION)
					StringCat(buf,"Professional ");
				if (osvi.wProductType==VER_NT_SERVER)
					StringCat(buf,"Server ");
			}
			else*/						
			dwBufLen=255;
			RegOpenKeyEx(HKEY_LOCAL_MACHINE,"SYSTEM\\CurrentControlSet\\Control\\ProductOptions",0,KEY_QUERY_VALUE,&hKey);
			RegQueryValueEx(hKey,"ProductType",NULL,NULL,(LPBYTE)szRegbuf,&dwBufLen);
			RegCloseKey(hKey);
			if (lstrcmpi("WINNT",szRegbuf)==0)
				StringCat(buf,"Workstation ");
			if (lstrcmpi("SERVERNT",szRegbuf)==0)
				StringCat(buf,"Server ");
			sprintf(tmpbuf,"version %d.%d %s (Build %d)",osvi.dwMajorVersion,osvi.dwMinorVersion,osvi.szCSDVersion,osvi.dwBuildNumber & 0xFFFF);
			StringCat(buf,tmpbuf);
			break;
		case VER_PLATFORM_WIN32_WINDOWS:
			if ((osvi.dwMajorVersion>4) || ((osvi.dwMajorVersion==4) && (osvi.dwMinorVersion>0)))
				StringCat(buf,"Microsoft Windows 98");
			else
				StringCat(buf,"Microsoft Windows 95");
			break;
		case VER_PLATFORM_WIN32s:
			StringCat(buf,"Microsoft Win32s");
			break;
	}
	StringCat(buf,"\n");
	return ERR_SUCCESS;					
}
#endif


void InitBSPLog(Int4 randseed,CharPtr sequence,Int4 tstruc)
{
	CharPtr decsequence;
	static Char pcUname[PATH_MAX]="";
	Char timedate[25];
#ifdef OS_UNIX
    struct utsname utsbuf;
#endif

    bspTempLog=BSNew(0);
    /* compute log file name */
	if (generation==0)
	    sprintf(LogOutName,"%sfold_%ld_%s_%ld_%s",CFG_local_datafilepath,(long)experiment,handle,(long)tstruc,fnamtrj);
	else
	    sprintf(LogOutName,"%sfold_%ld_%s_%ld_%s_%s_%d",CFG_local_datafilepath,(long)experiment,handle,(long)tstruc,handle,fnamtrj,generation);
    BSprintf(bspTempLog,"\nFoldtraj v%s log report",FOLDTRAJ_VERSION);
    BSprintf(bspTempLog,"\nTrajectory File: %s%s_%d",CFG_local_datafilepath,fnamtrj,generation);
    BSprintf(bspTempLog,"\tStructure File Base Name: %s_%s_%d",handle,fnamtrj,generation);
    BSprintf(bspTempLog,"\tRandom Seed: %ld\t# Generated: %ld\tStart Numbering at: %ld",randseed,STRUCS_PER_RUN,1);
    BSprintf(bspTempLog,"\tCompared to Native Structure: 1CDZ");

    decsequence=DecodeSequence(sequence,EXTAA_PARENT);
    BSprintf(bspTempLog,"\nSequence: %s",decsequence);
    decsequence=MemFree(decsequence);
    BSprintf(bspTempLog,"\nFolding Conditions:\t%ldx%ld Trajectory Distributions\t%s Random Walk\tAverage Timeout: %ld\tBackbone Error Tolerance: %4.2f\tBackbone Precision: %6.4f\tBackbone Atom Bounciness: %5.2f%%\tSidechain Atom Bounciness: %5.2f %%\tHydrogen Bumpchecking: %s\t# Rotamer Tries: %d/chi angle\tMarkov Scale Factor: %4.2f",
            TRAJDIV,TRAJDIV,WALKTYPE==WALK_CA?"C-Alpha":"Phi-Psi",TIMEOUT,BACKBONE_ERROR_TOLERANCE,
            BACKBONE_PRECISION,ATOM_BOUNCINESS_BB*100.0,ATOM_BOUNCINESS_SC*100.0,
            BUMPCHECK_HYDROGEN?"On":"Off",NUM_ROT_TRIES,MARKOV_SCALE_FACTOR);
    if (WALKTYPE!=WALK_CA && MARKOV_SCALE_FACTOR>0.0)
        ErrPostEx(SEV_ERROR,2,2,"Markov Scale Factor should be zero for non-Ca walk -- ignored");
    if (TRAJTYPE!=TRAJ_NA) {
        BSprintf(bspTempLog,"\tTrajectory Distribution: ");
        switch (TRAJTYPE) {
            case TRAJ_UNIFORM:
                BSprintf(bspTempLog,"Uniform");
                break;
            case TRAJ_STANDARD:
                BSprintf(bspTempLog,"Amino-Acid Based");
                break;
            case TRAJ_SSTRU:
                BSprintf(bspTempLog,"1-State Secondary Structure Prediction");
                break;
            case TRAJ_GOR:
                BSprintf(bspTempLog,"3-State Secondary Structure Prediction");
                break;
            default:
                BSprintf(bspTempLog,"Unknown");
        }
    }
    if (StringCmp(CONSTRAINT_FILE,""))
            BSprintf(bspTempLog,"\tConstraint file: %s",CONSTRAINT_FILE);
    BSprintf(bspTempLog,"\nSystem Information: ");
	
	if (StringLen(pcUname)==0) {
	    StringCpy(pcUname,"Not Available\n");
#ifdef OS_UNIX
	    /* use uname to get sysinfo */
		if (uname(&utsbuf)>=0) {
	        sprintf(pcUname,"%s %s %s %s %s\n",utsbuf.sysname,utsbuf.nodename,utsbuf.release,utsbuf.version,utsbuf.machine);
		}
#else
#ifdef OS_MSWIN
	    /* assume ver is available */
		if (GetMSWinVerInfo(pcUname)!=ERR_SUCCESS) {
            StringCpy(pcUname,"Not Available\n");
	
	    }       
#endif
#endif
	}
	BSWrite(bspTempLog,pcUname,StringLen(pcUname));
    DayTimeStr(timedate,TRUE,TRUE);
    BSprintf(bspTempLog,"Job started: %s\n",timedate);
    BSprintf(bspTempLog,"\nStructure\tTime(s)\tTries\tBad Backbone Tries\tCrashes\tDistant Constraint Tries (# Violations)\tChain Length\tRadius of Gyration(A)\tHydrophobics Radius of Gyration(A)\tEnd-to-end Distance(A)\tRn\tCn");
    BSprintf(bspTempLog,"\tSurface Accessibility (A^2)\tExposed Hydrophobics (A^2)");
    BSprintf(bspTempLog,"\t# Helical Residues (DSSP)\t# Extended Residues (DSSP)");
    BSprintf(bspTempLog,"\t# Extended Residues (CA-CA dist.)");
    BSprintf(bspTempLog,"\tRMSD from Native(A)");
    BSprintf(bspTempLog,"\tZhang Potential(Exclusive Window Size %d)\tBryant-Lawrence Potential(Exclusive Window Size %d)\tCrease Energy\n",ZHANG_WINDOWSIZE,BRYANT_WINDOWSIZE);
}


static TrajErr SetConfigFileValues(Int2 wt, FloatLo errtol, FloatLo abbb, FloatLo absc, FloatLo tunnelprob) {
	Char str[PATH_MAX];

	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","TRAJDIV","400"))
		return ERR_FAIL;
	if (wt==WALK_CA) {
		if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","WALKTYPE","CA"))
			return ERR_FAIL;
	}
	else {
		if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","WALKTYPE","PHIPSI"))
			return ERR_FAIL;
	}
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","TIMEOUT","100"))
		return ERR_FAIL;
	sprintf(str,"%1.3f",errtol);
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","BACKBONE_ERROR_TOLERANCE",str))
		return ERR_FAIL;
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","BACKBONE_PRECISION","0.050"))
		return ERR_FAIL;
	sprintf(str,"%1.3f",abbb);
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","ATOM_BOUNCINESS_BB",str))
		return ERR_FAIL;
	sprintf(str,"%1.3f",absc);
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","ATOM_BOUNCINESS_SC",str))
		return ERR_FAIL;
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","BUMPCHECK_HYDROGEN","FALSE"))
		return ERR_FAIL;
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","INCSIZE","10.0"))
		return ERR_FAIL;
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","START_BACKBONE","250.0"))
		return ERR_FAIL;
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","NUM_ROT_TRIES","7"))
		return ERR_FAIL;
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","MARKOV_SCALE_FACTOR","0.0"))
		return ERR_FAIL;
	sprintf(str,"%1.3f",tunnelprob);
	if (!SetAppParam(TASK_CFGFILE,"INITTRAJ","TUNNEL_PROB",str))
		return ERR_FAIL;
	return ERR_SUCCESS;
}


static TrajErr TRADESBenchmark(void)
{
	/* protein is 1CDZ, same as beta protein - mixed alpha/beta protein 96 AA */
	static Char sequence[]="ELPDFFQGKHFFLYGEFPGDERRKLIRYVTAFNGELEDYMSDRVQFVITAQEWDPSFEEALMDNPSLAFVRPRWIYSCNEKQKLLPHQLYGVVPQA";
	CharPtr pc;
	FILE *f;
	Char fnamseq[PATH_MAX];
	Char fnam[PATH_MAX];
	Char fnamtrj[PATH_MAX];
        Char fnuturndb[PATH_MAX];
	Int2 err,numAA,cnt;
	Int4 randseed,dbsize;
	pMakeTrjParamBlock mkdata;
	pFoldTrajParamBlock foldtrajprm; 
	CPUTimePtr pTime1,pTime2,pTime3,pTime4;
	BiostrucPtr bspBiostruc;
	PDNMS pdnmsModelstruc;
	PMSD pmsdHead;
	PDNML pdnmlModel;
	PMLD pmldThis;
	BiostrucIdPtr           pbi=NULL;
	BiostrucDescrPtr        pbd=NULL;


	printf("%d is the CODEBASE version here.\n6502:UNIX x86_32, 6503003:UNIX x86_64, 6500:Win32, 6401:UNIX PPC_32\n\n", S4VERSION);
        sprintf(fnuturndb,"UTURNDB"); /* to clean up errant U-Turn files CWVH 2010 */

	if (!OPENMMDBAPI(0,NULL)) {
	    ErrPostEx(SEV_ERROR,2,1,"Unable to open MMDBAPI");
	    return ERR_FAIL;
	}
	printf("One moment, opening rotamer library...\n");
	if (LoadRotLib()!=ERR_SUCCESS) {
		ErrPostEx(SEV_ERROR,1,8,"Cannot open rotamer library (%s%s) - if the file is missing, please re-install the software",CFG_local_datafilepath,NEWSCWRL_FNAM);
		return ERR_FAIL;
	}
	pc=DFPTmpNam(FALSE,FALSE);
	if ((f=FileOpen(pc,"w"))==NULL) {
		CLOSEMMDBAPI();
		FreeRotLib();
		ErrPostEx(SEV_ERROR,errno,1,"Failed to open %s file",pc);
		return ERR_FAIL;
	}
	fprintf(f,"%s\n",sequence);
	/* Make sure disk space was sufficient to write to file */
    if (fflush(f)==EOF) {
        FileClose(f);
		CLOSEMMDBAPI();
		FreeRotLib();
		ErrPostEx(SEV_ERROR,errno,1,"Failed to write to %s file",pc);
		return ERR_FAIL;
	}
	FileClose(f);
	StringCpy(fnamseq,pc);
	sprintf(fnam,"%s%s",CFG_local_datafilepath,BENCH_TRAJ);
	sprintf(fnamtrj,"%s%s%s",CFG_local_datafilepath,BENCH_TRAJ,ASN_EXT);

	mkdata = (pMakeTrjParamBlock)MemNew(sizeof(MakeTrjParamBlock));	
	mkdata->TrajMethod = 2; /* from sequence */
	mkdata->valfnamin = NULL;
	mkdata->pcChain = NULL;
	mkdata->modelnum = 1;
	mkdata->startres = 1;
	mkdata->endres = 0;
	mkdata->peakheight = 100;
	mkdata->noise = 0;
	mkdata->savechi = 0;
	mkdata->sigma_x = 0.0;
	mkdata->sigma_y = 0.0;
	mkdata->temperature = 0;
	mkdata->timestep = 0.0;
	mkdata->seqfnamin = fnamseq;
	mkdata->tgtype = 4;
	mkdata->sstrufname = NULL;
	mkdata->DoRPS = FALSE;
	mkdata->uturn = FALSE;
	mkdata->SStruinput = FALSE;
	mkdata->trjfnamout = fnam;
	mkdata->comprtype = USE_RLE;
	mkdata->constrfnamin = NULL;
	mkdata->units=UNITS_ARBITRARY;
	mkdata->templat=NULL;
	mkdata->zhangwindsize = 0;
	mkdata->alignfile=NULL;
	mkdata->ssmask=NULL;
	mkdata->dumpcsv=0;
	mkdata->benchmark=1;
	mkdata->all_coil = (Boolean)(FALSE);
	mkdata->all_beta = (Boolean)(FALSE);
	/* ensure config file is set right, even if user may have fudged it */
	if (SetConfigFileValues(WALK_PHIPSI,50.0/*errtol*/, 0.25/*abbb*/, 0.50/*absc*/, 0.0/*tunnelprob*/)!=ERR_SUCCESS) {
		CLOSEMMDBAPI();
		FreeRotLib();
		mkdata=MemFree(mkdata);
		FileRemove(pc);
#ifdef OS_UNIX
		ErrPostEx(SEV_ERROR,0,0,"Unable to write config file .foldtrajrc, cannot continue");
#else
		ErrPostEx(SEV_FATAL,0,0,"Unable to write config file foldtraj.ini, cannot continue");
#endif
		return ERR_FAIL;
	}
	printf("Predicting secondary structure and generating trajectory distribution...\n");
	pTime1=CPUTimeMeasure();
	MakeTrj((VoidPtr)mkdata);
        printf("MadeTrj complete, unpacking ASN.1 Trajectory Graph\n");
	pTime2=CPUTimeMeasure();
	mkdata=MemFree(mkdata);
	FileRemove(fnamseq);

	if (UnPackAsnTrajGraph(fnam,&numAA,sequence,NULL,&pbi,&pbd,NULL)==NULL) {
		CLOSEMMDBAPI();
		FreeRotLib();
		FileRemove(fnamtrj);
		ErrPostEx(SEV_ERROR,2,3,"Unable to read trajectory distribution %s, please create a new one",fnam);
		return ERR_FAIL;
	}
	randseed=54374;
	RandomSeed(randseed);
	StringCpy(tmpskelfname,DFPTmpNam(FALSE,FALSE));
	/* make the ASN.1 file for this protein's chemical graph */
	BuildSkelASN(sequence,tmpskelfname);
	InitBSPLog(randseed,sequence,totalstruc);
	TGInit(tmpdbasename,DB_READ,&dbsize);
	printf("Folding protein...\n");
	pTime3=CPUTimeMeasure();
	for (cnt=0;cnt<100;cnt++) {
		bspBiostruc=NULL;
		bspBiostruc=MIMEBiostrucAsnGet(tmpskelfname,"r",NULL);
		if (bspBiostruc==NULL) {
			CLOSEMMDBAPI();
			FreeRotLib();
			FileRemove(fnamtrj);
			BSFree(bspTempLog);
			FileRemove(tmpskelfname);
			TGClose();
			CleanUpDB(tmpdbasename);
                        CleanUpDB(fnuturndb);        
			ErrPostEx(SEV_ERROR,3,1,"Unable to fetch Biostruc");
			return ERR_FAIL;
		}
		if (pbi!=NULL) {
			if (bspBiostruc->id!=NULL)
				AsnGenericChoiceSeqOfFree(bspBiostruc->id,(AsnOptFreeFunc)BiostrucIdFree);
			bspBiostruc->id=pbi;
			pbi=NULL;
		}
		if (pbd!=NULL) {
			if (bspBiostruc->descr!=NULL)
				AsnGenericChoiceSeqOfFree(bspBiostruc->descr,(AsnOptFreeFunc)BiostrucDescrFree);
			bspBiostruc->descr=pbd;
			pbd=NULL;
		}
		/* now the chemical graph is correctly in memory as a modelstruc, with
			co-ordinates (and hence PALDs) assigned only to a-carbons */
		pdnmsModelstruc=MakeAModelstruc(bspBiostruc);
		if (pdnmsModelstruc==NULL) {
			CLOSEMMDBAPI();
			FreeRotLib();
			FileRemove(fnamtrj);
			BSFree(bspTempLog);
			FileRemove(tmpskelfname);
			TGClose();
			CleanUpDB(tmpdbasename);
                        CleanUpDB(fnuturndb);        
			ErrPostEx(SEV_ERROR,4,1,"Unable to convert Biostruc to Modelstruc. Please ensure your protein only contains the 20 standard amino acids.");
			return ERR_FAIL;
		}
		pmsdHead=(PMSD)(pdnmsModelstruc->data.ptrvalue);

		/* remove ppAsnOrder for all models so it will be rebuilt
			in WriteAsnAllModel */
		pdnmlModel=pmsdHead->pdnmlModels;
		while (pdnmlModel) {
			pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
			if (pmldThis->ppAsnOrder) {
				PTRVectorFree(pmldThis->ppAsnOrder,0);
				pmldThis->ppAsnOrder=NULL;
			}
			pdnmlModel=pdnmlModel->next;
		}
		if (dbsize!=((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue))->iResCount) {
			CLOSEMMDBAPI();
			FreeRotLib();
			FreeAModelstruc(pdnmsModelstruc);
			BSFree(bspTempLog);
			FileRemove(fnamtrj);
			FileRemove(tmpskelfname);
			TGClose();
			CleanUpDB(tmpdbasename);
                       	CleanUpDB(fnuturndb);        
			ErrPostEx(SEV_ERROR,1,1,"protein length inconsistency error, expect: %d actual: %d, possibly due to corrupt trajectory file; aborting",((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue))->iResCount,dbsize);
			return ERR_FAIL;
		}
		foldtrajprm=(pFoldTrajParamBlock)MemNew(sizeof(FoldTrajParamBlock));
		foldtrajprm->pmsdRoot=pmsdHead;
		foldtrajprm->Model=1;
		foldtrajprm->err=0;
		foldtrajprm->gen=0;
		foldtrajprm->errorfile=NULL;
		/* do the folding here */
		TRADEProtein((VoidPtr)foldtrajprm,1);
		err=foldtrajprm->err;
		foldtrajprm=MemFree(foldtrajprm);
		FreeAModelstruc(pdnmsModelstruc);
        if (err!=ERR_SUCCESS && err!=ERR_INCOMPLETE) {
			FileRemove(fnamtrj);
			FileRemove(tmpskelfname);
			BSFree(bspTempLog);
			TGClose();
			CleanUpDB(tmpdbasename);
  	                CleanUpDB(fnuturndb);        
			CLOSEMMDBAPI();
			FreeRotLib();
			ErrPostEx(SEV_ERROR,1,1,"Benchtraj returned a folding error %d",err);
			return ERR_FAIL;
		}
	}
	pTime4=CPUTimeMeasure();
	FileRemove(fnamtrj);
	FileRemove(tmpskelfname);
	BSFree(bspTempLog);
	TGClose();
	CleanUpDB(tmpdbasename);        
        CleanUpDB(fnuturndb);
	CLOSEMMDBAPI();
	FreeRotLib();
	printf("Benchmark complete.\n\nSummary\n-------\n          Usr time  Sys time\n          --------  --------\n");
	printf("Maketrj  %9.3f %9.3f\nFoldtraj %9.3f %9.3f\n\n",CPUTimeGetUser(pTime2)-CPUTimeGetUser(pTime1),CPUTimeGetSys(pTime2)-CPUTimeGetSys(pTime1),CPUTimeGetUser(pTime4)-CPUTimeGetUser(pTime3),CPUTimeGetSys(pTime4)-CPUTimeGetSys(pTime3));
        printf("Benchtraj successful\n");
	return ERR_SUCCESS;
}



Int2 Main()
{
static Char pcUname[2000]="";
/* Arg -l to keep the log file */  

Int2 retval = 0;
#ifdef OS_UNIX
    struct utsname utsbuf;
#endif

#define NUMARGS 1

static Nlm_Arg myargs[NUMARGS] = {
/*0*/ { "Keep error.log file?", "F", NULL, NULL, TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL} };


        if ( !Nlm_GetArgs("Benchtraj", NUMARGS, myargs) ) return 1;
 	
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);

#ifdef OS_UNIX
	    /* use uname to get sysinfo */
		if (uname(&utsbuf)>=0) {
	        sprintf(pcUname,"%s %s %s %s %s\n",utsbuf.sysname,utsbuf.nodename,utsbuf.release,utsbuf.version,utsbuf.machine);
		}
#else
#ifdef OS_MSWIN
	    /* assume ver is available */
		if (GetMSWinVerInfo(pcUname)!=ERR_SUCCESS) {
            StringCpy(pcUname,"Not Available\n");
	
	    }       
#endif
#endif

		printf("OS Reported:\n%s\n",pcUname);

        retval = TRADESBenchmark();
        if (!myargs[0].intvalue) FileRemove("error.log");

        return retval;
}

