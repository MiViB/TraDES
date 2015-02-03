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


#include <foldtrajlib.h>

#ifdef OS_MSWIN
#define TASK_CFGFILE ".\\foldtraj"
#else
#define TASK_CFGFILE "foldtraj"
#endif

/* load in optional parameters from disk */

void LoadParams(void)
{
	Char param[PATH_MAX];
	Char value[PATH_MAX];
	long lTIMEOUT = 100;
	long lTRAJDIV = 400;
	double fBACKBONE_ERROR_TOLERANCE = 50.0;
	double fBACKBONE_PRECISION = 0.5;
	double fINCSIZE = 20.0;
	double fTUNNEL_PROB = 0.0;
	double fSTART_BACKBONE = 250.0;
	int iNUM_ROT_TRIES = 7;
	double fATOM_BOUNCINESS_BB = 0.25;
	double fATOM_BOUNCINESS_SC = 0.50;
	double fMARKOV_SCALE_FACTOR = 0.0;

	/* set to defaults first */
	TIMEOUT=100;
	BACKBONE_ERROR_TOLERANCE=60.0;
	BACKBONE_PRECISION=0.5;
	INCSIZE=20.0;
	TUNNEL_PROB=0.0;
	START_BACKBONE=250.0;
	TRAJDIV=400;
	NUM_ROT_TRIES=7;
	ATOM_BOUNCINESS_BB=0.25;
	ATOM_BOUNCINESS_SC=0.50;
	MARKOV_SCALE_FACTOR=0.0;
	BUMPCHECK_HYDROGEN=1;
	WALKTYPE=WALK_UNKNOWN;
	value[0]='\0';
	StringCpy(param,"datafilepath");
	FindPath(TASK_CFGFILE,"LOCAL",param,value,PATH_MAX);
	if (value[0]!='\0') {
		StringCpy(CFG_local_datafilepath,value);
	}
	else {
		if (StringLen(CFG_local_datafilepath)==0)
#ifdef OS_UNIX
			StringCpy(CFG_local_datafilepath,"./");
#else
			StringCpy(CFG_local_datafilepath,".\\");
#endif
	}
	StringCpy(TRAJ3FNAME,CFG_local_datafilepath);
      
	StringCpy(param,"TIMEOUT");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		lTIMEOUT=atol(value);
		if (lTIMEOUT<4L || lTIMEOUT>250L) goto badparm;
	        TIMEOUT = (Int4) lTIMEOUT;
	}
	StringCpy(param,"BACKBONE_ERROR_TOLERANCE");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fBACKBONE_ERROR_TOLERANCE=atof(value);
		if (fBACKBONE_ERROR_TOLERANCE<10.0) goto badparm;
	        BACKBONE_ERROR_TOLERANCE = (FloatLo) fBACKBONE_ERROR_TOLERANCE;
	}
	StringCpy(param,"BACKBONE_PRECISION");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fBACKBONE_PRECISION=atof(value);
		if (fBACKBONE_PRECISION<0.0005 || fBACKBONE_PRECISION>5.0) goto badparm;
	        BACKBONE_PRECISION = (FloatLo) fBACKBONE_PRECISION;
	}
	StringCpy(param,"MARKOV_SCALE_FACTOR");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fMARKOV_SCALE_FACTOR=atof(value);
		if (fMARKOV_SCALE_FACTOR<0.0 || fMARKOV_SCALE_FACTOR>1.0) goto badparm;
                MARKOV_SCALE_FACTOR = (FloatLo) fMARKOV_SCALE_FACTOR;
        }
	StringCpy(param,"INCSIZE");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fINCSIZE=atof(value);
		if (fINCSIZE<1.0 || fINCSIZE>89.0) goto badparm;
                INCSIZE = (FloatLo) fINCSIZE;
        }
	StringCpy(param,"TUNNEL_PROB");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fTUNNEL_PROB=atof(value);
		if (fTUNNEL_PROB<0.0 || fTUNNEL_PROB>1.0) goto badparm;
                TUNNEL_PROB = (FloatLo) fTUNNEL_PROB;
        }
	StringCpy(param,"START_BACKBONE");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fSTART_BACKBONE=atof(value);
		if (fSTART_BACKBONE<=0.0 || fSTART_BACKBONE>360.0) goto badparm;
                START_BACKBONE = (FloatLo) fSTART_BACKBONE;
        }	
	StringCpy(param,"NUM_ROT_TRIES");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		iNUM_ROT_TRIES=atoi(value);
		if (iNUM_ROT_TRIES<1 || iNUM_ROT_TRIES>100) goto badparm;
                NUM_ROT_TRIES = (Int2) iNUM_ROT_TRIES;
        }	
	StringCpy(param,"BUMPCHECK_HYDROGEN");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		if (!StringCmp(value,"TRUE") || !StringCmp(value,"True") || !StringCmp(value,"true"))
			BUMPCHECK_HYDROGEN=1;
		else if (!StringCmp(value,"FALSE") || !StringCmp(value,"False") || !StringCmp(value,"false"))
			BUMPCHECK_HYDROGEN=0;
		else goto badparm;
        }	
	StringCpy(param,"WALKTYPE");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		if (!StringCmp(value,"CA")) {
			WALKTYPE=WALK_CA;
			StringCat(TRAJ3FNAME,"cawalk_dict");
		}
		else if (!StringCmp(value,"PHIPSI")) {
			WALKTYPE=WALK_PHIPSI;
			StringCat(TRAJ3FNAME,"phipsiwalk_dict");
		}
		else goto badparm;
        }	
	StringCpy(param,"ATOM_BOUNCINESS_BB");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fATOM_BOUNCINESS_BB=atof(value);
		if (fATOM_BOUNCINESS_BB<0.0) goto badparm;
                ATOM_BOUNCINESS_BB = (FloatLo) fATOM_BOUNCINESS_BB;
        }	
	StringCpy(param,"ATOM_BOUNCINESS_SC");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		fATOM_BOUNCINESS_SC=atof(value);
		if (fATOM_BOUNCINESS_SC<0.0) goto badparm;
                ATOM_BOUNCINESS_SC = (FloatLo) fATOM_BOUNCINESS_SC;
        }	
	StringCpy(param,"TRAJDIV");
	GetAppParam(TASK_CFGFILE,"INITTRAJ",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		lTRAJDIV=atol(value);
		if (lTRAJDIV<100L) goto badparm;
                TRAJDIV = (Int4) lTRAJDIV;
        }
	StringCpy(param,"CDDdbname");
	GetAppParam(TASK_CFGFILE,"LOCAL",param,"",value,PATH_MAX);
	if (value[0]!='\0') {
		StringCpy(CFG_local_cddname,value);
  }
  else
		StringCpy(CFG_local_cddname,"All");
	StringCpy(param,"CDDpath");
	FindPath(TASK_CFGFILE,"LOCAL",param,value,PATH_MAX);
	if (value[0]!='\0') {
		StringCpy(CFG_local_cddpath,value);
  }
  else
#ifdef OS_UNIX
		StringCpy(CFG_local_cddpath,"./");
#else
		StringCpy(CFG_local_cddpath,".\\");
#endif
	return;
badparm:
	PurgeGlobs();
	ErrPostEx(SEV_FATAL,9,8,"%s out of valid range.  Please correct and try again.",param);
}

void printgraph(PTGS ptgsHere, FILE *f)
  {
    Int4 cnt; 

  for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
        fprintf(f,"%d, ",(int)(ptgsHere->TrajGraph)[cnt]);
    }
  fprintf(f,"\n");
  }



/* generate initial trajectory graphs */
void MakeTG(CharPtr kbsequence, pMakeTrjParamBlock paramblock, Char*fnam)
{	Int4 clength,cnt,c1=0,cnt2; /*,c2,res;*/
	/*Int4 cntr;*/
	Int4 dictsize;
	/*FloatLo csum,ctemp=0.0,remain;*/
	Char cres[2];
	Char buf[255];
	/*FloatLo ftmp1=1.0,ftmp2=1.0,ftmp3=1.0;*/
	PTGS ptgsHere=NULL,ptgsTmp;
	PTGS ptgsHereH=NULL;
	/*PTGS ptgsHereE=NULL;
	PTGS ptgsHereC=NULL;*/
	Int2 trajwidth,numturns;
	/* store trajectory graphs */
	Int4 PNTR TrajGraph;
	Int4 PNTR TrajCisProGraph;
	Int4 TrajCisProIntegral;
	/* used for GOR structure prediction scheme */
	FloatLo arwProbMatrix[MAXRES][MAXSTRUCT];
	Char arwFiltStruct[MAXRES];
	CharPtr tmpfnam;
	time_t tm1, tm2;
	FILE *f, *csv = NULL, *ara = NULL;
	Char sequence[MAXSIZE];
	Char encseq[MAXSIZE];
	Char xxsequence[MAXSIZE];
	FloatLo cutoff,p;
	ByteStorePtr bsSeq;
	ValNodePtr r14vnp=NULL,FTvnp=NULL,FiltFTvnp=NULL,newTDvnp=NULL,turnResultvnp=NULL,vnpCur;
	CharPtr pcR14calls=NULL;
	FloatLoPtr fpR14score=NULL;
	Int2 uturnpos;
/*	UturnStrucPtr usp; */
	PFDS pfdsHere;
	NN nn;
	PNN pnnHere,pnnNext;
	Boolean DoUturn=FALSE;
	Boolean ulr;
	Int2 TrajType;
	CharPtr GORfnam = NULL;
	Int4 ctype;
	Boolean OverrideGOR;
	Boolean Uturn;
	CharPtr fragfile;
	Boolean Benchmark;
	Int4 i;
	
/* get rid of the long list of parameters passed and just pass the parameter block... */
	if (paramblock == NULL) return;
	TrajType=paramblock->tgtype;
	GORfnam=paramblock->sstrufname;
	ctype=paramblock->comprtype;
	OverrideGOR=(Boolean) paramblock->SStruinput;  /* this is TRUE if ss file input, FALSE if ss file to be output */
/*	Uturn= (Boolean) paramblock->uturn; */
	Uturn = FALSE; /*CWVH 2012 Disabled */
	fragfile=paramblock->fragfile;
	Benchmark=paramblock->benchmark;
	
	
	
	/* initialize database and get graph size from first dictionary entry */
	if ((TrajType!=TRAJ_UNIFORM)) {
		ulr=USE_LOTS_RAM;
		USE_LOTS_RAM=FALSE;
		TGInit(TRAJ3FNAME,DB_READ,&dictsize);
		if (dictsize!=60) {
			PurgeGlobs();
			ErrPostEx(SEV_FATAL,1,1,"Invalid trajectory distribution dictionary file: %s",TRAJ3FNAME);
		}
        	ptgsHereH=TrajGraphRead(1);
		trajwidth=ptgsHereH->dim;
		ptgsHereH=FreeTraj(ptgsHereH);
		TGClose();
		USE_LOTS_RAM=ulr;
	}
	else {
		trajwidth=TRAJDIV;
	}
	/* use parents of modified AAs for secondary structure prediction */
	clength=ConvertExtSequence(kbsequence,sequence,EXTAA_PARENT);
	/* will need to distinguish special amino acids in some cases
	   though */ 
	clength=ConvertExtSequence(kbsequence,xxsequence,EXTAA_X);
	/* allocate trajectory graph for actual protein */
	TrajGraph=(Int4 PNTR)MemNew(sizeof(Int4)*trajwidth*trajwidth);
	/* allocate an extra graph for CisPro */
	TrajCisProGraph=(Int4 PNTR)MemNew(sizeof(Int4)*trajwidth*trajwidth);
	/* stores trajectory graph integrals for actual protein */
	TrajCisProIntegral=0;

	if (paramblock->all_coil == TRUE) {
		TrajType = TRAJ_GOR; /* switched from TRAJ_SSTRU to TRAJ_GOR as it matches *.ss input in Peak height output  CWVH 12 June 2012 */
		paramblock->tgtype=TrajType;
		OverrideGOR = FALSE;
	}
	if (paramblock->all_beta == TRUE) {
		TrajType = TRAJ_GOR;
		paramblock->tgtype=TrajType;
		OverrideGOR = FALSE;
	}
	
	if ((TrajType==TRAJ_GOR) || (TrajType==TRAJ_SSTRU)) {
		if (OverrideGOR==TRUE) {  /* read in ss file */
  			if ((f=FileOpen(GORfnam,"r"))==NULL) {
				PurgeGlobs();
		                ErrPostEx(SEV_FATAL,1,1,"Error - unable to open input SS file %s",GORfnam);
                		return;
        		}
			c1=0;
			while (FileGets(buf,255,f)!= NULL) {
				if (buf[0]!='#') {
					if (StringLen(buf) > 5) {
						sscanf(buf,"%*c %c %f %f %f",&(arwFiltStruct[c1]),&(arwProbMatrix[c1][0]),&(arwProbMatrix[c1][1]),&(arwProbMatrix[c1][2]));
						if (floor(arwProbMatrix[c1][0]+arwProbMatrix[c1][1]+arwProbMatrix[c1][2]+0.5)!=100.0) {
							FileClose(f);
							PurgeGlobs();
							ErrPostEx(SEV_FATAL,1,1,"Probabilities in SS file do not sum to 100%% in %s, residue %ld",GORfnam,c1+1);
							return;
						}
					} /* terminal cr/lf issue CWVH June 12 2012 */
					c1++;
				}
			}
			FileClose(f);
		} /* Input SS prediction */
		else {
/* no prediction cases */
			  if ((paramblock->all_coil == TRUE) || (paramblock->all_beta == TRUE)) {
/* loop to end of sequence fill in values */
				  for (i=0;i<clength;i++) {
				  if ((paramblock->all_coil == TRUE) && (paramblock->all_beta == TRUE)) {
				    arwProbMatrix[i][0] = 0.0;
					arwProbMatrix[i][1] = 50.0;
					arwProbMatrix[i][2] = 50.0;
					arwFiltStruct[i] = 'C';
				  }
				  else if (paramblock->all_coil == TRUE) { /* all coil */
				    arwProbMatrix[i][0] = 0.0;
					arwProbMatrix[i][1] = 0.0;
					arwProbMatrix[i][2] = 100.0;
					arwFiltStruct[i] = 'C';

				  } else { /* all beta */
				    arwProbMatrix[i][0] = 0.0;
					arwProbMatrix[i][1] = 100.0;
					arwProbMatrix[i][2] = 0.0;
					arwFiltStruct[i] = 'E';				  
			      }
				  }
/*CWVH always write out the SS conditions used for recordkeeping */
  				if (GORfnam!=NULL) {
	  				if ((f=FileOpen(GORfnam,"w"))==NULL) {
						PurgeGlobs();
				             ErrPostEx(SEV_FATAL,1,1,"Error - unable to open output file %s",GORfnam);
				     		return;
						}
					fprintf(f,"#Trades Secondary Structure Sampling Method Selected:");
                        if (paramblock->all_coil == TRUE)  fprintf(f,"All-coil -c T ");	
						if (paramblock->all_beta == TRUE)  fprintf(f,"All-beta -b T ");
						fprintf(f,"\n#\n");
	
						fprintf(f,"#Seq.\tPred.\t H\t E\t C\n");
					for (c1=0;c1<clength;c1++)
						fprintf(f,"%c\t  %c\t%3.0f\t%3.0f\t%3.0f\n",sequence[c1],arwFiltStruct[c1],arwProbMatrix[c1][0],arwProbMatrix[c1][1],arwProbMatrix[c1][2]);
					FileClose(f);
				}
			   OverrideGOR = TRUE;
			   paramblock->SStruinput = OverrideGOR; /* make it appear like we read in a *.ss file 12 June 2012 CWVH */
			  } /* no prediction cases */
			  else { 
/* call GOR */
				tm1=GetSecs();
				if (traj_quiet==VERBOSITY_VERBOSE)
					printf("Predicting secondary structure via GOR method..\n");
				/* predict secondary structure with GOR method */
				if (fnGilGOR(sequence,arwProbMatrix,arwFiltStruct)!=ERR_SUCCESS) {
					PurgeGlobs();
						    ErrPostEx(SEV_FATAL,1,2,"Error occurred during GOR prediction");
	               			return;
				}
				tm2=GetSecs();
				if (traj_quiet==VERBOSITY_VERBOSE)
					printf("GOR prediction: %ld seconds\n",(long int)(tm2-tm1));
				/* FiltStruct is a char array giving the ss prediction, while
					probmatrix is a Nx3 array giving probabilities of H, E, C */
				if (GORfnam!=NULL) {
	  				if ((f=FileOpen(GORfnam,"w"))==NULL) {
						PurgeGlobs();
				             ErrPostEx(SEV_FATAL,1,1,"Error - unable to open output file %s",GORfnam);
				     		return;
						}
					fprintf(f,"#Predicted Secondary Structure by GOR method\n#\n");
					fprintf(f,"#Seq.\tPred.\t H\t E\t C\n");
					for (c1=0;c1<clength;c1++)
						fprintf(f,"%c\t  %c\t%3.0f\t%3.0f\t%3.0f\n",sequence[c1],arwFiltStruct[c1],100.0*arwProbMatrix[c1][0],100.0*arwProbMatrix[c1][1],100.0*arwProbMatrix[c1][2]);
					FileClose(f);
				}
			  }	
		}
	}

	/* to make procedure faster, would be better to integrate before
	   saving initial trajectory graph, wouldn't need to integrate
	   every time then - not needed though, it's fast already */
	cres[1]=0;
	/* begin storing actual trajectory graphs */
	/* initialize new CodeBase Database */
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Creating new trajectory graph database..\n");
	/* don't use NewTraj, don't want to allocate a tg buffer */
	ptgsHere=(PTGS)MemNew(sizeof(TGS));
	/* write out trajectories */
	
	if (traj_quiet==VERBOSITY_QUIET || traj_quiet==VERBOSITY_VERBOSE)
		printf("Writing out trajectory distributions..\nComputing trajectory distribution          ");
	tm1=GetSecs();
	
/* added by CWVH to dump CSV file of TG for R analysis */
    if(paramblock->dumpcsv == TRUE) {
    	tmpfnam=(CharPtr)MemNew(sizeof(Char)*(StringLen(CSV_EXT)+StringLen(fnam)+2));
		StringCpy(tmpfnam,fnam);
		StringCat(tmpfnam,CSV_EXT);
		if ((csv=FileOpen(tmpfnam,"w"))==NULL) {
	        ErrPostEx(SEV_ERROR,1,1,"Error - unable to open CSV output file %s",tmpfnam);
	        return;
    	}
	  MemFree(tmpfnam);
	}
	
/* added by CWVH - ARA file output is default together with TRJ file */
	if(!Benchmark) {	
			tmpfnam=(CharPtr)MemNew(sizeof(Char)*(StringLen(ARA_EXT)+StringLen(fnam)+2));
			StringCpy(tmpfnam,fnam);
			StringCat(tmpfnam,ARA_EXT);
			if ((ara=FileOpen(tmpfnam,"w"))==NULL) {
		       ErrPostEx(SEV_ERROR,1,1,"Error - unable to open ARA output file %s",tmpfnam);
		       return;
			}
			MemFree(tmpfnam);
	}

	
	for (cnt=0;cnt<clength;cnt++) {
		/* make a separate trajectory graph for cis-Pro, any SStype */
		/* integration destroys this so must rebuild each time a
		   proline is encountered */
		if (sequence[cnt]=='P')
			TrajCisProIntegral=FillTG(TrajType,"P",0,NULL,trajwidth,TrajCisProGraph,TRUE);
		if (traj_quiet==VERBOSITY_QUIET || traj_quiet==VERBOSITY_VERBOSE) {
			printf("\b\b\b\b\b\b\b\b\b%4ld/%-4ld",(long int)(cnt+1),(long int)clength);
			fflush(stdout);
		}
		if (traj_quiet==VERBOSITY_SILENT) {
			ProgramProgress=cnt+1;
		}
		/* use trajectory graph of parent residues too */
		cres[0]=sequence[cnt];
		ptgsHere->TrajIntegral=FillTG(TrajType,cres,arwFiltStruct[cnt],arwProbMatrix[cnt],trajwidth,TrajGraph,FALSE);
		/* re-use same structure since only writing one at a time */
		ptgsHere->resnum=cnt+1;
		/* use pointer arithmetic to point to cnt+1'th
		   trajectory graph */
		ptgsHere->TrajGraph=TrajGraph;
		if (sequence[cnt]=='P') {
			ptgsHere->CisTrajGraph=TrajCisProGraph;
			ptgsHere->CisTrajIntegral=TrajCisProIntegral;
		}
		else {
			ptgsHere->CisTrajGraph=NULL;
			ptgsHere->CisTrajIntegral=0;
		}
		ptgsHere->pCis=0.0;
		if (sequence[cnt+1]=='P') {
			ptgsHere->pCis=P_CISP;
		}
		/* assume modified cysteines don't bridge */
		if (xxsequence[cnt]=='C') {
			ptgsHere->pSS=(InDisulphide(cnt+1));
			if (ptgsHere->pSS==0.0)
				ptgsHere->pSS=P_SS;
		}
		else
			ptgsHere->pSS=0.0;
		ptgsHere->ChiWMean=CHI_W;
		ptgsHere->ChiWSD=CHISD_W;
		/* store X's for special amino acids */
		ptgsHere->AA=xxsequence[cnt];
		ptgsHere->dim=trajwidth;
		/* fill in remaining fields */
/* CWVH peak is assigned in TrajCalcSparsity */		
		/* Added by CWVH to dump csv files */
		if(paramblock->dumpcsv == TRUE) {
		   printgraph(ptgsHere,csv);
		}
		TrajCalcSparsity(ptgsHere,0,ara);
		/* make timeout depend on trajectory graph sparsity */
		TrajCalcTout(ptgsHere);
		/* integrate them */
		/***** IMPORTANT ****/
		/* integrate them before calculating NZ and after sparsity */
		TrajGraphIntegrate(ptgsHere);
		TrajCalcNZ(ptgsHere);
		/* rotamer info is reserved for future use */
		ptgsHere->rotid=0;
		ptgsHere->markovsf=MARKOV_SCALE_FACTOR;
	
		/* create/open database here to avoid interferece with database being opened
		   above in the FillTG call */
		ulr=USE_LOTS_RAM;
		USE_LOTS_RAM=FALSE;
		if (TGInit(tmpdbasename,DB_CREATE,NULL)!=ERR_SUCCESS) {
			PurgeGlobs();
			ErrPostEx(SEV_FATAL,11,1,"Unable to open database, cannot continue");
			return;
		}
		/* add to database */
		if ((TrajType==TRAJ_UNIFORM) && (ctype==USE_BZ))
			TrajGraphWrite(ptgsHere,USE_RLE,FALSE);
		else
			TrajGraphWrite(ptgsHere,ctype,FALSE);
		
		TGClose();
		USE_LOTS_RAM=ulr;
	}
	ptgsHere=MemFree(ptgsHere);
	/* close off all trajectory database files */
	fflush(csv);
	FileClose(csv);
/* Close the ARA file */	
	if (ara != NULL){
		fflush(ara);
	    FileClose(ara);
	}
	tm2=GetSecs();
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("\nTrajectory distribution creation: %ld seconds \n\n",(long int)(tm2-tm1));
	/* add fragments if any were given */
	if (fragfile) {
		if (traj_quiet==VERBOSITY_VERBOSE)
			printf("Adding fragments from file %s..\n",fragfile);
		if (ImportFrags(fragfile,xxsequence)!=ERR_SUCCESS)
			ErrPostEx(SEV_ERROR,1,1,"Unable to find fragment file %s",fragfile);
	}
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Storing archive in ASN.1 format..\n");
	ConvertExtSequence(kbsequence,encseq,EXTAA_ENCODE);
	if (ctype==USE_NONE)
		PackAsnTrajGraph(fnam,encseq,NULL,NULL,NULL,TRUE);
	else
		PackAsnTrajGraph(fnam,encseq,NULL,NULL,NULL,FALSE);
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Removing database files..\n");
	tmpfnam=(CharPtr)MemNew(sizeof(Char)*(StringLen(DB_EXT1)+StringLen(fnam)+2));
	CleanUpDB(tmpdbasename);
	tmpfnam=MemFree(tmpfnam);
	/* free trajectory graphs */
	TrajGraph=MemFree(TrajGraph);
	TrajCisProGraph=MemFree(TrajCisProGraph);
/*	FreeResUsedTable();*/
	
}


/*  
$Log: maketg.c,v $
Revision 1.35  2004/09/24 20:31:37  hfeldman
Print error message if fragment file not found

Revision 1.34  2004/09/24 19:09:36  hfeldman
Added import fragment option to maketrj

Revision 1.33  2004/07/06 20:15:58  hfeldman
Removed some old unneeded dependencies

Revision 1.32  2003/11/07 18:11:51  feldman
Fixed magic numbers

Revision 1.31  2003/10/17 16:46:46  ksnyder
Reset CFG_local_datafilepath back to its original value

Revision 1.30  2003/10/17 16:13:12  ksnyder
Fixed error in reading the output path from the foldtraj config file

Revision 1.29  2003/08/01 18:05:42  feldman
Changed maketrj defaults to more reasonable

Revision 1.28  2003/07/14 20:01:01  egarderm
Changed TASK_CFGFILE to .\foldtraj on Windows to allow local INI file

Revision 1.27  2003/03/14 21:06:01  feldman
Now turning off USE_LOTS_RAM during writing on TDs to avoid loading useless stuff into RAM

Revision 1.26  2003/01/09 20:00:25  feldman
Allow bailing from maketrj on signal

Revision 1.25  2002/07/25 16:31:55  feldman
Added tunnel prob

Revision 1.24  2001/09/13 16:22:42  feldman
Minor fragment related bugfixes

Revision 1.23  2001/09/10 20:32:48  feldman
Minor fix

Revision 1.19  2001/08/31 17:42:53  feldman
Added objective U-turn prediction when no template is provided

Revision 1.18  2001/06/18 16:05:49  feldman
Don't set data path if already set elsewhere

Revision 1.17  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.16  2001/04/04 21:26:01  feldman
-Removed my BSPrintf, now using the one in SLRIlib - thus needed to add
	this library to makefiles
-changed all fgets to FileGets

Revision 1.15  2001/03/30 22:22:27  feldman
Minor spelling and grammar fixes, bug fixes, and added
foldtraj ASCII screensaver ability

Revision 1.14  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.13  2001/03/23 18:44:17  feldman
Integrated foldtraj into vistraj (!)
and cleaned up a few minor bugs

Revision 1.12  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.11  2001/02/26 22:21:15  feldman
-Changed fragments to allow multiple possible fragments at the same
residue
-altered random walk to make use of fragments when present (a first
attempt at least...)

Revision 1.10  2001/02/15 20:27:57  feldman
Reworked maektrj to take a parameter block instead of a bunch
of individual parameters.

Revision 1.9  2001/02/06 18:40:39  feldman
Added a few functions for dealing with distance constraints and
tidied up so maketrj could join the library (foldtrajlib) without
conflicts

Revision 1.8  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.7  2001/01/12 20:01:51  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.6  2000/08/09 19:12:19  feldman
-minor bugfix update and fixed up makefiles removing USEDSSP
-added bioseq/seq-entry to trajectory graph header when using
 unfoldtraj or val2trj

Revision 1.5  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.4  2000/07/10 15:40:58  feldman
Updated TrajgraphWrite to not call TGInit and TGClose, thus
now just call TGInit once at start of program, TGClose at end
(also removed these calls from updatePcis and alterresiduesequence
)

Revision 1.3  2000/06/20 16:40:23  feldman
Incorporated sstru extended structure calculation into C code,
no longer need an external program

Revision 1.2  2000/06/15 17:11:03  feldman
Added Replace option when TrajGraphWrite-ing

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

