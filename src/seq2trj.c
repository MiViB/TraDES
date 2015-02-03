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
#include "seq2trj.h"

/* Global Variables */

Args myargs [9] =
{
  
/* Sequence Input - maketrj has been carved into two calling programs, one for sequence, the other for structure and unfoldtraj */
/*0*/  {"REQUIRED Input FASTA Amino Acid Sequence File Name\n      (First line must start with >, Protein sequence on second line)\n ",NULL,NULL,NULL,FALSE, 'f', ARG_FILE_IN,0.0,0,NULL},
/*1*/  {"Trajectory Distribution Type:\n      Default is 4 for 3-State Sec. Str. Prediction (GOR Method)\n      1 = Uniform, 2 = Standard, 3 = One-State Sec. Str. Prediction,\n","4","1","4",TRUE,'t',ARG_INT,0.0,0,NULL},
/*2*/  {"Trajectory Distribution ALL-COIL sampling.\nUse this for Intrinsic Disorder Simulation\n      (T/F) (T overrides -t)","FALSE",NULL, NULL, TRUE,'c',ARG_BOOLEAN,0.0,0,NULL},
/*3*/  {"Trajectory Distribution ALL-BETA sampling.\nUse this for Urea/GdHCL Denatured Simulation\n      (T/F) (overrides -t)\n      NOTE: -c T -b T results in 50% COIL, 50% BETA.","FALSE",NULL, NULL, TRUE,'b',ARG_BOOLEAN,0.0,0,NULL},
   
/*4*/  {"Input 3-State Secondary Structure File .ss : ",NULL,NULL,NULL,TRUE,'s',ARG_FILE_IN,0.0,0,NULL},
/*5*/  {"Input Distance Constraint File Name : ",NULL,NULL,NULL,TRUE,'d',ARG_FILE_IN,0.0,0,NULL},

/* Output */ 
/*6*/  {"REQUIRED Output TRJ, SS and ARA File Name (No extension)",NULL,NULL,NULL,FALSE,'o',ARG_FILE_OUT,0.0,0,NULL},
/*.trj is the output Trajectory File for structrj sampling,\n     .ss is any Secondary Structure weighting created interally,\n     .ara is full-area-at-half-max at each residue." */
/*7*/  {"Output .CSV Formatted Trajectory File (T/F) ", "FALSE", NULL, NULL, TRUE,'v',ARG_BOOLEAN,0.0,0,NULL},

/* General */
/*8*/  {"Quiet Operation (T/F)","TRUE",NULL, NULL, TRUE, 'q', ARG_BOOLEAN,0.0,0,NULL}
};


Int2 Nlm_Main() {
  Char ftmp[PATH_MAX];
  FILE *fp,*RWStatus;
  Int2 err;
  Int4 IndexArgs=-1;
  CharPtr pcThis=NULL;
  time_t time_now;
  CharPtr TimeNowStr;
  pMakeTrjParamBlock maketrjprm;
  Char errorfile[PATH_MAX];
  Char value[PATH_MAX];
  Char ssfile[PATH_MAX];

  ErrSetLogfile("error_seq2trj.log", ELOG_APPEND|ELOG_BANNER);
  ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
  ErrSetMessageLevel(SEV_FATAL);
  ErrSetLogLevel(SEV_ERROR);

  RandomSeed(GetSecs());


  if (!GetArgs("seq2trj: FASTA sequence to TraDES Trajectory Distribution *.trj.\nUse this to set up the trades input file from a sequence.\n\n",9,myargs)) {
	  printf("Use seq2trj - for args list.\n");
    return 1;       /* insufficient arguments entered */
  }
  
  if (myargs[8].intvalue==FALSE)
      traj_quiet = (Byte)VERBOSITY_VERBOSE;
  if (myargs[8].intvalue == TRUE)
      traj_quiet = (Byte)VERBOSITY_QUIET;


  if (!OPENMMDBAPI(0,NULL)){
   	ErrPostEx(SEV_FATAL,2,1, "Unable to open MMDBAPI");
    return 1;
  }
  /* call the main function */
  /* after filling in parameter block */

    maketrjprm=(pMakeTrjParamBlock)MemNew(sizeof(MakeTrjParamBlock));
	maketrjprm->TrajMethod=2; /* from sequence */
	maketrjprm->valfnamin=NULL; /* structure values are zeroed out */
	maketrjprm->pcChain=0;
	maketrjprm->modelnum=0;
	maketrjprm->startres=0;
	maketrjprm->endres=0;
	maketrjprm->peakheight=100;
	maketrjprm->noise=0;
	maketrjprm->sigma_x=0.0;
	maketrjprm->sigma_y=0.0;
	maketrjprm->savechi=0;
	maketrjprm->temperature=0;
	maketrjprm->timestep=0;
	maketrjprm->seqfnamin=myargs[0].strvalue; /* input sequence file */ 
	maketrjprm->tgtype=myargs[1].intvalue; /* handle override condtions */
	maketrjprm->DoRPS=0;
	if (myargs[4].strvalue == NULL) {
	  maketrjprm->SStruinput=FALSE; /* this is TRUE if ss file input, FALSE if ss file to be output */
      /* construct output filename */
      sprintf(ssfile,"%s.ss",myargs[6].strvalue);
	  maketrjprm->sstrufname=ssfile; /* file name to write out */
	}
	else {
	  maketrjprm->SStruinput=TRUE; 
	  maketrjprm->sstrufname=myargs[4].strvalue; /* file name to read in */
	}
    maketrjprm->trjfnamout=myargs[6].strvalue; /* non optional parameter now - must be given */
	maketrjprm->comprtype=USE_RLE; /* RLE compression */
	maketrjprm->constrfnamin=myargs[5].strvalue;
	maketrjprm->uturn=0;
	maketrjprm->units=UNITS_ARBITRARY; /* the =j reserved in original maketrj */
	maketrjprm->templat=NULL;
	maketrjprm->alignfile=NULL;
	maketrjprm->zhangwindsize = 0;
	maketrjprm->ssmask=NULL;
	maketrjprm->shiftgaps=(Boolean)(0);
	maketrjprm->fragfile=NULL;
	maketrjprm->dumpcsv=(Boolean)(myargs[7].intvalue);
	maketrjprm->benchmark=0;
	maketrjprm->all_coil = (Boolean)(myargs[2].intvalue);
	maketrjprm->all_beta = (Boolean)(myargs[3].intvalue);
		
	MakeTrj((VoidPtr)maketrjprm);
        maketrjprm=MemFree(maketrjprm);

	CLOSEMMDBAPI();
	     return 0;

}

