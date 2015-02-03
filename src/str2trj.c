#include "hfprogs.h"

/* Global Variables */


Args myargs [12] =
{
/* 3D Structure Input - maketrj has been carved into two calling programs, one for sequence, the other for structure and unfoldtraj */
/*0*/  {"REQUIRED \nInput 3D Asn.1 Structure File Name: (if no extension - looks for .val, .cn3, .prt)",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
/*1*/  {"PDB Chain Name (default = first one found)",NULL,NULL,NULL,TRUE,'c',ARG_STRING,0.0,0,NULL},
/*2*/  {"Model Number (default = 1)","1","1","9999",TRUE,'m',ARG_INT,0.0,0,NULL},
/*3*/  {"Start Residue (default = 1)","1","1","9999",TRUE,'s',ARG_INT,0.0,0,NULL},
/*4*/  {"End Residue (default = 0 represents last residue).","0","0","9999",TRUE,'e',ARG_INT,0.0,0,NULL},
/*5*/  {"Save Sidechain Rotamer angles\n      (Default 0=none; 1=all, 2=buried)","0","0","2",TRUE,'r',ARG_INT,0.0,0,NULL},

/* Output */ 
/*6*/  {"Output TRJ File Name: (default = {inputfile}.trj) ",NULL,NULL,NULL,TRUE,'o',ARG_FILE_OUT,0.0,0,NULL},

/* UnfoldTraj */
/*7*/  {"Unfolding MODE - Temperature (Kelvin)","0","0","9999.0",TRUE, 't', ARG_FLOAT,0.0,0,NULL},
/*8*/  {"Unfolding MODE - Time step (fs)\nUNFOLDING MODE Example: -t 375.0 -d 100.0\n","0","0","9999.0",TRUE, 'd',ARG_FLOAT,0.0,0,NULL},
/*9*/  {"Peak Width MODE - Standard deviation for x (Degrees)","0","0","45.0",TRUE,'x',ARG_FLOAT,0.0,0,NULL},
/*10*/ {"Peak Width MODE - Standard deviation for y (Degrees)\nPEAK WIDTH MODE Example: -x 5 -y 5\n","0","0","22.0",TRUE,'y',ARG_FLOAT,0.0,0,NULL},

/* General */
/*11*/  {"Quiet Operation (T/F)","TRUE",NULL, NULL, TRUE, 'q', ARG_BOOLEAN,0.0,0,NULL}
};


Int2 Nlm_Main() {
  Char fnamin[PATH_MAX];
  Char outfile[PATH_MAX];
  pMakeTrjParamBlock maketrjprm;
  Int4 iPeakHeight = 100;
  Int2    iDotLen = 0;

  ErrSetLogfile("error_str2trj.log", ELOG_APPEND|ELOG_BANNER);
  ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
  ErrSetMessageLevel(SEV_FATAL);
  ErrSetLogLevel(SEV_ERROR);

  RandomSeed(GetSecs());
 

  if (!GetArgs("str2trj: ASN.1 protein 3D structure converted to TraDES *.trj file\n",12,myargs)) {
	  printf("Use str2trj - for args list.\n");
    return 1;       /* insufficient arguments entered */
  }
  
 
  if (myargs[11].intvalue==FALSE)
      traj_quiet = (Byte)VERBOSITY_VERBOSE;
  if (myargs[11].intvalue == TRUE)
      traj_quiet = (Byte)VERBOSITY_QUIET;

  if (!OPENMMDBAPI(0,NULL)){
   	ErrPostEx(SEV_FATAL,2,1, "Unable to open MMDBAPI");
    return 1;
  }



     /* keep any input extension typed in */
     iDotLen = 0;
     fnamin[0] = '\0';
     iDotLen = (Int2) StringLen(myargs[0].strvalue);
     if (iDotLen > 5) {
		if (myargs[0].strvalue[iDotLen - 4] == '.') { /* extension typed, use it first one we hit */
		     StringCpy(fnamin,myargs[0].strvalue);
		     if (FileLength(fnamin) == 0)  {
		      		 ErrPostEx (SEV_FATAL, 1, 10, "Unable to find file  %s",myargs[0].strvalue);
	   			 return 10;
		    }	              
		}
     } /* but if not typed in or too short, try and add... */
     if (fnamin[0] == '\0') {
        StringCpy(fnamin,myargs[0].strvalue);
	StringCat(fnamin,".val");
	if (FileLength(fnamin) == 0) {
		     StringCpy(fnamin,myargs[0].strvalue);	              
             	     StringCat(fnamin,".cn3");
			if (FileLength(fnamin) == 0) {
		     	   StringCpy(fnamin,myargs[0].strvalue);	              
             	           StringCat(fnamin,".prt");
			  if (FileLength(fnamin) == 0) {
		      		 ErrPostEx (SEV_FATAL, 1, 11, "Unable to find file with .val .prt .cn3 extension %s",myargs[0].strvalue);
	   			 return 11;
			  }
           		}
	}
    }
    


  iDotLen = 0;
  if (myargs[6].strvalue == NULL)  /* no output file supplied, use input filename */
    StringCpy(outfile,myargs[0].strvalue);
  else 
    StringCpy(outfile,myargs[6].strvalue);
 
  iDotLen = (Int2) StringLen(outfile);
  if (iDotLen > 4) {
      if (outfile[iDotLen - 4] == '.') { /* IF extension, we remove it, interior code will add .trj  */
	     outfile[iDotLen -4] = '\0';      
	}
  }



/* automatically sets the integer peak height to produce smooth gaussian values 
   for small angles */
  if ((myargs[9].intvalue > 0) || (myargs[10].intvalue > 0))
   {
     iPeakHeight = 1000;
     if ((myargs[9].intvalue > 10) || (myargs[10].intvalue > 10))
     iPeakHeight = 500;
   }
  if ((myargs[7].floatvalue > 0.0) || (myargs[8].floatvalue > 0.0))
   {
     iPeakHeight = 1000;
     if ((myargs[9].floatvalue > 475.0) || (myargs[10].floatvalue > 1000.0))
     iPeakHeight = 500;
   } 
/* default peak height is 100 */

  /* call the main function now */
  /* after filling in parameter block */
  maketrjprm=(pMakeTrjParamBlock)MemNew(sizeof(MakeTrjParamBlock));
	maketrjprm->TrajMethod=1; /* from 3D structure val file */
	maketrjprm->valfnamin=fnamin;
	maketrjprm->pcChain=myargs[1].strvalue;
	maketrjprm->modelnum=(Int2) myargs[2].intvalue;
	maketrjprm->startres=(Int2)myargs[3].intvalue;
	maketrjprm->endres=(Int2)myargs[4].intvalue;
	maketrjprm->peakheight=iPeakHeight;
	maketrjprm->noise=0;
	maketrjprm->sigma_x=(FloatLo)myargs[9].floatvalue;
	maketrjprm->sigma_y=(FloatLo)myargs[10].floatvalue;
	maketrjprm->savechi=(Int2)myargs[5].intvalue;
	maketrjprm->temperature=(FloatLo)myargs[7].floatvalue;
	maketrjprm->timestep=(FloatLo)myargs[8].floatvalue;
	maketrjprm->seqfnamin=NULL;
	maketrjprm->tgtype=4; /* default value from maketrj for val file input*/
	maketrjprm->sstrufname=NULL;
	maketrjprm->DoRPS=0;
	maketrjprm->SStruinput=0;
	maketrjprm->trjfnamout=outfile;
	maketrjprm->comprtype=USE_RLE;
	maketrjprm->constrfnamin=NULL;
	maketrjprm->uturn=0;
	maketrjprm->units=UNITS_ARBITRARY;
	maketrjprm->templat=NULL;
	maketrjprm->alignfile=NULL;
	maketrjprm->zhangwindsize = 0;
	maketrjprm->ssmask=NULL;
	maketrjprm->shiftgaps=(Boolean)(FALSE);
	maketrjprm->fragfile=NULL;
	maketrjprm->dumpcsv=(Boolean)(FALSE);
	maketrjprm->benchmark=0;
    maketrjprm->all_coil = (Boolean)(FALSE);
	maketrjprm->all_beta = (Boolean)(FALSE);
		
	MakeTrj((VoidPtr)maketrjprm);
        maketrjprm=MemFree(maketrjprm);
	CLOSEMMDBAPI();
	     return 0;

}

