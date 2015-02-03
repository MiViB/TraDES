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
#include "hfprogs.h"
#include <string.h>

/* Global Variables */
Args Rargs[3] = {{"Input Asn.1 3D structure as .prt .val or .cn3 File Name.",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
                       {"Model Number: (0=Best Coordinate Set, 9999 = All)","0","0","9999",TRUE,'m',ARG_INT,0.0,0,NULL},
		       {"Output PDB File Name?",NULL,NULL,NULL,TRUE,'p',ARG_FILE_OUT,0.0,0,NULL}};
                       

Int2 Main()
{
    PMSD  pmsdRoot = NULL;

/* Variables for clipping status and pdb output file name*/
    static char pdbfileout[256];
    static char valfilein[PATH_MAX];
    Int2    iTest = 0;
    FILE    *pFile = NULL;
    Int2    iDotLen = 0;
	AsnIoPtr aip=NULL;
	NcbiMimeAsn1Ptr nmap=NULL;
	BiostrucPtr bsp=NULL;
	BiostrucSeqPtr bssp=NULL;
	PDNMS pdnmsModelstruc = NULL;
	Int2 modelnum = 0;
	PMLD pmldHere;
	ErrSev esMsg,esLog,esFatal;
	PDNML pdnmlModel;
	PMLD pmldThis;


    
/* Initialize MMDB-API */
        ErrSetLogfile("error_str2pdb.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);

	if (!GetArgs("str2pdb  Converts Asn.1 3d .val .prt or .cn3 file to PDB format.\nCaution - will not work on Biounits downloaded from NCBI, use Asymmetric Units\n",3,Rargs))
		return 1;
/*	printf("args0 :%s , args2: %s\n",Rargs[0].strvalue,Rargs[2].strvalue);*/

        

	
/* be forgiving about adding extensions 
  CWVH  June 2012
    No-extension is user frustration! offer forgiveness for typing the extention , punch a null overtop of the dot  
*/

        iDotLen = 0;
	if (Rargs[0].strvalue != NULL) {
           iDotLen = (Int2) StringLen(Rargs[0].strvalue);
           if (iDotLen > 5) {
		if (Rargs[0].strvalue[iDotLen - 4] == '.') {
	        Rargs[0].strvalue[iDotLen - 4] = '\0';      
		}
           }
	}


        if(Rargs[2].strvalue==NULL) {
        	StringCpy(pdbfileout,Rargs[0].strvalue);
		StringCat(pdbfileout,".pdb");
	}
	else
	 StringCpy(pdbfileout,Rargs[2].strvalue);


	printf("Output pdb file name: %s \n",pdbfileout);
	
        StringCpy(valfilein, Rargs[0].strvalue);
       	StringCat(valfilein,".val");
        if (FileLength(valfilein) == 0) 
	  {
              StringCpy(valfilein, Rargs[0].strvalue);
              StringCat(valfilein, ".cn3");
              if (FileLength(valfilein) == 0) {
	              StringCpy(valfilein, Rargs[0].strvalue);
	              StringCat(valfilein, ".prt");
	              if (FileLength(valfilein) == 0) {
			ErrPostEx(SEV_FATAL,13,1,"Unable to find input %s.val or %s.cn3 or %s.prt file",Rargs[0].strvalue,Rargs[0].strvalue,Rargs[0].strvalue);
			return 13;
		     }
		}
       	 }

	if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,12,1,"Unable to open MMDBAPI, check for missing bstdt.val dictionary file.");
		return 12;
	}
	
/*	printf("\nStructure loading ...filename: %s", valfilein); */


/* load an ASN.1 file *.val or *.c3d or *.prt write out Biostruc only */

	aip=AsnIoOpen(valfilein,"rb");
	if (aip==NULL) {
		
			ErrPostEx(SEV_FATAL,11,1,"Unable open binary or ascii ASN.1 stream in file %s",valfilein);
			return 11;
	}

/* first try biostruc load */
/* 	printf("try biostruc\n");  */
	esMsg=ErrGetMessageLevel(); esLog=ErrGetLogLevel(); esFatal=ErrGetFatalLevel();
        ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
 	bsp=BiostrucAsnRead(aip,NULL);
	AsnIoClose(aip);
  	ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);

	if (bsp == NULL) { /* try ascii */
		aip = AsnIoOpen(valfilein,"r");
		if (aip==NULL) {
			ErrPostEx(SEV_FATAL,10,1,"Unable open binary or ascii ASN.1 stream in file %s",valfilein);
			return 10;
		}
     	 	ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
	 	bsp=BiostrucAsnRead(aip,NULL);
		AsnIoClose(aip);
	  	ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);


	}


	if (bsp==NULL) {
/*		printf("try mime binary\n"); */
		/* then try NCBIMime load */
		aip=NULL;
		aip=AsnIoOpen(valfilein,"rb");
  		ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
		nmap=NcbiMimeAsn1AsnRead(aip,NULL);
  		ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);
		AsnIoClose(aip);

		if (nmap == NULL) {
		/*	printf("try mime ascii\n"); */
			aip=AsnIoOpen(valfilein,"r");
  			ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
			nmap=NcbiMimeAsn1AsnRead(aip,NULL);
  			ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);
			AsnIoClose(aip);
		}

        	if (nmap!=NULL) {
		/* got an NCBI mime */
			if (nmap->choice!=NcbiMimeAsn1_strucseq && nmap->choice!=NcbiMimeAsn1_strucseqs) {
				/* wrong MIME type */
				nmap=NcbiMimeAsn1Free(nmap);
		    		ErrPostEx(SEV_ERROR,9,1,"MIME-type wrapper is not a strucseq or strucseqs - no Biostruc to convert to PDB %s",valfilein);
				return 9;
			}
			/* unwrap the mime and leave the bsp for PDB conversion */
				bssp=(BiostrucSeqPtr)(nmap->data.ptrvalue);
				bsp=bssp->structure;
                		bssp->structure = NULL;
				nmap=NcbiMimeAsn1Free(nmap); /* discard the wrapper */
				nmap = NULL;
				bssp = NULL;
			}
		}

	if (bsp == NULL) {
	    	ErrPostEx(SEV_ERROR,8,1,"No Biostruc in files to convert to PDB %s",valfilein);
				return 8;
	}
	

        pdnmsModelstruc=MakeAModelstruc(bsp);
	if (pdnmsModelstruc==NULL) {
		ErrPostEx(SEV_ERROR,7,1,"Unable to convert Biostruc to Modelstruc");
		return 7;
	}

 	pmsdRoot=(PMSD)(pdnmsModelstruc->data.ptrvalue);
 

	if (pmsdRoot==NULL) {
		ErrPostEx(SEV_ERROR,6,1,"Internal error...");
		return 6;
	}
       if (pmsdRoot->pcPDBName != NULL) printf("PDB code for the structure is: %s\n",pmsdRoot->pcPDBName);


	/* get number of models */
        if (Rargs[1].intvalue == 0) { /* find the best all atom model */
	
	pdnmlModel=pmsdRoot->pdnmlModels;
	while (pdnmlModel) {
		pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
		if (pmldThis->iType==Model_type_ncbi_all_atom || pmldThis->iType==Model_type_pdb_model) {
			modelnum = pdnmlModel->choice;
		}
		pdnmlModel=pdnmlModel->next;
	}


/*	printf("Set modelnum = %d\n",modelnum); */

		if (!modelnum) {
			ErrPostEx(SEV_ERROR,5,1,"No all-atom models found");
			return 5;
			}	
	}
	if ( (Rargs[1].intvalue > 0) && (Rargs[1].intvalue < 9999) ) { /* find the user-specified model */
		pmldHere=GetModelN(pmsdRoot,(Int2) Rargs[1].intvalue);
		if (pmldHere==NULL){
			ErrPostEx(SEV_ERROR,3,1,"No model %d found in structure!!",Rargs[1].intvalue);
			return 3;
		}
		modelnum = (Int2) Rargs[1].intvalue; /* validated requested model actually in there.. */
	}



/* write out the pdb structure as a pdb file to inspect as test */
       if ((pFile=FileOpen(pdbfileout,"w")) !=NULL) {    
        
		if (Rargs[1].intvalue == 9999) {
			/* write out all models */
			iTest = WritePDBAllModel(pmsdRoot->pdnmsLink, pFile); 
		} 
		else {  /* write out the specified model */
 			iTest = WritePDBOneModel(pmsdRoot->pdnmsLink,  pFile,  modelnum);
		}         
        	fflush(pFile);
         	FileClose(pFile);
       }
	else {
	ErrPostEx(SEV_ERROR,2,1,"Cannot write to output PDB file %s",pdbfileout);
			return 2;
}
      


/* Shut Down MMDB-API */
/* All Modelstrucs remaining are freed in CloseMMDB-API() */
	CloseMMDBAPI();	
 	return 0;
}







/*  
$Log: 

*/

