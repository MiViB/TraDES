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

/* Program to merge one chain from one val file into another */
/* CWVH 2008 initial version */
/* Added mime-biostruc bioseq generator system from VAST server code March 2011 */
/* more cleanup in 2012 - the mkbioseq code has a mess of warnings that still need to be resolved */
/* This is the program formerly known as valmerge - has rigorous logging now due to some flaky libraries */

/* Both the nuschncpy and vast library with MakeBioseqs still 
  need to be purged of gcc -Wall compiler warnings and checked with Valgrind */

#include <hfprogs.h>
#include "nuschncpy.h"
#include "objmime.h"

#ifdef OS_MSWIN
#include "mkbioseq.h"
#endif

#ifdef OS_UNIX
/* gcc compiler compaining about the required function prototype in mkbioseq.h ... not sure why... so just force it here */
extern SeqEntryPtr MakeBioseqs (BiostrucPtr bsp, BiostrucResidueGraphSetPtr stdDictionary);
#endif

#define NUMARGS 8

Args Rargs[NUMARGS] = {
/*0*/  {"Input 3d Asn.1 file 1 (destination - To)",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
/*1*/  {"Input 3d Asn.1 file 2 (source - From)",NULL,NULL,NULL,FALSE,'g',ARG_FILE_IN,0.0,0,NULL},
/*2*/  {"Output file Name (extension ignored, default='strMerge_out.val')",NULL,NULL,NULL,TRUE,'o',ARG_FILE_OUT,0.0,0,NULL},
/*3*/  {"Model Number for file 1 (default - all-atom model)","0","0","9999",TRUE,'m',ARG_INT,0.0,0,NULL},
/*4*/  {"Model Number for file 2 (default - all-atom model)","0","0","9999",TRUE,'n',ARG_INT,0.0,0,NULL},
/*5*/  {"Chain to be copied from file 2 to file 1 (default: all chains)","-",NULL,NULL,TRUE,'c',ARG_STRING,0.0,0,NULL},
/*6*/  {"Output Asn.1 Binary (default 0=.val) (1=.cn3) or Asn.1 ASCII (2=.prt)","0","0","2",TRUE,'a',ARG_INT,0.0,0,NULL},
/*7*/  {"Verbose mode (T/f)","TRUE",NULL,NULL,TRUE,'q',ARG_BOOLEAN,0.0,0,NULL} };





Int2 Main()
{


	PMSD pmsd1,pmsd2;
	Int2 Model1,Model2;
	PDNMM pdnmm1,pdnmm2;
	PMMD pmmd1,pmmd2;
	char chainFrom;
	static char valfileout[256];
	BiostrucPtr bsp = NULL;
	SeqEntryPtr sep = NULL;
	PRGD prgdDictionary = NULL;
	NcbiMimeAsn1Ptr nmap;
	BiostrucSeqPtr bssp;
	AsnIoPtr aip = NULL;
	/*For copying local residue-graph*/
	ResidueGraphPtr pDict1=NULL,pDict2=NULL,pDict1Head=NULL,pDict1Last=NULL;
	Int2 iDotLen = 0;
	Int2 iAscii = 0;
	Int2 iBestModel1 = 0, iBestModel2 = 0;
	PDNML pdnmlModel = NULL;
	PDNMS pdnmsModelstruc = NULL;
	ErrSev esMsg,esLog,esFatal;
	PMSD pmsdRoot = NULL;
	PMLD pmldThis = NULL;

	
	

        ErrSetLogfile("error_strMerge.log", ELOG_APPEND|ELOG_BANNER);
        ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
        if (!GetArgs("strMerge   Merge One or All chains from one Asn.1 3D files into another Asn.1 file.\n",NUMARGS,Rargs))
                return 1;

/* set the output filename */
	if(Rargs[2].strvalue==NULL) {
        	StringCpy(valfileout,"strMerge_out.val");
        }	
        else {

      		iDotLen = 0;
	   	iDotLen = (Int2) StringLen(Rargs[2].strvalue);
           	if (iDotLen > 5) {
			if (Rargs[2].strvalue[iDotLen - 4] == '.') {
	        	Rargs[2].strvalue[iDotLen - 4] = '\0';      
			}
           	}
/* append the user-chosen extension */
		StringCpy(valfileout,Rargs[2].strvalue);
		if (Rargs[6].intvalue == 0) 
	        	StringCat(valfileout,".val");
		if (Rargs[6].intvalue == 1) 
	        	StringCat(valfileout,".cn3");
		if (Rargs[6].intvalue == 2) {
	        	StringCat(valfileout,".prt");
			iAscii = 1;
		}
	}

     
/* Initialize MMDB-API */


        if (!OpenMMDBAPI(0,NULL)) {
                ErrPostEx(SEV_FATAL,2,2,"Unable to open MMDBAPI");
                return 2;
        }

/*Only works for single character chain names PDB style - NCBI Biounits are not supported yet.*/
	chainFrom=toupper((Rargs[5].strvalue)[0]);


if (Rargs[7].intvalue == TRUE)
		printf("Loading Structures..\n" );

/* see if files exist */
	if (FileLength(Rargs[0].strvalue) == 0) 
	  {
		ErrPostEx(SEV_FATAL,2,3,"Unable to find 3d Asn.1 file 1: %s",Rargs[0].strvalue);
		return 3;
       	 }
	if (FileLength(Rargs[1].strvalue) == 0) 
	  {
		ErrPostEx(SEV_FATAL,2,4,"Unable to find 3d Asn.1 file 2: %s",Rargs[0].strvalue);
		return 4;
       	 }



/*Load the two structures with only one model in memory... */

/* Load in the destination-TO structure */
    	if (Rargs[3].intvalue == 0) { /* If model not set, find the best models - load in the entire structure the hard way and go looking, then reload just that model ... */

		if (Rargs[7].intvalue == TRUE)
			printf("Looking for the best model in (destination - To) file 1: %s  ..\n",Rargs[0].strvalue);
		aip=AsnIoOpen(Rargs[0].strvalue,"rb");
		if (aip==NULL) {
			ErrPostEx(SEV_FATAL,3,1,"Unable open binary or ascii ASN.1 stream in file 1 %s",Rargs[0].strvalue);
			return 301;
		}
		/* try binary Biostruc */
		esMsg=ErrGetMessageLevel(); esLog=ErrGetLogLevel(); esFatal=ErrGetFatalLevel();
        	ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
	 	bsp=BiostrucAsnRead(aip,NULL);
		AsnIoClose(aip);
  		ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);
		if (bsp == NULL) { /* try ascii Biostruc */
			aip = AsnIoOpen(Rargs[0].strvalue,"r");
			if (aip==NULL) {
				ErrPostEx(SEV_FATAL,3,2,"Unable open binary or ascii ASN.1 stream in file %s",Rargs[0].strvalue);
			return 302;
			}
    		 	ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
		 	bsp=BiostrucAsnRead(aip,NULL);
			AsnIoClose(aip);
		  	ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);	
		}
		if (bsp==NULL) {
			/*  try binary NCBIMime load */
			aip=NULL;
			aip=AsnIoOpen(Rargs[0].strvalue,"rb");
  			ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);	
			nmap=NcbiMimeAsn1AsnRead(aip,NULL);
  			ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);
			AsnIoClose(aip);
			if (nmap == NULL) {
				/* try ascii NCBImime load */
				aip=AsnIoOpen(Rargs[0].strvalue,"r");
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
			    		ErrPostEx(SEV_ERROR,3,3,"MIME-type wrapper is not a strucseq or strucseqs - no Biostruc to convert to PDB %s",Rargs[0].strvalue);
					return 303;
				}
				bssp=(BiostrucSeqPtr)(nmap->data.ptrvalue);
				bsp=bssp->structure;
	       	         	bssp->structure = NULL;
				nmap=NcbiMimeAsn1Free(nmap); /* discard the wrapper */
				nmap = NULL;
				bssp = NULL;
			}
		}
		if (bsp == NULL) {
	    		ErrPostEx(SEV_ERROR,3,4,"No Biostruc in file1 (destination - To)  %s",Rargs[0].strvalue);
				return 304;
		}	
      		pdnmsModelstruc=MakeAModelstruc(bsp);
		if (pdnmsModelstruc==NULL) {
			ErrPostEx(SEV_ERROR,3,5,"Unable to convert Biostruc to Modelstruc %s",Rargs[0].strvalue);
			return 305;
		}
 		pmsdRoot=(PMSD)(pdnmsModelstruc->data.ptrvalue);
		pdnmlModel = NULL;
		pdnmlModel=pmsdRoot->pdnmlModels;
		while (pdnmlModel) {
			pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
			if (pmldThis->iType==Model_type_ncbi_all_atom || pmldThis->iType==Model_type_pdb_model) {
				iBestModel1 = pdnmlModel->choice;
			}
			pdnmlModel=pdnmlModel->next;
		}
		ClearStructures();
		pmsdRoot = NULL;
		pdnmsModelstruc = NULL;
		/* now we got the model to load - reload just with that model */
		Model1=iBestModel1;
	}
	else  /* user-specified model, easy */
	   Model1=(Int2) Rargs[3].intvalue; 
	

/* Load in the source-From structure */
	if (Rargs[4].intvalue ==0) { /* If model not set, find the best models - load in the entire structure the hard way and go looking, then reload just that model ... */
		if (Rargs[7].intvalue == TRUE)
			printf("Looking for the best model in (source - From) file 1: %s  ..\n",Rargs[1].strvalue);
		aip=AsnIoOpen(Rargs[1].strvalue,"rb");
		if (aip==NULL) {
		
			ErrPostEx(SEV_FATAL,4,1,"Unable open binary or ascii ASN.1 stream in file %s",Rargs[1].strvalue);
			return 401;
		}
		/* try binary Biostruc */
		esMsg=ErrGetMessageLevel(); esLog=ErrGetLogLevel(); esFatal=ErrGetFatalLevel();
        	ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
	 	bsp=BiostrucAsnRead(aip,NULL);
		AsnIoClose(aip);
  		ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);
		if (bsp == NULL) { /* try ascii Biostruc */
			aip = AsnIoOpen(Rargs[1].strvalue,"r");
			if (aip==NULL) {
				ErrPostEx(SEV_FATAL,4,2,"Unable open binary or ascii ASN.1 stream in file %s",Rargs[1].strvalue);
			return 402;
			}
    		 	ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);
		 	bsp=BiostrucAsnRead(aip,NULL);
			AsnIoClose(aip);
		  	ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);	
		}
		if (bsp==NULL) {
			/*  try binary NCBIMime load */
			aip=NULL;
			aip=AsnIoOpen(Rargs[1].strvalue,"rb");
  			ErrSetMessageLevel(SEV_MAX); ErrSetLogLevel(SEV_MAX); ErrSetFatalLevel(SEV_MAX);	
			nmap=NcbiMimeAsn1AsnRead(aip,NULL);
  			ErrSetMessageLevel(esMsg); ErrSetLogLevel(esLog); ErrSetFatalLevel(esFatal);
			AsnIoClose(aip);
			if (nmap == NULL) {
				/* try ascii NCBImime load */
				aip=AsnIoOpen(Rargs[1].strvalue,"r");
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
			    		ErrPostEx(SEV_ERROR,4,3,"MIME-type wrapper is not a strucseq or strucseqs - no Biostruc to convert to PDB %s",Rargs[1].strvalue);
					return 403;
				}
				bssp=(BiostrucSeqPtr)(nmap->data.ptrvalue);
				bsp=bssp->structure;
	       	         	bssp->structure = NULL;
				nmap=NcbiMimeAsn1Free(nmap); /* discard the wrapper */
				nmap = NULL;
				bssp = NULL;
			}
		}
		if (bsp == NULL) {
	    		ErrPostEx(SEV_ERROR,4,4,"No Biostruc in files 2 (source - From) %s",Rargs[1].strvalue);
				return 404;
		}	
      		pdnmsModelstruc=MakeAModelstruc(bsp);
		if (pdnmsModelstruc==NULL) {
			ErrPostEx(SEV_ERROR,4,5,"Unable to convert Biostruc to Modelstruc %s",Rargs[1].strvalue);
			return 405;
		}
 		pmsdRoot=(PMSD)(pdnmsModelstruc->data.ptrvalue);
 		pdnmlModel = NULL;
		pdnmlModel=pmsdRoot->pdnmlModels;
		while (pdnmlModel) {
			pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
			if (pmldThis->iType==Model_type_ncbi_all_atom || pmldThis->iType==Model_type_pdb_model) {
				iBestModel2 = pdnmlModel->choice;
			}
			pdnmlModel=pdnmlModel->next;
		}
		ClearStructures();
		pmsdRoot = NULL;
		pdnmsModelstruc = NULL;
		/* now we got the model to load - reload just with that model */
		Model2=iBestModel2; 
	}
	else /* user-specified model, easy */
    		Model2=(Int2) Rargs[4].intvalue; 
  	
/* Wow, that was a lot of work, but now that Model1 and Model2 are set, load just ONLY those models. Defaults are wonderful :) */

	if (Rargs[7].intvalue == TRUE)
			printf("Loading Model %d (destination - To) file 1: %s  ..\n",Model1,Rargs[0].strvalue);
  	pmsd1=LoadABiostruc(Rargs[0].strvalue,0,0,&Model1);
 	if (pmsd1==NULL) {
			ErrPostEx(SEV_FATAL,2,5,"Unable to load 3d structure Model %d from Asn.1 file 1 (destination - To): %s",Model1,Rargs[0].strvalue);
			return 5;
	} 
	if (Rargs[7].intvalue == TRUE)
			printf("Loading Model %d (source - From) file 2: %s  ..\n",Model2,Rargs[1].strvalue);
	pmsd2=LoadABiostruc(Rargs[1].strvalue,0,0,&Model2);
 	if (pmsd2==NULL) {
		ErrPostEx(SEV_FATAL,2,6,"Unable to load 3d structure Model %d from Asn.1 file 2 (source - From): %s",Model2,Rargs[1].strvalue);
		return 6;
	}

/*

This has to be fixed here and in solvate ! Or do it in a separate desolvate utility first 
	 
	if (Rargs[7].intvalue == TRUE)
			printf("Removing solvent from (destination - To) file 1: %s  ..\n",Rargs[0].strvalue);
*/

	if (Rargs[7].intvalue == TRUE) {
 		pdnmm1=pmsd1->pdnmmHead;
		printf("List of file 1 (destination - To) molecules Before copy :\n");
		while(pdnmm1 != NULL){	
			pmmd1=(PMMD)pdnmm1->data.ptrvalue;
			if (pmmd1==NULL) break;
			if ((IsNAorProtein((PFB)pmmd1))==TRUE) { 
				printf("%d:",pdnmm1->choice);
				printf("%s, ", pmmd1->pcMolName);
			} 
			pdnmm1=pdnmm1->next;
		}
		printf("\n");
	} 

/*Find the chains for the molecule first*/	

	pdnmm2=pmsd2->pdnmmHead;
	if (pdnmm2==NULL){
		ErrPostEx(SEV_FATAL,2,7,"Apparently no list of molecules in (source - From) Asn.1 file 2: %s",Rargs[0].strvalue);
		return 7;
	}


	if (Rargs[7].intvalue == TRUE)
		printf("The next available chain of 3d file 1 (destination - To) for appending is [%s].\n", NextUniqueChainName(pmsd1));
	while(pdnmm2!=NULL){
		pmmd2=(PMMD)pdnmm2->data.ptrvalue;
		if(chainFrom=='-' || toupper(pmmd2->pcMolName[0]) == chainFrom){
				/* skip solvent */
				if ((IsSolvent((PFB)pmmd2))!=TRUE) {				
					if (Rargs[7].intvalue == TRUE)
						printf("Found chain %s (or HET) in file 2 (source - From) %s matching requested %c.\nCopying Molecule to file %s ..\n",
					pmmd2->pcMolName,Rargs[1].strvalue,chainFrom,Rargs[0].strvalue);
				
					if(!CopyBiomolecule(pmsd1,pdnmm2)){
						ErrPostEx(SEV_FATAL,2,8,"Copy Molecule Internal Code Failure. \n");
					return 8;
				}
			}
		}
		pdnmm2=pdnmm2->next;
	}

	if (Rargs[7].intvalue == TRUE) {
 		pdnmm1=pmsd1->pdnmmHead;
		printf("After copy - list of MERGED molecules to be written  :\n");
		while(pdnmm1 != NULL){	
			pmmd1=(PMMD)pdnmm1->data.ptrvalue;
			if (pmmd1==NULL) break;
			if ((IsNAorProtein((PFB)pmmd1))==TRUE) {
				printf("%d:",pdnmm1->choice);	
				printf("%s, ", pmmd1->pcMolName); 
			}
			pdnmm1=pdnmm1->next;
		}
		printf("\n");
	} 
 


/* For standard proteins, skipping the dictionary copy works sometimes - fails other times */

if (pmsd1->pDictLocal && pmsd2->pDictLocal) {
	  ErrPostEx(SEV_WARNING,3,1,"Warning - both structures have local residue graph dictionaries - assuming copied (source - From) chains(s) use global dictionary.\n");
	  if (Rargs[7].intvalue == TRUE) {
		  printf("Warning - both structures have local residue graph dictionaries\n");
		  printf("          Output keeps only the local dictionary in (destination - To) file 1 %s\n",Rargs[0].strvalue);
		  printf("          If the copied chain(s) originating in (source - From) file 2 %s\n", Rargs[1].strvalue);
		  printf("          only refer to standard residues from global dictionary bstdt.val\n");
		  printf("          this will be ok.\n");
	}
}


/* Local Dictionary Moving */
/* Currently only one of the structures can have a local dictionary */

/*Simple copy of local residue graph from 2 to 1 */
/* only works when file 1 has no local dictionary */
   if (pmsd1->pDictLocal==NULL &&  pmsd2->pDictLocal) {
	pDict2=pmsd2->pDictLocal;
	while(pDict2!=NULL){
		pDict1=(ResidueGraphPtr)AsnIoMemCopy(pDict2, (AsnReadFunc)ResidueGraphAsnRead, (AsnWriteFunc)ResidueGraphAsnWrite);
		if(pDict1==NULL) {
			ErrPostEx(SEV_FATAL,2,9,"Copy Local Residue Graph Dictionary Memory Alloction Failure\n");
			return 9;
		}
		if(pDict1Head==NULL){
			pDict1Head=pDict1;
			pDict1Last=pDict1;
		}
		else{
			pDict1Last->next=pDict1;
			pDict1Last=pDict1;
		}
		pDict2=pDict2->next;
	}
	pmsd1->pDictLocal=pDict1Head;
	if((pmsd1->pDictLocal) && (Rargs[7].intvalue == TRUE))
		printf("Local Residue Graph Dictionary copied from %s to %s\n",Rargs[1].strvalue,Rargs[0].strvalue);
   }


/*8888888888888888---Automating Local Dictionary Copies - To Do--8888888888888888888888888*/	

/* Only works for pmsd1 that don't have local dict and all pmsd2s have same localdict */
/* TraDES structures should never have local dictionaries - all in bstdt.val */
/* so if Structure 1 is a PDB file, this will keep its local dictionary */
/* IF structure 2 is from PDB file, it will move its local dictionary to the final structure */
/* It MAY fail if both files have local dictionaries AND the copied chains refer to it in structure 2 */


/* For the case when both structures have local dictionaries - complexity ensues..*/

/* Must append Local Dictionary 2 to Local Dictionary 1 and update all the dictionary numbering in sequence.
  Must update the new molecule's references to match the new numbering.

Some of this was done for the solvation code (i.e. append a dictionary and find the last number at the end...).
Must traverse the new molecule and update the numbering..

*/


/*888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/	
	
/* WRITE ASN.1 Biostruc (pre-mime) TO A FILE  */
	
    if (Rargs[7].intvalue == TRUE)
	printf("Writing preliminary merged Biostruc to %s..\n",valfileout);

	if (iAscii==1)  {
		if (WriteAsnOneModel(pmsd1->pdnmsLink,Model1,valfileout,SAVE_ASCII) == FALSE) {
			ErrPostEx(SEV_FATAL,2,10,"Error Writing ASCII Asn.1 file %s\n",valfileout);
			return 10;
		}
        }
        else	
		if (WriteAsnOneModel(pmsd1->pdnmsLink,Model1,valfileout,SAVE_BINARY) == FALSE) {
			ErrPostEx(SEV_FATAL,2,11,"Error Writing BINARY Asn.1 file %s\n",valfileout);
			return 11;
	}	
	
/* at this stage the file is complete but non MIMEd */
/* Note.. cannot write to a bytestore unless go in and change mmdbapi4.c */

    ClearStructures();
    if (Rargs[7].intvalue == TRUE)
	printf("Biostruc merged and written to %s\n Validating Asn.1 input parse..\n",valfileout);
	
/* re-read in ASN.1 biostruc - non mimed - from file - no longer MMDBAPI bound */
/* mucho validation... and if MakeBioseqs dies, you still have a Biostruc :) */

	bsp = NULL;
	if (iAscii==1)  {
		aip=AsnIoOpen(valfileout,"r");
	}
	else {        
		aip=AsnIoOpen(valfileout,"rb");
	}
	bsp = BiostrucAsnRead(aip,NULL);
	AsnIoClose(aip);
	if (bsp == NULL) {
		if (Rargs[7].intvalue == TRUE)
			printf("Asn.1 file %s will not load back in, terminating\n",valfileout);
		ErrPostEx(SEV_FATAL,2,12,"Merged Asn.1 Biostruc Failed Parse Validation %s\n",valfileout);
			return 12;
	}

	if (Rargs[7].intvalue == TRUE)
		printf("Biostruc Asn.1 File %s Validated...\n",valfileout);

/* create MIME Biostruct using VAST service MakeBioseqs code :o  */
	
/* shape this bsp into an MIME-wrapped Asn.1 structure nmap  */
			
	prgdDictionary = GetPRGDDictionary(); /* retrieves the bstdt dictionary pointer from MMDBAPI */
	if (prgdDictionary == NULL) { 
		ErrPostEx(SEV_ERROR,2,13,"Lost the dictionary internally...");
			return 13;
	}


/* call VAST code to make complete set of new sequences (protein only???) with some  annotation */
	sep = (SeqEntryPtr) MakeBioseqs (bsp, prgdDictionary);
	if (sep == NULL) { 
		ErrPostEx(SEV_ERROR,2,14,"MakeBioSeqs finds no sequences in Biostruc - Cannot add MIME-wrapper to Biostruc: %s",valfileout);
		return 14;
	}
 
		
/* Allocate and attach sequence and structure objects to nmap */
	
	bssp=BiostrucSeqNew();
	bssp->structure=bsp; /* attach Biostruc */
	ValNodeLink(&(bssp->sequences),sep); /* attach sequences */
	nmap=(NcbiMimeAsn1Ptr)ValNodeNew(NULL);
	nmap->choice=NcbiMimeAsn1_strucseq;
	nmap->next=NULL;
	nmap->data.ptrvalue=(VoidPtr)bssp;

		
	if (nmap != NULL) { /* write out Mime-type as read in with no changes  */
		if (Rargs[7].intvalue == TRUE)
			printf("Opening output file for MIME-Biostruc %s..\n",valfileout);
		if (iAscii == 1)
			aip=AsnIoOpen(valfileout,"w");
		else
			aip=AsnIoOpen(valfileout,"wb");
		if (aip==NULL) {
			ErrPostEx(SEV_FATAL,2,15,"Unable open file for writing MIME-wrapped Asn.1 %s",valfileout);
			return 15;
		}
		if (Rargs[7].intvalue == TRUE)
			printf("Writing MIME-wrapped Asn.1 file %s ..\n",valfileout);
                if (NcbiMimeAsn1AsnWrite(nmap,aip,NULL) == FALSE) {
			ErrPostEx(SEV_FATAL,2,16,"NcbiMimeAsn1AsnWrite failed to write Asn.1 file %s",valfileout);
			return 16;
		}
		AsnIoClose(aip);
		if (Rargs[7].intvalue == TRUE)
			printf("Asn.1 file %s completed.\n",valfileout);
   		/* NcbiMimeAsn1Free(nmap); */
     	}
	
   	CloseMMDBAPI(); 
	return 0;
}

