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

#include "valmerge.h"
#include "objmime.h"
#include "mkbioseq.h"

#define NUMARGS 7

Args Rargs[NUMARGS] = {{"Input val file 1",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
                        {"Input val file 2",NULL,NULL,NULL,FALSE,'g',ARG_FILE_IN,0.0,0,NULL},
			{"Output val file",NULL,NULL,NULL,TRUE,'o',ARG_FILE_OUT,0.0,0,NULL},
                       {"Model Number for file 1",NULL,"1","9999",FALSE,'m',ARG_INT,0.0,0,NULL},
			{"Model Number for file 2",NULL,"1","9999",FALSE,'n',ARG_INT,0.0,0,NULL},
                     {"Chain to be copied to val file 1 (default: all chains)","-",NULL,NULL,TRUE,'c',ARG_STRING,0.0,0,NULL},
			{"Output ascii  files?","FALSE",NULL,NULL,TRUE,'p',ARG_BOOLEAN,0.0,0,NULL}};


Int2 Main()
{

	
	
/*Find the chains for the molecule first*/
PMSD pmsd1,pmsd2;
Int2 Model1,Model2;
PDNMM pdnmm1,pdnmm2;
PMMD pmmd1,pmmd2;
char chainFrom;
static char valfileout[256];
Int2    iTest = 0;
FILE    *pFile = NULL;
PDNMS pdnms1=NULL;
BiostrucPtr pbsBS = NULL;
SeqEntryPtr sep = NULL;
PRGD prgdDictionary = NULL;
NcbiMimeAsn1Ptr nmap;
BiostrucSeqPtr bsqp;
AsnIoPtr aip = NULL;
/*For copying local residue-graph*/
ResidueGraphPtr pDict1=NULL,pDict2=NULL,pDict1Head=NULL,pDict1Last=NULL;

	
	
/* Initialize MMDB-API */
        ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);


        if (!GetArgs("ValMerge 1.0: Merge One or All chains from two mmdb val files to form a single 3D structure.\n",NUMARGS,Rargs))
                return 1;

	if(Rargs[2].strvalue==NULL) {
        	StringCpy(valfileout,"valmerge_out.val");
        }	
        else {
        	StringCpy(valfileout,Rargs[2].strvalue);
        }  

        if (!OpenMMDBAPI(0,NULL)) {
                ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI");
                return 2;
        }
/*Only works for single character chain names for now*/
	chainFrom=toupper((Rargs[5].strvalue)[0]);

/*load the two structures*/

    	Model1=Rargs[3].intvalue;
    	Model2=Rargs[4].intvalue;	
  	pmsd1=LoadABiostruc(Rargs[0].strvalue,0,0,&Model1);
  	pmsd2=LoadABiostruc(Rargs[1].strvalue,0,0,&Model2);
 	if (pmsd1==NULL || pmsd2==NULL) {
		/* error occurred */
		return 3;
	}
	
/*hack for 1QC6 structure, just drag the pDictLocal of pmsd2 to pmsd1*/	
/*Only work for pmsd1 that don't have local dict and all pmsd2s have same localdict*/


	pdnmm2=pmsd2->pdnmmHead;
	if (pdnmm2==NULL){
		printf("Unable to get list of molecules for structure 2\n");
		return 3;
	}
	
	printf("The next chain of Structure 1 to append is %s\n", NextUniqueChainName(pmsd1));
	while(pdnmm2!=NULL){
		
		pmmd2=(PMMD)pdnmm2->data.ptrvalue;
		/*printf("the chain now is %s\n",pmmd2->pcMolName);
		printf("the chain wanted is %c\n",chainFrom);*/
		if(chainFrom=='-' || toupper(pmmd2->pcMolName[0]) == chainFrom){
				if(!CopyBiomolecule(pmsd1,pdnmm2)){
				printf("Copy Biomolecule Failed.\n");
				return 3;
			}
		}
		pdnmm2=pdnmm2->next;
	 /*   printf("FYI - The next chain of Structure 1 to append would be %s\n", NextUniqueChainName(pmsd1)); */
	}

/*	pdnmm1=pmsd1->pdnmmHead;
	printf("New list of molecules including HETs :\n");
	while(pdnmm1 != NULL){
		printf("%d:",pdnmm1->choice);		
		pmmd1=(PMMD)pdnmm1->data.ptrvalue;
		if (pmmd1==NULL) break;
		printf("%s, ", pmmd1->pcMolName); 
		pdnmm1=pdnmm1->next;
	}
	printf("\n");

*/ 
 
/* Only work for pmsd1 that don't have local dict and all pmsd2s have same localdict */
/* TraDES structures should never have local dictionaries - all in bstdt.val */
/* so if Structure 1 is a PDB file, this will keep its local dictionary */
/* if structure 2 is a PDB file, it will move it to the final structure */
/* this may fail if both files have local dictionaries */


	/* LOCAL DICtioARY 2 id references must be renumbered to the end of the id rage of 1, and its references updated */



/*88888888888888888888888888888888888888888*/	

/*Adding local residue graphs*/
   if (pmsd1->pDictLocal==NULL &&  pmsd2->pDictLocal) {
	pDict2=pmsd2->pDictLocal;
	while(pDict2!=NULL){
		pDict1=(ResidueGraphPtr)AsnIoMemCopy(pDict2, (AsnReadFunc)ResidueGraphAsnRead, (AsnWriteFunc)ResidueGraphAsnWrite);
		if(pDict1==NULL) printf("pDict1 is NULL!!\n");
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
	if(pDict1Head==NULL) printf("pDict1Head is NULL!!\n");
	pmsd1->pDictLocal=pDict1Head;
	if(pmsd1->pDictLocal) printf("Local Residue Library added\n");
   }

/*88888888888888888888888888888888888888888*/	
	
	/* WRITE ASN.1 Biostruc (pre-mime) TO A FILE  */
	
	WriteAsnOneModel(pmsd1->pdnmsLink,Model1,valfileout,SAVE_BINARY);	
	
	/* at this stage the file is complete but non MIMEd */
	/* Note.. cannot write to a bytestore unless go in and change mmdbapi4.c */

    printf("Biostruc merged, now converting to Mime-Biostrucseq\n");
	
/* reread in ASN.1 biostruc - non mimed - from file - no longer MMDBAPI bound */
	aip=AsnIoOpen(valfileout,"rb");
	pbsBS = BiostrucAsnRead(aip,NULL);
	AsnIoClose(aip);

	if (pbsBS == NULL) printf("Failed to read in new Biostruc - pre-MIME stage\n");
/* create MIME Biostruct using VAST service MakeBioseqs code */
			
    prgdDictionary = GetPRGDDictionary(); /* retrieves the bstdt dictionary pointer from MMDBAPI */
	if (prgdDictionary == NULL) { printf("\nNo dictionary ??\n"); goto donehere; } else ; /* printf("\nGot dictionary\n"); */
       

/* call VAST code to make complete set of new sequences with annotation */
	sep = (SeqEntryPtr) MakeBioseqs (pbsBS, prgdDictionary);
    if (sep == NULL) { printf("No Seq Entry from Biostruc... "); goto donehere; }
		
/* Allocate and attach sequence and structure objects to nmap */
	
	bsqp=BiostrucSeqNew();
	bsqp->structure=pbsBS; /* attach Biostruc */
	ValNodeLink(&(bsqp->sequences),sep); /* attach sequences */
	nmap=(NcbiMimeAsn1Ptr)ValNodeNew(NULL);
	nmap->choice=NcbiMimeAsn1_strucseq;
	nmap->next=NULL;
	nmap->data.ptrvalue=(VoidPtr)bsqp;

	
	/* Write out MIME Biostruct */

	if (Rargs[6].intvalue==TRUE)
	  aip=AsnIoOpen(valfileout,SAVE_ASCII_STRING);
	else aip=AsnIoOpen(valfileout,SAVE_BINARY_STRING);
	NcbiMimeAsn1AsnWrite(nmap,aip,NULL);
	AsnIoClose(aip);
	printf("Completed\n");	
	
    CloseMMDBAPI();
    /* NcbiMimeAsn1Free(nmap); */
	/* there is a doubly freed pointer here dangling still, need to resolve it */

	return TRUE;
 
donehere:


CloseMMDBAPI();

/* was double-deleting some objects */	
	
return FALSE;
}

