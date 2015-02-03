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
Args Rargs[3] = {{"Input Asn.1 3D structure as .prt .val or .cn3  File Name.",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
               {"Sequences only (t/F):","FALSE",NULL,NULL,TRUE,'s',ARG_BOOLEAN,0.0,0,NULL},
               {"Include Ramachandran Angles (t/F):","FALSE",NULL,NULL,TRUE,'r',ARG_BOOLEAN,0.0,0,NULL}
};
                       



void WriteStdoutRamaOneModel(PMSD pmsdRoot, Int2 ModelNum ){

ValNodePtr vnpRama = NULL;
PRS prsHead = NULL;
PRS prsHere = NULL;
float phi,psi;
int num;
Int4 linelen;
Char aa[3]; 
Char ctemp[30];
Char chain[30];
PMGD pmgdAA = NULL;
CharPtr NCBIstdaaUC = "-ABCDEFGHIKLMNPQRSTVWXYZU*";

vnpRama=ConvertNode((PFB)pmsdRoot,AM_MGD);
if (vnpRama==NULL) return;
prsHead=Rama(vnpRama,ModelNum);
prsHere=prsHead;
while(prsHere) {
    phi = 0.0;
	psi = 0.0;
	phi=(float)prsHere->Phi;
	psi=(float)prsHere->Psi;
	pmgdAA = (PMGD) prsHere->pfbThis;
	num = (int) (pmgdAA->pdnmgLink->choice);
	StringCpy(ctemp, StringChr(NCBIstdaaUC,pmgdAA->pcIUPAC[0]));
        StringCpy(chain, ParentMolName((PFB) pmgdAA));
	aa[0] = ctemp[0];
	aa[1] = '\0';
	printf("%s, %s, %d, %f, %f\n", chain, aa , (int) num, phi,psi);
  	prsHere=prsHere->next;
	}
 freeRS(prsHead);
 return;
}




Int2 Main()
{
    PMSD  pmsdRoot = NULL;

    static char valfilein[PATH_MAX];
    Int2    iDotLen = 0;
	AsnIoPtr aip=NULL;
	NcbiMimeAsn1Ptr nmap=NULL;
	BiostrucPtr bsp=NULL;
	BiostrucSeqPtr bssp=NULL;
	PDNMS pdnmsModelstruc = NULL;
	ErrSev esMsg,esLog,esFatal;
	PDNML pdnmlModel;
	PMLD pmldThis;
	Boolean isMime = FALSE;
	Int4 iModels = 0, iNCBImodel = 0, iPDBmodel = 0, iAscii = 0;


    
/* Initialize MMDB-API */
        ErrSetLogfile("error_strSummary.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);

	if (!GetArgs("strSummary  reports on what is inside .val or .cn3 or .prt file to stdout.\nCaution - will not work on Biounits downloaded from NCBI, use Asymmetric Units\n",3,Rargs))
		return 1;


        /* we use exact filenames and extensions to avoid confusion */

	
        StringCpy(valfilein, Rargs[0].strvalue);
        if (FileLength(valfilein) == 0) 
	  {
			ErrPostEx(SEV_FATAL,13,1,"Unable to find input file %s",Rargs[0].strvalue);
			return 13;	
       	 }

	if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,12,1,"Unable to open MMDBAPI, check for missing bstdt.val dictionary file.");
		return 12;
	}
	



/* load an ASN.1 file *.val or *.c3d or *.prt write out Biostruc only */

	aip=AsnIoOpen(valfilein,"rb");
	if (aip==NULL) {
		
		ErrPostEx(SEV_FATAL,11,1,"Unable open ASN.1 stream in file %s",valfilein);
			return 11;
	}

/* first try biostruc load */
/* 	printf("try biostruc\n");  */
	esMsg=ErrGetMessageLevel();
	esLog=ErrGetLogLevel();
	esFatal=ErrGetFatalLevel();
        ErrSetMessageLevel(SEV_MAX);
        ErrSetLogLevel(SEV_MAX);
  	ErrSetFatalLevel(SEV_MAX);
	
 	bsp=BiostrucAsnRead(aip,NULL);
	AsnIoClose(aip);

  	ErrSetMessageLevel(esMsg);
  	ErrSetLogLevel(esLog);
  	ErrSetFatalLevel(esFatal);

	if (bsp == NULL) { /* try ascii */
		aip = AsnIoOpen(valfilein,"r");
		if (aip==NULL) {
			ErrPostEx(SEV_FATAL,10,1,"Unable open binary or ascii ASN.1 from file %s",valfilein);
			return 10;
		}

		esMsg=ErrGetMessageLevel();
		esLog=ErrGetLogLevel();
		esFatal=ErrGetFatalLevel();
     	 	ErrSetMessageLevel(SEV_MAX);
	        ErrSetLogLevel(SEV_MAX);
	  	ErrSetFatalLevel(SEV_MAX);
	
	 	bsp=BiostrucAsnRead(aip,NULL);
		AsnIoClose(aip);

	  	ErrSetMessageLevel(esMsg);
  		ErrSetLogLevel(esLog);
  		ErrSetFatalLevel(esFatal);
		if (bsp) iAscii = 1;


	}



	if (bsp==NULL) {
/*		printf("try mime binary\n"); */
		/* then try NCBIMime load */
		aip=NULL;
		aip=AsnIoOpen(valfilein,"rb");
  		ErrSetMessageLevel(SEV_MAX);
  		ErrSetLogLevel(SEV_MAX);
  		ErrSetFatalLevel(SEV_MAX);

		nmap=NcbiMimeAsn1AsnRead(aip,NULL);

  		ErrSetMessageLevel(esMsg);
  		ErrSetLogLevel(esLog);
  		ErrSetFatalLevel(esFatal);
		AsnIoClose(aip);

		if (nmap == NULL) {
/*			printf("try mime ascii\n"); */
			aip=AsnIoOpen(valfilein,"r");
  			ErrSetMessageLevel(SEV_MAX);
  			ErrSetLogLevel(SEV_MAX);
  			ErrSetFatalLevel(SEV_MAX);

			nmap=NcbiMimeAsn1AsnRead(aip,NULL);

  			ErrSetMessageLevel(esMsg);
  			ErrSetLogLevel(esLog);
  			ErrSetFatalLevel(esFatal);
			AsnIoClose(aip);
			if (nmap) iAscii = 1;
		}

        	if (nmap!=NULL) {
		/* got an NCBI mime */
			if (nmap->choice!=NcbiMimeAsn1_strucseq && nmap->choice!=NcbiMimeAsn1_strucseqs) {
				/* wrong MIME type */
				nmap=NcbiMimeAsn1Free(nmap);
		    		ErrPostEx(SEV_ERROR,9,1,"MIME-type wrapper is not a strucseq or strucseqs - no Biostruc to report in:  %s",valfilein);
				return 9;
			}
				isMime = TRUE;
			/* unwrap the mime and leave the bsp for PDB conversion */
				bssp=(BiostrucSeqPtr)(nmap->data.ptrvalue);
				bsp=bssp->structure;
                		bssp->structure = NULL;
/* may want to report on contents of wrapper ... */
				nmap=NcbiMimeAsn1Free(nmap); /* discard the wrapper */
				nmap = NULL;
				bssp = NULL;
			}
		}

	if (bsp == NULL) {
	    	ErrPostEx(SEV_ERROR,8,1,"No Biostruc in files to report in: %s",valfilein);
				return 8;
	}
	

        pdnmsModelstruc=MakeAModelstruc(bsp);
	if (pdnmsModelstruc==NULL) {
		ErrPostEx(SEV_ERROR,7,1,"Unable to convert Biostruc to Modelstruc");
		return 7;
	}

 	pmsdRoot=(PMSD)(pdnmsModelstruc->data.ptrvalue);

	if (pmsdRoot==NULL) {
		ErrPostEx(SEV_ERROR,6,1,"Internal error.");
		return 6;
	}

/* OK, Ready to Start Reporting */

	printf("Filename: %s\n", valfilein);
if (Rargs[1].intvalue == FALSE) {

        if (iAscii) printf("ASCII Asn.1 (.prt) ");
	else printf("Binary Asn.1 (.val) ");
	if (isMime) printf("Entry Point: NCBIMime\n");
	else printf("Entry Point: Biostruc\n");


        if (pmsdRoot->pcPDBName != NULL) printf("PDB code: %s\n",pmsdRoot->pcPDBName);
        printf("MMDB Id: %d\n", (int) pmsdRoot->iMMDBid);

	printf("Name: %s\n",GetStrucStrings(pdnmsModelstruc, LONG_NAME));
	printf("PDB Class: %s\n",GetStrucStrings(pdnmsModelstruc, PDB_CLASS));
	printf("PDB Source: %s\n",GetStrucStrings(pdnmsModelstruc, PDB_SOURCE)); 
	
	
	pdnmlModel=pmsdRoot->pdnmlModels;
	while (pdnmlModel) {
		iModels++;
		pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
		if (pmldThis->iType==Model_type_ncbi_all_atom) {
		 iNCBImodel = pdnmlModel->choice;
		}
		if (pmldThis->iType==Model_type_pdb_model) {
		 iPDBmodel = pdnmlModel->choice;
		}
		pdnmlModel=pdnmlModel->next;
	}

	printf("Number of Models: %d\n",iModels);
	if(iNCBImodel) printf("NCBI All Atom Model: %d\n",iNCBImodel); 
	if (iPDBmodel) printf("PDB All Atom Model: %d\n",iPDBmodel);
	if (pmsdRoot->pDictLocal) printf("Has Local Chemical Graph Dictionary\n");
        printf("Number of Molecules: %d\n", (int) pmsdRoot->iMolCount); /* number of molecules */
  	printf("Number of Objects: %d\n",(int)  pmsdRoot->iObjCount); /* number of objects */
   	printf("Number of Density Sets: %d\n",(int) pmsdRoot->iDensCount); /* number of densities */
        printf("Number of Inter-Molecule Bonds: %d\n", (int) pmsdRoot->iIMBCount); 
	if (Rargs[2].intvalue == TRUE) {
       		printf("\n\nRAMACHANDRAN Angles for all Proteins:\n");
       		if (iNCBImodel) {
       			printf("NCBI All Atom Coordinate Model Number %d\n\nChain, AA, Num, Phi, Psi\n",iNCBImodel); 
			WriteStdoutRamaOneModel(pmsdRoot, iNCBImodel );
		}
       		else if (iPDBmodel) {
       			printf("PDB All Atom Coordinate Model Number %d:\n\nChain, AA, Num, Phi, Psi\n",iPDBmodel); 
			WriteStdoutRamaOneModel(pmsdRoot, iPDBmodel );
		}
	}
}
        printf("\n\n");
/* Output FASTA for each chain and dump het list - oldie but goodie code ! */
	WriteFASTASeqHet(pdnmsModelstruc, stdout);
      

/* Shut Down MMDB-API */
/* All Modelstrucs remaining are freed in CloseMMDB-API() */
	CloseMMDBAPI();	
 	return 0;
}

