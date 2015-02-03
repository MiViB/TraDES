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
#include "val2prt.h"

/* Global Variables */
Args Rargs[NUMARGS] = {{"Input *.c3d Filename from VAST alignment (No Extension)",NULL,NULL,
		NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL}};


/* Program loads an ncbi-mime alignment from pairwise VAST */
/* Ignores the master */
/* Rotates each slave into the master coordinate space by applying transform */
/* Writes out the biostruc for each rotated slave - valmerge can be used to merge chains */



Int2 Main()
{
	Char fnamin[255];
	Char fnamout[255];
	PMSD pmsd;
	Int2 model=1;
	AsnIoPtr aip = NULL;
	AsnTypePtr atp = NULL, atpAliStruc = NULL, atpSlave = NULL;
	NcbiMimeAsn1Ptr pmime = NULL;
	

	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
	if (!GetArgs("Alignstruc 2 Biostruc Program",NUMARGS,Rargs))
		return 1;
	OpenMMDBAPI(0,NULL);
	/* load an ASN.1 Biostruc */
	StringCpy(fnamin,Rargs[0].strvalue);
	StringCpy(fnamout,Rargs[0].strvalue);
	StringCat(fnamin,".c3d");
	StringCat(fnamout,".prt");
	if ((aip = AsnIoOpen(fnamin,"r") == NULL)){
		printf("File %s could not be opened\n",fnamin);
		return 1;
	}
	pmime = NcbiMimeAsn1AsnRead(aip, atp);
	atpAliStruc = AsnFind("Biostruc-align");
	if (atpAliStruc == NULL) {
		printf("No Biostruc-align in file %s \n",fnamin);
		return 1;
	}
/*	atpSlave = AsnFind("");
	pFeatureSet = AsnFind("");
	pTransform = AsnFind("");

	pmsd=LoadABiostruc(fnamin,0,2,&model);
*/
	if (pmsd==NULL) {
		ErrPostEx(SEV_FATAL,3,1,"Unable to fetch Biostruc");
		return 3;
	}
/*	WriteAsnAllModel(pmsd->pdnmsLink,fnamout,SAVE_ASCII); */
	NcbiMimeAsn1AsnFree(pmime);
	CloseMMDBAPI();
 	return TRUE;
}


