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


#include <hfprogs.h>
#include <mmdbapi.h>

/* store Ramachandran plot information */
/*
typedef struct nlm_ramastruct {
	struct nlm_ramastruct PNTR next;
	PFB pfbThis;
	FloatHi Phi;
	FloatHi Psi;
	FloatHi Omega;
	FloatHi Mag;
	FloatHi Chi1;
	FloatHi Chi2;
	FloatHi Chi3;
	FloatHi Chi4;
} RS, *PRS, **PPRS;

*/

const int BS_STARTSIZE = 50000;


/* Global Variables */
Args Rargs[3] = {
{"Input VAL File Name (NO EXTENSION).",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
{"Foldtraj Range Start Number (optinal)","0","1","9999999",TRUE,'s',ARG_INT,0.0,0,NULL},
{"Foldtraj Range (optional)","0","1","50000",TRUE,'r',ARG_INT,0.0,0,NULL}};
             
ByteStorePtr bspA;
ByteStorePtr bspC;
ByteStorePtr bspD;
ByteStorePtr bspE;
ByteStorePtr bspF;
ByteStorePtr bspG;
ByteStorePtr bspH;
ByteStorePtr bspI;
ByteStorePtr bspK;
ByteStorePtr bspL;
ByteStorePtr bspM;
ByteStorePtr bspN;
ByteStorePtr bspP;
ByteStorePtr bspQ;
ByteStorePtr bspR;
ByteStorePtr bspS;
ByteStorePtr bspT;
ByteStorePtr bspV;
ByteStorePtr bspW;
ByteStorePtr bspY;
ByteStorePtr bspX;


void WriteRamaOneModel(PMSD pmsdRoot, Int2 ModelNum, Int4 structureNo ){

ValNodePtr vnpRama = NULL;
PRS prsHead = NULL;
PRS prsHere = NULL;
float phi,psi;
int num;
Int4 linelen;
Char aa[3]; 
Char ctemp[30];
Char line[255];
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
	aa[0] = ctemp[0];
	aa[1] = '\0';
	sprintf(line, "%d, %d, %f, %f\n",(int)structureNo, (int) num, phi,psi);
    linelen = StringLen(line);
	/* store into the appropriate ByteStore */
switch (aa[0]) {
case 'A':
   BSWrite(bspA, line, linelen);
   break;
case 'C':
   BSWrite(bspC, line, linelen);
   break;
case 'D':
   BSWrite(bspD, line, linelen);
   break;
case 'E':
   BSWrite(bspE, line, linelen);
   break;
case 'F':
   BSWrite(bspF, line, linelen);
   break;
case 'G':
   BSWrite(bspG, line, linelen);
   break;
case 'H':
   BSWrite(bspH, line, linelen);
   break;
case 'I':
   BSWrite(bspI, line, linelen);
   break;
case 'K':
   BSWrite(bspK, line, linelen);
   break;
case 'L':
   BSWrite(bspL, line, linelen);
   break;
case 'M':
   BSWrite(bspM, line, linelen);
   break;
case 'N':
   BSWrite(bspN, line, linelen);
   break;
case 'P':
   BSWrite(bspP, line, linelen);
   break;
case 'Q':
   BSWrite(bspQ, line, linelen);
   break;
case 'R':
   BSWrite(bspR, line, linelen);
   break;
case 'S':
   BSWrite(bspS, line, linelen);
   break;
case 'T':
   BSWrite(bspT, line, linelen);
   break;
case 'V':
   BSWrite(bspV, line, linelen);
   break;
case 'W':
   BSWrite(bspW, line, linelen);
   break;
case 'Y':
   BSWrite(bspY, line, linelen);
   break;
case 'X':
   BSWrite(bspX, line, linelen);
   break;
case '-': 
case 'Z':
case 'B':
case 'U':
case '*':
default:
break;

}

	prsHere=prsHere->next;
	}
 freeRS(prsHead);
 return;
}



void ByteStoreFileOut(ByteStorePtr bspHere, CharPtr filename){

CharPtr pcTempBuf;
Int4 numbytes;
FILE *fp;

 /*  write arbitrary bytestore to a file in one chunk of memory, if you got it*/ 
    numbytes=BSLen(bspHere);
	if (numbytes==0)
	    return;
	pcTempBuf=(CharPtr)MemNew((size_t)numbytes);
    BSSeek(bspHere,0L,SEEK_SET);
    BSRead(bspHere,pcTempBuf,numbytes);
	if ((fp=FileOpen(filename,"a"))==NULL) {
        ErrPostEx(SEV_ERROR,1,8,"Unable to open file %s for writing!",filename);
        return;
    }	
	if ( (FileWrite(pcTempBuf,sizeof(Char),numbytes,fp)) != numbytes ){
	    FileClose(fp);
	    pcTempBuf=MemFree(pcTempBuf);
        ErrPostEx(SEV_ERROR,1,9,"Unable to write to log file %s",filename);
        return;
    }
    FileClose(fp);
	pcTempBuf=MemFree(pcTempBuf);
}




Int2 Main()
{
    PMSD  pmsdRoot = NULL;
    Int2  ModelNum = 1;

/* Variables for clipping status and pdb output file name*/
    static char fileout[PATH_MAX];
	static char filein[PATH_MAX];
	static char txtnum[10];
    Int2    iTest = 0;
    FILE    *pFile = NULL;
	int     start;
	int     finish;
	int     thisone;
	int     writecount;

bspA = BSNew(BS_STARTSIZE);
bspC = BSNew(BS_STARTSIZE);
bspD = BSNew(BS_STARTSIZE);
bspE = BSNew(BS_STARTSIZE);
bspF = BSNew(BS_STARTSIZE);
bspG = BSNew(BS_STARTSIZE);
bspH = BSNew(BS_STARTSIZE);
bspI = BSNew(BS_STARTSIZE);
bspK = BSNew(BS_STARTSIZE);
bspL = BSNew(BS_STARTSIZE);
bspM = BSNew(BS_STARTSIZE);
bspN = BSNew(BS_STARTSIZE);
bspP = BSNew(BS_STARTSIZE);
bspQ = BSNew(BS_STARTSIZE);
bspR = BSNew(BS_STARTSIZE);
bspS = BSNew(BS_STARTSIZE);
bspT = BSNew(BS_STARTSIZE);
bspV = BSNew(BS_STARTSIZE);
bspW = BSNew(BS_STARTSIZE);
bspX = BSNew(BS_STARTSIZE);
bspY = BSNew(BS_STARTSIZE);


    
/* Initialize MMDB-API */
        ErrSetLogfile("error_ramangL.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);

    if (!GetArgs("ramangL: \nReport Phi, Psi from TraDES *.val file\n Creates up to 20 *.csv files for R plotting.\n",3,Rargs))
		return 1;
    StringCpy(fileout, "__");
	StringCat(fileout,Rargs[0].strvalue);
	if(StringLen(fileout) == 2) {
	  ErrPostEx(SEV_FATAL,2,1,"No filename provided");
	  return 3;
	}
    StringCat(fileout,".csv");
	start = Rargs[1].intvalue;
	writecount = 1;
	finish = (Rargs[1].intvalue + Rargs[2].intvalue);
	thisone = start;
	printf("Processing files from %0.7d - %0.7d \n",start,finish-1);
	
	if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI, check for bstdt.val");
		return 4;
	}
	
	/* loopy bit here */
	
	while (thisone < finish) {
	StringCpy(filein,Rargs[0].strvalue);
/* if this sprintf_s gives you trouble comipling, google for a new library - it is the windows safe form of sprintf 
please let me know if any platform builds hang on this one  - cwvh feb 2012 */
 /*   sprintf(txtnum,"%07d", thisone);   this will work if you need it */
	sprintf(txtnum, "%07d", thisone);
    StringCat(filein, txtnum);
	StringCat(filein, ".val");

	if (FileLength(filein) != 0) {
	pmsdRoot=LoadABiostruc(filein,FALSE,0,&ModelNum);
	if (pmsdRoot==NULL) {
		ErrPostEx(SEV_INFO,2,1,"Unable to load %s, skipping",filein);
		ClearStructures();
		continue;
	}    
        else
		ErrPostEx(SEV_INFO,2,1,"File %s does not exist, skipping",filein);



	printf("%s Loaded\n",filein);
	if (!(pmsdRoot->bWhat & AM_PROT)) {
			ErrPostEx(SEV_INFO,6,1,"No protein in structure, skipping");
			ClearStructures();
			continue;
		}

    WriteRamaOneModel(pmsdRoot, ModelNum, thisone);
	}
    ClearStructures();
    thisone++;
	writecount++;
	if (writecount == 1000) {
	printf("Writing batch of 1000\n");
    writecount = 1;

fileout[0] = 'A';
ByteStoreFileOut( bspA, fileout);
BSFree(bspA );
fileout[0] = 'C';
ByteStoreFileOut( bspC, fileout);
BSFree(bspC );
fileout[0] = 'D';
ByteStoreFileOut( bspD, fileout);
BSFree(bspD );
fileout[0] = 'E';
ByteStoreFileOut( bspE, fileout);
BSFree(bspE );
fileout[0] = 'F';
ByteStoreFileOut( bspF, fileout);
BSFree(bspF );
fileout[0] = 'G';
ByteStoreFileOut( bspG, fileout);
BSFree(bspG );
fileout[0] = 'H';
ByteStoreFileOut( bspH, fileout);
BSFree(bspH );
fileout[0] = 'I';
ByteStoreFileOut( bspI, fileout);
BSFree(bspI );
fileout[0] = 'K';
ByteStoreFileOut( bspK, fileout);
BSFree(bspK );
fileout[0] = 'L';
ByteStoreFileOut( bspL, fileout);
BSFree(bspL );
fileout[0] = 'M';
ByteStoreFileOut( bspM, fileout);
BSFree(bspM );
fileout[0] = 'N';
ByteStoreFileOut( bspN, fileout);
BSFree(bspN );
fileout[0] = 'P';
ByteStoreFileOut( bspP, fileout);
BSFree(bspP );
fileout[0] = 'Q';
ByteStoreFileOut( bspQ, fileout);
BSFree(bspQ );
fileout[0] = 'R';
ByteStoreFileOut( bspR, fileout);
BSFree(bspR );
fileout[0] = 'S';
ByteStoreFileOut( bspS, fileout);
BSFree(bspS );
fileout[0] = 'T';
ByteStoreFileOut( bspT, fileout);
BSFree(bspT );
fileout[0] = 'V';
ByteStoreFileOut( bspV, fileout);
BSFree(bspV );
fileout[0] = 'W';
ByteStoreFileOut( bspW, fileout);
BSFree(bspW );
fileout[0] = 'Y';
ByteStoreFileOut( bspY, fileout);
BSFree(bspY );
fileout[0] = 'X';
ByteStoreFileOut( bspX, fileout);
BSFree(bspX);

bspA = BSNew(BS_STARTSIZE);
bspC = BSNew(BS_STARTSIZE);
bspD = BSNew(BS_STARTSIZE);
bspE = BSNew(BS_STARTSIZE);
bspF = BSNew(BS_STARTSIZE);
bspG = BSNew(BS_STARTSIZE);
bspH = BSNew(BS_STARTSIZE);
bspI = BSNew(BS_STARTSIZE);
bspK = BSNew(BS_STARTSIZE);
bspL = BSNew(BS_STARTSIZE);
bspM = BSNew(BS_STARTSIZE);
bspN = BSNew(BS_STARTSIZE);
bspP = BSNew(BS_STARTSIZE);
bspQ = BSNew(BS_STARTSIZE);
bspR = BSNew(BS_STARTSIZE);
bspS = BSNew(BS_STARTSIZE);
bspT = BSNew(BS_STARTSIZE);
bspV = BSNew(BS_STARTSIZE);
bspW = BSNew(BS_STARTSIZE);
bspX = BSNew(BS_STARTSIZE);
bspY = BSNew(BS_STARTSIZE);


	}



	}
	
/* Shut Down MMDB-API */
/* All Modelstrucs remaining are freed in CloseMMDB-API() */
	CloseMMDBAPI();	

/* Write out all the files remaining */

fileout[0] = 'A';
ByteStoreFileOut( bspA, fileout);
BSFree(bspA );
fileout[0] = 'C';
ByteStoreFileOut( bspC, fileout);
BSFree(bspC );
fileout[0] = 'D';
ByteStoreFileOut( bspD, fileout);
BSFree(bspD );
fileout[0] = 'E';
ByteStoreFileOut( bspE, fileout);
BSFree(bspE );
fileout[0] = 'F';
ByteStoreFileOut( bspF, fileout);
BSFree(bspF );
fileout[0] = 'G';
ByteStoreFileOut( bspG, fileout);
BSFree(bspG );
fileout[0] = 'H';
ByteStoreFileOut( bspH, fileout);
BSFree(bspH );
fileout[0] = 'I';
ByteStoreFileOut( bspI, fileout);
BSFree(bspI );
fileout[0] = 'K';
ByteStoreFileOut( bspK, fileout);
BSFree(bspK );
fileout[0] = 'L';
ByteStoreFileOut( bspL, fileout);
BSFree(bspL );
fileout[0] = 'M';
ByteStoreFileOut( bspM, fileout);
BSFree(bspM );
fileout[0] = 'N';
ByteStoreFileOut( bspN, fileout);
BSFree(bspN );
fileout[0] = 'P';
ByteStoreFileOut( bspP, fileout);
BSFree(bspP );
fileout[0] = 'Q';
ByteStoreFileOut( bspQ, fileout);
BSFree(bspQ );
fileout[0] = 'R';
ByteStoreFileOut( bspR, fileout);
BSFree(bspR );
fileout[0] = 'S';
ByteStoreFileOut( bspS, fileout);
BSFree(bspS );
fileout[0] = 'T';
ByteStoreFileOut( bspT, fileout);
BSFree(bspT );
fileout[0] = 'V';
ByteStoreFileOut( bspV, fileout);
BSFree(bspV );
fileout[0] = 'W';
ByteStoreFileOut( bspW, fileout);
BSFree(bspW );
fileout[0] = 'Y';
ByteStoreFileOut( bspY, fileout);
BSFree(bspY );
fileout[0] = 'X';
ByteStoreFileOut( bspX, fileout);
BSFree(bspX);

    printf("Done!\n");
	fileout[0]='*';
	printf("Output file names are %s where * is amino acid ACDEFGHIKLMNPQRSTVWY or X\nReminder: Make Blank *.csv files for missing aa before using R fns",fileout);

 	return TRUE;
}
