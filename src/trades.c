/*Portions Copyright (c) 2007-2012 Hogue Laboratory, National University of Singapore
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

/* Protein Folding program using TRADES algorithm
This program takes an arbitrary primary amino acid sequence as
its input and attempts to determine the folded, 3-D structure
of the protein in aqueous solution.

Programmed by: Christopher Hogue and Howard Feldman
Designed by: Christopher Hogue and Howard Feldman



v 2.0 
Started May 15 2012

*/

#include "trades.h"
#include <tslri.h>

#ifdef OS_UNIX
#include <sys/utsname.h>
#endif

#define ZHANG_WINDOWSIZE 1   /* Reset to 1 CWVH 2012 - Small Peptide scores are off. Entropy Calculation Calibration */ 
#define BRYANT_WINDOWSIZE 3  /* set to 3 to match crease entropy window size so same counts are used */
#define VSCORE_WINDOWSIZE 1


/* Global Variables */


Args StrucTrjArgs[20] = {
/*0*/	{"Input Trajectory Distribution File (extension ignored, looks for .trj)\nREQUIRED:",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
/*1*/	{"Output filename to store structure (extension ignored, truncated)\n      (default = same as input trajectory distribution):",NULL,NULL,NULL,TRUE,'o',ARG_FILE_OUT,0.0,0,NULL},
/*2*/	{"Reproducible Random Seed? (0 = Truly Random):","0","-2147483647","2147483647",TRUE,'n',ARG_INT,0.0,0,NULL},
/*3*/	{"Keep input file when done?:","TRUE",NULL,NULL,TRUE,'d',ARG_BOOLEAN,0.0,0,NULL},
/*4*/	{"Output PDB files?:","FALSE",NULL,NULL,TRUE,'p',ARG_BOOLEAN,0.0,0,NULL},
/*5*/	{"Output ASN.1 files?:","TRUE",NULL,NULL,TRUE,'a',ARG_BOOLEAN,0.0,0,NULL},
/*6*/	{"Number of structures to build:","1","1","300000",TRUE,'b',ARG_INT,0.0,0,NULL},
/*7*/	{"Quiet operation?:","TRUE",NULL,NULL,TRUE,'q',ARG_BOOLEAN,0.0,0,NULL},
/*8*/   {"Comparison Asn.1 structure file (extension ignored, looks for: .val .cn3 .prt):",NULL,NULL,NULL,TRUE,'c',ARG_FILE_IN,0.0,0,NULL},
/*9*/	{"Comparison structure Chain\n      (default=first molecule in file):","-",NULL,NULL,TRUE,'h',ARG_STRING,0.0,0,NULL},
/*10*/	{"Experiment (default=none):",NULL,NULL,NULL,TRUE,'x',ARG_STRING,0.0,0,NULL},
/*11*/	{"Structure Numbering start at: ","1","1","10000000",TRUE,'s',ARG_INT,0.0,0,NULL},
/*12*/	{"Build only fragment covered range:","FALSE",NULL,NULL,TRUE,'F',ARG_BOOLEAN,0.0,0,NULL},	 
/*13*/	{"Write 20 Ramachandran Angle Files:","FALSE",NULL,NULL,TRUE,'r',ARG_BOOLEAN,0.0,0,NULL},
/*14*/	{"ASN.1 MIME BINARY .val Output (FALSE = ASCII .prt):","TRUE",NULL,NULL,TRUE,'v',ARG_BOOLEAN,0.0,0,NULL},
/*15*/	{"No. of H20 layers for Ellipsoid and PCA 3D axis Transform\n      (Default = 1, 0 for no solvation):","1","0","10",TRUE,'l',ARG_INT,0.0,0,NULL},
/*16*/	{"Keep H20 in Structure Output:","FALSE",NULL,NULL,TRUE,'w',ARG_BOOLEAN,0.0,0,NULL},
/*17*/  {"Stream Progress Update to stdout","FALSE",NULL,NULL,TRUE,'t',ARG_BOOLEAN,0.0,0,NULL},
/*18*/  {"Skip PCA rotation after Ellipsoid calculations:","FALSE",NULL,NULL,TRUE,'k',ARG_BOOLEAN,0.0,0,NULL},
/*19*/  {"Allow Atom Bounciness to Increase when stuck: (Only set true with input from str2trj!)","FALSE",NULL,NULL,TRUE,'z',ARG_BOOLEAN,0.0,0,NULL}
  	};

static char *ComplexityDescriptions[] = {"Virtual Bond Model","All Atom Model","Up to 5 Models","Up to 10 Models","All Models","Cn3D Subset"};
static char *ColorDescriptions[] = {"","Molecule Number","Secondary Structure","Thermal Factor","Element"};
static Int4 incompletecount=0;
Boolean timetoquit=FALSE;

const int BS_STARTSIZE = 50000;

ValNodePtr pvnSolventList = NULL;
ValNodePtr pvnCur = NULL;

/* CWVH Solvation Globals */
static int iAtomCount = 0;
static int iGraphCount = 0;
static int iWaterCount = 0;
static DValNodePtr pdnmmEnd = NULL;
static Int4 iLocalDict = 0;
static ResidueGraphPtr prgWat = NULL; /* pointer to the Water Chemical Graph */

/* Bytestores for Ramachandran Output */

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





void LIBCALLBACK DoCount(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr)
{
    if (pfbThis == NULL) return;  
    if (pfbThis->bMe != (Byte) AM_MAD) return;
 	if (GetAtomLocs((PMAD) pfbThis, iModel) != NULL)
			iAtomCount++;
}



void LIBCALLBACK DoCountGraphs(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr)
{
    if (pfbThis == NULL) return;  
    if (pfbThis->bMe != (Byte) AM_MGD) return;
			iGraphCount++;
}




void RemoveWaterDictionary(PMSD pmsdThis){
ResidueGraphPtr prgHere;
ResidueGraphPtr prgPrev;

prgHere = pmsdThis->pDictLocal;

	if (prgHere == NULL) {
		return;
	}
 /* simple case only one entry in local dictionary */
  if ((prgHere->next == NULL) && (prgHere->id == iLocalDict)) {
		  ResidueGraphFree(prgHere);
		  pmsdThis->pDictLocal = NULL;
		  prgWat = NULL;
          iLocalDict = 0;
		  return;
   }
  /* local dictionary already in front, seek to the end */
  while (prgHere->next != NULL) {
	        prgPrev = prgHere;
			prgHere = prgHere->next;	
  }
  if (prgHere->id == iLocalDict) {
       ResidueGraphFree(prgHere);
	   prgPrev->next = NULL;
  }
	/* Reset Library Globals */
  prgWat = NULL;
  iLocalDict = 0;
  
 
	
}



Int4 SolvateIt(PMAD pmadAtom, ValNodePtr pvnAtomList, DValNodePtr pdnListPmmds, Int2 iModel, PWS pwsHead, Boolean bSkipInternal )
{

 
 ValNodePtr pvnHere = NULL;
 ValNodePtr pvnStrings = NULL;
 ValNodePtr pvnDescr1 = NULL;
 ValNodePtr pvnDescr2 = NULL;
 vec vCenter, vTest, vDist;
 PMMD pmmdHere = NULL;
 CharPtr pcResidue = NULL;
 Char cNearRes[200];
 PMSD pmsdThis = NULL;
 Int2 iNeighbors = 0;
 Int2 iNearWaters = 0;
 char *pj;
 int i;


 Int2 iWater = 0;

 	  FloatLoPtr pflvWat = NULL;
	  PALD paldWat= NULL;
	  ValNodePtr pvnalWat = NULL;
	  PMAD pmadWat = NULL;
	  ValNodePtr pvnmaWat = NULL;
      PMGD pmgdWat = NULL;
	  DValNodePtr pdnmgWat = NULL;
	  PMMD pmmdWat = NULL;
	  DValNodePtr pdnmmWat = NULL;

 PALD pald1 = NULL, pald2 = NULL;
 FloatLo dist = 0.0;


 /* 10 water coordinates surrounding the target atom */
/*  some h-bonded together hence the spacing */
/* centred at 0 0 0 (empty - where target atom goes) with a maximum distance of 3.65, min of 2.71 A to centre */ 
vec water1[] = {
	 {2.087,   0.222,  -2.180}, 
	 {-0.189,   1.429,  -2.946}, 
	 {-1.615,  -1.035,  -1.964}, 
	 {-0.156,   2.813,  -0.272}, 
	 {2.348,  -1.947,  -0.717}, 
	 {-2.935,   1.523,   1.022}, 
	 {0.137,   0.292,   2.690}, 
	 {-1.239,  -2.175,   2.659}, 
	 {2.864,   1.510,   0.862}, 
	 {-0.267,  -2.729,  -0.012}}; 

/* the second set are rotated 45 degrees about each axis so that they might pack better on alternate calls on neighboring structure atoms */

vec water2[] = {
	{-2.180,   2.087,   0.222},
	{-2.946,  -0.189,   1.429},
	{-1.964,  -1.615,  -1.035},
	{-0.272,  -0.156,   2.813},
	{-0.717,   2.348,  -1.947},
	{1.022,  -2.935,   1.523},
	{2.690,   0.137,   0.292},
	{2.659,  -1.239,  -2.175},
	{0.862,   2.864,   1.510},
	{-0.012,  -0.267,  -2.729}};


static char keep[11];

vec Twater[] = {
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0},
	 {0.0,   0.0,  0.0}};





if (pmadAtom == NULL) return 0;
if (pvnAtomList == NULL) return 0;
if (pdnListPmmds == NULL) return 0;
if (pwsHead == NULL) return 0;

pvnHere = pvnAtomList;

cNearRes[0] = '\0';
i = 0;
pcResidue = ParentGraphIUPAC1((PFB) pmadAtom);
pmmdHere= GetParentMol((PFB)pmadAtom);
if ((pcResidue != NULL) && (pmmdHere->bWhat == (Byte) AM_PROT)) {
	  i++;
	  cNearRes[i-1] = pcResidue[0]; /* The IUPAC code of the amino acid residue */
	  cNearRes[i] = '\0';
}

if (bSkipInternal) { /* skip internal solvent molecules */
 while (pvnHere != NULL) {
	iNeighbors++;
	pald2 = (PALD)(pvnHere->data.ptrvalue);
	pmmdHere = GetParentMol((PFB) pald2); 
	pcResidue = ParentGraphIUPAC1((PFB) pald2);
	if ((pcResidue != NULL) && (pmmdHere->bWhat == (Byte) AM_PROT) && (pald2->pfbParent->bWhat != (Byte) AM_BACKBONE)) {
	  i++;
	  cNearRes[i-1] = pcResidue[0]; /* The IUPAC code of the amino acid residue  - sidechain atoms only */
	  cNearRes[i] = '\0'; /* list contains codes of nearby sidechain residues */ 
	}
    if (pmmdHere->bWhat & (Byte) AM_SOL)
	   iNearWaters++;
	pvnHere=pvnHere->next;
 }
 pmmdHere = NULL;
 pald2 = NULL;


i=0;
pj = strpbrk(cNearRes , "FMILVWA" );
while (pj != NULL) {
	i++;
	pj++;
	pj = strpbrk(pj, "FMILVWA");
}

 /* This neighbor list - in a dense region will contain about 60 atoms  and about 30 along a smooth surface */
/*printf("[ %d NN : %d nHOH : %d of FMILVW : %s]", (int) iNeighbors, (int) iNearWaters, i, cNearRes);  */
if ((iNeighbors >= 48) && (iNearWaters <= 12)) return 0;   
/* Ad-hoc rules to discourage buried atom solvation - tested on 1RRO & 1CC3 - this just skips dense regions of the protein   */

/* also ad-hoc as some internal core residues have very low density of neighbors, so based on the count of hydrophobic residues in low-density regions
, we skip these too */
/* value of 32 is close to the minimum coordination number for a packed atom in a protein, according to the paper of the Zhang-Delisi Potential */
/* when using a 6 angstrom cutoff distance */

if ((iNeighbors >=32) && (iNearWaters <= 10) && (i >= 17)) return 0; 

/* note that this also tends not to solvate narrow crevices on the protein surface as well. */


} /* bSkipInternal */

pmmdHere = (PMMD) pdnListPmmds->data.ptrvalue;
if (pdnmmEnd == NULL) {
  pdnmmEnd = pdnListPmmds;
  while (pdnmmEnd->next != NULL) {
		  pdnmmEnd = pdnmmEnd->next;
	  } /* get to end of molecule list */
}
    


/* pdmmEnd is a global  - keep track of the end of the molecule list so we don't traverse to the ends for each of thousands of atoms added */
pmsdThis = ToMSDParent((PFB)pmmdHere);

if (iAtomCount == 0) {  /* sets the global count of the number of atoms in the structure for adding water serial numbers */
  TraverseAtoms(pmsdThis->pdnmsLink,iModel,0, NULL, (pNodeFunc) (DoCount));
/*  printf("In solvate it with %d\n", (int) iAtomCount); */
}

if (iGraphCount == 0) { /* likewise, this counts the number of graphs in the structure for numbering the water mol (as graphs) */
	TraverseGraphs(pmsdThis->pdnmsLink,iModel,0,NULL, (pNodeFunc) (DoCountGraphs));
}

/* printf("Atoms: %d Graphs: %d \n",(int) iAtomCount, (int) iGraphCount); */

/* translate W1-W10 to PMAD coordinate  use a matrix of vecs */
pald1 = NULL;
pald1 = GetAtomLocs(pmadAtom, iModel);
if (pald1 == NULL) return 0;
vCenter[0] = AtomLocX(pald1);
vCenter[1] = AtomLocY(pald1);
vCenter[2] = AtomLocZ(pald1);
/* printf("Atom at %8.3f %8.3f %8.3f\n",(float) vCenter[0], (float)vCenter[1],(float)vCenter[2]); */
for (i = 0; i<10; i++) {
/*	printf("%d-",i); */
  Twater[i][0] = 0.0;
  Twater[i][1] = 0.0;
  Twater[i][2] = 0.0;
  keep[i] = '1';
  keep[i+1] = '\n';
  if ((iAtomCount % 2) == 0)
    VecAdd(Twater[i],vCenter,water2[i]); /* even calls */
  else
    VecAdd(Twater[i],vCenter,water1[i]); /* odd calls */
 /* 	printf("water translated %8.3f %8.3f %8.3f\n",(float) Twater[i][0], (float)Twater[i][1],(float)Twater[i][2]); */

}

pvnHere = pvnAtomList;
while (pvnHere != NULL) {
/* test distances to the pvnAtomList */
	 pald2 = (PALD)(pvnHere->data.ptrvalue);
     vTest[0] = AtomLocX(pald2);
	 vTest[1] = AtomLocY(pald2);
	 vTest[2] = AtomLocZ(pald2);
	 for (i = 0; i<10; i++) {
/* printf("In Distance Test %d %8.3f %8.3f %8.3f\n",i, (float) vTest[0],(float) vTest[1],(float) vTest[2]); */
		 if (keep[i] == '1') {
			VecSub(vDist,vTest,Twater[i]);  
			dist=getMag(vDist);
			if (dist < 2.65) {  
/* less than 2.65 to any atom - we mark to remove the water sometimes this may be a wee close, but never to the atom being solvated
because all the atom locations in the array are fixed , so it behaves randomly bouncy */
				keep[i] = '0';
			/*	printf("Dist crash %8.3f\n", (float) dist); */
			}
		 } /* skip waters already removed */
	 } /* for loop over waters */
	pvnHere = pvnHere -> next;
} /* while loop over structure atoms */
/* only keep waters have keep marks of 1 */
iWater = 0;
for (i = 0; i<10; i++) {
	if (keep[i] != '0')  {
      iWater++;
	  iWaterCount++;


/* Output for distance distribution histogram */
/*
pvnHere = pvnAtomList;
while (pvnHere != NULL) { */
/* Print out distances from the saved water to the pvnAtomList */
/*	 pald2 = (PALD)(pvnHere->data.ptrvalue);
     vTest[0] = AtomLocX(pald2);
	 vTest[1] = AtomLocY(pald2);
	 vTest[2] = AtomLocZ(pald2);
	 VecSub(vDist,vTest,Twater[i]);  
	 dist=getMag(vDist);
     printf("%8.3f\n", (float) dist);  
     pvnHere = pvnHere -> next;
} 
*/

/* printf("Adding water to Modelstruc\n"); */
     /* keep all the ones that fit */
     /* these need to be made into Het nodes and added to the Structure */
	 /* bottom up - allocate pvnma, pald, pmad, pmgd, pmmd */
     /* attach to end of the pdnListPmmds */
	 


	  pflvWat = FLVector(0, (Int4) 4); /* x, y, z, occupancy and b factor */
      paldWat = NewALD();
	  pvnalWat = ValNodeNew(NULL);
	  pmadWat = NewMAD();
	  pvnmaWat = ValNodeNew(NULL);
	  pmgdWat = NewMGD();
	  pdnmgWat = DValNodeNew(NULL);
	  pmmdWat = NewMMD();
	  pdnmmWat = DValNodeNew(NULL);
	  pvnDescr1 = ValNodeNew(NULL);
	  pvnDescr2 = ValNodeNew(NULL);

	  if ((pflvWat == NULL) ||
		  (paldWat == NULL) ||
		  (pvnalWat == NULL) ||
		  (pmadWat == NULL) ||
		  (pvnmaWat == NULL) ||
		  (pmgdWat == NULL) ||
		  (pdnmgWat == NULL) ||
		  (pmmdWat == NULL) ||
		  (pvnDescr1 == NULL) ||
		  (pvnDescr2 == NULL) ||
		  (pdnmmWat == NULL) ) { 
			  ErrPostEx(SEV_FATAL,0,10,"Out of Memory in SolvateIt");
			  return 0;}
	  pflvWat[0] = Twater[i][0];
	  pflvWat[1] = Twater[i][1];
	  pflvWat[2] = Twater[i][2];
	  pflvWat[3] = (float) 1.0; /* occupancy */
	  pflvWat[4] = (float) 99.99; /* maximum B factor */
	  iAtomCount++;
	  
	  /* test by writing out each keeper water to a PDB file as a HETATM line */
/*
fprintf(fp,"HETATM%5d  O   HOH 1%4d    %8.3f%8.3f%8.3f  1.00  0.00           O  \n",(int) iAtomCount, (int) iWaterCount, (float) pflvWat[0], 
		(float) pflvWat[1], (float) pflvWat[2]);
*/

	  paldWat->pflvData = pflvWat;
	  paldWat->iFloatNo = (Int1) 4; /* we add on max B factor to show error on these atoms */
      paldWat->pfbParent = (PFB) pmadWat;
	  paldWat->bMe = (Byte) AM_ALD;
	  paldWat->bWhat =(Byte) 0;
      paldWat->pvnalLink = pvnalWat; /* link to get model id */ 
      paldWat->iUniqueId = iAtomCount;  
      paldWat->cAltConf = ' '; 
    /*  paldWat->iCoordSet (Int1) 0 ;   index from ASN.1 model coord set  */
	  pvnalWat->choice = (Uint1) iModel;
	  pvnalWat->data.ptrvalue = (Pointer) paldWat;
	  pvnalWat->next = NULL;
	  pmadWat->pfbParent = (PFB) pmgdWat;
      pmadWat->bMe = (Byte) AM_MAD;
	  pmadWat->bWhat = (Byte) 0;
      pmadWat->pvnmaLink = pvnmaWat; /* element number is pvnmaLink->choice */
	  pmadWat->pcAName = prgWat->atoms->name;   /*   CharPtr pcAName points to PDB name in the Local Dictionary - not a StringSave alloced string */
      pmadWat->iIndex = prgWat->atoms->id; /* MMDB index  */
      pmadWat->pvnalLocate = pvnalWat; 
	  pvnmaWat->choice = 8; /* Atomic number of Oxygen */
	  pvnmaWat->data.ptrvalue = (Pointer) pmadWat;
	  pmgdWat->pfbParent = (PFB) pmmdWat;
      pmgdWat->bMe = (Byte) AM_MGD; 
      pmgdWat->bWhat = (Byte) DICT_LOCAL;  /* This refers to  LOCAL DICTIONARY FOR WATER */  
      pmgdWat->pdnmgLink = pdnmgWat;
      pmgdWat->pcGraphName = StringSave("HOH");  /* PDB 3-letter code */

   /*   pmgdWat->pcGraphNum  unnecessary - the PDB numbering string e.g . 38A */
      pmgdWat->iIDict  = iLocalDict;  /* index into local dictionary graph number used for added water  */
	  pvnStrings = NULL;
	  pvnStrings = prgWat->iupac_code;
      if (pvnStrings != NULL) 
         pmgdWat->pcIUPAC = StringSave((CharPtr) pvnStrings->data.ptrvalue); /* IUPAC name X = invalid graph code */  
      pmgdWat->pcPDBComment = NULL ; /* pointer into name of heterogen in dict */
	  pmgdWat->iAtomCount = (Int4) 1;
      pmgdWat->pvnmaAHead = pvnmaWat; /* the atom list (children) */
      pdnmgWat->choice = 1; /* residue id */
	  pdnmgWat->data.ptrvalue = (Pointer) pmgdWat;
	  pmmdWat->pfbParent = pmmdHere->pfbParent;
	  pmmdWat->bMe = (Byte) AM_MMD;
	  pmmdWat->bWhat = ((Byte) AM_SOL );
      pmmdWat->pdnmmLink = pdnmmWat;


      /* ValNodePtr pMolDescr; *//* the ASN.1 molecule descr */
	  /* this needs a new copy of the data structure */
	  /* in this case it is a list of ValNodes */

	  pvnDescr1->choice = BiomolDescr_name;
	  pvnDescr1->data.ptrvalue = (Pointer) StringSave("1");	           
	  pvnDescr2->choice = BiomolDescr_molecule_type;
	  pvnDescr2->data.intvalue = Molecule_type_solvent;
	  pvnDescr1->next = pvnDescr2;
	  pvnDescr2->next = NULL;
	  pmmdWat->pMolDescr = pvnDescr1; /* attach the molecule description */

/* Example of ASN.1 chemical graph entry for a specific water molecule */  
	  /****
{ id 108, descr { name "1", molecule-type solvent }, residue-sequence { { id 1, name " 215 ", residue-graph local 1 } } } },
     ***/

	  pmmdWat->pcMolName = StringSave(ParentMolName((PFB) pmadAtom ));
	  /* PDB solvent uses the chain name of the current atom being solvated */
      pmmdWat->iResCount = 1; /* number of graphs */
      pmmdWat->pdnmgHead = pdnmgWat;  /* the list of model graphs (children)  */
	  iGraphCount++;
	  pdnmmWat->choice = iGraphCount; /****/
	  pdnmmWat->data.ptrvalue = (Pointer) pmmdWat;
	  pdnmmEnd->next = pdnmmWat; /* finally added to the structure ! */ 
	  pdnmmWat->last = pdnmmEnd; 
	  pdnmmEnd = pdnmmWat; 
     /* then insert these in to  the BD tree */
	  /* so that new waters will avoid these placed ones */

 	  AddToBDTree((PFB)pmadWat,iModel,&(pwsHead->pbdTreeHead));
	  

	} /* if this water has coords to add */
} /* for loop over all waters in box */

/*printf("Added [%d] waters\n", (int) iWater); */
return iWater;
} /* SolvateIt */



void AddWaterLocalDictionary(PMSD pmsdThis) {

/*  add in the ASN.1 local dictionary  AND remove the vector of asn.1 oder information to clear for writers */    
/* caveat - difficult to use HOH in movie... unless all match to the same maximal chemical graph */
/*
/* the ASN.1 Residue-graph (text string) goes onto the end of the Biostruc-graph 
as "residue-graphs" */ 

AsnIoBSPtr    aibp = NULL;
BiostrucGraphPtr pbsgThis = NULL;
ResidueGraphPtr prgHere = NULL;
Char cDictLine[300];
ByteStorePtr pBSDict = NULL;
PMLD pmldTo = NULL;
PDNML pdnmlTo = NULL;


pBSDict = BSNew(1000);
if (pBSDict != NULL) {
	sprintf(cDictLine,"Residue-graph ::= { id 1, descr { name \"HOH\", other-comment \"TraDES COMPUTED WATERS CWVH 2012\" }, residue-type other, iupac-code { \"X\" }, atoms { { id 1, name \" O  \", iupac-code { \" O  \"}, element o } }, bonds { } } ");
	BSWrite(pBSDict, cDictLine, StringLen(cDictLine)+1 );
/* printf("%s\n",cDictLine);	*/
	aibp = AsnIoBSOpen ("r", pBSDict);
	if (aibp == NULL || aibp->aip == NULL) return;
	prgWat = ResidueGraphAsnRead(aibp->aip, NULL);
	prgHere = pmsdThis->pDictLocal;
	if (prgHere == NULL) {
		pmsdThis->pDictLocal = prgWat;
		iLocalDict = 1;
	}
	else
	{
		while (prgHere->next != NULL) {
			prgHere = prgHere->next;
		}
	/* add it to the end of the current local dictionary */
		prgHere->next = prgWat;
		prgWat->id = prgHere->id + 1;
		iLocalDict = prgWat->id;
	}
/* printf("Dictionary Index %d", (int) iLocalDict); */
/* RebuildChemGraph in mmdbapi4.c does this - just ensure all the right information is in the data structure */
	AsnIoBSClose (aibp);
	BSFree(pBSDict);
}

/* we also remove the vectors of ppAsnOrder to force mmdbapi to reconstruct these on writing out the file */
	pdnmlTo = pmsdThis->pdnmlModels;
		if (pdnmlTo != NULL) {
		      pmldTo = (PMLD)(pdnmlTo->data.ptrvalue);
			  if (pmldTo != NULL)
		         if(pmldTo->ppAsnOrder) {
			       PTRVectorFree(pmldTo->ppAsnOrder,0);
			       pmldTo->ppAsnOrder = NULL;
				 } 		
		}


}


Int4 ComputeSolvation (DValNodePtr pdnListPmmds, Int2 iModelNum, Boolean bSkipInternal, Int4 iLayers) {

ValNodePtr  pvnAtomList;
DValNodePtr pdnHere;
DValNodePtr pdnFirstSolvent = NULL;
DValNodePtr pdnLastSolvent = NULL;



PWS pwsHead = NULL;
PWS pwsThis = NULL;
vec middle;
Int2  numAtoms;
Int4 iTotalWater = 0;
int i;

PMMD pmmd1;
PMGD pmgd1;
PMAD pmad1;
PALD pald1;

PDNMG pdnmgHere;
PVNMA pvnmaHere;

    /* For each molecule in the list */

if (pdnListPmmds == NULL) return 0;
pdnHere = pdnListPmmds;
pdnmmEnd = NULL; 
iAtomCount = 0;
iGraphCount = 0;
iWaterCount = 0;/* zero out the global pointer and variables for any repeat calls */
/* printf("filling up world list in ComputeSolvation \n"); */

pmmd1 = (PMMD) pdnHere->data.ptrvalue;
if (IsSolvent((PFB) pmmd1)) return 0;  /* stop if first molecule is solvent  */

numAtoms = 0;
while (pdnHere != NULL) {
	pmmd1 = (PMMD)(pdnHere->data.ptrvalue);
/* We keep waters that are already present on database structures in initial BD-tree */ 
		pdnmgHere = pmmd1->pdnmgHead;
		/* For each residue in the molecule */
		while (pdnmgHere != NULL) {
	  	    pmgd1 = (PMGD)(pdnmgHere->data.ptrvalue);
		    pvnmaHere = pmgd1->pvnmaAHead;
		    /* For each atom */
	  		while (pvnmaHere != NULL) {
				pmad1 = (PMAD)(pvnmaHere->data.ptrvalue);
				if ((pald1 = GetAtomLocs (pmad1, iModelNum)) != NULL) {
				    /* If the atom is not a hydrogen */
					if ((((PVNMA)(pmad1->pvnmaLink))->choice != 1 /* ATOMIC_NO_H */)  &&
							(pald1->pflvData != NULL)) {
							/* ... Add it to the world of atoms */
							pwsThis = AddtoWorld (pwsHead,iModelNum,(PFB)(pmad1));
							numAtoms++;
							if (pwsHead == NULL) pwsHead = pwsThis;
					}
				}
				pvnmaHere = pvnmaHere->next;
			}
		pdnmgHere = pdnmgHere->next;
		}
	pdnHere = pdnHere->next;
}

InstantiateWorld (1, pwsHead);

/* printf("World Instantiated with %d atoms\n",numAtoms); */
AddWaterLocalDictionary(ToMSDParent((PFB) pmmd1));
pdnFirstSolvent = NULL;
pdnLastSolvent = NULL;

for (i = 1; i <= iLayers; i++) {

	if (i == 1) {  /* first case set solvation to be all non-water molecules in the molecule list, assuming all waters are at the end of the list */
		pdnFirstSolvent = pdnListPmmds;
		pdnHere = pdnListPmmds; /* walk the molecule list to find the first x-ray water */
		while (pdnHere != NULL) {
            pmmd1 = (PMMD) pdnHere->data.ptrvalue;
			if (!IsSolvent((PFB) pmmd1)){
				pdnLastSolvent = pdnHere;
			}
			pdnHere = pdnHere->next;
		}/* sets up pointers to the beginning and end of non-solvent molecules in x-ray structure layer just added */
	}

	if (i > 1) { 
		pdnFirstSolvent = pdnLastSolvent->next; /* everything here to the end is newly added solvent from the previous layer */
		pdnHere = pdnFirstSolvent;
        while (pdnHere != NULL) {
			pdnLastSolvent = pdnHere;
			pdnHere = pdnHere->next;
		}/* sets up pointers to the beginning and end of previous solvent layer - for the next layer of solvation */
	}




/* Do solvate for each molecule in the list of molecules to solvate (pdnFirstSolvent to pdnLastSolvent) */
pdnHere = pdnFirstSolvent;
if (pdnHere != NULL)
do {
	pmmd1 = (PMMD)(pdnHere->data.ptrvalue);
		pdnmgHere = pmmd1->pdnmgHead;
		/* For each residue in the molecule */
			while (pdnmgHere != NULL) {
	  			pmgd1 = (PMGD)(pdnmgHere->data.ptrvalue);
	    		pvnmaHere = pmgd1->pvnmaAHead;
				/* For each atom */
				while (pvnmaHere != NULL) {
					pmad1 = (PMAD)(pvnmaHere->data.ptrvalue);
					if ((pald1 = GetAtomLocs (pmad1, iModelNum)) != NULL) {
		  		  /* If the atom is not a hydrogen and it does have coordinates */
		 				if ((((PVNMA)(pmad1->pvnmaLink))->choice != 1 /* ATOMIC_NO_H*/) &&
							(pald1->pflvData != NULL)) {
								middle[0] = pald1->pflvData[0];
								middle[1] = pald1->pflvData[1];
								middle[2] = pald1->pflvData[2];
								/* get the 6 Anstrom Neighbor list */
								pvnAtomList = FindAtomsIn (pwsHead, middle, 6.0);
							    iTotalWater += SolvateIt(pmad1, pvnAtomList, pdnListPmmds, iModelNum, pwsHead, bSkipInternal );
								if (pvnAtomList != NULL) ValNodeFree(pvnAtomList);	
						} /* skip not suitable  H or no atom location */				
					} /* if atom 1 has a location */
					pvnmaHere = pvnmaHere->next;
	    		} /* while pvnmaHere */
				pdnmgHere = pdnmgHere->next;
			} /* while pdnmgHere */
	pdnHere = pdnHere->next;
} while (pdnHere->last != pdnLastSolvent);  
/* test condition is - did we just complete the last molecule to be solvated and are one molecule into the newly added layer ? */
 

} /* for each solvation layer to be added */



FreeAllWorlds ();



return iTotalWater;
}






void DoWipeSolvent(PMSD pmsdThis)
{
	
    PDNMM pdnmmThis = NULL;
	PDNMM pdnmmNext = NULL;
    PMMD pmmdThis = NULL;

    pdnmmThis = pmsdThis->pdnmmHead;
	while(pdnmmThis != NULL)
	    {   /*molecule*/
			pdnmmNext = pdnmmThis->next; /* save the next pointer */
/*			printf("In WipeSolvent\n"); */
			pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
			if (IsSolvent((PFB)pmmdThis)) {
		  	  if (pdnmmThis->last != NULL) /* first on the chain */
			    pdnmmThis->last->next = pdnmmThis->next;
			  else {  /* fix the PMSD link */
			   pmsdThis->pdnmmHead = pdnmmThis->next;
			  }
			  if (pdnmmThis->next != NULL) /* last on the chain */
			    pdnmmThis->next->last = pdnmmThis->last;
			  pdnmmThis->next = NULL;
			  pdnmmThis->last = NULL;
			  DValNodeFreeData(pdnmmThis, (pFreeFunc)FreeMMD);
/* 			  printf("Data Freed\n"); */ /* only called once! */  
		   }
			pdnmmThis = pdnmmNext; /* next molecule */
	    }   /*molecule*/
     
}


void LIBCALLBACK DoPCACenter(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr) {

FloatLoPtr pfmean = NULL;
PALD paldLocations;
int j;

    if (pfbThis == NULL) return;  
    if (pfbThis->bMe != (Byte) AM_MAD) return;
 	paldLocations = GetAtomLocs((PMAD) pfbThis, (Int2) iModel);
	if (paldLocations == NULL) return;
	if (paldLocations->iFloatNo == 0) return;
	if (paldLocations->pflvData == NULL) return;
    pfmean = (FloatLoPtr) ptr;
	if (pfmean == NULL) return;
	for (j = 0; j < 3; j++) {
		paldLocations->pflvData[j] -= pfmean[j];
	}
}



void LIBCALLBACK DoPCARotate(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr) {

FloatLoPtr PNTR ppfSTM = NULL;
PALD paldLocations;
vec pfInter = {0.0, 0.0, 0.0};
int j, k,l;

    if (pfbThis == NULL) return;  
    if (pfbThis->bMe != (Byte) AM_MAD) return;
 	paldLocations = GetAtomLocs((PMAD) pfbThis, (Int2) iModel);
	if (paldLocations == NULL) return;
	if (paldLocations->iFloatNo == 0) return;
	if (paldLocations->pflvData == NULL) return;

    ppfSTM = (FloatLoPtrPtr) ptr;
	if (ppfSTM == NULL) return;
	for (j = 0; j < 3; j++) {
			pfInter[j] = paldLocations->pflvData[j]; 
		}   
	for (k = 0; k < 3; k++) {
		paldLocations->pflvData[k] = 0.0;
	    for (l = 0; l < 3; l++) {
				paldLocations->pflvData[k] += pfInter[l] * ppfSTM[l][3-k-1]; 
			}
		}
}





void LIBCALLBACK DoPCACenterSolids(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr) {

FloatLoPtr pfmean = NULL;
int i, j;
PMOD pmodThis;

	if (pfbThis == NULL) return;  
    if (pfbThis->bMe != (Byte) AM_MOD) return;
 	pmodThis = (PMOD) pfbThis;
    if (pmodThis->ppflObject == NULL) return;
	pfmean = (FloatLoPtr) ptr;
	if (pfmean == NULL) return;
	for (i = 0; i < pmodThis->iCoordNo; i++) {	
      for (j = 0; j < 3; j++) {
		pmodThis->ppflObject[i][j] -= pfmean[j];
		}
	}
}



void LIBCALLBACK DoPCARotateSolids(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr) {

FloatLoPtr PNTR ppfSTM = NULL;
vec pfInter = {0.0, 0.0, 0.0};
int i, j, k,l;
PMOD pmodThis;

	if (pfbThis == NULL) return;  
    if (pfbThis->bMe != (Byte) AM_MOD) return;
 	pmodThis = (PMOD) pfbThis;
    if (!pmodThis->ppflObject) return;
	ppfSTM = (FloatLoPtrPtr) ptr;
	if (ppfSTM == NULL) return;

	for (i = 0; i < pmodThis->iCoordNo; i++) {	
		for (j = 0; j < 3; j++) {
			pfInter[j] = pmodThis->ppflObject[i][j]; 
		}   
		for (k = 0; k < 3; k++) {
			pmodThis->ppflObject[i][k] = 0.0;
		    for (l = 0; l < 3; l++) {
				pmodThis->ppflObject[i][k] += pfInter[l] * ppfSTM[l][3-k-1]; 
			}
		}
	}
}





void LIBCALLBACK DoAddSolvent(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr) {

	ValNodePtr pvnmaAtom;
	PMAD pmadThis;
    PMGD pmgdThis;
	PALD paldCoords;
  
   /* callback function that adds solvent molecules to a global linked list of PMAD nodes */
    if (pfbThis == NULL) return;  
 /*   printf("in DoAddsolvent\n"); */
	if (IsAtomNode(pfbThis)) {
/*	printf("Found an atom node/n"); */
    pmadThis = (PMAD) pfbThis;
    paldCoords=GetAtomLocs(pmadThis,(Int2) iModel/* model num*/);
	paldCoords = GetLocation(paldCoords, 0, ' ');
	if (paldCoords == NULL) return; /* no location */
    if (paldCoords->pflvData == NULL) return; /* no coordinate data */
	  
 /*   if (pmadThis->pcAName != NULL) return;  these are the ones we want! No dictionary name - this is fixed  */ 
    pmgdThis = GetParentGraph(pfbThis);
/*	printf("%s\n",pmgdThis->pcGraphName); */
/*
	if (GetParentMol(pfbThis) == NULL) printf("Cannot get parent molecule \n"); */

	if (pmgdThis->pcGraphName[0] != 'H') return;
	if (pmgdThis->pcGraphName[1] != 'O') return;
	if (pmgdThis->pcGraphName[2] != 'H') return;
/*	printf("Found one...that i added... \n"); */
	pvnmaAtom = ValNodeNew(NULL);
	if (pvnmaAtom == NULL) {
      	ErrPostEx(SEV_ERROR,0,10,"Memory Allocation Error on making linked list of solvent Atoms for PCA analysis.");
				 return;
    }
	pvnmaAtom->data.ptrvalue = (Pointer) pmadThis;
	if (pvnSolventList == NULL) {
         pvnSolventList = pvnmaAtom;
	     pvnCur = pvnmaAtom;
     }	
	else {
		pvnCur->next = pvnmaAtom;
   	    pvnCur = pvnmaAtom;
	 }

  }
}






void SolvateAndPCA(PMSD pmsdHead, Int2 ModelNum, Int4 structureNo, Int4 SkipPCA , Int4 iLayers) {

  ValNodePtr pvnSolvatedAtoms = NULL;
  ValNodePtr pvnHere = NULL;
  Int4 iSolvAtoms = 0, iM = 0;
  /* pointer to Solvated atom vector */
  FloatLoPtr *ppfVUA = NULL;
  FloatLoPtr pfVUA = NULL;
  Int4 iRows = 0;
  Int4 iWater = 0;
  int iBig, iMed, iSmall, i, j; 
/*   int i; */

  PMAD pmadAtom = NULL;
  PALD paldCoords = NULL;
  FloatLoPtr *ppfCov = NULL;  /* covariance matrix */
  FloatLoPtr pfEigen = NULL, pfInter = NULL;
  FloatLoPtr  *ppfSymm = NULL, *ppfSTM = NULL, *ppfSwap = NULL;
/*  FloatLoPtr *ppfRowProj = NULL, *ppfColProj = NULL;  */
  FloatLoPtr pfmean = NULL, pfstddev = NULL;
  Int4 linelen = 0;
  Char line[2000];
  
  /* variables for removal of ppAsnOrder */
  PMLD pmldTo = NULL;
  PDNML pdnmlTo = NULL;

  if (traj_quiet==VERBOSITY_VERBOSE)
	printf("Adding %d layers of water molecules..\n",iLayers); 
  iWater = ComputeSolvation(pmsdHead->pdnmmHead, ModelNum,TRUE,iLayers); /* TRUE is to skip internal solvation */
  if (traj_quiet==VERBOSITY_VERBOSE)
	printf("Solvated structure with %d water molecules.\n",iWater);
  /* traverse the structure and pull out all the AM_SOL PMADs into a linked list with PMADs replacing the Voronoi list with Solvated water list */
  /* this adds to the global pointer pvnSolventList */
  /* compute PCA on Solvated Atom List */
  if (traj_quiet==VERBOSITY_VERBOSE)
	printf("Computing PCA Axis Transformation with solvent layer..\n");
  /* printf("Gathering Solvent for PCA\n"); */
  pvnSolventList = NULL;
  TraverseAtoms(pmsdHead->pdnmmHead, ModelNum,0, NULL, (pNodeFunc) (DoAddSolvent));
  pvnSolvatedAtoms = pvnSolventList;  
  if (pvnSolventList != NULL) { 
  iSolvAtoms = 0;
  pvnHere = pvnSolvatedAtoms;
  while(pvnHere != NULL)
  {
	  iSolvAtoms++;
	  pvnHere=pvnHere->next;
  }
  /* printf("Counted %d Solvent Atoms \n",(int) iSolvAtoms);  this includes crystallographic waters */
  /* fill matrix rows */
  iM = 0;
  pvnHere = pvnSolvatedAtoms;
  while(pvnHere != NULL)
  {
	  pmadAtom = (PMAD) pvnHere->data.ptrvalue;
	  paldCoords=GetAtomLocs(pmadAtom,(Int2) 1/* model num*/);
	  /* Now extract the actual coordinates */
	  paldCoords = GetLocation(paldCoords, 0, ' ');
	  pfVUA = (FloatLoPtr) MemNew((size_t)sizeof(FloatLo) * 3);
	  if (pfVUA == NULL) {
			ErrPostEx(SEV_ERROR,0,10,"Memory Allocation Error on array to store unsolvated Atoms for PCA analysis.");
			return;
		}
	  pfVUA[0] = (FloatLo) paldCoords->pflvData[0]; /* x coordinate */
	  pfVUA[1] = (FloatLo) paldCoords->pflvData[1];   /* y coordinate */
      pfVUA[2] = (FloatLo) paldCoords->pflvData[2];   /* z coordinate */	 
	  iM++;
/* printf("Row %d packed: %f %f %f\n", (int) iM , (float) pfVUA[0], (float) pfVUA[1], (float) pfVUA[2]);  */
	  pvnHere->data.ptrvalue = (Pointer) pfVUA; /* hang the new vector here */	   
	  pvnHere=pvnHere->next;
  }
/* printf("Done packed: %f %f %f\n",(float) pfVUA[0], (float) pfVUA[1], (float) pfVUA[2]); */ 
  /* allocate matrix */
  ppfVUA = (FloatLoPtr *) MemNew((size_t) (sizeof(FloatLoPtr) * (iSolvAtoms +1)) );
  if (ppfVUA == NULL) {
	   ErrPostEx(SEV_ERROR,0,10,"Memory Allocation Error on array to store unsolvated Atoms for PCA analysis.");
	   return;
   }
   iM = 0;
/*   printf("ppfVUA allocated\n"); */
   pvnHere = pvnSolvatedAtoms;
   while(pvnHere != NULL) {
	   if (pvnHere->data.ptrvalue != NULL) {
		ppfVUA[iM] = (FloatLoPtr) pvnHere->data.ptrvalue;
		iM++;
	   }
		pvnHere = pvnHere->next;
	}
 /* printf("ppfVUA filled\n"); */
   iRows = iM;
  /* Call components of slri_stats.c */
	STAT_ColumnMean(ppfVUA,iRows,3,&pfmean);
	STAT_ColumnStdDev(ppfVUA,iRows,3,pfmean,&pfstddev);
	STAT_CenterAndReduce(ppfVUA,iRows,3,pfmean,pfstddev,FALSE);
	STAT_CovariationMatrix(ppfVUA,iRows,3,&ppfCov);	

/*	STAT_DumpMeanAndStdDev(pfmean,pfstddev,3,CovarianceFile); */
   /* 3 projections, iVecNum/Size/Col is 3, iVecLength/Row is iSolvAtoms */  
	PCA_TDR (ppfCov, 3, &pfEigen, &pfInter, &ppfSTM);  /* Triangular decomposition */
	PCA_TDQL(ppfSTM, 3, pfEigen,  pfInter);             /* Reduction of sym. trid. matrix */
    ValNodeFreeData(pvnSolventList); /* done with this list */
    MemFree(ppfVUA); 
 

/* Ellipsoid axis is scaled so that a = sqrt(3 * PCA) in paper Yaroslav E. Ryabov,† Charles Geraghty,† Amitabh Varshney,‡ and
David Fushman JACS 2006 ; PCA is  ( Eigen score / N ) for an entire sample .  Endpoints of a, b, c are now in pfEigen - checked with RasMol distances */
	pfEigen[0] = sqrt(3* (pfEigen[0] / iRows));  
	pfEigen[1] = sqrt(3* (pfEigen[1] / iRows));
	pfEigen[2] = sqrt(3* (pfEigen[2] / iRows));

	/* stupid sort! */
	iBig = 0;
	iMed = 0;
	iSmall = 0;
	for (i = 1; i <3; i++) {
		if (pfEigen[i] >= pfEigen[iBig]) 
		iBig = i;
	}
	for (i = 1; i <3; i++) {
		if (pfEigen[i] <= pfEigen[iSmall]) 
		iSmall = i;
	}
	for (i = 0; i < 3; i++) {
		if ((i != iBig) && (i != iSmall))
			iMed = i;  /* this handles the case of a tie */
	}
	/* copy and swap order of eigenvalues */
    pfInter[0] = pfEigen[0];
	pfInter[1] = pfEigen[1];
	pfInter[2] = pfEigen[2];
	pfEigen[0] = pfInter[iSmall];
	pfEigen[1] = pfInter[iMed];
	pfEigen[2] = pfInter[iBig];
	/* now reorg the eigenvectors in same order */
	ppfSwap = MATH_MatrixNew(3, 3);
	for (i = 0; i<3; i++) {
		for (j= 0; j<3; j++) { /* copy to holding matrix */
			ppfSwap[i][j] = ppfSTM[i][j];
		}
	}
	for (i = 0; i<3; i++) { /* transpose the eigenvector values by column */
		ppfSTM[i][0] = ppfSwap[i][iSmall];
		ppfSTM[i][1] = ppfSwap[i][iMed];
		ppfSTM[i][2] = ppfSwap[i][iBig];
	}

    MATH_VectorFree(&pfInter);
	MATH_MatrixFree(&ppfSwap,3);
	if (SkipPCA != TRUE) { 
		if (traj_quiet==VERBOSITY_VERBOSE)
			printf("Computing Coordinate Transform to Principal Axes.\n",iWater);

		/* Now that we have the mean and the rotation matrix, transform the entire coordinate set */
		TraverseAtoms(pmsdHead->pdnmmHead, ModelNum, 0, (Pointer) pfmean, (pNodeFunc) (DoPCACenter));
		TraverseAtoms(pmsdHead->pdnmmHead, ModelNum, 0, (Pointer) ppfSTM, (pNodeFunc) (DoPCARotate));
		/*	printf("Calling TraverseSolids\n"); */
		TraverseSolids(pmsdHead->pdnmsLink, ModelNum, 0, (Pointer) pfmean, (pNodeFunc) (DoPCACenterSolids));
		TraverseSolids(pmsdHead->pdnmsLink, ModelNum, 0, (Pointer) ppfSTM, (pNodeFunc) (DoPCARotateSolids));
	}
	pvnHere = NULL;
	sprintf(line,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n",  (int) iSolvAtoms,
	 (float) pfstddev[0],(float) pfstddev[1], (float) pfstddev[2],
     (float) pfEigen[2], (float) pfEigen[1], (float) pfEigen[0], LOG_TERM_CHAR); 

	linelen = StringLen(line);
	BSWrite(bspTempLog, line, linelen);

	MATH_VectorFree(&pfmean);
	MATH_VectorFree(&pfstddev);
	MATH_MatrixFree(&ppfCov,3);
	MATH_VectorFree(&pfEigen);
	MATH_MatrixFree(&ppfSymm,3);
	MATH_MatrixFree(&ppfSTM,3);

 } /* if pvnSolventList */
}


void RemoveSolvent (PMSD pmsdHead) {

  PMLD pmldTo = NULL;
  PDNML pdnmlTo = NULL;

	pdnmlTo = pmsdHead->pdnmlModels;
		if (pdnmlTo != NULL) {
		      pmldTo = (PMLD)(pdnmlTo->data.ptrvalue);
			  if (pmldTo != NULL)
		         if(pmldTo->ppAsnOrder) {
			       PTRVectorFree(pmldTo->ppAsnOrder,0);
			       pmldTo->ppAsnOrder = NULL;
				 } 		
		}
		DoWipeSolvent(pmsdHead);  
		RemoveWaterDictionary(pmsdHead); 
}


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
	sprintf(line,"%d, %d, %f, %f\n",(int)structureNo, (int) num, phi,psi);
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
Char ftmp2[PATH_MAX];

 /*  write arbitrary bytestore to a file in one chunk of memory, if you got it*/ 
    numbytes=BSLen(bspHere);
	if (numbytes==0)
	    return;
	pcTempBuf=(CharPtr)MemNew((size_t)numbytes);
    BSSeek(bspHere,0L,SEEK_SET);
    BSRead(bspHere,pcTempBuf,numbytes);
    StringCpy(ftmp2,filename);
  /*  StringCat(ftmp2,LOG_EXT); */
	if ((fp=FileOpen(ftmp2,"a"))==NULL) {
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






Int2 Main ()
{
	FILE *fp;
	Char sequence[MAXSIZE];
 	Char fnam[PATH_MAX];
 	Char fnamtrj[PATH_MAX];
 	Char fnamtrjnopath[PATH_MAX];
 	Char fnam2[PATH_MAX];
 	Char ftmp[PATH_MAX];
	Char writepath[PATH_MAX];
	Char pcUname[PATH_MAX];
	Char fileout[PATH_MAX];
	Int2 numAA;
	Int4 numtofold;
	TrajErr err;
	PMSD pmsdHead;
	BiostrucPtr bspBiostruc;
        PMLD pmldThis;
        PDNML pdnmlModel;
	PDNMS pdnmsModelstruc=NULL;
	PMMD pmmdNative=NULL;
	PMSD pmsd1;
	Int4 countstruc=0;
	FloatHi rmsdval=0.0, potzh=0.0, potbryant=0.0, potVoronoi = 0.0, potcrease=0.0;
	Int4 vPDBatomN = 0;
	Int4 processid,randseed;
	DValNodePtr pdnZhangAdjList = NULL;
	DValNodePtr pdnBryantAdjList = NULL;
	DValNodePtr pdnVoronoiAdjList = NULL;
	DValNodePtr pdnZhangAtmList = NULL;
	DValNodePtr pdnListPmmds;
	Int4Ptr piBryantPotential = NULL;
	Int4Ptr piZhangPtnl = NULL;
	Char chainid='-',chainhere='-';
        ValNodePtr vnpMolList,vnpHere;
	FILE *fpzh, *fpbl, *fpAtmInfo;
    Char timedate[25];			
	CharPtr decsequence;
	int numextres;
	
	CharPtr pc,pc2;
  	Int2 CreaseIncl=0,CreaseExcl=0,ModelNum,NativeModelNum=1;
	Boolean		      GlobalNonHtmlOutput=FALSE;
	BiostrucPtr		pbs;
	BiostrucIdPtr		pbi;
	BiostrucDescrPtr	pbd;
	ValNodePtr	vnpBioseq=NULL;
	BiostrucPtr bspTemp;
	AsnIoPtr aip;
	NcbiMimeAsn1Ptr nmap;
	Char progname[PATH_MAX];
	pFoldTrajParamBlock foldtrajprm;
	ValNodePtr vnpHead,vnpTemp;
	BioseqSetPtr bssp;
	BioseqPtr bsp;
	Int4 dbsize;
	Int4 writecount = 1;
	Int4 trystructure = 0;
	Int4 ramalimit = 0;
   	Int2    iDotLen = 0;

#ifdef OS_UNIX
	struct utsname utsbuf;
#endif
        
        ErrSetLogfile("error_trades.log", ELOG_APPEND|ELOG_BANNER);
        ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
        ErrSetMessageLevel(SEV_FATAL);
        ErrSetLogLevel(SEV_WARNING);

	/* we always build Model #1 */
	ModelNum=1;
	/* Method is set here */
	Method = 0; /* Dummy out Mobi Method - command line only */
	sprintf(progname,"TraDES v %s",FOLDTRAJ_VERSION);
	if (!GetArgs("TraDES - Trajectory Directed Ensemble Sampling Engine\nUsed to generate protein 3D structures by conformational space sampling.\n\n",20,StrucTrjArgs)) {
		printf("TO START: You need a *.trj file first (use seq2trj or str2trj). \n     For usage try: seq2trj - OR str2trj - \n\n");
		return 1;	/* insufficient arguments entered */
	}

	
/* be forgiving about adding extensions 
  CWVH  June 2012
    No-extension is user frustration! offer forgiveness for typing the extention , punch a null overtop of the dot  
   Int2    iDotLen = 0;
*/

	iDotLen = 0;
	if (StrucTrjArgs[0].strvalue != NULL) {
           iDotLen = StringLen(StrucTrjArgs[0].strvalue);
           if (iDotLen > 5) {
		if (StrucTrjArgs[0].strvalue[iDotLen - 4] == '.') {
		 ErrPostEx(SEV_WARNING,2,1,"Ignoring argument -f 3 letter extension after . on %s",StrucTrjArgs[0].strvalue);
                StrucTrjArgs[0].strvalue[iDotLen - 4] = '\0';      
		}
           }
	}
	iDotLen = 0;
	if (StrucTrjArgs[1].strvalue != NULL) {
           iDotLen = StringLen(StrucTrjArgs[1].strvalue);
           if (iDotLen > 5) {
		if (StrucTrjArgs[1].strvalue[iDotLen - 4] == '.') {
		 ErrPostEx(SEV_WARNING,2,1,"Ignoring argument -o 3 letter extension after . on %s",StrucTrjArgs[1].strvalue);
                StrucTrjArgs[1].strvalue[iDotLen - 4] = '\0';      
		}
           }
	}
	iDotLen = 0;
	if (StrucTrjArgs[8].strvalue != NULL) {
           iDotLen = StringLen(StrucTrjArgs[8].strvalue);
           if (iDotLen > 5) {
		if (StrucTrjArgs[8].strvalue[iDotLen - 4] == '.') {
		 ErrPostEx(SEV_WARNING,2,1,"Ignoring argument -c 3 letter extension after . on %s",StrucTrjArgs[8].strvalue);
                StrucTrjArgs[8].strvalue[iDotLen - 4] = '\0';      
		}
           }
	}


	StringCpy(fnamtrj,StrucTrjArgs[0].strvalue);
	USE_LOTS_RAM=1;
	BUILD_FRAGS_ONLY=StrucTrjArgs[12].intvalue;
	numtofold=StrucTrjArgs[6].intvalue; 

	if (StrucTrjArgs[7].intvalue==TRUE)
			traj_quiet=(Byte)VERBOSITY_QUIET;
		else
			traj_quiet=(Byte)VERBOSITY_VERBOSE;

	if (StrucTrjArgs[17].intvalue == TRUE)
		    traj_quiet=(Byte)VERBOSITY_STREAM;


	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Generating random seed..\n");

	processid = GetAppProcessID();
	
	if (StrucTrjArgs[2].intvalue)
		randseed=StrucTrjArgs[2].intvalue;	/* deterministic program */
	else
		randseed=GetSecs()*processid; 	   /* non-deterministic program if seed=0 */
	RandomSeed(randseed);
	StringCpy(ftmp,fnamtrj);
	pc=ftmp;
	do {
		pc2=StringStr(pc,DIRDELIMSTR);
		if (pc2)
			pc=pc2+1;
	} while (pc2);
	StringCpy(fnamtrjnopath,pc);

	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Opening MMDB API..\n");
    if (!OPENMMDBAPI(0,NULL)) {
	        ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI");
	        return 2;
        }


	if (StrucTrjArgs[13].intvalue == TRUE) {
		/* Initialize Ramachandran Output Bytestores */
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


        /* load native or comparison ASN.1 Structure */
 	if (StrucTrjArgs[8].strvalue!=NULL) {
		if (traj_quiet==VERBOSITY_VERBOSE)
			printf("Loading native structure..\n");
		StringCpy(ftmp,CFG_local_datafilepath);
		StringCat(ftmp,StrucTrjArgs[8].strvalue);
		StringCat(ftmp,MMDB_EXT);
		NativeModelNum=1;
		pmsd1=LoadABiostruc(ftmp,FALSE,ALLMDL,&NativeModelNum);
		if (pmsd1==NULL) {
			/* error occurred */
                        StringCpy(ftmp,CFG_local_datafilepath);
			StringCat(ftmp,StrucTrjArgs[8].strvalue);
			StringCat(ftmp,".cn3");
			pmsd1=LoadABiostruc(ftmp,FALSE,ALLMDL,&NativeModelNum);
			if (pmsd1 == NULL) {
	                        StringCpy(ftmp,CFG_local_datafilepath);
				StringCat(ftmp,StrucTrjArgs[8].strvalue);
				StringCat(ftmp,".prt");
				pmsd1=LoadABiostruc(ftmp,FALSE,ALLMDL,&NativeModelNum);
				if (pmsd1 == NULL) {
					ErrPostEx(SEV_FATAL,1,9,"Unable to load native structure file as %s.val or %s.cn3 or %s.prt ",StrucTrjArgs[8].strvalue,StrucTrjArgs[8].strvalue,StrucTrjArgs[8].strvalue);
				return 3;
				}
			}
		}
	        vnpMolList=ConvertNode((PFB)pmsd1,AM_MMD);
                vnpHere=vnpMolList;
                while(vnpHere) {
                        pmmdNative=(PMMD)(vnpHere->data.ptrvalue);
                        chainhere=(pmmdNative->pcMolName)[0];
                        /* ignore heterogens */
                        if (((pmmdNative->bWhat) & AM_PROT) && ((chainid=='-') || (chainhere==chainid)))
				break;		
                        vnpHere=vnpHere->next;
                }
		if (!pmmdNative || vnpHere==NULL) {
                  ValNodeFree(vnpMolList);
		  ErrPostEx(SEV_FATAL,0,0,"Unable to load Native molecular model data (correct chain given?)");
		  return 6;
		}
        ValNodeFree(vnpMolList);
	}
	else
		pmmdNative=NULL;

	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Reading Zhang-DeLisi Potential..\n");
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,ZHANG_POTENTIAL);
	if ((fpzh=FileOpen(ftmp,"r"))==NULL) {
	        	ErrPostEx(SEV_FATAL,1,7,"File not found %s",ftmp);
            	return 7;
        }
    
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,ZHANG_ATOMS);
    if ((fpAtmInfo=FileOpen(ftmp,"r"))==NULL) {
		    	ErrPostEx(SEV_FATAL,1,7,"File not found %s",ftmp);
            	return 7;
	}
	if (LoadZhangPotential (&piZhangPtnl, fpzh)!=ERR_SUCCESS) {
	        	ErrPostEx(SEV_FATAL,1,9,"Error in LoadZhangPotential");
            	return 9;
	}
	LoadZhangAtmList (&pdnZhangAtmList, fpAtmInfo);
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Reading Bryant-Lawrence Potential..\n");
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,BL_POTENTIAL);
	if ((fpbl=FileOpen(ftmp,"r"))==NULL) {
	        	ErrPostEx(SEV_FATAL,1,8,"File not found %s",ftmp);
            	return 8;
	}
	if (LoadBryantPotential (&piBryantPotential, fpbl, TRUE /*bUsingPep*/, FALSE /*bNegateInput*/)!=ERR_SUCCESS) {
	        	ErrPostEx(SEV_FATAL,1,10,"Error in LoadBryantPotential");
            	return 10;
	}
	FileClose (fpzh);
	FileClose (fpbl);
	FileClose (fpAtmInfo);
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Unpacking ASN.1 file..\n");

	if (UnPackAsnTrajGraph(fnamtrj,&numAA,sequence,NULL,&pbi,&pbd,&vnpBioseq)==NULL) {
			ErrPostEx(SEV_FATAL,2,3,"Unable to read trajectory distribution or file not found: %s.trj",fnamtrj);
	}
	if (StrucTrjArgs[3].intvalue==FALSE) {
		if (traj_quiet==VERBOSITY_VERBOSE)
			printf("Removing ASN.1 file..\n");
		StringCpy(fnam,fnamtrj);
		StringCat(fnam,ASN_EXT);
		FileRemove(fnam);
	}
	/* set Crease Energy window */
	CreaseExcl=3;
	CreaseIncl=numAA;
	
	/* Build the peptide with sequence as a Biostruc with no co-ords
	   yet or with extended strand co-ordinates */
    if (traj_quiet==VERBOSITY_VERBOSE)
        printf("Opening rotamer library..\n");
    if (LoadRotLib()!=ERR_SUCCESS) {
	 	ErrPostEx(SEV_FATAL,1,8,"Cannot open rotamer library");
       	return 8;
	}
	TmpNam(tmpskelfname);
    if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Creating chemical graph...\n");
	/* make the ASN.1 file for this protein's chemical graph */
	/* extended amino acids are in encoded format here EXTAA_ENCODE */
	BuildSkelASN(sequence,tmpskelfname);

/* echo all this information - this is all logged. */ 
    if (traj_quiet==VERBOSITY_VERBOSE) {
		decsequence=DecodeSequence(sequence,EXTAA_PARENT);
		printf("Preparing to fold sequence %s\n",decsequence);
		decsequence=MemFree(decsequence);
		printf("%ldx%ld Trajectory Distributions\n",(long int)TRAJDIV,(long int)TRAJDIV);
		if (StringCmp(CONSTRAINT_FILE,""))
			printf("Constraint file: %s\n",CONSTRAINT_FILE);
		if (TRAJTYPE!=TRAJ_NA) {
			printf("Trajectory Distribution: ");
			switch (TRAJTYPE) {
				case TRAJ_UNIFORM:
					printf("Uniform\n");
					break;
				case TRAJ_STANDARD:
					printf("Amino-Acid Based\n");
					break;
				case TRAJ_SSTRU:
					printf("1-State Secondary Structure Prediction\n");
					break;
				case TRAJ_GOR:
					printf("3-State Secondary Structure Prediction\n");
					break;
				default:
					printf("Unknown\n");
			}
		}
		printf("Average timeout: %ld\n",(long int)TIMEOUT);
		printf("Backbone error tolerance: %6.2f degrees squared\n",BACKBONE_ERROR_TOLERANCE);
		printf("Backbone precision: %6.4f degrees\n",BACKBONE_PRECISION);
		printf("Number of Rotamer Tries: %d/chi angle\n",NUM_ROT_TRIES);
		printf("Atom bounciness (Backbone): %6.2f%%\n",ATOM_BOUNCINESS_BB*100.0);
		printf("Atom bounciness (Sidechains): %6.2f%%\n",ATOM_BOUNCINESS_SC*100.0);
                printf("Bumpcheck Hydrogens: ");
        if (BUMPCHECK_HYDROGEN)
                printf("TRUE\n");
        else
                printf("FALSE\n");
        printf("Increment size: %6.2f degrees\n",INCSIZE);  /* not logged */
        printf("Fragment Tunnelling Probability: %1.3f\n",TUNNEL_PROB);  /* not logged */
        printf("Initial backbone test angle: %6.2f degrees\n",START_BACKBONE); /* not logged */
        if (WALKTYPE==WALK_CA)
             printf("Beginning CA random walk...\n");
        else if (WALKTYPE==WALK_PHIPSI)
             printf("Beginning Phi-Psi random walk...\n");
        if (TGUNITS==UNITS_ARBITRARY)
             printf("(Trajectory distributions created with arbitrary freqeuncy units)\n");
		printf("Note: these settings are stored in the .trj trajectory file.\nChanging their values in the configuration file (.structrjrc or structrj.ini)\nwill only have an effect on newly created trajectory ASN.1 files\n(created with seq2trj or str2trj); see .structrjrc/structrj.ini for details on parameters.\n");
	}


	bspTempLog=BSNew(BS_STARTSIZE);
/* compute log file name */
    sprintf(LogOutName,"%s_%s.log",StrucTrjArgs[10].strvalue ? StrucTrjArgs[10].strvalue : "",StrucTrjArgs[1].strvalue ? StrucTrjArgs[1].strvalue:fnamtrjnopath);
	sprintf(ftmp,"%s_%s.csv",StrucTrjArgs[10].strvalue ? StrucTrjArgs[10].strvalue : "" ,StrucTrjArgs[1].strvalue ? StrucTrjArgs[1].strvalue:fnamtrjnopath);
	StringCpy(fileout,"__");
	StringCat(fileout,ftmp);
				
/*	sprintf(LogOutName,"%s%s_%s",
			StrucTrjArgs[10].strvalue ? "_" : "",
			StrucTrjArgs[10].strvalue ? StrucTrjArgs[10].strvalue : "" ,
			fnamtrjnopath); */
	BSprintf(bspTempLog,"\nTraDES v %s log report",FOLDTRAJ_VERSION);
	BSprintf(bspTempLog,"\nTrajectory File: %s%s",CFG_local_datafilepath,fnamtrj);
	if (StrucTrjArgs[4].intvalue==TRUE || StrucTrjArgs[5].intvalue==TRUE)
		BSprintf(bspTempLog,"\tStructure File Base Name: %s%s%s",StrucTrjArgs[10].strvalue ? StrucTrjArgs[10].strvalue : "",StrucTrjArgs[10].strvalue ? "_" : "",StrucTrjArgs[1].strvalue ? StrucTrjArgs[1].strvalue:fnamtrj);
	BSprintf(bspTempLog,"\tRandom Seed: %ld\t# Generated: %ld\tStart Numbering at: %ld",randseed,StrucTrjArgs[6].intvalue,StrucTrjArgs[11].intvalue);
	if (StrucTrjArgs[8].strvalue!=NULL) {
		BSprintf(bspTempLog,"\tCompared to Native Structure: %s",StrucTrjArgs[8].strvalue);
		if (chainid!='-')
			BSprintf(bspTempLog," chain %c",chainid);
	}
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
#ifdef OS_UNIX
	/* use uname to get sysinfo */
	if (uname(&utsbuf)>=0) {
		sprintf(pcUname,"%s %s %s %s %s\n",utsbuf.sysname,utsbuf.nodename,utsbuf.release,utsbuf.version,utsbuf.machine);
		BSWrite(bspTempLog,pcUname,StringLen(pcUname));
	}
	else
#else
#ifdef OS_MSWIN
	/* assume ver is available */
	if (GetMSWinVerInfo(pcUname)==ERR_SUCCESS) {
		BSWrite(bspTempLog,pcUname,StringLen(pcUname));
	}	
	else
#endif
#endif


	BSprintf(bspTempLog,"Not Available\n");
	DayTimeStr(timedate,TRUE,TRUE);
	BSprintf(bspTempLog,"Job started: %s\n",timedate);

/* LOG - TITLE BAR */
/* Structure,Time,Tries,BadBB,Crashes,ViolatedConstr,N,Rgyr,HRgyr,NCdist,Rn,Cn,[ASA,HASA,Helix,Edssp,]Ecaca,[RMSD,]Zhang1,VSCORE1,Bryant3,Crease3,VPDBatomN,SolvAtoms1,SdX,SdY,SdZ,Ea,Eb,Ec; */

	BSprintf(bspTempLog,"\nStructure\tTime\tTries\tBadBB\tCrashes\tViolatedConstr\tN\tRgyr\tHRgyr\tNCdist\tRn\tCn");
#ifdef USE_DSSP
	BSprintf(bspTempLog,"\tASA\tHASA");
	BSprintf(bspTempLog,"\tHelix\tEdssp");
#endif
	BSprintf(bspTempLog,"\tEcaca");
	if (pmmdNative!=NULL)
		BSprintf(bspTempLog,"\tRMSD");
	BSprintf(bspTempLog,"\tZhang%d\tVSCORE%d\tBryant%d\tCrease%d\tVPDBatomN",
		ZHANG_WINDOWSIZE,VSCORE_WINDOWSIZE, BRYANT_WINDOWSIZE, BRYANT_WINDOWSIZE);
	if (StrucTrjArgs[15].intvalue > 0) { /* 0 means skip solvation steps and PCA */
	   BSprintf(bspTempLog,"\tSolvAtoms%d\tSdX\tSdY\tSdZ\tEa\tEb\tEc\n",StrucTrjArgs[15].intvalue); /*records number of water layers */
	}

/* LOG - END OF TITLE BAR */

   

	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Opening trajectory distribution database..\n");
	/* sanity check */
	TGInit(tmpdbasename,DB_READ,&dbsize);


/* STRUCTURE MAKING LOOP */
	while (numtofold>0) {
		countstruc++;
		numtofold--;
	        /* load the ASN.1 Biostruc as ASCII ASN.1 file */
	    bspBiostruc=NULL;
	    if (traj_quiet==VERBOSITY_VERBOSE)
			printf("\n_________________\nBeginning Structure #%d..\nReading in chemical graph..\n",countstruc);
		bspBiostruc=MIMEBiostrucAsnGet(tmpskelfname,"r",NULL);
	    if (bspBiostruc==NULL) {
		  	FileRemove(tmpskelfname);
			TGClose();	
	        ErrPostEx(SEV_FATAL,3,1,"Unable to fetch Biostruc");
		        return 3;
	    }
		/* now the chemical graph is correctly in memory as a modelstruc, with
		   co-ordinates (and hence PALDs) assigned only to a-carbons */
	    if (traj_quiet==VERBOSITY_VERBOSE)
			printf("Converting biostruc to modelstruc..\n");
	    pdnmsModelstruc=MakeAModelstruc(bspBiostruc);
	    if (pdnmsModelstruc==NULL) {
	  		CLOSEMMDBAPI();
	        ErrPostEx(SEV_FATAL,4,1,"Unable to convert Biostruc to Modelstruc.");
	        return 4;
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
		/* remove id and descr if necessary */
		/* assign new id and descr */
	    pbs=pmsdHead->pbsBS;
	    if (pbi!=NULL) {
			if (pbs->id!=NULL)
				AsnGenericChoiceSeqOfFree(pbs->id,(AsnOptFreeFunc)BiostrucIdFree);
			pbs->id=pbi;
		}
		if (pbd!=NULL) {
			if (pbs->descr!=NULL)
				AsnGenericChoiceSeqOfFree(pbs->descr,(AsnOptFreeFunc)BiostrucDescrFree);
			pbs->descr=pbd;
		}
		if (countstruc==1)
			firstiter=TRUE;
		else
			firstiter=FALSE;
		if (dbsize!=((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue))->iResCount) {
			PurgeGlobs();
			ErrPostEx(SEV_FATAL,1,1,"protein length inconsistency error, expect: %d actual: %d, possibly due to corrupt trajectory file; aborting",((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue))->iResCount,dbsize);
		}

/*ErrLogPrintf("foldtraj: fold begin\n");*/
		foldtrajprm=(pFoldTrajParamBlock)MemNew(sizeof(FoldTrajParamBlock));
		foldtrajprm->pmsdRoot=pmsdHead;
		foldtrajprm->Model=ModelNum;
		foldtrajprm->err=0;
		foldtrajprm->gen=0;
		foldtrajprm->errorfile= NULL; /* this was for MobiDick reporting , not real erorr reporting */
		trystructure = StrucTrjArgs[6].intvalue-numtofold+StrucTrjArgs[11].intvalue-1;
/*******************now call to randwalk.c *******/
		TRADEProtein((VoidPtr)foldtrajprm,trystructure);
		  /* no longer logging anything about unfinished structures */
		/* if walk completes randwalk.c logs these values: [excluded if DSSP code not linked with USE_DSSP] */
		/* Structure\tTime\tTriest\tBadBB\tCrashes\tViolatedConstr\tN\tRgyr\tHRgyr\tNCdist\tRn\tCn[\tASA\tHASA\tHelix\tEdssp] */
/*************************/
		err=foldtrajprm->err;
		foldtrajprm=MemFree(foldtrajprm);
/*ErrLogPrintf("foldtraj: fold end\n");*/

		if (err==ERR_INCOMPLETE) {
			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("\nTimed out, trying again.\n");
			numtofold++;
			countstruc--;
			incompletecount++;
			/* if get stuck a lot, increae laxness slightly, until it works */
			/* This is causing skew in Rgyr in IDP distributions CWVH June 12 2012 */
			/* continued increase in bounciness causing Rgyr to decrease over large samples */
			/* removing it from default trades behavior - relevant to folding/movie making and even then it is somewhat harmful! */
			if (StrucTrjArgs[19].intvalue == TRUE) {
				if (incompletecount%5==0) {
						ATOM_BOUNCINESS_BB+=0.05;
						ATOM_BOUNCINESS_SC+=0.10;
						BACKBONE_ERROR_TOLERANCE*=1.15;
				}
			}
			pbs->id=NULL;
			pbs->descr=NULL;
			FreeAModelstruc(pdnmsModelstruc);
		}
		else { 

/* REPORT on this structure */

			
			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("Calculating Ecaca extended residues..\n");
			numextres=(int)CalcExtendedResidues((PMMD)(((PMSD)(pdnmsModelstruc->data.ptrvalue))->pdnmmHead->data.ptrvalue),ModelNum);
			BSprintf(bspTempLog,"\t%d",numextres);

			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("Computing Zhang-DeLisi potential..\n");

			pdnListPmmds = CreateListPmmds (NULL, pmsdHead);
			ComputeZhangPotential (&pdnZhangAdjList, pdnListPmmds,
                                   1, FALSE/*bInclusiveWindow*/,
                                   ZHANG_WINDOWSIZE/*iWindowSize*/, piZhangPtnl,
                                   pdnZhangAtmList,FALSE);
			potzh = GetTotalPotential();
			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("Computing Bryant-Lawrence potential..\n");
			ComputeBryantPotential (&pdnBryantAdjList, pdnListPmmds,
                                    1, FALSE/*bInclusiveWindow*/,
                                    BRYANT_WINDOWSIZE/*iWindowSize*/, piBryantPotential,
                                    TRUE/*bUsingPep*/, FALSE);
			potbryant = GetTotalPotential();
			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("Computing VSCORE potential..\n");
			ComputeVscorePotential(((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue)),&pdnVoronoiAdjList,VSCORE_WINDOWSIZE);
			/* returns negative of correct value */
			potVoronoi = GetTotalVPotential(); /* this is old problematic Vscore computed as  -(vTotTerm + vSolvTer ) / vPDBAtomN */
			vPDBatomN = GetPDBatomsVPotential(); /* number of atoms */
			FreeListPmmds (&pdnListPmmds);
			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("Computing Crease Entropy..\n");
			if (CalcCreaseEnergy((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue),UNITS_BRYANT,CreaseIncl,CreaseExcl, FALSE ,FALSE,ModelNum)!=ERR_SUCCESS) {
				ErrPostEx(SEV_ERROR,0,10,"An error occurred in CalcCreaseEnergy");
				potcrease=9999999.9;
			}
			else
	            potcrease=GetCreaseEnergy();
		


/* align to native or example - compute RMSD */
			if (pmmdNative!=NULL) {
				if (traj_quiet==VERBOSITY_VERBOSE)
					printf("Aligning structure to native..\n");
				Align2StrucSVD(pmmdNative,(PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue),
			       ALIGN_CA, 1, pmmdNative->iResCount, NativeModelNum,ModelNum);

				rmsdval = GetRMSD(pmmdNative, (PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue), 
			       ALIGN_CA, 1, pmmdNative->iResCount, NativeModelNum,ModelNum);
			}



/* LOG Values of Scoring Functions */

			if (pmmdNative!=NULL) {
				BSprintf(bspTempLog,"\t%f\t%f\t%f\t%f\t%f\t%d", (float) rmsdval, (float) potzh, (float) potVoronoi, (float) potbryant, (float) potcrease, (int) vPDBatomN);
			}
			else { /* no RMSD */
				BSprintf(bspTempLog,"\t%f\t%f\t%f\t%f\t%d", (float) potzh, (float) potVoronoi, (float) potbryant, (float) potcrease,(int) vPDBatomN);
			}

/* Ramachandran Angle Reporting */ 
			if ((StrucTrjArgs[13].intvalue == TRUE) && (ramalimit <= 30000)) {
				ramalimit++;
				WriteRamaOneModel(pmsdHead,1,trystructure);
			}

/* Calculate & Add Solvent and Rotate via PCA and output Ellipsoidal parameters */
			/* one layer of solvent is needed to calculate the ellipsoidal parameters */
			if (StrucTrjArgs[15].intvalue > 0) {
				SolvateAndPCA(pmsdHead, 1, 1, StrucTrjArgs[18].intvalue, StrucTrjArgs[15].intvalue);
				/* this logs "SolvAtoms1\tSdX\tSdY\tSdZ\tEa\tEb\tEc\t;\n" */
			}
			else { /* finish writing log line without Ellipsoid and Solvent terms */
				BSprintf(bspTempLog,"\t%c\n",LOG_TERM_CHAR);
			}


/* If Unwanted in output, Remove any added solvent from Structures */
			if (((Boolean) StrucTrjArgs[16].intvalue == FALSE) && (StrucTrjArgs[15].intvalue > 0)) {
				RemoveSolvent (pmsdHead);
			}
            

/* WRITE ASN.1 FILE */
			if (StrucTrjArgs[5].intvalue==TRUE) {  
			    if (traj_quiet==VERBOSITY_VERBOSE)
					printf("Writing out ASN.1 structure file..\n");
				StringCpy(fnam2,StrucTrjArgs[10].strvalue ? StrucTrjArgs[10].strvalue : "");
				sprintf(ftmp,"%s%s_%07ld",StrucTrjArgs[10].strvalue ? "_" : "",StrucTrjArgs[1].strvalue ? StrucTrjArgs[1].strvalue:fnamtrjnopath,(long int)(StrucTrjArgs[6].intvalue-numtofold+StrucTrjArgs[11].intvalue-1));
				StringCat(fnam2,ftmp);
				StringCat(fnam2,MMDB_EXT);
				StringCpy(writepath, fnam2);
				/* write out as MIME Biostruc type */
/* This should be edited / replaced with a write to at Bytestore - change to mmdbapi */
				if (!(WriteAsnAllModel(pdnmsModelstruc,writepath,SAVE_BINARY))) {
					CLOSEMMDBAPI();
					TGClose();	
				  	ErrPostEx(SEV_FATAL,1,0,"ASNWrite failed");
					return 6;
				}
				bspTemp=NULL;
				bspTemp=MIMEBiostrucAsnGet(writepath,"rb",NULL);
				FileRemove(writepath);
				/* BINARY ASN.1 */
				if (StrucTrjArgs[14].intvalue == TRUE)
				   aip=AsnIoOpen(writepath,"wb");
                /* ASCII prt goes here*/
				else {
				   StringCpy(fnam2,ftmp);
				   StringCat(fnam2,".prt");
				   StringCpy(writepath,fnam2);
				   aip=AsnIoOpen(writepath,"w");
				}
				decsequence=DecodeSequence(sequence,EXTAA_X);
				nmap=BuildMIMEBiostruc(bspTemp,decsequence,vnpBioseq);
				decsequence=MemFree(decsequence);
				NcbiMimeAsn1AsnWrite(nmap,aip,NULL);
				AsnIoClose(aip);
				/* frees vnpBioseq as a side effect */
				vnpTemp=vnpBioseq;
				vnpHead=NULL;
				while (vnpTemp) {
					if (IS_Bioseq_set(vnpTemp)) {
						bssp=AsnIoMemCopy((BioseqSetPtr)(vnpTemp->data.ptrvalue),(AsnReadFunc)BioseqSetAsnRead,(AsnWriteFunc)BioseqSetAsnWrite);
						ValNodeAddPointer(&vnpHead,2,bssp);
					}
					else if (IS_Bioseq(vnpTemp)) {
						bsp=AsnIoMemCopy((BioseqPtr)(vnpTemp->data.ptrvalue),(AsnReadFunc)BioseqAsnRead,(AsnWriteFunc)BioseqAsnWrite);
						ValNodeAddPointer(&vnpHead,1,bsp);
					}
					vnpTemp=vnpTemp->next;
				}
				NcbiMimeAsn1Free(nmap);
				vnpBioseq=vnpHead;
			}

/* WRITE PDB FILE */

			if (StrucTrjArgs[4].intvalue==TRUE) {
			    if (traj_quiet==VERBOSITY_VERBOSE)
					printf("Writing out PDB structure file..\n");
				StringCpy(fnam2,StrucTrjArgs[10].strvalue ? StrucTrjArgs[10].strvalue : "");
				sprintf(ftmp,"%s%s_%07ld",StrucTrjArgs[10].strvalue ? "_" : "",StrucTrjArgs[1].strvalue ? StrucTrjArgs[1].strvalue:fnamtrjnopath,(long int)(StrucTrjArgs[6].intvalue-numtofold+StrucTrjArgs[11].intvalue-1));
				StringCat(fnam2,ftmp);
				StringCat(fnam2,PDB_EXT);
				StringCpy(writepath, fnam2);
				fp=FileOpen(writepath,"w");
			    if (!(WritePDBAllModel(pdnmsModelstruc,fp))) {
 					FileClose(fp);
					CLOSEMMDBAPI();
					TGClose();	
			        ErrPostEx(SEV_FATAL,1,0,"PDBWrite failed");
					return 6;
				}
				FileClose(fp);
			}

/* Free Good  Structures */
			/* trick freer so pbi and pbd are saved for next run */
			pbs->id=NULL;
			pbs->descr=NULL;
			FreeAModelstruc(pdnmsModelstruc);
			FreeAdjList (&pdnZhangAdjList);
			FreeAdjList (&pdnBryantAdjList);
			FreeAdjList (&pdnVoronoiAdjList);
			FreeCreaseEnergy();

/* flush bytestore every 1000 structures */
			writecount++;
			if (writecount == 1000) {
				if (traj_quiet==VERBOSITY_VERBOSE)
					printf("Writing cached log lines - batch of 1000\n");
				writecount = 1;
				ByteStoreFileOut( bspTempLog, LogOutName);
				BSFree(bspTempLog);
				bspTempLog = BSNew(BS_STARTSIZE);
				if ((StrucTrjArgs[13].intvalue == TRUE) && (ramalimit <= 30000)) {
					/* flush the Ramachandran bytestores too */
					
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

		} /* End of reporting section for good structures */	 



	}  
/* END OF STRUCTURE MAKING LOOP */

	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Finishing Log file %s\n",LogOutName);
	BSprintf(bspTempLog,"----------------------------------------\n\n");
	ByteStoreFileOut( bspTempLog, LogOutName);
	BSFree(bspTempLog);

	if (StrucTrjArgs[13].intvalue == TRUE) {
		/* flush the remaining Ramachandran bytestores and free them */
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
	}			

/*	DumpLog(LogOutName,bspTempLog,FALSE); */
	
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Final Cleanup..\n");
	if (pbi!=NULL)
		AsnGenericChoiceSeqOfFree(pbi,(AsnOptFreeFunc)BiostrucIdFree);
	if (pbd!=NULL)
		AsnGenericChoiceSeqOfFree(pbd,(AsnOptFreeFunc)BiostrucDescrFree);
	FileRemove(tmpskelfname);
	TGClose();	
	CleanUpDB(tmpdbasename);
	if (vnpBioseq!=NULL) {
		vnpTemp=vnpBioseq;
		while (vnpTemp) {
			if (IS_Bioseq_set(vnpTemp))
				BioseqSetFree((BioseqSetPtr)(vnpTemp->data.ptrvalue));
			else if (IS_Bioseq(vnpTemp))
				BioseqFree((BioseqPtr)(vnpTemp->data.ptrvalue));
			vnpTemp=vnpTemp->next;
		}
		ValNodeFree(vnpBioseq);
	}	
	 
	FreePNNList(); /* distance constraint list */
 	FreeRotLib();
    FreeZhangAtmList (&pdnZhangAtmList);
    FreeZhangPotential (&piZhangPtnl);
	FreeBryantPotential (&piBryantPotential);
	
    
	CLOSEMMDBAPI();
 	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Done.\n");
	return 0; /* successful program completion */

}

