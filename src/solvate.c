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
#include <tslri.h>
#include <mmdbapi.h>



const int BS_STARTSIZE = 5000;


/* Global Variables */
Args Rargs[10] = {
/*0*/ {"Input Asn.1 3d Structure File Name:",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
/*1*/ {"Scoring Function, Solvation, PCA Output (T/f):","TRUE",NULL,NULL,TRUE,'s',ARG_BOOLEAN,0.0,0,NULL},
/*2*/ {"Ramachandran Angle Output (T/f):","TRUE",NULL,NULL,TRUE,'r',ARG_BOOLEAN,0.0,0,NULL},
/*3*/ {"Write Solvated ASN.1 BINARY Biostruc (t/F):","FALSE",NULL,NULL,TRUE,'b',ARG_BOOLEAN,0.0,0,NULL},
/*4*/ {"Write Solvated ASN.1 ASCII Biostruc (t/F):","FALSE",NULL,NULL,TRUE,'a',ARG_BOOLEAN,0.0,0,NULL},
/*5*/ {"Write Solvated PDB File (t/F):","FALSE",NULL,NULL,TRUE,'p',ARG_BOOLEAN,0.0,0,NULL},
/*6*/ {"Omit Hydrophobic Core Solvent:","TRUE",NULL,NULL,TRUE,'h',ARG_BOOLEAN,0.0,0,NULL},
/*7*/ {"(Reserved for future use)","FALSE",NULL,NULL,TRUE,'n',ARG_BOOLEAN,0.0,0,NULL},
/*8*/ {"Number of Water Layers to Add (1-20)","1","1","20",TRUE,'l',ARG_INT,0.0,0,NULL},
/*9*/ {"Do PCA coordinate transformation: (T/f) (false leaves coordinates intact) ","TRUE",NULL,NULL,TRUE,'m',ARG_BOOLEAN,0.0,0,NULL},
};
             
/*7  Output files with NO solvent molecules (strips crystallographic water) */

ByteStorePtr bspA;
ByteStorePtr bspPotential;


Int4Ptr piBryantPotential = NULL;
Int4Ptr piZhangPtnl = NULL;
DValNodePtr pdnZhangAtmList = NULL;
ValNodePtr pvnSolventList = NULL;
ValNodePtr pvnCur = NULL;


#define ZHANG_WINDOWSIZE 1   /* Reset to 1 CWVH 2012 - Small Peptide scores are off. Entropy Calculation Calibration */ 
#define BRYANT_WINDOWSIZE 3
#define VSCORE_WINDOWSIZE 1


/* CWVH Solvation Globals */
static int iAtomCount = 0;
static int iGraphCount = 0;
static int iWaterCount = 0;
static DValNodePtr pdnmmEnd = NULL;
static Int4 iLocalDict = 0;
static ResidueGraphPtr prgWat = NULL; /* pointer to the Water Chemical Graph */


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



Int4 SolvateIt(PMAD pmadAtom, ValNodePtr pvnAtomList, DValNodePtr pdnListPmmds, Int2 iModel, PWS pwsHead, Boolean bSkipInternal, FloatLo flBFactor )
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
	  pflvWat[4] = (float) flBFactor; /* maximum B factor */
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


Int4Ptr ComputeSolvation (DValNodePtr pdnListPmmds, Int2 iModelNum, Boolean bSkipInternal, Int4 iLayers) {

ValNodePtr  pvnAtomList;
DValNodePtr pdnHere;
DValNodePtr pdnFirstSolvent = NULL;
DValNodePtr pdnLastSolvent = NULL;



PWS pwsHead = NULL;
PWS pwsThis = NULL;
vec middle;
Int2  numAtoms;
Int4 iTotalWater = 0;
Int4Ptr pI4LayerArray;
int i;

PMMD pmmd1;
PMGD pmgd1;
PMAD pmad1;
PALD pald1;

PDNMG pdnmgHere;
PVNMA pvnmaHere;

    /* For each molecule in the list */
Int4 iThisLayer = 0;
FloatLo flBFactorS = 99.00;
FloatLo flBFactor = 0.0;
FloatLo flLayer = 0.0;

pI4LayerArray = (Int4Ptr) MemNew((size_t)sizeof(Int4) * iLayers+1); /* total goes in 0 element */

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
	iThisLayer = 0;
	flLayer = (FloatLo) i;
	flBFactor = flBFactorS + (flLayer/100);

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
							    iThisLayer += SolvateIt(pmad1, pvnAtomList, pdnListPmmds, iModelNum, pwsHead, bSkipInternal, flBFactor );
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

iTotalWater += iThisLayer;
pI4LayerArray[i]=iThisLayer;
iThisLayer = 0;
pI4LayerArray[0] = iTotalWater;

} /* for each solvation layer to be added */



FreeAllWorlds ();



return pI4LayerArray;
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




void WritePotentialsOneModel(PMSD pmsdHead, Int2 ModelNum, Int4 structureNo, CharPtr filein, Boolean DumpSolvatedASNBinary, 
							 Boolean DumpSolvatedASNAscii, Boolean DumpSolvatedPDB, Boolean SkipCoreSolvent, 
							 Boolean NoSolventIn3DFile, Int4 iLayers, Boolean DoTransform) {
  DValNodePtr pdnListPmmds = NULL;
  DValNodePtr pdnZhangAdjList = NULL;
  DValNodePtr pdnBryantAdjList = NULL;
  DValNodePtr pdnVoronoiAdjList = NULL;
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

  PMMD pmmdHere;
  FloatHi potzh1=0.0,  potbryant3 = 0.0, potcrease3=0.0, potVoronoiF = 0.0,  vTotTerm = 0.0, vSolvTerm = 0.0, vSSTerm = 0.0, vASATotal = 0.0;
  Int4 vPDBatomN = 0;
  Int4 linelen = 0;
  Int2 CreaseIncl = 0;
  Char line[2000];
  Char tempstring[200];
  FILE *PDBFile = NULL;

  /* variables for removal of ppAsnOrder */
  PMLD pmldTo = NULL;
  PDNML pdnmlTo = NULL;
  Int4Ptr pI4WaterArray = NULL;

  /* printf("Computing Potentials..\n"); */
  pdnListPmmds = CreateListPmmds (NULL, pmsdHead);
  /* this is not the linked list for the entire structure - it is just a linked list with the FIRST molecule in it ! */
  if (pdnZhangAtmList == NULL) printf("NULL pdnZhangAtmList\n");
  printf("Zhang Potential Start \n");
  ComputeZhangPotential (&pdnZhangAdjList, pdnListPmmds,
                                   1, FALSE,
                                   1 , piZhangPtnl,
                                   pdnZhangAtmList,FALSE);
  potzh1 = GetTotalPotential();
  FreeAdjList (&pdnZhangAdjList);


  printf("Bryant Potential Start \n");
/* using BRYANT_WINDOWSIZE 3 */
  ComputeBryantPotential (&pdnBryantAdjList, pdnListPmmds,
                                    1, FALSE/*bInclusiveWindow*/,
                                    /*BRYANT_WINDOWSIZE*/3, piBryantPotential,
                                    TRUE/*bUsingPep*/, FALSE);
  potbryant3 = GetTotalPotential();
  FreeAdjList (&pdnBryantAdjList);
  ComputeVscorePotential((PMMD) pdnListPmmds->data.ptrvalue,&pdnVoronoiAdjList,VSCORE_WINDOWSIZE);
			/* returns negative of correct value */
  /* potVoronoi = GetTotalVPotential(); */ /* this is old problematic Vscore computed as  -(vTotTerm + vSolvTer ) / vPDBAtomN */
  vTotTerm =  GetTotTermVPotential();  /* Get the 3 raw score terms contact-contact */
  vSolvTerm = GetSolvTermVPotential(); /* contact-solvent */ 
  vSSTerm = GetSSTermVPotential();  /* solvent - solvent */
  vPDBatomN = GetPDBatomsVPotential(); /* number of atoms */
  vASATotal = GetASATotal(); /* Voronoi calculated ASA  - can use this to replace the DSSP calculation... */
  potVoronoiF = GetTotalVPotential(); /* Terms as per original TRADES implementation */
  FreeAdjList (&pdnVoronoiAdjList); /* Voronoi adjacency list does not include pairwise res-res contacts... */


/* cutoff 15, exclusive window size 3 */
/* FoldTraj CreaseExcl = 3  - differs from BL-Score so that it includes helical interactions!  */ 
/* CreaseIncl = NUMBER of AMINO ACIDS */
/* Using CREASE_WINDOWSIZE 3*/
  pmmdHere = (PMMD) pdnListPmmds->data.ptrvalue;
  CreaseIncl = pmmdHere->iResCount;
	  if (CalcCreaseEnergy((PMMD) pdnListPmmds->data.ptrvalue,UNITS_BRYANT,CreaseIncl,3/*CREASE_WINDOWSIZE*/,FALSE,FALSE,ModelNum)!=ERR_SUCCESS) {
		ErrPostEx(SEV_ERROR,0,10,"An error occurred in CalcCreaseEnergy");
		potcrease3=9999999.9;
  }
	else
	    potcrease3=GetCreaseEnergy();
  FreeCreaseEnergy();
  FreeListPmmds(&pdnListPmmds);


  pI4WaterArray = ComputeSolvation(pmsdHead->pdnmmHead, ModelNum,SkipCoreSolvent,iLayers); /* TRUE is to skip internal solvation */
  iWater = pI4WaterArray[0];
  printf("Solvated structure with %d water molecules\n",iWater);
  /* traverse the structure and pull out all the AM_SOL PMADs into a linked list with PMADs replacing the Voronoi list with Solvated water list */
  /* this adds to the global pointer pvnSolventList */
  /* compute PCA on Solvated Atom List */
  printf("PCA Start \n");
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
  printf("Counted %d Solvent Atoms \n",(int) iSolvAtoms); 
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


if (NoSolventIn3DFile) {
 /* clean up code - MMDBAPI should do this by itself, but checking to see if any leakage... */
		/* prune internal waters  - Remove ppAsnOrder, Molecule nodes pdnmm and dictionary entry */
		/* printf("Start of Solvent Removal\n"); */
		/* Remove any vectors of ppAsnOrder made from writing out ASN.1  */
  	 pdnmlTo = pmsdHead->pdnmlModels;
		if (pdnmlTo != NULL) {
		      pmldTo = (PMLD)(pdnmlTo->data.ptrvalue);
			  if (pmldTo != NULL)
		         if(pmldTo->ppAsnOrder) {
			       PTRVectorFree(pmldTo->ppAsnOrder,0);
			       pmldTo->ppAsnOrder = NULL;
				 } 		
		}
		DoWipeSolvent(pmsdHead);  /* If present, this wipes out x-ray but leaves structure in an unusable state for
			internal indexing problem in  MMDBAPI but Cn3D is OK.
			 {mmdbapi1.c, line 5092} Biostruc ASN.1 Internal Indexing Failure  */
		RemoveWaterDictionary(pmsdHead); 
		/* these calls must be done together to remove locally added water */ 	
		/* printf("WaterDictionary Removed\n"); */
	}


if ((DumpSolvatedPDB == TRUE) || (DumpSolvatedASNAscii == TRUE) || (DumpSolvatedASNBinary)) { /* */
/* Now that we have the mean and the rotation matrix, transform the entire coordinate set */
	if (DoTransform == TRUE) {
    TraverseAtoms(pmsdHead->pdnmmHead, ModelNum, 0, (Pointer) pfmean, (pNodeFunc) (DoPCACenter));
    TraverseAtoms(pmsdHead->pdnmmHead, ModelNum, 0, (Pointer) ppfSTM, (pNodeFunc) (DoPCARotate));
/*	printf("Calling TraverseSolids\n"); */
	TraverseSolids(pmsdHead->pdnmsLink, ModelNum, 0, (Pointer) pfmean, (pNodeFunc) (DoPCACenterSolids));
	TraverseSolids(pmsdHead->pdnmsLink, ModelNum, 0, (Pointer) ppfSTM, (pNodeFunc) (DoPCARotateSolids));
	}
/* Write out the - now axis aligned - structure */
	if (DumpSolvatedPDB == TRUE) {
		StringCpy(line, filein);
		StringCat(line,"_solv_r.pdb");	
		PDBFile = FileOpen(line,"w");
		printf("Writing Rotated PDB: %s\n",line);
		WritePDBOneModel(pmsdHead->pdnmsLink, PDBFile, ModelNum); 
		fflush(PDBFile);
		FileClose(PDBFile);
	}
	if (DumpSolvatedASNAscii == TRUE) {
	StringCpy(line, filein);
	StringCat(line,"_solv_r.prt");
	printf("Writing ASN.1 ASCII with solvent: %s\n",line); 
	WriteAsnOneModel(pmsdHead->pdnmsLink,ModelNum,line,SAVE_ASCII);
	}
	if (DumpSolvatedASNBinary == TRUE) {
	StringCpy(line, filein);
	StringCat(line,"_solv_r.val");
	printf("Writing ASN.1 Binary with solvent: %s\n",line); 
	WriteAsnOneModel(pmsdHead->pdnmsLink,ModelNum,line,SAVE_ASCII);
	}

}




pvnHere = NULL;

line[0] = '\0';
StringCpy(line,",H20_Added");
for (i = 1; i<=iLayers; i++) {
	sprintf(tempstring,",H20_L%d",i);
	StringCat(line,tempstring);
}

StringCat(line,"\n");
linelen = StringLen(line);
BSWrite(bspPotential, line, linelen);

sprintf(line,"\n%s.val,%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f", 
	 filein, (float) potzh1, (float) potbryant3, (float) potcrease3, (float) potVoronoiF, 
	 (float) vTotTerm, (float) vSolvTerm, (float) vSSTerm, (int) vPDBatomN, (float) vASATotal, (int) iSolvAtoms,
	 (float) pfstddev[0],(float) pfstddev[1], (float) pfstddev[2],
     (float) pfEigen[2], (float) pfEigen[1], (float) pfEigen[0]); /*, LOG_TERM_CHAR); */

linelen = StringLen(line);
BSWrite(bspPotential, line, linelen);

line[0] = '\0';
for (i = 0; i<=iLayers; i++) {
	sprintf(tempstring,",%d",pI4WaterArray[i]);
	StringCat(line,tempstring);
}
StringCat(line,"\n");
linelen = StringLen(line);
BSWrite(bspPotential, line, linelen);

MemFree(pI4WaterArray);



/*
sprintf(line,"%s.val,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%c\n", 
	 filein, (float) potzh1, (float) potbryant4, (float) potbryant3, (float) potcrease4, (float) potVoronoiF, 
	 (float) vTotTerm, (float) vSolvTerm, (float) vSSTerm, (int) vPDBatomN, (float) vASATotal, (int) iSolvAtoms,
	 (float) pfmean[0], (float) pfstddev[0],(float) pfmean[1],(float) pfstddev[1],(float) pfmean[2],(float) pfstddev[2],
     (float) pfEigen[2], (float) pfEigen[1], (float) pfEigen[0], 
     (float) ppfSTM[0][2], (float) ppfSTM[0][1], (float) ppfSTM[0][0], 
     (float) ppfSTM[1][2], (float) ppfSTM[1][1], (float) ppfSTM[1][0],
     (float) ppfSTM[2][2], (float) ppfSTM[2][1], (float) ppfSTM[2][0], LOG_TERM_CHAR);
*/




/*    printf("Start of PCA data structure removal \n"); */

MATH_VectorFree(&pfmean);
MATH_VectorFree(&pfstddev);
MATH_MatrixFree(&ppfCov,3);
MATH_VectorFree(&pfEigen);
MATH_MatrixFree(&ppfSymm,3);
MATH_MatrixFree(&ppfSTM,3);

 } /* if pvnSolventList */
} /* Write Potentials One Model */




void WriteRamaOneModel(PMSD pmsdRoot, Int2 ModelNum ){

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
	sprintf(line,"%s, %d, %f, %f\n", aa , (int) num, phi,psi);
    linelen = StringLen(line);
    BSWrite(bspA, line, linelen);
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
  Int2    iDotLen = 0;

/* Variables for clipping status and pdb output file name*/
    static char fileout[PATH_MAX];
	static char filein[PATH_MAX];
	static char txtnum[10];
    Int2    iTest = 0;
    FILE    *pFile = NULL;

/* Variables for potentials */
	FILE *fp, *fpAtmInfo=NULL;
	static char header[] = "Filename,Zhang1,Bryant3,Crease3,VoronoiF,VTotTerm,VSolvTerm,VSSTerm,VPDBatomN,VASATotal,SolvAtoms,SdX,SdY,SdZ,Ea,Eb,Ec";
    int retval = 0;

    
/* Initialize MMDB-API */
        ErrSetLogfile("error_solvate.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);

	if (!GetArgs("solvate  ASN.1 Structure Utility that Computes Phi, Psi, & Potential Scores,\n Adds solvent layers and rotates ASN.1 structures along ellipsoid axes\n",10,Rargs))
		return 1;
    StringCpy(fileout, "__");
	StringCat(fileout,Rargs[0].strvalue);
	if(StringLen(fileout) == 2) {
	  ErrPostEx(SEV_FATAL,2,1,"No filename provided");
	  return 3;
	}
    StringCat(fileout,".csv");
	
	if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI, check for bstdt.val");
		return 4;
	}
	
/* Intitialize potentials */


	if ((fp = FileOpen (ZHANG_POTENTIAL, "r")) == NULL) {
	    ErrPostEx (SEV_FATAL, 1, 7, "File not found (Zhang Potential) %s.",ZHANG_POTENTIAL);
	    return 7;
	}

	if ((fpAtmInfo = FileOpen (ZHANG_ATOMS, "r")) == NULL) {
	    ErrPostEx (SEV_FATAL, 1, 7, "File not found (Zhang AtmList) %s.",ZHANG_ATOMS);
	    return 8;
	}

	if (LoadZhangPotential (&piZhangPtnl, fp)!=ERR_SUCCESS) {
	    ErrPostEx (SEV_FATAL, 1, 9, "Error in LoadZhangPotential");
	    return 9;
	}
	if (LoadZhangAtmList (&pdnZhangAtmList, fpAtmInfo)!=ERR_SUCCESS) {
	    ErrPostEx (SEV_FATAL, 1, 9, "Error in LoadZhangAtomList");
	    return 9;
	}
    FileClose (fp);
 	FileClose (fpAtmInfo);
 

	if ((fp = FileOpen (BL_POTENTIAL, "r")) == NULL) {
	    ErrPostEx (SEV_FATAL, 1, 8, "File not found (BL potential)%s.", BL_POTENTIAL);
	    return 10;
	}

	if (LoadBryantPotential (&piBryantPotential, fp, TRUE /*bUsingPep*/, FALSE /*bNegateInput*/)!=ERR_SUCCESS) {
	    ErrPostEx (SEV_FATAL, 1, 10, "Error in LoadBryantPotential");
	    return 11;
	}
    FileClose (fp);

	bspA = BSNew(BS_STARTSIZE);
	bspPotential = BSNew(BS_STARTSIZE);

	if (Rargs[3].intvalue == TRUE) {
		 printf("Computing Potentials..\n"); 
		/* CSV file header */
		BSWrite(bspPotential,header,StringLen(header));
	}


		
     /* keep any extension typed in */
     iDotLen = 0;
     filein[0] = '\0';
     iDotLen = (Int2) StringLen(Rargs[0].strvalue);
     if (iDotLen > 5) {
		if (Rargs[0].strvalue[iDotLen - 4] == '.') { /* extension typed, use it first one we hit */
		     StringCpy(filein,Rargs[0].strvalue);
		     if (FileLength(filein) == 0)  {
		      		 ErrPostEx (SEV_FATAL, 1, 10, "Unable to find file  %s",Rargs[0].strvalue);
	   			 return 10;
		    }	              
		}
     }
     if (filein[0] == '\0') {
        StringCpy(filein,Rargs[0].strvalue);
	StringCat(filein,".val");
	if (FileLength(filein) == 0) {
		     StringCpy(filein,Rargs[0].strvalue);	              
             	     StringCat(filein,".cn3");
			if (FileLength(filein) == 0) {
		     	   StringCpy(filein,Rargs[0].strvalue);	              
             	           StringCat(filein,".prt");
			  if (FileLength(filein) == 0) {
		      		 ErrPostEx (SEV_FATAL, 1, 11, "Unable to find file with .val .prt .cn3 extension %s",Rargs[0].strvalue);
	   			 return 11;
			  }
           		}
	}
    }
    


	
	   pmsdRoot=LoadABiostruc(filein,FALSE,0,&ModelNum);
	   if (pmsdRoot==NULL) {
	   	 printf("Unable to load %s, exiting\n",filein);
		 ClearStructures();
		 BSFree(bspA);
  		 BSFree(bspPotential);
		 ErrPostEx (SEV_FATAL, 1, 12, "Unable to Load Biostruc named: %s",filein);
	   			 return 12;
		 return 1;
	   }
	   if (!(pmsdRoot->bWhat & AM_PROT)) {
			ErrPostEx (SEV_WARNING, 1, 12, "No protein proceeding with reckless abandon on file: %s",filein);
	   }
	  printf("%s Loaded\n",filein);
	  
      StringCpy(filein,Rargs[0].strvalue);
	
      if (Rargs[2].intvalue == TRUE) WriteRamaOneModel(pmsdRoot, ModelNum);
	  if (Rargs[1].intvalue == TRUE) WritePotentialsOneModel(pmsdRoot, ModelNum, 1, filein, (Boolean) Rargs[3].intvalue, (Boolean) Rargs[4].intvalue, 
		  (Boolean) Rargs[5].intvalue, (Boolean) Rargs[6].intvalue, (Boolean) Rargs[7].intvalue, (Int4) Rargs[8].intvalue, (Boolean) Rargs[9].intvalue);
	  
	  fileout[0] = 'A';
	  ByteStoreFileOut( bspA, fileout);
      printf("Output file %s has Ramachandran angles\n",fileout);

	  fileout[0] = 'Z';
	  ByteStoreFileOut( bspPotential, fileout);
	  printf("Output file %s has Potential, Solvation and PCA terms.\n",fileout);



	
/* Shut Down MMDB-API */
/* All Modelstrucs remaining are freed in CloseMMDB-API() */
	CloseMMDBAPI();	


/* Write out all the files remaining */


	BSFree(bspA);
	BSFree(bspPotential);

/* Free Potentials */
	FreeZhangAtmList (&pdnZhangAtmList);
	FreeZhangPotential (&piZhangPtnl);
 	FreeBryantPotential (&piBryantPotential);

 	return retval;
}
