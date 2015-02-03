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


/**
 *
 *     FILE NAME: potential.c
 *        AUTHOR: Mark Kotowycz
 * CREATION DATE: May 25, 1999
 *
 * DESCRIPTION -- Contains functions that load a potential table into
 *                memory and perform various operations on it.
 *
 */

#include <mmdbtraj.h>
#include "potential.h"


/* For getting a 2D, 3D OR 4D array value.  Pass in variables, and the last (Z) is an
   array with the dimensions of the multidimentional array */
#define Get2DArrayValue(A,B,Z)     (Int4) ((B)*(Z[0])+(A))
#define Get3DArrayValue(A,B,C,Z)   (Int4) ((C)*(Z[1])*(Z[0])+(B)*(Z[0])+(A))
#define Get4DArrayValue(A,B,C,D,Z) (Int4) ((D)*(Z[2])*(Z[1])*(Z[0])+(C)*(Z[1])*(Z[0])+(B)*(Z[0])+(A))


/* distance bins in C-Beta lookup table */
#define NUMBINS 5
/* index into length and angle arrays for Cis-Proline */
#define CISPRO 27
/* index into several length and angle arrays for Cystine */
#define CYSSYC 26		
/* index into length and angle arrays for Cis-Proline */
#define CISPRO 27

static CharPtr NCBIstdaaUC = "-ABCDEFGHIKLMNPQRSTVWXYZU*";

/* global look-up tables */
static Boolean gCBTable_IsLoaded = FALSE;
static FloatLo cbdir[30][NUMBINS+1][3];
static FloatLo cbcisdir[30][3];
static vec vZero={0.0,0.0,0.0};
static FloatLo bl_cacb[]={1.530F,0.0F,1.528F,1.533F,1.531F,1.534F,1.070F,1.542F,1.554F,1.528F,1.536F,1.528F,1.534F,1.527F,1.529F,1.532F,1.530F,1.560F,1.540F,1.534F,0.0F,1.541F,0.0F,0.0F,0.0F,0.0F,1.536F};
/*static FloatLo coszeta[30];
static FloatLo coseta[30];
static FloatLo bl_cao[30];
static FloatLo ba_ocaca[30];*/

/* Static global variables - used internally */
static FloatHi TotalPotential;
static DValNodePtr pdnListPepMGDs;
static DValNodePtr pdnListPepMGDsBak;

void BackupPepMGDs(void)
{
	pdnListPepMGDsBak=pdnListPepMGDs;
}

void RestorePepMGDs(void)
{
	pdnListPepMGDs=pdnListPepMGDsBak;
}

/**
 *
 * GetIndexEx
 *
 * Takes a string "str" and returns the corresponding value that
 * can be used to index the BL potential function array. Returns -1 if
 * "str" is not a valid representation of an amino acid.
 *
 */

Int2 GetIndexEx (CharPtr str, Boolean bAlphaIndex) {

    Int2 i, n = (Int2) StringLen (str);
	CharPtr pcAA1 = NULL;
	CharPtr pcAA3 = NULL;

	if(bAlphaIndex) {
		pcAA1 = BL_AA_1LETTER;
		pcAA3 = BL_AA_3LETTER;
	} else {
		pcAA1 = MK_AA_1LETTER;
		pcAA3 = MK_AA_3LETTER;
	}

	if (n == 1) {
		for (i=0; i < NUM_AA; i++) {
			if (StringNICmp (str, &(pcAA1[i]), 1) == 0) {
				return i;
			}
		}
    } else if (n == 3) {
		for (i=0; i < NUM_AA; i++) {
			if (StringNICmp (str, &(pcAA3[i*3]), 3) == 0) {
				return i;
			}
		}
    }

    return -1;

}

/* Get the index using the hydrophobic-hydrophilic arrangement of letters */
Int2 GetIndex (CharPtr str) {
	Boolean bAlphaIndex = FALSE;
	return GetIndexEx (str, bAlphaIndex);
}

/* Gets the index using alphabetical arrangement of letters */
Int2 GetAlphaIndex (CharPtr str) {
	Boolean bAlphaIndex = TRUE;
	return GetIndexEx (str, bAlphaIndex);
}



static Int2 GetBLFile(CharPtr pcID, FILE* PNTR pfp, Boolean bOpen)
{
	FILE *fp = NULL;
	CharPtr pcFileName = NULL;
	
	pcFileName = MemNew((size_t) sizeof(Char) *(StringLen(pcID)+1)*(StringLen(CFG_local_datafilepath)+1)*PATH_MAX);
	if(pcFileName == NULL) {
		ErrPostEx (SEV_FATAL, 0, 1, "Could not open file %s\n.", pcFileName);
		return 1;
	}

	sprintf(pcFileName,"%s%s",CFG_local_datafilepath,pcID);
	if(bOpen) {
		if ((fp = FileOpen (pcFileName, "r")) == NULL) {
			ErrPostEx (SEV_FATAL, 0, 1, "Could not open file %s\n.", pcFileName);
			pcFileName = MemFree(pcFileName);
			return 1;
		}
	} else { 
		if ((fp = FileOpen (pcFileName, "w+")) == NULL) {
			ErrPostEx (SEV_FATAL, 0, 1, "Could not create file %s\n.", pcFileName);
			pcFileName = MemFree(pcFileName);
			return 1;
		}
	}
	MemFree(pcFileName);
	*pfp = fp;
	return 0;
}


/* reads in table of Cbeta locations */
static TrajErr LoadCBTable(void)
{
	FILE *f;
	Char buf[255];
	Int2 res=0;
	Char rname[255];
	float  cbtemp [NUMBINS+1][3];
	float  costemptze;
	float  costempet;
	float  bl_caotemp;
	float  ba_ocacatemp;
	float  cbcisdirtemp[3];

	sprintf(buf,"%s%s",CFG_local_datafilepath,CBFNAME);
	if ((f=FileOpen(buf,"r"))==NULL) {
			ErrPostEx(SEV_ERROR, 0, 0, "Could not load in cbdata file.\n");
            return ERR_FAIL;
    }
    while (FileGets(buf,255,f)!=NULL) {
		sscanf(buf,"%s",rname);
		if (rname[0]=='#') continue;
		/* place these after IUPAC definition */
		if (!StringCmp(rname,"CX")) res=CYSSYC;
		else if (!StringCmp(rname,"PC")) res=CISPRO;
		else if (StringLen(rname)>1) continue;
		if (StringLen(rname)==1) res= (Int2) StringCSpn(NCBIstdaaUC,rname);
        
		sscanf(buf,"%*s %*s %f %f %f %f %f %f %f %f %f %f %f %f %f%f %f %f %f %f %f %f %f %f %f %f %f",
	&cbtemp[0][0],&cbtemp[1][0],&cbtemp[2][0],&cbtemp[3][0],&cbtemp[4][0],&cbtemp[5][0],&cbtemp[0][1],
	&cbtemp[1][1],&cbtemp[2][1],&cbtemp[3][1],&cbtemp[4][1],&cbtemp[5][1],&cbtemp[0][2],
	&cbtemp[1][2],&cbtemp[2][2],&cbtemp[3][2],&cbtemp[4][2],&cbtemp[5][2],
            &costemptze,&costempet,&bl_caotemp,&ba_ocacatemp,   
            &cbcisdirtemp[0],&cbcisdirtemp[1],&cbcisdirtemp[2]);
               cbdir[res][0][0] = (FloatLo) cbtemp[0][0];  
               cbdir[res][1][0] = (FloatLo) cbtemp[1][0];
               cbdir[res][2][0] = (FloatLo) cbtemp[2][0];
               cbdir[res][3][0] = (FloatLo) cbtemp[3][0];
               cbdir[res][4][0] = (FloatLo) cbtemp[4][0];
               cbdir[res][5][0] = (FloatLo) cbtemp[5][0];
               cbdir[res][0][1] = (FloatLo) cbtemp[0][1];
               cbdir[res][1][1] = (FloatLo) cbtemp[1][1];
               cbdir[res][2][1] = (FloatLo) cbtemp[2][1];
               cbdir[res][3][1] = (FloatLo) cbtemp[3][1];
               cbdir[res][4][1] = (FloatLo) cbtemp[4][1];
               cbdir[res][5][1] = (FloatLo) cbtemp[5][1];
               cbdir[res][0][2] = (FloatLo) cbtemp[0][2];   
               cbdir[res][1][2] = (FloatLo) cbtemp[1][2];   
               cbdir[res][2][2] = (FloatLo) cbtemp[2][2];   
               cbdir[res][3][2] = (FloatLo) cbtemp[3][2];
               cbdir[res][4][2] = (FloatLo) cbtemp[4][2];
               cbdir[res][5][2] = (FloatLo) cbtemp[5][2];
/*               coszeta[res]     = (FloatLo) costemptze;
               coseta[res]      = (FloatLo) costempet;
               bl_cao[res]      = (FloatLo) bl_caotemp;
               ba_ocaca[res]    = (FloatLo) ba_ocacatemp;*/
               cbcisdir[res][0] = (FloatLo) cbcisdirtemp[0];
               cbcisdir[res][1] = (FloatLo) cbcisdirtemp[1];
               cbcisdir[res][2] = (FloatLo) cbcisdirtemp[2];
	}
    FileClose(f);
	return ERR_SUCCESS;
}


/* given U basis vectors and two adjacent CAs and residue type, returns CB direction
   in vCBRef (relative to U basis vectors) and vCBHere (relative to molecular
   co-ordinate system) */
/************** requires that LoadCBTable has successfully executed ****************/
static void FindCBDir(vec uone,vec utwo,vec uthree,vec vCALast, vec vCANext,Int2 res,Int2 resNext,vec vCBRef, vec vCBHere, FloatLo CBOffset)
{
	vec vDist,vTmp,vNoise;
	FloatLo dist;  /*,noisemag,noisephi,noisetheta;*/
	Int2 distbin,cnt;
	FloatHi dbin[NUMBINS]={5.1,5.6,6.1,6.6,7.0};

        /* dist==d in Table I of Rey and Skolnick, J. Comp. Chem.,  */
	VecSub(vDist,vCANext,vCALast);  
	dist=getMag(vDist);
		/* distance bins range from 0 to NUMBINS */
	distbin=NUMBINS;
	for (cnt=0;cnt<NUMBINS;cnt++) {
		if (dist<dbin[cnt]) {
			distbin=cnt;
			break;
		}
	}
        /* use look-up table to get cosd1,cosd2,cosd3 */
	/* this strange mapping accounts for the error in Table I of the
	   Rey and Skolnick paper */
	if (resNext!=CISPRO) {
		vCBRef[0]=-cbdir[res][distbin][2];
		vCBRef[1]=cbdir[res][distbin][0];
		vCBRef[2]=-cbdir[res][distbin][1];
	}
	else {
		vCBRef[0]=cbcisdir[res][0];
		vCBRef[1]=cbcisdir[res][1];
		vCBRef[2]=cbcisdir[res][2];
	}
		/* add some random noise */
	vNoise[0]=(FloatLo)0.05*Rand1Distrib();
	vNoise[1]=(FloatLo)0.05*Rand1Distrib();
	vNoise[2]=(FloatLo)0.05*Rand1Distrib()+CBOffset;
	VecAdd(vCBRef,vCBRef,vNoise);
	Normalize(vCBRef,vCBRef);
        /* then true CiBeta dir - U*(vCBRef) where U=[uone utwo uthree] */
	VecScale(vCBHere,uone,vCBRef[0]);
	VecScale(vTmp,utwo,vCBRef[1]);
	VecAdd(vCBHere,vCBHere,vTmp);
	VecScale(vTmp,uthree,vCBRef[2]);
	VecAdd(vCBHere,vCBHere,vTmp);
}

/* given three consecutive CAs, calculates and returns R+, R- and U basis vectors */
static void SetUpRefAxes(vec uone,vec utwo,vec uthree,vec vRMinus,vec vRPlus,vec vCALast, vec vCAHere, vec vCANext)
{
        VecSub(vRMinus,vCALast,vCAHere);
        VecSub(vRPlus,vCANext,vCAHere);
        /* to ensure vectors are same length when we later add them,
           normalize both vectors */
        Normalize(vRMinus,vRMinus);
        Normalize(vRPlus,vRPlus);
        Cross(uone,vRPlus,vRMinus);
        Normalize(uone,uone);
        VecAdd(utwo,vRPlus,vRMinus);
        Normalize(utwo,utwo);
        NegateVec(utwo,utwo);
        Cross(uthree,uone,utwo);
        Normalize(uthree,uthree);
}



/**
 *
 * LoadBryantPotential (Int4Ptr *ppiBryantPotential, FILE *fp,
 *                      Boolean bUsingPep, Boolean bNegateInput):
 *
 * Creates an array to store a Bryant-Lawrence potential in memory,
 * fills this array with information read in from "fp" and makes
 * "*ppiBryantPotential" point to this array. Returns 0 on success.
 *
 * The potential function is represented as a 3-D table of integers: each
 * entry in the table represents the potential between two residues,
 * multiplied by POTENTIAL_TABLE_SCALE_FACTOR to eliminate fractions.
 * For simplicity, the table is stored in memory as a 1-dimensional array.
 * The potential between residues "res1" and "res2" in bin number "bin" is
 * stored at position (bin * NUM_AA * NUM_AA + res1 * NUM_AA + res2) of
 * the array.
 *
 * Assumptions:
 *
 *   - The first line of the input file is either "PEP=TRUE" or
 *     "PEP=FALSE", depending on whether the data in the input file has
 *     row and column entries for PEPs.
 *   - There are NUM_BINS bins in the input file and each of these is a
 *     table with NUM_AA x NUM_AA entries if PEP=TRUE, or (NUM_AA - 1) x
 *     (NUM_AA - 1) entries if PEP=FALSE.
 *   - The amino acid row and column titles can be given in any order
 *     as long as they are all present.
 *
 * If "bNegateInput" is TRUE, each input value is multiplied by -1 before
 * being stored in the potential function array.
 *
 */
TrajErr LoadBryantPotentialEx (Int4Ptr *ppiBryantPotential, FILE *fp, Boolean bPotential,
                          Boolean bUsingPep, Boolean bNegateInput, Boolean bAlphaIndex)
{
    Int2 col [NUM_AA];
    Int2 binCounter, rowCounter, colCounter, bin, row, num_aa=0,cnt;
    double tempFloat;

    Char buffer [MAX_STRING_LENGTH + 1];
    CharPtr str = buffer;

    Int4Ptr potential;

    fscanf (fp, "%s\n", str);

    if (!(StringCmp (str, "PEP=TRUE"))) {
		num_aa = NUM_AA;
    } else if ((!(StringCmp (str, "PEP=FALSE"))) && (!bUsingPep)) {
		num_aa = NUM_AA - 1;
    } else {
		ErrPostEx (SEV_ERROR, 3, 1, "Potential file corrupt.");
		return ERR_FAIL;
    }

    potential = (Int4Ptr) MemNew (NUM_BINS * NUM_AA * NUM_AA * sizeof (Int4));

    for (binCounter=0; binCounter < NUM_BINS; binCounter++) {	
		fscanf (fp, "%s", str);

		if (!(StringCmp (str, "0-5"))) {
			bin = BIN_0_5;
		} else if (!(StringCmp (str, "5-6"))) {
			bin = BIN_5_6;
		} else if (!(StringCmp (str, "6-7"))) {
			bin = BIN_6_7;
		} else if (!(StringCmp (str, "7-8"))) {
			bin = BIN_7_8;
		} else if (!(StringCmp (str, "8-9"))) {
			bin = BIN_8_9;
		} else if (!(StringCmp (str, "9-10"))) {
			bin = BIN_9_10;
		} else {
			return ERR_FAIL;
		}

		for (colCounter=0; colCounter < num_aa; colCounter++) {
			fscanf (fp, "%s", str);
			col [colCounter] = GetIndexEx (str, bAlphaIndex);
			if (col[colCounter] == -1) {
				ErrPostEx (SEV_ERROR, 3, 2, "Potential file corrupt.");
				return ERR_FAIL;
			}

		}
		fscanf (fp, "\n");

		for (rowCounter=0; rowCounter < num_aa; rowCounter++) {
			fscanf (fp, "%s", str);
			row = GetIndexEx (str,bAlphaIndex);
			if (row == -1) {
				ErrPostEx (SEV_ERROR, 3, 3, "Potential file corrupt.");
				return ERR_FAIL;
			}

			for (colCounter=0; colCounter < num_aa; colCounter++) {
				fscanf (fp, "%lf", &tempFloat);
				if (bNegateInput) {
					tempFloat = tempFloat * -1.0;
				}

				potential [bin*NUM_AA*NUM_AA + row*NUM_AA + col[colCounter]]
					= (Int4) floor (tempFloat * POTENTIAL_TABLE_SCALE_FACTOR);
			}

			fscanf (fp, "\n");

		}

    }
	if(bPotential == FALSE) {
		*ppiBryantPotential = potential;
		return ERR_SUCCESS;
	}

    /* extra code inserted here to add hydrophobicity terms */
    /* read and ignore bin names */
    for (binCounter=0; binCounter < NUM_BINS; binCounter++)
		if (fscanf (fp, "%s", str)!=1) {
			/* missing values */
			ErrPostEx (SEV_ERROR, 3, 4, "Potential file missing hydrophobicity terms.");
			return ERR_FAIL;		
		}
		rowCounter = 0;
		for (rowCounter=0; rowCounter < num_aa; rowCounter++) {
			fscanf (fp, "%s", str);
			row = GetIndexEx(str,bAlphaIndex);
			if (row == -1) {
				ErrPostEx (SEV_ERROR, 3, 3, "Potential file corrupt.");
				return ERR_FAIL;
			}
			for (colCounter=0; colCounter < NUM_BINS; colCounter++) {
				fscanf (fp, "%lf", &tempFloat);
				if (bNegateInput) {
					tempFloat = tempFloat * -1.0;
				}
				/* tempFloat (mu rd) to approproate rows and colums */
				/* row contains residue index */
				for (cnt=0;cnt<NUM_AA;cnt++)
					potential [colCounter*NUM_AA*NUM_AA + row*NUM_AA + cnt]
						+= (Int4) floor (tempFloat * POTENTIAL_TABLE_SCALE_FACTOR);
				for (cnt=0;cnt<NUM_AA;cnt++)
					potential [colCounter*NUM_AA*NUM_AA + cnt*NUM_AA + row]
						+= (Int4) floor (tempFloat * POTENTIAL_TABLE_SCALE_FACTOR);
			}
			fscanf (fp, "\n");
		}
    if (num_aa == NUM_AA) {
	/* ensure pep-pep =0 */
		for (cnt=0;cnt<NUM_BINS;cnt++) {
			potential[cnt*NUM_AA*NUM_AA + (NUM_AA-1)*NUM_AA + (NUM_AA-1)]=0;
		}
    }
    /* extra code ends here */
    *ppiBryantPotential = potential;
    return ERR_SUCCESS;

}

/* Load the Bryant Potential from file using the hydrophobic-hydrophilic index */
TrajErr LoadBryantPotential (Int4Ptr PNTR ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput) 
{
	Boolean bPotential = TRUE;
	Boolean bAlphaIndex = FALSE;
	return LoadBryantPotentialEx(&(*ppiBryantPotential), fp, bPotential, bUsingPep, bNegateInput, bAlphaIndex);
}

/* Load the Bryant Potential from file using the alphabetic index */
TrajErr LoadBryantPotential_AlphaIndex (Int4Ptr PNTR ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput) 
{
	Boolean bPotential = TRUE;
	Boolean bAlphaIndex = TRUE;
	return LoadBryantPotentialEx(&(*ppiBryantPotential), fp, bPotential, bUsingPep, bNegateInput, bAlphaIndex);
}

/* Load in a 21x21x6 or 20x20x6 table of counts */
TrajErr LoadBryantTable (Int4Ptr PNTR ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput) 
{ 
	Boolean bPotential = FALSE;
	Boolean bAlphaIndex = FALSE;
	return LoadBryantPotentialEx(&(*ppiBryantPotential), fp, bPotential, bUsingPep, bNegateInput, bAlphaIndex);
}

/* Open a table of counts by filename */
TrajErr LoadBryantTableByName(Int4Ptr PNTR ppiBryantPotential, CharPtr filename,
                          Boolean bUsingPep, Boolean bNegateInput) 
{
	FILE *fp = NULL;
	Boolean bOpen = TRUE;
	if(GetBLFile(filename,&fp,bOpen)) return ERR_FAIL;
	return LoadBryantTable (&(*ppiBryantPotential), fp,bUsingPep, bNegateInput);
}


/**
 *
 * Load4DPotential (Int4Ptr *ppiBryantPotential, FILE *fp,
 *                      Boolean bUsingPep, Boolean bNegateInput):
 *
 * Creates an array to store a 4D potential in memory,
 * fills this array with information read in from "fp" and makes
 * "*ppiBryantPotential" point to this array. Returns 0 on success.
 *
 * The potential function is represented as a 4-D table of integers: each
 * entry in the table represents the potential between two residues,
 * multiplied by POTENTIAL_TABLE_SCALE_FACTOR to eliminate fractions.
 * For simplicity, the table is stored in memory as a 1-dimensional array.
 * The potential between residues "res1" and "res2" in 3D bin number "bin" is
 * stored at position  of
 * (bin * NUM_BINS * NUM_AA * NUM_AA + bin * NUM_AA * NUM_AA + res2 * NUM_AA + res1)
 * the array.
 *
 * Assumptions:
 *
 *   - The first line of the input file is either "PEP=TRUE" or
 *     "PEP=FALSE", depending on whether the data in the input file has
 *     row and column entries for PEPs.
 *   - There are (NUM_BINS * NUM_BINS4D) bins in the input file and each of these is a
 *     table with NUM_AA x NUM_AA entries if PEP=TRUE, or (NUM_AA - 1) x
 *     (NUM_AA - 1) entries if PEP=FALSE.
 *   - The amino acid row and column titles can be given in any order
 *     as long as they are all present.
 *
 * If "bNegateInput" is TRUE, each input value is multiplied by -1 before
 * being stored in the potential function array.
 *
 */
TrajErr Load4DPotentialEx (Int4Ptr *ppiBryantPotential, FILE *fp, Boolean bPotential,
                          Boolean bUsingPep, Boolean bNegateInput, Boolean bAlphaIndex)
{
    Int2 col [NUM_AA];
    Int2 binCounter, bin4Dcounter, rowCounter, colCounter, bin, ssbin, row, num_aa=0;
    Int4 index;
    double tempFloat;

    Char buffer [MAX_STRING_LENGTH + 1];
    CharPtr str = buffer;
	Char buffer2[MAX_STRING_LENGTH + 1];
	CharPtr str2 = buffer2;
	Int4 piSSDim[] = {21,21,6,6};
    Int4Ptr potential;

    fscanf (fp, "%s\n", str);

    if (!(StringCmp (str, "PEP=TRUE"))) {
		num_aa = NUM_AA;
    } else if ((!(StringCmp (str, "PEP=FALSE"))) && (!bUsingPep)) {
		num_aa = NUM_AA - 1;
    } else {
		ErrPostEx (SEV_ERROR, 3, 1, "Potential file corrupt.");
		return ERR_FAIL;
    }

    potential = (Int4Ptr) MemNew (NUM_BINS * NUM_BINS4D * NUM_AA * NUM_AA * sizeof (Int4));
	for (bin4Dcounter=0; bin4Dcounter < NUM_BINS4D; bin4Dcounter++) {
		for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
			fscanf (fp, "%s %s", str2, str);

			if (!(StringCmp (str, "0-5"))) {
				bin = BIN_0_5;
			} else if (!(StringCmp (str, "5-6"))) {
				bin = BIN_5_6;
			} else if (!(StringCmp (str, "6-7"))) {
				bin = BIN_6_7;
			} else if (!(StringCmp (str, "7-8"))) {
				bin = BIN_7_8;
			} else if (!(StringCmp (str, "8-9"))) {
				bin = BIN_8_9;
			} else if (!(StringCmp (str, "9-10"))) {
				bin = BIN_9_10;
			} else {
				return ERR_FAIL;
			}


			if (!(StringCmp (str2, "a-a"))) {
				ssbin = HELIX_HELIX;
			} else if (!(StringCmp (str2, "a-b"))) {
				ssbin = HELIX_STRAND;
			} else if (!(StringCmp (str2, "a-c"))) {
				ssbin = HELIX_COIL;
			} else if (!(StringCmp (str2, "b-b"))) {
				ssbin = STRAND_STRAND;
			} else if (!(StringCmp (str2, "b-c"))) {
				ssbin = STRAND_COIL;
			} else if (!(StringCmp (str2, "c-c"))) {
				ssbin = COIL_COIL;
			} else {
				return ERR_FAIL;
			}

			for (colCounter=0; colCounter < num_aa; colCounter++) {
				fscanf (fp, "%s", str);
				col [colCounter] = GetIndexEx (str, bAlphaIndex);
				if (col[colCounter] == -1) {
					ErrPostEx (SEV_ERROR, 3, 2, "Potential file corrupt.");
					return ERR_FAIL;
				}

			}
			fscanf (fp, "\n");

			for (rowCounter=0; rowCounter < num_aa; rowCounter++) {
				fscanf (fp, "%s", str);
				row = GetIndexEx (str,bAlphaIndex);
				if (row == -1) {
					ErrPostEx (SEV_ERROR, 3, 3, "Potential file corrupt.");
					return ERR_FAIL;
				}

				for (colCounter=0; colCounter < num_aa; colCounter++) {
					fscanf (fp, "%lf", &tempFloat);
					if (bNegateInput) {
						tempFloat = tempFloat * -1.0;
					}

     index = Get4DArrayValue(col[colCounter], row, bin, ssbin, piSSDim);
					potential [index] = (Int4) floor (tempFloat * POTENTIAL_TABLE_SCALE_FACTOR);
				}

				fscanf (fp, "\n");

			}
		}
	}
	if(bPotential == FALSE) {
		*ppiBryantPotential = potential;
		return ERR_SUCCESS;
	}

    *ppiBryantPotential = potential;
    return ERR_SUCCESS;

}


TrajErr Load4DTable (Int4Ptr PNTR ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput)
{
	Boolean bPotential = FALSE;
	Boolean bAlphaIndex = FALSE;
	return Load4DPotentialEx(&(*ppiBryantPotential), fp, bPotential, bUsingPep, bNegateInput, bAlphaIndex);
}


/* Open a table of counts by filename */
TrajErr Load4DTableByName(Int4Ptr PNTR ppiBryantPotential, CharPtr filename,
                          Boolean bUsingPep, Boolean bNegateInput)
{
	FILE *fp = NULL;
	Boolean bOpen = TRUE;
	if(GetBLFile(filename,&fp,bOpen)) return ERR_FAIL;
	return Load4DTable (&(*ppiBryantPotential), fp,bUsingPep, bNegateInput);
}



/**
 *
 * GetBryantPotential (Int4Ptr piBryantPotential, Int2 res1, Int2 res2,
 *                     Int2 bin):
 *
 * Returns the potential between residues "res1" and "res2" in bin number
 * "bin" as given by "piBryantPotential".
 *
 */

FloatHi GetBryantPotential (Int4Ptr piBryantPotential, Int2 res1, Int2 res2,
                            Int2 bin) {

    return piBryantPotential [bin * NUM_AA * NUM_AA + res1 * NUM_AA + res2]
	/ POTENTIAL_TABLE_SCALE_FACTOR;

}



/**
 *
 * FindCAlpha2 (PVNMA pvnmaHead):
 * Same as FindCAlpha, but returns NULL if no C-Alpha was found.
 *
 */

/*
PMAD FindCAlpha2 (PVNMA pvnmaHead) {

    PVNMA pvnmaHere = pvnmaHead;

    while (pvnmaHere !=  NULL) {

	if (IsAtomCAlpha ((PFB)(pvnmaHere->data.ptrvalue))) {
	    return (PMAD) (pvnmaHere->data.ptrvalue);
	}

	pvnmaHere = pvnmaHere->next;
    }

    return NULL;

}
*/


/**
 *
 * FindCBeta2 (PVNMA pvnmaHead):
 * Same as FindCBeta, but returns NULL if no C-Beta was found.
 *
 */

/*
PMAD FindCBeta2 (PVNMA pvnmaHead) {

    PVNMA pvnmaHere = pvnmaHead;

    while (pvnmaHere !=  NULL) {

	if (IsAtomCBeta ((PFB)(pvnmaHere->data.ptrvalue))) {
	    return (PMAD) (pvnmaHere->data.ptrvalue);
	}

	pvnmaHere = pvnmaHere->next;
    }

    return NULL;

}
*/


/**
 *
 * FindCBackbone (PVNMA pvnmaHead) {
 *
 * Returns a PMAD corresponding to the first C-backbone atom found that is
 * not a C-Alpha, or NULL if no such atom exists.
 *
 */

PMAD FindCBackbone (PVNMA pvnmaHead) {

    PVNMA pvnmaHere = pvnmaHead;

    while (pvnmaHere != NULL) {

	if ((IsAtomBackBone ((PFB)(pvnmaHere->data.ptrvalue)))
	    && (AtomicNumber ((PMAD)(pvnmaHere->data.ptrvalue)) == ATOMIC_NO_C)
	    && !(IsAtomCAlpha ((PFB)(pvnmaHere->data.ptrvalue)))) {
	    return (PMAD) (pvnmaHere->data.ptrvalue);
	}

	pvnmaHere = pvnmaHere->next;
    }

    return NULL;

}



/**
 *
 * FindNBackbone (PVNMA pvnmaHead) {
 *
 * Returns a PMAD corresponding to the first backbone Nitrogen found,
 * or NULL if no such atom exists.
 *
 */

PMAD FindNBackbone (PVNMA pvnmaHead) {

    PVNMA pvnmaHere = pvnmaHead;

    while (pvnmaHere != NULL) {

	if ((IsAtomBackBone ((PFB)(pvnmaHere->data.ptrvalue)))
	    && (AtomicNumber ((PMAD)(pvnmaHere->data.ptrvalue)) == ATOMIC_NO_N)) {
	    return (PMAD) (pvnmaHere->data.ptrvalue);
	}

	pvnmaHere = pvnmaHere->next;
    }

    return NULL;

}



/**
 *
 * GetBinNumber (PALD pald1, PALD pald2):
 *
 * Takes two PALDs, computes the distance between them and returns the
 * appropriate MK_BIN value, or -1 if the distance doesn't fall into a
 * valid bin.
 *
 */

Int2 GetBinNumber (PALD pald1, PALD pald2) {

    FloatLo distance;
    vec vec1, vec2;

    vec1[0] = pald1->pflvData[0];
    vec1[1] = pald1->pflvData[1];
    vec1[2] = pald1->pflvData[2];

    vec2[0] = pald2->pflvData[0];
    vec2[1] = pald2->pflvData[1];
    vec2[2] = pald2->pflvData[2];

    VecSub (vec1, vec1, vec2);
    distance = getMag (vec1);

    if ((distance > 0.0) && (distance <= 5.0)) {
	return BIN_0_5;
    } else if ((distance > 5.0) && (distance <= 6.0)) {
	return BIN_5_6;
    } else if ((distance > 6.0) && (distance <= 7.0)) {
	return BIN_6_7;
    } else if ((distance > 7.0) && (distance <= 8.0)) {
	return BIN_7_8;
    } else if ((distance > 8.0) && (distance <= 9.0)) {
	return BIN_8_9;
    } else if ((distance > 9.0) && (distance <= 10.0)) {
	return BIN_9_10;
    } else {
/*		ErrPostEx(SEV_FATAL,1,1,"Bad Distance passed to GetBinNumber: %f",distance);*/
	return -1;
    }

}

/*
*
*GetSSBin(PMGD pmgd1, PMGD pmgd2)
*
*Looks up VAST secondary structure assignments in molecule graphs
*returns integer specifying the combination of ss structures elements
*returns -1 if combination of elements cannot be found !!!
*
*/


Int2 GetSSBin(PMGD pmgd1, PMGD pmgd2)
{
   if((pmgd1 == NULL) || (pmgd2 == NULL))
      return -1;
   /*printf("%d %d\n", pmgd1->bNCBISecStru, pmgd2->bNCBISecStru);
   fflush(stdout);*/

   if ((pmgd1->bNCBISecStru & (Byte) SS_HELIX) &&
      (pmgd2->bNCBISecStru & (Byte) SS_HELIX))
      return HELIX_HELIX;
   else if ((pmgd1->bNCBISecStru & (Byte) SS_HELIX) &&
           ((pmgd2->bNCBISecStru & (Byte) SS_STRAND) || (pmgd2->bNCBISecStru & (Byte) SS_SHEET)))
      return HELIX_STRAND;
   else if (((pmgd1->bNCBISecStru & (Byte) SS_STRAND) || (pmgd1->bNCBISecStru & (Byte) SS_SHEET)) &&
           (pmgd2->bNCBISecStru & (Byte) SS_HELIX))
      return HELIX_STRAND;
   else if (((pmgd1->bNCBISecStru & (Byte) SS_STRAND) || (pmgd1->bNCBISecStru & (Byte) SS_SHEET)) &&
           ((pmgd2->bNCBISecStru & (Byte) SS_STRAND) || (pmgd2->bNCBISecStru & (Byte) SS_SHEET)))
      return STRAND_STRAND;
   else if ((pmgd1->bNCBISecStru & (Byte) SS_HELIX) &&
           ((pmgd2->bNCBISecStru & (Byte) SS_TURN) || (pmgd2->bNCBISecStru == 0)))
      return HELIX_COIL;
   else if (((pmgd1->bNCBISecStru & (Byte) SS_TURN) || (pmgd1->bNCBISecStru == 0)) &&
           (pmgd2->bNCBISecStru & (Byte) SS_HELIX))
      return HELIX_COIL;
   else if (((pmgd1->bNCBISecStru & (Byte) SS_STRAND) || (pmgd1->bNCBISecStru & (Byte) SS_SHEET)) &&
           ((pmgd2->bNCBISecStru & (Byte) SS_TURN) || (pmgd2->bNCBISecStru == 0)))
      return STRAND_COIL;
   else if (((pmgd1->bNCBISecStru & (Byte) SS_TURN) || (pmgd1->bNCBISecStru == 0)) &&
           ((pmgd2->bNCBISecStru & (Byte) SS_STRAND) || (pmgd2->bNCBISecStru & (Byte) SS_SHEET)))
      return STRAND_COIL;
   else if (((pmgd1->bNCBISecStru & (Byte) SS_TURN) || (pmgd1->bNCBISecStru == 0)) &&
           ((pmgd2->bNCBISecStru & (Byte) SS_TURN) || (pmgd2->bNCBISecStru == 0)))
      return COIL_COIL;
   else
      return -1;
}

/*
* GetSBBin(Int2 bin, Int2 ssbin)
* takes Bryant-Lawrence potential and secondary structure potential bins and
* combines them into Bryant-Secondary potential bin number
* see defines in potential.h
*/


Int2 GetSBBin(Int2 bin, Int2 ssbin)
{
	if ((ssbin == HELIX_COIL) || (ssbin == STRAND_COIL) || (ssbin == COIL_COIL))
		return COIL_ALL;
	if (ssbin == HELIX_STRAND)
		return HELIX_STRAND_ALL;
	if (ssbin == HELIX_HELIX) {
		if ((bin == BIN_0_5) || (bin == BIN_5_6) || (bin == BIN_6_7))
			return HELIX_NEAR;
		if ((bin == BIN_7_8) || (bin == BIN_8_9) || (bin == BIN_9_10))
			return HELIX_FAR;
	}
	if (ssbin == STRAND_STRAND) {
		if ((bin == BIN_0_5) || (bin == BIN_5_6) || (bin == BIN_6_7))
			return STRAND_NEAR;
		if ((bin == BIN_7_8) || (bin == BIN_8_9) || (bin == BIN_9_10))
			return STRAND_FAR;
	}
        return -1;
}

/*
* GetSBBin2(Int2 bin, Int2 ssbin)
* takes Bryant-Lawrence potential and secondary structure potential bins and
* combines them into another Bryant-Secondary potential bin number
* see defines in potential.h
*/


Int2 GetSBBin2(Int2 bin, Int2 ssbin)
{
	if ((ssbin == HELIX_COIL) || (ssbin == STRAND_COIL))
		return HELIX_STRAND_COIL;
	if (ssbin == HELIX_STRAND)
		return HELIX_STRAND;
	if (ssbin == HELIX_HELIX)
		return HELIX_HELIX;
	if (ssbin == STRAND_STRAND) {
		if ((bin == BIN_0_5) || (bin == BIN_5_6) || (bin == BIN_6_7))
			return STRAND_NEAR;
		if ((bin == BIN_7_8) || (bin == BIN_8_9) || (bin == BIN_9_10))
			return STRAND_FAR;
	}
	if (ssbin == COIL_COIL)
		return COIL_COIL;
	return -1;
}


/**
 *
 * CreatePsBList (DValNodePtr *ppdnmgResult, DValNodePtr pdnListPmmds,
 *                Int2 iModelNum, Boolean bUsingPep, Boolean usecaonly):
 *
 * Takes a list of molecules "pdnListPmmds" and for each residue in each
 * molecule creates a pseudo-beta carbon MAD/MGD pair and adds it to the
 * list of MGDs pointed to by "*ppdnmgResult". The choice field of each
 * DValNode, the iDomain field of each MGD and the iIndex field of each
 * MAD store the enumerated residue number. The virtual PEP unit belonging
 * to residue i is enumerated as residue -i.
 *
 * Returns 0 on success.
 *
 */

Int2 CreatePsBList (DValNodePtr *ppdnmgResult, DValNodePtr pdnListPmmds,
		    Int2 iModelNum, Boolean bUsingPep, Boolean bModelling, Boolean usecaonly) {

    Boolean bSkipThis = FALSE;

    PDNMG pdnmg = NULL;
    PMMD pmmd = NULL;
    PMAD pmadA = NULL, pmadB = NULL, pmadN = NULL, pmadC = NULL, pmadAnext = NULL, pmadAlast = NULL, pmadPsB = NULL;
    PALD paldA = NULL, paldB = NULL, paldN = NULL, paldC = NULL, paldAnext = NULL, paldPsB = NULL;
    PMGD pmgdHere = NULL, pmgdNext = NULL, pmgdLast = NULL, pmgdPsB = NULL;

    vec vecA, vecB, vecN, vecC, vecAnext, bis, crs, vecPsB;
    vec vCAHere, vCALast, vCANext;
    vec uone, utwo, uthree, vRMinus, vRPlus, vCBHereRef, vCBHere;
    Int2 resHere = 0, resNext = 0;
    Int2 count = 0;
	
	if(gCBTable_IsLoaded == FALSE) { 
		if(LoadCBTable() == ERR_FAIL) return 1;
		gCBTable_IsLoaded = TRUE;
	}

    /* For each molecule in the list */
    while (pdnListPmmds != NULL) {
		pmmd = (PMMD)(pdnListPmmds->data.ptrvalue);
		pdnmg = pmmd->pdnmgHead;
		count = 1;

		/* For each residue in the molecule */
		while (pdnmg != NULL) {
			bSkipThis = FALSE;

			pmgdHere = (PMGD)(pdnmg->data.ptrvalue);
			
			if(bModelling && ( (pmgdHere->bReserved == 'X') || (pmgdHere->pcIUPAC[0] == 'X') || !(pmgdHere->bWhat&DICT_GLOBAL && pmgdHere->bWhat&RES_AA))) {
				bSkipThis = TRUE;
			}

			if ((pmadA = FindCAlpha (pmgdHere->pvnmaAHead)) == NULL) {
				ErrPostEx (SEV_WARNING, 4, 1, "No alpha carbon was found.");
				bSkipThis = TRUE;
			}

			if ((!bSkipThis) && ((paldA = GetAtomLocs (pmadA, iModelNum)) == NULL)) {
				ErrPostEx (SEV_WARNING, 4, 2, "Alpha carbon has no coordinates.");
				bSkipThis = TRUE;
			}

			if (!bSkipThis) {
				if (pdnmg->next != NULL) {
					pmgdNext = (PMGD)(pdnmg->next->data.ptrvalue);
					pmadAnext = FindCAlpha (pmgdNext->pvnmaAHead);
				} else {
					pmadAnext = NULL;
				}

				if (pdnmg->last != NULL) {
					pmgdLast = (PMGD)(pdnmg->last->data.ptrvalue);
					pmadAlast = FindCAlpha (pmgdLast->pvnmaAHead);
				} else {
					pmadAlast = NULL;
				}

				/* Case 1: The residue has a beta carbon */
				if (((pmadB = FindCBeta (pmgdHere->pvnmaAHead)) != NULL) &&
					((paldB = GetAtomLocs (pmadB, iModelNum)) != NULL) && !usecaonly) {

 					vecA[0] = paldA->pflvData[0];
					vecA[1] = paldA->pflvData[1];
					vecA[2] = paldA->pflvData[2];

					vecB[0] = paldB->pflvData[0];
					vecB[1] = paldB->pflvData[1];
					vecB[2] = paldB->pflvData[2];

					VecSub (vecPsB, vecB, vecA);
					Normalize (vecPsB, vecPsB);
					VecScale (vecPsB, vecPsB, (FloatLo) PSB_SCALE_FACTOR);
					VecAdd (vecPsB, vecA, vecPsB);

				/* Case 2: The residue has non-C-Alpha backbone atoms */
				} else if (((pmadN = FindNBackbone (pmgdHere->pvnmaAHead)) != NULL) &&
					((pmadC = FindCBackbone (pmgdHere->pvnmaAHead)) != NULL) &&
					((paldN = GetAtomLocs (pmadN, iModelNum)) != NULL) &&
					((paldC = GetAtomLocs (pmadC, iModelNum)) != NULL) && !usecaonly) {

					vecA[0] = paldA->pflvData[0];
					vecA[1] = paldA->pflvData[1];
					vecA[2] = paldA->pflvData[2];

					vecN[0] = paldN->pflvData[0];
					vecN[1] = paldN->pflvData[1];
					vecN[2] = paldN->pflvData[2];

					vecC[0] = paldC->pflvData[0];
					vecC[1] = paldC->pflvData[1];
					vecC[2] = paldC->pflvData[2];

					VecSub (vecN, vecN, vecA);
					VecSub (vecC, vecC, vecA);

					VecAdd (bis, vecN, vecC);
					Normalize (bis, bis);
					Cross (crs, vecN, vecC);
					Normalize (crs, crs);

					VecScale (crs, crs, (FloatLo) (PSB_SCALE_FACTOR * 0.816497));
					VecScale (bis, bis, (FloatLo) (PSB_SCALE_FACTOR * 0.57735));
					VecSub (vecPsB, crs, bis);
					VecAdd (vecPsB, vecA, vecPsB);

				/* Case 3: C-Alpha backbone only */
				} else if ((pmadAnext != NULL) && (pmadAlast != NULL) &&
					((paldAnext = GetAtomLocs (pmadAnext, iModelNum)) != NULL) &&
					(GetAtomLocs (pmadAlast, iModelNum) != NULL)) {

					resHere=GetResIdxFromMGD(NCBIstdaaUC,pmgdHere);
					resNext=GetResIdxFromMGD(NCBIstdaaUC,pmgdNext);

					GetCoOrds(pmgdLast," CA ",vZero,vCALast,iModelNum);
					GetCoOrds(pmgdHere," CA ",vZero,vCAHere,iModelNum);
					GetCoOrds(pmgdNext," CA ",vZero,vCANext,iModelNum);
					/* construct local axis system */
	
					SetUpRefAxes(uone,utwo,uthree,vRMinus,vRPlus,vCALast,vCAHere,vCANext);
					/* get CB dir unit vector in vCBTrue, relative to CA */
					/* Ref indicates with respect to the reference U axes */

					FindCBDir(uone,utwo,uthree,vCALast,vCANext,resHere,resNext,vCBHereRef,vCBHere,0);
					/* and assign it for use with Place Rotamer later */
					VecScale(vecPsB,vCBHere,bl_cacb[resHere-1]);
					VecAdd(vecPsB,vecPsB,vCAHere);

				} else {
/* 					ErrPostEx (SEV_WARNING, 4, 3, "Insufficient data to compute pseudo-beta carbon.");*/
					bSkipThis = TRUE;
				}
			}

			if (!bSkipThis) {
				pmgdPsB = NewMGD ();

				pmgdPsB->pfbParent = pmgdHere->pfbParent;
				pmgdPsB->bWhat = '\0';
				pmgdPsB->bReserved = pmgdHere->bReserved;
				pmgdPsB->iDomain = count;
				pmgdPsB->pcGraphName = StringSave (pmgdHere->pcGraphName);
				pmgdPsB->pcGraphNum = StringSave (pmgdHere->pcGraphNum);
				pmgdPsB->pcIUPAC = StringSave (pmgdHere->pcIUPAC);
				pmgdPsB->pvnContainedBy = NULL;
				pmgdPsB->pvnmbBHead = NULL;
				pmgdPsB->ppflBoundBox = NULL;
				pmgdPsB->iAtomCount = 1;

				pmadPsB = NewMAD ();

				pmadPsB->pfbParent = pmadA->pfbParent;
				pmadPsB->bWhat = '\0';
				pmadPsB->pcAName = StringSave("VCB ");
				pmadPsB->pvnBonds = NULL;
				pmadPsB->pvnContainedBy = NULL;
				pmadPsB->iIndex = count;

				paldPsB = NewALD ();

				paldPsB->pfbParent = (PFB) pmadPsB;
				paldPsB->bWhat = '\0';
				paldPsB->next = NULL;
				paldPsB->iCoordSet = (Char) iModelNum;
				paldPsB->iFloatNo = 4;

				paldPsB->pflvData = FLVector (0, (Int4) paldPsB->iFloatNo);

				paldPsB->pflvData[0] = vecPsB[0];
				paldPsB->pflvData[1] = vecPsB[1];
				paldPsB->pflvData[2] = vecPsB[2];

				ValNodeAddPointer (&(pmadPsB->pvnalLocate), iModelNum, paldPsB);
				ValNodeAddPointer (&(pmgdPsB->pvnmaAHead), 0, pmadPsB);
				DValNodeAddPointer (ppdnmgResult, count, pmgdPsB);

				if (bUsingPep && (pmadAnext != NULL) &&
					((paldAnext = GetAtomLocs (pmadAnext, iModelNum)) != NULL)) {

					vecA[0] = paldA->pflvData[0];
					vecA[1] = paldA->pflvData[1];
					vecA[2] = paldA->pflvData[2];

					vecAnext[0] = paldAnext->pflvData[0];
					vecAnext[1] = paldAnext->pflvData[1];
					vecAnext[2] = paldAnext->pflvData[2];

					VecAdd (vecPsB, vecA, vecAnext);
					VecScale (vecPsB, vecPsB, 0.5);

					pmgdPsB = NewMGD ();

					pmgdPsB->pfbParent = pmgdHere->pfbParent;
					pmgdPsB->bWhat = '\0';
					pmgdPsB->bReserved = pmgdHere->bReserved;
					pmgdPsB->iDomain = (Int2) -count;
					pmgdPsB->pcGraphName = StringSave ("Pep");
					pmgdPsB->pcGraphNum = StringSave (pmgdHere->pcGraphNum);
					pmgdPsB->pcIUPAC = StringSave ("-");
					pmgdPsB->pvnContainedBy = NULL;
					pmgdPsB->pvnmbBHead = NULL;
					pmgdPsB->ppflBoundBox = NULL;
					pmgdPsB->iAtomCount = 1;

					pmadPsB = NewMAD ();

					pmadPsB->pfbParent = pmadA->pfbParent;
					pmadPsB->bWhat = '\0';
					pmadPsB->pcAName = StringSave("VCB ");
					pmadPsB->pvnBonds = NULL;
					pmadPsB->pvnContainedBy = NULL;
					pmadPsB->iIndex = (Int2) -count;

					paldPsB = NewALD ();

					paldPsB->pfbParent = (PFB) pmadPsB;
					paldPsB->bWhat = '\0';
					paldPsB->next = NULL;
					paldPsB->iCoordSet = (Char) iModelNum;
					paldPsB->iFloatNo = 4;

					paldPsB->pflvData = FLVector (0, (Int4) paldPsB->iFloatNo);

					paldPsB->pflvData[0] = vecPsB[0];
					paldPsB->pflvData[1] = vecPsB[1];
					paldPsB->pflvData[2] = vecPsB[2];

					ValNodeAddPointer (&(pmadPsB->pvnalLocate), iModelNum, paldPsB);
					ValNodeAddPointer (&(pmgdPsB->pvnmaAHead), 0, pmadPsB);
					DValNodeAddPointer (ppdnmgResult, (Int2) (-count), pmgdPsB);
					DValNodeAddPointer (&pdnListPepMGDs, 0, pmgdPsB);
				}
			}
			pdnmg = pdnmg->next;
			count++;
		}
		pdnListPmmds = pdnListPmmds->next;
    }

    return 0;
}



void LIBCALLBACK FreeMAD2(Pointer ptr)
{    
	PMAD pmadThis;

	pmadThis=(PMAD)ptr;
	if (pmadThis) {
		if (pmadThis->pvnContainedBy) pmadThis->pvnContainedBy = ValNodeFree(pmadThis->pvnContainedBy);
		if (pmadThis->pvnBonds) pmadThis->pvnBonds = ValNodeFree(pmadThis->pvnBonds);
		if (pmadThis->pvnalLocate) FreeListVNAL(pmadThis->pvnalLocate);
		if (pmadThis->pcAName) pmadThis->pcAName = MemFree(pmadThis->pcAName);
		pmadThis = MemFree(pmadThis);
	}
}

/**
 *
 * FreeMGD2 (PMGD pmgdThis):
 *
 * Same as FreeMGD except that MGDs with negative ResId's are not
 * freed. For use with PsB lists and AdjLists.
 *
 */

void LIBCALLBACK FreeMGD2 (Pointer ptr) {
	PMGD pmgdThis = (PMGD)ptr;
    if (pmgdThis && (pmgdThis->iDomain >= 0)) {
        if (pmgdThis->pvnContainedBy) ValNodeFree (pmgdThis->pvnContainedBy);
        if (pmgdThis->pcGraphName) pmgdThis->pcGraphName = MemFree (pmgdThis->pcGraphName);
        if (pmgdThis->pvnmaAHead)
	    FreeVNDataFn(pmgdThis->pvnmaAHead,(pFreeFunc)FreeMAD2); 
        if (pmgdThis->pvnmbBHead) FreeListVNMB (pmgdThis->pvnmbBHead);
        if (pmgdThis->ppflBoundBox) FLMatrixFree (pmgdThis->ppflBoundBox, 0, 0);
        if (pmgdThis->pcIUPAC) pmgdThis->pcIUPAC = MemFree (pmgdThis->pcIUPAC);
        if (pmgdThis->pcGraphNum) pmgdThis->pcGraphNum = MemFree (pmgdThis->pcGraphNum);
		pmgdThis = MemFree(pmgdThis);
    }
}

void LIBCALLBACK FreeMGD3 (Pointer ptr) {
	PMGD pmgdThis;

	pmgdThis=(PMGD)ptr;
    if (pmgdThis) {
        if (pmgdThis->pvnContainedBy)
            ValNodeFree (pmgdThis->pvnContainedBy);
        if (pmgdThis->pcGraphName)
            MemFree (pmgdThis->pcGraphName);
        if (pmgdThis->pvnmaAHead)
	    FreeVNDataFn(pmgdThis->pvnmaAHead,(pFreeFunc)FreeMAD2); 
        if (pmgdThis->pvnmbBHead)
            FreeListVNMB (pmgdThis->pvnmbBHead);
        if (pmgdThis->ppflBoundBox)
            FLMatrixFree (pmgdThis->ppflBoundBox, 0, 0);
        if (pmgdThis->pcIUPAC)
            MemFree (pmgdThis->pcIUPAC);
        if (pmgdThis->pcGraphNum)
            MemFree (pmgdThis->pcGraphNum);
	MemFree(pmgdThis);
    }
}

/**
 *
 * FreePsBList (DValNodePtr *ppdnmgHead):
 *
 * Frees the memory allocated to *ppdnmgHead EXCEPT for the MGDs belonging
 * to pseudo-beta carbons that are PEPs. These are freed by FreeAdjList.
 * Sets *ppdnmgHead to NULL.
 *
 */

void FreePsBList (DValNodePtr *ppdnmgHead) {

    if (*ppdnmgHead != NULL) {
	DValNodeFreeData (*ppdnmgHead, (pFreeFunc) FreeMGD2);
    }

    *ppdnmgHead = NULL;

}

/* Allocates Int4 array for given POT_TYPE */
void NewContactCount(POT_TYPE pot_type, Int4Ptr *ppiCount, Int4 *piCount)
{
	if(ppiCount == NULL) {
		ErrPostEx(SEV_FATAL,0,0,"NewContactCount: You must pass a reference to a valid pointer.");
		return;
	}
	if(pot_type == POT_BL) {
		*piCount = NUM_BINS * NUM_AA * NUM_AA;
		*ppiCount = (Int4Ptr) MemNew ((size_t) (*piCount * sizeof(Int4)));
	} else if ((pot_type == POT_SS)||(pot_type == POT_SB) || (pot_type == POT_SB2)) {
		*piCount = NUM_BINS_SS * NUM_AA * NUM_AA;
		*ppiCount = (Int4Ptr) MemNew ((size_t) (*piCount * sizeof(Int4)));
	} else if (pot_type == POT_4D) {
		*piCount = NUM_BINS * NUM_BINS4D * NUM_AA * NUM_AA;
		*ppiCount = (Int4Ptr) MemNew ((size_t) (*piCount * sizeof(Int4)));
	}
	return;
}


/***********************************************************************************
* NewBLPotential                                                                       *
*                                                                                  *
* Creates an array to store a Bryant-Lawrence potential in memory.                 *
* The potential function is represented as a 3-D table of integers:                *
* For simplicity, the table is stored in memory as a 1-dimensional array.          *
* The potential between residues "res1" and "res2" in bin number "bin" is          *
* stored at position (bin * NUM_AA * NUM_AA + res1 * NUM_AA + res2) of             *
* the array.  NUM_AA typically represents 20 amino acids plus the peptide backbone.*
* If not counting the peptide backbone, set bUsingPep to FALSE.                    *
* There are NUM_BINS bins                                                          *
***********************************************************************************/
Int4 NewBLPotential (Int4Ptr *ppiBryantPotential)
{
	Int4 iCount = 0;
	NewContactCount(POT_BL,&(*ppiBryantPotential),&iCount);
	return 0;
}

/***********************************************************************************
* New4DPotential                                                                     *
*                                                                                  *
* Creates an array to store a potential in memory.                                 *
* The potential function is represented as a 4-D table of integers:                *
* For simplicity, the table is stored in memory as a 1-dimensional array.          *
* the conversion is done as follows:                                               *
* Get4DArrayValue(A,B,C,D,Z) (Int4) ((D)*(Z[2])*(Z[1])*(Z[0])+(C)*(Z[1])*(Z[0])+(B)*(Z[0])+(A))
* where A and B are NUM_AA, C is the thirs dimension, D is fourth dimension,       *
* and Z aray of dimensions of the potential i.e.NUM_AA, NUM_AA, NUM_BINS, NUM_BINS4D *
***********************************************************************************/
Int4 New4DPotential (Int4Ptr *ppiBryantPotential)
{
	Int4 iCount = 0;
	NewContactCount(POT_4D,&(*ppiBryantPotential),&iCount);
	return 0;
}


/**
 *
 * FreeBryantPotential (Int4Ptr *ppiBryantPotential):
 *
 * Frees memory allocated to the potential function array and makes
 * "*ppiBryantPotential" point to NULL.
 * 
 */
void FreeContactCount(Int4Ptr *ppiCount) 
{
    if (*ppiCount != NULL) {
		MemFree (*ppiCount);
    }
    *ppiCount = NULL;
}
void FreeBryantPotential (Int4Ptr *ppiBryantPotential) {FreeContactCount(&(*ppiBryantPotential));}




/**************************************************************************************
* CountStructureContacts                                                              *
*                                                                                     *
* This function counts the contacts between protein residues for a variety of methods * 
* Requires a list of models with the model number specified.  You must specify 
* whether there is an exclusive window of specified size and whether to count peptide 
* contacts.  The 
* parameter indicates whether to count the contacts between residues identified in *
* the bReserved field of the pmgd -> which can be used in a Threaded model.  This  *
* function will also ignore residues in the bReserved field that are marked 'X'    *
***********************************************************************************/
Int4 CountStructureContacts (Int4Ptr *ppiCount, DValNodePtr pdnListPmmds, Int2 iModelNum, 
		Boolean bInclusiveWindow, Int2 iWindowSize, Boolean bUsingPep, Boolean bModelling, Boolean bModel, POT_TYPE pot_type)
{
	Boolean bExclusiveWindow = (Boolean) !bInclusiveWindow;

	DValNodePtr pdnmgHeadPsBList = NULL;
	DValNodePtr pdnmgHere = NULL;
	ValNodePtr pvnAtomList = NULL, pvnHere = NULL;

	PWS pwsHead = NULL;
	PWS pwsThis = NULL;
	vec middle;
	Int2 iResId1, iResId2, iMolId1, iMolId2, res1, res2, bin, ssbin, sbbin;

	PMMD pmmd1, pmmd2;
	PMGD pmgd1, pmgd2;
	PMAD pmad1, pmad2;
	PALD pald1, pald2;

	Int4Ptr piCount = *ppiCount;
	Char myresname[2];
	Int4 icount = 0, itotalcounts = 0;
	Int4 piDim[4]; 
	Int4 value = 0;

	piDim[0] = NUM_AA;
	piDim[1] = NUM_AA;
	piDim[2] = NUM_BINS;
	piDim[3] = 0;
	if((pot_type == POT_SS)||(pot_type == POT_SB) || (pot_type == POT_SB2)) piDim[2] = NUM_BINS_SS;
	if(pot_type == POT_4D) piDim[3] = NUM_BINS4D;


	pdnListPepMGDs = NULL;
	if(CreatePsBList (&pdnmgHeadPsBList, pdnListPmmds, iModelNum, bUsingPep, bModelling, FALSE)) return 1;
	pdnmgHere = pdnmgHeadPsBList;

    /* Create a world containing all of the pseudo-beta carbons */
    while (pdnmgHere != NULL) {
		pmgd1 = (PMGD) pdnmgHere->data.ptrvalue;
			pwsThis = AddtoWorld (pwsHead, iModelNum, (PFB) pmgd1);
			if (pwsHead == NULL) {
				pwsHead = pwsThis;
			}
		pdnmgHere = pdnmgHere->next;
    }
	
	InstantiateWorld (1, pwsHead);
    pdnmgHere = pdnmgHeadPsBList;

    /* For each pseudo-beta carbon */
    while (pdnmgHere != NULL) {
		itotalcounts++;
	    pmgd1 = (PMGD)(pdnmgHere->data.ptrvalue);
		pmad1 = (PMAD)(pmgd1->pvnmaAHead->data.ptrvalue);
	    pald1 = GetAtomLocs (pmad1, iModelNum);

	    pmgd1 = (PMGD)(pmad1->pfbParent);
	    pmmd1 = (PMMD)(pmgd1->pfbParent);

		iResId1=pdnmgHere->choice;
	    iMolId1 = ((PDNMM)(pmmd1->pdnmmLink))->choice;

	    middle[0] = pald1->pflvData[0];
	    middle[1] = pald1->pflvData[1];
	    middle[2] = pald1->pflvData[2];

	    pvnAtomList = FindAtomsIn (pwsHead, middle, 10.0);
	    pvnHere = pvnAtomList;
		/* For each atom */
	    while (pvnHere != NULL) {
			itotalcounts++;
			pald2 = (PALD)(pvnHere->data.ptrvalue);
			pmad2 = (PMAD)(pald2->pfbParent);
			pmgd2 = (PMGD)(pmad2->pfbParent);
			pmmd2 = (PMMD)(pmgd2->pfbParent);

			iMolId2 = ((PDNMM)(pmmd2->pdnmmLink))->choice;
			iResId2 = pmad2->iIndex;

			if (((bExclusiveWindow) &&
				  ((iMolId2 != iMolId1) || ((iMolId2 == iMolId1) &&
				  (abs(abs(iResId1) - abs(iResId2)) > iWindowSize)))) ||
				((bInclusiveWindow) &&
				  (iMolId1 == iMolId2) && (iResId1 != iResId2) &&
				   (abs(abs(iResId1) - abs(iResId2)) >= iWindowSize))) {

				if (iResId1 > 0) {
					if(bModel == TRUE) {
						myresname[0] = (Char) toupper(pmgd1->bReserved);
					} else {
						myresname[0] = GetAAFromIDict(pmgd1);
					}
					myresname[1]='\0';
					res1 = GetAlphaIndex (myresname);
				} else {
					res1 = MK_PEP;
				}

				if (iResId2 > 0) {
					if(bModel == TRUE) {
						myresname[0] = (Char) toupper(pmgd2->bReserved);
					} else {
						myresname[0] = GetAAFromIDict(pmgd2);
					}
					myresname[1]='\0';
					res2 = GetAlphaIndex (myresname);
				} else {
					res2 = MK_PEP;
				}
				if(res1 < 0 || res2 < 0) {
					ErrPostEx(SEV_ERROR,0,0,"Invalid Residue.");
				} else {
					if((bin = GetBinNumber (pald1, pald2)) != -1) {
						if (pot_type == POT_BL) {
							value = Get3DArrayValue(res1,res2,bin,piDim);
							piCount [value]++;
							icount++;
						} else if (pot_type == POT_4D || pot_type == POT_SS) {
							if((ssbin = GetSSBin(pmgd1,pmgd2)) != -1) {
								if(pot_type == POT_4D) {
									value = Get4DArrayValue(res1,res2,bin,ssbin,piDim);
								} else if (pot_type == POT_SS) {
									value = Get3DArrayValue(res1,res2,ssbin,piDim);
								}
								piCount [value]++;
								icount++;
							}
						} else if (pot_type == POT_SB) {
							if((ssbin = GetSSBin(pmgd1,pmgd2)) != -1) {
								if((sbbin = GetSBBin(bin, ssbin)) != -1) {
									value = Get3DArrayValue(res1,res2,sbbin,piDim);
									piCount [value]++;
									icount++;
								}
							}
						} else if (pot_type == POT_SB2){
							if((ssbin = GetSSBin(pmgd1,pmgd2)) != -1) {
								if((sbbin = GetSBBin2(bin, ssbin)) != -1) {
									value = Get3DArrayValue(res1,res2,sbbin,piDim);
									piCount [value]++;
									icount++;
								}
							}
						}
					}
     				}
			}
			pvnHere = pvnHere->next;

	    }
	    ValNodeFree(pvnAtomList);
	    pdnmgHere = pdnmgHere->next;
    }
	ErrPostEx(SEV_INFO,0,0, "Possible Counts: %ld  Good Counts: %ld\n",(long) itotalcounts, (long)icount);
	*ppiCount = piCount;
    FreePsBList (&pdnmgHeadPsBList);
	pdnmgHeadPsBList = NULL;
	/* pdnListPepMGDs = DValNodeFreeData (pdnListPepMGDs, (pFreeFunc) FreeMGD3); */
	/* FreeAdjList(&pdnmgHeadPsBList); */
	FreeAllWorlds ();

	return 0;
}

/***********************************************************************************
* AddtoBLPotential                                                                 *
*                                                                                  *
* This function counts the contacts between residues within 6 different distance   *
* bins as specified by the Bryant-Lawrence potential.  Requires a list of models   *
* with the model number specified.  You must specify whether there is an exclusive *
* window of specified size and whether to count peptide contacts.  The last        *
* parameter indicates whether to count the contacts between residues identified in *
* the bReserved field of the pmgd -> which can be used in a Threaded model.  This  *
* function will also ignore residues in the bReserved field that are marked 'X'    *
***********************************************************************************/
Int4  AddtoBLPotential (Int4Ptr *ppiBryantPotential, DValNodePtr pdnListPmmds,
						Int2 iModelNum, Boolean bInclusiveWindow,
						Int2 iWindowSize, Boolean bUsingPep,
						Boolean bModelling, Boolean bModel)
{
	return AddtoBLPotentialEx(ppiBryantPotential, pdnListPmmds,iModelNum,
                            bInclusiveWindow,iWindowSize, bUsingPep,
                            bModelling, bModel, POT_BL);

}


Int4 AddtoBLPotentialEx (Int4Ptr *ppiBryantPotential, DValNodePtr pdnListPmmds,
						Int2 iModelNum, Boolean bInclusiveWindow,
						Int2 iWindowSize, Boolean bUsingPep,
						Boolean bModelling, Boolean bModel, POT_TYPE pot_type)
{
	return CountStructureContacts(&(*ppiBryantPotential),pdnListPmmds,iModelNum, bInclusiveWindow,iWindowSize,bUsingPep,bModelling,bModel,pot_type);
}

/***********************************************************************************
* BLPotential                                                                  *
*                                                                                  *
* This function writes the BL contacts from memory to a file.                      *
***********************************************************************************/
void DumpBLPotentialEx(Int4Ptr piBryantPotential, FloatLoPtr piHydrophobicity, FILE* fp, Boolean bUsingPep, Boolean bTable, Boolean bAlphaIndex)
{
	Int2 binCounter, rowCounter, colCounter, num_aa;
	Int2 bin = 0, col = 0, row = 0;
	
	CharPtr D_Bin[] = {"0-5","5-6","6-7","7-8","8-9","9-10"};
	CharPtr Alpha_A3_Bin[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","PEP"};
	CharPtr Hydro_A3_Bin[] = {"ALA","VAL","ILE","LEU","PRO","PHE","MET","TRP","GLY","SER","THR","CYS","TYR","GLN","ASN","ASP","GLU","HIS","LYS","ARG","PEP"};
	CharPtr *A3_Bin = NULL;
	Int4 piDim[] = {21,21,6};
	FloatHi value =0.0;

	if(piBryantPotential == NULL) {
		ErrPostEx (SEV_FATAL, 0, 1, "DumpBLContacts:  Pointer to BLContacts is NULL");
		return;
	} else if (fp == NULL) {
		ErrPostEx (SEV_FATAL, 0, 1, "DumpBLContacts:  Pointer to FILE is NULL");
		return;
	}

	/* Write Header to file */
	if(bUsingPep) {
		fprintf(fp,"PEP=TRUE\n\n");
		num_aa = NUM_AA;
	} else {
		fprintf(fp,"PEP=FALSE\n\n");
		num_aa = NUM_AA - 1;
		piDim[0] = piDim[1] = num_aa;
	}

	if(bAlphaIndex) {
		A3_Bin = Alpha_A3_Bin;
	} else {
		A3_Bin = Hydro_A3_Bin;
	}


	/* Write Potential to File */
	for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
		/* First Write the Bin number */
		fprintf (fp, "%-5s", D_Bin[binCounter]);
		bin = binCounter;
		/* Next write the amino acid column header */
		for (colCounter=0; colCounter < num_aa; colCounter++) {
			 fprintf (fp, "%10s", A3_Bin[colCounter]); 
		}
		fprintf (fp, "\n");
		/* Now write the row header */
		for (rowCounter=0; rowCounter < num_aa; rowCounter++) {
			fprintf (fp, "%-5s", A3_Bin[rowCounter]);
			row = GetIndexEx (A3_Bin[rowCounter], bAlphaIndex);
			if(row == -1) { 
				ErrPostEx (SEV_ERROR, 0, 0, "Potential file corrupt.");
				return;
			}
			/* and finally the amino acid contact values */
			for (colCounter=0; colCounter < num_aa; colCounter++) {
				col = GetIndexEx (A3_Bin[colCounter], bAlphaIndex);
				if (col == -1) { 
					ErrPostEx (SEV_ERROR, 0, 0, "Potential file corrupt.");
					return;
				}
				value = (FloatHi) piBryantPotential[Get3DArrayValue(row,col,bin,piDim)];
				if(!bTable) {
					value = value / POTENTIAL_TABLE_SCALE_FACTOR;
				}
					
				fprintf (fp, "%10.6f", (double) value);
			}
			fprintf (fp, "\n");
		}
		fprintf(fp,"\n");
    }

	/* Now write HYDROPHOBICITY terms */
	/* Must be of size bin * row */
	if(piHydrophobicity == NULL) {
		return;
	}
	/* Write the BIN Header */
	fprintf (fp, "\t");
	for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
		/* First Write the Bin number */
		fprintf (fp, " %s\t\t", D_Bin[binCounter]);
	}
	fprintf(fp,"\n");

	/* Write the row header */
	for (rowCounter=0; rowCounter < NUM_AA; rowCounter++) {
		fprintf (fp, "%s\t", A3_Bin[rowCounter]);
		row = GetIndexEx(A3_Bin[rowCounter],bAlphaIndex);
		if (row == -1) {
			ErrPostEx (SEV_ERROR, 3, 3, "Potential file corrupt.");
			return;
		}
		for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
			bin = binCounter;
			value = piHydrophobicity[bin*NUM_AA+row];
			fprintf(fp,"%6.3f\t\t", (double) value);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

/* Dumps a table of counts */
Int2 DumpBLTable(Int4Ptr piBryantPotential, CharPtr pcID, Boolean bUsingPep)
{
	FILE *fp = NULL;
	Boolean bOpen = FALSE;
	Boolean bTable = TRUE;
	Boolean bAlphaIndex = TRUE;

	if(GetBLFile(pcID, &fp, bOpen)) return 1;
	DumpBLPotentialEx(piBryantPotential, NULL, fp, bUsingPep, bTable, bAlphaIndex);
	FileClose(fp);

	return 0;
}

/* Dumps a potential table */
Int2 DumpBLPotential (Int4Ptr piBryantPotential, FloatLoPtr pfHydrophobicity, CharPtr pcID, Boolean bUsingPep)
{
	FILE *fp = NULL;
	Boolean bOpen = FALSE;
	Boolean bTable = FALSE;
	Boolean bAlphaIndex = TRUE;

	if(GetBLFile(pcID, &fp, bOpen)) return 1;
	DumpBLPotentialEx(piBryantPotential, pfHydrophobicity, fp, bUsingPep, bTable, bAlphaIndex);
	FileClose(fp);

	return 0;
}

/***********************************************************************************
* Dump4DPotential                                                                  *
*                                                                                  *
* This function writes the BL contacts from memory to a file.                      *
***********************************************************************************/
void Dump4DPotentialEx(Int4Ptr piBryantPotential, FloatLoPtr piHydrophobicity, FILE* fp, Boolean bUsingPep, Boolean bTable, Boolean bAlphaIndex)
{
	Int2 binCounter, rowCounter, colCounter, num_aa, bin4Dcounter;
	Int2 bin = 0, col = 0, row = 0, ssbin = 0;

	CharPtr D_Bin[] = {"0-5","5-6","6-7","7-8","8-9","9-10"};
	CharPtr SS_Bin[] = {"a-a","a-b","a-c","b-b","b-c","c-c"};
	CharPtr Alpha_A3_Bin[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","PEP"};
	CharPtr Hydro_A3_Bin[] = {"ALA","VAL","ILE","LEU","PRO","PHE","MET","TRP","GLY","SER","THR","CYS","TYR","GLN","ASN","ASP","GLU","HIS","LYS","ARG","PEP"};
	CharPtr *A3_Bin = NULL;
	Int4 piDim[] = {21,21,6,6};
	FloatHi value =0.0;

	if(piBryantPotential == NULL) {
		ErrPostEx (SEV_FATAL, 0, 1, "DumpBLContacts:  Pointer to BLContacts is NULL");
		return;
	} else if (fp == NULL) {
		ErrPostEx (SEV_FATAL, 0, 1, "DumpBLContacts:  Pointer to FILE is NULL");
		return;
	}

	/* Write Header to file */
	if(bUsingPep) {
		fprintf(fp,"PEP=TRUE\n\n");
		num_aa = NUM_AA;
	} else {
		fprintf(fp,"PEP=FALSE\n\n");
		num_aa = NUM_AA - 1;
		piDim[0] = piDim[1] = num_aa;
	}

	if(bAlphaIndex) {
		A3_Bin = Alpha_A3_Bin;
	} else {
		A3_Bin = Hydro_A3_Bin;
	}


	/* Write Potential to File */
	for (bin4Dcounter=0; bin4Dcounter < NUM_BINS4D; bin4Dcounter++) {
		ssbin = bin4Dcounter;
		for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
			/* First Write the ids for bins */
			fprintf (fp, "%-4s %-4s\t", SS_Bin[bin4Dcounter], D_Bin[binCounter]);
			bin = binCounter;
			/* Next write the amino acid column header */
			for (colCounter=0; colCounter < num_aa; colCounter++) {
				fprintf (fp, "%12s\t", A3_Bin[colCounter]);
			}
			fprintf (fp, "\n");
			/* Now write the row header */
			for (rowCounter=0; rowCounter < num_aa; rowCounter++) {
				fprintf (fp, "%-4s\t", A3_Bin[rowCounter]);
				row = GetIndexEx (A3_Bin[rowCounter], bAlphaIndex);
				if(row == -1) {
					ErrPostEx (SEV_ERROR, 0, 0, "Potential file corrupt.");
					return;
				}
				/* and finally the amino acid contact values */
				for (colCounter=0; colCounter < num_aa; colCounter++) {
					col = GetIndexEx (A3_Bin[colCounter], bAlphaIndex);
					if (col == -1) {
						ErrPostEx (SEV_ERROR, 0, 0, "Potential file corrupt.");
						return;
					}
					value = (FloatHi) piBryantPotential[Get4DArrayValue(row,col,bin,ssbin,piDim)];
					if(!bTable) {
						value = value / POTENTIAL_TABLE_SCALE_FACTOR;
					}
					fprintf (fp, "%12.9f\t", (double) value);
				}
			fprintf (fp, "\n");
			}
		fprintf(fp,"\n");
       }
	}

	/* Now write HYDROPHOBICITY terms */
	/* Must be of size bin * row */
	if(piHydrophobicity == NULL) {
		return;
	}
	/* Write the BIN Header */
	fprintf (fp, "\t");
	for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
		/* First Write the Bin number */
		fprintf (fp, " %s\t\t", D_Bin[binCounter]);
	}
	fprintf(fp,"\n");

	/* Write the row header */
	for (rowCounter=0; rowCounter < NUM_AA; rowCounter++) {
		fprintf (fp, "%s\t", A3_Bin[rowCounter]);
		row = GetIndexEx(A3_Bin[rowCounter],bAlphaIndex);
		if (row == -1) {
			ErrPostEx (SEV_ERROR, 3, 3, "Potential file corrupt.");
			return;
		}
		for (binCounter=0; binCounter < NUM_BINS; binCounter++) {
			bin = binCounter;
			value = piHydrophobicity[bin*NUM_AA+row];
			fprintf(fp,"%6.3f\t\t", (double) value);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}


/* Dumps a table of counts */
Int2 Dump4DTable(Int4Ptr piBryantPotential, CharPtr pcID, Boolean bUsingPep)
{
	FILE *fp = NULL;
	Boolean bOpen = FALSE;
	Boolean bTable = TRUE;
	Boolean bAlphaIndex = TRUE;

	if(GetBLFile(pcID, &fp, bOpen)) return 1;
	Dump4DPotentialEx(piBryantPotential, NULL, fp, bUsingPep, bTable, bAlphaIndex);
	FileClose(fp);

	return 0;
}




/**
 *
 * ComputeBryantPotential (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds,
 *                         Int2 iModelNum, Boolean bInclusiveWindow,
 *                         Int2 iWindowSize, Int4Ptr piBryantPotential,
 *                         Boolean bUsingPep, Boolean bDetailedList) {
 *
 * Takes a list of molecules specified by "pdnListPmmds" and "iModelNum",
 * and computes a BL potential for each residue in each molecule. The
 * "bInclusiveWindow" parameter determines whether an inclusive or
 * exclusive window is used and the window size is specified using the 
 * variable "iWindowSize". "piBryantPotential" is a pointer to the BL
 * potential function array and "bUsingPep" specifies whether "PEPs" are
 * to be included in the calculation.
 *
 * A new adjacency list is created to store the results of these
 * calculations and "*ppdnResult" is set to point to the head of this
 * list. The adjacency list is essentially a list of lists. The main
 * list contains a DValNode for every residue. The enumerated residue
 * number is stored in the choice field of the DValNode, whereas the
 * data.ptrvalue field points to a sublist of DValNodes that point to
 * AdjListNodes. The first AdjListNode in each sublist points to the PMGD
 * of the current residue and stores the total BL potential computed for
 * that residue. If "bDetailedList" is FALSE, each sublist will only
 * contain this first AdjListNode. If "bDetailedList" is TRUE, an
 * AdjListNode is also created for each residue that interacts with the
 * current residue. The potentials stored in these subsequent AdjListNodes
 * are the individual interaction potentials that sum to give the total
 * potential for a given residue.
 *
 * The total Bryant-Lawrence potential for the list of molecules is also
 * computed and is stored in the global variable "TotalPotential".
 *
 */


void ComputeBryantPotentialEx (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds, Int2 iModelNum, 
							 Boolean bInclusiveWindow, Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList, Boolean bModelling, Boolean bModel) {

    Boolean bExclusiveWindow = (Boolean) !bInclusiveWindow;

    DValNodePtr pdnmgHeadPsBList = NULL;
    DValNodePtr pdnmgHere = NULL;
    DValNodePtr pdnHeadSublist;
    DValNodePtr pdnTest = NULL;
    PALN palnSub = NULL;
    ValNodePtr pvnAtomList, pvnHere;
    AdjListNodePtr paln;

    PWS pwsHead = NULL;
    PWS pwsThis = NULL;
    vec middle;
    Int2 iResId1, iResId2, iMolId1, iMolId2, res1, res2, bin;

    PMMD pmmd1, pmmd2;
    PMGD pmgd1, pmgd2;
    PMGD pmgdSub;
    PMAD pmad1, pmad2;
    PALD pald1, pald2;

    FloatHi potential, ptnl;
    Char myresname[2];
    TotalPotential = 0.0;

    pdnListPepMGDs = NULL;

    CreatePsBList (&pdnmgHeadPsBList, pdnListPmmds, iModelNum, bUsingPep, bModelling, FALSE);
    pdnmgHere = pdnmgHeadPsBList;

    /* Create a world containing all of the pseudo-beta carbons */
    while (pdnmgHere != NULL) {
		pmgd1 = (PMGD) pdnmgHere->data.ptrvalue;
		pwsThis = AddtoWorld (pwsHead, iModelNum, (PFB) pmgd1);
		if (pwsHead == NULL) {
			pwsHead = pwsThis;
		}
		pdnmgHere = pdnmgHere->next;
    }
	InstantiateWorld (1, pwsHead);

    pdnmgHere = pdnmgHeadPsBList;

    /* For each pseudo-beta carbon */
    while (pdnmgHere != NULL) {
	    iResId1=pdnmgHere->choice;
	    pmgd1 = (PMGD)(pdnmgHere->data.ptrvalue);
		pmad1 = (PMAD)(pmgd1->pvnmaAHead->data.ptrvalue);
	    pald1 = GetAtomLocs (pmad1, iModelNum);
	    pmgd1 = (PMGD)(pmad1->pfbParent);
	    pmmd1 = (PMMD)(pmgd1->pfbParent);

	    iMolId1 = ((PDNMM)(pmmd1->pdnmmLink))->choice;

		palnSub=NULL;
		/* check if already in the list or what */
		pdnTest=*ppdnResult;
		while (pdnTest!=NULL) {
			if (pdnTest->choice==abs(iResId1)) {
				palnSub=(PALN)(((DValNodePtr)(pdnTest->data.ptrvalue))->data.ptrvalue);
				pmgdSub=palnSub->pmgd;
				if ((PMMD)(pmgdSub->pfbParent)==pmmd1) {
					/* same residue, same molecule */
					/* add to total potential */
					break;
				} else {
					palnSub=NULL;
				}
			}
			pdnTest=pdnTest->next;
		}  

		if (palnSub==NULL) {
			/* Add the residue to the head of the sublist */
			paln = (AdjListNodePtr) MemNew (sizeof (AdjListNode));
			paln->pmgd = pmgd1;
			pdnHeadSublist = NULL;
			DValNodeAddPointer (&pdnHeadSublist, 0, paln);
		} else {
			/* use existing head sub list */
			pdnHeadSublist=(DValNodePtr)(pdnTest->data.ptrvalue);
		}

	    potential = 0.0;

	    middle[0] = pald1->pflvData[0];
	    middle[1] = pald1->pflvData[1];
	    middle[2] = pald1->pflvData[2];

	    pvnAtomList = FindAtomsIn (pwsHead, middle, 10.0);
	    pvnHere = pvnAtomList;

	    while (pvnHere != NULL) {

			pald2 = (PALD)(pvnHere->data.ptrvalue);
			pmad2 = (PMAD)(pald2->pfbParent);
			pmgd2 = (PMGD)(pmad2->pfbParent);
			pmmd2 = (PMMD)(pmgd2->pfbParent);

			iMolId2 = ((PDNMM)(pmmd2->pdnmmLink))->choice;
			iResId2 = pmad2->iIndex;
			
			if (
				(bExclusiveWindow &&
				  ((iMolId2 != iMolId1) || 
				  ((iMolId2 == iMolId1) &&
				  ((abs (iResId2) > (abs(iResId1) + iWindowSize)) || 
				   (abs (iResId2) < (abs(iResId1) - iWindowSize)))))) 
				   ||
				(bInclusiveWindow &&
				  (iMolId1 == iMolId2) && 
				  (iResId1 != iResId2) &&
				  (abs (iResId2) >= (abs(iResId1) - iWindowSize)) &&
				  (abs (iResId2) <= (abs(iResId1) + iWindowSize)))
				  ) {

				if (iResId1 > 0) {
					if(!bModel) {
						myresname[0]=GetAAFromIDict(pmgd1);
					} else {
						myresname[0]=(Char)toupper(pmgd1->bReserved);
					}
					myresname[1]='\0';
					res1 = GetIndex (myresname);
				} else {
					res1 = MK_PEP;
				}

				if (iResId2 > 0) {
					if(!bModel) {
						myresname[0]=GetAAFromIDict(pmgd2);
					} else {
						myresname[0]=(Char)toupper(pmgd2->bReserved);
					}
					myresname[1]='\0';
					res2 = GetIndex (myresname);
				} else {
					res2 = MK_PEP;
				}

				bin = GetBinNumber (pald1, pald2);
				if (bin!=-1)
					ptnl = GetBryantPotential (piBryantPotential, res1, res2, bin);
				else
					ptnl = 0.0;
			   
				if (bDetailedList) {
					paln = (AdjListNodePtr) MemNew (sizeof (AdjListNode));
					paln->pmgd = pmgd2;
					paln->potential = ptnl/2.0;
					DValNodeAddPointer (&pdnHeadSublist, 0, paln);
				}

				potential += ptnl/2.0;
				TotalPotential += ptnl/2.0;
			}
			pvnHere = pvnHere->next;
	    }
	    ValNodeFree(pvnAtomList);

	    /* if didn't find match, new residue */
		if (palnSub==NULL) {
	    	((PALN)(pdnHeadSublist->data.ptrvalue))->potential = potential;
			DValNodeAddPointer (ppdnResult, (Int2) abs(iResId1), pdnHeadSublist);
		} else {
	    	palnSub->potential+=potential;
		}	
		pdnmgHere = pdnmgHere->next;
    }

    FreePsBList (&pdnmgHeadPsBList);
	FreeAllWorlds ();
}


void ComputeBryantPotential (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds, Int2 iModelNum,
							 Boolean bInclusiveWindow, Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList) 
{
	Boolean bModelling = FALSE, bModel = FALSE;
	ComputeBryantPotentialEx (&(*ppdnResult), pdnListPmmds, iModelNum, bInclusiveWindow,
                             iWindowSize, piBryantPotential, bUsingPep, bDetailedList, bModelling, bModel); 
}

void ComputeBryantPotentialModel (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds, Int2 iModelNum,
							 Boolean bInclusiveWindow, Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList) 
{
	Boolean bModelling = TRUE, bModel = TRUE;
	ComputeBryantPotentialEx (&(*ppdnResult), pdnListPmmds, iModelNum, bInclusiveWindow,
                             iWindowSize, piBryantPotential, bUsingPep, bDetailedList, bModelling, bModel);
}

void ComputeBryantPotentialTemplate (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds, Int2 iModelNum,
							 Boolean bInclusiveWindow, Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList) 
{
	Boolean bModelling = TRUE, bModel = FALSE;
	ComputeBryantPotentialEx (&(*ppdnResult), pdnListPmmds, iModelNum, bInclusiveWindow,
                             iWindowSize, piBryantPotential, bUsingPep, bDetailedList, bModelling, bModel);
}



/* only for use by Crease Energy */
void FreePepMGDs(void)
{
    if (pdnListPepMGDs != NULL) {
	/* corrected by howie -- MemFree -> FreeMGD */
	pdnListPepMGDs = DValNodeFreeData (pdnListPepMGDs, (pFreeFunc) FreeMGD3);
    }
}

/**
 *
 * FreeAdjList (DValNodePtr *ppdnHead):
 * Frees memory allocated to *ppdnHead and sets it to NULL.
 *
 */

void FreeAdjList (DValNodePtr *ppdnHead) {

    DValNodePtr pdnHere = *ppdnHead;
    DValNodePtr pdnNext, pdnSublistHere, pdnSublistNext;

    while (pdnHere != NULL) {

	pdnNext = pdnHere->next;
	pdnSublistHere = (DValNodePtr) pdnHere->data.ptrvalue;
	MemFree (pdnHere);

	while (pdnSublistHere != NULL) {

	    pdnSublistNext = pdnSublistHere->next;
	    MemFree ((PALN)(pdnSublistHere->data.ptrvalue));
	    MemFree (pdnSublistHere);
	    pdnSublistHere = pdnSublistNext;

	}

	pdnHere = pdnNext;

    }

    FreePepMGDs();

    *ppdnHead = NULL;

}



/**
 *
 * GetTotalPotential (void):
 * Returns the total potential of the current molecule.
 *
 */

FloatHi GetTotalPotential (void) {
    return TotalPotential;
}



/**
 *
 * PrintAdjList (DValNodePtr pdnHead, Int2 iModelNum, Boolean bVerbose
 *               Boolean bTemp):
 *
 * Prints out the contents of the adjacency list pointed to by "pdnHead".
 * The "bTemp" parameter specifies whether temperature factor information
 * should be printed along with the potentials.
 *
 */

void PrintAdjList (DValNodePtr pdnHead, Int2 iModelNum, Boolean bVerbose,
                   Boolean bTemp) {

    DValNodePtr pdnHere = pdnHead;
    DValNodePtr pdnSublist,pdnTest;
    PVNMA pvnmaHere;

    PMSD pmsd;
    PMMD pmmd;
    PMGD pmgd;
    PMAD pmad;
    PALD pald=NULL;

    FloatHi ptnl,TestSum;
    Int2 resNum, numAtoms;
    FloatLo b, MW;
    Char ss, data[4];
    CharPtr str = data;
    Char aa3 [] = MK_AA_3LETTER;
    Char myresname[2];

    if (pdnHead == NULL) {
	ErrPostEx (SEV_INFO, 7, 1, "Empty Adjacency list.");
	return;
    }

    pdnSublist = (DValNodePtr)(pdnHead->data.ptrvalue);
    pmgd = ((AdjListNodePtr)(pdnSublist->data.ptrvalue))->pmgd;
    pmmd = (PMMD)(pmgd->pfbParent);
    pmsd = (PMSD)(pmmd->pfbParent);

    if (bVerbose || bTemp) {
	printf ("> %s, model #%d\n", pmsd->pcPDBName, iModelNum);
	if (bTemp) {
	    printf ("Chain\tResNum\tAA\tMW\t\tPotential\tB factor\tSecStru\n");
	} else {
	    printf ("Chain\tResNum\tAA\tMW\t\tPotential\tSecStru\n");
	}
    } else {
	printf ("Model #%d", iModelNum);
    }

    while (pdnHere != NULL) {

	pdnSublist = (DValNodePtr)(pdnHere->data.ptrvalue);
	pmgd = ((AdjListNodePtr)(pdnSublist->data.ptrvalue))->pmgd;
	ptnl = ((AdjListNodePtr)(pdnSublist->data.ptrvalue))->potential;

	TestSum=0.0;
	pdnTest=pdnSublist->next;
	if (pdnTest!=NULL) {
		while (pdnTest) {
			TestSum += ((AdjListNodePtr)(pdnTest->data.ptrvalue))->potential;
			pdnTest=pdnTest->next;
		}
		if (TestSum-ptnl>0.00001)
			ErrPostEx(SEV_ERROR,1,0,"Potential Integrity Error");
	}

	if (bVerbose || bTemp) {

	    pmmd = (PMMD)(pmgd->pfbParent);
	    resNum = ((PDNMG)(pmgd->pdnmgLink))->choice;
	    myresname[0]=GetAAFromIDict(pmgd);
	    myresname[1]='\0';

	    StringNCpy (str, &(aa3[GetIndex(myresname)*3]), 3);
	    str[3] = '\0';

	    MW = MWaalist [GetResIdxFromMGD(aalist,pmgd)];

	    if (pmgd->bNCBISecStru & (Byte) SS_HELIX) {
		ss = 'H';
	    } else if (pmgd->bNCBISecStru & (Byte) SS_STRAND) {
		ss = 'E';
	    } else if (pmgd->bNCBISecStru & (Byte) SS_TURN) {
		ss = 'T';
	    } else {
		ss = 'C';
	    }

	    if (bTemp) {

		numAtoms = 0;
		b = 0.0;

		pvnmaHere = pmgd->pvnmaAHead;

		while (pvnmaHere != NULL) {

		    if ((pmad = (PMAD)(pvnmaHere->data.ptrvalue)) != NULL) {
			pald = GetAtomLocs (pmad, iModelNum);
		    }

		    if ((pald != NULL) && (pald->iFloatNo == 9)) {

			ErrPostEx (SEV_ERROR, 8, 1, "Formula not implemented yet.");

			/* The funky crystal formula goes here... */
			b += 0;

			numAtoms++;

		    } else if ((pald != NULL) && (pald->iFloatNo >= 4)) {

			b += pald->pflvData[4];
			numAtoms++;

		    }

		    pvnmaHere = pvnmaHere->next;

		}

		if (numAtoms != 0) {
		    b = b / numAtoms;
		}

		printf ("%c\t%d\t%s\t%f\t%f\t%f\t%c\n",
		    pmmd->pcMolName[0], resNum, str, MW, ptnl, b, ss);

	    } else {

		printf ("%c\t%d\t%s\t%f\t%f\t%c\n",
		    pmmd->pcMolName[0], resNum, str, MW, ptnl, ss);

	    }

	} else {
	    printf ("\t%f", ptnl);
	}
	pdnHere = pdnHere->next;

    }

    if (bVerbose) {
	printf ("\nTotal Potential: %f\n\n", GetTotalPotential());
    } else {
	printf ("\n");
    }

}



/**
 *
 * CreateListPmmds (CharPtr chainList, PMSD pmsd):
 *
 * Returns a pointer to a list of DValNodes, where the data.ptrvalue field
 * of each each node is a PMMD pointing to a protein chain specified by
 * "chainList".
 *
 */

DValNodePtr CreateListPmmds (CharPtr chainList, PMSD pmsd) {

    DValNodePtr pdnResult = NULL;
    Boolean bFoundChain;

    PDNMM pdnmmHere; 
    PMMD pmmd=NULL;

    Int2 choice=0;

    if (chainList == NULL) {
 
			/* Default case: Find first molecule in structure */
		pmmd = (PMMD) ((pmsd->pdnmmHead)->data.ptrvalue);
		choice = (pmsd->pdnmmHead)->choice;
		DValNodeAddPointer (&pdnResult, choice, pmmd);

    } else {

    		/* Mark all mmd's in structure as unvisited */
		pdnmmHere = pmsd->pdnmmHead;
		while (pdnmmHere != NULL) {
			((PMMD)(pdnmmHere->data.ptrvalue))->bReserved = (Byte) FALSE;
			pdnmmHere = pdnmmHere->next;
		}

		while (chainList[0] != '\0') {

			pdnmmHere = pmsd->pdnmmHead;
			bFoundChain = FALSE;

			/* Case-sensitive test */
			while ((pdnmmHere != NULL) && (!bFoundChain)) {

				pmmd = (PMMD)(pdnmmHere->data.ptrvalue);
				choice = pdnmmHere->choice;

				if (((StringNCmp (chainList, pmmd->pcMolName, 1)) == 0) &&
					(IsProtein ((PFB)(pmmd)))) {
					bFoundChain = TRUE;
				}

				pdnmmHere = pdnmmHere->next;

			}

			pdnmmHere = pmsd->pdnmmHead;

			/* Case-insensitive test */
			while ((pdnmmHere != NULL) && (!bFoundChain)) {

				pmmd = (PMMD)(pdnmmHere->data.ptrvalue);
				choice = pdnmmHere->choice;

				if (((StringNICmp (chainList, pmmd->pcMolName, 1)) == 0) &&
					(IsProtein ((PFB)(pmmd)))) {
					bFoundChain = TRUE;
				}

				pdnmmHere = pdnmmHere->next;

			}

			if ((bFoundChain) && (pmmd->bReserved)) {
				ErrPostEx (SEV_ERROR, 2, 1, "Multiple references to same chain.");
 			} else if ((bFoundChain) && !(pmmd->bReserved)) {
				DValNodeAddPointer (&pdnResult, choice, pmmd);
				pmmd->bReserved = (Byte) TRUE;
			} else {
				ErrPostEx (SEV_ERROR, 2, 2, "Chain not found.");
			}

			chainList++;

		}
    }

    return pdnResult;

}



/**
 *
 * FreeListPmmds (DValNodePtr *ppdnHead):
 *
 * Frees the list of DValNodes created by "CreateListPmmds" and sets
 * "*ppdnHead" to NULL;
 *
 */

void FreeListPmmds (DValNodePtr *ppdnHead) {

    if (*ppdnHead != NULL) {
	DValNodeFree (*ppdnHead);
    }

    *ppdnHead = NULL;

}



/**
 *
 * IsValidModel (PMSD pmsd, Int2 iModelNum):
 *
 * Takes a PMSD and a model number and returns TRUE/FALSE depending on 
 * whether the structure pointed to by "pmsd" has a model with "iModelNum"
 * as its ModelId.
 *
 */

Boolean IsValidModel (PMSD pmsd, Int2 iModelNum) {

    PDNML pdnmlHere = pmsd->pdnmlModels;

    while (pdnmlHere != NULL) {

	if (iModelNum == pdnmlHere->choice) {
	    return TRUE;
	}

	pdnmlHere = pdnmlHere->next;

    }

    return FALSE;

}



/**
 *
 * CreateListOfModels (CharPtr pcModelList, PMSD pmsd):
 *
 * Takes a string of comma-separated integers and parses it. If each of
 * the integers in the list is a valid model in structure "pmsd", the
 * function returns a pointer to a list of ValNodes that stores these
 * model numbers.
 *
 */

ValNodePtr CreateListOfModels (CharPtr pcModelList, PMSD pmsd) {

    ValNodePtr pvnResult = NULL;
    Char buffer [MAX_STRING_LENGTH];
    Int4 iModelNum, index = 0;

    while (pcModelList[0] != '\0') {

	if ((pcModelList[0] == ',') || (pcModelList[1] == '\0')) {

	    if ((pcModelList[1] == '\0') && (pcModelList[0] != ',')) {
		buffer [index] = pcModelList[0];
		index++;
	    }

	    buffer [index] = '\0';
	    iModelNum = (Int4) strtol (buffer, NULL, 10);
	    index = 0;

	    if (IsValidModel (pmsd, (Int2)iModelNum)) {
		ValNodeAddInt (&pvnResult, -1, iModelNum);
	    } else {
		ErrPostEx (SEV_ERROR, 6, 1,
		    "Invalid model number (%d).", iModelNum);
	    }

	} else {

	    buffer [index] = pcModelList[0];
	    index++;

	}

	pcModelList++;

    }

    return pvnResult;

}



/**
 *
 * FreeListOfModels (ValNodePtr *ppvnHead):
 *
 * Frees the list of ValNodes created by "CreateListOfModels" and sets
 * "*ppvnHead" to NULL;
 *
 */

void FreeListOfModels (ValNodePtr *ppvnHead) {

    if (*ppvnHead != NULL) {
	ValNodeFree (*ppvnHead);
    }

    *ppvnHead = NULL;

}



/**
 *
 * LoadZhangPotential (Int4Ptr *ppiZhangPtnl, FILE *fp):
 *
 * Creates an array to store a Zhang potential in memory, fills this array
 * with information read in from "fp" and makes "*ppiZhangPtnl" point to
 * this array. Returns 0 on success.
 *
 * The 2-D table of floats that is the Zhang potential is stored in memory
 * as a 1-D array of integers. POTENTIAL_TABLE_SCALE_FACTOR is used to
 * convert the floats to integers and vice versa.
 *
 */

TrajErr LoadZhangPotential (Int4Ptr *ppiZhangPtnl, FILE *fp) {

    int temp;
    double tempFloat;
    Int2 row, col, tempInt;
    Char buffer [MAX_STRING_LENGTH + 1];
    CharPtr str = buffer;
    Int4Ptr potential;

    fscanf (fp, "%s", str);

    if (StringICmp (str, "eij") != 0) {
	ErrPostEx (SEV_ERROR, 5, 1, "Potential file corrupt.");
	return ERR_FAIL;
    }

    for (col=1; col <= NUM_ZHANG; col++) {

	fscanf (fp, "%d", &temp);
	tempInt = (Int2) temp;
	if (tempInt != col) {
	    ErrPostEx (SEV_ERROR, 5, 2, "Potential file corrupt.");
	    return ERR_FAIL;
	}

    }

    potential = (Int4Ptr) MemNew ((NUM_ZHANG + 1) * (NUM_ZHANG + 1) * sizeof (Int4));

    for (row=1; row <= NUM_ZHANG; row++) {

	fscanf (fp, "%d", &temp);
	tempInt = (Int2) temp;
	if (tempInt != row) {
	    ErrPostEx (SEV_ERROR, 5, 3, "Potential file corrupt.");
	    return ERR_FAIL;
	}

	for (col=1; col <= NUM_ZHANG; col++) {

	    fscanf (fp, "%lf", &tempFloat);

	    potential [row * NUM_ZHANG + col]
		= (Int4) floor (tempFloat * POTENTIAL_TABLE_SCALE_FACTOR);

	}

    }

    *ppiZhangPtnl = potential;
    return ERR_SUCCESS;

}



/**
 *
 * FreeZhangPotential (Int4Ptr *ppiZhangPtnl):
 *
 * Frees memory allocated to the Zhang potential function array and sets
 * "*ppiZhangPtnl" to NULL.
 *
 */

void FreeZhangPotential (Int4Ptr *ppiZhangPtnl) {

    if (*ppiZhangPtnl != NULL) {
	MemFree (*ppiZhangPtnl);
    }

    *ppiZhangPtnl = NULL;

}



/**
 *
 * LoadZhangAtmList (DValNodePtr *ppdnResult, FILE *fp):
 *
 * Creates a list that is used to store indexing information for a Zhang
 * potential, fills this array with information read in from "fp" and
 * makes "*ppdnResult" point to this list. Returns 0 on success.
 *
 */

Int2 LoadZhangAtmList (DValNodePtr *ppdnResult, FILE *fp) {

    Char buffer [MAX_STRING_LENGTH + 1];
    Char temp [3];

    CharPtr str = buffer;
    CharPtr pcTemp = temp;

    ZhangAtmNodePtr pzanThis;

    while (!feof (fp)) {

	if (FileGets (str, MAX_STRING_LENGTH, fp)) {

	    pzanThis = (ZhangAtmNodePtr) MemNew (sizeof (ZhangAtmNode));

	    pzanThis->cResName = str[0];
	    StringNCpy (pzanThis->cAtomName, &str[2], 4);
	    pzanThis->cAtomName[4] = '\0';
	    StringNCpy (pcTemp, &str[7], 2);
	    pcTemp[2] = '\0';
	    pzanThis->iIndex = (Int2) atoi (pcTemp);

	    DValNodeAddPointer (ppdnResult, 0, pzanThis);

	}

    }

    return 0;

}



/**
 *
 * FreeZhangAtmNode (ZhangAtmNodePtr pzan):
 * Frees memory allocated to "pzan".
 *
 */

void LIBCALLBACK FreeZhangAtmNode (Pointer ptr) {
	ZhangAtmNodePtr pzan;

	pzan=(ZhangAtmNodePtr)ptr;
    if (pzan != NULL) {
	MemFree (pzan);
    }

}



/**
 *
 * FreeZhangAtmList (DValNodePtr *ppdnHead):
 * Frees memory allocated to *ppdnHead and sets it to NULL.
 *
 */

void FreeZhangAtmList (DValNodePtr *ppdnHead) {

    if (*ppdnHead != NULL) {
	DValNodeFreeData (*ppdnHead, (pFreeFunc) FreeZhangAtmNode);
    }

    *ppdnHead = NULL;

}



/**
 *
 * GetZhangIndex (Char cResName, CharPtr pcAtomName, DValNodePtr pdnAtmList):
 *
 * Takes atom type information and returns an integer value that can be
 * used to index a Zhang potential function array.
 * 
 * cResName:  One-letter amino acid code
 * pcAtomName:  4 character PDB atom name
 * pdnAtmList:  Pointer to a Zhang AtmList
 *
 */

Int2 GetZhangIndex (Char cResName, CharPtr pcAtomName, DValNodePtr pdnAtmList) {

    ZhangAtmNodePtr pzanHere;

    while (pdnAtmList != NULL) {

	pzanHere = (ZhangAtmNodePtr)(pdnAtmList->data.ptrvalue);

	if ((pzanHere->cResName == cResName) &&
	    (StringCmp (pzanHere->cAtomName, pcAtomName) == 0)) {
	    return pzanHere->iIndex;
	}

	pdnAtmList = pdnAtmList->next;

    }

    return -1;

}



/**
 *
 * GetZhangPotential (PMAD pmad1, PMAD pmad2,
 *                    Int4Ptr piZhangPtnl, DValNodePtr pdnAtmList):
 *
 * Takes two PMADs, a Zhang potential and a Zhang AtmList and returns
 * the value of the Zhang potential between atoms pmad1 and pmad2.
 *
 */

FloatHi GetZhangPotential (PMAD pmad1, PMAD pmad2,
                           Int4Ptr piZhangPtnl, DValNodePtr pdnAtmList) {

    Char cRes1, cRes2;
    CharPtr pcAtomName1, pcAtomName2;
    Int2 iIndex1, iIndex2;

    cRes1 = GetAAFromIDict((PMGD)(pmad1->pfbParent));
    cRes2 = GetAAFromIDict((PMGD)(pmad2->pfbParent));

    pcAtomName1 = pmad1->pcAName;
    pcAtomName2 = pmad2->pcAName;

    iIndex1 = GetZhangIndex (cRes1, pcAtomName1, pdnAtmList);
    iIndex2 = GetZhangIndex (cRes2, pcAtomName2, pdnAtmList);

    if (iIndex1 <= iIndex2) {
	return piZhangPtnl [iIndex1 * NUM_ZHANG + iIndex2]
	    / POTENTIAL_TABLE_SCALE_FACTOR;
    } else {
	return piZhangPtnl [iIndex2 * NUM_ZHANG + iIndex1]
	    / POTENTIAL_TABLE_SCALE_FACTOR;
    }

}



/**
 *
 * ComputeZhangPotential (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds,
 *                        Int2 iModelNum, Boolean bInclusiveWindow,
 *                        Int2 iWindowSize, Int4Ptr piZhangPtnl,
 *                        DValNodePtr pdnAtmList, Boolean bDetailedList):
 *
 * Computes the total Zhang potential for a list of molecules and stores
 * this value in the global variable "TotalPotential". Also creates an
 * adjacency list of the type created by "ComputeBryantPotential".
 *
 * ppdnResult:  Pointer to the newly created adjacency list
 * pdnListPmmds:  List of input molecules
 * iModelNum:  Model number
 * bInclusiveWindow:  True if using an inclusive window, false if exclusive
 * iWindowSize:  Size of include/exclude window
 * piZhangPtnl:  Pointer to a Zhang potential
 * pdnAtmList:  Pointer to a Zhang AtmList
 *
 */

void ComputeZhangPotential (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds,
                            Int2 iModelNum, Boolean bInclusiveWindow,
                            Int2 iWindowSize, Int4Ptr piZhangPtnl,
                            DValNodePtr pdnAtmList, Boolean bDetailedList) {

    ValNodePtr /*pvnWorldHead,*/ pvnAtomList, pvnHere;
    DValNodePtr pdnHere, pdnHeadSublist;

    AdjListNodePtr paln;

    PWS pwsHead = NULL;
    PWS pwsThis = NULL;
    vec middle;
    Int2 iMolId1, iMolId2, iResId1, iResId2, numAtoms;

    PMMD pmmd1, pmmd2;
    PMGD pmgd1, pmgd2;
    PMAD pmad1, pmad2;
    PALD pald1, pald2;

    PDNMG pdnmgHere;
    PVNMA pvnmaHere;

    FloatHi potential,ptnl;

    /* For each molecule in the list */
    pdnHere = pdnListPmmds;
    while (pdnHere != NULL) {
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
				    /* If the atom is not a hydrogen and it has a recognized
				     * PDB atom name and it has location data, then... */
		 		    if ((((PVNMA)(pmad1->pvnmaLink))->choice != ATOMIC_NO_H) &&
							(GetZhangIndex (GetAAFromIDict(pmgd1), pmad1->pcAName,
			  		  pdnAtmList) != -1) &&
							(pald1->pflvData != NULL)) {
							/* ... Add it to the world of atoms */
							pwsThis = AddtoWorld (pwsHead,iModelNum,(PFB)(pmad1));
							if (pwsHead == NULL) pwsHead = pwsThis;
		    		}
					}
					pvnmaHere = pvnmaHere->next;
		    }
		    pdnmgHere = pdnmgHere->next;
			}
			pdnHere = pdnHere->next;
    }
/*    pvnWorldHead =*/ InstantiateWorld (1, pwsHead);
    TotalPotential = 0.0;
    /* For each molecule in the list */
    while (pdnListPmmds != NULL) {
			pmmd1 = (PMMD)(pdnListPmmds->data.ptrvalue);
			pdnmgHere = pmmd1->pdnmgHead;
			iMolId1 = pdnListPmmds->choice;
			/* For each residue in the molecule */
			while (pdnmgHere != NULL) {
	  	  pmgd1 = (PMGD)(pdnmgHere->data.ptrvalue);
	    	pvnmaHere = pmgd1->pvnmaAHead;
		    iResId1 = pdnmgHere->choice;
		    potential = 0.0;
	  	  numAtoms = 0;
		    /* Add the residue to the head of the sublist */
        paln = (AdjListNodePtr) MemNew (sizeof (AdjListNode));
        paln->pmgd = pmgd1;
        pdnHeadSublist = NULL;
        DValNodeAddPointer (&pdnHeadSublist, 0, paln);
	  	  /* For each atom */
		    while (pvnmaHere != NULL) {
					pmad1 = (PMAD)(pvnmaHere->data.ptrvalue);
					if ((pald1 = GetAtomLocs (pmad1, iModelNum)) != NULL) {
		  		  /* If the atom is not a hydrogen and it has a recognized
		  		   * PDB atom name and it has location data, then... */
		 		    if ((((PVNMA)(pmad1->pvnmaLink))->choice != ATOMIC_NO_H) &&
							(GetZhangIndex (GetAAFromIDict(pmgd1), pmad1->pcAName,
					    pdnAtmList) != -1) &&
							(pald1->pflvData != NULL)) {
							numAtoms++;
							middle[0] = pald1->pflvData[0];
							middle[1] = pald1->pflvData[1];
							middle[2] = pald1->pflvData[2];
							pvnAtomList = FindAtomsIn (pwsHead, middle, 6.0);
							pvnHere = pvnAtomList;
							while (pvnHere != NULL) {
						    pald2 = (PALD)(pvnHere->data.ptrvalue);
						    pmad2 = (PMAD)(pald2->pfbParent);
						    pmgd2 = (PMGD)(pmad2->pfbParent);
						    pmmd2 = (PMMD)(pmgd2->pfbParent);
						    iMolId2 = ((PDNMM)(pmmd2->pdnmmLink))->choice;
						    iResId2 = ((PDNMG)(pmgd2->pdnmgLink))->choice;
						    ptnl=0.0;
						    if (bInclusiveWindow &&
									(iMolId1 == iMolId2) &&
									(iResId1 != iResId2) &&
									(iResId2 >= (iResId1 - iWindowSize)) &&
									(iResId2 <= (iResId1 + iWindowSize))) {
									ptnl = GetZhangPotential (pmad1,
								    pmad2, piZhangPtnl, pdnAtmList);
			  			  }
			  			  else if (!bInclusiveWindow &&
									((iMolId1 != iMolId2) ||
									((iMolId1 == iMolId2) &&
									((iResId2 < (iResId1 - iWindowSize)) ||
									 (iResId2 > (iResId1 + iWindowSize)))))) {
									ptnl = GetZhangPotential (pmad1,
								    pmad2, piZhangPtnl, pdnAtmList);
			    			}
				        potential+=ptnl;					
								if (bDetailedList && ptnl!=0.0) {
									paln = (AdjListNodePtr) MemNew (sizeof (AdjListNode));
									paln->pmgd = pmgd2;
									paln->potential = ptnl/2.0;
									DValNodeAddPointer (&pdnHeadSublist, 0, paln);
								}
						    pvnHere = pvnHere->next;
							}
							ValNodeFree(pvnAtomList);
		    		}
					}
					pvnmaHere = pvnmaHere->next;
		    }
/* do not normalize by number of atoms */
/*	    if (numAtoms != 0) {
		potential = potential / numAtoms;
	    }
*/
	  	  TotalPotential += potential;
	 	   /* Add the sublist to the adjacency list */
	  	  DValNodeAddPointer (ppdnResult, iResId1, pdnHeadSublist);
        ((AdjListNodePtr)(pdnHeadSublist->data.ptrvalue))->potential = potential/2.0;
	    	pdnmgHere = pdnmgHere->next;
			}
			pdnListPmmds = pdnListPmmds->next;
    }
    TotalPotential = TotalPotential / 2.0;
/*    pvnWorldHead =*/ FreeAllWorlds ();
}



/**
 *
 * PrintResNumbers (DValNodePtr pdnListPmmds):
 * Prints out residue numbers for a given list of PMMDs.
 *
 */

void PrintResNumbers (DValNodePtr pdnListPmmds) {

    PMMD pmmd;
    PDNMG pdnmgHere;

    printf ("ResNum");

    while (pdnListPmmds != NULL) {

	pmmd = (PMMD)(pdnListPmmds->data.ptrvalue);
	pdnmgHere = pmmd->pdnmgHead;

	while (pdnmgHere != NULL) {
	    printf ("\t%d", pdnmgHere->choice);
	    pdnmgHere = pdnmgHere->next;
	}

	pdnListPmmds = pdnListPmmds->next;
    }

    printf ("\n");

}



/**
 *
 * PrintSecStrucInfo (DValNodePtr pdnListPmmds):
 * Prints out secondary structure information for a given list of PMMDs.
 *
 */

void PrintSecStrucInfo (DValNodePtr pdnListPmmds) {

    PMMD pmmd;
    PDNMG pdnmgHere;
    PMGD pmgd;

    printf ("SecStru");

    while (pdnListPmmds != NULL) {

	pmmd = (PMMD)(pdnListPmmds->data.ptrvalue);
	pdnmgHere = pmmd->pdnmgHead;

	while (pdnmgHere != NULL) {

	    pmgd = (PMGD)(pdnmgHere->data.ptrvalue);

	    if (pmgd->bNCBISecStru & (Byte) SS_HELIX) {
		printf ("\tH");
	    } else if (pmgd->bNCBISecStru & (Byte) SS_STRAND) {
		printf ("\tE");
	    } else if (pmgd->bNCBISecStru & (Byte) SS_TURN) {
		printf ("\tT");
	    } else {
		printf ("\tC");
	    }

	    pdnmgHere = pdnmgHere->next;
	}

	pdnListPmmds = pdnListPmmds->next;
    }

    printf ("\n");

}



/**
 *
 * ReplaceTempWithPotential (DValNodePtr pdnAdjList, Int2 iModelNum):
 *
 * Takes an adjacency list and replaces the temperature information for
 * each alpha carbon residue with the potential score for that residue.
 *
 */

void ReplaceTempWithPotential (DValNodePtr pdnAdjList, Int2 iModelNum) {

    DValNodePtr pdnHere = pdnAdjList;
    DValNodePtr pdnSublist;

    PMGD pmgd;
    PMAD pmad;
    PALD pald;

    FloatHi ptnl;

    if (pdnAdjList == NULL) {
	ErrPostEx (SEV_INFO, 9, 1, "Empty Adjacency list.");
	return;
    }

    while (pdnHere != NULL) {

	pdnSublist = (DValNodePtr)(pdnHere->data.ptrvalue);
	pmgd = ((AdjListNodePtr)(pdnSublist->data.ptrvalue))->pmgd;
	ptnl = ((AdjListNodePtr)(pdnSublist->data.ptrvalue))->potential;

	if ((pmad = FindCAlpha (pmgd->pvnmaAHead)) != NULL) {
	    if (((pald = GetAtomLocs (pmad, iModelNum)) != NULL) &&
		(pald->iFloatNo >= 4)) {

		pald->pflvData[4] = (FloatLo) ((ptnl * -1.0) + 15.0);

	    }
	}

	pdnHere = pdnHere->next;

    }

}




/*************************************************************************
Kaca - Allocates one-dimensional array for contact map
is used with genome specific secondary structure potentials
NOTE: iSequenceLen is a sum of all lengths of protein chains
**************************************************************************/
void NewContactMap(BoolPtr *ppbContactMap, Int4 iSequenceLen)
{

	if (iSequenceLen == 0) {
		ErrPostEx (SEV_ERROR, 5, 3, "NewContactMap: Invalid param.");
		return;
	}
	*ppbContactMap = MemNew((size_t) (iSequenceLen * iSequenceLen * NUM_BINS * sizeof(Boolean)));
	return;

}

/****************************************************************************
Kaca - Prints the contacts map into a file; 6 bins of 0s and 1s
NOTE: iSequenceLen is a sum of all length of protein chains (get it using SumChainLen)
pcSequence is concatenated sequence for all chains
(in the same order as they appear in Biostruc - get it using ConcProtSequence)
*****************************************************************************/

void DumpContactMap(BoolPtr pbContactMap, Int4 iSequenceLen, FILE* pout, CharPtr pcSequence)
{
	Int2 binCounter = 0, rowCounter = 0, colCounter = 0;
	Int4 piDim[3];
	Int4 i = 0;

	if ((iSequenceLen == 0) || (pout == NULL) || (pcSequence == NULL)){
		ErrPostEx (SEV_ERROR, 5, 3, "DumpContactMap: Invalid param.");
		return;
	}

	piDim[0] = iSequenceLen;
	piDim[1] = iSequenceLen;
	piDim[2] = NUM_BINS;

	/* Iterate through bins */
	for (binCounter = 0; binCounter < NUM_BINS; binCounter++) {
		/* print out bin heading */
		fprintf(pout, "bin %d\n", (int) binCounter);
		/* print sequence column header */
		fprintf(pout, " ");
		for(i = 0; i<iSequenceLen; i++){
			fprintf(pout, "%c", pcSequence[i]);
		}
		fprintf(pout, "\n");
		/* iterate through rows */
		for (rowCounter = 0; rowCounter < iSequenceLen; rowCounter++) {
			/* print out row header aa */
			fprintf(pout, "%c", pcSequence[rowCounter]);
			/* iterate through columns */
			for (colCounter = 0; colCounter < iSequenceLen; colCounter++){
				if(rowCounter <= colCounter){
					if(pbContactMap[Get3DArrayValue(rowCounter, colCounter, binCounter, piDim)] == TRUE) fprintf(pout, "1");
					else fprintf(pout, ".");
				}
				else fprintf(pout, ".");
			}
			fprintf(pout, "\n"); /* end of column iteration */
		}
		fprintf(pout, "\n"); /* end of row iteration */
	} /* end of bin iteration */
	return;
}

/**********************************************************************************
Kaca - get sum of all protein residues in the structure
***********************************************************************************/
Int4 SumChainLen(PDNMM  pdnmm)
{
	Int4 itotalres = 0;

	if(pdnmm == NULL){
		ErrPostEx (SEV_ERROR, 5, 3, "SumChainLen: Invalid param.");
		return 0;
	}

	while (pdnmm != NULL){
		if(IsProtein((PFB)(pdnmm->data.ptrvalue)))
			itotalres = itotalres + ((PMMD)(pdnmm->data.ptrvalue))->iResCount;
		pdnmm = pdnmm->next;
	}
	return itotalres;

}


/**********************************************************************************
Kaca - get concatenated sequence of all protein chains withing biostruc in order of
appearance
ALLOCATES MEMORY!
***********************************************************************************/
CharPtr ConcProtSequence(PDNMM  pdnmm)
{
	Int4 itotalres = 0;
	CharPtr pseq = NULL;
	PDNMG pdnmg = NULL;
	PMGD pmgd = NULL;
	Int4 ipos = 0;
	Int2 protcount = 1;

	if(pdnmm == NULL){
		ErrPostEx (SEV_ERROR, 5, 3, "ConcProtSequence: Invalid param.");
		return NULL;
	}

	if((itotalres = SumChainLen(pdnmm)) == 0){
		ErrPostEx (SEV_ERROR, 5, 3, "ConcProtSequence: Zero sequence length.");
		return NULL;
	}
	pseq = MemNew((size_t) (itotalres+1) * sizeof(Char));

	/* iterate through molecules */
	while (pdnmm != NULL){
		if(IsProtein((PFB)(pdnmm->data.ptrvalue))){
			printf("protein %d\n", (int) protcount++);
			pdnmg = ((PMMD)(pdnmm->data.ptrvalue))->pdnmgHead;
			/* iterate through residues */
			while(pdnmg != NULL){
				pmgd = pdnmg->data.ptrvalue;
				pseq[ipos] = (pmgd->pcIUPAC)[0];
				ipos++;
				pdnmg = pdnmg->next;
			}
		}
		pdnmm = pdnmm->next;
	}
	pseq[ipos] = '\0';
	return pseq;

}



/**********************************************************************************
Kaca - given a valnode of protein chain length, molecule id and residue id
compute the position of the residue starting from chain 1 and residue 1
***********************************************************************************/

Int4 TranslateResNum(ValNodePtr pvnMolLen, Int4 iMolId, Int4 iResId)
{
	Int2 i = 0;
	Int4 iresnum = 0;

	if((pvnMolLen == NULL) || (iMolId == 0) || (iResId == 0)){
		ErrPostEx (SEV_ERROR, 5, 3, "TranslateResNum: Invalid param.");
		return 0;
	}
	if(ValNodeLen(pvnMolLen) < iMolId) {
		ErrPostEx (SEV_ERROR, 5, 3, "TranslateResNum: Invalid molecule ID.");
		return 0;
	}

	for (i = 1; i < iMolId; i++){
		iresnum = iresnum + pvnMolLen->data.intvalue;
		if(pvnMolLen != NULL) pvnMolLen = pvnMolLen->next;
	}
	iresnum = iresnum + iResId;
	return iresnum;
}


/**********************************************************************************
Kaca - Fills out N x N (N is length of the sequence, handles multiple chains)
contact map which is implemented as one-dimensional array of
Booleans (the function is based on CountStructureContacts)
the function can be used with POT_BL, POT_SS, POT_SB, POT_SB2
The Boolean map has dimensions of sequence_len x sequence_len x number_of_ bins but
is implemented as one-dimensional array
***********************************************************************************/
#define CONTACT_DIST 7
Int4 CreateContactMap(BoolPtr *ppbContactMap,
			DValNodePtr pdnListPmmds,
			Int2 iModelNum,
			Boolean bInclusiveWindow,
			Int2 iWindowSize,
			Boolean bModelling,
			POT_TYPE pot_type)
{
	Boolean bExclusiveWindow = (Boolean) !bInclusiveWindow;

	DValNodePtr pdnmgHeadPsBList = NULL;
	DValNodePtr pdnmgHere = NULL;
	DValNodePtr pdnPmmdsTemp = NULL;
	ValNodePtr pvnAtomList = NULL, pvnHere = NULL;
	ValNodePtr pvnMolLen = NULL;

	PWS pwsHead = NULL;
	PWS pwsThis = NULL;
	vec middle;
	Int2 iResId1, iResId2, iResIdTrue2, iMolId1, iMolId2, bin, ssbin, sbbin;
	Int4 iResNum1, iResNum2;

	PMMD pmmd1, pmmd2;
	PMGD pmgd1, pmgd2;
	PMAD pmad1, pmad2;
	PALD pald1, pald2;

	BoolPtr pbContactMap = *ppbContactMap;
	Int4 icount = 0, itotalcounts = 0;
	Int4 piDim[3];
	Int4 value = 0;
	Int4 itotalres = 0;

	/* figure out how many protein residues are in the structure and record chain lenghts in a valnode */
	pdnPmmdsTemp = pdnListPmmds;
	while (pdnPmmdsTemp != NULL){
		if(IsProtein((PFB)(pdnPmmdsTemp->data.ptrvalue))){
			ValNodeAddInt(&pvnMolLen, 0, ((PMMD)(pdnPmmdsTemp->data.ptrvalue))->iResCount);
			itotalres = itotalres + ((PMMD)(pdnPmmdsTemp->data.ptrvalue))->iResCount;
		}
		else ValNodeAddInt(&pvnMolLen, 0, 0);
		pdnPmmdsTemp = pdnPmmdsTemp->next;
	}


	piDim[0] = itotalres;
	piDim[1] = itotalres;
	piDim[2] = NUM_BINS;

	/* make pseudo-beta carbons */
	pdnListPepMGDs = NULL;
	if(CreatePsBList (&pdnmgHeadPsBList, pdnListPmmds, iModelNum, FALSE, bModelling, FALSE)) return 1;
	pdnmgHere = pdnmgHeadPsBList;

	/* Create a world containing all of the pseudo-beta carbons */
	while (pdnmgHere != NULL) {
		pmgd1 = (PMGD) pdnmgHere->data.ptrvalue;
			pwsThis = AddtoWorld (pwsHead, iModelNum, (PFB) pmgd1);
			if (pwsHead == NULL) {
				pwsHead = pwsThis;
			}
		pdnmgHere = pdnmgHere->next;
	}

	InstantiateWorld (1, pwsHead);
	pdnmgHere = pdnmgHeadPsBList;

	/* For each pseudo-beta carbon */
	while (pdnmgHere != NULL) {
	    itotalcounts++;
	    pmgd1 = (PMGD)(pdnmgHere->data.ptrvalue);
	    pmad1 = (PMAD)(pmgd1->pvnmaAHead->data.ptrvalue);
	    pald1 = GetAtomLocs (pmad1, iModelNum);

	    pmgd1 = (PMGD)(pmad1->pfbParent);
	    pmmd1 = (PMMD)(pmgd1->pfbParent);

	    iResId1 = pdnmgHere->choice;
	    iMolId1 = ((PDNMM)(pmmd1->pdnmmLink))->choice;

	    middle[0] = pald1->pflvData[0];
	    middle[1] = pald1->pflvData[1];
	    middle[2] = pald1->pflvData[2];

	    pvnAtomList = FindAtomsIn (pwsHead, middle, CONTACT_DIST);
	    pvnHere = pvnAtomList;
            /* For each atom within 10 angstroms sphere */
	    while (pvnHere != NULL) {
			itotalcounts++;
			pald2 = (PALD)(pvnHere->data.ptrvalue);
			pmad2 = (PMAD)(pald2->pfbParent);
			pmgd2 = (PMGD)(pmad2->pfbParent);
			pmmd2 = (PMMD)(pmgd2->pfbParent);

			iMolId2 = ((PDNMM)(pmmd2->pdnmmLink))->choice;
			iResId2 = pmad2->iIndex; /* what the heck is this number */
			iResIdTrue2 = ((PDNMG)(pmgd2->pdnmgLink))->choice;

			if (((bExclusiveWindow) && ((iMolId2 != iMolId1) || ((iMolId2 == iMolId1) && (abs(abs(iResId1) - abs(iResId2)) > iWindowSize)))) ||
			((bInclusiveWindow) && (iMolId1 == iMolId2) && (iResId1 != iResId2) && (abs(abs(iResId1) - abs(iResId2)) >= iWindowSize))) {

				if(iResId2 < 0 || iResIdTrue2 < 0) {
					ErrPostEx(SEV_ERROR,0,0,"Invalid Residue.");
				} else {
					/* get residue numbers counting from molecule1 and residue1 and substract one; it goes to an array*/
					iResNum1 = TranslateResNum(pvnMolLen, iMolId1, iResId1) - 1;
					iResNum2 = TranslateResNum(pvnMolLen, iMolId2, iResIdTrue2) - 1;
					if((bin = GetBinNumber (pald1, pald2)) != -1) {
						if (pot_type == POT_BL) {
							value = Get3DArrayValue(iResNum1,iResNum2,bin,piDim);
							pbContactMap[value] = TRUE;
							icount++;
						} else if (pot_type == POT_SS) {
							if((ssbin = GetSSBin(pmgd1,pmgd2)) != -1) {
								if (pot_type == POT_SS) {
									value = Get3DArrayValue(iResNum1,iResNum2,ssbin,piDim);
								}
								pbContactMap[value] = TRUE;
								icount++;
							}
						} else if (pot_type == POT_SB) {
							if((ssbin = GetSSBin(pmgd1,pmgd2)) != -1) {
								if((sbbin = GetSBBin(bin, ssbin)) != -1) {
									value = Get3DArrayValue(iResNum1,iResNum2,sbbin,piDim);
									pbContactMap[value] = TRUE;
									icount++;
								}
							}
						} else if (pot_type == POT_SB2){
							if((ssbin = GetSSBin(pmgd1,pmgd2)) != -1) {
								if((sbbin = GetSBBin2(bin, ssbin)) != -1) {
									value = Get3DArrayValue(iResNum1,iResNum2,sbbin,piDim);
									pbContactMap[value] = TRUE;
									icount++;
								}
							}
						}
					}
     				}
			}
			pvnHere = pvnHere->next;

	    }
	    ValNodeFree(pvnAtomList);
	    pdnmgHere = pdnmgHere->next;
	}
	ErrPostEx(SEV_INFO,0,0, "Possible Counts: %ld  Good Counts: %ld\n",(long) itotalcounts, (long)icount);
	FreePsBList (&pdnmgHeadPsBList);
	FreeAllWorlds ();

	return 0;
}




/*
$Log: potential.c,v $
Revision 1.39  2004/02/02 20:47:32  mjdumont
Minor changes

Revision 1.38  2003/09/14 02:38:23  feldman
Fixed unused variables and other minor compiler warnings

Revision 1.37  2002/10/26 16:34:32  michel
Updated CreatePsBList function call

Revision 1.36  2002/10/23 17:58:51  kaca
added contact map functionality


Revision 1.35  2002/09/25 19:56:28  feldman
Moved bryant minimization to foldtrajlib

Revision 1.34  2002/09/20 14:56:32  feldman
Added Bryant energy minimization

Revision 1.33  2002/08/13 19:46:45  kaca
added POT_SB2 Bryant/Secondary potential

Revision 1.32  2002/07/08 15:18:01  michel
fixed compiler warning

Revision 1.31  2002/07/02 20:43:15  kaca
added POT_SB potential

Revision 1.30  2002/06/18 19:18:11  michel
fixed very stupid bug

Revision 1.29  2002/06/18 16:07:54  michel
Abstracted structure contact counts

Revision 1.28  2002/06/14 22:12:26  feldman
Added Zhang bDetailedList option, and Zhang crease energy option

Revision 1.27  2002/03/26 16:54:27  kaca
added 4D potential functions, AddtoBLPotentialEx and GetSSBin

Revision 1.26  2002/02/19 15:42:10  feldman
Fixed hard-to-find bug in optimized version of potenital.c

Revision 1.25  2002/01/09 20:17:08  michel
Changed EGORDJKAdded name to EGORDoubleJacknife and added extra function parameter for detailed jk scores

Revision 1.24  2001/11/19 16:00:13  michel
minor change

Revision 1.23  2001/11/09 22:20:46  michel
Made GetIndex generic, added modelling feature to CreatePsBlist, fixed AddtoBryantPotential, other fixes

Revision 1.22  2001/10/26 20:24:25  michel
minor bug fix

Revision 1.21  2001/10/18 22:17:16  michel
Added functionality for BL table writing

Revision 1.20  2001/09/10 20:30:38  feldman
removed unused arrays

Revision 1.19  2001/09/10 03:40:48  michel
addressed compiler warnings

Revision 1.18  2001/08/31 18:48:28  feldman
Made LoadCBTable static

Revision 1.17  2001/08/30 14:45:22  michel
Changed file identifier to string for DumpBLPotential

Revision 1.16  2001/08/28 17:11:10  hogue
Added PsCB determination using CA (from randwalk.c)
Fixed bug in AddToBLPotential

Revision 1.15  2001/08/20 23:07:03  michel
compiler warning removal

Revision 1.14  2001/08/20 21:04:30  michel
fixed compiler warnings

Revision 1.13  2001/08/02 19:04:21  michel
bug fixes in BL potential functions

Revision 1.12  2001/07/26 19:19:51  feldman
Changed way fragments are stored on disk to be more platform independent,
and removed id parameter from AddtoWorld

Revision 1.11  2001/07/26 19:11:26  michel
Added new BL potential functions

Revision 1.10  2001/06/27 01:32:25  michel
Added new functions to generate Bryant-Lawrence potential

Revision 1.9  2001/04/04 21:26:01  feldman
-Removed my BSPrintf, now using the one in SLRIlib - thus needed to add
	this library to makefiles
-changed all fgets to FileGets

Revision 1.8  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.7  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.6  2001/02/07 22:31:42  feldman
fixed remaining code for unused variables and unreachable code, etc.

Revision 1.5  2001/02/07 18:46:45  feldman
Removed a few unused variables

Revision 1.4  2000/08/14 20:24:18  feldman
Added LIBCALLBACK for function pointers where needed

Revision 1.3  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.2  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

