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

#ifndef POTENTIAL_H
#define POTENTIAL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Mark's hydrophobic to hydrophilic arrangement of amino acids */
#define MK_AA_3LETTER "AlaValIleLeuProPheMetTrpGlySerThrCysTyrGlnAsnAspGluHisLysArgPep"
#define MK_AA_1LETTER "AVILPFMWGSTCYQNDEHKR-"
/* Bryant Lawrence style -> alphabetical arrangement of amino acids */
#define BL_AA_3LETTER "AlaArgAsnAspCysGlnGluGlyHisIleLeuLysMetPheProSerThrTrpTyrValPep"
#define BL_AA_1LETTER "ARNDCQEGHILKMFPSTWYV-"


/* for parsing text files */
#define MAX_STRING_LENGTH 256
/* scale for potential table */
#define POTENTIAL_TABLE_SCALE_FACTOR 1000000.0
/* average distance from Ca to Cb over all residues */
#define PSB_SCALE_FACTOR 2.4
#define BRYANT_SCALE_FACTOR 10.0 /* relative weight of Bryant energy in energy function */

/* atomic numbers */
#define ATOMIC_NO_H 1
#define ATOMIC_NO_C 6
#define ATOMIC_NO_N 7

/* distance bins for energy function */ 
#define NUM_BINS 6
#define NUM_BINS4D 6

/* zhang atoms */
#define NUM_ZHANG 18
#define NUM_AA 21

typedef enum {MK_ALA, MK_VAL, MK_ILE, MK_LEU, MK_PRO,
              MK_PHE, MK_MET, MK_TRP, MK_GLY, MK_SER,
              MK_THR, MK_CYS, MK_TYR, MK_GLN, MK_ASN,
              MK_ASP, MK_GLU, MK_HIS, MK_LYS, MK_ARG, MK_PEP} MK_AA;
typedef enum {BIN_0_5, BIN_5_6, BIN_6_7, BIN_7_8, BIN_8_9, BIN_9_10} MK_BIN;

/* secondary structure contact potential defines */
#define NUM_BINS_SS 6

#define HELIX_HELIX 0
#define HELIX_STRAND 1
#define HELIX_COIL 2
#define STRAND_STRAND 3
#define STRAND_COIL 4
#define COIL_COIL 5

#define HELIX_NEAR 0
#define HELIX_FAR 1
#define HELIX_STRAND_ALL 2
#define STRAND_NEAR 3
#define STRAND_FAR 4
#define COIL_ALL 5


/*#define HELIX_HELIX 0
#define HELIX_STRAND 1
#define COIL_COIL 5*/
#define HELIX_STRAND_COIL 2
/*#define STRAND_NEAR 3
#define STRAND_FAR 4*/

Int2 CreatePsBList (DValNodePtr *ppdnmgResult, DValNodePtr pdnListPmmds, Int2 iModelNum, Boolean bUsingPep, Boolean bModelling, Boolean usecaonly);
Int2 GetBinNumber (PALD pald1, PALD pald2);
FloatHi GetBryantPotential (Int4Ptr piBryantPotential, Int2 res1, Int2 res2, Int2 bin);
void FreePsBList (DValNodePtr *ppdnmgHead);

#ifdef __cplusplus
}
#endif

#endif

/*  
$Log: potential.h,v $
Revision 1.10  2002/09/25 19:56:29  feldman
Moved bryant minimization to foldtrajlib

Revision 1.9  2002/08/13 19:46:45  kaca
added POT_SB2 Bryant/Secondary potential

Revision 1.8  2002/07/02 20:43:51  kaca
added POT_SB potential defines

Revision 1.7  2002/06/18 16:07:54  michel
Abstracted structure contact counts

Revision 1.6  2002/03/26 16:55:17  kaca
added definitions for sec. str. bins

Revision 1.5  2001/11/09 22:20:46  michel
Made GetIndex generic, added modelling feature to CreatePsBlist, fixed AddtoBryantPotential, other fixes

Revision 1.4  2001/10/26 20:24:25  michel
minor bug fix

Revision 1.3  2001/10/18 22:17:16  michel
Added functionality for BL table writing

Revision 1.2  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.1  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies


*/
