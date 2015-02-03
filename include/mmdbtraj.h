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


#ifndef MMDBTRAJ_H
#define MMDBTRAJ_H

/* includes */
#include <mmdbtraj_pub.h>

#ifdef __cplusplus
extern "C" {
#endif

/* frequency cutoff for U-turn calling (see r14.c) */
#define DEFAULT_CUTOFF 0.20

/* floating and fixed point 3-D vector definitions respectively */
typedef Int4 ivec[3];

/* SCWRL Rotamer Library record type */
typedef struct nlm_rotlibstruct {
  Char res[4];
  Int4 phi;   
  Int4 psi;   
  Int4 N;
  Int4 r1; 
  Int4 r2;
  Int4 r3;
  Int4 r4;
  FloatLo p;
  FloatLo chi1;
  FloatLo chi2;
  FloatLo chi3;
  FloatLo chi4;
} rotlibrecord, *protlibrecord;

/* corresponding vector functions for fixed-point vectors */
/* warning - not fully tested!!! */
Int4 iDot PROTO((ivec a,ivec b));
Int4 igetMag PROTO((ivec a));
void iNormalize PROTO((ivec res,ivec a));
void iVecAdd PROTO((ivec res,ivec a,ivec b));
void iVecSub PROTO((ivec res,ivec a,ivec b));
void iVecScale PROTO((ivec res, ivec a,Int4 scale));
void iNegateVec PROTO((ivec res, ivec a));
void iPrintVec PROTO((ivec a));
void iTranslate PROTO((ivec a,ivec trans));
void iRotatecos PROTO((Int2 axis, Int4 sina, Int4 cosa, ivec a,ivec res));
/* with A, B and C having known co-ordinates, this finds the unique D such that CD has the
   given magnitude, angle BCD=Beta and Chi is the dihedral angle formed by ABCD */
void iCalcNextCoOrd PROTO((vec vA,vec vB,vec vC,vec vD,FloatLo MagCD,FloatLo Beta,FloatLo Chi));

/* molecular alignment functions */
/* center molecule about origin */
Int2 CenterMol(PMMD pmmdHere,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model);
/* rotate molecule about a given axis */
Int2 RotateMol(PMMD pmmdHere,FloatLo theta,Int2 Axis,Int2 Model);
/* single value decomposition (SVD) function, "borrowed" from MOLMOL source */
Int2 SVD(FloatLo **a,Int2 m,Int2 n,FloatLo *w,FloatLo **v);

/* Zhang/Bryant potential related functions */
Int2 GetIndex (CharPtr str);
void BackupPepMGDs(void);
void RestorePepMGDs(void);

#ifdef __cplusplus
}
#endif

#endif /* MMDBTRAJ_H */

/*  
$Log: mmdbtraj.h,v $
Revision 1.18  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.17  2000/08/18 19:00:17  adrian
changed r14seqfilecalc to r14bscalc
added DSSP function definitions to slriaccsurf.h

Revision 1.16  2000/07/13 20:31:52  feldman
Tidied up header further

Revision 1.14  2000/07/07 21:30:26  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.12  2000/06/23 19:47:04  feldman
Went back to using bzip2 for compression and no
longer use a bytestore but dump directly to a file

- added surface accessibility output to fold logs

Revision 1.11  2000/06/22 19:01:01  adrian
fixed headers for r14

Revision 1.10  2000/06/22 16:33:08  feldman
minor bugfix and added VDW radii to hfprogs

Revision 1.9  2000/06/22 13:05:01  feldman
Fixed minor bugs with log file nameing,
fixed error in extended residue calculation and DSSP usage
added alterresidue function

Revision 1.8  2000/06/21 19:01:59  john
Removed conflict indicators

Revision 1.7  2000/06/21 18:59:44  john
Finally fixed rotation
Added Proline & Cystine Procedures
Added FillTG functionality under Add Noise

Revision 1.6  2000/06/20 18:59:11  feldman
Removed dependence of r14.h, no longer needed, merged to mmdbtraj

Revision 1.5  2000/06/19 19:27:11  adrian
reverted to 1.3 code

Revision 1.3  2000/06/15 17:11:02  feldman
Added Replace option when TrajGraphWrite-ing

Revision 1.2  2000/06/15 15:52:10  feldman
Corrected Bzip2 log dumping and improved H-bond treatment (fixed
bug that was crashing program)

Revision 1.1.1.1  2000/06/09 18:12:47  feldman
TraDES include files

*/
