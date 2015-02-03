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

#ifndef SLRIACCSURF_H
#define SLRIACCSURF_H

#ifdef __cplusplus
extern "C" {
#endif

/* CalcDSSPAccSurf
 * Calculates surface accessibility
 *
 *
 */
Int4 *CalcDSSPAccSurf(PMMD pmmdHere,Int2 order,Int2 modelnum);
Int4 *CalcMSAAccSurf(PMMD pmmdHere,Int2 order,Int2 modelnum);

/* CalcDSSPAssign
   returns an array of characters representing the basic DSSP
   secondary structure assignment 
 */
void CalcDSSPAssign(PMMD pmmdHere, Int2 modelnum,CharPtr pcSec);
void CalcDSSPAssignEx(PMMD pmmdHere, Int2 modelnum,CharPtr pcSec,Boolean AddFeatures);
CharPtr ObjectiveUTurns(PMMD pmmdHere, Int2 iModelNum);
double ComputeResidueASP(int resnum,int order);
FloatHi *CalcDSSPASP(PMMD pmmdHere,Int2 order,Int2 modelnum);

#ifdef __cplusplus
}
#endif

#endif

/*

$Log: slriaccsurf.h,v $
Revision 1.10  2003/03/07 19:59:22  feldman
Added AddFeatures option to CalcDSSPAssign

Revision 1.9  2002/06/12 22:14:54  feldman
Added MSA functions

Revision 1.8  2001/11/08 16:55:41  feldman
updated function headers

Revision 1.7  2001/10/10 19:53:16  feldman
Added GetASP function

Revision 1.6  2001/10/02 21:50:43  feldman
Added double-precision vector functions and used with Energy computations

Revision 1.5  2001/10/01 14:28:37  feldman
Added ASP functions

Revision 1.4  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

*/
