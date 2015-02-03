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


#ifndef SLRI_MMDBLOCL_H
#define SLRI_MMDBLOCL_H

#include <mmdbapi1.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DocUid
#define DocUid Int4
#endif

Boolean LIBCALL IterateMMDBId(void);
void LIBCALL ResetMMDBIterate(void);
Int4 LIBCALL GetMMDBIdCurrent(void);
Boolean LIBCALL MMDBInit (void);
void LIBCALL MMDBFini (void);
BiostrucPtr LIBCALL MMDBBiostrucGet (DocUid uid, Int4 mdlLvl, Int4 maxModels);
DocUid  LIBCALL MMDBEvalPDB (CharPtr str);

#ifdef __cplusplus
}
#endif

#endif /* SLRI_MMDBLOCL_H */

/*  
$Log: slri_mmdblocl.h,v $
Revision 1.1  2001/03/15 18:54:51  feldman
removed mmdblocl from all makefiles which dont need it and added
my own mmdblocl.h to avoid using NCBI's which is different


*/
