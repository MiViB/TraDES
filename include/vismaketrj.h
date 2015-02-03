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

#ifndef VISMAKETRJ_H
#define VISMAKETRJ_H

#include <foldtrajlib.h>
#include <slri_mmdblocl.h>
#include <clustlib_pub.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BLOSUM_FNAME "BLOSUM62"

TrajErr BuildTrjGraph(PMMD pmmd, Int2 model, Int2 start, Int2 end, Int4 mag, FloatLo sigma_x, FloatLo sigma_y,CharPtr fnam,Int2 percent, FloatLo tmp, FloatLo delta_tmp, DValNodePtr pdnIncList, Int2 SaveChis, CharPtr ssmask);
void MakeTG(CharPtr kbsequence, pMakeTrjParamBlock paramblock, Char*fnam);

#ifdef __cplusplus
}
#endif

#endif

/*
$Log: vismaketrj.h,v $
Revision 1.9  2004/09/24 19:09:36  hfeldman
Added import fragment option to maketrj

Revision 1.8  2003/04/04 21:54:04  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.7  2001/09/10 20:29:21  feldman
Added Uturn option to maketrj

Revision 1.6  2001/06/26 20:19:10  feldman
Moved stuff to clustal library

Revision 1.5  2001/03/15 18:54:52  feldman
removed mmdblocl from all makefiles which dont need it and added
my own mmdblocl.h to avoid using NCBI's which is different

Revision 1.4  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.3  2001/03/13 15:08:17  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.2  2001/02/15 20:27:57  feldman
Reworked maektrj to take a parameter block instead of a bunch
of individual parameters.

Revision 1.1  2001/02/01 21:53:37  feldman
added header file


*/

