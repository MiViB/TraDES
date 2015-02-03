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

#ifndef ROTLIB_H
#define ROTLIB_H

#include <slri_misc.h>

#ifdef __cplusplus
extern "C" {
#endif

/* SCWRL library filenames and sizes */
#define SCWRL_BZ_SIZE 3332417L
#define SCWRL_SIZE 22780160L
#define SCWRL_FNAM "scwrlbin34.lib.bz2"
#define NEWSCWRL_BZ_SIZE 2803905L
#define NEWSCWRL_SIZE 12266240L
#define NEWSCWRL_FNAM "rotlib.bin.bz2"
#define BBDEP_LIBNAME "bbdep00.Nov.r1lib"
#define BBINDEP_LIBNAME "bbind00.Nov.lib"

#ifdef __cplusplus
}
#endif

#endif

/*
$Log: rotlib.h,v $
Revision 1.4  2001/07/24 17:27:13  feldman
Updated rotamer library version and fixed bug - order of building
for Ramachandran space had prevented rotamer library from being
used properly

Revision 1.3  2001/04/06 15:26:01  feldman
moved flipping functions to SLRI lib

Revision 1.2  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.1  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies


*/

