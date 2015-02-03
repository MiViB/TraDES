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

#ifndef TRADES_H
#define TRADES_H

/* includes */
#include <hfprogs.h>
#ifdef OS_UNIX
#include <signal.h>
/* not used anymore.. was used for web-based version only */
/*#define CPUTIME_MAX 6000*/
#endif
#include <asn.h>
#include <objloc.h>
#include <slri_misc.h>


#ifdef __cplusplus
extern "C" {
#endif

/* name of .rc file */
#ifdef OS_MSWIN
#define TASK_CFGFILE ".\\trades"
#else
#define TASK_CFGFILE "trades"
#endif

#define ZHANG_WINDOWSIZE 1   /* Reset to 1 CWVH 2012 - Small Peptide scores are off. Entropy Calculation Calibration */ 
#define BRYANT_WINDOWSIZE 3
#define VSCORE_WINDOWSIZE 1

/* write out log on interrupt */
Int2 Method;  /*CWVH added to replace Mobi Value */

#ifdef __cplusplus
}
#endif

#endif

/*
$Log: structrj.h,v $

Revision 2.0  2012/05/15  cwvhogue
Modified Scoring Functions for downstream IDP energy/entropy calculations, added solvation layer code, principal components
for elliposoidal parameters, and output is now normalized on x,y,z axes, extended log files.

Revision 1.13  2004/07/16 21:20:42  egarderm
Changed Bryant printing info to vscore, added window size define for vscore calls, replaced
Bryant energy computation with vscore

Revision 1.12  2003/07/14 20:01:00  egarderm
Changed TASK_CFGFILE to .\foldtraj on Windows to allow local INI file

Revision 1.11  2001/04/06 14:29:00  feldman
Completed fold-at-home Client-Server v1.0 now complete and ready for testing

Revision 1.10  2001/04/04 21:23:37  feldman
Completed file uploading segment of foldtrajlite

Revision 1.9  2001/03/30 22:21:39  feldman
Added foldtraj screensaver program

Revision 1.8  2001/03/29 02:52:24  feldman
fixed minor compiler warnings

Revision 1.7  2001/03/14 16:25:53  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.6  2000/10/16 15:52:41  feldman
moved ifdef OS_UNIX after hfprogs.h as needed

Revision 1.5  2000/08/23 21:19:39  feldman
Removed all remaining system calls to DSSP, now all internal
calls - thus you don't need dssp on your system for that part
to work

Also rearranged format of fold*.log slightly due to this - note
extended residues in DSSP are taken as 'B' or 'E' and helices
are 'H' only (not 'G' or 'I') - this is a slightly different
and more correct measure of sheet content than before (helix is same)

Revision 1.4  2000/08/14 20:42:01  feldman
change order of includes for Windows

Revision 1.3  2000/08/14 19:24:50  feldman
Corrected several Windows compatibility issues (all minor)
Improved aligntraj for homology modelling - now works correctly and only adds
  noise when absolutely necessary
Added AMBER potential functions

Revision 1.2  2000/07/14 20:09:24  feldman
Added parameter for MIMEBiostrucAsnGet for Bioseq

Revision 1.1  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies


*/
