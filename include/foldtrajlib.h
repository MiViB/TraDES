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


#ifndef FOLDTRAJLIB_H
#define FOLDTRAJLIB_H

/* includes */
#ifdef OS_MSWIN
#include <D4ALL.H>
#endif

#ifndef OS_MSWIN
#include <d4all.h>
#endif

#include <mmdbtraj_pub.h>
#include <objtraj.h>


#ifdef __cplusplus
extern "C" {
#endif

/* database averages and standard deviations for trans- */
#define BL_CACA 3.806		
/* and cis- residues */
#define BL_CACA_CISP 2.955	
#define BLSD_CACA 0.035
#define BLSD_CACA_CISP 0.097

/* for easy grepping of fold log */
#define LOG_TERM_CHAR ';' 

/* global variables */
/* parameter file variables */
extern paramblock prmblock;
/* quiet operation on if non-zero */
extern Byte traj_quiet;
extern Boolean firstiter;
extern Char tmpdbasename[PATH_MAX];
extern Char tmpskelfname[PATH_MAX];
/*extern Char loglock[PATH_MAX];*/
extern Char LogOutName[PATH_MAX];
extern Char ExtraOutName[PATH_MAX];
extern ByteStorePtr bspTempLog;
/*extern ByteStorePtr bspExtraLog;*/
extern FILE *fExtraLog;
extern Char CFG_local_cddname[PATH_MAX];
extern Char CFG_local_cddpath[PATH_MAX];

/* functions needed for web server stuff */
Int4 GetAlarmSecs(void);
CharPtr GetCGIpath(void);

/* various functions dealing with proteins at the atomic co-ordinate level */
/* set atom co-ordinates to vec */
TrajErr SetCoOrds PROTO((PALD paldHere,vec vHere));
/* makes a new PALD */
PALD AllocNewLoc PROTO((PMAD pmadThis,Int2 modelnum));
/* checks if atom has been placed */
Int2 IsAtomPlaced(PMAD pmadThis,Int2 Model);
/* remove pending H-bond status from atom */
void ClearHBondStatus(PMAD pmad,PHBS *pphbsHBondsToCheck);
/* conform time out */
void PrintFileTimedOut(void);

PEAS GetExtAAInfoEnc(Int2 idx);
Int2 AlreadyMadeTG(Int2 cnt);
void FreeResUsedTable(void);

/* distance constraint functions */
TrajErr CheckEndingDistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol);
TrajErr CheckEndCADistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol);
TrajErr CheckEndCDistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol);
TrajErr CheckDistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol);
TrajErr CheckAtomDistConstraints(CharPtr AtomName,Int2 res,vec vHere,Int2 Model,Int4 *viol);

#ifdef __cplusplus
}
#endif

#endif /* FOLDTRAJLIB_H */

/*  
$Log: foldtrajlib.h,v $
Revision 1.16  2003/01/24 16:44:37  feldman
Moved some defines around, made CHARMM callable as a thread

Revision 1.15  2001/10/10 19:53:54  feldman
Make H-bonds a bit more flexible (120 degrees rather than 150) and fixed up coding
for H-bonds in randwalk.c

Revision 1.14  2001/04/06 15:26:01  feldman
moved flipping functions to SLRI lib

Revision 1.13  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.12  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.11  2001/03/07 21:49:47  feldman
Made many important fixes to VisTraj such as bounds checking of
user-entered information, and tidied up interface further.
Fixed a few minor bugs as well, and tested RPSTraj feature of
VisTraj.

Revision 1.10  2001/02/09 20:18:33  feldman
Added fragments to PTGS

Revision 1.9  2001/02/06 18:37:39  feldman
-Moved around a few function headers
-Split up conform header into private and public sections

Revision 1.8  2001/01/26 20:37:38  feldman
Added getconstraints function

Revision 1.7  2001/01/23 18:02:09  feldman
Added headers for distance constraints functions

Revision 1.6  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.5  2001/01/12 20:01:50  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.4  2000/09/15 20:36:54  feldman
Added angle and dihedral to ASN trajgraph dist. constraints
Now counting distance constraint violations

Revision 1.3  2000/09/06 15:33:26  feldman
removed globals in MOBI from ensemb.c in library, now passed as
parameters instead

Revision 1.2  2000/08/14 20:32:44  feldman
Added some amber functions headers, and changed <> to "" for some
includes

Revision 1.1  2000/07/13 20:31:52  feldman
Tidied up header further

*/
