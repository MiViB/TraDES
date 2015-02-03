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


#ifndef HFPROGS_H
#define HFPROGS_H

/* include file for all executables */
#include <mmdbtraj_pub.h>
#include <ncbi.h>
#include "version.h"

#ifdef __cplusplus
extern "C" {
#endif

/* for easy grepping of fold log */
#define LOG_TERM_CHAR ';' 

/* The name of the HomTraj config file */
#define HOMTRAJ_CONF "HomTraj.conf"
#define HANDLE_LENGTH 8


/* Amino Acid Molecular Weight in AMU */
FloatLo MWaalist[] = {71.03712F,156.10112F,114.04293F,115.02695F,103.00919F,128.05858F,129.04260F,57.02147F,137.05891F,113.08407F,113.08407F,128.09497F,131.04049F,147.06842F,97.05277F,87.03203F,101.04768F,186.07932F,163.06333F,99.06842F};
/* All 20 amino acids */
CharPtr aalist = "ARNDCQEGHILKMFPSTWYV";

/* global parameters used by all foldtraj related executables */
paramblock prmblock;
/* quiet operation on if non-zero */
Byte traj_quiet;
Boolean firstiter;
Char tmpdbasename[PATH_MAX]="";
Char tmpskelfname[PATH_MAX]="";
/*Char loglock[PATH_MAX]="";*/
Char LogOutName[PATH_MAX]="";
Char ExtraOutName[PATH_MAX]="";
ByteStorePtr bspTempLog=NULL;
/*ByteStorePtr bspExtraLog=NULL;*/
FILE *fExtraLog=NULL;
Char CFG_local_cddname[PATH_MAX]="";
Char CFG_local_cddpath[PATH_MAX]="";
Char CFG_local_datafilepath[PATH_MAX]="";
/* for Vibrant monitors */
Int4 volatile ProgramProgress=0;
Int4 volatile ProgramProgressMax=0;

/* extended AA naming related functions */
TrajErr AdjustResidueNumbers(int numAA,int res_num,int Adjust_Type);
CharPtr AlterResidueInSequence(CharPtr seq,int resnum,CharPtr newaa,Int2 dowhat);

/* returns an array of containing predicted U turns
 * given the sequence 
 */
CharPtr CalcR14Assign(ByteStorePtr bsSequence);

/* r14 turn prediction-from-sequence functions */
Int2 r14calc(ByteStorePtr seqbsp, ValNodePtr *headvnpp);

/* Trajectory graph filter generating functions */
void MakeLPFilter(TrajFilter filt,FloatLo magnitude,FloatLo sd);
void MakeGaussFilter(TrajFilter filt,FloatLo magnitude,FloatLo sd);
void MakeSmoothFilter(TrajFilter filt);
/* trajectory graph operators */
void TrajScale(PTGS ptgsHere,FloatLo ScaleFactor);
TrajErr TrajAdd(PTGS ptgsDest,PTGS ptgsAdd);
void TrajCopy(Int4 *src,Int4 *dest,Int4 graphwidth);

TrajErr LoadCBTable(void);

/* main Maketrj function */
VoidPtr MakeTrj(VoidPtr ptr);

/* interrupt handler */
void MyInterruptHandler(int sig);

#ifdef __cplusplus
}
#endif

#endif /* HFPROGS_H */

/*  
$Log: hfprogs.h,v $
Revision 1.37  2003/11/07 19:25:02  ksnyder
Added HomTraj define

Revision 1.36  2003/03/07 23:33:13  feldman
MAde programprogress an Int4 and added wwwget paramblock

Revision 1.35  2002/09/26 13:23:22  michel
Moved BuildMIMEBiostruc to mmdbtrajlib

Revision 1.34  2002/06/13 13:09:53  michel
moved CalcExtendedResidues to mmdbtraj_pub.h

Revision 1.33  2001/08/31 18:31:14  feldman
Added U-turn ASN.1 header and U-turn database related function prototypes

Revision 1.32  2001/05/25 21:50:00  feldman
Added datafile path global variable

Revision 1.31  2001/04/06 14:29:00  feldman
Completed fold-at-home Client-Server v1.0 now complete and ready for testing

Revision 1.30  2001/03/29 20:56:58  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.29  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.28  2001/03/23 18:44:17  feldman
Integrated foldtraj into vistraj (!)
and cleaned up a few minor bugs

Revision 1.27  2001/03/15 15:59:56  feldman
revamped makefiles so they are a lot smaller now, all common stuff
is in .mk files - all have been tested and most still compile (Except
one or two with known problems)

Revision 1.26  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.25  2001/02/23 18:03:54  feldman
Changed MakeTrj to void*, void*

Revision 1.24  2001/02/15 20:29:03  feldman
Changed maketrj to take a parameter block

Revision 1.23  2001/02/06 18:37:39  feldman
-Moved around a few function headers
-Split up conform header into private and public sections

Revision 1.22  2001/01/23 18:02:09  feldman
Added headers for distance constraints functions

Revision 1.21  2001/01/12 20:01:50  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.20  2000/10/25 15:15:13  feldman
Made further updates, multiple model support is correct now
and relocated to a single function for loading correct model
and extensive error handling

Revision 1.19  2000/10/24 20:57:24  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.18  2000/10/06 18:14:41  feldman
Change Sarray so element value indicates residue it bonds with

Revision 1.17  2000/08/31 14:53:32  feldman
Changed filter functions to take filt, not *filt

Revision 1.16  2000/08/18 19:00:17  adrian
changed r14seqfilecalc to r14bscalc
added DSSP function definitions to slriaccsurf.h

Revision 1.15  2000/08/14 20:32:44  feldman
Added some amber functions headers, and changed <> to "" for some
includes

Revision 1.14  2000/08/09 19:17:07  feldman
-minor update/bugfixes added
-version.h contains default version number

Revision 1.13  2000/07/13 20:51:32  adrian
added cutoff parameter to r14seqfilecalc

Revision 1.12  2000/07/13 20:31:52  feldman
Tidied up header further

Revision 1.11  2000/07/07 21:30:26  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.9  2000/06/23 19:47:04  feldman
Went back to using bzip2 for compression and no
longer use a bytestore but dump directly to a file

- added surface accessibility output to fold logs

Revision 1.8  2000/06/22 19:01:01  adrian
fixed headers for r14

Revision 1.7  2000/06/22 16:33:08  feldman
minor bugfix and added VDW radii to hfprogs

Revision 1.6  2000/06/22 13:05:01  feldman
Fixed minor bugs with log file nameing,
fixed error in extended residue calculation and DSSP usage
added alterresidue function

Revision 1.5  2000/06/21 18:59:44  john
Finally fixed rotation
Added Proline & Cystine Procedures
Added FillTG functionality under Add Noise

Revision 1.4  2000/06/20 18:59:11  feldman
Removed dependence of r14.h, no longer needed, merged to mmdbtraj

Revision 1.3  2000/06/16 13:15:45  john
added rotatecos function

Revision 1.2  2000/06/15 17:11:02  feldman
Added Replace option when TrajGraphWrite-ing

Revision 1.1.1.1  2000/06/09 18:12:47  feldman
TraDES include files

*/
