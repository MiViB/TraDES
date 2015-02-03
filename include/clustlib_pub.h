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
 
#ifndef _CLUSTLIB_
#define _CLUSTLIB_

#include <ncbi.h>
#include <mmdbapi.h>  /* Structure manipulation */


#ifdef __cplusplus
extern "C" { /* } */
#endif


/*******************************************************************
* include <seqhound.h> in your main program for SeqHound Functions *
* use -z muldefs in linking application to use first defined symbol*
* for MMDBBiostrucGet                                              *
* include <accentr.h> in your main program for Entrez Functions    *
*******************************************************************/

typedef struct csanode {
	struct csanode PNTR next;
	CharPtr pcSeqName;
	Int4    iGi;
	Int4    iMMDB;
	ValNodePtr pvnSeqAlnPart;
	CharPtr pcSeqAln;
	Int4    iLen;
	PDNMS   pdnmsStructure;
	Boolean bStrucMask; /*  TRUE if this is just the sec-structure mask */
	Boolean bPosGapPen; /* TRUE if this is postition specific gap penalties */
    Pointer pExtra;     /* hanger for whatever... - used for hanging secondary structure info on sequencetodomainstructurealignment */
} CSAN, *PCSAN;

PCSAN LIBCALL NewCSAN(void);
void FreeCSAN(PCSAN pcsanThis);
PCSAN LIBCALL ReadCSAN(CharPtr fnam);
Int4 LIBCALL WriteCSAN(PCSAN pcsanThis, FILE *pFile);

/* SeqHound can be used to MakeStrucMask, simply link library */
PCSAN CSANStructMask(Int4 iMMDBid, Int4 iGiStr);
Int4 FillCSANWithStru(PCSAN pcsanThis, PMMD pmmdThis, Int4 iLen);
Int4 FillCSANWithMask(PCSAN pcsanThis, PMMD pmmdThis, Int4 iLen);

Int2 SSMasktoDisk(Int4 iMMDBid, Int4 iGiStr);

PCSAN NetValidateCSAN(PCSAN pcsanThis, Boolean bCheckSeq);


/* iType definitions */
#define CLUSTSEQ 1
#define CLUSTSTR 2
#define CLUSTSS  3

/* to write out the alignments in a file */
#define S2S_HEAD "alignments"
#define S2S_TAIL ".txt"

/* Generate a sequence to structure alignment with identity and domain occupancy thresholds */
Int4 MakeSequenceToDomainStructureAlignment
		(PCSAN pcsanAlign, Boolean bSecStruc, 
		FloatLo fDom_Id, FloatLo fDom_Occ,
		Int4 PNTR pidomcount, Int4 PNTR pidomthresh,
		ValNodePtr PNTR pvndomids, ValNodePtr PNTR pvndomocc, ValNodePtr PNTR pvndomlen, ValNodePtr PNTR pvndomok);

/* Get sequence alignment info with an already modeled sequence in the bReserved field */
Int2 CLUSTAL_GetSequenceAlignmentInfo(PDNMG pdnmgHead, ByteStorePtr bsp, FloatLo fDomIDThreshold, FloatLo fDomOccThreshold, 
			ValNodePtr *pvndomids, ValNodePtr *pvndomocc, ValNodePtr *pvndomlen, ValNodePtr *pvndomok, Int4 *piDomCount, Int4 *piDomOverThreshold);


/* Use this function to get either the sequence or the secondary structure information from the structure */
/* If bUseRealSeq is set to FALSE, then it looks in the bReserved position of the residue graph */
/* MakeSequenceToDomainStructureAlignment calls GetSeqInfoFromStruc and places the sequences in the pvnSeqAlnPart */
/* of the respective sequence and structure csanode.  The Secondary Structure sequence is placed in the pExtra part */

#define CLUSTLIB_INFO_SEQ     0 /* Get the sequence annotation */
#define CLUSTLIB_INFO_VASTSS  1 /* Get the VAST annotation for secondary structure */
#define CLUSTLIB_INFO_PDBSS   2 /* Get the PDB annotation for secondary structure */
#define CLUSTLIB_INFO_DOM     3 /* Get a comma delimited list of domain numbers for each residue */
void GetSeqInfoFromStruc(PMMD pmmdThis, Uint2 CLUSTLIB_INFO_TYPE, Boolean bUseRealSeq, CharPtr PNTR ppcSeq);


/* This writes a file named iSeqGI_iStrucGI.[seq|str|ss] and writes >pcDefline followed by pcSeq.*/
/* It also forces a column size of 60 characters */
void DumpPKBFASTA(Int4 iSeqGI, Int4 iStrucGI, CharPtr pcDefline, CharPtr pcSeq, Int2 iType);

/* This creates a database by appending successive sequences, removes X's from sequences and writes
to file named database<iFileID.[seq|str|ss] with a 50 character column... it also adds ! for the header 
and @ for the footer for each sequence */
void DumpGORFASTA(Int4 iSeqGI, Int4 iStrucGI, Int4 iFileID, CharPtr pcSeq, Int2 iType, Boolean bFirst);

/* if you want to supply your own file pointer */
void DumpPKBFASTAFile(CharPtr pcDefline, CharPtr pcSeq, FILE *fp);
void DumpGORFASTAFile(CharPtr pcDefline, CharPtr pcSeq, FILE *fp);

/* the generic function - fields cannot be null */
void DumpSequenceToFile(CharPtr pcHeader, CharPtr pcDefline, CharPtr pcFooter, CharPtr pcSeq, FILE *fp, Int4 columnsize);





/* Entrez is required for these functions */
#ifdef _ENTREZAPI_
Int4 FillCSANWithSeq(PCSAN pcsanThis, BioseqPtr pbsq, Int4 iLen);
void WriteFASTAfromGI(Int4 iGiSeq, FILE *pOut);
#endif 

PDNMG GetResidueListFromPCSAN(PCSAN pcsan);
PMMD MMDBGetMoleculeByGI(PDNMS pdnmsThis, Int4 iGi);

Int2 LoadBLOSUM(CharPtr fnam);
Int2 GetBLOSUM(Char x, Char y);
void FreeBLOSUM(void);

#ifdef __cplusplus
/* { */ }
#endif
#endif

/*
$Log: clustlib_pub.h,v $
Revision 1.15  2002/10/29 14:08:23  michel
moved BSFree out of conditional statement

Revision 1.14  2002/10/26 16:20:22  michel
minor fix

Revision 1.13  2002/10/21 17:48:38  michel
Changes to alignment logging - added loop & domain info

Revision 1.12  2002/06/25 23:03:32  michel
Added mkmdl and usemdl to windows gcpp workspace
Added function to clustal to apply new %id and %occ thresholds and
	extract alignment info (%id,%occ,domlen,etc) from structure models
Fixed gcpp/usemdl counting bugs and taxonomy calls

Revision 1.11  2002/06/14 22:47:56  michel
Addition of parameter

Revision 1.10  2002/03/30 22:02:16  michel
modified GetSeqInfoFromStruc for other types of residue info

Revision 1.9  2001/11/26 15:48:39  feldman
Fixed missing comment closer

Revision 1.8  2001/11/23 23:40:07  michel
added gilgor database formatting, modification to sequence to structure alignment function

Revision 1.7  2001/11/09 21:50:42  michel
Added domain ok list return on sequence to structure alignment

Revision 1.6  2001/09/09 05:38:13  michel
Added more returns for Sequencetostructure function
library can be linked by seqhound with linker flag to suppress multiply defined symbols

Revision 1.5  2001/08/28 21:48:57  michel
Removed printfs for ErrPostEx (set your message/log level appropriately)
Added paths for library if compiled with seqhound

Revision 1.4  2001/08/28 16:59:20  hogue
Fixed memory bug and changed function parameters

Revision 1.3  2001/06/27 01:37:01  michel
Added old slapamol functions to library with modifications and compile time
flags for seqhound or entrez access to sequences/structures

Revision 1.2  2001/06/26 20:40:09  feldman
Added makefile for UNIX,
and removed some functions that should not be there

Revision 1.1  2001/06/26 20:18:30  feldman
CLUSTAL/SLRI CSAN library

*/
