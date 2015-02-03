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


/* 
   CWVH June 2012 re-written to use pure ASN.1 operators and bypass MMDAB-API 
   to ensure there is a way of extracting ASN.1 into ascii for any hand editing necessary 
   also allows us to chain together val2prt and prt2val to strip out extra coordinates
   or convert Biounits on the prt2val side.
*/

/* includes */
#include "hfprogs.h"


Args Rargs[4] = { 
      {"Input 3d Asn.1 Structure File Name (extension ignored, looks for .val .cn3):",NULL,NULL, NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
      {"Unwrap MIME and output bare Biostruc (t/F)","FALSE",NULL,NULL,TRUE,'u',ARG_BOOLEAN,0.0,0,NULL},
      {"Use Cn3D 4.3 naming: {file}_prt.c3d (t/F)","FALSE",NULL,NULL,TRUE,'c',ARG_BOOLEAN,0.0,0,NULL},
      {"Strip Secondary Structure Annotation from Biostruc (t/F)","FALSE",NULL,NULL,TRUE,'s',ARG_BOOLEAN,0.0,0,NULL} };

Int2 Main()
{
	Char fnamin[255];
	Char fnamout[255];
	Int2    iDotLen = 0;
	AsnIoPtr aip=NULL;
	NcbiMimeAsn1Ptr nmap=NULL;
	BiostrucPtr bsp=NULL;
	BiostrucSeqPtr bssp=NULL;
	ErrSev esMsg,esLog,esFatal;
	Byte bSaveMode;


        ErrSetLogfile("error_bin2prt.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
	if (!GetArgs("bin2prt Program: Converts Binary Asn.1 files into ASCII",4,Rargs))
		return 1;


	
/* be forgiving about adding extensions 
  CWVH  June 2012
    No-extension is user frustration! offer forgiveness for typing the extention , punch a null overtop of the dot  
*/

        iDotLen = 0;
	if (Rargs[0].strvalue != NULL) {
           iDotLen = (Int2) StringLen(Rargs[0].strvalue);
           if (iDotLen > 5) {
		if (Rargs[0].strvalue[iDotLen - 4] == '.') {
	        Rargs[0].strvalue[iDotLen - 4] = '\0';      
		}
           }
	}
/* Check input file exist */
       StringCpy(fnamin, Rargs[0].strvalue);
       StringCat(fnamin,".val");
        if (FileLength(fnamin) == 0) 
	  {
              StringCpy(fnamin, Rargs[0].strvalue);
              StringCat(fnamin, ".cn3");
              if (FileLength(fnamin) == 0) {
		ErrPostEx(SEV_FATAL,2,1,"Unable to find input %s.val or %s.cn3 file",Rargs[0].strvalue,Rargs[0].strvalue);
		return 3;
	     }
       	 }
/* name output file with extension */
	StringCpy(fnamout,Rargs[0].strvalue);
	if (Rargs[2].intvalue == TRUE) 
	   StringCat(fnamout,"_prt.cn3");
	else
	  StringCat(fnamout,".prt");
	
	 if (Rargs[3].intvalue == 0)
		bSaveMode = (Byte) SAVE_ASCII;
	 else 
		bSaveMode = (Byte) ((Byte) SAVE_ASCII | (Byte) NOT_FEATURES);

if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI");
		return 2;
	}


/* load an ASN.1 file *.val or *.c3d, write out Biostruc only */

	aip=AsnIoOpen(fnamin,"rb");
	if (aip==NULL) {
		ErrPostEx(SEV_FATAL,2,1,"Unable open binary ASN.1 stream in file %s",fnamin);
		return 2;
	}

/* first try biostruc load */
 /*	printf("try biostruc\n");  */
	esMsg=ErrGetMessageLevel();
	esLog=ErrGetLogLevel();
	esFatal=ErrGetFatalLevel();
        ErrSetMessageLevel(SEV_MAX);
        ErrSetLogLevel(SEV_MAX);
  	ErrSetFatalLevel(SEV_MAX);
	
 	bsp=BiostrucAsnRead(aip,NULL);
	AsnIoClose(aip);

  	ErrSetMessageLevel(esMsg);
  	ErrSetLogLevel(esLog);
  	ErrSetFatalLevel(esFatal);

	if (bsp==NULL) {
/*		printf("try mime\n"); */
		/* then try NCBIMime load */
		aip=NULL;
		aip=AsnIoOpen(fnamin,"rb");
  		ErrSetMessageLevel(SEV_MAX);
  		ErrSetLogLevel(SEV_MAX);
  		ErrSetFatalLevel(SEV_MAX);

		nmap=NcbiMimeAsn1AsnRead(aip,NULL);

  		ErrSetMessageLevel(esMsg);
  		ErrSetLogLevel(esLog);
  		ErrSetFatalLevel(esFatal);
		AsnIoClose(aip);


        	if (nmap!=NULL) {
		/* got an NCBI mime */
			if (nmap->choice!=NcbiMimeAsn1_strucseq && nmap->choice!=NcbiMimeAsn1_strucseqs) {
				/* wrong MIME type */
				nmap=NcbiMimeAsn1Free(nmap);
		    		ErrPostEx(SEV_ERROR,6,9,"MIME-type wrapper is not a strucseq or strucseqs - no Biostruc to extract from %s",fnamin);
				return 1;
			}
			if (Rargs[1].intvalue == TRUE) { 
		/* we want to unwrap the mime and leave the bsp */
				bssp=(BiostrucSeqPtr)(nmap->data.ptrvalue);
				bsp=bssp->structure;
                		bssp->structure = NULL;
				nmap=NcbiMimeAsn1Free(nmap); /* discard the wrapper */
				nmap = NULL;
				bssp = NULL;
			}
		}
	}
        if (bsp != NULL) {
		printf("Writing Biostruc: %s\n",fnamout);
		/* write to prt using call in mmdbapi4.c */
		if (WriteOutBiostruc(bsp, fnamout,  bSaveMode) != FALSE);
			return 0;
                return 1;
	}
	if (nmap != NULL) { /* write out Mime-type as read in with no changes  */
		printf("Writing MIME-Biostruc: %s\n",fnamout);

		aip=AsnIoOpen(fnamout,"w");
		if (aip==NULL) {
			ErrPostEx(SEV_FATAL,2,1,"Unable open binary ASN.1 stream in file %s",fnamin);
		return 2;
		}
                NcbiMimeAsn1AsnWrite(nmap,aip,NULL);
		AsnIoClose(aip);
		nmap=NcbiMimeAsn1Free(nmap); 
	}
 


	CloseMMDBAPI(); 
 	return 0;

}




/*  
$Log: val2prt.c,v $
Revision 1.6  2002/12/09 16:00:28  feldman
Made work with more types of Biostrucs

Revision 1.5  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.4  2000/12/15 00:06:41  feldman
Now Val2Trj and unfoldtraj works correctly for Phi-Psi walk (I think)

Revision 1.3  2000/07/14 20:09:25  feldman
Added parameter for MIMEBiostrucAsnGet for Bioseq

Revision 1.2  2000/07/06 15:29:37  feldman
-- Updated old makefiles to compile with newest toolkit and
directory structure
-- added needed functions to hfprogs.h
-- update some executables to include hfprogs.h instead of mmdbtraj.h
-- replaced all instances of BiostrucAsnGet with MIMEBiostrucAsnGet
which can read with MIME biostrucs (v2.0) or normal ones (v1.0), and
as a result, foldtraj uses a new skel.prt, to result in a MIME biostruc
as its initial input

Revision 1.1.1.1  2000/06/09 18:14:05  feldman
TraDES project

*/

