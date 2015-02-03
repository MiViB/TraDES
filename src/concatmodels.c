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


/* includes */
#include <hfprogs.h>

/* Global Variables */
Args foldargs[] = {{"Input MMDB File Base Name (NO EXTENSION)","protein",NULL,NULL,
         TRUE,'i',ARG_FILE_IN,0.0,0,NULL},{"Number of steps to concatenate","10",
	 "1","2147483647",TRUE,'n',ARG_INT,0.0,0,NULL}};

Int2 Main()
{
	Int2 cnt,cnt2,cnt3;
        Char buf[MAXCOL];
        Char tmpbuf[MAXCOL];
        Char remain[MAXCOL];
	char ch;
	Char infile[FILENAME_MAX];
	Char outfile[FILENAME_MAX];
	Char curfile[FILENAME_MAX];
	FILE *f, *f1, *g;

	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
        ErrSetMessageLevel(SEV_FATAL);
        ErrSetLogLevel(SEV_ERROR);
        if (!GetArgs("concat models",DIM(foldargs),foldargs))
                return 1;       /* insufficient arguments entered */
	StringCpy(infile,foldargs[0].strvalue);
	StringCpy(outfile,foldargs[0].strvalue);
	StringCat(infile,".prt");
	StringCat(outfile,"_movie.prt");
	if ((f=FileOpen(infile,"r"))==NULL) {
		ErrPostEx(SEV_FATAL,1,1,"Unable to open input file %s",infile);
		return 1;
	}	
	if ((g=FileOpen(outfile,"w"))==NULL) {
		ErrPostEx(SEV_FATAL,1,1,"Unable to open output file %s",outfile);
		return 1;
	}	
	do {
    if ((FileGets(buf,MAXCOL,f))==NULL) {
			ErrPostEx(SEV_FATAL,1,1,"No models in input file");
			return 1;
		}
		FilePuts(buf,g);
		sscanf(buf,"%s %c",tmpbuf,&ch);
	} while ((ch!='{') || StringCmp(tmpbuf,"model"));
       	FileGets(buf,MAXCOL,f);
	FilePuts(buf,g);
	cnt=1;
	while (cnt>0) {
        	FileGets(buf,MAXCOL,f);
		cnt2=0;
		while (buf[cnt2]) {
			if (buf[cnt2]=='{')
				cnt++;
			if (buf[cnt2]=='}')
				cnt--;
			cnt2++;
			if (cnt<1) {
				StringCpy(remain,&buf[cnt2]);
				buf[cnt2]=' ';
				buf[cnt2+1]=',';
				buf[cnt2+2]='\n';
				buf[cnt2+3]=0;
				cnt2+=3;
			}
		}
		FilePuts(buf,g);
	}
	for (cnt3=1;cnt3<=foldargs[1].intvalue;cnt3++) {
		sprintf(curfile,"%s_%d.prt",foldargs[0].strvalue,cnt3);
		if ((f1=FileOpen(curfile,"r"))==NULL) {
			ErrPostEx(SEV_FATAL,1,1,"Unable to open input file %s",curfile);
			return 1;
		}	
		do {
	        	if ((FileGets(buf,MAXCOL,f1))==NULL) {
				ErrPostEx(SEV_FATAL,1,1,"No models in input file");
				return 1;
			}
			sscanf(buf,"%s %c",tmpbuf,&ch);
		} while ((ch!='{') || StringCmp(tmpbuf,"model"));
	       	FileGets(buf,MAXCOL,f1);
		FilePuts(buf,g);
		cnt=1;
		while (cnt>0) {
	        	FileGets(buf,MAXCOL,f1);
			sscanf(buf,"%s %c",tmpbuf,&ch);
			if ((ch=='4') && !StringCmp(tmpbuf,"id")) {
				cnt2=0;
				while (buf[cnt2]!='4') cnt2++;
				buf[cnt2]=0;
				FilePuts(buf,g);
				fprintf(g,"%d ,\n",cnt3+4);
				continue;
			}
			cnt2=0;
			while (buf[cnt2]) {
				if (buf[cnt2]=='{')
					cnt++;
				if (buf[cnt2]=='}')
					cnt--;
				cnt2++;
				if (cnt<1) {
					if (cnt3<foldargs[1].intvalue) {
						buf[cnt2]=' ';
						buf[cnt2+1]=',';
						buf[cnt2+2]='\n';
						buf[cnt2+3]=0;
						cnt2+=3;
					}
					else {
						buf[cnt2]=0;
					}
				}
			}
			FilePuts(buf,g);
		}
		FileClose(f1);
	}
	FilePuts(remain,g);
       	while ((FileGets(buf,MAXCOL,f))!=NULL)
	FilePuts(buf,g);
	FileClose(f);
	FileClose(g);
	return 1;
}


/*  
$Log: concatmodels.c,v $
Revision 1.5  2001/03/14 16:25:53  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.4  2000/12/01 19:05:09  feldman
oops, typo

Revision 1.3  2000/12/01 19:04:14  feldman
Fixed potential bug caused by using fprintf - use fputs instead

Revision 1.2  2000/07/06 15:29:35  feldman
-- Updated old makefiles to compile with newest toolkit and
directory structure
-- added needed functions to hfprogs.h
-- update some executables to include hfprogs.h instead of mmdbtraj.h
-- replaced all instances of BiostrucAsnGet with MIMEBiostrucAsnGet
which can read with MIME biostrucs (v2.0) or normal ones (v1.0), and
as a result, foldtraj uses a new skel.prt, to result in a MIME biostruc
as its initial input

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

