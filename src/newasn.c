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
#include <foldtrajlib.h>

/* creates an ASN.1 file (named fnam) containing the a-carbon skeleton for the
   protein sequence described by seq */
void BuildSkelASN(CharPtr encseq,CharPtr fnam)
{
	FILE *f,*g;
	Char buf[MAXCOL];
	Char tmpbuf[80];
	CharPtr seqbuf,seqptr,endptr;
	CharPtr seq;
	CharPtr xxseq;
	struct tm PNTR curTime;
	Char ftmp[PATH_MAX];
	PEAS peas;
	Int2 cnt,seqlen,encidx;
	Int2 dictindex[25]={1,0,13,10,19,40,22,25,28,0,34,31,37,7,0,43,16,4,46,49,0,58,52,0,55};
	Int2 dictNID[25]={6,0,7,9,10,17,4,13,10,0,9,10,8,7,0,7,8,10,7,9,0,9,20,0,17};
	Int2 dictCID[25]={1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1};

	seq=DecodeSequence(encseq,EXTAA_PARENT);
	xxseq=DecodeSequence(encseq,EXTAA_X);
	seqlen=StringLen(seq);
	if (seqlen<2) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Sequence must be at least 2 residues long");
		return;
	}
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,"skel.prt");
	if ((CheckMD5(ftmp,"4f140d77a7c0d2c391636702e18130d6")==FALSE) && (CheckMD5(ftmp,"a4d094f6c817826b7e63682741573e51")==FALSE)) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Error - file failed checksum %s",ftmp);
		return;
	}
	if ((f=FileOpen(ftmp,"r"))==NULL) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Error - unable to open %s",ftmp);
		return;
	}
	curTime=(struct tm PNTR)MemNew(sizeof(struct tm));
	if ((g=FileOpen(fnam,"w"))==NULL) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Error - unable to open temporary file %s for output",fnam);
		return;
	}
	/* use current date */
	GetDayTime(curTime);
	while (FileGets(buf,MAXCOL,f)!=NULL) {
		if (!StringNCmp(buf,"***insert",9)) {
			if (!StringNCmp(buf,"***insert date } here",21)) {
				fprintf(g,"              year %d ,\n", (int) 1900+curTime->tm_year);
				fprintf(g,"              month %d ,\n", (int) 1+curTime->tm_mon);
				fprintf(g,"              day %d } ,\n", (int) curTime->tm_mday);
			}
			else if (!StringNCmp(buf,"***insert date } } } here",25)) {
				fprintf(g,"              year %d ,\n", (int) 1900+curTime->tm_year);
				fprintf(g,"              month %d ,\n",(int) 1+curTime->tm_mon);
				fprintf(g,"              day %d } } } ,\n", (int) curTime->tm_mday);
			}
			else if (!StringNCmp(buf,"***insert date } } } } here",27)) {
				fprintf(g,"              year %d ,\n", (int) 1900+curTime->tm_year);
				fprintf(g,"              month %d ,\n",(int) 1+curTime->tm_mon);
				fprintf(g,"              day %d } } } } ,\n",(int) curTime->tm_mday);
			}
			else if (!StringNCmp(buf,"***insert residues here",23)) {
				encidx=0;
				for (cnt=0;cnt<seqlen;cnt++) {
					fprintf(g,"          {\n            id %d ,\n",(int) cnt+1);
					fprintf(g,"            name \"%4d \" ,\n            residue-graph\n",(int) cnt+1);
					fprintf(g,"              standard {\n                biostruc-residue-graph-set-id\n");
					fprintf(g,"                  other-database {\n                    db \"Standard residue dictionary\" ,\n");
					fprintf(g,"                    tag\n                      id 1 } ,\n");
					if (encseq[encidx]=='*') {
						peas=GetExtAAInfoEnc(1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]-1111*'0');
						/* special amino acid */
/*						if (cnt==0 && peas->modlocation=='-')
							fprintf(g,"                residue-graph-id %d } }",(int) (1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]+2-1111*'0'));
						else if (cnt==seqlen-1 && peas->modlocation=='-')
							fprintf(g,"                residue-graph-id %d } } }",(int) (1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]+1-1111*'0'));
						else*/
						if (cnt==seqlen-1)
							fprintf(g,"                residue-graph-id %d } } }",(int) (1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]-1111*'0'));
						else
							fprintf(g,"                residue-graph-id %d } }",(int) (1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]-1111*'0'));
						peas=MemFree(peas);
						encidx+=6;
					}
					else {
						/* use N/C-terminal residue if required */
						if (cnt==0)
							fprintf(g,"                residue-graph-id %d } }",(int) dictindex[seq[cnt]-'A']+2);
						else if (cnt==seqlen-1)
							fprintf(g,"                residue-graph-id %d } } }", (int) dictindex[seq[cnt]-'A']+1);
						else
							fprintf(g,"                residue-graph-id %d } }", (int) dictindex[seq[cnt]-'A']);
						encidx++;
					}
					fprintf(g," ,\n");
				}
				fprintf(g,"        inter-residue-bonds {\n");
				encidx=0;
				for (cnt=0;cnt<seqlen-1;cnt++) {
					fprintf(g,"          {\n            atom-id-1 {\n");
					fprintf(g,"              molecule-id 1 ,\n              residue-id %d ,\n",(int) cnt+1);
					if (encseq[encidx]=='*') {
						peas=GetExtAAInfoEnc(1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]-1111*'0');
						if ((cnt==0) && (peas->modlocation=='-'))
							fprintf(g,"              atom-id %d } ,\n            atom-id-2 {\n", (int)(peas->cidx[2]));
						else
							fprintf(g,"              atom-id %d } ,\n            atom-id-2 {\n", (int)(peas->cidx[0]));
						peas=MemFree(peas);
						encidx+=6;					
					}
					else {
						fprintf(g,"              atom-id %d } ,\n            atom-id-2 {\n", (int) dictCID[seq[cnt]-'A']);
						encidx++;
					}
					fprintf(g,"              molecule-id 1 ,\n              residue-id %d ,\n", (int) cnt+2);
					/* advanced encidx already otherwise would look at encidx+1 or whatever */
					if (encseq[encidx]=='*') {
						peas=GetExtAAInfoEnc(1000*encseq[encidx+1]+100*encseq[encidx+2]+10*encseq[encidx+3]+encseq[encidx+4]-1111*'0');
						/* last peptide bond */
						if ((cnt==seqlen-2) && (peas->modlocation=='-'))
							fprintf(g,"              atom-id %d } }",(int)(peas->nidx[1]));
						else
							fprintf(g,"              atom-id %d } }",(int)(peas->nidx[0]));
						peas=MemFree(peas);
					}
					else {
						if ((cnt==seqlen-2) && (seq[cnt+1]=='P'))
							fprintf(g,"              atom-id %d } }",(int) dictNID[seq[cnt+1]-'A']+1);
						else
							fprintf(g,"              atom-id %d } }",(int) dictNID[seq[cnt+1]-'A']);
					}
					if (cnt==seqlen-2) fprintf(g," } } } }");
					fprintf(g," ,\n");
				}
			}
			/* give a-carbons arbitrary co-ordinates along x-axis; these
			   will be changed later but ensure the graph is instantiated
			   correctly */
			else if (!StringNCmp(buf,"***insert co-ordinates here",27)) {
				fprintf(g,"                  number-of-points %d ,\n                  atoms {\n", (int) seqlen);
				fprintf(g,"                    number-of-ptrs %d ,\n                    molecule-ids {\n", (int) seqlen);
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      1 ,\n");
				fprintf(g,"                      1 } ,\n             residue-ids {\n");
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      %d ,\n", (int) cnt+1);
				fprintf(g,"                      %d } ,\n                    atom-ids {\n", (int) seqlen);
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      2 ,\n");
				fprintf(g,"                      2 } } ,\n                    sites {\n");
				fprintf(g,"                    scale-factor 1000 ,\n                    x {\n");
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      %ld , \n", (long) (5000*(Int4)cnt));
				fprintf(g,"                      %ld } ,\n                    y {\n",(long) (5000*(Int4)(seqlen-1)));
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      0 ,\n");
				fprintf(g,"                      0 } ,\n                    z {\n");
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      0 ,\n");
				fprintf(g,"                      0 } } ,\n                  temperature-factors\n");
				fprintf(g,"                    isotropic {\n                    scale-factor 1000 ,\n");
				fprintf(g,"                      b {\n");
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      0 ,\n");
				fprintf(g,"                      0 } } ,\n                  occupancies {\n");
				fprintf(g,"                    scale-factor 1000 ,\n                    o {\n");
				for (cnt=0;cnt<seqlen-1;cnt++) 
					fprintf(g,"                      1000 ,\n");
				/* attempt to update to ASN biostrucs 2.0 */
				fprintf(g,"                      1000 } } } } } } } } ,\n");
				/* end attempt */
/*				fprintf(g,"                      1000 } } } } } } } }\n");
*/			
			}
			else if (!StringNCmp(buf,"***insert residue numbering here",32)) {
				fprintf(g,"            num %d ,\n            names {\n",(int)seqlen);
				for (cnt=1;cnt<seqlen;cnt++)
					fprintf(g,"              \"%d\" ,\n",(int)cnt);
					fprintf(g,"              \"%d\" } } ,\n",(int)seqlen);
			}
			else if (!StringNCmp(buf,"***insert sequence here",23)) {
					fprintf(g,"        length %d ,\n        seq-data\n",(int)seqlen);
					seqbuf=(CharPtr)MemNew(sizeof(Char)*(seqlen+50));
					sprintf(seqbuf,"          ncbieaa \"%s\" } } } }",xxseq);
					seqptr=seqbuf;
					endptr=seqbuf+StringLen(seqbuf);
					do {
						StringNCpy(tmpbuf,seqptr,78);
						tmpbuf[78]='\0';
						fprintf(g,"%s\n",tmpbuf);
						seqptr+=78;
					} while (seqptr<endptr);
					seqbuf=MemFree(seqbuf);
			}
		}
		else
			fprintf(g,buf);
	}
	FileClose(f);
	FileClose(g);
	curTime=MemFree(curTime);
	seq=MemFree(seq);
	xxseq=MemFree(xxseq);
}

void BuildSkelASNDict(Int2 numrot,CharPtr fnam)
{
	FILE *f,*g;
	Char buf[MAXCOL];
	struct tm PNTR curTime;
	Int2 cnt,cnt2;

	if ((f=FileOpen("mmdbrdic.prt","r"))==NULL) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Error - unable to open mmdbrdic.prt");
		return;
	}
	if ((g=FileOpen(fnam,"w"))==NULL) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Error - unable to open temporary file for output");
		return;
	}
	curTime=(struct tm PNTR)MemNew(sizeof(struct tm));
	/* use current date */
	GetDayTime(curTime);
	while (FileGets(buf,MAXCOL,f)!=NULL) {
		if (!StringNCmp(buf,"***insert",9)) {
			for (cnt=1;cnt<=numrot;cnt++) {
				fprintf(g,"      {\n        id %d ,\n",cnt+3);
				fprintf(g,"        type pdb-model ,\n        descr {\n");
				fprintf(g,"          name \"Model %d\" ,\n          pdb-method \"Handmade\" ,\n",cnt);
				fprintf(g,"          pdb-comment \"created %d/%d/%d\" } ,\n        model-space {\n",(int)1+curTime->tm_mon,(int)curTime->tm_mday,(int)1900+curTime->tm_year);
				fprintf(g,"          coordinate-units angstroms ,\n          thermal-factor-units b ,\n");
				fprintf(g,"          occupancy-factor-units fractional } ,\n");
				fprintf(g,"        model-coordinates {\n          {\n");
				fprintf(g,"            id 1 ,\n            coordinates\n");
				fprintf(g,"              literal\n                atomic {\n");
				fprintf(g,"                  number-of-points 20 ,\n                  atoms {\n");
				fprintf(g,"                    number-of-ptrs 20 ,\n                    molecule-ids {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      %d ,\n",cnt2);
				fprintf(g,"                      20 } ,\n             residue-ids {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      1 ,\n");
				fprintf(g,"                      1 } ,\n                    atom-ids {\n");
				for (cnt2=1;cnt2<21;cnt2++) 
					fprintf(g,"                      2 ,\n");
				fprintf(g,"                      2 } } ,\n                    sites {\n");
				fprintf(g,"                    scale-factor 1000 ,\n                    x {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      %ld , \n",(long)cnt2*10);
				fprintf(g,"                      200 } ,\n                    y {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      1000 ,\n");
				fprintf(g,"                      1000 } ,\n                    z {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      1000 ,\n");
				fprintf(g,"                      1000 } } ,\n                  temperature-factors\n");
				fprintf(g,"                    isotropic {\n");
				fprintf(g,"                      scale-factor 1000 ,\n                      b {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      0 ,\n");
				fprintf(g,"                      0 } } ,\n                  occupancies {\n");
				fprintf(g,"                    scale-factor 1000 ,\n                    o {\n");
				for (cnt2=1;cnt2<20;cnt2++) 
					fprintf(g,"                      1000 ,\n");
				if (cnt<numrot)
					fprintf(g,"                      1000 } } } } } } ,\n");
				else
					fprintf(g,"                      1000 } } } } } } } }\n");
			}
		}
		else
			fprintf(g,buf);
	}
	FileClose(f);
	FileClose(g);
	curTime=MemFree(curTime);
}



/*  
$Log: newasn.c,v $
Revision 1.18  2002/09/26 13:23:22  michel
Moved BuildMIMEBiostruc to mmdbtrajlib

Revision 1.17  2002/02/25 22:08:01  feldman
Added checksumming of text files potentials, cbdata and skel.prt
If checksum fails, get error
Changed bailing from foldtrajlite to exit curses first if in it

Revision 1.16  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.15  2001/05/03 22:02:36  feldman
Fixed bioseq highlighting in Cn3D

Revision 1.14  2001/04/26 19:49:09  feldman
Fixed some minor Bioseq bugs and potential bugs

Revision 1.13  2001/04/26 16:01:07  feldman
Try to keep Bioseq descriptors whenever possible, but renumber
bioseq residues when making an NCBI MIME ASN1 structure

Revision 1.12  2001/04/15 18:40:32  feldman
fixed memory leak

Revision 1.11  2001/04/15 06:14:25  feldman
fixed memory leak in nmap building

Revision 1.10  2001/03/23 18:44:18  feldman
Integrated foldtraj into vistraj (!)
and cleaned up a few minor bugs

Revision 1.9  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.8  2001/02/08 22:28:11  feldman
Fixed important inconsistency in modified amino acids - now for
encoded sequence (*xxxxA) the xxxx is always the dictionary index
of the residue, thus we no longer need to add 1 or 2 to it in
cases where it is found at the N- or C-terminus; adding and
deleting residues with AlterResiduesInSequence will update any
affected residues correctly as well

Revision 1.7  2000/10/13 16:53:50  feldman
Minor error message change

Revision 1.6  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.5  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.4  2000/07/06 15:29:38  feldman
-- Updated old makefiles to compile with newest toolkit and
directory structure
-- added needed functions to hfprogs.h
-- update some executables to include hfprogs.h instead of mmdbtraj.h
-- replaced all instances of BiostrucAsnGet with MIMEBiostrucAsnGet
which can read with MIME biostrucs (v2.0) or normal ones (v1.0), and
as a result, foldtraj uses a new skel.prt, to result in a MIME biostruc
as its initial input

Revision 1.3  2000/06/22 16:33:08  feldman
minor bugfix and added VDW radii to hfprogs

Revision 1.2  2000/06/22 13:05:01  feldman
Fixed minor bugs with log file nameing,
fixed error in extended residue calculation and DSSP usage
added alterresidue function

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

