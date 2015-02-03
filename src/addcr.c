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

/* converts a UNIX text file to a DOS text file */

#include <stdio.h>

#define BUFSZ 16384

char tmpl[1024];
char buf[BUFSZ];
FILE *f,*g;
int len;

int main(int argc, char **argv)
{
	if (argc<2) {
		printf("Usage: %s [text file]\n",argv[0]);
		return 1;
	}
	if ((f=fopen(argv[1],"rb"))==NULL) {
		printf("Input file not found\n");
		return 1;
	}
	strcpy(tmpl,"./tmppoopXXXXXX");
	mkstemp(tmpl);
	g=fopen(tmpl,"wb");
	while (fgets(buf,BUFSZ,f)!=NULL) {
		len=strlen(buf)-1;
		if (buf[len]=='\n') {
			if (buf[len-1]!='\r') {
				buf[len]='\r';
				buf[len+1]='\n';
				buf[len+2]='\0';
			}
		}
		fputs(buf,g);
	}
	fclose(f);
	fclose(g);
	remove(argv[1]);
	f=fopen(tmpl,"rb");
	g=fopen(argv[1],"wb");
	while (fgets(buf,BUFSZ,f)!=NULL) {
		fputs(buf,g);
	}
	fclose(f);
	fclose(g);
	remove(tmpl);
	return 0;
}

/*
$Log: addcr.c,v $
Revision 1.1  2001/03/07 22:11:40  feldman
Convert UNIX text to DOS text.. WOW!


*/
