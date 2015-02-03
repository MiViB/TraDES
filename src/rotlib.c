/* Portions Copyright (c) 2007-2012 Hogue Laboratory, National University of Singapore
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



	    Rotamer library module 

  AUTHORS:  
            Howard J. Feldman           (feldman@mshri.on.ca)
            and Christopher W.V. Hogue  (hogue@mshri.on.ca)
	    Adapted from code by Roland Dunbrack used in SCWRL
*/


#include <mmdbtraj.h>
#include "rotlib.h"

static Int4 numrot[]={0,81,9,9,3,27,27,0,6,9,9,81,27,6,2,3,3,9,6,3};
static CharPtr SCWRLRotLib=NULL;

void FreeRotLib(void)
{
  if (SCWRLRotLib!=NULL)
    SCWRLRotLib=MemFree(SCWRLRotLib);
}

/* bunzip the SCWRL library into a global memory buffer */
TrajErr LoadRotLib(void)
{
  CharPtr bzbuf;
  Int2 err;
  Int4 bzerr;
  FILE *f;
  Char buf[PATH_MAX];
  unsigned int actsize;

  sprintf(buf,"%s%s",CFG_local_datafilepath,NEWSCWRL_FNAM);
  if ((f=FileOpen(buf,"rb"))==NULL) {
    ErrPostEx(SEV_ERROR,1,2,"Unable to open input file %s",buf);
    return ERR_FAIL;
  }
  bzbuf=(CharPtr)MemNew(sizeof(Char)*NEWSCWRL_BZ_SIZE);
  bzerr=FileRead(bzbuf,sizeof(Char),NEWSCWRL_BZ_SIZE,f);
  FileClose(f);
  if (bzerr!=NEWSCWRL_BZ_SIZE) {
    ErrPostEx(SEV_ERROR,2,2,"Input file %s possibly corrupt or in use, only %ld bytes could be read.",buf,(long int)bzerr);
    MemFree(bzbuf);
    return ERR_FAIL;
  }
  SCWRLRotLib=(CharPtr)MemNew(sizeof(Char)*NEWSCWRL_SIZE);
  actsize=(unsigned int)NEWSCWRL_SIZE;
  if ((err=BZ2_bzBuffToBuffDecompress((char *)SCWRLRotLib,&actsize,(char *)bzbuf,(unsigned int)NEWSCWRL_BZ_SIZE,0,0))!=BZ_OK) {
    ErrPostEx(SEV_ERROR,err,0,"Bzip decompression of SCWRL library error occurred");
    MemFree(bzbuf);
    MemFree(SCWRLRotLib);
    return ERR_FAIL;
  }
  bzbuf=MemFree(bzbuf);
  if (actsize!=(unsigned long int)NEWSCWRL_SIZE) {
    ErrPostEx(SEV_ERROR,1,3,"SCWRL library wrong size, integrity error");
    MemFree(SCWRLRotLib);
    return ERR_FAIL;
  }
  return ERR_SUCCESS;
}

TrajErr LoadOldRotLib(void)
{
  CharPtr bzbuf;
  Int2 err;
  FILE *f;
  unsigned int actsize;

  if ((f=FileOpen(SCWRL_FNAM,"rb"))==NULL) {
    ErrPostEx(SEV_ERROR,1,2,"Unable to open input file %s",SCWRL_FNAM);
    return ERR_FAIL;
  }
  bzbuf=(CharPtr)MemNew(sizeof(Char)*SCWRL_BZ_SIZE);
  FileRead(bzbuf,sizeof(Char),SCWRL_BZ_SIZE,f);
  FileClose(f);
  SCWRLRotLib=(CharPtr)MemNew(sizeof(Char)*SCWRL_SIZE);
  actsize=(unsigned long int)SCWRL_SIZE;
  if ((err=BZ2_bzBuffToBuffDecompress((char *)SCWRLRotLib,&actsize,(char *)bzbuf,(unsigned int)SCWRL_BZ_SIZE,0,0))!=BZ_OK) {
    ErrPostEx(SEV_ERROR,err,0,"Bzip decompression of SCWRL library error occurred");
    MemFree(bzbuf);
    MemFree(SCWRLRotLib);
    return ERR_FAIL;
  }
  bzbuf=MemFree(bzbuf);
  if (actsize!=(unsigned long int)SCWRL_SIZE) {
    ErrPostEx(SEV_ERROR,1,3,"SCWRL library wrong size, integrity error");
    MemFree(SCWRLRotLib);
    return ERR_FAIL;
  }
  return ERR_SUCCESS;
}

TrajErr FlipRotRec(pnewrotlibrecord rec)
{
  if (sizeof(FloatLo)!=4) {
    ErrPostEx(SEV_ERROR,48,0,"%s%s format incompatible with this system, unable to continue",CFG_local_datafilepath,NEWSCWRL_FNAM);
    return ERR_FAIL;
  }
/*  rec->phi=FlipI4(&(rec->phi));
  rec->psi=FlipI4(&(rec->psi));
  rec->N=FlipI4(&(rec->N));
  rec->r1=FlipI4(&(rec->r1));
  rec->r2=FlipI4(&(rec->r2));
  rec->r3=FlipI4(&(rec->r3));
  rec->r4=FlipI4(&(rec->r4));*/
  rec->p=FlipF(&(rec->p));
  rec->chi1=FlipF(&(rec->chi1));
  rec->chi2=FlipF(&(rec->chi2));
  rec->chi3=FlipF(&(rec->chi3));
  rec->chi4=FlipF(&(rec->chi4));
  rec->chi1sd=FlipF(&(rec->chi1sd));
  rec->chi2sd=FlipF(&(rec->chi2sd));
  return ERR_SUCCESS;
}

TrajErr get_rotamers(PMGD pmgdHere,FloatLo phi,FloatLo psi,FloatLo arandnum,FloatLo PNTR chi1,FloatLo PNTR chi2,FloatLo PNTR chi3,FloatLo PNTR chi4,FloatLo PNTR chi1sd,FloatLo PNTR chi2sd)
{
  Int2 k,res;
  Int4 cnt,lphi,lpsi,aaseek;
  CharPtr rechere;
  newrotlibrecord record;
	
  /* number of rotamers in scwrlbin34.lib file */
/* temp hack */
/*Char tc;
k=(Int2)(arandnum*9.0);
tc=(pmgdHere->pcIUPAC)[0];
lphi=k/3;
lpsi=k%3;
*chi1=(FloatLo)(60.0+120.0*lphi);
if ((tc=='F') || (tc=='W') || (tc=='H') || (tc=='Y'))
*chi2=(FloatLo)(90.0+180.0*lpsi);
else
*chi2=(FloatLo)(60.0+120.0*lpsi);
while ((*chi1)>180.0)
(*chi1)-=360.0;
while ((*chi2)>180.0)
(*chi2)-=360.0;
*chi3=0.0;
*chi4=0.0;
if (tc=='R') {
*chi3=175.076;
*chi4=163.155;
}
if (tc=='Q') {
*chi3=13.2519;
}
if (tc=='M') {
*chi3=173.635;
}
if (tc=='K') {
*chi3=171.131;
*chi4=179.21;
}
if (tc=='E') {
*chi3=18.7903;
}
return;
*/

  /* Location of first record for each amino acid type in binary rotamer library
     RLD 8/6/97 */
  res=GetResIdxFromMGD(aalist,pmgdHere);
  aaseek=0L;
  for (cnt=0;cnt<res;cnt++)
    aaseek=aaseek+numrot[cnt]*1369L*sizeof(newrotlibrecord);
    /* 1369 = 37 * 37 Ramachandran space divisions */
  /* For N-terminal residues (no phi), we use phi=-60
     part of bbdep library, where N-terminal neighbors
     have little effect on rotamer distribution. RLD 7/8/97 */
  /* For C-terminal residues (no psi), we use psi=+60
     part of bbdep library, where C-terminal neighbors
     have little effect on rotamer distribution. RLD 7/8/97 */
  /* but only residues having both phi and psi defined will 
     be passed to this function */
  /* shift "zero" to -180, -180 */
  phi+=180;
  psi+=180;
  /* convert phi and psi to 0-36 range each */
  while (phi<0.0)
    phi+=360.0;
  while (psi<0.0)
    psi+=360.0;
  while (phi>360.0)
    phi-=360.0;
  while (psi>360.0)
    psi-=360.0;
  phi=floor((phi+5.0)/10.0);
  psi=floor((psi+5.0)/10.0);
  lphi=(Int4)phi;
  lpsi=(Int4)psi;
  /* Location in binary library is calculated from xxxseek (1st record 
     for residue xxx) + (37*lphi + lpsi)*(size for each phi,psi bin 
     for residue type xxx, xxxsize); E.g. there are cyssize=3 records 
     for each cysteine residue; for phi=-120, psi=60,
     this record is the (cysseek + (37*6 + 24)*3) th record.
     RLD 8/6/97 */
  rechere=SCWRLRotLib+aaseek+(37L*lphi+lpsi)*numrot[res]*sizeof(newrotlibrecord);
  k=0;
  do {
    if (k>=numrot[res]) {
      ErrPostEx(SEV_ERROR,14,1,"Rotamer probabilities sum to less than 1, randnum=%f, k=%d, res=%s, phi=%f, psi=%f",arandnum,k,pmgdHere->pcIUPAC,phi,psi);
      return ERR_FAIL;
    }
    MemCopy(&record,rechere,sizeof(newrotlibrecord));
    if (GetEndian()==ENDIAN_UNKNOWN)
	return ERR_FAIL;
    if (GetEndian()==ENDIAN_LITTLE) {
      if (FlipRotRec(&record)!=ERR_SUCCESS)
    	return ERR_FAIL;
    }
    arandnum-=record.p;
    k++;
    /* advance record pointer */
    rechere+=sizeof(newrotlibrecord);
  /* read records until probability is arandnum - thus it acts as a
     rotamer pdf of sorts - but allow this to be slightly above zero
     to account for possible rounding errors */
  } while (arandnum>PRECISION);
  *chi1=record.chi1;
  *chi2=record.chi2;
  *chi3=record.chi3;
  *chi4=record.chi4;
  *chi1sd=record.chi1sd;
  *chi2sd=record.chi2sd;
/*printf("%s phi:%f psi:%f 1:%f 2:%f 3:%f 4:%f p:%f sd:%f sd:%f\n",pmgdHere->pcIUPAC,phi,psi,record.chi1,record.chi2,record.chi3,record.chi4,record.p,record.chi1sd,record.chi2sd);*/
  /* ensure chi2 of Phe, Tyr and Asp all lie between -90 and +90,
     which is a standard convention */
  if ((res==StringCSpn(aalist,"Y")) || (res==StringCSpn(aalist,"F")) || (res==StringCSpn(aalist,"D"))) {
    while ((*chi2)<-90.0)
      (*chi2)+=180.0;
    while ((*chi2)>90.0)
      (*chi2)-=180.0;
  }
  /* close files */
  return ERR_SUCCESS;
}

/* perhaps not graceful, but should do the job */
TrajErr AddSDtoSCWRLLIB(void)
{
	CharPtr NewSCWRL;
	Int4 posold,posnew;
	FILE *f,*g;
	Char readbuf[MAXCOL];
	Char outname[PATH_MAX];
	Char residue[4];
	Char lastresidue[4];
	FloatLo sd11=0.0,sd12=0.0,sd13=0.0,psi,lastpsi,chi1sd,chi2sd;
	float sd11tmp,sd12tmp,sd13tmp,sdtmp;
	FloatLo sd2[320];
	Int4 r1,r2,r3,r4,cnt,resnum=-1;
	Int2 rotidx[]={0,81,90,99,102,129,156,162,171,180,261,288,294,296,299,302,311,317};
	
  NewSCWRL=(CharPtr)MemNew(sizeof(Char)*NEWSCWRL_SIZE);
  printf("Opening Old Rotamer Library..\n");
	if (LoadOldRotLib()!=ERR_SUCCESS) {
		MemFree(NewSCWRL);
		return ERR_FAIL;
	}
	lastpsi=5000.0;
	printf("%-8.3f %% done",0.0);
	fflush(stdout);
	StringCpy(lastresidue,"XXX");
	if ((g=FileOpen(BBINDEP_LIBNAME ,"r"))==NULL) {
		ErrPostEx(SEV_ERROR,1,4,"Unable to open input file %s",BBINDEP_LIBNAME);
		MemFree(NewSCWRL);
	        return ERR_FAIL;
	}
	/* seek to correct strarting point in this file */	
	do {
		FileGets(readbuf,MAXCOL,g);
	} while (StringStr(readbuf,"Res Rotamer")==NULL);
	FileGets(readbuf,MAXCOL,g);
	for (cnt=0;cnt<320;cnt++) {
		FileGets(readbuf,MAXCOL,g);
		while (sscanf(readbuf,"%s %*d %*d %*d %*d %*d %*d %*f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f",residue,&sdtmp)==EOF)
			FileGets(readbuf,MAXCOL,g);
		sd2[cnt]=(FloatLo)sdtmp;
	}
 	FileClose(g);
	if ((g=FileOpen(BBDEP_LIBNAME ,"r"))==NULL) {
		ErrPostEx(SEV_ERROR,1,4,"Unable to open input file %s",BBDEP_LIBNAME);
		MemFree(NewSCWRL);
	        return ERR_FAIL;
	}
	posold=0L;
	posnew=0L;
	do {
 		r1=((protlibrecord)(SCWRLRotLib+posold))->r1;
 		r2=((protlibrecord)(SCWRLRotLib+posold))->r2;
 		r3=((protlibrecord)(SCWRLRotLib+posold))->r3;
 		r4=((protlibrecord)(SCWRLRotLib+posold))->r4;
		if (GetEndian()==ENDIAN_UNKNOWN) {
			MemFree(NewSCWRL);
		        return ERR_FAIL;
		}
		if (GetEndian()==ENDIAN_LITTLE) {
		  r1=FlipI4(&r1);
		  r2=FlipI4(&r2);
		  r3=FlipI4(&r3);
		  r4=FlipI4(&r4);
		}
 		psi=(FloatLo)((protlibrecord)(SCWRLRotLib+posold))->psi;
 		/* check if new line in file to be read */
 		if (psi!=lastpsi) {
 			FileGets(readbuf,MAXCOL,g);
			sscanf(readbuf,"%s %*d %*d %*d %*f %*f %f %*f %*f %f %*f %*f %f",residue,&sd11tmp,&sd12tmp,&sd13tmp);
			sd11=(FloatLo)sd11tmp;
			sd12=(FloatLo)sd12tmp;
			sd13=(FloatLo)sd13tmp;
		}
		lastpsi=psi;
		if (StringCmp(residue,lastresidue)) resnum++;
		StringCpy(lastresidue,residue);
		((pnewrotlibrecord)(NewSCWRL+posnew))->p=((protlibrecord)(SCWRLRotLib+posold))->p;
 		((pnewrotlibrecord)(NewSCWRL+posnew))->chi1=((protlibrecord)(SCWRLRotLib+posold))->chi1;
 		((pnewrotlibrecord)(NewSCWRL+posnew))->chi2=((protlibrecord)(SCWRLRotLib+posold))->chi2;
 		((pnewrotlibrecord)(NewSCWRL+posnew))->chi3=((protlibrecord)(SCWRLRotLib+posold))->chi3;
 		((pnewrotlibrecord)(NewSCWRL+posnew))->chi4=((protlibrecord)(SCWRLRotLib+posold))->chi4;
 		if (r1==1)
	 		chi1sd=sd11;
	 	else if (r1==2)
	 		chi1sd=sd12;
	 	else /* rotamer 3 */
	 		chi1sd=sd13;
		chi2sd=0.0;
		if (r2!=0) {
			/* HIS, PHE have 6 rather than 9 rotamers just to complicate things,
			   (and PRO has 2, not 3) */
			if (r4!=0)
				chi2sd=sd2[rotidx[resnum]+27*(r1-1)+9*(r2-1)+3*(r3-1)+(r4-1)];
			else if (r3!=0)
				chi2sd=sd2[rotidx[resnum]+9*(r1-1)+3*(r2-1)+(r3-1)];
			else {
				if ((!StringCmp(residue,"HIS")) || (!StringCmp(residue,"PHE")))
					chi2sd=sd2[rotidx[resnum]+2*(r1-1)+(r2-1)];
				else if (!StringCmp(residue,"PRO"))
					chi2sd=sd2[rotidx[resnum]+(r1-1)+(r2-1)];
				else
					chi2sd=sd2[rotidx[resnum]+3*(r1-1)+(r2-1)];
			}
		}
		if (GetEndian()==ENDIAN_UNKNOWN) {
			MemFree(NewSCWRL);
		        return ERR_FAIL;
		}
    		if (GetEndian()==ENDIAN_LITTLE) {
			chi1sd=FlipF(&chi1sd);
			chi2sd=FlipF(&chi2sd);
		}
		((pnewrotlibrecord)(NewSCWRL+posnew))->chi1sd=chi1sd;
		((pnewrotlibrecord)(NewSCWRL+posnew))->chi2sd=chi2sd;
		posold+=sizeof(rotlibrecord);
		posnew+=sizeof(newrotlibrecord); 		
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%-8.3f %% done",100.0*(FloatLo)posold/SCWRL_SIZE);
		fflush(stdout);
 	} while (posold<SCWRL_SIZE);
 	FileClose(g);
 	printf("\n\n");
 	StringCpy(outname,CFG_local_datafilepath);
	StringCat(outname,NEWSCWRL_FNAM);
 	/* remove extension */
 	outname[StringLen(outname)-4]=0;
  if ((f=FileOpen(outname,"wb"))==NULL) {
    ErrPostEx(SEV_ERROR,1,2,"Unable to open output file %s",outname);
    MemFree(NewSCWRL);
    return ERR_FAIL;
  }
  FileWrite(NewSCWRL,sizeof(Char),NEWSCWRL_SIZE,f);
  FileClose(f);
  NewSCWRL=MemFree(NewSCWRL);
  return ERR_SUCCESS;
}



/*  
$Log: rotlib.c,v $
Revision 1.8  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.7  2001/04/06 15:26:01  feldman
moved flipping functions to SLRI lib

Revision 1.6  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.5  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.4  2001/01/18 18:10:38  feldman
Added bzip2 1.01 compatibility

Revision 1.3  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.2  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

