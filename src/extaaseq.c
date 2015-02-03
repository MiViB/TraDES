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
#include <mmdbtraj.h>


/* returns 1 if an amino acid, 0 otherwise */
Int2 isAA(Char ch)
{
        Char chlow[2];
        
        chlow[0]=tolower(ch);
        chlow[1]=0;
        if (StringPBrk(chlow,"bjouxz")!=NULL) return 0;
        if (chlow[0]<'a' || chlow[0]>'z') return 0;
        return 1;
}

/* returns 1 if big endian, 0 if little endian and 2 if anything else (the
   latter is incorrectly treated as big endian in most code though! */
Int2 GetEndian(void)
{
  CharPtr TestPtr;
  Int4 TestVal;
  Int2 Endy;

  /* panic if not 4 bytes long!! */
  if (sizeof(Int4)!=4) {
    ErrPostEx(SEV_ERROR,99,97,"Unknown byte system on this machine");
    return ENDIAN_UNKNOWN;
  }
  TestPtr=(CharPtr)MemNew(sizeof(Char)*4);
  TestPtr[0]=0x01;
  TestPtr[1]=0x02;
  TestPtr[2]=0x03;
  TestPtr[3]=0x04;
  TestVal=*((Int4 *)TestPtr);
  Endy=ENDIAN_UNKNOWN;
  if (TestVal==0x04030201)
    Endy=ENDIAN_LITTLE;
  if (TestVal==0x01020304)
    Endy=ENDIAN_BIG;
  TestPtr=MemFree(TestPtr);
  if (Endy==ENDIAN_UNKNOWN) {
    ErrPostEx(SEV_ERROR,99,98,"Unknown Endian system on this machine");
  }
  return Endy;
}

/* allocate and returns a sequence obtained in the EXTAA_ENCODE mode
   converting to the equivalent of EXTAA_PARENT or EXTAA_X mode (see
   below) */
CharPtr DecodeSequence(CharPtr seq,Int2 extend)
{
	CharPtr dest;
	Int2 srcpos,destpos;

	dest=MemNew(sizeof(Char)*(StringLen(seq)+1));
	srcpos=0;
	destpos=0;
	while (seq[srcpos]!='\0') {
		if (seq[srcpos]=='*') {
			if (extend==EXTAA_PARENT)
				dest[destpos]=seq[srcpos+5];
			else /* extend==EXTAA_X */
				dest[destpos]='X';
			srcpos+=6;
		}
		else {
			dest[destpos]=seq[srcpos];
			srcpos++;
		}
		destpos++;
	}	
	/* null terminate string */
	dest[destpos]='\0';
	return dest;
}

PEAS GetExtAAInfoEnc(Int2 idx)
{
	FILE *f;
	Char inbuf[MAXCOL],fnam[PATH_MAX];
/*	Int2 numentries;*/
	PEAS peas;
	int idictidx,inidx0,inidx1,inidx2,icidx0,icidx1,icidx2;

	sprintf(fnam,"%s%s",CFG_local_datafilepath,EXTAA_FILENAME);
	if ((f=FileOpen(fnam,"r"))==NULL) {
		ErrPostEx(SEV_ERROR,1,2,"Unable to open AA translation table %s",fnam);
		return 0;
	}
	peas=(PEAS)MemNew(sizeof(EAS));
	while (FileGets(inbuf,MAXCOL,f)!=NULL) {
		if(inbuf[0] == '#') continue;
		sscanf(inbuf,"%*s %*s %*d %*s %c",&(peas->modlocation));
		if (peas->modlocation=='-') {
			sscanf(inbuf,"%s %s %d %4s %c %d %d %d %d %d %d",peas->keybname,peas->descr,&idictidx,peas->PDBName,&(peas->modlocation),&inidx0,&inidx1,&inidx2,&icidx0,&icidx1,&icidx2);
			peas->dictidx=(Int2)idictidx;
			peas->nidx[0]=(Int2)inidx0;
			peas->nidx[1]=(Int2)inidx1;
			peas->nidx[2]=(Int2)inidx2;
			peas->cidx[0]=(Int2)icidx0;
			peas->cidx[1]=(Int2)icidx1;
			peas->cidx[2]=(Int2)icidx2;
			/*numentries=3;*/
		}
		else {
			sscanf(inbuf,"%s %s %d %4s %c %d %d",peas->keybname,peas->descr,&idictidx,peas->PDBName,&(peas->modlocation),&inidx0,&icidx0);
			peas->dictidx=(Int2)idictidx;
			peas->nidx[0]=(Int2)inidx0;
			peas->cidx[0]=(Int2)icidx0;
			/*numentries=1;*/
		}
		/* compare ignoring case */
		if (peas->dictidx==idx) {
			/* we found the match */
			FileClose(f);
			return peas;
		}
		if ((peas->modlocation=='-') && (idx-peas->dictidx<3) && (idx-peas->dictidx>0)) {
			/* might be on ends */
			FileClose(f);
			return peas;
		}
	}
	FileClose(f);
	/* no match */
	peas=MemFree(peas);
	return NULL;
}

PEAS GetExtAAInfo(CharPtr buf)
{
	FILE *f;
	Char inbuf[MAXCOL],fnam[PATH_MAX];
/*	Int2 numentries;*/
	PEAS peas;
	int idictidx,inidx0,inidx1,inidx2,icidx0,icidx1,icidx2;

	sprintf(fnam,"%s%s",CFG_local_datafilepath,EXTAA_FILENAME);
	if ((f=FileOpen(fnam,"r"))==NULL) {
		ErrPostEx(SEV_ERROR,1,2,"Unable to open AA translation table %s",fnam);
		return 0;
	}
	peas=(PEAS)MemNew(sizeof(EAS));
	while (FileGets(inbuf,MAXCOL,f)!=NULL) {
			if(inbuf[0] == '#') continue;
			sscanf(inbuf,"%*s %*s %*d %*s %c",&(peas->modlocation));
			if (peas->modlocation=='-') {
				sscanf(inbuf,"%s %s %d %4s %c %d %d %d %d %d %d",peas->keybname,peas->descr,&idictidx,peas->PDBName,&(peas->modlocation),&inidx0,&inidx1,&inidx2,&icidx0,&icidx1,&icidx2);
				peas->dictidx=(Int2)idictidx;
				peas->nidx[0]=(Int2)inidx0;
				peas->nidx[1]=(Int2)inidx1;
				peas->nidx[2]=(Int2)inidx2;
				peas->cidx[0]=(Int2)icidx0;
				peas->cidx[1]=(Int2)icidx1;
				peas->cidx[2]=(Int2)icidx2;
/*				numentries=3;*/
			}
			else {
				sscanf(inbuf,"%s %s %d %4s %c %d %d",peas->keybname,peas->descr, &idictidx,peas->PDBName,&(peas->modlocation),&inidx0,&icidx0);
				peas->dictidx=(Int2)idictidx;
				peas->nidx[0]=(Int2)inidx0;
				peas->cidx[0]=(Int2)icidx0;
/*				numentries=1;*/
			}
			/* compare ignoring case */
			if (!StringICmp(buf,peas->keybname)) {
				/* we found the match */
				FileClose(f);
				return peas;
			}
		
	}
	FileClose(f);
	/* no match */
	peas=MemFree(peas);
	return NULL;
}

/* returns TRUE if PMAD is a S in a disulphide bridge */
PMAD IsDisulfideS(PMAD pmadHere)
{
        PMBD pmbdHere;
        PMAD pmad;
        ValNodePtr vnpHere;

        if (StringCmp(pmadHere->pcAName," SG "))
        	return FALSE;
        vnpHere=pmadHere->pvnBonds;
        while (vnpHere!=NULL) {
	        do {
	        	if (vnpHere==NULL)
	        		return FALSE;
  	      	pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
  	      	vnpHere=vnpHere->next;
  	      } while (((pmbdHere->bWhat) & BOND_SINGLE)==0);
	        if (pmbdHere->pmadFrom==pmadHere)
  	      	pmad=pmbdHere->pmadTo;
  	      else
  	      	pmad=pmbdHere->pmadFrom;
  	      if (!StringCmp(pmad->pcAName," SG "))
  	      	return pmad;
        }
        return NULL;
}

/* tells us if N-terminus is modified or standard (post-translational
modifications */
Boolean IsNTermModified(PMGD pmgdThis)
{
	PEAS peas;
	Boolean retval;

	if (pmgdThis->pcIUPAC[0]!='X')
		return FALSE;
	/* check for hydroxy-proline */
	if (pmgdThis->iIDict==212)
		return FALSE;
	/* check for pyro-glutamate */
	if (pmgdThis->iIDict==326)
		return FALSE;
	peas=GetExtAAInfoEnc(pmgdThis->iIDict);
	if (peas->modlocation=='N')
		retval=TRUE;
	else
		retval=FALSE;
	peas=MemFree(peas);
	return retval;
}

Boolean IsCTermModified(PMGD pmgdThis)
{
	PEAS peas;
	Boolean retval;

	if (pmgdThis->pcIUPAC[0]!='X')
		return FALSE;
	/* check for hydroxy-proline */
	if (pmgdThis->iIDict==211)
		return FALSE;
	peas=GetExtAAInfoEnc(pmgdThis->iIDict);
	if (peas->modlocation=='C')
		retval=TRUE;
	else
		retval=FALSE;
	peas=MemFree(peas);
	return retval;
}

/* convert an entered sequence in src into a one-letter-coded
amino acid string in dest.  There are 3 possible values for extend.

extend == EXTAA_X		replace any modified amino acids
				with X in the returned sequence

extend == EXTAA_PARENT		replace any modified amino acids
				with their unmodified parent in
				the returned sequence

extend == EXTAA_ENCODE		Any non-standard amino acids are
				represented by a *, a 4-digit integer
				(dictionary index) and then the unmodified
				amino acid 1-letter code 
*/
/* proc ensures N- and C-terminal modification are only entered at
the proper ends of the sequence (returns 0 if not) */
Int2 ConvertExtSequence(CharPtr src,CharPtr dest,Int2 extend)
{
	Int2 numAA;
	CharPtr destpos,pos,endpos;
	Char extbuf[255];
	PEAS peas;
	Int2 extcode;
	Boolean lastAA=FALSE;

	numAA=0;
	pos=src;
	destpos=dest;
	while (pos[0]!='\0') {
		if (pos[0]=='[') {
			/* special amino acid */
			endpos=StringChr(pos,']');
			/* check for close bracket - error if none */
			if (endpos==NULL)
				return 0;
			StringNCpy(extbuf,pos,(size_t)(endpos-pos+1));
			extbuf[endpos-pos+1]='\0';
			/* now extbuf holds extended aa code */
			pos=endpos;
			if (StringLen(extbuf)>MAX_KBAA_NAMELENGTH) {
				/* error, not in translation table */
				dest[0]='\0';
				return 0;
			}
			else if (extend==EXTAA_ENCODE) {
				/* get dictionary index of extended AA */
				peas=GetExtAAInfo(extbuf);
				if (peas==NULL) {
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				if ((peas->modlocation=='N' && numAA!=0) || lastAA) {
					peas=MemFree(peas);
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				if (peas->modlocation=='C') {
					/* cant add more now */
					lastAA=TRUE;
				}
				extcode=peas->dictidx;
				if (peas->modlocation=='-') {
					if (pos[1]=='\0')
						/* end of sequence */
						extcode++;
					if (numAA==0)
						extcode+=2;
				}
				destpos[0]='*';
				destpos[1]=extcode/1000+'0';
				destpos[2]=(extcode%1000)/100+'0';
				destpos[3]=(extcode%100)/10+'0';
				destpos[4]=(extcode%10)+'0';
				/* actual parent AA */
				destpos[5]=extbuf[1];
				destpos+=5;
				/* free up when done with this */
				MemFree(peas);
			}
			else if (extend==EXTAA_PARENT) {
				/* check if in table */
				peas=GetExtAAInfo(extbuf);
				if (peas==NULL) {
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				if ((peas->modlocation=='N' && numAA!=0) || lastAA) {
					peas=MemFree(peas);
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				if (peas->modlocation=='C') {
					/* cant add more now */
					lastAA=TRUE;
				}
				MemFree(peas);
				destpos[0]=toupper(extbuf[1]);
			}
			else { /* extend==EXTAA_X */
				/* check if in table */
				peas=GetExtAAInfo(extbuf);
				if (peas==NULL) {
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				if ((peas->modlocation=='N' && numAA!=0) || lastAA) {
					peas=MemFree(peas);
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				if (peas->modlocation=='C') {
					/* cant add more now */
					lastAA=TRUE;
				}
				MemFree(peas);
				destpos[0]='X';
			}
      numAA++;
			destpos++;
		}
		else {
			/* normal amino acid or junk */
      if (isAA(pos[0])) {
				if (lastAA) {
					/* error, return null */
					dest[0]='\0';
					return 0;
				}
				destpos[0]=toupper(pos[0]);
        numAA++;
				destpos++;
			}
		}
		pos++;
 	}
	destpos[0]='\0';
	return numAA;
}

Char GetAAFromIDict(PMGD pmgd)
{
	Char retval;
	PEAS peas;

	/* if local dictionary, still X */
	if ((pmgd->bWhat)&(Byte)DICT_LOCAL)
		return 'X';
	if (pmgd->pcIUPAC[0]=='X') {
		peas=GetExtAAInfoEnc(pmgd->iIDict);
		retval=peas->keybname[1];
		peas=MemFree(peas);
		/* pyro glutamate is proline */
		if (pmgd->iIDict==326)
			return 'P';
		return retval;	
	}
	return pmgd->pcIUPAC[0];
}

Int2 GetResIdxFromMGD(CharPtr AAlist,PMGD pmgdHere)
{
	Char resname[2];

	/* if local dictionary, still X */
	if ((pmgdHere->bWhat)&(Byte)DICT_LOCAL)
		return 0;
	if ((pmgdHere->pcIUPAC[0])=='X') {
		/* use iIDict to get parent amino acid */
		resname[0]=GetAAFromIDict(pmgdHere);
		resname[1]='\0';
		/* pyro-glutamate is like proline */
		if (pmgdHere->iIDict==326)
			return (Int2)StringCSpn(AAlist,"P");
		return (Int2)StringCSpn(AAlist,resname);
	}
	return (Int2)StringCSpn(AAlist,pmgdHere->pcIUPAC);
}

/* changes residue resnum in the sequence seq, which must be in EXTAA_ENCODE
   format.  What happens depend on the value of dowhat:

   EXTAASEQ_REPLACE: replaces the residue with newaa, which must be in keyboard
   entry format (see ConvertExtSequence)

   EXTAASEQ_INSERT: inserts the residue BEFORE residue resnum, so for example
   if resnum==2, newaa will be the new 2nd residue

   EXTAASEQ_DELETE: remove the residue (here newaa is ignored of course)

   returns the changed sequence (a piece of allocated memory)
*/
CharPtr AlterResidueInSequence(CharPtr seq,int resnum,CharPtr newaa,Int2 dowhat)
{
	CharPtr tmpseq,dest,pctmp;
	Int2 cnt,srcidx,destidx,tmpptr,length,extcode;
	PEAS peas,peas1;
	Char theAA[10];

	if (resnum<=0 || seq==NULL || newaa==NULL)
		return NULL;
	/* get seqlength */
	dest=(CharPtr)MemNew(sizeof(Char)*(StringLen(seq)+1));
	length=ConvertExtSequence(seq,dest,EXTAA_X);
	dest=MemFree(dest);
	/* convert newaa to proper code */
	if (StringLen(newaa)==1 && !isAA(newaa[0]))
		return NULL;
	if (StringLen(newaa)>1) {
		peas=GetExtAAInfo(newaa);
		if (peas==NULL) {
			/* no match, error, return null */	
			return NULL;
		}
		if (peas->modlocation=='N' && resnum!=1) {
			peas=MemFree(peas);
			/* no match, error, return null */	
			return NULL;
		}
		if (peas->modlocation=='C' && (resnum!=length || dowhat==EXTAASEQ_INSERT)) {
			peas=MemFree(peas);
			/* no match, error, return null */	
			return NULL;
		}
		extcode=peas->dictidx;
		if (peas->modlocation=='-' && resnum==length && dowhat!=EXTAASEQ_INSERT)
			extcode++;
		if (peas->modlocation=='-' && resnum==1)
			extcode+=2;
		theAA[0]='*';
		theAA[1]=extcode/1000+'0';
		theAA[2]=(extcode%1000)/100+'0';
		theAA[3]=(extcode%100)/10+'0';
		theAA[4]=(extcode%10)+'0';
		/* actual parent AA */
		theAA[5]=newaa[1];
		theAA[6]='\0';
		/* free up when done with this */
		MemFree(peas);
	}
	else
		StringCpy(theAA,newaa);
	tmpseq=(CharPtr)MemNew(sizeof(Char)*(StringLen(seq)+StringLen(theAA)+1));
	cnt=0;
	srcidx=0;
	destidx=0;
	while (cnt+1<resnum) {
		if (!seq[srcidx]) 
			break;
		if (isAA(seq[srcidx]))
			cnt++;
		tmpseq[destidx]=seq[srcidx];
		srcidx++;
		destidx++;
	}
	if (seq[srcidx]=='\0') {
		/* resnum is too high */
		ErrPostEx(SEV_ERROR,1,23,"Less than %d residues found in %s",resnum,seq);
		tmpseq=MemFree(tmpseq);
		return NULL;
	}
	if (dowhat==EXTAASEQ_DELETE) {
		if (resnum==1) {
			if (seq[0]=='*')
				pctmp=&seq[6];
			else
				pctmp=&seq[1];
			if (pctmp[0]=='*') {
				extcode=1000*(pctmp[1]-'0')+100*(pctmp[2]-'0')+10*(pctmp[3]-'0')+(pctmp[4]-'0');
				extcode+=2;			
				pctmp[1]=extcode/1000+'0';
				pctmp[2]=(extcode%1000)/100+'0';
				pctmp[3]=(extcode%100)/10+'0';
				pctmp[4]=(extcode%10)+'0';
			}
		}
		pctmp=&seq[srcidx];
		if (pctmp[0]=='*')
			pctmp+=6;
		else
			pctmp++;
		if (pctmp[0]=='\0') {
			/* delete last residue */
			pctmp=&tmpseq[destidx-6];
			/* check previous residue - if modified AA must correct iDict */
			if (pctmp[0]=='*') {
				extcode=1000*(pctmp[1]-'0')+100*(pctmp[2]-'0')+10*(pctmp[3]-'0')+(pctmp[4]-'0');
				extcode++;			
				pctmp[1]=extcode/1000+'0';
				pctmp[2]=(extcode%1000)/100+'0';
				pctmp[3]=(extcode%100)/10+'0';
				pctmp[4]=(extcode%10)+'0';
			}
		}	
	}
	if (resnum==1 && seq[0]=='*' && dowhat==EXTAASEQ_INSERT) {
		/* if inserting at beginning of sequence and special AA there... */
		extcode=1000*(seq[1]-'0')+100*(seq[2]-'0')+10*(seq[3]-'0')+(seq[4]-'0');		
		peas1=GetExtAAInfoEnc(extcode);
		/* check if will insert before an N-only residue */
		if (peas1->modlocation=='N') {
			ErrPostEx(SEV_ERROR,1,24,"Cannot insert residue - current first amino acid only ever found at N-terminus!");
			tmpseq=MemFree(tmpseq);
			return NULL;
		}
		/* else modlocation=='-' */
		/* so change from N-terminal form to "standard" form */
		extcode-=2;
		seq[1]=extcode/1000+'0';
		seq[2]=(extcode%1000)/100+'0';
		seq[3]=(extcode%100)/10+'0';
		seq[4]=(extcode%10)+'0';
		peas1=MemFree(peas1);
	}
	/* idx points to start of residue resnum now */
	switch (dowhat) {
		case EXTAASEQ_REPLACE:
			/* skip over current residue */
			while (!isAA(seq[srcidx])) 
				srcidx++;
			srcidx++;
		/* replace is like insert except for extra code above */
		case EXTAASEQ_INSERT:
			tmpptr=0;
			while (theAA[tmpptr]!='\0') {
				tmpseq[destidx]=theAA[tmpptr];
				destidx++;
				tmpptr++;
			}
			break;	
		case EXTAASEQ_DELETE:
			/* skip over current residue */
			while (!isAA(seq[srcidx])) 
				srcidx++;
			srcidx++;
			break;	
		default:;
	}
	while (seq[srcidx]!='\0') {
		tmpseq[destidx]=seq[srcidx];
		srcidx++;
		destidx++;
	}
	return tmpseq;
}

BiostrucPtr MIMEBiostrucAsnGet(CharPtr fnam,const CharPtr mode,ValNodePtr *BioseqKeep)
{
  AsnIoPtr aip=NULL;
  NcbiMimeAsn1Ptr nmap=NULL;
  BiostrucSeqPtr bssp=NULL;
  BiostrucPtr bsp=NULL;
  ErrSev esMsg,esLog,esFatal;
  if (BioseqKeep!=NULL)
	*BioseqKeep=NULL;



  aip=AsnIoOpen(fnam,mode);
  if (aip==NULL)
		return NULL;



/* first try biostruc */

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
	if (bsp!=NULL)
		return bsp;
/* then try NCBIMime */
	aip=NULL;
	aip=AsnIoOpen(fnam,mode);
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
		    ErrPostEx(SEV_ERROR,6,9,"Apparently wrong MIME type set in Asn.1 wrapper... Convert to ASN.1 text and fix it.");
			return NULL;
		}

		/* must be v2.0 file */
		bssp=(BiostrucSeqPtr)(nmap->data.ptrvalue);
		bsp=bssp->structure;

		if (BioseqKeep!=NULL) {
			*BioseqKeep=bssp->sequences;
			bssp->sequences=NULL;
/* caller must delete the sequences with SeqEntryFree() 
 whenever LoadABiostrucEx is called directly ! */
		}

		/* free the MIME wrapper now */
		bssp->structure=NULL;
		nmap=NcbiMimeAsn1Free(nmap);
		return bsp;
	}
	return NULL;
}



NcbiMimeAsn1Ptr BuildMIMEBiostruc(BiostrucPtr bsp,CharPtr seq,ValNodePtr sequences)
{
	NcbiMimeAsn1Ptr nmap;
	BiostrucSeqPtr bsqp;
	ValNodePtr seqlist,seqent;
	BioseqSetPtr bssp;
	SeqEntryPtr sep;
	BioseqPtr bp;
	SeqIdPtr sip;
	ValNodePtr vnpDescr,vnpTemp,SeqDescr,vnpLoc;
	PDBSeqIdPtr psip;
	ByteStorePtr pbs;
	NumberingPtr np;
	NumEnumPtr nep;
	Char tmpbuf[MAXCOL];
	Int4 cnt,lencnt;
	CharPtr pidx;
	SeqAnnotPtr sap;
	BiostrucFeatureSetPtr pbsfs;
	BiostrucFeaturePtr pbsf;
	ResidueIntervalPntrPtr pri;
	SeqIntPtr psi;
	SeqFeatPtr sfp,sfpHere,sfpHead;
	PRGD prgdDictionary = NULL;



		
	if (sequences==NULL) {
		nep=NumEnumNew();
		nep->num=StringLen(seq);
		nep->names=(CharPtr *)MemNew((nep->num)*sizeof(CharPtr));
		lencnt=0;
		for (cnt=0;cnt<nep->num;cnt++) {
			sprintf(tmpbuf,"%d",cnt+1);
			lencnt+=StringLen(tmpbuf)+1;
		}
		nep->buf=(CharPtr)MemNew(sizeof(Char)*lencnt);
		pidx=nep->buf;
		for (cnt=0;cnt<nep->num;cnt++) {
			sprintf(tmpbuf,"%d",cnt+1);
			StringMove(pidx,tmpbuf);
			(nep->names)[cnt]=pidx;
			pidx=pidx+StringLen(tmpbuf)+1;
		}
		np=(NumberingPtr)ValNodeNew(NULL);
		np->choice=Numbering_enum;
		np->data.ptrvalue=(VoidPtr)nep;
		vnpDescr=ValNodeNew(NULL);
		vnpDescr->choice=Seq_descr_num;
		vnpDescr->data.ptrvalue=(VoidPtr)np;
		psip=PDBSeqIdNew();
		vnpTemp=bsp->descr;
		while (vnpTemp->choice!=BiostrucDescr_name) {
			if (vnpTemp==NULL)
				break;
			vnpTemp=vnpTemp->next;
		}
		if (vnpTemp!=NULL)
			psip->mol=StringSave((CharPtr)(vnpTemp->data.ptrvalue));
		else
			psip->mol=StringSave("RAND");
		/* chain doesn't make sense anymore, always one chain */
		psip->chain=(Uint1)' ';
		psip->rel=DateCurr();
		sip=(SeqIdPtr)ValNodeNew(NULL);
		sip->choice=SEQID_PDB;
		sip->data.ptrvalue=(VoidPtr)psip;
		pbs=BSNew(1);
		BSWrite(pbs,seq,StringLen(seq));
		/* assign seq-id to molecule-graphs in biostruc next */
		if (bsp->chemical_graph->molecule_graphs->seq_id!=NULL)
			AsnGenericChoiceSeqOfFree((Pointer) bsp->chemical_graph->molecule_graphs->seq_id, (AsnOptFreeFunc) SeqIdFree);
		bsp->chemical_graph->molecule_graphs->seq_id=(SeqIdPtr)AsnIoMemCopy((Pointer)sip,(AsnReadFunc)SeqIdAsnRead,(AsnWriteFunc)SeqIdAsnWrite);
		/* add secondary structure as seqannots */
		pbsfs=bsp->features;
		while (pbsfs!=NULL) {
			if (pbsfs->descr->choice==BiostrucFeatureSetDescr_name && !StringCmp((CharPtr)(pbsfs->descr->data.ptrvalue),STRING_NCBI_SS))
				break;
			pbsfs=pbsfs->next;
		}
		sap=NULL;
		if (pbsfs!=NULL) {
			sap=SeqAnnotNew();
			sap->type=1; /* ftable */
			pbsf=pbsfs->features;
			sfpHead=NULL;
			while (pbsf!=NULL) {
				/* add one feature at a time */
				psi=SeqIntNew();
				pri=(ResidueIntervalPntrPtr)(((ValNodePtr)(((ValNodePtr)(pbsf->Location_location->data.ptrvalue))->data.ptrvalue))->data.ptrvalue);
				psi->from=pri->from-1;
				psi->to=pri->to-1;
				psi->id=SeqIdDup(sip); /* seqidptr */
				vnpLoc=ValNodeNew(NULL);
				vnpLoc->choice=SEQLOC_INT;
				vnpLoc->data.ptrvalue=(VoidPtr)psi;
				sfp=SeqFeatNew();
				sfp->data.choice=SEQFEAT_PSEC_STR;
				if (pbsf->type==Feature_type_helix)
					sfp->data.value.intvalue=1; /* 1 for helix, 2 for sheet */
				else if (pbsf->type==Feature_type_strand)
					sfp->data.value.intvalue=2;
				sfp->comment=StringSave(pbsf->name); /* helix 1, etc. */
				sfp->location=vnpLoc;
				if (sfpHead==NULL)
					sfpHead=sfp;
				else {
					sfpHere=sfpHead;
					while (sfpHere->next!=NULL)
						sfpHere=sfpHere->next;
					sfpHere->next=sfp;
				}
				pbsf=pbsf->next;
			}
			sap->data=(Pointer)sfpHead;
		}
		bp=BioseqNew();
		bp->id=sip;
		bp->descr=vnpDescr;
		bp->repr=Seq_repr_raw;
		bp->mol=Seq_mol_aa;
		bp->length=StringLen(seq);
		bp->seq_data=pbs;
		bp->seq_data_type=(Uint1)Seq_code_ncbieaa;
		bp->annot=sap;
		sep=SeqEntryNew();
		sep->choice=1; /* Bioseq */
		sep->data.ptrvalue=(VoidPtr)bp;
		bssp=BioseqSetNew();
		bssp->_class=BioseqseqSet_class_pdb_entry;
		bssp->seq_set=sep;
		seqlist=ValNodeNew(NULL);
		seqlist->choice=2; /* Bioseq-set */
		seqlist->data.ptrvalue=(VoidPtr)bssp;
	}
	else {
		seqlist=sequences;
		seqent=seqlist;
		while(IS_Bioseq_set(seqent)) {
			bssp=(BioseqSetPtr)(seqent->data.ptrvalue);
			seqent=bssp->seq_set;
		}
		if (IS_Bioseq(seqent)) {
			bp=(BioseqPtr)(seqent->data.ptrvalue);
			bp->length=StringLen(seq);
			/* remove secondary structure annontation - no longer applicable, remake it */
			bp->annot=SeqAnnotFree(bp->annot);
			pbsfs=bsp->features;
			while (pbsfs!=NULL) {
				if (pbsfs->descr->choice==BiostrucFeatureSetDescr_name && !StringCmp((CharPtr)(pbsfs->descr->data.ptrvalue),STRING_NCBI_SS))
					break;
				pbsfs=pbsfs->next;
			}
			sap=NULL;
			if (pbsfs!=NULL) {
				sap=SeqAnnotNew();
				sap->type=1; /* ftable */
				pbsf=pbsfs->features;
				sfpHead=NULL;
				while (pbsf!=NULL) {
					/* add one feature at a time */
					psi=SeqIntNew();
					pri=(ResidueIntervalPntrPtr)(((ValNodePtr)(((ValNodePtr)(pbsf->Location_location->data.ptrvalue))->data.ptrvalue))->data.ptrvalue);
					psi->from=pri->from-1;
					psi->to=pri->to-1;
					psi->id=SeqIdDup(bp->id); /* seqidptr */
					vnpLoc=ValNodeNew(NULL);
					vnpLoc->choice=SEQLOC_INT;
					vnpLoc->data.ptrvalue=(VoidPtr)psi;
					sfp=SeqFeatNew();
					sfp->data.choice=SEQFEAT_PSEC_STR;
					if (pbsf->type==Feature_type_helix)
						sfp->data.value.intvalue=1; /* 1 for helix, 2 for sheet */
					else if (pbsf->type==Feature_type_strand)
						sfp->data.value.intvalue=2;
					sfp->comment=StringSave(pbsf->name); /* helix 1, etc. */
					sfp->location=vnpLoc;
					if (sfpHead==NULL)
						sfpHead=sfp;
					else {
						sfpHere=sfpHead;
						while (sfpHere->next!=NULL)
							sfpHere=sfpHere->next;
						sfpHere->next=sfp;
					}
					pbsf=pbsf->next;
				}
				sap->data=(Pointer)sfpHead;
			}
			bp->annot=sap;
			pbs=BSNew(1);
			BSWrite(pbs,seq,StringLen(seq));
			BSFree(bp->seq_data);
			bp->seq_data=pbs;
			/* assign seq-id to molecule-graphs in biostruc next - necessary so Cn3D highlights properly */
			if (bsp->chemical_graph->molecule_graphs->seq_id!=NULL)
				AsnGenericChoiceSeqOfFree((Pointer) bsp->chemical_graph->molecule_graphs->seq_id, (AsnOptFreeFunc) SeqIdFree);
			bsp->chemical_graph->molecule_graphs->seq_id=(SeqIdPtr)AsnIoMemCopy((Pointer)(bp->id),(AsnReadFunc)SeqIdAsnRead,(AsnWriteFunc)SeqIdAsnWrite);
			SeqDescr=bp->descr;
			while (SeqDescr->choice!=Seq_descr_num && SeqDescr->next!=NULL)
				SeqDescr=SeqDescr->next;
			if (SeqDescr->choice==Seq_descr_num) {
				np=(NumberingPtr)(SeqDescr->data.ptrvalue);
				nep=(NumEnumPtr)(np->data.ptrvalue);
				nep->num=bp->length;
				MemFree(nep->names);
				MemFree(nep->buf);
				nep->names=(CharPtr *)MemNew((nep->num)*sizeof(CharPtr));
				lencnt=0;
				for (cnt=0;cnt<nep->num;cnt++) {
					sprintf(tmpbuf,"%d",cnt+1);
					lencnt+=StringLen(tmpbuf)+1;
				}
				nep->buf=(CharPtr)MemNew(sizeof(Char)*lencnt);
				pidx=nep->buf;
				for (cnt=0;cnt<nep->num;cnt++) {
					sprintf(tmpbuf,"%d",cnt+1);
					StringMove(pidx,tmpbuf);
					(nep->names)[cnt]=pidx;
					pidx=pidx+StringLen(tmpbuf)+1;
				}
			}
		}
	}			
	bsqp=BiostrucSeqNew();
	bsqp->structure=bsp;
	/* free features */
/*	AsnGenericChoiceSeqOfFree(bsp->features,(AsnOptFreeFunc)BiostrucFeatureSetFree);
	bsp->features=NULL;
*/
	bsqp->sequences=seqlist;
	nmap=(NcbiMimeAsn1Ptr)ValNodeNew(NULL);
	nmap->choice=NcbiMimeAsn1_strucseq;
	nmap->next=NULL;
	nmap->data.ptrvalue=(VoidPtr)bsqp;
	return nmap;
}
Int2 GetFirstFullModel(PMSD pmsdHere)
{
	PDNML pdnmlModel;
	PMLD pmldThis;
	
	pdnmlModel=pmsdHere->pdnmlModels;
	while (pdnmlModel) {
		pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
		if (pmldThis->iType==Model_type_ncbi_all_atom || pmldThis->iType==Model_type_pdb_model) {
			return pdnmlModel->choice;
		}
		pdnmlModel=pdnmlModel->next;
	}
	/* no good model found */
	return 0;
}

PMLD GetModelN(PMSD pmsdHere,Int2 n)
{
	PDNML pdnmlModel;
	
	pdnmlModel=pmsdHere->pdnmlModels;
	while (pdnmlModel) {
		if (pdnmlModel->choice==n)
			return (PMLD)(pdnmlModel->data.ptrvalue);
		pdnmlModel=pdnmlModel->next;
	}
	/* no good model found */
	return NULL;
}

PMSD LoadABiostrucEx(CharPtr fnam,Boolean remote,Int2 mdllvl,Int2 *reqmodel,ValNodePtr *ppvnBioseq)
{
	BiostrucPtr bsp;
	PDNMS pdnmsModelstruc;
	PMSD pmsdRoot;
	Int2 modelnum;
	PMLD pmldHere;

	bsp=NULL;	



	if (remote == FALSE) {

		bsp=MIMEBiostrucAsnGet(fnam,"rb",ppvnBioseq);
		if (bsp==NULL)
			bsp=MIMEBiostrucAsnGet(fnam,"r",ppvnBioseq);
	}
	else
		bsp=FetchBiostrucPDB(fnam,mdllvl,100);
	if (bsp==NULL) {
		ErrPostEx(SEV_ERROR,3,1,"Unable to fetch Biostruc %s",fnam);
		return NULL;
	}
	pdnmsModelstruc=MakeAModelstruc(bsp);
	if (pdnmsModelstruc==NULL) {
		ErrPostEx(SEV_ERROR,4,1,"Unable to convert Biostruc to Modelstruc");
		return NULL;
	}
 	pmsdRoot=(PMSD)(pdnmsModelstruc->data.ptrvalue);
	if (pmsdRoot==NULL) {
		ErrPostEx(SEV_ERROR,5,1,"Unable to determine number of models");
		return NULL;
	}
	/* get number of models */
	modelnum=GetFirstFullModel(pmsdRoot);
	if (!modelnum) {
		ErrPostEx(SEV_ERROR,5,3,"No all-atom models found");
		return NULL;
	}
	pmldHere=GetModelN(pmsdRoot,*reqmodel);
	if (pmldHere==NULL){
		ErrPostEx(SEV_ERROR,6,1,"No model %d found in structure!!",*reqmodel);
		return NULL;
	}
	if (pmldHere->iType!=Model_type_ncbi_all_atom && pmldHere->iType!=Model_type_pdb_model)
		ErrPostEx(SEV_ERROR,5,2,"Using model %ld instead of %ld since the latter is not all-atom",(long)modelnum,(long)(*reqmodel));
	else
		modelnum=*reqmodel;
	*reqmodel=modelnum;
	return pmsdRoot;
}

PMSD LoadABiostruc(CharPtr fnam,Boolean remote,Int2 mdllvl,Int2 *reqmodel)
{


	return LoadABiostrucEx(fnam,remote,mdllvl,reqmodel,NULL);

}

/*  
$Log: extaaseq.c,v $
Revision 1.18  2003/09/14 02:38:23  feldman
Fixed unused variables and other minor compiler warnings

Revision 1.17  2003/03/07 20:00:14  feldman
Add SS features when building MIME Biostruc

Revision 1.16  2003/01/24 16:52:25  feldman
Improved error msg

Revision 1.15  2003/01/08 16:52:33  feldman
maketrj now reads both ascii and binary val input files

Revision 1.14  2002/09/26 13:23:22  michel
Moved BuildMIMEBiostruc to mmdbtrajlib

Revision 1.13  2002/07/05 21:19:00  feldman
Made LoadABiostruc handle ASCII and binary both

Revision 1.12  2001/12/03 15:52:42  feldman
Change IsBondDisulfide return value

Revision 1.11  2001/09/24 14:51:10  feldman
Added IsDisulfide function

Revision 1.10  2001/06/18 17:34:05  feldman
Added datafile path to extaa file and amber parameter files

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

Revision 1.7  2001/02/08 19:46:36  feldman
Added code to improve SubstituteNames and make it do a lot
for checking of atomnames when entering constraints

Revision 1.6  2001/02/06 18:40:39  feldman
Added a few functions for dealing with distance constraints and
tidied up so maketrj could join the library (foldtrajlib) without
conflicts

Revision 1.5  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.4  2000/10/25 15:15:13  feldman
Made further updates, multiple model support is correct now
and relocated to a single function for loading correct model
and extensive error handling

Revision 1.3  2000/10/24 20:57:25  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.2  2000/07/14 20:09:25  feldman
Added parameter for MIMEBiostrucAsnGet for Bioseq

Revision 1.1  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now


*/

