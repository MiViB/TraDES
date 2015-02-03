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


#include <foldtrajlib.h>
#include <slri_misc.h>

#define RLEMAX 0x20000000UL
/* maximum number of degenerate atoms in a single NOE */
#define MAX_DEGEN 9		
/* times to try to satisfy an NOE before giving up on it */
#define CONSTR_TRIES_MAX 5000	
#define CONSTR_ANGLE_TOL 45.0
#define CONSTR_DIHED_TOL 45.0

static PNN pnnDistConstHead=NULL;
static PNN pnnDistConstTrueHead=NULL;
static Boolean mmdbclosed=TRUE;
static Boolean LogDumped=FALSE;
static vec vZero={0,0,0};
static Boolean largeram_loaded=FALSE;
static PTGS *ptgsarr=NULL;
static Int2 ptgsarrsize=0;

/* return TRUE if checksum matches file checksum */
Boolean CheckMD5(CharPtr fnam,CharPtr cksum)
{
	FILE *fpbl;
 	Uint4 mdlen;
	UcharPtr mdbuf;
	Uchar mdoutbuf[32];
	Char md5out[33];
	MD5Context context; 
	Int2 cnt;
	Char pctmp[3];

  if ((fpbl=FileOpen(fnam,"rb"))==NULL) {
		ErrPostEx(SEV_ERROR,1,8,"CheckMD5: File not found %s",fnam);
		return FALSE;
  }
	mdlen=FileLength(fnam);
	mdbuf=MemNew(mdlen*sizeof(Char));
	FileRead(mdbuf,sizeof(Char),(size_t)mdlen,fpbl);
	FileClose(fpbl);
	MD5Init(&context);
	MD5Update(&context,(UcharPtr)mdbuf,(Uint4)mdlen);
	MD5Final(&context,(UcharPtr)mdoutbuf);
	mdbuf=MemFree(mdbuf);
	StringCpy(md5out,"");
	for (cnt=0;cnt<16;cnt++) {
		sprintf(pctmp,"%02x",mdoutbuf[cnt]);
		StringCat(md5out,pctmp);
	}
	if (!StringCmp(md5out,cksum))
		return TRUE;
	return FALSE;
}

void CleanUpDB(CharPtr fnam)
{
	Char fnamtmp[PATH_MAX];

        StringCpy(fnamtmp,fnam);
        StringCat(fnamtmp,DB_EXT1);
        FileRemove(fnamtmp);
        StringCpy(fnamtmp,fnam);
        StringCat(fnamtmp,DB_EXT2);
        FileRemove(fnamtmp);
        StringCpy(fnamtmp,fnam);
        StringCat(fnamtmp,DB_EXT3);
        FileRemove(fnamtmp);
}

/* bzip is TRUE ONLY for distributedfolding and will also now strip out all but the best structure from the file
   as well as other funky stuff */
void DumpLog(CharPtr logname,ByteStorePtr bspHere,Boolean bzipit,CharPtr extralines)
{
    Char ftmp[PATH_MAX];
    Char ftmp2[PATH_MAX];
    Char tmpfnam[PATH_MAX];
    Char tmpout[PATH_MAX];
    Char buf[PATH_MAX];
    Char bestbuf[PATH_MAX]="";
    Char numbuf[PATH_MAX];
    Char mdinbuf[PATH_MAX];
	Int4 numbytes;
	CharPtr pcTempBuf,pcHere/*,chptr*/;
/*	Int2 gen=0;*/
	FILE *fp,*fp2,*g;
	int err,iprotver,icacaext;
	long lstrucnum,ldssphelix,ldsspsheet;
	float frgyr,fhprgyr,frnc,fsurfacc,fhpsurfacc,frmsd,fzhang,fbryant,fcrease;
	char cterm;
	Int4 struccnt=0;
	FloatHi bestenergy=9999999.9;
	MD5Context context;
	Uchar mdoutbuf[32];
	Char md5out[33];
	Char pctmp[3];
	Int2 cnt;

	if (bspTempLog==NULL)
	    return;
    numbytes=BSLen(bspHere);
	if (numbytes==0)
	    return;
	
	/* print separator */
	BSprintf(bspHere,"----------------------------------------\n\n");

/*	sprintf(loglock,".lock_%s",logname);*/    
    /* wait until clear to add to log */
/*  do {} while (CheckLock(loglock)==FALSE);*/
    
    LogDumped=TRUE;
    /* lock log file */
/*  fp=FileOpen(loglock,"w");
    fprintf(fp,"%ld\n",GetAppProcessID());
    FileClose(fp);*/
    
    StringCpy(ftmp2,logname);
    StringCat(ftmp2,LOG_EXT);
    if (traj_quiet==VERBOSITY_VERBOSE)
            printf("\nWriting out log file %s\n",ftmp2);
	
	/* get numbytes again since we just added separator */
    numbytes=BSLen(bspHere);
    pcTempBuf=(CharPtr)MemNew((size_t)numbytes);
    BSSeek(bspHere,0L,SEEK_SET);
    BSRead(bspHere,pcTempBuf,numbytes);
	pcHere=&pcTempBuf[numbytes-44];
		
	while (pcHere>pcTempBuf && *pcHere!=LOG_TERM_CHAR)
		pcHere--;
	if (pcHere==pcTempBuf) {
		/* no complete record to dump in log so abort */
		return;
	}
	if (pcHere<&pcTempBuf[numbytes-44]) {
		pcHere+=2;
		sprintf(pcHere,"----------------------------------------\n\n");
		numbytes=pcHere-pcTempBuf+42;
	}
	
    if ((fp=FileOpen(ftmp2,"ab"))==NULL) {
        ErrPostEx(SEV_ERROR,1,8,"Unable to open log file %s",ftmp);
        return;
    }
    
    if ( (FileWrite(pcTempBuf,sizeof(Char),numbytes,fp)) != numbytes ){
        FileClose(fp);
        ErrPostEx(SEV_ERROR,1,9,"Unable to write to log file %s",ftmp);
        return;
    }
    FileClose(fp);
  
    pcTempBuf=MemFree(pcTempBuf);

	
}

Int2 OPENMMDBAPI(Byte bExtent,CharPtr pcDictName)
{
	Int2 res;

	/* check if already open */
	if (mmdbclosed==FALSE)
		return 1;
	res=OpenMMDBAPI(bExtent,pcDictName);
	if (res)
		mmdbclosed=FALSE;
	return res;
}

void CLOSEMMDBAPI(void)
{
	if (mmdbclosed==TRUE)
		return;
	CloseMMDBAPI();
	mmdbclosed=TRUE;
}

void PurgeGlobs(void) {
  FILE *f,*g;
  Char ftmp[PATH_MAX];

  if (tmpdbasename[0]) {
  	CleanUpDB(tmpdbasename);
	tmpdbasename[0]='\0';
  }
  if (tmpskelfname[0]) {
  	FileRemove(tmpskelfname);
	tmpskelfname[0]='\0';
  }
/*  if (loglock[0]) {
  	FileRemove(loglock);
	loglock[0]='\0';
  }
*/
	ProgramProgress=0;
	ProgramProgressMax=0;
    if (!LogDumped) {
  	    if (bspTempLog!=NULL)
		    BSprintf(bspTempLog,"\n");
	    if (StringCmp(LogOutName,""))
	  	    DumpLog(LogOutName,bspTempLog,FALSE,NULL);
	    if (fExtraLog!=NULL)
		    FileClose(fExtraLog);
	    if (StringCmp(ExtraOutName,"")) {
	        f=FileOpen(ExtraOutName,"rb");
	        StringCpy(ftmp,ExtraOutName);
	        StringCat(ftmp,BZ_EXT);
	        if ((g=FileOpen(ftmp,"wb"))!=NULL){
                /* note streams are closed in compressStream */
                compressStream(f,g);
            }
	        FileRemove(ExtraOutName);
	    }
        /*DumpLog(ExtraOutName,bspExtraLog,TRUE);*/
    }
    /*if (info!=NULL) {
	    WWWInfoFree(info);
	    info=NULL;
    }*/
    CLOSEMMDBAPI();
}

/* check residue dictionary to see if atom exists */
TrajErr ValidateAtomName(Char atomname[5],Int4 DictIdx)
{
	PRGD prgdDict;
	ResidueGraphPtr prg;
	AtomPtr atom;
	Boolean matched=FALSE;
	
	prgdDict=LoadDict(DICT_DEFAULT);
	prg=prgdDict->residue_graphs;
	while (prg!=NULL) {
		if (prg->id==DictIdx) {
			/* found right residue */
			atom=prg->atoms;
			while (atom!=NULL) {
				/* compare name (or IUPAC-code?) */
				if (WildCardMatch(atom->name,atomname)) {
					matched=TRUE;
					break;
				}			
				atom=atom->next;
			}
			break;
		}	
		prg=prg->next;
	}
	prgdDict=BiostrucResidueGraphSetFree(prgdDict);
	if (matched) return ERR_SUCCESS;
	return ERR_FAIL;
}

/* convert X-PLOR/CNS atom names to PDB atom names */
/* but hopefully do nothing if already converted... */
/* dictidx should have +1 or +2 added for N- and C- termini if needed */
TrajErr SubstituteNames(Char nam[5],Int4 DictIdx)
{
	Char buf[5]="    ";
	Int2 l;
	Char ctmp;
	Char resid;
	PEAS peas;

	if (!StringCmp(nam,"HN") || !StringCmp(nam,"H")) {
		StringCpy(nam," H  ");
		return ValidateAtomName(nam,DictIdx);
	}
	l=StringLen(nam);
	/* remove trailing bracket */
	if (nam[l-1]==')') {
		l--;
		nam[l]=0;
	}
	if (l==4 && (isdigit((int)nam[0]) || nam[0]==' '))
		/* assume named correctly already */
		return ValidateAtomName(nam,DictIdx);
	/* test if already in correct format */
	if (isdigit((int)nam[0])) {
		/* and if so mess it up to X-PLOR name so rest of function will
			 work as expected */
		ctmp=nam[0];
		nam[0]=nam[1];
		nam[1]=nam[2];
		nam[2]=nam[3];
		nam[3]=ctmp;
	}
	if (DictIdx<=60)
		/* normal dictionary */
		resid=aalist[(DictIdx-1)/3];
	else {
		/* extended dictionary */
		peas=GetExtAAInfoEnc(DictIdx);
		if (peas==NULL) {
			ErrPostEx(SEV_ERROR,1,1,"Invalid dictionary index %ld passed to SubstituteNames",(long)DictIdx);
			return ERR_FAIL;
		}
		resid=(peas->keybname)[1];
		peas=MemFree(peas);
	}
	/* start messing around with string here */
	buf[1]=nam[0];
	if (l>1)
		buf[2]=nam[1];
	if (l==3) {
		if (nam[0]=='H') {
			/* must figure out if 3rd digit refers to hydrogen number or number of
			   atom to which it is attached */
			if ((resid=='F' && (nam[1]=='D' || nam[1]=='E')) || 
			 (resid=='H' && (nam[1]=='D' || nam[1]=='E')) || 
			 (resid=='T' && nam[1]=='G' && nam[2]=='1') || 
			 (resid=='W' && (nam[1]=='D' || nam[1]=='E' || nam[1]=='Z' || nam[1]=='H')) || 
			 (resid=='Y' && (nam[1]=='D' || nam[1]=='E'))) 
				buf[3]=nam[2];
			else
				buf[0]=nam[2];
		}
		else
			buf[3]=nam[2];
	}
	else if (l==4) {
		buf[3]=nam[2];
		buf[0]=nam[3];
	}
	/* convert wildcard characters; * in X-PLOR mean char 0 and
	   char 3 are anything */
	if (buf[0]=='*') buf[3]='*';
	if (buf[3]=='*') buf[0]='*';
	if (buf[0]=='#') buf[0]='*';
	if (buf[3]=='#') buf[3]='*';
	StringCpy(nam,buf);
	return ValidateAtomName(nam,DictIdx);
}

/* includes the +1 or +2 for start and end residues */
Int4 GetResIdxFromSeq(CharPtr encseq, Int4 iresnum)
{
	Int4 pos,residx;
	Char pctmp[2];
	
	pos=0;
	while (pos<(iresnum-1)) {
		if (encseq[pos]=='*')
			pos+=5;
		pos++;
	}
	if (encseq[pos]=='*') {
		residx=1000*(encseq[pos+1]-'0')+100*(encseq[pos+2]-'0')+10*(encseq[pos+3]-'0')+(encseq[pos+4]-'0');
	}
	else {
		pctmp[0]=encseq[pos];
		pctmp[1]='\0';
		residx=StringCSpn(aalist,pctmp)*3+1;
		if (iresnum==1)
			/* N-terminus */
			residx+=2;
		if (encseq[pos+1]=='\0')
			/* C-terminus */
			residx++;	
	}
	return residx;
}

TrajErr LoadDistConstraints(CharPtr fnam,CharPtr kbsequence)
{
	FILE *f;
	Char buf[255];
	int ires1,ires2;
	Int2 pos,nameidx;
	float fmeandist,fmindelta,fmaxdelta,fprob,fangl,fang2,fdihed01,fdihed12,fdihed23;
	Char atomname1[5];
	Char atomname2[5];
	Char cmd[255];
	Int4 residx1,residx2;
	CharPtr curpos;
	PNN pnnHere;
	Char encseq[MAXSIZE];

	/* use sequence with X's for modified AAs */
	ConvertExtSequence(kbsequence,encseq,EXTAA_ENCODE);
	if ((f=FileOpen(fnam,"r"))==NULL) return ERR_FAIL;
 	while (FileGets(buf,255,f)!=NULL) {
		/* check for comments */
		if ((buf[0]=='!') || (sscanf(buf,"%s",cmd)==EOF)) continue;
		/* clear command */
		cmd[0]=0;
		if (sscanf(buf,"%s",cmd)!=1) goto readerr;
		if (!StringCmp(cmd,"assign")) {
     	pos=6;  /* length of "assign" */
     	while (buf[pos]!='(') {
     		if (!buf[pos]) {
     			/* read next line */
     			if (FileGets(buf,255,f)==NULL) goto readerr;
     			pos=0;
     		}
     		else
     			pos++;
     	}
     	pos++;
     	/* point to first char after open bracket */
     	curpos=&(buf[pos]);
			while (sscanf(curpos,"%s",cmd)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			if (StringCmp(cmd,"resid")) goto readerr;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%d",&ires1)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
			}
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%s",cmd)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			if (StringCmp(cmd,"and")) goto readerr;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%s",cmd)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			if (StringCmp(cmd,"name")) goto readerr;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%4s",atomname1)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			/* atom name plus possibly junk trailing, now in cmd */
			nameidx=0;
			/* removing trailing crap */
			while (atomname1[nameidx] && (!isspace((int)(atomname1[nameidx]))) && (atomname1[nameidx]!=')'))
				nameidx++;
			atomname1[nameidx]=0;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip atom name and bracket */
      while ((*curpos)&&((*curpos)!=')')) curpos++;
      /* on ')' or end of line */
      while (!(*curpos)) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
        while ((*curpos)&&((*curpos)!=')')) curpos++;
      }
      /* scan for '(' now */
      while ((*curpos)&&((*curpos)!='(')) curpos++;
      while (!(*curpos)) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
        while ((*curpos)&&((*curpos)!='(')) curpos++;
      }
      curpos++;
			while (sscanf(curpos,"%s",cmd)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			if (StringCmp(cmd,"resid")) goto readerr;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%d",&ires2)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%s",cmd)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			if (StringCmp(cmd,"and")) goto readerr;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%s",cmd)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			if (StringCmp(cmd,"name")) goto readerr;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%4s",atomname2)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
			/* atom name plus possibly junk trailing, now in cmd */
			nameidx=0;
			/* removing trailing crap */
			while (atomname2[nameidx] && (!isspace((int)(atomname2[nameidx]))) && (atomname2[nameidx]!=')'))
				nameidx++;
			atomname2[nameidx]=0;
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip atom name and bracket */
      while ((*curpos)&&((*curpos)!=')')) curpos++;
      /* on ')' or end of line */
      while (!(*curpos)) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
        while ((*curpos)&&((*curpos)!=')')) curpos++;
      }
      curpos++;
			while (sscanf(curpos,"%f",&fmeandist)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%f",&fmindelta)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
			while (sscanf(curpos,"%f",&fmaxdelta)!=1) {
				curpos=buf;
				if (FileGets(buf,255,f)==NULL) goto readerr;
      }
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
      /* optional probability at end */
      fprob=0.0;
			sscanf(curpos,"%f",&fprob);
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
      /* optional probability at end */
      fangl=CONSTR_INFINITY;
			sscanf(curpos,"%f",&fangl);
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
      /* optional probability at end */
      fang2=CONSTR_INFINITY;
			sscanf(curpos,"%f",&fang2);
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
      /* optional probability at end */
      fdihed01=CONSTR_INFINITY;
			sscanf(curpos,"%f",&fdihed01);
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
      /* optional probability at end */
      fdihed12=CONSTR_INFINITY;
			sscanf(curpos,"%f",&fdihed12);
      /* skip WS */
      while ((*curpos)&&(isspace((int)(*curpos)))) curpos++;
      /* skip "word" read above */
      while ((*curpos)&&(!isspace((int)(*curpos)))) curpos++;
      /* optional probability at end */
      fdihed23=CONSTR_INFINITY;
			sscanf(curpos,"%f",&fdihed23);
			/* phew - one constraint read successfully */				
			residx1=GetResIdxFromSeq(encseq,ires1);
			residx2=GetResIdxFromSeq(encseq,ires2);
			if ((SubstituteNames(atomname1,residx1)!=ERR_SUCCESS) || (SubstituteNames(atomname2,residx2)!=ERR_SUCCESS)) {
				FileClose(f);
				PurgeGlobs();
				ErrPostEx(SEV_ERROR,1,15,"Error reading constraints, check atom names for correctness.");
				return ERR_FAIL;
			}					
			pnnHere=(PNN)MemNew(sizeof(NN));
			if (ires1<ires2) {
				pnnHere->res1=(Int2)ires1;
				pnnHere->res2=(Int2)ires2;
				StringCpy(pnnHere->AtomName1,atomname1);
				StringCpy(pnnHere->AtomName2,atomname2);
			}
			else {
				pnnHere->res1=(Int2)ires2;
				pnnHere->res2=(Int2)ires1;
				StringCpy(pnnHere->AtomName2,atomname2);
				StringCpy(pnnHere->AtomName1,atomname1);
			}
			pnnHere->MeanDist=(FloatLo)fmeandist;
			pnnHere->MinDelta=(FloatLo)fmindelta;
			pnnHere->MaxDelta=(FloatLo)fmaxdelta;
			pnnHere->Angle1=(FloatLo)fangl;
			pnnHere->Angle2=(FloatLo)fang2;
			pnnHere->Dihedral01=(FloatLo)fdihed01;
			pnnHere->Dihedral12=(FloatLo)fdihed12;
			pnnHere->Dihedral23=(FloatLo)fdihed23;
			if (fprob==0.0)
				pnnHere->prob=NOE_PROB;
			else
				pnnHere->prob=(FloatLo)fprob;
			pnnHere->tries=0;
			if (AddDistConstraint(pnnHere)!=ERR_SUCCESS) {
				FileClose(f);
				PurgeGlobs();
				ErrPostEx(SEV_ERROR,1,14,"Error loading constraints.");
				return ERR_FAIL;
			}
		}
	}	
	FileClose(f);
	return ERR_SUCCESS;
readerr:
	FileClose(f);
	PurgeGlobs();
	ErrPostEx(SEV_FATAL,1,13,"Error parsing file %s near:\n%s\nPlease correct and try again.",fnam,buf);
	return ERR_FAIL;
}

/* returns 1 if a PMADs co-ordinates are valid, 0 otherwise.
   invalid co-ordinates arise when a b-d tree node is removed but the
   PALD is left intact for efficiency.  This prevents accidental
   usage of these "ghost" co-ordinates */
Int2 IsAtomPlaced(PMAD pmadThis,Int2 Model)
{
  PALD paldThis;

  if (pmadThis==NULL) return 0;
  paldThis=GetAtomLocs(pmadThis,Model);
  if (paldThis==NULL) return 0;  
  if (((pmadThis->bReserved)&0x01)==0x01) return 1;
  return 0;
}

/* returns approximate distance from any atom to its Calpha */
FloatLo GetDistFromCAlpha(CharPtr AtomName,PMAD pmadHere)
{
  Char res;
  PMGD pmgdHere;
  PMBD pmbdHere;
  PMAD pmadNext;

  /* these distances are all approximate, and expirically determined */
  if (AtomName[1]=='H') {
	/* H only can have one bond */
	pmbdHere=(PMBD)((pmadHere->pvnBonds)->data.ptrvalue);
	pmadNext=pmbdHere->pmadTo;
	if (pmadNext==pmadHere)
		pmadNext=pmbdHere->pmadFrom;
	return 1.000+GetDistFromCAlpha(pmadNext->pcAName,pmadNext);
  }
  if (!StringCmp(AtomName," C  ")) return 1.527F;
  if (!StringCmp(AtomName," N  ")) return 1.464F;
  if (!StringCmp(AtomName," O  ")) return 2.401F;
  if (!StringCmp(AtomName," CA ")) return 0.000F;
  if (!StringCmp(AtomName," CB ")) return 1.534F;
  if (!StringCmp(AtomName," CG ")) return 2.560F;
  /* not equal to " CG " in this case */
  if (!StringNCmp(AtomName," CG",3)) return 2.686F;
  pmgdHere=(PMGD)(pmadHere->pfbParent);
  res=GetAAFromIDict(pmgdHere);
  if (!StringCmp(AtomName," CD ")) {
	if (res=='P') return 2.467F;
	else return  3.97F;
  }
  if (!StringNCmp(AtomName," CD",3)) {
	if (res=='H') return 3.726F;
	else if (res=='L' || res=='I') return 4.019F;
	else return 3.88F;
  }
  if (!StringCmp(AtomName," CE ")) return 5.237F;
  if (!StringCmp(AtomName," CE1")) {
	if (res=='H') return 4.623F;
	else  return 4.998F;
  }
  if (!StringCmp(AtomName," CE2")) {
	if (res=='W') return 4.740F;
	else  return 4.998F;
  }
  if (!StringCmp(AtomName," CE3")) return 4.735F;
  if (!StringCmp(AtomName," CZ ")) {
	if (res=='R') return 6.24F;
	else return 5.42F;
  }
  if (!StringNCmp(AtomName," CZ",3)) return 6.094F;
  if (!StringCmp(AtomName," CH2")) return 6.632F;
  if (!StringCmp(AtomName," ND1")) return 3.718F;
  if (!StringCmp(AtomName," ND2")) return 3.804F;
  if (!StringCmp(AtomName," NE ")) return 5.024F;
  if (!StringNCmp(AtomName," NE",3)) {
	if (res=='H') return 4.646F;
	else return 4.947F;
  }
  if (!StringCmp(AtomName," NZ ")) return 6.347F;
  if (!StringNCmp(AtomName," NH",3)) return 7.26F;
  if (!StringCmp(AtomName," OG1")) return 2.488F;
  if (!StringNCmp(AtomName," OD",3)) return 3.752F;
  if (!StringNCmp(AtomName," OE",3)) return 5.045F;
  if (!StringCmp(AtomName," OH ")) return 6.771F;
  if (!StringCmp(AtomName," SG ")) return 2.756F;
  if (!StringCmp(AtomName," SD ")) return 4.122F;
  PurgeGlobs();
  ErrPostEx(SEV_FATAL,23,23,"Unknown atom type %s encountered, please contact authors",AtomName);
  return 0.0;
}

/* only should be used by PHIPSI walk */
TrajErr CheckEndCADistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol)
{
  PNN pnnHere;
  PALD pald1,pald2;
  Int2 resHere;
  vec vAtom0,vAtom1,vAtom2;
  FloatLo bl,ba1,ba2,dihed,excessdihed,excessang;

  resHere=(pmgdHere->pdnmgLink)->choice;
  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    /* all further constraints apply later on in the protein chain */
    /* -2 because protein built 3 residues at a time, so make
       sure side chain has been placed */
    if (pnnHere->res1>resHere) return ERR_SUCCESS;
		/* atom after constraint, need to check Dihedral12 and Angle2 */
    if (pnnHere->res2==resHere && IsAtomPlaced(pnnHere->pmad1,Model)) {
    	if (pnnHere->Dihedral12!=CONSTR_INFINITY || pnnHere->Angle2!=CONSTR_INFINITY) {
	      /* constraint not yet passed by */
	      pald1=GetAtomLocs(pnnHere->pmad1,Model);
	      if (pald1==NULL)
	      	ErrPostEx(SEV_FATAL,0,0,"[CheckEndCADistConstraints]: An unexpected error occurred, terminating.");
	      pald2=GetAtomLocs(pnnHere->pmad2,Model);
	      if (pald2==NULL)
	      	ErrPostEx(SEV_FATAL,0,1,"[CheckEndCADistConstraints]: An unexpected error occurred, terminating.");
	      vAtom1[0]=AtomLocX(pald1);
	      vAtom1[1]=AtomLocY(pald1);
	      vAtom1[2]=AtomLocZ(pald1);
	      vAtom2[0]=AtomLocX(pald2);
	      vAtom2[1]=AtomLocY(pald2);
	      vAtom2[2]=AtomLocZ(pald2);
				GetCoOrds((PMGD)(pnnHere->pmad1->pfbParent)," CA ",vZero,vAtom0,Model);
				/* now co-ords are in vAtom0, vAtom1, (gap), vAtom2, vHere */
				excessang=0.0;
				excessdihed=0.0;
				GetDihedral(vAtom0,vAtom1,vAtom2,vHere,-180,&dihed,&ba1,&ba2,&bl,&bl,&bl);
				if (pnnHere->Angle2<CONSTR_INFINITY)
					excessang=fabs(pnnHere->Angle2-ba2);
				if (pnnHere->Dihedral12<CONSTR_INFINITY)
					excessdihed=fabs(pnnHere->Dihedral12-dihed);
				if (excessdihed>180.0)
					excessdihed=360.0-excessdihed;
				if (excessang>CONSTR_ANGLE_TOL || excessdihed>CONSTR_DIHED_TOL) {
					if ((pnnHere->tries)<CONSTR_TRIES_MAX) {
						(pnnHere->tries)++;
						if (pnnHere->tries==CONSTR_TRIES_MAX)
							(*viol)++;
						return ERR_DISTCONST;
					}
				}
/*printf("ang: %f,dihed: %f\n",excessang,excessdihed);*/
			}		  		
    }
    pnnHere=pnnHere->next;
  }
  /* no distance constraints violated */
  return ERR_SUCCESS;
}

/* only should be used by PHIPSI walk */
TrajErr CheckEndCDistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol)
{
  PNN pnnHere;
  PALD pald1,pald2;
  Int2 resHere;
	vec vAtom0,vAtom1,vAtom2;
  FloatLo bl,ba1,ba2,dihed,excessdihed;

  resHere=(pmgdHere->pdnmgLink)->choice;
  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    /* all further constraints apply later on in the protein chain */
    /* -2 because protein built 3 residues at a time, so make
       sure side chain has been placed */
    if (pnnHere->res1>resHere) return ERR_SUCCESS;
		/* atom after constraint, need to check Dihedral12 and Angle2 */
		/* atom 2nd after constraint, need to check Dihedral23 */
		if (pnnHere->res2==resHere && IsAtomPlaced(pnnHere->pmad1,Model)) {
    	if (pnnHere->Dihedral23!=CONSTR_INFINITY) {
	      /* constraint not yet passed by */
	      pald1=GetAtomLocs(pnnHere->pmad1,Model);
	      if (pald1==NULL)
	      	ErrPostEx(SEV_FATAL,0,3,"[CheckEndCDistConstraints]: An unexpected error occurred, terminating.");
	      pald2=GetAtomLocs(pnnHere->pmad2,Model);
	      if (pald2==NULL)
	      	ErrPostEx(SEV_FATAL,0,4,"[CheckEndCDistConstraints]: An unexpected error occurred, terminating.");
	      vAtom0[0]=AtomLocX(pald1);
	      vAtom0[1]=AtomLocY(pald1);
	      vAtom0[2]=AtomLocZ(pald1);
	      vAtom1[0]=AtomLocX(pald2);
	      vAtom1[1]=AtomLocY(pald2);
	      vAtom1[2]=AtomLocZ(pald2);
				GetCoOrds((PMGD)(pnnHere->pmad2->pfbParent)," CA ",vZero,vAtom2,Model);												
				/* now co-ords are in vAtom0, (gap), vAtom1, vAtom2, vHere */
				excessdihed=0.0;
				GetDihedral(vAtom0,vAtom1,vAtom2,vHere,-180,&dihed,&ba1,&ba2,&bl,&bl,&bl);					
				if (pnnHere->Dihedral23<CONSTR_INFINITY)
					excessdihed=fabs(pnnHere->Dihedral23-dihed);
				if (excessdihed>180.0)
					excessdihed=360.0-excessdihed;
				if (excessdihed>CONSTR_DIHED_TOL) {
					if ((pnnHere->tries)<CONSTR_TRIES_MAX) {
						(pnnHere->tries)++;
						if (pnnHere->tries==CONSTR_TRIES_MAX)
							(*viol)++;
						return ERR_DISTCONST;
					}
				}
/*printf("dihed: %f\n",excessdihed);*/
	  	}		
		}
    pnnHere=pnnHere->next;
  }
  /* no distance constraints violated */
  return ERR_SUCCESS;
}

TrajErr CheckEndingDistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol)
{
  PNN pnnHere;
  PALD pald1,pald2;
  Int2 resHere;
  vec vAtom0,vAtom1,vAtom2;
  FloatLo bl,ba1,ba2,dihed,excessdihed,excessang;
  PDNMG pdnmg0,pdnmg1;

  resHere=(pmgdHere->pdnmgLink)->choice;
  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    /* all further constraints apply later on in the protein chain */
    /* -2 because protein built 3 residues at a time, so make
       sure side chain has been placed */
    if (pnnHere->res1>=resHere) return ERR_SUCCESS;
		/* atom after constraint, need to check Dihedral12 and Angle2 */
    if (pnnHere->res2==(resHere-1) && IsAtomPlaced(pnnHere->pmad1,Model)) {
    	if (pnnHere->Dihedral12!=CONSTR_INFINITY || pnnHere->Angle2!=CONSTR_INFINITY) {
	      /* constraint not yet passed by */
	      pald1=GetAtomLocs(pnnHere->pmad1,Model);
	      if (pald1==NULL)
	      	ErrPostEx(SEV_FATAL,0,0,"[CheckEndingDistConstraints]: An unexpected error occurred, terminating.");
	      pald2=GetAtomLocs(pnnHere->pmad2,Model);
	      if (pald2==NULL)
	      	ErrPostEx(SEV_FATAL,0,1,"[CheckEndingDistConstraints]: An unexpected error occurred, terminating.");
	      vAtom1[0]=AtomLocX(pald1);
	      vAtom1[1]=AtomLocY(pald1);
	      vAtom1[2]=AtomLocZ(pald1);
	      vAtom2[0]=AtomLocX(pald2);
	      vAtom2[1]=AtomLocY(pald2);
	      vAtom2[2]=AtomLocZ(pald2);
	      pdnmg1=((PMGD)(pnnHere->pmad1->pfbParent))->pdnmgLink;
      	pdnmg0=pdnmg1->last;
	      if (pdnmg0!=NULL) {
					if (pdnmg0->choice!=(pdnmg1->choice-1))
		      	ErrPostEx(SEV_FATAL,0,2,"[CheckEndingDistConstraints]: An unexpected error occurred, terminating.");
					GetCoOrds((PMGD)(pdnmg0->data.ptrvalue)," CA ",vZero,vAtom0,Model);
					/* now co-ords are in vAtom0, vAtom1, (gap), vAtom2, vHere */
					excessang=0.0;
					excessdihed=0.0;
					GetDihedral(vAtom0,vAtom1,vAtom2,vHere,-180,&dihed,&ba1,&ba2,&bl,&bl,&bl);
					if (pnnHere->Angle2<CONSTR_INFINITY)
						excessang=fabs(pnnHere->Angle2-ba2);
					if (pnnHere->Dihedral12<CONSTR_INFINITY)
						excessdihed=fabs(pnnHere->Dihedral12-dihed);
					if (excessdihed>180.0)
						excessdihed=360.0-excessdihed;
					if (excessang>CONSTR_ANGLE_TOL || excessdihed>CONSTR_DIHED_TOL) {
						if ((pnnHere->tries)<CONSTR_TRIES_MAX) {
							(pnnHere->tries)++;
							if (pnnHere->tries==CONSTR_TRIES_MAX)
								(*viol)++;
							return ERR_DISTCONST;
						}
					}
/*printf("ang: %f,dihed: %f\n",excessang,excessdihed);*/
				}	
	  	}	
    }
		/* atom 2nd after constraint, need to check Dihedral23 */
		else if (pnnHere->res2==(resHere-2) && IsAtomPlaced(pnnHere->pmad1,Model)) {
    	if (pnnHere->Dihedral23!=CONSTR_INFINITY) {
	      /* constraint not yet passed by */
	      pald1=GetAtomLocs(pnnHere->pmad1,Model);
	      if (pald1==NULL)
	      	ErrPostEx(SEV_FATAL,0,3,"[CheckEndingDistConstraints]: An unexpected error occurred, terminating.");
	      pald2=GetAtomLocs(pnnHere->pmad2,Model);
	      if (pald2==NULL)
	      	ErrPostEx(SEV_FATAL,0,4,"[CheckEndingDistConstraints]: An unexpected error occurred, terminating.");
	      vAtom0[0]=AtomLocX(pald1);
	      vAtom0[1]=AtomLocY(pald1);
	      vAtom0[2]=AtomLocZ(pald1);
	      vAtom1[0]=AtomLocX(pald2);
	      vAtom1[1]=AtomLocY(pald2);
	      vAtom1[2]=AtomLocZ(pald2);
	      pdnmg1=((PMGD)(pnnHere->pmad1->pfbParent))->pdnmgLink;
      	pdnmg0=pdnmg1->next;
	      if (pdnmg0!=NULL) {
					if (pdnmg0->choice!=(pdnmg1->choice+1))
		      	ErrPostEx(SEV_FATAL,0,5,"[CheckEndingDistConstraints]: An unexpected error occurred, terminating.");
					GetCoOrds((PMGD)(pdnmg0->data.ptrvalue)," CA ",vZero,vAtom2,Model);												
					/* now co-ords are in vAtom0, (gap), vAtom1, vAtom2, vHere */
					excessdihed=0.0;
					GetDihedral(vAtom0,vAtom1,vAtom2,vHere,-180,&dihed,&ba1,&ba2,&bl,&bl,&bl);					
					if (pnnHere->Dihedral23<CONSTR_INFINITY)
						excessdihed=fabs(pnnHere->Dihedral23-dihed);
					if (excessdihed>180.0)
						excessdihed=360.0-excessdihed;
					if (excessdihed>CONSTR_DIHED_TOL) {
						if ((pnnHere->tries)<CONSTR_TRIES_MAX) {
							(pnnHere->tries)++;
							if (pnnHere->tries==CONSTR_TRIES_MAX)
								(*viol)++;
							return ERR_DISTCONST;
						}
					}
/*printf("dihed: %f\n",excessdihed);*/
				}	
	  	}		
		}
    pnnHere=pnnHere->next;
  }
  /* no distance constraints violated */
  return ERR_SUCCESS;
}

/* this is used by both WALKTYPE=WALK_CA and WALK_PHIPSI but is
   probably more effective for the latter */
TrajErr CheckDistConstraints(PMGD pmgdHere,vec vHere,Int2 Model,Int4 *viol)
{
  PNN pnnHere;
  PALD pald1;
  Int2 resHere;
  vec vAtom1,vTmp;
  FloatLo OptimalDist,DistFromCAlpha,MaxDistToGo;

  resHere=(pmgdHere->pdnmgLink)->choice;
  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    /* all further constraints apply later on in the protein chain */
    /* -2 because protein built 3 residues at a time, so make
       sure side chain has been placed */
    if (pnnHere->res1>=resHere) return ERR_SUCCESS;
    if (pnnHere->res2>=resHere && IsAtomPlaced(pnnHere->pmad1,Model)) {
      /* constraint not yet passed by */
      pald1=GetAtomLocs(pnnHere->pmad1,Model);
      /* maybe atom not placed yet */
      vAtom1[0]=AtomLocX(pald1);      
      vAtom1[1]=AtomLocY(pald1);      
      vAtom1[2]=AtomLocZ(pald1);
      VecSub(vTmp,vHere,vAtom1);
      /* shortest distance possible from current Ca to atom #2 in the constraint */
      OptimalDist=getMag(vTmp)-(pnnHere->MeanDist+pnnHere->MaxDelta);
      DistFromCAlpha=GetDistFromCAlpha(pnnHere->AtomName2,pnnHere->pmad2);
      MaxDistToGo=BL_CACA*(pnnHere->res2-resHere)+DistFromCAlpha;
      /* if true, can't possibly satisfy the distance constraint so backtrack */
      if (OptimalDist>MaxDistToGo) {
      	/* for disulphide bridge, always backtrack */
	if (((pnnHere->tries)<CONSTR_TRIES_MAX) || ((pnnHere->AtomName1)[1]=='S' && (pnnHere->AtomName2)[1]=='S')) {
		    	(pnnHere->tries)++;
				if (pnnHere->tries==CONSTR_TRIES_MAX)
					(*viol)++;
	      	return ERR_DISTCONST;
	      }
      }
      /* otherwise, it is possible so continue along */
    }
    pnnHere=pnnHere->next;
  }
  /* no distance constraints violated */
  return ERR_SUCCESS;
}

/* return 1 if atom names match, with wildcards, 0 else */			
Int2 WildCardMatch(CharPtr test,CharPtr pattern) {
	if ((test[1]!=pattern[1]) || (test[2]!=pattern[2]))
		return 0;
	if ((pattern[0]!='*') && (test[0]!=pattern[0]))
		return 0;
	if ((pattern[3]!='*') && (test[3]!=pattern[3]))
		return 0;
	return 1;
}
			
TrajErr CheckAtomDistConstraints(CharPtr AtomName,Int2 res,vec vHere,Int2 Model,Int4 *viol)
{
  PNN pnnHere;
  PALD paldTest;
  vec vTmp,vCofMHere,vCofMPlaced,vCAPrev,vCAPrevPrev;
  vec vAtomsPlaced[MAX_DEGEN];
  vec vAtomsHere[MAX_DEGEN];
  FloatLo d,mind,maxd,excessd,excessang,excessdihed;
  FloatLo dihed,ba1,ba2,bl;
  CharPtr CurName,PlacedName;
  PMAD pmadTest=NULL;
  PVNMA pvnmaTest;
  PDNMG pdnmgPrev,pdnmgPrevPrev;
  Int2 numhere,numplaced,cnt,cur;

  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    /* since list is ordered, we can save a bit of time here */
    if (pnnHere->res1>res) return ERR_SUCCESS;
    /* only find restraints involving this atom */
    if ((pnnHere->res1!=res) && (pnnHere->res2!=res)) goto nextconst;
    /* constraint not yet passed by */
		/* special case if intra-residue restraint */
/*    if ((pnnHere->res1==res) && (pnnHere->res2==res)) {
			ErrPostEx(SEV_WARNING,2,12,"Intra-residue restraints not supported (yet), restraint ignored");
			goto nextconst;
    }*/
		if (pnnHere->res1==res) {
      CurName=pnnHere->AtomName1;
      PlacedName=pnnHere->AtomName2;
      cur=1;
		}
		else { /* pnnHere->res2==res */
 	    CurName=pnnHere->AtomName2;
      PlacedName=pnnHere->AtomName1;
      cur=2;
		}
		/* we're in the right residue - now check atom names */
		if ((pnnHere->res1==res) && (pnnHere->res2==res)) {
			/* if intra-residue, figure out which matches which */
			if (WildCardMatch(AtomName,pnnHere->AtomName1)) {
	      CurName=pnnHere->AtomName1;
  	    PlacedName=pnnHere->AtomName2;
    	  cur=1;
			}
			else if (WildCardMatch(AtomName,pnnHere->AtomName2)) {
	      CurName=pnnHere->AtomName2;
  	    PlacedName=pnnHere->AtomName1;
    	  cur=2;
			}
			else goto nextconst;
		}
		else
			if (!WildCardMatch(AtomName,CurName)) goto nextconst;
	  /* now we have the right atom for sure */
		/* check if atom we are placing is last in a wildcard group */
		if (cur==1)
			pvnmaTest=(pnnHere->pmad1)->pvnmaLink;
		else
			pvnmaTest=(pnnHere->pmad2)->pvnmaLink;
		numhere=0;
		while (pvnmaTest!=NULL) {
			pmadTest=(PMAD)(pvnmaTest->data.ptrvalue);
			if (WildCardMatch(pmadTest->pcAName,CurName)) {
				if (!IsAtomPlaced(pmadTest,Model) && StringCmp(pmadTest->pcAName,AtomName)) goto nextconst;
				/* more to go in group still */
				else {
					if (numhere>=MAX_DEGEN) {
						PurgeGlobs();
						ErrPostEx(SEV_FATAL,4,44,"Cannot handle more than %d degenrate atoms, please contact authors",MAX_DEGEN);
						return ERR_DISTCONST;
					}
					/* record co-ords */
					if (StringCmp(pmadTest->pcAName,AtomName)) {
					  paldTest=GetAtomLocs(pmadTest,Model);
 	 			  	vAtomsHere[numhere][0]=AtomLocX(paldTest);
  	 	  		vAtomsHere[numhere][1]=AtomLocY(paldTest);
	  	 	  	vAtomsHere[numhere][2]=AtomLocZ(paldTest);
	  	 	  }
	  	 	  else {
	  	 	  	/* current atom is last in group */
	  	 	    vAtomsHere[numhere][0]=vHere[0];
	  	 	    vAtomsHere[numhere][1]=vHere[1];
	  	 	    vAtomsHere[numhere][2]=vHere[2];
	  	 	  }
					numhere++;
				}
			}
			pvnmaTest=pvnmaTest->next;
		}
	  /* now check if any atoms are not placed yet/ get their co-ords */
		if (cur==1)
			pvnmaTest=(pnnHere->pmad2)->pvnmaLink;
		else
			pvnmaTest=(pnnHere->pmad1)->pvnmaLink;
		numplaced=0;
		while (pvnmaTest!=NULL) {
			pmadTest=(PMAD)(pvnmaTest->data.ptrvalue);
			if (WildCardMatch(pmadTest->pcAName,PlacedName)) {
				if (!IsAtomPlaced(pmadTest,Model)) goto nextconst;
				else {
					if (numplaced>=MAX_DEGEN) {
						PurgeGlobs();
						ErrPostEx(SEV_FATAL,4,44,"Cannot handle more than %d degenerate atoms, please contact authors",MAX_DEGEN);
						return ERR_DISTCONST;
					}
					/* record co-ords */
				  paldTest=GetAtomLocs(pmadTest,Model);
 			  	vAtomsPlaced[numplaced][0]=AtomLocX(paldTest);
   	  		vAtomsPlaced[numplaced][1]=AtomLocY(paldTest);
	   	  	vAtomsPlaced[numplaced][2]=AtomLocZ(paldTest);
					numplaced++;
				}
			}
			pvnmaTest=pvnmaTest->next;
		}
		/* now we are assured all previous atoms are placed, and last atom of
		   current group is being placed.  All co-ordinates are recorded at this
		   point as well */
    mind=pnnHere->MeanDist-pnnHere->MinDelta;
    maxd=pnnHere->MeanDist+pnnHere->MaxDelta;
		/* now calculate distances between centres of mass */
		VecScale(vCofMPlaced,vZero,1.0);
		for (cnt=0;cnt<numplaced;cnt++)
			VecAdd(vCofMPlaced,vCofMPlaced,vAtomsPlaced[cnt]);
		VecScale(vCofMPlaced,vCofMPlaced,1.0/(FloatLo)numplaced);
		VecScale(vCofMHere,vZero,1.0);
		for (cnt=0;cnt<numhere;cnt++)
			VecAdd(vCofMHere,vCofMHere,vAtomsHere[cnt]);
		VecScale(vCofMHere,vCofMHere,1.0/(FloatLo)numhere);
		/* get distance vector */
		VecSub(vTmp,vCofMHere,vCofMPlaced);
		d=getMag(vTmp);
		if (d>maxd)
			excessd=d-maxd;
		else if (d<mind)
			excessd=mind-d;
		else
			excessd=0.0;
		excessang=0.0;
		excessdihed=0.0;
		if ((pnnHere->Angle1<CONSTR_INFINITY || pnnHere->Dihedral01<CONSTR_INFINITY) && pmadTest!=NULL) {
			if (WALKTYPE==WALK_CA) {
				pdnmgPrev=((PMGD)(pmadTest->pfbParent))->pdnmgLink->last;
				if (pdnmgPrev!=NULL) {
					pdnmgPrevPrev=pdnmgPrev->last;
					GetCoOrds((PMGD)(pdnmgPrev->data.ptrvalue)," CA ",vZero,vCAPrev,Model);
					if (pdnmgPrevPrev!=NULL) {
						GetCoOrds((PMGD)(pdnmgPrevPrev->data.ptrvalue)," CA ",vZero,vCAPrevPrev,Model);
						GetDihedral(vCAPrevPrev,vCAPrev,vCofMPlaced,vCofMHere,-180,&dihed,&ba1,&ba2,&bl,&bl,&bl);					
						if (pnnHere->Angle1<CONSTR_INFINITY)
							excessang=fabs(pnnHere->Angle1-ba2);
						if (pnnHere->Dihedral01<CONSTR_INFINITY)
							excessdihed=fabs(pnnHere->Dihedral01-dihed);
					}
				}
			}
			else if (WALKTYPE==WALK_PHIPSI) {
				GetCoOrds((PMGD)(pmadTest->pfbParent)," CA ",vZero,vCAPrev,Model);
				GetCoOrds((PMGD)(pmadTest->pfbParent)," N  ",vZero,vCAPrevPrev,Model);
				GetDihedral(vCAPrevPrev,vCAPrev,vCofMPlaced,vCofMHere,-180,&dihed,&ba1,&ba2,&bl,&bl,&bl);
				if (pnnHere->Angle1<CONSTR_INFINITY)
					excessang=fabs(pnnHere->Angle1-ba2);
				if (pnnHere->Dihedral01<CONSTR_INFINITY)
					excessdihed=fabs(pnnHere->Dihedral01-dihed);
			}
		}
		/* if too close or too far */
		/* using thus formula, there is a 33% chance a violation of 0.5 A
		   will be accepted, and no violations > 0.613 A will occur */
		if (fabs(Rand1())>1.0-excessd*excessd*8.0/3.0) {
/*printf("dist: %f\n",excessd);*/
     	/* for disulphide bridge, always backtrack */
			if (((pnnHere->tries)<CONSTR_TRIES_MAX) || ((pnnHere->AtomName1)[1]=='S' && (pnnHere->AtomName2)[1]=='S')) {
				(pnnHere->tries)++;
				if (pnnHere->tries==CONSTR_TRIES_MAX)
					(*viol)++;
				return ERR_DISTCONST;
			}
		}
		if (excessdihed>180.0)
			excessdihed=360.0-excessdihed;
		if (excessang>CONSTR_ANGLE_TOL || excessdihed>CONSTR_DIHED_TOL) {
			if (((pnnHere->tries)<CONSTR_TRIES_MAX) || ((pnnHere->AtomName1)[1]=='S' && (pnnHere->AtomName2)[1]=='S')) {
				(pnnHere->tries)++;
				if (pnnHere->tries==CONSTR_TRIES_MAX)
					(*viol)++;
				return ERR_DISTCONST;
			}
		}
/*printf("ang: %f,dihed: %f\n",excessang,excessdihed);*/
nextconst:
    pnnHere=pnnHere->next;
  }
  /* no distance constraints violated */
  return ERR_SUCCESS;
}

TrajErr OldCheckAtomDistConstraints(CharPtr AtomName,Int2 res,vec vHere,Int2 Model)
{
  PNN pnnHere;
  PALD pald1,paldCur,paldCur2;
  vec vTmp,vAtom1,vCur,vCur2;
  FloatLo d,mind,maxd;
  CharPtr CurName,TestName,TestName2,PlacedName;
  PMAD pmadCur,pmadHere,pmadCur2;
  PVNMA pvnmaCur,pvnmaCur2;
  Boolean ok;

  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    /* since list is ordered, we can save a bit of time here */
    if (pnnHere->res1>res) return ERR_SUCCESS;
    /* only find restraints involving this atom */
    if ((pnnHere->res1==res) || (pnnHere->res2==res)) {
      /* constraint not yet passed by */
      if ((pnnHere->res1==res) && (pnnHere->res2==res)) {
				/* special case if intra-residue restraint */
				ErrPostEx(SEV_WARNING,2,12,"Intra-residue restraints not supported (yet), restraint ignored");
      }
      else {
				if (pnnHere->res1==res) {
				  /* if placing res 1 now but res 2 not placed yet, nothing to do */
				  if (!IsAtomPlaced(pnnHere->pmad2,Model)) {
	  			  pnnHere=pnnHere->next;
	  			  continue;
	 			  }
          CurName=pnnHere->AtomName1;
          PlacedName=pnnHere->AtomName2;
				}
				else { /* pnnHere->res2==res */
				  /* if placing res 2 now but res 1 not placed yet, nothing to do */
				  if (!IsAtomPlaced(pnnHere->pmad1,Model)) {
				    pnnHere=pnnHere->next;
				    continue;
				  }
    	    CurName=pnnHere->AtomName2;
          PlacedName=pnnHere->AtomName1;
				}
				if ((AtomName[1]==CurName[1]) && (AtomName[2]==CurName[2])) {
          mind=pnnHere->MeanDist-pnnHere->MinDelta;
          maxd=pnnHere->MeanDist+pnnHere->MaxDelta;
  			  /* this atom may be in a constraint with some previous residue */
				  /* but we must account for degenerate atoms - only 1 of them need
				     be within the specified distance for now */
				  if ((CurName[0]!='*') && (CurName[0]!=AtomName[0])) {
				    pnnHere=pnnHere->next;
				    continue;
				  }
				  if ((CurName[3]!='*') && (CurName[3]!=AtomName[3])) {
	  			  pnnHere=pnnHere->next;
			    	continue;
				  }
				  /* now we have the right atom for sure */
	  			/* get co-ords of already placed atom */
	  			/* but what if already placed atom is degenerate??? - must handle */
				  if (pnnHere->res2==res)
	  			  pmadHere=pnnHere->pmad1;
				  else
		  		  pmadHere=pnnHere->pmad2;
		  		ok=FALSE;
				  pald1=GetAtomLocs(pmadHere,Model);
   			  vAtom1[0]=AtomLocX(pald1);
	   	  	vAtom1[1]=AtomLocY(pald1);
		   	  vAtom1[2]=AtomLocZ(pald1);
		   	  VecSub(vTmp,vHere,vAtom1);
				  d=getMag(vTmp);
				  /* d is distance between pmad1 and atom under consideration */
			  	if ((d>mind) && (d<maxd)) ok=TRUE;
			  	else if ((PlacedName[0]=='*') || (PlacedName[3]=='*')) {
			  		/* else check if another of the degenerate atoms will do */
				  	pvnmaCur2=pmadHere->pvnmaLink;
			    	while ((pvnmaCur2=pvnmaCur2->next)!=NULL) {
					    /* ->pmad1/2 always points to first wildcard match so no need to
					    	 start at the beginning of the list */
			      	pmadCur2=(PMAD)(pvnmaCur2->data.ptrvalue);
				      TestName2=pmadCur2->pcAName;
				      if ((TestName2[1]==PlacedName[1]) && (TestName2[2]==PlacedName[2])) {
				        if ((PlacedName[0]!='*') && (TestName2[0]!=PlacedName[0])) continue;
			        	if ((PlacedName[3]!='*') && (TestName2[3]!=PlacedName[3])) continue;
				        /* a definite match */
				        if (IsAtomPlaced(pmadCur2,Model)) {
								  paldCur2=GetAtomLocs(pmadCur2,Model);
						   	  vCur2[0]=AtomLocX(paldCur2);
					  	 	  vCur2[1]=AtomLocY(paldCur2);
					  	 	  vCur2[2]=AtomLocZ(paldCur2);
					  	 	  VecSub(vTmp,vHere,vCur2);
								  d=getMag(vTmp);
								  if ((d<mind) || (d>maxd));
								  else {
								    /* constraint already met elsewhere */
							  	  ok=TRUE;
							    	break;
								  }
				        }
	  			      else {
								  /* a match which hasnt been placed yet is good enough */
								  ok=TRUE;
								  break;
			      	  }
			      	}
			  		}
			  	}
			  	if (ok==FALSE) {
			    	/* constraint not satisfied... but may be OK if degenrate case */
				    /* OK if and only if another atom matches & (satisfies or not
				       placed yet) - check list of atoms at current residue */
				    if (CurName[0]!='*' && CurName[3]!='*') return ERR_DISTCONST;
				    /* ->pmad1/2 always points to first wildcard match so no need to
				    	 start at the beginning of the list */
					  if (pnnHere->res2==res)
				    	pvnmaCur=(pnnHere->pmad2)->pvnmaLink;
				  	else
				    	pvnmaCur=(pnnHere->pmad1)->pvnmaLink;
			  	  ok=FALSE;
				/* here we need to be more careful since this is the atom we
				   are currently placing, and not necessarily the first in the
				   list as in ->pmad 1/2 */
			    	while (pvnmaCur!=NULL) {
			      	pmadCur=(PMAD)(pvnmaCur->data.ptrvalue);
				      TestName=pmadCur->pcAName;
				      if ((TestName[1]==CurName[1]) && (TestName[2]==CurName[2]) && (StringCmp(TestName,AtomName))) {
				        if ((CurName[0]!='*') && (TestName[0]!=CurName[0])) {
								  pvnmaCur=pvnmaCur->next;
								  continue;
								}
			        	if ((CurName[3]!='*') && (TestName[3]!=CurName[3])) {
								  pvnmaCur=pvnmaCur->next;
								  continue;
								}
				        /* a definite match */
				        if (IsAtomPlaced(pmadCur,Model)) {
								  paldCur=GetAtomLocs(pmadCur,Model);
						   	  vCur[0]=AtomLocX(paldCur);
					  	 	  vCur[1]=AtomLocY(paldCur);
					  	 	  vCur[2]=AtomLocZ(paldCur);
					  	 	  VecSub(vTmp,vAtom1,vCur);
								  d=getMag(vTmp);
									/* must compare to each of degenrate Placed Atoms */
									if ((d>mind) && (d<maxd)) {
							  		ok=TRUE;
							  		break;
							  	}
							  	else if ((PlacedName[0]=='*') || (PlacedName[3]=='*')) {
							  		/* else check if another of the degenerate atoms will do */
								  	pvnmaCur2=pmadHere->pvnmaLink;
							    	while ((pvnmaCur2=pvnmaCur2->next)!=NULL) {
									    /* ->pmad1/2 always points to first wildcard match so no need to
									    	 start at the beginning of the list */
							      	pmadCur2=(PMAD)(pvnmaCur2->data.ptrvalue);
					  			    TestName2=pmadCur2->pcAName;
							  	    if ((TestName2[1]==PlacedName[1]) && (TestName2[2]==PlacedName[2])) {
								        if ((PlacedName[0]!='*') && (TestName2[0]!=PlacedName[0])) continue;
			        					if ((PlacedName[3]!='*') && (TestName2[3]!=PlacedName[3])) continue;
								        /* a definite match */
								        if (IsAtomPlaced(pmadCur2,Model)) {
												  paldCur2=GetAtomLocs(pmadCur2,Model);
										   	  vCur2[0]=AtomLocX(paldCur2);
									  	 	  vCur2[1]=AtomLocY(paldCur2);
					  	 					  vCur2[2]=AtomLocZ(paldCur2);
									  	 	  VecSub(vTmp,vCur,vCur2);
												  d=getMag(vTmp);
												  if ((d<mind) || (d>maxd));
												  else {
												    /* constraint already met elsewhere */
							  					  ok=TRUE;
											    	break;
												  }
				        				}
					  			      else {
												  /* a match which hasnt been placed yet is good enough */
												  ok=TRUE;
												  break;
			  				    	  }
							      	} 	
			  						}
							  	}
									if (ok==TRUE) break;									
				        }
	  			      else {
								  /* a match which hasnt been placed yet is good enough */
								  ok=TRUE;
								  break;
			      	  }
		  		    }
		  		    pvnmaCur=pvnmaCur->next;
				    }
	  			  if (ok==FALSE) return ERR_DISTCONST;
			  	}
				  /* constraint is satisfied */
        }
      }
    }
    pnnHere=pnnHere->next;
  }
  /* no distance constraints violated */
  return ERR_SUCCESS;
}

/* returns probabilitiy residue is in a disulfide bridge, else 0 */
FloatLo InDisulphide(Int2 resnum)
{
  PNN pnnHere;

  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    if ((!StringCmp(pnnHere->AtomName1," SG ")) && (!StringCmp(pnnHere->AtomName2," SG "))) {
      if ((pnnHere->res1==resnum) || (pnnHere->res2==resnum))
        return pnnHere->prob;
    }
    pnnHere=pnnHere->next;
  }
  return 0.0;
}

/* fills in PMAD values for PNN List after structure has been loaded */
/* also removes things from list based on probabilities */
void InitPNNList(PMMD pmmdHere)
{
  PDNMG pdnmgHere;
  PMGD pmgdHere;
  PVNMA pvnmaHere;
  PNN pnnHere,pnnNext,pnnLast=NULL;
  PMMD pmmdParent;
  PVNMB pvnmbNew,pvnmbNext,pvnmbLast;
  PMBD pmbdNew;
  ValNodePtr pvnBondHere,pvnBondLast;
  CharPtr NameHere,NameThere;

  /* move appropriate True Dist. Constraints to Dist. Constraint */
  /* first clear any old stuff left over */
  pnnHere=pnnDistConstHead;
  while (pnnHere!=NULL) {
    pnnNext=pnnHere->next;
    MemFree(pnnHere);
    pnnHere=pnnNext;
  }
  pnnDistConstHead=NULL;
  pnnHere=pnnDistConstTrueHead;
  /* undo all disulfide bridges PMAD bond pointers in chemical graph if present */
  while (pnnHere!=NULL) {
  	if ((pnnHere->AtomName1)[1]=='S' && (pnnHere->AtomName2)[1]=='S') {
  		if ((pnnHere->pmad1)==NULL) {
			pnnHere=pnnHere->next;
			continue;
		}
      pvnBondHere=(pnnHere->pmad1)->pvnBonds;
      pvnBondLast=NULL;
      while (pvnBondHere!=NULL) {
				if ((!StringCmp(" SG ",(((PMBD)(pvnBondHere->data.ptrvalue))->pmadFrom)->pcAName)) && (!StringCmp(" SG ",(((PMBD)(pvnBondHere->data.ptrvalue))->pmadTo)->pcAName))) {
					pvnBondHere=ValNodeFree(pvnBondHere);
					if (pvnBondLast!=NULL)
						pvnBondLast->next=NULL;
				}
				pvnBondLast=pvnBondHere;
				if (pvnBondHere!=NULL)
					pvnBondHere=pvnBondHere->next;					
			}
      pvnBondHere=(pnnHere->pmad2)->pvnBonds;
      pvnBondLast=NULL;
      while (pvnBondHere!=NULL) {
				if ((!StringCmp(" SG ",(((PMBD)(pvnBondHere->data.ptrvalue))->pmadFrom)->pcAName)) && (!StringCmp(" SG ",(((PMBD)(pvnBondHere->data.ptrvalue))->pmadTo)->pcAName))) {
					pvnBondHere=ValNodeFree(pvnBondHere);
					if (pvnBondLast!=NULL)
						pvnBondLast->next=NULL;
				}
				pvnBondLast=pvnBondHere;
				if (pvnBondHere!=NULL)
					pvnBondHere=pvnBondHere->next;					
			}
    }
    pnnHere=pnnHere->next;
  }
  pnnHere=pnnDistConstTrueHead;
  if (pnnHere!=NULL) {
	  if (pnnHere->pmad1!=NULL) {
		  pmmdParent=GetParentMol((PFB)(pnnHere->pmad1));
	  	pvnmbNew=pmmdParent->pvnmbIRBHead;
		  /* undo all disulfide bridge PMBDs in chemical graph if present */
		  while (pvnmbNew!=NULL) {
  		  pmbdNew=(PMBD)(pvnmbNew->data.ptrvalue);
    		if ((!StringCmp((pmbdNew->pmadFrom)->pcAName," SG ")) && (!StringCmp((pmbdNew->pmadTo)->pcAName," SG "))) {
	 	  		/* kill this bond */
	  	 		MemFree(pmbdNew);
  	  	 	(pmmdParent->iIRBCount)--;
    	 		/* mark for removal */
	    	  pvnmbNew->data.ptrvalue=NULL;
	  	  }
  	  	pvnmbNew=pvnmbNew->next;
	 		}
		 	/* and free the corresponding PVNMBs now too */
		  pvnmbNew=pmmdParent->pvnmbIRBHead;
			pvnmbLast=NULL;
		  while (pvnmbNew!=NULL) {
		    pvnmbNext=pvnmbNew->next;
	  	 	if (pvnmbNew->data.ptrvalue==NULL) {
	   			MemFree(pvnmbNew);
	   			if (pvnmbLast==NULL)
		   			pmmdParent->pvnmbIRBHead=pvnmbNext;
		   	}
	  	 	else {
	   			if (pvnmbLast!=NULL) 	
			   		pvnmbLast->next=pvnmbNew;
		   		pvnmbLast=pvnmbNew;
	  	 	}
	   		pvnmbNew=pvnmbNext;
		  }
		  if (pvnmbLast!=NULL)
			 	pvnmbLast->next=NULL;
		}
	} 	
	/* now copy true list to current list based on probabilities */
  pnnHere=pnnDistConstTrueHead;
  while (pnnHere!=NULL) {
		if (fabs(Rand1())<=pnnHere->prob) {
  		/* if pass probability test add it to the current constraint list */
			pnnNext=(PNN)MemNew(sizeof(NN));
			MemCopy(pnnNext,pnnHere,sizeof(NN));
			pnnNext->next=NULL;
 			if (pnnDistConstHead==NULL) {
 				/* head of list */
 				pnnDistConstHead=pnnNext;
 				pnnNext->prev=NULL;
  		}
 			else {
 				pnnNext->prev=pnnLast;
 				pnnLast->next=pnnNext;
  		}
 			pnnLast=pnnNext;
 		}
  	pnnHere=pnnHere->next;
 	}
 	/* add random S-S bridges */
	/* possibly in future */
/*  pdnmgHere=pmmdHere->pdnmgHead;
  while (pdnmgHere!=NULL) {
    pdnmgHere=pdnmgHere->next;
    pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
    if (pmgdHere->pcIUPAC
    pdnmgHere=pdnmgHere->next;
  }*/
 	/* now deal with current constraint list only, to fill in PMAD values */
  pnnHere=pnnDistConstHead;
  while (pnnHere) {
    pdnmgHere=pmmdHere->pdnmgHead;
    /* do 1st residue */
    while (pdnmgHere->choice!=pnnHere->res1) {
      pdnmgHere=pdnmgHere->next;
      if (pdnmgHere==NULL) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,2,21,"Res %d not found in this structure!",pnnHere->res1);
		return;
      }
    }
    pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
    pvnmaHere=pmgdHere->pvnmaAHead;
    NameThere=pnnHere->AtomName1;
    do {
      NameHere=((PMAD)(pvnmaHere->data.ptrvalue))->pcAName;
      if (!StringCmp(NameThere,NameHere)) break;
      /* take care of atom degeneracies */
      /* handle by pointing pmad1, pmad2 to first acceptable match found */
      if ((NameThere[0]=='*') && (NameThere[3]=='*')) {
				if (NameThere[1]==NameHere[1] && NameThere[2]==NameHere[2])
	  			break;
      }
      else if (NameThere[0]=='*') {
				if (NameThere[1]==NameHere[1] && NameThere[2]==NameHere[2] && NameThere[3]==NameHere[3])
	  			break;
      }
      else if (NameThere[3]=='*') {
				if (NameThere[0]==NameHere[0] && NameThere[1]==NameHere[1] && NameThere[2]==NameHere[2])
	  			break;
      }
      pvnmaHere=pvnmaHere->next;
      if (pvnmaHere==NULL) {
				PurgeGlobs();
				ErrPostEx(SEV_FATAL,2,23,"Atomic constraint res %d %s doesn't correspond to this structure!",pnnHere->res1,pnnHere->AtomName1);
				return;
      }
    } while (1);
    pnnHere->pmad1=(PMAD)(pvnmaHere->data.ptrvalue);
    /* do second atom now, note res2>=res1 always */
    while (pdnmgHere->choice!=pnnHere->res2) {
      pdnmgHere=pdnmgHere->next;
      if (pdnmgHere==NULL) {
				PurgeGlobs();
				ErrPostEx(SEV_FATAL,2,24,"Res %d not found in this structure!",pnnHere->res2);
				return;
      }
    }
    pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
    pvnmaHere=pmgdHere->pvnmaAHead;
    NameThere=pnnHere->AtomName2;
    do {
      NameHere=((PMAD)(pvnmaHere->data.ptrvalue))->pcAName;
      if (!StringCmp(NameThere,NameHere)) break;
      /* take care of atom degeneracies */
      /* handle by pointing pmad1, pmad2 to first acceptable match found */
      if ((NameThere[0]=='*') && (NameThere[3]=='*')) {
				if (NameThere[1]==NameHere[1] && NameThere[2]==NameHere[2])
				  break;
      }
      else if (NameThere[0]=='*') {
				if (NameThere[1]==NameHere[1] && NameThere[2]==NameHere[2] && NameThere[3]==NameHere[3])
				  break;
      }
      else if (NameThere[3]=='*') {
				if (NameThere[0]==NameHere[0] && NameThere[1]==NameHere[1] && NameThere[2]==NameHere[2])
				  break;
      }
      pvnmaHere=pvnmaHere->next;
      if (pvnmaHere==NULL) {
				PurgeGlobs();
				ErrPostEx(SEV_FATAL,2,23,"Atomic constraint res %d %s doesn't correspond to this structure!",pnnHere->res2,pnnHere->AtomName2);
				return;
      }
    } while (1);
    pnnHere->pmad2=(PMAD)(pvnmaHere->data.ptrvalue);
    if ((pnnHere->AtomName1)[1]=='S' && (pnnHere->AtomName2)[1]=='S') {
      /* make a disulfide bridge in chemical graph */
      pmmdParent=GetParentMol((PFB)(pnnHere->pmad1));
      (pmmdParent->iIRBCount)++;
      pvnmbNew=NewVNMB(&(pmmdParent->pvnmbIRBHead),Bond_order_single);
      pmbdNew=(PMBD)(pvnmbNew->data.ptrvalue);
      pmbdNew->pmadFrom=pnnHere->pmad1;
      pmbdNew->pmadTo=pnnHere->pmad2;
      pmbdNew->pfbParent=(PFB)pmmdParent;
      pmbdNew->bMe=(Byte)AM_MBD;
      pmbdNew->bWhat=(Byte)BOND_SINGLE;
      ValNodeAddPointer(&((pnnHere->pmad1)->pvnBonds),Bond_order_single,(Nlm_VoidPtr)pmbdNew);
      ValNodeAddPointer(&((pnnHere->pmad2)->pvnBonds),Bond_order_single,(Nlm_VoidPtr)pmbdNew);
    }
    pnnHere=pnnHere->next;
  }
}

void PrintDistConstTries(void)
{
  PNN pnnhere;

  pnnhere=pnnDistConstHead; 
  while(pnnhere) {
    printf("%4d %4s %4d %4s %5.3f %5.3f %5.3f %5.3f %4ld\n",pnnhere->res1,pnnhere->AtomName1,pnnhere->res2,pnnhere->AtomName2,pnnhere->MeanDist,pnnhere->MinDelta,pnnhere->MaxDelta,pnnhere->prob,(long  int)pnnhere->tries);
    pnnhere=pnnhere->next;
  }
}

void PrintDistConstraints(void)
{
  PNN pnnhere;

  pnnhere=pnnDistConstHead;
  while(pnnhere) {
    printf("%4d %4s %4d %4s %5.3f %5.3f %5.3f %5.3f\n",pnnhere->res1,pnnhere->AtomName1,pnnhere->res2,pnnhere->AtomName2,pnnhere->MeanDist,pnnhere->MinDelta,pnnhere->MaxDelta,pnnhere->prob);
    printf("\t%5.3f %5.3f %5.3f %5.3f %5.3f\n",pnnhere->Angle1,pnnhere->Angle2,pnnhere->Dihedral01,pnnhere->Dihedral12,pnnhere->Dihedral23);
    pnnhere=pnnhere->next;
  }
}

void PrintTrueDistConstraints(void)
{
  PNN pnnhere;

  pnnhere=pnnDistConstTrueHead;
  while(pnnhere) {
    printf("%4d %4s %4d %4s %5.3f %5.3f %5.3f %5.3f\n",pnnhere->res1,pnnhere->AtomName1,pnnhere->res2,pnnhere->AtomName2,pnnhere->MeanDist,pnnhere->MinDelta,pnnhere->MaxDelta,pnnhere->prob);
    printf("\t%5.3f %5.3f %5.3f %5.3f %5.3f\n",pnnhere->Angle1,pnnhere->Angle2,pnnhere->Dihedral01,pnnhere->Dihedral12,pnnhere->Dihedral23);
    pnnhere=pnnhere->next;
  }
}

/* in case it is needed - but use as read only please! */
PNN GetTrueDistConstraints(void)
{
  return pnnDistConstTrueHead;
}

PNN DeleteDistConstraint(PNN pnnNode)
{
	if (pnnNode==pnnDistConstTrueHead) {
		pnnDistConstTrueHead=pnnNode->next;
		if (pnnDistConstTrueHead!=NULL)
			pnnDistConstTrueHead->prev=NULL;
	}
	else {
		pnnNode->prev->next=pnnNode->next;
		if (pnnNode->next!=NULL)
			pnnNode->next->prev=pnnNode->prev;
	}
	MemFree(pnnNode);
	return NULL;
}

TrajErr AddDistConstraint(PNN pnnNode)
{
  Int2 r1,r2,endlist;
  PNN pnnHere;

  if (pnnNode==NULL) return ERR_FAIL;
  r1=pnnNode->res1;
  r2=pnnNode->res2;
  if ((r1<0) || (r2<0)) return ERR_FAIL;
  if (r1>r2) return ERR_FAIL;
  endlist=0;
  pnnHere=pnnDistConstTrueHead;
  if (pnnHere==NULL) {
    pnnDistConstTrueHead=pnnNode;
    pnnNode->next=NULL;
    pnnNode->prev=NULL;
    return ERR_SUCCESS;
  }
  while ((!endlist) && (pnnHere->res1<r1)) {
    if (pnnHere->next==NULL)
      endlist=1;
    else
      pnnHere=pnnHere->next;
  }
  while ((!endlist) && (pnnHere->res1==r1) && (pnnHere->res2<r2)) {
    if (pnnHere->next==NULL)
      endlist=1;
    else
      pnnHere=pnnHere->next;
  }
  /* append to end of list */
  if (endlist) {
    pnnHere->next=pnnNode;
    pnnNode->prev=pnnHere;
    pnnNode->next=NULL;
    return ERR_SUCCESS;
  }
  if (pnnHere==pnnDistConstTrueHead) {  /* && endlist == 0 */
    /* adding to head of list */
    pnnHere->prev=pnnNode;
    pnnNode->next=pnnHere;
    pnnNode->prev=NULL;
    pnnDistConstTrueHead=pnnNode;
    return ERR_SUCCESS;
  }
  /* insert between pnnHere and pnnHere->prev */
  pnnNode->next=pnnHere;
  pnnNode->prev=pnnHere->prev;
  (pnnHere->prev)->next=pnnNode;
  pnnHere->prev=pnnNode;
  return ERR_SUCCESS;
}

/* CodeBase database variables */
static CODE4 trajcb;       /* main database structure */
/* fields in each trajectory graph record (1 per residue) */
static FIELD4INFO f4iTrajDB[]=
{
  {"RES",r4num,6,0},        /* Residue number (primary key) */
  {"TG",r4memo,10,0},       /* actual trajectory graph binary data */
  {"INTEGRAL",r4num,10,0},  /* Trajectory Graph integral */
  {"CTYPE",r4num,6,0},      /* compression type used */
  {"BUFSZ",r4num,10,0},     /* BZ buffer size */
  {"PEAK",r4num,10,0},      /* Trajectory Graph largest value */
  {"PCIS",r4float,7,4},     /* probabaility residue is cis- (%) */
  {"CISTG",r4memo,10,0},    /* actual cis-trajectory graph binary data */
  {"CISINT",r4num,10,0},    /* cis-Trajectory Graph integral */
  {"CISBUFSZ",r4num,10,0},  /* cis-BZ buffer size */
  {"CISPEAK",r4num,10,0},   /* cis-Trajectory Graph largest value */
  {"PSS",r4float,7,4},      /* probabaility residue is S-S (%) */
  {"CHIWMEAN",r4float,7,2}, /* average value of omega for this residue (if cis, 180 is added) */
  {"CHIWSD",r4float,5,2},   /* standard deviation for omega (between residue i */
			    /* and i-1) */
  {"AA",r4str,1,0},         /* 1-letter code for amino acid */
  {"DIM",r4num,6,0},        /* equals TRAJDIV, dimension of TG */
  {"NZ1",r4num,6,0},        /* first non-zero row in TG */
  {"NUMNZ",r4num,6,0},      /* number of non-zero rows in TG from NZ1 to end */
  /* these four indicate sparseness of the TG, how close to done we are */
  {"NELT0",r4num,10,0},     /* # elements of graph =< 0 */
  {"NELT5",r4num,10,0},     /* # elements of graph =< 5% max */
  {"NELT10",r4num,10,0},    /* # elements of graph =< 10% max */
  {"NELT15",r4num,10,0},    /* # elements of graph =< 15% max */
  {"TOUT",r4num,6,0},	    /* timeout value */
  {"ENDIAN",r4num,6,0},     /* endianness flag */
  {"MARKOVSF",r4float,7,2}, /* Markovian scale factor */
  {"ROT",r4num,10,0},       /* Rotamer id (?) */
  {"NUMFRAG",r4num,6,0},       /* # fragments starting at this residue */
  {"FRAG",r4memo,10,0},       /* actual fragment binary data */
  {0,0,0,0}
};
/* searchable tags */
static TAG4INFO t4iTrajDB[]=
{
  {"RES","RES",0,r4unique,0},    /* Index by primary key */
  {0,0,0,0,0}
};
/* Data file pointer */
static DATA4 *pd4Traj=NULL;
/* Field pointers */
static FIELD4 *pf4res=NULL;
static FIELD4 *pf4comptype=NULL;
static FIELD4 *pf4tg=NULL;
static FIELD4 *pf4integral=NULL;
static FIELD4 *pf4BZbufsize=NULL;
static FIELD4 *pf4peak=NULL;
static FIELD4 *pf4pcis=NULL;
static FIELD4 *pf4cistg=NULL;
static FIELD4 *pf4cisintegral=NULL;
static FIELD4 *pf4cisBZbufsize=NULL;
static FIELD4 *pf4cispeak=NULL;
static FIELD4 *pf4pss=NULL;
static FIELD4 *pf4chiwmean=NULL;
static FIELD4 *pf4chiwsd=NULL;
static FIELD4 *pf4aaname=NULL;
static FIELD4 *pf4trajdiv=NULL;
static FIELD4 *pf4firstnzrow=NULL;
static FIELD4 *pf4numnzrows=NULL;
static FIELD4 *pf4nelt0=NULL;
static FIELD4 *pf4nelt5=NULL;
static FIELD4 *pf4nelt10=NULL;
static FIELD4 *pf4nelt15=NULL;
static FIELD4 *pf4tout=NULL;
static FIELD4 *pf4markovsf=NULL;
static FIELD4 *pf4endian=NULL;
static FIELD4 *pf4rot=NULL;
static FIELD4 *pf4numfrag=NULL;
static FIELD4 *pf4frag=NULL;
/* Tag pointer */
static TAG4 *pt4res = NULL;

/* performs a little endian to big endian flip on the bytes pointed to by
   buf; the size, bufsize, is given in Int4's, so this procedure assumes
   the buffer is really a series of bufsize Int4's, only cast as a CharPtr
   for the purpose of byte flipping */
void FlipBytes(CharPtr buf,Int4 bufsize)
{
  Int4 cnt;
  Char c1,c2,c3,c4;
  CharPtr bufhere;

  bufhere=buf;
  for (cnt=0;cnt<bufsize;cnt++) {
    c1=bufhere[3];
    c2=bufhere[2];
    c3=bufhere[1];
    c4=bufhere[0];
    bufhere[0]=c1;
    bufhere[1]=c2;
    bufhere[2]=c3;
    bufhere[3]=c4;
    bufhere+=4;
  }
}

TrajErr TGInit(CharPtr dbfnam,Int2 makeit,Int4 PNTR numrec) {
  /* initialize the database */
  Int2 cnt;
  Int4 reccount;
  INDEX4 *idx4;
  Boolean newdb=FALSE;

  if (code4init(&trajcb) != r4success)              /* cb is a global CODE4 structure */
    return ERR_FAIL;
  trajcb.accessMode         = OPEN4DENY_WRITE;      /* deny simultaneous writing to data file */
  trajcb.errDefaultUnique   = r4unique;             /* disallow duplicate keys & generate error */
  trajcb.errOpen            = 0;                    /* generate error if file cannot be opened */
  trajcb.optimize           = OPT4OFF;              /* optimize file when opened/created */
  trajcb.optimizeWrite      = OPT4OFF;              /* optimize file when writing */
  trajcb.autoOpen			= 0;					/* do not automatically load index file on d4open */
  trajcb.errOff				= 1;					/* do not output error messages to screen */
  /* now open the database file */
  pd4Traj=d4open(&trajcb,dbfnam);
  if (pd4Traj==NULL) {
	  if (makeit==DB_CREATE) {
		pd4Traj=d4create(&trajcb,dbfnam,f4iTrajDB,t4iTrajDB);
		newdb=TRUE;
	  }
  }
  if (pd4Traj==NULL) {
	ErrPostEx(SEV_ERROR,10,0,"Could not open/create data file %s",dbfnam);
    PurgeGlobs();
	ErrPostEx(SEV_FATAL,10,0,"Unable to continue");
    return ERR_FAIL;
  }
  if (newdb==FALSE) {
	  idx4=i4open(pd4Traj,NULL);
	  if (idx4==NULL) {
		  if (trajcb.errorCode==-64 || trajcb.errorCode==-490) {
			d4close(pd4Traj);
			trajcb.errorCode=0;
			trajcb.safety=0;
			trajcb.errOff=0;
			trajcb.accessMode=OPEN4DENY_RW;
			pd4Traj=d4open(&trajcb,dbfnam);
			if (i4create(pd4Traj,NULL,t4iTrajDB)==NULL) {
				ErrPostEx(SEV_ERROR,trajcb.errorCode,0,"Problem while reindexing data file %s",dbfnam);
				PurgeGlobs();
				ErrPostEx(SEV_FATAL,10,2,"Unable to continue");
				return ERR_FAIL;
			}
			d4close(pd4Traj);
			trajcb.accessMode=OPEN4DENY_WRITE;
			trajcb.safety=1;
			trajcb.autoOpen=1;
			pd4Traj=d4open(&trajcb,dbfnam);
			if (pd4Traj==NULL) {
				ErrPostEx(SEV_ERROR,trajcb.errorCode,0,"Reindexing failed for data file %s",dbfnam);
				PurgeGlobs();
				ErrPostEx(SEV_FATAL,10,3,"Unable to continue");
				return ERR_FAIL;
			}
			ErrPostEx(SEV_ERROR,1,0,"Warning, corrupt index detected in data file %s, file has been reindexed successfully",dbfnam);
		  }
		  else {
			ErrPostEx(SEV_ERROR,trajcb.errorCode,0,"Could not open/create index for data file %s",dbfnam);
			PurgeGlobs();
			ErrPostEx(SEV_FATAL,10,1,"Unable to continue");
			return ERR_FAIL;	
		  }
	  }
  }
  trajcb.errOff=0;
  reccount=(Int4)d4recCount(pd4Traj);
  if (makeit==DB_READ) {
    if (numrec!=NULL) {
      *numrec=reccount;
    }
  }
  /* Assign global field/tag pointers to database fields/tags */   
  pf4res=d4field(pd4Traj,"RES");
  pf4tg=d4field(pd4Traj,"TG");
  pf4integral=d4field(pd4Traj,"INTEGRAL");
  pf4comptype=d4field(pd4Traj,"CTYPE");
  pf4BZbufsize=d4field(pd4Traj,"BUFSZ");
  pf4peak=d4field(pd4Traj,"PEAK");
  pf4pcis=d4field(pd4Traj,"PCIS");
  pf4cistg=d4field(pd4Traj,"CISTG");
  pf4cisintegral=d4field(pd4Traj,"CISINT");
  pf4cisBZbufsize=d4field(pd4Traj,"CISBUFSZ");
  pf4cispeak=d4field(pd4Traj,"CISPEAK");
  pf4pss=d4field(pd4Traj,"PSS");
  pf4chiwmean=d4field(pd4Traj,"CHIWMEAN");
  pf4chiwsd=d4field(pd4Traj,"CHIWSD");
  pf4aaname=d4field(pd4Traj,"AA");
  pf4trajdiv=d4field(pd4Traj,"DIM");
  pf4firstnzrow=d4field(pd4Traj,"NZ1");
  pf4numnzrows=d4field(pd4Traj,"NUMNZ");
  pf4nelt0=d4field(pd4Traj,"NELT0");
  pf4nelt5=d4field(pd4Traj,"NELT5");
  pf4nelt10=d4field(pd4Traj,"NELT10");
  pf4nelt15=d4field(pd4Traj,"NELT15");
  pf4tout=d4field(pd4Traj,"TOUT");
  pf4markovsf=d4field(pd4Traj,"MARKOVSF");
  pf4endian=d4field(pd4Traj,"ENDIAN");
  pf4rot=d4field(pd4Traj,"ROT");
  pf4numfrag=d4field(pd4Traj,"NUMFRAG");
  pf4frag=d4field(pd4Traj,"FRAG");
  pt4res=d4tag(pd4Traj,"RES");
  /* select default tag to use */
  d4tagSelect(pd4Traj,pt4res);
  if (USE_LOTS_RAM) {
	/* load all trajectories in memory */
	ptgsarrsize=reccount;
	ptgsarr=MemNew(ptgsarrsize*sizeof(PTGS));
	for (cnt=1;cnt<=ptgsarrsize;cnt++)
		ptgsarr[cnt-1]=TrajGraphRead(cnt);
	largeram_loaded=TRUE;
  }
  return ERR_SUCCESS;
}

void TGClose(void) {
	Int2 cnt;

	/* Close the trajectory graph database, free memory */
	if (USE_LOTS_RAM) {
		largeram_loaded=FALSE;
		for (cnt=0;cnt<ptgsarrsize;cnt++)
			ptgsarr[cnt]=FreeTraj(ptgsarr[cnt]);
		ptgsarr=MemFree(ptgsarr);
	}
	code4close(&trajcb);
	code4initUndo(&trajcb);
}

void TGPack(void) {
	Int4 numrecords;
	Boolean ulr;

	ulr=USE_LOTS_RAM;
	USE_LOTS_RAM=FALSE;
	TGInit(tmpdbasename,DB_READ,&numrecords);
	d4pack(pd4Traj);
	d4memoCompress(pd4Traj);
	TGClose();
	USE_LOTS_RAM=ulr;
}

Int4 RLEPack(Uint4 *dest,Uint4 *src,Int4 numbytes)
{
	Uint4 runlength,srcpos,destpos;
	Uint4 lastbyte,nextbyte;

	lastbyte=src[0];
	runlength=1;
	srcpos=1;
	destpos=0;
	do {
		nextbyte=src[srcpos];
		if ((nextbyte!=lastbyte) || (runlength>RLEMAX-1UL)) {
			if ((runlength==1) && (lastbyte<0xFFFFFFFFUL-RLEMAX+1UL)) {
				dest[destpos]=lastbyte;
				destpos++;
			}
			else {
				dest[destpos]=0xFFFFFFFFUL-RLEMAX+runlength;
				dest[destpos+1]=lastbyte;
				destpos+=2;
			}
			lastbyte=nextbyte;
			runlength=1;
		}
		else {
			runlength++;
		}
		srcpos++;
	} while (srcpos<numbytes);
	if ((runlength==1) && (lastbyte<0xFFFFFFFFUL-RLEMAX+1UL)) {
		dest[destpos]=lastbyte;
		destpos++;
	}
	else {
		dest[destpos]=0xFFFFFFFFUL-RLEMAX+runlength;
		dest[destpos+1]=lastbyte;
		destpos+=2;
	}
	return destpos;
}

Int4 RLEUnPack(Uint4 *dest,Uint4 *src,Int4 numbytes,Int4 data_endian)
{
	Uint4 *srcpos,*destpos,runlength;
	Uint4 checkbyte,databyte;
	Int4 mustflip=0;
	Uint4 *uihere,*uiend,*bufend;

	if (GetEndian()==ENDIAN_UNKNOWN)
		return 0;
	if (GetEndian()!=data_endian) mustflip=1;
	if (numbytes==0)
		return 0;
	srcpos=src;
	destpos=dest;
	bufend=src+numbytes;
	do {
		checkbyte=*srcpos;
		if (mustflip)
			FlipBytes((CharPtr)(&checkbyte),1);
		if (checkbyte<0xFFFFFFFFUL-RLEMAX+1UL) {
			runlength=1;
			databyte=checkbyte;
			srcpos++;
		}
		else {
			runlength=checkbyte-0xFFFFFFFFUL+RLEMAX;
			databyte=*(srcpos+1);
			if (mustflip)
				FlipBytes((CharPtr)(&databyte),1);
			srcpos+=2;
		}
		/* an attempt to optimize speed here */
		uihere=destpos;
		uiend=destpos+runlength;
		for (;uihere<uiend;uihere++)
			*uihere=databyte;
		destpos+=runlength;
	} while (srcpos<bufend);
	return (Int4)(destpos-dest);
}



void UpdatePCis(Int4 resnum,FloatLo cis)
{
	d4seekDouble(pd4Traj,(double)resnum);
	f4assignDouble(pf4pcis,(double)cis);
}

FloatLo GetPCis(Int4 resnum)
{
	d4seekDouble(pd4Traj,(double)resnum);
	return (FloatLo)(f4double(pf4pcis));
}

TrajErr FlipFDS(PFDS pfds)
{
	if (sizeof(FloatLo)!=4) {
		ErrPostEx(SEV_ERROR,8,0,"Sizeof FloatLo is unexpected, cannot continue");
		return ERR_FAIL;
	}
	if (sizeof(Int4)!=4) {
		ErrPostEx(SEV_ERROR,8,0,"Sizeof Int4 is unexpected, cannot continue");
		return ERR_FAIL;
	}
	pfds->resnum=FlipI4(&(pfds->resnum));
	pfds->pSS=FlipF(&(pfds->pSS));
	pfds->Angle1=FlipF(&(pfds->Angle1));
	pfds->Angle2=FlipF(&(pfds->Angle2));
	pfds->AngleSD=FlipF(&(pfds->AngleSD));
	pfds->ChiWMean=FlipF(&(pfds->ChiWMean));
	pfds->ChiWSD=FlipF(&(pfds->ChiWSD));
	pfds->prob=FlipF(&(pfds->prob));
	pfds->tout=FlipI4(&(pfds->tout));
	pfds->rotid=FlipUI4(&(pfds->rotid));
	pfds->length=FlipI4(&(pfds->length));
	pfds->reserved=FlipI4(&(pfds->reserved));
	return ERR_SUCCESS;
}

TrajErr TrajGraphWrite(PTGS ptgsThis,Int4 ctype,Boolean Replace)
{
    Int2 rec;
    Int4 err,totalwritten=0,numfrag;
    CharPtr FromHere;
    CharPtr bzBuffer=NULL,RLEBuffer=NULL;
    Uint4 rlebufsize;
    unsigned int actbufsize;
    PFDS pfdsHere,pfdsNextFrag;
    CharPtr pcHere;
    Boolean isnewrow;

    if (!Replace) {
	    /* go to end of current file */
	    d4bottom(pd4Traj);
	    /* Append a new record to end of datafile & make it blank */
	    d4appendStart(pd4Traj,0);
	    d4blank(pd4Traj);
    }
    else {
	    /* replace current record */
	    rec=d4seekDouble(pd4Traj,(double)(ptgsThis->resnum));
	    if (rec!=r4success) {
	        PurgeGlobs();
  	        ErrPostEx(SEV_FATAL,12,4,"Attempt to replace non-existent residue number into database");
  	        return ERR_FAIL;
  	    }
    }
    /* Assign values to the new record */
    f4assignInt(pf4res,ptgsThis->resnum);
    FromHere=(CharPtr)(ptgsThis->TrajGraph)+(ptgsThis->dim)*((ptgsThis->firstnzrow)-1)*sizeof(Int4);
    
    /* use native endian format for data when writing */
/*  if (GetEndian()==ENDIAN_BIG)
    FlipBytes(FromHere,(Int4)((ptgsThis->dim)*(ptgsThis->numnzrows)));*/
    if (ctype==USE_BZ) {
        /* worst case buffer size allocated (from bzip2 manual) */
        actbufsize=(unsigned int)(sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim)*1.01)+600;
        bzBuffer=(CharPtr)MemNew((size_t)actbufsize);
        if ((err=BZ2_bzBuffToBuffCompress((char *)bzBuffer,&actbufsize,(char *)FromHere,(unsigned int)(sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->numnzrows)),9,0,50))!=BZ_OK) {
            PurgeGlobs();
            ErrPostEx(SEV_FATAL,err,0,"Bzip compression fatal error occurred");
        }
        totalwritten=(Int4)actbufsize;
        
        /* Write to codebase */
        if ( (f4memoAssignN(pf4tg,bzBuffer,totalwritten)) != r4success ){
            PurgeGlobs();
  	        ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	        return ERR_FAIL;
        }
    }
    
    if (ctype==USE_RLE) {
        rlebufsize=sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim)*2;
        RLEBuffer=(CharPtr)MemNew(rlebufsize);
        totalwritten=sizeof(Int4)*RLEPack((Uint4 *)RLEBuffer,(Uint4 *)FromHere,(ptgsThis->dim)*(ptgsThis->numnzrows));
        if ( (f4memoAssignN(pf4tg,RLEBuffer,totalwritten)) != r4success ){
            PurgeGlobs();
  	        ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	        return ERR_FAIL;
        }
    }
    if (ctype==USE_NONE) {
        totalwritten=sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim);
        if ( (f4memoAssignN(pf4tg,FromHere,totalwritten)) != r4success ){
            PurgeGlobs();
  	        ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	        return ERR_FAIL;
        }
    }
    
    /* flip back after writing in case value is still used later on */
/*  if (GetEndian()==ENDIAN_BIG)
    FlipBytes(FromHere,(Int4)((ptgsThis->dim)*(ptgsThis->numnzrows)));*/
  
    f4assignLong(pf4comptype,ctype);
    f4assignLong(pf4BZbufsize,totalwritten);
    f4assignLong(pf4integral,ptgsThis->TrajIntegral);
    f4assignLong(pf4peak,ptgsThis->Peak);
    f4assignDouble(pf4pcis,(double)(ptgsThis->pCis));
    totalwritten=0;
  
    if ((ptgsThis->CisTrajGraph)!=NULL) {
        /*if (GetEndian()==ENDIAN_BIG)
        FlipBytes((CharPtr)(ptgsThis->CisTrajGraph),(Int4)((ptgsThis->dim)*(ptgsThis->dim)));*/
        if (ctype==USE_BZ) {
            /* worst case buffer size allocated (from bzip2 manual) */
            /* reuse same bz buffer */
            actbufsize=((unsigned int)(sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim)*1.01))+600;
            if ((err=BZ2_bzBuffToBuffCompress((char *)bzBuffer,&actbufsize,(char *)(ptgsThis->CisTrajGraph),(unsigned int)(sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim)),9,0,50))!=BZ_OK) {
                PurgeGlobs();
                ErrPostEx(SEV_FATAL,err,1,"Bzip compression fatal error occurred");
            }
            totalwritten=(Int4)actbufsize;
            if ( (f4memoAssignN(pf4cistg,bzBuffer,totalwritten)) != r4success ){
                PurgeGlobs();
  	            ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	            return ERR_FAIL;
            }
        }
        if (ctype==USE_RLE) {
            rlebufsize=sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim)*2;
            totalwritten=sizeof(Int4)*RLEPack((Uint4 *)RLEBuffer,(Uint4 *)(ptgsThis->CisTrajGraph),(ptgsThis->dim)*(ptgsThis->dim));
            if ( (f4memoAssignN(pf4cistg,RLEBuffer,totalwritten)) != r4success ){
                PurgeGlobs();
  	            ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	            return ERR_FAIL;
            }
        }
        if (ctype==USE_NONE) {
            totalwritten=sizeof(Int4)*(ptgsThis->dim)*(ptgsThis->dim);
            if ( (f4memoAssignN(pf4cistg,(CharPtr)(ptgsThis->CisTrajGraph),totalwritten)) != r4success ){
                PurgeGlobs();
  	            ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	            return ERR_FAIL;
            }
        }
        
        /* put memory back the way it was to be safe */
    /*  if (GetEndian()==ENDIAN_BIG)
        FlipBytes((CharPtr)(ptgsThis->CisTrajGraph),(Int4)((ptgsThis->dim)*(ptgsThis->dim)));*/
    }
    if (ctype==USE_BZ)
        bzBuffer=MemFree(bzBuffer);
    if (ctype==USE_RLE)
        RLEBuffer=MemFree(RLEBuffer);
        
    f4assignLong(pf4cisBZbufsize,totalwritten);
    f4assignLong(pf4cisintegral,ptgsThis->CisTrajIntegral);
    f4assignLong(pf4cispeak,ptgsThis->CisPeak);
    f4assignDouble(pf4pss,(double)(ptgsThis->pSS));
    f4assignDouble(pf4chiwmean,(double)(ptgsThis->ChiWMean));
    f4assignDouble(pf4chiwsd,(double)(ptgsThis->ChiWSD));
    f4assignChar(pf4aaname,ptgsThis->AA);
    f4assignInt(pf4trajdiv,ptgsThis->dim);
    f4assignInt(pf4firstnzrow,ptgsThis->firstnzrow);
    f4assignInt(pf4numnzrows,ptgsThis->numnzrows);
    f4assignLong(pf4nelt0,ptgsThis->nelt0);
    f4assignLong(pf4nelt5,ptgsThis->nelt5);
    f4assignLong(pf4nelt10,ptgsThis->nelt10);
    f4assignLong(pf4nelt15,ptgsThis->nelt15);
    f4assignLong(pf4tout,ptgsThis->tout);
    f4assignDouble(pf4markovsf,(double)(ptgsThis->markovsf));
    f4assignLong(pf4endian,GetEndian());
    f4assignLong(pf4rot,ptgsThis->rotid);
  
    /* count fragments if any */
    numfrag=0;
    pfdsHere=ptgsThis->pfdsFragmentHead;
    if (pfdsHere!=NULL) do {
	    pfdsNextFrag=pfdsHere->nextfrag;
	    
	    while (pfdsHere!=NULL) {
  		    pfdsHere=pfdsHere->next;
	  	    numfrag++;	  	
	    }
        pfdsHere=pfdsNextFrag;
        
	} while (pfdsHere!=NULL);
    
    f4assignLong(pf4numfrag,numfrag);
    
    /* build up array of fragments in memory to write as BLOB */
    if (numfrag) {
        RLEBuffer=(CharPtr)MemNew(sizeof(FDTS)*numfrag);
	    pcHere=RLEBuffer;
	    pfdsHere=ptgsThis->pfdsFragmentHead;
	    do {
            pfdsNextFrag=pfdsHere->nextfrag;
            isnewrow=TRUE;
            while (pfdsHere!=NULL) {
			    MemCopy(pcHere,(VoidPtr)pfdsHere,sizeof(FDTS));
			    if (isnewrow) {
				    ((PFDTS)pcHere)->newrow=1;
				    isnewrow=FALSE;
			    }
			    else
				    ((PFDTS)pcHere)->newrow=0;
			    pcHere+=sizeof(FDTS);
	  		    pfdsHere=pfdsHere->next;
		    }	
		    pfdsHere=pfdsNextFrag;
	    } while (pfdsHere!=NULL);
	    
        if ( (f4memoAssignN(pf4frag,RLEBuffer,sizeof(FDTS)*numfrag)) != r4success ){
            PurgeGlobs();
  	        ErrPostEx(SEV_FATAL,errno,0,"Unable to write to database");
  	        return ERR_FAIL;
        }
        RLEBuffer=MemFree(RLEBuffer);	
	}
	
    /* Write new record to datafile */
    if (!Replace) {
	    rec=d4append(pd4Traj);
	    if (rec==r4success)
	        return ERR_SUCCESS;
	    if (rec==r4unique) {
	        PurgeGlobs();
	        ErrPostEx(SEV_FATAL,12,2,"Attempt to insert duplicate residue number into database");
	        return ERR_FAIL;
	    }
	    if (rec<0) {
	        PurgeGlobs();
	        ErrPostEx(SEV_FATAL,13,3,"Problem adding new record, CodeBase error %d",rec);
	        return ERR_FAIL;
	    }
    }
    
    return ERR_SUCCESS;
}

/* set cis=0 if no trajectory graph for Cis- possibility is required */
PTGS NewTraj(Int2 cis,Int2 graphwidth)
{
  PTGS ptgsHere;

  ptgsHere=(PTGS)MemNew(sizeof(TGS));
  ptgsHere->TrajGraph=(Int4 *)MemNew(graphwidth*graphwidth*sizeof(Int4));
  if (cis)
    ptgsHere->CisTrajGraph=(Int4 *)MemNew(graphwidth*graphwidth*sizeof(Int4));
  else
    ptgsHere->CisTrajGraph=NULL;
  return ptgsHere;
}

PFDS FreeFragmentList(PFDS pfdsHead)
{
	PFDS pfdsNext,pfdsHere,pfdsNextFrag;
	
	if (pfdsHead==NULL)
		return NULL;
	pfdsHere=pfdsHead;
	do {
		pfdsNextFrag=pfdsHere->nextfrag;
		while (pfdsHere!=NULL) {
			pfdsNext=pfdsHere->next;
			MemFree(pfdsHere);
			pfdsHere=pfdsNext;
		}
		pfdsHere=pfdsNextFrag;
	} while (pfdsHere!=NULL);
	return NULL;
}

PTGS FreeTraj(PTGS ptgsThis)
{
  if (largeram_loaded==TRUE)
	  return NULL;
  ptgsThis->TrajGraph=MemFree(ptgsThis->TrajGraph);
  if ((ptgsThis->CisTrajGraph)!=NULL)
    ptgsThis->CisTrajGraph=MemFree(ptgsThis->CisTrajGraph);
  if (ptgsThis->pfdsFragmentHead!=NULL)
  	FreeFragmentList(ptgsThis->pfdsFragmentHead);
  ptgsThis=MemFree(ptgsThis);
  return NULL;
}

PTGS TrajGraphRead(Int2 resnum)
{
  Int2 rec,graphwidth;
  PTGS ptgsHere=NULL;
  Char aahere;
  CharPtr FillHere;
  Int4 bzbufsize,ctype,data_endian,err,numfrag;
  unsigned int actbufsize;
  CharPtr pcHere=NULL;
  PFDS pfdsLast=NULL,pfdsHere,pfdsLastFrag;
 	Boolean isnewrow;
 	
  if (largeram_loaded)
	  return ptgsarr[resnum-1];
  d4top(pd4Traj);
  rec=d4seekDouble(pd4Traj,(double)resnum);    
  if (rec!=r4success) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,14,4,"Unable to find trajectory distribution for residue %d,d4seek return code %d",resnum,rec);
    return NULL;
  }
  /* if residue is Pro, may need cis-TG at this residue */
  graphwidth=f4int(pf4trajdiv);
  aahere=(Char)f4char(pf4aaname);
  if (aahere=='P')
    ptgsHere=NewTraj(1,graphwidth);
  else
    ptgsHere=NewTraj(0,graphwidth);
  ptgsHere->resnum=f4int(pf4res);
  ptgsHere->dim=graphwidth;
  ctype=f4int(pf4comptype);
  ptgsHere->firstnzrow=f4int(pf4firstnzrow);
  ptgsHere->numnzrows=f4int(pf4numnzrows);
  data_endian=f4long(pf4endian);
  FillHere=(CharPtr)(ptgsHere->TrajGraph);
  /* since ptgsHere->TrajGraph was allocated with MemNew, it should
     contain all zeroes now; thus we need only load in non-zero
     values */
  if ((ptgsHere->firstnzrow)>1)
    FillHere=FillHere+(ptgsHere->dim)*((ptgsHere->firstnzrow)-1)*sizeof(Int4);
  bzbufsize=f4long(pf4BZbufsize);    
  actbufsize=(unsigned int)(ptgsHere->dim)*(ptgsHere->numnzrows)*sizeof(Int4);
  if (ctype==USE_BZ) {
    if ((err=BZ2_bzBuffToBuffDecompress((char *)FillHere,&actbufsize,(char *)f4memoPtr(pf4tg),(unsigned int)bzbufsize,0,0))!=BZ_OK) { 
      PurgeGlobs();
      ErrPostEx(SEV_FATAL,err,0,"Bzip decompression error occurred");
    }
    if (GetEndian()==ENDIAN_UNKNOWN)
	ErrPostEx(SEV_FATAL,1,24,"Unknown endian system, cannot read in Trajectory Distribution");
    if (GetEndian()!=data_endian)
      FlipBytes((CharPtr)(ptgsHere->TrajGraph),(Int4)((ptgsHere->dim)*(ptgsHere->dim)));
  }
  if (ctype==USE_RLE) {
    /* bzbufsize holds num chars in RLE file written - flip bytes first for RLE */
    if ((err=RLEUnPack((Uint4 *)FillHere,(Uint4 *)f4memoPtr(pf4tg),bzbufsize/sizeof(Int4),data_endian))!=ptgsHere->dim*ptgsHere->numnzrows)
	ErrPostEx(SEV_FATAL,23,24,"RLEUnPack failed, size=%ld, should=%ld*%ld - likely this is caused by overclocked or faulty RAM chips, please test your RAM",(long)err,(long)ptgsHere->dim,ptgsHere->numnzrows);
  }
  if (ctype==USE_NONE) {
    MemCopy((CharPtr)(FillHere),(CharPtr)(f4memoPtr(pf4tg)),bzbufsize);
    if (GetEndian()==ENDIAN_UNKNOWN)
	ErrPostEx(SEV_FATAL,1,24,"Unknown endian system, cannot read in Trajectory Distribution");
    if (GetEndian()!=data_endian)
      FlipBytes((CharPtr)(ptgsHere->TrajGraph),(Int4)((ptgsHere->dim)*(ptgsHere->dim)));
  }
  /* flip bytes only if opposite from stored format from machine of origin */
  ptgsHere->TrajIntegral=f4long(pf4integral);
  ptgsHere->Peak=f4long(pf4peak);
  ptgsHere->pCis=(FloatLo)(f4double(pf4pcis));
  if (aahere=='P') {
    bzbufsize=f4long(pf4cisBZbufsize);  
    FillHere=(CharPtr)(ptgsHere->CisTrajGraph);
    actbufsize=(unsigned int)(ptgsHere->dim)*(ptgsHere->dim)*sizeof(Int4);
    if (ctype==USE_BZ) {
      if ((err=BZ2_bzBuffToBuffDecompress((char *)FillHere,&actbufsize,(char *)f4memoPtr(pf4cistg),(unsigned int)bzbufsize,0,0))!=BZ_OK) { 
        PurgeGlobs();
        ErrPostEx(SEV_FATAL,err,1,"Bzip decompression error occurred");
      }
      if (GetEndian()==ENDIAN_UNKNOWN)
	ErrPostEx(SEV_FATAL,1,24,"Unknown endian system, cannot read in Trajectory Distribution");
      if (GetEndian()!=data_endian)
        FlipBytes((CharPtr)(ptgsHere->CisTrajGraph),(Int4)((ptgsHere->dim)*(ptgsHere->dim)));
    }
    if (ctype==USE_RLE) {
      /* bzbufsize holds num chars in RLE file written */
      if ((err=RLEUnPack((Uint4 *)FillHere,(Uint4 *)f4memoPtr(pf4cistg),bzbufsize/sizeof(Int4),data_endian))!=ptgsHere->dim*ptgsHere->dim)
	ErrPostEx(SEV_FATAL,23,24,"RLEUnPack failed, size=%ld, should=%ld - likely this is caused by overclocked or faulty RAM chips, please test your RAM",(long)err,(long)ptgsHere->dim*ptgsHere->dim);
    }
    if (ctype==USE_NONE) {
      MemCopy((CharPtr)(FillHere),(CharPtr)(f4memoPtr(pf4cistg)),bzbufsize);
      if (GetEndian()==ENDIAN_UNKNOWN)
	ErrPostEx(SEV_FATAL,1,24,"Unknown endian system, cannot read in Trajectory Distribution");
      if (GetEndian()!=data_endian)
        FlipBytes((CharPtr)(ptgsHere->CisTrajGraph),(Int4)((ptgsHere->dim)*(ptgsHere->dim)));
    }
  }
  else
    ptgsHere->CisTrajGraph=NULL;
  ptgsHere->CisTrajIntegral=f4long(pf4cisintegral);
  ptgsHere->CisPeak=f4long(pf4cispeak);
  ptgsHere->pSS=(FloatLo)(f4double(pf4pss));
  ptgsHere->ChiWMean=(FloatLo)(f4double(pf4chiwmean));
  ptgsHere->ChiWSD=(FloatLo)(f4double(pf4chiwsd));
  ptgsHere->AA=(Char)f4char(pf4aaname);
  ptgsHere->nelt0=f4long(pf4nelt0);
  ptgsHere->nelt5=f4long(pf4nelt5);
  ptgsHere->nelt10=f4long(pf4nelt10);
  ptgsHere->nelt15=f4long(pf4nelt15);
  ptgsHere->tout=(Int2)f4long(pf4tout);
  ptgsHere->markovsf=(FloatLo)(f4double(pf4markovsf));
  ptgsHere->rotid=f4long(pf4rot);
  numfrag=f4long(pf4numfrag);
  if (numfrag)
  	pcHere=(CharPtr)(f4memoPtr(pf4frag));
	pfdsLastFrag=NULL;
	/* read in 2-D fragment array - switch columns by checking for prev==NULL */
  while (numfrag>0) {
  	pfdsHere=MemNew(sizeof(FDS));
    MemCopy((VoidPtr)pfdsHere,pcHere,sizeof(FDTS));
    isnewrow=((PFDTS)pcHere)->newrow;
    ((PFDTS)pcHere)->newrow=0;
    if (ptgsHere->pfdsFragmentHead==NULL) {
    	ptgsHere->pfdsFragmentHead=pfdsHere;
			pfdsLastFrag=pfdsHere;
			pfdsHere->prev=NULL;
    }
    else if (!isnewrow) {
/*    pfdsHere->prev!=NULL) {*/
    	pfdsLast->next=pfdsHere;
    	pfdsHere->prev=pfdsLast;
    }
    else { /* isnewrow */
/*pfdsHere->prev==NULL;*/
    	pfdsLast->next=NULL;
    	pfdsLastFrag->nextfrag=pfdsHere;
    	pfdsLastFrag=pfdsHere;
    }
    pfdsHere->nextfrag=NULL;
   	pfdsLast=pfdsHere;
    if (GetEndian()==ENDIAN_UNKNOWN)
			ErrPostEx(SEV_FATAL,1,24,"Unknown endian system, cannot read in Trajectory Distribution");
    if (GetEndian()!=data_endian) {
      if (FlipFDS(pfdsHere)!=ERR_SUCCESS)
				ErrPostEx(SEV_FATAL,1,25,"Cannot read in Trajectory Distribution");
		}
   	pcHere+=sizeof(FDTS);
   	numfrag--;
   	if (!numfrag)
			pfdsHere->next=NULL;
  }
  return ptgsHere;
}

void TrajGraphIntegrate(PTGS ptgsHere)
{
  Int4 cnt,cntmax;
  Int4 *tg,*ctg;
  Int4 itg,citg;

  cntmax=(ptgsHere->dim)*(ptgsHere->dim);
  tg=ptgsHere->TrajGraph;
  ctg=ptgsHere->CisTrajGraph;
  itg=ptgsHere->TrajIntegral;
  tg[cntmax-1]=itg-tg[cntmax-1];
  for (cnt=cntmax-2;cnt>=0;cnt--) {
    tg[cnt]=tg[cnt+1]-tg[cnt];
  }
  /* integrate cis- TG if present as well */
  if (ctg!=NULL) {
    citg=ptgsHere->CisTrajIntegral;
    ctg[cntmax-1]=citg-ctg[cntmax-1];
    for (cnt=cntmax-2;cnt>=0;cnt--) {
      ctg[cnt]=ctg[cnt+1]-ctg[cnt];
    }
  }
}

void TrajGraphDifferentiate(PTGS ptgsHere)
{
  Int4 cnt,cntmax;
  Int4 *tg,*ctg;

  cntmax=(ptgsHere->dim)*(ptgsHere->dim);
  tg=ptgsHere->TrajGraph;
  ctg=ptgsHere->CisTrajGraph;
  for (cnt=0;cnt<cntmax-1;cnt++) {
    tg[cnt]=tg[cnt+1]-tg[cnt];
  }
  tg[cntmax-1]=ptgsHere->TrajIntegral-tg[cntmax-1];
  /* integrate cis- TG if present as well */
  if (ctg!=NULL) {
    for (cnt=0;cnt<cntmax-1;cnt++) {
      ctg[cnt]=ctg[cnt+1]-ctg[cnt];
    }
    ctg[cntmax-1]=ptgsHere->CisTrajIntegral-ctg[cntmax-1];
  }
}

/* scales trajectory graph values by ScaleFactor and recalculates integral exactly */
/* other values like peak, numnzrows, etc. should be recalculated after calling this */
void TrajScale(PTGS ptgsHere,FloatLo ScaleFactor)
{
  Int4 csum=0,cnt;

  for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
    (ptgsHere->TrajGraph)[cnt]=(Int4)(ScaleFactor*(FloatLo)((ptgsHere->TrajGraph)[cnt]));
    csum+=((ptgsHere->TrajGraph)[cnt]);
  }
  ptgsHere->TrajIntegral=csum;
  if ((ptgsHere->CisTrajGraph)!=NULL) {
    csum=0;
    for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
      (ptgsHere->CisTrajGraph)[cnt]=(Int4)(ScaleFactor*(FloatLo)((ptgsHere->CisTrajGraph)[cnt]));
      csum+=((ptgsHere->CisTrajGraph)[cnt]);
    }
  }
  ptgsHere->CisTrajIntegral=csum;
}

/* produce a smoothing/blur filter */
void MakeSmoothFilter(TrajFilter filt)
{
	Int4 cnt;

	for (cnt=0;cnt<(2*TRAJFILTDIM+1)*(2*TRAJFILTDIM+1);cnt++)
		filt[cnt]=1;
}

/* returns a Gaussian filter (for trajectory graphs) with given
magnitude and standard deviation.  An s.d. of about 1 or a bit
higher is typical, while magnitude should typically be anywhere from
10 to 100 */
void MakeGaussFilter(TrajFilter filt,FloatLo magnitude,FloatLo sd)
{
	Int4 row,col,cnt;
	FloatLo dist;

	col=-TRAJFILTDIM;
	row=-TRAJFILTDIM;
	cnt=0;
	while (row<TRAJFILTDIM+1) {
		while (col<TRAJFILTDIM+1) {
			dist=(FloatLo)(-(col*col+row*row));
			filt[cnt]=(Int4)(magnitude*exp(dist/(2.0*sd*sd)));
			cnt++;
			col++;
		}
		col=-TRAJFILTDIM;
		row++;
	}
}

/* returns a sinc (low pass) filter (for trajectory graphs) with given
magnitude.  Here sd is a scale factor so that the first zero of the sinc
function will be at sd (in units of distance from center of filter) */
void MakeLPFilter(TrajFilter filt,FloatLo magnitude,FloatLo sd)
{
	Int4 row,col,cnt;
	FloatLo dist;

	col=-TRAJFILTDIM;
	row=-TRAJFILTDIM;
	cnt=0;
	while (row<TRAJFILTDIM+1) {
		while (col<TRAJFILTDIM+1) {
			dist=(FloatLo)(-(col*col+row*row));
			if (col==0 && row==0)
				filt[cnt]=(Int4)magnitude;
			else
				filt[cnt]=(Int4)(magnitude*sin(dist*PI/sd)/(dist*PI/sd));
			cnt++;
			col++;
		}
		col=-TRAJFILTDIM;
		row++;
	}
}

/* takes trajectory distribution in src (of dimension graphwidth) and
applies digital filter filt, storing the result in dest.  src and dest
should be distinct buffers or the result will be undefined (i.e. chaos) */ 
void DoFilt(Int4 *src,Int4 *dest,Int4 graphwidth,TrajFilter filt)
{
  Int4 c1,c2,filtint,actualdim,ivalue,rowabs,colabs;
  Int4 curpos,lastpos,filtsum,row,col,lastrow,lastcol,filtidx,filtpos;
  FloatLo remain,fvalue;

  filtint=0;
  actualdim=2*TRAJFILTDIM+1;
  for (c1=0;c1<actualdim;c1++)
	  for (c2=0;c2<actualdim;c2++)
		filtint+=filt[c1*actualdim+c2];
  if (filtint<=0) {
	ErrPostEx(SEV_ERROR,23,2,"Filter Integral non-positive: %d",filtint);
	return;
  }
  /* set position to upper left */
  if (WALKTYPE==WALK_PHIPSI)
	curpos=0;
  else
  	curpos=graphwidth*TRAJFILTDIM;
  /* when to stop... */
  if (WALKTYPE==WALK_PHIPSI)
	lastpos=graphwidth*graphwidth;
  else
	lastpos=graphwidth*(graphwidth-TRAJFILTDIM);
  remain=0.0;
  while (curpos<lastpos) {
	filtsum=0;
	row=(curpos/graphwidth-TRAJFILTDIM);
	col=curpos-TRAJFILTDIM;
	lastrow=curpos/graphwidth+TRAJFILTDIM+1;
	lastcol=curpos+TRAJFILTDIM+1;
	filtidx=0;
	/* set filtpos to upper left of filter area */
	while (row<lastrow) {
		/* row % graphwidth takes care of -ve rows, etc */
		rowabs=row%graphwidth;
		if (rowabs<0)
			rowabs+=graphwidth;
		colabs=col%graphwidth;
		if (colabs<0)
			colabs+=graphwidth;
		filtpos=colabs+rowabs*graphwidth;
		filtsum+=(src[filtpos]*filt[filtidx]);
		col++;
		filtidx++;
		if (col==lastcol) {
			col=curpos-TRAJFILTDIM;
			row++;
		}
	}
	/* ensure no negative values on trajectory graph */
	if (filtsum<0)
		filtsum=0;
	fvalue=(FloatLo)filtsum/(FloatLo)filtint;
	ivalue=(Int4)(floor(fvalue));
	remain=remain+fvalue-floor(fvalue);
	while (remain>0.5) {
		remain-=1.0;
		ivalue++;
	}
	dest[curpos]=ivalue;
  	curpos++;
  }
  if (WALKTYPE==WALK_CA) {
	  /* copy top and bottom edges */
	  for (curpos=0;curpos<graphwidth*TRAJFILTDIM;curpos++)
	    dest[curpos]=src[curpos];
	  for (curpos=graphwidth*(graphwidth-TRAJFILTDIM);curpos<graphwidth*graphwidth;curpos++)
	    dest[curpos]=src[curpos];
  }
}

void TrajCopy(Int4 *src,Int4 *dest,Int4 graphwidth)
{
  Int4 *end;

  end=src+graphwidth*graphwidth;
  while (src<end) {
  	*dest=*src;
	src++;
	dest++;
  }
}

/* applies a digital filter to a trajectory graph, wrapping around
correctly in the theta direction; since we are really on a sphere, it
becomes difficult to deal with the poles.  Luckily the "allowed" regions
of trajectory space lie in a relatively narrow band near the equator.
Note that our "filter" will get somewhat distorted as it nears the poles
but otherwise should work fine... */
/* Note the top few rows (near the poles) are NOT adjusted by the filter
but since they are usually zero, this should have little or no effect */
/* applies filter to Cis trajectory distribution if present as well */
void TrajFilt(PTGS ptgsHere,TrajFilter filt)
{
  Int4 *dest,*tmp;
  Int4 cnt,integ;

  dest=(Int4 *)MemNew(sizeof(Int4)*ptgsHere->dim*ptgsHere->dim);
  DoFilt(ptgsHere->TrajGraph,dest,ptgsHere->dim,filt);
  /* in case filter had negative values, recalc integrals */
  tmp=dest;
  integ=0;
  for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
    integ+=*tmp;
    tmp++;
  }  
  ptgsHere->TrajIntegral=integ;
  /* copy back from buffer */
  TrajCopy(dest,ptgsHere->TrajGraph,ptgsHere->dim);
  if ((ptgsHere->CisTrajGraph)!=NULL) {
  	DoFilt(ptgsHere->CisTrajGraph,dest,ptgsHere->dim,filt);
	tmp=dest;
  	integ=0;
  	for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
    		integ+=*tmp;
    		tmp++;
  	}  
 	ptgsHere->CisTrajIntegral=integ;
  	/* copy back from buffer */
  	TrajCopy(dest,ptgsHere->CisTrajGraph,ptgsHere->dim);
  }
  MemFree(dest);
  /* recalculate affected values */
  /* nelts will be for trans now even if should be for cis */
  TrajCalcSparsity(ptgsHere,0, NULL);
  TrajCalcTout(ptgsHere);
}

/* adds Add trajectory graph to Dest trajectory graph */
TrajErr TrajAdd(PTGS ptgsDest,PTGS ptgsAdd)
{
  Int4 csum=0,cnt;

  /* ensure they are the same size */
  if (ptgsDest->dim != ptgsAdd->dim) return ERR_FAIL;
  for (cnt=0;cnt<(ptgsDest->dim)*(ptgsDest->dim);cnt++) {
    (ptgsDest->TrajGraph)[cnt]=(ptgsDest->TrajGraph)[cnt]+(ptgsAdd->TrajGraph)[cnt];
    csum+=((ptgsDest->TrajGraph)[cnt]);
  }
  ptgsDest->TrajIntegral=csum;
  if ((ptgsDest->CisTrajGraph)!=NULL && (ptgsAdd->CisTrajGraph)!=NULL) {
    csum=0;
    for (cnt=0;cnt<(ptgsDest->dim)*(ptgsDest->dim);cnt++) {
      (ptgsDest->CisTrajGraph)[cnt]=(ptgsDest->CisTrajGraph)[cnt]+(ptgsAdd->CisTrajGraph)[cnt];
    csum+=((ptgsDest->CisTrajGraph)[cnt]);
    }
    ptgsDest->CisTrajIntegral=csum;
  }
  return ERR_SUCCESS;
}

void TrajCalcTout(PTGS ptgsHere)
{
	FloatLo fltmp;
	Int2 trajwidth;

	trajwidth=ptgsHere->dim;
        fltmp=6.0*(1.0-(FloatLo)(ptgsHere->nelt0)/((FloatLo)(trajwidth*trajwidth)));
        fltmp+=2.0*(1.0-(FloatLo)(ptgsHere->nelt5)/((FloatLo)(trajwidth*trajwidth)));
        fltmp-=1.0*(1.0-(FloatLo)(ptgsHere->nelt10)/((FloatLo)(trajwidth*trajwidth)));
        fltmp-=3.0*(1.0-(FloatLo)(ptgsHere->nelt15)/((FloatLo)(trajwidth*trajwidth))); 
        fltmp*=4.0;
        ptgsHere->tout=(Int4)(TIMEOUT*fltmp);
        /* ensure arbitrary minimum and maximum tries/residue */
        if ((Int4)(TIMEOUT*fltmp)<BACKTRACK_TRIES_MIN)   
        	ptgsHere->tout=(Int4)BACKTRACK_TRIES_MIN;
        if ((Int4)(TIMEOUT*fltmp)>BACKTRACK_TRIES_MAX)   
        	ptgsHere->tout=(Int4)BACKTRACK_TRIES_MAX;
           

}

/* fills in nelt entries, plus Peak value - uses trans trajectory
   graph unless cis=1 is given */



void TrajCalcSparsity(PTGS ptgsHere, Int2 cis, FILE * ara)
{
  Int4 num0=0,num5=0,num10=0,num15=0,num50=0,pk,pk5,pk10,pk15,pk50,cnt,tmp;

  /* start by finding peak values */
  pk=0;
  for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
    if ((ptgsHere->TrajGraph)[cnt]>pk)
      pk=(ptgsHere->TrajGraph[cnt]);
  }
  ptgsHere->Peak=pk;
  if (ptgsHere->CisTrajGraph!=NULL) {
    pk=0;
    for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
      if ((ptgsHere->CisTrajGraph)[cnt]>pk)
        pk=(ptgsHere->CisTrajGraph[cnt]);
    }
    ptgsHere->CisPeak=pk;
  }
  else
    ptgsHere->CisPeak=0;
  if (cis)
    pk=ptgsHere->CisPeak;
  else
    pk=ptgsHere->Peak;
  pk5=pk/20;
  pk10=pk/10;
  pk15=pk*3/20;
  pk50=pk/2;
  /* avoid crash from improperly detected cis here */	 
  if (cis && ptgsHere->CisTrajGraph==NULL) {	 
    ErrPostEx(SEV_ERROR,7,5,"Warning: TrajCalcSparsity called with cis on non-Cis residue, using trans");	 
    cis=0;	 
  }
  if (cis) {
    for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
      /* only apply to cis- TG if P(cis)>0.5 */
      tmp=(ptgsHere->CisTrajGraph)[cnt];
      /* check if < 15% max value, 10%, 5% and 0 */
      if (tmp>pk50) {
    	  num50++;
	  }
      if (tmp<=pk15) {
        num15++;
      if (tmp<=pk10) {
	    num10++;
	  if (tmp<=pk5)  {
	    num5++;
	  if (tmp<1)
	      num0++;
	  }}}
	}  
  }
  else {
    for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere->dim);cnt++) {
      /* only apply to cis- TG if P(cis)>0.5 */
      tmp=(ptgsHere->TrajGraph)[cnt];
      /* check if < 15% max value, 10%, 5% and 0 */
      if (tmp>pk50) {
    	  num50++;
	  }
      if (tmp<=pk15) {
        num15++;
      if (tmp<=pk10) {
	    num10++;
	  if (tmp<=pk5)  {
	    num5++;
	  if (tmp<1)
	      num0++;
	  }}}
	}
  }
  ptgsHere->nelt0=num0;
  ptgsHere->nelt5=num5;
  ptgsHere->nelt10=num10;
  ptgsHere->nelt15=num15;

  /* Calculate FAHM, write to ARA  */  
  if (ara) fprintf(ara, "%c\t%d\n",ptgsHere->AA ,num50);
  
}

void TrajCalcNZ(PTGS ptgsHere)
{
  Int4 row,col,firstrow=1,lastrow=ptgsHere->dim;
  
  for (row=0;row<ptgsHere->dim;row++) {
    col=0;
    while ((col<ptgsHere->dim) && ((ptgsHere->TrajGraph)[col+row*(ptgsHere->dim)]==0))
      col++;
    if (col<ptgsHere->dim) {
      firstrow=row+1;    
      break;
    }
  }
  /* entire TG is zero */
  if (row==ptgsHere->dim) {
    ptgsHere->numnzrows=0;
    ptgsHere->firstnzrow=(ptgsHere->dim)+1;    
    return;
  }
  for (row=(ptgsHere->dim)-1;row>=0;row--) {
    col=0;
    while ((col<ptgsHere->dim) && ((ptgsHere->TrajGraph)[col+row*(ptgsHere->dim)]==0))
      col++;
    if (col<ptgsHere->dim) {
      lastrow=row+1;    
      break;
    }
  }
  if (lastrow<firstrow) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,30,6,"Unexpected logic error - this should never happen");
  }
  ptgsHere->numnzrows=1+lastrow-firstrow;
  ptgsHere->firstnzrow=firstrow;  
  
  return;
}

/* free NOE list in memory */
PNN FreePNNList(void)
{
  PNN pnnNext,pnnHead;

  pnnHead=pnnDistConstTrueHead;
  while (pnnHead!=NULL) {
    pnnNext=pnnHead->next;
    MemFree(pnnHead);
    pnnHead=pnnNext;
  }
  pnnDistConstTrueHead=NULL;
  pnnHead=pnnDistConstHead;
  while (pnnHead!=NULL) {
    pnnNext=pnnHead->next;
    MemFree(pnnHead);
    pnnHead=pnnNext;
  }
  pnnDistConstHead=NULL;
  return NULL;
}

void BuildBiostrucStub(BiostrucPtr *bspstub)
{
	ValNodePtr pub=NULL,pubequiv=NULL,pubsub=NULL;
	CitArtPtr cap=NULL;
	AuthorPtr ap=NULL;
	CitJourPtr cjp=NULL;
	CitSubPtr csp=NULL;
	ImprintPtr imp=NULL;
	AffilPtr affil=NULL;
	BiostrucHistoryPtr bhp=NULL;
	DbtagPtr dp=NULL;
 
	*bspstub=BiostrucNew();
	(*bspstub)->id=NULL;	
	/* default to mmdb id -9999 to avoid conflicts */
	ValNodeAddInt(&((*bspstub)->id),(Int2)BiostrucId_mmdb_id,(Int4)-9999);	
	(*bspstub)->descr=NULL;
	ValNodeCopyStr(&((*bspstub)->descr),(Int2)BiostrucDescr_name,"RAND");	
	ValNodeCopyStr(&((*bspstub)->descr),(Int2)BiostrucDescr_pdb_comment,"remark  3: Ab Initio Random Protein Structure");
	ValNodeCopyStr(&((*bspstub)->descr),(Int2)BiostrucDescr_pdb_comment,"remark  3: generated by the TRADES algorithm");
	bhp=BiostrucHistoryNew();
	bhp->replaces=NULL;
	bhp->replaced_by=NULL;
	bhp->data_source=BiostrucSourceNew();
	(bhp->data_source)->OBbits__=(Int4)0;
	(bhp->data_source)->name_of_database=StringSave("Ab Initio Structure");
	(bhp->data_source)->database_entry_id=NULL;
	dp=DbtagNew();
	dp->db=StringSave("PDB");
	dp->tag=ObjectIdNew();
	(dp->tag)->str=StringSave("RAND");
	ValNodeAddPointer(&((bhp->data_source)->database_entry_id),(Int2)BiostrucId_other_database,(Pointer)dp);	
	/* today's date */
	(bhp->data_source)->database_entry_date=DateCurr();
	ValNodeAddPointer(&((*bspstub)->descr),(Int2)BiostrucDescr_history,(Pointer)bhp);
	csp=CitSubNew();
	csp->authors=AuthListNew();
	(csp->authors)->choice=(Uint1)1; /* std */
	ap=AuthorNew();
	ap->name=PersonIdNew();
	(ap->name)->choice=(Uint1)2; /* name-std */
	(ap->name)->data=(Pointer)NameStdNew();
	(((NameStdPtr)((ap->name)->data))->names)[0]=StringSave("Feldman");
	(((NameStdPtr)((ap->name)->data))->names)[1]=StringSave("Howard");
	(((NameStdPtr)((ap->name)->data))->names)[2]=StringSave("J.");
	(((NameStdPtr)((ap->name)->data))->names)[3]=StringSave("Howard J. Feldman");  /* me */
	(((NameStdPtr)((ap->name)->data))->names)[4]=StringSave("H.J.");
	(ap->lr)[0]=(Uint1)1;  /* primary author */
	(ap->lr)[1]=(Uint1)1;  /* compiler */
	ap->is_corr=(Uint1)0;  /* not corresponding author */
	ap->affil=NULL;
	ValNodeAddPointer(&((csp->authors)->names),(Int2)0,(Pointer)ap);
	ap=AuthorNew();
	ap->name=PersonIdNew();
	(ap->name)->choice=(Uint1)2; /* name-std */
	(ap->name)->data=(Pointer)NameStdNew();
	(((NameStdPtr)((ap->name)->data))->names)[0]=StringSave("Hogue");
	(((NameStdPtr)((ap->name)->data))->names)[1]=StringSave("Christopher");
	(((NameStdPtr)((ap->name)->data))->names)[2]=StringSave("W.V.");
	(((NameStdPtr)((ap->name)->data))->names)[3]=StringSave("Christopher W.V. Hogue");
	(((NameStdPtr)((ap->name)->data))->names)[4]=StringSave("C.W.V.");
	(ap->lr)[0]=(Uint1)2;  /* secondary author */
	(ap->lr)[1]=(Uint1)1;  /* compiler */
	ap->is_corr=(Uint1)1;  /* corresponding author */
	ap->affil=NULL;
	ValNodeAddPointer(&((csp->authors)->names),(Int2)0,(Pointer)ap);
	affil=AffilNew();
	/* indicate parsed string */
	affil->choice=2;
	affil->affil=StringSave("National University of Singapore");
	affil->div=StringSave("Department of Biological Sciences/Mechanobiology Institute of Singapore");
	affil->city=StringSave("Singapore");
	affil->sub=StringSave("Singapore");
	affil->country=StringSave("Singapore");
	affil->street=StringSave("5A Engineering Drive 1");
	affil->email=StringSave("chogue@blueprint.org");
	affil->fax=StringSave("+65 6774 2540");
	affil->phone=StringSave("+65 9487 5537");
	affil->postal_code=StringSave("117411");
	(csp->authors)->affil=affil;
	csp->date=DateCurr();
	ValNodeAddPointer(&pubsub,(Int2)PUB_Sub,(Pointer)csp);
	ValNodeAddPointer(&((*bspstub)->descr),(Int2)BiostrucDescr_attribution,(Pointer)pubsub);
	/* build publication reference */
	cap=CitArtNew();
	ValNodeCopyStr(&(cap->title),(Int2)Cit_title_name,"A Fast Method to Sample Real Protein Conformational Space");
	cap->authors=(AuthListPtr)AsnIoMemCopy((Pointer)csp->authors,(AsnReadFunc)AuthListAsnRead,(AsnWriteFunc)AuthListAsnWrite);
	cap->from=(Uint1)1; /* journal */
	cjp=CitJourNew();
	ValNodeCopyStr(&(cjp->title),(Int2)Cit_title_name,"Proteins: Structure, Function, and Genetics");
	ValNodeCopyStr(&(cjp->title),(Int2)Cit_title_iso_jta,"Proteins");
	ValNodeCopyStr(&(cjp->title),(Int2)Cit_title_ml_jta,"Proteins");
	ValNodeCopyStr(&(cjp->title),(Int2)Cit_title_issn,"0887-3585");
	ValNodeCopyStr(&(cjp->title),(Int2)Cit_title_jta,"PTS");
	imp=ImprintNew();
	imp->date=DateNew();
	((imp->date)->data)[0]=(Uint1)1;
	((imp->date)->data)[1]=(Uint1)100;
	((imp->date)->data)[2]=(Uint1)5;
	((imp->date)->data)[3]=(Uint1)1;
	((imp->date)->data)[4]=(Uint1)NOT_SET;
	((imp->date)->data)[5]=(Uint1)NOT_SET;
	((imp->date)->data)[6]=(Uint1)NOT_SET;
	(imp->date)->str=NULL;
	imp->volume=StringSave("39");
	imp->issue=StringSave("2");
	imp->pages=StringSave("112-131");
	imp->language=StringSave("ENG"); 
	cjp->imp=imp;
	cap->fromptr=(Pointer)cjp;
	cap->ids=NULL;	
/*	ValNodeCopyStr(&(cap->ids),(Int2)ARTICLEID_PMPID,"10.1002/(SICI)1097-0134(20000501)39:2<112::AID-PROT2>3.0.CO;2-B");
	ValNodeAddInt(&(cap->ids),(Int2)ARTICLEID_PUBMED,(Int4)10737933);
*/ /* can't write out .prt if these are uncommented, don't know why... */	
	ValNodeAddPointer(&pubequiv,(Int2)PUB_Article,(Pointer)cap);
	ValNodeAddInt(&pubequiv,(Int2)PUB_PMid,(Int4)10737933);
	ValNodeAddPointer(&pub,(Int2)PUB_Equiv,(Pointer)pubequiv);
	ValNodeAddPointer(&((*bspstub)->descr),(Int2)BiostrucDescr_attribution,(Pointer)pub);
	((*bspstub)->descr)->choice=BiostrucDescr_name;
	(*bspstub)->chemical_graph=BiostrucGraphNew();
	((*bspstub)->chemical_graph)->molecule_graphs=MoleculeGraphNew();
	(((*bspstub)->chemical_graph)->molecule_graphs)->residue_sequence=ResidueNew();
	((((*bspstub)->chemical_graph)->molecule_graphs)->residue_sequence)->residue_graph=ValNodeNew(NULL);
	((((*bspstub)->chemical_graph)->molecule_graphs)->residue_sequence)->residue_graph->choice=ResidueGraphPntr_local;
}

/* fills in fields of ASN.1 trajectory structure, including the .HA file
   as a raw octet stream, and writes out the ASN.1 file */
/* frees up pnnDistConstHead as well, removing stored NOE constraints
   from memory */
/* sequence must be in EXTAA_ENCODE format */
void PackAsnTrajGraph(CharPtr fnam,CharPtr sequence,ByteStorePtr id,ByteStorePtr descr,ValNodePtr bioseq,Boolean bzipit)
{
  TGraphPtr tgp=NULL;
  TGraphHeaderPtr tghp=NULL;	
  TGConstraintDataPtr cdpHere=NULL,cdpLast=NULL,cdpHead=NULL;
  PNN pnnHere;
  FILE *f,*g;
  AsnIoPtr paio;
  ByteStorePtr bsfile1,bsfile2,bsfile3;
  CharPtr cfile;
  Int4 filesize,cnt;
  Char tmpfnam[PATH_MAX];
  AsnIoBSPtr aibp=NULL;
  ValNodePtr vnpLast=NULL,vnpTemp,vnpTest;
  Boolean updHist=FALSE;
  BiostrucHistoryPtr bhp=NULL;
  Int4 updID=0;
  BiostrucReplacePtr brp=NULL;
  Int4 seqlength;

  TGPack();
  /* encode spoecial amino acids */
  seqlength=0;
  cnt=0;
  while (sequence[cnt]) {
    if (isAA(sequence[cnt]))
      seqlength++;
    cnt++;
  }
  tgp=TGraphNew();
  tghp=TGraphHeaderNew();
  tgp->tgheader=tghp;
  tghp->seq=sequence;
  tghp->bseq=bioseq;
  tghp->seqlen=seqlength;
  tghp->bberrtol=BACKBONE_ERROR_TOLERANCE;
  tghp->bbacc=BACKBONE_PRECISION;
  tghp->numrot=NUM_ROT_TRIES;
  tghp->atmbncbb=ATOM_BOUNCINESS_BB;
  tghp->atmbncsc=ATOM_BOUNCINESS_SC;
  tghp->bumph=BUMPCHECK_HYDROGEN;
  tghp->incsize=INCSIZE;
  tghp->tunnel=TUNNEL_PROB;
  tghp->startbb=START_BACKBONE;
  tghp->walktype=WALKTYPE;
  tghp->tgunits=TGUNITS;
  tghp->constrfile=CONSTRAINT_FILE;
  tghp->trajtype=TRAJTYPE;
  tghp->trajdiv=TRAJDIV;
  tghp->timeout=TIMEOUT;
  /* fill in biostruc stub */
  BuildBiostrucStub(&(tghp->bsstub));
  if (id!=NULL) {
    /* memory leak patched? */
    AsnGenericChoiceSeqOfFree((tghp->bsstub)->id,(AsnOptFreeFunc)BiostrucIdFree);
    /*BiostrucIdFree((tghp->bsstub)->id);*/
    (tghp->bsstub)->id=NULL;
    BSSeek(id,0L,SEEK_SET);
    aibp=AsnIoBSOpen("rb",id);
    while (BSTell(id)!=BSLen(id)) {
	vnpTemp=BiostrucIdAsnRead(aibp->aip,NULL);
	vnpTemp->next=NULL;
	if ((tghp->bsstub)->id==NULL) {
		/* first node */
		(tghp->bsstub)->id=vnpTemp;
	}
	else
		vnpLast->next=vnpTemp;
	/* negate MMDB id to avoid confusing with real structure */
	if (vnpTemp->choice==BiostrucId_mmdb_id && (vnpTemp->data.intvalue>0)) {
		/* must be starting from a real structure so we will update Biostruc history too */
		updHist=TRUE;
		updID=vnpTemp->data.intvalue;
		vnpTemp->data.intvalue=-(vnpTemp->data.intvalue);
	}
	vnpLast=vnpTemp;
    }
    aibp=AsnIoBSClose(aibp);
  }
  if (descr!=NULL) {
    /* memory leak patched? */
    AsnGenericChoiceSeqOfFree((tghp->bsstub)->descr,(AsnOptFreeFunc)BiostrucDescrFree);
/*    BiostrucDescrFree((tghp->bsstub)->descr);*/
    (tghp->bsstub)->descr=NULL;
    BSSeek(descr,0L,SEEK_SET);
    aibp=AsnIoBSOpen("rb",descr);
    while (BSTell(descr)!=BSLen(descr)) {
	    vnpTemp=BiostrucDescrAsnRead(aibp->aip,NULL);
	    vnpTemp->next=NULL;
	    if ((tghp->bsstub)->descr==NULL) {
		/* first node */
		(tghp->bsstub)->descr=vnpTemp;
	    }
	    else
		vnpLast->next=vnpTemp;
	    /* force name to RAND */
	    if (vnpTemp->choice==BiostrucDescr_name) {
		MemFree(vnpTemp->data.ptrvalue);
		vnpTemp->data.ptrvalue=StringSave("RAND");
	    }
	    else if (vnpTemp->choice==BiostrucDescr_history) {
		bhp=(BiostrucHistoryPtr)(vnpTemp->data.ptrvalue);
		if (updHist==TRUE) {
			/* update replaces field */
			brp=BiostrucReplaceNew();
			brp->id=NULL;
			/* give old MMDB id */
			ValNodeAddInt(&(brp->id),(Int2)BiostrucId_mmdb_id,(Int4)updID);	
			brp->date=DateCurr();
			if (bhp->replaces!=NULL)
				BiostrucReplaceFree(bhp->replaces);
			bhp->replaces=brp;
		}
		vnpTest=(bhp->data_source)->database_entry_id;
		while (vnpTest!=NULL) {
			if (vnpTest->choice==BiostrucId_other_database) {
				MemFree((((DbtagPtr)(vnpTest->data.ptrvalue))->tag)->str);
				(((DbtagPtr)(vnpTest->data.ptrvalue))->tag)->str=StringSave("RAND");
			}
			vnpTest=vnpTest->next;
		}
	    }
	    vnpLast=vnpTemp;
    }
    aibp=AsnIoBSClose(aibp);
  }
  /* parse distance constraints now */
  pnnHere=pnnDistConstTrueHead;
  while (pnnHere) {
    cdpHere=TGConstraintDataNew();
    cdpHere->res1=(Int4)(pnnHere->res1);
    cdpHere->res2=(Int4)(pnnHere->res2);
    cdpHere->atomname1=StringSave(pnnHere->AtomName1);
    cdpHere->atomname2=StringSave(pnnHere->AtomName2);
    cdpHere->meandist=(FloatHi)(pnnHere->MeanDist);
    cdpHere->mindelta=(FloatHi)(pnnHere->MinDelta);
    cdpHere->maxdelta=(FloatHi)(pnnHere->MaxDelta);
    cdpHere->angle1=(FloatHi)(pnnHere->Angle1);
    cdpHere->angle2=(FloatHi)(pnnHere->Angle2);
    cdpHere->dihedral01=(FloatHi)(pnnHere->Dihedral01);
    cdpHere->dihedral12=(FloatHi)(pnnHere->Dihedral12);
    cdpHere->dihedral23=(FloatHi)(pnnHere->Dihedral23);
    cdpHere->prob=(FloatHi)(pnnHere->prob);
    if (cdpLast!=NULL)
      cdpLast->next=cdpHere;
    else
      cdpHead=cdpHere;
    cdpLast=cdpHere;
    pnnHere=pnnHere->next;
  }
  if (cdpHere!=NULL)
    cdpHere->next=NULL;
  tghp->distconstr=cdpHead;
  StringCpy(tmpfnam,tmpdbasename);
  StringCat(tmpfnam,DB_EXT1);
  filesize=FileLength(tmpfnam);
  if ((f=FileOpen(tmpfnam,"rb"))==NULL) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,99,99,"Cannot open database file %s",tmpfnam);
  }
  cfile=(CharPtr)MemNew(sizeof(Char)*filesize);
  FileRead(cfile,sizeof(Char),filesize,f);
  FileClose(f);
  bsfile1=BSNew(sizeof(Char)*filesize);
  BSSeek(bsfile1,0L,SEEK_SET);
  BSWrite(bsfile1,cfile,filesize);
  cfile=MemFree(cfile);
  tgp->datacdx=bsfile1;
  StringCpy(tmpfnam,tmpdbasename);
  StringCat(tmpfnam,DB_EXT2);
  filesize=FileLength(tmpfnam);
  if ((f=FileOpen(tmpfnam,"rb"))==NULL) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,99,99,"Cannot open database file %s",tmpfnam);
  }
  cfile=(CharPtr)MemNew(sizeof(Char)*filesize);
  FileRead(cfile,sizeof(Char),filesize,f);
  FileClose(f);
  bsfile2=BSNew(sizeof(Char)*filesize);
  BSSeek(bsfile2,0L,SEEK_SET);
  BSWrite(bsfile2,cfile,filesize);
  cfile=MemFree(cfile);
  tgp->datadbf=bsfile2;
  StringCpy(tmpfnam,tmpdbasename);
  StringCat(tmpfnam,DB_EXT3);
  filesize=FileLength(tmpfnam);
  if ((f=FileOpen(tmpfnam,"rb"))==NULL) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,99,99,"Cannot open database file %s",tmpfnam);
  }
  cfile=(CharPtr)MemNew(sizeof(Char)*filesize);
  FileRead(cfile,sizeof(Char),filesize,f);
  FileClose(f);
  bsfile3=BSNew(sizeof(Char)*filesize);
  BSSeek(bsfile3,0L,SEEK_SET);
  BSWrite(bsfile3,cfile,filesize);
  cfile=MemFree(cfile);
  tgp->datafpt=bsfile3;
  StringCpy(tmpfnam,fnam);
  StringCat(tmpfnam,ASN_EXT);
  paio=AsnIoOpen(tmpfnam,"wb");
  if (TGraphAsnWrite(tgp,paio,NULL)==FALSE) {
    tgp=TGraphFree(tgp);
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,97,97,"Unable to write out ASN.1");
  }
  AsnIoClose(paio);
  /* don't want free-er to free the string constants, allocated elsewhere! */
  tghp->seq=NULL;
  tghp->constrfile=NULL;
  tghp->bseq=NULL; /* cwvh */
  /* only free these if user says so, may have been allocated elsewhere */
  tgp=TGraphFree(tgp);
  /* free NOE list from memory now too */
  pnnDistConstTrueHead=FreePNNList();
  /* bzip compress file if requested */
  if (bzipit==TRUE) {
    StringCpy(tmpfnam,fnam);
    StringCat(tmpfnam,ASN_EXT);
    if ((f=FileOpen(tmpfnam,"rb"))==NULL) {
      PurgeGlobs();
      ErrPostEx(SEV_FATAL,99,99,"Cannot open database file %s",tmpfnam);
    }
    StringCpy(tmpfnam,fnam);
    StringCat(tmpfnam,ASN_EXT);
    StringCat(tmpfnam,BZ_EXT);
    if ((g=FileOpen(tmpfnam,"wb"))==NULL) {
      PurgeGlobs();
      ErrPostEx(SEV_FATAL,99,99,"Cannot open output file %s",tmpfnam);
    }
    compressStream(f,g);
    StringCpy(tmpfnam,fnam);
    StringCat(tmpfnam,ASN_EXT);
    FileRemove(tmpfnam);
  }
}

/* robust tmpnam function - checks for existing files with same name, allows for use
   of DFPTEMP to specify temp file locations, ensures no '.' in filename (nodots=TRUE)
   and check if .dbf or .trj file exists of same name (adddbext=TRUE) */
CharPtr DFPTmpNam(Boolean nodots,Boolean adddbext)
{
	static Char tmpbuf[PATH_MAX];
	Char ftmp[PATH_MAX];
	CharPtr tmpdir,pc,pc2;
	Int4 fs1,fs2;

	do {
  		TmpNam(tmpbuf);
		tmpdir=getenv("DFPTEMP");
		/* for Distributed Folding Project, allow temp files to go in user-specified dir */
		if (tmpdir!=NULL) {
	  		pc=tmpbuf+StringLen(tmpbuf)-1;
  			while (pc>=tmpbuf) {
				if (*pc==DIRDELIMCHR) {
					pc++;
  					break;
				}
  				pc--;
  			}
			if (pc<tmpbuf)
				pc=tmpbuf;
			pc2=StringSave(pc);
			if (tmpdir[StringLen(tmpdir)-1]==DIRDELIMCHR)
				sprintf(tmpbuf,"%s%s",tmpdir,pc2);
			else
				sprintf(tmpbuf,"%s%c%s",tmpdir,DIRDELIMCHR,pc2);
			pc2=MemFree(pc2);
		}
		if (nodots) {
  			/* ensure no dots in the name or Codebase cannot create database - in path part is OK though */
  			pc=tmpbuf+StringLen(tmpbuf)-1;
  			while (pc>=tmpbuf) {
  				if (*pc=='.')
  					*pc='_';
  				if (*pc==DIRDELIMCHR)
  					break;
  				pc--;
  			}
		}
		/* add more randomness for stupid Microsoft TmpNams */
		sprintf(ftmp,"_%04d_%03d",(int)GetAppProcessID(),(int)fabs(999.0*Rand1()));
		StringCat(tmpbuf,ftmp);
		/* filename is complete */
		if (adddbext) {
		    StringCpy(ftmp,tmpbuf);
		    StringCat(ftmp,ASN_EXT);
			fs1=FileLength(ftmp);
			StringCpy(ftmp,tmpbuf);
			StringCat(ftmp,DB_EXT2);
			fs2=FileLength(ftmp);
		}
		else {
			fs1=FileLength(tmpbuf);
			fs2=fs1;
		}
	} while ((fs1!=0) || (fs2!=0));
	return tmpbuf;
}

/* extracts archive file and other fields from an ASN.1 trajectory file */
/* builds NOE storing data structure in memory off pnnDistConstHead */
/* seq will be in EXTAA_ENCODE format */
CharPtr UnPackAsnTrajGraph(CharPtr ASNfnam,Int2 PNTR seqlength,CharPtr seq,CharPtr DumpHeaderName,BiostrucIdPtr *ppbi,BiostrucDescrPtr *ppbd,ValNodePtr *bioseq)
{
  TGraphPtr tgp=NULL;
  TGraphHeaderPtr tghp=NULL;
  TGConstraintDataPtr cdp=NULL;
  PNN pnnHere=NULL,pnnLast=NULL;
  FILE *f,*g;
  AsnIoPtr paio;
  CharPtr cfile;
  Int4 filesize;
  Char tmpfnam[PATH_MAX];

  /* ensure unique name is chosen since we are adding extensions
     to the base name */
  StringCpy(tmpdbasename,DFPTmpNam(TRUE,TRUE));
  /* first check if a bzip compressed file exists, and use it if it 
     does, otherwise assume a non-compressed file should be used */
  StringCpy(tmpfnam,ASNfnam);
  StringCat(tmpfnam,ASN_EXT);
  StringCat(tmpfnam,BZ_EXT);
  if ((f=FileOpen(tmpfnam,"rb"))!=NULL) {
    /* compressed file found so decompress it */
    StringCpy(tmpfnam,ASNfnam);
    StringCat(tmpfnam,ASN_EXT);
    if ((g=FileOpen(tmpfnam,"wb"))==NULL) {
	  FileClose(f);
      PurgeGlobs();
      ErrPostEx(SEV_FATAL,1,1,"Unable to uncompress ASN.1 file %s",tmpfnam);
    }
    uncompressStream(f,g);
    StringCpy(tmpfnam,ASNfnam);
    StringCat(tmpfnam,ASN_EXT);
    StringCat(tmpfnam,BZ_EXT);
    FileRemove(tmpfnam);
  }
  StringCpy(tmpfnam,ASNfnam);
  StringCat(tmpfnam,ASN_EXT);
  if ((paio=AsnIoOpen(tmpfnam,"rb"))==NULL) {
    ErrPostEx(SEV_ERROR,1,1,"Unable to open trajectory distribution file %s",tmpfnam);
    return NULL;
  }
  if ((tgp=TGraphAsnRead(paio,NULL))==NULL) {
	AsnIoClose(paio);
    ErrPostEx(SEV_ERROR,1,1,"Unable to read trajectory distribution file %s",tmpfnam);
    return NULL;
  }
  AsnIoClose(paio);
  tghp=tgp->tgheader;
  if (DumpHeaderName!=NULL) {
          sprintf(tmpfnam,"%s",DumpHeaderName);
 	  /* dump ASN header on unpacking */
	  paio=AsnIoOpen(tmpfnam,"wb");
	  if (TGraphHeaderAsnWrite(tghp,paio,NULL)==FALSE) {
		ErrPostEx(SEV_ERROR,2,3,"Unable to write ASN.1 header");
	  }
	  AsnIoClose(paio);
  }
  StringCpy(seq,tghp->seq);
  *seqlength=(Int2)(tghp->seqlen);
  BACKBONE_ERROR_TOLERANCE=tghp->bberrtol;
  BACKBONE_PRECISION=tghp->bbacc;
  NUM_ROT_TRIES=tghp->numrot;
  ATOM_BOUNCINESS_BB=tghp->atmbncbb;
  ATOM_BOUNCINESS_SC=tghp->atmbncsc;
  BUMPCHECK_HYDROGEN=tghp->bumph;
  INCSIZE=tghp->incsize;
  START_BACKBONE=tghp->startbb;
  WALKTYPE=tghp->walktype;
  TGUNITS=tghp->tgunits;
  StringCpy(CONSTRAINT_FILE,tghp->constrfile);
  TRAJTYPE=tghp->trajtype;
  TRAJDIV=tghp->trajdiv;
  TIMEOUT=tghp->timeout;
  TUNNEL_PROB=tghp->tunnel;
  if (bioseq!=NULL) {
    *bioseq=tghp->bseq;
    tghp->bseq=NULL;
  }
  *ppbi=(tghp->bsstub)->id;
  /* so won't get freed */
  (tghp->bsstub)->id=NULL;
  *ppbd=(tghp->bsstub)->descr;
  /* so won't get freed */
  (tghp->bsstub)->descr=NULL;
  cdp=tghp->distconstr;
  /* free any potential junk in the list before reading in */
  pnnDistConstTrueHead=FreePNNList();
  /* parse distance constraints now */
  while (cdp) {
    pnnHere=(PNN)MemNew(sizeof(NN));
    pnnHere->res1=(Int2)(cdp->res1);
    pnnHere->res2=(Int2)(cdp->res2);
    StringNCpy(pnnHere->AtomName1,cdp->atomname1,5);
    (pnnHere->AtomName1)[4]=0;
    StringNCpy(pnnHere->AtomName2,cdp->atomname2,5);
    (pnnHere->AtomName2)[4]=0;
    pnnHere->MeanDist=(FloatLo)(cdp->meandist);
    pnnHere->MinDelta=(FloatLo)(cdp->mindelta);
    pnnHere->MaxDelta=(FloatLo)(cdp->maxdelta);
    pnnHere->Angle1=(FloatLo)(cdp->angle1);
    pnnHere->Angle2=(FloatLo)(cdp->angle2);
    pnnHere->Dihedral01=(FloatLo)(cdp->dihedral01);
    pnnHere->Dihedral12=(FloatLo)(cdp->dihedral12);
    pnnHere->Dihedral23=(FloatLo)(cdp->dihedral23);
    pnnHere->prob=(FloatLo)(cdp->prob);
    pnnHere->tries=0;
    pnnHere->next=NULL;
    if (pnnLast!=NULL) {
      pnnLast->next=pnnHere;
      pnnHere->prev=pnnLast;
    }
    else
      pnnHere->prev=NULL;
    /* true for first node */
    if (pnnDistConstTrueHead==NULL)
      pnnDistConstTrueHead=pnnHere;
    pnnLast=pnnHere;
    cdp=cdp->next;
  }
  filesize=BSLen((ByteStorePtr)(tgp->datacdx));
  BSSeek((ByteStorePtr)(tgp->datacdx),0L,SEEK_SET);
  cfile=(CharPtr)MemNew(sizeof(Char)*filesize);
  BSRead((ByteStorePtr)(tgp->datacdx),cfile,filesize);
  /* store extracted ASN.1 files in /tmp directory */
  StringCpy(tmpfnam,tmpdbasename);
  StringCat(tmpfnam,DB_EXT1);
  if ((f=FileOpen(tmpfnam,"wb"))==NULL) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,99,99,"Cannot write database file %s to disk",tmpfnam);
  }
  if (FileWrite(cfile,sizeof(Char),filesize,f)!=filesize) {
	FileClose(f);
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,97,99,"Cannot write database file %s to disk, disk full?",tmpfnam);
  }
  FileClose(f);
  cfile=MemFree(cfile);
  filesize=BSLen((ByteStorePtr)(tgp->datadbf));
  BSSeek((ByteStorePtr)(tgp->datadbf),0L,SEEK_SET);
  cfile=(CharPtr)MemNew(sizeof(Char)*filesize);
  BSRead((ByteStorePtr)(tgp->datadbf),cfile,filesize);
  StringCpy(tmpfnam,tmpdbasename);
  StringCat(tmpfnam,DB_EXT2);
  if ((f=FileOpen(tmpfnam,"wb"))==NULL) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,99,99,"Cannot write database file %s to disk",tmpfnam);
  }
  if (FileWrite(cfile,sizeof(Char),filesize,f)!=filesize) {
	FileClose(f);
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,97,99,"Cannot write database file %s to disk, disk full?",tmpfnam);
  }
  FileClose(f);
  cfile=MemFree(cfile);
  filesize=BSLen((ByteStorePtr)(tgp->datafpt));
  BSSeek((ByteStorePtr)(tgp->datafpt),0L,SEEK_SET);
  cfile=(CharPtr)MemNew(sizeof(Char)*filesize);
  BSRead((ByteStorePtr)(tgp->datafpt),cfile,filesize);
  StringCpy(tmpfnam,tmpdbasename);
  StringCat(tmpfnam,DB_EXT3);
  if ((f=FileOpen(tmpfnam,"wb"))==NULL) {
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,99,99,"Cannot write database file %s to disk",tmpfnam);
  }
  if (FileWrite(cfile,sizeof(Char),filesize,f)!=filesize) {
	FileClose(f);
    PurgeGlobs();
    ErrPostEx(SEV_FATAL,97,99,"Cannot write database file %s to disk, disk full?",tmpfnam);
  }
  FileClose(f);
  cfile=MemFree(cfile);
  tgp=TGraphFree(tgp);
  /* global variable */
  return tmpdbasename;
}

/* Adjusts the residue numbers in the temporary database according to an
	 insertion or deletion of a residue record */
TrajErr AdjustResidueNumbers(int numAA, int res_num, int Adjust_Type)
{
	Int4 j;
	
	if (Adjust_Type==EXTAASEQ_INSERT) {
    		for (j=numAA; j>=res_num; j--) {
	    		d4seekDouble(pd4Traj,(double)j);
	    		f4assignInt(pf4res,(j+1));
	 	}
	  return ERR_SUCCESS;
	}	
	else if (Adjust_Type==EXTAASEQ_DELETE) {
		d4seekDouble(pd4Traj,(double) res_num);
    		d4delete(pd4Traj);
    		d4pack(pd4Traj);
    		for (j=res_num+1; j<=numAA; j++) {
	    		d4seekDouble(pd4Traj,(double)j);
	    		f4assignInt(pf4res,(j-1));
	  	}
	  	return ERR_SUCCESS;
	}	
	return ERR_FAIL;	
}

/*void BZCompressFile(FILE *fin, FILE *fout)
{
	BZFILE* bzf=NULL;
	Uchar ibuf[5000];
	Int4 nIbuf;
	unsigned int nbytes_in,nbytes_out;
	int bzerr,bzerr_dummy,ret;
        
        bzf=bzWriteOpen(&bzerr,fout,9,0,50);
        if (bzerr!=BZ_OK) { 
                bzWriteClose(&bzerr_dummy,bzf,1,&nbytes_in,&nbytes_out );
                goto errh;;
        }
        while (1) {   
                if (feof(fin)) break;
                nIbuf=FileRead(ibuf,sizeof(Uchar),5000,fin);
                if (ferror(fin)) goto errh;
                if (nIbuf>0) bzWrite(&bzerr,bzf,(void*)ibuf,nIbuf);
                if (bzerr!=BZ_OK) {
			bzWriteClose(&bzerr_dummy,bzf,1,&nbytes_in,&nbytes_out);
			goto errh;
                }     
        }
        bzWriteClose(&bzerr,bzf,0,&nbytes_in,&nbytes_out);
        if (bzerr!=BZ_OK) goto errh;
        ret=fflush(fout);
        if (ret==EOF) goto errh;
        FileClose(fout);   
        FileClose(fin);
        return;
errh:
        PurgeGlobs();
        ErrPostEx(SEV_FATAL,1,1,"Problem occurred while writing compressed file");
 	return;   
}*/

/*void BZDecompressFile(FILE *fin, FILE *fout)
{
	BZFILE* bzf=NULL;
	int bzerr,bzerr_dummy,ret,nread,i;
	Uchar obuf[5000];
	Uchar unused[5000];
   	int nUnused;
   	Uchar* unusedTmp;
        
	nUnused = 0;
	while (1) {
		bzf=bzReadOpen(&bzerr,fin,0,0,unused,nUnused);
	      	if (bzf==NULL || bzerr!=BZ_OK) {
			bzReadClose(&bzerr_dummy,bzf);
			goto errh2;
		}
	      	while (bzerr==BZ_OK) {
	        	nread=bzRead(&bzerr,bzf,obuf,5000);
	         	if (bzerr==BZ_DATA_ERROR_MAGIC) {
				bzReadClose(&bzerr_dummy,bzf);
				goto errh2;
			}
			if ((bzerr==BZ_OK || bzerr==BZ_STREAM_END) && nread>0)
	            		FileWrite(obuf,sizeof(Uchar),nread,fout);
	         	if (ferror(fout)) goto errh2;
	      	}
	      	if (bzerr!=BZ_STREAM_END) {
		   	bzReadClose(&bzerr_dummy,bzf);
			goto errh2;
		}
	      	bzReadGetUnused(&bzerr,bzf,(void**)(&unusedTmp),&nUnused);
	      	for (i=0;i<nUnused;i++)
			unused[i]=unusedTmp[i];
	        bzReadClose(&bzerr,bzf);
	      	if (nUnused==0 && feof(fin)) break; 
	}
	if (ferror(fin)) goto errh2;
	FileClose(fin);
	if (ferror(fout)) goto errh2;
	ret=fflush(fout);
	if (ret!=0) goto errh2;
	FileClose(fout);
	return;   
errh2:
        PurgeGlobs();
        ErrPostEx(SEV_FATAL,1,1,"Problem occurred while writing decompressed file");
 	return;   
}*/

Int4 FillTG(Int2 TrajType,CharPtr cres,Char SStru,FloatLo prob[3],Int2 trajwidth,Int4 PNTR TrajGraph,Boolean UseCis)
{
	Int4 c1,c2,cntr,res;
	Int4 dictsize;
	FloatLo csum,ctemp=0.0,remain;
	FloatLo ftmp1=1.0,ftmp2=1.0,ftmp3=1.0;
	PTGS ptgsHereH=NULL;
	PTGS ptgsHereE=NULL;
	PTGS ptgsHereC=NULL;
	/* store trajectory graph */
	Int4 TrajIntegral;

	/* go through protein 1 residue at a time */
        res=StringCSpn(aalist,cres);
	if (TrajType!=TRAJ_UNIFORM) {
		csum=0.0;
		remain=0.0;
		TGInit(TRAJ3FNAME,DB_READ,&dictsize);
		if (dictsize!=60) {
			PurgeGlobs();
   			ErrPostEx(SEV_FATAL,1,1,"Invalid trajectory distribution dictionary file: %s",TRAJ3FNAME);
		}
        	ptgsHereH=TrajGraphRead(res*3+1);
        	ptgsHereE=TrajGraphRead(res*3+2);
		ptgsHereC=TrajGraphRead(res*3+3);
		if ((TrajType!=TRAJ_STANDARD) && (UseCis==FALSE)) {
			/* normalize area to size of traj graphs, TRAJDIV*TRAJDIV */
			ftmp1=prob[0]*(FloatLo)trajwidth*(FloatLo)trajwidth/(FloatLo)(ptgsHereH->TrajIntegral);
			ftmp2=prob[1]*(FloatLo)trajwidth*(FloatLo)trajwidth/(FloatLo)(ptgsHereE->TrajIntegral);
			ftmp3=prob[2]*(FloatLo)trajwidth*(FloatLo)trajwidth/(FloatLo)(ptgsHereC->TrajIntegral);
		}
		cntr=0;
		for (c1=0;c1<trajwidth;c1++)	
			for (c2=0;c2<trajwidth;c2++) {
				if (UseCis==TRUE) {
					/* just add 3 ss traj graphs here */
					ctemp=(FloatLo)((ptgsHereH->CisTrajGraph)[cntr])+(FloatLo)((ptgsHereE->CisTrajGraph)[cntr])+(FloatLo)((ptgsHereC->CisTrajGraph)[cntr]);
				}
				else {
					if (TrajType==TRAJ_GOR) {
						/* convolve 3 ss traj graphs */
						/* and normalize E and C trajectory
						   graphs to have same integral as H one
						   before doing so */
						ctemp=(FloatLo)((ptgsHereH->TrajGraph)[cntr]*ftmp1)+(FloatLo)((ptgsHereE->TrajGraph)[cntr]*ftmp2)+(FloatLo)((ptgsHereC->TrajGraph)[cntr]*ftmp3);
						/*
						printf("residue %c%d Htraj: %ld Etraj: %ld Ctraj: %ld pH: %f pE: %f pC: %f intH: %ld intE: %ld intC: %ld ctemp: %ld\n",aalist[res],cnt,TrajGraphH[TRAJ_EL(res,c1,c2)],TrajGraphE[TRAJ_EL(res,c1,c2)],TrajGraphC[TRAJ_EL(res,c1,c2)],
						 arwProbMatrix[cnt][0],arwProbMatrix[cnt][1],arwProbMatrix[cnt][2],aaintegrals[res][0],aaintegrals[res][1],aaintegrals[res][2],ctemp);
						*/
					}
					else if (TrajType==TRAJ_SSTRU) {
						/* just pick predicted SS type */
						switch (SStru) {
							case 'H':
								ctemp=(FloatLo)((ptgsHereH->TrajGraph)[cntr])*ftmp1/prob[0];
								break;
							case 'E':
								ctemp=(FloatLo)((ptgsHereE->TrajGraph)[cntr])*ftmp2/prob[1];
								break;
							case 'C':
								ctemp=(FloatLo)((ptgsHereC->TrajGraph)[cntr])*ftmp3/prob[2];
						}
					}
					else if (TrajType==TRAJ_STANDARD) {
						/* TRAJ_STANDARD */
						ctemp=(FloatLo)((ptgsHereH->TrajGraph)[cntr])+(FloatLo)((ptgsHereE->TrajGraph)[cntr])+(FloatLo)((ptgsHereC->TrajGraph)[cntr]);
					}
				}
				/* leave trajectory graph as PDF */	
				TrajGraph[cntr]=(Int4)floor(ctemp);
				/* try to remove effect of rounding a bit */
				remain+=ctemp-floor(ctemp);
				if (remain>1.0) {
					remain-=1.0;
					TrajGraph[cntr]++;
				}
				csum+=TrajGraph[cntr];
				cntr++;
			}
		TrajIntegral=(Int4)csum;
		ptgsHereH=FreeTraj(ptgsHereH);
		ptgsHereE=FreeTraj(ptgsHereE);
		ptgsHereC=FreeTraj(ptgsHereC);
		TGClose();
	}
	else if (TrajType==TRAJ_UNIFORM) {
		/* apply uniform trajectories of 1 */
		for (cntr=0;cntr<trajwidth*trajwidth;cntr++) {
			TrajGraph[cntr]=(Int4)1;
		}
	       	TrajIntegral=(Int4)(trajwidth*trajwidth);
	}
	return TrajIntegral;
}

/* add a pfds to end of a list, or start new list if null */
/* resnum is ignored except if head is NULL */
PFDS AddFragmentResidue(PFDS *pfdsHead,Int4 res,Char aa,FloatLo a1,FloatLo a2,FloatLo asd,FloatLo chiw,FloatLo chiwsd,Uint4 rotid,FloatLo p)
{
	PFDS pfdsHere,pfds;
	
	if (pfdsHead==NULL)
		return NULL;
	pfdsHere=(PFDS)MemNew(sizeof(FDS));
	pfdsHere->next=NULL;
	pfdsHere->prev=NULL;
	pfdsHere->nextfrag=NULL;
	/* no preservation of disulphides in fragments */
	if (aa!='C')
		pfdsHere->pSS=0.0;
	else
		pfdsHere->pSS=P_SS;
	pfdsHere->Angle1=a1;
	pfdsHere->Angle2=a2;
	pfdsHere->AngleSD=asd;
	pfdsHere->ChiWMean=chiw;
	pfdsHere->ChiWSD=chiwsd;
	pfdsHere->AA=aa;
	pfdsHere->reserved=0;	
	pfdsHere->length=0;
	/* # tries depends on standard deviations given */
	pfdsHere->tout=(Int4)(2.0*sqrt((FloatLo)chiwsd)*(FloatLo)asd);	
	if (pfdsHere->tout>BACKTRACK_TRIES_MAX)
		pfdsHere->tout=BACKTRACK_TRIES_MAX;
	if (pfdsHere->tout<BACKTRACK_TRIES_MIN)
		pfdsHere->tout=BACKTRACK_TRIES_MIN;
	pfdsHere->rotid=rotid;
	pfdsHere->prob=p;
	if (*pfdsHead==NULL) {
		*pfdsHead=pfdsHere;
		pfdsHere->resnum=res;
		return pfdsHere;
	}
	pfds=*pfdsHead;
	while (pfds->next!=NULL) {
		(pfds->length)++;
		pfds=pfds->next;
	}
	(pfds->length)++;
	pfds->next=pfdsHere;
	pfdsHere->prev=pfds;
	pfdsHere->resnum=(pfds->resnum)+1;
	return pfdsHere;
}

/* add the given fragment list to the given trajectory distribution */
TrajErr AddFragmentList(PTGS ptgsHere,PFDS pfdsHead)
{
	PFDS pfdsHere;
	
	if (pfdsHead==NULL)
		return ERR_FAIL;
	if (ptgsHere==NULL)
		return ERR_FAIL;
	if (ptgsHere->resnum!=pfdsHead->resnum) {
		ErrPostEx(SEV_ERROR,1,1,"Residue number in fragment does not match trajectory distribution");
		return ERR_FAIL;
	}
	if (WALKTYPE==WALK_CA && ptgsHere->resnum==1) {
		ErrPostEx(SEV_ERROR,1,1,"You may not add fragments to the first residue for a CA-walk trajectory distribution");
		return ERR_FAIL;
	}
	if (ptgsHere->pfdsFragmentHead==NULL) {
		ptgsHere->pfdsFragmentHead=pfdsHead;
		return ERR_SUCCESS;
	}
	pfdsHere=ptgsHere->pfdsFragmentHead;
	while (pfdsHere->nextfrag!=NULL)
		pfdsHere=pfdsHere->nextfrag;
	pfdsHere->nextfrag=pfdsHead;
	return ERR_SUCCESS;
}

TrajErr RemoveFragmentResidue(PTGS ptgsHere,PFDS pfdsHead,Int4 resnum)
{
	PFDS pfdsHere,pfdsNext;
	
	if (pfdsHead==NULL)
		return ERR_FAIL;
	if (ptgsHere==NULL)
		return ERR_FAIL;
	pfdsHere=pfdsHead;
	while (pfdsHere!=NULL && pfdsHere->resnum!=resnum) {
		(pfdsHere->length)--;
		pfdsHere=pfdsHere->next;
	}
	if (pfdsHere==NULL) {
		ErrPostEx(SEV_ERROR,1,1,"Residue %d not found in fragment",(int)resnum);
		pfdsHere=pfdsHead;
		/* undo length change */
		while (pfdsHere!=NULL && pfdsHere->resnum!=resnum) {
			(pfdsHere->length)++;
			pfdsHere=pfdsHere->next;
		}
		return ERR_FAIL;
	}
	if (pfdsHere->next==NULL && pfdsHere->prev==NULL) {
		/* remove whole list, last residue */
		return RemoveFragmentList(ptgsHere,pfdsHere);
	}
	if (pfdsHere->next!=NULL)
		pfdsHere->next->prev=pfdsHere->prev;
	if (pfdsHere->prev!=NULL)
		pfdsHere->prev->next=pfdsHere->next;
	pfdsNext=pfdsHere->next;
	MemFree(pfdsHere);
	/* renumber residues automagically */
	while (pfdsNext!=NULL) {
		(pfdsNext->resnum)--;
		pfdsNext=pfdsNext->next;
	}
	return ERR_SUCCESS;			
}

/* remove specified fragment from the set associated with the given
   trajectory distribution */
TrajErr RemoveFragmentList(PTGS ptgsHere,PFDS pfdsHead)
{
	PFDS pfdsHere;
	
	if (pfdsHead==NULL)
		return ERR_FAIL;
	if (ptgsHere==NULL)
		return ERR_FAIL;
	if (ptgsHere->pfdsFragmentHead==NULL)
		return ERR_FAIL;
	if (ptgsHere->pfdsFragmentHead==pfdsHead) {
		ptgsHere->pfdsFragmentHead=pfdsHead->nextfrag;
		pfdsHead->nextfrag=NULL;
		FreeFragmentList(pfdsHead);
		return ERR_SUCCESS;
	}				
	pfdsHere=ptgsHere->pfdsFragmentHead;
	while (pfdsHere->nextfrag!=NULL && pfdsHere->nextfrag!=pfdsHead)
		pfdsHere=pfdsHere->nextfrag;
	if (pfdsHere->nextfrag==NULL) {
		ErrPostEx(SEV_ERROR,1,1,"Fragment not found to delete");
		return ERR_FAIL;
	}
	pfdsHere->nextfrag=pfdsHead->nextfrag;
	pfdsHead->nextfrag=NULL;
	FreeFragmentList(pfdsHead);
	return ERR_SUCCESS;
}
		
FloatLo GetTotalFragmentProb(PTGS ptgsHere)
{
	PFDS pfds;
	FloatLo prob=0.0;

	if (ptgsHere==NULL)
		return prob;
	pfds=ptgsHere->pfdsFragmentHead;
	while (pfds) {
		prob=prob+pfds->prob;
		pfds=pfds->nextfrag;
	}
	return prob;
}

/* adds fragments from file fnam to currently open trajectory distribution
(you must provide the sequence as well, with X's for non-standard AAs) */
TrajErr ImportFrags(CharPtr fnam,CharPtr seq)
{
	Char path[PATH_MAX];
	Char buf[MAXCOL];
	Char resname[2];
	int istart,ilen;
	float fprob,phi,psi,sd,omega,omegasd,chi1,chi2,chi3,chi4;
	PFDS pfdsHere;
	PTGS ptgs;
	Int2 err,res,numc;
	Int4 line=0,rotid=0,numrecords;
	PRS prsHere;
	FILE *fp;
	Boolean ulr;
	static Int2 numchi[20]={0,4,2,2,1,3,3,0,2,2,2,4,3,2,2,1,1,2,2,1};

	if ((fp=FileOpen(fnam,"r"))==NULL)
		return ERR_FAIL;
	while (FileGets(buf,MAXCOL,fp)!=NULL) {
		line++;
		if (sscanf(buf,"%d %d %f",&istart,&ilen,&fprob)==3) {
			if (fprob<0.0) {
				fprob=0.0;
				ErrPostEx(SEV_ERROR,0,0,"Warning: Probability at line %d truncated to zero",line);
			}
			if (fprob>1.0) {
				fprob=1.0;
				ErrPostEx(SEV_ERROR,0,0,"Warning: Probability at line %d truncated to one",line);
			}
			res=istart-1;
			pfdsHere=NULL;
			while (res<istart+ilen-1) {
				if (FileGets(buf,MAXCOL,fp)!=NULL) {
					line++;
					chi1=0.0;
					chi2=0.0;
					chi3=0.0;
					chi4=0.0;
					err=sscanf(buf,"%f %f %f %f %f %f %f %f %f",&phi,&psi,&sd,&omega,&omegasd,&chi1,&chi2,&chi3,&chi4);
					if (err>4) {
						res++;
						/* found a valid fragment residue record - add it to fragment */
						/* validate data */
						if (phi<=-180.0 || phi>180.0)
						{
							while (phi>180.0)
								phi-=360.0;
							while (phi<=-180.0)
								phi+=360.0;
						}
						if (psi<=-180.0 || psi>180.0)
						{
							while (psi>180.0)
								psi-=360.0;
							while (psi<=-180.0)
								psi+=360.0;
						}
						if (sd<0.0) {
							sd=0.0;
							ErrPostEx(SEV_ERROR,0,1,"Warning: Standard Deviation at line %d truncated to 0",line);
						}
						if (sd>45.0) {
							sd=45.0;
							ErrPostEx(SEV_ERROR,0,1,"Warning: Standard Deviation at line %d truncated to 45",line);
						}
						if (omega<=-180.0 || omega>180.0)
						{
							while (omega>180.0)
								omega-=360.0;
							while (omega<=-180.0)
								omega+=360.0;
						}
						if (omegasd<0.0) {
							omegasd=0.0;
							ErrPostEx(SEV_ERROR,0,2,"Warning: Omega Standard Deviation at line %d truncated to 0",line);
						}
						if (omegasd>45.0) {
							omegasd=45.0;
							ErrPostEx(SEV_ERROR,0,2,"Warning: Omega Standard Deviation at line %d truncated to 45",line);
						}
						resname[0]=seq[res-1];
						resname[1]='\0';
						numc=numchi[StringCSpn(aalist,resname)];
						if (numc>0) {
							if (chi1<=-180.0 || chi1>180.0)
							{
								while (chi1>180.0)
									chi1-=360.0;
								while (chi1<=-180.0)
									chi1+=360.0;
							}
						}
						if (numc>1) {
							if (chi2<=-180.0 || chi2>180.0)
							{
								while (chi2>180.0)
									chi2-=360.0;
								while (chi2<=-180.0)
									chi2+=360.0;
							}
							if ((resname[0]=='D' || resname[0]=='F' || resname[0]=='Y') && (chi2<-90.0 || chi2>90.0)) {
								if (chi2<-90.0)
									chi2+=180.0;
								if (chi2>90.0)
									chi2-=180.0;
							}
						}
						if (numc>2) {
							if (chi3<=-180.0 || chi3>180.0)
							{
								while (chi3>180.0)
									chi3-=360.0;
								while (chi3<=-180.0)
									chi3+=360.0;
							}
						}
						if (numc>3) {
							if (chi4<=-180.0 || chi4>180.0)
							{
								while (chi4>180.0)
									chi4-=360.0;
								while (chi4<=-180.0)
									chi4+=360.0;
							}
						}

						prsHere=(PRS)MemNew(sizeof(RS));
						prsHere->Chi1=chi1;
						prsHere->Chi2=chi2;
						prsHere->Chi3=chi3;
						prsHere->Chi4=chi4;
						rotid=ComputeRotid(prsHere);
						prsHere=MemFree(prsHere);
						if (AddFragmentResidue(&pfdsHere,res,resname[0],phi,psi,sd,omega,omegasd,rotid,(FloatLo)fprob)==NULL) {
							ErrPostEx(SEV_ERROR,0,3,"The new residue of the fragment was not correctly added.");
						}
					}
				}
				else
					break;
			}
			if (pfdsHere) {
                ulr=USE_LOTS_RAM;
                USE_LOTS_RAM=FALSE;
				TGInit(tmpdbasename, DB_READ, &numrecords);
				ptgs =  TrajGraphRead((Int2) istart);
				if (AddFragmentList(ptgs,pfdsHere)!=ERR_SUCCESS)
					ErrPostEx(SEV_ERROR,0,4,"The new fragment was not correctly added.");
				TrajGraphWrite (ptgs, USE_RLE, TRUE);
				TGClose();
                USE_LOTS_RAM=ulr;
                ptgs=FreeTraj(ptgs);
			}
		}
	}
	FileClose(fp);
	return ERR_SUCCESS;
}

/*
$Log: trajtools.c,v $
Revision 1.79  2004/09/24 19:09:36  hfeldman
Added import fragment option to maketrj

Revision 1.78  2004/06/29 20:31:26  egarderm
Added CASP defines to correctly use the energy when evaluating structures

Revision 1.77  2003/11/18 19:52:59  egarderm
Fixed bracket typo

Revision 1.76  2003/11/18 19:45:12  egarderm
Changed fprintf error-checking to checking the fflfush return value (should not be EOF) before closing the file.

Revision 1.75  2003/11/13 22:26:57  feldman
Added a few extra error exit paths

Revision 1.74  2003/11/13 21:33:26  egarderm
Added error-handling for problematic file opening and writing when writing the trajectory distribution to the database or when dumping the log to disk

Revision 1.73  2003/10/23 17:30:35  feldman
Score by best crease energy column now

Revision 1.72  2003/09/18 18:28:58  feldman
Built in workaround for codebase -490 corrupt index problem - will rebuild index automagically when this happens now
(thanks to codebase tech support)

Revision 1.71  2003/09/14 02:38:23  feldman
Fixed unused variables and other minor compiler warnings

Revision 1.70  2003/08/05 16:12:17  feldman
Arrgh, cvs killed my last change.  here it is again

Revision 1.69  2003/08/05 16:08:45  feldman
made dfptmpnam more rigorous still

Revision 1.66  2003/07/17 13:31:49  feldman
Made universal DFPTmpNam function

Revision 1.65  2003/07/15 19:35:50  feldman
Made temp file naming more robust, especially if using DFPTEMP

Revision 1.64  2003/04/17 16:33:09  feldman
checks for gen 0 in DumpLog (for possible future use)

Revision 1.63  2003/03/14 21:06:01  feldman
Now turning off USE_LOTS_RAM during writing on TDs to avoid loading useless stuff into RAM

Revision 1.62  2003/03/11 19:31:24  feldman
changed scoring f'n back to RMS for now

Revision 1.61  2003/02/28 20:03:32  feldman
Log best energy not RMSD now

Revision 1.60  2003/02/27 21:20:33  feldman
Added parameter to dumplog

Revision 1.59  2003/01/24 16:52:07  feldman
Gets rid of junk at end of log in case of exiting midway

Revision 1.58  2003/01/09 20:01:04  feldman
Fixed bug for uselotsram

Revision 1.57  2002/10/16 23:05:40  feldman
Fixed stupid UNIX crashing bug

Revision 1.56  2002/10/10 20:32:32  feldman
New DFP clinet log format (only applies when bzipit is TRUE)

Revision 1.55  2002/08/30 01:42:51  feldman
Fixed memory violation with DFPTEMP

Revision 1.54  2002/08/11 15:42:37  feldman
Added DFPTEMP variable for changing temp dir in DFP

Revision 1.53  2002/07/25 16:31:54  feldman
Added tunnel prob

Revision 1.52  2002/07/20 16:42:43  feldman
Added LOTS_RAM option for trajgraphreading

Revision 1.51  2002/07/14 21:48:29  feldman
Added get cis function

Revision 1.50  2002/07/05 16:10:05  feldman
Added f'n to count total frag probability

Revision 1.49  2002/06/28 16:28:08  feldman
More Windows tmpnam fixes

Revision 1.48  2002/06/14 20:47:17  feldman
Added fragment length to data structures

Revision 1.47  2002/02/25 22:08:01  feldman
Added checksumming of text files potentials, cbdata and skel.prt
If checksum fails, get error
Changed bailing from foldtrajlite to exit curses first if in it

Revision 1.46  2002/02/07 21:21:07  feldman
Fixed Darwin fix

Revision 1.45  2001/12/06 23:31:07  feldman
Fix in tmpnaming for Darwin

Revision 1.44  2001/08/31 17:43:40  feldman
Now we pack the memo field of trajectory distributions before closing,
to remove wasted space

Revision 1.43  2001/07/31 19:45:51  feldman
minor changes

Revision 1.42  2001/07/26 19:19:52  feldman
Changed way fragments are stored on disk to be more platform independent,
and removed id parameter from AddtoWorld

Revision 1.41  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.40  2001/04/17 21:32:06  feldman
Added fragment editing to Vistraj - not quite fully functional yet
Also can now give trajectory graph name as argument #1 when running
and fixed potential sscanf problems so entering garbage into numeric
entry fields results in an error being shown (instead of accepting it)

Revision 1.39  2001/04/06 14:29:00  feldman
Completed fold-at-home Client-Server v1.0 now complete and ready for testing

Revision 1.38  2001/04/05 14:02:39  feldman
Fixed so curses should work on both HP and Alpha now, I hope.
Also added a few more error codes to trajstore

Revision 1.37  2001/04/04 21:26:01  feldman
-Removed my BSPrintf, now using the one in SLRIlib - thus needed to add
	this library to makefiles
-changed all fgets to FileGets

Revision 1.36  2001/03/30 22:22:27  feldman
Minor spelling and grammar fixes, bug fixes, and added
foldtraj ASCII screensaver ability

Revision 1.35  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.34  2001/03/29 02:52:24  feldman
fixed minor compiler warnings

Revision 1.33  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.32  2001/03/23 21:08:07  feldman
Added monitor for foldtraj inside vistraj

Revision 1.31  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.30  2001/03/13 15:08:17  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.29  2001/03/12 17:11:14  feldman
fixed a compiler warning

Revision 1.28  2001/02/27 20:52:32  feldman
minor bugfix

Revision 1.27  2001/02/26 22:21:15  feldman
-Changed fragments to allow multiple possible fragments at the same
residue
-altered random walk to make use of fragments when present (a first
attempt at least...)

Revision 1.26  2001/02/09 20:18:13  feldman
Added Fragments to trajectory graph structure

Revision 1.25  2001/02/07 18:46:45  feldman
Removed a few unused variables

Revision 1.24  2001/02/06 18:40:39  feldman
Added a few functions for dealing with distance constraints and
tidied up so maketrj could join the library (foldtrajlib) without
conflicts

Revision 1.23  2001/01/26 20:37:38  feldman
Added getconstraints function

Revision 1.22  2001/01/24 18:32:03  feldman
Some constraint parameter tweaking and removed bz_internal_error
function (put in bzlib_private.h)

Revision 1.21  2001/01/18 18:10:38  feldman
Added bzip2 1.01 compatibility

Revision 1.20  2001/01/16 22:01:36  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.19  2001/01/12 20:01:51  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.18  2000/09/15 20:33:22  feldman
Count distance constraint violations now
Allow constraint file to have angles and dihedrals along
with distances - directed distance constraint

Revision 1.17  2000/08/31 14:54:40  feldman
fixed filter functions to take filt rather than *filt

Revision 1.16  2000/08/25 18:11:09  feldman
Corrected an error message

Revision 1.15  2000/08/18 19:06:05  adrian
moved r14seqfilecalc to r14bscalc and moved back into r14.c where it belongs
fixed sequence numbering off-by one in output (still needs to be checked elsewhere)
more debug information added to RLEunpack error message in trajtools

Revision 1.14  2000/08/09 19:12:19  feldman
-minor bugfix update and fixed up makefiles removing USEDSSP
-added bioseq/seq-entry to trajectory graph header when using
 unfoldtraj or val2trj

Revision 1.13  2000/07/27 18:24:47  feldman
Removed all MOBI references in makefiles
Added new method to ensembletrj - straight summation of all positive
peaks

Revision 1.12  2000/07/14 20:09:25  feldman
Added parameter for MIMEBiostrucAsnGet for Bioseq

Revision 1.11  2000/07/13 20:41:06  feldman
Added inclusive/exclusive windows for crease energy and corrected
makeTG wrapping for Phi-Psi space

Revision 1.10  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.9  2000/07/10 15:40:58  feldman
Updated TrajgraphWrite to not call TGInit and TGClose, thus
now just call TGInit once at start of program, TGClose at end
(also removed these calls from updatePcis and alterresiduesequence
)

Revision 1.8  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.7  2000/07/06 15:29:38  feldman
-- Updated old makefiles to compile with newest toolkit and
directory structure
-- added needed functions to hfprogs.h
-- update some executables to include hfprogs.h instead of mmdbtraj.h
-- replaced all instances of BiostrucAsnGet with MIMEBiostrucAsnGet
which can read with MIME biostrucs (v2.0) or normal ones (v1.0), and
as a result, foldtraj uses a new skel.prt, to result in a MIME biostruc
as its initial input

Revision 1.6  2000/07/04 17:00:27  feldman
Fixed bug in maketrj which made residue one have a zero trajectory
graph

Revision 1.5  2000/06/20 18:59:46  feldman
moved to library directory and changed include files appropriately

Revision 1.4  2000/06/20 16:40:23  feldman
Incorporated sstru extended structure calculation into C code,
no longer need an external program

Revision 1.3  2000/06/15 17:11:03  feldman
Added Replace option when TrajGraphWrite-ing

Revision 1.2  2000/06/15 15:52:11  feldman
Corrected Bzip2 log dumping and improved H-bond treatment (fixed
bug that was crashing program)

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

