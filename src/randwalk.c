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
#ifdef USE_DSSP
#include <slriaccsurf.h>
#endif
#include <geometry.h>
#include <slri_misc.h>


/* do not change -- used for position Cb */
#define CB_PRECISION 0.00001  
/* can be used to give less (or more) importance to backbone
   error involving the C-Beta atoms */
#define FUDGEFACT 1.0	
/* distance bins in C-Beta lookup table */
#define NUMBINS 5		
/* index into several length and angle arrays for Cystine */
#define CYSSYC 26		
/* index into length and angle arrays for Cis-Proline */
#define CISPRO 27		
/* avoids possible but rare endless loop when searching for minima in backbone error surface */
#define TRIES_MAX 10000		
#define HYP_DICT_IDX 210
/* after doing this fraction of tries to place a rotamer with given
chi angles, and failing, it will switch to the rotamer library; the
only this that changes between tries before this is the exact chi due
to the standard deviation, and the random acceptance probabilities
for van der waals collisions */
#define CHANCE_TO_USE_FIXED_ROTAMER 0.5
/* standard deviation of chi1 and chi2 when given explicitly, to allow
a bit of looseness */
#define FIXED_ROTAMER_SD 5.0 /* degrees */


#define CHECKHBONDS					\
	if (IsAtomPlaced(pmad,1)) {	\
		pvnal=pmad->pvnalLocate;	\
		paldHbond=(PALD)(pvnal->data.ptrvalue);	\
		vHbond[0]=(paldHbond->pflvData)[0];	\
		vHbond[1]=(paldHbond->pflvData)[1];	\
		vHbond[2]=(paldHbond->pflvData)[2];	\
		VecSub(vTmp1,vAtom1,vHbond);		\
		VecSub(vTmp2,vAtom2,vHbond);		\
		Normalize(vTmp1,vTmp1);			\
		Normalize(vTmp2,vTmp2);			\
		HbondAngleCos=Dot(vTmp1,vTmp2);		\
		if (HbondAngleCos<HBOND_ANGLE_MAX)	\
			return TRUE;			\
	}						\
	else {						\
		if (!PlacedOnly)			\
			return TRUE;			\
	}

#define CLEARFRAGS						\
		if (inafrag==TRUE) {				\
			inafrag=FALSE;				\
			pfdsHead=FreeFragmentList(pfdsHead);	\
			pfdsHere=NULL;				\
		}

/* global look-up tables */
static FloatLo cbdir[30][NUMBINS+1][3];
static FloatLo cbcisdir[30][3];
static FloatLo coszeta[30];
static FloatLo coseta[30];
static FloatLo bl_cao[30];
static FloatLo ba_ocaca[30];
/* ordinal AA list */
extern CharPtr NCBIstdaaUC;
/* = "-ABCDEFGHIKLMNPQRSTVWXYZU*" */
/* more look-up tables */
extern FloatLo ba_ncacb[];
extern FloatLo ba_ccacb[];
extern FloatLo bl_cacb[];
extern Boolean volatile timetoquit;
static PWS pwsThis;
static Int2 lrx;
static Int4 lrxtries,lastprt=0;
/* store trajectory graphs */
static PTGS ptgsHere;
/* parameters read from parameter file */
static Int4 tries,fpabad,crashcnt,dcbad,dcviol;
/* used for GOR structure prediction scheme */
static Int2 IsCis[MAXRES+1];
/* some global co-ordinate vectors to minimize parameter passing */
static vec vCALast,vCAHere,vCBHere,vCLast,vNHere;
/* maximum side chain length 30 */
static vec vRef[30];
static vec vZero={0.0,0.0,0.0};
static Int2 cisLast,cisHere,cisNext;
/* timer */
static time_t tm1,tm2;
/* 0! up to 17! */
static const FloatHi factorial[18]={1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0,3628800.0,39916800.0,479001600.0,6227020800.0,87178291200.0,1307674368000.0,20922789888000.0,355687428096000.0};
static Char scpos[]="BGDEZHQILMXRTUFYWJKV";
static Char sccollide[]="GDEZHQILMXRTUFYWJKV";
static Char smallatoms[]="EZHQILMXRTUFYWJKV";
static PHBS phbsHBondsToCheck=NULL;
static Uint4 rotidHere=0,rotidNext=0;
static Int2 piclastx,piclasty;
static Int2 pictop,picbot,piclft,picrgt;
static Boolean picdir_n,picdir_s,picdir_w,picdir_e;
static Char tmp_str[PATH_MAX];
static Int2 lastresdrawn;
#ifdef OS_UNIX
static Boolean pichascolor=FALSE;
#endif
static FloatHi gencbrt=0.0;
static Int2 *ressrcres=NULL;
static Int2 *ressrcfrag=NULL;
static Char errorfile[PATH_MAX];





ValNodePtr BackupCoords(PMMD pmmdHere,Int2 Model)
{
	ValNodePtr vnpHead=NULL,vnpLast;
	PDNMG pdnmgHere;
	PVNMA pvnmaHere;
	PMAD pmadHere;
	PALD paldHere;
	FloatLoPtr fpAtom;
	
	pdnmgHere=pmmdHere->pdnmgHead;
	while (pdnmgHere!=NULL) {
		pvnmaHere=(((PMGD)(pdnmgHere->data.ptrvalue))->pvnmaAHead);
		while (pvnmaHere!=NULL) {
			pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
			paldHere=GetAtomLocs(pmadHere,Model);
			if (paldHere!=NULL) {
			  fpAtom=(FloatLoPtr)MemNew(3*sizeof(FloatLo));
				fpAtom[0]=AtomLocX(paldHere);
				fpAtom[1]=AtomLocY(paldHere);
				fpAtom[2]=AtomLocZ(paldHere);
				if (vnpHead==NULL)
					vnpLast=ValNodeAddPointer(&vnpHead,0,fpAtom);
				else
					vnpLast=ValNodeAddPointer(&vnpLast,0,fpAtom);								
			}
			pvnmaHere=pvnmaHere->next;
		}	
		pdnmgHere=pdnmgHere->next;
	}	
	return vnpHead;
}

void RestoreCoords(PMMD pmmdHere,Int2 Model,ValNodePtr vnpCoords)
{
	ValNodePtr vnpHere;
	PDNMG pdnmgHere;
	PVNMA pvnmaHere;
	PMAD pmadHere;
	PALD paldHere;
	FloatLoPtr fpAtom;
	FILE *fp;
	
	pdnmgHere=pmmdHere->pdnmgHead;
	vnpHere=vnpCoords;
	while (pdnmgHere!=NULL) {
		pvnmaHere=(((PMGD)(pdnmgHere->data.ptrvalue))->pvnmaAHead);
		while (pvnmaHere!=NULL) {
			pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
			paldHere=GetAtomLocs(pmadHere,Model);
			if (paldHere!=NULL) {
				if (vnpHere==NULL) {
					ErrPostEx(SEV_FATAL,0,0,"Error restoring atom co-ordinates, cannot continue");
					return;
				}
			  fpAtom=(FloatLoPtr)(vnpHere->data.ptrvalue);
			  (paldHere->pflvData)[0]=fpAtom[0];
			  (paldHere->pflvData)[1]=fpAtom[1];
			  (paldHere->pflvData)[2]=fpAtom[2];
			  vnpHere=vnpHere->next;
			}
			pvnmaHere=pvnmaHere->next;
		}	
		pdnmgHere=pdnmgHere->next;
	}	
}

/* finds neighbour PMAD to an H atom, used in H-bonding */
PMAD FindHNeighbour(PMAD pmadHere)
{
	PMBD pmbdHere;
	ValNodePtr vnpHere;

	/* only one single bond from an H */
	vnpHere=pmadHere->pvnBonds;
	do {
		pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
		vnpHere=vnpHere->next;
	} while (((pmbdHere->bWhat) & BOND_SINGLE)==0);
	if (pmbdHere->pmadFrom==pmadHere)
		return pmbdHere->pmadTo;
	return pmbdHere->pmadFrom; 	
}

/* returns TRUE if a H-bond may exist between two colliding atoms */
Boolean FindHBond(PMAD pmad1,PMAD pmad2,Boolean PlacedOnly,Boolean at1new)
{
	vec vHbond,vAtom1,vAtom2,vTmp1,vTmp2;
	PALD	pald1,pald2,paldHbond;
	PVNAL pvnal;
	FloatLo HbondAngleCos;
	ValNodePtr vnpHere;
	PMAD pmad;
	PMBD pmbdHere;
	Char res1,res2;
	Char at1[5],at2[5];
	
	if (pmad1==pmad2)
		return FALSE;
	if ((pmad1->pvnalLocate)==NULL)
		return FALSE;
	if ((pmad2->pvnalLocate)==NULL)
		return FALSE;
	if (!IsAtomPlaced(pmad2,1))
		return FALSE;
	if (!IsAtomPlaced(pmad1,1) && !at1new)
		return FALSE;
	res1=(((PMGD)(pmad1->pfbParent))->pcIUPAC)[0];
	res2=(((PMGD)(pmad2->pfbParent))->pcIUPAC)[0];
	StringCpy(at1,pmad1->pcAName);
	StringCpy(at2,pmad2->pcAName);
	pald1=(PALD)((pmad1->pvnalLocate)->data.ptrvalue);
	vAtom1[0]=(pald1->pflvData)[0];
	vAtom1[1]=(pald1->pflvData)[1];
	vAtom1[2]=(pald1->pflvData)[2];
	pald2=(PALD)((pmad2->pvnalLocate)->data.ptrvalue);
	vAtom2[0]=(pald2->pflvData)[0];
	vAtom2[1]=(pald2->pflvData)[1];
	vAtom2[2]=(pald2->pflvData)[2];
	vnpHere=pmad1->pvnBonds;
	while (vnpHere) {
		pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
		/* remember H on Glu and Asp are never placed so ignore them */
		if (pmbdHere->pmadFrom->pcAName[1]=='H' && ((res1!='D' && res1!='E') || at1[1]!='O' || (at1[2]!='D' && at1[2]!='E')) && (res1!='H' || StringCmp(at1," ND1")) && StringCmp(at1," OXT")) {
			/* check if h-bond possible */
			pmad=pmbdHere->pmadFrom;
			/* save some typing here */
			CHECKHBONDS				
		}
		if (pmbdHere->pmadTo->pcAName[1]=='H' && ((res1!='D' && res1!='E') || at1[1]!='O' || (at1[2]!='D' && at1[2]!='E')) && (res1!='H' || StringCmp(at1," ND1")) && StringCmp(at1," OXT")) {
			/* check if h-bond possible */
			pmad=pmbdHere->pmadTo;
			CHECKHBONDS
		}
		vnpHere=vnpHere->next;
	}
	vnpHere=pmad2->pvnBonds;
	while (vnpHere) {
		pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
		if (pmbdHere->pmadFrom->pcAName[1]=='H' && ((res2!='D' && res2!='E') || at2[1]!='O' || (at2[2]!='D' && at2[2]!='E')) && (res2!='H' || StringCmp(at2," ND1")) && StringCmp(at2," OXT")) {
			/* check if h-bond possible */
			pmad=pmbdHere->pmadFrom;
			CHECKHBONDS
		}
		if (pmbdHere->pmadTo->pcAName[1]=='H' && ((res2!='D' && res2!='E') || at2[1]!='O' || (at2[2]!='D' && at2[2]!='E')) && (res2!='H' || StringCmp(at2," ND1")) && StringCmp(at2," OXT")) {
			/* check if h-bond possible */
			pmad=pmbdHere->pmadTo;
			CHECKHBONDS
		}
		vnpHere=vnpHere->next;
	}
	return FALSE;
}

/* counts neighbour PMAD which are H atom, used in H-bonding */
/* returns two-digit number, 1st is # bonded atoms and 2nd is
   # bonded hydrogen atoms */
Int2 CountHNeighbour(PMAD pmadHere)
{
	PMBD pmbdHere;
	ValNodePtr vnpHere;
	Int2 neighbours,numbonds;
	Char res;

	neighbours=0;
	numbonds=0;
	/* only one single bond from an H */
	if (pmadHere->pcAName[1]=='H')
		return 10;
	res=(((PMGD)(pmadHere->pfbParent))->pcIUPAC)[0];
	if (!StringCmp(pmadHere->pcAName," OXT"))
		return 10;
	if (res=='H' && !StringCmp(pmadHere->pcAName," ND1"))
		return 20;
	if (pmadHere->pcAName[1]=='O' && (pmadHere->pcAName[2]=='E' || pmadHere->pcAName[2]=='D') && (res=='D' || res=='E'))		
		return 10;	
	vnpHere=pmadHere->pvnBonds;	
	while (vnpHere) {
		pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
		if (pmbdHere->pmadFrom->pcAName[1]=='H')
			neighbours++;
		if (pmbdHere->pmadTo->pcAName[1]=='H')
			neighbours++;
		numbonds++;
		vnpHere=vnpHere->next;
	}
	return neighbours+numbonds*10;
}

PHBS PHBSFree(PHBS phbs)
{
	PHBS phbsNext;
	
	while (phbs!=NULL) {
		phbsNext=phbs->next;
		MemFree(phbs);
		phbs=phbsNext;
	}
	return NULL;
}

void AddHBondToCheck(PMAD pmadDonor,PMAD pmadAcceptor,Boolean certain)
{
	PHBS phbs;
	
	phbs=(PHBS)MemNew(sizeof(HBS));
	phbs->pmadDonor=pmadDonor;
	phbs->pmadAcceptor=pmadAcceptor;
	phbs->certain=certain;
	phbs->next=phbsHBondsToCheck;
	if (phbs->next!=NULL)
		phbs->next->prev=phbs;
	phbs->prev=NULL;
	phbsHBondsToCheck=phbs;
}

static void LIBCALLBACK FreeUnwantedLocs(PFB curpfb,Int4 Model,Int4 dum1,Pointer dum2)
{
	PMGD pmgd;

	/* destroy co-ordinates */
	pmgd=(PMGD)(curpfb->pfbParent);
	if (!ressrcres[pmgd->pdnmgLink->choice-1]) {
		FreeListVNAL(((PMAD)curpfb)->pvnalLocate);
		((PMAD)curpfb)->pvnalLocate=NULL;
	}
}

void RemoveHBondToCheck(PMAD pmadDonor,PMAD pmadAcceptor)
{
	PHBS phbs,phbsNext;
	
	phbs=phbsHBondsToCheck;
	while (phbs!=NULL) {
		if (phbs->pmadDonor==pmadDonor && phbs->pmadAcceptor==pmadAcceptor) {
			if (phbs->prev==NULL) {
				phbsHBondsToCheck=phbs->next;
				if (phbsHBondsToCheck!=NULL)
					phbsHBondsToCheck->prev=NULL;
			}
			else if (phbs->next==NULL) {
				phbs->prev->next=NULL;
			}
			else {
				phbs->prev->next=phbs->next;
				phbs->next->prev=phbs->prev;
			}
			phbsNext=phbs->next;
			MemFree(phbs);
			phbs=phbsNext;
		}
		else
			phbs=phbs->next;
	}
}

void MakeHBondToCheckCertain(PMAD pmadDonor,PMAD pmadAcceptor)
{
	PHBS phbs;
	
	phbs=phbsHBondsToCheck;
	while (phbs!=NULL) {
		if (phbs->pmadDonor==pmadDonor && phbs->pmadAcceptor==pmadAcceptor)
			phbs->certain=TRUE;
		phbs=phbs->next;
	}
}

FloatLo GetRadNewAtom(Char res,CharPtr atnam,Int2 athn)
{
	FloatLo radNewAtom=0.0;
	
	/* get atom VDW radii for AtomName */
	switch (atnam[1]) {
		case 'C':
			if (athn==30)
				radNewAtom=VDWRAD_AROCH0;
			else if (athn==31)
				radNewAtom=VDWRAD_AROCH1;
			else
				radNewAtom=VDWRAD_C;
			break;
		case 'N':
			radNewAtom=VDWRAD_N;
			break;
		case 'O':
			radNewAtom=VDWRAD_O;
			break;
		case 'S':
			radNewAtom=VDWRAD_S;
			break;
		case 'P':
			radNewAtom=VDWRAD_P;
			break;
		case 'E':
			/* Selenium */
			if (atnam[0]=='S')
				radNewAtom=VDWRAD_SE;
			else
				radNewAtom=VDWRAD_C;
			break;
		case 'H':
			radNewAtom=VDWRAD_H;
			break;
		default:
			radNewAtom=VDWRAD_C;
	}
	if ((StringRChr(smallatoms,(int)(atnam[2]))!=NULL) && !(atnam[2]=='X' && atnam[3]=='T') && ((res=='E') || (res=='Q') || (res=='K') || (res=='R'))) {
		if (radNewAtom>VDWRAD_X)
			radNewAtom=VDWRAD_X;
	}
	if ((atnam[2]=='D') && (res=='P')) {
		if (radNewAtom>VDWRAD_X)
			radNewAtom=VDWRAD_X;
	}
	return radNewAtom;
}

/* give co-ordinates vCoord to AtomName in pmgdHere, Model; this routine tests for
   van der Waals clashes and H-bonds as well */
TrajErr AssignCoords(CharPtr AtomName,PMGD pmgdHere,vec vCoord,Int2 Model)
{
  PMAD pmadHere,pmadHBond,pmadTemp;
  PALD paldNew,paldHere;
	PVNAL pvnalHere;
	PMGD pmgdCrash;
	ValNodePtr vnpAtomList,vnpHere;
	PHBS phbsHere,phbsNext;
	FloatLo radNewAtom=0.0,radHere=0.0,tooclose=0.0,bumpdist,randbump,bounciness;
	CharPtr AtomHere;
	vec vHere,vDist;
	Int2 resdist,doit,isProNew,isProCrash;
	Char resHere,resCrash;
	Int2 at1hn,at2hn;
	Boolean HBondFound,marked,balreadyadded;
	FILE *fp;

/*printf("Adding %s %d...\n",AtomName,(pmgdHere->pdnmgLink)->choice);*/
	isProNew=0;
	resHere=GetAAFromIDict(pmgdHere);
	pmadHere=FindAtomName(pmgdHere->pvnmaAHead,AtomName);
	if (pmadHere==NULL) {
		ErrPostEx(SEV_FATAL,1,10,"Unable to find atom %s in dictionary entry %d",AtomName,pmgdHere->iIDict);
	}
	/* check distance contraints before testing for collisions */
	at1hn=CountHNeighbour(pmadHere);
	if (CheckAtomDistConstraints(AtomName,(pmgdHere->pdnmgLink)->choice,vCoord,Model,&dcviol)!=ERR_SUCCESS) {
		dcbad++;
		return ERR_DISTCONST;
	}
	/* allocate space for co-ordinates */
	pvnalHere=pmadHere->pvnalLocate;
	if (pvnalHere!=NULL)
		paldNew=(PALD)(pvnalHere->data.ptrvalue);
	else
		paldNew=AllocNewLoc(pmadHere,Model);
  if (SetCoOrds(paldNew,vCoord)==ERR_FAIL) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,0,0,"Out of memory while building?");
	}
	if (resHere=='P')
		isProNew=1;
	radNewAtom=GetRadNewAtom(resHere,AtomName,at1hn);
	if (AtomName[1]=='H' && ATOM_BOUNCINESS_BB<1.0 && ATOM_BOUNCINESS_SC<3.0) {
		/* check if we should be making H bonds */
		/* don't bother checking H-bonds if bounciness is set high - doesn't really make sense */
		pmadHBond=FindHNeighbour(pmadHere);
		phbsHere=phbsHBondsToCheck;
		marked=FALSE;
		while (phbsHere!=NULL) {
			phbsNext=phbsHere->next;
			if (pmadHBond==phbsHere->pmadDonor) {
				if (FindHBond(pmadHBond,phbsHere->pmadAcceptor,TRUE,FALSE)==TRUE) {
					/* H-bond constraint met */
					if (phbsHere->certain==FALSE) {
						/* donor then acceptor - i.e. we remove the node where acceptor here is the donor and vice-versa */
						RemoveHBondToCheck(phbsHere->pmadAcceptor,phbsHere->pmadDonor);
						/* in case this changes */
						phbsNext=phbsHere->next;
					}
					if (phbsHere->prev==NULL) {
						phbsHBondsToCheck=phbsHere->next;
						if (phbsHBondsToCheck!=NULL)
							phbsHBondsToCheck->prev=NULL;
					}
					else if (phbsHere->next==NULL) {
						phbsHere->prev->next=NULL;
					}
					else {
						phbsHere->prev->next=phbsHere->next;
						phbsHere->next->prev=phbsHere->prev;
					}
		      MemFree(phbsHere);
				}
				else {
					/* no h-bond, but this is supposed to be an h-bond donor */
					if (((pmadHBond->bReserved)&0x08) && !marked) {
						/* clear bit to indicate 1 less hydrogen to come */
						(pmadHBond->bReserved)&=0xF7;
						marked=TRUE;
					}
					else if (((pmadHBond->bReserved)&0x04) && !marked) {
						/* clear bit to indicate 1 less hydrogen to come */
						(pmadHBond->bReserved)&=0xF3;
						marked=TRUE;
					}
					else if ((((pmadHBond->bReserved)&0x02) && !marked) || (((pmadHBond->bReserved)&0x0E)==0)) {
						/* note could get here if donating to two acceptors - latter case means no more H's coming */
						/* clear bit to indicate 1 less hydrogen to come */
						(pmadHBond->bReserved)&=0xF1;
						marked=TRUE;
						/* last H on this pmad */
						if (phbsHere->certain==FALSE) {
							/* we'll remove this one and make the other requirement certain now */
							/* donor first, then acceptor - this will affect OTHER node, not this one */
				      MakeHBondToCheckCertain(phbsHere->pmadAcceptor,phbsHere->pmadDonor);
							if (phbsHere->prev==NULL) {
								phbsHBondsToCheck=phbsHere->next;
								if (phbsHBondsToCheck!=NULL)
									phbsHBondsToCheck->prev=NULL;
							}
							else if (phbsHere->next==NULL) {
								phbsHere->prev->next=NULL;
							}
							else {
								phbsHere->prev->next=phbsHere->next;
								phbsHere->next->prev=phbsHere->prev;
							}
				      MemFree(phbsHere);
			      }
			      else {
			      	/* twas certain, so we must backtrack */
							crashcnt++;
							return ERR_CRASH;
			      }
					}
				}
	   	}
	   	phbsHere=phbsNext;
		}
	}
	/* check for collisions whenever a new atom is placed */
	/* assume no atoms larger than radius of VDWRAD_MAX */
  vnpAtomList=FindAtomsIn(pwsThis,vCoord,radNewAtom+VDWRAD_MAX);
  vnpHere=vnpAtomList;
  /* walk through list to output results */
  while(vnpHere) {
  	/* recompute radNewAtom in case H-bonded with something which reduced its effective radius */
		radNewAtom=GetRadNewAtom(resHere,AtomName,at1hn);
	    paldHere=(PALD)(vnpHere->data.ptrvalue);
		pmgdCrash=GetParentGraph((PFB)paldHere);
		isProCrash=0;
		resCrash=GetAAFromIDict(pmgdCrash);
		if (resCrash=='P')
			isProCrash=1;
		/* check if different residues between colliding atoms */
		resdist=(((pmgdCrash->pdnmgLink)->choice)-((pmgdHere->pdnmgLink)->choice));
		if (abs(resdist)>0) {
			AtomHere=((PMAD)(paldHere->pfbParent))->pcAName;
			at2hn=CountHNeighbour((PMAD)(paldHere->pfbParent));
			/* NH and CO of peptide bond are allowed to "clash" though */
			/* certain atoms in adjacent residues are allowed to crash */
			doit=0;
			if (resdist==-1) {
				if (((!StringCmp(AtomHere," CB ")) || (!StringCmp(AtomHere,"2HA "))) && (
/*				(!StringCmp(AtomName," H  ")) ||*/
				(!StringCmp(AtomName," N  "))));
				else if ((!StringCmp(AtomHere," CA ")) &&
				((!StringCmp(AtomName," H  ")) || ((
				(!StringCmp(AtomName," CD ")) ||
				(!StringCmp(AtomName,"1HD ")) ||
				(!StringCmp(AtomName,"2HD "))) && isProNew) ||
				(!StringCmp(AtomName," N  "))));
				else if (((!StringCmp(AtomHere," HA ")) || (!StringCmp(AtomHere,"1HA "))) &&
				isProNew &&
				((!StringCmp(AtomName," CD ")) ||
				(!StringCmp(AtomName,"1HD ")) ||
				(!StringCmp(AtomName,"2HD "))));
				else if (((!StringCmp(AtomHere," HB ")) ||
				(!StringCmp(AtomHere,"1HB ")) ||
				(!StringCmp(AtomHere,"2HB "))) &&
				(!StringCmp(AtomName," H  ")));
				else if ((!StringCmp(AtomHere," C  ")) && (
				(!StringCmp(AtomName," H  ")) ||
				(((!StringCmp(AtomName,"1HD ")) ||
				(!StringCmp(AtomName,"2HD ")) ||
				(!StringCmp(AtomName," CD "))) && isProNew) ||
				(!StringCmp(AtomName," C  ")) ||
				(!StringCmp(AtomName," CA ")) ||
				(!StringCmp(AtomName," CB ")) ||
				(!StringCmp(AtomName,"2HA ")) ||
				(!StringCmp(AtomName," HA ")) ||
				(!StringCmp(AtomName,"1HA ")) ||
				(!StringCmp(AtomName," N  "))));
				else if ((!StringCmp(AtomHere," O  ")) && (
/*				(!StringCmp(AtomName," CB ")) ||*/
				(!StringCmp(AtomName," CA ")) ||
				(!StringCmp(AtomName," N  ")) ||
				(!StringCmp(AtomName," HA ")) ||
/*				(!StringCmp(AtomName,"2HA ")) ||*/
/*				(!StringCmp(AtomName," C  ")) ||*/
			  (!StringCmp(AtomName,"1HA "))));
				else if (isProCrash && (
				(!StringCmp(AtomName," N  ")) ||
				(!StringCmp(AtomName," N  "))
)) /*printf("excused\n")*/;
				else if ((!StringCmp(AtomHere," CA ")) && (!StringCmp(AtomName," CA ")));
				else doit=1;
			}
			else if (resdist==1) {
				if (((!StringCmp(AtomName," CB ")) ||
				(!StringCmp(AtomName,"2HA "))) &&
				(/*(!StringCmp(AtomHere," H  ")) ||*/
				(!StringCmp(AtomHere," N  "))));
				else if ((!StringCmp(AtomName," CA ")) &&
				((!StringCmp(AtomHere," H  ")) ||
				(((!StringCmp(AtomHere," CD ")) ||
				(!StringCmp(AtomHere,"1HD ")) ||
				(!StringCmp(AtomHere,"2HD "))) &&
				isProCrash) ||
				(!StringCmp(AtomHere," N  "))));
				else if (((!StringCmp(AtomName," HA ")) ||
				(!StringCmp(AtomName,"1HA "))) &&
				isProCrash &&
				((!StringCmp(AtomHere," CD ")) ||
				(!StringCmp(AtomHere,"1HD ")) ||
				(!StringCmp(AtomHere,"2HD "))));
				else if (((!StringCmp(AtomName," HB ")) ||
				(!StringCmp(AtomName,"1HB ")) ||
				(!StringCmp(AtomName,"2HB "))) &&
				(!StringCmp(AtomHere," H  ")));
				else if ((!StringCmp(AtomName," C  ")) &&
				((!StringCmp(AtomHere," H  ")) ||
				(((!StringCmp(AtomHere,"1HD ")) ||
				(!StringCmp(AtomHere,"2HD ")) ||
				(!StringCmp(AtomHere," CD "))) &&
				isProCrash) ||
				(!StringCmp(AtomHere," C  ")) ||
				(!StringCmp(AtomHere," CA ")) ||
				(!StringCmp(AtomHere," CB ")) ||
				(!StringCmp(AtomHere," HA ")) ||
				(!StringCmp(AtomHere,"1HA ")) ||
				(!StringCmp(AtomHere,"2HA ")) ||
				(!StringCmp(AtomHere," N  "))));
				else if ((!StringCmp(AtomName," O  ")) &&
				(/*(!StringCmp(AtomHere," CB ")) ||*/
				(!StringCmp(AtomHere," CA ")) ||
				(!StringCmp(AtomHere," N  ")) ||
				(!StringCmp(AtomHere," HA ")) ||
/*				(!StringCmp(AtomHere,"2HA ")) ||*/
/*				(!StringCmp(AtomHere," C  ")) ||*/
				 (!StringCmp(AtomHere,"1HA "))));
				else if (isProNew && (
				(!StringCmp(AtomHere," N  ")) ||
				(!StringCmp(AtomHere," N  "))
)) /*printf("excused\n")*/;
				else if ((!StringCmp(AtomHere," CA ")) && (!StringCmp(AtomName," CA ")));
				else doit=1;
			}
			else doit=1;
			if (doit) {
				/* get atom VDW radii for paldHere and AtomName */
				switch (AtomHere[1]) {
					case 'C': 
						if (at2hn==30)
							radHere=VDWRAD_AROCH0;
						else if (at2hn==31)
							radHere=VDWRAD_AROCH1;
						else
							radHere=VDWRAD_C;
						break;
					case 'N': 
						radHere=VDWRAD_N;
						break;
					case 'O': 
						radHere=VDWRAD_O;
						break;
					case 'S': 
						radHere=VDWRAD_S;
						break;
					case 'P':
						radHere=VDWRAD_P;
						break;
					case 'E': 
						if (AtomHere[0]=='S')
							radHere=VDWRAD_SE;
						else
							radHere=VDWRAD_C;
						break;
					case 'H': 
						radHere=VDWRAD_H;
						break;
					default: 
						radHere=VDWRAD_C;
				}
				if ((StringRChr(smallatoms,(int)(AtomHere[2]))!=NULL) && !(AtomHere[2]=='X' && AtomHere[3]=='T') && ((resCrash=='E') || (resCrash=='Q') || (resCrash=='K') || (resCrash=='R'))) {
					if (radHere>VDWRAD_X)
						radHere=VDWRAD_X;
				}
				if ((AtomHere[2]=='D') && (resCrash=='P')) {
					if (radHere>VDWRAD_X)
						radHere=VDWRAD_X;
				}
				vHere[0]=(paldHere->pflvData)[0];
				vHere[1]=(paldHere->pflvData)[1];
				vHere[2]=(paldHere->pflvData)[2];
				VecSub(vDist,vHere,vCoord);
/*if ((!StringCmp(AtomHere," CA ")) && (!StringCmp(AtomName," CA ")) &&
(getMag(vDist)<BL_CACA)) {*/
				bumpdist=radHere+radNewAtom;
				tooclose=getMag(vDist)-bumpdist;
				/* soft sphere collision, linear
				   probability of allowed collision down
				   to relevant ATOM_BOUNCINESS A distance */
				if ((StringRChr(sccollide,(int)(AtomName[2]))!=NULL) || (StringRChr(sccollide,(int)(AtomHere[2]))!=NULL))
					bounciness=ATOM_BOUNCINESS_SC;
				else
					bounciness=ATOM_BOUNCINESS_BB;
				if (AtomName[1]=='H' || AtomHere[1]=='H') {
					if (!BUMPCHECK_HYDROGEN)
						bounciness=999999.99;
				}
				if (BUILD_FRAGS_ONLY && !ressrcres[(pmgdHere->pdnmgLink)->choice-1])
						bounciness=999999.99;
				/* just in case atoms got added that are not from frags, such as 1st three CAs */
				if (BUILD_FRAGS_ONLY && !ressrcres[(pmgdCrash->pdnmgLink)->choice-1])
						bounciness=999999.99;
				if (bounciness>999.9 && tooclose<0.0) {
					if ((!StringCmp(AtomHere," CA ")) && (!StringCmp(AtomName," CA "))) {
						tooclose=getMag(vDist)-BL_CACA;
						bounciness=0.25;
					}
					else
						tooclose=0.1;
				}
				if (tooclose<0) {
					/* atoms too close, unless... check for H-bonds */
					if ((((((AtomName[1]=='O') || (AtomName[1]=='N')) && ((AtomHere[1]=='O') || (AtomHere[1]=='N'))))) && (FindHBond(pmadHere,(PMAD)(paldHere->pfbParent),FALSE,TRUE)==TRUE)) {
						/* possible H-bond so get smaller radii */
						if (AtomName[1]=='N') {
							if (at1hn==31)
								radNewAtom=VDWRAD_N3H1;
							else if (at1hn==32)
								radNewAtom=VDWRAD_N3H2;
							else /* if (at1hn==43) */
								radNewAtom=VDWRAD_N4H3;
						}
						if (AtomName[1]=='O') {
							if (at1hn==10)
								radNewAtom=VDWRAD_O1H0;
							else /* if (at1hn==21) */
								radNewAtom=VDWRAD_O2H1;
						}
						if (AtomHere[1]=='N') {
							if (at2hn==31)
								radHere=VDWRAD_N3H1;
							else if (at2hn==32)
								radHere=VDWRAD_N3H2;
							else /* if (at2hn==43) */
								radHere=VDWRAD_N4H3;
						}
						if (AtomHere[1]=='O') {
							if (at2hn==10)
								radHere=VDWRAD_O1H0;
							else /* if (at2hn==21) */
								radHere=VDWRAD_O2H1;
						}														
						tooclose=getMag(vDist)-(radHere+radNewAtom);
						if (tooclose<0) {
							/* still too close even if an H-bond so use bounciness to check acceptance */
							/* use non-H-bond van der Waals for bump checking */
							randbump=bumpdist*bounciness*fabs(Rand1());
							/* H's never cause a crash now */
							if ((randbump<fabs(tooclose)) && (AtomName[1]!='H') && (AtomHere[1]!='H')) {
								/* crash */
								crashcnt++;
					      ValNodeFree(vnpAtomList);
								return ERR_CRASH;
							}
						}
						/* note special case for N since may not have its H yet */
						else if ((at1hn%10!=0) || (!isProCrash && !StringCmp(AtomHere," N  "))) {
								HBondFound=FALSE;
								if (FindHBond(pmadHere,(PMAD)(paldHere->pfbParent),TRUE,TRUE)==FALSE) {
									/* there IS a hydrogen(s) attached to this atom so it could still hydrogen-bond */
									/* so mark atom for later checking */
									(pmadHere->bReserved)&=0xf1;
									if (at1hn%10>0) {
										(pmadHere->bReserved)|=0x02;
										HBondFound=TRUE;
									}
									if (at1hn%10>1)
										(pmadHere->bReserved)|=0x04;
									if (at1hn%10>2)
										(pmadHere->bReserved)|=0x08;
									/* record thing to H-bond with */
									/* if H on N not placed yet, add new node to Hbond list too */
									balreadyadded=FALSE;
									if (!StringCmp(AtomHere," N  ")) {
										pmadTemp=FindAtomName(pmgdCrash->pvnmaAHead," H  ");
										if (pmadTemp!=NULL) {
											if (!IsAtomPlaced(pmadTemp,Model)) {												
												(((PMAD)(paldHere->pfbParent))->bReserved)&=0xf1;
												(((PMAD)(paldHere->pfbParent))->bReserved)|=0x02;
												/* donor, followed by acceptor */
												AddHBondToCheck((PMAD)(paldHere->pfbParent),pmadHere,!HBondFound);
												balreadyadded=TRUE;
												if (HBondFound) {
													/* donor, followed by acceptor */
													AddHBondToCheck(pmadHere,(PMAD)(paldHere->pfbParent),FALSE);
												}
												HBondFound=TRUE;
											}
										}
									}
									if (HBondFound && !balreadyadded) {
											AddHBondToCheck(pmadHere,(PMAD)(paldHere->pfbParent),TRUE);
									}
/*{
PHBS phbs;
PMAD pmad;
PMGD pmgd;
ErrPostEx(SEV_ERROR,1,1,"Bond list");
phbs=phbsHBondsToCheck;
while (phbs!=NULL) {
pmad=phbs->pmadDonor;
pmgd=(PMGD)(pmad->pfbParent);
ErrPostEx(SEV_ERROR,1,pmgd->pdnmgLink->choice,"Donor: %s",pmad->pcAName);
pmad=phbs->pmadAcceptor;
pmgd=(PMGD)(pmad->pfbParent);
ErrPostEx(SEV_ERROR,1,pmgd->pdnmgLink->choice,"Acceptor: %s",pmad->pcAName);
phbs=phbs->next;
}
}	*/
								/* there IS a hydrogen(s) attached to this atom so it could still hydrogen-bond */
								/* so mark atom for later checking */
/*								if (at2hn%10>0)
									(((PMAD)(paldHere->pfbParent))->bReserved)|=0x02;
								if (at2hn%10>1)
									(((PMAD)(paldHere->pfbParent))->bReserved)|=0x04;
								if (at2hn%10>2)
									(((PMAD)(paldHere->pfbParent))->bReserved)|=0x08;
	*/							/* record thing to H-bond with */
		/*						ValNodeAddPointer(&(vnpHBondsToCheck),0,pmadHere);
			*/
									if (!HBondFound) {
										crashcnt++;
						     		ValNodeFree(vnpAtomList);
										return ERR_CRASH;
									}
								}
						}
						else {
							/* all hydrogens are placed so check for H-bond */
							if (FindHBond(pmadHere,(PMAD)(paldHere->pfbParent),TRUE,TRUE)==FALSE) {
								/* no H-bond so check normally */
								randbump=bumpdist*bounciness*fabs(Rand1());
								if ((randbump<fabs(tooclose)) && (AtomName[1]!='H') && (AtomHere[1]!='H')) {
									/* crash if Monte Carlo soft collision
									   test fails */
									crashcnt++;
								  ValNodeFree(vnpAtomList);
									return ERR_CRASH;
								}
							}
						}
					}
					else {
						/* no H-bond possible */
						randbump=bumpdist*bounciness*fabs(Rand1());
						if ((randbump<fabs(tooclose)) && (AtomName[1]!='H') && (AtomHere[1]!='H')) {
								/* crash if Monte Carlo soft collision
								   test fails */
/*printf("hard Crash %s %d %s %d
%f<%f\n",AtomHere,(pmgdCrash->pdnmgLink)->choice,AtomName,(pmgdHere->pdnmgLink)->choice,randbump,tooclose);*/
							crashcnt++;
						  ValNodeFree(vnpAtomList);
							return ERR_CRASH;
						}
/*else printf("soft Crash %s %d %s %d
%f>%f\n",AtomHere,(pmgdCrash->pdnmgLink)->choice,AtomName,(pmgdHere->pdnmgLink)->choice,randbump,tooclose);*/
					}
				}
			}
		}
    vnpHere=vnpHere->next;
  }
  /* free the atom list created */
  ValNodeFree(vnpAtomList);
	/* mark atom as having valid co-ordinates */
	(pmadHere->bReserved)|=0x01;
	/* only insert node if haven't already */
/*	if (AtomName[1]!='H') */
	/* if building only frags, only add to BD tree if from a frag */
/*	if (!BUILD_FRAGS_ONLY || ressrcres[(pmgdHere->pdnmgLink)->choice-1]) {*/
	AddToBDTree((PFB)pmadHere,Model,&(pwsThis->pbdTreeHead));
/*	}*/
	return ERR_SUCCESS;
}

/* given U basis vectors and CB relative to them, and vDir (either R+ or R- usually),
   calculates atom position of N, for example, returning it in vDest */
TrajErr FindAdjacentPeptideAtom(vec uone,vec utwo,vec uthree,vec vDir,vec vCBHereRef,FloatLo ba,FloatLo bl,FloatLo cosba,vec vDest)
{
	vec vTmp,vTest,vRef,vTrue;
	FloatLo testhi,testlo,errhere=0.0,errhi,errlo,thedot;
	
	vRef[0]=Dot(vDir,uone);
	vRef[1]=Dot(vDir,utwo);
	vRef[2]=Dot(vDir,uthree);
	VecAdd(vTmp,vCBHereRef,vRef);
	testlo=-180.0;
	testhi=0.0;
	CalcNextCoOrd(vTmp,vCBHereRef,vZero,vTest,1.0,ba,testlo);
	errlo=fabs(Dot(vTest,vRef)-cosba);
	CalcNextCoOrd(vTmp,vCBHereRef,vZero,vTest,1.0,ba,testhi);
	errhi=fabs(Dot(vTest,vRef)-cosba);
	while ((testhi-testlo)>BACKBONE_PRECISION) {
		CalcNextCoOrd(vTmp,vCBHereRef,vZero,vTest,1.0,ba,(testlo+testhi)/2.0);
		errhere=fabs(Dot(vTest,vRef)-cosba);
		if (errhi<errlo) {
			errlo=errhere;
 			testlo=(testlo+testhi)/2.0;
		}
		else {
			errhi=errhere;
			testhi=(testlo+testhi)/2.0;
		}
	}
	/* approx. 0.5 degrees error when 15-20 degrees */
	thedot=Dot(vTest,vRef);
	if (thedot>1.0) thedot=1.0;
	if (thedot<-1.0) thedot=-1.0;
	errhere=RADTODEG*fabs(acos(thedot)-acos(cosba));
	/* this function is called for residue 2 only! */
	if (BUILD_FRAGS_ONLY && !ressrcres[2-1])
		errhere=0;
	if (errhere*errhere>BACKBONE_ERROR_TOLERANCE/5.0) {
                ErrPostEx(SEV_INFO,1,15,"Unable to place atom, backing up, error is %f",errhere);
		fpabad++;
		return ERR_BADBACK;
	}
	CalcNextCoOrd(vTmp,vCBHereRef,vZero,vTest,1.0,ba,(testlo+testhi)/2.0);
	VecScale(vTrue,uone,vTest[0]);
	VecScale(vTmp,utwo,vTest[1]);
	VecAdd(vTrue,vTrue,vTmp);
	VecScale(vTmp,uthree,vTest[2]);
	VecAdd(vTrue,vTrue,vTmp);
	VecScale(vDest,vTrue,bl);
	return ERR_SUCCESS;
}

/* given U basis vectors and two adjacent CAs and residue type, returns CB direction
   in vCBRef (relative to U basis vectors) and vCBHere (relative to molecular
   co-ordinate system) */
void FindCBDir(vec uone,vec utwo,vec uthree,vec vCANext,Int2 res,Int2 resNext,vec vCBRef,FloatLo CBOffset)
{
	vec vDist,vTmp,vNoise;
	FloatLo dist;  /*,noisemag,noisephi,noisetheta;*/
	Int2 distbin,cnt;
	FloatLo dbin[NUMBINS]={5.1,5.6,6.1,6.6,7.0};

        /* dist==d in Table I of Rey and Skolnick, J. Comp. Chem.,  */
        VecSub(vDist,vCANext,vCALast);  
        dist=getMag(vDist);
	/* distance bins range from 0 to NUMBINS */
	distbin=NUMBINS;
	for (cnt=0;cnt<NUMBINS;cnt++) {
		if (dist<dbin[cnt]) {
			distbin=cnt;
			break;
		}
	}
        /* use look-up table to get cosd1,cosd2,cosd3 */
	/* this strange mapping accounts for the error in Table I of the
	   Rey and Skolnick paper */
	if (resNext!=CISPRO) {
		vCBRef[0]=-cbdir[res][distbin][2];
		vCBRef[1]=cbdir[res][distbin][0];
		vCBRef[2]=-cbdir[res][distbin][1];
	}
	else {
		vCBRef[0]=cbcisdir[res][0];
		vCBRef[1]=cbcisdir[res][1];
		vCBRef[2]=cbcisdir[res][2];
	}
/* add some random noise */
vNoise[0]=0.05*Rand1Distrib();
vNoise[1]=0.05*Rand1Distrib();
vNoise[2]=0.05*Rand1Distrib()+CBOffset;
VecAdd(vCBRef,vCBRef,vNoise);
Normalize(vCBRef,vCBRef);
        /* then true CiBeta dir - U*(vCBRef) where U=[uone utwo uthree] */
	VecScale(vCBHere,uone,vCBRef[0]);
	VecScale(vTmp,utwo,vCBRef[1]);
	VecAdd(vCBHere,vCBHere,vTmp);
	VecScale(vTmp,uthree,vCBRef[2]);
	VecAdd(vCBHere,vCBHere,vTmp);
}

/* given three consecutive CAs, calculates and returns R+, R- and U basis vectors */
void SetUpRefAxes(vec uone,vec utwo,vec uthree,vec vRMinus,vec vRPlus,vec vCANext)
{
        VecSub(vRMinus,vCALast,vCAHere);
        VecSub(vRPlus,vCANext,vCAHere);
        /* to ensure vectors are same length when we later add them,
           normalize both vectors */
        Normalize(vRMinus,vRMinus);
        Normalize(vRPlus,vRPlus);
        Cross(uone,vRPlus,vRMinus);
        Normalize(uone,uone);
        VecAdd(utwo,vRPlus,vRMinus);
        Normalize(utwo,utwo);
        NegateVec(utwo,utwo);
        Cross(uthree,uone,utwo);
        Normalize(uthree,uthree);
}

/* Ref is usually R+ or R-, and given CA one atom name, and eta or zeta, finds the
   other peptide atom or carbonyl O (depending on what ba and bl are given) saving
   it in vDest */
void FindAtomInResMinusOne(CharPtr AtomName,PMGD pmgdHere,vec vRef,FloatLo ba,FloatLo bl,Int2 Model,vec vDest)
{
	vec vTmp,vTrue,vHere;

	GetCoOrds(pmgdHere,AtomName,vCAHere,vTrue,Model);
	Cross(vTmp,vRef,vTrue);
	VecAdd(vHere,vTmp,vRef);
	CalcNextCoOrd(vHere,vTmp,vZero,vDest,bl,90,-ba);
	NegateVec(vDest,vDest);
}

/* given CAi and CAi-1, Ni and Ci-1, finds Hi and puts it in vDest */
void FindHi(vec vDest)
{
	vec vNHereAbs,vTmp,vTmp2,vCLastAbs;

	VecAdd(vNHereAbs,vNHere,vCAHere);
	VecAdd(vCLastAbs,vCLast,vCALast);
	VecSub(vTmp,vCLastAbs,vNHereAbs);
	VecSub(vTmp2,vCAHere,vNHereAbs);
	Normalize(vTmp2,vTmp2);
	Normalize(vTmp,vTmp);
	VecAdd(vTmp,vTmp,vTmp2);
	NegateVec(vTmp,vTmp);
	Normalize(vTmp,vTmp);
	VecScale(vTmp,vTmp,BL_NH);
	VecAdd(vDest,vTmp,vNHereAbs);
}

FloatLo GetProError(vec vCB,vec vCG,vec vCD,vec vN,vec vCLast,vec vA3)
{
	FloatLo tmp,proerror,dih,thedot;
	vec vA1,vA2,vA4;

	VecSub(vA1,vN,vCD);
	/* s.d. of NCD in A is about 100 times smaller than that for angles in degrees */
	tmp=(getMag(vA1)-BL_P_NCD)*100.0;  /* 0.014 sd */
	proerror=tmp*tmp;
	VecSub(vA2,vCG,vCD);               /* 1.5 sd */
	thedot=Dot(vA1,vA2)/(getMag(vA1)*getMag(vA2));
	if (thedot>1.00) thedot=1.00;
    if (thedot<-1.00) thedot=-1.00;
	tmp=RADTODEG*acos(thedot)-BA_P_CGCDN;
	proerror+=tmp*tmp;
	NegateVec(vA1,vA1);
	/* vA3 defined above as CA-CB */
	VecSub(vA4,vCB,vN);
	/* get Ca - N in vA4 */
	VecAdd(vA4,vA4,vA3);               /* 1.4 sd */
	thedot=Dot(vA1,vA4)/(getMag(vA1)*getMag(vA4));
	if (thedot>1.00) thedot=1.00;
    if (thedot<-1.00) thedot=-1.00;
	tmp=RADTODEG*acos(thedot)-BA_P_CANCD;
	proerror+=tmp*tmp;                 /* about 7.5 s.d. */
	/* try to make C-N-CA-CD planar i.e. N is planar */
	/* except at cis Pro residues */
	/* vCALast is equal to vCA in this context */
	if (((!cisLast) && (vCLast[0]!=0.0)) || (vCLast[1]!=0.0) || (vCLast[2]!=0.0)) {
		GetDihedral(vCLast,vN,vCALast,vCD,0,&dih,&tmp,&tmp,&tmp,&tmp,&tmp);
		tmp=(dih-180.0)/5.0;
		proerror+=tmp*tmp;
	}
/*printf("proerror: %f\n",proerror);*/
	return proerror;
}

TrajErr PlaceRotAtom(CharPtr AtomName,PMGD pmgdThis,Int2 ref,FloatLo bondlength,FloatLo bondangle,FloatLo chi)
{
  TrajErr err;

	/* always use model 1 for now */
        CalcNextCoOrd(vRef[ref],vRef[ref+1],vRef[ref+2],vRef[ref+3],bondlength,bondangle,chi);
	if ((err=AssignCoords(AtomName,pmgdThis,vRef[ref+3],1))!=ERR_SUCCESS)
		return err;
	return ERR_SUCCESS;
}

TrajErr BuildSC(PMGD pmgdThis,vec vCBLast,vec vNLast,FloatLo chi1,FloatLo chi2,FloatLo chi3,FloatLo chi4)
{
	/*vec vNew,vNew2,vTmp1,vTmp2;*/
	/*vec vHAdir;*/
	/*vec vCBdir={0,0,1};*/
	TrajErr err;
	/*FloatLo cosb,sinb;*/
	/*FloatLo blnca,blcac,bancac;*/
	Char ResId;

	/* maintain X here if present */
	ResId=pmgdThis->pcIUPAC[0];
	/* put normal sidechain on N- or C-terminus modified residues */
	if (IsNTermModified(pmgdThis) || IsCTermModified(pmgdThis))
		ResId=GetAAFromIDict(pmgdThis);
	VecScale(vRef[1],vCALast,1.0);
	VecScale(vRef[2],vCBLast,1.0);
	VecScale(vRef[0],vNLast,1.0);
/*        if (ResID!='G') {
		VecScale(vNew,vCBdir,bl_cacb[ResIUPAC]);
	        for (cnt=0;cnt<3;cnt++)
			vNew2[cnt]=vNew[0]*vXDir[cnt]+vNew[1]*vYDir[cnt]+vNew[2]*vZDir[cnt];
		VecAdd(vRef[2],vNew2,vRef[1]);
		if (!AssignCoords(" CB ",pmgdThis,vRef[2],1))
			return 0;
	}
	else {
		VecScale(vNew,vCBdir,BL_G_CA2HA);
	        for (cnt=0;cnt<3;cnt++)
			vNew2[cnt]=vNew[0]*vXDir[cnt]+vNew[1]*vYDir[cnt]+vNew[2]*vZDir[cnt];
	        VecAdd(vRef[2],vNew2,vCALast);
	       	if (!AssignCoords("2HA ",pmgdThis,vRef[2],1))
			return 0;
	}
*/	/* note: all Hydrogen dihedrals and bond lengths, angles, were
	   chosen somewhat arbitrarily - H.F. */
	/* it is crucual that atoms with bonded hydrogens be placed prior
	   to actual hydrogen atoms which are attached, to allow proper
	   H-bond testing */
	switch (ResId) {
		case 'A':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HB ",pmgdThis,0,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			break;
		case 'C':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" SG ",pmgdThis,0,BL_C_CBSG,BA_C_CACBSG,chi1))!=ERR_SUCCESS) return err;
			/* C-S-H angles only approximate */
			/* only put HG on if not a disulfide bridge! */
			if (InDisulphide((pmgdThis->pdnmgLink)->choice)==0.0)
				if ((err=PlaceRotAtom(" HG ",pmgdThis,1,BL_SH,108,100))!=ERR_SUCCESS) return err;
			break;
		case 'D':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_D_CBCG,BA_D_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OD1",pmgdThis,1,BL_D_CGOD1,BA_D_CBCGOD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OD2",pmgdThis,1,BL_D_CGOD2,BA_D_CBCGOD2,chi2+180))!=ERR_SUCCESS) return err;
			break;
		case 'E':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_E_CBCG,BA_E_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_E_CGCD,BA_E_CBCGCD,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_E_CDOE1,BA_E_CGCDOE1,chi3))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OE2",pmgdThis,2,BL_E_CDOE2,BA_E_CGCDOE2,chi3+180))!=ERR_SUCCESS) return err;
			break;
		case 'F':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_F_CBCG,BA_F_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_F_CGCD1,BA_F_CBCGCD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HD1",pmgdThis,2,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_F_CD1CE1,BA_F_CGCD1CE1,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_F_CGCD2,BA_F_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE2",pmgdThis,2,BL_F_CD2CE2,BA_F_CGCD2CE2,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_F_CE2CZ,BA_F_CD2CE2CZ,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HZ ",pmgdThis,4,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'H':
			/* Histidine epsilon-H ("HIE") */
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_H_CBCG,BA_H_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" ND1",pmgdThis,1,BL_H_CGND1,BA_H_CBCGND1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_H_ND1CE1,BA_H_CGND1CE1,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_CH,BA_HISCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_H_CGCD2,BA_H_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_HISCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_H_CD2NE2,BA_H_CGCD2NE2,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_NH,BA_HISCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'I':
			if ((err=PlaceRotAtom(" HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG1",pmgdThis,0,BL_I_CBCG1,BA_I_CACBCG1,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG1",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG1",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_I_CG1CD1,BA_I_CBCG1CD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HD1",pmgdThis,2,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HD1",pmgdThis,2,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HD1",pmgdThis,2,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			/* +/- 120 determines which enantionmer of Ile */
			if ((err=PlaceRotAtom(" CG2",pmgdThis,0,BL_I_CBCG2,BA_I_CACBCG2,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG2",pmgdThis,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG2",pmgdThis,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HG2",pmgdThis,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'K':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HZ ",pmgdThis,4,BL_NH3,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HZ ",pmgdThis,4,BL_NH3,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HZ ",pmgdThis,4,BL_NH3,BA_XCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'L':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_L_CBCG,BA_L_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_L_CGCD1,BA_L_CBCGCD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HD1",pmgdThis,2,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HD1",pmgdThis,2,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HD1",pmgdThis,2,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_L_CGCD2,BA_L_CBCGCD2,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HD2",pmgdThis,2,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HD2",pmgdThis,2,BL_CH,BA_XCH,120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HD2",pmgdThis,2,BL_CH,BA_XCH,-120))!=ERR_SUCCESS) return err;
			break;
		case 'M':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_M_CBCG,BA_M_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" SD ",pmgdThis,1,BL_M_CGSD,BA_M_CBCGSD,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_M_SDCE,BA_M_CGSDCE,chi3))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HE ",pmgdThis,3,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'N':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_N_CBCG,BA_N_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OD1",pmgdThis,1,BL_N_CGOD1,BA_N_CBCGOD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" ND2",pmgdThis,1,BL_N_CGND2,BA_N_CBCGND2,chi2+180))!=ERR_SUCCESS) return err;
			/* planar NH2 ? */
			if ((err=PlaceRotAtom("1HD2",pmgdThis,2,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HD2",pmgdThis,2,BL_NH,BA_XNH,180))!=ERR_SUCCESS) return err;
			break;
		case 'P':
			/* built elsewhere */
			break;
		case 'Q':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_Q_CBCG,BA_Q_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_Q_CGCD,BA_Q_CBCGCD,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_Q_CDOE1,BA_Q_CGCDOE1,chi3))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_Q_CDNE2,BA_Q_CGCDNE2,chi3+180))!=ERR_SUCCESS) return err;
			/* planar NH2? */
			if ((err=PlaceRotAtom("1HE2",pmgdThis,3,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HE2",pmgdThis,3,BL_NH,BA_XNH,180))!=ERR_SUCCESS) return err;
			break;
		case 'R':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_R_CBCG,BA_R_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_R_CGCD,BA_R_CBCGCD,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" NE ",pmgdThis,2,BL_R_CDNE,BA_R_CGCDNE,chi3))!=ERR_SUCCESS) return err;
			/* planar N */
			if ((err=PlaceRotAtom(" HE ",pmgdThis,3,BL_NH,BA_XNH,chi4+180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_R_NECZ,BA_R_CDNECZ,chi4))!=ERR_SUCCESS) return err;
			/* planar NH2 groups */
			if ((err=PlaceRotAtom(" NH2",pmgdThis,4,BL_R_CZNH2,BA_R_NECZNH2,CHI5_R))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_GUANNH,BA_XNH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_GUANNH,BA_XNH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" NH1",pmgdThis,4,BL_R_CZNH1,BA_R_NECZNH1,180+CHI5_R))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HH1",pmgdThis,5,BL_GUANNH,BA_XNH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HH1",pmgdThis,5,BL_GUANNH,BA_XNH,180))!=ERR_SUCCESS) return err;
			break;
		case 'S':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OG ",pmgdThis,0,BL_S_CBOG,BA_S_CACBOG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HG ",pmgdThis,1,BL_OH,BA_XOH,60))!=ERR_SUCCESS) return err;
			break;
		case 'T':
			if ((err=PlaceRotAtom(" HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OG1",pmgdThis,0,BL_T_CBOG1,BA_T_CACBOG1,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HG1",pmgdThis,1,BL_OH,BA_XOH,120))!=ERR_SUCCESS) return err;
			/* +/- 120 determines which enantionmer of Thr */
			if ((err=PlaceRotAtom(" CG2",pmgdThis,0,BL_T_CBCG2,BA_T_CACBCG2,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG2",pmgdThis,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG2",pmgdThis,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HG2",pmgdThis,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'V':
			if ((err=PlaceRotAtom(" HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG1",pmgdThis,0,BL_V_CBCG1,BA_V_CACBCG1,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG1",pmgdThis,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG1",pmgdThis,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HG1",pmgdThis,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG2",pmgdThis,0,BL_V_CBCG2,BA_V_CACBCG2,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("1HG2",pmgdThis,1,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HG2",pmgdThis,1,BL_CH,BA_XCH,-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("3HG2",pmgdThis,1,BL_CH,BA_XCH,120))!=ERR_SUCCESS) return err;
			break;
		case 'W':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_W_CBCG,BA_W_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_W_CGCD1,BA_W_CBCGCD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HD1",pmgdThis,2,BL_CH,BA_HISCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" NE1",pmgdThis,2,BL_W_CD1NE1,BA_W_CGCD1NE1,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_NH,BA_HISCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_W_CGCD2,BA_W_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE2",pmgdThis,2,BL_W_CD2CE2,BA_W_CGCD2CE2,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CZ2",pmgdThis,3,BL_W_CE2CZ2,BA_W_CD2CE2CZ2,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HZ2",pmgdThis,4,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CH2",pmgdThis,4,BL_W_CZ2CH2,BA_W_CE2CZ2CH2,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HH2",pmgdThis,5,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE3",pmgdThis,2,BL_W_CD2CE3,BA_W_CGCD2CE3,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE3",pmgdThis,3,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CZ3",pmgdThis,3,BL_W_CE3CZ3,BA_W_CD2CE3CZ3,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HZ3",pmgdThis,4,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			break;
		case 'X':
			/* modified amino acid */
			switch (pmgdThis->iIDict) {
                case 117:
                case 118:
                case 119:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        /* planar N */
                        if ((err=PlaceRotAtom(" HZ ",pmgdThis,4,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_NCA,BA_NCAC,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ1",pmgdThis,5,BL_X_CC,BA_X_XCC,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ1",pmgdThis,6,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ1",pmgdThis,6,BL_CH,BA_XCH,120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HQ1",pmgdThis,6,BL_CH,BA_XCH,-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OQ2",pmgdThis,5,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                  break;
                case 216:
                case 217:
                case 218:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" SG ",pmgdThis,0,BL_C_CBSG,BA_C_CACBSG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_M_SDCE,108,100))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,30))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,-90))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_X_CC,BA_CACO,150))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE ",pmgdThis,3,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH2",pmgdThis,4,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH2",pmgdThis,5,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH1",pmgdThis,4,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH1",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH1",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ ",pmgdThis,5,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ ",pmgdThis,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ ",pmgdThis,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI ",pmgdThis,6,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI ",pmgdThis,7,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL ",pmgdThis,7,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM2",pmgdThis,8,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM2",pmgdThis,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM2",pmgdThis,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HM2",pmgdThis,9,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM1",pmgdThis,8,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM1",pmgdThis,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM1",pmgdThis,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HX ",pmgdThis,10,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HX ",pmgdThis,10,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR ",pmgdThis,10,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HR ",pmgdThis,11,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT ",pmgdThis,11,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU2",pmgdThis,12,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HU2",pmgdThis,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HU2",pmgdThis,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HU2",pmgdThis,13,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU1",pmgdThis,12,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HU1",pmgdThis,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HU1",pmgdThis,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HU1",pmgdThis,13,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 219:
                case 220:
                case 221:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" SG ",pmgdThis,0,BL_C_CBSG,BA_C_CACBSG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_M_SDCE,108,100))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,30))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,-90))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_X_CC,BA_CACO,150))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE ",pmgdThis,3,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH2",pmgdThis,4,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH2",pmgdThis,5,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH1",pmgdThis,4,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH1",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH1",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ ",pmgdThis,5,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ ",pmgdThis,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ ",pmgdThis,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI ",pmgdThis,6,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI ",pmgdThis,7,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL ",pmgdThis,7,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM2",pmgdThis,8,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM2",pmgdThis,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM2",pmgdThis,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HM2",pmgdThis,9,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM1",pmgdThis,8,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM1",pmgdThis,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM1",pmgdThis,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HX ",pmgdThis,10,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HX ",pmgdThis,10,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR ",pmgdThis,10,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HR ",pmgdThis,11,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT ",pmgdThis,11,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU2",pmgdThis,12,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HU2",pmgdThis,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HU2",pmgdThis,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HU2",pmgdThis,13,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU1",pmgdThis,12,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HU1",pmgdThis,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HU1",pmgdThis,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CF ",pmgdThis,13,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HF ",pmgdThis,14,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HF ",pmgdThis,14,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CY ",pmgdThis,14,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HY ",pmgdThis,15,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW ",pmgdThis,15,BL_X_CCDBL,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CJ2",pmgdThis,16,BL_X_CC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HJ2",pmgdThis,17,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HJ2",pmgdThis,17,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HJ2",pmgdThis,17,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CJ1",pmgdThis,16,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HJ1",pmgdThis,17,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HJ1",pmgdThis,17,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HJ1",pmgdThis,17,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 223:
                case 224:
                case 225:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" SG ",pmgdThis,0,BL_C_CBSG,BA_C_CACBSG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_M_SDCE,108,100))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE2",pmgdThis,2,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE1",pmgdThis,3,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE1",pmgdThis,3,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HZ ",pmgdThis,4,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HZ ",pmgdThis,4,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH ",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH ",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ ",pmgdThis,5,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ ",pmgdThis,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ ",pmgdThis,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI ",pmgdThis,6,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HI ",pmgdThis,7,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HI ",pmgdThis,7,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL ",pmgdThis,7,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL ",pmgdThis,8,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL ",pmgdThis,8,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM ",pmgdThis,8,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM ",pmgdThis,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM ",pmgdThis,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HX ",pmgdThis,10,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HX ",pmgdThis,10,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR ",pmgdThis,10,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HR ",pmgdThis,11,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HR ",pmgdThis,11,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT ",pmgdThis,11,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HT ",pmgdThis,12,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HT ",pmgdThis,12,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU ",pmgdThis,12,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HU ",pmgdThis,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HU ",pmgdThis,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CF ",pmgdThis,13,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HF ",pmgdThis,14,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HF ",pmgdThis,14,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CY ",pmgdThis,14,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HY ",pmgdThis,15,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HY ",pmgdThis,15,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW ",pmgdThis,15,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW ",pmgdThis,16,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW ",pmgdThis,16,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CJ ",pmgdThis,16,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HJ ",pmgdThis,17,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HJ ",pmgdThis,17,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HJ ",pmgdThis,17,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 230:
                case 231:
                case 232:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        /* planar N */
                        if ((err=PlaceRotAtom(" HZ ",pmgdThis,4,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_NCA,BA_NCAC,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OQ1",pmgdThis,5,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ2",pmgdThis,5,BL_X_CC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ2",pmgdThis,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ2",pmgdThis,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI ",pmgdThis,6,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HI ",pmgdThis,7,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HI ",pmgdThis,7,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL ",pmgdThis,7,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL ",pmgdThis,8,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL ",pmgdThis,8,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM ",pmgdThis,8,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM ",pmgdThis,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM ",pmgdThis,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HX ",pmgdThis,10,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HX ",pmgdThis,10,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR ",pmgdThis,10,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HR ",pmgdThis,11,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HR ",pmgdThis,11,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT ",pmgdThis,11,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HT ",pmgdThis,12,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HT ",pmgdThis,12,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU ",pmgdThis,12,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HU ",pmgdThis,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HU ",pmgdThis,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CF ",pmgdThis,13,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HF ",pmgdThis,14,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HF ",pmgdThis,14,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CY ",pmgdThis,14,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HY ",pmgdThis,15,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HY ",pmgdThis,15,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW ",pmgdThis,15,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW ",pmgdThis,16,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW ",pmgdThis,16,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CJ ",pmgdThis,16,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HJ ",pmgdThis,17,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HJ ",pmgdThis,17,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CK ",pmgdThis,17,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HK ",pmgdThis,18,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HK ",pmgdThis,18,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HK ",pmgdThis,18,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 245:
                case 246:
                case 247:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_R_CBCG,BA_R_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_R_CGCD,BA_R_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE ",pmgdThis,2,BL_NH3CA,BA_NCAC,chi3))!=ERR_SUCCESS) return err;
                        /* planar N */
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_X_NC,BA_X_CNC,-120))!=ERR_SUCCESS) return err;
                        /* planar NH2 groups */
                        if ((err=PlaceRotAtom(" NH2",pmgdThis,4,BL_NCA,BA_HISCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_NH3,BA_XNH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_NH3,BA_XNH,120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NH1",pmgdThis,4,BL_NCA,BA_HISCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ1",pmgdThis,5,BL_NCA,BA_X_CNC,-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ1",pmgdThis,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ1",pmgdThis,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HQ1",pmgdThis,6,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ2",pmgdThis,5,BL_NCA,BA_X_CNC,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ2",pmgdThis,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ2",pmgdThis,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HQ2",pmgdThis,6,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 252:
                case 253:
                case 254:
                        if ((err=PlaceRotAtom(" HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG2",pmgdThis,0,BL_D_CBCG,BA_D_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OD2",pmgdThis,1,BL_D_CGOD1,BA_D_CBCGOD1,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OD3",pmgdThis,1,BL_D_CGOD2,BA_D_CBCGOD2,chi2+180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD3",pmgdThis,2,BL_OH,BA_XOH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" SG1",pmgdThis,0,BL_C_CBSG,BA_C_CACBSG,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_M_SDCE,BA_M_CGSDCE,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD1",pmgdThis,2,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD1",pmgdThis,2,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HD1",pmgdThis,2,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                  break;
                case 257:
                case 258:
                case 259:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_Q_CBCG,BA_Q_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_Q_CGCD,BA_Q_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_Q_CDOE1,BA_Q_CGCDOE1,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_Q_CDNE2,BA_Q_CGCDNE2,chi3+180))!=ERR_SUCCESS) return err;
                        /* planar NH2? */
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_NCA,BA_CACOXT,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HZ ",pmgdThis,4,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HZ ",pmgdThis,4,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HZ ",pmgdThis,4,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                  break;
                case 261:
                case 262:
                case 263:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_E_CBCG,BA_E_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_E_CGCD,BA_E_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_E_CDOE1,BA_E_CGCDOE1,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_COXT,BA_XOH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HZ ",pmgdThis,4,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HZ ",pmgdThis,4,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HZ ",pmgdThis,4,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE2",pmgdThis,2,BL_E_CDOE2,BA_E_CGCDOE2,chi3+180))!=ERR_SUCCESS) return err;
                  break;
                case 266:
                case 267:
                case 268:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_H_CBCG,BA_H_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" ND1",pmgdThis,1,BL_H_CGND1,BA_H_CBCGND1,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE3",pmgdThis,2,BL_NC,BA_NCAC,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE3",pmgdThis,3,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE3",pmgdThis,3,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HE3",pmgdThis,3,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_H_ND1CE1,BA_H_CGND1CE1,180)) !=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_CH,BA_HISCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_H_CGCD2,BA_H_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_HISCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_H_CD2NE2,BA_H_CGCD2NE2,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_NH,BA_HISCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 272:
                case 273:
                case 274:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ ",pmgdThis,4,BL_NH,BA_XNH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_NC,BA_NCAC,-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH ",pmgdThis,5,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH ",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH ",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                  break;
                case 275:
                case 276:
                case 277:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH1",pmgdThis,4,BL_NCA,BA_W_CGCD2CE2,-175.5))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH1",pmgdThis,5,BL_CH,BA_XCH,176.4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH1",pmgdThis,5,BL_CH,BA_XCH,58.5))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH1",pmgdThis,5,BL_CH,BA_XCH,-67.7))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH2",pmgdThis,4,BL_NCA,BA_W_CGCD2CE2,60.9))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_CH,BA_XCH,177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_CH,BA_XCH,56))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH2",pmgdThis,5,BL_CH,BA_XCH,-66))!=ERR_SUCCESS) return err;
                  break;
                case 278:
                case 279:
                case 280:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH1",pmgdThis,4,BL_CAC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH1",pmgdThis,5,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH1",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH1",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH2",pmgdThis,4,BL_CAC,BA_XCH,61))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_CH,BA_XCH,174))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_CH,BA_XCH,54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH2",pmgdThis,5,BL_CH,BA_XCH,-66))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH3",pmgdThis,4,BL_CAC,BA_XCH,-61))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH3",pmgdThis,5,BL_CH,BA_XCH,-174))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH3",pmgdThis,5,BL_CH,BA_XCH,66))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HH3",pmgdThis,5,BL_CH,BA_XCH,-54))!=ERR_SUCCESS) return err;
                  break;
                case 290:
                case 291:
                case 292:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_R_CBCG,BA_R_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_R_CGCD,BA_R_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE ",pmgdThis,2,BL_X_NC,BA_NCAC,chi3))!=ERR_SUCCESS) return err;
                        /* planar N */
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_X_NC,BA_X_CNC,180))!=ERR_SUCCESS) return err;
                        /* planar NH2 groups */
                        if ((err=PlaceRotAtom(" NH2",pmgdThis,4,BL_NC,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH2",pmgdThis,5,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH2",pmgdThis,5,BL_NH,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NH1",pmgdThis,4,BL_NC,BA_CACO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HH1",pmgdThis,5,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PQ ",pmgdThis,5,BL_X_NP,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI1",pmgdThis,6,BL_X_PODBL,BA_X_XPODBL,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI2",pmgdThis,6,BL_X_PO,BA_X_XPO,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI2",pmgdThis,7,BL_OH,BA_XOH,-30))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI3",pmgdThis,6,BL_X_PO,BA_X_XPO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI3",pmgdThis,7,BL_OH,BA_XOH,-60))!=ERR_SUCCESS) return err;
                  break;
                case 293:
                case 294:
                case 295:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_D_CBCG,BA_D_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OD1",pmgdThis,1,BL_COXT,BA_CACOXT,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PE ",pmgdThis,2,BL_X_PO,BA_XOH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OZ1",pmgdThis,3,BL_X_PODBL,BA_XOH,64))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OZ2",pmgdThis,3,BL_X_PO,BA_X_XPO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ2",pmgdThis,4,BL_OH,BA_X_POH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OZ3",pmgdThis,3,BL_X_PO,BA_X_XPO,-74))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ3",pmgdThis,4,BL_OH,BA_X_POH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OD2",pmgdThis,1,BL_D_CGOD2,BA_D_CBCGOD2,chi2+180))!=ERR_SUCCESS) return err;
                  break;
                case 296:
                case 297:
                case 298:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" SG ",pmgdThis,0,BL_C_CBSG,BA_C_CACBSG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PD ",pmgdThis,1,BL_X_PS,BA_XOH,119))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_X_PODBL,BA_X_XPODBL,-41))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE2",pmgdThis,2,BL_X_PO,BA_X_XPO,78))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_OH,BA_XOH,30))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE3",pmgdThis,2,BL_X_PO,BA_X_XPO,-161))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE3",pmgdThis,3,BL_OH,BA_XOH,-30))!=ERR_SUCCESS) return err;
                  break;
                case 299:
                case 300:
                case 301:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_H_CBCG,BA_H_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" ND1",pmgdThis,1,BL_H_CGND1,BA_H_CBCGND1,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_H_ND1CE1,BA_H_CGND1CE1,180)) !=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_CH,BA_HISCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_H_CGCD2,BA_H_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_HISCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_H_CD2NE2,BA_H_CGCD2NE2,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PZ ",pmgdThis,3,BL_X_NP,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OH1",pmgdThis,4,BL_X_PODBL,BA_X_XPODBL,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OH2",pmgdThis,4,BL_X_PO,BA_X_XPO,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HH2",pmgdThis,5,BL_OH,BA_XOH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OH3",pmgdThis,4,BL_X_PO,BA_X_XPO,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HH3",pmgdThis,5,BL_OH,BA_XOH,-60))!=ERR_SUCCESS) return err;
                  break;
                case 302:
                case 303:
                case 304:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_H_CBCG,BA_H_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" ND1",pmgdThis,1,BL_Y_CE2CZ,BA_H_CBCGND1,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PE1",pmgdThis,2,BL_X_NP,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OZ1",pmgdThis,3,BL_X_PODBL,BA_X_XPODBL,-80))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OZ2",pmgdThis,3,BL_X_PO,BA_X_XPO,20))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ2",pmgdThis,4,BL_OH,BA_XOH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OZ3",pmgdThis,3,BL_X_PO,BA_X_XPO,150))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ3",pmgdThis,4,BL_OH,BA_XOH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE3",pmgdThis,2,BL_X_ND1CE3,BA_H_CGND1CE1,180)) !=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE3",pmgdThis,3,BL_CH,BA_HISCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_H_CGCD2,BA_H_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_HISCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_H_CGND1,BA_H_CGCD2NE2,180))!=ERR_SUCCESS) return err;
                  break;
                case 305:
                case 306:
                case 307:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OG ",pmgdThis,0,BL_S_CBOG,BA_S_CACBOG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PD ",pmgdThis,1,BL_X_PO,BA_CACN,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_X_PODBL,BA_X_XPODBL,-173))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE2",pmgdThis,2,BL_X_PO,BA_X_XPO,-46))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_OH,BA_X_POH,-166))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE3",pmgdThis,2,BL_X_PO,BA_X_XPO,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE3",pmgdThis,3,BL_OH,BA_X_POH,-73))!=ERR_SUCCESS) return err;
                  break;
                case 308:
                case 309:
                case 310:
                        if ((err=PlaceRotAtom(" HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OG1",pmgdThis,0,BL_T_CBOG1,BA_T_CACBOG1,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PD ",pmgdThis,1,BL_X_PO,BA_CACN,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_X_PODBL,BA_X_XPODBL,-173))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE2",pmgdThis,2,BL_X_PO,BA_X_XPO,-46))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_OH,BA_X_POH,-166))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE3",pmgdThis,2,BL_X_PO,BA_X_XPO,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE3",pmgdThis,3,BL_OH,BA_X_POH,-73))!=ERR_SUCCESS) return err;
                        /* +/- 120 determines which enantionmer of Thr */
                        if ((err=PlaceRotAtom(" CG2",pmgdThis,0,BL_T_CBCG2,BA_T_CACBCG2,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG2",pmgdThis,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG2",pmgdThis,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HG2",pmgdThis,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 311:
                case 312:
                case 313:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_Y_CBCG,BA_Y_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_Y_CGCD1,BA_Y_CBCGCD1,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD1",pmgdThis,2,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_Y_CD1CE1,BA_Y_CGCD1CE1,180)) !=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_Y_CGCD2,BA_Y_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE2",pmgdThis,2,BL_Y_CD2CE2,BA_Y_CGCD2CE2,180)) !=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_Y_CE2CZ,BA_Y_CD2CE2CZ,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OH ",pmgdThis,4,BL_Y_CZOH,BA_Y_CE2CZOH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PQ ",pmgdThis,5,BL_X_PO,BA_CACOXT,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI1",pmgdThis,6,BL_X_PODBL,BA_XOH,64))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI2",pmgdThis,6,BL_X_PO,BA_X_XPO,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI2",pmgdThis,7,BL_OH,BA_X_POH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI3",pmgdThis,6,BL_X_PO,BA_X_XPO,-74))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI3",pmgdThis,7,BL_OH,BA_X_POH,180))!=ERR_SUCCESS) return err;
                  break;
                case 317:
                case 318:
                case 319:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("SEG ",pmgdThis,0,BL_X_CBSEG,BA_C_CACBSG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HG ",pmgdThis,1,BL_SH,108,100))!=ERR_SUCCESS) return err;
                  break;
                case 320:
                case 321:
                case 322:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_M_CBCG,BA_M_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("SED ",pmgdThis,1,BL_X_CBSEG,BA_M_CBCGSD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_M_SDCE,BA_M_CGSDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HE ",pmgdThis,3,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 323:
                case 324:
                case 325:
                        if ((err=PlaceRotAtom(" HB ",pmgdThis,0,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OG ",pmgdThis,0,BL_CO,BA_CACO,120))!=ERR_SUCCESS) return err;
                  break;
                case 330:
                case 331:
                case 332:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_Q_CBCG,BA_Q_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_Q_CGCD,BA_Q_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OE1",pmgdThis,2,BL_Q_CDOE1,BA_Q_CGCDOE1,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_Q_CDNE2,BA_Q_CGCDNE2,chi3+180))!=ERR_SUCCESS) return err;
                        /* planar NH2? */
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_NCA,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HZ ",pmgdThis,4,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HZ ",pmgdThis,4,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH ",pmgdThis,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH ",pmgdThis,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OQ ",pmgdThis,5,BL_S_CBOG,BA_XCH,179))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" PI ",pmgdThis,6,BL_X_PO,BA_CACOXT,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OL2",pmgdThis,7,BL_X_PODBL,BA_X_XCC,-1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OL3",pmgdThis,7,BL_X_PO,BA_X_POH,-125))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HL3",pmgdThis,8,BL_OH,BA_X_POH,166))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OL1",pmgdThis,7,BL_X_PO,BA_X_POH,123))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM ",pmgdThis,8,BL_S_CBOG,BA_CACOXT,49))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM ",pmgdThis,9,BL_CH,BA_XCH,-50))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM ",pmgdThis,9,BL_CH,BA_XCH,70))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_XCH,-170))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HX ",pmgdThis,10,BL_CH,BA_XCH,-175))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OR1",pmgdThis,10,BL_S_CBOG,BA_XOH,67))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HR1",pmgdThis,11,BL_OH,BA_X_POH,-177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR2",pmgdThis,10,BL_X_CC,BA_X_XCC,-52))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HR2",pmgdThis,11,BL_CH,BA_X_XCC,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HR2",pmgdThis,11,BL_CH,BA_X_XCC,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OT ",pmgdThis,11,BL_S_CBOG,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HT ",pmgdThis,12,BL_OH,BA_XOH,180))!=ERR_SUCCESS) return err;
                  break;
                case 333:
                case 334:
                case 335:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_H_CBCG,BA_H_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" ND1",pmgdThis,1,BL_H_CGND1,BA_H_CBCGND1,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_H_ND1CE1,BA_H_CGND1CE1,180)) !=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_X_CC,BA_HISCH,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HZ ",pmgdThis,4,BL_CH,BA_XCH,-179))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HZ ",pmgdThis,4,BL_CH,BA_XCH,61))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_X_CC,BA_XCH,-59))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH ",pmgdThis,5,BL_CH,BA_XCH,61))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH ",pmgdThis,5,BL_CH,BA_XCH,-59))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ ",pmgdThis,5,BL_X_CC,BA_XCH,-178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HQ ",pmgdThis,6,BL_CH,BA_XCH,172))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI2",pmgdThis,6,BL_X_CC,BA_XCH,-65))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OL4",pmgdThis,7,BL_CO,BA_CACO,3.5))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NL5",pmgdThis,7,BL_NC,BA_C_CACBSG,-179))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL5",pmgdThis,8,BL_NH,BA_XNH,-1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL5",pmgdThis,8,BL_NH,BA_XNH,-161))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NI1",pmgdThis,6,BL_NH3CA,BA_X_XPO,54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL1",pmgdThis,7,BL_NH3CA,BA_XCH,170))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL1",pmgdThis,8,BL_CH,BA_XCH,-177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL1",pmgdThis,8,BL_CH,BA_XCH,63))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HL1",pmgdThis,8,BL_CH,BA_XCH,-57))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL2",pmgdThis,7,BL_NH3CA,BA_XCH,-69))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL2",pmgdThis,8,BL_CH,BA_XCH,-178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL2",pmgdThis,8,BL_CH,BA_XCH,62))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HL2",pmgdThis,8,BL_CH,BA_XCH,-58))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL3",pmgdThis,7,BL_NH3CA,BA_XCH,50))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL3",pmgdThis,8,BL_CH,BA_XCH,174))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL3",pmgdThis,8,BL_CH,BA_XCH,54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HL3",pmgdThis,8,BL_CH,BA_XCH,-66))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_H_CGCD2,BA_H_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_HISCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NE2",pmgdThis,2,BL_H_CD2NE2,BA_H_CGCD2NE2,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_NH,BA_HISCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 336:
                case 337:
                case 338:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ ",pmgdThis,4,BL_GUANNH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_NC,BA_CACOXT,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OQ1",pmgdThis,5,BL_CO,BA_CACO,.5))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ2",pmgdThis,5,BL_X_CC,BA_C_CACBSG,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HQ2",pmgdThis,6,BL_CH,BA_XCH,-61))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HQ2",pmgdThis,6,BL_CH,BA_XCH,58))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI ",pmgdThis,6,BL_X_CC,BA_X_XCC,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HI ",pmgdThis,7,BL_CH,BA_XCH,-56))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HI ",pmgdThis,7,BL_CH,BA_XCH,63))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL ",pmgdThis,7,BL_X_CC,BA_XCH,-176))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL ",pmgdThis,8,BL_CH,BA_XCH,-65))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL ",pmgdThis,8,BL_CH,BA_XCH,54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM ",pmgdThis,8,BL_X_CC,BA_XCH,175))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM ",pmgdThis,9,BL_CH,BA_XCH,-76))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM ",pmgdThis,9,BL_CH,BA_XCH,43))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_X_XCC,163))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HX ",pmgdThis,10,BL_CH,BA_XCH,69))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" SR1",pmgdThis,10,BL_M_CGSD,BA_W_CGCD2CE2,-55))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT1",pmgdThis,11,BL_M_SDCE,BA_M_CGSDCE,-82))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HT1",pmgdThis,12,BL_CH,BA_H_CACBCG,88))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HT1",pmgdThis,12,BL_CH,BA_H_CACBCG,-150))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT3",pmgdThis,12,BL_X_CC,BA_H_CGCD2NE2,-31))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HT3",pmgdThis,13,BL_CH,BA_XCH,-103))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NU3",pmgdThis,13,BL_NCA,BA_NCAC,134))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HU3",pmgdThis,14,BL_GUANNH,BA_XNH,76))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR2",pmgdThis,10,BL_X_CC,BA_X_XCC,-169))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HR2",pmgdThis,11,BL_CH,BA_XCH,-151))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NT2",pmgdThis,11,BL_NCA,BA_NCAC,-28))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HT2",pmgdThis,12,BL_GUANNH,BA_XNH,-45))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU2",pmgdThis,12,BL_NC,BA_XNH,115))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OF ",pmgdThis,13,BL_CO,BA_CACO,177))!=ERR_SUCCESS) return err;
                  break;
                case 342:
                case 343:
                case 344:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HZ ",pmgdThis,4,BL_GUANNH,BA_XCH,-175))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_NCA,BA_W_CGCD2CE2,64))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HH ",pmgdThis,5,BL_CH,BA_XCH,57))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HH ",pmgdThis,5,BL_CH,BA_XCH,-62))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ ",pmgdThis,5,BL_X_CC,BA_X_XCC,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HQ ",pmgdThis,6,BL_CH,BA_XCH,59))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" OI1",pmgdThis,6,BL_COXT,BA_W_CGCD2CE2,174))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HI1",pmgdThis,7,BL_OH,BA_XOH,-170))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI2",pmgdThis,6,BL_X_CC,BA_X_XCC,-64))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HI2",pmgdThis,7,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HI2",pmgdThis,7,BL_CH,BA_XCH,-57))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL ",pmgdThis,7,BL_X_CC,BA_X_XCC,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL ",pmgdThis,8,BL_CH,BA_XCH,-58))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL ",pmgdThis,8,BL_CH,BA_XCH,62))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NM ",pmgdThis,8,BL_NCA,BA_NCAC,-178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HM ",pmgdThis,9,BL_GUANNH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HM ",pmgdThis,9,BL_GUANNH,BA_XCH,177))!=ERR_SUCCESS) return err;
                  break;
                case 345:
                case 346:
                case 347:
                        if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_K_CBCG,BA_K_CACBCG,chi1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HG ",pmgdThis,1,BL_CH,BA_XCH,chi2+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HG ",pmgdThis,1,BL_CH,BA_XCH,chi2-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CD ",pmgdThis,1,BL_K_CGCD,BA_K_CBCGCD,chi2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HD ",pmgdThis,2,BL_CH,BA_XCH,chi3+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HD ",pmgdThis,2,BL_CH,BA_XCH,chi3-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CE ",pmgdThis,2,BL_K_CDCE,BA_K_CGCDCE,chi3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HE ",pmgdThis,3,BL_CH,BA_XCH,chi4+120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HE ",pmgdThis,3,BL_CH,BA_XCH,chi4-120))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" NZ ",pmgdThis,3,BL_K_CENZ,BA_K_CDCENZ,chi4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CH ",pmgdThis,4,BL_X_NC,BA_X_CNC,81.4))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HH ",pmgdThis,5,BL_CH,BA_CACOXT,-1.7))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CQ ",pmgdThis,5,BL_NCA,BA_R_CDNECZ,177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HQ ",pmgdThis,6,BL_CH,BA_XNH,1.3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CI ",pmgdThis,6,BL_X_CCDBL,BA_XNH,-177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL2",pmgdThis,7,BL_X_CC,BA_CNCA,-0.5))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HL2",pmgdThis,8,BL_CH,BA_XCH,-115))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HL2",pmgdThis,8,BL_CH,BA_XCH,124))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HL2",pmgdThis,8,BL_CH,BA_XCH,5))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CL1",pmgdThis,7,BL_X_CC,BA_XNH,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HL1",pmgdThis,8,BL_CH,BA_XNH,-2))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CM ",pmgdThis,8,BL_X_CCDBL,BA_XNH,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HM ",pmgdThis,9,BL_CH,BA_XNH,-1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CX ",pmgdThis,9,BL_X_CC,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HX ",pmgdThis,10,BL_CH,BA_XNH,.8))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CR ",pmgdThis,10,BL_X_CCDBL,BA_XNH,179))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT2",pmgdThis,11,BL_X_CC,BA_CNCA,3))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HT2",pmgdThis,12,BL_CH,BA_XCH,43))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HT2",pmgdThis,12,BL_CH,BA_XCH,-79))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HT2",pmgdThis,12,BL_CH,BA_XCH,159))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CT1",pmgdThis,11,BL_X_CC,BA_XNH,-179))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HT1",pmgdThis,12,BL_CH,BA_XNH,-19))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CU ",pmgdThis,12,BL_X_CCDBL,BA_XNH,159))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" HU ",pmgdThis,13,BL_CH,BA_XNH,1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CF ",pmgdThis,13,BL_X_CC,BA_XNH,178))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CY1",pmgdThis,14,BL_X_CC,BA_CNCA,-132))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW2",pmgdThis,15,BL_X_CC,BA_XCH,42))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW2",pmgdThis,16,BL_CH,BA_XCH,57))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW2",pmgdThis,16,BL_CH,BA_XCH,-64))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HW2",pmgdThis,16,BL_CH,BA_XCH,177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW3",pmgdThis,15,BL_X_CC,BA_XCH,-77))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW3",pmgdThis,16,BL_CH,BA_XCH,-176))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW3",pmgdThis,16,BL_CH,BA_XCH,65))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HW3",pmgdThis,16,BL_CH,BA_XCH,-55))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW1",pmgdThis,15,BL_X_CC,BA_XCH,161))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW1",pmgdThis,16,BL_CH,BA_XCH,172))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW1",pmgdThis,16,BL_CH,BA_XCH,-70))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CJ ",pmgdThis,16,BL_X_CC,BA_XCH,50))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HJ ",pmgdThis,17,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HJ ",pmgdThis,17,BL_CH,BA_XCH,177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CY2",pmgdThis,14,BL_X_CCDBL,BA_XNH,51))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW5",pmgdThis,15,BL_X_CC,BA_CNCA,1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW5",pmgdThis,16,BL_CH,BA_XCH,16))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW5",pmgdThis,16,BL_CH,BA_XCH,134))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3HW5",pmgdThis,16,BL_CH,BA_XCH,-107))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" CW4",pmgdThis,15,BL_X_CC,BA_CNCA,-177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1HW4",pmgdThis,16,BL_CH,BA_XCH,-140))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2HW4",pmgdThis,16,BL_CH,BA_XCH,104))!=ERR_SUCCESS) return err;
                  break;				
				default:;
			}
			break;
		case 'Y':
			if ((err=PlaceRotAtom("1HB ",pmgdThis,0,BL_CH,BA_XCH,chi1+120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom("2HB ",pmgdThis,0,BL_CH,BA_XCH,chi1-120))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CG ",pmgdThis,0,BL_Y_CBCG,BA_Y_CACBCG,chi1))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD1",pmgdThis,1,BL_Y_CGCD1,BA_Y_CBCGCD1,chi2))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HD1",pmgdThis,2,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE1",pmgdThis,2,BL_Y_CD1CE1,BA_Y_CGCD1CE1,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE1",pmgdThis,3,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CD2",pmgdThis,1,BL_Y_CGCD2,BA_Y_CBCGCD2,chi2+180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HD2",pmgdThis,2,BL_CH,BA_ARCH,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CE2",pmgdThis,2,BL_Y_CD2CE2,BA_Y_CGCD2CE2,180)) !=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HE2",pmgdThis,3,BL_CH,BA_ARCH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" CZ ",pmgdThis,3,BL_Y_CE2CZ,BA_Y_CD2CE2CZ,0))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" OH ",pmgdThis,4,BL_Y_CZOH,BA_Y_CE2CZOH,180))!=ERR_SUCCESS) return err;
			if ((err=PlaceRotAtom(" HH ",pmgdThis,5,BL_OH,BA_XOH,60))!=ERR_SUCCESS) return err;
			break;
		default:;
	}
	return ERR_SUCCESS;
}

/* unset all those possibly set by Place Rotamer */
void UndoBDRotamer(PMGD pmgdHere,Int2 Model)
{
	PVNMA pvnmaHere;
	PMAD pmadHere;
	vec vHere;
	CharPtr AtomName;

	pvnmaHere=pmgdHere->pvnmaAHead;
	while (pvnmaHere) {
		pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
		AtomName=pmadHere->pcAName;
		/* exceptions are CB and ..XT atoms */
		if ((StringRChr(scpos,(int)(AtomName[2]))!=NULL) && (StringCmp(AtomName," CB ")) && !(AtomName[2]=='X' && AtomName[3]=='T')) {
			if (GetCoOrds(pmgdHere,AtomName,vZero,vHere,Model)!=NULL) {
				if (IsAtomPlaced(pmadHere,Model))
					BDRemove(vHere,&(pwsThis->pbdTreeHead),&phbsHBondsToCheck);
			}
		}
		pvnmaHere=pvnmaHere->next;
	}
}

/* given CA and C on previous residue, places rotamer and HA */
TrajErr PlaceRotamer(PMGD pmgdLast,vec vHALast,Int2 Model)
{
	/*PVNMA pvnmaHere;*/
	PDNMG pdnmgLastLast;
	PMGD pmgdLastLast;
	/*PMAD pmadThis;*/
	/*PALD paldHere;*/
	PRS prsHere=NULL;
	PMAD pmadCA;
	vec vTmp,vTmp2,vNLast,vNLastNorm,vCLastNorm,vCBLast;
	Char residue;
	Int2 cnt,protries,rtries;
	TrajErr err;
	vec vTmpCG,vTmp3,vA3,vNLastAct,vCLastLast;
	FloatLo dihedmin=0.0,dihedmax=0.0,dihedhere,errhere=0.0,errmax,errmin;
	FloatLo chi1=0.0,chi2=0.0,chi3=0.0,chi4=0.0,chi1sd,chi2sd,arandnum,phihere=0.0,psihere=0.0;
	FloatLo dih,tmp;
	Boolean ishyp,ispyg;
	Char AtomNameHere[5];
	static Int2 numchis[20]={0,4,2,2,1,3,3,0,2,2,2,4,3,2,2,1,1,2,2,1};
/*vec vTmp4*/
/*FloatLo chithis1,chithis2,a1,a2,b1,b2,c1,d1,e1,f1,m1,m2,m3,m4,m5,tx,ty,absmag,del,err1,err2;*/
	PRS prsRotamerHere;
	FILE *fp;
	
	residue=GetAAFromIDict(pmgdLast);
	rtries=NUM_ROT_TRIES*numchis[GetResIdxFromMGD(aalist,pmgdLast)];
	if (rtries<1) rtries=1;
	/* choose rotamers here */
	if (residue!='G') 
		GetCoOrds(pmgdLast," CB ",vCALast,vCBLast,Model);
	else 
		GetCoOrds(pmgdLast,"2HA ",vCALast,vCBLast,Model);
	VecAdd(vTmp,vCBLast,vCALast);
	GetCoOrds(pmgdLast," N  ",vCALast,vNLast,Model);
	/* get rotamers for all but Ala and Gly */
	if ((residue!='G') && (residue!='A')) {
		/* get local phi, psi */
		pmadCA=FindAtomName(pmgdLast->pvnmaAHead," CA ");
		if (getR(pmadCA,&prsHere,Model)!=ERR_FAIL) {
			/* should only be NULL for first and last residue */
			if (prsHere!=NULL) {
				phihere=prsHere->Phi;
				psihere=prsHere->Psi;
				prsHere=MemFree(prsHere);
			}
			/* default to beta-sheet region if can't get phi, psi */
			else {
				phihere=-120;
				psihere=-120;
			}
		}
		else {
			phihere=-120;
			psihere=-120;
		}
	}
	/* check for proline and pyroglutamate */
	if (residue=='P' || pmgdLast->iIDict==326) {
		/* ishyp and ispyg are mutually exclusive */
		if ((pmgdLast->iIDict)/3==HYP_DICT_IDX/3)
			ishyp=TRUE;
		else
			ishyp=FALSE;
		if (pmgdLast->iIDict==326)
			ispyg=TRUE;
		else
			ispyg=FALSE;
/*VecAdd(vNLast,vNLast,vCALast);
CalcNextCoOrd(vNLast,vCALast,vTmp,vTmp2,BL_CH,BA_XCH,117.5);
if (!AssignCoords("1HB ",pmgdLast,vTmp2,Model)) return 0;
CalcNextCoOrd(vNLast,vCALast,vTmp,vTmp2,BL_CH,BA_XCH,-121.9);
if (!AssignCoords("2HB ",pmgdLast,vTmp2,Model)) return 0;
CalcNextCoOrd(vNLast,vCALast,vTmp,vTmpCG,BL_P_CBCG,BA_P_CACBCG,28);
if (!AssignCoords(" CG ",pmgdLast,vTmpCG,Model)) return 0;
CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_CH,BA_XCH,117.4);
if (!AssignCoords("1HG ",pmgdLast,vTmp2,Model)) return 0;
CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_CH,BA_XCH,-122);
if (!AssignCoords("2HG ",pmgdLast,vTmp2,Model)) return 0;
a1=BL_P_NCD*cos(DEGTORAD*BA_P_CANCD)*cos(DEGTORAD*(90-ba_ncacb[13]))-cos(DEGTORAD*(ba_ncacb[13]-90))*BL_P_NCA;
a2=BL_P_NCD*sin(DEGTORAD*BA_P_CANCD)*sin(DEGTORAD*(90-ba_ncacb[13]));
b1=BL_P_NCD*cos(DEGTORAD*BA_P_CANCD)*sin(DEGTORAD*(90-ba_ncacb[13]))+bl_cacb[13]+sin(DEGTORAD*(ba_ncacb[13]-90))*BL_P_NCA;
b2=-BL_P_NCD*sin(DEGTORAD*BA_P_CANCD)*cos(DEGTORAD*(90-ba_ncacb[13]));
c1=-BL_P_NCD*sin(DEGTORAD*BA_P_CANCD);
d1=-BL_P_CBCG*cos(DEGTORAD*(BA_P_CACBCG-90));
e1=-BL_P_CBCG*sin(DEGTORAD*(BA_P_CACBCG-90));
f1=-d1;
m1=-2*b2*(b1-e1)-2*a1*a2;
m2=2*a2*d1;
m3=2*c1*f1;
m4=a1*a1+(b1-e1)*(b1-e1)-BL_P_CGCD*BL_P_CGCD+d1*d1+c1*c1;
m5=-2*a1*d1;
tx=m1+m2*cos(DEGTORAD*28);
ty=m3*sin(-DEGTORAD*28);
absmag=sqrt(tx*tx+ty*ty);
del=acos(ty/absmag);
if ((sin(del)*(tx/absmag))<0)
	del=2*PI-del;
if (m4+m5*cos(DEGTORAD*28)>absmag) {
	PurgeGlobs();
	ErrPostEx(SEV_FATAL,1,1,"Cannot build PRO with given Chi1 angle (%f degrees)",28);
	return 0;
}				
chithis1=asin((m4+m5*cos(DEGTORAD*28))/absmag);
chithis2=PI-chithis1;
chithis1-=del;
chithis2-=del;
chithis1*=-RADTODEG;
chithis2*=-RADTODEG;
CalcNextCoOrd(vTmp,vCALast,vNLast,vTmp2,BL_P_NCD,BA_P_CANCD,chithis1);
VecSub(vTmp3,vTmp,vTmpCG);
VecSub(vTmp4,vTmp2,vTmpCG);
err1=fabs(BA_P_CBCGCD-RADTODEG*acos(Dot(vTmp3,vTmp4)/(getMag(vTmp3)*getMag(vTmp4))));
VecSub(vTmp3,vTmpCG,vTmp2);
VecSub(vTmp4,vNLast,vTmp2);
err1+=fabs(BA_P_CGCDN-RADTODEG*acos(Dot(vTmp3,vTmp4)/(getMag(vTmp3)*getMag(vTmp4))));
CalcNextCoOrd(vTmp,vCALast,vNLast,vTmp2,BL_P_NCD,BA_P_CANCD,chithis2);
VecSub(vTmp3,vTmp,vTmpCG);
VecSub(vTmp4,vTmp2,vTmpCG);
err2=fabs(BA_P_CBCGCD-RADTODEG*acos(Dot(vTmp3,vTmp4)/(getMag(vTmp3)*getMag(vTmp4))));
VecSub(vTmp3,vTmpCG,vTmp2);
VecSub(vTmp4,vNLast,vTmp2);
err2+=fabs(BA_P_CGCDN-RADTODEG*acos(Dot(vTmp3,vTmp4)/(getMag(vTmp3)*getMag(vTmp4))));
if (err2>err1)
CalcNextCoOrd(vTmp,vCALast,vNLast,vTmp2,BL_P_NCD,BA_P_CANCD,chithis1);
if (!AssignCoords(" CD ",pmgdLast,vTmp2,Model)) return 0;
CalcNextCoOrd(vTmp,vTmpCG,vTmp2,vTmp3,BL_CH,BA_XCH,125.2);
if (!AssignCoords("1HD ",pmgdLast,vTmp3,Model)) return 0;
CalcNextCoOrd(vTmp,vTmpCG,vTmp2,vTmp3,BL_CH,BA_XCH,-112.5);
if (!AssignCoords("2HD ",pmgdLast,vTmp3,Model)) return 0;
*/
		/* begin method 2 */
		protries=0;
		/* if cis- mean is -16, same s.d. */
		VecSub(vA3,vCALast,vTmp);
		VecAdd(vNLastAct,vNLast,vCALast);
		pdnmgLastLast=(pmgdLast->pdnmgLink)->last;
                if (pdnmgLastLast!=NULL) {
			pmgdLastLast=(PMGD)(pdnmgLastLast->data.ptrvalue);
                        GetCoOrds(pmgdLastLast," C  ",vZero,vCLastLast,Model);
                }
                else {
                        vCLastLast[0]=0.0;
                        vCLastLast[1]=0.0;
                        vCLastLast[2]=0.0;
		}
proagain:
		do {
			protries++;
			UndoBDRotamer(pmgdLast,Model);
			arandnum=fabs(Rand1());
			if (get_rotamers(pmgdLast,phihere,psihere,arandnum,&chi1,&chi2,&chi3,&chi4,&chi1sd,&chi2sd)!=ERR_SUCCESS) {
				PurgeGlobs();
				ErrPostEx(SEV_FATAL,2,13,"A serious error has occured with rotamer library, bailing..");
			}
			if (rotidHere!=0 && protries<=(rtries*CHANCE_TO_USE_FIXED_ROTAMER)) {
				prsRotamerHere=(PRS)MemNew(sizeof(RS));
			  	GetChiFromRotid(&prsRotamerHere,rotidHere);
				/* allow some chis to be chosen randomly */
				if (prsRotamerHere->Chi1!=0.0) {
				  	chi1=prsRotamerHere->Chi1;
			  		chi1sd=FIXED_ROTAMER_SD;
				}
				if (prsRotamerHere->Chi2!=0.0) {
				  	chi2=prsRotamerHere->Chi2;
			  		chi2sd=FIXED_ROTAMER_SD;
				}
				if (prsRotamerHere->Chi3!=0.0)
			  		chi3=prsRotamerHere->Chi3;
				if (prsRotamerHere->Chi4!=0.0)
				  	chi4=prsRotamerHere->Chi4;
				prsRotamerHere=MemFree(prsRotamerHere);
			}
			chi1=chi1+chi1sd*Rand1Distrib();
			CalcNextCoOrd(vNLastAct,vCALast,vTmp,vTmpCG,BL_P_CBCG,BA_P_CACBCG,chi1);
			if ((err=AssignCoords(" CG ",pmgdLast,vTmpCG,Model))!=ERR_SUCCESS) {
				errhere=BACKBONE_ERROR_TOLERANCE;
				if (err==ERR_CRASH)
					crashcnt--;
				else if (err==ERR_DISTCONST)
					dcbad--;
				continue;
			}
			dihedmax=chi2+3.0*chi2sd;
			dihedmin=chi2-3.0*chi2sd;
			CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_P_CGCD,BA_P_CBCGCD,dihedmax);
			errmax=GetProError(vTmp,vTmpCG,vTmp2,vNLastAct,vCLastLast,vA3);
			CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_P_CGCD,BA_P_CBCGCD,dihedmin);
			errmin=GetProError(vTmp,vTmpCG,vTmp2,vNLastAct,vCLastLast,vA3);
			while ((dihedmax-dihedmin)>BACKBONE_PRECISION) {
				dihedhere=0.5*(dihedmax+dihedmin);
				CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_P_CGCD,BA_P_CBCGCD,dihedhere);
				errhere=GetProError(vTmp,vTmpCG,vTmp2,vNLastAct,vCLastLast,vA3);
				if (errmax>errmin) {
					errmax=errhere;
					dihedmax=dihedhere;
				}
				else {
					errmin=errhere;
					dihedmin=dihedhere;
				}
			}
			/* a small hack to allow more tolerance for Pro-1 */
			if (pdnmgLastLast==NULL)
				errhere/=4.0;
			/* lots of tolerance if we dont care about residue */
			if (BUILD_FRAGS_ONLY && !ressrcres[(pmgdLast->pdnmgLink)->choice-1])
				errhere=0.0;
			/* error should be less than backbone error */
		} while ((protries<rtries) && (errhere>(BACKBONE_ERROR_TOLERANCE/4.0)));
		if (errhere>(BACKBONE_ERROR_TOLERANCE/4.0)) {
	       	        ErrPostEx(SEV_INFO,1,15,"Unable to place PRO atom, backing up, error is %f (degrees squared)",errhere);
			fpabad++;
			return ERR_BADBACK;
		}
		/* after rtries, since not large errhere, must be crashing */
		/* or not meeting a distant constraint */
		if (protries>rtries) {
			ErrPostEx(SEV_INFO,1,15,"Unable to place PRO atom, backing up");
			crashcnt++;
			return ERR_CRASH;
		}
		CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_P_CGCD,BA_P_CBCGCD,(dihedmin+dihedmax)/2.0);
		if (ishyp)
			StringCpy(AtomNameHere," CD2");
		else
			StringCpy(AtomNameHere," CD ");
		if ((err=AssignCoords(AtomNameHere,pmgdLast,vTmp2,Model))!=ERR_SUCCESS) {
			if (err==ERR_CRASH)
				crashcnt--;
			else if (err==ERR_DISTCONST)
				dcbad--;
			goto proagain;
		}
		/* put H's last so can form/detect H-bonds properly */
		GetDihedral(vTmp,vTmpCG,vTmp2,vNLastAct,-160,&dih,&tmp,&tmp,&tmp,&tmp,&tmp);
		if (!ispyg) {
			CalcNextCoOrd(vTmp,vTmpCG,vTmp2,vTmp3,BL_CH,BA_XCH,dih+120.0);
			if (ishyp)
				StringCpy(AtomNameHere,"1HD2");
			else
				StringCpy(AtomNameHere,"1HD ");
			if ((err=AssignCoords(AtomNameHere,pmgdLast,vTmp3,Model))!=ERR_SUCCESS) {
				if (err==ERR_CRASH)
					crashcnt--;
				else if (err==ERR_DISTCONST)
					dcbad--;
				goto proagain;
			}
		}
		if (ispyg)
			CalcNextCoOrd(vTmp,vTmpCG,vTmp2,vTmp3,BL_CO,BA_CACO,dih-180.0);
		else
			CalcNextCoOrd(vTmp,vTmpCG,vTmp2,vTmp3,BL_CH,BA_XCH,dih-120.0);
		if (ishyp)
			StringCpy(AtomNameHere,"2HD2");
		else if (ispyg)
			StringCpy(AtomNameHere," OE ");
		else
			StringCpy(AtomNameHere,"2HD ");
		if ((err=AssignCoords(AtomNameHere,pmgdLast,vTmp3,Model))!=ERR_SUCCESS) {
			if (err==ERR_CRASH)
				crashcnt--;
			else if (err==ERR_DISTCONST)
				dcbad--;
			goto proagain;
		}
		CalcNextCoOrd(vNLastAct,vCALast,vTmp,vTmp2,BL_CH,BA_XCH,chi1+120.0);
		if ((err=AssignCoords("1HB ",pmgdLast,vTmp2,Model))!=ERR_SUCCESS) {
			if (err==ERR_CRASH)
				crashcnt--;
			else if (err==ERR_DISTCONST)
				dcbad--;
			goto proagain;
		}
		CalcNextCoOrd(vNLastAct,vCALast,vTmp,vTmp2,BL_CH,BA_XCH,chi1-120.0);
		if ((err=AssignCoords("2HB ",pmgdLast,vTmp2,Model))!=ERR_SUCCESS) {
			if (err==ERR_CRASH)
				crashcnt--;
			else if (err==ERR_DISTCONST)
				dcbad--;
			goto proagain;
		}
		if (ishyp) {
			StringCpy(AtomNameHere," HG ");
			CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_COXT,BA_XCH,(dihedmin+dihedmax)/2.0+120.0);
		}
		else {
			StringCpy(AtomNameHere,"1HG ");
			CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_CH,BA_XCH,(dihedmin+dihedmax)/2.0+120.0);
		}
		if ((err=AssignCoords(AtomNameHere,pmgdLast,vTmp2,Model))!=ERR_SUCCESS) {
			if (err==ERR_CRASH)
				crashcnt--;
			else if (err==ERR_DISTCONST)
				dcbad--;
			goto proagain;
		}
		CalcNextCoOrd(vCALast,vTmp,vTmpCG,vTmp2,BL_CH,BA_XCH,(dihedmin+dihedmax)/2.0-120.0);
		if (ishyp)
			StringCpy(AtomNameHere," OD1");
		else
			StringCpy(AtomNameHere,"2HG ");
		if ((err=AssignCoords(AtomNameHere,pmgdLast,vTmp2,Model))!=ERR_SUCCESS) {
			if (err==ERR_CRASH)
				crashcnt--;
			else if (err==ERR_DISTCONST)
				dcbad--;
			goto proagain;
		}
		if (ishyp) {
			CalcNextCoOrd(vTmp,vTmpCG,vTmp2,vTmp3,BL_OH,BA_XOH,180.0);
			if ((err=AssignCoords(" HD1",pmgdLast,vTmp3,Model))!=ERR_SUCCESS) {
				if (err==ERR_CRASH)
					crashcnt--;
				else if (err==ERR_DISTCONST)
					dcbad--;
				goto proagain;
			}
		}
		/* now it's safe to normalize these vectors */
		Normalize(vCBLast,vCBLast);
		Normalize(vNLastNorm,vNLast);
	}
	else {
		Normalize(vCBLast,vCBLast);
		Normalize(vNLastNorm,vNLast);
		Normalize(vCLastNorm,vCLast);
/*		VecAdd(vTmp,vNLastNorm,vCLastNorm);
		VecScale(vZDir,vCBLast,1.0);
		Cross(vXDir,vZDir,vTmp);
		Normalize(vXDir,vXDir);
		Cross(vYDir,vZDir,vXDir);*/
		err=ERR_FAIL;
		if ((residue=='G') || (residue=='A')) {
			err=BuildSC(pmgdLast,vTmp,vNLast,chi1,chi2,chi3,chi4);
		}
		else {
			for (cnt=0;cnt<rtries;cnt++) {
				arandnum=fabs(Rand1());
				if (get_rotamers(pmgdLast,phihere,psihere,arandnum,&chi1,&chi2,&chi3,&chi4,&chi1sd,&chi2sd)!=ERR_SUCCESS) {
					PurgeGlobs();
					ErrPostEx(SEV_FATAL,2,13,"A serious error has occured with rotamer library, bailing..");
 				}
				if (rotidHere!=0 && cnt<=(rtries*CHANCE_TO_USE_FIXED_ROTAMER)) {
					prsRotamerHere=(PRS)MemNew(sizeof(RS));
				  	GetChiFromRotid(&prsRotamerHere,rotidHere);
					/* allow some chis to be chosen randomly */
					if (prsRotamerHere->Chi1!=0.0) {
					  	chi1=prsRotamerHere->Chi1;
				  		chi1sd=FIXED_ROTAMER_SD;
					}
					if (prsRotamerHere->Chi2!=0.0) {
					  	chi2=prsRotamerHere->Chi2;
				  		chi2sd=FIXED_ROTAMER_SD;
					}
					if (prsRotamerHere->Chi3!=0.0)
			  			chi3=prsRotamerHere->Chi3;
					if (prsRotamerHere->Chi4!=0.0)
					  	chi4=prsRotamerHere->Chi4;
					prsRotamerHere=MemFree(prsRotamerHere);
				}
				/* don't care what chi values are if Ala or Gly */
				/* remove "false" CB co-ordinates from the BD-tree before assigning true value */
/*				BDRemove(vTmp,&(pwsThis->pbdTreeHead),NULL);
*/				/* build side chain with any chi1, chi2, chi3, chi4 we choose */
				chi1=chi1+chi1sd*Rand1Distrib();
				chi2=chi2+chi2sd*Rand1Distrib();
				chi3=chi3+CHI34SD*Rand1Distrib();
				chi4=chi4+CHI34SD*Rand1Distrib();
				/* ensure chi2 of Phe, Tyr and Asp all lie
				   between -90 and +90, which is a standard
				   convention */
  				if ((residue=='Y') || (residue=='F') || (residue=='D')) {
	  				while (chi2<-90.0)
				        	chi2+=180.0;   
	    				while (chi2>90.0)
	      					chi2-=180.0;
				}
				err=BuildSC(pmgdLast,vTmp,vNLast,chi1,chi2,chi3,chi4);
				/* exit on success */
				if (err==ERR_SUCCESS) break;
				/* otherwise remove rotamer and try again */
				UndoBDRotamer(pmgdLast,Model);
				if (err==ERR_CRASH)
					crashcnt--;
				else if (err==ERR_DISTCONST)
					dcbad--;
			}
			if (err==ERR_CRASH)
				crashcnt++;
			else if (err==ERR_DISTCONST)
				dcbad++;
		}
		if (err!=ERR_SUCCESS) return err;
	}
/* if (phbsHBondsToCheck!=NULL) {*/
		/* unsatisfied H-bonds, must go back */
/*		phbsHBondsToCheck=PHBSFree(phbsHBondsToCheck);
		crashcnt++;
		return ERR_CRASH;							
	}*/
	/* find location of HA too */
	VecAdd(vTmp2,vNLastNorm,vCBLast);
	/* get correct C co-ordinate */
	GetCoOrds(pmgdLast," C  ",vZero,vCLast,Model);
	VecSub(vTmp,vCLast,vCALast);
	Normalize(vTmp,vTmp);
	VecAdd(vTmp2,vTmp,vTmp2);
	Normalize(vTmp2,vTmp2);
	NegateVec(vTmp2,vTmp2);
	VecScale(vHALast,vTmp2,BL_CH);	
	return ERR_SUCCESS;
}

TrajErr PutSpecNTerm(PMGD pmgdLast)
{
	TrajErr err;
	
	switch (pmgdLast->iIDict) {
                case 105:
                case 106:
                case 107:
                case 108:
                case 109:
                case 110:
                case 111:
                case 112:
                case 113:
                case 114:
                case 115:
                case 116:
                case 120:
                case 121:
                case 123:
                case 124:
                case 126:
                case 127:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH,BA_XNH,-18))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NC,BA_XNH,-160))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C2 ",pmgdLast,1,BL_X_CC,BA_E_CACBCG,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H2 ",pmgdLast,2,BL_CH,BA_XCH,-173))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H2 ",pmgdLast,2,BL_CH,BA_XCH,-54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H2 ",pmgdLast,2,BL_CH,BA_XCH,66))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" O2 ",pmgdLast,1,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                  break;
                case 122:
			if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NC,BA_XNH,-48))!=ERR_SUCCESS) return err;
   			if ((err=PlaceRotAtom(" C2 ",pmgdLast,1,BL_X_CC,BA_E_CACBCG,0))!=ERR_SUCCESS) return err;
   			if ((err=PlaceRotAtom("1H2 ",pmgdLast,2,BL_CH,BA_XCH,-175))!=ERR_SUCCESS) return err;
   			if ((err=PlaceRotAtom("2H2 ",pmgdLast,2,BL_CH,BA_XCH,-56))!=ERR_SUCCESS) return err;
   			if ((err=PlaceRotAtom("3H2 ",pmgdLast,2,BL_CH,BA_XCH,66))!=ERR_SUCCESS) return err;
   			if ((err=PlaceRotAtom(" O2 ",pmgdLast,1,BL_CO,BA_CACO,180))!=ERR_SUCCESS) return err;
                  break;
                case 125:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH,BA_XNH,-160))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NC,BA_XNH,-18))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C2 ",pmgdLast,1,BL_X_CC,BA_E_CACBCG,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H2 ",pmgdLast,2,BL_CH,BA_XCH,-173))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H2 ",pmgdLast,2,BL_CH,BA_XCH,-54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H2 ",pmgdLast,2,BL_CH,BA_XCH,66))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" O2 ",pmgdLast,1,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                  break;
                case 172:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_GUANNH,BA_XNH,-26))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NC,BA_XNH,154))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" H1 ",pmgdLast,1,BL_CH,BA_XCH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" O2 ",pmgdLast,1,BL_CO,BA_CACO,180))!=ERR_SUCCESS) return err;
                  break;
                case 222:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NC,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" O2 ",pmgdLast,1,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C2 ",pmgdLast,1,BL_X_CC,BA_E_CACBCG,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H2 ",pmgdLast,2,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H2 ",pmgdLast,2,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C3 ",pmgdLast,2,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H3 ",pmgdLast,3,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H3 ",pmgdLast,3,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C4 ",pmgdLast,3,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H4 ",pmgdLast,4,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H4 ",pmgdLast,4,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C5 ",pmgdLast,4,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H5 ",pmgdLast,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H5 ",pmgdLast,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C6 ",pmgdLast,5,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H6 ",pmgdLast,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H6 ",pmgdLast,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C7 ",pmgdLast,6,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H7 ",pmgdLast,7,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H7 ",pmgdLast,7,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C8 ",pmgdLast,7,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H8 ",pmgdLast,8,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H8 ",pmgdLast,8,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C9 ",pmgdLast,8,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H9 ",pmgdLast,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H9 ",pmgdLast,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C! ",pmgdLast,9,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H! ",pmgdLast,10,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H! ",pmgdLast,10,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C@ ",pmgdLast,10,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H@ ",pmgdLast,11,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H@ ",pmgdLast,11,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C# ",pmgdLast,11,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H# ",pmgdLast,12,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H# ",pmgdLast,12,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C$ ",pmgdLast,12,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H$ ",pmgdLast,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H$ ",pmgdLast,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C% ",pmgdLast,13,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H% ",pmgdLast,14,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H% ",pmgdLast,14,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C^ ",pmgdLast,14,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H^ ",pmgdLast,15,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H^ ",pmgdLast,15,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C& ",pmgdLast,15,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H& ",pmgdLast,16,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H& ",pmgdLast,16,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H& ",pmgdLast,16,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 229:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NC,BA_XNH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" O2 ",pmgdLast,1,BL_CO,BA_CACO,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C2 ",pmgdLast,1,BL_X_CC,BA_E_CACBCG,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H2 ",pmgdLast,2,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H2 ",pmgdLast,2,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C3 ",pmgdLast,2,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H3 ",pmgdLast,3,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H3 ",pmgdLast,3,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C4 ",pmgdLast,3,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H4 ",pmgdLast,4,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H4 ",pmgdLast,4,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C5 ",pmgdLast,4,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H5 ",pmgdLast,5,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H5 ",pmgdLast,5,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C6 ",pmgdLast,5,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H6 ",pmgdLast,6,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H6 ",pmgdLast,6,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C7 ",pmgdLast,6,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H7 ",pmgdLast,7,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H7 ",pmgdLast,7,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C8 ",pmgdLast,7,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H8 ",pmgdLast,8,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H8 ",pmgdLast,8,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C9 ",pmgdLast,8,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H9 ",pmgdLast,9,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H9 ",pmgdLast,9,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C! ",pmgdLast,9,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H! ",pmgdLast,10,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H! ",pmgdLast,10,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C@ ",pmgdLast,10,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H@ ",pmgdLast,11,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H@ ",pmgdLast,11,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C# ",pmgdLast,11,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H# ",pmgdLast,12,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H# ",pmgdLast,12,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C$ ",pmgdLast,12,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H$ ",pmgdLast,13,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H$ ",pmgdLast,13,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C% ",pmgdLast,13,BL_X_CC,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H% ",pmgdLast,14,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H% ",pmgdLast,14,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H% ",pmgdLast,14,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 242:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH,BA_XCH,-55))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NCA,BA_XCH,-176))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H1 ",pmgdLast,1,BL_CH,BA_XCH,173))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H1 ",pmgdLast,1,BL_CH,BA_XCH,54))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H1 ",pmgdLast,1,BL_CH,BA_XCH,-68))!=ERR_SUCCESS) return err;
                  break;
                case 243:
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NCA,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H1 ",pmgdLast,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H1 ",pmgdLast,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H1 ",pmgdLast,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C2 ",pmgdLast,0,BL_NCA,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H2 ",pmgdLast,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H2 ",pmgdLast,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H2 ",pmgdLast,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C3 ",pmgdLast,0,BL_NCA,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H3 ",pmgdLast,1,BL_CH,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H3 ",pmgdLast,1,BL_CH,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H3 ",pmgdLast,1,BL_CH,BA_XCH,180))!=ERR_SUCCESS) return err;
                  break;
                case 282:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH3,BA_XCH,60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NCA,BA_XCH,-60))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H1 ",pmgdLast,1,BL_CH,BA_XCH,176))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H1 ",pmgdLast,1,BL_CH,BA_XCH,-57))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H1 ",pmgdLast,1,BL_CH,BA_XCH,64))!=ERR_SUCCESS) return err;
                  break;
                case 283:
                        if ((err=PlaceRotAtom(" H  ",pmgdLast,0,BL_NH3,BA_XCH,151))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom(" C1 ",pmgdLast,0,BL_NCA,BA_XCH,-87))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H1 ",pmgdLast,1,BL_CH,BA_XCH,-177))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H1 ",pmgdLast,1,BL_CH,BA_XCH,-57))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("3H1 ",pmgdLast,1,BL_CH,BA_XCH,64))!=ERR_SUCCESS) return err;
                  break;
    default:;
	}	
	return ERR_SUCCESS;
}

/* given N2, C1 and CA1, fills in remainder of N-terminus */
TrajErr FindNTerm(PMGD pmgdLast,Int2 Model)
{
	FloatLo psi;
	Char residue;
	vec vNLast,vHALast,vCBLast,vHLast,vNHereAbs,vCLastAbs,vTmp;
	Int2 res;
	TrajErr err;

	residue=GetAAFromIDict(pmgdLast);
	res=GetResIdxFromMGD(NCBIstdaaUC,pmgdLast);
	psi=180.0*Rand1();
	VecAdd(vNHereAbs,vNHere,vCAHere);
	VecAdd(vCLastAbs,vCLast,vCALast);
	/* BL different but BA same for NH3 */
	CalcNextCoOrd(vNHereAbs,vCLastAbs,vCALast,vNLast,BL_NH3CA+BLSD_NH3CA*Rand1Distrib(),BA_NCAC+BASD_NCAC*Rand1Distrib(),psi);
	if ((err=AssignCoords(" N  ",pmgdLast,vNLast,Model))!=ERR_SUCCESS) return err;
	VecSub(vNLast,vNLast,vCALast);
	CalcNextCoOrd(vNHereAbs,vCLastAbs,vCALast,vCBLast,bl_cacb[res-1],ba_ccacb[res-1],psi+120);
	/* temporarily give it co-ords for use by Place Rotamer below */
	if (residue!='G') {
		if ((err=AssignCoords(" CB ",pmgdLast,vCBLast,Model))!=ERR_SUCCESS) return err;
	}
	else {
		if ((err=AssignCoords("2HA ",pmgdLast,vCBLast,Model))!=ERR_SUCCESS) return err;
	}
	VecSub(vCBLast,vCBLast,vCALast);
	Normalize(vCBLast,vCBLast);
	/* length of CB unimportant */
	/* temporary vCLast=vCLastAbs */
	VecScale(vTmp,vCLast,1.0);
	VecScale(vCLast,vCLastAbs,1.0);
	if ((err=PlaceRotamer(pmgdLast,vHALast,Model))!=ERR_SUCCESS) return err;
	VecScale(vCLast,vTmp,1.0);
	Translate(vHALast,vCALast);
	if (residue!='G') {
		if ((err=AssignCoords(" HA ",pmgdLast,vHALast,Model))!=ERR_SUCCESS) return err;
	}
	else {
		if ((err=AssignCoords("1HA ",pmgdLast,vHALast,Model))!=ERR_SUCCESS) return err;
	}
	Translate(vNLast,vCALast);
	if (IsNTermModified(pmgdLast)) {
		/* here we handle N-terminal modifications */
		VecScale(vRef[0],vCLastAbs,1.0);
		VecScale(vRef[1],vCALast,1.0);
		VecScale(vRef[2],vNLast,1.0);
		if ((err=PutSpecNTerm(pmgdLast))!=ERR_SUCCESS) return err;
	}
	else {
		/* put H3 part of NH3 */
		if (residue!='P' && (pmgdLast->iIDict!=326))
			CalcNextCoOrd(vCLastAbs,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,180);
		else
			CalcNextCoOrd(vCLastAbs,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,0);
		if ((err=AssignCoords("3H  ",pmgdLast,vHLast,Model))!=ERR_SUCCESS) return err;
		if (residue!='P' && (pmgdLast->iIDict!=326))
			CalcNextCoOrd(vCLastAbs,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,-60);
		else
			CalcNextCoOrd(vCLastAbs,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,-120);
		if ((err=AssignCoords("2H  ",pmgdLast,vHLast,Model))!=ERR_SUCCESS) return err;
		if (residue!='P' && (pmgdLast->iIDict!=326)) {
			CalcNextCoOrd(vCLastAbs,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,60);
			if ((err=AssignCoords("1H  ",pmgdLast,vHLast,Model))!=ERR_SUCCESS) return err;
		}
	}
/*	if (phbsHBondsToCheck!=NULL) {*/
		/* unsatisfied H-bonds, must go back */
/*		phbsHBondsToCheck=PHBSFree(phbsHBondsToCheck);
		crashcnt++;
		return ERR_CRASH;
	}*/
	return ERR_SUCCESS;
}

/* given Nn, Cn-1, CAn, fills in remainder of C-terminus */
TrajErr FindCTerm(PMGD pmgdNext,vec vNNext,vec vCHere,vec vCANext,Int2 Model)
{
	FloatLo phi,chi;
	Char residue;
	Int2 res;
	TrajErr err;
	vec vCNext,vHANext,vCBNext,vONext,vCHereAbs,vNNextAbs,vTmp1,vTmp2;

	residue=GetAAFromIDict(pmgdNext);
	res=GetResIdxFromMGD(NCBIstdaaUC,pmgdNext);
	/* last phi random */
	phi=180.0*Rand1();
	VecAdd(vNNextAbs,vNNext,vCANext);
	VecAdd(vCHereAbs,vCHere,vCAHere);
	CalcNextCoOrd(vCHereAbs,vNNextAbs,vCANext,vCNext,BL_CAC+BLSD_CAC*Rand1Distrib(),BA_NCAC+BASD_NCAC*Rand1Distrib(),phi);
	if ((err=AssignCoords(" C  ",pmgdNext,vCNext,Model))!=ERR_SUCCESS) return err;
	VecSub(vCNext,vCNext,vCANext);
	CalcNextCoOrd(vCHereAbs,vNNextAbs,vCANext,vCBNext,bl_cacb[res-1],ba_ncacb[res-1],phi-120);
	/* temporarily give it co-ords for use by Place Rotamer below */
	if (residue!='G') {
		if ((err=AssignCoords(" CB ",pmgdNext,vCBNext,Model))!=ERR_SUCCESS) return err;
	}
	else {
		if ((err=AssignCoords("2HA ",pmgdNext,vCBNext,Model))!=ERR_SUCCESS) return err;
	}
	VecSub(vCBNext,vCBNext,vCANext);
	Normalize(vCBNext,vCBNext);
	/* length of CB unimportant */
	/* temporary vCALast=vCANext, vCLast=vCNext */
	VecScale(vTmp1,vCALast,1.0);
	VecScale(vTmp2,vCLast,1.0);
	VecScale(vCALast,vCANext,1.0);
	VecScale(vCLast,vCNext,1.0);
	if ((err=PlaceRotamer(pmgdNext,vHANext,Model))!=ERR_SUCCESS) return err;
	VecScale(vCALast,vTmp1,1.0);
	VecScale(vCLast,vTmp2,1.0);
	Translate(vHANext,vCANext);
	if (residue!='G') {
		if ((err=AssignCoords(" HA ",pmgdNext,vHANext,Model))!=ERR_SUCCESS) return err;
	}
	else {
		if ((err=AssignCoords("1HA ",pmgdNext,vHANext,Model))!=ERR_SUCCESS) return err;
	}
	chi=180.0*Rand1();
	Translate(vCNext,vCANext);
	/* put OO part of COO */
	CalcNextCoOrd(vNNextAbs,vCANext,vCNext,vONext,BL_CO,BA_CACO,chi);
	if ((err=AssignCoords(" O  ",pmgdNext,vONext,Model))!=ERR_SUCCESS) return err;
	if (IsCTermModified(pmgdNext)) {
		/* here we handle C-terminal modifications */
		VecScale(vRef[0],vNNextAbs,1.0);
		VecScale(vRef[1],vCANext,1.0);
		VecScale(vRef[2],vCNext,1.0);
		switch (pmgdNext->iIDict) {		
                case 140:
                case 141:
                case 142:
                case 143:
                case 144:
                case 145:
                case 146:
                case 147:
                case 148:
                case 149:
                case 150:
                case 151:
                case 152:
                case 153:
                case 154:
                case 155:
                case 156:
                case 157:
                case 158:
                case 159:
                        CalcNextCoOrd(vONext,vRef[1],vRef[2],vRef[3],BL_NC,BA_CACOXT,180);
                        if ((err=AssignCoords(" N1 ",pmgdNext,vRef[3],1))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("1H1 ",pmgdNext,1,BL_GUANNH,BA_XNH,0))!=ERR_SUCCESS) return err;
                        if ((err=PlaceRotAtom("2H1 ",pmgdNext,1,BL_GUANNH,BA_XNH,180))!=ERR_SUCCESS) return err;
                  break;                 		                  		
  		default:;
		}	
	}
	else {
		CalcNextCoOrd(vNNextAbs,vCANext,vCNext,vONext,BL_COXT,BA_CACOXT,chi+180);
		if ((err=AssignCoords(" OXT",pmgdNext,vONext,Model))!=ERR_SUCCESS) return err;
	}
/*	if (phbsHBondsToCheck!=NULL) {*/
		/* unsatisfied H-bonds, must go back */
/*		phbsHBondsToCheck=PHBSFree(phbsHBondsToCheck);
		crashcnt++;
		return ERR_CRASH;							
	}*/
	return ERR_SUCCESS;
}

/* value if vTmp must not be altered!!! */
FloatLo GetError(vec vNLast,vec vCBLast,Int2 resLast,Int2 resHere,FloatLo testang,FloatLo r1,FloatLo r2,FloatLo Omega,FloatLo r4,FloatLo r5,Byte *isd) {
	vec vTmp,vCTmp;
	FloatLo tmp,angerror,thedot;

	VecAdd(vTmp,vCBLast,vCAHere);
	/* correct R&S - for cis-Pro */
	if (resHere==CISPRO)
		CalcNextCoOrd(vTmp,vCAHere,vCALast,vCLast,BL_CAC+r5,CIS_ETA+r1,testang);
	else
		CalcNextCoOrd(vTmp,vCAHere,vCALast,vCLast,BL_CAC+r5,RADTODEG*acos(coseta[resLast])+r1,testang);
	VecSub(vCTmp,vCLast,vCALast);
	if (NCBIstdaaUC[resLast]!='G') {
		thedot=Dot(vCTmp,vCBLast)/(BL_CAC+r5);
		if (thedot>1.00) thedot=1.00;
		if (thedot<-1.00) thedot=-1.00;
		tmp=RADTODEG*acos(thedot)-ba_ccacb[resLast-1];
	}
	else
		tmp=0.0;
	/* fudge factor accounts for uncertainty in CB */
	angerror=tmp*tmp*FUDGEFACT;
	thedot=Dot(vCTmp,vNLast)/(BL_CAC+r5);
	if (thedot>1.00) thedot=1.00;
	if (thedot<-1.00) thedot=-1.00;
	tmp=RADTODEG*acos(thedot)-BA_NCAC;
	angerror=angerror+tmp*tmp;
	if (resHere==CISPRO)
		CalcNextCoOrd(vCLast,vCALast,vCAHere,vNHere,BL_NCA+r4,CIS_ZETA+r2,Omega);
	else
		CalcNextCoOrd(vCLast,vCALast,vCAHere,vNHere,BL_NCA+r4,RADTODEG*acos(coszeta[resHere])+r2,Omega);
	VecSub(vTmp,vNHere,vCAHere);
	/* only NULL for last residue */
	if (((vCBHere[0]!=0.0) || (vCBHere[1]!=0.0) || (vCBHere[2]!=0.0)) && (NCBIstdaaUC[resHere]!='G')) {
		thedot=Dot(vTmp,vCBHere)/(BL_NCA+r4);
		if (thedot>1.00) thedot=1.00;
		if (thedot<-1.00) thedot=-1.00;
		tmp=RADTODEG*acos(thedot)-ba_ncacb[resHere-1];
		angerror=angerror+tmp*tmp*FUDGEFACT;
	}
	/* error term for CN bond length */
	VecSub(vTmp,vNHere,vCLast);
	/* SD of bond length in Angstroms should be about 1/150 that of angles in degrees */
	if ((NCBIstdaaUC[resHere]=='P') || (resHere==CISPRO))  /* use Pro NC length if appropriate */
		tmp=(getMag(vTmp)-BL_P_NC)*150.0;
	else
		tmp=(getMag(vTmp)-BL_NC)*150.0;
	angerror=angerror+tmp*tmp;
	Cross(vTmp,vCBLast,vCTmp);
	/* forbid D-amino acids */
	*isd=0;
	if ((NCBIstdaaUC[resLast]!='G') && (Dot(vTmp,vNLast)>0)) {
		*isd=1;
	}
	return angerror;
}

/* value if vTmp must not be altered!!! */
/*FloatLo iGetError(vec vNLast,vec vCBLast,Int2 resLast,Int2
resHere,FloatLo testang,FloatLo r1,FloatLo r2,FloatLo r3,Byte *isd)
{
	vec vTmp,vCTmp;
	FloatLo tmp,angerror;

	VecAdd(vTmp,vCBLast,vCAHere);
	if (resHere==CISPRO)
		iCalcNextCoOrd(vTmp,vCAHere,vCALast,vCLast,BL_CAC,CIS_ETA+r1,testang);
	else
		iCalcNextCoOrd(vTmp,vCAHere,vCALast,vCLast,BL_CAC,RADTODEG*acos(coseta[resLast])+r1,testang);
	VecSub(vCTmp,vCLast,vCALast);
	if (NCBIstdaaUC[resLast]!='G')
		tmp=RADTODEG*acos(Dot(vCTmp,vCBLast)/BL_CAC)-ba_ccacb[resLast-1];
	else
		tmp=0.0;
	angerror=tmp*tmp*FUDGEFACT;
	tmp=RADTODEG*acos(Dot(vCTmp,vNLast)/BL_CAC)-BA_NCAC;
	angerror=angerror+tmp*tmp;
	if (resHere==CISPRO)
		iCalcNextCoOrd(vCLast,vCALast,vCAHere,vNHere,BL_NCA,CIS_ZETA+r2,r3);
	else
		iCalcNextCoOrd(vCLast,vCALast,vCAHere,vNHere,BL_NCA,RADTODEG*acos(coszeta[resHere])+r2,CHI_W+r3);
	VecSub(vTmp,vNHere,vCAHere);
	if ((vCBHere!=NULL) && (NCBIstdaaUC[resHere]!='G')) {
		tmp=RADTODEG*acos(Dot(vTmp,vCBHere)/BL_NCA)-ba_ncacb[resHere-1];
		angerror=angerror+tmp*tmp*FUDGEFACT;
	}
	VecSub(vTmp,vNHere,vCLast);
	if ((NCBIstdaaUC[resHere]=='P') || (resHere==CISPRO))
		tmp=(getMag(vTmp)-BL_P_NC)*150.0;
	else
		tmp=(getMag(vTmp)-BL_NC)*150.0;
	angerror=angerror+tmp*tmp;
	Cross(vTmp,vCBLast,vCTmp);
	*isd=0;
	if ((NCBIstdaaUC[resLast]!='G') && (Dot(vTmp,vNLast)>0)) {
		*isd=1;
	}
	return angerror;
}
*/

/* attempt to find optimal solution for backbone C and N atoms using other nearby atom
   locations */
TrajErr FindPeptideAtoms(PMGD pmgdLast,PMGD pmgdHere,Int2 Model, FloatLo Omega)
{
	vec vCBLast,vNLast,vTmp,vOLast;
	Char residueLast,residueHere;
	FloatLo minerr,minerrang,testhere,testhi,testlo,errhere,errhi,errlo,incsize,testtmp,r1,r2,r4,r5;
	Int2 resLast,resHere,looptries;
	Byte isd;
	TrajErr err;

	looptries=0;
	residueLast=GetAAFromIDict(pmgdLast);
	residueHere=GetAAFromIDict(pmgdHere);
	resLast=GetResIdxFromMGD(NCBIstdaaUC,pmgdLast);
	resHere=GetResIdxFromMGD(NCBIstdaaUC,pmgdHere);
	/* check for cis-Pro */
	if ((residueLast=='P') && (cisLast))
		resLast=CISPRO;
	if ((residueHere=='P') && (cisHere))
		resHere=CISPRO;
	if (residueLast!='G') 
		GetCoOrds(pmgdLast," CB ",vCALast,vCBLast,Model);
	else 
		GetCoOrds(pmgdLast,"2HA ",vCALast,vCBLast,Model);
	Normalize(vCBLast,vCBLast);
	GetCoOrds(pmgdLast," N  ",vCALast,vNLast,Model);
	Normalize(vNLast,vNLast);
        /* generate uncertainty in eta, 1.5 degrees standard deviation */
        r1=1.5*Rand1Distrib();
        /* uncertainty in zeta, 1.5 degrees s.d. */
        r2=1.5*Rand1Distrib();
        /* this dihedral (Ci-CAi-CA(i+1)-N(i+1)) actually varies about 3 times as
	   much as omega */
	r4=BLSD_NCA*Rand1Distrib();
	r5=BLSD_CAC*Rand1Distrib();
	/* perform binary search for max and min. */
	incsize=INCSIZE;
	testhere=START_BACKBONE;
	testhi=testhere+incsize;
	testlo=testhere-incsize;
	errhere=GetError(vNLast,vCBLast,resLast,resHere,testhere,r1,r2,Omega,r4,r5,&isd);	
	errhi=GetError(vNLast,vCBLast,resLast,resHere,testhi,r1,r2,Omega,r4,r5,&isd); 
	errlo=GetError(vNLast,vCBLast,resLast,resHere,testlo,r1,r2,Omega,r4,r5,&isd); 
recalc1:
	looptries++;
	/* prevent endless loop */
	if (looptries==TRIES_MAX)
                ErrPostEx(SEV_WARNING,1,15,"Avoiding endless loop");
	if ((errhi<errhere) && (errhere<errlo)) {
		/* errlo OK, increase errhi */
		do { 
			testhi+=incsize;
			testtmp=errhi;
			errhi=GetError(vNLast,vCBLast,resLast,resHere,testhi,r1,r2,Omega,r4,r5,&isd); 
		} while (errhi<testtmp); 
	}
	else if ((errhi>errhere) && (errhere>errlo)) {
		/* errhi OK, lower errlo */
		do { 
			testlo-=incsize;
			testtmp=errlo;
			errlo=GetError(vNLast,vCBLast,resLast,resHere,testlo,r1,r2,Omega,r4,r5,&isd); 
		} while (errlo<testtmp); 
	}
	else if ((errhi<errhere) && (errlo<errhere) && (looptries<TRIES_MAX)) {
		/* errhere near a maximum, not minimum!! */
		if (errhi>errlo) {
			/* probably want minimum closest to errlo */
			testhi=testlo;
			testhere=testhi-incsize;
			testlo=testhere-incsize;
			errhi=errlo;
			errhere=GetError(vNLast,vCBLast,resLast,resHere,testhere,r1,r2,Omega,r4,r5,&isd); 
			errlo=GetError(vNLast,vCBLast,resLast,resHere,testlo,r1,r2,Omega,r4,r5,&isd); 
			goto recalc1;
		}
		else {
			/* probably want minimum closest to errhi */
			testlo=testhi;
			testhere=testlo+incsize;
			testhi=testhere+incsize;
			errlo=errhi;
			errhere=GetError(vNLast,vCBLast,resLast,resHere,testhere,r1,r2,Omega,r4,r5,&isd); 
			errhi=GetError(vNLast,vCBLast,resLast,resHere,testhi,r1,r2,Omega,r4,r5,&isd); 
			goto recalc1;
		}
	}
	/* now errhere<errhi and errhere<errlo */
	while ((testhi-testlo)>BACKBONE_PRECISION) {
		errhere=GetError(vNLast,vCBLast,resLast,resHere,(testlo+testhi)/2.0,r1,r2,Omega,r4,r5,&isd);
		if (errhi<errlo) {
			errlo=errhere;
 			testlo=(testlo+testhi)/2.0;
		}
		else {
			errhi=errhere;
			testhi=(testlo+testhi)/2.0;
		}
	}	
	minerr=errhere;
	/* forbid D- amino acids */
	if (isd) minerr+=999999.0;
	minerrang=(testlo+testhi)/2.0;
	/* start far away from first minimum */
	testhere=minerrang+180.0;
	testhi=testhere+incsize;
	testlo=testhere-incsize;
	looptries=0;
	errhere=GetError(vNLast,vCBLast,resLast,resHere,testhere,r1,r2,Omega,r4,r5,&isd);
	errhi=GetError(vNLast,vCBLast,resLast,resHere,testhi,r1,r2,Omega,r4,r5,&isd); 
	errlo=GetError(vNLast,vCBLast,resLast,resHere,testlo,r1,r2,Omega,r4,r5,&isd); 
recalc2:
	looptries++;
	/* prevent endless loop */
	if (looptries==TRIES_MAX)
                ErrPostEx(SEV_WARNING,1,15,"Avoiding endless loop");
	if ((errhi<errhere) && (errhere<errlo)) {
		/* errlo OK, increase errhi */
		do { 
			testhi+=incsize;
			testtmp=errhi;
			errhi=GetError(vNLast,vCBLast,resLast,resHere,testhi,r1,r2,Omega,r4,r5,&isd); 
		} while (errhi<testtmp); 
	}
	else if ((errhi>errhere) && (errhere>errlo)) {
		/* errhi OK, lower errlo */
		do { 
			testlo-=incsize;
			testtmp=errlo;
			errlo=GetError(vNLast,vCBLast,resLast,resHere,testlo,r1,r2,Omega,r4,r5,&isd); 
		} while (errlo<testtmp); 
	}
	else if ((errhi<errhere) && (errlo<errhere) && (looptries<TRIES_MAX)) {
		/* errhere near a maximum, not minimum!! */
		if (errhi>errlo) {
			/* probably want minimum closest to errlo */
			testhi=testlo;
			testhere=testhi-incsize;
			testlo=testhere-incsize;
			errhi=errlo;
			errhere=GetError(vNLast,vCBLast,resLast,resHere,testhere,r1,r2,Omega,r4,r5,&isd); 
			errlo=GetError(vNLast,vCBLast,resLast,resHere,testlo,r1,r2,Omega,r4,r5,&isd); 
			goto recalc2;
		}
		else {
			/* probably want minimum closest to errhi */
			testlo=testhi;
			testhere=testlo+incsize;
			testhi=testhere+incsize;
			errlo=errhi;
			errhere=GetError(vNLast,vCBLast,resLast,resHere,testhere,r1,r2,Omega,r4,r5,&isd); 
			errhi=GetError(vNLast,vCBLast,resLast,resHere,testhi,r1,r2,Omega,r4,r5,&isd); 
			goto recalc2;
		}
	}
	/* now errhere<errhi and errhere<errlo */
	while ((testhi-testlo)>BACKBONE_PRECISION) {
		errhere=GetError(vNLast,vCBLast,resLast,resHere,(testlo+testhi)/2.0,r1,r2,Omega,r4,r5,&isd);
		if (errhi<errlo) {
			errlo=errhere;
 			testlo=(testlo+testhi)/2.0;
		}
		else {
			errhi=errhere;
			testhi=(testlo+testhi)/2.0;
		}
	}
	/* forbid D- amino acids */	
	if (isd) errhere+=999999.0;
	if (errhere<minerr) {
		minerr=errhere;
		minerrang=(testlo+testhi)/2.0;
	}
	if (!BUILD_FRAGS_ONLY || ressrcres[(pmgdHere->pdnmgLink)->choice-1]) {
		if (minerr>BACKBONE_ERROR_TOLERANCE) {
					ErrPostEx(SEV_INFO,1,15,"Unable to place atom, backing up, error is %f (degrees squared)",minerr);
			fpabad++;
	/*printf("backbone error: %f\n",minerr);*/
			return ERR_BADBACK;
		}
	}
	VecAdd(vTmp,vCBLast,vCAHere);
	if (resHere==CISPRO)
		CalcNextCoOrd(vTmp,vCAHere,vCALast,vCLast,BL_CAC+r5,CIS_ETA+r1,minerrang);
	else
		CalcNextCoOrd(vTmp,vCAHere,vCALast,vCLast,BL_CAC+r5,RADTODEG*acos(coseta[resLast])+r1,minerrang);
	if ((err=AssignCoords(" C  ",pmgdLast,vCLast,Model))!=ERR_SUCCESS) return err;
	/* try slightly random dihedrals */
	if (resHere==CISPRO)
		CalcNextCoOrd(vCLast,vCALast,vCAHere,vNHere,BL_NCA+r4,CIS_ZETA+r2,Omega);
	else
		CalcNextCoOrd(vCLast,vCALast,vCAHere,vNHere,BL_NCA+r4,RADTODEG*acos(coszeta[resHere])+r2,Omega);
	if ((err=AssignCoords(" N  ",pmgdHere,vNHere,Model))!=ERR_SUCCESS) return err;
	if (resHere==CISPRO)
		CalcNextCoOrd(vNHere,vCAHere,vCALast,vOLast,CIS_OCA,CIS_OCACA,Omega);
	else
		CalcNextCoOrd(vNHere,vCAHere,vCALast,vOLast,bl_cao[resLast],ba_ocaca[resLast],Omega);
	if ((err=AssignCoords(" O  ",pmgdLast,vOLast,Model))!=ERR_SUCCESS) return err;
	return ERR_SUCCESS;
}

void FindCBDirPhiPsi(vec vNLast,vec vCLast,vec vCALast,vec vCBLast,FloatLo ncacb,FloatLo ccacb)
{
	vec vN,vC,r1,r2,vTmp;
	FloatHi a1,a1sq,a2,a3,a4,a5,a6,a,b,c,d,cosa,cosb;

	cosa=cos(DEGTORAD*ncacb);
	cosb=cos(DEGTORAD*ccacb);
	VecSub(vN,vNLast,vCALast);
	VecSub(vC,vCLast,vCALast);
	Normalize(vN,vN);
	Normalize(vC,vC);
	a1=vC[2]*vN[0]-vC[0]*vN[2];
	a2=vC[2]*vN[1]-vC[1]*vN[2];
	a3=vC[2]*cosa-vN[2]*cosb;
	if (fabs(a1)<CB_PRECISION) {
		/* only 1 root */
		if (fabs(a2)<CB_PRECISION) {
			/* here we know (C x N)z must be non-zero (assigned to a2
			   here) since N and C are not perpendicular, so safe to
			   divide by it */
			a2=vC[0]*vN[1]-vC[1]*vN[0];
			a3=vC[0]*cosa-vN[0]*cosb;
			r1[1]=a3/a2;
			r2[1]=a3/a2;
			/* also Nz and Cz must be zero */
			/* one of Nx or Cx is non-zero, else parallel! */
			if (fabs(vN[0])<CB_PRECISION) {
				/* use Cx */
				r1[0]=(cosb-vC[1]*r1[1])/vC[0];
				r2[0]=r1[0];
			}
			else {
				r1[0]=(cosa-vN[1]*r1[1])/vN[0];
				r2[0]=r1[0];
			}
			/* get right sign for rz */
			Cross(vTmp,vN,vC);
			r1[2]=sqrt(1.0-r1[0]*r1[0]-r1[1]*r1[1]);
			r2[2]=sqrt(1.0-r2[0]*r2[0]-r2[1]*r2[1]);
			if (vTmp[2]<0.0) {
				r1[2]=-r1[2];
				r2[2]=-r2[2];
			}
		}
		else {
			/* if a1 = 0, a2 != 0 */
			r1[1]=a3/a2;
			r2[1]=a3/a2;
			/* need x and z=1-x*x-y*y */
			a=vN[0]*vN[0]+vN[2]*vN[2];
			if (fabs(a)<CB_PRECISION) {
				/* a is only 0 iff N along y axis and C not in x-y plane */
				/* thus Cz if non-zero and this a is non-zero */
				a=vC[0]*vC[0]+vC[2]*vC[2];
				b=(2.0*vC[0]*(vC[1]*r1[1]-cosb))/a;
				c=(cosb*cosb-vC[2]*vC[2]-r1[1]*(2.0*vC[1]*cosb)+r1[1]*r1[1]*(vC[1]*vC[1]+vC[2]*vC[2]))/a;
			}
			else {
				b=(2.0*vN[0]*(vN[1]*r1[1]-cosa))/a;
				c=(cosa*cosa-vN[2]*vN[2]-r1[1]*(2.0*vN[1]*cosa)+r1[1]*r1[1]*(vN[1]*vN[1]+vN[2]*vN[2]))/a;
			}
			d=b*b-4.0*c;
			/* should never happen since there is always at least one solution */
			if (d<0.0) {
				ErrPostEx(SEV_ERROR,1,13,"Discriminant <0 (%f) while calculating CB",d);
				d=0.0;
			}
			/* r1 is always the larger root in absolute magnitude */
			if (b>0.0)
				r1[0]=(-b-sqrt(d))/(2.0);
			else
				r1[0]=(-b+sqrt(d))/(2.0);
			/* should be impossible but just in case */
			if (fabs(r1[0])<CB_PRECISION) 
				r2[0]=r1[0];
			else
				r2[0]=c/r1[0];	
			/* check if Nz=0 */
			if (fabs(vN[2])<CB_PRECISION) {
				/* Cz not 0 or else a2=0 up above */
				r1[2]=(cosb-vC[0]*r1[0]-vC[1]*r1[1])/vC[2];
				r2[2]=(cosb-vC[0]*r2[0]-vC[1]*r2[1])/vC[2];
			}
			else {
				r1[2]=(cosa-vN[0]*r1[0]-vN[1]*r1[1])/vN[2];
				r2[2]=(cosa-vN[0]*r2[0]-vN[1]*r2[1])/vN[2];
			}
		}
	}
	else {
		a1sq=a1*a1;
		a4=vN[2]*(vN[0]*cosb-vC[0]*cosa)/a1;
		a5=vN[2]*(vC[0]*vN[1]-vC[1]*vN[0])/a1;
		a6=vN[2]*vN[2];
		if (fabs(vN[2])<CB_PRECISION) {
			/* since a1!=0, but Nz=0, we know Nx!=0 */
			/* use Nx */
			a1=-a1;
			a2=vC[0]*vN[1]-vC[1]*vN[0];
			a3=vC[0]*cosa-vN[0]*cosb;
			a4=-cosa;
			a5=vN[1];
			a6=vN[0]*vN[0];
		}
		a=a5*a5+a6*(1.0+a2*a2/a1sq);
		/* normalize quadratic */
		/* a is never zero */
		b=(2.0*(a4*a5-a6*a3*a2/a1sq))/a;
		c=(a6*(a3*a3/a1sq-1.0)+a4*a4)/a;
		d=b*b-4.0*c;
		/* should never happen since there is always at least one solution */
		if (d<0.0) {
			ErrPostEx(SEV_ERROR,1,13,"Discriminant =%f while calculating CB",d);
			d=0.0;
		}
		/* r1 is always the larger root in absolute magnitude */
		if (b>0.0)
			r1[1]=(-b-sqrt(d))/(2.0);
		else
			r1[1]=(-b+sqrt(d))/(2.0);
		/* should be impossible but just in case */
		if (fabs(r1[1])<CB_PRECISION) 
			r2[1]=r1[1];
		else
			r2[1]=c/r1[1];	
		if (fabs(vN[2])<CB_PRECISION) {
			/* since a1!=0, but Nz=0, we know Nx!=0 */
			r1[2]=(a3-a2*r1[1])/a1;
			r2[2]=(a3-a2*r2[1])/a1;
			r1[0]=(cosa-vN[2]*r1[2]-vN[1]*r1[1])/vN[0];
			r2[0]=(cosa-vN[2]*r2[2]-vN[1]*r2[1])/vN[0];
		}
		else {
			r1[0]=(a3-a2*r1[1])/a1;
			r2[0]=(a3-a2*r2[1])/a1;
			r1[2]=(cosa-vN[0]*r1[0]-vN[1]*r1[1])/vN[2];
			r2[2]=(cosa-vN[0]*r2[0]-vN[1]*r2[1])/vN[2];
		}
	}
	Normalize(r1,r1);
	Normalize(r2,r2);
	Cross(vTmp,vN,vC);
	Normalize(vTmp,vTmp);
	if (Dot(vTmp,r1)>0) {
		vCBLast[0]=r1[0];
		vCBLast[1]=r1[1];
		vCBLast[2]=r1[2];
	}
	else {
		vCBLast[0]=r2[0];
		vCBLast[1]=r2[1];
		vCBLast[2]=r2[2];
	}
}

TrajErr SetAllCoordsPhiPsi(PMGD pmgdHere,Int2 Model, FloatLo Phi, FloatLo Psi,FloatLo Omega)
{
        PMGD pmgdLast,pmgdLastLast=NULL;
        PDNMG pdnmgLastLast,pdnmgLast,pdnmgHere;
	Char residueHere,residueLast;
	Int2 resHere,resLast,resnum;
	Int4 clength;
	PMMD pmmdParent;
	vec vNLast,vCLastLast,vHALast,vTmp,vOLast,vHHere,vHLast,vCBLast;
	FloatLo psihere;
	
	pmmdParent=(PMMD)(pmgdHere->pfbParent);
	clength=pmmdParent->iResCount;
        pdnmgHere=pmgdHere->pdnmgLink;
	resnum=pdnmgHere->choice;
        pdnmgLast=pdnmgHere->last;
	if (resnum>2) {
	        pdnmgLastLast=pdnmgLast->last;
	        if ((pdnmgLast->choice-pdnmgLastLast->choice!=1)) {
	                ErrPostEx(SEV_ERROR,1,1,"Residues out of order, %d %d",pdnmgLast->choice,pdnmgLastLast->choice);
	                return ERR_FAIL;
	        }
	        pmgdLastLast=(PMGD)(pdnmgLastLast->data.ptrvalue);
	}
        if (pdnmgLast==NULL) {
                ErrPostEx(SEV_ERROR,1,1,"Unable to find C beta for first residue");
                return ERR_FAIL;
        }
        /* residues must be ordered consecutively in the linked list */
        if ((pdnmgHere->choice-pdnmgLast->choice!=1)) {
                ErrPostEx(SEV_ERROR,1,1,"Residues out of order, %d %d",pdnmgHere->choice,pdnmgLast->choice);
                return ERR_FAIL;
        }
        pmgdLast=(PMGD)(pdnmgLast->data.ptrvalue);
	residueHere=GetAAFromIDict(pmgdHere);
	residueLast=GetAAFromIDict(pmgdLast);
	resHere=GetResIdxFromMGD(NCBIstdaaUC,pmgdHere);
	resLast=GetResIdxFromMGD(NCBIstdaaUC,pmgdLast);
	/* take care of new cis-Proline residues */
	if ((residueHere=='P') && (cisHere))
		resHere=CISPRO;
	if ((residueLast=='P') && (cisLast))
		resLast=CISPRO;
	/* find C using Phi */
	GetCoOrds(pmgdLast," N  ",vZero,vNLast,Model);
	GetCoOrds(pmgdLast," CA ",vZero,vCALast,Model);
	if (resnum>2) {
		GetCoOrds(pmgdLastLast," C  ",vZero,vCLastLast,Model);
		CalcNextCoOrd(vCLastLast,vNLast,vCALast,vCLast,BL_CAC+Rand1Distrib()*BLSD_CAC,BA_NCAC+Rand1Distrib()*BASD_NCAC,Phi);		
		if (CheckEndCDistConstraints(pmgdLast,vCLast,Model,&dcviol)!=ERR_SUCCESS)
			return ERR_FAIL;			
		if (AssignCoords(" C  ",pmgdLast,vCLast,Model)!=ERR_SUCCESS) return ERR_FAIL;
	}
	else {
		/* complete N-terminus atoms */
		GetCoOrds(pmgdLast," C  ",vZero,vCLast,Model);
		if (IsNTermModified(pmgdLast)) {
			/* here we handle N-terminal modifications */
			VecScale(vRef[0],vCLast,1.0);
			VecScale(vRef[1],vCALast,1.0);
			VecScale(vRef[2],vNLast,1.0);
			if (PutSpecNTerm(pmgdLast)!=ERR_SUCCESS) return ERR_FAIL;
		}
		else {
			/* put H3 part of NH3 */
			if (residueLast!='P' && (pmgdLast->iIDict!=326))
				CalcNextCoOrd(vCLast,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,180.0);
			else
				CalcNextCoOrd(vCLast,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,0.0);
			if (AssignCoords("3H  ",pmgdLast,vHLast,Model)!=ERR_SUCCESS) return ERR_FAIL;
			if (residueLast!='P' && (pmgdLast->iIDict!=326))
				CalcNextCoOrd(vCLast,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,-60.0);
			else	
				CalcNextCoOrd(vCLast,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,-120.0);
			if (AssignCoords("2H  ",pmgdLast,vHLast,Model)!=ERR_SUCCESS) return ERR_FAIL;
			if (residueLast!='P' && (pmgdLast->iIDict!=326)) {
				CalcNextCoOrd(vCLast,vCALast,vNLast,vHLast,BL_NH3H,BA_CANH3H,60.0);
				if (AssignCoords("1H  ",pmgdLast,vHLast,Model)!=ERR_SUCCESS) return ERR_FAIL;
			}
		}
/*		if (phbsHBondsToCheck!=NULL) {*/
			/* unsatisfied H-bonds, must go back */
/*			phbsHBondsToCheck=PHBSFree(phbsHBondsToCheck);
			crashcnt++;
			return ERR_FAIL;							
		}*/
	}
	/* calculate CB direction vector */
	FindCBDirPhiPsi(vNLast,vCLast,vCALast,vCBLast,ba_ncacb[resLast-1],ba_ccacb[resLast-1]);
	/* and assign it for use with Place Rotamer later */
	VecScale(vTmp,vCBLast,bl_cacb[resLast-1]);
	VecAdd(vTmp,vTmp,vCALast);
	if (residueLast!='G') {
		if ((AssignCoords(" CB ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	else {
		if ((AssignCoords("2HA ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	/* find N using Psi */
	psihere=Psi;
	if (resnum==2)
		psihere=180.0*Rand1();
	if (residueHere=='P')
		CalcNextCoOrd(vNLast,vCALast,vCLast,vNHere,BL_P_NC+Rand1Distrib()*BLSD_NC,BA_P_CACN+Rand1Distrib()*BASD_P_CACN,psihere);
	else
		CalcNextCoOrd(vNLast,vCALast,vCLast,vNHere,BL_NC+Rand1Distrib()*BLSD_NC,BA_CACN+Rand1Distrib()*BASD_CACN,psihere);		
	/* if coming out of a gap distance constraint, AssignCoords checks D, Angle1 and dihedral01 */	
	if (AssignCoords(" N  ",pmgdHere,vNHere,Model)!=ERR_SUCCESS) return ERR_FAIL;
	/* place rotamer */
	if ((PlaceRotamer(pmgdLast,vHALast,Model))!=ERR_SUCCESS) return ERR_FAIL;
	/* place HA i-1 */
	VecAdd(vTmp,vHALast,vCALast);
	if (residueLast!='G') {
		if ((AssignCoords(" HA ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	else {
		if ((AssignCoords("1HA ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	/* find CAHere using omega */
	CalcNextCoOrd(vCALast,vCLast,vNHere,vCAHere,BL_NCA+Rand1Distrib()*BLSD_NCA,BA_CNCA+Rand1Distrib()*BASD_CNCA,Omega);	
	if (CheckEndCADistConstraints(pmgdHere,vCAHere,Model,&dcviol)!=ERR_SUCCESS)
		return ERR_FAIL;
	/* may help a bit - check if CA position forbids us from our goal */
	if (CheckDistConstraints(pmgdHere,vCAHere,Model,&dcviol)!=ERR_SUCCESS)
		return ERR_FAIL;
	if (AssignCoords(" CA ",pmgdHere,vCAHere,Model)!=ERR_SUCCESS) return ERR_FAIL;
	/* find OLast */
	if (resHere==CISPRO)
		CalcNextCoOrd(vNHere,vCALast,vCLast,vOLast,BL_CO,BA_CACO,180.0+Omega);
	else
		CalcNextCoOrd(vNHere,vCALast,vCLast,vOLast,BL_CO,BA_CACO,Omega);
	if (AssignCoords(" O  ",pmgdLast,vOLast,Model)!=ERR_SUCCESS) return ERR_FAIL;
	/* find HHere */	
	/* place H if not proline */
	if (residueHere!='P') {
		/* vNHere set for carbonyl O above */
		VecSub(vTmp,vCLast,vNHere);
		VecSub(vHHere,vCAHere,vNHere);
		Normalize(vHHere,vHHere);
		Normalize(vTmp,vTmp);
		VecAdd(vTmp,vTmp,vHHere);
		NegateVec(vTmp,vTmp);
		Normalize(vTmp,vTmp);
		VecScale(vTmp,vTmp,BL_NH);
		VecAdd(vHHere,vTmp,vNHere);
		if ((AssignCoords(" H  ",pmgdHere,vHHere,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	if (resnum==clength) {
		/* place C-terminus */
		/* N and C co-ordinates need to be relative to residue CA for
		   this function */
		VecSub(vNHere,vNHere,vCAHere);
		VecSub(vCLast,vCLast,vCALast);
		if (FindCTerm(pmgdHere,vNHere,vCLast,vCAHere,Model)!=ERR_SUCCESS) return ERR_FAIL;
	}
	return ERR_SUCCESS;
}

/* given a Calpha with co-ordinates defined for it and Calpha i+/-1, this
   function sets the co-ordinates of the corresponding CBeta - uses first
   PALD in the list if more than one exist - and all other residue atoms 
   cis will be set to 1 when the current residue is cis-, 0 otherwise */
TrajErr SetAllCoords(PMGD pmgdHere,Int2 Model,FloatLo CBOffset,FloatLo Psi,FloatLo Omega)
{
        PMGD pmgdLast,pmgdNext;
        PDNMG pdnmgLast,pdnmgNext,pdnmgHere;
        vec vRMinus,vRPlus,uone,utwo,uthree,vCANext,vCBHereRef;
	vec vTmp,vCHere,vHALast,vOLast,vNNext;
	vec vTmp1,vTmp2,vTmp3,vTmp4,vTmp5;
	Int2 cistmp1,cistmp2,cistmp3;
	Char residueHere,residueLast,residueNext;
	Int2 resHere,resNext,resLast,resnum;
	Int4 clength;
	PMMD pmmdParent;

	/* use separate build routine for phi-psi walk */
	if (WALKTYPE==WALK_PHIPSI) {
		/* note Phi is stored in CBOffset in this case */
		return SetAllCoordsPhiPsi(pmgdHere,Model,CBOffset,Psi,Omega);
	}
	/* ca walk */
	pmmdParent=(PMMD)(pmgdHere->pfbParent);
	clength=pmmdParent->iResCount;
        pdnmgHere=pmgdHere->pdnmgLink;
        pdnmgLast=pdnmgHere->last;
        pdnmgNext=pdnmgHere->next;
	resnum=pdnmgHere->choice;
        if (pdnmgLast==NULL) {
                ErrPostEx(SEV_ERROR,1,1,"Unable to find C beta for first residue");
                return ERR_FAIL;
        }
        if (pdnmgNext==NULL) {
                ErrPostEx(SEV_ERROR,1,1,"Unable to find C beta for last residue");
                return ERR_FAIL;
        }
        /* residues must be ordered consecutively in the linked list */
        if ((pdnmgNext->choice-pdnmgHere->choice!=1) || (pdnmgHere->choice-pdnmgLast->choice!=1)) {
                ErrPostEx(SEV_ERROR,1,1,"Residues out of order");
                return ERR_FAIL;
        }
        pmgdLast=(PMGD)(pdnmgLast->data.ptrvalue);
        pmgdNext=(PMGD)(pdnmgNext->data.ptrvalue);
	residueHere=GetAAFromIDict(pmgdHere);
	residueLast=GetAAFromIDict(pmgdLast);
	residueNext=GetAAFromIDict(pmgdNext);
	resHere=GetResIdxFromMGD(NCBIstdaaUC,pmgdHere);
	resLast=GetResIdxFromMGD(NCBIstdaaUC,pmgdLast);
	resNext=GetResIdxFromMGD(NCBIstdaaUC,pmgdNext);
	/* take care of new cis-Proline residues */
	if ((residueNext=='P') && (cisNext))
		resNext=CISPRO;
	if ((residueHere=='P') && (cisHere))
		resHere=CISPRO;
	if ((residueLast=='P') && (cisLast))
		resLast=CISPRO;
	GetCoOrds(pmgdHere," CA ",vZero,vCAHere,Model);
	GetCoOrds(pmgdLast," CA ",vZero,vCALast,Model);
	GetCoOrds(pmgdNext," CA ",vZero,vCANext,Model);
        /* construct local axis system */
	SetUpRefAxes(uone,utwo,uthree,vRMinus,vRPlus,vCANext);
	/* get CB dir unit vector in vCBTrue, relative to CA */
	/* Ref indicates with respect to the reference U axes */
	FindCBDir(uone,utwo,uthree,vCANext,resHere,resNext,vCBHereRef,CBOffset);
	/* and assign it for use with Place Rotamer later */
	VecScale(vTmp,vCBHere,bl_cacb[resHere-1]);
	VecAdd(vTmp,vTmp,vCAHere);
	if (residueHere!='G') {
		if ((AssignCoords(" CB ",pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	else {
		if ((AssignCoords("2HA ",pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	if (resnum==2) {
		/* next find approx. N dir in vNHere */
		if (resHere==CISPRO) {
			if (FindAdjacentPeptideAtom(uone,utwo,uthree,vRMinus,vCBHereRef,ba_ncacb[resHere-1],BL_NCA+BLSD_NCA*Rand1Distrib(),cos(DEGTORAD*CIS_ZETA),vNHere)!=ERR_SUCCESS)
				return ERR_FAIL;
		}
		else {
			if (FindAdjacentPeptideAtom(uone,utwo,uthree,vRMinus,vCBHereRef,ba_ncacb[resHere-1],BL_NCA+BLSD_NCA*Rand1Distrib(),coszeta[resHere],vNHere)!=ERR_SUCCESS)
				return ERR_FAIL;
		}
		VecAdd(vTmp,vNHere,vCAHere);
		if ((AssignCoords(" N  ",pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		/* vNHere relative to CA */
		/* do carbonyl O of previous residue */
		/* vOLast not a unit vector, relative to CAi-1 */
		if (resHere==CISPRO)
			FindAtomInResMinusOne(" N  ",pmgdHere,vRMinus,-CIS_OCACA,CIS_OCA,Model,vOLast);
		else
			FindAtomInResMinusOne(" N  ",pmgdHere,vRMinus,ba_ocaca[resLast],bl_cao[resLast],Model,vOLast);
		VecAdd(vTmp,vOLast,vCALast);
		if ((AssignCoords(" O  ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		/* do carbonyl C of previous residue in vCLast, relative to CAi-1 */
		if (resHere==CISPRO)
			FindAtomInResMinusOne(" N  ",pmgdHere,vRMinus,-CIS_ETA,BL_CAC+BLSD_CAC*Rand1Distrib(),Model,vCLast);
		else
			FindAtomInResMinusOne(" N  ",pmgdHere,vRMinus,RADTODEG*acos(coseta[resLast]),BL_CAC+BLSD_CAC*Rand1Distrib(),Model,vCLast);
		VecAdd(vTmp,vCLast,vCALast);
		if ((AssignCoords(" C  ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
/*DumpVal((PFB)pmgdHere);
exit(0);*/
	}
	else {
		if (FindPeptideAtoms(pmgdLast,pmgdHere,Model,Omega)!=ERR_SUCCESS) {
			return ERR_FAIL;
		}
		VecSub(vCLast,vCLast,vCALast);
		VecSub(vNHere,vNHere,vCAHere);
	}
	/* place H if not proline */
	if (residueHere!='P') {
		/* vNHere set for carbonyl O above */
		FindHi(vTmp);
		if ((AssignCoords(" H  ",pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	/* place sidechain on residue i-1 */
	if (resnum>2) {
		if ((PlaceRotamer(pmgdLast,vHALast,Model))!=ERR_SUCCESS) return ERR_FAIL;
		/* place HA i-1 */
		VecAdd(vTmp,vHALast,vCALast);
		if (residueLast!='G') {
			if ((AssignCoords(" HA ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		}
		else {
			if ((AssignCoords("1HA ",pmgdLast,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		}
	}
	else {
		/* resnum == 2 */
		/* place N-terminus */
		if ((FindNTerm(pmgdLast,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	if (resnum==clength-1) {
		/* need 2nd last C first */
		/* vCB is ignored in this case */
		/* temporary vCALast=vCAHere, vCAHere=vCANext, vCBHere=vZero, vCLast=vCHere,
		   vNHere=vNNext, cisLast=cisHere, cisHere=cisNext, cisNext=0 */
		VecScale(vTmp1,vCALast,1.0);
		VecScale(vTmp2,vCAHere,1.0);
		VecScale(vTmp3,vCBHere,1.0);
		VecScale(vTmp4,vCLast,1.0);
		VecScale(vTmp5,vNHere,1.0);
		cistmp1=cisLast;
		cistmp2=cisHere;
		cistmp3=cisNext;
		VecScale(vCALast,vCAHere,1.0);
		VecScale(vCAHere,vCANext,1.0);
		VecScale(vCBHere,vZero,1.0);
		/* must correct omega temporarily too! - the 9.8 factor accounts */
		/* for the 2.94 and 0.3 scalings applied above to Omega */
		if (cisHere!=cisNext)
			Omega=(180.0-Omega)/9.8;
		cisLast=cisHere;
		cisHere=cisNext;
		cisNext=0;
		if ((FindPeptideAtoms(pmgdHere,pmgdNext,Model,Omega))!=ERR_SUCCESS) return ERR_FAIL;
		/* give vNNext, vCHere their real value */
		VecScale(vCHere,vCLast,1.0);
		VecScale(vNNext,vNHere,1.0);
		/* vCALast is saved value of vCAHere */
		VecSub(vCHere,vCHere,vCALast);
		VecSub(vNNext,vNNext,vCANext);
		/* update changes to temporary global vars */
		VecScale(vCLast,vCHere,1.0);
		VecScale(vNHere,vNNext,1.0);
		/* last H */
		if (residueNext!='P') {
			/* temporary vCALast=vCAHere, vCAHere=vCANext, vNHere=vNNext, vCLast=vCHere */
			FindHi(vTmp);
			if ((AssignCoords(" H  ",pmgdNext,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		}
		/* 2nd last sidechain */
		/* temporary vCALast=vCAHere, vCLast=vCHere */
		if ((PlaceRotamer(pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		/* restore saved global variables */
		VecScale(vCALast,vTmp1,1.0);
		VecScale(vCAHere,vTmp2,1.0);
		VecScale(vCBHere,vTmp3,1.0);
		VecScale(vCLast,vTmp4,1.0);
		VecScale(vNHere,vTmp5,1.0);
		cisLast=cistmp1;
		cisHere=cistmp2;
		cisNext=cistmp3;
		/* put Omega back */
		if (cisHere!=cisNext)
			Omega=180.0-9.8*Omega;
		/* 2nd last HA */
		Translate(vTmp,vCAHere);
		if (residueHere!='G') {
			if ((AssignCoords(" HA ",pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		}
		else {
			if ((AssignCoords("1HA ",pmgdHere,vTmp,Model))!=ERR_SUCCESS) return ERR_FAIL;
		}
		/* and finally C-terminal atoms */
		if ((FindCTerm(pmgdNext,vNNext,vCHere,vCANext,Model))!=ERR_SUCCESS) return ERR_FAIL;
	}
	return ERR_SUCCESS;
}

/* unset all those possibly set by SetAllCoordsPhiPsi */
void UndoBDPhiPsi(PMGD pmgdHere,Int2 Model)
{
	PVNMA pvnmaHere,pvnmaLast;
	vec vHere;
	PDNMG pdnmgHere;
	PMGD pmgdLast;
	CharPtr AtomName;
	PMMD pmmdParent;
	PMAD pmadHere;
	Int2 clength,resnum;

	pmmdParent=(PMMD)(pmgdHere->pfbParent);
	clength=(Int2)(pmmdParent->iResCount);
	pdnmgHere=pmgdHere->pdnmgLink;
  pmgdLast=(PMGD)((pdnmgHere->last)->data.ptrvalue);
	resnum=pdnmgHere->choice;
	pvnmaHere=pmgdHere->pvnmaAHead;
	pvnmaLast=pmgdLast->pvnmaAHead;
	while (pvnmaHere) {
		pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
		AtomName=pmadHere->pcAName;	
		if (((!StringCmp(AtomName," CA ")) || (!StringCmp(AtomName," N  ")) || (!StringCmp(AtomName," H  "))) || (resnum==clength)) {
			if (GetCoOrds(pmgdHere,AtomName,vZero,vHere,Model)!=NULL) {
				if (IsAtomPlaced(pmadHere,Model))
					BDRemove(vHere,&(pwsThis->pbdTreeHead),&phbsHBondsToCheck);
			}
		}
		pvnmaHere=pvnmaHere->next;
	}
	while (pvnmaLast) {
		pmadHere=(PMAD)(pvnmaLast->data.ptrvalue);
		AtomName=pmadHere->pcAName;	
		if ((((StringCmp(AtomName," CA ")) && (StringCmp(AtomName," N  ")) && (StringCmp(AtomName," H  ")) && !((resnum==2) && (!StringCmp(AtomName," C  "))))) || ((resnum==2) && (!StringCmp(AtomName," H  ")))) {
			if (GetCoOrds(pmgdLast,AtomName,vZero,vHere,Model)!=NULL) {
				if (IsAtomPlaced(pmadHere,Model))
					BDRemove(vHere,&(pwsThis->pbdTreeHead),&phbsHBondsToCheck);
			} 
		}
		pvnmaLast=pvnmaLast->next;
	}
}

/* unset all those possibly set by SetAllCoords */
void UndoBD(PMGD pmgdHere,Int2 Model)
{
	PVNMA pvnmaHere,pvnmaNext,pvnmaLast;
	vec vHere;
	PDNMG pdnmgHere;
	PMGD pmgdNext,pmgdLast;
	CharPtr AtomName;
	PMMD pmmdParent;
	PMAD pmadHere;
	Int2 clength,resnum;

	if (WALKTYPE==WALK_PHIPSI) {
		UndoBDPhiPsi(pmgdHere,Model);
		return;
	}
	/* else WALKTYPE==WALK_CA */
	pmmdParent=(PMMD)(pmgdHere->pfbParent);
	clength=(Int2)(pmmdParent->iResCount);
	pdnmgHere=pmgdHere->pdnmgLink;
  pmgdNext=(PMGD)((pdnmgHere->next)->data.ptrvalue);
  pmgdLast=(PMGD)((pdnmgHere->last)->data.ptrvalue);
	resnum=pdnmgHere->choice;
	pvnmaHere=pmgdHere->pvnmaAHead;
	pvnmaNext=pmgdNext->pvnmaAHead;
	pvnmaLast=pmgdLast->pvnmaAHead;
	while (pvnmaHere) {
		pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
		AtomName=pmadHere->pcAName;	
		if (((!StringCmp(AtomName," CB ")) || (!StringCmp(AtomName,"2HA ")) || (!StringCmp(AtomName," N  ")) || (!StringCmp(AtomName," H  "))) || ((resnum==clength-1) && (StringCmp(AtomName," CA ")))) {
			if (GetCoOrds(pmgdHere,AtomName,vZero,vHere,Model)!=NULL) {
				if (IsAtomPlaced(pmadHere,Model))
					BDRemove(vHere,&(pwsThis->pbdTreeHead),&phbsHBondsToCheck);
			}
		}
		pvnmaHere=pvnmaHere->next;
	}
	while (pvnmaNext) {
		pmadHere=(PMAD)(pvnmaNext->data.ptrvalue);
		AtomName=pmadHere->pcAName;	
		if ((!StringCmp(AtomName," CA ")) || (resnum==clength-1)) { 
			if (GetCoOrds(pmgdNext,AtomName,vZero,vHere,Model)!=NULL) {
				if (IsAtomPlaced(pmadHere,Model))
					BDRemove(vHere,&(pwsThis->pbdTreeHead),&phbsHBondsToCheck);
			} 
		}
		pvnmaNext=pvnmaNext->next;
	}
	while (pvnmaLast) {
		pmadHere=(PMAD)(pvnmaLast->data.ptrvalue);
		AtomName=pmadHere->pcAName;	
		if (((StringCmp(AtomName," N  ")) && (StringCmp(AtomName," H  ")) && (StringCmp(AtomName," CA ")) && (StringCmp(AtomName," CB "))) || ((resnum==2) && (StringCmp(AtomName," CA ")))) {
			if (GetCoOrds(pmgdLast,AtomName,vZero,vHere,Model)!=NULL) {
				if (IsAtomPlaced(pmadHere,Model))
					BDRemove(vHere,&(pwsThis->pbdTreeHead),&phbsHBondsToCheck);
			}
		}
		pvnmaLast=pvnmaLast->next;
	}
}

/* return FALSE if cis- non-Pro residue found */
Boolean CheckFragCis(PFDS pfdsThis,PDNMG pdnmgThis)
{
	PFDS pfdsHere;
	Int4 OmegaTemp,resHere;
	PDNMG pdnmgHere;
	PMGD pmgdHere;
	
	pfdsHere=pfdsThis;
	pdnmgHere=pdnmgThis;
	while (pfdsHere!=NULL) {
		OmegaTemp=((Int4)(pfdsHere->ChiWMean))%360;
		while (OmegaTemp<0)
			OmegaTemp+=360;			
		if (OmegaTemp<90 || OmegaTemp>270) {
			/* cis residue in fragment */
			if (WALKTYPE==WALK_CA)
				resHere=pfdsHere->resnum;
			else /* WALKTYPE==WALK_PHIPSI */
				resHere=pfdsHere->resnum+1;
			while (pdnmgHere->choice>resHere)
				pdnmgHere=pdnmgHere->last;
			while (pdnmgHere->choice<resHere) {
				pdnmgHere=pdnmgHere->next;
				if (pdnmgHere==NULL)
					return FALSE;
			}
			pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
			if (GetAAFromIDict(pmgdHere)!='P')
				return FALSE;
		}			
		pfdsHere=pfdsHere->next;
	}
	return TRUE;
}

float AngleDistortion(float angle1, float angle2, float std)
{
	float angledif;

	if (std < 0.5) std=0.5;  /* The smaller this value, the less likely tunneling will occur
							    omegas tend to have small stds */

	angledif=fabs(angle1-angle2);
	if ( angledif < (360.0-angledif) )
		return (angledif / std);
	else
		return ((360.0-angledif) / std);
}

/* 
The following function estimates the likelihood that a given fragment can grow ("tunnel")
out of another as-of-yet unfinished fragment.  It is based on congruency between
their respective angles scaled by their standard deviations.  The model is developed
in "A Model For Estimating Protein Fragment Tunneling Probability in Ab Initio Protein Folding" [Michael Brougham,
Howard Feldman. 2004, Unpublished as of July 2004]
*/
float ProbabilityLikelihood(PFDS pfdsResInFrag1, PFDS pfdsResInFrag2, Int4 TotalFragmentLength1)
{
	/* FILE* fout; */
	Int4 length1,length2,minlength;
	int i;
	FloatLo Phi1,Psi1,Phi2,Psi2,AngleSTD1,AngleSTD2,Omega1,Omega2,OmegaSTD1,OmegaSTD2;
	FloatLo ResEuclideanDist;
	FloatLo PhiDist, PsiDist, OmegaDist;
	FloatLo Prob1,Prob2;
	FloatLo averageAngleDistortion;
	FloatLo overallAngleDistortion, tunnelingProb;
	FloatLo exponent;
	Int4 numOverlappingResidues;

	FloatLo denominator;
	FloatLo numerator;

	overallAngleDistortion=0.0;
	averageAngleDistortion=0.0;

	tunnelingProb=0.0;

	/*	printf("\n\nCALCULATING PROBABILITY LIKELIHOOD"); */

	if ( (pfdsResInFrag1 != NULL) && (pfdsResInFrag2 != NULL) )
	{

	/* determine the length of the first of the two fragments */

		Prob1 = pfdsResInFrag1->prob;
		Prob2 = pfdsResInFrag2->prob;

	/* ascertain which of the PFDS is shorter (ie has fewer remaining residues) */
	
		length1=pfdsResInFrag1->length;
		length2=pfdsResInFrag2->length;

		/*pfdsResInFrag2->resnum*/
		
		if (length1 < length2) 
			 minlength = length1;
		else minlength = length2;

		/*
		fout=fopen("NEWCODETEST.txt","a+");
		fprintf(fout,"\n\n\n\n-----------------------------------------------");
		fprintf(fout,"\n%d %d",length1,length2);
		*/

		numOverlappingResidues = minlength + 1;
						
		/* for, that length, traverse the remaining residues and compare their angle values */
		
		for (i = 0; i < minlength + 1; i++)
		{
			/* fprintf(fout,"\n\n\ni = %d", i); */
			Phi1 = pfdsResInFrag1->Angle1;
			Psi1 = pfdsResInFrag1->Angle2;
			Omega1 = pfdsResInFrag1->ChiWMean;
			AngleSTD1 = pfdsResInFrag1->AngleSD;
			OmegaSTD1 = pfdsResInFrag1->ChiWSD;

			Phi2 = pfdsResInFrag2->Angle1;
			Psi2 = pfdsResInFrag2->Angle2;
			Omega2 = pfdsResInFrag2->ChiWMean;
			AngleSTD2 = pfdsResInFrag2->AngleSD;
			OmegaSTD2 = pfdsResInFrag2->ChiWSD;

			PhiDist = (FloatLo)AngleDistortion(Phi1,Phi2,AngleSTD2);
			PsiDist = (FloatLo)AngleDistortion(Psi1,Psi2,AngleSTD2);
			OmegaDist = (FloatLo)AngleDistortion(Omega1,Omega2,OmegaSTD2);

			overallAngleDistortion= (PhiDist + PsiDist + OmegaDist) / 3;
			averageAngleDistortion+=overallAngleDistortion;

			/*
			fprintf(fout,"\nPhi1 = %f Phi2 = %f \nPsi1 = %f Psi2 = %f",Phi1,Phi2,Psi1,Psi2);
			fprintf(fout,"\nOmega1 = %f Omega2 = %f", Omega1, Omega2);
			
			fprintf(fout, "\nAngleSTD1 = %f AngleSTD2 = %f OmegaSTD1 = %f OmegaSTD2 = %f",
					AngleSTD1, AngleSTD2, OmegaSTD1, OmegaSTD2);

			fprintf(fout,"\nPhiDist = %f",PhiDist);
			fprintf(fout,"\nPsiDist = %f",PsiDist);
			fprintf(fout,"\nOmegaDist = %f",OmegaDist);

			fprintf(fout,"\noverallAngleDistortion: %f",overallAngleDistortion);
			fprintf(fout,"\naverageAngleDistortion: %f",averageAngleDistortion);
			*/
			/* calculate the Euclidean distance scaled by the above length */
		
		pfdsResInFrag1 = pfdsResInFrag1->next;
		pfdsResInFrag2 = pfdsResInFrag2->next;
		}
		
		averageAngleDistortion=averageAngleDistortion / (numOverlappingResidues);

		numerator = exp(-averageAngleDistortion);

		exponent = 2.0 - ((double)TotalFragmentLength1 / (double)numOverlappingResidues);

		denominator = pow(Prob1, exponent);

		tunnelingProb = numerator * Prob2 / denominator;

		/*
		fprintf(fout,"\naverageAngleDistortion: %f",averageAngleDistortion);
		fprintf(fout,"\nnumerator: %10.2ef",numerator);
		fprintf(fout,"\nTotalFragmentLength1: %f",TotalFragmentLength1);
		fprintf(fout,"\nProb1: %10.2ef", Prob1);
		fprintf(fout,"\nProb2: %10.2ef", Prob2);
		fprintf(fout,"\nTotalFragmentLength1: %d", TotalFragmentLength1);
		fprintf(fout,"\nnumOverlappingResidues: %d", numOverlappingResidues);
		fprintf(fout,"\nexponent: ", exponent);
		fprintf(fout,"\ndenominator: %10.2ef", denominator);
		fprintf(fout,"\ntunnelingProb: %10.2ef", tunnelingProb);
		*/

		/* fclose(fout); */
	}

	/* printf("\nTunneling Likelihood= %10.2ef", tunnelingProb);  */
	return tunnelingProb;
}

/* main recursive procedure for building the protein, one residue at a time */
TrajErr NextResidue(PMGD pmgdHere,vec vLastLast,vec vLast,vec vHere,FloatLo PrevCosPhi,FloatLo PrevTheta,Int2 Model)
{
	PDNMG pdnmgHere;
	PMGD pmgdNext;
	PMMD pmmdHere;
	PFDS pfdsThis,pfdsLast;
	Char tmpbuf[PATH_MAX];
	vec vNext,vTmp1,vTmp2;
	FloatLo bl,cosPhi=0.0,Theta=0.0,arandnum,LastcosPhi,Phi=0.0,Psi=0.0,Omega=0.0,PrevPhi,CurPhi;
	FloatLo PhiDist,ThetaDist,S2Dist=0.0,markovprob=0.0,fct,nfactorial,fragprob,curfragprob,tunnelprob;
	Int4 valtest,OmegaTemp;
	Int2 timeout=0,lr,tmpflag,proceedflag=0,cnt,fudge,fragnum;
	TrajErr err,CAisOK;
	Int4 boundhi,boundlo,boundtest;
	static Boolean inafrag=FALSE;
	static PFDS pfdsHere=NULL,pfdsHead=NULL;
	static FloatLo nlogn=-1.0;
	static FloatLo rtlrx=1.0;

	pmmdHere=GetParentMol((PFB)pmgdHere);
	/* escape the recursion with success */
	if (pmgdHere==NULL && WALKTYPE==WALK_PHIPSI) {
		CLEARFRAGS
		return ERR_SUCCESS;
	}
	pdnmgHere=pmgdHere->pdnmgLink;
	pdnmgHere=pdnmgHere->next;
	/* number of retries before backing up further */
	/* escape the recursion with success */
	if (pdnmgHere==NULL && WALKTYPE==WALK_CA) {
		CLEARFRAGS
		return ERR_SUCCESS;
	}
/*	tm2=GetSecs();*/

	if (nlogn<0.0)
		nlogn=(FloatLo)(pmmdHere->iResCount)*LOG10((FloatLo)(pmmdHere->iResCount));
	/* times out after 5000 + 10 N log N * gen ^1/3 * res ^1/2 tries */
	/* where gen=1 for foldtraj, generation for DFP, and res is current residue number */
	if ((FloatLo)tries>(5000.0+2.0*nlogn*gencbrt*rtlrx)) {
		CLEARFRAGS
		return ERR_INCOMPLETE;
	}
	if (ProgramProgressMax==-1) {
		CLEARFRAGS
		return ERR_CANCELLED;
	}

	if (pdnmgHere==NULL)
		/* last residue in PHIPSI walk */
		lr=(pmmdHere->iResCount)+1;
	else
		lr=pdnmgHere->choice;

		if (lr!=lrx) {
			lrx=lr;
			rtlrx=sqrt((double)lrx);
			if (traj_quiet==VERBOSITY_QUIET || traj_quiet==VERBOSITY_VERBOSE) {
				if (WALKTYPE==WALK_CA)
					printf("\b\b\b\b\b\b\b\b\b%4d/%-4ld",lrx,(long int)pmmdHere->iResCount);
				else
					printf("\b\b\b\b\b\b\b\b\b%4d/%-4ld",lrx-1,(long int)pmmdHere->iResCount);
				fflush(stdout);
			}
			/* check if file exists before writing */
			if (WALKTYPE==WALK_CA)
				ProgramProgress=lrx;
			else
				ProgramProgress=lrx-1;

		}
	if (pdnmgHere!=NULL)
		pmgdNext=(PMGD)(pdnmgHere->data.ptrvalue);
	else
		pmgdNext=NULL;
	tmpflag=0;
	err=ERR_SUCCESS;
	do {
		if (inafrag) {
			if (tmpflag==0)
				pfdsHere=pfdsHere->next;
			if (pfdsHere==NULL) {
				/* end of fragment reached, presumably */
				inafrag=FALSE;
				/* free it */			
				pfdsHead=FreeFragmentList(pfdsHead);
			}
		}
		/* read trajectory graph from database */
		/* could make a bit more efficient by not reading if in a frag... */
		fudge=0;
		/* phi-psi walk is off by one residue relative to CA walk in terms of trajectory distribution-residue
		   number correspondence */
		if (WALKTYPE==WALK_PHIPSI)
			fudge=-1;
		ptgsHere=TrajGraphRead(lr-1+fudge);
/*printf("read %d to assign N of residue %d\n",lr-2,pmgdHere->pdnmgLink->choice);*/
	/* integral stored directly now */
	/*	TrajGraphIntegrate(ptgsHere);*/
		
		/* attempt to place next residue, until success or timeout */
		/* now check if any fragments exist that we want to follow */
		/* before start timeout loop */
		/* tmpflag indicates first time through do..while */

		/* We now want to check for fragments at every residue */
		/* We will either enter a new fragment or tunnel */
		if (tmpflag==0) {

			fragprob=fabs(Rand1());
			curfragprob=0.0;
			tunnelprob=fabs(Rand1());
			pfdsThis=ptgsHere->pfdsFragmentHead;	/* Fragment we might tunnel to/start */
			
			fragnum=0;
			/* Iterate through all the fragments on the trajectory graph */
			while (pfdsThis!=NULL) {
				curfragprob+=pfdsThis->prob;
				
				/* We have selected a fragment starting at the current residue */
				/* CheckFragCis checks if the fragment has any cis- residues aligned to non-Prolines
				   if so, the fragment will be skipped regardless of probability etc */
				if (curfragprob>fragprob && CheckFragCis(pfdsThis,pdnmgHere)) {

					/*ErrPostEx(SEV_ERROR,1,1,"Fragment with probability %f selected for use", curfragprob);
					ErrPostEx(SEV_ERROR,1,1,"Fragment of length %u starting at residue %u", pfdsThis->length, pfdsThis->resnum);*/
					
					/* Check if we actually tunnel */
					if(inafrag==TRUE) {
						/* tunnel length 0 automatically for CA walk only */
						/*if(pfdsHere->length==0 && WALKTYPE==WALK_CA || (pfdsHere->length > 0 && tunnelprob < (TUNNEL_PROB/(FloatLo)pfdsHere->length))){ */
						if(pfdsHere->length==0 && WALKTYPE==WALK_CA || (pfdsHere->length > 0 && tunnelprob < TUNNEL_PROB*ProbabilityLikelihood(pfdsHere, pfdsThis,pfdsHead->length + 1))){
							/*ErrPostEx(SEV_ERROR,1,1,"Tunneling...");
							ErrPostEx(SEV_ERROR,1,1,"Remaining length is: %u", pfdsHere->length);
							if (pfdsHere->length>0)
							ErrPostEx(SEV_ERROR,1,1,"Current tunnel probability is %f", (TUNNEL_PROB/(FloatLo)pfdsHere->length));*/
							/* Since we are already in a fragment, free that memory */
							CLEARFRAGS;
						}	/* If we decided not to tunnel, break out of the loop */
						else
							break;
					}
					/* At this point we either enter a new fragment or tunnel */
					inafrag=TRUE;
/* printf("adding fragment %d from res %d (%d-%d)\n",fragnum,lr-1+fudge,pfdsThis->resnum,pfdsThis->resnum+pfdsThis->length);					*/
					/* attempt to follow this fragment path */
					/* make a local copy of fragment path */
					ressrcres[lr-2+fudge]=lr-1+fudge;
					ressrcfrag[lr-2+fudge]=fragnum;
					pfdsLast=NULL;
					while (pfdsThis!=NULL) {
						pfdsHere=MemNew(sizeof(FDS));
						MemCopy(pfdsHere,pfdsThis,sizeof(FDS));
						pfdsHere->prev=NULL;
						pfdsHere->next=NULL;
						pfdsHere->nextfrag=NULL;
						if (pfdsHead==NULL) {
							pfdsHead=pfdsHere;
						}
						else {
							pfdsLast->next=pfdsHere;
							pfdsHere->prev=pfdsLast;
							ressrcres[pfdsHere->resnum-1]=ressrcres[pfdsHere->resnum-2];
							ressrcfrag[pfdsHere->resnum-1]=ressrcfrag[pfdsHere->resnum-2];
						}
						pfdsLast=pfdsHere;
						pfdsThis=pfdsThis->next;
					}
					pfdsHere=pfdsHead;	/* Get ahold of the fragment head */
					break;
				}
				pfdsThis=pfdsThis->nextfrag;
				fragnum++;
			}

			/* At this point, if inafrag==FALSE,
			then no fragment chosen */
		}
		/* now if in a frag, get necessary values */
		if (inafrag) {				
			TIMEOUT=pfdsHere->tout;
			if (tmpflag==0) {
				tmpflag=1;
				timeout=pfdsHere->tout;
			}
			/* integral stored directly now */
		/*	TrajGraphIntegrate(ptgsHere);*/
			/* attempt to place next residue, until success or timeout */
			/* executed on the order of 5-10 times per residue, so this is a
			   good time to decide whether bond will be cis- or trans- */
			IsCis[lr]=0;
			OmegaTemp=((Int4)(pfdsHere->ChiWMean))%360;
			while (OmegaTemp<0)
				OmegaTemp+=360;			
			if (WALKTYPE==WALK_CA) {
				if (OmegaTemp<90 || OmegaTemp>270)
					IsCis[lr]=1;
				else
					IsCis[lr]=0;
			}
			else if (WALKTYPE==WALK_PHIPSI) {
				if (OmegaTemp<90 || OmegaTemp>270)
					IsCis[lr-1]=1;
				else
					IsCis[lr-1]=0;
			}
			if (WALKTYPE==WALK_CA) {
				rotidHere=rotidNext;
				rotidNext=pfdsHere->rotid;
			}
			else if (WALKTYPE==WALK_PHIPSI)
				rotidHere=pfdsHere->rotid;
		}		
		else {	/* not in a frag */
			ressrcres[lr-2+fudge]=0;
			ressrcfrag[lr-2+fudge]=0;

			/* first time function is called, get timeout from record
			   but since function is recursive, don't reset each time
			   we get here since timeout is local to the recursive
		  	   function while TIMEOUT is global */
			TIMEOUT=ptgsHere->tout;
			if (tmpflag==0) {
				tmpflag=1;
				timeout=ptgsHere->tout;
			}
			/* executed on the order of 5-10 times per residue, so this is a
		  	   good time to decide whether bond will be cis- or trans- */
			IsCis[lr]=0;
			if (WALKTYPE==WALK_CA && pmgdNext!=NULL) {
				if (GetAAFromIDict(pmgdNext)=='P') {
					if (fabs(Rand1())<(ptgsHere->pCis)) {
						IsCis[lr]=1;
					}
				}
			}
			else if (WALKTYPE==WALK_PHIPSI) {
				if (GetAAFromIDict(pmgdHere)=='P') {
					if (fabs(Rand1())<(ptgsHere->pCis)) {
						IsCis[lr-1]=1;
					}
					else
						IsCis[lr-1]=0;
				}
			}
			if (WALKTYPE==WALK_CA) {
				rotidHere=rotidNext;
				rotidNext=ptgsHere->rotid;
			}
			else if (WALKTYPE==WALK_PHIPSI)
				rotidHere=ptgsHere->rotid;
		}
		do {
			if (err==ERR_FAIL) UndoBD(pmgdHere,Model);
			timeout--;
			/* insert trajectories here */
/*printf("phi\tpsi\n");
for (cnt=0;cnt<10000;cnt++) {*/	
			/* skip getting phi1, psi1 since not used */
			if ((lr!=3) || (WALKTYPE!=WALK_PHIPSI)) {
				/* now choose phi,psi or phi,theta */
				if (inafrag) {
					if (WALKTYPE==WALK_PHIPSI) {
						Phi=pfdsHere->Angle1+Rand1Distrib()*(pfdsHere->AngleSD);
						Psi=pfdsHere->Angle2+Rand1Distrib()*(pfdsHere->AngleSD);
					}
					else if (WALKTYPE==WALK_CA) {
						/* note: please ensure Phi, not cos Phi */
						cosPhi=cos(DEGTORAD*(pfdsHere->Angle1+Rand1Distrib()*(pfdsHere->AngleSD)));
						Theta=pfdsHere->Angle2+Rand1Distrib()*(pfdsHere->AngleSD);												
					}
				}
				else do {
					if ((IsCis[lr-1] && WALKTYPE==WALK_CA) || (IsCis[lr-2] && WALKTYPE==WALK_PHIPSI))
						arandnum=fabs(Rand1())*(FloatLo)(ptgsHere->CisTrajIntegral);
					else
						arandnum=fabs(Rand1())*(FloatLo)(ptgsHere->TrajIntegral);
					/* do binary search to find right place in graph */
					boundlo=0;
					boundhi=(ptgsHere->dim)*(ptgsHere->dim)-1;
/*					if ((IsCis[lr-1] && WALKTYPE==WALK_CA) || (IsCis[lr-2] && WALKTYPE==WALK_PHIPSI))
						valhi=ptgsHere->CisTrajIntegral;
					else
						valhi=ptgsHere->TrajIntegral;*/
					while (boundhi-boundlo>1) {
						boundtest=(boundhi+boundlo)/2;
						if ((IsCis[lr-1] && WALKTYPE==WALK_CA) || (IsCis[lr-2] && WALKTYPE==WALK_PHIPSI))
							valtest=(ptgsHere->CisTrajGraph)[boundtest];
						else
							valtest=(ptgsHere->TrajGraph)[boundtest];
						if ((FloatLo)valtest>arandnum) {
							boundhi=boundtest;
						}
						else {
							boundlo=boundtest;
						}
					}
					/* horizontal & vertical axis convention for trajectory
					   graph defined here */
					/* uses boundlo, i.e. rounds down */
					/* add random factor to remove discretization
					   effect */
					/* perhaps use different trajectory graph for cis-Pro?? */
					proceedflag=1;
					if (WALKTYPE==WALK_CA) {
						cosPhi=(((FloatLo)(boundlo/(ptgsHere->dim))+0.5)/((FloatLo)(ptgsHere->dim)))*2.0-1.0;
						Theta=(((FloatLo)(boundlo%(ptgsHere->dim))+0.5)/((FloatLo)(ptgsHere->dim)))*360.0-180.0;
						/* skip first time and anywhere near Cis */
						if ((ptgsHere->markovsf) && PrevTheta!=999.0 && !(IsCis[lr]) && !(IsCis[lr-1])) {
							PrevPhi=acos(PrevCosPhi);
							CurPhi=acos(cosPhi);
							PhiDist=PrevPhi-CurPhi;
							ThetaDist=PrevTheta-Theta;
							if (fabs(ThetaDist)>180.0)
								ThetaDist=360.0-fabs(ThetaDist);
							S2Dist=fabs(RADTODEG*acos(cos(PhiDist)+sin(CurPhi)*sin(PrevPhi)*(cos(DEGTORAD*ThetaDist)-1.0)));
							/* ensure S2Dist in [0,180] range */
							while (S2Dist>360.0)
								S2Dist-=360.0;
							if (S2Dist>180.0)
								S2Dist=360.0-S2Dist;
							if (PrevTheta>-50.0 && PrevTheta<110.0) {
								/* helical area of space */
/*	
markovprob=6.0*(0.2572+S2Dist*(-0.01394+S2Dist*(0.0002991+S2Dist*(-0.000003053+S2Dist*(0.00000001492-S2Dist*0.00000000002801)))));
*/
								/* new updated helical formula = polynomial */
								markovprob=0.447814+S2Dist*(-0.0190036+S2Dist*(0.000362334+S2Dist*(-0.00000317002+S2Dist*(0.0000000127676-S2Dist*0.0000000000194259))));
								/* + poisson curve which is negligible beyond S2Dist=40 degrees*/
								if (S2Dist<40.0) {
									fct=floor(S2Dist/2.5);
									/* linear approximation to factorial for any float number */
									nfactorial=(FloatLo)(factorial[(int)fct+1]-factorial[(int)fct])*(S2Dist/2.5-fct)+factorial[(int)fct];
									markovprob=markovprob+0.128556*pow(3.15,S2Dist/2.5)/nfactorial;
								}
							}
							else
								/* new updated sheet formula */							
/*	
markovprob=6.0*(-0.004682+S2Dist*(0.005227+S2Dist*(-0.0001569+S2Dist*(0.000001836+S2Dist*(-0.000000009397+S2Dist*0.00000000001752)))));
*/
							markovprob=0.0120254+S2Dist*(0.113435+S2Dist*(-0.00456418+S2Dist*(0.0000762228+S2Dist*(-0.000000630383+S2Dist*(0.00000000255312-S2Dist*0.0000000000040486)))));
							/* prevent potential problems with negative/small probabilities */
							if (markovprob<0.001) markovprob=0.001;
							if (markovprob>1.0) markovprob=1.0;
							if (fabs(Rand1())>(1.0-markovprob)*(1.0-ptgsHere->markovsf)+markovprob)
/*if
(fabs(Rand1())>(1.0-S2Dist/360.0)*(1.0-ptgsHere->markovsf)+markovprob)
*/
								proceedflag=0;
						}
					}
					else {	
						Psi=(((FloatLo)(boundlo/(ptgsHere->dim))+0.5)/((FloatLo)(ptgsHere->dim)))*360.0-180.0;
						Phi=(((FloatLo)(boundlo%(ptgsHere->dim))+0.5)/((FloatLo)(ptgsHere->dim)))*360.0-180.0;

/*printf("%d\t%f\t%f\n",lr-1,Phi,Psi);*/
					}
				} while (!proceedflag);
			}
			tries++;
			/* choose omega */
			if (inafrag) {
				Omega=pfdsHere->ChiWMean;
				if (WALKTYPE==WALK_CA) {
					/* actually inverted improper dihedral,
					   not omega -- has more variation */
					if (IsCis[lr-1])
						Omega=Omega+0.30*Rand1Distrib()*(pfdsHere->ChiWSD);
					else
						Omega=Omega+2.94*Rand1Distrib()*(pfdsHere->ChiWSD);
				}
				else
					Omega=Omega+Rand1Distrib()*(pfdsHere->ChiWSD);
			}
			else { /* not in a frag */
				Omega=ptgsHere->ChiWMean;
				if (IsCis[lr-1]) {
					if ((fabs(ptgsHere->ChiWMean)>90.0) && (ptgsHere->ChiWSD>0.0))
						Omega=180.0-(ptgsHere->ChiWMean);
				}
				else {
					if ((fabs(ptgsHere->ChiWMean)<=90.0) && (ptgsHere->ChiWSD>0.0))
						Omega=180.0-(ptgsHere->ChiWMean);
				}
				if (WALKTYPE==WALK_CA) {
					/* actually inverted improper dihedral,
					   not omega -- has more variation */
					if (IsCis[lr-1])
						Omega=Omega+0.30*Rand1Distrib()*(ptgsHere->ChiWSD);
					else
						Omega=Omega+2.94*Rand1Distrib()*(ptgsHere->ChiWSD);
				}
				else
					Omega=Omega+Rand1Distrib()*(ptgsHere->ChiWSD);
			}			
			if (WALKTYPE==WALK_CA) {				
				/* compute co-ords of next C-alpha */
				if (!IsCis[lr]) {
					bl=BL_CACA+Rand1Distrib()*BLSD_CACA;
				}
				else
					bl=BL_CACA_CISP+Rand1Distrib()*BLSD_CACA_CISP;
				/* check if just coming out of a special distance constraint and if so,
					 use specified angles to place the CA */
				CAisOK=CheckEndingDistConstraints(pmgdNext,vNext,Model,&dcviol);					
				if (CAisOK==ERR_SUCCESS) {
					CalcNextCoOrd(vLastLast,vLast,vHere,vNext,bl,180-(RADTODEG*acos(cosPhi)),Theta);
					CAisOK=CheckDistConstraints(pmgdNext,vNext,Model,&dcviol);
				}
				if (CAisOK==ERR_SUCCESS)
					CAisOK=AssignCoords(" CA ",pmgdNext,vNext,Model);
				else
					dcbad++;
				if (CAisOK!=ERR_SUCCESS)
					err=ERR_FAIL;
				else {
					/* compute co-ords of residue */
					cisLast=IsCis[lr-2];
					cisHere=IsCis[lr-1];
					cisNext=IsCis[lr];
					VecSub(vTmp1,vLastLast,vLast);
					VecSub(vTmp2,vHere,vLast);
					LastcosPhi=-(Dot(vTmp1,vTmp2)/(getMag(vTmp1)*getMag(vTmp2)));
					/* detect turns and place Cb differently in these cases */
					if ((cosPhi<0.20) && ((Theta>=-150.0) && (Theta<-25.0)) && (LastcosPhi<0.20))
						err=SetAllCoords(pmgdHere,Model,-0.45,0.0,Omega);
					else
						err=SetAllCoords(pmgdHere,Model,0.0,0.0,Omega);
				}
			}
			else { /* WALK_PHIPSI */
				/* compute co-ords of residue */
				cisLast=IsCis[lr-2];
				cisHere=IsCis[lr-1];
				cisNext=IsCis[lr];
				err=SetAllCoords(pmgdHere,Model,Phi,Psi,Omega);
			}
/*			if (phbsHBondsToCheck!=NULL) {*/
				/* unsatisfied H-bonds, must go back */
/*				phbsHBondsToCheck=PHBSFree(phbsHBondsToCheck);
			}                                                 */
			/* backtrack a bit earlier than before, probabilistically */
			if (inafrag) {
				if (fabs(Rand1())-0.75>((FloatLo)timeout)/((FloatLo)(pfdsHere->tout)))
					timeout=0;
			}
			else if (fabs(Rand1())-0.75>((FloatLo)timeout)/((FloatLo)(ptgsHere->tout)))
				timeout=0;		
		} while ((timeout>0) && (err==ERR_FAIL));
		if (err==ERR_GIVEUPRES) timeout=0;
		ptgsHere=FreeTraj(ptgsHere);
		if (!timeout) {
			/* backtracking here */
			UndoBD(pmgdHere,Model);
			/* check for backing out of a fragment */
			if (inafrag) {
				pfdsHere=pfdsHere->prev;
				if (pfdsHere==NULL) {
					/* start of fragment reached, presumably */
					inafrag=FALSE;
					/* free it */			
					pfdsHead=FreeFragmentList(pfdsHead);
					/*if (WALKTYPE==WALK_CA) {	*/
						/* current residue is lr-1, and subtract one for array index */
						if (lr+fudge>2 && ressrcres[lr-3+fudge]>0) {
							/* previous residue was in a fragment */
							inafrag=TRUE;
							/* load up the fragment */
							ptgsHere=TrajGraphRead(ressrcres[lr-3+fudge]);
							pfdsThis=ptgsHere->pfdsFragmentHead;	/* Fragment we might tunnel to/start */
							/* Iterate through all the fragments on the trajectory graph */
							cnt=ressrcfrag[lr-3+fudge];
							while (cnt) {
								pfdsThis=pfdsThis->nextfrag;
								cnt--;
							}
/* printf("backtracking to fragment %d from res %d (%d-%d)\n",ressrcfrag[lr-3+fudge],ressrcres[lr-3+fudge],pfdsThis->resnum,pfdsThis->resnum+pfdsThis->length); */
							pfdsLast=NULL;
							while (pfdsThis!=NULL) {
								pfdsHere=MemNew(sizeof(FDS));
								MemCopy(pfdsHere,pfdsThis,sizeof(FDS));
								pfdsHere->prev=NULL;
								pfdsHere->next=NULL;
								pfdsHere->nextfrag=NULL;
								if (pfdsHead==NULL) {
									pfdsHead=pfdsHere;
								}
								else {
									pfdsLast->next=pfdsHere;
									pfdsHere->prev=pfdsLast;
								/*	ressrcres[pfdsHere->resnum-1]=ressrcres[pfdsHere->resnum-2];
									ressrcfrag[pfdsHere->resnum-1]=ressrcfrag[pfdsHere->resnum-2];*/
								}
								pfdsLast=pfdsHere;
								pfdsThis=pfdsThis->next;
							}
							ptgsHere=FreeTraj(ptgsHere);
							pfdsHere=pfdsHead;	/* Get ahold of the fragment head */
							/* go to previous residue */
							while (pfdsHere->resnum<lr-2+fudge)
								pfdsHere=pfdsHere->next;
							/* add back fragment source after tunnelling */
							ressrcres[lr-2+fudge]=ressrcres[lr-3+fudge];
							ressrcfrag[lr-2+fudge]=ressrcfrag[lr-3+fudge];
							cnt=pfdsHere->length;
							while (cnt) {
								ressrcres[lr-2+fudge+cnt]=ressrcres[lr-3+fudge];
								ressrcfrag[lr-2+fudge+cnt]=ressrcfrag[lr-3+fudge];
								cnt--;
							}
						}
/*else
printf("backtracking to no-fragment at res %d\n",lr-2+fudge);
*/
					/*}*/
				}
			}			
			return ERR_FAIL;
		}
		err=NextResidue(pmgdNext,vLast,vHere,vNext,cosPhi,Theta,Model);
	} while (err==ERR_FAIL);
	return err;
}

/* reads in table of Cbeta locations */
TrajErr LoadCBTable(void)
{
	FILE *f;
        Char buf[255];
	Int2 res=0;
	Char rname[255];
        float  cbtemp [NUMBINS+1][3];
        float  costemptze;
        float  costempet;
        float  bl_caotemp;
        float  ba_ocacatemp;
        float  cbcisdirtemp[3];

		sprintf(buf,"%s%s",CFG_local_datafilepath,CBFNAME);
		if ((CheckMD5(buf,"68753e20539e235f44a41bdcdcd62d04")==FALSE) && (CheckMD5(buf,"cac5372a34fc5255ffb0fc626bfdfe15")==FALSE)) {
			ErrPostEx(SEV_ERROR,1,1,"%s failed checksum",CBFNAME);
			return ERR_FAIL;
		}
		if ((f=FileOpen(buf,"r"))==NULL) {
                return ERR_FAIL;
        }
        while (FileGets(buf,255,f)!=NULL) {
		sscanf(buf,"%s",rname);
		if (rname[0]=='#') continue;
		/* place these after IUPAC definition */
		if (!StringCmp(rname,"CX")) res=CYSSYC;
		else if (!StringCmp(rname,"PC")) res=CISPRO;
		else
			if (StringLen(rname)>1) continue;
		if (StringLen(rname)==1)
	 		res=StringCSpn(NCBIstdaaUC,rname);
               sscanf(buf,"%*s %*s %f %f %f %f %f %f %f %f %f %f %f %f %f%f %f %f %f %f %f %f %f %f %f %f %f",
		&cbtemp[0][0],&cbtemp[1][0],&cbtemp[2][0],&cbtemp[3][0],&cbtemp[4][0],&cbtemp[5][0],&cbtemp[0][1],
		&cbtemp[1][1],&cbtemp[2][1],&cbtemp[3][1],&cbtemp[4][1],&cbtemp[5][1],&cbtemp[0][2],
		&cbtemp[1][2],&cbtemp[2][2],&cbtemp[3][2],&cbtemp[4][2],&cbtemp[5][2],
                &costemptze,&costempet,&bl_caotemp,&ba_ocacatemp,   
                &cbcisdirtemp[0],&cbcisdirtemp[1],&cbcisdirtemp[2]);
               cbdir[res][0][0] = (FloatLo) cbtemp[0][0];  
               cbdir[res][1][0] = (FloatLo) cbtemp[1][0];
               cbdir[res][2][0] = (FloatLo) cbtemp[2][0];
               cbdir[res][3][0] = (FloatLo) cbtemp[3][0];
               cbdir[res][4][0] = (FloatLo) cbtemp[4][0];
               cbdir[res][5][0] = (FloatLo) cbtemp[5][0];
               cbdir[res][0][1] = (FloatLo) cbtemp[0][1];
               cbdir[res][1][1] = (FloatLo) cbtemp[1][1];
               cbdir[res][2][1] = (FloatLo) cbtemp[2][1];
               cbdir[res][3][1] = (FloatLo) cbtemp[3][1];
               cbdir[res][4][1] = (FloatLo) cbtemp[4][1];
               cbdir[res][5][1] = (FloatLo) cbtemp[5][1];
               cbdir[res][0][2] = (FloatLo) cbtemp[0][2];   
               cbdir[res][1][2] = (FloatLo) cbtemp[1][2];   
               cbdir[res][2][2] = (FloatLo) cbtemp[2][2];   
               cbdir[res][3][2] = (FloatLo) cbtemp[3][2];
               cbdir[res][4][2] = (FloatLo) cbtemp[4][2];
               cbdir[res][5][2] = (FloatLo) cbtemp[5][2];
               coszeta[res]     = (FloatLo) costemptze;
               coseta[res]      = (FloatLo) costempet;
               bl_cao[res]      = (FloatLo) bl_caotemp;
               ba_ocaca[res]    = (FloatLo) ba_ocacatemp;
               cbcisdir[res][0] = (FloatLo) cbcisdirtemp[0];
               cbcisdir[res][1] = (FloatLo) cbcisdirtemp[1];
               cbcisdir[res][2] = (FloatLo) cbcisdirtemp[2];
	}
        FileClose(f);
	return ERR_SUCCESS;
}

/* main calling routine to fold the protein, after chemical graph has been
   instantiated  */
TrajErr TRADEProteinEx(PMSD pmsdRoot,Int2 Model, Int4 StructureNumber)
{
	PDNMG pdnmgLast,pdnmgHere;
	PMGD pmgdLast,pmgdHere;
	vec vTmp,vTmp2,vTmp3,vCANTerm,vCACTerm,vCATerm;
	PMMD pmmdParent;
/*	PMAD pmadHere,pmad2;*/
	ValNodePtr /*vnpHead,*/vnpHere,vnpChis;
/*	PHBS phbsHere;*/
	Int2 clength,efflength;
	TrajErr err;
	FloatLo Rgyr,RgyrHP,Rn,Cn;
	Int2 cnt;
	PRS prsHere,prsHead;
	pnewrotlibrecord prlr;
	FloatLo bl,ba,lastPhi,firstPsi,firstOmega,OmegaHere;
	vec v1,v2,v3,v4;
	FILE *fp;
#ifdef USE_DSSP
	Char cResHere;
	Char hpresidues[]="MILVYCFGAW";
	Int4 *surfacc;
  FloatHi Surface,HPSurface;
  Char pcDSSPAssign[MAXRES+1];
	Int4 NumDSSPHelix,NumDSSPSheet;
#endif

	fpabad=0;
	crashcnt=0;
	dcbad=0;
	dcviol=0;

	ProgramProgress=2;
	ProgramProgressMax=((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue))->iResCount;

	pwsThis=NULL;
	if (traj_quiet==VERBOSITY_VERBOSE)
	  	printf("Initializing distance constraint list..\n");
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Initializing BD-tree and worldstrucs..\n");
        pwsThis=AddtoWorld(pwsThis,Model,(PFB)pmsdRoot);
        /* instantiate the world */
/*        vnpHead=*/InstantiateWorld(1,pwsThis);
	FreeBDTree(pwsThis->pbdTreeHead);
	pwsThis->pbdTreeHead=NULL;
	pdnmgLast=((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue))->pdnmgHead;
	pdnmgHere=pdnmgLast->next;
	pmgdLast=(PMGD)(pdnmgLast->data.ptrvalue);
	pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
	pmmdParent=(PMMD)(pmgdHere->pfbParent);

	lastresdrawn=0;
	/* atoms will be marked as atoms are placed */
	ClearUpdate(pmsdRoot->pdnmmHead,Model);
	tries=2;
	clength=pmmdParent->iResCount;
	ressrcres=MemNew(sizeof(VoidPtr)*clength);
	ressrcfrag=MemNew(sizeof(VoidPtr)*clength);
	vTmp3[0]=0;
	vTmp3[1]=0;
	vTmp3[2]=1;
	vTmp[0]=0;	
	vTmp[1]=0;
	/* so WHATIF doesn't cry */	
	vTmp[2]=0.001;
	if (AssignCoords(" CA ",pmgdLast,vTmp,Model)!=ERR_SUCCESS) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,1,1,"Unable to place first atom! (This should never happen)");
	}
	IsCis[0]=0;
	IsCis[1]=0;
	IsCis[2]=0;
	if (traj_quiet==VERBOSITY_VERBOSE)
		printf("Loading CB table..\n");
	if (LoadCBTable()!=ERR_SUCCESS) {
		PurgeGlobs();
		ErrPostEx(SEV_FATAL,99,99,"Unable to load B-carbon table %s%s.",CFG_local_datafilepath,CBFNAME);
	}
        /* initialize CodeBase Database */
	if (WALKTYPE==WALK_CA) {
		/* check if the first peptide bond should be cis */
		ptgsHere=TrajGraphRead(1);
		if (GetAAFromIDict(pmgdHere)=='P') {
			if (fabs(Rand1())<(ptgsHere->pCis)) {
				IsCis[2]=1;
			}
		}
		vTmp2[1]=0;	
		vTmp2[2]=0;
		do {
			if (IsCis[2])
				vTmp2[0]=BL_CACA_CISP+Rand1Distrib()*BLSD_CACA_CISP;
			else
				vTmp2[0]=BL_CACA+Rand1Distrib()*BLSD_CACA;
		} while (AssignCoords(" CA ",pmgdHere,vTmp2,Model)!=ERR_SUCCESS);
		/* ensure ca1 and ca2 do not have van der waals collision */
		ptgsHere=FreeTraj(ptgsHere);
	}
	else { /* WALKTYPE==WALK_PHIPSI */
		/* place N and C in canonical position about origin (Ca) */
		vTmp2[0]=-BL_NCA;
		vTmp2[1]=0;
		vTmp2[2]=0;
		if (AssignCoords(" N  ",pmgdLast,vTmp2,Model)!=ERR_SUCCESS) {
			PurgeGlobs();
			ErrPostEx(SEV_FATAL,1,1,"Unable to place second atom!");
		}
		vTmp3[0]=-BL_CAC*cos(DEGTORAD*BA_NCAC);
		vTmp3[1]=BL_CAC*sin(DEGTORAD*BA_NCAC);
		vTmp3[2]=0;
		if (AssignCoords(" C  ",pmgdLast,vTmp3,Model)!=ERR_SUCCESS) {
			PurgeGlobs();
			ErrPostEx(SEV_FATAL,1,1,"Unable to place third atom!");
		}
	}
	lrx=0;
	lrxtries=0;
	lastprt=0;
	if (traj_quiet==VERBOSITY_QUIET || traj_quiet==VERBOSITY_VERBOSE) {
		printf("\nFolding..\nCalculating residue    2/%-4d",clength);
		fflush(stdout);
	}

/* Random walk */
	tm1=GetSecs();
	do {}
	while ((err=NextResidue(pmgdHere,vTmp3,vTmp,vTmp2,999.0,999.0,Model))==ERR_FAIL);
        /* free all the worlds we created */
	tm2=GetSecs();
/* Random walk ended */
	if (traj_quiet == VERBOSITY_STREAM) {
		printf("%07ld\n",StructureNumber);
	}

/*PrintDistConstTries();*/
	if (traj_quiet==VERBOSITY_QUIET)
		printf("\n");
	else if (traj_quiet==VERBOSITY_VERBOSE)
		printf("\nCleaning up..\n");
	/* remove garbage atoms */
	if (BUILD_FRAGS_ONLY) {
		TraverseOneModel(pmsdRoot->pdnmmHead,TRAVERSE_ATOM,Model,0,NULL,(pNodeFunc)(FreeUnwantedLocs));
	}
/*for (cnt=0;cnt<clength;cnt++) {
printf("%d: %d %d\n",cnt+1,ressrcres[cnt],ressrcfrag[cnt]);
}*/
	ressrcres=MemFree(ressrcres);
	ressrcfrag=MemFree(ressrcfrag);

    FreeAllWorlds();
	/* clean up temporary "dirty" bits in bReserved of PMADs */
	ClearUpdate(pmsdRoot->pdnmmHead,Model);
	if (err!=ERR_INCOMPLETE) {
		/* log it */
		BSprintf(bspTempLog,"%07ld\t",StructureNumber);		
		/* calculate NC dist and radius of gyration */
        pdnmgHere=((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue))->pdnmgHead;
 		while (pdnmgHere->next)
			pdnmgHere=pdnmgHere->next;
		do {
			pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
			if (GetCoOrds(pmgdHere," CA ",vZero,vCACTerm,Model))
				break;
			pdnmgHere=pdnmgHere->last;
		} while (pdnmgHere);
		/* effective length of sequence - with co-ordinates */
		efflength=pdnmgHere->choice;
		pdnmgHere=((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue))->pdnmgHead;
		while (pdnmgHere) {
			pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
			if (GetCoOrds(pmgdHere," CA ",vZero,vCANTerm,Model))
				break;
			pdnmgHere=pdnmgHere->next;
		}
		efflength=efflength-(pdnmgHere->choice-1);

		VecSub(vCATerm,vCACTerm,vCANTerm);
		Rgyr=GetRgyr((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue),Model);
		RgyrHP=GetHPRgyr((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue),Model);

		Rn=Rgyr*Rgyr/((FloatLo)efflength*BL_CACA*BL_CACA);
		Cn=getMag(vCATerm)*getMag(vCATerm)/((FloatLo)efflength*BL_CACA*BL_CACA);
#ifdef USE_DSSP
		if (traj_quiet==VERBOSITY_VERBOSE)
			printf("Calling DSSP to compute secondary structure and accessibility..\n");	
		surfacc=CalcDSSPAccSurf((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue),2,Model);
		Surface=0.0;
		for (cnt=0;cnt<clength;cnt++)
			Surface+=(FloatHi)surfacc[cnt];
		pdnmgHere=((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue))->pdnmgHead;
		HPSurface=0.0;
		for (cnt=0;cnt<clength;cnt++) {
			cResHere=(((PMGD)(pdnmgHere->data.ptrvalue))->pcIUPAC)[0];
  			if (StringRChr(hpresidues,(int)cResHere)!=NULL)
				HPSurface+=(FloatHi)surfacc[cnt];
			pdnmgHere=pdnmgHere->next;
		}
		MemFree(surfacc);
		CalcDSSPAssignEx((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue),Model,pcDSSPAssign,TRUE);
		NumDSSPHelix=0;
		NumDSSPSheet=0;
		pdnmgHere=((PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue))->pdnmgHead;
		/* shove DSSP identification in bPDBSecStru of PMGDs in case needed for later reference */
		for (cnt=0;cnt<clength;cnt++) {
			if (pcDSSPAssign[cnt]=='H') {
				NumDSSPHelix++;
				((PMGD)(pdnmgHere->data.ptrvalue))->bPDBSecStru=(Byte)SS_HELIX;
			}
			else if (pcDSSPAssign[cnt]=='B' || pcDSSPAssign[cnt]=='E') {
				NumDSSPSheet++;
				((PMGD)(pdnmgHere->data.ptrvalue))->bPDBSecStru=(Byte)SS_STRAND;
			}
			else
				((PMGD)(pdnmgHere->data.ptrvalue))->bPDBSecStru=(Byte)0;
			pdnmgHere=pdnmgHere->next;
		}
	
#endif
		/* for DF client, store trajectory version in length field - its not used beyond this point */
		/* identify client because TGUNITS will be three-digits or more */
		if (TGUNITS>99) {
			clength=TGUNITS;
		}
		if (traj_quiet==VERBOSITY_VERBOSE) {
			printf("\n(Structure): %07ld\n",StructureNumber);		
	        printf("(Time): %ld seconds\n(Tries): %ld\n",(long int)(tm2-tm1),(long int)tries);
			printf("No backbone solution tries (BadBB): %ld\n(Crashes): %ld\n",(long int)fpabad,(long int)crashcnt);
			printf("Tries due to Distance Constraints (ViolatedConstr): %ld (%ld)\n",(long int)dcbad,(long int)dcviol);
			printf("Sequence length (N): %d\n(Rgyr): %f\n",clength,Rgyr);
			printf("NCTerm distance (NCdist): %f\n(Rn): %f\n(Cn): %f\n",getMag(vCATerm),Rn,Cn);
#ifdef USE_DSSP
			printf("DSSP Surface Accessibility (ASA): %8.1f\nDSSP Exposed Hydrophobics (HASA): %8.1f\n",Surface,HPSurface);
			printf("DSSP Helical Residues (Helix): %ld\nDSSP Extended Residues (Edssp): %ld\n",(long int)NumDSSPHelix,(long int)NumDSSPSheet);
#endif
		}
#ifdef USE_DSSP
		if (traj_quiet==VERBOSITY_QUIET) {
			printf("\n(Structure): %07ld\n",StructureNumber);		
	        printf("(Time): %ld seconds\n(Tries): %ld\n",(long int)(tm2-tm1),(long int)tries);
		}		BSprintf(bspTempLog,"%ld\t%ld\t%ld\t%ld\t%ld(%ld)\t%d\t%f\t%f\t%f\t%f\t%f\t%8.1lf\t%8.1lf\t%ld\t%ld",(long int)(tm2-tm1),(long int)tries,(long int)fpabad,(long int)crashcnt,(long int)dcbad,(long int)dcviol,clength,Rgyr,RgyrHP,getMag(vCATerm),Rn,Cn,Surface,HPSurface,(long int)NumDSSPHelix,(long int)NumDSSPSheet);
#else
		BSprintf(bspTempLog,"%ld\t%ld\t%ld\t%ld\t%ld(%ld)\t%d\t%f\t%f\t%f\t%f\t%f",(long int)(tm2-tm1),(long int)tries,(long int)fpabad,(long int)crashcnt,(long int)dcbad,(long int)dcviol,clength,Rgyr,RgyrHP,getMag(vCATerm),Rn,Cn);
#endif
		if (traj_quiet==VERBOSITY_QUIET || traj_quiet==VERBOSITY_VERBOSE) {
			for (cnt=1;cnt<=clength;cnt++)
				if (IsCis[cnt]) {
					if (WALKTYPE==WALK_CA)
						printf("Residue %d is cis-Pro!\n",cnt-1);
					else
						printf("Residue %d is cis-Pro!\n",cnt);
				}
		}
		if (fExtraLog!=NULL) {
			prsHead=GetTrjAngle(pmmdParent,0,Model,0.0,NULL);
			vnpChis=MeasureChis(pmmdParent,Model);
			prsHere=prsHead;
			vnpHere=vnpChis;
			/* may want to print phi/psi/theta 0 and N here explicitly too */
			if (WALKTYPE==WALK_CA) {
				firstPsi=0.0;
				lastPhi=0.0;
				pdnmgHere=pmmdParent->pdnmgHead;
				pmgdLast=(PMGD)(pdnmgHere->data.ptrvalue);
				pmgdHere=(PMGD)((pdnmgHere->next)->data.ptrvalue);
				GetCoOrds(pmgdLast," CA ",vZero,v1,Model);
				GetCoOrds(pmgdLast," C  ",vZero,v2,Model);
				GetCoOrds(pmgdHere," N  ",vZero,v3,Model);
				GetCoOrds(pmgdHere," CA ",vZero,v4,Model);
				GetDihedral(v1,v2,v3,v4,-180,&firstOmega,&ba,&ba,&bl,&bl,&bl);
				pdnmgHere=pdnmgHere->next;
			}
			else if (WALKTYPE==WALK_PHIPSI) {
				/* calculate first and last phi, psi */
				pdnmgHere=pmmdParent->pdnmgHead;
				pmgdLast=(PMGD)(pdnmgHere->data.ptrvalue);
				pmgdHere=(PMGD)((pdnmgHere->next)->data.ptrvalue);
				GetCoOrds(pmgdLast," N  ",vZero,v1,Model);
				GetCoOrds(pmgdLast," CA ",vZero,v2,Model);
				GetCoOrds(pmgdLast," C  ",vZero,v3,Model);
				GetCoOrds(pmgdHere," N  ",vZero,v4,Model);
				GetDihedral(v1,v2,v3,v4,-180,&firstPsi,&ba,&ba,&bl,&bl,&bl);
				VecAdd(v1,v2,vZero);  /* v1=v2 */
				VecAdd(v2,v3,vZero);  /* v2=v3 */
				VecAdd(v3,v4,vZero);  /* v3=v4 */
				GetCoOrds(pmgdHere," CA ",vZero,v4,Model);
				GetDihedral(v1,v2,v3,v4,-180,&firstOmega,&ba,&ba,&bl,&bl,&bl);
 				/* now go to last residue */
				while (pdnmgHere->next)
					pdnmgHere=pdnmgHere->next;
				/* now point pmgdHere to last residue */
				pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
				pmgdLast=(PMGD)((pdnmgHere->last)->data.ptrvalue);
				GetCoOrds(pmgdLast," C  ",vZero,v1,Model);
				GetCoOrds(pmgdHere," N  ",vZero,v2,Model);
				GetCoOrds(pmgdHere," CA ",vZero,v3,Model);
				GetCoOrds(pmgdHere," C  ",vZero,v4,Model);
				GetDihedral(v1,v2,v3,v4,-180,&lastPhi,&ba,&ba,&bl,&bl,&bl);
			}
			fprintf(fExtraLog,"--------------------------------------\n");
			if (firstiter==TRUE)
				fprintf(fExtraLog,"%d %d\n",clength,WALKTYPE);
			if (WALKTYPE==WALK_CA)
				fprintf(fExtraLog,"%6.3f %6.1f %6.1f\n",0.0,firstPsi,firstOmega);
			else
				fprintf(fExtraLog,"%6.1f %6.1f %6.1f\n",0.0,firstPsi,firstOmega);
			while (prsHere)	{
				if (WALKTYPE==WALK_CA) {
					if (pdnmgHere->next!=NULL) {
						pmgdLast=(PMGD)(pdnmgHere->data.ptrvalue);
						pmgdHere=(PMGD)((pdnmgHere->next)->data.ptrvalue);
						GetCoOrds(pmgdLast," CA ",vZero,v1,Model);
						GetCoOrds(pmgdLast," C  ",vZero,v2,Model);
						GetCoOrds(pmgdHere," N  ",vZero,v3,Model);
						GetCoOrds(pmgdHere," CA ",vZero,v4,Model);
						GetDihedral(v1,v2,v3,v4,-180,&OmegaHere,&ba,&ba,&bl,&bl,&bl);
						pdnmgHere=pdnmgHere->next;
					}
					else
						OmegaHere=CHI_W;
					/* extra precision on phi since it is sometimes very small */
					fprintf(fExtraLog,"%6.3f %6.1f %6.1f\n",prsHere->Phi,prsHere->Psi,OmegaHere);
				}
				else
					fprintf(fExtraLog,"%6.1f %6.1f %6.1f\n",prsHere->Phi,prsHere->Psi,prsHere->Omega);
				prsHere=prsHere->next;
			}
			if (WALKTYPE==WALK_CA)
				fprintf(fExtraLog,"%6.3f %6.1f %6.1f\n",lastPhi,0.0,0.0);
			else
				fprintf(fExtraLog,"%6.1f %6.1f %6.1f\n",lastPhi,0.0,0.0);
			/* keep track of residue number */
			cnt=1;
			pdnmgHere=pmmdParent->pdnmgHead;
			while (vnpHere)	{
				prlr=(pnewrotlibrecord)(vnpHere->data.ptrvalue);
				if (firstiter==TRUE) {
					ptgsHere=TrajGraphRead(cnt);
					fprintf(fExtraLog,"%6.1f %6.1f %6.1f %6.1f %c %5.3f %d %d\n",prlr->chi1,prlr->chi2,prlr->chi3,prlr->chi4,ptgsHere->AA,ptgsHere->markovsf,ptgsHere->dim,ptgsHere->tout);
					ptgsHere=FreeTraj(ptgsHere);
					pdnmgHere=pdnmgHere->next;
				}
				else
					fprintf(fExtraLog,"%6.1f %6.1f %6.1f %6.1f\n",prlr->chi1,prlr->chi2,prlr->chi3,prlr->chi4);
				vnpHere=vnpHere->next;
				cnt++;
			}		
			freeRS(prsHead);
			FreeChis(vnpChis);
		} 
		/* code altered not to report any bad parameters to the log file from incomplete structures CWVH 2012 */
	}
	if (phbsHBondsToCheck!=NULL)
			phbsHBondsToCheck=PHBSFree(phbsHBondsToCheck);
	return err;
}

VoidPtr TRADEProtein(VoidPtr ptr, Int4 number)
{
	pFoldTrajParamBlock paramblock;
	Int2 gen;

	paramblock=(pFoldTrajParamBlock)ptr;
	if (paramblock==NULL)
		return NULL;
	gen=paramblock->gen;
	if (gen<1)
		gen=1;
	gencbrt=(FloatHi)pow((double)gen,(double)1.0/3.0);
	StringCpy(errorfile,paramblock->errorfile);
  	paramblock->err=TRADEProteinEx(paramblock->pmsdRoot,paramblock->Model, number);
	ProgramProgress=0;
	if (ProgramProgressMax>0)
		ProgramProgressMax=0;
	return NULL;
}

/*  
$Log: randwalk.c,v $
Revision 1.130  2004/10/04 17:18:02  egarderm
Fixed wrong comment

Revision 1.129  2004/09/08 20:58:30  mbrougham
A method was added to calculate the probability of tunneling.  This is scaled by TUNNEL_PROB via Vistraj, which can result in no tunneling if that parameter is set to 0.

Revision 1.128  2004/08/30 16:44:01  hfeldman
Added support for backtracking fragment in phi-psi walk

Revision 1.127  2004/07/15 21:50:27  ksnyder
Write errors, causing HomBail() to execute, to Maketrj Error log

Revision 1.126  2004/07/13 20:53:55  hfeldman
No backtracking memory of fragments for phi-psi walk

Revision 1.125  2004/06/28 14:08:59  hfeldman
Add all atoms to BD-tree, even if not bumpchecking them (for non-fragment atoms)

Revision 1.124  2004/06/23 17:40:05  hfeldman
Added option to only keep residues which come from fragments - co-ords are freed from the rest of residues

Revision 1.123  2003/12/11 17:32:21  hfeldman
Reduced tries before timeout by a factor of 5, seems to help fold things faster

Revision 1.122  2003/11/05 18:35:07  feldman
Include order fix for Darwin

Revision 1.121  2003/10/23 17:30:35  feldman
Score by best crease energy column now

Revision 1.120  2003/09/23 19:13:03  feldman
Added errorfile option to foldtraj for homtraj
Made foldtraj gradually increase laxness if gets stuck a lot, like in DFP

Revision 1.119  2003/09/22 16:05:10  feldman
Fixed compiler warnings

Revision 1.118  2003/09/16 21:23:10  feldman
Added code to track sources of residues in a structure - from a fragment or from trajectory distribution
This allows fragments to be re-loaded when backtracking so they are rebuilt after it goes forwards again
This has not been thoroughly tested for phi-psi walk but appears to work on a simple test case;  for ca
walk it seems to be working properly, generating better quality homology models

Revision 1.117  2003/08/23 17:03:27  feldman
Added constant to timeout - will be better for shorter peptides

Revision 1.116  2003/07/18 21:29:51  feldman
Decreased timeout time by factor of 2.5

Revision 1.115  2003/06/20 21:55:41  feldman
Print 'tight spot' message much less frequently

Revision 1.114  2003/04/07 22:09:12  feldman
Turned back off 3D objects - dont need them here for DFP

Revision 1.113  2003/04/06 21:01:03  feldman
Turn SS 3D objects back on for DFP beta

Revision 1.112  2003/04/04 21:54:04  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.111  2003/03/25 17:46:56  feldman
changed 'give up' timeout strategy

Revision 1.110  2003/03/13 21:34:27  feldman
Fixed bug with 'stuck' message

Revision 1.109  2003/03/11 22:19:03  feldman
Remove 3D objects from DFP movies/structures

Revision 1.108  2003/03/11 19:31:11  feldman
Added 'tight spot' message for DFP

Revision 1.107  2003/02/27 22:23:15  feldman
Restart now in 1500nlogn tries or nlogn seconds, whichever occurs first
Also for DFP, no more crazy ASCII art when it gets stuck

Revision 1.106  2003/02/19 23:27:53  feldman
Reduced endless loop message to a warning

Revision 1.105  2003/02/19 18:38:00  feldman
Made modification for exiting at ANY time from an external signal

Revision 1.1  2003/02/10 23:37:21  feldman
New folding algorithm

Revision 1.104  2003/01/09 20:00:43  feldman
Added extra space in progress message

Revision 1.103  2003/01/08 16:14:53  feldman
disable h-bond testing when bounciness is high

Revision 1.102  2002/12/19 20:14:50  feldman
Fixed include paths

Revision 1.101  2002/12/19 20:04:27  feldman
New screensaver

Revision 1.98  2002/09/05 17:32:08  feldman
Fixed cis-Pro in homology modelling

Revision 1.97  2002/08/11 16:40:46  feldman
Fixed logic ambiguity

Revision 1.96  2002/08/11 15:42:14  feldman
Bounciness > 1000 means disable bumpchecking completely

Revision 1.95  2002/08/01 21:55:39  feldman
Fixed bug - now will allow tunnelling from last residue of fragment

Revision 1.94  2002/07/29 13:59:57  feldman
removed TUNNEL_PROB

Revision 1.93  2002/07/20 16:42:42  feldman
Added LOTS_RAM option for trajgraphreading

Revision 1.92  2002/07/17 01:10:15  feldman
Fixed DFP bug for protein version in log in quiet mode is on

Revision 1.91  2002/07/15 21:56:28  feldman
Fixed fragment cis-proline bugs

Revision 1.90  2002/07/15 20:10:09  elena
Cleaned up debug info.

Revision 1.89  2002/07/15 19:59:29  elena
Added tunneling info the the log.

Revision 1.88  2002/07/11 13:47:43  elena
Modified tunneling probability to depend on remaining fragment length

Revision 1.87  2002/07/10 21:09:47  feldman
store protein version in protein length field of log for DFP

Revision 1.86  2002/06/13 20:06:01  elena
Added code to allow tunneling between fragments.

Revision 1.85  2002/04/02 15:27:43  feldman
Fix for hppa linux

Revision 1.84  2002/02/25 22:08:01  feldman
Added checksumming of text files potentials, cbdata and skel.prt
If checksum fails, get error
Changed bailing from foldtrajlite to exit curses first if in it

Revision 1.83  2002/02/22 19:51:21  feldman
Turn off annoying warnings for incomplete structures

Revision 1.82  2002/02/21 17:44:19  feldman
Made error msg a bit more verbose

Revision 1.81  2002/02/14 18:59:06  feldman
Removed extraneous variable

Revision 1.80  2002/02/12 16:28:08  feldman
Fixed broken foldtraj on windows

Revision 1.79  2002/02/10 03:46:27  feldman
Added NT Service capability

Revision 1.78  2002/02/05 17:53:09  feldman
Got rid of refresh call for quiet mode

Revision 1.77  2002/02/01 19:43:48  feldman
Added quiet option to client

Revision 1.76  2002/01/18 21:10:43  feldman
Changed Discriminant < 0 to non-fatal error, since it happens on some rare occasions on non-Intel hardware

Revision 1.75  2001/12/10 22:19:50  feldman
removed hack - should work on Darwin now

Revision 1.74  2001/12/10 21:08:26  feldman
Added brackets for compiler warning

Revision 1.73  2001/12/07 23:05:03  feldman
Fixed very subtle and rare H-bond bug

Revision 1.72  2001/12/06 22:35:05  feldman
Fixed compiler warnings
should compile on Darwin now
made default resolution 1024x768x16

Revision 1.71  2001/11/19 22:58:08  feldman
Turned off surface area computation for screensaver

Revision 1.70  2001/10/30 21:34:41  feldman
Fixed text-mode message

Revision 1.69  2001/10/23 15:24:02  feldman
Added backup and restore co-ords functions

Revision 1.68  2001/10/10 19:48:19  feldman
Revamped hydrogen bonding code - now much more robust and most likely correct

Revision 1.67  2001/09/13 16:22:42  feldman
Minor fragment related bugfixes

Revision 1.66  2001/09/10 20:32:00  feldman
Fixed bug - check if a fragment contains a cis- residue and if so, check
that it maps onto a proline - if not, pick the next fragment or no fragment
at all

Revision 1.65  2001/09/04 18:04:22  feldman
Fixed H-bond related bug

Revision 1.64  2001/08/27 22:07:09  feldman
Fixed flashing near end of protein

Revision 1.63  2001/08/21 16:03:58  feldman
fixed bug: When structure is completed, resets state to "not in a fragment" now

Revision 1.62  2001/08/15 20:12:48  feldman
Removed unused variable

Revision 1.61  2001/08/13 16:48:24  feldman
Add fragments in a different way and probability is the more intuitive
meaning now

Revision 1.60  2001/08/01 14:49:09  feldman
Fixed bug - Ha must be placed after rotamer for phi-psi

Revision 1.59  2001/07/31 19:44:43  feldman
Removed NextResidue parameters and fixed fragment bugs for setting
Cis residues, setting phi/psi and prevent from jumping ahead on fragments bug

Revision 1.58  2001/07/26 19:19:51  feldman
Changed way fragments are stored on disk to be more platform independent,
and removed id parameter from AddtoWorld

Revision 1.57  2001/07/24 17:27:12  feldman
Updated rotamer library version and fixed bug - order of building
for Ramachandran space had prevented rotamer library from being
used properly

Revision 1.56  2001/07/13 19:01:14  feldman
Id field removed from bd functions

Revision 1.55  2001/07/12 21:06:48  feldman
started changing to allow deletion of root node

Revision 1.54  2001/06/19 18:42:54  phan
changed rgyr functions to take pmmd

Revision 1.53  2001/06/14 19:13:47  feldman
Changed screen refresh function to avoid potential endless loop in
screensaver

Revision 1.52  2001/05/31 15:03:24  feldman
Changed a few variables around, made them volatile even though extern
(not sure if need to do this)

Revision 1.51  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.50  2001/05/04 16:18:59  feldman
Added Gly, Ala and Trp as hydrophobic residues when getting rgyr and
surface area of hydrophobics, based on Eisenberg ranking

Revision 1.49  2001/04/17 21:32:06  feldman
Added fragment editing to Vistraj - not quite fully functional yet
Also can now give trajectory graph name as argument #1 when running
and fixed potential sscanf problems so entering garbage into numeric
entry fields results in an error being shown (instead of accepting it)

Revision 1.48  2001/04/12 22:27:47  feldman
added IF CURSES where needed

Revision 1.47  2001/04/11 14:54:16  feldman
removed reference to "screensaver"

Revision 1.46  2001/04/09 15:43:14  feldman
Leave OSF_SOURCE defined for alpha processors

Revision 1.45  2001/04/08 22:17:52  feldman
fixed curses-related compiler warnings, and a small bug

Revision 1.44  2001/04/06 17:14:50  feldman
made colors brighter on WINDOWS and fixed a bug

Revision 1.43  2001/04/06 16:31:31  feldman
Made colors brighter

Revision 1.42  2001/04/06 14:29:00  feldman
Completed fold-at-home Client-Server v1.0 now complete and ready for testing

Revision 1.41  2001/04/05 14:02:39  feldman
Fixed so curses should work on both HP and Alpha now, I hope.
Also added a few more error codes to trajstore

Revision 1.40  2001/04/04 22:46:06  feldman
Added messages when data sent to server, success and failure

Revision 1.39  2001/04/04 21:26:01  feldman
-Removed my BSPrintf, now using the one in SLRIlib - thus needed to add
	this library to makefiles
-changed all fgets to FileGets

Revision 1.38  2001/04/03 16:27:13  feldman
Should work for ALL UNIXes this time...

Revision 1.37  2001/04/03 16:05:17  feldman
Fixed Windows compatibility

Revision 1.36  2001/04/03 14:31:20  feldman
changed curses to be first and added HP flag

Revision 1.35  2001/04/02 15:11:13  feldman
fixed minor color and keypress bug for UNIX

Revision 1.34  2001/04/02 00:17:05  feldman
added HP-compatibility

Revision 1.33  2001/04/01 19:44:42  feldman
made foldtrajlite Windows port using conio.h, began network portion

Revision 1.32  2001/03/30 22:22:27  feldman
Minor spelling and grammar fixes, bug fixes, and added
foldtraj ASCII screensaver ability

Revision 1.31  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.30  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.29  2001/03/23 21:08:07  feldman
Added monitor for foldtraj inside vistraj

Revision 1.28  2001/03/15 15:59:57  feldman
revamped makefiles so they are a lot smaller now, all common stuff
is in .mk files - all have been tested and most still compile (Except
one or two with known problems)

Revision 1.27  2001/03/13 17:48:54  feldman
fixed rotid bug

Revision 1.26  2001/03/13 15:08:17  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.25  2001/03/09 17:33:58  feldman
Added code to deal with rotamers.  Now rotamer chi angles can
be encoded into a 4-byte integer and used by foldtraj when
building the structures.  These can be edited by VisTraj as well

Revision 1.24  2001/03/07 21:49:47  feldman
Made many important fixes to VisTraj such as bounds checking of
user-entered information, and tidied up interface further.
Fixed a few minor bugs as well, and tested RPSTraj feature of
VisTraj.

Revision 1.23  2001/02/26 22:21:15  feldman
-Changed fragments to allow multiple possible fragments at the same
residue
-altered random walk to make use of fragments when present (a first
attempt at least...)

Revision 1.22  2001/02/08 22:28:11  feldman
Fixed important inconsistency in modified amino acids - now for
encoded sequence (*xxxxA) the xxxx is always the dictionary index
of the residue, thus we no longer need to add 1 or 2 to it in
cases where it is found at the N- or C-terminus; adding and
deleting residues with AlterResiduesInSequence will update any
affected residues correctly as well

Revision 1.21  2001/02/08 19:46:36  feldman
Added code to improve SubstituteNames and make it do a lot
for checking of atomnames when entering constraints

Revision 1.20  2001/02/07 18:46:45  feldman
Removed a few unused variables

Revision 1.19  2001/01/18 18:10:38  feldman
Added bzip2 1.01 compatibility

Revision 1.18  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.17  2001/01/12 20:01:51  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.16  2000/11/17 22:28:59  feldman
Minor bugfixes

Revision 1.15  2000/10/24 20:57:25  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.14  2000/09/26 18:32:54  feldman
Changed backtrack probability again

Revision 1.13  2000/09/15 20:32:08  feldman
Added printout of distance constraint violation count
Made backtracking random now; could backtrack anytime
so it will on average backtrack a bit more
(will test effect of this over the next few weeks)

Revision 1.12  2000/09/06 15:32:08  feldman
changed warnings to SEV_INFO to avoid excessive log files

Revision 1.11  2000/08/31 14:55:11  feldman
Minor corrections/changes

Revision 1.10  2000/08/23 21:40:59  feldman
Print new DSSP info using internal calls rather than
external ones

Revision 1.9  2000/08/09 19:12:19  feldman
-minor bugfix update and fixed up makefiles removing USEDSSP
-added bioseq/seq-entry to trajectory graph header when using
 unfoldtraj or val2trj

Revision 1.8  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.7  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.6  2000/07/06 15:26:14  feldman
Fixed bug - in no native file, don't try to save best RMSD structures!

Revision 1.5  2000/06/23 19:47:04  feldman
Went back to using bzip2 for compression and no
longer use a bytestore but dump directly to a file

- added surface accessibility output to fold logs

Revision 1.4  2000/06/22 13:05:01  feldman
Fixed minor bugs with log file nameing,
fixed error in extended residue calculation and DSSP usage
added alterresidue function

Revision 1.3  2000/06/20 16:40:23  feldman
Incorporated sstru extended structure calculation into C code,
no longer need an external program

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

