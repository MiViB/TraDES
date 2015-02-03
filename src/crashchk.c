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
#include "crashchk.h"



#define MAXDATA 10000
#define RNG -160

/* Global Variables */
Args Rargs[NUMARGS] = {
  /*0*/               {"Input MMDB filename",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
  /*1*/               {"Model Number","1","1","9999",TRUE,'m',ARG_INT,0.0,0,NULL},
  /*2*/	     	      {"Model Level (remote files only): 0 = Atoms; 1=BB; 2=all PDB; 3=VECT;","0","0","3",TRUE,'l',ARG_INT,0.0,0,NULL},
  /*3*/               {"Supress same-chain collision reporting (default=FALSE)",NULL,NULL,NULL,TRUE,'n',ARG_BOOLEAN,0.0,0,NULL}
                 };

static Char smallatoms[]="EZHQILMXRTUFYWJKV";

#define CHECKHBONDS					\
	if (pvnal!=NULL) {				\
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

static Boolean FindHBond(PMAD pmad1,PMAD pmad2,Boolean PlacedOnly)
{
	vec vHbond,vAtom1,vAtom2,vTmp1,vTmp2;
	PALD	pald1,pald2,paldHbond;
	PVNAL pvnal;
	FloatLo HbondAngleCos;
	ValNodePtr vnpHere;
	PMBD pmbdHere;
	Char res1,res2;
	Char at1[5],at2[5];
	
	if (pmad1==pmad2)
		return FALSE;
	if ((pmad1->pvnalLocate)==NULL)
		return FALSE;
	if ((pmad2->pvnalLocate)==NULL)
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
		if (pmbdHere->pmadFrom->pcAName[1]=='H' && ((res1!='D' && res1!='E') || at1[1]!='O' || (at1[2]!='D' && at1[2]!='E')) && (res1!='H' || StringCmp(at1," ND1")) && StringCmp(at1," OXT")) {
			/* check if h-bond possible */
			pvnal=pmbdHere->pmadFrom->pvnalLocate;
			/* save some typing here */
			CHECKHBONDS
		}
		if (pmbdHere->pmadTo->pcAName[1]=='H' && ((res1!='D' && res1!='E') || at1[1]!='O' || (at1[2]!='D' && at1[2]!='E')) && (res1!='H' || StringCmp(at1," ND1")) && StringCmp(at1," OXT")) {
			/* check if h-bond possible */
			pvnal=pmbdHere->pmadTo->pvnalLocate;
			CHECKHBONDS
		}
		vnpHere=vnpHere->next;
	}
	vnpHere=pmad2->pvnBonds;
	while (vnpHere) {
		pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
		if (pmbdHere->pmadFrom->pcAName[1]=='H' && ((res2!='D' && res2!='E') || at2[1]!='O' || (at2[2]!='D' && at2[2]!='E')) && (res2!='H' || StringCmp(at2," ND1")) && StringCmp(at2," OXT")) {
			/* check if h-bond possible */
			pvnal=pmbdHere->pmadFrom->pvnalLocate;
			CHECKHBONDS
		}
		if (pmbdHere->pmadTo->pcAName[1]=='H' && ((res2!='D' && res2!='E') || at2[1]!='O' || (at2[2]!='D' && at2[2]!='E')) && (res2!='H' || StringCmp(at2," ND1")) && StringCmp(at2," OXT")) {
			/* check if h-bond possible */
			pvnal=pmbdHere->pmadTo->pvnalLocate;
			CHECKHBONDS
		}
		vnpHere=vnpHere->next;
	}
	return FALSE;
}

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

Int2 Main()
{
	PMSD pmsdRoot;
	Int2 ModelNum,res1,res2,resdist,doit;
  PDNMG pdnmgHere,pdnmgHead;
  PDNMM	pdnmmHere;
  PMGD pmgdHere,pmgdCrash;
	PVNMA pvnmaHere;
  vec vHere,vTmp;
	PMAD pmadHere;
	PALD paldHere;
	FloatLo radAtom1,radAtom2,tooclose;
	CharPtr AtomName,AtomHere;
	ValNodePtr vnpHere,vnpHead,vnpAtomList;
	PWS pwsThis=NULL;
  Int2 at1hn,at2hn;
  Char reshere,rescrash;
  CharPtr chain1,chain2;
  Boolean bNoIntrachain; 

      ErrSetLogfile("error_crashchk.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
	if (!GetArgs("CrashChk 1.0 - finds atomic collisions in ASN.1 structure files with b-d tree",NUMARGS,Rargs))
		return 1;
	if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI");
		return 2;
	}

       
	/* load an ASN.1 Biostruc */
	ModelNum=Rargs[1].intvalue;
	pmsdRoot=LoadABiostruc(Rargs[0].strvalue,FALSE,Rargs[2].intvalue,&ModelNum);
	if (pmsdRoot==NULL) {
		/* error occurred */
		return 3;
	}
	/* note pmsdRoot contains all H2O, P, other heterogens when
	   determining the bounding box */
	printf("Struc loaded\n");

/* Supress intrachain bumpchecking */
       bNoIntrachain = 0; /* all chain bumpchecking by default */
       if (Rargs[3].intvalue) bNoIntrachain = 1; 
	


  pwsThis=AddtoWorld(pwsThis,ModelNum,(PFB)pmsdRoot);
  /* instantiate the world */
  vnpHead=InstantiateWorld(1,pwsThis);
  pdnmmHere=pmsdRoot->pdnmmHead;
  while(pdnmmHere){
	pdnmgHead=((PMMD)(pdnmmHere->data.ptrvalue))->pdnmgHead;
	pdnmgHere=pdnmgHead;
	while (pdnmgHere) {
		pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
		reshere=(pmgdHere->pcIUPAC)[0];
		pvnmaHere=pmgdHere->pvnmaAHead;
		while (pvnmaHere) {
			pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
			AtomName=pmadHere->pcAName;
			at1hn=CountHNeighbour(pmadHere);
	    /* get atom VDW radii for AtomName */
			radAtom1=GetRadNewAtom(reshere,AtomName,at1hn);	
    	paldHere=GetAtomLocs(pmadHere,ModelNum);
			if (paldHere!=NULL) {
				vTmp[0]=(paldHere->pflvData)[0];
			  vTmp[1]=(paldHere->pflvData)[1];
      	vTmp[2]=(paldHere->pflvData)[2];
				vnpAtomList=FindAtomsIn(pwsThis,vTmp,radAtom1+VDWRAD_MAX);
			  vnpHere=vnpAtomList;

 /* walk through list to output results */
			  while(vnpHere) {
					radAtom1=GetRadNewAtom(reshere,AtomName,at1hn);	
 			    paldHere=(PALD)(vnpHere->data.ptrvalue);
         	            pmgdCrash=GetParentGraph((PFB)paldHere);
	  		    rescrash=(pmgdCrash->pcIUPAC)[0];
			    /* check if different residues between colliding atoms */
		res1=(pmgdCrash->pdnmgLink)->choice;
		chain1=((PMMD)pmgdCrash->pfbParent)->pcMolName;
		res2=(pmgdHere->pdnmgLink)->choice;
		chain2=((PMMD)pmgdHere->pfbParent)->pcMolName;
		    	AtomHere=((PMAD)(paldHere->pfbParent))->pcAName;
					at2hn=CountHNeighbour((PMAD)(paldHere->pfbParent));
	        /* get atom VDW radii for paldHere */
          switch (AtomHere[1]) {
            case 'C':
							if (at2hn==30)
								radAtom2=VDWRAD_AROCH0;
							else if (at2hn==31)
								radAtom2=VDWRAD_AROCH1;
							else
	              radAtom2=VDWRAD_C;
              break;
            case 'N':
              radAtom2=VDWRAD_N;
              break;
            case 'O':
              radAtom2=VDWRAD_O;
              break;
            case 'S':
              radAtom2=VDWRAD_S;
              break;
            case 'H':
              radAtom2=VDWRAD_H;
              break;
         		default:
             	radAtom2=VDWRAD_C;
          }
					if ((StringRChr(smallatoms,(int)(AtomHere[2]))!=NULL) && !(AtomHere[2]=='X' && AtomHere[3]=='T') && ((rescrash=='E') || (rescrash=='Q') || (rescrash=='K') || (rescrash=='R'))) {
						if (radAtom2>VDWRAD_X)
							radAtom2=VDWRAD_X;
					}
					if ((AtomHere[2]=='D') && (rescrash=='P')) {
						if (radAtom2>VDWRAD_X)
							radAtom2=VDWRAD_X;
					}
		vHere[0]=(paldHere->pflvData)[0];
		vHere[1]=(paldHere->pflvData)[1];
		vHere[2]=(paldHere->pflvData)[2];
		VecSub(vHere,vHere,vTmp);
		tooclose=(getMag(vHere))-(radAtom1+radAtom2);
		if ((bNoIntrachain==1) && (chain1==chain2)) ; /* same chain so supress output */
		else if ((res1<res2) && (tooclose<0) ) {
			
				resdist=(((pmgdCrash->pdnmgLink)->choice)-((pmgdHere->pdnmgLink)->choice));
				doit=0;
				if (resdist==-1) {
					if (((!StringCmp(AtomHere," CB ")) || (!StringCmp(AtomHere,"2HA "))) && ((!StringCmp(AtomName," N  "))));
/*				(!StringCmp(AtomName," H  ")) ||*/				
					else if ((!StringCmp(AtomHere," CA ")) &&
				((!StringCmp(AtomName," H  ")) || ((
				(!StringCmp(AtomName," CD ")) ||
				(!StringCmp(AtomName,"1HD ")) ||
				(!StringCmp(AtomName,"2HD "))) && reshere=='P') ||
				(!StringCmp(AtomName," N  "))));
							else if (((!StringCmp(AtomHere," HA ")) || (!StringCmp(AtomHere,"1HA "))) &&
				reshere=='P' &&
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
				(!StringCmp(AtomName," CD "))) && reshere=='P') ||
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
							else if (rescrash=='P' && (
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
				rescrash=='P') ||
				(!StringCmp(AtomHere," N  "))));
							else if (((!StringCmp(AtomName," HA ")) ||
				(!StringCmp(AtomName,"1HA "))) &&
				rescrash=='P' &&
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
				rescrash=='P') ||
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
							else if (reshere=='P' && (
				(!StringCmp(AtomHere," N  ")) ||
				(!StringCmp(AtomHere," N  "))
)) /*printf("excused\n")*/;
							else if ((!StringCmp(AtomHere," CA ")) && (!StringCmp(AtomName," CA ")));
							else doit=1;
						}
						else doit=1;
					
						if (doit==1) {						
							/* atoms too close, unless... check for H-bonds */
							if ((((((AtomName[1]=='O') || (AtomName[1]=='N')) && ((AtomHere[1]=='O') || (AtomHere[1]=='N'))))) && (FindHBond(pmadHere,(PMAD)(paldHere->pfbParent),TRUE)==TRUE)) {
								/* possible H-bond so get smaller radii */
								if (AtomName[1]=='N') {
									if (at1hn==31)
										radAtom1=VDWRAD_N3H1;
									else if (at1hn==32)
										radAtom1=VDWRAD_N3H2;
									else /* if (at1hn==43) */
										radAtom1=VDWRAD_N4H3;
								}
								if (AtomName[1]=='O') {
									if (at1hn==10)
										radAtom1=VDWRAD_O1H0;
									else /* if (at1hn==21) */
										radAtom1=VDWRAD_O2H1;
								}
								if (AtomHere[1]=='N') {
									if (at2hn==31)
										radAtom2=VDWRAD_N3H1;
									else if (at2hn==32)
										radAtom2=VDWRAD_N3H2;
									else /* if (at2hn==43) */
										radAtom2=VDWRAD_N4H3;
								}
								if (AtomHere[1]=='O') {
									if (at2hn==10)
										radAtom2=VDWRAD_O1H0;
									else /* if (at2hn==21) */
										radAtom2=VDWRAD_O2H1;
								}														
								tooclose=getMag(vHere)-(radAtom1+radAtom2);
								if (tooclose<0) {
									/* still too close even if an H-bond so use bounciness to check acceptance */
									/* use non-H-bond van der Waals for bump checking */						
									printf("Crash between mol %s res %d - %s and mol %s res %d - %s, too close by %f Angstroms, dist=%f\n",chain1,res1,AtomHere,chain2,res2,AtomName,radAtom1+radAtom2-getMag(vHere),getMag(vHere));												
								}
								/* all hydrogens are placed so check for H-bond */
								if (FindHBond(pmadHere,(PMAD)(paldHere->pfbParent),TRUE)==FALSE) {
									printf("Crash between mol %s res %d - %s and mol %s res %d - %s, too close by %f Angstroms, dist=%f\n",chain1,res1,AtomHere,chain2,res2,AtomName,radAtom1+radAtom2-getMag(vHere),getMag(vHere));												
								}
								else
									printf("Crash between mol %s res %d - %s and mol %s res %d - %s, too close by %f Angstroms, dist=%f\n",chain1,res1,AtomHere,chain2,res2,AtomName,radAtom1+radAtom2-getMag(vHere),getMag(vHere));												
							}
							else {						
								printf("Crash between mol %s res %d - %s and mol %s res %d - %s, too close by %f Angstroms, dist=%f\n",chain1,res1,AtomHere,chain2,res2,AtomName,radAtom1+radAtom2-getMag(vHere),getMag(vHere));
							}
						}	          				
			
		}
                   /* iterate through neighbor list */
			  	vnpHere=vnpHere->next;
     		}
				/* free the atom list created */
	      ValNodeFree(vnpAtomList);
	   	}
			pvnmaHere=pvnmaHere->next;
		}
		pdnmgHere=pdnmgHere->next;
	}
	pdnmmHere=pdnmmHere->next;
}
	/* free the worlds */
	vnpHead=FreeAllWorlds();
	/* Shut Down MMDB-API */
	CloseMMDBAPI();	
	printf("Completed\n");
 	return TRUE;
}


/*  
$Log: bumpcheck.c,v $
Revision 1.10  2001/10/10 19:52:02  feldman
Should now match Foldtraj bumpchecking code exactly

Revision 1.9  2001/10/02 21:57:06  feldman
Updated bumpcheck so it matches FOLDTRAJ almost exactly now, and changed
argument to getRgyr

Revision 1.8  2001/07/26 19:20:45  feldman
Updated for new addtoworld

Revision 1.7  2001/03/14 16:25:53  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.6  2001/01/12 20:01:50  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.5  2000/10/25 15:15:12  feldman
Made further updates, multiple model support is correct now
and relocated to a single function for loading correct model
and extensive error handling

Revision 1.4  2000/10/24 20:57:24  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.3  2000/07/14 20:09:24  feldman
Added parameter for MIMEBiostrucAsnGet for Bioseq

Revision 1.2  2000/07/06 15:29:36  feldman
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

