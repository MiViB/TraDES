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


/* routines for generating a Ramachandran plot of a pfb using vibrant */

/* includes */
#include <mmdbtraj.h>

PMAD FindAtomName(PVNMA pvnmaHead,CharPtr AtomName)
{
        PVNMA pvnmaHere;
	PMAD pmadtmp;

        pvnmaHere=pvnmaHead;
        while (pvnmaHere) {
		pmadtmp=(PMAD)(pvnmaHere->data.ptrvalue);
	        if (!StringCmp(AtomName,pmadtmp->pcAName))
			return pmadtmp;
		pvnmaHere=pvnmaHere->next;
	}
	return NULL;
}

/* this takes a ValNodePtr to a singly linked list containing PMSDs, PMMDs
   or PMGDs and creates a new linked list containing their children
   PMMDs, PMGDs or PMADs and returns a pointer to its head; the original
   linked list is then freed */
ValNodePtr getChild(ValNodePtr vnpTop,Byte ParentType)
{
	ValNodePtr vnpHere,vnptarget=NULL,vnplast=NULL;
	PMSD pmsdtarget;
	PDNMM pdnmmtarget;
	PMMD pmmdtarget;
	PDNMG pdnmgtarget;
	PMGD pmgdtarget;
	PVNMA pvnmatarget;
	PMAD pmadtarget;

	/* remember top, then traverse list */
	vnpHere=vnpTop;
	while(vnpHere) {
		switch(ParentType) {
			/* if finding MSD children ... */
			case AM_MSD:
				pmsdtarget=(PMSD)(vnpHere->data.ptrvalue);
				/* add nodes for all MMDs in pmsdtarget */
				if (pmsdtarget != NULL) {
					pdnmmtarget=pmsdtarget->pdnmmHead;
					while (pdnmmtarget) {
						pmmdtarget=(PMMD)(pdnmmtarget->data.ptrvalue);
						if (pmmdtarget!=NULL) {
							vnplast=ValNodeAddPointer(&vnplast,0,(PFB)(pmmdtarget));
							if (vnptarget==NULL)
								vnptarget=vnplast;
						}
						pdnmmtarget=pdnmmtarget->next;
					}		
				}
				break;
			/* if finding MMD children ... */
			case AM_MMD:
				pmmdtarget=(PMMD)(vnpHere->data.ptrvalue);
				/* add nodes for all MGDs in pmmdtarget */
				if (pmmdtarget != NULL) {
					pdnmgtarget=pmmdtarget->pdnmgHead;
					while (pdnmgtarget) {
						pmgdtarget=(PMGD)(pdnmgtarget->data.ptrvalue);
						if (pmgdtarget!=NULL) {
							vnplast=ValNodeAddPointer(&vnplast,0,(PFB)(pmgdtarget));
							if (vnptarget==NULL)
								vnptarget=vnplast;
						}
						pdnmgtarget=pdnmgtarget->next;
					}		
				}
				break;
			/* if finding MGD children ... */
			case AM_MGD:
				pmgdtarget=(PMGD)(vnpHere->data.ptrvalue);
				/* add nodes for all MADs in pmsdtarget */
				if (pmgdtarget != NULL) {
					pvnmatarget=pmgdtarget->pvnmaAHead;
					while (pvnmatarget) {
						pmadtarget=(PMAD)(pvnmatarget->data.ptrvalue);
						if (pmadtarget!=NULL) {	
							vnplast=ValNodeAddPointer(&vnplast,0,(PFB)(pmadtarget));
							if (vnptarget==NULL)
								vnptarget=vnplast;
						}	
						pvnmatarget=pvnmatarget->next;
					}
				}
		}
		vnpHere=vnpHere->next;
	} 
	/* free original list and return new one */
	vnpTop=ValNodeFree(vnpTop);
	return vnptarget;	
}

/* this take any PFB of type PMSD, PMMD, PMGD PMBD or PMAD and generates a
   linked list, pointed to by the returned ValNodePtr, containing
   PMSDs, PMMDs, PMGDs PMBDs or PMADs (speicified by toNodeType) related
   to that input PFB; e.g. if pfbfrom is a PMMD and toNodeType is
   AM_MAD then the return list will contain all PMADs within the
   molecule pointed to by pfbfrom */ 
ValNodePtr ConvertNode(PFB pfbfrom,Byte toNodeType)
{
	Byte fromNodeType;
	ValNodePtr vnpHead=NULL,vnpBondHead=NULL,vnpHere,vnplast=NULL,vnpBlast=NULL;
	PMSD pmsdtarget;
	PMMD pmmdtarget;
	PMGD pmgdtarget;
	PMBD pmbdtarget;
	PVNMB pvnmbList;
	PFB pfb1,pfb2;
	

	if (pfbfrom==NULL) return NULL;
	fromNodeType=pfbfrom->bMe;
	/* resolve trivial case, same types */
	if (fromNodeType & toNodeType) {
		ValNodeAddPointer(&vnpHead,0,pfbfrom);
		return vnpHead;		
	}
	/* take care of converting from bonds */
	if (fromNodeType & AM_MBD) {
		pfb1=(PFB)(((PMBD)pfbfrom)->pmadFrom);
		pfb2=(PFB)(((PMBD)pfbfrom)->pmadTo);
		/* getparentgraph, etc cal handle NULL argument */
		switch (toNodeType) {
			case AM_MSD:
				pfb1=(PFB)ToMSDParent(pfb1);
				pfb2=(PFB)ToMSDParent(pfb2);
			case AM_MMD:
				pfb1=(PFB)GetParentMol(pfb1);
				pfb2=(PFB)GetParentMol(pfb2);
			case AM_MGD:
				pfb1=(PFB)GetParentGraph(pfb1);
				pfb2=(PFB)GetParentGraph(pfb2);
		}
		if (pfb1!=NULL)
			ValNodeAddPointer(&vnpHead,0,(PFB)pfb1);
		/* if both atoms belong to same molecule or whatever,
		   only return the first instance */
		if ((pfb2!=NULL) && (pfb1!=pfb2))
			ValNodeAddPointer(&vnpHead,0,(PFB)pfb2);
		return vnpHead;
	}
	/* note converting from bonds is now taken care of */
	/* if want to convert to MSD, can use toMSDParent */
	if (toNodeType & AM_MSD) {
		pmsdtarget=ToMSDParent(pfbfrom);
		ValNodeAddPointer(&vnpHead,0,(PFB)pmsdtarget);
		return vnpHead;		
	}
	/* if want to convert to MMD and is lower in hierarchy, use 
	   GetParentMol */	
	if ((toNodeType & AM_MMD) && (fromNodeType!=AM_MSD)) {
		pmmdtarget=GetParentMol(pfbfrom);
		ValNodeAddPointer(&vnpHead,0,(PFB)pmmdtarget);
		return vnpHead;		
	}
	/* if want to convert to MGD and lower in hierarchy, use
	   GetParentGraph */
	if ((toNodeType & AM_MGD) && (fromNodeType & AM_MAD)) {
		pmgdtarget=GetParentGraph(pfbfrom);
		ValNodeAddPointer(&vnpHead,0,(PFB)pmgdtarget);
		return vnpHead;		
	}
	/* otherwise make a node for the parent */
	ValNodeAddPointer(&vnpHead,0,pfbfrom);
	/* and use getChild to "recursively" replace each level with
	   the contents one level down until the right level is reached */
	switch (fromNodeType) {
		case AM_MSD:
			/* if converting to bonds, add inter-molecular
	   		bonds to vnpBondHead list */
			if (toNodeType & AM_MBD) {
				/* go through MSD list */
				vnpHere=vnpHead;
				while(vnpHere) {
					pmsdtarget=(PMSD)(vnpHere->data.ptrvalue);
					if (pmsdtarget!=NULL) {
						/* and add bonds in bond
						   list */
						pvnmbList=pmsdtarget->pvnmbIMBHead;
						while(pvnmbList) {
							pmbdtarget=(PMBD)(pvnmbList->data.ptrvalue);
							if (pmbdtarget!=NULL) {
								vnpBlast=ValNodeAddPointer(&vnpBlast,0,(PFB)(pmbdtarget));
								if (vnpBondHead==NULL)
									vnpBondHead=vnpBlast;
							}
							pvnmbList=pvnmbList->next;
						}
					}
					vnpHere=vnpHere->next;
				}
			}
			vnpHead=getChild(vnpHead,AM_MSD);
			if (toNodeType & AM_MMD) return vnpHead;
		case AM_MMD:
			/* if converting to bonds, add inter-residue
			   bonds to vnpBondHead list */
			if (toNodeType & AM_MBD) {
				/* go through MMD list */
				vnpHere=vnpHead;
				while(vnpHere) {
					pmmdtarget=(PMMD)(vnpHere->data.ptrvalue);
					if (pmmdtarget!=NULL) {
						/* and add bonds in bond
						   list */
						pvnmbList=pmmdtarget->pvnmbIRBHead;
						while(pvnmbList) {
							pmbdtarget=(PMBD)(pvnmbList->data.ptrvalue);
							if (pmbdtarget!=NULL) {
								vnpBlast=ValNodeAddPointer(&vnpBlast,0,(PFB)(pmbdtarget));
								if (vnpBondHead==NULL)
									vnpBondHead=vnpBlast;
							}
							pvnmbList=pvnmbList->next;
						}
					}
					vnpHere=vnpHere->next;
				}
			}
			vnpHead=getChild(vnpHead,AM_MMD);
			if (toNodeType & AM_MGD) return vnpHead;
		case AM_MGD:
			/* if converting to bonds, add graph
			   bonds to vnpHead list */
			if (toNodeType & AM_MBD) {
				/* go through MGD list */
				vnpHere=vnpHead;
				while(vnpHere) {
					pmgdtarget=(PMGD)(vnpHere->data.ptrvalue);
					if (pmgdtarget!=NULL) {
						/* and add bonds in bond
						   list */
						pvnmbList=pmgdtarget->pvnmbBHead;
						while(pvnmbList) {
							pmbdtarget=(PMBD)(pvnmbList->data.ptrvalue);
							if (pmbdtarget!=NULL) {
								vnplast=ValNodeAddPointer(&vnplast,0,(PFB)(pmbdtarget));
								if (vnpHead==NULL)
									vnpHead=vnplast;
							}
							pvnmbList=pvnmbList->next;
						}
					}
					vnpHere=vnpHere->next;
				}
			}
			/* only find atoms if to type isn't bonds */
			else vnpHead=getChild(vnpHead,AM_MGD);
			/* append vnpBondHead to vnpHead - ok if NULL */
			vnpHere=vnpHead;
			if (vnpHere!=NULL) {
				while(vnpHere->next)
					vnpHere=vnpHere->next;
				vnpHere->next=vnpBondHead;
			}
			if (toNodeType & (AM_MAD | AM_MBD)) return vnpHead;
	}
	/* in case of problems, return NULL */
	return NULL;
}

/* gives curvec the value of the vector between fromAtom and ToAtom in Model
   in Euclidean space, returning 1 if successful, 0 otherwise */
Int2 setVec(vec curvec,PMAD fromAtom,PMAD toAtom,Int2 Model)
{
	PALD Loc,Loc2;
	PMGD pmgdErr;

	Loc=GetAtomLocs(toAtom,Model);
	Loc2=GetAtomLocs(fromAtom,Model);
	/* ensure atoms have co-ordinates */
	if (Loc==NULL) {
		pmgdErr=GetParentGraph((PFB)toAtom);
		ErrPostEx(SEV_INFO,1,1,"Warning: unable to find co-ordinates for atoms in %s %d, model %d",
		 pmgdErr->pcGraphName,(pmgdErr->pdnmgLink)->choice,Model);
 		return 0;
	}	  
	if (Loc2==NULL) {
		pmgdErr=GetParentGraph((PFB)fromAtom);
		ErrPostEx(SEV_INFO,1,1,"Warning: unable to find co-ordinates for atoms in %s %d, model %d",
		 pmgdErr->pcGraphName,(pmgdErr->pdnmgLink)->choice,Model);
 		return 0;
	}
	/* set value */
	curvec[0]=AtomLocX(Loc2)-AtomLocX(Loc);
	curvec[1]=AtomLocY(Loc2)-AtomLocY(Loc);
	curvec[2]=AtomLocZ(Loc2)-AtomLocZ(Loc);
	/* return different value if multiple occupancies */
	if(Loc->next || Loc2->next)
		return 2;
	return 1;
}

/* Find a non-alpha-carbon bonded neighbour to curAtom, with
	atomic number anum, on or off the Backbone */
PMAD findNeighbour(PMAD curAtom,Uint1 anum,Boolean Backbone)
{
	PVNMB curBond;
	PMAD fromAtom,toAtom;
	PMBD pmbdBond;
	PVNMA pvnma;

	curBond=curAtom->pvnBonds;
	/* if curBond is NULL, will skip while loop */
	while (curBond) {
		pmbdBond=(PMBD)(curBond->data.ptrvalue);
		/* break out of loop if error occurs at a given bond */
		if (pmbdBond==NULL)
			continue;
		fromAtom=pmbdBond->pmadFrom;
		toAtom=pmbdBond->pmadTo;
		/* check if fromAtom is the one we want */
		if (fromAtom!=NULL) {
			pvnma=fromAtom->pvnmaLink;
			if (pvnma!=NULL)
				if (((pvnma->choice)==anum) &&
		 		 (IsAtomBackBone((PFB)fromAtom)==Backbone) &&
		 		 (!IsAtomCAlpha((PFB)fromAtom))) 
				 return(fromAtom);
		}
		/* if not, maybe toAtom is the one */
		if (toAtom!=NULL) {
			pvnma=toAtom->pvnmaLink;
			if (pvnma!=NULL)
				if (((pvnma->choice)==anum) &&
		 		 (IsAtomBackBone((PFB)toAtom)==Backbone) &&
		 		 (!IsAtomCAlpha((PFB)toAtom)))
				 return(toAtom);
		}
		curBond=curBond->next;
	}
	/* nothing found :( */
	return (NULL);
}

/* calculates phi and psi angles at a given alpha carbon curpmad, in a
   given Model and stores at end of linked list pointed to by prsHead
   NOTE: only returns values for first atom location if multiple
   occupancies exist; returns INCOMPLETE if unable to determine phi, psi */
TrajErr getR(PMAD curpmad,PRS PNTR prsHead,Int2 Model)
{
	PMAD CAlpha,Ci,NiPlusOne=NULL,Ni,CiMinusOne=NULL,CAlphaPlusOne;
	vec CaC2,CaN1,N1C1,N2Ca,N1C2,N1Ca,C2Ca;
	FloatLo Phi,Psi,Omega,thedot;
	Int2 num;
	Int4 chainlength;
	Char ResID[4];
	vec U,V,u,v,MainAxis; 
	PMGD pmgdParent,pmgdNext;
	PDNMG pdnmgParent;
	PMMD pmmdParent;
	PRS prsNode,prsHere;

    	/* we want only ModelAtomData or MAD nodes, specifically a-carbon */
	if (IsAtomCAlpha((PFB)curpmad)) {
		/* num = residue number, or 0 if not found */
		pmgdParent=GetParentGraph((PFB)curpmad);
		if (pmgdParent==NULL) num=0;
		else {
			pdnmgParent=pmgdParent->pdnmgLink;
			if (pdnmgParent==NULL) num=0;
			else num=pdnmgParent->choice;
		}
		if (!num) {
			ErrPostEx(SEV_ERROR,1,1,"Warning: Unable to determine residue number in Model %d",Model);
			return ERR_FAIL;
		}
		pmmdParent=GetParentMol((PFB)pmgdParent);
		/* find length of entire peptide */
		if (pmmdParent==NULL) {
			ErrPostEx(SEV_ERROR,2,1,"Warning: Unable to determine chain length for Model %d",Model);
			return ERR_FAIL;
		}
		else chainlength=pmmdParent->iResCount;
		/* ensure it's a protein and get residue name */
		if (!IsProtein((PFB)pmmdParent)) return ERR_INCOMPLETE;
		StringCpy(ResID,"XXX");
		/* pmgdParent can't be NULL if we got here */
		StringNCpy(ResID,pmgdParent->pcGraphName,3);
		/* ensure ane is null terminated */
		ResID[3]=0;
		if (num==1 || num==chainlength) {
		/*	ErrPostEx(SEV_INFO,2,1,"Warning: No Phi and Psi at protein ends"); CWVH 2012 dont warn the obvious and clog the log file... */ 
			return ERR_ENDOFPROTEIN;
		}	
		CAlpha=curpmad;
		/* find the other four atoms! */
		Ci=findNeighbour(CAlpha,6,TRUE);
		Ni=findNeighbour(CAlpha,7,TRUE);
		/* skip phi for 1st residue */
		if (Ni==NULL) {
			ErrPostEx(SEV_INFO,3,1,"Warning: unable to find Ni neighbour to residue %d, Model %d",num,Model);
			return ERR_INCOMPLETE;
		}
		else if (num>1) CiMinusOne=findNeighbour(Ni,6,TRUE);
		/* skip psi for last residue */
		if (Ci==NULL) {
			ErrPostEx(SEV_INFO,3,2,"Warning: unable to find Ci neighbour to residue %d, Model %d",num,Model);
			return ERR_INCOMPLETE;
		}
		else if (num<chainlength) NiPlusOne=findNeighbour(Ci,7,TRUE);
		/* set vectors required for calculation and ensure they
		   are valid */
		if (!setVec(CaC2,CAlpha,Ci,Model)) return ERR_INCOMPLETE;
		NegateVec(C2Ca,CaC2);
		if (!setVec(CaN1,CAlpha,Ni,Model)) return ERR_INCOMPLETE;
/* need c1cai-1*/
		if (CiMinusOne==NULL) {
			ErrPostEx(SEV_INFO,3,3,"Warning: unable to find Ci-1 neighbour to residue %d, Model %d",num,Model);
			return ERR_INCOMPLETE;
		}
		else if (num>1) {
			if (!setVec(N1C1,Ni,CiMinusOne,Model)) return ERR_INCOMPLETE;
		}
		if (NiPlusOne==NULL) {
			ErrPostEx(SEV_INFO,3,4,"Warning: unable to find Ni+1 neighbour to residue %d, Model %d",num,Model);
			return ERR_INCOMPLETE;
		}
		else if (num<chainlength) {
			if (!setVec(N1C2,NiPlusOne,Ci,Model)) return ERR_INCOMPLETE;
			/* get last Ca too */
			pmgdNext=(PMGD)(NiPlusOne->pfbParent);
			CAlphaPlusOne=FindCAlpha(pmgdNext->pvnmaAHead);
			if (CAlphaPlusOne==NULL) return ERR_INCOMPLETE;
			if (!setVec(N1Ca,NiPlusOne,CAlphaPlusOne,Model)) return ERR_INCOMPLETE;
		}
/*if (!setVec(N2Ca,NiPlusOne,CAlpha,Model)) return;*/
		if (!setVec(N2Ca,Ni,CAlpha,Model)) return ERR_INCOMPLETE;
/*if (!setVec(N1C2,Ni,Ci,Model)) return;*/
		Phi=0;
		/* calculate phi unless 1st residue */
		if (num>1) { 
			/* front vector for phi, u=projection onto the
			   plane with the main axis as its normal */
			Normalize(MainAxis,CaN1);
			VecScale(U,MainAxis,Dot(CaC2,MainAxis));
			NegateVec(U,U);
			VecAdd(U,CaC2,U);
			Normalize(u,U);
			/* back vector for phi, v=projection onto the 
			   same plane as in u */
			VecScale(V,MainAxis,Dot(N1C1,MainAxis));
			NegateVec(V,V);
			VecAdd(V,N1C1,V);
			Normalize(v,V);
			/* between zero and Pi by default */
			thedot=Dot(u,v);
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			Phi=RADTODEG*acos(thedot);
			/* Determine sign of Phi */
			Cross(U,u,v);
			if (Dot(U,MainAxis)>0)
				Phi=-Phi;
		}
		Omega=0;
		/* calculate Omega (between i and i+1) unless last residue */
		if (num<chainlength) { 
			Normalize(MainAxis,N1C2);
			VecScale(U,MainAxis,Dot(N1Ca,MainAxis));
			NegateVec(U,U);
			VecAdd(U,N1Ca,U);
			Normalize(u,U);
			VecScale(V,MainAxis,Dot(C2Ca,MainAxis));
			NegateVec(V,V);
			VecAdd(V,C2Ca,V);
			Normalize(v,V);
			thedot=Dot(u,v);
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			Omega=RADTODEG*acos(thedot);
			Cross(U,u,v);
			if (Dot(U,MainAxis)>0)
				Omega=-Omega;
		}
		Psi=0;
		/* find psi except for last residue */
		if (num<chainlength) {
			/* front vector for psi -- see above for phi */
			Normalize(MainAxis,CaC2);
			VecScale(U,MainAxis,Dot(N2Ca,MainAxis));
			NegateVec(U,U);
			VecAdd(U,N2Ca,U);
			NegateVec(U,U);
			Normalize(u,U);
			/* back vector for psi -- see above for phi */
			VecScale(V,MainAxis,Dot(N1C2,MainAxis));
			NegateVec(V,V);
			VecAdd(V,N1C2,V);
			NegateVec(V,V);
			Normalize(v,V);
			/* between zero and Pi by default */
			thedot=Dot(u,v);
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			Psi=RADTODEG*acos(thedot);
			/* Determine sign of Psi */
			Cross(U,u,v);
			if (Dot(U,MainAxis)>0)
				Psi=-Psi;
		}
		if ((num>1) && (num<chainlength)) {
			/* allocate a new RS node and assign values */
			prsNode=(PRS)MemNew(sizeof(RS));
			prsNode->next=NULL;
			prsNode->pfbThis=(PFB)pmgdParent;
			prsNode->Phi=Phi;
			prsNode->Psi=Psi;
			prsNode->Omega=Omega;
			/* if first node, save in prsHead */
			if (*prsHead==NULL)
				*prsHead=prsNode;
			/* otherwise go to end of list and append it */
			else {
				prsHere=*prsHead;
				while(prsHere->next)
					prsHere=prsHere->next;
				prsHere->next=prsNode;
			}
		}
	}
	return ERR_SUCCESS;
}

/* this frees a RS structure pointed to by prsHead and returns NULL */
PRS freeRS(PRS prsHead)
{
        PRS prsNext;
   
	/* walk down list deleting one node at a time */
        while (prsHead) {
                prsNext=prsHead->next; 
                MemFree(prsHead);
                prsHead=prsNext;
        }
        return NULL;
}

/* Rama takes a singly linked list pointed to by vnpHere, and returns
   a pointer to a data structure (PRS) containing the phi and psi angles
   for all a-carbons contained within all structures/molecules/residues/atoms
   within the linked-list, in a given Model (see Rplot.h for typedef) */
PRS Rama(ValNodePtr vnpHere,Int2 Model)
{
	PFB pfbHere;
	ValNodePtr vnpRamaHere,vnpRamaHead;
	PMGD pmgdRamaIt;
	PVNMA pvnmaAtoms;
	PMAD pmadAtom;
	PRS prsHead=NULL,prsNode,prsHere;
	TrajErr err;

	/* go through list of objects */
	while(vnpHere) {
		pfbHere=(PFB)(vnpHere->data.ptrvalue);
		if (pfbHere!=NULL) {
			/* convert current object to MGD list if not already */
			vnpRamaHead=ConvertNode(pfbHere,AM_MGD);
			vnpRamaHere=vnpRamaHead;
			/* step through MGD list one at a time */
			while(vnpRamaHere) {
				pmgdRamaIt=(PMGD)(vnpRamaHere->data.ptrvalue);
				if (pmgdRamaIt!=NULL) {
					pvnmaAtoms=pmgdRamaIt->pvnmaAHead;
					/* go through atoms in each MGD to pick
					   out a-carbon(s) */
					while(pvnmaAtoms) {
						pmadAtom=(PMAD)(pvnmaAtoms->data.ptrvalue);
						/* get phi, psi for a-carbons */
						if (pmadAtom!=NULL)
							if (IsAtomCAlpha((PFB)pmadAtom)) {
								err=getR(pmadAtom,&prsHead,Model);			
								if (err==ERR_INCOMPLETE) {			
									prsNode=(PRS)MemNew(sizeof(RS));
									prsNode->next=NULL;
									prsNode->pfbThis=(PFB)pmgdRamaIt;
									prsNode->Phi=0.0;
									prsNode->Psi=0.0;
									prsNode->Omega=0.0;
									prsNode->Mag=VL_MISSCALPHA;
									/* if first node, save in prsHead */
									if (prsHead==NULL)
										prsHead=prsNode;
									/* otherwise go to end of list and append it */
									else {
										prsHere=prsHead;
										while (prsHere->next)
											prsHere=prsHere->next;
										prsHere->next=prsNode;
									}
						    }
							}
						pvnmaAtoms=pvnmaAtoms->next;
					}
				}
				vnpRamaHere=vnpRamaHere->next;
			}
			/* free up the temporary MGD list */
			ValNodeFree(vnpRamaHead);
		}
		vnpHere=vnpHere->next;
	}
	return prsHead;
}


/*  
$Log: Rplot.c,v $
Revision 1.6  2001/07/24 17:27:12  feldman
Updated rotamer library version and fixed bug - order of building
for Ramachandran space had prevented rotamer library from being
used properly

Revision 1.5  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.4  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.3  2000/12/15 00:06:41  feldman
Now Val2Trj and unfoldtraj works correctly for Phi-Psi walk (I think)

Revision 1.2  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

