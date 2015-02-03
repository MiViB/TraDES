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


/* note that bit 0 of bReserved in PMADs is used by these functions to
   keep track of which atoms have been previously visited */ 

#include <mmdbtraj.h>

#define ROTATION_STEPSIZE 40.0          /* degrees */
#define ALIGN_PRECISION 0.5             /* degrees */
#define MAX_STEP 30
#define PYTHAG(a, b) ((at = fabs(a)) > (bt = fabs(b)) ? \
    (ct = bt / at, at * sqrt(1.0 + ct * ct)) : \
    (bt ? (ct = at / bt, bt * sqrt(1.0 + ct * ct)) : 0.0))
#define SV_MAX(a, b) (maxarg1 = (a), maxarg2 = (b), (maxarg1)>(maxarg2) ? \
    (maxarg1) : (maxarg2))
#define SV_SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* global variables */
static vec vTrans,vnegTrans;
static FloatLo d,sina,cosa,cosTheta,sinTheta,axisx;
static FloatLo maxarg1,maxarg2,at,bt,ct;
static vec vCofM_PutMeBack;

/* finds co-ordinates of AtomName in pmgdHere relative to vCAlpha, and stores
   result in vDest; normally vCAlpha is the CA co-ordinate in pmgdHere */
PALD GetCoOrds(PMGD pmgdHere,CharPtr AtomName,vec vCAlpha,vec vDest,Int2 Model)
{
        PMAD pmadThis;
        PALD paldHere;
        vec vTmp;

        pmadThis=FindAtomName(pmgdHere->pvnmaAHead,AtomName);
        if (pmadThis==NULL) return NULL;
        paldHere=GetAtomLocs(pmadThis,Model);
        if (paldHere==NULL) return NULL;
        vTmp[0]=AtomLocX(paldHere);
        vTmp[1]=AtomLocY(paldHere);
        vTmp[2]=AtomLocZ(paldHere);
        VecSub(vDest,vTmp,vCAlpha);
        return paldHere;
}

static void LIBCALLBACK clearbit(PFB curpfb,Int4 Model,Int4 dum1,Pointer dum2)
{
        (curpfb->bReserved)&=0xfe;
}

/* clears bReserved bit status in given modelnum for given pdnmmThis */
void ClearUpdate(PDNMM pdnmmThis,Int2 modelnum)
{
        TraverseOneModel(pdnmmThis,TRAVERSE_ATOM,modelnum,0,NULL,(pNodeFunc)(clearbit));
}

static void LIBCALLBACK Nlm_setbit(PFB curpfb,Int4 Model,Int4 dum1,Pointer dum2)
{
        (curpfb->bReserved)|=0x01;
}

/* clears bReserved bit status in given modelnum for given pdnmmThis */
void SetUpdate(PDNMM pdnmmThis,Int2 modelnum)
{
        TraverseOneModel(pdnmmThis,TRAVERSE_ATOM,modelnum,0,NULL,(pNodeFunc)(Nlm_setbit));
}

/* determines whether pmad1 or pmad2 is closer to pmadTmp, where all three
   atoms must lie on the same residue and one of pmad1 and pmad2 must lie
   between the other and pmadTmp on the chemical graph of the residue */
PMAD FindClosest(PMAD pmadTmp,PMAD pmad1,PMAD pmad2) 
{
        ValNodePtr pvnBList;
        PMAD pmadFound,pmadNew;
        PMBD pmbdHere;

        pvnBList=pmadTmp->pvnBonds;
        while (pvnBList) {
                pmbdHere=(PMBD)(pvnBList->data.ptrvalue);
                pmadFound=pmbdHere->pmadFrom;
                if (pmadFound==pmadTmp)
                        pmadFound=pmbdHere->pmadTo;
                if (pmadFound==pmad1) return pmad1;             
                if (pmadFound==pmad2) return pmad2;
                pmadNew=FindClosest(pmadFound,pmad1,pmad2);
                if (pmadNew!=NULL) return pmadNew;              
                pvnBList=pvnBList->next;
        }
        return NULL;
}

/* rotates the specified pmadHere (all its PALD's) in a given model by the
   parameters specified in the global variables cosTheta, sinTheta, d,
   etc. */
void RotateAtom(PMAD pmadHere,Int2 Model)
{
        PALD paldHere;
        vec vCur;

        /* do for all PALD's, if any exist */
        paldHere=GetAtomLocs(pmadHere,Model);
        while(paldHere) {
                vCur[0]=AtomLocX(paldHere);
                vCur[1]=AtomLocY(paldHere);
                vCur[2]=AtomLocZ(paldHere);
                /* must translate -vTrans to origin */
                Translate(vCur,vnegTrans);
                /* then map vector to the z-axis */
                Rotatecos(X_AXIS,sina,cosa,vCur,vCur);
                Rotatecos(Y_AXIS,-axisx,d,vCur,vCur);
                /* rotate Theta about the new z-axis */
                Rotatecos(Z_AXIS,sinTheta,cosTheta,vCur,vCur);
                /* return vector to its original position */
                Rotatecos(Y_AXIS,axisx,d,vCur,vCur);
                Rotatecos(X_AXIS,-sina,cosa,vCur,vCur);
                /* and undo the translation */
                Translate(vCur,vTrans);
                /* save the new atom co-ordinates */
                (paldHere->pflvData)[0]=vCur[0];
                (paldHere->pflvData)[1]=vCur[1];
                (paldHere->pflvData)[2]=vCur[2];
                paldHere=paldHere->next;
        }
}

/* calls RotateAtom on all atoms in a given model on the side of pmbdStart
   which contains pmadStart; note the global variable d, cosa, etc. must
   be set before calling this function since they are the same for each
   recursice call of this function */
void RotateAtoms(PMBD pmbdStart,PMAD pmadStart,Int2 Model)
{
        PMAD pmadNext;
        ValNodePtr vnpHere;
        PMBD pmbdHere;
        Byte btype;

        /* rotate the current atom and set if to "moved" */
        RotateAtom(pmadStart,Model);
        (pmadStart->bReserved)|=0x01;
        vnpHere=pmadStart->pvnBonds;
        /* get list of bonds at current location */
        while(vnpHere) {
                pmbdHere=(PMBD)(vnpHere->data.ptrvalue);
                btype=pmbdHere->bWhat;
                /* do not follow bond we just came on, and only do single,
                   double, triple and aromatic bonds */
                if ((pmbdHere!=pmbdStart) && (btype & (BOND_SINGLE|BOND_DOUBLE|BOND_PDOUBLE|BOND_TRIPLE))) {
                        /* find atom at other end of the bond */
                        pmadNext=pmbdHere->pmadTo;
                        if (pmadNext==pmadStart)
                                pmadNext=pmbdHere->pmadFrom;
                        /* and if we haven't "moved" it yet, do so
                           recursively */
                        if (((pmadNext->bReserved)&0x01)==0)
                                RotateAtoms(pmbdHere,pmadNext,Model);
                }
                /* otherwise continue for remaining bonds */
                vnpHere=vnpHere->next;
        }
}

/* rotates a molecule about the bond axis of the given PMBD (1st
   co-ordinate set if more than one exists for that bond) by an
   angle Theta in radians in a given model; pmadDir, which must
   be one of the atoms at the end of the bond, designates the
   positive sense of rotation about the bond as per the right hand
   rule; the co-ordinates on one side of the bond remain fixed while the
   other will rotate; all co-ordinates for each atom are rotated; if the
   specified bond lies on a closed ring, the function will not work
   (since it makes no sense to rotate about such a bond); the function
   automatically chooses which half to keep fixed and which half to
   actually change the co-ordinates of (the smaller half) */
void RotateMolecule(PMBD pmbdAxis,PMAD pmadDir,FloatLo Theta,Int2 Model)
{
        /* Note right now, can't tell it rotates +Theta or -Theta */
        PMAD pmadFrom,pmadTo,pmadTmp;
        PALD paldFrom,paldTo;
        PMMD pmmdFrom,pmmdTo;
        PVNMA pvnmaList;
        Int4 lengFrom,lengTo,curPos;
        vec vAxis;
        Uint1 ElemFrom,ElemTo;
        PDNMM pdnmmHere;

        /* we want pmadTo to be pmadDir and pmadFrom the other atom */
        pmadFrom=pmbdAxis->pmadFrom;
        pmadTo=pmbdAxis->pmadTo;
        if (pmadFrom==pmadDir) {
                pmadTo=pmadFrom;
                pmadFrom=pmbdAxis->pmadTo;
        }
        ElemTo=(pmadTo->pvnmaLink)->choice;
        ElemFrom=(pmadFrom->pvnmaLink)->choice;
        /* nothing to rotate if bond contains a hydrogen */
        if ((ElemTo==1) || (ElemFrom==1))
                return;
        pmmdFrom=GetParentMol((PFB)pmadFrom);
        pmmdTo=GetParentMol((PFB)pmadTo);
        lengFrom=pmmdFrom->iResCount;
        lengTo=pmmdTo->iResCount;
        curPos=2*(((GetParentGraph((PFB)pmadTo))->pdnmgLink)->choice);
        /* if intermolecular bond, want toAtom on smaller molecule */
        if (lengTo>lengFrom) {
                pmadTmp=pmadFrom;
                pmadFrom=pmadTo;
                pmadTo=pmadTmp;
                Theta=-Theta;
        }
        /* must be backbone c=o or c-o(-h) so want o as toAtom */
        else if (ElemFrom==8) {
                pmadTmp=pmadFrom;
                pmadFrom=pmadTo;
                pmadTo=pmadTmp;
                Theta=-Theta;
        }
        /* do nothing if to=O */
        else if (ElemTo==8);
        else if (IsAtomBackBone((PFB)pmadFrom)) {
                if (IsAtomBackBone((PFB)pmadTo)) {
                        /* 2 backbone atoms, non-hydrogen/oxygen */
                        /* peptide or n-c-alpha bond */
                        if (ElemFrom==7) {
                                /* n-c-alpha bond */
                                if (IsAtomCAlpha((PFB)pmadTo)) {
                                        if (curPos<lengFrom) {
                                                pmadTmp=pmadFrom;
                                                pmadFrom=pmadTo;
                                                pmadTo=pmadTmp;
                                                Theta=-Theta;
                                        }
                                }                               
                                /* peptide bond, from=N, to=C */
                                else {
                                        if (curPos>lengFrom) {
                                                pmadTmp=pmadFrom;
                                                pmadFrom=pmadTo;
                                                pmadTo=pmadTmp;
                                                Theta=-Theta;
                                        }
                                }
                        }
                        /* peptide or n-c-alpha bond */
                        else if (ElemTo==7) {
                                /* n-c-alpha bond */
                                if (IsAtomCAlpha((PFB)pmadTo)) {
                                        if (curPos>lengFrom) {
                                                pmadTmp=pmadFrom;
                                                pmadFrom=pmadTo;
                                                pmadTo=pmadTmp;
                                                Theta=-Theta;
                                        }
                                }
                                /* peptide bond */
                                else {
                                        if (curPos<lengFrom) {
                                                pmadTmp=pmadFrom;
                                                pmadFrom=pmadTo;
                                                pmadTo=pmadTmp;
                                                Theta=-Theta;
                                        }
                                }
                        }
                        /* c-alpha c bond */
                        else {
                                if (IsAtomCAlpha((PFB)pmadTo)) {
                                        if (curPos>lengFrom) {
                                                pmadTmp=pmadFrom;
                                                pmadFrom=pmadTo;
                                                pmadTo=pmadTmp;
                                                Theta=-Theta;
                                        }
                                }
                                else {
                                        if (curPos<lengFrom) {
                                                pmadTmp=pmadFrom;
                                                pmadFrom=pmadTo;
                                                pmadTo=pmadTmp;
                                                Theta=-Theta;
                                        }
                                }
                        }
                }
                /* else fromAtom is backbone, toAtom is not, so do not
                   swap */
        }
        else {
                /* from residue to backbone, so swap necessary */
                if (IsAtomBackBone((PFB)pmadTo)) {
                        pmadTmp=pmadFrom;
                        pmadFrom=pmadTo;
                        pmadTo=pmadTmp;
                        Theta=-Theta;
                }       
                /* both atoms are in a residue */
                else {
                        /* follow to see if need to swap or not */ 
                        pvnmaList=GetParentGraph((PFB)pmadTo)->pvnmaAHead;
                        while (!(IsAtomCBeta((PFB)(pvnmaList->data.ptrvalue))))
                                pvnmaList=pvnmaList->next;
                        pmadTmp=(PMAD)(pvnmaList->data.ptrvalue);
                        if (FindClosest(pmadTmp,pmadFrom,pmadTo)==pmadTo) {
                                pmadTmp=pmadFrom;
                                pmadFrom=pmadTo;
                                pmadTo=pmadTmp;
                                Theta=-Theta;
                        }
                }
        }
        cosTheta=cos(Theta);
        sinTheta=sin(Theta);
        paldFrom=GetAtomLocs(pmadFrom,Model);   
        paldTo=GetAtomLocs(pmadTo,Model);
        /* translation required to origin */
        vTrans[0]=AtomLocX(paldFrom);
        vTrans[1]=AtomLocY(paldFrom);
        vTrans[2]=AtomLocZ(paldFrom);
        NegateVec(vnegTrans,vTrans);
        /* find the bond vector (using 1st set of co-ordinates, if more
           than one exist) */
        vAxis[0]=AtomLocX(paldTo)-AtomLocX(paldFrom);
        vAxis[1]=AtomLocY(paldTo)-AtomLocY(paldFrom);
        vAxis[2]=AtomLocZ(paldTo)-AtomLocZ(paldFrom);
        Normalize(vAxis,vAxis);
        /* calculate fixed parameters needed to rotate atoms */
        d=sqrt(vAxis[1]*vAxis[1]+vAxis[2]*vAxis[2]);
        sina=vAxis[1]/d;
        cosa=vAxis[2]/d;
        axisx=vAxis[0];
        /* rotate the atoms now */
        RotateAtoms(pmbdAxis,pmadTo,Model);
        /* clear updated bReserved bit 0 status for all molecules so
           procedure can be safely called multiple times */
        pdnmmHere=(ToMSDParent((PFB)pmadTo))->pdnmmHead;
        ClearUpdate(pdnmmHere,Model);
}

PMAD FindCAlpha(PVNMA pvnmaHead)
{
        PVNMA pvnmaHere;

        if (pvnmaHead==NULL)
                return NULL;
        pvnmaHere=pvnmaHead;
        while (!(IsAtomCAlpha((PFB)(pvnmaHere->data.ptrvalue)))) {
                pvnmaHere=pvnmaHere->next;
                if (pvnmaHere==NULL)
                        return NULL;
        }
        return (PMAD)(pvnmaHere->data.ptrvalue);
}

PMAD FindCBeta(PVNMA pvnmaHead)
{
        PVNMA pvnmaHere;

        if (pvnmaHead==NULL)
                return NULL;
        pvnmaHere=pvnmaHead;
        while (!(IsAtomCBeta((PFB)(pvnmaHere->data.ptrvalue)))) {
                pvnmaHere=pvnmaHere->next;
                if (pvnmaHere==NULL)
                        return NULL;
        }
        return (PMAD)(pvnmaHere->data.ptrvalue);
}

/* translate center of mass of atoms indicated by atomtypes to
   origin; atomtypes are ALIGN_XXX where XXX=CA, BACKBONE or
   ALLATOM (calculates real center of mass, doesn't assume
   all atom types are equal weight!!) - ignores residues outside
   the closed interval [res1,res2] */
   
/* CWVH Dec 08 - modified to put vector values into  vCofM_PutMeBack. */


Int2 CenterMol(PMMD pmmdHere,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model)
{
        PMGD pmgdHere;
        PMAD pmadHere;
        PALD paldHere;
        PMMD pmmdThis;
        PMSD pmsdHere;
        PDNMM pdnmmHere;
        PDNMG pdnmgHere;
        PVNMA pvnmaHere;
        vec vHere,vCofM;
        vec vZero={0,0,0};
        Int2 attype;
        FloatLo totmass;
	PVNMO pvnmoThis;
	PMOD pmodThis;
	int i, j, k; 

        pdnmgHere=pmmdHere->pdnmgHead;
        VecScale(vCofM,vZero,1.0);
        totmass=0.0;
        if (pdnmgHere==NULL) return 0;
        while (pdnmgHere) {
                if ((pdnmgHere->choice>=res1) && (pdnmgHere->choice<=res2)) {
                        pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                        if (pmgdHere==NULL) return 0;
                        pvnmaHere=pmgdHere->pvnmaAHead;
                        if (pvnmaHere==NULL) return 0;
                        while (pvnmaHere) {
                                pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
                                if (pmadHere==NULL) return 0;
                                attype=pmadHere->bWhat;
                                if (((atomtypes&attype)!=0) || (atomtypes==ALIGN_ALLATOM)) {
                                        paldHere=GetAtomLocs(pmadHere,Model);
                                        /* some atoms may have no co-ords,
                                           e.g. HE1 on Glu */
                                        if (paldHere!=NULL) {
                                                vHere[0]=(paldHere->pflvData)[0];
                                                vHere[1]=(paldHere->pflvData)[1];
                                                vHere[2]=(paldHere->pflvData)[2];
                                                switch (pvnmaHere->choice) {
                                                        case 1:
                                                                totmass+=AW_H;
                                                                VecScale(vHere,vHere,AW_H);
                                                                break;
                                                        case 6:
                                                                totmass+=AW_C;
                                                                VecScale(vHere,vHere,AW_C);
                                                                break;
                                                        case 7:
                                                                totmass+=AW_N;
                                                                VecScale(vHere,vHere,AW_N);
                                                                break;
                                                        case 8:
                                                                totmass+=AW_O;
                                                                VecScale(vHere,vHere,AW_O);
                                                                break;
                                                        case 15:
                                                                totmass+=AW_P;
                                                                VecScale(vHere,vHere,AW_P);
                                                                break;
                                                        case 16:
                                                                totmass+=AW_S;
                                                                VecScale(vHere,vHere,AW_S);
                                                                break;
                                                        case 34:
                                                                totmass+=AW_SE;
                                                                VecScale(vHere,vHere,AW_SE);
                                                                break;
                                                        default:
                                                                totmass+=AW_C;
                                                                VecScale(vHere,vHere,AW_C);
                                                }
                                                VecAdd(vCofM,vCofM,vHere);
                                        }
                                }
                                pvnmaHere=pvnmaHere->next;
                        }
                }
                pdnmgHere=pdnmgHere->next;
        }
        VecScale(vCofM,vCofM,-1.0/totmass);
        /* now vCofM contains negative of centre of mass vector */
        vCofM_PutMeBack[0]= -(vCofM[0]); /* save this for putting everything back where it came from or so help me. */
        vCofM_PutMeBack[1]= -(vCofM[1]);
        vCofM_PutMeBack[2]= -(vCofM[2]);
         /* so traverse all atoms and translate appropriately */
         /*CWVH modified to translate all atoms in all molecules in the structure */
         /* this is so when we align a peptide, the bound globular protein is moved as well */
        pmsdHere = ToMSDParent( (PFB) pmmdHere);
        pdnmmHere = pmsdHere->pdnmmHead;
        while (pdnmmHere) { 
        pmmdThis=(PMMD)(pdnmmHere->data.ptrvalue);
        pdnmgHere=pmmdThis->pdnmgHead;
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                pvnmaHere=pmgdHere->pvnmaAHead;
                while (pvnmaHere) {
                        pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
                        paldHere=GetAtomLocs(pmadHere,Model);
                        if (paldHere!=NULL) {
                                /* first 3 elements of pflvdata are a vec */
                                vHere[0]=(paldHere->pflvData)[0];
                                vHere[1]=(paldHere->pflvData)[1];
                                vHere[2]=(paldHere->pflvData)[2];
                                Translate(vHere,vCofM);
                                (paldHere->pflvData)[0]=vHere[0];
                                (paldHere->pflvData)[1]=vHere[1];
                                (paldHere->pflvData)[2]=vHere[2];
                        }
                        pvnmaHere=pvnmaHere->next;
                }
                pdnmgHere=pdnmgHere->next;
        }
        pdnmmHere=pdnmmHere->next;
        }

/* CWVH 2012 also do the cylinders, bricks, cones.. */

  
        pvnmoThis = pmsdHere->pvnmoHead;
	  while (pvnmoThis != NULL)
	    {
 		pmodThis = (PMOD) pvnmoThis->data.ptrvalue;
		if (pmodThis != NULL) {
 			if (pmodThis->ppflObject != NULL) {
				for (i = 0; i < pmodThis->iCoordNo; i++) {
	                               for (j = 0; j < 3; j++) {
						vHere[j] = (pmodThis->ppflObject)[i][j]; 
					} 
					Translate(vHere,vCofM);
					for (k = 0; k < 3; k++) {  /* copy back translated coordinates */
						(pmodThis->ppflObject)[i][k] = vHere[k]; 
					} 
				}
			}
		}
		pvnmoThis = pvnmoThis->next;
	    }
        



        /* done */
        return 1;
}


Int2 PutMolBackWhereItCameFrom(PMMD pmmdHere,Int2 Model)
/* or so help me.. puts entire structure back - not just molecule  */
{
        PMGD pmgdHere;
        PMAD pmadHere;
        PALD paldHere;
        PMSD pmsdHere;
        PMMD pmmdThis;
        PDNMG pdnmgHere;
        PVNMA pvnmaHere;
        PDNMM pdnmmHere;
        vec vHere;
	PVNMO pvnmoThis;
	PMOD pmodThis;
	int i, j, k; 
        
        if (pmmdHere == NULL) return 0;
        
        pmsdHere = ToMSDParent( (PFB) pmmdHere);
        pdnmmHere = pmsdHere->pdnmmHead;
        while (pdnmmHere) { 
        pmmdThis=(PMMD)(pdnmmHere->data.ptrvalue);
        pdnmgHere=pmmdThis->pdnmgHead;
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                pvnmaHere=pmgdHere->pvnmaAHead;
                while (pvnmaHere) {
                        pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
                        paldHere=GetAtomLocs(pmadHere,Model);
                        if (paldHere!=NULL) {
                                /* first 3 elements of pflvdata are a vec */
                                vHere[0]=(paldHere->pflvData)[0];
                                vHere[1]=(paldHere->pflvData)[1];
                                vHere[2]=(paldHere->pflvData)[2];
                                Translate(vHere,vCofM_PutMeBack);
                                (paldHere->pflvData)[0]=vHere[0];
                                (paldHere->pflvData)[1]=vHere[1];
                                (paldHere->pflvData)[2]=vHere[2];
                        }
                        pvnmaHere=pvnmaHere->next;
                }
                pdnmgHere=pdnmgHere->next;
        }
        pdnmmHere = pdnmmHere->next;
        }
/* CWVH 2012 also do the cylinders, bricks, cones.. */

  
        pvnmoThis = pmsdHere->pvnmoHead;
	  while (pvnmoThis != NULL)
	    {
 		pmodThis = (PMOD) pvnmoThis->data.ptrvalue;
		if (pmodThis != NULL) {
 			if (pmodThis->ppflObject != NULL) {
				for (i = 0; i < pmodThis->iCoordNo; i++) {	/* each i is a vec, fill v1 */
					for (j = 0; j < 3; j++) {
						vHere[j] = (pmodThis->ppflObject)[i][j]; 
					} 
					Translate(vHere,vCofM_PutMeBack);
					for (k = 0; k < 3; k++) {  /* copy back translated coordinates */
						(pmodThis->ppflObject)[i][k] = vHere[k]; 
					} 
				}
			}
		}
		pvnmoThis = pvnmoThis->next;
	    }
        


    return 1;
}

/* rotates about the given axis by theta degrees */
Int2 RotateMol(PMMD pmmdHere,FloatLo theta,Int2 Axis,Int2 Model)
{
        PMGD pmgdHere;
        PMAD pmadHere;
        PALD paldHere;
        PDNMG pdnmgHere;
        PVNMA pvnmaHere;
        vec vHere;

        pdnmgHere=pmmdHere->pdnmgHead;
        if (pdnmgHere==NULL) return 0;
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                if (pmgdHere==NULL) return 0;
                pvnmaHere=pmgdHere->pvnmaAHead;
                if (pvnmaHere==NULL) return 0;
                while (pvnmaHere) {
                        pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
                        if (pmadHere==NULL) return 0;
                        paldHere=GetAtomLocs(pmadHere,Model);
                        if (paldHere!=NULL) {
                                /* first 3 elements of pflvdata are a vec */
                                vHere[0]=(paldHere->pflvData)[0];
                                vHere[1]=(paldHere->pflvData)[1];
                                vHere[2]=(paldHere->pflvData)[2];
                                Rotatecos(Axis,sin(DEGTORAD*theta),cos(DEGTORAD*theta),vHere,vHere);
                                (paldHere->pflvData)[0]=vHere[0];
                                (paldHere->pflvData)[1]=vHere[1];
                                (paldHere->pflvData)[2]=vHere[2];
                        }
                        pvnmaHere=pvnmaHere->next;
                }
                pdnmgHere=pdnmgHere->next;
        }
        /* done */
        return 1;
}

/* returns -1 on error */
static FloatLo GetSS(PMAD pmad1,PMAD pmad2,Int2 Model1,Int2 Model2)
{
	PALD pald1,pald2;
	vec v1,v2;

	if (!pmad1 || !pmad2)
		return -1.0;
	pald1=GetAtomLocs(pmad1,Model1);
	pald2=GetAtomLocs(pmad2,Model2);
	if ((pald1!=NULL) && (pald2!=NULL)) {
		/* first 3 elements of pflvdata are a vec */
		v1[0]=(pald1->pflvData)[0];
		v1[1]=(pald1->pflvData)[1];
		v1[2]=(pald1->pflvData)[2];
		v2[0]=(pald2->pflvData)[0];
		v2[1]=(pald2->pflvData)[1];
		v2[2]=(pald2->pflvData)[2];
		VecSub(v1,v1,v2);
		return Dot(v1,v1);
	}
	return -1.0;
}

/* returns RMSD between two molecules, and assumes they have
   identical chemical graphs; only considers atoms of type
   atomtypes and residues from resa1 to resb1 on pmmd1 and resa2 to resb2 on pmmd2 inclusive */
   
   /* NUS 2008 CWVH Changed to allow RMSD between different chemical graphs when atomtypes is ALIGN_BACKBONE only */
FloatLo GetRMSDEx(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 resa1,Int2 resb1,Int2 resa2,Int2 resb2,Int2 Model1,Int2 Model2)
{
        Int4 numatoms;
        FloatLo SS,RMS,SStmp;
        PMGD pmgd1,pmgd2;
        PMAD pmad1,pmad2;
        PALD pald1,pald2;
        PDNMG pdnmg1,pdnmg2;
        PVNMA pvnma1,pvnma2,pvnmaStep;
        vec v1,v2;
        Char errbuf[PATH_MAX];
        Int2 tmp;

        if (resb1-resa1!=resb2-resa2) {
                sprintf(errbuf,"Residue ranges do not match %d %d",resb1-resa1,resb2-resa2);
                goto rmsd_bail;
        }
        if (resa1>resb1) {
                tmp=resa1;
                resa1=resb1;
                resb1=tmp;
        }
        if (resa2>resb2) {
                tmp=resa2;
                resa2=resb2;
                resb2=tmp;
        }
        SS=0.0;
        numatoms=0;
        pdnmg1=pmmd1->pdnmgHead;
        pdnmg2=pmmd2->pdnmgHead;
        if (pdnmg1==NULL) {
                StringCpy(errbuf,"pdnmg1 is NULL");
                goto rmsd_bail;
        }
        if (pdnmg2==NULL) {
                StringCpy(errbuf,"pdnmg2 is NULL");
                goto rmsd_bail;
        }
        while (pdnmg1 && pdnmg1->choice<resa1)
                pdnmg1=pdnmg1->next;
        while (pdnmg2 && pdnmg2->choice<resa2)
                pdnmg2=pdnmg2->next;
        while (pdnmg1 && pdnmg2) {
                if (pdnmg1->choice>resb1 || pdnmg2->choice>resb2)
                        break;
                pmgd1=(PMGD)(pdnmg1->data.ptrvalue);
                pmgd2=(PMGD)(pdnmg2->data.ptrvalue);
                if (pmgd1==NULL) {
                        StringCpy(errbuf,"pmgd1 is NULL");
                        goto rmsd_bail;
                }
                if (pmgd2==NULL) {
                        StringCpy(errbuf,"pmgd2 is NULL");
                        goto rmsd_bail;
                }
                if ((atomtypes != ALIGN_BACKBONE) && (StringNCmp(pmgd1->pcGraphName,pmgd2->pcGraphName,3))) {
					if (!StringNCmp(pmgd1->pcGraphName,"MET",3) && !StringNCmp(pmgd2->pcGraphName,"MSE",3)) {
					}
					else if (!StringNCmp(pmgd1->pcGraphName,"MSE",3) && !StringNCmp(pmgd2->pcGraphName,"MET",3)) {
					}
					else if (!StringNCmp(pmgd1->pcGraphName,"CEA",3) && !StringNCmp(pmgd2->pcGraphName,"CYS",3)) {
					}
					else if (!StringNCmp(pmgd1->pcGraphName,"CYS",3) && !StringNCmp(pmgd2->pcGraphName,"CEA",3)) {
					}
					else {
                        sprintf(errbuf,"residue names do not match: %s %s",pmgd1->pcGraphName,pmgd2->pcGraphName);
                        goto rmsd_bail;
					}
                }
                pvnma1=pmgd1->pvnmaAHead;
                pvnma2=pmgd2->pvnmaAHead;
                if (pvnma1==NULL) {
                        StringCpy(errbuf,"pvnma1 is NULL");
                        goto rmsd_bail;
                }
                if (pvnma2==NULL) {
                        StringCpy(errbuf,"pvnma2 is NULL");
                        goto rmsd_bail;
                }
		if (atomtypes==ALIGN_CA) {
			pmad1=FindCAlpha(pvnma1);
			pmad2=FindCAlpha(pvnma2);
			SStmp=GetSS(pmad1,pmad2,Model1,Model2);
			if (SStmp>=0.0) {
				SS+=SStmp;
				numatoms++;
			}
		}
		else if (atomtypes == ALIGN_ALLATOM) {
		        while (pvnma1) {
	                        pmad1=(PMAD)(pvnma1->data.ptrvalue);
	                        pmad2=(PMAD)(pvnma2->data.ptrvalue);
	                        if (pmad1==NULL) {
	                                StringCpy(errbuf,"pmad1 is NULL");
	                                goto rmsd_bail;
	                        }
	                        if (pmad2==NULL) {
	                                StringCpy(errbuf,"pmad2 is NULL");
	                                goto rmsd_bail;
	                        }
        	                if ((pmad1->bWhat)!=(pmad2->bWhat)) {
                	                sprintf(errbuf,"bWhats do not match %d %d",pmad1->bWhat,pmad2->bWhat);
                        	        goto rmsd_bail;
	                        }
        	                if (StringCmp(pmad1->pcAName,pmad2->pcAName)) {
                	                sprintf(errbuf,"atom names do not match: %s %s",pmad1->pcAName,pmad2->pcAName);
                        	        goto rmsd_bail;
	                        }
        	                if (((atomtypes&(pmad1->bWhat))!=0) || (atomtypes==ALIGN_ALLATOM)) {
					SStmp=GetSS(pmad1,pmad2,Model1,Model2);
					if (SStmp>=0.0) {
						SS+=SStmp;
						numatoms++;
					}
                        	}
	                        pvnma1=pvnma1->next;
        	                pvnma2=pvnma2->next;
	                }
	                /* sanity check */
	                if (pvnma2!=NULL) {
	                        StringCpy(errbuf,"pvnma2 is not NULL");
	                        goto rmsd_bail;
	                }
		}
		else if (atomtypes == ALIGN_BACKBONE) {
	    	   while (pvnma1) {  /* match up equivalent backbone atoms */
	    	   	      pmad1=(PMAD)(pvnma1->data.ptrvalue);
	                  if (pmad1==NULL) return 0;
	                  pvnmaStep = pvnma2;
	    	   	      while (pvnmaStep) { /* go through all atoms */
	    	   	            pmad2=(PMAD)(pvnmaStep->data.ptrvalue);
  	    	   	            if (pmad2==NULL) return 0;
  	    	   	            if (((IsAtomBackBone((PFB)pmad1)) && (IsAtomBackBone((PFB)pmad2)))   	   	      	
  	 	   	                  && (!StringCmp(pmad1->pcAName,pmad2->pcAName))
  	 	   	                   && ((pmad1->bWhat) == (pmad2->bWhat))) {
  	 	   	                     SStmp=GetSS(pmad1,pmad2,Model1,Model2);
					               if (SStmp>=0.0) {
						              SS+=SStmp;
						              numatoms++;	 	   	                   
  	 	   	                       }
  	 	   	                   }
      	   	      	       pvnmaStep = pvnmaStep->next;
	    	   	      }
	                  pvnma1=pvnma1->next;            
	    	   }
		}
		
		
                pdnmg1=pdnmg1->next;
                pdnmg2=pdnmg2->next;
        }
        RMS=sqrt(SS/((FloatLo)numatoms));
        if (RMS<0.0) {
                ErrPostEx(SEV_ERROR,1,1,"Invalid RMSD computed (SS=%f numat=%d)",SS,numatoms);
                RMS=0.0;
        }
        return RMS;

rmsd_bail:
        ErrPostEx(SEV_ERROR,1,23,"%s",errbuf);
        return 0.0;
}

FloatLo GetRMSD(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model1,Int2 Model2)
{
        return GetRMSDEx(pmmd1,pmmd2,atomtypes,res1,res2,res1,res2,Model1,Model2);
}

/* finds optimal structural alignment between pmmd1 and pmmd1, using only
   atoms specified by atomtypes, and ignoring all but residues between
   res1 and res2 inclusive ( ={1, N} for whole protein) */
Int2 AlignStructures(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model1,Int2 Model2)
{
        Int2 rotx,roty,rotz;
        Int2 bestx,besty,bestz;
        FloatLo rmsdHere,rmsdBest,cur_precis,sweep;

        /* align their centres of mass first */
        if (!CenterMol(pmmd1,atomtypes,res1,res2,Model1)) return 0;
        if (!CenterMol(pmmd2,atomtypes,res1,res2,Model2)) return 0;
        cur_precis=ROTATION_STEPSIZE;
        sweep=360.0;
        do {
                /* now three rotational parameters to determine */
                /* hold pmmd1 fixed and rotate pmmd2 about until an
                   optimal fit is found */
                rmsdBest=GetRMSD(pmmd1,pmmd2,atomtypes,res1,res2,Model1,Model2);
                do {
                        bestx=0;
                        besty=0;
                        bestz=0;
                        if (!RotateMol(pmmd2,-0.5*sweep,Z_AXIS,Model2)) return 0;
                        for (rotz=0;rotz<(Int2)(sweep/cur_precis);rotz++) {
                                if (!RotateMol(pmmd2,cur_precis,Z_AXIS,Model2)) return 0;
                                rmsdHere=GetRMSD(pmmd1,pmmd2,atomtypes,res1,res2,Model1,Model2);
                                if (rmsdHere<rmsdBest) {
                                        rmsdBest=rmsdHere;
                                        bestz=rotz+1-(Int2)((sweep/cur_precis)/2.0);
                                }
                        }
                        if (!RotateMol(pmmd2,(FloatLo)(bestz-rotz)*cur_precis+sweep/2.0,Z_AXIS,Model2)) return 0;
                        if (!RotateMol(pmmd2,-0.5*sweep,Y_AXIS,Model2)) return 0;
                        for (roty=0;roty<(Int2)(sweep/cur_precis);roty++) {
                                if  (!RotateMol(pmmd2,cur_precis,Y_AXIS,Model2)) return 0;
                                rmsdHere=GetRMSD(pmmd1,pmmd2,atomtypes,res1,res2,Model1,Model2);
                                if (rmsdHere<rmsdBest) {
                                        rmsdBest=rmsdHere;
                                        besty=roty+1-(Int2)((sweep/cur_precis)/2.0);
                                }
                        }
                        if (!RotateMol(pmmd2,(FloatLo)(besty-roty)*cur_precis+sweep/2.0,Y_AXIS,Model2)) return 0;
                        if (!RotateMol(pmmd2,-0.5*sweep,X_AXIS,Model2)) return 0;
                        for (rotx=0;rotx<(Int2)(sweep/cur_precis);rotx++) {
                                if (!RotateMol(pmmd2,cur_precis,X_AXIS,Model2)) return 0;
                                rmsdHere=GetRMSD(pmmd1,pmmd2,atomtypes,res1,res2,Model1,Model2);
                                if (rmsdHere<rmsdBest) {
                                        rmsdBest=rmsdHere;
                                        bestx=rotx+1-(Int2)((sweep/cur_precis)/2.0);
                                }
                        }
                        if (!RotateMol(pmmd2,(FloatLo)(bestx-rotx)*cur_precis+sweep/2.0,X_AXIS,Model2)) return 0;
/*printf("Optimal rotation for structure 2: x=%f y=%f z=%f RMSD=%f\n",((FloatLo)bestx)*cur_precis,((FloatLo)besty)*cur_precis,((FloatLo)bestz)*cur_precis,rmsdBest);*/
                } while (bestx || besty || bestz);
                /* temporary storage */
                rmsdHere=sweep;
                sweep=2.0*cur_precis;
                cur_precis=sweep*cur_precis/rmsdHere;
        } while (sweep>ALIGN_PRECISION);
        /* and we are done */
        return 1;
}

Int2 AlignMultipleStructuresEx(PDNMM pdnmmHead1,PDNMM pdnmmHead2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model,Int4 offset1,Int4 offset2, FILE *fp)
{
        PMMD *ppmmdMolList1;
        PMMD *ppmmdMolList2;
        PDNMM pdnmmHere1=NULL,pdnmmHere2=NULL;
        Int2 NumMols1,NumMols2,cnt,leng,m1,m2;

        NumMols1=0;
        NumMols2=0;
        /* get # molecules and store PMMDs in a MolList */
        pdnmmHere1=pdnmmHead1;
        pdnmmHere2=pdnmmHead2;
        leng=((PMMD)(pdnmmHead1->data.ptrvalue))->iResCount;
        if (res1<1 || res1>leng || res2<1 || res2>leng || res1>res2) {
                ErrPostEx(SEV_ERROR,1,11,"Invalid residue range given for aligning");
                return 0;
        }
        while (pdnmmHere1) {
                NumMols1++;
                pdnmmHere1=pdnmmHere1->next;
        }
        while (pdnmmHere2) {
                NumMols2++;
                pdnmmHere2=pdnmmHere2->next;
        }
        if (NumMols1<2 && NumMols2<2) {
                ErrPostEx(SEV_ERROR,1,13,"Need at least 2 molecules for aligning");
                return 0;
        }
        ppmmdMolList1=(PMMD *)MemNew(sizeof(PMMD)*NumMols1);
        ppmmdMolList2=(PMMD *)MemNew(sizeof(PMMD)*NumMols2);
        cnt=0;
        pdnmmHere1=pdnmmHead1;
        pdnmmHere2=pdnmmHead2;
        while (pdnmmHere1) {
                ppmmdMolList1[cnt]=(PMMD)(pdnmmHere1->data.ptrvalue);
                if (ppmmdMolList1[cnt]->iResCount!=leng) {
                        ppmmdMolList1=MemFree(ppmmdMolList1);           
                        ErrPostEx(SEV_ERROR,1,14,"Proteins must be same length for aligning");
                        return 0;
                }
                cnt++;
                pdnmmHere1=pdnmmHere1->next;
        }
        cnt=0;
        pdnmmHere1=pdnmmHead1;
        pdnmmHere2=pdnmmHead2;
        while (pdnmmHere2) {
                ppmmdMolList2[cnt]=(PMMD)(pdnmmHere2->data.ptrvalue);
                if (ppmmdMolList2[cnt]->iResCount!=leng) {
                        ppmmdMolList2=MemFree(ppmmdMolList2);           
                        ErrPostEx(SEV_ERROR,1,14,"Proteins must be same length for aligning");
                        return 0;
                }
                cnt++;
                pdnmmHere2=pdnmmHere2->next;
        }
        for (m1=0;m1<NumMols1;m1++) {
                for (m2=0;m2<NumMols2;m2++) {
                        if (offset1+m1+1>offset2+m2+1) {
                                if (!Align2StrucSVD(ppmmdMolList1[m1],ppmmdMolList2[m2],atomtypes,res1,res2,Model,Model))
                                        return 0;
                                fprintf(fp, "%ld\t\t%ld\t%f\n",(long int)(offset1+m1+1),(long int)(offset2+m2+1),GetRMSD(ppmmdMolList1[m1],ppmmdMolList2[m2],atomtypes,res1,res2,Model,Model));
                        }
                }
        }
        ppmmdMolList1=MemFree(ppmmdMolList1);           
        ppmmdMolList2=MemFree(ppmmdMolList2);           
        return 1;
}

Int2 AlignMultipleStructures(PDNMM pdnmmHead1,PDNMM pdnmmHead2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model,Int4 offset1,Int4 offset2)
{
        return AlignMultipleStructuresEx(pdnmmHead1,pdnmmHead2,atomtypes,res1,res2,Model,offset1,offset2, stdout);
}

Int2 SVD(FloatLo **a,Int2 m,Int2 n,FloatLo *w,FloatLo **v)
{
        Int2 i,its,j,jj,k,l=0,nm=0,flag;
        FloatLo c,f,h,s,x,y,z;
        FloatLo anorm=0.0,g=0.0,scale=0.0;
        FloatLo *rv1;
        
        if (m<n) return 0;
        rv1=(FloatLo *)MemNew(n*sizeof(FloatLo));

  /* Housholder reduction to bidiagonal form. */
  for (i = 0; i < n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m) {
      for (k = i; k < m; k++)
        scale += fabs(a[k][i]);

      if (scale) {
        for (k = i; k < m; k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }

        f = a[i][i];
        g = - SV_SIGN(sqrt(s), f);
        h = f * g - s;
        a[i][i] = f - g;

        if (i != n - 1) {
          for (j = l; j < n; j++) {
                for (s = 0.0, k = i; k < m; k++)
              s += a[k][i] * a[k][j];
            f = s / h;
            for (k = i; k < m; k++)
              a[k][j] += f * a[k][i];
          }
        }

        for (k = i; k < m; k++)
          a[k][i] *= scale;
      }
    }

    w[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m && i != n - 1) {
      for (k = l; k < n; k++)
        scale += fabs(a[i][k]);

      if (scale) {
        for (k = l; k < n; k++) {   
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }

        f = a[i][l];
        g = - SV_SIGN(sqrt(s), f);
        h = f * g - s;
        a[i][l] = f - g;

        for (k = l; k < n; k++)
          rv1[k] = a[i][k] / h;

        if (i != m - 1)
          for (j = l; j < m; j++) {
            for (s = 0.0, k = l; k < n; k++)
              s += a[j][k] * a[i][k];
            for (k = l; k < n; k++)
              a[j][k] += s * rv1[k];
          }

        for (k = l;k < n; k++)
          a[i][k] *= scale;
      }
    }

    anorm = SV_MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

   /* Accumulation of right hand transformations. */
  for (i = n - 1; i >= 0; i--) {
    if (i < n - 1) {
      if (g != 0.0) {
        for (j = l; j < n; j++)
          /* Double division to avoid possible underflow. */
          v[j][i] = (a[i][j] / a[i][l]) / g;

        for (j = l; j < n; j++) {
          for (s = 0.0, k = l; k < n; k++)
            s += a[i][k] * v[k][j];
          for (k = l; k < n; k++)
            v[k][j] += s * v[k][i];
        }
      }

      for (j = l; j < n; j++)
        v[i][j] = v[j][i] = 0.0;
    }
              
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  /* Accumulation of left hand transformations. */
  for (i = n - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];

    if (i < n)
      for (j = l; j < n; j++)
        a[i][j] = 0.0;

    if (g != 0.0) {
      g = 1.0 / g;

      if (i != n - 1) {
        for (j = l; j < n; j++) {
          for (s = 0.0, k = l; k < m; k++)
            s += a[k][i] * a[k][j];
          f = (s / a[i][i]) * g;
          for (k = i; k < m; k++)
            a[k][j] += f * a[k][i];
        }
      }

      for (j = i; j < m; j++)
        a[j][i] *= g;
    } else {
      for (j = i; j < m; j++)
        a[j][i] = 0.0;
    }

    a[i][i] += 1.0;
  }

  /* Diagonalization of the bidiagonal form */
  for (k = n - 1; k >= 0; k--) {           /* Loop over singular values.    */
    for (its = 0; its < MAX_STEP; its++) { /* Loop over allowed iterations. */
      flag = 1;

      for (l = k; l >= 0; l--) {           /* Test for splitting.           */
        nm = l - 1;
        if ((fabs(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }

        if ((fabs(w[nm]) + anorm) == anorm)
          break;
      }

      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          if ((fabs(f) + anorm) == anorm)
            continue;

          g = w[i];
          h = PYTHAG(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g * h;
          s = - f * h;

          for (j = 0; j < m; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z * s;
            a[j][i] = z * c - y * s;
          }
        }
      }

      z = w[k];

      if (l == k) {                        /* Convergence.                  */
        if (z < 0.0) {              /* Singular value is made non negative. */
          w[k] = - z;
          for (j = 0; j < n; j++)
            v[j][k] = - v[j][k];
        }

        break;
      }

      if (its == MAX_STEP) {
        free(rv1);
        return 0;
      }

      /* Shift from bottom 2-by-2 minor. */
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SV_SIGN(g, f))) - h)) / x;

      /* Next QR transformation. */
      c = s = 1.0;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;

        for (jj = 0;jj < n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }

        z = PYTHAG(f, h);
        w[j] = z;
        if (z != 0.0) {             /* Rotation can be arbitrary if z == 0 */
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }

        f = c * g + s * y;
        x = c * y - s * g;

        for (jj = 0; jj < m; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }

      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  rv1=MemFree(rv1);
  return 1;
}

static TrajErr AddAtomToTot(PMAD pmad1,PMAD pmad2,Int2 Model1,Int2 Model2,FloatLo *a[3])
{
	PALD pald1,pald2;
	vec v1,v2;
	Int2 i1,i2;

	if (!pmad1 || !pmad2)
		return ERR_FAIL;
	pald1=GetAtomLocs(pmad1,Model1);
	pald2=GetAtomLocs(pmad2,Model2);
	if ((pald1!=NULL) && (pald2!=NULL)) {
		/* first 3 elements of pflvdata are a vec */
		v1[0]=(pald1->pflvData)[0];
		v1[1]=(pald1->pflvData)[1];
		v1[2]=(pald1->pflvData)[2];
		v2[0]=(pald2->pflvData)[0];
		v2[1]=(pald2->pflvData)[1];
		v2[2]=(pald2->pflvData)[2];
		for (i1=0;i1<3;i1++)
			for (i2=0;i2<3;i2++)
				a[i1][i2]+=v1[i1]*v2[i2];
	}
	return ERR_SUCCESS;
}

/* the resulting structures will be centres, and pmmd1 rotated to align with pmmd2,
   not the other way around (so pmmd2 is unchanged except for the centring */
/* aligns a1 to b1 in pmmd1 to a2 to b2 in pmmd2 */

/*CWVH Dec 2008 - so NOW a call to PutMolBackWhereItCameFrom will put a molecule back to last centered mol position */
/* in this case pmmd1 is last to be centered in Align2StrucSVDEx.  */
/* 2 calls to PutMolBackWhereItCameFrom with pmmd1, pmmd2 will move both back */
/* to original position and orientation of pmmd1 */
/* other models like the vector models are not rotated */

Int2 Align2StrucSVDEx(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 resa1,Int2 resb1,Int2 resa2,Int2 resb2,Int2 Model1,Int2 Model2)
{
        Int2 cnt,i1,i2;
        FloatLo *a[3],w[3],*v[3],RotMat[3][3],am[3*3],vm[3*3],Det;
        PMGD pmgd1,pmgd2;
        PMAD pmad1,pmad2;
        PALD pald1,pald2;
        PDNMG pdnmg1,pdnmg2;
        PVNMA pvnma1,pvnma2, pvnmaStep;
        PMSD pmsdHere;
        PDNMM pdnmmHere;
        PMMD pmmdThis;
        vec v1,v2;
/* CWVH 2012 Rotate the Cyl, Bricks, Solids */	
	PVNMO pvnmoThis = NULL;
	PMOD pmodThis = NULL;
	int i, j, k;

        /* check alignment range is the same length */
        if (resb2-resa2!=resb1-resa1)
                return 0;
        /* ensure res A < res B */
        if (resa1>resb1) {
                cnt=resa1;
                resa1=resb1;
                resb1=cnt;
        }
        if (resa2>resb2) {
                cnt=resa2;
                resa2=resb2;
                resb2=cnt;
        }
        for (cnt=0;cnt<3;cnt++) {
                a[cnt]=am+cnt*3;
                v[cnt]=vm+cnt*3;
        }
        /* align their centres of mass first */
        for (i1=0;i1<3;i1++)
                for (i2=0;i2<3;i2++)
                        a[i1][i2]=0.0;
        if (!CenterMol(pmmd2,atomtypes,resa2,resb2,Model2)) return 0;
        if (!CenterMol(pmmd1,atomtypes,resa1,resb1,Model1)) return 0;
        pdnmg1=pmmd1->pdnmgHead;
        pdnmg2=pmmd2->pdnmgHead;
        if (pdnmg1==NULL) return 0;
        if (pdnmg2==NULL) return 0;
        while (pdnmg1 && (pdnmg1->choice<resa1)) {
                pdnmg1=pdnmg1->next;
        }
        while (pdnmg2 && (pdnmg2->choice<resa2)) {
                pdnmg2=pdnmg2->next;
        }
        while (pdnmg1 && pdnmg2) {
                if (pdnmg2->choice>resb2 || pdnmg1->choice>resb1)
                        break;
                pmgd1=(PMGD)(pdnmg1->data.ptrvalue);
                pmgd2=(PMGD)(pdnmg2->data.ptrvalue);
                if (pmgd1==NULL) return 0;
                if (pmgd2==NULL) return 0;
		if (atomtypes==ALIGN_CA) {
			if (!StringNCmp(pmgd1->pcGraphName,"MET",3) && !StringNCmp(pmgd2->pcGraphName,"MSE",3)) {
			}
			else if (!StringNCmp(pmgd1->pcGraphName,"MSE",3) && !StringNCmp(pmgd2->pcGraphName,"MET",3)) {
			}
			else if (!StringNCmp(pmgd1->pcGraphName,"CEA",3) && !StringNCmp(pmgd2->pcGraphName,"CYS",3)) {
			}
			else if (!StringNCmp(pmgd1->pcGraphName,"CYS",3) && !StringNCmp(pmgd2->pcGraphName,"CEA",3)) {
			}
	                else if (StringNCmp(pmgd1->pcGraphName,pmgd2->pcGraphName,3)) return 0;
		}
		else {
	        if (atomtypes==ALIGN_BACKBONE) { ; }
	        else if (StringCmp(pmgd1->pcGraphName,pmgd2->pcGraphName)) return 0;
	        /* this allows us to do backbone based alignments for non-matching chemical graphs */
		}
        pvnma1=pmgd1->pvnmaAHead;
        pvnma2=pmgd2->pvnmaAHead;
        if (pvnma1==NULL) return 0;
        if (pvnma2==NULL) return 0;
		if (atomtypes==ALIGN_CA) {
			pmad1=FindCAlpha(pvnma1);
			pmad2=FindCAlpha(pvnma2);
			AddAtomToTot(pmad2,pmad1,Model2,Model1,a);
		}
		else if (atomtypes == ALIGN_ALLATOM) {
	                while (pvnma1) {
	                        pmad1=(PMAD)(pvnma1->data.ptrvalue);
	                        pmad2=(PMAD)(pvnma2->data.ptrvalue);
	                        if (pmad1==NULL) return 0;
	                        if (pmad2==NULL) return 0;
        	                if (StringCmp(pmad1->pcAName,pmad2->pcAName)) return 0;
                	        if ((pmad1->bWhat)!=(pmad2->bWhat)) return 0;
                	        AddAtomToTot(pmad2, pmad1, Model2, Model1, a); 
        	                pvnma1=pvnma1->next;
	                        pvnma2=pvnma2->next;
	                }   
	          }      /* sanity check * if (pvnma2!=NULL) return 0; */
	          
	    else if (atomtypes == ALIGN_BACKBONE) {
	    	   while (pvnma1) {  /* match up equivalent backbone atoms */
	    	   	      pmad1=(PMAD)(pvnma1->data.ptrvalue);
	                  if (pmad1==NULL) return 0;
	                  pvnmaStep = pvnma2;
	    	   	      while (pvnmaStep) { /* go through all atoms */
	    	   	            pmad2=(PMAD)(pvnmaStep->data.ptrvalue);
  	    	   	            if (pmad2==NULL) return 0;
  	    	   	            if (((IsAtomBackBone((PFB)pmad1)) && (IsAtomBackBone((PFB)pmad2)))   	   	      	
  	 	   	                  && (!StringCmp(pmad1->pcAName,pmad2->pcAName))
  	 	   	                   && ((pmad1->bWhat) == (pmad2->bWhat))) 
  	 	   	                   AddAtomToTot(pmad2, pmad1, Model2, Model1, a);
      	    	   	      	pvnmaStep = pvnmaStep->next;
	    	   	      }
	                  pvnma1=pvnma1->next;            
	    	   }
	    	   
	    	 }
            pdnmg1=pdnmg1->next;
            pdnmg2=pdnmg2->next;
        }
        SVD(a,3,3,w,v);
        for (i1=0;i1<3;i1++)
                for (i2=0;i2<3;i2++)
                        RotMat[i1][i2]=a[i2][0]*v[i1][0]+a[i2][1]*v[i1][1]+a[i2][2]*v[i1][2];
        Det=RotMat[0][0]*RotMat[1][1]*RotMat[2][2]+RotMat[0][1]*RotMat[1][2]*RotMat[2][0]+RotMat[0][2]*RotMat[1][0]*RotMat[2][1]-RotMat[2][0]*RotMat[1][1]*RotMat[0][2]-RotMat[2][1]*RotMat[1][2]*RotMat[0][0]-RotMat[2][2]*RotMat[1][0]*RotMat[0][1];
        if (Det<0.0) {
                if (w[0]<w[1] && w[0]<w[2])
                        for (i1=0;i1<3;i1++)
                                for (i2=0;i2<3;i2++)
                                        RotMat[i1][i2]=-a[i2][0]*v[i1][0]+a[i2][1]*v[i1][1]+a[i2][2]*v[i1][2];
                else if (w[1]<w[2])
                        for (i1=0;i1<3;i1++)
                                for (i2=0;i2<3;i2++)
                                        RotMat[i1][i2]=a[i2][0]*v[i1][0]-a[i2][1]*v[i1][1]+a[i2][2]*v[i1][2];
                else
                        for (i1=0;i1<3;i1++)
                                for (i2=0;i2<3;i2++)
                                        RotMat[i1][i2]=a[i2][0]*v[i1][0]+a[i2][1]*v[i1][1]-a[i2][2]*v[i1][2];
        }
        Normalize(RotMat[0],RotMat[0]);
        for (i1=1;i1<3;i1++) {
                v1[0]=0.0;
                v1[1]=0.0;
                v1[2]=0.0;
                for (i2=0;i2<i1;i2++) {
                        VecScale(v2,RotMat[i2],Dot(RotMat[i1],RotMat[i2]));
                        VecAdd(v1,v1,v2);
                }
                VecSub(RotMat[i1],RotMat[i1],v1);
                Normalize(RotMat[i1],RotMat[i1]);                       
        }       
        /* CWVH 2008 Dec - changed to rotate all molecules in the structure - but not all the Models..*/
        /* CWVH 2009 Feb - fixed rotate bug - we have to rotate pmmd2 into frame of pmmd1, not vice versa */
        pmsdHere = ToMSDParent( (PFB) pmmd2);
        pdnmmHere = pmsdHere->pdnmmHead;
        while (pdnmmHere) { 
        pmmdThis=(PMMD)(pdnmmHere->data.ptrvalue);
        pdnmg2=pmmdThis->pdnmgHead;
        while (pdnmg2) {
                pmgd2=(PMGD)(pdnmg2->data.ptrvalue);
                pvnma2=pmgd2->pvnmaAHead;
                while (pvnma2) {
                        pmad2=(PMAD)(pvnma2->data.ptrvalue);
                        pald2=GetAtomLocs(pmad2,Model2);
                        if (pald2!=NULL) {
                                /* first 3 elements of pflvdata are a vec */
                                v1[0]=(pald2->pflvData)[0];
                                v1[1]=(pald2->pflvData)[1];
                                v1[2]=(pald2->pflvData)[2];
                                for (i1=0;i1<3;i1++)
                                        v2[i1]=Dot(RotMat[i1],v1);
                                (pald2->pflvData)[0]=v2[0];
                                (pald2->pflvData)[1]=v2[1];
                                (pald2->pflvData)[2]=v2[2];
                        }
                        pvnma2=pvnma2->next;
                }
                pdnmg2=pdnmg2->next;
        }
        pdnmmHere = pdnmmHere->next;
        }

 /* CWVH 2012 Rotate the Solids as well  */

 	pmsdHere = ToMSDParent( (PFB) pmmd2);
        pvnmoThis = pmsdHere->pvnmoHead;
	  while (pvnmoThis != NULL)
	    {
 		pmodThis = (PMOD) pvnmoThis->data.ptrvalue;
		if (pmodThis != NULL) {
 			if (pmodThis->ppflObject != NULL) {
				for (i = 0; i < pmodThis->iCoordNo; i++) {	 
					for (j = 0; j < 3; j++) {
						v1[j] = (pmodThis->ppflObject)[i][j]; 
					} 
					for (i1=0;i1<3;i1++)   
                                        	v2[i1]=Dot(RotMat[i1],v1);
					for (k = 0; k < 3; k++) {  
						(pmodThis->ppflObject)[i][k] = v2[k]; 
					} 
				}
			}
		}
		pvnmoThis = pvnmoThis->next;
	    }
        
  
        return 1;
}

Int2 Align2StrucSVD(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model1,Int2 Model2)
{
        return Align2StrucSVDEx(pmmd1,pmmd2,atomtypes,res1,res2,res1,res2,Model1,Model2);
}

FloatLo GetRgyr(PMMD pmmdRoot,Int2 Model)
{
        PDNMG pdnmgHere;
        PMGD pmgdHere;
        vec vCofM,vCAHere,vTmp;
        vec vZero={0.0,0.0,0.0};
        Int2 clength;
        FloatLo Rgyr;

        if (pmmdRoot==NULL)
                return 0.0;
        clength=0;
        pdnmgHere=pmmdRoot->pdnmgHead;
        VecScale(vCofM,vZero,1.0);
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                if (GetCoOrds(pmgdHere," CA ",vZero,vCAHere,Model)) {
					VecAdd(vCofM,vCofM,vCAHere);
					clength++;
				}
                pdnmgHere=pdnmgHere->next;
        }
		if (clength)
	        VecScale(vCofM,vCofM,1.0/(FloatLo)clength);
        Rgyr=0.0;
        pdnmgHere=pmmdRoot->pdnmgHead;
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                if (GetCoOrds(pmgdHere," CA ",vZero,vCAHere,Model)) {
					VecSub(vTmp,vCAHere,vCofM);
					Rgyr+=Dot(vTmp,vTmp);
				}
                pdnmgHere=pdnmgHere->next;
        }
        Rgyr/=(FloatLo)clength;
        Rgyr=sqrt(Rgyr);
        return Rgyr;
}

FloatLo GetHPRgyr(PMMD pmmdRoot,Int2 Model)
{
        PDNMG pdnmgHere;
        PMGD pmgdHere;
		PALD pald;
        vec vCofM,vCBHere,vTmp;
        vec vZero={0.0,0.0,0.0};
        Int2 clength,numadd;
        FloatLo Rgyr;
        Char resHere;
        Char hpresidues[]="MILVYCFGAW";

        if (pmmdRoot==NULL)
                return 0.0;
        numadd=0;
        clength=0;
        pdnmgHere=pmmdRoot->pdnmgHead;
        VecScale(vCofM,vZero,1.0);
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                resHere=GetAAFromIDict(pmgdHere);
                if (resHere!='G') {
					if (GetCoOrds(pmgdHere," CB ",vZero,vCBHere,Model)) {
		                VecAdd(vCofM,vCofM,vCBHere);
						clength++;
					}
				}
                else {
					if (GetCoOrds(pmgdHere,"2HA ",vZero,vCBHere,Model)) {
		                VecAdd(vCofM,vCofM,vCBHere);
						clength++;
					}
				}
                pdnmgHere=pdnmgHere->next;
        }
        VecScale(vCofM,vCofM,1.0/(FloatLo)clength);
        Rgyr=0.0;
        pdnmgHere=pmmdRoot->pdnmgHead;
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                resHere=GetAAFromIDict(pmgdHere);
                if (resHere!='G')
                        pald=GetCoOrds(pmgdHere," CB ",vZero,vCBHere,Model);
                else
                        pald=GetCoOrds(pmgdHere,"2HA ",vZero,vCBHere,Model);
                if (pald && StringRChr(hpresidues,(int)resHere)!=NULL) {
                        VecSub(vTmp,vCBHere,vCofM);
                        Rgyr+=Dot(vTmp,vTmp);
                        numadd++;
                }
                pdnmgHere=pdnmgHere->next;
        }
        if (numadd==0)
                return 0.0;
        Rgyr/=(FloatLo)numadd;
        Rgyr=sqrt(Rgyr);
        return Rgyr;
}

FloatLo GetExactRgyr(PMMD pmmdRoot,Int2 Model,Boolean HPonly)
{
        PDNMG pdnmgHere;
        PMGD pmgdHere;
        PVNMA pvnmaHere;
        PMAD pmadHere;
        PALD paldHere;
        vec vCofM,vHere,vTmp;
        vec vZero={0,0,0};
        FloatLo Rgyr,totmass,mass;
        Uint1 element;
        Char hpresidues[]="MILVYCFGAW";
        Char resHere;
        
        if (pmmdRoot==NULL)
                return 0.0;
        pdnmgHere=pmmdRoot->pdnmgHead;
        totmass=0.0;
        VecScale(vCofM,vZero,1.0);
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                pvnmaHere=pmgdHere->pvnmaAHead;
                while (pvnmaHere) {
                        pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
                        element=pvnmaHere->choice;
                        switch (element) {
                                case 1:
                                        mass=AW_H;
                                        break;
                                case 6:
                                        mass=AW_C;
                                        break;
                                case 7:
                                        mass=AW_N;
                                        break;
                                case 8:
                                        mass=AW_O;
                                        break;
                                case 15:
                                        mass=AW_P;
                                        break;
                                case 16:
                                        mass=AW_S;
                                        break;
                                case 34:
                                        mass=AW_SE;
                                        break;
                                default:
                                        mass=AW_C;
                        }
                        paldHere=GetAtomLocs(pmadHere,Model);
                        if (paldHere!=NULL) {
                                vHere[0]=(paldHere->pflvData)[0];
                                vHere[1]=(paldHere->pflvData)[1];
                                vHere[2]=(paldHere->pflvData)[2];
                                VecScale(vHere,vHere,mass);
                                VecAdd(vCofM,vCofM,vHere);
                                totmass+=mass;
                        }
                        pvnmaHere=pvnmaHere->next;
                }
                pdnmgHere=pdnmgHere->next;
        }
        VecScale(vCofM,vCofM,1.0/totmass);
        Rgyr=0.0;
        totmass=0.0;
        pdnmgHere=pmmdRoot->pdnmgHead;
        while (pdnmgHere) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
                resHere=GetAAFromIDict(pmgdHere);
                if ((!HPonly) || (StringRChr(hpresidues,(int)resHere)!=NULL)) {
                        pvnmaHere=pmgdHere->pvnmaAHead;
                        while (pvnmaHere) {
                                pmadHere=(PMAD)(pvnmaHere->data.ptrvalue);
                                element=pvnmaHere->choice;
                                switch (element) {
                                        case 1:
                                                mass=AW_H;
                                                break;
                                        case 6:
                                                mass=AW_C;
                                                break;
                                        case 7:
                                                mass=AW_N;
                                                break;
                                        case 8:
                                                mass=AW_O;
                                                break;
                                        case 15:
                                                mass=AW_P;
                                                break;
                                        case 16:
                                                mass=AW_S;
                                                break;
                                        case 34:
                                                mass=AW_SE;
                                                break;
                                        default:
                                                mass=AW_C;
                                }
                                paldHere=GetAtomLocs(pmadHere,Model);
                                if (paldHere!=NULL) {
                                        vHere[0]=(paldHere->pflvData)[0];
                                        vHere[1]=(paldHere->pflvData)[1];
                                        vHere[2]=(paldHere->pflvData)[2];
                                        VecSub(vTmp,vHere,vCofM);
                                        Rgyr=Rgyr+mass*Dot(vTmp,vTmp);
                                        totmass+=mass;
                                }
                                pvnmaHere=pvnmaHere->next;
                        }
                }
                pdnmgHere=pdnmgHere->next;
        }
        Rgyr/=totmass;
        Rgyr=sqrt(Rgyr);
        return Rgyr;
}

Int4 CalcExtendedResidues(PMMD pmmdHere,Int2 Model)
{
        PDNMG pdnmgHere;
        PMGD pmgdHere;
        Int2 cnt,retval=0;
        vec *vecarray,vZero={0,0,0},vTmp;
		Boolean *missing;
        Boolean *resstraight;

        cnt=0;
        vecarray=(vec *)MemNew((pmmdHere->iResCount)*sizeof(vec));
		missing=(Boolean *)MemNew((pmmdHere->iResCount)*sizeof(Boolean));
        /* zero this to all FALSE */
        resstraight=(Boolean *)MemNew((pmmdHere->iResCount)*sizeof(Boolean));
        pdnmgHere=pmmdHere->pdnmgHead;
        /* read Ca co-ordinates into an array */
        while (pdnmgHere!=NULL) {
                pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);
				if (GetCoOrds(pmgdHere," CA ",vZero,vecarray[cnt],Model)==NULL) {
					missing[cnt]=TRUE;
				}
                cnt++;
                pdnmgHere=pdnmgHere->next;
        }
        /* and analyze */
        cnt=0;
        while ((cnt+4)<(pmmdHere->iResCount)) {
			if (!missing[cnt] && !missing[cnt+4]) {
                VecSub(vTmp,vecarray[cnt+4],vecarray[cnt]);
                if (getMag(vTmp)>13.25) {  /* 13.25 ANGSTROMS */
                        resstraight[cnt]=TRUE;
                        resstraight[cnt+1]=TRUE;
                        resstraight[cnt+2]=TRUE;
                        resstraight[cnt+3]=TRUE;
                        resstraight[cnt+4]=TRUE;
                }
			}
			cnt++;
        }
        for (cnt=0;cnt<pmmdHere->iResCount;cnt++)
                if (resstraight[cnt]==TRUE)
                        retval++;
        resstraight=MemFree(resstraight);
        vecarray=MemFree(vecarray);
		missing=MemFree(missing);
        return retval;
}

static TrajErr GetBackboneAtomCoords(PDNMG pdnmg, vec v0, vec v1, vec v2,Int2 Model)
{
        PMGD pmgdHere;
        vec vZero={0.0,0.0,0.0};
        
        pmgdHere=(PMGD)(pdnmg->data.ptrvalue);
        if (GetCoOrds(pmgdHere," N  ",vZero,v0,Model)==NULL)
                return ERR_FAIL;
        if (GetCoOrds(pmgdHere," CA ",vZero,v1,Model)==NULL)
                return ERR_FAIL;
        if (GetCoOrds(pmgdHere," C  ",vZero,v2,Model)==NULL)
                return ERR_FAIL;
        return ERR_SUCCESS;
}

static TrajErr GetThreeCACoords(PDNMG pdnmg0, PDNMG pdnmg1, PDNMG pdnmg2, vec v0, vec v1, vec v2,Int2 Model)
{
        PMGD pmgd0,pmgd1,pmgd2;
        PVNMA pvnma0,pvnma1,pvnma2;
        PMAD pmad0,pmad1,pmad2;
        PALD pald0,pald1,pald2;
        
        pmgd0=(PMGD)(pdnmg0->data.ptrvalue);
        pmgd1=(PMGD)(pdnmg1->data.ptrvalue);
        pmgd2=(PMGD)(pdnmg2->data.ptrvalue);
        pvnma0=pmgd0->pvnmaAHead;
        pvnma1=pmgd1->pvnmaAHead;
        pvnma2=pmgd2->pvnmaAHead;
        if ((pmad0=FindCAlpha(pvnma0))==NULL)
                return ERR_FAIL;
        if ((pmad1=FindCAlpha(pvnma1))==NULL)
                return ERR_FAIL;
        if ((pmad2=FindCAlpha(pvnma2))==NULL)
                return ERR_FAIL;
        pald0=GetAtomLocs(pmad0,Model);
        if (pald0==NULL)
                return ERR_FAIL;
        pald1=GetAtomLocs(pmad1,Model);
        if (pald1==NULL)
                return ERR_FAIL;
        pald2=GetAtomLocs(pmad2,Model);
        if (pald2==NULL)
                return ERR_FAIL;
        v0[0]=AtomLocX(pald0);
        v0[1]=AtomLocY(pald0);
        v0[2]=AtomLocZ(pald0);
        v1[0]=AtomLocX(pald1);
        v1[1]=AtomLocY(pald1);
        v1[2]=AtomLocZ(pald1);
        v2[0]=AtomLocX(pald2);
        v2[1]=AtomLocY(pald2);
        v2[2]=AtomLocZ(pald2);  
        return ERR_SUCCESS;
}

PNN ComputeTakeoffAngles(PMMD pmmd,Int2 Model,Int4 res1,Int4 res2,Int2 walk)
{
        PDNMG pdnmgHere,pdnmgLast,pdnmgLastLast;
/*      PMGD pmgdHere;
        PMAD pmad1,pmad2;*/
        vec v0,v1,v2,v3,v4,v5;
        FloatLo dh12,dh23,dh34,ba1,ba2,ba3,ba4,bl01,bl12,bl23,bl34,bl45;
        PNN pnnHere;
        
        pdnmgHere=pmmd->pdnmgHead;
        while (pdnmgHere->choice<res1) {
                pdnmgHere=pdnmgHere->next;
                if (pdnmgHere==NULL)
                        return NULL;
        }
        if (pdnmgHere->choice!=res1)
                return NULL;
/*      pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);*/
        if (walk==WALK_CA) {
                pdnmgLast=pdnmgHere->last;
                if (pdnmgLast==NULL)
                        return NULL;
                pdnmgLastLast=pdnmgLast->last;
                if (pdnmgLastLast==NULL)
                        return NULL;
                /* use model 1, assumed to be most accurate */
                if (GetThreeCACoords(pdnmgLastLast,pdnmgLast,pdnmgHere,v0,v1,v2,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        else if (walk==WALK_PHIPSI) {
                if (GetBackboneAtomCoords(pdnmgHere,v0,v1,v2,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        while (pdnmgHere->choice<res2) {
                pdnmgHere=pdnmgHere->next;
                if (pdnmgHere==NULL)
                        return NULL;
        }
        if (pdnmgHere->choice!=res2)
                return NULL;
/*      pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);*/
        if (walk==WALK_CA) {    
                pdnmgLast=pdnmgHere->next;
                if (pdnmgLast==NULL)
                        return NULL;    
                pdnmgLastLast=pdnmgLast->next;
                if (pdnmgLastLast==NULL)
                        return NULL;
                if (GetThreeCACoords(pdnmgHere,pdnmgLast,pdnmgLastLast,v3,v4,v5,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        else if (walk==WALK_PHIPSI) {
                if (GetBackboneAtomCoords(pdnmgHere,v3,v4,v5,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        /* compute restraint parameters now */
        /* should be same for PHIPSI and CA walk... (?) */
        GetDihedral(v0,v1,v2,v3,-180.0,&dh12,&ba1,&ba2,&bl01,&bl12,&bl23);
        GetDihedral(v1,v2,v3,v4,-180.0,&dh23,&ba2,&ba3,&bl12,&bl23,&bl34);
        GetDihedral(v2,v3,v4,v5,-180.0,&dh34,&ba3,&ba4,&bl23,&bl34,&bl45);
        /* required pieces of information are dh12,dh23,dh34,ba2,ba3,ba4,bl23
                 but ba4 is phi of Trajectory Distribution at atom 5 taken from
                 crystal structure, so already correct */
        /* thus 1 distance, 2 angles, 3 dihedrals needed to uniquely fix, 6 DOF */      
        /* now make the actual restraint */
        pnnHere=(PNN)MemNew(sizeof(NN));
        if (walk==WALK_CA) {
                StringCpy(pnnHere->AtomName1," CA ");
                StringCpy(pnnHere->AtomName2," CA ");
        }
        else if (walk==WALK_PHIPSI) {
                StringCpy(pnnHere->AtomName1," C  ");
                StringCpy(pnnHere->AtomName2," N  ");
        }
        pnnHere->MeanDist=bl23;
        pnnHere->Angle0=ba1;
        pnnHere->Angle1=ba2;
        pnnHere->Angle2=ba3;
        pnnHere->Angle3=ba4;
        pnnHere->Dihedral01=dh12;
        pnnHere->Dihedral12=dh23;
        pnnHere->Dihedral23=dh34;
        return pnnHere;
}

vec *ComputeTakeoffCoords(PMMD pmmd,Int2 Model,Int4 res1,Int4 res2,Int2 walk)
{
        PDNMG pdnmgHere,pdnmgLast,pdnmgLastLast;
/*      PMGD pmgdHere;
        PMAD pmad1,pmad2;*/
        vec v0,v1,v2,v3,v4,v5,*v;
        
        pdnmgHere=pmmd->pdnmgHead;
        while (pdnmgHere->choice<res1) {
                pdnmgHere=pdnmgHere->next;
                if (pdnmgHere==NULL)
                        return NULL;
        }
        if (pdnmgHere->choice!=res1)
                return NULL;
/*      pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);*/
        if (walk==WALK_CA) {
                pdnmgLast=pdnmgHere->last;
                if (pdnmgLast==NULL)
                        return NULL;
                pdnmgLastLast=pdnmgLast->last;
                if (pdnmgLastLast==NULL)
                        return NULL;
                /* use model 1, assumed to be most accurate */
                if (GetThreeCACoords(pdnmgLastLast,pdnmgLast,pdnmgHere,v0,v1,v2,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        else if (walk==WALK_PHIPSI) {
                if (GetBackboneAtomCoords(pdnmgHere,v0,v1,v2,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        while (pdnmgHere->choice<res2) {
                pdnmgHere=pdnmgHere->next;
                if (pdnmgHere==NULL)
                        return NULL;
        }
        if (pdnmgHere->choice!=res2)
                return NULL;
/*      pmgdHere=(PMGD)(pdnmgHere->data.ptrvalue);*/
        if (walk==WALK_CA) {    
                pdnmgLast=pdnmgHere->next;
                if (pdnmgLast==NULL)
                        return NULL;    
                pdnmgLastLast=pdnmgLast->next;
                if (pdnmgLastLast==NULL)
                        return NULL;
                if (GetThreeCACoords(pdnmgHere,pdnmgLast,pdnmgLastLast,v3,v4,v5,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        else if (walk==WALK_PHIPSI) {
                if (GetBackboneAtomCoords(pdnmgHere,v3,v4,v5,Model)!=ERR_SUCCESS)
                        return NULL;
        }
        v=(vec *)MemNew(sizeof(vec)*6);
        VecScale(v[0],v0,1.0);
        VecScale(v[1],v1,1.0);
        VecScale(v[2],v2,1.0);
        VecScale(v[3],v3,1.0);
        VecScale(v[4],v4,1.0);
        VecScale(v[5],v5,1.0);
        return v;
}

/*
$Log: rotate.c,v $
Revision 1.41  2009/02/16 17:49:52  chogue
fixed code to rotate second molecule passed as argument

Revision 1.40  2009/01/20 16:12:35  chogue

Rotate now moves all chains in a structure - not just the first one.

Revision 1.39  2008/12/11 17:18:54  chogue
modified rotate.c to store previous location of centered molecules and put them back where they came from with a new function.  All molecules in a structure are rotated and translated to allow alignment function to dock moleucles via peptides.

Revision 1.38  2008/12/11 15:22:19  chogue
allowed rotate library to superimpose non-identical chem graphs with backbone coordinates ALIGN_BACKBONE.

Revision 1.37  2004/07/19 21:56:48  hfeldman
Fixed typo

Revision 1.36  2004/07/16 15:08:02  hfeldman
Allow fuzzy name matching in SVD align for CA only

Revision 1.35  2004/06/29 18:07:07  hfeldman
Added special case exceptions for RMSD residue names being the same eg selenocystine=cystine

Revision 1.34  2004/06/23 20:27:26  hfeldman
resolved name conflict with windows

Revision 1.33  2004/06/23 20:01:11  fwu
only check if 1st 3 letters of residue name match for RMSD computation

Revision 1.32  2004/06/23 17:37:38  hfeldman
Correct deal with residues missing CA co-ordinates (by not including them in calculations)

Revision 1.31  2004/06/21 16:32:04  hfeldman
for CA RMSD, now look ONLY at CA atoms - other atoms need not match..

Revision 1.30  2004/06/17 22:08:42  hfeldman
Fixed invalid chars in file

Revision 1.29  2004/06/16 22:08:06  hfeldman
Added support for RMSD between things of different length e.g. domain vs. full protein

Revision 1.28  2003/08/03 20:49:44  feldman
Added extra error messages for problems with GetRMSD

Revision 1.27  2003/01/27 15:49:04  feldman
Fixed possible naming conflict on some OSes

Revision 1.26  2003/01/24 16:52:53  feldman
Added setbit to match clearbit

Revision 1.25  2002/10/25 14:44:07  feldman
Added bad RMS warning

Revision 1.24  2002/08/22 21:06:51  feldman
Added better bumpchecking for loop fragments

Revision 1.23  2002/08/20 19:28:25  michel
multalign results output to command line parameter file while keeping track of progresss

Revision 1.22  2002/08/19 21:15:16  michel
Multalign reads wildcard files from directory. Clusterrmsd requires multalign result filename as input

Revision 1.21  2002/07/31 22:22:59  feldman
Added two outside angles to ComputeTakeoffAngles
Fixed some bugs in fragment minimization for HM

Revision 1.20  2001/09/24 14:50:40  feldman
minor bugfix

Revision 1.19  2001/09/04 18:03:57  feldman
Fixed compiler warning

Revision 1.18  2001/08/31 17:45:26  feldman
Moved takeoff angle calculation code from clust.c to rotate.c

Revision 1.17  2001/07/09 14:26:30  feldman
Changed Atomic mass calculation by Element rather than atom name

Revision 1.16  2001/06/19 18:42:54  phan
changed rgyr functions to take pmmd

Revision 1.15  2001/05/04 16:18:59  feldman
Added Gly, Ala and Trp as hydrophobic residues when getting rgyr and
surface area of hydrophobics, based on Eisenberg ranking

Revision 1.14  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.13  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.12  2001/03/09 17:33:59  feldman
Added code to deal with rotamers.  Now rotamer chi angles can
be encoded into a 4-byte integer and used by foldtraj when
building the structures.  These can be edited by VisTraj as well

Revision 1.11  2001/02/07 22:31:42  feldman
fixed remaining code for unused variables and unreachable code, etc.

Revision 1.10  2001/01/18 18:10:38  feldman
Added bzip2 1.01 compatibility

Revision 1.9  2001/01/16 22:01:36  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.8  2000/11/17 22:28:59  feldman
Minor bugfixes

Revision 1.7  2000/10/24 20:57:25  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.6  2000/08/14 20:24:19  feldman
Added LIBCALLBACK for function pointers where needed

Revision 1.5  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.4  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.3  2000/06/22 13:05:02  feldman
Fixed minor bugs with log file nameing,
fixed error in extended residue calculation and DSSP usage
added alterresidue function

Revision 1.2  2000/06/20 16:40:23  feldman
Incorporated sstru extended structure calculation into C code,
no longer need an external program

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

