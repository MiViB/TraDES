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


#include <mmdbtraj.h>
#include <geometry.h>

/* global variables */
static PVNMA pvnmaHead;
static vec vRef[9];
extern CharPtr NCBIstdaaUC;
/* I made these up for Gly */
/* note that the last two elements correspond to CYX and cis-PRO respectively */
/* these values are from Rey and Skolnick as well, not Engh and Huber */
FloatLo bl_cacb[]={1.530F,0.0F,1.528F,1.533F,1.531F,1.534F,1.070F,1.542F,1.554F,1.528F,1.536F,1.528F,1.534F,1.527F,1.529F,1.532F,1.530F,1.560F,1.540F,1.534F,0.0F,1.541F,0.0F,0.0F,0.0F,0.0F,1.536F};
FloatLo ba_ncacb[]={109.9F,0.0F,110.3F,110.7F,110.9F,111.1F,109.9F,110.9F,111.1F,109.9F,109.4F,110.9F,110.1F,103.2F,110.7F,110.9F,110.2F,110.9F,110.8F,110.8F,0.0F,110.3F,0.0F,0.0F,0.0F,0.0F,103.2F};
FloatLo ba_ccacb[]={110.2F,0.0F,110.5F,111.1F,109.2F,110.8F,110.2F,110.1F,111.6F,109.5F,111.2F,110.6F,111.4F,111.4F,110.4F,109.9F,110.0F,110.9F,111.9F,110.5F,0.0F,110.3F,0.0F,0.0F,0.0F,0.0F,110.8F};

/* allocates a new PALD in a given model for a given PMAD */
PALD AllocNewLoc(PMAD pmadThis,Int2 modelnum)
{
	PVNAL pvnalThis=NULL,pvnalTemp=NULL;
	PALD paldThis=NULL,paldTemp=NULL;
/*	BiostrucModelPtr pbsmThis=NULL;
	PMSD pmsdThis=NULL;
	PMLD pmldThis=NULL;
	ModelCoordinateSetPtr pmcsThis=NULL;
	BiostrucPtr pbsBS=NULL;
	PDNML pdnmlModel=NULL; 
 
	pmsdThis=ToMSDParent((PFB)pmadThis);
	pbsBS=pmsdThis->pbsBS;
	pbsmThis=pbsBS->model;
	pdnmlModel=pmsdThis->pdnmlModels;
	while(pdnmlModel->choice!=modelnum) {
		pdnmlModel=pdnmlModel->next;
	}
	pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
   	pmcsThis=pbsmThis->model_coordinates; 	
	pmldThis->iNoCoordSet=pmcsThis->id;  */
	/* modelnum (2nd arg) is actually rotamer number */
	pvnalThis=NewVNAL(NULL,modelnum);
        if (!pvnalThis) goto allocerr;
        paldThis=(PALD)pvnalThis->data.ptrvalue;
	/* back-link to pmad */
        paldThis->pfbParent=(PFB)pmadThis;
	/* alloc the data vector */  
        paldThis->pflvData=FLVector(0,(Int4)4);
        if (!paldThis->pflvData) goto allocerr;
	paldThis->iCoordSet=1; /*(Int1)((pbsmThis->model_coordinates)->id);*/
        paldThis->iFloatNo=(Int1)4;
        if (pmadThis->pvnalLocate) {
		pvnalTemp=(PVNAL)ValNodeFindNext(pmadThis->pvnalLocate,NULL,modelnum);
		/* the case where THIS model has a location already */
                if (pvnalTemp) { 
                	paldTemp=(PALD)pvnalTemp->data.ptrvalue;
                        while (paldTemp->next) 
				/* find end of pald list */
                                paldTemp = paldTemp->next;
                        paldTemp->next=(PALD)pvnalThis->data.ptrvalue;
			/* attach the new loc */
                        pvnalThis->data.ptrvalue=NULL;
			/* lose the ValNode baggage */
                        pvnalThis=ValNodeFree(pvnalThis); 
        	}
                else
			/* attach to pmadThis */
                	ValNodeListCat(pmadThis->pvnalLocate, pvnalThis);
        }
        else
        	pmadThis->pvnalLocate=pvnalThis;                
	/* any co-ordinates for now */
        paldThis->pflvData[0]=1;
        paldThis->pflvData[1]=2;
        paldThis->pflvData[2]=3;
	return paldThis;
allocerr:
	ErrPostEx(SEV_FATAL,0,0,"Out of memory While building.");
	return NULL;
}

/* calculates D given A, B, C, magnitude of CD, angle BCD (beta) and the
   dihedral angle between ABCD (Chi) and stores in vD */
void CalcNextCoOrd(vec vA,vec vB,vec vC,vec vD,FloatLo MagCD,FloatLo Beta,FloatLo Chi)
{
        vec vAxis,vTmp,vBAA,vAB,vBC,vAA,vnegB;
        FloatLo sinChi,cosChi,d,sina,cosa,axisx;

        VecSub(vAB,vB,vA);
        VecSub(vBC,vC,vB);
        Normalize(vAxis,vBC);
        VecScale(vTmp,vAxis,Dot(vAB,vAxis));
        VecSub(vBAA,vAB,vTmp);
        /* vBAA is really BA' */
        NegateVec(vBAA,vBAA);
        /* now get absolute co-ords of A' */  
        VecAdd(vAA,vBAA,vB);
        /* and rotate CCW about BC by Chi to get D' */
        /* stolen from rotate.c */
        sinChi=sin(Chi*DEGTORAD);
        cosChi=cos(Chi*DEGTORAD);
        d=sqrt(vAxis[1]*vAxis[1]+vAxis[2]*vAxis[2]);
        /* if d=0, already on x-axis */
        if (d==0) {
                sina=0;
                cosa=1;
        }
        else {
                sina=vAxis[1]/d;  
                cosa=vAxis[2]/d; 
        }
        axisx=vAxis[0];
        NegateVec(vnegB,vB);
        /* must translate -vB to origin */
        Translate(vAA,vnegB);
        Normalize(vAA,vAA);
        /* then map vector to the z-axis */
        Rotatecos(X_AXIS,sina,cosa,vAA,vAA);
        Rotatecos(Y_AXIS,-axisx,d,vAA,vAA);
        /* rotate Chi about the new z-axis */
        Rotatecos(Z_AXIS,sinChi,cosChi,vAA,vAA);
        /* return vector to its original position */
        Rotatecos(Y_AXIS,axisx,d,vAA,vAA);
        Rotatecos(X_AXIS,-sina,cosa,vAA,vAA);
        /* determine correct length for D' position */
        VecScale(vAA,vAA,MagCD*sin(Beta*DEGTORAD)); 
        /* and undo the translation now placing the origin at C */
        Translate(vAA,vC);
        /* now vAA is at D' so to get D, just add (-cos Beta)*|CD|vAxis */
        VecScale(vTmp,vAxis,-MagCD*cos(Beta*DEGTORAD));
        VecAdd(vD,vAA,vTmp);
}

/* calculates D given A, B, C, magnitude of CD, angle BCD (beta) and the
   dihedral angle between ABCD (Chi) and stores in vD */
/*void iCalcNextCoOrd(vec vA,vec vB,vec vC,vec vD,FloatLo MagCD,FloatLo
Beta,FloatLo Chi)
{
	ivec ivAxis,ivTmp,ivBAA,ivAB,ivBC,ivAA,ivnegB;
	Int4 sinChi,cosChi,d,sina,cosa,axisx;
	Int4 iChi,iBeta,iCD;
	ivec ivA,ivB,ivC,ivD;

	ivA[0]=vA[0]*0x10000L;
	ivA[1]=vA[1]*0x10000L;
	ivA[2]=vA[2]*0x10000L;
	ivB[0]=vB[0]*0x10000L;
	ivB[1]=vB[1]*0x10000L;
	ivB[2]=vB[2]*0x10000L;
	ivC[0]=vC[0]*0x10000L;
	ivC[1]=vC[1]*0x10000L;
	ivC[2]=vC[2]*0x10000L;
	ivD[0]=vD[0]*0x10000L;
	ivD[1]=vD[1]*0x10000L;
	ivD[2]=vD[2]*0x10000L;
	iChi=Chi*0x10000;
	iBeta=Beta*0x10000;
	iCD=MagCD*0x10000;
	iVecSub(ivAB,ivB,ivA);
	iVecSub(ivBC,ivC,ivB);
	iNormalize(ivAxis,ivBC);
	iVecScale(ivTmp,ivAxis,iDot(ivAB,ivAxis));
	iVecSub(ivBAA,ivAB,ivTmp);
	iNegateVec(ivBAA,ivBAA);
	iVecAdd(ivAA,ivBAA,ivB);
	sinChi=(sin((FloatLo)(iChi/0x10000)*DEGTORAD))*0x10000L;
	cosChi=(cos((FloatLo)(iChi/0x10000)*DEGTORAD))*0x10000L;
	d=((Int4)sqrt((ivAxis[1]/0x100)*(ivAxis[1]/0x100)+(ivAxis[2]/0x100)*(ivAxis[2]/0x100)))*0x100;
	if (d==0) {
		sina=0;
		cosa=1<<16;
	}
	else {
		sina=(ivAxis[1]*0x100L/d)*0x100;
		cosa=(ivAxis[2]*0x100L/d)*0x100;
	}
	axisx=ivAxis[0];
	iNegateVec(ivnegB,ivB);
        iTranslate(ivAA,ivnegB);
	iNormalize(ivAA,ivAA);
        iRotatecos(X_AXIS,sina,cosa,ivAA,ivAA);
        iRotatecos(Y_AXIS,-axisx,d,ivAA,ivAA);
        iRotatecos(Z_AXIS,sinChi,cosChi,ivAA,ivAA);
        iRotatecos(Y_AXIS,axisx,d,ivAA,ivAA);
        iRotatecos(X_AXIS,-sina,cosa,ivAA,ivAA);
	iVecScale(ivAA,ivAA,iCD*sin((iBeta/0x10000)*DEGTORAD));
        iTranslate(ivAA,ivC);
	iVecScale(ivTmp,ivAxis,-iCD*cos((iBeta/0x10000)*DEGTORAD));
	iVecAdd(ivD,ivAA,ivTmp);
	vD[0]=(FloatLo)ivD[0]/0x10000L;
	vD[1]=(FloatLo)ivD[1]/0x10000L;
	vD[2]=(FloatLo)ivD[2]/0x10000L;
}
*/

TrajErr SetCoOrds(PALD paldHere,vec vHere)
{
	if (paldHere==NULL) {
		ErrPostEx(SEV_ERROR,2,2,"Cannot set co-ordinates\n");
		return ERR_FAIL;
	}
	(paldHere->pflvData)[0]=vHere[0];
	(paldHere->pflvData)[1]=vHere[1];
	(paldHere->pflvData)[2]=vHere[2];
	(paldHere->pflvData)[3]=1.0;
	(paldHere->pflvData)[4]=0.0;
/*printf("%s ",((PMAD)(paldHere->pfbParent))->pcAName);
PrintVec(vHere);*/
	return ERR_SUCCESS;
}

void PlaceAtom(CharPtr AtomName,Byte rotnum,Int2 ref,FloatLo bondlength,FloatLo bondangle,FloatLo chi)
{
	PMAD pmadHere;
	PALD paldNew;

	pmadHere=FindAtomName(pvnmaHead,AtomName);
	paldNew=GetAtomLocs(pmadHere,rotnum);
	if (paldNew==NULL)
		paldNew=AllocNewLoc(pmadHere,rotnum);
	CalcNextCoOrd(vRef[ref],vRef[ref+1],vRef[ref+2],vRef[ref+3],bondlength,bondangle,chi);
	SetCoOrds(paldNew,vRef[ref+3]);
}

/* takes a PMGD of residue type ResID and allocates a new PALD, assigning it
   co-ordinates using standard bond angles, etc. and given chi1 and chi2
   dihedral angles.  The desired CAlpha position is given as well */
void Buildit(PMGD pmgdThis, vec vCAlpha, Char ResID, FloatLo chi1,
FloatLo chi2, Byte rotnum)
{
	PMAD pmadHere;
	PALD paldNew;
	vec vNew,vTmp1,vTmp2;
	vec vHAdir;
	vec vCBdir={0,0,1};
	vec vC;
	Int2 ResIUPAC;
	Char ResIDString[2];
	FloatLo cosb,sinb,thedot;
	FloatLo blnca,blcac,bancac;
	FloatLo tx,ty,absmag,del,chithis1,chithis2,err1,err2;
	FloatLo a1,a2,b1,b2,c1,d1,e1,f1,m1,m2,m3,m4,m5;

	VecScale(vRef[1],vCAlpha,1.0);
	ResIDString[0]=ResID;
	ResIDString[1]=0;
	ResIUPAC=StringCSpn(NCBIstdaaUC,ResIDString)-1;
	pvnmaHead=pmgdThis->pvnmaAHead;
	blnca=BL_NCA;
	blcac=BL_CAC;
	bancac=BA_NCAC;
	if (ResID=='G') {
		bancac=BA_G_NCAC;
		blcac=BL_G_CAC;
		blnca=BL_G_NCA;
	}
	if (ResID=='P') {
		bancac=BA_P_NCAC;
		blnca=BL_P_NCA;
	}
	pmadHere=FindAtomName(pvnmaHead," CA ");
	paldNew=GetAtomLocs(pmadHere,rotnum);
        SetCoOrds(paldNew,vCAlpha);
	vRef[0][0]=0;
	vRef[0][1]=(-blnca)*sin(ba_ncacb[ResIUPAC]*DEGTORAD);
	vRef[0][2]=blnca*cos(ba_ncacb[ResIUPAC]*DEGTORAD);
	vC[2]=blcac*cos(ba_ccacb[ResIUPAC]*DEGTORAD);
	vC[1]=cos(bancac*DEGTORAD)*blcac-vC[2]*sin(ba_ncacb[ResIUPAC]*DEGTORAD);
	vC[1]/=cos(ba_ncacb[ResIUPAC]*DEGTORAD);
	vC[0]=sqrt(blcac*blcac-vC[1]*vC[1]-vC[2]*vC[2]);	
	Normalize(vTmp1,vC);
	Normalize(vTmp2,vRef[0]);
	VecAdd(vNew,vTmp1,vTmp2);
	Normalize(vNew,vNew);
	Cross(vTmp1,vCBdir,vNew);
	Normalize(vNew,vTmp1);
	cosb=vNew[0];
	sinb=-sqrt(1-cosb*cosb);
        Rotatecos(Z_AXIS,sinb,cosb,vRef[0],vRef[0]);
        Rotatecos(Z_AXIS,sinb,cosb,vC,vC);
	Normalize(vTmp1,vC);
	Normalize(vTmp2,vRef[0]);
	VecAdd(vHAdir,vTmp1,vTmp2);
	Translate(vRef[0],vRef[1]);
	pmadHere=FindAtomName(pvnmaHead," N  ");
	paldNew=AllocNewLoc(pmadHere,rotnum);
        SetCoOrds(paldNew,vRef[0]);
	Translate(vC,vRef[1]);
	pmadHere=FindAtomName(pvnmaHead," C  ");
	paldNew=AllocNewLoc(pmadHere,rotnum);
        SetCoOrds(paldNew,vC);
        if (ResID!='G') {
		pmadHere=FindAtomName(pvnmaHead," CB ");
		paldNew=AllocNewLoc(pmadHere,rotnum);
		VecScale(vNew,vCBdir,bl_cacb[ResIUPAC]);
	        VecAdd(vRef[2],vNew,vRef[1]);
        	SetCoOrds(paldNew,vRef[2]);
	}
	else {
		pmadHere=FindAtomName(pvnmaHead,"2HA ");
		paldNew=AllocNewLoc(pmadHere,rotnum);
		VecScale(vNew,vCBdir,BL_G_CA2HA);
	        VecAdd(vRef[2],vNew,vCAlpha);
        	SetCoOrds(paldNew,vRef[2]);
	}
	Normalize(vTmp1,vNew);
	VecAdd(vHAdir,vHAdir,vTmp1);
	NegateVec(vHAdir,vHAdir);
	Normalize(vHAdir,vHAdir);
	if ((pmadHere=FindAtomName(pvnmaHead," HA "))==NULL)
		pmadHere=FindAtomName(pvnmaHead,"1HA ");
        paldNew=AllocNewLoc(pmadHere,rotnum);
       	VecScale(vNew,vHAdir,BL_CH);
        VecAdd(vNew,vNew,vRef[1]);
        SetCoOrds(paldNew,vNew);
	switch (ResID) {
		case 'A':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,180);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,60);
			PlaceAtom("3HB ",rotnum,0,BL_CH,BA_XCH,-60);
			break;
		case 'C':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" SG ",rotnum,0,BL_C_CBSG,BA_C_CACBSG,chi1);
			PlaceAtom(" HG ",rotnum,1,BL_SH,108,100);
			break;
		case 'D':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" CG ",rotnum,0,BL_D_CBCG,BA_D_CACBCG,chi1);
			PlaceAtom(" OD1",rotnum,1,BL_D_CGOD1,BA_D_CBCGOD1,chi2);
			PlaceAtom(" OD2",rotnum,1,BL_D_CGOD2,BA_D_CBCGOD2,chi2+180);
			break;
		case 'E':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" CG ",rotnum,0,BL_E_CBCG,BA_E_CACBCG,chi1);
			PlaceAtom("1HG ",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom("2HG ",rotnum,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom(" CD ",rotnum,1,BL_E_CGCD,BA_E_CBCGCD,chi2);
			PlaceAtom(" OE1",rotnum,2,BL_E_CDOE1,BA_E_CGCDOE1,CHI3_E);
			PlaceAtom(" OE2",rotnum,2,BL_E_CDOE2,BA_E_CGCDOE2,CHI3_E+180);
			break;
		case 'F':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" CG ",rotnum,0,BL_F_CBCG,BA_F_CACBCG,chi1);
			PlaceAtom(" CD1",rotnum,1,BL_F_CGCD1,BA_F_CBCGCD1,chi2);
			PlaceAtom(" HD1",rotnum,2,BL_CH,BA_ARCH,0);
			PlaceAtom(" CE1",rotnum,2,BL_F_CD1CE1,BA_F_CGCD1CE1,180);
			PlaceAtom(" HE1",rotnum,3,BL_CH,BA_ARCH,180);
			PlaceAtom(" CD2",rotnum,1,BL_F_CGCD2,BA_F_CBCGCD2,chi2+180);
			PlaceAtom(" HD2",rotnum,2,BL_CH,BA_ARCH,0);
			PlaceAtom(" CE2",rotnum,2,BL_F_CD2CE2,BA_F_CGCD2CE2,180);
			PlaceAtom(" HE2",rotnum,3,BL_CH,BA_ARCH,180);
			PlaceAtom(" CZ ",rotnum,3,BL_F_CE2CZ,BA_F_CD2CE2CZ,0);
			PlaceAtom(" HZ ",rotnum,4,BL_CH,BA_ARCH,180);
			break;
		case 'H':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" CG ",rotnum,0,BL_H_CBCG,BA_H_CACBCG,chi1);
			PlaceAtom(" ND1",rotnum,1,BL_H_CGND1,BA_H_CBCGND1,chi2);
			PlaceAtom(" CE1",rotnum,2,BL_H_ND1CE1,BA_H_CGND1CE1,180);
			PlaceAtom(" HE1",rotnum,3,BL_CH,BA_HISCH,180);
			PlaceAtom(" CD2",rotnum,1,BL_H_CGCD2,BA_H_CBCGCD2,chi2+180);
			PlaceAtom(" HD2",rotnum,2,BL_CH,BA_HISCH,0);
			PlaceAtom(" NE2",rotnum,2,BL_H_CD2NE2,BA_H_CGCD2NE2,180);
			PlaceAtom(" HE2",rotnum,3,BL_NH,BA_HISCH,180);
			break;
		case 'I':
			PlaceAtom(" HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" CG1",rotnum,0,BL_I_CBCG1,BA_I_CACBCG1,chi1);
			PlaceAtom("1HG1",rotnum,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG1",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom(" CD1",rotnum,1,BL_I_CG1CD1,BA_I_CBCG1CD1,chi2);
			PlaceAtom("1HD1",rotnum,2,BL_CH,BA_XCH,60);
			PlaceAtom("2HD1",rotnum,2,BL_CH,BA_XCH,-60);
			PlaceAtom("3HD1",rotnum,2,BL_CH,BA_XCH,180);
			PlaceAtom(" CG2",rotnum,0,BL_I_CBCG2,BA_I_CACBCG2,chi1-120);
			PlaceAtom("1HG2",rotnum,1,BL_CH,BA_XCH,60);
			PlaceAtom("2HG2",rotnum,1,BL_CH,BA_XCH,-60);
			PlaceAtom("3HG2",rotnum,1,BL_CH,BA_XCH,180);
			break;
		case 'K':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_K_CBCG,BA_K_CACBCG,chi1);
			PlaceAtom("1HG ",rotnum,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom(" CD ",rotnum,1,BL_K_CGCD,BA_K_CBCGCD,chi2);
			PlaceAtom("1HD ",rotnum,2,BL_CH,BA_XCH,CHI3_K+120);
			PlaceAtom("2HD ",rotnum,2,BL_CH,BA_XCH,CHI3_K-120);
			PlaceAtom(" CE ",rotnum,2,BL_K_CDCE,BA_K_CGCDCE,CHI3_K);
			PlaceAtom("1HE ",rotnum,3,BL_CH,BA_XCH,CHI4_K+120);
			PlaceAtom("2HE ",rotnum,3,BL_CH,BA_XCH,CHI4_K-120);
			PlaceAtom(" NZ ",rotnum,3,BL_K_CENZ,BA_K_CDCENZ,CHI4_K);
			PlaceAtom("1HZ ",rotnum,4,BL_NH3,BA_XCH,60);
			PlaceAtom("2HZ ",rotnum,4,BL_NH3,BA_XCH,-60);
			PlaceAtom("3HZ ",rotnum,4,BL_NH3,BA_XCH,180);
			break;
		case 'L':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_L_CBCG,BA_L_CACBCG,chi1);
			PlaceAtom(" HG ",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom(" CD1",rotnum,1,BL_L_CGCD1,BA_L_CBCGCD1,chi2);
			PlaceAtom("1HD1",rotnum,2,BL_CH,BA_XCH,60);
			PlaceAtom("2HD1",rotnum,2,BL_CH,BA_XCH,-60);
			PlaceAtom("3HD1",rotnum,2,BL_CH,BA_XCH,180);
			PlaceAtom(" CD2",rotnum,1,BL_L_CGCD2,BA_L_CBCGCD2,chi2+120);
			PlaceAtom("1HD2",rotnum,2,BL_CH,BA_XCH,0);
			PlaceAtom("2HD2",rotnum,2,BL_CH,BA_XCH,120);
			PlaceAtom("3HD2",rotnum,2,BL_CH,BA_XCH,-120);
			break;
		case 'M':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_M_CBCG,BA_M_CACBCG,chi1);
			PlaceAtom("1HG ",rotnum,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom(" SD ",rotnum,1,BL_M_CGSD,BA_M_CBCGSD,chi2);
			PlaceAtom(" CE ",rotnum,2,BL_M_SDCE,BA_M_CGSDCE,CHI3_M);
			PlaceAtom("1HE ",rotnum,3,BL_CH,BA_XCH,60);
			PlaceAtom("2HE ",rotnum,3,BL_CH,BA_XCH,-60);
			PlaceAtom("3HE ",rotnum,3,BL_CH,BA_XCH,180);
			break;
		case 'N':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_N_CBCG,BA_N_CACBCG,chi1);
			PlaceAtom(" OD1",rotnum,1,BL_N_CGOD1,BA_N_CBCGOD1,chi2);
			PlaceAtom(" ND2",rotnum,1,BL_N_CGND2,BA_N_CBCGND2,chi2+180);
			PlaceAtom("1HD2",rotnum,2,BL_NH,BA_XNH,0);
			PlaceAtom("2HD2",rotnum,2,BL_NH,BA_XNH,180);
			break;
		case 'P':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,117.5F);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,-121.9F);
			PlaceAtom(" CG ",rotnum,0,BL_P_CBCG,BA_P_CACBCG,chi1);
			PlaceAtom("1HG ",rotnum,1,BL_CH,BA_XCH,117.4F);
			PlaceAtom("2HG ",rotnum,1,BL_CH,BA_XCH,-122.0F);
			a1=BL_P_NCD*cos(DEGTORAD*BA_P_CANCD)*cos(DEGTORAD*(90-ba_ncacb[ResIUPAC]))-cos(DEGTORAD*(ba_ncacb[ResIUPAC]-90))*BL_P_NCA;
			a2=BL_P_NCD*sin(DEGTORAD*BA_P_CANCD)*sin(DEGTORAD*(90-ba_ncacb[ResIUPAC]));
			b1=BL_P_NCD*cos(DEGTORAD*BA_P_CANCD)*sin(DEGTORAD*(90-ba_ncacb[ResIUPAC]))+bl_cacb[ResIUPAC]+sin(DEGTORAD*(ba_ncacb[ResIUPAC]-90))*BL_P_NCA;
			b2=-BL_P_NCD*sin(DEGTORAD*BA_P_CANCD)*cos(DEGTORAD*(90-ba_ncacb[ResIUPAC]));
			c1=-BL_P_NCD*sin(DEGTORAD*BA_P_CANCD);
			d1=-BL_P_CBCG*cos(DEGTORAD*(BA_P_CACBCG-90));
			e1=-BL_P_CBCG*sin(DEGTORAD*(BA_P_CACBCG-90));
			f1=-d1;
			m1=-2*b2*(b1-e1)-2*a1*a2;
			m2=2*a2*d1;
			m3=2*c1*f1;
			m4=a1*a1+(b1-e1)*(b1-e1)-BL_P_CGCD*BL_P_CGCD+d1*d1+c1*c1;
			m5=-2*a1*d1;
			tx=m1+m2*cos(DEGTORAD*chi1);
			ty=m3*sin(-DEGTORAD*chi1);
			absmag=sqrt(tx*tx+ty*ty);
			thedot=ty/absmag;
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			del=acos(thedot);
			if ((sin(del)*(tx/absmag))<0)
				del=2*PI-del;
			if (m4+m5*cos(DEGTORAD*chi1)>absmag) {
				ErrPostEx(SEV_FATAL,1,1,"Cannot build PRO with given Chi 1 angle (%f degrees)",chi1);
				return;
			}				
			chithis1=asin((m4+m5*cos(DEGTORAD*chi1))/absmag);
			chithis2=PI-chithis1;
			chithis1-=del;
			chithis2-=del;
			chithis1*=-RADTODEG;
			chithis2*=-RADTODEG;
			pmadHere=FindAtomName(pvnmaHead," CD ");
			paldNew=AllocNewLoc(pmadHere,rotnum);
			CalcNextCoOrd(vRef[2],vRef[1],vRef[0],vRef[4],BL_P_NCD,BA_P_CANCD,chithis1);
			VecSub(vTmp1,vRef[2],vRef[3]);
			VecSub(vTmp2,vRef[4],vRef[3]);
			thedot=Dot(vTmp1,vTmp2)/(getMag(vTmp1)*getMag(vTmp2));
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			err1=fabs(BA_P_CBCGCD-RADTODEG*acos(thedot));
			VecSub(vTmp1,vRef[3],vRef[4]);
			VecSub(vTmp2,vRef[0],vRef[4]);
			thedot=Dot(vTmp1,vTmp2)/(getMag(vTmp1)*getMag(vTmp2));
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			err1+=fabs(BA_P_CGCDN-RADTODEG*acos(thedot));
			CalcNextCoOrd(vRef[2],vRef[1],vRef[0],vRef[4],BL_P_NCD,BA_P_CANCD,chithis2);
			VecSub(vTmp1,vRef[2],vRef[3]);
			VecSub(vTmp2,vRef[4],vRef[3]);
			thedot=Dot(vTmp1,vTmp2)/(getMag(vTmp1)*getMag(vTmp2));
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			err2=fabs(BA_P_CBCGCD-RADTODEG*acos(thedot));
			VecSub(vTmp1,vRef[3],vRef[4]);
			VecSub(vTmp2,vRef[0],vRef[4]);
			thedot=Dot(vTmp1,vTmp2)/(getMag(vTmp1)*getMag(vTmp2));
			if (thedot>1.0) thedot=1.0;
			if (thedot<-1.0) thedot=-1.0;
			err2+=fabs(BA_P_CGCDN-RADTODEG*acos(thedot));
			printf("PRO rotamer %d total ring interior angle error: ",rotnum);
			if (err2>err1) {
				CalcNextCoOrd(vRef[2],vRef[1],vRef[0],vRef[4],BL_P_NCD,BA_P_CANCD,chithis1);
				printf("%f degrees\n",err1);
			}
			else
				printf("%f degrees\n",err2);
			SetCoOrds(paldNew,vRef[4]);
			PlaceAtom("1HD ",rotnum,2,BL_CH,BA_XCH,125.2F);
			PlaceAtom("2HD ",rotnum,2,BL_CH,BA_XCH,-112.5F);
			break;
		case 'Q':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_Q_CBCG,BA_Q_CACBCG,chi1);
			PlaceAtom("1HG ",rotnum,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom(" CD ",rotnum,1,BL_Q_CGCD,BA_Q_CBCGCD,chi2);
			PlaceAtom(" OE1",rotnum,2,BL_Q_CDOE1,BA_Q_CGCDOE1,CHI3_Q);
			PlaceAtom(" NE2",rotnum,2,BL_Q_CDNE2,BA_Q_CGCDNE2,CHI3_Q+180);
			PlaceAtom("1HE2",rotnum,3,BL_NH,BA_XNH,0);
			PlaceAtom("2HE2",rotnum,3,BL_NH,BA_XNH,180);
			break;
		case 'R':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_R_CBCG,BA_R_CACBCG,chi1);
			PlaceAtom("1HG ",rotnum,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",rotnum,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom(" CD ",rotnum,1,BL_R_CGCD,BA_R_CBCGCD,chi2);
			PlaceAtom("1HD ",rotnum,2,BL_CH,BA_XCH,CHI3_R+120);
			PlaceAtom("2HD ",rotnum,2,BL_CH,BA_XCH,CHI3_R-120);
			PlaceAtom(" NE ",rotnum,2,BL_R_CDNE,BA_R_CGCDNE,CHI3_R);
			PlaceAtom(" HE ",rotnum,3,BL_NH,BA_XNH,CHI4_R+180);
			PlaceAtom(" CZ ",rotnum,3,BL_R_NECZ,BA_R_CDNECZ,CHI4_R);
			PlaceAtom(" NH2",rotnum,4,BL_R_CZNH2,BA_R_NECZNH2,CHI5_R);
			PlaceAtom("1HH2",rotnum,5,BL_GUANNH,BA_XNH,0);
			PlaceAtom("2HH2",rotnum,5,BL_GUANNH,BA_XNH,180);
			PlaceAtom(" NH1",rotnum,4,BL_R_CZNH1,BA_R_NECZNH1,180+CHI5_R);
			PlaceAtom("1HH1",rotnum,5,BL_GUANNH,BA_XNH,0);
			PlaceAtom("2HH1",rotnum,5,BL_GUANNH,BA_XNH,180);
			break;
		case 'S':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" OG ",rotnum,0,BL_S_CBOG,BA_S_CACBOG,chi1);
			PlaceAtom(" HG ",rotnum,1,BL_OH,BA_XOH,60);
			break;
		case 'T':
			PlaceAtom(" HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom(" OG1",rotnum,0,BL_T_CBOG1,BA_T_CACBOG1,chi1);
			PlaceAtom(" HG1",rotnum,1,BL_OH,BA_XOH,120);
			PlaceAtom(" CG2",rotnum,0,BL_T_CBCG2,BA_T_CACBCG2,chi1-120);
			PlaceAtom("1HG2",rotnum,1,BL_CH,BA_XCH,60);
			PlaceAtom("2HG2",rotnum,1,BL_CH,BA_XCH,-60);
			PlaceAtom("3HG2",rotnum,1,BL_CH,BA_XCH,180);
			break;
		case 'V':
			PlaceAtom(" HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG1",rotnum,0,BL_V_CBCG1,BA_V_CACBCG1,chi1);
			PlaceAtom("1HG1",rotnum,1,BL_CH,BA_XCH,60);
			PlaceAtom("2HG1",rotnum,1,BL_CH,BA_XCH,-60);
			PlaceAtom("3HG1",rotnum,1,BL_CH,BA_XCH,180);
			PlaceAtom(" CG2",rotnum,0,BL_V_CBCG2,BA_V_CACBCG2,chi1+120);
			PlaceAtom("1HG2",rotnum,1,BL_CH,BA_XCH,0);
			PlaceAtom("2HG2",rotnum,1,BL_CH,BA_XCH,-120);
			PlaceAtom("3HG2",rotnum,1,BL_CH,BA_XCH,120);
			break;
		case 'W':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_W_CBCG,BA_W_CACBCG,chi1);
			PlaceAtom(" CD1",rotnum,1,BL_W_CGCD1,BA_W_CBCGCD1,chi2);
			PlaceAtom(" HD1",rotnum,2,BL_CH,BA_HISCH,0);
			PlaceAtom(" NE1",rotnum,2,BL_W_CD1NE1,BA_W_CGCD1NE1,180);
			PlaceAtom(" HE1",rotnum,3,BL_NH,BA_HISCH,180);
			PlaceAtom(" CD2",rotnum,1,BL_W_CGCD2,BA_W_CBCGCD2,chi2+180);
			PlaceAtom(" CE2",rotnum,2,BL_W_CD2CE2,BA_W_CGCD2CE2,180);
			PlaceAtom(" CZ2",rotnum,3,BL_W_CE2CZ2,BA_W_CD2CE2CZ2,180);
			PlaceAtom(" HZ2",rotnum,4,BL_CH,BA_ARCH,180);
			PlaceAtom(" CH2",rotnum,4,BL_W_CZ2CH2,BA_W_CE2CZ2CH2,0);
			PlaceAtom(" HH2",rotnum,5,BL_CH,BA_ARCH,180);
			PlaceAtom(" CE3",rotnum,2,BL_W_CD2CE3,BA_W_CGCD2CE3,0);
			PlaceAtom(" HE3",rotnum,3,BL_CH,BA_ARCH,0);
			PlaceAtom(" CZ3",rotnum,3,BL_W_CE3CZ3,BA_W_CD2CE3CZ3,180);
			PlaceAtom(" HZ3",rotnum,4,BL_CH,BA_ARCH,180);
			break;
		case 'Y':
			PlaceAtom("1HB ",rotnum,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",rotnum,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom(" CG ",rotnum,0,BL_Y_CBCG,BA_Y_CACBCG,chi1);
			PlaceAtom(" CD1",rotnum,1,BL_Y_CGCD1,BA_Y_CBCGCD1,chi2);
			PlaceAtom(" HD1",rotnum,2,BL_CH,BA_ARCH,0);
			PlaceAtom(" CE1",rotnum,2,BL_Y_CD1CE1,BA_Y_CGCD1CE1,180);
			PlaceAtom(" HE1",rotnum,3,BL_CH,BA_ARCH,180);
			PlaceAtom(" CD2",rotnum,1,BL_Y_CGCD2,BA_Y_CBCGCD2,chi2+180);
			PlaceAtom(" HD2",rotnum,2,BL_CH,BA_ARCH,0);
			PlaceAtom(" CE2",rotnum,2,BL_Y_CD2CE2,BA_Y_CGCD2CE2,180);
			PlaceAtom(" HE2",rotnum,3,BL_CH,BA_ARCH,180);
			PlaceAtom(" CZ ",rotnum,3,BL_Y_CE2CZ,BA_Y_CD2CE2CZ,0);
			PlaceAtom(" OH ",rotnum,4,BL_Y_CZOH,BA_Y_CE2CZOH,180);
			PlaceAtom(" HH ",rotnum,5,BL_OH,BA_XOH,60);
			break;
		default:;
	}
}

/* takes a PMGD of residue type ResID and recomputes H locations using
   ideal geometry for a given model */
TrajErr RebuildH(PMGD pmgdThis,Int2 Model)
{
	PMAD pmadHere;
	PALD paldNew;
	vec vNew,vNew1,vNew2;
	vec vHAdir;
	vec vCBdir;
	vec vZero={0,0,0};
	vec vC,vCrel,vNrel,vH;
	Char ResID;
	FloatLo chi1,chi2,chi3,chi4;
	PRS prshere;

	ResID=pmgdThis->pcIUPAC[0];
	pvnmaHead=pmgdThis->pvnmaAHead;
	prshere=(PRS)MemNew(sizeof(RS));
	if (MeasureResidueChis(pmgdThis,&prshere,Model,FALSE,0,NULL)!=ERR_SUCCESS) {
		prshere=MemFree(prshere);
		return ERR_FAIL;
	}
	chi1=prshere->Chi1;
	chi2=prshere->Chi2;
	chi3=prshere->Chi3;
	chi4=prshere->Chi4;
	prshere=MemFree(prshere);
	if (GetCoOrds(pmgdThis," N  ",vZero,vRef[0],Model)==NULL)
		return ERR_FAIL;
	if (GetCoOrds(pmgdThis," CA ",vZero,vRef[1],Model)==NULL)
		return ERR_FAIL;
	if (GetCoOrds(pmgdThis," C  ",vZero,vC,Model)==NULL)
		return ERR_FAIL;
	if (GetCoOrds(pmgdThis," N  ",vRef[1],vNrel,Model)==NULL)
		return ERR_FAIL;
	if (GetCoOrds(pmgdThis," C  ",vRef[1],vCrel,Model)==NULL)
		return ERR_FAIL;
	Normalize(vNrel,vNrel);
	Normalize(vCrel,vCrel);
    if (ResID!='G') {
		if (GetCoOrds(pmgdThis," CB ",vZero,vRef[2],Model)==NULL)
			return ERR_FAIL;
		if (GetCoOrds(pmgdThis," CB ",vRef[1],vCBdir,Model)==NULL)
			return ERR_FAIL;
		Normalize(vCBdir,vCBdir);
	}
	else {
		Cross(vNew1,vNrel,vCrel);
		Normalize(vNew1,vNew1);
		VecAdd(vNew2,vNrel,vCrel);
		Normalize(vNew2,vNew2);
		VecScale(vNew2,vNew2,-1.0);
		VecScale(vNew1,vNew1,sqrt(8.0));		
		VecAdd(vCBdir,vNew1,vNew2);
		Normalize(vCBdir,vCBdir);
		/* set 2HA for Gly */
		pmadHere=FindAtomName(pvnmaHead,"2HA ");
		paldNew=GetAtomLocs(pmadHere,Model);
		VecScale(vNew,vCBdir,BL_G_CA2HA);
	    VecAdd(vRef[2],vNew,vRef[1]);
        SetCoOrds(paldNew,vRef[2]);
	}
	Cross(vNew1,vNrel,vCrel);
	Normalize(vNew1,vNew1);
	VecAdd(vNew2,vNrel,vCrel);
	Normalize(vNew2,vNew2);
	VecScale(vNew2,vNew2,-sqrt(8.0));
	VecScale(vNew1,vNew1,-1.0);
	VecAdd(vHAdir,vNew1,vNew2);
	Normalize(vHAdir,vHAdir);
	if ((pmadHere=FindAtomName(pvnmaHead," HA "))==NULL)
		pmadHere=FindAtomName(pvnmaHead,"1HA ");
    paldNew=GetAtomLocs(pmadHere,Model);
   	VecScale(vNew,vHAdir,BL_CH);
    VecAdd(vNew,vNew,vRef[1]);
    SetCoOrds(paldNew,vNew);
	switch (ResID) {
		case 'A':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,180);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,60);
			PlaceAtom("3HB ",Model,0,BL_CH,BA_XCH,-60);
			break;
		case 'C':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			if (GetCoOrds(pmgdThis," SG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HG ",Model,1,BL_SH,108,100);
			break;
		case 'D':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			break;
		case 'E':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG ",Model,1,BL_CH,BA_XCH,chi2-120);
			PlaceAtom("2HG ",Model,1,BL_CH,BA_XCH,chi2+120);
			break;
		case 'F':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CD1",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HD1",Model,2,BL_CH,BA_ARCH,0);
			if (GetCoOrds(pmgdThis," CE1",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE1",Model,3,BL_CH,BA_ARCH,180);
			if (GetCoOrds(pmgdThis," CD2",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HD2",Model,2,BL_CH,BA_ARCH,0);
			if (GetCoOrds(pmgdThis," CE2",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE2",Model,3,BL_CH,BA_ARCH,180);
			if (GetCoOrds(pmgdThis," CZ ",vZero,vRef[6],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HZ ",Model,4,BL_CH,BA_ARCH,180);
			break;
		case 'H':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," ND1",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CE1",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE1",Model,3,BL_CH,BA_HISCH,180);
			if (GetCoOrds(pmgdThis," CD2",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HD2",Model,2,BL_CH,BA_HISCH,0);
			if (GetCoOrds(pmgdThis," NE2",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE2",Model,3,BL_NH,BA_HISCH,180);
			break;
		case 'I':
			PlaceAtom(" HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			if (GetCoOrds(pmgdThis," CG1",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG1",Model,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG1",Model,1,BL_CH,BA_XCH,chi2-120);
			if (GetCoOrds(pmgdThis," CD1",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD1",Model,2,BL_CH,BA_XCH,60);
			PlaceAtom("2HD1",Model,2,BL_CH,BA_XCH,-60);
			PlaceAtom("3HD1",Model,2,BL_CH,BA_XCH,180);
			if (GetCoOrds(pmgdThis," CG2",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG2",Model,1,BL_CH,BA_XCH,60);
			PlaceAtom("2HG2",Model,1,BL_CH,BA_XCH,-60);
			PlaceAtom("3HG2",Model,1,BL_CH,BA_XCH,180);
			break;
		case 'K':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG ",Model,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",Model,1,BL_CH,BA_XCH,chi2-120);
			if (GetCoOrds(pmgdThis," CD ",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD ",Model,2,BL_CH,BA_XCH,chi3+120);
			PlaceAtom("2HD ",Model,2,BL_CH,BA_XCH,chi3-120);
			if (GetCoOrds(pmgdThis," CE ",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HE ",Model,3,BL_CH,BA_XCH,chi4+120);
			PlaceAtom("2HE ",Model,3,BL_CH,BA_XCH,chi4-120);
			if (GetCoOrds(pmgdThis," NZ ",vZero,vRef[6],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HZ ",Model,4,BL_NH3,BA_XCH,60);
			PlaceAtom("2HZ ",Model,4,BL_NH3,BA_XCH,-60);
			PlaceAtom("3HZ ",Model,4,BL_NH3,BA_XCH,180);
			break;
		case 'L':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HG ",Model,1,BL_CH,BA_XCH,chi2-120);
			if (GetCoOrds(pmgdThis," CD1",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD1",Model,2,BL_CH,BA_XCH,60);
			PlaceAtom("2HD1",Model,2,BL_CH,BA_XCH,-60);
			PlaceAtom("3HD1",Model,2,BL_CH,BA_XCH,180);
			if (GetCoOrds(pmgdThis," CD2",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD2",Model,2,BL_CH,BA_XCH,0);
			PlaceAtom("2HD2",Model,2,BL_CH,BA_XCH,120);
			PlaceAtom("3HD2",Model,2,BL_CH,BA_XCH,-120);
			break;
		case 'M':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG ",Model,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",Model,1,BL_CH,BA_XCH,chi2-120);
			if (GetCoOrds(pmgdThis," SD ",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CE ",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HE ",Model,3,BL_CH,BA_XCH,60);
			PlaceAtom("2HE ",Model,3,BL_CH,BA_XCH,-60);
			PlaceAtom("3HE ",Model,3,BL_CH,BA_XCH,180);
			break;
		case 'N':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," ND2",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD2",Model,2,BL_NH,BA_XNH,0);
			PlaceAtom("2HD2",Model,2,BL_NH,BA_XNH,180);
			break;
		case 'P':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,117.5F);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,-121.9F);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG ",Model,1,BL_CH,BA_XCH,117.4F);
			PlaceAtom("2HG ",Model,1,BL_CH,BA_XCH,-122.0F);
			if (GetCoOrds(pmgdThis," CD ",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD ",Model,2,BL_CH,BA_XCH,125.2F);
			PlaceAtom("2HD ",Model,2,BL_CH,BA_XCH,-112.5F);
			break;
		case 'Q':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG ",Model,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",Model,1,BL_CH,BA_XCH,chi2-120);
			if (GetCoOrds(pmgdThis," CD ",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," NE2",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HE2",Model,3,BL_NH,BA_XNH,0);
			PlaceAtom("2HE2",Model,3,BL_NH,BA_XNH,180);
			break;
		case 'R':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG ",Model,1,BL_CH,BA_XCH,chi2+120);
			PlaceAtom("2HG ",Model,1,BL_CH,BA_XCH,chi2-120);
			if (GetCoOrds(pmgdThis," CD ",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HD ",Model,2,BL_CH,BA_XCH,chi3+120);
			PlaceAtom("2HD ",Model,2,BL_CH,BA_XCH,chi3-120);
			if (GetCoOrds(pmgdThis," NE ",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE ",Model,3,BL_NH,BA_XNH,chi4+180);
			if (GetCoOrds(pmgdThis," CZ ",vZero,vRef[6],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," NH2",vZero,vRef[7],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HH2",Model,5,BL_GUANNH,BA_XNH,0);
			PlaceAtom("2HH2",Model,5,BL_GUANNH,BA_XNH,180);
			if (GetCoOrds(pmgdThis," NH1",vZero,vRef[7],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HH1",Model,5,BL_GUANNH,BA_XNH,0);
			PlaceAtom("2HH1",Model,5,BL_GUANNH,BA_XNH,180);
			break;
		case 'S':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," OG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HG ",Model,1,BL_OH,BA_XOH,60);
			break;
		case 'T':
			PlaceAtom(" HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			if (GetCoOrds(pmgdThis," OG1",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HG1",Model,1,BL_OH,BA_XOH,120);
			if (GetCoOrds(pmgdThis," CG2",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG2",Model,1,BL_CH,BA_XCH,60);
			PlaceAtom("2HG2",Model,1,BL_CH,BA_XCH,-60);
			PlaceAtom("3HG2",Model,1,BL_CH,BA_XCH,180);
			break;
		case 'V':
			PlaceAtom(" HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG1",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG1",Model,1,BL_CH,BA_XCH,60);
			PlaceAtom("2HG1",Model,1,BL_CH,BA_XCH,-60);
			PlaceAtom("3HG1",Model,1,BL_CH,BA_XCH,180);
			if (GetCoOrds(pmgdThis," CG2",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom("1HG2",Model,1,BL_CH,BA_XCH,0);
			PlaceAtom("2HG2",Model,1,BL_CH,BA_XCH,-120);
			PlaceAtom("3HG2",Model,1,BL_CH,BA_XCH,120);
			break;
		case 'W':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CD1",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HD1",Model,2,BL_CH,BA_HISCH,0);
			if (GetCoOrds(pmgdThis," NE1",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE1",Model,3,BL_NH,BA_HISCH,180);
			if (GetCoOrds(pmgdThis," CD2",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CE2",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CZ2",vZero,vRef[6],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HZ2",Model,4,BL_CH,BA_ARCH,180);
			if (GetCoOrds(pmgdThis," CH2",vZero,vRef[7],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HH2",Model,5,BL_CH,BA_ARCH,180);
			if (GetCoOrds(pmgdThis," CE3",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE3",Model,3,BL_CH,BA_ARCH,0);
			if (GetCoOrds(pmgdThis," CZ3",vZero,vRef[6],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HZ3",Model,4,BL_CH,BA_ARCH,180);
			break;
		case 'Y':
			PlaceAtom("1HB ",Model,0,BL_CH,BA_XCH,chi1+120);
			PlaceAtom("2HB ",Model,0,BL_CH,BA_XCH,chi1-120);
			if (GetCoOrds(pmgdThis," CG ",vZero,vRef[3],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," CD1",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HD1",Model,2,BL_CH,BA_ARCH,0);
			if (GetCoOrds(pmgdThis," CE1",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE1",Model,3,BL_CH,BA_ARCH,180);
			if (GetCoOrds(pmgdThis," CD2",vZero,vRef[4],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HD2",Model,2,BL_CH,BA_ARCH,0);
			if (GetCoOrds(pmgdThis," CE2",vZero,vRef[5],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HE2",Model,3,BL_CH,BA_ARCH,180);
			if (GetCoOrds(pmgdThis," CZ ",vZero,vRef[6],Model)==NULL)
				return ERR_FAIL;
			if (GetCoOrds(pmgdThis," OH ",vZero,vRef[7],Model)==NULL)
				return ERR_FAIL;
			PlaceAtom(" HH ",Model,5,BL_OH,BA_XOH,60);
			break;
		default:;
	}
	if (pmgdThis->pdnmgLink->choice==1) {
		/* put H3 part of NH3 */
		if (ResID!='P')
			CalcNextCoOrd(vC,vRef[1],vRef[0],vH,BL_NH3H,BA_CANH3H,180);
		else
			CalcNextCoOrd(vC,vRef[1],vRef[0],vH,BL_NH3H,BA_CANH3H,0);
		pmadHere=FindAtomName(pvnmaHead,"3H  ");
		paldNew=GetAtomLocs(pmadHere,Model);
		SetCoOrds(paldNew,vH);
		if (ResID!='P')
			CalcNextCoOrd(vC,vRef[1],vRef[0],vH,BL_NH3H,BA_CANH3H,-60);
		else
			CalcNextCoOrd(vC,vRef[1],vRef[0],vH,BL_NH3H,BA_CANH3H,-120);
		pmadHere=FindAtomName(pvnmaHead,"2H  ");
		paldNew=GetAtomLocs(pmadHere,Model);
		SetCoOrds(paldNew,vH);
		if (ResID!='P') {
			CalcNextCoOrd(vC,vRef[1],vRef[0],vH,BL_NH3H,BA_CANH3H,60);
			pmadHere=FindAtomName(pvnmaHead,"1H  ");
			paldNew=GetAtomLocs(pmadHere,Model);
			SetCoOrds(paldNew,vH);
		}
	}
	return ERR_SUCCESS;
}

TrajErr MeasureResidueChis(PMGD pmgd,PRS *pprshere,Int2 Model,Boolean onlyburied,FloatLo temp,DValNodePtr dvnpPotential)
{
	Char ResID;
	Int2 numchi,cnt;
	Char ChiAtoms[7][7];
	vec vZero={0.0,0.0,0.0};
	vec v1,v2,v3,v4;
	FloatLo ba,bl,dihed,v;
	FloatHi x;
	Int2 res;
	DValNodePtr pdnSubList;

	StringCpy(ChiAtoms[0]," N  ");
	StringCpy(ChiAtoms[1]," CA ");
	StringCpy(ChiAtoms[2]," CB ");
	/* gets parent AA if modified AA, only 'X' if unknown */
	ResID=GetAAFromIDict(pmgd);
	/* clear all to 0 by default */
	(*pprshere)->Chi1=0.0;
	(*pprshere)->Chi2=0.0;
	(*pprshere)->Chi3=0.0;
	(*pprshere)->Chi4=0.0;
	if (onlyburied) {
		if (dvnpPotential!=NULL) {
			/* in case sigmas were given, assume room temp. */
			if (temp<=0.0)
				temp=298.15;
			/* get residue potential in x */
			res=pmgd->pdnmgLink->choice;
			while (res>1) {
				dvnpPotential=dvnpPotential->next;
				res--;
			}
			pdnSubList = (DValNodePtr)(dvnpPotential->data.ptrvalue);
			x = ((AdjListNodePtr)(pdnSubList->data.ptrvalue))->potential;
			v = sqrt((0.5+atan(74.896*((FloatLo)x)/temp+PI/2.0)/PI)); /* 74.896 = 2 PI * T0 / qr, T0=298K, qr=25 */
			if (v>0.5) /* not buried */
				return ERR_SUCCESS;
		}
		else {
			ErrPostEx(SEV_FATAL,1,32,"Potential not computed for savechis=buried only, this is a bug, please contact authors at trades@mshri.on.ca");
			return ERR_FAIL;
		}
	}
	switch (ResID) {
		case 'A':
			numchi=0;
			break;
		case 'C':
			numchi=1;
			StringCpy(ChiAtoms[3]," SG ");
			break;
		case 'D':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," OD1");
			break;
		case 'E':
			numchi=3;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD ");
			StringCpy(ChiAtoms[5]," OE1");
			break;
		case 'F':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD1");
			break;
		case 'G':
			numchi=0;
			break;
		case 'H':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," ND1");
			break;
		case 'I':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG1");
			StringCpy(ChiAtoms[4]," CD1");
			break;
		case 'K':
			numchi=4;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD ");
			StringCpy(ChiAtoms[5]," CE ");
			StringCpy(ChiAtoms[6]," NZ ");
			break;
		case 'L':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD1");
			break;
		case 'M':
			numchi=3;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," SD ");
			StringCpy(ChiAtoms[5]," CE ");
			break;
		case 'N':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," OD1");
			break;
		case 'P':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD ");
			break;
		case 'Q':
			numchi=3;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD ");
			StringCpy(ChiAtoms[5]," OE1");
			break;
		case 'R':
			numchi=4;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD ");
			StringCpy(ChiAtoms[5]," NE ");
			StringCpy(ChiAtoms[6]," CZ ");
			break;
		case 'S':
			numchi=1;
			StringCpy(ChiAtoms[3]," OG ");
			break;
		case 'T':
			numchi=1;
			StringCpy(ChiAtoms[3]," OG1");
			break;
		case 'V':
			numchi=1;
			StringCpy(ChiAtoms[3]," CG1");
			break;
		case 'W':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD1");
			break;
		case 'Y':
			numchi=2;
			StringCpy(ChiAtoms[3]," CG ");
			StringCpy(ChiAtoms[4]," CD1");
			break;
		default:  /* includes 'X' == unknown */
			return ERR_FAIL;
	}
	if (GetCoOrds(pmgd,ChiAtoms[0],vZero,v1,Model)==NULL)
		return ERR_FAIL;
	if (GetCoOrds(pmgd,ChiAtoms[1],vZero,v2,Model)==NULL)
		return ERR_FAIL;
	if (GetCoOrds(pmgd,ChiAtoms[2],vZero,v3,Model)==NULL)
		return ERR_FAIL;
	for (cnt=0;cnt<numchi;cnt++) {
		if (GetCoOrds(pmgd,ChiAtoms[cnt+3],vZero,v4,Model)==NULL)
			return ERR_FAIL;
		/* use dummy ba and bl variables */
		GetDihedral(v1,v2,v3,v4,-180,&dihed,&ba,&ba,&bl,&bl,&bl);
		/* ensure not exactly zero */
		if (dihed==0.0)
			dihed=0.01;
		switch (cnt) {
			case 0:
				(*pprshere)->Chi1=dihed;
				break;
			case 1:
				if ((ResID=='Y') || (ResID=='F') || (ResID=='D')) {
					/* ensure convention is met */
					while (dihed<-90.0)
						dihed+=180.0;
					while (dihed>90.0)
						dihed-=180.0;
				}
				if (dihed==0.0)
					dihed=0.01;
				(*pprshere)->Chi2=dihed;
				break;
			case 2:
				(*pprshere)->Chi3=dihed;
				break;
			case 3:
				(*pprshere)->Chi4=dihed;
				break;
			default:;
		}
		/* shift vectors */
		VecAdd(v1,v2,vZero);  /* v1 = v2 +0 */
		VecAdd(v2,v3,vZero);  /* v2 = v3 +0 */
		VecAdd(v3,v4,vZero);  /* v3 = v4 +0 */
	}
	return ERR_SUCCESS;
}

ValNodePtr MeasureChis(PMMD pmmd,Int2 Model)
{
	PDNMG pdnmg;
	PMGD pmgd;
	PRS prstemp;
	pnewrotlibrecord prlrHere=NULL;
	ValNodePtr vnpHead=NULL,vnpLast=NULL,vnpHere;

	prstemp=(PRS)MemNew(sizeof(RS));
	pdnmg=pmmd->pdnmgHead;
	while (pdnmg) {
		pmgd=(PMGD)(pdnmg->data.ptrvalue);
		prlrHere=MemNew(sizeof(newrotlibrecord));
		/* build up valnode list of rotamer records */
		vnpHere=ValNodeAddPointer(&vnpLast,pdnmg->choice,(Pointer)prlrHere);
		if (vnpHead==NULL)
			vnpHead=vnpHere;
		vnpLast=vnpHere;
		/* now fill in prlr chi values */
		prlrHere->p=0.0;
		prlrHere->chi1sd=0.0;
		prlrHere->chi2sd=0.0;
		if (MeasureResidueChis(pmgd,&prstemp,Model,FALSE,0.0,NULL)==ERR_SUCCESS) {
			prlrHere->chi1=prstemp->Chi1;
			prlrHere->chi2=prstemp->Chi2;
			prlrHere->chi3=prstemp->Chi3;
			prlrHere->chi4=prstemp->Chi4;
		}			
		pdnmg=pdnmg->next;
	}
	prstemp=MemFree(prstemp);
	return vnpHead;
}

void GetChiFromRotid(PRS *pprsHere,Uint4 rotid)
{
	FloatHi fchi1=0.0,fchi2=0.0,fchi3=0.0,fchi4=0.0;
	Uint1 ichi1,ichi2,ichi3,ichi4;
	
	if (rotid==0)
		return;
	ichi1=(rotid%256);
	rotid>>=8;
	ichi2=(rotid%256);
	rotid>>=8;
	ichi3=(rotid%256);
	rotid>>=8;
	ichi4=(rotid%256);
	if (ichi1!=0)
		fchi1=(((FloatHi)(ichi1-1))/254.0)*360.01-180.0;
	if (ichi2!=0)
		fchi2=(((FloatHi)(ichi2-1))/254.0)*360.01-180.0;
	if (ichi3!=0)
		fchi3=(((FloatHi)(ichi3-1))/254.0)*360.01-180.0;
	if (ichi4!=0)
		fchi4=(((FloatHi)(ichi4-1))/254.0)*360.01-180.0;		
	(*pprsHere)->Chi1=fchi1;
	(*pprsHere)->Chi2=fchi2;
	(*pprsHere)->Chi3=fchi3;
	(*pprsHere)->Chi4=fchi4;
}

Uint4 ComputeRotid(PRS prsHere)
{
	FloatHi fchi1,fchi2,fchi3,fchi4;
	Uint1 ichi1,ichi2,ichi3,ichi4;
	Uint4 rotid;
	
	if (prsHere==NULL)
		return 0;
	fchi1=(FloatHi)prsHere->Chi1;
	fchi2=(FloatHi)prsHere->Chi2;
	fchi3=(FloatHi)prsHere->Chi3;
	fchi4=(FloatHi)prsHere->Chi4;
	/* convert to range [0,1), the [0,254) and add 1.5 to get [1,256),
		 then truncate */
	ichi1=(Uint1)(((fchi1+180.0)/360.01)*254.0+1.5);
	/* just in case */
	if (ichi1==0)
		ichi1=1;
	ichi2=(Uint1)(((fchi2+180.0)/360.01)*254.0+1.5);
	if (ichi2==0)
		ichi2=1;
	ichi3=(Uint1)(((fchi3+180.0)/360.01)*254.0+1.5);
	if (ichi3==0)
		ichi3=1;
	ichi4=(Uint1)(((fchi4+180.0)/360.01)*254.0+1.5);
	if (ichi4==0)
		ichi4=1;
	rotid=0;
	if (fchi4!=0.0)
		rotid+=ichi4;
	rotid<<=8;
	if (fchi3!=0.0)
		rotid+=ichi3;
	rotid<<=8;
	if (fchi2!=0.0)
		rotid+=ichi2;
	rotid<<=8;
	if (fchi1!=0.0)
		rotid+=ichi1;
	return rotid;		
}

ValNodePtr FreeChis(ValNodePtr vnp)
{
	ValNodeFreeData(vnp);
	return NULL;
}


/*  
$Log: buildit.c,v $
Revision 1.11  2003/11/13 20:05:24  feldman
Moved buildit to mmdbtraj lib

Revision 1.10  2003/05/02 16:43:47  feldman
Removed unneded purgeglobs dependency

Revision 1.9  2003/05/01 20:33:50  feldman
Added 'rebuild hydrogens' function

Revision 1.8  2003/04/04 21:54:04  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.7  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.6  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.5  2001/03/09 17:33:58  feldman
Added code to deal with rotamers.  Now rotamer chi angles can
be encoded into a 4-byte integer and used by foldtraj when
building the structures.  These can be edited by VisTraj as well

Revision 1.4  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.3  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.2  2000/07/06 15:29:37  feldman
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


