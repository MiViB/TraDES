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

#define MAXRAND 2147483648.0    /* 2^31 */

/* return a pseudo-random number from -1 to 1, uniform distribution */
FloatLo Rand1(void)
{
	return 2.0*((FloatLo)RandomNum()/MAXRAND)-1.0;
}

/* approximates Weibull distribution, parameters a and b */
FloatLo RandWeibull(FloatLo a, FloatLo b)
{
	FloatLo tmp;

	tmp=0.0-log(fabs(Rand1()));
	tmp=pow(tmp,1.0/b);
	return tmp*a;
}

/* approximates Gaussian distribution, mean 0.00, std. deviation 1.0 */
FloatLo Rand1Distrib(void)
{
	/* note 1.1547=4*sqrt(1/12), where 1/12 is the theoretical variance of
	   the Rand1() random variable, so this guarantees a theoretical
	   standard deviation of 1.0 */
	return (Rand1()+Rand1()+Rand1()+Rand1())/1.1547;
}

void Cross(vec res,vec a, vec b)
{
	vec ans;

	ans[0]=a[1]*b[2]-a[2]*b[1];
	ans[1]=a[2]*b[0]-a[0]*b[2];
	ans[2]=a[0]*b[1]-a[1]*b[0];
	VecScale(res,ans,1.0);
}

FloatLo Dot(vec a,vec b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

FloatLo getMag(vec a)
{
	return sqrt(Dot(a,a));
}

void Normalize(vec res,vec a)
{
	float b;

	b=getMag(a);
	res[0]=a[0]/b;
	res[1]=a[1]/b;
	res[2]=a[2]/b;
}

void VecAdd(vec res, vec a,vec b)
{
	res[0]=a[0]+b[0];
	res[1]=a[1]+b[1];
	res[2]=a[2]+b[2];
}

void VecSub(vec res, vec a,vec b)
{
	res[0]=a[0]-b[0];
	res[1]=a[1]-b[1];
	res[2]=a[2]-b[2];
}

void VecScale(vec res, vec a,FloatLo scale)
{
	res[0]=a[0]*scale;
	res[1]=a[1]*scale;
	res[2]=a[2]*scale;
}

void NegateVec(vec res, vec a)
{
	res[0]=-a[0];
	res[1]=-a[1];
	res[2]=-a[2];
}

void PrintVec(vec a)
{
	printf("X=%8.3f Y=%8.3f Z=%8.3f \n",(float) a[0], (float) a[1], (float) a[2]);
}

void Translate(vec a,vec trans)
{
	VecAdd(a,a,trans);
}

void Rotatecos(Int2 axis, FloatLo sina, FloatLo cosa, vec a,vec res)
{
/*	vec matrx[3];*/
	static vec ans;
	static Int2 nextaxis,prevaxis;

	prevaxis=(axis+2)%3;
	nextaxis=(axis+1)%3;
	ans[axis]=a[axis];
	ans[nextaxis]=a[nextaxis]*cosa-a[prevaxis]*sina;
	ans[prevaxis]=a[prevaxis]*cosa+a[nextaxis]*sina;
	/* use different variable just in case a and res are same */
	res[0]=ans[0];
	res[1]=ans[1];
	res[2]=ans[2];
/*	matrx[axis][axis]=1;
	matrx[(axis+1)%3][axis]=0;
	matrx[(axis+2)%3][axis]=0;
	matrx[axis][(axis+1)%3]=0;
	matrx[axis][(axis+2)%3]=0;
	matrx[(axis+1)%3][(axis+1)%3]=cosa;
	matrx[(axis+2)%3][(axis+2)%3]=cosa;
	matrx[(axis+2)%3][(axis+1)%3]=sina;
	matrx[(axis+1)%3][(axis+2)%3]=-sina;*/
	/* use different variable just in case a and res are same */
/*	ans[0]=Dot(a,matrx[0]);
	ans[1]=Dot(a,matrx[1]);
	ans[2]=Dot(a,matrx[2]);
	VecScale(res,ans,1.0);*/
}

void GetDihedral(vec v1,vec v2,vec v3,vec v4,FloatLo rng,FloatLo *dihed,FloatLo *ba1,FloatLo *ba2,FloatLo *bl1,FloatLo *bl2,FloatLo *bl3)
{
        vec v12,v23,v34,MainAxis,u,v,U,V;
        FloatLo angle;   
                
        VecSub(v34,v4,v3);
        VecSub(v23,v3,v2);
        VecSub(v12,v2,v1);
        Normalize(MainAxis,v23);
        VecScale(U,MainAxis,Dot(v34,MainAxis));
        NegateVec(U,U);
        VecAdd(U,v34,U);
        Normalize(u,U);
        VecScale(V,MainAxis,Dot(v12,MainAxis));
        VecSub(V,v12,V); 
        NegateVec(V,V);  
        Normalize(v,V);   
        /* between zero and Pi by default */
        angle=Dot(u,v);
        if (angle>1.00)
                angle=1.00;
        if (angle<-1.00)
                angle=-1.00;
        angle=RADTODEG*acos(angle);
        Cross(U,u,v);    
        if (Dot(U,MainAxis)>0)
                angle=-angle;
        while (angle<rng)
                angle+=360.0;
        while (angle>rng+360)
                angle-=360.0;
        *dihed=angle;   
        *bl1=getMag(v12);   
        *bl2=getMag(v23);
        *bl3=getMag(v34);
        Normalize(v12,v12);   
        Normalize(v23,v23);  
        Normalize(v34,v34);
		angle=Dot(v12,v23);
		if (angle>1.00)
                angle=1.00;
        if (angle<-1.00)
                angle=-1.00;
        angle=180.0-RADTODEG*acos(angle);
        while (angle<rng)
                angle+=360.0;
        while (angle>rng+360)
                angle-=360.0;
        *ba1=angle;
		angle=Dot(v23,v34);
		if (angle>1.00)
                angle=1.00;
        if (angle<-1.00)
                angle=-1.00;
        angle=180.0-RADTODEG*acos(angle);
        while (angle<rng)     
                angle+=360.0;
        while (angle>=rng+360.0)
                angle-=360.0;
        *ba2=angle;
}

void VechitoVec(vec res,vechi a)
{
	res[0]=(FloatLo)a[0];
	res[1]=(FloatLo)a[1];
	res[2]=(FloatLo)a[2];
}

void VectoVechi(vechi res,vec a)
{
	res[0]=(FloatHi)a[0];
	res[1]=(FloatHi)a[1];
	res[2]=(FloatHi)a[2];
}

void Crosshi(vechi res,vechi a, vechi b)
{
	vechi ans;

	ans[0]=a[1]*b[2]-a[2]*b[1];
	ans[1]=a[2]*b[0]-a[0]*b[2];
	ans[2]=a[0]*b[1]-a[1]*b[0];
	VecScalehi(res,ans,1.0);
}

FloatHi Dothi(vechi a,vechi b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

FloatHi getMaghi(vechi a)
{
	return sqrt(Dothi(a,a));
}

void Normalizehi(vechi res,vechi a)
{
	FloatHi b;

	b=getMaghi(a);
	res[0]=a[0]/b;
	res[1]=a[1]/b;
	res[2]=a[2]/b;
}

void VecAddhi(vechi res, vechi a,vechi b)
{
	res[0]=a[0]+b[0];
	res[1]=a[1]+b[1];
	res[2]=a[2]+b[2];
}

void VecSubhi(vechi res, vechi a,vechi b)
{
	res[0]=a[0]-b[0];
	res[1]=a[1]-b[1];
	res[2]=a[2]-b[2];
}

void VecScalehi(vechi res, vechi a,FloatHi scale)
{
	res[0]=a[0]*scale;
	res[1]=a[1]*scale;
	res[2]=a[2]*scale;
}

void NegateVechi(vechi res, vechi a)
{
	res[0]=-a[0];
	res[1]=-a[1];
	res[2]=-a[2];
}

void PrintVechi(vechi a)
{
	printf("X=%8.3f Y=%8.3f Z=%8.3f \n",(double) a[0], (double) a[1], (double) a[2]);
}

void Translatehi(vechi a,vechi trans)
{
	VecAddhi(a,a,trans);
}

void GetDihedralhi(vechi v1,vechi v2,vechi v3,vechi v4,FloatHi rng,FloatHi *dihed,FloatHi *ba1,FloatHi *ba2,FloatHi *bl1,FloatHi *bl2,FloatHi *bl3)
{
	vec v1lo,v2lo,v3lo,v4lo;
	FloatLo rnglo,dihedlo,ba1lo,ba2lo,bl1lo,bl2lo,bl3lo;
	
	VechitoVec(v1lo,v1);
	VechitoVec(v2lo,v2);
	VechitoVec(v3lo,v3);
	VechitoVec(v4lo,v4);
	rnglo=(FloatLo)rng;
	GetDihedral(v1lo,v2lo,v3lo,v4lo,rnglo,&dihedlo,&ba1lo,&ba2lo,&bl1lo,&bl2lo,&bl3lo);
	*dihed=(FloatHi)dihedlo;
	*ba1=(FloatHi)ba1lo;
	*ba2=(FloatHi)ba2lo;
	*bl1=(FloatHi)bl1lo;
	*bl2=(FloatHi)bl2lo;
	*bl3=(FloatHi)bl3lo;
}

Int4 iDot(ivec a,ivec b)
{
	return ((a[0]/0x100)*(b[0]/0x100)+(a[1]/0x100)*(b[1]/0x100)+(a[2]/0x100)*(b[2]/0x100));
}

Int4 igetMag(ivec a)
{
	return (Int4)(sqrt((double)((a[0]/0x100)*(a[0]/0x100)+(a[1]/0x100)*(a[1]/0x100)+(a[2]/0x100)*(a[2]/0x100))))*0x100;
}

void iNormalize(ivec res,ivec a)
{
	Int4 b;

	b=igetMag(a);
	res[0]=(a[0]*0x100/b)*0x100;
	res[1]=(a[1]*0x100/b)*0x100;
	res[2]=(a[2]*0x100/b)*0x100;
}

void iVecAdd(ivec res,ivec a,ivec b)
{
	res[0]=a[0]+b[0];
	res[1]=a[1]+b[1];
	res[2]=a[2]+b[2];
}

void iVecSub(ivec res,ivec a,ivec b)
{
	res[0]=a[0]-b[0];
	res[1]=a[1]-b[1];
	res[2]=a[2]-b[2];
}

void iVecScale(ivec res, ivec a,Int4 scale)
{
	res[0]=(a[0]/0x100)*(scale/0x100);
	res[1]=(a[1]/0x100)*(scale/0x100);
	res[2]=(a[2]/0x100)*(scale/0x100);
}

void iNegateVec(ivec res, ivec a)
{
	res[0]=-a[0];
	res[1]=-a[1];
	res[2]=-a[2];
}

void iPrintVec(ivec a)
{
	printf("X=%d.%d Y=%d.%d Z=%d.%d \n",a[0]/0x10000,a[0]&0xffff,a[1]/0x10000,a[1]&0xffff,a[2]/0x10000,a[2]&0xffff);
}

void iTranslate(ivec a,ivec trans)
{
	iVecAdd(a,a,trans);
}

void iRotatecos(Int2 axis, Int4 sina, Int4 cosa,ivec a,ivec res)
{
	ivec imatrx[3];
	ivec ians;
	
	imatrx[axis][axis]=1<<16;
	imatrx[(axis+1)%3][axis]=0;
	imatrx[(axis+2)%3][axis]=0;
	imatrx[axis][(axis+1)%3]=0;
	imatrx[axis][(axis+2)%3]=0;
	imatrx[(axis+1)%3][(axis+1)%3]=cosa;
	imatrx[(axis+2)%3][(axis+2)%3]=cosa;
	imatrx[(axis+2)%3][(axis+1)%3]=sina;
	imatrx[(axis+1)%3][(axis+2)%3]=-sina;
	/* use different variable just in case a and res are same */
	ians[0]=iDot(a,imatrx[0]);
	ians[1]=iDot(a,imatrx[1]);
	ians[2]=iDot(a,imatrx[2]);
	res[0]=ians[0];
	res[1]=ians[1];
	res[2]=ians[2];
}


/*  
$Log: vecmath.c,v $
Revision 1.6  2002/08/01 21:54:30  feldman
Fixed some potential rounding errors

Revision 1.5  2001/10/02 21:53:02  feldman
added

Revision 1.4  2001/03/29 02:52:24  feldman
fixed minor compiler warnings

Revision 1.3  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.2  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

