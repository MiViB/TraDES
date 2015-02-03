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
#include <geometry.h>

/* The function walks through atom information returns a pointer to trajectory angles, which
contains information on angle. Angles can not be calculated for the first residue.
Phi represents the cos of Phi and Psi, which stores theta is in degrees, in the case of a
CA walk, otherwise the normal phi and psi are stored here */ 
/* for a phi-psi walk, no values are returned for the first and last residues,
and "omega i" is actually omega (i+1) */
/* for the ca walk, all values are stored with the first element containing only
cos phi at residue 2 (the angle between the first 3 residues).  The final node
(for a total N-2 of them) contains cos phi and theta at residue N-1 (used to
place residue N).  It may be desirable to make up false values for residues 1 and
N for the sake of completeness but this is left up to the user */
PRS GetTrjAngle(PMMD pmmd, Int2 DoChis, Int2 Model, FloatLo tmp, DValNodePtr pdnIncList) {
  PDNMG pdnmg,pdnmgNext;
  PRS prshead=NULL, prshere=NULL, prslast=NULL;
  PMGD pmgd,pmgdNext;
  PVNMA pvnma;
  PMAD pmad;
  PALD pald;
  Int2 PNTR missCalpha;
  FloatLoPtr * flpp;
  Int2 j=0,i=0;
  vec v1={0.0,0.0,0.0},v2={0.0,0.0,0.0},v3={0.0,0.0,0.0},v4={0.0,0.0,0.0},vCANext,vCAHere,vTmp;
  vec vZero={0.0,0.0,0.0};
  FloatLo dummy[3];
  FloatLo dihed, ba1, ba2, bl1, bl2,bl3,ftmp;
  ValNodePtr vnpRama;
 
  if (pmmd == NULL){
    ErrPostEx(SEV_ERROR,3,1, "Bad molecule");
    return NULL; 
  }
  if (WALKTYPE==WALK_PHIPSI) {
    vnpRama=ConvertNode((PFB)pmmd,AM_MGD);
    if (vnpRama!=NULL)
      prshead=Rama(vnpRama,Model);
    else return NULL;
    if (prshead==NULL)
      return NULL;
    /* first prs has phi2, psi2, and so on */
    prshere=prshead;
    while (prshere) {
/*printf("%f %f\n",prshere->Phi,prshere->Psi);*/
      if (prshere->Mag!=VL_MISSCALPHA)
        prshere->Mag=VL_TRANSRESIDUE;
      pmgd=(PMGD)(prshere->pfbThis);
      pdnmg=pmgd->pdnmgLink;
      pdnmgNext=pdnmg->next;
      if (pdnmgNext!=NULL) {
        pmgdNext=(PMGD)(pdnmgNext->data.ptrvalue);
        if (GetCoOrds(pmgdNext," CA ",vZero,vCANext,Model) && GetCoOrds(pmgd," CA ",vZero,vCAHere,Model)) {
 	       VecSub(vTmp,vCANext,vCAHere);
 	       if ((getMag(vTmp)<3.2) && (prshere->Mag!=VL_MISSCALPHA))
		  prshere->Mag=VL_CISRESIDUE;
	}
      }
      if (DoChis>0) {
		  if (MeasureResidueChis(pmgd,&prshere,Model,DoChis==1?FALSE:TRUE,tmp,pdnIncList)!=ERR_SUCCESS)
	      	ErrPostEx(SEV_WARNING,1,2,"Unable to get chi angles for residue %ld",(long int)(pdnmg->choice));
	    }
      prshere=prshere->next;
    }
    vnpRama=ValNodeFree(vnpRama);
    return prshead;
  }
  /* else WALKTYPE==WALK_CA */
  flpp = (FloatLoPtr *) MemNew(sizeof(FloatLoPtr)*(pmmd -> iResCount));
  missCalpha = (Int2 *) MemNew(sizeof(Int2)*(pmmd -> iResCount)); 

  pdnmg = pmmd -> pdnmgHead;
  while (pdnmg) {
    pmgd = (PMGD) pdnmg -> data.ptrvalue;
    pvnma = pmgd -> pvnmaAHead;

    /* deal with gap in sequence are represented in missCalpha */
    if ((pmad = FindCAlpha(pvnma))==NULL){
      flpp[i]=dummy;
      missCalpha[i] = 1;
      i++;
    }
    else {
      pald = GetAtomLocs (pmad, Model);    
      if (pald != NULL){
	flpp[i] = pald -> pflvData;
	missCalpha[i] = 0;
	i++;
      }
      else{ 
	flpp[i]=dummy;
	missCalpha[i] = 1;
	i++;
      }
    }
    pdnmg = pdnmg->next;   
  }
  pdnmg = ((pmmd->pdnmgHead)->next);
  
  for (j=1; j< ((pmmd -> iResCount)-1);j++){ 
    if (j>1){
      for (i=0; i<3; i++)
	v1[i] = flpp [j-2][i];} 
    for (i=0; i<3; i++)
      v2[i] = flpp [j-1][i]; 
    for (i=0;  i<3; i++)
      v3[i] = flpp [j-0][i]; 
    for (i=0; i<3; i++)
      v4[i] = flpp [j+1][i];  
    prshere =(PRS)MemNew(sizeof(RS));
  
    if (prshead==NULL)
      prshead = prshere;
    else
      prslast -> next = prshere;
    
    GetDihedral(v1,v2,v3,v4,0,&dihed,&ba1,&ba2,&bl1,&bl2,&bl3);
    
    prshere -> pfbThis = (PFB) (pdnmg -> data.ptrvalue);
        
    if (((j==1) && (missCalpha[j-1]==1 || missCalpha[j]==1 || missCalpha[j+1]==1)) || ((j>1) && ((missCalpha[j-2]==1 || missCalpha[j-1]==1 || missCalpha[j]==1 || missCalpha[j+1]==1)))) {
      prshere -> Mag = VL_MISSCALPHA;
    }
    else{
		ftmp=floor(0.5 +(-cos(ba2*DEGTORAD)*(FloatLo)(TRAJDIV/2.0))+(FloatLo) TRAJDIV/2.0 )-(FloatLo) (TRAJDIV/2.0);
		if (ftmp>=TRAJDIV)
			ftmp-=TRAJDIV;
      prshere -> Phi = ftmp/((FloatLo)(TRAJDIV/2.0)); 
	  /* cheat here to ensure we never reach the poles which are invalid and singular points */
	  if (prshere->Phi>0.99)
		  prshere->Phi=0.99;
	  if (prshere->Phi<-0.99)
		  prshere->Phi=-0.99;
      if (j>1) {
	prshere -> Psi = dihed+180.0;
	while (prshere->Psi>=360.0)
	  prshere->Psi -= 360.0;
	ftmp=floor(0.5+prshere -> Psi * (FloatLo)(TRAJDIV/360.0));
	if (ftmp>=TRAJDIV)
		ftmp-=TRAJDIV;
	prshere -> Psi = ftmp*360.0/(FloatLo)TRAJDIV;
      }
      else 
	prshere -> Psi = 0.0;	
      if (bl3<3.2) 
	prshere ->Mag = VL_CISRESIDUE;
      else 
	prshere ->Mag = VL_TRANSRESIDUE; 
    }
    if (DoChis>0) {
      if (MeasureResidueChis((PMGD)(pdnmg->data.ptrvalue),&prshere,Model,DoChis==1?FALSE:TRUE,tmp,pdnIncList)!=ERR_SUCCESS)
      	ErrPostEx(SEV_WARNING,1,2,"Unable to get chi angles for residue %ld",(long int)(pdnmg->choice));
    }
    prslast = prshere;
/*printf("%d ph: %f ps: %f\n",j+1,prshere->Phi,prshere->Psi); */
    pdnmg = pdnmg -> next;    
  }  
  
  missCalpha = MemFree(missCalpha);
  prshere->next=NULL;
  flpp = MemFree(flpp);
  return prshead;
}

/* Calculates the area of peak value of the Gaussian */
Int4 CalcMagnitude(FloatLo sigma_x, FloatLo sigma_y, FloatLo yDeg, Int2 percent,PTGS ptgs,Int2 IsCis){
  Int4 peakG;
  Int4 integral;
  FloatLo yInc;
  
  if (IsCis == 1)
    integral = ptgs -> CisTrajIntegral;
  else
    integral = ptgs -> TrajIntegral;

  if (sigma_y == 0 && sigma_x == 0)    
    peakG =(Int4) (((FloatLo)integral/(FloatLo)percent)*(100.0-(FloatLo)percent));
  else
    if (sigma_y == 0)
      peakG =(Int4) ((((FloatLo)integral/(FloatLo)percent)*(100.0-(FloatLo)percent))/ (((FloatLo)(ptgs -> dim))*(sqrt(2.0*(PI)) *sigma_x/360.0)));
    else 
      if (sigma_x == 0){
	yInc = 1.0/fabs(yDeg - RADTODEG * acos(cos(yDeg*DEGTORAD) -2.0/(FloatLo)(ptgs -> dim))); 
	peakG =(Int4) ((((FloatLo)integral/(FloatLo)percent)*(100.0-(FloatLo)percent))/((yInc * sqrt(2.0*(PI)) * sigma_y)));
      }
	
      else{
	/*  */
	yInc = 1.0/fabs(yDeg - RADTODEG * acos(cos(yDeg*DEGTORAD)-2.0/(FloatLo)(ptgs -> dim))); 
	peakG =(Int4) ((((FloatLo)integral/(FloatLo)percent)*(100.0-(FloatLo)percent))/ (2.0 * PI * yInc *((FloatLo)(ptgs ->dim))*sigma_x*sigma_y/360.0));
      }
  return peakG;
}


/* Calculates a Gaussian distribution given: deviation in x and y, and x y cordinates which are all in degrees, magnitude which will be an interger. */    
Int4 CalcTGSElement(FloatLo sigma_x,FloatLo sigma_y,FloatLo xCord,FloatLo yCord,FloatLo x,FloatLo y,Int4 mag){
  Int4 tmpvalue=0,tmpvalue2=0,value=0;
  FloatLo i,j;

  /* special cases */
  if (sigma_y == 0 && sigma_x == 0){ 
    if (fabs(x-xCord)<PRECISION && fabs(y-yCord)<PRECISION)
      value = mag;   
    else 
      value = 0;
  }
  else if (sigma_x == 0){
    if (fabs(x - xCord)<PRECISION) {
      value = (Int4)floor((FloatLo)mag*exp(-( ((y-yCord)*(y-yCord))/(2*sigma_y*sigma_y) )));
	    /* wrap around y-axis in Phi-Psi space only */
 	    if (WALKTYPE==WALK_PHIPSI) {
	      j = 360.0;
	      do {
		      tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((y+j-yCord)*(y+j-yCord))/(2*sigma_y*sigma_y) )));
					value += tmpvalue2;
					j += 360.0;
    	  } while (tmpvalue2 >(mag/1000));
	
  	    j = 360.0;
    	  do {
		      tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((y-j-yCord)*(y-j-yCord))/(2*sigma_y*sigma_y) )));
					value += tmpvalue2;
					j += 360.0;
    	  } while (tmpvalue2 >(mag/1000));	
			}
    }
    else
      value = 0;
  }
  else if (sigma_y == 0){
    if (fabs(y - yCord)<PRECISION){
      i = 360.0;
      value = (Int4)floor((FloatLo)mag*exp(-( ((x-xCord)*(x-xCord))/(2*sigma_x*sigma_x) )));
      do {
				tmpvalue = (Int4)floor((FloatLo)mag*exp(-( ((x+i-xCord)*(x+i-xCord))/(2*sigma_x*sigma_x) )));
				value += tmpvalue;
				i += 360.0;
      } while (tmpvalue >(mag/1000));
      
      i = 360.0;
      do {
				tmpvalue = (Int4)floor((FloatLo)mag*exp(-( ((x-i-xCord)*(x-i-xCord))/(2*sigma_x*sigma_x) )));
				value += tmpvalue;
				i += 360.0;
      } while (tmpvalue >(mag/1000));	
    }
    else
      value=0;
  }
  else{ 
    value = (Int4)floor((FloatLo)mag*exp(-( ((x-xCord)*(x-xCord)/(2*sigma_x*sigma_x))+ ((y-yCord)*(y-yCord)/(2*sigma_y*sigma_y)) )));
    i = 360.0;	    
    
    do {
      tmpvalue = (Int4)floor((FloatLo)mag*exp(-( ((x+i-xCord)*(x+i-xCord)/(2*sigma_x*sigma_x))+ ((y-yCord)*(y-yCord)/(2*sigma_y*sigma_y)) )));
      value += tmpvalue;
	    /* wrap around y-axis in Phi-Psi space only */
	    if (WALKTYPE==WALK_PHIPSI) {
		    j = 360.0;	
		    do {
	  	    tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((x+i-xCord)*(x+i-xCord)/(2*sigma_x*sigma_x))+ ((y+j-yCord)*(y+j-yCord)/(2*sigma_y*sigma_y)) )));
	    	  value += tmpvalue2;
	      	j += 360.0;
		    } while (tmpvalue2 >(mag/1000));
		
	  	  j = 360.0;
	    	do {
	      	tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((x+i-xCord)*(x+i-xCord)/(2*sigma_x*sigma_x))+( (y-j-yCord)*(y-j-yCord)/(2*sigma_y*sigma_y)) )));
		      value += tmpvalue2;
		      j += 360.0;
  		  } while (tmpvalue2 >(mag/1000));	
	    }
      i += 360.0;
    } while (tmpvalue >(mag/1000));
    
    i = 360.0;
    do {
      tmpvalue = (Int4)floor((FloatLo)mag*exp(-( ((x-i-xCord)*(x-i-xCord)/(2*sigma_x*sigma_x))+( (y-yCord)*(y-yCord)/(2*sigma_y*sigma_y)) )));
      value += tmpvalue;
	    /* wrap around y-axis in Phi-Psi space only */
	    if (WALKTYPE==WALK_PHIPSI) {
		    j = 360.0;	
		    do {
	  	    tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((x-i-xCord)*(x-i-xCord)/(2*sigma_x*sigma_x))+ ((y+j-yCord)*(y+j-yCord)/(2*sigma_y*sigma_y)) )));
	    	  value += tmpvalue2;
	      	j += 360.0;
		    } while (tmpvalue2 >(mag/1000));
		
	  	  j = 360.0;
	    	do {
	      	tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((x-i-xCord)*(x-i-xCord)/(2*sigma_x*sigma_x))+( (y-j-yCord)*(y-j-yCord)/(2*sigma_y*sigma_y)) )));
		      value += tmpvalue2;
		      j += 360.0;
  		  } while (tmpvalue2 >(mag/1000));	
	    }
      i += 360.0;
    } while (tmpvalue >(mag/1000));	
    if (WALKTYPE==WALK_PHIPSI) {
      j = 360.0;
	    do {
  	    tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((x-xCord)*(x-xCord)/(2*sigma_x*sigma_x))+ ((y+j-yCord)*(y+j-yCord)/(2*sigma_y*sigma_y)) )));
				value += tmpvalue2;
				j += 360.0;
   	  } while (tmpvalue2 >(mag/1000));
	    j = 360.0;
   	  do {
  	    tmpvalue2 = (Int4)floor((FloatLo)mag*exp(-( ((x-xCord)*(x-xCord)/(2*sigma_x*sigma_x))+ ((y-j-yCord)*(y-j-yCord)/(2*sigma_y*sigma_y)) )));
				value += tmpvalue2;
				j += 360.0;
   	  } while (tmpvalue2 >(mag/1000));	
		}
  }
  return value; 
}


/* Calculates the following values for the trajectory graph: residue#,amino acid letter, graph dimension, peak value, rotomer info, presents of a sulfide bond, standard deviations, integral,trajectory graph data, first non-zero rows, number of non-zero rows, number of elements (0%,5%,10%,15%), time out and cis probability */
/* NOTE: this function calls FillTG so make sure databases are closed! */
TrajErr FillInTGS(PTGS ptgsHere, PDNMG pdnmg, PRS prsHere, PRS prsLast, Int2 start, FloatLo sigma_x, FloatLo sigma_y, Int4 mag, Int4 *Sarray,Int2 percent){
  Int4 csum,cnt,xCord,yCord,c2;
  Char cres[2];
  Int4 height;
  FloatLo xDeg, yDeg,x,y; 
  Char lastRes = ' ';
  PMGD pmgd;

  height = mag;
  if (ptgsHere == NULL)
    return ERR_FAIL;
  if (pdnmg == NULL)
    return ERR_FAIL;
  pmgd = (PMGD) (pdnmg -> data.ptrvalue);
  if (pdnmg ->next != NULL)
    lastRes = ((PMGD)((pdnmg ->next)->data.ptrvalue))->pcIUPAC[0];
  ptgsHere -> resnum = (pdnmg -> choice)-start+1;
/*printf("%d....%d.. ",pdnmg->choice,ptgsHere->resnum);  */
  ptgsHere -> AA = pmgd -> pcIUPAC[0];

  ptgsHere -> dim = TRAJDIV;  
  ptgsHere -> rotid = ComputeRotid(prsHere);
  ptgsHere -> markovsf = MARKOV_SCALE_FACTOR;  
  ptgsHere -> pSS = 0.0;
  if (Sarray[(ptgsHere -> resnum)-1]) 
    ptgsHere -> pSS = 1.0;
  /* default W -- CA walk does not record omega angles, just cis or not */
  ptgsHere -> ChiWMean=CHI_W;
  ptgsHere -> ChiWSD=CHISD_W;
  sprintf(cres,"%c", pmgd ->pcIUPAC[0]); 
  
  /* all but first/last residue */
  if (prsHere!=NULL){
    if (WALKTYPE==WALK_CA) {
      xCord = (Int4)floor((((prsHere -> Psi))/360.0)*TRAJDIV);
      yCord = (Int4)floor((((prsHere -> Phi)+1.0)/2.0)*TRAJDIV);
    }
    else { /* WALKTYPE==WALK_PHIPSI */
      xCord = (Int4)floor((((prsHere -> Phi)+180.0)/360.0)*TRAJDIV);
      yCord = (Int4)floor((((prsHere -> Psi)+180.0)/360.0)*TRAJDIV);
/*printf("%ld %ld\n",xCord,yCord);*/
      ptgsHere->ChiWMean=prsHere->Omega;
      /* take geometric mean of deviations here */
      ptgsHere->ChiWSD=2.0*sqrt(sigma_x*sigma_y);
      if (ptgsHere->ChiWSD>CHISD_W)
        ptgsHere->ChiWSD=CHISD_W;
    }
    if (WALKTYPE==WALK_CA) {
      xDeg = prsHere -> Psi;
      yDeg = acos(prsHere -> Phi)*RADTODEG;
    }
    else {
      xDeg = prsHere->Phi+180.0;
      yDeg = prsHere->Psi+180.0;
    }
    if ((prsHere ->Mag)==VL_CISRESIDUE && lastRes=='P'){
      ptgsHere -> pCis = 1.0;
    }
    else 
      if ((prsHere -> Mag) == VL_MISSCALPHA){
	if (lastRes == 'P')
	  ptgsHere -> pCis = P_CISP;
	else 
	  ptgsHere -> pCis = 0;
      }
    else
      ptgsHere -> pCis =0.0;
    
    /* Call Howie's function, when missing a CAlpha. Assumes no cis trj graph */   
    if (prsHere ->Mag == VL_MISSCALPHA){
      ptgsHere -> TrajIntegral = FillTG(TRAJ_STANDARD, cres, 0, NULL,TRAJDIV,ptgsHere->TrajGraph, FALSE); 
      if (ptgsHere -> AA == 'P')	
	ptgsHere -> CisTrajIntegral = FillTG(TRAJ_STANDARD, cres, 0, NULL,TRAJDIV,ptgsHere->CisTrajGraph, TRUE)  ;
      else
	ptgsHere ->CisTrajIntegral = 0;
    }

    else if ((prsLast != NULL && (prsLast ->Mag)==VL_CISRESIDUE)) {
      if (percent) {
	ptgsHere -> CisTrajIntegral = FillTG(TRAJ_STANDARD, cres, 0, NULL,TRAJDIV,ptgsHere->CisTrajGraph,TRUE);  
	height =  CalcMagnitude(sigma_x, sigma_y,yDeg,percent, ptgsHere,1); 
      }
      /* Calculation to noise */
      for (c2=0;c2<(TRAJDIV*TRAJDIV);c2++){
	x = ((FloatLo)(c2%TRAJDIV)/(FloatLo)TRAJDIV)*360.0;
	if (WALKTYPE==WALK_CA)
	  y = (acos(((FloatLo)(c2/TRAJDIV)/(FloatLo)TRAJDIV)*2.0-1.0)*RADTODEG);
	else
	  y = ((FloatLo)(c2/TRAJDIV)/(FloatLo)TRAJDIV)*360.0;
	/* Force x to equal xDeg to avoid rounding errors */
	if (sigma_x ==0 && c2%TRAJDIV == xCord)
	  x = xDeg;
	if (sigma_y ==0 && c2/TRAJDIV == yCord)
	  y = yDeg;
      if (percent) 
	ptgsHere -> CisTrajGraph[c2] += CalcTGSElement(sigma_x,sigma_y,xDeg,yDeg,x,y,height);
      else
	ptgsHere -> CisTrajGraph[c2] = CalcTGSElement(sigma_x,sigma_y,xDeg,yDeg,x,y,height);
      }
      csum =0;
      for (cnt=0;cnt<(ptgsHere->dim)*(ptgsHere -> dim);cnt++)
	csum += (ptgsHere ->CisTrajGraph[cnt]);
      ptgsHere -> CisTrajIntegral = csum;
      ptgsHere->TrajIntegral=1;
      ptgsHere->TrajGraph[0]=1;
    }
    else{
     
      /* trans peptide */
      if (percent) {
	ptgsHere -> TrajIntegral = FillTG(TRAJ_STANDARD, cres, 0, NULL,TRAJDIV,ptgsHere->TrajGraph,FALSE); 
     
	/* Calculates the area of peak value of the Gaussian,if percent specified on command line */
	height =  CalcMagnitude(sigma_x,sigma_y,yDeg,percent, ptgsHere,0);
      }
      for (c2=0;c2<(TRAJDIV*TRAJDIV);c2++){
	x = ((FloatLo)(c2%TRAJDIV)/(FloatLo)TRAJDIV)*360.0;
	if (WALKTYPE==WALK_CA)
	  y = (acos(((FloatLo)(c2/TRAJDIV)/(FloatLo)TRAJDIV)*2.0-1.0)*RADTODEG);
	else
	y = ((FloatLo)(c2/TRAJDIV)/(FloatLo)TRAJDIV)*360.0;
	/* Force x to equal xDeg to avoid rounding errors */
	if (sigma_x == 0) {
		if (c2%TRAJDIV == xCord)
			x = xDeg;
		else {
			/* x != xDeg (arbitrarily add 10) */
			x = xDeg + 10.0;
			if (x>360.0)
				x-=360.0;
		}
	}
	if (sigma_y == 0) {
		if (c2/TRAJDIV == yCord)
			y = yDeg;
		else {
			y = yDeg + 10.0;
			if (y>360.0)
				y-=360.0;
		}
	}
	if (percent) 	
	  ptgsHere -> TrajGraph[c2] += CalcTGSElement(sigma_x,sigma_y,xDeg,yDeg,x,y,height);	
	else 
	  ptgsHere -> TrajGraph[c2] = CalcTGSElement(sigma_x,sigma_y,xDeg,yDeg,x,y,height);
      }    
      
      csum =0;
      for (cnt=0;cnt<(ptgsHere ->dim)*(ptgsHere -> dim);cnt++)
	csum += (ptgsHere ->TrajGraph[cnt]);
      ptgsHere ->TrajIntegral = csum;     
      ptgsHere -> CisTrajIntegral = 0;
    }  
  }   
  else {
      ptgsHere->ChiWMean=CHI_W;
      ptgsHere->ChiWSD=CHISD_W;
	if (lastRes == 'P')
	  ptgsHere -> pCis = P_CISP;
	else 
	  ptgsHere -> pCis = 0.0;
      ptgsHere->TrajIntegral=1;
      ptgsHere->TrajGraph[0]=1;
      if (ptgsHere->AA == 'P') {	
	ptgsHere->CisTrajIntegral=1;
        ptgsHere->CisTrajGraph[0]=1;
      }
      else
	ptgsHere->CisTrajIntegral=0;
  }
  if (prsLast == NULL || (Int2) (prsLast -> Mag == VL_MISSCALPHA))
    TrajCalcSparsity(ptgsHere,0,NULL);
  else
    TrajCalcSparsity(ptgsHere,(Int2)(prsLast -> Mag),NULL);
  TrajCalcTout(ptgsHere);
  return ERR_SUCCESS;
}  

TrajErr AddDisulfideRestraint(Int4 a,Int4 b)
{
	PNN pnnhere;

	pnnhere = (PNN) MemNew(sizeof(NN));
	if (b<a){	  
		pnnhere -> res1 = b;
		pnnhere -> res2 = a;
	}
	else {	  
		pnnhere -> res1 = a;
		pnnhere -> res2 = b;
	}
	  
	StringCpy(pnnhere -> AtomName1," SG ");  
	StringCpy(pnnhere -> AtomName2," SG ");  
	pnnhere -> MeanDist = BL_C_SGSG;
	pnnhere -> MinDelta = BLSD_C_SGSG;
	pnnhere -> MaxDelta = BLSD_C_SGSG;
	/* set these to INFINITY - i.e. ignore */
	pnnhere -> Angle1 = CONSTR_INFINITY;
	pnnhere -> Angle2 = CONSTR_INFINITY;
	pnnhere -> Dihedral01 = CONSTR_INFINITY;
	pnnhere -> Dihedral12 = CONSTR_INFINITY;
	pnnhere -> Dihedral23 = CONSTR_INFINITY;
	pnnhere -> prob = 1.0;
	if (AddDistConstraint(pnnhere) != ERR_SUCCESS) {
		ErrPostEx(SEV_ERROR,8,1,"Attempted to add invalid disulfide bridge");
		return ERR_FAIL;
	}
	return ERR_SUCCESS;
}

/* Determines disulfide bond */
TrajErr IsBondSulfide(PMGD pmgd,Int4 *Sarray,Int2 AddRestraint){
    /*PVNMA pvnmaHead;*/
    PVNMB pvnmb;
    PMMD pmmd;
    PMAD pmadAfter, pmadBefore;
    Int4 numres,a,b,i;
    Char atomb,atoma;
    PMBD pmbd;
    /*PDNMG pdnmg;*/
    PMGD pmgdAfter, pmgdBefore;     

    pmmd = (PMMD) pmgd -> pfbParent;
    pvnmb = pmmd -> pvnmbIRBHead;
    numres = (pmmd ->iResCount);
   
    
    /* Fills in array with zeros */
    for (i =0;i<numres; i++)
      Sarray[i]= 0;
    
    while (pvnmb != NULL ){
      pmbd = (PMBD) pvnmb -> data.ptrvalue;
      pmadBefore = (pmbd->pmadFrom);
      pmadAfter = (pmbd ->pmadTo);
      atoma = pmadBefore -> pcAName[1];
      atomb = pmadAfter -> pcAName[1];
      
      if (atoma=='S' &&  atomb=='S'){
	pmgdAfter = (PMGD)(pmadAfter -> pfbParent); 
	a = (pmgdAfter -> pdnmgLink) -> choice;
	pmgdBefore = (PMGD)(pmadBefore -> pfbParent); 
	b = (pmgdBefore -> pdnmgLink) -> choice;       	
	
	
	Sarray[a-1]= b;
	Sarray[b-1]= a;
	if (AddRestraint){
		if (AddDisulfideRestraint(a,b)!=ERR_SUCCESS)
			return ERR_FAIL;
	}      
      }
      pvnmb = pvnmb -> next;
	}
	return ERR_SUCCESS;
}

/* If temperature is specified then calculate new sigma values
	x = Zhang-DeLisi potential at this residue
	tmp = temperature in Kelvin
	delta_tmp = timestep in fs */
FloatLo CalcSigma(PDNMG pdnmg,FloatLo tmp, FloatLo delta_tmp,FloatHi x){
  PMGD pmgd; 
  /*Int2 i;*/
  FloatLo v,molW, sigma;
 
  pmgd = (PMGD) (pdnmg -> data.ptrvalue);
  
  molW = MWaalist[GetResIdxFromMGD(aalist,pmgd)];
  
  v = sqrt((tmp/ molW)*(8314.51+5293.18*atan(74.896*((FloatLo)0.0/*x*/)/tmp+1.5708))); /* NULL out Z-D potential for now - seems to work better? */
  if (WALKTYPE==WALK_CA)
        sigma = delta_tmp * v / 9404 /* 8334 was old, wrong number */;         /* Is it for both sigma? */
  else
	sigma=delta_tmp*v/4534.5;
/*  printf ("v = %f \n",v);*/
  return sigma; 
}

/* To direct the construction of a trajectory graph */
TrajErr BuildTrjGraph(PMMD pmmd, Int2 model, Int2 start, Int2 end, Int4 mag, FloatLo sigma_x, FloatLo sigma_y,CharPtr fnam,Int2 percent, FloatLo tmp, FloatLo  delta_tmp, DValNodePtr pdnIncList, Int2 SaveChis, CharPtr ssmask){
  PDNMG pdnmg;
  PMGD pmgd;
  PTGS ptgsHere=NULL;
  PRS prsHere, prsHead, prsLast;
  Int4 *Sarray;
  FloatHi potential;
  Int2 cnt,a;
  Boolean ulr;
  DValNodePtr pdnSubList;

  pdnmg = (pmmd -> pdnmgHead); 
  pmgd = (PMGD) (pdnmg -> data.ptrvalue);
  Sarray = (Int4 *) MemNew(sizeof(Int4)*(pmmd -> iResCount));
  
  if (IsBondSulfide(pmgd,Sarray,0)==ERR_FAIL)
		return ERR_FAIL;
  if (tmp == 0) {
  	/* add disulphides now */
  	for (cnt=start;cnt<=end;cnt++) {
  		a=Sarray[cnt-1];
  		if (a<cnt && a>0 && a>=start) {
  			if (AddDisulfideRestraint(a-start+1,cnt-start+1)!=ERR_SUCCESS)
  				return ERR_FAIL;
  		}
  	}
  }
  /* Fills in TGS for first residue */
  if (start == 1) {
	pdnSubList = (DValNodePtr)(pdnIncList->data.ptrvalue);
    potential = ((AdjListNodePtr)(pdnSubList->data.ptrvalue))->potential;
    /*printf ("potential %f \n",potential);*/
    if (((pmgd->pcIUPAC)[0]) == 'P')
      ptgsHere = NewTraj(1, TRAJDIV);
    else
      ptgsHere = NewTraj(0, TRAJDIV);

    
    if (tmp != 0){
      sigma_x = CalcSigma(pdnmg,tmp,delta_tmp,potential);
      sigma_y = sigma_x;
    }
	if (ssmask!=NULL) {
		if (ssmask[0]=='H' && WALKTYPE==WALK_PHIPSI) {
			sigma_x=0.0;
			sigma_y=0.0;
		}
	}

    /*********** change sigma's to array form ******************/
    if (FillInTGS(ptgsHere,pdnmg,NULL,NULL,start,sigma_x,sigma_y,mag,Sarray,percent)!=ERR_SUCCESS){
      ptgsHere = FreeTraj(ptgsHere);
      Sarray = MemFree(Sarray);
      ErrPostEx(SEV_ERROR,4,1,"Unable to create trajectory distribution");
      return ERR_FAIL;
    }     
    TrajGraphIntegrate(ptgsHere);
    TrajCalcNZ(ptgsHere);
	ulr=USE_LOTS_RAM;
	USE_LOTS_RAM=FALSE;
    if (TGInit(fnam,DB_CREATE,NULL)!=ERR_SUCCESS) {
  	PurgeGlobs();
  	ErrPostEx(SEV_FATAL,11,1,"Unable to open database, cannot continue");
  	return ERR_FAIL;
    }
    TrajGraphWrite(ptgsHere,USE_RLE,FALSE);
    TGClose();
	USE_LOTS_RAM=ulr;
    ptgsHere = FreeTraj(ptgsHere);
  }
  pdnmg = pdnmg ->next;
  /* pdnIncList has not been moved to ->next yet, still 'clean' */
  prsHead = GetTrjAngle(pmmd, SaveChis, model, tmp, pdnIncList);
  if (prsHead == NULL){
    prsHead = freeRS (prsHead);
    ErrPostEx(SEV_ERROR,4,1,"Unable to retrieve angles.");
    return ERR_FAIL;
  }
  prsHere= prsHead;
  prsLast = NULL;
  pdnIncList = pdnIncList->next; 
  while (pdnmg->next){
    
    if ((pdnmg ->choice >= start) && (pdnmg ->choice <= end)){ 
	  pdnSubList = (DValNodePtr)(pdnIncList->data.ptrvalue);
      potential = ((AdjListNodePtr)(pdnSubList->data.ptrvalue))->potential;
      /*printf ("res %d potential %f \n",pdnmg->choice,potential);*/
      pmgd = (PMGD) (pdnmg -> data.ptrvalue);
      if (((pmgd->pcIUPAC)[0]) == 'P')
        ptgsHere = NewTraj(1, TRAJDIV);
      else {
        ptgsHere = NewTraj(0, TRAJDIV);
		/* avoid cis problems */
		if ((prsLast != NULL && (prsLast ->Mag)==VL_CISRESIDUE))
			prsLast ->Mag=VL_TRANSRESIDUE;
      }
      
      if (traj_quiet == VERBOSITY_QUIET || traj_quiet == VERBOSITY_VERBOSE){
	printf("\b\b\b\b\b\b\b\b\b\b\b %4d/%-4d ",(pdnmg -> choice)-start+1,end-start+1); 
	fflush(stdout);
      } 
			if (traj_quiet==VERBOSITY_SILENT) {
				if (ProgramProgress<0) {
					/* abort signal from calling process */
					ptgsHere = FreeTraj (ptgsHere);
					prsHead = freeRS(prsHead);
				    CleanUpDB(tmpdbasename);
					return ERR_FAIL;					
				}
				ProgramProgress=(pdnmg->choice)-start+1;
				ProgramProgressMax=end-start+1;
			}

      if (tmp != 0){
	sigma_x = CalcSigma(pdnmg,tmp,delta_tmp,potential);
	sigma_y = sigma_x;
      }
	  if (ssmask!=NULL) {
		if (ssmask[pdnmg->choice-1]=='H') {
			if (WALKTYPE==WALK_PHIPSI || ssmask[pdnmg->choice-2]=='H') {
				sigma_x=0.0;
				sigma_y=0.0;
			}
		}
	  }
      

      if (FillInTGS(ptgsHere,pdnmg,prsHere,prsLast,start,sigma_x,sigma_y,mag,Sarray,percent)!=ERR_SUCCESS){
	ptgsHere = FreeTraj (ptgsHere);
	prsHead = freeRS(prsHead);
	ErrPostEx(SEV_ERROR,4,1,"Unable to create trajectory distribution");
	return ERR_FAIL;
      }
      TrajGraphIntegrate(ptgsHere);
      TrajCalcNZ(ptgsHere);
	  ulr=USE_LOTS_RAM;
	  USE_LOTS_RAM=FALSE;
      if (TGInit(fnam,DB_CREATE,NULL)!=ERR_SUCCESS) {
  	PurgeGlobs();
  	ErrPostEx(SEV_FATAL,11,1,"Unable to open database, cannot continue");
  	return ERR_FAIL;
      }
      TrajGraphWrite(ptgsHere,USE_RLE,FALSE);
      TGClose();
	  USE_LOTS_RAM=ulr;
      ptgsHere = FreeTraj(ptgsHere);
    }
    pdnmg = pdnmg -> next;
    prsLast = prsHere;
    prsHere = (prsHere -> next);    
    pdnIncList = pdnIncList->next; 
  }  
 
  /* For last residue */
  if (end == pmmd ->iResCount){

	pdnSubList = (DValNodePtr)(pdnIncList->data.ptrvalue);
    potential = ((AdjListNodePtr)(pdnIncList->data.ptrvalue))->potential;
    /*printf ("potential %f \n",potential);*/
   
 
    pmgd = (PMGD) (pdnmg -> data.ptrvalue);
    if (((pmgd->pcIUPAC)[0]) == 'P')
      ptgsHere = NewTraj(1,TRAJDIV);   
    else {
		ptgsHere = NewTraj(0,TRAJDIV);
		/* avoid cis problems */
		if ((prsLast != NULL && (prsLast ->Mag)==VL_CISRESIDUE))
			prsLast ->Mag=VL_TRANSRESIDUE;
	}
    
    if (tmp != 0){
      sigma_x = CalcSigma(pdnmg,tmp,delta_tmp,potential);
      sigma_y = sigma_x;
    }
	if (ssmask!=NULL) {
		if (ssmask[pdnmg->choice-1]=='H') {
			if (WALKTYPE==WALK_PHIPSI || ssmask[pdnmg->choice-2]=='H') {
				sigma_x=0.0;
				sigma_y=0.0;
			}
		}
	}

    if (FillInTGS(ptgsHere,pdnmg,NULL,prsLast,start,sigma_x,sigma_y,mag,Sarray,percent)!=ERR_SUCCESS){
      ptgsHere = FreeTraj (ptgsHere);
      prsHead = freeRS(prsHead);
      ErrPostEx(SEV_ERROR,4,1,"Unable to create trajectory distribution");
      return ERR_FAIL;
    }
    TrajGraphIntegrate(ptgsHere);
    TrajCalcNZ(ptgsHere);
    
    if (traj_quiet == VERBOSITY_QUIET || traj_quiet == VERBOSITY_VERBOSE){
      printf("\b\b\b\b\b\b\b\b\b\b\b %4d/%-4d \n",(pdnmg -> choice)-start+1,end-start+1); 
    }

	/* prevent TGInit from doing unnecessary TD loading */			
	ulr=USE_LOTS_RAM;
	USE_LOTS_RAM=FALSE;
    if (TGInit(fnam,DB_CREATE,NULL)!=ERR_SUCCESS) {
  	PurgeGlobs();
  	ErrPostEx(SEV_FATAL,11,1,"Unable to open database, cannot continue");
  	return ERR_FAIL;
    }
    TrajGraphWrite(ptgsHere,USE_RLE,FALSE);
    TGClose();
	USE_LOTS_RAM=ulr;
    ptgsHere = FreeTraj(ptgsHere);
  }
  freeRS(prsHead);
  Sarray = MemFree(Sarray); 
  return ERR_SUCCESS;
}


/*  
$Log: GetTrjAngle.c,v $
Revision 1.37  2004/07/08 20:38:31  hfeldman
Fixed potential memory violation when CAs are missing

Revision 1.36  2003/08/04 22:11:38  feldman
Fixed missing brace

Revision 1.35  2003/08/04 22:08:37  feldman
Added warning message for cis and fixed missing cis-check at last residue

Revision 1.34  2003/08/04 19:14:26  feldman
Fixed potential crash when minimization messes up geometry enough to turn trans into cis CA-CA distance

Revision 1.33  2003/07/18 14:39:09  feldman
Fixed uninitialized variable bug

Revision 1.32  2003/04/28 02:23:01  feldman
Fixed bug in helix fixing - for Ca walk, dont fix first residue in helix
Fixed bug when single pt in trajectory space - now prevent from 'bleeding' to adjacent cells

Revision 1.31  2003/04/04 21:53:37  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj
Fixed bug in getting Zhang potential for computing s.d. with unfoldtraj

Revision 1.30  2003/03/14 21:05:27  feldman
GetTrjAngle returns max abs. value of 0.99 for cos phi in Calpha walk to prevent singularity at the poles
Turn off USE_LOTS_RAM when writing TDs

Revision 1.29  2003/01/09 20:00:25  feldman
Allow bailing from maketrj on signal

Revision 1.28  2002/12/20 23:07:51  feldman
Correction for dt computation

Revision 1.27  2002/12/09 16:10:36  feldman
Avoid seg fault if try to homology model a non Pro onto a Cis-Pro in template, by changing residue to trans-

Revision 1.26  2002/08/01 21:54:30  feldman
Fixed some potential rounding errors

Revision 1.25  2001/04/26 19:49:09  feldman
Fixed some minor Bioseq bugs and potential bugs

Revision 1.24  2001/04/26 15:58:45  feldman
fixed disulfide placement bug

Revision 1.23  2001/03/30 22:22:27  feldman
Minor spelling and grammar fixes, bug fixes, and added
foldtraj ASCII screensaver ability

Revision 1.22  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.21  2001/03/27 20:46:00  feldman
removed unused variables

Revision 1.20  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.19  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.18  2001/03/13 15:08:16  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.17  2001/03/09 17:33:58  feldman
Added code to deal with rotamers.  Now rotamer chi angles can
be encoded into a 4-byte integer and used by foldtraj when
building the structures.  These can be edited by VisTraj as well

Revision 1.16  2001/02/15 20:27:57  feldman
Reworked maektrj to take a parameter block instead of a bunch
of individual parameters.

Revision 1.15  2001/02/06 18:40:39  feldman
Added a few functions for dealing with distance constraints and
tidied up so maketrj could join the library (foldtrajlib) without
conflicts

Revision 1.14  2001/01/26 15:41:47  feldman
Removed quiet parameter from BuildTrjGraph

Revision 1.13  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.12  2000/12/15 00:06:41  feldman
Now Val2Trj and unfoldtraj works correctly for Phi-Psi walk (I think)

Revision 1.11  2000/11/17 22:28:59  feldman
Minor bugfixes

Revision 1.10  2000/10/23 20:38:25  feldman
fixed potential crash

Revision 1.9  2000/10/06 17:33:25  feldman
Changed disulfide bridge addition slightly

Revision 1.8  2000/09/15 20:34:39  feldman
Angle and Dihedral now stored with constraints

Revision 1.7  2000/08/09 19:12:19  feldman
-minor bugfix update and fixed up makefiles removing USEDSSP
-added bioseq/seq-entry to trajectory graph header when using
 unfoldtraj or val2trj

Revision 1.6  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.5  2000/07/10 15:40:58  feldman
Updated TrajgraphWrite to not call TGInit and TGClose, thus
now just call TGInit once at start of program, TGClose at end
(also removed these calls from updatePcis and alterresiduesequence
)

Revision 1.4  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.3  2000/07/04 17:00:27  feldman
Fixed bug in maketrj which made residue one have a zero trajectory
graph

Revision 1.2  2000/06/15 17:11:03  feldman
Added Replace option when TrajGraphWrite-ing

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

