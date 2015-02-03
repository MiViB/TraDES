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


#include <stdio.h>
#include <math.h>

#include <stddef.h>
#include <stdlib.h>
#include <ctype.h> 
#include <mmdbtraj.h>

#define LSIZE 500
#define MAXLINE 100
#define DISLOCATION 8
#define OFFSET 9        /* DISLOCATION+1 */
#define WINSIZ 17       /* 2*DISLOCATION+1 */
#define NPAIRS 136      /* WINSIZ*(WINSIZ-1)/2 */
#define BLANK 21        /* The index of the Character '^' in array amino */
#define interpol_coeff 0.75
#define MINFREQ 10.
#define SKIP_BLANK 0
#define Nterm 3
#define Cterm 2
#define NR_END 1

#define ExpInv ((FloatLo) 1.0 / (FloatLo) WINSIZ);

Int4 seq_indx(Int4 c);
Int4 obs_indx(Char c);
TrajErr read_file(CharPtr fname, Int4 nprot, CharPtr *obs, CharPtr *title, Int4Ptr pnter);
TrajErr LoadGOR_DB(CharPtr fname, Int4 PNTR pnprot, CharPtr **pseq, CharPtr **ptitle, Int4Ptr *pires);
TrajErr Indices(Int4 np, Int4Ptr dis1, Int4Ptr dis2);
void predic(Int4 nres, CharPtr seq, CharPtr pred, FloatLoPtr *proba);
void First_Pass(Int4 nres, FloatLoPtr *proba, CharPtr pred);
void Second_Pass(Int4 nres, FloatLoPtr *proba, CharPtr pred);
void printout(Int4 nres, CharPtr seq, CharPtr predi, CharPtr title, FloatLoPtr *proba, FILE *fp);


/*
 * External variables
 */

static Char conf[5] = {' ','H','E','C','S'};
static FloatHi infopair[3][NPAIRS+1][23][23];
static FloatHi infodir[3][WINSIZ+1][23];
static FloatLo nS[4]/*, pS[4]*/;
static FloatLo Singlet[4][WINSIZ+1][23];
static FloatLo Doublet[4][NPAIRS+1][23][23];




/***************************************************************************/
/*                                                                         */
/* GOR secondary structure prediction method                               */
/* J. Garnier, J.-F. Gibrat, B. Robson, Methods in Enzymology,             */
/* R.F. Doolittle Ed. In press (1995)                                      */
/*                                                                         */
/* Modified by: Gil Alterovitz                                             */
/* Position: Fulbright at Dr. C. W. Hogue's Lab                            */
/* Samuel Lunenfeld Research Institute, Mt. Sinai Hospital						*/
/* Email Contact: bp837@cleveland.freenet.edu                              */
/*							                                                      */
/*                                                                         */
/* October 20, 1998                                                        */
/* Libraries: main NCBI library                                            */
/* Parameters:                                                             */
/* 1. String of AA residue characters (X is ignored).                      */
/* 2. address of 2D Array which will hold residue number x 3 probability   */
/*	values.                                                                 */
/* E.g. array[0][0] = first residue, helix probability.                    */
/*      array[0][1] = first residue, sheet probability.                     */
/*      array[0][2] = first residue, coil probability.								*/
/*                                                                         */
/*                                                                         */
/* This program gets its input from the file ggor.input                     */
/* ggor.input contains:                                                     */
/* the number of proteins in the database                                  */
/* File names must be in current directory so portable accross platforms:  */
/* the name of the file containing the sequences for the database proteins */
/* the name of the file containing the observed 2nd structures (Kabsch &   */
/* Sander classification) for the database proteins.                       */
/* H= Helix, E= Sheet, C = Coi                                             */
/*                                                                         */
/***************************************************************************/

Int4 INDMAXVAL(FloatLo val[], Int4 i1, Int4 i2)
{
/*
 * Return the index of the greatest element of val between positions i1 and i2
 */
  Int4 i, ini;

  ini = i1;
  for(i = i1+1; i <= i2; i++) {
    if(val[ini] < val[i])
      ini = i;
  }
  return(ini);
}



Int4 INDMINVAL(FloatLo val[], Int4 i1, Int4 i2)
{
/*
 * Return the index of the smallest element of val between positions i1 and i2
 */
  Int4 i, ini;

  /* printf("val[0]= %f\n",val[0]); */
  ini = i1;
  for(i = i1+1; i <= i2; i++) {
    if(val[ini] > val[i])
      ini = i;
  }
  return(ini);
}


UcharPtr cvector(Int4 nl,  Int4 nh)
{
    return UCVector(nl, nh);
}

void free_cvector(CharPtr v,  Int4 nl,  Int4 nh)
{
    UCVectorFree((UcharPtr) v, nl);
}

void CmatrixFree(CharPtr *m, Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch)
/* Free a Char matrix allocated by cmatrix() */
{
  MemFree((CharPtr) (m[nrl]+ncl-NR_END));
  MemFree((CharPtr) (m+nrl-NR_END));
}

CharPtr *Cmatrix(Int4 nrl, Int4 nrh, Int4 ncl, Int4 nch)
/* Allocate a Char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	Int4 i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	CharPtr *m;

	/* allocate pointers to rows */
	m = (CharPtr *) MemNew((size_t) ((nrow+NR_END) * sizeof(CharPtr)));
	if(!m) {
		ErrPostEx(SEV_ERROR, 0, 0, "Out of Memory 1 in cmatrix");
		ErrShow();
		return NULL;
	} 
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (CharPtr) MemNew((size_t) ((nrow*ncol+NR_END) * sizeof(Char)));
	if(!m[nrl]) if(!m) {
		ErrPostEx(SEV_ERROR, 0, 0, "Out of Memory 2 in cmatrix");
		ErrShow();
		return NULL;
	}  
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i = nrl + 1; i <= nrh; i++) {
		m[i] = m[i-1] + ncol;
	}
	/* return pointer to array of pointers to rows */
	return m;
}


/*****************************************************************************/
/*                                                                           */
/* This routine reads the sequence and observed secondary structures for all */
/* the proteins in the data base.                                            */
/*                                                                           */
/*****************************************************************************/
TrajErr LoadGOR_DB(CharPtr fname, Int4 PNTR pnprot, CharPtr **pseq, CharPtr **ptitle, Int4Ptr *pires)
{
	FILE *fp;
	Int4 ip = 0, nres = 0, nprot = 0, i = 0;
	Char c;
	Char buf[MAXLINE];
	Char keep[MAXRES];

	CharPtr *seq = NULL, *title = NULL;
	Int4Ptr ires = NULL;


	if((fp = FileOpen(fname,"r")) == NULL) {
		ErrPostEx(SEV_ERROR, 1,5,"LoadGOR_DB: Could not find file %s\n",fname);
		return ERR_FAIL;
	}

	/* Get the number of proteins */
	while(FileGets(buf,MAXLINE,fp) != NULL) {
		while((c = (Char)getc(fp)) != '@') ;
		while((c = (Char)getc(fp)) != '\n');
		nprot++;
	}
	FileClose(fp);


	/* Allocate Memory */
	if((seq = Cmatrix(1,nprot,1,MAXRES)) == NULL) return ERR_FAIL;
	if((title = (CharPtr *) MemNew((size_t)sizeof(CharPtr)*(nprot+1))) == NULL) return ERR_FAIL;
	if((ires = I4Vector(1,nprot)) == NULL) return ERR_FAIL;

	if((fp = FileOpen(fname,"r")) == NULL) return ERR_FAIL;
	
	for(ip = 1; ip <= nprot; ip++) {
		FileGets(buf,MAXLINE,fp);
		title[ip] = StringSave(buf);

		nres = 0;
		while((c = (Char)getc(fp)) != '@') {
			if(c == '\n' || c == ' ' || c =='\t') continue;
				nres++;
			if(nres > MAXRES) {
				ErrPostEx(SEV_ERROR, 1,6,"The value of MAXRES should be increased: %ld", (long) MAXRES);
				FileClose(fp);
				return ERR_FAIL;
			}
			if((c >= 'A' && c < 'Z') && c != 'J' && c != 'O' && c != 'U') {
				keep[nres] = c;
			} else {
				ErrPostEx(SEV_ERROR, 1,7,"Error in protein: %ld residue: %ld\nInvalid amino acid type or secondary structure state in character: %c\n",
				(long) ip, (long) nres, c);
				keep[nres] = 'X';
			}
		}
		while((c = (Char)getc(fp)) != '\n')
		;
		for(i = 1; i <= nres; i++) {
			seq[ip][i] = keep[i];
		}
		ires[ip] = nres;
	}
	FileClose(fp);

	*pnprot = nprot;
	*pseq = seq;
	*ptitle = title;
	*pires = ires;
	return ERR_SUCCESS;
}

void FreeGOR_DB(Int4 nprot, CharPtr **pseq, CharPtr **ptitle, Int4Ptr *pires)
{
	CmatrixFree(*pseq,1,nprot,1,MAXRES);
	CmatrixFree(*ptitle,1,nprot,1,MAXLINE);
	I4VectorFree(*pires,1);
	*pseq = NULL;
	*ptitle = NULL;
	*pires = NULL;
}


/*****************************************************************************/
/*                                                                           */
/* This routine reads the sequence and observed secondary structures for all */
/* the proteins in the data base.                                            */
/*                                                                           */
/*****************************************************************************/
TrajErr read_file(CharPtr fname, Int4 nprot, CharPtr *obs, CharPtr *title, Int4Ptr pnter)
{
  FILE *fp;
  Int4 ip, nres, i;
  Char c;
  CharPtr keep;

  fp = FileOpen(fname,"r");
  if(fp == NULL) {
    ErrPostEx(SEV_ERROR, 1,5,"Could not find file %s\n",fname);
    return ERR_FAIL;
  }

  keep = (CharPtr) MemNew((size_t) MAXRES*sizeof(Char));

  for(ip = 1; ip <= nprot; ip++) {
    FileGets(title[ip],MAXLINE,fp);
    nres = 0;
    while((c = (Char)getc(fp)) != '@') {
      if(c == '\n' || c == ' ' || c =='\t') continue;
      nres++;
      if(nres > MAXRES) {
		ErrPostEx(SEV_ERROR, 1,6,"The value of MAXRES should be increased: %ld", (long) MAXRES);
		MemFree(keep);
		FileClose(fp);
		return ERR_FAIL;
      }
      if((c >= 'A' && c < 'Z') && c != 'B' && c != 'J' && c != 'O' && c != 'U') {
		keep[nres] = c;
      }
      else {
		ErrPostEx(SEV_ERROR, 1,7,"Error in protein: %ld residue: %ld\nInvalid amino acid type or secondary structure state in character: %c\n",
		 (long) ip, (long) nres, c);
		MemFree(keep);
		FileClose(fp);
        return ERR_FAIL;
      }
    }
    while((c = (Char)getc(fp)) != '\n')
      ;
    for(i = 1; i <= nres; i++)
      obs[ip][i] = keep[i];
    pnter[ip] = nres;
  }

  MemFree(keep);

  FileClose(fp);
  return ERR_SUCCESS;

}
/*****************************************************************************/
/*                                                                           */
/* This function returns an integer for each amino acid type.                */
/*                                                                           */
/*****************************************************************************/
Int4 seq_indx(Int4 c)
{
  switch(c) {
  case 'A': return(1);
  case 'C': return(2);
  case 'D': return(3);
  case 'E': return(4);
  case 'F': return(5);
  case 'G': return(6);
  case 'H': return(7);
  case 'I': return(8);
  case 'K': return(9);
  case 'L': return(10);
  case 'M': return(11);
  case 'N': return(12);
  case 'P': return(13);
  case 'Q': return(14);
  case 'R': return(15);
  case 'S': return(16);
  case 'T': return(17);
  case 'V': return(18);
  case 'W': return(19);
  case 'Y': return(20);
  case '^': return(21);
  case '-': return(22);
  default : return(23);
  }

}
/*****************************************************************************/
/*                                                                           */
/* This function returns an integer for each secondary structure type.       */
/*                                                                           */
/*****************************************************************************/
Int4 obs_indx(Char c)
{

  switch(c) {
  case 'H': return(1);
  case 'E': return(2);
  case 'C': return(3);
  case 'X': return(0);
  default:
	return(0);
  }

}


/*****************************************************************************/
/*                                                                           */
/* Determine indices dis1 dis2 as a function of np                           */
/*                                                                           */
/*****************************************************************************/
TrajErr Indices(Int4 np, Int4Ptr dis1, Int4Ptr dis2)
{
  Int4 i, j, k;

  k = 0;
  for(i = -DISLOCATION; i <= DISLOCATION; i++) {
    for(j= i+1; j <= DISLOCATION; j++) {
      k++;
      if(k == np) {
	*dis1 = i;
	*dis2 = j;
	return ERR_SUCCESS;
      }
    }
  }
  ErrPostEx(SEV_ERROR, 1,8,"Error invalid value of np= %ld\n", (long) np);
  return ERR_FAIL;
}

/*********************************************************************************/
/*                                                                               */
/*                          Normalize the probabilities                          */
/*                                                                               */
/*********************************************************************************/
void GORNormalize(FloatLoPtr proba, FloatHiPtr v)
{
  FloatHi denom;
  
  denom = (FloatHi) (1.0 / (1.0 + exp(v[1]) + exp(v[2])));
  proba[1] = (FloatLo) (exp(v[1]) * denom);
  proba[2] = (FloatLo) (exp(v[2]) * denom);
  proba[3] = (FloatLo) denom;
}

/***********************************************************************************************************/
/*                                                                                                         */
/* This routine performs the prediction of the current protein                                             */
/*                                                                                                         */
/***********************************************************************************************************/
void predic(Int4 nres, CharPtr seq, CharPtr pred, FloatLoPtr *proba)
{
	FloatHi it[3];
	Int4 aa1, aa2;
	Int4 konf, ires;
	Int4 dis1, dis2, np;

	/*
	* Calculate sum of information values for each secondary structure type (konf)
	*/

	for(ires = 1; ires <= nres; ires++) {  
		it[1] = it[2] = (FloatHi) 0.0;
		for(dis1 = -DISLOCATION; dis1 <= +DISLOCATION; dis1++) {
			if(ires+dis1 < 1 || ires +dis1 > nres) {
				if(SKIP_BLANK) continue;                 /* If SKIP_BLANK "amino acid" of type ' ', i.e., */
				aa1 = 21;                                /* aa1 = 21 are not included in the calculation  */
			} else {
				aa1 = seq_indx(seq[ires+dis1]);
			}
			for(dis2 = dis1+1; dis2 <= +DISLOCATION; dis2++) {
				if(ires+dis2 < 1 || ires +dis2 > nres) {
					if(SKIP_BLANK) continue;                 
					aa2 = 21;                                
				} else {
					aa2 = seq_indx(seq[ires+dis2]);
				}
				np = (dis1+8) * (WINSIZ-1) - ((dis1+8)*(dis1+9)/2) + (dis2+8);
				for(konf = 1; konf <= 2; konf++)  {
					it[konf] = it[konf] + infopair[konf][np][aa1][aa2];
				}
			}
		}
		for(dis1 = -DISLOCATION; dis1 <= +DISLOCATION; dis1++) {
			if(ires+dis1 < 1 || ires +dis1 > nres) {
				if(SKIP_BLANK) continue;                 
				aa1 = 21;                                
			} else {
				aa1 = seq_indx(seq[ires+dis1]);
			}
			for(konf = 1; konf <= 2; konf++) {
				it[konf] = it[konf] + infodir[konf][dis1+9][aa1];
			}
		}
		GORNormalize((FloatLoPtr) proba[ires],it);
		pred[ires] = conf[INDMAXVAL(proba[ires],1,3)];
	}

	/*
	* If "blank residues" are not included the first Nterm and the last Cterm residues are predicted as coils
	*/

	if(SKIP_BLANK) {                              
		for(ires = 1; ires <= Nterm; ires++) {
			pred[ires] = 'C';
		}
		for(ires = nres-Cterm+1; ires <= nres; ires++) {
			pred[ires] = 'C';
		}
	}

}
  

/***************************************************************************/
/*                                                                         */
/*                           Routine First_Pass                            */
/*                                                                         */
/***************************************************************************/
void First_Pass(Int4 nres, FloatLo **proba, Char *pred)
{
/*
 * 1) Look for areas that are a mixture of Es and Hs.
 * 2) When such an area is isolated check whether Es and Hs occurs in two blocks.
 * If yes and number of Hs > 4 and number of Es > 3 do nothing
 * In all other cases compute the product of probabilities for all residues in the area
 * and assign to this area the conformation having the highest probability over the area.
 *
 */

  Int4 ires;
  Int4 lim1=0, lim2;
  Int4 open;
  Int4 kk;
  Int4 type;
  Int4 block[3];
  Int4 nseg;
  Int4 size[3] = {0,4,3};
  FloatHi ptot[3];

  pred[1] = pred[nres] = 'C';
  open = 0;
  for(ires = 1; ires <= nres; ires++) {
    if(pred[ires] != 'C') {
      if(!open) {
	open = 1;
	lim1 = ires;
      }
    } else {
      if(open) {
	open = 0;
	lim2 = ires - 1;
	type = obs_indx(pred[lim1]);
	block[1] = block[2] = 0;
	nseg = 1;
	block[nseg]++;
	for(kk = lim1+1; kk <= lim2; kk++) {
	  if(obs_indx(pred[kk]) != type)
	    nseg++;
	  if(nseg <= 2) block[nseg]++;
	  type = obs_indx(pred[kk]);
	}
	if(nseg > 2 || block[1] < size[obs_indx(pred[lim1])] || block[2] < size[obs_indx(pred[lim2])]) {
	  ptot[1] = ptot[2] = (FloatHi) 1.0;
	  for(kk = lim1; kk <= lim2; kk++) {
	    ptot[1] = ptot[1] * proba[kk][1];
	    ptot[2] = ptot[2] * proba[kk][2];
	  }
	  if(ptot[1] > ptot[2]) {
	    for(kk = lim1; kk <= lim2; kk++)
	      pred[kk] = 'H';
	  } else {
	    for(kk = lim1; kk <= lim2; kk++)
	      pred[kk] = 'E';
	  }
	}
      }
    }
  }

}
/***************************************************************************/
/*                                                                         */
/*                           Routine Second_Pass                           */
/*                                                                         */
/***************************************************************************/
void Second_Pass(Int4 nres, FloatLo **proba, Char *pred)
{
/*
 * Correct strands having less than 2 and helices having less than 4 residues.
 * Either the secondary structure element is suppressed or additional
 * residues are recruted to reach the required number.
 * 
 */

  Int4 ires, ires1;
  Int4 len;
  Int4 standard[4] = {0,4,2,0};
  Int4 missing;
  Int4 k;
  Int4 lim1, lim2, lim3, lim4, Lim1=0, Lim2=0, Lim3=0, Lim4=0, KeepNterm=0, KeepCterm=0;
  FloatLo cost, costmax;
  Int4 type;
  Int4 type_Cterm, type_Nterm;
 
  len = 0;
  type = obs_indx(pred[1]);
  for(ires = 2; ires <= nres; ires++) {
    if(type != obs_indx(pred[ires])) {
      if(len < standard[type]) { /* Check all possibilities */
	costmax = (FloatLo) 0.0;
	missing = standard[type] - len;
/*
 * Check the cost of increasing the secondary structure element
 */
	lim1 = ires - len - missing;
	for(k = 1; k <= missing+1; k++) {
	  lim2 = lim1 + standard[type] - 1;
	  if(lim1 < 1 || lim2 > nres) {
	    lim1++;
	    continue;
	  }
	  cost = (FloatLo) 1.0;
	  for(ires1 = lim1; ires1 <= lim2; ires1++)
	    cost *= proba[ires1][type];
	  if(cost > costmax) {
	    costmax = cost;
	    Lim1 = lim1;
	    Lim2 = lim2;
	    KeepNterm = type;
	    Lim3 = 0;
	    Lim4 = -1;
	  }
	  lim1++;
	}
/*
 * Check the cost of suppressing the secondary structure element using the same segments as previously
 */
	type_Nterm = obs_indx(pred[ires-len-1]);
	type_Cterm = obs_indx(pred[ires]);
	lim1 = ires - len - missing;
	for(k = 1; k <= missing+1; k++) {
	  lim4 = lim1 + standard[type] - 1;
	  if(lim1 < 1 || lim4 > nres) {
	    lim1++;
	    continue;
	  }
	  lim2 = ires - 1;
	  lim3 = lim2 + 1;
	  while(lim3 >= ires - len) {
	    cost = (FloatLo) 1.0;
	    for(ires1 = lim1; ires1 <= lim2; ires1++)
	      cost *= proba[ires1][type_Nterm];
	    for(ires1 = lim3; ires1 <= lim4; ires1++)
	      cost *= proba[ires][type_Cterm];
	    if(cost > costmax) {
	      costmax = cost;
	      Lim1 = lim1;
	      Lim2 = lim2;
	      Lim3 = lim3;
	      Lim4 = lim4;
	      KeepNterm = type_Nterm;
	      KeepCterm = type_Cterm;
	    }
	    lim2--;
	    lim3--;
	  }
	  lim1++;
	}
/*
 * Modify pred accordingly
 */
	for(ires1 = Lim1; ires1 <= Lim2; ires1++)
	  pred[ires1] = conf[KeepNterm];
	for(ires1 = Lim3; ires1 <= Lim4; ires1++)
	  pred[ires1] = conf[KeepCterm];
/*
 * Move to the end of the modified segment if necessary
 */
	if(Lim2 > ires || Lim4 > ires) {
	  if(Lim2 > Lim4)
	    ires = Lim2;
	  else 
	    ires = Lim4;
	}

      } /* End of segment correction */
 
      len = 1;
    } else {
      len++;
    }
    type = obs_indx(pred[ires]);
  }

}



/* Compute the frequencies from proteins in the data base. */
TrajErr GOR_Parameters(Int4 nprot_dbase, Int4Ptr nres, CharPtr *obs, CharPtr *seq)
{
	Int4 pro;
	Int4 ires;
	Int4 konf, dis, aa1, aa2, np;
	Int4 dis1, dis2;
	FloatHi C1, C2;
	FloatLo f1, f2, f3;

	C1 = (FloatHi) 2. * (FloatHi) ExpInv;
	C2 = (FloatHi) 1. - C1;


	/*
	* Initialisation
	*/

	for(konf = 0; konf < 4; konf++) {
		for(dis = 0; dis < WINSIZ+1; dis++){
			for(aa1 = 0; aa1 < 23; aa1++) {
				Singlet[konf][dis][aa1] = (FloatLo) 0.0;
	}}}

	for(konf = 0; konf < 4; konf++) {
		for(np = 0; np < NPAIRS+1; np++) {
			for(aa1 = 0; aa1 < 23; aa1++) {
				for(aa2 = 0; aa2 < 23; aa2++) {
					Doublet[konf][np][aa1][aa2] = (FloatLo) 0.0;
	}}}}

	nS[0] = nS[1] = nS[2] = nS[3] = 0;

	/*
	* Loop over all the proteins of the data base. 
	*/

	for(pro = 1; pro <= nprot_dbase; pro++) {
		/* Determine frequencies related to the sequence of the query protein (the 1st row in the alignment) */
		for(ires = 1; ires <= nres[pro]; ires++) {
			konf = obs_indx(obs[pro][ires]);
			if(konf == 0) {                    /* Skip X conformations, i.e., residues for */
				continue;  /* which the secondary structure is unknown */
			} 
			nS[konf]++;   

			for(dis = -DISLOCATION; dis <= DISLOCATION; dis++) {
				if(ires+dis < 1 || ires+dis > nres[pro]) {
					aa1 = BLANK;
				} else {
					aa1 = seq_indx(seq[pro][ires+dis]);
					Singlet[konf][dis+OFFSET][aa1] += (FloatLo) 1.0;
				}
			}

			np = 0;
			for(dis1 = -DISLOCATION; dis1 <= DISLOCATION; dis1++) {
				if(ires+dis1 < 1 || ires+dis1 > nres[pro]) {
					aa1 = BLANK;
				} else {
					aa1 = seq_indx(seq[pro][ires+dis1]);
				}
				for(dis2 = dis1+1; dis2 <= DISLOCATION; dis2++) {
					if(ires+dis2 < 1 || ires+dis2 > nres[pro]) {
						aa2 = BLANK; 
					} else {
						aa2 = seq_indx(seq[pro][ires+dis2]); 
					}
					np++;
					Doublet[konf][np][aa1][aa2] += (FloatLo) 1.0;
	}}}} /* End of loop over the proteins in the data base index pro */

	/*
	* Calculate probabilities for the 3 secondary structures, H, E and C.
	*/

	nS[0] = nS[1] + nS[2] + nS[3];

	/*  for(konf = 1; konf <= 3; konf++)
	pS[konf] = (FloatLo) nS[konf] / (FloatLo) nS[0];
	*/
	/*
	* Calculate information parameters (sort of)
	*/

	for(konf = 1; konf <= 2; konf++) {
		for(np = 1; np <= NPAIRS; np++) {
			for(aa1 = 1; aa1 <= 21; aa1++) {
				for(aa2 = 1; aa2 <= 21; aa2++) {
					f1 = Doublet[konf][np][aa1][aa2];
					f2 = Doublet[3][np][aa1][aa2];
					if(f1 < MINFREQ) {
						if (Indices(np,&dis1,&dis2)!=ERR_SUCCESS) {
							return ERR_FAIL;
						}
						f3 = Singlet[konf][dis1][aa1] * Singlet[konf][dis2][aa1] / (FloatLo) nS[konf];
						f1 = (f3 - f1) * (FloatLo) interpol_coeff + f1;
						if(f1 < 1.e-6) {
							f1 = (FloatLo) 1.0;
						}
					}
					if(f2 < MINFREQ) {
						if (Indices(np,&dis1,&dis2)!=ERR_SUCCESS) {
							return ERR_FAIL;
						}
						f3 = Singlet[3][dis1][aa1] * Singlet[3][dis2][aa1] / (FloatLo) nS[3];
						f2 = (f3 - f2) * (FloatLo) interpol_coeff + f2;
						if(f2 < 1.e-6) {
							f2 = (FloatLo) 1.0;
						}
					}
					infopair[konf][np][aa1][aa2] =  C1 * (FloatHi) (log(f1)-log(f2));
	}}}}


	for(konf = 1; konf <= 2; konf++) {
		for(dis = 1; dis <= WINSIZ; dis++) {
			for(aa1 = 1; aa1 <= 21; aa1++) {
				f1 = Singlet[konf][dis][aa1];
				f2 = Singlet[3][dis][aa1];
				if(f1 < 1.e-6) { 
					f1 = (FloatLo) 1.0;
				}
				if(f2 < 1.e-6) {
					f2 = (FloatLo) 1.0;
				}
				infodir[konf][dis][aa1] = C2 * (FloatHi) (log(f2)- log(f1));
	}}}
	
	return ERR_SUCCESS;
}



/*******************************************************************************************************
* EGOR_DJK                                                                                             *
********************************************************************************************************
  This function creates a probability matrix based on the GOR method for secondary structure prediction. 
  It uses a given protein sequence and its corresponding secondary structure to create this matrix.  
  The secondary structure of another sequence that is aligned to the first sequence and the SS can be 
  predicted and the success rate estimated using the double jackknife statistic. 
  Up to 3 "database" files must be supplied, each of which must have an identical number of records,
  each of identical length:  1) a secondary structure file (database_*.ss), containing the correct
  secondary structure for each amino acid position denoted by H (Helix), E(Sheet), and C (coil), 2) 
  a FASTA-like sequence database (database_*.str), that corresponds to the sequence underlying the 
  secondary structure in 1) and finally, 3) some other, aligned sequence to 2) and 1).  Each file is 
  formatted as follows: a header starting with '!' followed by a definition line, followed by a
  capitalized protein sequence  or secondary structure and the 'record' ends with '@'.  None of the
  sequences should contain masked out 'X' or gapped '-' residues due to some alignment process. 
  This function uses one sequence database along with the secondary structure database as the base
  for the probability matrix.  It then takes out each sequence out of the database, one at a time,
  and predicts either that sequence (normal jackknife) or the other sequence.  The predictive success rate 
  is tabulated for all entries and written to a file "gor[d]jk_[str|seq]db_%ld.txt" dependending 
  on which was chosen for the jackknife test.
******************************************************************************************************/
TrajErr EGORDoubleJacknife(Int4 dbID, Boolean bDJK, Boolean bUseStrDB, FILE *fpdetailedresults)
{
/*	FloatLo arwProbMatrix[MAXRES][MAXSTRUCT];
	Char arwFiltStruct[MAXRES];
*/	
	Int4 nprot = 0;
	Int4 nprot_dbase = 0;
	Int4 i = 0, j = 0, k = 0;
	Int4 pro = 0;
	Int4 iw2 = 0,iPred = 0,wA = 0,wB = 0,wC = 0;	
	
	Int4 wCorrectNRES = 0;
	Int4 wNumXs = 0;
	Int4 wTotalRes = 0;
	Int4 wTotalPred = 0;
	
	Int4 wHHObsPred = 0;
	Int4 wHEObsPred = 0;
	Int4 wHCObsPred = 0;
	Int4 wEHObsPred = 0;
	Int4 wEEObsPred = 0;
	Int4 wECObsPred = 0;
	Int4 wCHObsPred = 0;
	Int4 wCEObsPred = 0;
	Int4 wCCObsPred = 0;
	Int4Ptr temp = NULL, temp_d = NULL, nres = NULL, nres_temp = NULL, nres_d = NULL, NRES = NULL;
	
	Char f1[PATH_MAX];
	Char f2[PATH_MAX];
	Char fss[PATH_MAX];
	Char fresults[PATH_MAX];
	CharPtr cPredOrig = NULL;
	CharPtr predi = NULL;
	CharPtr *obs = NULL, *seq = NULL, *OBS = NULL, *SEQ = NULL;
	CharPtr *obs_d = NULL, *seq_d = NULL, *obs_temp = NULL, *seq_temp = NULL;
	CharPtr *title_obs = NULL,  *title_seq = NULL,  *TITLE = NULL;
	CharPtr *title_obs_d = NULL, *title_seq_d = NULL;
	FloatLo fAval= 0.0,fBval= 0.0,fCval= 0.0,fDval= 0.0;
	FloatLo fprob = 0.0;
	FloatLo fTotalProb = 0.0;
	FloatLoPtr *probai = NULL;
	
	FILE *fp = NULL;


	/* so, we get to choose whether we want to create a database from the structure sequences or from the genome sequence */

	if(bDJK && bUseStrDB) {
		sprintf(f1,"database_%ld.str",(long) dbID);
		sprintf(f2,"database_%ld.seq",(long) dbID);
		sprintf(fresults,"gordjk_strdb_%ld.txt",(long) dbID);
	} else if (bDJK && !bUseStrDB) {
		sprintf(f1,"database_%ld.seq",(long) dbID);
		sprintf(f2,"database_%ld.str",(long) dbID);
		sprintf(fresults,"gordjk_seqdb_%ld.txt",(long) dbID);
	} else if (!bDJK && bUseStrDB) {
		sprintf(f1,"database_%ld.str",(long) dbID);
		sprintf(f2,"database_%ld.str",(long) dbID);
		sprintf(fresults,"gorjk_strdb_%ld.txt",(long) dbID);
	} else if (!bDJK && !bUseStrDB) {
		sprintf(f1,"database_%ld.seq",(long) dbID);
		sprintf(f2,"database_%ld.seq",(long) dbID);
		sprintf(fresults,"gorjk_seqdb_%ld.txt",(long) dbID);
	}
	sprintf(fss,"database_%ld.ss",(long) dbID);
	
	if((fp = FileOpen(fresults, "w")) == NULL) return ERR_FAIL;
	if(fpdetailedresults != NULL) fprintf(fpdetailedresults, "Protein | Length  |  Number Predicted  |  Structure Prob. Sum\n");

/*
 * Input the sequences and secondary structures for the data base
 */
	if (LoadGOR_DB(f1,&nprot_dbase, &seq, &title_seq, &temp)!=ERR_SUCCESS) return ERR_FAIL;
	if (LoadGOR_DB(fss,&nprot_dbase, &obs, &title_obs, &nres)!=ERR_SUCCESS) return ERR_FAIL;
	if (LoadGOR_DB(f2,&nprot_dbase, &seq_d, &title_seq_d, &temp_d)!=ERR_SUCCESS) return ERR_FAIL;
	if (LoadGOR_DB(fss,&nprot_dbase, &obs_d, &title_obs_d, &nres_d)!=ERR_SUCCESS) return ERR_FAIL;

	if((seq_temp = Cmatrix(1,nprot_dbase,1,MAXRES)) == NULL) return ERR_FAIL;
	if((obs_temp = Cmatrix(1,nprot_dbase,1,MAXRES)) == NULL) return ERR_FAIL;
	if((nres_temp = I4Vector(1,nprot_dbase)) == NULL) return ERR_FAIL;

	/* for JK */
	nprot = 1;
	if((SEQ = Cmatrix(1,nprot,1,MAXRES)) == NULL) return ERR_FAIL;
	if((OBS = Cmatrix(1,nprot,1,MAXRES)) == NULL) return ERR_FAIL;
	if((TITLE = Cmatrix(1,nprot,1,MAXLINE)) == NULL) return ERR_FAIL;
	if((NRES = I4Vector(1,nprot)) == NULL) return ERR_FAIL;

	
	if((predi = (CharPtr) cvector(1,MAXRES)) == NULL) return ERR_FAIL;
	if((cPredOrig = (CharPtr) cvector(1,MAXRES)) == NULL) return ERR_FAIL;
	if((probai = FLMatrix(1,MAXRES,1,3)) == NULL) return ERR_FAIL;

	for (i = 1;i <= nprot_dbase;i++) {
		for(j = 1; j <= nres_d[i]; j++) {
			OBS[1][j] = obs_d[i][j];
			SEQ[1][j] = seq_d[i][j];
		}
		NRES[1] = nres_d[i];
				
		for (pro = 1, k = 1; pro <= nprot_dbase; pro++) {
			if(pro == i) continue; 
			nres_temp[k] = nres[pro];
			for (iw2=1;iw2<=nres[pro];iw2++) {
				seq_temp[k][iw2] = seq[pro][iw2];
				obs_temp[k][iw2] = obs[pro][iw2];
			}
			k++;
		}

		if (GOR_Parameters(nprot_dbase-1,nres_temp,obs_temp,seq_temp)!=ERR_SUCCESS) {
			return ERR_FAIL;
		}

		/*
		* Predict the secondary structure of protein pro.
		* Carry out the prediction for the sequence alone 
		*/
		for(pro = 1; pro <= nprot; pro++) {
			predic(NRES[pro],SEQ[pro],predi,probai);
			for (iPred=1;iPred<=NRES[pro];iPred++) {
				cPredOrig[iPred]=predi[iPred];
			}
			First_Pass(NRES[pro],probai,predi);
			Second_Pass(NRES[pro],probai,predi);
			for (iPred=1;iPred<=NRES[pro];iPred++) {
				if (cPredOrig[iPred]!=predi[iPred]) {
					wA = obs_indx(cPredOrig[iPred]);
					wB = obs_indx(predi[iPred]);
					wC = 6-wA-wB;
					fDval = probai[iPred][wA]-probai[iPred][wB];
					fBval = probai[iPred][wA];	/* B = B+D = A */				                
					fAval = probai[iPred][wA] - fDval*(probai[iPred][wC]/(probai[iPred][wA]  + probai[iPred][wC]));	/* A=A-D*(C/(A+C)) */   
					fCval = probai[iPred][wC] - fDval*(probai[iPred][wA]/(probai[iPred][wA] + probai[iPred][wC]));	/* C=C-D*(A/(A+C)) */
					/* if C<0, then { A=A+C; C=0 } (note we're adding a negative number to A) */
					if (fCval < 0) {
						fAval = fAval + fCval;
						fCval = 0;               
					}
					probai[iPred][wA] = fAval;
					probai[iPred][wB] = fBval;
					probai[iPred][wC] = fCval;
				}    
			} 

			wTotalRes+=NRES[pro];
			wCorrectNRES=0;
			fprob=0.0;

			for(iw2 = 1; iw2 <= NRES[pro]; iw2++) {
				wNumXs=0;

				if (OBS[pro][iw2] == 'H') {
					fprob+=probai[iw2][1];
				}
				if (OBS[pro][iw2] == 'E') {
					fprob+=probai[iw2][2];
				}
				if (OBS[pro][iw2] == 'C') {
					fprob+=probai[iw2][3];
				}


				if ((predi[iw2] == 'H') && (OBS[pro][iw2] == 'H')) {
					wHHObsPred++;
					wCorrectNRES++;
				}
				if ((predi[iw2] == 'H') && (OBS[pro][iw2] == 'E')) {
					wHEObsPred++;
				}
				if ((predi[iw2] == 'H') && (OBS[pro][iw2] == 'C')) {
					wHCObsPred++;
				}

				if ((predi[iw2] == 'E') && (OBS[pro][iw2] == 'H')) {
					wEHObsPred++;
				}
				if ((predi[iw2] == 'E') && (OBS[pro][iw2] == 'E')) {
					wEEObsPred++;
					wCorrectNRES++;
				}
				if ((predi[iw2] == 'E') && (OBS[pro][iw2] == 'C')) {
					wECObsPred++;
				}

				if ((predi[iw2] == 'C') && (OBS[pro][iw2] == 'H')) {
					wCHObsPred++;
				}
				if ((predi[iw2] == 'C') && (OBS[pro][iw2] == 'E')) {
					wCEObsPred++;
				}
				if ((predi[iw2] == 'C') && (OBS[pro][iw2] == 'C')) {
					wCCObsPred++;
					wCorrectNRES++;
				}

				if (predi[iw2] == OBS[pro][iw2]) {
   	      			wTotalPred++;
				}
	
				if (OBS[pro][iw2]=='X') {
					wNumXs++;
					printf("!!!! ***** %ld", (long) wTotalRes);
				}
				wTotalRes-=wNumXs;
			}
			if(fpdetailedresults) fprintf(fpdetailedresults, "%ld %ld %ld %f\n", (long) i, (long) (NRES[pro] - wNumXs), (long) wCorrectNRES, (float) fprob);
			fTotalProb+=fprob;
		}

/*
		for (iw2=0;iw2<NRES[1];iw2++) {
			arwFiltStruct[iw2] = predi[iw2+1];
			arwProbMatrix[iw2][0] = probai[iw2+1][1];
  			arwProbMatrix[iw2][1] = probai[iw2+1][2];
  			arwProbMatrix[iw2][2] = probai[iw2+1][3];
		}
*/
	}

	printf("\nTotal Pred. Residues: %ld, Total Residues: %ld, Q3 = %f\n", (long) wTotalPred, (long) wTotalRes, (float) ((FloatLo) wTotalPred / (FloatLo) wTotalRes));
	printf("Information Function: %f\n", (float) ((FloatLo) fTotalProb / (FloatLo) wTotalRes));
	printf("Qh = %f; Qe = %f; Qc = %f\n",(float) ((FloatLo) wHHObsPred / (FloatLo) (wHHObsPred+wEHObsPred+wCHObsPred)), (float) ((FloatLo) wEEObsPred / (FloatLo) (wHEObsPred+wEEObsPred+wCEObsPred)), (float) ((FloatLo) wCCObsPred / (FloatLo) (wHCObsPred+wECObsPred+wCCObsPred)) );

	printf("HH: %ld, HE: %ld, HC: %ld\n", (long) wHHObsPred,  (long) wHEObsPred, (long) wHCObsPred);
	printf("EH: %ld, EE: %ld, EC: %ld\n", (long) wEHObsPred,  (long) wEEObsPred, (long) wECObsPred);
	printf("CH: %ld, CE: %ld, CC: %ld\n", (long) wCHObsPred,  (long) wCEObsPred, (long) wCCObsPred);

	if(fp) {
		fprintf(fp,"Total Pred. Residues: %ld, Total Residues: %ld, Q3 = %f\n", (long) wTotalPred, (long) wTotalRes, (float) wTotalPred / (float) wTotalRes);
		fprintf(fp,"Information Function: %f\n",  ((FloatLo) fTotalProb / (FloatLo) wTotalRes));
		fprintf(fp,"Qh = %f; Qe = %f; Qc = %f\n", (float) ((FloatLo) wHHObsPred / (FloatLo) (wHHObsPred+wEHObsPred+wCHObsPred)), (float) ((FloatLo) wEEObsPred / (FloatLo) (wHEObsPred+wEEObsPred+wCEObsPred)), (float) ((FloatLo) wCCObsPred / (FloatLo) (wHCObsPred+wECObsPred+wCCObsPred)) );
		fprintf(fp,"HH: %ld, HE: %ld, HC: %ld\n", (long) wHHObsPred,  (long) wHEObsPred, (long) wHCObsPred);
		fprintf(fp,"EH: %ld, EE: %ld, EC: %ld\n", (long) wEHObsPred,  (long) wEEObsPred, (long) wECObsPred);
		fprintf(fp,"CH: %ld, CE: %ld, CC: %ld\n", (long) wCHObsPred,  (long) wCEObsPred, (long) wCCObsPred);

		FileClose(fp);
	}


	/*
	* Free memory
	*/

	CmatrixFree(seq,1,nprot_dbase,1,MAXRES);
	CmatrixFree(seq_d,1,nprot_dbase,1,MAXRES);
	CmatrixFree(seq_temp,1,nprot_dbase,1,MAXRES);

	CmatrixFree(obs,1,nprot_dbase,1,MAXRES);
	CmatrixFree(obs_d,1,nprot_dbase,1,MAXRES);
	CmatrixFree(obs_temp,1,nprot_dbase,1,MAXRES);

	CmatrixFree(title_obs,1,nprot_dbase,1,MAXLINE);
	CmatrixFree(title_obs_d,1,nprot_dbase,1,MAXLINE);
	CmatrixFree(title_seq,1,nprot_dbase,1,MAXLINE);
	CmatrixFree(title_seq_d,1,nprot_dbase,1,MAXLINE);

	CmatrixFree(SEQ,1,nprot,1,MAXRES);
	CmatrixFree(OBS,1,nprot,1,MAXRES);
	CmatrixFree(TITLE,1,nprot,1,MAXLINE);

	I4VectorFree(temp,1);
	I4VectorFree(nres,1);
	I4VectorFree(nres_d,1);
	I4VectorFree(NRES,1);
	
	free_cvector(predi,1,MAXRES);
	free_cvector(cPredOrig,1,MAXRES);
	FLMatrixFree(probai,1,1);

	return ERR_SUCCESS;
}


TrajErr fnGilGOR(Char *szInputSeq, FloatLo arwProbMatrix[MAXRES][MAXSTRUCT], Char arwFiltStruct[MAXRES])
{


  Int4 nprot =  0;
  Int4 nprot_dbase = 0;
  Int4Ptr temp, nres, NRES;
  Int4 i;
  Int4 pro;
  Int4 nerr;
  CharPtr predi;
  Char fname1[PATH_MAX];
  Char fname2[PATH_MAX];
  FloatLoPtr *probai;
  CharPtr cPredOrig;
  Int4 iw2,iPred,wA,wB,wC;
  CharPtr *obs, *seq, *SEQ;
  CharPtr *title_obs,  *title_seq,  *TITLE;
  FloatLo fAval,fBval,fCval,fDval;

/* added by Howie */
nprot_dbase=(Int4)834;
sprintf(fname1,"%s%s",CFG_local_datafilepath,"database.seq");
sprintf(fname2,"%s%s",CFG_local_datafilepath,"database.obs");
/* to here */
  nprot = 1;
  
/*
 * Memory allocations
 */

  seq = Cmatrix(1,nprot_dbase,1,MAXRES);
  if (seq==NULL)
	return ERR_FAIL;
  SEQ = Cmatrix(1,nprot,1,MAXRES);
  if (SEQ==NULL)
	return ERR_FAIL;
  obs = Cmatrix(1,nprot_dbase,1,MAXRES);
  if (obs==NULL)
	return ERR_FAIL;
  title_obs = Cmatrix(1,nprot_dbase,1,MAXLINE);
  if (title_obs==NULL)
	return ERR_FAIL;
  title_seq = Cmatrix(1,nprot_dbase,1,MAXLINE);
  if (title_seq==NULL)
	return ERR_FAIL;
  TITLE = Cmatrix(1,nprot,1,MAXLINE);
  if (TITLE==NULL)
	return ERR_FAIL;
  temp = I4Vector(1,nprot_dbase);
  nres = I4Vector(1,nprot_dbase);
  NRES = I4Vector(1,nprot);
  predi = (CharPtr) cvector(1,MAXRES);
  cPredOrig = (CharPtr) cvector(1,MAXRES);
  probai = FLMatrix(1,MAXRES,1,3);

  for (i=1;i<=(Int4) StringLen(szInputSeq);i++)
  {
		SEQ[1][i] = szInputSeq[i-1];
  }
  NRES[1] = (Int4) StringLen(szInputSeq);
  TITLE[1] = StringCpy(TITLE[1], "Input Sequence");

/*
 * Input the sequences and observed secondary structures for the data base
 */


  if (read_file(fname1,nprot_dbase,seq,title_seq,temp)!=ERR_SUCCESS)
	return ERR_FAIL;


  if (read_file(fname2,nprot_dbase,obs,title_obs,nres)!=ERR_SUCCESS)
	return ERR_FAIL;

/*
 * Check that the data are consistent in the two files
 */

  nerr = 0;
  for(i = 1; i <= nprot_dbase; i++)
    if(temp[i] != nres[i]) {
      ErrPostEx(SEV_ERROR, 1,2,"Error in %ldth protein database consistancy.  Number of residues in one file= %ld, in other file= %ld\n%s\n%s\n\n",(long) i,
       (long) temp[i], (long) nres[i],title_seq[i],title_obs[i]);
      nerr++;
      return ERR_FAIL;
    }

  for(i = 1; i <= nprot_dbase; i++)
    if(StringNCmp(title_seq[i],title_obs[i],100) != 0) {
      ErrPostEx(SEV_ERROR, 1,3,"Error in %ldth protein database consistancy.\n %s \n %s \n",(long) i,title_seq[i],title_obs[i]);
      nerr++;
      return ERR_FAIL;
    }

  if(nerr > 0) {
    ErrPostEx(SEV_ERROR, 1,4,"Database had %ld errors\n",(long) nerr); 
    return ERR_FAIL;    
  }
  
 
/*
 * Calculate the parameters
 */
  if (GOR_Parameters(nprot_dbase,nres,obs,seq)!=ERR_SUCCESS) {
	return ERR_FAIL;
  }
  
/*
 * Predict the secondary structure of protein pro.
 */

  for(pro = 1; pro <= nprot; pro++) {

/*
 * Carry out the prediction for the sequence alone
 */

    predic(NRES[pro],SEQ[pro],predi,probai);
    for (iPred=1;iPred<=NRES[pro];iPred++)
	cPredOrig[iPred]=predi[iPred];
    First_Pass(NRES[pro],probai,predi);
    Second_Pass(NRES[pro],probai,predi);
    for (iPred=1;iPred<=NRES[pro];iPred++) {
	if (cPredOrig[iPred]!=predi[iPred]) {
                wA = obs_indx(cPredOrig[iPred]);
                wB = obs_indx(predi[iPred]);
		wC = 6-wA-wB;
		fDval = probai[iPred][wA]-probai[iPred][wB];
		/* B = B+D = A */
		fBval = probai[iPred][wA];
		/* A=A-D*(C/(A+C)) */                   
		fAval = probai[iPred][wA] - fDval*(probai[iPred][wC]/(probai[iPred][wA]  + probai[iPred][wC]));
		/* C=C-D*(A/(A+C)) */
		fCval = probai[iPred][wC] - fDval*(probai[iPred][wA]/(probai[iPred][wA] + probai[iPred][wC]));
		/* if C<0, then { A=A+C; C=0 } (note we're adding a negative number to A) */
		if (fCval < 0)
		{
		        fAval = fAval + fCval;
		        fCval = 0;               
		}
                probai[iPred][wA] = fAval;
		probai[iPred][wB] = fBval;
		probai[iPred][wC] = fCval;
        }    
    } 

/*
 * Print the results for the protein
 */


	for (iw2=0;iw2<NRES[1];iw2++) {
		arwFiltStruct[iw2] = predi[iw2+1];
		arwProbMatrix[iw2][0] = probai[iw2+1][1];
  		arwProbMatrix[iw2][1] = probai[iw2+1][2];
  		arwProbMatrix[iw2][2] = probai[iw2+1][3];
	}
  }

/*
 * Free memory
 */

  CmatrixFree(seq,1,nprot_dbase,1,MAXRES);
  CmatrixFree(SEQ,1,nprot,1,MAXRES);
  CmatrixFree(obs,1,nprot_dbase,1,MAXRES);
  CmatrixFree(title_obs,1,nprot_dbase,1,MAXLINE);
  CmatrixFree(title_seq,1,nprot_dbase,1,MAXLINE);
  CmatrixFree(TITLE,1,nprot,1,MAXLINE);
  I4VectorFree(temp,1);
  I4VectorFree(nres,1);
  I4VectorFree(NRES,1);
  free_cvector(predi,1,MAXRES);
  free_cvector(cPredOrig,1,MAXRES);
  FLMatrixFree(probai,1,1);

  return ERR_SUCCESS;
}




/*  
$Log: gil_gor.c,v $
Revision 1.12  2003/11/07 18:11:31  feldman
Added file path to database file names

Revision 1.11  2002/01/09 20:17:08  michel
Changed EGORDJKAdded name to EGORDoubleJacknife and added extra function parameter for detailed jk scores

Revision 1.10  2002/01/03 23:02:54  michel
fixed compiler warning

Revision 1.9  2001/12/21 22:45:58  michel
Added (Secondary structure) Enhanced GOR prediction with double jacknifing

Revision 1.8  2001/12/07 22:52:39  michel
added egor_function

Revision 1.7  2001/04/04 21:26:01  feldman
-Removed my BSPrintf, now using the one in SLRIlib - thus needed to add
	this library to makefiles
-changed all fgets to FileGets

Revision 1.6  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.5  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.4  2000/08/14 19:24:50  feldman
Corrected several Windows compatibility issues (all minor)
Improved aligntraj for homology modelling - now works correctly and only adds
  noise when absolutely necessary
Added AMBER potential functions

Revision 1.3  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.2  2000/06/10 16:15:05  feldman
Uodated bzip library with modifications and
corrected Makefiles

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

