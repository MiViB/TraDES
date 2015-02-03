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


/* utilities for dealing with crease energy */
#include <mmdbtraj.h>

static DValNodePtr pdnZhangAtmList = NULL;
static Int4Ptr piBryantPotential = NULL;
static Int4Ptr piZhangPtnl = NULL;
static FloatLoPtr CreaseEnergy = NULL;
static Int4 NumCreaseResidues = 0;
static FloatLo TotalCreaseEnergy = 0.0;

TrajErr LoadPotential(Int2 Units)
{
	FILE *fpzh,*fpAtmInfo,*fpbl;
	Char buf[PATH_MAX];

	if (Units==UNITS_ZHANG) {
			sprintf(buf,"%s%s",CFG_local_datafilepath,ZHANG_POTENTIAL);
	        if ((fpzh=FileOpen(buf,"r"))==NULL) {
	                ErrPostEx(SEV_ERROR,1,7,"File not found %s",buf);
	                return ERR_FAIL;
	         }
			sprintf(buf,"%s%s",CFG_local_datafilepath,ZHANG_ATOMS);
	         if ((fpAtmInfo=FileOpen(buf,"r"))==NULL) {
	                ErrPostEx(SEV_ERROR,1,7,"File not found %s",buf);
			FileClose(fpzh);
	                return ERR_FAIL;
	        }
		/* load potentials */
	        if (LoadZhangPotential (&piZhangPtnl, fpzh)!=ERR_SUCCESS) {
	                ErrPostEx(SEV_ERROR,2,7,"Error in LoadZhangPotential");
			FileClose(fpzh);
			FileClose(fpAtmInfo);
	                return ERR_FAIL;
		}
	        LoadZhangAtmList (&pdnZhangAtmList, fpAtmInfo);
	        FileClose (fpzh);
	        FileClose (fpAtmInfo);
	}
	if (Units==UNITS_BRYANT) {
			sprintf(buf,"%s%s",CFG_local_datafilepath,BL_POTENTIAL);
	        if ((fpbl=FileOpen(buf,"r"))==NULL) {
	                ErrPostEx(SEV_ERROR,1,8,"File not found %s",buf);
	                return ERR_FAIL;
	        }
	        if (LoadBryantPotential (&piBryantPotential, fpbl, TRUE /*bUsingPep*/, FALSE /*bNegateInput*/)!=ERR_SUCCESS) {
	                ErrPostEx(SEV_ERROR,2,8,"Error in LoadBryantPotential");
			FileClose(fpbl);
	                return ERR_FAIL;
		}
	        FileClose (fpbl);
	}
	return ERR_SUCCESS;
}

void UnLoadPotentials(void)
{
        FreeZhangAtmList (&pdnZhangAtmList);
        FreeZhangPotential (&piZhangPtnl);
        FreeBryantPotential (&piBryantPotential);
}

/* units specify the potential to be used here */
/* operates on first molecule */
/* cutoff specifies furthest distance between residues in contact to consider, so
   that interacting N and C termini will not affect the whole protein for example */
TrajErr CalcCreaseEnergy(PMMD pmmd,Int2 Units,Int2 cutoff,Int2 ExcludeWindow,Boolean Decay,Boolean Color,Int2 ModelNum)
{
	DValNodePtr pdnListPmmds=NULL;
	DValNodePtr pdnHead = NULL,pdnRes,pdnInteract;
	PMGD pmgd2;
	PMAD pmad;
	PALD pald;
	PALN palnHere;
    	Int2 choice=0;
	Int4 numres,cnt,res1,res2,midresx2,distx2;
	FloatHi pot,potHere;

	if (ExcludeWindow>=cutoff) {
		ErrPostEx(SEV_ERROR,1,15,"Window should be smaller than cutoff to CalcCreaseEnergy");
		return ERR_FAIL;
	}
	if (pmmd==NULL) return ERR_FAIL;
	if (LoadPotential(Units)!=ERR_SUCCESS) return ERR_FAIL;
	choice=(pmmd->pdnmmLink)->choice;
	DValNodeAddPointer(&pdnListPmmds,choice,pmmd);
	numres=pmmd->iResCount;
	NumCreaseResidues=numres;
	CreaseEnergy=(FloatLoPtr)MemNew((numres+1)*sizeof(FloatLo));
	for (cnt=0;cnt<=numres;cnt++)
		CreaseEnergy[cnt]=0.0;
	if (Units==UNITS_ZHANG)
        	ComputeZhangPotential (&pdnHead, pdnListPmmds, ModelNum, TRUE/*bInclusiveWindow*/,
        		cutoff/*iWindowSize*/, piZhangPtnl, pdnZhangAtmList, TRUE/*bDetailedList*/);
	else if (Units==UNITS_BRYANT) {
		/* prevent dangling pointer problem by saving any previous
		   Pep MGDs */
		BackupPepMGDs();
	        ComputeBryantPotential (&pdnHead, pdnListPmmds, ModelNum, TRUE/*bInclusiveWindow*/,
	        	cutoff/*iWindowSize*/, piBryantPotential, TRUE/*bUsingPep*/, TRUE/*bDetailedList*/);
	}
	else {
	        FreeListPmmds (&pdnListPmmds);
		UnLoadPotentials();
		return ERR_FAIL;
	}
	pdnRes=pdnHead;
	while (pdnRes!=NULL) {
		res1=pdnRes->choice;
		pdnInteract=(DValNodePtr)(pdnRes->data.ptrvalue);
		if (pdnInteract==NULL) {
			UnLoadPotentials();
			return ERR_FAIL;
		}
		/* skip first node which is just summary info */
		pdnInteract=pdnInteract->next;
		while (pdnInteract!=NULL) {
			palnHere=(PALN)(pdnInteract->data.ptrvalue);
			pmgd2=palnHere->pmgd;
			res2=(pmgd2->pdnmgLink)->choice;
			/* add contribution from each residue */
			pot=(palnHere->potential)*2.0;
			/* interaction between res1 and res2 of value pot */
			/* only treat each interaction once */
			if (res2>res1) {
				if (abs(res2-res1)>cutoff) {
					ErrPostEx(SEV_ERROR,1,0,"Internal inconsistency during crease energy calculation %d %d %d",res1,res2,cutoff);
					UnLoadPotentials();
					return ERR_FAIL;
				}
				/* skip i to i+window size */
				if (abs(res2-res1)>ExcludeWindow) {
					/* middle residue number * 2 */
					midresx2=res1+res2;
					for (cnt=res1;cnt<=res2;cnt++) {
						distx2=abs(2*cnt-midresx2);
						/* dist is # residues from current residue to middle of
						   interaction, in units of 0.5 residues */
						/* adjust potential by distance */
						potHere=pot;
						if (Decay) {
/*							potHere=potHere/((FloatLo)cnt*(FloatLo)(NumCreaseResidues+1-cnt));*/
							potHere=2.0*potHere/((FloatLo)distx2+1.0);
						}
						CreaseEnergy[cnt]+=potHere;
					}
				}
			}
			pdnInteract=pdnInteract->next;
		}	
		pdnRes=pdnRes->next;
	}
	if (Color==TRUE) {
		pdnRes=pdnHead;
		/* set B factor to crease energy */
		while (pdnRes!=NULL) {
			res1=pdnRes->choice;
			pdnInteract=(DValNodePtr)(pdnRes->data.ptrvalue);
			palnHere=(PALN)(pdnInteract->data.ptrvalue);
			pmgd2=palnHere->pmgd;
	       		if ((pmad = FindCAlpha (pmgd2->pvnmaAHead)) != NULL) {
	            		if (((pald = GetAtomLocs (pmad, ModelNum)) != NULL) && (pald->iFloatNo >= 4))
	                		pald->pflvData[4] = (FloatLo) (CreaseEnergy[res1] * (-1.0)) + 15.0;
	 		}
			pdnRes=pdnRes->next;
		}
	}
        FreeAdjList (&pdnHead);
        FreeListPmmds (&pdnListPmmds);
	UnLoadPotentials();
	if (Units==UNITS_BRYANT) {
		/* restore any previously saved Pep MGDs */
		RestorePepMGDs();
	}
	TotalCreaseEnergy=0.0;
	for (cnt=1;cnt<=NumCreaseResidues;cnt++) {
		TotalCreaseEnergy+=CreaseEnergy[cnt];
	}
	return ERR_SUCCESS;
}

void FreeCreaseEnergy(void)
{
	MemFree(CreaseEnergy);
	CreaseEnergy=NULL;
}

FloatLo GetCreaseEnergy(void)
{
	return TotalCreaseEnergy;
}

FloatLo GetCreaseRes(Int4 res)
{
	if (CreaseEnergy==NULL) return 0.0;
	if (res<1 || res>NumCreaseResidues) return 0.0;
	return CreaseEnergy[res];
}

void PrintCreaseEnergy(void)
{
	Int4 cnt;

	if (CreaseEnergy==NULL)
		return;
	printf("Residue\tCrease Energy\n-------\t-------------\n");
	for (cnt=1;cnt<=NumCreaseResidues;cnt++)
		printf("  %3ld\t%f\n",(long int)cnt,CreaseEnergy[cnt]);
	printf("\nTotal Crease Energy: %f\n",GetCreaseEnergy());
	return;
}



/*  
$Log: crease.c,v $
Revision 1.9  2003/03/07 19:59:54  feldman
Fixed ambiguous expression

Revision 1.8  2002/08/11 15:41:23  feldman
Improved error message

Revision 1.7  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.6  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.5  2000/10/24 20:57:25  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.4  2000/07/13 20:41:06  feldman
Added inclusive/exclusive windows for crease energy and corrected
makeTG wrapping for Phi-Psi space

Revision 1.3  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.2  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

