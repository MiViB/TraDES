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


  AUTHORS:  Kevin Snyder (ksnyder@blueprint.org)
            Howard Feldman (hfeldman@blueprint.org)
            and Christopher W.V. Hogue  (chogue@blueprint.org)


***************************************************************************
*/

#include "nuschncpy.h"


BiomolDescrPtr CpyBiomolDescrList(BiomolDescrPtr bdp) {
	BiomolDescrPtr bdpNew = NULL, bdpNewHead = NULL, bdpLast = NULL;

	if(bdp == NULL) {
		return bdpNewHead;
	}

	while(bdp != NULL) {
		bdpNew = (BiomolDescrPtr)AsnIoMemCopy(bdp, (AsnReadFunc)BiomolDescrAsnRead, (AsnWriteFunc)BiomolDescrAsnWrite);
		if(bdpNewHead == NULL) {
			bdpNewHead = bdpNew;
			bdpLast = bdpNew;
		}
		else {
			bdpLast->next = bdpNew;
			bdpLast = bdpNew;
		}
		bdp = bdp->next;
	}

	return bdpNewHead;
}

SeqIdPtr CpySeqIdList(SeqIdPtr sip) {
	SeqIdPtr sipNew = NULL, sipNewHead = NULL, sipLast = NULL;

	if(!sip) {
		return sipNewHead;
	}

	while(sip) {
		sipNew = (SeqIdPtr)AsnIoMemCopy(sip, (AsnReadFunc)SeqIdAsnRead, (AsnWriteFunc)SeqIdAsnWrite);
		if(sipNewHead == NULL) {
			sipNewHead = sipNew;
			sipLast = sipNew;
		}
		else {
			sipLast->next = sipNew;
			sipLast = sipNew;
		}
		sip = sip->next;
	}

	return sipNewHead;
}

PALD CpyPALD(PALD pald, PMAD pmadNew, PVNAL pvnalNew) {
	PALD paldNew = NULL;
	FloatLoPtr pflvData = NULL;
	Int2 x;

	paldNew = NewALD();
	
	paldNew->pfbParent = (PFB)pmadNew;
	if(pald->next) {
		paldNew->next = CpyPALD(pald->next, pmadNew, pvnalNew);
	}
	paldNew->pvnalLink = pvnalNew;
	paldNew->iUniqueId = pald->iUniqueId;
	paldNew->cAltConf = pald->cAltConf;
	paldNew->iCoordSet = pald->iCoordSet;
	paldNew->iFloatNo = pald->iFloatNo;
	
	pflvData = FLVector(0,((Int4)(paldNew->iFloatNo)));					

	/* Copy values  */
	/* flatland bug fix < to <= in copy matrix numbers */
	for(x = 0; x <= (Int2)(paldNew->iFloatNo); x++) {
		pflvData[x] = pald->pflvData[x];
	}
	paldNew->pflvData = pflvData;
	paldNew->pGraphic = pald->pGraphic;

	return paldNew;
}

PVNAL CpyPVNALList(PVNAL pvnal, PMAD pmadNew) {
	PVNAL pvnalNew = NULL, pvnalNewHead = NULL, pvnalLast = NULL;
	PALD pald, paldNew = NULL;

	/* One atom location for each model */
	while(pvnal) {
		pald = (PALD)(pvnal->data.ptrvalue);
		pvnalNew = ValNodeNew(NULL);
		pvnalNew->choice = pvnal->choice;

		paldNew = CpyPALD(pald, pmadNew, pvnalNew);
		
		pvnalNew->data.ptrvalue = (VoidPtr)paldNew;
	
		if(pvnalNewHead == NULL) {
			pvnalNewHead = pvnalNew;
			pvnalLast = pvnalNew;
		}	
		else {
			pvnalLast->next = pvnalNew;
			pvnalLast = pvnalNew;
		}

		pvnal = pvnal->next;
	}

	return pvnalNewHead;
}

PVNMA CpyPVNMAList(PVNMA pvnma, PMGD pmgdNew) {
	PVNMA pvnmaNew = NULL, pvnmaNewHead = NULL, pvnmaLast = NULL;
	PMAD pmad, pmadNew = NULL;

	while(pvnma) {
		pmad = (PMAD)(pvnma->data.ptrvalue);
		pvnmaNew = ValNodeNew(NULL);
		pvnmaNew->choice = pvnma->choice;
		pmadNew = NewMAD();

		pmadNew->pfbParent = (PFB)pmgdNew;
		pmadNew->bWhat = pmad->bWhat;
		pmadNew->pvnmaLink = pvnmaNew;
		pmadNew->pcAName = pmad->pcAName;
		
		pmadNew->iIndex = pmad->iIndex;
		pmadNew->pvnalLocate = CpyPVNALList(pmad->pvnalLocate, pmadNew);

		pvnmaNew->data.ptrvalue = (VoidPtr)pmadNew;

		if(pvnmaNewHead == NULL) {
			pvnmaNewHead = pvnmaNew;
			pvnmaLast = pvnmaNew;
		}
		else {
			pvnmaLast->next = pvnmaNew;
			pvnmaLast = pvnmaNew;
		}

		pvnma = pvnma->next;
	}

	return pvnmaNewHead;
}

PDNMG CpyPDNMGList(PDNMG pdnmg, PMMD pmmdNew) {
	PDNMG pdnmgNew = NULL, pdnmgNewHead = NULL, pdnmgLast = NULL;
	PMGD pmgdNew = NULL, pmgdParent = NULL;

	while(pdnmg) {
		pmgdParent = (PMGD)(pdnmg->data.ptrvalue);
		pdnmgNew = DValNodeNew(NULL);
		pdnmgNew->choice = pdnmg->choice;
		pmgdNew = NewMGD();

		pmgdNew->pfbParent = (PFB)pmmdNew;
		pmgdNew->bWhat = pmgdParent->bWhat;
		pmgdNew->bNCBISecStru = pmgdParent->bNCBISecStru;
		pmgdNew->bPDBSecStru = pmgdParent->bPDBSecStru;
		pmgdNew->pdnmgLink = pdnmgNew;
		pmgdNew->iDomain = pmgdParent->iDomain;
		if(pmgdParent->pcGraphName) {
			pmgdNew->pcGraphName = StringSave(pmgdParent->pcGraphName);
		}
		if(pmgdParent->pcGraphNum) {
			pmgdNew->pcGraphNum = StringSave(pmgdParent->pcGraphNum);
		}
		pmgdNew->iIDict = pmgdParent->iIDict;
		if(pmgdParent->pcIUPAC) {
			pmgdNew->pcIUPAC = StringSave(pmgdParent->pcIUPAC);	
		}
		pmgdNew->pcNCBISS = pmgdParent->pcNCBISS;
		pmgdNew->pcPDBSS = pmgdParent->pcPDBSS;
		
		if(pmgdParent->pcPDBComment) {
			pmgdNew->pcPDBComment = StringSave(pmgdParent->pcPDBComment);
		}
		pmgdNew->iAtomCount = pmgdParent->iAtomCount;
		pmgdNew->iBondCount = pmgdParent->iBondCount;
		pmgdNew->pvnmaAHead = CpyPVNMAList(pmgdParent->pvnmaAHead, pmgdNew);

		pdnmgNew->data.ptrvalue = (VoidPtr)pmgdNew;

		if(pdnmgNewHead == NULL) {
			pdnmgNewHead = pdnmgNew;
			pdnmgLast = pdnmgNew;
		}
		else {
			pdnmgLast->next = pdnmgNew;
			pdnmgNew->last = pdnmgLast;
			pdnmgLast = pdnmgNew;
		}
		pdnmg = pdnmg->next;
	}

	return pdnmgNewHead;
}

PMMD CpyPMMD(PMMD pmmd) {
	PMMD pmmdNew;
	PMSD pmsd;

	pmsd = (PMSD)(pmmd->pfbParent);

	pmmdNew = NewMMD();
	pmmdNew->bWhat = pmmd->bWhat;
        pmmdNew->pfbParent = (PFB)pmsd;
        if(pmmd->pcSeqId) {
                pmmdNew->pcSeqId = StringSave(pmmd->pcSeqId);
        }
        if(pmmd->pMolDescr) {
                pmmdNew->pMolDescr = CpyBiomolDescrList(pmmd->pMolDescr);
        }
        if(pmmd->pSeqId) {
                pmmdNew->pSeqId = CpySeqIdList(pmmd->pSeqId);
        }
        pmmdNew->iResCount = pmmd->iResCount;
        pmmdNew->iIRBCount = pmmd->iIRBCount;
        pmmdNew->iGi = pmmd->iGi;
        if(pmmd->pcSequence) {
                pmmdNew->pcSequence = StringSave(pmmd->pcSequence);
        }
	pmmdNew->pdnmgHead = CpyPDNMGList(pmmd->pdnmgHead, pmmdNew);

	return pmmdNew;
}

PDNMG GetPDNMGFromIndex(PDNMG pdnmg, Int2 index) {
	PDNMG pdnmgHere = NULL;

	pdnmgHere = pdnmg;

	while(pdnmgHere) {
		if(pdnmgHere->choice == index) {
			return pdnmgHere;
		}
		pdnmgHere = pdnmgHere->next;
	}

	return NULL;
}

PMAD GetPMADFromName(PVNMA pvnmaA, CharPtr pcAName) {
	PVNMA pvnmaHere = NULL;
	PMAD pmad;

	pvnmaHere = pvnmaA;
	while(pvnmaHere) {
		pmad = (PMAD)(pvnmaHere->data.ptrvalue);
		if(!StringCmp(pmad->pcAName, pcAName)) {
			return pmad;
		}

		pvnmaHere = pvnmaHere->next;
	}

	return NULL;
}

PVNMB CpyPVNMBList(PVNMB pvnmb, PDNMG pdnmgNewHead, Pointer parentPtr) {
	PDNMG pdnmg, pdnmgNew, pdnmgPmadFrom = NULL, pdnmgPmadTo = NULL;
	PMGD pmgd, pmgdNew, pmgdTmp = NULL;
	PVNMB pvnmbNew, pvnmbNewHead = NULL, pvnmbLast = NULL;
	PMBD pmbd, pmbdNew;

	while(pvnmb) {
		pmbd = (PMBD)(pvnmb->data.ptrvalue);	
		pvnmbNew = ValNodeNew(NULL);
		pvnmbNew->choice = pvnmb->choice;

		pmbdNew = NewMBD();

		pmbdNew->pfbParent = (PFB)parentPtr;
		pmbdNew->pvnmbLink = pvnmbNew;
		pmbdNew->bWhat = pmbd->bWhat;
		
		pmgdTmp = (PMGD)(pmbd->pmadFrom->pfbParent);
		pdnmgPmadFrom = GetPDNMGFromIndex(pdnmgNewHead, ((PDNMG)(pmgdTmp->pdnmgLink))->choice);
		pmbdNew->pmadFrom = GetPMADFromName(((PMGD)(pdnmgPmadFrom->data.ptrvalue))->pvnmaAHead, pmbd->pmadFrom->pcAName);

		pmgdTmp = (PMGD)(pmbd->pmadTo->pfbParent);
		pdnmgPmadTo = GetPDNMGFromIndex(pdnmgNewHead, ((PDNMG)(pmgdTmp->pdnmgLink))->choice);
		pmbdNew->pmadTo = GetPMADFromName(((PMGD)(pdnmgPmadTo->data.ptrvalue))->pvnmaAHead, pmbd->pmadTo->pcAName);

		/* Add PMBD to the appropriate PMADs */
		ValNodeAddPointer(&pmbdNew->pmadFrom->pvnBonds, 0, (VoidPtr)pmbdNew);			
		ValNodeAddPointer(&pmbdNew->pmadTo->pvnBonds, 0, (VoidPtr)pmbdNew);			
		
		pvnmbNew->data.ptrvalue = (VoidPtr)pmbdNew;

		if(pvnmbNewHead == NULL) {
                        pvnmbNewHead = pvnmbNew;
                        pvnmbLast = pvnmbNew;
                }
                else {
                        pvnmbLast->next = pvnmbNew;
                        pvnmbLast = pvnmbNew;
                }

		pvnmb = pvnmb->next;
	}

	return pvnmbNewHead;
}

void CpyBonds(PMMD pmmd, PMMD pmmdNew) {
	PDNMG pdnmg, pdnmgNew;
	PMGD pmgd, pmgdNew;
	PVNMB pvnmb;

	/* Copy Graph bonds */
	pdnmg = pmmd->pdnmgHead;
	pdnmgNew = pmmdNew->pdnmgHead;
	
	while(pdnmg) {
		pmgd = (PMGD)(pdnmg->data.ptrvalue);
		pmgdNew = (PMGD)(pdnmgNew->data.ptrvalue);

		pmgdNew->pvnmbBHead = CpyPVNMBList(pmgd->pvnmbBHead, pmmdNew->pdnmgHead, pmgdNew); 

		pdnmg = pdnmg->next;
		pdnmgNew = pdnmgNew->next;
	}

	/* Copy inter-residue bonds */
	pmmdNew->pvnmbIRBHead = CpyPVNMBList(pmmd->pvnmbIRBHead, pmmdNew->pdnmgHead, pmmdNew);
}

void AddPALDsToVector(PointerPtr ppAsnOrder, PDNMM pdnmmHead) {
	PDNMM pdnmm;
	PMMD pmmd;
	PDNMG pdnmg;
        PMGD pmgd;
        PVNMA pvnma;
        PMAD pmad;
        PVNAL pvnal;
	PALD pald;
	Int4 paldCount = 1;
	
	pdnmm = pdnmmHead;
	while(pdnmm) {
		pmmd = (PMMD)(pdnmm->data.ptrvalue);
                pdnmg = pmmd->pdnmgHead;
                while(pdnmg) {
                        pmgd = (PMGD)(pdnmg->data.ptrvalue);
                        pvnma = pmgd->pvnmaAHead;
                        while(pvnma) {
                                pmad = (PMAD)(pvnma->data.ptrvalue);
                                pvnal = pmad->pvnalLocate;
                                while(pvnal) {
					pald = (PALD)(pvnal->data.ptrvalue);
					ppAsnOrder[paldCount] = (Pointer) pald;
					paldCount++;

					/* Add locations for occupancy */
					while(pald->next) {
						pald = pald->next;
						ppAsnOrder[paldCount] = (Pointer) pald;
						paldCount++;
					}
					pvnal = pvnal->next;
				}
				pvnma = pvnma->next;
			}
			pdnmg = pdnmg->next;
		}
		pdnmm = pdnmm->next;
	}

	/* Set the last element to NULL */
	ppAsnOrder[paldCount] = NULL;
}



BiostrucFeaturePtr CpyPBSF(BiostrucFeaturePtr pbsfParent, Int4 newMolNum, Int4 lastFeatId, ResidueIntervalPntrPtr pri) {
	BiostrucFeaturePtr pbsfNew;
	ResidueIntervalPntrPtr priNew;
	ValNodePtr vnp, vnp2, vnp3;

	pbsfNew = BiostrucFeatureNew();
	pbsfNew->id = lastFeatId + 1;
	pbsfNew->name = StringSave(pbsfParent->name);
	pbsfNew->type = pbsfParent->type;
	
	priNew = ResidueIntervalPntrNew();
	priNew->molecule_id = newMolNum;
	priNew->from = pri->from;
	priNew->to = pri->to;

	vnp3=ValNodeNew(NULL);
        vnp3->choice=ResiduePntrs_interval;
        vnp3->data.ptrvalue=(VoidPtr)priNew;
        vnp2=ValNodeNew(NULL);
        vnp2->choice=ChemGraphPntrs_residues;
        vnp2->data.ptrvalue=(VoidPtr)vnp3;
        vnp=ValNodeNew(NULL); 
        vnp->choice=Location_location_subgraph;
        vnp->data.ptrvalue=(VoidPtr)vnp2;

	pbsfNew->Location_location=vnp;

	return pbsfNew;
}

void AddNewChainFeatures(PMSD pmsd, Int4 parentMolNum, Int4 newMolNum) {
	BiostrucFeatureSetPtr pbsfs;
	BiostrucFeaturePtr pbsfLast, pbsfHere, pbsfNew;
	ValNodePtr vnp, vnp2, vnp3;
	ResidueIntervalPntrPtr pri;

	pbsfs = pmsd->pbsBS->features;
	/* Iterate through each feature type */
	while(pbsfs) {
		pbsfHere= pbsfs->features;
		/* Grab the pointer to the last feature */
		while(pbsfHere) {
			if(!pbsfHere->next) {
				pbsfLast = pbsfHere;
			}
			pbsfHere = pbsfHere->next;
		}

		pbsfHere = pbsfs->features;
		while(pbsfHere) {
			vnp = pbsfHere->Location_location;
			if(vnp->choice == Location_location_subgraph) {
				vnp2 = (ValNodePtr)(vnp->data.ptrvalue);
				if(vnp2->choice == ChemGraphPntrs_residues) {
					vnp3 = (ValNodePtr)(vnp2->data.ptrvalue);
					if(vnp3->choice == ResiduePntrs_interval) {
						pri = (ResidueIntervalPntrPtr)(vnp3->data.ptrvalue);
						if(pri->molecule_id == parentMolNum) {
							pbsfNew = CpyPBSF(pbsfHere,newMolNum,pbsfLast->id, pri);
							pbsfLast->next = pbsfNew;
							pbsfLast = pbsfNew;
						}		
					}
				}
			}
			pbsfHere = pbsfHere->next;
		}

		pbsfs = pbsfs->next;
	}
}


CharPtr Misc_TrimSpacesAroundString ( CharPtr str ){
      char *ptr, *dst, *revPtr;
      int spaceCounter = 0;
  
      ptr = dst = revPtr = str;
      
      if ( !str || str[0] == '\0' )
          return str;
  
      while ( *revPtr != '\0' )
          if ( *revPtr++ <= ' ' )
              spaceCounter++;
      
      if ( (revPtr - str) <= spaceCounter )
      {
         *str = '\0';
          return str;
      }
  
      while ( revPtr > str && *revPtr <= ' ' ) revPtr--;
      while ( ptr < revPtr && *ptr <= ' ' ) ptr++;
      while ( ptr <= revPtr ) *dst++ = *ptr++;
      *dst = '\0';
      
      return str;
}


CharPtr NextUniqueChainName(PMSD pmsd) {
/* Must return a new letter for the added chain name */
	PDNMM pdnmmHere;
	PMMD pmmdHere;
	char *ptr,*chainptr;
	CharPtr pcChainName=NULL;
	CharPtr newChainName=NULL;
	char Alphabet[]="123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	Int2 i=0;	
	
	chainptr=Alphabet;
	pdnmmHere=pmsd->pdnmmHead;
	
	while(pdnmmHere != NULL){
		
		/*printf("pdnmmHere is %d\n",pdnmmHere->choice);*/
		pmmdHere=(PMMD)pdnmmHere->data.ptrvalue;
		if (pmmdHere==NULL) break;
		pcChainName=pmmdHere->pcMolName;
		if (pcChainName==NULL) break;
		Misc_TrimSpacesAroundString(pcChainName);
		/* printf("pcChainame here is %s\n",pcChainName); */
		if (StringLen(pcChainName)==1 ) {
			ptr=Alphabet;
			while (*ptr != '\0'){
				if (*ptr == pcChainName[0]) {
					
					*ptr = '0';  /* marks the ones in use with 0 */
					break;
				}
				ptr++;
			}	
		
		}
		
		pdnmmHere=pdnmmHere->next;
	}
	
	newChainName=(CharPtr)MemNew(2*sizeof(Char));
	
	while(*chainptr !='\0'){
		if(*chainptr != '0') {
			newChainName[0] = *chainptr;
			break;
		}
		chainptr++;
	}
	return newChainName;
/* See original code ...  beware of non biopolymer names */

}


Boolean CopyNewChain(PDNMM pdnmmFrom, PDNMM pdnmmTo,PMSD pmsdTo) {
	PDNMM pdnmmNew = NULL;
	PMMD pmmdNew, pmmdTo, pmmdFrom;
	PDNMG pdnmg;
	PMGD pmgd;
	PVNMA pvnma;
	PMAD pmad;
	PVNAL pvnal;
	PALD pald;
	FloatLo delta = 0.001;
	ValNodePtr vnp;
	SeqIdPtr sip;
	PDBSeqIdPtr psp;
	
	pmmdTo = (PMMD)(pdnmmTo->data.ptrvalue);
	/*pmsdTo = (PMSD)pmmdTo->pfbParent;*/

	pmmdFrom = (PMMD)(pdnmmFrom->data.ptrvalue);

	(pmsdTo->iMolCount)++;

	/* Create a new molecule */
	pdnmmNew = DValNodeNew(NULL);
	pdnmmNew->choice = pdnmmTo->choice + 1;
	pdnmmTo->next = pdnmmNew;
	pdnmmNew->last = pdnmmTo;

	pmmdNew = CpyPMMD(pmmdFrom);

	/* Alter values in PMMD to make new chain unique */
	pmmdNew->iChainId = pmmdTo->iChainId + 1;
	
	if (pmmdNew->bWhat & AM_PROT)
		pmmdNew->pcMolName = NextUniqueChainName(pmsdTo); 

	pmmdNew->pdnmmLink = pdnmmNew;
	
	/*CONSIDER THIS CODE..TEST ITS EFFECT ON GENERATED ASN.1 ... */
	vnp = pmmdNew->pMolDescr;
	/* Change name in molecule description */
	while(vnp) {
		if(vnp->choice == BiomolDescr_name) {
			MemFree(vnp->data.ptrvalue);
			vnp->data.ptrvalue = (VoidPtr)(StringSave(pmmdNew->pcMolName));
			break;
		}
		vnp = vnp->next;
	}
	
	/* Need to go back and fill in bond info */
	CpyBonds(pmmdFrom, pmmdNew);

	pdnmmNew->data.ptrvalue = (VoidPtr)pmmdNew;

	/* Copy features of parent molecule to new molecule */
	/* this will require separate pmsdFrom, pmsdTo parameters */
	/*****OR  JUST COPY CHAIN AND DO NO ADD ITS FEATURES ? */
/*	AddNewChainFeatures(pmsdTo, pmmdFrom->iChainId, pmmdNew->iChainId);	*/
	
	
	
	/* change the added PDB id code to be unique from RAND to RNDx */
	
    sip = pmmdNew->pSeqId;	
	if (sip != NULL)
	if (sip->choice == SEQID_PDB)
	{
		psp = sip->data.ptrvalue;
		if (!StringCmp(psp->mol,"RAND")) {
		/*	printf("RAND CASE ALTERING IN CHEMICAL GRAPH\n"); */
			psp->mol[0] = '\0';
			StringCpy(psp->mol,"RND");
			StringCat(psp->mol,pmmdNew->pcMolName);
			psp->chain = (Uint1) *(pmmdNew->pcMolName);
			psp->rel = DateCurr();	
			
			
		}
		
	}
		
	
	 	
	return TRUE;
}


Boolean CopyBiomolecule(PMSD pmsdTo, PDNMM pdnmmThisMolecule ) {

	PDNMM pdnmmHeadTo = NULL, pdnmmLastTo = NULL, pdnmmParentTo = NULL, pdnmmHereTo;
	PMMD pmmdTo = NULL, pmmdHereTo = NULL,pmmdChain=NULL;
	PMLD pmldTo = NULL;
	PDNML pdnmlTo = NULL;

	Boolean addedChain = FALSE;
	SeqIdPtr sip;
	PDBSeqIdPtr psp;
	ValNodePtr vnp = NULL;
	CharPtr chnid = NULL;

/*	ErrPostEx(SEV_WARNING,0,0,"Copy chain %s To structure %ld",  NextUniqueChainName(pmsdTo), pmsdTo-iMMDBid);
*/
	pdnmmHeadTo = pmsdTo->pdnmmHead;
	pmmdTo = (PMMD)(pdnmmHeadTo->data.ptrvalue);
/*For trades generated structures, give every chain a chain name*/
			
    /* Copy the molecule pdnmmThisMolecule in its entirety to pmsdTo */
	
	/* Grab the last molecule in the TO list.  This is where new chains will be appended */
	pdnmmLastTo = pdnmmHeadTo;
	while(pdnmmLastTo) {
		pmmdChain=(PMMD)pdnmmLastTo->data.ptrvalue;
                /*pmmdChain->pcMolName = NextUniqueChainName(pmsdTo);*/
/*		printf("pmmdchain->pcMolname here is a%da\n",pmmdChain->pcMolName[0]); */
		
		
		/* change the original PDB id codes to be unique from RAND to RNDx for multi-chain Cn3D lookups*/
		
		sip = pmmdChain->pSeqId;	
		if (sip != NULL)
		if (sip->choice == SEQID_PDB)
		{
			psp = sip->data.ptrvalue;
			if (!StringCmp(psp->mol,"RAND")) {
				/*	printf("RAND CASE ALTERING IN CHEMICAL GRAPH\n"); */
				if (pmmdChain->pcMolName[0]==(Uint1) '\0'){
					pmmdChain->pcMolName[0]= (Uint1) 'A';
				/*	printf("changed to a\n"); */
				}
				psp->mol[0] = '\0';
				StringCpy(psp->mol,"RND");
				psp->mol[3]= (Uint1) *(pmmdChain->pcMolName);
				psp->mol[4] = '\0';
				psp->chain = (Uint1) *(pmmdChain->pcMolName);
				psp->rel = DateCurr();	
				
				/* must also change the chain code in the ASN.1 stub */
				
				if (vnp = ValNodeFindNext(pmmdChain->pMolDescr, NULL, BiomolDescr_name))
					chnid = vnp->data.ptrvalue;
				chnid[0] = (Uint1) pmmdChain->pcMolName[0];
								
			}
			
		}
		
		
		if(!pdnmmLastTo->next) {
			break;
		}
		pdnmmLastTo = pdnmmLastTo->next;
	}	
	
/*	printf("Original Changes Phase Done \n"); */

	/* Add chain to modelstruc */
	
    if(CopyNewChain(pdnmmThisMolecule, pdnmmLastTo,pmsdTo)) {
				addedChain = TRUE;
				pdnmmLastTo = pdnmmLastTo->next;
			}			
	if(addedChain) {
			ErrPostEx(SEV_WARNING,0,0,"Created new chain");
		}
/* this message appears in ribosome processing bug as last warning before it goes off endless loop*/
	
	
	/*88888888888888888888888888888888888888888*/	

	/* free-err double free-ing perhaps?  - or a dictionary bug 
	   modelstruc loaders wipe out ppAsnOrder if the are too large in memory...? still true? 
	   if so this fixes the freeing so it is pointer safe */

	if(addedChain) {
		/* Usage of ppAsnOrder appears to be optional.  For purposes of speed, it is freed instead of updated */
		pdnmlTo = pmsdTo->pdnmlModels;
		if (pdnmlTo) {
		      pmldTo = (PMLD)(pdnmlTo->data.ptrvalue);
			  if (pmldTo)
		         if(pmldTo->ppAsnOrder) {
			       PTRVectorFree(pmldTo->ppAsnOrder,0);
			       pmldTo->ppAsnOrder = NULL;
				 } 		
		}
	}



	return TRUE;
}
/*
$Log: nuschncpy.c,v $
Revision 1.3  2011/01/04 09:47:00  mingxi
Commit by mingxi

Revision 1.2  2009/02/16 17:53:00  chogue
File conversion utility for TraDES compatibility with PDB file output

Revision 1.1  2009/01/27 17:47:23  chogue
initial revision of nuschncpy derived from GenerateBiounit

 

*/
