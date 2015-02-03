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
 
#include "clustlib.h"
#include <slri_misc.h>
#include <seqhound.h>

#ifdef _ENTREZAPI_
#include <accentr.h>  /* EntrezSeqEntryGet in WriteFASTAfromGI */
                      /* EntrezLinkUidList in NetValidateCSAN */
#include <seqport.h>  /* SeqPortPtr in CompareSeqWithAlign for NetValidateCSAN */
#include <objacces.h> /* For LinkSetPtr in NetValidateCSAN*/
#include <tofasta.h>  /* BioseqtoFasta in WriteFASTAfromGI*/
#endif

#define BIGGEROF(bach1,bach2) (((bach1) > (bach2)) ? (bach1) : (bach2))


/* Turns on code (ie printf) for debuging */
/*#define DEBUGCLUST*/
static Int2 *Blosum=NULL;
static Char blocol[PATH_MAX];

PCSAN LIBCALL NewCSAN(void) 
{
	PCSAN pcsanThis = NULL;   
	pcsanThis = (PCSAN) MemNew((size_t)sizeof(CSAN));
	return pcsanThis;
}

void FreeCSAN(PCSAN pcsanThis)
{
	if (pcsanThis->next) 
		FreeCSAN(pcsanThis->next);
		pcsanThis->next = NULL;
	if (pcsanThis->pcSeqName) 
		MemFree(pcsanThis->pcSeqName);
	if (pcsanThis->pcSeqAln) 
		MemFree(pcsanThis->pcSeqAln);
	if (pcsanThis->pExtra)
		MemFree(pcsanThis->pExtra);
	if (pcsanThis->pvnSeqAlnPart) 
		ValNodeFreeData(pcsanThis->pvnSeqAlnPart);
	if (pcsanThis) {
		MemFree(pcsanThis);
	}
/*  if (pcsanThis->pdnmsStructure)
  	FreeAModelstruc(pcsanThis->pdnmsStructure);
****pdnmsStructure is to be freed with mmdbapi's ClearStructure (Done in main part of program) *****/
}
 

Int4 LIBCALL WriteCSAN(PCSAN pcsanThis, FILE *pFile){
  Int4 iLine = 0;
  Int4 iLineCt = 0;
  Int4 iLineIn = 0;
  Int4 iPad = 0;
  PCSAN pcsanTemp = NULL;
  Char cHold = '-';
  CharPtr pcTemp = NULL;
  CharPtr pcStart = NULL;
  CharPtr pcEnd = NULL;
  Boolean bSameLen = TRUE;
  ValNodePtr pvnSeqAlnTemp;
  Int4    iTotal;  
  Int4    iLen, iSlen; 


  if (!pcsanThis) return 0;
  if (!pFile) return 0;
    
  fprintf(pFile,"CLUSTAL X   SLRI I/O Module\n\n\n");
  pcsanTemp = pcsanThis;
  iLen = pcsanTemp->iLen;
  while (pcsanTemp)  /* Find longest string */
    {

#ifdef DEBUGCLUST
      printf("%s[%s]\n\n",  pcsanTemp->pcSeqName,  pcsanTemp->pcSeqAln);
#endif

      /* chop up the sequence into a valnode linked list of bits to print out */
      if (iLen < pcsanTemp->iLen) {
	iLen = pcsanTemp->iLen;  /* find length of longest sequence */
	bSameLen = FALSE;
      }
      pcsanTemp = pcsanTemp->next;
    }
  pcsanTemp = pcsanThis;
  if (!bSameLen)
    
    /* normalize string lengths if necessary */
    while (pcsanTemp) {
      /* reallocate all to match the largest iLen */
      pcTemp = (CharPtr) MemNew((size_t) ((iLen * sizeof(char)) + 2));
      StringCpy(pcTemp, pcsanTemp->pcSeqAln);
      MemFree(pcsanTemp->pcSeqAln);
      pcsanTemp->pcSeqAln = pcTemp;
      iPad = iLen - StringLen(pcTemp);
      if (iPad) {
	StringCat(pcTemp,"-");
      }
      pcsanTemp = pcsanTemp->next;
    }
  
   iLine = 0;
   iLineCt = iLen / 60;
   if (iLen % 60) iLineCt++;

#ifdef DEBUGCLUST
   printf("Line[%ld], LineCt[%ld]\n", (long) iLine, (long) iLineCt); 
#endif   
   
   pcsanTemp = pcsanThis;
   while (pcsanTemp)  /* Chop and allocate pieces into linked list */
     {
       pcStart = pcsanTemp->pcSeqAln;
       pcEnd = pcStart;
       iTotal = 0;
       while (iTotal < iLen)
	 {
          iSlen=0;
          do {
	    iTotal++;
	    iSlen++;
	    pcEnd++;
	    
#ifdef DEBUGCLUST
	    printf("(%ld)",(long) iSlen); 
#endif
	    
	  } while ((*pcEnd != '\0') && (iSlen < 60));
          
          if (*pcEnd == '\0') 
	    { 
	      pcTemp = StringSave(pcStart);
              ValNodeAddStr(&pcsanTemp->pvnSeqAlnPart, 0,  pcTemp);
	      break;
	    }
	  else
	    {
	      cHold = *pcEnd;
	      *pcEnd = '\0';
	      pcTemp = StringSave(pcStart);
	      ValNodeAddStr(&pcsanTemp->pvnSeqAlnPart, 0,  pcTemp);
	      *pcEnd = cHold;
	      pcStart = pcEnd;
	    }
      	}
       pcsanTemp = pcsanTemp->next;
     }
   iLine = 0;
   iLineCt = iLen / 60;
   if (iLen % 60) iLineCt++;
   
   while (iLine < iLineCt)
     {
       pcsanTemp = pcsanThis;
       while (pcsanTemp)
	 {
	   pvnSeqAlnTemp = pcsanTemp->pvnSeqAlnPart;
	   iLineIn = 0;
	   while (iLineIn < iLine)
	     {
	       pvnSeqAlnTemp = pvnSeqAlnTemp->next;
	       iLineIn++;
	     }
	   fprintf(pFile, "%-23.23s   %s\n",pcsanTemp->pcSeqName,(CharPtr) pvnSeqAlnTemp->data.ptrvalue);
	   fflush(pFile);
	   pcsanTemp = pcsanTemp->next;
	 }
       iLine++;
       fprintf(pFile, "\n");
     }
   fflush(pFile);
   return iLen;
}

/* Takes a clustal file and parses it */
PCSAN LIBCALL ReadCSAN(CharPtr fnam)
{
  FILE *pFile; 
  BiostrucPtr bspBiostruc;
  PDNMS pdnms;
  PCSAN pcsanLast = NULL;
  PCSAN pcsanHead = NULL;
  PCSAN pcsanTemp = NULL;
  Char pcBuf[250];
  Char pcTmpFnam[PATH_MAX];
  CharPtr ppcParse[125];
  CharPtr pcTemp = NULL;
  CharPtr pcTest = NULL;
  CharPtr pcSeqPart = NULL;
  ValNodePtr pvn = NULL;
  Int4 i = 0;
  Int4 iString = 0;
#ifdef DEBUGCLUST
  Int4 iLen =0;
#endif
  Int4 iWC = 0;
  Int4 iWL = 0;
  Int4 iGi = 0;
  long int lGi=0;
  Int4 firstGI = 0;
  Int4 iNumSeq = 0;
  Int4 iAliLen = 0;
  Int2 flag =0;
	AsnIoPtr aip;

  if (!fnam) 
    return(NULL);
  
  if ((pFile = FileOpen(fnam, "r"))==NULL)
    return(NULL);
  
  /* MUST SAVE pcSeqName portion of the string in the node...*/   

  do   /* get each line, tokenize into ppcParse */
    {
      pcBuf[0] = '\0';
      pcTest = FileGets(pcBuf,(size_t)250,pFile);

#ifdef DEBUGCLUST
      printf("[%s]\n",pcTest);
#endif

      if (pcTest) {
#ifdef DEBUGCLUST
	iLen = StringLen(pcTest);
	/* Break into words delimited by whitespace or "|" */
	/* gi|123456| */
	
	  printf("LineLength %d\n",(int) iLen);
#endif
	
	  /* makes array 0 */
	  for (i=0; i < 125; i++)
	    ppcParse[i] = NULL;
	
	  pcTemp = pcTest;
	  ppcParse[0] = pcTest;
	  iWC = 1;
	  while ((*pcTemp != '\0') && (iWC < 125)) {
	    if ((*pcTemp == '|') || (*pcTemp == ' ') || (*pcTemp == '\t') || (*pcTemp == '\n')) {
	      *pcTemp = '\0';  /* terminate current word */
	      do {
		pcTemp++;  /* skip whitespace to next word */
	      }
	
	      while ((*pcTemp == '|') || (*pcTemp == ' ') || (*pcTemp == '\t') || (*pcTemp == '\n'));
	      ppcParse[iWC] = pcTemp;  /* start next word */
	      iWC++;
	    }
	    else
	      pcTemp++;
	  }
	
	
#ifdef DEBUGCLUST	
	  for (i=0; i<125; i++)
	    if (ppcParse[i] != NULL)
	      printf("%ld/%ld [%s]\n",i,iWC, ppcParse[i]);	
#endif
	
	  /* If ppcParse[0] is "gi" then this is a VALID sequence line and the last
	     word parsed is the sequence */
	  if (ppcParse[0][0] != '\0') {
	
	
	    /* WHACK THIS TO SCAN FOR GI-NUMBER PAIR */
	
	
	    iWL = StringLen(ppcParse[0]);
	    if (iWL == 2)
	      if (((ppcParse[0][0] == 'g') ||
		   (ppcParse[0][0] == 'G')) &&
		  ((ppcParse[0][1] == 'i') ||
		   (ppcParse[0][1] == 'I')))
		{
		  /* next word should be GI number */
		  iGi = 0;

		  if (sscanf(ppcParse[1],"%ld",&lGi)) {
			iGi=(Int4)lGi;
		    /*printf("GI is (%ld)\n",iGi); */
		
		    if (flag == 0){
		      flag = 1;
		      firstGI = iGi;
		    }
		
		    /* find the 2nd last string = sequence */
		    /* last string is stripped off line terminator, by the way... */
		
		    /* start after GI */
		    iString = 2;
		
		    /* iString is the LAST sting - line terminator */
		    while (ppcParse[iString+1] != NULL){
		      iString++;
		    }
		
		    pcsanTemp = pcsanHead;
		
		    /* look for existing pcsan */
		    while (pcsanTemp){
		
		      /* WHACK match by GI or by identical leader string */
		      if (iGi == pcsanTemp->iGi){
			break;
			}
		      pcsanTemp = pcsanTemp->next;
		    }
		
		    /* otherwise add a new one */
		    if (!pcsanTemp) {
		      pcsanTemp = NewCSAN();
		      pcsanTemp->next=NULL;
		      if (pcsanHead==NULL) {
			pcsanHead=pcsanTemp;
			pcsanLast=pcsanTemp;
		      }
		      else {
			pcsanLast->next=pcsanTemp;
			pcsanLast=pcsanTemp;
		      }
		      pcsanTemp->iGi = iGi;
		
		      /* Sequence name for unknown sequence */
		      if (iGi == firstGI)			
			pcsanTemp->pcSeqName = StringSave(ppcParse[2]);
		      else
			  /* Checks if there is a pdb file name available */
			if (!((ppcParse[2][0] == 'p') ||
			      (ppcParse[2][0] == 'P')) &&
			    ((ppcParse[2][1] == 'd') ||
			     (ppcParse[2][1] == 'D')) &&
			    ((ppcParse[2][2] == 'b') ||
			     (ppcParse[2][2] == 'B')))
			  {
			    ErrPostEx(SEV_ERROR,1,1, "Requires a pdb file");
			    return NULL;
			  }
		
		      /* Seq name for sequence with a specified chain */
		      if ((iWC == 7) && (iGi != firstGI)){
			pcsanTemp->pcSeqName = (CharPtr)MemNew(sizeof(Char)*(1+StringLen(ppcParse[3])+StringLen(ppcParse[4])));
			StringCpy(pcsanTemp->pcSeqName,ppcParse[3]);
			StringCat(pcsanTemp->pcSeqName,ppcParse[4]);
		      }			
		
		      else  			
			  /* Seq name for sequence without a specified chain */
			if ((iWC == 6) && (iGi != firstGI)){
			  pcsanTemp->pcSeqName = StringSave(ppcParse[3]);
			  /*pcsanTemp->pcSeqName = ppcParse[3];*/
			}
		    }
		
		    /* Initialize fields */
		    pcSeqPart = NULL;
		    pcSeqPart = StringSave(ppcParse[iString-1]);
		
#ifdef DEBUGCLUST
		    printf("SEQ [%ld]= [%s]\n",(long) iString-1, pcSeqPart);
#endif
		
		    ValNodeAddStr(&pcsanTemp->pvnSeqAlnPart, 0,  pcSeqPart);
		  }  /* if GI */
		}  /* if line starts with GI */
	  }  /* if line starts with a word */
      } /* if pcTest */
    } while (pcTest);


#ifdef DEBUGCLUST
  printf("REPORT\n");
  pcsanTemp = pcsanHead;

  while (pcsanTemp) {
    printf("GI = [%ld]\n", (long) pcsanTemp->iGi);
    printf("Sequence = [%s]\n", pcsanTemp->pcSeqName);
    pvn = pcsanTemp->pvnSeqAlnPart;

    while (pvn) {
      printf(">%s<\n",(char *) pvn->data.ptrvalue);
      pvn = pvn->next;
    }
    pcsanTemp = pcsanTemp->next;
  }
#endif


  FileClose(pFile);

  /*CONSOLIDATE THE PIECES*/

  /*count the lengths of each one*/
  pcsanTemp = pcsanHead;
  while (pcsanTemp) {
    iNumSeq++;     	
    pvn = pcsanTemp->pvnSeqAlnPart;

    /* Count lengths - same? */
    while (pvn) {
      pcsanTemp->iLen += strlen((char *) pvn->data.ptrvalue);
      pvn = pvn->next;
    }
    pcsanTemp = pcsanTemp->next;
  }

  /*compare - same sequence lengths ?*/
  pcsanTemp = pcsanHead;
  iAliLen = pcsanTemp->iLen;


#ifdef DEBUGCLUST
  printf("Alignment Length %ld\n",(long) iAliLen);
#endif


  while (pcsanTemp) {
    if (pcsanTemp->iLen != iAliLen) {
      ErrPostEx(SEV_ERROR,1,1,"Validation Error: Sequence Alignment Lengths Differ %d vs %d\n", pcsanTemp->iLen, iAliLen);
      return NULL;
    }                                   	
    pcsanTemp = pcsanTemp->next;
  }

  /*make room for full sequences*/	
  pcsanTemp = pcsanHead;
  while (pcsanTemp) {
    pcsanTemp->pcSeqAln = (CharPtr) MemNew((size_t) iAliLen+1 * sizeof(char));
    if (!pcsanTemp->pcSeqAln){
      ErrPostEx(SEV_ERROR,0,0,"Memory Alloc Error\n");
      FreeCSAN(pcsanTemp);
      return NULL;
    }
    pcsanTemp = pcsanTemp->next;
  }

  /*concatenate them together*/
  pcsanTemp = pcsanHead;
  while (pcsanTemp) {
    pvn = pcsanTemp->pvnSeqAlnPart;
    StringCpy(pcsanTemp->pcSeqAln,(char *) pvn->data.ptrvalue);
    pvn = pvn->next;
    while (pvn) {
      StringCat(pcsanTemp->pcSeqAln, (char *) pvn->data.ptrvalue);
      pvn = pvn->next;
    }
    pcsanTemp = pcsanTemp->next;
  }

  /*lose the fragments*/
  pcsanTemp = pcsanHead;
  while (pcsanTemp) {
    ValNodeFreeData(pcsanTemp->pvnSeqAlnPart);
    pcsanTemp->pvnSeqAlnPart = NULL;
    pcsanTemp = pcsanTemp->next;
  }

  /* Load pdnms structure */
  pcsanTemp = pcsanHead;
  while (pcsanTemp->next){
    pcsanTemp = pcsanTemp -> next ;
    StringNCpy(pcTmpFnam, pcsanTemp -> pcSeqName,4);
    pcTmpFnam[4] = '\0';
		StringCat(pcTmpFnam,MMDB_EXT);
    ErrPostEx(SEV_INFO,0,0,"Loading %s\n",pcTmpFnam);
    bspBiostruc = NULL;
	aip=AsnIoOpen(pcTmpFnam,"rb");
    bspBiostruc = BiostrucAsnRead(aip,NULL);
	AsnIoClose(aip);
/*    FetchBiostrucPDB(pcTmpFnam, 0,100);*/
    if (bspBiostruc == NULL){
      ErrPostEx(SEV_ERROR,3,1,"Unable to fetch Biostruc %s", pcTmpFnam);
      return NULL;
    }


    pdnms = MakeAModelstruc (bspBiostruc);

    if (pdnms == NULL){
      ErrPostEx(SEV_ERROR,4,1, "Unable to convert Biostruc to Modelstruc");
      return NULL;
    }
    pcsanTemp -> iMMDB = ((PMSD)(pdnms->data.ptrvalue))->iMMDBid;
    pcsanTemp -> pdnmsStructure = pdnms;
  }
  return pcsanHead;

/* parserr:
  FileClose(pFile);
  FreeCSAN(pcsanHead);
  ErrPostEx(SEV_ERROR,0,0, "Error in input format \n");
  return NULL;*/
}

#ifdef _ENTREZAPI_
void WriteFASTAfromGI(Int4 iGiSeq, FILE *pOut)
{
    SeqEntryPtr  pSeqEntry = NULL;
    BioseqPtr pBioseq = NULL;
    Int2 retcode = 0;
    SeqIdPtr pSeqId = NULL;


    retcode = 0;
    pSeqEntry = EntrezSeqEntryGet(iGiSeq, retcode);
    if (pSeqEntry == NULL) {
       ErrPostEx(SEV_ERROR,0,0,"No sequence for GI [%ld]\n", (long) iGiSeq);
	   return;
	}
    pSeqId = ValNodeNew(NULL);
    pSeqId->choice = SEQID_GI;
    pSeqId->data.intvalue = iGiSeq;
    pBioseq = BioseqFindInSeqEntry(pSeqId, pSeqEntry);

    if (BioseqToFasta(pBioseq, pOut, FALSE) != TRUE) {
        ErrPostEx(SEV_ERROR,0,0,"BioseqToFasta Failed on GI=%ld\n",(long)iGiSeq);  
        return;
    }    
    SeqEntryFree(pSeqEntry);
    MemFree(pSeqId);

}
#endif

PMMD MMDBGetMoleculeByGI(PDNMS pdnmsThis, Int4 iGi)
{
  PMSD pmsdThis = NULL;
  PDNMM pdnmmThis = NULL;
  PMMD pmmdThis = NULL;

  if (!pdnmsThis) return NULL;
  if (!iGi) return NULL;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
  pdnmmThis = pmsdThis->pdnmmHead;
  while (pdnmmThis) {
      pmmdThis = (PMMD) pdnmmThis->data.ptrvalue;
      if (iGi == pmmdThis->iGi) return pmmdThis;
      pdnmmThis = pdnmmThis->next;
  }
  return NULL;
}

PDNMG GetResidueListFromPCSAN(PCSAN pcsan)
{
  PDNMM pdnmmHere;
  PMMD pmmdHere=NULL;
  Char chainid;
	
	if (pcsan==NULL) {
		ErrPostEx(SEV_FATAL,1,22,"Inconsistency check failure, contact authors");
		return NULL;
	}
	pdnmmHere=((PMSD)((pcsan->pdnmsStructure)->data.ptrvalue))->pdnmmHead;
		if (StringLen(pcsan->pcSeqName)==4)
		  pmmdHere=(PMMD)(pdnmmHere->data.ptrvalue);
		else {
		  chainid=(pcsan->pcSeqName)[4];
		  while (pdnmmHere!=NULL) {
		    pmmdHere=(PMMD)(pdnmmHere->data.ptrvalue);
		    if (chainid==(pmmdHere->pcMolName)[0])
		      break;
		    pdnmmHere=pdnmmHere->next;
		  }				
		  if (pdnmmHere==NULL)
		    ErrPostEx(SEV_FATAL,1,21,"Unable to find chain in molecule");
		}
		/* now pmmdHere points to the right molecule */
		return pmmdHere->pdnmgHead;
}

/* Fills CSAN structre with the structure's sequence */
Int4  FillCSANWithStru(PCSAN pcsanThis, PMMD pmmdThis, Int4 iLen)
{
     Int4 iCount = 0;
     CharPtr pcA;
     PDNMG pdnmgThis = NULL;
     PMGD pmgdThis = NULL;

     if (!pcsanThis) return 0;
     if (!pmmdThis) return 0;
     if (!IsProtein(pmmdThis)) return 0;
     if (!iLen) return 0;
 
     pcsanThis->pcSeqAln = (CharPtr)MemNew((size_t) (1+ sizeof(char) * iLen));
     pcA = pcsanThis->pcSeqAln;
     pdnmgThis = pmmdThis->pdnmgHead;
     iCount = 0;
     while ((pdnmgThis != NULL)  && (iCount < iLen))
      {
            iCount++;
            pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
            *pcA = (char) pmgdThis->pcIUPAC[0] ;  
            pcA++;
            pdnmgThis = pdnmgThis->next;
      }
     while (iCount < iLen)
       {
            *pcA = '-';
            pcA++;  
            iCount++;
       }
    pcsanThis->pcSeqAln[iLen] = '\0';  
    return iCount;
}


/* Fills the CSAN structure with secondary structure information */
Int4 FillCSANWithMask(PCSAN pcsanThis, PMMD pmmdThis, Int4 iLen)
{
     Int4 iCount = 0;
     CharPtr pcA;
     PDNMG pdnmgThis = NULL;
     PMGD pmgdThis = NULL;
     Char cCode = ' ';

     
     if (!pcsanThis) return 0;
     if (!pmmdThis) return 0;
     if (!IsProtein(pmmdThis)) return 0;
     if (!iLen) return 0;
 
 
     pcsanThis->pcSeqAln = (CharPtr)MemNew((size_t) (1+ sizeof(char) * iLen));
     pcA = pcsanThis->pcSeqAln;
     pdnmgThis = pmmdThis->pdnmgHead;
     iCount = 0;
     while ((pdnmgThis != NULL)  && (iCount < iLen))
      {
           iCount++;
           /* Traverse graph data and get secondary structure information */
		   /* A for alpha helix, B for beta sheet */
		   /* If Author SS and VAST SS definitions agree, capitalize */

		   pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
           cCode = '.';  /* default loop character */
           if (pmgdThis->bPDBSecStru & (Byte) SS_HELIX)  /*author described */
              cCode = 'a';
           if (pmgdThis->bNCBISecStru & (Byte) SS_HELIX)  /* VAST ss definition */
              cCode = 'a';
           if ((pmgdThis->bPDBSecStru & (Byte) SS_HELIX) && (pmgdThis->bNCBISecStru & (Byte) SS_HELIX))
              cCode = 'A';  /* consensus !*/
           if (pmgdThis->bPDBSecStru & (Byte) SS_STRAND)  /* author described */
              cCode = 'b';
           if  (pmgdThis->bNCBISecStru & (Byte) SS_STRAND)  /* VAST ss definition */
              cCode = 'b';
           if ((pmgdThis->bPDBSecStru & (Byte) SS_STRAND) && (pmgdThis->bNCBISecStru & (Byte) SS_STRAND))
              cCode = 'B';  /* consensus */
            *pcA = cCode ;  
            pcA++;
            pdnmgThis = pdnmgThis->next;
      }
     while (iCount < iLen)  /* pad to the end if this is short */
       {
            *pcA = '-';
            pcA++;  
            iCount++;
       }
    pcsanThis->pcSeqAln[iLen] = '\0';  
    return iCount;
}


/*  Fills CSAN with the sequence from the bioseq, padded to length iLen*/
/*  WARNING not called and not tested... */

#ifdef _ENTREZAPI_
Int4 FillCSANWithSeq(PCSAN pcsanThis, BioseqPtr pbsq, Int4 iLen)
{
     SeqPortPtr spp = NULL;
     Uint1 code = Seq_code_ncbieaa;
     Uint1 residue;
     Int4 iCount = 0;
     CharPtr pcA;
     
     if (!pcsanThis) return 0;
     if (!ISA_aa(pbsq->mol)) return 0;
     if (!iLen) return 0;
     
     spp = SeqPortNew(pbsq, 0, -1, 0, code);
     if (!spp) return 0;
     SeqPortSeek(spp, 0, SEEK_SET);
     pcsanThis->pcSeqAln = (CharPtr)MemNew((size_t) (1+ sizeof(char) * iLen));
     pcA = pcsanThis->pcSeqAln;
     residue = SeqPortGetResidue(spp);
     iCount = 0;
     while ((residue != SEQPORT_EOF) && (residue != '\0') && (iLen < iCount))
      {
            iCount++;
            *pcA = (char) residue;  
            pcA++;
            residue = SeqPortGetResidue(spp);
      }
     while (iCount < iLen)
       {
            *pcA = '-';
            pcA++;  
            iCount++;
      }
    pcsanThis->pcSeqAln[iLen+1] = '\0';  
    SeqPortFree(spp);
    return iCount;
}
#endif


PCSAN CSANStructMask(Int4 iMMDBid, Int4 iGiStr)
{
    BiostrucPtr pbsThis = NULL;
    PDNMS pdnmsStructure = NULL;
    PMSD  pmsdThis = NULL;
    PMMD  pmmdMol = NULL;
    PCSAN pcsanStruc = NULL;
    PCSAN pcsanMask = NULL;
	Char pcName[PATH_MAX];

    if (iGiStr == 0 || iMMDBid == 0) return NULL;
	
    /* get a pdnms to the structure MMDBid, then a pdnmm to the embedded Gi */
	pbsThis = SHoundGet3DEx(iMMDBid,  ALLMDL, 100);
	if(pbsThis == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Unable to get biostruc for mmdbid %ld.\n",(long) iMMDBid);
		return NULL;
	}
	if((pdnmsStructure = MakeAModelstruc(pbsThis)) == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Fatal Error during MMDB structure retrieval %ld\n", (long) iMMDBid);
		BiostrucFree(pbsThis);
		return NULL;
	}
     
    if ((pmmdMol = MMDBGetMoleculeByGI(pdnmsStructure, iGiStr)) == NULL) { 
		ErrPostEx(SEV_ERROR,0,0,"Fatal Error in Molecule retrieval\n");
		return NULL;
	}
    
	pmsdThis = (PMSD) pdnmsStructure->data.ptrvalue;
    
    /* make node for Mask */
    pcsanMask = NewCSAN();
    pcsanMask->iGi = iGiStr;
	pcsanMask->iMMDB = iMMDBid;
    pcsanMask->iLen = pmmdMol->iResCount;
	pcsanMask->next = NewCSAN();
    
	/* Make node for Structure */
	pcsanStruc= pcsanMask->next;
    pcsanStruc->iGi = iGiStr;
	pcsanStruc->iMMDB = iMMDBid;
    pcsanStruc->iLen = pmmdMol->iResCount;
	pcsanStruc->pdnmsStructure = pdnmsStructure;

	/* for each pcsan, copy ID headers into SeqName field */
    pcName[0] = '\0';
    sprintf(pcName,"!SS_gi|%ld|pdb|%4s|%1s",(long) iGiStr, pmsdThis->pcPDBName, pmmdMol->pcMolName);
    pcsanMask->pcSeqName = StringSave(pcName);
    if (StringLen(pcName) > 23) pcsanMask->pcSeqName[23] = '\0';
    
	sprintf(pcName,"gi|%ld|pdb|%4s|%1s",(long) iGiStr, pmsdThis->pcPDBName, pmmdMol->pcMolName);
    pcsanStruc->pcSeqName = StringSave(pcName);
    if (StringLen(pcName) > 23) pcsanMask->pcSeqName[23] = '\0';
    
    /* fill the sec-structure mask with the structure, pad the remainder */ 
    FillCSANWithMask(pcsanMask, pmmdMol, pmmdMol->iResCount);

    /* fill the CSAN node with the structure, pad the remainder if shorter */
    FillCSANWithStru(pcsanStruc, pmmdMol, pmmdMol->iResCount);

    return pcsanMask;

}


Int2 SSMasktoDisk(Int4 iMMDBid, Int4 iGiStr)
{
	FILE *pOut;
	Char path[PATH_MAX];
	PCSAN pcsanMask = NULL;

	pcsanMask =  CSANStructMask(iMMDBid, iGiStr);

	/* write out clustal file */
	sprintf(path,"%ld.aln",(long) iGiStr);
	if ((pOut = FileOpen(path, "w")) == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Cannot open output file, aborting: %s\n", path);
		return 0;
	}
	WriteCSAN(pcsanMask, pOut); 
	fflush(pOut);
	FileClose(pOut);
	FreeCSAN(pcsanMask);

	ErrPostEx(SEV_INFO,0,0,"Done\n");
	return 0;
}


#ifdef _ENTREZAPI_
static Boolean CompareSeqWithAlign(BioseqPtr pbsq,  CharPtr pcAlign)
{
     SeqPortPtr spp = NULL;
     Uint1 code = Seq_code_ncbieaa;
     Uint1 residue;
     Int4 ctr = 0;
     CharPtr pcA;
     
     if (!pcAlign) return FALSE;
     if (!pbsq) return FALSE;
     if (!ISA_aa(pbsq->mol)) return FALSE;
     
     spp = SeqPortNew(pbsq, 0, -1, 0, code);
     if (!spp) return FALSE;
     SeqPortSeek(spp, 0, SEEK_SET);
     pcA = StringUpper(pcAlign);
     residue = SeqPortGetResidue(spp);
     while ((residue != SEQPORT_EOF) && (residue != '\0') && (*pcA != '\0'))
      {
          /*  printf("%c = ",(char) residue); */
            while ((*pcA == '-'))
              {
                pcA++;  
              }
            if (*pcA == '\0') goto out;
          /*  printf("%c\n",*pcA); */
            if (*pcA != (char) residue)
              {
            /*   printf("No Match\n"); */   
               return FALSE;
              }
            pcA++;
            residue = SeqPortGetResidue(spp);
      }
out: 
    SeqPortFree(spp);
    return TRUE;
}


PCSAN NetValidateCSAN(PCSAN pcsanThis, Boolean bCheckSeq)
{
    PCSAN pcsanTemp = NULL, pcsanFrom = NULL, pcsanTo = NULL;
    ValNodePtr pvn = NULL;
    Int4 iNumSeq = 0;
    Int4 iAliLen = 0;
    long GetGi = 0;
    LinkSetPtr lsp1 = NULL;
    BiostrucPtr pbsThis = NULL;
    SeqEntryPtr  sep;
    SeqIdPtr pSeqId = NULL;
    BioseqPtr pbsq;
    Int2 retcode = 0;
    Int4 i = 0;
    Int4 l = 0;
    Int4  seqnum = 0;
    Boolean bSeqMatch = FALSE;
  
	if (!pcsanThis) return NULL;

	pcsanTemp = pcsanThis;
	while (pcsanTemp) {
		iNumSeq++;
		pvn = pcsanTemp->pvnSeqAlnPart;
		while (pvn) {
			/* Count lengths - same? */
			pcsanTemp->iLen += strlen((char *) pvn->data.ptrvalue);
			pvn = pvn->next;
		}
		pcsanTemp = pcsanTemp->next;
	}
	/* Pairwise??  */ 
	if (iNumSeq != 2) {
		FreeCSAN(pcsanThis);
		ErrPostEx(SEV_ERROR,0,0,"Validation Error: Not a pairwise alignment file - %ld sequences found\n", (long) iNumSeq);
		return NULL;
	}
	/* Same sequence lengths ?*/
	pcsanTemp = pcsanThis;
	iAliLen = pcsanTemp->iLen;
	ErrPostEx(SEV_INFO,0,0,"Alignment Length %ld\n",(long) iAliLen);
	if (iAliLen != pcsanTemp->next->iLen) {
		ErrPostEx(SEV_ERROR,0,0,"Validation Error: Sequence Alignment Lengths Differ/nGi#1=%ld Len=%ld Gi#2=%ld Len=%ld\n",
		(long) pcsanTemp->iGi, (long) iAliLen, (long) pcsanTemp->next->iGi, (long) pcsanTemp->next->iLen);
		FreeCSAN(pcsanThis);
		return NULL;
	}
	/* allocate complete sequences*/
	pcsanThis->pcSeqAln = (CharPtr) MemNew((size_t) iAliLen+1 * sizeof(char));
	if (!pcsanThis->pcSeqAln) { 
		ErrPostEx(SEV_ERROR,0,0,"Memory Alloc Error\n"); 
		FreeCSAN(pcsanThis); 
		return NULL;
	}
	pcsanThis->next->pcSeqAln = (CharPtr) MemNew((size_t) iAliLen+1 * sizeof(char));
	if (!pcsanThis->next->pcSeqAln) { 
		ErrPostEx(SEV_ERROR,0,0,"Memory Alloc Error\n"); 
		FreeCSAN(pcsanThis); 
		return NULL;
	}
	/* copy over the sequence bits into one item */
	pvn = pcsanThis->pvnSeqAlnPart;
	StringCpy(pcsanThis->pcSeqAln,(char *) pvn->data.ptrvalue);
	pvn = pvn->next;
    while (pvn) {
		StringCat(pcsanThis->pcSeqAln, (char *) pvn->data.ptrvalue);
		pvn = pvn->next;
	}

	pvn = pcsanThis->next->pvnSeqAlnPart;
	StringCpy(pcsanThis->next->pcSeqAln,(char *) pvn->data.ptrvalue);
	pvn = pvn->next;
	while (pvn)	{
		StringCat(pcsanThis->next->pcSeqAln, (char *) pvn->data.ptrvalue);
		pvn = pvn->next;
	}

	/* printf("Seq 1 %ld char [%s]\n",(long) strlen(pcsanThis->pcSeqAln),pcsanThis->pcSeqAln);
	printf("Seq 2 %ld char [%s]\n",(long) strlen(pcsanThis->next->pcSeqAln),pcsanThis->next->pcSeqAln);
	*/
	/* test each gi for structure in MMDB */

	GetGi = pcsanThis->iGi;
	EntrezLinkUidList(&lsp1, TYP_AA, TYP_ST, 1, &GetGi, FALSE);
	if (lsp1 != NULL) {
		if (lsp1->num > 0)	{
			pcsanThis->iMMDB = lsp1->uids[0];
			ErrPostEx(SEV_INFO,0,0,"Sequence [%ld] is the structure [%ld]\n",(long) GetGi, 
			(long) pcsanThis->iMMDB);
			pcsanFrom = pcsanThis;
			pcsanTo = pcsanThis->next;
		}
	}
	GetGi = pcsanThis->next->iGi;
	EntrezLinkUidList(&lsp1, TYP_AA, TYP_ST, 1, &GetGi, FALSE);
	if (lsp1 != NULL) {
		if (lsp1->num > 0) {
			pcsanThis->next->iMMDB = lsp1->uids[0];
			ErrPostEx(SEV_INFO,0,0,"Sequence [%ld] is the structure [%ld]\n",(long) GetGi, 
			(long) pcsanThis->next->iMMDB);
			pcsanFrom = pcsanThis->next;
			pcsanTo = pcsanThis;
		}
	}
	/* should be ONE structure */

	if ((pcsanThis->iMMDB == 0) && (pcsanThis->next->iMMDB == 0)) {
		ErrPostEx(SEV_ERROR,0,0,"Neither sequence matches a structure in MMDB at NCBI by GI \n");
		return NULL;
	}
	if ((pcsanThis->iMMDB != 0) && (pcsanThis->next->iMMDB != 0)) {
		ErrPostEx(SEV_ERROR,0,0,"Ambiguous File: Both sequences match a structure in MMDB at NCBI by GI\n");
		return NULL;
	}
	if (pcsanThis->iMMDB != 0) {
		if ((pbsThis = MMDBBiostrucGet( pcsanThis->iMMDB,  0, 1)) != NULL) {
			pcsanThis->pdnmsStructure = MakeAModelstruc(pbsThis);
		}
		if (pcsanThis->pdnmsStructure) {
			ErrPostEx(SEV_INFO,0,0,"Got structure 1\n"); */
	/* WriteFASTASeqHet(pcsanThis->pdnmsStructure,  stdout); */
		} else {
			ErrPostEx(SEV_ERROR,0,0,"Fatal Error during MMDB structure retrieval %ld\n", (long) pcsanThis->iMMDB);
		}
	}
	if (bCheckSeq) { /* get the sequence */
		retcode = 0;
		sep = EntrezSeqEntryGet(pcsanThis->iGi, retcode);
		if (sep == NULL) {
			ErrPostEx(SEV_ERROR,0,0,"No sequence for GI [%ld]\n", (long) pcsanThis->iGi);
		}
		pSeqId = ValNodeNew(NULL);
		pSeqId->choice = SEQID_GI;
		pSeqId->data.intvalue = pcsanThis->iGi;
		pbsq = BioseqFindInSeqEntry(pSeqId, sep);
		bSeqMatch = CompareSeqWithAlign(pbsq, pcsanThis->pcSeqAln);
		SeqEntryFree(sep);
		ValNodeFree(pSeqId);
		if ((bSeqMatch == FALSE) && (bCheckSeq == TRUE)) {
			ErrPostEx(SEV_ERROR,0,0,"Failed to match with NCBI's version of [%ld]\n",(long) pcsanThis->iGi);
			return NULL; /* fail code here ... */
		}
	}
	if (pcsanThis->next->iMMDB != 0) {
		if ((pbsThis = MMDBBiostrucGet( pcsanThis->next->iMMDB,  0, 1)) != NULL) {
			pcsanThis->next->pdnmsStructure = MakeAModelstruc(pbsThis);
		}
		if (pcsanThis->next->pdnmsStructure)  {
			ErrPostEx(SEV_INFO,0,0,"Got structure 2\n"); 
		/*   WriteFASTASeqHet(pcsanThis->next->pdnmsStructure,  stdout); */
		} else {
			ErrPostEx(SEV_ERROR,0,0,"Fatal Error during MMDB structure retrieval %ld\n", (long) pcsanThis->next->iMMDB);
		}
	}
	if (bCheckSeq) {
	/* get the sequence */
		retcode = 0;
		sep = EntrezSeqEntryGet(pcsanThis->next->iGi, retcode);
		if (sep == NULL) 
			ErrPostEx(SEV_ERROR,0,0,"No sequence for GI [%ld]\n", (long) pcsanThis->next->iGi);
		pSeqId = ValNodeNew(NULL);
		pSeqId->choice = SEQID_GI;
		pSeqId->data.intvalue = pcsanThis->next->iGi;
		pbsq = BioseqFindInSeqEntry(pSeqId, sep);
		bSeqMatch = CompareSeqWithAlign(pbsq, pcsanThis->next->pcSeqAln);
		SeqEntryFree(sep);
		ValNodeFree(pSeqId);
		if ((bSeqMatch == FALSE) && (bCheckSeq == TRUE)) {
			ErrPostEx(SEV_ERROR,0,0,"Failed to match with NCBI's version of [%ld]\n",(long) pcsanThis->iGi);
			return NULL; /* fail code here ... */
		}
	}
	/*  printf("returning pcsanThis\n"); */
	return pcsanThis;
}
#endif


/* 
	GetSeqInfoFromStruc

	Retrieves the  sequence or the secondary structure sequence from a protein structure
	or from a sequence in the bReserved Field (modelling)
*/


void GetSeqInfoFromStruc(PMMD pmmdThis, Uint2 CLUSTLIB_INFO_TYPE, Boolean bUseRealSeq, CharPtr PNTR ppcSeq)
{
	Boolean bFirst = TRUE;
	PDNMG pdnmgThis = NULL;
	PMGD pmgdThis = NULL;
	Int4 iResCt = 0;
	CharPtr pcTemp = NULL;
	Char pcHere[PATH_MAX];

	if ((Int2) pmmdThis->bWhat == AM_PROT) {
		pdnmgThis = pmmdThis->pdnmgHead;
	} else {
		pdnmgThis = NULL;
		return;
	}
	while(pdnmgThis) {
		iResCt++;
		pdnmgThis = pdnmgThis->next;
	}
	if((pcTemp = MemNew((size_t) sizeof(Char) * (iResCt*3)+1)) == NULL) {
		return;
	}


	iResCt = 0;
	pdnmgThis = pmmdThis->pdnmgHead;
	while (pdnmgThis) {
		iResCt++;
		pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
		/* NCBI secondary Structure */
		if (CLUSTLIB_INFO_TYPE == CLUSTLIB_INFO_VASTSS) { 
			if (pmgdThis->bReserved != (Byte) 'X') {
				if (pmgdThis->bNCBISecStru & (Byte) SS_HELIX) {
					StringCat(pcTemp, "H");
				} else if (pmgdThis->bNCBISecStru & (Byte) SS_STRAND) {
					StringCat(pcTemp, "E");
				} else {		
					StringCat(pcTemp, "C"); /* default coil character */
				}
			} else {
				StringCat(pcTemp, "X");
			}
		/* PDB SS */
		} else if (CLUSTLIB_INFO_TYPE == CLUSTLIB_INFO_PDBSS) {
			if (pmgdThis->bReserved != (Byte) 'X') {
				if (pmgdThis->bPDBSecStru & (Byte) SS_HELIX) {
					StringCat(pcTemp, "H");	
				} else if (pmgdThis->bPDBSecStru & (Byte) SS_STRAND) {
					StringCat(pcTemp, "E");
				} else {		
					StringCat(pcTemp, "C");/* default coil character */
				}
			} else {
				StringCat(pcTemp, "X");
			}
		/* For domain info, comma delimited */
		} else if (CLUSTLIB_INFO_TYPE == CLUSTLIB_INFO_DOM) {
			if (pmgdThis->bReserved != (Byte) 'X') {
				if(bFirst == TRUE) {
					sprintf(pcHere,"%d",(int) pmgdThis->iDomain);
					bFirst = FALSE;
				} else {
					sprintf(pcHere,",%d",(int) pmgdThis->iDomain);
				}
				StringCat(pcTemp, pcHere);
			} else {
				StringCat(pcTemp, "X,");
			}	
		/* For sequence */
		} else if(CLUSTLIB_INFO_TYPE == CLUSTLIB_INFO_SEQ) {
			/* prints 1 letter code in the pcIUPAC field */
			if ((bUseRealSeq) && (pmgdThis->bReserved != (Byte) 'X')) {
				if((pmgdThis->pcIUPAC)) StringCat(pcTemp, &pmgdThis->pcIUPAC[0]);
			} else { /* prints the 1 letter code in the bReserved field */
				sprintf(pcHere, "%c", toupper((char) pmgdThis->bReserved));
				StringCat(pcTemp, pcHere);
			}
		} 
		pdnmgThis = pdnmgThis->next;
	}
	
	*ppcSeq = pcTemp;
	return;
}


void DumpSequenceToFile(CharPtr pcHeader, CharPtr pcDefline, CharPtr pcFooter, CharPtr pcSeq, FILE *fp, Int4 columnsize)
{
	Int4 i = 0, len = 0, iResCt = 0;

	if( (pcDefline == NULL) || (pcSeq == NULL) || (fp == NULL) ) {
		ErrPostEx(SEV_ERROR,0,0,"Invalid Parameters.");
		return;
	}

	/* Header */
	fprintf(fp,"%s%s",pcHeader,pcDefline);

	/* Body */
	len = StringLen(pcSeq);
	for(i = 0; i < len; i++) {
		if(iResCt % columnsize == 0) {
			fprintf(fp,"\n");
		}
		fprintf(fp,"%c",pcSeq[i]);
		iResCt++;
	}

	/* Footer */
	fprintf(fp,"%s", pcFooter);
	return;
}

void DumpPKBFASTAFile(CharPtr pcDefline,CharPtr pcSeq, FILE *fp)
{
	Int4 columnsize = 60;
	DumpSequenceToFile(">",pcDefline,"\n", pcSeq,fp,columnsize);
}

void DumpGORFASTAFile(CharPtr pcDefline, CharPtr pcSeq, FILE *fp)
{
	Int4 columnsize = 50;
	DumpSequenceToFile("!",pcDefline,"@\n", pcSeq, fp, columnsize);
}

void DumpPKBFASTA(Int4 iSeqGI, Int4 iStrucGI, CharPtr pcDefline, CharPtr pcSeq, Int2 iType) 
{
	FILE *fp = NULL;
	Char path[PATH_MAX];
	CharPtr end[] = {"seq","str","ss"};
	CharPtr pcend;

	if(iType == CLUSTSEQ) {
		pcend = end[0];
	} else if (iType == CLUSTSTR) {
		pcend = end[1];
	} else if (iType == CLUSTSS) {
		pcend = end[2];
	} else {
		ErrPostEx(SEV_ERROR,0,0,"Invalid Type");
		return;
	}

	sprintf(path,"%ld_%ld.%s",(long) iSeqGI,(long) iStrucGI, pcend); 
	fp = FileOpen(path, "w");
	DumpPKBFASTAFile(pcDefline,pcSeq,fp);
	FileClose(fp);
}

Int2 CompressSequence(CharPtr pcSeq, CharPtr PNTR ppNewSeq)
{
	CharPtr pcSeqHere = NULL;
	Int4 i = 0, iCount = 0, len = 0;
	
	if((len = StringLen(pcSeq)) == 0) return 1;
	if((pcSeqHere = (CharPtr) MemNew((size_t)sizeof(Char)*(len +1))) == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Could not allocate new CharPtr");
		return 1;
	}

	for(i = 0; i < len; i++) {	
		if(pcSeq[i] != 'X') {
			pcSeqHere[iCount] = pcSeq[i];
			iCount++;
		}	
	}
	*ppNewSeq = pcSeqHere;
	return 0;
}

void DumpGORFASTA(Int4 iSeqGI, Int4 iStrucGI, Int4 iFileID, CharPtr pcSeq, Int2 iType, Boolean bFirst)
{
	FILE *fp = NULL;
	Char path[PATH_MAX];
	Char pcDefline[PATH_MAX];
	CharPtr pcCompressedSeq = NULL;
	CharPtr end[] = {"seq","str","ss"};
	CharPtr pcend = NULL;

	if(iType == CLUSTSEQ) {
		pcend = end[0];
	} else if (iType == CLUSTSTR) {
		pcend = end[1];
	} else if (iType == CLUSTSS) {
		pcend = end[2];
	} else {
		ErrPostEx(SEV_ERROR,0,0,"Invalid Type");
		return;
	}

	sprintf(path,"database_%ld.%s",(long) iFileID, pcend);
	if(bFirst) fp = FileOpen(path,"w");
	else fp = FileOpen(path,"a+");
	if(fp == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Could not open file %s",path);
		return;
	}
	sprintf(pcDefline,"%ld_%ld.",(long) iSeqGI, (long) iStrucGI);
	CompressSequence(pcSeq,&pcCompressedSeq);
	DumpGORFASTAFile(pcDefline, pcCompressedSeq, fp);
	FileClose(fp);
}





static Boolean ResidueHasLoc(PMGD pmgdThis)
{
	PVNMA pvnmaThis;
	PMAD pmadThis;
	Boolean bFound;

	if (!pmgdThis) return FALSE;
	bFound = FALSE;
	pvnmaThis = pmgdThis->pvnmaAHead;
	while (pvnmaThis) {
		pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
		if (pmadThis->pvnalLocate) {
			bFound = TRUE;
		}

		pvnmaThis = pvnmaThis->next;
	}
	return bFound;
}

Char GetDomainLetterFromNumber(Uint1 num)
{
	Char alpha[] = "ABCDEFGHIJKLMNOPQRSTUVXYZ";
	if(num < 24) return alpha[num];
	else if (num < 48) return tolower(alpha[num]);
	else return '-';
}


Int2 CLUSTAL_GetSequenceAlignmentInfo(PDNMG pdnmgHead, ByteStorePtr bsp, FloatLo fDomIDThreshold, FloatLo fDomOccThreshold, 
			ValNodePtr *pvndomids, ValNodePtr *pvndomocc, ValNodePtr *pvndomlen, ValNodePtr *pvndomok, Int4 *piDomCount, Int4 *piDomOverThreshold)
{
	Boolean bIsOk = FALSE;
	Char ch;
	Int4 iDom = 0, iDomRes = 0, iDomainCt = 0;
	FloatLo flDomID = 0.0, flDomOcc = 0.0;
	PDNMG pdnmgThis = NULL, pdnmgHere = NULL;
	PMGD pmgdThis = NULL, pmgdHere = NULL;

	
	/* get the number of domains */
	pdnmgThis = pdnmgHead;
	while (pdnmgThis) {
		if((pmgdThis = (PMGD) pdnmgThis->data.ptrvalue) != NULL) {
			if ((Int4) pmgdThis->iDomain > iDomainCt) {
				iDomainCt = (Int4) pmgdThis->iDomain;
			}
		}
		pdnmgThis = pdnmgThis->next;
	}
	if(piDomCount) *piDomCount = iDomainCt;
	if(piDomOverThreshold) *piDomOverThreshold = 0;

	for (iDom = 0; iDom <= iDomainCt; iDom++) {
		flDomID = flDomOcc = 0.0;
		iDomRes = 0;
		pdnmgThis = pdnmgHead;
		while (pdnmgThis) {
			pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
			if ((Int4) pmgdThis->iDomain == iDom) {
				/* accumulate number of identities and domain occupancy per domain too */
				if (((char) pmgdThis->bReserved != '\0') && ((char) pmgdThis->bReserved !='X') ) {
					flDomOcc++;
				}
				if ((pmgdThis->pcIUPAC) && ( toupper((char) pmgdThis->bReserved) == (char) pmgdThis->pcIUPAC[0])) {
					flDomID++;
				}
				iDomRes++;
			} else if ((Int4) pmgdThis->iDomain == -(iDom)) {
				iDomRes++;
			}
			pdnmgThis = pdnmgThis->next;
		}
	
		if (iDomRes > 0) {
			flDomID  = flDomID / (FloatLo)iDomRes;
			flDomOcc = flDomOcc / (FloatLo)iDomRes;
		}
		ch = GetDomainLetterFromNumber((Uint1)iDom);
		if (( flDomOcc < fDomOccThreshold) || (flDomID < fDomIDThreshold)) {
			if(bsp && iDomRes > 0) BSprintf(bsp,"Domain [%c] masked out- [%6.2f]occupied [%6.2f]identical over [%5ld] res.\n",
				(char) ch, 
				(float) flDomOcc * 100,
				(float) flDomID * 100,
				(long)  iDomRes);
			pdnmgHere = pdnmgHead;
			while (pdnmgHere) {
				pmgdHere = (PMGD) pdnmgHere->data.ptrvalue;
				if ( ((Int4) pmgdHere->iDomain == iDom) && ((char) pmgdHere->bReserved != '\0') && 
					((char) pmgdHere->bReserved !='X') ) {
					pmgdHere->bReserved = 'X';
				}
				pdnmgHere = pdnmgHere->next;
			}
		} else {
			if(bsp && iDomRes > 0) BSprintf(bsp,"Domain [%c] is OK     - [%6.2f]occupied [%6.2f]identical over [%5ld] res.\n",
				(char) ch,  
				(float) flDomOcc * 100, 
				(float) flDomID * 100, 
				(long)  iDomRes);
			if(piDomOverThreshold) (*piDomOverThreshold)++;
			bIsOk = TRUE;
		}
		if(pvndomids) ValNodeAddFloat(&(*pvndomids),0,flDomID*100);
		if(pvndomocc) ValNodeAddFloat(&(*pvndomocc),0,flDomOcc*100);
		if(pvndomlen) ValNodeAddInt(&(*pvndomlen),0,iDomRes);
		if(pvndomok)  ValNodeAddInt(&(*pvndomok),0,(Int4) bIsOk);
		
		bIsOk = FALSE;
	}  
	return 0;
}
/* This function takes in a 2 node pcsan, containing 1- the aligned structure
sequence and structure and 2 - the aligned sequence to be mapped to the structure
 
The function traverses domains in target molecule 
  traverse graphs - count number of residues in domain "d"
  count number of occupied, non-'X' bReserved values in domain "d" 
  
  if bReserved/number of residues is under 75 % 
  or if under 25% identities within the domain
  X out the bReserved values in domain - this domain fails

can write out the fasta formatted sequence corresponding to the PKB sequence...
  traverse molecules - if bReserved, output its value, otherwise output X. 
    (extremophile total fractional sequence)
  traverse molecules - if bReserved, output the 1-letter code of the molecule, otherwise X (control)
    (original structure control fractional sequence)

  traverse molecule - if bReserved !=X
    traverse the atoms of the backbone - output PDB atom record lines;
    if bReserved = (lowercase) traverse all the atoms of the graph - output PDB atom record lines

pcSeq   = NVAFSKNFVSLTS-----KSSAAVGYEKLPISLVSS-----PSKKYLNEYLEITKEILNLANYNVH---------
pcStr   = NVVQRAEIRRMTVIEYDPKAKQADEYRALARKVVDNKLLVIPNPITM-----LLMEFGIMEVEDESIVGKTAEEV
pcMerge = nvAFSKNFVSLTSXXXXXkSSAaVGyEKlPISLvSSXXXXXpSKKYL-----ITKeILNLANYNVHXXXXXXXXX
pcCompr = nvAFSKNFVSLtSXXXXXkSSAaVGyEKlPISLvSSXXXXXpSKKYLITKeILNLANYNVHXXXXXXXXX

  NVkkkkkkkkkkk-----kkkkAkkYkkkkkkkVkk-----PkkkkkkkkEkkkkkkkkkk---------

  SCWRL data PDB output (k = backbone coordinates)
NVkkkkkkkkkkkkkkkAkkYkkkkkkkVkkPkkkkkkkkEkkkkkkkkkk

  SCWRL target sequence:  (remove X'es)
nvAFSKNFVSLTSKSSAaVGyEKLPISLvSSpSKKYLITKeILNLANYNVH


Allocate a sequence of same size length;
copy over sequence onto new one
  - SEQUENCE[i] (uppercase) i
  - if SEQUENCE[i] == STRUCTURE[i] (lowercase) i
  - if STRUCTURE '-' copy '-'
  - if SEQUENCE  '-' copy 'X'

make version of resulting sequence lacking '-'  (should match original structure in length)
poke values into bReserved along pdnmg list.

*/

Int4 MakeSequenceToDomainStructureAlignment
				(PCSAN pcsanAlign, Boolean bSecStruc,
				FloatLo fDomIDThreshold, FloatLo fDomOccThreshold,
				Int4 PNTR pidomcount, Int4 PNTR pidomthresh,
				ValNodePtr PNTR pvndomids, ValNodePtr PNTR pvndomocc, ValNodePtr PNTR pvndomlen, ValNodePtr PNTR pvndomok)
{
	PCSAN pcsanStructure = NULL, pcsanSequence = NULL;
	Char pcFileName[PATH_MAX];
	CharPtr pcStr = NULL,   pcSeq = NULL;
	CharPtr pcMerge = NULL, pcMtemp = NULL;
	CharPtr pcCompr = NULL, pcCtemp = NULL;
	CharPtr pcSecStruc = NULL;
	CharPtr pcLoop=NULL, pcLtemp = NULL;
	CharPtr pcDomain=NULL, pcDtemp = NULL;
	PDNMG pdnmgThis = NULL;
	PMMD pmmdThis = NULL;
	PMGD pmgdThis = NULL;
	Int4 i = 0, err;
	Int4 iStrLen = 0, iSeqLen = 0, iLnum = 0, iLtemp = 0;
	Boolean bUseRealSeq = FALSE;
	CharPtr pcControl = NULL, pcQuery = NULL;
	FILE *fp = NULL;
	ByteStorePtr bsp;
	Boolean bFirst = TRUE;

	
	if (pcsanAlign->pdnmsStructure) {
		pcsanStructure = pcsanAlign;
		pcsanSequence = pcsanAlign->next;
	} else {
		pcsanStructure = pcsanAlign->next;
		pcsanSequence = pcsanAlign;
	}

	pcStr = pcsanStructure->pcSeqAln;
	pcSeq = pcsanSequence->pcSeqAln;
	if((iStrLen = StringLen(pcStr)) != (iSeqLen = StringLen(pcSeq))) {
		ErrPostEx(SEV_ERROR,0,0,"Template Structure and Model Sequence are not the same length!");
		return 1;
	}

	if( ((pcMerge = (CharPtr) MemNew ((size_t) sizeof(Char) * (pcsanStructure->iLen + 1))) == NULL) ||
		((pcCompr = (CharPtr) MemNew ((size_t) sizeof(Char) * (pcsanStructure->iLen + 1))) == NULL) ||
		((pcDomain =(CharPtr) MemNew ((size_t) sizeof(Char) * (pcsanStructure->iLen + 1))) == NULL) ||
		((pcLoop  = (CharPtr) MemNew ((size_t) sizeof(Char) * (pcsanStructure->iLen + 1))) == NULL) ) {
			ErrPostEx(SEV_ERROR,0,0,"Unable to allocate memory");
			return 1;
	}

	bsp = BSNew(10);
	BSSeek(bsp, 0L, SEEK_SET);
	BSprintf(bsp,"SEQ GI: %-6ld\n",(long)pcsanSequence->iGi);
	BSprintf(bsp,"STR GI: %-6ld\tMMDBID %-ld\n",(long)pcsanStructure->iGi,(long)pcsanStructure->iMMDB);

	/* Mark an X where there is a gap or no alignment */
	pcMtemp = pcMerge;
	pcCtemp = pcCompr;
	pcLtemp = pcLoop;
	while (*pcSeq != '\0') {
		if (*pcSeq == '-') {                 /* sequence gap */
			*pcCtemp = *pcMtemp = 'X';       /* mask out residue */
			pcCtemp++;
			
			*pcLtemp = *pcSeq;
			iLtemp = 0;
			pcLtemp++;
		} else {
			if (*pcStr == '-') {             /* structure gap */
				*pcMtemp = '-';	             /* keep gap */
				if(!bFirst) {*pcLtemp = *pcSeq;pcLtemp++;iLtemp++;iLnum++;}
			
			} else {
				bFirst = FALSE;
				*pcLtemp = '-';
				iLtemp = 0;
				pcLtemp++;
				/* aligned region */
				if (*pcSeq == *pcStr) {      /* identity */
					*pcMtemp = (Char) tolower(*pcSeq); 
				} else {
					*pcMtemp = *pcSeq;           /* non-identity */
				}
				*pcCtemp = *pcMtemp;             /* compress to structure */
				pcCtemp++;
				
			}
		}
		pcSeq++;
		pcStr++;
		pcMtemp++;		
	}
	*pcMtemp = '\0';
	*pcCtemp = '\0';
	pcLtemp[-iLtemp] = '\0';

	/* Write the compressed alignment residue into the reserved block of the structure */
	if((pmmdThis = MMDBGetMoleculeByGI(pcsanStructure->pdnmsStructure, pcsanStructure->iGi)) == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Unable to get mmdb molecule by GI %ld", (long) pcsanStructure->iGi);
		return 1;
	}
	pdnmgThis = pmmdThis->pdnmgHead;

	/* Invert the domain number if no coordinates found */
	pcCtemp = pcCompr;
	pcDtemp = pcDomain;
	while ((pdnmgThis) && (*pcCtemp != '\0')) {
		if((pmgdThis = (PMGD) pdnmgThis->data.ptrvalue) != NULL) {
			*pcDtemp = (Char) GetDomainLetterFromNumber((Uint1)pmgdThis->iDomain);
			if (!ResidueHasLoc(pmgdThis)) {
				pmgdThis->bReserved = 'X';
				pmgdThis->iDomain = (Int2) -pmgdThis->iDomain; 
				i++;
			} else {
				pmgdThis->bReserved = (Byte) *pcCtemp;
			}
			pcDtemp++;
			pcCtemp++;
		}
		pdnmgThis = pdnmgThis->next;
	}
	
	BSprintf(bsp,"Initial Seq = [%s]\n",pcsanSequence->pcSeqAln);
	BSprintf(bsp,"Initial Str = [%s]\n",pcsanStructure->pcSeqAln);
	BSprintf(bsp,"Merged  Seq = [%s]\n",pcMerge);
	BSprintf(bsp,"Comprs  Seq = [%s]\n",pcCompr);
	BSprintf(bsp,"Domain      = [%s]\n",pcDomain);
	if((iLnum - iLtemp) > 0) {
		BSprintf(bsp,"Loops       = [%s]\n",pcLoop);
	}
	pdnmgThis = pmmdThis->pdnmgHead;

	/* evaluate each domain according to threshold criteria */
	/* remember that domains can span multiple sequence segments */
	if((err = CLUSTAL_GetSequenceAlignmentInfo(pmmdThis->pdnmgHead,bsp,fDomIDThreshold,fDomOccThreshold,
		&(*pvndomids),&(*pvndomocc),&(*pvndomlen),&(*pvndomok),&(*pidomcount),&(*pidomthresh))) > 0) {
		ErrPostEx(SEV_ERROR,0,0,"Unable to get sequence alignment info");
		return err;
	}
	BSprintf(bsp,"Missing Coordinates: %ld\n", (long) i);

	if((*pidomthresh > 0) && (bsp != NULL)) {
		sprintf(pcFileName,"%s%s",S2S_HEAD,/*(long)itaxid,*/S2S_TAIL);
		if((fp = FileOpen(pcFileName,"a+")) == NULL) return 1;
		BSprintf(bsp,";\n");
		PrintBS(bsp,fp);
		FileClose(fp);
		
		/* Now Get the aligned sequence */
		bUseRealSeq = FALSE;
		GetSeqInfoFromStruc(pmmdThis, FALSE, bUseRealSeq, &pcQuery);
		ValNodeAddStr(&pcsanSequence->pvnSeqAlnPart, 0,  pcQuery);
		
		bUseRealSeq = TRUE;
		GetSeqInfoFromStruc(pmmdThis, FALSE, bUseRealSeq, &pcControl);
		ValNodeAddStr(&pcsanStructure->pvnSeqAlnPart, 0,  pcControl);

		iStrLen = StringLen(pcControl);

		/* and secondary structure, if requested */
		if(bSecStruc) {
			GetSeqInfoFromStruc(pmmdThis, bSecStruc, bUseRealSeq, &pcSecStruc);
			if(pcsanStructure->pExtra == NULL) {
				pcsanStructure->pExtra = pcSecStruc;
			} else {
				ErrPostEx(SEV_ERROR,0,0,"pscanStructure does not have a free pExtra field for SecStruc");
				return 1;
			}
		}
	}
	BSFree(bsp);


	pcMerge = MemFree(pcMerge);
	pcCompr = MemFree(pcCompr);
	pcLoop = MemFree(pcLoop);
	pcDomain = MemFree(pcDomain);
	return 0;
}



/* BLOSUM MATRIX parsing */
Int2 LoadBLOSUM(CharPtr fnam)
{
	FILE *f;
	Char buf[PATH_MAX];
	Boolean first=TRUE;
	Char ch[2];
	CharPtr pbuf;
	Int2 *pihere=NULL;
	int num;
	
	StringCpy(blocol,"");
	if ((f=FileOpen(fnam,"r"))==NULL) {
		ErrPostEx(SEV_ERROR,1,2,"Unable to open file %s",fnam);
		return 1;
	}
	while (FileGets(buf,PATH_MAX,f)!=NULL) {
		if (buf[0]=='#')
			continue;
		if (first) {
			first=FALSE;
			pbuf=buf;
			while (sscanf(pbuf,"%s",ch)==1) {
				StringCat(blocol,ch);
				pbuf=StringStr(buf,ch)+1;
			}
			Blosum=(Int2 *)MemNew(sizeof(Int2)*StringLen(blocol)*StringLen(blocol));
			pihere=Blosum;
		}
		else {
			pbuf=buf+1;
			while (sscanf(pbuf,"%d",&num)==1) {
				*pihere=(Int2)num;
				pihere++;
				while (pbuf[0]==' ')
					pbuf++;
				while (pbuf[0]=='-' || (pbuf[0]>='0' && pbuf[0]<='9'))
					pbuf++;
			}
		}
	}
	FileClose(f);
	return 0;
}

Int2 GetBLOSUM(Char x, Char y)
{
	Char xx[2],yy[2];
	Int2 idx1,idx2;

	if (Blosum==NULL)
		return 0;
	xx[0]=x;
	xx[1]='\0';
	yy[0]=y;
	yy[1]='\0';
	idx1=(Int2)StringCSpn(blocol,xx);
	idx2=(Int2)StringCSpn(blocol,yy);
	if (idx1==(Int2)StringLen(blocol))
		return 0;
	if (idx2==(Int2)StringLen(blocol))
		return 0;
	return Blosum[idx1*StringLen(blocol)+idx2];
}


void FreeBLOSUM(void)
{
	if (Blosum!=NULL)
		Blosum=MemFree(Blosum);
}




/*
$Log: clustlib.c,v $
Revision 1.29  2004/11/30 16:23:08  hfeldman
Changed structure fetching to shound function

Revision 1.28  2004/06/14 14:49:07  hfeldman
Added missing fileclose

Revision 1.27  2003/12/01 23:02:47  mjdumont
Added extra logging

Revision 1.26  2003/12/01 22:44:51  mjdumont
Fixed broken call to SHoundGet3D.

Revision 1.25  2003/10/27 15:31:26  michel
Removed compiler warnings, updated projects for msvc dependency to NCBI libraries

Revision 1.24  2002/10/29 14:28:54  michel
removed another memory leak

Revision 1.23  2002/10/29 14:08:23  michel
moved BSFree out of conditional statement

Revision 1.22  2002/10/26 16:20:22  michel
minor fix

Revision 1.21  2002/10/21 17:48:38  michel
Changes to alignment logging - added loop & domain info

Revision 1.20  2002/06/26 21:51:41  michel
Amended compression step and alignment info

Revision 1.19  2002/06/25 23:03:31  michel
Added mkmdl and usemdl to windows gcpp workspace
Added function to clustal to apply new %id and %occ thresholds and
	extract alignment info (%id,%occ,domlen,etc) from structure models
Fixed gcpp/usemdl counting bugs and taxonomy calls

Revision 1.18  2002/06/14 22:47:56  michel
Addition of parameter

Revision 1.17  2002/06/10 17:23:41  rong
Fixed memory allocation bug

Revision 1.16  2002/03/30 22:02:16  michel
modified GetSeqInfoFromStruc for other types of residue info

Revision 1.15  2001/11/27 23:28:34  michel
removed compiler warning

Revision 1.14  2001/11/23 23:40:07  michel
added gilgor database formatting, modification to sequence to structure alignment function

Revision 1.13  2001/11/09 21:50:42  michel
Added domain ok list return on sequence to structure alignment

Revision 1.12  2001/09/10 16:49:38  michel
added SLRI_MMDBBiostrucGet to call MMDBGetBiostrucGet (for clustal library)

Revision 1.11  2001/09/09 05:38:13  michel
Added more returns for Sequencetostructure function
library can be linked by seqhound with linker flag to suppress multiply defined symbols

Revision 1.10  2001/09/04 21:25:44  michel
Moved mmdblocl.c out of seqhound local library so that it can be directly compiled with the cbmmdb parser
This allows clustlib to use MMDBBiostrucGet with SeqHound (which calls SHoundGet3D)
or without SeqHound (in which case mmdblocl.c must be compiled with application)
SeqHound's MMDBBiostrucGet in mmdblocl.c has been renamed to MMDBBiostrucGetHere

Revision 1.9  2001/08/29 17:45:08  michel
minor change

Revision 1.8  2001/08/28 21:48:57  michel
Removed printfs for ErrPostEx (set your message/log level appropriately)
Added paths for library if compiled with seqhound

Revision 1.7  2001/08/28 16:59:20  hogue
Fixed memory bug and changed function parameters

Revision 1.6  2001/08/27 17:35:19  michel
minor changes

Revision 1.5  2001/08/20 20:48:55  michel
memory Bug fixes

Revision 1.4  2001/08/13 17:39:30  michel
minor fix to keep crosscomp quiet

Revision 1.3  2001/06/27 01:37:01  michel
Added old slapamol functions to library with modifications and compile time
flags for seqhound or entrez access to sequences/structures

Revision 1.2  2001/06/26 20:40:09  feldman
Added makefile for UNIX,
and removed some functions that should not be there

Revision 1.1  2001/06/26 20:18:30  feldman
CLUSTAL/SLRI CSAN library


*/
