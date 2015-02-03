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


#include <vismaketrj.h>
/* seqhound dependency removed */

VoidPtr MakeTrj(VoidPtr ptr) {
	pMakeTrjParamBlock paramblock;
	PMSD pmsd; 
	PMMD pmmd=NULL;
	PMGD pmgd;
	PDNMG pdnmg;
	PDNMM pdnmm;
	Char sequence[MAXSIZE];
	Char encseq[MAXSIZE];
	Char xxsequence[MAXSIZE];
	Char dummysequence[MAXSIZE];
	Char tmpfnamout[PATH_MAX];
	Char filteredfnamout[PATH_MAX];
 	Char valfileext[PATH_MAX];
	CharPtr pfnamout;
	Char ftmp[PATH_MAX];
	Int4 numres;
	Int2 i,model,idx;
	FILE *fp;
	FILE *fpAtmInfo;
	Char buf[MAXCOL];
	Int2 numAA=0;
	Char chain[2];
	DValNodePtr pdnIncList = NULL;
	DValNodePtr pdnZhangAtmList = NULL;
	DValNodePtr pdnListPmmds;
	Int4Ptr ptnlFunc = NULL;
	BiostrucPtr		pbs;
	AsnIoBSPtr		aibp=NULL;
	ByteStorePtr		bsid=NULL,bsdescr=NULL;
	ValNodePtr		vnpDescr,vnpId;
	ValNodePtr		vnpBioseq=NULL;
	PMLD pmldHere;
	Int2 iDotLen  = 0;
#ifdef OS_MSWIN
	Char tmpfnam[PATH_MAX];
#endif

	paramblock=(pMakeTrjParamBlock)ptr;
	if (paramblock==NULL) return NULL;
	/*  RandomSeed(GetSecs());*/
	if (paramblock->trjfnamout==NULL) {
		if (paramblock->valfnamin==NULL) {
			StringCpy(tmpfnamout,DEFAULT_TRJ_NAME);
		} else {
			StringCpy(tmpfnamout,paramblock->valfnamin);
			/* remove paths from output filename when copying input filename */
			pfnamout=tmpfnamout+StringLen(tmpfnamout)-1;
			while (pfnamout>tmpfnamout) {
				if (pfnamout[0]=='/' || pfnamout[0]=='\\') {
					pfnamout++;
					break;
				}
				pfnamout--;
			}
			if (pfnamout>tmpfnamout) {
				idx=0;
				do {
					tmpfnamout[idx]=pfnamout[idx];
					idx++;
				} while (pfnamout[idx-1]);
			}					
		}
	} else {
		StringCpy(tmpfnamout, paramblock->trjfnamout);
	}

	StringCpy(filteredfnamout,tmpfnamout);
	if (paramblock->valfnamin!=NULL) {
		StringCpy(valfileext,paramblock->valfnamin);
  
		/* add extension only if one not already there! we may pass in a .prt or .cn3 file */
		iDotLen = 0;
            	iDotLen = (Int2) StringLen(valfileext);
           	if (iDotLen > 5) {
		   if (valfileext[iDotLen - 4] != '.') {
		        StringCat(valfileext,MMDB_EXT);	         
		   }
                }
	}


        paramblock->valfnamin = valfileext;  /*VALGRIND overwrite error StringCat(paramblock->valfnamin,MMDB_EXT);*/
	


	chain[0]=0;
	chain[1]=0;
	if (paramblock->pcChain) {
		chain[0]= paramblock->pcChain[0];
	}
	model=paramblock->modelnum;
	if (paramblock->constrfnamin!=NULL) {
		StringCpy(CONSTRAINT_FILE,paramblock->constrfnamin);
	} else {
		StringCpy(CONSTRAINT_FILE,"");
	}
	if (paramblock->TrajMethod == 2) { /* 2 is from sequence */
		TRAJTYPE=paramblock->tgtype;
	} else {
		TRAJTYPE=TRAJ_NA;  /*  from structure */ 
	}

	/* all trajectory graphs created by default have arbitrary units */
	if (paramblock->units>0) {
		TGUNITS=paramblock->units;
	} else {
		TGUNITS=UNITS_ARBITRARY;
	}

	/* ensure unique name is chosen since we are adding extensions
	to the base name */
	StringCpy(tmpdbasename,DFPTmpNam(TRUE,TRUE));  
	if (traj_quiet == VERBOSITY_VERBOSE)
	printf("Reading parameter file..\n");
	/* read .inittrajrc */
	LoadParams();
	/* section for InitTraj */
	if (paramblock->TrajMethod == 2) {
			if(paramblock->seqfnamin == NULL && paramblock->alignfile == NULL) {
				ErrPostEx(SEV_ERROR,3,2,"No sequence filename passed in input!");
				ProgramProgressMax=-999;
				return NULL;
			}
			if(paramblock->seqfnamin && paramblock->alignfile == NULL) {
				/* Make sure the file is set */
				if((fp=FileOpen(paramblock->seqfnamin,"r"))==NULL) {
					ErrPostEx(SEV_ERROR,1,1,"Error - unable to find input file %s",paramblock->seqfnamin);
					ProgramProgressMax=-999;
					return NULL;
				}
				else {
					if (traj_quiet == VERBOSITY_VERBOSE) printf("Reading FASTA sequence file %s\n",paramblock->seqfnamin);
					StringCpy(sequence,"");
					numAA=0;
					while (FileGets(buf,MAXCOL,fp)!=NULL) {
						if (buf[0]!='>') {
							StringCat(sequence,buf);
						} else {
							while (buf[StringLen(buf)-1]!='\n') {
								if (FileGets(buf,MAXCOL,fp)==NULL) {
									ErrPostEx(SEV_ERROR,3,2,"No sequence found in FASTA sequence file %s!",paramblock->seqfnamin);
									ProgramProgressMax=-999;
									return NULL;
								}
							}
						}
					}
					FileClose(fp);
	   				numAA=ConvertExtSequence(sequence,dummysequence,EXTAA_X);
					if (numAA==0) {
						ErrPostEx(SEV_ERROR,3,2,"Sequence FASTA file parsing reports 0 length of amino acids in sequence file %s!",paramblock->seqfnamin);
						ProgramProgressMax=-999;
						return NULL;
					}
					if (numAA<3) {
						ErrPostEx(SEV_ERROR,3,2,"Sequence in FASTA must be at least 3 residues long! %s",paramblock->seqfnamin);
						ProgramProgressMax=-999;
						return NULL;
					}
				}
			}
			/* otherwise will get sequence from alignment file */

			if (paramblock->constrfnamin!=NULL) {
				if ((LoadDistConstraints(paramblock->constrfnamin,sequence)) != ERR_SUCCESS) {
					ErrPostEx(SEV_ERROR,3,1,"Unable to load constraint file %s", paramblock->constrfnamin);
					ProgramProgressMax=-999;
					return NULL;
				}
			}
		
				if (traj_quiet==VERBOSITY_SILENT) {
					ProgramProgress=1;
					if (numAA)
						ProgramProgressMax=numAA;
					else {
						/* should never get here.. */
						ErrPostEx(SEV_ERROR,3,1,"Unable to determine protein length, aborting.");
						ProgramProgressMax=-999;
						return NULL;
					}
				}
/* CWVH Shortened call, passing parameter block, Model is not used */				
			MakeTG(sequence,paramblock,tmpfnamout);
			if (traj_quiet == VERBOSITY_VERBOSE)
				printf("Done.\n");
			if (traj_quiet == VERBOSITY_QUIET)
				printf("\n");
	} else {
		/* below is only for Val2Traj/ Unfoldtraj */
		model=1;
		pmsd=LoadABiostrucEx(paramblock->valfnamin,FALSE,ALLMDL,&model,&vnpBioseq);
		if (pmsd == NULL){
			ErrPostEx(SEV_ERROR,5,1,"Unable to determine number of models");
			ProgramProgressMax=-999;
			return NULL;
		}
		model=GetFirstFullModel(pmsd);
			if (!model){
			ErrPostEx(SEV_ERROR,5,2,"No all-atom models found");
			ProgramProgressMax=-999;
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		pmldHere=GetModelN(pmsd,paramblock->modelnum);
		if (pmldHere==NULL) {
			ErrPostEx(SEV_ERROR,5,3,"Model %d not found in structure!",paramblock->modelnum);
			ProgramProgressMax=-999;
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		if (pmldHere->iType!=Model_type_ncbi_all_atom && pmldHere->iType!=Model_type_pdb_model) {
			ErrPostEx(SEV_ERROR,6,1,"Model %ld not all atom, overriding and using model %ld",(long)paramblock->modelnum,(long)model);
		} else {
			model=paramblock->modelnum;
		}
		/* write id and descr to temporary bytestores */
		pbs = pmsd -> pbsBS;
		/* arbitrary sizes */
		bsid=BSNew(sizeof(Char)*1);
		bsdescr=BSNew(sizeof(Char)*1);
		BSSeek(bsid,0L,SEEK_SET);
		BSSeek(bsdescr,0L,SEEK_SET);
		vnpId=pbs->id;
		aibp=AsnIoBSOpen("wb",bsid);
		while (vnpId!=NULL) {
			BiostrucIdAsnWrite(vnpId,aibp->aip,NULL);
			vnpId=vnpId->next;
		}
		aibp=AsnIoBSClose(aibp);
		vnpDescr=pbs->descr;
		aibp=AsnIoBSOpen("wb",bsdescr);
		while (vnpDescr!=NULL) {
			BiostrucDescrAsnWrite(vnpDescr,aibp->aip,NULL);
			vnpDescr=vnpDescr->next;
		}
		aibp=AsnIoBSClose(aibp);
		pdnmm = (PDNMM) (pmsd -> pdnmmHead);
     
		/* Selects chain to construct trajgraph.  Added case insensitivity */
		if (StringStr(VL_ALPHABET_LOWER, chain) || paramblock->pcChain == NULL){
			if (paramblock->pcChain) {
				while (pdnmm) {
					pmmd = (PMMD) (pdnmm -> data.ptrvalue);
					if (!StringCmp(chain, pmmd -> pcMolName)) {
						break;
					} else {
						pdnmm = pdnmm -> next;
					}
				}
			} else {
				pmmd = (PMMD) (pdnmm -> data.ptrvalue);
			}
			if (pdnmm == NULL) {
				chain[0] = toupper(chain[0]);	
				pdnmm = (PDNMM) (pmsd -> pdnmmHead);
			}
		}
		if (StringStr(VL_ALPHABET_UPPER, chain)) {
			if (paramblock->pcChain) {
				while (pdnmm) {
					pmmd = (PMMD) (pdnmm -> data.ptrvalue);
					if (!StringCmp(chain, pmmd -> pcMolName)) {
						break;
					} else {
						pdnmm = pdnmm -> next;
					}
				}
			} else  {
				pmmd = (PMMD) (pdnmm -> data.ptrvalue);
			}
		}	
		if (pdnmm == NULL) {
			ErrPostEx(SEV_ERROR,5,1,"No  %s  chain exists for this protein",chain);
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		} 
		pdnmg = (pmmd -> pdnmgHead);
		numres = (pmmd -> iResCount);
		pmgd = (PMGD) (pdnmg -> data.ptrvalue);
		if ((pmmd -> bWhat & AM_PROT) ==0){
			ErrPostEx(SEV_ERROR,6,1,"Not a protein, unable to continue");
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		if (paramblock->startres>numres){
			ErrPostEx(SEV_ERROR,6,1,"Start residue invalid");
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		if (paramblock->endres==0)
			paramblock->endres = (Int2) numres;
		if (paramblock->endres > numres){
			ErrPostEx(SEV_ERROR,7,1,"End residue invalid");
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		if (paramblock->startres>paramblock->endres){
			ErrPostEx(SEV_ERROR,8,1,"Start residue invalid");
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		if (paramblock->endres-paramblock->startres<2){
			ErrPostEx(SEV_ERROR,8,1,"Must have at least 3 resdidues");
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		i=0;
		/* fills in the sequence */
		while (pdnmg){
			if ((pdnmg->choice) >= paramblock->startres && (pdnmg ->choice) <= paramblock->endres) {
				pmgd = (PMGD) (pdnmg -> data.ptrvalue);
				sequence[i]= pmgd -> pcIUPAC[0];
				i++;
			}
			pdnmg = pdnmg -> next;
		}
		sequence[i]='\0';
		if (traj_quiet==VERBOSITY_QUIET || traj_quiet==VERBOSITY_VERBOSE) {
			printf("Progress:      %4d/%-4d ",paramblock->startres,paramblock->endres-paramblock->startres+1);
			fflush(stdout); 
		}  
		sprintf(ftmp,"%s%s",CFG_local_datafilepath,ZHANG_POTENTIAL);
		if ((fp = FileOpen (ftmp, "r")) == NULL) {
			ErrPostEx (SEV_ERROR, 1, 6, "File not found: %s",ftmp);
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		sprintf(ftmp,"%s%s",CFG_local_datafilepath,ZHANG_ATOMS);
		if ((fpAtmInfo = FileOpen (ftmp, "r")) == NULL) {
			ErrPostEx (SEV_ERROR, 1, 7, "File not found: %s",ftmp);
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}

		if (LoadZhangPotential (&ptnlFunc, fp)!=ERR_SUCCESS) {
			ErrPostEx (SEV_ERROR, 1, 8, "Error loading Zhang potential");
			ProgramProgressMax=-999;
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAModelstruc(pmsd->pdnmsLink);
			return NULL;
		}
		LoadZhangAtmList (&pdnZhangAtmList, fpAtmInfo);
		FileClose(fp);
		FileClose(fpAtmInfo);
		pdnListPmmds = CreateListPmmds (NULL, pmsd);
		/* DValNodePtr pdnmm is of chain selected */
		ComputeZhangPotential (&pdnIncList, (DValNodePtr) (pdnmm), model,
		/*bInclusiveWindow*/ FALSE, /*iWindowSize*/ paramblock->zhangwindsize,
		ptnlFunc, pdnZhangAtmList, FALSE /*detailed list?*/);
		if (BuildTrjGraph(pmmd,model,paramblock->startres,paramblock->endres,paramblock->peakheight,paramblock->sigma_x,paramblock->sigma_y,tmpdbasename,paramblock->noise,paramblock->temperature,paramblock->timestep,pdnIncList,paramblock->savechi,paramblock->ssmask) != ERR_SUCCESS){
			if (ProgramProgress>=0) {
				ErrPostEx(SEV_ERROR,9,1,"Unable to build trajectory distribution");
			}
			BSFree(bsid);
			BSFree(bsdescr);
			FreeAdjList (&pdnIncList);
			FreeZhangAtmList (&pdnZhangAtmList);
			FreeZhangPotential (&ptnlFunc);  
			FreeListPmmds (&pdnListPmmds);
			FreeAModelstruc(pmsd->pdnmsLink);
			ProgramProgressMax=-999;
			return NULL;
		}
		if (paramblock->fragfile) {
			ConvertExtSequence(sequence,xxsequence,EXTAA_X);
			if (traj_quiet==VERBOSITY_VERBOSE)
				printf("Adding fragments from file %s..\n",paramblock->fragfile);
			if (ImportFrags(paramblock->fragfile,xxsequence)!=ERR_SUCCESS)
				ErrPostEx(SEV_ERROR,1,1,"Unable to find fragment file %s",paramblock->fragfile);
		}
		ConvertExtSequence(sequence,encseq,EXTAA_ENCODE);
		if (paramblock->comprtype==USE_NONE) {
			PackAsnTrajGraph(tmpfnamout,encseq,bsid,bsdescr,vnpBioseq,TRUE); 
		} else {
			PackAsnTrajGraph(tmpfnamout,encseq,bsid,bsdescr,vnpBioseq,FALSE); 
		}
		BSFree(bsid);
		BSFree(bsdescr);
		CleanUpDB(tmpdbasename);
		FreeAModelstruc(pmsd->pdnmsLink);
		FreeAdjList (&pdnIncList);
		FreeZhangAtmList (&pdnZhangAtmList);
		FreeZhangPotential (&ptnlFunc);  
		FreeListPmmds (&pdnListPmmds);
                SeqEntryFree(vnpBioseq); /* now freed outside of PackAsnTrajGraph by caller of creator LoadABiostrucEx */
	}
	fflush(stdout);
	ProgramProgress=0;
	ProgramProgressMax=0;
	return NULL;
}

/*  
$Log: vismaketrj.c,v $
Revision 1.49  2004/09/24 20:31:37  hfeldman
Print error message if fragment file not found

Revision 1.48  2004/09/24 19:09:36  hfeldman
Added import fragment option to maketrj

Revision 1.47  2004/07/27 22:46:25  ksnyder
Pass value of gap shifting command line argument to HomTraj function

Revision 1.46  2004/07/15 22:38:16  ksnyder
Write errors, causing 'goto hombail', to Maketraj error log

Revision 1.45  2004/07/13 15:45:32  hfeldman
Changed behaviour to ignore seq file if aln file given

Revision 1.44  2004/07/12 16:57:59  mjdumont
HomTraj doesn't require query FASTA - obtained from alignments.\nLow complexity filtering in

Revision 1.43  2004/07/09 22:10:14  hfeldman
Change X to A in maketrj sequences whenever -A param is given (homology modelling)

Revision 1.42  2004/07/06 20:25:24  hfeldman
Fixed missing parameter

Revision 1.41  2004/07/06 20:19:32  hfeldman
Changed arguments to HomTraj call

Revision 1.40  2004/06/29 15:59:26  egarderm
Changed all FreeDNMS call to FreeAModelstruc calls - otherwise memory was being freed incorrectly, and was being corrupted.

Revision 1.39  2004/06/17 22:03:21  hfeldman
Added argument to homtraj

Revision 1.38  2004/06/09 17:23:55  mjdumont
Broke out RPSBlastTraj mega-function to HomTraj and Hom_* functions

Revision 1.37  2003/12/19 21:38:38  hfeldman
convert Xs to As for homtraj sequence

Revision 1.36  2003/11/07 18:08:53  ksnyder
added path to blosum matrix

Revision 1.35  2003/09/23 17:54:48  feldman
Fixed up error messages for homsrv

Revision 1.34  2003/09/23 17:09:30  feldman
Added error message output for homtraj

Revision 1.33  2003/09/14 02:38:23  feldman
Fixed unused variables and other minor compiler warnings

Revision 1.32  2003/07/17 13:31:50  feldman
Made universal DFPTmpNam function

Revision 1.31  2003/07/15 19:35:51  feldman
Made temp file naming more robust, especially if using DFPTEMP

Revision 1.30  2003/07/15 19:06:12  feldman
Added some extra freeing statements on certain errors

Revision 1.29  2003/07/15 13:56:46  feldman
Changed field name of structure

Revision 1.28  2003/07/14 20:01:01  egarderm
Changed TASK_CFGFILE to .\foldtraj on Windows to allow local INI file

Revision 1.26  2003/04/04 21:54:04  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.25  2003/03/25 17:46:39  feldman
added zhang window param to maketrj

Revision 1.24  2003/02/27 21:20:09  feldman
Fixed tmp name bug on windows

Revision 1.23  2003/01/24 16:50:54  feldman
Dont filter filenames except as CGI and made minor change to error msg

Revision 1.22  2003/01/08 16:52:33  feldman
maketrj now reads both ascii and binary val input files

Revision 1.21  2002/07/29 19:42:21  phan
Template option added
E-value back to 10.0

Revision 1.20  2002/07/08 20:32:20  phan
Added creation of fragments. Fixed bugs.

Revision 1.19  2002/07/07 21:12:02  feldman
Finished TG building section

Revision 1.18  2002/02/07 21:21:07  feldman
Fixed Darwin fix

Revision 1.17  2001/12/06 23:31:07  feldman
Fix in tmpnaming for Darwin

Revision 1.16  2001/12/03 15:54:04  feldman
Can pass units to function now

Revision 1.15  2001/09/10 20:29:21  feldman
Added Uturn option to maketrj

Revision 1.14  2001/06/26 20:19:10  feldman
Moved stuff to clustal library

Revision 1.13  2001/06/21 16:19:38  feldman
Fixed small bug and Windows temp file issue

Revision 1.12  2001/05/25 21:48:04  feldman
Added functionality to allow most Foldtraj data files to be in a directory
separate from the executable (given in the config file for example)
Also added screensaver stuff to randwalk

Revision 1.11  2001/04/26 19:49:09  feldman
Fixed some minor Bioseq bugs and potential bugs

Revision 1.10  2001/03/30 22:22:27  feldman
Minor spelling and grammar fixes, bug fixes, and added
foldtraj ASCII screensaver ability

Revision 1.9  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.8  2001/03/27 20:24:22  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.7  2001/03/13 15:08:17  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.6  2001/03/07 21:49:47  feldman
Made many important fixes to VisTraj such as bounds checking of
user-entered information, and tidied up interface further.
Fixed a few minor bugs as well, and tested RPSTraj feature of
VisTraj.

Revision 1.5  2001/02/27 20:52:32  feldman
minor bugfix

Revision 1.4  2001/02/23 18:03:54  feldman
Changed MakeTrj to void*, void*

Revision 1.3  2001/02/15 20:27:57  feldman
Reworked maektrj to take a parameter block instead of a bunch
of individual parameters.

Revision 1.2  2001/02/09 20:17:40  feldman
Fixed closemmdbapi bug and default filenames

Revision 1.1  2001/02/06 20:41:13  feldman
Moved to library directory, file modified to fit in foldtrajlib library

Revision 1.2  2001/01/25 22:25:26  feldman
Integrated Maketrj, changed lots of text messages, fixed
many minor bugs, improved overall interface

Revision 1.1  2001/01/25 14:36:08  feldman
First attempt to merge maketrj into vistraj


*/

