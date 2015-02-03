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

#include <hfprogs.h>
#include <slri_misc.h>

#define NUMARGS 18


/* Global Variables */
Args Rargs[NUMARGS] = {
	    {"Input MMDB filename 1",NULL,NULL,NULL,TRUE,'f',ARG_FILE_IN,0.0,0,NULL},
		{"Input MMDB filename 2",NULL,NULL,NULL,TRUE,'g',ARG_FILE_IN,0.0,0,NULL},
		{"Chain 1 (default = first in file)","-",NULL,NULL,TRUE,'h',ARG_STRING,0.0,0,NULL},
		{"Chain 2 (default = first in file)","-",NULL,NULL,TRUE,'i',ARG_STRING,0.0,0,NULL},
		{"Model Number 1","1","1","9999",TRUE,'m',ARG_INT,0.0,0,NULL},
		{"Model Number 2","1","1","9999",TRUE,'n',ARG_INT,0.0,0,NULL},
		{"Window size (0=whole protein)","0","0","32767",TRUE,'w',ARG_INT,0.0,0,NULL},
		{"Align: 0 = All Atoms; 1=BB; 2=CA","0","0","2",TRUE,'a',ARG_INT,0.0,0,NULL},
		{"res 1 mol 1 (0=use -w option instead)","0","0","32767",TRUE,'b',ARG_INT,0.0,0,NULL},
		{"res 2 mol 1 (0=use -w option instead)","0","0","32767",TRUE,'c',ARG_INT,0.0,0,NULL},
		{"res 1 mol 2 (0=use -w option instead)","0","0","32767",TRUE,'d',ARG_INT,0.0,0,NULL},
		{"res 2 mol 2 (0=use -w option instead)","0","0","32767",TRUE,'e',ARG_INT,0.0,0,NULL},
		{"Output supressed for R stdin (t/F)","0","0","0",TRUE,'o',ARG_BOOLEAN,0.0,0,NULL},
        	{"Input File - List of Filenames to Align",NULL,NULL,NULL,TRUE,'l',ARG_FILE_IN,0.0,0,NULL},
        	{"File List - Start at (use with -l)","0","0","10000",TRUE,'j',ARG_INT,0.0,0,NULL},
		{"Create File Queue from directory","0","0","0",TRUE,'q',ARG_BOOLEAN,0.0,0,NULL},
		{"Write PDB output file (t/F)","0","0","0",TRUE,'p',ARG_BOOLEAN,0.0,0,NULL},
	        {"Write .prt ASCII text Asn.1 output file (t/F)","0","0","0",TRUE,'t',ARG_BOOLEAN,0.0,0,NULL}
        };

/*ed
#Row 1 salign parameter is -j for starting position
# -l for list of files
# and parameter -o T supresses writing out rotated structure

*/



#define TEMP_FILE_NAME "val_file_list.txt"
void GenerateFileNamesFromDirectory(CharPtr pcWildCard)
{
	Char pcCommand [PATH_MAX];

#ifdef WIN32
	sprintf(pcCommand,"dir /b %s > %s",pcWildCard, TEMP_FILE_NAME);
#else 
	sprintf(pcCommand,"ls -1 %s > %s",pcWildCard, TEMP_FILE_NAME);
#endif
	system(pcCommand);
}



ValNodePtr GetFileNamesFromFile(CharPtr pcFileName)
{
	ValNodePtr vnp = NULL, vnpHead = NULL;
	Char buf[PATH_MAX];
	CharPtr pcTemp = NULL;
	FILE *fp = NULL;

	if((fp = FileOpen(pcFileName,"r")) == NULL) {
		ErrPostEx(SEV_ERROR,0,0,"Invalid File Name");
		return NULL;
	}

	while(FileGets(buf,PATH_MAX,fp)) {
		pcTemp = StringChr(buf,'\r');
		if(pcTemp != NULL) *pcTemp = '\0';
		pcTemp = StringChr(buf,'\n');
		if(pcTemp != NULL) *pcTemp = '\0';

		vnp = ValNodeCopyStr(&vnp, 0, buf);
		if(vnpHead == NULL) vnpHead = vnp;
	}
	FileClose(fp);

	return vnpHead;
}




Int2 Main()
{
	PMSD pmsd1,pmsd2, pmsd3;
	PMMD pmmd1,pmmd2;
	PDNMM pdnmm;
	Int2 ModelNum1,ModelNum2,atomtypes,pstart,wsize,a1,a2,b1,b2;
	Char chain1,chain2;
	Int4 start_i = 0;
    CharPtr pcStartstruc1 = NULL; /* first structure column on the diagonal */
    CharPtr pcStartstruc2 = NULL; /* second structure column on the diagonal */
	FILE *fp = NULL, *pFile = NULL;
	CharPtr pcFileName = NULL;	
	Char buf[255];
	Char fname[255];
	Int4 ilen = 0;
	CharPtr *ppcDirList = NULL;
	Int4 numtocompare, numleft,i,j;
	ValNodePtr vnp = NULL , vnpHead = NULL;
	Int4 iArray = 0;

	PDNML pdnmlTo = NULL;
        PMLD pmldTo  = NULL; 

    

        ErrSetLogfile("error_salign.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
	if (!GetArgs("salign - ASN.1 Structure Alignment Program. ",NUMARGS,Rargs))
		return 1;
	
    if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_ERROR,2,1,"Unable to open MMDBAPI");
		return 2;
	}
	atomtypes=ALIGN_ALLATOM;
	if (Rargs[7].intvalue==1)
			atomtypes=ALIGN_BACKBONE;
	if (Rargs[7].intvalue==2)
			atomtypes=ALIGN_CA;
	a1=Rargs[8].intvalue;
	b1=Rargs[9].intvalue;
	a2=Rargs[10].intvalue;
	b2=Rargs[11].intvalue;
	wsize=Rargs[6].intvalue;
	if (Rargs[13].strvalue == NULL) {   
		
		/* The usual salign - skip this if input file list specified */
		printf("Please note: Output file salign_out.val\ncorresponds to rotated/translated MMDB file 2\nMMDB file 1 is unmoved.\n");
	
		/* load an ASN.1 Biostruc */
		ModelNum1=Rargs[4].intvalue;
		pmsd1=LoadABiostruc(Rargs[0].strvalue,FALSE,ALLMDL,&ModelNum1);
		if (pmsd1==NULL) {
			/* error occurred */
			ErrPostEx(SEV_ERROR,3,1,"Unable to load biostruc 1");
			return 3;
		}
		ModelNum2=Rargs[5].intvalue;
		pmsd2=LoadABiostruc(Rargs[1].strvalue,FALSE,ALLMDL,&ModelNum2);
		if (pmsd2==NULL) {
			/* error occurred */
			ErrPostEx(SEV_ERROR,3,2,"Unable to load biostruc 2");
		return 4;
		}
		chain2=toupper((Rargs[3].strvalue)[0]);
		chain1=toupper((Rargs[2].strvalue)[0]);
		pdnmm=pmsd1->pdnmmHead;
		while (pdnmm) {
			pmmd1=(PMMD)(pdnmm->data.ptrvalue);
			if (chain1=='-')
				break;
			if (toupper((pmmd1->pcMolName)[0])==chain1)
				break;
			pdnmm=pdnmm->next;
		}
		if (!pdnmm) {
			/* error occurred */
			ErrPostEx(SEV_ERROR,4,1,"Unable to find chain %c in structure 1",chain1);
			return 5;
		}
		pdnmm=pmsd2->pdnmmHead;
		while (pdnmm) {
			pmmd2=(PMMD)(pdnmm->data.ptrvalue);
			if (chain2=='-')
				break;
			if (toupper((pmmd2->pcMolName)[0])==chain2)
				break;
			pdnmm=pdnmm->next;
		}
		if (!pdnmm) {
			/* error occurred */
			ErrPostEx(SEV_ERROR,4,2,"Unable to find chain %c in structure 2",chain2);
			return 6;
		}
	
		if (wsize>pmmd1->iResCount || wsize>pmmd2->iResCount) {
			ErrPostEx(SEV_FATAL,8,1,"Window size bigger than protein length!");
			return 7;
		}
		
		if (!wsize) wsize=pmmd1->iResCount;
	
		if (!a1 || !b1 || !a2 || !b2) {
			printf("Window Size: %d\n",wsize);
			printf("\nresidue\tRMSD (A)\n");
			for (pstart=1;pstart<(pmmd1->iResCount-wsize+2);pstart++) {
				if (!Align2StrucSVD(pmmd1,pmmd2,atomtypes,pstart,pstart+wsize-1,ModelNum1,ModelNum2)) {
					ErrPostEx(SEV_FATAL,9,1,"Structure alignment failed!");
					return 8;
				}			
				printf("%d\t%f\n",pstart+wsize/2,GetRMSD(pmmd1,pmmd2,atomtypes,pstart,pstart+wsize-1,ModelNum1,ModelNum2));
			}
		}
		else {
			/* in this mode, align only part of structure as specified */
			printf("Range1=%d-%d, Range2=%d-%d, M1=%d, M2=%d/n",a1,b1,a2,b2,ModelNum1,ModelNum2);	

			printf("\nresidue ranges\tRMSD (A)\n");	

			if (!Align2StrucSVDEx(pmmd1,pmmd2,atomtypes,a1,b1,a2,b2,ModelNum1,ModelNum2)) {
				ErrPostEx(SEV_FATAL,9,1,"Structure alignment failed!");
				return 9;
			}			
			printf("%d-%d,%d-%d\t%f\n",a1,b1,a2,b2,GetRMSDEx(pmmd1,pmmd2,atomtypes,a1,b1,a2,b2,ModelNum1,ModelNum2));
		}
		if (Rargs[12].intvalue == 0) {
			/* fast mode ignores output -just wants RMSD calculation */
		/*  PutMolBackWhereItCameFrom(pmmd1,ModelNum1); not necessary as we are not writing it out anymore */
			PutMolBackWhereItCameFrom(pmmd2,ModelNum2);
		/* write rotated molecules out to disk as ASN.1 binary */
		/* WriteAsnOneModel(pmsd1->pdnmsLink,ModelNum1,"align1.val",SAVE_BINARY); now unnecessary as equivalent to input file */ 

 		pdnmlTo = pmsd2->pdnmlModels;
		if (pdnmlTo != NULL) {
		      pmldTo = (PMLD)(pdnmlTo->data.ptrvalue);
			  if (pmldTo != NULL)
		         	if(pmldTo->ppAsnOrder) {
				       PTRVectorFree(pmldTo->ppAsnOrder,0);
				       pmldTo->ppAsnOrder = NULL;
				} 		
		}
	
			WriteAsnOneModel(pmsd2->pdnmsLink,ModelNum2,"salign_out.val",SAVE_BINARY); 
		if (Rargs[17].intvalue != 0) 
			WriteAsnOneModel(pmsd2->pdnmsLink,ModelNum2,"salign_out.prt",SAVE_ASCII); 		
		if (Rargs[16].intvalue != 0)    			
			if ((pFile=FileOpen("salign_out.pdb","w")) !=NULL) {    
         		   	WritePDBOneModel(pmsd2->pdnmsLink,  pFile,  ModelNum2);        
        	           	fflush(pFile);
         			FileClose(pFile);
		        }
       		}
	 
	/* Not supporting output of MIME Biostrucs... as most of the output is destined for strMerge, which does... */

	
	} /* Done the usual salign two-structure alignment */
	else
	{  /* we are doing a list of files  */
		
		if (Rargs[15].intvalue != 0){
        GenerateFileNamesFromDirectory("*.val");
		if(Rargs[12].intvalue == 0) printf("File created from local directory: %s\n",TEMP_FILE_NAME);
		vnpHead = GetFileNamesFromFile(TEMP_FILE_NAME);
		}
		else {
	/* load the file list */	
		vnpHead = GetFileNamesFromFile((CharPtr) Rargs[13].strvalue);
		}

        if (vnpHead == NULL) {
					/* error occurred */
					ErrPostEx(SEV_ERROR,3,1,"Bad list filename or list file is empty.");
					return 12;
		} 
		Misc_ValNodeListToCharPtrArray(vnpHead, &ppcDirList, &iArray);
		ValNodeFreeData(vnpHead);
        

		start_i = Rargs[14].intvalue;
		if (start_i == 0) start_i = 1;		

				/* load an ASN.1 Biostruc */
				ModelNum1=Rargs[4].intvalue;
				pmsd1=LoadABiostruc(ppcDirList[(start_i - 1)],FALSE,ALLMDL,&ModelNum1);
				if (pmsd1==NULL) {
					/* error occurred */
					ErrPostEx(SEV_ERROR,3,1,"Unable to load first biostruc in list");
					return 3;
				}
				/* set the usual salign stuff up */

				chain1=toupper((Rargs[2].strvalue)[0]);
				pdnmm=pmsd1->pdnmmHead;
				while (pdnmm) {
					pmmd1=(PMMD)(pdnmm->data.ptrvalue);
					if (chain1=='-')
						break;
					if (toupper((pmmd1->pcMolName)[0])==chain1)
						break;
					pdnmm=pdnmm->next;
				}
				if (!pdnmm) {
				/* error occurred */
					ErrPostEx(SEV_ERROR,4,1,"Unable to find chain %c in structure 1",chain1);
					return 5;
				}
				if (!wsize) wsize=pmmd1->iResCount;
	            if(Rargs[12].intvalue == 0) printf("Window Size: %d\n",wsize);
			    if(Rargs[12].intvalue == 0) printf("\nresidue\tRMSD (A)\n");
                ModelNum2=Rargs[5].intvalue;
				chain2=toupper((Rargs[3].strvalue)[0]);
	
				
			
    		/* Process structure 1 loop list */
	            for (i=start_i;i<numtocompare;i++){
					pmsd2=LoadABiostruc(ppcDirList[i],FALSE,ALLMDL,&ModelNum2);
				    if (pmsd2==NULL) {
					  ErrPostEx(SEV_ERROR,3,2,"Unable to load biostruc in list loop");
				      return 4;
				    }
					pdnmm=pmsd2->pdnmmHead;
					while (pdnmm) {
						pmmd2=(PMMD)(pdnmm->data.ptrvalue);
						if (chain2=='-')
							break;
						if (toupper((pmmd2->pcMolName)[0])==chain2)
							break;
						pdnmm=pdnmm->next;
					}
					if (!pdnmm) {
						/* error occurred */
						ErrPostEx(SEV_ERROR,4,2,"Unable to find chain %c in loop structure 2",chain2);
						return 6;
					}
	
					if (wsize>pmmd1->iResCount || wsize>pmmd2->iResCount) {
						ErrPostEx(SEV_FATAL,8,1,"Window size bigger than protein length!");
						return 7;
					}
		
					if (!wsize) wsize=pmmd1->iResCount;
	
					if (!a1 || !b1 || !a2 || !b2) {
						
						for (pstart=1;pstart<(pmmd1->iResCount-wsize+2);pstart++) {
							if (!Align2StrucSVD(pmmd1,pmmd2,atomtypes,pstart,pstart+wsize-1,ModelNum1,ModelNum2)) {
								ErrPostEx(SEV_FATAL,9,1,"Structure alignment failed!");
								return 8;
							}			
							printf("%d\t%f\n",pstart+wsize/2,GetRMSD(pmmd1,pmmd2,atomtypes,pstart,pstart+wsize-1,ModelNum1,ModelNum2));
						}
					}
					else {
						/* in this mode, align only part of structure as specified */
						if(Rargs[12].intvalue == 0) printf("Range1=%d-%d, Range2=%d-%d, M1=%d, M2=%d/n",a1,b1,a2,b2,ModelNum1,ModelNum2);	

						if(Rargs[12].intvalue == 0) printf("\nresidue ranges\tRMSD (A)\n");	

						if (!Align2StrucSVDEx(pmmd1,pmmd2,atomtypes,a1,b1,a2,b2,ModelNum1,ModelNum2)) {
							ErrPostEx(SEV_FATAL,9,1,"Structure alignment failed!");
							return 9;
						}			
						printf("%d-%d,%d-%d\t%f\n",a1,b1,a2,b2,GetRMSDEx(pmmd1,pmmd2,atomtypes,a1,b1,a2,b2,ModelNum1,ModelNum2));
					}
					if (Rargs[12].intvalue == 0) {
					/* fast mode ignores output -just wants RMSD calculation */
						PutMolBackWhereItCameFrom(pmmd2,ModelNum2);
					/* write rotated molecules out to disk as ASN.1 binary biostruc with no Bioseq for concatmodels */
						sprintf(fname,"saln_%d_%s",start_i,ppcDirList[i]);
						WriteAsnOneModel(pmsd2->pdnmsLink,ModelNum2,fname,SAVE_BINARY); 
					}
					
					if (i != numtocompare - start_i - 1 ) {
   					  FreeAModelstruc(pmsd2->pdnmsLink);
					} else /* keep it for next column */
                    pmsd3 = pmsd2;
					pmsd2 = NULL;
				} /* for */
               			
/*  DO 2ND COLUMN HERE  NBlast style two column process in NxN matrix */

	if (start_i > 1) {
                ClearStructures();
				pmsd1 = NULL;
				pdnmm = NULL;
                pmmd1 = NULL;
				pmmd2 = NULL;
			    
		        start_i = numtocompare - start_i;	
				pmsd1=pmsd3;
				if (pmsd1==NULL) {
					/* error occurred */
					ErrPostEx(SEV_ERROR,3,1,"Missed the starting biostruc in second column pass..");
					return 3;
				}
				/* set the usual salign stuff up */

				pdnmm=pmsd1->pdnmmHead;
				while (pdnmm) {
					pmmd1=(PMMD)(pdnmm->data.ptrvalue);
					if (chain1=='-')
						break;
					if (toupper((pmmd1->pcMolName)[0])==chain1)
						break;
					pdnmm=pdnmm->next;
				}
				if (!pdnmm) {
				/* error occurred */
					ErrPostEx(SEV_ERROR,4,1,"Unable to find chain %c in structure 1",chain1);
					return 5;
				}
			
    		/* Process second column loop list */
	            for (i=start_i;i<numtocompare;i++){
					pmsd2=LoadABiostruc(ppcDirList[i],FALSE,ALLMDL,&ModelNum2);
				    if (pmsd2==NULL) {
					  ErrPostEx(SEV_ERROR,3,2,"Unable to load biostruc in list loop");
				      return 4;
				    }
					pdnmm=pmsd2->pdnmmHead;
					while (pdnmm) {
						pmmd2=(PMMD)(pdnmm->data.ptrvalue);
						if (chain2=='-')
							break;
						if (toupper((pmmd2->pcMolName)[0])==chain2)
							break;
						pdnmm=pdnmm->next;
					}
					if (!pdnmm) {
						/* error occurred */
						ErrPostEx(SEV_ERROR,4,2,"Unable to find chain %c in loop structure 2",chain2);
						return 6;
					}
	
					if (wsize>pmmd1->iResCount || wsize>pmmd2->iResCount) {
						ErrPostEx(SEV_FATAL,8,1,"Window size bigger than protein length!");
						return 7;
					}
		
					if (!wsize) wsize=pmmd1->iResCount;
	
					if (!a1 || !b1 || !a2 || !b2) {
						
						for (pstart=1;pstart<(pmmd1->iResCount-wsize+2);pstart++) {
							if (!Align2StrucSVD(pmmd1,pmmd2,atomtypes,pstart,pstart+wsize-1,ModelNum1,ModelNum2)) {
								ErrPostEx(SEV_FATAL,9,1,"Structure alignment failed!");
								return 8;
							}			
							printf("%d\t%f\n",pstart+wsize/2,GetRMSD(pmmd1,pmmd2,atomtypes,pstart,pstart+wsize-1,ModelNum1,ModelNum2));
						}
					}
					else {
						/* in this mode, align only part of structure as specified */
						if(Rargs[12].intvalue == 0) printf("Range1=%d-%d, Range2=%d-%d, M1=%d, M2=%d/n",a1,b1,a2,b2,ModelNum1,ModelNum2);	

						if(Rargs[12].intvalue == 0) printf("\nresidue ranges\tRMSD (A)\n");	

						if (!Align2StrucSVDEx(pmmd1,pmmd2,atomtypes,a1,b1,a2,b2,ModelNum1,ModelNum2)) {
							ErrPostEx(SEV_FATAL,9,1,"Structure alignment failed!");
							return 9;
						}			
						printf("%d-%d,%d-%d\t%f\n",a1,b1,a2,b2,GetRMSDEx(pmmd1,pmmd2,atomtypes,a1,b1,a2,b2,ModelNum1,ModelNum2));
					}
					if (Rargs[12].intvalue == 0) {
					/* fast mode ignores output -just wants RMSD calculation */
						PutMolBackWhereItCameFrom(pmmd2,ModelNum2);
					/* write rotated molecules out to disk as ASN.1 binary biostruc with no Bioseq for concatmodels */
						sprintf(fname,"saln_%d_%s",start_i,ppcDirList[i]);

						/***do a better naming job because you have both filenames >>>>>>>>>>>>>>>>>>>>>>*/

						WriteAsnOneModel(pmsd2->pdnmsLink,ModelNum2,fname,SAVE_BINARY); 
					}
					
					FreeAModelstruc(pmsd2->pdnmsLink);
					pmsd2 = NULL;
				} /* for */
               			
	} /* if second column to pass */

          for(i = 0; i < iArray; i++) ppcDirList[i] = MemFree(ppcDirList[i]);
	 
       



	} /* list of files processing */
    
	/* Shut Down MMDB-API */
	CloseMMDBAPI();	
 	return 0;
}



/*
$Log: salign.c,v $
Revision 1.14  2009/01/19 16:06:24  chogue
salign modified to output only 2nd rotated/translated file, 1st file unmoved

Revision 1.13  2008/12/11 17:17:16  chogue
modified salign.c and rotate.c to move structures back where they came from for docking, rotate and translate entire structure - all molecules.

Revision 1.12  2008/12/11 15:20:34  chogue
removed comments from salign.c output blocking file writes for trades pkg

Revision 1.11  2004/07/16 15:07:13  hfeldman
Added additional error checking

Revision 1.10  2004/06/23 19:55:49  fwu
made 0 valid value for sequence range

Revision 1.9  2004/06/21 15:01:48  hfeldman
fixed typo

Revision 1.8  2004/06/17 16:40:27  hfeldman
Added chain arguments

Revision 1.7  2004/06/16 22:08:52  hfeldman
Added options for testing alignment of molecules of different length

Revision 1.6  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.5  2000/10/25 15:15:13  feldman
Made further updates, multiple model support is correct now
and relocated to a single function for loading correct model
and extensive error handling

Revision 1.4  2000/10/24 20:57:25  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.3  2000/07/14 20:09:25  feldman
Added parameter for MIMEBiostrucAsnGet for Bioseq

Revision 1.2  2000/07/06 15:29:37  feldman
-- Updated old makefiles to compile with newest toolkit and
directory structure
-- added needed functions to hfprogs.h
-- update some executables to include hfprogs.h instead of mmdbtraj.h
-- replaced all instances of BiostrucAsnGet with MIMEBiostrucAsnGet
which can read with MIME biostrucs (v2.0) or normal ones (v1.0), and
as a result, foldtraj uses a new skel.prt, to result in a MIME biostruc
as its initial input

Revision 1.1.1.1  2000/06/09 18:14:05  feldman
TraDES project

*/

