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

#ifdef _WIN32
int * __cdecl errno(void){return 0;}
#endif

#include "vis_traj.h"


#ifdef OS_UNIX_DARWIN
Char CFG_local_writepath[PATH_MAX];
Int2 Method;
#endif

/*****************************************************************
		         NCBI OGL DUMMY FUNCTIONS
*****************************************************************/

/*void Vtrj_OGL_DrawViewer3D(TOGL_Data * OGL_Data)
{
}*/

/*void Nlm_GetCurrentOGLFontMenuSettings(Nlm_Int2 *font, Nlm_Int2 *size,
    Nlm_Boolean *isBold, Nlm_Boolean *isItalic, Nlm_Boolean *isUnderlined)
{
} */

TOGL_Data * Cn3D_GetCurrentOGLData(void)
{
    return myOGL_data; /*Vtrj_ColorData.OGL_Data;*/
}

#ifdef WIN_MOTIF
void ** Cn3D_GetCurrentOGLDisplayHndl(void)
{
	return &(Vtrj_ColorData.OGL_Data->display);
}

void ** Cn3D_GetCurrentOGLVisinfoHndl(void)
{
	return &(Vtrj_ColorData.OGL_Data->visinfo);
}
#endif /* WIN_MOTIF */

Boolean timetoquit=FALSE;
Boolean CheckForAbort(void) {
	return FALSE;
}

/*****************************************************************
					WEBSITE FUNCTIONS
******************************************************************/

static void LaunchTradesSite(IteM i)
{
	Char       str [PATH_MAX];
#ifdef WIN_MOTIF
	NS_Window  window = NULL;
#endif

	sprintf (str, "http://www.blueprint.org/");
#ifdef WIN_MSWIN
	Nlm_MSWin_OpenDocument (str);
#endif
#ifdef WIN_MOTIF
	NS_OpenURL (&window, str, NULL, TRUE);
	NS_WindowFree (window);
#endif
}

/***************************************************************
										NEW TRAJGRAPH FUNCTIONS
***************************************************************/


static void GenerateModAAList(void)
{

	FILE *f;
	Char inbuf[MAXCOL], tmp_modAA[30],fnam[PATH_MAX];
	PEAS peas;

	tmp_modAA[0] = '\0';
	inbuf[0] = '\0';

	sprintf(fnam,"%s%s",CFG_local_datafilepath,EXTAA_FILENAME);
	if ((f=FileOpen(fnam,"r"))==NULL)
		ErrPostEx(SEV_ERROR,1,1,"Unable to open AA translation table %s",fnam);
	
	while (FileGets(inbuf, MAXCOL,f)!=NULL)
	{
		if (inbuf[0] == '[')
		{
			sscanf(inbuf, "%s", tmp_modAA);
			
			peas = GetExtAAInfo(tmp_modAA);

			if(peas->modlocation=='-') ValNodeAddPointer(&MmodAAlist, 0,(Nlm_VoidPtr) peas);
			if(peas->modlocation=='N') ValNodeAddPointer(&NmodAAlist, 0,(Nlm_VoidPtr) peas);
			if(peas->modlocation=='C') ValNodeAddPointer(&CmodAAlist, 0,(Nlm_VoidPtr) peas);

		}

	}
	FileClose(f);
}

/*--------------------------------------------------------------------------------------------------*/


static void DefaultVisMakeTrajVariables(void)
{

VMkTrj_valin[0] = '\0';
VMkTrj_trjout[0] = '\0';
VMkTrj_chain[0] = '\0';
VMkTrj_model = 1;
VMkTrj_startres = 1; 
VMkTrj_endres = 0;
VMkTrj_peakheight = 100; 
VMkTrj_savechi = 0;
VMkTrj_sdx = 0.0;
VMkTrj_sdy = 0.0;
VMkTrj_temp = 0.0;
VMkTrj_timestep = 0.0;
VMkTrj_quiteop = TRUE;
VMkTrj_noise = 0;
VMkTrj_method = 1;
StringCpy(VMkTrj_seqfile, "seq.in");
StringCpy(VMkTrj_ssout, "ssout");
VMkTrj_tgtype = 4;
VMkTrj_comprtype = 1; 
VMkTrj_constraints[0] = '\0';
VMkTrj_goroverride = FALSE;
VMkTrj_seqtype = 1;
VMkTrj_AAseq[0] = '\0';
VMkTrj_dorps = FALSE;
VMkTrj_douturn = FALSE;
VMkTrj_SSorTT = 1;
/*traj_quiet = (Byte)VERBOSITY_VERBOSE;*/
traj_quiet = (Byte)VERBOSITY_SILENT;

}

static void ClearSorVVariables(void)
{
	VMkTrj_comprtype = 1; 
	VMkTrj_constraints[0] = '\0';
	VMkTrj_trjout[0] = '\0';
}

/*------------------------------------------------------------------------*/

static void RefreshAAlist(void)
{
		Reset(newAAlist);
		isNew = TRUE;

		res_num=NewAAlen+1;
		if (res_num >= 1) numAA = res_num;
	
		LoadAAList(GetAAPos());
	
		isNew = FALSE;
	
		res_num=0;
		numAA=0;


}

static void GetAAseqlen(void)
{
	Char currentAAseq[MAXSIZE], testmod[20];
	static Char prevAAseq[MAXSIZE];
	CharPtr seqlenstring=NULL/*, nogood=NULL*/;
	Int2 currentAAseqlen=0, x=0, newAAseqlen=0, beginmod=0, posmod=0, posEmod;
	static Boolean isNewMod;
	Boolean /*isOk = TRUE,*/ shouldcopy = TRUE;
	PEAS peas=NULL, tmppeas=NULL;
	ValNodePtr tmpvalnode=NULL;

	isNewMod=FALSE;

	currentAAseq[0] = '\0';

	GetTitle(NewWiz2_text, currentAAseq, sizeof(currentAAseq));
	

	currentAAseqlen = StringLen(currentAAseq);
	if(currentAAseqlen==0)
	{
		SetTitle(seqlen_title, "0");
		NewAAlen=0;
		return;
	}


	for(x=0; x<=currentAAseqlen; x++)
	{
		if (currentAAseq[x] == '[')
		{
			posmod=x;
			isNewMod = TRUE;
			newAAseqlen++;
			shouldcopy = FALSE;
		}
		
		if (isNewMod==FALSE)
		{
			newAAseqlen++;
			shouldcopy = TRUE;
		}

		if(currentAAseq[x] ==  ']' && isNewMod==TRUE)
		{
			isNewMod = FALSE;
			testmod[beginmod] = ']';
			testmod[beginmod+1] = '\0';
			posEmod=x;
			peas = GetExtAAInfo(testmod);
			beginmod=0;

			if(peas==NULL) 
			{
/*					Message (MSG_OK, "The amino acid %s, is not a valid modified amino acid. \n Use the \"Insert A.A.\" button for assistance.", testmod);
					StringCpy(currentAAseq, prevAAseq);
					SetTitle(NewWiz2_text, currentAAseq);
					testmod[0] = '\0';
					return;
*/			}
			else
			{
				tmpvalnode = NmodAAlist;
				while(tmpvalnode!=NULL)
				{
					tmppeas = (PEAS) tmpvalnode->data.ptrvalue;

					if(!StringICmp(peas->keybname,tmppeas->keybname) && posmod!=0)
					{
						Message (MSG_OK, "You may not enter an amino acid before the N-terminus modified amino acid.");
						StringCpy(currentAAseq, prevAAseq);
						SetTitle(NewWiz2_text, currentAAseq);
						return;
					}
					tmpvalnode = tmpvalnode->next;
				}

				tmpvalnode = CmodAAlist;
				while(tmpvalnode!=NULL)
				{
					tmppeas = (PEAS) tmpvalnode->data.ptrvalue;

					if(!StringICmp(peas->keybname,tmppeas->keybname) && currentAAseq[posEmod+1]!='\0')
					{
						Message (MSG_OK, "You may not enter an amino acid after the C-terminus modified amino acid.");
						StringCpy(currentAAseq, prevAAseq);
						SetTitle(NewWiz2_text, currentAAseq);
						return;
					}
					tmpvalnode = tmpvalnode->next;
				}

			}

			shouldcopy = TRUE;
		}

		if(isNewMod==TRUE)
		{
			testmod[beginmod] = currentAAseq[x];
			beginmod++;
		}
	}

	NewAAlen=newAAseqlen-1;
	
	seqlenstring = Ultostr(NewAAlen,(Int2) NULL);
	SetTitle(seqlen_title, seqlenstring);
	if (shouldcopy == TRUE) StringCpy(prevAAseq, currentAAseq);
}

/*--------------------------------------------------------------------------------*/

static void DoneNewAAProc(IteM i)
{
	Hide(newAA_win);
	isAAused=FALSE;
}

static void EnabOkAAProc(IteM i)
{
	Int2 whichAA=0;

	whichAA = GetValue(newAAlist);

	SetTitle(newAA_info, descrAA[whichAA]);

	Enable(okAA_bttn);
}

static void InsertNewAAProc(IteM i)
{
	Int2 whichAA=0;
	Char AAinput[MAXSIZE];

	whichAA = (GetValue(newAAlist)-1);

	GetTitle(NewWiz2_text, AAinput, sizeof(AAinput));

	StringCat(AAinput, modAA[whichAA]);

	SetTitle(NewWiz2_text, AAinput);

	GetAAseqlen();

	if(isAAused==TRUE) RefreshAAlist();

	Enable(Wiz2Next_bttn);
}


static void DisplayAAList(void)
{

	GrouP g, Vtrj_grp1;
	
	isAAused=TRUE;

	newAA_win = FixedWindow(-50,-33,-10,-10, "Insert Amino Acid Selection", (WndActnProc) DoneNewAAProc);

	g = HiddenGroup(newAA_win,-1,-1, (GrpActnProc) NULL);
		
	SetGroupMargins(g, 17, 3);
	SetGroupSpacing(g, 3, 3);	
	
	newAAlist = SingleList (g,9,9, (LstActnProc) EnabOkAAProc);

	isNew = TRUE;
	
	Reset(newAAlist);

	res_num=NewAAlen+1;
	if (res_num >= 1) numAA = res_num+1;
	
	LoadAAList(GetAAPos());
	
	isNew = FALSE;
	
	res_num=0;
	numAA=0;

	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);
	
	SetGroupSpacing(Vtrj_grp1, 15, 3);

	okAA_bttn = PushButton(Vtrj_grp1, "Insert", (BtnActnProc) InsertNewAAProc);
	Disable(okAA_bttn);
	DefaultButton(Vtrj_grp1, "Close", (BtnActnProc) DoneNewAAProc);
	
	AlignObjects(ALIGN_CENTER, Vtrj_grp1, g, NULL);
	
	Break(g);
	
	newAA_info = ScrollText (g, 14, 3, systemFont, TRUE, NULL);

	AlignObjects(ALIGN_CENTER, newAA_info, g, NULL);

	Show(newAA_win);
	ProcessEvents();


}


static void SetTGType(IteM i)
{
	Int2 which_type=0;
	
	which_type =GetValue(TGtype_list);

	if(which_type==1) Disable(Gor_Group);
	if(which_type==2) Disable(Gor_Group);
	if(which_type==3) Enable(Gor_Group);
	if(which_type==4) Enable(Gor_Group);
	
}

static void SetSSorTT(IteM i)
{
	
	Char Vtrj_temp[10], Vtrj_timestep[10], Vtrj_sdx[10], Vtrj_sdy[10];

	Vtrj_temp[0] = '\0';
	Vtrj_timestep[0] = '\0';
	Vtrj_sdy[0] = '\0';
	Vtrj_sdx[0] = '\0';

	VMkTrj_SSorTT = GetValue(UseSSorTT);

	if (VMkTrj_SSorTT == 1)
	{
		Enable(SS_Group);
		Disable(TT_Group);
		VMkTrj_temp = 0.0;
		VMkTrj_timestep = 0.0;
		sprintf(Vtrj_temp, "%3.1f", VMkTrj_temp);
		sprintf(Vtrj_timestep, "%3.1f", VMkTrj_timestep);
		SetTitle(timestep_text, Vtrj_timestep);
		SetTitle(temp_text, Vtrj_temp);

	}
	if (VMkTrj_SSorTT ==2)
	{
		Enable(TT_Group);
		Disable(SS_Group);
		VMkTrj_sdx = 0.0;
		VMkTrj_sdy = 0.0;
		sprintf(Vtrj_sdx, "%3.1f", VMkTrj_sdx);
		sprintf(Vtrj_sdy, "%3.1f", VMkTrj_sdy);
		SetTitle(sdx_text, Vtrj_sdx);
		SetTitle(sdy_text, Vtrj_sdy);
	}

}


static CharPtr CreateSeqTmpFile(void)
{
	CharPtr path=NULL;
	FILE *f;	

	path=TmpNam(NULL);

	if ((f=FileOpen(path,"w"))==NULL)
	{
		ErrPostEx(SEV_ERROR,1,1,"Unable to create temporary sequence file %s", path);
		return NULL;
	}
	else
	{
		FilePuts(VMkTrj_AAseq,f);
		FileClose(f);
		return path;
	}
}

/*-------------------------------------------------------------*/

void MonitorMkTrj(void)
{	
	while (1) {	
#ifdef OS_MSWIN
		Sleep(1000);
#else
		sleep(1);
#endif
		if (ProgramProgress<=0)
			return;
		MonitorIntValue(Vtrj_mon, (Nlm_Int4)ProgramProgress);
	}	
}

void MonitorFoldTraj(void)
{	
	Boolean retval=TRUE;
	
	while (1) {	
#ifdef OS_MSWIN
		Sleep(1000);
#else
		sleep(1);
#endif
		if (ProgramProgress<=0)
			return;
		/* retval stores cancel value */
		retval=MonitorIntValue(Vtrj_mon, (Nlm_Int4)ProgramProgress);
		if (retval==FALSE) {
			/* cancelled by user */
			/* send signal to foldtraj */
			ProgramProgressMax=-1;
			return;
		}
	}	
}

static void DoneWiz4Proc(IteM i)
{
	Hide(NewWiz4_win);
	isNew=FALSE;
	lastg=0;
	needApply=FALSE;
	EnableDisableProc();
}


static void DoneNewWizProc(void)
{
	IteM i=NULL;
/*	Int2 len=0, x, y=0;
	Int2 clength=0;
	Char trjout_tmp[PATH_MAX];*/
	CharPtr fnamein=NULL, fnameout=NULL, chain=NULL, seqfile=NULL,
		constfile=NULL, strufile=NULL, hasext=NULL;
	CharPtr tmp_seqfile=NULL;
	pMakeTrjParamBlock mkdata=NULL;
	Boolean freeme=FALSE;

/*	clength=NewAAlen;*/

	mkdata = (pMakeTrjParamBlock)MemNew(sizeof(MakeTrjParamBlock));
	
	if (VMkTrj_seqtype==1 && VMkTrj_method == 1) {
			tmp_seqfile = CreateSeqTmpFile();
			StringCpy(VMkTrj_seqfile, tmp_seqfile);
	}

	if(VMkTrj_valin[0] != '\0')
	{
		hasext = StringStr (VMkTrj_valin, MMDB_EXT);
		if (hasext !=NULL) hasext[0] = '\0';

		fnamein =VMkTrj_valin;

		/*if (VMkTrj_trjout[0]=='\0') StringCpy(VMkTrj_trjout, VMkTrj_valin);*/
	}
	if(VMkTrj_seqfile[0] != '\0')
	{
		seqfile = VMkTrj_seqfile;
/*		if (VMkTrj_trjout[0]=='\0') StringCpy(VMkTrj_trjout, "protein");*/
	}
	if(VMkTrj_trjout[0] != '\0')
	{
		hasext = StringStr (VMkTrj_trjout, ASN_EXT);           	
		if (hasext !=NULL) hasext[0] = '\0';

		fnameout = VMkTrj_trjout;
	}

	if(VMkTrj_chain[0] != '\0') chain = VMkTrj_chain;
	if(VMkTrj_constraints[0] != '\0') constfile = VMkTrj_constraints;
	if(VMkTrj_ssout[0] != '\0') strufile=VMkTrj_ssout;

	
	DoneWiz4Proc(i);
	isNew=TRUE;
	EnableDisableProc();
		
	if (VMkTrj_method==1)
		mkdata->TrajMethod = 2;
	else if (VMkTrj_method==2)
		mkdata->TrajMethod = 1;
	mkdata->valfnamin = fnamein;
	mkdata->pcChain = chain;
	mkdata->modelnum = VMkTrj_model;
	mkdata->startres = VMkTrj_startres;
	mkdata->endres = VMkTrj_endres;
	mkdata->peakheight = VMkTrj_peakheight;
	mkdata->noise = VMkTrj_noise;
	mkdata->savechi = VMkTrj_savechi;
	mkdata->sigma_x = VMkTrj_sdx;
	mkdata->sigma_y = VMkTrj_sdy;
	mkdata->temperature = VMkTrj_temp;
	mkdata->timestep = VMkTrj_timestep;
	mkdata->seqfnamin = seqfile;
	mkdata->tgtype = VMkTrj_tgtype;
	mkdata->sstrufname = strufile;
	mkdata->DoRPS = VMkTrj_dorps;
	mkdata->uturn = VMkTrj_douturn;
	mkdata->SStruinput = VMkTrj_goroverride;
	mkdata->trjfnamout = fnameout;
	mkdata->comprtype = VMkTrj_comprtype;
	mkdata->constrfnamin = constfile;
	mkdata->zhangwindsize = 0;
	mkdata->ssmask=NULL;
	mkdata->errorfile=NULL;

	Disable(Vtrj_Exit);
	Disable(Vtrj_Help);
	Disable(Vtrj_About);
	Disable(Vtrj_WebSite);
	Disable(Vtrj_BgColor);

	if(NlmThreadsAvailable()==TRUE)
	{
		NlmThreadCreate(MakeTrj, (VoidPtr) mkdata);
		do {	
		} while (ProgramProgressMax<=0 && ProgramProgressMax!=-999);
		

		if (ProgramProgressMax==-999) {
			/* bail */
			NlmThreadJoinAll();
			ProgramProgress=0;
			ProgramProgressMax=0;
			isNew=FALSE;
	
			Enable(Vtrj_Exit);
			Enable(Vtrj_Help);
			Enable(Vtrj_About);
			Enable(Vtrj_WebSite);
			Enable(Vtrj_BgColor);
			EnableDisableProc();
			return;

		}
		Vtrj_mon = MonitorIntNewEx("Creating New Trajectory Distribution...", 0, (Nlm_Int4)ProgramProgressMax,FALSE);

		MonitorMkTrj();

		NlmThreadJoinAll();
		MonitorFree(Vtrj_mon);
	}
	else
	{
		Vtrj_mon = Nlm_MonitorStrNewEx("VisTraj", 50, FALSE);
		MonitorStrValue(Vtrj_mon, "Creating New Trajectory Distribution...");
		MakeTrj((VoidPtr)mkdata);
		MonitorFree(Vtrj_mon);
	}
	
	mkdata=MemFree(mkdata);
	ProgramProgress=0;
	ProgramProgressMax=0;
	
	Enable(Vtrj_Exit);
	Enable(Vtrj_Help);
	Enable(Vtrj_About);
	Enable(Vtrj_WebSite);
	Enable(Vtrj_BgColor);

	isNew = TRUE;
	if (fnameout==NULL) {
		fnameout=fnamein;
		if (fnameout==NULL) {
			freeme=TRUE;
			fnameout=StringSave(DEFAULT_TRJ_NAME);
		}
		else {
			/* maketrj adds .val */
			hasext = StringStr (fnameout, MMDB_EXT);           	
			if (hasext !=NULL) hasext[0] = '\0';	
		}
	}
	
	LoadProc(fnameout);

	Hide(aaname_panel);
	Show(aaname_panel);

	if (VMkTrj_trjout[0]=='\0' && !StringCmp(fnameout,DEFAULT_TRJ_NAME))
	{
		Message (MSG_OK, "It is advised that you change the file name of your new trajectory\ndistribution from the default name given.");
	}
	if (VMkTrj_seqtype==1 && VMkTrj_method == 1)
		FileRemove(VMkTrj_seqfile);
	if (freeme)
		fnameout=MemFree(fnameout);
	DefaultVisMakeTrajVariables();

	isNew = FALSE;
  NewAAlen = 0;
}

static void DoneWiz3VProc(IteM i)
{
	Hide(NewWiz3V_win);
	isNew=FALSE;
	EnableDisableProc();
}

static void DoneWiz3SProc(IteM i)
{
	Hide(NewWiz3S_win);
	isNew=FALSE;
	EnableDisableProc();
}

static void DoneWiz2Proc(IteM i)
{
	Hide(NewWiz2_win);
	Hide(newAA_win);
	isNew=FALSE;
	EnableDisableProc();
}

static void DoneWiz1Proc(IteM i)
{
	Hide(NewWiz1_win);
	isNew=FALSE;
	EnableDisableProc();
}

/*------------------------------------------------------------*/

static void BackWiz2Proc(IteM i)
{
	DoneWiz2Proc(i);
	DisplayNewWiz1();
}

static void BackWiz3VProc(IteM i)
{
	DoneWiz3VProc(i);
	DisplayNewWiz2(i);
}


static void BackWiz3SProc(IteM i)
{
	DoneWiz3SProc(i);
	DisplayNewWiz2(i);
}

static void BackWiz4Proc(IteM i)
{

	DoneWiz4Proc(i);
	if(VMkTrj_method==2) DisplayNewWiz3Val(i);
	if(VMkTrj_method==1) DisplayNewWiz3Seq(i);

}

/*------------------------------------------------------------*/

static Boolean CheckNewAASequence(void)
{

	Char tmp_seq[MAXSIZE], new_seq[MAXSIZE], tmp_mod[MAX_KBAA_NAMELENGTH+1];
	Int2 len=0, x, y, z, n=0, newlen=0, Vtrj_ans=ANS_NO;
	PEAS peas=NULL;
	Boolean ignore=FALSE, isMod=FALSE;

	GetTitle(NewWiz2_text, tmp_seq, sizeof(tmp_seq));

	len = StringLen(tmp_seq);
	
	for(x=0; x<=len; x++)
	{
		if(tmp_seq[x]=='[')
		{
			isMod = TRUE;
			z=0;
			y=x;
			while(tmp_seq[y]!=']')
			{
				tmp_mod[z] = tmp_seq[y]; 
				y++;
				z++;
			}			
			tmp_mod[z]= ']';
			tmp_mod[z+1] = '\0';
			peas = GetExtAAInfo(tmp_mod);
			if(peas==NULL) ignore= TRUE;
			else ignore=FALSE;
		}

		else if(tmp_seq[x]==']') isMod=FALSE;

		else if(isMod==FALSE && tmp_seq[x] != ']')
		{
			if(isAA(tmp_seq[x])==0) ignore=TRUE;
			else ignore=FALSE;
		}

		if(ignore==FALSE)
		{
			new_seq[n] = tmp_seq[x];
			n++;
		}
	}
	new_seq[n] = '\0';

	newlen = StringLen(new_seq);

	if (len!=newlen)
	{
		Vtrj_ans = Message (MSG_YN, "The sequence you have entered has been altered to conform to the standard amino acid naming convension. Would you like to continue?");	
		if(Vtrj_ans==ANS_YES) return TRUE;
		
		if(Vtrj_ans==ANS_NO)
		{
			SetTitle(NewWiz2_text, new_seq);
			Select(NewWiz2_text);
			GetAAseqlen();
			return FALSE;
		}
	}
	/*else*/ return TRUE;

}

static void CheckWiz2Info(IteM i)
{
	Char wiz2input[MAXSIZE];
	Boolean isOk=FALSE;

	wiz2input[0] = '\0';

	GetTitle(NewWiz2_text, wiz2input, sizeof(wiz2input));
	

	if(VMkTrj_method == 2)
	{
		if(FileLengthEx(wiz2input) != -1)
		{
			StringCpy(VMkTrj_valin, wiz2input);
			isOk=TRUE;
		}
		else
		{
			Message (MSG_OK, "The file %s does not exist.", wiz2input);
			SetTitle(NewWiz2_text, "");
			isOk=FALSE;
			return;
		}

	}
	if(VMkTrj_method == 1)
	{
		if(VMkTrj_seqtype==2)
		{
			if(FileLengthEx(wiz2input) != -1)
			{
				StringCpy(VMkTrj_seqfile, wiz2input);
				isOk=TRUE;
			}
			else
			{

				Message (MSG_OK, "The file %s, does not exist.", wiz2input);
				SetTitle(NewWiz2_text, "");
				isOk=FALSE;
				return;
			}
		}
		if(VMkTrj_seqtype==1)
		{
			isOk = CheckNewAASequence();
			if (isOk==FALSE) return;
			StringCpy(VMkTrj_AAseq, wiz2input);
		}
	}

	if (isOk==TRUE)
	{
		DoneWiz2Proc(i);
		if(VMkTrj_method == 1) DisplayNewWiz3Seq(i);
		if(VMkTrj_method == 2) DisplayNewWiz3Val(i);
	}

}

static void CheckWiz3VInfo(IteM i)
{

	Char tmpinput[MAXSIZE];	
	int tmpint=0,err;
	float tmpfloat=0.0;

	tmpinput[0] = '\0';

	GetTitle(chain_text, tmpinput, sizeof(tmpinput));
	StringCpy(VMkTrj_chain, tmpinput);

	GetTitle(model_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%d", &tmpint);
	if (tmpint<1 || tmpint>9999 || err==0)
	{
		Message(MSG_OK, "Model number must be from 1 to 9999");
		Select(model_text);
		return;
	}
	VMkTrj_model = (Int2)tmpint;

	GetTitle(startres_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%d", &tmpint);
	if (tmpint<1 || tmpint>9999 || err==0)
	{
		Message(MSG_OK, "Start Residue number must be from 1 to 9999");
		Select(startres_text);
		return;
	}
	VMkTrj_startres = (Int2)tmpint;

	GetTitle(endres_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%d", &tmpint);
	if (tmpint<0 || tmpint>9999 || err==0)
	{
		Message(MSG_OK, "Finish Residue number must be from 0 to 9999");
		Select(endres_text);
		return;
	}
	VMkTrj_endres = (Int2)tmpint;

	if (VMkTrj_endres-VMkTrj_startres<2 && VMkTrj_endres!=0)
	{
		Message(MSG_OK, "Sequence must be at least three residues long");
		Select(endres_text);
		return;
	}

	GetTitle(peakheight_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%d", &tmpint);
	if (tmpint<1 || tmpint>999999 || err==0)
	{
		Message(MSG_OK, "Gaussian Peak Height must be from 1 to 999999");
		Select(peakheight_text);
		return;
	}
	VMkTrj_peakheight = (Int2)tmpint;

	GetTitle(sdx_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%f", &tmpfloat);
	if (tmpfloat<0.0 || tmpfloat>45.0 || err==0)
	{
		Message(MSG_OK, "X standard deviation must be between 0 and 45 degrees");
		Select(sdx_text);
		return;
	}
	VMkTrj_sdx = (FloatLo)tmpfloat;

	GetTitle(sdy_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%f", &tmpfloat);
	if (tmpfloat<0.0 || tmpfloat>22.0 || err==0)
	{
		Message(MSG_OK, "Y standard deviation must be between 0 and 22 degrees");
		Select(sdy_text);
		return;
	}
	VMkTrj_sdy = (FloatLo)tmpfloat;

	GetTitle(temp_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%f", &tmpfloat);
	if (tmpfloat<0.0 || tmpfloat>9999.0 || err==0)
	{
		Message(MSG_OK, "Temperature must be between 0K and 9999K");
		Select(temp_text);
		return;
	}
	VMkTrj_temp = (FloatLo)tmpfloat;

	GetTitle(timestep_text, tmpinput, sizeof(tmpinput));
	err=sscanf(tmpinput, "%f", &tmpfloat);
	if (tmpfloat<0.0 || tmpfloat>9999.0 || err==0)
	{
		Message(MSG_OK, "Timestep must be between 0 fs and 9999 fs");
		Select(timestep_text);
		return;
	}
	VMkTrj_timestep = (FloatLo)tmpfloat;

	VMkTrj_comprtype = (GetValue(comptypeV_list)-1);
	VMkTrj_savechi = (GetValue(savechiV_list)-1);
	
	GetTitle(NewWiz3V_text, tmpinput, sizeof(tmpinput));
	StringCpy(VMkTrj_constraints, tmpinput);

	GetTitle(NewWiz3V2_text, tmpinput, sizeof(tmpinput));
	StringCpy(VMkTrj_trjout, tmpinput);

/*	tmpint = GetValue(VCallRPS);
	if (tmpint==1) VMkTrj_dorps = TRUE;
	if (tmpint==2) VMkTrj_dorps = FALSE;
*/

	DoneWiz3VProc(i);
	DisplayNewWiz4(i);

}

static void CheckWiz3SInfo(IteM i)
{
	Char tmpinput[PATH_MAX];	
	Int2 tmpint=0;

	tmpinput[0] = '\0';


	GetTitle(NewWiz3S3_text, tmpinput, sizeof(tmpinput));
	StringCpy(VMkTrj_ssout, tmpinput);

	GetTitle(NewWiz3S_text, tmpinput, sizeof(tmpinput));
	StringCpy(VMkTrj_constraints, tmpinput);

	GetTitle(NewWiz3S2_text, tmpinput, sizeof(tmpinput));
	StringCpy(VMkTrj_trjout, tmpinput);

	VMkTrj_comprtype = (GetValue(comptypeS_list)-1);
	VMkTrj_tgtype = GetValue(TGtype_list);
	
	tmpint = GetValue(UseSSFile);
	if (tmpint==1) VMkTrj_goroverride = TRUE;
	if (tmpint==2) VMkTrj_goroverride = FALSE;

	tmpint = GetValue(SCallRPS);
	if (tmpint==1) VMkTrj_dorps = TRUE;
	if (tmpint==2) VMkTrj_dorps = FALSE;

	tmpint = GetValue(SCallUturn);
	if (tmpint==1) VMkTrj_douturn = TRUE;
	if (tmpint==2) VMkTrj_douturn = FALSE;

	DoneWiz3SProc(i);
	DisplayNewWiz4(i);

}


static void CheckWiz4Info(IteM i)
{
	Int2 ans;
	
	if (needApply==TRUE) {
		ans=Message (MSG_YN, "Apply unsaved changes?");
		if (ans==ANS_YES) {
			ApplyGNewEditProc(i);
			/* check if apply failed */
			if (needApply==TRUE)
				return;		
		}
	}
	DoneNewWizProc();
}


/*------------------------------------------------------------*/


static void EnableWiz2Next(IteM i)
{
	Enable(Wiz2Next_bttn);

	if(VMkTrj_method==1 && VMkTrj_seqtype==1) GetAAseqlen();
	if(isAAused==TRUE) RefreshAAlist();
}

/*------------------------------------------------------------*/

static void DoWiz3BrowseProc(Int2 VorS)
{
	
	Char path[PATH_MAX];

	path[0] = '\0';
	
	
	if (GetInputFileName(path, sizeof(path), "*", NULL))
	{
		if (VorS==1) SetTitle(NewWiz3V_text, path);
		if (VorS==2) SetTitle(NewWiz3S_text, path);
		
		if (VorS==3) 
		{
			SetTitle(NewWiz3V2_text, path);
			Enable(Wiz3VNext_bttn);
		}
		
		if (VorS==4) 
		{
			SetTitle(NewWiz3S2_text, path);
			Enable(Wiz3SNext_bttn);
		}

		if (VorS==5) SetTitle(NewWiz3S3_text, path);
	}
}

static void Wiz3VBrowseProc(IteM i)
{
	DoWiz3BrowseProc(1);
	return;
}

static void Wiz3V2BrowseProc(IteM i)
{
	DoWiz3BrowseProc(3);
	return;
}

static void Wiz3SBrowseProc(IteM i)
{
	DoWiz3BrowseProc(2);
	return;
}

static void Wiz3S2BrowseProc(IteM i)
{
	DoWiz3BrowseProc(4);
	return;
}

static void Wiz3S3BrowseProc(IteM i)
{
	DoWiz3BrowseProc(5);
	return;
}

/*------------------------------------------------------------*/


static void DoWiz2BrowseProc(void)
{
	
	IteM i=NULL;
	Char ext[4];
	Char path[PATH_MAX];

	path[0] = '\0';
	ext[0] = '\0';


	if (VMkTrj_method==1) StringCpy(ext, "*");
	if (VMkTrj_method==2) StringCpy(ext, "val");
		
	
	if (GetInputFileName(path, sizeof(path), ext, NULL))
	{
		SetTitle(NewWiz2_text, path);
		EnableWiz2Next(i);
	}
}

static void Wiz2BrowseProc(IteM i)
{
	if(VMkTrj_method==1 && VMkTrj_seqtype==1)
	{
		DisplayAAList();
		return;
	}

	DoWiz2BrowseProc();
	return;
}

/*------------------------------------------------------------*/

static void LoadGNewParams(void)
{
	Char value[PATH_MAX];
	int tmp_int=0;
	float tmp_flt=0.0;

	if (NewGlobs!=NULL) NewGlobs = MemFree(NewGlobs);
	NewGlobs = (VNGPtr) MemNew(sizeof(VNG));

	GetAppParam(TASK_CFGFILE,"INITTRAJ","WALKTYPE","",value,PATH_MAX);
	StringCpy(NewGlobs->wt, value);

	GetAppParam(TASK_CFGFILE,"INITTRAJ","TIMEOUT","",value,PATH_MAX);
	sscanf(value, "%d", &tmp_int);
	NewGlobs->to = (Int4) tmp_int;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","BACKBONE_ERROR_TOLERANCE","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->bbet = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","BACKBONE_PRECISION","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->bbp = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","ATOM_BOUNCINESS_BB","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->abbb = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","ATOM_BOUNCINESS_SC","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->absc = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","BUMPCHECK_HYDROGEN","",value,PATH_MAX);
	StringCpy(NewGlobs->bch, value);

	GetAppParam(TASK_CFGFILE,"INITTRAJ","INCSIZE","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->is = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","TUNNEL_PROB","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->tp = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","START_BACKBONE","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->sbb = (FloatHi) tmp_flt;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","NUM_ROT_TRIES","",value,PATH_MAX);
	sscanf(value, "%d", &tmp_int);
	NewGlobs->nrt = (Int2) tmp_int;

	GetAppParam(TASK_CFGFILE,"INITTRAJ","MARKOV_SCALE_FACTOR","",value,PATH_MAX);
	sscanf(value, "%f", &tmp_flt);
	NewGlobs->msf = (FloatHi) tmp_flt;
}

static void LoadGNewEditProc(IteM i)
{
	Int2 which_newG=0;
	Char printwhat[100];


	which_newG = GetValue(editGNewlist);

	
	switch (which_newG)
	{
		case 1:	
			sprintf(printwhat, "%s", NewGlobs->wt);
			break;
		case 2:	
			sprintf(printwhat, "%1.1f", NewGlobs->bbet);
			break;
		case 3:	
			sprintf(printwhat, "%1.4f", NewGlobs->bbp);
			break;
		case 4:	
			sprintf(printwhat, "%d", NewGlobs->nrt);
			break;
		case 5:	
			sprintf(printwhat, "%1.2f", NewGlobs->abbb);
			break;
		case 6:	
			sprintf(printwhat, "%1.2f", NewGlobs->absc);
			break;
		case 7:	
			sprintf(printwhat, "%s", NewGlobs->bch);
			break;	
		case 8:	
			sprintf(printwhat, "%1.1f", NewGlobs->is);
			break;
		case 9:	
			sprintf(printwhat, "%1.1f", NewGlobs->sbb);
			break;
		case 10:	
			sprintf(printwhat, "%d", NewGlobs->to);
			break;
		case 11:	
			sprintf(printwhat, "%1.2f", NewGlobs->msf);
			break;
		case 12:	
			sprintf(printwhat, "%1.2f", NewGlobs->tp);
			break;
	}
	
	
	if(which_newG==1)
	{
		if(printwhat[0]=='P') SetValue(Newwalktype_popup, 1);
		if(printwhat[0]=='C') SetValue(Newwalktype_popup, 2);
		Hide(currentGNewvalue_text);
		Hide(Newbchydrogen_popup);
		Show(Newwalktype_popup);
	}
	else if(which_newG==7)
	{
		if(printwhat[0]=='T') SetValue(Newbchydrogen_popup, 1);
		if(printwhat[0]=='F') SetValue(Newbchydrogen_popup, 2);
		Hide(currentGNewvalue_text);
		Hide(Newwalktype_popup);
		Show(Newbchydrogen_popup);
	}
	else
	{
		Hide(Newbchydrogen_popup);
		Hide(Newwalktype_popup);
		Show(currentGNewvalue_text);
		if (lastg!=which_newG) {
			SetTitle(currentGNewvalue_text, printwhat);
			lastg=which_newG;
		}
	}

	LoadInfoTextProc(which_newG, 5);
/*	Disable(NewGApply_bttn);*/
}

/*------------------------------------------------------------*/

static void EnableGNewSaveEditProc(IteM i)
{
	Int2 which_newG=0, whichone=0;
	int itemp=0;
	float ftemp=0.0;
	Char str[30];

	
	which_newG = GetValue(editGNewlist);
	
	GetTitle(currentGNewvalue_text, str, sizeof(str));

/*	Enable(editGNewlist);*/
	switch (which_newG)
	{
	
		case 1:	
			whichone = GetValue(Newwalktype_popup);
			if(whichone==1) StringCpy(NewGlobs->wt, "PHIPSI");
			if(whichone==2) StringCpy(NewGlobs->wt, "CA");
			break;
		case 2:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->bbet = (FloatHi) ftemp;
			break;
		case 3:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->bbp = (FloatHi) ftemp;
			break;
		case 4:	
			sscanf(str, "%d", &itemp);
			NewGlobs->nrt = (Int2) itemp;
			break;
		case 5:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->abbb = (FloatHi) ftemp;
			break;
		case 6:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->absc = (FloatHi) ftemp;
			break;
		case 7:	
			whichone = GetValue(Newbchydrogen_popup);
			if(whichone==1) StringCpy(NewGlobs->bch, "TRUE");
			if(whichone==2) StringCpy(NewGlobs->bch, "FALSE");
			break;	
		case 8:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->is = (FloatHi) ftemp;
			break;
		case 9:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->sbb = (FloatHi) ftemp;
			break;
		case 10:	
			sscanf(str, "%d", &itemp);
			NewGlobs->to = (Int4) itemp;
			break;
		case 11:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->msf = (FloatHi) ftemp;
			break;
		case 12:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->tp = (FloatHi) ftemp;
			break;
	}
	
	LoadGNewEditProc(i);
	
	Enable(NewGApply_bttn);
	needApply=TRUE;
	/*Disable(editGNewlist);*/
}

static void ApplyGNewEditProc(IteM i)
{
	Int2 whichone=0;
	int itemp=0;
	float ftemp=0.0;
	Char str[30];

	
/*	Enable(editGNewlist);*/
	ftemp=NewGlobs->bbet;
	if (ftemp<10.0) {
		Message (MSG_OK, "Backbone Error Tolerance must be 10 or greater");
		SetValue(editGNewlist,2);
		LoadGNewEditProc(i);
		return;
	}	
	ftemp=NewGlobs->bbp;
	if (ftemp<0.0005 || ftemp>5.0) {
		Message (MSG_OK, "Backbone Precision must be between 0.0005 and 5");
		SetValue(editGNewlist,3);
		LoadGNewEditProc(i);
		return;
	}
	itemp=NewGlobs->nrt;
	if (itemp<1 || itemp>100) {
		Message(MSG_OK,"Rotamer Tries must be between 1 and 100");
		SetValue(editGNewlist,4);
		LoadGNewEditProc(i);
		return;
	}
	ftemp=NewGlobs->abbb;
	if (ftemp<0.0) {
		Message(MSG_OK,"Atom Bounciness must be non-negative");
		SetValue(editGNewlist,5);
		LoadGNewEditProc(i);
		return;
	}
	ftemp=NewGlobs->absc;
	if (ftemp<0.0) {
		Message(MSG_OK,"Atom Bounciness must be non-negative");
		SetValue(editGNewlist,6);
		LoadGNewEditProc(i);
		return;
	}
	ftemp=NewGlobs->is;
	if (ftemp<1.0 || ftemp>89.0) {
		Message(MSG_OK,"Increment Size must be between 1 and 89");
		SetValue(editGNewlist,8);
		LoadGNewEditProc(i);
		return;
	}
	ftemp=NewGlobs->tp;
	if (ftemp<0.0 || ftemp>1.0) {
		Message(MSG_OK,"Tunnelling Probability must be between 0.0 and 1.0");
		SetValue(editGNewlist,12);
		LoadGNewEditProc(i);
		return;
	}
	ftemp=NewGlobs->sbb;
	if (ftemp<=0.0 || ftemp>360.0) {
		Message(MSG_OK,"Backbone Start must be greater than 0 and less than or equal to 360");
		SetValue(editGNewlist,9);
		LoadGNewEditProc(i);
		return;
	}
	itemp=NewGlobs->to;
	if (itemp<BACKTRACK_TRIES_MIN || itemp>BACKTRACK_TRIES_MAX) {
		Message(MSG_OK,"Mean Tries per Residue must be between %s and %s inclusive",BACKTRACK_TRIES_MIN,BACKTRACK_TRIES_MAX);
		SetValue(editGNewlist,10);
		LoadGNewEditProc(i);
		return;
	}
	ftemp=NewGlobs->msf;
	if (ftemp<0.0 || ftemp>1.0) {
		Message(MSG_OK,"Markov Scale Factor must be between zero and one");
		SetValue(editGNewlist,11);
		LoadGNewEditProc(i);
		return;
	}

	whichone = GetValue(Newwalktype_popup);
	if(whichone==1) StringCpy(NewGlobs->wt, "PHIPSI");
	if(whichone==2) StringCpy(NewGlobs->wt, "CA");
	StringCpy(str, NewGlobs->wt);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","WALKTYPE",str);
	sprintf(str, "%f",   NewGlobs->bbet);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","BACKBONE_ERROR_TOLERANCE",str);
	sprintf(str, "%f",   NewGlobs->bbp);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","BACKBONE_PRECISION",str);
	sprintf(str, "%d",   NewGlobs->nrt);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","NUM_ROT_TRIES",str);
	sprintf(str, "%f",   NewGlobs->abbb);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","ATOM_BOUNCINESS_BB",str);
	sprintf(str, "%f",   NewGlobs->absc);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","ATOM_BOUNCINESS_SC",str);
	whichone = GetValue(Newbchydrogen_popup);
	if(whichone==1) StringCpy(NewGlobs->bch, "TRUE");
	if(whichone==2) StringCpy(NewGlobs->bch, "FALSE");
	StringCpy(str, NewGlobs->bch);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","BUMPCHECK_HYDROGEN",str);
	sprintf(str, "%f",   NewGlobs->is);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","INCSIZE",str);
	sprintf(str, "%f",   NewGlobs->tp);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","TUNNEL_PROB",str);
	sprintf(str, "%f",   NewGlobs->sbb);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","START_BACKBONE",str);
	sprintf(str, "%d",   NewGlobs->to);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","TIMEOUT",str);
	sprintf(str, "%f", NewGlobs->msf);
	SetAppParam(TASK_CFGFILE,"INITTRAJ","MARKOV_SCALE_FACTOR",str);
	
	needApply=FALSE;
	Disable(NewGApply_bttn);
	LoadGNewEditProc(i);
}

static void DisplayNewWiz4(IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3;
	
	LoadGNewParams();

	isNew = TRUE;
	EnableDisableProc();
	needApply=FALSE;

	NewWiz4_win = MovableModalWindow(-50,-33,-15,-10, "New Trajectory Distribution Wizard - Step 4", (WndActnProc) DoneWiz4Proc);

	g = HiddenGroup(NewWiz4_win,-1,-1, (GrpActnProc) NULL);
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 10, 10);			
	


	Vtrj_grp1 =  HiddenGroup(g,-1,-1,(GrpActnProc) NULL);
	
	editGNewlist = SingleList (Vtrj_grp1,16,8, (LstActnProc) LoadGNewEditProc);
		
	ListItem (editGNewlist,"Trajectory Space");
	ListItem (editGNewlist,"Backbone Error Tolerance");
	ListItem (editGNewlist,"Backbone Precision");
	ListItem (editGNewlist,"Rotamer Tries");
	ListItem (editGNewlist,"Backbone Atom Bounciness");
	ListItem (editGNewlist,"Sidechain Atom Bounciness");
	ListItem (editGNewlist,"Bumpcheck Hydrogen");
	ListItem (editGNewlist,"Increment Size");
	ListItem (editGNewlist,"Backbone Start");
	ListItem (editGNewlist,"Mean Tries per Residue");
	ListItem (editGNewlist,"Markov Scale Factor");
	ListItem (editGNewlist,"Tunnelling Probability");
		
	Vtrj_grp3 = HiddenGroup(Vtrj_grp1, 0, 0, (GrpActnProc) NULL);

	currentGNewvalue_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) EnableGNewSaveEditProc);
	Newwalktype_popup = PopupList(Vtrj_grp3, TRUE, (PupActnProc) EnableGNewSaveEditProc);
	Newbchydrogen_popup = PopupList(Vtrj_grp3, TRUE, (PupActnProc) EnableGNewSaveEditProc);

	PopupItem(Newwalktype_popup, "RAMACHANDRAN");
	PopupItem(Newwalktype_popup, "ALPHA_CARBON");	

	PopupItem(Newbchydrogen_popup, "TRUE");
	PopupItem(Newbchydrogen_popup, "FALSE");

	Hide(Newwalktype_popup);
	Hide(Newbchydrogen_popup);

	Select(currentGNewvalue_text);

	Vtrj_grp2 = HiddenGroup(g, 4, 0, (GrpActnProc) NULL);
	NewGApply_bttn = PushButton(Vtrj_grp2, "Apply", (BtnActnProc) ApplyGNewEditProc);
	Disable(NewGApply_bttn);

	PushButton(Vtrj_grp2, "Cancel", (BtnActnProc) DoneWiz4Proc);
	PushButton(Vtrj_grp2, "<< Back", (BtnActnProc) BackWiz4Proc);
	Wiz3VNext_bttn = PushButton(Vtrj_grp2, "Begin", (BtnActnProc) CheckWiz4Info);

	Break(g);
	
	GNewEdit_Info = ScrollText (g, 24, 8, systemFont, TRUE, NULL);

	AlignObjects(ALIGN_JUSTIFY, GNewEdit_Info, editGNewlist, NULL);
	
	AlignObjects(ALIGN_JUSTIFY, currentGNewvalue_text, editGNewlist, NULL);
	AlignObjects(ALIGN_JUSTIFY, Newwalktype_popup, editGNewlist, NULL);
	AlignObjects(ALIGN_JUSTIFY, Newbchydrogen_popup, editGNewlist, NULL);
	
	AlignObjects(ALIGN_CENTER, Vtrj_grp2, GNewEdit_Info, NULL);

	Show(NewWiz4_win);
	ProcessEvents();

}

static void DisplayNewWiz3Val(IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3, Vtrj_grp5,
		Vtrj_grp6, Vtrj_grp7, Vtrj_grp8, Vtrj_grp9, Vtrj_grp11;

	Char Vtrj_ChainName[100], Vtrj_sdx[10], Vtrj_sdy[10], Vtrj_temp[10], Vtrj_timestep[10],
		Vtrj_model[10], Vtrj_endres[10], Vtrj_peakheight[10], Vtrj_startres[10];

	isNew=TRUE;
	EnableDisableProc();

	Vtrj_ChainName[0] = '\0';
	Vtrj_sdx[0] = '\0';
	Vtrj_sdy[0] = '\0';
	Vtrj_temp[0] = '\0';
	Vtrj_timestep[0] = '\0';
	Vtrj_model[0] = '\0';
	Vtrj_startres[0] = '\0';
	Vtrj_endres[0] = '\0';
	Vtrj_peakheight[0] = '\0';


	sprintf(Vtrj_endres, "%d", VMkTrj_endres);
	sprintf(Vtrj_sdx, "%3.1f", VMkTrj_sdx);
	sprintf(Vtrj_sdy, "%3.1f", VMkTrj_sdy);
	sprintf(Vtrj_temp, "%3.1f", VMkTrj_temp);
	sprintf(Vtrj_timestep, "%3.1f", VMkTrj_timestep);
	sprintf(Vtrj_model, "%d", VMkTrj_model);
	sprintf(Vtrj_peakheight, "%d", VMkTrj_peakheight);
	sprintf(Vtrj_startres, "%d", VMkTrj_startres);

	StringCpy(Vtrj_ChainName, VMkTrj_chain);

	

	NewWiz3V_win = MovableModalWindow(-50,-33,-10,-10, "New Trajectory Distribution Wizard - Step 3 (MMDB Structure File)", (WndActnProc) DoneWiz3VProc);

	g = HiddenGroup(NewWiz3V_win,0,1, (GrpActnProc) NULL);	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 15, 15);

	Vtrj_grp1 = NormalGroup(g, 2, 0, "General",systemFont, (GrpActnProc) NULL);
/*	StaticPrompt (Vtrj_grp1, "            ", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp1, "(Default: First One Found)", 0, dialogTextHeight, systemFont, 'l');*/
	StaticPrompt (Vtrj_grp1, "Chain Name (Default: 1st One Found)", 0, dialogTextHeight, systemFont, 'l');
	chain_text = DialogText (Vtrj_grp1, Vtrj_ChainName,  10, (TxtActnProc) NULL);
	StaticPrompt (Vtrj_grp1, "Model Number", 0, dialogTextHeight, systemFont, 'l');
	model_text = DialogText (Vtrj_grp1, Vtrj_model,  5, (TxtActnProc) NULL);	
/*	StaticPrompt (Vtrj_grp1, "Call RPSBlast? (requires CDD database)", 0, dialogTextHeight, systemFont, 'l');

	VCallRPS = HiddenGroup(Vtrj_grp1, 2, 0, (GrpActnProc) NULL);
	RadioButton (VCallRPS, "Yes");
	RadioButton (VCallRPS, "No");
	SetValue(VCallRPS, Vtrj_rps);
*/
	Vtrj_grp2 = NormalGroup(g, 2, 0, "Residue", systemFont, (GrpActnProc) NULL);
/*	StaticPrompt (Vtrj_grp2, "       ", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp2, "       ", 0, dialogTextHeight, systemFont, 'l');*/
	StaticPrompt (Vtrj_grp2, "Start #", 0, dialogTextHeight, systemFont, 'l');
	startres_text = DialogText (Vtrj_grp2, Vtrj_startres,  5, (TxtActnProc) NULL);
	StaticPrompt (Vtrj_grp2, "Finish #", 0, dialogTextHeight, systemFont, 'l');
	endres_text = DialogText (Vtrj_grp2, Vtrj_endres,  10, (TxtActnProc) NULL);
	StaticPrompt (Vtrj_grp2, "       ", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp2, "       (0: Represents Last Residue)", 0, dialogTextHeight, systemFont, 'l');
	
	
	SS_Group = NormalGroup(g, 2, 0, "Standard Deviation", systemFont, (GrpActnProc) NULL);
/*	StaticPrompt (SS_Group, "       ", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (SS_Group, "       ", 0, dialogTextHeight, systemFont, 'l');*/
	StaticPrompt (SS_Group, "For x dir.", 0, dialogTextHeight, systemFont, 'l');
	sdx_text = DialogText (SS_Group, Vtrj_sdx,  3, (TxtActnProc) NULL);
	StaticPrompt (SS_Group, "For y dir.", 0, dialogTextHeight, systemFont, 'l');
	sdy_text = DialogText (SS_Group, Vtrj_sdy,  3, (TxtActnProc) NULL);
	StaticPrompt (SS_Group, "       ", 0, dialogTextHeight, systemFont, 'l');
	

	Break(g);

	Vtrj_grp5 = NormalGroup(g, 0, 1, "Details", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp5, 50, 10);

	Vtrj_grp11 = HiddenGroup(Vtrj_grp5, 2, 0, (GrpActnProc) NULL);

	Vtrj_grp3 = HiddenGroup(Vtrj_grp11, 0, 3, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp3, 3, 10);

	StaticPrompt (Vtrj_grp3, "Gaussian Peak Height", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp3, "Save Sidechain Chi Angles", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp3, "Data Compression Type (Experts only)", 0, dialogTextHeight, systemFont, 'l');
	peakheight_text = DialogText (Vtrj_grp3, Vtrj_peakheight,  5, (TxtActnProc) NULL);
	
	savechiV_list = PopupList(Vtrj_grp3, TRUE, (PupActnProc) NULL);

	PopupItem(savechiV_list, "None");
	PopupItem(savechiV_list, "All");
	PopupItem(savechiV_list, "Buried Only");
	SetValue(savechiV_list, VMkTrj_savechi+1);
	
	comptypeV_list = PopupList(Vtrj_grp3, TRUE, (PupActnProc) NULL);

	PopupItem(comptypeV_list, "None");
	PopupItem(comptypeV_list, "RLE");
	PopupItem(comptypeV_list, "BZip");
	SetValue(comptypeV_list, VMkTrj_comprtype+1);
	
	TT_Group = HiddenGroup(Vtrj_grp5, 0, 2, (GrpActnProc) NULL);
	SetGroupSpacing(TT_Group, 3, 10);

	StaticPrompt (TT_Group, "Time Step (fs)", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (TT_Group, "Temperature (K)", 0, dialogTextHeight, systemFont, 'l');
	timestep_text = DialogText (TT_Group, Vtrj_timestep,  5, (TxtActnProc) NULL);
	temp_text = DialogText (TT_Group, Vtrj_temp,  5, (TxtActnProc) NULL);


	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp5, g, NULL);
	AlignObjects(ALIGN_RIGHT, TT_Group, Vtrj_grp5, NULL);

	Break(g);

	Vtrj_grp9 = NormalGroup(g, 2, 0, "Decide", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp9, 20, 3);

	StaticPrompt (Vtrj_grp9, "Use Standard Deviation or TimeStep/Temp.?", 0, dialogTextHeight, systemFont, 'l');
	UseSSorTT = HiddenGroup(Vtrj_grp9, 2, 0, (GrpActnProc) SetSSorTT);
	RadioButton (UseSSorTT, "Standard Deviation");
	RadioButton (UseSSorTT, "TimeStep/Temp.");
	SetValue(UseSSorTT, VMkTrj_SSorTT);

	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp9, g, NULL);

	Break(g);

	Vtrj_grp6 = NormalGroup(g, 2, 0, "Constraint File", systemFont,(GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp6, 10, 10);
	SetGroupMargins(Vtrj_grp6, 5, 5);

	NewWiz3V_text = DialogText (Vtrj_grp6, VMkTrj_constraints,  20, (TxtActnProc) NULL);
	PushButton(Vtrj_grp6, "Browse", (BtnActnProc) Wiz3VBrowseProc);

	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp6, g, NULL);

	Break(g);

	Vtrj_grp8 = NormalGroup(g, 2, 0, "Output File", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp8, 10, 0);
	SetGroupMargins(Vtrj_grp8, 5, 5);

	NewWiz3V2_text = DialogText (Vtrj_grp8, VMkTrj_trjout,  20, (TxtActnProc) NULL);
	PushButton(Vtrj_grp8, "Browse", (BtnActnProc) Wiz3V2BrowseProc);
	StaticPrompt (Vtrj_grp8, "(Default: Input File Name)", 0, dialogTextHeight, systemFont, 'l');

	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp8, g, NULL);

	Break(g);

	Vtrj_grp7 = HiddenGroup(g, 3, 0, (GrpActnProc) NULL);
	PushButton(Vtrj_grp7, "Cancel", (BtnActnProc) DoneWiz3VProc);
	PushButton(Vtrj_grp7, "<< Back", (BtnActnProc) BackWiz3VProc);
	Wiz3VNext_bttn = PushButton(Vtrj_grp7, "Next >>", (BtnActnProc) CheckWiz3VInfo);

	AlignObjects(ALIGN_CENTER, Vtrj_grp7, g, NULL);


	if(VMkTrj_SSorTT == 1)
	{
		Disable(TT_Group);
		Enable(SS_Group);
	}
	if(VMkTrj_SSorTT ==2)
	{
		Disable(SS_Group);
		Enable(TT_Group);
	}

	Show(NewWiz3V_win);
	Select(NewWiz3V2_text);
	ProcessEvents();

}

static void DisplayNewWiz3Seq(IteM i)
{
	
	Int2 Vtrj_goroverride=0, Vtrj_rps=0, Vtrj_uturn=0;

	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp4,Vtrj_grp7, Vtrj_grp9,
		Vtrj_grp10, Vtrj_grp11, Vtrj_grp12, Vtrj_grp13;

	PrompT tmp_prmt;
	Char pcdflt[PATH_MAX];

	isNew=TRUE;
	EnableDisableProc();

	if(VMkTrj_goroverride) Vtrj_goroverride=1;
	else Vtrj_goroverride=2;

	if(VMkTrj_dorps) Vtrj_rps=1;
	else Vtrj_rps=2;

	if(VMkTrj_douturn) Vtrj_uturn=1;
	else Vtrj_uturn=2;

	NewWiz3S_win = MovableModalWindow(-50,-33,-10,-10, "New Trajectory Distribution Wizard - Step 3 (Sequence File)", (WndActnProc) DoneWiz3SProc);

	g = HiddenGroup(NewWiz3S_win,0,1, (GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 15, 15);

	Vtrj_grp1 = NormalGroup(g, 1, 0, "General", systemFont, (GrpActnProc) NULL);

	Vtrj_grp9 = HiddenGroup(Vtrj_grp1, 3, 0, (GrpActnProc) NULL);
	SetGroupMargins(Vtrj_grp9, 3, 3);
	SetGroupSpacing(Vtrj_grp9, 3, 3);

	StaticPrompt (Vtrj_grp9, "Data Compresion Type", 0, dialogTextHeight, systemFont, 'l');
	comptypeS_list = PopupList(Vtrj_grp9, TRUE, (PupActnProc) NULL);
	
	Vtrj_grp12 = HiddenGroup(Vtrj_grp9, 2, 0, (GrpActnProc) NULL);
	SetGroupMargins(Vtrj_grp12, 3, 3);
	SetGroupSpacing(Vtrj_grp12, 3, 3);

	tmp_prmt = StaticPrompt (Vtrj_grp12, "        ", 0, dialogTextHeight, systemFont, 'r');
	StaticPrompt (Vtrj_grp12, "(Expert Use Only)", 0, dialogTextHeight, systemFont, 'r');
	Hide(tmp_prmt);

	Vtrj_grp10 = HiddenGroup(Vtrj_grp1, 2, 0, (GrpActnProc) NULL);
	SetGroupMargins(Vtrj_grp10, 3, 3);
	SetGroupSpacing(Vtrj_grp10, 3, 3);

	StaticPrompt (Vtrj_grp10, "Trajectory Distribution Creation Method", 0, dialogTextHeight, systemFont, 'l');
	TGtype_list = PopupList(Vtrj_grp10, TRUE, (PupActnProc) SetTGType);
	
	Vtrj_grp11 = HiddenGroup(Vtrj_grp1, 2, 0, (GrpActnProc) NULL);
	SetGroupMargins(Vtrj_grp11, 3, 3);
	SetGroupSpacing(Vtrj_grp11, 3, 3);

	StaticPrompt (Vtrj_grp11, "Call RPSBlast? (Requires CDD database, see Help for details)", 0, dialogTextHeight, systemFont, 'l');
	SCallRPS = HiddenGroup(Vtrj_grp11, 2, 0, (GrpActnProc) NULL);
	RadioButton (SCallRPS, "Yes");
	RadioButton (SCallRPS, "No");
	SetValue(SCallRPS, Vtrj_rps);

	Vtrj_grp13 = HiddenGroup(Vtrj_grp1, 2, 0, (GrpActnProc) NULL);
	SetGroupMargins(Vtrj_grp13, 3, 3);
	SetGroupSpacing(Vtrj_grp13, 3, 3);

	StaticPrompt (Vtrj_grp13, "Predict U-turns? (may be slow, requires UTURN library)", 0, dialogTextHeight, systemFont, 'l');
	SCallUturn = HiddenGroup(Vtrj_grp13, 2, 0, (GrpActnProc) NULL);
	RadioButton (SCallUturn, "Yes");
	RadioButton (SCallUturn, "No");
	SetValue(SCallUturn, Vtrj_uturn);

	PopupItem(comptypeS_list, "None");
	PopupItem(comptypeS_list, "RLE");
	PopupItem(comptypeS_list, "BZip");
	SetValue(comptypeS_list, VMkTrj_comprtype+1);

	PopupItem(TGtype_list, "Uniform");
	PopupItem(TGtype_list, "Standard");
	PopupItem(TGtype_list, "1 State Secondary Structrue");
	PopupItem(TGtype_list, "3 State Secondary Structure");
	SetValue(TGtype_list, VMkTrj_tgtype);


	Break(g);

	Vtrj_grp2 = NormalGroup(g, 2, 0, "Constraint File", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp2, 10, 10);
	SetGroupMargins(Vtrj_grp2, 5, 5);

	NewWiz3S_text = DialogText (Vtrj_grp2, VMkTrj_constraints,  20, (TxtActnProc) NULL);
	PushButton(Vtrj_grp2, "Browse", (BtnActnProc) Wiz3SBrowseProc);

	Break(g);

	Gor_Group = NormalGroup(g, 2, 0, "Secondary Structure Prediction File", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Gor_Group, 10, 10);
	SetGroupMargins(Gor_Group, 5, 5);

	NewWiz3S3_text = DialogText (Gor_Group, VMkTrj_ssout,  20, (TxtActnProc) NULL);
	PushButton(Gor_Group, "Browse", (BtnActnProc) Wiz3S3BrowseProc);
	StaticPrompt (Gor_Group, "Write Prediction to Above File or Read from the File?", 0, dialogTextHeight, systemFont, 'l');
	UseSSFile = HiddenGroup(Gor_Group, 2, 0, (GrpActnProc) NULL);
	RadioButton (UseSSFile, "Read");
	RadioButton (UseSSFile, "Write");
	SetValue(UseSSFile, Vtrj_goroverride);

	AlignObjects(ALIGN_JUSTIFY, Gor_Group, g, NULL);

	Break(g);

	Vtrj_grp4 = NormalGroup(g, 2, 0, "Output File", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp4, 10, 0);
	SetGroupMargins(Vtrj_grp4, 5, 5);

	NewWiz3S2_text = DialogText (Vtrj_grp4, VMkTrj_trjout,  20, (TxtActnProc) NULL);
	PushButton(Vtrj_grp4, "Browse", (BtnActnProc) Wiz3S2BrowseProc);
	sprintf(pcdflt,"(Default: %s)",DEFAULT_TRJ_NAME);
	StaticPrompt (Vtrj_grp4, pcdflt, 0, dialogTextHeight, systemFont, 'l');
	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp4, g, NULL);

	
	Break(g);


	Vtrj_grp7 = HiddenGroup(g, 3, 0, (GrpActnProc) NULL);
	PushButton(Vtrj_grp7, "Cancel", (BtnActnProc) DoneWiz3SProc);
	PushButton(Vtrj_grp7, "<< Back", (BtnActnProc) BackWiz3SProc);
	Wiz3SNext_bttn = PushButton(Vtrj_grp7, "Next >>", (BtnActnProc) CheckWiz3SInfo);


	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp1, g, NULL);
	AlignObjects(ALIGN_JUSTIFY, Vtrj_grp2, g, NULL);
	AlignObjects(ALIGN_CENTER, Vtrj_grp7, g, NULL);

	if(VMkTrj_tgtype==1) Disable(Gor_Group);
	if(VMkTrj_tgtype==2) Disable(Gor_Group);
	if(VMkTrj_tgtype==3) Enable(Gor_Group);
	if(VMkTrj_tgtype==4) Enable(Gor_Group);

	Show(NewWiz3S_win);
	Select(NewWiz3S2_text);
	ProcessEvents();

}

static void DisplayNewWiz2(IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3, Vtrj_grp4;
	PrompT Val_prompt, Seq_prompt, AAS_prompt;
	ButtoN Wiz2browse_bttn;
	Char Vtrj_FileIn[PATH_MAX];

	Vtrj_FileIn[0] = '\0';
	DoneWiz1Proc(i);

	isNew=TRUE;
	EnableDisableProc();

	if(VMkTrj_method ==1)
	{
		if (VMkTrj_seqtype==1) StringCpy(Vtrj_FileIn, VMkTrj_AAseq);
		if (VMkTrj_seqtype==2) StringCpy(Vtrj_FileIn, VMkTrj_seqfile);
	}
	if(VMkTrj_method ==2) StringCpy(Vtrj_FileIn, VMkTrj_valin);

	NewWiz2_win = MovableModalWindow(-50,-33,-10,-10, "New Trajectory Distribution Wizard - Step 2", (WndActnProc) DoneWiz2Proc);

	g = HiddenGroup(NewWiz2_win,0,1, (GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 10, 20);

	Vtrj_grp4 = HiddenGroup(g, 0,1, (GrpActnProc) NULL);

	StaticPrompt (Vtrj_grp4, "                                            ", 0, dialogTextHeight, systemFont, 'r');
	StaticPrompt (Vtrj_grp4, "Sequence Length = ", 0, dialogTextHeight, systemFont, 'r');
	seqlen_title = StaticPrompt (Vtrj_grp4, "    0", 0,dialogTextHeight, programFont, 'r');
	
	Break(g);

	Vtrj_grp1 = HiddenGroup(g, 0, 0, (GrpActnProc) NULL);
	
	Val_prompt = StaticPrompt (Vtrj_grp1, "NCBI-MMDB Structure File", 0, dialogTextHeight, systemFont, 'r');
	Seq_prompt = StaticPrompt (Vtrj_grp1, "Sequence File Name", 0, dialogTextHeight, systemFont, 'r');
	AAS_prompt = StaticPrompt (Vtrj_grp1, "Input Amino Acid Sequence", 0, dialogTextHeight, systemFont, 'r');
	

	Vtrj_grp2 = HiddenGroup(g, 2, 0, (GrpActnProc) NULL);

	NewWiz2_text = DialogText (Vtrj_grp2, Vtrj_FileIn,  20, (TxtActnProc) EnableWiz2Next);
	Wiz2browse_bttn = PushButton(Vtrj_grp2, "Insert A.A.", (BtnActnProc) Wiz2BrowseProc);

	Break(g);

	Vtrj_grp3 = HiddenGroup(g, 3, 0, (GrpActnProc) NULL);
	
	PushButton(Vtrj_grp3, "Cancel", (BtnActnProc) DoneWiz2Proc);
	PushButton(Vtrj_grp3, "<< Back", (BtnActnProc) BackWiz2Proc);
	Wiz2Next_bttn = PushButton(Vtrj_grp3, "Next >>", (BtnActnProc) CheckWiz2Info);

	AlignObjects(ALIGN_CENTER, Vtrj_grp3, g, NULL);

	Disable(Wiz2Next_bttn);

	if (VMkTrj_method ==2)
	{
		SetTitle(Wiz2browse_bttn, "Browse");
		Show(Val_prompt);
		Hide(Seq_prompt);
		Hide(AAS_prompt);
		Hide(Vtrj_grp4);
		if (VMkTrj_valin[0] != '\0') Enable(Wiz2Next_bttn);
	}
	if (VMkTrj_method ==1)
	{
		if (VMkTrj_seqtype==2)
		{
			SetTitle(Wiz2browse_bttn, "Browse");
			Show(Seq_prompt);
			Hide(Val_prompt);
			Hide(AAS_prompt);
			Hide(Vtrj_grp4);
			Enable(Wiz2Next_bttn);
		}
		if (VMkTrj_seqtype==1)
		{
			SetTitle(Wiz2browse_bttn, "Insert A.A.");
			Show(AAS_prompt);
			Show(Vtrj_grp4);
			Hide(Val_prompt);
			Hide(Seq_prompt);
			if (VMkTrj_AAseq[0] != '\0') Enable(Wiz2Next_bttn);
		}
	}


	Show(NewWiz2_win);
	Select(NewWiz2_text);
	ProcessEvents();
}


static void DisplayNewWiz1(void)
{
	GrouP g, Vtrj_grp1;
		
	isNew=TRUE;
	EnableDisableProc();

	NewWiz1_win = MovableModalWindow(-50,-33,-10,-10, "New Trajectory Distribution Wizard - Step 1", (WndActnProc) DoneWiz1Proc);

	g = HiddenGroup(NewWiz1_win,0,1, (GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 50, 20);			

	
	Wiz1InputType = NormalGroup (g, 2, 0, "Input Type", systemFont, (GrpActnProc) CheckInputTypeProc);
	RadioButton (Wiz1InputType, "Sequence");
	RadioButton (Wiz1InputType, "MMDB Structure");
	SetValue(Wiz1InputType, VMkTrj_method);
	
	SequenceType = NormalGroup(g,2,0, "Sequence Type", systemFont, (GrpActnProc) CheckSequenceTypeProc);
	RadioButton (SequenceType, "Amino Acid Input");
	RadioButton (SequenceType, "From File");
	SetValue(SequenceType, VMkTrj_seqtype);
	Hide(SequenceType);


	Break(g);

	Vtrj_grp1 = HiddenGroup (g, 2, 0, (GrpActnProc) NULL);
	PushButton(Vtrj_grp1, "Cancel", (BtnActnProc) DoneWiz1Proc);
	DefaultButton(Vtrj_grp1, "Next >>", (BtnActnProc) DisplayNewWiz2);

	AlignObjects(ALIGN_CENTER, Vtrj_grp1, g, NULL);
	
	if(VMkTrj_method==1) Show(SequenceType);
	
	Show(NewWiz1_win);
	ProcessEvents();
}

/*------------------------------------------------------------*/

static void CheckSequenceTypeProc(IteM i)
{
	ClearSorVVariables();

	VMkTrj_seqtype = GetValue(SequenceType);
}

static void CheckInputTypeProc(IteM i)
{

	ClearSorVVariables();

	VMkTrj_method = GetValue(Wiz1InputType);

	if (VMkTrj_method ==1) Show(SequenceType);
	if (VMkTrj_method ==2) Hide(SequenceType);

}

/*------------------------------------------------------------*/

static void CreateNewTrajProc(IteM i)
{
	Int2 Vtrj_ans = ANS_NO;
	Boolean is_OK = TRUE;


	isNew=TRUE;

	if (isLoaded==TRUE)
	{
		Vtrj_ans = Message (MSG_YN, "The current trajectory distribution must be closed. Would you like to close it now?");	
	
		if (Vtrj_ans==ANS_NO) return;

		if (Vtrj_ans==ANS_YES) 
		{
		
			is_OK = CheckIfSaved();
			if (is_OK==FALSE) return;
			CloseTrajProc(i);
		}
	}

	Hide(aaname_panel);
	Show(aaname_panel);
	
	EnableDisableProc();
	DefaultVisMakeTrajVariables();
	DisplayNewWiz1();

}

/***************************************************************
										FOLDTRAJ FUNCTIONS
****************************************************************/

TrajErr VisFoldTraj(Int2 dumptype)
{
	static Boolean rotlibloaded=FALSE;	/* only load rotamer library once */
	FILE *fp;
	AsnIoPtr aip;
	Char fnamout[PATH_MAX];
	NcbiMimeAsn1Ptr nmap;
	BiostrucPtr bspBiostruc,bspTemp,pbs;
	PDNMS pdnmsModelstruc=NULL;
	PDNML pdnmlModel;
	PMSD pmsdHead;
	PMLD pmldThis;
	TrajErr err;
	Boolean dumppdb,dumpval;
	Int4 dbsize;
  pFoldTrajParamBlock foldtrajprm;
#ifdef WIN_MOTIF
	int childpid;
#endif	
	
	if (dumptype==DUMPTYPE_VAL) {
		dumpval=TRUE;
		dumppdb=FALSE;
	}
	else if (dumptype==DUMPTYPE_PDB) {
		dumpval=FALSE;
		dumppdb=TRUE;
	}
	else {
		Message (MSG_OK,"VisFoldTraj: Invalid Dump Type");
		/* invalid call */
		return ERR_FAIL;
	}
	if (HasBeenFiltered==TRUE)
	{		
		Vtrj_SaveResidue();
		HasBeenFiltered = FALSE;
	}
	if (!rotlibloaded) {
		Vtrj_mon = Nlm_MonitorStrNewEx("FoldTraj", 30, FALSE);
		MonitorStrValue(Vtrj_mon, "Loading Rotamer Library...");
		if (LoadRotLib()!=ERR_SUCCESS) {
			MonitorFree(Vtrj_mon);
			Message (MSG_OK, "Cannot open rotamer library");
			return ERR_FAIL;
		}
		MonitorFree(Vtrj_mon);
		rotlibloaded=TRUE;
	}
	TmpNam(tmpskelfname);
	BuildSkelASN(sequence,tmpskelfname);
	bspTempLog=BSNew(0);
	bspBiostruc=NULL;
	bspBiostruc=MIMEBiostrucAsnGet(tmpskelfname,"r",NULL);
  if (bspBiostruc==NULL) {
		FileRemove(tmpskelfname);
    Message(MSG_OK,"Unable to fetch Biostruc");
	  return ERR_FAIL;
  }
  pdnmsModelstruc=MakeAModelstruc(bspBiostruc);
  if (pdnmsModelstruc==NULL) {
    Message(MSG_OK,"Unable to convert Biostruc to Modelstruc.");
    return ERR_FAIL;
  }
	pmsdHead=(PMSD)(pdnmsModelstruc->data.ptrvalue);
	/* remove ppAsnOrder for all models so it will be rebuilt
	   in WriteAsnAllModel */
	pdnmlModel=pmsdHead->pdnmlModels;
	while (pdnmlModel) {
	 	pmldThis=(PMLD)(pdnmlModel->data.ptrvalue);
	  if (pmldThis->ppAsnOrder) {
	    PTRVectorFree(pmldThis->ppAsnOrder,0);
	    pmldThis->ppAsnOrder=NULL;
	  }
	  pdnmlModel=pdnmlModel->next;
	}
	/* remove id and descr if necessary */
	/* assign new id and descr */
	pbs=pmsdHead->pbsBS;
	if (pbi!=NULL) {
		if (pbs->id!=NULL)
			AsnGenericChoiceSeqOfFree(pbs->id,(AsnOptFreeFunc)BiostrucIdFree);
		pbs->id=pbi;
	}
	if (pbd!=NULL) {
		if (pbs->descr!=NULL)
			AsnGenericChoiceSeqOfFree(pbs->descr,(AsnOptFreeFunc)BiostrucDescrFree);
		pbs->descr=pbd;
	}
	/* make a progress monitor for foldtraj now */
	TGInit(tmpdbasename,DB_READ,&dbsize);
	if (dbsize!=((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue))->iResCount) {
		Message(MSG_OK,"protein length inconsistency error, expect: %d actual: %d, possibly due to corrupt trajectory file; aborting",((PMMD)((pmsdHead->pdnmmHead)->data.ptrvalue))->iResCount,dbsize);
		TGClose();	
		return ERR_FAIL;		
	}	
	foldtrajprm=(pFoldTrajParamBlock)MemNew(sizeof(FoldTrajParamBlock));
	foldtrajprm->pmsdRoot=pmsdHead;
	foldtrajprm->Model=(Int2)1;
	foldtrajprm->err=0;
	foldtrajprm->gen=0;
	foldtrajprm->errorfile=NULL;
	DisableAllItems();
	if (NlmThreadsAvailable()==TRUE) {
		NlmThreadCreate(TRADEProtein, (VoidPtr)foldtrajprm);
		do {
		} while (ProgramProgressMax<=0);
		Vtrj_mon = MonitorIntNew("Folding Protein...", 0, (Nlm_Int4)ProgramProgressMax);
		MonitorFoldTraj();
		NlmThreadJoinAll();
		MonitorFree(Vtrj_mon);
	}
	else {
		Vtrj_mon = Nlm_MonitorStrNewEx("FoldTraj", 20, FALSE);
		MonitorStrValue(Vtrj_mon, "Folding Protein...");
		TRADEProtein((VoidPtr)foldtrajprm);
		MonitorFree(Vtrj_mon);
	}
	err=foldtrajprm->err;
	foldtrajprm=MemFree(foldtrajprm);
	TGClose();
	ProgramProgress=0;
	ProgramProgressMax=0;

	EnableAllItems();
	
 	StringCpy(fnamout,GENERATE_TMP_FILENAME);
	if (dumpval && err!=ERR_CANCELLED) {
		StringCat(fnamout,MMDB_EXT);				
		/* write out as MIME Biostruc type */
    if (!(WriteAsnAllModel(pdnmsModelstruc,fnamout,SAVE_BINARY))) {
			Message(MSG_OK,"ASNWrite failed");
			return ERR_FAIL;
		}
		bspTemp=NULL;
		bspTemp=MIMEBiostrucAsnGet(fnamout,"rb",NULL);
		FileRemove(fnamout);
		aip=AsnIoOpen(fnamout,"wb");
		nmap=BuildMIMEBiostruc(bspTemp,Xsequence,vnpBioseq);
		NcbiMimeAsn1AsnWrite(nmap,aip,NULL);
		AsnIoClose(aip);
		/* remove Bioseq from nmap since free-er frees it */
		((BiostrucSeqPtr)(nmap->data.ptrvalue))->sequences=NULL;
		NcbiMimeAsn1Free(nmap);
	}
	if (dumppdb && err!=ERR_CANCELLED) {
		StringCat(fnamout,PDB_EXT);
		fp=FileOpen(fnamout,"w");
		if (!(WritePDBAllModel(pdnmsModelstruc,fp))) {
			FileClose(fp);
			Message(MSG_OK,"PDBWrite failed");
			return ERR_FAIL;
		}
		FileClose(fp);
	}
	if (err==ERR_INCOMPLETE)
 		Message(MSG_OK,"Warning: Timed out, your structure file is incomplete.");
	/* trick freer so pbi and pbd are saved for next run */
	pbs->id=NULL;
	pbs->descr=NULL;
	FreeAModelstruc(pdnmsModelstruc);
	/* Free useless log now */	
	bspTempLog=BSFree(bspTempLog);
	FileRemove(tmpskelfname);
	if (err==ERR_CANCELLED)
		return err;
	/* launch viewer, file = fnamout */
#ifdef WIN_MOTIF	
	childpid=fork();
	if (childpid==0) {
		/* child process */
		if (dumpval) {
			if (execlp("Cn3D","Cn3D",fnamout,NULL)<0) {
				Message (MSG_OK,"Unable to load Cn3D");
				exit(1);
			}
		}
		if (dumppdb) {
			if (execlp("rasmol","rasmol",fnamout,NULL)<0) {
				Message (MSG_OK,"Unable to load Rasmol");
				exit(1);
			}
		}
	}
#endif
#ifdef WIN_MSWIN
	if ((int)ShellExecute(0,"open",fnamout,NULL,NULL,SW_SHOWNORMAL)<=32) {
		if (dumpval)
			Message (MSG_OK,"Unable to load Cn3D");
		if (dumppdb)
			Message (MSG_OK,"Unable to load PDB viewer");
	}
#endif
	/* parent continues */
	return ERR_SUCCESS;
}		
	
static void GeneratePdbProc(IteM i)
{
	if (VisFoldTraj(DUMPTYPE_PDB)!=ERR_SUCCESS) {
		/* an error occurred during foldtraj */
	}
}

static void GenerateValProc(IteM i)
{
	if (VisFoldTraj(DUMPTYPE_VAL)!=ERR_SUCCESS) {
		/* an error occurred during foldtraj */
	}
}

/***************************************************************
										MULTIPLE VIEW FUNCTIONS
****************************************************************/
/*static void DoneMultView(IteM i)
{
	Hide(MultView_win);
}

static void EnableMultOK(void)
{
	Enable(MultOK_bttn);	
}*/

/*static void GetMultNum(IteM i)
{
	Int2 numofres;
	
	numofres = GetValue(multnum_list);
	
	DoneMultView(i);

	Message (MSG_OK, "Feature to come in version 2");
	
}*/


/*static void MultipleViewProc (IteM i)
{

	GrouP g, Vtrj_grp1, Vtrj_grp2;
	
	MultView_win = MovableModalWindow(-50,-33,-10,-10, "Multiple View",  (WndActnProc) DoneMultView);
	
	g = HiddenGroup(MultView_win, -1, -1, (GrpActnProc) NULL);

	SetGroupSpacing(g, 3, 15);
	SetGroupMargins(g, 3, 3);

	Vtrj_grp1 = HiddenGroup(g, 0, 2, (GrpActnProc) NULL);

	StaticPrompt (Vtrj_grp1, "Display how many residues?", 0, dialogTextHeight, systemFont, 'l');
        	
    multnum_list = PopupList(Vtrj_grp1, TRUE, (PupActnProc) EnableMultOK);

	PopupItem(multnum_list, "One");
	PopupItem(multnum_list, "Two");
	PopupItem(multnum_list, "Three");
	PopupItem(multnum_list, "Four");
	PopupItem(multnum_list, "Five");
	PopupItem(multnum_list, "Six");

	Vtrj_grp2 = HiddenGroup(g, 2, -1, (GrpActnProc) NULL);
	
	MultOK_bttn = PushButton(Vtrj_grp2, "OK", (BtnActnProc) GetMultNum);
	Disable(MultOK_bttn);

	DefaultButton(Vtrj_grp2, "Cancel", (BtnActnProc) DoneMultView);
               	
	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

    Show(MultView_win);
	
	ProcessEvents();

}*/

/*****************************************************************
										DISTANCE CONSTRAINTS FUNCTIONS
******************************************************************/

static void DoneAddDistProc(IteM i)
{
	Hide(adddist_win);
	constrwinopen=FALSE;
	DisableAllItems();
}

static void DoneDistProc(IteM i)
{
	if (constrwinopen==TRUE)
		return;
	Hide(editdist_win);
	EnableAllItems();
}

static void ApplyDistProc(IteM i)
{
	PNN pnnHere=NULL;
	Char tmp_char[20], atomname1[5], atomname2[5];
	int itemp=0,err;
	float ftemp=0.0;
	Int2 ires1=0, ires2=0, whichdist=0, x=0;
	Int4 residx1=0, residx2=0;
	FloatLo fmeandist=0.0, fmindelta=0.0, fmaxdelta=0.0, fang1=0.0, fang2=0.0,
		fdihed01=0.0, fdihed12=0.0, fdihed23=0.0, fprob=0.0;


	tmp_char[0] = '\0';

	
	GetTitle(resnum1_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for Residue #1.");
		Select(resnum1_text);
		return;
	}
	err=sscanf(tmp_char, "%d", &itemp);
	if (itemp < 1 || itemp > numAA || err==0)
	{
		Message (MSG_OK, "Please enter a valid residue number for Residue #1.");
		Select(resnum1_text);
		return;
	}

	ires1= (Int2) itemp;
	
	GetTitle(resnum2_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for Residue #2.");
		Select(resnum2_text);
		return;
	}
	err=sscanf(tmp_char, "%d", &itemp);
	if (itemp < 1 || itemp > numAA || err==0)
	{
		Message (MSG_OK, "Please enter a valid residue number for Residue #2.");
		Select(resnum2_text);
		return;
	}
	if (itemp == (int)ires1)
	{
		Message (MSG_OK, "Residue #2 must me different from Residue #1.");
		Select(resnum2_text);
		return;
	}
	ires2= (Int2) itemp;

	GetTitle(resname1_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a name for Atom #1.");
		Select(resname1_text);
		return;
	}
	if (StringLen(tmp_char)>4)
	{
		Message (MSG_OK, "The name for Atom #1 cannot be longer than 4 characters.");
		Select(resname1_text);
		return;
	}
	StringCpy(atomname1, tmp_char);
	StrUpper(atomname1);
	atomname1[4] = '\0';

	GetTitle(resname2_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a name for Atom #2.");
		Select(resname2_text);
		return;
	}
	if (StringLen(tmp_char)>4)
	{
		Message (MSG_OK, "The name for Atom #2 cannot be longer than 4 characters.");
		Select(resname2_text);
		return;
	}
	StringCpy(atomname2, tmp_char);
	StrUpper(atomname2);
	atomname2[4] = '\0';

	GetTitle(meandist_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for the Mean Distance");
		Select(meandist_text);
		return;
	}
	err=sscanf(tmp_char, "%f", &ftemp);
	if (ftemp<=0.0 || err==0)
	{
		Message (MSG_OK, "Mean Distance must be positive.");
		Select(meandist_text);
		return;
	}
	fmeandist= (FloatLo) ftemp;

	GetTitle(mindelta_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for the Minimum Delta");
		Select(mindelta_text);
		return;
	}
	err=sscanf(tmp_char, "%f", &ftemp);
	if (ftemp<0.0 || ftemp>=(float)fmeandist || err==0)
	{
		Message (MSG_OK, "Minimum Delta must be non-negative and less than Mean Distance");
		Select(mindelta_text);
		return;
	}
	fmindelta= (FloatLo) ftemp;

	GetTitle(maxdelta_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for the Maximum Delta");
		Select(maxdelta_text);
		return;
	}
	err=sscanf(tmp_char, "%f", &ftemp);
	if (ftemp<0.0 || err==0)
	{
		Message (MSG_OK, "Maximum Delta must be non-negative");
		Select(maxdelta_text);
		return;
	}
	fmaxdelta= (FloatLo) ftemp;

	GetTitle(angle1_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
	    fang1=CONSTR_INFINITY;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (((ftemp<=0.0 || ftemp>180.0) && (ftemp!=CONSTR_INFINITY)) || err==0)
		{
			Message (MSG_OK, "Angle #1 must be greater than 0 and less than or equal to 180 degrees");
			Select(angle1_text);
			return;
		}
		fang1= (FloatLo) ftemp;
	}

	GetTitle(angle2_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
	    fang2=CONSTR_INFINITY;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (((ftemp<=0.0 || ftemp>180.0) && (ftemp!=CONSTR_INFINITY)) || err==0)
		{
			Message (MSG_OK, "Angle #2 must be greater than 0 and less than or equal to 180 degrees");
			Select(angle2_text);
			return;
		}
		fang2= (FloatLo) ftemp;
	}

	GetTitle(dihedral01_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
	    fdihed01=CONSTR_INFINITY;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (((ftemp<=0.0 || ftemp>360.0) && (ftemp!=CONSTR_INFINITY)) || err==0)
		{
			Message (MSG_OK, "Dihedral 0-1 must be greater than 0 and less than or equal to 360 degrees");
			Select(dihedral01_text);
			return;
		}
		fdihed01= (FloatLo) ftemp;
	}

	GetTitle(dihedral12_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
	    fdihed12=CONSTR_INFINITY;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (((ftemp<=0.0 || ftemp>360.0) && (ftemp!=CONSTR_INFINITY)) || err==0)
		{
			Message (MSG_OK, "Dihedral 1-2 must be greater than 0 and less than or equal to 360 degrees");
			Select(dihedral12_text);
			return;
		}
		fdihed12= (FloatLo) ftemp;
	}

	GetTitle(dihedral23_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
	    fdihed23=CONSTR_INFINITY;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (((ftemp<=0.0 || ftemp>360.0) && (ftemp!=CONSTR_INFINITY)) || err==0)
		{
			Message (MSG_OK, "Dihedral 2-3 must be greater than 0 and less than or equal to 360 degrees");
			Select(dihedral23_text);
			return;
		}
		fdihed23= (FloatLo) ftemp;
	}

	GetTitle(prob_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
	    fprob=0.0;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<0.0 || ftemp>1.0 || err==0)
		{
			Message (MSG_OK, "Probability must be between zero and one");
			Select(prob_text);
			return;
		}
		fprob= (FloatLo) ftemp;
	}

	residx1=GetResIdxFromSeq(sequence,ires1);
	residx2=GetResIdxFromSeq(sequence,ires2);
	if(SubstituteNames(atomname1,residx1)!=ERR_SUCCESS)
	{
		Message (MSG_OK, "The first atom name entered is invalid.");
		Select(resname1_text);
		return;
	}
		
	if(SubstituteNames(atomname2,residx2) != ERR_SUCCESS)
	{
		Message (MSG_OK, "The second atom name entered is invalid.");
		Select(resname2_text);
		return;
	}
	
	if(editdist==FALSE) pnnHere=(PNN)MemNew(sizeof(NN));
	
	if(editdist==TRUE)
	{
		whichdist = GetValue(distlist);
		pnnHere = Vtrj_pnn;
		for(x=1; x<whichdist; x++)	pnnHere = pnnHere->next;
	}
	

	if (ires1<ires2)
	{
		pnnHere->res1=(Int2)ires1;
		pnnHere->res2=(Int2)ires2;
		StringCpy(pnnHere->AtomName1,atomname1);
		StringCpy(pnnHere->AtomName2,atomname2);
	}
	else
	{
		pnnHere->res1=(Int2)ires2;
		pnnHere->res2=(Int2)ires1;
		StringCpy(pnnHere->AtomName2,atomname2);
		StringCpy(pnnHere->AtomName1,atomname1);
	}
	pnnHere->MeanDist=(FloatLo)fmeandist;
	pnnHere->MinDelta=(FloatLo)fmindelta;
	pnnHere->MaxDelta=(FloatLo)fmaxdelta;
	pnnHere->Angle1=(FloatLo)fang1;
	pnnHere->Angle2=(FloatLo)fang2;
	pnnHere->Dihedral01=(FloatLo)fdihed01;
	pnnHere->Dihedral12=(FloatLo)fdihed12;
	pnnHere->Dihedral23=(FloatLo)fdihed23;
	if (fprob==0.0)
		pnnHere->prob= NOE_PROB;
	else
		pnnHere->prob=(FloatLo)fprob;
	pnnHere->tries=0;
	
	if (editdist==FALSE)
	{
		if(AddDistConstraint(pnnHere) != ERR_SUCCESS)
			Message (MSG_OK, "The new distance constraint was not correctly added. \n Please refer to the VisTraj Help File for further information.");
	}

	DistBeenEdited=TRUE;

	Vtrj_pnn = GetTrueDistConstraints(); 

	if (Vtrj_pnn->next==NULL)
		SetValue(distlist,1);
	Enable(dist_doc);
	Enable(distlist);

	LoadDistProc(i);
	LoadDistList();
	DoneAddDistProc(i);
}


static void DoDistProc(IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3, Vtrj_grp4, Vtrj_grp5, Vtrj_grp6;
	PNN tmp_pnn=NULL;
	Int2 whichdist=0, x=0;
	Char temp_str[10];
	Char tmpprob[PATH_MAX];	

	adddist_win = MovableModalWindow(-50,-33,-10,-10, "", (WndActnProc) DoneAddDistProc);

	Reset(adddist_win);

	if(editdist==TRUE)
	{
		SetTitle(adddist_win, "Edit Distance Constraints");
		temp_str[0] = '\0';
		tmp_pnn = Vtrj_pnn;
		whichdist = GetValue(distlist);
		for(x=1; x<whichdist; x++)	tmp_pnn = tmp_pnn->next;
	}

	if(editdist==FALSE) SetTitle(adddist_win, "Add Distance Constriants");

	g = HiddenGroup(adddist_win,-1,-1, (GrpActnProc) NULL);
	
	Vtrj_grp1 =  HiddenGroup(g,1,0,(GrpActnProc) NULL);
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 5, 5);			
	
	Vtrj_grp3 = NormalGroup(Vtrj_grp1, 4, 0, "Constrained Residues", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp3, 10, 10);			


	StaticPrompt (Vtrj_grp3, "1st Residue #", 0, dialogTextHeight, systemFont, 'l');    	
    resnum1_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) NULL);
	if(editdist==TRUE) Disable(resnum1_text);
	if(editdist==FALSE) Enable(resnum1_text);
	
	StaticPrompt (Vtrj_grp3, "1st Atom Name", 0, dialogTextHeight, systemFont, 'l');    	
    resname1_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) NULL);
	
	StaticPrompt (Vtrj_grp3, "2nd Residue #", 0, dialogTextHeight, systemFont, 'l');    	
    resnum2_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) NULL);
	if(editdist==TRUE) Disable(resnum2_text);
	if(editdist==FALSE) Enable(resnum2_text);

	StaticPrompt (Vtrj_grp3, "2nd Atom Name", 0, dialogTextHeight, systemFont, 'l');    	
    resname2_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) NULL);




	Vtrj_grp4 = NormalGroup(Vtrj_grp1, 1, 0, "Probabilistic Data", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp4, 3, 3);			
	
	Vtrj_grp6 = HiddenGroup(Vtrj_grp4, 4, 0, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp6, 10, 10);			

	StaticPrompt (Vtrj_grp6, "Mean Distance", 0, dialogTextHeight, systemFont, 'l');    	
    meandist_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp6, "Minimum Delta", 0, dialogTextHeight, systemFont, 'l');    	
    mindelta_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp6, "Maximum Delta", 0, dialogTextHeight, systemFont, 'l');    	
    maxdelta_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp6, "Probability", 0, dialogTextHeight, systemFont, 'l');    	
    prob_text = DialogText (Vtrj_grp6, "0.0",  5, (TxtActnProc) NULL);

	sprintf(tmpprob,"                                                       Put 0.0 to use default probability of %4.2f",NOE_PROB);
	StaticPrompt (Vtrj_grp4, tmpprob, 0, dialogTextHeight, systemFont, 'l');    	

	sprintf(tmpprob,"Angles and Dihedrals (optional) (%1.2f or blank = ignored field)",CONSTR_INFINITY);
	Vtrj_grp5 = NormalGroup(Vtrj_grp1, 4, 0, tmpprob, systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp5, 10, 10);			

	StaticPrompt (Vtrj_grp5, "Dihedral 0-1", 0, dialogTextHeight, systemFont, 'l');    	
    dihedral01_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp5, "Angle #1", 0, dialogTextHeight, systemFont, 'l');    	
    angle1_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp5, "Dihedral 1-2", 0, dialogTextHeight, systemFont, 'l');    	
    dihedral12_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp5, "Angle #2", 0, dialogTextHeight, systemFont, 'l');    	
    angle2_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp5, "Dihedral 2-3", 0, dialogTextHeight, systemFont, 'l');    	
    dihedral23_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);


	Vtrj_grp2 = HiddenGroup(g,2,-1,(GrpActnProc) NULL);
	
	if(editdist==FALSE)
		PushButton(Vtrj_grp2, "Add", (BtnActnProc) ApplyDistProc);
	else
		PushButton(Vtrj_grp2, "Apply", (BtnActnProc) ApplyDistProc);
	PushButton(Vtrj_grp2, "Cancel", (BtnActnProc) DoneAddDistProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);


	if(editdist==TRUE)
	{
		sprintf(temp_str, "%d", tmp_pnn->res1);
		SetTitle(resnum1_text, temp_str);
		
		sprintf(temp_str, "%d", tmp_pnn->res2);
		SetTitle(resnum2_text, temp_str);
		
		SetTitle(resname1_text, tmp_pnn->AtomName1);
		
		SetTitle(resname2_text, tmp_pnn->AtomName2);
		
		sprintf(temp_str, "%6.2f", tmp_pnn->Angle1);
		SetTitle(angle1_text, temp_str);
		
		sprintf(temp_str, "%6.2f", tmp_pnn->Angle2);
		SetTitle(angle2_text, temp_str);
		
		sprintf(temp_str, "%6.2f", tmp_pnn->Dihedral01);
		SetTitle(dihedral01_text, temp_str);
		
		sprintf(temp_str, "%6.2f", tmp_pnn->Dihedral12);
		SetTitle(dihedral12_text, temp_str);
		
		sprintf(temp_str, "%6.2f", tmp_pnn->Dihedral23);
		SetTitle(dihedral23_text, temp_str);
		
		sprintf(temp_str, "%5.3f", tmp_pnn->prob);
		SetTitle(prob_text, temp_str);

		sprintf(temp_str, "%5.2f", tmp_pnn->MeanDist);
		SetTitle(meandist_text, temp_str);

		sprintf(temp_str, "%5.2f", tmp_pnn->MaxDelta);
		SetTitle(maxdelta_text, temp_str);

		sprintf(temp_str, "%5.2f", tmp_pnn->MinDelta);
		SetTitle(mindelta_text, temp_str);

	}

	constrwinopen=TRUE;
	Show(adddist_win);
	ProcessEvents();	

}


static void AddDistProc(IteM i)
{
	editdist=FALSE;
	DoDistProc(i);
}

static void EditDistProc(IteM i)
{
	editdist=TRUE;
	DoDistProc(i);
}


static void DeleteDistProc(IteM i)
{

	PNN tmp_pnn=NULL;
	Int2 whichdist=0, x=0;
	Int2 Vtrj_ans=ANS_NO;
	Boolean listend=FALSE;

	
	Vtrj_ans = Message (MSG_YN, "Are you sure you would like to delete this distance constraint?");	

	if(Vtrj_ans==ANS_NO) return;

	tmp_pnn = Vtrj_pnn;
	whichdist = GetValue(distlist);
	for(x=1; x<whichdist; x++)	tmp_pnn = tmp_pnn->next;
	if (tmp_pnn->next==NULL)
		listend=TRUE;
	DeleteDistConstraint(tmp_pnn);

	Vtrj_pnn = GetTrueDistConstraints();
	
	DistBeenEdited=TRUE;

	Reset(dist_doc);
	if (listend==TRUE)
		whichdist--;
	if (whichdist>0)
		SetValue(distlist,whichdist);
	LoadDistList();
	if (whichdist>0)
		LoadDistProc(i);

}

static void LoadDistProc(IteM i)
{
	PNN tmp_pnn=NULL;
	Int2 whichdist=0, x=0;
	Char tmp_str[10];
	

	tmp_pnn = Vtrj_pnn;
	if (tmp_pnn==NULL)
		return;
	tmp_str[0] = '\0';

	Enable(distE_bttn);
	Enable(distD_bttn);

	Reset(dist_doc);


	whichdist = GetValue(distlist);

	for(x=1; x<whichdist; x++)	tmp_pnn = tmp_pnn->next;

	sprintf(tmp_str, "%d", tmp_pnn->res1);
	AppendText (dist_doc, "Residue #1:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%d", tmp_pnn->res2);
	AppendText (dist_doc, "Residue #2:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);
			
	StringCpy(tmp_str, tmp_pnn->AtomName1);
	AppendText (dist_doc, "Atom Name #1:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);
		
	StringCpy(tmp_str, tmp_pnn->AtomName2);
	AppendText (dist_doc, "Atom Name #2:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%5.2f", tmp_pnn->MeanDist);
	AppendText (dist_doc, "Mean Distance:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%5.2f", tmp_pnn->MinDelta);
	AppendText (dist_doc, "Minimum Delta:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%5.2f", tmp_pnn->MaxDelta);
	AppendText (dist_doc, "Maximum Delta:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%6.2f", tmp_pnn->Angle1);
	AppendText (dist_doc, "Angle #1:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%6.2f", tmp_pnn->Angle2);
	AppendText (dist_doc, "Angle #2:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%6.2f", tmp_pnn->Dihedral01);
	AppendText (dist_doc, "Dihedral 0-1:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%6.2f", tmp_pnn->Dihedral12);
	AppendText (dist_doc, "Dihedral 1-2:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%6.2f", tmp_pnn->Dihedral23);
	AppendText (dist_doc, "Dihedral 2-3:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%5.3f", tmp_pnn->prob);
	AppendText (dist_doc, "Probability:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (dist_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);


	UpdateDocument (dist_doc, 0, 0);
	
}

static void LoadDistList(void)
{

	PNN tmp_pnn=NULL;
	Char AtN1[5], AtN2[5], res1_str[10], res2_str[10], dist_title[40];
	Int2 res1=0, res2=0, whichdist=0, numdist=0;

/*	Disable(distE_bttn);
	Disable(distD_bttn);*/

	DishdgColFmt.font = ParseFont ("Times,14,b");
	DissubColFmt.font = ParseFont ("Helvetica,11,b");
	DistxtColFmt.font = ParseFont ("fixed,10");
/*	lstColFmt.font = programFont;
	tblColFmt.font = programFont;*/

	whichdist=GetValue(distlist);
	if (whichdist==0)
		whichdist=1;
	Reset(distlist);

	tmp_pnn = Vtrj_pnn;

	if(tmp_pnn==NULL)
	{
		Reset(dist_doc);
		AppendText (dist_doc, "There are no distance constraints.", &hdgParFmt, &DishdgColFmt, programFont);
		UpdateDocument (dist_doc, 0, 0);
		Disable(distlist);
		Disable(dist_doc);
		Disable(distE_bttn);
		Disable(distD_bttn);
	}

	while(tmp_pnn !=NULL)
	{
		dist_title[0] = '\0';
		res1 = tmp_pnn->res1;
		res2 = tmp_pnn->res2;
		StringCpy(AtN1, tmp_pnn->AtomName1);
		StringCpy(AtN2, tmp_pnn->AtomName2);

		sprintf(res1_str, "%d", res1);
		sprintf(res2_str, "%d", res2);
		
		StringCat(dist_title, res1_str);
		StringCat(dist_title, " - ");
		StringCat(dist_title, res2_str);
		StringCat(dist_title, " (");
		StringCat(dist_title, AtN1);
		StringCat(dist_title, " - ");
		StringCat(dist_title, AtN2);
		StringCat(dist_title, ")");

		ListItem(distlist, dist_title);
		numdist++;

		tmp_pnn = tmp_pnn->next;
	}
	
	if (whichdist>numdist)
		whichdist=numdist;
	SetValue(distlist,whichdist);

}

static void DisplayDistConsProc (IteM i)
{

	GrouP g, Vtrj_grp1, Vtrj_grp2;

	editdist_win = MovableModalWindow(-50,-33,-10,-10, "Distance Constraints", (WndActnProc) DoneDistProc);

	g = HiddenGroup(editdist_win,-1,-1, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 20, 20);			
	
	distlist = SingleList (Vtrj_grp1,10,15, (LstActnProc) LoadDistProc);

	dist_doc = DocumentPanel (Vtrj_grp1, stdCharWidth * 20, stdLineHeight * 15);
	
	Vtrj_grp2 = HiddenGroup(g,4,-1,(GrpActnProc) NULL);
	
	SetGroupSpacing(g, 15, 15);

	distE_bttn = PushButton(Vtrj_grp2, "Edit", (BtnActnProc) EditDistProc);
	distA_bttn = PushButton(Vtrj_grp2, "Add", (BtnActnProc) AddDistProc);
	distD_bttn = PushButton(Vtrj_grp2, "Delete", (BtnActnProc) DeleteDistProc);
	PushButton(Vtrj_grp2, "Close", (BtnActnProc) DoneDistProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

	LoadDistList();
	LoadDistProc(i);

	Show(editdist_win);
	ProcessEvents();	

}

/*****************************************************************
										FRAGMENT LIST FUNCTIONS
******************************************************************/

static void DoneAddFragProc(IteM i)
{
	Hide(addfrag_win);
	fragwinopenadd=FALSE;
	if (fragwinopenres==FALSE)
		Enable(fraglist);
	DisableAllItems();
}

static void DoneFragListProc(IteM i)
{
	if (fragwinopenres==TRUE || fragwinopenadd==TRUE)
		return;
	Hide(editfraglist_win);
	EnableAllItems();
}

static void DoneFragProc(IteM i)
{
	if (fragwinopenadd==TRUE)
		return;
	Hide(editfrag_win);
	fragwinopenres=FALSE;
	Enable(fraglist);
	DisableAllItems();
}

static void ApplyFragProc(IteM i)
{
	PFDS pfdsHere,pfdstmp;
	Char tmp_char[20];
	int itemp=0,err;
	float ftemp=0.0;
	Int2 whichfrag=0, x=0;
	Int4 res;
	Char aa;
	FloatLo p,a1,a2,asd,chiw,chiwsd,chi1=0.0,chi2=0.0,chi3=0.0,chi4=0.0;
	Uint4 rotid;
	PRS prsHere;
	Char resname[2];
	Int2 numc;

	tmp_char[0] = '\0';
	
	GetTitle(fragres_text, tmp_char, sizeof(tmp_char));
	sscanf(tmp_char, "%d", &itemp);
	res=(Int4)itemp;
	
	GetTitle(fragaa_text, tmp_char, sizeof(tmp_char));
	aa=tmp_char[0];

	GetTitle(fragprob_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for the probability.");
		Select(fragprob_text);
		return;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<0.0 || ftemp>1.0 || err==0)
		{
			Message (MSG_OK, "Probability must be between zero and one");
			Select(fragprob_text);
			return;
		}
		p=(FloatLo)ftemp;
	}
	
	GetTitle(fraga1_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for Phi.");
		Select(fraga1_text);
		return;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<=-180.0 || ftemp>180.0 || err==0)
		{
			while (ftemp>180.0)
				ftemp-=360.0;
			while (ftemp<=-180.0)
				ftemp+=360.0;
		}
		a1=(FloatLo)ftemp;
	}

	GetTitle(fraga2_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		if (WALKTYPE==WALK_PHIPSI)
			Message (MSG_OK, "Please enter a value for Psi.");
		else
			Message (MSG_OK, "Please enter a value for Theta.");
		Select(fraga2_text);
		return;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<=-180.0 || ftemp>180.0 || err==0)
		{
			while (ftemp>180.0)
				ftemp-=360.0;
			while (ftemp<=-180.0)
				ftemp+=360.0;
		}
		a2=(FloatLo)ftemp;
	}

	GetTitle(fragasd_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		if (WALKTYPE==WALK_PHIPSI)
			Message (MSG_OK, "Please enter a value for Phi/Psi (SD).");
		else
			Message (MSG_OK, "Please enter a value for Phi/Theta (SD).");
		Select(fragasd_text);
		return;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<0.0 || ftemp>45.0 || err==0) {
			if (WALKTYPE==WALK_PHIPSI)
				Message(MSG_OK,"Please use a value between 0 and 45 degrees for Phi/Psi (SD)");
			else
				Message(MSG_OK,"Please use a value between 0 and 45 degrees for Phi/Theta (SD)");
			Select(fragasd_text);
			return;
		}
		asd=(FloatLo)ftemp;
	}

	GetTitle(fragchiw_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for Omega.");
		Select(fragchiw_text);
		return;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<=-180.0 || ftemp>180.0 || err==0)
		{
			while (ftemp>180.0)
				ftemp-=360.0;
			while (ftemp<=-180.0)
				ftemp+=360.0;
		}
		chiw=(FloatLo)ftemp;
	}

	GetTitle(fragchiwsd_text, tmp_char, sizeof(tmp_char));
	if (StringLen(tmp_char)==0)
	{
		Message (MSG_OK, "Please enter a value for Omega (SD).");
		Select(fragchiwsd_text);
		return;
	}
	else {
		err=sscanf(tmp_char, "%f", &ftemp);
		if (ftemp<0.0 || ftemp>45.0 || err==0) {
			Message(MSG_OK,"Please use a value between 0 and 45 degrees for Omega (SD)");
			Select(fragchiwsd_text);
			return;
		}
		chiwsd=(FloatLo)ftemp;
	}

	if (frag_restoadd==0)
		resname[0]=(resinfo->AA);
	else
		resname[0]=Xsequence[frag_restoadd-1];
	resname[1]='\0';
	numc=numchi[StringCSpn(aalist,resname)];
	if (numc>0) {
		GetTitle(fragchi1_text, tmp_char, sizeof(tmp_char));
		if (StringLen(tmp_char)==0)
		{
			Message (MSG_OK, "Please enter a value for Chi 1.");
			Select(fragchi1_text);
			return;
		}
		else {
			if (!StringCmp(tmp_char,"Rotamer Library")) {
				ftemp=0.0;
				err=1;
			}
			else
				err=sscanf(tmp_char, "%f", &ftemp);
			if (ftemp<-180.0 || ftemp>180.0 || err==0) {
				Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi 1");
				Select(fragchi1_text);
				return;
			}
			chi1=(FloatLo)ftemp;
		}
	}
	if (numc>1) {
		GetTitle(fragchi2_text, tmp_char, sizeof(tmp_char));
		if (StringLen(tmp_char)==0)
		{
			Message (MSG_OK, "Please enter a value for Chi 2.");
			Select(fragchi2_text);
			return;
		}
		else {
			if (!StringCmp(tmp_char,"Rotamer Library")) {
				ftemp=0.0;
				err=1;
			}
			else
				err=sscanf(tmp_char, "%f", &ftemp);
			if (ftemp<-180.0 || ftemp>180.0 || err==0) {
				Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi 2");
				Select(fragchi2_text);
				return;
			}
			else if ((resinfo->AA=='D' || resinfo->AA=='F' || resinfo->AA=='Y') && (ftemp<-90.0 || ftemp>90.0)) {
				Message(MSG_OK,"Please use a value between -90 and 90 degrees for chi 2 of D, F or Y");
				Select(fragchi2_text);
				return;
			}
			chi2=(FloatLo)ftemp;
		}
	}
	if (numc>2) {
		GetTitle(fragchi3_text, tmp_char, sizeof(tmp_char));
		if (StringLen(tmp_char)==0)
		{
			Message (MSG_OK, "Please enter a value for Chi 3.");
			Select(fragchi3_text);
			return;
		}
		else {
			if (!StringCmp(tmp_char,"Rotamer Library")) {
				ftemp=0.0;
				err=1;
			}
			else
				err=sscanf(tmp_char, "%f", &ftemp);
			if (ftemp<-180.0 || ftemp>180.0 || err==0) {
				Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi 3");
				Select(fragchi3_text);
				return;
			}
			chi3=(FloatLo)ftemp;
		}
	}
	if (numc>3) {
		GetTitle(fragchi4_text, tmp_char, sizeof(tmp_char));
		if (StringLen(tmp_char)==0)
		{
			Message (MSG_OK, "Please enter a value for Chi 4.");
			Select(fragchi4_text);
			return;
		}
		else {
			if (!StringCmp(tmp_char,"Rotamer Library")) {
				ftemp=0.0;
				err=1;
			}
			else
				err=sscanf(tmp_char, "%f", &ftemp);
			if (ftemp<-180.0 || ftemp>180.0 || err==0) {
				Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi 4");
				Select(fragchi4_text);
				return;
			}
			chi4=(FloatLo)ftemp;
		}
	}
								
	prsHere=(PRS)MemNew(sizeof(RS));
	prsHere->Chi1=chi1;
	prsHere->Chi2=chi2;
	prsHere->Chi3=chi3;
	prsHere->Chi4=chi4;
	rotid=ComputeRotid(prsHere);
	prsHere=MemFree(prsHere);
	
	if (editfrag==TRUE)
	{
		whichfrag = GetValue(fraglist);
		pfdsHere = resinfo->pfdsFragmentHead;
		for (x=1; x<whichfrag; x++)
			pfdsHere = pfdsHere->nextfrag;
		whichfrag = GetValue(fragreslist);
		for (x=1; x<whichfrag; x++)
			pfdsHere = pfdsHere->next;
		pfdsHere->Angle1=a1;
		pfdsHere->Angle2=a2;
		pfdsHere->AngleSD=asd;
		pfdsHere->ChiWMean=chiw;
		pfdsHere->ChiWSD=chiwsd;
		pfdsHere->tout=(Int4)(2.0*sqrt((FloatLo)chiwsd)*(FloatLo)asd);
		if (pfdsHere->tout>BACKTRACK_TRIES_MAX)
		  pfdsHere->tout=BACKTRACK_TRIES_MAX;
		if (pfdsHere->tout<BACKTRACK_TRIES_MIN)
		  pfdsHere->tout=BACKTRACK_TRIES_MIN;
		pfdsHere->rotid=rotid;
		pfdsHere->prob=p;
		if (GetValue(fragreslist)==1) {
			/* propogate prob down list */
			while (pfdsHere!=NULL) {
				pfdsHere->prob=p;
				pfdsHere = pfdsHere->next;
			}
		}	
	}
	
	else
	{
		if (frag_restoadd==0) {
			pfdsHere=NULL;
			if (AddFragmentResidue(&pfdsHere,res,aa,a1,a2,asd,chiw,chiwsd,rotid,p)==NULL)
				Message (MSG_OK, "The new residue of the fragment was not correctly added. \n Please refer to the VisTraj Help File for further information.");
			else {
				if (AddFragmentList(resinfo,pfdsHere)!=ERR_SUCCESS)
					Message (MSG_OK, "The new fragment was not correctly added. \n Please refer to the VisTraj Help File for further information.");
			}
		}
		else {
			pfdstmp = resinfo->pfdsFragmentHead;
			whichfrag = GetValue(fraglist);
			for(x=1; x<whichfrag; x++)
				pfdstmp = pfdstmp->nextfrag;
			if (AddFragmentResidue(&pfdstmp,res,aa,a1,a2,asd,chiw,chiwsd,rotid,p)==NULL)
				Message (MSG_OK, "The new residue of the fragment was not correctly added. \n Please refer to the VisTraj Help File for further information.");		
		}		
	}

	Vtrj_SaveResidue();
	FragBeenEdited=TRUE;
	
	if (frag_restoadd==0 && editfrag==FALSE) {
		Enable(fraglist_doc);
		Enable(fraglist);
	}
	else {
		Enable(frag_doc);
		Enable(fragreslist);
	}

	LoadFragList();
	LoadFragListProc(i);
	LoadFragHere();
	LoadFragProc(i);
	DoneAddFragProc(i);
}


static void DoFragProc(IteM i,Boolean NewList)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3, Vtrj_grp4, Vtrj_grp5, Vtrj_grp6;
	PFDS pfdstmp=NULL;
	Int2 whichfrag=0, x=0;
	Char temp_str[20];
	Int2 numc;
	Char resname[2];
	PRS prsHere;
	
	if (WALKTYPE==WALK_CA && resinfo->resnum==1) {
		Message(MSG_OK,"Sorry, you may not add fragments to the first residue");
		return;
	}
		
	addfrag_win = MovableModalWindow(-50,-33,-10,-10, "", (WndActnProc) DoneAddFragProc);

	Reset(addfrag_win);
	
	frag_restoadd=0;
	if (editfrag==TRUE)
	{
		SetTitle(addfrag_win, "Edit Residue Fragment");
		temp_str[0] = '\0';
		pfdstmp = resinfo->pfdsFragmentHead;
		whichfrag = GetValue(fraglist);
		for(x=1; x<whichfrag; x++)
			pfdstmp = pfdstmp->nextfrag;
		whichfrag = GetValue(fragreslist);
		for(x=1; x<whichfrag; x++)
			pfdstmp = pfdstmp->next;
	}

	if (editfrag==FALSE)
	{
		SetTitle(addfrag_win, "Add New Residue Fragment");
	}

	prsHere=(PRS)MemNew(sizeof(RS));
	if (editfrag==TRUE)
	  GetChiFromRotid(&prsHere,pfdstmp->rotid);
	else
	  GetChiFromRotid(&prsHere,0);
	
	g = HiddenGroup(addfrag_win,-1,-1, (GrpActnProc) NULL);
	
	Vtrj_grp1 =  HiddenGroup(g,1,0,(GrpActnProc) NULL);
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 5, 5);			
	
	Vtrj_grp3 = NormalGroup(Vtrj_grp1, 4, 0, "Residue Properties", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp3, 10, 10);			
	
	if (editfrag==TRUE)
		sprintf(temp_str, "%d", pfdstmp->resnum);
	else {
		if (NewList==TRUE)
			sprintf(temp_str, "%d", resinfo->resnum);
		else {
			pfdstmp = resinfo->pfdsFragmentHead;
			whichfrag = GetValue(fraglist);
			for(x=1; x<whichfrag; x++)
				pfdstmp = pfdstmp->nextfrag;
			while (pfdstmp->next!=NULL)
				pfdstmp=pfdstmp->next;
			frag_restoadd=(pfdstmp->resnum)+1;
			sprintf(temp_str, "%d", frag_restoadd);			
		}
	}
	if (editfrag==TRUE)
		resname[0]=pfdstmp->AA;
	else {
		if (NewList)
			resname[0]=(resinfo->AA);
		else
			resname[0]=Psequence[frag_restoadd-1];
	}
	resname[1]='\0';
	numc=numchi[StringCSpn(aalist,resname)];
	
	StaticPrompt (Vtrj_grp3, "Residue #", 0, dialogTextHeight, systemFont, 'l');    	
	fragres_text = DialogText (Vtrj_grp3, temp_str,  5, (TxtActnProc) NULL);
	Disable(fragres_text);
	
	if (editfrag==TRUE)
		sprintf(temp_str, "%d", pfdstmp->length);
	else
		sprintf(temp_str, "0");
	StaticPrompt (Vtrj_grp3, "Length Remaining", 0, dialogTextHeight, systemFont, 'l');    	
	fraglen_text = DialogText (Vtrj_grp3, temp_str,  5, (TxtActnProc) NULL);
	Disable(fraglen_text);
	
	if (editfrag==TRUE)
		sprintf(temp_str, "%c", pfdstmp->AA);
	else {
		if (NewList==TRUE)
			sprintf(temp_str, "%c", resinfo->AA);
		else
			sprintf(temp_str, "%c", Xsequence[frag_restoadd-1]);
	}
	StaticPrompt (Vtrj_grp3, "Amino Acid", 0, dialogTextHeight, systemFont, 'l');    	
	fragaa_text = DialogText (Vtrj_grp3, temp_str,  5, (TxtActnProc) NULL);
	Disable(fragaa_text);
	
	if (editfrag==TRUE)
		sprintf(temp_str, "%1.2f",pfdstmp->prob);
	else {
		if (NewList==FALSE)
			sprintf(temp_str, "%1.2f",pfdstmp->prob);
		else
			sprintf(temp_str, "1.0");
	}
	if (!StringCmp(temp_str,"0.00") && pfdstmp->prob>0.0) {
		sprintf(temp_str, "%1.2e",pfdstmp->prob);
	}
	StaticPrompt (Vtrj_grp3, "Probability", 0, dialogTextHeight, systemFont, 'l');    	
	fragprob_text = DialogText (Vtrj_grp3, temp_str,  5, (TxtActnProc) NULL);
	if (NewList==FALSE)
		Disable(fragprob_text);
	if (editfrag==TRUE && GetValue(fragreslist)==1)
		Enable(fragprob_text);

	Vtrj_grp4 = NormalGroup(Vtrj_grp1, 1, 0, "Backbone Dihedrals", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp4, 3, 3);			
	
	Vtrj_grp6 = HiddenGroup(Vtrj_grp4, 4, 0, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp6, 10, 10);			

	StaticPrompt (Vtrj_grp6, "Phi", 0, dialogTextHeight, systemFont, 'l');    	
	fraga1_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);
	
 if (WALKTYPE==WALK_PHIPSI)
		StaticPrompt (Vtrj_grp6, "Psi", 0, dialogTextHeight, systemFont, 'l');
	else
		StaticPrompt (Vtrj_grp6, "Theta", 0, dialogTextHeight, systemFont, 'l');
	fraga2_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	if (WALKTYPE==WALK_PHIPSI)
		StaticPrompt (Vtrj_grp6, "Phi/Psi (SD)", 0, dialogTextHeight, systemFont, 'l');    	
	else
		StaticPrompt (Vtrj_grp6, "Phi/Theta (SD)", 0, dialogTextHeight, systemFont, 'l');    	
	fragasd_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp6, "Omega", 0, dialogTextHeight, systemFont, 'l');    	
	fragchiw_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	StaticPrompt (Vtrj_grp6, "Omega (SD)", 0, dialogTextHeight, systemFont, 'l');    	
	fragchiwsd_text = DialogText (Vtrj_grp6, "",  5, (TxtActnProc) NULL);

	Vtrj_grp5 = NormalGroup(Vtrj_grp1, 4, 0, "Sidechain Dihedrals", systemFont, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp5, 10, 10);			

	if (numc>0) {
		StaticPrompt (Vtrj_grp5, "Chi 1", 0, dialogTextHeight, systemFont, 'l');    	
		fragchi1_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);
	}
	if (numc>1) {
		StaticPrompt (Vtrj_grp5, "Chi 2", 0, dialogTextHeight, systemFont, 'l');    	
		fragchi2_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);
	}
	if (numc>2) {
		StaticPrompt (Vtrj_grp5, "Chi 3", 0, dialogTextHeight, systemFont, 'l');    	
		fragchi3_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);
	}
	if (numc>3) {
		StaticPrompt (Vtrj_grp5, "Chi 4", 0, dialogTextHeight, systemFont, 'l');    	
		fragchi4_text = DialogText (Vtrj_grp5, "",  5, (TxtActnProc) NULL);
	}
	if (numc>0) {
		StaticPrompt (Vtrj_grp5, " (0 = ", 0, dialogTextHeight, systemFont, 'l');    	
		StaticPrompt (Vtrj_grp5, "Rotamer Library)", 0, dialogTextHeight, systemFont, 'l');    	
	}
		
	Vtrj_grp2 = HiddenGroup(g,2,-1,(GrpActnProc) NULL);
	
	if (editfrag==FALSE)
		PushButton(Vtrj_grp2, "Add", (BtnActnProc) ApplyFragProc);
	else
		PushButton(Vtrj_grp2, "Apply", (BtnActnProc) ApplyFragProc);
	PushButton(Vtrj_grp2, "Cancel", (BtnActnProc) DoneAddFragProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

	if (editfrag==TRUE)
	{
		sprintf(temp_str, "%d", pfdstmp->resnum);
		SetTitle(fragres_text, temp_str);
		
		sprintf(temp_str, "%d", pfdstmp->length);
		SetTitle(fraglen_text, temp_str);
		
		sprintf(temp_str, "%6.2f", pfdstmp->Angle1);
		SetTitle(fraga1_text, temp_str);
		
		sprintf(temp_str, "%6.2f", pfdstmp->Angle2);
		SetTitle(fraga2_text, temp_str);
		
		sprintf(temp_str, "%5.2f", pfdstmp->AngleSD);
		SetTitle(fragasd_text, temp_str);
		
		sprintf(temp_str, "%6.2f", pfdstmp->ChiWMean);		
		SetTitle(fragchiw_text, temp_str);
		
		sprintf(temp_str, "%5.2f", pfdstmp->ChiWSD);
		SetTitle(fragchiwsd_text, temp_str);
		
		sprintf(temp_str, "%c", pfdstmp->AA);
		SetTitle(fragaa_text, temp_str);
		
		sprintf(temp_str, "%5.3f", pfdstmp->prob);
		SetTitle(fragprob_text, temp_str);
		
		if (numc>0) {
			if (prsHere->Chi1==0.0)
				sprintf(temp_str, "Rotamer Library");
			else
				sprintf(temp_str, "%6.2f", prsHere->Chi1);
			SetTitle(fragchi1_text, temp_str);
		}		
		if (numc>1) {
			if (prsHere->Chi2==0.0)
				sprintf(temp_str, "Rotamer Library");
			else
				sprintf(temp_str, "%6.2f", prsHere->Chi2);
			SetTitle(fragchi2_text, temp_str);
		}
		if (numc>2) {
			if (prsHere->Chi3==0.0)
				sprintf(temp_str, "Rotamer Library");
			else		
				sprintf(temp_str, "%6.2f", prsHere->Chi3);
			SetTitle(fragchi3_text, temp_str);
		}
		if (numc>3) {
			if (prsHere->Chi4==0.0)
				sprintf(temp_str, "Rotamer Library");
			else		
				sprintf(temp_str, "%6.2f", prsHere->Chi4);
			SetTitle(fragchi4_text, temp_str);
		}
	}

	prsHere=MemFree(prsHere);
#ifndef OS_MSWIN
	Disable(fraglist);
#endif
	fragwinopenadd=TRUE;
	Show(addfrag_win);
	ProcessEvents();	
}

static void AddFragProc(IteM i)
{
	editfrag=FALSE;
	DoFragProc(i,FALSE);
}

static void AddFragListProc(IteM i)
{
	editfrag=FALSE;
	DoFragProc(i,TRUE);
}

static void EditFragProc(IteM i)
{
	editfrag=TRUE;
	DoFragProc(i,FALSE);
}

static void DeleteFragListProc(IteM i)
{

	PFDS pfdstmp=NULL;
	Int2 whichfrag=0, x=0;
	Int2 Vtrj_ans=ANS_NO;
	
	Vtrj_ans = Message (MSG_YN, "Are you sure you would like to delete this fragment?");	

	if (Vtrj_ans==ANS_NO) return;

	pfdstmp = resinfo->pfdsFragmentHead;
	whichfrag = GetValue(fraglist);
	for(x=1; x<whichfrag; x++)
		pfdstmp = pfdstmp->nextfrag;
	
	RemoveFragmentList(resinfo,pfdstmp);

	FragBeenEdited=TRUE;
	Vtrj_SaveResidue();

	Reset(fraglist_doc);
	LoadFragList();
	if (resinfo->pfdsFragmentHead!=NULL)
		LoadFragListProc(i);
}

static void LoadFragHere(void)
{

	PFDS pfdstmp=NULL;
	Char res1_str[10], frag_title[40];
	Int2 res1=0, whichfrag=0, x, whichres;

	Disable(fragE_bttn);
	Disable(fragD_bttn);

	DishdgColFmt.font = ParseFont ("Times,14,b");
	DissubColFmt.font = ParseFont ("Helvetica,11,b");
	DistxtColFmt.font = ParseFont ("fixed,10");
/*	lstColFmt.font = programFont;
	tblColFmt.font = programFont;*/

	whichres=GetValue(fragreslist);
	if (whichres==0)
			whichres=1;
	Reset(fragreslist);

	pfdstmp = resinfo->pfdsFragmentHead;

	whichfrag = GetValue(fraglist);

	for(x=1; x<whichfrag; x++)
		pfdstmp = pfdstmp->nextfrag;
		
	while(pfdstmp!=NULL)
	{
		frag_title[0] = '\0';
		res1 = pfdstmp->resnum;

		sprintf(res1_str, "%d", res1);
		
		StringCat(frag_title, res1_str);

		ListItem(fragreslist, frag_title);

		pfdstmp = pfdstmp->next;
	}
	SetValue(fragreslist,whichres);
}

static void DeleteFragProc(IteM i)
{

	PFDS pfdstmp=NULL,pfdstmp2;
	Int2 whichfrag=0, x=0;
	Int2 Vtrj_ans=ANS_NO;

	
	Vtrj_ans = Message (MSG_YN, "Are you sure you would like to delete this fragment element?");	

	if (Vtrj_ans==ANS_NO) return;

	pfdstmp = resinfo->pfdsFragmentHead;
	whichfrag = GetValue(fraglist);
	for(x=1; x<whichfrag; x++)
		pfdstmp = pfdstmp->nextfrag;
	pfdstmp2=pfdstmp;
	whichfrag = GetValue(fragreslist);
	for(x=1; x<whichfrag; x++)
		pfdstmp2 = pfdstmp2->next;
	if (pfdstmp2->next!=NULL) {
		Message(MSG_OK, "Sorry, you may only delete the last residue in the fragment");
		return;
	}

	RemoveFragmentResidue(resinfo,pfdstmp,pfdstmp2->resnum);

	FragBeenEdited=TRUE;
	Vtrj_SaveResidue();

	Reset(frag_doc);
	SetValue(fragreslist,whichfrag-1);
	LoadFragList();
	if (resinfo->pfdsFragmentHead!=NULL)
		LoadFragListProc(i);
	if (whichfrag==1) {		
		DoneFragProc(i);
	}
	else {
		LoadFragHere();
		LoadFragProc(i);
	}
}

static void LoadFragProc(IteM i)
{
	PFDS pfdstmp=NULL;
	Int2 whichfrag=0, x=0;
	Char tmp_str[20];
	PRS prsRotamerHere;
	Char resname[2];
	Int2 numc;

	tmp_str[0] = '\0';

	Enable(fragE_bttn);
	Enable(fragD_bttn);

	Reset(frag_doc);


	pfdstmp = resinfo->pfdsFragmentHead;

	whichfrag = GetValue(fraglist);
	for(x=1; x<whichfrag; x++)
		pfdstmp = pfdstmp->nextfrag;
	
	whichfrag = GetValue(fragreslist);
	for(x=1; x<whichfrag; x++)
		pfdstmp = pfdstmp->next;
	if (pfdstmp->next!=NULL)
		Disable(fragD_bttn);		
	prsRotamerHere=(PRS)MemNew(sizeof(RS));
		
	sprintf(tmp_str, "%d", pfdstmp->resnum);
	AppendText (frag_doc, "Residue:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%d", pfdstmp->length);
	AppendText (frag_doc, "Length Remaining:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);
	
	sprintf(tmp_str, "%c", pfdstmp->AA);
	AppendText (frag_doc, "Amino Acid:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%1.1f", pfdstmp->Angle1);
	AppendText (frag_doc, "Phi:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%1.1f", pfdstmp->Angle2);
	if (WALKTYPE==WALK_PHIPSI)
		AppendText (frag_doc, "Psi:", &DishdgParFmt, &DissubColFmt, programFont);
	else {
		AppendText (frag_doc, "Theta:", &DishdgParFmt, &DissubColFmt, programFont);
	}
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%1.1f", pfdstmp->AngleSD);
	if (WALKTYPE==WALK_PHIPSI)
		AppendText (frag_doc, "Phi/Psi (SD):", &DishdgParFmt, &DissubColFmt, programFont);
	else
		AppendText (frag_doc, "Phi/Theta (SD):", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%1.1f", pfdstmp->ChiWMean);
	AppendText (frag_doc, "Omega:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%1.1f", pfdstmp->ChiWSD);
	AppendText (frag_doc, "Omega (SD):", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	resname[0]=pfdstmp->AA;
	resname[1]='\0';
	numc=numchi[StringCSpn(aalist,resname)];
  GetChiFromRotid(&prsRotamerHere,pfdstmp->rotid);
  if (numc>0) {
	  if (prsRotamerHere->Chi1!=0.0)
			sprintf(tmp_str, "%1.1f", prsRotamerHere->Chi1);
		else
			StringCpy(tmp_str,"Rotamer Library");
		AppendText (frag_doc, "Chi 1:", &DishdgParFmt, &DissubColFmt, programFont);
		AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);	
	}
  if (numc>1) {
	  if (prsRotamerHere->Chi2!=0.0)
			sprintf(tmp_str, "%1.1f", prsRotamerHere->Chi2);
		else
			StringCpy(tmp_str,"Rotamer Library");
		AppendText (frag_doc, "Chi 2:", &DishdgParFmt, &DissubColFmt, programFont);
		AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);	
	}
  if (numc>2) {
	  if (prsRotamerHere->Chi3!=0.0)
			sprintf(tmp_str, "%1.1f", prsRotamerHere->Chi3);
		else
			StringCpy(tmp_str,"Rotamer Library");
		AppendText (frag_doc, "Chi 3:", &DishdgParFmt, &DissubColFmt, programFont);
		AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);	
	}
  if (numc>3) {
	  if (prsRotamerHere->Chi4!=0.0)
			sprintf(tmp_str, "%1.1f", prsRotamerHere->Chi4);
		else
			StringCpy(tmp_str,"Rotamer Library");
		AppendText (frag_doc, "Chi 4:", &DishdgParFmt, &DissubColFmt, programFont);
		AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);	
	}
		
	sprintf(tmp_str, "%d", pfdstmp->tout);
	AppendText (frag_doc, "# Tries:", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);

	sprintf(tmp_str, "%1.2f", pfdstmp->pSS);
	AppendText (frag_doc, "p(SS):", &DishdgParFmt, &DissubColFmt, programFont);
	AppendText (frag_doc, tmp_str, &DishdgParFmt, &DistxtColFmt, programFont);			
		
  prsRotamerHere=MemFree(prsRotamerHere);
	UpdateDocument (frag_doc, 0, 0);
	
}

static void LoadFragListProc(IteM i)
{
	PFDS pfdstmp=NULL;
	Int2 whichfrag=0, x=0;
	Char tmp_str[20];
	
	pfdstmp = resinfo->pfdsFragmentHead;
	if (pfdstmp==NULL)
		return;

	tmp_str[0] = '\0';

	Enable(fraglistE_bttn);
	Enable(fraglistD_bttn);

	Reset(fraglist_doc);

	whichfrag = GetValue(fraglist);

	for (x=1; x<whichfrag; x++)
		pfdstmp = pfdstmp->nextfrag;
	
	sprintf(tmp_str, "%1.2f", pfdstmp->prob);
	if (!StringCmp(tmp_str,"0.00") && pfdstmp->prob>0.0) {
		sprintf(tmp_str, "%1.2e",pfdstmp->prob);
	}
/*	AppendText (fraglist_doc, "Probability:", &DishdgParFmt, &DissubColFmt, programFont);*/
	AppendText (fraglist_doc, tmp_str, &DishdgParFmt, &DissubColFmt, programFont);
	
	UpdateDocument (fraglist_doc, 0, 0);
	
}

static void LoadFragList(void)
{

	PFDS pfdstmp=NULL,pfdstmp2=NULL;
	Char buf[PATH_MAX];
	Char res1_str[10], res2_str[10], frag_title[40];
	Int2 res1=0, res2=0;
	Int2 whichfrag,listsize=0;

	Disable(fraglistE_bttn);
	Disable(fraglistD_bttn);

	DishdgColFmt.font = ParseFont ("Times,14,b");
	DissubColFmt.font = ParseFont ("Helvetica,11,b");
	DistxtColFmt.font = ParseFont ("fixed,10");
/*	lstColFmt.font = programFont;
	tblColFmt.font = programFont;*/

	whichfrag = GetValue(fraglist);
	if (whichfrag==0)
		whichfrag=1;
	Reset(fraglist);

	pfdstmp = resinfo->pfdsFragmentHead;
	
	if (pfdstmp==NULL)
		hasafrag[res_num-1]=FALSE;
	else
		hasafrag[res_num-1]=TRUE;
	Hide(aaname_panel);
	Show(aaname_panel);
	if (pfdstmp==NULL)
	{
		Reset(fraglist_doc);
		sprintf(buf,"There are no fragments attached to residue %d.",resinfo->resnum);
		AppendText (fraglist_doc, buf, &hdgParFmt, &FraghdgColFmt, programFont);
		UpdateDocument (fraglist_doc, 0, 0);
		Disable(fraglist);
		Disable(fraglist_doc);
		Disable(fraglistE_bttn);
		Disable(fraglistD_bttn);
	}

	while(pfdstmp!=NULL)
	{
		frag_title[0] = '\0';
		res1 = pfdstmp->resnum;
		pfdstmp2=pfdstmp;
		while (pfdstmp2->next!=NULL)
			pfdstmp2=pfdstmp2->next;
		res2 = pfdstmp2->resnum;

		sprintf(res1_str, "%d", res1);
		sprintf(res2_str, "%d", res2);
		
		StringCat(frag_title, res1_str);
		StringCat(frag_title, " - ");
		StringCat(frag_title, res2_str);

		ListItem(fraglist, frag_title);

		pfdstmp = pfdstmp->nextfrag;
		listsize++;
	}
	if (whichfrag>listsize)
		whichfrag=listsize;
	SetValue(fraglist,whichfrag);
}

static void EditFragListProc(IteM i)
{

	GrouP g, Vtrj_grp1, Vtrj_grp2;

	editfrag_win = MovableModalWindow(-50,-33,-10,-10, "Edit Residue Fragment", (WndActnProc) DoneFragProc);

	g = HiddenGroup(editfrag_win,-1,-1, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);

	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 20, 20);
	StaticPrompt (Vtrj_grp1, "Residue", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp1, "Dihedrals", 0, dialogTextHeight, systemFont, 'l');
	fragreslist = SingleList (Vtrj_grp1,10,15, (LstActnProc) LoadFragProc);

	frag_doc = DocumentPanel (Vtrj_grp1, stdCharWidth * 20, stdLineHeight * 15);

	Vtrj_grp2 = HiddenGroup(g,4,-1,(GrpActnProc) NULL);

	SetGroupSpacing(g, 15, 15);

	fragE_bttn = PushButton(Vtrj_grp2, "Edit", (BtnActnProc) EditFragProc);
	fragA_bttn = PushButton(Vtrj_grp2, "Append", (BtnActnProc) AddFragProc);
	fragD_bttn = PushButton(Vtrj_grp2, "Delete", (BtnActnProc) DeleteFragProc);
	PushButton(Vtrj_grp2, "Close", (BtnActnProc) DoneFragProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

	LoadFragHere();
	LoadFragProc(i);

	fragwinopenres=TRUE;
#ifndef OS_MSWIN
	Disable(fraglist);
#endif
	Show(editfrag_win);
	ProcessEvents();

}

static void DisplayFragListProc (IteM i)
{

	GrouP g, Vtrj_grp1, Vtrj_grp2;

	editfraglist_win = MovableModalWindow(-50,-33,-10,-10, "Residue Fragments", (WndActnProc) DoneFragListProc);

	g = HiddenGroup(editfraglist_win,-1,-1, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);

	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 20, 20);

	StaticPrompt (Vtrj_grp1, "Residue Range", 0, dialogTextHeight, systemFont, 'l');
	StaticPrompt (Vtrj_grp1, "Probability", 0, dialogTextHeight, systemFont, 'l');
	fraglist = SingleList (Vtrj_grp1,10,15, (LstActnProc) LoadFragListProc);

	fraglist_doc = DocumentPanel (Vtrj_grp1, stdCharWidth * 8, stdLineHeight * 15);

	Vtrj_grp2 = HiddenGroup(g,4,-1,(GrpActnProc) NULL);

	SetGroupSpacing(g, 15, 15);

	fraglistE_bttn = PushButton(Vtrj_grp2, "Edit", (BtnActnProc) EditFragListProc);
	fraglistA_bttn = PushButton(Vtrj_grp2, "Add", (BtnActnProc) AddFragListProc);
	fraglistD_bttn = PushButton(Vtrj_grp2, "Delete", (BtnActnProc) DeleteFragListProc);
	PushButton(Vtrj_grp2, "Close", (BtnActnProc) DoneFragListProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

	LoadFragList();
	LoadFragListProc(i);

	Show(editfraglist_win);
	ProcessEvents();
}

static void ImportFragListProc (IteM i)
{
	Char path[PATH_MAX];
	Char buf[MAXCOL];
	Char resname[2];
	int istart,ilen;
	float fprob,phi,psi,sd,omega,omegasd,chi1,chi2,chi3,chi4;
	PFDS pfdsHere;
	PTGS ptgs;
	Int2 err,res,numc;
	Int4 line=0,rotid=0,numrecords;
	PRS prsHere;
	FILE *fp;
	Boolean ulr;

	path[0] = '\0';

	if (!GetInputFileName(path, sizeof(path), "txt", "Text Files"))
		return;

	if ((fp=FileOpen(path,"r"))==NULL)
		return;
	while (FileGets(buf,MAXCOL,fp)!=NULL) {
		line++;
		if (sscanf(buf,"%d %d %f",&istart,&ilen,&fprob)==3) {
			if (fprob<0.0) {
				fprob=0.0;
				Message (MSG_OK, "Warning: Probability at line %d truncated to zero",line);
			}
			if (fprob>1.0) {
				fprob=1.0;
				Message (MSG_OK, "Warning: Probability at line %d truncated to one",line);
			}
			res=istart-1;
			pfdsHere=NULL;
			while (res<istart+ilen-1) {
				if (FileGets(buf,MAXCOL,fp)!=NULL) {
					line++;
					chi1=0.0;
					chi2=0.0;
					chi3=0.0;
					chi4=0.0;
					err=sscanf(buf,"%f %f %f %f %f %f %f %f %f",&phi,&psi,&sd,&omega,&omegasd,&chi1,&chi2,&chi3,&chi4);
					if (err>4) {
						res++;
						/* found a valid fragment residue record - add it to fragment */
						/* validate data */
						if (phi<=-180.0 || phi>180.0)
						{
							while (phi>180.0)
								phi-=360.0;
							while (phi<=-180.0)
								phi+=360.0;
						}
						if (psi<=-180.0 || psi>180.0)
						{
							while (psi>180.0)
								psi-=360.0;
							while (psi<=-180.0)
								psi+=360.0;
						}
						if (sd<0.0) {
							sd=0.0;
							Message (MSG_OK, "Warning: Standard Deviation at line %d truncated to 0",line);
						}
						if (sd>45.0) {
							sd=45.0;
							Message (MSG_OK, "Warning: Standard Deviation at line %d truncated to 45",line);
						}
						if (omega<=-180.0 || omega>180.0)
						{
							while (omega>180.0)
								omega-=360.0;
							while (omega<=-180.0)
								omega+=360.0;
						}
						if (omegasd<0.0) {
							omegasd=0.0;
							Message (MSG_OK, "Warning: Omega Standard Deviation at line %d truncated to 0",line);
						}
						if (omegasd>45.0) {
							omegasd=45.0;
							Message (MSG_OK, "Warning: Omega Standard Deviation at line %d truncated to 45",line);
						}
						resname[0]=Xsequence[res-1];
						resname[1]='\0';
						numc=numchi[StringCSpn(aalist,resname)];
						if (numc>0) {
							if (chi1<=-180.0 || chi1>180.0)
							{
								while (chi1>180.0)
									chi1-=360.0;
								while (chi1<=-180.0)
									chi1+=360.0;
							}
						}
						if (numc>1) {
							if (chi2<=-180.0 || chi2>180.0)
							{
								while (chi2>180.0)
									chi2-=360.0;
								while (chi2<=-180.0)
									chi2+=360.0;
							}
							if ((resname[0]=='D' || resname[0]=='F' || resname[0]=='Y') && (chi2<-90.0 || chi2>90.0)) {
								if (chi2<-90.0)
									chi2+=180.0;
								if (chi2>90.0)
									chi2-=180.0;
							}
						}
						if (numc>2) {
							if (chi3<=-180.0 || chi3>180.0)
							{
								while (chi3>180.0)
									chi3-=360.0;
								while (chi3<=-180.0)
									chi3+=360.0;
							}
						}
						if (numc>3) {
							if (chi4<=-180.0 || chi4>180.0)
							{
								while (chi4>180.0)
									chi4-=360.0;
								while (chi4<=-180.0)
									chi4+=360.0;
							}
						}

						prsHere=(PRS)MemNew(sizeof(RS));
						prsHere->Chi1=chi1;
						prsHere->Chi2=chi2;
						prsHere->Chi3=chi3;
						prsHere->Chi4=chi4;
						rotid=ComputeRotid(prsHere);
						prsHere=MemFree(prsHere);
						if (AddFragmentResidue(&pfdsHere,res,resname[0],phi,psi,sd,omega,omegasd,rotid,(FloatLo)fprob)==NULL) {
							Message (MSG_OK, "The new residue of the fragment was not correctly added. \n Please refer to the VisTraj Help File for further information.");
						}
					}
				}
				else
					break;
			}
			if (pfdsHere) {
				if ((Int2)istart==resinfo->resnum) {
					ptgs=resinfo;
				}
				else {
	                ulr=USE_LOTS_RAM;
	                USE_LOTS_RAM=FALSE;
					TGInit(tmpdbasename, DB_READ, &numrecords);
					ptgs =  TrajGraphRead((Int2) istart);
				}
				if (AddFragmentList(ptgs,pfdsHere)!=ERR_SUCCESS)
					Message (MSG_OK, "The new fragment was not correctly added. \n Please refer to the VisTraj Help File for further information.");
				else
					hasafrag[(Int2)istart-1]=TRUE;
				if ((Int2)istart!=resinfo->resnum) {
					TrajGraphWrite (ptgs, USE_RLE, TRUE);
					TGClose();
	                USE_LOTS_RAM=ulr;
	                ptgs=FreeTraj(ptgs);
				}
				else
					Vtrj_SaveResidue();
			}
		}
	}

	FileClose(fp);
	FragBeenEdited=TRUE;
	Hide(aaname_panel);
	Show(aaname_panel);

/*	Vtrj_SaveResidue();
	LoadFragList();
	LoadFragListProc(i);
	LoadFragHere();
	LoadFragProc(i);
	DoneAddFragProc(i);*/

}


/*****************************************************************
										MAKESHIFT TOOLTIPS
******************************************************************/

static void LoadInfoTextProc (Int2 Vtrj_topic, Int2 Vtrj_heading)
{

	Char printedstring[5000];
    Char inbuf[MAXCOL];
	Char fnam[PATH_MAX];
    Int2 topicnum=0, headingnum=0, Vtrj_len=0;
	FILE *f;
	Boolean shouldloadit=FALSE;
	Char onespace[]=" ";

	printedstring[0]='\0';
	inbuf[0]='\0';


	sprintf(fnam,"%s%s",CFG_local_datafilepath,VTRJ_HELP_FILE);
	if((f=FileOpen(fnam,"r"))==NULL)
		ErrPostEx(SEV_ERROR,1,1,"Unable to open VisTraj Help file %s", fnam);

	while (FileGets(inbuf, MAXCOL,f)!=NULL)
	{
		if (shouldloadit==TRUE)
		{
			if (inbuf[0] =='#')
			{
				topicnum++;
			}
			else if(inbuf[0]!='!' && topicnum==Vtrj_topic) 
			{
				Vtrj_len = StringLen(inbuf);
				StringNCat(printedstring,inbuf,Vtrj_len-1);
				StringCat(printedstring,onespace);
			}
		}
		else if (inbuf[0] == '!')
		{
				headingnum++;
					
				if (headingnum==Vtrj_heading) shouldloadit = TRUE;
		}
	}
	FileClose(f);
	if (Vtrj_heading==2) StringCpy(descrAA[tmp_Glob], printedstring);
/*	if (Vtrj_heading==2) SetTitle(newAA_info, printedstring);*/
	if (Vtrj_heading==4) SetTitle(REdit_Info, printedstring);
	if (Vtrj_heading==5)
	{
		if (isNew==FALSE) SetTitle(GEdit_Info, printedstring);
		if (isNew==TRUE) SetTitle(GNewEdit_Info, printedstring);
	}

}

/*****************************************************************
                                        GLOBAL VARIABLE EDIT FUNCTIONS
*****************************************************************/

static void SaveOrigGlobs (void)
{
	if (OrigGlobs!=NULL) OrigGlobs = MemFree(OrigGlobs);
	OrigGlobs = (VOGPtr) MemNew(sizeof(VOG));
	
	OrigGlobs->bbet = BACKBONE_ERROR_TOLERANCE;
	OrigGlobs->bbp = BACKBONE_PRECISION;
	OrigGlobs->nrt = NUM_ROT_TRIES;
	OrigGlobs->abbb = ATOM_BOUNCINESS_BB;
	OrigGlobs->absc = ATOM_BOUNCINESS_SC;
	OrigGlobs->is = INCSIZE;
	OrigGlobs->tp = TUNNEL_PROB;
	OrigGlobs->sbb = START_BACKBONE;
	OrigGlobs->tgu = TGUNITS;
	OrigGlobs->wt = WALKTYPE;
	strcpy(OrigGlobs->cf, CONSTRAINT_FILE);
	OrigGlobs->tt = TRAJTYPE;
	OrigGlobs->td = TRAJDIV;
	OrigGlobs->to = TIMEOUT;
	OrigGlobs->bch = BUMPCHECK_HYDROGEN;
	OrigGlobs->msf = MARKOV_SCALE_FACTOR;
}

static void LoadGEditProc (IteM i)
{	

	Boolean isenab = FALSE;
		
	Char printwhat[100];
	
/*	SetTitle(currentGvalue_text, "");*/
	
	editwhatG = GetValue(editGlist);
	
	switch (editwhatG)
	{
		case 1:
			if(!StringCmp(NewGlobs->wt,"PHIPSI")) SetValue(walktype_popup, 2);
			else if(!StringCmp(NewGlobs->wt,"CA")) SetValue(walktype_popup, 3);
			else SetValue(walktype_popup, 1);
			Disable(walktype_popup);
			break;
		case 2:	
			sprintf(printwhat, "%1.1f", NewGlobs->bbet);
			isenab = TRUE;
			break;
		case 3:	
			sprintf(printwhat, "%1.4f", NewGlobs->bbp);
			isenab = TRUE;
			break;
		case 4:	
			sprintf(printwhat, "%d", NewGlobs->nrt);
			isenab = TRUE;
			break;
		case 5:	
			sprintf(printwhat, "%1.2f", NewGlobs->abbb);
			isenab = TRUE;
			break;
		case 6:	
			sprintf(printwhat, "%1.2f", NewGlobs->absc);
			isenab = TRUE;
			break;
		case 7:	
			if (!StringCmp(NewGlobs->bch,"TRUE")) SetValue(bchydrogen_popup, 1);
			else SetValue(bchydrogen_popup, 2);
			break;	
		case 8:	
			sprintf(printwhat, "%1.1f", NewGlobs->is);
			isenab = TRUE;
			break;
		case 9:	
			sprintf(printwhat, "%1.1f", NewGlobs->sbb);
			isenab = TRUE;
			break;
		case 10:	
			sprintf(printwhat, "%d", NewGlobs->to);
			isenab = FALSE;
			break;
		case 11:	
			sprintf(printwhat, "%1.2f", NewGlobs->msf);
			isenab = FALSE;
			break;
		case 12:	
			sprintf(printwhat, "%1.2f", NewGlobs->tp);
			isenab = TRUE;
			break;
		case 13:
			if(TGUNITS==UNITS_ARBITRARY) SetValue(tgunits_popup, 2);
			else if(TGUNITS==UNITS_CREASE) SetValue(tgunits_popup, 3);
			else if(TGUNITS==UNITS_ZHANG) SetValue(tgunits_popup, 4);
			else if(TGUNITS==UNITS_BRYANT) SetValue(tgunits_popup, 5);
			else SetValue(tgunits_popup, 1);
			Disable(tgunits_popup);
			break;
		case 14:	
			sprintf(printwhat, "%s",   CONSTRAINT_FILE);
			isenab = FALSE;
			break;
		case 15:	
			if(TRAJTYPE==TRAJ_NA) SetValue(ttype_popup,1);
			if(TRAJTYPE==TRAJ_UNIFORM) SetValue(ttype_popup,2);
			if(TRAJTYPE==TRAJ_STANDARD) SetValue(ttype_popup,3);
			if(TRAJTYPE==TRAJ_SSTRU) SetValue(ttype_popup,4);
			if(TRAJTYPE==TRAJ_GOR) SetValue(ttype_popup,5);
			Disable(ttype_popup);
			break;
		case 16:	
			sprintf(printwhat, "%d",   TRAJDIV);
			isenab = FALSE;
			break;
	}
	
	if (editwhatG==1)
	{
		Show(walktype_popup);
		Hide(bchydrogen_popup);
		Hide(currentGvalue_text);
		Hide(tgunits_popup);
		Hide(ttype_popup);
	}
	else if (editwhatG==7)
	{
		Hide(walktype_popup);
		Hide(tgunits_popup);
		Hide(ttype_popup);
		Show(bchydrogen_popup);
		Hide(currentGvalue_text);
	}
	else if (editwhatG==13)
	{
		Show(tgunits_popup);
		Hide(bchydrogen_popup);
		Hide(walktype_popup);
		Hide(ttype_popup);
		Hide(currentGvalue_text);
	}
	else if (editwhatG==15)
	{
		Show(ttype_popup);
		Hide(bchydrogen_popup);
		Hide(walktype_popup);
		Hide(tgunits_popup);
		Hide(currentGvalue_text);
	}
	else
	{
		Hide(walktype_popup);
		Hide(ttype_popup);
		Hide(tgunits_popup);
		Hide(bchydrogen_popup);
		Show(currentGvalue_text);
		if (lastg!=editwhatG) {
			SetTitle(currentGvalue_text, printwhat);
			lastg=editwhatG;
		}
		if (isenab == TRUE) Enable(currentGvalue_text);
		if (isenab == FALSE) Disable(currentGvalue_text);

	}
/*	Select(currentGvalue_text);	*/
	LoadInfoTextProc(editwhatG, 5);
}

static void EnableGSaveEditProc (IteM i)
{	
	Int2 editvalue=0,whichone=0;
  Char str[30];
	int itemp;
	float ftemp;

	editvalue = GetValue(editGlist);
  GetTitle(currentGvalue_text, str, sizeof(str));
              	
	switch (editvalue)
	{
	
		case 2:	
			sscanf(str, "%f",  &ftemp);
			NewGlobs->bbet = (FloatHi) ftemp;
			break;
		case 3:	
			sscanf(str, "%f",  &ftemp);
			NewGlobs->bbp = (FloatHi) ftemp;
			break;
		case 4:	
			sscanf(str, "%d",   &itemp);
			NewGlobs->nrt = (Int2) itemp;
			break;
		case 5:	
			sscanf(str, "%f",  &ftemp);
			NewGlobs->abbb = (FloatHi) ftemp;
			break;
		case 6:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->absc = (FloatHi) ftemp;
			break;
		case 7:	
			whichone=GetValue(bchydrogen_popup);
			if (whichone==1) StringCpy(NewGlobs->bch,"TRUE");
			if (whichone==2) StringCpy(NewGlobs->bch,"FALSE");
			break;	
		case 8:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->is = (FloatHi) ftemp;
			break;
		case 9:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->sbb = (FloatHi) ftemp;
			break;
/*		case 14:	
			sscanf(str, "%d",   &itemp);
			TIMEOUT = (Int4) itemp;
			break;*/
		case 12:	
			sscanf(str, "%f", &ftemp);
			NewGlobs->tp = (FloatHi) ftemp;
			break;
	}
	
	LoadGEditProc(i);
	
	Enable (applyG_bttn);
	needApply=TRUE;
}

static void ApplyGEditProc (IteM i)
{

	int itemp = 0;
	float ftemp = 0.00;
	Int2 whichone=0;
		
	Ghasbeenedited = TRUE;

	ftemp=NewGlobs->bbet;
	if (ftemp<10.0) {
		Message(MSG_OK,"A minimum value of 10 is required for Backbone Error Tolerance");
		SetValue(editGlist, 2);
		LoadGEditProc(i);
		return;
	}
	ftemp=NewGlobs->bbp;
	if (ftemp<0.0005 || ftemp>5.0) {
		Message(MSG_OK,"A Backbone Precision must be between 0.0005 and 5.0");
		SetValue(editGlist, 3);
		LoadGEditProc(i);
		return;
	}
	itemp=NewGlobs->nrt;
	if (itemp<1 || itemp>100) {
		Message(MSG_OK,"Rotamer Tries must be between 1 and 100");
		SetValue(editGlist, 4);
		LoadGEditProc(i);
		return;
	}
	ftemp=NewGlobs->abbb;
	if (ftemp<0.0) {
		Message(MSG_OK,"Atom Bounciness must be non-negative");
		SetValue(editGlist, 5);
		LoadGEditProc(i);
		return;
	}
	ftemp=NewGlobs->absc;
	if (ftemp<0.0) {
		Message(MSG_OK,"Atom Bounciness must be non-negative");
		SetValue(editGlist, 6);
		LoadGEditProc(i);
		return;
	}
	ftemp=NewGlobs->is;
	if (ftemp<1.0 || ftemp>89.0) {
		Message(MSG_OK,"Increment Size must be between 1 and 89");
		SetValue(editGlist, 8);
		LoadGEditProc(i);
		return;
	}
	ftemp=NewGlobs->tp;
	if (ftemp<0.0 || ftemp>1.0) {
		Message(MSG_OK,"Tunnelling Probability must be between 0.0 and 1.0");
		SetValue(editGlist, 12);
		LoadGEditProc(i);
		return;
	}
	ftemp=NewGlobs->sbb;
	if (ftemp<=0.0 || ftemp>360.0) {
		Message(MSG_OK,"Backbone Start must be greater than 0 and less than or equal to 360");				
		SetValue(editGlist, 9);
		LoadGEditProc(i);
		return;
	}
					
	BACKBONE_ERROR_TOLERANCE = (FloatHi) NewGlobs->bbet;
	BACKBONE_PRECISION = (FloatHi) NewGlobs->bbp;
	NUM_ROT_TRIES = (Int2) NewGlobs->nrt;
	whichone=GetValue(bchydrogen_popup);
	if(whichone==1) BUMPCHECK_HYDROGEN = TRUE;
	if(whichone==2) BUMPCHECK_HYDROGEN = FALSE;
	ATOM_BOUNCINESS_BB = (FloatHi) NewGlobs->abbb;
	ATOM_BOUNCINESS_SC = (FloatHi) NewGlobs->absc;
	INCSIZE = (FloatHi) NewGlobs->is;
	TUNNEL_PROB = (FloatHi) NewGlobs->tp;
	START_BACKBONE = (FloatHi) NewGlobs->sbb;
	needApply=FALSE;

/*		case 14:	
			sscanf(str, "%d",   &itemp);
			TIMEOUT = (Int4) itemp;
			break;*/
	
	Disable(applyG_bttn);
	LoadGEditProc(i);

}

static void DoneEditGInfoProc (IteM i)
{
	Int2 ans;
	
	if (needApply==TRUE) {
		ans=Message (MSG_YN, "Apply unsaved changes before closing?");
		if (ans==ANS_YES) {
			ApplyGEditProc (i);
			/* check if apply failed */
			if (needApply==TRUE)
				return;		
		}
	}
	lastg=0;
	needApply=FALSE;
	Hide(editGinfo_win);
}

static void EditGInfoProc (IteM i)
{
	
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3;
/*	ButtoN cancelG_bttn;*/
	
	if (NewGlobs!=NULL) NewGlobs = MemFree(NewGlobs);
	NewGlobs = (VNGPtr) MemNew(sizeof(VNG));
	if (WALKTYPE==WALK_CA)
		StringCpy(NewGlobs->wt,"CA");
	else if (WALKTYPE==WALK_PHIPSI)
		StringCpy(NewGlobs->wt,"PHIPSI");
	else
		StringCpy(NewGlobs->wt,"UNKNOWN");
	NewGlobs->to = (Int4)TIMEOUT;
	NewGlobs->bbet = (FloatHi)BACKBONE_ERROR_TOLERANCE;
	NewGlobs->bbp = (FloatHi)BACKBONE_PRECISION;
	NewGlobs->abbb = (FloatHi)ATOM_BOUNCINESS_BB;
	NewGlobs->absc = (FloatHi)ATOM_BOUNCINESS_SC;
	if (BUMPCHECK_HYDROGEN==TRUE)
		StringCpy(NewGlobs->bch,"TRUE");
	else
		StringCpy(NewGlobs->bch,"FALSE");
	NewGlobs->tp = (FloatHi)TUNNEL_PROB;
	NewGlobs->is = (FloatHi)INCSIZE;
	NewGlobs->sbb = (FloatHi)START_BACKBONE;
	NewGlobs->nrt = (Int2)NUM_ROT_TRIES;
	NewGlobs->msf = (FloatHi)MARKOV_SCALE_FACTOR;	
	
	isNew = FALSE;
	needApply = FALSE;

	editGinfo_win = MovableModalWindow(-50,-33,-10,-10, "Edit Global Trajectory Properties", (WndActnProc) DoneEditGInfoProc);

	g = HiddenGroup(editGinfo_win,-1,-1, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,-1,-1,(GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 10, 10);			
	
	editGlist = SingleList (Vtrj_grp1,17,8, (LstActnProc) LoadGEditProc);
		
	ListItem (editGlist,"Trajectory Space");
	ListItem (editGlist,"Backbone Error Tolerance");
	ListItem (editGlist,"Backbone Precision");
	ListItem (editGlist,"Rotamer Tries");
	ListItem (editGlist,"Backbone Atom Bounciness");
	ListItem (editGlist,"Sidechain Atom Bounciness");
	ListItem (editGlist,"Bumpcheck Hydrogen");
	ListItem (editGlist,"Increment Size");
	ListItem (editGlist,"Backbone Start");	
	ListItem (editGlist,"Mean Tries per Residue");
	ListItem (editGlist,"Markov Scale Factor");
	ListItem (editGlist,"Tunnelling Probability");			
	ListItem (editGlist,"Distribution Units");
	ListItem (editGlist,"Constraint File");
	ListItem (editGlist,"Trajectory Distribution Creation Method");	
	ListItem (editGlist,"Discretization Size");		
	
	Vtrj_grp3 = HiddenGroup(Vtrj_grp1, 0, 0, (GrpActnProc) NULL);

	currentGvalue_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) EnableGSaveEditProc);
	Disable(currentGvalue_text);
	walktype_popup = PopupList(Vtrj_grp3, TRUE, (PupActnProc) NULL);
	bchydrogen_popup = PopupList(Vtrj_grp3, TRUE, (PupActnProc) EnableGSaveEditProc);
	ttype_popup = PopupList(Vtrj_grp3, TRUE, (PupActnProc) NULL);
	tgunits_popup = PopupList(Vtrj_grp3, TRUE, (PupActnProc) NULL);

	PopupItem(walktype_popup, "UNKNOWN");
	PopupItem(walktype_popup, "RAMACHANDRAN");
	PopupItem(walktype_popup, "ALPHA-CARBON");	

	PopupItem(bchydrogen_popup, "TRUE");
	PopupItem(bchydrogen_popup, "FALSE");

	PopupItem(ttype_popup, "N\\A");
	PopupItem(ttype_popup, "UNIFORM");
	PopupItem(ttype_popup, "STANDARD");
	PopupItem(ttype_popup, "ONE STATE SECONDARY STRUCTURE");
	PopupItem(ttype_popup, "THREE STATE SECONDARY STRUCTURE");

	PopupItem(tgunits_popup, "UNKNOWN");
	PopupItem(tgunits_popup, "ARBITRARY");
	PopupItem(tgunits_popup, "CREASE");
	PopupItem(tgunits_popup, "ZHANG-DELISI");
	PopupItem(tgunits_popup, "BRYANT-LAWRENCE");

	Hide(walktype_popup);
	Hide(bchydrogen_popup);
	Hide(ttype_popup);
	Hide(tgunits_popup);

/*	Select(currentGvalue_text);*/

	AlignObjects(ALIGN_JUSTIFY, currentGvalue_text, editGlist, NULL);
	AlignObjects(ALIGN_JUSTIFY, walktype_popup, editGlist, NULL);
	AlignObjects(ALIGN_JUSTIFY, bchydrogen_popup, editGlist, NULL);
	AlignObjects(ALIGN_JUSTIFY, ttype_popup, editGlist, NULL);
	AlignObjects(ALIGN_JUSTIFY, tgunits_popup, editGlist, NULL);
	
	
	Vtrj_grp2 = HiddenGroup(g, 2, -1, (GrpActnProc) NULL);
	SetGroupSpacing(Vtrj_grp2, 10, 10);			

	applyG_bttn = PushButton(Vtrj_grp2, "Apply", (BtnActnProc) ApplyGEditProc);
	Disable(applyG_bttn);

	/*cancelG_bttn =*/ DefaultButton(Vtrj_grp2, "Close", (BtnActnProc) DoneEditGInfoProc);
		

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, Vtrj_grp1, NULL);

	Break(g);
	
	GEdit_Info = ScrollText (g, 24, 8, systemFont, TRUE, NULL);

	Show(editGinfo_win);
	ProcessEvents();

}



/*****************************************************************
                                       CYSTINE & PROLINE EDIT FUNCTIONS
*****************************************************************/


static void Vtrj_ChangeFromP(void)
{
	
	Int4 numrecords=0;

	resinfo->CisTrajGraph = MemFree (resinfo->CisTrajGraph);
	resinfo->CisPeak = 0;
	resinfo->CisTrajIntegral = 0;
	
	if (resinfo->ChiWMean > -90.0 && resinfo->ChiWMean < 90.0)
		resinfo->ChiWMean += 180.0;
	
	if (res_num>1)
	{
		TGInit(tmpdbasename, DB_READ, &numrecords);
		UpdatePCis(res_num-1, 0.0);
		TGClose();
	}
		
}

static void Vtrj_ChangeToP(void)
{
	Int4 numrecords=0;

	if (res_num>1)
	{
		TGInit(tmpdbasename, DB_READ, &numrecords);
		UpdatePCis(res_num-1, P_CISP);
		TGClose();
	}
	
	resinfo->CisTrajGraph = (Int4 *) MemNew(sizeof(Int4)*resinfo->dim*resinfo->dim);	
	resinfo->CisTrajIntegral = FillTG(TRAJ_STANDARD, "P", 0, NULL, resinfo->dim, resinfo->CisTrajGraph, TRUE);

}


static void Vtrj_ChangeFromC(void)
{
	resinfo->pSS = 0.0;
	origresinfo->pSS = 0.0;
                reseditinfo->pSS = 0.0;
}

static void Vtrj_ChangeToC(void)
{
	resinfo->pSS = P_SS;
	origresinfo->pSS = P_SS;
                reseditinfo->pSS = P_SS;
	
}


/*****************************************************************
                                       AMINO ACID FUNCTIONS
*****************************************************************/


static void DoneModifyRProc (IteM i)
{
	Hide(modifyR_win);
	ResReplace=TRUE;
	isAAused=FALSE;
}

static void OkmodRProc (IteM i)
{
	Int2 whichAAlist = 0, dowhat=-1;
	Char the_aa= '\0';
/*	Char the_aa_s[2];*/
	
	whichAAlist = GetValue(modifyRlist);	
 	
	if (isAA(modAA[whichAAlist-1][0])==0)
	{
		the_aa = 'X';
		modresParent=modAA[whichAAlist-1][1];
	}
	
	else the_aa = modAA[whichAAlist-1][0];
	
/*	the_aa_s[0] = the_aa;
	the_aa_s[1] = '\0';*/

 	if (reseditinfo->AA=='C') Vtrj_ChangeFromC();
	else if (reseditinfo->AA=='P') Vtrj_ChangeFromP();
 	else if (the_aa=='C') Vtrj_ChangeToC();
 	else if (the_aa=='P') Vtrj_ChangeToP();
 	
 	reseditinfo->AA = the_aa;
 	resinfo->AA=reseditinfo->AA;
	
 	if (ResReplace==FALSE) dowhat = EXTAASEQ_INSERT;
 	else dowhat = EXTAASEQ_REPLACE;
 	
 	newsequence = AlterResidueInSequence(sequence, (int) res_num, modAA[whichAAlist-1],  dowhat);
 	
 	
 	if (newsequence == NULL)
 		ErrPostEx(SEV_ERROR,1,1,"Error in creating new amino acid sequence");
	
	AAmodified=TRUE;

 	if (ResReplace == FALSE)
 	{
 		Vtrj_InsertResidue(the_aa);
 		DoneModifyRProc(i);
		return;
 	}
 	
	StringCpy(sequence, newsequence);
	newsequence = MemFree(newsequence);
	Vtrj_LoadSequence();
	Vtrj_SaveResidue();
	SwitchResProc(FALSE);
	Hide(aaname_panel);
	Show(aaname_panel);
	AAmodified=FALSE;
	lastg=0;
	needApply=FALSE;
	Disable(applyR_bttn);

	DoneModifyRProc(i);
}

static Int2 GetAAPos(void)
{

	if(res_num==1) return MODAA_N;
	if(res_num > 1 && res_num < numAA) return MODAA_M;
	if(res_num==numAA) return MODAA_C;

	return 0;
}

static void LoadAAList(Int2 WhereModAA)
{
	int i=0, isok;
	Char tempaa = 'A';
	Char aa[2];
	ValNodePtr tmp_modAA=NULL;
	Boolean doagain=TRUE, shouldcontinue=TRUE;
	PEAS peashere = NULL;
	
	tmp_Glob=0;

	while (tempaa != 'Z')
	{
		 isok = isAA(tempaa);
		 if (isok==1)
		 {
		 	aa[0] = tempaa;
		 	aa[1] = '\0';
		 	modAA[i][0] = tempaa;
		 	modAA[i][1] = '\0';
			tmp_Glob=i+1;
			LoadInfoTextProc(tmp_Glob, 2);
		 	if (isNew==FALSE) ListItem (modifyRlist, aa);
			if (isNew==TRUE) ListItem(newAAlist, aa);
		 	++i;
		 }
		 ++tempaa;
	}

	tmp_Glob=0;

	if(WhereModAA==MODAA_C) tmp_modAA = CmodAAlist;
	if(WhereModAA==MODAA_M) tmp_modAA = MmodAAlist;
	if(WhereModAA==MODAA_N) tmp_modAA = NmodAAlist;


	while (doagain==TRUE)
	{

		while (tmp_modAA !=NULL)	
		{
			
			peashere = (PEAS) tmp_modAA->data.ptrvalue;

			StringCpy(modAA[i], peashere->keybname);
			StringCpy(descrAA[i+1], peashere->descr);

			if(isNew==FALSE) 
			{
				ListItem (modifyRlist, modAA[i]);

			}
			if(isNew==TRUE)
			{
				ListItem(newAAlist, modAA[i]);

			}

			tmp_modAA = tmp_modAA->next;
			++i;
		}
		
		if (shouldcontinue==FALSE) break;

		if(WhereModAA==MODAA_M)
		{
			if(isNew==FALSE) doagain=FALSE;
			if(isNew==TRUE)
			{
				tmp_modAA=CmodAAlist;
				shouldcontinue = FALSE;
			}
		}
		
		if(WhereModAA==MODAA_C)
		{
			tmp_modAA=MmodAAlist;
			shouldcontinue = FALSE;
		}

		if(WhereModAA==MODAA_N)
		{
			tmp_modAA=MmodAAlist;
			shouldcontinue = FALSE;
		}
	
	}

}

static void EnabOkRProc (IteM i)
{
	Enable(okR_bttn);
}

static void LoadModifyRProc(IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2;
/*	ButtoN cancelAA_bttn;*/
	PrompT Mutating_prompt;
	Char printamino[5];
/*	Int2 Whereseq=0;*/
	Int2 ans;
	
	if (needApply==TRUE) {
		ans=Message (MSG_YN, "Apply unsaved changes before mutating?");
		if (ans==ANS_YES) {
			ApplyREditProc (i);
			/* check if apply failed */
			if (needApply==TRUE)
				return;		
		}
	}
	
	isAAused=TRUE;
	/* reload info for editing */
	reseditinfo->pCis=resinfo->pCis;
	reseditinfo->pSS=resinfo->pSS;
	reseditinfo->ChiWMean=resinfo->ChiWMean;
	reseditinfo->ChiWSD=resinfo->ChiWSD;
	reseditinfo->rotid=resinfo->rotid;
	reseditinfo->tout=resinfo->tout;
	reseditinfo->markovsf=resinfo->markovsf;
	lastg=0;
	needApply=FALSE;
	Disable(applyR_bttn);

	modifyR_win = MovableModalWindow(-50,-33,-10,-10, "Amino Acid Selection", (WndActnProc) DoneModifyRProc);

	g = HiddenGroup(modifyR_win,-1,-1, (GrpActnProc) NULL);
		
	SetGroupMargins(g, 17, 3);
	SetGroupSpacing(g, 3, 3);	
	
	Vtrj_grp2 = HiddenGroup(g, 2, 0, (GrpActnProc) NULL);
	StaticPrompt (Vtrj_grp2, "Current Amino Acid:", 0, dialogTextHeight, systemFont, 'l');
	Mutating_prompt = StaticPrompt (Vtrj_grp2, " ", 10, dialogTextHeight, programFont, 'l');
	sprintf(printamino, "%c",   reseditinfo->AA);
  SetTitle(Mutating_prompt, printamino);
	
	modifyRlist = SingleList (g,9,9, (LstActnProc) EnabOkRProc);


	LoadAAList(GetAAPos());	
	
	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);
	
	SetGroupSpacing(Vtrj_grp1, 15, 3);

	if (ResReplace==TRUE)
		okR_bttn = PushButton(Vtrj_grp1, "Mutate", (BtnActnProc) OkmodRProc);
	else
		okR_bttn = PushButton(Vtrj_grp1, "Insert", (BtnActnProc) OkmodRProc);
	Disable(okR_bttn);
/*	cancelAA_bttn =*/ DefaultButton(Vtrj_grp1, "Cancel", (BtnActnProc) DoneModifyRProc);
	
	AlignObjects(ALIGN_CENTER, Vtrj_grp1, g, NULL);
	
	Show(modifyR_win);
	ProcessEvents();

}

/*****************************************************************
                                       RESIDUE PROPERTIES EDIT FUNCTIONS
*****************************************************************/

static void LoadREditProc (IteM i)
{	
	Boolean  isenab = FALSE;
  Char printwhat[100], printamino[5];
  PRS prsthis;
	Int2 numc;
	Char resname[2];
	
/*	Disable (applyR_bttn);*/
        
	printwhat[0] = '\0';

/*	SetTitle(currentRvalue_text, "");*/
	
	prsthis=(PRS)MemNew(sizeof(RS));
	GetChiFromRotid(&prsthis,reseditinfo->rotid);	
	resname[0]=(Psequence[res_num-1]);
	resname[1]='\0';
	numc=numchi[StringCSpn(aalist,resname)];
	
	editwhatR = GetValue(editRlist);
				
	switch (editwhatR)
	{
		case 1:	
			sprintf(printwhat, "%d", reseditinfo->Peak);
			break;
		case 2:	
			sprintf(printwhat, "%d", reseditinfo->TrajIntegral);
			break;
		case 3:	
			sprintf(printwhat, "%d", reseditinfo->CisPeak);
			break;
		case 4:	
			sprintf(printwhat, "%d", reseditinfo->CisTrajIntegral);
			break;
		case 5:	
			sprintf(printwhat, "%1.2f", reseditinfo->pCis);
			isenab = TRUE;
			break;
		case 6:	
			sprintf(printwhat, "%1.2f", reseditinfo->pSS);
			isenab = TRUE;
			break;
		case 7:	
			sprintf(printwhat, "%1.1f", reseditinfo->ChiWMean);
			isenab = TRUE;
			break;
		case 8:	
			sprintf(printwhat, "%1.1f", reseditinfo->ChiWSD);
			isenab = TRUE;
			break;
		case 9:	
			if (numc>0) {
				if (prsthis->Chi1==0.0)
					sprintf(printwhat, "Rotamer Library");
				else
					sprintf(printwhat, "%1.1f", prsthis->Chi1);
				isenab = TRUE;
			}
			else
				isenab = FALSE;	
			break;
		case 10:	
			if (numc>1) {
				if (prsthis->Chi2==0.0)
					sprintf(printwhat, "Rotamer Library");
				else
					sprintf(printwhat, "%1.1f", prsthis->Chi2);
				isenab = TRUE;
			}
			else
				isenab = FALSE;	
			break;
		case 11:	
			if (numc>2) {
				if (prsthis->Chi3==0.0)
					sprintf(printwhat, "Rotamer Library");
				else
					sprintf(printwhat, "%1.1f", prsthis->Chi3);
				isenab = TRUE;
			}
			else
				isenab = FALSE;	
			break;
		case 12:	
			if (numc>3) {
				if (prsthis->Chi4==0.0)
					sprintf(printwhat, "Rotamer Library");
				else
					sprintf(printwhat, "%1.1f", prsthis->Chi4);
				isenab = TRUE;
			}
			else
				isenab = FALSE;	
			break;				
		case 13:	
			sprintf(printwhat, "%d", reseditinfo->dim);
			break;
		case 14:	
			sprintf(printwhat, "%d", reseditinfo->firstnzrow);
			break;
		case 15:	
			sprintf(printwhat, "%d", reseditinfo->numnzrows);
			break;
		case 16:	
			sprintf(printwhat, "%d", reseditinfo->nelt0);
			break;
		case 17:	
			sprintf(printwhat, "%d", reseditinfo->nelt5);
			break;
		case 18:	
			sprintf(printwhat, "%d", reseditinfo->nelt10);
			break;
		case 19:	
			sprintf(printwhat, "%d", reseditinfo->nelt15);
			break;
		case 20:	
			sprintf(printwhat, "%d", reseditinfo->tout);
			isenab = TRUE;
			break;
		case 21:	
			sprintf(printwhat, "%1.2f", reseditinfo->markovsf);
			isenab = TRUE;
			break;
	
	}
	
	prsthis=MemFree(prsthis);
	
	sprintf(printamino, "%c", reseditinfo->AA);

	SetTitle(AminoEdit_title, printamino);
	
	/*if (ismod == TRUE)  Show(modifyR_bttn);
	else if (ismod == FALSE)  Hide(modifyR_bttn);*/
	if (lastg!=editwhatR) {
		SetTitle(currentRvalue_text, printwhat);
		lastg=editwhatR;
	}
	
	if (isenab == TRUE) Enable(currentRvalue_text);
	else if (isenab == FALSE)  Disable(currentRvalue_text);
	
	LoadInfoTextProc(editwhatR, 4);
	
}

static void EnableRSaveEditProc (IteM i)
{
	Int2 editvalue=0;
	Char str[30];
	int temp_tout;
	float ftemp;
	PRS prseditinfo;
	
	if (AAmodified==FALSE)
	{		
		editvalue = GetValue(editRlist);
    GetTitle(currentRvalue_text, str, sizeof(str));

  	prseditinfo=(PRS)MemNew(sizeof(RS));
  	GetChiFromRotid(&prseditinfo,reseditinfo->rotid);
		switch (editvalue)
		{
			
			case 5:
				sscanf(str, "%f", &ftemp);
				reseditinfo->pCis=(FloatLo)ftemp;
				break;
			case 6:
				sscanf(str, "%f", &ftemp);
				reseditinfo->pSS=(FloatLo)ftemp;
				break;
			case 7:
				sscanf(str, "%f", &ftemp);
				reseditinfo->ChiWMean=(FloatLo)ftemp;
				break;			
			case 8:
				sscanf(str, "%f", &ftemp);
				reseditinfo->ChiWSD=(FloatLo)ftemp;
				break;
			case 9:
				sscanf(str, "%f", &ftemp);
				prseditinfo->Chi1=(FloatHi)ftemp;
				reseditinfo->rotid=ComputeRotid(prseditinfo);
				break;
			case 10:
				sscanf(str, "%f", &ftemp);
				prseditinfo->Chi2=(FloatHi)ftemp;
				reseditinfo->rotid=ComputeRotid(prseditinfo);
				break;
			case 11:
				sscanf(str, "%f", &ftemp);
				prseditinfo->Chi3=(FloatHi)ftemp;
				reseditinfo->rotid=ComputeRotid(prseditinfo);
				break;
			case 12:
				sscanf(str, "%f", &ftemp);
				prseditinfo->Chi4=(FloatHi)ftemp;
				reseditinfo->rotid=ComputeRotid(prseditinfo);
				break;											
			case 20:
				sscanf(str, "%d", &temp_tout);
				reseditinfo->tout=(Int2)temp_tout;
				break;
			case 21:
				sscanf(str, "%f", &ftemp);
				reseditinfo->markovsf=(FloatLo)ftemp;
				break;
		}
		reseditinfo->rotid=ComputeRotid(prseditinfo);
		prseditinfo=MemFree(prseditinfo);
	}
	
	LoadREditProc(i);
	Enable(applyR_bttn);
	needApply=TRUE;
}

static void ApplyREditProc (IteM i)
{
			
	int temp_tout;
	float ftemp;
	PRS prsinfo,prseditinfo;
		
	Rhasbeenedited = TRUE;
  		
	if (AAmodified==FALSE)
	{		
  	prsinfo=(PRS)MemNew(sizeof(RS));
  	prseditinfo=(PRS)MemNew(sizeof(RS));
  	GetChiFromRotid(&prsinfo,resinfo->rotid);
  	GetChiFromRotid(&prseditinfo,reseditinfo->rotid);
		
		ftemp=reseditinfo->pCis;
		if (ftemp<0.0 || ftemp>1.0) {			
			Message(MSG_OK,"p(cis) must be between zero and one");
			SetValue(editRlist,5);
			LoadREditProc(i);
			return;
		}
		ftemp=reseditinfo->pSS;
		if (ftemp<0.0 || ftemp>1.0) {
			Message(MSG_OK,"p(SS) must be between zero and one");
			SetValue(editRlist,6);
			LoadREditProc(i);
			return;
		}
		ftemp=reseditinfo->ChiWMean;
		if (ftemp<=-180.0 || ftemp>180.0) {
			Message(MSG_OK,"Please use a value greater than -180 and less than or equal to 180 degrees for Omega (Mean)");
			SetValue(editRlist,7);
			LoadREditProc(i);
			return;
		}
		ftemp=reseditinfo->ChiWSD;
		if (ftemp<0.0 || ftemp>45.0) {
			Message(MSG_OK,"Please use a value between 0 and 45 degrees for Omega (Standard Deviation)");
			SetValue(editRlist,8);
			LoadREditProc(i);
			return;
		}
		ftemp=prseditinfo->Chi1;
		if (ftemp<-180.0 || ftemp>180.0) {
			Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi1");
			SetValue(editRlist,9);
			LoadREditProc(i);
			return;
		}			
		ftemp=prseditinfo->Chi2;
		if (ftemp<-180.0 || ftemp>180.0) {
			Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi2");
			SetValue(editRlist,10);
			LoadREditProc(i);
			return;
		}
		else if ((Psequence[res_num-1]=='D' || Psequence[res_num-1]=='F' || Psequence[res_num-1]=='Y') && (ftemp<-90.0 || ftemp>90.0)) {
			Message(MSG_OK,"Please use a value between -90 and 90 degrees for chi2 of D, F or Y");
			SetValue(editRlist,10);
			LoadREditProc(i);
			return;
		}
		ftemp=prseditinfo->Chi3;
		if (ftemp<-180.0 || ftemp>180.0) {
			Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi3");
			SetValue(editRlist,11);
			LoadREditProc(i);
			return;
		}
		ftemp=prseditinfo->Chi4;
		if (ftemp<-180.0 || ftemp>180.0) {
			Message(MSG_OK,"Please use a value between -180 and 180 degrees for chi4");
			SetValue(editRlist,12);
			LoadREditProc(i);
			return;
		}
		temp_tout=reseditinfo->tout;
		if (temp_tout<BACKTRACK_TRIES_MIN || temp_tout>BACKTRACK_TRIES_MAX) {
			Message(MSG_OK,"Tries per Residue must be between %d and %d inclusive",BACKTRACK_TRIES_MIN,BACKTRACK_TRIES_MAX);
			SetValue(editRlist,20);
			LoadREditProc(i);
			return;
		}
		ftemp=reseditinfo->markovsf;
		if (ftemp<0.0 || ftemp>1.0) {
			Message(MSG_OK,"Markov Scale Factor must be between zero and one");
			SetValue(editRlist,21);
			LoadREditProc(i);
			return;
		}

		resinfo->pCis = reseditinfo->pCis;
		resinfo->pSS = reseditinfo->pSS;
		resinfo->ChiWMean = reseditinfo->ChiWMean;
		resinfo->ChiWSD = reseditinfo->ChiWSD;
		prsinfo->Chi1=prseditinfo->Chi1;		
		prsinfo->Chi2=prseditinfo->Chi2;		
		prsinfo->Chi3=prseditinfo->Chi3;		
		prsinfo->Chi4=prseditinfo->Chi4;		
		resinfo->rotid = reseditinfo->rotid;
		resinfo->tout = reseditinfo->tout;
		resinfo->markovsf = reseditinfo->markovsf;				
		resinfo->rotid=ComputeRotid(prsinfo);
		reseditinfo->rotid=ComputeRotid(prseditinfo);
		prsinfo=MemFree(prsinfo);
		prseditinfo=MemFree(prseditinfo);
		
		needApply=FALSE;
		Disable(applyR_bttn);
	}
	
	LoadREditProc(i);
	Vtrj_SaveResidue();
	
}

static void DoneEditRInfoProc (IteM i)
{
	Int2 ans;
	
	if (needApply==TRUE) {
		ans=Message (MSG_YN, "Apply unsaved changes before closing?");
		if (ans==ANS_YES) {
			ApplyREditProc (i);
			/* check if apply failed */
			if (needApply==TRUE)
				return;		
		}
	}
	lastg=0;
	needApply=FALSE;
	Hide(modifyR_win);
	Hide(editRinfo_win);
}



static void EditRInfoProc (IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3, Vtrj_grp4;
/*	ButtoN cancelR_bttn;
	PrompT AAEdit_title;*/
	
	reseditinfo->pCis=resinfo->pCis;
	reseditinfo->pSS=resinfo->pSS;
	reseditinfo->ChiWMean=resinfo->ChiWMean;
	reseditinfo->ChiWSD=resinfo->ChiWSD;
	reseditinfo->rotid=resinfo->rotid;
	reseditinfo->tout=resinfo->tout;
	reseditinfo->markovsf=resinfo->markovsf;
	
	needApply=FALSE;
	
	editRinfo_win = MovableModalWindow(-50,-33,460,-10, "Edit Residue Properties",  (WndActnProc) DoneEditRInfoProc);
	
	g = HiddenGroup(editRinfo_win,2,0, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,0,2,(GrpActnProc) NULL);

	Vtrj_grp2 = HiddenGroup(Vtrj_grp1, -1, -1, (GrpActnProc) NULL);
	
	SetGroupMargins(g, 10, 10);
	SetGroupSpacing(g, 5, 15);			
	
	editRlist = SingleList (Vtrj_grp2,16,8, (LstActnProc) LoadREditProc);

	/*ListItem (editRlist,"AA");	*/
	ListItem (editRlist, "Peak");
	ListItem (editRlist, "Integral");
	ListItem (editRlist,"Cis-Peak");
	ListItem (editRlist,"Cis-Integral");
	ListItem (editRlist, "p(cis)");
	ListItem (editRlist,"p(SS)");
	ListItem (editRlist,"Omega (Mean)");
	ListItem (editRlist,"Omega (SD)");
	ListItem (editRlist,"Chi1");
	ListItem (editRlist,"Chi2");
	ListItem (editRlist,"Chi3");
	ListItem (editRlist,"Chi4");	
	ListItem (editRlist,"Discretization Size");
	ListItem (editRlist,"First Row");
	ListItem (editRlist,"Number of Rows");
	ListItem (editRlist,"# Elements <= 0% peak");
	ListItem (editRlist,"# Elements <= 5% peak");
	ListItem (editRlist,"# Elements <= 10% peak");	
	ListItem (editRlist,"# Elements <= 15% peak");	
	ListItem (editRlist,"Tries per Residue");
	ListItem (editRlist,"Markov Scale Factor");
		
	
	currentRvalue_text = DialogText (Vtrj_grp1, "",  5, (TxtActnProc) EnableRSaveEditProc);
	Disable(currentRvalue_text);
	
	Select(currentRvalue_text);

	
	Vtrj_grp3 = HiddenGroup(g, 0, 3, (GrpActnProc) NULL);
	
/*	AAEdit_title =*/ StaticPrompt (Vtrj_grp3, "Amino Acid:", 0, dialogTextHeight, systemFont, 'c');

	AminoEdit_title = StaticPrompt (Vtrj_grp3, " ", 0, dialogTextHeight, programFont, 'c');

	modifyR_bttn = PushButton(Vtrj_grp3, "Mutate Residue", (BtnActnProc) LoadModifyRProc);
		
	Vtrj_grp4 = HiddenGroup(g, 3, 0, (GrpActnProc) NULL);

	SetGroupSpacing(Vtrj_grp4, 10, 0);

	applyR_bttn = PushButton(Vtrj_grp4, "Apply", (BtnActnProc) ApplyREditProc);
	Disable(applyR_bttn);

/*	cancelR_bttn =*/ DefaultButton(Vtrj_grp4, "Close", (BtnActnProc) DoneEditRInfoProc);


	AlignObjects(ALIGN_CENTER, g, Vtrj_grp4, NULL);

	Break(g);
	
	REdit_Info = ScrollText (g, 24, 8, systemFont, TRUE, NULL);
	
	/*Sets which function is called when a key is pressed*/
	SetSlateChar ((SlatE) editRinfo_win, Vtrj_ModualKeyPressed);

	
	LoadREditProc(i);
	
	Show(editRinfo_win);
	
	ProcessEvents();

}



/*****************************************************************
			NOISE FUNCTIONS
*****************************************************************/

static void ClearGorText (void)
{
        	SetTitle(gor_e_text, "");
        	SetTitle(gor_h_text,"");
        	SetTitle(gor_c_text,"");

}

static Boolean CheckGor(void)
{
	int h=0, e=0, c=0,err;
    Char strh[5];
    Char stre[5];
    Char strc[5];
		
    GetTitle(gor_h_text, strh, sizeof(strh));
    GetTitle(gor_e_text, stre, sizeof(stre));
    GetTitle(gor_c_text, strc, sizeof(strc));
              	
	 err=sscanf(strh, "%d", &h);
	 if (err==0)
	 	h=0;
	 err=sscanf(stre, "%d", &e);
	 if (err==0)
	 	e=0;
	 err=sscanf(strc, "%d", &c);
	 if (err==0)
	 	c=0;
	
	 if (h+e+c!=100)
	 {
		Message (MSG_OK, "You must enter percentages that add up to 100.");
		ClearGorText();
		return FALSE;
	 }
	 
	 if (h<0 || e<0 || c<0)
	 {
		Message (MSG_OK, "You must enter positive percentages!");
		ClearGorText();
		return FALSE;
	 }
		
	Vtrj_prob[0] = (FloatLo) h/100.0;
	Vtrj_prob[1] = (FloatLo) e/100.0;
	Vtrj_prob[2] = (FloatLo) c/100.0;
		
	return TRUE;

}


static Boolean CheckSstru(void)
{
	Int2 whichstru=0;

	whichstru = GetValue(sstru_list);
	
	
  	if (whichstru==3)
  	{
  		Vtrj_prob[0] = 0.00;
  		Vtrj_prob[1] = 0.00;
  		Vtrj_prob[2] = 100.0;
		Vtrj_sstru = 'C';
  	}
  	
  	if (whichstru==2)
  	{
  		Vtrj_prob[0] = 0.00;
  		Vtrj_prob[1] = 100.0;
  		Vtrj_prob[2] = 0.00;
		Vtrj_sstru = 'E';
  	}
  	
  	if (whichstru==1)
  	{
  		Vtrj_prob[0] = 100.0;
  		Vtrj_prob[1] = 0.00;
  		Vtrj_prob[2] = 0.00;
		Vtrj_sstru = 'H';
  	}
  	
	return TRUE;
}

static Boolean CheckSP(void)
{
	float x=0.0, y=0.0, z=0.0;
	Char strx[10];
	Char stry[10];
	Char strz[10];
	int err;

	GetTitle(gor_h_text, strx, sizeof(strx));
	GetTitle(gor_e_text, stry, sizeof(stry));
	GetTitle(gor_c_text, strz, sizeof(strz));

	err=sscanf(strx, "%f", &x);
	if (err==0)
		x=0;
	err=sscanf(stry, "%f", &y);
	if (err==0)
		y=0;
	err=sscanf(strz, "%f", &z);
	if (err==0)
	 	z=0;
	
	if (z<=0)
	{
		Message(MSG_OK, "Magnitude must be one or greater!");
		SetTitle(gor_c_text, "");
		return FALSE;
	}
	if (WALKTYPE==WALK_CA)
	{ 		
	 	if (x<-180.0 || x >=180.0)
	 	{
	 		Message(MSG_OK, "Your Theta Value is out of bounds. \n Please enter a number between 0 and 360.");
	 		SetTitle(gor_h_text, "");
	 		return FALSE;
	 	}
	 	if (y<-1.0 || y>=1.0)
	 	{
	 		Message(MSG_OK, "Your Cos(Phi) Value is out of bounds. \n Please enter a number between -1 and 1.");
	 		SetTitle(gor_e_text, "");
	 		return FALSE;
	 	}
	 	
	 	x += 180.0;
	 	x /=  360.0;
	 	x *= (float) resinfo->dim;
	 	
	 	y += 1.0;
	 	y /= 2.0;
	 	y *= (float) resinfo->dim;
	}
	
	if (WALKTYPE==WALK_PHIPSI)
	{ 	
		if (x<-180.0 || x>=180.0)
		{
			Message(MSG_OK, "Your Phi Value is out of bounds. \n Please enter a number between -180 and 180.");
			SetTitle(gor_h_text, "");
			return FALSE;
		}
		if (y<-180.0 || y>=180.0)
		{
			Message(MSG_OK, "Your Psi Value is out of bounds. \n Please enter a number between -180 and 180.");
			SetTitle(gor_e_text, "");
			return FALSE;
		}
		
		x += 180.0;
		x /= 360.0;
		x *= (float) resinfo->dim;
		
		y += 180.0;
		y /= 360.0;
		y *= (float) resinfo->dim;
		
	 }
	 		
	
	Vtrj_SPposition = (int) resinfo->dim*(int)y+(int)x;
	Vtrj_SPpeak = (FloatLo) z;		
	 			
		
	return TRUE;


}

static void EnableGorOK (IteM i)
{
	Enable(noise_ok_bttn);
}

static void DoneNoiseProc (IteM i)
{
  	Hide(addnoise_win);
}

static void DoGor(void)
{
	Disable(noise_ok_bttn);
	
    SetTitle(gor_h_title, "Helix [%]");
    SetTitle(gor_e_title, "Strand [%]");
    SetTitle(gor_c_title, "Coil [%]");
    SetTitle(prop_title, "Composition");
    ClearGorText();
}

static void DoSstru(void)
{

	Enable(noise_ok_bttn);
    
	SetTitle(gor_h_title, "");
    SetTitle(gor_e_title, "");
    SetTitle(gor_c_title, "");
    SetTitle(prop_title, "Secondary Structure");
    
	ClearGorText();
	
	Hide(gor_h_text);
	Hide(gor_e_text);
	Hide(gor_c_text);
	Hide(gor_h_title);
	Hide(gor_e_title);
	Hide(gor_c_title);

	Show(sstru_list);

}


static void DoSP(void)
{

	Disable(noise_ok_bttn);
	
	if (WALKTYPE==WALK_CA)
	{
		SetTitle(gor_h_title, "Theta [Degrees]");
		SetTitle(gor_e_title, "Cos(Phi)");
	}
	if (WALKTYPE==WALK_PHIPSI)
	{
	    SetTitle(gor_h_title, "Phi [Degrees]");
	    SetTitle(gor_e_title, "Psi [Degrees]");
	}
		        	
    SetTitle(gor_c_title, "Peak Magnitude");
    SetTitle(prop_title, "Single Peak Placement");
    ClearGorText();

}


static void InVisGor (IteM i)
{
	Hide(gor_h_text);
	Hide(gor_e_text);
	Hide(gor_c_text);
	Hide(gor_h_title);
	Hide(gor_e_title);
	Hide(gor_c_title);
	Hide(prop_title);
	Hide(sstru_list);
}


static void VisGor (IteM i)
{
	Show(gor_h_text);
  Select(gor_h_text);
	Show(gor_e_text);
	Show(gor_c_text);
	Show(gor_h_title);
	Show(gor_e_title);
	Show(gor_c_title);
	Show(prop_title);

	Hide(sstru_list);
}


static void MakeNoiseProc (IteM i)
{

	Boolean isok = FALSE;
	Int2 whichnoise = 0, whichadd = 0, trajtype = 0;
	Char cres[2];
	Int2 Cis=0;
	
	cres[0] = resinfo->AA;
	
	if (cres[0] == 'X')
	{
		GetModResInSequence();
		cres[0] = modresParent;
	}	
		
	cres[1] = '\0';

	if (trajnoise != NULL) FreeTraj(trajnoise);
	if (resinfo->CisTrajGraph!=NULL) Cis = 1;
	
	trajnoise = NewTraj(Cis, resinfo->dim);
		
	whichnoise = GetValue(noise_type);
	whichadd = GetValue(add_or_rep);
	
	switch (whichnoise)
	{
		case 1:
			isok = TRUE;
			trajtype = TRAJ_UNIFORM;
			break;
		case 2:
			isok = TRUE;
			trajtype = TRAJ_STANDARD;
			break;			
		case 3:
			trajtype = TRAJ_SSTRU;
			isok = CheckSstru();
			break;
		case 4:
			trajtype = TRAJ_GOR;
			isok = CheckGor();
			break;
		case 5:
			trajtype = -1;
			isok = CheckSP();
			break;
	}
	
	if (isok==FALSE)
		return;
	trajnoise->dim = resinfo->dim;
	
	switch (whichadd)
	{
		case 1:
			if (trajtype >=0)
			{
				if (isCis==FALSE) trajnoise->TrajIntegral = FillTG(trajtype, cres, Vtrj_sstru, Vtrj_prob, resinfo->dim, trajnoise->TrajGraph, FALSE);
				if (isCis==TRUE) trajnoise->CisTrajIntegral = FillTG(trajtype, cres, Vtrj_sstru, Vtrj_prob, resinfo->dim, trajnoise->CisTrajGraph, TRUE);
			}
			else if (trajtype <0)
			{
				if (isCis==FALSE)
				{
					trajnoise->TrajGraph[Vtrj_SPposition] += Vtrj_SPpeak;
					trajnoise->TrajIntegral += Vtrj_SPpeak;
				}
				
				if (isCis==TRUE)
				{
					trajnoise->CisTrajGraph[Vtrj_SPposition] += Vtrj_SPpeak;
					trajnoise->CisTrajIntegral += Vtrj_SPpeak;
				}
			}
			
			TrajAdd(resinfo, trajnoise);
			break;
			
		case 2:
			if (trajtype >=0)
			{
				if (isCis==FALSE) resinfo->TrajIntegral = FillTG(trajtype, cres, Vtrj_sstru, Vtrj_prob, resinfo->dim, resinfo->TrajGraph, FALSE);
				if (isCis==TRUE) resinfo->CisTrajIntegral = FillTG(trajtype, cres, Vtrj_sstru, Vtrj_prob, resinfo->dim, resinfo->CisTrajGraph, TRUE);
			}
			else if (trajtype <0)
			{
				if (isCis==FALSE)
				{
					resinfo->TrajGraph = (Int4 *) MemGet(sizeof(Int4)*resinfo->dim*resinfo->dim, MGET_CLEAR);	
				 	resinfo->TrajGraph[Vtrj_SPposition] = Vtrj_SPpeak;
					resinfo->TrajIntegral = Vtrj_SPpeak;
				 	
				} 	
			  	if (isCis==TRUE)
			  	{
					resinfo->CisTrajGraph = (Int4 *) MemGet(sizeof(Int4)*resinfo->dim*resinfo->dim, MGET_CLEAR);	
				 	resinfo->CisTrajGraph[Vtrj_SPposition] = Vtrj_SPpeak;
					resinfo->CisTrajIntegral = Vtrj_SPpeak;
			  	}
			}
			break;
	}	

	Vtrj_SaveResidue();

	Nhasbeenedited = TRUE;
	
	Vtrj_sstru = '\0';
	Vtrj_prob[0] = 0.0;
	Vtrj_prob[1] = 0.0;
	Vtrj_prob[2] = 0.0;
	Vtrj_SPpeak = 0.0;
	Vtrj_SPposition = 0;
	
	DoneNoiseProc(i);	
	
}

static void CheckNoiseProc (IteM  i)
{
	Int2 whichnoise = 0 ;
	
	whichnoise = GetValue(noise_type);
	
	switch (whichnoise)
	{
		case 1:
			EnableGorOK(i);
			InVisGor(i);
			break;
		case 2:
			EnableGorOK(i);
			InVisGor(i);
			break;
		case 3:
			VisGor(i);
			DoSstru();
			break;
	    case 4:
	       	VisGor(i);
	       	DoGor();
	       	break;
		case 5:
			VisGor(i);
			DoSP();
			break;	                	
	}
		

}


static void AddNoiseProc (IteM i)
{

	GrouP g, Vtrj_grp1, Vtrj_grp2, Vtrj_grp3, Vtrj_grp4;
/*	ButtoN NoiseDone_bttn;*/
		
	addnoise_win = MovableModalWindow(-50,-33,-10,-10, "Add/Replace Noise", (WndActnProc) DoneNoiseProc);

	g = HiddenGroup(addnoise_win,0,1, (GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 70, 3);			
	
	noise_type = HiddenGroup (g, 0, 5, (GrpActnProc) CheckNoiseProc);
	RadioButton (noise_type, "Uniform");
	RadioButton (noise_type, "Standard");
	RadioButton (noise_type, "One State Secondary Structure");
	RadioButton (noise_type, "Three State Secondary Structure");
	RadioButton (noise_type, "Single Peak");
	SetValue(noise_type, 2);
	
	Vtrj_grp1=  HiddenGroup(g,0,4,(GrpActnProc) NULL);
    StaticPrompt (Vtrj_grp1, "", 0, dialogTextHeight, systemFont, 'l');
    gor_h_title = StaticPrompt (Vtrj_grp1, "Theta [Degrees]", 0, dialogTextHeight, systemFont, 'l');
    gor_e_title = StaticPrompt (Vtrj_grp1, "Cos(Phi)", 0, dialogTextHeight, systemFont, 'l');
    gor_c_title = StaticPrompt  (Vtrj_grp1, "Peak Magnitude", 0, dialogTextHeight, systemFont, 'l');
    prop_title = StaticPrompt (Vtrj_grp1, "Single Peak Placement", 0, dialogTextHeight, systemFont, 'l');
	Vtrj_grp3 = HiddenGroup(Vtrj_grp1, 0,0,(GrpActnProc) NULL);
	gor_h_text = DialogText (Vtrj_grp3, "",  5, (TxtActnProc) EnableGorOK);
	gor_e_text = DialogText (Vtrj_grp1, "",  5, (TxtActnProc) EnableGorOK);
	gor_c_text = DialogText (Vtrj_grp1, "",  5, (TxtActnProc) EnableGorOK);
	sstru_list = PopupList(Vtrj_grp3, TRUE, (PupActnProc) EnableGorOK);

	AlignObjects (ALIGN_JUSTIFY, gor_h_text, gor_e_text, NULL);

	
	Vtrj_grp2 = HiddenGroup(g, 0, 2, (GrpActnProc) NULL);	
	add_or_rep = HiddenGroup(Vtrj_grp2, 0, 2, (GrpActnProc) NULL);
	RadioButton (add_or_rep, "Add");
	RadioButton (add_or_rep, "Replace");
	SetValue(add_or_rep, 1);

	Break(g);
	
	Vtrj_grp4 = HiddenGroup(g, 2, -1, (GrpActnProc) NULL);

	noise_ok_bttn = DefaultButton(Vtrj_grp4, "OK", (BtnActnProc) MakeNoiseProc);
	
/*	NoiseDone_bttn =*/ PushButton(Vtrj_grp4, "Cancel", (BtnActnProc) DoneNoiseProc);
	
	AlignObjects (ALIGN_CENTER, Vtrj_grp4, g, NULL);

	PopupItem(sstru_list, "Helix            ");
	PopupItem(sstru_list, "Sheet            ");
	PopupItem(sstru_list, "Coil             ");
	SetValue(sstru_list, 1);

	InVisGor(i);
	
	Show(addnoise_win);
	ProcessEvents();
	
}

/*****************************************************************
			RESIDUE ADD & REMOVE FUNCTIONS
*****************************************************************/

static void Vtrj_InsertResidue(Char the_aa)
{

    Char the_aa_s[2];
    Int4 numrecords = 0;


	if (the_aa=='X') the_aa=modresParent;	
	the_aa_s[0] = the_aa;
    the_aa_s[1] = '\0';


	StringCpy(sequence, newsequence);
	newsequence = MemFree(newsequence);
	
	TGInit(tmpdbasename, DB_READ, &numrecords);

	if(AdjustResidueNumbers(numAA, res_num, EXTAASEQ_INSERT)==ERR_FAIL)
 		ErrPostEx(SEV_ERROR,1,1,"Could not properly adjust residue numbers");
	
	TGClose();
 		
	resinfo->TrajGraph = MemFree (resinfo->TrajGraph);
  	resinfo->TrajGraph = (Int4 *) MemNew(sizeof(Int4)*resinfo->dim*resinfo->dim);	
	resinfo->TrajIntegral = FillTG(TRAJ_STANDARD, the_aa_s, 0, NULL, resinfo->dim, resinfo->TrajGraph, FALSE);
	if (resinfo->CisTrajGraph!=NULL)
		resinfo->CisTrajGraph = MemFree (resinfo->CisTrajGraph);
	if (the_aa=='P')
		Vtrj_ChangeToP();
	else
		Vtrj_ChangeFromP();

	Vtrj_SaveResidue();

	numAA+= 1;
		
	ResReplace = TRUE;
	
	SetSwitchParams (resid_switch, 1, numAA);
	SetValue(resid_switch, res_num);
	
	Vtrj_LoadSequence();
	SwitchResProc(FALSE);
	Hide(aaname_panel);
	Show(aaname_panel);
	
	Rhasbeenedited=TRUE;
}

static void AddResidueProc(IteM i)
{

	Char tmp_str[MAX_KBAA_NAMELENGTH+1];
	PEAS tmppeas=NULL;
	ValNodePtr tmpvalnode=NULL;

	tmp_str[0] = '\0';

	if(Xsequence[res_num-1]=='X')
	{
		GetTitle(modAA_title, tmp_str, sizeof(tmp_str));
		
		tmpvalnode = NmodAAlist;
		
		while(tmpvalnode!=NULL)
		{
			tmppeas = (PEAS) tmpvalnode->data.ptrvalue;
			if(!StringICmp(tmp_str,tmppeas->keybname))
			{
				Message (MSG_OK, "You may not insert an amino acid before this N-terminus modified amino acid.");
				return;
			}
			tmpvalnode = tmpvalnode->next;
		}
	}

	ResReplace = FALSE;
	LoadModifyRProc(i);
	ResReplace = TRUE;
}

static void RemoveResidueProc (IteM i)
{
	Int2 Vtrj_ans = ANS_NO;
	Int4 numrecords = 0;

	if(StringLen(Xsequence)<=3)
	{
		Message (MSG_OK, "You cannot have a sequence length of less than 3.");
		return;
	}

	
	Vtrj_ans = Message (MSG_YN, "Are you sure you would like to remove amino acid \"%c\" from the sequence?", Xsequence[res_num-1]);	
	
	
	if (Vtrj_ans==ANS_NO) return;
	
 	newsequence = AlterResidueInSequence(sequence, (int) res_num, "A",  EXTAASEQ_DELETE);
 	
 	if (newsequence == NULL)
 		ErrPostEx(SEV_ERROR,1,1,"Error in creating new amino acid sequence");
 	
 	StringCpy(sequence, newsequence);
	newsequence = MemFree(newsequence);

	TGInit(tmpdbasename, DB_READ, &numrecords);

	if(AdjustResidueNumbers(numAA, res_num, EXTAASEQ_DELETE)==ERR_FAIL)
 		ErrPostEx(SEV_ERROR,1,1,"Could not properly adjust residue numbers");
	
	TGClose();
	
	numAA -= 1;

	if(res_num>numAA) res_num = numAA;

	SetSwitchParams (resid_switch, 1, numAA);
	SetValue(resid_switch, res_num);
	
	
	Vtrj_LoadSequence();
	SwitchResProc(FALSE);  /*Not called to switch residues but rather to reload all
				variables with the new information produced
				by FillTG*/

	Rhasbeenedited=TRUE;					
}

/*****************************************************************
			RESIDUE POSITION FUNCTIONS
*****************************************************************/

/*static void Vtrj_DrawString (RectPtr rptr, CharPtr text, FonT fnt)

{
  if (fnt != NULL) SelectFont (fnt);


  rptr->bottom = rptr->top + LineHeight();

#ifdef WIN_MSWIN
	rptr->right = rptr->right + MaxCharWidth();
#endif

  DrawString (rptr, text, 'l', FALSE);
  rptr->left = rptr ->left + (MaxCharWidth() * 4);

#ifdef WIN_MOTIF  
  rptr->right = rptr->left + MaxCharWidth();
#endif

#ifdef WIN_MSWIN
	rptr->right = rptr->left;
#endif

  rptr->top = rptr->bottom - LineHeight();


}*/

static void DrawAANameProc (PaneL p)
{
	RecT r;
	static FonT  titleFont = NULL;
	Char pc[2];
	Int2 cnt,reshere;

#ifdef WIN_MOTIF
    titleFont = GetFont ("Courier", 14, FALSE, FALSE, FALSE, "");
#endif

#ifdef WIN_MSWIN
    titleFont = GetFont ("Courier", 14, TRUE, FALSE, FALSE, "");
#endif
	if (titleFont!=NULL)
	 	SelectFont (titleFont);
	ObjectRect (aaname_panel, &r);
	InsetRect (&r, 4, 2);
	pc[1]='\0';
	for (cnt=0;cnt<9;cnt++) {
		pc[0]=shortsequence[cnt];
		reshere=res_num-1+cnt-4;
		if (reshere>=0 && reshere<numAA && hasafrag!=NULL && hasafrag[reshere]==TRUE) {
			if (cnt==4)
				SelectColor(255,164,0);
			else
				SelectColor(0,164,0);
		}
		else if (cnt==4)
			SelectColor(255,0,0);
		else
			SelectColor(0,0,0);
	 	DrawString (&r, pc, 'l', FALSE);	  	
	  r.left+=MaxCharWidth();
	}
}


static void DisplayModResidue()
{

	PEAS tmppeas=NULL;
	ValNodePtr tmpvalnode=NULL;
	Int2 tmp_modresid=modresid;

	if(res_num==1)
	{
		tmpvalnode = NmodAAlist;
		
		while(tmpvalnode!=NULL)
		{
			tmppeas = (PEAS) tmpvalnode->data.ptrvalue;
			if(tmp_modresid==tmppeas->dictidx)
			{
				SetTitle (modAA_title, tmppeas->keybname);
				return;
			}
			tmpvalnode = tmpvalnode->next;
		}

		tmpvalnode = MmodAAlist;
		tmp_modresid-=2;

		while(tmpvalnode!=NULL)
		{
			tmppeas = (PEAS) tmpvalnode->data.ptrvalue;
			if(tmp_modresid==tmppeas->dictidx)
			{
				SetTitle (modAA_title, tmppeas->keybname);
				return;
			}
			tmpvalnode = tmpvalnode->next;
		}

	}


	else if(res_num==numAA)
	{
		tmpvalnode = CmodAAlist;
		
		while(tmpvalnode!=NULL)
		{
			tmppeas = (PEAS) tmpvalnode->data.ptrvalue;
			if(tmp_modresid==tmppeas->dictidx)
			{
				SetTitle (modAA_title, tmppeas->keybname);
				return;
			}
			tmpvalnode = tmpvalnode->next;
		}

		tmpvalnode = MmodAAlist;
		tmp_modresid--;

		while(tmpvalnode!=NULL)
		{
			tmppeas = (PEAS) tmpvalnode->data.ptrvalue;
			if(tmp_modresid==tmppeas->dictidx)
			{
				SetTitle (modAA_title, tmppeas->keybname);
				return;
			}
			tmpvalnode = tmpvalnode->next;
		}

	}

	else
	{

		tmpvalnode = MmodAAlist;

		while(tmpvalnode!=NULL)
		{
			tmppeas = (PEAS) tmpvalnode->data.ptrvalue;
			if(tmp_modresid==tmppeas->dictidx)
			{
				SetTitle (modAA_title, tmppeas->keybname);
				return;
			}
			tmpvalnode = tmpvalnode->next;
		}

	}

}


static void GetModResInSequence(void)
{
	Char modres[3];
	int nthX=0, tempX = 0;
	int i, modresidtemp;
	
	for (i=0; i<=res_num; ++i)
	{
		if (Xsequence[i] == 'X')
		{
			tempX += 1;
			if (i == (int) res_num-1) nthX = tempX;
		}
	}
	
	
		
	for (i=0; i <=((numAA)+(nthX*5)); ++i)
	{	
		if (sequence[i]=='*')
		{
			if ((int) res_num - 1 == i-(nthX-1)*5)
			{
				modres[0] = sequence[i+2];
				modres[1] = sequence[i+3];
				modres[2] = sequence[i+4];
				modresParent = sequence[i+5];
				sscanf(modres, "%d",  &modresidtemp);
				modresid = (Int2) modresidtemp;
				return;
			}
		}
	}
	
		
	
}

static void ChangeAATitle(void)
{
	

	if (res_num == 1)
	{
		shortsequence[0] = ' ';
		shortsequence[1] = ' ';
		shortsequence[2] = ' ';
		shortsequence[3] = ' ';	
	}		
	
	if (res_num == 2)
	{
		shortsequence[0] = ' ';
		shortsequence[1] = ' ';
		shortsequence[2] = ' ';
		shortsequence[3] = Xsequence[res_num-2];	
		
	}

	if (res_num == 3)
	{
		shortsequence[0] = ' ';
		shortsequence[1] = ' ';
		shortsequence[2] = Xsequence[res_num-3];
		shortsequence[3] = Xsequence[res_num-2];	
	}
	
	if (res_num == 4)
	{
		shortsequence[0] = ' ';
		shortsequence[1] = Xsequence[res_num-4];
		shortsequence[2] = Xsequence[res_num-3];
		shortsequence[3] = Xsequence[res_num-2];	
	}
	
	if (res_num > 4)
	{
		shortsequence[0] = Xsequence[res_num-5];
		shortsequence[1] = Xsequence[res_num-4];
		shortsequence[2] = Xsequence[res_num-3];
		shortsequence[3] = Xsequence[res_num-2];	
	}
	
			
	shortsequence[4] = Xsequence[res_num - 1];
	if (res_num<numAA)
		shortsequence[5] = Xsequence[res_num];
	else
		shortsequence[5] = ' ';
	if (res_num<numAA-1)
		shortsequence[6] = Xsequence[res_num+1];
	else
		shortsequence[6] = ' ';
	if (res_num<numAA-2)
		shortsequence[7] = Xsequence[res_num+2];
	else
		shortsequence[7] = ' ';
	if (res_num<numAA-3)
		shortsequence[8] = Xsequence[res_num+3];
	else
		shortsequence[8] = ' ';
	shortsequence[9] = '\0';
	
	pointedAA[0] = shortsequence[4];
	pointedAA[1] = '\0';	
	
	SetTitle (resnum_title, AAposition);
	
	DrawAANameProc(aaname_panel);	
    /*	if (resinfo!=NULL) Update();*/

	if (Xsequence[res_num-1] == 'X')
	{
		Show(modAA_temptitle);
		Show(modAA_title);
		GetModResInSequence();
		DisplayModResidue();
	}
	else
	{
		Hide(modAA_temptitle);
		Hide(modAA_title);
	}
}

static void SetPeaks(void)
{
	
	if (isCis==TRUE) ThePeak = (FloatLo) resinfo->CisPeak;
	else if (isCis==FALSE) ThePeak = (FloatLo) resinfo->Peak;

       	ScaleZ = ThePeak/SCALEZ;
        	ScaleSphere = ThePeak/SCALESPHERE;

}

static void Vtrj_StoreOrigResVariables(void)
{
	origresinfo->Peak = resinfo->Peak;
	origresinfo->CisTrajIntegral = resinfo->CisTrajIntegral;
	origresinfo->CisPeak = resinfo->CisPeak;
	origresinfo->pCis = resinfo->pCis;
	origresinfo->pSS = resinfo->pSS;
	origresinfo->ChiWMean = resinfo->ChiWMean;
	origresinfo->ChiWSD = resinfo->ChiWSD;
	origresinfo->AA = resinfo->AA;
	origresinfo->dim = resinfo->dim;
	origresinfo->firstnzrow = resinfo->firstnzrow;
	origresinfo->numnzrows = resinfo->numnzrows;
	origresinfo->nelt0 = resinfo->nelt0;
	origresinfo->nelt5 = resinfo->nelt5;
	origresinfo->nelt10 = resinfo->nelt10;
	origresinfo->nelt15 = resinfo->nelt15;
	origresinfo->tout = resinfo->tout;
	origresinfo->markovsf = resinfo->markovsf;
	origresinfo->rotid = resinfo->rotid;
	origresinfo->TrajIntegral = resinfo->TrajIntegral;
}

static void Vtrj_StoreResVariables(void)
{
	reseditinfo->Peak = resinfo->Peak;
	reseditinfo->CisTrajIntegral = resinfo->CisTrajIntegral;
	reseditinfo->CisPeak = resinfo->CisPeak;
	reseditinfo->pCis = resinfo->pCis;
	reseditinfo->pSS = resinfo->pSS;
	reseditinfo->ChiWMean = resinfo->ChiWMean;
	reseditinfo->ChiWSD = resinfo->ChiWSD;
	reseditinfo->AA = resinfo->AA;
	reseditinfo->dim = resinfo->dim;
	reseditinfo->firstnzrow = resinfo->firstnzrow;
	reseditinfo->numnzrows = resinfo->numnzrows;
	reseditinfo->nelt0 = resinfo->nelt0;
	reseditinfo->nelt5 = resinfo->nelt5;
	reseditinfo->nelt10 = resinfo->nelt10;
	reseditinfo->nelt15 = resinfo->nelt15;
	reseditinfo->tout = resinfo->tout;
	reseditinfo->markovsf = resinfo->markovsf;
	reseditinfo->rotid = resinfo->rotid;
	reseditinfo->TrajIntegral = resinfo->TrajIntegral;
}

static Int4 SwitchResProc (Boolean IsSwitching)
{
	  	
	Int2 Cis=0;
	Int4 numrecords;
	
	if(IsSwitching==TRUE)
	{
		isCis = FALSE;
		SetValue(which_trans,1);
	}
	  	
/*	if (JumpFromStart==TRUE)
	{
		res_num = (Int2) Vtrj_args[1].intvalue;
		
		if (res_num > numAA) res_num = numAA;
		if (res_num < 1) res_num = 1;
		
		SetValue(resid_switch, Vtrj_args[1].intvalue);
		JumpFromStart=FALSE;
	}
	else 
*/
	res_num = GetValue(resid_switch);
	
        	AAposition =  Ultostr(res_num,(Int2) NULL);
	
  ChangeAATitle();
	
	if (HasBeenFiltered==TRUE && IsSwitching==TRUE)
	{		
		Vtrj_SaveResidue();
		HasBeenFiltered = FALSE;
	}
	        	
  	if (resinfo != NULL) FreeTraj(resinfo);
	if (origresinfo != NULL) FreeTraj(origresinfo);
	if (reseditinfo != NULL) reseditinfo = MemFree(reseditinfo);
	

	TGInit(tmpdbasename, DB_READ, &numrecords);
	
        	resinfo =  TrajGraphRead((Int2) res_num);
	TGClose();        	
        	if (resinfo==NULL)
        	{
		ErrPostEx(SEV_ERROR,1,1,"No Traj Graph Loaded");
		exit(1);        		
	}
		
        	
        	if (resinfo->CisTrajGraph!=NULL) Cis = 1;
        	
        	origresinfo = NewTraj(Cis, resinfo->dim);
        	reseditinfo = (PTGS) MemNew (sizeof(TGS));
	
	TrajGraphDifferentiate(resinfo);
        	TrajCopy(resinfo->TrajGraph, origresinfo->TrajGraph, resinfo->dim);
                if (Cis==1) TrajCopy(resinfo->CisTrajGraph, origresinfo->CisTrajGraph, resinfo->dim);
                        	
        	Vtrj_StoreOrigResVariables();
        	Vtrj_StoreResVariables();
        	
	SetPeaks();
        	
	EnableDisableProc();
		 	
		
	RangeProc();

    Vtrj_LoadVertexes();

	if (IsSwitching==TRUE) {
		Vtrj_scale = 10.0;
		Vtrj_rotateX = 90;
		Vtrj_rotateY = -35;
		Vtrj_rotateZ = 0;
	}
       	
	isSwitching=TRUE;
    Vtrj_DrawAll();
	isSwitching=FALSE;
        
        /*Vtrj_ClearXYZ();*/

        	
	isRused = Visible(editRinfo_win);
	if (IsSwitching==TRUE) AAmodified=FALSE;
	
	if (isRused==TRUE) LoadREditProc(NULL);
	
	if (isAAused==TRUE)
	{
		Reset(modifyRlist);
		LoadAAList(GetAAPos());	
	}
		
	Hide(aaname_panel);
	Show(aaname_panel);
	

return numrecords;
	
}

static void PreSwitchResProc(IteM i)
{
	SwitchResProc(TRUE);
}

static void TextEnteredProc (IteM i)
{
	Char str[5];
	GetTitle (jumpto_text, str, sizeof(str));
	
	if (StringLen(str) == 0)
	{
		Disable (go_bttn);
	}
	else
	{
		Enable (go_bttn);
	}
	
	
}


static void GoBttnPressedProc (IteM i)
{
	int ijump_num,err;
	Int4 jump_num;
	Char jump_to_where[5];
	Int2 count;
	Boolean isok;
	
	isok = TRUE;
	
	GetTitle (jumpto_text, jump_to_where, sizeof(jump_to_where));
	
	for (count = 0; count <= StringLen(jump_to_where); count++)
	{
		if (jump_to_where[count] > '9')
		{
			Message (MSG_OK, "You have entered an invalid residue number.");
			SetTitle (jumpto_text, "");
			Disable (go_bttn);
			isok = FALSE;
			break;
		}
	 }
	
	 err=sscanf(jump_to_where, "%d", &ijump_num);
	 if (err==0)
	 	ijump_num=0;
	 jump_num = (Int4) ijump_num;
	
	 if (isok==TRUE && (jump_num > numAA || jump_num < 1))
	 {
	 	Message (MSG_OK, "No such residue number in this sequence.");
	 	SetTitle (jumpto_text, "");
	 	Disable (go_bttn);
	 	isok = FALSE;
	 }
	
	 if (isok == TRUE)
	 {	
	 	res_num = jump_num;
	 	SetTitle (jumpto_text, "");
	 	Disable(go_bttn);	 	
	 	SetValue (resid_switch, jump_num);
	 	SwitchResProc(TRUE);
	
	 }
}




/*****************************************************************
                                          TRAJ-GRAPH LOAD & SAVE FUNCTIONS
*****************************************************************/

static void Vtrj_LoadSequence(void)
{
	Xsequence = DecodeSequence(sequence, EXTAA_X);
	Psequence = DecodeSequence(sequence, EXTAA_PARENT);
}

static void Vtrj_UnLoadSequence(void)
{
	Xsequence = MemFree(Xsequence);
	Psequence = MemFree(Psequence);
	shortsequence[0] = '\0';
	sequence[0] = '\0';
	newsequence = MemFree(newsequence);
}

static void GetNewFileProc (IteM i)
{
    Char path[PATH_MAX];
	path[0] = '\0';
	
	isNew = FALSE;

	if (GetInputFileName(path, sizeof(path), "trj*", NULL))
	{
		StringCpy(GlobalPath, path);
		LoadProc(path);
	}
}

static void LoadProc (CharPtr path)
{
	Int4 numrecords;
    CharPtr main_title, hastrj=NULL, tmp_fnamtrj;
	Char tmp_progname[PATH_MAX];
	Boolean is_OK=FALSE;
	Int2 cnt;

	is_OK = CheckIfSaved();
	if (is_OK==FALSE) return;
	if (hasafrag!=NULL)
		hasafrag=MemFree(hasafrag);
	StringCpy(tmp_progname, prog_name);
/*	tmp_progname[20] = '\0';*/

	tmp_fnamtrj = StringSave(path);
	
	if (fnamtrj!=NULL)
		fnamtrj=MemFree(fnamtrj);
	if (isNew==TRUE)
	{
		fnamtrj=(CharPtr)MemNew(5+StringLen(tmp_fnamtrj));
		StringCpy(fnamtrj,tmp_fnamtrj);
		StringCat(fnamtrj, ASN_EXT);
		StringCpy(GlobalPath, path);
	}
	else if (isNew==FALSE) fnamtrj = StringSave(tmp_fnamtrj);

	hastrj = StringStr (fnamtrj, ASN_EXT);
               	
	if (hastrj==NULL) {
		ErrPostEx(SEV_ERROR,1,1,"Invalid file name");
		return;
	}
	else if (hastrj !=NULL) hastrj[0] = '\0';
		
	main_title=(CharPtr)MemNew(StringLen(tmp_progname)+StringLen(fnamtrj)+1);
	StringCpy(main_title,tmp_progname);
	StringCat(main_title,fnamtrj);
	
	if (UnPackAsnTrajGraph(fnamtrj,&numAA,sequence,NULL,&pbi,&pbd,&vnpBioseq)==NULL) {
   	return;
 	}
 	if (hasafrag==NULL)
 		hasafrag=(Boolean *)MemNew(numAA*sizeof(Boolean));
	TGInit(tmpdbasename, DB_READ, &numrecords);
 	for (cnt=0;cnt<numAA;cnt++) {
  	resinfo=TrajGraphRead(cnt+1);
  	if (resinfo->pfdsFragmentHead!=NULL)
  		hasafrag[cnt]=TRUE;
  	resinfo=FreeTraj(resinfo);
  }
	TGClose();        	
	Vtrj_LoadSequence();

	if (WALKTYPE==WALK_CA) StringCpy(TRAJ3FNAME,"cawalk_dict");
  else if (WALKTYPE==WALK_PHIPSI) StringCpy(TRAJ3FNAME,"phipsiwalk_dict");
	
	SetTitle(main_win, main_title);
	
	isLoaded = TRUE;
	
	
	SetValue(which_view, 2);
	
	res_num = 1;

	SetSwitchParams (resid_switch, 1, numAA);
	
	numrecords=SwitchResProc (TRUE);

	if (numrecords != numAA)
		ErrPostEx(SEV_ERROR,1,1,"Number of records does not match number of amino acids");
		
	SaveOrigGlobs();
	Enable(jumpto_text);

	Vtrj_pnn = GetTrueDistConstraints();
	main_title=MemFree(main_title);
	tmp_fnamtrj=MemFree(tmp_fnamtrj);	
}

static Boolean SaveTrajProc (void)
{
	Char tmp_progname[PATH_MAX];
	Char path[PATH_MAX];
	CharPtr hasstr=NULL;
	CharPtr main_title;
	Boolean isOk=FALSE;
	ByteStorePtr bsid=NULL, bsdescr=NULL;
	ValNodePtr vnpDescr, vnpId;
	AsnIoBSPtr aibp;	

	StringCpy(tmp_progname,prog_name);
	path[0] = '\0';


	if (ShouldSaveAs == TRUE)
	{
		if (GetOutputFileName(path, sizeof(path), NULL)) isOk = TRUE;
		else return FALSE;
	}

	if (ShouldSaveAs == FALSE)
	{
		StringCpy(path, GlobalPath);
		isOk = TRUE;
	}

	if (isOk == TRUE)
	{
        hasstr = StringStr (path, ASN_EXT);

        if (hasstr !=NULL) hasstr[0] = '\0';

        hasstr = NULL;

		if (HasBeenFiltered) {
			Vtrj_SaveResidue();
			Disable(Vtrj_None);
		}
		/* arbitrary sizes */
		bsid=BSNew(sizeof(Char)*1);
		bsdescr=BSNew(sizeof(Char)*1);
		BSSeek(bsid,0L,SEEK_SET);
		BSSeek(bsdescr,0L,SEEK_SET);
		vnpId=pbi;
		aibp=AsnIoBSOpen("wb",bsid);
		while (vnpId!=NULL)
		{
			BiostrucIdAsnWrite(vnpId,aibp->aip,NULL);
			vnpId=vnpId->next;
		}
		aibp=AsnIoBSClose(aibp);
		vnpDescr=pbd;
		aibp=AsnIoBSOpen("wb",bsdescr);
		while (vnpDescr!=NULL)
		{
			BiostrucDescrAsnWrite(vnpDescr,aibp->aip,NULL);
			vnpDescr=vnpDescr->next;
		}
		aibp=AsnIoBSClose(aibp);


		PackAsnTrajGraph (path, sequence, bsid, bsdescr, vnpBioseq, FALSE);
		
/*        hasstr = StringStr (tmp_progname, " - ");

        if (hasstr !=NULL) hasstr[3] = '\0';
*/	
		main_title = StringSave(StringCat(tmp_progname, path));
	
	SetTitle(main_win, main_title);

	StringCpy(GlobalPath, path);

	Rhasbeenedited=FALSE;
	Ghasbeenedited=FALSE;
	Nhasbeenedited=FALSE;
	HasBeenFiltered=FALSE;
	DistBeenEdited=FALSE;
	FragBeenEdited=FALSE;

	}

	return TRUE;
}


static void Vtrj_SaveAsProc (IteM i)
{
	ShouldSaveAs = TRUE;
	SaveTrajProc();

}

static void Vtrj_SaveProc (IteM i)
{
	ShouldSaveAs = FALSE;
	SaveTrajProc();
}

/*static void SavePNGProc(IteM i)
{

          Nlm_SaveImagePNG("screenshotOGL.png");

          Message (MSG_OK, "Feature to come in version 2");

}*/

static void Vtrj_CalculateInfo(Boolean CalcTout)
{
	Int2 num = 0,nextnum;
	Int4 numrecords=0;
	FloatLo pcis;
	PTGS restemp;
	
/*	if (resinfo->pCis > 0.5 && (WALKTYPE==WALK_PHIPSI))
			num =1;*/
	if (/*WALKTYPE==WALK_CA &&*/ (res_num<numAA)) {
		if (resinfo->pCis > 0.5)
			nextnum=1;
		else
			nextnum=0;
		TGInit(tmpdbasename, DB_READ, &numrecords);	
        restemp =  TrajGraphRead((Int2) res_num+1);
		TGClose();
		TrajGraphDifferentiate(restemp);
		TrajCalcSparsity(restemp,nextnum);
		TrajCalcTout(restemp);
		/***** IMPORTANT ****/
		/* integrate them before calculating NZ and after sparsity */
		TrajGraphIntegrate(restemp);
		TrajCalcNZ(restemp);
		TGInit(tmpdbasename, DB_READ, &numrecords);
		TrajGraphWrite (restemp, USE_RLE, TRUE);
		TGClose();
		restemp=FreeTraj(restemp);
	}
	if (/*WALKTYPE==WALK_CA &&*/ (res_num>1)) {
		TGInit(tmpdbasename, DB_READ, &numrecords);
		pcis=GetPCis(res_num-1);
		TGClose();
		if (pcis>0.5)
			num=1;
		else
			num=0;
	}
	/* fill in remaining fields */
	TrajCalcSparsity(resinfo,num);
	/* make timeout depend on trajectory distribution sparsity */
	if (CalcTout==TRUE)
		TrajCalcTout(resinfo);

}

static void Vtrj_SaveResidue(void)
{
	Int4 numrecords=0;
	
	Vtrj_CalculateInfo(FALSE);
	/* integrate them */
	/***** IMPORTANT ****/
	/* integrate them before calculating NZ and after sparsity */
	TrajGraphIntegrate(resinfo);
	TrajCalcNZ(resinfo);
	/* rotamer info is reserved for future use */
/*	resinfo->rotid=0; used now */
	
	if (resinfo==NULL)
		ErrPostEx(SEV_ERROR,1,1,"Error in writing new trajectory distribution");

	TGInit(tmpdbasename, DB_READ, &numrecords);

	if (numrecords != numAA)
		ErrPostEx(SEV_ERROR,1,1,"Number of records does not match number of amino acids");
		
	
	TrajGraphWrite (resinfo, USE_RLE, ResReplace);
	TGClose();
	

	Rhasbeenedited=TRUE;
		
	SwitchResProc(FALSE);  /*Not called to switch residues but rather to reload all
				variables with the new information produced
				by FillTG*/

	Disable(Vtrj_None);

}


/*****************************************************************
			OTHER FUNCTIONS
*****************************************************************/

static void LoadPInfo(void)
{
	Boolean MmdbIn=FALSE, NameIn=FALSE, AttribIn=FALSE,
		/*HistoryIn= FALSE,*/ OthercIn=FALSE, PdbcIn=FALSE;

	BiostrucDescrPtr	tmp_pbd;
	BiostrucIdPtr	tmp_pbi;
	
	ObjectIdPtr objid;
	DbtagPtr dbt;
	ValNodePtr pvnPub=NULL;
	CharPtr pcNameString=NULL;
	CitSubPtr pcsThis=NULL;
	DataVal dvAuthors;

  Char pinfo[MAXCOL+1];

  hdgColFmt.font = ParseFont ("Times,18,b");
	subColFmt.font = ParseFont ("Helvetica,11,b");
	txtColFmt.font = ParseFont ("fixed,12");
	lstColFmt.font = programFont;
	tblColFmt.font = programFont;

		
	tmp_pbi = pbi;	

    		
	while (tmp_pbi!=NULL)
	{
		switch (tmp_pbi->choice)
		{	
			case BiostrucId_mmdb_id:
			    if (MmdbIn==FALSE) AppendText (pinfo_doc, "MMDB Id", &hdgParFmt, &subColFmt, programFont);
			    MmdbIn = TRUE;
				sprintf(pinfo, "%d",   tmp_pbi->data.intvalue);			
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				break;
			
			case BiostrucId_other_database:
				dbt = tmp_pbi->data.ptrvalue;
			
			    AppendText (pinfo_doc, "Other Database", &hdgParFmt, &subColFmt, programFont);
				StringCpy(pinfo,(CharPtr) dbt->db);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				
				objid = dbt->tag;
				
				sprintf(pinfo, "%d",   objid->id);			
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				
				StringCpy(pinfo,(CharPtr) objid->str);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				
				break;
			
			case BiostrucId_local_id:	
				objid = tmp_pbi->data.ptrvalue;
			
			    AppendText (pinfo_doc, "Local Id", &hdgParFmt, &subColFmt, programFont);
				
				sprintf(pinfo, "%d",   objid->id);			
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				
				StringCpy(pinfo,(CharPtr) objid->str);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
		}
		tmp_pbi = tmp_pbi->next; 	
	}	

	
	tmp_pbd = pbd;	

	while (tmp_pbd!=NULL)
	{
		switch (tmp_pbd->choice)
		{
			case BiostrucDescr_name:
			    if (NameIn==FALSE) AppendText (pinfo_doc, "Name", &hdgParFmt, &subColFmt, programFont);
				NameIn = TRUE;			
				StringCpy(pinfo,(CharPtr) tmp_pbd->data.ptrvalue);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				break;
				
			case BiostrucDescr_pdb_comment:
			    if (PdbcIn==FALSE) AppendText (pinfo_doc, "Pdb Comment", &hdgParFmt, &subColFmt, programFont);
			    PdbcIn = TRUE;
				StringCpy(pinfo,(CharPtr) tmp_pbd->data.ptrvalue);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
		        break;
			
			case BiostrucDescr_other_comment:
			    if (OthercIn==FALSE) AppendText (pinfo_doc, "Other Comment", &hdgParFmt, &subColFmt, programFont);
				OthercIn = TRUE;
				StringCpy(pinfo,(CharPtr) tmp_pbd->data.ptrvalue);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
			    break;
			
			case BiostrucDescr_history:
				/*if (HistoryIn==FALSE) AppendText (pinfo_doc, "History", &hdgParFmt, &subColFmt, programFont);
			    HistoryIn = TRUE;
			    StringCpy(pinfo,(CharPtr) tmp_pbd->data.ptrvalue);
				AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);*/
			    break;
			
			case BiostrucDescr_attribution:
				pvnPub = (ValNodePtr) tmp_pbd->data.ptrvalue;
				if (pvnPub->choice==PUB_Sub) {
					pcsThis=(CitSubPtr)pvnPub->data.ptrvalue;
					dvAuthors.ptrvalue=pcsThis->authors; /* get the pointer to the author id's */
					pcNameString=StdAuthListNamesPrint(&dvAuthors);  /* pretty print them */
					sprintf(pinfo,"%s",pcNameString);
					pcNameString=MemFree(pcNameString);
					if (AttribIn==FALSE)
						AppendText (pinfo_doc, "Attribution", &hdgParFmt, &subColFmt, programFont);
					AttribIn = TRUE;
					AppendText (pinfo_doc, pinfo, &txtParFmt, &txtColFmt, programFont);
				}
				break;
		}
		tmp_pbd = tmp_pbd->next;
	 }
		
}


static void DonePInfoProc(void)
{
	Hide(pinfo_win);
}


static void ProteinInfoProc (IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2;

	pinfo_win = MovableModalWindow(-50,-33,-30,-30, "General Protein Information", (WndActnProc) DonePInfoProc);

	g = HiddenGroup(pinfo_win,-1,-1, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);
		
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 3, 3);			
	
	pinfo_doc = DocumentPanel (Vtrj_grp1, stdCharWidth * 50, stdLineHeight * 20);

	Vtrj_grp2 = HiddenGroup(g,-1,-1,(GrpActnProc) NULL);
	
	DefaultButton(Vtrj_grp2, "     Close     ", (BtnActnProc) DonePInfoProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

	LoadPInfo();	
		
	Show(pinfo_win);
	
	ProcessEvents();

}



/*****************************************************************
       			      NCBI OGL FUNCTIONS
*****************************************************************/

/*void Vtrj_OGL_DrawViewer3D(TOGL_Data * OGL_Data);*/

/*Nlm_PaneL  Nlm_Autonomous3DPanel PROTO((Nlm_GrouP prnt, Nlm_Int2 pixwidth, Nlm_Int2 pixheight, Nlm_PnlActnProc draw, Nlm_SltScrlProc vscrl, Nlm_SltScrlProc hscrl, Nlm_Int2 extra, Nlm_PnlActnProc reset, Nlm_GphPrcsPtr classPtr, Nlm_Boolean *IndexMode, void **display, void **visinfo));
void Nlm_Set3DColorMap PROTO((Nlm_PaneL w, Nlm_Uint2 totalColors,
                             Nlm_Uint1Ptr red, Nlm_Uint1Ptr green,
                             Nlm_Uint1Ptr blue, void **display));

static TOGL_Data *PanelToOGL(Nlm_PaneL panel);
*/
/*static const GLfloat Color_Off[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
static const GLfloat Color_MostlyOff[4] = { 0.05f, 0.05f, 0.05f, 1.0f };
static const GLfloat Color_MostlyOn[4] = { 0.95f, 0.95f, 0.95f, 1.0f };
static const GLfloat Color_On[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
*/

static void OGL_DeleteViewer3D(TOGL_Data * OGL_Data)
/* delete OGL_Data */
{
    if (OGL_Data == NULL)
        return;

    MA_Destroy(OGL_Data->ma);

    /* free items on the heap */
    MemFree(OGL_Data->Layers);
    MemFree(OGL_Data->ModelMatrix);
    MemFree(OGL_Data);
}


                              /**********************************
                              		ScrollBar Procedures
                              **********************************/


static void OGL_ViewerVScrollProc(Nlm_BaR sb, Nlm_SlatE viewer,
                                  Nlm_Int2 newval, Nlm_Int2 oldval)
{
	if (alldisabled==TRUE)
		return;
	if (isLoaded==TRUE && ProgramProgress<=0)
	{
	    	Vtrj_rotateY +=  (newval-oldval);	
    		Vtrj_rotateY = Vtrj_rotateY % 360;
			if (QuickRotate == TRUE && Drawit == FALSE) {
				Drawit=TRUE;
	    		Vtrj_DrawAll();
				Drawit=FALSE;
			}
			else
	    		Vtrj_DrawAll();
	}
}

static void OGL_ViewerHScrollProc(Nlm_BaR sb, Nlm_SlatE viewer,
                                  Nlm_Int2 newval, Nlm_Int2 oldval)
{
	if (alldisabled==TRUE)
		return;
	if (isLoaded==TRUE && ProgramProgress<=0)
	{
	    	Vtrj_rotateX +=  -(newval-oldval);
		Vtrj_rotateX = Vtrj_rotateX % 360;    	
		if (QuickRotate == TRUE && Drawit == FALSE) {
			Drawit=TRUE;
    		Vtrj_DrawAll();
			Drawit=FALSE;
		}
		Vtrj_DrawAll();
	}
}

/**********************************
    		Mouse Procedures
**********************************/

static TOGL_Data *MAToOGL(MAPtr ma)
/* extracts OGL_Data out of the extra pointer in the mouse data */
{
    return (TOGL_Data *) MA_GetExtra(ma);
}


static TOGL_Data *PanelToOGL(Nlm_PaneL panel)
/* extract OGL_Data out of the panel data */
{
    MAPtr ma;
    Nlm_GetPanelExtra(panel, &ma);

    return MAToOGL(ma);
}


/**********************************
		CallBack Procedures
**********************************/


static void OGL_ResetViewerProc_CB(Nlm_PaneL panel)
{
    TOGL_Data *OGL_Data = PanelToOGL(panel);
    if (!OGL_Data)
        return;

    OGL_Data = NULL;
    OGL_DeleteViewer3D(PanelToOGL(panel));
}



static void Vtrj_OGL_DrawViewer3D_CB(Nlm_PaneL panel)
/* callback */
{
	/*Exposed Events*/
    if (isLoaded==FALSE) Vtrj_ClearScreen();
    if (isLoaded==TRUE) 
	{
		Vtrj_redraw = TRUE;
		Vtrj_DrawAll();
		Vtrj_redraw = FALSE;
	}
/*    Vtrj_OGL_DrawViewer3D(PanelToOGL(panel));*/
}

TOGL_Data *Vtrj_OGL_CreateViewer(Nlm_GrouP prnt,
                            Uint2Ptr width, Uint2 height,
                            Int4 flags,
                            Nlm_MenU ma_group_menu,
                            Nlm_MenU ma_action_menu,
                            Nlm_MAInitOGLFunc ma_init_func,
                            VoidPtr ma_init_data)
                             /* initialize the OpenGL library */
{
    TOGL_Data *OGL_Data;
    Nlm_Uint2 x_width;

    OGL_Data = (TOGL_Data *) MemNew(sizeof(TOGL_Data));
    if (OGL_Data == NULL)
    	return NULL;

    OGL_Data->ModelMatrix = (Nlm_VoidPtr) MemNew(16 * sizeof(GLdouble));
    if (OGL_Data->ModelMatrix == NULL)
        return NULL;

    OGL_Data->Layers = (TOGL_Layers *) MemNew(sizeof(TOGL_Layers));
    if (OGL_Data->Layers == NULL)
        return NULL;
    	
    	    	
    OGL_Data->ParentWindow = Nlm_ParentWindow((Nlm_Handle) prnt);

   OGL_Data->Panel = Nlm_Autonomous3DPanel(prnt,
                                            (Int2) * width, (Int2) height,
                                            Vtrj_OGL_DrawViewer3D_CB,
                                            ((flags & Y_ROTATE_SBAR) ?
                                             OGL_ViewerVScrollProc : NULL),
                                            ((flags & X_ROTATE_SBAR) ?
                                             OGL_ViewerHScrollProc : NULL),
                                            sizeof(MAPtr),
                                            OGL_ResetViewerProc_CB, NULL,
                                            &OGL_Data->IndexMode, &OGL_Data->display, &OGL_Data->visinfo);

     {
        Nlm_RecT rect;
        Nlm_GetPosition(OGL_Data->Panel, &rect);
        rect.right = (Int2) (rect.left + *width);
        rect.bottom = (Int2) (rect.top + height);
        OGL_SetPosition3D(OGL_Data, &rect);
        x_width = (Uint2) (rect.right - rect.left);
    }

    if (flags & X_ROTATE_SBAR) {
        Nlm_BaR sb = Nlm_GetSlateHScrollBar((Nlm_SlatE) OGL_Data->Panel);
        Nlm_CorrectBarValue(sb, 0);
        Nlm_SetRange(sb, 10, 10, 360);
    }
    if (flags & Y_ROTATE_SBAR) {
        Nlm_BaR sb = Nlm_GetSlateVScrollBar((Nlm_SlatE) OGL_Data->Panel);
        Nlm_CorrectBarValue(sb, 0);
        Nlm_SetRange(sb, 10, 10, 360);
    }
		
    *width=x_width;

    return OGL_Data;
}


void Cn3DResizeProc(WindoW w)
{
    RecT r;
    ObjectRect(w, &r);

    OffsetRect(&r, (Int2) (-r.left), (Int2) (-r.top+40));

	InsetRect(&r, 5, 50);
	OGL_SetPosition3D(myOGL_data, &r);

}


/****************************************************************************************
 							BACKGROUND COLOR FUNCTIONS
 *****************************************************************************************/

static void DoneBG(IteM i)
{
 	Hide(Background_win);
}

static void EnableBGOK(IteM i)
{
 	Enable(BGOK_bttn);
}

static void ChangeBackground(IteM i)
{
 	Int2 which_color=0;

 	which_color = GetValue(background_popup);

 	switch (which_color)
 	{
 		case 1:
 			Vtrj_bgcolor=BLACK;
 			break;
 		case 2:
 			Vtrj_bgcolor=WHITE;
 			break;
 	}


 	glClearColor(Vtrj_RGB[Vtrj_bgcolor][0], Vtrj_RGB[Vtrj_bgcolor][1], Vtrj_RGB[Vtrj_bgcolor][2], 0.0);

 	if(isLoaded==TRUE) Vtrj_DrawAll();
 	else Vtrj_ClearScreen();

 	DoneBG(i);
 }


static void DisplayBackGroundProc(IteM i)
{
 	GrouP g, Vtrj_grp1, Vtrj_grp2;
 	
 	Background_win = MovableModalWindow(-50,-33,-10,-10, "Background Color",  (WndActnProc) DoneBG);
 	
 	g = HiddenGroup(Background_win, -1, -1, (GrpActnProc) NULL);

 	SetGroupSpacing(g, 3, 15);
 	SetGroupMargins(g, 3, 3);

 	Vtrj_grp1 = HiddenGroup(g, 0, 2, (GrpActnProc) NULL);

 	StaticPrompt (Vtrj_grp1, "Choose a background color: ", 0, dialogTextHeight, systemFont, 'l');
         	
  background_popup = PopupList(Vtrj_grp1, TRUE, (PupActnProc) EnableBGOK);

 	PopupItem(background_popup, "Black");
 	PopupItem(background_popup, "White");

 	switch (Vtrj_bgcolor)
 	{
 		case BLACK:
 			SetValue(background_popup, 1);
 			break;
 		case WHITE:
 			SetValue(background_popup, 2);
 			break;
 	}
 	Vtrj_grp2 = HiddenGroup(g, 2, -1, (GrpActnProc) NULL);
 	
 	BGOK_bttn = PushButton(Vtrj_grp2, "Apply", (BtnActnProc) ChangeBackground);
 	Disable(BGOK_bttn);

 	DefaultButton(Vtrj_grp2, "Cancel", (BtnActnProc) DoneBG);
                	
 	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

     Show(Background_win);
 	
 	ProcessEvents();

}

/*****************************************************************
                         MOUSE FUNCTIONS
*****************************************************************/
void Vtrj_HoldProc (PaneL p, PoinT pt)
{
	/*Empty Function for now*/ 	
}

void Vtrj_ReleaseProc (PaneL p, PoinT pt)
{
  if (alldisabled==TRUE)
  	return;
	if (isLoaded==FALSE || ProgramProgress>0)
		return;
	if (dragged==FALSE)
	{
		Vtrj_GetPoints();
	}
	
	if (dragged==TRUE && QuickRotate==TRUE)
	{
		Drawit = TRUE;
		Vtrj_DrawAll();	
	}

}


void Vtrj_ClickProc (PaneL p, PoinT pt)
{
  if (alldisabled==TRUE)
  	return;
	whereclicked.y = pt.y;
	whereclicked.x = pt.x;
	clicked = TRUE;
	dragged = FALSE;
	if (QuickRotate==TRUE) Drawit = FALSE;
}

void Vtrj_DragProc (PaneL p, PoinT pt)
{
  if (alldisabled==TRUE)
  	return;
	wheredragged.y = pt.y;
	wheredragged.x = pt.x;

	
if (isLoaded == TRUE && ProgramProgress<=0)
	{
		if (shftKey)
		{
			Vtrj_Move();
		}	
	
		if (ctrlKey)
		{
			Vtrj_Zoom();
		}
	
		if (!ctrlKey && !shftKey)
		{
			Vtrj_Rotate();

		}
			
		Vtrj_DrawAll();	
	}
	dragged = TRUE;
}

/*****************************************************************
                 OPEN GL TRANSFORMATION FUNCTIONS
*****************************************************************/
/*static void Vtrj_ClearXYZ (void)
{
               	
               	SetTitle (xval_title, "");
               	SetTitle (yval_title, "");
               	SetTitle (zval_title, "");

}*/

/*static void Vtrj_PrintXYZ (int i)
{

	Ppp PVp = GPVp;
	GLfloat x= 0.00, y = 0.00, z = 0.00;
	CharPtr xtext, ytext, ztext;
	
	if (i > reseditinfo->dim*reseditinfo->dim)
               		ErrPostEx(SEV_ERROR,1,1,"Point selected outside array boundaries");
	
	if (PVp[i].z1 <= PVp[i].z2)
	{
	 	x = PVp[i].x2;
	 	y = PVp[i].y2;
	 	z = PVp[i].z2;
	 	
	}
	else
	{
		x = PVp[i].x1;
		y = PVp[i].y1;
		z = PVp[i].z1;
	}
	if (PVp[i].z3 >= z)
	{
		x = PVp[i].x3;
		y = PVp[i].y3;
		z = PVp[i].z3;
	}
	if (PVp[i].z4 >= z)
	{
		x = PVp[i].x4;
		y = PVp[i].y4;
		z = PVp[i].z4;
	}
	
	if (WALKTYPE==WALK_CA)
	{ 		
	 	x *=  360.0;
	 	x /= (float) resinfo->dim;
	 	
	 	y -= 1.0;
	 	y *= 2.0;
	 	y /= (float) resinfo->dim;
	}
	
	if (WALKTYPE==WALK_PHIPSI)
	{ 	
		x -= 180.0;
		x *= 360.0;
		x /= (float) resinfo->dim;
		
		y -= 180.0;
		y *= 360.0;
		y /= (float) resinfo->dim;
	}
	
	z *= ScaleZ;
		
        	xtext =  Ltostr((Int4) x,(Int2) NULL);
        	ytext =  Ltostr((Int4) y,(Int2) NULL);
        	ztext =  Ltostr((Int4) z,(Int2) NULL);
        	
               	SetTitle (xval_title, xtext);
               	SetTitle (yval_title, ytext);
               	SetTitle (zval_title, ztext);

} */


static void Vtrj_GetPoints(void)
{

	Boolean isSphere = FALSE;
	GLuint hit/*, hits*/;
	GLuint selectBuf[100];
                Int2 viewtype;
	
	viewtype = GetValue(which_view);
	
	if (viewtype==1) isSphere = TRUE;
		
	glGetIntegerv(GL_VIEWPORT, viewport_properties);
	
	
	glSelectBuffer(100, selectBuf);	
	glRenderMode(GL_SELECT);
  	glInitNames();
	glPushName(~0);

  	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPickMatrix((GLdouble) whereclicked.x, (GLdouble)
	(viewport_properties[3] - whereclicked.y), 2.0, 2.0,
	viewport_properties);
	glOrtho(-Vtrj_scale, (Vtrj_scale*Vtrj_Xratio), (-Vtrj_scale*Vtrj_Yratio), Vtrj_scale, -50.0, 50.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	gluLookAt(sin(DEGTORAD*Vtrj_rotateY)*cos(DEGTORAD*Vtrj_rotateX),sin(DEGTORAD*Vtrj_rotateY)*sin(DEGTORAD*Vtrj_rotateX),cos(DEGTORAD*Vtrj_rotateY),
	    0.0, 0.0, 0.0, sin(DEGTORAD*(Vtrj_rotateY+90.0))*cos(DEGTORAD*Vtrj_rotateX),sin(DEGTORAD*(Vtrj_rotateY+90.0))*sin(DEGTORAD*Vtrj_rotateX),cos(DEGTORAD*(Vtrj_rotateY+90.0)));
	
  	Vtrj_Draw_Poly(GL_SELECT, isSphere);

	hit = glRenderMode(GL_RENDER);
	
	glPopMatrix();
 	if ( (int) hit <=0)
 	{
 		/*Vtrj_ClearXYZ();*/
 		return;
 	}
 	
/* 	hits = selectBuf[( (int) hit - 1) * 4 + 3];
 	
	Vtrj_PrintXYZ((int) hits);*/
}



static void Vtrj_Move(void)
{
	static Int2 disx= 0, disy = 0;
	static PoinT lastposition;

	
	if (clicked==TRUE)
	{
		lastposition = whereclicked;
		clicked = FALSE;
	}  	
	
	
	disx =  (wheredragged.x - lastposition.x)/MOUSE_SENSITIVITY;
	disy =  (wheredragged.y - lastposition.y)/MOUSE_SENSITIVITY;
	
	
		glPopMatrix();
        glTranslatef(disx, -disy, 0.0);
		glPushMatrix();
	
        lastposition = wheredragged;


}

static void Vtrj_Zoom(void)
{

	static Int2 disx= 0;
	static PoinT lastposition;

	
	if (clicked==TRUE)
	{
		lastposition = whereclicked;
		clicked = FALSE;
	}  	
	
	
	disx =  (wheredragged.x - lastposition.x);

	if (disx >= 0) Vtrj_scale *= 1.1;
	
	if (disx < 0) Vtrj_scale *= 0.9;
	
                                lastposition = wheredragged;
	
}


static void Vtrj_Rotate(void)
{
	static PoinT lastposition;
	
	if (clicked==TRUE)
	{
		lastposition = whereclicked;
		clicked = FALSE;
	}  	

	
	
	Vtrj_rotateX += (lastposition.x - wheredragged.x);
	Vtrj_rotateY += -(lastposition.y - wheredragged.y);

	Vtrj_rotateY = Vtrj_rotateY % 360;
	Vtrj_rotateX = Vtrj_rotateX % 360;	

                lastposition = wheredragged;
}


/*****************************************************************
                  SHADING FUNCTIONS
*****************************************************************/

static void InitializeShading(void)
{


	/*glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);*/
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lightModel);


	glLightfv(GL_LIGHT0, GL_AMBIENT, lightColorAmb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColorDif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightColorSpec);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
		
/*	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor2);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition2);*/

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
/*	glEnable(GL_LIGHT1);*/


}


/*****************************************************************
                     OPEN GL DRAWING FUNCTIONS
*****************************************************************/


static void Vtrj_LoadVertexes(void)
{

	Int2 viewtype=2;
	
	viewtype = GetValue(which_view);
	
  	if (viewtype == 2) Vtrj_LoadVertexesPlane();
  	if (viewtype == 1) Vtrj_LoadVertexesSphere();

}


static void Vtrj_LoadVertexesSphere(void)
{

	Ppp PVp;

                Int2 count;
	Int2 x, y;
	Int4 *TrajHere=NULL;
	FloatLo maxpoint=0.00;
	FloatLo theta, theta2, cosphi, cosphi2,  sinphi, sinphi2,
		radius1, radius2, radius3, radius4, radius = 5.0;
		
	count = resinfo->dim;
	
	if (GPVp!=NULL) GPVp = MemFree(GPVp);
	GPVp = (Ppp) MemNew(count*count*sizeof(PolyPoints));
	PVp = GPVp;
		
                if (isCis==FALSE)	TrajHere = resinfo->TrajGraph;
                if (isCis==TRUE) TrajHere = resinfo->CisTrajGraph;

	/*For both of the for loops, count-2 is used to remove 2 large residual spikes
	at the south pole (appears if count-1 is used)*/

	for (y = 0; y<= count-1; y++)
	{
		cosphi =  ((FloatLo) y/ (FloatLo) count)*2.00-1.0;
		cosphi2 = ((FloatLo) (y+1)/ (FloatLo) count)*2.00-1.0;


		sinphi  = sin(acos(cosphi));
		sinphi2 = sin(acos(cosphi2));
		
		for (x = 0;  x<= count-1; x++)
		{
			
			theta =  ((FloatLo) x/(FloatLo) count)*2.00*PI-PI;
			theta2 = ((FloatLo) (x+1)/(FloatLo) count)*2.00*PI-PI;

			
			radius1 = radius + (FloatLo) *(TrajHere)/ScaleSphere;
			if(y==count-1) 
				radius2 = radius + (FloatLo) *(TrajHere-count*y)/ScaleSphere;
			else 
			    radius2 = radius + (FloatLo) *(TrajHere+count)/ScaleSphere;
			if(x==count-1)
				radius3 = radius + (FloatLo) *(TrajHere-x)/ScaleSphere;
			else
				radius3 = radius + (FloatLo) *(TrajHere+1)/ScaleSphere;
			if(x==count-1 && y==count-1)
				radius4 = radius + (FloatLo) *(TrajHere-count*y-x)/ScaleSphere;
			else if(x==count-1)
				radius4 = radius + (FloatLo) *(TrajHere+1)/ScaleSphere;
			else if(y==count-1)
				radius4 = radius + (FloatLo) *(TrajHere-count*y+1)/ScaleSphere;
			else
				radius4 = radius + (FloatLo) *(TrajHere+count+1)/ScaleSphere;

  		PVp->x1 = radius1*sinphi*cos(theta);
			PVp->y1 = radius1*sinphi*sin(theta);
			PVp->z1 = radius1*cosphi;
			
 		  PVp->x2 = radius2*sinphi2*cos(theta);
			PVp->y2 = radius2*sinphi2*sin(theta);
			PVp->z2 = radius2*cosphi2;
			
			PVp->x3 = radius3*sinphi*cos(theta2);
			PVp->y3 = radius3*sinphi*sin(theta2);
			PVp->z3 = radius3*cosphi;
			
			PVp->x4 = radius4*sinphi2*cos(theta2);
			PVp->y4 = radius4*sinphi2*sin(theta2);
			PVp->z4 = radius4*cosphi2;			
			
			
			if (radius1 <= radius2) maxpoint = radius2;
			if (radius1 >= radius2) maxpoint = radius1;
			if (radius3 >= maxpoint) maxpoint = radius3;
			if (radius4 >= maxpoint) maxpoint = radius4;

			maxpoint = (maxpoint-radius)*ScaleSphere;
			
			PVp->maxpoint = maxpoint;
			
			PVp->color = RED;			
			if (maxpoint <= (FloatLo) ThePeak*0.75) PVp->color = RED;
			if (maxpoint <= (FloatLo) ThePeak*0.50) PVp->color = YELLOW;
			if (maxpoint <= (FloatLo) ThePeak*0.25) PVp->color = GREEN;
			
			if (radius1==radius && radius2==radius && radius3==radius && radius4==radius) PVp->iszero = TRUE;
			else PVp->iszero = FALSE;
			
		TrajHere++;		
		PVp++;
		}
	}
}


static void Vtrj_LoadVertexesPlane(void)
{
	
	Ppp PVp;
	
                Int2 count, position;
                GLfloat X, X1, Y, Y1;
	Int2 x, y;
	Int4 *TrajHere=NULL;
	GLfloat maxpoint=0.00;

	count = resinfo->dim;
	position = count/2;
	
	if (GPVp!=NULL) GPVp = MemFree(GPVp);
	GPVp = (Ppp) MemNew(count*count*sizeof(PolyPoints));
	PVp = GPVp;
	
	if (isCis==FALSE)	TrajHere = resinfo->TrajGraph;	
     	if (isCis==TRUE)	TrajHere = resinfo->CisTrajGraph;
	
	for (y = 0; y<= count-2; y++)
	{
		Y = (GLfloat) (y-position)/SCALE_SIZE_XY;
		Y1 = (GLfloat) (y+1-position)/SCALE_SIZE_XY;
		
		for (x = 0;  x<= count-2; x++)
		{
			
			X = (GLfloat) (x-position)/SCALE_SIZE_XY;
			X1 = (GLfloat) (x+1-position)/SCALE_SIZE_XY;
						
			PVp->x1 = X;
			PVp->y1 = Y;
			PVp->z1 = (GLfloat) *(TrajHere)/ScaleZ;
			
			PVp->x2 = X;
			PVp->y2 = Y1;
			PVp->z2 = (GLfloat) *(TrajHere+count)/ScaleZ;
			
			PVp->x3 = X1;
			PVp->y3 = Y;
			PVp->z3 = (GLfloat) *(TrajHere+1)/ScaleZ;
			
			PVp->x4 = X1;
			PVp->y4 = Y1;
			PVp->z4 = (GLfloat) *(TrajHere+count+1)/ScaleZ;			
			
			
			if (PVp->z1 <= PVp->z2) maxpoint = PVp->z2;
			else maxpoint = PVp->z1;
			if (PVp->z3 >= maxpoint) maxpoint = PVp->z3;
			if (PVp->z4 >= maxpoint) maxpoint = PVp->z4;
			
			maxpoint = maxpoint*ScaleZ;
			
			PVp->maxpoint = maxpoint;
			
			PVp->color = RED;			
			if (maxpoint <= (GLfloat) ThePeak*0.75) PVp->color = RED;
			if (maxpoint <= (GLfloat) ThePeak*0.50) PVp->color = YELLOW;
 			if (maxpoint <= (GLfloat) ThePeak*0.25) PVp->color = GREEN;
 			
 			if (PVp->z1 + PVp->z2 + PVp->z3 + PVp->z4==0.00) PVp->iszero = TRUE;
 			else PVp->iszero = FALSE;
 						
			TrajHere++;
			PVp++;
		}
		TrajHere++;	
	}
}

static void InitializeOGLDraw(void)
{

			glPushMatrix();
			
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			
			glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_LESS);
			
			glCullFace(GL_BACK);
/*			glCullFace(GL_FRONT);*/
			
			if (ShadingOn==TRUE) InitializeShading();

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho(-Vtrj_scale, (Vtrj_scale*Vtrj_Xratio), (-Vtrj_scale*Vtrj_Yratio), Vtrj_scale, -50.0, 50.0);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glTranslatef(2.0, 0.0, 0.0);
			gluLookAt(sin(DEGTORAD*Vtrj_rotateY)*cos(DEGTORAD*Vtrj_rotateX),sin(DEGTORAD*Vtrj_rotateY)*sin(DEGTORAD*Vtrj_rotateX),cos(DEGTORAD*Vtrj_rotateY),
/* eye is at (0,0,30) */
	    0.0, 0.0, 0.0,      /* center is at (0,0,0) */

sin(DEGTORAD*(Vtrj_rotateY+90.0))*cos(DEGTORAD*Vtrj_rotateX),sin(DEGTORAD*(Vtrj_rotateY+90.0))*sin(DEGTORAD*Vtrj_rotateX),cos(DEGTORAD*(Vtrj_rotateY+90.0)));
/* up is in positive Y direction */
}             	

static void Vtrj_DrawNorthPole(void)
{
	Int2 color = MAGENTA;
	GLUquadricObj *quad;
	
	if (ShadingOn==FALSE) glColor3f(Vtrj_RGB[color][0], Vtrj_RGB[color][1], Vtrj_RGB[color][2]);
	else
	{
		glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbPurple);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifPurple);
		glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecPurple);
		glMaterialfv(GL_FRONT, GL_SHININESS, shinePurple);
    }
	
	quad = gluNewQuadric();
	
	gluQuadricDrawStyle(quad, filltype);
	
	glTranslatef(0.0, 0.0, 5.0);	
	gluCylinder(quad, 0.1, 0.1, 2, 10, 1);
	glTranslatef(0.0, 0.0,  2.0);
	gluCylinder(quad, 0.2, 0.0, 1, 5, 1);
	glTranslatef(0.0, 0.0,-2.0);
	glTranslatef(0.0, 0.0, -5.0);
}

static void Vtrj_NumberAxis(Int2 color)
{

	if (WALKTYPE != WALK_CA)
	{		
				/*Y-AXIS*/
		Vtrj_PutText("Psi", color, -8.2, 0.0, 0.003, 1.0, 0.0);
		Vtrj_PutText("-180", color, -8.0, -6.8, 0.002, 0.4, 0.0);
		Vtrj_PutText("-144", color, -8.0, -5.48, 0.002, 0.4, 0.0);
		Vtrj_PutText("-108", color, -8.0, -4.16, 0.002, 0.4, 0.0);
		Vtrj_PutText("-72", color, -7.85, -2.84, 0.002, 0.4, 0.0);
		Vtrj_PutText("-36", color, -7.85, -1.54, 0.002, 0.4, 0.0);
		Vtrj_PutText("0", color, -7.5, -0.2, 0.002, 0.4, 0.0);
		Vtrj_PutText("36", color, -7.7, 1.12, 0.002, 0.4, 0.0);
		Vtrj_PutText("72", color, -7.7, 2.44, 0.002, 0.4, 0.0);
		Vtrj_PutText("108", color, -7.8, 3.76, 0.002, 0.4, 0.0);
		Vtrj_PutText("144", color, -7.8, 5.08, 0.002, 0.4, 0.0);
		Vtrj_PutText("180", color, -7.8, 6.4, 0.002, 0.4, 0.0);
	
		Vtrj_PutText("-", color, -7.25, -6.8, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -5.48, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -4.16, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -2.84, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -1.54, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -0.2, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 1.12, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 2.44, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 3.76, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 5.08, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 6.4, 0.003, 1.4, 0.0);


				/*X-AXIS*/
		Vtrj_PutText("Phi", color, 0.0, -8.2, 0.003, 1.0, 0.0);
		Vtrj_PutText("-180", color, -7.0, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("-144", color, -5.66, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("-108", color, -4.32, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("-72", color, -2.98, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("-36", color, -1.64, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("0", color, -0.3, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("36", color, 1.04, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("72", color, 2.38, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("108", color, 3.72, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("144", color, 5.06, -7.7, 0.002, 0.4, 0.0);
		Vtrj_PutText("180", color, 6.4, -7.7, 0.002, 0.4, 0.0);
	
		Vtrj_PutText("|", color, -6.55, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -5.21, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -3.87, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -2.7, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -1.3, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -0.25, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 1.2, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 2.5, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 3.9, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 5.3, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 6.6, -7.25, 0.002, 1.4, 0.0);
		
	
	
	}

	if (WALKTYPE == WALK_CA)
	{
				/*Y-AXIS*/
		Vtrj_PutText("Cos (Phi)", color, -7.9, 0.0, 0.003, 1.0, 90.0);
		Vtrj_PutText("-1.0", color, -7.9, -6.8, 0.002, 0.4, 0.0);
		Vtrj_PutText("-0.8", color, -7.9, -5.48, 0.002, 0.4, 0.0);
		Vtrj_PutText("-0.6", color, -7.9, -4.16, 0.002, 0.4, 0.0);
		Vtrj_PutText("-0.4", color, -7.9, -2.84, 0.002, 0.4, 0.0);
		Vtrj_PutText("-0.2", color, -7.9, -1.54, 0.002, 0.4, 0.0);
		Vtrj_PutText("0.0", color, -7.7, -0.2, 0.002, 0.4, 0.0);
		Vtrj_PutText("0.2", color, -7.7, 1.12, 0.002, 0.4, 0.0);
		Vtrj_PutText("0.4", color, -7.7, 2.44, 0.002, 0.4, 0.0);
		Vtrj_PutText("0.6", color, -7.7, 3.76, 0.002, 0.4, 0.0);
		Vtrj_PutText("0.8", color, -7.7, 5.08, 0.002, 0.4, 0.0);
		Vtrj_PutText("1.0", color, -7.7, 6.4, 0.002, 0.4, 0.0);

		Vtrj_PutText("-", color, -7.25, -6.8, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -5.48, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -4.16, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -2.84, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -1.54, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, -0.2, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 1.12, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 2.44, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 3.76, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 5.08, 0.003, 1.4, 0.0);
		Vtrj_PutText("-", color, -7.25, 6.4, 0.003, 1.4, 0.0);
			

				/*X-AXIS*/
		Vtrj_PutText("Theta", color, 0.0, -8.6, 0.003, 1.0, 0.0);
		Vtrj_PutText("-180", color, -6.7, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("-144", color, -5.66, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("-108", color, -4.32, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("-72", color, -2.98, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("-36", color, -1.64, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("0", color, -0.3, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("36", color, 1.04, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("72", color, 2.38, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("108", color, 3.72, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("144", color, 5.06, -7.65, 0.002, 0.4, 0.0);
		Vtrj_PutText("180", color, 6.4, -7.65, 0.002, 0.4, 0.0);
	
		Vtrj_PutText("|", color, -6.65, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -5.5, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -4.15, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -2.8, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -1.4, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, -0.15, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 1.2, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 2.5, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 3.9, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 5.3, -7.25, 0.002, 1.4, 0.0);
		Vtrj_PutText("|", color, 6.6, -7.25, 0.002, 1.4, 0.0);
	
	}


}

static void Vtrj_DrawAxis(FloatLo length)
{
	Int2 color = 0;
	GLUquadricObj *quad;
	
	if(Vtrj_bgcolor==WHITE) color = BLACK;
 	else if(Vtrj_bgcolor==BLACK) color = WHITE;
	length *= 2.00;
	
	if (ShadingOn==FALSE) glColor3f(Vtrj_RGB[color][0], Vtrj_RGB[color][1], Vtrj_RGB[color][2]);
	else
	{
		if (color==BLACK) {
			glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbBlack);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifBlack);
			glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecBlack);
			glMaterialfv(GL_FRONT, GL_SHININESS, shineBlack);
		}		
		else { /* color==WHITE */
			glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbWhite);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifWhite);
			glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecWhite);
			glMaterialfv(GL_FRONT, GL_SHININESS, shineWhite);
		}		
    }
	
	quad = gluNewQuadric();
	
	gluQuadricDrawStyle(quad, filltype);
	
	glTranslatef(-7.0, -7.0, 0.0);	
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	
	gluCylinder(quad, 0.05, 0.05, (length + 1), 50, 1);
	glTranslatef(0.0, 0.0, (GLfloat) (length+1));
	gluCylinder(quad, 0.1, 0.0, (length + 1) /18, 50, 1);
	glTranslatef(0.0, 0.0, (GLfloat) -(length+1));
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glRotatef(90.0, 0.0, 1.0, 0.0);
	
	gluCylinder(quad, 0.05, 0.05, (length+1), 50, 1);
	glTranslatef(0.0, 0.0, (GLfloat) (length+1));
	gluCylinder(quad, 0.1, 0.0, (length+1)/18, 50, 1);
	glTranslatef(0.0, 0.0, (GLfloat) -(length+1));
	glRotatef(-90.0, 0.0, 1.0, 0.0);
               	glTranslatef(7.0, 7.0, 0.0);	

    

	Vtrj_NumberAxis(color);				

}

static void Vtrj_DrawBase(FloatLo pos)
{
	Int2 color = BLUE;
	

                if (ShadingOn==FALSE) glColor3f(Vtrj_RGB[color][0], Vtrj_RGB[color][1], Vtrj_RGB[color][2]);
                else
                {
				    glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbBlue);
					glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifBlue);
				    glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecBlue);
				    glMaterialfv(GL_FRONT, GL_SHININESS, shineBlue);
                }

	glBegin(GL_QUADS);
	glVertex3f((GLfloat) (-pos), (GLfloat) (pos), 0.00);
	glVertex3f((GLfloat) (-pos), (GLfloat) (-pos), 0.00);
	glVertex3f((GLfloat) (pos), (GLfloat) (-pos), 0.00);
	glVertex3f((GLfloat) (pos), (GLfloat) (pos), 0.00);
	glEnd();

	glBegin(GL_QUADS);
	glVertex3f((GLfloat) (-pos), (GLfloat) (pos), -0.4);
	glVertex3f((GLfloat) (-pos), (GLfloat) (-pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (-pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (pos), -0.4);
	glEnd();

	glBegin(GL_QUADS);
	glVertex3f((GLfloat) (-pos), (GLfloat) (pos), -0.4);
	glVertex3f((GLfloat) (-pos), (GLfloat) (-pos), -0.4);
	glVertex3f((GLfloat) (-pos), (GLfloat) (-pos), 0.0);
	glVertex3f((GLfloat) (-pos), (GLfloat) (pos), 0.0);
	glEnd();
	
	glBegin(GL_QUADS);
	glVertex3f((GLfloat) (pos), (GLfloat) (-pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (pos), 0.0);
	glVertex3f((GLfloat) (pos), (GLfloat) (-pos), 0.0);
	glEnd();

		
	glBegin(GL_QUADS);
	glVertex3f((GLfloat) (-pos), (GLfloat) (pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (pos), 0.0);
	glVertex3f((GLfloat) (-pos), (GLfloat) (pos), 0.0);
	glEnd();

	
	glBegin(GL_QUADS);
	glVertex3f((GLfloat) (-pos), (GLfloat) (-pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (-pos), -0.4);
	glVertex3f((GLfloat) (pos), (GLfloat) (-pos), 0.0);
	glVertex3f((GLfloat) (-pos), (GLfloat) (-pos), 0.0);
	glEnd();
	
		
}

static void Vtrj_DrawPrimeMeridian(void)
{
	Int2 color =MAGENTA;
	GLUquadricObj *quad3;
	
	if (ShadingOn==FALSE) glColor3f(Vtrj_RGB[color][0], Vtrj_RGB[color][1], Vtrj_RGB[color][2]);
	else
	{
		glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbPurple);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifPurple);
		glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecPurple);
		glMaterialfv(GL_FRONT, GL_SHININESS, shinePurple);
    }
	
	quad3 = gluNewQuadric();
	
	gluQuadricDrawStyle(quad3, filltype);
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	gluPartialDisk(quad3,4.8,5.2,50,10,0,180);
	glRotatef(90.0, 1.0, 0.0, 0.0);
}


static void Vtrj_DrawSphere(void)
{
	Int2 color =BLUE;
	GLUquadricObj *quad2;
	
	if(ShadingOn==FALSE) glColor3f(Vtrj_RGB[color][0], Vtrj_RGB[color][1], Vtrj_RGB[color][2]);
	else
	{
	    glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbBlue);
	    glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifBlue);
	    glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecBlue);
	    glMaterialfv(GL_FRONT, GL_SHININESS, shineBlue);
	}
	
	quad2 = gluNewQuadric();
	
	gluQuadricDrawStyle(quad2, filltype);
	
	gluSphere(quad2,5, 50, 50);


}

static void Vtrj_DrawAsPlane(void)
{
	
                Vtrj_DrawBase((FloatLo) (resinfo->dim - (resinfo->dim)/2)/SCALE_SIZE_XY);
                if (Drawit==TRUE || isSwitching==TRUE || isChanging==TRUE || Vtrj_redraw==TRUE) Vtrj_Draw_Poly(GL_RENDER, FALSE);
				Vtrj_DrawAxis((FloatLo) (resinfo->dim - (resinfo->dim)/2)/SCALE_SIZE_XY);
                Vtrj_PostPicture();
}


static void Vtrj_DrawAsSphere(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	Vtrj_DrawSphere();
	Vtrj_DrawPrimeMeridian();
	Vtrj_DrawNorthPole();
	if (Drawit==TRUE || isSwitching==TRUE || isChanging==TRUE || Vtrj_redraw==TRUE) Vtrj_Draw_Poly(GL_RENDER, TRUE);
	Vtrj_PostPicture();
}

static void Vtrj_ChangeView(IteM i)
{
	isChanging=TRUE;
	Vtrj_DrawAll();
	isChanging=FALSE;

}

static void Vtrj_DrawAll(void)
{
	Int2 viewtype=2;
	static Int2 prev_view=0;
	
	viewtype = GetValue(which_view);
		
  	if (viewtype == 2)
  	{
  		if (prev_view ==1) Vtrj_LoadVertexes();
  		InitializeOGLDraw();
	  	Vtrj_DrawAsPlane();
	 }
  	if (viewtype == 1)
  	{
  		if (prev_view == 2) Vtrj_LoadVertexes();
  		InitializeOGLDraw();
  		Vtrj_DrawAsSphere();
  	}
  	
  	prev_view = viewtype;
	
}

static void Vtrj_Draw_Poly (GLenum mode, Boolean Sphere)
{

	int x;
	int max= 0 ;
	Ppp PVp = GPVp;

                	
                if (Sphere == TRUE) max = (resinfo->dim)*(resinfo->dim)-1;
                else max = ((resinfo->dim)-1)*((resinfo->dim)-1)-1;

/*	glEnable(GL_BLEND);*/
	glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
/*	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);*/
	
                for (x=0; x<=max; x++)
                {
                	
                	if (PVp->iszero==FALSE)
						if (PVp->maxpoint > rangemin && PVp->maxpoint <= rangemax)

                                {
                                	if (mode==GL_SELECT) glLoadName(x);
		
			
			if (ShadingOn==FALSE) glColor3f(Vtrj_RGB[PVp->color][0], Vtrj_RGB[PVp->color][1], Vtrj_RGB[PVp->color][2]);
			else
			{
				switch (PVp->color) {
					case RED:
						glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbRed);
						glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifRed);
						glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecRed);
					    glMaterialfv(GL_FRONT, GL_SHININESS, shineRed);
						break;
					case YELLOW:
						glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbYellow);
						glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifYellow);
						glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecYellow);
					    glMaterialfv(GL_FRONT, GL_SHININESS, shineYellow);
						break;
					case GREEN:
						glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbGreen);
						glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifGreen);
						glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecGreen);
					    glMaterialfv(GL_FRONT, GL_SHININESS, shineGreen);
						break;
					default:
						glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbBlue);
						glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifBlue);
						glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecBlue);
					    glMaterialfv(GL_FRONT, GL_SHININESS, shineBlue);
				}
            }
            
			glBegin(drawtype);
			glVertex3f(PVp->x2, PVp->y2, PVp->z2);
			glVertex3f(PVp->x4, PVp->y4, PVp->z4);
			glVertex3f(PVp->x3, PVp->y3, PVp->z3);
            glVertex3f(PVp->x1, PVp->y1, PVp->z1);
			glEnd();
		}
		PVp++;
	}
	
	
 }


/*****************************************************************
                        FILTER FUNCTIONS
*****************************************************************/

/*Gaussian Filter*/

static void Vtrj_GausFilter(void)
{
	TrajFilter filt;
	
	MakeGaussFilter(filt, gaus_mag, gaus_sd);

	TrajFilt(resinfo, filt);

	Vtrj_CalculateInfo(TRUE);
	
	Vtrj_StoreResVariables();
	SetPeaks();	

	Vtrj_LoadVertexes();
	RangeProc();
	Vtrj_redraw=TRUE;
	Vtrj_DrawAll();
	Vtrj_redraw=FALSE;
	/*Vtrj_ClearXYZ();	*/
	HasBeenFiltered = TRUE;
	EnableDisableProc();
}


static void DoneGausAttr(IteM i)
{
	Hide(GausAttr_win);
 }


static void GetGausValues(IteM i)
{

    float fmag_num, fsd_num;
    int err;
	Char mag[10], sd[10];
	
		GetTitle (gaus_mag_text, mag, sizeof(mag));
		GetTitle (gaus_sd_text, sd, sizeof(sd));
	
		err=sscanf(mag, "%f", &fmag_num);
		if (err==0)
			fmag_num=0.0;
		err=sscanf(sd, "%f", &fsd_num);
		if (err==0)
			fsd_num=0.0;
		if (fmag_num<10.0)
		{
			Message(MSG_OK, "Magnitude must be 10 or greater!");
			SetTitle(gaus_mag_text, "");
			Select(gaus_mag_text);
			 return;
		}
		else if (fsd_num<0.5)
		{
			Message(MSG_OK, "Standard deviation must be 0.5 or greater!");
			SetTitle(gaus_sd_text, "");
			Select(gaus_sd_text);
			 return;
		}
	
	 gaus_mag = (FloatLo) fmag_num;
	 gaus_sd = (FloatLo) fsd_num;
	
	 DoneGausAttr(i);	
	
	 Vtrj_GausFilter();
	
}



static void ShowGausAttr(IteM i)
{

	GrouP g, Vtrj_grp2, Vtrj_grp3;
	
	GausAttr_win = MovableModalWindow(-50,-33,-10,-10, "Enter Gaussian Filter Parameters",  (WndActnProc) DoneGausAttr);
	
	g = HiddenGroup(GausAttr_win, -1, -1, (GrpActnProc) NULL);

	SetGroupSpacing(g, 3, 15);
	SetGroupMargins(g, 3, 3);

	Vtrj_grp2 = HiddenGroup(g, 2, 0, (GrpActnProc) NULL);

    StaticPrompt (Vtrj_grp2, "Magnitude", 0, dialogTextHeight, systemFont, 'l');
        	
    gaus_mag_text = DialogText (Vtrj_grp2, "100.00",  10, (TxtActnProc) NULL);
	
	StaticPrompt (Vtrj_grp2, "Standard Deviation", 0, dialogTextHeight, systemFont, 'l');
        	
    gaus_sd_text = DialogText(Vtrj_grp2, "1.00", 5, (TxtActnProc) NULL);
	
	Vtrj_grp3 = HiddenGroup(g, 2, -1, (GrpActnProc) NULL);
	
	DefaultButton(Vtrj_grp3, "OK", (BtnActnProc) GetGausValues);
	PushButton(Vtrj_grp3, "Cancel", (BtnActnProc) DoneGausAttr);
               	
	AlignObjects(ALIGN_CENTER, Vtrj_grp3, g, NULL);

    Show(GausAttr_win);
	
	ProcessEvents();
}


				/*No Filter*/

				
static void Vtrj_NoneFilter(IteM i)
{

	
	if (isCis==FALSE)
		TrajCopy(origresinfo->TrajGraph, resinfo->TrajGraph, resinfo->dim);
	
	if (isCis==TRUE)
		TrajCopy(origresinfo->CisTrajGraph, resinfo->CisTrajGraph, resinfo->dim);

	resinfo->Peak = origresinfo->Peak;
	resinfo->CisPeak = origresinfo->CisPeak;
		
	Vtrj_CalculateInfo(TRUE);
	SetPeaks();
			
	Vtrj_LoadVertexes();
	RangeProc();
	Vtrj_redraw=TRUE;
    Vtrj_DrawAll();
	Vtrj_redraw=FALSE;
    /*Vtrj_ClearXYZ();*/
    HasBeenFiltered = FALSE;
	EnableDisableProc();


}


				/*Smooth Filter*/
				
static void Vtrj_SmoothFilter(IteM i)
{

	TrajFilter filt;
		  	
                MakeSmoothFilter(filt);

	TrajFilt(resinfo, filt);

	Vtrj_CalculateInfo(TRUE);
	Vtrj_StoreResVariables();
	
	SetPeaks();
		
	Vtrj_LoadVertexes();
	RangeProc();
	Vtrj_redraw=TRUE;
    Vtrj_DrawAll();
	Vtrj_redraw=FALSE;
    /*Vtrj_ClearXYZ();	*/
    HasBeenFiltered = TRUE;
	EnableDisableProc();

}


				/*Low Pass Filter*/

				
static void Vtrj_LPFilter(IteM i)
{
	TrajFilter filt;
	
    MakeLPFilter(filt, lp_mag, lp_sd);

	TrajFilt(resinfo, filt);
	Vtrj_CalculateInfo(TRUE);
	Vtrj_StoreResVariables();
	
    SetPeaks();
	
	Vtrj_LoadVertexes();
	RangeProc();
	Vtrj_redraw=TRUE;
    Vtrj_DrawAll();
	Vtrj_redraw=FALSE;
    /*Vtrj_ClearXYZ();	*/
    HasBeenFiltered = TRUE;
	EnableDisableProc();
        	
}


static void DoneLPAttr(IteM i)
{
	Hide(LPAttr_win);
 }


static void GetLPValues(IteM i)
{
	int err;
 	float fmag_num, fsd_num;
	Char mag[10], sd[10];
	
	GetTitle (lp_mag_text, mag, sizeof(mag));
	GetTitle (lp_sd_text, sd, sizeof(sd));
	
	 err=sscanf(mag, "%f", &fmag_num);
	 if (err==0)
	 	fmag_num=0.0;
	 err=sscanf(sd, "%f", &fsd_num);
	 if (err==0)
	 	fsd_num=0.0;
	
	if (fmag_num<10.0)
	{
		Message(MSG_OK, "Magnitude must be 10 or greater!");
		SetTitle(lp_mag_text, "");
		Select(lp_mag_text);
		return;
	}
	else if (fsd_num<0.5)
	{
		Message(MSG_OK, "Standard deviation must be 0.5 or greater!");
		SetTitle(lp_sd_text, "");
		Select(lp_mag_text);
		return;
	}
	lp_mag = (FloatLo) fmag_num;
	lp_sd = (FloatLo) fsd_num;
	
	DoneLPAttr(i);	
	
	Vtrj_LPFilter(i);
	
}



static void ShowLPAttr(IteM i)
{
	
	GrouP g, Vtrj_grp1, Vtrj_grp2;

	LPAttr_win = MovableModalWindow(-50,-33,-10,-10, "Enter Low Pass Filter Parameters", (WndActnProc) DoneLPAttr);
	
    g = HiddenGroup(LPAttr_win, -1, -1, (GrpActnProc) NULL);

	SetGroupSpacing(g, 3, 15);
	SetGroupMargins(g, 3, 3);

	Vtrj_grp1 = HiddenGroup(g, 2, 0, (GrpActnProc) NULL);
        	
	StaticPrompt (Vtrj_grp1, "Magnitude", 0, dialogTextHeight, systemFont, 'l');
	
	lp_mag_text = DialogText (Vtrj_grp1, "100.00",  10, (TxtActnProc) NULL);
	
	StaticPrompt (Vtrj_grp1, "Standard Deviation", 0, dialogTextHeight, systemFont, 'l');
        	
    lp_sd_text = DialogText(Vtrj_grp1, "1.50", 5, (TxtActnProc) NULL);
	
	Vtrj_grp2 = HiddenGroup(g,2, -1, (GrpActnProc) NULL);

	DefaultButton(Vtrj_grp2, "OK", (BtnActnProc) GetLPValues);

	PushButton(Vtrj_grp2, "Cancel", (BtnActnProc) DoneLPAttr);
               	
	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);
    
			
	Show(LPAttr_win);
	
	ProcessEvents();
}


/*****************************************************************
                                               TEXT FUNCTIONS
*****************************************************************/


void stroke_output(GLfloat size, char *format,...)
{
	va_list args;
	char buffer[200], *p;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
	glPushMatrix();
	glScalef(size, size, size);
	for (p = buffer; *p; p++)
	glutStrokeCharacter(GLUT_STROKE_ROMAN, *p);
	glPopMatrix();
}


void Vtrj_PutText (CharPtr string, int color, GLfloat x, GLfloat y, GLfloat size, GLfloat width, GLfloat rotz)
{
	glEnable(GL_LINE_SMOOTH);
  	glEnable(GL_BLEND);
	glLineWidth(width);	
	
	if (ShadingOn==FALSE) glColor3f(Vtrj_RGB[color][0], Vtrj_RGB[color][1], Vtrj_RGB[color][2]);
	else
	{
		if (color==BLACK) {
			glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbBlack);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifBlack);
			glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecBlack);
			glMaterialfv(GL_FRONT, GL_SHININESS, shineBlack);
		}		
		else { /* color==WHITE */
			glMaterialfv(GL_FRONT, GL_AMBIENT, matColorAmbWhite);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, matColorDifWhite);
			glMaterialfv(GL_FRONT, GL_SPECULAR, matColorSpecWhite);
			glMaterialfv(GL_FRONT, GL_SHININESS, shineWhite);
		}		
    }
			
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-Vtrj_scale, (Vtrj_scale*Vtrj_Xratio), (-Vtrj_scale*Vtrj_Yratio), Vtrj_scale, -50.0, 50.0);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	
	glPushMatrix();
	glTranslatef(x, y, 0.0);
	glRotatef(rotz, 0.00, 0.00, 1.00);
	stroke_output(size, string);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);
	
}

/*****************************************************************
                                               DRAW CHOICE FUNCTIONS
*****************************************************************/

static void QuickRotProc(IteM i)
{
    Int2 rot_choice;
	rot_choice = GetValue(which_rot);
	
	switch (rot_choice)
	{
		case 1:
			QuickRotate = TRUE;
			Drawit = FALSE;
			break;
		case 2:
			QuickRotate = FALSE;
			Drawit = TRUE;
			break;
	}
}

static void ShadeProc(IteM i)
{
        	Int2 shade_choice;
	shade_choice = GetValue(which_shade);
	
	switch (shade_choice)
	{
		case 1:
			ShadingOn = TRUE;
			break;
		case 2:
			ShadingOn = FALSE;
			break;
	}

	
        	Vtrj_DrawAll();
/*        	Vtrj_ClearXYZ();*/

}

static void Vtrj_TransCisProc(IteM i)
{
        	Int2 trans_choice;
	trans_choice = GetValue(which_trans);
	
	switch (trans_choice)
	{
		case 1:
			isCis = FALSE;
			break;
		case 2:
			isCis = TRUE;
			break;
	}

	SetPeaks();
	RangeProc();
	
	Vtrj_LoadVertexes();
    Vtrj_DrawAll();
    /*Vtrj_ClearXYZ();*/
}



static void DrawProc(IteM i)
{
	Int2 draw_choice;
	draw_choice = GetValue(which_draw);
	
	switch (draw_choice)
	{
		case 1:
			drawtype = GL_QUADS;
			filltype = GLU_FILL;
			Vtrj_DrawAll();
			break;
		case 2:
			drawtype = GL_LINE_LOOP;
			filltype = GLU_LINE;
			Vtrj_DrawAll();			
			break;
	}
	
}

static void VisHeightProc(IteM i)
{
	RangeProc();
	Vtrj_redraw=TRUE;
	Vtrj_DrawAll();
	Vtrj_redraw=FALSE;
}

static void RangeProc(void)
{
	Int2 range_choice;

	range_choice = GetValue(which_range);
	
	switch (range_choice)
	{
		case 1:
			/* set to -1 so 0 > rangemin... */
			rangemin = -1;
			rangemax = ThePeak;
			break;
		case 2:
			rangemin = -1;
			rangemax = ThePeak*0.5;
			break;
		case 3:
			rangemin = -1;
			rangemax = ThePeak*0.25; 	
			break;
		case 4:	
			rangemin = ThePeak*0.50;
			rangemax = ThePeak;
			break;
		case 5:
			rangemin = ThePeak*0.25;
			rangemax = ThePeak;
			break;
		case 6:
			rangemin = ThePeak*0.25;
			rangemax = ThePeak*0.50;
			break;
	}			
	
}



static void TopViewProc (IteM i)
{
	
	Vtrj_scale = 10.0;
	Vtrj_rotateX = 90;
	Vtrj_rotateY = 0.0;
	Vtrj_rotateZ = 0.0;
	
	VisHeightProc(i);

}

static void BottomViewProc (IteM i)
{
		
	Vtrj_scale = 10.0;
	Vtrj_rotateX = 90.0;
	Vtrj_rotateY = 180.0;
	Vtrj_rotateZ = 0.0;
	
	VisHeightProc(i);
}


static void XViewProc (IteM i)
{
	
	Vtrj_scale = 10.0;
	Vtrj_rotateX = 90;
	Vtrj_rotateY = -90;
	Vtrj_rotateZ = 0.0;
	
	VisHeightProc(i);

}


static void YViewProc (IteM i)
{
	
	Vtrj_scale = 10.0;
	Vtrj_rotateX = -90;
	Vtrj_rotateY = -90;
	Vtrj_rotateZ = 0.0;
	
	VisHeightProc(i);

}


static void NormViewProc (IteM i)
{

	Vtrj_scale = 10.0;
	Vtrj_rotateX = 90;
	Vtrj_rotateY = -35;
	Vtrj_rotateZ = 0;
	
	VisHeightProc(i);

}



/*****************************************************************
                                               HELP WINDOW FUNCTIONS
*****************************************************************/


static void LoadHelpMenu(void)
{
	Char menuchoice[MAXCOL+1];
	Char inbuf[MAXCOL];
	Char fnam[PATH_MAX];
	int x, j=1;
	FILE *f;
	
	inbuf[0] = '\0';
	menuchoice[0] = '\0';
	
	sprintf(fnam,"%s%s",CFG_local_datafilepath,VTRJ_HELP_FILE);
	if ((f=FileOpen(fnam,"r"))==NULL)
		ErrPostEx(SEV_ERROR,1,1,"Unable to open VisTraj Help file %s", fnam);
	
	while (FileGets(inbuf, MAXCOL,f)!=NULL)
	{
		if (inbuf[0] == '!')
		{
			for (x=0; x<StringLen(inbuf); x++)
				menuchoice[x] = inbuf[x+1];
			
			menuchoice[x-2] = '\0';
			StringCpy(Vtrj_Help_Menu[j], menuchoice);
	 	  	ListItem (helpmenulist, menuchoice);
			j++;
		}
	}
	FileClose(f);
		
}

static void DoneHelpProc(IteM i)
{
	Hide(help_win);	
	inHelp=FALSE;
}

static void LoadHelpTopicProc(IteM i)
{

	Char printedstring[MAXCOL+1];
    Char inbuf[MAXCOL];
    Char menuchoice[MAXCOL+1];
	Char fnam[PATH_MAX];
	Int2 loadwhat = 0, x;
	FILE *f;
	Boolean shouldloadit=FALSE;
	
	loadwhat = GetValue(helpmenulist);
	Reset(help_doc);
	
	hdgColFmt.font = ParseFont ("Times,18,b");
	subColFmt.font = ParseFont ("Helvetica,12,b");
#ifdef WIN_MSWIN
	txtColFmt.font = ParseFont ("Courier,12");
#else
	txtColFmt.font = ParseFont ("fixed,12");
#endif
	lstColFmt.font = programFont;
	tblColFmt.font = programFont;

	sprintf(fnam,"%s%s",CFG_local_datafilepath,VTRJ_HELP_FILE);
	if ((f=FileOpen(fnam,"r"))==NULL)
		ErrPostEx(SEV_ERROR,1,1,"Unable to open VisTraj Help file %s", fnam);
	
	while (FileGets(inbuf, MAXCOL,f)!=NULL)
	{
		while ((inbuf[StringLen(inbuf)-1]=='\n') || (inbuf[StringLen(inbuf)-1]=='\r'))
			inbuf[StringLen(inbuf)-1]='\0';
		if (shouldloadit==TRUE)
		{
			if (inbuf[0] =='#')
			{
				for (x=0; x<StringLen(inbuf); x++) printedstring[x] = inbuf[x+1];
 			  	AppendText (help_doc, printedstring, &hdgParFmt, &subColFmt, programFont);
			}
			else if(inbuf[0]!='!') {
				AppendText (help_doc, inbuf, &txtParFmt, &txtColFmt, programFont);
			}
			else {
				FileClose(f);
				UpdateDocument (help_doc, 0, 0);
				Show(help_doc);
				return;
			}
		}
		
		else if (inbuf[0] == '!')
		{
			for (x=0; x<StringLen(inbuf); x++)
			   menuchoice[x] = inbuf[x+1];
			
			if (StringCmp(menuchoice,Vtrj_Help_Menu[loadwhat])==0)
			{
				
				AppendText (help_doc, menuchoice, &hdgParFmt, &hdgColFmt, programFont);
				shouldloadit = TRUE;
			}
		}
	}
	UpdateDocument (help_doc, 0, 0);
	FileClose(f);

}


static void LoadHelpProc (IteM i)
{
	GrouP g, Vtrj_grp1, Vtrj_grp2;

	if (inHelp) return;
	inHelp=TRUE;
	help_win = FixedWindow(-50,-33,-30,-30, "VisTraj Help", (WndActnProc) DoneHelpProc);

	g = HiddenGroup(help_win,-1,-1, (GrpActnProc) NULL);
	Vtrj_grp1 =  HiddenGroup(g,2,0,(GrpActnProc) NULL);
	
	SetGroupMargins(g, 3, 3);
	SetGroupSpacing(g, 3, 3);			
	
	helpmenulist = SingleList (Vtrj_grp1,18,10, (LstActnProc) LoadHelpTopicProc);
#ifdef WIN_MSWIN
	help_doc = DocumentPanel (Vtrj_grp1, stdCharWidth * 46, stdLineHeight * 20);
#else
	help_doc = DocumentPanel (Vtrj_grp1, stdCharWidth * 33, stdLineHeight * 20);
#endif
	Vtrj_grp2 = HiddenGroup(g,-1,-1,(GrpActnProc) NULL);
	
	DefaultButton(Vtrj_grp2, "     Exit     ", (BtnActnProc) DoneHelpProc);

	AlignObjects(ALIGN_CENTER, Vtrj_grp2, g, NULL);

	LoadHelpMenu();	
		
	Show(help_win);
	
	ProcessEvents();
}


/*****************************************************************
                                               MAIN WINDOW FUNCTIONS
*****************************************************************/


static Boolean CheckIfSaved (void)
{
	Int2 Vtrj_ans = ANS_NO;
	Boolean is_OK=TRUE;
	
	if (Rhasbeenedited==TRUE || Ghasbeenedited==TRUE || Nhasbeenedited==TRUE || HasBeenFiltered==TRUE || DistBeenEdited==TRUE || FragBeenEdited==TRUE)
		Vtrj_ans = Message (MSG_YNC, "You have unsaved changes.  Would you like to save them now?");	
	
	if (Vtrj_ans==ANS_CANCEL) return FALSE;
	
	if (Vtrj_ans==ANS_YES) 
	{
		ShouldSaveAs = FALSE;
		is_OK = SaveTrajProc();
	}
	
	if (Vtrj_ans==ANS_NO)
	{
		Rhasbeenedited=FALSE;
		Ghasbeenedited=FALSE;
		Nhasbeenedited=FALSE;
		HasBeenFiltered=FALSE;
		DistBeenEdited=FALSE;
		FragBeenEdited=FALSE;
	}
	
	return is_OK;
}

static void CleanUpTrajProc(void)
{
	
	SetSwitchParams (resid_switch, 0, 0);
	SetTitle(resnum_title,"");
	Vtrj_UnLoadSequence();
	pointedAA[0] = '\0';
	DrawAANameProc(aaname_panel);
	Disable(jumpto_text);
	SetTitle(main_win, prog_name);
	/*Vtrj_ClearXYZ();*/
	Hide(modAA_temptitle);
	Hide(modAA_title);
	Vtrj_ClearScreen();
	isLoaded=FALSE;
	if (pbi!=NULL)
		AsnGenericChoiceSeqOfFree(pbi,(AsnOptFreeFunc)BiostrucIdFree);
	if (pbd!=NULL)
		AsnGenericChoiceSeqOfFree(pbd,(AsnOptFreeFunc)BiostrucDescrFree);
	if (vnpBioseq!=NULL)
		AsnGenericChoiceSeqOfFree(vnpBioseq,(AsnOptFreeFunc)SeqEntryFree);
}

static Boolean QuitTrajProc (void)
{
	Boolean is_OK=FALSE;

	is_OK = CheckIfSaved();	
	
	if (is_OK==FALSE) return is_OK;

	if (isLoaded == TRUE)
	{
		CleanUpDB(tmpdbasename);
	 	FreePNNList();
		isLoaded=FALSE;
	}

	return is_OK;
}


static void ExitProc (IteM i)
{
	Boolean is_OK=FALSE;
	
  if (alldisabled==TRUE)
  	return;
	if (ProgramProgress>0)
		return;
	is_OK = QuitTrajProc();
	if(is_OK==FALSE) return;

	if (hasafrag!=NULL)	hasafrag=MemFree(hasafrag);
	if (resinfo != NULL) FreeTraj(resinfo);
	if (origresinfo != NULL) FreeTraj(origresinfo);
	if (reseditinfo != NULL) reseditinfo = MemFree(reseditinfo);
	if (OrigGlobs!=NULL) OrigGlobs = MemFree(OrigGlobs);
	if (GPVp!=NULL) GPVp = MemFree(GPVp);
	if (NmodAAlist!=NULL) ValNodeFreeData(NmodAAlist);
	if (CmodAAlist!=NULL) ValNodeFreeData(CmodAAlist);
	if (MmodAAlist!=NULL) ValNodeFreeData(MmodAAlist);
	if (fnamtrj!=NULL) fnamtrj=MemFree(fnamtrj);
 	FreeRotLib();
	CLOSEMMDBAPI();
	QuitProgram ();

}

static void CloseTrajProc (IteM i)
{
	Boolean is_OK=FALSE;

	
	is_OK = QuitTrajProc();
	if(is_OK==FALSE) return;

	CleanUpTrajProc();

	Hide(aaname_panel);
	Show(aaname_panel);

	EnableDisableProc();
}


static void AboutProc (IteM i)
{
	Char buf[1024];

	if (!StringCmp(FOLDTRAJ_VERSION," Unknown"))
		sprintf(buf,"OpenGL Trajectory Distribution Visualizer V%s\n\nAUTHORS:  John J. Salama\n                     Howard J. Feldman\n                    Christopher W.V. Hogue\n\nWEB SITE: http://www.blueprint.org/\nE-mail bug reports, praise and suggestions to: \hogue@nus.edu.sg",VISTRAJ_VERSION);
	else
		sprintf(buf,"OpenGL Trajectory Distribution Visualizer V%s\nBuilt on %s\n\nAUTHORS:  John J. Salama\n                    Howard J. Feldman\n                    Christopher W.V. Hogue\n\nWEB SITE: http://www.blueprint.org/\nE-mail bug reports, praise and suggestions to:\hogue@nus.edu.sg",VISTRAJ_VERSION,FOLDTRAJ_VERSION);
	Message (MSG_OK,buf); 
}

static void Vtrj_GetWinSize(void)
{

	FloatLo width=0.0, height=0.0;

	width = (FloatLo) screenRect.right-screenRect.left;
	height = (FloatLo) screenRect.bottom-screenRect.top;

	if ((width/height) > 1.0)
	{
		Vtrj_Xratio = (GLfloat) width/height;
		Vtrj_Yratio = (GLfloat) 1.0;	
	}

	else
	{
		Vtrj_Xratio = (GLfloat) 1.0;
		Vtrj_Yratio = (GLfloat) height/width;
	}
}

void Vtrj_ResizeProc(WindoW w)
{

	Cn3DResizeProc(w);
	Vtrj_GetWinSize();
	if (isLoaded == TRUE) Vtrj_DrawAll();
}


static void ResizeProc(IteM i)
{
	if (ProgramProgress>0)
		return;
/* Dummy Function in order to integrate with NCBI OGL resize procedure*/
	Vtrj_ResizeProc(main_win);
}



static void SetupWindow (void)
{
	GrouP main_group, scnd_group;
		
	numAA = 1;
	res_num = 1;
	
	main_group = HiddenGroup (main_win,2,0, (GrpActnProc) NULL);
	scnd_group = HiddenGroup (main_win, 4,0, (GrpActnProc) NULL);
	
	SetGroupMargins(scnd_group, 2, 15);
	SetGroupSpacing(scnd_group, 10, 0);
			
	resid_switch = LeftRightSwitch (main_group,0, (SwtChngProc) PreSwitchResProc);
        	SetSwitchParams (resid_switch, 1, numAA);

	Advance (main_win);
	
        	StaticPrompt (main_win, "Residue #: ", 0, dialogTextHeight, systemFont, 'l');

        	Advance (main_win);

        	resnum_title = StaticPrompt (main_win," ", 50, dialogTextHeight, systemFont, 'l');
     	
        	Advance (main_win);

        	StaticPrompt (main_win,"           Amino Acid Sequence: ",0,dialogTextHeight, systemFont,'l');

        	Advance(main_win);
	aaname_panel = SimplePanel (main_win, 100, 20, DrawAANameProc);        	
        	
       	Advance(main_win);
		
/*	StaticPrompt(main_win, "     X=",0,dialogTextHeight, programFont, 'r');
		
        	Advance(main_win);
        	
        	xval_title = StaticPrompt (main_win, "", 30, dialogTextHeight, programFont, 'r');
        	
        	Advance(main_win);
	
	StaticPrompt(main_win, "Y=",0,dialogTextHeight, programFont, 'r');
	
        	Advance(main_win);
        	
        	yval_title = StaticPrompt (main_win, "", 30, dialogTextHeight, programFont, 'r');
        	
        	Advance(main_win);
	
	StaticPrompt(main_win, "Z=",0,dialogTextHeight, programFont, 'r');

	Advance(main_win);
        	
        	zval_title = StaticPrompt (main_win, "", 30, dialogTextHeight, programFont, 'r');
        	
        	Advance(main_win);
*/     	       	
        	Break (main_win);

        	StaticPrompt (scnd_group, "Jump to Residue:", 0, dialogTextHeight, systemFont, 'l');

        	Advance (main_win);

       	jumpto_text = DialogText (scnd_group,"", 5,(TxtActnProc) TextEnteredProc);
		Disable(jumpto_text);
        Select(jumpto_text);

        	Advance (main_win);

	go_bttn = DefaultButton (scnd_group, "Go", (BtnActnProc) GoBttnPressedProc);

	Advance (main_win);
	
	Disable (go_bttn);
	
	Advance (main_win);
	
	modAA_temptitle = StaticPrompt(main_win, "                     Modified A.A.:", 0, dialogTextHeight, systemFont, 'r');
	
	Hide(modAA_temptitle);
	
	Advance (main_win);
	
	modAA_title = StaticPrompt(main_win, "", 100, dialogTextHeight, programFont, 'l');
	
	Break(main_win);
}



static void SetupMenus (void)
{
	
MenU Vtrj_file_menu, Vtrj_edit_menu, Vtrj_view_menu, Vtrj_help_menu,
		 Vtrj_filter_menu, Vtrj_draw_menu, Vtrj_type_menu,
		 Vtrj_range_menu, Vtrj_pos_menu, Vtrj_shade_menu,
		 Vtrj_quikrot_menu, Vtrj_generate_menu;
	

	Vtrj_file_menu = PulldownMenu (main_win, "File/F");
    Vtrj_edit_menu = PulldownMenu (main_win, "Edit/E");
    Vtrj_filter_menu = PulldownMenu (main_win, "Filter/T");
	Vtrj_view_menu = PulldownMenu (main_win, "View/V");
	Vtrj_draw_menu = PulldownMenu (main_win, "Draw/D");
    Vtrj_help_menu = PulldownMenu (main_win, "Help/H");
        	
    Vtrj_type_menu = SubMenu (Vtrj_draw_menu, "Rendering");
    Vtrj_shade_menu = SubMenu(Vtrj_draw_menu, "Shading");
    Vtrj_range_menu = SubMenu (Vtrj_draw_menu, "Visible Height");
    Vtrj_pos_menu = SubMenu (Vtrj_draw_menu, "Camera Position");
    Vtrj_quikrot_menu = SubMenu(Vtrj_draw_menu, "Quick Rotation");

	Vtrj_New = CommandItem (Vtrj_file_menu, "New.../N", CreateNewTrajProc);
	SeparatorItem (Vtrj_file_menu);
    Vtrj_Load = CommandItem (Vtrj_file_menu, "Open.../O", GetNewFileProc);
    Vtrj_Save = CommandItem (Vtrj_file_menu, "Save/S", Vtrj_SaveProc);
	Vtrj_SaveAs = CommandItem (Vtrj_file_menu, "Save As.../A", Vtrj_SaveAsProc);
	Vtrj_Close = CommandItem (Vtrj_file_menu, "Close/C", CloseTrajProc);
	SeparatorItem (Vtrj_file_menu);
    /*Vtrj_SavePNG = CommandItem (Vtrj_file_menu, "Save as PNG", SavePNGProc);*/
    Vtrj_generate_menu = SubMenu(Vtrj_file_menu, "Generate Conformer");
	Vtrj_Protein = CommandItem (Vtrj_file_menu, "Structure Info", ProteinInfoProc);
	SeparatorItem (Vtrj_file_menu);
	Vtrj_Exit = CommandItem (Vtrj_file_menu, "Exit", ExitProc);

	
	Vtrj_Global = CommandItem (Vtrj_edit_menu, "Global Properties", EditGInfoProc);
	Vtrj_Residue = CommandItem (Vtrj_edit_menu, "Residue Properties", EditRInfoProc);
	
	SeparatorItem(Vtrj_edit_menu);
	
	Vtrj_Dist_Cons = CommandItem (Vtrj_edit_menu, "Distance Contraints", DisplayDistConsProc);
	Vtrj_Frag = CommandItem (Vtrj_edit_menu, "Residue Fragments", DisplayFragListProc);
	Vtrj_ImportFrag = CommandItem (Vtrj_edit_menu, "Import Fragments from File", ImportFragListProc);

	SeparatorItem(Vtrj_edit_menu);

	Vtrj_Noise = CommandItem (Vtrj_edit_menu, "Add\\Replace Noise", AddNoiseProc);

	SeparatorItem(Vtrj_edit_menu);

	Vtrj_AddRes = CommandItem (Vtrj_edit_menu, "Insert Residue", AddResidueProc);
	Vtrj_RemoveRes = CommandItem (Vtrj_edit_menu, "Remove Residue", RemoveResidueProc);
	Vtrj_Mutate = CommandItem (Vtrj_edit_menu, "Mutate Residue", LoadModifyRProc);

	Vtrj_None = CommandItem (Vtrj_filter_menu, "Revert to Original", Vtrj_NoneFilter);
	Vtrj_Smooth = CommandItem (Vtrj_filter_menu, "Smooth", Vtrj_SmoothFilter);
	Vtrj_Gaussian = CommandItem (Vtrj_filter_menu, "Gaussian", ShowGausAttr);
	Vtrj_LP = CommandItem (Vtrj_filter_menu, "Low Pass", ShowLPAttr);

	which_view = ChoiceGroup(Vtrj_view_menu, (ChsActnProc) Vtrj_ChangeView);
	Vtrj_Sphere = ChoiceItem (which_view, "Sphere");
	Vtrj_Plane = ChoiceItem (which_view, "Plane");
	SetValue (which_view, 2);

	SeparatorItem(Vtrj_view_menu);

	which_trans = ChoiceGroup(Vtrj_view_menu, (ChsActnProc) Vtrj_TransCisProc);
	Vtrj_Trans = ChoiceItem (which_trans, "Trans");
	Vtrj_Cis = ChoiceItem (which_trans, "Cis");
	SetValue (which_trans, 1);
	
	which_draw = ChoiceGroup(Vtrj_type_menu, (ChsActnProc) DrawProc);
	Vtrj_Filled = ChoiceItem (which_draw, "Solid");
	Vtrj_Wire = ChoiceItem (which_draw, "Wire Frame");
	SetValue (which_draw, 1);
	
	which_shade = ChoiceGroup(Vtrj_shade_menu, (ChsActnProc) ShadeProc);
	Vtrj_Shade = ChoiceItem (which_shade, "On");
	Vtrj_NoShade = ChoiceItem (which_shade, "Off");
	SetValue (which_shade, 2);
	
	which_range = ChoiceGroup(Vtrj_range_menu, (ChsActnProc) VisHeightProc);
	Vtrj_Range1 = ChoiceItem (which_range, "0 - 100%");
	Vtrj_Range2 = ChoiceItem (which_range, "0 - 50%");
	Vtrj_Range3 = ChoiceItem (which_range, "0 - 25%");
	Vtrj_Range4 = ChoiceItem (which_range, "50 - 100%");
	Vtrj_Range5 = ChoiceItem (which_range, "25 - 100%");
	Vtrj_Range6 = ChoiceItem (which_range, "25 - 50%");
	SetValue(which_range, 1);
	
	Vtrj_PosNorm = CommandItem (Vtrj_pos_menu, "Default", NormViewProc);
	Vtrj_PosXSide = CommandItem (Vtrj_pos_menu, "X Axis Side" , XViewProc);
	Vtrj_PosYSide = CommandItem (Vtrj_pos_menu, "Y Axis Side", YViewProc);
	Vtrj_PosTop = CommandItem (Vtrj_pos_menu, "Top",TopViewProc);
	Vtrj_PosBottom = CommandItem (Vtrj_pos_menu, "Bottom", BottomViewProc);
	
	Vtrj_GenerateVal = CommandItem(Vtrj_generate_menu, "PDB (Launch Rasmol)",GeneratePdbProc);
	Vtrj_GeneratePdb = CommandItem(Vtrj_generate_menu, "MMDB (Launch Cn3D)",GenerateValProc);

	which_rot = ChoiceGroup(Vtrj_quikrot_menu, (ChsActnProc) QuickRotProc);
	Vtrj_QuickOn = ChoiceItem(which_rot, "On");
	Vtrj_QuickOff = ChoiceItem(which_rot, "Off");
	SetValue(which_rot, 2);
	
	SeparatorItem(Vtrj_draw_menu);
	Vtrj_BgColor=CommandItem(Vtrj_draw_menu, "Background Color", DisplayBackGroundProc);
	/*Vtrj_Multiple = CommandItem (Vtrj_draw_menu, "Multiple Residues", MultipleViewProc);*/
	                              	
	Vtrj_Help = CommandItem (Vtrj_help_menu, "Contents", LoadHelpProc);
	Vtrj_About = CommandItem (Vtrj_help_menu, "About", AboutProc);
	
	/* SeparatorItem(Vtrj_help_menu);
	Vtrj_WebSite=CommandItem (Vtrj_help_menu, "TRADES Website", LaunchTradesSite);	*/

}	


static void EnableDisableProc(void)
{
	
		
	if (isLoaded ==FALSE)
	{
		Disable(Vtrj_Global);
		Disable(Vtrj_Residue);
		Disable(Vtrj_None);
		Disable(Vtrj_Smooth);
		Disable(Vtrj_Gaussian);
		Disable(Vtrj_LP);
		Disable(Vtrj_Sphere);
		Disable(Vtrj_Plane);
		Disable(Vtrj_Filled);
		Disable(Vtrj_Wire);
		Disable(Vtrj_Range1);
		Disable(Vtrj_Range2);
		Disable(Vtrj_Range3);
		Disable(Vtrj_Range4);
		Disable(Vtrj_Range5);
		Disable(Vtrj_Range6);
		Disable(Vtrj_GeneratePdb);
		Disable(Vtrj_GenerateVal);
		Disable(Vtrj_Protein);
		Disable(Vtrj_PosXSide);
		Disable(Vtrj_PosYSide);
		Disable(Vtrj_PosTop);
		Disable(Vtrj_PosBottom);
		Disable(Vtrj_PosNorm);
		Disable(Vtrj_Trans);
		Disable(Vtrj_Cis);
		Disable(Vtrj_SavePNG);
		Disable(Vtrj_Save);
		Disable(Vtrj_Noise);
		Disable(Vtrj_AddRes);
		Disable(Vtrj_RemoveRes);
		Disable(Vtrj_Shade);
		Disable(Vtrj_NoShade);
		Disable(Vtrj_QuickOn);
		Disable(Vtrj_QuickOff);
		/*Disable(Vtrj_Multiple);*/
		Disable(Vtrj_Dist_Cons);
		Disable(Vtrj_Frag);
		Disable(Vtrj_ImportFrag);
		Disable(Vtrj_Mutate);
		Disable(Vtrj_Close);
		Disable(Vtrj_SaveAs);

		if (isNew==TRUE) {
			Disable(Vtrj_Load);
			Disable(Vtrj_New);
		}
		if (isNew==FALSE) {
			Enable(Vtrj_Load);
			Enable(Vtrj_New);
		}

	}
	
                
	if (isLoaded ==TRUE)
	{
		
		Enable(Vtrj_Global);
		Enable(Vtrj_Residue);
		Enable(Vtrj_None);
		Enable(Vtrj_Smooth);
		Enable(Vtrj_Gaussian);
		Enable(Vtrj_LP);
		Enable(Vtrj_Sphere);
		Enable(Vtrj_Plane);
		Enable(Vtrj_Filled);
		Enable(Vtrj_Wire);
		Enable(Vtrj_Range1);
		Enable(Vtrj_Range2);
		Enable(Vtrj_Range3);
		Enable(Vtrj_Range4);
		Enable(Vtrj_Range5);
		Enable(Vtrj_Range6);
		Enable(Vtrj_GeneratePdb);
		Enable(Vtrj_GenerateVal);
		Enable(Vtrj_Protein);
		Enable(Vtrj_PosXSide);
		Enable(Vtrj_PosYSide);
		Enable(Vtrj_PosTop);
		Enable(Vtrj_PosBottom);
		Enable(Vtrj_PosNorm);
		Enable(Vtrj_Trans);
		Enable(Vtrj_SavePNG);
		Enable(Vtrj_Save);
		Enable(Vtrj_Noise);
		Enable(Vtrj_AddRes);
		Enable(Vtrj_RemoveRes);
		Enable(Vtrj_Shade);
		Enable(Vtrj_NoShade);
		Enable(Vtrj_QuickOn);
		Enable(Vtrj_QuickOff);
		/*Enable(Vtrj_Multiple);*/
		Enable(Vtrj_Dist_Cons);
		Enable(Vtrj_Frag);
		Enable(Vtrj_ImportFrag);
		Enable(Vtrj_Mutate);
		Enable(Vtrj_SaveAs);
		Enable(Vtrj_Close);
		Enable(Vtrj_Load);
		Enable(Vtrj_New);

        if (WALKTYPE!=WALK_CA) Disable(Vtrj_Sphere);
		if (WALKTYPE==WALK_CA) Enable(Vtrj_Sphere);
		
		if (resinfo->CisTrajGraph==NULL) Disable(Vtrj_Cis);
		if(resinfo->CisTrajGraph!=NULL)	Enable(Vtrj_Cis);

		if (HasBeenFiltered==FALSE) Disable(Vtrj_None);
		if (HasBeenFiltered==TRUE) Enable(Vtrj_None);

	}

}

static void DisableAllItems(void)
{
	Boolean wasNew,wasLoaded;
	
	/* disable stuff */
	Disable(Vtrj_Exit);
	Disable(Vtrj_Help);
	Disable(Vtrj_About);
	Disable(resid_switch);
	Disable(jumpto_text);
	Disable(Vtrj_WebSite);
	Disable(Vtrj_BgColor);
	wasLoaded=isLoaded;
	wasNew=isNew;
	isLoaded=FALSE;
	isNew=TRUE;
	EnableDisableProc();
	isLoaded=wasLoaded;
	isNew=wasNew;
	alldisabled=TRUE;
}

static void EnableAllItems(void)
{
	EnableDisableProc();
	Enable(Vtrj_Exit);
	Enable(Vtrj_Help);
	Enable(Vtrj_About);
	Enable(Vtrj_WebSite);
	Enable(Vtrj_BgColor);
	Enable(resid_switch);
	Enable(jumpto_text);
	alldisabled=FALSE;
}

/*****************************************************************
			                KEYBOARD FUNCTIONS
*****************************************************************/

static void Vtrj_ModualKeyPressed(SlatE s, Char key)
{

}

static void Vtrj_MainKeyPressed(SlatE s, Char key)
{
	
	if (isLoaded ==TRUE && ProgramProgress<=0)
	{
	  	
		
		/*Right Arrow Key*/
		if ((int) key == 29)
		{
			res_num += 1;
			if (res_num>numAA) res_num=1;
			SetValue(resid_switch, res_num);
			SwitchResProc(TRUE);
		}
	
		/*Left Arrow Key*/
		else if ((int) key == 28)
		{
			res_num -= 1;
			if (res_num < 1) res_num = numAA;
			SetValue(resid_switch, res_num);
			SwitchResProc(TRUE);
		}
	}
}


/*****************************************************************
			                MAIN
*****************************************************************/


extern Int2 Main(void)
{
	Int4 processid,randseed;
	Char progname[PATH_MAX];
	Char fname[PATH_MAX];
	Int4 argc;
	Int2 winx,winy;
	Char **argv;
	Char cHere;

	/*Corrects a runtime error*/
	glutInit(&argc,argv);
	
	ErrSetFatalLevel(SEV_FATAL);
	ErrSetMessageLevel (SEV_FATAL);
	ErrSetLogLevel (SEV_ERROR);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
	ErrSetLogfile ("vtrj_error.log", ELOG_BANNER | ELOG_APPEND);
  /* grab Ctrl-C */
#ifdef OS_UNIX
  signal(SIGINT, MyInterruptHandler);
  signal(SIGSEGV, MyInterruptHandler);
  signal(SIGKILL, MyInterruptHandler);
  signal(SIGTERM, MyInterruptHandler);
  signal(SIGALRM, MyInterruptHandler);
#endif                                                                                                    		
	if (!OPENMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_ERROR,2,1,"Unable to open MMDBAPI");
		return 1;
	}
	USE_LOTS_RAM=FALSE;
	BUILD_FRAGS_ONLY=FALSE;
	/* save working directory */
	getcwd(CFG_local_datafilepath,PATH_MAX);
	/* append slash */
	cHere=CFG_local_datafilepath[StringLen(CFG_local_datafilepath)-1];
	if (cHere!=DIRDELIMCHR) {
		StringCat(CFG_local_datafilepath,DIRDELIMSTR);
	}
	processid = GetAppProcessID();
	randseed=GetSecs()*processid;
	RandomSeed(randseed);
	sprintf(progname,"%s V%s",VISTRAJ_EXE_NAME,VISTRAJ_VERSION);
	SetAppProperty ("AppName", progname);
	winy=screenRect.bottom-screenRect.top-100;
	if (winy>800)
		winy=800;
	winx=winy;
	if (winx<640)
		winx=640;
	main_win = DocumentWindow (0,0,winx,winy, progname, (WndActnProc) ExitProc, (WndActnProc)ResizeProc);
	SetupMenus ();
	SetupWindow ();
	SetupOGL();

	/*Sets which function is called when a key is pressed*/
	SetSlateChar ((SlatE) myOGL_data->Panel, Vtrj_MainKeyPressed);
		
	EnableDisableProc();
	
	StringCpy(prog_name, VISTRAJ_EXE_NAME);
	StringCat(prog_name, " V");
	StringCat(prog_name, VISTRAJ_VERSION);
	StringCat(prog_name, " - ");
	Show (main_win);
	
	/*Sets which function is called when a mouse button is pressed*/
	SetPanelClick(myOGL_data->Panel, Vtrj_ClickProc, Vtrj_DragProc, Vtrj_HoldProc, Vtrj_ReleaseProc);	
/*	SetPanelClick(aaname_panel, DrawAANameProc, NULL, NULL, NULL);*/

  argc=GetArgc();
  if (argc>1) {
  	argv=GetArgv();
  	StringCpy(fname,argv[1]);
  	/* append .trj if not given */
  	if (StringStr(fname,ASN_EXT)==NULL)
  		StringCat(fname,ASN_EXT);
	StringCpy(GlobalPath, fname);
  	LoadProc(fname);
  	SetValue(resid_switch,1);
  }
	/*GetArgs ("vis_traj", DIM(Vtrj_args),Vtrj_args);*/  /*Very Annoying Feature*/
/*	if (Vtrj_args[0].strvalue!=NULL)
	{
		if (Vtrj_args[1].intvalue!=0)  JumpFromStart=TRUE;
		SetValue(resid_switch, (Int2) Vtrj_args[1].intvalue);
		LoadProc(Vtrj_args[0].strvalue);
	}
*/
	GenerateModAAList();
	
	ProcessEvents ();
	
	return 0;
}


/*  
$Log: vis_traj.c,v $
Revision 1.88  2008/10/08 17:18:13  mingxi

Corrects a runtime error concerning "GlutInit" in vis_traj.c

Revision 1.87  2004/07/13 19:36:01  hfeldman
Show full fragment probability now if < 0.01, and fixed memory issue with importing fragments

Revision 1.86  2004/06/23 17:35:37  hfeldman
Added global frags define

Revision 1.85  2003/12/08 22:49:30  hfeldman
Added preliminary support for import of fragments from text file

Revision 1.84  2003/09/23 19:13:03  feldman
Added errorfile option to foldtraj for homtraj
Made foldtraj gradually increase laxness if gets stuck a lot, like in DFP

Revision 1.83  2003/09/23 17:09:30  feldman
Added error message output for homtraj

Revision 1.82  2003/08/21 18:12:49  feldman
changed pauses to sleeps

Revision 1.81  2003/04/04 21:51:05  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.80  2003/03/25 17:51:16  feldman
Added generation param to foldtraj

Revision 1.79  2003/03/25 17:46:19  feldman
Added zhang window param to maketrj

Revision 1.78  2003/02/27 21:32:28  feldman
Added needed dummy functions for DFP

Revision 1.77  2003/01/16 14:43:21  feldman
Fixed PHIPSI walk cis pro editing bug

Revision 1.76  2002/07/25 16:30:21  feldman
Added tunnel prob

Revision 1.75  2002/07/20 16:40:50  feldman
Changed location of tginit

Revision 1.74  2002/07/14 21:49:30  feldman
Fixed CalculateInfo function for WALK_CA

Revision 1.73  2002/06/14 20:47:17  feldman
Added fragment length to data structures

Revision 1.72  2002/03/14 21:55:49  feldman
Removed unneeded functions

Revision 1.71  2002/01/24 23:36:11  feldman
Added fix for OSX

Revision 1.70  2001/12/03 15:54:55  feldman
UNITS_UNKNOWN is default units now

Revision 1.69  2001/10/19 14:42:42  feldman
Fixed bug - didnt save residue after removing residue fragments

Revision 1.68  2001/09/10 20:25:36  feldman
Added U-turns to trajectory creation

Revision 1.67  2001/09/07 21:47:35  feldman
Fixed bug with single peak placement

Revision 1.66  2001/08/16 21:23:04  feldman
fixed seg. fault

Revision 1.65  2001/08/16 21:17:45  feldman
Change color when edit fragments too

Revision 1.64  2001/08/16 20:54:18  feldman
Color sequence characters different is have fragments attached

Revision 1.63  2001/07/12 21:09:12  feldman
Fixed command line file loading bug

Revision 1.62  2001/06/19 16:39:13  feldman
Changed error level, and corrected minor spelling mistake bug

Revision 1.61  2001/06/18 17:33:30  feldman
added file path to extaa file

Revision 1.60  2001/06/18 16:06:29  feldman
Search for help file in Data Path

Revision 1.59  2001/06/14 19:19:00  feldman
changed vis_traj to size itself according to your current display resolution

Revision 1.58  2001/05/25 21:49:36  feldman
Added lighting and materials

Revision 1.57  2001/05/03 22:02:14  feldman
removed stray 0

Revision 1.56  2001/05/03 18:19:58  feldman
chnaged behaviour of Apply slightly for Residue Properties

Revision 1.55  2001/05/02 16:22:18  feldman
Changed Apply button to behave like Windows Apply buttons (i.e. only
hit once when all changes are done)

Revision 1.54  2001/04/26 19:49:09  feldman
Fixed some minor Bioseq bugs and potential bugs

Revision 1.53  2001/04/26 16:00:17  feldman
Ensure menus are disabled on UNIX when opening multiple modal windows

Revision 1.52  2001/04/18 22:33:51  feldman
Windows-speific fixes, and made so distance constraints and fragments are
auto-selected (you dont have to click on one at the start)

Revision 1.51  2001/04/18 21:24:18  feldman
Added fragment editing menu

Revision 1.49  2001/04/02 15:26:34  feldman
fixed aboutproc message

Revision 1.48  2001/04/02 15:10:20  feldman
fixed more memory leaks

Revision 1.47  2001/04/01 19:44:42  feldman
made foldtrajlite Windows port using conio.h, began network portion

Revision 1.46  2001/03/29 16:11:55  feldman
Added date and version to about box

Revision 1.45  2001/03/29 02:52:24  feldman
fixed minor compiler warnings

Revision 1.44  2001/03/27 20:46:00  feldman
removed unused variables

Revision 1.43  2001/03/27 20:24:22  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.42  2001/03/23 21:08:07  feldman
Added monitor for foldtraj inside vistraj

Revision 1.40  2001/03/22 22:14:34  feldman
cleaned up a bit more garbage of unneeded functions

Revision 1.39  2001/03/22 18:59:17  feldman
Removed RPS Blast option when input is MMDB structure (didn't make sense)

Revision 1.38  2001/03/19 19:03:44  feldman
Fixed window size on win 98

Revision 1.37  2001/03/19 16:44:22  feldman
Fixed a few compiler warnings and a small memory leak

Revision 1.36  2001/03/19 15:13:14  feldman
fixed slight peak coloring bug

Revision 1.35  2001/03/19 14:30:18  feldman
fixed a few more minor bugs

Revision 1.34  2001/03/18 18:25:20  feldman
put back mmdblocl where needed.. alternative is LIB3 and LIB5 which make bulkier executable

Revision 1.33  2001/03/16 16:43:07  feldman
minor cosmetic improvements

Revision 1.32  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.31  2001/03/13 17:49:24  feldman
added to about box, allow distance constraints to be edited always now

Revision 1.30  2001/03/13 15:08:17  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.29  2001/03/12 17:10:50  feldman
Fixed compiler warnings (again)

Revision 1.28  2001/03/09 17:33:59  feldman
Added code to deal with rotamers.  Now rotamer chi angles can
be encoded into a 4-byte integer and used by foldtraj when
building the structures.  These can be edited by VisTraj as well

Revision 1.27  2001/03/07 21:49:48  feldman
Made many important fixes to VisTraj such as bounds checking of
user-entered information, and tidied up interface further.
Fixed a few minor bugs as well, and tested RPSTraj feature of
VisTraj.

Revision 1.23  2001/02/27 20:51:14  feldman
-minor bugfixes
-added background color changing

Revision 1.22  2001/02/23 18:03:54  feldman
Changed MakeTrj to void*, void*

Revision 1.21  2001/02/23 16:41:58  feldman
Fixed minor warning

Revision 1.20  2001/02/15 21:10:00  feldman
Fixed amino acid font

Revision 1.19  2001/02/15 20:43:05  feldman
Split vis_traj into two files, to separated MMDBAPI from X stuff,
and added minor tweaks, progress monitor and bugfixes

Revision 1.16  2001/02/06 18:39:43  feldman
Major overhaul of UI, added distance constraints, maketrj functionality and more

Revision 1.14  2000/08/24 22:21:47  feldman
Fixed up sphere stuff, now correct.
Note that previous residue is at sphere centre, and
previous previous residue is at SOUTH pole (not north
as stated earlier), nand 3rd last is along prime meridian.
Now theta is indeed the dihedral of the 4 Ca's (right sign
for sure) and phi is the supplement of the angle between the
last three atoms as stated earlier.

Revision 1.13  2000/08/16 19:46:22  john
Updated code to work with new NCBI toolkit

Revision 1.12  2000/08/14 20:36:29  john
Made windows compatible.  Fixed plane view rotation in windows environment

Revision 1.7  2000/07/10 18:12:07  john
Updated code to work with Howards new headers and altered functions

Revision 1.1.1.1  2000/06/12 15:50:27  feldman
OpenGL Trajectory Distribution Visualizer

*/
