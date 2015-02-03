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

#ifndef VIS_TRAJ_H
#define VIS_TRAJ_H

#ifdef WIN_MSWIN
#include <windows.h>
#include <direct.h>
#include <gl/glut.h>
#endif /*WIN_MSWIN*/

#ifdef WIN_MOTIF
#include <GL/glut.h>
#include <netscape.h>  /* defines OS_UNIX */
#endif  /*WIN_MOTIF*/

#ifdef WIN_MAC
#include <GLUT/glut.h>
#endif

#ifdef OS_UNIX
#include <signal.h>
#endif

#ifdef WIN_MSWIN
#undef WIN_MSWIN
#endif

#include <stdarg.h>
#include <hfprogs.h>
#include <conform.h>
#include <document.h>
#include <prtutil.h>
#include "vtrj_ncbi.h"

#include <ncbithr.h>

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************
			DEFINES
*****************************************************************/
#define VISTRAJ_VERSION "1.1.2"
#define VISTRAJ_EXE_NAME "VisTraj"
#define VTRJ_HELP_FILE "vtrj.hlp"
#define GENERATE_TMP_FILENAME "last_generated"
#define MODAA_C 1
#define MODAA_M 2
#define MODAA_N 3
#define MAX_HELP_LIST 200
#define MAX_AA_LINES 200
#define SCALESPHERE 5.00
#define SCALE_SIZE_XY 30.000
#define SCALEZ 12.00
#define MOUSE_SENSITIVITY 5
#define VTRJ_COLOR_MAX 37
#ifdef OS_MSWIN
#define TASK_CFGFILE ".\\foldtraj"
#else
#define TASK_CFGFILE "foldtraj"
#endif
#define DUMPTYPE_PDB 1
#define DUMPTYPE_VAL 2

/*****************************************************************
			COLOR DEFINITIONS
*****************************************************************/
#define MAGENTA 2
#define BLUE    4
#define GREEN   8
#define YELLOW  9
#define RED     12
#define WHITE   15
#define BLACK   16
/* colours below are currently unused */
#define PURPLE  3
#define SKY     5
#define GOLD    10
#define ORANGE  11
#define PINK    13
#define BROWN   21


/*****************************************************************
			TYPEDEFS
*****************************************************************/

typedef struct _TOGL_Layers
/* this struct contains the information used to manage the different layers of the display */
{
    GLuint FirstLayer;
    GLuint LastLayer;
    GLuint SelectedLayer;
    Nlm_Boolean IsOn[OGLMAXLAYERS];
} TOGL_Layers;

typedef struct Vtrj_PolyPoints
{
	GLfloat x1;
	GLfloat y1;
	GLfloat z1;
	GLfloat x2;
	GLfloat y2;
	GLfloat z2;
	GLfloat x3;
	GLfloat y3;
	GLfloat z3;
	GLfloat x4;
	GLfloat y4;
	GLfloat z4;
	int color;
	Boolean iszero;
	GLfloat maxpoint;
}  PolyPoints, *Ppp;

typedef struct Vtrj_globals
{
	FloatHi bbet;
	FloatHi bbp;
	Int2 nrt;
	FloatHi abbb;
	FloatHi absc;
	FloatHi is;
	FloatHi sbb;
	Int2 tgu;
	Int2 wt;
	Char cf[PATH_MAX];
	Int2 tt;
	Int4 td;
	Int4 to;
	Int2 bch;
	FloatHi msf;
	FloatHi tp;
} VOG, *VOGPtr;
	
typedef struct Vtrj_Newglobals
{
	FloatHi bbet;
	FloatHi bbp;
	Int2 nrt;
	FloatHi abbb;
	FloatHi absc;
	FloatHi is;
	FloatHi sbb;
	Int2 tgu;
	Char wt[PATH_MAX];
	Char cf[PATH_MAX];
	Int2 tt;
	Int4 td;
	Int4 to;
	Char bch[10];
	FloatHi msf;
	FloatHi tp;
} VNG, *VNGPtr;

/*****************************************************************
			GLOBAL VARIABLES
*****************************************************************/

/*Args Vtrj_args[] = {
			{"Input Trajectory Distribution filename",NULL,NULL,NULL,TRUE,'f', ARG_STRING,0.0,0,NULL},
			{"Input Residue Number",NULL,NULL,NULL,TRUE,'r', ARG_INT,0.0,1,NULL}
		};
*/

Int2 Vtrj_bgcolor=BLACK;
static Boolean fragwinopenres=FALSE,fragwinopenadd=FALSE,constrwinopen=FALSE,alldisabled=FALSE;

WindoW main_win;
static TexT jumpto_text, gaus_mag_text, gaus_sd_text,
	lp_mag_text, lp_sd_text, REdit_Info, GEdit_Info, GNewEdit_Info, newAA_info;

static Int2 numchi[20]={0,4,2,2,1,3,3,0,2,2,2,4,3,2,2,1,1,2,2,1};
static Int2 lastg=0;
static Char prog_name[PATH_MAX];		
static ButtoN go_bttn;
static SwitcH resid_switch;
static PaneL aaname_panel;
static LisT editGlist, editRlist, modifyRlist, helpmenulist, editGNewlist, newAAlist, distlist, fraglist, fragreslist;
static Int2 res_num;
static CharPtr AAposition;

static PrompT resnum_title, /*xval_title, yval_title, zval_title,*/
	modAA_title, modAA_temptitle;

static Int2 editwhatG, editwhatR;
static PoinT whereclicked, wheredragged;
static Boolean clicked = FALSE;
static Boolean dragged = FALSE;
static Boolean inHelp=FALSE;
static Boolean *hasafrag=NULL;

static ChoicE which_view, which_draw, which_range, which_trans,
	which_shade, which_rot;

static FloatLo ThePeak = 0.00;
static FloatLo rangemin = 0.00;
static FloatLo rangemax = 0.00;
static GLfloat ScaleZ, ScaleSphere;
static int Vtrj_rotateX= 90, Vtrj_rotateY=0, Vtrj_rotateZ=0;
static GLfloat Vtrj_scale = 10.0;
static GLfloat Vtrj_Xratio = 1.0, Vtrj_Yratio = 1.0;

static FloatLo gaus_mag=10.00, gaus_sd=1.00;
static FloatLo lp_mag=10.00, lp_sd=1.50;

static GLenum drawtype=GL_QUADS;
static GLenum filltype=GLU_FILL;
static GLint viewport_properties[4];

/*static GLdouble mod_matrix[16];
static GLdouble proj_matrix[16];*/

static IteM Vtrj_Load, Vtrj_Save, Vtrj_Exit, Vtrj_Global, Vtrj_Residue,
	Vtrj_None, Vtrj_Smooth, Vtrj_Gaussian, Vtrj_LP, Vtrj_Sphere,
	Vtrj_Plane, Vtrj_Filled, Vtrj_Wire, Vtrj_Range1, Vtrj_Range2,
	Vtrj_Range3, Vtrj_Range4, Vtrj_Range5, Vtrj_Range6,
	Vtrj_GenerateVal, Vtrj_GeneratePdb, Vtrj_Protein, Vtrj_About, Vtrj_Help, Vtrj_PosXSide,
	Vtrj_PosNorm, Vtrj_PosYSide, Vtrj_PosTop, Vtrj_PosBottom,
	Vtrj_Trans, Vtrj_Cis, Vtrj_SavePNG, Vtrj_Noise, Vtrj_AddRes,
	Vtrj_RemoveRes, Vtrj_QuickOn, Vtrj_Frag, Vtrj_ImportFrag,
	Vtrj_QuickOff, Vtrj_Dist_Cons, Vtrj_Mutate, Vtrj_SaveAs,
	Vtrj_Close, Vtrj_New, Vtrj_WebSite, Vtrj_BgColor, Vtrj_Shade, Vtrj_NoShade/*, Vtrj_Multiple*/;

/*shading global variables*/
/* light is at infinity */
static GLfloat lightPosition[] = {1.0, 1.0, 1.0, 0.0};
/*static GLfloat lightPosition2[] = {20.0, 0.0, 0.0, 1.0};*/
static GLfloat lightModel[] = {1.0, 1.0, 1.0, 1.0};
static GLfloat matColorAmbGreen[] = {0.0615, 0.3745, 0.0615, 1.0};
static GLfloat matColorDifGreen[] = {0.02568, 0.21424, 0.02568, 1.0};
static GLfloat matColorSpecGreen[] = {0.633, 0.727811, 0.633, 1.0};
static GLfloat shineGreen[] = {76.8};
static GLfloat matColorAmbBlue[] = {0.0, 0.0, 0.65, 1.0};
static GLfloat matColorDifBlue[] = {0.0, 0.0, 0.8, 1.0};
static GLfloat matColorSpecBlue[] = {0.6, 0.6, 0.7, 1.0};
static GLfloat shineBlue[] = {32.0};
static GLfloat matColorAmbPurple[] = {0.35, 0.0, 0.35, 1.0};
static GLfloat matColorDifPurple[] = {0.8, 0.0, 0.8, 1.0};
static GLfloat matColorSpecPurple[] = {0.7, 0.6, 0.7, 1.0};
static GLfloat shinePurple[] = {32.0};
static GLfloat matColorAmbRed[] = {0.6, 0.0, 0.0, 1.0};
static GLfloat matColorDifRed[] = {0.25, 0.0, 0.0, 1.0};
static GLfloat matColorSpecRed[] = {0.7, 0.6, 0.6, 1.0};
static GLfloat shineRed[] = {32.0};
static GLfloat matColorAmbYellow[] = {0.66, 0.44, 0.06, 1.0};
static GLfloat matColorDifYellow[] = {0.28, 0.17, 0.03, 1.0};
static GLfloat matColorSpecYellow[] = {0.99, 0.91, 0.81, 1.0};
static GLfloat shineYellow[] = {27.8};
static GLfloat matColorAmbWhite[] = {0.9, 0.9, 0.9, 1.0};
static GLfloat matColorDifWhite[] = {0.78, 0.78, 0.78, 1.0};
static GLfloat matColorSpecWhite[] = {0.91, 0.91, 0.91, 1.0};
static GLfloat shineWhite[] = {27.8};
static GLfloat matColorAmbBlack[] = {0.13, 0.13, 0.13, 1.0};
static GLfloat matColorDifBlack[] = {0.11, 0.11, 0.11, 1.0};
static GLfloat matColorSpecBlack[] = {0.5, 0.5, 0.5, 1.0};
static GLfloat shineBlack[] = {27.8};
static GLfloat lightColorAmb[] = {1.0, 1.0, 1.0, 1.0};
static GLfloat lightColorDif[] = {1.0, 0.9, 0.5, 1.0};
static GLfloat lightColorSpec[] = {1.0, 0.9, 0.5, 1.0};


/*help option global variables*/

DoC help_doc, pinfo_doc, dist_doc, frag_doc, fraglist_doc;

/*static ParData lstParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
*/
static ColData lstColFmt = {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

static ParData hdgParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ParData DishdgParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData hdgColFmt = {600, 0, 80, 0, NULL, 'c', TRUE, FALSE, FALSE, FALSE, TRUE};

/*static ParData subParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
*/
static ColData subColFmt = {750, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

/*static ParData tblParFmt = {TRUE, FALSE, FALSE, FALSE, TRUE, 0, 0};
*/
static ColData tblColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static ParData txtParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData txtColFmt = {750, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

/*distance constraints and fragment global variables*/

static ColData FraghdgColFmt = {100, 0, 35, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};
static ColData DishdgColFmt = {250, 0, 80, 0, NULL, 'c', TRUE, FALSE, FALSE, FALSE, TRUE};
static ColData DissubColFmt = {250, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};
static ColData DistxtColFmt = {500, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

Char Vtrj_Help_Menu[MAX_HELP_LIST][MAXCOL];

PNN Vtrj_pnn=NULL;


WindoW editGinfo_win, editRinfo_win, GausAttr_win, LPAttr_win,
	addnoise_win, modifyR_win, help_win, pinfo_win, MultView_win, editdist_win, adddist_win,
	addfrag_win, editfrag_win, editfraglist_win;

WindoW NewWiz1_win, NewWiz2_win, NewWiz3S_win, NewWiz3V_win, NewWiz4_win, newAA_win,Background_win;
	
TexT currentRvalue_text, currentGvalue_text, currentGNewvalue_text, gor_h_text, gor_c_text, gor_e_text,
	NewWiz2_text, NewWiz3V_text, NewWiz3S_text, NewWiz3V2_text, NewWiz3S2_text,
	NewWiz3S3_text;

TexT chain_text, endres_text, startres_text, temp_text, sdx_text, sdy_text, peakheight_text,
	timestep_text, model_text;
	
PopuP multnum_list;

PrompT 	prop_title, gor_h_title, gor_c_title, gor_e_title, AminoEdit_title, seqlen_title;

static TexT fraglen_text, fragres_text, fragaa_text, fraga1_text, fraga2_text, fragasd_text, fragchiw_text,
	fragchiwsd_text, fragchi1_text, fragchi2_text, fragchi3_text, fragchi4_text, fragprob_text;

static TexT resnum1_text, resnum2_text, resname1_text, resname2_text, meandist_text, mindelta_text,
	maxdelta_text, angle1_text, angle2_text, dihedral01_text, dihedral12_text, dihedral23_text,
	prob_text;

GrouP noise_type, add_or_rep;

GrouP Wiz1InputType, SequenceType, UseSSFile, UseSSorTT, SS_Group, TT_Group, Gor_Group, SCallRPS,
	VCallRPS, SCallUturn;
	
ButtoN applyR_bttn, applyG_bttn, noise_ok_bttn, modifyR_bttn, okR_bttn, MultOK_bttn, NewGApply_bttn;

ButtoN Wiz2Next_bttn, Wiz3VNext_bttn, Wiz3SNext_bttn, okAA_bttn, distE_bttn, distA_bttn,
	 distD_bttn, fragE_bttn, fragA_bttn, fragD_bttn, fraglistE_bttn, fraglistA_bttn, fraglistD_bttn, BGOK_bttn;

PopuP savechiV_list, comptypeV_list, comptypeS_list, TGtype_list, Newwalktype_popup, Newbchydrogen_popup,
	walktype_popup, bchydrogen_popup, tgunits_popup, ttype_popup, background_popup;

Int2 NewAAlen=0;
Int4 frag_restoadd=0;

Boolean isLoaded = FALSE, isCis=FALSE, isRused=FALSE, ResReplace=TRUE, isAAused=FALSE,
	isFirstX=FALSE, isSwitching=FALSE, isChanging=FALSE, Vtrj_redraw=FALSE;
Boolean Rhasbeenedited = FALSE, Ghasbeenedited = FALSE, Nhasbeenedited = FALSE;
Boolean ShadingOn = FALSE/*, JumpFromStart=FALSE*/;
Boolean HasBeenFiltered=FALSE, QuickRotate=FALSE, Drawit=TRUE,needApply=FALSE;
Boolean AAmodified=FALSE, isNew=FALSE, editdist=FALSE, DistBeenEdited=FALSE, editfrag=FALSE, FragBeenEdited=FALSE;

ValNodePtr CmodAAlist = NULL;
ValNodePtr NmodAAlist = NULL;
ValNodePtr MmodAAlist = NULL;

Char shortsequence[10];	
Char pointedAA[2];

MonitorPtr Vtrj_mon=NULL;

TOGL_Data  *myOGL_data;

Uint1 Vtrj_RGB[VTRJ_COLOR_MAX][3] =
{
  {255, 255, 255}, /* default     0 */
  {255,  20, 147}, /* hotpink     1 */
  {255,   0, 255}, /* magenta     2 */
  {155,	48, 255}, /* purple      3 */
   {0,   0, 255}, /* blue        4 */
   {30, 144, 255}, /* sky         5 */
    {0, 255, 255}, /* cyan        6 */
    {0, 255, 127}, /* sea         7 */
    {0, 255,   0}, /* green       8 */
  {255, 255,   0}, /* yellow      9 */
  {255, 165,   0}, /* gold       10 */
  {255,  69,   0}, /* orange     11 */
  {255,   0,   0}, /* red        12 */
  {255, 114,  86}, /* pink       13 */
  {255, 174, 185}, /* pinktint   14 */
  {255, 255, 255}, /* white      15 */
    {0,   0,   0}, /* black      16 */
  {176, 226, 255}, /* bluetint   17 */
  {154, 255, 154}, /* greentint  18 */
  {255, 236, 139}, /* yellowtint 19 */
  {125, 125, 125}, /* gray       20 */
  {139,  87,  66}, /* brown      21 */
  {255,   0,   0}, /* perm colors 22: red */
    {0, 255,   0}, /* perm colors 23: green */
  {255,   0, 255}, /* perm colors 24: magenta */
   {30, 144, 255}, /* perm colors 25: sky */
  {155,	48, 255}, /* perm colors 26: purple */
   { 0, 255,   0}, /* SS colors 27: helix, green */
  {255, 165,   0}, /* SS colors 28: strand, gold */
  {255,  69,   0}, /* SS colors 29: turn, orange */
   { 0, 255, 255}, /* SS colors 30: coil, cyan */
  {255, 255,   0}, /* highlight colors 31: yellow */
   { 0,   0,   0}, /* background colors 32: black */
  {160,  82,  45}, /* sienna 33 */
  {240, 230, 140}, /* khaki 34 */
  {171, 130, 255}, /* light purple 35 */
  {255, 200,   0}  /* orange yellow 36 */
};


static Ppp GPVp = NULL;
static VOGPtr OrigGlobs = NULL;
static VNGPtr NewGlobs = NULL;

static PopuP sstru_list;

static Char GlobalPath[PATH_MAX];
Boolean ShouldSaveAs = TRUE;

CharPtr fnamtrj=NULL;
Int2 numAA;
Char sequence[MAXSIZE];
CharPtr Xsequence;
CharPtr Psequence;
CharPtr newsequence;
Int2 modresid;
Char modresParent = '\0';

BiostrucIdPtr	pbi;
BiostrucDescrPtr	pbd;
ValNodePtr vnpBioseq;

static PTGS resinfo = NULL;
static PTGS origresinfo = NULL;
static PTGS reseditinfo = NULL;
static PTGS trajnoise = NULL;

static Char Vtrj_sstru = '\0';
static FloatLo Vtrj_prob[3] = {0.0, 0.0, 0.0};
static int Vtrj_SPposition = 0;
static FloatLo Vtrj_SPpeak = 0.0;
static Char modAA[MAX_AA_LINES] [MAX_KBAA_NAMELENGTH+1];
static Char descrAA[MAX_AA_LINES] [PATH_MAX];

Int2 tmp_Glob;

/* Visual Make Traj Global Variables */

Char VMkTrj_valin[PATH_MAX];
Char VMkTrj_trjout[PATH_MAX];
Char VMkTrj_chain[100];
Int2 VMkTrj_model;
Int2 VMkTrj_startres;
Int2 VMkTrj_endres;
Int2 VMkTrj_peakheight;
Int2 VMkTrj_savechi;
FloatLo VMkTrj_sdx;
FloatLo VMkTrj_sdy;
FloatLo VMkTrj_temp;
FloatLo VMkTrj_timestep;
Boolean VMkTrj_quiteop;
Int2 VMkTrj_noise;
Int2 VMkTrj_method;
Char VMkTrj_seqfile[PATH_MAX];
Char VMkTrj_ssout[PATH_MAX];
Int2 VMkTrj_tgtype;
Int2 VMkTrj_comprtype;
Char VMkTrj_constraints[PATH_MAX];
Boolean VMkTrj_goroverride;
Int2 VMkTrj_seqtype;
Char VMkTrj_AAseq[MAXSIZE];
Boolean VMkTrj_dorps;
Boolean VMkTrj_douturn;
Int2 VMkTrj_SSorTT;

TCn3D_ColorData Vtrj_ColorData;

/*****************************************************************
			FUNCTION DECLARATIONS
*****************************************************************/
static void Vtrj_GetPoints(void);
static void RangeProc(void);
static Int4 SwitchResProc(Boolean IsSwitching);
static void ApplyREditProc (IteM i);
static void ApplyGNewEditProc(IteM i);
static void Vtrj_Move(void);
static void Vtrj_Zoom(void);
static void Vtrj_Rotate(void);
static void Vtrj_DrawAll(void);
static void EnableDisableProc(void);
static void DisableAllItems(void);
static void EnableAllItems(void);
static void Vtrj_Draw_Poly(GLenum hit, Boolean Sphere);
static void LoadProc (CharPtr path);
static void GetModResInSequence(void);
static void Vtrj_LoadSequence(void);
static void Vtrj_SaveResidue(void);
static void Vtrj_InsertResidue(Char the_aa);
static void CloseTrajProc (IteM i);
static Boolean CheckIfSaved(void);
static void CheckInputTypeProc(IteM i);
static void CheckSequenceTypeProc(IteM i);
static void Vtrj_LoadVertexesSphere(void);
static void Vtrj_LoadVertexesPlane(void);
static void Vtrj_LoadVertexes(void);
static void Vtrj_ModualKeyPressed(SlatE s, Char key);
static Int2 GetAAPos(void);
static void LoadAAList(Int2 WhereModAA);
static void DisplayNewWiz1(void);
static void DisplayNewWiz2(IteM i);
static void DisplayNewWiz3Val(IteM i);
static void DisplayNewWiz3Seq(IteM i);
static void DisplayNewWiz4(IteM i);
static void LoadInfoTextProc (Int2 Vtrj_topic, Int2 Vtrj_heading);
static void LoadDistList(void);
static void LoadDistProc(IteM i);
static void LoadFragList(void);
static void LoadFragHere(void);
static void LoadFragProc(IteM i);
static void LoadFragListProc(IteM i);

void SetupOGL(void);
void Vtrj_PostPicture(void);
void Vtrj_ClearScreen(void);

void Vtrj_PutText(CharPtr string, int color, GLfloat x, GLfloat y, GLfloat size,
	GLfloat width, GLfloat rotz);
Nlm_PaneL  Nlm_Autonomous3DPanel PROTO((Nlm_GrouP prnt, Nlm_Int2 pixwidth,
                                        Nlm_Int2 pixheight, Nlm_PnlActnProc draw,
                                        Nlm_SltScrlProc vscrl, Nlm_SltScrlProc hscrl,
                                        Nlm_Int2 extra, Nlm_PnlActnProc reset,
                                        Nlm_GphPrcsPtr classPtr, Nlm_Boolean *IndexMode,
                                        void **display, void **visinfo));
/*
static void VisHeightProc(IteM i);
static void ApplyREditProc (IteM i);
static void InitializeOGLDraw(void);
void Vtrj_ResizeProc(WindoW w);
static void NormViewProc (IteM i);
static void XViewProc (IteM i);
static void YViewProc (IteM i);
static void TopViewProc (IteM i);
static void BottomViewProc (IteM i);
static void Vtrj_ClearXYZ(void);
static void LoadREditProc (IteM i);
static void Vtrj_CalculateInfo(Boolean CalcTout);

void Nlm_SaveImagePNG(Nlm_Char *fname);

NLM_EXTERN void LIBCALL Cn3D_Redraw(Boolean New);
*/

#ifdef __cplusplus
}
#endif

#endif

/*  
$Log: vis_traj.h,v $
Revision 1.52  2003/12/08 22:49:31  hfeldman
Added preliminary support for import of fragments from text file

Revision 1.51  2003/11/10 18:59:41  feldman
Added includes for Mac OS X

Revision 1.50  2003/07/14 20:01:01  egarderm
Changed TASK_CFGFILE to .\foldtraj on Windows to allow local INI file

Revision 1.49  2003/04/04 21:51:06  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.48  2002/07/25 16:30:22  feldman
Added tunnel prob

Revision 1.47  2002/06/14 20:47:17  feldman
Added fragment length to data structures

Revision 1.46  2001/09/10 20:25:36  feldman
Added U-turns to trajectory creation

Revision 1.45  2001/08/16 20:54:19  feldman
Color sequence characters different is have fragments attached

Revision 1.44  2001/06/22 17:32:11  feldman
Updated to v1.1.1

Revision 1.43  2001/06/18 16:06:30  feldman
Search for help file in Data Path

Revision 1.42  2001/05/25 21:49:36  feldman
Added lighting and materials

Revision 1.41  2001/05/02 16:22:19  feldman
Changed Apply button to behave like Windows Apply buttons (i.e. only
hit once when all changes are done)

Revision 1.40  2001/04/26 16:00:18  feldman
Ensure menus are disabled on UNIX when opening multiple modal windows

Revision 1.39  2001/04/18 21:24:19  feldman
Added fragment editing menu

Revision 1.38  2001/04/17 21:32:06  feldman
Added fragment editing to Vistraj - not quite fully functional yet
Also can now give trajectory graph name as argument #1 when running
and fixed potential sscanf problems so entering garbage into numeric
entry fields results in an error being shown (instead of accepting it)

Revision 1.37  2001/04/02 15:10:21  feldman
fixed more memory leaks

Revision 1.36  2001/04/01 19:44:42  feldman
made foldtrajlite Windows port using conio.h, began network portion

Revision 1.35  2001/03/29 20:56:59  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.34  2001/03/29 16:11:55  feldman
Added date and version to about box

Revision 1.33  2001/03/28 17:29:28  feldman
removed errno from mobidick

Revision 1.32  2001/03/27 20:24:22  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.31  2001/03/23 21:08:07  feldman
Added monitor for foldtraj inside vistraj

Revision 1.30  2001/03/23 18:44:18  feldman
Integrated foldtraj into vistraj (!)
and cleaned up a few minor bugs

Revision 1.29  2001/03/22 22:14:35  feldman
cleaned up a bit more garbage of unneeded functions

Revision 1.28  2001/03/19 19:05:43  feldman
a few more unused variables

Revision 1.27  2001/03/19 16:44:22  feldman
Fixed a few compiler warnings and a small memory leak

Revision 1.26  2001/03/19 14:30:18  feldman
fixed a few more minor bugs

Revision 1.25  2001/03/16 16:43:07  feldman
minor cosmetic improvements

Revision 1.24  2001/03/14 19:16:12  feldman
removed more stuff from headers for windows

Revision 1.23  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.22  2001/03/13 15:08:17  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.21  2001/03/07 21:49:48  feldman
Made many important fixes to VisTraj such as bounds checking of
user-entered information, and tidied up interface further.
Fixed a few minor bugs as well, and tested RPSTraj feature of
VisTraj.

Revision 1.20  2001/03/02 16:43:35  feldman
More minor bugfixes and tweaks to improve usability

Revision 1.19  2001/02/27 20:51:15  feldman
-minor bugfixes
-added background color changing

Revision 1.18  2001/02/15 21:02:44  feldman
Fixed minor Windows compatibility issue

Revision 1.17  2001/02/15 20:43:05  feldman
Split vis_traj into two files, to separated MMDBAPI from X stuff,
and added minor tweaks, progress monitor and bugfixes

Revision 1.14  2001/02/06 18:39:44  feldman
Major overhaul of UI, added distance constraints, maketrj functionality and more

Revision 1.12  2000/10/27 21:09:54  feldman
Update vis_traj project file to work properly

Revision 1.11  2000/08/11 17:27:39  john
Added ifdefs for X window dependant code
Replaced all "grp" with "Vtrj_grp" due to Windows grp constants

Revision 1.7  2000/07/10 18:12:07  john
Updated code to work with Howards new headers and altered functions

Revision 1.5  2000/07/06 18:21:48  john
Added error handeling routine through ErrPostEx

Revision 1.1.1.1  2000/06/12 15:50:27  feldman
OpenGL Trajectory Distribution Visualizer

*/
