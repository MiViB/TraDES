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



#ifndef VTRJ_NCBI_H
#define VTRJ_NCBI_H

#include <vibmouse.h>
#include <ddvcolor.h>
/*
  Do not include gl.h here.  On Windows, gl.h is dependent on windows.h.
  However, some functions in vibrant have the same name as some of the
  functions in windows.h, so files that include both can have conflicts.
*/

#ifdef __cplusplus
extern "C" {
#endif

/* defines */
/* from shim3d.h */

#ifndef X_ROTATE_SBAR
#define X_ROTATE_SBAR 0x1
#endif

#ifndef Y_ROTATE_SBAR
#define Y_ROTATE_SBAR 0x2
#endif

#ifndef Z_ROTATE_SBAR
#define Z_ROTATE_SBAR 0x4
#endif

#define OGLMAXLAYERS 128

typedef enum {
    MouseOGL_DoNothing = 0,
    MouseOGL_RotateYX,
    MouseOGL_RotateZX,
    MouseOGL_RotateYZ,
    MouseOGL_Move,
    MouseOGL_Zoom,
    MouseOGL_NumStd
} Nlm_enumStdMAOGL;

typedef Boolean(*Nlm_MAInitOGLFunc) (MAPtr ma, Nlm_VoidPtr data);

typedef struct _OGL_BoundBox
/* bounds a volume */
{
    Nlm_FloatLo x[2];
    Nlm_FloatLo y[2];
    Nlm_FloatLo z[2];
    Nlm_Boolean set;
} TOGL_BoundBox;

struct _TOGL_Layers;

typedef struct _TOGL_Layers TOGL_Layers_;

typedef struct _OGL_Data
/* general runtime information */
{
    /* the current viewpoint. used for zoom and move */
    Nlm_FloatLo CameraDistance; /* distance (on Z-axis) from origin */
    Nlm_FloatLo CameraAngle;    /* in radians */
    Nlm_FloatLo CameraDirection[2]; /* point in Z=0 plane camera points at */
    Nlm_Boolean NeedCameraSetup;

    Nlm_FloatLo MaxSize;        /* biggest side of the bound box */
    TOGL_BoundBox BoundBox;     /* the containing box of the molecule */
    Nlm_WindoW ParentWindow;
    Nlm_PaneL Panel;            /* needed to set palette */
    ValNodePtr PaletteExpanded; /* the palette itself */
    ValNodePtr PaletteIndex;    /* palette index. type is TOGL_PaletteIndex */
    Nlm_BaR Z_rotate;           /* z rotation scroll bar */
    MAPtr ma;
    MA_GroupPtr ma_std_group[MouseOGL_NumStd];
    Nlm_VoidPtr ModelMatrix;    /* temporary copy of modelview matrix */
    Nlm_Boolean IsPlaying;      /* is the animation running? */
    Nlm_FloatHi Tick;
    TOGL_Layers_ *Layers;       /* the layers and their state */
    Nlm_Int4 SpaceWidth;        /* width of space character */
    Nlm_Int4 SpaceHeight;       /* height of space character */
    Nlm_Boolean IndexMode;      /* number of bits per pixel.  If < 16, used color index mode */
    DDV_ColorCell Background;   /* background color */
    /* note that the highlight color is the responsibility of the application */
    Nlm_Boolean SelectMode;     /* are we doing a selection? */
    Nlm_VoidPtr SelectBuffer;   /* buffer where the selection are dumped */
    Nlm_PoinT SelectPoint;      /* the point on the screen that was clicked for selection */
    Nlm_Uint4 SelectHits;       /* the number of hits */

    /* various info on X stuff related to OpenGL rendering context */
    /* these are void pointers, because including the X headers above this
       causes all sorts of name conflicts in various modules. Big pain! */
    void *display;              /* is actually a Display*      */
    void *visinfo;              /* is actually an XVisualInfo* */
} TOGL_Data;

/* from cn3dshim.h */

typedef struct _Cn3D_AnimateDlg {
    WindoW Cn3D_wAnimate;
} Cn3D_AnimateDlg;


typedef struct _Cn3D_ColorData {
    DDV_ColorGlobal *pDDVColorGlobal;
    Boolean IsUserData;     /* is the DDV_ColorGlobal userdata? */
    SeqAnnot *sap;          /* the current seqalign */
    ValNode *pvnsep;        /* the current seqentry */
    Boolean StandAlone;     /* is Cn3D running standalone? */
    Uint2 sapprocid, sepprocid, userkey;
    WindoW Cn3D_w;
    Cn3D_AnimateDlg AnimateDlg;
    Boolean UseEntrez;  /* turn on entrez use */
    Boolean EntrezOn;  /* is entrez on? */
		Int4 rows;
		Boolean AlignMode;
		Boolean IBM;
		IteM BlastMany;
		TOGL_Data *OGL_Data;
} TCn3D_ColorData;

extern Nlm_Boolean OGL_SetPosition3D(TOGL_Data * OGL_Data,
                                     Nlm_RectPtr rect);

#ifdef __cplusplus
}
#endif

#endif

/*
$Log: vtrj_ncbi.h,v $
Revision 1.7  2003/11/03 20:05:26  feldman
Added comment

Revision 1.6  2003/04/28 17:52:31  feldman
Updated for new toolkit

Revision 1.5  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.4  2001/02/15 20:43:05  feldman
Split vis_traj into two files, to separated MMDBAPI from X stuff,
and added minor tweaks, progress monitor and bugfixes

Revision 1.3  2001/01/25 22:25:26  feldman
Integrated Maketrj, changed lots of text messages, fixed
many minor bugs, improved overall interface

Revision 1.2  2000/08/16 15:56:13  feldman
Updated typedefs to match new NCBI typedefs...this caused the resizing on windows platform to behave properly.

Revision 1.1.1.1  2000/06/12 15:50:27  feldman
OpenGL Trajectory Distribution Visualizer

*/

