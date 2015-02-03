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

#ifndef VIS_TRAJGL_H
#define VIS_TRAJGL_H

#ifdef WIN_MSWIN
#include <windows.h>
#include <gl/gl.h>
#endif

#ifdef WIN_MOTIF
#include <GL/glx.h>
#include <Xm/Xm.h>
#endif

#ifdef WIN_MAC
#include <OpenGL/gl.h>
#endif

#include <vibrant.h>
#include "vtrj_ncbi.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN_MOTIF
#define Nlm_SlateTool  Widget
#endif /* WIN_MOTIF */

/*#ifdef WIN_MSWIN
#define Nlm_SlateTool HWND
#endif *//* WIN_MSWIN */

#if defined(WIN_MSWIN)
    HDC current_hdc = NULL;
	HGLRC current_hglrc = NULL;
/*	HDC hdc = NULL;
    HGDIOBJ current_hgdiobj = NULL;
    HBITMAP hbm = NULL;
    PIXELFORMATDESCRIPTOR pfd;
    int nPixelFormat;
    HGLRC hglrc = NULL;*/
#endif

extern TOGL_Data  *myOGL_data;
extern TCn3D_ColorData Vtrj_ColorData;
extern WindoW main_win;

#ifdef WIN_MOTIF
extern Display *Nlm_currentXDisplay;
GLXContext currentCtx;
GLXDrawable currentXdrw;
#endif

#ifdef WIN_MOTIF
NLM_EXTERN void Nlm_SetOGLContext (Nlm_SlateTool a, Nlm_Boolean *im, Display **d, XVisualInfo **v);
#endif  /*WIN_MOTIF*/

void Vtrj_ResizeProc(WindoW w);
TOGL_Data *Vtrj_OGL_CreateViewer(Nlm_GrouP prnt,
                            Uint2Ptr width, Uint2 height,
                            Int4 flags,
                            Nlm_MenU ma_group_menu,
                            Nlm_MenU ma_action_menu,
                            Nlm_MAInitOGLFunc ma_init_func,
                            VoidPtr ma_init_data);

#endif

#ifdef __cplusplus
}
#endif

/*

$Log: vis_trajgl.h,v $
Revision 1.6  2003/11/10 18:59:40  feldman
Added includes for Mac OS X

Revision 1.5  2001/03/14 19:16:12  feldman
removed more stuff from headers for windows

Revision 1.4  2001/03/14 16:25:56  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.3  2001/02/27 20:51:15  feldman
-minor bugfixes
-added background color changing

Revision 1.2  2001/02/15 21:02:44  feldman
Fixed minor Windows compatibility issue

Revision 1.1  2001/02/15 20:43:05  feldman
Split vis_traj into two files, to separated MMDBAPI from X stuff,
and added minor tweaks, progress monitor and bugfixes

*/

