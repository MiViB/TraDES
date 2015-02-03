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


#include "vis_trajgl.h"

static Boolean Cn3D_InitMA(MAPtr ma, VoidPtr data)
{
	return TRUE;
}

void SetupOGL(void)
{
	Int2 win_size;
	Uint2 win_usize;
	
	win_size = (Int2) MIN(screenRect.right, screenRect.bottom);

	win_size -= 128;
	
	if (win_size < 200)  win_size = 200;
            	else if (win_size > 800) win_size = 750;            	
            	win_usize = (Uint2) win_size;         		
            	myOGL_data =
                Vtrj_OGL_CreateViewer(main_win, &win_usize, win_usize,  X_ROTATE_SBAR |
                	Y_ROTATE_SBAR /*|  Z_ROTATE_SBAR*/, NULL, NULL,  Cn3D_InitMA, NULL);            		

	Vtrj_ColorData.OGL_Data = myOGL_data; 	
	RealizeWindow(main_win);

#ifdef WIN_MOTIF
    /* now that all windows are realized, set X OpenGL context */
    Vtrj_ColorData.OGL_Data->display = Nlm_currentXDisplay;

    Nlm_SetOGLContext(NULL, NULL, NULL, NULL);
    currentCtx = glXGetCurrentContext();
    currentXdrw = glXGetCurrentDrawable();
#endif

#ifdef WIN_MSWIN
         current_hglrc = wglGetCurrentContext();
         current_hdc = wglGetCurrentDC();

#endif /*WIN_MSWIN*/
	
/*    Cn3D_Redraw(TRUE);
*/
    SetResize(main_win, Vtrj_ResizeProc);
    Vtrj_ResizeProc(main_win);
}

void Vtrj_PostPicture(void)
{
#ifdef WIN_MOTIF
	glXSwapBuffers(Nlm_currentXDisplay, currentXdrw);
#endif /*WIN_MOTIF*/

#ifdef WIN_MSWIN
               SwapBuffers(wglGetCurrentDC());
#endif /*WIN_MSWIN*/

	glDisable(GL_NORMALIZE);
/*	glDisable(GL_BLEND);*/
	glPopMatrix();
}

void Vtrj_ClearScreen(void)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#ifdef WIN_MOTIF	
	glXSwapBuffers(Nlm_currentXDisplay, currentXdrw);
#endif /*WIN_MOTIF*/

#ifdef WIN_MSWIN
               SwapBuffers(wglGetCurrentDC());
#endif /*WIN_MSWIN*/

}

/*

$Log: vis_trajgl.c,v $
Revision 1.2  2001/03/04 19:43:30  feldman
fixed compiler warnings

Revision 1.1  2001/02/15 20:43:05  feldman
Split vis_traj into two files, to separated MMDBAPI from X stuff,
and added minor tweaks, progress monitor and bugfixes

*/

