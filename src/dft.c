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

#include <dft.h>

/* dft
 * discrete fourier transform
 *
 *
 */
void dft(ValNodePtr vnp, ValNodePtr PNTR transform_vnpp) {
  Int4 numpoints;
  Int4 n, k;
  complexptr curptr = NULL;
  complexptr curFptr = NULL;
  ValNodePtr curnode = NULL;

  numpoints = ValNodeLen(vnp);
  /*  printf("there are %d points\n", numpoints); */

  if (numpoints > 0) {

  /* for each frequency component */
  for (n = 0; n < numpoints; n++) { 
    curFptr = MemNew(sizeof(complex));
    curFptr->real = 0; curFptr->imag = 0;
    
    /* start again at beginning of list of points*/
    curnode=vnp;
    /* take sum over all points */
    for (k = 0; k < numpoints; k++) { 

    curptr = curnode->data.ptrvalue;

    curFptr->real +=
      (1.0 / (FloatHi) numpoints) * 
      ((curptr->real * cos(2.0 * PI * (FloatHi) n * (FloatHi) k / (FloatHi) numpoints)) -
	(curptr->imag * sin(2.0 * PI * (FloatHi) n * (FloatHi) k / (FloatHi) numpoints)));

    curFptr->imag += 
      (1.0 / (FloatHi) numpoints) * 
      ((curptr->imag * cos(2.0 * PI * (FloatHi) n * (FloatHi) k / (FloatHi) numpoints)) + 
       (curptr->real * sin(2.0 * PI * (FloatHi) n * (FloatHi) k / (FloatHi) numpoints)));			     

    curnode = curnode->next;
    }
    
    ValNodeAddPointer(transform_vnpp, 1, curFptr);
    /*    printf("freq: %f real: %8f, imag: %8f \n",((FloatLo) n / (FloatLo) numpoints), curFptr->real, curFptr->imag);*/
  }
  } 
} /* dft */



/* idft
 * 
 * inverse discrete fourier transform
 *
 */
void idft(ValNodePtr freqdom_vnp, ValNodePtr PNTR timedom_vnpp) {
  Int4 numpoints;
  Int4 n, k;
  complexptr curptr = NULL;
  complexptr curFptr = NULL;
  ValNodePtr curFnode = NULL; 
  
  numpoints = ValNodeLen(freqdom_vnp);
  curFnode = ValNodeFindNext(freqdom_vnp, curFnode, 1);

  if (numpoints > 0) {
 /* create as many data points as components */
  for (k = 0; k < numpoints; k++) { 
    curptr = MemNew(sizeof(complex));
    curptr->real = 0; curptr->imag = 0;

    /* start at first frequency */
    curFnode = freqdom_vnp;      

    /* sum over all frequencies */
    for (n = 0; n < numpoints; n++) { 
      curFptr = curFnode->data.ptrvalue;

      curptr->real +=
	((curFptr->real * cos(2.0 * PI * n * k / numpoints)) +
	 (curFptr->imag * sin(2.0 * PI * n * k / numpoints)));

      curptr->imag += 
	((curFptr->imag * cos(2.0 * PI * n * k / numpoints)) - 
	 (curFptr->real * sin(2.0 * PI * n * k / numpoints)));	
		     
      curFnode = curFnode->next;
    } /* sum over frequencies */
    
    ValNodeAddPointer(timedom_vnpp, 0, curptr);
    /* printf("IFT real: %f, imag: %f \n", curptr->real, curptr->imag); */ 
    
  } /* data points */
  } 
}

/* lowPassFilter 
 *
 * filters out all frequencies HIGHER than the cutoff
 * frequencies are implicit
 * frequency ranges from 0 to 1, in steps of n/N
 * 
 *  so we cut out frequencies between cutoff and 1-cutoff
 * 
 */
void lowPassFilter(ValNodePtr FDom_vnp, ValNodePtr PNTR FDomFiltered_vnpp, FloatLo cutoff_fh) {
  Int2 numpoints = 0;  
  Int2 n = 0;
  complexptr curptr = NULL;

  if (*FDomFiltered_vnpp != NULL) {
    /*   ErrPostEX(SEV_INFO,"Overwriting existing transform."); */
   *FDomFiltered_vnpp = NodeListFree(*FDomFiltered_vnpp);
  }

  numpoints = ValNodeLen(FDom_vnp);

  for(n=0; n < numpoints; n++) {  
   curptr = MemNew(sizeof(complex));
   curptr->real = 0; curptr->imag = 0;

   if ((((FloatHi) n / (FloatHi) numpoints) < cutoff_fh) ||
       (((FloatHi) n / (FloatHi) numpoints) > (1 - cutoff_fh))) {

     curptr->real = ((complexptr) (FDom_vnp->data.ptrvalue))->real;
     curptr->imag = ((complexptr) (FDom_vnp->data.ptrvalue))->imag;     
   }

   else {
     curptr->real = 0.0;
     curptr->imag = 0.0;
   }
   
   ValNodeAddPointer(FDomFiltered_vnpp, 1, curptr);
   FDom_vnp = FDom_vnp->next;
     }
}



/*
$Log: dft.c,v $
Revision 1.5  2001/08/15 20:15:13  feldman
Changed Fourier transform to opposite sign

Revision 1.4  2001/08/03 20:52:15  feldman
Fixed crashing bug and output

Revision 1.3  2000/08/18 19:06:05  adrian
moved r14seqfilecalc to r14bscalc and moved back into r14.c where it belongs
fixed sequence numbering off-by one in output (still needs to be checked elsewhere)
more debug information added to RLEunpack error message in trajtools

Revision 1.2  2000/06/22 19:25:44  feldman
Fixed remaining warnings

Revision 1.1  2000/06/20 18:59:46  feldman
moved to library directory and changed include files appropriately

Revision 1.2  2000/06/20 17:31:54  adrian
added license headers
 dft
*/
