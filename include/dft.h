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

#ifndef DFT_H
#define DFT_H

#include <ncbi.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PI
#define PI 3.14159
#endif

typedef struct nlm_complex {
  FloatHi real;
  FloatHi imag;
} complex, PNTR complexptr;

/* sequences of complex numbers should be stored in ValNodes
   set the choice to 0 if the data is in the time or spatial domains
   set the choice to 1 if the data is in the frequency domain
*/

/* dft  -  discrete fourier transform
 *
 * Input:  A sequence of N complex datapoints in the time or space domain
 * Output: The discrete fourier transform of the input data,
 *         as a list of N complex frequency components
 *
 */
void dft(ValNodePtr vnp, ValNodePtr PNTR transform_vnpp);

/* idft - inverse discrete fourier transform
 *
 *
 * Input:  A series of frequency components.
 * Output: A series of points in the time or space domain.
 */
void idft(ValNodePtr freqdom_vnp, ValNodePtr PNTR timedom_vnp);

/* lowPassFilter
 * 
 * Input:  A series of frequency components.
 * Output: A series of frequency components, with all frequencies
 *         above the cutoff set to zero.
 */ 
void lowPassFilter(ValNodePtr FDom_vnp, ValNodePtr PNTR FDomFiltered_vnpp, FloatLo cutoff_fh);

/*
 * initialize a result array
 * returns an ID for the array
 */
Int2 arrayInit(Int4 width, Int4 height);

/* arrayWriteFromList
 * writes a valnode list into a column of the array
 *
 */
void arrayWriteFromList(ValNodePtr vnp, Int4 column);



/* arrayRead
 *
 * read a column from a file into an array
 *
 */
void arrayRead(CharPtr File, Int2 arrayID, Int2 column);

void arrayPrint(FloatHi PNTR ar, Int2 length);

/* arrayfromList
 * copy a sequence of numbers from a valnode list into a 
 * dynamically allocated array 
 * 
 * returns: length of the list that has been copied
 */
Int2 arrayfromList(ValNodePtr vnp, FloatHi PNTR PNTR araddress);

#ifdef __cplusplus
}
#endif

#endif 

/*
$Log: dft.h,v $
Revision 1.6  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.5  2000/08/14 20:31:47  feldman
Corrected complex structure


*/


