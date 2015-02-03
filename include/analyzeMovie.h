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

#ifndef ANALYZE_H
#define ANALYZE_H

/* includes */
#include <hfprogs.h>
#include <gifgen.h>
/*include the functions retrieving the radius of gyration*/
#include <mmdbtraj_pub.h>
#include <slriaccsurf.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

#define NUMARGS			5
#define NUM_OF_COLORS	10			/*number of colors availaible to the 
									graphs and defined in define color */
#define NUM_COLOR_RAINBOW  10		/*Number of colors in the gradient used to 
									 represent the diffferent kind of gradients*/
#define LENGTH_SUFFIX_MOVIE 12		/*length of the suffixe that we want to retrieve from the
									  name of the movie to obtain the name of the protein we want*/ 
#define NUM_MODEL		/*100*/ MAXGEN

/*modes for graphing in graphArray*/
#define MODE_MODELTIME	1

#define ARRAY_VS_TIME				1	/*graphing like: f(x) = ax+b for x in [c;d]*/			
#define DOT_ASA_VS_TIME				2	/*graph the ASA*/
#define DOT_2NDRY_STUC_VS_TIME		3	/*plot the 2ndry structures kinds*/
#define DOT_AVG_2NDRY_STUC_VS_TIME	5	/*plot the 2ndry struct. when averaged*/
#define CONTACT_MAP					6   
#define ENERGY_PLOTTING				7	/*plot residue number VS time with plot representing Energy level*/
#define RGYR_VS_RMSD				8	/*plot using the bryant energy level*/
#define ARRAY_VS_TIME_FLOAT			9	/*plot the vertical scale using float numbers*/
#define ARRAY_VS_RESIDUE			10  /*when plotting data per residue*/
#define	RECTANGLE_VS_CONTACT		11  /*trace rectangles which height is proportionnal to the occupancy*/

#define NUM_KIND_OF_2NDRY_STRUCT	7
#define NUM_GROUP_2NDRY_STRUCT		3	/*helix, sheet, (coil and others)*/ 
#define RGYR_VS_RMSD_SQUARE_SIZE	20	/*plot using the bryant energy level*/


/*Some useful functions to travel arrays of int, float, etc ... and print their 
values to the screen as well as to a file with row tab seperated*/
/*to print the radius of gyration using hydrophobic forces only or not*/
void PrintRGyr(FloatLo * flRadiusHP, FloatLo * flRadiusNoHP, Int4 nbOfModel, CharPtr filename);
/*to print only a one dimension array of floats*/
void Printfl(FloatLo * flArray, Int4 nbOfModel, CharPtr data, CharPtr filename);
/*to print a 2 dimension array of floats*/
void Printppfl(FloatLo ** pflToPrint, Int2 xaxis, Int2 yaxis, CharPtr xdata, CharPtr ydata, CharPtr filename);
/*to print a 3 dimension array of floats*/
void Printpppfl(FloatLo *** pppflToPrint,Int4 iNumRes,Int2 numPartition,CharPtr filename);

/*to print only a one dimension array of integers*/
void Printpi4(Int4 * pi4ToPrint, Int4 nbOfModel, CharPtr data, CharPtr filename);
/*to print a 2 dimension array of  integers*/
void Printppi4(Int4 ** pi4ToPrint, Int2 xaxis, Int2 yaxis,  CharPtr xdata, CharPtr ydata,CharPtr filename);
/*to print a 3 dimension array of  integers*/
void PrintppChar(CharPtr * arrayCharToPrint, Int2 xaxis, Int2 yaxis,  CharPtr xdata, CharPtr ydata, CharPtr filename);

/*given arrays of float or ints, returns the largest/smallest number*/
Int4 findSmallestInt(Int4 iArray[], Int4 arraySize);
Int4 findGreatestInt(Int4 iArray[], Int4 arraySize);
Int4 findSmallestIntInppInt(Int4 ** iArray, Int4 numElmtDim1, Int4 numElmtDim2);
Int4 findGreatestIntInppInt(Int4 ** iArray, Int4 numElmtDim1, Int4 numElmtDim2);
FloatLo findSmallestFloat(FloatLo floatArray[], Int2 arraySize);
FloatLo findGreatestFloat(FloatLo floatArray[], Int2 arraySize);
FloatLo findSmallestFloatInppFloat(FloatLo ** flArray, Int4 numElmtDim1, Int4 numElmtDim2);
FloatLo findGreatestFloatInppFloat(FloatLo ** flArray, Int4 numElmtDim1, Int4 numElmtDim2);

/*Given a PMAD, returns if the atom has an aromatic bond or not*/
Boolean IsAromatic(PMAD pmadThis);

/*GraphArray:
arg. description
  -CharPtr gifFileName: name given to the gif file generated by the function
  -FloatLo * flArray1 & FloatLo * flArray2: Array of data we want to plot in mode ARRAY_VS_TIME
											Set to NULL if not in this mode and not needed.
  -Int2 numModels:	Number of models to plot (on x axis)
						--> good in modes: basically all the mode .._VS_TIME and ENERGY_PLOTTING
					in CONTACT_MAP  mode: represent the partition number that is drawn 
					in RGYR_VS_RMSD mode: represents the maximum RMSD (which is on the x axis) 
					in ARRAY_VS_RESIDUE MODE: represents the number of residue 

  -Int2 graphNum :	 set to "1" if you want to generate a new picture
				else set to "n" (n>0) and this is the nth addition to the gif pictures
  -FloatLo flLow1, FloatLo flHigh1: determies the limits on the y axis, and determines 
									the plotting range on the y axis
									Note: in DOT_RESIDUE_VS_TIME mode, 
											flHigh1 is the number of residue in the protein
  -CharPtr pcLegend1, CharPtr pcLegend2:legend corresponding to flArray1 & flArray2
  -CharPtr pcyAxis, CharPtr pcxAxis: Title of each of those axis.
  -Int2 xSize, Int2 ySize:	xSize determines the Width of the gif picture.
							ySize determines the height of the gif picture.
  -Int2 tracingMode: determines the tracing mode.
  -void * option1: in DOT_RESIDUE_VS_TIME should containt the 2-d array containing the wanted informations.
*/
void GraphArray(CharPtr gifFileName,FloatLo * flArray1 ,FloatLo * flArray2, Int2 numModels,Int2 graphNum, 
				FloatLo flLow1, FloatLo flHigh1, CharPtr pcLegend1, CharPtr pcLegend2,CharPtr pcyAxis, CharPtr pcxAxis, 
				Int2 xSize, Int2 ySize , Int2 tracingMode, VoidPtr option1);



/*definition of contacts*/
#define MAX_DIST_CONTACT_ALIPH_C	5.4	/*maximum distance to considerate a contact with an aliphatic C*/
#define MAX_DIST_CONTACT_OTHER		4.6	/*maximum distance to considerate a contact for other atoms*/
#define NUM_DIV_IN_SIMULATION		4	/*Number of partitions of the contact map that we want to see
										  eg: if NUM_DIV_IN_SIMULATION = 4 and we have 100 frames, 
										  then we will compute the contacts on 1-25, 25-50, 50-75, 75-100*/
#define ZHANG_WINDOWSIZE			4
#define BRYANT_WINDOWSIZE			4


/*data used for the NOEs*/
#define NOE_CUTOFF					10	/* distance in amstrong below which 2 atoms are considered as 
										   being in contact when computing the NOE's */
#define NUM_CONTACTS_IN_DRK			51 /* 59 */

/*function added to complete/replace the function gdImageGif.
  Correct some uninitialyzed data.
*/
static void advGdImageGif(gdImagePtr  imagePtr, FILE * fileOut);


#ifdef __cplusplus
}
#endif

#endif

