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
#include "analyzeMovie.h"

/*	Global Variables */
Args Rargs[NUMARGS] =  
	{{"Input MMDB filename (enter name of the protein (PATH\\X); to open multiple files, the filenames must be of the format PATH\\X_movie_#.val)",NULL,NULL,NULL,FALSE,'f',ARG_FILE_IN,0.0,0,NULL},
	{"Native structure filename",NULL,NULL,NULL,TRUE,'g',ARG_FILE_IN,0.0,0,NULL},
	{"# files to process","1","1","20",TRUE,'n',ARG_INT,0.0,0,NULL},
	{"Average data (if no, all data from separate simulations are plotted on the same graph)?","FALSE",NULL,NULL,TRUE,'a',ARG_BOOLEAN,0.0,0,NULL},
	{"# partitions of total number of frames when making contact maps","4","1","50",TRUE,'p',ARG_INT,0.0,0,NULL}
/*	{"Width of the graphs","700","100","50000",TRUE,'w',ARG_INT,0.0,0,NULL},
	{"Height of the graphs","520","100","50000",TRUE,'h',ARG_INT,0.0,0,NULL} */
	};


/*Given two arrays of radius of gyration, prints the radius of gyration */
void PrintRGyr(FloatLo * flRadiusHP, FloatLo * flRadiusNoHP, Int4 nbOfModel, CharPtr filename)
{
	Int4 model;
	FILE * fp;
	fp = FileOpen(filename,"w");
/*	printf("\n\nModel #\t| rGyrHP\t| rGyr\t");*/
	fprintf(fp,"\n\nModel #\trGyrHP\trGyr\t");

	for(model = 1; model <= nbOfModel ; model ++){
/*		printf("\n%d\t|%3.3f\t\t|%3.3f\t\t|", (int)model, (float)flRadiusHP[model -1], (float)flRadiusNoHP[model -1]);*/
		fprintf(fp,"\n%d\t%3.3f\t%3.3f\t", (int)model, (float)flRadiusHP[model -1], (float)flRadiusNoHP[model -1]);
	}
	FileClose(fp);
}

/*Given an array of rmsd, prints the RMSD */
void Printfl(FloatLo * flArray, Int4 nbOfModel, CharPtr data, CharPtr filename)
{
	Int2 model;
	FILE * fp;
	fp = FileOpen(filename,"w");

/*	printf("\n\nModel #\t| %s\t|", data);*/
	fprintf(fp,"\n\nModel #\t%s\t", data);
	for(model = 1; model <= nbOfModel ; model ++){
/*		printf("\n%d\t|%3.3f\t|", (int)model, (float)flArray[model -1]);*/
		fprintf(fp,"\n%d\t%3.3f\t", (int)model, (float)flArray[model -1]);
	}
	FileClose(fp);
}

void Printfl3(FloatLo * flArray[3], Int4 nbOfModel, Char data[3][25], CharPtr filename)
{
	Int2 model;
	FILE * fp;
	fp = FileOpen(filename,"w");

/*	printf("\n\nModel #\t| %s\t| %s\t| %s\t|", data[0],data[1],data[2]);*/
	fprintf(fp,"\n\nModel #\t%s\t%s\t%s\t", data[0],data[1],data[2]);
	for(model = 1; model <= nbOfModel ; model ++){
/*		printf("\n%d\t|%3.3f\t|%3.3f\t|%3.3f\t|", (int)model, (float)flArray[0][model -1], (float)flArray[1][model -1], (float)flArray[2][model -1]);*/
		fprintf(fp,"\n%d\t%3.3f\t%3.3f\t%3.3f\t", (int)model, (float)flArray[0][model -1], (float)flArray[1][model -1], (float)flArray[2][model -1]);
	}
	FileClose(fp);
}

/*Given a 2-d array of Float */
void Printppfl(FloatLo ** pflToPrint, Int2 xaxis, Int2 yaxis, CharPtr xdata, CharPtr ydata, CharPtr filename)
{
	Int4 i,j;
	FILE * fp;
	fp = FileOpen(filename,"w");
/*	printf(    "\n\n%d %s(yaxis) \\ %d %s(xaxis) \n",yaxis,ydata,xaxis, xdata);*/
	fprintf(fp,"\n\n%d %s(yaxis) \\ %d %s(xaxis) \n",yaxis,ydata,xaxis, xdata);
	/*x legend*/
/*	printf("|\t|");*/
	fprintf(fp,"\t");
	for (i=1; i<=xaxis; i++){
/*		printf("%d\t|",i);*/
		fprintf(fp,"%d\t|",i);
	}
	for (i = 0; i < yaxis ; i ++){
/*		printf(    "\n%d\t|",i+1);*/
		fprintf(fp,"\n%d\t" ,i+1);
		for(j =0; j < xaxis; j++){
/*			printf(    "%5.3f\t|", (float)pflToPrint[j][i]);*/
			fprintf(fp,"%5.3f\t",  (float)pflToPrint[j][i]);
		}
	}
	FileClose(fp);
}

/*Given an array of int4 */
void Printpi4(Int4 * pi4ToPrint, Int4 nbOfModel, CharPtr data, CharPtr filename)
{
	Int4 model;
	FILE * fp;
	fp = FileOpen(filename,"w");

/*	printf("\n\nModel #\t| %s\t|", data);*/
	for(model = 1; model <= nbOfModel ; model ++){
/*		printf ("\n%d\t|%d\t|", (int)model, (int)pi4ToPrint[model -1]);*/
		fprintf(fp,"\n%d\t%d\t", (int)model, (int)pi4ToPrint[model -1]);
	}

	FileClose(fp);
}

/*Given an array of int4 */
void Printppi4(Int4 ** pi4ToPrint, Int2 xaxis, Int2 yaxis,  CharPtr xdata, CharPtr ydata,CharPtr filename)
{
	Int4 i,j;
	FILE * fp;
	fp = FileOpen(filename,"w");
/*	printf(    "\n\n%d %s(yaxis) \\ %d %s(xaxis) \n",yaxis,ydata,xaxis, xdata);*/
	fprintf(fp,"\n\n%d %s(yaxis) \\ %d %s(xaxis) \n",yaxis,ydata,xaxis, xdata);
	/*x legend*/
/*	printf("|\t|");*/
	fprintf(fp,"\t");
	for (i=1; i<=xaxis; i++){
/*		printf("%d\t|",i);*/
		fprintf(fp,"%d\t|",i);		
	}
	for (i = 0; i < yaxis ; i ++){
/*		printf(    "\n%d\t|",i+1);*/
		fprintf(fp,"\n%d\t",i+1);
		for(j =0; j < xaxis; j++){
/*			printf(    "%d\t|", (int)pi4ToPrint[j][i]);*/
			fprintf(fp,"%d\t", (int)pi4ToPrint[j][i]);
		}
	}
	FileClose(fp);
}


/*Given an array of int4 */
void PrintppChar(CharPtr * arrayCharToPrint, Int2 xaxis, Int2 yaxis,  CharPtr xdata, CharPtr ydata, CharPtr filename)
{
	Int4 i,j;
	FILE * fp;
	fp = FileOpen(filename,"w");
/*	printf("\n\n%d %s(yaxis) \\ %d %s(xaxis) \n",(int)yaxis,ydata,(int)xaxis, xdata);*/
	fprintf(fp,"\n\n%d %s(yaxis) \\ %d %s(xaxis) \n",(int)yaxis,ydata,(int)xaxis, xdata);
	/*x legend*/
/*	printf("|\t|");*/
	fprintf(fp,"\t");
	for (i=1; i<=xaxis; i++){
/*		printf("%ld\t|",(long)i); */
		fprintf(fp,"%ld\t",(long)i); 
	}
	for (i = 0; i < yaxis ; i ++){
/*		printf("\n%ld\t|",(long)i+1);*/
		fprintf(fp,"\n%ld\t",(long)i+1);
		for(j =0; j < xaxis; j++){
/*			printf("%c\t|", arrayCharToPrint[j][i]);*/
			fprintf(fp,"%c\t", arrayCharToPrint[j][i]);
		}
	}
	FileClose(fp);
}

void Printpppfl(FloatLo *** pppflToPrint,Int4 iNumRes,Int2 numPartition,CharPtr filename){
	Int2 i,j,k;

	FILE * fp;
	fp = FileOpen(filename,"w");

/*	printf("\t|"); */
	fprintf(fp,"\t");

	for (i=1; i<=iNumRes; i++){
/*		printf("%d\t|",(int)i); */
		fprintf(fp,"%d\t",(int)i);
	}

	for(i=0; i<numPartition; i++){
/*		printf("\n Partition: %d\n", i);*/
		fprintf(fp,"\n\nPartition: %d %s\n",i,i==0?"(native)":"");
		for(j=0; j< iNumRes; j++){
/*			printf("\n %d\t|",(int)j+1);*/
			fprintf(fp,"\n %d\t",(int) j+1);
			for(k=0; k<=j; k++){
/*				printf("%2.4f\t|",    (float)pppflToPrint[i][j][k]);*/
				fprintf(fp,"%2.4f\t", (float)pppflToPrint[i][j][k]);
			}
		}
	}
	FileClose(fp);
}

/*given an array of floats, returns the smallest number*/
Int4 findSmallestInt(Int4 iArray[], Int4 arraySize)
{
	Int4	iMinSoFar;
	Int4	cnt;
	iMinSoFar = iArray[0];
	for(cnt=1; cnt < arraySize; cnt ++)
	{
		if( (iArray[cnt] < iMinSoFar) )
		{
			iMinSoFar = iArray[cnt];
		}
	}
	return iMinSoFar;
}
/*given an array of ints, returns the greatest number*/
Int4 findGreatestInt(Int4 iArray[], Int4 arraySize)
{
	Int4 iMaxSoFar;
	Int4 cnt;
	iMaxSoFar = iArray[0];
	for(cnt=1; cnt < arraySize; cnt ++)
	{
		if(iArray[cnt] > iMaxSoFar)
		{
			iMaxSoFar = iArray[cnt];
		}
	}
	return iMaxSoFar;
}

/*given 2-d array of ints, returns the minimum Int */
/*given an array of Ints, returns the smallest number*/
Int4 findSmallestIntInppInt(Int4 ** iArray, Int4 numElmtDim1, Int4 numElmtDim2)
{
	Int4	iMinSoFar;
	Int4	cnt, temp;
	iMinSoFar = findSmallestInt(iArray[0],numElmtDim2);
	for(cnt=1; cnt < numElmtDim1; cnt ++)
	{
		temp = findSmallestInt(iArray[cnt],numElmtDim2);
		if(  temp < iMinSoFar )
		{
			iMinSoFar = temp;
		}
	}
	return iMinSoFar;
}

/*given 2-d array of ints, returns the maximum Int */
/*given an array of ints, returns the smallest number*/
Int4 findGreatestIntInppInt(Int4 ** iArray, Int4 numElmtDim1, Int4 numElmtDim2)
{
	Int4	iMaxSoFar;
	Int4	cnt, temp;
	iMaxSoFar = findGreatestInt(iArray[0],numElmtDim2);
	for(cnt=1; cnt < numElmtDim1; cnt ++)
	{
		temp = findGreatestInt(iArray[cnt],numElmtDim2);
		if(  temp > iMaxSoFar )
		{
			iMaxSoFar = temp;
		}
	}
	return iMaxSoFar;
}

/*given an array of floats, returns the smallest number*/
FloatLo findSmallestFloat(FloatLo floatArray[], Int2 arraySize)
{
	FloatLo flMinSoFar;
	Int2	cnt;
	flMinSoFar = floatArray[0];
	for(cnt=1; cnt < arraySize; cnt ++)
	{
		if( (floatArray[cnt] < flMinSoFar) )
		{
			flMinSoFar = floatArray[cnt];
		}
	}
	return flMinSoFar;
}

/*given an array of floats, returns the greatest number*/
FloatLo findGreatestFloat(FloatLo floatArray[], Int2 arraySize)
{
	FloatLo flMaxSoFar;
	Int2	cnt;
	flMaxSoFar = floatArray[0];
	for(cnt=1; cnt < arraySize; cnt ++)
	{
		if(floatArray[cnt] > flMaxSoFar)
		{
			flMaxSoFar = floatArray[cnt];
		}
	}
	return flMaxSoFar;
}

/*given 2-d array of ints, returns the minimum Int */
/*given an array of Ints, returns the smallest number*/
FloatLo findSmallestFloatInppFloat(FloatLo ** flArray, Int4 numElmtDim1, Int4 numElmtDim2)
{
	FloatLo flMinSoFar, temp;
	Int4	cnt;
	flMinSoFar = findSmallestFloat(flArray[0],numElmtDim2);
	for(cnt=1; cnt < numElmtDim1; cnt ++)
	{
		temp = findSmallestFloat(flArray[cnt],numElmtDim2);
		if(  temp < flMinSoFar )
		{
			flMinSoFar = temp;
		}
	}
	return flMinSoFar;
}

/*given 2-d array of ints, returns the maximum Int */
/*given an array of ints, returns the smallest number*/
FloatLo findGreatestFloatInppFloat(FloatLo ** flArray, Int4 numElmtDim1, Int4 numElmtDim2)
{
	FloatLo	flMaxSoFar, temp;
	Int4	cnt;
	flMaxSoFar = findGreatestFloat(flArray[0],numElmtDim2);
	for(cnt=1; cnt < numElmtDim1; cnt ++)
	{
		temp = findGreatestFloat(flArray[cnt],numElmtDim2);
		if(  temp > flMaxSoFar )
		{
			flMaxSoFar = temp;
		}
	}
	return flMaxSoFar;
}

/*note: update the number of colors in analyzeUnfold.h when adding colors*/
void defineColor(gdImagePtr imagePtr, Int4 * arrayOfColor)
{
	Int4	/*white,*/ black , sky , red, green, lblue , dblue, redMix;
	Int4    greenLight,redLight,greenDark;

	/*define the colors*/
	/*White allocated first to set the background color*/
	gdImageColorAllocate(imagePtr, 255, 255, 255);

    redMix   = gdImageColorAllocate(imagePtr, 255, 145, 104);
	arrayOfColor[0] = redMix;

	black = gdImageColorAllocate(imagePtr,   0,   0,   0);
	arrayOfColor[1] = black;

    red   = gdImageColorAllocate(imagePtr, 255, 0, 0);
	arrayOfColor[2] = red;

	sky   = gdImageColorAllocate(imagePtr,   0, 200, 170);
	arrayOfColor[3] = sky;

	greenLight= gdImageColorAllocate(imagePtr,   190, 232,148);
	arrayOfColor[4] = greenLight;

	lblue = gdImageColorAllocate(imagePtr,  51, 102, 255);     
	arrayOfColor[5] = lblue;

	green = gdImageColorAllocate(imagePtr, 0, 155, 0);
	arrayOfColor[6] = green;

    redLight = gdImageColorAllocate(imagePtr, 225,117,117);
	arrayOfColor[7] = redLight;

    greenDark = gdImageColorAllocate(imagePtr, 106, 170, 17);
	arrayOfColor[8] = greenDark;

	dblue = gdImageColorAllocate(imagePtr, 0, 0, 100);
	arrayOfColor[NUM_OF_COLORS -1] = dblue;

	return;
}

void defineRainbowColor(gdImagePtr imagePtr, Int4 * arrayOfColor)
{

	Int4 orange,greenLight,green,greenDark,blue;
	Int4 blueDark,blueDarkest,redDark, redLight	,red,yellow;
	
	redDark= gdImageColorAllocate(imagePtr, 161,36,50);
	arrayOfColor[0] = redDark;

	red= gdImageColorAllocate(imagePtr, 227,64,33);
	arrayOfColor[1] = red;

	redLight= gdImageColorAllocate(imagePtr, 239,141,123);
	arrayOfColor[2] = redLight;

	orange   = gdImageColorAllocate(imagePtr, 232,152,49);
	arrayOfColor[3] = orange;

    yellow = gdImageColorAllocate(imagePtr, 237,239,15);
	arrayOfColor[4] = yellow;

	greenLight = gdImageColorAllocate(imagePtr,104,233,153);
	arrayOfColor[5] = greenLight;

	green = gdImageColorAllocate(imagePtr,51,192,104);     
	arrayOfColor[6] = green;

	greenDark   = gdImageColorAllocate(imagePtr, 17,150,67);
	arrayOfColor[7] = greenDark;

    blue = gdImageColorAllocate(imagePtr, 102,159,212);
	arrayOfColor[8] = blue;

	blueDark = gdImageColorAllocate(imagePtr, 27,91,150);
	arrayOfColor[9] = blueDark;

	blueDarkest = gdImageColorAllocate(imagePtr, 27,91,148);
	arrayOfColor[10] = blueDarkest;
}

void calculatePercent(FloatLo * arrayToStore, CharPtr * arrayOf2ndryStruct, CharPtr strucToLookFor, Int2 numModels, Int2 numResidues){
	FloatLo temp;
	Int2 i,j;
	for (i=0; i<numModels; i++){
		temp = 0.0;
		for(j=0; j< numResidues; j++){
			if(StringChr(strucToLookFor, arrayOf2ndryStruct[i][j]) != NULL){
				temp ++; }
		}
		arrayToStore[i]= (temp / numResidues) * 100;
	}
}

void Calc2ndryStrucRepart(FloatLo ** ppfl2ndryStruct, CharPtr *ppc2ndrystruct, Int4 numModels, Int4 numRes)
{
		calculatePercent( ppfl2ndryStruct[0] ,ppc2ndrystruct, "HI", numModels, numRes);
		calculatePercent( ppfl2ndryStruct[1] ,ppc2ndrystruct, "BE", numModels, numRes);
		calculatePercent( ppfl2ndryStruct[2] ,ppc2ndrystruct, "TSG ", numModels, numRes);
		/*calculatePercent( ppfl2ndryStruct[3] ,ppc2ndrystruct, " ",  numModels, numRes);*/
}

void Calc2ndryStrucRepartAvg(FloatLo ** ppfl2ndryStruct, CharPtr * ppc2ndryStruct,
							 CharPtr ** pppc2ndryStruct, Int4 numModels, Int4 numRes, Int4 numFiles)
{
	FloatLo ** temp;	/*Array that will temporarly contain the percents for calculation of the average*/
	Int2 i,j,k;
	Int4 *** pppi4ComputAvgStruct; /*array to store the predominant strucuture when averaging*/

	/*allocate memory for the temp array*/
	temp = (FloatLo**) MemNew( sizeof(FloatLo) * NUM_GROUP_2NDRY_STRUCT);
	for(i=0; i<NUM_GROUP_2NDRY_STRUCT; i++){
		temp[i] = (FloatLo*) MemNew( numModels * sizeof(FloatLo) );
	}

	/*allocate memory for pppi4ComputAvgStruct*/
	pppi4ComputAvgStruct = (Int4***) MemNew(sizeof(Int4*) * numModels);
	for (i=0; i< numModels; i++){
		pppi4ComputAvgStruct[i] = (Int4**) MemNew(sizeof(Int4*) * numRes);
		for(j =0 ; j< numRes; j++){
			pppi4ComputAvgStruct[i][j] = (Int4 *) MemNew(sizeof(Int4*) * NUM_GROUP_2NDRY_STRUCT);
		}
	}


	/*COMPUTING OF ppfl2ndryStruct*/
	/*compute the average of the secondary structures*/
	for(i=0; i<numFiles; i++){
		calculatePercent( temp[0] ,pppc2ndryStruct[i], "HI", numModels, numRes);
		calculatePercent( temp[1] ,pppc2ndryStruct[i], "BE", numModels, numRes);
		calculatePercent( temp[2] ,pppc2ndryStruct[i], "TSG ", numModels, numRes);

		/*"merge the result from temp with ppfl2ndrystruct*/
		for(j=0; j<NUM_GROUP_2NDRY_STRUCT; j++){
			for (k=0; k< numModels; k++){
				ppfl2ndryStruct[j][k] = (	(temp[j][k] / (FloatLo)(i+1)) +  
											(ppfl2ndryStruct[j][k]* i /(i+1) ))  ;
			}
		}
	}

	/*COMPUTING OF ppc2ndryStruct*/
	for (i=0; i< numFiles; i++){
		for(j =0 ; j< numModels; j++){
			for (k=0; k< numRes; k++){
				if (pppc2ndryStruct[i][j][k] == 'H' || pppc2ndryStruct[i][j][k] == 'I') 
					pppi4ComputAvgStruct[j][k][0] ++;
				else if (pppc2ndryStruct[i][j][k] == 'B' || pppc2ndryStruct[i][j][k] == 'E') 
					pppi4ComputAvgStruct[j][k][1] ++;
				else if (pppc2ndryStruct[i][j][k] == 'T' || pppc2ndryStruct[i][j][k] == 'S' 
					     ||pppc2ndryStruct[i][j][k] == 'G' || pppc2ndryStruct[i][j][k] == ' ' ) 
					pppi4ComputAvgStruct[j][k][2] ++;
			}
		}
	}

	/*fill the 2ndry structure average array*/
	for(i=0; i < numModels; i++){
		for(j=0; j< numRes; j++){
			for (k=0; k< NUM_GROUP_2NDRY_STRUCT; k++){
				if ( ((FloatLo)pppi4ComputAvgStruct[i][j][k]/(FloatLo)numFiles) > 0.5 ){
					if(k==0) 
						ppc2ndryStruct[i][j] = 'H';
					else if (k==1) 
						ppc2ndryStruct[i][j] = 'E';						
					else if (k==2) 
						ppc2ndryStruct[i][j] = 'C';													
				}
			}
			if (ppc2ndryStruct[i][j] == '\0')
				ppc2ndryStruct[i][j] = 'C';	
		}
	}

	for (i=0; i< numModels; i++){
		for(j =0 ; j< numRes; j++){
			MemFree(pppi4ComputAvgStruct[i][j]);
		}
		MemFree(pppi4ComputAvgStruct[i] );
	}
	MemFree( pppi4ComputAvgStruct);

	for(i=0; i<NUM_GROUP_2NDRY_STRUCT; i++){
		MemFree(temp[i]);
	}
	MemFree(temp);
}

static void advGdImageGif(gdImagePtr  imagePtr, FILE * fileOut){
	Int2 i;
	for(i = (Int2)imagePtr->colorsTotal; i< (Int2)gdMaxColors; i++){
		imagePtr->red[i]=0;
		imagePtr->green[i]=0;
		imagePtr->blue[i]=0;
	}
	gdImageGif(imagePtr, fileOut);
}

/*function used to graph the radius of gyration in function of the time*/
void GraphArray(CharPtr gifFileName,FloatLo * flArray1 ,FloatLo * flArray2, Int2 numModels,Int2 graphNum,
				FloatLo flLow1, FloatLo flHigh1, CharPtr pcLegend1, CharPtr pcLegend2,CharPtr pcyAxis, CharPtr pcxAxis,
				Int2 xSize, Int2 ySize , Int2 tracingMode, VoidPtr option1)
{
	gdImagePtr	imagePtr;
	FILE		*fileOut, *fileIn;
	Int4		arrayOfColor[NUM_OF_COLORS],arrayRainbowColors[NUM_COLOR_RAINBOW+1];
	Int2		cnt=0,i=0,j,iMinOnVScale, istep, minASA, maxASA,iMaxRmsd=0,cur,gap;
	FloatLo		xTopLeft, yTopLeft, xBotLeft, yBotLeft, xBotRight, yBotRight, xWide, yWide;
	FloatLo		xPrevious , yPrevious1=0.0,  yPrevious2=0.0, xUnit, yUnit,xTemp, yTemp, yMin, flStep;
	Char		pcSmall[50];
	Char		pcTemp[PATH_MAX];
	FloatLo		flLow, flHigh,graphMarginLow, graphMarginHigh, flvstep=0.0,flxstep=0.0, flvmargin;
	FloatLo		flMinEn, flMaxEn;

	graphMarginLow  = 0.0;
	graphMarginHigh = 0.0;

	flLow	= flLow1  - graphMarginLow;
	flHigh	= flHigh1 + graphMarginHigh;

	xTopLeft = (FloatLo)50.0 ;
	yTopLeft = (FloatLo)40.0 ;

	xBotLeft = (FloatLo)xTopLeft;
	yBotLeft = (FloatLo)ySize - 50;

	xBotRight = (FloatLo)xSize - (FloatLo)130 ;
	yBotRight = (FloatLo)yBotLeft ;

	xWide =	xBotRight -xBotLeft - 15;
	yWide = yBotLeft - yTopLeft - 15;

	yUnit	= (yWide / (flHigh - flLow));
	xUnit	= ((xWide) / NUM_MODEL);

	/*Create the file and trace legend and graph*/
	if( graphNum ==1 )
	{
		printf("\nCreating the gif: %s\n", gifFileName);
		imagePtr = gdImageCreate(xSize, ySize);
		/*define the colors*/
		defineColor(imagePtr, arrayOfColor);

		/*draw the bases of the diagram*/
		/*vertical bar + arrows*/
		gdImageLine(imagePtr, xTopLeft,  yBotLeft , xTopLeft	,  yTopLeft,	arrayOfColor[1]);
		gdImageLine(imagePtr, xTopLeft,  yTopLeft ,	xTopLeft +5	,  yTopLeft +5, arrayOfColor[1]);
		gdImageLine(imagePtr, xTopLeft,  yTopLeft ,	xTopLeft -5	,  yTopLeft +5, arrayOfColor[1]);

		/*horizontal bar + arrows*/
		gdImageLine(imagePtr, xBotLeft,  yBotLeft,  xBotRight,     yBotRight   , arrayOfColor[1]);
		gdImageLine(imagePtr, xBotRight, yBotRight, xBotRight -5,  yBotRight -5, arrayOfColor[1]);
		gdImageLine(imagePtr, xBotRight, yBotRight, xBotRight -5,  yBotRight +5, arrayOfColor[1]);

		/*write horizontal and vertical indices*/
		gdImageString(imagePtr, gdFont6X12,xBotRight + 15 , yBotRight, pcxAxis, arrayOfColor[1]);
		gdImageString(imagePtr, gdFont6X12,30, 25 , pcyAxis, arrayOfColor[1]);

		if( StringCmp(pcLegend1, "") != 0) {
			sprintf(pcTemp, "%s%s", "_", pcLegend1);
			gdImageString(imagePtr, gdFont6X12,10 , 5, pcTemp, arrayOfColor[2]);
		}
		if(  StringCmp(pcLegend2,"") != 0) {
			sprintf(pcTemp, "%s%s", "..", pcLegend2);
			gdImageString(imagePtr, gdFont6X12,10 , 15, pcTemp, arrayOfColor[2]);
		}

		/*draw the legend verticaly on the diagram and insert the indices*/
		iMinOnVScale = (Int4)flLow;		/*make sure that we start from an integer*/
		istep  = (Int4)( (flHigh - flLow) /10) ;
		if(istep==0){istep =1;}
		flvstep = yWide /(flHigh - flLow);
		flvmargin =istep -  (flLow - (FloatLo)iMinOnVScale);


		/*print intermediates values*/

		/*case where the range is larger than a single unit*/
		if( (abs(flHigh - flLow) ) > 1.00 ){

			/*print extreme values*/
			/*print the lowest value*/
			if ((fabs(flLow)-floor(fabs(flLow)>0.1)) && (fabs(flLow)-floor(fabs(flLow)<0.9)) && flLow>1.0 && flLow<99.9)
				sprintf(pcSmall,"%3.1f", flLow ) ;
			else
				sprintf(pcSmall,"%d", (int)flLow ) ;
			if((tracingMode == CONTACT_MAP ) || (tracingMode == DOT_ASA_VS_TIME) || (tracingMode == ENERGY_PLOTTING)
				|| (tracingMode ==DOT_2NDRY_STUC_VS_TIME) || (tracingMode == DOT_AVG_2NDRY_STUC_VS_TIME) )
				/*sprintf(pcSmall,"%d", (int)flLow +1) ;			*/;
				else {

			gdImageString(	imagePtr, gdFont6X12, 15 ,
								(Int2)(yBotLeft -5 ),
								pcSmall, arrayOfColor[1]);
			gdImageLine(	imagePtr, xBotLeft -5,
								(Int2)(yBotLeft ),
								xBotLeft +5,(Int2)(yBotLeft ),
								arrayOfColor[1]);
			}
			/*print the highest value*/
/*			sprintf(pcSmall,"%d", (int)flHigh ) ;

			gdImageString(	imagePtr, gdFont6X12, 15 ,
								(Int2)(yTopLeft+ 12),
								pcSmall, arrayOfColor[1]);
			gdImageLine(	imagePtr, xBotLeft -5,
								(Int2)(yTopLeft +15),
								xBotLeft +5,(Int2)(yTopLeft + 15),
								arrayOfColor[1]);*/

			for (cnt= 1 ;(Int2)(yBotLeft -( (flvmargin + ((FloatLo)istep * (FloatLo)(cnt -1)))* flvstep)) > yTopLeft + 5; cnt++)
			{
				sprintf(pcSmall,"%d",(iMinOnVScale + (istep  * cnt) ) ) ;
				if((tracingMode == CONTACT_MAP ) || (tracingMode == DOT_ASA_VS_TIME) || (tracingMode == ENERGY_PLOTTING)
					|| (tracingMode ==DOT_2NDRY_STUC_VS_TIME) || (tracingMode == DOT_AVG_2NDRY_STUC_VS_TIME) ) {
					gap=(flHigh-flLow+8)/10;
					cur=gap*(cnt);
					if (cnt==10)
						cur=flHigh-flLow;
					sprintf(pcSmall,"%d",cur) ;
					gdImageString(	imagePtr, gdFont6X12, 15 ,
									(Int2)(yBotLeft -( cur* flvstep) -5),
									pcSmall, arrayOfColor[1]);
					gdImageLine(	imagePtr, xBotLeft -5,
									(Int2)(yBotLeft -( cur* flvstep)),
									xBotLeft +5,(Int2)(yBotLeft -cur* flvstep),
									arrayOfColor[1]);
					if (cnt==10)
						break;
				}
				else  {
					if((Int2)(yBotLeft -( (flvmargin + ((FloatLo)istep * (FloatLo)(cnt -1)))* flvstep) -5) < yBotLeft -15) {
						gdImageString(	imagePtr, gdFont6X12, 15 ,
										(Int2)(yBotLeft -( (flvmargin + ((FloatLo)istep * (FloatLo)(cnt -1)))* flvstep) -5),
										pcSmall, arrayOfColor[1]);
						gdImageLine(	imagePtr, xBotLeft -5,
										(Int2)(yBotLeft -( (flvmargin + ((FloatLo)istep * (FloatLo)(cnt -1)))* flvstep)),
										xBotLeft +5,(Int2)(yBotLeft -( (flvmargin + ((FloatLo)istep * (FloatLo)(cnt -1)))* flvstep)),
										arrayOfColor[1]);
					}
				}
			}
		}

		/*If the range is small, then we print float intermediate values*/
		else{
			/*print extreme values*/
			/*print the lowest value*/
			sprintf(pcSmall,"%5.4f", flLow ) ;

			gdImageString(	imagePtr, gdFont6X12, 5 ,
								(Int2)(yBotLeft -5 ),
								pcSmall, arrayOfColor[1]);
			gdImageLine(	imagePtr, xBotLeft -5,
								(Int2)(yBotLeft ),
								xBotLeft +5,(Int2)(yBotLeft ),
								arrayOfColor[1]);

			/*print the highest value*/
			/*sprintf(pcSmall,"%5.4f", flHigh ) ;

			gdImageString(	imagePtr, gdFont6X12, 5 ,
								(Int2)(yTopLeft+ 12),
								pcSmall, arrayOfColor[1]);
			gdImageLine(	imagePtr, xBotLeft -5,
								(Int2)(yTopLeft +15),
								xBotLeft +5,(Int2)(yTopLeft + 15),
								arrayOfColor[1]);*/

			flStep  = (flHigh - flLow) /10.0 ;
			flvstep = (FloatLo)yWide /(flHigh - flLow);

			for (cnt=1;(yBotLeft -(( (flStep * (FloatLo)(cnt -1)))* flvstep)) > yTopLeft + 15 ; cnt++){
				sprintf(pcSmall,"%5.4f",(FloatLo) flLow+(flStep  * (FloatLo)(cnt-1)) ) ;
				if((Int2)(yBotLeft - ( ((flStep * (FloatLo)(cnt -1)))* flvstep)  -5) < yBotLeft -15) {
					gdImageString(	imagePtr, gdFont6X12, 5 ,
									(Int2)(yBotLeft -( ( ((FloatLo)flStep * (FloatLo)(cnt -1)))* flvstep) -5),
									pcSmall, arrayOfColor[1]);
					gdImageLine(	imagePtr, xBotLeft -5,
									(Int2)(yBotLeft -( ( ((FloatLo)flStep * (FloatLo)(cnt -1)))* flvstep)),
									xBotLeft +5,(Int2)(yBotLeft -( ( ((FloatLo)flStep * (FloatLo)(cnt -1)))* flvstep)),
									arrayOfColor[1]);
				}
			}
		}


		/*horizontally adds the model numbers*/
		if( (tracingMode == ARRAY_VS_TIME) || (tracingMode == DOT_ASA_VS_TIME) || (tracingMode == ENERGY_PLOTTING)
			|| (tracingMode ==DOT_2NDRY_STUC_VS_TIME) || (tracingMode == DOT_AVG_2NDRY_STUC_VS_TIME) )
		{
			flxstep = xWide / NUM_MODEL;
			/*draw the legend horizontaly on the diagram and insert the indices*/
			for (cnt= 0 ;cnt< ((numModels+9)/10+1) ; cnt++){
				sprintf(pcSmall,"%d",(Int2)(cnt * 10) );
				gdImageString(	imagePtr, gdFont6X12, ( (xTopLeft-10) + cnt * (xWide / (NUM_MODEL/10)) ),
								yBotLeft + 10 , pcSmall, arrayOfColor[1]);
				gdImageLine(	imagePtr,  (xBotLeft) + (cnt * xWide / (NUM_MODEL/10) ) ,
								yBotLeft -5,  (xBotLeft) + (cnt * xWide / (NUM_MODEL/10)) ,
								yBotLeft +5, arrayOfColor[1]);
			}
		}

		/*horizontally adds the residue numbers*/
		else if ( tracingMode ==ARRAY_VS_RESIDUE){
			flxstep = xWide / numModels;
			/*draw the legend horizontaly on the diagram and insert the indices*/
			for (cnt= 1 ;cnt<= 10; cnt++){
				gap=(numModels+8)/10;
				cur=gap*cnt;
				if (cnt==10)
					cur=numModels;
				sprintf(pcSmall,"%d",cur );
				gdImageString(	imagePtr, gdFont6X12, ( (xTopLeft-10) + cur * flxstep ),
								yBotLeft + 10 , pcSmall, arrayOfColor[1]);
				gdImageLine(	imagePtr,  xBotLeft  + cur * flxstep,
								yBotLeft -5,  xBotLeft  + cur * flxstep,
								yBotLeft +5, arrayOfColor[1]);
			}
		}
		else if (tracingMode == CONTACT_MAP ){
			flxstep = xWide / flHigh1;
			/*draw the legend horizontaly on the diagram and insert the indices*/
			for (cnt= 1 ;cnt<= 10; cnt++){
				gap=(flHigh1+8)/10;
				cur=gap*cnt;
				if (cnt==10)
					cur=flHigh1;
				sprintf(pcSmall,"%d",cur);
				gdImageString(	imagePtr, gdFont6X12, ( (xTopLeft-10) + cur * flxstep ),
								yBotLeft + 10 , pcSmall, arrayOfColor[1]);
				gdImageLine(	imagePtr,  (xBotLeft ) + cur * flxstep ,
								yBotLeft -5,  (xBotLeft ) + cur * flxstep ,
								yBotLeft +5, arrayOfColor[1]);
			}
		}
		else if(tracingMode == RGYR_VS_RMSD){
			iMaxRmsd = (findGreatestFloat(flArray2, numModels));
			flxstep  = (FloatLo)xWide / (FloatLo)iMaxRmsd ;
			/*draw the legend horizontaly on the diagram and insert the indices*/
			for (cnt= 1 ;cnt<= 10; cnt++){
				sprintf(pcSmall,"%3.1f",(float)((FloatLo)cnt * (FloatLo)((FloatLo)iMaxRmsd/10.0)) );
				gdImageString(	imagePtr, gdFont6X12, ( (xTopLeft-10) + cnt * (xWide / 10) ),
								yBotLeft + 10 , pcSmall, arrayOfColor[1]);
				gdImageLine(	imagePtr,  (xBotLeft ) + (cnt * xWide / 10 ) ,
								yBotLeft -5,  (xBotLeft ) + (cnt * xWide / 10 ) ,
								yBotLeft +5, arrayOfColor[1]);
			}
		}

		/*Horizontaly adds the contact number*/
		else if(tracingMode == RECTANGLE_VS_CONTACT){
			flxstep  = (FloatLo)xWide / (FloatLo)numModels ;
			/*draw the legend horizontaly on the diagram and insert the indices*/
			for (cnt= 1 ;cnt<= 10; cnt++){
				gap=(numModels+8)/10;
				cur=gap*cnt;
				if (cnt==10)
					cur=numModels;
				sprintf(pcSmall,"%d",cur);

				gdImageString(	imagePtr, gdFont6X12, ( (xTopLeft-5) + cur * flxstep ),
								yBotLeft + 10 , pcSmall, arrayOfColor[1]);

				gdImageLine(	imagePtr,  ((FloatLo)xBotLeft ) + ((FloatLo)(cur+1) * flxstep ) ,
								yBotLeft -5,  ((FloatLo)xBotLeft ) + ((FloatLo)(cur+1) * flxstep ) ,
								yBotLeft +5, arrayOfColor[1]);
			}
		}
	}

	/*add the additional plots to the graph if this is not a new graph*/
	else{
		printf("\nReopening the gif: %s\n", gifFileName);
		fileIn = FileOpen(gifFileName, "rb");
		if(fileIn == NULL){
			ErrPostEx(SEV_FATAL,0,0,"Unable to open file %s",gifFileName);
			exit(1);
		}
		imagePtr = 	gdImageCreateFromGif(fileIn);
		/*define the colors*/
		defineColor(imagePtr, arrayOfColor);
		FileClose(fileIn);

		if (tracingMode == ARRAY_VS_TIME){
			if( flArray1 != NULL){
				sprintf(pcTemp,"%s%s","_",pcLegend1 );
				gdImageString(imagePtr, gdFont6X12,(60 * graphNum) , 5, pcTemp, arrayOfColor[(graphNum + 1) % NUM_OF_COLORS]);
			}
			if( flArray2 != NULL){
				sprintf(pcTemp,"%s%s","..",pcLegend2);
				gdImageString(imagePtr, gdFont6X12,(60 * graphNum) , 15, pcTemp, arrayOfColor[(graphNum + 1) % NUM_OF_COLORS]);
			}
		}
	}


	/*NOW PLOTTING IN THE DIFFERENT MODES*/
	if (tracingMode == ARRAY_VS_TIME){
		if( flArray1 != NULL)
			yPrevious1 =   (yBotLeft -  ( flArray1[0] - flLow) * ( yUnit ));
		if( flArray2 != NULL)
			yPrevious2 =   (yBotLeft -  ( flArray2[0] - flLow) * ( yUnit ));
		xPrevious = xBotLeft + xUnit;

		for (cnt=1; cnt<= numModels ;cnt++)
		{
			if( flArray1 != NULL){
				gdImageLine	(imagePtr, xPrevious, yPrevious1 ,
							(xBotLeft + (cnt) * xUnit ), ( yBotLeft -( (flArray1[cnt-1] -flLow) * ( yUnit )) ) ,
							 arrayOfColor[(graphNum + 1) % NUM_OF_COLORS ]);
				yPrevious1   = ( yBotLeft - ( (flArray1[cnt-1]   - flLow) * ( yUnit )));
			}

			if (flArray2 != NULL){
				gdImageDashedLine(imagePtr, xPrevious, yPrevious2,
								 (xBotLeft + (cnt) * xUnit ),( yBotLeft -( ( flArray2[cnt-1] - flLow) * ( yUnit )) ) ,
								  arrayOfColor[(graphNum +1) % NUM_OF_COLORS]);
				yPrevious2 = ( yBotLeft - ( (flArray2[cnt-1] - flLow) * ( yUnit )));
			}
			xPrevious = (xBotLeft + (cnt) * xUnit );
		}

	}else if (tracingMode == RECTANGLE_VS_CONTACT){
		if( flArray1 != NULL)
			yPrevious1 =   (yBotLeft -  (( flArray1[0] - flLow) * ( yUnit )) );
		xPrevious = xBotLeft + (((FloatLo)xWide) / (FloatLo)numModels);

		for (cnt=1; cnt<= numModels ;cnt++)
		{
			/*gdImageLine	(imagePtr, xPrevious, yPrevious1 ,
						(xBotLeft + (cnt +1) * xUnit ), ( yBotLeft -( (flArray1[cnt] -flLow) * ( yUnit )) ) ,
						 arrayOfColor[(graphNum + 1) % NUM_OF_COLORS ]);*/
			gdImageFilledRectangle(imagePtr, 
					((FloatLo)xBotLeft + ((FloatLo)cnt     * (((FloatLo)xWide) / (FloatLo)numModels))), yPrevious1,
					((FloatLo)xBotLeft + ((FloatLo)(cnt+ 1.0) * (((FloatLo)xWide) / (FloatLo)numModels)))-1, yBotLeft-1, 
					arrayOfColor[cnt%NUM_OF_COLORS]);				

			yPrevious1   = ( yBotLeft - ( (flArray1[cnt-1]   - flLow) * ( yUnit )));

			xPrevious = (xBotLeft + (FloatLo)cnt + 1.0) * (((FloatLo)xWide) / (FloatLo)numModels) ;
		}

	}else if (tracingMode == ARRAY_VS_RESIDUE ){
		xUnit=xWide / numModels;
		if( flArray1 != NULL) 
			yPrevious1 =   (yBotLeft -  ( flArray1[0] - flLow) * ( yUnit ));
		if( flArray2 != NULL) 
			yPrevious2 =   (yBotLeft -  ( flArray2[0] - flLow) * ( yUnit ));
		xPrevious = xBotLeft + xUnit;

		for (cnt=1; cnt<= numModels ;cnt++)
		{
			if( flArray1 != NULL){ 
				gdImageLine	(imagePtr, xPrevious, yPrevious1 ,
							(xBotLeft + (cnt) * xUnit ), ( yBotLeft -( (flArray1[cnt-1] -flLow) * ( yUnit )) ) ,
							 arrayOfColor[(graphNum + 1) % NUM_OF_COLORS ]);
				yPrevious1   = ( yBotLeft - ( (flArray1[cnt-1]   - flLow) * ( yUnit )));
			}

			if (flArray2 != NULL){
				gdImageDashedLine(imagePtr, xPrevious, yPrevious2, 
								 (xBotLeft + (cnt) * xUnit ),( yBotLeft -( ( flArray2[cnt-1] - flLow) * ( yUnit )) ) ,
								  arrayOfColor[(graphNum +1) % NUM_OF_COLORS]);
				yPrevious2 = ( yBotLeft - ( (flArray2[cnt-1] - flLow) * ( yUnit )));
			}
			xPrevious = (xBotLeft + (cnt) * xUnit );
		}

	}else if(tracingMode == DOT_ASA_VS_TIME){
		minASA = findSmallestIntInppInt((Int4**)option1,numModels, (Int4)flHigh1 );
		maxASA = findGreatestIntInppInt((Int4**)option1,numModels, (Int4)flHigh1 );
		printf ("\nRange of ASA: \n-min:%d \n-max:%d", minASA, maxASA);
		defineRainbowColor(imagePtr, arrayRainbowColors);
		
		/*graph the legend*/
		gdImageString(imagePtr, gdFont7X13b, xSize - 120 , 15 , "Legend: in A^2", arrayOfColor[1]);
		for(cnt=1; cnt<=NUM_COLOR_RAINBOW; cnt++){
				i = minASA + ((cnt-1) * ((maxASA - minASA)/NUM_COLOR_RAINBOW));
				j = minASA + ( cnt * ((maxASA - minASA)/NUM_COLOR_RAINBOW));
				sprintf(pcTemp,"%d to %d ",i,j);
				gdImageFilledRectangle(imagePtr,
						xSize - 107, 30+ (cnt * 15),
						xSize - 102, 40+ (cnt * 15),
						arrayRainbowColors[cnt-1]);

				gdImageString(imagePtr, gdFont6X12,
							  xSize - 100 , 30 + (cnt * 15),
							  pcTemp,arrayOfColor[1]);
		}

		/*insert the colored squared*/
		for (i=0; i< numModels; i++){
			for(j=0; j<(Int2)flHigh; j++){
				sprintf(pcTemp, "%d", (Int4)( ((Int4**)option1)[i][j] ));
				cnt =(NUM_COLOR_RAINBOW) *(Int4)(((Int4**)option1)[i][j] - minASA+1)/(maxASA - minASA);
				gdImageFilledRectangle(imagePtr,
						(xBotLeft + (i * flxstep))+1, (yBotLeft - (flvstep * (j+1))),
						(xBotLeft + ((i+1) * flxstep)), (yBotLeft - (flvstep * j  )) -1,
						arrayRainbowColors[cnt]);
				/*gdImageString(imagePtr, gdFont6X12, (xBotLeft + (i * flxstep)),
								(yBotLeft - (flvstep * (j+1))) , pcTemp, arrayOfColor[1]);*/
			}
		}

	}else if(tracingMode == ENERGY_PLOTTING){

		flMinEn = findSmallestFloatInppFloat((FloatLo**)option1,numModels, flHigh1 );
		flMaxEn = findGreatestFloatInppFloat((FloatLo**)option1,numModels, flHigh1 );
		printf ("\nRange of Energy: \n-min:%5.3f \n-max:%5.3f", flMinEn, flMaxEn);
		defineRainbowColor(imagePtr, arrayRainbowColors);

		/*graph the legend*/
		gdImageString(imagePtr, gdFont7X13b, xSize - 120 , 15 , "Energy", arrayOfColor[1]);
		for(cnt=1; cnt<=NUM_COLOR_RAINBOW; cnt++){
				i = flMinEn + ((cnt-1) * ((flMaxEn - flMinEn)/NUM_COLOR_RAINBOW));
				j = flMinEn + ( cnt * ((flMaxEn - flMinEn)/NUM_COLOR_RAINBOW));
				sprintf(pcTemp,"%d to %d ",i,j);
				gdImageFilledRectangle(imagePtr, 
						xSize - 107, 30+ (cnt * 15),
						xSize - 102, 40+ (cnt * 15), 
						arrayRainbowColors[cnt-1]);					

				gdImageString(imagePtr, gdFont6X12, 
							  xSize - 100 , 30 + (cnt * 15), 
							  pcTemp,arrayOfColor[1]);
		}

		/*insert the colored squared*/
		for (i=0; i< numModels; i++){
			for(j=0; j<(Int2)flHigh; j++){
				sprintf(pcTemp, "%f", (FloatLo)( ((FloatLo**)option1)[i][j] ));
				cnt =(NUM_COLOR_RAINBOW) *(FloatLo)(((FloatLo**)option1)[i][j] - flMinEn+1)/(flMaxEn - flMinEn);
				gdImageFilledRectangle(imagePtr, 
						(xBotLeft + (i * flxstep))+1, (yBotLeft - (flvstep * (j+1))),
						(xBotLeft + ((i+1) * flxstep)), (yBotLeft - (flvstep * j  )) -1, 
						arrayRainbowColors[cnt]);				
				/*gdImageString(imagePtr, gdFont6X12, (xBotLeft + (i * flxstep)), 
								(yBotLeft - (flvstep * (j+1))) , pcTemp, arrayOfColor[1]);*/
			}
		}

	}else if( (tracingMode == DOT_2NDRY_STUC_VS_TIME) || (tracingMode == DOT_AVG_2NDRY_STUC_VS_TIME) ){
		defineRainbowColor(imagePtr, arrayRainbowColors);
		
		/*graph the legend*/
		gdImageString(imagePtr, gdFont7X13b, xSize - 120 , 15 , "Legend:", arrayOfColor[1]);
		j=30;
		if(tracingMode == DOT_2NDRY_STUC_VS_TIME ){
			for(cnt=1; cnt<= (NUM_KIND_OF_2NDRY_STRUCT +1); cnt++){
				if(cnt ==1){StringCpy(pcTemp,"H: alpha-helix"); i= 1;}
				if(cnt ==2){StringCpy(pcTemp,"I: pi-helix");    i= 0;}
				if(cnt ==3){StringCpy(pcTemp,"G: 3_10 helix");  i= 2;}
				if(cnt ==4){StringCpy(pcTemp, "B: Isolated");
							StringCpy(pcSmall,"Beta bridge");
							i= 9;}
				if(cnt ==5){StringCpy(pcTemp,  "E: Ext. strand");
							StringCpy(pcSmall, "found in Beta-ladder");
							i=8 ;}
				if(cnt ==6){StringCpy(pcTemp,"T: H-bonded turn"); i= 7;}
				if(cnt ==7){StringCpy(pcTemp,"S: bend"); i= 6;}
				if(cnt ==8){StringCpy(pcTemp,"Not in 2ndry struc."); i= 4;}

				gdImageFilledRectangle(imagePtr, 
						xSize - 127, j+ (cnt * 15),
						xSize - 122, j +10+ (cnt * 15), 
						arrayRainbowColors[i]);

				gdImageString(imagePtr, gdFont6X12, 
							  xSize - 120 , j + (cnt * 15), 
							  pcTemp,arrayOfColor[1]);
				if( (cnt==4) || (cnt==5) ){
					j= j+10;
					gdImageString(	imagePtr, gdFont6X12, 
									xSize - 120 , j + (cnt * 15), 
									pcSmall,arrayOfColor[1]);
				}
			}
		}
		else if(tracingMode == DOT_AVG_2NDRY_STUC_VS_TIME){
			for(cnt=1; cnt<= (NUM_GROUP_2NDRY_STRUCT); cnt++){
				if(cnt ==1){StringCpy(pcTemp," Helix"); i= 1;}
				if(cnt ==2){StringCpy(pcTemp," Sheet"); i= 9;}
				if(cnt ==3){StringCpy(pcTemp," Coils and Others"); i= 7;}

				gdImageFilledRectangle(imagePtr, 
						xSize - 127, j+ (cnt * 15),
						xSize - 122, j +10+ (cnt * 15), 
						arrayRainbowColors[i]);					

				gdImageString(imagePtr, gdFont6X12, 
							  xSize - 120 , j + (cnt * 15), 
							  pcTemp,arrayOfColor[1]);
			}
		}
		/*insert the colored squared*/
		for (i=0; i< numModels; i++){
			for(j=0; j<(Int2)flHigh; j++){
				sprintf(pcTemp, "%c", (char)( ((CharPtr*)option1)[i][j] ));
				if (tracingMode == DOT_2NDRY_STUC_VS_TIME){
					if((char)( ((CharPtr*)option1)[i][j]) == 'H'){cnt=1;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'I'){cnt=0;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'G'){cnt=2;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'B'){cnt=9;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'E'){cnt=8;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'T'){cnt=7;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'D'){cnt=6;}
					else if((char)( ((CharPtr*)option1)[i][j]) == ' '){cnt=4;}
				}
				else if(tracingMode == DOT_AVG_2NDRY_STUC_VS_TIME){
					if     ((char)( ((CharPtr*)option1)[i][j]) == 'H'){cnt=1;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'E'){cnt=9;}
					else if((char)( ((CharPtr*)option1)[i][j]) == 'C'){cnt=7;}
				}
				gdImageFilledRectangle(imagePtr, 
						(xBotLeft + (i * flxstep))+1, (yBotLeft - (flvstep * (j+1))),
						(xBotLeft + ((i+1) * flxstep)), (yBotLeft - (flvstep * j  )) -1, 
						arrayRainbowColors[cnt]);

				/*gdImageString(imagePtr, gdFont6X12, (xBotLeft + (i * flxstep)),
				(yBotLeft - (flvstep * (j+1))) , pcTemp, arrayOfColor[1]);*/
			}
		}
	}
	else if (tracingMode == CONTACT_MAP){
		defineRainbowColor(imagePtr, arrayRainbowColors);
		/*graphs the legend*/
		for(cnt=5; cnt<=NUM_COLOR_RAINBOW; cnt++){
			i =  (int)((((FloatLo)cnt-5.0)/ ((FloatLo)NUM_COLOR_RAINBOW - 4.0) * 100.0 )+0.5);
			j =  (int)((((FloatLo)cnt-4.0)/ ((FloatLo)NUM_COLOR_RAINBOW - 4.0) * 100.0 )+0.5);
			sprintf(pcTemp,"%d%% to %d%% ",i,j);
			gdImageFilledRectangle(imagePtr, 
					xSize - 107, 30+ (cnt * 15),
					xSize - 102, 40+ (cnt * 15), 
					arrayRainbowColors[cnt-1]);					

			gdImageString(imagePtr, gdFont6X12, 
						  xSize - 95 , 30 + (cnt * 15), 
							  pcTemp,arrayOfColor[1]);
		}

		/*iterate throught the contact arrays*/
		for(i=0; i<(Int2)flHigh; i++){
			for(j=0; j<=(Int2)i; j++){
				/*choose the color accordingly to the density*/
				cnt = ((FloatLo)NUM_COLOR_RAINBOW - 5.0 ) *  (FloatLo)(((FloatLo***)option1)[numModels][i][j]);
				cnt = cnt + 4;
				/*graph the native contacts*/
				/*case where the contact is present in the native structure*/
				/*the plot is then traced in the upper left side*/
				if((FloatLo)(((FloatLo***)option1)[numModels][i][j]) != 0.0 && (FloatLo)(((FloatLo***)option1)[0][i][j]) != 0.0){
						gdImageFilledRectangle(imagePtr, 
								(xBotLeft + (j * flxstep))+1, (yBotLeft - (flvstep * (i+1))),
								(xBotLeft + ((j+1) * flxstep)), (yBotLeft - (flvstep * i  )) -1, 
								arrayRainbowColors[cnt]);
				}
				/*if it's a new contact, graph is traced on the bottom right side*/
				if( ((FloatLo)(((FloatLo***)option1))[numModels][i][j] != 0.0) && ((FloatLo)(( (FloatLo***)option1)[0][i][j]) == 0.0) ){
						gdImageFilledRectangle(imagePtr, 
								(xBotLeft + (i     * flxstep))+1, (yBotLeft - (flvstep * (j+1))),
								(xBotLeft + ((i+1) * flxstep)), (yBotLeft - (flvstep *  j)) -1, 
								arrayRainbowColors[cnt]);
					}

				/*else{
						gdImageFilledRectangle(imagePtr, 
								(xBotLeft + (i     * flxstep))+1, (yBotLeft - (flvstep * (j+1))),
								(xBotLeft + ((i+1) * flxstep))+1, (yBotLeft - (flvstep *  j)) -1, 
								arrayRainbowColors[cnt]);
				}*/
			}
		}
	}else if(tracingMode == RGYR_VS_RMSD){

		flMinEn = findSmallestFloat((FloatLo*)option1,numModels );
		flMaxEn = findGreatestFloat((FloatLo*)option1,numModels );
		printf ("\nRange of Energy: \n-min:%5.3f \n-max:%5.3f", flMinEn, flMaxEn);
		defineRainbowColor(imagePtr, arrayRainbowColors);
		
		/*graph the legend*/
		gdImageString(imagePtr, gdFont7X13b, xSize - 120 , 15 , "Bryant Energy", arrayOfColor[1]);
		for(cnt=1; cnt<=NUM_COLOR_RAINBOW; cnt++){
				i = flMinEn + ((cnt-1) * ((flMaxEn - flMinEn)/NUM_COLOR_RAINBOW));
				j = flMinEn + ( cnt * ((flMaxEn - flMinEn)/NUM_COLOR_RAINBOW));
				sprintf(pcTemp,"%d to %d ",i,j);
				gdImageFilledRectangle(imagePtr, 
						xSize - 107, 30+ (cnt * 15),
						xSize - 102, 40+ (cnt * 15), 
						arrayRainbowColors[cnt-1]);					

				gdImageString(imagePtr, gdFont6X12, 
							  xSize - 100 , 30 + (cnt * 15), 
							  pcTemp,arrayOfColor[1]);
		}

		/*plot each dot the rmsd(horyzontal axis) and rgyr(vertical) define the position of the dot while
		  the color of the dot is defined by the Bryant energy level*/
		yMin = findSmallestFloat(flArray1, numModels);
		for(i=0; i< numModels; i++){
				cnt =(NUM_COLOR_RAINBOW) *( ( ((FloatLo*)option1)[i] - (flMinEn+1))/(flMaxEn - flMinEn));

				xTemp =  flArray2[i];
				yTemp =  flArray1[i] - yMin;

				/*we don't want the square to go over the vertical limit*/
				if( ((yBotLeft - (flvstep * (yTemp))) -RGYR_VS_RMSD_SQUARE_SIZE) <  yTopLeft + 15){
					/*we don't want it to go over the right limit either*/
					if ( (xBotLeft + (xTemp * flxstep)) > (xBotRight - (15 + RGYR_VS_RMSD_SQUARE_SIZE))){
					gdImageFilledRectangle(imagePtr, 
							(xBotLeft + (xTemp * flxstep))- (RGYR_VS_RMSD_SQUARE_SIZE + 15), (yBotLeft - (flvstep * (yTemp /*+ flvstep*/))),
							(xBotLeft + (xTemp * flxstep)) -15, (yBotLeft - (flvstep * yTemp ))+RGYR_VS_RMSD_SQUARE_SIZE, 
							arrayRainbowColors[cnt]);	
					
					}else{
					gdImageFilledRectangle(imagePtr, 
							(xBotLeft + (xTemp * flxstep)), (yBotLeft - (flvstep * (yTemp /*+ flvstep*/))),
							(xBotLeft + (xTemp  * flxstep))+RGYR_VS_RMSD_SQUARE_SIZE, (yBotLeft - (flvstep * yTemp ))+RGYR_VS_RMSD_SQUARE_SIZE, 
							arrayRainbowColors[cnt]);	
					}
				}else{
					/*we don't want it to go over the right limit either*/
					if ( (xBotLeft + (xTemp * flxstep)) > (xBotRight - (15+ RGYR_VS_RMSD_SQUARE_SIZE))){
					gdImageFilledRectangle(imagePtr, 
							(xBotLeft + (xTemp * flxstep))- (RGYR_VS_RMSD_SQUARE_SIZE + 15), (yBotLeft - (flvstep * yTemp )) -RGYR_VS_RMSD_SQUARE_SIZE,
							(xBotLeft + (xTemp * flxstep)) - 15, (yBotLeft - (flvstep * yTemp )), 
							arrayRainbowColors[cnt]);
					
					}else{
					gdImageFilledRectangle(imagePtr, 
							(xBotLeft + (xTemp * flxstep))+1.0, (yBotLeft - (flvstep * yTemp )) -RGYR_VS_RMSD_SQUARE_SIZE,
							(xBotLeft + (xTemp  * flxstep))+RGYR_VS_RMSD_SQUARE_SIZE+1.0, (yBotLeft - (flvstep * yTemp )),
							arrayRainbowColors[cnt]);	
					}
				}
			}
		}

	fileOut = FileOpen(gifFileName, "wb");
	advGdImageGif(imagePtr, fileOut);
	FileClose(fileOut);
	gdImageDestroy(imagePtr);
	return;
}

Boolean CreateWorld(PWS ** pppwsThis,  PMSD * pmsdArray, Int2 iNumSimulation, Int4 iNumRes, Int4 iNumModels,PMMD pmmdNative){
	PMAD	pmadThis;
	PMMD	pmmdThis;
	PDNMG	pdnmgThis;
	PMGD	pmgdThis;
	PVNMA	pvnmaThis;

	Int2 i,j,k;

	/*	Criteria of atoms included in the world:
		Heavy atoms between non-neighbooring residues within 5.4 A for alephatic C
		& 4.6A for all other atoms*/
	if(pmsdArray[0] == NULL){return FALSE;}

	for(i=0; i< iNumSimulation; i++){
		pmmdThis	= ((PMMD)((PDNMM)pmsdArray[i]->pdnmmHead)->data.ptrvalue);

		/*create the world for each model*/
		for(j=0; j< iNumModels; j++){
			pdnmgThis	= pmmdThis->pdnmgHead;
            if (j==0 && pmmdNative!=NULL){
                pdnmgThis=pmmdNative->pdnmgHead;
            }

			/*go through each residue*/
			for(k=0; k< iNumRes; k++){
  				pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
 				pvnmaThis = pmgdThis->pvnmaAHead;

				/*add the heavy atoms to the world*/
				while(pvnmaThis){
					pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
					if( (pmadThis ->pvnmaLink->choice ) != 1){
						if( pppwsThis[i][j] == NULL){
							pppwsThis[i][j] = AddtoWorld(pppwsThis[i][j], j, (PFB)pmadThis);
						}
						else{
							AddtoWorld(pppwsThis[i][j], j, (PFB)pmadThis);
						}

					}
					pvnmaThis = pvnmaThis->next;
				}
				pdnmgThis = pdnmgThis->next;
			}
			InstantiateWorld(1, pppwsThis[i][j]);
		}
	}
	return TRUE;
}

Boolean IsAromatic(PMAD pmadThis){
	ValNodePtr pvnHere;

	pvnHere = pmadThis->pvnBonds;

	while (pvnHere){
		if( IsBondAromatic((PFB)pvnHere->data.ptrvalue)){
			return TRUE;
		}	
		pvnHere=pvnHere->next;
	}
	return FALSE;
}


Boolean RetrieveContacts(PMSD * pmsdArray, PWS ** pppwsThis, Boolean **** bContactMap, Int2 iNumSimulation, Int2 iNumRes, Int2 iNumModels,PMMD pmmdNative){
	ValNodePtr vnpList,vnpListHere;
	PMAD pmadThis;
	PALD paldLoc,paldRowAtom;
	PMGD pmgdThis;
	PMMD pmmdThis;
	PVNMA pvnmaThis;
	PDNMG pdnmgThis;
	Boolean bIsAliphC1/*, bIsAliphC2*/;
	vec vdist/*, vdist2, res*/;
	/*FloatLo distance;*/
	Int2 resNumAtom,resNumContact,/*molnum,*/i,j;

	if (pppwsThis==NULL) return FALSE;

	for(i=0; i< iNumSimulation; i++){

		pmmdThis= ((PMMD)((PDNMM)pmsdArray[i]->pdnmmHead)->data.ptrvalue);

		/*create the world for each model*/
        for(j=0; j< iNumModels; j++){
            if (j==0){
                if(pmmdNative!=NULL){
                    pmmdThis=pmmdNative;
                }
                else {
                    continue;
                }
            }
            
            pdnmgThis	= pmmdThis->pdnmgHead;

			/*go through each residue*/
			while(pdnmgThis){
  				pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
 				pvnmaThis = pmgdThis->pvnmaAHead;

				/*add the heavy atoms to the world*/
				while(pvnmaThis){
					pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
					if( (pmadThis ->pvnmaLink->choice ) != 1){

           				paldLoc=GetAtomLocs(pmadThis,pppwsThis[i][j]->Model);
						if (paldLoc!=NULL) {

							vdist[0]=AtomLocX(paldLoc);
							vdist[1]=AtomLocY(paldLoc);
        					vdist[2]=AtomLocZ(paldLoc);

							/* saves residue number */
							pmgdThis=GetParentGraph((PFB)paldLoc);
							resNumAtom=(pmgdThis->pdnmgLink)->choice;

							/*check if the current atom is an aromatic atom*/
							bIsAliphC1	= (!IsAromatic(pmadThis)) && (pmadThis->pvnmaLink->choice  == 12);

							/* get list of PALDs "near" this atom */
							if(!bIsAliphC1){
								vnpList=FindAtomsIn(pppwsThis[i][j],vdist, MAX_DIST_CONTACT_OTHER);
							}
							else{
								vnpList=FindAtomsIn(pppwsThis[i][j],vdist, MAX_DIST_CONTACT_ALIPH_C);
							}
							vnpListHere=vnpList;

							/*retrieve the neighboors*/
							while(vnpListHere) {
								paldRowAtom=(PALD)(vnpListHere->data.ptrvalue);

								pmgdThis=GetParentGraph((PFB)paldRowAtom);
								pmmdThis=GetParentMol((PFB)paldRowAtom);

								resNumContact=(pmgdThis->pdnmgLink)->choice;
								/*molnum=(pmmdThis->pdnmmLink)->choice;

								bIsAliphC2 = (!IsAromatic((PMAD)(paldRowAtom->pfbParent))) && (((PMAD)(paldRowAtom->pfbParent))->pvnmaLink->choice  == 12);*/

								/*if both carbons are aliphatic carbons then we have to use a smaller distance MAX_DIST_CONTACT_ALIPH_C*/
								/*if (bIsAliphC1 && !bIsAliphC2){
									if (paldRowAtom!=NULL) {
										vdist2[0]=AtomLocX(paldRowAtom);
										vdist2[1]=AtomLocY(paldRowAtom);
        								vdist2[2]=AtomLocZ(paldRowAtom);
									}
									VecSub(res, vdist,vdist2);
									distance = getMag(res);
									if(distance < MAX_DIST_CONTACT_OTHER){
										bContactMap[i][j][resNumAtom-1][resNumContact-1] = TRUE;

								}
								else{*/
									bContactMap[i][j][resNumAtom-1][resNumContact-1]= TRUE;

								vnpListHere=vnpListHere->next;
							}					
							ValNodeFree(vnpList);
						}
					}
					pvnmaThis = pvnmaThis->next;
				}
				pdnmgThis = pdnmgThis->next;
			}
		}/*end of iteration through all the models*/
	}/*end of iterations through all the simulations*/
	return TRUE;
}

void AverageContacts(Boolean **** bAllContactMap, FloatLo *** pppflFinalContactMap,Int2 iNumSimulation, Int2 iNumRes, Int2 iNumModels, Int2 numPartition){
	Int2 i,j,k,l,temp;

	/*iterate through the contact map and count the number of contacts*/
	/*store the number of contacts in the various partitions*/
	for(i=0; i< iNumSimulation; i++){
		for(j=1; j<iNumModels; j++){  /* start at 1 so we skip over the native contacts in partition 1 */
			for(k=0; k< iNumRes; k++){
				for(l=0; l<= k; l++){
					if(bAllContactMap[i][j][k][l]){
						temp = (Int2)((FloatLo)(j) / ceil((FloatLo) iNumModels/ (FloatLo)numPartition ));
						pppflFinalContactMap[temp+1][k][l] ++; /*temp +1 because the first array contains the native structure*/
					}
				}
			}
		}
	}

	/*average the number of contacts*/
	for(j=0; j<=numPartition; j++){
		for (k=0; k<iNumRes; k++){
			for (l=0; l<=k; l++){
					/*if it's the last partition, we have to consider the remainder of the number
					  of models divided by the number of partitions*/
				if(j == numPartition){
					pppflFinalContactMap[j][k][l] = (FloatLo)pppflFinalContactMap[j][k][l] / ((FloatLo)(iNumModels - ((numPartition -1 ) * (ceil((FloatLo) iNumModels/ (FloatLo)numPartition ))))* (FloatLo)iNumSimulation);
				}else{
					pppflFinalContactMap[j][k][l] = (FloatLo)pppflFinalContactMap[j][k][l] / ((FloatLo)(ceil((FloatLo) iNumModels/ (FloatLo)numPartition )) * iNumSimulation) ;
				}
				if(l==k || k-l ==1 /*|| k-l ==2*/){
					pppflFinalContactMap[j][k][l] = 0.0;
				}
			}
		}
	}

	/*we store the contacts in the native structure*/
	for(k=iNumRes-1; k >= 0; k--){
		for(l=iNumRes-1; l>= k; l--){
			if(bAllContactMap[0][0][k][l]){
				/*pppflFinalContactMap[0][k][l] = 1.0;*/
				pppflFinalContactMap[0][l][k] = 1.0;
			}
			if (l==k || l-k ==1 /*|| l-k == 2*/){
				pppflFinalContactMap[0][l][k] =0.0;
			}
		}
	}
}

void ComputeContactOrder(FloatLo * flpContactOrder, FloatLo * flpContactOrderNoCloseContact ,FloatLo ** flpContactOrderPerRes,
						 FloatLo ** flpCONoCloseContactPerRes,Boolean **** bAllContactMap, Int4 iNumSimulation, Int4 iNumModels,
						 Int4 iNumRes, Int2 iNumPartition){
	FloatLo flSumSeqSep, flSumSeqSepNoClose, * flSumSeqSepPerRes, * flSumSeqSepNoClosePerRes;
	Int4 i4TotalNumContacts, i4TotalNumContactsNoClose, * i4TotalNumContactsPerRes, * i4TotalNumContactsNoClosePerRes;
	Int2 i,j,k,l,m, temp, *iNumModelInPart;

	flSumSeqSepPerRes = (FloatLo*) MemNew( iNumRes * sizeof(FloatLo));
	flSumSeqSepNoClosePerRes= (FloatLo*) MemNew( iNumRes * sizeof(FloatLo));
	i4TotalNumContactsPerRes= (Int4*) MemNew( iNumRes * sizeof(Int4));
	i4TotalNumContactsNoClosePerRes= (Int4*) MemNew( iNumRes * sizeof(Int4));
	iNumModelInPart = (Int2*) MemNew ( iNumPartition * sizeof(Int2));

	for(i=0; i< iNumSimulation; i++){
		/*iterate through each residue and compute the contact order.*/
		for(j=0; j< iNumModels ; j++){
			i4TotalNumContacts = 0;		/*reset to 0 to compute CO of new structure*/
			i4TotalNumContactsNoClose =0;
			flSumSeqSepNoClose = 0.0;
			flSumSeqSep = 0.0;

			for(m=0; m<iNumRes; m++){
				flSumSeqSepPerRes[m]=0.0;
				flSumSeqSepNoClosePerRes[m]=0.0;
				i4TotalNumContactsPerRes[m]=0;
				i4TotalNumContactsNoClosePerRes[m]=0;
			}
			/*iterate through the structure computing the number of contacts as it goes*/
			for(k=0; k< iNumRes; k++){
				for (l=0; l<=k; l++){
					/*If there is a contact, we add it to the number of contact order*/
					if(bAllContactMap[i][j][k][l]){
						i4TotalNumContacts ++;
						flSumSeqSep += ( (FloatLo)k-(FloatLo)l);

						i4TotalNumContactsPerRes[k] ++;
						i4TotalNumContactsPerRes[l] ++;
						flSumSeqSepPerRes[k] += ( (FloatLo)k-(FloatLo)l);
						flSumSeqSepPerRes[l] += ( (FloatLo)k-(FloatLo)l);

						if (k-l > 1){
							i4TotalNumContactsNoClose ++;
							flSumSeqSepNoClose += ( (FloatLo)k-(FloatLo)l);

							i4TotalNumContactsNoClosePerRes[k] ++;
							i4TotalNumContactsNoClosePerRes[l] ++;
							flSumSeqSepNoClosePerRes[k] += ( (FloatLo)k-(FloatLo)l);
							flSumSeqSepNoClosePerRes[l] += ( (FloatLo)k-(FloatLo)l);
						}
					}
				
				}
			}
			/*compute the CO per residue*/
			temp = (Int2)((FloatLo)(j) / ceil((FloatLo) iNumModels/ (FloatLo)iNumPartition ));
			iNumModelInPart[temp] ++;
			for(m=0; m <iNumRes; m++){
				if (i4TotalNumContactsPerRes[m] !=0)
					flpContactOrderPerRes[temp][m]	  +=(( 1.0/ ((FloatLo)iNumRes * (FloatLo)i4TotalNumContactsPerRes[m]) ) * flSumSeqSepPerRes[m]);
				if (i4TotalNumContactsNoClosePerRes[m] !=0)
					flpCONoCloseContactPerRes[temp][m]+=(( 1.0/ ((FloatLo)iNumRes * (FloatLo)i4TotalNumContactsNoClosePerRes[m]) ) * flSumSeqSepNoClosePerRes[m]);
			}

			/*compute the contact order and store it*/
			if (i4TotalNumContacts != 0)
				flpContactOrder[j] += ( ( 1.0/ ((FloatLo)iNumRes * (FloatLo)i4TotalNumContacts) ) * flSumSeqSep);
			if (i4TotalNumContactsNoClose !=0)
				flpContactOrderNoCloseContact[j] += ( ( 1.0/ ((FloatLo)iNumRes * (FloatLo)i4TotalNumContactsNoClose) ) * flSumSeqSepNoClose);
		}
	}
		
	/*go through all the simulations and average the contact orders*/
	for(j=0; j< iNumModels ; j++){
		flpContactOrder[j] = flpContactOrder[j] / (FloatLo)iNumSimulation;
		flpContactOrderNoCloseContact[j] = flpContactOrderNoCloseContact[j] / (FloatLo)iNumSimulation;
	}

	for(m=0; m< iNumPartition;m++){
		for(i=0; i<iNumRes; i++){
			flpContactOrderPerRes[m][i]		= flpContactOrderPerRes[m][i]	/ (iNumModelInPart[m]);
			flpCONoCloseContactPerRes[m][i]	= flpCONoCloseContactPerRes[m][i]/ (iNumModelInPart[m]);
		}
	}

	MemFree (flSumSeqSepPerRes); 
	MemFree (flSumSeqSepNoClosePerRes);
	MemFree (i4TotalNumContactsPerRes);
	MemFree (i4TotalNumContactsNoClosePerRes);
	MemFree (iNumModelInPart);
}	


void ComputeContactMap(PMSD * pmsdArray, Int2 iNumSimulation, Int4 iNumRes, Int4 iNumModels, Int2 numPartition,CharPtr fileName,PMMD pmmdNative){
	Boolean **** bAllContactMap;		/*contains all the contacts computer*/
										/*2 levels before contacts: # of pmsd, # of models*/
	FloatLo *** pppflFinalContactMap;	/*final contact map that will be drawn*/
	PWS ** pppwsThis;
	FloatLo * flpContactOrder, *flpContactOrderNoCloseContact;
	FloatLo ** flpContactOrderPerRes, **flpCONoCloseContactPerRes;
	FloatLo flLowTmp,flHighTmp;
	Int2 i,j,k;
	char pcTemp[PATH_MAX];

	/*Allocate memory for all the arrays*/
	bAllContactMap = (Boolean ****) MemNew(iNumSimulation * sizeof( Boolean ***));
	for(i=0; i< iNumSimulation; i++){
		bAllContactMap[i] = (Boolean ***) MemNew( iNumModels * sizeof (Boolean**));
		for(j=0; j < iNumModels; j++){
			bAllContactMap[i][j] = (Boolean **) MemNew( iNumRes * sizeof (Boolean*));
			for(k=0; k< iNumRes; k++){
				bAllContactMap[i][j][k] = (Boolean *) MemNew( iNumRes * sizeof (Boolean));
			}
		}
	}

	/*NOTE: numPartition +1 patitions are allocated to store the contacts because the first
	 partition will contain the contacts in the native structure*/
	pppflFinalContactMap = (FloatLo ***) MemNew((numPartition+1) * sizeof(FloatLo **));
	for(i=0; i<= numPartition; i++){
		pppflFinalContactMap[i] = (FloatLo **) MemNew(iNumRes * sizeof(FloatLo *));
		for(j=0; j < iNumRes; j++){
			pppflFinalContactMap[i][j] = (FloatLo *) MemNew (iNumRes *sizeof (FloatLo));
		}
	}

	pppwsThis = (PWS **)MemNew(iNumSimulation * sizeof(PWS *) );
	for(i=0; i< iNumSimulation; i++){
		pppwsThis[i] = (PWS*) MemNew( iNumModels * sizeof(PWS) );
	}

	flpContactOrder = (FloatLo *) MemNew (iNumModels * sizeof(FloatLo));
	flpContactOrderNoCloseContact =(FloatLo *) MemNew (iNumModels * sizeof(FloatLo));

	flpContactOrderPerRes= (FloatLo **) MemNew (numPartition * sizeof(FloatLo *));
	for(i=0; i< numPartition; i++){
		flpContactOrderPerRes[i] = (FloatLo *) MemNew (iNumRes *sizeof(FloatLo));
	}

	flpCONoCloseContactPerRes= (FloatLo **) MemNew (numPartition * sizeof(FloatLo *));
	for(i=0; i< numPartition; i++){
		flpCONoCloseContactPerRes[i] = (FloatLo *) MemNew (iNumRes *sizeof(FloatLo));
	}

	/*Create the world(s) */
	CreateWorld(pppwsThis, pmsdArray, iNumSimulation,iNumRes,iNumModels,pmmdNative);

	/*Retrieve the necessary informations to draw the contact map*/
	RetrieveContacts(pmsdArray, pppwsThis, bAllContactMap, iNumSimulation,iNumRes,iNumModels,pmmdNative);

	/*Average the contacts*/
	AverageContacts(bAllContactMap, pppflFinalContactMap ,iNumSimulation,iNumRes,iNumModels, numPartition);

	/*compute the conctact order*/
	ComputeContactOrder(flpContactOrder,flpContactOrderNoCloseContact,flpContactOrderPerRes ,
						flpCONoCloseContactPerRes,bAllContactMap, iNumSimulation,
						iNumModels, iNumRes, numPartition);

	/*printout the contact map to the screen and file*/
	if(iNumSimulation ==1){
		sprintf(pcTemp, "%s_contact_map_part%d%s", fileName, numPartition, ".dat");
		Printpppfl(pppflFinalContactMap, iNumRes,numPartition+1, pcTemp);

		sprintf(pcTemp, "%s_CO%s", fileName, ".dat");
		Printfl(flpContactOrder, iNumModels, "Contact order", pcTemp);
		sprintf(pcTemp, "%s_CO_no_close_contact%s", fileName, ".dat");
		Printfl(flpContactOrderNoCloseContact, iNumModels, "Contact order without close contact", pcTemp);

		for(i=0; i< numPartition; i++){
			sprintf(pcTemp, "%s_CO_per_res_part%d%s", fileName,i ,".dat");
			Printfl(flpContactOrderPerRes[i], iNumRes, "Contact order with close contact per residue", pcTemp);

			sprintf(pcTemp, "%s_CO_no_close_contact_per_res_part%d%s", fileName,i, ".dat");
			Printfl(flpCONoCloseContactPerRes[i], iNumRes, "Contact order without close contact per residue", pcTemp);
		}

	}
	else{
		sprintf(pcTemp, "%s_Avrg_contact_map%s%d%s%d%s",fileName,"_upto",iNumSimulation,"_part_",numPartition,".dat");
		Printpppfl(pppflFinalContactMap, iNumRes,numPartition+1, pcTemp);

		sprintf(pcTemp, "%s_Avrg_CO%s%d%s", fileName,"_upto",iNumSimulation, ".dat");
		Printfl(flpContactOrder, iNumModels, "Contact order", pcTemp);
		sprintf(pcTemp, "%s_Avrg_CO_no_close_contact%s%d%s", fileName,"_upto",iNumSimulation, ".dat");
		Printfl(flpContactOrderNoCloseContact, iNumModels, "Contact order without close contact", pcTemp);

		for(i=0; i< numPartition; i++){
			sprintf(pcTemp, "%s_Avrg_CO_per_res%s%d_part%d%s", fileName,"_upto",iNumSimulation,i ,".dat");
			Printfl(flpContactOrderPerRes[i], iNumRes, "Contact order with close contact per residue", pcTemp);

			sprintf(pcTemp, "%s_Avrg_CO_no_close_contact_per_res%s%d_part%d%s", fileName,"_upto",iNumSimulation,i, ".dat");
			Printfl(flpCONoCloseContactPerRes[i], iNumRes, "Contact order without close contact per residue", pcTemp);
		}


	}

	/*Draw the contact maps and the contact- order graphs*/
	if(iNumSimulation ==1){
		flLowTmp  = 0.0;
		flHighTmp = findGreatestFloat(flpContactOrder, iNumModels);
		sprintf(pcTemp, "%s_CO%s", fileName, ".gif");
		GraphArray(pcTemp, flpContactOrder, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Contact Order", "", "Contact Order",
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		flHighTmp = findGreatestFloat(flpContactOrderNoCloseContact, iNumModels);
		sprintf(pcTemp, "%s_CO_no_close_contact%s", fileName, ".gif");
		GraphArray(pcTemp, flpContactOrderNoCloseContact, NULL, iNumModels, 1, flLowTmp, flHighTmp,
					"Contact Order without considering close contacts", "", "Contact Order",
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		for(i=0; i<=numPartition; i++){

			if (i<numPartition) {
				flHighTmp = findGreatestFloat(flpContactOrderPerRes[i], iNumRes);
				sprintf(pcTemp, "%s_CO_per_res_part%d%s", fileName,i ,".gif");
				GraphArray(pcTemp, flpContactOrderPerRes[i], NULL, iNumRes, 1, flLowTmp, flHighTmp,
						"Contact Order considering close contacts and per residue", "", "Contact Order",
						"Residue #", 645, 500, ARRAY_VS_RESIDUE, NULL);

				flHighTmp = findGreatestFloat(flpCONoCloseContactPerRes[i], iNumRes);
				sprintf(pcTemp, "%s_CO_no_close_contact_per_res_part%d%s", fileName,i, ".gif");
				GraphArray(pcTemp, flpCONoCloseContactPerRes[i], NULL, iNumRes, 1, flLowTmp, flHighTmp,
						"Contact Order without considering close contacts per residue", "", "Contact Order",
						"Residue #", 645, 500, ARRAY_VS_RESIDUE, NULL);
      }
			sprintf(pcTemp, "%s_contact_map_partition%d%s",fileName,i,".gif");
			GraphArray(pcTemp, NULL, NULL, i, 1, 0, iNumRes, "", "", "RESIDUE #", "RESIDUE #",
						645, 500, CONTACT_MAP,pppflFinalContactMap);
		}
	}
	else{
		flLowTmp  = 0.0;
		flHighTmp = findGreatestFloat(flpContactOrder, iNumModels);
		sprintf(pcTemp, "%s_Avrg_CO%s%d%s", fileName,"_upto",iNumSimulation,  ".gif");
		GraphArray(pcTemp, flpContactOrder, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Contact Order", "", "Contact Order",
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		flHighTmp = findGreatestFloat(flpContactOrderNoCloseContact, iNumModels);
		sprintf(pcTemp, "%s_Avrg_CO_no_close_contact%s%d%s", fileName,"_upto",iNumSimulation,  ".gif");
		GraphArray(pcTemp, flpContactOrderNoCloseContact, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Contact Order", "", "Contact Order",
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		for(i=0; i<=numPartition; i++){
			if (i<numPartition) {
				flHighTmp = findGreatestFloat(flpContactOrderPerRes[i], iNumRes);
				sprintf(pcTemp, "%s_Avrg_CO_per_res%s%d_part%d%s", fileName,"_upto",iNumSimulation,i ,".gif");
				GraphArray(pcTemp, flpContactOrderPerRes[i], NULL, iNumRes, 1, flLowTmp, flHighTmp,
						"Contact Order considering close contacts and per residue", "", "Contact Order",
						"Residue #", 645, 500, ARRAY_VS_RESIDUE, NULL);

				flHighTmp = findGreatestFloat(flpCONoCloseContactPerRes[i], iNumRes);
				sprintf(pcTemp, "%s_Avrg_CO_no_close_contact_per_res%s%d_part%d%s", fileName,"_upto",iNumSimulation,i, ".gif");
				GraphArray(pcTemp, flpCONoCloseContactPerRes[i], NULL, iNumRes, 1, flLowTmp, flHighTmp,
						"Contact Order without considering close contacts per residue", "", "Contact Order",
						"Residue #", 645, 500, ARRAY_VS_RESIDUE, NULL);
      }
			sprintf(pcTemp, "%s_Avrg_contact_map%s%d%s%d%s",fileName,"_upto",iNumSimulation,"_parition",i,".gif");
			GraphArray(pcTemp, NULL, NULL, i, 1, 0, iNumRes, "", "", "RESIDUE #", "RESIDUE #", 645, 500, CONTACT_MAP,pppflFinalContactMap);
		}
	}
	
	/*Memory freeing*/
	for(i=0; i< iNumSimulation; i++){
		for(j=0; j < iNumModels; j++){
			for(k=0; k< iNumRes; k++){
				MemFree (bAllContactMap[i][j][k]);
			}
			MemFree (bAllContactMap[i][j] );
		}
		MemFree(bAllContactMap[i]);
	}
	MemFree(bAllContactMap);

	/*for(i=0; i<= numPartition; i++){
		for(j=0; j < iNumRes; j++){
			MemFree(pppflFinalContactMap[i][j]);
		}
		MemFree (pppflFinalContactMap[i]); 
	}
	MemFree (pppflFinalContactMap);*/

	for(i=0; i< iNumSimulation; i++){
		MemFree(pppwsThis[i]);
	}
	MemFree(pppwsThis);

	MemFree(flpContactOrder);
	MemFree(flpContactOrderNoCloseContact);

	for(i= 0; i< numPartition; i++){
		MemFree(flpContactOrderPerRes[i]);
	}
	MemFree(flpContactOrderPerRes);

	for(i= 0; i< numPartition; i++){
		MemFree(flpCONoCloseContactPerRes[i]);
	}
	MemFree(flpCONoCloseContactPerRes);


}

void ComputeEnergyGraphs(PMSD * pmsdArray,Int4 iNumSimulation, Int4 iNumModels, Int4 creaseIncl,
						 Int4 creaseExcl, Boolean creaseDecay,CharPtr fileName, 
						 FloatLo * pflRmsdArray, FloatLo * pflRadiusHP, FloatLo * pflRadiusNoHP){

	DValNodePtr ** pdnZhangAdjList ;
	DValNodePtr ** pdnBryantAdjList ;
	Int4Ptr piBryantPotential = NULL;
	Int4Ptr piZhangPtnl = NULL;
	DValNodePtr pdnZhangAtmList = NULL, pdnListPmmds;
	DValNodePtr  pdnZHere,pdnBHere,pdnZPot,pdnBPot;
	PALN palnZ,palnB;
	FILE * fpzh, *fpbl, *fpAtmInfo;
	Int2	i,j, numResidue, cnt, resZ,resB;
	FloatLo	* pflPotCrease, * pflPotZh, * pflPotBryant;
	FloatLo  ** ppflPotCrease, ** ppflPotZh, **ppflPotBryant;
	FloatLo flHighTmp, flLowTmp;
	char ftmp[PATH_MAX];

	numResidue	= ((PMMD)(pmsdArray[0]->pdnmmHead->data.ptrvalue))->iResCount;

	/*Memory allocation for this "bunch" of arrays*/
	pdnZhangAdjList		= (DValNodePtr **)MemNew(iNumSimulation * sizeof (DValNodePtr *));
	for(i=0; i< iNumSimulation; i++){
		pdnZhangAdjList[i]= (DValNodePtr*)MemNew(iNumModels * sizeof(DValNodePtr));
	}

	pdnBryantAdjList		= (DValNodePtr **)MemNew(iNumSimulation * sizeof (DValNodePtr *));
	for(i=0; i< iNumSimulation; i++){
		pdnBryantAdjList[i]= (DValNodePtr*)MemNew(iNumModels * sizeof(DValNodePtr));
	}

	ppflPotZh	= (FloatLo **)MemNew(iNumModels * sizeof(FloatLo *));
	for (i=0; i< iNumModels; i++){
		ppflPotZh[i] = (FloatLo *)MemNew(numResidue * sizeof(FloatLo));
	}
	
	ppflPotBryant	= (FloatLo **)MemNew(iNumModels * sizeof(FloatLo *));
	for (i=0; i< iNumModels; i++){
		ppflPotBryant[i] = (FloatLo *)MemNew(numResidue * sizeof(FloatLo));
	}

	ppflPotCrease	= (FloatLo **)MemNew(iNumModels * sizeof(FloatLo *));
	for (i=0; i< iNumModels; i++){
		ppflPotCrease[i] = (FloatLo *)MemNew(numResidue * sizeof(FloatLo));
	}

	pflPotCrease		= (FloatLo   *)MemNew(iNumModels * sizeof (FloatLo));
	pflPotZh			= (FloatLo   *)MemNew(iNumModels * sizeof (FloatLo));
	pflPotBryant		= (FloatLo   *)MemNew(iNumModels * sizeof (FloatLo));

	/*retrieves Zhang potential energy information from file*/
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,ZHANG_POTENTIAL);
	fpzh=FileOpen(ftmp,"r");
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,ZHANG_ATOMS);
	fpAtmInfo=FileOpen(ftmp,"r");

	LoadZhangPotential (&piZhangPtnl, fpzh);
	LoadZhangAtmList (&pdnZhangAtmList, fpAtmInfo);

	/*retrieves Bryant-Lawrence Potential energy information from file*/
	sprintf(ftmp,"%s%s",CFG_local_datafilepath,BL_POTENTIAL);
	fpbl=FileOpen(ftmp,"r");
	LoadBryantPotential (&piBryantPotential, fpbl, TRUE /*bUsingPep*/, FALSE /*bNegateInput*/);

	FileClose (fpzh);
	FileClose (fpbl);
	FileClose (fpAtmInfo);

	/*set the proper informations in order to set the crease potential*/
	if(creaseIncl ==0){
		creaseIncl = numResidue;
	}
	if(creaseExcl > (numResidue -2)){
		creaseExcl = numResidue -2;
	}

	for(i=0; i< iNumSimulation; i++){
		for(j=0; j< iNumModels; j++){
			pdnListPmmds = CreateListPmmds (NULL, pmsdArray[i]);

			ComputeZhangPotential (&pdnZhangAdjList[i][j], pdnListPmmds,j+1, FALSE/*bInclusiveWindow*/,ZHANG_WINDOWSIZE/*iWindowSize*/, piZhangPtnl,pdnZhangAtmList,FALSE);
			pflPotZh[j] += GetTotalPotential();

			ComputeBryantPotential (&pdnBryantAdjList[i][j], pdnListPmmds,j+1, FALSE/*bInclusiveWindow*/, BRYANT_WINDOWSIZE/*iWindowSize*/, piBryantPotential,TRUE/*bUsingPep*/, FALSE);
			/* returns negative of correct value */
			pflPotBryant[j] += GetTotalPotential();

			if (CalcCreaseEnergy((PMMD)((pmsdArray[i]->pdnmmHead)->data.ptrvalue),UNITS_BRYANT,creaseIncl,creaseExcl,creaseDecay,FALSE,j+1)!=ERR_SUCCESS) {
				ErrPostEx(SEV_ERROR,0,10,"An error occurred in CalcCreaseEnergy");
				pflPotCrease[j]=9999999.9;}
			else{
				pflPotCrease[j]+=GetCreaseEnergy();
				for (cnt=0;cnt<numResidue;cnt++) {
					ppflPotCrease[j][cnt]	+= GetCreaseRes(cnt+1);
				}
			}
			FreeListPmmds (&pdnListPmmds);
		}
	}

	/*retrieve the potentials from the adjList */
	for(i=0; i< iNumSimulation; i++){
		for(j=0; j<iNumModels; j++){
			pdnZHere=pdnZhangAdjList[i][j];
			pdnBHere=pdnBryantAdjList[i][j];
			while (pdnZHere!=NULL) {
				/* get potentials from pdnAdjLists */
				pdnZPot=(DValNodePtr)(pdnZHere->data.ptrvalue);
				palnZ=(PALN)(pdnZPot->data.ptrvalue);
				resZ=palnZ->pmgd->pdnmgLink->choice-1;
				ppflPotZh[j][resZ]		+= palnZ->potential;					
				pdnZHere=pdnZHere->next;
			}
			while (pdnBHere!=NULL) {
				pdnBPot=(DValNodePtr)(pdnBHere->data.ptrvalue);
				palnB=(PALN)(pdnBPot->data.ptrvalue);
				resB=palnB->pmgd->pdnmgLink->choice-1;
				ppflPotBryant[j][resB]	+= palnB->potential;
				pdnBHere=pdnBHere->next;
			}		
		}
	}

	/*average all the total energies*/
	for(j=0; j< iNumModels; j++){
		pflPotZh[j]		= pflPotZh[j]     / (FloatLo) iNumSimulation;
		pflPotBryant[j]	= pflPotBryant[j] / (FloatLo) iNumSimulation; 
		pflPotCrease[j]	= pflPotCrease[j] / (FloatLo) iNumSimulation;
	}
	/*average all the residue energies*/
	for(i= 0; i< iNumModels; i++){
		for(cnt =0; cnt < numResidue; cnt++){
			ppflPotZh[i][cnt]		= ppflPotZh    [i][cnt] / (FloatLo) iNumSimulation;
			ppflPotBryant[i][cnt]	= ppflPotBryant[i][cnt] / (FloatLo) iNumSimulation;
			ppflPotCrease[i][cnt]	= ppflPotCrease[i][cnt] / (FloatLo) iNumSimulation;		
		}
	}

	
	/*graphing and plotting*/
	/*printout the Energy files to the screen and file*/
	if(iNumSimulation ==1){
		sprintf(ftmp, "%s_PotZh_eachRes%s", fileName,".dat");	
		Printppfl(ppflPotZh, iNumModels, numResidue,  "model #", "residue #", ftmp);
		sprintf(ftmp, "%s_PotBryant_eachRes%s", fileName,".dat");	
		Printppfl(ppflPotBryant, iNumModels, numResidue,  "model #", "residue #",ftmp);
		sprintf(ftmp, "%s_PotCrease_eachRes%s", fileName,".dat");	
		Printppfl(ppflPotCrease, iNumModels, numResidue,  "model #", "residue #",ftmp);

		sprintf(ftmp, "%s_PotZh_totStruct%s", fileName,".dat");	
		Printfl(pflPotZh, iNumModels,"PotZh_totStruct", ftmp);
		sprintf(ftmp, "%s_PotBryant_totStruct%s", fileName,".dat");	
		Printfl(pflPotBryant, iNumModels, "PotBryant_totStruct",ftmp);
		sprintf(ftmp, "%s_PotCrease_totStruct%s", fileName,".dat");	
		Printfl(pflPotCrease, iNumModels, "PotCrease_totStruct",ftmp);
	}
	else{
		sprintf(ftmp, "%s_Avrg_PotZh_eachRes_upto_%d%s", fileName,iNumSimulation,".dat");	
		Printppfl(ppflPotZh, iNumModels, numResidue,  "model #", "residue #",ftmp);
		sprintf(ftmp, "%s_Avrg_PotBryant_eachRes_upto_%d%s", fileName,iNumSimulation,".dat");	
		Printppfl(ppflPotBryant, iNumModels, numResidue,  "model #", "residue #", ftmp);
		sprintf(ftmp, "%s_Avrg_PotCrease_eachRes_upto_%d%s", fileName,iNumSimulation,".dat");	
		Printppfl(ppflPotCrease, iNumModels, numResidue,  "model #", "residue #", ftmp);

		sprintf(ftmp, "%s_Avrg_PotZh_totStruct_upto_%d%s", fileName,iNumSimulation,".dat");	
		Printfl(pflPotZh, iNumModels, "Avrg_PotZh_totStruct",ftmp);
		sprintf(ftmp, "%s_Avrg_PotBryant_totStruct_upto_%d%s", fileName,iNumSimulation,".dat");	
		Printfl(pflPotBryant, iNumModels,"Avrg_PotBryant_totStruct" ,ftmp);
		sprintf(ftmp, "%s_Avrg_PotCrease_totStruct_upto_%d%s", fileName,iNumSimulation,".dat");	
		Printfl(pflPotCrease, iNumModels,"Avrg_PotCrease_totStruct" ,ftmp);
	}

	/*Draw the Energy graphs*/
	if(iNumSimulation ==1){

		sprintf(ftmp, "%s_PotZh_eachRes%s", fileName,".gif");	
		GraphArray(ftmp, NULL, NULL, i, 1, 0, numResidue, "Zhang Potential", "", "RESIDUE #", 
					"Time", 645, 500, ENERGY_PLOTTING,ppflPotZh);
		sprintf(ftmp, "%s_PotBryant_eachRes%s", fileName,".gif");	
		GraphArray(ftmp, NULL, NULL, i, 1, 0, numResidue, "Bryant Potential", "", "RESIDUE #", 
					"Time", 645, 500, ENERGY_PLOTTING,ppflPotBryant);
		sprintf(ftmp, "%s_PotCrease_eachRes%s", fileName,".gif");	
		GraphArray(ftmp, NULL, NULL, i, 1, 0, numResidue, "Crease Potential", "", "RESIDUE #", 
					"Time", 645, 500, ENERGY_PLOTTING,ppflPotCrease);

		flLowTmp  = findSmallestFloat(pflPotZh, iNumModels);
		flHighTmp = findGreatestFloat(pflPotZh, iNumModels);
		sprintf(ftmp, "%s_PotZh_totStruct%s", fileName,".gif");	
		GraphArray(ftmp, pflPotZh, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Z-D potential", "", "Z-D potential", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		flLowTmp  = findSmallestFloat(pflPotBryant, iNumModels);
		flHighTmp = findGreatestFloat(pflPotBryant, iNumModels);
		sprintf(ftmp, "%s_PotBryant_totStruct%s", fileName,".gif");	
		GraphArray(ftmp, pflPotBryant, NULL, iNumModels, 1, flLowTmp, flHighTmp, "B-L potential", "", "B-L potential", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		flLowTmp  = findSmallestFloat(pflPotCrease, iNumModels);
		flHighTmp = findGreatestFloat(pflPotCrease, iNumModels);
		sprintf(ftmp, "%s_PotCrease_totStruct%s", fileName,".gif");	
		GraphArray(ftmp, pflPotCrease, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Crease Energy", "", "Crease Energy", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);


		if(pflRmsdArray){
		    /*draw the radius of gyration against the RMSD*/
		    /*First using the HP only*/
		    flLowTmp  = findSmallestFloat(pflRadiusHP, iNumModels);
		    flHighTmp = findGreatestFloat(pflRadiusHP, iNumModels);
		    sprintf(ftmp, "%s_RGYR_HP_VS_RMSD_%s", fileName,".gif");	
		    GraphArray(ftmp, pflRadiusHP, pflRmsdArray, iNumModels, 1, flLowTmp, flHighTmp, "RGYR_VS_RMSD_HP", "", "Rgyr (HP)", 
					    "RMSD", 645, 500, RGYR_VS_RMSD, pflPotBryant);
		    /*Then using All the forces*/
		    flLowTmp  = findSmallestFloat(pflRadiusNoHP, iNumModels);
		    flHighTmp = findGreatestFloat(pflRadiusNoHP, iNumModels);
		    sprintf(ftmp, "%s_RGYR_NOHP_VS_RMSD_%s", fileName,".gif");	
		    GraphArray(ftmp, pflRadiusNoHP, pflRmsdArray, iNumModels, 1, flLowTmp, flHighTmp, "RGYR_VS_RMSD_NoHP", "", "Rgyr", 
					    "RMSD", 645, 500, RGYR_VS_RMSD, pflPotBryant);
		}

	}

	else{
		sprintf(ftmp, "%s_Avrg_PotZh_eachRes_upto_%d%s", fileName,iNumSimulation,".gif");	
		GraphArray(ftmp, NULL, NULL, i, 1, 0, numResidue, "Zhang Potential", "", "RESIDUE #", 
					"Time", 645, 500, ENERGY_PLOTTING,ppflPotZh);
		sprintf(ftmp, "%s_Avrg_PotBryant_eachRes_upto_%d%s", fileName,iNumSimulation,".gif");	
		GraphArray(ftmp, NULL, NULL, i, 1, 0, numResidue, "Bryant Potential", "", "RESIDUE #", 
					"Time", 645, 500, ENERGY_PLOTTING,ppflPotBryant);
		sprintf(ftmp, "%s_Avrg_PotCrease_eachRes_upto_%d%s", fileName,iNumSimulation,".gif");	
		GraphArray(ftmp, NULL, NULL, i, 1, 0, numResidue, "Crease Potential", "", "RESIDUE #", 
					"Time", 645, 500, ENERGY_PLOTTING,ppflPotCrease);

		flLowTmp  = findSmallestFloat(pflPotZh, iNumModels);
		flHighTmp = findGreatestFloat(pflPotZh, iNumModels);
		sprintf(ftmp, "%s_Avrg_PotZh_totStruct_upto_%d%s", fileName,iNumSimulation,".gif");	
		GraphArray(ftmp, pflPotZh, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Z-D potential", "", "Z-D potential", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		flLowTmp  = findSmallestFloat(pflPotBryant, iNumModels);
		flHighTmp = findGreatestFloat(pflPotBryant, iNumModels);
		sprintf(ftmp, "%s_Avrg_PotBryant_totStruct_upto_%d%s", fileName,iNumSimulation,".gif");	
		GraphArray(ftmp, pflPotBryant, NULL, iNumModels, 1, flLowTmp, flHighTmp, "B-L potential", "", "B-L potential", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		flLowTmp  = findSmallestFloat(pflPotCrease, iNumModels);
		flHighTmp = findGreatestFloat(pflPotCrease, iNumModels);
		sprintf(ftmp, "%s_Avrg_PotCrease_totStruct_upto_%d%s", fileName,iNumSimulation,".gif");	
		GraphArray(ftmp, pflPotCrease, NULL, iNumModels, 1, flLowTmp, flHighTmp, "Crease Energy", "", "Crease Energy", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);
		
		if(pflRmsdArray){
		    /*draw the radius of gyration against the RMSD*/
		    /*First using the HP only*/
		    flLowTmp  = findSmallestFloat(pflRadiusHP, iNumModels);
		    flHighTmp = findGreatestFloat(pflRadiusHP, iNumModels);
		    sprintf(ftmp, "%s_RGYR_HP_VS_RMSD_upto_%d%s", fileName,iNumSimulation,".gif");	
		    GraphArray(ftmp, pflRadiusHP, pflRmsdArray, iNumModels, 1, flLowTmp, flHighTmp, "RGYR_VS_RMSD_HP", "", "Rgyr (HP)", 
					    "RMSD", 645, 500, RGYR_VS_RMSD, pflPotBryant);
		    /*Then using All the forces*/
		    flLowTmp  = findSmallestFloat(pflRadiusNoHP, iNumModels);
		    flHighTmp = findGreatestFloat(pflRadiusNoHP, iNumModels);
		    sprintf(ftmp, "%s_RGYR_NOHP_VS_RMSD_upto_%d%s", fileName,iNumSimulation,".gif");	
		    GraphArray(ftmp, pflRadiusNoHP, pflRmsdArray, iNumModels, 1, flLowTmp, flHighTmp, "RGYR_VS_RMSD_NoHP", "", "Rgyr", 
					    "RMSD", 645, 500, RGYR_VS_RMSD, pflPotBryant);
        }
	}

	/*Memory freeing*/
	/*free memory of structure retrieving residue energy*/
	for(i=0; i< iNumSimulation; i++){
		/*for(j=0; j< iNumModels; j++){
			FreeAdjList (pdnZhangAdjList[i][j]);
			FreeAdjList (pdnBryantAdjList[i][j]);
		}*/
		MemFree(pdnZhangAdjList[i]);
		MemFree(pdnBryantAdjList[i]);
	}
	MemFree(pdnZhangAdjList);
	MemFree(pdnBryantAdjList);
	FreeCreaseEnergy();

	for(i=0; i<iNumModels; i++){
		MemFree(ppflPotZh[i]);
		MemFree(ppflPotBryant[i]);
		MemFree(ppflPotCrease[i]);
	}
	MemFree(ppflPotZh);
	MemFree(ppflPotBryant);
	MemFree(ppflPotCrease);

	/*free energy of structure retrieving total energy on structure*/
	MemFree(pflPotCrease);
	MemFree(pflPotZh);
	MemFree(pflPotBryant);

	FreeZhangPotential (&piZhangPtnl);
	FreeBryantPotential (&piBryantPotential);

	FreeZhangAtmList (&pdnZhangAtmList);
	FreeZhangPotential (&piZhangPtnl);
	FreeBryantPotential (&piBryantPotential);
}

/*Helper function of ComputeNOES.
	returns true if two residues considered in contact.
	return false if not.
*/
Boolean IsThereAContact ( Int2 iResidue1, Int2 iResidue2, Int2 iAtom1, Int2 iAtom2, PMSD currentPmsd,Int2 iModel){
	vec vecAtom1, vecAtom2, vecAtom1alt, vecAtom2alt, vecRes, vecZero;
	PDNMG	pdnmgTemp;
	Int2 i,j,fnd=0;
	FloatLo distance;
	Boolean alt1=FALSE,alt2=FALSE;
/*	Boolean exception1, exception2;

	exception1 = FALSE;
	exception2 = FALSE;*/

	vecZero[0]=0;
	vecZero[1]=0;
	vecZero[2]=0;

	/*go to the first residue in the list*/
	pdnmgTemp = ((PMMD)((PDNMM)currentPmsd->pdnmmHead)->data.ptrvalue)->pdnmgHead;


	j=MAX(iResidue1, iResidue2);
/*	if (j==100){
		if (iResidue1 ==100 ){
				exception1 = TRUE;
				iResidue1=36;
				j=MAX(iResidue1, iResidue2);
		}
		else{
				exception2 = TRUE;
				iResidue2=36;
				j=MAX(iResidue1, iResidue2);
		}
	}*/

	for(i=1; i<= j; i++){

		if (pdnmgTemp->choice== iResidue1){
			fnd++;
			alt1=FALSE;
			if (iAtom1==0)
				GetCoOrds(pdnmgTemp->data.ptrvalue, " H  ", vecZero,vecAtom1, iModel);
			else if (iAtom1==4) {
				GetCoOrds(pdnmgTemp->data.ptrvalue, "1HD2", vecZero,vecAtom1, iModel);  /* ASN delta H */
				GetCoOrds(pdnmgTemp->data.ptrvalue, "2HD2", vecZero,vecAtom1alt, iModel);  /* ASN delta H */
				alt1=TRUE;
			}
			else if (iAtom1==5) {
				if ((((PMGD)(pdnmgTemp->data.ptrvalue))->pcIUPAC)[0]=='W')
					GetCoOrds(pdnmgTemp->data.ptrvalue, " HE1", vecZero,vecAtom1, iModel); /* TRP epsilon H */
				else {
					GetCoOrds(pdnmgTemp->data.ptrvalue, "1HE2", vecZero,vecAtom1, iModel); /* GLN epsilon H */
					GetCoOrds(pdnmgTemp->data.ptrvalue, "2HE2", vecZero,vecAtom1alt, iModel); /* GLN epsilon H */
					alt1=TRUE;
				}
			}
/*			if(!exception1)
				GetCoOrds(pdnmgTemp->data.ptrvalue, " H  ", vecZero,vecAtom1, iModel);
			else
				GetCoOrds(pdnmgTemp->data.ptrvalue, " HE1", vecZero,vecAtom1, iModel);*/
		}
		else if (pdnmgTemp->choice==iResidue2){
			fnd++;
			alt2=FALSE;
			if (iAtom2==0)
				GetCoOrds(pdnmgTemp->data.ptrvalue, " H  ", vecZero,vecAtom2, iModel);
			else if (iAtom2==4) {
				GetCoOrds(pdnmgTemp->data.ptrvalue, "1HD2", vecZero,vecAtom2, iModel);  /* ASN delta H */
				GetCoOrds(pdnmgTemp->data.ptrvalue, "2HD2", vecZero,vecAtom2alt, iModel);  /* ASN delta H */
				alt2=TRUE;
			}
			else if (iAtom2==5) {
				if ((((PMGD)(pdnmgTemp->data.ptrvalue))->pcIUPAC)[0]=='W')
					GetCoOrds(pdnmgTemp->data.ptrvalue, " HE1", vecZero,vecAtom2, iModel); /* TRP epsilon H */
				else {
					GetCoOrds(pdnmgTemp->data.ptrvalue, "1HE2", vecZero,vecAtom2, iModel); /* GLN epsilon H */
					GetCoOrds(pdnmgTemp->data.ptrvalue, "2HE2", vecZero,vecAtom2alt, iModel); /* GLN epsilon H */
					alt2=TRUE;
				}
			}
/*			if(!exception2)
				GetCoOrds(pdnmgTemp->data.ptrvalue," H  ", vecZero,vecAtom2, iModel);
			else
				GetCoOrds(pdnmgTemp->data.ptrvalue," HE1", vecZero,vecAtom2, iModel);*/
		}
		pdnmgTemp = pdnmgTemp->next;
	}
	if (fnd<2)
		return FALSE;

	VecSub(vecRes, vecAtom1, vecAtom2);
	distance = getMag(vecRes);

	if (distance < NOE_CUTOFF)
		return TRUE;
	/* check alternate atom pairs if exist */
	if (alt1) {
		VecSub(vecRes, vecAtom1alt, vecAtom2);
		distance = getMag(vecRes);
		if (distance < NOE_CUTOFF)
			return TRUE;
	}
	if (alt2) {
		VecSub(vecRes, vecAtom1, vecAtom2alt);
		distance = getMag(vecRes);
		if (distance < NOE_CUTOFF)
			return TRUE;
	}
	if (alt1 && alt2) {
		VecSub(vecRes, vecAtom1alt, vecAtom2alt);
		distance = getMag(vecRes);
		if (distance < NOE_CUTOFF)
			return TRUE;
	}
	return FALSE;		
}


void ComputeNOES(PMSD * pmsdArray,Int2 iNumSimulation, Int2 iNumModels,Int2 iNumPartition, CharPtr fileName){
	/*arrays that contain all the residues that we want to compare*/
/*	Int2 pi2Residue1[NUM_CONTACTS_IN_DRK] = {	10,14,18,24,27,34,42,46,50,13,
												38,43,45,14,29,42,40, 3,40,34,
												38,27,37,27,24,36,38,14,42,22,
												100,12,37, 9,18,20,34,40,18,34,
												37,18,10,14,18,29, 9, 8,11,14,
												18,20,10,12,18,12,12,10,12
												};	

	Int2 pi2Residue2[NUM_CONTACTS_IN_DRK] = {	14,18,22,28,31,38,46,50,54,18,
												43,48,50,20,35,48,46,10,47,42,
												46,36,46,38,100,48,50,27,55,100,
												50,28,53,27,100,38,52,58,37,55,
												58,40,33,38,40,53,100,100,40,45,
												50,53,45,48,55,51,53,53,55
											};*/
	Int2 pi2Residue1[NUM_CONTACTS_IN_DRK] = {	8,9,12,13,14,15,19,22,29,29,31,
												32,33,34,35,35,36,36,40,48,50,
												51,54,55,12,13,14,25,29,31,32,
												34,35,36,39,42,50,51,13,29,
												30,31,36,29,29,29,39,37,22,
												23,18
											};	

	Int2 pi2Residue2[NUM_CONTACTS_IN_DRK] = {	11,12,15,16,17,18,22,25,32,32,34,
												35,36,37,38,38,39,39,43,51,53,
												54,57,58,16,17,18,29,33,35,
												36,38,39,40,43,46,54,55,18,34,
												35,36,41,35,35,36,46,46,36,
												37,36
											};
	/* 0 = backbone amide, 4=delta N, 5=epsilon N */
	Int2 pi2Residue1Atom[NUM_CONTACTS_IN_DRK] = {	0,0,0,0,0,0,0,0,0,4,0,
													0,0,0,0,4,0,5,0,0,0,
													4,0,0,0,0,0,0,4,0,0,
													0,4,5,0,0,0,4,0,4,
													0,0,5,4,0,0,0,0,0,
													5,0
												};	

	Int2 pi2Residue2Atom[NUM_CONTACTS_IN_DRK] = {	0,0,0,0,0,0,0,0,0,0,0,
													4,5,0,0,0,0,0,0,4,0,
													0,4,0,0,0,0,4,0,4,
													5,0,0,0,0,0,0,0,0,0,
													4,5,0,0,4,0,0,0,5,
													0,5
												};

	Boolean *** bNOEContacts;
	FloatLo	**	ppflOccupancy;	/*Contains the contact occupancy for each partition*/
	FloatLo	*	pflNumNOEPerFrame;

	char pcTemp[PATH_MAX];
	Int2 i,j,k, temp, temp2;
	FloatLo flLowTmp, flHighTmp;
	
	/*Memory allocation*/
	bNOEContacts = (Boolean ***) MemNew (iNumSimulation * sizeof (Boolean **));
	for (i=0; i<iNumSimulation; i++){
		bNOEContacts[i] = (Boolean **) MemNew (iNumModels * sizeof ( Boolean *));
		for (j=0; j < iNumModels; j++){
			bNOEContacts[i][j] = (Boolean*) MemNew (NUM_CONTACTS_IN_DRK * sizeof (Boolean));
		}
	}
	
	ppflOccupancy = (FloatLo **) MemNew ( iNumPartition * sizeof(FloatLo *));
	for (i=0; i< iNumPartition; i++){
		ppflOccupancy[i] = (FloatLo *) MemNew ( NUM_CONTACTS_IN_DRK * sizeof (FloatLo));
	}

	pflNumNOEPerFrame = (FloatLo *) MemNew ( iNumModels * sizeof (FloatLo));

	/*set temp to the number of frames per partition temp2 is for the last partition*/
	temp  = ceil((FloatLo) iNumModels/ (FloatLo)iNumPartition );
	temp2 = (FloatLo)(iNumModels - ((iNumPartition -1 ) * (ceil((FloatLo) iNumModels/ (FloatLo)iNumPartition ))));
	
	/*go through the models and simulations in order to retrieve the contacts*/
	for(i=0; i<iNumSimulation; i++){
		for(j=0; j<iNumModels; j++){
			for(k=0; k<NUM_CONTACTS_IN_DRK; k++){
				bNOEContacts[i][j][k] = IsThereAContact(pi2Residue1[k], pi2Residue2[k], pi2Residue1Atom[k], pi2Residue2Atom[k], pmsdArray[i],j+1);
				if (bNOEContacts[i][j][k] == TRUE){
					pflNumNOEPerFrame[j] ++;
					if( (j / temp) < iNumPartition -1)
						ppflOccupancy[j / temp][k] ++;
					else
						ppflOccupancy[iNumPartition -1][k] ++;
				}
			}
		}
	}

	/*Average the total number of NOE's*/
	for(i=0; i< iNumModels; i++){
		pflNumNOEPerFrame[i] = pflNumNOEPerFrame[i] / iNumSimulation;
	}

	for(i=0; i< iNumPartition; i++){
		for (j=0; j< NUM_CONTACTS_IN_DRK; j++){
			if(i !=iNumPartition -1){
				ppflOccupancy[i][j] = ppflOccupancy[i][j] / (iNumSimulation * temp)* 100;
			}
			else{
				ppflOccupancy[i][j] = ppflOccupancy[i][j] / (iNumSimulation * temp2)* 100;
			}
		}
	}

	/*editing on file/screen and graphing*/
	if(iNumSimulation ==1){
		sprintf(pcTemp, "%s_NOE_Per_Frame_%s", fileName,  ".dat");	
		Printfl(pflNumNOEPerFrame, iNumModels,"NOE per FRAME", pcTemp);

		for(i=0; i< iNumPartition; i++){
			sprintf(pcTemp, "%s_NOE_Occupancy_part%d%s", fileName,i ,".dat");	
			Printfl(ppflOccupancy[i], NUM_CONTACTS_IN_DRK, "Occupancy per contact in this partition", pcTemp);
		}
	}
	else{
		sprintf(pcTemp, "%s_Avrg_NOE_Per_Frame%s%d%s",fileName,"_upto",iNumSimulation,".dat");
		Printfl(pflNumNOEPerFrame, iNumModels,"NOE per FRAME", pcTemp);

		for(i=0; i< iNumPartition; i++){
			sprintf(pcTemp, "%s_Avrg_NOE_Occupancy%s%d_part%d%s", fileName,"_upto",iNumSimulation,i ,".dat");	
			Printfl(ppflOccupancy[i], NUM_CONTACTS_IN_DRK, "Occupancy per contact in this partition", pcTemp);
		}
	}

	/*Draw the #NOE per frame as well as the Occupancy of each contact over a certain partition*/
	if(iNumSimulation ==1){
		flLowTmp  = 0.0;
		flHighTmp = findGreatestFloat(pflNumNOEPerFrame, iNumModels);
		sprintf(pcTemp, "%s_NOE_per_Frame%s", fileName, ".gif");	
		GraphArray(pcTemp, pflNumNOEPerFrame, NULL, iNumModels, 1, flLowTmp, flHighTmp, 
			"NOE/Frame", "", "NOE/Frame", 
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		for(i=0; i<iNumPartition; i++){
			flHighTmp = 100;
			sprintf(pcTemp, "%s_NOE_Occupancy_part%d%s", fileName,i ,".gif");	
			GraphArray(pcTemp, ppflOccupancy[i], NULL, (Int2)NUM_CONTACTS_IN_DRK, 1, flLowTmp, flHighTmp, 
						"Occupancy of each contact", "", "Occupancy (%)", 
						"Contact Number", 645, 500, RECTANGLE_VS_CONTACT, NULL);
		}
	}

	else{
		flLowTmp  = 0.0;
		flHighTmp = findGreatestFloat(pflNumNOEPerFrame, iNumModels);
		sprintf(pcTemp, "%s_Avrg_NOE_per_Frame%s%d%s", fileName,"_upto",iNumSimulation,  ".gif");	
		GraphArray(pcTemp, pflNumNOEPerFrame, NULL, iNumModels, 1, flLowTmp, flHighTmp, "NOE/Frame", "", "NOE per Frame",
					"Time", 645, 500, ARRAY_VS_TIME, NULL);

		for(i=0; i<iNumPartition; i++){
			flHighTmp = 100;
			sprintf(pcTemp, "%s_Avrg_NOE_Occupancy%s%d_part%d%s", fileName,"_upto",iNumSimulation,i ,".gif");	
			GraphArray(pcTemp, ppflOccupancy[i], NULL, NUM_CONTACTS_IN_DRK, 1, flLowTmp, flHighTmp, 
						"Occupancy of each contact", "", "Occupancy (%)", 
						"Contact Number", 645, 500, RECTANGLE_VS_CONTACT, NULL);
		}
	}
	
	/*Memory freeing*/
	for(i=0; i< iNumSimulation; i++){
		for(j=0; j < iNumModels; j++){
			MemFree (bNOEContacts[i][j] );
		}
		MemFree(bNOEContacts[i]);
	}
	MemFree(bNOEContacts);

	for(i=0; i<iNumPartition; i++){
		MemFree(ppflOccupancy[i]);
	}
	MemFree(ppflOccupancy);

	MemFree(pflNumNOEPerFrame);
}


static Int2 RetrieveInfo(CharPtr fileName,Boolean remoteFile, Int2 mLevel,
						 Int2 * modelNum,  Boolean output, Boolean createGif,
						 int numFile, Boolean average, Int2 numPartition,
						 Int4 graphWidth, Int4 graphHeight, Int4 creaseIncl,
						 Int4 creaseExcl, Boolean creaseDecay, CharPtr pcnatfnam)
{
	PMSD pmsdRoot, pmsdNative;
	PMSD * pmsdArray;
	PMMD pmmdThis;
	PMMD pmmdNative = NULL;
	PDNMM pdnmmHere;
	Char pcTemp[PATH_MAX];

	/*Array used to store the data before plotting*/
	FloatLo *	pflRadiusHP ;
	FloatLo *	pflRadiusNoHP;
	FloatLo *	pflRmsdArray = NULL;
	Int4	**	ppi4Surfacc;
	FloatLo	*	pflSurfaccTot;
	Int4	*	pi4Temp;
	CharPtr *   ppc2ndrystruct;


	/*2-d array used to plot multiple movie's data on the same diagram*/	
	FloatLo **	pppflRadiusHP;	
	FloatLo **	pppflRadiusNoHP;
	FloatLo **  pppflRmsdArray = NULL;
	CharPtr **  pppc2ndryStruct;

	/*2-d array to store the %age of each secondary structure*/
	FloatLo **	ppfl2ndryStruct;

	/* array used to retrieve the structure containing the minimum 
	 * number of models (used when calculating averages)*/
    Int4	*	iArrayNumModels;
    FloatLo flLow1,  flHigh1, flLow2, flHigh2;
    Int4 models,i,j,k, iMinNumModels;
    Char pNames1[3][25]={"% helix (avg)","% sheet (avg)","% coil and other (avg)"};
    Char pNames2[3][25]={"% helix","% sheet","% coil"};

    if (pcnatfnam != NULL){
        /*opens the native structure required to trace the rmsd*/

	if (FileLength(pcnatfnam) == 0) 
        {
		printf("\nInput native structure file not found: %s", pcnatfnam);	
		/* error occurred */	
		return 3;
        }
        printf("\nNative structure loading ...filename: %s\n", pcnatfnam);
        pmsdNative=LoadABiostruc(pcnatfnam,remoteFile,mLevel,modelNum);
        if (pmsdNative==NULL) {
            /* error occurred */
            printf("\nNative structure loading error from file: %s\n", pcnatfnam);
            return 3;
        }
        printf("Native structure loaded: %s.\n", pcnatfnam);

        if (!IsStructureNode(pmsdNative)){
            printf("\nBad internal data type returned from native structure load %s\n",pcnatfnam);
            return ERR_FAIL;
        }
        pmmdNative = pmsdNative->pdnmmHead->data.ptrvalue;
    }

	/*we want to trace the average of many models*/
	if(average){

			pmsdArray = (PMSD*)MemNew( numFile * sizeof(PMSD));

			/*We need to load a biostruct first to retrieve the number of models per structure*/
			sprintf(pcTemp, "%s%s%d%s",  fileName,"_movie_",1, ".val");
			if (FileLength(pcTemp) == 0) 
                        {
				printf("\nInput sequential movie structure file not found: %s", pcTemp);	
				/* error occurred */
				return 3;
	                 }

			printf("\nStructure loading ...filename: %s\n", pcTemp);
			pmsdArray[0]=LoadABiostruc(pcTemp,remoteFile,mLevel,modelNum);
			if (pmsdArray[0]==NULL) {
			printf("\nStructure failed to load: %s\n", pcTemp);
				/* error occurred */
				return 3;
			}
			printf("Structure # %d loaded.\n", 1);
			if (!IsStructureNode(pmsdArray[0])){
				printf("\nBad internal data type returned from structure file load: %s\n",pcTemp);
				return ERR_FAIL;
			}

			pdnmmHere = pmsdArray[0]->pdnmmHead;
			pmmdThis = (PMMD)pdnmmHere->data.ptrvalue;

			/*memory allocation for all the data arrays*/
			pflRadiusHP		= (FloatLo*) MemNew( (NUM_MODEL+1) * sizeof(FloatLo));
			pflRadiusNoHP	= (FloatLo*) MemNew( (NUM_MODEL+1) * sizeof(FloatLo));
			if(pmmdNative){
			    pflRmsdArray	= (FloatLo*) MemNew( (NUM_MODEL+1) * sizeof(FloatLo));
			}
			ppi4Surfacc		= (Int4 **)  MemNew( (NUM_MODEL+1) * sizeof (Int4*) );  

			pppc2ndryStruct = (CharPtr**)MemNew( numFile   * sizeof(CharPtr*));
			ppfl2ndryStruct = (FloatLo**)MemNew( NUM_GROUP_2NDRY_STRUCT * sizeof(FloatLo*) );
			ppc2ndrystruct	= (CharPtr *)MemNew( (NUM_MODEL+1) * sizeof(CharPtr));
			iArrayNumModels = (Int4*) MemNew(numFile * sizeof(Int4));

			for(k=0; k <numFile ; k++){
				pppc2ndryStruct[k] = (CharPtr*)MemNew((NUM_MODEL+1) * sizeof(CharPtr));
				for(i=0; i<=NUM_MODEL; i++){
					pppc2ndryStruct[k][i] = (CharPtr) MemNew( (pmmdThis->iResCount * sizeof(char) ) +1 );
				}
			}

			for(i=0; i<NUM_GROUP_2NDRY_STRUCT; i++){
				ppfl2ndryStruct[i] = (FloatLo*) MemNew( pmsdArray[0]->iModels * sizeof(FloatLo) );
			}

			for(i=0; i<=NUM_MODEL; i++){
				ppc2ndrystruct[i] = (CharPtr) MemNew( (pmmdThis->iResCount * sizeof(char) ) +1 );
			}
			
			/*fill the array with the first movie*/
			for (models=1;models<=pmsdArray[0]->iModels;models++)
			{
				pflRadiusHP  [models -1] = GetExactRgyr(pmmdThis, models, TRUE);
				pflRadiusNoHP[models -1] = GetExactRgyr(pmmdThis, models, FALSE);
				
				if(pmmdNative){
				    Align2StrucSVD(	pmmdNative,(PMMD)((pmsdArray[0]->pdnmmHead)->data.ptrvalue),
									ALIGN_CA, 1, pmmdNative->iResCount, 1,models);
				    pflRmsdArray[models -1]	= GetRMSD(pmmdNative, (PMMD)((pmsdArray[0]->pdnmmHead)->data.ptrvalue), 
												    ALIGN_CA, 1, pmmdNative->iResCount, 1,models);
                }
				ppi4Surfacc[models -1] = CalcDSSPAccSurf((PMMD)pmmdThis, 2, models);
				CalcDSSPAssignEx(pmmdThis, models , pppc2ndryStruct[0][models -1],FALSE);
			}

			/*fill the array with the rest of the movies*/
			for(i=2; i <= numFile; i++)
			{
				/* load an ASN.1 Biostruc */
				sprintf(pcTemp, "%s%s%d%s",  fileName,"_movie_",i, ".val");

				if (FileLength(pcTemp) == 0) 
	                        {
					printf("\nInput sequential movie structure file not found: %s", pcTemp);	
					/* error occurred */
					return 3;
       		                 }
				printf("\nStructure loading ...filename: %s", pcTemp);

				pmsdArray[i -1]=LoadABiostruc(pcTemp,remoteFile,mLevel,modelNum);
				if (pmsdArray[i -1]==NULL) {
          				printf("\nStructure FAILED to load with assumed sequential filename: %s", pcTemp);
  				/* error occurred */
					return 3;
				}
				printf("\nStructure # %d loaded.\n", i);
			
				if (!IsStructureNode(pmsdArray[i -1])){
					printf("\nInternal pointer failure, cannot continue.\n");
					return ERR_FAIL;
				}

				pdnmmHere = pmsdArray[i -1]->pdnmmHead;
				pmmdThis = (PMMD)pdnmmHere->data.ptrvalue;

				for (models=1;models<=pmsdArray[i -1]->iModels;models++)
				{
					/*fill the array for the radius of gyration and add up the radius*/
					pflRadiusHP  [models -1] += GetExactRgyr(pmmdThis, models, TRUE);
					pflRadiusNoHP[models -1] += GetExactRgyr(pmmdThis, models, FALSE);
					
					if(pmmdNative){			
					/*fill the array for the RMSD and add up the rmsds*/
					    Align2StrucSVD(	pmmdNative,(PMMD)((pmsdArray[i-1]->pdnmmHead)->data.ptrvalue),
										ALIGN_CA, 1, pmmdNative->iResCount, 1,models );
					    pflRmsdArray[models -1] += GetRMSD(pmmdNative, (PMMD)((pmsdArray[i-1]->pdnmmHead)->data.ptrvalue), 
                                                            ALIGN_CA, 1, pmmdNative->iResCount, 1,models);
                    }
                    
					CalcDSSPAssignEx(pmmdThis, models , pppc2ndryStruct[i-1][models -1],FALSE);
					/*retreive the surface accessibility array and add it*/
					pi4Temp = CalcDSSPAccSurf((PMMD)pmmdThis, 2, models);
					if (ppi4Surfacc[models -1]!=NULL) {
						for(k=0; k < pmmdThis->iResCount; k++)
							ppi4Surfacc[models -1][k] += pi4Temp[k];
					}
					MemFree(pi4Temp);
				}
			}

			for(i=0; i< numFile; i++){
				iArrayNumModels[i] = pmsdArray[i]->iModels;
			}

			iMinNumModels = findSmallestInt(iArrayNumModels, numFile);
			if (iMinNumModels>100)
				iMinNumModels=100;
			printf("\nThe minimum number of models is: %d\n",  iMinNumModels);

			/*calculate the average*/
			for (models=1;models<=iMinNumModels;models++)
			{
				pflRadiusHP  [models -1] = pflRadiusHP  [models -1]/numFile;
				pflRadiusNoHP[models -1] = pflRadiusNoHP  [models -1]/numFile;
				
				if(pmmdNative){
				    pflRmsdArray[models -1] = pflRmsdArray[models -1]/numFile;
				}
				
				for(k=0; k < pmmdThis->iResCount; k++){
					ppi4Surfacc[models -1][k] = ppi4Surfacc[models -1][k]/numFile;}
			}

			pflSurfaccTot = MemNew(sizeof(FloatLo) * iMinNumModels);
			/*fill the array containing the total surface accessibility*/
			for (models=0;models<iMinNumModels;models++){
				for(k=0; k < pmmdThis->iResCount; k++){
					pflSurfaccTot[models] += ppi4Surfacc[models][k];
				}
			}

			/*calculate the average secondary structures*/
			Calc2ndryStrucRepartAvg(ppfl2ndryStruct,ppc2ndrystruct,pppc2ndryStruct,  
									iMinNumModels, pmmdThis->iResCount, numFile);

			/*Energy graphs*/
            ComputeEnergyGraphs(pmsdArray, numFile,iMinNumModels, creaseIncl ,
                                creaseExcl, creaseDecay,fileName, pflRmsdArray,pflRadiusHP,pflRadiusNoHP);

			/*compute the contact map*/
			ComputeContactMap(pmsdArray, numFile, pmmdThis->iResCount, iMinNumModels, numPartition, fileName, pmmdNative);

			/*We only trace the NOES for the SH3 domain*/ 
			if(StringStr(fileName,"drk") != NULL){
				ComputeNOES(pmsdArray, numFile, iMinNumModels,numPartition,fileName);
			}


			if (output){
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_rGyr_upto",numFile,".dat");
				PrintRGyr(pflRadiusHP, pflRadiusNoHP ,iMinNumModels, pcTemp);

				if(pmmdNative){
				    sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_RMSD_upto",numFile,".dat");
				    Printfl(pflRmsdArray, iMinNumModels, "RMSD",pcTemp);
				}

				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_ASA_tot_upto",numFile,".dat");
				Printfl(pflSurfaccTot, iMinNumModels, "Acc. surf",pcTemp);
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_ASA_Plot_upto",numFile,".dat");
				Printppi4(ppi4Surfacc, iMinNumModels, pmmdThis->iResCount, "models", "residue", pcTemp);
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_2ndry_struct_percent_upto",numFile,".dat");
				Printfl3(ppfl2ndryStruct, pmsdArray[0]->iModels, pNames1,pcTemp);
				sprintf(pcTemp, "%s%s%d%s", fileName,"__Avrg_2ndrystruct_Plot_upto",numFile,".dat");
				PrintppChar(ppc2ndrystruct, pmsdArray[0]->iModels, pmmdThis->iResCount, "avg models", "residue",pcTemp);
			}

			if (createGif){

				flLow1	= findSmallestFloat(pflRadiusHP, iMinNumModels);
				if (flLow1 > findSmallestFloat(pflRadiusNoHP, iMinNumModels) ){
					flLow1 = findSmallestFloat(pflRadiusNoHP, iMinNumModels);
				}
				flHigh1	= findGreatestFloat(pflRadiusHP, iMinNumModels);
				if (flHigh1 <findGreatestFloat(pflRadiusNoHP, iMinNumModels) ){
					flHigh1 = findGreatestFloat(pflRadiusNoHP, iMinNumModels);
				}
				
				printf("\nRange of rGyr\n -greatest:%5.2f \n -smallest:%5.2f \n", flHigh1,flLow1);
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_rGyr_upto",numFile,".gif");	
				GraphArray(pcTemp, pflRadiusHP, pflRadiusNoHP, iMinNumModels, 
							1,flLow1,flHigh1, "HPonly", "All", "Rgyr", 
							"Time", 600,400, ARRAY_VS_TIME, NULL);	

				if(pmmdNative){
				    flLow1	= findSmallestFloat(pflRmsdArray, iMinNumModels);
				    flHigh1	= findGreatestFloat(pflRmsdArray, iMinNumModels);
				    sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_RMSD_upto",numFile,".gif");	
				    GraphArray(pcTemp, pflRmsdArray, NULL, iMinNumModels, 
						        1, flLow1,flHigh1, "RMSD", "", "RMSD", "Time",
						        600,400, ARRAY_VS_TIME, NULL);
                }

				flLow1	= findSmallestFloat(pflSurfaccTot, iMinNumModels);
				flHigh1	= findGreatestFloat(pflSurfaccTot, iMinNumModels);
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_ASA_tot_upto",numFile,".gif");	
				GraphArray(pcTemp,pflSurfaccTot, NULL, iMinNumModels, 
							1, flLow1,(FloatLo)flHigh1, "ASA", "", "ASA",
							"Time",600,400, ARRAY_VS_TIME, NULL);
				
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_ASA_Plot_upto",numFile,".gif");	
				GraphArray(pcTemp,NULL, NULL, iMinNumModels, 
					1, 0,pmmdThis->iResCount, "AREA OF SURFACE ACCESSIBILITY", "", "Residue Number", "Time",
					800,400, DOT_ASA_VS_TIME, ppi4Surfacc);	
				
				flLow1	= findSmallestFloatInppFloat(ppfl2ndryStruct, NUM_GROUP_2NDRY_STRUCT,pmsdArray[0]->iModels);
				flHigh1	= findGreatestFloatInppFloat(ppfl2ndryStruct, NUM_GROUP_2NDRY_STRUCT,pmsdArray[0]->iModels);
				sprintf(pcTemp, "%s%s%d%s", fileName,"_Avrg_2ndry_struct_percent_upto",numFile,".gif");	
				/*plot the progression of the various kinds of structures*/
				GraphArray(pcTemp,ppfl2ndryStruct[0], NULL, iMinNumModels, 
							1, 0.0,100.0, "helices", "", "Percentage of the residues", "Time",
							800,400, ARRAY_VS_TIME, NULL);						
				GraphArray(pcTemp,ppfl2ndryStruct[1], NULL, iMinNumModels, 
							2, 0.0,100.0, "sheets", "", "Percentage of the residues", "Time",
							800,400, ARRAY_VS_TIME, NULL);						
				GraphArray(pcTemp,ppfl2ndryStruct[2], NULL, iMinNumModels, 
							3, 0.0,100.0, "coils and others", "", "Percentage of the residues", "Time",
							800,400, ARRAY_VS_TIME, NULL);						
				/*graph the 2ndry structure averaged*/
				sprintf(pcTemp, "%s%s%d%s", fileName,"__Avrg_2ndrystruct_Plot_upto",numFile,".gif");	
				GraphArray(pcTemp,NULL, NULL, iMinNumModels, 
					1, 0,pmmdThis->iResCount, "Secondary structure assignment", "", "Residue Number", "Time",
					800,400, DOT_AVG_2NDRY_STUC_VS_TIME, ppc2ndrystruct);			
			}

			/*Free up the memory of all the arrays used*/
			MemFree(pmsdArray);

			MemFree(pflRadiusHP);
			MemFree(pflRadiusNoHP);
            
            if(pmmdNative){
			    MemFree(pflRmsdArray);
			}
			
			for (i=0; i <=NUM_MODEL; i++){
				MemFree(ppi4Surfacc[i]);
			}
			MemFree(ppi4Surfacc);
			MemFree(pflSurfaccTot);

			for(i=0; i<numFile;i++){
				for(j=0; j<=NUM_MODEL; j++){
					MemFree(pppc2ndryStruct[i][j]);
				}
				MemFree(pppc2ndryStruct[i]);
			}
			MemFree(pppc2ndryStruct);

			for(i=0;i<NUM_GROUP_2NDRY_STRUCT; i++){
				MemFree(ppfl2ndryStruct[i]);
			}
			MemFree(ppfl2ndryStruct);

			for(i=0; i<=NUM_MODEL; i++){
				MemFree(ppc2ndrystruct[i]);
			}
			MemFree	(ppc2ndrystruct);
			MemFree (iArrayNumModels);
		}

	    /*We trace only one graph*/
	    else if( numFile ==1 ){

			/* load an ASN.1 Biostruc */
			if (FileLength(fileName) == 0) 
                        {
				printf("\nInput movie structure file not found: %s", fileName);
				/* error occurred */
				return 3;
                        }
			printf("\nMovie structure loading ...filename: %s", fileName);
			pmsdRoot=LoadABiostruc(fileName,remoteFile,mLevel,modelNum);
			if (pmsdRoot==NULL) {
				printf("\nMovie structure failed to load with filename: %s", fileName);
				/* error occurred */
				return 3;
			}
			printf("\nMovie Structure loaded. %s\n", fileName);

			if (!IsStructureNode(pmsdRoot)){
				printf("\nThe pointer passed to RepresentRGyr must be a PMSD\n");
				return ERR_FAIL;
			}
			pdnmmHere = pmsdRoot->pdnmmHead;
			pmmdThis = (PMMD)pdnmmHere->data.ptrvalue;

			pflRadiusHP   = (FloatLo*) MemNew( pmsdRoot->iModels * sizeof(FloatLo));
			pflRadiusNoHP = (FloatLo*) MemNew( pmsdRoot->iModels * sizeof(FloatLo));

			if(pmmdNative){
			    pflRmsdArray  = (FloatLo*) MemNew( pmsdRoot->iModels * sizeof(FloatLo));
			}
			pflSurfaccTot = MemNew(sizeof(FloatLo) * pmsdRoot->iModels);
			ppi4Surfacc	  = (Int4 **)  MemNew( pmsdRoot->iModels * sizeof(Int4*) );  

			ppc2ndrystruct= (char **)  MemNew( pmsdRoot->iModels * sizeof(char*));
			ppfl2ndryStruct=(FloatLo**)MemNew(NUM_GROUP_2NDRY_STRUCT * sizeof(FloatLo*) );
			
			for(i=0; i<pmsdRoot->iModels; i++){
				ppc2ndrystruct[i] = (CharPtr) MemNew( (pmmdThis->iResCount * sizeof(char) ) +1 );
			}

			for(i=0; i<NUM_GROUP_2NDRY_STRUCT; i++){
				ppfl2ndryStruct[i] = (FloatLo*) MemNew( pmsdRoot->iModels * sizeof(FloatLo) );
			}

			for (models=1;models<=pmsdRoot->iModels;models++){
				/*retrieve the radius of gyration*/
				pflRadiusHP  [models -1] = GetExactRgyr(pmmdThis, models, TRUE);
				pflRadiusNoHP[models -1] = GetExactRgyr(pmmdThis, models, FALSE);

                if(pmmdNative){
				    /*retrieve the RMSD*/ 
				    Align2StrucSVD(	pmmdNative,(PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue),
									ALIGN_CA, 1, pmmdNative->iResCount, 1,models );
				    pflRmsdArray[models -1]= GetRMSD(pmmdNative, (PMMD)((pmsdRoot->pdnmmHead)->data.ptrvalue), 
									                    ALIGN_CA, 1, pmmdNative->iResCount, 1,models);
                }

				ppi4Surfacc[models -1] = CalcDSSPAccSurf((PMMD)pmmdThis, 2, models);	
				CalcDSSPAssignEx(pmmdThis, models , ppc2ndrystruct[models -1],FALSE);
			}

			/*fill the array containing the total surface accessibility*/
			for (models=0;models<pmsdRoot->iModels;models++){
				for(k=0; k < pmmdThis->iResCount; k++){
					pflSurfaccTot[models] += ppi4Surfacc[models][k];
				}
			}

			/*retrieve the 2ndry structure repartition*/
			Calc2ndryStrucRepart(ppfl2ndryStruct,ppc2ndrystruct,pmsdRoot->iModels,pmmdThis->iResCount);

			/*compute the contact map*/
			ComputeContactMap(&pmsdRoot, 1, pmmdThis->iResCount, pmsdRoot->iModels, numPartition, fileName, pmmdNative);

			/*Energy graphs*/
            ComputeEnergyGraphs(&pmsdRoot, 1,pmsdRoot->iModels, creaseIncl ,
                                creaseExcl, creaseDecay ,fileName,pflRmsdArray,pflRadiusHP,pflRadiusNoHP);
			
			/*We only trace the NOES for the SH3 domain*/ 
			if(StringStr(fileName,"drk") != NULL){
				ComputeNOES(&pmsdRoot, 1, pmsdRoot->iModels,numPartition,fileName);
			}		

			if (output){
				sprintf(pcTemp, "%s%s", fileName, "_rGyr.dat");	
				PrintRGyr(pflRadiusHP, pflRadiusNoHP ,pmsdRoot->iModels, pcTemp);

				if(pmmdNative){
				    sprintf(pcTemp, "%s%s", fileName, "_RMSD.dat");	
				    Printfl(pflRmsdArray,pmsdRoot->iModels,"RMSD",pcTemp);
				}

				sprintf(pcTemp, "%s%s%s", fileName,"_ASAtot",".dat");	
				Printfl(pflSurfaccTot, pmsdRoot->iModels, "Acc. surf",pcTemp);
				sprintf(pcTemp, "%s%s%s", fileName,"_ASA_Plot",".dat");	
				Printppi4(ppi4Surfacc, pmsdRoot->iModels, pmmdThis->iResCount, "models", "residue",pcTemp);
				sprintf(pcTemp, "%s%s%s", fileName,"_2ndrystruct_Plot",".dat");
				PrintppChar(ppc2ndrystruct, pmsdRoot->iModels, pmmdThis->iResCount, "models", "residue",pcTemp);
				sprintf(pcTemp, "%s%s%s", fileName,"_2ndry_struct_percent",".dat");
				Printfl3(ppfl2ndryStruct, pmsdRoot->iModels, pNames2,pcTemp);
				/*Printfl(ppfl2ndryStruct[3], pmsdRoot->iModels, "%age of non structure");*/
			}
			if (createGif){	
				flLow1	= findSmallestFloat(pflRadiusHP, pmsdRoot->iModels);
				if (flLow1 > findSmallestFloat(pflRadiusNoHP, pmsdRoot->iModels) ){
					flLow1 = findSmallestFloat(pflRadiusNoHP, pmsdRoot->iModels);
				}
				flHigh1	= findGreatestFloat(pflRadiusHP, pmsdRoot->iModels);
				if (flHigh1 <findGreatestFloat(pflRadiusNoHP,  pmsdRoot->iModels) ){
					flHigh1 = findGreatestFloat(pflRadiusNoHP, pmsdRoot->iModels);
				}
				printf("\nRange of rGyr\n -greatest:%5.2f \n -smallest:%5.2f \n", flHigh1,flLow1);
				sprintf(pcTemp, "%s%s", fileName, "_rGyr.gif");	
				GraphArray(pcTemp, pflRadiusHP, pflRadiusNoHP, pmsdRoot->iModels, 
							1, flLow1,flHigh1, "HPonly", "All", "Rgyr",
							"Time",600,400, ARRAY_VS_TIME, NULL);		


                if(pmmdNative){
				    flLow2	= findSmallestFloat(pflRmsdArray, pmsdRoot->iModels);
				    flHigh2	= findGreatestFloat(pflRmsdArray, pmsdRoot->iModels);
				    printf("\nRange of RMSD\n -greatest:%5.2f \n -smallest:%5.2f \n", flHigh2,flLow2);
				    sprintf(pcTemp, "%s%s", fileName, "_RMSD.gif");	
				    GraphArray(pcTemp, pflRmsdArray, NULL, pmsdRoot->iModels, 
						        1, flLow2,flHigh2, "RMSD", "", "RMSD", "Time",
						        600,400, ARRAY_VS_TIME, NULL);
                }

				flLow1	= findSmallestFloat(pflSurfaccTot, pmsdRoot->iModels);
				flHigh1	= findGreatestFloat(pflSurfaccTot, pmsdRoot->iModels);
				sprintf(pcTemp, "%s%s%s", fileName,"_ASAtot",".gif");	
				GraphArray(pcTemp,pflSurfaccTot, NULL, pmsdRoot->iModels, 
					1, flLow1,(FloatLo)flHigh1, "ASA", "", "ASA", "Time",
					600,400, ARRAY_VS_TIME, NULL);						

				sprintf(pcTemp, "%s%s%s", fileName,"_ASA_Plot",".gif");	
				GraphArray(pcTemp,NULL, NULL, pmsdRoot->iModels, 
					1, 0,pmmdThis->iResCount, "AREA OF SURFACE ACCESSIBILITY", "", "Residue Number", "Time",
					800,400, DOT_ASA_VS_TIME, ppi4Surfacc);						

				sprintf(pcTemp, "%s%s%s", fileName,"_2ndrystruct_Plot",".gif");
				GraphArray(pcTemp,NULL, NULL, pmsdRoot->iModels, 
					1, 0,pmmdThis->iResCount, "Secondary structure assignment", "", "Residue Number", "Time",
					800,400, DOT_2NDRY_STUC_VS_TIME, ppc2ndrystruct);			

				flLow1	= findSmallestFloatInppFloat(ppfl2ndryStruct, NUM_GROUP_2NDRY_STRUCT,pmsdRoot->iModels);
				flHigh1	= findGreatestFloatInppFloat(ppfl2ndryStruct, NUM_GROUP_2NDRY_STRUCT,pmsdRoot->iModels);
				sprintf(pcTemp, "%s%s%s", fileName,"_2ndry_struct_percent_",".gif");	
				GraphArray(pcTemp,ppfl2ndryStruct[0], NULL, pmsdRoot->iModels, 
					1, 0.0,100.0, "helices", "", "Percentage of the residues", "Time",
					800,400, ARRAY_VS_TIME, ppc2ndrystruct);						
				GraphArray(pcTemp,ppfl2ndryStruct[1], NULL, pmsdRoot->iModels, 
					2, 0.0,100.0, "sheets", "", "Percentage of the residues", "Time",
					800,400, ARRAY_VS_TIME, ppc2ndrystruct);						
				GraphArray(pcTemp,ppfl2ndryStruct[2], NULL, pmsdRoot->iModels, 
					3, 0.0,100.0, "Coils and others", "", "Percentage of the residues", "Time",
					800,400, ARRAY_VS_TIME, ppc2ndrystruct);						
				/*GraphArray(pcTemp,ppfl2ndryStruct[3], NULL, pmsdRoot->iModels, 
					4, flLow1,flHigh1, "not in a structure", "", "Percentage of the residues", "Time",
					800,400, ARRAY_VS_TIME, ppc2ndrystruct);						*/

			}

			/*Memory freeing*/
			MemFree(pflRadiusHP);
			MemFree(pflRadiusNoHP);
			if(pmmdNative){
			    MemFree(pflRmsdArray);
			}
			MemFree(pflSurfaccTot);

			for (i=0; i <pmsdRoot->iModels; i++){
				MemFree(ppi4Surfacc[i]);
			}
			MemFree(ppi4Surfacc);

			for(i=0; i<pmsdRoot->iModels; i++){
				MemFree(ppc2ndrystruct[i]);
			}
			MemFree(ppc2ndrystruct);

			for(i=0; i<NUM_GROUP_2NDRY_STRUCT; i++){
				MemFree(ppfl2ndryStruct[i]);
			}
			MemFree(ppfl2ndryStruct);
		}

	    /*print all the curves of the various structures on the same gif picture*/
	    else {
            pmsdArray=(PMSD*) MemNew(numFile * sizeof(PMSD));

			/*Allocate memory for the two dimensionnal array*/
			pppflRadiusHP   = (FloatLo **) MemNew (sizeof (FloatLo*) * numFile);  
			pppflRadiusNoHP = (FloatLo **) MemNew (sizeof (FloatLo*) * numFile);

			for (i=0; i< numFile; i++){
				pppflRadiusHP[i]   = (FloatLo*)MemNew( sizeof(FloatLo)* (NUM_MODEL+1) );   
				pppflRadiusNoHP[i] = (FloatLo*)MemNew( sizeof(FloatLo)* (NUM_MODEL+1) ); 
			}

			/*allocate memory for the rmsd array*/

			if(pmmdNative){
			    pppflRmsdArray = (FloatLo **) MemNew (sizeof (FloatLo*) * numFile);
			    for (i=0; i< numFile; i++){
				    pppflRmsdArray[i] = (FloatLo*)MemNew( sizeof(FloatLo)* (NUM_MODEL+1) );
			    }
			}

			for(i=1; i <= numFile; i++)
			{
				/* load an ASN.1 Biostruc */
				sprintf(pcTemp, "%s%s%d%s",  fileName,"_movie_",i, ".val");

				if (FileLength(pcTemp) == 0) 
	                        {
					printf("\nInput sequential movie structure file not found: %s", pcTemp);	
					/* error occurred */
					return 3;
       		                 }
				printf("\nStructure loading ...filename: %s", pcTemp);

				pmsdRoot=LoadABiostruc(pcTemp,remoteFile,mLevel,modelNum);
				if (pmsdRoot==NULL) {
					printf("\nError loading sequential movie structure file: %s", pcTemp);
					/* error occurred */ 
					return 3;
				}
				printf("\nSequential Movie Structure # %d loaded.\n", i);
				pmsdArray[i-1] = pmsdRoot;
				pdnmmHere = pmsdRoot->pdnmmHead;
				pmmdThis = (PMMD)pdnmmHere->data.ptrvalue;

				for (models=1;models<=pmsdRoot->iModels;models++)
				{
					pppflRadiusHP  [i-1 ][models -1] = GetExactRgyr(pmmdThis, models, TRUE);
					pppflRadiusNoHP[i-1 ][models -1] = GetExactRgyr(pmmdThis, models, FALSE);
				}
				

				if(pmmdNative){
				    for (models=1;models<=pmsdRoot->iModels;models++)
				    {
					    Align2StrucSVD(	pmmdNative,(PMMD)((pmsdArray[i-1]->pdnmmHead)->data.ptrvalue),
									    ALIGN_CA, 1, pmmdNative->iResCount, 1,models );
					    pppflRmsdArray [i-1 ][models -1] = GetRMSD(pmmdNative, (PMMD)((pmsdArray[i-1]->pdnmmHead)->data.ptrvalue), 
																    ALIGN_CA, 1, pmmdNative->iResCount, 1,models);
				    }
                }
			}

			/*retrieve the range of radius of gyration*/
			flLow1 = findSmallestFloat ( pppflRadiusHP[0], pmsdArray[0]->iModels);
			flHigh1= findGreatestFloat ( pppflRadiusHP[0], pmsdArray[0]->iModels);
			for(i=0; i<numFile; i++){
				if(flLow1 > findSmallestFloat ( pppflRadiusHP[i], pmsdArray[i]->iModels) ){
					flLow1 = findSmallestFloat ( pppflRadiusHP[i], pmsdArray[i]->iModels);
				}
				if(flHigh1< findGreatestFloat ( pppflRadiusHP[i], pmsdArray[i]->iModels) ){
					flHigh1 = findGreatestFloat ( pppflRadiusHP[i], pmsdArray[i]->iModels);
				}
				if(flLow1 > findSmallestFloat ( pppflRadiusNoHP[i], pmsdArray[i]->iModels) ){
					flLow1 = findSmallestFloat ( pppflRadiusNoHP[i], pmsdArray[i]->iModels);
				}
				if(flHigh1< findGreatestFloat ( pppflRadiusNoHP[i], pmsdArray[i]->iModels) ){
					flHigh1 = findGreatestFloat ( pppflRadiusNoHP[i], pmsdArray[i]->iModels);
				}
			}
			printf("rGyr:\n -lowest: %5.2f\n -highest: %5.2f\n", flLow1,flHigh1);

			if(pmmdNative){
			    /*retrieve the range of rmsd*/
			    flLow2 = findSmallestFloat ( pppflRmsdArray[0], pmsdArray[0]->iModels);
			    flHigh2= findGreatestFloat ( pppflRmsdArray[0], pmsdArray[0]->iModels);
			    for(i=0; i<numFile; i++){
				    if(flLow2 > findSmallestFloat ( pppflRmsdArray[i], pmsdArray[i]->iModels) ){
					    flLow2 = findSmallestFloat ( pppflRmsdArray[i], pmsdArray[i]->iModels);
				    }
				    if(flHigh2< findGreatestFloat ( pppflRmsdArray[i], pmsdArray[i]->iModels) ){
					    flHigh2 = findGreatestFloat ( pppflRmsdArray[i], pmsdArray[i]->iModels);
				    }
			    }
			    printf("RMSD:\n -lowest: %5.2f\n -highest: %5.2f\n", flLow2,flHigh2);
			}

			for(i=1; i <= numFile; i++)
			{
				if (output){
					sprintf(pcTemp, "%s%s%d%s", fileName,"_Upto",numFile, "_rGyr.dat");	
					PrintRGyr(pppflRadiusHP[i-1], pppflRadiusNoHP[i-1] ,pmsdArray[i-1 ]->iModels, pcTemp);

					if(pmmdNative){
					    sprintf(pcTemp, "%s%s%d%s", fileName,"_Upto",numFile, "_RMSD.dat");
					    Printfl(pppflRmsdArray[i-1], pmsdArray[i-1 ]->iModels,"RMSD",pcTemp);
                    }

				}
				if (createGif){		
					sprintf(pcTemp, "%s%s%d%s", fileName,"_Upto",numFile, "_rGyr.gif");	
					GraphArray(pcTemp, pppflRadiusHP[i-1], pppflRadiusNoHP[i-1], 
							pmsdArray[i -1]->iModels, i, flLow1, flHigh1, "HPonly", 
							"All", "Rgyr", "Time",
							1200,800, ARRAY_VS_TIME, NULL);		

					if(pmmdNative){
					    sprintf(pcTemp, "%s%s%d%s", fileName,"_Upto",numFile, "_RMSD.gif");	
					    GraphArray(pcTemp, pppflRmsdArray[i-1], NULL, 
								    pmsdArray[i -1]->iModels, i, flLow2, flHigh2, "RMSD",
								    "", "RMSD", "Time",1200,800, ARRAY_VS_TIME, NULL);
                    }
				}
			}

			/*Memory freeing*/
			MemFree(pmsdArray);

			/* free the 2 dimensionnal array*/
			for (i=0; i <numFile; i++){
				MemFree(pppflRadiusHP[i]);
				MemFree(pppflRadiusNoHP[i]);
			}
			MemFree(pppflRadiusHP);
			MemFree(pppflRadiusNoHP);

			if(pmmdNative){
			    for (i=0; i <numFile; i++){
				    MemFree(pppflRmsdArray[i]);
			    }
			    MemFree(pppflRmsdArray);
			}

	}
	return 0;
}


Int2 Main()
{
	Int2 ModelNum;

        ErrSetLogfile("error_analyzeMovie.log", ELOG_APPEND|ELOG_BANNER);
	ErrSetOptFlags(EO_SHOW_SEVERITY|EO_SHOW_CODES|EO_LOG_FILELINE|EO_SHOW_USERSTR|EO_SHOW_ERRTEXT|EO_BEEP|EO_WAIT_KEY|EO_LOGTO_USRFILE);
	if (!GetArgs("Analyze Movie Program",NUMARGS,Rargs))
		return 1;
	if (!OpenMMDBAPI(0,NULL)) {
		ErrPostEx(SEV_FATAL,2,1,"Unable to open MMDBAPI");
		return 2;
	}

	ModelNum=1;
	RetrieveInfo(	Rargs[0].strvalue,FALSE,ALLMDL,&ModelNum,
					TRUE, TRUE, Rargs[2].intvalue,
					Rargs[3].intvalue, Rargs[4].intvalue, 700, 520,
					0, 3, FALSE, Rargs[1].strvalue);

	/* Free the Modelstruc (and its enclosed Biostruc) */
	/* FreeAModelstruc(PDNMS pdnmsThis); not necessary */
	/* This can be done individually - but all Modelstrucs */
   	/* remaining are freed in CloseMMDB-API() */
	/* free the atom list created */
	/* Shut Down MMDB-API */
	printf("\nNow closing the structures.\n\n");
	CloseMMDBAPI();	
 	return TRUE;
}


/*  
$Log: analyzeMovie.c,v $
Revision 1.5  2004/07/08 19:05:01  egarderm
Changed code and added explicit conditions to allow code to run without a native structure. Also changed the argument for native structure to be optional

Revision 1.4  2003/09/14 02:25:49  feldman
Removed unused variables

Revision 1.3  2003/04/06 20:02:54  feldman
Removed SS objects side effect

Revision 1.2  2003/03/21 23:26:07  feldman
Checked in wrong version before

Revision 1.33  2003/03/20 21:59:45  lewis
Tidied up graphs and naming a bit more

Revision 1.32  2003/03/19 16:45:48  feldman
Removed extraneous options

Revision 1.31  2003/03/14 20:57:03  feldman
Updated NOEs to correct ones for drk

Revision 1.30  2003/03/12 20:33:29  lewis
Fixed x-axis scaling for when num frames is not 100

Revision 1.29  2003/03/11 21:42:55  lewis
Made detection of native structure filename much smarter and also changed movie length to MAXGEN

Revision 1.28  2003/03/11 19:35:16  lewis
Native contact map is now of actual native structure, and excluded from partition 1 map

Revision 1.27  2003/03/11 16:49:27  feldman
Fixed native contact map to actual use native structure file

Revision 1.26  2002/11/13 04:13:44  feldman
Minor bug fix

Revision 1.25  2002/11/12 21:03:58  feldman
One more minor fix

Revision 1.24  2002/11/12 20:53:50  feldman
Tidied up graphs and fixed a number of tricky bugs

Revision 1.23  2001/10/23 15:27:28  feldman
Fixed several minor bugs, plots should be mostly correct now...

Revision 1.22  2001/10/16 20:10:37  feldman
Fixed a bug with model numbering

Revision 1.21  2001/10/11 22:05:24  feldman
Fixed rainbow coloring problem due to round-off error

Revision 1.20  2001/10/11 18:48:30  feldman
Fixed a few minor bugs, should work with all movies now

Revision 1.19  2001/08/31 19:26:47  feldman
Removed unused variables

Revision 1.18  2001/08/30 19:02:39  phan
Added NOE's graphs.

Revision 1.17  2001/08/20 15:57:13  feldman
removed unused variables and libraries

Revision 1.16  2001/08/13 19:26:56  phan
Added the contact order per residue

Revision 1.15  2001/08/13 01:10:57  phan
Added RGYR vs RMSD

Revision 1.14  2001/08/03 17:47:44  feldman
Fixed a few bugs

Revision 1.13  2001/08/03 17:37:47  phan
Added Energy calculations

Revision 1.12  2001/08/01 17:48:07  phan
Corrected compiling warnings due to mismatch in printf statements

Revision 1.11  2001/08/01 17:13:51  feldman
Fixed compiler warnings

Revision 1.10  2001/08/01 16:49:59  phan
Added contact maps.

Revision 1.8  2001/07/19 15:19:33  phan
Corrected good scaling when drawing rGyr.

Revision 1.7  2001/07/13 16:00:47  phan
Added creation of tab seperated files with graph datas

Revision 1.6  2001/07/12 14:37:45  phan
Added plotting of secondery structure when averaging movies.

Revision 1.5  2001/07/05 17:08:27  phan
Update the represent-it function that now can calculate the average for ASA

Revision 1.4  2001/06/28 14:49:15  phan
Fixed the vertical scale in the "graphArray function"

Revision 1.3  2001/06/27 20:56:40  phan
Reviewed allocation memory(use of MemNew instead of Malloc)
Optimized function to plot graph.

Revision 1.2  2001/06/25 15:34:28  phan
Fixed a problem with the scale on the gif pictures.

Revision 1.1  2001/06/22 19:17:19  phan
Program to trace diagrams of rGyr, RMSD .... of giiven MMDB structures.

*/

