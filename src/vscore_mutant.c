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
/**********************************************************************
* PROGRAM Vscore-dec2003.c
* calculates the contact score for a specified protein, based
* on a Voronoi tesselation procedure.
*
**********************************************************************/

#include "vscoredata.h"
#include <mmdbtraj.h>


#define CELLSIZE 7.0  /* set to (maximum contact radius x 2) */
#define Rw 1.4        /* radius of water */
#define PAUSE ch=getchar()  /* for debugging */

#define ATOMNAME_LEN 5
#define RESNAME_LEN 4
#define MAX_ATOM_CLASS 168

#define MAX_UNKNOWN_TOLERANCE 30  /* increased to 30 CWVH 2012 */

/* ------------------------- structure definitions ---------------------- */

struct plane {
    double Ai[4];      /* parameters A,B,C,D of contact plane Ax+By+Cz+D=0  */
    double dist;       /* distance from plane to origin */
    int    index;      /* index to which record in PDB or ligand array. */
    double area;       /* contact area in square angstroms */
    char   flag;       /* 'X' if no contact, 'E' if an engulfing atom. */
};

struct atom {
    int   atomnum;      /* record number from PDB */
    float coor[3];      /* xyz coordinates of atom */
    char  atomname[ATOMNAME_LEN];  /* name of atom, CA, etc */
    char  res[RESNAME_LEN];       /* residue name from PDB */
    int   resnum;
    float radius;       /* radius of atom */
    int   boxnum;       /* the box number atom is assigned to. */
    char  source;       /* source of atom, protein, ligand,  */
    float SAS;          /* solvent exposed surface area */
    float vol;          /* atom volume */
    char  done;         /* flag if atom contacts have already been calculated */
    int   class;        /* atom class, 0 to 167. */
    float score;        /* atom score */
    float solvscore;    /* solvent score */
/*	PMAD  pmadAtom;      pointer back to MMDBAPI atom */
};

struct atomindex {
    int   nument;      /* number of entries in box  */
    int   first;       /* location of first entry in PDBlist */
};

struct contactlist {
    int     index;   /* index to PDB atom */
    double  area;    /* contact area, square angstroms */
    double  dist;    /* distance to atom zero */
    char    flag;    /* to keep or not. 'X' == omit. */
};

struct contactinfo {
    int  atomcount;  /* number of contributing atoms */
    char  res[RESNAME_LEN];    /* residue name */
    char name[ATOMNAME_LEN];    /* atom name */
    float areasum;   /* sum of areas for contributing atoms (for weighted average) */
    int  ctctot;     /* total number of recorded contacts */
};

struct site {
    double xi[3];
    int atom[20];
    int maxnum;
};


/* -------------- structure definitions for subroutines ----------------- */

struct vertex {
    double xi[3];      /* x,y,z coordinates (x1,x2,x3) */
    double dist;       /* distance to origin */
    int    plane[3];   /* identification of intersecting planes. -1 = sphere. */
};
struct ptindex {
    int  numpts;       /* number of points defining face */
    /*int  pt[40];*/     /*   index to polyhedron points */
    int  pt[50];       /* index to polyhedron points */
};
struct edgevector {
    double V[3];       /* vector for edge */
    int startpt;       /* initial vertex */
    int endpt;         /* final vertex */
    int plane[2];      /* planes defining edge (using single atoms for now) */
    int startplane;    /* third plane at start point */
    int endplane;      /* third plane at end point */
    char arc;          /* flag for arc point calculations */
};
struct ca_struct {
    long int prev;     /* previous contact location in ca_index */
    long int atom;     /* PDBarray number of current contact (NOT PDB record number) */
    float area;        /* contact area */
    float dist;        /* distance between atoms */
};


/* ----------------- Global variables ----------------- */

struct atom *PDB;                 /* pointer to PDB array (dynamically allocated) */
struct atomindex *box;            /* index to PDB atoms within cubic grid */
int      *PDBlist;                /* list of atoms ordered by box number */
long int *seed;                   /* seed vertices for new polyhedr */
char      planedef;               /* = X, R, or B */
struct ca_struct *ca_rec;         /* array - contact area records */
long int *ca_index;               /* array - index to first ca_recs for each atom. */
long int  numcarec = 0;           /* number of contact atom records   */
long int  ca_recsize;             /* current max size of array */
int       dim;                    /* size of side of box, in units CELLSIZE. */
float     globalmin[3];           /* minimum coordinate value for protein atoms */
float     globalmax[3];           /* maximum coordinate value for protein atoms */
struct plane       cont[200];     /* atom and contact plane information */
struct contactlist contlist[200]; /* list of possible atom contacts */
char      ch;                     /* for debugging  */
char      showbonded;             /* = Y or N. */
char      normalize;              /* = Y or N. normalize areas to area of sphere. */
struct contactinfo continfo[170]; /* contact information for elements in array continfo */
float contfreq[170][170];         /* array of contact frequencies */

int       numsites = 0;
struct site ionsite[100];

double    totscore;    /* atom-atom contact score */
double    solvscore;   /* atom-solvent contact score */
double    SSscore;     /* solvent-solvent contact score */
char     FILENAME[60];    /* input file name */

/* Static global variable - used internally */
static FloatHi TotalVPotential;
static FloatHi TotScore;
static FloatHi SolvScore;
static FloatHi SSScore;
static FloatHi ASATotal; 
static Int4 NumPDBatoms;
/* CWVH 2012 testing 
static ValNodePtr pvnHead = NULL;
static ValNodePtr pvnCur = NULL;
static FILE *fp;
static Int4 CapCount = 0;
*/

/* ----------------- Global variables for subroutines ----------------- */

struct vertex     poly[400];     /* polyhedron vertices */
struct edgevector vedge[400];
struct ptindex    ptorder[200];  /* for ordering vertices around each face */
struct vertex     centerpt[200]; /* center points for each contact face */


/* --------------------- function prototypes ---------------------- */

int    calc_freqarray(struct atom PDB[], int PDBtot, int atomlist[], int numatoms);
void   init_contfreq(struct contactinfo continfo[], float *contfreq);
TrajErr assign_class(struct atom PDB[], long int PDBtot);
struct atom *read_PDB2(char *FILE_ptr, struct atom *PDB, long int *numPDBatoms, char HETATMs);
int    index_protein(int totatoms);
double calc_score(struct atom PDB[], long int PDBtot, int atomlist[],int windowSize);
void   ion_sites(struct plane cont[], struct vertex poly[], int NV, int atomzero);

/* ------- prototypes from subroutines -------------- */
double spherical_arc(struct vertex ptAo, struct vertex ptB, struct vertex ptC, float rado);
double cosPQR( double ptP[], double ptQ[], double ptR[]);
char   test_point(double ptX[], struct plane cont[], int NC, float rado, int planeA, int planeB, int planeC);
char   order_faces(int atomzero, struct vertex poly[], struct vertex centerpt[], float rado, int NC, int NV, struct plane cont[], struct ptindex ptorder[]);
void   project_points(struct vertex poly[], struct vertex centerpt[], float rado, int NC, int NV, struct plane cont[]);
int    voronoi_poly2(struct atom PDB[], int atomzero, struct plane cont[], float rado, int NC, struct contactlist contlist[]);
int    add_vertex(struct vertex poly[], int vn, double coor[], int planeA, int planeB, int planeC);
void   add_vedge(struct edgevector vedge[], int edgenum, struct plane cont[], int plane0, int plane1, int testplane, struct vertex poly[], int startpt);
int    solve_3x3(double eq0[], double eq1[], double eq2[], double pt[]);
int    solve_2xS(struct plane eq0, struct plane eq1, float rado, double pt0[], double pt1[]);
void   save_seeds(struct plane cont[], struct vertex poly[], int NV, int atomzero);
void   get_firstvert(struct plane cont[], int *planeA, int *planeB, int *planeC, int NC, int atomzero);
void   calc_areas(struct vertex poly[], struct vertex centerpt[], float rado, int NC, int NV, struct plane cont[], struct ptindex ptorder[], int atomzero);
int    index_protein(int PDBtot);
void   assign_radii(struct atom *PDB, long int PDBtot);
void   save_areas(struct plane cont[], struct contactlist contlist[], int NC, int atomzero);
int    get_contlist4(struct atom PDB[], int atomzero, struct contactlist contlist[], int PDBtot, float rado, int dim);

/* ======================================================================== */

/**
 * Function parsePMMD
 * A function to parse a pmmd structure into the data strcutures required
 * for further score calculation
 * The function populates the global PDB array with parsed data
 * Note: hydrogen atoms are skipped during parsing
 * Note: parser reads only one of multiple possible structures - if the AltLoc
 * record is different than ' ' or 'A', atom is omitted
 * @param molecule - the PMMD structure to be parsed
 * @param numAtoms - the number of atoms in the structure, get filled in when
 * the function completes
 * @return an array of atom ranges for each residue
 */

Int4Ptr* parsePMMD(PMMD molecule, long int *numAtoms)
{
    long int PDBarraysize = 0;
    long int atomCount = 0; /* number of atoms read from the array */
    Int4 lastAtomCount = atomCount+1;
    Int2 resIndex = 0;
    Int2 first = 0;
    Int2 last = 1;
    Int2 i;
    
    Int4 numRes;
    PDNMG pdnmgHeadRes;
    PMGD pmgdResX;
    PVNMA pvnmaResAtoms;
    PMAD pmadCurAtom;
    PALD paldCoords;
    
    Int4Ptr* resBreakdown;
    
    /* Traverse molecule to get the atom count */
    pdnmgHeadRes = molecule->pdnmgHead;
    while(pdnmgHeadRes!=NULL){
        pmgdResX = (PMGD)pdnmgHeadRes->data.ptrvalue;       
        PDBarraysize += pmgdResX->iAtomCount;
        pdnmgHeadRes = pdnmgHeadRes->next;
    }
        
    /* Initialize array to total number of atoms - size may be larger than # atoms w/known coords */    
    PDB = (struct atom *)(MemNew(PDBarraysize*sizeof(struct atom)));
    if(!PDB) {
        ErrPostEx(SEV_ERROR,1,1,"Memory allocation for atom array failed while calculating v-score\n");
        return NULL;
    }
    
    numRes = molecule->iResCount;
    /* Initialize 2D array to store atom index boundaries (only with coords) for each residue */
    resBreakdown = (Int4Ptr *)MemNew((numRes+1) * sizeof(Int4Ptr));
    for(i = 0; i <= numRes; i++){
        resBreakdown[i] = (Int4Ptr) MemNew(2 * sizeof(Int4Ptr));
        resBreakdown[i][first]=0;
        resBreakdown[i][last]=0;
    }
    
    pdnmgHeadRes = molecule->pdnmgHead;
    /* Iterate through each residue in the molecule */
    while(pdnmgHeadRes!=NULL){ 

        /* Extract current residue node */
        pmgdResX = (PMGD)(pdnmgHeadRes->data.ptrvalue);
        
        /* ------- Skip HOH records ---------- */
        if (!StringCmp(pmgdResX->pcGraphName,"HOH")){
	        pdnmgHeadRes = pdnmgHeadRes->next;
            continue;
        }
        
        /* Should skip alternate residue conformation - taken care of by biostruc */
        
        /* Now get the atom list for this residue */
        pvnmaResAtoms = pmgdResX->pvnmaAHead;
        
        /* Iterate through the atoms */
        while (pvnmaResAtoms != NULL){
            
            pmadCurAtom = (PMAD)(pvnmaResAtoms->data.ptrvalue);
            
            /* If atom has no coordinates, it is not useful for the score calculation - ignore */
            paldCoords=GetAtomLocs(pmadCurAtom,1/* model num*/);
            if (paldCoords==NULL) {
                pvnmaResAtoms = pvnmaResAtoms->next;
                continue;
            }
        
            /* ------- Skip alternate location records ---------- */
            if ((paldCoords->cAltConf != ' ') && (paldCoords->cAltConf != 'A')){
                pvnmaResAtoms = pvnmaResAtoms->next;
                continue;
            }
                        
            /* ------- Skip explicit hydrogens ------- */
            if (pvnmaResAtoms->choice == 1){
                pvnmaResAtoms = pvnmaResAtoms->next;
                continue;
            }

            /* Now that all conditions were checked, store the relevant info in the respective atom struc */
            
            /* atom number */
            PDB[atomCount].atomnum = atomCount;
            
            /* atom name */
            StringNCpy(PDB[atomCount].atomname, pmadCurAtom->pcAName, ATOMNAME_LEN-1);
            PDB[atomCount].atomname[ATOMNAME_LEN-1] = '\0';
            
            /* residue name */
            StringNCpy(PDB[atomCount].res, pmgdResX->pcGraphName, RESNAME_LEN-1);
            PDB[atomCount].res[RESNAME_LEN-1] = '\0';
            
            /* residue number */
            sscanf(pmgdResX->pcGraphNum,"%d",&(PDB[atomCount].resnum));
           
            /* Inititalize default values */
            PDB[atomCount].done = 'N';
            PDB[atomCount].score = 0.0;
            PDB[atomCount].solvscore = 0.0;
        
            /* Now extract the actual coordinates */
            PDB[atomCount].coor[0] = AtomLocX(paldCoords);   /* x coordinate */
            PDB[atomCount].coor[1] = AtomLocY(paldCoords);   /* y coordinate */
            PDB[atomCount].coor[2] = AtomLocZ(paldCoords);   /* z coordinate */

		/*	PDB[atomCount].pmadAtom = pmadCurAtom;  pointer back to MMDBAPI data structure CWVH 2012 */ 
            
            atomCount++;
            
            pvnmaResAtoms = pvnmaResAtoms->next;
        }
        resBreakdown[resIndex][first] = lastAtomCount;
        resBreakdown[resIndex][last] = atomCount;
        lastAtomCount = atomCount+1;
        resIndex++;
        
        pdnmgHeadRes = pdnmgHeadRes->next;
    }
    *numAtoms = atomCount;
    assign_radii(PDB, atomCount);
    return resBreakdown;
}


TrajErr ComputeVscorePotential (PMMD pmmdToScore, DValNodePtr *ppdnResult, Int2 iWindowSize)
{
    int      atomi;           /* atom number of atomzero */
    long int PDBtot;          /* number of PDB atoms */

    /* ion site variables. */
    int sitei;    /* site counter. */
    int sai;      /* site contact atom counter (0-3) */
    int siteatom; /* index to site atom. */
    int cnt;
    
    BiostrucPtr bspBiostruc = NULL;
    
    PMMD pmmdMol;
 	PDNMG pdnmgHere;
 	PMGD pmgdHere;
	Int4 totNumAtoms = 0;
	Int4Ptr * atomsPerRes;
	Int4 i, j;
	
	/* Variable for list storing results */
	AdjListNodePtr paln;
	DValNodePtr pdnHeadSublist;
	int resNum;
	
	showbonded = 'N';  /* default */
    planedef = 'R';    /* new default */
    normalize = 'N';   /* default */
    
    TotalVPotential = 0.0;
	
	pmmdMol = pmmdToScore;
	
	if ( !IsProtein((PFB)pmmdMol) ){
        ErrPostEx(SEV_ERROR,4,1,"Structure is not a protein, aborting\n");
        return ERR_FAIL;
    }
        
    numsites = 0;
    numcarec = 0;
    totNumAtoms = 0;
    
	atomsPerRes = parsePMMD(pmmdMol, &PDBtot);
	if (!atomsPerRes){
	    ErrPostEx(SEV_ERROR,1,1,"Structure parsing failed\n");
	    return ERR_FAIL;
	}
	if(PDBtot <= 0){
	    ErrPostEx(SEV_ERROR,1,1,"Could not establish atom count for structure\n");
	    return ERR_FAIL;
	}
	else
	{
		NumPDBatoms = (Int4) PDBtot;
	}

    /* Initialize contact atom index */
    ca_index = (long int *) MemNew(PDBtot * sizeof(long int));
    ca_recsize = 5*PDBtot;
    ca_rec = (struct ca_struct *) MemNew(ca_recsize*sizeof(struct ca_struct));
    seed = (long int *) MemNew(3*PDBtot*sizeof(long int));

    if((!ca_rec) || (!ca_index) || (!seed)) {
        ErrPostEx(SEV_FATAL,1,1,"Memory allocation error while calculating v-score");
    }
    for(atomi=0; atomi<PDBtot; ++atomi) {
        PDB[atomi].vol = 0.0;
        ca_index[atomi] = -1;   /*initialize pointer array */
        seed[atomi*3] = -1;     /* initialize seed array */
    }
        
    /* Initialize global min and max arrays */
    for (i=0;i<3;i++){
        globalmin[i] = 0.0;
        globalmax[i] = 0.0;
    }

    /* assign protein atoms to boxes in cubic grid */
    dim = index_protein(PDBtot);
    /* initialize contact array and assign atom types */
    init_contfreq(continfo, contfreq[0]);
    if (assign_class(PDB, PDBtot)==ERR_FAIL){
        /* Abort */
        ErrPostEx(SEV_ERROR,1,1,"Could not assign structure to class\n");
        return ERR_FAIL;
    }

    /* calculate protein score */
    /* last argument is residue adjacency window size to be ignored in score */
    totscore = calc_score(PDB, PDBtot, PDBlist,iWindowSize);

    /* Check for ion sites - they would otherwise affect the score */
    if(numsites == 0) {
        /* No action required */
    } else {
        for(sitei=0; sitei<numsites; ++sitei) {
            for(sai=0; sai<ionsite[sitei].maxnum; ++sai) {
                siteatom = ionsite[sitei].atom[sai];
                if(PDB[siteatom].solvscore < 0.0) {
                    PDB[siteatom].solvscore = 0.0;
                }
                if(PDB[siteatom].score < 0.0) {
                    PDB[siteatom].score = 0.0;
                }
            }
        }
    }
    
    /* Traverse the residues and extract the scores */
    pdnmgHere = pmmdMol->pdnmgHead;
    for (i=0; i<pmmdMol->iResCount; i++){
        totscore = 0.0;
        solvscore = 0.0;
        /* Extract the residue struct */
        pmgdHere = (PMGD)(pdnmgHere->data.ptrvalue);
        
        /* Calculate the score for the current residue */
        for (j=atomsPerRes[i][0]; j<=atomsPerRes[i][1]; j++){
            if(PDB[j-1].class < MAX_ATOM_CLASS) {
                totscore += PDB[j-1].score;
                solvscore += PDB[j-1].solvscore;
                totNumAtoms++;
            }
        }
        
        /* Now build list component from the bottom up */
    	
		/* Add the residue to the head of the sublist */
		paln = (AdjListNodePtr) MemNew (sizeof (AdjListNode));
		paln->pmgd = pmgdHere; /* pointer to residue */
		/* this is also troubling CWVH cannot crease Voronoi without breakout of terms */
		/* we don't have the tracking informaiton to know what other residue is in the contact area ... ?*/
		paln->potential = -(totscore+solvscore)/(float)PDBtot; /* residue score */
		paln->potentialT1 = totscore;  /* break out Voronoi Total Term */
	    paln->potentialT2 = solvscore; /* break out Voronoi Solvent term */
	    paln->potentialCount  = atomsPerRes[i][1]- atomsPerRes[i][0] + 1; /* this should be the Voronoi useful atom count for residue ?? */
    	
		/* Create clean sublist structure and add to it */
		pdnHeadSublist = NULL;
		DValNodeAddPointer (&pdnHeadSublist, 0, paln); /* store in sublist */
		/* this adds a single "anonymous" potential value instead of a list of pairing residues CWVH */
        /* to CREASE Voronoi would require storing a list of pairwise residue contacts here */
    	
		/* Extract residue number in correct format */
		sscanf(pmgdHere->pcGraphNum,"%d",&resNum);
    	
		/* Add residue component to result list */
		DValNodeAddPointer (ppdnResult, (Int2) abs(resNum), pdnHeadSublist);
            
        pdnmgHere = pdnmgHere->next;
    }
    
    /* Now free the atom range array */
    for(cnt = 0; cnt <= pmmdMol->iResCount; cnt++){
        atomsPerRes[cnt] = MemFree(atomsPerRes[cnt]);
    }
    atomsPerRes = MemFree(atomsPerRes);
        
    totscore = 0.0;
    solvscore = 0.0;

    for(atomi=0; atomi<PDBtot; ++atomi) {
        if(PDB[atomi].class < MAX_ATOM_CLASS) {
            totscore += PDB[atomi].score;
            solvscore += PDB[atomi].solvscore;
        }
    }
    /* MUST break out the parameters and log them */

    TotalVPotential = -(totscore+solvscore)/(float)PDBtot; 
	TotScore = totscore;
    SolvScore = solvscore;
        
    /* Free all allocated memory */
    ca_index = MemFree(ca_index);
    ca_rec = MemFree(ca_rec);
    seed = MemFree(seed);
    box = MemFree(box);
    PDBlist = MemFree(PDBlist);
    PDB = MemFree(PDB);

    return ERR_SUCCESS;
}

FloatHi GetTotalVPotential()
{
    return  (FloatHi) TotalVPotential; /* is doubly solvent substracted */
}

FloatHi GetTotTermVPotential()
{
	return TotScore;
}
FloatHi GetSolvTermVPotential()
{
	return SolvScore;
}

FloatHi GetSSTermVPotential()
{
	return SSScore;
}

Int4 GetPDBatomsVPotential()
{
	return NumPDBatoms;
}

FloatHi GetASATotal()
{
	return ASATotal;
}


/******************************
* subroutine init_contfreq
******************************/

/* initialize contact matrix names and values */
/* same order as in subroutine assign_class */

void init_contfreq(struct contactinfo continfo[], float *contfreq) 
{
    int i,j;
    char resnames[170][4] = {"SAS",
        "ALA","ALA","ALA","ALA","ALA",
        "ARG","ARG","ARG","ARG","ARG","ARG","ARG","ARG","ARG","ARG","ARG",
        "ASN","ASN","ASN","ASN","ASN","ASN","ASN","ASN",
        "ASP","ASP","ASP","ASP","ASP","ASP","ASP","ASP",
        "CYS","CYS","CYS","CYS","CYS","CYS",
        "GLN","GLN","GLN","GLN","GLN","GLN","GLN","GLN","GLN",
        "GLU","GLU","GLU","GLU","GLU","GLU","GLU","GLU","GLU",
        "GLY","GLY","GLY","GLY",
        "HIS","HIS","HIS","HIS","HIS","HIS","HIS","HIS","HIS","HIS",
        "ILE","ILE","ILE","ILE","ILE","ILE","ILE","ILE",
        "LEU","LEU","LEU","LEU","LEU","LEU","LEU","LEU",
        "LYS","LYS","LYS","LYS","LYS","LYS","LYS","LYS","LYS",
        "MET","MET","MET","MET","MET","MET","MET","MET",
        "PHE","PHE","PHE","PHE","PHE","PHE","PHE","PHE","PHE","PHE","PHE",
        "PRO","PRO","PRO","PRO","PRO","PRO","PRO",
        "SER","SER","SER","SER","SER","SER",
        "THR","THR","THR","THR","THR","THR","THR",
        "TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP","TRP",
        "TYR","TYR","TYR","TYR","TYR","TYR","TYR","TYR","TYR","TYR","TYR","TYR",
        "VAL","VAL","VAL","VAL","VAL","VAL","VAL",
        "OXT","UNK"};
    char atomnames[170][ATOMNAME_LEN] = {" SAS",
        " N  "," CA "," C  "," O  "," CB ",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",
        " N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",
        " N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",
        " N  "," CA "," C  "," O  "," CB "," SG ",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",
        " N  "," CA "," C  "," O  ",
        " N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",
        " N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",
        " N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD ",
        " N  "," CA "," C  "," O  "," CB "," OG ",
        " N  "," CA "," C  "," O  "," CB "," OG1"," CG2",
        " N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2", 
        " N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",
        " N  "," CA "," C  "," O  "," CB "," CG1"," CG2",
        " OXT"," UNK"};			    

    /* assign names and set contacts to zero */
    for(i=0; i<170; ++i) {
        strncpy(continfo[i].res, resnames[i],3);
        continfo[i].res[3] = '\0';
        strncpy(continfo[i].name, atomnames[i],4);
        continfo[i].name[4] = '\0';
        continfo[i].atomcount = 0;
        continfo[i].areasum = 0.0;
        for(j=0; j<170; ++j) {
            *(contfreq+i*170+j) = 0;
        }
    }
    return;
}

/******************************
* subroutine  assign_class
* created 24/08/2001  BJM
******************************/

/* assigns an atom class for protein atoms. */
/* here, 170 contact types (0-169), SAS + each protein atom type, including terminal oxygen. */
/* class #169 is unknown contact type (ignore) */

/* alphabetical order: */
/* SAS:   0 */
/* ALA:   1-  5 N  CA C  O  CB */
/* ARG:   6- 16 N  CA C  O  CB CG  CD  NE  CZ  NH1 NH2 */
/* ASN:  17- 24 N  CA C  O  CB CG  OD1 ND2 */
/* ASP:  25- 32 N  CA C  O  CB CG  OD1 OD2 */
/* CYS:  33- 38 N  CA C  O  CB SG */
/* GLN:  39- 47 N  CA C  O  CB CG  CD  OE1 NE2 */
/* GLU:  48- 56 N  CA C  O  CB CG  CD  OE1 OE2 */
/* GLY:  57- 60 N  CA C  O */
/* HIS:  61- 70 N  CA C  O  CB CG  ND1 CD2 CE1 NE2 */
/* ILE:  71- 78 N  CA C  O  CB CG1 CG2 CD1 */
/* LEU:  79- 86 N  CA C  O  CB CG  CD1 CD2 */
/* LYS:  87- 95 N  CA C  O  CB CG  CD  CE  NZ */
/* MET:  96-103 N  CA C  O  CB CG  SD  CE */
/* PHE: 104-114 N  CA C  O  CB CG  CD1 CD2 CE1 CE2 CZ */
/* PRO: 115-121 N  CA C  O  CB CG  CD */
/* SER: 122-127 N  CA C  O  CB OG */
/* THR: 128-134 N  CA C  O  CB OG1 CG2 */
/* TRP: 135-148 N  CA C  O  CB CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2  */
/* TYR: 149-160 N  CA C  O  CB CG  CD1 CD2 CE1 CE2 CZ  OH */
/* VAL: 161-167 N  CA C  O  CB CG1 CG2 */

TrajErr assign_class(struct atom PDB[], long int PDBtot)
{
    long int atomi;  /* atom counter */
    int classnum;

    int numCA = 0;  /* for checking if side chains included */
    int numCB = 0;
    int numCG = 0;
    int numUNK = 0; /* number of unassigned atoms */

    for(atomi=0; atomi<PDBtot; ++atomi) {

        /* base number for amino acids */
        if(strncmp(PDB[atomi].res, "ALA", RESNAME_LEN-1) == 0) classnum = 1;
        else if(strncmp(PDB[atomi].res, "ARG", RESNAME_LEN-1) == 0) classnum = 6;
        else if(strncmp(PDB[atomi].res, "ASN", RESNAME_LEN-1) == 0) classnum = 17;
        else if(strncmp(PDB[atomi].res, "ASP", RESNAME_LEN-1) == 0) classnum = 25;
        else if(strncmp(PDB[atomi].res, "CYS", RESNAME_LEN-1) == 0) classnum = 33;
        else if(strncmp(PDB[atomi].res, "GLN", RESNAME_LEN-1) == 0) classnum = 39;
        else if(strncmp(PDB[atomi].res, "GLU", RESNAME_LEN-1) == 0) classnum = 48;
        else if(strncmp(PDB[atomi].res, "GLY", RESNAME_LEN-1) == 0) classnum = 57;
        else if(strncmp(PDB[atomi].res, "HIS", RESNAME_LEN-1) == 0) classnum = 61;
        else if(strncmp(PDB[atomi].res, "ILE", RESNAME_LEN-1) == 0) classnum = 71;
        else if(strncmp(PDB[atomi].res, "LEU", RESNAME_LEN-1) == 0) classnum = 79;
        else if(strncmp(PDB[atomi].res, "LYS", RESNAME_LEN-1) == 0) classnum = 87;
        else if(strncmp(PDB[atomi].res, "MET", RESNAME_LEN-1) == 0) classnum = 96;
        else if(strncmp(PDB[atomi].res, "PHE", RESNAME_LEN-1) == 0) classnum = 104;
        else if(strncmp(PDB[atomi].res, "PRO", RESNAME_LEN-1) == 0) classnum = 115;
        else if(strncmp(PDB[atomi].res, "SER", RESNAME_LEN-1) == 0) classnum = 122;
        else if(strncmp(PDB[atomi].res, "THR", RESNAME_LEN-1) == 0) classnum = 128;
        else if(strncmp(PDB[atomi].res, "TRP", RESNAME_LEN-1) == 0) classnum = 135;
        else if(strncmp(PDB[atomi].res, "TYR", RESNAME_LEN-1) == 0) classnum = 149;
        else if(strncmp(PDB[atomi].res, "VAL", RESNAME_LEN-1) == 0) classnum = 161;
        else {
            ErrPostEx(SEV_ERROR,4,1,"Residue for atom %ld not assigned, using class 0",atomi);
            PDB[atomi].class = 0;
            continue;
        }

        /* add values for atom position in list */
        if(strncmp(PDB[atomi].atomname, " N  ", ATOMNAME_LEN-1) == 0) {
            /*classnum += 0;  do nothing */
        } else if(strncmp(PDB[atomi].atomname, " CA ", ATOMNAME_LEN-1) == 0) {
            classnum += 1;
            ++numCA;
        } else if(strncmp(PDB[atomi].atomname, " C  ", ATOMNAME_LEN-1) == 0) {
            classnum += 2;
        } else if(strncmp(PDB[atomi].atomname, " O  ", ATOMNAME_LEN-1) == 0) {
            classnum += 3;
        } else if(strncmp(PDB[atomi].atomname, " CB ", ATOMNAME_LEN-1) == 0) {
            classnum += 4;
            ++numCB;
        } else if((strncmp(PDB[atomi].atomname, " CG ", ATOMNAME_LEN-1) == 0) || 
            (strncmp(PDB[atomi].atomname, " CG1", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " SG ", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " OG ", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " OG1", ATOMNAME_LEN-1) == 0)) {
                classnum += 5;
                ++numCG;
        } else if((strncmp(PDB[atomi].atomname, " CG2", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " CD ", ATOMNAME_LEN-1) == 0) ||
            ((strncmp(PDB[atomi].atomname," CD1", ATOMNAME_LEN-1) == 0) && (strncmp(PDB[atomi].res, "ILE", 3) != 0)) ||
            (strncmp(PDB[atomi].atomname, " OD1", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " ND1", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " SD ", ATOMNAME_LEN-1) == 0)) {
                classnum += 6;
        } else if((strncmp(PDB[atomi].atomname, " OD2", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " ND2", ATOMNAME_LEN-1) == 0) ||
            ((strncmp(PDB[atomi].atomname," CD1", ATOMNAME_LEN-1) == 0) && (strncmp(PDB[atomi].res, "ILE", 3) == 0)) ||
            (strncmp(PDB[atomi].atomname, " CD2", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " NE ", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " OE1", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " CE ", ATOMNAME_LEN-1) == 0)) {
                classnum += 7;
        } else if(((strncmp(PDB[atomi].atomname," CZ ", ATOMNAME_LEN-1) == 0) && (strncmp(PDB[atomi].res, "ARG", 3) == 0)) || 
            (strncmp(PDB[atomi].atomname, " NE1", ATOMNAME_LEN-1) == 0) ||
            ((strncmp(PDB[atomi].atomname," NE2", ATOMNAME_LEN-1) == 0) && (strncmp(PDB[atomi].res, "GLN", 3) == 0)) ||
            (strncmp(PDB[atomi].atomname, " OE2", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " CE1", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " NZ ", ATOMNAME_LEN-1) == 0)) {
                classnum += 8;
        } else if((strncmp(PDB[atomi].atomname, " NH1", ATOMNAME_LEN-1) == 0) ||
            ((strncmp(PDB[atomi].atomname," NE2", ATOMNAME_LEN-1) == 0) && (strncmp(PDB[atomi].res, "HIS", 3) == 0)) ||
            (strncmp(PDB[atomi].atomname, " CE2", ATOMNAME_LEN-1) == 0)) {
                classnum += 9;
        } else if((strncmp(PDB[atomi].atomname, " NH2", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " CZ ", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " CE3", ATOMNAME_LEN-1) == 0)) {
                classnum += 10;
        } else if((strncmp(PDB[atomi].atomname, " CZ2", ATOMNAME_LEN-1) == 0) ||
            (strncmp(PDB[atomi].atomname, " OH ", ATOMNAME_LEN-1) == 0)) {
                classnum += 11;
        } else if(strncmp(PDB[atomi].atomname, " CZ3", ATOMNAME_LEN-1) == 0) {
            classnum += 12;
        } else if(strncmp(PDB[atomi].atomname, " CH2", ATOMNAME_LEN-1) == 0) {
            classnum += 13;
        } else if(strncmp(PDB[atomi].atomname, " OXT", ATOMNAME_LEN-1) == 0) {
            classnum = MAX_ATOM_CLASS;  /* terminal oxygen of amino acid chain */
        } else {
            /*printf("atom type for atom %ld not assigned\n", atomi);*/
            classnum = MAX_ATOM_CLASS+1;  /* error class number */
            ++numUNK;
        }

        PDB[atomi].class = classnum;
    }

    if(numUNK > MAX_UNKNOWN_TOLERANCE) {
        ErrPostEx(SEV_ERROR,4,1,"Encountered %d unknown atoms, aborting current structure",numUNK);
        return ERR_FAIL;
    } 

   /* if(!numCG) {   TAKEN OUT AS IT CAUSES ERROR ON poly-GLY, poly-ALA CWVH 
        ErrPostEx(SEV_ERROR,4,1,"Structure has no side chains, aborting current structure",numUNK);
        return ERR_FAIL;
    } */ 

    return ERR_SUCCESS;
}



/*

ValNodePtr GetSolvatedAtoms()
{
 return pvnHead;
}

void FreeSolvatedAtoms()
{
	ValNodeFree(pvnHead);
	pvnHead = NULL;
	pvnCur = NULL;
}

*/
/*  this is intended to creat a linked list of solvated atoms without contacts for 
solvation profile PCA analysis */
/*
void GatherSolvatedAtoms(int AtomI) {

	if (pvnHead == NULL) {
        pvnHead = ValNodeAddPointer(&pvnHead, 0, PDB[AtomI].pmadAtom);
		pvnCur = pvnHead;
	}
	else  
		ValNodeAddPointer(&pvnCur, 0, (Pointer) PDB[AtomI].pmadAtom);
}

*/
/*
void SetWriteFile(FILE *fptr)
{
 fp = fptr;
}
This function writes out the unsolvated atoms as a PDB file specified in SetWriteFile() for testing

void WriteCapCenters(int NC, int AtomI) {

int cai;
float x, y , z;

x = (float) (PDB[AtomI].coor[0]);
y = (float) (PDB[AtomI].coor[1]);
z = (float) (PDB[AtomI].coor[2]);

fprintf(fp,"ATOM  %5d %s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %s\n",(int) CapCount,PDB[AtomI].atomname, PDB[AtomI].res,(int) AtomI,x,y,z,PDB[AtomI].atomname);
CapCount++;
*/




/****************************
* subroutine calc_score
****************************/

/* this subroutine calculates the atom-atom contact areas and SAS for a given protein. */
/* needs global variables 'score' and 'dim'. */

double calc_score(struct atom PDB[], long int PDBtot, int atomlist[], int windowSize)
{
    int    atomi;    /* atom counter */
    int    atomzero; /* current center atom */
    int    cai;      /* contact atom counter for atomzero */
    int    NC;       /* number of contacts around atomzero */
    int    NV;       /* number of vertices in polyhedron around atomzero */
    double areao;    /* area of atom zero */
    double areatot;  /* total contact area */
    float  rado;     /* radius of atomzero PLUS radius of water */
    float  SAS;      /* solvent exposed surface in square angstroms */
    char   surfatom; /* atom type, 'I' internal, 'S' surface */
    double atomscore;       /* score of the current atom */
    double SAStot = 0.0;    /* total SAS of protein */
    double AAtot  = 0.0;    /* total atom-atom contact area of protein */

    solvscore = 0.0;
    totscore = 0.0;

    for(atomi=0; atomi<PDBtot; ++atomi) {

        /* omit unclassed atoms (OXT, UNK) */
        if(PDB[atomi].class > (MAX_ATOM_CLASS-1)){
            continue;
        }

        atomscore = 0.0;
        /* ============= atom contact calculations ============= */
        atomzero = atomlist[atomi];
        rado = PDB[atomzero].radius + Rw;
        NC = get_contlist4(PDB, atomzero, contlist, PDBtot, rado, dim);
        NV = voronoi_poly2(PDB, atomzero, cont, rado, NC, contlist);
        /*
        if(NV == -1) { *//* error catching of subroutine */
        /* for scoring only, assign invalid atoms score of zero.*/ 
        /*PDB[atomzero].SAS = 0.0;
        PDB[atomzero].score = 0.0;
        continue;
        }
        */
        surfatom = order_faces(atomzero, poly, centerpt, rado, NC, NV, cont, ptorder);
		
        calc_areas(poly, centerpt, rado, NC, NV, cont, ptorder, atomzero);
        /*save_areas(cont, contlist, NC, atomzero); */

/******************CWVH **************************************************************************************************/
	 /*   if (surfatom == 'S')
            GatherSolvatedAtoms(atomi);
	 */
/******/

        /* check vertices to see if good ion site. */
        ion_sites(cont, poly, NV, atomzero);

        /* ==========  SAS contact area calculations =========== */
        areao = 4.0*PI*rado*rado;
        areatot = 0.0;
        for(cai=0; cai<NC; ++cai) {
            if(cont[cai].flag != 'X') {
                areatot += cont[cai].area;
            }
        }
		
        SAS = areao - areatot;
        if( SAS > 0.0) {  /* remove fractional negative values (~ -0.0000002) */
            PDB[atomzero].SAS = SAS;
            SAStot += SAS;
        } else {
            PDB[atomzero].SAS = 0.0;
        }

        /* store solvent score. */
        PDB[atomzero].solvscore =  score[PDB[atomzero].class][0]* PDB[atomzero].SAS;

        /* DEBUG */
        solvscore +=PDB[atomzero].solvscore;

        /* add atom contact scores to total score, excluding 'windowSize' adjacent residues */
        for(cai=0; cai<NC; ++cai) {
            if(cont[cai].flag != 'X') {
                /*if((abs(PDB[atomzero].resnum - PDB[cont[cai].index].resnum) > windowSize) 
                    && (PDB[atomzero].class < MAX_ATOM_CLASS)) {*/
                if((abs(PDB[atomzero].resnum - PDB[cont[cai].index].resnum) > windowSize) 
                    && (PDB[atomzero].class < MAX_ATOM_CLASS)
                    && (PDB[cont[cai].index].class < MAX_ATOM_CLASS)) {
                        atomscore += cont[cai].area * score[PDB[atomzero].class][PDB[cont[cai].index].class];
                        AAtot += cont[cai].area;
                        /*printf("atom %3d  %s%s gr%3d + contact %s%s gr%3d, area= %5.2f, score= %7.4f\n", atomzero,  */
                        /* PDB[atomzero].res, PDB[atomzero].atomname, PDB[atomzero].class, PDB[cont[cai].index].res,  */
                        /* PDB[cont[cai].index].atomname, PDB[cont[cai].index].class, cont[cai].area, */
                        /*cont[cai].area*score[PDB[atomzero].class][PDB[cont[cai].index].class]); */

                    }
            }
        }

        /*
        printf("atom%5d %s%s gr%3d,  AAscore= %7.3f  SAscore= %7.3f\n",atomzero,
        PDB[atomzero].res, PDB[atomzero].atomname, PDB[atomzero].class, 
        atomscore, PDB[atomzero].solvscore);
        */

        totscore += atomscore;
        PDB[atomzero].score = atomscore;
        /*PAUSE; */
    }

    /* total SAS lost on folding == total non-bonded contact area. */
    /*printf(" total non-bonded contact area/ #atoms = %6.2f\n", AAtot/PDBtot); */
	ASATotal = (FloatHi) AAtot;
    SSscore = score[0][0]*AAtot;
	SSScore = (FloatHi) SSscore;
    /*printf("  %5ld  %6.3f %6.3f %6.3f\n", PDBtot, totscore/(float)PDBtot,  */
    /* solvscore/(float)PDBtot, SSscore/(float)PDBtot); */
    return(totscore);
}




/****************************
* subroutine ion_sites
****************************/

/* uses global variables numsites and ionsites[]. */

void ion_sites(struct plane cont[], struct vertex poly[], int NV, int atomzero)
{
    int vi;
    int ai;
    char flag;
    int oxygens = 0; 
    int ode_nde_s = 0;  /* ASP or GLU OD/OE, HIS ND/NE, or CYS SG  */
    char *atomname;
    char *res;
    int initoxy = 0;
    int initONS = 0;
    int msi; /* metal (ion) site counter */
    int sai; /* site contact atom counter */
    char present; /* flag if contact already defined */
    double sqrdist;

    atomname = PDB[atomzero].atomname;
    res =  PDB[atomzero].res;
    /* backbone oxygen or ASP/GLU terminal oxygen */
    if((*(atomname+1) == 'O') && !(*(atomname+2) == 'G') && !(*(atomname+2) == 'H'))  {
        initoxy = 1;
        if((*(atomname+2) != ' ') && (!strncmp(res, "ASP", 3) || !strncmp(res, "GLU", 3))) { 
            initONS = 1;
        }
        /* return on any N's other than HIS ND/NE */
    } else if((*(atomname+1) == 'N') && (*(atomname+2) != ' ') && !strncmp(res, "HIS", 3)) {
        initONS = 1;
        /* keep if CYS SG */
    } else if((*(atomname+1) == 'S') && (*(atomname+2) == 'G')) {
        initONS = 1;
    } else {
        return;
    }

    /* Find possible metal site. */
    /* criteria:  3+ atoms of OD/E, ND/E, SG; or 4 atoms O. */

    for(vi=0; vi<NV; ++vi) {
        oxygens = initoxy;
        ode_nde_s = initONS;
        flag = 'Y';
        for(ai=0; ai<3; ++ai) {
            if(poly[vi].plane[ai] == -1) continue;
            atomname = PDB[cont[poly[vi].plane[ai]].index].atomname;
            res =  PDB[cont[poly[vi].plane[ai]].index].res;

            if((*(atomname+1) == 'O') && !(*(atomname+2) == 'G') && !(*(atomname+2) == 'H'))  {
                ++oxygens;
                if((*(atomname+2) != ' ') && (!strncmp(res, "ASP", 3) || !strncmp(res, "GLU", 3))) {
                    ++ode_nde_s;
                }
            } else if(*(atomname+1) == 'N') {
                if((*(atomname+2) != ' ') && !strncmp(res, "HIS", 3)) {
                    ++ode_nde_s;
                }
            } else if((*(atomname+1) == 'S') && (*(atomname+2) == 'G')) {
                ++ode_nde_s;
            }
        }

        if((oxygens == 4) || (ode_nde_s == 4)) {
            /* tentative ion site. Record coordinates and 4 defining atoms.  */
            ionsite[numsites].xi[0] = PDB[atomzero].coor[0] + poly[vi].xi[0];
            ionsite[numsites].xi[1] = PDB[atomzero].coor[1] + poly[vi].xi[1];
            ionsite[numsites].xi[2] = PDB[atomzero].coor[2] + poly[vi].xi[2];
            ionsite[numsites].atom[0] = atomzero;
            ionsite[numsites].atom[1] = cont[poly[vi].plane[0]].index;
            ionsite[numsites].atom[2] = cont[poly[vi].plane[1]].index;
            ionsite[numsites].atom[3] = cont[poly[vi].plane[2]].index;
            ionsite[numsites].maxnum = 4;
            /*compare with previous sites */
            for(msi=0; msi<numsites; ++msi) {
                sqrdist = (ionsite[numsites].xi[0]-ionsite[msi].xi[0])*(ionsite[numsites].xi[0]-ionsite[msi].xi[0]) 
                    + (ionsite[numsites].xi[1]-ionsite[msi].xi[1])*(ionsite[numsites].xi[1]-ionsite[msi].xi[1])
                    + (ionsite[numsites].xi[2]-ionsite[msi].xi[2])*(ionsite[numsites].xi[2]-ionsite[msi].xi[2]);
                if(sqrdist < 2.56) {
                    flag = 'N';
                    break;
                }
            }
            if(flag == 'Y') {
                /*DEBUG PRINT */
                /**
                printf("Site %3d atoms %5d%5d%5d%5d C:(%5.1f,%5.1f,%5.1f) D: %4.2f", 
                numsites, atomzero, 
                cont[poly[vi].plane[0]].index, cont[poly[vi].plane[1]].index, cont[poly[vi].plane[2]].index, 
                ionsite[numsites].xi[0], ionsite[numsites].xi[1], ionsite[numsites].xi[2],
                poly[vi].dist);
                printf(" A: %4s%4s %4s%4s %4s%4s %4s%4s.\n",
                PDB[atomzero].res, PDB[atomzero].atomname, 
                PDB[cont[poly[vi].plane[0]].index].res, PDB[cont[poly[vi].plane[0]].index].atomname, 
                PDB[cont[poly[vi].plane[1]].index].res, PDB[cont[poly[vi].plane[1]].index].atomname, 
                PDB[cont[poly[vi].plane[2]].index].res, PDB[cont[poly[vi].plane[2]].index].atomname); 	
                **/
                ++numsites;
            } else {
                /* else add contact atoms to previous ion site. */
                for(ai=0; ai<4; ++ai) {
                    present = 'N';
                    for(sai=0; sai<ionsite[msi].maxnum; ++sai) {
                        if(ionsite[msi].atom[sai] == ionsite[numsites].atom[ai]) {
                            present = 'Y';
                        }
                    }
                    if(present == 'N') {
                        ionsite[msi].atom[ionsite[msi].maxnum] = ionsite[numsites].atom[ai];
                        ++ionsite[msi].maxnum;
                        /* average site location */
                        ionsite[msi].xi[0] = (ionsite[msi].xi[0] + ionsite[numsites].xi[0])/2.0;
                        ionsite[msi].xi[1] = (ionsite[msi].xi[1] + ionsite[numsites].xi[1])/2.0;
                        ionsite[msi].xi[2] = (ionsite[msi].xi[2] + ionsite[numsites].xi[2])/2.0;
                    }
                }
            }
        }
    }
    return;
}



/****************************
* subroutine get_firstvert
****************************/

void get_firstvert(struct plane cont[], int *planeA, int *planeB, int *planeC, int NC, int atomzero)
{
    long int seedi;
    int cai;
    double mindist;
    double ptA[3];
    double vectA[4];     /* 4 values so it can be used as a plane equation as well */
    double vt;           /* vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt */
    double vtmin;        /* minimum value for vt  */
    double vtdiv;        /* check for division by zero in vt calculation  */
    double temppt[3];
    double ptdist;

    /* NEW - get previous seed vertex */
    seedi = atomzero*3;
    *planeA = -1;
    *planeB = -1;
    *planeC = -1;
    if(seed[seedi] != -1) {
        for(cai=0; cai<NC; ++cai) {
            if(cont[cai].index == seed[seedi]) {
                *planeA = cai;
            } else if(cont[cai].index == seed[seedi+1]) {
                *planeB = cai;
            } else if(cont[cai].index == seed[seedi+2]) {
                *planeC = cai;
            }
        }
    }
    if((*planeA != -1)&&(*planeB != -1)&&(*planeC != -1)) {
        return;
    } else {

        /* -------------  find initial edge, on plane closest to origin  ------------- */
        mindist = 9.9e+9;  /* dummy value */
        for(cai=0; cai<NC; ++cai) {
            if(cont[cai].dist < mindist) {
                mindist = cont[cai].dist;
                *planeA = cai;
            }
        }

        mindist = 9.9e+9;  /* dummy value */
        for(cai=0; cai<NC+4; ++cai) {
            if(cai != *planeA) {
                vectA[0] = cont[*planeA].Ai[1]*cont[cai].Ai[2] - cont[*planeA].Ai[2]*cont[cai].Ai[1];
                vectA[1] = cont[*planeA].Ai[2]*cont[cai].Ai[0] - cont[*planeA].Ai[0]*cont[cai].Ai[2];
                vectA[2] = cont[*planeA].Ai[0]*cont[cai].Ai[1] - cont[*planeA].Ai[1]*cont[cai].Ai[0];
                vectA[3] = 0.0;
                if(solve_3x3(cont[*planeA].Ai, cont[cai].Ai, vectA, temppt) != -1) {
                    ptdist = sqrt(temppt[0]*temppt[0] + temppt[1]*temppt[1] + temppt[2]*temppt[2]);
                    if(ptdist < mindist) {
                        *planeB = cai;
                        mindist = ptdist;
                        ptA[0] = temppt[0];
                        ptA[1] = temppt[1];
                        ptA[2] = temppt[2];
                    }
                }
            }
        }
        /* recalc vector normal to planes A and B */
        vectA[0] = cont[*planeA].Ai[1]*cont[*planeB].Ai[2] - cont[*planeA].Ai[2]*cont[*planeB].Ai[1];
        vectA[1] = cont[*planeA].Ai[2]*cont[*planeB].Ai[0] - cont[*planeA].Ai[0]*cont[*planeB].Ai[2];
        vectA[2] = cont[*planeA].Ai[0]*cont[*planeB].Ai[1] - cont[*planeA].Ai[1]*cont[*planeB].Ai[0];
        vectA[3] = 0.0;

        /* get starting vertex on polyhedron */
        vtmin = 9.9e+9;   /* dummy value */
        for(cai=0; cai<NC+4; ++cai) {
            if((cai != *planeA) && (cai != *planeB)) {
                vtdiv = (cont[cai].Ai[0]*vectA[0] +cont[cai].Ai[1]*vectA[1] +cont[cai].Ai[2]*vectA[2]);
                if(vtdiv != 0.0) {
                    vt = -(cont[cai].Ai[0]*ptA[0] +cont[cai].Ai[1]*ptA[1] +cont[cai].Ai[2]*ptA[2] +cont[cai].Ai[3])/vtdiv;
                    if(fabs(vt) < vtmin) {
                        vtmin = fabs(vt);
                        *planeC = cai;
                    }
                }
            }
        }
    }
    return;
}


/*************************
* subroutine save_seeds
*************************/

/* saves starting vertices for atoms to be done */

void save_seeds(struct plane cont[], struct vertex poly[], int NV, int atomzero)
{
    int vi;
    int seedi;

    for(vi=0; vi<NV; ++vi) {
        if(poly[vi].plane[2] != -1) {
            seedi = 3*cont[poly[vi].plane[0]].index;
            if(seed[seedi] == -1) {
                seed[seedi] = atomzero;
                seed[seedi+1] = cont[poly[vi].plane[1]].index;
                seed[seedi+2] = cont[poly[vi].plane[2]].index;
            }
            seedi = 3*cont[poly[vi].plane[1]].index;
            if(seed[seedi] == -1) {
                seed[seedi] = atomzero;
                seed[seedi+1] = cont[poly[vi].plane[0]].index;
                seed[seedi+2] = cont[poly[vi].plane[2]].index;
            }
            seedi = 3*cont[poly[vi].plane[2]].index;
            if(seed[seedi] == -1) {
                seed[seedi] = atomzero;
                seed[seedi+1] = cont[poly[vi].plane[0]].index;
                seed[seedi+2] = cont[poly[vi].plane[1]].index;
            }
        }
    }  
    return;
}



/*******************************
* subroutine add_vedge
*******************************/

/* adds a new edge to the edge list. Direction of vector is tested so that */
/* it points away from the startpt, along the body of the polyhedron. */

/* stores information in variable vedge[edgenum]. */

void add_vedge(struct edgevector vedge[], int edgenum, struct plane cont[], int plane0, int plane1, 
               int testplane, struct vertex poly[], int startpt)
{
    double testpt[3];
    double testval;

    vedge[edgenum].V[0] = cont[plane0].Ai[1]*cont[plane1].Ai[2] - cont[plane0].Ai[2]*cont[plane1].Ai[1];
    vedge[edgenum].V[1] = cont[plane0].Ai[2]*cont[plane1].Ai[0] - cont[plane0].Ai[0]*cont[plane1].Ai[2];
    vedge[edgenum].V[2] = cont[plane0].Ai[0]*cont[plane1].Ai[1] - cont[plane0].Ai[1]*cont[plane1].Ai[0];
    vedge[edgenum].startpt = startpt;
    vedge[edgenum].endpt = -1;  /* flag, edge not completed. */
    vedge[edgenum].plane[0] = plane0;
    vedge[edgenum].plane[1] = plane1;
    vedge[edgenum].startplane = testplane;
    vedge[edgenum].arc = '.'; /* dummy value. */

    /* test direction of vector */
    testpt[0] = poly[startpt].xi[0] + vedge[edgenum].V[0];
    testpt[1] = poly[startpt].xi[1] + vedge[edgenum].V[1];
    testpt[2] = poly[startpt].xi[2] + vedge[edgenum].V[2];

    testval = cont[testplane].Ai[0]*testpt[0] +cont[testplane].Ai[1]*testpt[1] 
    + cont[testplane].Ai[2]*testpt[2] +cont[testplane].Ai[3];
    if(testval > 0.0) { /* vector is in wrong direction */
        vedge[edgenum].V[0] = -vedge[edgenum].V[0];
        vedge[edgenum].V[1] = -vedge[edgenum].V[1];
        vedge[edgenum].V[2] = -vedge[edgenum].V[2];
    }
    return;
}

/*******************************
* subroutine add_vertex
*******************************/

/* add polyhedron vertex to local list for current atom contacts */

int add_vertex(struct vertex poly[], int vn, double coor[], int planeA, int planeB, int planeC)
{ 
    poly[vn].xi[0] = coor[0];
    poly[vn].xi[1] = coor[1];
    poly[vn].xi[2] = coor[2];
    poly[vn].plane[0] = planeA;
    poly[vn].plane[1] = planeB;
    poly[vn].plane[2] = planeC;
    poly[vn].dist = sqrt(coor[0]*coor[0] +coor[1]*coor[1] +coor[2]*coor[2]);
    return(0);
}

/*********************************
* function order_faces          
*********************************/

/* output is array ptorder, the order of points around each face. */
/* return values are 'I' for internal atom, 'S' for surface atom. */

char order_faces(int atomzero,struct vertex poly[], struct vertex centerpt[], float rado,
                 int NC, int NV, struct plane cont[], struct ptindex ptorder[])
{
    int planeX;             /* current plane, 1 to NC */
    int planeY;             /* second plane, to get adjacent points on polygon */
    int surfcount;          /* number of points defining a given plane of contact */
    int vi, vi2;            /* vertices counter */
    int tempsi;             /* temporary storage for exchanging indices */
    double tempcos;         /* for exchanging cosines */
    double cos10X[50];      /* cos of angle between points poly[1],poly[0],and poly[X] */
    double temppt[3];       /* temp coordinates of new arc point */
    char surfatom;          /* return value: 'I' internal atom, 'S' surface atom */

    surfatom = 'I'; /* internal atom */
    for(vi=0; vi<NV; ++vi) {
        if(poly[vi].plane[2] == -1) {
            surfatom = 'S'; /* surface atom  */
            break;
        }
    }
    
    /* Initialize the cos angles array to zeroes */
    for(vi=0;vi<50;vi++){
        cos10X[vi] = 0.0;
    }
    
    /* for surface calculation only */
    /* if(surfatom == 'I') return('I'); */

    for(planeX=0; planeX < NC; ++planeX) {
        if(cont[planeX].flag == 'X') { /* hidden */
            ptorder[planeX].numpts = 0;
            continue;
        }

        surfcount=0;
        for(vi=0; vi<NV; ++vi) {
            /* index all points comprising surface for planeX */
            if((poly[vi].plane[0]==planeX) || (poly[vi].plane[1]==planeX) || (poly[vi].plane[2]==planeX)) {
                ptorder[planeX].pt[surfcount] = vi;
                ++surfcount;
            }
        }

        if(surfcount > 2) {

            /*-------------------------------------------------------------------*/
            /* three or more surface points, need additional point so no arcs    */
            /* are greater than 180 deg. Make point opposite first pt on an arc. */
            /* (if no points on arc, extra point isn't needed.)                  */
            /*-------------------------------------------------------------------*/

            for(vi=0; vi<surfcount; ++vi) {
                if(poly[ptorder[planeX].pt[vi]].plane[2] == -1) { /* on arc, calc pt. opposite */

                    /* get coordinates of point */
                    temppt[0] = 2.0*centerpt[planeX].xi[0] - poly[ptorder[planeX].pt[vi]].xi[0];
                    temppt[1] = 2.0*centerpt[planeX].xi[1] - poly[ptorder[planeX].pt[vi]].xi[1];
                    temppt[2] = 2.0*centerpt[planeX].xi[2] - poly[ptorder[planeX].pt[vi]].xi[2];

                    /*keep point if it's valid */
                    if(test_point(temppt, cont, NC, rado, planeX, -1, -1) == 'Y') {
                        poly[NV].xi[0] = temppt[0];
                        poly[NV].xi[1] = temppt[1];
                        poly[NV].xi[2] = temppt[2];
                        poly[NV].plane[0] = planeX;
                        poly[NV].plane[1] = -1;
                        poly[NV].plane[2] = -1;
                        poly[NV].dist = rado;
                        ptorder[planeX].pt[surfcount] = NV;
                        ++surfcount;
                        ++NV;
                    }
                    break;
                }
            }
        }
        ptorder[planeX].numpts = surfcount;

        if(surfcount > 3) {
            /* get two points on same line (two common planes). */
            /* all points already share one plane (planeX), find another. */
            if(poly[ptorder[planeX].pt[0]].plane[0] == planeX) {
                planeY = poly[ptorder[planeX].pt[0]].plane[1];
            } else {
                planeY = poly[ptorder[planeX].pt[0]].plane[0];
            }

            /*find another point on the same line (2 common planes) */
            for(vi=1; vi<surfcount; ++vi) {
                if((poly[ptorder[planeX].pt[vi]].plane[0]==planeY) || (poly[ptorder[planeX].pt[vi]].plane[1]==planeY) 
                    || (poly[ptorder[planeX].pt[vi]].plane[2]==planeY)) {
                        break;
                    }
            }    

            /*swap index for pt[1] and pt[vi], so points 0 and 1 are on same line */
            tempsi = ptorder[planeX].pt[vi];
            ptorder[planeX].pt[vi] = ptorder[planeX].pt[1];
            ptorder[planeX].pt[1] = tempsi;

            /* calculate cosine between points indexed 1,0,X */
            for(vi=2; vi<surfcount; ++vi) {
                if (vi < 50){
                    cos10X[vi] = cosPQR(poly[ptorder[planeX].pt[1]].xi, poly[ptorder[planeX].pt[0]].xi, poly[ptorder[planeX].pt[vi]].xi);
                }
            }

            /* order by cosines, decreasing order */
            for(vi=2; vi<surfcount-1; ++vi) {
                for(vi2=vi+1; vi2<surfcount; ++vi2) {
                    /*if(cos10X[vi] < cos10X[vi2]) {*/
                    /* Check if indices need to be swapped, and the validity of all array boundaries */
                    if((cos10X[vi] < cos10X[vi2]) && (vi < 50) && (vi2 < 50) && (planeX < 200)) {
                        /* swap indices if points in wrong order */
                        tempsi = ptorder[planeX].pt[vi];
                        ptorder[planeX].pt[vi] = ptorder[planeX].pt[vi2];
                        ptorder[planeX].pt[vi2] = tempsi;
                        tempcos = cos10X[vi];
                        cos10X[vi] = cos10X[vi2];
                        cos10X[vi2] = tempcos;
                    }
                }
            }
        }
    }
    return(surfatom);
}


/**************************
* subroutine calc_areas
**************************/

void calc_areas(struct vertex poly[], struct vertex centerpt[], float rado, int NC, int NV, 
struct plane cont[], struct ptindex ptorder[], int atomzero)
{
    char   engflag;         /* ='Y' if an engulfing plane is present */
    int    planeX;          /* current plane, 1 to NC */
    int    NP;              /* number of points on face */
    int    vi;              /* vertices counter */
    double area;            /* area of polygon */
    int    pa1, pa2;        /* possible common planes */
    int    commplane;       /* plane shared by adjacent points (not planeX) */
    int    epi;             /* engulfing plane counter */
    int    engplane[4];     /* index to engulfed planes */
    double maxSAS;          /* maximum solvent exposed surface, no atoms other than engulfing. */
    struct vertex ptB, ptC; /* arc point intersections */
    int    currpt, nextpt;
    double cosNN1[40];      /* angle between vertex N and vertex N+1  */
    double cosNzero[40];    /* angle between vertex N and vertex N+1  */
    double tanprod;         /* product of tangents */
    int    v0, va, vb;      /* vertices for arc calculation */

    double U,V,W,X;
    double tansqrS, tansqrSA, tansqrSB, tansqrSC;

    engflag = 'N';
    epi = 0;

    /*RESET AREAS TO ZERO */
    for(planeX=0; planeX < NC; ++planeX) {
        cont[planeX].area = 0.0;
        if(cont[planeX].flag == 'E') {
            engflag = 'Y';
            if (epi < 4){
                engplane[epi] = planeX;
            }
            ++epi;
        }
    }

    if(engflag == 'Y') { /* engulfing plane correction - project points onto sphere surface. */
        project_points(poly, centerpt, rado, NC, NV, cont);
    }

    /* ---------------------------- */
    /* calculate area for each face */
    /* ---------------------------- */

    for(planeX=0; planeX < NC; ++planeX) {
        NP = ptorder[planeX].numpts;
        area = 0.0;
        if(cont[planeX].flag == 'X') {
            continue;
        }

        /* if there are no points on a valid contact, area is spherical cap */
        if(NP == 0) {
            if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') {
                cont[planeX].area = 2.0*PI*rado*(rado-centerpt[planeX].dist);
            }
        } else if(NP == 2) {  /* only two contact points, check which part of arc */
            if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') { /* area is (cap - arc) */
                cont[planeX].area = 2.0*PI*rado*(rado-centerpt[planeX].dist)
                    - spherical_arc(centerpt[planeX], poly[ptorder[planeX].pt[0]], 
                    poly[ptorder[planeX].pt[1]], rado);
            } else {            /* area is arc. */
                cont[planeX].area = spherical_arc(centerpt[planeX], poly[ptorder[planeX].pt[0]], 
                    poly[ptorder[planeX].pt[1]], rado);
            }
        } else {

            /* ------ calculate cosines and angles ------ */
            for(vi=0; vi<NP; ++vi) {
                v0 = ptorder[planeX].pt[0];
                va = ptorder[planeX].pt[vi];
                vb = ptorder[planeX].pt[(vi+1)%NP];

                /* calculate cosines between adjacent vertices */
                if (vi < 40){
                    cosNN1[vi] = (( poly[va].xi[0]*poly[vb].xi[0] + poly[va].xi[1]*poly[vb].xi[1] 
                    + poly[va].xi[2]*poly[vb].xi[2] ) / (poly[va].dist*poly[vb].dist));

                    /* calculate cosines between vertex zero and vertex 'vi' */
                    if(vi != 0) {
                        cosNzero[vi] = ((poly[v0].xi[0]*poly[va].xi[0] + poly[v0].xi[1]*poly[va].xi[1] 
                        + poly[v0].xi[2]*poly[va].xi[2] ) / (poly[v0].dist*poly[va].dist));
                    }
                }
            }

            /* ----- calculate area of triangles in face ----- */
            for(vi=1; vi<(NP-1); ++vi) {
                U = sqrt((1+cosNzero[vi])*(1+cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
                V = sqrt((1-cosNzero[vi])*(1-cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
                W = sqrt((1-cosNzero[vi])*(1+cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
                X = sqrt((1+cosNzero[vi])*(1-cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
                tansqrS  = (1-U+V+W+X)/(1+U-V-W-X);
                tansqrSA = (1-U-V-W+X)/(1+U+V+W-X);
                tansqrSB = (1-U-V+W-X)/(1+U+V-W+X);
                tansqrSC = (1-U+V-W-X)/(1+U-V+W+X);
                tanprod = sqrt(tansqrS*tansqrSA*tansqrSB*tansqrSC);
                if(tanprod > 0.0) {
                    area += 4.0*rado*rado*atan(sqrt(tanprod));
                }
            }

            /* ----- add area of arc segments  ----- */
            for(vi=0; vi<NP; ++vi) {
                va = ptorder[planeX].pt[vi];
                vb = ptorder[planeX].pt[(vi+1)%NP];

                /*check if adjacent points are arc segments */
                if((poly[va].plane[2] == -1) && (poly[vb].plane[2] == -1)) {
                    /* if on same two planes, no arc. */
                    if((poly[va].plane[0]+poly[va].plane[1]) != (poly[vb].plane[0]+poly[vb].plane[1])) {
                        area += spherical_arc(centerpt[planeX], poly[va], poly[vb], rado);
                    }
                }
            }
            cont[planeX].area = area;
        }
    }

    /* -------------------------------------------------------- */
    /*  add correction terms for engulfing planes, if required */
    /* -------------------------------------------------------- */

    if(engflag == 'Y') {
        for(planeX=0; planeX < NC; ++planeX) {
            if(cont[planeX].flag != 'E') {
                continue;
            }

            NP = ptorder[planeX].numpts;
            for(vi=0; vi<NP; ++vi) {
                currpt = ptorder[planeX].pt[vi];
                nextpt = ptorder[planeX].pt[(vi+1)%NP];

                /* find common second plane, if any. */
                if(poly[currpt].plane[0] == planeX) {
                    pa1 = poly[currpt].plane[1];
                    pa2 = poly[currpt].plane[2];
                } else {
                    pa1 = poly[currpt].plane[0];
                    if(poly[currpt].plane[1] == planeX) {
                        pa2 = poly[currpt].plane[2];
                    } else {
                        pa2 = poly[currpt].plane[1];
                    }
                }

                if((pa1 == poly[nextpt].plane[0]) || (pa1 == poly[nextpt].plane[1])
                    || (pa1 == poly[nextpt].plane[2])) {
                        commplane = pa1;
                    } else if((pa2 == poly[nextpt].plane[0]) || (pa2 == poly[nextpt].plane[1])
                        || (pa2 == poly[nextpt].plane[2])) { 
                            commplane = pa2;
                        } else {
                            continue;
                        }
                        if((commplane != -1) && (cont[commplane].flag != 'E')) {
                            /* add correction to commplane area. here centerpt is from engulfing plane. */
                            cont[commplane].area += spherical_arc(centerpt[planeX], poly[currpt], poly[nextpt], rado);
                            if(NP == 2) break;  /* otherwise would repeat adding area */
                        }
            }
        }

        /* ----------------------------------------------- */
        /* ------ calculate engulfed contact areas ------- */
        /* ----------------------------------------------- */

        if(epi == 1) {
            cont[engplane[0]].area = 2.0*PI*rado*(rado+cont[engplane[0]].dist);
        } else if(epi == 2) { 
            if(solve_2xS(cont[engplane[0]],cont[engplane[1]], rado, ptB.xi, ptC.xi)== -1) {
                cont[engplane[0]].area = 2.0*PI*rado*rado;
                cont[engplane[1]].area = 2.0*PI*rado*rado;
            } else {
                ptB.dist = rado;
                ptC.dist = rado;
                maxSAS = spherical_arc(centerpt[engplane[0]], ptB, ptC, rado);
                maxSAS += spherical_arc(centerpt[engplane[1]], ptB, ptC, rado);
                cont[engplane[0]].area = 2.0*PI*rado*rado - 0.5*maxSAS;
                cont[engplane[1]].area = 2.0*PI*rado*rado - 0.5*maxSAS;
            }
        } else if(epi>=3) {
            /* no exposed surface if there are three or more engulfing contacts  */
            for(planeX=0; planeX<NC; ++planeX) {
                if(cont[planeX].flag == 'E') {
                    cont[planeX].area = 4.0*PI*rado*rado/epi;
                } else {
                    cont[planeX].area = 0.0;
                }
            }
        }
    }
    return;
}


/***************************
* subroutine save_areas
***************************/

void save_areas(struct plane cont[], struct contactlist contlist[], int NC, int atomzero)
{
    int cai;
    long int currindex, previndex, nextindex;

    if(numcarec > (ca_recsize-100)) {
        ca_recsize += 10000;
        ca_rec = MemMore(ca_rec, ca_recsize*sizeof(struct ca_struct));
        if(!ca_rec) {
            ErrPostEx(SEV_FATAL,1,1,"Memory re-allocation error for ca_rec");
        }
    }

    /* first overwrite previous records */
    cai = 0;
    currindex = ca_index[atomzero];
    previndex = -1;
    while((currindex != -1) && (cai<NC)) {
        if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
            ca_rec[currindex].atom = cont[cai].index;
            ca_rec[currindex].area = cont[cai].area;
            ca_rec[currindex].dist = contlist[cai].dist;
            nextindex = ca_rec[currindex].prev; /* next index is prev record from old list */
            ca_rec[currindex].prev = previndex;
            previndex = currindex;
            currindex = nextindex;
        }
        ++cai;
    }
    ca_index[atomzero] = previndex;

    /* then add new records */
    while(cai<NC) {
        if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
            ca_rec[numcarec].atom = cont[cai].index;
            ca_rec[numcarec].area = cont[cai].area;
            ca_rec[numcarec].dist = contlist[cai].dist;
            ca_rec[numcarec].prev = ca_index[atomzero];
            ca_index[atomzero] = numcarec;
            ++numcarec;
        }
        ++cai;
    }
    PDB[atomzero].done = 'Y';

    /* index contacts for atoms not yet done */
    for(cai=0; cai<NC; ++cai) {
        if(PDB[cont[cai].index].done != 'Y') {
            if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
                ca_rec[numcarec].atom = atomzero;
                ca_rec[numcarec].prev = ca_index[cont[cai].index];
                ca_index[cont[cai].index] = numcarec;
                ++numcarec;
            }
        }
    }
    return;
}

/*****************************
* subroutine project_points
*****************************/

/* this subroutine corrects for engulfing atoms. The center of the engulfing plane */
/* contact is used as the center of projection instead of the center of the atom. */
/* points already on the surface are not moved, preserving the SAS. */

void project_points(struct vertex poly[], struct vertex centerpt[], float rado, int NC, 
                    int NV, struct plane cont[])
{
    int epi;               /* engulfing plane counter */
    int cai;               /* contact atom counter */
    int engplane[4];       /* index to engulfing planes */
    double projpt[3];      /* center point for projecting internal pts onto surface of sphere */
    double pt0[3], pt1[3]; /* temporary points for intersection solutions */
    double V[3];           /* vector from projection point to vertex */
    double *P;             /* pointer for vertex coordinates to be projected */
    double a,b,c,k;        /* for solving quadratic eqn for point projection */
    int vi;                /* vertex counter */

    /* count and mark engulfing planes */
    epi = 0;
    for(cai=0; cai<NC; ++cai) {
        if(cont[cai].flag == 'E') {
            engplane[epi] = cai;
            ++epi;
        }
    }

    /* get projpt[] for projecting points to surface */
    if(epi == 1) {
        projpt[0] = centerpt[engplane[0]].xi[0];
        projpt[1] = centerpt[engplane[0]].xi[1];
        projpt[2] = centerpt[engplane[0]].xi[2];
    } else if(epi == 2) {
        solve_2xS(cont[engplane[0]],cont[engplane[1]],rado, pt0, pt1);
        projpt[0] = (pt0[0]+pt1[0])/2;
        projpt[1] = (pt0[1]+pt1[1])/2;
        projpt[2] = (pt0[2]+pt1[2])/2;
    } else {
        solve_3x3(cont[engplane[0]].Ai, cont[engplane[1]].Ai, cont[engplane[2]].Ai, pt0);
        projpt[0] = pt0[0];
        projpt[1] = pt0[1];
        projpt[2] = pt0[2];
    }

    for(vi=0; vi<NV; ++vi) {
        if(poly[vi].plane[2] != -1) {
            /* project point to intersection of engulfing plane(s) and surface of sphere */
            P = poly[vi].xi;
            V[0] = P[0] - projpt[0];
            V[1] = P[1] - projpt[1];
            V[2] = P[2] - projpt[2];
            a = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
            b = 2*(P[0]*V[0] +P[1]*V[1] +P[2]*V[2]);
            c = P[0]*P[0] + P[1]*P[1] + P[2]*P[2] - rado*rado;  /* c is < 0  */
            k = (sqrt(b*b - 4.0*a*c) - b)/(2*a);                /* k is > 0 */
            P[0] += k*V[0];
            P[1] += k*V[1];
            P[2] += k*V[2];      
            poly[vi].dist = rado;
        }
    }
    return;
}


/************************
* function test_point     
************************/

/* this function tests a given point versus a set of plane constraints */
/* it returns a value of 'Y' if the point is OK, and 'N' if there */
/* was a violation of the plane inequality. */
/* ptX =(xo, yo, zo); if Axo+Byo+Czo+D > 0, point is behind plane (hidden) */
/* planes A,B,C are the planes that the point lies on, don't test. */

char test_point(double ptX[], struct plane cont[], int NC, float rado, int planeA, int planeB, int planeC)
{
    int cp;      /* counter, number of planes */

    /* if point is not behind any plane, keep point. */
    for(cp=0; cp<NC; ++cp) {
        /*if pt is on current plane, get next plane */
        if((cp != planeA) && (cp != planeB) && (cp != planeC) && (cont[cp].flag != 'X')) {
            if((cont[cp].Ai[0]*ptX[0] + cont[cp].Ai[1]*ptX[1] + cont[cp].Ai[2]*ptX[2] + cont[cp].Ai[3]) > 0.0) {
                /*point is behind plane, not on polyhedron. */
                return('N');
            }
        }
    }
    return('Y');
}


/****************************
* function solve_3x3   
****************************/

/* determines the intersection of three planes */
/* (solves a system of three linear equations and 3 unknowns) */
/* input is 3 four element arrays, representing Ax+By+Cz+D=0 */
/* output is a three element array, (xi,xj,xk). */

int solve_3x3(double eq0[], double eq1[], double eq2[], double pt[])
{
    double cof00, cof01, cof02;   /* matrix cofactors */
    double cof10, cof11, cof12;
    double cof20, cof21, cof22;
    double det;                   /* determinant of matrix */

    cof00 =  eq1[1]*eq2[2] - eq2[1]*eq1[2];
    cof01 = -eq1[0]*eq2[2] + eq2[0]*eq1[2];
    cof02 =  eq1[0]*eq2[1] - eq2[0]*eq1[1];

    cof10 = -eq0[1]*eq2[2] + eq2[1]*eq0[2];
    cof11 =  eq0[0]*eq2[2] - eq2[0]*eq0[2];
    cof12 = -eq0[0]*eq2[1] + eq2[0]*eq0[1];

    cof20 =  eq0[1]*eq1[2] - eq1[1]*eq0[2];
    cof21 = -eq0[0]*eq1[2] + eq1[0]*eq0[2];
    cof22 =  eq0[0]*eq1[1] - eq1[0]*eq0[1];

    det = eq0[0]*cof00 + eq0[1]*cof01 + eq0[2]*cof02;
    if(det == 0.0) {
        /*printf("no solution for equation set, determinant is zero\n"); */
        return(-1);
    } else {
        pt[0] = -(eq0[3]*cof00 + eq1[3]*cof10 + eq2[3]*cof20)/det;
        pt[1] = -(eq0[3]*cof01 + eq1[3]*cof11 + eq2[3]*cof21)/det;
        pt[2] = -(eq0[3]*cof02 + eq1[3]*cof12 + eq2[3]*cof22)/det;
        return(0);
    }
}

/**********************************
* Function solve_2xS       
* revised 19/02/2001  BJM
***********************************/

/* determines the intersection of two planes and a sphere radius 'rad' */
/* input is 2 four element arrays, representing Ax+By+Cz+D=0 */
/* plus the radius of a sphere centered on the origin. */
/* output is two three element arrays pt0 and pt1, (xi,xj,xk). */
/* return value is -1 if no real solution exists. */

int solve_2xS(struct plane eq0, struct plane eq1, float rado, double pt0[], double pt1[])
{
    double eq2[3];              /* eqn of plane through (0,0,0) and perp. to other two */
    double cof00, cof01, cof02; /* matrix cofactors */
    double cof10, cof11, cof12; /* (don't need cof20, cof21, cof22) */
    double det;                 /* determinant of matrix */
    double avgpt[3];            /* average of two solution points */
    double t;                   /* parameter in eqn of line: x=xo+At, y=yo+Bt, z=zo+Ct.  */
    int xi;                     /* coordinate counter */

    eq2[0] = eq0.Ai[1]*eq1.Ai[2] - eq0.Ai[2]*eq1.Ai[1];
    eq2[1] = eq0.Ai[2]*eq1.Ai[0] - eq0.Ai[0]*eq1.Ai[2];
    eq2[2] = eq0.Ai[0]*eq1.Ai[1] - eq0.Ai[1]*eq1.Ai[0];

    cof00 =  eq1.Ai[1]*eq2[2] - eq2[1]*eq1.Ai[2];
    cof01 = -eq1.Ai[0]*eq2[2] + eq2[0]*eq1.Ai[2];
    cof02 =  eq1.Ai[0]*eq2[1] - eq2[0]*eq1.Ai[1];

    cof10 = -eq0.Ai[1]*eq2[2] + eq2[1]*eq0.Ai[2];
    cof11 =  eq0.Ai[0]*eq2[2] - eq2[0]*eq0.Ai[2];
    cof12 = -eq0.Ai[0]*eq2[1] + eq2[0]*eq0.Ai[1];

    det = eq2[0]*eq2[0] + eq2[1]*eq2[1] + eq2[2]*eq2[2];
    if(det == 0) {
        /*printf("no solution in solve_2xS\n"); */
        return(-1);
    }

    avgpt[0] = -(eq0.Ai[3]*cof00 + eq1.Ai[3]*cof10)/det;
    avgpt[1] = -(eq0.Ai[3]*cof01 + eq1.Ai[3]*cof11)/det;
    avgpt[2] = -(eq0.Ai[3]*cof02 + eq1.Ai[3]*cof12)/det;

    t = (rado*rado-avgpt[0]*avgpt[0]-avgpt[1]*avgpt[1]-avgpt[2]*avgpt[2])/det;
    if(t<0) {
        return(-1);
    } else {
        t = sqrt(t);
    }
    for(xi=0; xi<3; ++xi) {
        pt0[xi] = avgpt[xi] + t*eq2[xi];
        pt1[xi] = avgpt[xi] - t*eq2[xi];
    }
    return(0);
}


/***************************
* function cosPQR         *
***************************/

/* this function returns the cosine of the angle between three points P,Q,R  */
/* with Q being the center point. */

double cosPQR( double ptP[], double ptQ[], double ptR[])
{
    double QP[3];    /* vector from Q to P */
    double QR[3];    /* vector from Q to R */
    double cosine;   /* cosine of angle PQR at Q. */

    /* calculate vectors */
    QP[0] = ptP[0] - ptQ[0];
    QP[1] = ptP[1] - ptQ[1];
    QP[2] = ptP[2] - ptQ[2];
    QR[0] = ptR[0] - ptQ[0]; 
    QR[1] = ptR[1] - ptQ[1]; 
    QR[2] = ptR[2] - ptQ[2]; 

    /*calculate cosine */
    cosine = (QP[0]*QR[0]+ QP[1]*QR[1]+ QP[2]*QR[2])
        /sqrt((QP[0]*QP[0]+ QP[1]*QP[1]+ QP[2]*QP[2]) * (QR[0]*QR[0]+ QR[1]*QR[1]+ QR[2]*QR[2]));

    return(cosine);
}


/********************************
* function spherical_arc       *	   
********************************/

/* given two points in cartesian space and the center of the spherical */
/* cap between atom A and the origin, plus the radius of the sphere, */
/* this function returns the area of an arc between points B and C. */
/* the sides of the arc are the great circle distance (shortest distance) */
/* and the arc of the spherical cap centered on line AO. */

double spherical_arc(struct vertex ptAo, struct vertex ptB, struct vertex ptC, float rado)
{
    /* here, cosAOC=cosAOB, AOB = AOC. */
    /* BAC is the angle at the point on the origin, on line AO */

    double cosAOB, cosBOC; /* cosines of angles centered on (0,0,0) */
    double cosBAC, angBAC;    /* angle and cosine at vertex of cap */
    double U, V, W;
    double tansqrS, tansqrSA, tansqrSB;
    double tanprod;        /* product of tangents for spherical triangle */
    double area;           /* the value to be returned */

    cosAOB = (ptAo.xi[0]*ptB.xi[0] + ptAo.xi[1]*ptB.xi[1] + ptAo.xi[2]*ptB.xi[2])/(ptAo.dist*ptB.dist);
    cosBOC = (ptB.xi[0]*ptC.xi[0] + ptB.xi[1]*ptC.xi[1] + ptB.xi[2]*ptC.xi[2])/(ptB.dist*ptC.dist);

    U = (1+cosAOB)*sqrt((1+cosBOC)/8.0);
    V = (1-cosAOB)*sqrt((1+cosBOC)/8.0);
    W = sqrt((1-cosAOB)*(1+cosAOB)*(1-cosBOC)/8.0);  /* W == X */
    tansqrS  = (1-U+V+W+W)/(1+U-V-W-W);
    tansqrSB = (1-U-V)/(1+U+V);
    tansqrSA = (1-U+V-W-W)/(1+U-V+W+W);
    if((tansqrS*tansqrSA) > 0.0) {
        tanprod = sqrt(tansqrSB*tansqrSB*tansqrS*tansqrSA);
    } else {
        tanprod = 0.0;
    }

    cosBAC = cosPQR(ptB.xi, ptAo.xi, ptC.xi);
    if(cosBAC>1.0) { 
        angBAC = 0.0; 
    } else if(cosBAC<-1.0) {
        angBAC = PI; 
    } else { 
        angBAC = acos(cosBAC);
    } 

    /* area is area of spherical cap segment minus area of spherical triangle */
    area = rado*(angBAC*(rado-ptAo.dist) - 4.0*rado*atan(sqrt(tanprod)));
    return(area);
}

/********************************
* function assign_radii
********************************/

/* uses radii from file radii.dat. */

void assign_radii(struct atom *PDB, long int PDBtot)
{
    long int atomi;
    float  radius[11]      = { 1.61,  1.76,  1.88,  1.64,  1.64,  1.64,  1.64,  1.42,  1.46,  1.77,  1.80};

    for(atomi=0; atomi<PDBtot; ++atomi) {
        /*   OXYGEN */
        if(PDB[atomi].atomname[1] == 'O') {
            if((PDB[atomi].atomname[2] == 'G')||(PDB[atomi].atomname[2] == 'H')) {
                PDB[atomi].radius = radius[8];  /* O2H1 */
            } else { 
                PDB[atomi].radius = radius[7];  /* O2H0 */
            }

            /* SULFUR */
        } else if(PDB[atomi].atomname[1] == 'S') {
            PDB[atomi].radius = radius[9]; /* S */

            /* NITROGEN */
        } else if(PDB[atomi].atomname[1] == 'N') {
            if(!strncmp(PDB[atomi].res, "PRO", 3)) {
                PDB[atomi].radius = radius[3];  /* N3H0 */
            } else if(PDB[atomi].atomname[2] == ' ') {
                PDB[atomi].radius = radius[4];  /* N3H1 */
            } else if(PDB[atomi].atomname[2] == 'E') {
                PDB[atomi].radius = radius[4];  /* N3H1 */
            } else if((PDB[atomi].atomname[2] == 'D') && (!strncmp(PDB[atomi].res, "HIS", 3))) {
                PDB[atomi].radius = radius[4];  /* N3H1 */
            } else if((PDB[atomi].atomname[2] == 'D') && (!strncmp(PDB[atomi].res, "ASN", 3))) {
                PDB[atomi].radius = radius[5];  /* N3H2 */
            } else if((PDB[atomi].atomname[2] == 'E') && (!strncmp(PDB[atomi].res, "GLN", 3))) {
                PDB[atomi].radius = radius[5];  /* N3H2 */
            } else if((PDB[atomi].atomname[2] == 'E') && (strncmp(PDB[atomi].res, "GLN", 3))) {
                PDB[atomi].radius = radius[4];  /* N3H1 */
            } else if(PDB[atomi].atomname[2] == 'H') {
                PDB[atomi].radius = radius[5];  /* N3H2 */
            } else if(PDB[atomi].atomname[2] == 'Z') {
                PDB[atomi].radius = radius[6];  /* N4 */
            } else {
                PDB[atomi].radius = radius[5];  /* N3H2,  N default */
            }

            /* CARBON */
        } else if(PDB[atomi].atomname[1] == 'C') {
            if(PDB[atomi].atomname[2] == ' ') {
                PDB[atomi].radius = radius[0];  /* C3H0 backbone carbon */
            } else if((!strncmp(PDB[atomi].res, "ASP", 3)) && (PDB[atomi].atomname[2] == 'G')) {
                PDB[atomi].radius = radius[0];  /* C3H0 */
            } else if((!strncmp(PDB[atomi].res, "GLU", 3)) && (PDB[atomi].atomname[2] == 'D')) {
                PDB[atomi].radius = radius[0];  /* C3H0 */
            } else if((!strncmp(PDB[atomi].res, "ASN", 3)) && (PDB[atomi].atomname[2] == 'G')) { 
                PDB[atomi].radius = radius[0];  /* C3H0 */
            } else if((!strncmp(PDB[atomi].res, "GLN", 3)) && (PDB[atomi].atomname[2] == 'D')) { 
                PDB[atomi].radius = radius[0];  /* C3H0 */
            } else if((!strncmp(PDB[atomi].res, "ARG", 3)) && (PDB[atomi].atomname[2] == 'Z')) { 
                PDB[atomi].radius = radius[0];  /* C3H0 */
            } else if((!strncmp(PDB[atomi].res, "PHE", 3)) || (!strncmp(PDB[atomi].res, "HIS", 3)) 
                || (!strncmp(PDB[atomi].res, "TYR", 3)) || (!strncmp(PDB[atomi].res, "TRP", 3))) {
                    if((PDB[atomi].atomname[2] == 'A') || (PDB[atomi].atomname[2] == 'B')) {
                        PDB[atomi].radius = radius[2];  /* C4 */
                    } else {
                        PDB[atomi].radius = radius[1];  /* C3H1 */
                    }
                } else { /* carbon is C4, aliphatic */
                    PDB[atomi].radius = radius[2];    /* aliphatic carbon */
                }
        } else {
            /* default radius */
            PDB[atomi].radius = radius[10];     /* default for unknown atom; */
        }
    }
    return;
}


/******************************
* subroutine index_protein
* created 16/02/2001
*
* assigns all protein atoms to boxes within a cubic grid
* returns cube dimesions as number of boxes per side
* assigns values to global array box[]
******************************/

int index_protein(int PDBtot)
{
    int   xi;              /* coordinate counter 0,1,2 */
    int   pdbi;            /* atom counter 0 to (PDBtot-1) */
    float cellsize = CELLSIZE;  /* box dimensions */
    float maxwidth;        /* maximum protein dimension */
    int   dim;             /* dimension of DxDxD cube. */
    int   boxi;            /* box index */
    int   startind;        /* start point of box in PDBlist */


    /* ------ find global minimum and maximum ------- */
    for(xi=0; xi<3; ++xi) {
        globalmin[xi] = PDB[0].coor[xi];
        globalmax[xi] = PDB[0].coor[xi];
    }
    for (pdbi=0; pdbi<PDBtot; ++pdbi) {
        for (xi=0; xi<3; ++xi) {
            if (PDB[pdbi].coor[xi] < globalmin[xi]) {
                globalmin[xi] = PDB[pdbi].coor[xi];
            }
            if (PDB[pdbi].coor[xi] > globalmax[xi]) {
                globalmax[xi] = PDB[pdbi].coor[xi];
            }
        }
    }

    /* ------ get largest dimension of protein ------- */
    maxwidth = 0.0;
    for(xi=0; xi<3; ++xi) {
        if((globalmax[xi]-globalmin[xi]) > maxwidth) {
            maxwidth = (globalmax[xi]-globalmin[xi]);
        }
    }
    /* expand dimensions of cube by one (for round-off) */
    dim = (int)(maxwidth/cellsize) + 1;

    /* allocate memory for atomindex and PDBlist arrays */
    box = (struct atomindex *) MemNew(dim*dim*dim*sizeof(struct atomindex));
    PDBlist = (int *) MemNew(PDBtot*sizeof(int));

    if((!box) || (!PDBlist)) {
        ErrPostEx(SEV_FATAL,1,1,"Memory allocation error for box or PDBlist while calculating v-score");
    }

    /* count entries per box, assign box number to atom */
    for (pdbi=0; pdbi<PDBtot; ++pdbi) {
        boxi = (int)((PDB[pdbi].coor[0]-globalmin[0])/cellsize)*dim*dim
            + (int)((PDB[pdbi].coor[1]-globalmin[1])/cellsize)*dim
            + (int)((PDB[pdbi].coor[2]-globalmin[2])/cellsize);
        ++box[boxi].nument;
        PDB[pdbi].boxnum = boxi;
    }

    /* assign start pointers for boxes in PDBlist */
    startind = 0;
    for (boxi=0; boxi<dim*dim*dim; ++boxi) {
        box[boxi].first = startind;
        startind += box[boxi].nument;
    }

    /* clear array (needed for recounting index) */
    for(boxi=0; boxi<dim*dim*dim; ++boxi) {
        box[boxi].nument = 0;
    }

    /* fill PDBlist index */
    for (pdbi=0; pdbi<PDBtot; ++pdbi) {
        boxi = PDB[pdbi].boxnum;
        PDBlist[box[boxi].first +box[boxi].nument] = pdbi;
        ++box[boxi].nument;
    }
    return(dim);
}


/*********************************
* subroutine get_contlist4
*********************************/

/* uses box index of protein atoms to find contacts in range of a given atom. */
/* requires global variable 'box[]'. */
/* checks previous atoms, keeps only those with non-zero contact area. */

int get_contlist4(struct atom PDB[], int atomzero, struct contactlist contlist[], 
                  int PDBtot, float rado, int dim) 
{
    double sqrdist;             /* distance squared between two points */
    double neardist;            /* max distance for contact between atom spheres */
    int NC;                     /* number of contacts */
    int bai;                    /* box atom counter */
    int boxi;                   /* current box number */
    int atomj;                  /* index number of atom in PDB */
    int boxzero;                /* box atomzero is in */
    int i;                      /* dummy counter for surrounding boxes */
    long int currindex;

    NC = 0;
    /* mark previously contacted atoms */
    currindex = ca_index[atomzero];
    while(currindex != -1) {
        PDB[ca_rec[currindex].atom].done = 'C'; /* makes contact */
        currindex = ca_rec[currindex].prev;
    }

    /* get pdb atom contacts from current and adjacent boxes */
    boxzero = PDB[atomzero].boxnum;
    for(i=0; i<27; ++i) {
        /* get up to 27 boxes surrounding current box */
        boxi = boxzero +dim*dim*((i/9)-1) +dim*(((i/3)%3)-1) +(i%3)-1;
        if((boxi < 0) || (boxi >= dim*dim*dim)) continue;  /* don't do boxes outside of range */
        bai = 0;
        while(bai<box[boxi].nument) {
            atomj = PDBlist[box[boxi].first+bai]; 
            /* check previous contacts */
            if(PDB[atomj].done == 'Y') {
                ++bai;
                continue;
            }
            sqrdist = (PDB[atomzero].coor[0]-PDB[atomj].coor[0])*(PDB[atomzero].coor[0]-PDB[atomj].coor[0]) 
                + (PDB[atomzero].coor[1]-PDB[atomj].coor[1])*(PDB[atomzero].coor[1]-PDB[atomj].coor[1])
                + (PDB[atomzero].coor[2]-PDB[atomj].coor[2])*(PDB[atomzero].coor[2]-PDB[atomj].coor[2]);
            neardist =  rado + PDB[atomj].radius +Rw;
            if((sqrdist < neardist*neardist) && (sqrdist != 0.0)) {
                /* add atoms to list */
                contlist[NC].index = atomj;
                contlist[NC].dist = sqrt(sqrdist);
                ++NC;
            }
            ++bai;
        }
    }

    /* reset atoms to 'done' */
    currindex = ca_index[atomzero];
    while(currindex != -1) {
        PDB[ca_rec[currindex].atom].done = 'Y'; /* re-mark as done */
        currindex = ca_rec[currindex].prev;
    }
    return(NC);
}



/******************************
* subroutine voronoi_poly2
* created 08/07/2001  BJM
******************************/

int voronoi_poly2(struct atom PDB[], int atomzero, struct plane cont[], float rado, 
                  int NC, struct contactlist contlist[])
{
    int    cai;          /* contact atom counter */
    struct atom *ca_ptr; /* pointer to pdb atom */
    double atomdist;     /* distance to atom */
    double planedist;    /* distance to plane */
    double mindist;      /* distance to closest plane */
    int    planeA;       /* closest plane to origin */
    int    planeB;       /* second plane, with planeA defines closest edge */
    int    planeC;       /* new intersection plane for edge (endpt) */
    int    oldplaneC;    /* old intersection plane for edge (startpt) */
    double vt;           /* vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt */
    double vtmin;        /* minimum value for vt  */
    double vtdiv;        /* check for division by zero in vt calculation  */
    double temppt[3];
    struct vertex *stp; /* pointer to start vertex, coordinates */
    double *V;          /* pointer to edge vector */
    int    startedge;
    int    edgenum;
    int    vn = 0;
    char   edgeflag;
    int    edgei;      /* edge counter */
    int    vi, vj;     /* vertices counters */
    double arcpt0[3], arcpt1[3];
    int    testpA, testpB;
    double testvalA, testvalB;
    char   arcflag = 'N';

    /* failsafe variables: */
    char   recalc;       /* flag if hull is being recalculated (orig. unbounded) */
    float  origcoor[3];  /* original pdb coordinates for atom.  */

    recalc = 'N';
    /* RESTART: */
    planeA = -1;
    planeB = -1;
    planeC = -1;

    /* generate planes of contact with D = planedist */
    mindist = 9.9e+9;
    for(cai=0; cai<NC; ++cai) {
        ca_ptr = &PDB[contlist[cai].index];
        atomdist = contlist[cai].dist;

        if(planedef == 'B') {  /* bisection - original Voronoi procedure */
            planedist = atomdist/2.0;
        } else if(planedef == 'X') { /* extended radical plane (McConkey et al). */
            planedist = (atomdist*atomdist + rado*rado - (Rw + ca_ptr->radius)*(Rw + ca_ptr->radius))/(2*atomdist);
        } else { /* radical plane (Gellatly and Finney) - default. */
            planedist = (atomdist*atomdist + (rado-Rw)*(rado-Rw) - (ca_ptr->radius)*(ca_ptr->radius))/(2*atomdist);
        } 

        cont[cai].Ai[0] = (ca_ptr->coor[0] - PDB[atomzero].coor[0])/atomdist;
        cont[cai].Ai[1] = (ca_ptr->coor[1] - PDB[atomzero].coor[1])/atomdist;
        cont[cai].Ai[2] = (ca_ptr->coor[2] - PDB[atomzero].coor[2])/atomdist;
        cont[cai].Ai[3] = -planedist;
        cont[cai].dist  = fabs(planedist);
        cont[cai].index = contlist[cai].index;
        cont[cai].flag = 'X'; /* initialize contact flags to 'no contact'  */

        /* set plane0 as closest plane */
        if(cont[cai].dist < mindist) {
            mindist = cont[cai].dist;
            planeA = cai;
        }
    }

    /* add four planes surrounding atom, outer limit for voronoi polyhedron */
    cont[NC].Ai[0] = 0.707;
    cont[NC].Ai[1] = 1.0;
    cont[NC].Ai[2] = 0.0;
    cont[NC].Ai[3] = -10.0;
    cont[NC+1].Ai[0] = 0.707;
    cont[NC+1].Ai[1] = -1.0;
    cont[NC+1].Ai[2] = 0.0;
    cont[NC+1].Ai[3] = -10.0;
    cont[NC+2].Ai[0] = -0.707;
    cont[NC+2].Ai[1] = 0.0;
    cont[NC+2].Ai[2] = 1.0;
    cont[NC+2].Ai[3] = -10.0;
    cont[NC+3].Ai[0] = -0.707;
    cont[NC+3].Ai[1] = 0.0;
    cont[NC+3].Ai[2] = -1.0;
    cont[NC+3].Ai[3] = -10.0;

    /* get starting vertex from seed or calc new vertex */
    get_firstvert(cont, &planeA, &planeB, &planeC, NC, atomzero);
    solve_3x3(cont[planeA].Ai, cont[planeB].Ai, cont[planeC].Ai, temppt);

    /* add first vertex to vertex list */
    add_vertex(poly, 0, temppt, planeA, planeB, planeC);

    /* flag contacts as present */
    cont[planeA].flag = 'Y';
    cont[planeB].flag = 'Y';
    cont[planeC].flag = 'Y';

    /* calculate edge vectors */
    add_vedge(vedge, 0, cont, planeA, planeB, planeC, poly, 0);
    add_vedge(vedge, 1, cont, planeB, planeC, planeA, poly, 0);
    add_vedge(vedge, 2, cont, planeC, planeA, planeB, poly, 0);

    startedge = 0;
    edgenum = 3;
    vn = 1;

    /* --------------------------------------------------- */
    /* Generate new polyhedron points from edge vectors    */
    /* --------------------------------------------------- */

    while(1) {
        /* get next unfinished vector = startedge */
        while((vedge[startedge].endpt >= 0) && ((edgenum-startedge) > 0)) {
            ++startedge;
        }
        if((edgenum-startedge) <= 0) {
            /* all edges are done, polyhedron complete. */
            break;
        }

        vtmin = 9.9e+9; /* dummy value */
        stp = &poly[vedge[startedge].startpt];
        V = vedge[startedge].V;
        planeA = vedge[startedge].plane[0];
        planeB = vedge[startedge].plane[1];
        oldplaneC = vedge[startedge].startplane;
        planeC = -1;

        /* get closest positive intersection point */
        for(cai=0; cai<NC; ++cai) {
            /* check if contact is to be done - for now, do all. */
            if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
                vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
                if(vtdiv != 0.0) {
                    vt = -(cont[cai].Ai[0]*stp->xi[0] +cont[cai].Ai[1]*stp->xi[1] 
                    +cont[cai].Ai[2]*stp->xi[2] +cont[cai].Ai[3])/vtdiv;
                    if((vt < vtmin) && (vt > 0)) {
                        vtmin = vt;
                        planeC = cai;
                    }
                }
            }
        }
        poly[vn].xi[0] = stp->xi[0] + vtmin*V[0];
        poly[vn].xi[1] = stp->xi[1] + vtmin*V[1];
        poly[vn].xi[2] = stp->xi[2] + vtmin*V[2];

        /* if point is outside sphere, check vs. external planes */
        if((poly[vn].xi[0]*poly[vn].xi[0] + poly[vn].xi[1]*poly[vn].xi[1] 
        + poly[vn].xi[2]*poly[vn].xi[2]) > rado*rado) {
            for(cai=NC; cai<NC+4; ++cai) {
                /* check if contact is to be done - for now, do all. */
                if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
                    vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
                    if(vtdiv != 0.0) {
                        vt = -(cont[cai].Ai[0]*stp->xi[0] +cont[cai].Ai[1]*stp->xi[1] 
                        +cont[cai].Ai[2]*stp->xi[2] +cont[cai].Ai[3])/vtdiv;
                        if((vt < vtmin) && (vt > 0)) {
                            vtmin = vt;
                            planeC = cai;
                        }
                    }
                }
            }
            poly[vn].xi[0] = stp->xi[0] + vtmin*V[0];
            poly[vn].xi[1] = stp->xi[1] + vtmin*V[1];
            poly[vn].xi[2] = stp->xi[2] + vtmin*V[2];
        }

        add_vertex(poly, vn, poly[vn].xi, planeA, planeB, planeC);
        vedge[startedge].endpt = vn;
        vedge[startedge].endplane = planeC;

        /*flag contact as present */
        cont[planeC].flag = 'Y';

        /* ========  ADD EDGES  ======== */

        /* check edge (planeA, planeC) */
        edgeflag = 'Y';
        edgei = startedge+1;
        while(edgei < edgenum) {
            if(((vedge[edgei].plane[0] == planeA)&&(vedge[edgei].plane[1] == planeC)) ||
                ((vedge[edgei].plane[0] == planeC)&&(vedge[edgei].plane[1] == planeA))) {
                    /* already on list, add current vertex as endpt */
                    vedge[edgei].endpt = vn;
                    vedge[edgei].endplane = planeB;
                    edgeflag = 'N';
                    break;
                }
                ++edgei;
        }
        if(edgeflag == 'Y') { /* add edge */
            if(edgenum >= 400)
                return(-1);

            add_vedge(vedge, edgenum, cont, planeA, planeC, planeB, poly, vn);
            ++edgenum;
        }

        /* check edge (planeB, planeC) */
        edgeflag = 'Y';
        edgei = startedge+1;
        while(edgei < edgenum) {
            if(((vedge[edgei].plane[0] == planeB)&&(vedge[edgei].plane[1] == planeC)) ||
                ((vedge[edgei].plane[0] == planeC)&&(vedge[edgei].plane[1] == planeB))) {
                    /* already on list, add current vertex as endpt */
                    vedge[edgei].endpt = vn;
                    vedge[edgei].endplane = planeA;
                    edgeflag = 'N';
                    break;
                }
                ++edgei;
        }
        if(edgeflag == 'Y') { /* add edge */
            /* ===== failsafe - if solution is not converging, perturb atom  ===== */
            /* ===== coordinates and recalculate.                            ===== */
            if(edgenum >= 400) {

                /* some decoy structures have multiple atoms within van der Waals  */
                /* contact range, and may cause the hull to be unbounded. */
                /* for purposes of scoring, these atoms are given a score of zero */
                /* and further calculations omitted. */

                return(-1);

                /*
                printf("********* invalid solution for hull, recalculating *********\n");
                seed[atomzero*3] = -1;*/  /* reset to no seed vertex */
/*                origcoor[0] = PDB[atomzero].coor[0];
                origcoor[1] = PDB[atomzero].coor[1];
                origcoor[2] = PDB[atomzero].coor[2];*/
                /* perturb atom coordinates*/
/*                PDB[atomzero].coor[0] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
                PDB[atomzero].coor[1] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
                PDB[atomzero].coor[2] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
                recalc = 'Y';
                goto RESTART;
                */
            }
            add_vedge(vedge, edgenum, cont, planeB, planeC, planeA, poly, vn);
            ++edgenum;
        }
        ++vn;
    }

    /*--------------------------------------------------*/
    /*  now have voronoi polyhedron around given atom.  */
    /*  remove vertices outside of sphere, and          */
    /*  calculate intersection points with sphere.      */ 
    /*--------------------------------------------------*/

    /* flag edges that may cross sphere boundary */
    for(edgei=0; edgei<edgenum; ++edgei) {
        if((rado < poly[vedge[edgei].startpt].dist) || (rado < poly[vedge[edgei].endpt].dist)) {
            /* one or both vertices fall outside of sphere */
            arcflag = 'Y';
            vedge[edgei].arc = '?';
        } else {
            vedge[edgei].arc = 'X';
        }
    }

    /* calculate new arc points */
    for(edgei=0; edgei<edgenum; ++edgei) {
        if(vedge[edgei].arc != 'X') {
            if(solve_2xS(cont[vedge[edgei].plane[0]], cont[vedge[edgei].plane[1]], rado, arcpt0, arcpt1) == -1) {
                vedge[edgei].arc = 'X';  /* mark edge as no associated arc point */
                continue;
            }
            /* test new arc points vs. adjacent planes, add if ok. */
            testpA = vedge[edgei].startplane;
            testpB = vedge[edgei].endplane;
            testvalA = cont[testpA].Ai[0]*arcpt0[0] + cont[testpA].Ai[1]*arcpt0[1]
            + cont[testpA].Ai[2]*arcpt0[2] + cont[testpA].Ai[3];
            testvalB = cont[testpB].Ai[0]*arcpt0[0] + cont[testpB].Ai[1]*arcpt0[1] 
            + cont[testpB].Ai[2]*arcpt0[2] + cont[testpB].Ai[3];
            if((testvalA < 0.0) && (testvalB < 0.0)) { /* point is good */
                add_vertex(poly, vn, arcpt0, vedge[edgei].plane[0], vedge[edgei].plane[1], -1);
                poly[vn].dist = rado;
                ++vn;
            }
            testvalA = cont[testpA].Ai[0]*arcpt1[0] + cont[testpA].Ai[1]*arcpt1[1] 
            + cont[testpA].Ai[2]*arcpt1[2] + cont[testpA].Ai[3];
            testvalB = cont[testpB].Ai[0]*arcpt1[0] + cont[testpB].Ai[1]*arcpt1[1] 
            + cont[testpB].Ai[2]*arcpt1[2] + cont[testpB].Ai[3];
            if((testvalA < 0.0) && (testvalB < 0.0)) { /* point is good */
                add_vertex(poly, vn, arcpt1, vedge[edgei].plane[0], vedge[edgei].plane[1], -1);
                poly[vn].dist = rado;
                ++vn;
            }
        }
    }

    /* reduce poly vertex list */
    vj=0;
    for(vi=0; vi<vn; ++vi) {
        if(poly[vi].dist <= rado) {
            poly[vj] = poly[vi];
            ++vj;
        }
    }
    vn = vj;

    /* ----- calculate center points and mark engulfing contacts ----- */
    for(cai=0; cai<NC; ++cai) {
        centerpt[cai].xi[0] = -cont[cai].Ai[0]*cont[cai].Ai[3];
        centerpt[cai].xi[1] = -cont[cai].Ai[1]*cont[cai].Ai[3];
        centerpt[cai].xi[2] = -cont[cai].Ai[2]*cont[cai].Ai[3];
        centerpt[cai].dist = fabs(cont[cai].Ai[3]);
        if(cont[cai].Ai[3] > 0.0) {
            cont[cai].flag = 'E';
        }
    }
  

    if(recalc == 'N') {
        save_seeds(cont, poly, vn, atomzero);
    } else {
        /* reset atom coordinates to original values */
        PDB[atomzero].coor[0] = origcoor[0];
        PDB[atomzero].coor[1] = origcoor[1];
        PDB[atomzero].coor[2] = origcoor[2];
    }
    return(vn);
}


