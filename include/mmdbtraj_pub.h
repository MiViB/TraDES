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

#ifndef MMDBTRAJ_PUB_H
#define MMDBTRAJ_PUB_H

/* includes */
#include <bzlib.h>
#include <asn.h>
#include <objmime.h>
/* #include <objensemb.h> */
/* #include <objuturn.h> */
#include <mmdbapi.h>
/* #include <connect/ncbi_connutil.h> */

#ifdef __cplusplus
extern "C" {
#endif

/* external file names */
#define CBFNAME "cbdata"
#define BL_POTENTIAL "blpotential.txt"
#define ZHANG_POTENTIAL "zhangeij.txt"
#define ZHANG_ATOMS "zhangatm.txt"
#define EXTAA_FILENAME "SLRIextaa"


/* file extensions used throughout */
#ifdef OS_MSWIN
#define CGI_EXT ".exe"
#else
#define CGI_EXT ""
#endif
#define DB_EXT1 ".cdx"
#define DB_EXT2 ".dbf"
#define DB_EXT3 ".fpt"
#define ASN_EXT ".trj"
#define MMDB_EXT ".val"
#define PDB_EXT ".pdb"
#define BZ_EXT ".bz2"
#define LOG_EXT ".log"
#define ENS_EXT ".ens"
#define CDD_EXT ".acd"
#define CSV_EXT ".csv"
#define ARA_EXT ".ara"

#define DEFAULT_TRJ_NAME "protein"

/* global program parameters */
#define TIMEOUT prmblock.prm_timeout
#define BACKBONE_ERROR_TOLERANCE prmblock.prm_backerrtol
#define BACKBONE_PRECISION prmblock.prm_backprec
#define NUM_ROT_TRIES prmblock.prm_numrottries
#define ATOM_BOUNCINESS_BB prmblock.prm_atombouncebb
#define ATOM_BOUNCINESS_SC prmblock.prm_atombouncesc
#define BUMPCHECK_HYDROGEN prmblock.prm_bumph
#define INCSIZE prmblock.prm_incsize
#define START_BACKBONE prmblock.prm_startback
#define TRAJ3FNAME prmblock.prm_traj3fname
#define TRAJDIV prmblock.prm_trajdiv
#define WALKTYPE prmblock.prm_walk
#define TGUNITS prmblock.prm_units
#define TRAJTYPE prmblock.prm_trajtype
#define CONSTRAINT_FILE prmblock.prm_constrfile
#define MARKOV_SCALE_FACTOR prmblock.prm_markovsf
#define TUNNEL_PROB prmblock.prm_tunnel
#define USE_LOTS_RAM prmblock.prm_lotsram
#define BUILD_FRAGS_ONLY prmblock.prm_fragsonly

/* used to TGInit */
#define DB_READ 0
#define DB_CREATE 1

/* trajectory type to use for parameter TrajType */
#define TRAJ_NA 0
#define TRAJ_UNIFORM 1
#define TRAJ_STANDARD 2
#define TRAJ_SSTRU 3
#define TRAJ_GOR 4

/* specify compression type to use in various functions */
#define USE_NONE 0
#define USE_RLE 1
#define USE_BZ 2

/* supported types of random walk */
#define WALK_UNKNOWN 0
#define WALK_CA 1
#define WALK_PHIPSI 2

/* define verbosity modes */
#define VERBOSITY_SILENT (Byte)0
#define VERBOSITY_QUIET (Byte)1
#define VERBOSITY_VERBOSE (Byte)2
#define VERBOSITY_STREAM (Byte)3 


/* supported units for trajectory graph heights */
/* not actually used anywhere, just for info */
	/* USED by crease functions, CWVH */
#define UNITS_UNKNOWN 0
#define UNITS_ARBITRARY 1
#define UNITS_CREASE 2
#define UNITS_ZHANG 3
#define UNITS_BRYANT 4
#define UNITS_VORONOI 5

/* sequence alteration modes */
#define EXTAASEQ_REPLACE 1
#define EXTAASEQ_INSERT 2
#define EXTAASEQ_DELETE 3

/* portability */
#define ENDIAN_LITTLE 0
#define ENDIAN_BIG 1
#define ENDIAN_UNKNOWN 2

/* sequence encoding methods */
#define EXTAA_X 	((Int2)0)
#define EXTAA_PARENT 	((Int2)1)
#define EXTAA_ENCODE 	((Int2)2)

/* error types */
#define ERR_SUCCESS     	((TrajErr)0)
#define ERR_FAIL        	((TrajErr)1)
#define ERR_GIVEUPRES   	((TrajErr)2)
#define ERR_INCOMPLETE  	((TrajErr)3)
#define ERR_BADBACK     	((TrajErr)4)
#define ERR_DISTCONST   	((TrajErr)5)
#define ERR_CRASH       	((TrajErr)6)
#define ERR_ENDOFPROTEIN	((TrajErr)7)
#define ERR_NOPDB			((TrajErr)8)
#define ERR_CANCELLED		((TrajErr)9)

/* maketrj cis and missing residue markers */
#define VL_ALPHABET_LOWER  "abcdefghijklmnopqrstuvwxyz"
#define VL_ALPHABET_UPPER  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define VL_TRANSRESIDUE 0.0F
#define VL_CISRESIDUE 1.0F
#define VL_MISSCALPHA 2.0F

/* probability a given NOE (read in from a file) is correct */
#define NOE_PROB 0.95F
/* larger than NOE angle or dihedral possible values */
#define CONSTR_INFINITY 9999.0F

/* absolute min. and max. number of tries before backtracking */
/* irregardless of trajectory graph area */
#define BACKTRACK_TRIES_MIN 4
#define BACKTRACK_TRIES_MAX 250

/* number of bins (from -180 to 180) for phi to make histogram */
#define PHIBINS 6
/* number of bins (from -180 to 180) for psi to make histogram */
#define PSIBINS 6

/* 21 graphs collected by gettraj (#21 is cis-Pro) */
#define NUMAA 21

/* 20 is cis-Pro for gettraj */
#define GETTRAJ_CISP 20
/* for parsing primary sequence file */
#define MAXCOL 200
#define MAXRES 5000
/* longest sequence the program can accept */
#define MAXSIZE MAXRES
#define MAXGEN 250
#define GENZEROSIZE 10000
/* Generation size MUST be >1, otherwise standard deviation is infinity */
#define GENXSIZE 100

/* trajectory graph filter size */
#define TRAJFILTDIM 3

/* number of GOR secondary structure types (H, E and C) */
#define MAXSTRUCT 3	

#define MAX_KBAA_NAMELENGTH 20

/* for detecting BD-tree collisions, in angstroms */
#define PRECISION 0.00001

/* alignment defines */
#define ALIGN_ALLATOM 255
#define ALIGN_BACKBONE (AM_BACKBONE|AM_CALPHA)
#define ALIGN_CA AM_CALPHA

/* pi */
#ifndef PI
#define PI 3.14159265358979323846
#endif
#define RADTODEG (180.0/PI)
#define DEGTORAD (PI/180.0)

/* tell RotateCos which axis you wish to rotate about */
#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

/* a few atomic masses from Sargent-Welch Scientific Company P.T. */
#define AW_H 1.0079F
#define AW_C 12.011F
#define AW_N 14.0067F
#define AW_O 15.9994F
#define AW_P 30.97376F
#define AW_S 32.06F
#define AW_SE 78.96F

/* database-derived probability any given proline will be cis- */
#define P_CISP 0.0463F
/* for future use with S-S bridges */
#define P_SS 0.00F
/* average omega dihedral angle */
#define CHI_W 179.8F
#define CHISD_W 4.5F

/* used to calculate index into 1-D trajectory array using 3-D indexing */
#define TRAJ_EL(A,B,C,D) (Int4)((A)*(D)*(D)+(B)*(D)+(C))

/* Amino Acid Molecular Weight in AMU */
extern FloatLo MWaalist[];
/* All 20 amino acids */
extern CharPtr aalist;

/* to distinguish error return values from other values */
typedef Int2 TrajErr;

/* vectors */
typedef FloatLo vec[3];
typedef FloatHi vechi[3];

/* for filtering trajectory graphs */
typedef Int4 TrajFilter[(2*TRAJFILTDIM+1)*(2*TRAJFILTDIM+1)];

/* global variable parameter block */
typedef struct nlm_paramblock {
	Int4 prm_trajdiv;
	Int4 prm_timeout;
	Int2 prm_numrottries;
	Int2 prm_bumph;
	Int2 prm_walk;
	Int2 prm_units;
	Int2 prm_trajtype;
	FloatHi prm_startback;
	FloatHi prm_incsize;
	FloatHi prm_backerrtol;
	FloatHi prm_backprec;
	FloatHi prm_atombouncebb;
	FloatHi prm_atombouncesc;
	FloatHi prm_markovsf;
	Char prm_traj3fname[PATH_MAX];
	Char prm_constrfile[PATH_MAX];
	Boolean prm_lotsram;
	Boolean prm_fragsonly;
	FloatLo prm_tunnel;
} paramblock, *paramblockptr;

/* store Ramachandran plot information */
typedef struct nlm_ramastruct {
	struct nlm_ramastruct PNTR next;
	PFB pfbThis;
	FloatHi Phi;
	FloatHi Psi;
	FloatHi Omega;
	FloatHi Mag;
	FloatHi Chi1;
	FloatHi Chi2;
	FloatHi Chi3;
	FloatHi Chi4;
} RS, *PRS, **PPRS;

typedef struct slri_ExtAAStruc {
	Char keybname[MAX_KBAA_NAMELENGTH+1];
	Char descr[PATH_MAX];
	Int2 dictidx;
	Char PDBName[5];
	Char modlocation;
	Int2 nidx[3];
	Int2 cidx[3];
} EAS, *PEAS;

/* SCWRL Rotamer Library record type */
typedef struct nlm_newrotlibstruct {
  FloatLo p;
  FloatLo chi1;
  FloatLo chi2;
  FloatLo chi3;
  FloatLo chi4;
	FloatLo chi1sd;
	FloatLo chi2sd;
} newrotlibrecord, *pnewrotlibrecord;

/* conform user database record */
typedef struct nlm_tgdbstruct {
	CharPtr id;
	CharPtr basename;
	CharPtr filename;
	Int4 crtdate;
	Int4 expdate;
	struct nlm_tgdbstruct PNTR next;
} TDR, *PTDR;

typedef struct {
    Char cResName;
    Char cAtomName [5];
    Int2 iIndex;
} ZhangAtmNode, *ZhangAtmNodePtr;

typedef struct {
    PMGD pmgd;
    FloatHi potential;
	FloatHi potentialT1;  /* break out Voronoi Total Term */
	FloatHi potentialT2; /* break out Voronoi Solvent term */
	Int4 potentialCount; /*CWVH atom count in residue */
} AdjListNode, *AdjListNodePtr, *PALN;

typedef struct nlm_ambervanderwaals {
	FloatHi Rstar;
	FloatHi epsilon;
	FloatHi Rstar14;	/* for CHARMM 1-4 */
	FloatHi epsilon14;
} AVDW, *PAVDW;

typedef struct nlm_amberatomtype {
	struct nlm_amberatomtype PNTR next;
	Char PDBAtomName[5];
	Char AmberAtomType[5];
	FloatHi charge;
	Uint1 group;	/* for CHARMM non-bonded groups */
} AAT, *PAAT;

typedef struct nlm_amberimpropertype {
	struct nlm_amberimpropertype PNTR next;
	Char atm1[5];
	Char atm2[5];
	Char atm3[5];
	Char atm4[5];
} AIT, *PAIT;

typedef struct nlm_amberbond {
	struct nlm_amberbond PNTR next;
	Int2 attyp;
	FloatHi springconst;
	FloatHi length;
} AB, *PAB;

typedef struct nlm_amberangle {
	struct nlm_amberangle PNTR next;
	Int2 attyp1;
	Int2 attyp2;
	FloatHi springconst;
	FloatHi angle;
	FloatHi UBspring;	/* Urey-Bradley, if needed */
	FloatHi UBdist;
} AAN, *PAAN;

typedef struct nlm_ambertorsion{
	struct nlm_ambertorsion PNTR next;
	Int2 attyp1;
	Int2 attyp2;
	Int2 attyp3;
	FloatHi divisor;
	FloatHi barrier;
	FloatHi phase;
	FloatHi period;
} AT, *PAT;

/* AMBER energy data structure */
typedef struct nlm_AMBEREnergy {
	FloatHi EBond;
	FloatHi EAngle;
	FloatHi EDihed;
	FloatHi Evdw;
	FloatHi Evdw14;
	FloatHi Eel;
	FloatHi Eel14;
	FloatHi Esolv;
} AmberEnergy,*pAmberEnergy;

/* CHARMM solvation table */
typedef struct nlm_CHARMMSolvation {
	FloatHi volume;
	FloatHi gref;
	FloatHi gfree;
	FloatHi lambda;
} CSolvate,*PCSolvate;

/* HomTraj Loop Info */
typedef struct homtraj_loop_info {
	Int4 qstart;
	Int4 tstart;
	Int4 arrayGapStart;
	Int4 seglen;
	CharPtr qseq;
	CharPtr tseq;
	FloatHi bondEnergy;
	FloatHi totalEnergy;
	Char type;
	struct homtraj_loop_info PNTR next;
} HLI, *PHLI;

/* NOE record */
typedef struct nlm_noelike_constraintnodestruct {
  Int2 res1;
  Int2 res2;
  Char AtomName1[5];
  Char AtomName2[5];
  FloatLo MeanDist;
  FloatLo MinDelta;
  FloatLo MaxDelta;
  FloatLo Angle0;
  FloatLo Angle1;
  FloatLo Angle2;
  FloatLo Angle3;
  FloatLo Dihedral01;
  FloatLo Dihedral12;
  FloatLo Dihedral23;
  PMAD pmad1;
  PMAD pmad2;
  FloatLo prob;
  Int4 tries;
  struct nlm_noelike_constraintnodestruct PNTR next;
  struct nlm_noelike_constraintnodestruct PNTR prev;
} NN, *PNN, **PPNN;

/* fragment data structure in memory */
typedef struct nlm_fragmentdatastruct {
  Int4 resnum;
  FloatLo pSS;
  FloatLo Angle1;	/* normally phi (do NOT store cos Phi) */
  FloatLo Angle2;	/* normally psi or theta*/
  FloatLo AngleSD;	/* applied to both angles */
  FloatLo ChiWMean;
  FloatLo ChiWSD;
  Char AA;
  Int4 tout;
  Uint4 rotid;
  FloatLo prob;
  Int4 length;
  Int4 reserved;
  struct nlm_fragmentdatastruct PNTR prev;
  struct nlm_fragmentdatastruct PNTR next;
  struct nlm_fragmentdatastruct PNTR nextfrag;	/* should be null except at head of each list */
} FDS, *PFDS;

/* fragment data structure on disk */
typedef struct nlm_fragmenttruncstruct {
  Int4 resnum;
  FloatLo pSS;
  FloatLo Angle1;	/* normally phi (do NOT store cos Phi) */
  FloatLo Angle2;	/* normally psi or theta*/
  FloatLo AngleSD;	/* applied to both angles */
  FloatLo ChiWMean;
  FloatLo ChiWSD;
  Char AA;
  Int4 tout;
  Uint4 rotid;
  FloatLo prob;
  Int4 length;
  Int4 reserved;
  Char newrow;
} FDTS, *PFDTS;

/* trajectory graph structure in memory */
typedef struct nlm_trajectorystruct {
  Int2 resnum;
  Int4 *TrajGraph;
  Int4 TrajIntegral;
  Int4 Peak;
  FloatLo pCis;
  Int4 *CisTrajGraph;
  Int4 CisTrajIntegral;
  Int4 CisPeak;
  FloatLo pSS;
  FloatLo ChiWMean;
  FloatLo ChiWSD;
  PFDS pfdsFragmentHead;
  Char AA;
  Int2 dim;
  Int2 firstnzrow;
  Int2 numnzrows;
  Int4 nelt0;
  Int4 nelt5;
  Int4 nelt10;
  Int4 nelt15;
  Int2 tout;
  FloatLo markovsf;
  Uint4 rotid;
  Boolean dumpcsv;
} TGS, *PTGS, **PPTGS;

/* dotplot node */
typedef struct nlm_dotstruct {
	struct nlm_dotstruct PNTR nextcol;
	struct nlm_dotstruct PNTR nextrow;
	vec dist;
	Int2 resID;
	PALD paldrow;
	PALD paldcol;
} DS, *PDS, **PPDS;

/* node of a bd-tree */
typedef struct nlm_BDNode {
	PALD paldThis;
	vec dim;
	struct nlm_BDNode PNTR way[6];
} BDNode, *PBDNode;

/* this type is used to define a world, a set of linked list of sets of
   molecules or other PFBs, each in a given model number and MSD */
typedef struct nlm_worldstruct {
	struct nlm_worldstruct PNTR next;
	Int2 Model;
	PFB pfbThis;
	PBDNode pbdTreeHead;
} WS, *PWS, **PPWS;

typedef struct nlm_alignedrangenode {
	struct nlm_alignedrangenode PNTR next;
	Int4 start;
	Int4 end;
	FloatLo score;
} ARN, *ARNP;

typedef struct nlm_alignedregions {
	struct nlm_alignedregions PNTR next;
	CharPtr seqname;
	ARNP ranges;
} AR, *ARP;

typedef struct nlm_MakeTrjParamBlock {
	/* from structure or from sequence file */
	Int2 TrajMethod;        /* 1=from structure, 2=from sequence file */

	/* if from structure file (TrajMethod=1), these are used */
	CharPtr valfnamin;      /* structure filename (no extension) */
	CharPtr pcChain;        /* chain in structure, default=first found */
	Int2 modelnum;          /* Model number in structure, normally 1 */
	Int2 startres;          /* start residue */
	Int2 endres;            /* end residue (put 0 for last residue) */
	Int4 peakheight;        /* maximum peak height */
	Int2 noise;             /* percent noise to add - normally zero */
	Int2 units;		/* units of graphs */
	Int2 savechi;	/* if TRUE, rotamer chis will be preserved too */
	/* peak uncertainty given either as standard deviations */	
	FloatLo sigma_x;        /* s.d. in x direction (degrees) */
	FloatLo sigma_y;        /* s.d. in y direction (degrees) */
	/* or calculated from temperature and timestep */
	FloatLo temperature;    /* Temperature (degrees K) */
	FloatLo timestep;       /* timestep (fs) */
	
	/* if sequence file, these are used */
	CharPtr seqfnamin;      /* Sequence filename */
	Int2 tgtype;            /* 1=Uniform, 2=AA-based, 3=1-state GOR, 4=3-state GOR */
	CharPtr sstrufname;     /* sec. structure prediciton output file, or if SStruinput
	                           is TRUE, name of sec. structure prediciton input file */
	Boolean DoRPS;          /* call RPSBlast to get alignments? */
	Boolean SStruinput;     /* if FALSE, sec. structure prediction is done with GOR
	                           if TRUE, it is read in from file "sstrufname"; note if
	                           tgtype<3 this setting and sstrufname are irrelevant */
	Boolean uturn;		/* do U-turn predictions and add fragments for them */

	/* these options are common to all cases */
	CharPtr trjfnamout;     /* output .trj name (no extension), same as tmpfnamin if
	                           not given and tmpfnamin has a value */
	Int2 comprtype;         /* compression scheme to use: 0=None, 1=RLE, 2=Bzip2 */
	CharPtr constrfnamin;   /* Constraint filename to read in, if any */
	CharPtr templat;		/* template structure to use */
	CharPtr alignfile;      /* a file containing an alignment between query and template sequences */
	Int2 zhangwindsize;		/* size of exclusize window for ZD energy calculation 
								when computing s.d. from temp. and timestep */
	CharPtr ssmask;			/* DSSP-style string - if given, will try to keep helices
								from moving */
	CharPtr errorfile;		/* filename to print homtraj errors to */
	Boolean shiftgaps;               /* Shift alignment gaps to optimal locations */
	CharPtr fragfile;		/* filename containing fragments to import and add to nascent trajectory distribution */
	Boolean dumpcsv;        /* dump a .csv file with the same name as the .trj file */
	Boolean benchmark;      /* don't dump ARA when benchmarking */	
    Boolean all_coil;     /* make all-coil sampling 3 state = 100% coil */
    Boolean all_beta;     /* make all-beta sampling 3 state = 100% extended */
	                     /* both all_coil + all_beta true = 50% coil, 50% extended */
} MakeTrjParamBlock, *pMakeTrjParamBlock;

typedef struct nlm_FoldTrajParamBlock {
	PMSD pmsdRoot;		/* struc to fold */
	Int2 Model;		/* model number to fold */
	TrajErr err;		/* filled in by foldtraj */
	Int2 gen;
	CharPtr errorfile;	/* file to print homtraj errors to */
} FoldTrajParamBlock, *pFoldTrajParamBlock;

typedef struct nlm_CharmmParamBlock {
	PMMD pmmd;				/* molecule to minimize */
	Int2 Model;				/* model number to minimize */
	Int2 numsteps;			/* maximum number of steps */
	Boolean ShowProgress;	/* verbose output? */
	FloatHi softatompad;	/* value to add to vdw radius for soft atoms */
} CharmmParamBlock, *pCharmmParamBlock;

/* dotplot node */
typedef struct nlm_hbondstruc {
	struct nlm_hbondstruc PNTR next;
	struct nlm_hbondstruc PNTR prev;
	PMAD pmadDonor;
	PMAD pmadAcceptor;
	Boolean certain;
} HBS, *PHBS, **PPHBS;

typedef struct nlm_charmmnbxstruc {
	struct nlm_charmmnbxstruc PNTR next;
	Char atom1[5];
	Char atom2[5];
} CNS, *PCNS, **PPCNS;

extern Char CFG_local_datafilepath[PATH_MAX];
/* for progress monitors */
extern Int4 volatile ProgramProgress;
extern Int4 volatile ProgramProgressMax;

/* create worlds for BD-trees */
PWS AddtoWorldEx PROTO((PWS pwsCurrent,Int2 Model,PFB pfbThis,Boolean AllConfs));
PWS AddtoWorld PROTO((PWS pwsCurrent,Int2 Model,PFB pfbThis));
ValNodePtr InstantiateWorld PROTO((Int2 Id,PWS pwsNew));
ValNodePtr FreeWorld PROTO((ValNodePtr vnpKill));
ValNodePtr FreeAllWorlds PROTO((void));
PWS FindWorld PROTO((Int2 ID));
/* return a list of all atoms with given radius of centre */
ValNodePtr FindAtomsIn PROTO((PWS pwsWorld,vec centre,FloatLo radius));
/* insert co-ordinates of all atoms in pfbThis to pbdTreeHead bd-tree;
   id is only used for error messages and is safe to set to zero */
void AddToBDTree PROTO((PFB pfbThis,Int2 Model,PBDNode *ppbdTreeHead));
/* remove a node from pbdTreeHead having vCoOrd as co-ordinates */
void BDRemove PROTO((vec vCoOrd,PBDNode *ppbdTreeHead,PHBS *pphbsHbond));
void FreeBDTree PROTO((PBDNode pbdTreeHead));
void TraverseBDTree(PBDNode pbdTreeHead, Int4 level);

/* clears bReserved bit status in given modelnum for given pdnmmThis */
void ClearUpdate(PDNMM pdnmmThis,Int2 modelnum);
/* sets bReserved */
void SetUpdate(PDNMM pdnmmThis,Int2 modelnum);

/* vector function declarations */
/* res=a x b */
void Cross PROTO((vec res,vec a,vec b));
/* returns inner product of a and b */
FloatLo Dot PROTO((vec a,vec b));
/* magnitude of a */
FloatLo getMag PROTO((vec a));
/* res=a/|a| */
void Normalize PROTO((vec res,vec a));
/* res=a + b */
void VecAdd PROTO((vec res,vec a,vec b));
/* res=a - b */
void VecSub PROTO((vec res,vec a,vec b));
/* res=scale*a where scale is a scalar, not a vector */
void VecScale PROTO((vec res, vec a,FloatLo scale));
/* res=-a */
void NegateVec PROTO((vec res, vec a));
/* print out vector */
void PrintVec PROTO((vec a));
/* a=a+trans */
void Translate PROTO((vec a,vec trans));
/* res becomes a, rotated about the given axis (see #defines above) by angle with
   sine sina and cosine cosa */
void Rotatecos PROTO((Int2 axis, FloatLo sina, FloatLo cosa, vec a,vec res));
/* given 4 consecutive co-ordinates (vectors), returns their dihedral angle, the two normal
   angles they form and the three bond lengths they form; rng is simply used such that
   all angles (in degrees) have return values in the interval [rng, rng+360) */
void GetDihedral PROTO((vec v1,vec v2,vec v3,vec v4,FloatLo rng,FloatLo *dihed,FloatLo *ba1,FloatLo *ba2,FloatLo *bl1,FloatLo *bl2,FloatLo *bl3));

void Crosshi PROTO((vechi res,vechi a,vechi b));
FloatHi Dothi PROTO((vechi a,vechi b));
FloatHi getMaghi PROTO((vechi a));
void Normalizehi PROTO((vechi res,vechi a));
void VecAddhi PROTO((vechi res,vechi a,vechi b));
void VecSubhi PROTO((vechi res,vechi a,vechi b));
void VecScalehi PROTO((vechi res, vechi a,FloatHi scale));
void NegateVechi PROTO((vechi res, vechi a));
void PrintVechi PROTO((vechi a));
void Translatehi PROTO((vechi a,vechi trans));
void GetDihedralhi PROTO((vechi v1,vechi v2,vechi v3,vechi v4,FloatHi rng,FloatHi *dihed,FloatHi *ba1,FloatHi *ba2,FloatHi *bl1,FloatHi *bl2,FloatHi *bl3));
void VechitoVec(vec res,vechi a);
void VectoVechi(vechi res,vec a);

/* various dotplot function declarations */
PDS DotPlot PROTO((PWS pwsWorld,Int2 NOEelement));
PDS DotPlotSparse PROTO((PWS pwsWorld,FloatLo maxRadius));
PDS GetCAlphaTrace PROTO((PMMD pmmdMol,Int2 Model));
PDS freeDS PROTO((PDS pdsHead));

/* low-level Ramachandran plot function for a pre-loaded Modelstruc, works
   on a single residue/atom */
TrajErr getR PROTO((PMAD curpmad,PRS PNTR prsHead,Int2 Model));
void GetTrace PROTO((PRS prsLast[][PHIBINS][PSIBINS],PRS prsHead[][PHIBINS][PSIBINS],PMMD pmmdHere,Int2 Model,Char pdbid[5]));
PRS Rama PROTO((ValNodePtr vnpHere,Int2 Model));
PRS freeRS PROTO((PRS prsHead));

/* CodeBase-related trajectory graph functions */
/* allocate and free trajectory graph structures */
PTGS NewTraj(Int2 cis,Int2 graphwidth);
PTGS FreeTraj(PTGS ptgsThis);
/* robust tmpnam for DFP */
CharPtr DFPTmpNam(Boolean nodots,Boolean adddbext);
/* pack/unpack ASN.1 trajectory database */
CharPtr UnPackAsnTrajGraph(CharPtr ASNfnam,Int2 PNTR seqlength,CharPtr seq,CharPtr DumpHeaderName,BiostrucIdPtr *ppbi,BiostrucDescrPtr *ppbd,ValNodePtr *bioseq);
void PackAsnTrajGraph(CharPtr fnam,CharPtr sequence,ByteStorePtr id,ByteStorePtr descr,ValNodePtr bioseq,Boolean bzipit);
/* erase database files from disk */
void CleanUpDB(CharPtr fnam);
/* open trajectory graph database for read/write access; if dbfnam does not exist, a
   new database will be created if makeit=DB_CREATE, or an error will be given if
   makeit=DB_READ */
TrajErr TGInit(CharPtr dbfnam,Int2 makeit,Int4 PNTR numrec);
/* shut down and free memory used by database */
void TGClose(void);
void TGPack(void);
TrajErr BuildEmptyTG(CharPtr seq);
/* read a trajectory graph from the database into a structure buffer */
PTGS TrajGraphRead(Int2 resnum);
/* write a trajectory graph record to the end of the database */
TrajErr TrajGraphWrite(PTGS ptgsThis,Int4 ctype,Boolean Replace);
/* integrate a trajectory graph in memory only */
void TrajGraphIntegrate(PTGS ptgsHere);
/* differentiate a trajectory graph in memory only */
void TrajGraphDifferentiate(PTGS ptgsHere);
/* apply a filter to trajectory graph */
void TrajFilt(PTGS ptgsHere,TrajFilter filt);
void TrajScale(PTGS ptgsHere,FloatLo ScaleFactor);
/* fill in sparsity values nelt0-nelt6 for a given TG */
void TrajCalcSparsity(PTGS ptgsHere,Int2 cis, FILE *ara);
/* fill in nz elements in a given PTGS */
void TrajCalcNZ(PTGS ptgsHere);
TrajErr TrajAdd(PTGS ptgsDest,PTGS ptgsAdd);
void TrajCalcTout(PTGS ptgsHere);
/* change pCis value of a residue directly in the unpacked temporary database */
void UpdatePCis(Int4 resnum,FloatLo cis);
FloatLo GetPCis(Int4 resnum);

/* check if bond is disulphide */
TrajErr IsBondSulfide(PMGD pmgd,Int4 *Sarray,Int2 AddRestraint);
PMAD IsDisulfideS(PMAD pmadHere);

/* get first useful model in a structure file */
Int2 GetFirstFullModel(PMSD pmsdHere);
PMLD GetModelN(PMSD pmsdHere,Int2 n);

/* Homology Modeling */
/*TrajErr HomTraj(CharPtr extsequence,Int2 numpass,FloatLo stddev,CharPtr fnam,CharPtr templat,CharPtr alignfile,CharPtr ssfile,CharPtr errfnam);*/
TrajErr HomTraj(CharPtr fnamout,CharPtr alignfile,CharPtr ssfile,CharPtr errfnam, Boolean shiftgaps);

/******************** GOR ******************/
/* Load a GOR formatted database from fname file to 3 arrays of size pnprot, with pseq sequences, ptitle deflines and pires residues */
TrajErr LoadGOR_DB(CharPtr fname, Int4 PNTR pnprot, CharPtr **pseq, CharPtr **ptitle, Int4Ptr *pires);

/* Free the contents from LoadGOR_DB, and set variables to NULL */
void FreeGOR_DB(Int4 nprot, CharPtr **pseq, CharPtr **ptitle, Int4Ptr *pires);

/* Perform Single or Double Jackknife statistical analysis using the Enhanced GOR prediction method */
TrajErr EGORDoubleJacknife(Int4 dbID, Boolean bDoubleJK, Boolean bUseStrDB, FILE *fpdetailedresults);

/* given a primary sequence szInputSeq, this returns the GOR secondary structure predictions,
   with arwFiltStruct containing the one-state predictions at each residue (H, E or C) and
   arwProMatrix the 3-state predictions (probabilities of H, E and C at each residue) */
TrajErr fnGilGOR(Char *szInputSeq,FloatLo arwProbMatrix[MAXRES][MAXSTRUCT],Char arwFiltStruct[MAXRES]);
/********************************************/

/* Zhang potential creating/computing/freeing */
TrajErr LoadZhangPotential (Int4Ptr *ppiZhangPtnl, FILE *fp);
void FreeZhangPotential (Int4Ptr *ppiZhangPtnl);
Int2 LoadZhangAtmList (DValNodePtr *ppdnResult, FILE *fp);
void FreeZhangAtmList (DValNodePtr *ppdnHead);
void ComputeZhangPotential (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds,
                            Int2 iModelNum, Boolean bInclusiveWindow,
                            Int2 iWindowSize, Int4Ptr piZhangPtnl,
                            DValNodePtr pdnAtmList, Boolean bDetailedList);
void LIBCALLBACK FreeZhangAtmNode (Pointer pzan);
void FreeAdjList (DValNodePtr *ppdnHead);

/* Vscore potential computing */
TrajErr ComputeVscorePotential (PMMD pmmdToScore, DValNodePtr *ppdnResult, Int2 iWindowSize);
FloatHi GetTotalVPotential();
FloatHi GetTotTermVPotential();
FloatHi GetSolvTermVPotential();
FloatHi GetSSTermVPotential();
FloatHi GetASATotal();
Int4 GetPDBatomsVPotential();

/*******************************************************
* Enumerated type for different 3D potentials derived
* as a cross between BL and secondary structure contacts
* used in AddtoBLPotential in Kaca's experimentation
* with different bins
********************************************************/
typedef enum {POT_BL = 1, POT_4D = 2, POT_SS = 3, POT_SB = 4, POT_SB2 = 5} POT_TYPE;
void NewContactCount(POT_TYPE pot_type, Int4Ptr *ppiCount, Int4 *piCount);
void FreeContactCount(Int4Ptr *ppiCount);
Int4 CountStructureContacts (Int4Ptr *ppiCount, DValNodePtr pdnListPmmds, Int2 iModelNum,
		Boolean bInclusiveWindow, Int2 iWindowSize, Boolean bUsingPep, Boolean bModelling, Boolean bModel, POT_TYPE pot_type);


		
/******************************
* Bryant Potential functions  *
******************************/
Int4 NewBLPotential (Int4Ptr *ppiBryantPotential);
void FreeBryantPotential (Int4Ptr *ppiBryantPotential);

/* Counts contacts for a structure or a modelled structure with residues in the bReserved position of the pgmd */
Int4 AddtoBLPotential (Int4Ptr *ppiBryantPotential, DValNodePtr pdnListPmmds,
					Int2 iModelNum, Boolean bInclusiveWindow,
					Int2 iWindowSize, Boolean bUsingPep,
					Boolean bModelling, Boolean bModel);


/* Dump either the BL Table of counts or the potential itself */
/* The first two will write to the file with the given pcID as filename */
TrajErr DumpBLTable(Int4Ptr piBryantPotential, CharPtr pcID, Boolean bUsingPep);
TrajErr LoadBryantTableByName(Int4Ptr PNTR ppiBryantPotential, CharPtr filename,Boolean bUsingPep, Boolean bNegateInput);

TrajErr DumpBLPotential (Int4Ptr piBryantPotential, FloatLoPtr pfHydrophobicity, CharPtr pcID, Boolean bUsingPep);


/* Load either the BL Table of counts or the potential itself */
TrajErr LoadBryantTable (Int4Ptr *ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput);
TrajErr LoadBryantTableByName(Int4Ptr PNTR ppiBryantPotential, CharPtr filename,
                          Boolean bUsingPep, Boolean bNegateInput);

TrajErr LoadBryantPotential (Int4Ptr *ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput);
TrajErr LoadBryantPotential_AlphaIndex (Int4Ptr PNTR ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput);



void ComputeBryantPotential (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds,
                             Int2 iModelNum, Boolean bInclusiveWindow,
                             Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList);

/* The next two functions are for use with sequence to structure models */
void ComputeBryantPotentialModel (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds, Int2 iModelNum,
							 Boolean bInclusiveWindow, Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList);
void ComputeBryantPotentialTemplate (DValNodePtr *ppdnResult, DValNodePtr pdnListPmmds, Int2 iModelNum,
							 Boolean bInclusiveWindow, Int2 iWindowSize, Int4Ptr piBryantPotential,
                             Boolean bUsingPep, Boolean bDetailedList);

void PrintAdjList (DValNodePtr pdnHead, Int2 iModelNum, Boolean bVerbose,
                   Boolean bTemp);
FloatHi GetTotalPotential (void);
TrajErr MinimizeBLPotential(PMMD pmmd,Int2 Model,Int2 numstep,Int2 tsize);

/********************************************************
* 4D Potential functions (taken out from BL potential)  *
********************************************************/
Int4 AddtoBLPotentialEx(Int4Ptr *ppiBryantPotential, DValNodePtr pdnListPmmds,
					Int2 iModelNum, Boolean bInclusiveWindow,
					Int2 iWindowSize, Boolean bUsingPep,
					Boolean bModelling, Boolean bModel, POT_TYPE pot_type);
Int4 New4DPotential (Int4Ptr *ppiBryantPotential);
TrajErr Dump4DTable(Int4Ptr piBryantPotential, CharPtr pcID, Boolean bUsingPep);
TrajErr Dump4DPotential (Int4Ptr piBryantPotential, FloatLoPtr pfHydrophobicity, CharPtr pcID, Boolean bUsingPep);
TrajErr Load4DPotentialEx (Int4Ptr *ppiBryantPotential, FILE *fp, Boolean bPotential,
                          Boolean bUsingPep, Boolean bNegateInput, Boolean bAlphaIndex);
TrajErr Load4DTable (Int4Ptr PNTR ppiBryantPotential, FILE *fp,
                          Boolean bUsingPep, Boolean bNegateInput);
TrajErr Load4DTableByName(Int4Ptr PNTR ppiBryantPotential, CharPtr filename,
                          Boolean bUsingPep, Boolean bNegateInput);

/* used by potential.c */
DValNodePtr CreateListPmmds (CharPtr chainList, PMSD pmsd);
void FreeListPmmds (DValNodePtr *ppdnHead);
ValNodePtr CreateListOfModels (CharPtr pcModelList, PMSD pmsd);
void FreeListOfModels (ValNodePtr *ppvnHead);
Boolean IsValidModel (PMSD pmsd, Int2 iModelNum);
void PrintResNumbers (DValNodePtr pdnListPmmds);
void PrintSecStrucInfo (DValNodePtr pdnListPmmds);

Int2 GetIndex (CharPtr str); /* Gets a numeric index according to the hydrophobic arrangement of letters */
Int2 GetAlphaIndex (CharPtr str);/* Gets a numeric index according to the alphabetical arrangement of letters */ 

/* potentials are stored in a 20 x 20 x 6 static array
                              aa x aa x bin
   These r14 potentials are supplied in a file */
/* r14init: load the r14 potentials to memory */
Int2 r14init(void);
/* r14Tstar: assign the E*(i) Turn Propensity Potential to each residue in the sequence. */
Int2 r14Tstar(Int2 AA0i, Int2 AA1i, Int2 AA2i, Int2 AA3i, FloatLo *Tstar_i);
Int2 r14turncall(ValNodePtr Evnp, ValNodePtr Efiltvnp, ValNodePtr PNTR turnListvnpp);
Int2 r14bscalc(ByteStorePtr seqbsp, FloatHi cutoff,
                    ValNodePtr PNTR vnppTD, ValNodePtr PNTR vnppFT,
                    ValNodePtr PNTR vnppFiltFT, ValNodePtr PNTR vnppFiltTD);

/*
TrajErr UTDBInit(Boolean count);
void UTDBClose(void);
void UTDBPack(void);
TrajErr WriteUturnDBRec(UturnStrucPtr usp,Int4 uid);
TrajErr ReadUturnDBRec(UturnStrucPtr *pusp,Int4 uid);
TrajErr FindSuitableUturn(UturnStrucPtr *pusp,Int4 seqlen,PNN pnn,FloatLoPtr dist);
PFDS BuildFragmentListFromUSP(UturnStrucPtr usp,Int2 resnum,CharPtr sequence,FloatLo prob);
*/

/* crease energy function */
TrajErr CalcCreaseEnergy(PMMD pmmd,Int2 Units,Int2 cutoff,Int2 ExcludeWindow,Boolean Decay,Boolean Color,Int2 ModelNum);
void FreeCreaseEnergy(void);
void PrintCreaseEnergy(void);
FloatLo GetCreaseEnergy(void);
/* give crease energy at a particular residue */
FloatLo GetCreaseRes(Int4 res);
void ReplaceTempWithPotential (DValNodePtr pdnAdjList, Int2 iModelNum);

/* build residue side chain */
void Buildit PROTO((PMGD pmgdThis, vec vCAlpha, Char ResID, FloatLo chi1, FloatLo chi2, Byte rotnum));

/* functions opening/makeing data files */
/* creates a valid ASCII ASN.1 file (named fnam) containing the a-carbon skeleton for the
   protein sequence described by seq */
void BuildSkelASN PROTO((CharPtr seq,CharPtr fnam));
void BuildSkelASNDict(Int2 numrot,CharPtr fnam);
/* open/close rotamer library */
TrajErr LoadRotLib(void);
void FreeRotLib(void);
/* read in parameters in config file */
void LoadParams(void);
/* read in constraint files */
TrajErr LoadDistConstraints(CharPtr fnam,CharPtr sequence);

PMSD LoadABiostruc(CharPtr fnam,Boolean remote,Int2 mdllvl,Int2 *reqmodel);
PMSD LoadABiostrucEx(CharPtr fnam,Boolean remote,Int2 mdllvl,Int2 *reqmodel,ValNodePtr *ppvnBioseq);
BiostrucPtr MIMEBiostrucAsnGet(CharPtr fnam,const CharPtr mode,ValNodePtr *BiostrucKeep);
/* MIME structure builder */
NcbiMimeAsn1Ptr BuildMIMEBiostruc(BiostrucPtr bsp,CharPtr seq,ValNodePtr sequences);

VoidPtr TRADEProtein(VoidPtr ptr, Int4 StructureNumber);
Int4 FillTG(Int2 TrajType,CharPtr cres,Char SStru,FloatLo prob[3],Int2 trajwidth,Int4 PNTR TrajGraph,Boolean UseCis);
TrajErr FillInTGS(PTGS ptgsHere,PDNMG pdnmg,PRS prsHere,PRS prsLast, Int2 start, FloatLo sigma_x, FloatLo sigma_y,Int4 mag,Int4 *Sarray, Int2 percent);

/* extended sequence conversion functions */
Int2 ConvertExtSequence(CharPtr src,CharPtr dest,Int2 extend);
CharPtr DecodeSequence(CharPtr seq,Int2 extend);
/* returns 0 if ch is not a valid 1-letter amino acid code */
Int2 isAA PROTO((Char ch));
Char GetAAFromIDict (PMGD pmgd);
Int2 GetResIdxFromMGD(CharPtr aalist,PMGD pmgdHere);
Int2 WildCardMatch(CharPtr test,CharPtr pattern);
Int4 GetResIdxFromSeq(CharPtr encseq, Int4 iresnum);

/* RMSD alignment functions */
FloatLo GetRMSD(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model1,Int2 Model2);
FloatLo GetRMSDEx(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 resa1,Int2 resb1,Int2 resa2,Int2 resb2,Int2 Model1,Int2 Model2);
Int2 PutMolBackWhereItCameFrom(PMMD pmmdHere,Int2 Model);
Int2 Align2StrucSVD(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model1,Int2 Model2);
Int2 Align2StrucSVDEx(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 resa1,Int2 resb1,Int2 resa2,Int2 resb2,Int2 Model1,Int2 Model2);
Int2 AlignMultipleStructures(PDNMM pdnmmHead1,PDNMM pdnmmHead2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model,Int4 offset1,Int4 offset2);
Int2 AlignMultipleStructuresEx(PDNMM pdnmmHead1,PDNMM pdnmmHead2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model,Int4 offset1,Int4 offset2, FILE *fp);
Int2 AlignStructures(PMMD pmmd1,PMMD pmmd2,Int2 atomtypes,Int2 res1,Int2 res2,Int2 Model1,Int2 Model2);
void RotateMolecule PROTO((PMBD pmbdAxis,PMAD pmadDir,FloatLo Theta,Int2 Model));
PNN ComputeTakeoffAngles(PMMD pmmd,Int2 Model,Int4 res1,Int4 res2,Int2 walk);
vec *ComputeTakeoffCoords(PMMD pmmd,Int2 Model,Int4 res1,Int4 res2,Int2 walk);

/* get radius of gyration */
FloatLo GetRgyr(PMMD pmmdRoot,Int2 Model);
FloatLo GetHPRgyr(PMMD pmmdRoot,Int2 Model);
FloatLo GetExactRgyr(PMMD pmmdRoot,Int2 Model,Boolean HPonly);

/* extended MMDBAPI opening/closing */
Int2 OPENMMDBAPI(Byte bExtent,CharPtr pcDictName);
void CLOSEMMDBAPI(void);

/* residue/atom co-ordinate related functions */
PALD GetCoOrds(PMGD pmgdHere,CharPtr AtomName,vec vCAlpha,vec vDest,Int2 Model);
void CalcNextCoOrd PROTO((vec vA,vec vB,vec vC,vec vD,FloatLo MagCD,FloatLo Beta,FloatLo Chi));
PRS GetTrjAngle(PMMD pmmd, Int2 DoChis, Int2 Model, FloatLo tmp, DValNodePtr pdnIncList);
Boolean IsNTermModified(PMGD pmgdThis);
Boolean IsCTermModified(PMGD pmgdThis);
PMAD FindCAlpha PROTO((PVNMA pvnmaHead));
PMAD FindCBeta PROTO((PVNMA pvnmaHead));
PMAD FindCBackbone (PVNMA pvnmaHead);
PMAD FindNBackbone (PVNMA pvnmaHead);
PMAD findNeighbour PROTO((PMAD curAtom,Uint1 anum,Boolean Backbone));
PMAD FindAtomName PROTO((PVNMA pvnmaHead,CharPtr AtomName));
/* convert pfbfrom into node type specified in toNodeType, which must be one of:
   AM_MSD, AM_MMD, AM_MGD, AM_MAD or AM_MBD; when traversing up the hierarchy,
   this behaves much like a GetParent function, and when traversing down,
   objects are expanded to all sub-objects */
ValNodePtr ConvertNode PROTO((PFB pfbfrom,Byte toNodeType));
ValNodePtr BackupCoords(PMMD pmmdHere,Int2 Model);
void RestoreCoords(PMMD pmmdHere,Int2 Model,ValNodePtr vnpCoords);
Int4 CalcExtendedResidues(PMMD pmmdHere,Int2 Model);

/* returns a bounding sphere for an entire World */
FloatLo BoundSphereWorld PROTO((PWS pwsWorld, vec centre));
/* returns bounding box for entire World */
Int2 BoundBoxWorld PROTO((PWS pwsWorld, vec *finalbox));
/* bounding box function declarations */
/* find bounding box */
Int2 BoundBox PROTO((PFB pfbThis, Int2 Model, vec *result));
/* find a bounding sphere, return value is radius and centre is returned in centre;
   note this may not always be the smallest possible bounding sphere */
FloatLo BoundSphere PROTO((PFB pfbThis, Int2 Model, vec centre));

/* printf to status file used by conform */
void StatusPrintf(CharPtr format, ...);

/* global bail procedure, does cleanup */
void PurgeGlobs(void);

/* distance constraint functions */
void InitPNNList(PMMD pmmdHere);
PNN FreePNNList(void);
/* extra distance constraint for S-S */
TrajErr AddDisulfideRestraint(Int4 a,Int4 b);
FloatLo InDisulphide(Int2 resnum);

/* TGDB related functions */
void ClearExpiredStatus(void);
void AddStatusTGDB(void);
Int4 MyWWWSendFile(CharPtr Name,FILE *File);
TrajErr ValidateUser(Int4 uid);

/* sidechain angle functions */
TrajErr get_rotamers(PMGD pmgdHere,FloatLo phi,FloatLo psi,FloatLo arandnum,FloatLo PNTR chi1,FloatLo PNTR chi2,FloatLo PNTR chi3,FloatLo PNTR chi4,FloatLo PNTR chi1sd,FloatLo PNTR chi2sd);
ValNodePtr MeasureChis(PMMD pmmd,Int2 Model);
TrajErr MeasureResidueChis(PMGD pmgd,PRS *pprshere,Int2 Model,Boolean onlyburied,FloatLo temp,DValNodePtr dvnpPotential);
ValNodePtr FreeChis(ValNodePtr vnp);
/* convert between PTGS->rotid and PRS->Chi# */
void GetChiFromRotid(PRS *pprsHere,Uint4 rotid);
Uint4 ComputeRotid(PRS prsHere);
/* hydrogen rebuilding function */
TrajErr RebuildH(PMGD pmgdThis,Int2 Model);

/* chooses a random number from -1 to 1 */
FloatLo Rand1 PROTO((void));
FloatLo Rand1Distrib(void);
FloatLo RandWeibull(FloatLo a, FloatLo b);
 
/* get endianness of machine */
Int2 GetEndian(void);

/* general purpose functions */
/* remove all but alphanumeric and _ characters from a filename */
void FilterFnam(CharPtr fnam);
/* like ANSI fgets function but for a ByteStore */
Int2 BSgets(CharPtr buf,Int2 maxsz,ByteStorePtr bs);
/* remove extra padding from database strings */
Int4 StripEndSpace(CharPtr fnam);
/* test file MD5 checksum */
Boolean CheckMD5(CharPtr fnam,CharPtr cksum);

PEAS GetExtAAInfo(CharPtr buf);
TrajErr SubstituteNames(Char nam[5],Int4 DictIdx);

/* distance constraint functions */
TrajErr AddDistConstraint(PNN pnnNode);
void PrintDistConstraints(void);
void PrintTrueDistConstraints(void);
void PrintDistConstTries(void);
PNN GetTrueDistConstraints(void);
PNN DeleteDistConstraint(PNN pnnNode);

/* fragment related functions */
PFDS FreeFragmentList(PFDS pfdsHead);
PFDS AddFragmentResidue(PFDS *pfdsHead,Int4 res,Char aa,FloatLo a1,FloatLo a2,FloatLo asd,FloatLo chiw,FloatLo chiwsd,Uint4 rotid,FloatLo p);
TrajErr AddFragmentList(PTGS ptgsHere,PFDS pfdsHead);
TrajErr RemoveFragmentList(PTGS ptgsHere,PFDS pfdsHead);
TrajErr RemoveFragmentResidue(PTGS ptgsHere,PFDS pfdsHead,Int4 resnum);
FloatLo GetTotalFragmentProb(PTGS ptgsHere);
TrajErr ImportFrags(CharPtr fnam,CharPtr seq);

/*******************************************************
*Kaca - functionality related to contact maps
********************************************************/
void NewContactMap(BoolPtr *ppbContactMap, Int4 iSequenceLen);
void DumpContactMap(BoolPtr pbContactMap, Int4 iSequenceLen, FILE* pout, CharPtr pcSequence);
Int4 SumChainLen(PDNMM  pdnmm);
CharPtr ConcProtSequence(PDNMM  pdnmm);
Int4 TranslateResNum(ValNodePtr pvnMolLen, Int4 iMolId, Int4 iResId);
Int4 CreateContactMap(BoolPtr *ppbContactMap,DValNodePtr pdnListPmmds,Int2 iModelNum,Boolean bInclusiveWindow,Int2 iWindowSize,Boolean bModelling,POT_TYPE pot_type);



#ifdef __cplusplus
}
#endif

#endif /* MMDBTRAJ_PUB_H */

/*
$Log: mmdbtraj_pub.h,v $
Revision 1.121  2008/12/11 17:19:26  chogue
added function in rotate.c - updated header file

Revision 1.120  2004/09/24 19:08:19  hfeldman
Added import fragment option to maketrj

Revision 1.119  2004/08/30 20:36:53  egarderm
Changed gen zero sized to 10,000 (smaller gen needed for non-CASP proteins)

Revision 1.118  2004/07/27 22:48:11  ksnyder
Added gap shift Boolean to maketrjparam structure; changed HomTraj function prototype; added PHLI structure

Revision 1.117  2004/07/14 18:03:17  egarderm
External definitions for vscore functions

Revision 1.116  2004/07/12 16:57:59  mjdumont
HomTraj doesn't require query FASTA - obtained from alignments.\nLow complexity filtering in

Revision 1.115  2004/07/06 20:21:30  hfeldman
Changed homtraj arguments

Revision 1.114  2004/07/06 14:43:46  egarderm
Smaller gen 0 size to ensure no overflow to any storage ints

Revision 1.113  2004/06/30 20:55:11  egarderm
Changed gen zero size to 50000 to provide a better starting point

Revision 1.112  2004/06/23 17:35:56  hfeldman
Added define to build from fragments onlyx

Revision 1.111  2004/06/16 22:07:22  hfeldman
Added params to RMSD functions

Revision 1.110  2004/06/09 17:23:55  mjdumont
Broke out RPSBlastTraj mega-function to HomTraj and Hom_* functions

Revision 1.109  2004/04/21 20:10:46  egarderm
Changed genX size back to 100

Revision 1.108  2004/03/23 19:03:51  egarderm
Doubled generation size to account for fast protein

Revision 1.107  2003/12/19 19:26:13  hfeldman
added new f'n header

Revision 1.106  2003/12/11 17:32:54  hfeldman
Increased generation size from 50 to 100

Revision 1.105  2003/10/01 17:04:37  feldman
Added occupancy redundancy as optional flag for adding to BD tree

Revision 1.104  2003/09/23 19:13:03  feldman
Added errorfile option to foldtraj for homtraj
Made foldtraj gradually increase laxness if gets stuck a lot, like in DFP

Revision 1.103  2003/09/23 17:09:30  feldman
Added error message output for homtraj

Revision 1.102  2003/08/29 22:13:22  ksnyder
Added new VERBOSITY setting HOMSRV

Revision 1.101  2003/07/17 13:31:59  feldman
Added f'n header

Revision 1.100  2003/07/14 19:05:18  ksnyder
*** empty log message ***

Revision 1.99  2003/05/23 19:26:32  feldman
Added soft atom padding to CHARMM

Revision 1.98  2003/05/01 20:34:29  feldman
Added f'n header

Revision 1.97  2003/04/04 21:49:16  feldman
Added buried only to savechis and fix helices option to maketrj/unfoldtraj

Revision 1.96  2003/03/25 17:46:18  feldman
Added zhang window param to maketrj

Revision 1.95  2003/03/11 20:44:46  lewis
Maxgen is 250 now

Revision 1.94  2003/03/11 19:27:44  feldman
New generation sizes

Revision 1.93  2003/03/07 23:33:13  feldman
MAde programprogress an Int4 and added wwwget paramblock

Revision 1.92  2003/02/19 17:39:29  feldman
Make generation sizes smaller for beta test

Revision 1.91  2003/01/24 16:44:37  feldman
Moved some defines around, made CHARMM callable as a thread

Revision 1.90  2003/01/24 16:34:16  feldman
Added generation sizes

Revision 1.89  2003/01/24 16:33:50  feldman
Added generation defines

Revision 1.88  2003/01/08 16:53:46  feldman
Added parameter for bioseq to loadabiostruc

Revision 1.87  2002/10/23 17:59:46  kaca
added declarations for contact maps

Revision 1.86  2002/09/26 13:23:22  michel
Moved BuildMIMEBiostruc to mmdbtrajlib

Revision 1.85  2002/09/20 14:54:40  feldman
New function headers

Revision 1.84  2002/08/22 21:07:05  feldman
Added a function header

Revision 1.83  2002/08/20 19:28:25  michel
multalign results output to command line parameter file while keeping track of progresss

Revision 1.82  2002/08/13 19:44:41  kaca
added POT_SB2 type

Revision 1.81  2002/07/31 22:22:59  feldman
Added two outside angles to ComputeTakeoffAngles
Fixed some bugs in fragment minimization for HM

Revision 1.80  2002/07/29 19:41:20  phan
Added template struc option

Revision 1.79  2002/07/25 16:30:05  feldman
Added Tunnel prob

Revision 1.78  2002/07/20 16:41:15  feldman
Added new parameter for LOTS_OF_RAM

Revision 1.77  2002/07/15 14:26:21  feldman
Added Getpcis function

Revision 1.76  2002/07/09 14:09:11  phan
Added parameter to rpsblast

Revision 1.75  2002/07/05 16:10:05  feldman
Added f'n to count total frag probability

Revision 1.74  2002/07/02 20:44:14  kaca
added POT_SB potential

Revision 1.73  2002/06/18 16:07:24  michel
Abstracted structure contact counts

Revision 1.72  2002/06/14 20:47:17  feldman
Added fragment length to data structures

Revision 1.71  2002/06/13 13:09:53  michel
moved CalcExtendedResidues to mmdbtraj_pub.h

Revision 1.70  2002/06/12 22:14:54  feldman
Added MSA functions

Revision 1.69  2002/03/26 16:52:47  kaca
added 4D potential function defs

Revision 1.68  2002/02/25 22:08:01  feldman
Added checksumming of text files potentials, cbdata and skel.prt
If checksum fails, get error
Changed bailing from foldtrajlite to exit curses first if in it

Revision 1.67  2002/01/09 20:17:08  michel
Changed EGORDJKAdded name to EGORDoubleJacknife and added extra function parameter for detailed jk scores

Revision 1.66  2001/12/21 22:45:58  michel
Added (Secondary structure) Enhanced GOR prediction with double jacknifing

Revision 1.65  2001/12/03 15:51:37  feldman
Added CHARMM NBX structure

Revision 1.64  2001/11/09 22:20:46  michel
Made GetIndex generic, added modelling feature to CreatePsBlist, fixed AddtoBryantPotential, other fixes

Revision 1.63  2001/11/08 16:55:41  feldman
updated function headers

Revision 1.62  2001/10/26 20:24:25  michel
minor bug fix

Revision 1.61  2001/10/23 15:09:31  feldman
Added coordinate backup and restore

Revision 1.60  2001/10/18 22:17:15  michel
Added functionality for BL table writing

Revision 1.59  2001/10/10 19:53:54  feldman
Make H-bonds a bit more flexible (120 degrees rather than 150) and fixed up coding
for H-bonds in randwalk.c

Revision 1.58  2001/10/02 21:50:43  feldman
Added double-precision vector functions and used with Energy computations

Revision 1.57  2001/10/01 14:43:03  feldman
Added BD Node counting function

Revision 1.56  2001/10/01 14:28:20  feldman
Added CHARMM group option

Revision 1.55  2001/09/26 15:31:11  feldman
Added Solvation term to energy

Revision 1.54  2001/09/24 14:53:42  feldman
Added CHARMM fucntions

Revision 1.53  2001/09/13 16:27:51  feldman
Had to add a flag to Uturn database init to indicate if should count entries

Revision 1.52  2001/09/10 20:27:46  feldman
Added Uturn option to maketrj

Revision 1.51  2001/09/07 21:50:03  feldman
Added probability as uturn creation parameter

Revision 1.50  2001/09/05 16:59:19  feldman
Changed Uturn function to take PNN

Revision 1.49  2001/08/31 18:31:14  feldman
Added U-turn ASN.1 header and U-turn database related function prototypes

Revision 1.48  2001/08/30 14:45:21  michel
Changed file identifier to string for DumpBLPotential

Revision 1.47  2001/08/02 19:04:21  michel
bug fixes in BL potential functions

Revision 1.46  2001/07/26 19:18:31  feldman
Added new data structure to stores fragments on disk with no platform-
dependent-sized pointers

Revision 1.45  2001/07/26 19:12:07  michel
Added new BL potential functions

Revision 1.44  2001/07/13 21:15:33  feldman
Made amber output more quiet

Revision 1.43  2001/07/13 17:40:17  feldman
Removed id parameter from BD functions

Revision 1.42  2001/07/12 21:08:03  feldman
Added more AMBER related functions

Revision 1.41  2001/06/27 01:32:25  michel
Added new functions to generate Bryant-Lawrence potential

Revision 1.40  2001/06/26 20:27:01  feldman
Moved clustal stuff to library

Revision 1.39  2001/06/19 18:41:47  phan
Changed rgyr f'ns to take pmmd

Revision 1.38  2001/05/25 21:50:00  feldman
Added datafile path global variable

Revision 1.37  2001/04/17 21:32:06  feldman
Added fragment editing to Vistraj - not quite fully functional yet
Also can now give trajectory graph name as argument #1 when running
and fixed potential sscanf problems so entering garbage into numeric
entry fields results in an error being shown (instead of accepting it)

Revision 1.36  2001/04/06 14:29:00  feldman
Completed fold-at-home Client-Server v1.0 now complete and ready for testing

Revision 1.35  2001/04/05 14:02:38  feldman
Fixed so curses should work on both HP and Alpha now, I hope.
Also added a few more error codes to trajstore

Revision 1.34  2001/04/04 21:26:00  feldman
-Removed my BSPrintf, now using the one in SLRIlib - thus needed to add
	this library to makefiles
-changed all fgets to FileGets

Revision 1.33  2001/03/30 22:22:26  feldman
Minor spelling and grammar fixes, bug fixes, and added
foldtraj ASCII screensaver ability

Revision 1.32  2001/03/29 20:56:58  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.31  2001/03/27 20:24:21  feldman
Tidied up thread communication and windows 3D viewer launching

Revision 1.30  2001/03/23 21:08:07  feldman
Added monitor for foldtraj inside vistraj

Revision 1.29  2001/03/14 16:25:54  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.28  2001/03/13 15:08:16  feldman
Added ability to save sidechain chi angles when making a new
trajectory distribution - foldtraj will then use this

Revision 1.27  2001/03/09 17:33:58  feldman
Added code to deal with rotamers.  Now rotamer chi angles can
be encoded into a 4-byte integer and used by foldtraj when
building the structures.  These can be edited by VisTraj as well

Revision 1.26  2001/03/07 21:49:47  feldman
Made many important fixes to VisTraj such as bounds checking of
user-entered information, and tidied up interface further.
Fixed a few minor bugs as well, and tested RPSTraj feature of
VisTraj.

Revision 1.25  2001/02/26 22:19:22  feldman
Changed fragments to allow for multiple possible fragments at
each residue

Revision 1.24  2001/02/15 20:29:03  feldman
Changed maketrj to take a parameter block

Revision 1.23  2001/02/13 02:18:54  feldman
added needed filename define

Revision 1.22  2001/02/09 20:18:33  feldman
Added fragments to PTGS

Revision 1.21  2001/02/08 22:28:11  feldman
Fixed important inconsistency in modified amino acids - now for
encoded sequence (*xxxxA) the xxxx is always the dictionary index
of the residue, thus we no longer need to add 1 or 2 to it in
cases where it is found at the N- or C-terminus; adding and
deleting residues with AlterResiduesInSequence will update any
affected residues correctly as well

Revision 1.20  2001/02/08 19:46:36  feldman
Added code to improve SubstituteNames and make it do a lot
for checking of atomnames when entering constraints

Revision 1.19  2001/02/06 18:37:39  feldman
-Moved around a few function headers
-Split up conform header into private and public sections

Revision 1.17  2001/01/23 18:02:09  feldman
Added headers for distance constraints functions

Revision 1.16  2001/01/16 22:01:35  feldman
Updated contraints to now have the proper 6 degrees of freedom.
Support still needs to be added for Phi-Psi walk

Revision 1.15  2001/01/12 20:01:50  feldman
Added functionality to maketrj to call RPSBlast and then look up
results in local CDD database - this functionality makes aligntraj
obsolete, so that program will not be developed further

Major functions were cut and paste from cddserver.c in toolkit and
adapted into clust.c.
Note this addition required 2 more parameters in the .foldtrajrc file
for RPSBlast database name and CDD path

Revision 1.14  2000/12/15 00:06:41  feldman
Now Val2Trj and unfoldtraj works correctly for Phi-Psi walk (I think)

Revision 1.13  2000/10/25 15:15:13  feldman
Made further updates, multiple model support is correct now
and relocated to a single function for loading correct model
and extensive error handling

Revision 1.12  2000/10/24 20:57:24  feldman
Added better support for models, now avoids Vector model and
one-coord per residue models when possible and always uses
All-atom models

Revision 1.11  2000/10/06 18:14:41  feldman
Change Sarray so element value indicates residue it bonds with

Revision 1.10  2000/09/15 20:36:54  feldman
Added angle and dihedral to ASN trajgraph dist. constraints
Now counting distance constraint violations

Revision 1.9  2000/08/23 21:36:51  feldman
Updated a few new functions

Revision 1.8  2000/08/18 19:00:17  adrian
changed r14seqfilecalc to r14bscalc
added DSSP function definitions to slriaccsurf.h

Revision 1.7  2000/08/17 17:53:36  feldman
tidied up makefile paths and includes

Revision 1.6  2000/08/14 20:32:44  feldman
Added some amber functions headers, and changed <> to "" for some
includes

Revision 1.5  2000/08/14 20:26:14  feldman
Added LIBCALLBACK to FreeZhangAtmNode

Revision 1.4  2000/08/09 19:17:07  feldman
-minor update/bugfixes added
-version.h contains default version number

Revision 1.3  2000/07/13 20:31:52  feldman
Tidied up header further

Revision 1.2  2000/07/10 15:40:58  feldman
Updated TrajgraphWrite to not call TGInit and TGClose, thus
now just call TGInit once at start of program, TGClose at end
(also removed these calls from updatePcis and alterresiduesequence
)

Revision 1.1  2000/07/07 21:30:26  feldman
Tidied up .h header files and removed some
inter-dependencies


*/

