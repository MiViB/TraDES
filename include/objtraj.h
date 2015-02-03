#ifndef _objtraj_ 
#define _objtraj_ 

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" { /* } */
#endif


/**************************************************
*
*    Generated objects for Module TrajGraph
*    Generated using ASNCODE Revision: 6.9 at Jun 29, 2004  6:03 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objtrajAsnLoad PROTO((void));


/**************************************************
*
*    TGraph
*
**************************************************/
typedef struct struct_TGraph {
   struct struct_TGraph_header PNTR   tgheader;
   ByteStorePtr   datacdx;
   ByteStorePtr   datadbf;
   ByteStorePtr   datafpt;
} TGraph, PNTR TGraphPtr;


NLM_EXTERN TGraphPtr LIBCALL TGraphFree PROTO ((TGraphPtr ));
NLM_EXTERN TGraphPtr LIBCALL TGraphNew PROTO (( void ));
NLM_EXTERN TGraphPtr LIBCALL TGraphAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL TGraphAsnWrite PROTO (( TGraphPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TGraphHeader
*
**************************************************/
typedef struct struct_TGraph_header {
   struct struct_Biostruc PNTR   bsstub;
   CharPtr   seq;
   ValNodePtr   bseq;
   Int4   seqlen;
   FloatHi   bberrtol;
   FloatHi   bbacc;
   Int4   numrot;
   FloatHi   incsize;
   FloatHi   startbb;
   FloatHi   atmbncbb;
   FloatHi   atmbncsc;
   Int4   bumph;
   struct struct_TGConstraint_Data PNTR   distconstr;
   CharPtr   constrfile;
   Int4   trajtype;
   Int4   trajdiv;
   Int4   timeout;
   Int4   walktype;
   Int4   tgunits;
   FloatHi   tunnel;
} TGraphHeader, PNTR TGraphHeaderPtr;


NLM_EXTERN TGraphHeaderPtr LIBCALL TGraphHeaderFree PROTO ((TGraphHeaderPtr ));
NLM_EXTERN TGraphHeaderPtr LIBCALL TGraphHeaderNew PROTO (( void ));
NLM_EXTERN TGraphHeaderPtr LIBCALL TGraphHeaderAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL TGraphHeaderAsnWrite PROTO (( TGraphHeaderPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TGConstraintData
*
**************************************************/
typedef struct struct_TGConstraint_Data {
   struct struct_TGConstraint_Data PNTR next;
   Int4   res1;
   Int4   res2;
   CharPtr   atomname1;
   CharPtr   atomname2;
   FloatHi   meandist;
   FloatHi   mindelta;
   FloatHi   maxdelta;
   FloatHi   angle1;
   FloatHi   angle2;
   FloatHi   dihedral01;
   FloatHi   dihedral12;
   FloatHi   dihedral23;
   FloatHi   prob;
} TGConstraintData, PNTR TGConstraintDataPtr;


NLM_EXTERN TGConstraintDataPtr LIBCALL TGConstraintDataFree PROTO ((TGConstraintDataPtr ));
NLM_EXTERN TGConstraintDataPtr LIBCALL TGConstraintDataNew PROTO (( void ));
NLM_EXTERN TGConstraintDataPtr LIBCALL TGConstraintDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL TGConstraintDataAsnWrite PROTO (( TGConstraintDataPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SSPredictionOneState
*
**************************************************/
typedef struct struct_SSPrediction_OneState {
   CharPtr   pred;
   CharPtr   cksum;
} SSPredictionOneState, PNTR SSPredictionOneStatePtr;


NLM_EXTERN SSPredictionOneStatePtr LIBCALL SSPredictionOneStateFree PROTO ((SSPredictionOneStatePtr ));
NLM_EXTERN SSPredictionOneStatePtr LIBCALL SSPredictionOneStateNew PROTO (( void ));
NLM_EXTERN SSPredictionOneStatePtr LIBCALL SSPredictionOneStateAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SSPredictionOneStateAsnWrite PROTO (( SSPredictionOneStatePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    HomTrajFrag
*
**************************************************/
typedef struct struct_HomTraj_Frag {
   struct struct_HomTraj_Frag PNTR next;
   struct struct_HomTraj_Frag PNTR   last;
   CharPtr   alnname;
   CharPtr   qseq;
   CharPtr   tseq;
   Int4   qstart;
   Int4   tstart;
   Int4   seglen;
   Uint1   isgap;
} HomTrajFrag, PNTR HomTrajFragPtr;


NLM_EXTERN HomTrajFragPtr LIBCALL HomTrajFragFree PROTO ((HomTrajFragPtr ));
NLM_EXTERN HomTrajFragPtr LIBCALL HomTrajFragNew PROTO (( void ));
NLM_EXTERN HomTrajFragPtr LIBCALL HomTrajFragAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL HomTrajFragAsnWrite PROTO (( HomTrajFragPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    HomTrajAlignment
*
**************************************************/
typedef struct struct_HomTraj_Alignment {
   CharPtr   qseq;
   Int4   qlen;
   FloatHi   cddeval;
   FloatHi   alneval;
   struct struct_HomTraj_Frag PNTR   frags_pri;
   struct struct_HomTraj_Frag PNTR   frags_sec;
} HomTrajAlignment, PNTR HomTrajAlignmentPtr;


NLM_EXTERN HomTrajAlignmentPtr LIBCALL HomTrajAlignmentFree PROTO ((HomTrajAlignmentPtr ));
NLM_EXTERN HomTrajAlignmentPtr LIBCALL HomTrajAlignmentNew PROTO (( void ));
NLM_EXTERN HomTrajAlignmentPtr LIBCALL HomTrajAlignmentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL HomTrajAlignmentAsnWrite PROTO (( HomTrajAlignmentPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objtraj_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

