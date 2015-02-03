#ifndef _objensemb_ 
#define _objensemb_ 

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
*    Generated objects for Module Ensemble-Struc
*    Generated using ASNCODE Revision: 6.9 at Jul 25, 2002 12:15 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objensembAsnLoad PROTO((void));


/**************************************************
*
*    EnsembleStruc
*
**************************************************/
typedef struct struct_EnsembleStruc {
   Int4   endianness;
   struct struct_TGraph_header PNTR   ensemble_header;
   struct struct_Ensemble_residue PNTR   resid;
} EnsembleStruc, PNTR EnsembleStrucPtr;


NLM_EXTERN EnsembleStrucPtr LIBCALL EnsembleStrucFree PROTO ((EnsembleStrucPtr ));
NLM_EXTERN EnsembleStrucPtr LIBCALL EnsembleStrucNew PROTO (( void ));
NLM_EXTERN EnsembleStrucPtr LIBCALL EnsembleStrucAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EnsembleStrucAsnWrite PROTO (( EnsembleStrucPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    EnsembleResidue
*
**************************************************/
typedef struct struct_Ensemble_residue {
   struct struct_Ensemble_residue PNTR next;
   Int4   resnum;
   CharPtr   aa;
   FloatHi   markovfact;
   Int4   dim;
   Int4   timeout;
   struct struct_Ensemble_residue_instance PNTR   instances;
} EnsembleResidue, PNTR EnsembleResiduePtr;


NLM_EXTERN EnsembleResiduePtr LIBCALL EnsembleResidueFree PROTO ((EnsembleResiduePtr ));
NLM_EXTERN EnsembleResiduePtr LIBCALL EnsembleResidueNew PROTO (( void ));
NLM_EXTERN EnsembleResiduePtr LIBCALL EnsembleResidueAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EnsembleResidueAsnWrite PROTO (( EnsembleResiduePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    EnsembleResidueInstance
*
**************************************************/
typedef struct struct_Ensemble_residue_instance {
   struct struct_Ensemble_residue_instance PNTR next;
   Int4   id;
   FloatHi   coord1;
   FloatHi   coord2;
   FloatHi   chiwmean;
   Uint1   iscis;
   FloatHi   chi1;
   FloatHi   chi2;
   FloatHi   chi3;
   FloatHi   chi4;
   FloatHi   energy1;
   FloatHi   energy2;
   FloatHi   energy3;
} EnsembleResidueInstance, PNTR EnsembleResidueInstancePtr;


NLM_EXTERN EnsembleResidueInstancePtr LIBCALL EnsembleResidueInstanceFree PROTO ((EnsembleResidueInstancePtr ));
NLM_EXTERN EnsembleResidueInstancePtr LIBCALL EnsembleResidueInstanceNew PROTO (( void ));
NLM_EXTERN EnsembleResidueInstancePtr LIBCALL EnsembleResidueInstanceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EnsembleResidueInstanceAsnWrite PROTO (( EnsembleResidueInstancePtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objensemb_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

