#include <objall.h>
#include <objmmdb1.h>
#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objtraj.h>

static Boolean loaded = FALSE;

#include <tgraph.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objtrajAsnLoad(void)
{

   if ( ! loaded) {
      NLM_EXTERN_LOADS

      if ( ! AsnLoad ())
      return FALSE;
      loaded = TRUE;
   }

   return TRUE;
}



/**************************************************
*    Generated object loaders for Module TrajGraph
*    Generated using ASNCODE Revision: 6.9 at Jun 29, 2004  6:03 PM
*
**************************************************/


/**************************************************
*
*    TGraphNew()
*
**************************************************/
NLM_EXTERN 
TGraphPtr LIBCALL
TGraphNew(void)
{
   TGraphPtr ptr = MemNew((size_t) sizeof(TGraph));

   return ptr;

}


/**************************************************
*
*    TGraphFree()
*
**************************************************/
NLM_EXTERN 
TGraphPtr LIBCALL
TGraphFree(TGraphPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   TGraphHeaderFree(ptr -> tgheader);
   BSFree(ptr -> datacdx);
   BSFree(ptr -> datadbf);
   BSFree(ptr -> datafpt);
   return MemFree(ptr);
}


/**************************************************
*
*    TGraphAsnRead()
*
**************************************************/
NLM_EXTERN 
TGraphPtr LIBCALL
TGraphAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TGraphPtr ptr;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TGraph ::= (self contained) */
      atp = AsnReadId(aip, amp, TGRAPH);
   } else {
      atp = AsnLinkType(orig, TGRAPH);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TGraphNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == TGRAPH_tgheader) {
      ptr -> tgheader = TGraphHeaderAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_datacdx) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> datacdx = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_datadbf) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> datadbf = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_datafpt) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> datafpt = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TGraphFree(ptr);
   goto ret;
}



/**************************************************
*
*    TGraphAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TGraphAsnWrite(TGraphPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, TGRAPH);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tgheader != NULL) {
      if ( ! TGraphHeaderAsnWrite(ptr -> tgheader, aip, TGRAPH_tgheader)) {
         goto erret;
      }
   }
   if (ptr -> datacdx != NULL) {
      av.ptrvalue = ptr -> datacdx;
      retval = AsnWrite(aip, TGRAPH_datacdx,  &av);
   }
   if (ptr -> datadbf != NULL) {
      av.ptrvalue = ptr -> datadbf;
      retval = AsnWrite(aip, TGRAPH_datadbf,  &av);
   }
   if (ptr -> datafpt != NULL) {
      av.ptrvalue = ptr -> datafpt;
      retval = AsnWrite(aip, TGRAPH_datafpt,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    TGraphHeaderNew()
*
**************************************************/
NLM_EXTERN 
TGraphHeaderPtr LIBCALL
TGraphHeaderNew(void)
{
   TGraphHeaderPtr ptr = MemNew((size_t) sizeof(TGraphHeader));

   return ptr;

}


/**************************************************
*
*    TGraphHeaderFree()
*
**************************************************/
NLM_EXTERN 
TGraphHeaderPtr LIBCALL
TGraphHeaderFree(TGraphHeaderPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   BiostrucFree(ptr -> bsstub);
   MemFree(ptr -> seq);
   AsnGenericChoiceSeqOfFree(ptr -> bseq, (AsnOptFreeFunc) SeqEntryFree);
   AsnGenericUserSeqOfFree(ptr -> distconstr, (AsnOptFreeFunc) TGConstraintDataFree);
   MemFree(ptr -> constrfile);
   return MemFree(ptr);
}


/**************************************************
*
*    TGraphHeaderAsnRead()
*
**************************************************/
NLM_EXTERN 
TGraphHeaderPtr LIBCALL
TGraphHeaderAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TGraphHeaderPtr ptr;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TGraphHeader ::= (self contained) */
      atp = AsnReadId(aip, amp, TGRAPH_HEADER);
   } else {
      atp = AsnLinkType(orig, TGRAPH_HEADER);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TGraphHeaderNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == TGRAPH_HEADER_bsstub) {
      ptr -> bsstub = BiostrucAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_seq) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seq = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_bseq) {
      ptr -> bseq = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqEntryAsnRead, (AsnOptFreeFunc) SeqEntryFree);
      if (isError && ptr -> bseq == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_seqlen) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seqlen = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_bberrtol) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> bberrtol = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_bbacc) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> bbacc = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_numrot) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> numrot = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_incsize) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> incsize = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_startbb) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> startbb = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_atmbncbb) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> atmbncbb = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_atmbncsc) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> atmbncsc = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_bumph) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> bumph = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_distconstr) {
      ptr -> distconstr = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) TGConstraintDataAsnRead, (AsnOptFreeFunc) TGConstraintDataFree);
      if (isError && ptr -> distconstr == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_constrfile) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> constrfile = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_trajtype) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> trajtype = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_trajdiv) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> trajdiv = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_timeout) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> timeout = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_walktype) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> walktype = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_tgunits) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> tgunits = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGRAPH_HEADER_tunnel) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> tunnel = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TGraphHeaderFree(ptr);
   goto ret;
}



/**************************************************
*
*    TGraphHeaderAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TGraphHeaderAsnWrite(TGraphHeaderPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, TGRAPH_HEADER);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> bsstub != NULL) {
      if ( ! BiostrucAsnWrite(ptr -> bsstub, aip, TGRAPH_HEADER_bsstub)) {
         goto erret;
      }
   }
   if (ptr -> seq != NULL) {
      av.ptrvalue = ptr -> seq;
      retval = AsnWrite(aip, TGRAPH_HEADER_seq,  &av);
   }
   AsnGenericChoiceSeqOfAsnWrite(ptr -> bseq, (AsnWriteFunc) SeqEntryAsnWrite, aip, TGRAPH_HEADER_bseq, TGRAPH_HEADER_bseq_E);
   av.intvalue = ptr -> seqlen;
   retval = AsnWrite(aip, TGRAPH_HEADER_seqlen,  &av);
   av.realvalue = ptr -> bberrtol;
   retval = AsnWrite(aip, TGRAPH_HEADER_bberrtol,  &av);
   av.realvalue = ptr -> bbacc;
   retval = AsnWrite(aip, TGRAPH_HEADER_bbacc,  &av);
   av.intvalue = ptr -> numrot;
   retval = AsnWrite(aip, TGRAPH_HEADER_numrot,  &av);
   av.realvalue = ptr -> incsize;
   retval = AsnWrite(aip, TGRAPH_HEADER_incsize,  &av);
   av.realvalue = ptr -> startbb;
   retval = AsnWrite(aip, TGRAPH_HEADER_startbb,  &av);
   av.realvalue = ptr -> atmbncbb;
   retval = AsnWrite(aip, TGRAPH_HEADER_atmbncbb,  &av);
   av.realvalue = ptr -> atmbncsc;
   retval = AsnWrite(aip, TGRAPH_HEADER_atmbncsc,  &av);
   av.intvalue = ptr -> bumph;
   retval = AsnWrite(aip, TGRAPH_HEADER_bumph,  &av);
   AsnGenericUserSeqOfAsnWrite(ptr -> distconstr, (AsnWriteFunc) TGConstraintDataAsnWrite, aip, TGRAPH_HEADER_distconstr, TGRAPH_HEADER_distconstr_E);
   if (ptr -> constrfile != NULL) {
      av.ptrvalue = ptr -> constrfile;
      retval = AsnWrite(aip, TGRAPH_HEADER_constrfile,  &av);
   }
   av.intvalue = ptr -> trajtype;
   retval = AsnWrite(aip, TGRAPH_HEADER_trajtype,  &av);
   av.intvalue = ptr -> trajdiv;
   retval = AsnWrite(aip, TGRAPH_HEADER_trajdiv,  &av);
   av.intvalue = ptr -> timeout;
   retval = AsnWrite(aip, TGRAPH_HEADER_timeout,  &av);
   av.intvalue = ptr -> walktype;
   retval = AsnWrite(aip, TGRAPH_HEADER_walktype,  &av);
   av.intvalue = ptr -> tgunits;
   retval = AsnWrite(aip, TGRAPH_HEADER_tgunits,  &av);
   av.realvalue = ptr -> tunnel;
   retval = AsnWrite(aip, TGRAPH_HEADER_tunnel,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    TGConstraintDataNew()
*
**************************************************/
NLM_EXTERN 
TGConstraintDataPtr LIBCALL
TGConstraintDataNew(void)
{
   TGConstraintDataPtr ptr = MemNew((size_t) sizeof(TGConstraintData));

   return ptr;

}


/**************************************************
*
*    TGConstraintDataFree()
*
**************************************************/
NLM_EXTERN 
TGConstraintDataPtr LIBCALL
TGConstraintDataFree(TGConstraintDataPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> atomname1);
   MemFree(ptr -> atomname2);
   return MemFree(ptr);
}


/**************************************************
*
*    TGConstraintDataAsnRead()
*
**************************************************/
NLM_EXTERN 
TGConstraintDataPtr LIBCALL
TGConstraintDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   TGConstraintDataPtr ptr;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* TGConstraintData ::= (self contained) */
      atp = AsnReadId(aip, amp, TGCONSTRAINT_DATA);
   } else {
      atp = AsnLinkType(orig, TGCONSTRAINT_DATA);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = TGConstraintDataNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == TGCONSTRAINT_DATA_res1) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> res1 = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_res2) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> res2 = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_atomname1) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> atomname1 = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_atomname2) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> atomname2 = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_meandist) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> meandist = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_mindelta) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> mindelta = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_maxdelta) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> maxdelta = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_angle1) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> angle1 = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_angle2) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> angle2 = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_dihedral01) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> dihedral01 = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_dihedral12) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> dihedral12 = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_dihedral23) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> dihedral23 = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == TGCONSTRAINT_DATA_prob) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> prob = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = TGConstraintDataFree(ptr);
   goto ret;
}



/**************************************************
*
*    TGConstraintDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
TGConstraintDataAsnWrite(TGConstraintDataPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, TGCONSTRAINT_DATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> res1;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_res1,  &av);
   av.intvalue = ptr -> res2;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_res2,  &av);
   if (ptr -> atomname1 != NULL) {
      av.ptrvalue = ptr -> atomname1;
      retval = AsnWrite(aip, TGCONSTRAINT_DATA_atomname1,  &av);
   }
   if (ptr -> atomname2 != NULL) {
      av.ptrvalue = ptr -> atomname2;
      retval = AsnWrite(aip, TGCONSTRAINT_DATA_atomname2,  &av);
   }
   av.realvalue = ptr -> meandist;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_meandist,  &av);
   av.realvalue = ptr -> mindelta;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_mindelta,  &av);
   av.realvalue = ptr -> maxdelta;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_maxdelta,  &av);
   av.realvalue = ptr -> angle1;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_angle1,  &av);
   av.realvalue = ptr -> angle2;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_angle2,  &av);
   av.realvalue = ptr -> dihedral01;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_dihedral01,  &av);
   av.realvalue = ptr -> dihedral12;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_dihedral12,  &av);
   av.realvalue = ptr -> dihedral23;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_dihedral23,  &av);
   av.realvalue = ptr -> prob;
   retval = AsnWrite(aip, TGCONSTRAINT_DATA_prob,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SSPredictionOneStateNew()
*
**************************************************/
NLM_EXTERN 
SSPredictionOneStatePtr LIBCALL
SSPredictionOneStateNew(void)
{
   SSPredictionOneStatePtr ptr = MemNew((size_t) sizeof(SSPredictionOneState));

   return ptr;

}


/**************************************************
*
*    SSPredictionOneStateFree()
*
**************************************************/
NLM_EXTERN 
SSPredictionOneStatePtr LIBCALL
SSPredictionOneStateFree(SSPredictionOneStatePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> pred);
   MemFree(ptr -> cksum);
   return MemFree(ptr);
}


/**************************************************
*
*    SSPredictionOneStateAsnRead()
*
**************************************************/
NLM_EXTERN 
SSPredictionOneStatePtr LIBCALL
SSPredictionOneStateAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SSPredictionOneStatePtr ptr;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SSPredictionOneState ::= (self contained) */
      atp = AsnReadId(aip, amp, SSPREDICTION_ONESTATE);
   } else {
      atp = AsnLinkType(orig, SSPREDICTION_ONESTATE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SSPredictionOneStateNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SSPREDICTION_ONESTATE_pred) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> pred = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SSPREDICTION_ONESTATE_cksum) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> cksum = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SSPredictionOneStateFree(ptr);
   goto ret;
}



/**************************************************
*
*    SSPredictionOneStateAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SSPredictionOneStateAsnWrite(SSPredictionOneStatePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SSPREDICTION_ONESTATE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> pred != NULL) {
      av.ptrvalue = ptr -> pred;
      retval = AsnWrite(aip, SSPREDICTION_ONESTATE_pred,  &av);
   }
   if (ptr -> cksum != NULL) {
      av.ptrvalue = ptr -> cksum;
      retval = AsnWrite(aip, SSPREDICTION_ONESTATE_cksum,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    HomTrajFragNew()
*
**************************************************/
NLM_EXTERN 
HomTrajFragPtr LIBCALL
HomTrajFragNew(void)
{
   HomTrajFragPtr ptr = MemNew((size_t) sizeof(HomTrajFrag));

   ptr -> isgap = 0;
   return ptr;

}


/**************************************************
*
*    HomTrajFragFree()
*
**************************************************/
NLM_EXTERN 
HomTrajFragPtr LIBCALL
HomTrajFragFree(HomTrajFragPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   HomTrajFragFree(ptr -> last);
   MemFree(ptr -> alnname);
   MemFree(ptr -> qseq);
   MemFree(ptr -> tseq);
   return MemFree(ptr);
}


/**************************************************
*
*    HomTrajFragAsnRead()
*
**************************************************/
NLM_EXTERN 
HomTrajFragPtr LIBCALL
HomTrajFragAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   HomTrajFragPtr ptr;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* HomTrajFrag ::= (self contained) */
      atp = AsnReadId(aip, amp, HOMTRAJ_FRAG);
   } else {
      atp = AsnLinkType(orig, HOMTRAJ_FRAG);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = HomTrajFragNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == HOMTRAJ_FRAG_last) {
      ptr -> last = HomTrajFragAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_alnname) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> alnname = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_qseq) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> qseq = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_tseq) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> tseq = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_qstart) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> qstart = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_tstart) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> tstart = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_seglen) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seglen = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_FRAG_isgap) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> isgap = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = HomTrajFragFree(ptr);
   goto ret;
}



/**************************************************
*
*    HomTrajFragAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
HomTrajFragAsnWrite(HomTrajFragPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, HOMTRAJ_FRAG);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> last != NULL) {
      if ( ! HomTrajFragAsnWrite(ptr -> last, aip, HOMTRAJ_FRAG_last)) {
         goto erret;
      }
   }
   if (ptr -> alnname != NULL) {
      av.ptrvalue = ptr -> alnname;
      retval = AsnWrite(aip, HOMTRAJ_FRAG_alnname,  &av);
   }
   if (ptr -> qseq != NULL) {
      av.ptrvalue = ptr -> qseq;
      retval = AsnWrite(aip, HOMTRAJ_FRAG_qseq,  &av);
   }
   if (ptr -> tseq != NULL) {
      av.ptrvalue = ptr -> tseq;
      retval = AsnWrite(aip, HOMTRAJ_FRAG_tseq,  &av);
   }
   av.intvalue = ptr -> qstart;
   retval = AsnWrite(aip, HOMTRAJ_FRAG_qstart,  &av);
   av.intvalue = ptr -> tstart;
   retval = AsnWrite(aip, HOMTRAJ_FRAG_tstart,  &av);
   av.intvalue = ptr -> seglen;
   retval = AsnWrite(aip, HOMTRAJ_FRAG_seglen,  &av);
   av.boolvalue = ptr -> isgap;
   retval = AsnWrite(aip, HOMTRAJ_FRAG_isgap,  &av);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    HomTrajAlignmentNew()
*
**************************************************/
NLM_EXTERN 
HomTrajAlignmentPtr LIBCALL
HomTrajAlignmentNew(void)
{
   HomTrajAlignmentPtr ptr = MemNew((size_t) sizeof(HomTrajAlignment));

   return ptr;

}


/**************************************************
*
*    HomTrajAlignmentFree()
*
**************************************************/
NLM_EXTERN 
HomTrajAlignmentPtr LIBCALL
HomTrajAlignmentFree(HomTrajAlignmentPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> qseq);
   AsnGenericUserSeqOfFree(ptr -> frags_pri, (AsnOptFreeFunc) HomTrajFragFree);
   AsnGenericUserSeqOfFree(ptr -> frags_sec, (AsnOptFreeFunc) HomTrajFragFree);
   return MemFree(ptr);
}


/**************************************************
*
*    HomTrajAlignmentAsnRead()
*
**************************************************/
NLM_EXTERN 
HomTrajAlignmentPtr LIBCALL
HomTrajAlignmentAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   HomTrajAlignmentPtr ptr;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* HomTrajAlignment ::= (self contained) */
      atp = AsnReadId(aip, amp, HOMTRAJ_ALIGNMENT);
   } else {
      atp = AsnLinkType(orig, HOMTRAJ_ALIGNMENT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = HomTrajAlignmentNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == HOMTRAJ_ALIGNMENT_qseq) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> qseq = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_ALIGNMENT_qlen) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> qlen = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_ALIGNMENT_cddeval) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> cddeval = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_ALIGNMENT_alneval) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> alneval = av.realvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_ALIGNMENT_frags_pri) {
      ptr -> frags_pri = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) HomTrajFragAsnRead, (AsnOptFreeFunc) HomTrajFragFree);
      if (isError && ptr -> frags_pri == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == HOMTRAJ_ALIGNMENT_frags_sec) {
      ptr -> frags_sec = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) HomTrajFragAsnRead, (AsnOptFreeFunc) HomTrajFragFree);
      if (isError && ptr -> frags_sec == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = HomTrajAlignmentFree(ptr);
   goto ret;
}



/**************************************************
*
*    HomTrajAlignmentAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
HomTrajAlignmentAsnWrite(HomTrajAlignmentPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtrajAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, HOMTRAJ_ALIGNMENT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> qseq != NULL) {
      av.ptrvalue = ptr -> qseq;
      retval = AsnWrite(aip, HOMTRAJ_ALIGNMENT_qseq,  &av);
   }
   av.intvalue = ptr -> qlen;
   retval = AsnWrite(aip, HOMTRAJ_ALIGNMENT_qlen,  &av);
   av.realvalue = ptr -> cddeval;
   retval = AsnWrite(aip, HOMTRAJ_ALIGNMENT_cddeval,  &av);
   av.realvalue = ptr -> alneval;
   retval = AsnWrite(aip, HOMTRAJ_ALIGNMENT_alneval,  &av);
   AsnGenericUserSeqOfAsnWrite(ptr -> frags_pri, (AsnWriteFunc) HomTrajFragAsnWrite, aip, HOMTRAJ_ALIGNMENT_frags_pri, HOMTRAJ_ALIGNMENT_frags_pri_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> frags_sec, (AsnWriteFunc) HomTrajFragAsnWrite, aip, HOMTRAJ_ALIGNMENT_frags_sec, HOMTRAJ_ALIGNMENT_frags_sec_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

