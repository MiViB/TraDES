/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "tgraph.h";
static AsnValxNode avnx[1] = {
    {2,NULL,0,0.0,NULL } };

static AsnType atx[80] = {
  {401, "TGraph" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[1],0,&atx[2]} ,
  {0, "tgheader" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[53]} ,
  {402, "TGraph-header" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[3],0,&atx[10]} ,
  {0, "bsstub" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {404, "Biostruc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[54]} ,
  {0, "seq" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[8]} ,
  {406, "TGraph-Seq" ,1,0,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[13]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "bseq" ,128,2,0,1,0,0,0,0,NULL,&atx[11],&atx[9],0,&atx[12]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {403, "Seq-entry" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[4]} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seqlen" ,128,3,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[15]} ,
  {407, "TGraph-SeqLen" ,1,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[26]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "bberrtol" ,128,4,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[17]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "bbacc" ,128,5,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[18]} ,
  {0, "numrot" ,128,6,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[19]} ,
  {0, "incsize" ,128,7,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[20]} ,
  {0, "startbb" ,128,8,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[21]} ,
  {0, "atmbncbb" ,128,9,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[22]} ,
  {0, "atmbncsc" ,128,10,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[23]} ,
  {0, "bumph" ,128,11,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[24]} ,
  {0, "distconstr" ,128,12,0,1,0,0,0,0,NULL,&atx[45],&atx[25],0,&atx[46]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {408, "TGConstraint-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[27],0,&atx[28]} ,
  {0, "res1" ,128,0,0,0,0,0,0,0,NULL,&atx[28],NULL,0,&atx[29]} ,
  {409, "TGConstraint-ResidueID" ,1,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[31]} ,
  {0, "res2" ,128,1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,&atx[30]} ,
  {0, "atomname1" ,128,2,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[32]} ,
  {410, "TGConstraint-AtomName" ,1,0,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[34]} ,
  {0, "atomname2" ,128,3,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[33]} ,
  {0, "meandist" ,128,4,0,0,0,0,0,0,NULL,&atx[34],NULL,0,&atx[35]} ,
  {411, "TGConstraint-Dist" ,1,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[38]} ,
  {0, "mindelta" ,128,5,0,0,0,0,0,0,NULL,&atx[34],NULL,0,&atx[36]} ,
  {0, "maxdelta" ,128,6,0,0,0,0,0,0,NULL,&atx[34],NULL,0,&atx[37]} ,
  {0, "angle1" ,128,7,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[39]} ,
  {412, "TGConstraint-Angle" ,1,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[58]} ,
  {0, "angle2" ,128,8,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[40]} ,
  {0, "dihedral01" ,128,9,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[41]} ,
  {0, "dihedral12" ,128,10,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[42]} ,
  {0, "dihedral23" ,128,11,0,0,0,0,0,0,NULL,&atx[38],NULL,0,&atx[43]} ,
  {0, "prob" ,128,12,0,0,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "constrfile" ,128,13,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[47]} ,
  {0, "trajtype" ,128,14,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[48]} ,
  {0, "trajdiv" ,128,15,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[49]} ,
  {0, "timeout" ,128,16,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[50]} ,
  {0, "walktype" ,128,17,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[51]} ,
  {0, "tgunits" ,128,18,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[52]} ,
  {0, "tunnel" ,128,19,0,1,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {0, "datacdx" ,128,1,0,0,0,0,0,0,NULL,&atx[54],NULL,0,&atx[56]} ,
  {405, "TGraph-trajdata" ,1,0,0,0,0,0,0,0,NULL,&atx[55],NULL,0,&atx[6]} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "datadbf" ,128,2,0,0,0,0,0,0,NULL,&atx[54],NULL,0,&atx[57]} ,
  {0, "datafpt" ,128,3,0,0,0,0,0,0,NULL,&atx[54],NULL,0,NULL} ,
  {413, "SSPrediction-OneState" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[59],0,&atx[61]} ,
  {0, "pred" ,128,0,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[60]} ,
  {0, "cksum" ,128,1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {414, "HomTraj-Frag" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[62],0,&atx[71]} ,
  {0, "last" ,128,0,0,1,0,0,0,0,NULL,&atx[61],NULL,0,&atx[63]} ,
  {0, "alnname" ,128,1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[64]} ,
  {0, "qseq" ,128,2,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[65]} ,
  {0, "tseq" ,128,3,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[66]} ,
  {0, "qstart" ,128,4,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[67]} ,
  {0, "tstart" ,128,5,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[68]} ,
  {0, "seglen" ,128,6,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[69]} ,
  {0, "isgap" ,128,7,0,0,1,0,0,0,&avnx[0],&atx[70],NULL,0,NULL} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {415, "HomTraj-Alignment" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[72],0,NULL} ,
  {0, "qseq" ,128,0,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[73]} ,
  {0, "qlen" ,128,1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[74]} ,
  {0, "cddeval" ,128,2,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[75]} ,
  {0, "alneval" ,128,3,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[76]} ,
  {0, "frags-pri" ,128,4,0,0,0,0,0,0,NULL,&atx[45],&atx[77],0,&atx[78]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[61],NULL,0,NULL} ,
  {0, "frags-sec" ,128,5,0,1,0,0,0,0,NULL,&atx[45],&atx[79],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[61],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "TrajGraph" , "tgraph.h",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module TrajGraph
*
**************************************************/

#define TGRAPH &at[0]
#define TGRAPH_tgheader &at[1]
#define TGRAPH_datacdx &at[53]
#define TGRAPH_datadbf &at[56]
#define TGRAPH_datafpt &at[57]

#define TGRAPH_HEADER &at[2]
#define TGRAPH_HEADER_bsstub &at[3]
#define TGRAPH_HEADER_seq &at[5]
#define TGRAPH_HEADER_bseq &at[8]
#define TGRAPH_HEADER_bseq_E &at[9]
#define TGRAPH_HEADER_seqlen &at[12]
#define TGRAPH_HEADER_bberrtol &at[15]
#define TGRAPH_HEADER_bbacc &at[17]
#define TGRAPH_HEADER_numrot &at[18]
#define TGRAPH_HEADER_incsize &at[19]
#define TGRAPH_HEADER_startbb &at[20]
#define TGRAPH_HEADER_atmbncbb &at[21]
#define TGRAPH_HEADER_atmbncsc &at[22]
#define TGRAPH_HEADER_bumph &at[23]
#define TGRAPH_HEADER_distconstr &at[24]
#define TGRAPH_HEADER_distconstr_E &at[25]
#define TGRAPH_HEADER_constrfile &at[46]
#define TGRAPH_HEADER_trajtype &at[47]
#define TGRAPH_HEADER_trajdiv &at[48]
#define TGRAPH_HEADER_timeout &at[49]
#define TGRAPH_HEADER_walktype &at[50]
#define TGRAPH_HEADER_tgunits &at[51]
#define TGRAPH_HEADER_tunnel &at[52]

#define TGRAPH_TRAJDATA &at[54]

#define TGRAPH_SEQ &at[6]

#define TGRAPH_SEQLEN &at[13]

#define TGCONSTRAINT_DATA &at[26]
#define TGCONSTRAINT_DATA_res1 &at[27]
#define TGCONSTRAINT_DATA_res2 &at[29]
#define TGCONSTRAINT_DATA_atomname1 &at[30]
#define TGCONSTRAINT_DATA_atomname2 &at[32]
#define TGCONSTRAINT_DATA_meandist &at[33]
#define TGCONSTRAINT_DATA_mindelta &at[35]
#define TGCONSTRAINT_DATA_maxdelta &at[36]
#define TGCONSTRAINT_DATA_angle1 &at[37]
#define TGCONSTRAINT_DATA_angle2 &at[39]
#define TGCONSTRAINT_DATA_dihedral01 &at[40]
#define TGCONSTRAINT_DATA_dihedral12 &at[41]
#define TGCONSTRAINT_DATA_dihedral23 &at[42]
#define TGCONSTRAINT_DATA_prob &at[43]

#define TGCONSTRAINT_RESIDUEID &at[28]

#define TGCONSTRAINT_ATOMNAME &at[31]

#define TGCONSTRAINT_DIST &at[34]

#define TGCONSTRAINT_ANGLE &at[38]

#define SSPREDICTION_ONESTATE &at[58]
#define SSPREDICTION_ONESTATE_pred &at[59]
#define SSPREDICTION_ONESTATE_cksum &at[60]

#define HOMTRAJ_FRAG &at[61]
#define HOMTRAJ_FRAG_last &at[62]
#define HOMTRAJ_FRAG_alnname &at[63]
#define HOMTRAJ_FRAG_qseq &at[64]
#define HOMTRAJ_FRAG_tseq &at[65]
#define HOMTRAJ_FRAG_qstart &at[66]
#define HOMTRAJ_FRAG_tstart &at[67]
#define HOMTRAJ_FRAG_seglen &at[68]
#define HOMTRAJ_FRAG_isgap &at[69]

#define HOMTRAJ_ALIGNMENT &at[71]
#define HOMTRAJ_ALIGNMENT_qseq &at[72]
#define HOMTRAJ_ALIGNMENT_qlen &at[73]
#define HOMTRAJ_ALIGNMENT_cddeval &at[74]
#define HOMTRAJ_ALIGNMENT_alneval &at[75]
#define HOMTRAJ_ALIGNMENT_frags_pri &at[76]
#define HOMTRAJ_ALIGNMENT_frags_pri_E &at[77]
#define HOMTRAJ_ALIGNMENT_frags_sec &at[78]
#define HOMTRAJ_ALIGNMENT_frags_sec_E &at[79]
