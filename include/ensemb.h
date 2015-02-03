/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "ensemb.h";
static AsnType atx[35] = {
  {401, "EnsembleStruc" ,1,0,0,0,0,1,0,0,NULL,&atx[33],&atx[1],0,&atx[4]} ,
  {0, "endianness" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ensemble-header" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {402, "TGraph-header" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[7]} ,
  {0, "resid" ,128,2,0,0,0,0,0,0,NULL,&atx[34],&atx[6],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {403, "Ensemble-residue" ,1,0,0,0,0,0,0,0,NULL,&atx[33],&atx[8],0,&atx[17]} ,
  {0, "resnum" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[9]} ,
  {0, "aa" ,128,1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[11]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "markovfact" ,128,2,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[13]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "dim" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[14]} ,
  {0, "timeout" ,128,4,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[15]} ,
  {0, "instances" ,128,5,0,0,0,0,0,0,NULL,&atx[34],&atx[16],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {404, "Ensemble-residue-instance" ,1,0,0,0,0,0,0,0,NULL,&atx[33],&atx[18],0,&atx[20]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[19]} ,
  {0, "coord1" ,128,1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[21]} ,
  {405, "Angular-Units" ,1,0,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[30]} ,
  {0, "coord2" ,128,2,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[22]} ,
  {0, "chiwmean" ,128,3,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[23]} ,
  {0, "iscis" ,128,4,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[25]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "chi1" ,128,5,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[26]} ,
  {0, "chi2" ,128,6,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[27]} ,
  {0, "chi3" ,128,7,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[28]} ,
  {0, "chi4" ,128,8,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[29]} ,
  {0, "energy1" ,128,9,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[31]} ,
  {406, "Energy-Units" ,1,0,0,0,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {0, "energy2" ,128,10,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[32]} ,
  {0, "energy3" ,128,11,0,0,0,0,0,0,NULL,&atx[30],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "Ensemble-Struc" , "ensemb.h",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = NULL;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module Ensemble-Struc
*
**************************************************/

#define ENSEMBLESTRUC &at[0]
#define ENSEMBLESTRUC_endianness &at[1]
#define ENSEMBLESTRUC_ensemble_header &at[3]
#define ENSEMBLESTRUC_resid &at[5]
#define ENSEMBLESTRUC_resid_E &at[6]

#define ENSEMBLE_RESIDUE &at[7]
#define ENSEMBLE_RESIDUE_resnum &at[8]
#define ENSEMBLE_RESIDUE_aa &at[9]
#define ENSEMBLE_RESIDUE_markovfact &at[11]
#define ENSEMBLE_RESIDUE_dim &at[13]
#define ENSEMBLE_RESIDUE_timeout &at[14]
#define ENSEMBLE_RESIDUE_instances &at[15]
#define ENSEMBLE_RESIDUE_instances_E &at[16]

#define ENSEMBLE_RESIDUE_INSTANCE &at[17]
#define ENSEMBLE_RESIDUE_INSTANCE_id &at[18]
#define ENSEMBLE_RESIDUE_INSTANCE_coord1 &at[19]
#define ENSEMBLE_RESIDUE_INSTANCE_coord2 &at[21]
#define ENSEMBLE_RESIDUE_INSTANCE_chiwmean &at[22]
#define ENSEMBLE_RESIDUE_INSTANCE_iscis &at[23]
#define ENSEMBLE_RESIDUE_INSTANCE_chi1 &at[25]
#define ENSEMBLE_RESIDUE_INSTANCE_chi2 &at[26]
#define ENSEMBLE_RESIDUE_INSTANCE_chi3 &at[27]
#define ENSEMBLE_RESIDUE_INSTANCE_chi4 &at[28]
#define ENSEMBLE_RESIDUE_INSTANCE_energy1 &at[29]
#define ENSEMBLE_RESIDUE_INSTANCE_energy2 &at[31]
#define ENSEMBLE_RESIDUE_INSTANCE_energy3 &at[32]

#define ANGULAR_UNITS &at[20]

#define ENERGY_UNITS &at[30]
