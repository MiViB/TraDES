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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Following are all the constant bond angles and lengths used, mostly in the rotamer
   dictionary. The format is BL for a bond length or BA for a bond angle or CHIX where
   X is a number, for a chi dihedral angle.  This is immediately followed by the one-
   letter code for the amino acid it applies to; if omitted, it is constant for all
   amino acids.  Lastly the two (for bonds) or three (for angles) atoms involved are
   given using standard PDB atom names.  All distances in Angstroms, all angles in degrees */
/* database average virtual angles and lengths for cis-Proline */
#define CIS_ETA 57.81F
#define CIS_ZETA 56.56F
#define CIS_OCACA 84.46F
#define CIS_OCA 2.388F
/* bond lengths for hydrogens from B. R. Galin, Thesis */
#define BL_OH 0.96F
#define BL_CH 1.08F
#define BL_NH 0.98F
/* NH3 on LYS */
#define BL_NH3 1.04F
#define BL_GUANNH 1.00F
/* an "educated" guess */
#define BL_SH 0.96F
/* chose BA_CACBH according to geometry only */
#define BA_XCH 109.5F
#define BA_XOH 109.5F
#define BA_ARCH 120.0F
#define BA_HISCH 126.0F
#define BA_XNH 120.0F
#define BL_NH3H 1.04F
#define BA_CANH3H 109.5F
/* all values below from Engh and Huber, Acta Cryst. (1991) A47, 392-400 */
/* except chi3, chi4 angles which are from PKB data */
#define BL_NCA 1.458F
#define BL_NH3CA 1.491F
#define BL_CAC 1.525F
#define BA_NCAC 111.2F
#define BLSD_NCA 0.019F
#define BLSD_CAC 0.021F
/* approximate s.d. for all sidechain chi angles except Pro which is
   built in a special way */
#define CHI34SD 15.0F
/* these s.d. only used for termini */
#define BLSD_NH3CA 0.021F
#define BASD_NCAC 2.8F
/* for C-terminus only */
#define BL_CO 1.231F
#define BL_COXT 1.249F
#define BA_CACO 120.8F
#define BA_CACOXT 117.0F
/* next three only used in determining N- and C-terminus */
#define BL_NC 1.329F
#define BL_P_NC 1.341F
#define BLSD_NC 0.015F
#define BA_CNCA 121.7F
#define BASD_CNCA 1.8F
#define BA_CACN 116.2F
#define BASD_CACN 2.0F
#define BA_P_CACN 116.9F
#define BASD_P_CACN 1.5F
#define BL_C_CBSG 1.808F
#define BA_C_CACBSG 114.4F
/* from pkb */
#define BL_C_SGSG 2.038F
/* a conservative guess */
#define BLSD_C_SGSG 0.05F
#define BL_D_CBCG 1.503F
#define BA_D_CACBCG 112.6F
#define BL_D_CGOD1 1.231F
#define BA_D_CBCGOD1 120.8F
#define BL_D_CGOD2 1.249F
#define BA_D_CBCGOD2 118.4F
#define BL_E_CBCG 1.520F
#define BA_E_CACBCG 114.1F
#define BL_E_CGCD 1.503F
#define BA_E_CBCGCD 112.6F
#define BL_E_CDOE1 1.231F
#define BA_E_CGCDOE1 120.8F
#define BL_E_CDOE2 1.249F
#define BA_E_CGCDOE2 118.4F
#define CHI3_E 18.7903F
#define BL_F_CBCG 1.502F
#define BA_F_CACBCG 113.8F
#define BL_F_CGCD1 1.384F
#define BA_F_CBCGCD1 120.7F
#define BL_F_CGCD2 1.384F
#define BA_F_CBCGCD2 120.7F
#define BL_F_CD1CE1 1.382F
#define BA_F_CGCD1CE1 120.7F
#define BL_F_CD2CE2 1.382F
#define BA_F_CGCD2CE2 120.7F
#define BL_F_CE2CZ 1.382F
#define BA_F_CD2CE2CZ 120.0F
/* 2nd hydrogen on GLY */
#define BL_G_CA2HA 1.08F
#define BL_G_NCA 1.451F
#define BL_G_CAC 1.516F
#define BA_G_NCAC 112.5F
/* HIS "E" used, protonated at NE2 */
#define BL_H_CBCG 1.497F
#define BA_H_CACBCG 113.8F
#define BL_H_CGND1 1.371F
#define BA_H_CBCGND1 121.6F
#define BL_H_CGCD2 1.356F
#define BA_H_CBCGCD2 129.1F
#define BL_H_ND1CE1 1.319F
#define BA_H_CGND1CE1 105.6F
#define BL_H_CD2NE2 1.374F
#define BA_H_CGCD2NE2 106.5F
#define BL_I_CBCG1 1.530F
#define BA_I_CACBCG1 110.4F
#define BL_I_CBCG2 1.521F
#define BA_I_CACBCG2 110.5F
#define BL_I_CG1CD1 1.513F
#define BA_I_CBCG1CD1 113.8F
#define BL_K_CBCG 1.520F
#define BA_K_CACBCG 114.1F
#define BL_K_CGCD 1.520F
#define BA_K_CBCGCD 111.3F
#define BL_K_CDCE 1.520F
#define BA_K_CGCDCE 111.3F
#define BL_K_CENZ 1.489F
#define BA_K_CDCENZ 111.9F
#define CHI3_K 171.131F
#define CHI4_K 179.21F
#define BL_L_CBCG 1.530F
#define BA_L_CACBCG 116.3F
#define BL_L_CGCD1 1.521F
#define BA_L_CBCGCD1 110.7F
#define BL_L_CGCD2 1.521F
#define BA_L_CBCGCD2 110.7F
#define BL_M_CBCG 1.520F
#define BA_M_CACBCG 114.1F
#define BL_M_CGSD 1.803F
#define BA_M_CBCGSD 112.7F
#define BL_M_SDCE 1.791F
#define BA_M_CGSDCE 100.9F
#define CHI3_M 173.635F
#define BL_N_CBCG 1.516F
#define BA_N_CACBCG 112.6F
#define BL_N_CGOD1 1.231F
#define BA_N_CBCGOD1 120.8F
#define BL_N_CGND2 1.328F
#define BA_N_CBCGND2 116.4F
#define BL_P_CBCG 1.492F
#define BA_P_CACBCG 104.5F
#define BA_P_NCAC 111.8F
#define BL_P_NCA 1.466F
#define BA_P_CANCD 112.0F
#define BL_P_NCD 1.473F
#define BA_P_CBCGCD 106.1F
#define BL_P_CGCD 1.503F
#define BA_P_CGCDN 103.2F
#define BL_Q_CBCG 1.520F
#define BA_Q_CACBCG 114.1F
#define BL_Q_CGCD 1.516F
#define BA_Q_CBCGCD 112.6F
#define BL_Q_CDOE1 1.231F
#define BA_Q_CGCDOE1 120.8F
#define BL_Q_CDNE2 1.328F
#define BA_Q_CGCDNE2 116.4F
#define CHI3_Q 13.2519F
#define BL_R_CBCG 1.520F
#define BA_R_CACBCG 114.1F
#define BL_R_CGCD 1.520F
#define BA_R_CBCGCD 111.3F
#define BL_R_CDNE 1.460F
#define BA_R_CGCDNE 112.0F
#define BL_R_NECZ 1.329F
#define BA_R_CDNECZ 124.2F
#define BL_R_CZNH1 1.326F
#define BA_R_NECZNH1 120.0F
#define BL_R_CZNH2 1.326F
#define BA_R_NECZNH2 120.0F
#define CHI3_R 175.076F
#define CHI4_R 163.155F
/* delocalized e- in guanido group ensures planar characteristic for chi5 of Arg */
#define CHI5_R 180.0F
#define BL_S_CBOG 1.417F
#define BA_S_CACBOG 111.1F
#define BL_T_CBOG1 1.433F
#define BA_T_CACBOG1 109.6F
#define BL_T_CBCG2 1.521F
#define BA_T_CACBCG2 110.5F
#define BL_V_CBCG1 1.521F
#define BA_V_CACBCG1 110.5F
#define BL_V_CBCG2 1.521F
#define BA_V_CACBCG2 110.5F
#define BL_W_CBCG 1.498F
#define BA_W_CACBCG 113.6F
#define BL_W_CGCD1 1.365F
#define BA_W_CBCGCD1 126.9F
#define BL_W_CGCD2 1.433F
#define BA_W_CBCGCD2 126.8F
#define BL_W_CD1NE1 1.374F
#define BA_W_CGCD1NE1 110.2F
#define BL_W_CD2CE2 1.409F
#define BA_W_CGCD2CE2 107.2F
#define BL_W_CE2CZ2 1.394F
#define BA_W_CD2CE2CZ2 122.4F
#define BL_W_CZ2CH2 1.368F
#define BA_W_CE2CZ2CH2 117.5F
#define BL_W_CD2CE3 1.398F
#define BA_W_CGCD2CE3 133.9F
#define BL_W_CE3CZ3 1.382F
#define BA_W_CD2CE3CZ3 118.6F
#define BL_Y_CBCG 1.512F
#define BA_Y_CACBCG 113.9F
#define BL_Y_CGCD1 1.389F
#define BA_Y_CBCGCD1 120.8F
#define BL_Y_CGCD2 1.389F
#define BA_Y_CBCGCD2 120.8F
#define BL_Y_CD1CE1 1.382F
#define BA_Y_CGCD1CE1 121.2F
#define BL_Y_CD2CE2 1.382F
#define BA_Y_CGCD2CE2 121.2F
#define BL_Y_CE2CZ 1.378F
#define BA_Y_CD2CE2CZ 119.6F
#define BL_Y_CZOH 1.376F
#define BA_Y_CE2CZOH 119.9F

/* non-standard AA values - very approximate */
#define BL_X_CBSEG 2.2F
#define BL_X_CC 1.520F
#define BA_X_XCC 111.3F
#define BL_X_CCDBL 1.35F
#define BA_X_CNC 108.0F
#define BL_X_NC 1.276F
#define BL_X_ND1CE3 1.345F
#define BL_X_NP 1.80F
#define BL_X_PO 1.615F
#define BA_X_XPO 102.6F
#define BA_X_POH 105.0F
#define BL_X_PODBL 1.734F
#define BA_X_XPODBL 127.3F
#define BL_X_PS 1.92F

/* van der Waal radii from Tsai et al (1999) J. Mol. Biol. 290, 253-266 */
/* commented set is from previous foldtraj versions */
/* largest VDW radius used */
#define VDWRAD_MAX 1.9F
#define VDWRAD_C 1.88F
/* 1.5 */
#define VDWRAD_AROCH0 1.61F
#define VDWRAD_AROCH1 1.76F
#define VDWRAD_N 1.64F
/* for h-bonds */
#define VDWRAD_N3H1 1.48F
#define VDWRAD_N3H2 1.54F
#define VDWRAD_N4H3 1.28F
#define VDWRAD_O1H0 1.28F
#define VDWRAD_O2H1 1.32F
/* 1.35 */
#define VDWRAD_O 1.42F
/* 1.35 */
#define VDWRAD_S 1.77F
/* 1.6 */
/* atoms far from the C-Beta, to account for current rotamer rigidness */
#define VDWRAD_X 0.8F
/* guesses for P and SE */
#define VDWRAD_P 1.6F
#define VDWRAD_SE 1.9F
/* very approximate */
#define VDWRAD_H 0.0F
/* = cos150, so H-bonds with bond angles less than 120 will be rejected and treated as a crash */
#define HBOND_ANGLE_MAX -0.5F

#ifdef __cplusplus
}
#endif

#endif /* GEOMETRY_H */

/*  
$Log: geometry.h,v $
Revision 1.4  2001/03/29 20:56:58  feldman
Fixed some Windows compiler warnings, pertaining to typecasting, mostly
doubles to floats

Revision 1.3  2001/03/08 15:03:40  feldman
Fixed license agreement for freely distributable files

Revision 1.2  2000/08/09 19:17:07  feldman
-minor update/bugfixes added
-version.h contains default version number

Revision 1.1  2000/07/07 21:30:26  feldman
Tidied up .h header files and removed some
inter-dependencies

*/
