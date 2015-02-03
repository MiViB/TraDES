TraDES File Dependencies.
June 5 2012 CWVHogue.

This document lists the files required by each command-line executable program.
It also indicates which *_error.log file each executable will log errors to.
All dependent files must be in the same directory/folder as the executable. 

_________
benchtraj - TraDES benchmarking and unit test.
	- logs errors in error.log - but saved only if - specified on command line.
REQUIRED FILES (same directory):
	bstdt.val - the ASN.1 chemical graph dictionary with added non-standard amino acids
	rotlib.bin.bz2 - the backbone-dependent rotamer library
	phipshiwalk_dict.cdx ,phipshiwalk_dict.dbf ,phipshiwalk_dict.fpt - the Ramachandran dictionary
	database.obs, database.seq - the GOR secondary structure prediciton database
	skel.prt - the stub ASN.1 structure for output
	cbdata - file with custom beta-carbon coordinate placement information
	
_______
seq2trj   - Make *.trj files from FASTA
          - logs errors in seq2trj_error.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	phipshiwalk_dict.cdx ,phipshiwalk_dict.dbf ,phipshiwalk_dict.fpt - the Ramachandran Space dictionary
	cawalk_dict.cdx, cawalk_dict.dbf, cawalk_dict.fpt - the Alpha Carbon space dictionary 
	database.obs, database.seq - the GOR secondary structure prediciton database
	skel.prt - the stub ASN.1 structure for output

_______	
str2trj   - Make *.trj files from 3D structures in ASN.1 format
          - logs errors in str2trj_error.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	phipshiwalk_dict.cdx ,phipshiwalk_dict.dbf ,phipshiwalk_dict.fpt - the Ramachandran Space dictionary
	cawalk_dict.cdx, cawalk_dict.dbf, cawalk_dict.fpt - the Alpha Carbon space dictionary 
	database.obs, database.seq - the GOR secondary structure prediciton database
	skel.prt - the stub ASN.1 structure for output
		
______	
trades    - The TraDES sampling engine (formerly foldtraj)
		  - logs errors in trades_error.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	rotlib.bin.bz2 - the backbone-dependent rotamer library
	skel.prt - the stub ASN.1 structure for output
	blpotential.txt - the Bryant-Lawrence potential
	cbdata - file with custom beta-carbon coordinate placement information
	zhangatm.txt  - Zhang-Delisi potential atom names
	zhangij.txt   - Zhang-Delisi potential atom scores


________	
str2pdb  - Convert ASN.1 files to PDB files (save space, use ASN.1!)
		  - logs errors in error_str2pdb.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids

________	
strSummary  - Report Contents, Sequences, Ramachandran angles of  ASN.1 3d files 
		  - logs errors in error_strSummary.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	
_______
ramangL   - Output Ramachandran Angles for a set of *.val samples
		  - logs errors in error_ramangL.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	
_______
solvate   - Add/remove solvent layers to an ASN.1 3D structure.
		- logs errors in error_solvate.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	blpotential.txt - the Bryant-Lawrence potential
	cbdata - file with custom beta-carbon coordinate placement information
	zhangatm.txt  - Zhang-Delisi potential atom names
	zhangij.txt   - Zhang-Delisi potential atom scores
	
________
solvateL  - Loop form of solvate for TraDES data postprocessing.
		  - logs errors in error_solvate.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	blpotential.txt - the Bryant-Lawrence potential
	cbdata - file with custom beta-carbon coordinate placement information
	zhangatm.txt  - Zhang-Delisi potential atom names
	zhangij.txt   - Zhang-Delisi potential atom scores

______
salign    - Align ASN.1 structure, for movies or dock-by-alignment
		  - logs errors in error_salign.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids

________	
strMerge  - Merges two ASN.1 files - for dock-by-alignment methods
		  - logs errors in error_strMerge.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids

________	
crashchk  - Post-process steric collision testing code 
		  - logs errors in error_crashchk.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids

_______
bin2prt   - Convert ASN.1 binary structure to ASN.1 ascii (for movie making)
		  - logs errors in error_bin2prt.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids

____________	
concatmodels  - Movie making - concatenates ASN.1 ascii structures into one file.
		  - logs errors in error.log
REQUIRED FILES (same directory):
		  - none.
_______	
prt2bin   - Convert ASN.1 ascii to binary (form movie making)
		  - logs errors in error_prt2bin.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids

____________
analyzeMovie  - Output a variety of graphs and collated data for movie files.
		  - logs errors in error_analyzeMovie.log
REQUIRED FILES (same directory):
	bstdt.val, - the ASN.1 chemical graph dictionary with added non-standard amino acids
	blpotential.txt - the Bryant-Lawrence potential
	cbdata - file with custom beta-carbon coordinate placement information
	zhangatm.txt  - Zhang-Delisi potential atom names
	zhangij.txt   - Zhang-Delisi potential atom scores
	
	