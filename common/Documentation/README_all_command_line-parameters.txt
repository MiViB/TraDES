TraDES Executable Command-Line Reference.
May 31 2012 CWVHogue.

Core Tools:
               benchtraj - TraDES benchmarking and unit test.
	seq2trj   - Make *.trj files from FASTA
	str2trj   - Make *.trj files from 3D structures in ASN.1 format
	trades    - The TraDES sampling engine (formerly foldtraj)
	str2pdb  - Convert ASN.1 files to PDB files (save space, use ASN.1!)
	strSummary - Reports on Asn.1 file contents, outputs FASTA, Ramachandran angles
	             NOTE - use strSummary for file compatibility inspection. 
	
Analysis Tools: (Input ASN.1 3D files only)
	ramangL   - Output Ramachandran Angles for a set of *.val samples
	solvate   - Add/remove solvent layers to an ASN.1 3D structure.
	solvateL  - Loop form of solvate for TraDES data postprocessing.

Higher-Order Assembly Tools: (For dock-by-alignment methods)
	salign    - Align ASN.1 structure, for movies or dock-by-alignment
	strMerge  - Merges two ASN.1 files - for dock-by-alignment methods
	crashchk  - Post-process steric collision testing code 

Ensemble Movie Tools: (For Unfolding Pathways or Solvated Disorder Sets):
	bin2prt   - Convert ASN.1 binary structure to ASN.1 ascii (for movie making)
	concatmodels  - Movie making - concatenates ASN.1 ascii structures into one file.
	prt2bin   - Convert ASN.1 ascii to binary (form movie making)
	analyzeMovie  - Output a variety of graphs and collated data for movie files.

__________________________________________________________________________________________
COMMAND REFERENCE
For detailed usage info, type the program name with no arguments or just a '-' at the end.

The following is the current list of command-line parameters you will find from TraDES 
executables.


-------benchtraj-------
 Benchtraj   arguments:

  -l  Keep error.log file? [T/F]  Optional
    default = F

-------seq2trj-------
 seq2trj: FASTA sequence to TraDES Trajectory Distribution *.trj.
Use this to set up the trades input file from a sequence.

   arguments:

  -f  REQUIRED Input FASTA Amino Acid Sequence File Name
      (First line must start with >, Protein sequence on second line)
  [File In]
  -t  Trajectory Distribution Type:
      Default is 4 for 3-State Sec. Str. Prediction (GOR Method)
      1 = Uniform, 2 = Standard, 3 = One-State Sec. Str. Prediction,
 [Integer]  Optional
    default = 4
    range from 1 to 4
  -c  Trajectory Distribution ALL-COIL sampling.
Use this for Intrinsic Disorder Simulation
      (T/F) (T overrides -t) [T/F]  Optional
    default = FALSE
  -b  Trajectory Distribution ALL-BETA sampling.
Use this for Urea/GdHCL Denatured Simulation
      (T/F) (overrides -t)
      NOTE: -c T -b T results in 50% COIL, 50% BETA. [T/F]  Optional
    default = FALSE
  -s  Input 3-State Secondary Structure File .ss :  [File In]  Optional
  -d  Input Distance Constraint File Name :  [File In]  Optional
  -o  REQUIRED Output TRJ, SS and ARA File Name (No extension) [File Out]
  -v  Output .CSV Formatted Trajectory File (T/F)  [T/F]  Optional
    default = FALSE
  -q  Quiet Operation (T/F) [T/F]  Optional
    default = TRUE

     Use seq2trj - for args list.
     
     
-------str2trj-------
str2trj: ASN.1 protein 3D structure converted to TraDES *.trj file
   arguments:

  -f  REQUIRED 
Input 3D Asn.1 Structure File Name: (if no extension - looks for .val, .cn3, .prt) [File In]
  -c  PDB Chain Name (default = first one found) [String]  Optional
  -m  Model Number (default = 1) [Integer]  Optional
    default = 1
    range from 1 to 9999
  -s  Start Residue (default = 1) [Integer]  Optional
    default = 1
    range from 1 to 9999
  -e  End Residue (default = 0 represents last residue). [Integer]  Optional
    default = 0
    range from 0 to 9999
  -r  Save Sidechain Rotamer angles
      (Default 0=none; 1=all, 2=buried) [Integer]  Optional
    default = 0
    range from 0 to 2
  -o  Output TRJ File Name: (default = {inputfile}.trj)  [File Out]  Optional
  -t  Unfolding MODE - Temperature (Kelvin) [Real]  Optional
    default = 0
    range from 0 to 9999.0
  -d  Unfolding MODE - Time step (fs)
UNFOLDING MODE Example: -t 375.0 -d 100.0
 [Real]  Optional
    default = 0
    range from 0 to 9999.0
  -x  Peak Width MODE - Standard deviation for x (Degrees) [Real]  Optional
    default = 0
    range from 0 to 45.0
  -y  Peak Width MODE - Standard deviation for y (Degrees)
PEAK WIDTH MODE Example: -x 5 -y 5
 [Real]  Optional
    default = 0
    range from 0 to 22.0
  -q  Quiet Operation (T/F) [T/F]  Optional
    default = TRUE

Use str2trj - for args list.



-------trades-------
TraDES - Trajectory Directed Ensemble Sampling Engine
Used to generate protein 3D structures by conformational space sampling.

   arguments:

  -f  Input Trajectory Distribution File (extension ignored, looks for .trj)
REQUIRED: [File In]
  -o  Output filename to store structure (extension ignored, truncated)
      (default = same as input trajectory distribution): [File Out]  Optional
  -n  Reproducible Random Seed? (0 = Truly Random): [Integer]  Optional
    default = 0
    range from -2147483647 to 2147483647
  -d  Keep input file when done?: [T/F]  Optional
    default = TRUE
  -p  Output PDB files?: [T/F]  Optional
    default = FALSE
  -a  Output ASN.1 files?: [T/F]  Optional
    default = TRUE
  -b  Number of structures to build: [Integer]  Optional
    default = 1
    range from 1 to 300000
  -q  Quiet operation?: [T/F]  Optional
    default = TRUE
  -c  Comparison Asn.1 structure file (extension ignored, looks for: .val .cn3 .prt): [File In]  Optional
  -h  Comparison structure Chain
      (default=first molecule in file): [String]  Optional
    default = -
  -x  Experiment (default=none): [String]  Optional
  -s  Structure Numbering start at:  [Integer]  Optional
    default = 1
    range from 1 to 10000000
  -F  Build only fragment covered range: [T/F]  Optional
    default = FALSE
  -r  Write 20 Ramachandran Angle Files: [T/F]  Optional
    default = FALSE
  -v  ASN.1 MIME BINARY .val Output (FALSE = ASCII .prt): [T/F]  Optional
    default = TRUE
  -l  No. of H20 layers for Ellipsoid and PCA 3D axis Transform
      (Default = 1, 0 for no solvation): [Integer]  Optional
    default = 1
    range from 0 to 10
  -w  Keep H20 in Structure Output: [T/F]  Optional
    default = FALSE
  -t  Stream Progress Update to stdout [T/F]  Optional
    default = FALSE
  -k  Skip PCA rotation after Ellipsoid calculations: [T/F]  Optional
    default = FALSE
  -z  Allow Atom Bounciness to Increase when stuck: (Only set true with input from str2trj!) [T/F]  Optional
    default = FALSE


TO START: You need a *.trj file first (use seq2trj or str2trj). 
     For usage try: seq2trj - OR str2trj - 



-------str2pdb-------
str2pdb  Converts Asn.1 3d .val .prt or .cn3 file to PDB format.
Caution - will not work on Biounits downloaded from NCBI, use Asymmetric Units
   arguments:

  -f  Input Asn.1 3D structure as .prt .val or .cn3 File Name. [File In]
  -m  Model Number: (0=Best Coordinate Set, 9999 = All) [Integer]  Optional
    default = 0
    range from 0 to 9999
  -p  Output PDB File Name? [File Out]  Optional


-------strSummary-----
strSummary  reports on what is inside .val or .cn3 or .prt file to stdout.
Caution - will not work on Biounits downloaded from NCBI, use Asymmetric Units
   arguments:

  -f  Input Asn.1 3D structure as .prt .val or .cn3  File Name. [File In]
  -s  Sequences only (t/F): [T/F]  Optional
    default = FALSE
  -r  Include Ramachandran Angles (t/F): [T/F]  Optional
    default = FALSE

-------ramangL--------
ramangL: 
Report Phi, Psi from TraDES *.val file
 Creates up to 20 *.csv files for R plotting.
   arguments:

  -f  Input VAL File Name (NO EXTENSION). [File In]
  -s  Foldtraj Range Start Number (optinal) [Integer]  Optional
    default = 0
    range from 1 to 9999999
  -r  Foldtraj Range (optional) [Integer]  Optional
    default = 0
    range from 1 to 50000


-------solvate-------
solvate  ASN.1 Structure Utility that Computes Phi, Psi, & Potential Scores,
 Adds solvent layers and rotates ASN.1 structures along ellipsoid axes
   arguments:

  -f  Input Asn.1 3d Structure File Name: [File In]
  -s  Scoring Function, Solvation, PCA Output (T/f): [T/F]  Optional
    default = TRUE
  -r  Ramachandran Angle Output (T/f): [T/F]  Optional
    default = TRUE
  -b  Write Solvated ASN.1 BINARY Biostruc (t/F): [T/F]  Optional
    default = FALSE
  -a  Write Solvated ASN.1 ASCII Biostruc (t/F): [T/F]  Optional
    default = FALSE
  -p  Write Solvated PDB File (t/F): [T/F]  Optional
    default = FALSE
  -h  Omit Hydrophobic Core Solvent: [T/F]  Optional
    default = TRUE
  -n  (Reserved for future use) [T/F]  Optional
    default = FALSE
  -l  Number of Water Layers to Add (1-20) [Integer]  Optional
    default = 1
    range from 1 to 20
  -m  Do PCA coordinate transformation: (T/f) (false leaves coordinates intact)  [T/F]  Optional
    default = TRUE

-------solvateL-------
solvateL
 Computes Solvation, Phi, Psi, & Alternate Potential Scores
 LOOPs over a set of numbered TraDES *.val files
 Creates up to 21 *.csv files for R plotting. Rotates structures along ellipsoid axes
   arguments:

  -f  Input Foldtraj VAL File Name PREFIX 'XXXXX_' (NO EXTENSION). [File In]
  -b  Foldtraj Range Begin Number (optinal) [Integer]  Optional
    default = 0
    range from 1 to 9999999
  -t  Foldtraj Range To Number (optional) [Integer]  Optional
    default = 0
    range from 1 to 9999999
  -s  Scoring Function, Solvation, PCA Output (T/f): [T/F]  Optional
    default = TRUE
  -r  Ramachandran Angle Output (30000 files max!) (T/f): [T/F]  Optional
    default = TRUE
  -v  Write Solvated ASN.1 BINARY Biostruc (t/F): [T/F]  Optional
    default = FALSE
  -a  Write Solvated ASN.1 ASCII Biostruc (t/F): [T/F]  Optional
    default = FALSE
  -p  Write Solvated PDB File (t/F): [T/F]  Optional
    default = FALSE
  -h  Omit Hydrophobic Core Solvent: [T/F]  Optional
    default = TRUE
  -n  Output files with NO solvent molecules (strips crystallographic water) (t/F): [T/F]  Optional
    default = FALSE
  -l  Number of Water Layers to Add (1-20) [Integer]  Optional
    default = 1
    range from 1 to 20
  -m  Do PCA coordinate transformation: (T/f) (false leaves coordinates intact)  [T/F]  Optional
    default = TRUE


-------salign-------
salign - ASN.1 Structure Alignment Program.    arguments:

  -f  Input MMDB filename 1 [File In]  Optional
  -g  Input MMDB filename 2 [File In]  Optional
  -h  Chain 1 (default = first in file) [String]  Optional
    default = -
  -i  Chain 2 (default = first in file) [String]  Optional
    default = -
  -m  Model Number 1 [Integer]  Optional
    default = 1
    range from 1 to 9999
  -n  Model Number 2 [Integer]  Optional
    default = 1
    range from 1 to 9999
  -w  Window size (0=whole protein) [Integer]  Optional
    default = 0
    range from 0 to 32767
  -a  Align: 0 = All Atoms; 1=BB; 2=CA [Integer]  Optional
    default = 0
    range from 0 to 2
  -b  res 1 mol 1 (0=use -w option instead) [Integer]  Optional
    default = 0
    range from 0 to 32767
  -c  res 2 mol 1 (0=use -w option instead) [Integer]  Optional
    default = 0
    range from 0 to 32767
  -d  res 1 mol 2 (0=use -w option instead) [Integer]  Optional
    default = 0
    range from 0 to 32767
  -e  res 2 mol 2 (0=use -w option instead) [Integer]  Optional
    default = 0
    range from 0 to 32767
  -o  Output supressed for R stdin (t/F) [T/F]  Optional
    default = 0
    range from 0 to 0
  -l  Input File - List of Filenames to Align [File In]  Optional
  -j  File List - Start at (use with -l) [Integer]  Optional
    default = 0
    range from 0 to 10000
  -q  Create File Queue from directory [T/F]  Optional
    default = 0
    range from 0 to 0
  -p  Write PDB output file (t/F) [T/F]  Optional
    default = 0
    range from 0 to 0
  -t  Write .prt ASCII text Asn.1 output file (t/F) [T/F]  Optional
    default = 0
    range from 0 to 0

-------strMerge-------
strMerge   Merge One or All chains from one Asn.1 3D files into another Asn.1 file.
   arguments:

  -f  Input 3d Asn.1 file 1 (destination - To) [File In]
  -g  Input 3d Asn.1 file 2 (source - From) [File In]
  -o  Output file Name (extension ignored, default='strMerge_out.val') [File Out]  Optional
  -m  Model Number for file 1 (default - all-atom model) [Integer]  Optional
    default = 0
    range from 0 to 9999
  -n  Model Number for file 2 (default - all-atom model) [Integer]  Optional
    default = 0
    range from 0 to 9999
  -c  Chain to be copied from file 2 to file 1 (default: all chains) [String]  Optional
    default = -
  -a  Output Asn.1 Binary (default 0=.val) (1=.cn3) or Asn.1 ASCII (2=.prt) [Integer]  Optional
    default = 0
    range from 0 to 2
  -q  Verbose mode (T/f) [T/F]  Optional
    default = TRUE

-------crashchk-------
CrashChk 1.0 - finds atomic collisions in ASN.1 structure files with b-d tree   arguments:

  -f  Input MMDB filename [File In]
  -m  Model Number [Integer]  Optional
    default = 1
    range from 1 to 9999
  -l  Model Level (remote files only): 0 = Atoms; 1=BB; 2=all PDB; 3=VECT; [Integer]  Optional
    default = 0
    range from 0 to 3
  -n  Supress same-chain collision reporting (default=FALSE) [T/F]  Optional

 

-------bin2prt-------
 bin2prt Program: Converts Binary Asn.1 files into ASCII   arguments:

  -f  Input 3d Asn.1 Structure File Name (extension ignored, looks for .val .cn3): [File In]
  -u  Unwrap MIME and output bare Biostruc (t/F) [T/F]  Optional
    default = FALSE
  -c  Use Cn3D 4.3 naming: {file}_prt.c3d (t/F) [T/F]  Optional
    default = FALSE
  -s  Strip Secondary Structure Annotation from Biostruc (t/F) [T/F]  Optional
    default = FALSE

-------prt2bin-------
prt2bin   Converts 3d Asn.1 files from ASCII to Binary   arguments:

  -f  Input ASCII 3d Asn.1 Structure File Name (extension ignored, looks for .prt .cn3): [File In]
  -m  MIME-wrap Output Biostruc (if missing, adds protein Bioseqs from internal sequences) (T/f) [T/F]  Optional
    default = TRUE
  -c  Use Cn3D 4.3 naming _val.cn3 (t/F) [T/F]  Optional
    default = FALSE  
  
-------concatmodels-------
concat models   arguments:

  -i  Input MMDB File Base Name (NO EXTENSION) [File In]  Optional
    default = protein
  -n  Number of steps to concatenate [Integer]  Optional
    default = 10
    range from 1 to 2147483647


-------analyzeMovie-------
Analyze Movie Program   arguments:

  -f  Input MMDB filename (enter name of the protein (PATH\X); to open multiple files, the filenames must be of the format PATH\X_movie_#.val) [File In]
  -g  Native structure filename [File In]  Optional
  -n  # files to process [Integer]  Optional
    default = 1
    range from 1 to 20
  -a  Average data (if no, all data from separate simulations are plotted on the same graph)? [T/F]  Optional
    default = FALSE
  -p  # partitions of total number of frames when making contact maps [Integer]  Optional
    default = 4
    range from 1 to 50

 
