TraDES-2 (c) The Hogue Laboratory, National University of Singapore

Release of 12 June 2012 with Source Code, (excepting commercial library source code). 
Refer to "TraDES_Source_Licenses.txt" for complete copyright and license details.

For more information:

Visit the TraDES website:
	http://trades.blueprint.org/ 

For code/executable downloads:
	ftp://ftp.blueprint.org/pub/trades/

Instructional VIDEOS are under construction here: 
	http://videos.blueprint.org
 

 Really Fast Start - Sample Commands:
 
 ./benchtraj
 ./seq2trj -f ala.faa -o mytrj
 ./trades -f mytrj -b 10 -p T
 ./seq2trj -f prot.faa -o myIDP_allcoil -c T 
 ./trades -f myID_allcoil -b 25 -p T -a F

 Load one of the created *.pdb files into a PDB viewer 
 or use Cn3D on the *.val output.
 
MORE DOCUMENTATION in the /Documentation Directory!
 
 
-----------------------
This file documents the distributable executables within the TraDES-2 package.

Contents of this file:

1. License and References to Cite.

2. TraDES-2 Installation, Testing, TroubleShooting.

3. Install OK, Examples of TraDES - The Basics.

4. Summary of TraDES-2 Executable Command-Line Programs.

5. Additional TraDES Data Files.

6. Getting/building/forking TraDES from Source Code.

7. Known Bugs and Roadmap.

8. Contributors and Acknowledgements.


_____________________________________________________________________________________________________
1. License and References to Cite.

Software contributions from the Hogue Laboratory are provided under a BSD license and are open source.
Other contributions in this package have license and copyright terms of their own. 

In particular the commercial library CodeBase is not open-source. However you can rebuild
the TraDES software with the provided precomipled libraries, subject to the CodeBase sub-licence.

READ the information in the file TraDES_Source_Licences.txt for complete license details.


------------------------
REFERENCES TO CITE

If you use this program for any research which is published, you must
reference the following publications or more recent ones from our laboratory.

1. Feldman, Howard J. and Hogue, Christopher W.V. (2000)  A Fast Method
to Sample Real Protein Conformational Space.  Proteins: Structure, 
Function and Genetics, 39, 112-131.

for the C-alpha walk process and/or

2. Feldman, Howard J. and Hogue Christopher W.V. (2002) Probabilistic Sampling
of Protein Conformations: New Hope for Brute Force? Proteins: Structure,
Function and Genetics, 46, 8-23.

for the default Ramachandran walk processs.

See the included Publications in the /Documentation folder:

fastmethod.pdf  - First TraDES publication 1 above.
finalnewhope.pdf - Second TraDES publication 2 above.
vistraj.pdf - Vistraj paper
6490532_Method_to_construct_protein_stru.pdf  - the U.S. Patent for the TraDES package

_________________________________________________________________________________________________
2. TraDES-2 Installation 
---------------------
TraDES depends on a number of files that must be in the same directory.

The list of file dependencies is in "README_file_dependencies". 

Unzip the archive File to the directory where you wish to exectute the programs.

If you choose to work in a new folder/directory, simply copy all the files 
in the /TraDES_2/ directory (excluding /Documentation and /Samples)
to the folder/directory location where they are to be used. 

Note. TraDES is intentionally a single-threaded program. 
Multiple instances of TraDES can be run,
but each TraDES process must live in a separate folder/directory.

TRADES execuatables are COMMAND-LINE programs, with the exception
of VisTraj. VisTraj is found in previous releases, and is not included in 
the current TraDES-2 release. TraDES-2 is backwards-compatible with VisTraj 1.1.1
which can run under Windows XP.  (We are working on an update of VisTraj).

OTHER THINGS TO INSTALL

---R---

We provide a comprehensive R analysis of TraDES sampling for Intrinsic
Disorder Studies with TRAP (TraDES R Analysis Package)

Install R if you wish to replicate out own publication quality graphs and analysis.
See:
http://cran.r-project.org/ for R downloading and information.
http://videos.blueprint.org/ for R tutorials relevant to TraDES analysis.

---Cn3D--- Visualization 

TraDES extensively utilizes the NCBI Asn.1 3D structure file format.
Movie files output by TraDES can be played in Cn3D and annotated as well, with
rendering settings saved inside the Asn.1 file.

Download Cn3D 4.x from NCBI's website for visualizing *.val files output by TraDES.

----A few words about 3D Structure File Formats---

Any *.val file can be converted to a PDB file - and *.val files are smaller, and
furthermore always compress better than PDB files. If you wish to conserve disk space
with large ensemble samples, convert to PDB only when necessary with the tool 'valtopdb'. 
If you wish to have faster downstream processing or pipeline speed, trades can output
pdb files directly with the -p T option.

Finally note that, downstream processing of TraDES output with programs in this package
all require saving the *.val files, as are output from TraDES by default.
This includes salign, ramangL, solvate, solvateL, concatmodels, anayzeMovie, and valmerge.
If you override Asn.1 output from the command-line currently there is only a very slow
web-based method to recreate Asn.1 files from PDB - that is to use the NCBI VAST search
server.  

If you need to make a *.val file from a PDB file, submit your PDB file here, and
the service will return a link to a fully parsed *.val file BEFORE it runs the full search.
http://www.ncbi.nlm.nih.gov/Structure/VAST/vastsearch.html

NOTE: Currently NCBI's web services send out "Biounits" to Cn3D 4.3. These 
are not compatible with TraDES. Click on the "(x) Assymetric Unit" box before you
download an ASN.1 file.  NCBI also began naming ASN.1 files as *.c3d . These
are equivalent to *.val files, and can be renamed to *.val for use in TraDES.

To obtain ASN.1 versions of PDB database files visit the NCBI ftp site:
   ftp://ftp.ncbi.nlm.nih.gov/mmdb/mmdbdata/mmdb.idx      - for the index of MMDB-ID to PDB-ID conversions
   ftp://ftp.ncbi.nlm.nih.gov/mmdb/mmdbdata/57035.val.gz  - is PDB-ID 1OMD 

One last thing.  ASN.1 files with structures can have different starting objects, 
and these have changed over time.  Cn3D 4.3 has improved its handing and can 
load old- or new- style Asn.1, so be sure you upgrade to the latest version 
of Cn3D to work with TraDES.


-----Testing your installation-----------------------------------------------------
First do a TEST RUN. Start by running a CMD shell from windows (type CMD in the windows run box) or
from a terminal in Unix/Linux:

./benchtraj 

(Windows executables are run without the "./" in front.)

benchtraj    

No command-line parameters needed.  

The benchtraj program exercises most of the underlying library code. 
It will indicate whether or not the software binaries work on your platform.

Benchtraj makes 100 copies of a protein sequence 1CDZ which is a mixed 
alpha-beta protein of 96 amino acids, using a series of standard settings 
that we have employed for over a decade. It always runs on a single core.



If this works it will print out a benchmark timing for
the creation and sampling phases of the code like this:

---example benchtraj output----

benchtraj
6503003 is the CODEBASE version number
 (6502= Intel UNIX 32bit, 6503003 = 64Bit Unix, 6500 = Win32, 6401= NonIntel Unix 32bit)
One moment, opening rotamer library...
Predicting secondary structure and generating trajectory distribution...
MadeTrj complete, unpacking ASN.1 Trajectory Graph
Folding protein...
Benchmark complete.

Summary
-------

          Usr time  Sys time
          --------  --------
Maketrj      1.630     0.220
Foldtraj    30.570     0.700

Benchtraj successful
---------------------------------

If benchtraj does not work on your platform, try the following
troubleshooting section, try a different TraDES distribution or request
a set of binaries for your platform by writing to chogue@blueprint.org



---TroubleShooting on Unix/Linux--------------------------------------
If the benchtraj program fails, it is either an architecture problem
or the platform is missing required libraries (Unix/Linux).

If you issue the command 
ldd benchtraj
you can find out the library dependencies:

e.g. 
[chogue@localhost Fedora14_x86_64]$ ldd benchtraj 
	linux-vdso.so.1 =>  (0x00007fffed9ff000)
	libm.so.6 => /lib64/libm.so.6 (0x0000003a9a600000)
	libpthread.so.0 => /lib64/libpthread.so.0 (0x0000003a99600000)
	libncurses.so.5 => /lib64/libncurses.so.5 (0x0000003aace00000)
	libtinfo.so.5 => /lib64/libtinfo.so.5 (0x0000003aa3e00000)
	libc.so.6 => /lib64/libc.so.6 (0x0000003a99200000)
	/lib64/ld-linux-x86-64.so.2 (0x0000003a98e00000)
	libdl.so.2 => /lib64/libdl.so.2 (0x0000003a99a00000)

This list of libraries must be installed on your system for the program to work.
Ususally these are installed by dependency based package installing systems:

yum on RedHat/Fedora/CentOS
zypper on SuSE
apt-get on Ubuntu
port on Mac OS X

For example if the system complains about a  missing library like libncurses.so.5
you can issue the command (as root or with sudo)

yum install ncurses  

to resolve the library dependency.
---------------------------------------------------------------------


__________________________________________________________________________________________
3. Install OK, Examples of TraDES - The Basics

There is a supplied FASTA formatted "fake" protein sequence file called
prot.faa

First run the sequence to trajectory distribution builder:

./seq2trj -f ala.faa -o polyala

This by default makes a "trajectory distribution file" or *.trj file.

Then run the structure sampling engine:

./trades -f polyala -b 10 -p T

This will make 10 3D structure samples of the protein sequence from ala.faa 
in Asn.1 and in PDB file formats (-t T):
polyala_0000001.val to polyala_0000010.val 
polyala_0000001.pdb to polyala_0000010.pdb
 
To see the .val files you can find the Cn3D structure viewer on the NCBI website. 
Take a look at the structures you have made in a 3D PDB viewer or Cn3D.

You will see a lot of alpha-helix in the resulting structures.


For a quick look to show what sequences are inside an Asn.1 structure use:
./strSummary -f polyala_0000001.val
Or try the example file
./strSummary -f 1YU5.val

For sampling from disordered space try this:

./seq2trj -f ala.faa -o ala_coil -c T

./trades -f ala_coil -b 10 -p T

Look at a few of the results.


-----Great, it works... but how do I study an Intrinsically Disorderd Protein?------

We provide an automated R system "TRAP" TraDES R Analysis Package 
that runs the analysis, makes graphs, assembles ensembles into movies.
Read "README_TraDES_R_Analysis_Package.txt" for instructions.

For now, here is an overview of the low-level process used inside TRAP
that you can follow along with on the command-line.

Simple Overview.

a) Use a disordered prediction package or do a 3D fold sequence analysis and
   remove the "folded" parts where possible. Make a FASTA
   file out of it like the minimal FASTA file "ala.faa", or use a standard 
   FASTA file from a sequence database.
   In this example we use the included file: prot.faa

b) Use seq2trj to make a trajectory distribution (*.trj) file 
   using an ALL-COIL sampling method from the command line:
   
./seq2trj -f prot.faa -o myIDP_allcoil -c T 

3 files are created:
	myIDP_allcoil.trj    This is the Trajectory Distribution File required to make 3d structures
	myIDP_allcoil.ss     This is the Secondary Structure model used (in this case 100% coil with '- c')
	myIDP_allcoil.ara    This is a the area-at-half-max of Ramachandran Space 
	                       at each amino acid (max = 160000 for "uniform" sampling)

d) Create some number (up to 300000) of 3D structure files with the trades engine.

./trades -f myID_allcoil -b 10000 -p T -a F

This will output 10001 files into the current directory, 10000 structures and 1 log file:
	myIDP_allcoil_0000001.pdb
	...
	myIDP_allcoil_0010000.pdb
	_myIDP_allcoil.log  

e) You can analyze the 10000 sample statistical data points contained in the 
_myIDP_allcoil.log file with Excel if you wish. The log file is designed to
be opened in Excel without modification.

For a detailed anaysis of a log file, you can use functions within 
the TraDES R Analysis Package (TRAP), which will be released shortly.

The TraDES R Analysis Package can run the entire analysis system on
a starting FASTA file, sample conformational space using multiple
processors, collect the results, collate the statistics and 
plot publication-quality graphs.

See the TRAP documentation file for details when it arrives in the next release. 
README_TraDES_R_Analysis_Package.txt

We provide a series of videos to get started in R at http://videos.blueprint.org
which use simple TraDES log file processing as examples.


_____________________________________________________________________________________________________
4. Summary of TraDES-2 Executable Command-Line Programs.

Here is a list of the major executable programs in the the TraDES package and their 
primary purpose. 

A detailed list of the command-line parameters is in the file:
README_all_command_line_parameters.txt

A detailed list of file dependencies for each program is in the file:
README_file_dependencies.txt

Core Tools:
	benchtraj - TraDES benchmarking and unit test.
		Essential reading - this file.

	seq2trj   - Make *.trj files from FASTA (formerly maketrj + scripts - now simplified!)
		Essential reading - this file and "README_Using_seq2trj.txt"

	str2trj   - Make *.trj files from 3D structures in ASN.1 format (formerly maketrj + scripts)
		Essential reading: "README_Using_str2trj.txt" and "README_Unfolding.txt"

	trades    - The TraDES sampling engine (formerly foldtraj)
		Essential reading: this file and "README_Using_trades_engine.txt"
	
	str2pdb  - Convert ASN.1 files to PDB files (one-way only!)
	
	strSummary - Shows file format and contents of a .val .prt or .cn3 Asn1 structure
	             dumps out FASTA sequences from the 3D structure file
	             can dump out ramachandran angles as well.

Analysis Tools: (Input ASN.1 3D files only)
	ramangL   - Output Ramachandran Angles for a set of *.val samples
	solvate   - Add/remove solvent layers to an ASN.1 3D structure.
	solvateL  - Loop form of solvate for TraDES data postprocessing.
		Essential reading: "README_output_and_analysis.txt"

Higher-Order Protein Complex Assembly Tools: (For dock-by-superposition methods)
	salign    - Align ASN.1 structure, for movies or dock-by-superposition
	strMerge  - Merges two ASN.1 files - for dock-by-superpositoin methods
	crashchk  - Post-process steric collision testing code 
		Essential reading: "README_dock_by_superposition"

Ensemble Movie Tools: (For Unfolding Pathways or Solvated Disorder Sets):
	bin2prt   - Convert ASN.1 binary structure to ASN.1 ascii (for movie making)
	concatmodels  - Movie making - concatenates ASN.1 ascii structures into one file.
	prt2bin   - Converts ASN.1 ascii to binary (for movie making)
	analyzeMovie  - Output a variety of graphs and collated data for movie files.
		Essential reading: "README_Unfolding.txt", "README_AnalyzeMovie.txt", and
		"README_TraDES_R_Analysis_Package.txt"




__________________________________________________________________________________________________
5. Additional TraDES Documentation and Data Files

**Note - some of the documentation is still under construction. 

/Documentation/

whatsnew.txt - list of recent and past changes to the TraDES software package

mapcommands.txt - converting from earlier verions of foldtraj to trades

README_seq2trj.txt - details on seq2trj options and output

README_str2trj.txt - details about how to use str2trj with trades for near-neighbor sampling

README_Unfolding.txt - details about how to unfold a protein and make a Cn3D compatible movie file.

README_AnalyzeMovie.txt - details about how to use analyzeMovie to look at unfolding or disorder ensembles.

README_TraDES_R_Analysis_Package.txt - A framework for running Intrinsic Disordered Protein 3D ensemble analysis.

README_dock_by_superposition.txt - A fast method to assemble higher-ordered protein complexes
   from X-ray structures and TraDES computed IDP ensembles and screen them for steric compatibility.

README_output_and_analysis.txt - What parameters are computed in the trades log files.

README_VisTraj.txt - *Where to get VisTraj for Windows XP*

README_phospho_ptm_other_amino_acids.txt - How to easily put phospho-Tyr, biotinyl-Lys or other nonstandard amino acids into TraDES.



CompatibilityNotes.txt - some information on the platforms which are
supported by this software and possible problems

SLRIextaa.html - instructions to users for entering non-standard amino acids into their sequence given to InitTraj


/Samples/
SLRIextaa - table of non-standard amino acids used by foldtraj

geometry.h - this is supplied solely for your reference, and may be deleted
	safely.  This file contains all the bond length, bond angles, and so on, used
	by foldtraj, and may be of interest or use to users

phd2ss.c - 
psipred.c  - sample C programs to convert the output of PHDsec or PsiPred
to the input format required by maketrj.  
Use and modify this if you want to use your
favorite secondary structure prediction method to bias trajectory
distribution creation, to get the prediction in the input format that
maketrj will accept.  See also sample.ss

sample.constr - A sample constraint file to demonstrate the format required
by maketrj.  Extensive comments within this file describe the format in
some detail.

sample.ss - A sample secondary structure prediction file, with comments
explaining the format in detail.  Use this format to give to maketrj to
override the internal GOR prediction that will be made otherwise.  Use a
short program like phd2ss.c (see above) to output a file in this format
from the output of PHDsec, for example



		
---Data Files Required by TraDES---
see the file "README_file_dependencies" for the exact dependency list.

blpotential.txt - contains data for the Bryant-Lawrence potential (see
structure summary section below)

bstdt.val - a binary file containing chemical graphs for all residue types,
including non-standard ones

cawalk_dict.* - trajectory distribution dictionary, in CodeBase/FoxPro format,
separated by secondary structure types, and specific for C-alpha random
walks 

cbdata - b-carbon data tables from Rey and Skolnick - don't modify this file

database.obs +
database.seq - the databases used by the GOR secondary structure method,
consisting of 834 proteins from a non-redundant PDB dataset.  Don't modify
these unless you know what you are doing!!

phipsiwalk_dict.* - trajectory distribution dictionary, in CodeBase/FoxPro 
format, separated by secondary structure types, and specific for Phi-Psi
random walks.  Modify .foldtrajrc to determine whether a C-alpha
or Phi-Psi random walk is performed 

rotlib.bin.bz2 - derived from Dunbrack's backbone dependent rotamer
library, this is a database containing probabilities of various
chi angle combinations in different regions of Ramachandran space
and is used to place rotamers probabilistically.  For details:

	Dunbrack, R.L.Jr. and Cohen, F.E.  (1997)  Bayesian statistical
	analysis of protein side-chain rotamer preferences.  Protein
	Science.  6, 1661-1681.

skel.prt - used to build output structure files - do not modify!!


zhangatm.txt +
zhangeij.txt - contains data for the Zhang-DeLisi potential (see
structure summary section below)



_________________________________________________________________________________________________
6. Getting/building/forking TraDES from Source Code.

Initial release is at:
ftp://ftp.blueprint.org/pub/trades/src/

Public access through a version control system is planned, but commits and 
changes to the base source code must go through chogue@blueprint.org .  

If you need training to understand the design principles behind this source code 
and the NCBI C toolkit visit
http://www.blueprint.org/Home/ncbi-toolkit-course-2008 
and go through lectures/exercises 1-9.  

See the file TraDES/makefiles/README_Debug_gdb_valgrind.txt 
for instructions to run with the GNU debugger or the Valgrind memory analyzer.

MS Windows DEBUG configurations are set for the NCBI libraries, but
not at the moment for any of the other libraries or TraDES executables.

Follow on twitter: @cwvhogue to see what I am working on...
_________________________________________________________________________________________________
7. Known Bugs and Roadmap.

The trades core sampling subroutine has been run over 32 BILLION times and is quite
thoroughly debugged. Occasionally proteins get stuck in random walks,
and the program is set to "give up" when it stops making progress, without crashing
and try again.  This happens most often when using *.trj files made with str2trj
and TraDES slows down when trying to reconstruct folded or -nearly folded- structures
as it encounters more "friction" from collisions. 
This is normal.    

Known Bugs:
1. bin2prt, bin2pdb, str2trj, solvate, strMerge
   New NCBI ASN.1 biounits fail - the ones downloaded from the website.
   Click the (x) Assymetric Units check box on NCBI structure summary web page to 
   download Asn.1 files in the original format compatible with TraDES, rename from *.c3d to *.val.
   Or obtain Asn.1 files from the NCBI ftp site, which work as-is, just gunzip.
   ftp://ftp.ncbi.nlm.nih.gov/mmdb/mmdbdata/mmdb.idx  - for the index of MMDB-ID to PDB-ID conversions
   ftp://ftp.ncbi.nlm.nih.gov/mmdb/mmdbdata/57035.val.gz - is PDB-ID 1OMD
   

RoadMap:
1. Continued improvements to TRAP R analysis - this will be an actively updating segment of the code.
2. Ressurection or replacement of VisTraj.
3. Refinements to change the default scoring functions in analyzeMovie.
4. More tools for constructing full-length proteins from folded segments spliced to disordered segments.
6. A new potential function is being built. 
7. New Dunbrack "smoothed" rotamer library to be added. (major rework required)
8. Removal of CodeBase dependencies. (major rework required)

_________________________________________________________________________________________________
8. Contributors and Acknowledgements.

Original Author: Howard Feldman (str2trj, seq2trj, trades, salign, and lots more)
& Christopher Hogue (mmdbapi, b-d tree, Cn3D, st2trj, seq2trj, solvate, trades). 

Contributions from John J. Salama (VisTraj) and Kevin Snyder (valmerge)
Phillipe Phan (analyzeMovie) Adrian Heilbut, Mark Kotowycz, Van Le (trj visualization),
Michel Dumontier (scoring function libraries), Elena Garderman (everywhere), Mingxi Yao (valtopdb),
Gil Alteroviz (GOR implementation),  Boris Steipe (fragment based construction)
Brendan McConkey & Michael Brougham(atom-atom contact scoring function)

Principal Investigator, Software Architect and Coder...:
Christopher W. V. Hogue, B.Sc. Ph.D
twitter: @cwvhogue
chogue@blueprint.org

For more information, visit the TraDES web site at:

http://trades.blueprint.org/

e-mail comments, questions, criticisms and praise to:
chogue@blueprint.org

The latest version of foldtraj will live here:

ftp://ftp.blueprint.org/pub/TraDES/

