AnalyzeMovie Instructions 
CWVHogue June 2012
-----------------------------------

Quick start:
------------

You can find 1YU5_movie.val.gz and 1YU5.val in the distribution directory
./TraDES-2/
(or in the ./TraDES/build directory after compiling)

Uncompress the included gzipped Asn.1 structure movie file 

1YU5_movie.val.gz 

(7zip works on Windows)

FIRST Look at the file in Cn3D 4.3 or higher. 
Use the down arrow key to go to the first frame.
Use the right arrow key to step frame by frame to see it unfold.
p key to play, s key to stop
Use the left arrow to see it go backwards - the direction of folding.

If you want to see how massive this is as a pdb file, give this a try:
./str2pdb -f 1YU5_movie.val -m 9999

This movie is made with the TraDES package with the unfold-it script (on Unix)
See - README_Unfolding.txt if you want to try it.

SO THEN...
AnalyzeMovie does an analysis on the finished movie file:

./analyzeMovie -f 1YU5_movie.val g 1YU5.val

- a number of .gif and .dat files will be produced in the working directory 
as described below.


What is it?
-----------
AnalyzeMovie is a protein structural analysis software package.  It is used
from the command line only to analyze protein dynamics simulations and output
numerous graphs of certain properties (summarized below).  Graphs are output
both as .GIFs for viewing with graphic software, and as ASCII text files,
which could then be imported into a spreadsheet program for further analysis
or plotting.

The input of the program takes the form of an NCBI ASN.1 MMDB structure file,
which may contain any number of protein structures as a series of 'models'.
There models may be multiple models from an NMR experiment, or may be the same
molecule as it varies in position and shape over time.  Structures in this
file format can be retrieved from the MMDB database at NCBI:

http://www.ncbi.nlm.nih.gov/Structure/

and viewed with NCBI's molecular viewing software, Cn3D, available here:

http://www.ncbi.nlm.nih.gov/Structure/CN3D/cn3d.shtml

The output will be a set of GIF files and a set of *.dat files which are the
text equivalents of the corresponding graphs.  To get a list of possible
arguments to the program, run it from the command line with no arguments.

Arguments
---------
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

	


Note in the following descriptions, an NCBI ASN.1 MMDB file, which normally is
given the .val extension, is a generic structure file format.  It contains a
single 'structure' (a set of one or more molecules).  But it may contain
multiple models.  A model is a set of co-ordinates for the structure.  Each
model is independent of another and there may be 250 or more models in a .val
file.  Multiple models may arise from an NMR experiment where multiple solution
coordinate sets are obtained.  They could also come from a dynamics simulation
where each model represents one 'frame' or timestep of the simulation.  The
coordinates of the atoms changes with each model or frame, but the molecules
themselves (techinically, the chemical graph, specify all the atoms and bonds
and their interconnections) remain the same from model to model.

If your structures are in PDB format, see below for a note on converting them
to the NCBI MMDB format.

IMPORTANT NOTE FOR WINDOWS USERS - due to a 'bug' in the software, when running
the Windows version you must include '.\' in front of the program name to avoid
an error message.  For example ".\analyzeMovie.exe -f blah_movie_1.val
-g blah.val"

-f  -- input filename - this may be a single file or a set of related structure
	files.  If more than one, they all must be in the same directory, and
	all must have the same file name, except they must end in '_movie_X.val'
	where X=1, 2, 3, ... for each succesive file.  If specifying a single
	input file, give the entire path and filename here, while if specifying
	multiple files, just give the path and first part of the filename - the
	_movie_X.val will automatically appended to each filename as they are
	loaded.  When a single file is used, graphs are created up to the
	highest model number represented in the input file.  When using multiple
	input files, the largest model number common to all input file is used.
	For example if three input files are used, one with 63 model, one with
	22 models and one with 35 models, then only the first 22 models of each
	file will be used and the rest ignored.

-g  -- native structure - this should be the path and filename of the native
	structure file (used for computing RMSD and native contact maps).  It
	must again be in NCBI's ASN.1 MMDB file format (.val extension normally)
	and can be obtained from NCBI.  If the native structure is not known,
	you may use one of the input files as the native structure, and then
	simply ignore the RMSD plots of course.

-n  -- number of files - this tells the program how many files are specified in
	the -f option.  By default this is 1, in which case just one file is
	specified in the -f option.  If this option is greater than 1, multiple
	files will be loaded (see -f option above for details).

-a  -- average data? - this only has an effect when the -n option is greater
	than 1.  If true, all data from each file will be averaged and plotted
	on a single plot where it makes sense to do so.  If false, all data will
	be kept separate and plotted as separate lines on the same graph.  The
	default is not to average the data.

-p  -- numer of partitions - rather than making a contact map for every single
	model in the structure file(s), contact maps are aggregated and averaged
	over a number of models, specified here.  With the default of four, four
	contact maps are made corresponding to the first quarter of the models,
	the second quarter, etc.

Graphs
------
Following are a summary of the plots produced by the software and their
interpretation.  Please note there differ depending on what choice you made for
the -a option, if you have input multiple files.  For simplicity lets assume you
gave 2 input files each with 100 models.  Then for radius of gyration vs. time,
for example, you would get a plot of the radius of gyration of each model, with
model number along the 'time' axis (it is assumed the models correspond to some
dynamic simulation of the protein).  With -at, the number plotted for each time
point will be the average radius of gyration between the two input files, at the
given model number.  If -af is used on the other hand, two distinct lines will
be plotted on the same graph, one for each file.  The -at graph would be the
average of these two latter lines.  Averaging (or lack thereof) is carried out
in a similar way for all the other plots.

Lastly, note that while most graphs have corresponding text files containing the
raw data, with the same filename as the graph but a different extension, in some
cases this will vary slightly.  If in doubt, inspect the files generated by the
program after it is run.

Secondary Structure over Time

This plot shows the secondary structure of the protein throughout the folding
trajectory. Each time unit on the x axis corresponds to one generation. Helical
content is shown in red, hydrogen-bonded beta sheet is in blue, and unstructured
residues appear green.

Secondary Structure per Residue over Time

This shows the same information as the previous plot, but on a
residue-by-residue basis, with residue number on the y-axis. The legend on the
right indicates which colours correspond to which secondary structure types.
This can be compared to the secondary structure of the native fold to obtain
some measure of suceess.

Total Surface Accessibility over Time

This is the total, in Angstroms^2, of solvent-accessible surface area of the
protein from one generation to the next. Typically a folded protein will have
half the exposed surface area of an unfolded one.

Surface Accessiblity per Residue over Time

This is the same information as in the previous plot, but split up on a
residue-by-resiude basis, so we can see, at any given time, which residues are
most exposed and most buried. Most often hydrophobic amino acids are found
buried at the core of the protein. The legend indicates the colour scaling.

B-L Potential over Time

This shows the Bryant-Lawrence (B-L) potential over the duration of the folding
trajectory. The B-L potential is a residue-based contact potential, and is
computed from contacts between residues less than 10A apart in space. Ideally
this should decrease as the protein folds. This potential has arbitrary units.

B-L Potential per Residue over Time

This is the same as the previous plot, but broken down on a residue-by-residue
basis, allowing us to see exactly which residues contribute most (either in a
favourable or unfavourable manor) to the total energy at any given time. The
legend shows the mapping of the colour scheme to energy value. This potential
has arbitrary units, but smaller (more negative) is better.

Crease Energy over Time

The Crease Energy is a slight modification of the B-L potential, which provides
an estimate of the steepness of the energy well around the structure being
scored. Again it is unitless but smaller is better. The Crease Energy tends to
be more sensitive to minor structural perturbations than the B-L and can more
reliably 'pick out' low RMSD structures.

Crease Energy per Residue over Time

This shows the Crease Energy at each residue (see Legend for colouring scheme)
over time, and allows us to see which residues are making the most favourable
contacts in the structure.

Z-D Potential over Time

The Zhang-DeLisi (Z-D) potential is an atom-based contact potential, computed by
counting contacts between atoms less than 6A away in space. Its units can be
converted to kcal/mol if some assumptions are made but for our purposes, we
again consider it unitless with smaller (more negative) being better.

Z-D Potential per Residue over Time

The Z-D potential (see previous graph description) plotted at each residue (see
Legend for colour scheme) over time.

Radius of Gyration over Time

Radius of Gyration (Rgyr) gives a measure of the compactness of a structure.
Generally folded proteins are quite compact and it is in fact difficult to find
structures more compact than the native state in many cases. An unfolded protein
on the other hand can have a very large Rgyr. This plot shows Rgyr (in A) over
time (dotted line) as well as the hydrophobic Rgyr (HPRgyr) (solid line). The
latter is the compactness of only the hydrophobic residues, ignoring the rest of
the protein. Since the hydrophobic residues are most often found at the protein
core, it is expected that in a real folded protein, the HPRgyr will be smaller
than Rgyr.

RMSD to native over Time

The RMSD gives a measure of the similarity of a structure to another. This plot
gives the RMSD (in A) between any given generation and the true native fold of
the protein (since we know the correct structure in this case). Ideally it
should go all the way down towards zero RMSD, indicating that the protein folded
into the correct shape.

Radius of Gyration vs. RMSD

This plot shows the Rgyr (in A) on the y-axis, RMSD (in A) on the x-axis, and
B-L potential by colour (see legend). See above plot descriptions for details on
what these quantities are. Exactly one square appears for each generation, but
their time order cannot be explicitly determined from the plot. The correct
sequence can however be deduced by looking at, for example, the RMSD vs. time
plot and comparing the two. This plot gives a good overall summary of what is
happening as the protein folds from one generation to the next.

Native structure contact map

This plot show the contact map of the native protein structure. See the
following description for an explanation of contact maps in general.

Contact Map

A contact map is a good way of summarizing graphically a protein structure.
Along both axes are the residue numbers, or protein sequence. Then, a dot is
placed on the map for every pair of residues that are in contact (within some
arbitrary distance, say 5A). For example if residues 10 and 50 are in contact, a
dot will be placed at x=10, y=50. Since (i,i+1) and (i,i+2) contacts are
ubiquitous, they are normally ignored (hence the missing diagonal). If a contact
is found in the native structure, it appears above the diagonal, while if it is
a non-native contact, it is shown below the diagonal. In our previous example,
if (10,50) were a native contact we would place the dot at x=10, y=50 but if it
were a non-native contact we would place it at x=50, y=10 instead. Lastly, it
would be cumbersome and confusing to look at contacts for each individual
generation, so instead we average the contacts over a number of frames.
Partition 1 shows the contacts from the first quarter of the set of generation,
partition 2 shows the second quarter, and so on. The legend indicates the
colouring scheme which represents the occupancy of each contact during that
quarter. For example, if there were a total of 50 generations, partition 1 would
show the average from the first 12 - a contact with an occupancy of 25% in
partition 1 would have been present in 3 of the first 12 generations. In these
plots a maximum contact distance of 4.6A was used, except for aliphatic carbons
which were allowed to be up to 5.4A away.

Ideally as the protein folds, it should acquire more and more native contacts,
and they should become more persistent (higher occupancy). However, at RMSDs
larger than about 6A, it is generally unlikely that a sturcture would have many
native-like contacts. We do not expect to see a lot of native contacts until we
are able to reach structures with RMSDs to native in this neighbourhood.

PDB format files
----------------
The dominant file format for protein structures amongst structural biologists
today is still primarily the PDB file.  In order to obtain files in the MMDB
file format, it is easiest to simply download them from the MMDB datbase
(address given above).  The MMDB contains all structures on the PDB except for
theoretical models.


