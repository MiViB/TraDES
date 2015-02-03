
Using trades engine

CWVH updated 5 Jun 2012

Start with str2trj or seq2trj and obtain a *.trj file.

A *.trj file is a compressed set of 3D probability distributions 
	in either Ramachandran Space, or Alpha Carbon space.
	See the original papers in the Documentation directory.

Trades simply needs one to start.

./seq2trj -f ala.faa -o test

./trades -f test

and it will make 1 structure.

Trades makes 3d structure files, a .log file for the run,
and if requested, ramachandran angle files in .csv format.

The log file contains extensive calculations, as described
in the README_output_and_analysis.txt file.







Details about the trades engine command-line parameters are:

trades 
--------------

-f : this is the name of the .trj input file created by InitTraj.
	 

-o :    the name of the output structure file; if none is given,
	the same base name as the input file is used.  Again, the
	appropriate extensions will be automatically appended, do
	not include these in the name you give.

-n :    This is used as the seed for the random 
	number generator.  If you wish to generate reproducible
	results, specify a seed.  The same structure will be
	generated whenever the same seed is used.  The default
	of 0 means to use the clock to get a truly random seed
	instead.  Note that due to differences in floating point
	rounding, using the same seed on different machine
	architechtures MAY result in different results.

-d : if this option is set, the input .trj file will be deleted
	once the program is done, otherwise, by default, it will
	remain on your machine.  This may be useful in scripts.
	Otherwise, the same input may be reused over and over
	for multiple folding "experiments" as it contains no
	random elements and will be the same each time for a
	given protein and trajectory distribution setting  

-p :    Output a .pdb structure file, named with
	the -o option above; use this to turn out PDB output.
        or use str2pdb to convert Asn.1 into PDB

-a : by default, only the MMDB ASN.1 .val structure file will be
	output.  If you require one, for use with Cn3D, use this
	option.  The filename is given in the -f option (see
	above). PDB files can be made with str2pdb at any time,
        or use the -P 

-b : use this to give the number of structures you want
	generated.  Structure files, if output, will be
	numbered sequentially, for example 'protein_0000001.pdb',
	'protein_0000002.pdb' and so on

-s : if you wish to start structure numbering at something other
	than 1 (e.g. continuing where you left off earlier) then
	specify here the initial value, otherwise, it starts at 1

-x : if you wish to give an "experiment name" here, this name will
	be included in all structure filenames output.  This is
	useful for distinguishing sets of structures generated 
	under different conditions or for different purposes

-c : if the native structure (in MMDB ASN.1 format, available at
	http://www.ncbi.nlm.nih.gov/Structure/) is known and in
	the current directory, provide it here and the RMSD to the
	native structure will be reported for all randomly
	generated structures.  Note that both structures must be
	atom-for-atom identical in order to obtain a meaningful
	RMSD

-h : if the native structure given in -c that you wish to align your
	random structures to is not the first molecule in the file,
	you must specify teh chain name here

-q : normally operation is in "quiet" mode.  For verbose
	operation, to see what is being done at all times, use
	this option


-F : Fragment file sampling - unused option at the moment

-r : Collect Ramachandran information for up to 30000 structures. 
	This option outputs up to 20 .csv files starting with the first 
	letter of the amino acid.

-v : Switch the default from binary to ASCII Asn.1 here.

-l : Number of water layers to add to the structure. 
	A Principal Components Analysis uses the water shell to calculate 
	the ellipsoid parameter, then to align the structure along 
	the X, Y, and Z axes.  This uses 1 layer of water
	in calculation, but you can add up to 10.

-w : default is false, but if you want solvated structures output, turn on.

-t : Stream progress - this is used to feed updates to R shells, but can be
	rigged for other systems to monitor trades jobs.

-k : this flag bypasses any PCA rotation of the structure after 
	PCA analysis. When set, the first alpha-carbon coordinate of the 
	structure will be 0,0,0. By default the protein is centered at 
	the PCA determined centre of its water layer at 0,0,0.

-z : this flag is used for unfolding mode with str2trj. It detects when a 
        trades job gets stuck a lot and slowly ramps up the bounciness 
        parameters. This is set FALSE by default for IDP sampling mode


