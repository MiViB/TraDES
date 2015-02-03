Using seq2trj - Updated June 5 2012 CWVHogue


Parameters:
--------------------------------------------------------------------------

-f : the name of the file containing the FASTA or raw sequence
	info (1-letter AA abbreviations) to be used for sampling.

-o : the name to give the trajectory file output by the program.
	.trj will automatically be appended to the filename
	you specify, which defaults to protein(.trj).  This file
	contains all the probability distributions used by
	foldtraj to generate your random structure, plus some
	other useful information and parameters.
	
	Other output files include a (.ss) file if str2trj used
	some method of combining alpha beta and coil probabilities.
	By default this uses the GOR secondary structre prediction 
	method, which we find does not overpredict secondary structure
	in intrinsic disordered proteins.
	
	The last output file is a (.ara) which is a single number
	representing the area of conformational space available at each
	amino acid. It is computed from the residue trajectory distribution 
	after it is made with whatever -t or -c or -b combination you 
	specify.  The integer value is the number of cells (out of 400x400)
	that are below 50% of the maximum probability value. 
	(Full Area at Half Max - a 3d version of the usual bandwidth measure FWHM)
	-t 1 - Uniform sampling sets these all to maximum 160000. 
	Alpha-helical signals under -t 3 will produce minimal values.

-s :  You can provide a secondary structure prediction file to be used
	for 3-state sampling method. See the /Samples directory 
	The input file has five columns,
	and one row per residue in the protein.  The values,
	from left to right, are residue ID, predicted structure
	(H=helix, E=extended/sheet, C=coil), % probability of
	Helix, % probability of sheet and % probability of coil.
	The last three columns should always sum to 100% aside
	from rounding error

-c : use this option along with a full path/filename to impose
	X-PLOR/CNS-like constraints on the random walk.  For a
	sample constraint file, please see 'sample.constr' included
	in the /Samples for details.  This file explains how to
	enter distance constraints and gives a few examples.  Other
	restraints (such as dihedrals) may be supported in a future
	release. This is very useful for disulphides.

-t : you may choose from four possible trajectory distributions construction methods
	(probability distributions for random walk):
	1: Uniform - all regions of space equally likely,
		flat probability distribution function (PDF)
		In testing, this is a waste of time but provided
		so that you can see for yourself that proteins 
		prefer certain Ramachandran angles.
	2: Amino Acid Based - a different PDF is chosen for each
		of the 20 possible amino acid types, based on a
		statistical analysis of a non-redundant set of the
		PDB
	3: Secondary structure - the GOR (or input .ss file) 1-state prediction is used
		in conjunction with amino acid type to limit the
		domain of the PDF; regions outside the predicted
		secondary structure type are forbidden.
	4: GOR Method - the 3-state GOR prediction (or input .ss file) is convolved with
		the 3 regions of conformational space, again 
		separated by amino acid type, to produce a biased
		PDF.  For example, if a given Ala residue is
		predicted to be 75% helix, 20% coil, 5% sheet, we
		take 75% of the normalized Ala helix PDF, 20% of the
		coil PDF and 5% of the sheet PDF, and add them up.
		This method leads to the most realistic results.

-c : this sets the trajectory distribution to use 100% coil sampling
-b : this sets the trajectory distribution to use 100% beta sampling
	NOTE: these override any of the -t settings. 
	and a combination of -c T -b T sets 50% coil, 50% beta
	
	Use -c for sampling IDP if you know there is no stable secondary structure
	in it, otherwise use the default -t 4 
		
-v : this outputs the entire Nx160000 trajectory graph as a .csv file, if you need 
	the computed data in plain text format. 
	
-q : this sets to quiet operation.
 
