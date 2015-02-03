This document requires revision. CWVH Jun 5 2012
Please refer to the updated file "README_TraDES_R_Analysis_Package.txt"



What is computed in the TraDES Log Files?

IMPORTANT: Do not try to read the log file with a text editor or the 
	columns will not line up, since they are tab separated.  Use
	a word processor or import it into a spreadsheet for further
	analysis.

Structure Summary
-----------------

When trades terminates, you will be given a summary of the 
random structure you have just generated, in the form of a
log file (as well as screen output), with the following 
information:

Time taken: from placement of the first residue to the final
	residue

Tries (total): the total number of attempts to place a residue.
	The minimum possible value is the length of your protein,
	but this is usually longer due to backtracking.  This
	is further broken up into the next three entries

No backbone solution tries: the number of tries to place a
	residue which failed because the backbone N and C
	peptide bond atoms could not be placed without
	introducing large errors into the structure

Crashes: the number of tries to place a residue which failed due
	to a van der Waals collision

Tries due to Distance Constraints (# Violations): if a constraint
	file was used with the walk, this is the number of atom
	placements which failed because they would violate a
	distance constraint.  Following this, in parentheses are
	the number of constraints actually violated in the final
	structure

	Note that the above three entries approximately sum to
	total tries.  The discrepancy is due to not counting the
	successful tries (!)

Sequence length: the length of your protein sequence in AA

Rgyr: the radius of gyration of the protein structure, in Angstroms

NCTerm distance: the distance, in Angstroms, from the N-terminal
	alpha carbon to the C-terminal alpha carbon

Rn: Flory's characteristic or "normalized" radius of gyration,
	calculated as Rgyr^2/(N*l^2) where Rgyr is radius
	of gyration, N is length of protein (in AA) and l is
	3.81A, the distance between C-alphas (consecutive
	polymer units)

Cn: Flory's characteristic or "normalized" end-to-end polymer
	distance, calculated as Rn above but replacing Rgyr with
	NCTerm distance.

	For large proteins, Cn/Rn approaches 6.0

Surface Accessibility: the approximate surface area of the
	structure, calculated with DSSP, in square Angstroms

Exposed Hydrophobics: the surface area of hydrophobic residues
	only, in square Angstroms

# Helical Residues (DSSP): directly from DSSP output 'H'

# Extended Residues (DSSP): directly from DSSP output 'E' or 'B'

In addition, to this information, the log file will also contain:

A header giving the file names and folding conditions, sequence,
system information and date.

Structure number: so you can easily look at the structure if you
	saved it

Hydrophobics Radius of Gyration: the radius of gyration of the
	hydrophobic residues in the structure, in Angstroms

# Extended Residues (CA-CA dist.): number of residues in the extended
	conformation based solely on geometric CA-CA distances rather
	than hydrogen bonding patterns

Zhang Potential: an atom-based statistical potential giving an
	idea of how many favourable contacts are made; see

	Zhang, C. et al.  (1997)  Determination of Atomic Desolvation 
	Energies From the Structures of Crystallized Proteins.  
	J. Mol. Biol.  267, 707-726.

#Bryant-Lawrence Potential: a residue-based threading potential.

	Bryant, S.H. and Lawrence, C.E.  (1993)  An Empirical Energy
	Function for Threading Protein Sequence Through the Folding 
	Motif.  Proteins.  16, 92-112.

  With both scoring functions as implemented here, local contacts 
  are ignored (i.e. contributions from residues less than 4 residues
  or 1 helix turn apart are not counted (Bryant & Crease))

Crease Energy: an experimental energy measure we are working with, as
	an extension to the Bryant-Lawrence potential. It is more
        sensitive than the Bryant-Lawrence potential. Unpublished data.

Voronoi Potential: an atom-contact based scoring function
       
        McConkey, B.J., Sobolev, V. and Edelman, M. (2003). Discrimination
        of native protein structures using atom-atom contact scoring.
        Proc. Natl Acad. Sci. U.S.A. 100:3215-3220 

	
	
	