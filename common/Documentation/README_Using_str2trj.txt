Using str2trj



Examples:

To obtain a trajectory distribution for the reconstruction of a structure
in the TraDES_2/Samples directory, copy 

1YU5.val 

into the TraDES_2 folder, then run:


str2trj -f 1YU5 -o 1YU5

which will produce 1YU5.trj as its output file.   

Then:

trades -f 1YU5 -b 10 -p T

will make 10 samples of near-neighbors of this structure.


As another example, to reconstruct the crystal structure 1YU5.val and add
noise with standard deviation 4 degrees in the theta direction and 6
degrees in the phi direction, and an output file called
1YU5_4x6y.trj, saving the sidechain rotamer angles in the output file as well,

str2trj -f 1YU5 -x 4 -y 6 -o 1YU5_4x6y -r 1

Then

trades -f 1YU5_4x6y -b 10 -p T

((((TO ADD back to str2trj)))) 
The increased peak height will serve to provide added precision in the
tails of the gaussian since trajectory distributions heights are integers,
so that once we fall below a height of 1, everything is treated as zero
probability.

 
Once you have the .trj file created, running trades will
attempt to reconstruct the protein backbone and sidechains.  

Each time it is run, a different structure will result, but all 
will be similar to the given crystal structure.  The RMSD will increase as 
bigger values are fed to the -x and -y options here.  

For very large -x and/or -y values, the
structures will no longer be any more 
similar that two random structures.




----------More about the arguments...----------------



str2trj: ASN.1 protein 3D structure converted to TraDES *.trj file
Use for protein near-neighbor or decoy sampling.

SEE ALSO README_Unfolding.txt for more information about using str2trj
for protein unfolding simulations.

   arguments:

  -f  REQUIRED 
Input 3D Asn.1 VAL File Name PREFIX (No extension) [File In]

  -c  PDB Chain Name (default = first one found) [String]  Optional
  -m  Model Number (default = 1) [Integer]  Optional
    default = 1
  -s  Start Residue (default = 1) [Integer]  Optional
    default = 1
  -e  End Residue (default = 0 represents last residue). [Integer]  Optional
    default = 0

	These settings allow you to select a single chain, and/or subsequence of that chain.
        The model number allows you to select one model from a movie or nmr structure.

  -r  Save Sidechain Rotamer angles
      (Default 0=none; 1=all, 2=buried) [Integer]  Optional
    default = 0
    range from 0 to 2

	This option tells str2Trj to "preserve" the sidechain dihedral (chi)
  	angles along with the backbone conformation, so that when the structure
  	is rebuilt with foldtraj, it will attempt to use all the same rotamers
  	as in the original structure, except where forbidden by collisions.
  	Without this option, all sidechains are placed using the bacbone-
  	dependent rotamer library of Dunbrack et al.

  -o  REQUIRED 
Output TRJ File Name Prefix (No extension):  [File Out]

        This is the output file name, to which .trj will automatically be
	appended.  

  -x  Peak Width MODE - Standard deviation for x (Degrees) [Real]  Optional
    default = 0
    range from 0 to 45.0

        Adds gaussian noise to the crystal structure, to
	allow the alpha carbons to wander somewhat from their orientations
	in the crystal structure then give a standard deviation in degrees
	here, for the x dimension (theta direction in trajectory space).
	A gaussian is placed in trajectory space centred about the true
	location of the residue, rather than a single sharp peak, and with
	given standard deviation. By default no noise is added.



  -y  Peak Width MODE - Standard deviation for y (Degrees)
PEAK WIDTH MODE Example: -x 5 -y 5
 [Real]  Optional
    default = 0
    range from 0 to 22.0
        
        As with -x, this adds noise in the y or "cos(phi)" direction of
	trajectory space.  Again, use units of degrees phi, and the
	default is 0.  In most cases when noise is added, you will want
	the -x and -y standard deviations to be equal.


  -q  Quiet Operation (T/F) [T/F]  Optional
    default = TRUE





