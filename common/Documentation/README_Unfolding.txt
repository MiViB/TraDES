


Unfolding with str2trj (see unfold-it below for even more unfolding magic!)
----------------

Simply give, in order, the name of the Asn.1 file containing the structure
to be unfolded, the temperature at which to simulate unfolding (-t degrees 
Kelvin) and the timestep (-d in femtoseconds) (10e-15) to determine how far it 
will move.  

You can try this example with the Asn.1 version of the 153L structure provided.

./str2trj -f mmdb_153L.cn3 -t 400 -d 50 -r 1

./trades -f mmdb_153L -p T

The trades engine will keep working on it till it produces a file, which
you can see in Cn3D, or also in a PDB viewer, as -pT outputs PDB files.

The -r parameter "preserves" the sidechain dihedral (chi) rotamer angles 
along with the backbone conformation, so that when the structure is rebuilt 
with the trades engine, it will attempt to use all the same rotamers as in 
the original structure, except where forbidden by collisions.  
Without this option, all sidechains are placed using the
bacbone-dependent rotamer library of Dunbrack et al.


Typically each residue can move a single degree at 298K in about 44fs, and
they move faster at higher temperatures, slower at lower temperatures.

Note that str2trj only creates the trajectory distribution file, you
must still run trades to get a possible structure.  Each run of trades
will produce a possible structure after the molecule has moved for the
given amount of time at the given temperature, starting from the crystal
structure given each time (in vaccuo, with no external forces).



str2trj also allows you you to use a model number other than model 1 
in a crystal structure for making a .trj file, or using a chain other 
than the first in a structure file.  

Try with the included file mmdb_1mcp.cn3

./strSummary -f mmdb_1mcp.cn3

---> output of FASTA shows there are 2 chains L and H in the Immunoglobin FAB fragment


The following command extracts the first 100 residues of the L chain for sampling
with preserved rotamers, at 300K over 50 time steps, and makes a file 1mcpL_N100.trj

./str2trj -f mmdb_1mcp.cn3 -c L -s 1 -e 100 -r 1 -t 300 -d 50 -o 1mcpL_N100

./trades -f 1mcpL_N100 

Also, a range of residues (other than the full range of a protein chain) 
may be specified if only 1 domain of a protein, for example, is to be modelled.



Using the unfold-it.sh script (Unix only)

unfold-it is a "Wrapper" around str2trj and trades.

It can unfold 
a) small single-domain folded proteins (150aa or less)
b) it can be used to do dynamics on an intrinsic disordered protein (big)

It cannot unfold larger folded proteins without excessive atomic collisions, and modifications to the str2trj program are needed to fix this.

How it works:
It carries out a full unfolding of a protein by repeatedly calling 
-str2trj to make the next trajectory distribution file
-trades to sample the structures
-salign to align the structures with the previous one
-bin2prt, concatmodels and prt2bin to join the unfolded structure into 
the finished "_movie.val" file.

Rename your protein to a *.val file before using this script.

After running unfold-it you should have a complete protein unfolding "movie" 
which you can view in Cn3D and analyze with the analyzeMovie program.  
See README_AnalyzeMovie.txt for details and examples.

Note that unfold-it script only works on the first chain in a structure, however
it can be modified in the first step to extract a segment by hand.
 

