To Debug with GDB or Valgrind.
CWV Hogue June 2012  


Unpack the entire archive into a new location.

BEFORE YOU RUN THE BUILD SCRIPT. Changes need to be made to
1. your platform-dependent makefile include that gets renamed to "ncbi.mk" 
2. the trades.mk optimization file that tweaks optimizations for TraDES executables
3. the "make.target" file that builds the executable and by default strips out debugging info.
4. the required runtime files need to be copied into the src directory for testing.


---------------------------------
PART 1
Check to see what platform-dependent makefile is listed in
your build script:

e.g. If using CentoS 6 x64 your build script is:
./build_script/build_centos6_x64.sh

Use grep:
grep hoguemake ./build_script/build_centos6_x64.sh

and see that it uses:
cp ./ncbi/hoguemake/linux_CentOS_RHEL_5_5_x64.ncbi.mk ./ncbi/platform/linux.ncbi.mk



Edit this platform dependent makefile in its location:
./ncbi/hoguemake/linux_CentOS_RHEL_5_5_x64.ncbi.mk 

ADD -g to NCBI_CFLAGS like so:

NCBI_CFLAGS1 = -c -g 

CHANGE NCBI_LDFLAGS1 and NCBI_OPTFLAG like so:

NCBI_LDFLAGS = -O0
NCBI_OPTFLAG = -O0


Save it in ./ncbi/hoguemake 

This sets the entire library build tree (excepting precompiled objects) to compiled in debug mode.

-------------------------------------
PART 2

Edit the file 
./nusi/TraDES/makefiles/trades.mk

(this unsets a lot of compiler optimization flags and replaces them with the -g symbol)

Change

LDFLAGS1 = $(NCBI_LDFLAGS1)
OPTFLAG = -O2 $(EXTRAOPT)

To the debug -g version:

LDFLAGS1 = -g $(NCBI_LDFLAGS1)
OPTFLAG = -g 



-----------------------------------
PART 3


THEN you need to change the executable makefiles, for example to debug str2trj
edit this file: 
./TraDES/makefile/make.str2trj

COMMENT OUT the following lines:
 	strip $(EXE)
 	mv $(EXE) $(BUILDDIR)

like so:

#	strip $(EXE)
#	mv $(EXE) $(BUILDDIR)

Save it.



-------------------------------
PART 4

Last thing - you will need the data files to debug executables:
copy the files in ./nusi/common to ./nusi/TraDES/src


---------------------

FINALLY - run the build script as usual - debuggable executables will
be placed in ./nusi/TraDES/src


To rebuild an executable you are working on 
from within the ./nusi/Trades/ directory (where the makefiles are moved to)
for an example makefile edited as above issue
make -f make.str2trj clean
make -f make.str2trj

If you want to debug a different executable after running the build script
without running the entire script again Apply step 3 above to make.target 
and manually copy it to the ./nusi/TraDES directory before running it. 

Refer to the GNU Debugger (GDB) and the ValGrind memory checker for more information.









 
