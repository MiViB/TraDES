########################################################################
#       Generic command-line makefile for NCBI applications
#	Adapted for compilation of TRADES projects
#
#   This assumes the following variables are set in a file called "ncbi.mk"
#     which is included by this one. In this case it resides in this directory
#
#    NCBI_CC = default compiler to use
#    NCBI_CFLAGS1 = default compiler switches (not optimization)
#    NCBI_LDFLAGS1 = default compiler switches when linking (not optimization)
#    NCBI_OPTFLAG = default optimization (-O)
#    NCBI_INCDIR = default toolkit include paths
#    NCBI_LIBDIR = default toolkit library paths
#    NCBI_ALTLIB  = debugging toolkit library paths
#    NCBI_VIBFLAG = additional include paths/defines for Vibrant apps
#    NCBI_VIBLIBS  = additional libraries needed to link Vibrant apps
#    NCBI_OTHERLIBS = additional libraries needed to link the system
#
#   This enables platform independent makefiles for ncbi toolkit apps. You
#   can change any default behavior by setting the variable on the command
#   line.
#
#   Use gcc complier.   "make NCBI_CC=gcc" or "make CC=gcc"
#   Change optimization to debug
#                "make NCBI_OPTFLAG="-g" NCBI_LIBPATH=$NCBI_ALTLIB"
#                or
#                "make OPTFLAG="-g" LIBPATH=$NCBI_ALTLIB"
#
#   You can also change the behavior permanently in the makefile by changing
#    the assignment of the internal variable to the one from the include
#    file:
#
#   Use gcc compiler.
#    Change   CC=$(NCBI_CC)
#    To       CC=gcc
#    in the makefile itself
#
#   Make debugging version
#    OPTFLAG=-g
#    LIBPATH=$(NCBI_ALTDIR)
#
#   You need to specify the EXE and SRC symbols
#
#######################################################################
#
# machine dependent paths, etc. for trajectory graph realted programs
#
#######################################################################

include ncbi.mk
include codebase.mk

BZINC = ../bzip2
DSSPINC = ../dssp
BZLIB = ../bzip2
DSSPLIB = ../dssp
SRCDIR = ./src
MMDB_SRCDIR = ./src/libmmdbtrj
FOLD_SRCDIR = ./src/libfoldtrj
MAKE_SRCDIR = ./src/libmaketrj
TRADESLIBS = ./lib
TRADESINC = ./include
TSLRILIB = ../slri/lib
TSLRIINC = ../slri/include
PUT_VERSION = -DFOLDTRAJ_VERSION=\"`date +%Y.%m.%d`\"
TRAJLIBOBJS = $(SRCDIR)/vecmath.o $(SRCDIR)/rotate.o $(SRCDIR)/bbox.o $(SRCDIR)/Rplot.o $(SRCDIR)/rotlib.o $(SRCDIR)/potential.o $(SRCDIR)/vscore.o $(SRCDIR)/crease.o $(SRCDIR)/dft.o $(SRCDIR)/buildit.o $(SRCDIR)/gil_gor.o $(SRCDIR)/extaaseq.o
FOLDTRAJLIBOBJS = $(SRCDIR)/randwalk.o $(SRCDIR)/GetTrjAngle.o $(SRCDIR)/trajtools.o $(SRCDIR)/newasn.o $(SRCDIR)/objtraj.o
MAKETRAJLIBOBJS = $(SRCDIR)/maketg.o $(SRCDIR)/vismaketrj.o
GIFSRC = $(NCBI)/vibrant/picture.c $(NCBI)/vibrant/mapping.c $(NCBI)/vibrant/drawing.c
TRAJLIB = libmmdbtraj.a
FOLDTRAJLIB = libfoldtraj.a
MAKETRAJLIB = libmaketrj.a
TRAJLIBNAME = $(TRADESLIBS)/$(TRAJLIB)
FOLDTRAJLIBNAME = $(TRADESLIBS)/$(FOLDTRAJLIB)
MAKETRAJLIBNAME = $(TRADESLIBS)/$(MAKETRAJLIB)
INCDIRS = -I$(TRADESINC) -I$(BZINC) -I$(DSSPINC) -I$(CBINC) -I$(TSLRIINC)
LIBDIRS = -L$(TRADESLIBS) -L$(BZLIB) -L$(DSSPLIB) -L$(CBLIB) -L$(TSLRILIB)
BUILDDIR = ./build
EXTRAOPT = -ffast-math -finline-functions -Wall -Wconversion -Wnested-externs -Winline

#Aliased NCBI libraries - WARNING, these numbers may change!!!
LIB1 = -lncbi
LIB2 = -lncbiobj
LIB3 = -lncbicdr
LIB4 = -lvibrant
LIB5 = -lncbiacc
LIB6 = -lnetcli
LIB7 = -lnetentr
LIB8 = -lncbiNacc
LIB9 = -lncbiCacc
# LIB10 is reserved for NCSA socket library
LIB10 =
LIB11 = -lncbimla
LIB12 = -lncbitax
LIB13 = -lncbiid0
#LIB14 = -lncbibls0
LIB15 = -lnetcliE
LIB16 = -lnetcliS
LIB17 = -lnetcliES
LIB19 = -lncbispel
# LIB20 is for the NCBI desktop utilities
LIB20 = -lncbidesk
LIB21 = -lncbibls2
LIB22 = -lncbimmdb
LIB23 = -lncbitool
LIB24 = -lncbisugg
LIB25 = -lncbiwww
LIB26 = -lncbitax1
LIB27 = -lncbimsc1
LIB28 = -lvibgif
LIB29 = -lncbitxc1
LIB30 = -lncbicn3d
LIB31 = -lvibnet
LIB33 = -lsmartnet
LIB35 = -lssparse
LIB36 = -lnetblast
LIB38 = -lncbiprsn
LIB40 = -lncbitxc2
LIB41 = -lncbiid1
LIB42 = -lctutils
LIB43 = -losutils
LIB44 = -lidload
LIB45 = -lddvlib
LIB400 = -lvibrantOGL
LIB3000 = -lncbicn3dOGL

# TraDES required libraries
LIB95 = -ldssp
LIB98 = -lbz
LIB99 = -lmmdbtraj
LIB100 = -lfoldtraj
LIB101 = -ltslri
LIB104 = -lmaketrj
LIB199 = -lcb

MISCLIBS = -lcb 
GLUT_LIB = -lglut

#######################################################################
#
# default flags for compiling and loading
#
#######################################################################

CC = $(NCBI_CC)
CFLAGS1 = $(NCBI_CFLAGS1)
LDFLAGS1 = $(NCBI_LDFLAGS1)
OPTFLAG = -O2 $(EXTRAOPT)
OTHERLIBS = $(NCBI_OTHERLIBS)
VIBLIBS = $(NCBI_VIBLIBS)
VIBFLAG = $(NCBI_VIBFLAG)
INCPATH = $(NCBI_INCDIR) $(INCDIRS)
LIBPATH = $(NCBI_LIBDIR) $(LIBDIRS)
DEF = -DUSE_DSSP
