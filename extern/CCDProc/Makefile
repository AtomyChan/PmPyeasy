#
# Makefile for the CCDProc routines
#   ccdproc, imcombine, and normalize
#
# Edit to set flags as appropriate for your installation.
#
# You need the path to the HEASARC CFITSIO libraries.  You should
# only need to change CCDBin if you are doing an installation where
# you wish the binary executable files to go somewhere other than
# right here (.).
#
# "make install" tag is provided to clean up once you've tested it.
# If you set CCDBin to something other than the current directory, 
# you will need to uncomment the part of install: that moves the 
# binaries to the CCDBin directory.  install cleans up the .o files.
#
# Type "make clean" to remove extranous/old .o files and executables
# before making the first time.
#
# R. Pogge, OSU Astronomy Dept.
# pogge@astronomy.ohio-state.edu
# 1999 May 27
#
# Last modification: 2005 Jan 21 [rwp/osu]
#
###########################################################################

# path to the CFITSIO libraries (no terminal /)

#CFITSIO=/home/ping/pmpy/extern/cfitsio-4.0.0
CFITSIO=/usr/local/lib

# path to the bin directory.  Default is ., see install: below for
# other related bits

CCDBin=/data/microfun/0/microfun/bin

# Sun cc
#CC=cc
#CFLAGS= -g $(CINCLUDES)
#COPTFLAGS= -O $(CINCLUDES)

# gcc
CC=g++
CINCLUDES= -I$(CFITSIO)
# CFLAGS= -g -traditional -W -Wunused $(CINCLUDES) # noisy version
CFLAGS= $(CINCLUDES)
COPTFLAGS= -O $(CFLAGS) 

LDFLAGS=
CLIBS= -L$(CFITSIO) -lcfitsio -lm

#
# library routines
#

OBJS= gparam.o fit2.o Utils.o

#
# compile all executables
#
all: ccdproc imcombine normalize flt2int dicer imstat expix immath

#
# individual routines
#
ccdproc: $(OBJS) ccdproc.o CCDProc.h
	$(CC) $(LDFLAGS) -o ccdproc ccdproc.o $(OBJS) $(CLIBS)
	strip ccdproc

imcombine: $(OBJS) imcombine.o CCDProc.h
	$(CC) $(LDFLAGS) -o imcombine imcombine.o $(OBJS) $(CLIBS)
	strip imcombine

normalize: $(OBJS) normalize.o CCDProc.h
	$(CC) $(LDFLAGS) -o normalize normalize.o $(OBJS) $(CLIBS)
	strip normalize

flt2int: $(OBJS) flt2int.o CCDProc.h
	$(CC) $(LDFLAGS) -o flt2int flt2int.o $(OBJS) $(CLIBS)
	strip flt2int

dicer: $(OBJS) dicer.o CCDProc.h
	$(CC) $(LDFLAGS) -o dicer dicer.o $(OBJS) $(CLIBS)
	strip dicer

imstat: $(OBJS) imstat.o CCDProc.h
	$(CC) $(LDFLAGS) -o imstat imstat.o $(OBJS) $(CLIBS)
	strip imstat

expix: $(OBJS) expix.o CCDProc.h
	$(CC) $(LDFLAGS) -o expix expix.o $(OBJS) $(CLIBS)
	strip expix 

immath: $(OBJS) immath.o CCDProc.h
	$(CC) $(LDFLAGS) -o immath immath.o $(OBJS) $(CLIBS)
	strip immath 

# 
# install the binaries & mop up.
#   uncomment the \mv's if CCDBin is other than . 
#
install:
	\rm *.o
	\mv ccdproc $(CCDBin)/ccdproc
	\mv imcombine $(CCDBin)/imcombine
	\mv normalize $(CCDBin)/normalize
	\mv flt2int $(CCDBin)/flt2int
	\mv imstat $(CCDBin)/imstat
	\mv expix  $(CCDBin)/expix 
	\mv immath  $(CCDBin)/immath
	\mv dicer  $(CCDBin)/dicer

#
# clean up after yourself
#
clean:
	\rm -f *.o ccdproc imcombine normalize flt2int dicer imstat expix immath
