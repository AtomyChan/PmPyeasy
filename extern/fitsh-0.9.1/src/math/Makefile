SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: dft fit intersec spline elliptic expint poly.o polyfit.o match.a splinefit.o

.PHONY: all clean fit elliptic intersec spline expint dft

###############################################################################

dft:
	$(MAKE) -C dft

elliptic:
	$(MAKE) -C elliptic

expint:
	$(MAKE) -C expint

fit:
	$(MAKE) -C fit

intersec:
	$(MAKE) -C intersec

spline:
	$(MAKE) -C spline

###############################################################################

tpoint.o: tpoint.c tpoint.h
	$(CC) $(CFLAGS) -c tpoint.c

cpmatch.o: cpmatch.c cpmatch.h
	$(CC) $(CFLAGS) -c cpmatch.c

spmatrix.o: spmatrix.c spmatrix.h
	$(CC) $(CFLAGS) -c spmatrix.c

delaunay.o: delaunay.c delaunay.h
	$(CC) $(CFLAGS) -c delaunay.c

poly.o: poly.c poly.h
	$(CC) $(CFLAGS) -c poly.c

polyfit.o: polyfit.c polyfit.h
	$(CC) $(CFLAGS) -c polyfit.c

splinefit.o: splinefit.c splinefit.h
	$(CC) $(CFLAGS) -c splinefit.c

trimatch.o: trimatch.c trimatch.h
	$(CC) $(CFLAGS) -c trimatch.c

match.a: tpoint.o cpmatch.o spmatrix.o delaunay.o trimatch.o
	$(AR) src match.a tpoint.o cpmatch.o spmatrix.o delaunay.o trimatch.o

###############################################################################

clean:
	$(MAKE) -C dft clean
	$(MAKE) -C elliptic clean
	$(MAKE) -C expint clean
	$(MAKE) -C fit clean
	$(MAKE) -C intersec clean
	$(MAKE) -C spline clean
	rm -f *.o *.a
