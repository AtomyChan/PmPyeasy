SHELL = /bin/sh

CC=gcc
AR=ar
LD=ld

INC=../include

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result -I$(INC)

all: libastro.a

.PHONY: all clean

ASTRO_MODULES= \
	astro.o \
	cmath.o \
	dtc.o \
	earth.o \
	easpec.o \
	planets.o \
	output.o 

###############################################################################

libastro.a: $(ASTRO_MODULES)
	$(AR) src libastro.a $(ASTRO_MODULES)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

astro.o: astro.c $(INC)/astro/astro.h
	$(CC) $(CFLAGS) -c astro.c

cmath.o: cmath.c $(INC)/astro/cmath.h
	$(CC) $(CFLAGS) -c cmath.c

dtc.o: dtc.c $(INC)/astro/dtc.h 
	$(CC) $(CFLAGS) -c dtc.c

earth.o: earth.c $(INC)/astro/earth.h $(INC)/astro/astro.h
	$(CC) $(CFLAGS) -c earth.c

easpec.o: easpec.c $(INC)/astro/easpec.h $(INC)/astro/astro.h
	$(CC) $(CFLAGS) -c easpec.c

planets.o: planets.c $(INC)/astro/planets.h $(INC)/astro/astro.h
	$(CC) $(CFLAGS) -c planets.c

output.o: output.c $(INC)/astro/output.h
	$(CC) $(CFLAGS) -c output.c

###############################################################################

clean:
	rm -f *.o *.a

###############################################################################
