SHELL = /bin/sh

CC=gcc
AR=ar
LD=ld
RANLIB=ranlib

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FITSH_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

DLEXT=so
DLSWC=-shared

INC=../include

.PHONY: all clean

all: libpsn.a

psn.o: psn.c $(INC)/psn/psn.h
	$(CC) $(CFLAGS) -I$(INC) -c psn.c

libpsn.a: psn.o
	$(AR) src libpsn.a psn.o
	$(RANLIB) libpsn.a

libpsn.$(DLEXT): psn.o
	$(LD) -o libpsn.$(DLEXT) $(DLSWC) psn.o
	chmod -x libpsn.$(DLEXT)

clean:
	rm -f *.o *.a *.$(DLEXT)

