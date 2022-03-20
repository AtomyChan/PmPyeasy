SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: intersec.o intersec-cri.o

.PHONY: all clean

intersec.o: intersec.c intersec.h
	$(CC) $(CFLAGS) -c intersec.c

intersec-cri.o: intersec-cri.c intersec-cri.h
	$(CC) $(CFLAGS) -c intersec-cri.c

clean:
	rm -f *.o *.a 
