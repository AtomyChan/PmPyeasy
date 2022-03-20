SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: lmfit.o downhill.o

.PHONY: all clean

lmfit.o: lmfit.c lmfit.h
	$(CC) $(CFLAGS) -c lmfit.c

downhill.o: downhill.c downhill.h
	$(CC) $(CFLAGS) -c downhill.c

clean:
	rm -f *.o *.a
