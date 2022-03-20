SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld
RANLIB=ranlib

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FITSH_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: expint.o 

.PHONY: all clean

expint.o: expint.c expint.h
	$(CC) $(CFLAGS) -c expint.c

clean:
	rm -f *.o *.a
