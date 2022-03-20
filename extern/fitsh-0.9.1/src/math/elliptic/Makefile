SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

.PHONY: all clean

TARGETS=elliptic.o ntiq.o

all: $(TARGETS)

elliptic.o: elliptic.c elliptic.h
	$(CC) $(CFLAGS) -c elliptic.c

ntiq.o: ntiq.c ntiq.h
	$(CC) $(CFLAGS) -c ntiq.c

clean:
	rm -f *.o $(TARGETS)

