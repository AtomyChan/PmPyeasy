SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: sort.o multiindex.o

.PHONY: all clean

sort.o: sort.c sort.h
	$(CC) $(CFLAGS) -c sort.c 

multiindex.o: multiindex.c multiindex.h
	$(CC) $(CFLAGS) -c multiindex.c

clean:
	rm -f *.o
