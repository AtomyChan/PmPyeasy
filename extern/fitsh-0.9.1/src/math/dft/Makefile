SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: pbfft.o 

.PHONY: all clean

pbfft.o: pbfft.c pbfft.h
	$(CC) $(CFLAGS) -c pbfft.c

clean:
	rm -f *.o *.a
