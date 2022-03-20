SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: bicubic.o biquad.o biquad-isc.o spline.o

.PHONY: all clean

bicubic.o: bicubic.c bicubic.h
	$(CC) $(CFLAGS) -c bicubic.c

biquad.o: biquad.c biquad.h
	$(CC) $(CFLAGS) -c biquad.c

biquad-isc.o: biquad-isc.c biquad-isc.h
	$(CC) $(CFLAGS) -c biquad-isc.c

spline.o: spline.c spline.h
	$(CC) $(CFLAGS) -c spline.c

clean:
	rm -f *.o *.a
