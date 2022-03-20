SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: iof.o tokenize.o scanarg.o format.o

.PHONY: all clean

iof.o: iof.c iof.h
	$(CC) $(CFLAGS) -c iof.c

scanarg.o: scanarg.c scanarg.h
	$(CC) $(CFLAGS) -c scanarg.c

tokenize.o: tokenize.c tokenize.h
	$(CC) $(CFLAGS) -c tokenize.c

format.o: format.c format.h
	$(CC) $(CFLAGS) -c format.c

clean:
	rm -f *.o *.a
