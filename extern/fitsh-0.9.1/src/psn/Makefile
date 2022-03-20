SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

INC=../../include

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result -I$(INC)

.PHONY: all clean

all: psn-general.o psn-general-ds.o

psn-general.o: psn-general.c psn-general.h
	$(CC) $(CFLAGS) -c psn-general.c

psn-general-ds.o: psn-general-ds.c psn-general-ds.h psn-general.h
	$(CC) $(CFLAGS) -c psn-general-ds.c

clean:
	rm -f *.o *.a $(TARGETS)


