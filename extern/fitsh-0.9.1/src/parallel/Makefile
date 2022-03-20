SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: parallel.a

.PHONY: all clean

none.o: none.c parallel.h
	$(CC) $(CFLAGS) -c none.c 

ipc.o: ipc.c parallel.h
	$(CC) $(CFLAGS) -c ipc.c

parallel.o: parallel.c parallel.h
	$(CC) $(CFLAGS) -c parallel.c

parallel.a: none.o ipc.o parallel.o 
	$(AR) src parallel.a none.o ipc.o parallel.o

clean:
	rm -f *.o *.a
