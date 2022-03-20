SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

all: link.a floodfill.o

.PHONY: all clean

link.a: linkpoint.o linkblock.o 
	$(AR) src link.a linkpoint.o linkblock.o 

floodfill.o: floodfill.c floodfill.h
	$(CC) $(CFLAGS) -c floodfill.c

linkpoint.o: linkpoint.c linkpoint.h link.h
	$(CC) $(CFLAGS) -c linkpoint.c

linkblock.o: linkblock.c linkblock.h link.h
	$(CC) $(CFLAGS) -c linkblock.c

clean:
	rm -f *.o *.a
