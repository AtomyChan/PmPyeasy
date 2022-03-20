SHELL = /bin/sh

CC=gcc
AR=ar
LD=ld
RANLIB=ranlib

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FITSH_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

DLEXT=so
DLSWC=-shared

INC=../include

.PHONY: all clean

all: librandom.a

random.o: random.c $(INC)/random/random.h
	$(CC) $(CFLAGS) -I$(INC) -c random.c

librandom.a: random.o
	$(AR) src librandom.a random.o
	$(RANLIB) librandom.a

librandom.$(DLEXT): random.o
	$(LD) -o librandom.$(DLEXT) $(DLSWC) random.o
	chmod -x librandom.$(DLEXT)

clean:
	rm -f *.o *.a *.$(DLEXT)

