SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FI_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result

DLEXT=dylib
DLSWC=-dylib

INC=../include

all: libfits.a

FITS_MODULES=\
	fits-bintable.o \
	fits-common.o \
	fits-draw.o \
	fits-header.o \
	fits-image.o \
	fits-table.o \
	fits-ui.o \
	fits-core.o

.PHONY: all clean

fits-common.o: fits-common.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-common.c

fits-header.o: fits-header.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-header.c

fits-image.o: fits-image.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-image.c

fits-table.o: fits-table.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-table.c

fits-bintable.o: fits-bintable.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-bintable.c

fits-draw.o: fits-draw.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-draw.c

fits-ui.o: fits-ui.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-ui.c

fits-core.o: fits-core.c fits-common.h $(INC)/fits/fits.h
	$(CC) $(CFLAGS) -I$(INC) -c fits-core.c

libfits.a: $(FITS_MODULES)
	$(AR) src libfits.a $(FITS_MODULES)

libfits.$(DLEXT): $(FITS_MODULES)
	$(CC) $(CFLAGS) -o libfits.$(DLEXT) $(DLSWC) $(FITS_MODULES)
	chmod -x libfits.$(DLEXT)

clean:
	rm -f *.o *.a *.$(DLEXT)
