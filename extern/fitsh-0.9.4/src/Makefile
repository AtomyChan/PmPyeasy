SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld
RANLIB=ranlib

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FITSH_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result -I$(INC)

DEPTH=..

INC=$(DEPTH)/include

CFLARGEFILE=

DLEXT=so
DLSWC=-shared
DLDYN=-rdynamic
DLLIB=-ldl

TARGETS  = fiarith ficalib ficombine ficonv fiheader fiign fiinfo fiphot
TARGETS += firandom fistar fitrans 
TARGETS += grcollect grmatch grtrans gropt
TARGETS += lfit 
TARGETS += linear.so

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

DEP_FIARITH=\
	history.o \
	fitsmask.o \
	statistics.o \
	imgtrans.o \
	common.o \
	tensor.o \
	ui.o \
	longhelp.o

MOD_FIARITH=\
	io/iof.o \
	io/scanarg.o \
	psn/psn-general.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/spline/biquad.o \
	math/dft/pbfft.o

LIB_FIARITH=$(DEPTH)/libfits/libfits.a $(DEPTH)/libpsn/libpsn.a $(DEPTH)/librandom/librandom.a -lm

DEP_FICOMBINE=\
	fitsmask.o \
	statistics.o \
	tensor.o \
	common.o \
	ui.o \
	longhelp.o \
	history.o \
	combine.o

MOD_FICOMBINE=\
	io/iof.o \
	io/scanarg.o \
	math/fit/lmfit.o \
	math/poly.o

LIB_FICOMBINE=$(DEPTH)/libfits/libfits.a -lm

DEP_FICALIB=\
	fitsmask.o \
	statistics.o \
	tensor.o \
	common.o \
	ui.o \
	longhelp.o \
	history.o \
	str.o \
	combine.o \
	fbase.o

MOD_FICALIB=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/spline/spline.o \
	math/splinefit.o 

LIB_FICALIB=$(DEPTH)/libfits/libfits.a -lm

DEP_FICONV=\
	history.o \
	fitsmask.o \
	statistics.o \
	tensor.o \
	common.o \
	ui.o \
	longhelp.o \
	kernel.a

MOD_FICONV=\
	math/fit/lmfit.o \
	math/poly.o \
	math/spline/spline.o \
	math/spline/bicubic.o \
	math/spline/biquad.o \
	io/iof.o \
	io/tokenize.o \
	io/scanarg.o

LIB_FICONV=$(DEPTH)/libfits/libfits.a -lm

DEP_FIHEADER=\
	longhelp.o \
	str.o \
	ui.o

MOD_FIHEADER=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_FIHEADER=$(DEPTH)/libfits/libfits.a -lm

DEP_FIIGN=\
	fitsmask.o \
	statistics.o \
	longhelp.o \
	ui.o \
	maskdraw.o \
	common.o \
	history.o \
	tensor.o

MOD_FIIGN=\
	math/fit/lmfit.o \
	math/poly.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_FIIGN=$(DEPTH)/libfits/libfits.a -lm

DEP_FIINFO=\
	fitsmask.o \
	statistics.o \
	fiinfo-pnm.o \
	fiinfo-image.o \
	ui.o \
	longhelp.o \
	tensor.o

MOD_FIINFO=\
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/splinefit.o \
	math/spline/biquad.o \
	math/spline/bicubic.o \
	math/spline/spline.o  \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	link/link.a

LIB_FIINFO=$(DEPTH)/libfits/libfits.a -lm

DEP_FIPHOT=\
	fiphot-io.o \
	fitsmask.o \
	statistics.o \
	magnitude.o \
	common.o \
	longhelp.o \
	ui.o \
	kernel.a \
	apphot.o \
	tensor.o \
	weight-io.o \
	weight-gen.o \
	str.o \

MOD_FIPHOT=\
	math/intersec/intersec.o \
	math/intersec/intersec-cri.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	math/spline/spline.o \
	math/polygon.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_FIPHOT=$(DEPTH)/libfits/libfits.a -lm

DEP_FIRANDOM=\
	history.o \
	longhelp.o \
	ui.o \
	tensor.o \
	psf-io.o \
	star-draw.o \
	firandom-eval.o \
	magnitude.o \
	common.o \
	fitsmask.o \
	maskdraw.o \
	weight-star.o \
	weight-io.o

MOD_FIRANDOM=\
	math/fit/lmfit.o \
	math/poly.o \
	math/expint/expint.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	psn/psn-general.o \
	math/intersec/intersec-cri.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_FIRANDOM=$(DEPTH)/libfits/libfits.a $(DEPTH)/libpsn/libpsn.a $(DEPTH)/librandom/librandom.a -lm

DEP_FISTAR=\
	fistar-io.o \
	fitsmask.o \
	statistics.o \
	tensor.o \
	stars.a \
	ui.o \
	longhelp.o \
	common.o \
	imgtrans.o \
	background.o \
	magnitude.o \
	psf.a

MOD_FISTAR=\
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	math/expint/expint.o \
	math/intersec/intersec-cri.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	index/sort.o \
	index/multiindex.o \
	link/floodfill.o \
	link/link.a \
	math/delaunay.o

LIB_FISTAR=$(DEPTH)/libfits/libfits.a -lm

DEP_FITRANS=\
	fitsmask.o \
	fbase.o \
	statistics.o \
	transform.o \
	common.o \
	history.o \
	ui.o \
	longhelp.o \
	tensor.o \
	weight-io.o \
	weight-gen.o 

MOD_FITRANS=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	io/format.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	math/spline/bicubic.o \
	math/spline/spline.o \
	math/fit/lmfit.o \
	math/poly.o

LIB_FITRANS=$(DEPTH)/libfits/libfits.a -lm

DEP_GRCOLLECT=\
	cache.o \
	statistics.o \
	longhelp.o \
	ui.o \
	str.o

MOD_GRCOLLECT=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_GRCOLLECT=-lm

DEP_GRMATCH=\
	longhelp.o \
	ui.o \
	transform.o

MOD_GRMATCH=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/convexhull.o \
	math/match.a

LIB_GRMATCH=-lm

DEP_GRTRANS=\
	statistics.o \
	tensor.o \
	ui.o \
	longhelp.o \
	projection.o \
	transform.o

MOD_GRTRANS=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o

LIB_GRTRANS=-lm

DEP_GROPT=\
	optcalc.o \
	longhelp.o \
	ui.o 

MOD_GROPT=\
	io/tokenize.o \
	io/scanarg.o \
	math/polygon.o \
	math/matrixvector.o 

LIB_GROPT=-lm

DEP_LFIT=\
	str.o \
	xfunct.o \
	longhelp.o \
	ui.o \
	lfit-builtin.o \
	lfit-info.o

MOD_LFIT=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/dft/pbfft.o \
	math/fit/downhill.o \
	math/elliptic/elliptic.o \
	math/elliptic/ntiq.o \
	math/spline/spline.o

LIB_LFIT=$(DEPTH)/libpsn/libpsn.a $(DEPTH)/libastro/libastro.a $(DEPTH)/librandom/librandom.a -lm $(DLLIB)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

DEP_FIC_MPSTACK=\
	fitsmask.o \
	fbase.o \
	statistics.o \
	transform.o \
	common.o \
	history.o \
	ui.o \
	longhelp.o \
	tensor.o

MOD_FIC_MPSTACK=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	io/format.o \
	math/spline/bicubic.o \
	math/spline/spline.o \
	math/poly.o

LIB_FIC_MPSTACK=$(DEPTH)/libfits/libfits.a -lm

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

DIRS=index io link math psn

all: dirs
	$(MAKE) targets

targets: $(TARGETS)

dirs: $(DIRS)

.PHONY: all clean dirs $(DIRS) targets

###############################################################################

index:
	$(MAKE) -C index

io:
	$(MAKE) -C io

link:
	$(MAKE) -C link

math:
	$(MAKE) -C math

psn:
	$(MAKE) -C psn

###############################################################################

background.o: background.c background.h
	$(CC) $(CFLAGS) -c background.c

cache.o: cache.c cache.h
	$(CC) $(CFLAGS) -c cache.c

combine.o: combine.c combine.h
	$(CC) $(CFLAGS) -c combine.c

fbase.o: fbase.c fbase.h
	$(CC) $(CFLAGS) -c fbase.c

imgtrans.o: imgtrans.c imgtrans.h
	$(CC) $(CFLAGS) -c imgtrans.c

magnitude.o: magnitude.c magnitude.h
	$(CC) $(CFLAGS) -c magnitude.c

pselect.o: pselect.c pselect.h
	$(CC) $(CFLAGS) -c pselect.c

statistics.o: statistics.c statistics.h
	$(CC) $(CFLAGS) -c statistics.c

str.o: str.c str.h
	$(CC) $(CFLAGS) -c str.c

tensor.o: tensor.c tensor.h
	$(CC) $(CFLAGS) -c tensor.c

transform.o: transform.c transform.h
	$(CC) $(CFLAGS) -c transform.c

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

fitsmask.o: fitsmask.c fitsmask.h
	$(CC) $(CFLAGS) -c fitsmask.c

apphot.o: apphot.c apphot.h fitsh.h
	$(CC) $(CFLAGS) -c apphot.c

common.o: common.c common.h fitsh.h
	$(CC) $(CFLAGS) -c common.c

psf-base.o: psf-base.c psf-base.h psf.h
	$(CC) $(CFLAGS) -c psf-base.c

psf-io.o: psf-io.c psf-io.h psf.h
	$(CC) $(CFLAGS) -c psf-io.c

psf-determine.o: psf-determine.c psf-determine.h psf.h stars.h
	$(CC) $(CFLAGS) -c psf-determine.c

ui.o: ui.c fitsh.h
	$(CC) $(CFLAGS) -c ui.c

history.o: history.c fitsh.h
	$(CC) $(CFLAGS) -c history.c

maskdraw.o: maskdraw.c maskdraw.h
	$(CC) $(CFLAGS) -c maskdraw.c

kernel-base.o: kernel-base.c kernel.h fitsh.h
	$(CC) $(CFLAGS) -c kernel-base.c

kernel-io.o: kernel-io.c kernel.h fitsh.h
	$(CC) $(CFLAGS) -c kernel-io.c

kernel.a: kernel-base.o kernel-io.o
	$(AR) src kernel.a kernel-base.o kernel-io.o
	$(RANLIB) kernel.a

star-base.o: star-base.c stars.h
	$(CC) $(CFLAGS) -c star-base.c

star-draw.o: star-draw.c stars.h
	$(CC) $(CFLAGS) -c star-draw.c

star-model.o: star-model.c star-model.h stars.h
	$(CC) $(CFLAGS) -c star-model.c

star-psf.o: star-psf.c stars.h
	$(CC) $(CFLAGS) -c star-psf.c

star-cand-pp.o: star-cand-pp.c stars.h
	$(CC) $(CFLAGS) -c star-cand-pp.c

star-cand-biq.o: star-cand-biq.c stars.h
	$(CC) $(CFLAGS) -c star-cand-biq.c

star-cand-trb.o: star-cand-trb.c stars.h
	$(CC) $(CFLAGS) -c star-cand-trb.c

star-cand-lnk.o: star-cand-lnk.c stars.h
	$(CC) $(CFLAGS) -c star-cand-lnk.c

MOD_STARS=star-base.o star-model.o star-psf.o star-draw.o \
          star-cand-pp.o star-cand-biq.o star-cand-trb.o star-cand-lnk.o 

MOD_PSF=psf-io.o psf-base.o psf-determine.o

stars.a: $(MOD_STARS)
	$(AR) src stars.a $(MOD_STARS)
	$(RANLIB) stars.a

psf.a: $(MOD_PSF)
	$(AR) src psf.a $(MOD_PSF)
	$(RANLIB) psf.a

projection.o: projection.c projection.h
	$(CC) $(CFLAGS) -c projection.c

weight-gen.o: weight-gen.c weight.h
	$(CC) $(CFLAGS) -c weight-gen.c

weight-io.o: weight-io.c weight.h
	$(CC) $(CFLAGS) -c weight-io.c

weight-star.o: weight-star.c weight.h
	$(CC) $(CFLAGS) -c weight-star.c

longhelp.o: longhelp.c longhelp.h
	$(CC) $(CFLAGS) -c longhelp.c

optcalc.o: optcalc.c optcalc.h
	$(CC) $(CFLAGS) -c optcalc.c

###############################################################################

fiarith: fiarith.c fitsh.h $(DEP_FIARITH) $(MOD_FIARITH)
	$(CC) $(CFLAGS) -o fiarith fiarith.c $(DEP_FIARITH) $(MOD_FIARITH) $(LIB_FIARITH)

ficalib: ficalib.c fitsh.h $(DEP_FICALIB) $(MOD_FICALIB)
	$(CC) $(CFLAGS) -o ficalib ficalib.c $(DEP_FICALIB) $(MOD_FICALIB) $(LIB_FICALIB)

ficombine: ficombine.c fitsh.h $(DEP_FICOMBINE) $(MOD_FICOMBINE)
	$(CC) $(CFLAGS) -o ficombine ficombine.c $(DEP_FICOMBINE) $(MOD_FICOMBINE) $(LIB_FICOMBINE)

ficonv: ficonv.c fitsh.h $(DEP_FICONV) $(MOD_FICONV)
	$(CC) $(CFLAGS) -o ficonv ficonv.c $(DEP_FICONV) $(MOD_FICONV) $(LIB_FICONV)

fiheader: fiheader.c fitsh.h $(DEP_FIHEADER) $(MOD_FIHEADER)
	$(CC) $(CFLAGS) -o fiheader fiheader.c $(DEP_FIHEADER) $(MOD_FIHEADER) $(LIB_FIHEADER)

fiign: fiign.c fitsh.h $(DEP_FIIGN) $(MOD_FIIGN)
	$(CC) $(CFLAGS) -o fiign fiign.c $(DEP_FIIGN) $(MOD_FIIGN) $(LIB_FIIGN)

fiinfo-image.o: fiinfo-image.c fiinfo.h 
	$(CC) $(CFLAGS) -c fiinfo-image.c 

fiinfo-pnm.o: fiinfo-pnm.c fiinfo.h 
	$(CC) $(CFLAGS) -c fiinfo-pnm.c 

fiinfo: fiinfo.c fiinfo.h fitsh.h $(DEP_FIINFO) $(MOD_FIINFO)
	$(CC) $(CFLAGS) -o fiinfo fiinfo.c $(DEP_FIINFO) $(MOD_FIINFO) $(LIB_FIINFO)

firandom-eval.o: firandom-eval.c firandom.h
	$(CC) $(CFLAGS) -c firandom-eval.c

firandom: firandom.c fitsh.h $(DEP_FIRANDOM) $(MOD_FIRANDOM)
	$(CC) $(CFLAGS) -o firandom firandom.c $(DEP_FIRANDOM) $(MOD_FIRANDOM) $(LIB_FIRANDOM)

fiphot-io.o: fiphot-io.c fiphot.h
	$(CC) $(CFLAGS) -c fiphot-io.c

fiphot: fiphot.c fitsh.h $(DEP_FIPHOT) $(MOD_FIPHOT)
	$(CC) $(CFLAGS) -o fiphot fiphot.c $(DEP_FIPHOT) $(MOD_FIPHOT) $(LIB_FIPHOT)

fistar-io.o: fistar-io.c fistar.h
	$(CC) $(CFLAGS) -c fistar-io.c

fistar: fistar.c fitsh.h $(DEP_FISTAR) $(MOD_FISTAR) fistar.h
	$(CC) $(CFLAGS) -o fistar fistar.c $(DEP_FISTAR) $(MOD_FISTAR) $(LIB_FISTAR)

fitrans: fitrans.c fitsh.h $(DEP_FITRANS) $(MOD_FITRANS)
	$(CC) $(CFLAGS) -o fitrans fitrans.c $(DEP_FITRANS) $(MOD_FITRANS) $(LIB_FITRANS)

grcollect: grcollect.c fitsh.h $(DEP_GRCOLLECT) $(MOD_GRCOLLECT)
	$(CC) $(CFLAGS) -o grcollect grcollect.c $(DEP_GRCOLLECT) $(MOD_GRCOLLECT) $(LIB_GRCOLLECT)

grmatch: grmatch.c fitsh.h $(DEP_GRMATCH) $(MOD_GRMATCH)
	$(CC) $(CFLAGS) -o grmatch grmatch.c $(DEP_GRMATCH) $(MOD_GRMATCH) $(LIB_GRMATCH)

grtrans: grtrans.c fitsh.h $(DEP_GRTRANS) $(MOD_GRTRANS)
	$(CC) $(CFLAGS) -o grtrans grtrans.c $(DEP_GRTRANS) $(MOD_GRTRANS) $(LIB_GRTRANS)

gropt: gropt.c fitsh.h $(DEP_GROPT) $(MOD_GROPT)
	$(CC) $(CFLAGS) -o gropt gropt.c $(DEP_GROPT) $(MOD_GROPT) $(LIB_GROPT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

xfunct.o: xfunct.c xfunct.h
	$(CC) $(CFLAGS) -c xfunct.c

lfit-builtin.o: lfit-builtin.c lfit-builtin.h 
	$(CC) $(CFLAGS) -c lfit-builtin.c

lfit-info.o: lfit-info.c lfit-info.h 
	$(CC) $(CFLAGS) -c lfit-info.c

lfit: lfit.c $(DEP_LFIT) $(MOD_LFIT)
	$(CC) $(CFLAGS) -o lfit lfit.c $(DEP_LFIT) $(MOD_LFIT) $(DLDYN) $(LIB_LFIT)

linear.o: linear.c 
	$(CC) $(CFLAGS) -c linear.c 

linear.$(DLEXT): linear.o
	$(LD) -o linear.$(DLEXT) linear.o $(DLSWC) -lc
	chmod 644 linear.$(DLEXT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

fic_mpstack: fic_mpstack.c fitsh.h $(DEP_FIC_MPSTACK)
	$(CC) $(CFLAGS) -o fic_mpstack fic_mpstack.c $(DEP_FIC_MPSTACK) $(MOD_FIC_MPSTACK) $(LIB_FIC_MPSTACK)

###############################################################################

clean:
	$(MAKE) -C index clean
	$(MAKE) -C io clean
	$(MAKE) -C link clean
	$(MAKE) -C math clean
	$(MAKE) -C psn clean
	rm -f *.o *.a *.exe $(TARGETS)

