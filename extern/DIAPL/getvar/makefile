CC   = gcc		# the C compiler
COPT = -O3 -Wall -s	# options for the C compiler
LIBS = -lm      	# libraries

OBJS = getvar.o im_convolve.o neighbor.o indexx.o errmess.o pfitsio1.o \
       read_psf.o init_psf.o psf_core.o read_kernel.o spatial_coeffs.o \
       make_vectors.o base_func.o make_kernel.o get_params.o make_psf.o \
       get_fwhm.o covar_sig.o centroid.o aperphot.o bkg.o center_stars.o \
       quick_sort.o get_repeaters.o find_stars.o cross_id.o get_phot.o \

DEFS = get_params.o bkg.o center_stars.o cross_id.o \
       find_stars.o get_fwhm.o get_params.o get_phot.o get_repeaters.o \
       init_psf.o getvar.o make_kernel.o make_psf.o make_vectors.o neighbor.o \
       psf_core.o read_kernel.o read_psf.o spatial_coeffs.o

FUNC = bkg.o center_stars.o get_repeaters.o getvar.o make_psf.o make_vectors.o

UTIL = getvar.o

.c.o:
	${CC} ${COPT} -c $<

getvar:	${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}

${OBJS} : makefile

${DEFS} : defs.h

${FUNC} : funcs.h
