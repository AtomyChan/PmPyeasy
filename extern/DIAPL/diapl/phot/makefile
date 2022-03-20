CC   = gcc		# the C compiler
COPT = -O3 -Wall -s	# options for the C compiler
LIBS = -lm      	# libraries

OBJS = indexx.o bkg.o im_convolve.o phot.o \
       read_psf.o init_psf.o psf_core.o read_kernel.o spatial_coeffs.o \
       make_vectors.o base_func.o make_kernel.o get_params.o make_psf.o \
       get_phot.o get_fwhm.o errmess.o pfitsin1.o

DEFS = bkg.o get_params.o get_phot.o init_psf.o phot.o make_kernel.o \
       make_psf.o make_vectors.o psf_core.o read_kernel.o read_psf.o \
       spatial_coeffs.o get_fwhm.o

FUNC = bkg.o phot.o make_psf.o

UTIL = phot.o

.c.o:
	${CC} ${COPT} -c $<

phot: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}


${OBJS}: makefile

${DEFS} : defs.h

${FUNC} : funcs.h
