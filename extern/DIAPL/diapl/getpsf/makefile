ORIG  = .		# where the sources lie
DEST  = .		# where to put the executable
CC    = gcc		# the C compiler
COPT  =  -O3 -s -Wall	# options for the C compiler
LIBS  =  -lm		# libraries


OBJS = getpsf.o bkg.o is_peak.o local_max.o psf_stars.o \
       silly_phot.o quick_sort.o get_params.o center.o is_isolated.o \
       refine_max.o copy_obj.o psf_bitmaps.o fitpsf.o ludcmp.o lubksb.o \
       get_gaussian.o init_psf.o psf_core.o local_fit.o clip_stars.o \
       write_psf.o errmess.o pfitsin1.o

DEFS = bkg.o center.o clip_stars.o copy_obj.o fitpsf.o get_gaussian.o \
       get_params.o init_psf.o is_isolated.o is_peak.o local_fit.o getpsf.o \
       psf_bitmaps.o psf_core.o psf_stars.o refine_max.o silly_phot.o \
       write_psf.o

FUNC = bkg.o clip_stars.o fitpsf.o get_gaussian.o local_fit.o local_max.o \
       getpsf.o psf_stars.o

.c.o :
	${CC} ${COPT} -c $<

getpsf:	${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile

${DEFS} : defs.h

${FUNC} : funcs.h

