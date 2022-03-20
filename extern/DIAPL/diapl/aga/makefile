CC   = gcc	# the C compiler
COPT = -O3 -Wall -s	# options for the C compiler
LIBS = -lm	# libraries

OBJS = base_func.o clone_domains.o init_difimages.o \
       ludcmp.o lubksb.o make_domains.o make_vectors.o get_headlen.o \
       get_params.o read_sector.o write_sector.o xy_convolve.o \
       mask_badpix.o clean_domains.o apply_norm.o expand_matrix.o pfitsio1.o \
       spatial_convolve.o local_clip.o local_solution.o clip_domains.o \
       quick_sort.o max.o get_domains.o write_kernel.o errmess.o aga.o

DEFS = clip_domains.o local_clip.o init_difimages.o make_domains.o \
       make_vectors.o mask_badpix.o clean_domains.o clone_domains.o \
       expand_matrix.o get_domains.o get_params.o local_solution.o \
       spatial_convolve.o write_kernel.o aga.o

FUNC = clip_domains.o get_domains.o init_difimages.o get_params.o \
       local_clip.o make_vectors.o aga.o

UTIL = write_sector.o read_sector.o init_difimages.o

.c.o :
	${CC} ${COPT} -c $<

aga: ${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile

${DEFS} : defs.h

${FUNC} : funcs.h
