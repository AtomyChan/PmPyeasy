CC   = gcc       # the C compiler
COPT = -O3 -Wall -s # options for the C compiler
LIBS = -lm       # libraries

OBJS = crf.o fourn.o get_max.o get_params.o pack_data.o \
       get_sector.o pfitsin1.o errmess.o cross.o

FUNC = crf.o cross.o


DEFS = crf.o get_params.o get_sector.o cross.o pack_data.o

UTIL = pfitsin1.o errmess.o

.c.o :
	${CC} ${COPT} -c $<
     

cross: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}

${OBJS} : makefile

${FUNC} : funcs.h

${DEFS} : defs.h
