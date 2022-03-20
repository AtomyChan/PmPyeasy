CC   = gcc		# the C compiler
COPT = -O2 -s -Wall 	# options for the C compiler
LIBS = -lm		# libraries

OBJS =	mstack.o histogram.o bin.o get_peak.o median.o \
        get_params.o read_sector.o write_sector.o pfitsio1.o errmess.o

FUNC =	mstack.o

DEFS =	get_params.o mstack.o

.c.o:
	${CC} ${COPT} -c $<

mstack: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}


${OBJS} : makefile

${FUNC}	: funcs.h

${DEFS}	: defs.h
