CC   = gcc		# the C compiler
COPT = -O3 -s -Wall 	# options for the C compiler
LIBS = -lm		# libraries

OBJS =	cutfitsim.o errmess.o swap.o pfitshead.o

.c.o:
	${CC} ${COPT} -c $<

cutfitsim: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}

${OBJS} : makefile
