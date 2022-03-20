CC   = gcc		# the C compiler
COPT = -O3 -s -Wall 	# options for the C compiler
LIBS = -lm		# libraries

OBJS =	im2float.o errmess.o pfitsio.o

.c.o:
	${CC} ${COPT} -c $<

im2float: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}

${OBJS} : makefile
