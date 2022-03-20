CC   = gcc		# the C compiler
COPT = -O3 -s -Wall	# options for the C compiler
LIBS = -lm		# libraries

OBJS = xygrid.o fitn.o lubksb.o ludcmp.o poly.o sigma.o errmess.o

.c.o :
	${CC} ${COPT} -c $<


xygrid: ${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile
