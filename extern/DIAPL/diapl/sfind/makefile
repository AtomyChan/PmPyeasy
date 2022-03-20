CC   = gcc		# the C compiler
COPT = -O3 -Wall	# options for the C compiler
LIBS = -lm		# libraries

OBJS = aperphot.o centroid.o indexx.o bkg.o covar.o pfitsin.o errmess.o sfind.o

UTIL = pfitsin.o

.c.o :
	${CC} ${COPT} -c $<


sfind: ${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile
