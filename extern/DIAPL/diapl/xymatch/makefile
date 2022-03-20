CC   = gcc		# the C compiler
COPT = -O3 -s -Wall 	# options for the C compiler
LIBS = -lm		# libraries

OBJS = tri_gen.o triangles.o bright_end.o xy_lin.o indexx.o refine.o \
       errmess.o xymatch.o

.c.o :
	${CC} ${COPT} -c $<


xymatch: ${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile
