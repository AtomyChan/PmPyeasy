CC   = gcc           # the C compiler
COPT = -O3 -Wall -s  # options for the C compiler
LIBS = -lm           # libraries

OBJS = hedit.o pfitsio1.o errmess.o template.o

UTIL = pfitsio1.o errmess.o

.c.o :
	${CC} ${COPT} -c $<

template: ${OBJS}
	${CC} -o $@ ${OBJS} ${LIBS} ${COPT}

${OBJS} : makefile 

${UTIL} : pfitsio1.h errmess.h
