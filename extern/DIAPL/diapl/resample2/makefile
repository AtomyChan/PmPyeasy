CC   = gcc		# the C compiler
COPT = -O3 -s -Wall 	# options for the C compiler
LIBS = -lm		# libraries

OBJS = poly.o spline.o splint.o bicspl.o hedit.o \
       indexx.o fix_cosmics.o resample2.o rm_bad_col.o res_bad_col.o \
       gasdev.o ran1.o errmess.o pfitsio1.o

DEFS = bicspl.c splint.c resample2.c

UTIL = rfits_u.o wfits_f.o errmess.o

.c.o :
	${CC} ${COPT} -c $<


resample2: ${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile

${DEFS}	: defs.h

${UTIL}	: errmess.h pfitsio1.h
