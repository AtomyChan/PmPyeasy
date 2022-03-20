CC    = gcc # the C compiler
COPTS =  -O2 -Wall -s # options for the C compiler
LIBS  =  -lm 

ALL   =  fwhm.o errmess.o pfitsin.o

fwhm: $(ALL)
	$(CC) $(COPTS) $(ALL) -o fwhm $(LIBS)

fwhm.o : $(STDH) fwhm.c
	$(CC) $(COPTS) -c fwhm.c

errmess.o : $(STDH) errmess.c
	$(CC) $(COPTS) -c errmess.c

pfitsin.o : $(STDH) pfitsin.c
	$(CC) $(COPTS) -c pfitsin.c
