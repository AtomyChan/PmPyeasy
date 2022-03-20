CC    = gcc # the C compiler
COPTS =  -Wall -s -O2 # options for the C compiler

ALL   =  fitshedit.o errmess.o pfitshead.o

fitshedit: $(ALL)
	$(CC) $(COPTS) $(ALL) -o fitshedit

fitshedit.o: $(STDH) fitshedit.c
	$(CC) $(COPTS) -c fitshedit.c

errmess.o : $(STDH) errmess.c
	$(CC) $(COPTS) -c errmess.c

pfitshead.o : $(STDH) pfitshead.c
	$(CC) $(COPTS) -c pfitshead.c
