FC=gfortran-7

FOBJECTS1 = dophot.f90 tuneup.f90 paramfile.f90 fitssubs.f90 \
	io_dophot.f90 stat.f90 mathsub.f90 \
	makenoise.f90 warmstart.f90 addstar.f90 thresholds.f90 \
	skypar.f90 skyfun.f90 stpar.f90 pseud2d.f90 surface.f90 \
	chisq.f90 ludcmp.f90 lubksb.f90 \
	makemask.f90 fillerup.f90 transmask.f90 offpic.f90 volcano.f90 \
	isearch.f90 cosmic.f90 toobright.f90 oblit.f90 oblit2.f90 toofaint.f90 \
	shape.f90 chisqshape.f90 twofit.f90 \
	improve.f90 aperpar.f90 probgal.f90 \
	vpsf.f90 ellip.f90 shapemag.f90 notifyend.f90 \
	choose_io.f90 inte_io.f90 comp_io.f90 bina_io.f90 \
        daoph.f90 incomp_out.f90 shadow_out.f90 apcorsat.f90 \
	apcorr2.f90 apcfit.f90 apapply.f90 findmatch2.f90 finish.f90 \
	/usr/local/lib/libcfitsio.so
 
dophot: $(FOBJECTS1)
	$(FC) $(FOBJECTS1) -o dophot
