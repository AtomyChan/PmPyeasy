SHELL=/bin/sh

FIRANDOM=../../../src/firandom

TARGETS=uniform.fits \
	globularcluster.fits \
	coma1.fits \
	coma2.fits \
	background.fits \
	vignetting.fits \
	grid1.fits \
	grid2.fits \
	grid3.fits 

all: $(TARGETS)

.PHONY: all clean

# Simple image: 2000 uniformly (both in mag and spatially) distributed 
# stars, all of them have a fwhm of 2.3. Background has gaussian 
# distribution with a mean of 500 (thus, sky) and stddev of 10 (skysigma).
uniform.fits:
	$(FIRANDOM) --size 512,512 --sky "g(500,10)" \
	--list "f=4.3, 2000*[x=r(-1,1),y=r(-1,1),m=r(11,15)]" \
	--mag-flux 10,100000 --bitpix -32 --output uniform.fits

# A globular cluster: 2000 stars with a gaussian spatial distribution
# (being located at the center and having an fwhm of 0.1 of the image size),
# and some field stars (200 brighter and 2000 fainter...).
globularcluster.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)" \
	-l "f=2.3, 2000*[x=g(0,0.2),y=g(0,0.2),m=r(11,15)],200*[x=r(-1,1),y=r(-1,1),m=r(10,16)],2000*[x=r(-1,1),y=r(-1,1),m=r(13,18)]" \
	--mag-flux 10,10000 -b -32 -o globularcluster.fits

# The coma-effect. The star profile is described by the parameters
# s, d and k. s is the sigma value of a gaussian profile (the f, fwhm
# is always 2.34 times larger than s). The d and k vales are the deviation
# momenta, so the profile can be described by the equation
# exp ( - 1/2 x M^-1 x ), where M is the matrix
#	[ s+d	k ] 
#	[ 	  ]
#	[ k   s-d ].
# The strength of the effect can be described with this 
COMA=0.2
# parameter. $(COMA) should be between 0 and 1.
# The ellipticity of the stars can be described by two set of parameters:
# this s, d and k _or_ by f (fwhm), e (ellipticity) and p (position angle)
# (see other examples below).

coma1.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)" \
	-l "1000*[x=r(-1,1),y=r(-1,1),m=r(10,15),s=1.5,d=+$(COMA)*s*(x^2-y^2),k=+2*$(COMA)*s*x*y]"\
	--mag-flux 10,10000 -b -32 -o coma1.fits

# Another type of coma effect: the sign is changed in d and k.
coma2.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)" \
	-l "1000*[x=r(-1,1),y=r(-1,1),m=r(10,15),s=1.5,d=-$(COMA)*s*(x^2-y^2),k=-2*$(COMA)*s*x*y]"\
	--mag-flux 10,10000 -b -32 -o coma2.fits

# The background is not constant...
background.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)+x*20-y*8" \
	-l "4000*[x=r(-1,1),y=r(-1,1),m=r(10,17)]"\
	--mag-flux 10,10000 -b -32 -o background.fits

# the simulation of the vignetting effect:
vignetting.fits:
	$(FIRANDOM) -s 512,512 -m "200+600*(1-0.4*(x^2+y^2))+g(0,10)" \
	-l "f=2.7,5000*[x=r(-1,1),y=r(-1,1),m=r(10,15)-2.5*lg(1-0.4*(x^2+y^2))]" \
	--mag-flux 10,10000 -b -32 -o vignetting.fits
# An equivalent way, the flux of the stars can be defined directly, 
# using the parameter 'i'. If the flux of the stars are specified by only 
# this parameter, the --max-flux switch can be omitted (however, it
# can always be omitted, the default value is --max-flux 10,10000 ;)).
# -l "f=2.7,5000*[x=r(-1,1),y=r(-1,1),i=(1-0.4*(x^2+y^2))*10000*10^(-0.4*r(0,5))]" \

# A perfect grid: in a block of N*[...], the global variable 'n' is changed
# between 0 and N-1 so it can be used for generating grids.
grid1.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)" \
	-l "f=2.5,100*[x=-0.5+0.1*div(n,10),y=-0.5+0.1*mod(n,10),m=r(10,12)]" \
	--mag-flux 10,10000 -b -32 -o grid1.fits

# Another grid which is not so perfect...
grid2.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)" \
	-l "f=2.5,100*[x=-0.7+0.14*div(n,10)+r(0,0.04),y=-0.7+0.14*mod(n,10)+r(0,0.04),m=r(10,12)]" \
	--mag-flux 10,10000 -b -32 -o grid2.fits

# Another grid where other parameters (fwhm, ellipticity) are also varying 
# definitely, and the flux is constant:
grid3.fits:
	$(FIRANDOM) -s 512,512 -m "g(500,10)" \
	-l "100*[x=-0.7+0.14*div(n,10),y=-0.7+0.14*mod(n,10),m=10,f=1.0+0.2*div(n,10),e=0.04*mod(n,10),p=30]" \
	--mag-flux 10,10000 -b -32 -o grid3.fits

clean:
	rm -f $(TARGETS)

