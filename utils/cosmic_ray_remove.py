#! /usr/bin/python

# The following two lines are only needed as cosmic.py is not in this directory nor in the python path.
# They would not be required if you copy cosmics.py in this directory.
import sys
sys.path.append("/home/supernova/software/cosmic_ray_removal/cosmics.py_0.4") # The directory that contains cosmic.py


import cosmics

def remove_CR(inimage, outimage, gain, readnoise, sigclip=3, maxiter=2):

	# Read the FITS :
	array, header = cosmics.fromfits(inimage)
	# array is a 2D numpy array
	
	# Build the object :
	c = cosmics.cosmicsimage(array, gain=gain, readnoise=readnoise, sigclip = sigclip, sigfrac = 0.3, objlim = 5.0)
	# There are other options, check the manual...
	
	# Run the full artillery :
	c.run(maxiter = maxiter)
	
	# Write the cleaned image into a new FITS file, conserving the original header :
	cosmics.tofits(outimage, c.cleanarray, header)
	
	# If you want the mask, here it is :
	cosmics.tofits("mask_CR.fits", c.mask, header)
	# (c.mask is a boolean numpy array, that gets converted here to an integer array)


if __name__ == "__main__":
	import sys
	
	#OSIRIS 
	gain = 0.95
	readnoise = 4.5

	inimage = sys.argv[1]
	outimage = sys.argv[2]


	remove_CR(inimage, outimage, gain, readnoise, sigclip=3, maxiter=2)
