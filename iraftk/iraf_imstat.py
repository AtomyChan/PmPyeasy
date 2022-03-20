#! /usr/bin/python

from collections import OrderedDict
from pyraf import iraf
from iraf import imstatistics

def imstat_iraf(images,fields = "image,npix,mode, midpt,mean,stddev,min,max",lower = 'INDEF',upper = 'INDEF',nclip = 0,lsigma = 3.0,usigma = 3.0,binwidth = 1, verbose=1):
	'''
	compute and print image pixel statistics

	INPUTS:
	    images: The  input  images  or image sections for which pixel statistics are to be computed.
    	    fields: The statistical quantities to be computed and printed. The available
  		    fields are the following.
         	    image - the image name
          	     npix - the number of pixels used to do the statistics
         	     mean - the mean of the pixel distribution
       		    midpt - estimate of the median of the pixel distribution
       		     mode - the mode of the pixel distribution
             	   stddev - the standard deviation of the pixel distribution
              	     skew - the skew of the pixel distribution
         	 kurtosis - the kurtosis of the pixel distribution
            	      min - the minimum pixel value
            	      max - the maximum pixel value
    
	    lower: The minimum good data limit.  All pixels are above  the  default value of INDEF.
    	    upper: The  maximum  good data limit.  All pixels are above the default value of INDEF.
    	    nclip: The maximum number of iterative clipping cycles. By default no clipping is performed.
    	    lsigma: The low side clipping factor in sigma.
    	    usigma: The high side clipping factor in sigma.
    	    binwidth: The  width of the histogram bins used for computing the midpoint (estimate of the median) 		     and the mode.  The units are in sigma.

	'''

	imstatistics.unlearn()

	out = imstatistics(images, fields = fields, lower = lower,upper = upper,nclip = nclip,lsigma = lsigma,usigma = usigma,binwidth = binwidth, Stdout=1)
	
	if verbose:
		print out

	stat_ret  = out[1].split()

	stat_retdict = OrderedDict()
	for i,field in enumerate(fields.split(',')):
		stat_retdict[field] = stat_ret[i]


	return stat_retdict



if __name__ == "__main__":
	import sys

	images = sys.argv[1]

	if images == '--help':
		print "iraf_imstat.py input_images [nclip] [output_fields]"
		sys.exit()

	if len(sys.argv) > 2:
		nclip = int(sys.argv[2])
	else:
		print "nclip=0 will be used, so no clipping will be performed"
		nclip = 0

	if len(sys.argv) > 3:
		fields = sys.argv[3]
	else:
		fields = "image,npix,mode, midpt,mean,stddev,min,max"
	


	out = imstat_iraf(images,fields, nclip=nclip)
	for field in fields.split(','):
		print "%s: %s"%(field, out[field])
