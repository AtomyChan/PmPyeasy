#! /anaconda/bin/python

import ccdproc
from astropy.io import fits

from astropy import units as u

def remove_CR_ccdproc_cosmicray_lacosmic(image_in, image_out, gain, readnoise, sigclip=None):
	'''
	see ccdproc.cosmicray_lacosmic for detailed parameters
	
	here:
	http://ccdproc.readthedocs.io/en/latest/api/ccdproc.cosmicray_lacosmic.html#ccdproc.cosmicray_lacosmic

	'''

	hdu = fits.open(image_in)
	data = hdu[0].data
	hdr = hdu[0].header
	ccddata = ccdproc.CCDData(data, unit=u.adu)
	data_with_deviation = ccdproc.create_deviation(ccddata, gain= gain* u.electron / u.adu, readnoise= readnoise * u.electron)

	if sigclip is None:
		sigclip = 4.0

	#clean_data = ccdproc.cosmicray_lacosmic(data_with_deviation, thresh=sigclip ) #old version and didn't work...
	clean_data = ccdproc.cosmicray_lacosmic(data_with_deviation, sigclip=sigclip )
	fits.writeto(image_out, data=clean_data.data, header=hdr)


if __name__ == "__main__":
	import sys
	
	#default
	def_gain = 1
	def_readnoise = 10
	def_sigclip = 4.0

	image_in = sys.argv[1]
	if image_in == '--help':
		print "ccdproc_cosmicray_removal_lacosmic.py image_in image_out [gain] [readnoise] [sigclip]"
		sys.exit()

	image_out = sys.argv[2]
	
	if len(sys.argv)>5:
		gain = float(sys.argv[3])
		readnoise = float(sys.argv[4])
		sigclip = float(sys.argv[5])
	else:
		gain = def_gain
		readnoise = def_readnoise
		sigclip = def_sigclip

	remove_CR_ccdproc_cosmicray_lacosmic(image_in, image_out, gain, readnoise, sigclip=sigclip)

