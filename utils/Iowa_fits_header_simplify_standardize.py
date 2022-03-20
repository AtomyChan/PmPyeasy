#! /usr/bin/python

from astropy.io import fits
import os

def simplify_and_standardize_fits_header(img, output, out_header_only = False, overwrite=False):
	'''
	This script is used to solve the non-standard wcs info in Iowa images
	'''

	hdu = fits.open(img)
	hdr_old = hdu[0].header

	hdu_new = fits.PrimaryHDU()

	hdu_new.header['SIMPLE'] = hdr_old['SIMPLE']
	hdu_new.header['BITPIX'] = hdr_old['BITPIX']
	hdu_new.header['NAXIS'] = hdr_old['NAXIS']
	hdu_new.header['NAXIS1'] = hdr_old['NAXIS1']
	hdu_new.header['NAXIS2'] = hdr_old['NAXIS2']

	hdu_new.header['CTYPE1'] = hdr_old['CTYPE1']
	hdu_new.header['CTYPE2'] = hdr_old['CTYPE2']
	hdu_new.header['CRPIX1'] = hdr_old['CRPIX1']
	hdu_new.header['CRPIX2'] = hdr_old['CRPIX2']
	hdu_new.header['CRVAL1'] = hdr_old['CRVAL1']
	hdu_new.header['CRVAL2'] = hdr_old['CRVAL2']
	hdu_new.header['CROTA1'] = hdr_old['CROTA1']
	hdu_new.header['CROTA2'] = hdr_old['CROTA2']
	hdu_new.header['CDELT1'] = hdr_old['CDELT1']
	hdu_new.header['CDELT2'] = hdr_old['CDELT2']

	if not out_header_only:
		hdu_new.data = hdu[0].data

	if os.path.exists(output):
		if overwrite:
			os.remove(output)
		else:
			raise IOError('%s already exists...')		

	hdu_out = fits.HDUList([hdu_new])
	hdu_out.writeto(output)


if __name__ == "__main__":
	import sys
	
	in_img = sys.argv[1]
	out_img = sys.argv[2]

	simplify_and_standardize_fits_header(in_img, out_img, out_header_only = False)
