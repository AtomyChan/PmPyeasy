#! /usr/bin/python

import os
from astropy.io import fits

def update_image_wcs_header_info(input_image, wcs_image, output_image = None, clobber=True):
	'''
	copy wcs info keywords from wcs_image to input_image
	'''

	string_header_keys_to_copy = [
	    "CTYPE1",
	    "CTYPE2",
	    "CUNIT1",
	    "CUNIT2"
	    ]
	
	float_header_keys_to_copy = [
	    "EQUINOX",
	    "LONPOLE",
	    "LATPOLE",
	    "CRVAL1",
	    "CRVAL2",
	    "CRPIX1",
	    "CRPIX2",
	    "CD1_1",
	    "CD1_2",
	    "CD2_1",
	    "CD2_2",
	    "IMAGEW",
	    "IMAGEH",
	    "A_ORDER",
	    "A_0_2",
	    "A_1_1",
	    "A_2_0",
	    "B_ORDER",
	    "B_0_2",
	    "B_1_1",
	    "B_2_0",
	    "AP_ORDER",
	    "AP_0_1",
	    "AP_0_2",
	    "AP_1_0",
	    "AP_1_1",
	    "AP_2_0",
	    "BP_ORDER",
	    "BP_0_1",
	    "BP_0_2",
	    "BP_1_0",
	    "BP_1_1",
	    "BP_2_0"
	    ]
	    
	
	wcs_hdu = fits.open(wcs_image)
	wcs_header = wcs_hdu[0].header.copy()
	wcs_hdu.close()
	
	
	input_hdu = fits.open(input_image)
	image_data = input_hdu[0].data
	updated_imagefile_header = input_hdu[0].header.copy()
	for hk in string_header_keys_to_copy:
	    try:
	        updated_imagefile_header.update(hk, wcs_header[hk])
	    except:
	        print hk, "string header update failed"
	for hk in float_header_keys_to_copy:
	    try:
	        updated_imagefile_header.update(hk, float(wcs_header[hk]))
	    except:
	        print hk, "float header update failed"
	
	
	output_hdu = fits.PrimaryHDU(image_data)
	output_hdu.header = updated_imagefile_header
	
	output_hdu.verify("fix")
	output_hdulist = fits.HDUList([output_hdu])
	if output_image is None:
		output_image =input_image.replace(".fits", ".wcs.fits") 

	if os.path.exists(output_image):
		if clobber:
			os.remove(output_image)	
		else:
			print "output image %s exists which will be overwrote"

	output_hdulist.writeto(output_image)
	input_hdu.close()


if __name__ == "__main__":
	import sys
		
	if len(sys.argv)<3:
		print "usage: %s input_image wcs_image [output_image]"%sys.argv[0]
		sys.exit()

	input_image = sys.argv[1]
	wcs_image = sys.argv[2]
	if len(sys.argv)>3:
		output_image = sys.argv[3]
	else:
		output_image = None

	update_image_wcs_header_info(input_image, wcs_image, output_image = output_image, clobber=True)
