#! /usr/bin/python

from pyraf import iraf
from iraf import imred
from iraf import bias
from iraf import colbias

import os

def colbias_iraf(input_img, output_img, biassec='[3:52, 1:2102]',  trimsec = '[220:1960, 185:1925]'):

	if os.path.exists(output_img):
		os.remove(output_img)

	colbias.unlearn()
	print input_img
	print output_img
	colbias(input_img, output_img, bias= biassec, trim=trimsec, interactive=False,)




if __name__ == "__main__":
	import sys
	
	input_img = sys.argv[1]
	output_img = sys.argv[2]

	colbias_iraf(input_img, output_img)
