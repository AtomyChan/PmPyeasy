#! /usr/bin/python

import os
import numpy as np
from pyraf import iraf
from iraf import imcombine

def imcombine_iraf(input_images, output, combine='median', reject='none'):
	'''
	INPUTS: 
		input_images:
		output:
		combine: 'median', 'average', 'lmedian', 'sum', 'quaduature', 'nmodel'
		reject: 'none', 'minmax', 'ccdclip', 'crreject', 'sigclip', 'avsigclip', 'pclip'
	'''

	imcombine.unlearn()
		
	combine_valid = ['median', 'average', 'lmedian', 'sum', 'quaduature', 'nmodel']
	reject_valid = ['none', 'minmax', 'ccdclip', 'crreject', 'sigclip', 'avsigclip', 'pclip']

	if combine not in combine_valid:
		raise IOError('invalid input for combine method')
	if reject not in reject_valid:
		raise IOError('invalid input for rejection method')

	imcombine(input_images, output=output,  combine=combine, reject=reject)	



if __name__ == "__main__":
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_output = 'imgcombine_temp.fits'
	parser.add_option('-o', '--output', dest='output', type='string', default=def_output, help='output file; default: %s'%def_inputimage)

	options, remainder = parser.parse_args()

	input_images = options.input_image
	outfile = options.output

	print input_images, outfile
	imcombine_iraf(input_images, outfile)
