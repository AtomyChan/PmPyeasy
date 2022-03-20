#! /usr/bin/python

#select candidate psf stars from a photometry file with iraf/daophot

import numpy as np
from pyraf import iraf
from iraf import digiphot,daophot
import os
from astropy.io import fits, ascii
from astropy.table import Table

def daophot_pstselect_iraf(image, photfile, pstfile, maxnpsf):
	'''
	INPUTS:
		image: image for which to build psf star list
		photfile: photometry file
		pstfile: output psf star list file
		maxnpsf: maximum number of psf stars
	'''

	daophot.pstselect.unlearn()
	daophot.pstselect.verify = 'no'	

	daophot.datapar.unlearn()
	daophot.daopars.unlearn()

	daophot.pstselect(image=image, photfile=photfile,pstfile=pstfile,maxnpsf=maxnpsf)	



if __name__ == "__main__":

	import optparse
	parser = optparse.OptionParser()

	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_photfile = 'default'
	parser.add_option('--photfile', dest='photfile', type='string', default=def_photfile, help='photometry file(s); default: %s.mag.?'%def_inputimage)

	def_pstfile = 'default'
	parser.add_option('--pstfile', dest='pstfile', type='string', default=def_pstfile, help='output psf star list(s); default: %s.pst.?'%def_inputimage)

	def_maxnpsf = 50
	parser.add_option('-n','--maxnpsf', dest='maxnpsf', type=int, default=def_maxnpsf, help='maximum number of psf stars; default: %s'%def_maxnpsf)

	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	pstfile = options.pstfile
	maxnpsf = options.maxnpsf

	daophot_pstselect_iraf(input_image, photfile, pstfile, maxnpsf)
