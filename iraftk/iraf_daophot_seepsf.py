#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import seepsf

def daophot_seepsf_iraf(psfimage, outputimg, xpsf='INDEF', ypsf='INDEF'):
	'''
	INPUTS:
		psfimage:
		outputimg:
	'''
	
	seepsf(psfimage, outputimg, xpsf=xpsf, ypsf=ypsf)



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_psfimage = ''
	parser.add_option('-i','--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='input psf image')

	def_outputimg = ''
	parser.add_option('-o','--outputimg', dest = 'outputimg', type= 'string', default = def_outputimg, help='output psf image')
	
	
	options, remainder = parser.parse_args()
	psfimage = options.psfimage
	outputimg = options.outputimg

	daophot_seepsf_iraf(psfimage, outputimg)
