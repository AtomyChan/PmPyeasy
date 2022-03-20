#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import substar

def daophot_substar_iraf(image, photfile, exfile, psfimage, subimage):
	'''
	INPUTS:
		image:
		photfile:
		exfile:
		psfimage:
		subimage:
	'''

	substar.datapars.unlearn()
	substar.daopars.unlearn()
	substar.unlearn()

	substar(image=image, photfile=photfile, exfile=exfile, psfimage=psfimage, subimage=subimage,  verify='no', update='no')	



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')
	def_photfile = 'default'
	parser.add_option('--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='output photometry file')

	def_exfile = 'default'
	parser.add_option('--exfile', dest = 'exfile', type= 'string', default = def_exfile, help='input exclude file')

	def_psfimage = 'default'
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the output psf model image')
	
	def_subimage = 'default'
	parser.add_option('--subimage', dest = 'subimage', type= 'string', default = def_subimage, help='input group file')

	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	exfile = options.exfile
	psfimage = options.psfimage
	subimage = options.subimage

	daophot_substar_iraf(input_image, photfile, exfile, psfimage, subimage)
