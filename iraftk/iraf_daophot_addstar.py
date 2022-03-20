#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import addstar

def daophot_addstar_iraf(image, addimg, psfimage, photfile, psfrad=None):
	'''
	INPUTS:
		image:
		groupfile:
		psfimage:
		nstarfile:
		rejfile:
	'''
	addstar.datapars.unlearn()
	addstar.daopars.unlearn()
	if psfrad is not None:
		addstar.daopars.psfrad = psfrad
	addstar.unlearn()

	addstar(image=image, photfile=photfile,  addimage=addimg, psfimage=psfimage, verify='no', update='no')



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_addimg = 'default'
	parser.add_option('-o', '--addimg', dest = 'addimage', type= 'string', default = def_addimg, help='output image')
	
	def_psfimage = ''
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the output psf model image')
	
	def_photfile = ''
	parser.add_option('--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='photometry file containing the source list to add')

	options, remainder = parser.parse_args()
	input_image = options.input_image
	addimg = options.addimage
	psfimage = options.psfimage
	photfile = options.photfile

	daophot_addstar_iraf(input_image, addimg, psfimage, photfile)
