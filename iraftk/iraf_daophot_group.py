#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import group

def daophot_group_iraf(image, photfile, psfimage, groupfile):
	'''
	INPUTS:
		image:
		photfile:
		psfimage:
		groupfile:
	'''

	group.datapars.unlearn()
	group.daopars.unlearn()
	group.unlearn()

	group(image=image, photfile=photfile, psfimage=psfimage, groupfile=groupfile,  verify='no', update='no')	



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')
	def_photfile = 'default'
	parser.add_option('--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='output photometry file')

	def_psfimage = 'default'
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the output psf model image')
	
	def_groupfile = 'default'
	parser.add_option('--groupfile', dest = 'groupfile', type= 'string', default = def_groupfile, help='input group file')

	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	psfimage = options.psfimage
	groupfile = options.groupfile

	daophot_group_iraf(input_image, photfile, psfimage, groupfile)
