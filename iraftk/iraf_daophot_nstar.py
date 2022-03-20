#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import nstar

def daophot_nstar_iraf(image, groupfile, psfimage, nstarfile, rejfile):
	'''
	INPUTS:
		image:
		groupfile:
		psfimage:
		nstarfile:
		rejfile:
	'''
	nstar.datapars.unlearn()
	nstar.daopars.unlearn()
	nstar.unlearn()

	nstar(image=image, groupfile=groupfile,  psfimage=psfimage, nstarfile=nstarfile, rejfile=rejfile, verify='no', update='no')	



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_groupfile = 'default'
	parser.add_option('--groupfile', dest = 'groupfile', type= 'string', default = def_groupfile, help='input group file')
	
	def_psfimage = 'default'
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the output psf model image')
	
	def_nstarfile = 'default'
	parser.add_option('--nstarfile', dest = 'nstarfile', type= 'string', default = def_nstarfile, help='output photometry file')

	def_rejfile = 'default'
	parser.add_option('--rejfile', dest = 'rejfile', type= 'string', default = def_rejfile, help='output rejection file')

	options, remainder = parser.parse_args()
	input_image = options.input_image
	groupfile = options.groupfile
	psfimage = options.psfimage
	nstarfile = options.nstarfile
	rejfile = options.rejfile

	daophot_nstar_iraf(input_image, groupfile, psfimage, nstarfile, rejfile)
