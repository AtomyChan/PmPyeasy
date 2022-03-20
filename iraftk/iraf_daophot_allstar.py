#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import allstar

def daophot_allstar_iraf(image, photfile, psfimage, allstarfile, rejfile, subimage, fitrad=5.0, psfrad=11, fitsky='yes', sannulus=8, wsannulus=10):
	'''
	INPUTS:
		image: 
		photfile: input photometry file; default image.psf.?
		psfimage: PSF image; default image.psf.?
		allstarfile: output photometry file; default image.als.?
		rejfile: output rejection file; default image.arj.?
		subimage: substracted image; default image.sub.?
		fitrad:
		psfrad:
		fitsky: recompute group sky value during the fit?
		sannulus: inner radius of the sky fit annulus
		wsannulus: width of the sky fitting annulus
	'''

	allstar.datapars.unlearn()
	allstar.daopars.unlearn()
	allstar.daopars.fitrad = fitrad
	allstar.daopars.psfrad = psfrad
	allstar.daopars.fitsky = fitsky
	allstar.daopars.sannulus = sannulus
	allstar.daopars.wsannulus = wsannulus
	allstar.unlearn()

	allstar(image=image, photfile=photfile,  psfimage=psfimage, allstarfile=allstarfile, rejfile=rejfile, subimage=subimage, verify='no', update='no')	



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_photfile = 'default'
	parser.add_option('--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='input photometry file')
	
	def_psfimage = 'default'
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the psf model image')
	
	def_allstarfile = 'default'
	parser.add_option('--allstarfile', dest = 'allstarfile', type= 'string', default = def_allstarfile, help='output photometry file')

	def_rejfile = 'default'
	parser.add_option('--rejfile', dest = 'rejfile', type= 'string', default = def_rejfile, help='output rejection file')

	def_fitrad = 5.0
	parser.add_option('-r','--fitrad', dest='fitrad', type=float, default=def_fitrad, help='fitting radius in scale units; default: %s'%def_fitrad)


	def_subimage = 'default'
	parser.add_option('--subimage', dest = 'subimage', type= 'string', default = def_subimage, help='the subtracted image')
	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	psfimage = options.psfimage
	allstarfile = options.allstarfile
	rejfile = options.rejfile
	fitrad = options.fitrad
	subimage =  options.subimage

	daophot_allstar_iraf(input_image, photfile, psfimage, allstarfile, rejfile, subimage, fitrad=fitrad)
