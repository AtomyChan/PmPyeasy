#! /usr/bin/python

from pyraf import iraf
from iraf import daophot
from iraf import peak
from iraf import substar


def daophot_peak_iraf(image, photfile, psfimage, peakfile, rejfile, fitrad=3, refitsky='yes', recenter='yes', sannulus=None, wsannulus=10, psfrad=11):
	'''
	INPUTS:
		image:
		photfile:
		psfimage:
		peakfile:
		rejfile:
	'''
	peak.unlearn()
	peak.datapars.unlearn()
	peak.daopars.unlearn()
	peak.daopars.fitsky = refitsky
	peak.daopars.recenter = recenter
	if sannulus is None:
		sannulus = fitrad
	peak.daopars.sannulus = sannulus
	peak.daopars.wsannulus = wsannulus
	peak.daopars.fitrad = fitrad
	peak.daopars.psfrad = psfrad

	peak(image=image, photfile=photfile,  psfimage=psfimage, peakfile=peakfile, rejfile=rejfile, verify='no', update='no')	


def daophot_substar_iraf(image, photfile, exfile, psfimage, subimage, fitrad=3, psfrad=11):
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
	substar.daopars.fitrad = fitrad
	substar.daopars.psfrad = psfrad
	substar.unlearn()

	substar(image=image, photfile=photfile, exfile=exfile, psfimage=psfimage, subimage=subimage,  verify='no', update='no')	


if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()

	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_photfile = 'default'
	parser.add_option('--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='input photometry file')
	
	def_psfimage = 'default'
	parser.add_option('--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the input psf model image')
	
	def_peakfile = 'default'
	parser.add_option('--peakfile', dest = 'peakfile', type= 'string', default = def_peakfile, help='output photometry file')

	def_rejfile = 'default'
	parser.add_option('--rejfile', dest = 'rejfile', type= 'string', default = def_rejfile, help='output rejection file')

	def_fitrad = 5.0
	parser.add_option('-r','--fitrad', dest='fitrad', type=float, default=def_fitrad, help='fitting radius in scale units; default: %s'%def_fitrad)

	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	psfimage = options.psfimage
	peakfile = options.peakfile
	rejfile = options.rejfile
	fitrad = options.fitrad

	daophot_peak_iraf(input_image, photfile, psfimage, peakfile, rejfile, fitrad=fitrad)
