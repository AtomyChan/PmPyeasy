from pyraf import iraf
from iraf import daophot
from iraf import psf
import os

def daophot_psf_iraf(image, photfile, pstfile, psfimage, opstfile=None, groupfile=None, psfrad=15.0, fitrad=9, varorder=0):
	'''
	INPUTS:
		image:
		photfile:
		pstfile:
		psfimage:
		opstfile:
		groupfile:
		psfrad: the radius of the circle in scale units within which the PSF model is dedined
	'''

	psf.datapars.unlearn()
	psf.daopars.unlearn()
	psf.daopars.varorder=varorder
	psf.daopars.psfrad = psfrad
	psf.daopars.fitrad = fitrad
	psf.unlearn()

	if opstfile is None:
		opstfile = 'default'

	if groupfile is None:	
		groupfile = 'default'

	if os.path.exists(psfimage):
		os.remove(psfimage)

	psf(image=image, photfile=photfile, pstfile=pstfile, psfimage=psfimage, opstfile=opstfile, groupfile=groupfile, matchbyid='yes', interactive='no', mkstars='no', showplots='no', verify='no', update='no')	



if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_photfile = 'default'
	parser.add_option('-p','--photfile', dest = 'photfile', type= 'string', default = def_photfile, help='input photometry file')
	
	def_pstfile = 'default'
	parser.add_option('-s','--pstfile', dest = 'pstfile', type= 'string', default = def_pstfile, help='input psf star photometry file')
	
	def_psfimage = 'default'
	parser.add_option('-f','--psfimage', dest = 'psfimage', type= 'string', default = def_psfimage, help='the output psf model image')
	
	def_psfrad = 15.0
	parser.add_option('-r','--psfrad', dest='psfrad', type=float, default=def_psfrad, help='the radius of the circle within which the psf model is defined; default: %s'%def_psfrad)

	def_varorder = 0
	parser.add_option('-v','--varorder', dest='varorder', type=int, default=def_varorder, help='the order of variability of the psf model computed by the DAOPHOT psf task; default: %s'%def_varorder)

	options, remainder = parser.parse_args()
	input_image = options.input_image
	photfile = options.photfile
	pstfile = options.pstfile
	psfimage = options.psfimage
	psfrad   = options.psfrad	
	varorder = options.varorder

	daophot_psf_iraf(input_image, photfile, pstfile, psfimage, opstfile=None, groupfile=None, psfrad =psfrad, varorder=varorder)
