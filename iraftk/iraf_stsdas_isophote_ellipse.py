#! /anaconda/bin/python

import numpy as np
import pyraf
from pyraf import iraf
from iraf import stsdas
from iraf import analysis
from iraf import isophote
from iraf import ellipse
from iraf import bmodel
try:
	from iraf import cmodel
except:
	print 'cmodel not available'

def stsdas_analysis_isophote_ellipse(input_image, output, x0, y0, minsma, maxsma, ellip0=0.5, pa0=65, sma0=10,step=0.05, linear=False,  recenter=False, xylearn=False, physical=True, conver=0.05, minit=20, maxit=80, hcenter=False, hellip=False,hpa=False, wander='INDEF', maxgerr=1.0,  olthresh=0.0, integrmode='bi-linear', usclip=3.0, lsclip=3.0, nclip=0, fflag=0.5, mag0=0.0, refer=1.0, zerolevel=0.0, interactive=False):
	'''
	Fit elliptical isophote to galaxy image
	    
	INPUTS:
		x0, y0, ellip0, pa0, sma0: initial ellipse parameters
		x0,y0:center coordinate X,Y; [>=1.0]
		ellip0: ellipticity (defined as 1-b/a), [0.05 to 1.0]
		pa0: position angle [-90 to 90]
		sma0: semi-major axis length, [>=5.0]

		minsma: minimum semi-major axis length to be measured, if set to zero, if central pixel intensity will be measured
		maxsma: maximum semi-major axis length to be measured
	'''
	
	ellipse.unlearn()
	#ellipse.x0 = x0
	#ellipse.y0 = y0
	#ellipse.ellip0 = ellip0
	#ellipse.pa0 = pa0  #-90 < PA <= 90
	#ellipse.sma0 = sma0
	#ellipse.minsma = minsma
	#ellipse.maxsma = maxsma
	#ellipse.setParam('step',  step)
	ellipse.step = step
	#ellipse.linear = linear
	#ellipse.recenter = recenter
	#ellipse.xylearn = xylearn
	#ellipse.physical = physical
	##ellipse.unlearn() #set the algorithm control parameters for the ellipse task
	#ellipse.conver = conver
	#ellipse.minit = minit
	#ellipse.maxit = maxit
	#ellipse.hcenter = "no"
	#ellipse.hellip = "no"
	#ellipse.hpa = "no"
	#ellipse.wander = wander
	#ellipse.maxgerr = maxgerr
	#ellipse.olthresh = olthresh #if set to zero, the x0,y0 value found in geompar are used without questioning
	##ellipse.unlearn() #set image sampling parameters for the ellipse task
	#ellipse.integrmode = integrmode
	#ellipse.usclip = usclip
	#ellipse.lsclip = lsclip
	#ellipse.nclip = nclip
	#ellipse.fflag = fflag
	##ellipse.harmonics = "2 3 4 6 8 10 12"
	##ellipse.unlearn() #pset with paramters that define the magnitude scale
	#ellipse.mag0 = mag0
	#ellipse.refer = refer
	#ellipse.zerolevel = zerolevel
	#print ellipse.hcenter
	#print ellipse.hellip
	#print ellipse.hpa
	#print ellipse.ellip0
	#print ellipse.pa0
	#print ellipse.maxsma
	#print ellipse.minit
	#print ellipse.maxit
	#print ellipse.harmonics
	print ellipse.step
	ellipse.geompar.unlearn() #set the geometri parameters for the ellipse task
	#ellipse.geompar.x0 = x0
	#ellipse.geompar.y0 = y0
	#ellipse.geompar.ellip0 = ellip0
	#ellipse.geompar.pa0 = pa0  #-90 < PA <= 90
	#ellipse.geompar.sma0 = sma0
	#ellipse.geompar.minsma = minsma
	#ellipse.geompar.maxsma = maxsma
	ellipse.geompar.step = step
	#ellipse.geompar.linear = linear
	#ellipse.geompar.recenter = recenter
	#ellipse.geompar.xylearn = xylearn
	#ellipse.geompar.physical = physical
	#ellipse.controlpar.unlearn() #set the algorithm control parameters for the ellipse task
	#ellipse.controlpar.conver = conver
	#ellipse.controlpar.minit = minit
	#ellipse.controlpar.maxit = maxit
	#ellipse.controlpar.hcenter = hcenter
	#ellipse.controlpar.hellip = hellip
	#ellipse.controlpar.hpa = hpa
	#ellipse.controlpar.wander = wander
	#ellipse.controlpar.maxgerr = maxgerr
	#ellipse.controlpar.olthresh = olthresh #if set to zero, the x0,y0 value found in geompar are used without questioning
	#ellipse.samplepar.unlearn() #set image sampling parameters for the ellipse task
	#ellipse.samplepar.integrmode = integrmode
	#ellipse.samplepar.usclip = usclip
	#ellipse.samplepar.lsclip = lsclip
	#ellipse.samplepar.nclip = nclip
	#ellipse.samplepar.fflag = fflag
	#ellipse.samplepar.harmonics = "2 3 4 6 8 10 12"
	#ellipse.magpar.unlearn() #pset with paramters that define the magnitude scale
	#ellipse.magpar.mag0 = mag0
	#ellipse.magpar.refer = refer
	#ellipse.magpar.zerolevel = zerolevel
	#print ellipse.harmonics
	#print ellipse.samplepar.harmonics	

	ellipse(input_image, output, x0=x0, y=y0, ellip0=ellip0, pa0=pa0, sma0=sma0, step=step, olthresh=olthresh, maxsma=maxsma, harmonics='2 3 4 6', interactive=interactive)


def stsdas_analysis_isophote_bmodel(table, output, parent_image=None, backgr=0):
	'''
	build a model image from the results of isophotal analysis
	'''
	bmodel.unlearn()
	if parent_image is not None:
		bmodel.parent = parent_image
	
	bmodel(table=table, output=output, backgr=backgr)


def stsdas_analysis_isophote_cmodel(table, output, parent_image=None, backgr=0):
	'''
	build a model image from the results of isophotal analysis
	'''
	cmodel.unlearn()
	if parent_image is not None:
		cmodel.parent = parent_image
	
	cmodel(table=table, output=output, backgr=backgr)


if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_isomodel = ''
	parser.add_option('-m','--isophote_model', dest = 'isophote_model', type= 'string', default = def_isomodel, help='isophote model')

	def_isoimage = ''
	parser.add_option('-o', '--isophote_image',  dest = 'isophote_image', type= 'string', default = def_isoimage, help='output isophote image')

	def_x0 = 661 
	parser.add_option('-x','--x0', dest='x0', type=float, default=def_x0, help='the initial guess on the galaxy center; default: %s'%def_x0)

	def_y0 = 775
	parser.add_option('-y','--y0', dest='y0', type=float, default=def_y0, help='the initial guess on the galaxy center; default: %s'%def_y0)

	def_minsma = 1
	parser.add_option('--minsma', dest='minsma', type=float, default=def_minsma, help='min of the galaxy semimajor exis length; default: %s'%def_minsma)

	def_maxsma = 120
	parser.add_option('--maxsma', dest='maxsma', type=float, default=def_maxsma, help='max of the galaxy center semimajor axix length; default: %s'%def_maxsma)

	def_sma0 = 20
	parser.add_option('--sma0', dest='sma0', type=float, default=def_sma0, help='the initial starting point on semimajor axis; default: %s'%def_sma0)

	def_step = 0.1
	parser.add_option('--step', dest='step', type=float, default=def_step, help='see the document of isophote-ellipse; default: %s'%def_step)

	def_bkg = 0
	parser.add_option('--bkg', dest='bkg', type=float, default=def_bkg, help='the bkg value; default: %s'%def_bkg)

	options, remainder = parser.parse_args()
	input_image = options.input_image
	isophote_model = options.isophote_model
	isophote_image = options.isophote_image
	x0 = options.x0
	y0 = options.y0
	minsma = options.minsma
	maxsma = options.maxsma
	sma0 = options.sma0
	step = options.step
	bkg  = options.bkg

	stsdas_analysis_isophote_ellipse(input_image, isophote_model, x0, y0, minsma, maxsma, ellip0=0.2, pa0=45, sma0=sma0, step=step)
	#stsdas_analysis_isophote_bmodel(isophote_model, isophote_image, parent_image=input_image, backgr=bkg) 
	stsdas_analysis_isophote_cmodel(isophote_model, isophote_image, parent_image=input_image, backgr=bkg) 
