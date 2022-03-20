#! /usr/bin/python

import sys
sys.path.append('/home/asassn/ASASSN/Scripts/fitting_collections/scipy')


from fitting_functions import superflat_model_order1, superflat_model_order1_residual
from fitting_functions import superflat_model_order2, superflat_model_order2_residual
from fitting_functions import superflat_model_distance, superflat_model_distance_residual
from leastsq_fitting import scipy_leastsq

from superflat_result_display import superflat_calibration_result_display

import numpy as np
import matplotlib.pylab as plt

def superflat_calibration(x,y,obs_mags, std_mags, obs_magerr, std_magerr, sn_x, sn_y, model='poly1', param0=None, display=False):
	'''
	calibration using superflat method which take spatial variance into account

	INPUTS:
		x:
		y:
		obs_mags:
		std_mags:
		obs_magerrs:
		std_magerrs:
		sn_x:
		sn_y: 
		model: 'poly1', 'poly2', 'distance'
		param0: initial guess for parameters
	'''

	if param0 is None:
		if model == 'poly1':
			param0 = [0,0,0]
		elif model == 'poly2':
			param0 = [0,0,0,0,0,0]
		elif model == 'distance':
			param0 = [0,0,0,0]
		else:
			raise IOError("model %s not supported"%model)
	
	
	errs = np.sqrt(obs_magerr**2+std_magerr**2)
	x = x[errs!=0]
	y = y[errs!=0]
	std_mags = std_mags[errs!=0]
	obs_mags = obs_mags[errs!=0]
	std_magerr = std_magerr[errs !=0]
	obs_magerr = obs_magerr[errs !=0]

	errs = errs[errs != 0]

	args = (x,y, std_mags-obs_mags,errs)


	if display:
		plt.errorbar(std_mags, std_mags-obs_mags, yerr=np.sqrt(obs_magerr**2+std_magerr**2), fmt='o')
		plt.xlabel('std_mags')
		plt.ylabel('std_mag - obs_mag')
		plt.show()

	
	if model == 'poly1':
		residual = superflat_model_order1_residual
		flatmodel = superflat_model_order1
	elif model == 'poly2':
		residual = superflat_model_order2_residual
		flatmodel = superflat_model_order2
	elif model == 'distance':
		residual = superflat_model_distance_residual
		flatmodel = superflat_model_distance
	else:
		raise IOError("model %s not supported"%model)


	fitting_out = scipy_leastsq(residual, param0, args)

	params_fit = fitting_out[0]
	print params_fit

#	c1 = params_fit[0]
#	c2 = params_fit[1]
#	c3 = params_fit[2]

	offsets = flatmodel(x, y, *params_fit)

	if display:
		plt.errorbar(std_mags, std_mags-obs_mags-offsets, yerr=np.sqrt(obs_magerr**2+std_magerr**2), fmt='o')
		plt.show()

		
		caldata = np.hstack((x.reshape(len(x),1),y.reshape(len(x),1),obs_mags.reshape(len(x),1),obs_magerr.reshape(len(x),1),std_mags.reshape(len(x),1),std_magerr.reshape(len(x),1), ))

		superflat_calibration_result_display(caldata, flatmodel, params_fit)


	offset_sn = flatmodel(sn_x, sn_y, *params_fit)

	return offset_sn


if __name__ == "__main__":

	import sys


	import optparse
	parser = optparse.OptionParser()

	def_model = 'poly1'
	parser.add_option('-m','--model', dest="model", type="string", default=def_model, help="the superflat model; the default is [%s]"%def_model)

	def_inputfile = ''
	parser.add_option('-i','--inputfile', dest="inputfile", type="string", default=def_inputfile, help="the input calibration file; the default is [%s]"%def_inputfile)

	def_cols = '0,1,2,3,6,7'
        parser.add_option('-c','--cols',dest="cols", type="string", default=def_cols, help="the column indexs in input calibration file which correspond to img_x,img_y, std_m, std_merr, obs_m, obs_merr [%s]"%def_cols)

	def_snx = 0
        parser.add_option('-x','--snx',dest="snx", type="float", default=def_snx, help="sn position x on the reference image [%s]"%def_snx)

	def_sny = 0
        parser.add_option('-y','--sny',dest="sny", type="float", default=def_sny, help="sn position y on the reference image [%s]"%def_sny)

	def_params0 = None
	parser.add_option('-p','--pms0', dest='params0', type="string", default=def_params0, help="the initiall guess for the model paramters [%s]"%def_params0)

	options, remainder = parser.parse_args()

	model = options.model
	inputfile = options.inputfile
	cols = options.cols

	cols_int = [int(col) for col in cols.split(',')]
	x_col = cols_int[0]
	y_col = cols_int[1]
	stdm_col = cols_int[2]
	stdmerr_col = cols_int[3]
	obsm_col = cols_int[4]
	obsmerr_col = cols_int[5]

	snx = options.snx
	sny = options.sny
	
	pms0 = options.params0
	if pms0 is not None:
		pms0 = [float(pm0) for pm0 in pms0.split(',')]
	
	data = np.loadtxt(inputfile)
	
	x = data[:,x_col]
	y = data[:,y_col]
	
	obs_mags = data[:,obsm_col]
	obs_magerr = data[:,obsmerr_col]

	std_mags = data[:,stdm_col]
	std_magerr = data[:,stdmerr_col]

#	sn_x = 492
#	sn_y = 498

	offset_sn = superflat_calibration(x,y,obs_mags, std_mags, obs_magerr, std_magerr, snx, sny, model=model, param0=pms0, display=True)
	print offset_sn
