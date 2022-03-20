#! /usr/bin/python

import sys
sys.path.append('/home/asassn/ASASSN/Scripts/fitting_collections/scipy')


from fitting_functions import superflat_model_order1, superflat_model_order1_residual
from fitting_functions import superflat_model_order2, superflat_model_order2_residual
from fitting_functions import superflat_model_distance, superflat_model_distance_residual
from leastsq_fitting import scipy_leastsq


import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def superflat_calibration_result_display(data, model, params):

	
	x = data[:,0]
	y = data[:,1]
	z = data[:,4] - data[:,2]
	
	# Set up a regular grid of interpolation points
	xi, yi = np.linspace(np.min(x), np.max(x), 200), np.linspace(np.min(y), np.max(y), 200)
	xi, yi = np.meshgrid(xi, yi)
	
	# Interpolate
	rbf = interpolate.Rbf(x, y, z, function='linear')
	zi = rbf(xi, yi)
	
	plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
	plt.scatter(x, y, c=z)
	plt.colorbar()
	plt.show()
	
	
	
	zi_fit = np.array([model(xgrid, ygrid, *params) for xgrid,ygrid in zip(np.ravel(xi),np.ravel(yi))])
	zi_display = zi_fit.reshape(xi.shape)
	plt.imshow(zi_display, vmin=z.min(), vmax=z.max(), origin='lower',extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
	plt.scatter(x, y, c=z)
	plt.colorbar()
	plt.show()




