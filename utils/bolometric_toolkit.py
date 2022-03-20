#! /usr/bin/python

#This file contains some function for bolometric construction


import numpy as np
from numpy import exp
from astropy.table import Table

try:
	import lmfit
except:
	print "lmfit module import error; functions using lmfit can not be used"

from collections import OrderedDict
import matplotlib.pylab as plt


def galactic_reddening_NED(target,flt):
	'''
	Get the Galactic Extinction from NED 'http://ned.ipac.caltech.edu/'
	'''

	A_lambda = {'ASASSN-14il':{'U':0.096,
				   'B':0.080,	
				   'V':0.061,
				   'R':0.048,
				   'I':0.033,
				   'u':0.094,
				   'g':0.073,
				   'r':0.051,	
				   'i':0.038,	
				   'z':0.028,	
				   'J':0.016,	
				   'H':0.010,	
				   'K':0.007,	
				   'L':0.003},
		     'ASASSN-15pz':{'U':0.078,
				  'B':0.065,
			          'V':0.049,
				  'R':0.039,
				  'I':0.027,
				  'u':0.076,
				  'g':0.059,
				  'r':0.041,
				  'i':0.031,
				  'z':0.023,
				  'J':0.013,
				  'H':0.008,
				  'K':0.005},
		     'CSS161010':{'U':0.415,
				  'B':0.347,
			          'V':0.263,
				  'R':0.208,
				  'I':0.144,
				  'u':0.406,
				  'g':0.317,
				  'r':0.219,
				  'i':0.163,
				  'z':0.121,
				  'J':0.068,
				  'H':0.043,
				  'K':0.029},
			'SN2017dpu':{'U':0.134,
				   'B':0.112,
				   'V':0.085,
				   'r':0.07,
				   'i':0.052}}
	
	return A_lambda[target][flt]


def galactic_dereddening_swift_bands(flt,AV,spec_type='BB_10000K'):
	'''
	
	The value for U B V r i from NED
	
	For swift bands the extinction are adopted from 
	'http://iopscience.iop.org/0004-637X/721/2/1608/suppdata/apj343574t14_mrt.txt'

	which is a table for R_X, where A_X = R_X*E(B-V) in filter X = W2, M2, W1, U, B, V

	'''
	
	RX = {'SN_1992A' : {'W2': 6.44, 'M2':8.06 , 'W1':5.45 , 'U':4.91 , 'B':4.16 , 'V':3.16},
	      'SN_1994I' : {'W2': 6.32, 'M2':8.30 , 'W1':5.35 , 'U':4.88 , 'B':4.10 , 'V':3.14}, 
	      'SN_1999em': {'W2': 7.63, 'M2':8.37 , 'W1':6.18 , 'U':5.00 , 'B':4.16 , 'V':3.16}, 
	      'BB_2000K' : {'W2': 2.65, 'M2':4.90 , 'W1':3.33 , 'U':4.65 , 'B':3.92 , 'V':3.09}, 
	      'BB_2200K' : {'W2': 2.83, 'M2':5.32 , 'W1':3.52 , 'U':4.71 , 'B':3.95 , 'V':3.10}, 
	      'BB_2400K' : {'W2': 3.03, 'M2':5.79 , 'W1':3.71 , 'U':4.75 , 'B':3.97 , 'V':3.11}, 
	      'BB_2600K' : {'W2': 3.23, 'M2':6.22 , 'W1':3.91 , 'U':4.79 , 'B':3.98 , 'V':3.11}, 
	      'BB_2800K' : {'W2': 3.45, 'M2':6.58 , 'W1':4.12 , 'U':4.81 , 'B':4.00 , 'V':3.12}, 
	      'BB_3000K' : {'W2': 3.70, 'M2':6.86 , 'W1':4.32 , 'U':4.83 , 'B':4.01 , 'V':3.12}, 
	      'BB_3200K' : {'W2': 3.97, 'M2':7.07 , 'W1':4.53 , 'U':4.85 , 'B':4.02 , 'V':3.13}, 
	      'BB_3400K' : {'W2': 4.26, 'M2':7.25 , 'W1':4.72 , 'U':4.86 , 'B':4.04 , 'V':3.13}, 
	      'BB_3600K' : {'W2': 4.57, 'M2':7.39 , 'W1':4.90 , 'U':4.87 , 'B':4.05 , 'V':3.13}, 
	      'BB_3800K' : {'W2': 4.87, 'M2':7.51 , 'W1':5.07 , 'U':4.88 , 'B':4.06 , 'V':3.13}, 
	      'BB_4000K' : {'W2': 5.18, 'M2':7.61 , 'W1':5.22 , 'U':4.89 , 'B':4.06 , 'V':3.14}, 
	      'BB_4200K' : {'W2': 5.47, 'M2':7.70 , 'W1':5.35 , 'U':4.90 , 'B':4.07 , 'V':3.14}, 
	      'BB_4400K' : {'W2': 5.74, 'M2':7.77 , 'W1':5.48 , 'U':4.91 , 'B':4.08 , 'V':3.14}, 
	      'BB_4600K' : {'W2': 6.00, 'M2':7.84 , 'W1':5.59 , 'U':4.92 , 'B':4.08 , 'V':3.14}, 
	      'BB_4800K' : {'W2': 6.23, 'M2':7.90 , 'W1':5.69 , 'U':4.92 , 'B':4.09 , 'V':3.14}, 
	      'BB_5000K' : {'W2': 6.44, 'M2':7.96 , 'W1':5.78 , 'U':4.93 , 'B':4.10 , 'V':3.15}, 
	      'BB_5200K' : {'W2': 6.63, 'M2':8.01 , 'W1':5.86 , 'U':4.93 , 'B':4.10 , 'V':3.15}, 
	      'BB_5600K' : {'W2': 6.95, 'M2':8.09 , 'W1':6.00 , 'U':4.94 , 'B':4.11 , 'V':3.15}, 
	      'BB_6000K' : {'W2': 7.21, 'M2':8.17 , 'W1':6.12 , 'U':4.95 , 'B':4.12 , 'V':3.15}, 
	      'BB_6600K' : {'W2': 7.50, 'M2':8.25 , 'W1':6.27 , 'U':4.96 , 'B':4.13 , 'V':3.15}, 
	      'BB_7000K' : {'W2': 7.64, 'M2':8.30 , 'W1':6.35 , 'U':4.97 , 'B':4.13 , 'V':3.15}, 
	      'BB_8000K' : {'W2': 7.89, 'M2':8.40 , 'W1':6.52 , 'U':4.98 , 'B':4.14 , 'V':3.16}, 
	      'BB_9000K' : {'W2': 8.05, 'M2':8.47 , 'W1':6.64 , 'U':4.99 , 'B':4.15 , 'V':3.16}, 
	      'BB_10000K': {'W2': 8.14, 'M2':8.53 , 'W1':6.74 , 'U':5.00 , 'B':4.16 , 'V':3.16}, 
	      'BB_12000K': {'W2': 8.25, 'M2':8.60 , 'W1':6.89 , 'U':5.01 , 'B':4.17 , 'V':3.16}, 
	      'BB_13000K': {'W2': 8.27, 'M2':8.63 , 'W1':6.94 , 'U':5.01 , 'B':4.17 , 'V':3.17}, 
	      'BB_14000K': {'W2': 8.29, 'M2':8.65 , 'W1':6.99 , 'U':5.01 , 'B':4.18 , 'V':3.17}, 
	      'BB_15000K': {'W2': 8.31, 'M2':8.67 , 'W1':7.03 , 'U':5.02 , 'B':4.18 , 'V':3.17}, 
	      'BB_20000K': {'W2': 8.33, 'M2':8.73 , 'W1':7.17 , 'U':5.03 , 'B':4.19 , 'V':3.17}, 
	      'BB_25000K': {'W2': 8.33, 'M2':8.76 , 'W1':7.25 , 'U':5.03 , 'B':4.19 , 'V':3.17}, 
	      'BB_30000K': {'W2': 8.33, 'M2':8.78 , 'W1':7.30 , 'U':5.04 , 'B':4.20 , 'V':3.17}, 
	      'BB_35000K': {'W2': 8.33, 'M2':8.79 , 'W1':7.33 , 'U':5.04 , 'B':4.20 , 'V':3.17}}



	A_X_mag = AV/RX[spec_type]['V']*RX[spec_type][flt]
	
	return A_X_mag
		

	

def Vega_AB_mag_convertion(mag,flt,mode,direction='Vega2AB'):
	'''
	transform the magnitude from vega system to AB magnitude system

	For mode == 'Swift', ref to "http://swift.gsfc.nasa.gov/analysis/uvot_digest/zeropts.html"

	For mode == 'Bessell', refer to 'http://www.astronomy.ohio-state.edu/~martini/usefuldata.html'
	These data are mostly from Blanton et al. (2007)
	

	INPUT:
	mode: 'Swift' or 'Bessell'
	direction: 'Vega2AB' or 'AB2Vega'

	'''

	if mode not in ['Swift','Bessell']:
		raise IOError('Only Swift and Bessell bands are supported now!')
	

	#Values for m_AB - m_Vega
	offset_dict = {'Swift':{'V':-0.01,
				'B':-0.13,
				'U':+1.02,
				'W1':+1.51,
				'M2':+1.69,
				'W2':+1.73,
				'White':+0.80},
		       'Bessell':{'U':+0.79,
				  'B':-0.09,
				  'V':+0.02,
				  'R':+0.21,
				  'I':+0.45,
				  'J':+0.91,
				  'H':+1.39,
				  'K':+1.85,
				  'u':+0.91,
				  'g':-0.08,
				  'r':+0.16,
				  'i':+0.37,
				  'z':+0.54,
				  'Y':+0.634}}
	if direction == 'Vega2AB':
		mag_new = mag + offset_dict[mode][flt]
	elif direction == 'AB2Vega':
		mag_new = mag - offset_dict[mode][flt]
	else:
		raise IOError('directions Vega2AB or AB2Vega are supported!')

	return mag_new

	


def mag2flux(mag,flt, bp_mp =None, mode='wavelength',I_unit='cgs',wavelength_units='Angstroms'):
	'''
	convert magnitude in AB system into flux

	if 'I_unit' == 'Jy' then F_nu in Jy unit; 
	if 'I_unit' == 'cgs' then F_nu in  unit of erg/s/cm^2/Hz; 

	F_lambda is in unit of erg/s/cm^2/A
	'''
	f0_nu = 3631	# Jy

	lambda_eff_dict = {'W2':2030.5,
			   'M2':2228.1,
			   'W1':2589.1,
			   'UVOT_U':3501.2,
			   'UVOT_B':4328.6,
			   'UVOT_V':5402.1,
        		   'U':3462.8,
        		   'B':4379.7,
        		   'V':5485.1,
        		   'R':6400,
        		   'I':7900,
        		   'J':12600,
        		   'H':16000,
        		   'K':22200,
        		   'g':5200,
        		   'r':6169.5,
        		   'i':7496.6,
        		   'z':9100,
        		   'unit':'Angstroms'}


	if mode == 'frequency':
		f0 = f0_nu*1e-23   #frequence F_\nu = 3631 Jy and 1 Jy = 10^-23 erg/s/cm^2/Hz
	elif mode == 'wavelength':

		if bp_mp is None:
			bp_mp = lambda_eff_dict[flt]

		f0_lambda = f0_nu /(3.34*10**4*bp_mp**2)	#erg/s/cm^2/A
		f0 = f0_lambda
	else:
		raise IOError('unrecognized mode for experssion of zero point of flux!')

	if bp_mp is None:
		lam  = lambda_eff_dict[flt]
	else:
		lam = bp_mp

	flux = f0*10**(-0.4*mag)

	return lam,flux




def flux2mag(flux,flt, bp_mp =None,  mode='wavelength',I_unit='cgs',wavelength_units='Angstroms'):
	'''
	convert magnitude in AB system into flux

	if 'I_unit' == 'Jy' then F_nu in Jy unit; 
	if 'I_unit' == 'cgs' then F_nu in  unit of erg/s/cm^2/Hz; 

	F_lambda is in unit of erg/s/cm^2/A
	'''
	f0_nu = 3631	# Jy

	if mode == 'flux':
		f0 = f0_nu*1e-23   #frequence F_\nu = 3631 Jy and 1 Jy = 10^-23 erg/s/cm^2/Hz
	elif mode == 'wavelength':
		lambda_eff_dict = {'W2':2030.5,
				   'M2':2228.1,
				   'W1':2589.1,
				   'UVOT_U':3501.2,
				   'UVOT_B':4328.6,
				   'UVOT_V':5402.1,
				   'U':3462.8,
				   'B':4379.7,
        		   	   'R':6400,
        		           'I':7900,
        		           'J':12600,
        		            'H':16000,
        		            'K':22200,
        		            'g':5200,
				   'V':5485.1,
				   'r':6169.5,
				   'i':7496.6,
				   'z':9100,
				   'unit':'Angstroms'}
		
		if bp_mp is None:
			bp_mp = lambda_eff_dict[flt]

		f0_lambda = f0_nu /(3.34*10**4*bp_mp**2)	#erg/s/cm^2/A
		f0 = f0_lambda

	else:
		raise IOError('unrecognized mode for experssion of zero point of flux!')

	mag = -2.5*np.log10(flux/f0)

	return mag





	

unitdict = {'cgs':{'h':6.626068e-27,
                   'k':1.3806503e-16,
                   'c':2.99792458e10,
                   'mh':1.67262158e-24 * 1.00794,
                   'length':'cm'},
            'mks':{'h':6.626068e-34,
                   'k':1.3806503e-23,
                   'c':2.99792458e8,
                   'mh':1.67262158e-27 * 1.00794,
                   'length':'m'}
            }

wavelength_dict = {'meters':1.0,'m':1.0,
                   'centimeters':1e-2,'cm':1e-2,
                   'millimeters':1e-3,'mm':1e-3,
                   'nanometers':1e-9,'nm':1e-9,
                   'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
                   'kilometers':1e3,'km':1e3,
                   'angstroms':1e-10,'A':1e-10,'Angstroms':1e-10,
                   }

def blackbody_wavelength(lam,temperature, scale=1.0,units='cgs',wavelength_units='Angstroms', normalize=max, beta=0):
	'''
	INPUTS:
		lam:
		temperature:

	'''
	# load constants in desired units
	h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']
	
	# converta lambd to cm/m
	lam = lam * wavelength_dict[wavelength_units] / (1e-2 if units=='cgs' else 1)
	
	I = 2*h*c**2 / lam**5 * (exp(h*c/(k*temperature*lam)) - 1)**-1
	
	if normalize and hasattr(I,'__len__'):
	    if len(I) > 1:
	        return I/normalize(I) * scale
	    else:
	        return I * scale
	else:
	    return I * scale




def fit_blackbody_lmfit(xdata, flux, guesses=(0,0), err=None,blackbody_function=blackbody_wavelength, quiet=True, **kwargs):
	"""
	Parameters
	----------
	xdata : array
	    Array of the X-values (frequency, wavelength) of the data
	flux : array
	    The fluxes corresponding to the xdata values
	guesses : (Temperature,Scale) or (Temperature,Beta,Scale)
	    The input guesses.  3 parameters are used for greybody
	    fitting, two for temperature fitting.
	blackbody_function: function
	    Must take x-axis (e.g. frequency), temperature, scale, and then
	    optionally beta args
	quiet : bool
	    quiet flag passed to mpfit
	
	kwargs are past to blackbody function
	
	Examples
	--------
	>>> wavelength = np.array([20,70,160,250,350,500,850,1100])
	>>> flux = modified_blackbody_wavelength(wavelength, 15, beta=1.75,
	        wavelength_units='microns', normalize=False, logN=22, logscale=16)
	>>> err = 0.1 * flux
	>>> flux += np.random.randn(len(wavelength)) * err
	>>> tguess, bguess, nguess = 20.,2.,21.5
	>>> lm = fit_blackbody_lmfit(wavelength, flux, err=err,
	         blackbody_function=modified_blackbody_wavelength, logscale=16,
	         guesses=(tguess,bguess,nguess),
	         wavelength_units='microns')
	>>> print lm.params
	
	>>> # If you want to fit for a fixed beta, do this:
	>>> parameters = lmfit.Parameters(OrderedDict([ (n,lmfit.Parameter(x)) for n,x
	        in zip(('T','beta','N'),(20.,2.,21.5)) ]))
	>>> import lmfit
	>>> parameters['beta'].vary = False
	>>> lm = fit_blackbody_lmfit(wavelength, flux, err=err,
	         blackbody_function=modified_blackbody_wavelength, logscale=16,
	         guesses=parameters,
	         wavelength_units='microns')
	>>> print lm.params
	"""
	
	def lmfitfun(x,y,err):
	    if err is None:
	        def f(p): return (y-blackbody_function(x, *[p[par].value for par in p],normalize=False, **kwargs))
	    else:
	        def f(p): return (y-blackbody_function(x, *[p[par].value for par in p],normalize=False, **kwargs))/err
	    return f
	
	if not isinstance(guesses,lmfit.Parameters):
	    guesspars = lmfit.Parameters( OrderedDict([ (n,lmfit.Parameter(value=x,name=n)) for n,x in zip(('T','beta','N'),guesses) ]))
	else:
	    guesspars = guesses
	
	minimizer = lmfit.minimize( lmfitfun(xdata,np.array(flux),err),
	        guesspars)
	
	return minimizer


def read_and_prepare_mags(data,flts):
	'''
	Read data into data dictionary 
	'data' format: JD mag1 magerr1 mag2 magerr2 mag3 magerr3 mag4 magerr4...
	flts =[flt1 flt2 flt3 flt4 ...]
	'''

	data_dict = {}
	data_err_dict = {}

	for i,flt in enumerate(flts):

		data_dict[flt] = data[2*i+1]
		data_err_dict[flt] = data[2*i+2]
	
	return data_dict,data_err_dict
		



def get_SED_AB_mag(sed_input_data_file, target_name  = None, A_V=None, correct_gal_extinction=False):
	'''
	INPUTS:
		sed_input_data_file: in the following format

			#flt    mag    magerr    sys
			 B      13.5   0.04      Vega
			 g      14.5   0.03      AB
			 r      14.0   0.00      AB
	'''

	lambda_eff_dict = {'W2':2056.6,
        		   'M2':2246.4,
        		   'W1':2581.0,
        		   'U':3462.8,
        		   'B':4379.7,
        		   'V':5485.1,
        		   'R':6400,
        		   'I':7900,
        		   'J':12600,
        		   'H':16000,
        		   'K':22200,
        		   'g':5200,
        		   'r':6169.5,
        		   'i':7496.6,
        		   'z':9100}


	data = Table.read(sed_input_data_file,format='ascii')

	flts = []
	lams = []
	mags = []
	magerrs = []

	swift_flts = ['W1','M2','W2']
	
	for fltdata in data:
		flt    = fltdata['flt']
		flts.append(flt)
		lams.append(lambda_eff_dict[flt])
		magsys = fltdata['sys']
		mag    = fltdata['mag']

		if correct_gal_extinction:

			if flt in swift_flts:
				if A_V is None:
					A_V = galactic_reddening_NED(target_name, 'V')

				A_mag = galactic_dereddening_swift_bands(flt, AV, spec_type='BB_10000K')
				mag = mag - A_mag
			else:
				if target_name is not None:
					try:
						A_mag = galactic_reddening_NED(target_name, flt)
						mag = mag - A_mag
					except:
						print "extinction for %s not implemented in bolometric_toolkit.py yet"%target_name

		magerr = fltdata['magerr']
 
		if magsys == 'AB':
			mags.append(mag)
			magerrs.append(magerr)
		else:
			if flt in swift_flts:
				mag_AB = Vega_AB_mag_convertion(mag,flt,'Swift',direction='Vega2AB')
			else:
				mag_AB = Vega_AB_mag_convertion(mag,flt,'Bessell',direction='Vega2AB')
			mags.append(mag)
			magerrs.append(magerr)


	output_data = [flts, lams, mags, magerrs]
	output_table = Table(output_data, names=('flt','lam','mag','magerr'))
				
			
	return output_table



def get_SED_flux(sed_mags_table, target_name  = None):
	'''
	INPUTS:
		sed_mags_table:
	
	'''
	lams     = []
	fluxs    = []
	fluxerrs = []

	for fltdata in sed_mags_table:
		flt  = fltdata['flt']
		mag  = fltdata['mag']
		magerr = fltdata['magerr']

		lamb, flux = mag2flux(mag,flt,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms')
		lamb, flux_upper = mag2flux(mag-magerr, flt)
		fluxerr = flux_upper - flux

		lams.append(lamb)
		fluxs.append(flux)
		fluxerrs.append(fluxerr)

	output_data = [lams, fluxs, fluxerrs]
	output_table = Table(output_data, names=('lam','flux', 'fluxerr'))

	return output_table




def show_SED(photdata_dict, ):
	'''
	
	for 

	mag2flux(mag,flt,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms')

	'''


if __name__ == '__main__':


	#Basic information known:
	SN_name = 'ASASSN-14il'
	z = 0.021989
	D = 86	#Mpc 	??
	Derr = 0
	AV = galactic_reddening_NED(SN_name,'V')
	Lsun = 3.83*10**33	#erg/s
	t0 = 2456931.500	#discovery date

	phot_rets = np.loadtxt('ASASSN-14il_bolometri_lc_input.txt')

	mag_host_dict=    {'W2':16.1,	'M2':16.11,  'W1':16.09, 'U':15.86,	'B':16.00,	'V':15.48,  'r':15.24,  'i':15.03}
	magerr_host_dict= {'W2':0.03,	'M2':0.03,   'W1':0.03,	 'U':0.03,	'B':0.05,	'V':0.05,   'r':0.08, 	'i':0.10}

	lambda_eff_dict = {'W2':2056.6,
        		   'M2':2246.4,
        		   'W1':2581.0,
        		   'U':3462.8,
        		   'B':4379.7,
        		   'V':5485.1,
        		   'r':6169.5,
        		   'i':7496.6,
        		   'unit':'Angstroms'}

	flts = ['W2','M2','W1','U','B','V','r','i']
	flts_Vega =  ['W2','M2','W1','U','B','V']
	flts_Swift = ['W2','M2','W1','U']
	flts_Bessell = ['B','V','r','i']
	
	colors = ['b','r','g','c','m','y','k']


	jds   = np.array([])
	Ts    = np.array([])
	Terrs = np.array([])	
	Rs    =	np.array([])	
	Rerrs = np.array([])	
	Ls    =	np.array([])	
	Lerrs =	np.array([])		
	chi2s = np.array([])

	
	lengends = []
		
	for i,phot_ret in enumerate(phot_rets):
	
		jd = phot_ret[0]
	
		print jd
		mag_tot_dict,magerr_tot_dict = read_and_prepare_mags(phot_ret,flts)
	
		wavelengths = np.array([])
		fluxs = np.array([])
		fluxerrs = np.array([])


		

		for flt in flts:
			lambda_flt = lambda_eff_dict[flt]

			if flt in flts_Swift:
				A_flt = galactic_dereddening_swift_bands(flt,AV,spec_type='BB_10000K')	# The spectra type is rough estimation
			elif flt in flts_Bessell:
				A_flt = galactic_reddening_NED(SN_name,flt)
			else:
				raise IOError('get galactic reddening wrong! Please Check ...')
	
			mag_tot  = mag_tot_dict[flt]
			if mag_tot == 99.99:
				continue
	
			magerr_tot = magerr_tot_dict[flt]
			mag_host = mag_host_dict[flt]
			magerr_host = magerr_host_dict[flt]

			if mag_tot > mag_host:
				continue			
	
			if flt in flts_Vega:
				if flt in flts_Swift:
					mag_tot  = Vega_AB_mag_convertion(mag_tot,flt,mode='Swift',direction='Vega2AB')
					mag_host = Vega_AB_mag_convertion(mag_host,flt,mode='Swift',direction='Vega2AB')	
				elif flt in flts_Bessell:
					mag_tot  = Vega_AB_mag_convertion(mag_tot,flt,mode='Bessell',direction='Vega2AB')         	
					mag_host = Vega_AB_mag_convertion(mag_host,flt,mode='Bessell',direction='Vega2AB')	
				else:
					raise ValueError('Get AB magnitude error, please check...')
	
			lamb_flt,flux_tot  = mag2flux(mag_tot-A_flt, flt)
			lamb_flt,flux_host = mag2flux(mag_host-A_flt,flt)
			flux_flt = flux_tot - flux_host
	
			fluxerr_flt = np.sqrt(magerr_tot**2 + magerr_host**2)*1.086*flux_flt
	
			

			wavelengths = np.append(wavelengths,lambda_flt)
			fluxs = np.append(fluxs,flux_flt)
			fluxerrs = np.append(fluxerrs, fluxerr_flt)
			
			
		if len(wavelengths) <3:
			continue
			
		lm_ret = fit_blackbody_lmfit(wavelengths, fluxs, guesses=(9000,1e-30), err=fluxerrs)
		params = lm_ret.params
		print lm_ret.params
		T = params['T'].value
		Terr = params['T'].stderr
		print T
	
		scale = params['beta'].value
		scale_err = params['beta'].stderr
		print scale
		
		scale_keep = scale	
		
		scale = scale*1e8	#the factor is from the difference between erg/s/cm^2/A and erg/s/cm^2/cm
		scale_err = scale_err*1e8
		Mpc2cm =  3.0857*10**24
		D_cm = D*Mpc2cm
		Derr_cm = Derr*Mpc2cm
		
		R = np.sqrt(scale/np.pi)*D_cm	# cm
		Rerr = np.sqrt(np.pi) * np.sqrt( (scale**(1./2)*Derr_cm)**2  + (D_cm*1./2*scale**(-1./2)*scale_err)**2 )
		

		sigma = 5.76*10**-5	#erg/s/cm^2/deg^-4
		L = 4*(np.pi)**2*R**2*sigma*T**4
		Lerr = 4*(np.pi)**2*sigma*np.sqrt( (2*R*T**4*Rerr)**2 + (4*T**3*R**2*Terr)**2 )


		wavelengths_plot = np.linspace(100,10000,100)
		fluxs_bb  = blackbody_wavelength(wavelengths_plot,T, scale=scale_keep, normalize=False)

		jds = np.append(jds,jd)
		Ts = np.append(Ts,T)
		Terrs = np.append(Terrs,Terr)
	
		Rs = np.append(Rs,R)
	        Rerrs = np.append(Rerrs,Rerr)

		Ls = np.append(Ls,L)
	        Lerrs = np.append(Lerrs,Lerr)
	
		#print fluxs_bb

#		if i>30:
#			break
#		
#		if np.mod(i,4) != 0:
#			continue
#		color = colors[np.mod(i,7)]
#		plt.errorbar(wavelengths,fluxs,yerr = fluxerrs,fmt=color+'o')
#		plt.plot(wavelengths_plot,fluxs_bb,color)
#		lengend_t = np.str(np.round(jd - t0))
#		lengends.append(lengend_t)	
#		lengends.append(lengend_t)	
	
#	plt.xscale('log')
#	plt.yscale('log')
#	plt.legend(lengends,loc=2)
#	plt.xlabel('$wavelength(\AA)$')
#	plt.ylabel('$F_\lambda(erg/s/cm^2/\AA)$')
#	plt.show()

	ret_out = np.hstack((jds.reshape(len(jds),1),Ts.reshape(len(jds),1),Terrs.reshape(len(jds),1),Rs.reshape(len(jds),1),Rerrs.reshape(len(jds),1),Ls.reshape(len(jds),1),Lerrs.reshape(len(jds),1)))
	print ret_out
	fid = open('bbfit_result_ping.txt','wt')
	np.savetxt(fid,ret_out,fmt='%12.5f %d %d %.2e %.2e %.2e %.2e')
	fid.close()


	plt.errorbar(jds,Ts,yerr=Terrs,fmt='o')
