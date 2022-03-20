#! /usr/bin/python


import numpy as np



def flux2mag(flux,flt,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms'):
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
		lambda_eff_dict = {'W2':2056.6,
				   'M2':2246.4,
				   'W1':2581.0,
				   'U':3462.8,
				   'B':4379.7,
				   'V':5485.1,
				   'g':5200,
				   'r':6169.5,
				   'i':7496.6,
				   'z':9100,
				   'unit':'Angstroms'}
		
		f0_lambda = f0_nu /(3.34*10**4*lambda_eff_dict[flt]**2)	#erg/s/cm^2/A
		f0 = f0_lambda

	else:
		raise IOError('unrecognized mode for experssion of zero point of flux!')

	mag = -2.5*np.log10(flux/f0)

	return mag


def mag2flux(mag,flt,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms'):
	'''
	convert magnitude in AB system into flux

	if 'I_unit' == 'Jy' then F_nu in Jy unit; 
	if 'I_unit' == 'cgs' then F_nu in  unit of erg/s/cm^2/Hz; 

	F_lambda is in unit of erg/s/cm^2/A
	'''
	f0_nu = 3631	# Jy

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
        		   'z':9100,
        		   'unit':'Angstroms'}


	if mode == 'frequency':
		f0 = f0_nu*1e-23   #frequence F_\nu = 3631 Jy and 1 Jy = 10^-23 erg/s/cm^2/Hz
	elif mode == 'wavelength':
		f0_lambda = f0_nu /(3.34*10**4*lambda_eff_dict[flt]**2)	#erg/s/cm^2/A
		f0 = f0_lambda
	else:
		raise IOError('unrecognized mode for experssion of zero point of flux!')

	lam  = lambda_eff_dict[flt]
	flux = f0*10**(-0.4*mag)

	return lam,flux

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

	
