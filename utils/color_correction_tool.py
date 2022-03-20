#! /usr/bin/python

import os
import ASASSN_photometry
from kmpfit_errorsinXandYPlot import kmp_linear_fitting
from common import sigma_clipping

import matplotlib.pylab as plt
import numpy as np


from matplotlib.pyplot import figure, show, rc
from scipy.odr import Data, Model, ODR, RealData, odr_stop



def compare_plot_offsetmag2color(SN1,img1,SN2,img2,flt1=None,flt2=None,img3=None, color_mode='color_S1', mag_offset1_flt1=0, mag_offset2_flt1=0,  mag_offset_colorflt=0, trim_img1 = None, trim_img2 = None):
	'''
	suppose you got two systems then two set magnitudes, S1 and S2, 
	one of which (let's say S1) is considered as standard system 
	To characterize the new system S2, you can plot 
	( m_flt1_S1 - m_flt1_S2 ) versus ( m_flt1_S1 - m_flt2_S1 )
	or 
	( m_flt1_S1 - m_flt1_S2 ) versus ( m_flt1_S2 - m_flt2_S2 )
	
	above flt1 is the filter of interest and flt2 is another filter 
	to derive color 


	Inputs:
		SN1 and SN2: 		ASASSN photometry object
		img1 and img2: 		image key for in SN1 and SN2
		flt1 and flt2: 		flt1 is the filter of interest and flt2 is another filter to derive color
		color_mode: 			'color_S1' or 'color_S2'
		      			'color_S1' uses color derived from system S1
		      			'color_S2' uses color derived from system S2
		mag_offset1_flt1:	add a vaule to S1 magnitude in flt1
		mag_offset2_flt1:	add a value to S2 magnitude in flt1
		mag_offset_colorflt:	add a vaule to magnitude in filter which is used to derive colors

		trim_img1: None or 4-value list [xl,xu,yl,yu]
		trim_img2:


		      
	'''


	snname1 = SN1.current_sn
	snname2 = SN2.current_sn
	if snname1 != snname2:
		raise IOError('working on two different target...')

	S1 = SN1.current_telescope
	S2 = SN2.current_telescope
	
	#-----------------------------------------------------------
	#prepare photometry results for stars in system S1
	imgkey1 = img1.split('.')[0]
	photmethod = SN1.photometry_method
	if photmethod == 'psfphot':
		photretfile = imgkey1+'.psfphot'
		photret = os.path.join(SN1.psf_photometry_dir, photretfile)
	else:
		photretfile = imgkey1+'.apphot'
		photret = os.path.join(SN1.aperture_photometry_dir, photretfile)
		
	xymags1 = np.loadtxt(photret)

	#trim the match source list if needed
	xs_S1 = xymags1[:,0]
	ys_S1 = xymags1[:,1]

	if trim_img1 is not None:
		xl = trim_img1[0]
		xu = trim_img1[1]
		yl = trim_img1[2]
		yu = trim_img1[3]

		x_trim_mask = np.logical_and(xs_S1>xl, xs_S1<xu)	
		y_trim_mask = np.logical_and(ys_S1>yl, ys_S1<yu)
		trim_mask   = np.logical_and(x_trim_mask, y_trim_mask)
		
		xymags1 = xymags1[trim_mask,:]

	#prepare magnitudes 
	xymags1[:,2] = xymags1[:,2] + mag_offset1_flt1


	#----------------------------------------------------------------
	#prepare photometry results for stars in system S2
	imgkey2 = img2.split('.')[0]
	photmethod = SN2.photometry_method
	if photmethod == 'psfphot':
		photretfile = imgkey2 +'.psfphot'
		photret = os.path.join(SN2.psf_photometry_dir, photretfile)
	else:
		photretfile = imgkey2 +'.apphot'
		photret = os.path.join(SN2.aperture_photometry_dir, photretfile)
		
	xymags2 = np.loadtxt(photret)

	#trim the match source list if needed
	xs_S2 = xymags2[:,0]
	ys_S2 = xymags2[:,1]

	if trim_img2 is not None:
		xl = trim_img2[0]
		xu = trim_img2[1]
		yl = trim_img2[2]
		yu = trim_img2[3]

		x_trim_mask = np.logical_and(xs_S2>xl, xs_S2<xu)	
		y_trim_mask = np.logical_and(ys_S2>yl, ys_S2<yu)
		trim_mask   = np.logical_and(x_trim_mask, y_trim_mask)
		
		xymags2 = xymags2[trim_mask,:]

	#prepare magnitudes 
	xymags2[:,2] = xymags2[:,2] + mag_offset2_flt1



	#---------------------------------------------------------------
	#prepare photometry results to derive colors of stars
	if color_mode == 'color_S1':
		SN3 = SN1
		trim_img3 = trim_img1
	elif color_mode == 'color_S2':
		SN3 = SN2
		trim_img3 = trim_img2
	else:
		raise IOError("invalid input for color_mode. only color_S1 and color_S2 are supported!")


	if img3 is not None:
		imgkey3 = img3.split('.')[0]	
		flt2 = SN3.photometry_info[imgkey3+'.fits']['flt']
	elif flt2 is not None:	
		tplimgs  = SN3._photometry__find_template_imagekey()
		tplimg_flt2 = tplimgs[flt2]
		imgkey3 = tplimg_flt2.split('.')[0]
	else:
		raise IOError('at least img3 or flt2 is not None...')

	
	photmethod = SN3.photometry_method
	if photmethod == 'psfphot':
		photretfile = imgkey3+'.psfphot'
		photret = os.path.join(SN3.psf_photometry_dir, photretfile)
	else:
		photretfile = imgkey3+'.apphot'
		photret = os.path.join(SN3.aperture_photometry_dir, photretfile)
		
	xymags3 = np.loadtxt(photret)

	#trim the match source list if needed
	xs_S3 = xymags3[:,0]
	ys_S3 = xymags3[:,1]


	if trim_img3 is not None:
		xl = trim_img3[0]
		xu = trim_img3[1]
		yl = trim_img3[2]
		yu = trim_img3[3]

		x_trim_mask = np.logical_and(xs_S3>xl, xs_S3<xu)	
		y_trim_mask = np.logical_and(ys_S3>yl, ys_S3<yu)
		trim_mask   = np.logical_and(x_trim_mask, y_trim_mask)
		
		xymags3 = xymags3[trim_mask,:]

	#prepare magnitudes
	 
	xymags3[:,2] = xymags3[:,2] + mag_offset_colorflt


	print xymags1
	print xymags2
	print xymags3


	plt.plot(xymags1[:,0],xymags1[:,1],'o')	
	plt.plot(xymags3[:,0],xymags3[:,1],'ro')
	plt.show()

	###=======================================================================###
	#match stars in two images from the same system
	if color_mode == 'color_S1':	
		reflist_filename = snname1+"_"+S1+"_"+imgkey1+'.txt'
		inputlist_filename = snname1+"_"+S1+"_"+imgkey3+'.txt'
	else:
		reflist_filename = snname2+"_"+S2+"_"+imgkey2+'.txt'
		inputlist_filename = snname2+"_"+S2+"_"+imgkey3+'.txt'
	
	
	np.savetxt(reflist_filename, xymags1, fmt="%8.3f %8.3f %8.3f %8.3f")
	np.savetxt(inputlist_filename, xymags3, fmt="%8.3f %8.3f %8.3f %8.3f")


	match_output = S1+snname1+imgkey1 + "_" + S1+snname1+imgkey3+'.match'
	match_coeff  = S1+snname1+imgkey1 + "_" + S1+snname1+imgkey3+'.coeff'

	SN1._photometry__fitsh_grmatch(reflist_filename, inputlist_filename, match_output, match_coeff)
	















	#-----------------------------------------------------------------------
	#match stars in two systems
