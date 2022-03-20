#! /usr/bin/python

#Description on the general principles of the photometry pipeline
#First, to protect the repository data under '/home/asassn/ASASSN/repository', symbolic links are created in workplace directory for each target, the default directory is '/home/asassn/ASASSN/$telescope/$sn', where $telescope is the telescope or insrument name. The linked image names are simple indexes(000.fits,001.fits...) and are used as a tracer of each image.
#Second, the photometry information is recorded in a dictory 'photometry_info' whose contents are dictories with the indexed image names as keywords. For example, xxx.photometry_info['000.fits'] is a dictionary which contain the datails of image '000.fits' as listed below:

# 'realimg', 	the source image name under '/home/asassn/ASASSN/repository/$telescope/$sn'
# 'flt',	filter
# 'obstime',	observation time in Julian date
# 'exptime',	exposure time in seconds
# 'bitpix',	data format in fits file
# 'flt2int',	float to integer transformation already? 1:yes; 0:no
# 'fixbad',	bad CCD pixes fixed already? 1:yes; 0:no
# 'rmCR',	cosmic ray removed? 1:yes; 0:no
# 'fwhm',	full width at half maximum of image profile
# 'bkg',	background value
# 'nstar',	source number detected in the field
# 'template',	is this image used as template? 1:yes; 0:no
# 'x',		
# 'y',		(x,y) is the image coordinate of your supernova
# 'instmag',	instrument magnitude of your target
# 'instmagerr',	
# 'relmag',	relative magnitude to the template image
# 'relmagerr',
# 'calmag',	magnitude after calibration
# 'calmagerr',
# 'drop'	is this image abandoned? non zero value assigned for abandoned image and different value indicate different 
#		reason, which can be accessed in 'xxx.drop_status'

#First created on Apr 6, 2016 by Ping Chen


#Notice: The following command can be executed in a interactive way in ipython, which is prefered for single target photometry


from ASASSN_photometry import photometry
import numpy as np

#You can use the following target as an example ASASSN-16bw
#Please delete the file '/home/asassn/ASASSN/LCOGT/$sn/workplace/*_phometry_info.txt' for a clean start if the file exists



#=================================== initiate ===========================================

sn = 'ASASSN-15nx'	#give the name of the supernova which will be used for the workplace directory for this target
telescope = 'LCOGT'

#initiate photometry instance with supernova name and telescope/instrument name
#you can either use the pipeline directory, for which you should put all your images in the direcotory under '/home/asassn/ASASSN/repository/$telescope/$sn'

phot = photometry(sn,telescope)

#or you can specify your prefered workspace directory 
#unblock the following line and comment the above line if you want this option

#$ img_dir = 		#the directory where your images are stored
#$ workspace_dir = 	#the directory to store the middle products and final results
#$ phot = photometry(sn,telescope,pipe=False,result_dir=workspace_dir,data_repository_dir=img_dir)



#================================= registration =========================================

phot.get_fwhm_bkg()			#get fwhm and background level
phot.get_star_number_rough()		#source extraction




#=============================== template image ==========================================
#select the "best" image as a reference/template image, basing on bkg,fwhm and star number found in the image field
#with a given relative weight in 'phot.template_factors' and the default value is [0.7,0.3,-0.5]
#unblock the following line and give your values

#$ phot.template_factors = []

phot.get_template_image()

#Or you can assign the templates in an abitrary way and this works in interactive way

#$ phot._photometry__dict2table()
#$ print np.unique(phot.result_table_newsaved['flt'])	#get filters in which the images are obtained

#check images information, for example you have images in B band

#$ Binfo = phot._photometry__select_rows_from_table(phot.result_table_newsaved,'flt','B')
#$ print Binfo	#Binfo is an astropy table and you can access each column like Binfo['fwhm']

#You can check the image quality by eye 

#$ d = phot.display_image_with_ds9(image_key)	#image_key is the linked image name like '001.fits'
						#this will open ds9 and load image named after 'image_key' in 
						#'/home/asassn/ASASSN/$telescope/$sn/data'
#$ d.set('exit')

#Once you get your 'best' template image, you can assign the template image as below

#$ phot.assign_template_image(flt,image_key)	$ eg. phot.assign_template_image('V','013.fits')


#================================ get target coordinate on image ==============================
#get the target coordinate on the template image and obtain the coordinates for other images by source match

phot._photometry__find_template_imagekey()
print phot.templates

phot.renew_template_xys = True		#use this before you want to clean the record of the target coordinate on template 
					#image and  redo the getting xy position task


#you have three options to get the supernova image coordinate 
phot.get_target_xys(tpl_method='ds9')		#the template image will be loaded into ds9 and you can select the target position by clicking on the wanted place
#$ phot.get_target_xys(tpl_method='manual')	#prompt come up and ask for xy position
#$ phot.get_target_xys(tpl_method='astrometry')	#If you want to use this option, astrometry solution is required. And you can put the astrometry solution (image file with correct wcs information or pure wcs file, the file name is 'cal_'+template_imagekey, for example 'cal_013.fits')



phot.save_results()	#save photometry work progress which can be reloaded when you have to run photometry for the same target and you don't want to repeat the above tasks


#===================================== aperture photometry ========================================
#a default parameter configuration file for aperture photometry with iraf/apphot task was put under '/home/asassn/ASASSN/$telescope/$sn/para'
#check the file for default value and change them to proper value for your images
#you may want to change the following paramters: app, aperture size in pixels; skyin,inner radius for background; skyout, outer radius for background

phot.get_apphot_iraf_parameters()


phot.aperture_photometry_apphot_iraf()


#If you  have a supernova near galaxy center, not isolated from the host galaxy center, and you still want
#to do aperture photometry for this target with a big aperture size

#phot.apphot_iraf_options['app'] = 18
#phot.apphot_iraf_options['skyin'] = 20
#phot.apphot_iraf_options['skyout'] = 35

#for flt in phot.templates.keys():
#	phot._photometry__aperture_photometry_apphot_iraf_target(flt)


phot.photometry_method = 'apphot'
phot.save_results()


#!!!!! psf photometry is still under optimization and will not be documented here !!!!!!!!!!!!!!!!


#======================================== calibration ====================================================
#calibration philosophy: relative calibration to the template image; standard  calibration of template image relative to standard reference stars provided

#phot.cal_plot = True	#uncomment this line if you want to see the calibration plot for each image

phot.get_relative_mags()


phot.get_standards()	#get standard magnitudes from APASS which have B,V,g',i',r' magnitudes 

#you can check standard stars in the working image field as below

#$ phot.load_std_stars_region_to_ds9(flt)	#for example phot.load_std_stars_region_to_ds9('V'), this will load the region file under '/home/asassn/ASASSN/$telescope/$sn/warehouse/std_ref'; the line number start from 0 is indicated in the ds9 images and you can modify the 'std_ref_$flt.reg' file to use only good reference stars


#phot.standard_calibration(tpl_has_wcs=False,all_tpl_astrometry=True,)
phot.standard_calibration(tpl_has_wcs=True)	#LCOGT images have wcs information 

phot.save_results()

phot.save_lc_data()	#save the light curve date under '/home/asassn/ASASSN/$telescope/$sn/workplace'


#phot.lc_plot('relmag')
phot.lc_plot('calmag')			#plot the light curve


