#! /usr/bin/env python

from astropy.io import fits
import os
import numpy.ma as ma
import numpy as np
import shutil
from numpy.lib import recfunctions
from numpy.core.records import recarray

from common import *

def setting_dir(master_dir):
	'''
	This sets up the directories to put all kinds of stuff
	
	Inputs:
		master_dir
	'''
	dirs = {}
	dirs['data_dir']       = os.path.join(master_dir, 'data')
	dirs['warehouse_dir']  = os.path.join(master_dir, 'warehouse')
	dirs['result_dir']     = os.path.join(master_dir, 'workplace')
	dirs['para_dir']       = os.path.join(master_dir, 'para')
	for dirkey in dirs.keys():
	    mkdir_p(dirs[dirkey])
	dirs['modified_image_dir']  = os.path.join(dirs['warehouse_dir'],'modified_images')
	dirs['stars_dir']      = os.path.join(dirs['warehouse_dir'],'xys')
	dirs['template_dir']   = os.path.join(dirs['warehouse_dir'],'template')
	dirs['std_ref_dir']    = os.path.join(dirs['warehouse_dir'],'std_ref')
	dirs['subtraction_dir']= os.path.join(dirs['warehouse_dir'],'image_subtraction')
	dirs['apphot_dir']     = os.path.join(dirs['warehouse_dir'], 'aperture_photometry')
	dirs['psfphot_dir']    = os.path.join(dirs['warehouse_dir'], 'psf_photometry')
	dirs['isophote_dir']   = os.path.join(dirs['warehouse_dir'], 'isophote')
	dirs['ds9reg_dir']     = os.path.join(dirs['warehouse_dir'], 'regions')
	dirs['psfmodel_dir']   = os.path.join(dirs['warehouse_dir'], 'psfmodel')
	dirs['astrometry_dir']   = os.path.join(dirs['warehouse_dir'], 'astrometry')
	dirs['log_dir']        = os.path.join(dirs['result_dir'] + 'log')	
	for dirkey in dirs.keys():
	     mkdir_p(dirs[dirkey])
	dirs['dophot_fortran'] = os.path.join(dirs['psfphot_dir'], 'dophot_fortran')
	dirs['dophot_c']       = os.path.join(dirs['psfphot_dir'], 'dophot_C')
	dirs['daophot_iraf']   = os.path.join(dirs['psfphot_dir'], 'daophot_iraf')
	for dirkey in dirs.keys():
	     mkdir_p(dirs[dirkey])
	
	
	return dirs

def image_preprocess(data_dir,filenames,middle_product_dir,log_dir):
    '''
	Preprocess task on images
    '''
    #replace the saturated pixes with a valid value which determined by 'fix_bad_pixels'
    fix_bad_pixels(data_dir,filenames,middle_product_dir,log_dir,mode= 'vicinity',progress_report = False)


    #transform the image date from float32(-32) to int(16) and the processed images will be stored in 'middle_product_dir'
    flt2int(data_dir,filenames,middle_product_dir,log_dir,progress_report = False)

    return 0


def get_image_info(data_dir,info_init,middle_product_dir,starfile_dir,log_dir,level = 0):
    '''
    Get the info of the images like 'filter','bkg','fwhm','nstar'
    INPUTS:
        data_dir,middle_product_dir,starfile_dir,log_dir
        level:
    '''

    if level == 0:
	 #better initial format needed
	filenames = info_init['name']
#	print filenames
        keywords = ['EXPTIME','FILTER']
        fitsinfo = get_fits_info(data_dir,filenames,keywords)

	exptime = fitsinfo['EXPTIME']
	flts    = fitsinfo['FILTER']
	#print flts
        image_info = recfunctions.append_fields(info_init, 'exptime',exptime,'f4',usemask = False, asrecarray = True)
        image_info = recfunctions.append_fields(image_info,'flts',flts,'S2',usemask = False, asrecarray = True)
	#print image_info

	info_dump_file = log_dir + 'image_middleinfo_dump.dat'
	save_results(info_dump_file, image_info)
	level = level+1

    if level ==1:
        info_dump_file = log_dir + 'image_middleinfo_dump.dat'
        image_info = restore_results(info_dump_file)
        filenames = image_info['name']
        bkg,fwhm,nstar=get_fwhm_bkg(data_dir,filenames,middle_product_dir,log_dir)

        image_info = recfunctions.append_fields(image_info,'bkg',bkg,float,usemask = False, asrecarray = True)
        image_info = recfunctions.append_fields(image_info,'fwhm',fwhm,float,usemask = False, asrecarray = True)
        image_info = recfunctions.append_fields(image_info,'nstar',nstar,int,usemask = False, asrecarray = True)

        info_dump_file = log_dir + 'image_middleinfo_dump.dat'
        save_results(info_dump_file, image_info)
        level = level+1

    if level ==2:
        info_dump_file = log_dir + 'image_middleinfo_dump.dat'
        image_info = restore_results(info_dump_file)
	filenames = image_info['name']
        starnum = get_star_number_rough(filenames,middle_product_dir,starfile_dir,log_dir)
        image_info = recfunctions.append_fields(image_info,'starnum',starnum,int,usemask = False, asrecarray = True)
    if os.path.isfile(info_dump_file):
        os.remove(info_dump_file)

    return image_info

def get_fwhm_bkg(data_dir,filenames,middle_product_dir,log_dir):
    '''
    data_dir:  the directory containing the image files
    filenames: without the absolute path to the data, just the fits file names under the data_dir' directory
    middle_product_dir: store the middle product like images with BITPIX = 16 and bad pixes fixed
    '''

    N = len(filenames)

    logfile = log_dir + 'get_fwhm_bkg.log'
    check_and_handle_file(logfile)
#    os.system('touch %s' %(logfile))
    fid = open(logfile,'a')

    nstar = []
    fwhm  = []
    bkg   = []
    for i,fi in enumerate(filenames):
        image_name = middle_product_dir+'int16_bpf_'+fi
	print image_name

	#the task will not stop for individual anomaly, the anomaly will be marked by -99, -99.99 and 9999.99
	try:
	    data = __get_fwhm_bkg(image_name)
	    nstar.append(data[-1])
	    fwhm.append(data[-2])
	    bkg.append(data[-3])

	    fid.write("%s succeed\n" % fi)

	except Exception,e:
	    nstar.append(-99)
	    fwhm.append(-99.99)
	    bkg.append(9999.99)
	    print "fwhm error on ",fi
	    fid.write("%s fail\n" % fi)

#        progress_report(i,N,'fwhm')
    fid.close()
    nstar = np.array(nstar)
    fwhm  = np.array(fwhm)
    bkg  = np.array(bkg)

    return bkg,fwhm,nstar


def __get_fwhm_bkg(image_name):
	'''
	'''
	templog    = 'fwhm.log'
	check_and_handle_file(templog)

        command = 'fwhm %s >> %s' %(image_name,templog)

	os.system(command)
	log = open(templog).readline()
	print log

	data = log.strip().split()[-5:]
	os.remove(templog)

	return data


def flt2int(data_dir,filenames,middle_product_dir,log_dir,progress_report = False):
	'''
	convert the data type from float32 to int to meet the requirement of some functions in diapl
	'''
	N = len(filenames)
	logfile = log_dir + 'flt2int.log'
	check_and_handle_file(logfile)
	fid = open(logfile,'a')
	
	for i,fi in enumerate(filenames):
	
	    input_filename = middle_product_dir+'bpf_'+fi
	    print input_filename
	
	    output_image = middle_product_dir+'int16_bpf_'+fi
	    check_and_handle_file(output_image)
	
	    try:
	        __flt2int(input_filename,output_image)
	        fid.write("%s succeed\n" % fi)
	
	    except:
	        fid.write("%s fail\n" % fi)
	        continue
	
	    if os.path.isfile(input_filename):
	        os.remove(input_filename)
	
	    if progress_report:
	        progress_report(i,N,'flt2in')
	
	fid.close()
	return 0

def __flt2int(input_image,output_image,):
	'''
	convert the data type from float32 to int to meet the requirement of some functions in diapl
	'''
	command = 'flt2int %s %s' %(input_image,output_image)
	os.system(command)
	
	return 0

def fix_bad_pixels(data_dir,filenames,middle_product_dir,log_dir,mode,progress_report = False):
    '''
    Replace the saturated pixes with normal value which determined by the 'mode'
    acceptable inputs for mode are'vicinity' or 'bkg'
    '''
    N = len(filenames)
    logfile = log_dir + 'fix_bad_pixels.log'
    check_and_handle_file(logfile)
    fid = open(logfile,'a')

    for i,fi in enumerate(filenames):

        #image = middle_product_dir+'int16_'+fi
        image = data_dir+fi
        print image
        fits_new = middle_product_dir+'bpf_'+fi

	try:
	    __fix_bad_pixels(image,fits_new,mode)
	    fid.write("%s succeed\n" % fi)
	except:
	    fid.write("%s fail\n" % fi)
	    continue

        if progress_report:
            progress_report(i,N,'fix_bad_pixels')
        #check_and_handle_file(image)
    fid.close()
    return 0

def __fix_bad_pixels(image,fits_new,mode):
	'''
	    Replace the saturated pixes with normal value which determined by the 'mode'
	'''
        hdu = fits.open(image)
        header = hdu[0].header
        data = hdu[0].data
        hdu.close()
        bad_pixes = get_bad_pixes(data,60000)
        data_new = fill_bad_pixes(data,bad_pixes,mode=mode)

        hdu_new = fits.PrimaryHDU(data_new,header)
        hdulist = fits.HDUList([hdu_new])

        check_and_handle_file(fits_new)
        hdulist.writeto(fits_new)

	return 0

def get_star_number_rough(filenames,input_data_dir,output_dir,log_dir):
	'''
	Get the approximate number of stars in the observation images
	'''
	N = len(filenames)
	
	logfile = log_dir + 'get_star_number_rough.log'
	check_and_handle_file(logfile)
	fid = open(logfile,'a')
	
	starnum = []
	for i,fi in enumerate(filenames):
	    filename = input_data_dir+'int16_bpf_'+fi
	
	    starfile =fi.split('.')[0]+ '_star.txt'
	    starfile_save = output_dir + starfile
	    check_and_handle_file(starfile)
	
	    try:
	        starnum_this = __get_star_number_rough(filename,starfile_save)
	        starnum.append(len(open(starfile_save).readlines()))
	        fid.write("%s succeed\n" % fi)
	    except:
	        fid.write("%s fail\n" % fi)
	        starnum.append(9999)
	        continue
	
	    progress_report(i,N,'sfind')
	
	fid.close()
	
	return np.array(starnum)

def __get_star_number_rough(filename,starfile_save):
	'''
	'''

        command = "sfind sfind.par instrument.par %s %s" %(filename,starfile_save)

        os.system(command)

        starnum_this = len(open(starfile_save).readlines())

	return starnum_this


def get_bad_pixes(image,threshold):
    '''
    bad pixes have value higher than bad value threshold
    '''
    mask = np.zeros(image.shape)
    mask[np.where(image > threshold)]=1
    return mask

def get_surrounding(radius = 1):
    '''
    get the pixes which are inside the circle with center at [0,0] and radius = radius
    '''
    shifts = []

    for x in np.arange(-int(radius),int(radius)+1):
        for y in np.arange(-int(radius),int(radius)+1):
            if x**2+y**2 <= radius**2:
                shifts.append([x,y])

    return shifts



def get_background(image):
    '''
    Get the background intensity of the given image
    '''
    imgbkg = image.reshape(1,np.size(image))
    cri = np.mean(imgbkg) + np.std(imgbkg)
    imgbkg = imgbkg[imgbkg<cri]
    #print imgbkg.shape
    background,bkgfilter,std = sigma_clipping(imgbkg, sig=3, meanfunc=np.median)
    #show_histogram(background,nbins = 100)
    bkg = np.mean(background)

    return bkg,std


def fill_bad_pixes(image,bad_pixes,mode='bkg'):
    '''
    fill the bad pixes with the desired value which is determined by option 'mode'
    mode = 'bkg' or 'vicinity'
    when mode = 'bkg' the bad_pix will be filled by the value of the image background
    when mode = 'vicinity' the bad_pix will be filled by the average value ofthe surrounding pixes(time consuming!!)
    '''
    M,N = image.shape
    X,Y = np.where(bad_pixes ==1)
    if mode=='vicinity':
        shifts = get_surrounding(radius = 1)
        #print shifts
        temp = ma.zeros((len(shifts),M,N))

        for i,[x,y] in enumerate(shifts):
            xlo = 0+x
            ylo = 0+y
            xho = M+x
            yho = N+y
            img_shift = cut_or_extend_imagedata(image,xlo,ylo,xho,yho)
            badpix_shift = cut_or_extend_imagedata(bad_pixes,xlo,ylo,xho,yho)
            masked_image = ma.array(img_shift,mask = badpix_shift)
            #print type(masked_image)
            temp[i] = masked_image
        #print type(temp)
        mean_image=np.mean(temp,axis = 0)
        for x,y in zip(X,Y):
            #print image[x,y],mean_image[x,y]
            image[x,y]=mean_image[x,y]
            #print i
    if mode=='bkg':
        bkg,std = get_background(image)
        image[X,Y]=bkg

    return image


def cut_or_extend_imagedata(imagedata,xlo,ylo,xho,yho):
    '''
    cut or extend the input image data
    we can get the size of the image from the image coordinate in which the left bottom pix has (0,0)
    and the right top pix has (M-1,N-1) where (M,N) = imagedata.shape
    the image coordinates of the left bottom pix and right top pix in the output image are (xlo,ylo),(xho,yho) respectively
    (xlo,ylo),(xho,yho) are give in the input image frame
    NOTE: the record of data on CCD is not the same as the normal array. the x axis is  along the bottom to top direction and
    the y axis is along the left to right direction
    '''
    data = imagedata
    temp = np.zeros((xho-xlo,yho-ylo),dtype=data.dtype.name)
    xtemp1,ytemp1 = temp.shape

    x0,y0=tuple([0,0])
    x1,y1=data.shape

    if xlo<=0:
        xf0 = x0
        xt0 = x0-xlo
    else:
        xf0 = xlo
        xt0 = x0

    if ylo<=0:
        yf0 = y0
        yt0 = y0 - ylo
    else:
        yf0 = ylo
        yt0 = y0

    if xho>=x1:
        xf1 = x1
        xt1 = x1-xlo
    else:
        xf1 = xho
        xt1 = xtemp1

    if yho>=y1:
        yf1 = y1
        yt1 = y1-ylo
    else:
        yf1 = yho
        yt1 = ytemp1

    #print xtemp1,ytemp1
    #print xf0,yf0,xf1,yf1
    #print xt0,yt0,xt1,yt1

    temp[xt0:xt1,yt0:yt1]=data[xf0:xf1,yf0:yf1]
    new_data = temp
    return new_data
