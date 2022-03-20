#!/usr/bin/python

##pipeline for photometry by Ping Chen

##modified by Ping on 2016-06-19: make aperture and PSF photometry compossible
##philosophy: the default is aperture photometry; PSF photometry needs specific flag

##modified by Ping on 2016-07-04: remove things about triming images

##SN position (x,y) comes from fits image in ds9 and iraf style, which treat the first pixel as (1,1)
##Many data process (sfind, dophot etc) treat the first poiont as (0,0)

##modification on 2016-10-11: (for aperture photometry,) aperture size used to be fixed to a value for all images;
##now adjust aperture size as well as background according to the image fwhm: aperture size 2*fwhm, sky_in 3*fwhm, sky_out 5*fwhm

##modification on 2016-10-28: aperture size 1.5*fwhm
##modify matplotlib related part under crontab job restriction(GUI not supported)
##self.standards type changed from ndarray to astropy.table

##modification on 2018-02-26: add functions to integrate photometry for WFCAM images
##add new image information element airmass to the photometry record table

##todo: estimate image background and improve get_fwhm_fistar
##todo: use astropy.stats.sigma_clipped_stats

###!!NOTES: self.psfphot_method default value is 'dophot' which will put the xxx.psfphot in self.psfphot_dophot_fortran_dir by default

##from __future__ import print_function

#import Packages
import numpy as np
import os
import datetime
import shutil
from astropy.table import Table,Column, vstack, hstack
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d,interp2d, griddata
from collections import OrderedDict
import ConfigParser
from photutils import fit_2dgaussian

try:
    import pyraf 
    print "pyraf loaded"
    from iraf import stsdas
    print "stsdas loaded"
except:
    print "pyraf loading has issue..."


try:
	import pyds9
except:
	print "pyds9 import failure..."

try:
	import matplotlib.pylab as plt
except:
	print "no display mode"

from extern import cosmics
from utils.ccdproc_cosmicray_removal_lacosmic import remove_CR_ccdproc_cosmicray_lacosmic
from utils.common import __get_fits_info as get_fits_info
from utils.common import sigma_clipping,get_best, radec_format_transformation

from utils.prepare_work import setting_dir
from utils.get_image_info_difftelescope import get_obstime,get_exptime,get_filter,get_bitpix, get_airmass_telescope
from utils.timeout_handle import Command
from utils.catalog_search import query_local_APASS,query_VO_SCS,query_General_MAST, query_local_ATLAS_Refcat2, gaia_query_mast_Catalogs
from utils.catalog_search import panstarrs_query_Vizier,apass_query_Vizier,twomass_query_Vizier, SDSS_DR12_query_Vizier, UKIDSS_query_Vizier, gaia2_query_Vizier
from utils.photometry_collection import apphot_iraf

from utils.mag_flux_convertion import Vega_AB_mag_convertion, mag2flux, flux2mag
from utils.get_airmass_obs import get_airmass_given_time_target_observatory
from utils.Iowa_fits_header_simplify_standardize import simplify_and_standardize_fits_header

from mpltk.event_handler import mouse_pick
from mpltk.plottask import *

from iraftk.iraf_daophot_daofind import daofind_iraf
from iraftk.iraf_imstat import imstat_iraf
from iraftk.iraf_imarith import imarith_iraf
from iraftk.iraf_imcopy import imcopy_iraf
from iraftk.iraf_imcombine import imcombine_iraf
from iraftk.iraf_xytran import geomap_iraf, geoxytran_iraf
from iraftk.iraf_utils import *
from iraftk.iraf_daophot_phot import daophot_phot_iraf
from iraftk.iraf_daophot_psf import daophot_psf_iraf
from iraftk.iraf_daophot_seepsf import daophot_seepsf_iraf
from iraftk.iraf_daophot_peak import daophot_peak_iraf, daophot_substar_iraf
from iraftk.iraf_daophot_allstar import daophot_allstar_iraf
from iraftk.iraf_stsdas_isophote_ellipse import stsdas_analysis_isophote_ellipse,stsdas_analysis_isophote_bmodel#,stsdas_analysis_isophote_cmodel

from pyds9tk.select_or_delete_from_given_list_ds9_display import get_xy_on_image_from_ds9, input_xys_through_ds9_get_mags
from pyds9tk.select_source_on_ds9 import get_lines_by_selecting_xys_on_image_from_ds9
from pyds9tk.ds9_display_sources import create_ds9_region_file, pylab_talk_to_ds9


import warnings
warnings.filterwarnings('ignore')

from photbasesetup import basesetup

import pickle

def load_phot_target(picklefile):
	'''
	Load pickled photometry object
	'''
	try:
		sn = pickle.load(open(picklefile,'r'))
	except:
		sn = None
		print "failed to load picklefile %s"%picklefile
	return sn

class photometry():
	'''
	photometry class
	'''
# Initial set up
	def __init__(self,sn,telescope,pipe=True, auto_init = True, result_dir=None,data_repository_dir=None,photometry_method = 'apphot', image_subtraction=False):
		'''
		inputs:
			sn: the target name
			telescope: the telescope name
			pipe: whether in pipeline mode
			auto_init:
			result_dir: the master workplace folder for the given target
			data_repository_dir: the input data repository directory
			photometry_method:
		'''
		self.base = basesetup()
		self.pipe = pipe
		self.auto_init = auto_init
		self.current_sn = sn
		self.current_telescope = telescope
		self.photometry_method = photometry_method
		self.image_subtraction = image_subtraction
		self.readme = "" #Please keep the notes

		self.ds9D = None #store the current instance of ds9
		#set up data repository and photometry directory system
		if not pipe:
			self.repository_dir_current_sn =  os.path.abspath(data_repository_dir)
			self.current_sn_dir = os.path.abspath(result_dir)
		else:
			self.repository_dir_current_telescope = os.path.join(self.base.data_repository_master_dir,self.current_telescope)
			self.repository_dir_current_sn = os.path.join(self.repository_dir_current_telescope,self.current_sn)
			self.current_sn_master_dir = os.path.join(self.base.workplace_master_dir,self.current_telescope)
			self.current_sn_dir = os.path.join(self.current_sn_master_dir, self.current_sn)
		if not os.path.exists(self.repository_dir_current_sn):
			raise IOError("directory %s not found"%self.repository_dir_current_sn)
		if not os.path.exists(self.current_sn_dir):
			os.mkdir(self.current_sn_dir)

		self.__init_settings__()
		self.__init_dophot_pm()
		self.dophot_select_stars_method = 'lt' #
		self.dophot_select_star_value = 8
		self.__init_renew_trigger()
		self.__init_fitsh_grmatch_options()
		self.__init_fitsh_fitrans_options()
		self.__init_hotpants_options()
		self.hotpants_par_hold_list = [] #parameters from this list will directly use those from self.hotpants_par without further modification
		self.__init_std_catalog_config()


		self.imagepreprocess_verbose = False
		self.cal_plot = False
		self.sigclip_first_preoffset = 1
		self.cal_offset_method='funcfit'#'median','mean' or 'funcfit';the method adopted to get the offset between instrumental mags and std mags
		self.cal_offset_funcfit_type = 'constant' #'constant', 'o1poly', 'o2poly'
		self.cal_offset_nsig_clip = 5.0
		self.cal_offset_verbose = False
		self.cal_offset_err_method = 'brightness_dependent_fit' #'gobal_std'
		self.cal_offset_err_funcfit_type = 'o2poly'

		self.cal_check_stds = False

		if self.photometry_method == 'apphot':
			methodcode = 'AP'
		elif self.photometry_method == 'psfphot':
			methodcode = 'PSF'
		else:
			raise IOError("Invalid input for photometry_method...")

		if self.image_subtraction:
			subflag = '_sub'
		else:
			subflag =  ''
		photfile_initversion =  self.current_sn+'_photometry_info_'+ methodcode + subflag +'.txt'

		self.result_table_file_init =  os.path.join(self.result_dir, photfile_initversion)
		photfile_modversions = [f for f in os.listdir(self.result_dir) if f.split('.')[-1].isdigit() and photfile_initversion=='.'.join(f.split('.')[:-1])]
		output_mvs =  np.sort([int(f.split('.')[-1]) for f in photfile_modversions]) #modified output versions
		print photfile_modversions
		print output_mvs
		if len(output_mvs)>0:
			self.currentversion = output_mvs[-1]
			self.result_table_file = self.result_table_file_init + '.' + str(output_mvs[-1])
		else:
			self.currentversion = 0
			self.result_table_file = self.result_table_file_init

		self.notesfile = os.path.join(self.result_dir, 'README.txt')

		self.new = True
		self.result_table = None
		self.photometry_record_table = None

		if os.path.isfile(self.result_table_file):
			print "Archived results from %s have been loaded"%self.result_table_file
			self.result_table =Table.read(self.result_table_file,format='ascii.fixed_width')
			self.new = False


		if not self.new:
			self.reln_image_to_datadir()

		self.ln_image_to_datadir(self.repository_dir_current_sn,self.raw_image_dir)

		self.images = self.get_images_dict(self.raw_image_dir)
		self.apphots  = {}
		self.psfphots = {}
		for img in self.images.keys():
			img_apphot  = os.path.join( self.aperture_photometry_dir, img.split('.')[0]+self.apphot_ret_file_suffix)
			psf_photometry_dir = self.__get_internal_psfphot_dir()
			img_psfphot = os.path.join(psf_photometry_dir, img.split('.')[0]+self.psfphot_ret_file_suffix)
			self.apphots[img] = img_apphot
			self.psfphots[img]= img_psfphot

		self.stdobs_match_files	= {}
		self.stdobs_transcoef_files = {}

		self.__init_photometry_info()
		self.__init_iraf_apphot_config()
		self.__init_fill_images_info_and_prepare_parfiles()
		self.__init_iraf_isophote_ellipse_pm()
		self.__init_stdref_filenames()
		self.psfrad_nfwhm = 4.0
		self.fitrad_nfwhm = 1.3 #around 3 gaussian sigma

	def __init_fill_images_info_and_prepare_parfiles(self):
		'''
		fill out the image information dictionary: realimg, filter name, obstime, exposure time, bitpix
		'''
		if self.new:
			self.__prepare_parameter_files()
			for key in self.photometry_info.keys():
				realimg_abs = os.path.realpath(self.images[key])
				self.photometry_info[key]['realimg'] = realimg_abs.split('/')[-1]
			if self.auto_init:
				for key in self.photometry_info.keys():
					self.__get_image_size_from_fitsheader(key)
				for key in self.photometry_info.keys():
					self.__renew_flt_info(key)
				for key in self.photometry_info.keys():
					self.__renew_obstime_info(key)
				for key in self.photometry_info.keys():
					self.__renew_exptime_info(key)
				for key in self.photometry_info.keys():
					self.__renew_bitpix_info(key)
		else:
			self.__load_old_results(self.result_table,self.photometry_info)
			for key in self.photometry_info.keys():
				if key not in self.result_table['name']:
					realimg_abs = os.path.realpath(self.images[key])
					self.photometry_info[key]['realimg'] = realimg_abs.split('/')[-1]
					if self.auto_init:
						self.__get_image_size_from_fitsheader(key)
						self.__renew_flt_info(key)
						self.__renew_obstime_info(key)
						self.__renew_exptime_info(key)
						self.__renew_bitpix_info(key)

	def __init_settings__(self):

		dirdict = setting_dir(self.current_sn_dir)
		self.raw_image_dir     = dirdict['data_dir']
		self.warehouse_dir     = dirdict['warehouse_dir']
		self.modified_image_dir= dirdict['modified_image_dir']
		self.stars_dir         = dirdict['stars_dir']
		self.template_dir      = dirdict['template_dir']
		self.result_dir        = dirdict['result_dir']
		self.log_dir           = dirdict['log_dir']
		self.std_ref_dir       = dirdict['std_ref_dir']
		self.subtraction_dir   = dirdict['subtraction_dir']
		self.parafile_dir      = dirdict['para_dir']
		self.psf_photometry_dir= dirdict['psfphot_dir']
		self.ds9_region_dir    = dirdict['ds9reg_dir']
		self.psfmodel_dir      = dirdict['psfmodel_dir']
		self.aperture_photometry_dir = dirdict['apphot_dir']
		self.isophote_dir      = dirdict['isophote_dir']
		self.psfphot_dophot_fortran_dir = dirdict['dophot_fortran']
		self.psfphot_dophot_c_dir       = dirdict['dophot_c']
		self.psfphot_daophot_iraf_dir   = dirdict['daophot_iraf']
		self.astrometry_dir = dirdict['astrometry_dir']

		self.parafile_template_dir = os.path.join(self.base.conffile_dir, 'para_templates')

		self.starfile_suffix = '.star'                        #This is the simplied version of source list with the format:
		self.iraf_daofind_output_suffix = '_iraf_daofind.out' #This is the output filename suffix for iraf/daofind, for example 000_iraf_daofind.out
		self.regfile_suffix = '.reg'
		self.apphot_ret_file_suffix = '.apphot'   	   #This is apphot photometry output modified/simplified version: x,y,mag,magerr
		self.phot_iraf_apphot_suffix = '_iraf_apphot.out'  #This is original apphot photometry output, refer to iraf for the detailed format
		self.apphot_iraf_parafile = 'apphot_iraf_param.par'
		self.psfphot_ret_file_suffix = '.psfphot'

		self.templates = {}
		self.templates_after_astrometry = {}
		self.template_factors = [0.7,0.3,-0.5]		#the weight in choosing the best image

		try:
			if self.current_sn in self.base.first_targets_info['Name']:
				RA_str  = self.base.first_targets_info[self.base.first_targets_info['Name']==self.current_sn]['RA'].data[0]
				Dec_str = self.base.first_targets_info[self.base.first_targets_info['Name']==self.current_sn]['Dec'].data[0]
				RA_deg,Dec_deg = radec_format_transformation(RA_str,Dec_str)
				print "%s (RA,DEC)=(%s,%s): (%s,%s)"%(self.current_sn, RA_str, Dec_str, RA_deg, Dec_deg)
			elif self.current_sn in self.base.second_targets_info['Name']:
				RA_str  = self.base.second_targets_info[self.base.second_targets_info['Name']==self.current_sn]['RA'].data[0]
				Dec_str = self.base.second_targets_info[self.base.second_targets_info['Name']==self.current_sn]['Dec'].data[0]
				RA_deg,Dec_deg = radec_format_transformation(RA_str,Dec_str)
				print "%s (RA,DEC)=(%s,%s): (%s,%s)"%(self.current_sn, RA_str, Dec_str, RA_deg, Dec_deg)
			elif self.current_sn in self.base.registerd_photometry_objects_info:
				RA_deg = float(self.base.registerd_photometry_objects_info[self.current_sn]['ra'])
				Dec_deg = float(self.base.registerd_photometry_objects_info[self.current_sn]['dec'])
				RA_str,Dec_str = radec_format_transformation(RA_deg,Dec_deg, mode='float2str')
				print "%s (RA,DEC)=(%s,%s): (%s,%s)"%(self.current_sn, RA_str, Dec_str, RA_deg, Dec_deg)
			else:
				RA_deg = None
				Dec_deg = None
				RA_str = None
				Dec_str = None

			self.sn_ra_world_deg = RA_deg
			self.sn_dec_world_deg = Dec_deg
			self.sn_ra_str = RA_str
			self.sn_dec_str = Dec_str
		except Exception as e:
			print e
			self.sn_ra_world_deg = None
			self.sn_dec_world_deg = None
			self.sn_ra_str = None
			self.sn_dec_str = None

		self.current_ccd_readnoise = None
		self.current_ccd_epadu = None

		self.pixscale = None
		self.std_region_radius = 0.33	#degree

		self.physical_image_offset_x = 0	#difference between images coordinate and physical coordinate
		self.physical_image_offset_y = 0

		self.stdcal_xytran_stdnum = 3	#how many stars to select when xytran method used to link std and obs

		self.trim_before_grmatch = False
		self.trim_left = None		# For example, self.trim_left = n1, self.trim_right=n2, self.trim_bottom = n3, self.trim_up =n4
		self.trim_right = None		# the sources within [n1:n2] in x direction and [n3:n4] in y direction survive
		self.trim_bottom = None 	# self.trim_right and self.trim_up can be negative value
		self.trim_up = None


		self.filter_after_fistar = False #select identidied sources within specfic region
		self.filter_after_daofind = True #select identidied sources within specfic region
		self.edge_xl = 0	# source detection image coordinate range
		self.edge_xu = 100000
		self.edge_yl = 0
		self.edge_yu = 100000

		self.magnitude_cut_before_grmatch = False
		self.reference_mag_min = None	#the faint end of reference magnitudes for calibration
		self.reference_mag_max = None	#the bright end of reference magnitudes for calibration
		self.reference_mag_min_note = 'the magnitude cut in the faint end'
		self.reference_mag_max_note = 'the magnitude cut in the bright end'

		self.reference_magerr_cut = None	#the threshold for filtering magnitudes with large uncertainties
		self.input_magerr_cut = None

		self.criteria_stdref_tpl_match = 10	# stdref and tpl match in 'surrouding_search' method
		self.criteria_tpl_obs_match = 10
		self.criteria_point_match_general = 5   # general souch match criteria 
		self.criteria_psfphot_locate_target = 2	# find corresponding target in psfphot results for SN
		#self.apphot_xy_centroid	= True		#

		self.maxtime_on_single_fwhm = 60        #maximum time allowed for obtaining fwhm from single image
		self.maxtime_on_astrometry = 300	#maxmum time for obtaining astrometry solution

		self.tpl_obs_before_match_display = False
		self.tpl_obs_match_result_display     = False	#display the match result of stars on template and input image
		self.stdref_tpl_match_result_display  = False   #display the match result of stars on template and standard calibration database


		#if os.path.isfile(os.path.join(self.std_ref_dir,self.std_ref_stars_file)):
		#	self.standards = np.loadtxt(os.path.join(self.std_ref_dir,self.std_ref_stars_file))

		self.drop_status = {}
		self.drop_status['0'] = 'safe'
		self.drop_status['1'] = 'bad image quality: background/bad flat'
		self.drop_status['2'] = 'bad image quality: sparse stars/low transparency/star finding error'
		self.drop_status['3'] = 'fwhm'
		self.drop_status['4'] = 'target not found on image'
		self.drop_status['5'] = 'grmatch failure'
		self.drop_status['6'] = 'standard calibration failure'
		self.drop_status['7'] = 'reference match failure'
		self.drop_status['8'] = 'photometry result anomaly/photometry result records less than 10'
		self.drop_status['9'] = 'more than one sources detected at given position'
		self.drop_status['10'] = 'target affccted by cosmic ray or other artifacts'
		self.drop_status['11'] = 'uncharted new phenomenon'
		self.drop_status['12'] = 'get filter information error'
		self.drop_status['13'] = 'image subtraction failure'

		self.psfphot_method = 'dophot' #alternative: daophot
		self.dophot_version = 'fortran' #the version of dophot used for PSF photometry: 'fortran' or 'C'
		self.dophot_verbose = False
		self.dophot_image_prepare = 'imcopy' #or softlink

		#light curves
		self.lcs_raw = {}

		self.host_mags_file = 'host_mags.txt'
		self.host_mags_dict = {}
		self.lcs_hostfree = {}

		#insulate the pmfile to prevent from deleting by another process when multiple psf photometry running
		self.insulate_pmfile = False

	def __init_std_catalog_config(self):
		'''
		standard catalog configuration
		'''
		#local APASS data directories
		self.local_apass_dir = os.path.join(self.base.stdcat_dir,'APASS')
		self.local_atlas_refcat2_dir = os.path.join(self.base.stdcat_dir, 'ATLAS_Refcat2')

		self.std_ref_stars_file = {}
		self.std_ref_stars_file['apass'] = 'standard_reference_star_apass.txt'
		self.std_ref_stars_file['2mass'] = 'standard_reference_star_2mass.txt'
		self.std_ref_stars_file['panstarrs'] = 'standard_reference_star_panstarrs.txt'
		self.std_ref_stars_file['refcat2'] = 'standard_reference_star_atlas_refcat2.txt'
		self.std_ref_stars_file['sdss'] = 'standard_reference_star_sdss.txt'
		self.std_ref_stars_file['ukidss'] = 'standard_reference_star_ukidss.txt'
		self.std_ref_stars_file['gaia'] = 'standard_reference_star_gaia.txt'

		self.flt_std_ref_stars_file = {}
		self.flt_std_ref_stars_file_big = {}
		self.flt_std_ref_and_obs_measurement_stars_file = {}

		self.standards = None	#all information for standard stars

		self.panstarrs_stds_gp = None #only (x,y, mag_flt, magerr_flt)
		self.panstarrs_stds_rp = None
		self.panstarrs_stds_ip = None
		self.panstarrs_stds_zp = None
		self.panstarrs_stds_yp = None

		self.apass_stds_B = None
		self.apass_stds_V = None
		self.apass_stds_rp= None
		self.apass_stds_ip= None

		self.twomass_stds_J  = None
		self.twomass_stds_H  = None
		self.twomass_stds_K  = None

		self.reference_flt = None	# if template images don't have wcs info and we don't want API astrometry for all template then only one image with this 'reference_flt' will be given wcs and this wcs info will be used for other templates

		self.Vizier_catalog_all_columns = False

		self.gaia_catalog_method = 1
		self.gaia_colnames_dict = {'RA':'RA_ICRS', 'Dec':'DE_ICRS','G':'Gmag', 'Gerr':'e_Gmag', 'BP':'BPmag', 'BPerr':'e_BPmag', 'RP':'RPmag', 'RPerr':'e_RPmag'}

		self.ukidss_catalog_method = 1
		self.ukidss_colnames_dict = {'RA':'RAJ2000', 'Dec':'DEJ2000','J':'Jmag', 'Jerr':'e_Jmag', 'H':'Hmag', 'Herr':'e_Hmag', 'K':'Kmag', 'Kerr':'e_Kmag'}

		self.sdss_catalog_method = 1
		self.sdss_colnames_dict =  {'RA':'RA_ICRS', 'Dec':'DE_ICRS','u':'umag', 'g':'gmag', 'r':'rmag', 'i':'imag', 'z':'zmag', 'gerr':'e_gmag', 'rerr':'e_rmag', 'ierr':'e_imag', 'zerr':'e_zmag'}
		for kw in self.sdss_colnames_dict.keys():
			if kw in ['RA', 'Dec']:
				continue
			if len(kw)==1:
				kwnew = kw+'p'
			else:
				kwnew = kw.replace('err', 'perr')
			self.sdss_colnames_dict[kwnew] = self.sdss_colnames_dict[kw]

		self.atlas_refcat2_method = 1
		self.refcat2_colnames_dict = {'RA':'RAJ2000', 'Dec':'DEJ2000','B':'Bmag', 'V':'Vmag', 'R':'Rmag', 'I':'Imag', 'gp':'gmag', 'rp':'rmag', 'ip':'imag', 'zp':'zmag','Berr':'e_Bmag', 'Verr':'e_Vmag', 'Rerr':'e_Rmag', 'Ierr':'e_Imag', 'gperr':'e_gmag', 'rperr':'e_rmag', 'iperr':'e_imag', 'zperr':'e_zmag', 'J':'Jmag', 'Jerr':'e_Jmag', 'H':'Hmag', 'Herr':'e_Hmag', 'K':'Kmag', 'Kerr':'e_Kmag'}

		self.apass_catalog_method = 2
		self.apass_catalog_method_notes = "multiple methods to retrive APASS catalog provided: [1. local data (the whole APASS distribution) extraction; 2, apass_query_Vizier]"
		self.apass_colnames_dict = {'RA':'RAJ2000', 'Dec':'DEJ2000', 'B':'Bmag', 'V':'Vmag', 'g':'g_mag', 'r':'r_mag', 'i':'i_mag', 'Berr':'e_Bmag', 'Verr':'e_Vmag', 'gerr':'e_g_mag', 'rerr':'e_r_mag', 'ierr':'e_i_mag'}

		#The followings are related with stds from PanSTARRS catalog
		self.panstarrs_catalog_method = 2
		self.panstarrs_catalog_method_notes = "multiple methods to retrive PS1 catalog provided: [1, query_General_MAST (http://gsss.stsci.edu); 2, panstarrs_query_Vizier (based on astroquery.Vizier)]"
		self.panstarrs_colnames_dict = {'RA':'RAJ2000', 'Dec':'DEJ2000', 'g':'gmag', 'r':'rmag', 'i':'imag', 'z':'zmag', 'y':'ymag', 'gerr':'e_gmag', 'rerr':'e_rmag', 'ierr':'e_imag', 'zerr':'e_zmag',  'yerr':'e_ymag', 'gp':'gmag', 'rp':'rmag', 'ip':'imag', 'zp':'zmag', 'yp':'ymag', 'gperr':'e_gmag', 'rperr':'e_rmag', 'iperr':'e_imag', 'zperr':'e_zmag',  'yperr':'e_ymag', 'B': 'Bmag', 'V': 'Vmag','R': 'Rmag','I': 'Imag', 'Berr': 'e_Bmag','Verr': 'e_Vmag','Rerr': 'e_Rmag','Ierr': 'e_Imag'}

		self.panstarrs_colnames_dict_method1 = {'RA':'raMean', 'Dec':'decMean','gp':'gMeanPSFMag', 'gperr':'gMeanPSFMagErr', 'rp':'rMeanPSFMag', 'rperr':'rMeanPSFMagErr', 'ip':'iMeanPSFMag','iperr':'iMeanPSFMagErr','zp':'zMeanPSFMag', 'zperr':'zMeanPSFMagErr', 'yp':'yMeanPSFMag', 'yperr':'yMeanPSFMagErr', 'R':'Rmag', 'Rerr':'Rmagerr', 'I':'Imag', 'Ierr':'Imagerr'}
		self.panstarrs_colnames_dict_method1_AP = {'RA':'raMean','Dec':'decMean', 'gp_AP':'gMeanApMag', 'gperr_AP':'gMeanApMagErr', 'rp_AP':'rMeanApMag', 'rperr_AP':'rMeanApMagErr', 'ip_AP':'iMeanApMag', 'iperr_AP':'iMeanApMagErr', 'zp_AP':'zMeanApMag', 'zperr_AP':'zMeanApMagErr', 'yp_AP':'yMeanApMag', 'yperr_AP':'yMeanApMagErr', 'R_AP':'Rmag', 'Rerr_AP':'Rmagerr', 'I_AP':'Imag', 'Ierr_AP':'Imagerr'}


		self.twomass_catalog_method = 2
		self.twomass_catalog_method_notes = "multiple methods to retrive 2MASS catalog provided: [1, query_VO_SCS (https://irsa.ipac.caltech.edu); 2, twomass_query_Vizier (based on astroquery.Vizier)]"
		self.twomass_colnames_dict_method1 = {'RA':'ra', 'Dec': 'dec','J':'j_m', 'Jerr':'j_msigcom', 'H':'h_m', 'Herr':'h_msigcom', 'K':'k_m', 'Kerr':'k_msigcom'}
		self.twomass_colnames_dict = {'RA':'RAJ2000', 'Dec':'DEJ2000', 'J':'Jmag', 'Jerr':'e_Jmag', 'H':'Hmag', 'Herr':'e_Hmag', 'K':'Kmag', 'Kerr':'e_Kmag'}

		self.panstarrs_mag_photometry_method = 'PSF' # PanSTARRS catalog provides magnitudes from two photometry kinds: PSF or AP, only works for method 1
		self.panstarrs_min_nx = None		     # colname 'ng', 'nr', 'ni', 'nz', 'ny'
		self.panstarrs_min_nstackDetections = None

		self.gaia_colnames_dict = {'RA':'RA_ICRS', 'Dec':'DE_ICRS', 'G':'Gmag', 'Gerr':'e_Gmag'}

		self.std_catalog_colnames = {'apass': self.apass_colnames_dict, 'panstarrs': self.panstarrs_colnames_dict, '2mass':self.twomass_colnames_dict, 'refcat2': self.refcat2_colnames_dict, 'gaia': self.gaia_colnames_dict}

		self.selfphot_mag_faint_cut = None	#the magnitudes from private photometry results used as for calibration
		self.selfphot_mag_saturate_cut = None
		self.selfphot_magerr_up_cut = None #upper limit for magnitude uncertainty
		self.selfphot_match_method_for_RI_stds = 'surrounding_search' #if cal std from selfphot, when prepare std for R,I band then r,i std are needed, then matching

		self.panstarrs_mag_faint_cut = None
		self.panstarrs_mag_saturate_cut = None

		self.sdss_mag_faint_cut = None
		self.sdss_mag_saturate_cut = None

		self.refcat2_mag_faint_cut = None
		self.refcat2_mag_saturate_cut = None

		self.apass_mag_faint_cut = None
		self.apass_mag_saturate_cut = None
		self.apass_nobs_min = None
		self.apass_mobs_min = None
		self.apass_remove_single_measurement = True

		self.twomass_mag_faint_cut = None
		self.twomass_mag_saturate_cut = None

	def __init_stdref_filenames(self):
		'''
		prepare some filename for standard catalog
		'''
		#self.__find_template_imagekey()
		flts = np.unique([self.photometry_info[imgkey]['flt'] for imgkey in self.images.keys()])
		for flt in flts:
			self.flt_std_ref_stars_file[flt] = 'std_ref_' +flt+ '.txt'
			self.flt_std_ref_stars_file_big[flt] = 'std_ref_whole_info_' +flt+ '.txt'
			if flt in self.templates.keys():
				tpl_imgkey = self.templates[flt]
				tpl_imgkey_s = tpl_imgkey.split('.')[0]
				self.flt_std_ref_and_obs_measurement_stars_file[flt] = tpl_imgkey_s +'_std_'+flt+'.match'

	def __init_fitsh_grmatch_options(self):
		'''
		See docs/fitsh_grmatch.help for details on the options
		'''
		self.fitsh_grmatch_type = 'point'  #'point'--point matching, 'coord'--coordinate matching, or 'id'--identifier matching

		self.fitsh_grmatch_pointmatch_pars = {}
		self.fitsh_grmatch_pointmatch_pars['--col-ref'] = '1,2' #The index of the first column is always 1
		self.fitsh_grmatch_pointmatch_pars['--col-inp'] = '1,2'
		self.fitsh_grmatch_pointmatch_pars['--order'] = 1 #If the order is A, >= (A+1)*(A+2)/2 valid points are needed to fit the transformation
		self.fitsh_grmatch_pointmatch_pars['--max-distance'] = 1
		self.fitsh_grmatch_pointmatch_pars['--triangulation'] = 'auto,mixed,maxnumber=200'
		self.fitsh_grmatch_pointmatch_pars['--col-ref-ordering'] = -3 #negative sign indicates ascending, small values first
		self.fitsh_grmatch_pointmatch_pars['--col-inp-ordering'] = -3
		self.fitsh_grmatch_pointmatch_pars['--fit'] = None  #iterations=<N>,firstrejection=<F>,sigma=<S>
		self.fitsh_grmatch_pointmatch_pars['--weight'] = None

		self.fitsh_grmatch_coordmatch_pars = {}
		self.fitsh_grmatch_coordmatch_pars['--col-ref'] = '1,2' #The index of the first column is always 1
		self.fitsh_grmatch_coordmatch_pars['--col-inp'] = '1,2'
		self.fitsh_grmatch_coordmatch_pars['--max-distance'] = 1

		self.fitsh_grmatch_idmatch_pars = {}

	def __init_fitsh_fitrans_options(self):
		'''
		See docs/fitsh_fitrans.help for details
		'''
		self.fitsh_fitrans_pars = {}
		self.fitsh_fitrans_pars['-m'] = False
		self.fitsh_fitrans_pars['-l'] = False
		self.fitsh_fitrans_pars['-c'] = False
		self.fitsh_fitrans_pars['-k'] = True
		self.fitsh_fitrans_pars['--reverse'] = True #--reverse  tranform the input image to the reference image

	def __init_renew_trigger(self):
		'''
		controller on whether renew some properties
		'''
		self.renew_template_xys = False			#
		self.renew_target_xys = False
		self.renew_aperture_photometry = False		#renew the instmag of self.photometry_info['xxx']['instmag']
		self.renew_aperture_photometry_retfile = False	#renew the aperture photometry result file
		self.renew_psf_photometry = False
		self.renew_psf_photometry_retfile = False
		self.renew_relative_mag = False
		self.renew_stdcal_mag = False
		self.renew_std_ref_match = False
		self.renew_stdcal_xytran = False
		self.renew_standards = False

		self.renew_apphot_parfile = True #overwrite existing iraf-apphot parameter file?

	def __init_iraf_apphot_config(self):
		'''
		parameters init for iraf.apphot
		'''
		self.apphot_iraf_options = {}
		self.apphot_iraf_calgorithm = 'centroid' #the default centering algorithm for apphot, other options 'gauss', 'ofilter'
		self.apphot_iraf_datamin = None
		self.apphot_iraf_datamax = None

		self.apphot_iraf_options['fwhmpsf'] = None
		self.apphot_iraf_options['app'] = 8
		self.apphot_iraf_options['skyin'] = 10
		self.apphot_iraf_options['skywidth'] = 20
		self.apphot_iraf_options['sky_sigma_down'] = 3
		self.apphot_iraf_options['sky_sigma_up'] = 3
		self.apphot_iraf_options['def_zeropt'] = 25

		self.apphot_iraf_autoscale_aperture = True
		self.apphot_iraf_app_nfwhm = 2
		self.apphot_iraf_skyin_nfwhm = 3
		self.apphot_iraf_skywidth_nfwhm = 3

	def __init_iraf_isophote_ellipse_pm(self):
		'''
		parameters for iraf-ststas-analysis-isophote-ellipse
		'''
		self.ellipse_pm = {}
		self.ellipse_pm['ellip0']=0.2
		self.ellipse_pm['pa0']=45
		self.ellipse_pm['sma0']=10
		self.ellipse_pm['step']=0.1
		self.ellipse_pm['linear']=False
		self.ellipse_pm['recenter']=False
		self.ellipse_pm['xylearn']=False
		self.ellipse_pm['physical']=True
		self.ellipse_pm['conver']=0.05
		self.ellipse_pm['minit']=10
		self.ellipse_pm['maxit']=50
		self.ellipse_pm['hcenter']=False
		self.ellipse_pm['hellip']=False
		self.ellipse_pm['hpa']=False
		self.ellipse_pm['wander'] = 'INDEF'
		self.ellipse_pm['maxgerr']=1.0
		self.ellipse_pm['olthresh']=0.0
		self.ellipse_pm['integrmode']='bi-linear'
		self.ellipse_pm['usclip']=3.0
		self.ellipse_pm['lsclip']=3.0
		self.ellipse_pm['nclip']=0
		self.ellipse_pm['fflag']=0.5
		self.ellipse_pm['mag0']=0.0
		self.ellipse_pm['refer']=1.0
		self.ellipse_pm['zerolevel']=0.0
		self.ellipse_pm['interactive'] = False

	def __init_dophot_pm(self):
		'''
		initiate the pm data for dophot_C and dophot_fortran
		'''
		self.dophot_C_pm = OrderedDict()

		self.dophot_C_pm['FWHM']    = None    #Approx FWHM of objects (pixels) along major axis.
		self.dophot_C_pm['SKY']     = None
		self.dophot_C_pm['EPERDN']  = None
		self.dophot_C_pm['RDNOISE'] = None
		self.dophot_C_pm['AXIS_RATIO']  = 1.0 #The initial guess of the ratio of the minor axis divided by the major axis for star objects.
		self.dophot_C_pm['TILT']        = 0.0 #The initial guess of the position angle of the major axis with respect to the positive x-axis for a typical unblended stellar object on an image

		self.dophot_C_pm['APBOX_X']     = 16.0#The sizes of the sides of the aperture photometry box in the x direction.
		self.dophot_C_pm['APBOX_Y']     = 16.0#The sizes of the sides of the aperture photometry box in the y direction.
		self.dophot_C_pm['MASKBOX_X']   = 5   # Size of mask box size in x.
		self.dophot_C_pm['MASKBOX_Y']   = 5   # Size of mask box size in y.
		self.dophot_C_pm['NFITBOX_X']   = 12.0#The sizes of the sides of the fitting box used for finding objects in the x direction.
		self.dophot_C_pm['NFITBOX_Y']   = 12.0#The sizes of the sides of the fitting box used for finding objects in the y direction.

		self.dophot_C_pm['IBOTTOM'] = -50        #Lowest allowed data value in data numbers.
		self.dophot_C_pm['ITOP'] = 65535         #Maximum allowed data value in data numbers.
		self.dophot_C_pm['THRESHMIN'] = 100.0    #The value above sky of the lowest threshold that DoPHOT uses to search for objects, in DN
		self.dophot_C_pm['THRESHMAX'] = 40000.0  #The maximum **possible** first threshold that DoPHOT uses to search for objects, in DN.
		self.dophot_C_pm['THRESHDEC'] = 1.0      #Threshold decrement in powers-of-2. For THRESHMAX=40000,THRESHMIN=100, the thresholds (relative to the sky value) would be
							 #25600,12800,6400,3200,1600,800,400,200,100
		self.dophot_C_pm['DOFINALFIT'] = 'YES'   #Do a last fitting iteration after all objects found  fully fits faintest objects of type 3, but slows DoPHOT
		self.dophot_C_pm['THRESHEMP'] = 0.0      #Fit empirical PSF at and below this value.
		self.dophot_C_pm['MAX_SOUGHT'] = 32768   #Quit after this number of improved stars.
		self.dophot_C_pm['MAX_PERF']  = 2000     #Only average up to this number of stars to get shape parameters.
		self.dophot_C_pm['RANGE_MAG'] = 30.      #Sets threshhold some magnitudes fainter than brightest fit.

		self.dophot_C_pm['AUTOSCALE'] = 'NO'     #Auto-scaling of sizes by FWHM.
		self.dophot_C_pm['AUTOTHRESH'] = 'NO'    #Auto-scaling of thresholds.
		self.dophot_C_pm['SCALEFITBOX'] = 3.0    #Size of fit box in units of FWHM.
		self.dophot_C_pm['FITBOXMIN'] = 5.0      #Smallest allowed fit box size.
		self.dophot_C_pm['SCALEAPBOX'] = 6.0     #Size of aperture phot box in units of FWHM.
		self.dophot_C_pm['APBOXMIN'] = 7.0       #Smallest allowed aperture phot box size.
		self.dophot_C_pm['SCALEMASKBOX'] = 1.5   #Size of mask box in units of FWHM.
		self.dophot_C_pm['AMASKBOXMIN'] = 5.0    #Smallest allowed mask box size.
		self.dophot_C_pm['SIGMAIBOTTOM'] = 10.0  #Level of IBOTTOM below sky in units of noise.
		self.dophot_C_pm['SIGMATHRESHMIN'] = 2.0 #Level of THRESHMIN above sky in units of noise.
		self.dophot_C_pm['FIXPOS'] = 'NO'        #Fix star positions?

		self.dophot_C_pm['PARAMS_DEFAULT'] = 'conf/paramdefault_dophot_C'
		self.dophot_C_pm['PARAMS_OUT']    = None            #Output parameters file name.
		self.dophot_C_pm['IMAGE_IN']      = None            #Input image name.
		self.dophot_C_pm['IMAGE_OUT']     = None            #Output image name.
		self.dophot_C_pm['EMP_SUBRAS_OUT']= None            #Empirical PSF subraster (most recent).
		self.dophot_C_pm['OBJECTS_IN']    = None            #Input object list file name.
		self.dophot_C_pm['OBJECTS_OUT']   = None            #Output object list file name.
		self.dophot_C_pm['SHADOWFILE_IN'] = None            #Input shadow file name.
		self.dophot_C_pm['SHADOWFILE_OUT']= None            #Output shadow file name.
		self.dophot_C_pm['ERRORS_OUT']    = None            #Errors of fit to be output if shadow file is requested
		self.dophot_C_pm['LOGFILE'] = None                #Log file name.  'TERM' for screen but failed in practice
		self.dophot_C_pm['LOGVERBOSITY'] = 1

		self.dophot_C_pm['ICRIT']      = 10        #Obliterate if # of pixels > ITOP exceeds this.
		self.dophot_C_pm['CENTINTMAX'] = 40000.0   #Obliterate if central intensity exceeds this.
		self.dophot_C_pm['CTPERSAT']   = 6.0e4     #Assumed intensity for saturated pixels.
		self.dophot_C_pm['NBADLEFT']   = 0         #Ignore pixels closer to the left edge than this.
		self.dophot_C_pm['NBADRIGHT']  = 0         #Ignore pixels closer to the right edge than this.
		self.dophot_C_pm['NBADTOP']    = 0         #Ignore pixels closer to the top edge than this.
		self.dophot_C_pm['NBADBOT']    = 0         #Ignore pixels closer to the bottom edge than this.

		self.dophot_C_pm['PSFTYPE']     = 'PGAUSS'   #PSF type: (PGAUSS, GAUSS, EXTPGAUSS)
		self.dophot_C_pm['SKYTYPE']     = 'PLANE'    #SKY type: (PLANE, HUBBLE, MEDIAN)
		self.dophot_C_pm['OBJTYPE_OUT'] = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE
		self.dophot_C_pm['OBJTYPE_IN']  = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE

		self.dophot_C_pm['EMP_STAR_X'] = 0         #X position of empirical template.
		self.dophot_C_pm['EMP_STAR_Y'] = 0         #Y position of empirical template.
		self.dophot_C_pm['EMP_STAR_Z'] = 0         #Central intensity of empirical template.

		self.dophot_fortran_pm = OrderedDict()
		self.dophot_fortran_pm['SKY']     =  None	#Approximate mean sky value in data numbers.
		self.dophot_fortran_pm['FWHM']    =  None	#Approx FWHM of objects (pixels) along major axis.
		self.dophot_fortran_pm['EPERDN']  =  None	#Electrons per data number.
		self.dophot_fortran_pm['RDNOISE'] =  None	#Readout noise in electrons.
		self.dophot_fortran_pm['TOP'] = 65535		#Maximum allowed data value in data numbers.
		self.dophot_fortran_pm['AUTOTHRESH'] = 'NO'		#psf photometry control parameters
		self.dophot_fortran_pm['THRESHMIN'] = 100.0    #The value above sky of the lowest threshold that DoPHOT uses to search for objects, in DN
		self.dophot_fortran_pm['THRESHMAX'] = 40000.0  #The maximum **possible** first threshold that DoPHOT uses to search for objects, in DN.
		self.dophot_fortran_pm['PARAMS_DEFAULT'] = 'conf/paramdefault_dophot_fortran'
		self.dophot_fortran_pm['PSFTYPE']        = 'PGAUSS' #for C version, there is a new model EXTPGAUSS
		self.dophot_fortran_pm['OBJTYPE_OUT']    = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE
		self.dophot_fortran_pm['OBJTYPE_IN']     = 'COMPLETE' #INTERNAL, COMPLETE, INCOMPLETE

	def __init_hotpants_options(self):
		'''
		All command line options for hotpants
		'''
		self.hotpants_pars = OrderedDict()
		self.hotpants_pars['-inim']  =None  # : comparison image to be differenced
		self.hotpants_pars['-tmplim']=None  # template image
		self.hotpants_pars['-outim'] =None  #output difference image
		self.hotpants_pars['-tu'] = 65000  # upper valid data count, template (25000)
		self.hotpants_pars['-tuk']= 65000  # upper valid data count for kernel, template (tuthresh)
		self.hotpants_pars['-tl'] = 0      # lower valid data count, template (0)
		self.hotpants_pars['-tg'] = 1      # gain in template (1)
		self.hotpants_pars['-tr'] = 0      # e- readnoise in template (0)
		self.hotpants_pars['-tp'] =0       # ADU pedestal in template (0)
		self.hotpants_pars['-tni']=None    # input template noise array (undef)
		self.hotpants_pars['-tmi']=None    # input template mask image (undef)
		self.hotpants_pars['-iu'] = 65000   # upper valid data count, image (25000)
		self.hotpants_pars['-iuk']= 65000   # upper valid data count for kernel, image (iuthresh)
		self.hotpants_pars['-il'] = 0       # lower valid data count, image (0)
		self.hotpants_pars['-ig'] = 1       # gain in image (1)
		self.hotpants_pars['-ir'] = 0       # e- readnoise in image (0)
		self.hotpants_pars['-ip'] = 0       # ADU pedestal in image (0)
		self.hotpants_pars['-ini']=None    # input image noise array (undef)
		self.hotpants_pars['-imi']=None    # input image mask image (undef)
		self.hotpants_pars['-ki'] =None    # use kernel table in image header (undef)
		self.hotpants_pars['-r']  = 10      # convolution kernel half width (10)
		self.hotpants_pars['-kcs']= 21      # size of step for spatial convolution (2 * rkernel + 1)
		self.hotpants_pars['-ft'] = 20.0    # RMS threshold for good centroid in kernel fit (20.0)
		self.hotpants_pars['-sft']= 0.5     # scale fitthresh by this fraction if... (0.5)
		self.hotpants_pars['-nft']= 0.1     # this fraction of stamps are not filled (0.1)
		self.hotpants_pars['-mins']=1.0    # Fraction of kernel half width to spread input mask (1.0)
		self.hotpants_pars['-mous']=1.0    # Ditto output mask, negative = no diffim masking (1.0)
		self.hotpants_pars['-omi']=None    # Output bad pixel mask (undef)
		self.hotpants_pars['-gd'] =None    # only use subsection of full image (full image); [xmin xmax ymi ymax]
		self.hotpants_pars['-nrx']=1       # number of image regions in x dimension (1)
		self.hotpants_pars['-nry']=1       # number of image regions in y dimension (1)
		self.hotpants_pars['-nsx']=10      # number of each region's stamps in x dimension (10)
		self.hotpants_pars['-nsy']=10      # number of each region's stamps in y dimension (10)
		self.hotpants_pars['-afssc'] =1    # autofind stamp centers so #=-nss when -ssf,-cmp (1)
		self.hotpants_pars['-nss'] =3       # number of centroids to use for each stamp (3)
		self.hotpants_pars['-rss'] = 15     # half width substamp to extract around each centroid (15)
		self.hotpants_pars['-savexy']=None  # save positions of stamps for convolution kernel (undef)
		self.hotpants_pars['-c'] =None      # force convolution on (t)emplate or (i)mage (undef)
		self.hotpants_pars['-n'] = 't'      # normalize to (t)emplate, (i)mage, or (u)nconvolved (t)
		self.hotpants_pars['-fom'] = 'v'    # (v)ariance, (s)igma or (h)istogram convolution merit (v)
		self.hotpants_pars['-sconv']= False # all regions convolved in same direction (False)
		self.hotpants_pars['-ko'] =2        # spatial order of kernel variation within region (2)
		self.hotpants_pars['-bgo']=1        # spatial order of background variation within region (1)
		self.hotpants_pars['-ssig'] =3      # threshold for sigma clipping statistics  (3.0)
		self.hotpants_pars['-ks'] =2        # high sigma rejection for bad stamps in kernel fit (2.0)
		self.hotpants_pars['-kfm']=0.990    # fraction of abs(kernel) sum for ok pixel (0.990)
		self.hotpants_pars['-okn']= False   # rescale noise for 'ok' pixels (False)
		self.hotpants_pars['-fi'] =1.0e-30  # value for invalid (bad) pixels (1.0e-30)
		self.hotpants_pars['-fin']=0.0e00   # noise image only fillvalue (0.0e+00)
		self.hotpants_pars['-convvar']=False# convolve variance not noise (False)
		self.hotpants_pars['-oni']=None     # output noise image (undef)
		self.hotpants_pars['-ond']=None     # output noise scaled difference image (undef)
		self.hotpants_pars['-nim']= False   # add noise image as layer to sub image (False)
		self.hotpants_pars['-ndm']= False   # add noise-scaled sub image as layer to sub image (False)
		self.hotpants_pars['-oci']=None     # output convolved image (undef)
		self.hotpants_pars['-cim']= False   # add convolved image as layer to sub image (False)
		self.hotpants_pars['-allm']=False   # output all possible image layers (False)
		self.hotpants_pars['-nc'] = False   # do not clobber output image (False)
		self.hotpants_pars['-hki']= False   # print extensive kernel info to output image header (False)
		self.hotpants_pars['-oki']=None     # new fitsfile with kernel info (under)
		self.hotpants_pars['-sht']= False   # output images 16 bitpix int, vs -32 bitpix float (False)
		self.hotpants_pars['-obs']=1.0      # if -sht, output image BSCALE, overrides -inim (1.0)
		self.hotpants_pars['-obz']=0.0      # if -sht, output image BZERO , overrides -inim (0.0)
		self.hotpants_pars['-nsht']= False  # output noise image 16 bitpix int, vs -32 bitpix float (False)
		self.hotpants_pars['-nbs']=1.0      # noise image only BSCALE, overrides -obs (1.0)
		self.hotpants_pars['-nbz']=0.0      # noise image only BZERO,  overrides -obz (0.0)
		self.hotpants_pars['-v']= 1         # level of verbosity, 0-2 (1)
		self.hotpants_pars['-ng'] = '3 6 0.70 4 1.50 2 3.00'
					#[ngauss degree0 sigma0 .. degreeN sigmaN]
					# ngauss = number of gaussians which compose kernel (3)
					# degree = degree of polynomial associated with gaussian # (6 4 2)
					# sigma  = width of gaussian # (0.70 1.50 3.00)

	def __init_photometry_info(self):
		'''
		initiate the photometry info dictionary
		'''
		self.photometry_info = OrderedDict()
		self.photometry_info_keys = ['name','realimg','flt','obstime','camera', 'exptime','bitpix', 'NX', 'NY', 'flt2int','fixbad','rmCR',\
					     'fwhm','bkg','airmass','nstar','template','x','y','instmag','instmagerr','relmag','relmagerr', 'magzpt','magzpterr', 'calmag','calmagerr','drop']
		self.photometry_info_dtypes = ['S10','S100','S10','f8','S30', 'f8','i4', 'i4', 'i4', 'i4','i4','i4','f8','f8','f8','i4','i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','i4']


		if self.image_subtraction:
			self.photometry_info_keys.append('geotranfac')
			self.photometry_info_keys.append('fwhmconv')
			self.photometry_info_keys.append('convon')
			self.photometry_info_keys.append('xsub')
			self.photometry_info_keys.append('ysub')
			self.photometry_info_keys.append('fitsky')
			for dtype in ['f8', 'f8', 'S10', 'f4', 'f4', 'S3']:
				self.photometry_info_dtypes.append(dtype)
			
			self.photometry_info_keys.append('magap1')  #aperture photometry
			self.photometry_info_keys.append('magerrap1')
			self.photometry_info_keys.append('mskyap1')
			self.photometry_info_keys.append('magap2')
			self.photometry_info_keys.append('magerrap2')
			self.photometry_info_keys.append('mskyap2')
			self.photometry_info_keys.append('mag1')  #use above aperture phot msky and peak fit 
			self.photometry_info_keys.append('magerr1')
			self.photometry_info_keys.append('mag2')
			self.photometry_info_keys.append('magerr2')
			self.photometry_info_keys.append('magfit1')
			self.photometry_info_keys.append('magerrfit1')
			self.photometry_info_keys.append('mskyfit1')
			self.photometry_info_keys.append('magfit2')
			self.photometry_info_keys.append('magerrfit2')
			self.photometry_info_keys.append('mskyfit2')
			N = 16
			for i in range(N): 
				self.photometry_info_dtypes.append("f8")

			self.photometry_info_keys.append('magdiff')
			self.photometry_info_dtypes.append("f8")

		self.photometry_info_init_values = {}
		self.photometry_info_init_values['realimg'] = '?'
		self.photometry_info_init_values['flt'] = '?'
		self.photometry_info_init_values['obstime'] = 0
		self.photometry_info_init_values['camera'] = '?'
		self.photometry_info_init_values['exptime'] = 0
		self.photometry_info_init_values['bitpix'] = 0
		self.photometry_info_init_values['NX']  = 0
		self.photometry_info_init_values['NY']  = 0
		self.photometry_info_init_values['flt2int'] = 0
		self.photometry_info_init_values['fixbad'] = 0
		self.photometry_info_init_values['rmCR'] = 0	#remove cosmic ray
		self.photometry_info_init_values['fwhm'] = 0
		self.photometry_info_init_values['bkg'] = 0
		self.photometry_info_init_values['airmass'] = 0
		self.photometry_info_init_values['nstar'] = 0
		self.photometry_info_init_values['template'] = 0
		self.photometry_info_init_values['x'] = 0.0
		self.photometry_info_init_values['y'] = 0.0
		self.photometry_info_init_values['instmag'] = 99.99
		self.photometry_info_init_values['instmagerr'] = 99.99
		self.photometry_info_init_values['relmag'] = 99.99
		self.photometry_info_init_values['relmagerr'] = 99.99
		self.photometry_info_init_values['magzpt'] = 99.99
		self.photometry_info_init_values['magzpterr'] = 99.99
		self.photometry_info_init_values['calmag'] = 99.99
		self.photometry_info_init_values['calmagerr'] = 99.99
		self.photometry_info_init_values['drop'] = 0
		self.photometry_info_init_values['geotranfac'] = 0
		self.photometry_info_init_values['fwhmconv'] = 0
		self.photometry_info_init_values['convon'] = '?'
		self.photometry_info_init_values['xsub'] = 0.0
		self.photometry_info_init_values['ysub'] = 0.0
		self.photometry_info_init_values['fitsky'] = '?'
		self.photometry_info_init_values['magap1'] = 99.99
                self.photometry_info_init_values['magerrap1'] = 99.99
                self.photometry_info_init_values['mskyap1'] = 0
                self.photometry_info_init_values['magap2'] = 99.99
                self.photometry_info_init_values['magerrap2'] = 99.99
                self.photometry_info_init_values['mskyap2'] = 0
                self.photometry_info_init_values['mag1'] = 99.99
                self.photometry_info_init_values['magerr1'] = 99.99
                self.photometry_info_init_values['mag2'] = 99.99
                self.photometry_info_init_values['magerr2'] = 99.99
                self.photometry_info_init_values['magfit1'] = 99.99
                self.photometry_info_init_values['magerrfit1'] = 99.99
                self.photometry_info_init_values['mskyfit1'] = 0
                self.photometry_info_init_values['magfit2'] = 99.99
                self.photometry_info_init_values['magerrfit2'] = 99.99
                self.photometry_info_init_values['mskyfit2'] = 0
		self.photometry_info_init_values['magdiff'] = 99.99



		for img in self.images:
			self.photometry_info[img] = OrderedDict()
			self.photometry_info[img]['name'] = img
			for infokey in self.photometry_info_keys:
				if infokey == 'name':
					continue
				self.photometry_info[img][infokey]   = self.photometry_info_init_values[infokey]

		self.anomaly_signs = {}
		self.anomaly_signs['bkg'] =  -9999.99
		self.anomaly_signs['fwhm'] = -99.99
		self.anomaly_signs['nstar'] = -9999

	def __prepare_parameter_files(self):
		'''
		prepare a copy of parameter files for the working target
		'''
		if not os.path.exists(self.parafile_dir):
			os.mkdir(self.parafile_dir)

		for parafile in os.listdir(self.parafile_template_dir):
			parafile_abs_from = os.path.join(self.parafile_template_dir,parafile)
			parafile_abs_to = os.path.join(self.parafile_dir,parafile)
			if not os.path.isfile(parafile_abs_to):
			#self.__delete_file_if_exist(parafile_abs_to)
				shutil.copy(parafile_abs_from,parafile_abs_to)

# Notes manage
	def check_notes(self):
		'''
		there is a notes file under result directory, try to look into that for reduction notes
		'''
		if not os.path.exists(self.notesfile):
			print "Notes file not exists..."
		else:
			for noteline in open(self.notesfile).readlines():
				print noteline

	def add_note(self, note):
		'''
		add note to self.readme
		'''
		self.readme = self.readme + note+'\n'

	def save_notes(self, outfile=None):
		'''
		save notes in self.readme to the outfile; default outfile: self.notesfile, reduction README file
		'''
		if outfile is None:
			outfile = self.notesfile
		nowtime = datetime.datetime.now().isoformat()
		timestamp = '\n[%s]\n'%nowtime
		fid = open(outfile, 'awt')
		fid.write(timestamp)
		fid.write(self.readme)
		fid.close()
		self.readme = ''

#Prepare images
	def ln_image_to_datadir(self,from_dir,to_dir):
		'''
		spaces in filenames can be annoying (cause unexpected errors in bash shell) and remove them first
		'''
		self.__filename_modification_remove_spaces(from_dir)

		ori_images = os.listdir(from_dir)
		lned_images = os.listdir(to_dir)
		num = len(lned_images)

		if num>9999:
			numdig = 5
		elif num > 999:
			numdig = 4
		else:
			numdig = 3

		lned_images_realpath = [os.path.realpath(os.path.join(to_dir,img)) for img in lned_images]
		for image in ori_images:
			image_abs = os.path.realpath(os.path.join(from_dir,image))
			if image_abs not in lned_images_realpath:
				obs_indicator = '0'*(numdig-len(str(num)))+str(num)
				new_lned_image = obs_indicator + '.fits'
				from_image = os.path.join(from_dir,image)
				to_image = os.path.join(to_dir,new_lned_image)

				print "ln -s %s %s"%(from_image,to_image)
				command = "ln -s %s %s"%(from_image,to_image)
				os.system(command)
				num = num+1

	def reln_image_to_datadir(self, verbose=0):
		'''
		re-link the registered images in photometry result table to the workplace folder
		'''
		photret_legacy = self.result_table
		realimgs = photret_legacy['realimg'].data
		dstimgs_retfile = photret_legacy['name'].data

		from_dir = self.repository_dir_current_sn
		to_dir = self.raw_image_dir
		dstimgs_old_all = os.listdir(to_dir)

		for img in dstimgs_old_all:
			dstimg_old_abs = os.path.join(to_dir,img)
                        os.remove(dstimg_old_abs)

		for srcimg,dstimg in zip(realimgs,dstimgs_retfile):
			img_to = os.path.join(to_dir,dstimg)
			img_from = os.path.join(from_dir,srcimg)
			command = "ln -s %s %s"%(img_from,img_to)
			if verbose:
				print command
			os.system(command)

	def __filename_modification_remove_spaces(self,file_dir):
		'''
		remove space in filename
		'''
		filenames = os.listdir(file_dir)
		filenames_new = [filename.replace(" ","") for filename in filenames]
		for fromfile,tofile in zip(filenames,filenames_new):
			if fromfile != tofile:
				src = os.path.join(file_dir,fromfile)
				dst = os.path.join(file_dir,tofile)
				os.rename(src,dst)

	def image_flt2int(self, which_dir = 'modified_image', nocheck=False):
		'''
		convert data type from float to integer for those with bitpix == -32
		'''
		try:
			for key in self.photometry_info.keys():
				if nocheck:
					self.__flt2int(key, which_dir = which_dir )
				else:
					if (self.photometry_info[key]['bitpix'] == -32 or self.photometry_info[key]['bitpix'] == -64) and (not self.photometry_info[key]['flt2int']):
						self.__flt2int(key, which_dir = which_dir )
		except:
			print "flt2int not available..."

		print "CCDProc flt2int... Done!"

	def __flt2int(self, image_key, which_dir='modified_image'):
        	'''
        	convert the data type from float32 to int to meet the requirement of some functions in diapl
        	'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		output_image_temp = os.path.join(self.modified_image_dir, 'int16_temp.fits')
		output_image = os.path.join(self.modified_image_dir,image_key)

		flt2int = os.path.join(self.base.ccdproc_dir, 'flt2int')
        	command = '%s %s %s'%(flt2int,input_image,output_image_temp)
		try:
			os.system(command)
			status = 1
			if os.path.isfile(output_image):
				os.remove(output_image)
			os.rename(output_image_temp, output_image)

		except:
			print "skip image %s in flt2int"%image_key
			status = 0

		self.photometry_info[image_key]['flt2int'] = status

#Get basic image info
	def __get_fits_image_size(self, image):
		'''
		get image size from fits header keywords 'NAXIS1' and 'NAXIS2'
		'''
		hdu = fits.open(image)
		hdr = hdu['PRIMARY'].header
		NX = hdr['NAXIS1']
		NY = hdr['NAXIS2']
		hdu.close()
		return NX,NY

	def __renew_bitpix_info(self,image_key):
		img_abs = self.images[image_key]
                bitpix_key = self.base.telescopes_info[self.current_telescope]['bitpixkey']
                self.photometry_info[image_key]['bitpix'] = get_bitpix(img_abs,bitpix_key,self.current_telescope)

	def __renew_flt_info(self,image_key):
		img_abs = self.images[image_key]
		flt_key = self.base.telescopes_info[self.current_telescope]['fltkey']
		try:
			flt = get_filter(img_abs,flt_key,self.current_telescope)
			self.photometry_info[image_key]['flt'] = flt
		except:
			self.photometry_info[image_key]['drop'] = 12

	def __renew_airmass_info(self, image_key):
		img_abs = self.images[image_key]
		airmass_key = self.base.telescopes_info[self.current_telescope]['airmasskey']
                self.photometry_info[image_key]['airmass'] = get_airmass_telescope(img_abs,airmass_key,self.current_telescope)

	def __get_airmass_astroobs_single(self, imgkey):
		obsJD = self.photometry_info[imgkey]['obstime']
		airmass_ret = get_airmass_given_time_target_observatory(obsJD, self.current_telescope, self.sn_ra_str, self.sn_dec_str)
		return airmass_ret

	def get_airmass(self):
		airmass_key = self.base.telescopes_info[self.current_telescope]['airmasskey']
		if airmass_key != '':
			for img in self.images.keys():
				img_abs = self.images[img]
				self.photometry_info[img]['airmass'] = get_airmass_telescope(img_abs,airmass_key,self.current_telescope)
		else:
			for img	in self.images.keys():
				airmass_ret = self.__get_airmass_astroobs_single(img)
				self.photometry_info[img]['airmass'] = airmass_ret

	def __renew_obstime_info(self,image_key):
		img_abs = self.images[image_key]
		obstime_key = self.base.telescopes_info[self.current_telescope]['obstimekey']
		self.photometry_info[image_key]['obstime'] = get_obstime(img_abs,obstime_key,self.current_telescope)

	def __renew_exptime_info(self,image_key):
		img_abs = self.images[image_key]
		exptime_key = self.base.telescopes_info[self.current_telescope]['exptimekey']
		self.photometry_info[image_key]['exptime'] = get_exptime(img_abs,exptime_key,self.current_telescope)

	def __get_image_size_from_fitsheader(self, image_key, ext=0, verbose=1):
		'''
		this is to get image size from the fits file header
		'''
		img_abs = self.images[image_key]
		hdulist = fits.open(img_abs)
		hdr0 = hdulist[0].header
		if 'NAXIS1' in hdr0 and 'NAXIS2' in hdr0:
			self.photometry_info[image_key]['NX'] = hdr0['NAXIS1']
			self.photometry_info[image_key]['NY'] = hdr0['NAXIS2']
		else:
			if len(hdulist)>1:
				hdr1 = hdulist[1].header
				if 'NAXIS1' in hdr1 and 'NAXIS2' in hdr1:
					self.photometry_info[image_key]['NX'] = hdr1['NAXIS1']
					self.photometry_info[image_key]['NY'] = hdr1['NAXIS2']
				else:
					if verbose:
						print "image size NAXIS1 and NAXIS2 not found in %s"%img_abs	

	def renew_image_info(self,image_key=None,renewall=False,):
		'''
		renew image information: filter,obstime and bixpix
		'''
		if renewall:
			for key in self.photometry_info.keys():
				self.__renew_flt_info(key)
				self.__renew_obstime_info(key)
				self.__renew_bitpix_info(key)
				self.__renew_exptime_info(key)
		elif image_key is not None and image_key in self.images.keys():
			self.__renew_flt_info(image_key)
			self.__renew_obstime_info(image_key)
			self.__renew_bitpix_info(image_key)
			self.__renew_exptime_info(key)
		else:
			print "no valid input for this action..."

#refine images
	def remove_cosmic_ray(self, which_dir = 'modified_image', method='ccdproc',sigclip=5, skipdroped=True):
		'''
		remove cosmic ray
		INPUTS:
			which_dir: the directory containing the input images
			method: 'ccdproc' or 'cosmics'
		'''

		for imgkey in self.photometry_info.keys():
			if self.photometry_info[imgkey]['drop'] != 0 and skipdroped:
				continue
			try:
				if not self.photometry_info[imgkey]['rmCR']:
					self.remove_cosmic_ray_single_image(imgkey, which_dir = which_dir, method=method, sigclip=sigclip)
					self.photometry_info[imgkey]['rmCR'] = 1
			except Exception as e:
				print e
				print "Remove cosmic ray failure on image %s" %imgkey
				self.photometry_info[imgkey]['rmCR'] = 0

	def remove_cosmic_ray_single_image(self,image_key,which_dir = 'modified_image', method='ccdproc', output_image = None,single_use = False, image_gain = 1.0, image_readnoise=10.0, sigclip =5):
		'''
		INPUTS:
			image_key:
			which_dir:
			method: cosmics or ccdproc
			output_image:
			single_use:
			image_gain:
			image_read_noise:
			sigclip:
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir,image_key)

		data = fits.open(input_image)[0].data
		if np.isnan(data[0,0]):
			self.photometry_info[image_key]['drop'] = 99
			return			

		# Build the object:
		if len(self.base.telescopes_info[self.current_telescope]['gain']):
			gain = float(self.base.telescopes_info[self.current_telescope]['gain'])
		else:
			gain = image_gain

		if len(self.base.telescopes_info[self.current_telescope]['readnoise']):
			readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
		else:
			readnoise = image_readnoise

		if method == 'cosmics':
			array, header = cosmics.fromfits(input_image) # Read the FITS; # array is a 2D numpy array
			c = cosmics.cosmicsimage(array, gain=gain, readnoise=readnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0) # There are other options, check the manual...
			c.run(maxiter = 4) # Run the full artillery :

			# Write the cleaned image into a new FITS file, conserving the original header :
			if os.path.isfile(output_image):
				os.remove(output_image)
			cosmics.tofits(output_image, c.cleanarray, header)
			# If you want the mask, here it is :
			#mask_filename = ''
			#cosmics.tofits(mask_filename, c.mask, header)
			# (c.mask is a boolean numpy array, that gets converted here to an integer array)
		elif method == 'ccdproc':
			print input_image, output_image, gain, readnoise
			if os.path.exists(output_image):
				output_temp = "temp_clean_image_%s.fits"%self.current_sn
				if os.path.exists(output_temp):
					os.remove(output_temp)
				movedata = 1
			else:
				output_temp = output_image
				movedata = 0
			remove_CR_ccdproc_cosmicray_lacosmic(input_image, output_temp, gain, readnoise, sigclip=sigclip)
			if movedata:
				shutil.move(output_temp, output_image)
		else:
			raise IOError("method of %s not supported"%method)
		if single_use:
			self.photometry_info[image_key]['rmCR'] = 1

	def fix_bad_pixels(self,mode='vicinity', which_dir='raw_image', bad_threshold = 60000, dolist=None, skiplist=None):
		'''
		replace the 'bad' pixels with 'good' values which determined by the 'fix_mode'

		Inputs:
			mode: how you want to fix the 'bad' pixels, available options are the followings,
				'vicinity', replace the 'bad' values with the values from the vicinity
				'bkg', replace the 'bad' values with the rough background
				'zero', replace the 'bad' values with zeros
				'interp', replace with 'bad' values with the interpolated values

		'''
		print "Fix bad pixels on image..."
		for imgkey in self.photometry_info.keys():
			print imgkey
			if dolist is not None:
				if imgkey not in dolist:
					continue
			if skiplist is not None:
				if imgkey in skiplist:
					continue
			try:
				self.__fix_bad_pixels(imgkey,mode=mode, which_dir = which_dir, bad_threshold = bad_threshold )
				self.photometry_info[imgkey]['fixbad'] = 1
			except:
				print "Fix bad pixels failure on image %s" %imgkey
				self.photometry_info[imgkey]['fixbad'] = 0
		print "Fix bad pixels on image... Done!"

	def fix_bad_pixels_with_given_mask_region(self, imgkey, which_dir='raw', fillmethod='bkg', fillvalue=None, outfile=None):
		'''
		if 
		'''
		maskfile = self._photometry__get_internal_image(imgkey, which_dir='mask')	
		maskdata = fits.open(maskfile)[0].data
		imgfile = self._photometry__get_internal_image(imgkey, which_dir=which_dir)
		imghdu = fits.open(imgfile)
		imgdata = imghdu[0].data
		if fillmethod == 'bkg':
			fillvalue = self.photometry_info[imgkey]['bkg']
		imgdata[np.where(maskdata==1)] = fillvalue
		imghdu[0].data = imgdata
		if outfile is None:
			outfile = self._photometry__get_internal_image(imgkey, which_dir='md')
		imghdu.writeto(outfile, overwrite=True)


	def fix_negative_values_pixels(self, which_dir ='raw_image', fillmethod='bkg', verbose=1):
		'''
		fix the negative pixel value issue
		'''
		for img in self.images.keys():
			if fillmethod == 'bkg':
				bkgvalue = self.photometry_info[img]['bkg']
			self.__replace_negative_value(img, which_dir=which_dir, fillvalue=bkgvalue)
			if verbose:
				print "negative values in %s will be replaced with %s"%(img, bkgvalue)

	def __replace_negative_value(self, image_key, which_dir = 'raw_image', fillvalue=None, fillmethod='bkg', output_image=None):
		'''
		replace the negative values with fillvalue or the value derived from fillmethod
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir,image_key)

		hdu = fits.open(input_image)
		header = hdu[0].header
		data = hdu[0].data
		hdu.close()
		bad_pixels = self.__get_bad_pixels(data,0, mode='negative')
		data_new = self.__fill_bad_pixels(data,bad_pixels, fill_value=fillvalue, mode=fillmethod)

		hdu_new = fits.PrimaryHDU(data_new,header)
		hdulist = fits.HDUList([hdu_new])

		if os.path.isfile(output_image):
			os.remove(output_image)
		hdulist.writeto(output_image)


	def __fix_bad_pixels_star_model(self, image_key, mode='bkg', which_dir = 'raw_image', bad_threshold = 60000, data_region=None, output_image = None,single_use = False, verbose=0):
		'''
		this is similar with self._pixels_star_model but have the ability to deal with situation where one bright star is within the interpolation region
		'''
		while True:
			tmpimg = 'temp_fixbad.fits'
			self.__fix_bad_pixels(image_key, mode=mode, which_dir = which_dir, bad_threshold = bad_threshold, data_region=data_region, output_image = tempimg)
			self.iraf_photometry_peak(image_key, image=tmpimg)


	def __fix_bad_pixels(self, image_key, mode='bkg', which_dir = 'raw_image', bad_threshold = 60000, data_region=None, output_image = None,single_use = False, verbose=0):
		'''
		Replace the saturated pixels with normal value which determined by the 'mode'
		data_region: xl,xu,yl,yu --> [xl:xu, yl:yu]
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		if output_image is None:
			output_image = os.path.join(self.modified_image_dir,image_key)

		hdu = fits.open(input_image)
		header = hdu[0].header
		data = hdu[0].data
		hdu.close()
		if data_region is not None:
			try:
				a1l,a1u,a2l,a2u = data_region
			except Exception as e:
				raise ValueError(e)
			datacut = data[a1l:a1u,a2l:a2u]
			bad_pixels = self.__get_bad_pixels(datacut,bad_threshold)
			print bad_pixels.shape
			print datacut.shape
			data_filled = self.__fill_bad_pixels(datacut,bad_pixels,mode=mode, verbose=verbose)
			print data_filled.shape
			data[a1l:a1u,a2l:a2u] = data_filled
			data_new = data
		else:	
			bad_pixels = self.__get_bad_pixels(data,bad_threshold)
			data_new = self.__fill_bad_pixels(data,bad_pixels,mode=mode, verbose=verbose)

		hdu_new = fits.PrimaryHDU(data_new,header)
		hdulist = fits.HDUList([hdu_new])
		if os.path.isfile(output_image):
			os.remove(output_image)
		hdulist.writeto(output_image)

		if single_use:
			self.photometry_info[image_key]['fixbad'] = 1

#FWHM and background
	def get_fwhm_bkg(self, method = 'diapl', which_dir = 'modified_image', time_out_max = 60, renew =False, skipdroped=False):
		'''
		get fwhm and background

		INPUTs:
			method: 'diapl' or 'fistar' or 'borrow'
			which_dir: 'modified_image' or 'raw_image'
			time_out_max: the maximum time before time out for single image
		'''
		for imgkey in self.photometry_info.keys():
			if self.photometry_info[imgkey]['drop'] != 0 and skipdroped:
				continue
			if self.photometry_info[imgkey]['fwhm'] == 0 or renew:
				self.__get_fwhm_bkg(imgkey,method =method, which_dir = which_dir, time_out_max=time_out_max)

	def renew_fwhm_bkg(self, method ='diapl', which_dir = 'modified_image', time_out_max = 60):
		'''
		get fwhm and background
		'''
		for key in self.photometry_info.keys():
                        self.__get_fwhm_bkg(key, method =method, which_dir = which_dir, time_out_max = time_out_max)

	def __get_fwhm_bkg(self,imgkey,method='diapl', which_dir = 'modified_image', time_out_max = 60):
		'''
		get the fwhm and background of given image with imagekey 'imgkey'

		INPUTS:
			imgkey:
			method: 'fwhm' from diapl, or 'fistar' from fitsh, 'borrow' is using results from external reductions for example from WFCAM CASU reduction pipeline which give SEEING and SKYLEVEL
			which_dir: 'raw_image' or 'modified_image'
			time_out_max: time out threshold for this task
		'''

		if method == 'diapl':
			success,result = self.__get_fwhm_bkg_single_diapl(imgkey, which_dir = which_dir, time_out_max = time_out_max)
			if success:
				fwhm  = round(np.float(result[-2]),2)
				bkg   = round(np.float(result[-3]),2)
				print imgkey,fwhm,bkg
			else:
				fwhm  = self.anomaly_signs['fwhm']
				bkg   = self.anomaly_signs['bkg']
				self.photometry_info[imgkey]['drop'] = 3
		elif method == 'fistar':
			success, result = self.__get_fwhm_bkg_single_fistar(imgkey, which_dir = which_dir, time_out_max = time_out_max)
			if success:
				bkg  = round(result[0],2)
				fwhm = round(result[1],2)
				print imgkey,fwhm,bkg
			else:
				fwhm  = self.anomaly_signs['fwhm']
				bkg   = self.anomaly_signs['bkg']
				self.photometry_info[imgkey]['drop'] = 3
		elif method == 'borrow':
			result = self.__get_fwhm_bkg_single_borrow(imgkey, which_dir =which_dir)
			bkg  = round(result[0],2)
			fwhm = round(result[1],2)
			print imgkey,fwhm,bkg
		else:
			raise IOError("method %s not supported"%method)

		self.photometry_info[imgkey]['fwhm'] = fwhm
		self.photometry_info[imgkey]['bkg'] = bkg

		if self.photometry_info[imgkey]['fwhm'] == -99.99:
			self.photometry_info[imgkey]['drop'] = 3

	def __get_fwhm_bkg_single_borrow(self, image_key, which_dir='modified_image'):
		'''
		get fwhm and background by extracting external reduction results
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		fitsinfo = get_fits_info(input_image,['SEEING', 'SKYLEVEL'], extension=None)
		fwhm = fitsinfo['SEEING']
		bkg = fitsinfo['SKYLEVEL']
		outdata = [bkg, fwhm]

		return outdata

	def __get_fwhm_bkg_single_diapl(self,image_key, which_dir = 'modified_image', time_out_max = 60, monitorplot=False, deletemfile=True, update_imginfo=False):
		'''
		get fwhm and background with function 'fwhm' from  diapl package
		'''
		snname = self.current_sn
		templog = image_key.split('.')[0] + '_fwhm_%s.log'%snname
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		input_image_int16 = 'tmpimage_delete_after_%s.fits'%snname
		if os.path.exists(input_image_int16):
			os.remove(input_image_int16)

		if self.photometry_info[image_key]['bitpix'] != 16:
			flt2int = os.path.join(self.base.ccdproc_dir, 'flt2int')
			command = "%s %s %s"%(flt2int, input_image, input_image_int16)
			os.system(command)
		else:
			shutil.copy(input_image, input_image_int16)

		print input_image

		telcode = self.current_telescope
		fwhm_parfile_telspec = os.path.join(self.parafile_dir, 'fwhm_%s.par'%telcode)
		if os.path.exists(fwhm_parfile_telspec):
			fwhm_parfile = fwhm_parfile_telspec
		else:
			fwhm_parfile = os.path.join(self.parafile_dir, 'fwhm.par')

		fwhmstarsfile= os.path.join(self.stars_dir, image_key.split('.')[0]+'_fwhm.stars')

		fwhm = os.path.join(self.base.diapl_dir, 'fwhm_ping2')
		command_line = '%s %s %s %s>>%s' %(fwhm, fwhm_parfile, input_image_int16, fwhmstarsfile, templog)
		print command_line
		command = Command(command_line)
		failure_signal = command.run(timeout=time_out_max)

		if os.path.exists(input_image_int16) and deletemfile:
			os.remove(input_image_int16)

		if failure_signal:
			print "Sorry friend, we can't get the result within %s s"%time_out_max
			success = False
			data = None
		else:
			loglines = open(templog).readlines()
			fwhmline = loglines[-1]
			print fwhmline
			data = fwhmline.strip().split()[-5:]
			success = True

			if monitorplot:
				fwhmdata = np.loadtxt(fwhmstarsfile)
				fig = plt.figure(figsize=(5,5))
				plt.hist(fwhmdata[:,4], bins=30)
				plt.show()

		if update_imginfo and success:
			print "fwhm and bkg in photometry_info will be updated for image %s"%image_key
			fwhm  = round(np.float(data[-2]),2)
			bkg   = round(np.float(data[-3]),2)
			self.photometry_info[image_key]['fwhm'] = fwhm
			self.photometry_info[image_key]['bkg'] = bkg

		if os.path.exists(templog) and deletemfile:
			os.remove(templog)

		return success, data

	def __get_fwhm_bkg_single_fistar(self,imgkey, which_dir = 'modified_image', time_out_max = 60):

		success, info = self.__fistar_action(imgkey, which_dir = which_dir, time_out_max = time_out_max)
		return success, info

	def __fistar_action(self,imgkey, which_dir = 'modified_image', time_out_max = 60):
		'''
		run fitsh/fistar task
		'''
		imgkey_s = imgkey.split('.')[0]

		if which_dir == 'modified_image':
			img_file = os.path.join(self.modified_image_dir, imgkey)
		else:
			img_file = self.images[imgkey]

		print img_file
		hdu = fits.open(img_file)
		img_data = hdu[0].data
		print img_data.shape
		N1,N2 = img_data.shape
		print N1,N2

		from estimate_background import estimate_bkg
		bkgmean, bkgmedian, bkgstd = estimate_bkg(img_data)
		peakthreshold = bkgmedian+3*bkgstd
		source_file = os.path.join(self.stars_dir,imgkey_s + '.temp')

		fistar = os.path.join(self.base.fitsh_dir, 'fistar')
		fistar_command = "%s -i %s -o %s -t %s -s x -F id,x,y,bg,amp,fwhm"%(fistar, img_file,source_file,peakthreshold)
		command = Command(fistar_command)

		failure_signal = command.run(timeout=time_out_max)
		if not failure_signal:
			source_info = np.loadtxt(source_file)
			bkg = np.mean(source_info[:,3])
			fwhm = np.median(source_info[:,5])
			info = [bkg,fwhm]
			success = True
		else:
			print "dear friend, we can't get the result within %s s"%time_out_max
			info = None
			success = False

		return success, info

	def get_fwhm_bkg_2dgaussian_photutils_with_input_star_xys(self, imgkey, xys, fitwidthx=20, fitheighty=20, which_dir='modified_image', update_phot_record=0, verbose=1, display_fit=0):
		'''
		get the fwhm and local bkg with given list of stars by fitting 2d gaussian to it/them.
		INPUTS:
			imgkey:
			xys: star coordinate array with shape of (N,2)
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		imagedata = fits.open(input_image)[0].data

		fwhm_list = []
		bkg_list = []
		for xy in xys:
			try:
				xfwhm, yfwhm, bkg = self.__estimate_fwhm_and_localbkg_single_star_photutils_2dgaussian(imagedata, xy[0], xy[1], fitwidthx, fitheighty, verbose=verbose, display_fit=display_fit)
				fwhm = (xfwhm+yfwhm)/2.0
				fwhm_list.append(fwhm)
				bkg_list.append(bkg)
			except:
				print "failed on star (%s,%s)"%(xy[0], xy[1])
		if len(fwhm_list)>5:
			fwhm = np.median(fwhm_list)
			bkg = np.median(bkg_list)
		else:
			fwhm = np.mean(fwhm_list)
			bkg  = np.mean(bkg_list)

		fwhm = np.round(fwhm, 1)
		bkg = np.round(bkg, 1)
		if update_phot_record:
			self.photometry_info[imgkey]['fwhm'] = fwhm
			self.photometry_info[imgkey]['bkg'] = bkg

		return fwhm, bkg

	def get_fwhm_bkg_2dgaussian_photutils_pick_on_ds9(self, imgkey, nstar=1, fitwidthx=20, fitheighty=20, which_dir='modified_image', update_phot_record=0, verbose=1, display_fit=0):
		'''
		get the fwhm and local bkg with selected star(s) by fitting 2d gaussian to it/them.
		INPUTS:
			imgkey:
			nstar: int, number of stars selected for fitting; or None, continue until no coordinate selection obtained
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		imagedata = fits.open(input_image)[0].data

		newds9 = True
		xy = self.__get_xy_on_image_from_ds9(input_image, regionfile=None, newds9=newds9)
		xfwhm, yfwhm, bkg = self.__estimate_fwhm_and_localbkg_single_star_photutils_2dgaussian(imagedata, xy[0], xy[1], fitwidthx, fitheighty, verbose=verbose, display_fit=display_fit)
		fwhm = (xfwhm+yfwhm)/2.0

		if nstar > 1 or nstar is None:
			newds9 = 0
			istar = 1
			fwhm_list = [fwhm]
			bkg_list = [bkg]
			continue_xy = True
			while continue_xy:
				xy = self.__get_xy_on_image_from_ds9(input_image, regionfile=None, newds9=newds9)
				if xy is None:
					continue_xy = False
				else:
					xfwhm, yfwhm, bkg = self.__estimate_fwhm_and_localbkg_single_star_photutils_2dgaussian(imagedata, xy[0], xy[1], fitwidthx, fitheighty, verbose=verbose, display_fit=display_fit)
					fwhm = (xfwhm+yfwhm)/2.0
					istar +=1
					fwhm_list.append(fwhm)
					bkg_list.append(bkg)
					continue_xy = istar <nstar
			fwhm = np.mean(fwhm_list)
			bkg = np.mean(bkg_list)


		if update_phot_record:
			self.photometry_info[imgkey]['fwhm'] = fwhm
			self.photometry_info[imgkey]['bkg'] = bkg

		return fwhm, bkg

	def __estimate_fwhm_and_localbkg_single_star_photutils_2dgaussian(self, imagedata, cx0, cy0, fitwidthx, fitheighty, verbose=1, display_fit=0, strict_region=False):
		'''
		Fit data within [(cy0-fitheighty/2):(cy0+fitheighty/2), (cx0-fitwidthx/2):(cx0+fitwidthx/2)] of imagedata with 2d gaussian model
		'''
		a1 = int(cx0-fitwidthx/2)
		a2 = int(cx0+fitwidthx/2)
		b1 = int(cy0-fitheighty/2)
		b2 = int(cy0+fitheighty/2)
		NY, NX = imagedata.shape

		if strict_region:
			if b1 < fitheighty/2.0 or a1 < fitwidthx/2.0 or b2 > NY-fitheighty/2.0 or a2 > NX-fitwidthx/2.0:
				raise ValueError("The star is too close to the image edge or outside of the image")
		else:
			if a1 < 0:
				a1 = 0
			if a2 > (NX-1):
				a2 = NX - 1
			if b1 < 0:
				b1 = 0
			if b2 > (NY-1):
				b2 = NY - 1

		fitdata = imagedata[b1:b2,a1:a2]
		gfit = fit_2dgaussian(fitdata)
		if verbose:
			for pname,pvalue in zip(gfit.param_names,gfit.parameters):
				print "%s: %s"%(pname, str(np.round(pvalue,2)))

		x_fwhm = 2*np.sqrt(2.0*np.log(2))*np.abs(gfit.x_stddev.value)
		y_fwhm = 2*np.sqrt(2.0*np.log(2))*np.abs(gfit.y_stddev.value)
		localbkg = gfit.constant.value

		if display_fit:
			show_2d_gaussian_fit_result(fitdata, gfit)

		return x_fwhm, y_fwhm, localbkg

#Source detection
	def source_detection(self,method='daofind', which_dir = None, fistar_ntime_above_bkg = 3, daofind_nsigma_thresh = 5, nmin=5, redo_daofind=True):
		'''
		get the approximate number of stars in the observation images
		INPUTS:
			method: 'sfind' or 'fistar' or 'daofind'
			which_dir:
			fistar_ntime_above_bkg:
		'''
		for imgkey in self.photometry_info.keys():
			if which_dir is None:
				if self.photometry_info[imgkey]['bitpix'] == -32 or self.photometry_info[key]['bitpix'] == -64:
					imgdir  = 'modified_image'
				else:
					imgdir = 'raw_image'
			else:
				imgdir = which_dir

			if self.photometry_info[imgkey]['nstar'] == 0 and self.photometry_info[imgkey]['drop'] == 0:
				self.__source_detection(imgkey,method=method,ntime_above_bkg = fistar_ntime_above_bkg,  daofind_nsigma_thresh = daofind_nsigma_thresh, which_dir = imgdir, drop_nstar_less_than_this=nmin, redo_daofind=redo_daofind)

	def source_detection_flt(self, flt, renew=False, method='daofind', which_dir = 'raw_image', fistar_ntime_above_bkg = 3, daofind_nsigma_thresh = 5, nmin=5):
		'''
		source detection
		'''
		for imgkey in self.images.keys():
			if self.photometry_info[imgkey]['flt'] != flt or self.photometry_info[imgkey]['drop'] != 0:
				continue
			if self.photometry_info[imgkey]['nstar'] != self.photometry_info_init_values['nstar'] and (not renew):
				continue
			self.__source_detection(imgkey,method=method, ntime_above_bkg = fistar_ntime_above_bkg,  daofind_nsigma_thresh = daofind_nsigma_thresh, which_dir = which_dir, drop_nstar_less_than_this=nmin)

	def renew_source_detection(self,method='sfind', which_dir = 'modified_image', fistar_ntime_above_bkg = 3,  daofind_nsigma_thresh = 5):
		'''
		get the approximate number of stars in the observation images
		'''
		for imgkey in self.photometry_info.keys():
			if self.photometry_info[imgkey]['drop'] == 0:
				self.__source_detection(imgkey,method=method,ntime_above_bkg = fistar_ntime_above_bkg,  daofind_nsigma_thresh = daofind_nsigma_thresh, which_dir = which_dir)

	def __source_detection(self,imgkey,method = 'sfind',ntime_above_bkg =3,  daofind_nsigma_thresh = 3, drop_nstar_less_than_this = 10, redo_daofind=True, which_dir = 'modified_image'):

		if method == 'sfind':
			try:
				print imgkey
				starnum = self.__source_detection_single_sfind_diapl(imgkey, which_dir = which_dir)
				self.photometry_info[imgkey]['nstar'] = starnum
				if starnum < drop_nstar_less_than_this:
					self.photometry_info[imgkey]['drop'] = 2
			except Exception,e:
				nstar = self.anomaly_signs['nstar']
				print "source finding or region file saving error on %s"%imgkey
				self.photometry_info[imgkey]['nstar'] = nstar
				self.photometry_info[imgkey]['drop'] = 2
		elif method == 'fistar':
			try:
				print imgkey
				starnum = self.__source_detection_single_fistar_fitsh(imgkey,ntime_above_bkg = ntime_above_bkg, which_dir = which_dir)
				self.photometry_info[imgkey]['nstar'] = starnum
				if starnum < drop_nstar_less_than_this:
					self.photometry_info[imgkey]['drop'] = 2
			except Exception,e:
				nstar = self.anomaly_signs['nstar']
				print "source finding or region file saving error on %s"%imgkey
				self.photometry_info[imgkey]['nstar'] = nstar
				self.photometry_info[imgkey]['drop'] = 2
		elif method == 'daofind':
			try:
				print imgkey
				starnum = self.source_detection_single_apphot_daofind(imgkey, which_dir = which_dir, threshold=daofind_nsigma_thresh, redo_daofind=redo_daofind)
				self.photometry_info[imgkey]['nstar'] = starnum
				if starnum < drop_nstar_less_than_this:
					print "%s stars detected on %s, and drop value 2 will be assigned"%(starnum, imgkey)
					self.photometry_info[imgkey]['drop'] = 2
			except Exception,e:
				print e
				nstar = self.anomaly_signs['nstar']
				print "source finding or region file saving error on %s"%imgkey
				self.photometry_info[imgkey]['nstar'] = nstar
				self.photometry_info[imgkey]['drop'] = 2
		else:
			raise IOError("Invalid input for method...")

	def __source_detection_single_apphot_daofind(self, image_key, daofind_outfile, which_dir='raw_image', fwhmpsf=None, emission=True, sigma=None, datamin=None, datamax=None, readnoise=None, epadu=None, threshold=3, nsigma=1.5,  psffwhm_hc=30):
		'''
		source detection using daofind from iraf apphot package

		'''
		imgkey= image_key.split('.')[0]
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		if fwhmpsf is None:
			fwhmpsf= self.photometry_info[image_key]['fwhm']

		if fwhmpsf <= 0 or fwhmpsf >psffwhm_hc: #check if value reasonable
			raise ValueError("unrealistic fwhm of stellar profile, please check...")

		if epadu is None:
			try:
				epadu = float(self.base.telescopes_info[self.current_telescope]['gain'])
			except:
				print "gain value of 1.0 will be used for %s"%image_key
				epadu = 1.0

		if sigma is None:
			stat_retdict = imstat_iraf(input_image,fields = "image,npix,mode, midpt,mean,stddev,min,max",lower = 'INDEF',upper = 'INDEF',nclip = 5,lsigma = 3.0,usigma = 3.0,binwidth = 1.0, verbose=0)
			sigma = stat_retdict['stddev']
		if datamin is None:
			datamin = 'INDEF'
		if datamax is None:
			datamax = 'INDEF'
		if readnoise is None:
			try:
				readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
			except:
				print "readnoise value of 10.0 will be used for %s"%image_key
				readnoise = 10.0

		#print fwhmpsf, sigma, datamax, threshold
		daofind_iraf(input_image, output=daofind_outfile, fwhmpsf=fwhmpsf, emission=emission, sigma=sigma, datamin=datamin, datamax=datamax, readnoise=readnoise, epadu=epadu, threshold=threshold, nsigma=nsigma)


	def source_detection_single_apphot_daofind(self, image_key, daofind_outfile=None, which_dir='raw_image', threshold=3, deleteINDEF=True, redo_daofind=True):
		'''
		source detection for single image
		'''
		imgkey_s = image_key.split('.')[0]
		if daofind_outfile is None:
			daofind_outfile = os.path.join(self.stars_dir, imgkey_s+'.coo')
		if not os.path.exists(daofind_outfile) or redo_daofind:
			print daofind_outfile
			self.__delete_file_if_exist(daofind_outfile)
			self.__source_detection_single_apphot_daofind(image_key, daofind_outfile, which_dir=which_dir, threshold=threshold)
		daofind_ret = Table.read(daofind_outfile, format='daophot')

		edge_xl = self.edge_xl
		edge_yl = self.edge_yl
		edge_xu = self.edge_xu
		edge_yu = self.edge_yu
		if self.photometry_info[image_key]['NX'] != 0 and edge_xu<0:
			edge_xu = self.photometry_info[image_key]['NX'] + edge_xu
		if self.photometry_info[image_key]['NY'] != 0 and edge_yu<0:
			edge_yu = self.photometry_info[image_key]['NY'] + edge_yu

		print "source detection region [%s:%s, %s:%s]"%(edge_xl, edge_xu, edge_yl, edge_yu)
		daofind_ret = daofind_ret[(daofind_ret['XCENTER']>edge_xl)*(daofind_ret['XCENTER']<edge_xu)*(daofind_ret['YCENTER']>edge_yl)*(daofind_ret['YCENTER']<edge_yu)]
		if deleteINDEF and daofind_ret.masked:  #table read with daophot results in a table with INDEF values masked
			for colname in ['MAG', 'SHARPNESS','SROUND','GROUND']:
				daofind_ret = daofind_ret[~daofind_ret[colname].mask]

		imgkey= image_key.split('.')[0]
		stars_filename = os.path.join(self.stars_dir,imgkey+self.starfile_suffix)
		self.__delete_file_if_exist(stars_filename)
		if os.path.exists(stars_filename):
			os.remove(stars_filename)
		daofind_ret.write(stars_filename, format='ascii.commented_header')

		#save ds9 region of found sources
		regfile = os.path.join(self.stars_dir,imgkey_s+self.regfile_suffix)
		stars = np.loadtxt(stars_filename)
		create_ds9_region_file(stars, regfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = 10)
		starnum_this = len(open(stars_filename).readlines())

		return starnum_this

	def __source_detection_single_fistar_fitsh(self,image_key,ntime_above_bkg=3, which_dir = 'modified_image'):
		'''
		get the source list with fistar task from fitsh package
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		stars_filename = os.path.join(self.stars_dir,image_key.split('.')[0]+self.starfile_suffix)
		self.__delete_file_if_exist(stars_filename)

		peakthres = self.photometry_info[image_key]['bkg'] + ntime_above_bkg*np.sqrt(self.photometry_info[image_key]['bkg'])
		output_psf = os.path.join(self.stars_dir,image_key.split('.')[0]+'.psf')

		fistar = os.path.join(self.base.fitsh_dir, 'fistar')
		command = "%s -i %s -o %s -t %s --model elliptic -F x,y,magnitude,bg"%(fistar, input_image,stars_filename,peakthres)
		try:
                	os.system(command)
                except:
                	raise IOError("source finding on image %s fails..."%image_key)

		#fistar source detection algorithm find false targets in the edge region of the image if the image has bright edge
		stars = np.loadtxt(stars_filename)
		if self.filter_after_fistar:
			if self.edge_xl is not None:
				xmin = self.edge_xl
				stars = stars[stars[:,0]> xmin]
			if self.edge_xu is not None:
				xmax = self.edge_xu
				stars = stars[stars[:,0]< xmax]
			if self.edge_yl is not None:
				ymin = self.edge_yl
				stars = stars[stars[:,1]> ymin]
			if self.edge_yu is not None:
				ymax = self.edge_yu
				stars = stars[stars[:,1]< ymax]
			np.savetxt(stars_filename, stars, fmt='%6.2f %6.2f %5.2f %6.2f')

		#save ds9 region of found sources
                regfile = os.path.join(self.stars_dir,image_key.split('.')[0]+self.regfile_suffix)
		stars = np.loadtxt(stars_filename)
		create_ds9_region_file(stars, regfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = 10)

                starnum_this = len(open(stars_filename).readlines())

		return starnum_this

	def __source_detection_single_sfind_diapl(self,image_key, which_dir = 'modified_image'):
		'''
		get star number through sfind and save the star list file
		'''
		sfind_par = os.path.join(self.parafile_dir,'sfind.par')
		instrument_par = os.path.join(self.parafile_dir,'instrument.par')
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)

		stars_filename = os.path.join(self.stars_dir,image_key.split('.')[0]+self.starfile_suffix)
		self.__delete_file_if_exist(stars_filename)
		sfind = os.path.join(self.base.diapl_dir, 'sfind')
		command = "%s %s %s %s %s" %(sfind, sfind_par,instrument_par,input_image,stars_filename)
		print command
		try:
			os.system(command)
			print "1"
		except:
			raise IOError("source finding on image %s fails..."%image_key)

		regfile = os.path.join(self.stars_dir,image_key.split('.')[0]+self.regfile_suffix)
		stars = np.loadtxt(stars_filename)
		create_ds9_region_file(stars, regfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = 10)

		print "2"

		starnum_this = len(open(stars_filename).readlines())

		return starnum_this


	def get_source_detections_table(self, imgkey, method='daofind'):
		'''
		Load the source detections to a table
		'''
		imgkey_s = imgkey.split('.')[0]
		sourcefile = os.path.join(self.stars_dir, imgkey_s+self.starfile_suffix)
		sources = Table.read(sourcefile, format='ascii')
		if method == 'daofind':
			renamecols_dict = {'col1': 'x', 'col2':'y', 'col3':'mag', 'col4':'sharpness', 'col5':'sround', 'col6':'ground', 'col7':'id'}
			for colname_old in renamecols_dict.keys():
				sources.rename_column(colname_old, renamecols_dict[colname_old])

		return sources

#Get reference image
	def get_template_image(self):
		'''
		Find the template image
		'''
		self.__get_template_image()

	def renew_template_image(self):
		'''
		re-find the template image
		'''
		for image_key in self.photometry_info.keys():
			self.photometry_info[image_key]['template'] = 0

		self.__get_template_image()

	def __get_template_image(self, updatetable=1):
		'''
		assign template image for each filter

		This function will not remove old template image
		refer to self.renew_template_image, if you want to assign new template image
		'''
		self.__dict2table()

		imgs_info = self.photometry_record_table
		imgs_info = self.__select_rows_from_table(imgs_info,'drop',0)
		names = imgs_info['name']
		flts = imgs_info['flt']
		flts_unique = np.unique(flts)

		print flts_unique

		self.__find_template_imagekey(updatetable=updatetable)
		flts_tpl_old = self.templates.keys()

		for flt in flts_unique:
			if flt in flts_tpl_old:
				continue

			flt_mask = flts==flt
			if np.sum(flt_mask)==1:
				image_key = names[flt_mask][0]
				self.photometry_info[image_key]['template'] = 1
			else:
				self.__get_template_image_single(flt)

	def __get_template_image_single(self,flt,verbose=1):
		'''
		if verbose==2 the template will be displayed and saved by aplpy
		'''
		self.__dict2table()
		info_sf = self.__select_rows_from_table(self.photometry_record_table,'flt',flt,)
		info_sf = self.__select_rows_from_table(info_sf,'drop',0)

		bkg     = info_sf['bkg']
		fwhm    = info_sf['fwhm']
		starnum = info_sf['nstar']

		bkg_mask = bkg != self.anomaly_signs['bkg']
		fwhm_mask = fwhm != self.anomaly_signs['fwhm']
		starnum_mask = starnum != self.anomaly_signs['nstar']
		valid_mask = bkg_mask*fwhm_mask*starnum_mask

		info_sf = info_sf[valid_mask]
		bkg     = info_sf['bkg']
		fwhm    = info_sf['fwhm']
		starnum = info_sf['nstar']
		coadd_template = 0 	# use the best single image as the template
		if coadd_template == 0:
		    num = 1		#how many best images you want for the given filter
		    inputpara = np.hstack((bkg.reshape(len(bkg),1),fwhm.reshape(len(fwhm),1),starnum.reshape(len(starnum),1)))
		    if verbose:
			print inputpara
		    ret,wantedindice=get_best(inputpara,self.template_factors,num)
		    image_key = info_sf['name'][wantedindice[0]]
		    self.photometry_info[image_key]['template'] = 1

		    if verbose == 2:
			import aplpy
			tpl_image_name = 'template_'+flt+'.eps'
			tpl_image_save = os.path.join(self.template_dir,tpl_image_name)
			self.__delete_file_if_exist(tpl_image_save)
			img = aplpy.FITSFigure(self.images[image_key])
			img.show_grayscale()
			img.save(tpl_image_save)

	def __find_template_imagekey(self, updatetable=1):
		'''
		find template images
		'''
		if updatetable:
			self.__dict2table()
		#info_notdrop = self.__select_rows_from_table(self.photometry_record_table,'drop',0)
		info_template  = self.__select_rows_from_table(self.photometry_record_table,'template',1)
		#print info_template

		for imgkey in info_template['name']:
			flt = self.photometry_info[imgkey]['flt']
			if flt in self.templates.keys():
				if self.templates[flt] != imgkey:
					raise ValueError('two template images found for the same band')
			self.templates[flt] = imgkey
			self.templates_after_astrometry[flt] = os.path.join(self.template_dir,'cal_' + imgkey)

	def assign_template_image(self,flt,image_key):
		'''
		delete old if exist and add new one assigned
		'''
		self.__dict2table()
		info_sf = self.__select_rows_from_table(self.photometry_record_table,'flt',flt,)
		names = info_sf['name']
		template_flags     = info_sf['template']

		for img in names:
			self.photometry_info[img]['template'] = 0
		if image_key not in names:
			raise ValueError("%s not in %s band"%(image_key,flt))
		self.photometry_info[image_key]['template'] = 1

#Source identification
	def get_target_xys(self,tpl_method='ds9', select_from_sources=True, which_dir = 'raw_image', updatetable=1):
		'''
		get the target position on image
		INPUTS:
			tpl_method: 'astrometry', 'ds9' or 'manual'
			which_dir: the directory containing the image to work on
		'''
		if updatetable:
			self.__dict2table()
		self.__find_template_imagekey(updatetable=updatetable)
		names = self.photometry_record_table['name']
		for flt in self.templates.keys():
			self.__get_templates_xy(flt,method = tpl_method, select_from_sources=select_from_sources, which_dir = which_dir)
		for image_key in names:
			if self.photometry_info[image_key]['drop']>0:
				continue
			if self.photometry_info[image_key]['template'] == 1:
				continue
			if self.photometry_info[image_key]['x'] == 0.0 or self.photometry_info[image_key]['y'] == 0:
				self.__get_target_xy_single(image_key, select_from_sources=select_from_sources, which_dir=which_dir)
				if self.photometry_info[image_key]['x'] == 0.0 or self.photometry_info[image_key]['y'] == 0:
					self.photometry_info[image_key]['drop'] = 5

	def __get_xys_on_image(self,image_key, xyfile=None):
		'''
		read the star file and get the xys of found stars
		return (N,2) array
		INPUTS:
			xyfile: the source file containing the x,y coordinate
		'''
		if xyfile is None:
			stars_file = image_key.split('.')[0]+self.starfile_suffix
			stars_file_abs = os.path.join(self.stars_dir,stars_file)
		else:
			stars_file_abs = xyfile

		if os.path.isfile(stars_file_abs):
			data = np.loadtxt(stars_file_abs)
			xys = data[:,0:2]
		else:
			print "The star file %s does not exist... you can perform sources finding through 'sfind' or 'daofind'"%stars_file_abs
			xys = None

		return xys

	def __get_templates_xy(self,flt,method = 'astrometry', which_dir='raw_image', select_from_sources=True, updatetable=1):
		'''
		Input:
			method: 'astrometry' or 'ds9' or 'manual'
		'''

		print "working on %s band template"%flt
		self.__find_template_imagekey(updatetable=updatetable)
		tpl_imgkey = self.templates[flt]
		tpl_image_cal = self.templates_after_astrometry[flt]

		if self.photometry_info[tpl_imgkey]['x'] == 0.0 or self.photometry_info[tpl_imgkey]['y'] == 0 or self.renew_template_xys:
			if method == 'astrometry':
				xy = self.__get_tpl_xy_on_image_from_astrometry(flt, which_dir=which_dir)
			elif method == 'ds9':
				starfile = os.path.join(self.stars_dir, tpl_imgkey.split('.')[0]+self.starfile_suffix)
				if os.path.exists(starfile) and select_from_sources:
					starxys = np.loadtxt(starfile)
					regionfile = os.path.join(self.stars_dir, tpl_imgkey.split('.')[0]+"_temp.reg")
					circle_radius = self.criteria_point_match_general
					create_ds9_region_file(starxys, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = circle_radius, color='red', width=2, load_text=False)
					input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)
					xy = self.__get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=True)
					self.__delete_file_if_exist(regionfile)					
					xy_target = [xy[0],xy[1]]
					match_criteria = self.criteria_point_match_general
					yesfound, xy_match,index = self.__find_corresponding(starxys[:,[0,1]], xy_target, match_criteria)
					if yesfound and len(index)==1:
						print "the corresponding object #%s found at (%s)"%(index[0], xy_match)
						matchfound = starxys[index[0]]
						xy = [matchfound[0], matchfound[1]]
					else:
						print "No exclusive match found for input ds9 cursor coordinate (x,y) in the input source list within circle radius %s"%match_criteria 
						print "the ds9 cursor coordinate (%s,%s) will be used"%(xy[0], xy[1])
				else:
					xy = self.get_xy_on_image_from_ds9(image_key = tpl_imgkey, which_dir = which_dir)
			elif method == 'manual':
				x = raw_input("Enter the x position of the current supernova on %s band template image:"%flt)
				y = raw_input("Enter the y position of the current supernova on %s band template image:"%flt)
				xy = [float(x),float(y)]
			else:
				print "Invalid input for method"
				raise IOError("Sorry...")

			x_image = round(xy[0],1)
			y_image = round(xy[1],1)
			self.photometry_info[tpl_imgkey]['x'] = x_image
			self.photometry_info[tpl_imgkey]['y'] = y_image

	def __get_target_xy_single(self,image_key,depend_on_tpl = True, non_tpl_method = 'ds9', select_from_sources=True, which_dir = 'raw_image'):
		'''
		Get the supernova position on images
		Input:
			depend_on_tpl: if True, the supernova position on the template image required and working-on image get its position relative the template image
					if False, the working on image gets its position independently
			non_tpl_method: if depend_on_tpl is False, then individual image get sn position with method assigned by this input
					'ds9', pick up on ds9 display
					'manual', input after prompt
					'wcs', convert from self.sn_ra_world_deg and self.sn_dec_world_deg with image wcs info
		'''
		print image_key
		imgkey_s = image_key.split('.')[0]
		if depend_on_tpl:
			flt = self.photometry_info[image_key]['flt']
			tpl_imgkey = self.templates[flt]
			tplimgkey_s = tpl_imgkey.split('.')[0]
			ref_list   = os.path.join(self.stars_dir,tplimgkey_s+self.starfile_suffix)
			input_list = os.path.join(self.stars_dir,imgkey_s+self.starfile_suffix)
			match_output = os.path.join(self.stars_dir,imgkey_s+ '_tpl.match')
			self.__delete_file_if_exist(match_output)
			trans_fitting = os.path.join(self.stars_dir,imgkey_s + '_tpl.coef')
			self.__delete_file_if_exist(trans_fitting)

			self.fitsh_grmatch_type = 'point'
			matched_table = self.__fitsh_grmatch(ref_list, input_list, match_output, trans_output=trans_fitting)
			if matched_table is None:
				self.photometry_info[image_key]['drop'] = 5
				x_image = 0
				y_image = 0
			else:
				trans_output = os.path.join(self.stars_dir,imgkey_s+'_tpl.trans')
				self.__delete_file_if_exist(trans_output)
				self.__fitsh_grtrans(match_output,trans_output,trans_fitting)
				porder,dxfit,dyfit = self.__extract_transformation_fitting(trans_fitting)
				tpl_xy = [self.photometry_info[tpl_imgkey]['x'],self.photometry_info[tpl_imgkey]['y']]
				target_xy = self.__transform_xys(dxfit,dyfit, porder, tpl_xy)
				print target_xy
				x_image = round(target_xy[0],1)
				y_image = round(target_xy[1],1)

				if x_image < 0 or y_image <0:
					self.photometry_info[image_key]['drop'] = 5
		else:
			if non_tpl_method == 'ds9':
				starfile = os.path.join(self.stars_dir, imgkey_s+self.starfile_suffix)
				if os.path.exists(starfile) and select_from_sources:
					starxys = np.loadtxt(starfile)
					regionfile = os.path.join(self.stars_dir, imgkey_s+"_temp.reg")
					circle_radius = self.criteria_point_match_general
					create_ds9_region_file(starxys, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = circle_radius, color='red', width=2, load_text=False)
					input_image = self.__get_internal_image(image_key, which_dir=which_dir)
					xy = self.__get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=True)
					self.__delete_file_if_exist(regionfile)					
					if xy is None:
						print "no input coordinate (x,y) for the target; drop this"
						self.photometry_info[image_key]['drop'] = 99
						return 
					xy_target = [xy[0],xy[1]]
					match_criteria = self.criteria_point_match_general
					yesfound, xy_match,index = self.__find_corresponding(starxys[:,[0,1]], xy_target, match_criteria)
					if yesfound and len(index)==1:
						print "the corresponding object #%s found at (%s)"%(index[0], xy_match)
						matchfound = starxys[index[0]]
						xy = [matchfound[0], matchfound[1]]
					else:
						print "No exclusive match found for input ds9 cursor coordinate (x,y) in the input source list within circle radius %s"%match_criteria 
						print "the ds9 cursor coordinate (%s,%s) will be used"%(xy[0], xy[1])
				else:
					xy = self.get_xy_on_image_from_ds9(image_key = image_key, regionfile=None, which_dir=which_dir)
			elif non_tpl_method == 'manual':
				x = raw_input("Enter the x position of the current supernova on image %s:"%image_key)
				y = raw_input("Enter the y position of the current supernova on image %s:"%image_key)
				xy = [float(x),float(y)]
			elif non_tpl_method == 'wcs':
				inputimage = self.__get_internal_image(image_key, which_dir=which_dir)
				xy = self.__get_xy_on_image_from_wcsinfo(inputimage)
			else:
				print "Invalid input for method"
				raise IOError("Sorry...")

			x_image = round(xy[0],1)
			y_image = round(xy[1],1)

		self.photometry_info[image_key]['x'] = x_image
		self.photometry_info[image_key]['y'] = y_image

	def __get_tpl_xy_on_image_from_astrometry(self, flt, verbose=1):
		'''
		Output:
			image_xy: array
		'''
		self.__find_template_imagekey()
		inputimage = self.templates_after_astrometry[flt]
		tplxy = self.__get_xy_on_image_from_wcsinfo(inputimage)
		return tplxy

	def __get_xy_on_image_from_wcsinfo(self, inputimage, verbose=0):
		'''
		world2image conversion: wcs info in fits header required
		'''
		if self.sn_ra_world_deg is None or self.sn_dec_world_deg is None:
			if self.current_sn in self.base.registerd_photometry_objects_info.keys():
				self.sn_ra_world_deg  = float(self.base.registerd_photometry_objects_info[self.current_sn]['ra'])
				self.sn_dec_world_deg = float(self.base.registerd_photometry_objects_info[self.current_sn]['dec'])

		if self.sn_ra_world_deg is None or self.sn_dec_world_deg is None:
			sn_ra = raw_input("Enter the RA of the current supernova:")
			sn_dec = raw_input("Enter the Dec of the current supernova:")
			self.sn_ra_world_deg = float(sn_ra)
			self.sn_dec_world_deg = float(sn_dec)

		world_radec = np.array([self.sn_ra_world_deg,self.sn_dec_world_deg])
		world_radec = world_radec.reshape(1,2)
		image_xy = self.__world2image(inputimage,world_radec)
		if verbose>0:
			print image_xy

		if verbose>1:
			image = self.images[self.templates[flt]]
			self.label_single_xy_region(image,image_xy)

		return image_xy[0]

	def get_xy_on_image_from_ds9(self,image_key=None,input_image=None, which_dir='raw_image', regionfile =None, newds9=True):
		'''
		pick up xy coordinate from ds9 image
		'''

		if input_image is None:
			if image_key is not None:
				input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			else:
				raise IOError("No valid input...")
		else:
			if os.path.exists(input_image):
				print input_image
			else:
				raise IOError("input image doesn't exist!")

		xy = self.__get_xy_on_image_from_ds9(input_image, regionfile=regionfile, newds9=newds9)

		return xy



#Aperture photometry
	def get_apphot_iraf_parameters(self, image_key, saveparfile = True):
		'''
		Get the default config parameters for apphot.daophot
		'''
		imgkey = image_key.split('.')[0]
		parfile = os.path.join(self.aperture_photometry_dir, imgkey+'.par')

		config = ConfigParser.ConfigParser()
		if os.path.exists(parfile) and (not self.renew_apphot_parfile):
			with open(parfile,'rw') as cfgfile:
				config.readfp(cfgfile)
				app     = config.getfloat('param','app')
				skyin   = config.getfloat('param','skyin')
				skywidth  = config.getfloat('param','skywidth')
				fwhmpsf = config.getfloat('param','fwhmpsf')
				sky_sigma_down  = config.getfloat('param','sky_sigma_down')
				sky_sigma_up    = config.getfloat('param','sky_sigma_up')
				def_zeropt      = config.getfloat('param','def_zeropt')

				self.apphot_iraf_options['fwhmpsf']       = fwhmpsf
				self.apphot_iraf_options['app']           = app
				self.apphot_iraf_options['skyin']         = skyin
				self.apphot_iraf_options['skywidth']      = skywidth
				self.apphot_iraf_options['sky_sigma_down']= sky_sigma_down
				self.apphot_iraf_options['sky_sigma_up']  = sky_sigma_up
				self.apphot_iraf_options['def_zeropt']    = def_zeropt
		else:
			fwhmpsf = self.photometry_info[image_key]['fwhm']
			self.apphot_iraf_options['fwhmpsf'] = fwhmpsf

			if self.apphot_iraf_autoscale_aperture:
				self.__autoscale_apphot_iraf_aperture_pars(fwhmpsf)

			if saveparfile:
				self.__delete_file_if_exist(parfile)
				config.add_section('param')
				config.set('param', 'fwhmpsf', self.apphot_iraf_options['fwhmpsf'])
				config.set('param', 'app', self.apphot_iraf_options['app'])
				config.set('param', 'skyin', self.apphot_iraf_options['skyin'])
				config.set('param', 'skywidth', self.apphot_iraf_options['skywidth'])
				config.set('param', 'sky_sigma_down', self.apphot_iraf_options['sky_sigma_down'])
				config.set('param', 'sky_sigma_up', self.apphot_iraf_options['sky_sigma_up'])
				config.set('param', 'def_zeropt', self.apphot_iraf_options['def_zeropt'])

				fid = open(parfile,'w')
				config.write(fid)
				fid.close()

	def __autoscale_apphot_iraf_aperture_pars(self, fwhmpsf):
		'''
		Auto-scale the aperture sizes (source aperture and background region)
		'''
		if fwhmpsf <=0 or fwhmpsf>50: #check to prevent radicuous results
			raise ValueError("The star profile with negative fwhm or fwhm value larger than 50 pixels. Really??")

		self.apphot_iraf_options['app']     = self.apphot_iraf_app_nfwhm * fwhmpsf
		self.apphot_iraf_options['skyin']   = self.apphot_iraf_skyin_nfwhm * fwhmpsf
		self.apphot_iraf_options['skywidth']= self.apphot_iraf_skywidth_nfwhm * fwhmpsf

	def aperture_photometry_apphot_iraf(self,flts = None, aperture_size_fixed = False, \
	centering_target = True, centering_ref_stars = True, updatetable=1, which_dir = 'raw_image'):
		'''
		Aperture photometry on images with filter given in 'flts', which default is all filters for
		which the template image is given.

		Two functions involved in for this purpose:
		self.aperture_photometry_apphot_iraf_flt(flt): aperture photometry for all given objects in the image
		self.aperture_photometry_apphot_iraf_target_flt(flt): aperture photometry for target


		Inputs:
			flts: default None; allowed input for example ['B','V','R','I']
			which_dir: perform aperture photometry on images stored in self.raw_image_dir or self.modified_image_dir, indicated by
			'raw_image' or 'modified_image'
		'''
		if updatetable:
			self.__dict2table()
		if len(self.templates) ==0:
			self.__find_template_imagekey(updatetable=updatetable)

		if flts is None:
			flts = self.templates.keys()

		for flt in flts:
			print "Working on %s band"%flt
			self.aperture_photometry_apphot_iraf_flt(flt, which_dir = which_dir, aperture_size_fixed=aperture_size_fixed, centering = centering_ref_stars)
			self.aperture_photometry_apphot_iraf_target_flt(flt, which_dir = which_dir, aperture_size_fixed=aperture_size_fixed, centering = centering_target, updatetable=updatetable)

	def aperture_photometry_apphot_iraf_target_flt(self,flt, which_dir = 'raw_image', \
	aperture_size_fixed = False, imagefwhm=None, centering=True, x=None, y=None, updatetable=1):
		'''
		pipeline script for aperture photometry on images with filter 'flt'
		which_dir: 'raw_image', 'modified_image', or 'subtracted_image'
		See also:
			self.aperture_photometry_apphot_iraf_target_single
		'''
		if updatetable:
			self.__dict2table()
		info_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		images = info_flt['name']

		for image_key in images:
			print image_key
			if self.photometry_info[image_key]['drop'] >0:
				continue
			if self.photometry_info[image_key]['instmag'] != 99.99 and not self.renew_aperture_photometry:
				print "already done"
				continue
			photret_singlept = self.aperture_photometry_apphot_iraf_target_single(image_key, which_dir =  which_dir, fwhm_img=imagefwhm, centering = centering, x=x, y=y)

			if len(photret_singlept) == 0:
				self.photometry_info[image_key]['drop'] = 4
				continue

			mag = photret_singlept[0,2]
			magerr = photret_singlept[0,3]

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr

	def aperture_photometry_apphot_iraf_target_single(self, image_key, x=None, y=None, \
	fwhm_img=None, ap = None, skyin=None, skywidth=None, centering= True, which_dir = 'raw_image', verbose=0):
		'''
		For each image: get target xy; get apphot input parameters; perform aperture photometry
		Involved functions:
			self.get_apphot_iraf_parameters
			self.aperture_photometry_apphot_iraf_single_image_xys_given

		INPUTS:
			image_key:
			x:
			y:
			fwhm_img: if given, then app, skyin and skywidth are based on this
			ap, skyin, skywidth: if given, then overwrite the autoscaled aperture radiu, skyin and skywidth

		'''
		print image_key
		if x is None:
			x = self.photometry_info[image_key]['x']
		if y is None:
			y = self.photometry_info[image_key]['y']

		if x == self.photometry_info_init_values['x'] or y==self.photometry_info_init_values['y']:
			raise ValueError("no valid souce position is given...")

		xy = np.array([x,y]).reshape((1,2))
		image = self.__get_internal_image(image_key, which_dir=which_dir)
		self.get_apphot_iraf_parameters(image_key, saveparfile=True)
		options = self.apphot_iraf_options.copy()
		if fwhm_img is not None:
			options['app'] = fwhm_img*self.apphot_iraf_app_nfwhm
			options['skyin'] = fwhm_img*self.apphot_iraf_skyin_nfwhm
			options['skywidth']= fwhm_img*self.apphot_iraf_skywidth_nfwhm
		if ap is not None:
			options['app'] = ap
		if skyin is not None:
			options['skyin'] = skyin
		if skywidth is not None:
			options['skywidth'] = skywidth
		output = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0]+'_%s_%s.mag'%(str(x),str(y)))
		photret_singlept = self.aperture_photometry_apphot_iraf_single_image_xys_given(image, xy, options, output=output, centering = centering)

		if verbose:
			daophotret = Table.read(output, format='daophot')
			print daophotret[daophotret.colnames[0:11]]
			print daophotret[daophotret.colnames[11:22]]
			print daophotret[daophotret.colnames[22:33]]

		return photret_singlept


	def __aperture_photometry_single_target_ds9_pick(self, image_key, match_criteria=3, which_dir='raw_image',  verbose=1 , update_xy=True, newds9=True, box_length=100, box_xc=None,box_yc=None):
		'''
		this is created for the initial purpose of SMARTS-NIR photometry where the source can be very sparse on the image. The PSF photometry is obtained first then select the target mag from the
		magnitude list.
		'''
		imgkey_s = image_key.split('.')[0]
		apphot_ret_file = os.path.join(self.aperture_photometry_dir,imgkey_s+self.apphot_ret_file_suffix)
		if not os.path.exists(apphot_ret_file):
			print "aperture photometry result file %s not exists"%apphot_ret_file
			self.photometry_info[image_key]['drop'] = 4
		else:

			apphot_ret = np.loadtxt(apphot_ret_file)
			xys = apphot_ret[:,0:2]
			mags = apphot_ret[:,2]
			magerrs = apphot_ret[:,3]

			input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			regionfile = os.path.join(self.aperture_photometry_dir, imgkey_s+"_xymag.reg")

			if (box_xc is not None) and (box_yc is not None):
				tmpmask = (apphot_ret[:,0]>(box_xc-box_length/2))*(apphot_ret[:,0]<(box_xc+box_length/2))*(apphot_ret[:,1]<(box_yc+box_length/2))*(apphot_ret[:,1]>(box_yc-box_length/2))
				apphot_reg = apphot_ret[tmpmask,:]
			else:
				apphot_reg = apphot_ret

			create_ds9_region_file(apphot_reg, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

			xy = self.get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=newds9)
			print "xy from mouse pick up:", xy
			x = xy[0]
			y = xy[1]
			xy_target = [x,y]
			yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, match_criteria)
			if (not yesfound):
				print "no object found within criteeia..."
				mag = 99.99
				magerr = 99.99
				sigdrop = 4
			elif len(index)>1:
				print 'more than one objects fall into the criteria region of reference star. This object is droped!'
				mag = 99.99
				magerr = 99.99
				sigdrop = 4
			else:
				print "the corresponding object found at (%s)"%xy_match
				if update_xy:
					snx = xy_match[0][0]
					sny = xy_match[0][1]
					print "the (x,y) coordinate will be updated to (%s,%s)"%(snx, sny)
					self.photometry_info[image_key]['x'] = snx
					self.photometry_info[image_key]['y'] = sny
				mag = apphot_ret[index,2]
				magerr = apphot_ret[index,3]
				sigdrop = 0
				if verbose:
					print "the mag: %s+/-%s"%(mag, magerr)

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr
			self.photometry_info[image_key]['drop'] = sigdrop

	def aperture_photometry_apphot_iraf_single_image(self, image_key, source_infile=None, phot_outfile=None, output_dir=None, overwrite=False, which_dir='modified_image', aperture_size_fixed=False, centering=True):
		'''
		Aperture photometry for individual image. This is not intended to be used in patch mode (not for pipeline usage)
		INPUTS:	
			source_infile: the input file containing the x,y coordites of the sources; requirement: readable by np.loadtxt and the first two columns are the x,y coordinate
			phot_outfile: the photometry output file
			output_dir: the default is to store the photometry file xxx.apphot in the self.aperture_photometry_dir; otherwise to the folder provided here
		'''

		imgkey_s = image_key.split('.')[0]
		photret_file = image_key.split('.')[0]  + self.apphot_ret_file_suffix
		if output_dir is None:
			output_dir = self.aperture_photometry_dir

		photret_file_abs = os.path.join(output_dir,photret_file)

		if os.path.exists(photret_file_abs):
			if not overwrite:
				print "photometry result file for %s exists: %s; please specify overwrite=True if you want to perform aperture photometry again"%(image_key, photret_file_abs)
				return

		xys = self.__get_xys_on_image(image_key, xyfile=source_infile)
		image = self.__get_internal_image(image_key, which_dir=which_dir)
		self.get_apphot_iraf_parameters(image_key, saveparfile=True)
		options = self.apphot_iraf_options.copy()
		if phot_outfile is None:
			photfile = os.path.join(self.aperture_photometry_dir, imgkey_s+'.mag')
		else:
			photfile = phot_outfile
		self.__delete_file_if_exist(photfile)
		photret = self.aperture_photometry_apphot_iraf_single_image_xys_given(image, xys, options, output=photfile, centering=centering)
		self.__delete_file_if_exist(photret_file_abs)
		np.savetxt(photret_file_abs,photret,fmt="%6.2f %6.2f %6.3f %6.3f")

		phot_source_reg =  os.path.join(self.aperture_photometry_dir, imgkey_s + '_apphot.reg')
		create_ds9_region_file(photret, phot_source_reg, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)


	def aperture_photometry_apphot_iraf_flt(self,flt, which_dir = 'raw_image', aperture_size_fixed = False, centering = True):
		'''
		pipeline script for aperture photometry on images with filter 'flt'
		For each image: get objects xys; get apphot input parameters; perform aperture photometry

		Involved functions:
			self.aperture_photometry_apphot_iraf_single_image
		'''
		info_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		images = info_flt['name']

		for image_key in images:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			imgkey_s = image_key.split('.')[0]
			photret_file = image_key.split('.')[0]  + self.apphot_ret_file_suffix
			photret_file_abs = os.path.join(self.aperture_photometry_dir,photret_file)
			if os.path.exists(photret_file_abs):
				if not self.renew_aperture_photometry_retfile:
					print "!!! Aperture photometry for %s %s already exists. Set self.renew_aperture_photometry_retfile to True to overwrite."%(image_key, photret_file_abs)
					continue
			self.aperture_photometry_apphot_iraf_single_image(image_key, overwrite=True, which_dir=which_dir, aperture_size_fixed=aperture_size_fixed, centering=centering)



	def aperture_photometry_apphot_iraf_single_point(self,image,xy,options=None, imagefwhm=None, centering = True, aperture_size_fixed=1):
		'''
		The same as self.aperture_photometry_apphot_iraf_single_image_xys_given but no region file will be saved
		Please see self.aperture_photometry_apphot_iraf_single_image_xys_given for details

		INPUTS:
			image: the input image not image key
			xy:
			options: iraf appot parameter dictory: app, def_zeropoint, fwhmpsf, sky_sigma_down, sky_sigma_up, skyin, skywidth
			imagefwhm: the fwhm of the image, if provided this will overwrite the one in options
			centering: if True, then centering algorithm in self.apphot_iraf_calgorithm will be used
		'''
		if options is None:
			options = self.apphot_iraf_options.copy()

		if imagefwhm is not None:
			options['fwhmpsf'] = imagefwhm
		else:
			print "Attention!!!: default fwhm value of %s will be used"%options['fwhmpsf']

		if not aperture_size_fixed:
			fwhm_img = options['fwhmpsf']
			if fwhm_img <0 or fwhm_img>50:
				raise ValueError("The star profile with negative fwhm or fwhm value larger than 50 pixels. Really??")
			options['app'] = 2*fwhm_img
			options['skyin'] = 3*fwhm_img
			options['skywidth']= 3*fwhm_img

		ret = self.aperture_photometry_apphot_iraf_single_image_xys_given(image,  xy, options, output=None, centering = centering, source_reg_filename=None, delete_ERR=True)
		print ret

		return ret

	def __aperture_photometry_apphot_iraf_single_image(self, image, coo_list, output, options, centering = True, which_dir='raw_image'):
		'''
		underlying aperture photometry task which utilize 'apphot_iraf'.
		Inputs:
			image:
			options: dict containing input parameters for aperture photometry
				{'app': xx,
				 'def_zeropt': xx,
				 'fwhmpsf': xx,
				 'sky_sigma_down': xx,
				 'sky_sigma_up': xx,
				 'skyin': xx,
				 'skywidth': xx}

			centroid_option: required by apphot_iraf function; default 'centroid'

		See photometry_collection.apphot_iraf for details of aperture photometry task
		'''
		#image = self.__get_internal_image(image_key, which_dir=which_dir)

		datamin = self.apphot_iraf_datamin
		datamax = self.apphot_iraf_datamax

		if centering:
			calgorithm = self.apphot_iraf_calgorithm
		else:
			calgorithm = 'none'

		epadu = float(self.base.telescopes_info[self.current_telescope]['gain'])
		readnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])

		if self.current_ccd_readnoise is not None:
			readnoise = self.current_ccd_readnoise
		if self.current_ccd_epadu is not None:
			epadu = self.current_ccd_epadu

		app = options['app']
		zeropt = options['def_zeropt']
		fwhmpsf= options['fwhmpsf']
		sky_sigma_down = options['sky_sigma_down']
		sky_sigma_up = options['sky_sigma_up']
		skyin = options['skyin']
		wsky = options['skywidth']

		hdu = fits.open(image)
		if len(hdu)>1:
			image = image+'[0]'
		hdu.close()

		daophot_phot_iraf(image, centroid_algorithm=calgorithm, coo_list=coo_list, output=output, fwhmpsf=fwhmpsf, sigma=10, app=app, skyin=skyin, wsky=wsky, sky_sigma_down=sky_sigma_down, sky_sigma_up=sky_sigma_up, readnoise=readnoise, epadu=epadu, datamin = datamin, datamax = datamax, emission='yes', zeropt=zeropt)


	def aperture_photometry_apphot_iraf_single_image_xys_given(self, image,  xys, options, output=None, centering = True,source_reg_filename=None, delete_ERR=True):
		'''
		INPUTS:
			xys:	the objects position on the images at which the aperture photometry will take place
		'''

		coo_list =  'temp_delete_after.coo'
		f = open(coo_list, 'w')
		if np.array(xys).shape == (2,): #find if there is only one set of coords to measure
			f.write(str(xys[0]) + ' ' + str(xys[1]) + '\n')
		else:
			for this_star in range(len(xys)):
				f.write(str(xys[this_star][0]) + ' ' + str(xys[this_star][1]) + '\n')
		f.flush()
		f.close

		if output is None:
			output = 'temp_delete_after.mag'
		self.__aperture_photometry_apphot_iraf_single_image(image, coo_list, output, options, centering = centering)

		photret = Table.read(output, format='daophot')
		if delete_ERR:
			photret = photret[(photret['CIER']==0).data*(photret['SIER']==0).data*(photret['PIER']==0).data]  #get rid of results with error in centering, skyfiting and photometry

		photret = photret[['XCENTER','YCENTER','MAG','MERR']]
		if photret.masked:
			photret = photret[~photret['MAG'].mask]
			photret = photret[~photret['MERR'].mask]
		photret = np.array([photret['XCENTER'].data,photret['YCENTER'].data,photret['MAG'].data,photret['MERR'].data]).transpose()

		return photret


#PSF photometry DoPHOT
	def psf_photometry_dophot(self, which_dir = 'raw_image'):
		'''
		See self.psf_photometry_dophot_flt
		'''
		self.__find_template_imagekey()
		for flt in self.templates.keys():
			self.psf_photometry_dophot_flt(flt, which_dir = which_dir)
			self.psf_photometry_dophot_flt_target(flt)

	def psf_photometry_dophot_flt(self,flt, which_dir = 'raw_image'):
		'''
		See self.psf_photometry_dophot_single
		'''
		self.__dict2table()
		flt_table = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		for image_key in flt_table['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			psfret_savefile = self.__get_internal_psfphot_mag_file(image_key)
			if not os.path.isfile(psfret_savefile) or  self.renew_psf_photometry_retfile:
				self.psf_photometry_dophot_single(image_key,psfret_savefile,output_residuals=True, which_dir=which_dir)

	def psf_photometry_dophot_flt_target(self,flt):
		'''
		See self.psf_photometry_dophot_single_target
		'''
		self.__dict2table()
		flt_table = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		for image_key in flt_table['name']:
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if self.photometry_info[image_key]['instmag'] != 99.99 and not self.renew_psf_photometry:
				continue
			match_criteria = self.criteria_psfphot_locate_target
			self.psf_photometry_dophot_single_target(image_key,match_criteria = match_criteria)

	def psf_photometry_dophot_single_target_ds9_pick(self, image_key, match_criteria=3, which_dir='raw_image',  verbose=1 , update_xy=True, newds9=True, box_length=100, box_xc=None,box_yc=None):
		'''
		this is created for the initial purpose of SMARTS-NIR photometry where the source can be very sparse on the image. The PSF photometry is obtained first then select the target mag from the
		magnitude list.
		'''
		imgkey_s = image_key.split('.')[0]
		psfphot_ret_file = self.__get_internal_psfphot_mag_file(image_key)
		if not os.path.exists(psfphot_ret_file):
			print "PSF photometry result file %s not exists"%psfphot_ret_file
			self.photometry_info[image_key]['drop'] = 4
		else:
			psfphot_ret = np.loadtxt(psfphot_ret_file)
			xys = psfphot_ret[:,0:2]
			mags = psfphot_ret[:,2]
			magerrs = psfphot_ret[:,3]

			input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			regionfile = os.path.join(os.path.dirname(psfphot_ret_file), imgkey_s+"_xymag.reg")

			if (box_xc is not None) and (box_yc is not None):
				tmpmask = (psfphot_ret[:,0]>(box_xc-box_length/2))*(psfphot_ret[:,0]<(box_xc+box_length/2))*(psfphot_ret[:,1]<(box_yc+box_length/2))*(psfphot_ret[:,1]>(box_yc-box_length/2))
				psfphot_reg = psfphot_ret[tmpmask,:]
			else:
				psfphot_reg = psfphot_ret

			create_ds9_region_file(psfphot_reg, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

			xy = self.get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=newds9)
			print "xy from mouse pick up:", xy

			x = xy[0]
			y = xy[1]
			xy_target = [x,y]
			yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, match_criteria)

			if (not yesfound):
				print "no object found within criteeia..."
				mag = 99.99
				magerr = 99.99
				sigdrop = 4
			elif len(index)>1:
				print 'more than one objects fall into the criteria region of reference star. This object is droped!'
				mag = 99.99
				magerr = 99.99
				sigdrop = 4
			else:
				print "the corresponding object found at (%s)"%xy_match
				if update_xy:
					snx = xy_match[0][0]
					sny = xy_match[0][1]
					print "the (x,y) coordinate will be updated to (%s,%s)"%(snx, sny)
					self.photometry_info[image_key]['x'] = snx
					self.photometry_info[image_key]['y'] = sny

				mag = psfphot_ret[index,2]
				magerr = psfphot_ret[index,3]
				sigdrop = 0
				if verbose:
					print "the mag: %s+/-%s"%(mag, magerr)

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr
			self.photometry_info[image_key]['drop'] = sigdrop

	def psf_photometry_dophot_single_target(self,image_key,ds9_display = False,match_criteria = 3, drop_nstar_less_than_this = 3):
		'''
		get the mag for the target from mag list obtained from PSF photometry
		'''

		psfphot_ret_file = self.__get_internal_psfphot_mag_file(image_key)
		psfphot_ret = np.loadtxt(psfphot_ret_file)

		if len(psfphot_ret) < drop_nstar_less_than_this:
			self.photometry_info[image_key]['drop'] = 8
			return

		xys = psfphot_ret[:,0:2]
		mags = psfphot_ret[:,2]
		magerrs = psfphot_ret[:,3]
		x = self.photometry_info[image_key]['x']
		y = self.photometry_info[image_key]['y']
		xy = [x-1, y-1]

		exist,xy_match,indice_match = self.__find_corresponding(xys, xy, match_criteria)
		if exist:
			if len(indice_match)>1:
				print "more than one source detected at the SN position... drop this one"
				self.photometry_info[image_key]['drop'] = 9
			else:
				self.photometry_info[image_key]['instmag'] = mags[indice_match]
				self.photometry_info[image_key]['instmagerr'] = magerrs[indice_match]

			if ds9_display:
				input_image = self.images[image_key]
				d = self.__display_image_with_ds9(input_image, newds9=True)
				print "Red circle indicate the corresponding position from PSF photometry results"
				if len(indice_match) >1:
					for xy_reg in xy_match:
						self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='red', width=1)
				else:
					xy_reg = xy_match.ravel()
					self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='red', width=1)
				print "Green circle indicate the SN position previously provided"
				xy_reg = [x,y]
				self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='green',width=1)
		else:
			if ds9_display:
				input_image = self.images[image_key]
				d = self.__display_image_with_ds9(input_image, newds9=True)
				print "Green circle indicate the SN position previously provided"
				xy_reg = [x,y]
				self.__load_source_regions_to_ds9(d, xy_reg, radius=5, color='green',width=1)
			print "photometry result not found for %s at position %s"%(image_key,xy)
			self.photometry_info[image_key]['drop'] = 4

	def psf_photometry_dophot_single(self,image_key,ret_save_filename=None, \
	output_residuals=False, which_dir ='raw_image', pmfile=None, renew_pmfile=True, verbose=1):
		'''
		perform PSF photometry with DoPhot package

		INPUTS:
			image_key:
			ret_save_filename: the output photometry file containing the well formated x,y,mag,magerr
			output_residuals: save the residual image after model subtraction?
			which_dir: where the input image is from
			pmfile: the input parameter file; None if not specified
			renew_pmfile: you can provide explicitly the parameter file under /path/to/sndir/warehouse/psfphotometry/xxx.pm, or the parameter file can be from last run. The above existing parameter file will be replaced if renew_pmfile is True

		'''
		if verbose:
			print image_key
		img_key_s  = image_key.split('.')[0]

		if pmfile is not None: #if parameter file is specified through input variable, the pmfile will be used
			if not os.path.exists(pmfile):
				raise ValueError("the input parameter file %s not exists"%pmfile)
			pmfile_external = True
		else:
			pmfile_tmp = img_key_s + '.pm'
			pmfile_external = False
			cwd = os.getcwd()
			cwdod = os.path.join(cwd, self.current_sn)#oject dir under cwd
			cwdodid = os.path.join(cwdod, self.current_telescope)#instrument dir under cwdod
			if self.insulate_pmfile:
				if not os.path.exists(cwdod):
					os.mkdir(cwdod)
				if not os.path.exists(cwdodid):
					os.mkdir(cwdodid)
				pmfile  = os.path.join(self.current_sn, os.path.join(self.current_telescope, pmfile_tmp))
				photimg = os.path.join(self.current_sn, os.path.join(self.current_telescope, image_key))
			else:
				pmfile = pmfile_tmp
				photimg = image_key

			pmfile_src = os.path.join(cwd,pmfile)  #this is the pm file which will be the direct input of dophot command execution
			if self.dophot_version == 'fortran':
				psfphotdir = self.psfphot_dophot_fortran_dir
			else:
				psfphotdir = self.psfphot_dophot_c_dir
			pmfile_dst = os.path.join(psfphotdir, pmfile_tmp) #this is the archived pm file for specified image stored under event specific directory

			if os.path.exists(pmfile_dst) and (not renew_pmfile):
				self.__delete_file_if_exist(pmfile_src)
				shutil.copy(pmfile_dst, pmfile_src)
				savepmfile = False #need to save/archive the newly created pmfile
			else:
				savepmfile = True
				if self.dophot_version == 'fortran':
					dophot_fortran_pm = self.__prepare_pm_data_dophot_fortran(image_key, photimg)
					self.__prepare_pm_file_dophot_fortran(pmfile, **dophot_fortran_pm)
				elif self.dophot_version == 'C':
					dophot_C_pm = self.__prepare_pm_data_dophot_C(image_key, photimg)
					self.__prepare_pm_file_dophot_C(pmfile, **dophot_C_pm)
				else:
					raise IOError("dophot version %s is not supported..."%self.dophot_version)
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		input_image_here = os.path.join(os.getcwd(), photimg)
		self.__delete_file_if_exist(input_image_here)

		if self.dophot_image_prepare == 'softlink':
			ln_command = "ln -s %s %s"%(input_image,input_image_here)
			if verbose:
				print ln_command
			os.system(ln_command)
		elif self.dophot_image_prepare == 'imcopy':
			imcopy_iraf(input_image+'[0]', input_image_here)
			if verbose:
				print "imcopy %s[0] %s"%(input_image, input_image_here)

		if self.dophot_version == 'fortran':
			dophot = os.path.join(self.base.dophot_fortran_dir, 'dophot')
			storagedir = self.psfphot_dophot_fortran_dir
		elif self.dophot_version == 'C':
			dophot = os.path.join(self.base.dophot_C_dir, 'dophot')
			storagedir = self.psfphot_dophot_c_dir
		else:
			raise IOError("dophot version %s is not supported..."%self.dophot_version)

		if verbose:
			psfphot_command = '%s %s'%(dophot,pmfile)
		else:
			psfphot_command = '%s %s >/dev/null'%(dophot,pmfile)

		print psfphot_command
		os.system(psfphot_command)

		if not pmfile_external:
			if savepmfile:
				shutil.move(pmfile_src,pmfile_dst)
		self.__clean_dophot_output_files(image_key, storagedir)
		os.remove(input_image_here)

		ret_good_small = self.__extract_dophot_phot_results(image_key)
		if ret_save_filename is None:
			ret_save_filename = self.__get_internal_psfphot_mag_file(image_key)
		np.savetxt(ret_save_filename,ret_good_small)

	def __get_internal_psfphot_dir(self):
		'''
		check self.psfphot_method and self.dophot_version
		'''
		if self.psfphot_method == 'daophot':
			psfphotdir = self.psfphot_daophot_iraf_dir
		elif self.psfphot_method == 'dophot':
			if self.dophot_version == 'fortran':
				psfphotdir = self.psfphot_dophot_fortran_dir
			elif self.dophot_version == 'C':
				psfphotdir = self.psfphot_dophot_c_dir
			else:
				raise ValueError('not supported version of dophot')
		else:
			raise ValueError('Sorry, not supported, can not get what you want')

		return psfphotdir

	def __get_internal_psfphot_mag_file(self, imgkey):
		'''
		See __get_internal_psfphot_dir
		'''
		psfphotdir = self.__get_internal_psfphot_dir()
		outfile = os.path.join(psfphotdir,imgkey.split('.')[0]+self.psfphot_ret_file_suffix)

		return outfile


	def __extract_dophot_phot_results(self, image_key):
		'''
		extract and output the clean result in format of x,y,mag,magerr
		'''
		img_key_s = image_key.split('.')[0]
		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		if self.dophot_version == 'fortran':
			psfphotdir = self.psfphot_dophot_fortran_dir
		else:
			psfphotdir = self.psfphot_dophot_c_dir

		psfphot_ret_dst = os.path.join(psfphotdir,output_obj)
		#output result file of dophot photometry with complete format has 15 columns; but elements in some rows stick toghther
		#we need to eliminate wrong rows
		temp_file = os.path.join(psfphotdir, img_key_s + '_temp.out')
		fid_out = open(temp_file, 'w')
		renew_required = False
		for line in open(psfphot_ret_dst).readlines():
			line_segs = line.split()
			if len(line_segs) < 15:
				renew_required = True
				continue
			fid_out.write(line)
		fid_out.close()

		if renew_required:
			shutil.move(temp_file, psfphot_ret_dst)

		psfret_whole = np.loadtxt(psfphot_ret_dst)
		mags = psfret_whole[:,4]
		mask_good = np.where(mags != -99.999)[0]
		ret_good  = psfret_whole[mask_good,:]

		if self.dophot_select_stars_method == 'eq': 	#source type identation 1(star); see docs/dophot_notes.txt for summary of object types
			ret_good_small = ret_good[ret_good[:,1]==self.dophot_select_star_value]	 #only x,y,mag,magerr
		elif self.dophot_select_stars_method == 'lt':
			ret_good_small = ret_good[ret_good[:,1]<self.dophot_select_star_value]
		elif self.dophot_select_stars_method == 'gt':
			ret_good_small = ret_good[ret_good[:,1]>self.dophot_select_star_value]
		else:
			raise ValueError("invalid input for self.dophot_select_stars_method")

		xymags = ret_good_small[:,[2,3,4,5]]
		#save the source regions as a ds9 region file,
		phot_source_reg =  os.path.join(psfphotdir, image_key.split('.')[0] + '_all.reg')
		create_ds9_region_file(ret_good, phot_source_reg, clobber = True,  x_col=2, y_col=3, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		#save the sources as a ds9  reg file, only sources meet the requirement of identification and photometry quality
		phot_source_reg =  os.path.join(psfphotdir, image_key.split('.')[0] + '.reg')
		create_ds9_region_file(xymags, phot_source_reg, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return xymags


	def __clean_dophot_output_files(self, image_key, storagedir):
		'''
		Clean the miscellaneous output from dophot task
		'''
		img_key_s = image_key.split('.')[0]
		finish_file   = img_key_s +'.finish'
		finish_src = os.path.join(os.getcwd(),finish_file)
		finish_dst = os.path.join(storagedir,finish_file)
		if os.path.exists(finish_src):
			shutil.move(finish_src,finish_dst)

		psf_file   = img_key_s +'_psf.fits'
		psffile_src = os.path.join(os.getcwd(),psf_file)
		psffile_dst = os.path.join(storagedir,psf_file)
		if os.path.exists(psffile_src):
			shutil.move(psffile_src,psffile_dst)

		logfile    = img_key_s +'_dophot.log'
		logfile_src = os.path.join(os.getcwd(),logfile)
		logfile_dst = os.path.join(storagedir,logfile)
		if os.path.exists(logfile_src):
			shutil.move(logfile_src,logfile_dst)

		pmfile_out = img_key_s +'out.pm'
		pmfile_out_src = os.path.join(os.getcwd(),pmfile_out)
		pmfile_out_dst = os.path.join(storagedir,pmfile_out)
		shutil.move(pmfile_out_src,pmfile_out_dst)

		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		output_image = 'out_' + image_key
		output_obj_src = os.path.join(os.getcwd(), output_obj)
		psfphot_ret_dst = os.path.join(storagedir,output_obj)
		shutil.move(output_obj_src,psfphot_ret_dst)

		output_img_src = os.path.join(os.getcwd(), output_image)
		psfphot_residual_img_dst = os.path.join(storagedir,output_image)
		if os.path.isfile(output_img_src):
			shutil.move(output_img_src,psfphot_residual_img_dst)

	def __prepare_pm_data_dophot_fortran(self, image_key, photimg):
		'''
		prepare the input data for the pm(parameter modification) file
		'''
		img_key_s = image_key.split('.')[0]

		fwhm = self.photometry_info[image_key]['fwhm']
		skyvalue = self.photometry_info[image_key]['bkg']
		eperdn   = float(self.base.telescopes_info[self.current_telescope]['gain'])
		rdnoise  = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		output_image = 'out_' + image_key
		pmfile_out = img_key_s +'out.pm'
		finishfile = img_key_s + '.finish'

		autopmdata_dict = {}
		autopmdata_dict['FWHM']    = fwhm
		autopmdata_dict['SKY']     = skyvalue
		autopmdata_dict['EPERDN']  = eperdn
		autopmdata_dict['RDNOISE'] = rdnoise
		autopmdata_dict['IMAGE_IN'] = photimg
		autopmdata_dict['IMAGE_OUT']=output_image
		autopmdata_dict['OBJECTS_OUT']=output_obj
		autopmdata_dict['PARAMS_OUT']=pmfile_out
		autopmdata_dict['FINISHFILE'] =finishfile

		dophot_fortran_pm = OrderedDict()
		for kw in self.dophot_fortran_pm.keys():
			dophot_fortran_pm[kw] = self.dophot_fortran_pm[kw]

		for kw in autopmdata_dict.keys():
			if kw not in dophot_fortran_pm:
				dophot_fortran_pm[kw] = autopmdata_dict[kw]
			elif dophot_fortran_pm[kw] is None:
				dophot_fortran_pm[kw] = autopmdata_dict[kw]
			else:
				print "%s explicitly set"%pm
		return dophot_fortran_pm

	def __prepare_pm_file_dophot_fortran(self,pm_file, **dophot_fortran_pm):
		'''
		prepare dophot input file
		'''
        	fid_out = open(pm_file,'wt')
		for kw in dophot_fortran_pm.keys():
			if dophot_fortran_pm[kw] is None:
				continue
			elif isinstance(dophot_fortran_pm[kw],str):
				fid_out.write("%s = '%s' \n"%(kw, dophot_fortran_pm[kw]))
			else:
				fid_out.write("%s = %s \n"%(kw, str(dophot_fortran_pm[kw])))
		fid_out.write('END')
		fid_out.close()

	def __prepare_pm_data_dophot_C(self, image_key, photimg):
		'''
		prepare the input data for the pm(parameter modification) file
		'''
		img_key_s = image_key.split('.')[0]

		fwhm = self.photometry_info[image_key]['fwhm']
		skyvalue = self.photometry_info[image_key]['bkg']
		eperdn = float(self.base.telescopes_info[self.current_telescope]['gain'])
		rdnoise = float(self.base.telescopes_info[self.current_telescope]['readnoise'])
		saturation_value = float(self.base.telescopes_info[self.current_telescope]['satval'])

		output_obj = img_key_s + '.out'	#the psf photometry result directly from the dophot  with COMPLETE data format
		output_image = 'out_' + image_key
		pmfile_out = img_key_s +'out.pm'
		logfile    = img_key_s +'_dophot.log'
		psf_file   = img_key_s +'_psf.fits'
		centintmax = saturation_value

		autopmdata_dict = {}
		autopmdata_dict['FWHM']    = fwhm
		autopmdata_dict['SKY']     = skyvalue
		autopmdata_dict['EPERDN']  = eperdn
		autopmdata_dict['RDNOISE'] = rdnoise
		autopmdata_dict['IMAGE_IN'] = photimg
		autopmdata_dict['IMAGE_OUT']=output_image
		autopmdata_dict['OBJECTS_IN'] = ' '
		autopmdata_dict['OBJECTS_OUT']=output_obj
		autopmdata_dict['PARAMS_OUT'] = pmfile_out
		autopmdata_dict['EMP_SUBRAS_OUT'] = psf_file
		autopmdata_dict['SHADOWFILE_OUT'] = ' '
		autopmdata_dict['ERRORS_OUT']     = ' '
		autopmdata_dict['LOGFILE'] = logfile

		dophot_C_pm = OrderedDict()
		for kw in self.dophot_C_pm.keys():
			dophot_C_pm[kw] = self.dophot_C_pm[kw]

		for kw in autopmdata_dict.keys():
			if kw not in dophot_C_pm:
				dophot_C_pm[kw] = autopmdata_dict[kw]
			elif dophot_C_pm[kw] is None:
				dophot_C_pm[kw] = autopmdata_dict[kw]
			else:
				print "%s explicitly set"%kw

		return dophot_C_pm

	def __prepare_pm_file_dophot_C(self,pm_file, **dophot_C_pm):
		'''
		prepare pm file
		'''
		fid_out = open(pm_file,'wt')
		for kw in dophot_C_pm.keys():
			if dophot_C_pm[kw] is None:
				continue
			elif isinstance(dophot_C_pm[kw],str):
				fid_out.write("%s = '%s' \n"%(kw, dophot_C_pm[kw]))
			else:
				fid_out.write("%s = %s \n"%(kw, str(dophot_C_pm[kw])))
		fid_out.write('END')
		fid_out.close()

# PSF photometry IRAF
	def compare_psf_photometry_dophot_and_irafdaophot(self, imgkey, dophot_version='fortran', match_criteria=2):
		'''
		check the offset
		'''
		imgkey_s = imgkey.split('.')[0]
		if dophot_version == 'fortran':
			dophotdir = self.psfphot_dophot_fortran_dir
		else:
			dophotdir = self.psfphot_dophot_c_dir
		dophot_xymags = np.loadtxt(os.path.join(dophotdir, imgkey_s+self.psfphot_ret_file_suffix))
		daophot_xymags = np.loadtxt(os.path.join(self.psfphot_daophot_iraf_dir, imgkey_s+self.psfphot_ret_file_suffix))

		xymags1,index1, xymags2, index2 = self.__crude_match(dophot_xymags, daophot_xymags, match_criteria)

		fig, ax = plt.subplots(figsize=(7,7))
		ax.errorbar(xymags1[:,2], xymags1[:,2] - xymags2[:,2], yerr=np.sqrt(xymags1[:,3]**2+xymags2[:,3]**2), fmt='ko')
		ax.set_xlabel('DoPhot')
		ax.set_ylabel('DoPhot - IrafDaophot')
		plt.show()

	def iraf_get_psf_model_single_image(self, imgkey, fitrad=None, psfrad=None, varorder=0, \
	peakmax=45000, snrmin=50, edge_width = 60, sround_fraction=0.68, use_tplbase_as_prior = False, tplimgkey=None, xytran_need=True, \
	eyecheck=False, display_tool='matplotlib', pick_mag_ref_pst = False, reuse_pst=True, overwrite=False, \
	Nmin_warn = 20, Nmax=120,apply_reldiff_filter = True, reldiff_cut1 = 0.05,reldiff_cut2=0.01, learning=False, image_dir='raw_image', verbose=0):
		'''
		Get the psf model for image with imgkey
		INPUTS:
			imgkey:
			image_dir:
			fitrad:
			psfrad:
			peakmax: sources with peak value higher than this will be excluded
			snrmin: only sources with SNR > snrmin will be used for PSF modelling
			edge_width: sources within the edge will be excluded from psf stars
			eyecheck:
			pick_mag_ref_pst:
			reuse_pst:
			overwrite:
			Nmin_warn: if psf stars less than this, warning message will generate
			Nmax: maximum number of psf stars will be used
		'''
		log_info = []

		imgkey_s = imgkey.split('.')[0]
		image = self.__get_internal_image(imgkey, which_dir=image_dir)
		print image

		fwhmpsf = self.photometry_info[imgkey]['fwhm']
		if fwhmpsf == self.photometry_info_init_values['fwhm']:
			raise ValueError('FWHM of image %s not available'%imgkey)

		if psfrad is None:
			psfrad = int(fwhmpsf*self.psfrad_nfwhm)
		if fitrad is None:
			fitrad = int(fwhmpsf*self.fitrad_nfwhm)  

		psffile = os.path.join(self.psfmodel_dir, imgkey_s + '_psf.fits')
		if os.path.exists(psffile) and (not overwrite):
			raise ValueError("PSF model for %s already exists. Please set overwrite=1 if you want that."%image)

		pstfile = os.path.join(self.psfmodel_dir, imgkey_s + '.pst.1')

		if os.path.exists(pstfile) and reuse_pst:
			photret2 = Table.read(pstfile, format='daophot')
			photfile_check = pstfile
		else:
			self.__delete_file_if_exist(pstfile)
			NX, NY = self.__get_fits_image_size(image)
			if verbose:
				print "image size: %s %s"%(NX, NY)
			xmax = NX - edge_width
			ymax = NY - edge_width

			self.get_apphot_iraf_parameters(imgkey, saveparfile=False)
			options = self.apphot_iraf_options.copy()
			print options

			if not use_tplbase_as_prior:
				daofind_outfile = os.path.join(self.psfmodel_dir, imgkey_s+'.coo.1')
				self.__delete_file_if_exist(daofind_outfile)
				self.__source_detection_single_apphot_daofind(imgkey, daofind_outfile=daofind_outfile,threshold=7, which_dir=image_dir, datamax=peakmax)
				srs1 = Table.read(daofind_outfile, format='daophot')
				if srs1.masked:
					print "delete rows with empty SHARPNESS and SROUND"
					srs1 = srs1[~srs1['SHARPNESS'].mask]
					srs1 = srs1[~srs1['SROUND'].mask]
				print "First Source detection done"
				if verbose:
					print "Current sources listed below"
					print srs1

				sround = srs1['SROUND'].data
				dist = np.absolute(sround - np.median(sround))
				distcut = dist[np.argsort(dist)][int(np.floor(len(sround)*sround_fraction))]
				srs1 = srs1[dist<=distcut]
				print "Source selectionn done"
				if verbose:
					print "Remaining sources listed below"
					print srs1

				if srs1.masked:
					srs1 = srs1[(srs1['XCENTER']>edge_width).data*(srs1['YCENTER']>edge_width).data*(srs1['XCENTER']<xmax).data*(srs1['YCENTER']<ymax).data]
				else:
					srs1 = srs1[(srs1['XCENTER']>edge_width)*(srs1['YCENTER']>edge_width)*(srs1['XCENTER']<xmax)*(srs1['YCENTER']<ymax)]

				coofile = os.path.join(self.psfmodel_dir, imgkey_s+'.coo.2')
				self.__delete_file_if_exist(coofile)
				iraf_coo_file_id_filter(daofind_outfile, coofile, srs1['ID'].data) #

				photfile = os.path.join(self.psfmodel_dir, imgkey_s+'.mag.1')
				self.__delete_file_if_exist(photfile)
				self.__aperture_photometry_apphot_iraf_single_image(image, coofile, photfile, options, centering = True)

				photret1 = Table.read(photfile, format='daophot')
				photret1 = photret1[(photret1['CIER']==0).data*(photret1['SIER']==0).data*(photret1['PIER']==0).data]  #get rid of results with error in centering, skyfiting and photometry
				photret1 = photret1[~photret1['MERR'].mask]
				if verbose:
					print "Remaining sources after photometry result check"
					print photret1


				photret1 = photret1[(1.087/photret1['MERR']).data>snrmin]
				if len(photret1)<Nmin_warn:
					log_info.append('less than %s psf star candidates available: %s'%(Nmin_warn, len(photret1)))
	 			while len(photret1)>Nmax:
					if snrmin>500:
						break
					snrs = (1.087/photret1['MERR']).data
					#fig, ax = plt.subplots(figsize=(6,6))
					#ax.hist(snrs,bins=20)
					#plt.show()
					snrmin  = snrmin*1.1
					photret1 = photret1[(1.087/photret1['MERR']).data>snrmin]
				if verbose:
					print "Remaining sources after SNR check"
					print photret1
				log_info.append('Final SNR cut: %s'%snrmin)

				photfile_check = os.path.join(self.psfmodel_dir, imgkey_s + '.mag.2') #This is for the purpose of psf star filtering
				self.__delete_file_if_exist(photfile_check)
				coofile_check = os.path.join(self.psfmodel_dir, imgkey_s+'.coo.3')
				self.__delete_file_if_exist(coofile_check)
				#iraf_mag_file_id_filter(photfile, coofile_check, photret1['ID'].data) #This is wrong because daophot.phot doesn't take daophot.phot ouptut as initial cooridnate input
				coos = photret1[['XCENTER', 'YCENTER', 'MAG']]
				coos.write(coofile_check, format='ascii.commented_header')
				options['app'] = options['app'] + fwhmpsf
				self.__aperture_photometry_apphot_iraf_single_image(image, coofile_check, photfile_check, options, centering = False)

				photfile_check_2 = os.path.join(self.psfmodel_dir, imgkey_s + '.mag.3')
				self.__delete_file_if_exist(photfile_check_2)
				options['app'] = options['app'] + fwhmpsf
				self.__aperture_photometry_apphot_iraf_single_image(image, coofile_check, photfile_check_2, options, centering = False)


				photret2 = Table.read(photfile_check, format='daophot')
				print len(photret1), len(photret2)
				if verbose:
					print "Current remaining sources:"
					print photret2

				if apply_reldiff_filter:
					reldiff_flux_1 = np.absolute(photret2['FLUX']/photret1['FLUX']-1)
	
					photret2_2 = Table.read(photfile_check_2, format='daophot')
					print len(photret2), len(photret2_2)
					reldiff_flux_2 = np.absolute(photret2_2['FLUX']/photret2['FLUX']-1)
	
					if learning:
						showdata = photret2.copy()
						showdata.add_column(Column(data=reldiff_flux_1, name='reldiff1'))
						showdata.add_column(Column(data=reldiff_flux_2, name='reldiff2'))
						for star in showdata[['XCENTER','YCENTER','reldiff1', 'reldiff2']]:
							print star
	
						fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6))
						ax1.hist(reldiff_flux_1, bins=80)
						ax2.hist(reldiff_flux_2, bins=80)
						plt.show()
	
					photret2 = photret2[(reldiff_flux_1<reldiff_cut1).data*(reldiff_flux_2<reldiff_cut2).data]
					#The first star is important to set the flux scale of the PSF model
				distances = np.sqrt((photret2['XCENTER']-NX/2.0)**2 + (photret2['YCENTER']-NY/2.0)**2)
				photret2 = photret2[np.argsort(distances)]

			else:
				if tplimgkey is None:
					raise ValueError('If want to use template image PSF star list as prior input, template image key is required')
				pstfile_prior = os.path.join(self.psfmodel_dir, tplimgkey.split('.')[0] + '.pst.1')
				if not os.path.exists(pstfile_prior):
					raise ValueError('The psf star file %s not exists'%pstfile_prior)
				pststars_prior = Table.read(pstfile_prior, format='daophot')

				if xytran_need:
					trans_fitting = os.path.join(self.stars_dir,imgkey_s + '_tpl.coef')
					if not os.path.exists(trans_fitting):
						raise ValueError('The tranformation file %s is required, please see self.get_target_xys'%trans_fitting)
					porder,dxfit,dyfit = self.__extract_transformation_fitting(trans_fitting)
					xys = self.__transform_xys(dxfit,dyfit, porder, pststars_prior[['XCENTER', 'YCENTER']])
					xys = xys[(xys[:,0]>edge_width)*(xys[:,0]<xmax)*(xys[:,1]>edge_width)*(xys[:,1]<ymax)]
				else:
					xys = np.array([pststars_prior['XCENTER'].data.data, pststars_prior['YCENTER'].data.data]).transpose()
					imagedata = fits.open(image)[0].data
					okindexs = []
					for i,xy in enumerate(xys):
						x = int(xy[0])
						y = int(xy[1])
						if x-psfrad<0 or x+psfrad>NX-1 or y-psfrad<0 or y+psfrad>NY-1:
							continue
						if imagedata[y,x-psfrad] == 0 or imagedata[y,x+psfrad] == 0 or imagedata[y-psfrad,x] == 0 or imagedata[y+psfrad,x] == 0:
							continue
						okindexs.append(i)
					xys = xys[okindexs]

				coofile = os.path.join(self.psfmodel_dir, imgkey_s+'.coo.2')
				self.__delete_file_if_exist(coofile)
				np.savetxt(coofile, xys)

				photfile_check = os.path.join(self.psfmodel_dir, imgkey_s + '.mag.2') #This is for the purpose of psf star filtering
				self.__delete_file_if_exist(photfile_check)
				self.__aperture_photometry_apphot_iraf_single_image(image, coofile, photfile_check, options, centering = False)
				photret2 = Table.read(photfile_check, format='daophot')
				if verbose:
					print "Current remaining source list"
					print photret2


		#add other constraints here
		if eyecheck:
			if display_tool=='ds9':
				print "Please select stars which to be removed from psf star list"
				photret2, selectedindexs = self.select_or_delete_rows_from_given_table_by_pick_on_ds9_display(image, photret2, xcol='XCENTER', ycol='YCENTER', mode='d', match_criteria = 5, coordinate_type='image', radec_deg=True, circle_radius=psfrad)
			elif display_tool == 'matplotlib':
				imagedata = fits.open(image)[0].data
				self.show_source_stamps(imagedata, photret2, 'XCENTER', 'YCENTER', halfwidth=psfrad, baseindex=1)
				print "Choose to select[s] from the list or delete[d] from the list:"
				mode = raw_input('to delete or to select (default: [d])') or 'd'
				inds = raw_input('input the IDs of the wanted sources:')
				if inds!='':
					inds = [int(ind) for ind in  inds.split(',')]
				else:
					inds = []
				if mode.startswith('s'):
					photret2 = photret2[inds]
				elif mode.startswith('d'):
					inds_use = list(np.arange(len(photret2)))
					for ind in inds:
						inds_use.remove(ind)
					photret2 = photret2[inds_use]
			else:
				raise ValueError('display tool %s not supported'%display_tool)

		if pick_mag_ref_pst:
			print "Please select the star to be used for PSF mag"
			photrettemp, selectedindexs = self.select_or_delete_rows_from_given_table_by_pick_on_ds9_display(image, photret2, xcol='XCENTER', ycol='YCENTER', mode='s', match_criteria = 5, coordinate_type='image', radec_deg=True, circle_radius=psfrad)
			idstemp = [ii for ii in np.arange(len(photret2)) if ii not in selectedindexs]
			idstemp.insert(0, selectedindexs[0])
			print idstemp
			photret2 = photret2[idstemp]

		#pstfile = os.path.join(self.psfmodel_dir, imgkey_s + '.pst.1')
		#self.__delete_file_if_exist(pstfile)
		if os.path.exists(pstfile):
			if os.path.samefile(photfile_check, pstfile):
				pstfile_copy = pstfile + '.copy'
				shutil.copy(pstfile, pstfile_copy)
				photfile_check = pstfile_copy
		self.__delete_file_if_exist(pstfile)
		ids = photret2['ID'].data   #!!!!!
		iraf_mag_file_id_filter(photfile_check, pstfile, ids, firstfirst=True)

		#psffile = os.path.join(self.psfmodel_dir, imgkey_s + '_psf.fits')
		self.__delete_file_if_exist(psffile)
		daophot_psf_iraf(image, photfile_check, pstfile, psffile, opstfile=None, groupfile=None, psfrad=psfrad, fitrad=fitrad, varorder=varorder)

		psfimage = os.path.join(self.psfmodel_dir, imgkey_s + '_psf_image.fits')
		self.__delete_file_if_exist(psfimage)
		if varorder >0:
			xpsf = NX/2
			ypsf = NY/2
		else:
			xpsf = 'INDEF'
			ypsf = 'INDEF'
		daophot_seepsf_iraf(psffile, psfimage,  xpsf=xpsf, ypsf=ypsf)

		return log_info

	def iraf_psf_photometry_peak(self, image_key, image =None, psffile=None, photfile=None, \
	peakfile=None, x=None, y=None, fitrad=None, psfrad=None, fitsky='yes', recenter='yes', \
	ap_phot_recenter_xy = True, apradius=None, skyin=None, skywidth=None, residual_image = True, subimage=None, \
	residual_smallcut=False, cutout_width=120, update_phot_record=False, which_dir='raw_image'):
		'''
		IRAF PSF photometry for sources with peak task
		INPUTS:
			image_key:
			psffile: if None, then imgkey_s + '_psf.fits' in self.psfmodel_dir will be used
			photfile: the filename for the output of the aperture photometry 
			peakfile: the filename for the ouptut of peak task
			x, y: source coordinate; scaler for single source or list-like input for multiple sources; \
				if None, then target coordinate for image_key will be used
			fitrad, psfrad: parameters for daophot peak task
			fitsky: refit sky during peak task
			recenter: recenter during dophot peak task
			ap_phot_recenter_xy: whether to recenter the source during daophot task
			skyin: radius of sky annulus during daophot task
			skywidth: width of sky annulus during daophot task
			residual_image: output residual image
			subimage: residual image filename
			residual_smallcut:
			cutout_width:
			update_phot_record:
			which_dir:
		'''
		imgkey_s = image_key.split('.')[0]
		if image is None:
			image = self.__get_internal_image(image_key, which_dir=which_dir)
		if psffile is None:
			psffile = os.path.join(self.psfmodel_dir, imgkey_s + '_psf.fits')
		if not os.path.exists(psffile):
			raise ValueError('PSF model for image %s not avaiable'%image_key)
		if x is None or y is None:
			x = self.photometry_info[image_key]['x']
			y = self.photometry_info[image_key]['y']
		coo_list =  os.path.join(self.psfphot_daophot_iraf_dir, imgkey_s+'_temp.coo')
		f = open(coo_list, 'w')
		if np.isscalar(x) and np.isscalar(y):
			xc = x
			yc = y
			f.write(str(x) + ' ' + str(y) + '\n')
		else:
			xc = x[0]
			yc = y[0]
			for xi, yi in zip(x,y):
				f.write(str(xi) + ' ' + str(yi) + '\n')
		f.flush()
		f.close

		self.get_apphot_iraf_parameters(image_key)
		options = self.apphot_iraf_options.copy()
		if apradius is not None:
			options['app'] = apradius
		if skyin is not None:
			options['skyin'] = skyin
		if skywidth is not None:
			options['skywidth'] = skywidth
		if photfile is None:
			photfile = os.path.join(self.aperture_photometry_dir, imgkey_s+'_xys.mag')
		self.__delete_file_if_exist(photfile)
		self.__aperture_photometry_apphot_iraf_single_image(image, coo_list, photfile, options, centering = ap_phot_recenter_xy)


		rejfile = ''
		if peakfile is None:
			peakfile = os.path.join(self.psfphot_daophot_iraf_dir, imgkey_s+'.peak')
		self.__delete_file_if_exist(peakfile)
		fwhm = self.photometry_info[image_key]['fwhm']

		hinfo = get_fits_info(psffile, ['FITRAD','PSFRAD'])
		if fitrad is None:
			fitrad = hinfo['FITRAD']
		if psfrad is None:
			psfrad = hinfo['PSFRAD']
		daophot_peak_iraf(image, photfile, psffile, peakfile, rejfile, fitrad=fitrad, psfrad=psfrad, refitsky=fitsky, recenter=recenter)

		if residual_image:
			exfile = ''
			if subimage is None:
				subimage = os.path.join(self.psfphot_daophot_iraf_dir, '%s_xys_model_sub.fits'%(imgkey_s))
			self.__delete_file_if_exist(subimage)
			daophot_substar_iraf(image, peakfile, exfile, psffile, subimage, fitrad=fitrad, psfrad=psfrad)
			if residual_smallcut:
				hw = int(cutout_width/2)
				inimage = subimage+'[%s:%s,%s:%s]'%(str(int(xc)-hw), str(int(xc)+hw), str(int(yc)-hw), str(int(yc)+hw))
				imcopy_iraf(inimage, subimage)

		peakret = Table.read(peakfile, format='daophot')
		if update_phot_record:
			mag = peakret['MAG'].data[0]
			magerr = peakret['MERR'].data[0]
			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr

		return peakret


	def psf_photometry_iraf_single_image(self, image_key,fitrad=None, psfrad=None, \
	redo_daofind=False, redo_phot=False, threshold=5, fitsky='no', sannulus=None, \
	wsannulus=None,  which_dir='raw_image'):
		'''
		IRAF PSF photometry for single image
		INPUTS:
			image_key:
			fitrad: fitting radius
			psfrad: psfmode radius
			redo_daofind: if coordinate list file exists, redo?
			redo_phot: if photfile exists, redo?
			threshold: if redo_daofind, the threshold
			fitsky: recompute group sky value during the fit in allstar task?
			sannulus: if fitsky, the inner radius of sky annulus
			wsannulus: if fitsky, the width of sky annulus
			which_dir:
		'''
		imgkey_s = image_key.split('.')[0]
		image = self.__get_internal_image(image_key, which_dir=which_dir)

		psffile = os.path.join(self.psfmodel_dir, imgkey_s + '_psf.fits')
		if not os.path.exists(psffile):
			raise ValueError('PSF model for image %s not avaiable'%image_key)

		photfile = os.path.join(self.aperture_photometry_dir, imgkey_s+'.mag')
		if (not os.path.exists(photfile)) or redo_phot:
			coofile = os.path.join(self.stars_dir, imgkey_s+'.coo')
			if (not os.path.exists(coofile)) or redo_daofind:
				self.__delete_file_if_exist(coofile)
				self.source_detection_single_apphot_daofind(image_key, coofile, which_dir=which_dir, threshold=threshold)
			self.aperture_photometry_apphot_iraf_single_image(image_key, overwrite=True, which_dir=which_dir, aperture_size_fixed=False, centering=True)

		allstarfile = os.path.join(self.psfphot_daophot_iraf_dir, imgkey_s+'.allstar')
		rejfile = os.path.join(self.psfphot_daophot_iraf_dir, imgkey_s+'.rej')
		subimage = os.path.join(self.psfphot_daophot_iraf_dir, imgkey_s+'_sub.fits')
		self.__delete_file_if_exist(allstarfile)
		self.__delete_file_if_exist(rejfile)
		self.__delete_file_if_exist(subimage)
		hinfo = get_fits_info(psffile, ['FITRAD','PSFRAD'])
		if fitrad is None:
			fitrad = hinfo['FITRAD']
		if psfrad is None:
			psfrad = hinfo['PSFRAD']

		if sannulus is None:
			sannulus = int(psfrad/2.0)
		if wsannulus is None:
			wsannulus = int(psfrad/2.0)
		daophot_allstar_iraf(image, photfile, psffile, allstarfile, rejfile, subimage, fitrad=fitrad, psfrad=psfrad, fitsky=fitsky, sannulus=sannulus, wsannulus=wsannulus)

		psfphotret = Table.read(allstarfile, format='daophot')
		out = psfphotret[['XCENTER','YCENTER','MAG','MERR']]
		psfphot_retfile = self.__get_internal_psfphot_mag_file(image_key)
		out.write(psfphot_retfile, format='ascii.commented_header')


	def psf_photometry_iraf_single_target(self, image_key, search_radius=3, which_dir='raw_image'):
		'''
		PSF photometry for the target of interest
		'''
		imgkey_s = image_key.split('.')[0]
		psfphot_ret_file = self.__get_internal_psfphot_mag_file(image_key)
		psfphot_ret = np.loadtxt(psfphot_ret_file)

		xys = psfphot_ret[:,0:2]
		mags = psfphot_ret[:,2]
		magerrs = psfphot_ret[:,3]
		x = self.photometry_info[image_key]['x']
		y = self.photometry_info[image_key]['y']
		xy = [x, y]

		exist,xy_match,indice_match = self.__find_corresponding(xys, xy, search_radius)
		if exist and len(indice_match)==1:
			self.photometry_info[image_key]['instmag'] = mags[indice_match]
			self.photometry_info[image_key]['instmagerr'] = magerrs[indice_match]
		else:
			photret = self.iraf_psf_photometry_peak(image_key, x=x, y=y, recenter='yes', residual_smallcut=True, update_phot_record=True, which_dir=which_dir)



	def psf_photometry_iraf_flt(self, flt, dolist = None, renew_psfmodel=False, renew_phot_image=False, \
	renew_phot_target=False, snrmin=50, edge_width=60, tplpsf_prior=False, tplimgkey=None, target_phot_indpt = True, \
	hybrid_model = False, niter_hybrid=5, maxsmas=None, galras=None, galdecs=None, galxs=None, galys=None, starxs=None, starys=None, peakfitsky='yes', \
	allstarfitsky='yes', recenter='yes', which_dir='raw_image', debug=0):
		'''
		See self.psf_photometry_iraf_single_image and self.iraf_psf_photometry_peak
		INPUTS:
			flt:
			dolist: if not None, only imgkey in this list will be done
			maxsmas: maximum semiaxis length for the isophote modelling; only required when hybr == True
			renew_psfmodel:
			renew_phot_image: renew the photometry for the whole image
			renew_phot_target: renew the photometry for the target
			snrmin, edge_width: see self.iraf_get_psf_model_single_image
			tplpsf_prior, tplimgkey: see self.iraf_get_psf_model_single_image
			target_phot_indpt: whether the photometry for the target is independent of the photometry for the whole image
			hybrid_model: whether use hybrid model (PSF for the target and isophote for the background galaxy) for the photometry of target; only matters when target_phot_indpt == True
			galras, galdecs: galaxy center in degree; only required when hybrid_model == True
			galxs, galys:
			starxs, starys: list of stars to be masked for hybrid model
			recenter: whether recenter during the peak task
		'''
		failed_list = []
		for imgkey in self.images.keys():
			if dolist is not None:
				if imgkey not in dolist:
					continue
			if self.photometry_info[imgkey]['flt'] != flt:
				continue
			print imgkey
			imgkey_s = imgkey.split('.')[0]
			if self.photometry_info[imgkey]['drop'] != 0:
				continue

			psffile = os.path.join(self.psfmodel_dir, imgkey_s+'_psf.fits')
			if (not os.path.exists(psffile)) or renew_psfmodel:
				if tplpsf_prior and (tplimgkey is not None) and tplimgkey!=imgkey:
					use_tplbase_as_prior = True
				else:
					use_tplbase_as_prior = False
				self.iraf_get_psf_model_single_image(imgkey, use_tplbase_as_prior = use_tplbase_as_prior, tplimgkey=tplimgkey, image_dir=which_dir, eyecheck=False, overwrite=True, snrmin=snrmin, edge_width=edge_width)

			psfphot_ret_file = self.__get_internal_psfphot_mag_file(imgkey)
			if (not os.path.exists(psfphot_ret_file)) or renew_phot_image:
				self.psf_photometry_iraf_single_image(imgkey, which_dir=which_dir, fitsky=allstarfitsky)

			if self.photometry_info[imgkey]['instmag'] == self.photometry_info_init_values['instmag'] or renew_phot_target:
				try:
					if target_phot_indpt:
						if hybrid_model:
							if maxsmas is None:
								raise ValueError('maxsmas is required for hybrid model')
							if (galras is None or galdecs is None) and (galxs is None or galys is None):
								raise ValueError('(galras, galdecs) or (galxs, galys) is required for hybrid model')
							if galras is not None and galdecs is not None:
								if starxs is not None:
									print "# WARNING: starlist only allowed when galx, galy are used"
								self.psf_photometry_iraf_single_target_hybrid(imgkey, galras, galdecs, maxsmas, gal_coord_image=False, niter=niter_hybrid, which_dir=which_dir, fitsky=peakfitsky)
							else:
								#if tplimgkey is None:
								#	raise ValueError('if want to used (galx, galy), the reference image is required')
								if imgkey == tplimgkey:
									self.psf_photometry_iraf_single_target_hybrid(imgkey, galxs, galys, maxsmas, gal_coord_image=True, xstars=starxs, ystars=starys, niter=niter_hybrid, which_dir=which_dir, fitsky=peakfitsky)
								else:
									trans_fitting = os.path.join(self.stars_dir,imgkey_s + '_tpl.coef')
									if not os.path.exists(trans_fitting):
										raise ValueError('The tranformation file %s is required, please see self.get_target_xys'%trans_fitting)
									porder,dxfit,dyfit = self.__extract_transformation_fitting(trans_fitting)
									galxcs = []
									galycs = []
									for galx, galy in zip(galxs, galys):
										tpl_galxy = [galx, galy]
										target_galxy = self.__transform_xys(dxfit,dyfit, porder, tpl_galxy)
										#print target_galxy
										galxc = round(target_galxy[0],1)
										galyc = round(target_galxy[1],1)
										galxcs.append(galxc)
										galycs.append(galyc)
									if starxs is not None and starys is not None:
										xstars = []
										ystars = []
										for x,y in zip(starxs, starys):
											xy  = target_galxy = self.__transform_xys(dxfit,dyfit, porder, [x,y])
											xstars.append(xy[0])
											ystars.append(xy[1])
									else:
										xstars = None
										ystars = None
									self.psf_photometry_iraf_single_target_hybrid(imgkey, galxcs, galycs, maxsmas, gal_coord_image=True, xstars=xstars, ystars=ystars, niter=niter_hybrid, which_dir=which_dir, fitsky=peakfitsky)
						else:
							x = self.photometry_info[imgkey]['x']
							y = self.photometry_info[imgkey]['y']
							self.iraf_psf_photometry_peak(imgkey, x=x, y=y, recenter=recenter, residual_smallcut=True, update_phot_record=True, which_dir=which_dir, fitsky=peakfitsky)
					else:
						self.psf_photometry_iraf_single_target(imgkey, which_dir=which_dir)
				except Exception as e:
					if debug:
						raise ValueError(e)
					else:
						failed_list.append(imgkey)
		return failed_list

	def iraf_get_isophote_model(self, imgkey,  x0, y0, minsma, maxsma, \
	input_image=None, isoimage=None, residual_image = False, subimage=None, bkgvalue=0, isokwargs=None, which_dir='raw_image'):
		'''
		galaxy isophote modelling
		INPUTS:
			imgkey:
			x0:
			y0:
			minsma:
			maxsma:
			bkgvalue:
		'''
		imgkey_s = imgkey.split('.')[0]
		if input_image is None:
			input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		isomodel = os.path.join(self.isophote_dir, imgkey_s+'_iso.dat')
		self.__delete_file_if_exist(isomodel)
		if isokwargs is None:
			kwargs = self.ellipse_pm.copy()
		else:
			kwargs = isokwargs
		stsdas_analysis_isophote_ellipse(input_image, isomodel, x0, y0, minsma, maxsma, **kwargs)

		if isoimage is None:
			isoimage = os.path.join(self.isophote_dir, imgkey_s+'_iso.fits')
		self.__delete_file_if_exist(isoimage)
		stsdas_analysis_isophote_bmodel(isomodel, isoimage, parent_image=input_image, backgr=bkgvalue)
		#stsdas_analysis_isophote_cmodel(isomodel, isoimage, parent_image=input_image, backgr=bkgvalue)

		if residual_image:
			if subimage is None:
				subimage = os.path.join(self.isophote_dir, imgkey_s + '_iso_sub.fits')
			imarith_iraf(input_image, '-', isoimage, subimage)


	def psf_photometry_iraf_single_target_hybrid(self, imgkey, xgals, ygals, maxsmas, minsmas=None, fitsky='yes', \
	xstars=None, ystars=None, use_archival_source_info=False, gal_coord_image=True, niter=1, gal_subimage_smcut=True, \
	isoimage_smcut=True, hybrid_subimage_smcut=True, peak_subimage_smcut=True, which_dir='raw_image', verbose=1):
		'''
		PSF photometry combined with isophote model for photometry of target embeded in host galaxy
		Updates 20220228: allow multiple galaxy models to be subtracted
		STEPS: (1) psf photometry of the stellar objects and subtract them; (2) model and subtract the galaxy from the image one after another starting from the biggest galaxy on the star-subtracted image; (3) subtract all galaxy models from the original image; (4) perform psf photometry on the residual image of step (3); (5) redo step (2),(3),(4). One difference here in (2) is that before modelling specific galaxy, not only subtracting the stars but also other galaxies; (6) subtract all galaxies and non-target stars from the original image and get the photometry of the target.
		INPUTS:
			imgkey:
			xgals:
			ygals:
			maxsmas:
			fitsky: do not control the fitsky parameter in the first peak task 
			xstars: x coordinates for point sources other than the target
			ystars:
			gal_coord_image: whether the input galaxy center coordinate (xgal, ygal) is in image coordinate
			niter:
			which_dir:
		'''

		image = self.__get_internal_image(imgkey, which_dir=which_dir)
		imgkey_s= imgkey.split('.')[0]
		xtarget = self.photometry_info[imgkey]['x']  #
		ytarget = self.photometry_info[imgkey]['y']

		galmodelfiles = {}
		starspsfdata = []

		galinfofile = os.path.join(self.isophote_dir, imgkey_s+'_gals_info.csv')
		starinfofile = os.path.join(self.isophote_dir, imgkey_s+'_bkgstars_info.csv')
		if use_archival_source_info:
			if not os.path.exists(galinfofile):
				raise ValueError('the galaxy info file %s not exists'%galinfofile)
			galsinfo = Table.read(galinfofile, format='ascii.csv')
			xgals = galsinfo['xc']
			ygals = galsinfo['yc']
			maxsmas = galsinfo['maxsma']
			minsmas = galsinfo['minsma']
			if os.path.exists(starinfofile):
				starsinfo = Table.read(starinfofile, format='ascii.csv')
				xstars = starsinfo['xc']
				ystars = starsinfo['yc']
			else:
				if verbose:
					print "No archival info for background stars"


		if xstars is not None and ystars is not None:
			if np.isscalar(xstars):
				xstars = [xstars]
			if np.isscalar(ystars):
				ystars = [ystars]
			xstars.insert(0, xtarget)
			ystars.insert(0, ytarget)
		else:
			xstars = xtarget
			ystars = ytarget

		if np.isscalar(xgals) and np.isscalar(ygals) and np.isscalar(maxsmas):
			if minsmas is None:
				minsmas = 5
			xgals = [xgals]
			ygals = [ygals]
			maxsmas = [maxsmas]
			minsmas = [minsmas]
			Ngal = 1
		else:
			if len(xgals) != len(xgals) or len(ygals) != len(maxsmas):
				raise ValueError("the number of galaxy (x,y,sma) does not match")
			Ngal = len(xgals)
			if minsmas is None:
				minsmas = np.ones(Ngal)*5

		if not gal_coord_image:
			world_radec = np.array([xgals, ygals]).reshape(Ngal,2)
			image_xy = self.__world2image(image, world_radec)
			xgals = image_xy[:,0]
			ygals = image_xy[:,1]
		
		if not use_archival_source_info:
			galsinfo = Table([xgals, ygals, maxsmas, minsmas], names =['xc','yc','maxsma','minsma'])
			starsinfo = Table([xstars, ystars], names=['xc', 'yc'])
			galsinfo.write(galinfofile, format='ascii.csv', overwrite=True)
			starsinfo.write(starinfofile, format='ascii.csv', overwrite=True)

		psffile = os.path.join(self.psfmodel_dir, imgkey_s+'_psf.fits')
		peak_subimage= os.path.join(self.isophote_dir, imgkey_s+'_peak_sub.fits')  #PSF modelled star subtracted image
		peakfile = os.path.join(self.isophote_dir, imgkey_s+'_xys.peak')
		spphot = self.iraf_psf_photometry_peak(imgkey, psffile=psffile, peakfile=peakfile, x=xstars, y=ystars, fitrad=None, psfrad=None, fitsky='yes', recenter='no', ap_phot_recenter_xy = False, residual_image = True, subimage=peak_subimage, residual_smallcut=False, update_phot_record=False, which_dir=which_dir)
		print spphot

		starfree_galsub_image = os.path.join(self.isophote_dir, imgkey_s+'_starfree_galsub.fits')
		shutil.copy(peak_subimage, starfree_galsub_image)
		gal_subimage = os.path.join(self.isophote_dir, imgkey_s+'_gal_sub.fits') #stars stay
		shutil.copy(image, gal_subimage)	
		for i, (xgal, ygal, minsma, maxsma) in enumerate(zip(xgals, ygals, maxsmas, minsmas)):
			isoimage = os.path.join(self.isophote_dir, imgkey_s+'_%s_%s_%s_iso.fits'%(str(xgal), str(ygal), str(maxsma)))
			self.__delete_file_if_exist(isoimage)
			if i==0:
				bkgvalue = self.photometry_info[imgkey]['bkg']
			else:
				bkgvalue = 0
			isokwargs = self.ellipse_pm.copy()
			isokwargs['sma0'] = np.max([maxsma*0.2, 5])
			self.iraf_get_isophote_model(imgkey,  xgal, ygal, minsma, maxsma, bkgvalue=bkgvalue, input_image= starfree_galsub_image, isoimage=isoimage, isokwargs=isokwargs, residual_image= False, which_dir=which_dir)
			galmodelfiles[str(i)] = isoimage
			starfree_galsub_image_temp = os.path.join(self.isophote_dir, imgkey_s+'_starfree_galsub_temp.fits')
			self.__delete_file_if_exist(starfree_galsub_image_temp)
			imarith_iraf(starfree_galsub_image, '-', isoimage, starfree_galsub_image_temp)
			shutil.copy(starfree_galsub_image_temp, starfree_galsub_image) 

			starsafe_galsub_image = os.path.join(self.isophote_dir, imgkey_s+'_gal_sub_temp.fits')
			self.__delete_file_if_exist(starsafe_galsub_image)
			imarith_iraf(gal_subimage, '-', isoimage, starsafe_galsub_image)
			shutil.copy(starsafe_galsub_image, gal_subimage)

		iteri = 1
		while iteri<niter:
			print "iteration #"+str(i)
			peakfile = os.path.join(self.isophote_dir, imgkey_s+'_xys.peak')
			spphot = self.iraf_psf_photometry_peak(imgkey, image=gal_subimage, psffile=psffile, peakfile=peakfile, x=xstars, y=ystars, fitrad=None, psfrad=None, fitsky=fitsky, recenter='yes', ap_phot_recenter_xy = False, residual_image = False, subimage=None, residual_smallcut=True, update_phot_record=False, which_dir=which_dir)
			print spphot
			exfile = ''
			peak_subimage= os.path.join(self.isophote_dir, imgkey_s+'_peak_sub.fits')
			self.__delete_file_if_exist(peak_subimage)
			hinfo = get_fits_info(psffile, ['FITRAD','PSFRAD'])
			fitrad = hinfo['FITRAD']
			psfrad = hinfo['PSFRAD']
			daophot_substar_iraf(image, peakfile, exfile, psffile, peak_subimage, fitrad=fitrad, psfrad=psfrad)

			starfree_galsub_image = os.path.join(self.isophote_dir, imgkey_s+'_starfree_galsub.fits')
			shutil.copy(peak_subimage, starfree_galsub_image)
			gal_subimage = os.path.join(self.isophote_dir, imgkey_s+'_gal_sub.fits') #stars stay
			shutil.copy(image, gal_subimage)	
			for i, (xgal, ygal, minsma, maxsma) in enumerate(zip(xgals, ygals, minsmas, maxsmas)):
				isoimage = os.path.join(self.isophote_dir, imgkey_s+'_%s_%s_%s_iso.fits'%(str(xgal), str(ygal), str(maxsma)))
				self.__delete_file_if_exist(isoimage)
				if i==0:
					bkgvalue = self.photometry_info[imgkey]['bkg']
				else:
					bkgvalue = 0

				starfree_galsub_image_fit =  os.path.join(self.isophote_dir, imgkey_s+'_gal%s.fits'%str(i))#
				shutil.copy(peak_subimage, starfree_galsub_image_fit)
				for j in np.arange(Ngal):
					if j == i:
						continue
					starfree_galsub_image_fit_temp = os.path.join(self.isophote_dir, imgkey_s+'_starfree_galsub_fit_temp.fits')
					self.__delete_file_if_exist(starfree_galsub_image_fit_temp)
					isoimage_old = galmodelfiles[str(j)]
					imarith_iraf(starfree_galsub_image_fit, '-', isoimage_old, starfree_galsub_image_fit_temp)
					shutil.copy(starfree_galsub_image_fit_temp, starfree_galsub_image_fit) 
					
				isokwargs = self.ellipse_pm.copy()
				isokwargs['sma0'] = np.max([maxsma*0.2, 5])
				self.iraf_get_isophote_model(imgkey,  xgal, ygal, minsma, maxsma, bkgvalue=bkgvalue, input_image= starfree_galsub_image_fit, isoimage=isoimage, residual_image= False, which_dir=which_dir)
				galmodelfiles[str(i)] = isoimage #check this 
			
				starfree_galsub_image_temp = os.path.join(self.isophote_dir, imgkey_s+'_starfree_galsub_temp.fits')
				self.__delete_file_if_exist(starfree_galsub_image_temp)
				imarith_iraf(starfree_galsub_image, '-', isoimage, starfree_galsub_image_temp)
				shutil.copy(starfree_galsub_image_temp, starfree_galsub_image) 

				starsafe_galsub_image = os.path.join(self.isophote_dir, imgkey_s+'_gal_sub_temp.fits')
				self.__delete_file_if_exist(starsafe_galsub_image)
				imarith_iraf(gal_subimage, '-', isoimage, starsafe_galsub_image)
				shutil.copy(starsafe_galsub_image, gal_subimage)

			iteri += 1

		subimage = os.path.join(self.isophote_dir, imgkey_s+'_hybrid_subimage.fits')
		spphot = self.iraf_psf_photometry_peak(imgkey, image=gal_subimage, psffile=psffile, x=xtarget, y=ytarget, fitrad=None, psfrad=None, fitsky=fitsky, recenter='yes', ap_phot_recenter_xy = False, residual_image = True, subimage=subimage, residual_smallcut=False, update_phot_record=True, which_dir=which_dir)

		xgal = xgals[0]
		ygal = ygals[0]	
		maxsma = maxsmas[0]
		if peak_subimage_smcut:
			inimage = peak_subimage+'[%s:%s,%s:%s]'%(str(int(xgal)-maxsma), str(int(xgal)+maxsma), str(int(ygal)-maxsma), str(int(ygal)+maxsma))
			imcopy_iraf(inimage, peak_subimage)

		if hybrid_subimage_smcut:
			inimage = subimage+'[%s:%s,%s:%s]'%(str(int(xgal)-maxsma), str(int(xgal)+maxsma), str(int(ygal)-maxsma), str(int(ygal)+maxsma))
			imcopy_iraf(inimage, subimage)

		if isoimage_smcut:
			#for i, (xgal, ygal, maxsma) in enumerate(zip(xgals, ygals, maxsmas)):
			#	isoimage = galmodelfiles[str(i)]
			#	inimage = isoimage+'[%s:%s,%s:%s]'%(str(int(xgal)-maxsma), str(int(xgal)+maxsma), str(int(ygal)-maxsma), str(int(ygal)+maxsma))
			for i in np.arange(Ngal):
				isoimage = galmodelfiles[str(i)]
				inimage = isoimage+'[%s:%s,%s:%s]'%(str(int(xgal)-maxsma), str(int(xgal)+maxsma), str(int(ygal)-maxsma), str(int(ygal)+maxsma))
				imcopy_iraf(inimage, isoimage)

		if gal_subimage_smcut:
			inimage = gal_subimage+'[%s:%s,%s:%s]'%(str(int(xgal)-maxsma), str(int(xgal)+maxsma), str(int(ygal)-maxsma), str(int(ygal)+maxsma))
			imcopy_iraf(inimage, gal_subimage)

		print spphot

#Photometry files
	def select_photometry_measurement_given_region(self, image_key,xl=None,xh=None,yl=None,yh=None, photometry_method='apphot', outfile=None):
		'''
		select photometry points within [xl:xh,yl:yh]

		The selection works on the phootometry output file (self.apphot_ret_file_suffix for aperture photometry and self.psfphot_ret_file_suffix for psf photometry)
		and the selected result will be written into 'outfile' or overwritten into original photometry output file if 'outfile' is not provided

		'''
		imgkey = image_key.split('.')[0]
		if photometry_method == 'apphot':
			infile = os.path.join(self.aperture_photometry_dir, imgkey+self.apphot_ret_file_suffix)
		elif photometry_method == 'psfphot':
			infile = self.__get_internal_psfphot_mag_file(image_key)
		else:
			raise ValueError("not support")
		data = np.loadtxt(infile)
		if xl is not None:
			data = data[data[:,0]>xl,:]
		if xh is not None:
			data = data[data[:,0]<xh, :]

		if yl is not None:
			data = data[data[:,1]>yl,:]
		if yh is not None:
			data = data[data[:,1]<yh, :]

		if outfile is None:
			outfile = infile
			print "orginal file %s will be overwrited"%infile
		np.savetxt(outfile, data, fmt='%6.2f %6.2f %5.2f %5.2f')

	def delete_photometry_measurement_given_region(self, image_key, x_center, y_center, radius, infile = None, xcol_infile= 0, ycol_infile = 1,  photometry_method='apphot', outfile=None):

		'''
		delete photometry points within circle at given center with given radius
		'''
		imgkey = image_key.split('.')[0]
		if infile is None:
			if photometry_method == 'apphot':
				infile = os.path.join(self.aperture_photometry_dir, imgkey+self.apphot_ret_file_suffix)
			elif photometry_method == 'psfphot':
				infile  = self.__get_internal_psfphot_mag_file(image_key)
			else:
				raise ValueError("not support")
		data = np.loadtxt(infile)
		x_image = data[:,xcol_infile]
		y_image = data[:,ycol_infile]

		dists = np.sqrt((x_image - x_center)**2 + (y_image - y_center)**2)
		mask = dists > radius
		data = data[mask,:]
		if outfile is None:
			outfile = infile
			print "orginal file %s will be overwrited"%infile

		M,N = data.shape
		np.savetxt(outfile, data, fmt='%6.2f '*(N-1)+'%6.2f')

#Relative calibration
	def get_relative_mags(self, tpl_obs_match_method = 'grmatch', updatetable=1):
		'''
		relative magnitude calibration
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		flts_unique = self.templates.keys()
		offset_method = self.cal_offset_method
		for flt in flts_unique:
			print "Now working on %s band"%flt
			self.get_relative_mags_flt_general(flt, tpl_obs_match_method = tpl_obs_match_method, offset_method=offset_method, updatetable=updatetable)

	def get_relative_mags_flt_general(self,flt, tpl_obs_match_method = 'grmatch', offset_method='median', updatetable=1):

		if updatetable:
			self.__dict2table()
		self.__find_template_imagekey(updatetable=updatetable)

		tpl_imgkey = self.templates[flt]
		images_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)

		for image_key in images_flt['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if self.photometry_info[image_key]['relmag'] != 99.99:
				if not self.renew_relative_mag:
					continue
			self.__get_relative_mag_single(image_key,tpl_imgkey, tpl_obs_match_method = tpl_obs_match_method, offset_method=offset_method)

	def get_relative_mags_flt_auto_shift_method(self, flt, update_matched_inmags=True, refnumin=5, verbose=1):
		'''
		Relative calibration by manually select the objects for the calibration process. The magnitude table of reference objects on the reference image should be prepared before this.

		See also: self.get_template_reference_mags_flt()
		'''
		self.__dict2table()
		self.__find_template_imagekey()

		tpl_imgkey = self.templates[flt]
		tpl_imgkey_s = tpl_imgkey.split('.')[0]
		images_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)

		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		if not os.path.exists(refmag_filename):
			raise ValueError("please get the reference mags first from function get_template_reference_mags_flt")
		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')
		if verbose:
			print "mag table for reference stars on template image:"
			print refmag_table

		refstarnum = len(refmag_table)

		for image_key in images_flt['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if self.photometry_info[image_key]['relmag'] != 99.99:
				if not self.renew_relative_mag:
					continue
			if image_key == tpl_imgkey:
				self.photometry_info[image_key]['relmag'] = self.photometry_info[image_key]['instmag']
				self.photometry_info[image_key]['relmagerr'] = 0.0
				continue

			relmag, relmagerr = self.__get_relative_mag_single_simple_shift_and_match(image_key, refmag_table=refmag_table, update_matched_inmags=update_matched_inmags,refnumin=refnumin, verbose=verbose)
			self.photometry_info[image_key]['relmag'] = relmag
			self.photometry_info[image_key]['relmagerr'] = relmagerr

	def get_relative_mags_flt_manual_match_method(self, flt, which_dir='raw_image', verbose=1, update_matched_inmags=True):
		'''
		Relative calibration by manually select the objects for the calibration process. The magnitude table of reference objects on the reference image should be prepared before this.

		See also: self.get_template_reference_mags_flt()
		'''
		self.__dict2table()
		self.__find_template_imagekey()

		tpl_imgkey = self.templates[flt]
		tpl_imgkey_s = tpl_imgkey.split('.')[0]
		images_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)

		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		if not os.path.exists(refmag_filename):
			raise ValueError("please get the reference mags first from function get_template_reference_mags_flt")
		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')
		if verbose:
			print "mag table for reference stars on template image:"
			print refmag_table

		refstarnum = len(refmag_table)

		for image_key in images_flt['name']:
			print image_key
			if self.photometry_info[image_key]['drop'] > 0:
				continue

			if self.photometry_info[image_key]['relmag'] != 99.99:
				if not self.renew_relative_mag:
					continue

			if image_key == tpl_imgkey:
				self.photometry_info[image_key]['relmag'] = self.photometry_info[image_key]['instmag']
				self.photometry_info[image_key]['relmagerr'] = 0.0
				continue

			relmag, relmagerr = self.__get_relative_mag_single_manual_match(image_key, flt, refmag_table=refmag_table, which_dir=which_dir, verbose=verbose, refstarnum=refstarnum, update_matched_inmags=update_matched_inmags)
			self.photometry_info[image_key]['relmag'] = relmag
			self.photometry_info[image_key]['relmagerr'] = relmagerr

	def __get_tpl_calibrated_mag_manual_method(self, tpl_imgkey, flt, which_dir ='raw_image'):
		'''
		select the standard stars which match the reference stars on the template image and perform the calibration

		What required before this:
			self.get_template_reference_mags_flt
			self.show_standards_for_wcs_not_available
		'''
		input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)
		tpl_skey = tpl_imgkey.split('.')[0]

		std_xymags_file = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.txt"%(flt, tpl_skey))
		if not os.path.exists(std_xymags_file):
			raise IOError("%s not exists; please get it ready first from function show_standards_for_wcs_not_available"%std_xymags_file)

		xysmags = np.loadtxt(std_xymags_file)
		regionfile = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.reg"%(flt, tpl_skey))
		stdcal_idxymag_table = self.input_xys_through_ds9_get_mags(input_image, xysmags, regionfile=regionfile)
		stdcalmag_filename  = os.path.join(self.template_dir, "%s_%s_std_idxymag.txt"%(flt, tpl_skey))
		stdcal_idxymag_table.write(stdcalmag_filename, format='ascii.fixed_width')

		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_skey))
		if not os.path.exists(refmag_filename):
			raise ValueError("please get the reference mags first; the function will be used get_template_reference_mags_flt ")
		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')

		stdmags	= []
		stdmagerrs = []
		tplmags = []
		tplmagerrs = []

		for stdstar in stdcal_idxymag_table:
			sid = stdstar['num']
			if sid not in refmag_table['num']:
				raise ValueError("reference star with id number %s not found in reference table"%sid)
			tplstar = refmag_table[refmag_table['num']==sid][0]

			stdmags.append(stdstar['mag'][0])
			stdmagerrs.append(stdstar['magerr'][0])
			tplmags.append(tplstar['mag'])
			tplmagerrs.append(tplstar['magerr'])

		print stdmags
		mags_ref_matched = np.transpose(np.array([stdmags, stdmagerrs]))
		mags_input_matched = np.transpose(np.array([tplmags, tplmagerrs]))
		print mags_ref_matched
		print mags_input_matched

		tpl_instmag = self.photometry_info[tpl_imgkey]['instmag']
		tpl_instmagerr = self.photometry_info[tpl_imgkey]['instmagerr']

		cal_mag,cal_magerr, temp_future_use = self.__relative_calibration_single(mags_ref_matched,mags_input_matched,tpl_instmag)
		self.photometry_info[tpl_imgkey]['calmag'] = cal_mag
		self.photometry_info[tpl_imgkey]['calmagerr'] =  np.sqrt( cal_magerr**2 + tpl_instmagerr**2 )

	def __get_offset_between_given_refmags_inmags(self, refmags, inmags, instmag, average_method='weighted', verbose=1):
		'''
		get the offset between refmags and inmags
		'''
		offsets = np.array(refmags[:,0]) - np.array(inmags[:,0])
		offsetserr = np.sqrt(np.array(refmags[:,1])**2 + np.array(inmags[:,1])**2)
		print "The offsets and uncertainties on offsets are:"
		print offsets
		print offsetserr
		offsetserr[np.where(offsetserr==0.0)] = 0.01
		print offsetserr
		if average_method == 'weighted':
			weights = 1./offsetserr**2
			offset_mean = np.sum(weights*offsets)/ np.sum(weights)
			offseterr_ave = 1./ np.sqrt(np.sum(weights))
		elif average_method == 'simple':
			offset_mean = np.mean(offsets)
			offseterr_ave = 0.0
		else:
			raise IOError("invalid input for average method")

		relmag = instmag + offset_mean
		relmagerr = np.max([np.std(offsets), offseterr_ave])  #To be conservative with the uncertainty. Attention here!!
		if verbose:
			print "offset=%s+/-%s"%(offset_mean,relmagerr)

		return relmag, relmagerr

	def __match_refmags_and_inmags_by_id(self, refmag_table, inmag_table):
		'''
		refmag_table and inmag_table have columns of 'num','mag','magerr'
		'''
		refmags	= []
		refmagerrs = []
		inmags = []
		inmagerrs = []
		for instar in inmag_table:
			sid = instar['num']
			if sid not in refmag_table['num']:
				raise ValueError("reference star with id number %s not found in reference table"%sid)
			refstar = refmag_table[refmag_table['num']==sid][0]

			if instar['mag'] == 99.99:
				continue
			inmags.append(instar['mag'])
			inmagerrs.append(instar['magerr'])
			refmags.append(refstar['mag'])
			refmagerrs.append(refstar['magerr'])

		refmags_match = np.transpose(np.array([refmags, refmagerrs]))
		inmags_match  = np.transpose(np.array([inmags, inmagerrs]))

		return refmags_match, inmags_match

	def __find_corresponding_refmags_for_given_imgkey(self, image_key=None, flt=None, verbose=1):
		'''
		For given imgkey, what it the reference mags table?
		First find the template image with tpl_imgkey_s in the same band with the input image, then the refmag_table is [flt]_[tpl_imgkey_s]_idxymag.txt
		'''
		if flt is None:
			if image_key is None:
				raise ValueError('None for both image_key and flt, really?! Give one of them pls')
			flt = self.photometry_info[image_key]['flt']
		try:
			if len(self.templates) == 0:
				self.__find_template_imagekey()
			if flt not in self.templates.keys():
				raise IOError("template image for %s not available"%flt)
			tpl_imgkey = self.templates[flt]
			tpl_imgkey_s = tpl_imgkey.split('.')[0]
			refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
			if not os.path.exists(refmag_filename):
				raise ValueError("please get the reference mags first from function get_template_reference_mags_flt")
			refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')
			if verbose:
				print "mag table for reference stars on template image:"
				print refmag_table
		except Exception as e:
			print "failed to find the magnitude table for the reference image"
			raise ValueError(e)

		return refmag_table


	def __get_idxys_corresponding_to_reference_image(self, image_key,refmag_table):
		'''
		INPUTS:
			refmag_table: columns with name 'num','x', and 'y' are Required
		'''
		flt = self.photometry_info[image_key]['flt']
		if len(self.templates) == 0:
			self.__find_template_imagekey()
		if flt not in self.templates.keys():
			raise IOError("template image for %s not available"%flt)
		tpl_imgkey = self.templates[flt]

		snx = self.photometry_info[image_key]['x']
		sny = self.photometry_info[image_key]['y']
		if snx == self.photometry_info_init_values['x'] or sny == self.photometry_info_init_values['y']:
			raise ValueError('To use this method, the position for the source of interest on the input image should be provided')

		snxref = self.photometry_info[tpl_imgkey]['x']
		snyref = self.photometry_info[tpl_imgkey]['y']
		if snxref == self.photometry_info_init_values['x'] or snyref == self.photometry_info_init_values['y']:
			raise ValueError('To use this method, the position for the source of interest on the template image should be provided')

		xshift = snx - snxref
		yshift = sny - snyref
		idxys = np.array([refmag_table['num'], refmag_table['x']+xshift, refmag_table['y']+yshift]).transpose()

		return idxys

	def __get_relative_mag_single_simple_shift_and_match(self, image_key,refmag_table=None, update_matched_inmags = False, refnumin=5, average_method='weighted', update_rettable=False, verbose=1):
		'''
		initilly write this for SMARTS NIR calibration when the images have the same orientation and simple shift can match the sources
		'''
		image_key_s = image_key.split('.')[0]
		if refmag_table is None:
			print "The refmag_table is not provided. Try to find it."
			refmag_table = self.__find_corresponding_refmags_for_given_imgkey(image_key=image_key, verbose=verbose)

		idxys = self.__get_idxys_corresponding_to_reference_image(image_key,refmag_table)
		inmag_table = self.__relative_calibration_prepare_autoshift_matching(image_key, idxys, clobber=update_matched_inmags, verbose=verbose)
		if self.photometry_method   == 'apphot':
			magtable_savedir = self.aperture_photometry_dir
		elif self.photometry_method == 'psfphot':
			magtable_savedir = self.__get_internal_psfphot_dir()
		else:
			raise IOError("Invalid input for photometry_method")
		savefile = os.path.join(magtable_savedir, 'prematch_'+image_key_s+'.txt')
		inmag_table.write(savefile, format='ascii.fixed_width')

		if verbose:
			print "mag table for reference stars on input image:"
			print inmag_table

		refmags_match, inmags_match = self.__match_refmags_and_inmags_by_id(refmag_table, inmag_table)
		instmag = self.photometry_info[image_key]['instmag']

		if len(refmags_match)<refnumin:
			if verbose:
				print "less than %s reference stars will be used for the relative calibration."%refnumin
				print "The %s average of the offsets will be applied and std of the offsets will be used as uncertainty"%average_method
			relmag, relmagerr = self.__get_offset_between_given_refmags_inmags(refmags_match, inmags_match, instmag, average_method=average_method, verbose=verbose)
		else:
			relmag, relmagerr, temp_future_use = self.__relative_calibration_single(refmags_match, inmags_match, instmag)

		if verbose:
			print "relmag=%s, relmagerr=%s"%(relmag, relmagerr)

		if update_rettable:
			self.photometry_info[image_key]['relmag'] = relmag
			self.photometry_info[image_key]['relmagerr'] = relmagerr

		return relmag, relmagerr

	def __get_relative_mag_single_manual_match(self, image_key, refmag_table=None, which_dir ='raw_image', verbose=1, update_matched_inmags = False, refstarnum=None,  refnumin=5, average_method='weighted', update_rettable=False):
		'''
		relative calibration beween template image and input image. the matching is done by interactively selecting
		INPUTS:
			image_key:
			flt:
			which_dir:
			verbose:
			update_matched_inmags: update the matched mag list if it already exists
			refstarnum: the number of manually seletected reference star
			refnumin: the minimum number of reference stars required for the method of self.__relative_calibration_single; otherwise simple weighted average will be used
			average_method: 'weighted' or 'simple'
		'''
		image_key_s = image_key.split('.')[0]

		if refmag_table is None:
			print "The refmag_table is not provided. Try to find it."
			refmag_table = self.__find_corresponding_refmags_for_given_imgkey(image_key=image_key)

		inmag_table = self.__relative_calibration_prepare_interactive_matching(image_key, which_dir=which_dir, clobber=update_matched_inmags, refstarnum=refstarnum)
		if self.photometry_method   == 'apphot':
			magtable_savedir = self.aperture_photometry_dir
		elif self.photometry_method == 'psfphot':
			magtable_savedir = self.__get_internal_psfphot_dir()
		else:
			raise IOError("Invalid input for photometry_method")
		savefile = os.path.join(magtable_savedir, 'prematch_'+image_key_s+'.txt')
		inmag_table.write(savefile, format='ascii.fixed_width')

		if verbose:
			print "mag table for reference stars on input image:"
			print inmag_table

		refmags_match, inmags_match = self.__match_refmags_and_inmags_by_id(refmag_table, inmag_table)
		instmag = self.photometry_info[image_key]['instmag']

		if len(refmags_match)<refnumin:
			if verbose:
				print "less than %s reference stars will be used for the relative calibration."%refnumin
				print "The %s average of the offsets will be applied and std of the offsets will be used as uncertainty"%average_method
			relmag, relmagerr = self.__get_offset_between_given_refmags_inmags(refmags_match, inmags_match, instmag, average_method=average_method, verbose=verbose)
		else:
			relmag, relmagerr, temp_future_use = self.__relative_calibration_single(refmags_match, inmags_match, instmag)

		if verbose:
			print "relmag=%s, relmagerr=%s"%(relmag, relmagerr)

		if update_rettable:
			self.photometry_info[image_key]['relmag'] = relmag
			self.photometry_info[image_key]['relmagerr'] = relmagerr

		return relmag, relmagerr

	def __get_relative_mag_single(self,image_key,tpl_imgkey,  tpl_obs_match_method = 'grmatch', offset_method='median'):
		'''
		relative calibration between template image and input image
		'''
		instmag = self.photometry_info[image_key]['instmag']
		instmagerr = self.photometry_info[image_key]['instmagerr']
		if image_key == tpl_imgkey:
			relmag = instmag
			relmagerr = 0
		else:
			mags_ref_matched,mags_input_matched = self.__relative_calibration_prepare(image_key,tpl_imgkey, tpl_obs_match_method = tpl_obs_match_method)

			if self.photometry_info[image_key]['drop']>0:
				print "image %s has been dropped... please check"%image_key
				return

			if self.photometry_method   == 'apphot':
				match_prereject_filename  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')
				match_caldata_filename  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl_caldata.match')
			elif self.photometry_method == 'psfphot':
				psf_photometry_dir = self.__get_internal_psfphot_dir()
				match_prereject_filename  = os.path.join(psf_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')
				match_caldata_filename  = os.path.join(psf_photometry_dir,image_key.split('.')[0]+ '_tpl_caldata.match')
			else:
				raise IOError("Invalid input for photometry_method")

			relmag,relmagerr, caldata_index  = self.__relative_calibration_single(mags_ref_matched,mags_input_matched,instmag, offset_method=offset_method)
			match_alldata = np.loadtxt(match_prereject_filename)
			match_caldata = match_alldata[caldata_index,:]
			np.savetxt(match_caldata_filename, match_caldata, fmt="%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f")

            	print relmag,relmagerr
                self.photometry_info[image_key]['relmag'] = relmag
                self.photometry_info[image_key]['relmagerr'] = relmagerr

	def __relative_calibration_prepare_autoshift_matching(self, image_key, idxys, clobber=True, verbose=1):
		'''
		prepare magnitudes list for relative calibration: input magnitudes with corresponding id against the reference mags

		INPUTs:
			image_key: input image key, eg. 000.fits
			idxys:

		PRODUCTS:
			matched list between input magnitudes and reference magnitudes, i.e 'xxx_tpl.match'
		'''
		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		image_key_s = image_key.split('.')[0]
		if self.photometry_method   == 'apphot':
			input_list_filename    = os.path.join(self.aperture_photometry_dir, image_key_s +self.apphot_ret_file_suffix)
			match_output_filename  = os.path.join(self.aperture_photometry_dir,image_key_s + '_tpl.match')
			photdir = self.aperture_photometry_dir
		elif self.photometry_method == 'psfphot':
			input_list_filename    = self.__get_internal_psfphot_mag_file(image_key)
			photdir = self.__get_internal_psfphot_dir()
			match_output_filename  = os.path.join(photdir,image_key_s + '_tpl.match')
		else:
			raise IOError("Invalid input for photometry_method")

		if clobber:
			self.__delete_file_if_exist(match_output_filename)

		if not os.path.exists(match_output_filename):
			if not os.path.exists(input_list_filename):
				raise IOError("Input photometry result file %s not exists"%input_list_filename)
			inputphot_ret = np.loadtxt(input_list_filename)
			calmatch_table = self.__input_xys_get_mags(inputphot_ret, idxys, verbose=verbose)
			calmatch_table.write(match_output_filename, format='ascii.fixed_width')
		else:
			if verbose:
				print "matched file between template image and input image %s already exist"%match_output_filename

		calmatch_table = Table.read(match_output_filename, format='ascii.fixed_width') #the table read from previouly save table

		return calmatch_table

	def __relative_calibration_prepare_interactive_matching(self,image_key, which_dir='raw_image', clobber=True, refstarnum=None, verbose=1):
		'''
		prepare magnitudes list for relative calibration: input magnitudes with corresponding id against the reference mags

		INPUTs:
			image_key: input image key, eg. 000.fits
			refstarnum: If prompt for input of id number is not expected, make sure to provide the refstarnum (for example N), which will give id number from 0 to N-1 by default.
		OUTPUTS:
			mags_input_matched:

		PRODUCTS:
			matched list between input magnitudes and reference magnitudes, i.e 'xxx_tpl.match'
		'''
		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		image_key_s = image_key.split('.')[0]

		if self.photometry_method   == 'apphot':
			input_list_filename    = os.path.join(self.aperture_photometry_dir, image_key_s +self.apphot_ret_file_suffix)
			match_output_filename  = os.path.join(self.aperture_photometry_dir,image_key_s + '_tpl.match')
			photdir = self.aperture_photometry_dir
		elif self.photometry_method == 'psfphot':
			input_list_filename    =  self.__get_internal_psfphot_mag_file(image_key)
			photdir = self.__get_internal_psfphot_dir()
			match_output_filename  = os.path.join(photdir, image_key_s + '_tpl.match')
		else:
			raise IOError("Invalid input for photometry_method")

		if clobber:
			self.__delete_file_if_exist(match_output_filename)

		if not os.path.exists(match_output_filename):
			if not os.path.exists(input_list_filename):
				raise IOError("Input photometry result file %s not exists"%input_list_filename)

			inputphot_ret = np.loadtxt(input_list_filename)
			xys = inputphot_ret[:,0:2]
			mags = inputphot_ret[:,2]
			magerrs = inputphot_ret[:,3]

			input_image = self.__get_internal_image(image_key, which_dir=which_dir)
			psf_photometry_dir = self.__get_internal_psfphot_dir()
			regionfile = os.path.join(psf_photometry_dir, image_key_s +"_xymag.reg")
			if not os.path.exists(regionfile):
				create_ds9_region_file(inputphot_ret, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)
			calmatch_table = self.input_xys_through_ds9_get_mags(input_image, inputphot_ret, regionfile=regionfile, wantnum=refstarnum)
			calmatch_table.write(match_output_filename, format='ascii.fixed_width')
		else:
			if verbose:
				print "matched file between template image and input image %s already exist"%match_output_filename

		calmatch_table = Table.read(match_output_filename, format='ascii.fixed_width') #the table read from previouly save table

		return calmatch_table

	def __relative_calibration_prepare(self,image_key,tpl_imgkey, remove_transient = True, tpl_obs_match_method = 'grmatch'):
		'''
		prepare magnitudes list for relative calibration: input magnitudes and reference magnitudes

		INPUTs:
			image_key: input image key, eg. 000.fits
			tpl_imgkey: reference image key, eg. 001.fits
			remove_transient: remove the transient from the matched result of reflist and input list
			tpl_obs_match_method: 'grmatch' or 'surrounding_search'

		OUTPUTS:
			mags_ref_matched:
			mags_input_matched:

		PRODUCTS:
			matched list between input magnitudes and reference magnitudes, i.e 'xxx_tpl.match'
			transform coefficient, i.e 'xxx_tpl.coef'

		Notes:
			the matched magnitudes will be filterd according to self.reference_mag_min and self.reference_mag_max

		'''
		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		if self.photometry_method   == 'apphot':
			ref_list_filename      = os.path.join(self.aperture_photometry_dir,tpl_imgkey.split('.')[0]+self.apphot_ret_file_suffix)
			input_list_filename    = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0]+self.apphot_ret_file_suffix)
			match_output_filename  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl.match')
			match_output_filename_tstrmd  = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')  # transient removed; too faint, too bright and too large uncertainty stars are rmeoved
			trans_fitting_filename = os.path.join(self.aperture_photometry_dir,image_key.split('.')[0] + '_tpl.coef')
		elif self.photometry_method == 'psfphot':
			ref_list_filename      = self.__get_internal_psfphot_mag_file(tpl_imgkey)
			input_list_filename    = self.__get_internal_psfphot_mag_file(image_key)
			psf_photometry_dir = self.__get_internal_psfphot_dir()
			match_output_filename  = os.path.join(psf_photometry_dir,image_key.split('.')[0]+ '_tpl.match')
			match_output_filename_tstrmd  = os.path.join(psf_photometry_dir,image_key.split('.')[0]+ '_tpl_tstrmd.match')  # transient removed; too faint, too bright and too large uncertainty stars are rmeoved
			trans_fitting_filename = os.path.join(psf_photometry_dir,image_key.split('.')[0] + '_tpl.coef')
		else:
			raise IOError("Invalid input for photometry_method")

		self.__delete_file_if_exist(match_output_filename)
		self.__delete_file_if_exist(trans_fitting_filename)

		if self.trim_before_grmatch:
			if self.trim_left is None or self.trim_right is None or self.trim_bottom is None or self.trim_up is None:
				raise IOError("trim boundary required!")

			ref_list_data = np.loadtxt(ref_list_filename)
			input_list_data = np.loadtxt(input_list_filename)

			print len(input_list_data)
			nl = self.trim_left
			nb = self.trim_bottom

			if self.trim_right <0 or self.trim_up <0:

				img_template = self.images[tpl_imgkey]
				hdu_tpl = fits.open(img_template)
				header_tpl = hdu_tpl['PRIMARY'].header
				try:
					NX_tpl = header_tpl['NAXIS1']
					NY_tpl = header_tpl['NAXIS2']
				except:
					print "image size not found in the tempalte image"


				img_input = self.images[image_key]
				hdu_input = fits.open(img_input)
				header_input = hdu_input['PRIMARY'].header
				try:
					NX_input = header_input['NAXIS1']
					NY_input = header_input['NAXIS2']
				except:
					print "image size not found in the input image"

				nr_tpl = NX_tpl + self.trim_right
				nu_tpl = NY_tpl + self.trim_up
				nr_input = NX_input + self.trim_right
				nu_input = NY_input + self.trim_up

			else:
				nr_tpl = self.trim_right
				nu_tpl = self.trim_up
				nr_input = self.trim_right
				nu_input = self.trim_up

			print "input image boundaries:",   nl,nr_input,nb,nu_input
			print "template image boundaries:",nl,nr_tpl,  nb,nu_tpl

			ref_list_data_trim_x = ref_list_data[np.logical_and(ref_list_data[:,0]>nl,ref_list_data[:,0]<nr_tpl),:]
			ref_list_data_trim = ref_list_data_trim_x[np.logical_and(ref_list_data_trim_x[:,1]>nb,ref_list_data_trim_x[:,1]<nu_tpl),:]

			input_list_data_trim_x = input_list_data[np.logical_and(input_list_data[:,0]>nl,input_list_data[:,0]<nr_input),:]
			input_list_data_trim = input_list_data_trim_x[np.logical_and(input_list_data_trim_x[:,1]>nb,input_list_data_trim_x[:,1]<nu_input),:]
			print len(input_list_data_trim)
			ref_list_filename = ref_list_filename + '.trim'
			input_list_filename = input_list_filename + '.trim'

			np.savetxt(ref_list_filename,  ref_list_data_trim)
			np.savetxt(input_list_filename,input_list_data_trim)

		if self.magnitude_cut_before_grmatch:
			if self.reference_mag_min is None and self.reference_mag_max is None:
				raise IOError("magnitude cut required... self.reference_mag_min or self.refence_mag_max or both are required")

			ref_list_data = np.loadtxt(ref_list_filename)
			input_list_data = np.loadtxt(input_list_filename)

			if self.reference_mag_min is not None:
				ref_list_data_magcut = ref_list_data[ref_list_data[:,2]<self.reference_mag_min]
				input_list_data_magcut = input_list_data[input_list_data[:,2]<self.reference_mag_min]
			else:
				ref_list_data_magcut = ref_list_data
				input_list_data_magcut = input_list_data

			if self.reference_mag_max is not None:
				ref_list_data_magcut = ref_list_data_magcut[ref_list_data_magcut[:,2]>self.reference_mag_max]
				input_list_data_magcut = input_list_data_magcut[input_list_data_magcut[:,2]>self.reference_mag_max]

			ref_list_filename = ref_list_filename + '.magcut'
			input_list_filename = input_list_filename + '.magcut'

			np.savetxt(ref_list_filename,  ref_list_data_magcut)
			np.savetxt(input_list_filename,input_list_data_magcut)

		visbfmatch = self.tpl_obs_before_match_display #show the data before matching
		if visbfmatch:
			refdata_plot = np.loadtxt(ref_list_filename)
			inputdata_plot = np.loadtxt(input_list_filename)
			xy_simple_scatter_plot_2sets(refdata_plot[:,0], refdata_plot[:,1], inputdata_plot[:,0], inputdata_plot[:,1], label1='reference', label2='input')

		if tpl_obs_match_method == 'grmatch':
			self.fitsh_grmatch_type = 'point'
			matched_table = self.__fitsh_grmatch(ref_list_filename, input_list_filename, match_output_filename, trans_output=trans_fitting_filename)
		elif tpl_obs_match_method == 'surrounding_search':
			input_list_data = np.loadtxt(input_list_filename)
			ref_list_data   = np.loadtxt(ref_list_filename)
			tg_xys  = input_list_data[:,0:2]
			ref_xys = ref_list_data[:,0:2]
			criteria = self.criteria_tpl_obs_match
			target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)

			if len(indexmat) != 0:
				match_ret = np.hstack((ref_list_data[indexref],input_list_data[indexmat]))
				np.savetxt(match_output_filename,match_ret)

		match_result_plot = self.tpl_obs_match_result_display
		if match_result_plot:
			if tpl_obs_match_method == 'grmatch':
				if not os.path.exists(match_output_filename):
					raise IOError("grmatch failure...")

				if self.photometry_method == 'apphot':
					trans_output = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0] + '_tpl.trans')
				elif self.photometry_method == 'psfphot':
					psf_photometry_dir = self.__get_internal_psfphot_dir()
					trans_output = os.path.join(psf_photometry_dir, image_key.split('.')[0] + '_tpl.trans')
				else:
					raise IOError("Invalid input for photometry_method")

				self.__fitsh_grtrans(match_output_filename, trans_output, trans_fitting_filename)
				trans_result = np.loadtxt(trans_output)
				ref_xys = trans_result[:,[0,1]]
				input_xys  = trans_result[:,[4,5]]
				xy_simple_scatter_plot(ref_xys[:,0]-input_xys[:,0], ref_xys[:,1]-input_xys[:,1], xlabel='x_ref - x_input', ylabel='y_ref - y_input')
			elif tpl_obs_match_method == 'surrounding_search':
				if os.path.exists(match_output_filename):
					match_result = np.loadtxt(match_output_filename)
					ref_xys   = match_result[:,[0,1]]
					input_xys = match_result[:,[4,5]]
					xy_simple_scatter_plot_2sets(ref_xys[:,0], ref_xys[:,1], input_xys[:,0], input_xys[:,1], label1='reference', label2='input')
				else:
					print "no valid match obtained..."
		try:
			matched_ret   = np.loadtxt(match_output_filename)

			if self.reference_mag_min is not None:
				mag_faint_cut = self.reference_mag_min
				matched_ret = matched_ret[matched_ret[:,2]<mag_faint_cut,:]

			if self.reference_mag_max is not None:
				mag_bright_cut  = self.reference_mag_max
				matched_ret = matched_ret[matched_ret[:,2]>mag_bright_cut,:]

			if self.reference_magerr_cut is not None:
				magerr_cut_ref = self.reference_magerr_cut
				matched_ret = matched_ret[matched_ret[:,3]<magerr_cut_ref,:]

			if self.input_magerr_cut is not None:
				magerr_cut_input = self.input_magerr_cut
				matched_ret = matched_ret[matched_ret[:,7]<magerr_cut_input,:]

			#get rid of the transient
			if remove_transient:
				x = self.photometry_info[tpl_imgkey]['x']
				y = self.photometry_info[tpl_imgkey]['y']

				transient_exist, corresp_xy, index_transient = self.__find_corresponding(matched_ret[:,0:2], [x,y], 2)
				if transient_exist:
					indexs = range(len(matched_ret))
					for index_remove in index_transient:
						indexs.remove(index_remove)
					matched_ret = matched_ret[indexs,:]
				np.savetxt(match_output_filename_tstrmd, matched_ret, fmt="%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f")
			ref_matched   = matched_ret[:,[0,1,2,3]]
			input_matched = matched_ret[:,[4,5,6,7]]
		except:
			self.photometry_info[image_key]['drop'] = 5
			return 0,0

		xys_ref_matched   = ref_matched[:,[0,1]]
		xys_input_matched = input_matched[:,[0,1]]
		mags_ref_matched = ref_matched[:,[2,3]]
		mags_input_matched = input_matched[:,[2,3]]

		return mags_ref_matched,mags_input_matched

	def __relative_calibration_single(self,ref_mags,input_mags,instrument_mag, xdata=None, offset_method='median',  offset_err_stat_binsize=0.1, offset_err_stat_xstepsize=0.5, sigclipping_first  = True, title_text=None):
		'''
		updates 2022-03-18: 1) improvement on the uncertainty, from single uncertainty for all brightness to target brightness dependent uncertainty
				    2) self.cal_offset_funcfit_type add new funcfit type of 'o2poly' to account for the non-linearity for example in LCOGT0.4m camera https://s3.us-west-2.amazonaws.com/www.lco.global/documents/Photometric_nonlinearity_SBIG_6303_cameras.20220310.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA6FT4CXR4ZJRYWHNN%2F20220318%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20220318T134518Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&X-Amz-Signature=a88fbbb7e902662fd520c59f80b1db273f25a05f37de9fc7ea8dbfcdce17fcfe
		INPUTS:
			ref_mags:
			input_mags:
			instrument_mag:
			xdata: default is ref_mags
			offset_method: median, mean, or funcfit.
			offset_err_stat_binsize: if <1 then means fraction of total number, otherwise the absolute value 
			offset_err_stat_xstepsize: 
			sigclipping_first: simple sigma clipping before any fitting of the offset
			title_text: if self.cal_plot, the title of the plot
		'''
		if xdata is None:
			xdata = ref_mags

		if sigclipping_first:
			offsets_first = ref_mags[:,0] - input_mags[:,0]
			survive_mask,discard_mask = self.__sigma_clipping(offsets_first, sig=5, variability_method='mad')
			mags_ref = ref_mags[survive_mask,0]
			mags_input = input_mags[survive_mask,0]
			magserr_ref = ref_mags[survive_mask,1]
			magserr_input = input_mags[survive_mask,1]
			xs = xdata[survive_mask,0]
			xerrs = xdata[survive_mask,1]
		else:
			mags_ref = ref_mags[:,0]
			mags_input = input_mags[:,0]
			magserr_ref = ref_mags[:,1]
			magserr_input = input_mags[:,1]
			survive_mask = np.arange(len(xdata))
			xs = xdata[:,0]
			xerrs = xdata[:,1]

		mags_offset = mags_ref - mags_input
		magserr_offset = np.sqrt(magserr_ref**2+magserr_input**2)

		sortids = np.argsort(xs)	
		xs = xs[sortids]
		xerrs = xerrs[sortids]
		mags_offset = mags_offset[sortids]
		magserr_offset = magserr_offset[sortids]


		if offset_method == "funcfit":
			magoffset0 = np.median(mags_offset)
			if self.cal_offset_funcfit_type == 'o1poly':
				def func(x,a,b):
					return a+b*x
				def yerr(x,aerr,berr):
					return np.sqrt((x*berr)**2+aerr**2)
				def func_str(popt): #give the function form to display on the final plot  a+b*x
					a=np.round(popt[0], 3)
					b=np.round(popt[1], 3)
					return '%s+%sx'%(a,b)
				
				p0 = np.array([magoffset0, 0])
			elif self.cal_offset_funcfit_type == 'o2poly':
				def func(x,a,b,c):
					return a+b*x+c*x**2
				def yerr(x,aerr,berr):
					return np.sqrt((x**2*cerr)**2+(x*berr)**2+aerr**2)
				def func_str(popt): #give the function form to display on the final plot  n(x-x0)**2+q
					a, b, c = popt
					n = str(np.round(c,5))
					x0 = str(np.round(-b/c/2, 3))
					q = str(np.round(a-b**2/c/4,3))
					return '%s(x-%s)^2+%s'%(n,x0,q)
				p0 = np.array([magoffset0, 0, 0])
			elif self.cal_offset_funcfit_type == 'constant':
				def func(x,c):
					return x*0+c
				def yerr(x,cerr):
					return x*0+cerr
				def func_str(popt):
					return "offset=%s"%str(np.round(popt[0],3))
				p0 = np.array([magoffset0])
			else:
				raise ValueError('%s for self.cal_offset_funcfit_type not supported yet'%self.cal_offset_funcfit_type)

			nsig_clip = self.cal_offset_nsig_clip
			fitdata, popt, perr = self.__funcfit_remove_outlier(xs,mags_offset,magserr_offset,func, p0, nsig=nsig_clip, rm_mode = 'all')
			offset = func(fitdata[:,0],*popt)
			x0 = instrument_mag
			xfinal = instrument_mag + func(x0, *popt)
			while np.abs(xfinal-x0)>0.001:
				x0 = xfinal
				xfinal = instrument_mag + func(x0, *popt)
			offset_ret = func(xfinal, *popt)
		elif offset_method in ['median', 'mean']: #this should be consider as deprecated after implenting the function type of 'constant' for the cal_offset_funcfit_type
			N = len(xs)
			fitdata = np.hstack((xs.reshape((N,1)), mags_offset.reshape((N,1)), magserr_offset.reshape((N,1))))
			def func(x,c):
				return x*0+c
			def yerr(x,cerr):
				return x*0+cerr
			def func_str(popt):
				return "offset=%s"%str(np.round(popt[0],3))
			if offset_method == 'median':
				c =  np.median(mags_offset)
			else:
				c =  np.mean(mags_offset)
			
			popt = [c]
			offset = func(fitdata[:,0], *popt)
			offset_ret = func(instrument_mag, *popt)

		else:
			raise IOError("invalid input for offset_method...")


		xmin = np.min(fitdata[:,0])
		xmax = np.max(fitdata[:,0])
		offset_dev = np.abs(fitdata[:,1]-offset) #deviation from the offset model
		if self.cal_offset_err_method == 'global_std':
			cut = int(np.floor((len(fitdata)*0.688)))
			indice = np.argsort(offset_dev)[cut]
			offset_err = offset_dev[indice]
			def err_model_func(x,c):
				return x*0+c
			err_stat_popt = [offset_err]
			plot_ylim = 10*offset_err
		elif self.cal_offset_err_method == 'brightness_dependent_fit':
			Ntotal = len(fitdata)
			if offset_err_stat_binsize<1:
				Nbin = int(Ntotal*offset_err_stat_binsize)
			else:
				Nbin = int(offset_err_stat_binsize)

			halfNbin= int(Nbin/2.0)
			xstart =  fitdata[halfNbin,0]
			xend = 	fitdata[-(halfNbin+1),0]
			xs_init = np.arange(xstart, xend, offset_err_stat_xstepsize) #first guess on the xs 
			if xend > xs_init[-1]:
				xs_init = np.append(xs_init, xend)
			xs_residual = []
			ys_residual = [] #offset residual 'std' 
			for x in xs_init:
				xind = np.argsort(np.abs(fitdata[:,0]-x))[0] 
				xs_temp = fitdata[(xind-halfNbin):(xind+halfNbin), 0]
				ys_temp = offset_dev[(xind-halfNbin):(xind+halfNbin)]
				xs_residual.append(np.mean(xs_temp))
				ys_residual.append(1.48*np.median(np.abs(ys_temp-np.median(ys_temp))))
			xs_residual = np.array(xs_residual)
			ys_residual = np.array(ys_residual)
			
			def err_model_func(x,a,b,c):
				return a+b*x+c*x**2
			err_stat_popt, pcov = curve_fit(err_model_func, xs_residual-xmin, ys_residual, p0=[0, 0, 0], sigma=None)
			offset_err = err_model_func(instrument_mag+offset_ret-xmin, *err_stat_popt)
			plot_ylim = err_model_func(np.max(fitdata[:,0])-xmin, *err_stat_popt)*6
		else:
			raise ValueError("%s not supported"%self.cal_offset_err_method)

		if self.cal_plot:
			from matplotlib import gridspec

			fig = plt.figure(figsize=(9, 6))
			gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
			ax0 = plt.subplot(gs[0])
			if title_text is not None:
				ax0.set_title(title_text, fontsize=12)
			plt.xlabel('reference magnitudes')
			plt.ylabel('magnitude offset (reference - input)')

			ax0.errorbar(xs, mags_offset, yerr=magserr_offset, fmt='gx', alpha=0.3, label='funcfit input data')

			if sigclipping_first:
				xs_discard = xdata[discard_mask,0]
				offset_discard_mag = ref_mags[discard_mask,0] - input_mags[discard_mask,0]
				offset_discard_magerr = np.sqrt(ref_mags[discard_mask,1]**2 + input_mags[discard_mask,1]**2)
				if len(xs_discard) != 0:
					ax0.errorbar(xs_discard, offset_discard_mag, yerr=offset_discard_magerr,fmt= 'g.', alpha=0.3, label='pre-funcfit discarded (5 sigma)')

			ax0.errorbar(fitdata[:,0],fitdata[:,1],yerr=fitdata[:,2], fmt='r.', alpha=0.3, label='final fit data')

			fitret_label = func_str(popt)
			ax0.plot(fitdata[:,0], offset,'k', label=fitret_label)
			if self.cal_offset_err_method == "brightness_dependent_fit":
				ax0.plot(xs_residual, func(xs_residual,*popt)+ys_residual, 'ko', ms=10)
				ax0.plot(xs_residual, func(xs_residual,*popt)-ys_residual, 'ko', ms=10)
			ax0.plot(fitdata[:,0], func(fitdata[:,0],*popt)+ err_model_func(fitdata[:,0]-xmin,*err_stat_popt), 'k')
			ax0.plot(fitdata[:,0], func(fitdata[:,0],*popt)- err_model_func(fitdata[:,0]-xmin,*err_stat_popt), 'k')
			
			ax1 = plt.subplot(gs[1])
			hist,edges = np.histogram(fitdata[:,1]-offset,bins=15)
			bin_centers = edges[:-1] + (edges[1:]-edges[:-1])/2.
			ax1.step(hist,bin_centers)

			ax0.legend(loc=0)

			if instrument_mag != 0:
				upper_limit = offset_ret + 5*offset_err
				lower_limit = offset_ret - 5*offset_err
				ax0.plot([instrument_mag+offset_ret,instrument_mag+offset_ret],[lower_limit,upper_limit])

			ax0.set_ylim([offset_ret-plot_ylim, offset_ret+plot_ylim])
			ax0.set_xlim([xmin, xmax])
			plt.show()

		relative_mag = instrument_mag + offset_ret
		relative_magerr = offset_err
		caldata_index = survive_mask  #the index of input data which survive the first sig_clipping

		return relative_mag,relative_magerr, caldata_index

	def __sigma_clipping(self,input_data, sig=3, meanfunc=np.median, variability_method="stddev", verbose=0):
		"""
		Identify outliers considering the mean (if meanfunc=np.mean) or median (if meanfunc=np.median) value and 3 sigma (3*stdev),
		iterating until convergence.

		20180318: add new parameter variability_method which is used to input the choice of how to measure the variability of a univariate sample of quantitative data.
			  Two options are supported. If 'stddev', stddev =  np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
						     If 'mad', stddev = 1.48*np.median(np.abs(data-np.median(data)))

		INPUTS:
			input_data:

		"""
		data = input_data.copy()
		last_total = len(data)

		if variability_method == "stddev":
			stdev = np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
		elif variability_method == "mad":
			stdev = 1.48*np.median(np.abs(data-np.median(data)))
		else:
			raise IOError("invalid input for variability_method")

		if verbose:
			print "The stddev in the first iteration:", stdev

		diff = data - meanfunc(data)
		sfilter = np.abs(diff) < sig*stdev
		current_total = len(data[sfilter])
		# Continue iterating until convergence (no more points are removed)
		i = 1
		while last_total > current_total:
			data = data[sfilter]
			last_total = current_total

			if variability_method == "stddev":
				stdev = np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
			elif variability_method == "mad":
				stdev = 1.48*np.median(np.abs(data-np.median(data)))
			else:
				raise IOError("invalid input for variability_method")

			if verbose:
				i = i+1
				print "The stddev in iteration %s:"%i, stdev

			diff = data - meanfunc(data)
			sfilter = np.abs(diff) < sig*stdev
			current_total = len(data[sfilter])

		if variability_method == "stddev":
			stdev = np.sqrt(np.sum((data-meanfunc(data))**2)/last_total)
		elif variability_method == "mad":
			stdev = 1.48*np.median(np.abs(data-np.median(data)))
		else:
			raise IOError("invalid input for variability_method")

		meanv = meanfunc(data)
		diff = input_data - meanv
		survive_mask = np.abs(diff) < sig*stdev
		discard_mask = 	np.abs(diff) > sig*stdev

		return survive_mask,discard_mask

	def __funcfit_remove_outlier(self, x,y,yerr,func,p0, nsig=3, criteria_mode= 'individual', rm_mode = 'all', stddev_method = 'mad'):
		"""
		fitting the func to input data (x, y+/-yerr) and get the best fitting parameter.
		During the fitting process, the iterative rejection on the input data will be applied

		INPUT:
			x:
			y:
			yerr:
			func: function object takes x and parameter list and output y
			p0:
			nsig:
			criteria_mode: 'whole': the fixed single criteria for the sample; 'indivual': yerr value used for each point in criteria
			rm_mode: 'all': remove all beyond the criteria for one iteration, 'worst': only remove the worst for one iteration
			stddev_method: 'direct' or 'mad';  np.std() for direct, and 1.48*MAD() for mad
		"""
		def get_stddev(y, stddev_method):
			if stddev_method == 'direct':
				stddev = np.std(y)
			elif stddev_method == 'mad':
				stddev = np.median(np.abs(y-np.median(y)))
			else:
				raise ValueError("invalid input for stddev_method...")
			return stddev

		def remove_outliers(x, y, yerr, abs_diffs, nsig, stddev, rm_mode, criteria_mode):
			if self.cal_offset_verbose:
				print "find those have large deviation"
			xyyerr = np.array([x,y,yerr]).transpose()
			if criteria_mode == 'whole':
				sigma = np.ones(len(x))*stddev
			elif criteria_mode == 'individual':
				sigma = yerr
			else:
				raise ValueError('invalid input for criteria_mode')
			goodcriteria = nsig*sigma
			if rm_mode == 'all':
				goodflt = abs_diffs < goodcriteria
				x    = x[goodflt]
				y    = y[goodflt]
				yerr = yerr[goodflt]
				if self.cal_offset_verbose:
					print "the following are removed from fitting"
					print xyyerr[~goodflt]
			elif rm_mode == 'worst':
				indice = np.argsort(abs_diffs/sigma)
				if abs_diffs[indice[-1]] > goodcriteria[-1]:
					goodflt = indice[:-1]
					x    = x[goodflt]
					y    = y[goodflt]
					yerr = yerr[goodflt]
					if self.cal_offset_verbose:
						print "the following are removed from fitting"
						print xyyerr[indice[-1]]
						
			else:
				raise ValueError('invalid input for rm_mode...')
			return x,y,yerr

		popt, pcov = curve_fit(func, x, y, p0=p0, sigma=yerr)
		p0 = popt
		perr = np.sqrt(np.diag(pcov))
		y_model = func(x,*popt)
		abs_diffs = np.abs(y-y_model)
		last_total = len(y)

		stddev = get_stddev(y-y_model, stddev_method)
		x,y,yerr = remove_outliers(x,y,yerr,abs_diffs,nsig,stddev, rm_mode, criteria_mode)
		current_total = len(y)

		while last_total > current_total: #Continue iterating until convergence (no more points are removed)
			last_total = current_total
			popt, pcov = curve_fit(func, x, y, p0=p0, sigma=yerr)
			p0 = popt
			perr = np.sqrt(np.diag(pcov))
			y_model = func(x,*popt)
			abs_diffs = np.abs(y-y_model)

			stddev = get_stddev(y-y_model, stddev_method)
			x,y,yerr = remove_outliers(x,y,yerr,abs_diffs,nsig,stddev, rm_mode, criteria_mode)
			current_total = len(y)

		gooddata = np.hstack((x.reshape((len(x),1)), y.reshape((len(y),1)), yerr.reshape((len(yerr),1))))
		return gooddata,popt,perr

	def __crude_match(self,array1,array2,criteria,verbose=False):
	    '''
	    Conduct crude match between two given list of positions of stars

	    INPUTS:
	    array1: the first star array with first and second column array1[:,0:2] are star position(x,y) and array1.shape = (N1,m)
	    array2: the second star array with first and second column array1[:,0:2] are star position(x,y)and array1.shape = (N2,m)
	            array1 and array have the same width(columns) and ONE of N1 and N2 can equal 1
	    criteria: = c(PIXES) the criteria for match is that one star is within c pixels from the corrresponding star's central position
	    '''

	    if len(array1) < len(array2):
	        refarray = array2
	        matarray = array1
	        alter = 1
	    else:
	        refarray = array1
	        matarray = array2
	        alter = 0

	    N1 = len(matarray)
	    N2 = len(refarray)

	    matmask=[]
	    refmask=[]
	    if matarray.shape == (2,) or matarray.shape == (1,2):
	        matarray = matarray.copy().ravel()
	        diffarray = refarray - matarray
	        temp = []
	        for j,diff in enumerate(diffarray):
	            if np.abs(diff[0])<criteria and np.abs(diff[1])< criteria:
	                temp.append(j)
	                #print diff

	        if len(temp)==1:
	            matched_ref = refarray[temp[0]]
	            matched_mat = matarray
	            matmask.append(1)
	            refmask.append(temp[0])
	        else:
	            if len(temp)>1:
	                print 'more than one objects fall into the criteria region of reference star. This object is droped!'
	            stars1 = []
	            stars2 = []
	            mask1 =[]
	            mask2 =[]
	            return stars1,mask1,stars2,mask2
	    else:
	        for i,current in enumerate(matarray):
	            diffarray = refarray - current
	            if verbose:
	                print diffarray
	            temp = []
	            for j,diff in enumerate(diffarray):
	                try:
	                    diff1 = diff[0]
	                    diff2 = diff[1]
	                except:
	                    diff1 = diff[0,0]
	                    diff2 = diff[0,1]

	                if np.abs(diff1)<criteria and np.abs(diff2)< criteria:
	                    temp.append(j)
	                    #print diff
	            if len(temp)==1:
	                matmask.append(i)
	                refmask.append(temp[0])
	            elif len(temp)>1:
	                print 'more than one objects fall into the criteria region of reference star. This object is droped!'

	        if verbose:
	            print refmask
	            print matmask

	        matched_ref = refarray[refmask]
	        matched_mat = matarray[matmask]

	    #make sure the output stars1,stars2 correspond to array1,array2 respectively

	    if alter:
	        stars1 = matched_mat
	        stars2 = matched_ref
	        mask1 = matmask
	        mask2 = refmask
	    else:
	        stars1 = matched_ref
	        stars2 = matched_mat
	        mask1 = refmask
	        mask2 = matmask

	    if verbose:
	        print "From array1:\n"
	        print stars1
	        print "From array2:\n"
	        print stars2

	    return stars1,mask1,stars2,mask2

	def __find_corresponding(self,base_array,target,criteria,verbose=False):
		'''
		find sources in 'base_array' within a distance of 'criteria' to the given 'target'
		Input:
			base_array: (N,2) array
			target:     two element array or list
			criteria: scaler, the distance to the target; match radius

		Output:
			exist: True or False
			base_array[ret_mask]: (N,2) shape array and N can be 1
			ret_mask: list, indice of rows meeting requirement
		'''
		if isinstance(target,list):
			target = np.array(target)

		target = target.reshape((1,2))
		diffs = base_array-target
		distances = np.sqrt(diffs[:,0]**2 + diffs[:,1]**2)
		ret_mask= []
		for i,dist in enumerate(distances):
			if dist<criteria:
				ret_mask.append(i)
		if ret_mask == []:
			exist = False
			return exist,0,0
		else:
			exist = True
			return exist,base_array[ret_mask],ret_mask

#calibration with color term
	def explore_color_term_effect_in_calibrtion(self, flt, flt2, std_catalog='apass', refimg_dir='raw'):
		'''
		study the color term coefficients
		'''
		ps  = np.array([])
		eps = np.array([])

		for img in self.images.keys():
			if self.photometry_info[img]['flt'] == flt and self.photometry_info[img]['drop']==0:
				print img
				figoutfile = os.path.join(self.std_ref_dir, flt+'_'+img.split('.')[0]+'_cal.pdf')
				popt,perr = self.photcal_with_colorterm_single(img, flt, flt2, std_catalog=std_catalog, refimg_dir=refimg_dir, figoutfile=figoutfile)
				ps = np.append(ps,popt)
				eps = np.append(eps, perr)

		return ps, eps

	def photcal_with_colorterm_single(self, image_key, flt, flt2, std_catalog='apass', refimg_dir='raw', auto_shift_xy = False, sigclipfirst=True, sig1=6, sig2=3, figoutfile=None, target_color=None, verbose=1):
		'''
		This is created for the purpose of exploring the potential color term effect in calibration against to given standard stars

		INPUTS:
			image_key:
			flt:
			flt2: another flt to set the color
			sigclipfirst: do sigma clipping first before the calibration process
			sig1: sigma parameter for self.__sigma_clipping
			sig2: sigma parameter for self.__funcfit_remove_outlier
		'''
		imgkey_s = image_key.split('.')[0]

		refimg_abs = self.__get_internal_image(image_key, which_dir=refimg_dir)
		#refimg_abs = self.self.images[image_key]
		#self.get_standards(std_catalog)
		#stds = self.standards
		if std_catalog == 'apass':
			stds = self.get_standard_reference_star_APASS(flt,wcs=True, single_use = True, refimg_key=image_key, refimg = refimg_abs)
		elif std_catalog == 'panstarrs':
			stds = self.get_standard_reference_star_PanSTARRS(flt, flt2=flt2, photmethod = 'PSF', reference_image =refimg_abs, single_use = True, savedatakey=image_key, save_std_flt_file = True, save_std_flt_reg = False, wcsinfo_hdu='PRIMARY')
		elif std_catalog == 'refcat2':
			stds = self.get_standard_reference_star_ATLASRefcat2(flt, flt2=flt2, reference_image =refimg_abs, single_use = True, savedatakey=image_key, save_std_flt_file = True, save_std_flt_reg = False, wcsinfo_hdu='PRIMARY')
		else:
			raise ValueError('not ready')

		ref_list_file_raw = os.path.join(self.std_ref_dir, 'std_ref_whole_info_' + imgkey_s + '.txt')

		stdcolnames = self.std_catalog_colnames[std_catalog]
		magcol = stdcolnames[flt]
		magcol2 = stdcolnames[flt2]
		magerrcol  = stdcolnames[flt+'err']
		magerrcol2 = stdcolnames[flt2+'err']

		wantcols = ['x','y',magcol, magerrcol, magcol2, magerrcol2]
		if std_catalog == 'apass':
			refdata = Table.read(ref_list_file_raw, format='ascii.csv')
			refdata = refdata[(refdata[magcol2]!=0)*(refdata[magcol]!=0)*(refdata[magcol]<50)*(refdata[magcol2]<50)]
		else:
			refdata = stds
		ref_list_file = os.path.join(self.std_ref_dir, 'std_ref_'+flt+flt2+'_'+imgkey_s + '.txt')
		self.__select_columns_from_table_to_ascii(refdata, wantcols, ref_list_file)

		
		if auto_shift_xy:
			target_xy_wcs = self.__get_xy_on_image_from_wcsinfo(refimg_abs)
			stds_shift_x =  self.photometry_info[image_key]['x']- target_xy_wcs[0]
			stds_shift_y =  self.photometry_info[image_key]['y']- target_xy_wcs[1]
			if verbose:
				print "(x, y) transformation from (RA, Dec) = (%s, %s)"%(target_xy_wcs[0], target_xy_wcs[1])
				print "(x, y) target coordinate on the image = (%s, %s)"%(self.photometry_info[image_key]['x'], self.photometry_info[image_key]['y'])
				print "%s stds shift x: %s; shift y:%s"%(image_key, stds_shift_x, stds_shift_y)
		else:
			stds_shift_x = 0
			stds_shift_y = 0

		match_output_file = os.path.join(self.std_ref_dir, imgkey_s +'_std_'+flt+flt2+'.match')
		data = self.__stdcal_link_std_obs_surrounding_search_single_image(image_key, flt, ref_list_file, stds_shift_x=stds_shift_x, stds_shift_y=stds_shift_y)

		if self.reference_mag_min is not None:
			data = data[data[:,2]<self.reference_mag_min]
		if self.reference_magerr_cut is not None:
			data = data[data[:,3]<self.reference_magerr_cut]

		mags_ref    = data[:,2]
		magserr_ref  = data[:,3]
		mags_input  = data[:,8]
		magserr_input= data[:,9]
		color = data[:,2] - data[:,4]
		colorerr =  np.sqrt(data[:,3]**2+data[:,5]**2)

		mags_offset = mags_ref - mags_input
		magserr_offset = np.sqrt(magserr_ref**2+magserr_input**2)

		magoffset0 = np.median(mags_offset)
		if self.cal_offset_funcfit_type == 'o1poly':
			def func(x,a,b):
				return a*x+b
			def yerr(x,aerr,berr):
				return np.sqrt((x*aerr)**2+berr**2)
			p0 = np.array([0,magoffset0])
		elif self.cal_offset_funcfit_type == 'constant':
			def func(x,c):
				return x*0+c
			def yerr(x,cerr):
				return x*0+cerr
			p0 = np.array([magoffset0])
		else:
			raise ValueError('%s for self.cal_offset_funcfit_type not supported yet'%self.cal_offset_funcfit_type)

		if sigclipfirst:
			survive_mask,discard_mask = self.__sigma_clipping(mags_offset,variability_method='mad', sig=sig1)
			color_rmd = color[discard_mask]
			mags_offset_rmd = mags_offset[discard_mask]
			magserr_offset_rmd = magserr_offset[discard_mask]
			color_svv = color[survive_mask]
			mags_offset_svv = mags_offset[survive_mask]
			magserr_offset_svv = magserr_offset[survive_mask]

		fitdata, popt, perr = self.__funcfit_remove_outlier(color_svv, mags_offset_svv, magserr_offset_svv,func, p0, nsig=sig2, rm_mode = 'all')
		offset = func(fitdata[:,0],*popt)
	
		fitdata_outputfile = os.path.join(self.std_ref_dir, imgkey_s+'_colorterm_fitdata.txt')
		np.savetxt(fitdata_outputfile, fitdata)

		cut = int(np.floor((len(fitdata)*0.688)))
		offset_dev = np.abs(fitdata[:,1]-offset)
		indice = np.argsort(offset_dev)[cut]
		offset_err = offset_dev[indice]

		if self.cal_plot:
			from matplotlib import gridspec
			fig = plt.figure(figsize=(9, 6))
			gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
			ax0 = plt.subplot(gs[0])
			plt.xlabel('X')
			plt.ylabel('magnitude offset (reference - input)')
			ax1 = plt.subplot(gs[1],sharey=ax0)

			ax0.errorbar(color_svv, mags_offset_svv, yerr=magserr_offset_svv, fmt='gx',label='raw')
			ax0.errorbar(fitdata[:,0],fitdata[:,1],yerr=fitdata[:,2],fmt='ro',label='fitdata')
			if sigclipfirst:
				ax0.errorbar(color_rmd, mags_offset_rmd, yerr=magserr_offset_rmd, fmt='bo',label='very raw', alpha=0.5)

			ax0.plot(fitdata[:,0], offset,'k')
			ax0.plot(fitdata[:,0], offset+offset_err,'b',alpha=0.4)
			ax0.plot(fitdata[:,0], offset-offset_err,'b',alpha=0.4)

			hist,edges = np.histogram(fitdata[:,1],bins=15)
			bin_centers = edges[:-1] + (edges[1:]-edges[:-1])/2.
			ax1.step(hist,bin_centers)

			if figoutfile is not None:
				plt.title(image_key)
				plt.savefig(figoutfile)
				plt.close()
			else:
				plt.show()

		if self.cal_check_stds:
			pylab_talk_to_ds9(color, mags_offset, data, refimg_abs, xcol=0,ycol=1, yerrplot=magserr_offset)


		if target_color is not None and self.cal_offset_funcfit_type=='o1poly':
			offset = popt[1] + popt[0]*target_color
			offseterr = np.sqrt(perr[1]**2+(perr[0]*target_color)**2)
			print offset, offseterr
			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['instmag'] + offset
			self.photometry_info[image_key]['calmagerr'] = np.sqrt(self.photometry_info[image_key]['instmagerr']**2 + offseterr**2)

		return popt,perr

#Standard calibration
	def standard_calibration_flt_manual_match_method(self, flt, which_dir='raw_image'):
		'''
		standard calibration for situation where automatic matching doesn't work
		'''
		if len(self.templates) == 0:
			self.__find_template_imagekey()

		if flt not in self.templates.keys():
			raise ValueError("template image in filter %s not specified..."%flt)

		tpl_imgkey = self.templates[flt]

		if self.photometry_info[tpl_imgkey]['calmag'] == 99.99 or self.renew_stdcal_mag:
			self.__get_tpl_calibrated_mag_manual_method(tpl_imgkey, flt, which_dir=which_dir)
		else:
			print "old calibration result will be used for the template image"

		self.__dict2table()
		table_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)

		for image_key in table_flt['name']:
			if self.photometry_info[image_key]['template'] == 1:
				continue

			if self.photometry_info[image_key]['drop'] > 0:
				continue

			offset = self.photometry_info[tpl_imgkey]['calmag'] - self.photometry_info[tpl_imgkey]['relmag']
			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['relmag'] + offset

			instmagerr = self.photometry_info[image_key]['instmagerr']
			relmagerr = self.photometry_info[image_key]['relmagerr']
			stdmagerr = np.sqrt(self.photometry_info[tpl_imgkey]['calmagerr']**2 - self.photometry_info[tpl_imgkey]['instmagerr']**2)

			#Improvement needed!!!
			if self.photometry_info[tpl_imgkey]['drop'] >0 or stdmagerr==99.99:
				stdmagerr = 2*relmagerr

			calmagerr = np.sqrt(instmagerr**2 + relmagerr**2 + stdmagerr**2)
			self.photometry_info[image_key]['calmagerr'] = calmagerr

	def get_std_obs_match_filenames(self):
		'''
		prepare the filenames for calibration
		'''
		for img_key in self.images.keys():
			flt = self.photometry_info[img_key]['flt']
			stdobs_matchfile = os.path.join(self.std_ref_dir, 'std_%s_%s.txt'%(flt, img_key.split('.')[0]))
			stdobs_trans_coeffile = os.path.join(self.std_ref_dir, 'std_%s_%s_trans_coef.txt'%(flt, img_key.split('.')[0]))
			self.stdobs_match_files[img_key] = stdobs_matchfile
			self.stdobs_transcoef_files[img_key] = stdobs_trans_coeffile

	def get_secondary_stds(self, flt, secondary_source, secondary_photmethod='apphot', radec_format = False, xy_format=True, cut_radius=None):
		'''
		When primary standard stars are not enough for calibration purpose, we can use images from other resources to link the primary standard stars with stars in input image

		INPUTS:
			flt:
			secondary_source: the telescope/instrument name from which the secondard standard stars are from
			photmethod: photometry method
			radec_format: the coordiante in ra, dec format
			xy_format: the coordinate in image x,y format
			cut_radius: if not None, only return stars with distance to the supernova on the secondary-source image within 'cut_radius'
		'''
		sn_current = self.current_sn
		tel_current = self.current_telescope
		apphot_dir = self.aperture_photometry_dir
		psfphot_dir = self.__get_internal_psfphot_dir()
		photmethod_current = self.photometry_method

		second_phot_infofile = self.result_table_file_init.replace(tel_current, secondary_source)

		if photmethod_current != secondary_photmethod:
			if secondary_photmethod == 'psfphot':
				second_phot_infofile = second_phot_infofile[:-4] + '_PSF' + second_phot_infofile[-4:]
			else:
				second_phot_infofile = second_phot_infofile.replace('_PSF', '')

		if not os.path.exists(second_phot_infofile):
			print "Sorry, you request secondary standard photometry %s not exist"%second_phot_infofile
			return

		sdphotinfo = Table.read(second_phot_infofile, format='ascii.fixed_width')
		sdtpls = sdphotinfo[sdphotinfo['template'] == 1]

		sdtpl_flt = sdtpls[sdtpls['flt'] == flt]
		if len(sdtpl_flt) != 1:
			print "template image not found in secondary source, or more than one found..."
			return

		sdtplflt = sdtpl_flt[0]
		sdimg = sdtplflt['name']
		if secondary_photmethod == 'psfphot':
			sdtplphot_dir = psfphot_dir.replace(tel_current, secondary_source)
			sdtplphot_file = os.path.join(sdtplphot_dir, sdimg[:-5]+self.psfphot_ret_file_suffix)
		elif secondary_photmethod == 'apphot':
			sdtplphot_dir = apphot_dir.replace(tel_current, secondary_source)
			sdtplphot_file = os.path.join(sdtplphot_dir, sdimg[:-5]+self.apphot_ret_file_suffix)
		else:
			raise ValueError("invalid input for secondary_photomethod")

		sdx = sdtplflt['x']
		sdy = sdtplflt['y']
		sdoffset = sdtplflt['calmag'] - sdtplflt['instmag']
		sdtplphot = np.loadtxt(sdtplphot_file)
		sdtplphot[:,2] = sdtplphot[:,2] + sdoffset
		if cut_radius is not None:
			mask = np.sqrt( (sdtplphot[:,0] - sdx)**2 + (sdtplphot[:,1] -sdy)**2 ) < cut_radius
			sdtplphot = sdtplphot[mask]
		if radec_format:
			print "this feature still under construction"
		if xy_format:
			sdstdfile = "%s_%s_stds.txt"%(secondary_source, flt)
			outfile = os.path.join( self.std_ref_dir, sdstdfile )
			print "secondary stds are saved here:%s"%outfile
			np.savetxt(outfile, sdtplphot)

	def get_calmag_instmag_offset_single_image(self, image_key, photmethod='apphot', aperture_size_fixed=False, renew_existing_photret=True, img_input_which_dir='raw_image', apphot_centering=True, link_stdref_obs_method='surrounding_search', std_catalog = 'apass', panstarrs_photmethod = 'PSF', ):
		'''
		get the offset between instrumental magnitudes and calibrated magnitudes.
		'''
		if photmethod == 'apphot':
			photret_file = image_key.split('.')[0]  + self.apphot_ret_file_suffix
			photret_file_abs = os.path.join(self.aperture_photometry_dir,photret_file)
			if os.path.exists(photret_file_abs) and (not renew_existing_photret):
				print "existing photometry file will be used"
			else:
				print image_key
				xys = self.__get_xys_on_image(image_key)
				image = self.__get_internal_image(image_key, which_dir=which_dir)
				self.get_apphot_iraf_parameters(image_key, saveparfile=False)
				options = self.apphot_iraf_options.copy()
				print options
				photret = self.aperture_photometry_apphot_iraf_single_image_xys_given(image, xys, options, output=None, centering=apphot_centering)
				#xpos, ypos,mag,mag_err
				self.__delete_file_if_exist(photret_file_abs)
				np.savetxt(photret_file_abs,photret,fmt="%6.2f %6.2f %6.3f %6.3f")
			flt = self.photometry_info[image_key]['flt']
			offset, offset_err = self.stdcal_single_image(image_key,flt, instmag=0, instmagerr=0, link_stdref_obs_method= link_stdref_obs_method, std_catalog = std_catalog, panstarrs_photmethod =  panstarrs_photmethod,  std_world2img_refimage=None,  stds_shift_x=0, stds_shift_y=0, std_obs_matched = False, vstdcal_single_imageerbose = 0, update_stdmag=False)
		elif photmethod == 'psfphot':
			print "not available yet..."
			offset = 99.99
			offset_err = 99.99
		else:
			raise ValueError("not recognised input photometry method: %s"%photmethod)
		return offset, offset_err

	def standard_calibration_direct(self, flts='all', std_catalog='apass', auto_shift_xy=True, verbose=0, refimg_which_dir='raw_image', update_photret_table=0):
		'''
		directly calibrate the science image to the standard catalogue without the relative calibration within the same band

		See self.stdcal_single_image for detail
		'''
		if flts == 'all':
			if len(self.templates)==0:
				self.__find_template_imagekey(updatetable=1)
			fltstocal = self.templates.keys()
		elif isinstance(flts,list) or isinstance(flts,np.ndarray):
			fltstocal = flts
		else:
			raise IOError('Invalid input for flts')

		for flt in fltstocal:
			if update_photret_table:
				self.__dict2table()
			fltdata = self.photometry_record_table[self.photometry_record_table['flt']==flt]
			for img in fltdata['name']:
				if self.photometry_info[img]['drop'] == 0 and (self.photometry_info[img]['calmag'] == 99.99 or self.renew_stdcal_mag):
					self.stdcal_single_image(img, flt, std_catalog=std_catalog, link_stdref_obs_method='surrounding_search', auto_shift_xy=auto_shift_xy, refimg_which_dir = refimg_which_dir, verbose=verbose)

	def standard_calibration(self,flts='all', std_catalog = None,  tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False,link_stdref_obs_method=None, std_obs_matched = False, stds_shift_x = 0, stds_shift_y = 0, updatetable=1):
		'''


		Inputs:
			flts: 'all', all filters whose template image found
			      list or np.ndarray, for example flts=['B','V','rp']

			std_catalog: 'apass', '2mass' or 'panstarrs'

			tpl_has_wcs: template images already have wcs information?
					If True, the wcs information will be applied to get the image positions
					of standard stars and use grmatch to get more precise match between
					stars on template image and stars from standard database;
				    	If False, astrometry will be needed to put standard stars on the template image;
					whether all templates need astrometry solution
					depend on 'all_tpl_astrometry'

			all_tpl_astrometry: get astrometry for all templates?
					If True, get astrometry for each template;
					If False, astrometry will be obtained for the best template, and other
					templates use that solution to get a rough image positions and
					 apply grmatch further.

			renew_tpl_astrometry:


			link_stdref_obs_method: the following options available
						'surrounding_search': this require the template image has correct wcs solution; and stars at the same position within given criteria are treated as iddentical stars
						'grmatch': match stdref stars and obs stars with 'grmatch' task; template image has valid wcs solution but not correct/precise
						'world_image_grmatch': the template image doesn't have wcs solution and difficult to obtain the wcs; then perform some simple world coordinate to image coordinate transformation based on knowledge of image size and CCD pixel scale and the [RA, Dec] of one pixel on the image, then perform grmatch the same as option 'grmatch'

			std_obs_matched:


		Please see self.standard_calibration_flt for details
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		flts_obs = self.templates.keys()

		if flts == 'all':
			for flt in flts_obs:
				print "working on %s band"%flt
				self.standard_calibration_flt(flt, std_catalog = std_catalog, tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry,std_obs_matched = std_obs_matched,link_stdref_obs_method=link_stdref_obs_method, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y)
		elif isinstance(flts,list) or isinstance(flts,np.ndarray):
			for flt in flts:
				print "working on %s band"%flt
				self.standard_calibration_flt(flt, std_catalog = std_catalog, tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry, std_obs_matched = std_obs_matched,link_stdref_obs_method = link_stdref_obs_method, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y)
		else:
			raise IOError('Invalid input for flts')

	def standard_calibration_flt_simple(self, flt, m_zpt, em_zpt):
		'''
		calmag = instmag + m_zpt
		'''
		self.__dict2table()
		table_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		for image_key in table_flt['name']:
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if (self.photometry_info[image_key]['calmag'] != 99.99) and (not self.renew_stdcal_mag):
				continue

			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['instmag'] + m_zpt
			calmagerr = np.sqrt(em_zpt**2 + self.photometry_info[image_key]['instmagerr']**2)
			self.photometry_info[image_key]['calmagerr'] = calmagerr

	def standard_calibration_flt(self,flt, std_catalog = None,  tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False,link_stdref_obs_method = None,std_obs_matched = False, stds_shift_x =0, stds_shift_y =0, updatetable=1):
		'''
		standard calibration for image with filter 'flt'
		The current philosophy here is to perform calibration for the reference image and apply the calibration result to all images

		check the self.photometry_info[template_imgkey]['calmag'] first. If the value is not 99.99, then skip the calibration for reference image

		Inputs:
			flt:

			tpl_has_wcs: same as self.standard_calibration
			all_tpl_astrometry: same as self.standard_calibration
			renew_tpl_astrometry: same as self.standard_calibration


		Involved functions:
			self.__stdcal_template_prepare_tplimg_flt
			self.__stdcal_template_prepare_stdfile_flt
			self.__stdcal_template_flt

		'''
		apass_flts = ['B','V','R','I','gp','rp','ip']
		twomass_flts = ['J', 'H', 'K']
		self.__find_template_imagekey(updatetable=updatetable)

		flt_tplkey = self.templates[flt]
		if self.photometry_info[flt_tplkey]['calmag'] == 99.99:
			if not std_obs_matched:
				if link_stdref_obs_method != 'world_image_grmatch':
					self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
				if std_catalog is not None:
					which_stdcatalog = std_catalog
				else:
					if flt in apass_flts:
						which_stdcatalog = 'apass'
					elif flt in twomass_flts:
						which_stdcatalog = '2mass'
					else:
						raise ValueError('%s band photometry calibration is not supported'%flt)

				self.__stdcal_template_prepare_stdfile_flt(flt,link_stdref_obs_method = link_stdref_obs_method, which_catalog=which_stdcatalog)
				if link_stdref_obs_method is None:
					if not tpl_has_wcs and all_tpl_astrometry:
						link_stdref_obs_method = 'surrounding_search'
					else:
						link_stdref_obs_method = 'grmatch'

			self.__stdcal_template_flt(flt, link_stdref_obs_method = link_stdref_obs_method, std_obs_matched = std_obs_matched, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y )
		elif self.renew_stdcal_mag:
			if not std_obs_matched:
				if link_stdref_obs_method != 'world_image_grmatch':
					self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
				if std_catalog is not None:
					which_stdcatalog = std_catalog
				else:
					if flt in apass_flts:
						which_stdcatalog = 'apass'
					elif flt in twomass_flts:
						which_stdcatalog = '2mass'
					else:
						raise ValueError('%s band photometry calibration is not supported'%flt)

				self.__stdcal_template_prepare_stdfile_flt(flt,link_stdref_obs_method = link_stdref_obs_method, which_catalog=which_stdcatalog)
				if link_stdref_obs_method is None:
					if not tpl_has_wcs and all_tpl_astrometry:
						link_stdref_obs_method = 'surrounding_search'
					else:
						link_stdref_obs_method = 'grmatch'

			self.__stdcal_template_flt(flt, link_stdref_obs_method = link_stdref_obs_method, std_obs_matched = std_obs_matched , stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y )

		else:
			print "old standard calibration result for images in filter %s will be used"%flt

		if updatetable:
			self.__dict2table()
		table_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)

		for image_key in table_flt['name']:
			tpl_imgkey = self.templates[flt]
			if self.photometry_info[image_key]['template'] == 1:
				continue
			if self.photometry_info[image_key]['drop'] > 0:
				continue
			if (self.photometry_info[image_key]['calmag'] != 99.99) and (not self.renew_stdcal_mag):
				continue

			offset = self.photometry_info[tpl_imgkey]['calmag'] - self.photometry_info[tpl_imgkey]['relmag']
			self.photometry_info[image_key]['calmag'] = self.photometry_info[image_key]['relmag'] + offset

			instmagerr = self.photometry_info[image_key]['instmagerr']
			relmagerr = self.photometry_info[image_key]['relmagerr']
			stdmagerr = np.sqrt(self.photometry_info[tpl_imgkey]['calmagerr']**2 - self.photometry_info[tpl_imgkey]['instmagerr']**2)

			#Improvement needed!!!
			if self.photometry_info[tpl_imgkey]['drop'] >0 or stdmagerr==99.99:
				stdmagerr = 2*relmagerr

			calmagerr = np.sqrt(instmagerr**2 + relmagerr**2 + stdmagerr**2)
			self.photometry_info[image_key]['calmagerr'] = calmagerr


	def get_astrometry_solution_ccmap(self, imgkey, stars_list_source, ):
		'''
		first get reference sources and input sources; match sources from two catalog; ccmap solution
		
		INPUTS:
			stars_list_source: 'apphot', 'psfphot', 'prephot'
		'''
		imgkey_s = imgkey.split('.')[0]
		if stars_list_source == 'apphot':
			input_list_file = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
			#column names in ref_list and input_list: xpos, ypos, mag,mag_err
		elif stars_list_source == 'psfphot':
			input_list_file = self.__get_internal_psfphot_mag_file(imgkey)
			# XCENTER YCENTER MAG MERR
		elif stars_list_source == 'prephot':
			input_list_file = os.path.join(self.stars_dir, imgkey_s+self.starfile_suffix)	
			# XCENTER YCENTER MAG SHARPNESS SROUND GROUND ID
		else:
			raise IOError("Invalid input for stars_list_source")

		inxymags_regionfile = os.path.join(os.path.dirname(input_list_file), imgkey_s + '.reg')

		input_stars = Table.read(input_list_file, format='ascii.fast_no_header')
		return input_stars


	def __get_astrometry_solution_grtrans(self, ref_catalog_input, ref_catalog_output, input_source_list, racol, deccol, magcol, raref=None, decref=None, verbose=0):
		'''
		get astrometry with programs grmatch and grtrans from fitsh
		INPUTS:
			ref_catalog_input: the input reference star catalog with ra,dec columns
			ref_catalog_output: the output reference star catalog with image x,y columns appended
			input_source_list: formated as x,y,mag,magerr
			racol: RA column index in input_catalog
			deccol: Dec column index in input_catalog
			magcol: magnitude column index in input_catalog
			raref: the reference RA for sky projection to image plane; around the center of FOV
			decref: the reference Dec for sky projection to image plane; around the center of FOV

		!!!UNDER DEVELOPMENT
		'''
		if raref is None:
			raref = self.sn_ra_world_deg
		if decref is None:
			decref = self.sn_dec_world_deg
		if self.pixscale is None:
			self.pixscale = 1.0
		scale = 1.0/(self.pixscale/3600.0/180*np.pi)
		grtrans = os.path.join(self.base.fitsh_dir, 'grtrans')
		temp = np.loadtxt(ref_catalog_input)
		M,N = temp.shape
		xrefcol = N+1
		yrefcol = N+2
		grtranscmd1 = "%s --input %s --wcs tan,scale=%s,ra=%s,dec=%s --col-radec %s,%s --col-out %s,%s --output %s"%(grtrans, ref_catalog,scale,raref,decref,racol, deccol,xrefcol,yrefcol,ref_catalog_output)
		try:
			if verbose:
				print grtranscmd1
			os.system(grtranscmd1)
		except:
			raise OSError("grtranscmd1 failure, pls check")

		self.fitsh_grmatch_pointmatch_pars['--col-ref'] = '%s,%s'%(xrefcol,yrefcol) #The index of the first column is always 1
		#self.fitsh_grmatch_pointmatch_pars['--col-inp'] = '1,2'
		#self.fitsh_grmatch_pointmatch_pars['--order'] = 1 #If the order is A, >= (A+1)*(A+2)/2 valid points are needed to fit the transformation
		#self.fitsh_grmatch_pointmatch_pars['--max-distance'] = 1
		#self.fitsh_grmatch_pointmatch_pars['--triangulation'] = 'auto,mixed,maxnumber=200'
		self.fitsh_grmatch_pointmatch_pars['--col-ref-ordering'] = -magcol #negative sign indicates ascending, small values first
		#self.fitsh_grmatch_pointmatch_pars['--col-inp-ordering'] = -3

		self.fitsh_grmatch_type = 'point'
		self.__fitsh_grmatch(ref_catalog_output, input_source_list, match_output, trans_output=trans_output, mode ='file')


	def __update_wcs(self, imgkey, wcsfile, input_img_dir = 'raw'):
		'''
		update the wcs of the input image and output to self.astrometry_dir
		'''
		inimage = self.__get_internal_image(imgkey, which_dir=input_img_dir)
		outimage = self.__get_internal_image(imgkey,which_dir='wcs')
		wcs = fits.open(wcsfile)[0].header
		hdu = fits.open(inimage)
		hdu[0].header.update(wcs)
		hdu.writeto(outimage)
		hdu.close()
	
	def __compute_distance_on_sphere(self, ras, decs, refra, refdec):
		'''
		all inputs are in degree
		'''
		distances  = np.arccos(np.sin(decs/180.0*np.pi)*np.sin(refdec/180.0*np.pi) + np.cos(decs/180.0*np.pi)*np.cos(refdec/180.0*np.pi)*np.cos((ras-refra)/180.0*np.pi))/np.pi*180
		return distances

	def __get_sources_within_given_region(self, sources, racol, deccol, refra, refdec, radius):
		'''
		get sources within a specified region in spherical coordinate
		'''
		distances = self.__compute_distance_on_sphere(sources[racol].data, sources[deccol].data, refra, refdec)
		sources_found = sources[distances<radius]
		return sources_found

	def __get_astrometry_solution_API_client(self,img,wcs_out, racenter=None, deccenter=None, pixelscale=None, radius =None, time_out_max = 300, apikey ='aubtqramhfmriufp',verbose=1):
		'''
		Get astrometry for image 'img' with nova.astrometry.net API and the wcs result will be saved in 'wcs_out'

		Please see 'astrometry_API/client.py' for details
		'''

		api_dir = os.path.join(self.base.base_dir, 'extern/nova_astrometry_api')
		API_client = os.path.join(api_dir, "client.py")
		if not os.path.exists(API_client):
			raise IOError("API client for astrometry not available")

		img_key = os.path.basename(img)
		print "working on astrometry for %s"%img

		if racenter is None:
			racenter = self.sn_ra_world_deg
		if deccenter is None:
			deccenter = self.sn_dec_world_deg
		
		if radius is None:
			imgfr = 0.07 
		if pixelscale is None:
			pixelscale = 0.2
		if verbose:
			monitor = '-w'
		else:
			monitor = ' >/dev/null'

		command_line = "python %s --server http://nova.astrometry.net/api/ --apikey %s --wcs %s --upload %s --ra %s --dec %s --radius %s --scale-units arcsecperpix --scale-est %s -p %s"%(API_client, apikey,  wcs_out, img, racenter, deccenter, imgfr, pixelscale,monitor)
		if verbose:
			print command_line
		command = Command(command_line)
		failure_signal = command.run(timeout=time_out_max)
		#failure_signal  non-zero value means failure

		if not failure_signal:
			success = True
		else:
			success = False

		return success

	def __stdcal_template_prepare_tplimg(self,flts='all',tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False, updatetable=1):
		'''
		Prepare template images
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		flts_obs = self.templates.keys()

		if flts == 'all':
			for flt in flts_obs:
				self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_wcs,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
		elif isinstance(flts,list) or isinstance(flts,np.ndarray):
			for flt in flts:
				self.__stdcal_template_prepare_tplimg_flt(flt,tpl_has_wcs=tpl_has_astrometry,all_tpl_astrometry=all_tpl_astrometry,renew_tpl_astrometry=renew_tpl_astrometry)
		else:
			raise IOError('Invalid input for flts')

	def __stdcal_template_prepare_tplimg_flt(self,flt,tpl_has_wcs=False,all_tpl_astrometry=True,renew_tpl_astrometry=False):
		'''
		prepare reference image which will be used to link the calibration reference stars
		from both standard datebase and new measurements.
		Inputs:
			refer to 'standard_calibration' for explanation
		'''
		maxtime_on_astrometry = self.maxtime_on_astrometry
		tpl_key = self.templates[flt]

        	if self.current_telescope == 'WFCAM':
			tpl_img_orig = os.path.join(self.modified_image_dir, tpl_key)
		else:
	        	tpl_img_orig = os.path.join(self.raw_image_dir,tpl_key)
	        tpl_img_for_stdcal = self.templates_after_astrometry[flt]

		finish_code = 1
                if tpl_has_wcs:
                	if not os.path.exists(tpl_img_for_stdcal) or renew_tpl_astrometry:
				shutil.copy(tpl_img_orig, tpl_img_for_stdcal)

                elif all_tpl_astrometry:
                	if not os.path.exists(tpl_img_for_stdcal) or renew_tpl_astrometry:
	        		success = self.__get_astrometry_solution_API_client(tpl_img_orig,tpl_img_for_stdcal, time_out_max = maxtime_on_astrometry)
				if not success:
					finish_code = 0
		else:
			if self.reference_flt is not None:
				ref_flt = self.reference_flt
				tpl_ref_flt_key  = self.templates[ref_flt]
				tpl_img_ref_flt_orig = os.path.join(self.raw_image_dir,tpl_ref_flt_key)
				tpl_img_for_stdcal_ref_flt = self.templates_after_astrometry[ref_flt]

				if not os.path.exists(tpl_img_for_stdcal_ref_flt) or renew_tpl_astrometry:
					success = self.__get_astrometry_solution_API_client(tpl_img_ref_flt_orig,tpl_img_for_stdcal_ref_flt, time_out_max = maxtime_on_astrometry)
					if not success:
						finish_code = 0
				if flt != ref_flt:
					shutil.copy(tpl_img_for_stdcal_ref_flt,tpl_img_for_stdcal)
			else:
				raise KeyError("Specify self.reference_flt please!")

		return finish_code

	def __stdcal_template_prepare_stdfile_flt(self,flt,link_stdref_obs_method=None, which_catalog = 'apass'):
		'''
		Get the reference stars for magnitude calibration

		INPUTs:
			flt:
			link_stdref_obs_method: see self.standard_calibration
			which_catalog: 'apass' or '2mass'

		'''
		if link_stdref_obs_method == "world_image_grmatch":
			wcs = False
		else:
			wcs = True

		if which_catalog == 'apass':
			self.get_standard_reference_star_APASS(flt, wcs=wcs)
		elif which_catalog == '2mass':
			self.get_standard_reference_star_2MASS(flt, wcs=wcs)
		elif which_catalog == 'panstarrs':
			magtype = self.panstarrs_mag_photometry_method
			self.get_standard_reference_star_PanSTARRS(flt, photmethod = magtype, single_use = False)
		elif which_catalog == 'refcat2':
			self.get_standard_reference_star_ATLASRefcat2(flt, single_use = False)
		else:
			raise ValueError('catalog %s is not supported...'%which_catalog)

	def __stdcal_template_link_std_obs_flt(self,flt, method = 'grmatch', stds_shift_x=0, stds_shift_y = 0):
		'''
		get the one-one mapping between standards and new measurements

		Inputs:
			methods: available options currently available:
				grmatch:
				selection_on_ds9:
				external:
				surrounding_search:
		'''
		if method == 'grmatch':
			self.__stdcal_template_link_std_obs_grmatch_flt(flt)
		elif method == 'world_image_grmatch':
			self.__stdcal_template_link_std_obs_grmatch_flt(flt)
		elif method == 'surrounding_search':
			self.__stdcal_template_link_std_obs_surrounding_search_flt(flt, stds_shift_x=stds_shift_x, stds_shift_y= stds_shift_y)
		elif method == 'selection_on_ds9':
			self.__stdcal_template_link_std_obs_selection_on_ds9_flt(flt)
		elif method == 'external':
			self.__stdcal_template_link_std_obs_selection_external_flt(flt)
		else:
			raise IOError("Invalid input for method, available options are grmatch, selection_on_ds9 and external")

	def __stdcal_link_std_obs_single_image(self,imgkey,flt,method = 'grmatch', single_use = False, stds_shift_x=0, stds_shift_y =0, external_refimg=None, image_dir='raw_image'):
		'''
		match sources on given image with standards

		INPUTS:
			image:
			flt:
			method: 'grmatch', 'surrounding_search', 'xytran', 'selection_on_ds9', 'external'
			single_use:
			stds_shift_x: shift standards in x direction by this value if method is 'surrounding_search'
			stds_shift_y: shift standards in y direction by this value if method is 'surrounding_search'
		'''
		imgkey_s = imgkey.split('.')[0]
		if not single_use:
			ref_list_file = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
		else:
			ref_list_file = os.path.join(self.std_ref_dir, 'std_ref_' + imgkey_s + '.txt')

		if not os.path.exists(ref_list_file):
			print "The required standards file %s not avaiable"%ref_list_file

		if method == 'grmatch':
			match_ret = self.__stdcal_template_link_std_obs_grmatch_single_image(imgkey,flt, single_use = single_use)
		elif method == 'surrounding_search':
			match_ret = self.__stdcal_link_std_obs_surrounding_search_single_image(imgkey,flt, ref_list_file, stds_shift_x=stds_shift_x, stds_shift_y= stds_shift_y)
		elif method == 'xytran':
			match_ret = self.__stdcal_template_link_std_obs_xytrans_single_image(imgkey, flt, single_use = single_use, verbose = 1, external_refimg=external_refimg, image_dir=image_dir)
		elif method == 'selection_on_ds9':
			match_ret = self.__stdcal_template_link_std_obs_selection_on_ds9_single_image(imgkey,flt)
		elif method == 'external':
			match_ret = self.__stdcal_template_link_std_obs_selection_external_single_image(imgkey,flt)
		else:
			raise IOError("Invalid input for method, available options are grmatch, surrouding_search, xytran, selection_on_ds9 and external")

		return match_ret

	def __stdcal_template_link_std_obs_surrounding_search_flt(self,flt, stds_shift_x=0, stds_shift_y = 0):
		'''
		Get one-one mapping between standard reference stars and sources detected in template image
		'''
		tpl_imgkey = self.templates[flt]
      		ref_list_file = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
		match_ret = self.__stdcal_link_std_obs_surrounding_search_single_image(tpl_imgkey,flt, ref_list_file)

		return match_ret

	def __stdcal_template_link_std_obs_xytrans_single_image(self, imgkey,flt, single_use = False, verbose = False, external_refimg=None, image_dir='raw_image'):
		'''
		First manually select 3 pair of sources from std catalog and new measurements; then use the selected pairs to do transformation for the whole catalog; finally do 		  the normal surrounding search
		'''
		imgkey_s = imgkey.split('.')[0]

		if not single_use:
			ref_list_file = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
			stdxymags_regionfile =  os.path.join(self.std_ref_dir, flt + '.reg')
		else:
			ref_list_file = os.path.join(self.std_ref_dir, 'std_ref_' + imgkey_s + '.txt')
			stdxymags_regionfile = os.path.join(self.std_ref_dir, imgkey_s + '.reg')

		if self.photometry_method == 'apphot':
			input_list_file = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
			#column names in ref_list and input_list: xpos, ypos, mag,mag_err
			inxymags_regionfile = os.path.join(self.aperture_photometry_dir, imgkey_s + '.reg')
		elif self.photometry_method == 'psfphot':
			input_list_file = self.__get_internal_psfphot_mag_file(imgkey)
			inxymags_regionfile = os.path.join(os.path.dirname(input_list_file), imgkey_s + '.reg')
		else:
			raise IOError("Invalid input for photometry_method")

		inxymags = np.loadtxt(input_list_file)
		stdxymags = np.loadtxt(ref_list_file)
		refimg = self.__get_internal_image(imgkey, which_dir=image_dir)

		match_listfile = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref_xytrans.input')
		match_database = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref_xytrans.coef')

		if os.path.exists(match_listfile) and os.path.exists(match_database) and (not self.renew_stdcal_xytran):
			print "%s and %s exist; if you want to renew, modify self.renew_stdcal_xytran"%(match_listfile,match_database)
		else:
			if not os.path.exists(stdxymags_regionfile):
				create_ds9_region_file(stdxymags, stdxymags_regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='green', width=1, load_text=False,  text_uoffset=25, text_loffset=25, textcol1=2,textcol2=3)
			if external_refimg is not None:
				dispimg = external_refimg
			else:
				dispimg = refimg

			select_num = self.stdcal_xytran_stdnum

			stdlist = self.input_xys_through_ds9_get_mags(dispimg, stdxymags, regionfile=stdxymags_regionfile, wantnum=select_num, verbose=1)#output table colnames: num, x, y, mag, magerr
			if not os.path.exists(inxymags_regionfile):
				create_ds9_region_file(inxymags, inxymags_regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='green', width=1, load_text=False,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)
			inlist = self.input_xys_through_ds9_get_mags(refimg, inxymags, regionfile=inxymags_regionfile, wantnum=select_num,  verbose=1)

			stdlist.write(os.path.join(self.std_ref_dir, 'std_ref_pickups_for_%s.txt'%imgkey_s), format='ascii.fixed_width')
			inlist.write(os.path.join(self.std_ref_dir, 'obs_pickups_on_%s.txt'%imgkey_s), format='ascii.fixed_width')

			reference_list = Table()
			reference_list.add_columns([stdlist['x'], stdlist['y']])
			reference_list.rename_column('x','xstd')
			reference_list.rename_column('y','ystd')
			reference_list.add_columns([inlist['x'], inlist['y']])
			reference_list.write(match_listfile, format='ascii.no_header')


		geomap_iraf(match_listfile, match_database)

		tg_xys = inxymags[:,0:2]
		ref_xys_beforetrans = stdxymags[:,0:2]
		ref_xys_beforetrans_file = os.path.join(self.std_ref_dir, 'std_ref_'+imgkey_s + '_xy_beforetran.txt')
		np.savetxt(ref_xys_beforetrans_file, ref_xys_beforetrans)
		ref_xys_aftertran_file   = os.path.join(self.std_ref_dir, 'std_ref_'+imgkey_s + '_xy_aftertran.txt')
		geoxytran_iraf(ref_xys_beforetrans_file, ref_xys_aftertran_file, match_database, match_listfile)

		ref_xys = np.loadtxt(ref_xys_aftertran_file)
		criteria = self.criteria_stdref_tpl_match
		target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)

		if verbose:
			print indexmat

		match_ret = np.hstack((stdxymags[indexref], inxymags[indexmat]))

        	std_ref_and_obs_measurement_stars_file = imgkey_s +'_std_'+flt+'.match'
        	match_output_file  = os.path.join(self.std_ref_dir, std_ref_and_obs_measurement_stars_file)
		np.savetxt(match_output_file,match_ret)

		match_result_plot = self.stdref_tpl_match_result_display

		if match_result_plot:
			tg_xys_matched  = tg_xys[indexmat]
			ref_xys_matched = ref_xys[indexref]
			plt.plot(tg_xys_matched[:,0]-ref_xys_matched[:,0],tg_xys_matched[:,1]-ref_xys_matched[:,1],'o')
			plt.show()

		return match_ret

	def __stdcal_link_std_obs_surrounding_search_single_image(self, imgkey, flt, ref_list_file, stds_shift_x=0, stds_shift_y=0, match_output_file=None, verbose = False):
		'''
		match two catalog according to the (x,y) coordinates.
		The two catalogs (std and obs) are obtained for image with imgkey.
		The first two columns of both catalogs are image coordinate x and y.

		INPUTS:
			imgkey:
			flt:
			single_use: if True, 'the ref_list_file' will be slightly different
			stdwhole:if True, the std file contain whole package of catalog information
			stds_shift_x:
			stds_shift_y:
		'''
		imgkey_s = imgkey.split('.')[0]
		if match_output_file is None:
			match_output_file  = os.path.join(self.std_ref_dir, imgkey_s +'_std_'+flt+'.match')

		if self.photometry_method == 'apphot':
			input_list_file = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
			#column names in ref_list and input_list: xpos, ypos, mag,mag_err
		elif self.photometry_method == 'psfphot':
			input_list_file = self.__get_internal_psfphot_mag_file(imgkey)
		else:
			raise IOError("Invalid input for photometry_method")
		if verbose:
			print input_list_file

		'''	
		input_list = np.loadtxt(input_list_file)
		tg_xys = input_list[:,0:2]
		ref_list = np.loadtxt(ref_list_file)
		ref_xys = ref_list[:,0:2]
		ref_xys[:,0] = ref_xys[:,0] + stds_shift_x
		ref_xys[:,1] = ref_xys[:,1] + stds_shift_y

		criteria = self.criteria_stdref_tpl_match
		target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)

		if verbose:
			print indexmat

		match_ret = np.hstack((ref_list[indexref],input_list[indexmat]))
		np.savetxt(match_output_file,match_ret)

		match_result_plot = self.stdref_tpl_match_result_display
		if match_result_plot:
			tg_xys_matched  = tg_xys[indexmat]
			ref_xys_matched = ref_xys[indexref]
			fig = plt.figure(figsize=(6,6))
			plt.plot(tg_xys_matched[:,0]-ref_xys_matched[:,0],tg_xys_matched[:,1]-ref_xys_matched[:,1],'o')
			plt.show()
		'''

		input_list = np.loadtxt(input_list_file)
		ref_list = np.loadtxt(ref_list_file)
		ref_list[:,0] = ref_list[:,0] + stds_shift_x
		ref_list[:,1] = ref_list[:,1] + stds_shift_y
		self.fitsh_grmatch_type = 'coord'
		self.fitsh_grmatch_coordmatch_pars['--max-distance'] = self.criteria_stdref_tpl_match
		#match_ret = self.__fitsh_grmatch(ref_list_file, input_list_file, match_output_file, mode ='file') #
			
		#xyid_matchedfile = os.path.join(self.std_ref_dir, imgkey_s +'_std_xyidxyid.match') 
		#matched_xyid = self.__fitsh_grmatch(ref_xys, tg_xys, xyid_matchedfile, mode ='data') #
		#matchedref = ref_list[map(int, matched_xyid[:,2])] #need to append index to original data
		#matchedinput = input_list[map(int, matched_xyid[:,5])]
		#match_ret = np.hstack((matchedref, matchedinput))
		#np.savetxt(match_output_file,match_ret)

		match_ret = self.__fitsh_grmatch(ref_list, input_list, match_output_file, mode ='data') #
		NX,NY = ref_list.shape

		match_result_plot = self.stdref_tpl_match_result_display
		if match_result_plot:
			fig = plt.figure(figsize=(6,6))
			plt.plot(match_ret[:,0]-match_ret[:,NY], match_ret[:,1]-match_ret[:,NY+1],'o')
			plt.show()

		return match_ret

	def __stdcal_template_link_std_obs_selection_on_ds9_flt(self,flt):
		tpl_imgkey = self.templates[flt]

		self.__stdcal_template_link_std_obs_selection_on_ds9_single_image(tpl_imgkey,flt)

		raise IOError("under construction...")

	def get_template_reference_mags_flt(self, flt,  which_dir='raw_image', renew_refmags=False):
		'''
		get the template mags table. The output format is: id, x, y, mag, magerr

		INPUTS:
			flt:
			which_dir:
			renew_refmags:
		'''
		if len(self.templates) == 0:
			self.__find_template_imagekey()
		if flt not in self.templates.keys():
			raise ValueError("no template image avaiable for filter %s"%flt)
		image_key = self.templates[flt]
		tpl_imgkey_s = image_key.split('.')[0]

		if self.photometry_method is None:
			self.photometry_method = raw_input("please input the photometry method used(apphot/psfphot):")

		if self.photometry_method   == 'apphot':
			ref_list_filename      = os.path.join(self.aperture_photometry_dir, tpl_imgkey_s + self.apphot_ret_file_suffix)
		elif self.photometry_method == 'psfphot':
			ref_list_filename      = self.__get_internal_psfphot_mag_file(image_key)
		else:
			raise IOError("Invalid input for photometry_method")

		if not os.path.exists(ref_list_filename):
			raise IOError("PSF photometry result file %s not exists"%ref_list_filename)

		psfphot_ret = np.loadtxt(ref_list_filename)
		xys = psfphot_ret[:,0:2]
		mags = psfphot_ret[:,2]
		magerrs = psfphot_ret[:,3]

		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		regionfile = os.path.join(self.template_dir, "%s_%s_xymag.reg"%(flt, tpl_imgkey_s))
		if (not os.path.exists(regionfile)) or renew_refmags:
			create_ds9_region_file(psfphot_ret, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='red', width=1, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

		refmag_table = self.input_xys_through_ds9_get_mags(input_image, psfphot_ret, regionfile=regionfile)
		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		refmag_table.write(refmag_filename, format='ascii.fixed_width')

		return refmag_table

	def __input_xys_get_mags(self, xysmags, idxys, verbose=1, search_radius=None):
		'''
		Extract mags from xysmags (x,y,mag,magerr) with input image coordinate (x,y) from idxys

		INPUTS:
			xysmags:
			idxys:
			search_radius: if None, then self.criteria_point_match_general will be used
		'''
		xys = xysmags[:,0:2]
		if search_radius is None:	
			criteria = self.criteria_point_match_general 
		else:
			criteria = search_radius 
		num_flags = []
		xs = []
		ys = []
		mags = []
		magerrs = []

		for idnum, x, y in idxys:
			xy_target = [x,y]
			yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, criteria)
			if (not yesfound):
				print "no object found within criteeia..."
				mag = 99.99
				magerr = 99.99
			elif len(index)>1:
				print 'more than one objects fall into the criteria region of reference star. This object is droped!'
				mag = 99.99
				magerr = 99.99
			else:
				print "the corresponding object found at (%s)"%xy_match
				mag = xysmags[index,2][0]
				magerr = xysmags[index,3][0]

			num_flags.append(idnum)
			xs.append(x)
			ys.append(y)
			mags.append(mag)
			magerrs.append(magerr)

		num_flags_col = Column(name='num',data=num_flags)
		xs_col = Column(name='x', data = xs)
		ys_col = Column(name='y', data = ys)
		mags_col = Column(name = 'mag', data = mags)
		magerrs_col = Column(name = 'magerr',data = magerrs)
		out_table = Table()
		out_table.add_columns([num_flags_col, xs_col, ys_col, mags_col, magerrs_col])
		if verbose>1:
			print out_table
		return out_table

	def input_xys_through_ds9_get_mags(self, input_image, xysmags, stopandask=False, regionfile=None, wantnum=None, newds9=True, verbose=1):
		'''
		pick up stars from ds9 display and extract mags from file xysmagfile

		INPUTS:
			input_image: the reference image to load the xysmags
			xysmags: format x,y,mag,magerr
			stopandask: whether to stop and ask for the id for each pick up in the loop
			regionfile: the file contains regions to be displayed while picking on image
			wantnum: the number of points to be picked. It can be a integer number which give the number of points to be selected or None which allow the pick process continue until no valid point obtained
			newds9: whether refresh the ds9 image
		'''
		if verbose:
			print "Now start clicking on ds9 to pick up sources"
			print "You can modify self.criteria_point_match_general  to define the searching cone radius"
		xys = xysmags[:,0:2]

		ids = []
		xs = []
		ys = []
		mags = []
		magerrs = []

		id = 0
		while True:
			if wantnum is not None:
				if id>(wantnum-1):
					break
			if id==0:
				if newds9:
					newds9_loop = True
				else:
					newds9_loop = False
				regionfile = regionfile
			else:
				newds9_loop = False
				regionfile = None

			if stopandask:
				id_input = raw_input('id for the next source (default: %s)'%str(id)) or id
				id = int(id_input)

			print "Please pick up the star of id =%s"%str(id)
			xy = self.get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=newds9_loop)
			if xy is None:
				break
			print "xy from mouse pick up:", xy
			x = xy[0]
			y = xy[1]
			xy_target = [x,y]

			criteria  = self.criteria_point_match_general
			yesfound, xy_match,index = self.__find_corresponding(xys, xy_target, criteria)
			if (not yesfound):
				print "No object found within criteria..."
				print "You can change self.criteria_point_match_general to give a new search distance"
				mag = 99.99
				magerr = 99.99
			elif len(index)>1:
				print 'more than one objects fall into the criteria region of reference star. This object is droped!'
				print "You can change self.criteria_point_match_general to give a new search distance"
				mag = 99.99
				magerr = 99.99
			else:
				print "the corresponding object found at (%s)"%xy_match
				mag = xysmags[index,2][0]
				magerr = xysmags[index,3][0]

			ids.append(id)
			xs.append(x)
			ys.append(y)
			mags.append(mag)
			magerrs.append(magerr)

			id += 1

		ids_col = Column(name='num',data=ids)
		xs_col = Column(name='x', data = xs)
		ys_col = Column(name='y', data = ys)
		mags_col = Column(name = 'mag', data = mags)
		magerrs_col = Column(name = 'magerr',data = magerrs)
		out_table = Table()
		out_table.add_columns([ids_col, xs_col, ys_col, mags_col, magerrs_col])
		if verbose>1:
			print out_table
		return out_table



	def __stdcal_template_link_std_obs_external_flt(self,flt):
		raise IOError("under construction...")

	def __stdcal_template_link_std_obs_grmatch_flt(self,flt):

		tpl_imgkey = self.templates[flt]
		match_ret = self.__stdcal_template_link_std_obs_grmatch_single_image(tpl_imgkey,flt)

		return match_ret

	def __stdcal_template_link_std_obs_grmatch_single_image(self, imgkey,flt,single_use= False):
		'''
		match standard stars with sources detected on given image with grmatch

		INPUTS:
			imgkey:
			flt:
		'''
		imgkey_s = imgkey.split('.')[0]

		if not single_use:
			ref_list = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file[flt])	#x y mag magerr
		else:
			ref_list = os.path.join(self.std_ref_dir, 'std_ref_' + imgkey_s + '.txt')

	        std_ref_and_obs_measurement_stars_file  = imgkey_s +'_std_'+flt+'.match'
		match_output  = os.path.join(self.std_ref_dir, std_ref_and_obs_measurement_stars_file)
		trans_fitting = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref.coef')

		if self.photometry_method == 'apphot':
			input_list = os.path.join(self.aperture_photometry_dir, imgkey_s + self.apphot_ret_file_suffix)
			#column names in ref_list and input_list: xpos, ypos, mag,mag_err
		elif self.photometry_method == 'psfphot':
			input_list = self.__get_internal_psfphot_mag_file(imgkey)
		else:
			raise IOError("Invalid input for photometry_method")

		if (not os.path.exists(match_output)) or self.renew_std_ref_match:
			self.__delete_file_if_exist(match_output)
			self.__delete_file_if_exist(trans_fitting)

			self.fitsh_grmatch_type = 'point'
			matched_table = self.__fitsh_grmatch(ref_list, input_list, match_output, trans_output=trans_fitting)

			match_result_plot = self.stdref_tpl_match_result_display
			if match_result_plot:
				if not os.path.exists(match_output):
					raise IOError("grmatch failure...")
				trans_output = os.path.join(self.std_ref_dir, imgkey_s + '_tpl_stdref.trans')
				self.__fitsh_grtrans(match_output, trans_output, trans_fitting)
				trans_result = np.loadtxt(trans_output)
				ref_xys = trans_result[:,[0,1]]
				input_xys  = trans_result[:,[4,5]]
				plt.plot(ref_xys[:,0]-input_xys[:,0],ref_xys[:,1]-input_xys[:,1],'o')
				plt.show()
		else:
			print "std and obs match result exists %s"%match_output
			matched_table = np.loadtxt(match_output)

		return matched_table

	def __stdcal_template_flt(self,flt,link_stdref_obs_method='grmatch', std_obs_matched = False, stds_shift_x = 0, stds_shift_y = 0):
		'''
		Find one-one corresponding between standard reference stars and measurements in template image
		perform relative calibration to get the mag offset between standards and new measurements.

		Main functions involved:
			self.__stdcal_template_link_std_obs_flt
			self.__relative_calibration_single

		Inputs:
			flt:
			link_stdref_obs_method: See self.__stdcal_template_link_std_obs_flt for allowed options
			std_obs_matched: if True, then one mapping file containing std mags and new measurements is given
					on which the calibration will be based

		'''
		tpl_imgkey = self.templates[flt]
		tpl_key = tpl_imgkey.split('.')[0]

		if self.photometry_info[tpl_imgkey]['calmag'] != 99.99:
			if not self.renew_stdcal_mag:
				return

		if std_obs_matched:
			self.get_std_obs_match_filenames()
			std_obs_match = self.stdobs_match_files[tpl_imgkey]
			if not os.path.exists(std_obs_match):
				raise IOError("the std_obs_match file %s not exist"%std_obs_match)
			matched_ret = np.loadtxt(std_obs_match)
			M,N = matched_ret.shape
			if N<5:
				mags_ref_matched   = matched_ret[:,[0,1]]
				mags_input_matched = matched_ret[:,[2,3]]
			else:
				mags_ref_matched   = matched_ret[:,[2,3]]
				mags_input_matched = matched_ret[:,[6,7]]

		else:
			print link_stdref_obs_method
			self.__stdcal_template_link_std_obs_flt(flt,method=link_stdref_obs_method, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y )
			std_obs_match = os.path.join(self.std_ref_dir,self.flt_std_ref_and_obs_measurement_stars_file[flt])
			matched_ret   = np.loadtxt(std_obs_match)

			if self.reference_mag_min is not None:
				mag_faint_cut = self.reference_mag_min
				matched_ret = matched_ret[matched_ret[:,2]<mag_faint_cut,:]

			if self.reference_mag_max is not None:
				mag_bright_cut  = self.reference_mag_max
				matched_ret = matched_ret[matched_ret[:,2]>mag_bright_cut,:]

			if self.reference_magerr_cut is not None:
				magerr_cut_ref = self.reference_magerr_cut
				matched_ret = matched_ret[matched_ret[:,3]<magerr_cut_ref,:]

			if self.input_magerr_cut is not None:
				magerr_cut_input = self.input_magerr_cut
				matched_ret = matched_ret[matched_ret[:,7]<magerr_cut_input,:]

			ref_matched   = matched_ret[:,[0,1,2,3]]
			input_matched = matched_ret[:,[4,5,6,7]]
			mags_ref_matched = ref_matched[:,[2,3]]
			mags_input_matched = input_matched[:,[2,3]]

		tpl_instmag = self.photometry_info[tpl_imgkey]['instmag']
		tpl_instmagerr = self.photometry_info[tpl_imgkey]['instmagerr']
		cal_mag,cal_magerr, temp_future_use = self.__relative_calibration_single(mags_ref_matched,mags_input_matched,tpl_instmag)

		self.photometry_info[tpl_imgkey]['calmag'] = cal_mag
		self.photometry_info[tpl_imgkey]['calmagerr'] =  np.sqrt( cal_magerr**2 + tpl_instmagerr**2 )

	def stdcal_single_image(self, imgkey,flt, instmag =None, instmagerr=None, eyecheck=False, \
	link_stdref_obs_method='grmatch', std_catalog = 'apass', panstarrs_photmethod = 'PSF', \
	std_world2img_refimage=None, refimg_which_dir='raw_image', stds_shift_x=0, stds_shift_y=0, \
	std_obs_matched = False, verbose = 0, update_stdmag=True, auto_shift_xy = False, \
	renew_selfphot_stds = False, secondary_std_photrettable=None, secondard_std_photfile=None):
		'''
		This function directly calibration given image to standard stars

		Find one-one corresponding between standard reference stars and measurements in given image
		perform relative calibration to get the mag offset between standards and new measurements.


		Main functions involved:
			self.__stdcal_template_link_std_obs_flt
			self.__relative_calibration_single

		Inputs:
			imgkey:
			flt:
			instmag:
			instmagerr:
			eyecheck: whether do checking on ds9 display
			link_stdref_obs_method: See self.__stdcal_template_link_std_obs_flt for allowed options
			std_obs_matched: if True, then one mapping file containing std mags and new measurements is given
					on which the calibration will be based

			std_catalog: 'apass', 'panstarrs', '2mass', or from private photometry result for example 'selfphot,LCOGT,001,psfphot' or 'selfphot,LCOGT,001,010,psfphot'
				    where LCOGT is telescope code, 001 is image index, and PSF is photometry method for which apphot is another option
			panstarrs_photmethod: 'PSF' or 'AP'

			stds_shift_x:
			stds_shift_y:

			std_world2img_refimage: the image containing WCS information acoording which transform the standard stars from world coordinate to image coordinate


		'''
		print imgkey
		key_s = imgkey.split('.')[0]

		if self.photometry_info[imgkey]['calmag'] != 99.99:
			if not self.renew_stdcal_mag:
				if verbose>0:
					print "calibration result exists for %s"%imgkey
				return
			else:
				if verbose > 0:
					print "renew calibration result for %s"%imgkey

		if std_obs_matched:
			self.get_std_obs_match_filenames()
			std_obs_match = self.stdobs_match_files[imgkey]
			matched_ret = np.loadtxt(std_obs_match)
			M,N = matched_ret.shape
			if N<5:  #only magnitude pairs
				mags_ref_matched   = matched_ret[:,[0,1]]
				mags_input_matched = matched_ret[:,[2,3]]
			else:   # x, y, mag, magerr, x, y, mag,magerr
				mags_ref_matched   = matched_ret[:,[2,3]]
				mags_input_matched = matched_ret[:,[6,7]]

		else:
			#prepare the standard stars
			print link_stdref_obs_method

			if std_world2img_refimage is None:
				if self.current_telescope == 'Iowa':
					outimg = os.path.join(self.template_dir, imgkey)
					refimg_abs = self.Iowa_prepare_template_image_wcs_single(imgkey, outimg=outimg)
				else:
					refimg_abs = self.__get_internal_image(imgkey, which_dir=refimg_which_dir)
				external_refimg = None
			else:
				refimg_abs = std_world2img_refimage
				external_refimg = std_world2img_refimage

			if std_catalog == 'apass':
				stds = self.get_standard_reference_star_APASS(flt,wcs=True, single_use = True, refimg_key=imgkey, refimg = refimg_abs, verbose=verbose)
				if verbose:
					print stds
			elif std_catalog == 'panstarrs':
				stds = self.get_standard_reference_star_PanSTARRS(flt, photmethod = panstarrs_photmethod, reference_image = refimg_abs, single_use=True, savedatakey=imgkey, external_refimg = external_refimg, save_std_flt_file = True, save_std_flt_reg = True)
			elif std_catalog == 'refcat2':
				stds = self.get_standard_reference_star_ATLASRefcat2(flt, reference_image = refimg_abs, single_use=True, savedatakey=imgkey, external_refimg = external_refimg, save_std_flt_file = True, save_std_flt_reg = True)
			elif std_catalog == 'sdss':
				stds = self.get_standard_reference_star_SDSS(flt, reference_image = refimg_abs, single_use=True, savedatakey=imgkey, external_refimg = external_refimg, save_std_flt_file = True, save_std_flt_reg = True)
			elif std_catalog == '2mass':
				stds = self.get_standard_reference_star_2MASS(flt,wcs=True, reference_image=refimg_abs, single_use=True, savedatakey=imgkey)
			elif std_catalog.split(',')[0] == 'selfphot':
				std_catalog_segs = std_catalog.split(',')
				if len(std_catalog_segs)<4:
					if secondary_std_photrettable is None or secondard_std_photfile is None:
						print "if std_catalog=selfphot, secondary_std_photrettable and secondard_std_photfile are necessary inputs"
						raise ValueError('invalid input for std_catalog')
				elif len(std_catalog_segs) == 4:
					stdcatcode, telcode, imgindex, photmethod = std_catalog_segs
					imgindex2 = None
				elif len(std_catalog_segs) == 5:
					stdcatcode, telcode, imgindex, imgindex2, photmethod = std_catalog_segs
					imgindex2 = imgindex2+'.fits'
				else:
					raise ValueError('invalid input for std_catalog')

				refimgkey = imgindex+'.fits'

				refxymagfile = os.path.join(self.std_ref_dir, 'xymag_'+telcode+refimgkey.replace('fits', 'txt'))
				refradecmagfile = os.path.join(self.std_ref_dir, 'radecmag_'+telcode+refimgkey.replace('fits', 'txt'))
				bridge_image = os.path.join(self.std_ref_dir, telcode+refimgkey)
				if os.path.exists(bridge_image) and os.path.exists(refradecmagfile) and (not renew_selfphot_stds):
					print "selfphot standards %s already exist and will be used; please set renew_selfphot_stds=1 if you want to rebuild"
				else:
					print "now go to build the standard catalog from internal photometry"
					if secondary_std_photrettable is None or secondard_std_photfile is None:
						secondary_std_photrettable, secondard_std_photfile = self.__build_bridge_data_for_stdcal_from_selfphot_prepare(refimgkey, flt, telcode, photmethod,  verbose=1)
					self.__build_bridge_data_for_stdcal_from_selfphot(refimgkey, flt, secondary_std_photrettable, secondard_std_photfile, telcode, secondary_imgkey=imgindex2, saverefimg=bridge_image, xymagfile=refxymagfile, radecmagfile=refradecmagfile, verbose=1)
				if std_world2img_refimage is None:
					refimg_abs = bridge_image
					external_refimg = bridge_image
				stds = self.__get_standard_reference_star_selfphot(refradecmagfile, refimg_abs, single_use=1, savedatakey=imgkey, title_text=imgkey)
			else:
				raise ValueError('std catalog %s not supported'%std_catalog)

			if auto_shift_xy:
				inputimage = self.images[imgkey]
				target_xy_wcs = self.__get_xy_on_image_from_wcsinfo(inputimage)
				stds_shift_x =  self.photometry_info[imgkey]['x']- target_xy_wcs[0]
				stds_shift_y =  self.photometry_info[imgkey]['y']- target_xy_wcs[1]
				if verbose:
					print "(x, y) transformation from (RA, Dec) = (%s, %s)"%(target_xy_wcs[0], target_xy_wcs[1])
					print "(x, y) target coordinate on the image = (%s, %s)"%(self.photometry_info[imgkey]['x'], self.photometry_info[imgkey]['y'])
					print "%s stds shift x: %s; shift y:%s"%(imgkey, stds_shift_x, stds_shift_y)

			matched_ret = self.__stdcal_link_std_obs_single_image(imgkey,flt,method=link_stdref_obs_method, single_use = True, stds_shift_x = stds_shift_x, stds_shift_y = stds_shift_y, external_refimg=external_refimg, image_dir =refimg_which_dir)

			if self.reference_mag_min is not None:
				mag_faint_cut = self.reference_mag_min
				matched_ret = matched_ret[matched_ret[:,2]<mag_faint_cut,:]
			if self.reference_mag_max is not None:
				mag_bright_cut  = self.reference_mag_max
				matched_ret = matched_ret[matched_ret[:,2]>mag_bright_cut,:]
			if self.reference_magerr_cut is not None:
				magerr_cut_ref = self.reference_magerr_cut
				matched_ret = matched_ret[matched_ret[:,3]<magerr_cut_ref,:]
			if self.input_magerr_cut is not None:
				magerr_cut_input = self.input_magerr_cut
				matched_ret = matched_ret[matched_ret[:,7]<magerr_cut_input,:]

			ref_matched   = matched_ret[:,[0,1,2,3]]
			input_matched = matched_ret[:,[4,5,6,7]]
			mags_ref_matched = ref_matched[:,[2,3]]
			mags_input_matched = input_matched[:,[2,3]]

		if instmag is None:
			instmag = self.photometry_info[imgkey]['instmag']
		if instmagerr is None:
			instmagerr = self.photometry_info[imgkey]['instmagerr']

		offset_method = self.cal_offset_method
		sigclipping = self.sigclip_first_preoffset
		if eyecheck:
			ref_matched_temp = Table(ref_matched, names=['x','y','m','merr'])
			ref_matched_temp, selectedindexs = self.select_or_delete_rows_from_given_table_by_pick_on_ds9_display(refimg_abs, ref_matched_temp, xcol='x', ycol='y', mode='d', match_criteria = 5, coordinate_type='image', radec_deg=True, circle_radius=10)
			saveindexs = [i for i in np.arange(len(ref_matched_temp)) if i not in selectedindexs]
			mags_ref_matched   = mags_ref_matched[saveindexs]
			mags_input_matched = mags_input_matched[saveindexs]
		if verbose:
			print mags_ref_matched
			print mags_input_matched
		cal_mag,cal_magerr, temp_future_use = self.__relative_calibration_single(mags_ref_matched,mags_input_matched, instmag, offset_method=offset_method, sigclipping_first = sigclipping, title_text=imgkey)

		if update_stdmag:
			self.photometry_info[imgkey]['calmag'] = cal_mag
			self.photometry_info[imgkey]['calmagerr'] = np.sqrt( cal_magerr**2 + instmagerr**2 )

		return cal_mag, cal_magerr

	def __get_tpl_image_size(self,flt, tplimg=None):
		'''
		get the image size
		'''
		if tplimg is not None:
			NX,NY = self.__get_fits_image_size(tplimg)
		else:
			cal_template = self.templates_after_astrometry[flt]
			image  = fits.open(cal_template)
			header = image['PRIMARY'].header
			data = image['PRIMARY'].data

			#LCOGT use new CCD with different NAXIS1 and NAXIS2
			NX = None
			NY = None

			self.__find_template_imagekey()
        		ori_template = os.path.join(self.raw_image_dir,self.templates[flt])
			hdu_ori = fits.open(ori_template)
			header_ori = hdu_ori['PRIMARY'].header

			try:
				NX = header_ori['NAXIS1']
				NY = header_ori['NAXIS2']
			except:
				print "image size not found in the original tempalte image"

			print NX,NY


			if NX is None or NY is None:
				try:
					NX = header['NAXIS1']
					NY = header['NAXIS2']
				except:
					print "image size not found in the astrometry solution"

			if NX is None or NY is None:
				try:
					NX = int(self.base.telescopes_info[self.current_telescope]['nx'])
	        	                NY = int(self.base.telescopes_info[self.current_telescope]['ny'])
				except:
					raise ValueError("NAXIS1 and NAXIS2 not available in image header and please specify it though telescope infomation file")

			print NX,NY

		return NX,NY

	
	def get_standard_reference_star_SDSS(self, flt, flt2=None, \
	reference_image = None, single_use = False, savedatakey=None, external_refimg=None, \
	save_std_flt_file = True, save_std_flt_reg = True, wcsinfo_hdu='PRIMARY'):
		'''
		prepare SDSS std stars and the expected product is a file containing img_x, img_y, mag_flt, magerr_flt where img_x and img_y are image
		coordinate related to the provided reference image and mag_flt is magnitudes in filter of interest. The reference image is provided in
		different ways under different circumstances, which are listed below.

		Reference image:
			(1) single_use --> False: the reference image is related to self.templates[flt]
			(2) single_use --> True:  'reference_image' is provided from internal data of current telescope+sn, 						  for example '001.fits'
			(3) 'external_refimg' is not None then will be used at higher priorities than (1), (2)

		INPUTS:
			flt:
			flt2:

		'''
		#check the input filter
		flts_acceptable = ['gp','rp','ip','zp']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable
		if (flt2 is not None) and flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable

		colnames_dict = self.sdss_colnames_dict
		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['sdss'])
		if os.path.exists(std_ref_stars_file_abs):
			refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')
			#print refstars
		else:
			print "====> ATTENTION!!! SDSS catalog not available; try to get transformed magnitudes from other catalogs"
			std_ref_stars_file_abs_borrow1 = os.path.join(self.std_ref_dir,self.std_ref_stars_file['panstarrs'])
			std_ref_stars_file_abs_borrow2 = os.path.join(self.std_ref_dir,self.std_ref_stars_file['refcat2'])
			if os.path.exists(std_ref_stars_file_abs_borrow1):
				refstars_borrow1 = Table.read(std_ref_stars_file_abs_borrow1, format='ascii.csv')
				print "PanSTARRS catalog loaded"
				refstars_borrow = self.add_transformed_flt_mags_PanSTARRS(refstars_borrow1)
			elif os.path.exists(std_ref_stars_file_abs_borrow2):
				refstars_borrow2 = Table.read(std_ref_stars_file_abs_borrow2, format='ascii.csv')
				print "ATLAS refcat2 catalog loaded"
				refstars_borrow = self.add_transformed_flt_mags_ATLASRefcat2(refstars_borrow2)
			else:
				raise ValueError("no catalog avaliable for the transformation")
			#print refstars_borrow
			try:
				refstars_borrow['RAJ2000'].name = colnames_dict['RA']
				refstars_borrow['DEJ2000'].name = colnames_dict['Dec']
			except:
				print refstars_borrow.colnames
				raise KeyError('RAJ2000 and DEJ2000 not found')

			#colnames_sdss = np.unique(self.sdss_colnames_dict.values())
			#print refstars_borrow.colnames
			colnames_sdss = ['RA_ICRS', 'DE_ICRS', 'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag','zmag', 'e_zmag']

			for colname in colnames_sdss:
				print colname
				if colname in [colnames_dict['RA'], colnames_dict['Dec'], 'umag', 'e_umag']:
					continue
				try:
					colname_from = colname.replace('mag','SDSSmag')
					if colname in refstars_borrow.colnames:
						refstars_borrow.remove_column(colname)
					refstars_borrow[colname_from].name=colname
				except:
					raise KeyError('please check the colnames')

			refstars = refstars_borrow[colnames_sdss]

			#check the get the wcs info of the reference image
		if external_refimg is not None:
			if not os.path.isfile(external_refimg):
				raise IOError('external reference image %s does not exist'%external_refimg)
			reference_image = external_refimg
		else:
			if single_use:
				if reference_image is None:
					raise ValueError('the reference image refimg is required...')
				if not os.path.isfile(reference_image):
					raise IOError('reference image %s does not exist'%reference_image)
			else:
				reference_image = self.templates_after_astrometry[flt]

		refstars = self.add_image_coordinate_to_standard_catalog(refstars, reference_image, colnames_dict, wcsinfo_hdu=wcsinfo_hdu)
		#print refstars

		magcol = colnames_dict[flt]
		magerrcol = colnames_dict[flt+'err']

		cols_wanted = ['x', 'y', magcol, magerrcol]
		if flt2 is not None:
			magcol2 = colnames_dict[flt2]
			magerrcol2 = colnames_dict[flt2+'err']
			cols_wanted.append(magcol2)
			cols_wanted.append(magerrcol2)

		std_table = refstars[cols_wanted]
		#print std_table

		#select stds in science image field
		if external_refimg is not None:
			NX = 9999
			NY = 9999
			lowcut = -9999

		if reference_image is None:
			NX,NY = self.__get_tpl_image_size(flt)
			lowcut = 0
		else:
			NX,NY = self.__get_fits_image_size(reference_image)
			lowcut = 0
		print NX,NY

		mask1 = np.logical_and(std_table['x']>lowcut, std_table['x']<NX)
		mask2 = np.logical_and(std_table['y']>lowcut, std_table['y']<NY)
		std_table = std_table[mask1*mask2]
		#print std_table

		if self.sdss_mag_faint_cut is not None:
			faint_mag_cut = self.sdss_mag_faint_cut
			std_table = std_table[std_table[colnames_dict[flt]]<faint_mag_cut]
		if self.sdss_mag_saturate_cut is not None:
			saturate_mag_cut = self.sdss_mag_saturate_cut
			std_table = std_table[std_table[colnames_dict[flt]]>saturate_mag_cut]

		if std_table.masked:
			mask = np.array([True]*len(std_table))
			for colname in std_table.colnames:
				mask = mask*(~std_table[colname].mask)
			std_table = std_table[mask]#remove rows with masked magnitudes

		if single_use:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key = reference_image.split('/')[-1]
				savedata_key_s = savedata_key.split('.')[0]

		if save_std_flt_file:
			if not single_use:
				std_ref_flt = self.flt_std_ref_stars_file[flt]
				std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
			else:
				std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

			self.__delete_file_if_exist(std_ref_filename_abs)
			std_table.write(std_ref_filename_abs, format='ascii.fast_no_header')

		if save_std_flt_reg: #save the sources as a ds9  reg file
			if not single_use:
				source_regfile =  os.path.join(self.std_ref_dir, flt + '.reg')
			else:
				source_regfile = os.path.join(self.std_ref_dir, savedata_key_s + '.reg')
			xys = np.array([std_table['x'].data, std_table['y'].data]).transpose()
			create_ds9_region_file(xys, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std_table




	def get_standard_reference_star_PanSTARRS_dev(self, reference_image, flt, photmethod = 'PSF'):
		'''
		INPUTS:
			reference_image:
			flt:
			wcs:
			photmethod:
		'''
		std_table = self.get_standard_reference_star_PanSTARRS(flt, photmethod = photmethod, reference_image = reference_image, save_std_flt_file = False, save_std_flt_reg = False)

		if flt == 'gp':
			self.panstarrs_stds_gp = std_table
		elif flt == 'rp':
			self.panstarrs_stds_rp = std_table
		elif flt == 'ip':
			self.panstarrs_stds_ip = std_table
		elif flt == 'zp':
			self.panstarrs_stds_zp = std_table
		elif flt == 'yp':
			self.panstarrs_stds_yp = std_table
		else:
			raise ValueError("flt %s not supported"%flt)

	def add_transformed_flt_mags_PanSTARRS(self, refstars, photmethod='PSF'):
		'''

		Transformation:
			PS1 results from https://iopscience.iop.org/article/10.1088/0004-637X/750/2/99/pdf
			B = g_P1 + 0.212 + 0.556* (g - r)_P1 + 0.034* (g - r)_P1**2  (+/-0.032)

			x | y | A_0 | A_1 | A_2 | +or- | B_0 | B_1 | +or-
			(g - r)_P1 (B - g_P1) 0.212 0.556 0.034 0.032 0.213 0.587 0.034
			(g - r)_P1 (V - r_P1) 0.005 0.462 0.013 0.012 0.006 0.474 0.012
			(g - r)_P1 (R_C - r_P1) -0.137 -0.108 -0.029 0.015 -0.138 -0.131 0.015
			(g - r)_P1 (I_C - i_P1) -0.366 -0.136 -0.018 0.017 -0.367 -0.149 0.016
			(g - r)_P1 (V - w_P1) -0.021 0.299 0.187 0.025 -0.011 0.439 0.035
			(g - r)_P1 (V - g_P1) 0.005 -0.536 0.011 0.012 0.006 -0.525 0.012

			(g - r)_P1 (gSDSS - gP1) 0.013 0.145 0.019 0.008 0.014 0.162 0.009
			(g - r)_P1 (rSDSS - rP1) -0.001 0.004 0.007 0.004 -0.001 0.011 0.004
			(g - r)_P1 (iSDSS - iP1) -0.005 0.011 0.010 0.004 -0.004 0.020 0.005
			(g - r)_P1 (zSDSS - zP1) 0.013 -0.039 -0.012 0.010 0.013 -0.050 0.010
			(g - r)_P1 (zSDSS - yP1) -0.031 0.111 0.004 0.024 -0.031 0.115 0.024
			(g - r)_P1 (rSDSS - wP1) -0.024 -0.149 0.155 0.018 -0.016 -0.029 0.031
		'''

		#Below are for catalog obtained from MAST

		if self.panstarrs_catalog_method == 1:
			if photmethod == 'PSF':
				gmags = refstars['gMeanPSFMag']
				gmagerrs = refstars['gMeanPSFMagErr']
				rmags = refstars['rMeanPSFMag']
				rmagerrs = refstars['rMeanPSFMagErr']
				imags = refstars['iMeanPSFMag']
				imagerrs = refstars['iMeanPSFMagErr']
				zmags = refstars['zMeanPSFMag']
				zmagerrs = refstars['zMeanPSFMagErr']
			elif photmethod == 'AP':
				gmags = refstars['gMeanApMag']
				gmagerrs = refstars['gMeanApMagErr']
				rmags = refstars['rMeanApMag']
				rmagerrs = refstars['rMeanApMagErr']
				imags = refstars['iMeanApMag']
				imagerrs = refstars['iMeanApMagErr']
				zmags = refstars['zMeanApMag']
				zmagerrs = refstars['zMeanApMagErr']
			else:
				raise ValueError("invalid input for photmethod...")
		elif self.panstarrs_catalog_method == 2:
			if photmethod == 'AP':
				raise ValueError('AP magnitudes not avaiable for self.panstarrs_catalog_method == 2')
			gmags = refstars['gmag']
			gmagerrs = refstars['e_gmag']
			rmags = refstars['rmag']
			rmagerrs = refstars['e_rmag']
			imags = refstars['imag']
			imagerrs = refstars['e_imag']
			zmags = refstars['zmag']
			zmagerrs = refstars['e_zmag']
		else:
			raise ValueError('self.panstarrs_catalog_method valid option: 1,2 ')

		#0.212 0.556 0.034 0.032
		Bmags = gmags + 0.212 + 0.556* (gmags - rmags) + 0.034* (gmags - rmags)**2
		Bmagerrs = np.sqrt(0.032**2 + (1+0.556**2+(2*0.034)**2)*gmagerrs**2 + (0.556**2+(2*0.034)**2)*rmagerrs**2)

		#0.005 0.462 0.013 0.012
		Vmags = rmags + 0.005 + 0.462* (gmags - rmags) + 0.013* (gmags - rmags)**2
		Vmagerrs = np.sqrt(0.012**2 + (1+0.462**2+(2*0.013)**2)*rmagerrs**2 + (0.462**2+(2*0.013)**2)*gmagerrs**2)

		#-0.137 -0.108 -0.029 0.015
		Rmags = rmags -0.137 - 0.108* (gmags - rmags) - 0.029* (gmags - rmags)**2
		Rmagerrs = np.sqrt(0.015**2 + (1+0.108**2+(2*0.029)**2)*rmagerrs**2 + (0.108**2+(2*0.029)**2)*gmagerrs**2)

		#-0.366 -0.136 -0.018 0.017
		Imags = imags -0.366 - 0.136* (gmags - rmags) - 0.018* (gmags - rmags)**2
		Imagerrs = np.sqrt(0.017**2 + imagerrs**2 + (0.136**2+(2*0.018)**2)*rmagerrs**2 + (0.136**2+(2*0.018)**2)*gmagerrs**2)
		
		#0.013 0.145 0.019 0.008
		gSDSSmags = gmags + 0.013 + 0.145* (gmags - rmags) + 0.019* (gmags - rmags)**2
		gSDSSmagerrs = np.sqrt(0.008**2 + (1+0.145**2+(2*0.019)**2)*gmagerrs**2 + (0.145**2+(2*0.019)**2)*rmagerrs**2)
		
		#-0.001 0.004 0.007 0.004
		rSDSSmags = rmags - 0.001 + 0.004* (gmags - rmags) + 0.007* (gmags - rmags)**2
		rSDSSmagerrs = np.sqrt(0.004**2 + (0.004**2+(2*0.007)**2)*gmagerrs**2 + (1+0.004**2+(2*0.007)**2)*rmagerrs**2)

		#-0.005 0.011 0.010 0.004
		iSDSSmags = imags - 0.005 + 0.011* (gmags - rmags) + 0.010* (gmags - rmags)**2
		iSDSSmagerrs = np.sqrt(0.004**2 + imagerrs**2 + (0.011**2+(2*0.010)**2)*gmagerrs**2 + (0.011**2+(2*0.010)**2)*rmagerrs**2)

		#0.013 -0.039 -0.012 0.010
		zSDSSmags = zmags + 0.013 - 0.039* (gmags - rmags) - 0.012* (gmags - rmags)**2
		zSDSSmagerrs = np.sqrt(0.010**2 + zmagerrs**2 + (0.039**2+(2*0.012)**2)*gmagerrs**2 + (0.039**2+(2*0.012)**2)*rmagerrs**2)


		#Bmag_col = Column(data=Bmags, name='Bmag')
		#Vmag_col = Column(data=Vmags, name='Vmag')
		#Rmag_col = Column(data=Rmags, name='Rmag')
		#Imag_col = Column(data=Imags, name='Imag')
		Bmags.name = 'Bmag'
		Vmags.name='Vmag'
		Rmags.name='Rmag'
		Imags.name='Imag'
		gSDSSmags.name = 'gSDSSmag'
		rSDSSmags.name = 'rSDSSmag'
		iSDSSmags.name = 'iSDSSmag'
		zSDSSmags.name = 'zSDSSmag'
		if self.panstarrs_catalog_method == 1:
			#Bmagerr_col = Column(data=Bmagerrs, name='Bmagerr')
			#Vmagerr_col = Column(data=Vmagerrs, name='Vmagerr')
			#Rmagerr_col = Column(data=Rmagerrs, name='Rmagerr')
			#Imagerr_col = Column(data=Imagerrs, name='Imagerr')
			Bmagerrs.name='Bmagerr'
			Vmagerrs.name='Vmagerr'
			Rmagerrs.name='Rmagerr'
			Imagerrs.name='Imagerr'
			gSDSSmagerrs.name = 'gSDSSmagerr'
			rSDSSmagerrs.name = 'rSDSSmagerr'
			iSDSSmagerrs.name = 'iSDSSmagerr'
			zSDSSmagerrs.name = 'zSDSSmagerr'
		else:
			#Bmagerr_col = Column(data=Bmagerrs, name='e_Bmag')
			#Vmagerr_col = Column(data=Vmagerrs, name='e_Vmag')
			#Rmagerr_col = Column(data=Rmagerrs, name='e_Rmag')
			#Imagerr_col = Column(data=Imagerrs, name='e_Imag')
			Bmagerrs.name='e_Bmag'
			Vmagerrs.name='e_Vmag'
			Rmagerrs.name='e_Rmag'
			Imagerrs.name='e_Imag'
			gSDSSmagerrs.name = 'e_gSDSSmag'
			rSDSSmagerrs.name = 'e_rSDSSmag'
			iSDSSmagerrs.name = 'e_iSDSSmag'
			zSDSSmagerrs.name = 'e_zSDSSmag'

		#refstars.add_columns([Bmag_col, Bmagerr_col, Vmag_col, Vmagrr_col, Rmag_col, Rmagerr_col, Imag_col, Imagerr_col])
		refstars.add_columns([Bmags, Bmagerrs, Vmags, Vmagerrs, Rmags, Rmagerrs, Imags, Imagerrs, gSDSSmags, gSDSSmagerrs,  rSDSSmags, rSDSSmagerrs, iSDSSmags, iSDSSmagerrs, zSDSSmags, zSDSSmagerrs])

		return refstars

	def add_image_coordinate_to_standard_catalog(self, refstars, image, colnames_dict, wcsinfo_hdu='PRIMARY'):
		'''
		add image coordinate x,y to PanSTARRS catalog
		'''
		hdu  = fits.open(image)
		header = hdu[wcsinfo_hdu].header

		w = WCS(header)
		ra  = refstars[colnames_dict['RA']].data
		dec = refstars[colnames_dict['Dec']].data
		world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
		pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)

		#fits image treat first pixel as (1,1)
		#to make the star pixel position consistent with that of the result from sfind (check here)
		#pix = pix-1  # uncomment this if source detection with sfind 
		xy_table = Table(pix,names=['x','y'])
		refstars_new = hstack((xy_table,refstars))

		return refstars_new


	def filter_sources_PanSTARRS(self, refstars, colnames_dict):
		'''
		select specific sources from PanSTARRS catalog
		'''
		#stds filtering to get stars with good magnitudes
		#the current values can be modified to have real filtering effect
		mag_up = 99
		magerr_up = 99
		mag_low = -100
		magerr_low = -100
		refstars = refstars[refstars[colnames_dict[flt]]>mag_low]
		refstars = refstars[refstars[colnames_dict[flt+'err']]>magerr_low]

		refstars = refstars[refstars[colnames_dict[flt]]<mag_up]
		refstars = refstars[refstars[colnames_dict[flt+'err']]<magerr_up]
		#print refstars

		if self.panstarrs_catalog_method ==1:
			if self.panstarrs_min_nx is not None:
				colname = 'n'+flt[0]
				refstars = refstars[refstars[colname]>self.panstarrs_min_nx]
			if self.panstarrs_min_nstackDetections is not None:
				refstars = refstars[refstars['nStackDetections'] > self.panstarrs_min_nstackDetections]

		return refstars

	def get_expanded_refined_filtered_stds_PanSTARRS(self, refimg, add_transformed_mags = True, add_image_coordinate=True, wcsinfo_hdu='PRIMARY'):
		'''
		here we go
		'''
		#get the related panstarrs catalog data
		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['panstarrs'])
		refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')
		refstars = refstars.filled(fill_value=99.99)

		if self.panstarrs_catalog_method == 1:
			if photmethod == 'PSF':
				colnames_dict = self.panstarrs_colnames_dict_method1
			elif photmethod == 'AP':
				colnames_dict = self.panstarrs_colnames_dict_method1_AP
			else:
				raise ValueError('photmethod has to be PSF or AP')
		elif self.panstarrs_catalog_method == 2:
			colnames_dict = self.panstarrs_colnames_dict
		else:
			raise ValueError('self.panstarrs_catalog_method valid option: 1, 2')

		refstars = self.add_transformed_flt_mags_PanSTARRS(refstars)
		refstars = self.add_image_coordinate_to_standard_catalog(refstars, refimg, colnames_dict, wcsinfo_hdu=wcsinfo_hdu)

		return refstars


	def get_standard_reference_star_PanSTARRS(self, flt, flt2=None, photmethod = 'PSF', \
	reference_image = None, single_use = False, savedatakey=None, external_refimg=None, \
	save_std_flt_file = True, save_std_flt_reg = True, wcsinfo_hdu='PRIMARY'):
		'''
		prepare PanSTARRS std stars and the expected product is a file containing img_x, img_y, mag_flt, magerr_flt where img_x and img_y are image
		coordinate related to the provided reference image and mag_flt is magnitudes in filter of interest. The reference image is provided in
		different ways under different circumstances, which are listed below.

		Reference image:
			(1) single_use --> False: the reference image is related to self.templates[flt]
			(2) single_use --> True:  'reference_image' is provided from internal data of current telescope+sn, 						  for example '001.fits'
			(3) 'external_refimg' is not None then will be used at higher priorities than (1), (2)

		INPUTS:
			flt:
			flt2:
			SDSS: whether or note transform the g,r,i,z magnitudes to SDSS system
			photmethod: 'PSF' or 'AP'; 'PSF' not supported for self.panstarrs_catalog_method=2

		'''
		#check the input filter
		flts_acceptable = ['gp','rp','ip','zp','yp', 'B', 'V', 'R', 'I']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable
		if (flt2 is not None) and flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable

		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['panstarrs'])
		refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')
		#print refstars

		if self.panstarrs_catalog_method == 1:
			if photmethod == 'PSF':
				colnames_dict = self.panstarrs_colnames_dict_method1
			elif photmethod == 'AP':
				colnames_dict = self.panstarrs_colnames_dict_method1_AP
			else:
				raise ValueError('photmethod has to be PSF or AP')
		elif self.panstarrs_catalog_method == 2:
			colnames_dict = self.panstarrs_colnames_dict
		else:
			raise ValueError('self.panstarrs_catalog_method valid option: 1, 2')

		refstars = self.add_transformed_flt_mags_PanSTARRS(refstars)
		#print refstars
		#check the get the wcs info of the reference image
		if external_refimg is not None:
			if not os.path.isfile(external_refimg):
				raise IOError('external reference image %s does not exist'%external_refimg)
			reference_image = external_refimg
		else:
			if single_use:
				if reference_image is None:
					raise ValueError('the reference image refimg is required...')
				if not os.path.isfile(reference_image):
					raise IOError('reference image %s does not exist'%reference_image)
			else:
				reference_image = self.templates_after_astrometry[flt]

		refstars = self.add_image_coordinate_to_standard_catalog(refstars, reference_image, colnames_dict, wcsinfo_hdu=wcsinfo_hdu)
		#print refstars

		magcol = colnames_dict[flt]
		magerrcol = colnames_dict[flt+'err']

		cols_wanted = ['x', 'y', magcol, magerrcol]
		if flt2 is not None:
			magcol2 = colnames_dict[flt2]
			magerrcol2 = colnames_dict[flt2+'err']
			cols_wanted.append(magcol2)
			cols_wanted.append(magerrcol2)

		std_table = refstars[cols_wanted]
		#print std_table

		#select stds in science image field
		if external_refimg is not None:
			NX = 9999
			NY = 9999
			lowcut = -9999

		if reference_image is None:
			NX,NY = self.__get_tpl_image_size(flt)
			lowcut = 0
		else:
			NX,NY = self.__get_fits_image_size(reference_image)
			lowcut = 0
		print NX,NY

		mask1 = np.logical_and(std_table['x']>lowcut, std_table['x']<NX)
		mask2 = np.logical_and(std_table['y']>lowcut, std_table['y']<NY)
		std_table = std_table[mask1*mask2]
		#print std_table

		if self.panstarrs_mag_faint_cut is not None:
			faint_mag_cut = self.panstarrs_mag_faint_cut
			std_table = std_table[std_table[colnames_dict[flt]]<faint_mag_cut]
		if self.panstarrs_mag_saturate_cut is not None:
			saturate_mag_cut = self.panstarrs_mag_saturate_cut
			std_table = std_table[std_table[colnames_dict[flt]]>saturate_mag_cut]

		if std_table.masked:
			mask = np.array([True]*len(std_table))
			for colname in std_table.colnames:
				mask = mask*(~std_table[colname].mask)
			std_table = std_table[mask]#remove rows with masked magnitudes

		if single_use:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key = reference_image.split('/')[-1]
				savedata_key_s = savedata_key.split('.')[0]

		if save_std_flt_file:
			if not single_use:
				std_ref_flt = self.flt_std_ref_stars_file[flt]
				std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
			else:
				std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

			self.__delete_file_if_exist(std_ref_filename_abs)
			std_table.write(std_ref_filename_abs, format='ascii.fast_no_header')

		if save_std_flt_reg: #save the sources as a ds9  reg file
			if not single_use:
				source_regfile =  os.path.join(self.std_ref_dir, flt + '.reg')
			else:
				source_regfile = os.path.join(self.std_ref_dir, savedata_key_s + '.reg')
			xys = np.array([std_table['x'].data, std_table['y'].data]).transpose()
			create_ds9_region_file(xys, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std_table




	def add_transformed_flt_mags_ATLASRefcat2(self, refstars):
		'''
		Transformation:
			PS1 results from https://iopscience.iop.org/article/10.1088/0004-637X/750/2/99/pdf
			B = g_P1 + 0.212 + 0.556* (g - r)_P1 + 0.034* (g - r)_P1**2  (+/-0.032)

			x | y | A_0 | A_1 | A_2 | +or- | B_0 | B_1 | +or-
			(g - r)_P1 (B - g_P1) 0.212 0.556 0.034 0.032 0.213 0.587 0.034
			(g - r)_P1 (V - r_P1) 0.005 0.462 0.013 0.012 0.006 0.474 0.012
			(g - r)_P1 (R_C - r_P1) -0.137 -0.108 -0.029 0.015 -0.138 -0.131 0.015
			(g - r)_P1 (I_C - i_P1) -0.366 -0.136 -0.018 0.017 -0.367 -0.149 0.016
			(g - r)_P1 (V - w_P1) -0.021 0.299 0.187 0.025 -0.011 0.439 0.035
			(g - r)_P1 (V - g_P1) 0.005 -0.536 0.011 0.012 0.006 -0.525 0.012
		'''
		gmags = refstars['gmag']
		gmagerrs = refstars['e_gmag']
		rmags = refstars['rmag']
		rmagerrs = refstars['e_rmag']
		imags = refstars['imag']
		imagerrs = refstars['e_imag']
		zmags = refstars['zmag']
		zmagerrs = refstars['e_zmag']

		Bmags = gmags + 0.212 + 0.556* (gmags - rmags) + 0.034* (gmags - rmags)**2
		Bmagerrs = np.sqrt(0.032**2 + (1+0.556**2+(2*0.034)**2)*gmagerrs**2 + (0.556**2+(2*0.034)**2)*rmagerrs**2)

		Vmags = rmags + 0.005 + 0.462* (gmags - rmags) + 0.013* (gmags - rmags)**2
		Vmagerrs = np.sqrt(0.012**2 + (1+0.462**2+(2*0.013)**2)*rmagerrs**2 + (0.462**2+(2*0.013)**2)*gmagerrs**2)

		Rmags = rmags -0.137 - 0.108* (gmags - rmags) - 0.029* (gmags - rmags)**2
		Rmagerrs = np.sqrt(0.015**2 + (1+0.108**2+(2*0.029)**2)*rmagerrs**2 + (0.108**2+(2*0.029)**2)*gmagerrs**2)

		Imags = imags -0.366 - 0.136* (gmags - rmags) - 0.018* (gmags - rmags)**2
		Imagerrs = np.sqrt(0.017**2 + imagerrs**2 + (0.136**2+(2*0.018)**2)*rmagerrs**2 + (0.136**2+(2*0.018)**2)*gmagerrs**2)

		#0.013 0.145 0.019 0.008
		gSDSSmags = gmags + 0.013 + 0.145* (gmags - rmags) + 0.019* (gmags - rmags)**2
		gSDSSmagerrs = np.sqrt(0.008**2 + (1+0.145**2+(2*0.019)**2)*gmagerrs**2 + (0.145**2+(2*0.019)**2)*rmagerrs**2)
		
		#-0.001 0.004 0.007 0.004
		rSDSSmags = rmags - 0.001 + 0.004* (gmags - rmags) + 0.007* (gmags - rmags)**2
		rSDSSmagerrs = np.sqrt(0.004**2 + (0.004**2+(2*0.007)**2)*gmagerrs**2 + (1+0.004**2+(2*0.007)**2)*rmagerrs**2)

		#-0.005 0.011 0.010 0.004
		iSDSSmags = imags - 0.005 + 0.011* (gmags - rmags) + 0.010* (gmags - rmags)**2
		iSDSSmagerrs = np.sqrt(0.004**2 + imagerrs**2 + (0.011**2+(2*0.010)**2)*gmagerrs**2 + (0.011**2+(2*0.010)**2)*rmagerrs**2)

		#0.013 -0.039 -0.012 0.010
		zSDSSmags = zmags + 0.013 - 0.039* (gmags - rmags) - 0.012* (gmags - rmags)**2
		zSDSSmagerrs = np.sqrt(0.010**2 + zmagerrs**2 + (0.039**2+(2*0.012)**2)*gmagerrs**2 + (0.039**2+(2*0.012)**2)*rmagerrs**2)

		Bmags.name = 'Bmag'
		Vmags.name='Vmag'
		Rmags.name='Rmag'
		Imags.name='Imag'
		Bmagerrs.name='e_Bmag'
		Vmagerrs.name='e_Vmag'
		Rmagerrs.name='e_Rmag'
		Imagerrs.name='e_Imag'

		gSDSSmags.name = 'gSDSSmag'
		rSDSSmags.name = 'rSDSSmag'
		iSDSSmags.name = 'iSDSSmag'
		zSDSSmags.name = 'zSDSSmag'
		gSDSSmagerrs.name = 'e_gSDSSmag'
		rSDSSmagerrs.name = 'e_rSDSSmag'
		iSDSSmagerrs.name = 'e_iSDSSmag'
		zSDSSmagerrs.name = 'e_zSDSSmag'

		refstars.add_columns([Bmags, Bmagerrs, Vmags, Vmagerrs, Rmags, Rmagerrs, Imags, Imagerrs, gSDSSmags, gSDSSmagerrs,  rSDSSmags, rSDSSmagerrs, iSDSSmags, iSDSSmagerrs, zSDSSmags, zSDSSmagerrs])

		return refstars


	def get_standard_reference_star_ATLASRefcat2(self, flt, flt2=None, \
	reference_image = None, single_use = False, savedatakey=None, external_refimg=None, \
	save_std_flt_file = True, save_std_flt_reg = True, wcsinfo_hdu='PRIMARY'):
		'''
		prepare ATLAS-Refcat2 std stars and the expected product is a file containing img_x, img_y, mag_flt, magerr_flt where img_x and img_y are image
		coordinate related to the provided reference image and mag_flt is magnitudes in filter of interest. The reference image is provided in
		different ways under different circumstances, which are listed below.

		Reference image:
			(1) single_use --> False: the reference image is related to self.templates[flt]
			(2) single_use --> True:  'reference_image' is provided from internal data of current telescope+sn, for example '001.fits'
			(3) 'external_refimg' is not None then will be used at higher priorities than (1), (2)

		INPUTS:
			flt:
			wcs: whether the input reference image has WCS info in header
			r2R_flts: gr or ri
			i2I_flts: ri or iz

		'''
		#check the input filter
		flts_acceptable = ['gp','rp','ip','zp','yp', 'B', 'V', 'R', 'I']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable
		if (flt2 is not None) and flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable

		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['refcat2'])
		refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')
		colnames_dict = self.refcat2_colnames_dict	

		refstars = self.add_transformed_flt_mags_ATLASRefcat2(refstars)
		#check and get the wcs info of the reference image
		if external_refimg is not None:
			if not os.path.isfile(external_refimg):
				raise IOError('external reference image %s does not exist'%external_refimg)
			reference_image = external_refimg
		else:
			if single_use:
				if reference_image is None:
					raise ValueError('the reference image refimg is required...')
				if not os.path.isfile(reference_image):
					raise IOError('reference image %s does not exist'%reference_image)
			else:
				reference_image = self.templates_after_astrometry[flt]

		refstars = self.add_image_coordinate_to_standard_catalog(refstars, reference_image, colnames_dict, wcsinfo_hdu=wcsinfo_hdu)
		#print refstars

		magcol = colnames_dict[flt]
		magerrcol = colnames_dict[flt+'err']
		cols_wanted = ['x', 'y', magcol, magerrcol]
		if flt2 is not None:
			magcol2 = colnames_dict[flt2]
			magerrcol2 = colnames_dict[flt2+'err']
			cols_wanted.append(magcol2)
			cols_wanted.append(magerrcol2)

		std_table = refstars[cols_wanted]
		#print std_table

		#select stds in science image field
		if external_refimg is not None:
			NX = 9999
			NY = 9999
			lowcut = -9999

		if reference_image is None:
			NX,NY = self.__get_tpl_image_size(flt)
			lowcut = 0
		else:
			NX,NY = self.__get_fits_image_size(reference_image)
			lowcut = 0
		print NX,NY

		mask1 = np.logical_and(std_table['x']>lowcut, std_table['x']<NX)
		mask2 = np.logical_and(std_table['y']>lowcut, std_table['y']<NY)
		std_table = std_table[mask1*mask2]
		#print std_table

		if self.refcat2_mag_faint_cut is not None:
			faint_mag_cut = self.refcat2_mag_faint_cut
			std_table = std_table[std_table[colnames_dict[flt]]<faint_mag_cut]
		if self.refcat2_mag_saturate_cut is not None:
			saturate_mag_cut = self.refcat2_mag_saturate_cut
			std_table = std_table[std_table[colnames_dict[flt]]>saturate_mag_cut]

		if std_table.masked:
			mask = np.array([True]*len(std_table))
			for colname in std_table.colnames:
				mask = mask*(~std_table[colname].mask)
			std_table = std_table[mask]#remove rows with masked magnitudes

		if single_use:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key = reference_image.split('/')[-1]
				savedata_key_s = savedata_key.split('.')[0]

		if save_std_flt_file:
			if not single_use:
				std_ref_flt = self.flt_std_ref_stars_file[flt]
				std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
			else:
				std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

			self.__delete_file_if_exist(std_ref_filename_abs)
			std_table.write(std_ref_filename_abs, format='ascii.fast_no_header')

		if save_std_flt_reg: #save the sources as a ds9  reg file
			if not single_use:
				source_regfile =  os.path.join(self.std_ref_dir, flt + '.reg')
			else:
				source_regfile = os.path.join(self.std_ref_dir, savedata_key_s + '.reg')
			xys = np.array([std_table['x'].data, std_table['y'].data]).transpose()
			create_ds9_region_file(xys, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std_table



	def get_standard_reference_star_2MASS(self,flt,wcs=True, eliminate =True, reference_image=None, single_use=False, savedatakey=None, verbose=1):
		'''
		self.standards have multiple colnames; This function is used to extract image position (x,y) and single band magnitude

		INPUTS:
			eliminate: eliminate stars which are outside the field of view of the image

		OUTPUTS:
			output file contains single band standard stars, with columns of 'x','y','mag','magerr'
		'''
		#check the unput filter
		flts_acceptable = ['J','H','K']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable

		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['2mass'])
		refstars = Table.read(std_ref_stars_file_abs, format='ascii.csv')
		std_ref_big_file = self.flt_std_ref_stars_file_big[flt]
		std_ref_big_absfile = os.path.join(self.std_ref_dir, std_ref_big_file)
		self.__delete_file_if_exist(std_ref_big_absfile)

		if single_use:
			if reference_image is None:
				raise ValueError('the reference image refimg is required...')
			if not os.path.exists(reference_image):
				raise IOError('reference image %s not exist'%reference_image)
			refimg_key = reference_image.split('/')[-1]
			refimg_key_s = refimg_key.split('.')[0]
			if verbose:
				print "image %s will be used to relate the 2MASS catalog and observed sources"%reference_image
		else:
			reference_image  = self.templates_after_astrometry[flt]

		image  = fits.open(reference_image)
		header = image['PRIMARY'].header
		if wcs:
			NX,NY = self.__get_tpl_image_size(flt, tplimg=reference_image)
			w = WCS(header)
			ra  = refstars[self.twomass_colnames_dict['RA']].data.data	#refstars['ra'] is masked table column
			dec = refstars[self.twomass_colnames_dict['Dec']].data.data
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)

			#fits image treat first pixel as (1,1)
			#to make the star pixel position consistent with that of the result from sfind
			pix = pix-1
			xy_table = Table(pix,names=['img_x','img_y'])
			refstars_new = hstack((xy_table,refstars))
			refstars_new.write(std_ref_big_absfile,format='ascii.csv')

			cols_wanted = [self.twomass_colnames_dict[flt],self.twomass_colnames_dict[flt+'err']]
			mags_table = refstars[cols_wanted]

			#empty values exist in the table, fill the empty with 99.99
			if mags_table.masked:
				mags_table = mags_table.filled(fill_value=99.99)
			std_table = hstack((xy_table, mags_table))

			if eliminate:
				mask1 = np.logical_and(xy_table['img_x']>0, xy_table['img_x']<NX)
				mask2 = np.logical_and(xy_table['img_y']>0, xy_table['img_y']<NY)
				mask = mask1*mask2
				std_table = std_table[mask]
		else:
			raise IOError('under construction...')



		if not single_use:
			std_ref_flt = self.flt_std_ref_stars_file[flt]
			std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
		else:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key_s = refimg_key_s
			std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_filename_abs)
		std_table.write(std_ref_filename_abs, format='ascii.fast_no_header')

		#save the sources as a ds9  reg file
		source_regfile =  os.path.join(self.std_ref_dir, flt + '.reg')
		create_ds9_region_file(pix, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std_table

	def __get_standard_reference_star_selfphot(self, radecmagfile, refimg, \
	single_use=1, wcs=1, external_refimg =None,  verbose=1, savedatakey=None):
		'''
		INPUTS:
			radecmagfile:
			refimg:
			single_use:
			wcs: whether the refimg has wcs info in header;
			verbose:
		'''

		if single_use:
			if savedatakey is not None:
				savedata_key_s = savedatakey.split('.')[0]
			else:
				savedata_key = refimg.split('/')[-1]
				savedata_key_s = savedata_key.split('.')[0]
			if verbose:
				print "image %s will be used to relate the std catalog and observed sources"%refimg

		refstars = np.loadtxt(radecmagfile, skiprows=0)
		if verbose:
			print "the primitive std stars list are below:"
			print refstars

		#you can exclude standard stars with too big uncertainties
		if self.selfphot_magerr_up_cut is not None:
			refstars = refstars[refstars[:,3]<self.selfphot_magerr_up_cut]
			if verbose:
				print "stds with uncertainties >%s removed"%self.selfphot_magerr_up_cut

		if single_use:
			if refimg is None:
				raise ValueError('the reference image refimg is required...')
			elif not os.path.exists(refimg):
				raise IOError('reference image %s not exist'%refimg)
			else:
				image = fits.open(refimg)
				header = image['PRIMARY'].header
		else: #dead code block
			std_ref_big_file = self.flt_std_ref_stars_file_big[flt]
			std_ref_big_absfile = os.path.join(self.std_ref_dir, std_ref_big_file)
			self.__delete_file_if_exist(std_ref_big_absfile)
			cal_template = self.templates_after_astrometry[flt]
			image  = fits.open(cal_template)
			header = image['PRIMARY'].header

		if wcs:
			NX,NY = self.__get_fits_image_size(refimg)
			w = WCS(header)
			ra  = refstars[:,0]
			dec = refstars[:,1]
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)

			#fits image treat first pixel as (1,1)
			#to make the star pixel position consistent with that of the result from sfind
			pix = pix-1
			refstars_new = np.concatenate((pix,refstars),axis=1)

			mask  = np.logical_and(refstars_new[:,4]<90,refstars_new[:,4]<90)
			std = refstars_new[mask,:]
			std = std[:,[0,1,4,5]]

			if verbose:
				print "updated std stars list are below:"
				print std

		if self.selfphot_mag_faint_cut is not None:
			faint_mag_cut = self.selfphot_mag_faint_cut
			std  = std[std[:,2]<faint_mag_cut]
			if verbose:
				print "brighter than %s sources survived"%faint_mag_cut

		if self.selfphot_mag_saturate_cut is not None:
			saturate_mag_cut = self.selfphot_mag_saturate_cut
			std = std[std[:,2]>saturate_mag_cut]
			if verbose:
				print "fainter than %s sources survived"%saturate_mag_cut

		if wcs:
			#only the standard stars within the input image region are needed
			mask1 = np.logical_and(std[:,0]>0,std[:,0]<NX)
			mask2 = np.logical_and(std[:,1]>0,std[:,1]<NY)
			mask = mask1*mask2
			std = std[mask,:]


		if not single_use:
			std_ref_flt = self.flt_std_ref_stars_file[flt]
			std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
		else:
			std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+savedata_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_filename_abs)
		fid = open(std_ref_filename_abs,'w')
		np.savetxt(fid,std,fmt="%8.4f %8.4f %8.3f %8.3f")
		fid.close()

		#save the sources as a ds9  reg file
		if not single_use:
			source_regfile = os.path.join(self.std_ref_dir, flt + '.reg')
		else:
			source_regfile = os.path.join(self.std_ref_dir, savedata_key_s + '.reg')

		create_ds9_region_file(std, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)
		return std

	def __build_bridge_data_for_stdcal_from_selfphot_prepare(self, imgkey, flt, telcode, photmethod, verbose=1):
		'''
		find the corresponding file;
		reply on the consistent file structure
		'''
		srcimgdir = self.raw_image_dir.replace(self.current_telescope, telcode)
		srcretdir = self.result_dir.replace(self.current_telescope, telcode)

		if photmethod == 'apphot':
			suffix1 = '_AP.txt'	# photometry results table file
			photretdir = self.aperture_photometry_dir.replace(self.current_telescope, telcode)
		elif photmethod == 'psfphot':
			suffix1 = '_PSF.txt'
			psf_photometry_dir = self.__get_internal_psfphot_dir()
			photretdir = psf_photometry_dir.replace(self.current_telescope, telcode)
		else:
			raise ValueError('invalid input for photmethod')

		photrettable = os.path.join(srcretdir, '%s_photometry_info%s'%(self.current_sn, suffix1))
		photfile = os.path.join(photretdir, imgkey.replace('fits', photmethod))

		if not os.path.exists(photrettable):
			raise ValueError('The default file %s not exists'%photrettable)

		if not os.path.exists(photfile):
			raise ValueError('The default file %s not exists'%photfile)

		return photrettable, photfile

	def __build_bridge_data_for_stdcal_from_selfphot(self, imgkey, flt, photrettable, photfile, telcode,  secondary_imgkey=None, saverefimg=None, xymagfile=None, radecmagfile=None, verbose=1):
		'''
		prepare standard stars for calibration from photometry result for images of other telescope
		For R, I band, we need conversion from r,i band for which the secondary image from another filter is required

		INPUTS:
			imgkey: for example, '001.fits'
			telcode: for example, 'LCOGT', 'LT'
			photmethod: 'apphot', 'psfphot'
			saverefimg: the filename of the image to be saved
			xymagfile: the std catalog with format x,y,mag,magerr
			radecmagfile: the std catalog with format ra, dec, mag, magerr
		'''
		photret = Table.read(photrettable, format='ascii.fixed_width')
		imgphot = photret[photret['name']==imgkey]
		zptmag = imgphot['magzpt'][0]
		zptmagerr = imgphot['magzpterr'][0]
		flt1 = imgphot['flt'][0]
		if zptmag == 99.99 or zptmagerr == 99.99:
			instmag = imgphot['instmag'][0]
			instmagerr = imgphot['instmagerr'][0]
			calmag = imgphot['calmag'][0]
			calmagerr = imgphot['calmagerr'][0]
			zptmag = calmag - instmag
			zptmagerr = np.sqrt(calmagerr**2-instmagerr**2)

		xymags1 = np.loadtxt(photfile)
		xymags1[:,2] = xymags1[:,2] + zptmag
		xymags1[:,3] = np.sqrt(xymags1[:,3]**2 + zptmagerr**2)

		if secondary_imgkey is not None:
			photfile2 = photfile.replace(imgkey, secondary_imgkey)
			imgphot = photret[photret['name']==secondary_imgkey]
			zptmag = imgphot['magzpt'][0]
			zptmagerr = imgphot['magzpterr'][0]
			flt2 = imgphot['flt'][0]

			if zptmag == 99.99 or zptmagerr == 99.99:
				instmag = imgphot['instmag'][0]
				instmagerr = imgphot['instmagerr'][0]
				calmag = imgphot['calmag'][0]
				calmagerr = imgphot['calmagerr'][0]
				zptmag = calmag - instmag
				zptmagerr = np.sqrt(calmagerr**2-instmagerr**2)

			xymags2 = np.loadtxt(photfile2)
			xymags2[:,2] = xymags2[:,2] + zptmag
			xymags2[:,3] = np.sqrt(xymags2[:,3]**2 + zptmagerr**2)


			if self.selfphot_match_method_for_RI_stds == 'grmatch':
				raise ValueError("grmatch method for matching %s and %s under development"%(imgkey,secondary_imgkey))
				#grmatch_option = self.grmatch_option_obsimg2refimg
				#matched_table = self.__fitsh_grmatch(photfile, photfile2, match_output_filename, trans_fitting_filename, option=grmatch_option)
			elif self.selfphot_match_method_for_RI_stds == 'surrounding_search':
				tg_xys  = xymags2[:,0:2]
	        	        ref_xys = xymags1[:,0:2]
	        	        criteria = self.criteria_tpl_obs_match
        		        target_stars,indexmat,ref_stars,indexref = self.__crude_match(tg_xys,ref_xys,criteria)
				if len(indexmat) != 0:
			                match_ret = np.hstack((xymags1[indexref],xymags2[indexmat]))
				else:
					print "No match result obtained !!! You can:"
					print "Try to modify self.criteria_tpl_obs_match and see if improve"
					print "Try to modify self.selfphot_match_method_for_RI_stds to grmatch and run agian"
			else:
				raise ValueError('invalid input for self.selfphot_match_method_for_RI_stds')

			if flt == 'R':
				if flt1 == 'rp' and flt2 == 'ip':
					xymags = match_ret[:,[0,1,2,3]]
					rmag = match_ret[:,2]
					rmagerr = match_ret[:,3]
					imag = match_ret[:,6]
					imagerr = match_ret[:,7]
					Rmag = rmag - 0.2936*(rmag - imag) - 0.1439
					sigma = 0.0072
					Rmagerr = np.sqrt(sigma**2 + np.sqrt((0.2936*imagerr)**2+(0.7064*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error
					xymags[:,2] = Rmag
					xymags[:,3] = Rmagerr
				else:
					raise ValueError("when prepare for R band std, r band should be main filter for transformation data and i band be secondary")
			elif flt == 'I':
				if flt1 == 'ip' and flt2 == 'rp':

					if verbose:
						print "SDSS-i band stds:"
						print match_ret[:, [0,1,2,3]]
						print "SDSS-r band stds:"
						print match_ret[:,[4,5,6,7]]

					xymags = match_ret[:,[0,1,2,3]]
					imag = match_ret[:,2]
					imagerr = match_ret[:,3]
					rmag = match_ret[:,6]
					rmagerr = match_ret[:,7]
					Imag = rmag - 1.2444*(rmag - imag) - 0.3820
					sigma = 0.0078
					Imagerr = np.sqrt(sigma**2 + np.sqrt((1.2444*imagerr)**2+(0.2444*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error
					xymags[:,2] = Imag
					xymags[:,3] = Imagerr
				else:
					raise ValueError("when prepare for I band std, i band should be main filter for transformation data and r band be secondary")
			else:
				raise ValueError("std transformation only works for R and I")
		else:
			if flt != flt1:
				raise ValueError("The obtained filter doesn't match the desired one")
			xymags = xymags1

		if not xymagfile:
			savefile = os.path.join(self.std_ref_dir, 'xymag_'+telcode+imgkey.replace('fits', 'txt'))
		else:
			savefile = os.path.join(self.std_ref_dir, xymagfile)
		self.__delete_file_if_exist(savefile)
		fid = open(savefile,'w')
		np.savetxt(fid, xymags, fmt="%8.4f %8.4f %8.3f %8.3f")
		fid.close()

		if not saverefimg:
			saveimg =  os.path.join(self.std_ref_dir, telcode+imgkey)
		else:
			saveimg = os.path.join(self.std_ref_dir, saverefimg)
		self.__delete_file_if_exist(saveimg)
		if verbose:
			print "Image %s from telescope %s will be saved as %s"%(imgkey, telcode, saveimg)

		srcimgdir = self.raw_image_dir.replace(self.current_telescope, telcode)
		shutil.copy(os.path.join(srcimgdir, imgkey), saveimg)

		radecmags = xymags.copy()
		radecmags[:,0:2] = self.__image2world(saveimg, xymags[:,0:2], hduindex=None)
		if not radecmagfile:
			savefile2 = os.path.join(self.std_ref_dir, 'radecmag_'+telcode+imgkey.replace('fits', 'txt'))
		else:
			savefile2 = os.path.join(self.std_ref_dir, radecmagfile)
		self.__delete_file_if_exist(savefile2)
		fid2 = open(savefile2,'w')
		np.savetxt(fid2, radecmags, fmt="%9.6f %9.6f %8.3f %8.3f")
		fid2.close()



	def get_standard_reference_star_APASS(self, flt, wcs=True, single_use = False, refimg_key=None, refimg = None, center_ra = None, center_dec = None, center_X = None, center_Y = None, verbose=0):
		'''
		Currently the calibrated star are from APASS DR9 which have Johson BV band and SDSS gri band photometry
		And we extend the use to Cousin R I band by transforming the SDSS r i band to I R according to
		'https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Lupton2005'

		    R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
		    R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072

		    I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
		    I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063        (no z band in APASS)

		This function will create convert APASS standards from world coordinate to image coordinate.
		you will need one reference image to determine the transformation and region of interest.
		If 'single_use' is not True, then the template image for filter 'flt' will be used as the reference image;
		otherwise, 'regimg' is required

		INPUTS:
			flt:
			wcs: if True, the wcs information exists in reference image header and will be used
			single_use:
			refimg:
			center_ra:
			center_dec:
			center_X:
			center_Y:
		OUTPUTS:

		'''
		if single_use:
			if refimg is None:
				raise ValueError('the reference image refimg is required...')
			if not os.path.exists(refimg):
				raise IOError('reference image %s not exist'%refimg)
			if refimg_key is None:
				refimg_key = refimg.split('/')[-1]
			refimg_key_s = refimg_key.split('.')[0]
			if verbose:
				print "image %s will be used to relate the APASS catalog and observed sources"%refimg
		else:
			refimg  = self.templates_after_astrometry[flt]

		#check the input filter
		flts_acceptable = ['B','V','I','R','gp','rp','ip','g','r','i']
		if flt not in flts_acceptable:
			raise IOError,'For now only %s acceptable !' % flts_acceptable

		std_ref_stars_file_abs = os.path.join(self.std_ref_dir,self.std_ref_stars_file['apass'])
		if self.apass_catalog_method == 1:
			refstars = Table.read(std_ref_stars_file_abs, format='ascii')
		elif self.apass_catalog_method == 2:
			refstars = Table.read(std_ref_stars_file_abs, format='csv')
		else:
			raise ValueError('check')

		#you can exclude standard stars with insufficient observation
		if self.apass_mobs_min is not None:
			refstars = refstars[refstars['mobs']>self.apass_mobs_min]
			if verbose:
				print "APASS stds with mobs >%s selected"%mobs_min
		if self.apass_nobs_min is not None:
			refstars = refstars[refstars['nobs']>self.apass_nobs_min]
			if verbose:
				print "APASS stds with nobs >%s selected"%nobs_min

		hdu  = fits.open(refimg)
		header = hdu['PRIMARY'].header
		NX,NY = self.__get_fits_image_size(refimg)

		refstarcolnames = refstars.colnames
		if ('RAJ2000' in refstarcolnames) and ('DEJ2000' in refstarcolnames):
			RA_colname  = 'RAJ2000'
			Dec_colname = 'DEJ2000'
		elif ('RA' in refstarcolnames) and ('Dec' in refstarcolnames):
			RA_colname  = 'RA'
			Dec_colname = 'Dec'
		else:
			raise ValueError("catalog colnames are not supported yet...")
		ra  = refstars[RA_colname]
		dec = refstars[Dec_colname]
		if wcs:
			w = WCS(header)
			world = np.concatenate((ra.reshape(len(ra),1),dec.reshape(len(dec),1)),axis=1)
			pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)
			pix = pix-1 #fits image treat first pixel as (1,1),to make the star pixel position consistent with that of the result from sfind
			xs = pix[:,0]
			ys = pix[:,1]
		else:
			pixscale_deg = self.pixscale/60./60.
			if center_ra is None or center_dec is None:
				center_ra  = self.sn_ra_world_deg
				center_dec = self.sn_dec_world_deg
				print "target coordinate used as image center coordinate"

                        if center_X is None or center_Y is None:
                                if not single_use:
                                        self.__find_template_imagekey()
                                        tpl_imgkey = self.templates[flt]
                                        center_X = self.photometry_info[tpl_imgkey]['x']
                                        center_Y = self.photometry_info[tpl_imgkey]['y']
                                else:
                                        center_X = self.photometry_info[refimg_key]['x']
                                        center_Y = self.photometry_info[refimg_key]['y']
				print "target position (x,y) on image used as reference point"

			xs = -(ra - center_ra)/pixscale_deg  + X_cen
			ys = (dec - center_dec)/pixscale_deg + Y_cen


		xcol = Column(data=xs, name='x') # in pixel
		ycol = Column(data=ys, name='y')
		refstars.add_columns([xcol,ycol], [0,0])

		if refstars.masked:
			refstars = refstars.filled(99.999)


		if flt in ['I','R']: #transformation needed
			refstars = refstars[(refstars['r_mag']<50)*(refstars['i_mag']<50)]
			refstars = refstars[(refstars['r_mag']!=0)*(refstars['i_mag']!=0)]
			rmag = refstars['r_mag'].data
			rmagerr = refstars['e_r_mag'].data
                        imag = refstars['i_mag'].data
                        imagerr = refstars['e_i_mag'].data

			if flt == 'I':
				mag = rmag - 1.2444*(rmag - imag) - 0.3820
				magerr = np.sqrt(0.0078**2 + np.sqrt((1.2444*imagerr)**2+(0.2444*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error (sigma = 0.0078)
			else:
				mag = rmag - 0.2936*(rmag - imag) - 0.1439
				magerr = np.sqrt(0.0072**2 + np.sqrt((0.2936*imagerr)**2+(0.7064*rmagerr)**2))  #uncertainty come from two parts, ri measurements and transformation systematic error (sigma = 0.0072)
		else:
			if flt in ['B','V']:
				magcolname = flt+'mag'
			elif flt in ['g','r','i']:
				magcolname = flt+'_mag'
			else:
				magcolname = flt[0]+'_mag'
			magerrcolname = 'e_'+magcolname

			refstars = refstars[refstars[magcolname]<50]
			refstars = refstars[refstars[magcolname]!=0]
			if self.apass_remove_single_measurement:
				refstars = refstars[refstars[magerrcolname]!=0]
			mag = refstars[magcolname].data
			magerr = refstars[magerrcolname].data

		if not single_use:
			std_ref_big_absfile = os.path.join(self.std_ref_dir, self.flt_std_ref_stars_file_big[flt])
		else:
			std_ref_big_absfile = os.path.join(self.std_ref_dir, 'std_ref_whole_info_'+refimg_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_big_absfile)
		refstars.write(std_ref_big_absfile,format='csv')


		x = refstars['x'].data
		y = refstars['y'].data
		std = np.array([x,y,mag,magerr]).transpose()

		if self.apass_mag_faint_cut is not None:
			std  = std[std[:,2]<self.apass_mag_faint_cut]
			if verbose:
				print "brighter than %s sources survived"%faint_mag_cut

		if self.apass_mag_saturate_cut is not None:
			std = std[std[:,2]>self.apass_mag_saturate_cut]
			if verbose:
				print "fainter than %s sources survived"%saturate_mag_cut

		if self.apass_remove_single_measurement:
			std = std[std[:,3]>0]
			if verbose:
				print "measurements with magnitude uncertainty = 0 removed"

		mask = ((std[:,0]>0)*(std[:,0]<NX))*((std[:,1]>0)*(std[:,1]<NY))
		std = std[mask,:] #only the standard stars within the input image region are needed

		if not single_use:
			std_ref_flt = self.flt_std_ref_stars_file[flt]
			std_ref_filename_abs = os.path.join(self.std_ref_dir,std_ref_flt)
		else:
			std_ref_filename_abs = os.path.join(self.std_ref_dir, 'std_ref_'+refimg_key_s+'.txt' )

		self.__delete_file_if_exist(std_ref_filename_abs)
		fid = open(std_ref_filename_abs,'w')
		np.savetxt(fid,std,fmt="%8.4f %8.4f %8.3f %8.3f")
		fid.close()

		#save the sources as a ds9  reg file
		if not single_use:
			source_regfile = os.path.join(self.std_ref_dir, flt + '.reg')
		else:
			source_regfile = os.path.join(self.std_ref_dir, refimg_key_s + '.reg')

		create_ds9_region_file(std, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		return std


#Photometry result output
	def save(self, outfile=None):
		'''
		pickle the current photometry instance to self.result_dir
		'''
		if outfile is None:
			outfile = 'photometry.pickle'
		savetofile = os.path.join(self.result_dir, outfile)

		if os.path.exists(savetofile):
			print "%s already exists; leave now and let you check!"
			return

		print "the photometry object will be pickled to %s"%savetofile
		f = open(savetofile, 'w')
		pickle.dump(self, f)
		f.close()

	def save_results(self,save_filename = None, overwrite_or_create_initversion=False, create_newversion=True):
		'''
		convert photometry_info dictionary to table and save.
		INPUTS:
			save_filename:
			overwrite_or_create_initversion: if yes, override the initial version (the version without suffix number)
			create_newversion: if yes, currentversion+1
		'''
		self.__dict2table()
		self.photometry_record_table.sort('obstime')

		if save_filename is None:
			if overwrite_or_create_initversion:
				save_filename = self.result_table_file_init
			else:
				if create_newversion:
					version = self.currentversion + 1
				else:
					version = self.currentversion
				self.result_table_file = self.result_table_file_init + '.' + str(version)
				save_filename = self.result_table_file
		self.__delete_file_if_exist(save_filename)
		Table.write(self.photometry_record_table, save_filename, format='ascii.fixed_width')
		self.add_note('Save the current result to file %s'%save_filename)
		self.save_notes()

	def save_lc_data(self,flts='all', versionflag='', binned =False, hostfree=False, bin_width=0.5,remove_inbin_outlier = False, updatetable=1):
		'''
		INPUTS:
			binned: binning the data first before save?
			hostfree: subtract host flux before save?
		'''
		self.__find_template_imagekey(updatetable=updatetable)

		if flts == 'all':
			flts_save = self.templates.keys()
		elif isinstance(flts,str):
			flts_save = [flts]
		elif isinstance(flts,list):
			flts_save = flts
		else:
			raise IOError('Invalid input for flts...')

		for flt in flts_save:
			self.save_lc_data_flt(flt, versionflag= versionflag, binned = binned, hostfree=hostfree, bin_width=bin_width,remove_inbin_outlier = remove_inbin_outlier)

	def __get_lc_flt(self, flt):
		'''
		get light curve (t,mag,magerr) from self.photometry_info
		'''
		self.__dict2table()
		tpl_table = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		tpl_nondrop_table = self.__select_rows_from_table(tpl_table,'drop',0)

		lc_flt  = Table()
		lc_flt.add_columns([tpl_nondrop_table['obstime'],tpl_nondrop_table['calmag'],tpl_nondrop_table['calmagerr']])
		lc_flt.sort('obstime')

		self.lcs_raw[flt] = lc_flt

		return lc_flt

	def __binning_data(self, obstimes, mags, magerrs, bin_width=0.5, remove_inbin_outlier=False):
		'''
		Data binning
		'''
                obstimes_bined = []
                mags_bined = []
                magerrs_bined = []

                bined_already_index = []
                for i,obstime in enumerate(obstimes):
                        if i in bined_already_index:
                                continue

                        distances = obstimes - obstime
                        #indexs = find(np.abs(distances)< bin_width)
                        indexs = np.where(np.abs(distances)< bin_width)[0]
                        bined_already_index = np.append(bined_already_index,indexs)

                        obstime_bined = np.mean(obstimes[indexs])
                        mags_tobin = mags[indexs]
                        magerrs_tobin = magerrs[indexs]

                        #remove outlier measurements
                        if remove_inbin_outlier and len(mags_tobin)>1:
                                mags_tobin_nonoutlier,index_keep,stdev = sigma_clipping(mags_tobin,sig=2)
                                magerrs_tobin_nonoutlier = magerrs_tobin[index_keep]
                        else:
                                mags_tobin_nonoutlier = mags_tobin
                                magerrs_tobin_nonoutlier = magerrs_tobin
                        weights = 1./magerrs_tobin_nonoutlier**2
                        mag_bined = np.sum(weights*mags_tobin_nonoutlier) / np.sum(weights)
                        magerr_bined = 1./ np.sqrt(np.sum(weights))

                        obstimes_bined.append(obstime_bined)
			mags_bined.append(mag_bined)
                        magerrs_bined.append(magerr_bined)

		N = len(obstimes_bined)
		obstimes_bined_array = np.array(obstimes_bined).reshape(N,1)
		mags_bined_array    = np.array(mags_bined).reshape(N,1)
		magerrs_bined_array = np.array(magerrs_bined).reshape(N,1)
		lc_flt = np.hstack((obstimes_bined_array,mags_bined_array,magerrs_bined_array))

		return lc_flt

	def save_lc_data_flt(self,flt, versionflag='', binned=False, hostfree=False,bin_width=0.5,remove_inbin_outlier=False):
		'''
		save light curve, filename format: sn_name+'_'+telescope_code+'-'+flt+'-PSF[/AP]' + ['-'+versionflag]+ '.txt'
		INPUTS:
			versionflag:
		'''
		if hostfree:
			self.__get_host_mags_dict()
			mag_host = self.host_mags_dict[flt][0]
			magerr_host = self.host_mags_dict[flt][1]
			lc_flt = self.subtract_host_magnitude_flt(flt, mag_host, magerr_host = magerr_host, A_flt=0, display_result = False)
			flt = flt + '-HFree'
		else:
			lc_flt = self.__get_lc_flt(flt)

		obstimes = lc_flt['obstime']
                mags = lc_flt['calmag']
                magerrs = lc_flt['calmagerr']

		if binned:
			lc_flt = self.__binning_data(obstimes, mags,magerrs,bin_width=bin_width, remove_inbin_outlier=remove_inbin_outlier)
			flt = flt+'-ave'

		sn_name = self.current_sn
		telescope_code =self.base.telescopes_info[self.current_telescope]['code']
		if self.photometry_method == 'apphot':
			methodcode = '-AP'
		elif self.photometry_method == 'psfphot':
			methodcode = '-PSF'
		else:
			raise IOError("Invalid input for photometry_method...")

		if versionflag != '':
			versionflag = '-'+versionflag
		if flt in ['gp', 'rp', 'ip', 'zp']:
			flt = flt[0]

		lc_filename = sn_name+'_'+telescope_code+'-'+flt+methodcode+versionflag+'.txt'

		outfile = os.path.join(self.result_dir,lc_filename)
		np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")

		if self.base.save_data_local_dir != '' and os.path.exists(self.base.save_data_local_dir):
			lc_dir_archive = os.path.join(self.base.save_data_local_dir,sn_name)
			if not os.path.exists(lc_dir_archive):
				os.mkdir(lc_dir_archive)

			outfile = os.path.join(lc_dir_archive, lc_filename)
			np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")

	def __get_host_mags_dict(self):
		'''
		Read in the magnitudes of host galaxy from the externally prepared and well stored file
		'''
		host_mags_file  = os.path.join( self.std_ref_dir , self.host_mags_file )

		if not os.path.exists(host_mags_file):
			raise IOError("host magnitudes not available: %s"%host_mags_file)

		host_mags = np.loadtxt(host_mags_file, dtype= {'names':('flt','mag','magerr'), 'formats':('S8','f4','f4')})

		for host_flt_mag in host_mags:
			flt = host_flt_mag[0]
			mag  = host_flt_mag[1]
			magerr = host_flt_mag[2]
			self.host_mags_dict[flt] = [mag, magerr]

	def subtract_host_magnitude(self, flts='all', display_result =False ):
		'''
		For nucleus targets, if archive host magnitude exist, we can subtract the host flux to get rough estimation on transient light curve

		'''
		self.__get_host_mags_dict()

		if flts == 'all':
			self.__find_template_imagekey()
			flts = self.templates.keys()

		for flt in flts:
			host_mag = self.host_mags_dict[flt][0]
			host_magerr = self.host_mags_dict[flt][1]

			lc_hostfree = self.subtract_host_magnitude_flt(flt, host_mag, magerr_host = host_magerr, display_result = display_result)

			self.lcs_hostfree[flt]  = lc_hostfree

	def subtract_host_magnitude_flt(self, flt, mag_host, magerr_host = 0, A_flt=0, display_result = True):
		'''
		subtract host flux

		INPUTS:
			flt:
			mag_host:
			magerr_host:
			A_flt: extinction in filter 'flt'
		'''
		obstime_new = []
		mag_new = []
		magerr_new = []

		flts_Vega =  ['U','B','V']
		lc_flt = self.__get_lc_flt(flt)

		for lcp in lc_flt:
			obstime = lcp['obstime']
			mag     = lcp['calmag']
			magerr  = lcp['calmagerr']

			print obstime, mag, magerr

			if mag > mag_host:
				continue

			if flt in flts_Vega:
				mag_tot_AB  = Vega_AB_mag_convertion(mag,flt,mode='Bessell',direction='Vega2AB')
	              		mag_host_AB = Vega_AB_mag_convertion(mag_host,flt,mode='Bessell',direction='Vega2AB')
			else:
				mag_tot_AB  = mag
				mag_host_AB = mag_host


			if flt == 'gp':
				flt_temp = 'g'
			elif flt == 'rp':
				flt_temp = 'r'
			elif flt == 'ip':
				flt_temp = 'i'
			elif flt == 'zp':
				flt_temp = 'z'
			elif flt == 'yp':
				flt_temp = 'y'
			else:
				flt_temp = flt

			lamb_flt,flux_tot  = mag2flux(mag_tot_AB - A_flt, flt_temp)
                        lamb_flt,flux_host = mag2flux(mag_host_AB - A_flt,flt_temp)
                        flux_transient = flux_tot - flux_host

			magerr_tot = magerr
                        fluxerr_transient   = np.sqrt(magerr_tot**2 + magerr_host**2)*1.086*flux_transient
			fluxupper_transient = flux_transient + fluxerr_transient

			mag_transient    = flux2mag(flux_transient,flt_temp,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms')
			magerr_transient = mag_transient - flux2mag(fluxupper_transient,flt_temp,mode='wavelength',I_unit='cgs',wavelength_units='Angstroms')
			obstime_new.append(obstime)
			mag_new.append(mag_transient)
			magerr_new.append(magerr_transient)

		if display_result:
			plt.errorbar(obstime_new, mag_new, yerr=magerr_new, fmt='o')
			plt.gca().invert_yaxis()
			plt.xlabel('JD')
			plt.ylabel('mag')
			plt.title('%s band transient light curve'%flt)
                	plt.show()

		lc_hostfree_data = [obstime_new, mag_new, magerr_new]
		lc_hostfree_table = Table(data=lc_hostfree_data, names=['obstime', 'calmag', 'calmagerr'])

		return lc_hostfree_table

#global uncertainty view tools
	def get_single_star_SNR(self,flt,photometry_function = 'iraf'):
		'''
		simple statistics on SNR of the photometry targets in given band
		'''
		obstimes = []
		snrs = []

		for img in self.images.keys():
			if self.photometry_info[img]['drop'] != 0:
				continue
			if self.photometry_info[img]['flt'] != flt:
				continue

			obstime = self.photometry_info[img]['obstime']
			snr = self.__single_star_SNR(img,photometry_function=photometry_function)

			obstimes.append(obstime)
			snrs.append(snr)

		SNR_table = Table()
		obstime_col = Column(data=obstimes,name='obstime')
		snr_col = Column(data=snrs,name='snr')
		SNR_table.add_columns([obstime_col,snr_col])

		return SNR_table

	def __single_star_SNR(self,image_key, x=None, y=None, photometry_function='iraf', which_dir = 'raw_image', fixed_app_size = False):
		'''
		Estimate the signal to noise ratio of photometry target
		The SNR calculation from IRAF measurement is done as below:
		IRAF photometry results give mag and magerr, where merr = 1.087*(N/S) where N is noise and S is signal
		IRAF calculation works in the level of real measurement of counts or data number (DN) for source and noise.
		For a measurement of N counts, the sigma = sqrt(N*epadu)/epadu = sqrt(N/epadu)
		N = sqrt(flux/epadu+ area*stdev**2 + [area*stdev/sqrt(nsky)]**2) where flux is total number of counts excluding sky in the aperture,
		stdev is the standard deviation of the background (this term absorb different noise contributions such as read noise, dark and sky)
		S = flux.
		So, S/N = 1.087/merr

		Please refer to /Users/chenping/Documents/instrument_data_analysis/SNR_CCD_optical.pdf for Signal-to-Noise in Optical Astronomy

		INPUTS:
			image_key:
			photometry_function:
			which_dir:
			fixed_app_size:
		'''
		if photometry_function == 'iraf':
			if x is None or y is None:
				x = self.photometry_info[image_key]['x']
				y = self.photometry_info[image_key]['y']
			xy = np.array([x,y]).reshape((1,2))
			image = self.__get_internal_image(image_key, which_dir=which_dir)
			self.get_apphot_iraf_parameters(image_key)
			options = self.apphot_iraf_options.copy()
			output = os.path.join(self.aperture_photometry_dir, '%s_%s_%s_singlepoint.mag'%(image_key.split('.')[0], str(x), str(y)))
			photret_singlept = self.aperture_photometry_apphot_iraf_single_image_xys_given(image,xy,options, output=output) #XCENTER, YCENTER, MAG, MERR
			print photret_singlept
			SNR = 1.087/photret_singlept[0,3]
		else:
			raise ValueError('sorry, not implemented yet')

		return SNR

	def get_photometry_uncertainty_flt(self,flt, phot_method = 'apphot', mode='relmag', updatetable=1):
		'''
		get photometry unertainties for source of different magnitudes

		This is used to estimate systematic uncertainties by studying the scatter of magnitudes from all images for each star

		How? Register all stars in the template image and then find all matches in other images for each star. Then study the scatter of stars with different magnitudes

		'''
		self.__find_template_imagekey(updatetable=updatetable)
		tpl_image = self.templates[flt]
		tpl_imgkey = tpl_image.split('.')[0]
		if updatetable:
               		self.__dict2table()
		photret_table_all = self.photometry_record_table

		photret_flt = self.__select_rows_from_table(photret_table_all,'flt',flt)
		photret_flt_nondrop = self.__select_rows_from_table(photret_flt,'drop',0)

		if phot_method == 'apphot':
			magfile = os.path.join(self.aperture_photometry_dir, tpl_imgkey+self.apphot_ret_file_suffix)
		elif phot_method == 'psfphot':
			magfile = self.__get_internal_psfphot_mag_file(tpl_image)
		else:
			raise IOError("Invalid input for phot_method...")
		photret_dir = os.path.dirname(magfile)

		xys_dict = {}
		mags_dict = {}

		offset_tpl = self.photometry_info[tpl_image][mode] - self.photometry_info[tpl_image]['instmag']
		tpl_data = np.loadtxt(magfile)

		print tpl_data

		for i,phot_single in enumerate(tpl_data):
			mags_dict[i] = [phot_single[2]+offset_tpl]
			xys_dict[i] = [phot_single[0],phot_single[1]]

		for image in photret_flt_nondrop['name']:
			if image == tpl_image:
				continue

			if mode != 'calmag' and mode !='relmag':
				raise IOError("Invlalid input for mode")

			offset_mag = self.photometry_info[image][mode] - self.photometry_info[image]['instmag']

			imgkey = image.split('.')[0]
			this_tpl_match_file = os.path.join(photret_dir, imgkey+'_tpl.match')
			this_phot_data = np.loadtxt(this_tpl_match_file)

			tpl_xys_match = this_phot_data[:,[0,1]]

			for dict_key in xys_dict.keys():
				xy = xys_dict[dict_key]
				exist_flag,xy_match,match_index = self.__find_corresponding(tpl_xys_match,xy,0.5)
				if not exist_flag:
					continue

				mag = this_phot_data[match_index,6] + offset_mag
				mags_dict[dict_key].append(mag)

		mags = []
		stds  = []
		for dict_key in mags_dict.keys():
			mags_this_star = mags_dict[dict_key]
			if len(mags_this_star)<4:
				continue
			mags.append(np.mean(mags_this_star))
			stds.append(np.std(mags_this_star))

		mags_array = np.array(mags).reshape(len(mags),1)
		stds_array = np.array(stds).reshape(len(stds),1)
		out = np.hstack((mags_array,stds_array))
		outsort = out[np.argsort(out[:,0]),:]

		out_file_name = os.path.join(self.result_dir,'mags_stds_'+flt+'.txt')
		self.__delete_file_if_exist(out_file_name)

		np.savetxt(out_file_name,outsort,fmt="%5.2f %5.2f")

		plt.plot(mags,stds,'o')
		plt.xlabel('mag')
		plt.ylabel('mag scatter')
		plt.show()

	def get_light_curves_ref_stars(self, flt, mag_min_ref=None, mag_max_ref=None, phot_method='apphot', lctype='calmag', only_caldata=True, verbose=0):
		'''
		extract light curves for all reference stars which meet the given criteria


		INPUTS:
			flt:
			mag_min_ref: faintest end
			mag_max_ref: brightest end
			phot_method: 'apphot' or 'psfphot'
			lctype: 'relmag' or 'calmag'
		'''

		if flt not in self.templates.keys():
			self.__find_template_imagekey(updatetable=updatetable)

		tpl_image = self.templates[flt]
		if verbose:
			print "the template image: %s %s"%(tpl_image, self.images[tpl_image])
		tpl_imgkey = tpl_image.split('.')[0]

		if phot_method == 'apphot':
			tpl_photfile = os.path.join(self.aperture_photometry_dir, tpl_imgkey+self.apphot_ret_file_suffix)
		elif phot_method == 'psfphot':
			tpl_photfile = self.__get_internal_psfphot_mag_file(tpl_image)
		else:
			raise IOError("Invalid input for phot_method...")
		photret_dir = os.path.dirname(tpl_photfile)
		if verbose:
			print "The photometry file of tempalte image: %s"%tpl_photfile

		tpl_photdata = np.loadtxt(tpl_photfile)
		if verbose:
			print "photometry data of tempalte image:\n"
			print tpl_photdata

		lcs_ref_dir = os.path.join(self.warehouse_dir, 'lcs_ref')
		if not os.path.exists(lcs_ref_dir):
			os.mkdir(lcs_ref_dir)

		lcs_ref_spm_dir = os.path.join(lcs_ref_dir, phot_method)
		if not os.path.exists(lcs_ref_spm_dir):
			os.mkdir(lcs_ref_spm_dir)

		for photline in tpl_photdata:
			x,y,mag,magerr=photline

			if mag_min_ref is not None:
				if mag > mag_min_ref:
					continue

			if mag_max_ref is not None:
				if mag < mag_max_ref:
					continue

			reflcfile =os.path.join( lcs_ref_spm_dir,  'ref_'+flt+'_' + str(np.round(x,1)).replace('.','p') + '_' + str(np.round(y,1)).replace('.', 'p') + '.txt')
			lcdata = self.get_light_curve_for_control_star(flt, x, y, phot_method=phot_method, lctype=lctype, only_caldata=only_caldata)

			calstar_index = np.where(lcdata['mag']!=99.99)[0]
			if len(calstar_index)>1:
				lcdata.write(reflcfile, format='ascii.commented_header')
			else:
				print "the star at image coordinate (%s,%s) only exist on template image"%(x,y)

	def get_light_curve_for_control_star(self,flt,tpl_x,tpl_y, phot_method='apphot', lctype='calmag', only_caldata=True,  verbose=0, updatetable=1):
		'''
		extract light curve for a given star with position on reference image (tpl_x, tpl_y)

		INPUTS:
			flt:
			tpl_x:
			tpl_y:
			phot_method: 'apphot' or 'psfphot'
			lctype: 'calmag' or 'relmag'
			only_caldata: only extract control which have been used in the relative calibration
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		tpl_image = self.templates[flt]
		tpl_imgkey = tpl_image.split('.')[0]
		if updatetable:
                	self.__dict2table()
		photret_table_all = self.photometry_record_table
		photret_flt = self.__select_rows_from_table(photret_table_all,'flt',flt)
		photret_flt_nondrop = self.__select_rows_from_table(photret_flt,'drop',0)

		if phot_method == 'apphot':
			tplphotfile = os.path.join(self.aperture_photometry_dir, tpl_imgkey+self.apphot_ret_file_suffix)
		elif phot_method == 'psfphot':
			tplphotfile = self.__get_internal_psfphot_mag_file(tpl_image)
		else:
			raise IOError("Invalid input for phot_method...")
		photret_dir = os.path.dirname(tplphotfile)

		tpl_data = np.loadtxt(tplphotfile)
		if verbose:
			print tpl_data

		mags = []
		magerrs = []
		obstimes = []
		for i,image in enumerate(photret_flt_nondrop['name']):
			obstimes.append(photret_flt_nondrop['obstime'][i])
			offset_mag = self.photometry_info[image][lctype] - self.photometry_info[image]['instmag']

			imgkey = image.split('.')[0]
			if imgkey == tpl_imgkey:
				if phot_method == 'apphot':
					tpl_phot_suffix = '.apphot'
				else:
					tpl_phot_suffix = '.psfphot'
				this_tpl_match_file = os.path.join(photret_dir,imgkey+tpl_phot_suffix)
			else:
				if only_caldata:
					this_tpl_match_file = os.path.join(photret_dir, imgkey+'_tpl_caldata.match')
				else:
					this_tpl_match_file = os.path.join(photret_dir, imgkey+'_tpl.match')

			this_phot_data = np.loadtxt(this_tpl_match_file)
			tpl_xys_match = this_phot_data[:,[0,1]]

			xy = [tpl_x,tpl_y]
			exist_flag,xy_match,match_index = self.__find_corresponding(tpl_xys_match,xy,1.5)
			if not exist_flag:
				mags.append(99.99)
				magerrs.append(99.99)
				continue
			if imgkey == tpl_imgkey:
				mag = this_phot_data[match_index,2] + offset_mag
				magerr_inst = this_phot_data[match_index,3]
			else:
				mag = this_phot_data[match_index,6] + offset_mag
				magerr_inst = this_phot_data[match_index,7]

			magerr_relcal = self.photometry_info[image]['relmagerr']
			magerr = np.sqrt(magerr_inst**2 + magerr_relcal**2)


			mags.append(mag[0])
			magerrs.append(magerr[0])

		ref_star_data  = [obstimes, mags, magerrs]
		control_star_table = Table(ref_star_data, names = ['obstime', 'mag', 'magerr'])

		#control_star_table = Table()
		#obstime_col = Column(data=obstimes,name='obstime')
		#mag_col     = Column(data=mags,name='mag')
		#magerr_col  = Column(data=magerrs,name='magerr')
		#control_star_table.add_columns([obstime_col,mag_col,magerr_col])

		return control_star_table

#global filtering
	def check_battle_damage(self,flts='all', updatetable=1):
		'''
		Check images which have been droped with non-zero drop flag
		'''
		if updatetable:
			self.__dict2table()

		if flts == 'all':
			self.__find_template_imagekey(updatetable=updatetable)
			flts_tocheck = self.templates.keys()

			if len(flts_tocheck) == 0:
				flts_tocheck = np.unique(self.photometry_record_table['flt'])

		elif isinstance(flts,str):
			flts_tocheck = [flts]
		elif isinstance(flts,list):
			flts_tocheck = flts
		else:
			raise IOError('Invalid input for flts...')


		droped_images = []
		for flt in flts_tocheck:

			print "Check %s band image:"%flt

			info_flt = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
			for img in info_flt['name']:
				drop_status =  self.photometry_info[img]['drop']
				if drop_status != 0:
					print img,drop_status
					droped_images.append(img)
		return droped_images

	def block_specific_image(self,colname,criteria,mode='eq', updatetable=1):
		'''
		block images which you don't want to process

		if mode == 'eq', image with specific value == criteria will be assigned 10 to drop flag
		if mode == 'gt', image with specific value > criteria will be assigned 10 to drop flag
		if mode == 'lt', image with specific value < criteria will be assigned 10 to drop flag
		'''
		if updatetable:
			self.__dict2table()
		colnames = self.photometry_record_table.colnames

		if colname not in colnames:
			raise IOError("%s not exist in table"%colname)

		for i,value in enumerate(self.photometry_record_table[colname]):
			imgkey = self.photometry_record_table['name'][i]
			print imgkey
			if mode == 'eq':
				if value == criteria:
					self.photometry_info[imgkey]['drop'] = 10
			elif mode == 'gt':
				if value > criteria:
					self.photometry_info[imgkey]['drop'] = 10
			elif mode == 'lt':
				if value < criteria:
					self.photometry_info[imgkey]['drop'] = 10
			else:
				print "what"

	def bring_dead_back_to_life(self,drop_flag='all', updatetable=1):
		'''
		renew drop flag for images with drop signal as 'drop_flag'

		Inputs:
			drop_flag: default value is 'all' which bring all dead to life
				   You can rescue specific images by specify drop_flags,
				   for example drop_flag = '3'
				   Please see self.drop_status for the details on drop flags

		'''
		drop_flags_allowed = self.drop_status.keys()
		if updatetable:
			self.__dict2table()
		imginfo_table = self.photometry_record_table

		if drop_flag == 'all':
			save_images = imginfo_table['name']
		elif str(drop_flag) in drop_flags_allowed:
			info_table_specific_drop = self.__select_rows_from_table(imginfo_table,'drop',drop_flag)
			N_save = len(info_table_specific_drop)
			save_images = info_table_specific_drop['name']
		else:
			raise IOError('Sorry, invalid input for drop_flag')

		if drop_flag != 'all' and N_save == 1:
			self.photometry_info[save_images]['drop'] = 0
		else:
			for img in save_images:
				self.photometry_info[img]['drop'] = 0

	def check_images_given_criteria(self, flt=None, obstimemin=None, obstimemax=None, selectcolname=None, colmatch=None, matchtype='eq', skipdroped=0, updatetable=0, image_dir='raw_image'):
		'''
		filter the images meeting given criteria, load them to ds9 and check their quality to designate drop flag.
		'''
		if updatetable:
			self.__dict2table()
		photret_table = self.photometry_record_table.copy()
		if skipdroped:
			photret_table = photret_table[photret_table['drop']==0]

		if flt is not None:
			photret_table = photret_table[photret_table['flt']==flt]
		
		if obstimemin is not None:
			photret_table = photret_table[photret_table['obstime']>obstimemin]
		if obstimemax is not None:
			photret_table = photret_table[photret_table['obstime']<obstimemax]

		if selectcolname is not None and colmatch is not None:
			if matchtype == 'eq':
				photret_table = photret_table[photret_table[selectcolname]==colmatch]
			elif matchtype == 'gt':
				photret_table = photret_table[photret_table[selectcolname]>colmatch]
			elif matchtype == 'lt':
				photret_table = photret_table[photret_table[selectcolname]<colmatch]
			else:
				raise ValueError('matchtype %s not supported'%matchtype)

		for imgkey in photret_table['name']:
			self.__load_check_and_decidedrop_single_image(imgkey, image_dir=image_dir)

	def check_images_picked_from_lc_flt(self,flt,lctype = 'calmag', image_dir = 'raw_image', updatetable=1):
		'''
		display the image corresponding to the point (obstime,mag) in light curve

		lctype:	'calmag', 'relmag', 'instmag'
		'''
		if updatetable:
			self.__dict2table()
		photret_table_all = self.photometry_record_table

		photret_table_nondrop = self.__select_rows_from_table(photret_table_all,'drop',0)
		flt_photret = self.__select_rows_from_table(photret_table_nondrop, 'flt', flt)
		obstimes = flt_photret['obstime']

		if lctype not in ['calmag','relmag','instmag']:
			raise IOError("Invalid input for lctype, calmag, relmag or instmag are allowed")

		mags = flt_photret[lctype]
		magerrs = flt_photret[lctype+'err']
		obstimes_selected,mags_selected,index_selected=mouse_pick(obstimes,mags,err=magerrs) #mouse_pick can pick up multiple points

		for obstime_wanted,mag_wanted in zip(obstimes_selected, mags_selected):
			table_selected_from_obstime = self.__select_rows_from_table(flt_photret,'obstime',obstime_wanted,mode='lgt',criteria=0.001)
			table_selected_from_mag = self.__select_rows_from_table(flt_photret, lctype, mag_wanted,mode='lgt',criteria=0.01)

			for name in table_selected_from_obstime['name']:
				if name in table_selected_from_mag['name']:
					final_selected_imgkey = name
			self.__load_check_and_decidedrop_single_image(final_selected_imgkey, image_dir=image_dir)

#matplotlib visualization
	def show_source_stamps(self, imagedata, sources, xcol, ycol, halfwidth=10, baseindex=1):
		'''
		IRAF, ds9 both index starting from 1
		'''
		from matplotlib import rcParams
		rcParams['figure.subplot.left'] = 0.03
		rcParams['figure.subplot.right'] = 0.97
		rcParams['figure.subplot.bottom'] = 0.03
		rcParams['figure.subplot.top'] = 0.97
		
		fig = plt.figure(figsize=(15,9))
		numrow = int(len(sources)/15+1)
		gs = gridspec.GridSpec(numrow, 15, wspace=0.0, hspace=0)
		for i,(x,y) in enumerate(sources[[xcol, ycol]]):
			indrow = int(i/15)
			indcol = np.mod(i,15)
			ax = plt.subplot(gs[indrow, indcol])
			#ax = fig.add_subplot(numrow, 10, i+1)			
			xint = int(x) - baseindex
			yint = int(y) - baseindex
			stampdata = imagedata[(yint-halfwidth):(yint+halfwidth), (xint-halfwidth):(xint+halfwidth)]
			ax.imshow(stampdata, origin='lower', cmap='gray_r', interpolation='none', vmin=np.min(stampdata), vmax=np.max(stampdata)*0.1)
			nx, ny = stampdata.shape
			linecut1 = stampdata[halfwidth,:]
			ax.plot(np.arange(ny), linecut1/np.max(linecut1)*0.5*halfwidth+halfwidth)
			linecut2 = stampdata[:,halfwidth]
			ax.plot( linecut2/np.max(linecut2)*0.5*halfwidth+halfwidth, np.arange(nx))
			ax.text(1.0*halfwidth, 1.6*halfwidth, str(i), fontsize=10)
			ax.text(0, 0, '(%s,%s)'%(xint, yint))
			ax.set_xlim([-1, 2*halfwidth+2])
			ax.set_ylim([-1, 2*halfwidth+2])
			ax.set_xticks([])
			ax.set_yticks([])
		plt.show(block=False)

	def xy_plot_photometry_table(self, xcol, ycol, flt=None, xerrcol=None, yerrcol=None, skipdroped=1, updatetable=1):
		'''
		X-Y plot of two given columns from photometry result table
		'''
		if updatetable:
			self.__dict2table()
		if flt is not None:
			table_flt =self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
		else:
			table_flt = self.photometry_record_table
		if skipdroped:
			table_flt = self.__select_rows_from_table(table_flt,'drop',0)
		xy_plot_general(table_flt, xcol, ycol,xerrcol=xerrcol, yerrcol=yerrcol)

	def lc_plot(self, magcolname,flts = 'all',bining =False, bin_width = 0.5, t0=0,remove_inbin_outliers=False, plot_droped = False,plot_bined_only = False, updatetable=1, extra_filter_table=None, extra_filter_colname=None, extra_filter_colvalue=None, xlabeltext =None, ylabeltext=None):
		'''
		Inputs:
			magcol: instmag; relmag; calmag; if image_subtraction == 1, there are more
			flts: 'all' or 'V' or ['B','V']
			bining: bin the data of not if multiple measurements obtained within the one day
			bin_width: the widt for the bin
			extra_filter_table, extra_filter_colname, extra_filter_colvalue: if not None, first filtering the self.photometry_record_table by extra_filter_table[extra_filter_colname]==extra_filter_colvalue
		'''
		if updatetable:
			self.__dict2table()
		plttable = self.photometry_record_table.copy()

		if (extra_filter_table is not None) and (extra_filter_colname is not None) and (extra_filter_colvalue is not None):
			plttable = plttable[extra_filter_table[extra_filter_colname]==extra_filter_colvalue]

		self.__find_template_imagekey(updatetable=updatetable)

		magerrcolname = magcolname.replace('mag', 'magerr')

		if magcolname not in plttable.colnames:
			raise ValueError("wrong input for magcol")

		if flts == 'all':
			#flts_plot = self.templates.keys()
			flts_plot = np.unique(plttable['flt'])
		else:
			if isinstance(flts,str):
				flts  = [flts]
			flts_plot = flts

		from matplotlib import rcParams
		rcParams['figure.subplot.left'] = 0.13
		rcParams['figure.subplot.right'] = 0.95
		rcParams['figure.subplot.bottom'] = 0.13
		rcParams['figure.subplot.top'] = 0.95
		fig, ax = plt.subplots(figsize=(9,6))

		for flt in flts_plot:
			print flt
			table_flt =self.__select_rows_from_table(plttable,'flt',flt)
			if not plot_droped:
				table_flt = self.__select_rows_from_table(table_flt,'drop',0)

			if len(table_flt) == 0:
				continue

			table_flt.sort('obstime')
			obstimes = table_flt['obstime']
			mags = table_flt[magcolname]
			magerrs = table_flt[magerrcolname]

			if not plot_bined_only:
				ax.errorbar(obstimes-t0,mags,yerr=magerrs,fmt='o', label=flt)

			if bining:
				obstimes_bined = []
				mags_bined = []
				magerrs_bined = []

				bined_already_index = []
				for i,obstime in enumerate(obstimes):
					if i in bined_already_index:
						continue

					distances = obstimes - obstime
					indexs = np.where(np.abs(distances)< bin_width)[0]
					print indexs
					bined_already_index = np.append(bined_already_index,indexs)

					obstime_bined = np.mean(obstimes[indexs])
					mags_tobin = mags[indexs]
					magerrs_tobin = magerrs[indexs]

					#remove outlier measurements
					if remove_inbin_outliers and len(mags_tobin)>1:
						mags_tobin_nonoutlier,index_keep,stdev = sigma_clipping(mags_tobin,sig=2)
						magerrs_tobin_nonoutlier = magerrs_tobin[index_keep]
					else:
						mags_tobin_nonoutlier = mags_tobin
						magerrs_tobin_nonoutlier = magerrs_tobin

					print mags_tobin_nonoutlier

					weights = 1./magerrs_tobin_nonoutlier**2
					mag_bined = np.sum(weights*mags_tobin_nonoutlier) / np.sum(weights)
					magerr_bined = 1./ np.sqrt(np.sum(weights))

					obstimes_bined.append(obstime_bined)
					mags_bined.append(mag_bined)
					magerrs_bined.append(magerr_bined)
				ax.errorbar(obstimes_bined-t0,mags_bined, yerr=magerrs_bined,fmt='s', label=flt+'_bin')
		if ylabeltext is None:
			ylabeltext = 'Magnitude'
		ax.set_ylabel(ylabeltext, fontsize=18)
		if xlabeltext is None:
			if t0==0:
				xlabeltext = 'JD'
			else:
				xlabeltext = 'JD - %s'%str(t0)
		ax.set_xlabel(xlabeltext, fontsize=18)
		ax.invert_yaxis()
		ax.legend(loc=0)
		plt.show()

#Photometry records management
	def shift_obstime_from_exposure_start_to_exposure_middle(self):
		'''
		when accurate obstimes are required, we need to make sure the reported obstime is at the middle of the exposure
		'''

		for img in self.images.keys():
			self.photometry_info[img]['obstime'] = self.photometry_info[img]['obstime'] + self.photometry_info[img]['exptime']/2.0/3600/24

	def fill_phot_info_specific_col(self,colname,flts='all', updatetable=1, fill_value=None, skiplist=None):
		'''
		Fill in the photometry record dictionary with fill_value or the default value in self.photometry_info_init_values

		Inputs:
			colname: check self.photometry_info_keys for valid inputs
			flts:	default 'all',
				non-default input should be a list or np.ndarray containing the
				ilters, eg flts= ['B'] or flts = ['B','V']
			updatetable:
			fill_value: if not None, will override the default value
		'''
		if updatetable:
			self.__dict2table()
		imginfo_table = self.photometry_record_table

		if flts == 'all':
			flts = np.unique(imginfo_table['flt'])

		if isinstance(flts,str):
			flts = [flts]

		for flt in flts:
			init_images = self.__select_rows_from_table(imginfo_table,'flt',flt)
			for img in init_images['name']:
				if skiplist is not None:
					if img in skiplist:
						continue
				if fill_value is None:
					self.photometry_info[img][colname] = self.photometry_info_init_values[colname]
				else:
					self.photometry_info[img][colname] = fill_value

	def update_photometry_record_table(self, colnames):
		if isinstance(colnames, str):
			colnames = [colnames]
		for i,imgkey in enumerate(self.photometry_record_table['name']):
			for colname in colnames:
				self.photometry_record_table[colname][i] = self.photometry_info[imgkey][colname]
			

	def __dict2table(self, sortcol='obstime'):
		'''
		self.photometry_info --> self.photometry_record_table
		'''
		self.photometry_record_table = Table(names=self.photometry_info_keys,dtype=self.photometry_info_dtypes)
		#print len(self.photometry_record_table.colnames)
		#print self.photometry_record_table.colnames
		for img in self.photometry_info.keys():
			#print len(self.photometry_info[img].values())
			#print self.photometry_info[img]
			self.photometry_record_table.add_row(self.photometry_info[img].values())
		for colname in ['relmag', 'relmagerr','calmag','calmagerr']:
			self.photometry_record_table[colname] = np.round(self.photometry_record_table[colname], 3)
		self.photometry_record_table.sort(sortcol)

	def __load_old_results(self,info_table,info_dict):
		'''
		load info into info_dict from info_table
		'''
		for key in info_dict.keys():
			for colname in info_table.colnames:
				if colname not in self.photometry_info_keys:
					continue
				name_key = np.where(info_table['name'] == key)[0]
				if len(name_key):
					name_key_indice = name_key[0]
					info_dict[key][colname] = info_table[name_key_indice][colname]

	def get_phot_zpts(self):
		'''
		Get photometric zero points
		'''
		for imgkey in self.images.keys():
			if self.photometry_info[imgkey]['instmag'] == 99.99 or self.photometry_info[imgkey]['calmag'] == 99.99:
				continue
			zpt, zpterr = self.get_phot_zpt_single_image(imgkey)
			self.photometry_info[imgkey]['magzpt'] = zpt
			self.photometry_info[imgkey]['magzpterr'] = zpterr

	def get_phot_zpt_single_image(self, imgkey):
		'''
		Get the photometric zero point of the image when the instru and calibrated magnitudes are avaiable.
		For self.photometry_method == 'apphot', instmag = -2.5*log10(fcount)+25
		For self.photometry_method == 'psfphot' (DoPhot), instmag = -2.5*log10(fcount) where fcount is total count

		The zero point is time averaged.
		'''
		instmag = self.photometry_info[imgkey]['instmag']
		instmagerr = self.photometry_info[imgkey]['instmagerr']
		calmag = self.photometry_info[imgkey]['calmag']
		calmagerr  = self.photometry_info[imgkey]['calmagerr']
		exptime = self.photometry_info[imgkey]['exptime']

		zpt = calmag - instmag - 2.5*np.log10(exptime)
		if self.photometry_method == 'apphot':
			zpt = zpt + 25.0
		zpterr = np.sqrt(calmagerr**2-instmagerr**2)

		return zpt, zpterr

#Image transformation
	def __transform_xys(self, dxfit_values, dyfit_values, porder, input_xys):
		'''
		allowd input_xys: (x,y) or (N,2) array
		dxfit_values = [a0,a1,a2]
		dyfit_values = [b0,b1,b2]
		(x,y) ?==> (x',y'): x'= a0+a1*x+a2*y; y'=b0+b1*x+b2*y
		'''
		def poly_trans(porder,dxfit,dyfit,x,y):
			if porder == 1:
				xout = dxfit[0] + dxfit[1]*x + dxfit[2]*y
				yout = dyfit[0] + dyfit[1]*x + dyfit[2]*y
			elif porder == 2:
				#!!! the following formula not checked
				xout = dxfit[0] + dxfit[1]*x + dxfit[2]*y  + dxfit[3]*x**2 + dxfit[4]*x*y + dxfit[5]*y**2
				yout = dyfit[0] + dyfit[1]*x + dyfit[2]*y  + dyfit[3]*x**2 + dyfit[4]*x*y + dyfit[5]*y**2
			elif porder == 3:
				#!!! the following formula not checked
				xout = dxfit[0] + dxfit[1]*x + dxfit[2]*y  + dxfit[3]*x**2 + dxfit[4]*x*y + dxfit[5]*y**2 + dxfit[6]*x**3 + dxfit[7]*x**2*y + dxfit[8]*x*y**2 + dxfit[9]*y**3
				yout = dyfit[0] + dyfit[1]*x + dyfit[2]*y  + dyfit[3]*x**2 + dyfit[4]*x*y + dyfit[5]*y**2 + dyfit[6]*x**3 + dyfit[7]*x**2*y + dyfit[8]*x*y**2 + dyfit[9]*y**3
			else:
				raise ValueError('poly order higher than 3 not supported yet')
			return xout,yout


		if isinstance(input_xys,list):
			x = input_xys[0]
			y = input_xys[1]
			xout,yout = poly_trans(porder, dxfit_values, dyfit_values, x,y)
			output_xys = [xout,yout]
		else:
			output_xys = None
			for xy in input_xys:
				x = xy[0]
				y = xy[1]
				xout,yout = poly_trans(porder, dxfit_values, dyfit_values, x,y)
				output_xy = np.array([xout,yout]).reshape((1,2))
				if output_xys is None:
					output_xys = output_xy
				else:
					output_xys = np.append(output_xys,output_xy,axis=0)

		return output_xys

	def __extract_transformation_fitting(self,trans_file):
		'''
		Read in data from transformation coefficient file
		'''
		fid = open(trans_file,'r')
		data= fid.readlines()
		order = int(data[1].strip().split('=')[1])
		dxfit = data[-2].strip()
		dyfit = data[-1].strip()
		dxfit_strs = dxfit.split('=')[1].split(',')
		dyfit_strs = dyfit.split('=')[1].split(',')
		dxfit_values = [float(xv) for xv in dxfit_strs]
		dyfit_values = [float(yv) for yv in dyfit_strs]

		return order, dxfit_values,dyfit_values

	def __fitsh_grtrans(self,input_list,output_list,input_trans,):
		'''
		For details on grtrans, please refer to http://fitsh.szofi.net/task/grtrans

		grtransh [options] <input> [-o <output>]
		'''
		self.__delete_file_if_exist(output_list)
		grtrans = os.path.join(self.base.fitsh_dir, 'grtrans')
		command = "%s -i %s -o %s --input-transformation %s --col-xy 1,2"%(grtrans,input_list,output_list,input_trans)
		try:
			os.system(command)
		except:
			print "grtrans failure..."

	def __fitsh_fitrans(self, input_image, output_image, trans_file, sx=None, sy=None):
		'''
		perform geometric transformatios on the input_image according to transformation in trans_file
		For details on fitrans refer to http://fitsh.szofi.net/task/fitrans

		INPUTS:
			input_image:
			output_image:
			trans_file:
			sx, sy: the size of the output_image if it should differ from the original image size
		'''
		self.__delete_file_if_exist(output_image)
		if not os.path.exists(trans_file):
			raise ValueError('transformation file %s not exist'%trans_file)
		fitrans = os.path.join(self.base.fitsh_dir, 'fitrans')
		command = "%s -i %s -o %s -T %s"%(fitrans, input_image, output_image, trans_file)
		if sx is not None and sy is not None:
			command = command + ' -s %s,%s'%(sx, sy)
		for kw in self.fitsh_fitrans_pars.keys():
			if self.fitsh_fitrans_pars[kw]:
				command = command + ' ' + kw
		try:
			os.system(command)
		except:
			print "fitrans failure"

	def __fitsh_grmatch(self,ref_list,input_list, match_output, trans_output=None, mode ='file'):
		'''
		For details on grmatch refer to http://fitsh.szofi.net/task/grmatch

		Basics for grmatch:
		grmatch [options] -r <reference> -i <input> [-o <output>]

		The program 'grmatch' matches lines read from two input files,
		namely from a reference and from an input file.
		All implemented algorithms are symmetric, in the manner that
		the result should be the same if these two files
		are swapped. The only case when the order of these files is
		important is when a geometrical transformation is
		also returned (see point matching below), in this case
		the swapping of the files results the inverse form of the
		original transformation. The lines (rows) can be matched using various criteria.


		If there is no negative sign before the column index, the data are sorted
		in descending(!) order, therefore the lines with the lines with the highest(!)
		values are selected for triangulation. If there is a negative sign before the index,
		the data are sorted in ascending order by these values, therefore the lines with
		the smallest(!) values are selected for triangulation.

		Inputs:
			ref_list: filename or data depending on 'mode'
			input_list:
			match_output: the table file containing match sources
			trans_output: the transformation coefficient file
			mode: 'file' or 'data'
		'''
		self.__delete_file_if_exist(match_output)
		if trans_output is not None:
			self.__delete_file_if_exist(trans_output)

		if mode == 'file':
			ref_list_filename = ref_list
			input_list_filename = input_list
		elif mode == 'data':
			ref_list_filename = "grmatch_ref_list_temp.txt"
			input_list_filename = "grmatch_input_list_temp.txt"
			np.savetxt(ref_list_filename,ref_list)
			np.savetxt(input_list_filename, input_list)

		grmatch = os.path.join(self.base.fitsh_dir, 'grmatch')
		command = "%s -r %s -i %s -o %s --output-transformation %s"%(grmatch,ref_list_filename,input_list_filename,match_output,trans_output)
		if self.fitsh_grmatch_type == 'point':
			command = command + ' --match-points'
			for kw in self.fitsh_grmatch_pointmatch_pars.keys():
				if self.fitsh_grmatch_pointmatch_pars[kw] is not None:
					command = command + ' ' + kw + ' ' + str(self.fitsh_grmatch_pointmatch_pars[kw])
		elif self.fitsh_grmatch_type == 'coord':
			command = command + ' --match-coord'
			for kw in self.fitsh_grmatch_coordmatch_pars.keys():
				if self.fitsh_grmatch_coordmatch_pars[kw] is not None:
					command = command + ' ' + kw + ' ' + str(self.fitsh_grmatch_coordmatch_pars[kw])
		elif self.fitsh_grmatch_type == 'id':
			command = command + ' --match-id'
			for kw in self.fitsh_grmatch_idmatch_pars.keys():
				if self.fitsh_grmatch_idmatch_pars[kw] is not None:
					command = command + ' ' + kw + ' ' + str(self.fitsh_grmatch_idmatch_pars[kw])
		else:
			raise ValueError('invalid input for self.fitsh_grmatch_type')

		try:
			os.system(command)
		except:
			print "grmatch failure"

		if mode == 'data':
			os.remove(ref_list_filename)
			os.remove(input_list_filename)

		if os.path.exists(match_output):
			matched_ret = np.loadtxt(match_output)
		else:
			matched_ret = None

		return matched_ret

	def __get_trans_matrix_diapl_cross(self,list1,list2):

		print "on the way"

	def __get_rough_cross_shift_by_imaga_correlation(self,ref_img,target_img):
		'''
		cross cross.par instrument.par output_file reference_image target_image

		the function give the relative shift of target image against the reference image
		suppose sources on refence image (x,y); sources on the target image (x',y')
		in the ideal circumstance, x = x'+dx; y = y'+dy

		'''

		#prepare the images(cut the image to 2048*2048 for cross correlation)
		in_img_name = ref_img
		cutted_ref_name = 'cutted_ref_img.fits'
		self.__cut_image_diapl(in_img_name,cutted_ref_name,1,1,2048,2048)

		in_img_name = target_img
		cutted_target_name = 'cutted_target_img.fits'
		self.__cut_image_diapl(in_img_name,cutted_target_name,1,1,2048,2048)

		shifts_file = "temp_shift.dat"
		self.__delete_file_if_exist(shifts_file)

		#cross correlate two images to get the shift amout;
		cross_par = os.path.join(self.parafile_dir,'cross.par')
		instrument_par = os.path.join(self.parafile_dir,'instrument.par')

		cross = os.path.join(self.base.diapl_dir, 'cross')
		cross_command = "%s %s %s %s %s %s >/dev/null" %(cross, cross_par,instrument_par,shifts_file, cutted_ref_name, cutted_target_name)
		if os.system(cross_command):
		    print "Error: sorry dude, cannot run", command
		    sys.exit(1)


		os.remove(cutted_ref_name)
		os.remove(cutted_target_name)

		fid = open(shifts_file)
		line = fid.read()
		prepared_target_image, dx, dy = tuple(line.split())

		os.remove(shifts_file)

		return dx,dy

#coordinate transformation
	def __world2image_fitshead(self, hdr, worldcoor):
		w = WCS(hdr)
		image_coor = w.wcs_world2pix(worldcoor,1) # Pixel coordinates of (RA, DEC)
		return image_coor

	def __world2image(self,image,world, hduindex=None):
		hdu  = fits.open(image)
		if hduindex is None:
			header = hdu['PRIMARY'].header
		else:
			header = hdu[hduindex].header

		image_coor = self.__world2image_fitshead(header, world)

		return image_coor

	def __image2world_fitshead(self, hdr, imagecoor):
		'''
		convert from image coordinate to physical world coordinate

		'''
		w = WCS(hdr)
		world_coor = w.wcs_pix2world(imagecoor,1)

		return world_coor

	def __image2world(self,image,image_xys, hduindex=None):
		hdu = fits.open(image)

		if hduindex is None:
			header = hdu['PRIMARY'].header
		else:
			header = hdu[hduindex].header

		world_coor = self.__image2world_fitshead(header, image_xys)

		return world_coor

#Standard star catalog
	def get_standards(self,catalog = 'apass', center_ra = None, center_dec = None, distance = None):
		'''
		search the standard star catalog to extract subset of the data for future match with stars in working image field
		currently supported catalogs: apass and 2mass

		INPUTS:
			catalog: 'apass', '2mass' or 'panstarrs', 'refcat2','sdss'
			center_ra: catalog search center ra; if None, self.sn_ra_world_deg will be used
			center_dec: catalog search center dec; if None, self.sn_dec_world_deg will be used
			distance: search radius in degree; if None, self.std_region_radius will be used
		'''

		std_ref_stars_file = os.path.join(self.std_ref_dir, self.std_ref_stars_file[catalog])
		renew_stds = self.renew_standards
		if os.path.exists(std_ref_stars_file) and renew_stds:
			os.remove(std_ref_stars_file)

		if not os.path.exists(std_ref_stars_file):
			if center_ra is None or center_dec is None:
				center_ra, center_dec = self.__get_standards_prepare()
			if catalog == 'panstarrs' and center_dec<-30:
				print "====> ATTENTION!!! you are asking for PS1 data beyond -30 degree of the south"
				#raise ValueError('PanSTARRS catalog does not cover region beyond -30 degree for the south')
			if distance is None:
				distance = self.std_region_radius
			self.__get_standards(catalog=catalog, center_ra = center_ra, center_dec =center_dec, distance = distance)
		else:
			self.standards = Table.read(std_ref_stars_file, format='ascii.csv')

	def __get_standards(self, catalog='apass', center_ra=None, center_dec=None, distance=None):
		'''
		standard reference data will be saved in self.std_ref_stars_file[catalog]
		and self.standards will be updated

		INPUTS:
			Refer to self.get_standards
		'''
		if center_ra is None or center_dec is None:
			raise ValueError('subset region center (center_ra, center_dec) required...')

		if distance is None:
			raise ValueError('standard region radius required...')

		std_ref_dir = self.std_ref_dir
		std_ref_stars_file = self.std_ref_stars_file[catalog]
		output_filename = os.path.join(std_ref_dir, std_ref_stars_file)
		print "The standard stars will be saved at %s"%output_filename

		if catalog == '2mass':
			twomassmethod = self.twomass_catalog_method
			if twomassmethod == 1:
				query_VO_SCS(center_ra, center_dec,distance,table='fp_psc',out_format='csv', output_file = output_filename)
				stdstemp = Table.read(output_filename, format='ascii.csv')
				for colname in self.twomass_colnames_dict.keys():
					colnamenew = self.twomass_colnames_dict[colname]
					colnameold = self.twomass_colnames_dict_method1[colname]
					stdstemp.rename_column(colnameold, colnamenew)
				stdstemp.write(output_filename, format='ascii.csv', overwrite=1)
			elif twomassmethod == 2:
				twomass_query_Vizier(center_ra, center_dec, distance, outfile=output_filename, allcolumns=self.Vizier_catalog_all_columns)
			else:
				raise IOError('2MASS catalog download method %s not available yet'%str(twomassmethod))
		elif catalog == 'apass':
			apassmethod  = self.apass_catalog_method
			if apassmethod == 1:
				stdref_database_dir  = self.local_apass_dir
				if not os.path.exists(stdref_database_dir):
					raise ValueError("the local apass data expected in %s not available"%stdref_database_dir)
				query_local_APASS(center_ra, center_dec,distance,output_file=output_filename, local_database_dir=stdref_database_dir)
			elif apassmethod == 2:
				apass_query_Vizier(center_ra, center_dec, distance, outfile=output_filename, allcolumns=self.Vizier_catalog_all_columns)
			else:
				raise IOError('PS1 catalog download method %s not available yet'%str(PS1method))
		elif catalog == 'panstarrs':
			PS1method = self.panstarrs_catalog_method
			if PS1method ==1:
				if self.std_region_radius > 0.25:
					raise ValueError('sorry, radius too large...')
				query_General_MAST(center_ra, center_dec, distance, FORMAT='csv', catalog='PS1V3OBJECTS', filename=output_filename)
			elif PS1method == 2:
				panstarrs_query_Vizier(center_ra, center_dec, distance, outfile=output_filename, allcolumns=self.Vizier_catalog_all_columns)
			else:
				raise IOError('PS1 catalog download method %s not available yet'%str(PS1method))
		elif catalog == 'refcat2':
			refcat2method = self.atlas_refcat2_method
			if refcat2method == 1:
				stdref_database_dir  = self.local_atlas_refcat2_dir
				if not os.path.exists(stdref_database_dir):
					raise ValueError("the local ATLAS refcat2 data expected in %s not available"%stdref_database_dir)
				query_local_ATLAS_Refcat2(center_ra, center_dec, distance, local_database_master_dir=stdref_database_dir, outfile=output_filename)
			else:
				raise IOError('ATLAS refcat2 catalog download method %s not available yet'%str(refcat2method))
		elif catalog == 'sdss':
			sdssmethod = self.sdss_catalog_method
			if sdssmethod == 1:
				SDSS_DR12_query_Vizier(center_ra, center_dec, distance, outfile=output_filename, allcolumns=self.Vizier_catalog_all_columns)
			else:
				raise IOError('SDSS catalog download method %s not available yet'%str(sdssmethod))
		elif catalog == 'ukidss':
			ukidssmethod = self.ukidss_catalog_method
			if ukidssmethod == 1:
				UKIDSS_query_Vizier(center_ra, center_dec, distance, outfile=output_filename, allcolumns=self.Vizier_catalog_all_columns)
			else:
				raise IOError('UKIDSS catalog download method %s not available yet'%str(ukidssmethod))
		elif catalog == 'gaia':
			gaiamethod = self.gaia_catalog_method
			if gaiamethod == 1:
				gaia2_query_Vizier(center_ra, center_dec, distance, outfile=output_filename, allcolumns=self.Vizier_catalog_all_columns)
			elif gaiamethod == 2:
				gaia_query_mast_Catalogs(center_ra, center_dec, distance, outfile=output_filename, version=2)
			else:
				raise IOError('UKIDSS catalog download method %s not available yet'%str(ukidssmethod))
		else:
			print "Warning: the catalog %s not supported yet"%catalog

		self.standards = Table.read(output_filename, format='ascii.csv')

	def __get_standards_prepare(self, verbose=1):
		'''
		get the coordinate of the target in interest
		'''
		if self.sn_ra_world_deg is not None and self.sn_dec_world_deg is not None:
			if verbose:
				print "target (RA,Dec): (%s, %s)"%(self.sn_ra_world_deg, self.sn_dec_world_deg)
		else:
			sn_ra = raw_input("Enter the RA of the current supernova:")
			sn_dec = raw_input("Enter the Dec of the current supernova:")
			self.sn_ra_world_deg = float(sn_ra)
			self.sn_dec_world_deg = float(sn_dec)

		ra = self.sn_ra_world_deg
		dec = self.sn_dec_world_deg

		return ra,dec

	def load_2mass_standards_to_ds9_image(self, input_image,text1col=None, text2col=None,newds9=True):
		'''
		load 2mass catalog in DS9 to given image
		The standard stars are loaded from catalog file standard_reference_star_2mass.txt

		INPUTS:
		'''
		twomassstdfile = os.path.join(self.std_ref_dir, self.std_ref_stars_file['2mass'])
		if not os.path.exists(twomassstdfile):
			raise ValueError('standard star catalog %s not exists... you can prepare that with self.get_standards()'%twomassstdfile)

		standards = Table.read(twomassstdfile, format='ascii.csv')
		self.__load_standards_region_to_ds9_image(input_image, standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)

	def load_panstarrs_standards_to_ds9_image(self, input_image, text1col=None, text2col=None, magcol=None, magbcut=None, magfcut=None, nbright=None, fbright=None, newds9=True):
		'''
		load PS1 catalog in DS9 to given image
		The standard stars are loaded from catalog file standard_reference_star_panstarrs.txt
		INPUTS:
			input_image:
			magcol: the magnitude column providing data for filtering
			magbcut: the bright end cut of the loaded sources; bright corresponds to small value of magnitudes
			magfcut: the faint end cut of the loaded sources
			nbright: the brightest number of targets to be loaded;
			fbright: the fraction of the brightest sources to be loaded
		'''

		ps1stdfile = os.path.join(self.std_ref_dir, self.std_ref_stars_file['panstarrs'])
		if not os.path.exists(ps1stdfile):
			raise ValueError('standard star catalog %s not exists... you can prepare that with self.get_standards()'%ps1stdfile)

		standards = Table.read(ps1stdfile, format='ascii.csv')
		if standards.masked:
			standards = standards.filled(99.99)

		if np.any(np.array([magbcut, magfcut, nbright, fbright])):
			if magcol is None:
				raise ValueError('the magnitude column required to provide data for filtering ')
			if magbcut is not None and magfcut is not None:
				if magbcut > magfcut:
					raise ValueError("bright end cut must brighter than the faint end cut")
				standards = standards[standards[magcol]>magbcut*standards[magcol]<magfcut]
			elif magbcut is not None:
				standards = standards[standards[magcol]>magbcut]
			elif magfcut is not None:
				standards = standards[standards[magcol]<magfcut]
			else:
				standards.sort(magcol)
				Nn = nbright
				if nbright is not None and fbright is not None:
					Nf = len(standards)*fbright
					if Nn > Nf:
						standards = standards[:Nf]
					else:
						standards = standards[:Nn]
				else:
					if nbright is not None:
						standards = standards[:Nn]
					if fbright is not None:
						Nf = len(standards)*fbright
						standards = standards[:Nf]

		self.__load_standards_region_to_ds9_image(input_image, standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)

	def load_apass_standards_to_ds9_image(self, img, text1col=None, text2col=None,newds9=True):
		'''
		The standard stars are loaded from catalog file standard_reference_star_apass.txt

		See self.__load_standards_region_to_ds9_image
		'''
		apassstdfile = os.path.join(self.std_ref_dir, self.std_ref_stars_file['apass'])
		if not os.path.exists(apassstdfile):
			raise ValueError('standard star catalog %s not exists... you can prepare that with self.get_standards()'%apassstdfile)

		standards = Table.read(apassstdfile, format='ascii.csv')
		self.__load_standards_region_to_ds9_image(img, standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)

	def stat_apass_standards(self, outfile=None):
		'''
		display the mag VS. e_mag plot
		'''
		if self.standards is None:
			raise ValueError('prepare the standard catalog first')

		if outfile is not None:
			outfile = os.path.join(self.std_ref_dir, outfile)
		data = self.standards
		columns = ['Bmag','Vmag','r_mag','i_mag']
		simple_plot_table_data_four_columns(data, columns, outfile=outfile)

	def load_current_standards_to_ds9_image(self, imgkey=None, img=None, text1col=None, text2col=None, newds9=True, which_dir='raw'):
		'''
		The standard stars are from self.standards
		INPUTS:
			img: the input image (with absolute path)
			text1col: the column name in the standard star table
			text2col: the column name in the standard star table
		'''
		if imgkey is None and img is None:
			raise ValueError('no input image')
		else:
			if img is None:
				img = self.__get_internal_image(imgkey, which_dir =which_dir)

		if not self.standards:
			raise ValueError("please load standard catalog to self.standards first")

		self.__load_standards_region_to_ds9_image(img, self.standards, color='red', text1col=text1col, text2col=text2col, newds9=newds9)

	def __load_standards_region_to_ds9_image(self, img, standards, color='red', text1col=None, text2col=None, newds9=True):
		'''
		load the input standards to ds9 image; WCS info required in the image header

		INPUTS:
			img:
			standards: the input standards catalog
			newds9: whether start a new pyds9 new instance
		'''
		if standards.masked:
			standards = standards.filled(99.99)

		region_filename = os.path.join(self.std_ref_dir, 'standards_temp.reg')
		stdcolnames = standards.colnames
		if ('RAJ2000' in stdcolnames) and ('DEJ2000' in stdcolnames):
			RA_colname  = 'RAJ2000'
			Dec_colname = 'DEJ2000'
		elif ('RA_ICRS' in stdcolnames) and ('DE_ICRS' in stdcolnames): #gaia
			RA_colname  = 'RA_ICRS'
			Dec_colname = 'DE_ICRS'
		elif ('RA' in stdcolnames) and ('Dec' in stdcolnames):
			RA_colname  = 'RA'
			Dec_colname = 'Dec'
		elif  ('ra' in stdcolnames) and ('dec' in stdcolnames):
			RA_colname  = 'ra'
			Dec_colname = 'dec'
		else:
			raise ValueError("catalog columns names are not supported yet...")

		try:
			RAs  = standards[RA_colname].data
			Decs = standards[Dec_colname].data
		except:
			print "please check whether input parameter std_catalog is correct..."
		if text1col is not None and text1col not in stdcolnames:
			raise ValueError("the required colume %s not exist"%text1col)
		if text2col is not None and text2col not in stdcolnames:
			raise ValueError("the required colume %s not exist"%text2col)

		if text1col is not None and text2col is not None:
			loadtext = True
			textcol1 = 2
			textcol2 = 3
			input_sources = np.array([RAs, Decs, standards[text1col].data, standards[text2col].data]).transpose()
		elif text1col is not None:
			loadtext = True
			textcol1 = 2
			textcol2 = None
			input_sources = np.array([RAs, Decs, standards[text1col].data]).transpose()
		elif text2col is not None:
			loadtext = True
			textcol1 = None
			textcol2 = 2
			input_sources = np.array([RAs, Decs, standards[text2col].data]).transpose()
		else:
			loadtext = False
			textcol1 = None
			textcol2 = None
			input_sources = np.array([RAs, Decs]).transpose()

		create_ds9_region_file(input_sources, region_filename, x_col=0, y_col=1, coordinate_type = 'fk5', radec_deg = True, circle_radius = 0.0042, color=color, load_text=loadtext, text_uoffset=0.005, text_loffset=0.005, textcol1=textcol1, textcol2=textcol2) #circle radius around 15 arcsec = APASS aperture size for photometry

		if not os.path.exists(img):
			raise IOError("image %s doesn't exist"%img)

		d = self.__display_image_with_ds9(img, newds9=newds9)
		self.__load_source_regions_to_ds9(d, region_filename)

	def show_reference_stars_with_id_on_template_image_flt(self, flt, which_dir='raw_image'):
		'''
		This is used to display the reference stars with labeled ID which aids selecting corresponding stars in input images or from standards stars
		'''
		if len(self.templates) == 0:
			self.__find_template_imagekey()

		if flt not in self.templates.keys():
			raise ValueError("no template image avaiable for filter %s"%flt)

		tpl_imgkey = self.templates[flt]
		tpl_imgkey_s = tpl_imgkey.split('.')[0]
		refmag_filename = os.path.join(self.template_dir, "%s_%s_idxymag.txt"%(flt, tpl_imgkey_s))
		refmag_table = Table.read(refmag_filename, format='ascii.fixed_width')

		xymagids = np.transpose(np.array([refmag_table['x'], refmag_table['y'], refmag_table['mag'], refmag_table['num']]))
		input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)
		regionfile = os.path.join(self.template_dir, "%s_%s_idxymag.reg"%(flt, tpl_imgkey_s))

		self.__load_stars_xy_mag_to_ds9(input_image, xymagids, regionfile, mag_offset = 0, coortype = 'image', color = 'red', width=1, radius=10)

	def show_standards_for_wcs_not_available(self,flt, std_catalog='apass', radius_arcmin = 16, rotate_deg = 0, brightest_N = None, invert_x=False, invert_y = False, x_mid = None, y_mid = None, display_mags = False, which_dir='raw_image'):
		'''
		prepare and display standard stars on the template image for the purpose of 'mannually' matching standard stars with calibration stars on the template image

		INPUTS:
			flt:
			std_catalog: 'apass', '2mass'
			radius_arcmin:
			rotate_deg:
			brightest_N:
			invert_x:
			invert_y:
			x_mid:
			y_mid:
			display_mags:
		'''

		if std_catalog == 'apass':
			corresmap = {'RA':"RA", 'Dec':"Dec", 'V':"V", 'B':"B",'gp':"gp", 'rp':"rp", 'ip':"ip", 'Verr':"Verr", 'Berr':"Berr", 'gperr':"gperr", 'rperr':"rperr", 'iperr':"iperr" }
		elif std_catalog == '2mass':
			corresmap = {'RA':"ra", 'Dec':"dec", 'J':"j_m", 'Jerr':"j_msig", 'H':"h_m", 'Herr':"h_msig",'K':"k_m", 'Kerr':"k_msig", }
		else:
			raise ValueError("%s not supported..."%std_catalog)


		racen  = self.sn_ra_world_deg
		deccen = self.sn_dec_world_deg

		if self.standards is None:
			self.get_standards(catalog=std_catalog)

		ras  = self.standards[corresmap['RA']]
		decs = self.standards[corresmap['Dec']]

		#
		if flt == 'I':
			flt = 'ip'
			print "SDSS-i band data used for I band"
		if flt == 'R':
			flt = 'rp'
			print "SDSS-r band data used for R band"

		#mag_col = corresmap[flt]
		mags    = self.standards[corresmap[flt]]
		magerrs = self.standards[corresmap[flt+'err']]

		RAfactor = np.cos(deccen/180.0*np.pi)
		distances = np.sqrt(((ras-racen)*RAfactor)**2 + (decs-deccen)**2)  #approximation 

		radius_deg = radius_arcmin/60.

		ra_yes    = ras[distances<radius_deg]
		dec_yes  = decs[distances<radius_deg]
		mags_yes = mags[distances<radius_deg]
		magerrs_yes = magerrs[distances<radius_deg]

		self.__find_template_imagekey()
		tpl_imgkey = self.templates[flt]
		x0 = self.photometry_info[tpl_imgkey]['x']
		y0 = self.photometry_info[tpl_imgkey]['y']

		if x0 == 0 or y0==0:
			raise ValueError("the image coordinate for the target is required")

		if self.pixscale is None:
			raise ValueError("please specify the pixel scale")

		pixscale_deg = self.pixscale / 3600.

		x = (ra_yes - racen)*RAfactor/pixscale_deg
		y = (dec_yes - deccen)/pixscale_deg

		rotate_rad = rotate_deg/180.*np.pi

		xp =  x*np.cos(rotate_rad) + y*np.sin(rotate_rad)
		yp = -x*np.sin(rotate_rad) + y*np.cos(rotate_rad)

		X = x0-xp	#left east and right west
		Y = y0+yp

		if invert_x:
			if x_mid is None:
				x_mid = (np.max(X) + np.min(X))/2
			X = 2*x_mid - X
		if invert_y:
			if y_mid is None:
				y_mid = (np.max(Y) + np.min(Y))/2
			Y = 2*y_mid - Y

		if brightest_N is None:
			N = len(mags)
		else:
			N = brightest_N

		mags_sort = np.argsort(mags_yes)

		X_plot = X[mags_sort[:N]]
		Y_plot = Y[mags_sort[:N]]
		mags_plot = mags_yes[mags_sort[:N]]
		magerrs_plot = magerrs_yes[mags_sort[:N]]

		xymags = np.transpose(np.array([X_plot, Y_plot, mags_plot, magerrs_plot]))
		input_image = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)

		tpl_skey = tpl_imgkey.split('.')[0]
		std_xymags_file = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.txt"%(flt, tpl_skey))
		np.savetxt(std_xymags_file, xymags, fmt="%8.2f %8.2f %8.3f %8.3f")

		regionfile = os.path.join(self.std_ref_dir, "%s_%s_std_xymag.reg"%(flt, tpl_skey))
		self.__load_stars_xy_mag_to_ds9(input_image, xymags, regionfile, mag_offset = 0, coortype = 'image', color = 'red', width=1, radius=10)


		return xymags

#DS9 Display
	def select_or_delete_rows_from_given_table_by_pick_on_ds9_display(self, image, intable, xcol=0, ycol=1, mode='s', match_criteria = 3, coordinate_type='image', radec_deg=True, circle_radius=10):
		'''
		INPUTS:
			image: the image on which source regions are displayed
			intable:
			xcol, ycol: column name or the column index
			mode: 's' or 'd': the picked up points will be selected from the input table or deleted from  input table
			match_criteria: if None, then self.criteria_point_match_general will be used
			coordinate_type:
			radec_deg: whether the RA,Dec (x and y column when the coordinates are in physical world) are in unit of degree
		'''
		try:
			if (xcol not in intable.colnames) or (ycol not in intable.colnames):
				xcol = intable.colnames[xcol]
				ycol = intable.colnames[ycol]

			xyids = np.array([intable[xcol].data, intable[ycol].data, np.arange(len(intable))]).transpose()
		except:
			raise ValueError("The input columns are not recognised")
		regionfile = 'temp_delete_after.reg'
		create_ds9_region_file(xyids, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = coordinate_type, radec_deg = radec_deg, circle_radius = circle_radius, color='green', width=1, load_text=True, textcol1=2)

		if not os.path.exists(image):
			raise ValueError("The input image %s not exists"%image)
		d = self.__display_image_with_ds9(image, newds9=True)
		self.__load_source_regions_to_ds9(d, regionfile)
		xy = xyids[0,[0,1]]
		self.__load_source_regions_to_ds9(d,xy,radius =20, color='red', width=2)

		indexs = []

		if match_criteria is None:
			match_criteria = self.criteria_point_match_general

		while True:
			xy = self.__get_xy_on_image_from_ds9(image, regionfile=None, newds9=False)
			if xy is None:
				break
			xy_target = [xy[0],xy[1]]
			yesfound, xy_match,index = self.__find_corresponding(xyids[:,[0,1]], xy_target, match_criteria)

			if (not yesfound):
				print "No object found within criteria..."
				print "You can change input parameter match_criteria or self.criteria_point_match_general to give a new search distance"
			elif len(index)>1:
				print 'more than one objects fall into the criteria region of reference star. This object is droped!'
				print "You can change input parameter match_criteria or self.criteria_point_match_general to give a new search distance"
			else:
				print "the corresponding object #%s found at (%s)"%(index[0], xy_match)
				indexs.append(index[0])

		if mode == 's':
			saveindexs = indexs
		elif mode == 'd':
			saveindexs = [i for i in np.arange(len(intable)) if i not in indexs]
		else:
			raise ValueError('invalid input for mode')
		outtable = intable[saveindexs]

		return outtable, indexs


	def __load_source_regions_to_ds9(self,ds9_object,input_reg,radius=10,color='green',width=1, verbose=1):
		'''
		Inputs:
			ds9_object:
			input_reg: ds9 region file or position list with two element (x,y)
			radius:
			color:
			width:
		'''
		if isinstance(input_reg,str) and os.path.isfile(input_reg):
			ds9_object.set('regions load %s'%input_reg)
		elif isinstance(input_reg,list):
			x = input_reg[0]
			y = input_reg[1]
			#region_set = 'regions command {circle %s %s %s #color=%s width=%s}'%(x,y,radius,color,width)
			#ds9_object.set(region_set)
			if verbose:
				print "circle %s %s %s"%(x,y,radius)
			ds9_object.set('regions', 'image; circle %s %s %s #color=%s width=%s'%(x,y,radius,color,width))
		else:
			try:
				xy = input_reg.ravel()
				x = xy[0]
				y = xy[1]
				ds9_object.set('regions', 'image; circle %s %s %s #color=%s width=%s'%(x,y,radius,color,width))
				#ds9_object.set('regions command {circle %s %s %s #color=%s width=%s}'%(x,y,radius,color,width))
			except:
				raise IOError("invalid input for ds9 region")

	def load_dophot_photret_to_ds9(self, image_key, xcol=2, ycol=3, text1col=4, text2col=1, \
	filtercol=1, filtertype='eq',filtervalue=1, coortype = 'image', color = 'green', width=1, radius=15, \
	newds9=True, which_dir='raw_image'):
		'''
		load phot object with different caption choice
		INPUTS:
			xcol, ycol: define the object coordinate
			text1col, text2col: the text columns to show
			filtercol: if None then no filtering
			filtertype: 'lt', less than; 'eq', equal to; 'gt', greater than
			filtervalue: the criteria for filtering
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		if self.dophot_version == 'fortran':
			photdir = self.psfphot_dophot_fortran_dir
		else:
			photdir = self.psfphot_dophot_c_dir
		imgkey = image_key.split('.')[0]
		photret_file = os.path.join(photdir, imgkey+'.out')
		photret = np.loadtxt(photret_file)

		if filtercol is not None:
			if filtertype == 'eq':
				mask = photret[:,filtercol] == filtervalue
			elif filtertype == 'lt':
				mask = photret[:,filtercol] < filtervalue
			elif filtertype == 'gt':
				mask = photret[:,filtercol] > filtervalue
			else:
				raise ValueError('not support filter type %s'%filtertype)
			photret = photret[mask]

		regionfile = os.path.join(photdir, imgkey+"_photret.reg")
		create_ds9_region_file(photret, regionfile, clobber = True,  x_col=xcol, y_col=ycol, x_offset=0, y_offset=0, coordinate_type = coortype, radec_deg = True, circle_radius = radius, color=color, width=width, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=text1col, textcol2=text2col)

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		d.set('regions load %s'%regionfile)


	def load_stars_xy_mag_to_ds9(self,image_key,photometry_method = 'apphot', \
	offset= 0, mag_cut = 1, xl = None, xh=None, yl=None, yh = None, \
	coortype = 'image', color = 'green', width=1, radius=15, \
	image_physical_diff_x = 0, image_physical_diff_y = 0, newds9=True, which_dir='raw_image'):
		'''
		load star region ans magnitude measurements to the image in DS9

		Inputs:
			image_key:
			photometry_method:
			offset: add an offset to magnitudes measurements
			mag_cut: bright %(mag_cut*100) will be displayed
			image_physical_diff_x:
			image_physical_diff_y:
		'''
		imgkey = image_key.split('.')[0]
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		if photometry_method == 'apphot':
			photret_file = os.path.join(self.aperture_photometry_dir, imgkey+ self.apphot_ret_file_suffix)
		elif photometry_method == 'psfphot':
			photret_file = self.__get_internal_psfphot_mag_file(image_key)
		else:
			raise IOError("input for photometry_method is not supported...")
		photdir = os.path.dirname(photret_file)
		xymags = np.loadtxt(photret_file)

		if xl is not None:
			xymags = xymags[xymags[:,0] > xl,:]
		if xh is not None:
			xymags = xymags[xymags[:,0] < xh,:]
		if yl is not None:
			xymags = xymags[xymags[:,1] > yl,:]
		if yh is not None:
			xymags = xymags[xymags[:,1] < yh,:]

		xymags = xymags[np.argsort(xymags[:,2]),:]
		N = len(xymags)
		N_want = int(np.ceil(N*mag_cut))
		mask = range(N_want)
		xymags_show = xymags[mask,:]

		regionfile = os.path.join(photdir, imgkey+"_xymag.reg")
		self.ds9D  = self.__load_stars_xy_mag_to_ds9(input_image, xymags_show, regionfile, mag_offset = offset, coortype = coortype, color =color, width=width, radius=radius, newds9=newds9)

	def __load_stars_xy_mag_to_ds9(self, input_image, xymags, regionfile, mag_offset  = 0, coortype = 'image', color = 'green', width=1, radius=15, newds9=True):
		'''
		load object list to ds9
		'''
		if mag_offset != 0:
			xymags[:,2] = xymags[:,2] + mag_offset
		create_ds9_region_file(xymags, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = coortype, radec_deg = True, circle_radius = radius, color=color, width=width, load_text=True,  text_uoffset=25, text_loffset=25, textcol1=2, textcol2=3)

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		d.set('regions load %s'%regionfile)

		return d

	def load_matched_stars_to_ds9(self, img, matched_source_file = None, ref_reg_file = None, input_reg_file = None):
		'''
		load two set sources which are matched one to one.

		INPUTS:
			matched_source_file:	first and second columns are images coordinates of ref sources; and fifth and sixth columns are images coordinates of input sources
			ref_reg_file:
			input_reg_file:
		'''

		#check inputs
		if ref_reg_file is not None and input_reg_file is not None:
			if os.path.exists(ref_reg_file) and os.path.exists(input_reg_file):
				regionfile_provided = True
			else:
				raise IOError("input region file %s and/or %s doesn't exist..."%(ref_reg_file, input_reg_file))
		else:
			regionfile_provided = False

		if matched_source_file is None and (not regionfile_provided):
			raise IOError("No valid region file provided...")

		#explicit ds9  region files have higher priority
		if regionfile_provided:
			self.__load_matched_stars_to_ds9(img, ref_reg_file, input_reg_file)
		else:
			ref_reg_file_temp = 'temp_ref_sources.reg'
			input_reg_file_temp = 'temp_input_source.reg'
			input_sources = np.loadtxt(matched_source_file)

			create_ds9_region_file(input_sources, ref_reg_file_temp,   clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 20)
			create_ds9_region_file(input_sources, input_reg_file_temp, clobber = True,  x_col=4, y_col=5, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

			self.__load_matched_stars_to_ds9(img, ref_reg_file_temp, input_reg_file_temp)

	def __load_matched_stars_to_ds9(self, input_image, ref_reg_file, input_reg_file, newds9=True):
		'''
		load ds9 region files 'ref_reg_file' and 'input_reg_file' to image 'img'
		'''
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		print "reference with large circle"
		self.__load_source_regions_to_ds9(d, ref_reg_file)
		self.__load_source_regions_to_ds9(d, input_reg_file)

	def load_std_stars_region_to_ds9_flt(self,flt,load_mags=False, on_image_after_astrometry = False, newds9=True, updatetable=1):
		'''
		load standard reference stars to the template image
		Here standard stars are all available ones before match with sources measured on the science image
		INPUTs:
			flt:
			load_mags: if True, then load the magnitudes
			on_image_after_astrometry: if True, then the image in ./warehouse/template/cal_xxx.fits will be used
		'''
		self.__find_template_imagekey(updatetable=updatetable)
		if not on_image_after_astrometry:
			image_key = self.templates[flt]
        		input_image = self.images[image_key]
		else:
			input_image = self.templates_after_astrometry[flt]

		d = self.__display_image_with_ds9(input_image, newds9=newds9)

		regfile = os.path.join(self.std_ref_dir,flt+self.regfile_suffix)
		if not os.path.exists(regfile):
			raise ValueError("region file %s not exist"%regfile)
		self.__load_source_regions_to_ds9(d,regfile)

		return d

	def __add_text_to_ds9(self,d,x,y,text,color='red',width=2, physical_image_offset = True):
		'''
		add 'text' to ds9 frame d at position (x,y)

		Notes:
			Use self.physical_image_offset_x and self.physical_image_offset_y to deal with the difference between image coordinate and physical coordinate
		'''
		if physical_image_offset:
			x_diff = self.physical_image_offset_x
			y_diff = self.physical_image_offset_y
		else:
			x_diff = 0
			y_diff = 0

		x_physical = x + x_diff
		y_physical = y + y_diff
		d.set('regions command {text %s %s #text="%s" color="%s" width=%s}'%(x_physical,y_physical,text,color,width))

	def load_fwhm_stars_to_ds9(self, image_key, which_dir='raw_image', newds9=True):
		'''
		load the star list from which FWHM value derived to ds9
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		fwhmstarsfile= os.path.join(self.stars_dir, image_key.split('.')[0]+'_fwhm.stars')
		source_regfile = os.path.join(self.stars_dir, image_key.split('.')[0] + '_fwhmstars.reg')
		fwhmstars = np.loadtxt(fwhmstarsfile)
		create_ds9_region_file(fwhmstars, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)

		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		self.__load_source_regions_to_ds9(d,source_regfile)

		return d


	def load_psf_stars_to_ds9(self, image_key, which_dir='raw_image', newds9=True):
		'''
		load the stars from which the PSF model is built
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		psfstarfile = os.path.join(self.psfmodel_dir, image_key.split('.')[0]+'.pst.1')
		psfstars = Table.read(psfstarfile, format='daophot')
		psfstars_xys = np.array([psfstars['XCENTER'].data, psfstars['YCENTER'].data]).transpose()
		source_regfile = os.path.join(self.psfmodel_dir, image_key.split('.')[0]+'_psfstars.reg')
		create_ds9_region_file(psfstars_xys, source_regfile, clobber = True,  x_col=0, y_col=1, x_offset=1, y_offset=1, coordinate_type = 'image', circle_radius = 10)
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		self.__load_source_regions_to_ds9(d,source_regfile)

		return d


	def load_source_regions_to_ds9(self,image_key,regtype = 'prephot',whichone='staronly', which_dir='raw_image', newds9=True):

		'''
		load source list to ds9

		Input:
			regtype: 'prephot' or 'postphot'
			whichone:'staronly' or 'all', only work for psfphotmetry
				 'staronly' mean identification restricted to stars
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		d = self.__display_image_with_ds9(input_image, newds9=newds9)

		if regtype == 'prephot':
			source_regfile = os.path.join(self.stars_dir, image_key.split('.')[0] + self.regfile_suffix)
		elif regtype == 'postphot':
			if self.photometry_method is None:
				self.photometry_method = raw_input('Please input the photometry method used(apphot/psfphot):')
			if self.photometry_method == 'apphot':
				source_regfile = os.path.join(self.aperture_photometry_dir, image_key.split('.')[0] + self.regfile_suffix)
				if not os.path.isfile(source_regfile):
					raise IOError("region file %s doesn't exist"%source_regfile)
			elif self.photometry_method == 'psfphot':
				psf_photometry_dir = self.__get_internal_psfphot_dir()
				if whichone == 'staronly':
					source_regfile = os.path.join(psf_photometry_dir, image_key.split('.')[0] + self.regfile_suffix)
				else:
					source_regfile = os.path.join(psf_photometry_dir, image_key.split('.')[0] + '_all'+ self.regfile_suffix)

				if not os.path.isfile(source_regfile):
					raise IOError("region file %s doesn't exist"%source_regfile)
			else:
				raise IOError('Invalid input for photometry_method')
		else:
			raise IOError('Invalid input for regtype')

		self.__load_source_regions_to_ds9(d,source_regfile)

		x = self.photometry_info[image_key]['x']
		y = self.photometry_info[image_key]['y']
		xy = [x,y]
		if x != 0.0 and y != 0.0:
			self.__load_source_regions_to_ds9(d,xy,radius =20, color='red', width=2)

		return d

	def display_multi_images_and_highlight_sn(self,flts = None,jd_start=None,jd_end = None, display_drop_with_flag = 0, updatetable=1, which_dir='raw'):
		'''

		display images with sn position labeled as well obstime and flt information shown
		Default is all registered images will be displayed unless flts/jd_start/jd_end is given
		'''
		if updatetable:
			self.__dict2table()
		flts_valid = np.unique(self.photometry_record_table['flt'])

		if flts is None:
			img_table_flts = self.photometry_record_table
		else:
			img_table_flts = Table()

			if isinstance(flts,str):
				flts = [flts]
			for flt in flts:
				if flt not in flts_valid:
					raise KeyError("no image taken with %s filter"%flt)
				table_toadd = self.__select_rows_from_table(self.photometry_record_table,'flt',flt)
				img_table_flts = vstack([img_table_flts, table_toadd])

		img_table_flts_jd = img_table_flts
		if jd_start is not None:
			img_table_flts_jd = self.__select_rows_from_table(img_table_flts_jd, 'obstime',jd_start,mode='gt' )
 		if jd_end is not None:
			img_table_flts_jd = self.__select_rows_from_table(img_table_flts_jd, 'obstime',jd_end,mode='lt' )

		img_table_flts_jd_dropflag = self.__select_rows_from_table(img_table_flts_jd, 'drop', display_drop_with_flag)
		img_keys = img_table_flts_jd_dropflag['name']

		self.__display_multi_images_and_highlight_sn(img_keys, which_dir=which_dir)

	def __display_multi_images_and_highlight_sn(self,img_keys, which_dir='raw'):
		'''
		display multiple images in ds9
		'''
		pioneer_img = img_keys[0]

		d = self.display_image_with_ds9(pioneer_img, which_dir=which_dir)
		x = self.photometry_info[pioneer_img]['x']
		y = self.photometry_info[pioneer_img]['y']
		xy_reg = [x,y]
		self._photometry__load_source_regions_to_ds9(d, xy_reg, radius=25, color='red', width=2)

		x_text = x + 100
		y_text = y + 100
		jd = self.photometry_info[pioneer_img]['obstime']
		t = Time(jd,scale='utc',format='jd')
		t_isot = t.isot
		fwhm_this = self.photometry_info[pioneer_img]['fwhm']
		flt_this = self.photometry_info[pioneer_img]['flt']
		text  = flt_this + '@' + t_isot + 'with fwhm=%s'%fwhm_this
		self._photometry__add_text_to_ds9(d,x_text,y_text,text,color='green',width=2,)

		for img in img_keys:
		        if img == pioneer_img:
		                continue
		        d.set('frame new')
		        #current_img = self.images[img]
			current_img = self.__get_internal_image(img, which_dir=which_dir)
		        d.set('file %s'%current_img)
		        d.set('zoom to fit')

		        x = self.photometry_info[img]['x']
		        y = self.photometry_info[img]['y']
		        xy_reg = [x,y]
		        self._photometry__load_source_regions_to_ds9(d, xy_reg, radius=25, color='red', width=2)

		        x_text = x + 100
		        y_text = y + 100
		        jd = self.photometry_info[img]['obstime']
		        t = Time(jd,scale='utc',format='jd')
		        t_isot = t.isot
			fwhm_this = self.photometry_info[img]['fwhm']
		        flt_this = self.photometry_info[img]['flt']
		        text  = flt_this + '@' + t_isot + ' with fwhm=%s'%fwhm_this
		        self._photometry__add_text_to_ds9(d,x_text,y_text,text,color='green',width=2,)

	def __display_image_with_ds9(self,input_image, newds9=True):
		'''
		display input_image in ds9
		'''
		if newds9:
             		d = pyds9.DS9()
               		d.set('fits %s'%input_image)
               		d.set('zoom to fit')
                	d.set('scale zscale')
			self.ds9D = d
		else:
			if not self.ds9D:
				raise ValueError("no available ds9 instance in self.ds9D")
			d = self.ds9D

        	return d

	def display_image_with_ds9(self,image_key, which_dir = 'raw_image', newds9=True):
        	'''
		display image with pyds9
		'''
		input_image = self.__get_internal_image(image_key, which_dir=which_dir)
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		return d

	def __load_check_and_decidedrop_single_image(self, imgkey, image_dir='raw_image'):
		'''
		load image with imgkey to ds9 and decide on the drop
		'''
		d = self.display_image_with_ds9(imgkey,which_dir=image_dir)
		x = self.photometry_info[imgkey]['x']
		y = self.photometry_info[imgkey]['y']
		xy_reg = [x,y]
		self.__load_source_regions_to_ds9(d, xy_reg, radius=25, color='red', width=2)
		drop_signal_default = 0
		drop_signal = raw_input("Enter the drop signal for image %s (default: %s):" %(imgkey,drop_signal_default)) or drop_signal_default

		if not str(drop_signal).isdigit():
			raise IOError("Invalid input for drop_signal")

		drop_sinal = int(drop_signal)
		#if drop_signal:
		self.photometry_info[imgkey]['drop'] = drop_signal

	def label_single_xy_region(self,image,xy):
		d = self.__display_image_with_ds9(image)
		d = self.__load_source_regions_to_ds9(d,xy,radius=15,color='red',width=2)
		return d

	def __get_xy_on_image_from_ds9(self,input_image, regionfile =None, newds9=True):
		'''
		pick up xy coordinate from ds9 image
		'''
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		if regionfile is not None:
			print regionfile
			d.set('regions load %s'%regionfile)

		xystr = d.get('imexam coordinate image')
		if xystr == '' or xystr == '0 0':
			xy = None
		else:
			x,y = xystr.split()
			xy = np.array([float(x),float(y)])
			print xy
			self.__load_source_regions_to_ds9(d,[float(x),float(y)],radius=20,color='red',width=2, verbose=0)

		return xy

	def display_sources_with_ds9(self, imgkey, sources, xcol='x', ycol='y', radius=10,color='green',width=1, which_dir='raw_image', newds9=True):
		'''
		Display sources in 'sources' on image with 'imgkey' from directory 'which_dir'


		'''
		if hasattr(sources, 'keys'):
			if xcol in sources.keys() and ycol in sources.keys():
				xdata = sources[xcol]
				ydata = sources[ycol]
			else:
				xdata = sources[sources.keys()[0]]
				ydata = sources[sources.keys()[1]]
		elif isinstance(sources, np.ndarray):
			if len(sources.shape) == 1:
				xdata = sources[0]
				ydata = sources[1]
			else:
				xdata = sources[:,0]
				ydata = sources[:,1]
		else:
			try:
				xdata = sources[0]
				ydata = sources[1]
			except:
				raise IOError("the input sources format not recognised")

		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		d = self.__display_image_with_ds9(input_image, newds9=newds9)
		for x,y in zip(xdata, ydata):
			xy = [x,y]
			self.__load_source_regions_to_ds9(d, xy, radius=radius,color=color,width=width, verbose=0)

#Common
	def __get_internal_image(self, image_key, which_dir='raw_image'):
		'''
		GET the corresponding filename of different kind defined by which_dir
		INPUTS:
			which_dir, 'raw_iamge' or 'raw':
				   'modified_image' or 'md':
				   'interp_image' or 'interp':
				   'conv_image' or 'conv':
				   'subtraced_image' or 'sub': image subtraction output image
				   'astrometry_image' or 'wcs':
				   'psfsub_dophot_fortran' or 'psdf': PSF photometry output residual image from DoPhot fortran version
				   'psfsub_dophot_c' or 'psdc': PSF photometry output residual image from DoPhot C version
				   'psfsub_daophot_iraf' or 'psdi': PSF photometry output residual image from IRAF
				   'psf_subras': emperical psf subtraster from dophot_C
				   'psf_model' or 'psf': PSF model image
				   'iso_image' or 'iso': isophote model image
				   'galsub_image' or 'galsub': galaxy model subtracted image
				   'mask_image' or 'mask':
				   'astrometry' or 'wcs': images with new or updated WCS solution
		'''
		imgkey_s = image_key.split('.')[0]
		if which_dir == 'raw_image' or which_dir == 'raw':
			input_image = self.images[image_key]
		elif which_dir == 'modified_image' or which_dir == 'md':
			input_image = os.path.join(self.modified_image_dir, image_key)
		elif which_dir == 'subtracted_image' or which_dir == 'sub':
			input_image = os.path.join(self.subtraction_dir, 'sub_'+image_key)
		elif which_dir == 'interp_image' or which_dir == 'interp':
			input_image = os.path.join(self.subtraction_dir, 'interp_'+image_key)
		elif which_dir == 'conv_image' or which_dir == 'conv':
			input_image = os.path.join(self.subtraction_dir, 'conv_'+image_key)
		elif which_dir == 'astrometry_image' or which_dir =='wcs':
			input_image = os.path.join(self.template_dir,'cal_'+image_key)
		elif which_dir == 'psf_model' or which_dir == 'psf':
			input_image = os.path.join(self.psfmodel_dir, imgkey_s+'_psf_image.fits')
		elif which_dir == 'psfsub_dophot_fortran' or which_dir =='psdf':
			input_image = os.path.join(self.psfphot_dophot_fortran_dir, 'out_'+image_key)
		elif which_dir == 'psfsub_dophot_c' or which_dir =='psdc':
			input_image = os.path.join(self.psfphot_dophot_c_dir, 'out_'+image_key)
		elif which_dir == 'psfsub_daophot_iraf' or which_dir == 'psdi':
			input_image = os.path.join(self.psfphot_daophot_iraf_dir , imgkey_s + '_sub.fits')
		elif which_dir == 'psf_subras':
			input_image = os.path.join(self.psfphot_dophot_c_dir, imgkey_s + '_psf.fits')
		elif which_dir == 'iso_image' or which_dir == 'iso':
			input_image = os.path.join(self.isophote_dir, imgkey_s + '_iso.fits')
		elif which_dir == 'galsub_image' or which_dir == 'galsub':
			input_image = os.path.join(self.isophote_dir, imgkey_s + '_gal_sub.fits')
		elif which_dir == 'mask_image' or which_dir == 'mask':
			input_image = os.path.join(self.ds9_region_dir, 'mask_'+image_key)
		elif which_dir == 'astrometry' or which_dir == 'wcs':
			input_image = os.path.join(self.astrometry_dir, image_key)
		else:
			raise IOError('Non valiad input for which_dir...')

		if not os.path.exists(input_image):
			print "warning: the request image does not exist"

		return input_image

	def __delete_file_if_exist(self,filename):
		'''
		if file 'filename' exist, remove the existing file
		'''
		if os.path.isfile(filename):
			os.remove(filename)

	def __select_rows_from_table(self,table,colname,value,mode='eq',criteria = 0):
		'''
		select rows in table which meet specific requirements given by mode and criteria
		available mode 'eq','lt','gt','lgt'

		'eq':X == value; 'lt': X<= value; 'gt': X>=value; 'lgt': value-criteria  <= X <= value+criteria
		if mode is 'lgt' then non-zero criteria is needed, otherwise 'lgt' is the same as 'eq' with default criteria

		return the table by rows which meet the criteria
		'''
		colvalues = table[colname]
		if mode   == 'eq':
			mask  = colvalues == value
		elif mode == 'lt':
			mask  = colvalues <= value
		elif mode == 'gt':
			mask  = colvalues >= value
		elif mode == 'lgt':
			if criteria == 0:
				print "Attention: if mode is 'lgt' then non-zero criteria is needed, otherwise 'lgt' is the same as 'eq' with default criteria!!!"
			mask = np.logical_and(colvalues>=value-criteria,colvalues<=value+criteria)
		else:
			raise KeyError("the criteria %s is not available"%criteria)

		selected = table[mask]

		return selected

	def __select_columns_from_table_to_ascii(self, table, colnames, outfile):
		'''
		output data in colnames to pure ascii data file
		'''
		outdata = table[colnames]
		self.__delete_file_if_exist(outfile)
		outdata.write(outfile, format='ascii.fast_commented_header')


	def get_image_table(self, updatetable=1):
		'''
		get image table: obsdate, B, V, ...
		'''
		if updatetable:
			self.__dict2table()

		flts = np.unique(self.photometry_record_table['flt'])

		if self.current_telescope == 'LCOGT':
			obsdates = np.unique([img.split('-')[2] for img in self.photometry_record_table['realimg']])
		elif self.current_telescope == 'SMARTS':
			obsdates = np.unique([img[4:10] for img in self.photometry_record_table['realimg']])
		else:
			print "current telescope not supported yet"
			return

	def create_empty_image(self, outimg, NX=1000, NY=1000, mean=0, std=1, refimg=None):
		'''
		create an "empty" image
		
		'''
		if refimg is None:
			data = np.random.randn(NX, NY)*std+mean
			hdu = fits.PrimaryHDU(data=data)	
			hdu.writeto(outimg)	
		else:
			try:
				refimghdu = fits.open(refimg)
			except Exception as e:
				raise IOError(e)
			
			refimgdata = refimghdu[0].data
			refimghdu[0].data =mean+std*np.random.randn(refimgdata.shape[0], refimgdata.shape[1])
			
			refimghdu.writeto(outimg)	
		

	def create_mask_image(self, imgkey, mask_image_dir=None, base_image_dir='raw_image', \
	region_type='square', square_half_length=10 ):
		'''
		star with simple, box region; can be sophisticated as polygon

		INPUTS:
			imgkey:
			mask_image_dir: the folder to search for and store the mask image; if None, save the output to self.warehouse_dir+'/regions'
			base_image_dir: if no existing mask image, then go to thid folder to find the base image
			region_type: square, rectangle, polygon
			square_half_length:
		'''
		from matplotlib import path as mplpath
		baseimagefile = self.__get_internal_image(imgkey, which_dir=base_image_dir)
		baseimage = fits.open(baseimagefile)[0].data

		maskimage = np.zeros(baseimage.shape)
		if mask_image_dir is not None:
			if not os.path.exists(mask_image_dir):
				mask_image_dir  = os.path.dirname(self.__get_internal_image(imgkey, which_dir=mask_image_dir))
			maskimagefile = os.path.join(mask_image_dir, 'mask_'+imgkey)
			if os.path.exists(maskimagefile):
				maskimage = fits.open(maskimagefile)[0].data
		else:
			maskimagefile = self.__get_internal_image(imgkey, which_dir='mask_image')

		if baseimage.shape != maskimage.shape:
			raise ValueError('the base image and the mask image have different shape...')

		NX, NY = maskimage.shape
		indxs,indys = np.meshgrid(np.arange(NX), np.arange(NY))
		points = np.array([indxs.flatten(), indys.flatten()]).transpose()

		pickpts = []
		firstdisplay = True
		while True:
			xy = self.__get_xy_on_image_from_ds9(baseimagefile, regionfile=None, newds9=firstdisplay)
			firstdisplay = False
			if xy is None:
				break
			xy_0based = [xy[0]-1,xy[1]-1]
			pickpts.append(xy_0based)

		if region_type == 'square':
			for xc,yc in pickpts:
				hl = square_half_length
				vertexs = [[yc-hl, xc-hl], [yc+hl, xc-hl], [yc+hl, xc+hl], [yc-hl, xc+hl]]
				square = mplpath.Path(vertexs)
				inside = square.contains_points(points)
				for idx,idy in points[inside]:
					maskimage[idx,idy] = 1
		elif region_type == 'rectangle':
			for i, (blcx, blcy) in enumerate(pickpts): #bottom left corner point
				if np.mod(i,2) == 0:
					trcx, trcy = pickpts[i+1] #top right corner point
					vertexs = [[blcy, blcx],[trcy, blcx],[trcy,trcx],[blcy,trcx]]
					rectangle = mplpath.Path(vertexs)
					inside = rectangle.contains_points(points)
					for idx,idy in points[inside]:
						maskimage[idx,idy] = 1
				else:
					continue
		elif region_type == 'polygon':
			vertexs = [[vty, vtx] for (vtx, vty) in pickpts]
			polygon = mplpath.Path(vertexs)
			inside = polygon.contains_points(points)
			for idx,idy in points[inside]:
				maskimage[idx,idy] = 1
		else:
			print 'sorry, region type %s not supported'%region_type
			return

		maskimage_hdu = fits.PrimaryHDU(data=maskimage)
		maskimage_hdu.writeto(maskimagefile, overwrite=True)



#Souce match
	def match_sources_from_two_catalog(self, incatname, refcatname):
		'''
		Match two catalogs by coordinate
		'''
		if incatname not in self.std_ref_stars_file.keys() or refcatname not in self.std_ref_stars_file.keys():
			raise ValueError('The requested catalog not registered in self.std_ref_stars_file')

		incat  = Table.read(os.path.join(self.std_ref_dir, self.std_ref_stars_file[incatname]), format='ascii.csv')
		refcat = Table.read(os.path.join(self.std_ref_dir, self.std_ref_stars_file[refcatname]), format='ascii.csv')

		inracol  = self.std_catalog_colnames[incatname]['RA']
		indeccol = self.std_catalog_colnames[incatname]['Dec']
		refracol = self.std_catalog_colnames[refcatname]['RA']
		refdeccol= self.std_catalog_colnames[refcatname]['Dec']
		outfile = os.path.join(self.std_ref_dir, "%s_%s_matched_radecid.txt"%(incatname, refcatname))
		incat_mtc, refcat_mtc = self.match_sources_from_two_catalog_general(incat, refcat, inracol, indeccol, refracol, refdeccol, outfile=outfile)

		return incat_mtc, refcat_mtc

	def match_sources_from_two_catalog_general(self, catalog1, catalog2, racol1, deccol1, racol2, deccol2, outfile=None, criteria=1.0):
		'''
		INPUTS:
			catalog1, catalog2:
			racol1, deccol1:
			racol2, deccol2:
			criteria: in arcsec
		'''
		incoordata  = np.array([catalog1[racol1].data, catalog1[deccol1].data, np.arange(len(catalog1))]).transpose()
		refcoordata = np.array([catalog2[racol2].data, catalog2[deccol2].data, np.arange(len(catalog2))]).transpose()
		if outfile is None:
			outfile = 'cat1_cat2_match_ret.txt'

		
		self.fitsh_grmatch_type = 'coord'
		self.fitsh_grmatch_coordmatch_pars['--max-distance'] = criteria/3600.0
		matched_radecid = self.__fitsh_grmatch(incoordata, refcoordata, outfile, mode ='data') #

		matchdata1 = catalog1[map(int, matched_radecid[:,2])]
		matchdata2 = catalog2[map(int, matched_radecid[:,5])]

		return matchdata1, matchdata2

	def match_sources_on_two_images(self, in_imgkey, ref_imgkey, option=None,mode='file', ref_col1=1, ref_col2=2, input_col1=1, input_col2=2, input_order_col=-3, ref_order_col=-3, max_distance=1, matched_xyxy_outfile=None, matched_trans_output=None,  match_result_display=False, show_match_stat=False):
		'''
		match two set of sources detected on input image and reference image

		INPUTS:
			matched_xyxy_outfile:  this is the simplified version of the matched source list with only x,y from two source list


		MIDDLE PRODUCTS:
			match_output: raw matched source list with whole information from input source list and reference source list
			match_trans_output: the matched and transformed source list
		'''
		in_imgkey_s = in_imgkey.split('.')[0]
		ref_imgkey_s = ref_imgkey.split('.')[0]
		ref_list = os.path.join(self.stars_dir, ref_imgkey_s + self.starfile_suffix)
		input_list = os.path.join(self.stars_dir, in_imgkey_s + self.starfile_suffix)
		match_output = os.path.join(self.stars_dir, "match_%s_%s.txt"%(in_imgkey_s, ref_imgkey_s))
		transcoeff_output = os.path.join(self.stars_dir, "match_%s_%s.coeff"%(in_imgkey_s, ref_imgkey_s))

		reflist_data = np.loadtxt(ref_list)
		Nr_ref, Nc_ref = reflist_data.shape

		#check here
		self.fitsh_grmatch_type = 'point'
		matched_data = self.__fitsh_grmatch(ref_list,input_list, match_output, trans_output=transcoeff_output, mode =mode)
		wanted_columns = [ref_col1-1, ref_col2-1, input_col1+Nc_ref-1, input_col2+Nc_ref-1]
		matched_xyxy = matched_data[:,wanted_columns]

		if matched_xyxy_outfile is None:
			matched_xyxy_outfile = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(in_imgkey_s, ref_imgkey_s))
		np.savetxt(matched_xyxy_outfile, matched_xyxy, fmt="%8.2f %8.2f %8.2f %8.2f")

		porder, dxfit_values, dyfit_values = self.__extract_transformation_fitting(transcoeff_output)
		ref_xys = matched_xyxy[:,[0,1]]
		input_xys = matched_xyxy[:,[2,3]]
		transformed_xys = self.__transform_xys(dxfit_values, dyfit_values, porder, ref_xys)

		matched_xyxy_trans = np.hstack((transformed_xys, input_xys))
		if matched_trans_output is None:
			matched_trans_output = os.path.join(self.stars_dir, "match_%s_%s_xyxy_trans.txt"%(in_imgkey_s, ref_imgkey_s))
		np.savetxt(matched_trans_output, matched_xyxy_trans, fmt="%8.2f %8.2f %8.2f %8.2f")

		if match_result_display:
			xy_simple_scatter_plot_2sets(input_xys[:,0], input_xys[:,1], transformed_xys[:,0], transformed_xys[:,1], label1='DATA1', label2='reference transformed', xlabel='X', ylabel='Y')

		if show_match_stat:
			sources_match_result_stat_plot(input_xys, transformed_xys)

#Image subtraction
	def image_subtraction_all(self, which_dir='raw_image', xstamp=8, ystamp=8, renew_sub=False):
		self.__dict2table()
		self.__find_template_imagekey()

		faildict = {}
		for flt,tpl_imgkey in self.templates.items():
			failimgs = self.image_subtraction_flt(flt, tpl_imgkey, which_dir=which_dir, xstamp=xstamp, ystamp=ystamp, renew_sub=renew_sub)
			for img in failimgs:
				self.photometry_info[img]['drop'] = 13
			faildict[flt]= failimgs
			print "The following images failed on image subtraction:"
			print failimgs

		return faildict

	def image_subtraction_flt(self, flt, tpl_imgkey, updatetable=0, xstamp=8, ystamp=8, \
	renew_sub=False, apply_mask=False, which_dir='raw_image',  onlydothese=None):
		'''
		image subtraction process
		'''

		if updatetable:
			self.__dict2table()

		self.hotpants_pars['-nsx']= xstamp
		self.hotpants_pars['-nsy']= ystamp


		fltimgs = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		faillist = []
		for imgkey in fltimgs['name']:
			print imgkey
			if onlydothese is not None:
				if imgkey not in onlydothese:
					print('skip this...')
					continue
			subimg = os.path.join(self.subtraction_dir, 'sub_'+imgkey)
			if os.path.exists(subimg) and (not renew_sub):
				continue

			if self.photometry_info[imgkey]['drop'] != 0:
				continue
			if imgkey == tpl_imgkey:
				continue

			try:
				self.match_sources_on_two_images(imgkey, tpl_imgkey, match_result_display=0, show_match_stat=0, option=None)
			except:
				faillist.append(imgkey)
				continue
			if self.photometry_info[imgkey]['bkg']<100:
				self.add_constant_to_image(imgkey, 100, input_which_dir=which_dir)
			#self.align_and_resample_image_given_matched_source_list(imgkey, tpl_imgkey, which_dir=which_dir)
			self.align_and_resample_image_given_transfile(imgkey, tpl_imgkey, refimg_dir=which_dir, inputimg_dir=which_dir)
			if apply_mask:
				maskimg = self.__get_internal_image(tpl_imgkey, which_dir='mask')
				if not os.path.exists(maskimg):
					raise ValueError('Sorry, please prepare mask image %s first'%maskimg)
				inputimg_mask = maskimg
				refimg_mask = maskimg
			else:
				inputimg_mask = None
				refimg_mask = None
			self.hotpants_image_subtraction_single(imgkey, tpl_imgkey=tpl_imgkey, inputimg_mask=inputimg_mask, refimg_mask=refimg_mask) #, xstamp=xstamp,ystamp=ystamp, ilthresh=1, tlthresh=1, iuthresh=60000,tuthresh=60000)

		return faillist


	def image_subtraction_stdcal_flt(self, flt, tpl_imgkey, skiplist=None):
		'''
		work for hotpants results
		'''
		if self.photometry_info[tpl_imgkey]['calmag'] == self.photometry_info_init_values['calmag']:
			raise ValueError('calibration for the template image not avaiable')
		offset = self.photometry_info[tpl_imgkey]['calmag'] - self.photometry_info[tpl_imgkey]['instmag']
		offset_err = np.sqrt(self.photometry_info[tpl_imgkey]['calmagerr']**2 - self.photometry_info[tpl_imgkey]['instmagerr']**2)

		for image_key in self.images.keys():
			if skiplist is not None:
				if image_key in skiplist:
					continue
			if image_key == tpl_imgkey:
				continue
			if self.photometry_info[image_key]['drop'] !=0 or self.photometry_info[image_key]['flt']!=flt:
				continue

			calmag = self.photometry_info[image_key]['instmag'] + offset
			calmagerr = np.sqrt(self.photometry_info[image_key]['instmagerr']**2 + offset_err**2)

			self.photometry_info[image_key]['calmag'] = calmag
			self.photometry_info[image_key]['calmagerr'] = calmagerr

	def image_subtraction_get_target_xy_on_subtracted_images(self, flt, tpl_imgkey, which_dir='raw_image', showplot=True):
		'''
		find the source coordinate on the template image
		'''
		tplimg = self.__get_internal_image(tpl_imgkey, which_dir=which_dir)
		xy = self.__get_xy_on_image_from_wcsinfo(tplimg)
		xc0 = np.round(xy[0])
		yc0 = np.round(xy[1])

		xs = []
		ys = []
		for imgkey in self.images.keys():
			if imgkey == tpl_imgkey or self.photometry_info[imgkey]['drop'] != 0 or self.photometry_info[imgkey]['flt']!=flt:
				continue
			photret_singlept = self.aperture_photometry_apphot_iraf_target_single(imgkey, which_dir = 'subtracted_image', centering = True, x=xc0, y=yc0)
			print imgkey, photret_singlept
			if len(photret_singlept) == 0:
				print "source not detected on image %s"%imgkey
				continue
			x = np.round(photret_singlept[0,0], 2)
			y = np.round(photret_singlept[0,1], 2)
			xs.append(x)
			ys.append(y)
			self.photometry_info[imgkey]['xsub'] = x
			self.photometry_info[imgkey]['ysub'] = y

		if showplot:
			fig,ax = plt.subplots(figsize=(6,6))
			plt.plot(xs, ys, 'ko')
			plt.show()

		xc = np.median(xs)
		yc = np.median(ys)

		return xc, yc


	def aperture_photometry_on_subtracted_images_flt(self, flt, tpl_imgkey, xc=None, yc=None, renewphot=False, centering=False, apsizefixed=False):
		'''
		See self.aperture_photometry_apphot_iraf_target_flt
		'''
		if xc is None or yc is None:
			xc,yc = self.image_subtraction_get_target_xy_on_subtracted_images(flt, tpl_imgkey, which_dir='raw_image', showplot=False)

		tpl_fwhm = self.photometry_info[tpl_imgkey]['fwhm']
		for image_key in self.images.keys():
			if image_key == tpl_imgkey or self.photometry_info[image_key]['drop'] !=0 or self.photometry_info[image_key]['flt']!=flt:
				continue
			if self.photometry_info[image_key]['instmag'] != self.photometry_info_init_values['instmag'] and (not renewphot):
				print "already done"
				continue

			subimg = self.__get_internal_image(image_key, which_dir='subtracted_image')
			hinfo = get_fits_info(subimg, 'CONVOL00')
			if hinfo['CONVOL00'] == 'IMAGE':
				imagefwhm = tpl_fwhm
			else:
				imagefwhm = self.photometry_info[image_key]['fwhm']

			photret_singlept = self.aperture_photometry_apphot_iraf_target_single(image_key, which_dir = 'subtracted_image', fwhm_img=imagefwhm, centering = centering, x=xc, y=yc)

			if len(photret_singlept) == 0:
				self.photometry_info[image_key]['drop'] = 4
				continue

			mag = photret_singlept[0,2]
			magerr = photret_singlept[0,3]

			self.photometry_info[image_key]['instmag'] = mag
			self.photometry_info[image_key]['instmagerr'] = magerr

	def psf_photometry_on_subtracted_images_flt(self, flt, tpl_imgkey, \
	xc=None, yc=None,  fitsky='yes', recenter='yes', ap_phot_recenter=False, renewphot=False, \
	renew_psfmodel=False, use_tplpst_prior = True,skiplist=None, which_dir='raw_image'):
		'''
		See self.iraf_psf_photometry_peak
		INPUTS:
		'''
		if xc is None or yc is None:
			xc,yc = self.image_subtraction_get_target_xy_on_subtracted_images(flt, tpl_imgkey, which_dir=which_dir, showplot=False)
			
		for imgkey in self.images.keys():
			if skiplist is not None:
				if imgkey in skiplist:
					continue
			if imgkey == tpl_imgkey or self.photometry_info[imgkey]['drop'] != 0 or self.photometry_info[imgkey]['flt']!=flt:
				continue
			if self.photometry_info[imgkey]['instmag'] != self.photometry_info_init_values['instmag'] and (not renewphot):
				continue
			subimg = self.__get_internal_image(imgkey, which_dir='subtracted_image')
			if not os.path.exists(subimg):
				print "No subtracted image for %s"%imgkey
				continue
			hinfo = get_fits_info(subimg, 'CONVOL00')
			if hinfo['CONVOL00'] == 'IMAGE':
				psffile = os.path.join(self.psfmodel_dir, tpl_imgkey.split('.')[0] + '_psf.fits')
				if not os.path.exists(psffile):
					self.iraf_get_psf_model_single_image(tpl_imgkey, image_dir=which_dir)
			else:
				psffile = os.path.join(self.psfmodel_dir, imgkey.split('.')[0] + '_psf.fits')
				if not os.path.exists(psffile) or renew_psfmodel:
					self.iraf_get_psf_model_single_image(imgkey, image_dir='interp_image', overwrite=True, reuse_pst=False, use_tplbase_as_prior=use_tplpst_prior, tplimgkey=tpl_imgkey, xytran_need=0)

			hinfo = get_fits_info(psffile, ['FITRAD','PSFRAD'])
			fitrad = hinfo['FITRAD']  #default self.fitrad_nfwhm* fwhm  where self.fitrad_nfwhm=1.3 by default
			psfrad = hinfo['PSFRAD']  #default self.psfrad_nfwhm* fwhm  where self.psfrad_nfwhm=4.0 by default
			print imgkey, fitrad, psfrad
			fwhm = psfrad/4.0

			imgkey_s = imgkey.split('.')[0]
			subimage = os.path.join(self.subtraction_dir, '%s_xys_model_sub.fits'%(imgkey_s))
			photfile = os.path.join(self.subtraction_dir, imgkey_s+'_xys.mag')
			peakfile = os.path.join(self.subtraction_dir, imgkey_s+'.peak')
			
			# fitsky, aperture radius=2*fwhm,
			apradius = int(fwhm*2.0)
			skyin = int(fwhm*2)
			skywidth = int(fwhm*2)
			fitsky1yes = 'yes'
			fitrad = int(psfrad/4.0*1.3)
			peakphot = self.iraf_psf_photometry_peak(imgkey, x=xc, y=yc, photfile=photfile, psffile=psffile,peakfile=peakfile, fitrad=fitrad, psfrad=psfrad, fitsky=fitsky1yes, recenter=recenter, ap_phot_recenter_xy= ap_phot_recenter, apradius=apradius, skyin=skyin, skywidth=skywidth, residual_image=True, subimage=subimage, residual_smallcut=True, update_phot_record=False, which_dir='subtracted_image')
			peakmsky = peakphot['MSKY'].data[0]
			peakmag = peakphot['MAG'].data[0]
			peakmagerr = peakphot['MERR'].data[0]

			apphot = Table.read(photfile, format='daophot')
			apmsky = apphot['MSKY'].data[0]
			apmag = apphot['MAG'].data[0]
			apmagerr = apphot['MERR'].data[0]			

			self.photometry_info[imgkey]['mskyfit1'] = peakmsky
			self.photometry_info[imgkey]['mskyap1'] = apmsky
			self.photometry_info[imgkey]['magap1'] = apmag
			self.photometry_info[imgkey]['magerrap1'] = apmagerr
			self.photometry_info[imgkey]['magfit1'] = peakmag
			self.photometry_info[imgkey]['magerrfit1'] = peakmagerr

			if fitsky == 'yes':
				self.photometry_info[imgkey]['instmag'] = peakmag
				self.photometry_info[imgkey]['instmagerr'] = peakmagerr

			
			# not fitsky, aperture radius=2*fwhm
			fitrad = int(psfrad/4.0*1.3)
			fitsky1no = 'no'
			peakphot = self.iraf_psf_photometry_peak(imgkey, x=xc, y=yc, photfile=photfile, psffile=psffile,peakfile=peakfile, fitrad=fitrad, psfrad=psfrad, fitsky=fitsky1no, recenter=recenter, ap_phot_recenter_xy= ap_phot_recenter, apradius=apradius, skyin=skyin, skywidth=skywidth, residual_image=True, subimage=subimage, residual_smallcut=True, update_phot_record=False, which_dir='subtracted_image')
			peakmsky = peakphot['MSKY'].data[0]
			peakmag = peakphot['MAG'].data[0]
			peakmagerr = peakphot['MERR'].data[0]
			
			self.photometry_info[imgkey]['mag1'] = peakmag
			self.photometry_info[imgkey]['magerr1'] = peakmagerr

			if fitsky == 'no':
				self.photometry_info[imgkey]['instmag'] = peakmag
				self.photometry_info[imgkey]['instmagerr'] = peakmagerr


			# fitsky, aperture radius 2*fwhm, 
			apradius = int(fwhm*2.0)
			skyin = int(fwhm*3)
			skywidth = int(fwhm*3)
			fitsky2yes = 'yes'
			fitrad = int(psfrad/4.0*1.0)
			peakphot = self.iraf_psf_photometry_peak(imgkey, x=xc, y=yc, photfile=photfile, psffile=psffile,peakfile=peakfile, fitrad=fitrad, psfrad=psfrad, fitsky=fitsky2yes, recenter=recenter, ap_phot_recenter_xy=ap_phot_recenter, apradius=apradius, skyin=skyin, skywidth=skywidth, residual_image=True, subimage=subimage, residual_smallcut=True, update_phot_record=False, which_dir='subtracted_image')
			peakmsky = peakphot['MSKY'].data[0]
			peakmag = peakphot['MAG'].data[0]
			peakmagerr = peakphot['MERR'].data[0]

			apphot = Table.read(photfile, format='daophot')
			apmsky = apphot['MSKY'].data[0]
			apmag = apphot['MAG'].data[0]
			apmagerr = apphot['MERR'].data[0]			

			self.photometry_info[imgkey]['mskyfit2'] = peakmsky
			self.photometry_info[imgkey]['mskyap2'] = apmsky
			self.photometry_info[imgkey]['magap2'] = apmag
			self.photometry_info[imgkey]['magerrap2'] = apmagerr
			self.photometry_info[imgkey]['magfit2'] = peakmag
			self.photometry_info[imgkey]['magerrfit2'] = peakmagerr


			# not fitsky
			fitsky2no = 'no'
			fitrad = int(psfrad/4.0*1.3)
			peakphot = self.iraf_psf_photometry_peak(imgkey, x=xc, y=yc, photfile=photfile, psffile=psffile,peakfile=peakfile, fitrad=fitrad, psfrad=psfrad, fitsky=fitsky2no, recenter=recenter, ap_phot_recenter_xy=ap_phot_recenter, apradius=apradius, skyin=skyin, skywidth=skywidth, residual_image=True, subimage=subimage, residual_smallcut=True, update_phot_record=False, which_dir='subtracted_image')
			peakmsky = peakphot['MSKY'].data[0]
			peakmag = peakphot['MAG'].data[0]
			peakmagerr = peakphot['MERR'].data[0]
			
			self.photometry_info[imgkey]['mag2'] = peakmag
			self.photometry_info[imgkey]['magerr2'] = peakmagerr


			#self.photometry_info[imgkey]['fitsky'] = fitsky
			#self.photometry_info[imgkey]['magdiff'] = peakmag-apmag


	def image_subtraction_residual_check_with_aperture_photometry(self, imgkey, phot_outfile=None):
		'''
		This is to check how the residual behave for all stars 
		'''
		
		xyfile = os.path.join(self.subtraction_dir, imgkey.split('.')[0]+'.coo')
		self.__source_detection_single_apphot_daofind(imgkey, xyfile, which_dir='interp', threshold=5, sigma=None)	
		xymags = Table.read(xyfile, format='daophot') #XCENTER     YCENTER      MAG     SHARPNESS      SROUND       GROUND      ID 
		xys = np.array([xymags['XCENTER'].data, xymags['YCENTER'].data]).transpose()

		photfile  = os.path.join(self.subtraction_dir, imgkey.split('.')[0]+'.apphot')	
		
		image = self.__get_internal_image(imgkey, which_dir='sub')
		self.get_apphot_iraf_parameters(imgkey, saveparfile=True)
		options = self.apphot_iraf_options.copy()
		if phot_outfile is None:
			photfile = os.path.join(self.aperture_photometry_dir, imgkey.split('.')[0]+'.mag')
		else:
			photfile = phot_outfile
		self.__delete_file_if_exist(photfile)
		xymags2 = self.aperture_photometry_apphot_iraf_single_image_xys_given(image, xys, options, output=photfile, centering=False) #xymags2 only have x,y,mag,magerr
		photret = Table.read(photfile, format='daophot')		

		fig, ax = plt.subplots(figsize=(9,6))
		mask = photret['FLUX']!=0
		ax.plot(xymags['MAG'][mask], photret['FLUX'][mask], 'ko')
		ax.hlines(0, xmin=np.min(xymags['MAG']), xmax=np.max(xymags['MAG']), colors='r')
		plt.show()

		#highlight those with high residual counts
		countcut = raw_input("cut on count:")
		magcut = raw_input("cut on mag:")
		if countcut  == '' and magcut == '':
			sources = photret
		elif countcut == '':
			countcut = float(countcut)
			sources = photret[np.absolute(photret['FLUX'])>countcut]
		elif magcut == '':
			magcut = float(magcut)
			sources = photret[np.absolute(xymags['MAG'])>magcut]
		else:
			magcut = float(magcut)
			countcut = float(countcut)
			sources = photret[(np.absolute(photret['FLUX'])>countcut)*(np.absolute(xymags['MAG'])>magcut)]

		self.display_sources_with_ds9(imgkey, sources, xcol='XCENTER', ycol='YCENTER', radius=10, color='green', width=1, which_dir='sub', newds9=True)




	def image_subtraction_psfphotret_diagnosis(self, flt, maxmagerr=0.5, maxstd=0.5, docolnames=None, skipcolnames=None):
		'''
		Check the photometry results on subtracted images with different photometry configurations 
		'''
		self.__dict2table()
		fltdata = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		fltdata = fltdata[fltdata['drop']==0]

		magcolnames = ['magap1','magfit1','mag1','magap2', 'magfit2', 'mag2']	
		magerrcolnames =['magerrap1','magerrfit1','magerr1','magerrap2', 'magerrfit2', 'magerr2']
		fig, [ax1, ax2, ax3, ax4] = plt.subplots(1,4,figsize=(18,12), sharex=True, gridspec_kw = {'left':0.05, 'right':0.95, 'bottom':0.5})
		ax5 = fig.add_axes([0.05,0.05,0.9,0.4])
		markers = ['>','s','s', '*','o','o']
		linetypes = {'1':'-', '2':'-.'}
		colors = ['k','b','g','k','m','c']

		for i, (magkey,magerrkey) in enumerate(zip(magcolnames, magerrcolnames)):
			if docolnames is not None:
				if magkey not in docolnames:
					continue
			if skipcolnames is not None:
				if magkey in skipcolnames:
					continue
			vnum = magkey[-1]
			ax1.errorbar(fltdata['obstime'], fltdata[magkey], yerr=fltdata[magerrkey], fmt=colors[i]+markers[i], alpha=0.5, label=magkey)
		ax1.errorbar(fltdata['obstime'], fltdata['instmag'], yerr=fltdata['instmagerr'], fmt='ro', label='adopted')
			
		for i, (magkey,magerrkey) in enumerate(zip(magcolnames, magerrcolnames)):
			if magkey in ['mag1', 'mag2']:
				continue
			mskykey = magkey.replace('mag', 'msky')
			vnum = magkey[-1]
			ax2.plot(fltdata['obstime'], fltdata[mskykey], linetypes[vnum], color=colors[i], label=mskykey)

		magcolnames_psf    =['magfit1',   'mag1',   'magfit2', 'mag2']	
		magerrcolnames_psf =['magerrfit1','magerr1','magerrfit2', 'magerr2']

		mags = fltdata[magcolnames_psf]
		magerrs = fltdata[magerrcolnames_psf]
		indexs  = list(range(len(mags)))
		mags_rowmean = []
		magerrs_rowstd = []
		mags_diff = []
		for i, (magrow, magerrrow) in enumerate(zip(mags, magerrs)):
			wholerow = fltdata[i]
			imgkey = wholerow['name']
			magerrdata = np.array([magerrrow[colname] for colname in magerrrow.colnames])
			magdata = np.array([magrow[colname] for colname in magrow.colnames])
			if np.any(magerrdata>maxmagerr) or np.any(magdata==0) or np.std(magdata)>0.3:
				print "romove row %s %s"%(str(i), imgkey)
				indexs.remove(i)
				continue
			psfmagmean = np.mean(magdata)
			psfmagstd  = np.std(magdata)
			mags_rowmean.append(psfmagmean)
			magerrs_rowstd.append(psfmagstd)

			instmag = wholerow['instmag']
			magdiff = instmag -psfmagmean	
			self.photometry_info[imgkey]['magdiff'] = magdiff
			mags_diff.append(magdiff)

		ax5.set_xlim([0,100])
		ax5.set_ylim([0,10])

		ax3.errorbar(fltdata['obstime'][indexs], mags_rowmean, yerr=magerrs_rowstd, fmt='ko', label='psfmag_mean')
		ax3.errorbar(fltdata['obstime'][indexs], fltdata['instmag'][indexs], yerr=fltdata['instmagerr'][indexs], fmt='ro', label='instmag')
		for i, (imgkey, obstime, instmag) in enumerate(zip(fltdata['name'], fltdata['obstime'], fltdata['instmag'])):
			ax3.plot([obstime, obstime],[instmag-0.4, instmag+0.4], 'k', alpha=0.3)
			if np.mod(i,2) == 0:	
				ax3.text(obstime, instmag+0.5, str(i), horizontalalignment='center')
			else:
				ax3.text(obstime, instmag-0.5, str(i), horizontalalignment='center')
	
			irow = i/10
			icol = np.mod(i,10)
			ax5.text(5+10*icol, 10-irow-1, '%s:%s'%(i, imgkey))

		

		ax4.plot(fltdata['obstime'][indexs], magerrs_rowstd, '-ko', label='psfmag_std')
		#ax4.plot(fltdata['obstime'][indexs], magerrs_rowstd, 'k')
		ax4.plot(fltdata['obstime'][indexs], mags_diff, '-ro', label='instmag - psfmag_mean')
		#ax4.plot(fltdata['obstime'][indexs], mags_diff, 'r')
		for i,magcolname in enumerate(magcolnames):
			ax4.plot(fltdata['obstime'][indexs], fltdata['instmag'][indexs]-fltdata[magcolname][indexs], '-'+markers[i], color=colors[i], alpha=0.5, label='instmag - %s'%magcolname)

		ax1.invert_yaxis()
		ax3.invert_yaxis()
		ax1.legend(loc=0)
		ax2.legend(loc=0)
		ax3.legend(loc=0)
		ax4.legend(loc=0)



		plt.savefig(os.path.join(self.result_dir, '%s_image_subtraction_psfphotret_diagnosis.png'%flt))
		plt.show()
		
	def image_subtraction_psfphot_renew(self, flt, magcolname, whichones='all'):
		'''
		renew instmag and instmagerr with alternatives: mag1, mag2, magap1, magap2, magfit1, magfit2
		'''
		if isinstance(whichones, str):
			if whichones == 'all':
				whichones = [imgkey for imgkey in self.images.keys() if self.photometry_info[imgkey]['flt']==flt]
				note_imgkeys = 'ALL'
			else:
				whichones = [whichones]
				note_imgkeys = whichones
		else:
			note_imgkeys = ' '.join(whichones)

		self.add_note("instmag(s) for image(s) [%s] are replaced by %s"%(note_imgkeys, magcolname))

		magerrcolname = magcolname.replace('mag', 'magerr')
		for imgkey in whichones:
			self.photometry_info[imgkey]['instmag'] = self.photometry_info[imgkey][magcolname]
			self.photometry_info[imgkey]['instmagerr'] = self.photometry_info[imgkey][magerrcolname]
		
		

	def hotpants_image_subtraction_single(self, input_imgkey, tpl_imgkey=None, tplimage=None, \
	fwhm_tplimg = None, fwhm_interpimg = None, inputimg_mask=None, refimg_mask=None, verbose=1):
		'''
		INPUTS:
			See self.__hotpants_image_subtraction for details
		'''
		in_imgkey_s = input_imgkey.split('.')[0]
		input_image = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)

		if tpl_imgkey is not None:
			tplimage = os.path.join(self.template_dir, 'tplsub_'+tpl_imgkey)
			fwhm_tplimg = self.photometry_info[tpl_imgkey]['fwhm']
			geotranfac = self.__get_image_geometric_transformation_scale_factor(input_imgkey, tpl_imgkey)
		elif tplimage is not None:
			if fwhm_tplimg is None:
				raise ValueError("fwhm of the template is required")
		else:
			raise ValueError("please provide valid template image")

		output_image = os.path.join(self.subtraction_dir, 'sub_'+input_imgkey)
		convolved_image = os.path.join(self.subtraction_dir, 'conv_'+input_imgkey)

		fwhm_inimg  = self.photometry_info[input_imgkey]['fwhm']
		if fwhm_interpimg is None:
			fwhm_interpimg = fwhm_inimg*geotranfac

		if verbose:
			print "The input image pre-sampling FWHM: %s"%fwhm_inimg
			print "The input image after-resampling FWHM: %s"%fwhm_interpimg
			print "The reference image FWHM: %s"%fwhm_tplimg

		if fwhm_inimg <= 0 or fwhm_inimg ==99.99:
			raise ValueError("non-realistic FWHM value for the input image...")

		fwhm = int(np.max([fwhm_interpimg, fwhm_tplimg]))

		telescope_gain = self.base.telescopes_info[self.current_telescope]['gain']
		telescope_rdnoise = self.base.telescopes_info[self.current_telescope]['readnoise']
		tgain = telescope_gain
		igain = telescope_gain
		trdnoise = telescope_rdnoise
		irdnoise = telescope_rdnoise

		sigdiff = np.sqrt(np.abs(fwhm_tplimg**2 - fwhm_interpimg**2))
		sigma1 = np.round(0.5*sigdiff, 1)
		sigma2 = np.round(1*sigdiff, 1)
		sigma3 = np.round(2.0*sigdiff, 1)
		#ATTENTIONS!!!
		#hotpants will very likely go wrong if the smallest gaussian width < 0.3
		#it's embarrassing to say that I DON'T know why...
		#check the sigma1, and use np.max([sigma1, 0.3]) as the first gaaussian width
		if sigma1 < 0.3:
			print "!!!! The calculated gaussian width from the difference of the FWHM of two images are too small... "
			sigma1 = 0.3
			sigma2 = 0.6
			sigma3 = 1.2

		def fill_hotpants_par(parname, parvalue):
			if parname not in self.hotpants_par_hold_list:
				pardict[parname] = parvalue
			else:
				print "You have hold %s fixed"%parname

		pardict = self.hotpants_pars.copy()

		if fwhm_tplimg > fwhm_interpimg:
			fill_hotpants_par('-c', 'i')
		else:
			fill_hotpants_par('-c', 't')

		fill_hotpants_par('-inim', input_image)
		fill_hotpants_par('-tmplim', tplimage)
		fill_hotpants_par('-outim', output_image)
		fill_hotpants_par('-oci', convolved_image)
		if inputimg_mask is not None:
			fill_hotpants_par('-imi', inputimg_mask)
		if refimg_mask is not None:
			fill_hotpants_par('-tmi', refimg_mask)
		fill_hotpants_par('-tg',tgain)
		fill_hotpants_par('-ig',igain)
		fill_hotpants_par('-tr',trdnoise)
		fill_hotpants_par('-ir',irdnoise)
		fill_hotpants_par('-r',2*fwhm)
		fill_hotpants_par('-kcs', 4*fwhm+1)
		fill_hotpants_par('-rss',3*fwhm)
		fill_hotpants_par('-ng', '3 6 %s 4 %s 2 %s'%(sigma1, sigma2, sigma3))

		#if fwhm_tplimg > fwhm_interpimg:
		#	pardict['-c'] = 'i'
		#else:
		#	pardict['-c'] = 't'

		#pardict['-inim'] = input_image
		#pardict['-tmplim'] = tplimage
		#pardict['-outim'] = output_image
		#pardict['-oci']  = convolved_image
		#if inputimg_mask is not None:
		#	pardict['-imi'] = inputimg_mask
		#if refimg_mask is not None:
		#	pardict['-tmi'] = refimg_mask
		#pardict['-tg']=tgain
		#pardict['-ig']=igain
		#pardict['-tr']=trdnoise
		#pardict['-ir']=irdnoise
		#pardict['-r'] = 2*fwhm
		#pardict['-kcs'] = 4*fwhm+1
		#pardict['-rss'] = 3*fwhm
		#pardict['-ng']= '3 6 %s 4 %s 2 %s'%(sigma1, sigma2, sigma3)

		self.__hotpants_image_subtraction(pardict=pardict)

	def __hotpants_image_subtraction(self, pardict=None, verbose=1):
		'''
		Execute hotpants with options from pardict if provided otherwise from self.hotpants_pars
		Check self.__init_hotpants_options for the detailed options.
		INPUTS:
		'''
		if pardict is None:
			pardict = self.hotpants_pars.copy()
		if pardict['-inim'] is None or pardict['-tmplim'] is None or pardict['-outim'] is None:
			raise ValueError('input image, template image and output image are required')

		options = ''
		for kw in pardict.keys():
			if pardict[kw] is not None:
				if isinstance(pardict[kw], bool):
					if pardict[kw]:
						options = options + ' ' + kw
				else:
					options = options + ' ' +kw + ' ' + str(pardict[kw])

		hotpants = os.path.join(self.base.hotpants_dir, 'hotpants')
		hotpants_command = hotpants + ' ' + options
		if verbose:
			print hotpants_command

		try:
			os.system(hotpants_command)
		except:
			print "hotpants fails..."

	def get_image_geometric_transformation_scale_factor_flt(self, flt, tplimgkey, dolist=None, skiplist=None):
		'''
		To check the plate scale difference with the template image
		'''
		for imgkey in self.images.keys():
			if dolist is not None:
				if imgkey not in dolist:
					continue
			if skiplist is not None:
				if imgkey in skiplist:
					continue
			if self.photometry_info[imgkey]['flt'] != flt:
				continue
			print imgkey
			transcoeff_file = os.path.join(self.stars_dir, "match_%s_%s.coeff"%(imgkey.split('.')[0], tplimgkey.split('.')[0]))
			if not os.path.exists(transcoeff_file):
				print "No transformation coefficient file for %s"%imgkey
				continue
			porder, dxfit_values, dyfit_values = self.__extract_transformation_fitting(transcoeff_file)
			self.photometry_info[imgkey]['geotranfac'] =  dxfit_values[1]


	def __get_image_geometric_transformation_scale_factor(self, imgkey, tplimgkey):
		transcoeff_file = os.path.join(self.stars_dir, "match_%s_%s.coeff"%(imgkey.split('.')[0], tplimgkey.split('.')[0]))
		porder, dxfit_values, dyfit_values = self.__extract_transformation_fitting(transcoeff_file)
		geotranfac = np.abs(1.0/dxfit_values[1])
		self.photometry_info[imgkey]['geotranfac'] =  geotranfac
		return geotranfac

	def get_convolution_direction_flt(self, flt):
		'''
		check the convolution direction, convolve on TEMPLATE or IMAGE
		'''
		for imgkey in self.images.keys():
			if self.photometry_info[imgkey]['flt'] != flt:
				continue
			subimage = self.__get_internal_image(imgkey, which_dir='sub')
			if not os.path.exists(subimage):
				continue
			print imgkey
			CONVOL00 = hinfo = get_fits_info(subimage, 'CONVOL00')['CONVOL00']
			self.photometry_info[imgkey]['convon'] = CONVOL00

#Image alignment and resampling
	def prepare_template_image(self, tplmethod=1, tplflt=None, tplimgkey=None, tplimg_bricks=None, reference_imgkey=None,  which_dir='raw_image'):
		'''
		This function is used to prepare template image for image subtraction. Different procedures are required in different situations.
			(1) use one single image for template image
			(2) co-adding multiple images, which is ok when the input images are well aligned
			(3) align and resample multiple images and co-add them to build the template image, which is ok when the images are taken under stable conditions consecutively
			(4) convolution with kernel is needed before co-adding in 2)


		INPUTS:
			tplmethod: 1, 2, 3, 4 corresponding to procedures required above
			tplimgkey: the selected template image when tplmethod is 1
			tplimg_bricks: the raw materials to build the template image, list of names of key for the input images for example ['001.fits','002.fits','003.fits']
			reference_img: the reference image for alignment, resampling and even convolution
			which_dir: 'raw_image', 'modified_image', 'subtraction'
		'''

		if tplmethod == 1:
			self.__assign_single_image_as_template(tplimgkey, tplflt, which_dir=which_dir)
		elif tplmethod == 2:
			if which_dir == 'raw_image':
				imgdir = self.raw_image_dir
			elif which_dir == 'modified_image':
				imgdir = self.modified_image_dir
			elif which_dir == 'subtraction':
				imgdir = self.subtraction_dir
			else:
				raise ValueError("input for which_dir is not recognised...")

			imgs_input = open('tplimg_build_bricks_temp.txt', 'awt')
			for img in tplimg_bricks:
				imginput = os.path.join(imgdir, img)
				imgs_input.write(imginput+'\n')

			imgs_input.close()
			if tplflt is None:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_temp.fits')
			else:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_%s.fits'%tplflt)

			imcombine_iraf("@%s"%imgs_input, tplimg_output)
			self.__delete_file_if_exist(imgs_input)
		elif tplmethod == 3:
			if tplflt is None:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_temp.fits')
			else:
				tplimg_output = os.path.join(self.template_dir, 'tplsub_%s.fits'%tplflt)


			self.align_and_stack_multiple_images(tplimg_bricks, reference_imgkey, tplimg_output, input_which_dir=which_dir)

		elif tplmethod == 4:
			print "under construction"

	def __assign_single_image_as_template(self, imgkey, flt, which_dir= 'raw_image'):
		'''
		just copy the given image to the destination directory
		'''
		self.__dict2table()
		info_sf = self.__select_rows_from_table(self.photometry_record_table,'flt',flt,)
		names = info_sf['name']
		if imgkey not in names:
			raise ValueError("%s not in %s band"%(imgkey,flt))
		tplimage = self.__get_internal_image(imgkey, which_dir=which_dir)
		dstimg = os.path.join(self.template_dir, "tplsub_"+imgkey)
		shutil.copy(tplimage, dstimg)

		for img in names:
			self.photometry_info[img]['template'] = 0
		self.photometry_info[imgkey]['template'] = 1

	def align_and_stack_multiple_images(self, imgkeylist, refimgkey, output_image, input_which_dir='raw_image', renew_matched_sourcelist = False, renew_resampled_img=False):
		'''
		This function is to align and co-add multiple images.

		How?
			(1) align and resample the input images to the frame of reference image, see self.align_and_resample_image_given_matched_source_list
			(2) combine reference image and resampled input images, see imcombine_iraf

		INPUTS:
			imgkeylist: input images list without the reference image
			refimgkey: the reference image which will be one of the images to be combined
		'''

		tplimg_listfile = 'tplimg_build_bricks_temp.txt'
		self.__delete_file_if_exist(tplimg_listfile)

		imgs_input = open(tplimg_listfile, 'awt')

		refimage = self.__get_internal_image(refimgkey, which_dir=input_which_dir)
		imgs_input.write(refimage + '\n')

		ref_imgkey_s  = refimgkey.split('.')[0]

		for imgkey in imgkeylist:
			print imgkey
			if imgkey == refimgkey:
				continue

			in_imgkey_s = imgkey.split('.')[0]

			matched_xyxy_outfile = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(in_imgkey_s, ref_imgkey_s))
			if (not os.path.exists(matched_xyxy_outfile)) or renew_matched_sourcelist:
				self.match_sources_on_two_images(imgkey, refimgkey, option=None, mode='file', ref_col1=1, ref_col2=2, input_col1=1, input_col2=2, input_order_col=-3, ref_order_col=-3, max_distance=1, matched_xyxy_outfile=None, match_result_display=False, show_match_stat=False)

			resample_outfile = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)
			if (not os.path.exists(resample_outfile)) or renew_resampled_img:
				outimg = self.align_and_resample_image_given_matched_source_list(imgkey, refimgkey, which_dir=input_which_dir, verbose=0)
			imgs_input.write(resample_outfile+'\n')

		imgs_input.close()
		input_imglist_file = "@%s"%tplimg_listfile
		imcombine_iraf(input_imglist_file, output_image)
		self.__delete_file_if_exist(tplimg_listfile)

	def align_and_resample_image_given_transfile(self, input_imgkey, ref_imgkey, trans_file=None, refimg_dir='raw_image', inputimg_dir='raw_image'):
		'''
		This function is to align the input image with the reference image and resample the input image to the same image coordinate of the reference image. The transformation file is used.
		INPUTS:
			input_imgkey:
			ref_imgkey:
			trans_file: the human readable transformation file with key=value entries
		'''
		in_imgkey_s = input_imgkey.split('.')[0]
		ref_imgkey_s = ref_imgkey.split('.')[0]

		if trans_file is None:
			trans_file = os.path.join(self.stars_dir, "match_%s_%s.coeff"%(in_imgkey_s, ref_imgkey_s))
			if not os.path.exists(trans_file):
				raise ValueError("The required transformation file not available. You can prepare it with self.match_sources_on_two_images")
		output_image = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)
		input_image = self.__get_internal_image(input_imgkey, which_dir=inputimg_dir)

		ref_image = self.__get_internal_image(ref_imgkey, which_dir=refimg_dir)
		Nx_ref,Ny_ref = self.__get_fits_image_size(ref_image)

		self.__fitsh_fitrans(input_image, output_image, trans_file, sx=Nx_ref, sy=Ny_ref)

	def align_and_resample_image_given_matched_source_list(self, input_imgkey, ref_imgkey, matched_source_file = None, which_dir='raw_image',  verbose=1):
		'''
		This function is to align the input image with the reference image and resample the input image to the same image coordinate of the reference image. The matched source list should be built in advance and given as input
		INPUTS:
			input_imgkey:
			ref_imgkey:
			matched_source_file: the file containing matched sources with the format of x,y,x,y
		'''
		in_imgkey_s = input_imgkey.split('.')[0]
		ref_imgkey_s = ref_imgkey.split('.')[0]

		if matched_source_file is None:
			matched_source_file = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(in_imgkey_s, ref_imgkey_s))
		if not os.path.exists(matched_source_file):
			raise IOError("matched source list file %s not exists"%matched_source_file)
		xygrid_para = os.path.join(self.parafile_dir, 'xygrid.par')
		xygrid_outfile = os.path.join(self.stars_dir, 'xygrid_%s_%s.txt'%(in_imgkey_s, ref_imgkey_s))
		xygrid = os.path.join(self.base.diapl_dir, 'xygrid')
		xygrid_command = "%s %s %s %s"%(xygrid, xygrid_para, matched_source_file, xygrid_outfile)
		try:
			if verbose:
				print xygrid_command

			os.system(xygrid_command)
		except:
			print "xygrid failure for image %s ..."%input_imgkey

		resample_para = os.path.join(self.parafile_dir, 'resample2.par')
		instrument_para = os.path.join(self.parafile_dir, 'instrument.par')
		resample_outfile = os.path.join(self.subtraction_dir, 'interp_%s.fits'%in_imgkey_s)
		input_image = self.__get_internal_image(input_imgkey, which_dir=which_dir)
		resample2 = os.path.join(self.base.diapl_dir, 'resample2')
		resample_command  = "%s %s %s %s %s %s"%(resample2, resample_para, instrument_para, xygrid_outfile, input_image, resample_outfile)

		try:
			if verbose:
				print resample_command
			os.system(resample_command)
		except:
			print "resample failure for image %s ..."%input_imgkey

		return resample_outfile

	def interpolate_resample_to_the_template_image(self, imgkey,  matched_source_listfile=None, tpl_imgkey=None, which_dir_input = 'modified_image', verbose=True):
		'''
		Interpolate and resample the input image to the frame of reference image. The matched souce list matched_source_listfile or the reference image tpl_imgkey is required

		INPUTS:
			imgkey:
			matched_source_listfile:
			tpl_imgkey:
			which_dir_input:
		'''
		imgkey_s = imgkey.split('.')[0]

		if matched_source_listfile is None:
			if tpl_imgkey is None:
				raise ValueError("the matched source list file between input image and reference image is not provided...")
			else:
				ref_imgkey_s = tpl_imgkey.split('.')[0]
				matched_source_listfile = os.path.join(self.stars_dir, "match_%s_%s_xyxy.txt"%(imgkey_s, ref_imgkey_s))

		xygrid_parafile = os.path.join(self.parafile_dir, 'xygrid.par')
		xygrid_outfile = os.path.join(self.stars_dir, "%s_xygrid.out"%imgkey_s)
		xygrid_command = "xygrid %s %s %s"%(xygrid_parafile, matched_source_listfile, xygrid_outfile)
		if os.path.exists(xygrid_outfile):
			os.remove(xygrid_outfile)

		if verbose:
			print xygrid_command
		os.system(xygrid_command)

		if not os.path.exists(xygrid_outfile):
			raise IOError("xygrid didn't get the expected output...")
		else:
			resample2_infile = xygrid_outfile

		resample2_parafile = os.path.join(self.parafile_dir, 'resample2.par')
		instrument_parafile = os.path.join(self.parafile_dir, 'instrument.par')
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		interp_image = os.path.join(self.subtraction_dir, 'interp_'+imgkey)
		if os.path.exists(interp_image):
			os.remove(interp_image)

		resample2 = os.path.join(self.base.diapl_dir, 'resample2')
		resample2_command = "%s %s %s %s %s %s"%(resample2, resample2_parafile, instrument_parafile, resample2_infile, input_image, interp_image)
		if verbose:
			print resample2_command
		os.system(resample2_command)

#Image manipulation
	def add_constant_to_image(self, imgkey, addvalue, input_which_dir='modified_image', output_image=None):
		'''
		add constant value addvalue to the inputimage
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=input_which_dir)
		if output_image is None:
			output_image = input_image

		imarith_iraf(input_image, '+', addvalue, output_image)

	def __subtract_images(self, img1, img2):
		'''
		subtract img1 by img2 and save the output image in self.modified_image_dir

		'''
		img_sci_abs = os.path.join(self.raw_image_dir, img1)
		img_sky_abs = os.path.join(self.raw_image_dir, img2)
		img_ret_abs = os.path.join(self.modified_image_dir, img1)
		self.__delete_file_if_exist(img_ret_abs)

		imarith_iraf(img_sci_abs, '-', img_sky_abs, img_ret_abs)

	def flip_image_data(self, imgkey, invertwhich, which_dir='raw_image', output_image=None):
		'''
		flip array in the left/right direction or up/down direction
		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		hdu_input = fits.open(input_image)
		data_input = hdu_input[0].data

		if invertwhich == 'x':
			data_new = np.fliplr(data_input)
		elif invertwhich == 'y':
			data_new = np.flipud(data_input)
		else:
			raise ValueError('invalid input for invertwhich')

		hdu_output = fits.PrimaryHDU(data=data_new)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir, imgkey)

		if os.path.exists(output_image):
			os.remove(output_image)

		hdu_output.writeto(output_image)

	def rotate_image_data(self, imgkey, rotdeg=90, which_dir = 'raw_image', output_image = None):
		'''

		INPUTS:
			imgkey:
			rotdeg: 90, 180, 270 clockwise direction

		'''
		input_image = self.__get_internal_image(imgkey, which_dir=which_dir)
		hdu_input = fits.open(input_image)
		data_input = hdu_input[0].data
		rotnum = rotdeg/90
		data_new = np.rot90(data_input, k=rotnum)

		hdu_output = fits.PrimaryHDU(data=data_new)
		if output_image is None:
			output_image = os.path.join(self.modified_image_dir, imgkey)

		if os.path.exists(output_image):
			os.remove(output_image)

		hdu_output.writeto(output_image)

	def image_cutout(self, xl=None, yl=None, xh=None, yh=None, edge_width=None, \
	input_image_dir='raw_image', output_image_dir='modified_image', skiplist=None, dolist=None):
		'''
		cut out the wanted portion
		'''
		for imgkey in self.images.keys():
			if dolist is not None:
				if imgkey not in dolist:
					continue
			if skiplist is not None:
				if imgkey in skiplist:
					continue
			self.image_cutout_single(imgkey, xl=xl,yl=yl,xh=xh,yh=yh, edge_width=edge_width, input_image_dir=input_image_dir, output_image_dir=output_image_dir)

	def image_cutout_single(self, imgkey, xl=None, yl=None, xh=None, yh=None, edge_width=None, \
	input_image_dir='raw_image', output_image_dir='modified_image'):
		'''
		image cut-out
		First, edge_width set the first cutout section; then xl, xh, yl, yh can individually override
		'''
		input_image  = self.__get_internal_image(imgkey, which_dir=input_image_dir)
		output_image = self.__get_internal_image(imgkey, which_dir=output_image_dir)

		temp_image = 'cuttemp.fits'
		if os.path.exists(temp_image):
			os.remove(temp_image)

		NX, NY = self.__get_fits_image_size(input_image)

		x1 = 1
		x2 = NX
		y1 = 1
		y2 = NY

		if edge_width is not None:
			x1 = edge_width+1
			y1 = edge_width+1
			x2 = NX - edge_width
			y2 = NY - edge_width
		if xl is not None:
			x1 = xl
		if yl is not None:
			y1 = yl
		if xh is not None:
			x2 = xh
		if yh is not None:
			y2 = yh

		self.__cut_image_imcopy(input_image, temp_image, x1, x2, y1, y2)
		shutil.move(temp_image, output_image)

	def __cut_image_imcopy(self, input_image, output_image, x1,x2,y1,y2):
		'''
		input_image[x1:x2,y1:y2] --> output_image
		'''
		input_image_section = input_image + '[%s:%s,%s:%s]'%(x1,x2,y1,y2)
		imcopy_iraf(input_image_section, output_image)


	def __cut_image_diapl(self,in_img_name,cutted_img_name,x0,y0,xe,ye):
		'''
		cut the rectangle region defined by low left corner(x0,y0) and upper right corner (xe,ye) of the 'in_img_name'
		and save as 'out_img_name'
		'''
		cut = os.path.join(self.base.diapl_dir, 'cutfitsim')
		cut_command = "%s %s %s %s %s %s %s" % (cut,in_img_name,cutted_img_name,x0,xe,y0,ye)

		try:
			os.system(cut_command)
		except:
			print "cut image failure on %s"%in_img_name



	def __get_bad_pixels(self,image,threshold, mode='saturation'):
		'''
		bad pixels have value higher than saturation threshold in saturation mode or have negative values in negative mode
		INPUTS:
			image:
			threshold:
			mode: saturation or negative
		'''
		mask = np.zeros(image.shape)
		if mode == 'saturation':
			mask[np.where(image > threshold)]=1
		elif mode == 'negative':
			mask[np.where(image<threshold)] =1
		else:
			raise ValueError("invalid input for mode...")

		return mask

	def __fill_bad_pixels(self,image,bad_pixels,fill_value=None, mode='bkg', verbose=0):
		'''
		fill the bad pixels with the desired value which is determined by option 'mode'
		mode = 'bkg' or 'vicinity'
		when mode = 'bkg' the bad_pix will be filled by the value of the image background
		when mode = 'vicinity' the bad_pix will be filled by the average value ofthe surrounding pixels(time consuming!!)
		when mode = 'zero' the bad value will be replaced with zeros
		when mode = 'interp' the bad values will be replaced with the interpolated value at that coordinate
		'''
		M,N = image.shape
		print M,N
		X,Y = np.where(bad_pixels ==1)

		if fill_value is not None:
			image[X,Y] = fill_value
		else:
			if mode=='vicinity':
				shifts = self.__get_surrounding(radius = 1)
				#print shifts
				temp = np.ma.zeros((len(shifts),M,N))

				for i,[x,y] in enumerate(shifts):
				    xlo = 0+x
				    ylo = 0+y
				    xho = M+x
				    yho = N+y
				    img_shift = self.__cut_or_extend_imagedata(image,xlo,ylo,xho,yho)
				    badpix_shift = self.__cut_or_extend_imagedata(bad_pixels,xlo,ylo,xho,yho)
				    masked_image = np.ma.array(img_shift,mask = badpix_shift)
				    #print type(masked_image)
				    temp[i] = masked_image
				#print type(temp)
				mean_image=np.mean(temp,axis = 0)
				for x,y in zip(X,Y):
				    #print image[x,y],mean_image[x,y]
				    image[x,y]=mean_image[x,y]
				    #print i
			elif mode=='bkg':
			    	bkg,std =self.__get_rough_background(image)
			        image[X,Y]=bkg
			elif mode == 'zero':
			    	image[X,Y] = 0
			elif mode == 'interp':
				fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))
				ax1.imshow(image, origin='lower', cmap='gray_r', interpolation='none', vmin=np.min(image),vmax=np.median(image))
				ax2.imshow(bad_pixels, origin='lower', cmap='gray_r', interpolation='none')
				plt.show()
				xs,ys = np.where(bad_pixels==0)	
				print xs, ys
				print M,N
				fig,ax = plt.subplots(figsize=(6,6))
				img_ma = np.ma.masked_array(image, bad_pixels)
				ax.imshow(img_ma, origin='lower', cmap='gray_r', interpolation='none', vmin=np.min(img_ma),vmax=np.median(img_ma))
				plt.show()
				
				#f = interp2d(xs, ys, image[xs,ys], kind='linear')
				#img2d_interp = f(np.arange(M),np.arange(N)).transpose()
				
				#grid_x, grid_y = np.mgrid[0:M,0:N]
				#img2d_interp = griddata(np.array([xs,ys]).reshape(len(xs),2), image[xs,ys], (grid_x, grid_y), method='linear')
				
				img2d_interp = image.copy()
				for i, (rowdata, rowmask) in enumerate(zip(image, bad_pixels)):
					if np.sum(rowmask)!=0:
						xfit = np.arange(len(rowdata))
						inds = np.where(rowmask==0)
						f = interp1d(xfit[inds], rowdata[inds], kind='linear')
						value_interp = f(xfit[np.where(rowmask==1)])
						img2d_interp[i, np.where(rowmask==1)] = value_interp
				
				fig,ax = plt.subplots(figsize=(6,6))
				ax.imshow(img2d_interp, origin='lower', cmap='gray_r', interpolation='none')
				plt.show()
				print img2d_interp.shape
				print X,Y
				values_fill = img2d_interp[X,Y]
				image[X,Y] = values_fill
				if verbose:
					fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))
					ax1.imshow(np.ma.masked_array(image, bad_pixels), origin='lower', cmap='gray_r', interpolation='none')
					ax2.imshow(image, origin='lower', cmap='gray_r', interpolation='none')
					plt.show()
			else:
				raise ValueError('not supported filling method %s'%mode)

		return image

	def __get_rough_background(self,image):
		'''
		get the background intensity of the given image
		'''
		imgbkg = image.reshape(1,np.size(image))
		cri = np.mean(imgbkg) + np.std(imgbkg)
		imgbkg = imgbkg[imgbkg<cri]
		#print imgbkg.shape
		background,bkgfilter,std = sigma_clipping(imgbkg, sig=3, meanfunc=np.median)
		#show_histogram(background,nbins = 100)
		bkg = np.mean(background)

		return bkg,std

	def __get_surrounding(self,radius = 1):
		'''
		get the pixels which are inside the circle with center at [0,0] and radius = radius
		'''
		shifts = []

		for x in np.arange(-int(radius),int(radius)+1):
		    for y in np.arange(-int(radius),int(radius)+1):
		        if x**2+y**2 <= radius**2:
		    		shifts.append([x,y])

		return shifts

	def __cut_or_extend_imagedata(self,imagedata,xlo,ylo,xho,yho):
		'''
		cut or extend the input image data
		we can get the size of the image from the image coordinate in which the left bottom pix has (0,0)
		and the right top pix has (M-1,N-1) where (M,N) = imagedata.shape
		the image coordinates of the left bottom pix and right top pix in the output image are (xlo,ylo),(xho,yho) respectively
		(xlo,ylo),(xho,yho) are give in the input image frame
		note: the record of data on CCD is not the same as the normal array. the x axis is  along the bottom to top direction and
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

		temp[xt0:xt1,yt0:yt1]=data[xf0:xf1,yf0:yf1]
		new_data = temp
		return new_data

	def __create_new_image_from_selected_extension(self, input_image,output_image, want_hduindex, hdr_kw_del = ['XTENSION'], output_verify='ignore'):
		'''
		get the wanted extension data and header, delete specified header keywords, and write to output image
		'''
		hdu = fits.open(input_image)
		hdudata = hdu[want_hduindex].data
		hduhdr  = hdu[want_hduindex].header

		for delkw in hdr_kw_del:
			hduhdr.pop(delkw)

		newimg = fits.PrimaryHDU(data=hdudata, header=hduhdr)
		newimg.writeto(output_image, output_verify=output_verify)

#WFCAM specific
	def wfcam_prepare_image(self, extension_want=None, verbose=0):
		'''
		WFCAM reduced image has 4 science images in 4 seperate extensions. Select the one containing the target and write to individual image
		'''
		for img in self.photometry_info.keys():
			input_image = self.images[img]
			output_image =  os.path.join(self.modified_image_dir, img)
			if os.path.exists(output_image):
				continue

			self.__get_wfcam_target_image(input_image, output_image,extension_want=extension_want, verbose=verbose)

	def __get_wfcam_target_image(self, input_image, output_image, extension_want=None, verbose=0):
		'''
		the image containing the target should be closest to the target coordinate
		'''
		if extension_want is not None:
			self.__create_new_image_from_selected_extension(input_image, output_image, extension_want)
		else:
			if self.sn_ra_world_deg is None or self.sn_dec_world_deg is None:
				print "the target coordinate is not available"

			snrarad = self.sn_ra_world_deg / 180.0 * np.pi
			sndecrad = self.sn_dec_world_deg / 180.0 * np.pi

			distances = np.array([])
			hdu = fits.open(input_image)
			for i in (np.arange(4)+1):
				hdr = hdu[i].header
				data = hdu[i].data
				NX,NY = data.shape
				imagecoor = np.array([NX/2, NY/2]).reshape(1,2)
				center_coor = self.__image2world_fitshead(hdr, imagecoor) #the sky coordinate of the central pixel in the current ccd chip
				cradeg,cdecdeg = center_coor[0]
				crarad = cradeg/180.0*np.pi
				cdecrad = cdecdeg/180.0*np.pi
				raraddiff = np.abs(snrarad - crarad)
				distance = np.arccos(np.sin(sndecrad)*np.sin(cdecrad) + np.cos(sndecrad)*np.cos(cdecrad)*np.cos(raraddiff))	#compute the great-circle distance
				distances = np.append(distances, distance)
				if verbose:
					print i, distance

			extension_want = np.where(distances==np.min(distances))[0][0] + 1
			if verbose:
				print "the image containing the target for %s is in extension %s"%(input_image, extension_want)
			self.__create_new_image_from_selected_extension(input_image, output_image, extension_want)

	def wfcam_photometry_calibration_casu_zpt(self):
		'''
		For photometry on WFCAM images, the CASU data reduction pipeline gives the photometry zero point which is derived from calibration with 2MASS, details here http://www3.interscience.wiley.com/cgi-bin/fulltext?ID=122210794&PLACEBO=IE.pdf&mode=pdf

		'''

		for img in self.photometry_info.keys():
			if self.photometry_info[img]['airmass'] == 0:
				self.__renew_airmass_info(img)

			if self.photometry_info[img]['magzpt'] == 99.99:
				self.__get_magzpt_wfcam(img)

			instmag = self.photometry_info[img]['instmag']
			instmagerr = self.photometry_info[img]['instmagerr']
			magzpt = self.photometry_info[img]['magzpt']
			magzpterr = self.photometry_info[img]['magzpterr']

			exptime = self.photometry_info[img]['exptime']
			photmethod = self.photometry_method

			airmass = self.photometry_info[img]['airmass']

			mag, magerr = self.__wfcam_stdcal_casu_magzpt(instmag, instmagerr, magzpt, magzpterr, exptime, photmethod, airmass)

			self.photometry_info[img]['calmag'] = mag
			self.photometry_info[img]['calmagerr'] = magerr

	def __get_magzpt_wfcam(self, imgkey):
		'''
		Get the WFCAM photometry zeropoint from the image header which has been obtained by CASU data reduction pipeline
		'''
		input_image = os.path.join(self.modified_image_dir, imgkey)
		fitsinfo = get_fits_info(input_image, ['MAGZPT','MAGZRR'])
		magzpt = fitsinfo['MAGZPT']
		magzpterr = fitsinfo['MAGZRR']

		self.photometry_info[imgkey]['magzpt'] = magzpt
		self.photometry_info[imgkey]['magzpterr'] = magzpterr

	def __wfcam_stdcal_casu_magzpt(self, instmag, instmagerr, zpt, zpterr, exptime, photmethod, airmass):
		'''
		Apply the CASU photometry zeropoint to instrumental magnitude
		'''
		if photmethod == 'apphot':
			mag = zpt + instmag - 25 + 2.5*np.log10(exptime) - 0.05*(airmass-1)
		elif photmethod == 'psfphot':
			mag = zpt + instmag + 2.5*np.log10(exptime) - 0.05*(airmass-1)
		else:
			raise ValueError("photometry method %s not supported..."%photmethod)

		magerr = np.sqrt(instmagerr**2 + zpterr**2)

		return mag, magerr

# SMARTS-specific
	def smarts_renew_fwhm_and_psfphot_same_night(self, flt, fwhm, night=None, phottable=None,  which_dir='modified_image'):
		'''
		Renew fwhm and redo psf photometry
		'''
		if phottable is None and night is None:
			raise ValueError("provide the phottable or night")
		if phottable is None:
			phottable = self.smarts_select_photret_given_night(flt, night)
		for img in phottable['name']:
			self.photometry_info[img]['fwhm'] = fwhm
			self.psf_photometry_dophot_single(img, output_residuals=True, which_dir=which_dir)
			self.psf_photometry_dophot_single_target_ds9_pick(img, which_dir=which_dir)

	def smarts_relative_calibration_same_night(self, flt, night, offset=None, offseterr=None, image_dir = 'modified_image', refstarnum=None, refnumin=10):
		'''
		Images observed with the same exposure time within a short time period should have close zero points.
		Use the most reliable zero point derived from the best image for all images.
		'''
		phottable = self.smarts_select_photret_given_night(flt, night)

		for img in phottable['name']:
			if self.photometry_info[img]['drop']:
				continue
			if offset is not None and offseterr is not None:
				self.photometry_info[img]['relmag'] = self.photometry_info[img]['instmag'] + offset
				self.photometry_info[img]['relmagerr'] = offseterr
			else:
				self.__get_relative_mag_single_manual_match(img, flt, which_dir=image_dir, refstarnum=refstarnum, refnumin=refnumin, update_matched_inmags=1,update_rettable=1)

	def smarts_select_photret_given_night_flt(self, flt, night, updatetable=0, verbose=0):
		'''
		flt, night
		'''
		nights_flt, nights_flt_unique = self.obs_nights_stat_flt(flt, updatetable=updatetable, verbose=verbose)
		ptable_flt = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		photret_night = ptable_flt[nights_flt==night]
		return photret_night

	def obs_nights_stat_flt(self, flt, updatetable=0, verbose=1):
		'''
		List the nights of smarts observation for the current target in give band
		'''
		if updatetable or self.photometry_record_table is None:
			self.__dict2table()

		ptable_flt = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		if self.current_telescope in ['SMARTS', 'SMARTS-IR']:
			nights_flt, nights_flt_unique = self.smarts_obs_nights_stat_given_phot_table(ptable_flt)
		else:
			raise ValueError('This function only works for SMARTS telescope for now')
		if verbose:
			print "%s observed in %s band on nights listed below:"%(self.current_sn, flt)
			print nights_flt_unique

		return nights_flt, nights_flt_unique

	def smarts_obs_nights_stat_given_phot_table(self, ptable):
		'''
		smarts image name contains night info
		'''
		if self.current_telescope == 'SMARTS':
			n1 = 4
			n2 = 10
		else:
			n1 = 5
			n2 = 11
		nights = np.array([img[n1:n2] for img in ptable['realimg']])
		unights = np.unique(nights)

		return nights, unights

	def smarts_nir_flt_image_reduction_sky_subtraction(self, flt, tstart=None, tend=None, updatetable=1):
		'''
		SMARTS NIR images data reduction: remove the sky by subtracting two different images.
		The latest image to the target image will be used as the sky template

		'''
		if updatetable or self.photometry_record_table is None:
			self.__dict2table()

		flttable = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		flttable.sort('obstime')
		if tstart is not None:
			flttable = flttable[flttable['obstime']>tstart]
		if tend is not None:
			flttable = flttable[flttable['obstime']<tend]
		self.smarts_nir_multigroups_image_reduction_sky_substraction(flttable)

	def smarts_nir_multigroups_image_reduction_sky_substraction(self, ptable, exptime_diff_tolerance_digit=0, verbose=1):
		'''
		INPUTS:
			ptable:
			exptime_diff_tolerance_digit:
			verbose:
		'''
		obsnights, unights = self.smarts_obs_nights_stat_given_phot_table(ptable)
		print unights
		for night in unights:
			nightimgs = ptable[obsnights==night]
			exptimes_rounded = map(lambda x: round(x, exptime_diff_tolerance_digit), nightimgs['exptime'])
			exptimes_unique = np.unique(exptimes_rounded) #difference tolerance of 1 second by default
			if verbose:
				print night, exptimes_unique
			for exptime in exptimes_unique:
				imgs = nightimgs[exptimes_rounded==exptime]
				self.smarts_nir_single_group_image_reduction_sky_subtraction(imgs, exptime_check_criteria=0.01)

	def smarts_nir_single_group_image_reduction_sky_subtraction(self, imgs, exptime_check_criteria=0.01, verbose=1):
		'''
		Method: A number (ngroup) of images as a group, the image_subtraction happen within the group.
				In the group, each image is subtracted by the image before it or after it.

		INPUTS:
			imgs: photometry record table with a group of images (usually taken on the same night)
		'''
		ngroup = len(imgs)
		if ngroup<2:
			raise ValueError("The input group has less than 2 images...")
		for i,obs in enumerate(imgs):
			img_sci = obs['name']
			img_sci_realname = obs['realimg']
			imod = np.mod(i,ngroup)
			if imod == (ngroup-1):
				isky = i - 1  #the last image in the group
			else:
				if i == (len(imgs)-1): # the last image in the table
					isky = i -1
				else:
					if np.mod(imod,2) == 0:
						isky = i + 1      #even index, sky is the one after it
					else:
						isky = i - 1      #odd index, sky is the one before it

			img_sky = imgs[isky]['name']
			img_sky_realname = imgs[isky]['realimg']

			exptime_sci = imgs[i]['exptime']
			exptime_sky = imgs[isky]['exptime']
			if (exptime_sci-exptime_sky)/exptime_sci>exptime_check_criteria:
				raise ValueError('science frame and sky frame do not have the same exposure time')

			if verbose:
				print "%s - %s"%(img_sci, img_sky)
				print "%s - %s"%(img_sci_realname, img_sky_realname)

			self.__subtract_images(img_sci, img_sky)

	def smarts_get_target_xys_flt(self, flt, image_dir='modified_image', updatetable=0):
		if updatetable or self.photometry_record_table is None:
			self.__dict2table()
		flttable = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		for imgkey in flttable['name']:
			if self.photometry_info[imgkey]['drop'] != 0 or self.photometry_info[imgkey]['x'] != 0.0:
				continue
			drop, x, y = self.smarts_get_target_xy_single(imgkey, image_dir=image_dir)
			if drop:
				self.photometry_info[imgkey]['drop'] = 41
			else:
				self.photometry_info[imgkey]['x'] = x
				self.photometry_info[imgkey]['y'] = y

	def smarts_get_target_xy_single(self, imgkey, image_dir = 'modified_image', select_from_sources=True, app=5, fwhmpsf=5, skywidth=10, update_photometry_record=0, verbose=0):
		'''
		First pick (x,y) through on DS9, then run single-point aperture photometry to determine whether the source is detected
		'''
		imgkey_s = imgkey.split('.')[0]
		input_image = self.__get_internal_image(imgkey, which_dir=image_dir)

		spphot_options = self.apphot_iraf_options.copy()
		spphot_options['app'] = app
		spphot_options['fwhmpsf'] = fwhmpsf
		spphot_options['skywidth'] = skywidth

		starfile = os.path.join(self.stars_dir, imgkey_s+self.starfile_suffix)
		if os.path.exists(starfile) and select_from_sources:
			starxys = np.loadtxt(starfile)
			regionfile = os.path.join(self.stars_dir, imgkey_s+"_temp.reg")
			circle_radius = self.criteria_point_match_general
			create_ds9_region_file(starxys, regionfile, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', circle_radius = circle_radius, color='red', width=2, load_text=False)
			xy = self.__get_xy_on_image_from_ds9(input_image=input_image, regionfile=regionfile, newds9=True)
			self.__delete_file_if_exist(regionfile)
			
			if xy is None:
				drop = 1
				xnew = 0.0
				ynew = 0.0
			else:
				xy_target = [xy[0],xy[1]]
				match_criteria = self.criteria_point_match_general
				yesfound, xy_match,index = self.__find_corresponding(starxys[:,[0,1]], xy_target, match_criteria)

				if yesfound and len(index)==1:
					print "the corresponding object #%s found at (%s)"%(index[0], xy_match)
					matchfound = starxys[index[0]]
					drop = 0
					xnew = matchfound[0]
					ynew = matchfound[1]
				else:
					if (not yesfound):
						print "No object found within criteria..."
						print "You can change input parameter match_criteria or self.criteria_point_match_general to give a new search distance"
					else:
						print 'more than one objects fall into the criteria region of reference star. This object is droped!'
						print "You can change input parameter match_criteria or self.criteria_point_match_general to give a new search distance"
					drop = 1
					xnew = 0.0
					ynew = 0.0
		else:
			xy = self.__get_xy_on_image_from_ds9(input_image=input_image, regionfile=None, newds9=True)
			if xy is None:
				drop = 1
				xnew = 0.0
				ynew = 0.0
			else:
				spphot = self.aperture_photometry_apphot_iraf_single_point(input_image, xy, options=spphot_options)
				xnew = spphot[0,0]
				ynew = spphot[0,1]
				mag_err = spphot[0,-1]
				if verbose:
					print "(%s, %s) --> (%s, %s)"%(xy[0], xy[1], xnew, ynew)
				if mag_err > 0.33:
					if verbose:
						print "%s mag_err > 0.33, will be dropped"%imgkey
					drop = 1
				else:
					drop = 0
					if update_photometry_record:
						self.photometry_info[imgkey]['x'] = xnew
						self.photometry_info[imgkey]['y'] = ynew

		return drop, xnew, ynew

	def smarts_ccd_remove_photometry_in_bad_region(self,xlow=105, ylow=30):
		'''
		SMARTS 1.3m ANDICAM CCD has some region which are not used
		'''
		photmethod = self.photometry_method
		for img in self.images.keys():
			if self.photometry_info[img]['drop'] == 0:
				self.select_photometry_measurement_given_region(img, xl=xlow, yl=ylow, photometry_method=photmethod)

	def smarts_get_imgs_fwhm_bkg_flt(self, flt, renew=0, verbose=0):
		'''
		find stars for measuring fwhm and bkg by shift the sources on reference image
		'''
		refmag_table = self.__find_corresponding_refmags_for_given_imgkey(flt=flt)
		if self.photometry_record_table is None:
			self.__dict2table()
		flttable = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		for image_key in flttable['name']:
			if self.photometry_info[image_key]['drop'] != 0:
				continue
			if self.photometry_info[image_key]['fwhm'] != 0 and (not renew):
				continue
			print image_key
			self.smarts_get_img_fwhm_bkg_single(image_key, refmag_table, update_phot_record=1, display_fit=0, verbose=verbose)

	def smarts_get_img_fwhm_bkg_single(self, image_key, refmag_table, update_phot_record=1, display_fit=1, verbose=0):
		idxys = self.__get_idxys_corresponding_to_reference_image(image_key,refmag_table)
		xys = idxys[:,[1,2]]
		fwhm, bkg = self.get_fwhm_bkg_2dgaussian_photutils_with_input_star_xys(image_key, xys, fitwidthx=20, fitheighty=20, which_dir='modified_image', update_phot_record=update_phot_record, verbose=verbose, display_fit=display_fit)
		print "fwhm:%s bkg:%s"%(fwhm, bkg)



#DEMONEXT-specific
	def demonext_coadd_images(self, input_which_dir='raw_image', output_dir=None, updatetable=1):
		'''
		try to co-add multiple images taken on the same night to get high SNR of the combined image
		'''
		if updatetable:
			self.__dict2table()
		dates_all  = [img[1:9]  for img in self.photometry_record_table['realimg']]
		dates_unique = np.unique(dates_all)

		flts_unique = np.unique(self.photometry_record_table['flt'])

		for date in dates_unique:
			for flt in flts_unique:
				self.__demonext_coadd_image_night_flt(date, flt, input_which_dir=input_which_dir, output_dir=output_dir)

	def __demonext_coadd_image_night_flt(self, date, flt, input_which_dir='raw_image', output_dir=None, clobber=False):


		if input_which_dir == 'raw_image':
			imgdir = self.raw_image_dir
		elif input_which_dir == 'modified_image':
			imgdir = self.modified_image_dir
		else:
			raise ValueError("the input image dir not recognised")

		input_imgkeylist = [img for img in self.images.keys() if date in self.photometry_info[img]['realimg'] and self.photometry_info[img]['flt']==flt]
		print input_imgkeylist
		refimgkey = input_imgkeylist[0]


                combined_image_file = '%s_%s.fits'%(flt, date)

		if output_dir is None:
			output_dir = os.path.join(self.repository_dir_current_telescope, "coadd_temp")
			print "the co-added images will be stored here: %s"%output_dir
			if os.path.exists(output_dir):
				os.remove(output_dir)

		if not os.path.exists(output_dir):
			os.mkdir(output_dir)

                combined_image = os.path.join(output_dir, combined_image_file)

                if os.path.exists(combined_image) and clobber:
                        os.remove(combined_image)


		self.align_and_stack_multiple_images(input_imgkeylist, refimgkey, combined_image, input_which_dir=input_which_dir, renew_matched_sourcelist = False, renew_resampled_img=False)

#Iowa telescope specific
	def Iowa_prepare_template_image_wcs(self):
		'''
		The Iowa images have the non-standard header issue which will cause the problems when processing by the astropy.wcs function
		To get rid of the problem, only nessary keywords are kept when preparing Iowa template images
		'''
		for flt in self.templates.keys():
			tplimg = self.templates[flt]
			tempimg = self.Iowa_prepare_template_image_wcs_single(tplimg)

	def Iowa_prepare_template_image_wcs_single(self, tplimg, outimg=None):
		'''
		solve the non-standard issue of Iowa image
		'''
		inimg = os.path.join(self.raw_image_dir, tplimg)
		if outimg is None:
			outimg = os.path.join(self.template_dir, 'cal_%s'%tplimg)
		simplify_and_standardize_fits_header(inimg, outimg, out_header_only=False, overwrite=True)
		return outimg

#LCOGT telescope 
	def lcogt_get_images_propid(self, update_phot_record_table=1):
		'''
		register image PROPID for the image table sorted by obstime 
		'''
		if update_phot_record_table:
			self.__dict2table()
		
		internal_proposal_IDs = ['CLN2014B-004','CON2014B-004','CLN2015A-002','NAOC2015A-001','CON2015A-004','ARI2015A-004','NAOC2016A-001','ARI2016A-002','CLN2016A-005','CON2016A-008','ARI2017AB-001','CLN2017AB-006','NAOC2017AB-001', 'NAOC2018A-001','ARI2018A-003', 'DDT2018A-002', 'TAU2018A-004', 'CON2018B-001', 'CLN2018B-006', 'NAOC2018B-001', 'NOAO2018B-013', 'NOAO2018B-020','CLN2018B-001', 'TAU2018B-008','NOAO2019A-019', 'NOAO2019A-023','NAOC2019A-002', 'SUPA2019A-002','TAU2019A-005', 'CON2019A-003','NAOC2019B-003','NOAO2019B-014','TAU2019B-001','NAOC2020A-001','TAU2020A-006','TAU2020A-013','NAOC2020B-002']
		
		propid_list = []
		internal_list = []
		for img in self.photometry_record_table['name']:
			image = self.images[img]
			hdu = fits.open(image)
			header = hdu[0].header
			propid = header['PROPID']
			propid_list.append(propid)
			if propid not in internal_proposal_IDs:
				internal_list.append('no')
			else:
				internal_list.append('yes')
		propinfotable = Table([self.photometry_record_table['name'].data, self.photometry_record_table['realimg'].data, propid_list, internal_list], names=['name','realimg','propid','internal'])
		return propinfotable

	def lcogt_seperate_lcs(self, flts='all', versionflag=''):
		'''
		seperate LC data from internal proposals and other proposals
		'''
		
		self.__find_template_imagekey(updatetable=1)
		if flts == 'all':
			flts_save = self.templates.keys()
		elif isinstance(flts,str):
			flts_save = [flts]
		elif isinstance(flts,list):
			flts_save = flts
		else:
			raise IOError('Invalid input for flts...')

		proptable = self.lcogt_get_images_propid(update_phot_record_table=1) 
		self.__dict2table()
		
		internaltable = self.photometry_record_table[proptable['internal']=='yes']
		pubtable = self.photometry_record_table[proptable['internal']=='no']

		sn_name = self.current_sn
		telescope_code =self.base.telescopes_info[self.current_telescope]['code']
		if self.photometry_method == 'apphot':
			methodcode = '-AP'
		elif self.photometry_method == 'psfphot':
			methodcode = '-PSF'
		else:
			raise IOError("Invalid input for photometry_method...")


		if self.base.save_data_local_dir != '' and os.path.exists(self.base.save_data_local_dir):
			lc_dir_archive = os.path.join(self.base.save_data_local_dir,sn_name)
			if not os.path.exists(lc_dir_archive):
				os.mkdir(lc_dir_archive)
			savetosniper = True
		else:
			savetosniper = False


		for flt in flts_save:
			tpl_table = self.__select_rows_from_table(internaltable,'flt',flt)
			if len(tpl_table)==0:
				continue
			tpl_nondrop_table = self.__select_rows_from_table(tpl_table,'drop',0)
			if len(tpl_nondrop_table) == 0:
				continue
			lc_flt  = Table()
			lc_flt.add_columns([tpl_nondrop_table['obstime'],tpl_nondrop_table['calmag'],tpl_nondrop_table['calmagerr']])
			lc_flt.sort('obstime')
			if versionflag != '':
				versionflagnew = '-'+versionflag+'-internal'
			else:
				versionflagnew = '-internal'
			if flt in ['gp', 'rp', 'ip', 'zp']:
				flt = flt[0]
			lc_filename = sn_name+'_'+telescope_code+'-'+flt+methodcode+versionflagnew+'.txt'
			outfile = os.path.join(self.result_dir,lc_filename)
			np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")
			if savetosniper:
				outfile = os.path.join(lc_dir_archive, lc_filename)
				np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")
			'''
			tpl_table = self.__select_rows_from_table(pubtable,'flt',flt)
			tpl_nondrop_table = self.__select_rows_from_table(tpl_table,'drop',0)
			lc_flt  = Table()
			lc_flt.add_columns([tpl_nondrop_table['obstime'],tpl_nondrop_table['calmag'],tpl_nondrop_table['calmagerr']])
			lc_flt.sort('obstime')
			if versionflag != '':
				versionflag = '-'+versionflag
			else:
				versionflag = '-pub'
			if flt in ['gp', 'rp', 'ip', 'zp']:
				flt = flt[0]
			lc_filename = sn_name+'_'+telescope_code+'-'+flt+methodcode+versionflag+'.txt'
			outfile = os.path.join(self.result_dir,lc_filename)
			np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")
			if savetosniper:
				outfile = os.path.join(lc_dir_archive, lc_filename)
				np.savetxt(outfile,lc_flt,fmt="%12.4f %8.3f %8.3f")

			'''

# handy scripts
	def one_by_one_check_and_photometry(self, flt, drop_signal=0, psfphot_threshold_default=50):
		'''
		one by one checking the image in filter 'flt' with given drop_signal and decide whether the photometry needed
		'''
		self.__dict2table()
		fltdata = self.photometry_record_table[self.photometry_record_table['flt']==flt]
		for img in fltdata['name']:
			if self.photometry_info[img]['drop'] == drop_signal:
				self.display_image_with_ds9(img, which_dir='raw_image')

				drop_default = 0
				drop = raw_input("Enter the drop signal for image %s (default: %s):" %(img,drop_default)) or drop_default
				if drop:
					self.photometry_info[img]['drop'] = int(drop)
				else:
        				bkg = float(raw_input("bkg:"))
        				fwhm = float(raw_input("fwhm:"))
        				self.photometry_info[img]['bkg'] = bkg
        				self.photometry_info[img]['fwhm'] = fwhm
					psfphot_threshold = raw_input('the threshold for peak for sources in psf photometry [default: %s]'%psfphot_threshold_default) or psfphot_threshold_default
					self.psfphot_thresmin = float(psfphot_threshold)
        				self.psf_photometry_dophot_single(img, output_residuals=1, which_dir='raw_image')
        				self.psf_photometry_dophot_single_target_ds9_pick(img)

	def get_images_dict(self,data_dir):
		'''
		return the images info dict containing the relative image names as keys and absolute names as values
		'''
		images_dict = OrderedDict()
		filelist = os.listdir(data_dir)
		for img in filelist:
			images_dict[img] = os.path.join(data_dir,img)
		return images_dict


	def perform_whole_image_subtraction_process(self, tpldict=None, cutimage=True, edge_width=100, fwhm_maxtime=15, nsx=10, nsy=10):
		'''
		remove CR, get fwhm and bkg, source detection, prepare tpl images, hotpants image subtraction
		'''
		if cutimage:
			self.image_cutout(edge_width=edge_width)
			self.remove_cosmic_ray(which_dir='md')
		else:
			self.remove_cosmic_ray(which_dir='raw')
		
		self.fill_phot_info_specific_col('bitpix', fill_value=-32)
		self.get_fwhm_bkg(which_dir='md', time_out_max=fwhm_maxtime)
		self.source_detection(which_dir='md', redo_daofind=1)
		self.save_results(overwrite_or_create_initversion=0, create_newversion=0)
		
		if tpldict is not None:
			for flt,img in tpldict.items():
				self.assign_template_image(flt, img)
		self.__find_template_imagekey()
		
		for flt,tplimg in self.templates.items():
			print flt, tplimg
			self.prepare_template_image(tplflt=flt, tplimgkey=tplimg, which_dir='md')
			self.photometry_info[tplimg]['drop'] = 100

		self.image_subtraction_all(which_dir='md', renew_sub=0, xstamp=nsx, ystamp=nsy)

