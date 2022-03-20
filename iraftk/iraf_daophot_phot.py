#! /usr/bin/python

#aperture photometry with iraf/apphot
import os
from pyraf import iraf
from iraf import daophot
from iraf import phot
from iraf_utils import iraf_get_default_filename

def daophot_phot_iraf(image, centroid_algorithm='centroid', coo_list='default', output='default', fwhmpsf=3, sigma=10, app=8, skyin=15, wsky=15, sky_sigma_down=3, sky_sigma_up=3, readnoise=5, epadu=1, datamin = 'INDEF', datamax = 'INDEF', emission='yes', zeropt=25):
	'''
	INPUTS:
	image:
	centroid_algorithm: 'none', 'centroid', 'gauss' or 'ofilter'
	coo_list: the input coordinate list(s); default: image.coo.?
	output: output photometry file(s); default: image.mag.?
	fwhmpsf: FWHM of the PSF in scale unit
	sigma: standard deiation of background in counts
	app:  aperture radius(radii) in scale unit
	skyin: inner sky annulus radius
	wsky: the width of sky annulus
	sky_sigma_down: lower k-sigma rejection limit in sky sigma
	sky_sigma_up: upper k-sigma rejection limit in sky sigma
	readnoise: CCD readout noise in electrons
	epadu: gain in electron per count
	datamin: minimum good data value
	datamax: maximum good data value

	OUTPUTS:
	The full photometry will be saved in aperture photometry file output, and for each star measured the following record is written :

	image  xinit  yinit  id  coords  lid
	xcenter  ycenter  xshift  yshift  xerr  yerr  cier cerror
	msky  stdev  sskew  nsky  nsrej  sier  serror
	itime  xairmass  ifilter  otime
	rapert  sum  area  mag  merr  pier  perr

	Image and coords are the name of the image and coordinate file respectively.
	Id and lid are the sequence numbers of stars in the output and coordinate files respectively.
	Sier and serror are the sky fitting error code and accompanying error message respectively.
	Msky, stdev and sskew are the best estimate of the sky value (per pixel), standard deviation and skew respectively.
	Nsky and nsrej are the number of sky pixels and the number of sky pixels rejected respectively.
	Itime is the exposure time, xairmass is self-evident, ifilter is an id string identifying the filter used in the observations, and otime is a string containing the time of the observation in whatever units the user has set up.
	Rapert, sum, area, and flux are the radius of the aperture in scale units, the total number of counts including sky in the aperture, the area of the aperture in square pixels, and the total number of counts excluding sky in the aperture.
	Mag and merr are the magnitude and error in the magnitude in the aperture (see below).

	        flux = sum - area * msky
       		mag = zmag - 2.5 * log10 (flux) + 2.5 * log10 (itime)
        	merr = 1.0857 * err / flux
         	err = sqrt (flux / epadu + area * stdev**2 + area**2 * stdev**2 / nsky)
	'''
	phot.datapars.unlearn()
	phot.datapars.scale = 1.0
	phot.datapars.fwhmpsf = fwhmpsf
	phot.datapars.emission = emission
	phot.datapars.sigma = sigma
	if datamin is None:
		datamin = 'INDEF'
	if datamax is None:
		datamax = 'INDEF'
	phot.datapars.datamin = datamin
	phot.datapars.datamax = datamax
	phot.datapars.noise = 'poisson'
	phot.datapars.ccdread = ''
	phot.datapars.readnoise = readnoise
	phot.datapars.itime = 1.0
	phot.datapars.epadu = epadu
	phot.datapars.xairmass = 'INDEF'
	phot.datapars.ifilter = 'INDEF'
	phot.datapars.otime = 'INDEF'

	
	phot.centerpars.unlearn()# iraf.digiphot.daophot.centerpars :
	phot.centerpars.calgorithm = centroid_algorithm
	phot.centerpars.cbox = 10.0
	phot.centerpars.cthreshold = 0.0
	phot.centerpars.minsnratio = 1.0
	phot.centerpars.cmaxiter = 10.0
	phot.centerpars.maxshift = fwhmpsf/2.0
	phot.centerpars.clean = False
	phot.centerpars.rclean = 1.0
	phot.centerpars.rclip = 2.0
	phot.centerpars.kclean = 3.0
	phot.centerpars.mkcenter = False

	
	phot.fitskypars.unlearn()#iraf.digiphot.daophot.fitskypars :
	phot.fitskypars.salgorithm = 'median'
	phot.fitskypars.annulus = skyin
	phot.fitskypars.dannulus = wsky
	phot.fitskypars.skyvalue = 0.0
	phot.fitskypars.smaxiter = 10.0
	phot.fitskypars.sloclip = 0.0
	phot.fitskypars.shiclip = 0.0
	phot.fitskypars.snreject = 50.0
	phot.fitskypars.sloreject = sky_sigma_down
	phot.fitskypars.shireject = sky_sigma_up
	phot.fitskypars.khist = 3.0
	phot.fitskypars.binsize = 0.1
	phot.fitskypars.smooth = False
	phot.fitskypars.rgrow = 0.0
	phot.fitskypars.mksky = False

	
	phot.photpars.unlearn() #iraf.digiphot.daophot.photpars :
	phot.photpars.weighting = 'constant'
	phot.photpars.apertures = app
	phot.photpars.zmag = zeropt
	phot.photpars.mkapert = False

	photparams = {'radplot':False}

	if os.path.exists(output):
	    os.remove(output)

	phot.unlearn()
	phot(image=image, skyfile='', output=output, coords=coo_list, verify='no',interactive='no',verbose=True, Stdout=1, **photparams)
	#'IMAGE', 'XINIT', 'YINIT', 'ID', 'COORDS', 'LID', 'XCENTER','YCENTER', 'XSHIFT', 'YSHIFT', 'XERR', 'YERR', 'CIER', 'CERROR',
	#'MSKY', 'STDEV', 'SSKEW', 'NSKY', 'NSREJ', 'SIER', 'SERROR','ITIME', 'XAIRMASS', 'IFILTER', 'OTIME', 'RAPERT', 'SUM', 'AREA',
	#'FLUX', 'MAG', 'MERR', 'PIER', 'PERROR'




if __name__ == "__main__":
	from astropy.io import fits,ascii
	import optparse
	parser = optparse.OptionParser()

	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_x = None
	parser.add_option('-x', dest = 'x', type = float, default= def_x, help='x position of target; default: %s'%def_x )

	def_y = None
	parser.add_option('-y', dest = 'y', type = float, default= def_x, help='y position of target; default: %s'%def_y )

	def_coor = 'default'
	parser.add_option('--coor', dest='coor', type='string', default=def_coor, help='input coordinate list(s); default: %s.coo.?'%def_inputimage)

	def_daophotoutput = 'default'
	parser.add_option('--daophot_output', dest='daophot_output', type='string', default=def_daophotoutput, help='iraf/daophot internal output file; default: %s.mag.?'%def_inputimage)

	def_centroid = 'centroid'
	parser.add_option('--centroid', dest='centroid', type='string', default=def_centroid, help='centroid algorithm, valid inputs are centroid, none, gauss and ofilter; default: %s'%def_centroid)

	def_app = 8
	parser.add_option('-a','--aperture_radius', dest='app', type=float, default=def_app, help='aperture radius(radii) in scale unit; default: %s'%def_app)


	def_skyin = 15
	parser.add_option('--skyin', dest='skyin', type = float, default=def_skyin, help='inner sky annulus radius; default: %s'%def_skyin)


	def_wsky = 15
	parser.add_option('--sky_width', dest='sky_width', type = float, default=def_wsky, help='the width of sky annulus; default: %s'%def_wsky)


	def_fwhmpsf = 5
	parser.add_option('--fwhmpsf', dest = 'fwhmpsf', type = float, default= def_fwhmpsf, help='FWHM of the PSF in scale unit; default: %s'%def_fwhmpsf )


	def_sigma = 10
	parser.add_option('--sigma', dest = 'sigma', type = float , default= def_sigma, help='standard of deviation of background in counts; default: %s'%def_sigma )

	def_datamin = 'INDEF'
	parser.add_option('--datamin', dest = 'datamin', type = "string", default= def_datamin, help='minimum good value; default: %s'%def_datamin )

	def_datamax = 'INDEF'
	parser.add_option('--datamax', dest = 'datamax', type = "string", default= def_datamax, help='maximum good value; default: %s'%def_datamax )


	def_readnoise = 5
	parser.add_option('--readnoise', dest = 'readnoise', type = float, default= def_readnoise, help='CCD readnoise in electrons; default: %s'%def_readnoise )

	def_epadu = 1
	parser.add_option('--epadu', dest = 'epadu', type = float, default= def_epadu, help='gain in electrons per count; default: %s'%def_epadu )

	def_okeys = 'XCENTER,YCENTER,MSKY,STDEV,NSKY,AREA,SUM,MAG,MERR'
	parser.add_option('--colnames', dest='record_names', type='string', default=def_okeys, help='the names of the record for which to save to photometry table; default:%s'%def_okeys)

	def_output = ''
	parser.add_option('-o','--output', dest='output', type='string', default=def_output, help='the filename of the output file where the selected record of the photometry will be saved; default: %s'%def_output)

	options, remainder = parser.parse_args()

	input_image = options.input_image
	centroid_algorithm = options.centroid
	coo_list = options.coor

	x = options.x
	y = options.y

	output = options.daophot_output
	fwhmpsf= options.fwhmpsf
	sigma= options.sigma

	app = options.app
	skyin = options.skyin
	wsky = options.sky_width

	outkeys = options.record_names
	outkeys = outkeys.split(',')
	if len(outkeys)==0:
		outkeys = None
	outfile = options.output

	datamin= options.datamin
	if datamin != 'INDEF':
		datamin = float(datamin)

	datamax= options.datamax
	if datamax != 'INDEF':
		datamax = float(datamax)

	readnoise= options.readnoise
	epadu= options.epadu

	emission = 'yes'
	sky_sigma_down = 3
	sky_sigma_up = 3
	zeropt= 25

	if not os.path.exists(input_image):
		raise ValueError("input image %s doesn't exist"%input_image)

	if x is not None and y is not None:
		coo_list = 'temp.coo'
		fid = open(coo_list,'w')
		fid.write(str(x)+' ')
		fid.write(str(y))
		fid.close()
	elif coo_list == 'default':
		print "the input coordinate list: %s.coo.?"%input_image
	else:
		print "the input coordinate list: %s"%coo_list


	if output == 'default':
		output = iraf_get_default_filename(input_image, 'mag')
	daophot_phot_iraf(input_image, centroid_algorithm=centroid_algorithm, coo_list=coo_list, output=output, fwhmpsf=fwhmpsf, sigma=sigma, app=app, skyin=skyin, wsky=wsky, sky_sigma_down=sky_sigma_down, sky_sigma_up=sky_sigma_up, readnoise=readnoise, epadu=epadu, datamin =datamin, datamax = datamax)

	photret = ascii.read(output, format='daophot')
	#get the photometric data from the daophot output
	if outkeys is not None:
		selectedphot = photret[outkeys]
	else:
		selectedphot = photret

	for colname in selectedphot.keys():
		selectedphot[colname]=selectedphot[colname].filled(fill_value=99.99)
	print selectedphot

	if outfile != '':
		selectedphot.write(outfile, format='ascii.commented_header')
