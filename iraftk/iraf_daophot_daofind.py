from pyraf import iraf
from iraf import daophot
from iraf import daofind

def daofind_iraf(input_image, output='default', starmap='',skymap='', fwhmpsf=2.5, emission=True, sigma=0, datamin='INDEF', datamax='INDEF', readnoise=0, epadu=1.0, threshold=4, nsigma=1.5, sharplo=0.2,sharphi=1.0, roundlo=-1.0, roundhi=1.0,  mkdetections=False):
	'''
	DAOFIND approximate the stellar point spread function with an elliptical Gaussian function, whose sigma along th semi-major axis is 0.42466*datapars.fwhmpsf/datapars.scale pixels, semi-minor to semi-major axis ratio is ratio, and major axis position angle is theta. Using this model, a convolution kernel, truncated at nsigma sigma, and normalized so as to sum to zero(?? from iraf/daofind help page), is contructed. Then the density enhancement image starmap is computed by convolving the input image with gaussian kernel. After image convolution, DAOFIND steps through starmap searching for density enhancements greater than findpars.threshold*datapars.sigma, and brighter than all other density enhancement within a semi-major axis of 0.42466*findpars.nsigma*datapars.fwhmpsf.

	INPUTS:
		input_image:  input image(s)
		output: output coordinate file(s), default: image.coo.?
		starmap: output density enhencement image
		skymap: output sky image
		fwhmpsf: FWHM of the PSF in scale unit
		emission: features are positive
		sigma: standard of deviation of background in counts
		datamin: minimum good value
		datamax: maximum good value
		readnoise: CCD readnoise in electrons
		epadu: gain in electrons per count
		threshold: threshold in sigma for feature detection
		nsigma: width of convolution in sigma
		mkdetections: mark detection on image display

	OUTPUTS:
		the output coordinate file contains:
		xcenter ycenter mag sharpness sround ground id
		where mag = -2.5*log10(peak density /detection threshold)
	'''
	daofind.datapars.unlearn()
	daofind.datapars.fwhmpsf = fwhmpsf # fwhm of the psf in scale units
	daofind.datapars.emission = emission #positive feature
	daofind.datapars.sigma = sigma  #standard deviation of backgroud in counts
	daofind.datapars.datamin = datamin
	daofind.datapars.datamax = datamax
	daofind.datapars.readnoise = readnoise #ccd readout noise in electrons
	daofind.datapars.epadu  = epadu #gain in electrons per count

	daofind.findpars.unlearn()
	daofind.findpars.threshold = threshold #threshold in sigma for feature detection
	daofind.findpars.nsigma = nsigma    #width of convolution kernel in sigma
	#daofind.findpars.ratio = 1.0
	#daofind.findpars.theta = 0.0
	daofind.findpars.sharplo = sharplo #lower bound on sharpness for feature detection
	daofind.findpars.sharphi = sharphi #upper bound on sharpness for feature detection
	daofind.findpars.roundlo = roundlo #lower bound on roundness for feature detection
	daofind.findpars.roundhi = roundhi #upper bound on roundness for feature detection
	#daofind.findpars.mkdetections = mkdetections #mark detected stars on the display


	daofind.unlearn()
	daofind.verify = 'no'
	daofind.update = 'no'
	daofind.verbose = 'no'
	daofind(image=input_image, output = output, starmap=starmap, skymap=skymap)




if __name__ == "__main__":

	from get_image_bkg_stddev_rough import estimate_bkg_stddev

	import optparse
	parser = optparse.OptionParser()


	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_output = 'default'
	parser.add_option('-o', '--output', dest='output', type='string', default=def_output, help='output file; default: %s'%def_output)

	def_starmap = ''
	parser.add_option('--starmap', dest='starmap', type='string', default=def_starmap, help='the name of the image prefix and/or directory where the density enhancement image will be stored: %s'%def_starmap)

	def_skymap = ''
	parser.add_option('--skymap', dest='skymap', type='string', default=def_skymap, help='the name of the image prefix and/or directory where the mean density image will be stored: %s'%def_skymap)

	def_fwhmpsf = 5
	parser.add_option('--fwhmpsf', dest = 'fwhmpsf', type = float, default= def_fwhmpsf, help='FWHM of the PSF in scale unit; default: %s'%def_fwhmpsf )

	def_sigma = None
	parser.add_option('--sigma', dest = 'sigma', type = float , default= def_sigma, help='standard of deviation of background in counts; default: %s'%def_sigma )

	def_datamin = None
	parser.add_option('--datamin', dest = 'datamin', type = float, default= def_datamin, help='minimum good value; default: %s'%def_datamin )

	def_datamax = None
	parser.add_option('--datamax', dest = 'datamax', type = float, default= def_datamax, help='maximum good value; default: %s'%def_datamax )


	def_readnoise = 5
	parser.add_option('--readnoise', dest = 'readnoise', type = float, default= def_readnoise, help='CCD readnoise in electrons; default: %s'%def_readnoise )

	def_epadu = 1
	parser.add_option('--epadu', dest = 'epadu', type = float, default= def_epadu, help='gain in electrons per count; default: %s'%def_epadu )

	def_threshold = 3
	parser.add_option('--threshold', dest = 'threshold', type = float, default= def_threshold, help='threshold in sigma for feature detection; default: %s'%def_threshold )

	def_sharphi = 1
	parser.add_option('--sharphi', dest = 'sharphi', type = float, default= def_sharphi, help='image sharpness; default: %s'%def_sharphi )

	def_sharplo = 0.2
	parser.add_option('--sharplo', dest = 'sharplo', type = float, default= def_sharplo, help='image sharpness; default: %s'%def_sharplo )

	def_deleteINDEF = False
        parser.add_option('-d','--deleteINDEF', dest="deleteINDEF",action="store_true", default=def_deleteINDEF, help="whether delete lines containing INDEF[%s]"%def_deleteINDEF)

	options, remainder = parser.parse_args()

	input_image = options.input_image
	output = options.output
	starmap = options.starmap
	skymap  = options.skymap
	fwhmpsf= options.fwhmpsf
	emission=True
	sigma= options.sigma
	datamin= options.datamin
	datamax= options.datamax
	readnoise= options.readnoise
	epadu= options.epadu
	threshold= options.threshold
	sharplo = options.sharplo
	sharphi  = options.sharphi
	deleteINDEF = options.deleteINDEF

	nsigma=1.5

	if sigma is None:
		stddev1, stddev2 = estimate_bkg_stddev(input_image, epadu=epadu)
		sigma = stddev2

	if datamin is None:
		datamin = 'INDEF'

	if datamax is None:
		datamax = 'INDEF'

	daofind_iraf(input_image, output=output, starmap=starmap, skymap=skymap, fwhmpsf=fwhmpsf, emission=emission, sigma=sigma, datamin=datamin, datamax=datamax, readnoise=readnoise, epadu=epadu, threshold=threshold, sharplo=sharplo, sharphi=sharphi,  nsigma=nsigma)
