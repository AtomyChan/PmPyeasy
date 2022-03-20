from pyraf import iraf
from iraf import imcoords
from iraf import ccmap
import os

def imcoords_ccmap_iraf(xyRADec, wcsfile, images, fitinfofile=None, xcolumn=1, ycolumn=2, lngcolumn=3,latcolumn=4, lngunits='degrees', latunits='degrees',insystem='icrs', update='yes', verbose='no', interactive='no'):
	'''
	compute plate solutions using matched pixel and celestial coordinate lists
	See: https://iraf.net/irafhelp.php?val=ccmap&help=Help+Page

	INPUTS:
		xyRADec: The input text files containing the pixel and celestial coordinates of points in the input images
		wcsfile: The text database file where the computed plate solutions are stored.
		images: The images associated with the input coordinate files.
		fitinfofile: Optional output files containing a summary of the results including a description of the plate geometry and a listing of the input coordinates, the fitted coordinates, and the fit residuals.
		xcolumn, ycolumn, lngcolumn, latcolumn: The input coordinate file columns containing the x, y, ra / longitude and dec / latitude values.
		lngunits, latunits: The units of the input ra / longitude and dec / latitude coordinates. The options are "hours", "degrees", and "radians" for ra / longitude, and "degrees" and "radians" for dec / latitude.
		insystem: The input celestial coordinate system. The systems of most interest to users are "icrs", "j2000", and "b1950" which stand for the ICRS J2000.0, FK5 J2000.0 and FK4 B1950.0 celestial coordinate systems respectively.
		update: Update the world coordinate system in the input image headers ? 'yes' or 'no'
		verbose: print message about progress of the task?
		interactive: fit the transformation interactively? 
	'''

	ccmap.unlearn()
	ccmap.images = images
	if fitinfofile is not None:
		ccmap.results = fitinfofile
	ccmap.xcolumn = xcolumn
	ccmap.ycolumn = ycolumn
	ccmap.lngcolumn = lngcolumn
	ccmap.latcolumn = latcolumn
	ccmap.lngunits = lngunits
	ccmap.latunits = latunits
	ccmap.insystem = insystem
	ccmap.update = update
	ccmap.verbose = verbose
	ccmap.interactive = interactive

	ccmap(input=xyRADec, database=wcsfile)	


if __name__ == "__main__":
	
	import optparse
	parser = optparse.OptionParser()

	def_xyradecfile = ''
	parser.add_option('-i','--xyradec', dest = 'xyradecfile', type= 'string', default = def_xyradecfile, help='input file containing the pixel and celestial coordinates of points in the input image')

	def_wcsfile = ''
	parser.add_option('-o','--wcsfile', dest = 'wcsfile', type= 'string', default = def_wcsfile, help='output file containing the plate solution')

	def_inputimage = ''
	parser.add_option('--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	
	options, remainder = parser.parse_args()
	xyradecfile = options.xyradecfile
	wcsfile = options.wcsfile
	input_image = options.input_image

	imcoords_ccmap_iraf(xyradecfile, wcsfile, input_image, fitinfofile=None, xcolumn=1, ycolumn=2, lngcolumn=3,latcolumn=4, lngunits='degrees', latunits='degrees',insystem='icrs', update='yes', verbose='no', interactive='no')
