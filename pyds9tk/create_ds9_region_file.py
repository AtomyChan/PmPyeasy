#! /usr/bin/python

import os
from common import radec_format_transformation

def create_ds9_region_file_ellipse(input_sources, region_filename, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, a_col=2, b_col=3, lengthscale=6, theta_col=4, color='green', width=1):
	'''
	create ds9 region file 
		
	INPUTS:
		input_source: ndarray
		x_col:
		y_col:
		a_col: major-axis length
		b_col: minor-axis length
		theta_col: angle 	
	'''
	if os.path.exists(region_filename):
		if clobber:
			os.remove(region_filename)
		else:
			print "region file %s already exist; quit now under non clobber mode..."%region_filename
			return
		
	fid = open(region_filename,'w')
	fid.write("global color=%s width=%s\n"%(color, width))

        fid.write('image\n')
	
	ellipses = input_sources[:,[x_col,y_col,a_col,b_col,theta_col]]
	
	for ellipse in ellipses:
        	fid.write("ellipse(%5.2f,%5.2f, %5.2f, %5.2f, %5.2f)\n"% (ellipse[0]+1+x_offset, ellipse[1]+1+y_offset, ellipse[2]*lengthscale, ellipse[3]*lengthscale, ellipse[4]/np.pi*180))
	
	fid.close()


def create_ds9_region_file(input_sources, region_filename, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10, color='green', width=1, load_text=False,  text_uoffset=25, text_loffset=25, mag_col=2, magerr_col=3):
	'''
	create ds9 region file 
		
	INPUTS:
		input_source: ndarray
		x_col:
		y_col:
		x_offset: only for coordinate type 'image' or 'physical'
		y_offset: only for coordinate type 'image' or 'physical'
		coordinate_type: 'image','physical','fk5'
		radec_deg:
		circle_radius:
		color:
		width:
		text_uoffset:
		text_loffset:
		mag_col:
		magerr_col:

	'''
	if os.path.exists(region_filename):
		if clobber:
			os.remove(region_filename)
		else:
			print "region file %s already exist; quit now under non clobber mode..."%region_filename
			return
		
	fid = open(region_filename,'w')
	fid.write("global color=%s width=%s\n"%(color, width))

        fid.write(coordinate_type+'\n')
	
	xys = input_sources[:,[x_col,y_col]]


	if coordinate_type == 'image' or coordinate_type == 'physical':
		for xy in xys:
        		fid.write("circle(%5.2f,%5.2f,%s)\n"% (xy[0] + x_offset, xy[1] + y_offset, circle_radius))
	elif coordinate_type == 'fk5':
		for xy in xys:
			if radec_deg:
				ra  = xy[0]
				dec = xy[1]
				ra_str,dec_str = radec_format_transformation(ra,dec,mode='float2str',delimiter=':')
			else:
				ra_str  = xy[0]
				dec_str = xy[1]

        		fid.write("circle(%s, %s, %s)\n"% (ra_str, dec_str, circle_radius))
	else:
		raise ValueError("%s not supported..."%coordinate_type)

	
	if load_text:
		xymags = input_sources[:,[x_col, y_col, mag_col, magerr_col]]	
		
		if coordinate_type == 'image' or coordinate_type == 'physical':
			for xymag in xymags:
				x,y,mag,magerr = xymag
				yu = y + text_uoffset
				yl = y - text_loffset
				fid.write("# text(%s, %s) text={%s}\n"%(x,yu,mag))
				fid.write("# text(%s, %s) text={%s}\n"%(x,yl,magerr))
		elif coordinate_type == 'fk5':
				print "loading text for FK5 coordinate is under construction..."
		else:
			raise ValueError("%s not supported..."%coordinate_type)

       	fid.close()

if __name__ == "__main__":
	import sys
	import numpy as np
	insourcefile = sys.argv[1]
	if insourcefile in '--help':
		print "Prepare ds9 source region file"
		print "Usage: %s insourcefile regionfile xcolindex ycolindex"%sys.argv[0]
		sys.exit()
	regionfile = sys.argv[2]
	xcolindex= int(sys.argv[3])
	ycolindex= int(sys.argv[4])
	input_sources = np.loadtxt(insourcefile)
	create_ds9_region_file(input_sources, regionfile, clobber = True,  x_col=xcolindex, y_col=ycolindex)

	#region_filename = sys.argv[2]
	#input_sources = np.loadtxt(insourcefile)
	#create_ds9_region_file_ellipse(input_sources, region_filename, clobber = True,  x_col=0, y_col=1, a_col=2, b_col=3, x_offset=0, y_offset=0,theta_col=4, color='green', width=1)
