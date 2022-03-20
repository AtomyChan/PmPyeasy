#! /usr/bin/python

import os
import sys
import numpy as np
import pyds9

base_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(os.path.join(base_dir,'utils'))
from common import radec_format_transformation

sys.path.append(os.path.join(base_dir,'mpltk'))
from event_handler import mouse_pick


def radec_format_transformation_list(ras, decs, mode='str2float', delimiter=':'):
	'''
	See radec_format_transformation
	'''
	ras_new = []
	decs_new = []
	for ra, dec in zip(ras, decs):
		ra_new, dec_new = radec_format_transformation(ra, dec, mode=mode, delimiter=delimiter)
		ras_new.append(ra_new)
		decs_new.append(dec_new)

	return np.array(ras_new), np.array(decs_new)



def create_ds9_region_file(input_sources, region_filename, clobber = True,  x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg=True, circle_radius = 10, color='green', width=1, load_text=False,  text_uoffset=25, text_loffset=25, textcol1=None, textcol2=None):
	'''
	create ds9 region file

	For the units, x_offset,y_offset,circle_radius,text_uoffset,text_loffset share the same convention:
	in unit of pixel for coordinate type 'image' or 'physical'; in unit of degree for fk5

	INPUTS:
		input_source: ndarray
		x_col: the column for the x coordinate
		y_col: the column for the y coordinate
		x_offset: shift the x coordinate by this value
		y_offset:
		coordinate_type: 'image','physical','fk5'
		radec_deg: RA,Dec coordinate in unit of degree?
		circle_radius:
		color:
		width:
		text_uoffset:
		text_loffset:
		textcol1: the first column of text to be loaded
		textcol2:

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
		if coordinate_type == 'image' or coordinate_type == 'physical':
			textxcol  = xys[:,0]
			textycol1 = xys[:,1] + text_uoffset
			textycol2 = xys[:,1] - text_loffset
		elif coordinate_type == 'fk5':
			if radec_deg:
				racol_degree  = xys[:,0]
				deccol_degree = xys[:,1]
			else:
				racol_degree,deccol_degree = radec_format_transformation_list(xys[:,0], xys[:,1],mode='str2float',delimiter=':')

			textycol1_degree = deccol_degree + text_uoffset
			textycol2_degree = deccol_degree - text_loffset

			textxcol, textycol1 = radec_format_transformation_list(racol_degree, textycol1_degree, mode='float2str',delimiter=':')
			textxcol, textycol2 = radec_format_transformation_list(racol_degree, textycol2_degree, mode='float2str',delimiter=':')
		else:
			raise ValueError("%s not supported..."%coordinate_type)

		if textcol1 is not None:
			for x,y1,text1 in zip(textxcol, textycol1, input_sources[:,textcol1]):
				fid.write("# text(%s, %s) text={%s}\n"%(x,y1, text1))

		if textcol2 is not None:
			for x,y2,text2 in zip(textxcol, textycol2, input_sources[:,textcol2]):
				fid.write("# text(%s, %s) text={%s}\n"%(x,y2, text2))

       	fid.close()


def add_text_to_ds9(d,x,y, text, physical_image_offset_x =0, physical_image_offset_y=0,  color='red',width=2):
        '''
        add 'text' to ds9 frame d at position (x,y)

        Notes:
                Use self.physical_image_offset_x and self.physical_image_offset_y to deal with the difference between image coordinate and physical coordinate
        '''
        x_diff = physical_image_offset_x
        y_diff = physical_image_offset_y

        x_physical = x + x_diff
        y_physical = y + y_diff

        d.set('regions command {text %s %s #text="%s" color="%s" width=%s coor="image"}'%(x_physical,y_physical,text,color,width))


def load_source_regions_to_ds9(ds9_object,input_reg,radius=10,color='green',width=1):
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
                region_set = 'regions command {circle %s %s %s #color=%s width=%s}'%(x,y,radius,color,width)
                print region_set

                ds9_object.set(region_set)
        else:
                try:
                        xy = input_reg.ravel()
                        x = xy[0]
                        y = xy[1]
                        ds9_object.set('regions command {circle %s %s %s #color=%s width=%s}'%(x,y,radius,color,width))
                except:
                        raise IOError("invalid input for ds9 region")




def ds9_show_source_regions(source_file, img, x_col=0, y_col=1, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10):

	input_sources = np.loadtxt(source_file, usecols=(x_col, y_col))
	region_filename = source_file + '.reg'

        create_ds9_region_file(input_sources, region_filename, clobber = True,  x_col=x_col, y_col=y_col, x_offset=x_offset, y_offset=y_offset, coordinate_type = coordinate_type, radec_deg = radec_deg, circle_radius = circle_radius)

        if not os.path.exists(img):
                raise IOError("image %s doesn't exist"%img)

        d = pyds9.DS9()
        d.set('fits %s'%img)
        d.set('zoom to fit')
        d.set('scale zscale')

        load_source_regions_to_ds9(d, region_filename)

	return d



def pylab_talk_to_ds9(xplot, yplot, data, input_image, xcol=0, ycol=1, yerrplot=None):
	'''
	load source list on both python matplotlib.pylab plot and ds9; select point on pylab plot and show on ds9

	idcol, ds9xcol, ds9ycol, pylabxcol, pylabycol, ds9textcol1=None, ds9textcol2=None
	'''
	xs,ys,indexs = mouse_pick(xplot, yplot, yerrplot, yreverse=False)

        d = pyds9.DS9()
        d.set('fits %s'%input_image)
        d.set('zoom to fit')
        d.set('scale zscale')

	for x,y in zip(xs,ys):
		linepicked = data[(xplot==x)*(yplot==y)][0]
		print linepicked
		ds9xy = [linepicked[0], linepicked[1]]
		load_source_regions_to_ds9(d, ds9xy, radius=10,color='green',width=1)




if __name__ == "__main__":
	'''
	Do you want to load a source list to an image on ds9?? Ok! the current task is born for this!!
	'''
	if len(sys.argv) == 1:
		show_help = True
		img = ''
	else:
		show_help = False
		img = sys.argv[1]

	if img == '--help' or show_help:
		print "usage: ds9_display_sources.py image source_list_file"
		print "where source_list_file is the source list file and by default the first and second column will be loaded\n"
		sys.exit()

	source_file = sys.argv[2]
	#ds9_show_source_regions(source_file, img)

	pylab_talk_to_ds9(source_file, img)
