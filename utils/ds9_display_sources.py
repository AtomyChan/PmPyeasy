#! /usr/bin/python

import os
import sys
#sys.path.append('/home/asassn/ASASSN/photometry_pipeline/')
from common import create_ds9_region_file
import numpy as np

import pyds9


def add_text_to_ds9(d,x,y, text, physical_image_offset_x =0, physical_image_offset_y=0,  color='red',width=2,):
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

	ds9_show_source_regions(source_file, img)
