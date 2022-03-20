#! /usr/bin/python

import numpy as np
import pyds9
import os
from ds9_display_sources import ds9_show_source_regions, load_source_regions_to_ds9, add_text_to_ds9

from astropy.table import Table

def get_lines_by_selecting_xys_on_image_from_ds9(input_image, output_list,  source_list_given = True, source_list_file=None, source_list_format= 'loadtxt',  x_col = 0, y_col=1, match_criteria = 3, nstar=None, verbose=1):
        '''
	INPUTS:
		input_image:
		output_list: 
		source_list_given:
		source_list_file:
		source_list_format: the method to load the source list, 'loadtxt' ---> np.loadtxt; 'table_daophot' ---> astropy.table format 'daophot'
		x_col: 
		y_col:
		match_criteria: 
		nstar: the number of lines expected to be selected
        '''

        if os.path.exists(input_image):
                print input_image
        else:
                raise IOError("input image doesn't exist!")


	if source_list_given:
		if source_list_file is None:
			raise ValueError("source list file required")

		if not os.path.exists(source_list_file):		
			raise ValueError("%s not exist")
	
		d = ds9_show_source_regions(source_list_file, input_image, x_col=x_col, y_col=y_col, x_offset=0, y_offset=0, coordinate_type = 'image', radec_deg = True, circle_radius = 10)

		source_list_lines = open(source_list_file).readlines()

		if source_list_format == 'loadtxt':
			source_list = np.loadtxt(source_list_file)
			xys = source_list[:,[x_col, y_col]]
		elif source_list_format== 'table_daophot':
			source_list = Table.read(source_list_file, format='daophot')
			xys = np.array([source_list['XCENTER'],source_list['YCENTER']]).transpose()
		else:
			raise ValueError("not supported format for source list: %s"%source_list_format)

	else:
    		d = pyds9.DS9()
       		d.set('fits %s'%input_image)
       		d.set('scale zscale')
       		d.set('zoom to fit')
	
	
	selected_sources_xys = []
	selected_data_rows = []
	selected_row_indexs = []
	i = 1

	continue_pick = True
	while continue_pick:
				
     	   	xys_pick = d.get('imexam coordinate image')
       		x,y = xys_pick.split()
     	 	xy = np.array([float(x),float(y)])
		
        	print "mouse pick (x,y) = (%s, %s)"%(x,y)

		if source_list_given:
#			print xys
			print xy

			exist,xy_match,indice_match = find_corresponding(xys, xy, match_criteria)

     			if exist:
          		      	if len(indice_match)>1:
                       			print "more than one source detected at the SN position... drop this one"
             			else:
					print "corresponding source in provided list found for seleted target!"
					xy = xy_match[0]
					selected_data_row = source_list[indice_match[0]]
					selected_row_index = indice_match[0]
			else:
				print "No corresponding source in provided list found for seleted target..."
		else:
			selected_data_row = xy

		if nstar is None:
			keep_this_default = 'yes'
              	 	keep_this  = raw_input("keep currently selected source: [yes or no; default:%s]:"%keep_this_default) or keep_this_default
		
		else: 
			keep_this = 'yes'		


		keep_this = keep_this.lower()

		if keep_this == 'yes' or keep_this == 'y':
			selected_sources_xys.append(xy)
			selected_data_rows.append(selected_data_row)
			selected_row_indexs.append(selected_row_index)
			
			input_reg = xy
			load_source_regions_to_ds9(d,input_reg,radius=20,color='red',width=1)

			x = xy[0] 
			y = xy[1]
			text = str(i)
			add_text_to_ds9(d,x,y, text, physical_image_offset_x =0, physical_image_offset_y=0,  color='red',width=2,)

			i = i + 1
		if nstar is None:
			continue_select_default = 'yes'
			continue_select  = raw_input("continue select source: [yes or no; default:%s]:"%continue_select_default) or continue_select_default
		else:
			if i > nstar:
				continue_select = 'no'
			else:
				continue_select = 'yes'
			
		continue_select = continue_select.lower()

                if continue_select == 'yes' or continue_select == 'y':
			continue_pick = True
		else:
			continue_pick = False

	if verbose:
		print selected_sources_xys
		print np.array(selected_data_rows)
	
	selected_table = source_list[selected_row_indexs] 	

	if source_list_format == 'table_daophot':
		wanted_lines = range(41) + [indexn+41 for indexn in selected_row_indexs]
		print wanted_lines
		
		if os.path.exists(output_list):
			os.remove(output_list)
	
		fid = open(output_list, 'a')
		for lineindex in wanted_lines:
			fid.write(source_list_lines[lineindex])
		fid.close()
	elif source_list_format == 'loadtxt':
		np.savetxt(output_list, np.array(selected_data_rows))
	else:
		raise ValueError("source list not supported")
	if verbose:
		print selected_table

#	selected_table.write(output_list, format='ascii.commented_header')



def find_corresponding(base_array,target,criteria,verbose=False):
        '''
        find the location of 'target' in 'base_array' when it meet the 'criteria'
        Input:
                base_array: (N,2) array         
                target:     two element array or list
                criteria: scaler, the distance to the target in both x and y direction

        Output:
                exist: True or False
                base_array[ret_mask]: (N,2) shape array and N can be 1
                ret_mask: list, indice of rows meeting requirement
        '''
        if isinstance(target,list):
                target = np.array(target)

        target = target.reshape((1,2))
        diffs = base_array-target
        ret_mask= []
        for i,diff in enumerate(diffs):
            if np.abs(diff[0])<criteria and np.abs(diff[1])<criteria:
                ret_mask.append(i)
        if ret_mask == []:
            exist = False
            return exist,0,0
        else:
            exist = True
            return exist,base_array[ret_mask],ret_mask



if __name__ == "__main__":
	'''
	Want to select a sample points from an image? This will make your life easier!!!
	'''
	import sys
		
	
	
	#if len(sys.argv) == 1:
	#	show_help = True
	#	input_image = ''
	#else:
	#	show_help = False
	#	input_image = sys.argv[1]

	#if input_image == '--help' or show_help:
	#	print "Usage: select_source_on_ds9.py input_image output_list [source_list_file]"
	#	print "where source_list_file is optional; if provided, only sources in the given list can be selected"
	#	sys.exit()	

	#output_list = sys.argv[2]
	#if len(sys.argv) >3:
	#	source_list_file = sys.argv[3]
	#	source_list_given  = True
	#else:
	#	source_list_file = None
	#	source_list_given = False

	import optparse
	parser = optparse.OptionParser()
		
	def_inputimage = ''
	parser.add_option('-i','--input_image', dest = 'input_image', type= 'string', default = def_inputimage, help='input image')

	def_output = 'default'
	parser.add_option('-o', '--output', dest='output', type='string', default=def_output, help='output file; default: %s'%def_output)

	def_sourcelist = ''
	parser.add_option('-s', '--slist', dest='inlist', type='string', default=def_sourcelist, help='input source list; default: %s'%def_sourcelist)	

	def_sourcelist_format = 'loadtxt'
	parser.add_option('-f', '--format', dest='fileformat', type='string', default=def_sourcelist_format, help='input source list format/loading method; default: %s, [loadtxt/table_daophot]'%def_sourcelist_format)	

	def_xcol = 0
	parser.add_option('-x', '--xcol', dest='xcol', type=int, default=def_xcol, help='column number of x coordinate; default: %s'%def_xcol)	

	def_ycol = 1
	parser.add_option('-y', '--ycol', dest='ycol', type=int, default=def_ycol, help='column number of y coordinate; default: %s'%def_ycol)	

	def_matchcriteria = 10
	parser.add_option('-m', '--mc', dest='match_criteria', type=float, default=def_matchcriteria, help='match criteria; default: %s'%def_matchcriteria)	

	def_selnum = None
	parser.add_option('-n','--nstar', dest='nstar', type=int, default=def_selnum, help='number of sources want to select; default:%s'%def_selnum)

	options, remainder = parser.parse_args()

	input_image = options.input_image
	output_list = options.output
	source_list_file = options.inlist
	source_list_format = options.fileformat
	xcol = options.xcol
	ycol = options.ycol
	match_criteria = options.match_criteria
	nstar = options.nstar

	if source_list_file == '':
		source_list_given = False
	else:
		source_list_given = True

	get_lines_by_selecting_xys_on_image_from_ds9(input_image, output_list,  source_list_given = source_list_given, source_list_file= source_list_file, source_list_format=source_list_format, x_col = xcol, y_col=ycol, match_criteria = match_criteria, nstar=nstar)
