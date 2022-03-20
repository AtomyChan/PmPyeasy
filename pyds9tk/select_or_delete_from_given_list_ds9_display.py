#! /usr/bin/python

import os
import pyds9
import numpy as np
from astropy.table import Table


def delete_entry_from_given_table_on_ds9_display(fitsimage, inputfile_table, table_format= 'ascii.fast_no_header', x_col =0, y_col = 1, coordinate_system = 'image', ):
	'''
	INPUTS:
		fitsimage:
		inputtable:
		table_format: common ones 'ascii', 'ascii.commented_header', 'ascii.csv', 'ascii.fixed_width', 
		x_col:
		y_col:
		coordinate_system:

	'''
	intable = Table.read(inputfile_table, format=table_format)

	xy_colnames = [intable.colnames[x_col],intable.colnames[y_col]]
	xy_cols = intable[xy_colnames]


	


def input_xys_through_ds9_get_mags(fitsimg,xysmagsfile,criteria):
        '''
        pick up stars from ds9 display and extract mags from file xysmagfile
                
        INPUTS:
                fitsimg:
                xysmagsfile:
        '''

        xysmags = np.loadtxt(xysmagsfile)
        xys = xysmags[:,0:2]

        print xys

        i = 0

        num_flags = []
        xs = []
        ys = []
        mags = []
        magerrs = []

        while True:
                num_default = i
                num_flag = raw_input("please input the num flag for the star (default:%s )"%num_default) or num_default

                num_int = int(num_flag)
                i = num_int + 1

                if num_int < 0:
                        break


                xy = get_xy_on_image_from_ds9(input_image=fitsimg)
                print "xy from mouse pick up:", xy


                x = xy[0] - 1
                y = xy[1] - 1

                xy_target = [x,y]

                yesfound, xy_match,index = find_corresponding(xys, xy_target, criteria)


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
                        mag = xysmags[index,2]
                        magerr = xysmags[index,3]

                num_flags.append(num_int)
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

        print out_table

        return out_table



def get_xy_on_image_from_ds9(input_image):
        '''
	
	mouse pick from ds9 display
        '''
        if os.path.exists(input_image):
                print input_image
        else:
                raise IOError("input image doesn't exist!")


        d = pyds9.DS9()
        d.set('fits %s'%input_image)
        d.set('scale zscale')
        d.set('zoom to fit')
        xys = d.get('imexam coordinate image')
        x,y = xys.split()
        xy = np.array([float(x),float(y)])
        print xy

        return xy




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
