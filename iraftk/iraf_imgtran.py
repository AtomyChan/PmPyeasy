#! /usr/bin/python

from pyraf import iraf
from iraf import geomap,geoxytran


def geomap_iraf(input_list, output_database, xmin=None, xmax=None, ymin=None, ymax=None, verbose=True, interactive=False):
	'''
	INPUTS:
		input_list: the input coordinate file
		output_database: the output database

	Most input parameters for geomap function use default values here
 	''' 
	
	geomap.unlearn()

	if xmin is not None:
		xmin = xmin
	else:
		xmin = 'INDEF'

	if xmax is not None:
		xmax = xmax
	else:
		xmax = 'INDEF'

	if ymin is not None:
		ymin = ymin
	else:
		ymin = 'INDEF'

	if ymax is not None:
		ymax = ymax
	else:
		ymax = 'INDEF'


	geomap(input=input_list, database=output_database, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, interactive=interactive, verbose=verbose)

	
def geoxytran_iraf(input_coor, output_coor, database, transforms, direction='forward'):
	'''
	INPUT:
		input_coor: input coordinate file to be transformed
		output_coor: output transformed coordiante file
		database: the geomap database file
		transforms: names of coordinate transforms in the database
	'''

	geoxytran.unlearn()
	geoxytran(input=input_coor, output=output_coor, database=database, transforms=transforms, direction=direction)



if __name__ == "__main__":
	import sys
	
	argv1 = sys.argv[1]
	if argv1 == '--help':
		print "usage: thisfile reference_list match_database input_coordinate output_coordinate"
		print "reference_list: transform reference targets position list, formated as x1,y1,x2,y2 \n match_database: the transformation database file \n input_coordinate: the file containing the input coordinate \n output_coordinate: the file containing output coordinate"
		sys.exit()
	
	else:
		reference_list = argv1

	match_database = sys.argv[2]

	input_coordinate  = sys.argv[3]
	output_coordinate = sys.argv[4]

	geomap_iraf(reference_list, match_database)
	geoxytran_iraf(input_coordinate, output_coordinate, match_database, reference_list)


#	print 

