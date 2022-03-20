#! /usr/bin/python

from pyraf import iraf 
from iraf import imcopy

def imcopy_iraf(inputimg, outputimg):
	'''
	INPUTS:
		inputimg: input image
		outputimg: output image
	'''
	imcopy.unlearn()
	imcopy(input=inputimg,output=outputimg)


if __name__ == "__main__":
	import sys

	inimg = sys.argv[1]
	if inimg == '--help':
		print "usage: %s input output"%sys.argv[0]
		sys.exit()
	outimg = sys.argv[2]
	imcopy_iraf(inimg, outimg)

