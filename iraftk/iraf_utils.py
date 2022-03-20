#! /anaconda/bin/python

import numpy as np
import os
import shutil
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


def iraf_get_default_filename(image, mode, checkdir=None):
	'''
	figure out the directory and get the filename  for the default setup
	'''
	image_segs = image.split('/')

	if checkdir is None:
		if len(image_segs) == 1:
			imagedir = './'
		else:
			if image_segs[0] not in ['','.']:
				imagedir = ''
			elif image_segs[0] == '':
				imagedir = '/'
				image_segs.pop(0)
			else:
				imagedir = './'
				image_segs.pop(0)

			if len(image_segs)>1:
				for dirseg in image_segs[:-1]:
					#print dirseg
					imagedir = os.path.join(imagedir, dirseg)
		checkdir = imagedir

	root = image_segs[-1]
	output_base = root + '.'+ mode
	outputs_old = [f for f in os.listdir(checkdir) if f.split('.')[-1].isdigit() and output_base=='.'.join(f.split('.')[:-1])]
	outputs_ov =  np.sort([int(f.split('.')[-1]) for f in outputs_old]) #old output version
	if len(outputs_ov)>0:
		ver = outputs_ov[-1]+1
	else:
		ver = 1

	output = output_base+'.'+str(ver)

	return output


def iraf_coo_file_id_filter(infile, output, idlist):
    '''
    select lines in infile xxx.coo.X file with ID in idlist
    '''
    fid = open(output, 'a')

    for line in open(infile,'r').readlines():
        if line[0] == '#':
            fid.write(line)
        else:
            id = int(line.split()[-1])
            if id in idlist:
                fid.write(line)
    fid.close()

def iraf_mag_file_id_filter(infile, output, idlist, firstfirst=True):
    '''
    select lines in infile xxx.mag.X file with ID in idlist
    '''
    iraf_id_filter(infile, output, idlist, linegroup=5, idindex=3, idlist_first_as_output_first=firstfirst)


def iraf_id_filter(infile, output, idlist, linegroup, idindex, idlist_first_as_output_first=True):
	'''
	INPUTS:
		infile:
		output:
		idlist:
		linegroup: how many lines for the same object
		idindex: the index of the ID entry in the row containing ID
	'''
	if os.path.exists(output):
		if os.path.samefile(infile, output):
			infile_copy = infile+'.copy'
			shutil.copy(infile, infile_copy)
			infile = infile_copy
		os.remove(output)
	fid = open(output, 'a')
	if idlist_first_as_output_first:
		idfirst = [idlist[0]]
		fidtemp = open(infile,'r')
		iraf_id_filter_kit(fidtemp, fid, idfirst, linegroup, idindex, write_headline=True)
		fidtemp.close()

		idlist = idlist[1:]
		infid = open(infile,'r')
		iraf_id_filter_kit(infid, fid, idlist, linegroup, idindex, write_headline=False)
		infid.close()
	else:
		infid = open(infile,'r')
		iraf_id_filter_kit(infid, fid, idlist, linegroup, idindex, write_headline=True)
		infid.close()

	fid.close()



def iraf_id_filter_kit(infid, outputfid, idlist, linegroup, idindex, write_headline=True):
	Nheadlines = 0
	yesnextline = 0
	for i,line in enumerate(infid.readlines()):
		if line[0] == '#':
			if write_headline:
				outputfid.write(line)
			Nheadlines += 1
		else:
			if np.mod(i-Nheadlines, linegroup)==0:
				id = int(line.split()[idindex])
				if id in idlist:
					outputfid.write(line)
					yesnextline = 1
				else:
					yesnextline = 0
			else:
				if yesnextline:
					outputfid.write(line)
