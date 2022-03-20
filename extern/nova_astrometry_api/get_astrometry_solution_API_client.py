#! /usr/bin/python

import os
import sys
import subprocess, threading


class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            print 'Thread started'
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            print 'Thread finished'

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print 'Terminating process'
            self.process.terminate()
            thread.join()

	returncode = self.process.returncode
        print returncode

	return returncode


def get_wcs_solution(img, wcs_out, ra, dec, radius, time_out_max = 300, apikey ='aubtqramhfmriufp' ):
	'''
	Get astrometry for image 'img' with nova.astrometry.net API and the wcs result will be saved in 'wcs_out'
	Please see client.py for details

	INPUTS:
		img:
		wcs_out:
		ra: in degree
		dec: in degree
		radius: the estimate for the field radius
	'''		

	API_client = "./nova_astrometry_api/client.py"
	if not os.path.exists(API_client):
		raise IOError("API client for astrometry not available")

	print "working on astrometry for %s"%img

	command_line = "python %s --server http://nova.astrometry.net/api/ --apikey %s --wcs %s --upload %s --ra %s --dec %s --radius %s -p >/dev/null"%(API_client, apikey,  wcs_out, img, ra, dec, radius)

	print command_line

	command = Command(command_line)

	failure_signal = command.run(timeout=time_out_max)
	#failure_signal  non-zero value means failure

	if not failure_signal:
		success = True
	else:
		success = False

	return success
	


if __name__ == "__main__":
	

	import optparse

	parser = optparse.OptionParser()

	def_inimg = ''
	parser.add_option('-i', '--input_image', dest='in_img',  type="string", default=def_inimg, help="input image for astrometry solution")

	def_outwcs = 'astrometry_result.wcs'
	parser.add_option('-w', '--wcs', dest='wcs',  type="string", default=def_outwcs, help="astrometry solution, the filename containing the wcs solution")


	def_maxtime = 300
	parser.add_option('-m', '--maxtime', dest='maxtime', type=float, default=def_maxtime, help="maximum time on the task")

	def_newimage = False
	parser.add_option('-n', '--newimage', dest='newimage', action="store_true", default=def_newimage, help="whether want the new image with wcs replaced for the input image")

	def_outimg = 'new_image_astrometry.fits'
	parser.add_option('-o', '--output_image', dest='out_image',  type="string", default=def_outimg, help="output image with new astrometry solution")
	

	options, remainder = parser.parse_args()

	input_image = options.in_img
	wcs_out     = options.wcs
	maxtime     = options.maxtime
	newimage    = options.newimage
	output_image = options.out_image

	if os.path.exists(wcs_out):
		os.remove(wcs_out)

	if os.path.exists(output_image):
		os.remove(output_image)
	

	if input_image == '': 
		raise ValueError("Please specify the image though -i --input_image")

	if not os.path.exists(input_image):
		raise ValueError('the input image %s not exists'%input_image)


	

	failure = get_astrometry_solution_API_client(input_image, wcs_out, time_out_max = maxtime, )
	print failure

	if failure:
		print "Sorry dear friend, we can't get the astrometry solution for you... maybe more time can help"
	else:
		if newimage:
			#easy way to go: replace the whole header with wcs information 
		
			from astropy.io import fits
			wcs_hdu  =  fits.open(wcs_out)
			wcs_hdr = wcs_hdu[0].header
		
			hdu = fits.open(input_image)
			hdu[0].header = wcs_hdr
	
			hdu.writeto(output_image)
	
			hdu.close()

