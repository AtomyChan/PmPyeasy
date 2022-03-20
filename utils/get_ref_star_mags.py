#! /usr/bin/python

from astropy.table import Table

import os
import ASASSN_photometry
sn = ASASSN_photometry.photometry('2017cbv', 'LCOGT', photometry_method='apphot')


# 2017cbv template 017.fits
xys = {}
xys['star1'] = [1751,1802]
xys['star2'] = [1873.5,1915]
xys['star3'] = [1747.5,1732]
xys['star3'] = [1982,1802]
xys['star3'] = [1747.5,1732]
xys['star4'] = [1982,1802]
xys['star5'] = [2304,1830]
xys['star6'] = [2621,2044]
xys['star7'] = [2548.5,2223]
xys['star8'] = [2580,2307]
xys['star9'] = [2135,2105]
xys['star10'] = [1474.5,2073]
xys['star11'] = [1443.5,2248]

save_dir = '/home/asassn/ASASSN/LCOGT/2017cbv/workplace'

for star in xys.keys():
	print star
	xy = xys[star]
	x = xy[0]
	y = xy[1]
	out = sn.get_light_curve_for_control_star('B',x,y)

	filename = os.path.join(save_dir, star+'_B_mags.txt')

	print out

	out.write(filename,format='ascii')

#plt.plot(star1['obstime'], star1['mag'],'o')
#plt.show()
