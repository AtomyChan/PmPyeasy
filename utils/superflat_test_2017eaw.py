import ASASSN_photometry
import os
import numpy as np
from superflat_calibrate import superflat_calibration

sn = ASASSN_photometry.photometry('2017eaw', 'RP', photometry_method='psfphot')
sn._photometry__dict2table()
#V = sn._photometry__select_rows_from_table(sn.result_table_newsaved,'flt','V')
#rp = sn._photometry__select_rows_from_table(sn.result_table_newsaved,'flt','rp')
ip = sn._photometry__select_rows_from_table(sn.result_table_newsaved,'flt','ip')


#for img in V['name']:
#	if sn.photometry_info[img]['realimg'][-8:-6]=='NM':
#		if sn.photometry_info[img]['drop'] == 0 and sn.photometry_info[img]['template'] == 0:
#			print img, sn.photometry_info[img]['realimg']
#			calfile = os.path.join(sn.psf_photometry_dir, img.split('.')[0]+'_tpl.match')
#			data = np.loadtxt(calfile)
#			offset_sn = superflat_calibration(data[:,0],data[:,1],data[:,6], data[:,2], data[:,7], data[:,3], 1031, 1192, model='poly2', param0=None, display=False)
#			offset_old = sn.photometry_info[img]['relmag'] - sn.photometry_info[img]['instmag']
#			print offset_old, offset_sn
#			sn.photometry_info[img]['relmag'] = sn.photometry_info[img]['instmag'] + offset_sn
#			sn.photometry_info[img]['relmagerr'] = 0.01
#
#sn.save_results()


#for img in rp['name']:
#	if sn.photometry_info[img]['realimg'][-8:-6]=='NM':
#		if sn.photometry_info[img]['drop'] == 0 and sn.photometry_info[img]['template'] == 0:
#			print img, sn.photometry_info[img]['realimg']
#			calfile = os.path.join(sn.psf_photometry_dir, img.split('.')[0]+'_tpl.match')
#			data = np.loadtxt(calfile)
#			offset_sn = superflat_calibration(data[:,0],data[:,1],data[:,6], data[:,2], data[:,7], data[:,3], 1031, 1192, model='poly2', param0=None, display=False)
#			offset_old = sn.photometry_info[img]['relmag'] - sn.photometry_info[img]['instmag']
#			print offset_old, offset_sn, offset_sn - offset_old
#			sn.photometry_info[img]['relmag'] = sn.photometry_info[img]['instmag'] + offset_sn
#			sn.photometry_info[img]['relmagerr'] = 0.01
#
#sn.save_results()


for img in ip['name']:
	if sn.photometry_info[img]['realimg'][-8:-6]=='NM':
		if sn.photometry_info[img]['drop'] == 0 and sn.photometry_info[img]['template'] == 0:
			print img, sn.photometry_info[img]['realimg']
			calfile = os.path.join(sn.psf_photometry_dir, img.split('.')[0]+'_tpl.match')
			data = np.loadtxt(calfile)
			offset_sn = superflat_calibration(data[:,0],data[:,1],data[:,6], data[:,2], data[:,7], data[:,3], 1031, 1192, model='poly2', param0=None, display=False)
			offset_old = sn.photometry_info[img]['relmag'] - sn.photometry_info[img]['instmag']
			print offset_old, offset_sn, offset_sn - offset_old
			sn.photometry_info[img]['relmag'] = sn.photometry_info[img]['instmag'] + offset_sn
			sn.photometry_info[img]['relmagerr'] = 0.01

sn.save_results()
