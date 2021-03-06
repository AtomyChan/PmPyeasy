from common import __get_fits_info
import os
from astropy.time import Time

#get image information
#

def get_obstime(img,obstime_key,telescope):

	if telescope == 'LCOGT':
		obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'LCOGT2m':
		obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'LCOGT0m4':
		obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'JB':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'JG':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'KA':
		obstime = get_obstime_OAUV(img,obstime_key)
	elif telescope == 'NJ':
		obstime = get_obstime_OAUV(img,obstime_key)
	elif telescope == 'GC':
		obstime = get_obstime_GC(img,obstime_key)
	elif telescope == 'OAUV':
		obstime = get_obstime_OAUV(img,obstime_key)
	elif telescope == 'RK':
        	obstime = get_obstime_RK(img,obstime_key)
	elif telescope == 'TV':
        	obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'WFCCD':
		obstime = get_obstime_RK(img,obstime_key)
	elif telescope == 'SITe2K':
		obstime = get_obstime_RK(img,obstime_key)
	elif telescope == 'LDSS3C':
		obstime = get_obstime_RK(img,obstime_key)
	elif telescope == 'SDSS':
		obstime = get_obstime_RK(img,obstime_key)
	elif telescope == 'Iowa':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'MDM':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'IMACS':
        	obstime = get_obstime_RK(img,obstime_key)
	elif telescope == 'weihai':
        	obstime = get_obstime_weihai(img,obstime_key)
	elif telescope == 'LTH2RG':
        	obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'OSIRIS':
        	obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'PanSTARRS':
		obstime = get_obstime_MJD(img, obstime_key)
	elif telescope == 'TJO':
                obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'LT':
                obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'DEMONEXT':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'SMARTS':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'SMARTS-IR':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'REM':
		obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'REM-NIR':
		obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'GS':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'RP':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'JN':
		obstime = get_obstime_JN(img, obstime_key)
	elif telescope == 'ALFOSC':
		obstime = get_obstime_OAUV(img, obstime_key)
	elif telescope == 'MODS1R':
		obstime = get_obstime_MJD(img, obstime_key)
	elif telescope == 'NOTCAM':
		obstime = get_obstime_OAUV(img, obstime_key)
	elif telescope == 'ARIES130cm':
		obstime = get_obstime_OAUV(img, obstime_key)
	elif telescope == 'EFOSC2':
		obstime = get_obstime_MJD(img, obstime_key)
	elif telescope == 'LOTIS':
		obstime = get_obstime_OAUV(img, obstime_key)
	elif telescope == 'SWIFTtest':
		obstime = get_obstime_OAUV(img, obstime_key, extension=1)
	elif telescope == 'IRIDA':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'Rozhen50':
		obstime = get_obstime_JN(img,obstime_key)
	elif telescope == 'Rozhen2m':
		obstime = get_obstime_JN(img,obstime_key)
	elif telescope == 'StuParker':
		obstime = get_obstime_JN(img,obstime_key)
	elif telescope == 'DT':
		obstime = get_obstime_STRAIGHT(img,obstime_key)
	elif telescope == 'WFCAM':
                obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'UFTI':
                obstime = get_obstime_MJD(img,obstime_key)
	elif telescope == 'FORS2':
                obstime = get_obstime_MJD(img,obstime_key)
	else:
		raise KeyError("Telescope %s is not supported yet"%telescope)

	return obstime





def get_obstime_OAUV(img,obstime_key, extension=None):

	fitsinfo = __get_fits_info(img,obstime_key, extension=extension)
	tt = Time(fitsinfo[obstime_key])

	obstime = tt.jd

	return obstime

def get_obstime_weihai(img, obstime_key):

	key1 = obstime_key.split()[0]
	key2 = obstime_key.split()[1]
	key3 = obstime_key.split()[2]

	fitsinfo = __get_fits_info(img, obstime_key.split())

	if fitsinfo[key1] is not None:
		obstime = fitsinfo[key1]
	else:
		tt = Time(fitsinfo[key2] + 'T' + fitsinfo[key3])
        	obstime = tt.jd

	return obstime

def get_obstime_GC(img, obstime_key):
	'''
	For GC images, if the image is co-added image then MIDPOINT exists which is the at the middle of the exposure and will be used as observation time
	otherwise the JD value in the header will be used directly
	'''
	key1 = obstime_key.split()[0]
	key2 = obstime_key.split()[1]

	fitsinfo = __get_fits_info(img,obstime_key.split())

	v1 = fitsinfo[key1]
	v2 = fitsinfo[key2]

	if v1 is None:
		obstime = v2
	else:
		tt = Time(v1)
		obstime = tt.jd

	return obstime


def get_obstime_RK(img,obstime_key):

	key1 = obstime_key.split()[0]
	key2 = obstime_key.split()[1]

	fitsinfo = __get_fits_info(img,obstime_key.split())
	tt = Time(fitsinfo[key1] + 'T' + fitsinfo[key2])

	obstime = tt.jd

	return obstime

def get_obstime_JN(img, obstime_key):

	fitsinfo = __get_fits_info(img, obstime_key)
	tt = Time(fitsinfo[obstime_key])
	obstime  = tt.jd

	return obstime




def get_obstime_MJD(img,obstime_key):

	fitsinfo = __get_fits_info(img,obstime_key)
	MJD = fitsinfo[obstime_key]

	obstime = MJD + 2400000.5

	return obstime

def get_obstime_STRAIGHT(img,obstime_key):

	fitsinfo = __get_fits_info(img,obstime_key)
	obstime = fitsinfo[obstime_key]

	return obstime





def get_exptime(img,exptime_key,telescope):

	if telescope == 'LCOGT':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'LCOGT2m':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'LCOGT0m4':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'JB':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'JG':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'GC':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'KA':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'NJ':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'OAUV':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'RK':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'TV':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'WFCCD':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'SITe2K':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'LDSS3C':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'SDSS':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'Iowa':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'MDM':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'IMACS':
        	exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'weihai':
        	exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'LTH2RG':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'OSIRIS':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'TJO':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'LT':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'DEMONEXT':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'SMARTS':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'SMARTS-IR':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'REM':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'REM-NIR':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'GS':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'RP':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'JN':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'ALFOSC':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'MODS1R':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'PanSTARRS':
		exptime = get_exptime_STRAIGHT(img,exptime_key)
	elif telescope == 'NOTCAM':
		exptime = get_exptime_NOTCAM(img, exptime_key)
	elif telescope == 'ARIES130cm':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'EFOSC2':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'LOTIS':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'SWIFTtest':
		exptime = get_exptime_STRAIGHT(img, exptime_key, extension=1)
	elif telescope == 'IRIDA':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'Rozhen50':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'Rozhen2m':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'StuParker':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'DT':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'WFCAM':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'UFTI':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	elif telescope == 'FORS2':
		exptime = get_exptime_STRAIGHT(img, exptime_key)
	else:
		raise KeyError("Telescope %s is not supported yet"%telescope)

	return exptime




def get_exptime_STRAIGHT(img,exptime_key, extension=None):

	fitsinfo = __get_fits_info(img,exptime_key, extension=extension)
	exptime = fitsinfo[exptime_key]

	return exptime


def get_exptime_NOTCAM(img, exptime_key):

	key1 = exptime_key.split()[0]
	key2 = exptime_key.split()[1]

	fitsinfo = __get_fits_info(img, exptime_key.split())

	sexptime = fitsinfo[key1]
	ncombine = fitsinfo[key2]

	exptime = sexptime * ncombine

	return exptime

def get_exptime_WFCAM_total(img, exptime_keys):
	'''
	Final exptime = EXP_TIME*NEXP*NJITTER*NUSTEP
	'''
	key1 = exptime_keys.split()[0]
	key2 = exptime_keys.split()[1]
	key3 = exptime_keys.split()[2]
	key4 = exptime_keys.split()[3]

	fitsinfo = __get_fits_info(img, exptime_keys.split())

	sexptime = fitsinfo[key1]
	nexp = fitsinfo[key2]
	njitter = fitsinfo[key3]
	nustep = fitsinfo[key4]

	totalexptime = sexptime * nexp * njitter*nustep

	return totalexptime



def get_filter(img,fltkey,telescope):

	if telescope == 'LCOGT':
		flt = get_filter_LCOGT(img,fltkey)
	elif telescope == 'LCOGT2m':
		flt = get_filter_LCOGT(img,fltkey)
	elif telescope == 'LCOGT0m4':
		flt = get_filter_LCOGT(img,fltkey)
	elif telescope == 'JB':
		flt = get_filter_JB(img,fltkey)
	elif telescope == 'GC':
		flt = get_filter_GC(img,fltkey)
	elif telescope == 'KA':
		flt = get_filter_KA(img, fltkey)
	elif telescope == 'NJ':
		flt = get_filter_NJ(img,fltkey)
	elif telescope == 'JG':
		flt = get_filter_JG(img,fltkey)
	elif telescope == 'OAUV':
		flt = get_filter_OAUV(img,fltkey)
	elif telescope == 'RK':
		flt = get_filter_RK(img,fltkey)
	elif telescope == 'TV':
		flt = get_filter_TV(img,fltkey)
	elif telescope == 'WFCCD':
		flt = get_filter_WFCCD(img,fltkey)
	elif telescope == 'SITe2K':
		flt = get_filter_SITe2K(img,fltkey)
	elif telescope == 'LDSS3C':
		flt = get_filter_LDSS3C(img,fltkey)
	elif telescope == 'SDSS':
		flt = get_filter_SDSS(img,fltkey)
	elif telescope == 'Iowa':
		flt = get_filter_iowa(img,fltkey)
	elif telescope == 'MDM':
		flt = get_filter_MDM(img,fltkey)
	elif telescope == 'IMACS':
        	flt = get_filter_IMACS(img,fltkey)
	elif telescope == 'weihai':
        	flt = get_filter_weihai(img,fltkey)
	elif telescope == 'LTH2RG':
        	flt = get_filter_LTH2RG(img,fltkey)
	elif telescope == 'OSIRIS':
        	flt = get_filter_OSIRIS(img,fltkey)
	elif telescope == 'TJO':
        	flt = get_filter_TJO(img,fltkey)
	elif telescope == 'LT':
        	flt = get_filter_LT(img,fltkey)
	elif telescope == 'DEMONEXT':
		flt = get_filter_DEMONEXT(img,fltkey)
	elif telescope == 'SMARTS':
		flt = get_filter_SMART(img,fltkey)
	elif telescope == 'SMARTS-IR':
		flt = get_filter_SMART_IR(img,fltkey)
	elif telescope == 'REM':
		flt = get_filter_REM(img,fltkey)
	elif telescope == 'REM-NIR':
		flt = get_filter_REM(img,fltkey)
	elif telescope == 'GS':
		flt = get_filter_GS(img, fltkey)
	elif telescope == 'RP':
		flt = get_filter_RP(img, fltkey)
	elif telescope == 'JN':
		flt = get_filter_JN(img, fltkey)
	elif telescope == 'ALFOSC':
		flt = get_filter_ALFOSC(img, fltkey)
	elif telescope == 'MODS1R':
		flt = get_filter_MODS1R(img,fltkey)
	elif telescope == 'PanSTARRS':
		flt = get_filter_PanSTARRS(img,fltkey)
	elif telescope == 'NOTCAM':
		flt = get_filter_NOTCAM(img, fltkey)
	elif telescope == 'ARIES130cm':
		flt = get_filter_ARIES130cm(img, fltkey)
	elif telescope == 'EFOSC2':
		flt = get_filter_EFOSC2(img, fltkey)
	elif telescope == 'LOTIS':
		flt = get_filter_LOTIS(img, fltkey)
	elif telescope == 'SWIFTtest':
		flt = get_filter_SWIFTtest(img, fltkey, extension=1)
	elif telescope == 'IRIDA':
		flt = get_filter_IRIDA(img, fltkey)
	elif telescope == 'Rozhen50':
		flt = get_filter_Rozhen50(img, fltkey)
	elif telescope == 'Rozhen2m':
		flt = get_filter_Rozhen2m(img, fltkey)
	elif telescope == 'StuParker':
		flt = get_filter_StuParker(img, fltkey)
	elif telescope == 'DT':
		flt = get_filter_DT(img, fltkey)
	elif telescope == 'WFCAM':
		flt = get_filter_WFCAM(img, fltkey)
	elif telescope == 'UFTI':
		flt = get_filter_UFTI(img, fltkey)
	elif telescope == 'FORS2':
		flt = get_filter_FORS2(img, fltkey)
	else:
		raise KeyError("Telescope %s is not supported yet"%telescope)

	return flt

def get_filter_DT(img, fltkey):
	fitsinfo = __get_fits_info(img, fltkey)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	return flt

def get_filter_UFTI(img, fltkey):
	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')

	if flt == 'J98':
		flt = 'J'

	if flt == 'H98':
		flt = 'H'

	if flt == 'K98':
		flt = 'K'

	return flt

def get_filter_WFCAM(img, fltkey):
	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	return flt

def get_filter_StuParker(img, fltkey):
	flt = 'C'
	return flt

def get_filter_FORS2(img, fltkey):
	flt = 'V'
	return flt


def get_filter_Rozhen50(img, fltkey, extension=None):
	fitsinfo = __get_fits_info(img, fltkey, extension=extension)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	return flt


def get_filter_Rozhen2m(img, fltkey, extension=None):
	fitsinfo = __get_fits_info(img, fltkey, extension=extension)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')

	if flt == 'Sloan_g':
		flt = 'gp'
	if flt =='Sloan_r':
		flt = 'rp'
	if flt == 'Sloan_i':
		flt = 'ip'
	if flt =='Sloan_z':
		flt = 'zp'


	return flt


def get_filter_IRIDA(img, fltkey):
	fitsinfo = __get_fits_info(img, fltkey)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	if flt == 'r':
		flt = 'rp'
	if flt =='i':
		flt = 'ip'

	return flt

def get_filter_TV(img, fltkey):
	fitsinfo = __get_fits_info(img, fltkey)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	if flt == 'Clear':
		flt = 'C'
	return flt

def get_filter_SWIFTtest(img, fltkey, extension=None):
	fitsinfo = __get_fits_info(img, fltkey, extension=extension)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	return flt



def get_filter_LOTIS(img, fltkey, extension=None):
	fitsinfo = __get_fits_info(img, fltkey, extension=extension)

	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	return flt



def get_filter_EFOSC2(img,fltkey):

	fltkey = 'HIERARCH ESO INS FILT1 NAME'
	fitsinfo = __get_fits_info(img, fltkey)

	flt = fitsinfo['HIERARCH ESO INS FILT1 NAME']

	flt = flt.split('#')[0]

	return flt







def get_filter_ARIES130cm(img, fltkey):

	img_real = os.path.realpath(img)

	flt_seg = img_real.split('_')[-2]
	flt = flt_seg[0]

	if flt == 'u':
		flt = 'U'
	elif flt == 'b':
		flt = 'B'
	elif flt == 'v':
		flt = 'V'
	elif flt == 'r':
		flt = 'R'
	elif flt == 'i':
		flt = 'I'
	else:
		raise ValueError("filter information not obtained for %s"%img)


	return flt

def get_filter_LCOGT(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

	if flt == 'zs':
		flt = 'zp'

	return flt


def get_filter_JB(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()
	return flt


def get_filter_GC(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()

	if flt == 'OG530':
		flt = 'R'

	flt = flt.replace("'",'p')

	return flt

def get_filter_KA(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	#flt = flt.strip()
	return flt

def get_filter_NJ(img, fltkey):
	flt = 'C'
	return flt


def get_filter_JG(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')
	return flt


def get_filter_OAUV(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()

	if flt == 'Johnson_B':
		flt = 'B'
	elif flt == 'Johnson_V':
		flt = 'V'
	elif flt == 'SDSS_r':
		flt = 'rp'
	elif flt == 'SDSS_i':
		flt = 'ip'
	else:
		raise KeyError("unrecongnized filter, please check...")

	return flt



def get_filter_RK(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()
	return flt

def get_filter_WFCCD(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()
	return flt

def get_filter_SITe2K(img, fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()




	if flt[0]== "R":
		flt = "R"
	elif flt[0] == 'B':
		flt = 'B'
	elif flt[0] == 'V':
		flt = 'V'
	elif flt[0] == 'r':
		flt = 'rp'
	elif flt[0] == 'i':
		flt = 'ip'
	else:
		raise ValueError("filter not understood... please check")

	return flt




def get_filter_SDSS(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.strip()

	if flt == 'u':
		flt = 'up'

	if flt == 'g':
		flt = 'gp'


	if flt == 'r':
		flt = 'rp'

	if flt == 'i':
		flt = 'ip'

	if flt == 'z':
		flt = 'zp'

	return flt


def get_filter_iowa(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

#	print flt

	flt = flt.replace(" ","")

#	print flt

	if flt == 'B-Blue' or flt == 'B-J-Cblue':
		flt = 'B'
	elif flt == 'G-Sloang' or flt=='G-Sloan_g':
		flt = 'gp'
	elif flt == 'G-Green' or flt == 'V-Visual' or flt == 'V-AstrodonVisual':
		flt = 'V'
	elif flt == 'R-Red' or flt == 'R-Sloan_r':
		flt = 'rp'
	elif flt == 'W-Longpass' or flt == 'W-LongpassIR' or flt == 'I-Sloan_i':
		flt = 'ip'
	elif flt == "G-GRISM":
		flt = 'V'
	elif flt == 'X-Astrodon_Red':
		flt = 'R'
	elif flt == "N-None":
		flt = 'N'

	elif flt == "T-Transmissiongrism":
		flt = 'T'

	else:
		print img
		raise KeyError("the filter information in fits header for iowa telescope is not understood, please check...")

	return flt

def get_filter_PanSTARRS(img, fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(" ","")

	if flt[0] == 'g':
		flt = 'gp'
	elif flt[0] == 'r':
		flt = 'rp'
	elif flt[0] == 'i':
		flt = 'ip'
	elif flt[0] == 'z':
		flt = 'zp'
	elif flt[0] == 'y':
		flt = 'yp'
	else:
		print img
		raise KeyError("the filter information in fits header for iowa telescope is not understood, please check...")

	return flt

def get_filter_MODS1R(img, fltkey):
	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(" ","")

	if flt == 'r_sdss':
		flt = 'rp'

	return flt

def get_filter_MDM(img,fltkey):


	fltkey1 = fltkey.split()[0]
	fltkey2 = fltkey.split()[1]


	fitsinfo = __get_fits_info(img,fltkey.split())

	flt1 = fitsinfo[fltkey1]
	flt1 = flt1.replace(' ','')

	flt2 = fitsinfo[fltkey2]
	flt2 = flt2.replace(' ','')

	if flt1 == 'Open':
		flt = flt2
	else:
		flt = flt1


	if flt == 'SDSS-r':
		flt = 'rp'

	if flt == 'SDSS-i':
		flt = 'ip'


	return flt


def get_filter_NOTCAM(img, fltkey):

	fltkey1 = fltkey.split()[0]
	fltkey2 = fltkey.split()[1]

	fitsinfo = __get_fits_info(img,fltkey.split())

	flt1 = fitsinfo[fltkey1]
	flt1 = flt1.replace(' ','')

	flt2 = fitsinfo[fltkey2]
	flt2 = flt2.replace(' ','')


	if flt1 == 'Open':
		flt = flt2
	else:
		flt = flt1

	print flt

	if flt[0] == 'J':
		flt = 'J'
	elif flt[0] == 'H':
		flt = 'H'
	elif flt[0] == 'K' or flt[0] == 'Ks':
		flt = 'K'
	else:
		raise ValueError("input filter %s not suppported"%flt)

	return flt




def get_filter_ALFOSC(img, fltkey):

	print img
	fltkey1 = fltkey.split()[0]
	fltkey2 = fltkey.split()[1]


	fitsinfo = __get_fits_info(img,fltkey.split())

	flt1 = fitsinfo[fltkey1]
	flt1 = flt1.replace(' ','')

	flt2 = fitsinfo[fltkey2]
	flt2 = flt2.replace(' ','')


	if flt1 == 'Open':
		flt = flt2
	else:
		flt = flt1

	print flt

	if flt[0] == 'B':
		flt = 'B'
	elif flt[0] == 'V':
		flt = 'V'
	elif flt[0] == 'u':
		flt = 'up'
	elif flt[0] == 'g':
		flt = 'gp'
	elif flt[0] == 'r':
		flt = 'rp'
	elif flt[0] == 'i':
		flt = 'ip'
	elif flt[0] == 'z':
		flt = 'zp'
	else:
		raise ValueError("input filter %s not suppported"%flt)


	return flt



def get_filter_RP(img, fltkey):

	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ", "")

	if flt == 'SG':
		flt = 'gp'

	if flt == 'SR':
		flt = 'rp'

	if flt == 'SI':
		flt = 'ip'

	return flt


def get_filter_JN(img, fltkey):

	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ", "")

	if flt == 'FiltreB':
		flt = 'B'

	if flt == 'FiltreV':
		flt = 'V'

	return flt



def get_filter_OSIRIS(img, fltkey):

	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ", "")

	if flt == 'Sloan_g':
		flt = 'gp'

	if flt == 'Sloan_r':
		flt = 'rp'

	if flt == 'Sloan_i':
		flt = 'ip'

	return flt

def get_filter_IMACS(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

	if flt == 'Bessell_V2':
		flt = 'V'
	if flt == 'g_Sloan':
		flt = 'gp'

	if flt == 'Sloan_r':
		flt = 'rp'

	flt = flt.strip()

	return flt



def get_filter_LDSS3C(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]
	flt = flt.replace(' ','')

	if flt == 'r_Sloan':
		flt = 'rp'
	if flt == 'g_Sloan':
		flt = 'gp'
	if flt == 'i_Sloan':
		flt = 'ip'
	if flt == 'z_Sloan':
		flt = 'zp'

	return flt




def get_filter_weihai(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ","")
	if flt == 'r1':
		flt='rp'
	if flt == 'i1':
		flt = 'ip'

	return flt


def get_filter_DEMONEXT(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ","")
	if flt == 'r':
		flt='rp'
	if flt == 'i':
		flt = 'ip'

	return flt


def get_filter_REM(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ","")
	if flt == 'g':
		flt='gp'
	if flt == 'r':
		flt = 'rp'
	if flt == 'i':
		flt = 'ip'

	return flt


def get_filter_GS(img,fltkey):
	fitsinfo = __get_fits_info(img, fltkey)
	flt = fitsinfo[fltkey]
	flt  = flt.replace(" ","")

	if flt == 'Blue':
		flt = 'B'
	if flt == 'Red':
		flt = 'R'

	return flt



def get_filter_LTH2RG(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]


	return flt

def get_filter_LT(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

#	print flt

	flt = flt.replace(" ","")

#	print flt

	if flt == 'Bessell-B':
		flt = 'B'
	elif flt == 'Bessell-V':
		flt = 'V'
	elif flt == 'SDSS-R':
		flt = 'rp'
	elif flt == 'SDSS-I':
		flt = 'ip'
	elif flt == 'SDSS-Z':
		flt = 'zp'
	elif flt == 'SDSS-U':
		flt = 'up'
	elif flt == 'SDSS-G':
		flt = 'gp'
	elif flt == 'MK-J':
		flt = 'J'
	elif flt == 'MK-H':
		flt = 'H'
	elif flt == 'MK-K':
		flt = 'K'
	else:
		print img
		raise KeyError("the filter information in fits header for LT telescope is not understood, please check...")

	return flt


def get_filter_TJO(img,fltkey):

	fitsinfo = __get_fits_info(img,fltkey)
	flt = fitsinfo[fltkey]

	flt = flt.replace(" ","")

	return flt


def get_filter_SMART(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
        flt = fitsinfo[fltkey]

        flt = flt.replace(" ","")

        return flt


def get_filter_SMART_IR(img,fltkey):
	fitsinfo = __get_fits_info(img,fltkey)
        flt = fitsinfo[fltkey]

        flt = flt.replace(" ","")

        return flt

def get_bitpix(img,bitpix_key,telescope):
	if telescope == 'WFCAM':
		fitsinfo = __get_fits_info(img, bitpix_key, extension =1)
		bitpix = fitsinfo[bitpix_key]
	else:
		fitsinfo = __get_fits_info(img,bitpix_key)
       	 	bitpix = fitsinfo[bitpix_key]
        return bitpix

def get_airmass_telescope(img, airmasskey, telescope):

	if telescope == "WFCAM":
		airmass = get_airmass_WFCAM(img, airmasskey)
	elif telescope == "UFTI":
		airmass = get_airmass_WFCAM(img, airmasskey)
	elif telescope == "LCOGT":
		airmass = get_airmass_direct(img, airmasskey)
	elif telescope == 'REM':
		airmass = get_airmass_direct(img, airmasskey)
	elif telescope == 'REM-NIR':
		airmass = get_airmass_direct(img, airmasskey)
	elif telescope == 'SMARTS':
		airmass = get_airmass_direct(img, airmasskey)
	elif telescope == 'SMARTS-IR':
		airmass = get_airmass_direct(img, airmasskey)
	else:
		raise KeyError("Telescope %s is not supported yet"%telescope)

	return airmass

def get_airmass_direct(img, amkey):

	fitsinfo = __get_fits_info(img, amkey)
	airmass = fitsinfo[amkey]

	return airmass

def get_airmass_WFCAM(img, amkeys):

	amkey1 = amkeys.split()[0]
	amkey2 = amkeys.split()[1]

	fitsinfo = __get_fits_info(img, amkeys.split())

	amstart = fitsinfo[amkey1]
	amend = fitsinfo[amkey2]

	airmass = (amstart + amend)/2.0

	return airmass
