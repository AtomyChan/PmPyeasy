#! /usr/bin/python 

import numpy as np
import astroobs.obs as obs
from astropy.time import Time

def get_observatory_info(tel_obs):
	
	#Magellan telescope
	LCO_info = {}
	LCO_info['latitude']  = "-29:00:12.0"
	LCO_info['longitude'] = "-70:42:06.0"
	LCO_info['elev']      = 2282
	LCO_info['tz']        = "Chile/Continental"
	
	#La Silla observatory
	LaSilla_info = {}
	LaSilla_info['latitude'] = "-29:15:24.0"	
	LaSilla_info['longitude']= "-70:43:48.0"	
	LaSilla_info['elev']     = 2347
	LaSilla_info['tz']       = "Chile/Continental"	

	#GTC telescope
	Lapalma_info = {}
	Lapalma_info['latitude']  = "28:45:30.0"
	Lapalma_info['longitude'] = "-17:52:48.0"
	Lapalma_info['elev']      = 2327
	Lapalma_info['tz']        = "Atlantic/Canary"
	
	
	#LBT telescope
	MGO_info  = {}
	MGO_info['latitude']  = "32:42:06.0"
	MGO_info['longitude'] = "-109:53:30.0"
	MGO_info['elev']      = 3181.0
	MGO_info['tz']        = 'US/Mountain'

	#MDM Observatory
	MDM_info  = {}
	MDM_info['latitude'] = "31:57:00.0"
	MDM_info['longitude'] = "-111:37:00.0"
	MDM_info['elev'] = 760
	MDM_info['tz'] = 'US/Mountain'

	#APO observatory
	APO_info = {}
	APO_info['latitude']  = "32:46:49" 
	APO_info['longitude'] = "-105:49:13"
	APO_info['elev']      = 2788
	APO_info['tz']        = 'US/Mountain'

	#DEMONEXT at Winer Observatory in Sonoita, Arizona
	#Iowa Gemini at Winer Observatory in Sonoita, Arizona
	WinerO_info = {}
        WinerO_info['latitude']  = "+31:39:56.08"
        WinerO_info['longitude'] = "-110:36:06.42"
        WinerO_info['elev']      = 1515.7
        WinerO_info['tz']        = 'US/Mountain'

	#FLWO observatory:   6.5m MMT, 1.5m Tillinghast optical spectroscopic telescope with FAST instrument
	FLWO_info = {}
        FLWO_info['latitude']  = "31:40:51.4"
        FLWO_info['longitude'] = "-110:52:39.0"
        FLWO_info['elev']      = 2320
        FLWO_info['tz']        = "US/Mountain"



	#VLT telescope 
	VLT_info = {}
	VLT_info['latitude']  = "-24:37:31.5"
	VLT_info['longitude'] = "-70:24:10.1"
	VLT_info['elev']      = 2648
	VLT_info['tz'] = 'Chile/Continental'
	
	#Palomar Observatory
	Palomar_info  = {}
	Palomar_info['latitude'] = "33:21:21.6"
	Palomar_info['longitude'] = "-116:51:46.8"
	Palomar_info['elev'] = 1706
	Palomar_info['tz'] = 'US/Pacific'


	#Mauna Kea observatory
	MK_info = {}
	MK_info['latitude'] = "19:49:28"
	MK_info['longitude']= "-155:28:24"
	MK_info['elev'] = 4205
	MK_info['tz'] = "US/Aleutian"

	#CTIO observatory
	CTIO_info = {}
	CTIO_info['latitude'] = "-30:09:55.0"
	CTIO_info['longitude'] = "-70:48:54.0"
	CTIO_info['elev'] = 2215.0
	CTIO_info['tz']	  = "Chile/Continental"



	#RP CA 
	RP_CA_info = {}
	RP_CA_info['latitude']  = "37:10:00.00"
	RP_CA_info['longitude'] = "-119:24:47.00"
	RP_CA_info['elev']      = 1405
	RP_CA_info['tz']        = 'US/Pacific'

	#RP NM 
	RP_NM_info = {}
	RP_NM_info['latitude']  = "32:54:46.68"
	RP_NM_info['longitude'] = "-105:31:29.00"
	RP_NM_info['elev']      = 2200
	RP_NM_info['tz']        = 'US/Mountain'

	
	#Xingming Observatory
	XM_info  = {}
	XM_info['latitude']    = "43:28:15.0" 
	XM_info['longitude']   = "87:10:39.6"
	XM_info['elev']        = 2080
	XM_info['tz']          = 'Asia/Urumqi'
	
	
	#Dorothy Hill Observatory   Person in charge David Trappett
	DT_info = {}
	DT_info['latitude']    = "-26:29:24.5" 
	DT_info['longitude']   = "152:35:57.42"
	DT_info['elev']        = 580
	DT_info['tz']          = 'Australia/Brisbane'


	#SAAO where SALT is located
	SAAO_info = {}
	SAAO_info['latitude']  = "-32:22:34"
	SAAO_info['longitude'] = "20:48:38"
	SAAO_info['elev']      = 1798
	SAAO_info['tz'] = 'Africa/Johannesburg'

	#CBA Belgium Observatory
	CBAB_info = {}
	CBAB_info['latitude']  = "50:43:10.0"
	CBAB_info['longitude'] = "05:06:00.1"
	CBAB_info['elev'] = 90
	CBAB_info['tz']  = 'Europe/Brussels'

	#CBA Estremadura Observatory
	


	import astroobs
	obslist = astroobs.obs.ObservatoryList()

	if tel_obs == 'LCO' or tel_obs == 'Magellan' or tel_obs == 'dupont':
		observatory_info = LCO_info
	elif tel_obs == 'LaSilla' or tel_obs == 'NTT' or tel_obs == 'REM':
		observatory_info = LaSilla_info
	elif tel_obs == 'Lapalma' or tel_obs == 'GTC' or tel_obs == 'LT' or tel_obs == 'NOT':
		observatory_info = Lapalma_info
	elif tel_obs == 'MGO' or tel_obs == 'LBT':
		observatory_info = MGO_info
	elif tel_obs == "FAST" or tel_obs == "MMT":
		observatory_info = FLWO_info
	elif tel_obs == 'CTIO' or tel_obs == 'GeminiSouth' or tel_obs == 'GS' or tel_obs == 'SOAR' or tel_obs == 'SMARTS':
		observatory_info = CTIO_info
	elif tel_obs == 'VLT':
		observatory_info = VLT_info
	elif tel_obs == 'SAAO' or tel_obs == 'SALT':
		observatory_info = SAAO_info
	elif tel_obs == 'APO':
		observatory_info = APO_info
	elif tel_obs == 'Palomar' or tel_obs == 'Hale':
		observatory_info = Palomar_info
	elif tel_obs == 'MDM' or tel_obs == 'mdm':
		observatory_info = MDM_info
	elif tel_obs in ['MaunaKea', 'MK', 'GeminiNorth', 'GN', 'Keck', 'Subaru', 'CFHT', 'UKIRT', 'IRTF']:
		observatory_info = MK_info
	elif tel_obs == 'DEMONEXT' or tel_obs == 'demonext' or tel_obs == 'Iowa':
		observatory_info = WinerO_info
	elif tel_obs == 'RP_CA' or tel_obs == 'RPCA':
		observatory_info = RP_CA_info
	elif tel_obs == 'RP_NM' or tel_obs == 'RPNM':
		observatory_info = RP_NM_info
	elif tel_obs == 'Xingming' or tel_obs == 'XM':
		observatory_info = XM_info
	elif tel_obs  == 'DT':
		observatory_info = DT_info
	elif tel_obs  == 'CBAB' or tel_obs == 'TV':
		observatory_info = CBAB_info
	elif tel_obs in obslist.obsids:
		
		obs_info_big = obslist[tel_obs]
		observatory_info = {}
		observatory_info['latitude']  = obs_info_big['lat']
		observatory_info['longitude'] = obs_info_big['long']
		observatory_info['elev']      = obs_info_big['elevation']
		observatory_info['tz']        = obs_info_big['timezone']

	else:
		raise IOError("sorry, %s not supported"%tel_obs)


	return observatory_info



def get_night_airmass_for_given_target_at_given_observatory(observatory_longitude, observatory_latitude, observatory_elevation, observatory_timezone, target_ra, target_dec, observation_date_UT, observatory_name = 'Custom_Observatory', target_name= 'this_target'):
	'''
	OUTPUTS:
		
	'''


	observatory = 	obs.Observatory(observatory_name,long=observatory_longitude,lat=observatory_latitude,elevation=observatory_elevation,timezone=observatory_timezone, ut_date=observation_date_UT)

	#observatory.dates is a vector of Dublin Julian Dates
	#Dublin JD epoch:12h Dec 31, 1899 calculation:JD - 2415020 Introduced by the IAU in 1955
	plan_times = observatory.dates + 2415020 
	t_astrosunset = observatory.sunsetastro + 2415020
	t_astrosunrise = observatory.sunriseastro + 2415020

	obstarget = obs.Target(target_ra, target_dec, target_name)

	obstarget.process(observatory)
	plan_airmass = obstarget.airmass

	#on night sky or not? airmass < 2 exists, Yes; airmass <2 not exists, No
	airmass_darknight = plan_airmass[np.logical_and(plan_times>t_astrosunset,plan_times<t_astrosunrise)]
	JD_times_darknight   = plan_times[np.logical_and(plan_times>t_astrosunset,plan_times<t_astrosunrise)]	


	return JD_times_darknight, airmass_darknight, t_astrosunset, t_astrosunrise


def get_airmass_given_time_target_observatory(obsJD, tel_obs, ra_str, dec_str):

	obsdate = Time(obsJD, format='jd', scale='utc').isot
	obsdate = obsdate.split('T')[0]
	year, month, day = obsdate.split('-')
	obs_plan_JD_date = (int(year),int(month),int(day),23,59,59)

	observatory_info = get_observatory_info(tel_obs)

	obsplan_times, obsplan_airmass, tastrosunset, tastrosunrise = get_night_airmass_for_given_target_at_given_observatory(observatory_info['longitude'], observatory_info['latitude'], observatory_info['elev'], observatory_info['tz'], ra_str, dec_str, obs_plan_JD_date)
	
	#for simplicity, 
	offsets = [0, 1, -1, 2, -2]
	add_ind = 0
	
	tastrosunset_cal = tastrosunset 
	tastrosunrise_cal = tastrosunrise
	while (tastrosunset_cal-obsJD)*(tastrosunrise_cal-obsJD)>0:
		add_ind += 1
		offset = offsets[add_ind]
		tastrosunset_cal = tastrosunset + offset	
		tastrosunrise_cal = tastrosunrise + offset

	obsplan_times = obsplan_times + offsets[add_ind]

	airmass_want = np.interp(obsJD, obsplan_times, obsplan_airmass)
	
	return airmass_want


if __name__ == "__main__":
	obsJD = 2458059.5
	tel_obs = 'CBAB'
	ra_str = '05:07:42.72'
	dec_str = '+24:47:56.58'
	
	ret = get_airmass_given_time_target_observatory(obsJD, tel_obs, ra_str, dec_str)
	print ret
