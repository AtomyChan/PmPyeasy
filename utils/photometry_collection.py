import numpy as np
from pyraf import iraf
from iraf import digiphot,apphot
import os
from astropy.io import fits

try:
	from aperture_photometry import aperphot
except:
	print "self customised aperture photometry function not loaded"

def apphot_selfcustomized(image,xy_pnts,coo_list,options):


    hdu = fits.open(image)
    data = hdu[0].data
    hdr  = hdu[0].header
    flt = hdr['FILTER']

    f = open(coo_list, 'w')
    if np.array(xy_pnts).shape == (2,): #find if there is only one set of coords to measure
        f.write(str(xy_pnts[0]) + ' ' + str(xy_pnts[1]) + '\n')
        N = 1
    else:
        for this_star in range(len(xy_pnts)):
            f.write(str(xy_pnts[this_star][0]) + ' ' + str(xy_pnts[this_star][1]) + '\n')
        N = len(xy_pnts)
    f.flush()
    f.close

    dap_tg = options['app']*2
    dap_insky = options['skyin']*2
    dap_outsky = options['skywidth']*2

    if N == 1:
        pos = xy_pnts
        thisstar = aperphot(data, timekey=None, pos=pos, dap=[dap_tg,dap_insky,dap_outsky])
        phot_res = [pos[0],pos[1],thisstar.phot,thisstar.bg,thisstar.ebg,thisstar.nsky,thisstar.ntarg]
    else:
        phot_res  = np.zeros((N,7))
        for i,star in enumerate(xy_pnts):
            pos = star[0:2]
            print i
            thisstar = aperphot(data, timekey=None, pos=pos, dap=[dap_tg,dap_insky,dap_outsky])
            phot_res[i] = [pos[0],pos[1],thisstar.phot,thisstar.bg,thisstar.ebg,thisstar.nsky,thisstar.ntarg]

    return phot_res

def apphot_iraf(image, xy_pnts, calgorithm, coo_list, ap_list,options, datamin = None, datamax = None, epadu=None, readnoise=None):
    '''
	INPUTS:
		image:
		xy_pnts:
		calgorithm: 'centroid', 'gauss', 'none', 'ofilter'
		coo_list:
		ap_list:
		options:
		datamin:
		datamax:
    '''
    f = open(coo_list, 'w')
    if np.array(xy_pnts).shape == (2,): #find if there is only one set of coords to measure
        f.write(str(xy_pnts[0]) + ' ' + str(xy_pnts[1]) + '\n')
    else:
        for this_star in range(len(xy_pnts)):
            f.write(str(xy_pnts[this_star][0]) + ' ' + str(xy_pnts[this_star][1]) + '\n')
    f.flush()
    f.close

    apphot.phot.unlearn()
    apphot.datapars.unlearn()
    apphot.datapars.scale = 1.0
    apphot.datapars.fwhmpsf = options['fwhmpsf']
    apphot.datapars.emission = True
    apphot.datapars.sigma = 'INDEF'

    if datamin is None:
        datamin = 'INDEF'
    apphot.datapars.datamin = datamin

    if datamax is None:
        datamax = 46000
    apphot.datapars.datamax = datamax

    apphot.datapars.noise = 'poisson'

    apphot.datapars.ccdread = ''
    if readnoise is not None:
	apphot.datapars.readnoise = readnoise
    else:
        apphot.datapars.readnoise = 13.5
    apphot.datapars.itime = 1.0
    if epadu is not None:
	apphot.datapars.epadu = epadu
    else:
        apphot.datapars.epadu = 1.4
    apphot.datapars.xairmass = 'INDEF'
    apphot.datapars.ifilter = 'INDEF'
    apphot.datapars.otime = 'INDEF'

    # iraf.digiphot.apphot.centerpars :
    if calgorithm not in ['centroid', 'gauss', 'none', 'ofilter']:
	raise ValueError("center algorithm %s not supported"%calgorithm)

    apphot.centerpars.calgorithm = calgorithm

    apphot.centerpars.cbox = 10.0
    apphot.centerpars.cthreshold = 0.0
    apphot.centerpars.minsnratio = 1.0
    apphot.centerpars.cmaxiter = 10.0
    apphot.centerpars.maxshift = 2.0
    apphot.centerpars.clean = False
    apphot.centerpars.rclean = 1.0
    apphot.centerpars.rclip = 2.0
    apphot.centerpars.kclean = 3.0
    apphot.centerpars.mkcenter = False

    # iraf.digiphot.apphot.fitskypars :
    apphot.fitskypars.unlearn()
    apphot.fitskypars.salgorithm = 'median'
    apphot.fitskypars.annulus = options['skyin']
    apphot.fitskypars.dannulus = options['skywidth']
    apphot.fitskypars.skyvalue = 0.0
    apphot.fitskypars.smaxiter = 10.0
    apphot.fitskypars.sloclip = 0.0
    apphot.fitskypars.shiclip = 0.0
    apphot.fitskypars.snreject = 50.0
    apphot.fitskypars.sloreject = options['sky_sigma_down']
    apphot.fitskypars.shireject = options['sky_sigma_up']
    apphot.fitskypars.khist = 3.0
    apphot.fitskypars.binsize = 0.1
    apphot.fitskypars.smooth = False
    apphot.fitskypars.rgrow = 0.0
    apphot.fitskypars.mksky = False

    # iraf.digiphot.apphot.photpars :
    apphot.photpars.unlearn()
    apphot.photpars.weighting = 'constant'
    apphot.photpars.apertures = options['app']
    apphot.photpars.zmag = options['def_zeropt']
    apphot.photpars.mkapert = False

    photparams = {
        'radplot':False,
        }

    if os.path.exists(ap_list):
        os.remove(ap_list)

    # run photometry using the newly created coxyfile for providing input coordinates
    apphot.phot(image=image, skyfile='', output=ap_list,coords=coo_list, verify='no',interactive='no',verbose=True, Stdout=1, **photparams)

    f = open(ap_list, 'r')
    maglines = f.readlines()
    f.close()

    #get the photometric data from the apphot output
    xpos = np.array([])
    ypos = np.array([])
    cent_err = np.array([])

    bkg = np.array([])
    bkg_stdev = np.array([])
    sky_counts = np.array([])

    target_area = np.array([])
    mag = np.array([])
    flux = np.array([])
    mag_err = np.array([])
    phot_err = np.array([])


    if np.array(xy_pnts).shape == (2,): #find if there is only one set of coords to measure
        number = 1
    else:
        number = len(xy_pnts)
    for i in range(number):

        xpos = np.append(xpos, float(maglines[75 + i * (4 + 1) + 1].split()[0]))
        ypos = np.append(ypos, float(maglines[75 + i * (4 + 1) + 1].split()[1]))
        cent_err = np.append(cent_err, maglines[75 + i * (4 + 1) + 1].split()[7])

        if maglines[75 + i * (4 + 1) + 2].split()[0]=='INDEF':
            #mag=float('nan') ## this would be preferable, but formatting below makes this impossible
            bkg = np.append(bkg, float('99.999'))
        else:
            bkg = np.append(bkg, float(maglines[75 + i * (4 + 1) + 2].split()[0]))

        if maglines[75 + i * (4 + 1) + 2].split()[1]=='INDEF':
            bkg_stdev = np.append(bkg_stdev, float('99.999'))
        else:
            bkg_stdev = np.append(bkg_stdev, float(maglines[75 + i * (4 + 1) + 2].split()[1]))

        sky_counts = np.append(sky_counts, float(maglines[75 + i * (4 + 1) + 4].split()[3]))
        target_area = np.append(target_area, float(maglines[75 + i * (4 + 1) + 4].split()[2]))


        if maglines[75 + i * (4 + 1) + 4].split()[4]=='INDEF':
            #mag=float('nan') ## this would be preferable, but formatting below makes this impossible
            mag = np.append(mag, float('99.999'))
        else:
            mag = np.append(mag, float(maglines[75 + i * (4 + 1) + 4].split()[4]))

        if maglines[75 + i * (4 + 1) + 4].split()[3]=='INDEF':
            #mag=float('nan') ## this would be preferable, but formatting below makes this impossible
            flux = np.append(flux, float('99.999'))
        else:
            flux = np.append(flux, float(maglines[75 + i * (4 + 1) + 4].split()[3]))


        if maglines[75 + i * (4 + 1) + 4].split()[5]=='INDEF':
            #mag_err=float('nan') ## this would be preferable, but formatting below requires something else
            mag_err = np.append(mag_err, float('99.999'))
        else:
            mag_err = np.append(mag_err, float(maglines[75 + i * (4 + 1) + 4].split()[5]))

        if len(maglines[75 + i * (4 + 1) + 4].split()) < 8:
            phot_err = np.append(phot_err, float('99.999'))
        else:
            if maglines[75 + i * (4 + 1) + 4].split()[7]=='INDEF':
                phot_err = np.append(phot_err, float('99.999'))
            else:
                phot_err = np.append(phot_err, maglines[75 + i * (4 + 1) + 4].split()[7])


    ret = np.hstack((xpos.reshape((len(xpos),1)), ypos.reshape((len(xpos),1)), flux.reshape((len(xpos),1)), bkg.reshape((len(xpos),1)), \
                 bkg_stdev.reshape((len(xpos),1)), sky_counts.reshape((len(xpos),1)), target_area.reshape((len(xpos),1)), \
                mag.reshape((len(xpos),1)),mag_err.reshape((len(xpos),1))))


    return ret
