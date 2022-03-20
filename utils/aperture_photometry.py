from numpy import meshgrid, median,isfinite,sort,ndarray,string_,sign, pi, nan, arange
from numpy import std, sum, mean,tile, concatenate, floor, Inf,sqrt, exp, nan, max,Inf
from numpy import array, isfinite, swapaxes,zeros, bool,histogram,sqrt, linspace, ones,pi

import numpy as np
from astropy.io import fits
from os import path
from scipy import interpolate
from pylab import find
from scipy import optimize


class phot:
    def __init__(self, **kw):
        '''
        Generate a photometry object w/keywords and kw+'str' descriptive keywords
        Inputs cannot be numpy arrays
        '''
        keylist = ['time', 'phot', 'ephot', 'bg', 'ebg', 'aper', 'position', 'eposition', \
            'ditherposition', 'object', 'filename', 'exptime', 'ra', 'dec']
        defaults = dict(photunit=None)
        for key in keylist:
            defaults[key]= None
            defaults[key+'str']= None
            keylist = keylist + [key+'str']
        for keyword in defaults:
            if (not kw.has_key(keyword)):
                kw[keyword] = defaults[keyword]

        for k in kw.keys():
            if kw[k].__class__==str:
                exec('self.'+k+'="'+str(kw[k]))+'"'
            else:
                exec('self.'+k+'='+str(kw[k]))        
        return 

    def __str__(self): # Return info
        lin= 'Phot: {} at pixel ({:.1f},{:.1f})\n'.format(self.object,self.position[0],self.position[1])
        lin+= '  Analyzed frame {}\n'.format(self.filename)
        lin+= '  Exposure time = {:.1f}s\n'.format(self.exptime)
        lin+= '  Aperture = {}, {}, {}\n'.format(self.aper[0],self.aper[1],self.aper[2])        
        lin+= '  Counts/s = {} +/-  {}\n'.format(self.phot,self.ephot)
        return lin


def aperphot(fn, timekey=None, pos=[0,0], dap=[2,4,6], mask=None, verbose=False, nanval=999, resamp=None, retfull=False):
    '''
    Do aperture photometry on a specified file.

    :INPUTS:
      pos : 2-sequence
        center of apertures (as if indexing: fn[pos[0], pos[1]])
  
      dap : 3-sequence
        Photometry aperture DIAMETERS:
           -- target aperture (within which to sum flux)
           -- inner sky aperture
           -- outer sky aperture
  
      resamp : int
        Factor by which to interpolate frame before measuring photometry
        (in essence, does partial-pixel aperture photometry)
  
      Aperture masking:
         If no mask is passed in, use the star's input position and
         aperture diameters to create binary pixel masks for stellar and
         background photometry.  If a mask is passed in, stellar and
         background aperture masks are generated from where all input
         mask elements equal 1 and 2, respectively.
  
      retfull: 
          Also return arrays of target mask, sky mask, and frame used.
          This option is a memory hog!
  
    :OUTPUTS:  
      :class:`phot` object.  
  
    '''

    thisobs = phot()
    x0, y0 = pos
    dap_targ, dap_skyinner, dap_skyouter = dap
    if resamp is None or resamp<1:
        resamp = 1
    else:
        resamp = float(resamp)
        
    # Determine size:
    if isinstance(fn,str):
        nx = fits.getval(fn, 'NAXIS1')
        ny = fits.getval(fn, 'NAXIS2')
    elif isinstance(fn,ndarray):
        nx,ny = fn.shape

    nx0, ny0 = nx, ny
    nx = ((nx - 1)*resamp + 1.)  # Avoid resampling at pixel locations
    ny = ((ny - 1)*resamp + 1.)  #   outside the original boundaries.

    # Generate or load masks:
    if mask==None:
        xx,yy = meshgrid(np.arange(nx)/resamp, np.arange(ny)/resamp)
        mask_targ = makemask(xx, yy, (x0, y0, dap_targ))
        mask_s1 = makemask(xx, yy, (x0,y0, dap_skyinner))
        mask_s2 = makemask(xx, yy, (x0,y0, dap_skyouter))
        mask_sky = mask_s2 - mask_s1
    else:
        mask_targ = mask==1
        mask_sky = mask==2
        if resamp>1:
            print "In aperphot, resamp>1 and user-specified mask passed in... beware!"

    # Load data frame:
    thisobs = phot()
    if isinstance(fn,ndarray):
        frame = fn
    elif isinstance(fn, str) or isinstance(fn,string_):
        if not path.isfile(fn):
            print "file %s not found! exiting..." % fn
            return thisobs
        frame = fits.getdata(fn)
        fixval(frame, nanval)

    # Resample data frame
    if resamp>1:
        frame0 = frame.copy()
        xx0 = range(nx0)
        yy0 = range(ny0)
        x1,y1 = np.arange(nx)/resamp, np.arange(ny)/resamp
        rectspline = interpolate.fitpack2.RectBivariateSpline(xx0, yy0, frame0, kx=1, ky=1, s=0)
        frame = rectspline(x1, y1)

    # Measure background and aperture photometry
    thisbg, thisebg = estbg(frame, bins = 20, mask=mask_sky, plotalot=verbose, rout=[3,99])
    thisphot = (mask_targ*(frame - thisbg)).sum() /resamp/resamp
    peak = frame.max()
    peak_targ = (mask_targ * frame).max()
    peak_annulus = (mask_sky * frame).max()

    thisobs.bg=thisbg
    thisobs.ebg=thisebg
    thisobs.bgstr='phot.estbg: SDOM on bg histogram mean & dispersion after outlier rejection'
    thisobs.phot=thisphot
    thisobs.photstr='by-hand background-subtracted aperture photometry'
    thisobs.ntarg = mask_targ.sum()/resamp/resamp
    thisobs.nsky = mask_sky.sum()/resamp/resamp

    thisobs.peak = peak
    thisobs.peak_targ = peak_targ
    thisobs.peak_annulus = peak_annulus
    thisobs.peakstr = 'peak pixel value in frame'
    thisobs.peak_targstr = 'peak pixel value in target aperture'
    thisobs.peak_annulusstr = 'peak pixel value in sky annulus'
    thisobs.position = pos
    thisobs.positionstr = 'user-specified, zero-indexed pixel coordinates.'
    if isinstance(fn, str):
        header = fits.getheader(fn)
        if not timekey==None:
            if timekey in header: 
                thisobs.time=header['timekey']
                thisobs.timestr='heliocentric modified julian date'
        if 'object' in header:  thisobs.object = header['object']
        if 'exptime' in header: thisobs.exptime = header['exptime']
    thisobs.aper = dap
    thisobs.aperstr = 'target, inner, outer aperture diameters, in pixels.'
    thisobs.filename=fn
    thisobs.resamp = resamp

    # Simple stats :: JXP 2014 May 6
    var = thisphot + np.sqrt(thisobs.nsky)*thisobs.bg 
    thisobs.ephot = np.sqrt(var)


    if retfull:
        thisobs.mask_targ = mask_targ
        thisobs.mask_sky  = mask_sky
        thisobs.frame = frame

    return thisobs


def estbg(im, mask=None, bins=None, plotalot=False, rout=(3,200), badval=nan, verbose=False):
    '''
    Estimate the background value of a masked image via histogram fitting.

    INPUTS:
      im -- numpy array.  Input image.

    OPTIONAL INPUTS:
      mask -- numpy array. logical mask, False/0 in regions to ignore
      bins -- sequence.  edges of bins to pass to HIST
      plotalot -- bool.  Plot the histogram and fit.
      rout -- 2-tuple of (nsigma, niter) for analysis.removeoutliers.
              Set to (Inf, 0) to not cut any outliers.
      badval -- value returned when things go wrong.

    OUTPUT:
      b, s_b -- tuple of (background, error on background) from gaussian fit.
                 Note that the error is analagous to the standard deviation on the mean

    COMMENTS:
      The fit parameters appear to be robust across a fairly wide range of bin sizes.
      '''

        
    def gaussianChiSquared(guess, x, y, err):
        return (egaussian(guess, x, y, e=err)**2).sum()

    if mask is None:
        mask = ones(im.shape)
    dat = im.ravel()[find(mask<>0)]
    if plotalot:        
        print "mean(img)=",mean(dat), "std(img)=",std(dat), "nsig*std=",rout[0]*std(dat)                
    dat = removeoutliers(dat, rout[0], remove='both', center='mean', niter=rout[1], verbose=plotalot)
    ndat = len(dat)
    
    if ndat==0:
        print "No data to work with!"
        return (badval, badval)
    if bins==None:
        if plotalot or verbose: print "no bins entered!"
        datmean = dat.mean()
        datstd = stdr(dat, nsigma=3)
        dobin = False
    else:
        dobin = True
        if isinstance(bins,int):
            bins = linspace(dat.min(), dat.max(), bins)
        elif isinstance(bins,list) or isinstance(bins,np.ndarray) and len(instance.shape)==1:
            bins = bins
        else:            
            raise IOError("wrong bins was given")

    if plotalot or verbose: 
        print "dat.mean, dat.std>>" + str((dat.mean(), dat.std()))        
    
    if dobin:
        binwidth = mean(bins[1::]-bins[:-1])
        bincenter = 0.5*(bins[1::]+bins[:-1])
        datIndex = (dat>=bins.min()) * (dat<=bins.max())
        hout = histogram(dat[datIndex], bins) #,new=True)
        gy = hout[0]
        erry = sqrt(gy)
        usableIndex = gy>0

        eff_binwidth = mean(bins[usableIndex][1::]-bins[usableIndex][:-1])
        guess = [gy.sum()*eff_binwidth, std(dat[datIndex]), median(dat[datIndex])]

        if 1.0*usableIndex.sum()/usableIndex.size < 0.5:
            out = guess
        else:
            out = optimize.fmin(gaussianChiSquared, guess,args=(bincenter[usableIndex],gy[usableIndex], erry[usableIndex]),disp=plotalot)

        if plotalot:
            from pylab import figure, errorbar, plot, colorbar, title,show
            print 'guess>>',guess
            print 'fit>>',out
            figure()
            errorbar(bincenter[usableIndex], gy[usableIndex], erry[usableIndex], fmt='ob')
            plot(bincenter, gaussian(out, bincenter),'-r', linewidth=2)
            title('Mean: %f, Std. Dev.: %f' % (out[2], out[1]))
            show()
        ret = out[2], out[1]/sqrt(ndat)
    else:
        ret = datmean, datstd/sqrt(ndat)

    return  ret

def makemask(x,y,params, shape='circle'):
    '''
    Generate a binary (logical) mask with given shape and location.

    INPUTS:
         x      = x-coodinate system (made with meshgrid)
         y      = y-coodinate system (made with meshgrid)
      params:
        shape='circle':
           params(1)  = x-offset
           params(2)  = y-offset
           params(3)  = x-diameter
           params(4)  = OPTIONAL y-diameter
        shape='quad':
           params: list of quadrants to include in the mask.  The
              midpoint for determining quadrants will be
              mid = (xmax+xmin)/2.  Quadrants are:
              0:   x<midX  and  y<midY
              1:   x<midX  and  y>=midY
              2:   x>=midX  and  y<midY
              3:   x>=midX  and  y>=midY

    OPTIONAL INPUTS:
        shape=:  desired mask shape.  Currently only 'circle' is valid.

    OUTPUTS:
      mask   = NxM grided representation of specified mask
                  where NxM are the size of the x,y input meshgrids
    '''
    
    if not x.shape==y.shape:
        print "x,y meshgrid coordinates must be the same shape! Exiting."
        return -1
    if shape=='circle':
        if len(params)<3:
            print "Must give at least 3 parameters to mask."
            return -1
        x0 = params[0]
        y0 = params[1]
        xdia =params[2]
        if len(params)==3:
            ydia = xdia
        else:
            ydia = params[3]
        mask = (  (((x-x0)/(xdia/2.))**2 + ((y-y0)/(ydia/2.))**2) < 1 )
    elif shape=='quad':
        midx = (x.max()+x.min())/2.
        midy = (y.max()+y.min())/2.
        mask = zeros(x.shape, bool)
        for ii in range(len(params)):
            if params[ii]==0:
                mask += (x<midx) * (y<midy)
            elif params[ii]==1:
                mask += (x<midx) * (y>=midy)
            elif params[ii]==2:
                mask += (x>=midx) * (y<midy)
            elif params[ii]==3:
                mask += (x>=midx) * (y>=midy)

    return mask



def fixval(arr, repval, retarr=False):
    '''
    Fix non-finite values in a numpy array, and replace them with 'repval'.
    :INPUT:
       arr -- numpy array, with values to be replaced.
       repval -- value to replace non-finite elements with

    :OPTIONAL INPUT:

       retarr -- if False, changes values in arr directly (more
       efficient).  if True, returns a fixed copy of the input array,
       which is left unchanged.
    '''
       
    if retarr:
        arr2 = arr.ravel().copy()
    else:
        arr2 = arr.ravel()

    finiteIndex = np.isfinite(arr2)
    if not finiteIndex.any():
        badIndex = find((1-finiteIndex))
        arr2[badIndex] = repval

    if retarr:
        return arr2.reshape(arr.shape)
    else:
        return arr2



def removeoutliers(data, nsigma, remove='both', center='mean', niter=Inf, retind=False, verbose=False):
    '''
    Strip outliers from a dataset, iterating until converged.
    :INPUT:
      data -- 1D numpy array.  data from which to remove outliers.
      nsigma -- positive number.  limit defining outliers: number of
                standard deviations from center of data.

    :OPTIONAL INPUTS:               
      remove -- ('min'|'max'|'both') respectively removes outliers
                 below, above, or on both sides of the limits set by
                 nsigma.

      center -- ('mean'|'median'|value) -- set central value, or
                 method to compute it.

      niter -- number of iterations before exit; defaults to Inf,
               which can occasionally result in empty arrays returned

      retind -- (bool) whether to return index of good values as
                second part of a 2-tuple.
    '''

    def getcen(data, method):
        "Get central value of a 1D array (helper function)"
        if method.__class__==str:
            if method=='median':
                cen = median(data)
            else:
                cen = mean(data)
        else:
            cen = method
        return cen

    def getgoodindex(data, nsigma, center, stdev, remove):
        "Get number of outliers (helper function!)"
        if stdev==0:
            distance = data*0.0
        else:
            distance = (data-center)/stdev
        if remove=='min':
            goodind = distance>-nsigma
        elif remove=='max':
            goodind = distance<nsigma
        else:
            goodind = abs(distance)<=nsigma
        return goodind

    data = data.ravel().copy()

    ndat0 = len(data)
    ndat = len(data)
    itern=0
    goodind = ones(data.shape,bool)
    goodind *= isfinite(data)
    while ((ndat0<>ndat) or (itern==0)) and (itern<niter) and (ndat>0) :
        ndat0 = len(data[goodind])
        cen = getcen(data[goodind], center)
        stdev = data[goodind].std()
        thisgoodind = getgoodindex(data[goodind], nsigma, cen, stdev, remove)
        goodind[find(goodind)] = thisgoodind
        if verbose:
            print "cen>>",cen
            print "std>>",stdev
        ndat = len(data[goodind])
        itern +=1
        if verbose:
            print ndat0, ndat
    if retind:
        ret = data[goodind], goodind
    else:
        ret = data[goodind]
    return ret
        
def egaussian(p, x, y, e=None):
    '''
    Compute the deviation between the values y and the gaussian defined by p, x:

    p is a three- or four-component array, list, or tuple
    Returns:   y - p3 - p0/(p1*sqrt(2pi)) * exp(-(x-p2)**2 / (2*p1**2))
    if an error array, e (typ. one-sigma) is entered, the returned value is divided by e.

    SEE ALSO:  :func:`gaussian`
    '''
    
    if e is None:
        e=ones(x.shape)
    fixval(e,y.max()*1e10)

    z = (y - gaussian(p, x))/e
    fixval(z,0)

    return z

def gaussian(p, x):    
    '''
    Compute a gaussian distribution at the points x.

    p is a three- or four-component array, list, or tuple:

    y =  [p3 +] p0/(p1*sqrt(2pi)) * exp(-(x-p2)**2 / (2*p1**2))

    p[0] -- Area of the gaussian
    p[1] -- one-sigma dispersion
    p[2] -- central offset (mean location)
    p[3] -- optional constant, vertical offset

    NOTE: FWHM = 2*sqrt(2*ln(2)) * p1  ~ 2.3548*p1

    SEE ALSO:  :func:`egaussian`
    '''
    
    
    if not isinstance(x, np.ndarray):
        x = array(x, dtype=float, copy=False)

    if len(p)==3:
        p = array(p, copy=True)
        p = concatenate((p, [0]))

    return  p[3] + p[0]/(p[1]*sqrt(2*pi)) * exp(-(x-p[2])**2 / (2*p[1]**2))


def stdr(x, nsigma=3, niter=Inf, finite=True, verbose=False, axis=None):
    '''
    Return the standard deviation of an array after removing outliers.    
    :INPUTS:
      x -- (array) data set to find std of

    :OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean    

    SEE ALSO: :func:`removeoutliers`, 
              :func:`numpy.isfinite`
    '''
    
    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=Inf, verbose=verbose)
        return x.std()
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the action along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):            
            ret[ii] = stdr(x[ii], nsigma=nsigma,niter=niter,finite=finite,verbose=verbose, axis=nextaxis)
        return ret
