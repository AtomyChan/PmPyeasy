#! /usr/bin/python

#some commonly used functions

from astropy.io import fits
from numpy.linalg import inv
import sys
import os
import numpy as np
from scipy.optimize import leastsq


def get_weighted_average(values, errs):
	'''
	get weighted average from multiple measurements with uncertainties as weights
	Refer to the following pages	

	http://www.colorado.edu/physics/phys2150/phys2150_sp14/phys2150_lec4.pdf
	http://labrad.fisica.edu.uy/docs/promedios_ponderados_taylor.pdf


	Inputs:
	values:
	errs: measurements uncertainty 

	Outputs:
	wtave: weighted average
	wtunc: uncertainty on averaged result
	'''
	if isinstance(values,list):
		values = np.array(values)	
	
	if isinstance(errs, list):
		errs = np.array(errs)

	weights = 1./errs**2
        wtave = np.sum(weights*values) / np.sum(weights)
        wtunc = 1./ np.sqrt(np.sum(weights))

	return wtave, wtunc


def get_xy_on_image_from_ds9(input_image, coordinate_system='image'):
        '''
	load the input_image to ds9 and pick the point you want with mouse
	INPUTS:
		input_image:
		coordinate_system: image or fk5(in degree)
        '''

	try:
		import pyds9
	except:
		print "pyds9 not available"

        d = pyds9.DS9()
        d.set('fits %s'%input_image)
        d.set('scale zscale')
        d.set('zoom to fit')

	#currently coordinate_system support 'image' and 'fk5'
        xys = d.get('imexam coordinate %s'%coordinate_system)
        x,y = xys.split()
        xy = np.array([float(x),float(y)])
        print xy

        return xy

def select_rows_from_table(table,colname,value,mode='eq',criteria = 0):
        '''
        select rows in table which meet specific requirements given by mode and criteria                
        available mode 'eq','lt','gt','lgt'

        'eq':X == value; 'lt': X<= value; 'gt': X>=value; 'lgt': value-criteria  <= X <= value+criteria  
        if mode is 'lgt' then non-zero criteria is needed, otherwise 'lgt' is the same as 'eq' with default criteria

        return the table by rows which meet the criteria
        '''
        colvalues = table[colname]
        if mode   == 'eq':
                mask  = colvalues == value
        elif mode == 'lt':
                mask  = colvalues <= value
        elif mode == 'gt':
                mask  = colvalues >= value
        elif mode == 'lgt':
                if criteria == 0:
                        print "Attention: if mode is 'lgt' then non-zero criteria is needed, otherwise 'lgt' is the same as 'eq' with default criteria!!!"
                mask = np.logical_and(colvalues>=value-criteria,colvalues<=value+criteria)
        else:
                raise KeyError("the criteria %s is not available"%criteria)

        selected = table[mask]

        return selected


def radec_format_transformation_list(ras, decs, mode='str2float', delimiter=':'):
	'''
	See radec_format_transformation
	'''
	ras_new = []
	decs_new = []
	for ra, dec in zip(ras, decs):
		ra_new, dec_new = radec_format_transformation(ra, dec, mode=mode, delimiter=delimiter)
		ras_new.append(ra_new)
		decs_new.append(dec_new)

	return np.array(ras_new), np.array(decs_new)

def radec_format_transformation(ra,dec,mode='str2float',delimiter=':'):
	'''
	xx:xx:xx.xx to xx.xxxxx
	xxhxxmxx.xxs to xx.xxxx
	
	mode:	'str2float'; 'float2str'
	delimiter: ':'; 'hdms'
	'''
	from astropy import coordinates as coor
	from astropy import units as u
	from astropy.coordinates import SkyCoord	

	if mode == 'str2float':

		if delimiter == ':':

			if isinstance(ra,str):
				ra_str = ra
			else:
				raise IOError("ra must be string when mode is str2float")

			if isinstance(dec,str):
				dec_str = dec
			else:
				raise IOError("dec must be string when mode is str2float")			

			ra_segs = ra_str.split(':')
			rahms = ra_segs[0]+'h'+ra_segs[1]+'m'+ra_segs[2]+'s'

			dec_segs = dec_str.split(':')
			decdms = dec_segs[0]+'d'+dec_segs[1]+'m'+dec_segs[2]+'s'

		elif delimiter == 'hdms':
			rahms = ra
			decdms = dec
		else:
			raise IOError('Invalid input for delimiter')
		

		sn_coor = rahms+' '+decdms
		sn_coor_degree =  SkyCoord(rahms, decdms, frame='fk5', unit=(u.hourangle, u.deg))
		ra_nf = sn_coor_degree.ra.value
		dec_nf = sn_coor_degree.dec.value
		
	elif mode == 'float2str':
		c  = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='fk5')
		radec = c.to_string('hmsdms')
		ra_nf  =  radec.split()[0]
		dec_nf =  radec.split()[1]

		if delimiter == ':':
			ra_nf = ra_nf.replace('h',':')
			ra_nf = ra_nf.replace('m',':')
			ra_nf = ra_nf.replace('s','')
			dec_nf = dec_nf.replace('d',':')
			dec_nf = dec_nf.replace('m',':')
			dec_nf = dec_nf.replace('s','')
	else:
		raise IOError("Invalid input for mode")
		
	return ra_nf,dec_nf


def save_results(dump_filename, data):	
	import gzip
	import cPickle as pickle
	pickle.dump(data, gzip.open(dump_filename, "wb", compresslevel=3), protocol=2)

def restore_results(dump_filename):	
	import gzip
	import cPickle as pickle
	return pickle.load(gzip.open(dump_filename, "rb"))

def remove_row_in_array(array,clm,rm_value,mode='eq'):
	'''
	Remove the specific row in array which match the 'rm_value' in 'clm'
	input:
		array: the array of interest
		clm: the column which want to be checked, 0 for the first column 
		rm_value: the critical value to conduct the check
		mode: acceptable inputs are 'eq'(equal),'lt'(less than) and 'gt'(greater than)
	'''
	from pylab import find

	N = len(array)
	indef_mask= range(N)

	if array.shape == (N,):
		try:
			mat_column = array[clm]
		except Exception as e:
			print e
	else:
		mat_column = array[:,clm]

	if mode =='eq':
		remove = find(mat_column == rm_value)
	elif mode == 'lt':
		remove = find(mat_column < rm_value)
	elif mode == 'gt':
		remove = find(mat_column > rm_value)
	else:
		raise TypeError,'Acceptable inputs for mode eq lt and gt'
#	print remove
	for rm in remove:
		indef_mask.remove(rm)
	indef_mask = np.array(indef_mask)
	#print indef_mask
	
	if array.shape == (N,):
		ret = array[indef_mask]
	else:
		ret = array[indef_mask,:]
	
	return ret,indef_mask


def find_duplicates(a, key):
    """
    Find duplicates in a column of a recarray. This is a simplified version of:
    ::

        import numpy.lib.recfunctions as rfn
        rfn.find_duplicates(...)
    """
    a = np.asanyarray(a).ravel()
    # Get the sorting data (by selecting the corresponding field)
    base = a[key]
    # Get the sorting indices and the sorted data
    sortidx = base.argsort()
    sorteddata = base[sortidx]
    # Compare the sorting data
    flag = (sorteddata[:-1] == sorteddata[1:])
    flag = np.concatenate(([False], flag))
    # We need to take the point on the left as well (else we're missing it)
    flag[:-1] = flag[:-1] + flag[1:]
    duplicates = a[sortidx][flag]
    duplicates_index = sortidx[flag]
    return (duplicates, duplicates_index)


def removeoutliers(data, nsigma, remove='both', center='mean', niter=np.Inf, retind=False, verbose=False):
	
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

	from numpy import median, ones, isfinite
	from pylab import find

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

def funcfit_remove_outlier(x,y,func,p0,err=None, nsig=3, rm_mode = 'worst',diff_mode = 'useerr'):
	"""
	Identify outliers considering the mean (if meanfunc=np.mean) or median (if meanfunc=np.median) value
	and 3 sigma (3*stdev),	iterating until convergence.
	INPUT:
		data:
		nsig:
		meanfunc:
		
	"""
	data_model = func(x,p0)
	if err is None:
		err = np.ones(x.shape)
		
	def residuals(p, y, x,yerr):			
		return (y - func(x, p))/yerr

	plsq = leastsq(residuals, p0, args=(y, x,err))
	
	last_total = len(y)
	if diff_mode != 'useerr':
		diff_err = np.ones(x.shape)
	else:
		diff_err = err
		
	stdev = np.std(residuals(plsq[0],y,x,diff_err))
	absdiff = np.abs(residuals(plsq[0],y,x,diff_err))
	if rm_mode == 'all':
		goodflt = absdiff < nsig*stdev
		x = x[goodflt]
		y = y[goodflt]
		err = err[goodflt]
		diff_err = diff_err[goodflt]
	elif rm_mode == 'worst':
		indice = np.argsort(absdiff)
		if absdiff[indice[-1]] > nsig*stdev:
			goodflt = indice[:-1]
			x = x[goodflt]
			y = y[goodflt]
			err = err[goodflt]
			diff_err = diff_err[goodflt]
		
			
	current_total = len(y)
	# Continue iterating until convergence (no more points are removed)
	while last_total > current_total:		
		#print current_total, stdev
		last_total = current_total
		plsq = leastsq(residuals, plsq[0], args=(y, x,err))
		stdev = np.std(residuals(plsq[0],y,x,diff_err))
		absdiff = np.abs(residuals(plsq[0],y,x,diff_err))
		if rm_mode == 'all':
			goodflt = absdiff < nsig*stdev
			x = x[goodflt]
			y = y[goodflt]
			err = err[goodflt]
			diff_err = diff_err[goodflt]
		elif rm_mode == 'worst':
			indice = np.argsort(absdiff)
			if absdiff[indice[-1]] > nsig*stdev:
				goodflt = indice[:-1]
				x = x[goodflt]
				y = y[goodflt]
				err = err[goodflt]
				diff_err = diff_err[goodflt]
		current_total = len(y)
	
	gooddata = np.hstack((x.reshape((len(x),1)),y.reshape((len(y),1)),err.reshape((len(err),1))))
	print gooddata
	pret= plsq[0]
	return gooddata,pret

def sigma_clipping(data, sig=3, meanfunc=np.median, display_result = False):
	"""
	Identify outliers considering the mean (if meanfunc=np.mean) or median (if meanfunc=np.median) value and 3 sigma (3*stdev),
	iterating until convergence.
	"""
	last_total = len(data)
	# First iteration
	stdev = np.std(data)
	diff = data - meanfunc(data)
	sfilter = np.abs(diff) < sig*stdev
	current_total = len(data[sfilter])
	# Continue iterating until convergence (no more points are removed)
	while last_total > current_total:		
		#print current_total, stdev
		last_total = current_total
		stdev = np.std(data[sfilter])
		diff = data - meanfunc(data[sfilter])
		sfilter = np.abs(diff) < sig*stdev
		current_total = len(data[sfilter])

	if display_result:
		import matplotlib.pylab as plt

		plt.plot(data, 'ko', label='Original')
		plt.plot(data[sfilter], 'ro', label='after clipping')
		plt.legend()
		plt.show()
	

	return data[sfilter],sfilter,stdev


def check_and_handle_file(filename):
	'''
	if file 'filename' exist, remove the existing file
	'''
	if os.path.isfile(filename):
		os.remove(filename)

def mkdir_p(path):
	"""
	Creates a directory. Same behaviour as 'mkdir -p'.
	"""
	import errno

	try:		
		os.makedirs(path)
	except OSError as exc: # Python >2.5		
		if exc.errno == errno.EEXIST:			
			pass
		else:			
			raise
def mkdir_rm(path):
	'''
	Create a directory. If the path exists, delete the directory and all files in it.
	'''
	import errno

	try:		
		os.makedirs(path)
	except OSError as exc: # Python >2.5		
		if exc.errno == errno.EEXIST:			
			delete_file_folder(path)
			os.makedirs(path)
		else:			
			raise
	

def delete_file_folder(src):
    '''
    delete files and folders
    '''
    if os.path.isfile(src):  
        try:  
            os.remove(src)  
        except:  
            pass
    elif os.path.isdir(src):  
        for item in os.listdir(src):  
            itemsrc=os.path.join(src,item)  
            delete_file_folder(itemsrc)  
        try:  
            os.rmdir(src)  
        except:  
            pass 



def get_fits_info(data_dir,files,keywords,verbose = False):
	'''
	Get the wanted keywords value in fits header
	
	INPUT:
		data_dir: the directory where the fits files are
		files: fits file names
		keywords
	'''
	info = {}
	message = ''
	
	for key in keywords:
		info[key] = []
		message += key + ' '
		
	print "Get fits info " + message
	N = len(files)
#	print files
	for i,filename in enumerate(files):
		image = os.path.join(data_dir, filename)
		if verbose:
			print image

		info_this = __get_fits_info(image,keywords)

		for key in keywords:
			info[key].append(info_this[key])

		if not verbose:
 		    progress_report(i,N,'get fits header info')
		
	return info

def __get_fits_info(fits_image,keywords, extension=None):

	info_this = {}
	
       	try:
	    hdu = fits.open(fits_image)
	except Exception,e:
	    print "IOError on ",fits_image
	    raise IOError,e

	if extension is None:
		header = hdu['PRIMARY'].header
	else:
		header = hdu[extension].header

	if isinstance(keywords,str):
		if keywords in header.keys():
			value = header[keywords]
		else:
			value = None
		info_this[keywords] = value
	elif isinstance(keywords,list) or isinstance(keywords,numpy.ndarray):
		for key in keywords:
			if key in header.keys():
				value = header[key]
			else:
				value = None
			info_this[key] = value
#		print value
	hdu.close()
	
	return info_this


def report_progress(current_work_progress, last_reported_progress):
	"""
	returns: True every 10% of progress.
	"""
	return (int(current_work_progress) % 10 == 0 and current_work_progress - last_reported_progress > 10) or last_reported_progress < 0 or current_work_progress == 100

def progress_report(current,total,message,bar_length=100):
	'''
	report the working progress
	'''
	percent = float(current) / total
	progress = '=' * (int(round(percent * bar_length))-1)+'>'
	spaces = ' ' * (bar_length - len(progress))
	if current == (total-1):
		sys.stdout.write("\r[{0}]: [{1}] {2}%\n".format(message,progress + spaces, 100))
	else:
		sys.stdout.write("\r[{0}]: [{1}] {2}%".format(message,progress + spaces, int(round(percent * 100))))
	sys.stdout.flush()
	
def optimize_over_two(var1,var2):
	'''
	When we want to select the images with the low background and small fwhm as the same time, we need to balance the fight between bkg and fwhm
	'''
	score = (var1-min(var1))/(max(var1)-min(var1))+(var2-min(var2))/(max(var2)-min(var2))
	indice = np.argsort(score)
	return indice

def get_best(inputpara,factors,num):
	'''
	np.max()-np.min()
	'''
	paras = inputpara.copy()

	N,M = paras.shape
	m   = len(factors)
	if M != m:
		raise OSError('The number of given criterias and the number of parameters to be considered are not consistent!')
	
	indice = np.arange(N)
	#print indice
	#for i in range(m):
	#	para = paras[indice,i]
	#	pararemain,sfilter,std = sigma_clipping(para,sig=2)
	#	print sfilter
	#	#print pararemain
	#	#print sfilter
	#	#print type(sfilter)
	#	indice = indice[sfilter]
	#	print indice

	#paras = paras[indice,:]
	grades = np.zeros((len(indice),))
	for j in range(m):
		para = paras[:,j]
		#print para
		#print factors[j]
		grade_this= factors[j]*(para-np.min(para))/(np.max(para)-np.min(para))
		#print grade_this
		grades = grades + grade_this
#		plt.plot(para)

	#print grades
#	plt.plot(grades)
#	plt.show()

	sortindice = np.argsort(grades)	
	small_indice = sortindice[range(num)]
	
	wantedindice  = indice[small_indice]
	ret = inputpara[wantedindice]
	print ret

	return ret,wantedindice
	

def get_matrix(xa,ya,xb,yb):
    '''
    (xa,ya) is the independent variables
    (xb,yb) is the dependent variables
    '''
    #check the input 
    if len(xa)!=len(ya) or len(xb)!=len(yb):
        raise OSError('measurements for independent variable or dependent variable have different length')
    if len(xa)==len(xb):
        N = len(xa)
    else:
        raise OSError('the length of dependent and independent variable differs')
    
    XA = xa.sum()
    YA = ya.sum()
    XB = xb.sum()
    YB = yb.sum()
    XAXA = (xa*xa).sum()
    XAYA = (xa*ya).sum()
    XAXB = (xa*xb).sum()
    XAYB = (xa*yb).sum()
    YAYA = (ya*ya).sum()
    YAYB = (ya*yb).sum()
    XBYA = (xb*ya).sum()
    #get the first matrix and vector
    C1 = np.matrix([[N,XA,YA],[XA,XAXA,XAYA],[YA,XAYA,YAYA]])
    D1 = np.matrix([[XB],[XAXB],[XBYA]])
    
    #get the second matrix and vector
    C2 = np.matrix([[N,XA,YA],[XA,XAXA,XAYA],[YA,XAYA,YAYA]])
    D2 = np.matrix([[YB],[XAYB],[YAYB]])
    
    return C1,D1,C2,D2
    
def least_square_solution(xa,ya,xb,yb):
    '''
    this is a python relization of least square of two variables
    measurements have the form (X_Ai,Y_Ai) (X_Bi,Y_Bi)
    the given model is  X_Bi = dx + a*X_Ai + b*Y_Ai
                        Y_Bi = dy + c*X_Ai + d*Y_Ai
    (xa,ya) is the independent variables
    (xb,yb) is the dependent variables
    '''
    C1,D1,C2,D2 = get_matrix(xa,ya,xb,yb)
    para1 = inv(C1)*D1
    para2 = inv(C2)*D2
    
    #X_Bi = para1[0] + para1[1]*X_Ai + para1[2]*Y_Ai
    #Y_Bi = para2[0] + para2[1]*X_Ai + para2[2]*Y_Ai
    return para1,para2

def show_histogram(x, xlabel='Units', nbins=50):
	"""
	Build a histogram from 'x' values and plot it.
	"""
	import matplotlib.pyplot as plt
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# the histogram of the data
	n, bins, patches = ax.hist(x, nbins, normed=0, facecolor='green', alpha=0.75)
	# Mean
	l = plt.axvline(x = np.mean(x), linewidth=1, color='red')
	ax.annotate('Mean', xy=(np.mean(x), np.max(n)),  xycoords='data',xytext=(10, 20), textcoords='offset points', \
	    size=8,bbox=dict(boxstyle="round", fc="0.8"),arrowprops=dict(arrowstyle="->", \
	    connectionstyle="angle,angleA=0,angleB=90,rad=10",edgecolor='black'), \
	    horizontalalignment='right', verticalalignment='top',)
	# Median
	l = plt.axvline(x = np.median(x), linewidth=1, color='orange')
	ax.annotate('Median', xy=(np.median(x), np.max(n)),  xycoords='data',xytext=(10, 35), textcoords='offset points', \
	    size=8,bbox=dict(boxstyle="round", fc="0.8"),arrowprops=dict(arrowstyle="->", \
	    connectionstyle="angle,angleA=0,angleB=90,rad=10",edgecolor='black'), \
	    horizontalalignment='right', verticalalignment='top',)
	ax.set_xlabel(xlabel)
	ax.set_ylabel('Counts')
	ax.grid(True)
	plt.show()
