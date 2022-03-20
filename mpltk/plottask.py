import matplotlib.pylab as plt
from matplotlib import gridspec
from photutils import fit_2dgaussian
import numpy as np

def show_2d_gaussian_fit_result(fitdata, gfit):
	'''
	INPUTS:
		fitdata: M x N array
		gfit: photutils.GaussianConst2D object with fitting result
	'''
	fig = plt.figure(figsize=(7,7))
	ax1 = fig.add_axes((0.1,0.1,0.35,0.35))
	ax1.contour(fitdata, cmap='gray')
	ax1.set_aspect('equal')
	modelcontour = ax1.imshow(fitdata, cmap='hsv')
	NY,NX = fitdata.shape
	X,Y = np.meshgrid(np.arange(NX), np.arange(NY))
	fitret = gfit(X,Y)
	ax1.contour(fitret, cmap='gray')
	fig.colorbar(modelcontour, ax=ax1)
	ax1.invert_yaxis()
	ax1.set_xlabel('X')
	ax1.set_ylabel('Y')


	ax2 = fig.add_axes((0.1,0.55,0.35,0.35))
	residual = ax2.imshow(fitdata-fitret, cmap='gray')
	fig.colorbar(residual, ax =ax2)
	ax2.invert_yaxis()
	ax2.set_xlabel('X')
	ax2.set_ylabel('Y')

	x_fwhm = 2*np.sqrt(2.0*np.log(2))*np.abs(gfit.x_stddev.value)
	y_fwhm = 2*np.sqrt(2.0*np.log(2))*np.abs(gfit.y_stddev.value)
	xc = gfit.x_mean.value
	yc = gfit.y_mean.value
	bkg = gfit.constant.value
	amp = gfit.amplitude.value

	xmax = gfit.x_mean.value
	ymax = gfit.y_mean.value
	ax3 = fig.add_axes((0.6,0.1,0.35,0.35))
	ax3.step(np.arange(NX), fitdata[int(ymax),:], where='mid')
	ax3.plot(np.arange(NX), gfit(np.arange(NX), np.zeros(NX)+ymax), 'k')
	ax3.vlines(xc, bkg-0.1*amp, bkg+1.1*amp, linestyles='dashed')
	ax3.hlines(bkg+amp*0.5, xc-y_fwhm*0.5, xc+y_fwhm*0.5, colors='r', linewidth=3, alpha=0.3)
	ax3.set_title('fwhm$_X$=%s'%str(np.round(y_fwhm,1)))
	ax3.set_xlabel('X')

	ax4 = fig.add_axes((0.6,0.55,0.35,0.35))
	ax4.step(np.arange(NY), fitdata[:,int(xmax)], where='mid')
	ax4.plot(np.arange(NY), gfit(np.zeros(NY)+xmax, np.arange(NY)), 'k')
	ax4.vlines(yc, bkg-0.1*amp, bkg+1.1*amp, linestyles='dashed')
	ax4.hlines(bkg+amp*0.5, yc-x_fwhm*0.5, yc+x_fwhm*0.5, colors='r', linewidth=3, alpha=0.3)
	ax4.set_title('fwhm$_Y$=%s'%str(np.round(x_fwhm,1)))
	ax4.set_xlabel('Y')

	plt.show()



def d1x1d1x2d1yd2y_plot_two_tables(data1, data2, x1col1, x1col2, y1col, y2col, x1errcol1=None, x1errcol2=None, y1errcol=None, y2errcol=None,color='k', alpha=1.0, newaxes=1):
	'''
	plot x1col1 (from data1) - x1col2 (from data1)  VS. y1col (from data1) - y2col (from data2)
	'''
	x1data = data1[x1col1]
	x2data = data1[x1col2]
	y1data = data1[y1col]
	y2data = data2[y2col]
	x1mx2data = x1data - x2data
	y1my2data = y1data - y2data
	if x1errcol1 is not None and x1errcol2 is not None:
		x1mx2errdata = np.sqrt(data1[x1errcol1]**2+data1[x1errcol2]**2)
	else:
		x1mx2errdata = None
	if y1errcol is not None and y2errcol is not None:
		y1my2errdata = np.sqrt(data1[y1errcol]**2+data2[y2errcol]**2)
	else:
		y1my2errdata = None
	xy_simple_scatter_plot(x1mx2data, y1my2data, xerrdata=x1mx2errdata, yerrdata=y1my2errdata, xlabel=x1col1+'-'+x1col2, ylabel=y1col+'-'+y2col, color=color, alpha=alpha, newaxes=newaxes)


def d1xd1yd2y_plot_two_tables(data1, data2, xcol, y1col, y2col, xerrcol=None, y1errcol=None, y2errcol=None):
	'''
	plot xcol from data1 VS. y1col(from data1) - y2col(from data2)
	'''
	xdata = data1[xcol]
	y1data = data1[y1col]
	y2data = data2[y2col]
	y1my2data  = y1data - y2data
	if xerrcol is not None:
		xerrdata = data1[xerrcol]
	else:
		xerrdata = None
	if y1errcol is not None and y2errcol is not None:
		y1my2errdata = np.sqrt(data1[y1errcol]**2+data2[y2errcol]**2)
	else:
		y1my2errdata = None
	xy_simple_scatter_plot(xdata, y1my2data, xerrdata=xerrdata, yerrdata=y1my2errdata, xlabel=xcol, ylabel=y1col+'-'+y2col)

def xy1y2_plot_single_table(data, xcol, y1col, y2col, xerrcol=None, y1errcol=None, y2errcol=None):
	'''
	plot xcol VS. y1col-y2col
	'''
	xdata = data[xcol]
	y1data = data[y1col]
	y2data = data[y2col]
	y1my2data  = y1data - y2data
	if xerrcol is not None:
		xerrdata = data[xerrcol]
	else:
		xerrdata = None
	if y1errcol is not None and y2errcol is not None:
		y1my2errdata = np.sqrt(data[y1errcol]**2+data[y2errcol]**2)
	else:
		y1my2errdata = None
	xy_simple_scatter_plot(xdata, y1my2data, xerrdata=xerrdata, yerrdata=y1my2errdata, xlabel=xcol, ylabel=y1col+'-'+y2col)





def xy_plot_general(data, xcol, ycol, xerrcol=None, yerrcol=None):
	'''
	plot xcol VS. ycol of data
	'''
	xdata = data[xcol]
	ydata = data[ycol]
	if xerrcol is not None:
		xerrdata = data[xerrcol]
	else:
		xerrdata = None
	if yerrcol is not None:
		yerrdata = data[yerrcol]
	else:
		yerrdata = None
	xy_simple_scatter_plot(xdata, ydata, xerrdata=xerrdata, yerrdata=yerrdata, xlabel=xcol, ylabel=ycol)


def xy_simple_scatter_plot(xdata, ydata, xerrdata=None, yerrdata=None, xlabel='X', ylabel='Y', color='k', alpha=1.0, newaxes=1):
	if newaxes:
		fig = plt.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		ax.errorbar(xdata, ydata, xerr=xerrdata, yerr=yerrdata, fmt='o', color=color, alpha=alpha)
		ax.set_xlabel(xlabel, fontsize=15)
		ax.set_ylabel(ylabel, fontsize=15)
		plt.show()
	else:
		plt.errorbar(xdata, ydata, xerr=xerrdata, yerr=yerrdata, fmt='o', color=color, alpha=alpha)
		plt.xlabel(xlabel, fontsize=15)
		plt.ylabel(ylabel, fontsize=15)

def xy_simple_scatter_plot_2sets(xdata1, ydata1, xdata2, ydata2, xerrdata1=None, yerrdata1=None, xerrdata2=None, yerrdata2=None, label1='DATA1', label2='DATA2', xlabel='X', ylabel='Y'):
	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)
	ax.errorbar(xdata1, ydata1, xerr=xerrdata1, yerr=yerrdata1, fmt='ro', label=label1)
	ax.errorbar(xdata2, ydata2, xerr=xerrdata2, yerr=yerrdata2, fmt='ko', label=label2)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	plt.legend()
	plt.show()


def simple_plot_table_data_four_columns(data, columns, outfile=None):

	fig = plt.figure(figsize=(8,7))
	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)

	for ax,magcol in zip([ax1,ax2,ax3,ax4],columns):
		ax.plot(self.standards[magcol], self.standards['e_'+magcol],'.')
		ax.set_xlabel(magcol)
		ax.set_ylabel('e_'+magcol)
		ax.set_ylim([-0.02, 0.3])

	if outfile is not None:
		plt.savefig(outfile)

	plt.show()

def sources_match_result_stat_plot(input_xys, transformed_xys):

	fig = plt.figure(figsize=(12, 6))
	gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
	ax1 = plt.subplot(gs[0])
	ax1.set_xlabel('$\Delta X$')
	ax1.set_ylabel('$\Delta Y$')
	ax2 = plt.subplot(gs[1])
	ax2.set_xlabel('$\Delta X or \Delta Y$')
	ax2.set_ylabel('N')
	xdiff = input_xys[:,0]-transformed_xys[:,0]
	ydiff = input_xys[:,1]-transformed_xys[:,1]
	ax1.plot(xdiff, ydiff, 'ro')
	ax2.hist(xdiff, bins=20,  alpha=0.5, color='r',  label='X')
	ax2.hist(ydiff, bins=20, alpha = 0.5, color='b', label='Y')
	ax2.legend()
	plt.show()


def __plot_light_curve(lc_info,t0 = 2457000,xylim=False):

	def jd2iso(jd):
		from astropy.time import Time
		t = Time(jd,scale='utc',format='jd')
		t_iso = t.iso
		t_iso_short = t_iso.split()[0]

		return t_iso_short


	colors = ['g','b','r','y','c','k','m']
	symbols = ['o','s','D','*','+','p','H']

	ax = host_subplot(111, axes_class=AA.Axes)

	target_name = lc_info['sn']
	lc_data_dir = lc_info['lc_spec_dir']
	lc_filenames = lc_info['lc']

	cm = {}
	ci = 0
	mi = 0

	t_min = None
	t_max = None
	m_min = None
	m_max = None

	for i,lc_file in enumerate(lc_filenames):
		instrument_band = lc_file.split('_')[1].split('.')[0]
		instrument  = instrument_band.split('-')[0]
		flt = instrument_band.split('-')[1]
		if flt == 'uu':
			flt = 'U'
		if flt == 'bb':
			flt = 'B'
		if flt == 'vv':
			flt = 'V'

		print flt

		if instrument not in cm.keys():
			cm[instrument] = mi
			mi = mi+1
		if flt not in cm.keys():

			cm[flt] = ci
			ci = ci+1

		lc = np.loadtxt(os.path.join(lc_data_dir,lc_file))

		try:
			t     = lc[:,0]-t0
			mag   = lc[:,1]
			magerr= lc[:,2]
		except:
			t = lc[0] - t0
			mag = lc[1]
			magerr = lc[2]

		t_min_this = np.min(t)
		t_max_this = np.max(t)
		m_min_this = np.min(mag)
		m_max_this = np.max(mag)

		if t_min is None or t_min_this < t_min:
			t_min = t_min_this
		if t_max is None or t_max_this > t_max:
			t_max = t_max_this
		if m_min is None or m_min_this < m_min:
			m_min = m_min_this
		if m_max is None or m_max_this > m_max:
			m_max = m_max_this

		plt_hand = __plot_light_curve_single(ax,t,mag,err = magerr,color = colors[np.mod(cm[flt],7)],marker = symbols[np.mod(cm[instrument],7)])


	print t_min,t_max,m_min,m_max

	if xylim:
		plt.xlim([t_min-0.8,t_max+0.8])
		plt.ylim([m_min-0.5,m_max+0.5])

	xticks = plt.xticks()
	ax2 = ax.twin() # ax2 is responsible for "top" axis and "right" axis
	ax2.set_xticks(xticks[0])

	ax2ticklabels = [jd2iso(xticklabel +t0 ) for xticklabel in xticks[0]]

	ax2ticklabels_new = []
	suffix_prior = ''
	for ticklabel in ax2ticklabels:
	 	ticklabel_coms = ticklabel.split('-')
		suffix =ticklabel_coms[0]
		if suffix == suffix_prior:
			ax2ticklabels_new.append(ticklabel_coms[1]+'-'+ticklabel_coms[2])
		else:
			suffix_prior = ticklabel_coms[0]
			ax2ticklabels_new.append(ticklabel_coms[0]+'-'+ticklabel_coms[1]+'-'+ticklabel_coms[2])

	ax2.set_xticklabels(ax2ticklabels_new)

	ax2.axis["right"].major_ticklabels.set_visible(False)

	plt.draw()

	ncols = np.round(len(lc_filenames)/4)
	legends = [lc_filename.split('_')[1].split('.')[0] for lc_filename in lc_filenames]
        plt.legend(legends,loc=0,numpoints=1)
        title_str = 'Light Curve ' + target_name
#        plt.title(title_str)
	plt.rcParams['axes.labelsize'] = 20
	plt.rcParams['font.size'] = 20
        plt.xlabel('JD (+'+str(t0)+')')
        ylabel_str = 'magnitude'
        plt.ylabel(ylabel_str)
        plt.gca().invert_yaxis()
#	ax = plt.gca()
#	print ax.axis['left']
#	ax.axis['left'].major_ticklabels.set_size(18)
#	ax.axis['bottom'].major_ticklabels.set_size(18)
#	ax2.axis['top'].major_ticklabels.set_size(14)

#	ax.xaxis.get_label().set_fontproperties(font)
#	ax.yaxis.get_label().set_fontproperties(font)
#	plt.tick_params(labelsize= 18)
#	ax.xaxis.set_tick_params(labelsize= 28)

	plt.grid()
        plt.show()


def __plot_light_curve_single(ax,t,mag,err=None,color='b'):
	'''

	'''
	if err is not None:
	        lc_plot  = ax.errorbar(t,mag,yerr = err,fmt= color+'o')
	else:
		lc_plot = ax.plot(t,mag,color+'o')

        return lc_plo
