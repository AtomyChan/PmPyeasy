import matplotlib.pylab as plt
from common import restore_results
import matplotlib.font_manager as fm


import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


font = fm.FontProperties()
font.set_size(18)





def plot_light_curve(result_dir,target_name,filters,mode,app,phot_method='app',t0 = 0):

	def jd2iso(jd):
		from astropy.time import Time
		t = Time(jd,scale='utc',format='jd')
		t_iso = t.iso
		t_iso_short = t_iso.split()[0]
	
		return t_iso_short


        colors = ['r','b','g','y','c','k','g']
#        ax = plt.subplot(111)
	ax = host_subplot(111, axes_class=AA.Axes)
        for i,flt in enumerate(filters):
	    if phot_method == 'psf':
	        datafile = result_dir + target_name + '_' + mode + '_lc_psf_' + flt +'.dat'	
    	    else:
                datafile = result_dir + target_name + '_' + mode + '_lc_'+str(app)+"_" + flt +'.dat'

            lc = restore_results(datafile)       
	    t     = lc[:,0]-t0
	    mag   = lc[:,1]
	    magerr= lc[:,2]

	    plt_hand = __plot_light_curve(ax,t,mag,err = magerr,color = colors[i])


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

        plt.legend(filters,loc=0,numpoints=1)
        title_str = 'Light Curve ' + target_name
#        plt.title(title_str)
	plt.rcParams['axes.labelsize'] = 20
	plt.rcParams['font.size'] = 20
        plt.xlabel('MJD (+'+str(t0)+')')
        ylabel_str = mode +' magnitude'
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

        return 0 


def __plot_light_curve(ax,t,mag,err=None,color='b'):
	'''

	'''
	if err is not None:
	        lc_plot  = ax.errorbar(t,mag,yerr = err,fmt= color+'o')
	else:
		lc_plot = ax.plot(t,mag,color+'o')

        return lc_plot

