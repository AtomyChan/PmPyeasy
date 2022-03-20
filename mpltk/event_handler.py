import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.image import AxesImage

def __mouse_pick(x,y,err=None,yreverse=True):

    fig, ax1 = plt.subplots(1,1)
    ax1.set_title('click on points to select the point in interest')
    ax1.set_ylabel('x', bbox=dict(facecolor='red'))
    ax1.set_xlabel('y', bbox=dict(facecolor='red'))
    if err is not None:        
    	ax1.errorbar(x,y,yerr=err,fmt= '.',picker=3)
    else:
	err = np.random.random(len(x))/10.
	ax1.errorbar(x,y,yerr=err,fmt= '.',picker=3)

    if yreverse:
	plt.gca().invert_yaxis()
    
    global coorxs,coorys,pick_index
    coorxs = []
    coorys = []
    pick_index = [] 


	
    def onpick1(event):	
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind

            print('onpick1 line:', ind, zip(np.take(xdata, ind), np.take(ydata, ind)))
            pick_index.append(ind[0])
            coorxs.append(np.take(xdata, ind))
	    coorys.append(np.take(ydata, ind))

        elif isinstance(event.artist, Rectangle):
            patch = event.artist
            print('onpick1 patch:', patch.get_path())
        elif isinstance(event.artist, Text):
            text = event.artist
            print('onpick1 text:', text.get_text())

    fig.canvas.mpl_connect('pick_event', onpick1)
    plt.show()

    return coorxs,coorys,pick_index


def mouse_pick(x,y,err=None,yreverse=True):
    xs,ys,indexs = __mouse_pick(x,y,err=err,yreverse=yreverse)
    return xs,ys,indexs

