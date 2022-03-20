from __future__ import print_function
import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.image import AxesImage


def plotdata(x,y,err=None,yreverse=True):

    fig = plt.figure(100)	
    ax1 = fig.add_subplot(111)
    ax1.set_title('click on points to select the point in interest')
    ax1.set_ylabel('x', bbox=dict(facecolor='red'))
    ax1.set_xlabel('y', bbox=dict(facecolor='red'))
    if err is not None:        
    	ax1.errorbar(x,y,yerr=err,fmt= 'o',picker=3)
    else:
	ax1.plot(x,y, 'o',picker=3)

    if yreverse:
	ax1.invert_yaxis()
    
    plt.show()
    return fig


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
    else:
	print('can not handle this event type')



def mouse_pick(x,y,err=None,yreverse=True):
    fig = plotdata(x,y,err=err,yreverse=yreverse)
    fig.canvas.mpl_connect('pick_event', onpick1)
    return coorxs,coorys,pick_index


if __name__ == "__main__":

	x = np.random.random(10)
	y = np.random.random(10)
	
	xs,ys,indexs = mouse_pick(x,y,err=None,yreverse=False)
	print(xs,ys,indexs)
