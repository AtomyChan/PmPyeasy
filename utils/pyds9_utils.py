from pyds9 import DS9 

def display_image_with_ds9(fitsfile):

	d = DS9()
	d.set('fits %s'%fitsfile)
	d.set('zoom to fit')
	d.set('scale zscale')

	return d
