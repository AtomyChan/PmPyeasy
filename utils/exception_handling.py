#! /usr/bin/env python
import numpy as np

def obstime2filename(mapping,obstime,distance=1):
    '''
	distance: default 1min
    '''
    distance = distance/60.0/24.0
    filenames = mapping['name']
    mjds = mapping['obstime']
    mask = np.logical_and(mjds>obstime-distance,mjds<obstime+distance)
    return filenames[mask]





