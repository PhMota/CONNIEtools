# coding: utf-8

import numpy as np
import scipy

def stdCovMed( x, y ):
    X = x - np.median(x)
    Y = y - np.median(y)
    return np.mean(X*Y)
    
def madCovSqr( data0, data1 ):
    return np.median( np.abs( (data0 - np.median(data0))*((data1 - np.median(data1))) ) )

def normFit( data, loc=None, scale=None ):
    return scipy.stats.norm.fit( data, floc=loc, fscale=scale )

def tiles( data, p ):
    m = np.median(data)
    if p - 1 == 0: return [m]
    return tiles(data[data<m], p-1 ) + [m] + tiles(data[data>m], p-1 )

def quartiles(x):
    return tiles(x,2)

def IQR( data, scale=1.34897 ):
    return scipy.stats.iqr( data )/scale

def MAD( data, axis = None, scale=1.4826 ):
    if axis == 0:
        median = np.median( data, axis = 0 )[None,:]
    elif axis == 1:
        median = np.median( data, axis = 1 )[:,None]
    else:
        median = np.median( data )
    return scale * np.median( np.abs( data - median ), axis = axis )


def mean_outliers( data, axis=None, n=2 ):
    med = np.median(data,axis)
    MAD_ = MAD(data,axis)
    return np.mean( data[ np.abs(data-med)<n*MAD_ ] )
