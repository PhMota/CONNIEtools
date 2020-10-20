#!/usr/bin/python
# coding: utf-8

from __future__ import print_function
import warnings
warnings.filterwarnings("ignore")
import os
import sys
try:
    import numpy as np
except ImportError:
    print('missing module, please run')
    print('module load softwares/python/2.7-gnu-5.3')
    exit(0)

from termcolor import colored
import astropy.io.fits
import scipy.stats
import scipy.ndimage
from glob import glob
import re

import matplotlib
matplotlib.use('gtk3agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable

import Statistics as stats
from Timer import Timer
#from TerminalColor import text
from termcolor import colored

from PrintVar import print_var
np.warnings.filterwarnings("ignore")
import argparse
from collections import OrderedDict


import Simulation
import constants


bias_from_width = { 9340: 450, 8540: 150 }
bin_from_width = { 9340: 5, 8540: 1 }
vbias_from_height = { 900: 70, 900-1: 70, 4*1055: 90, 4*1055-4: 90 }

def MAD( data, axis = None, scale=1.4826 ):
    if axis == 0:
        median = np.nanmedian( data, axis = 0 )[None,:]
    elif axis == 1:
        median = np.nanmedian( data, axis = 1 )[:,None]
    elif axis is None:
        median = np.nanmedian( data, axis = None )
    else:
        raise Exception( 'axis not found', axis )
    return scale * np.nanmedian( np.abs( data - median ), axis = axis )

def crop_var( x, axis = None, mode = None ):
    x2 = crop_mean(x**2, axis = axis, mode=mode)
    m = crop_mean(x, axis = axis, mode=mode)**2
    return x2 - m**2

def crop_std( x, axis = None, mode = None): return np.sqrt( crop_var( x, axis, mode ) )

def crop_mean( x, axis = None ):
    a = np.apply_along_axis( lambda x: crop_mean_1d(x, abs_err=.25), axis, x )
    return a

def crop_mean_1d( x, step = 1, abs_err = .5 ):
    x = np.sort( x.flatten() )
    left, right = 0, len(x)

    mean = np.nanmean(x)
    median = np.nanmedian(x)

    while True:
        if mean > median: right -= step
        else: left += step

        mean = np.nanmean( x[left:right] )
        median = np.nanmedian( x[left:right] )
        if abs(mean - median) < abs_err: break
    return mean, mean

def norm_mean_1d( x, step = 1, abs_err = .25 ):
    y = np.array(x)
    mean = np.nanmean(y)
    median = np.nanmedian(y)
    while True:
        if mean > median: y[ y == np.nanmax(y) ] = np.nan
        else: y[ y == np.nanmin(y) ] = np.nan

        mean = np.nanmean(y)
        median = np.nanmedian(y)
        if abs( mean - median ) < abs_err: break
    return mean, y

def norm_mean_1d2( x, step = 1, abs_err = 0 ):
    from matplotlib import pylab as plt
    mean = [ np.nanmean(x) ]
    median = [ np.nanmedian(x) ]
    diff = [ mean[-1] - median[-1] ]
    ddiff = []
    count = 0
    while True:
        if mean[-1] > median[-1]: x[ x == np.nanmax(x) ] = np.nan
        else: x[ x == np.nanmin(x) ] = np.nan

        mean.append( np.nanmean(x) )
        median.append( np.nanmedian(x) )
        diff.append( mean[-1] - median[-1] )
        ddiff.append( diff[-2] - diff[-1] )
        #print( np.sum(np.isnan(x)), mean[-1], median[-1], diff[-1], diff[-2]-diff[-1] )
        #if abs(diff[-1]) > abs(diff[-2]): break
        if abs( diff[-1] ) < abs_err: break
        if count > 1100: break
        count += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.plot( mean, label='mean' )
    #ax.plot( median, label='median' )
    #ax.plot( diff, label='diff' )
    ax.plot( ddiff, label='ddiff' )
    ax.legend()
    plt.show()
    return mean[-1], x

def norm_mean( x, axis=None, step = 1 ):
    np.apply_along_axis( norm_mean_1d, axis, x )

def outliers2nan_1d( x, step=1, abs_err = .5 ):
    y = np.array(x)
    mean = np.nanmean(y)
    median = np.nanmedian(y)
    while abs( mean - median ) >= abs_err:
        if mean > median: y[ y == np.nanmax(y) ] = np.nan
        else: y[ y == np.nanmin(y) ] = np.nan

        mean = np.nanmean(y)
        median = np.nanmedian(y)
    return y

def outliers2nan( x, axis=None, abs_err = .25 ):
    y = np.apply_along_axis( lambda x: outliers2nan_1d(x, abs_err=abs_err), axis, x )
    #print( 'outliers2nan' )
    #print( 'x', x.shape, np.sum(np.isnan(x)) )
    #print( 'y', y.shape, np.sum(np.isnan(y)) )
    return y

def outliers2nanp_1d( x, pmin = 1e-2 ):
    y = np.array(x)
    mean = np.nanmean(y)
    std = np.nanstd(y)
    median = np.nanmedian(y)
    p = scipy.stats.norm.pdf( y, loc=mean, scale=std )
    while np.nanmin(p) <= pmin:
        y[ np.argwhere( p < pmin ) ] = np.nan
        mean = np.nanmean(y)
        median = np.nanmedian(y)
        std = np.nanstd(y)
        p = scipy.stats.norm.pdf( y, loc=mean, scale=std )
        #print( 'minp', np.sum(np.isnan(y)), np.sum(np.isnan(p)), np.nanmin(p) )
    return y

def outliers2nanp( x, axis=None, pmin = 1e-2 ):
    y = np.apply_along_axis( lambda x: outliers2nanp_1d(x, pmin=pmin), axis, x )
    return y

def outliers2nanPoissonNorm( x, sigma, pmin = 1e-2 ):
    y = np.array(x)

    count = 0
    while True:
        mean = np.nanmean(y)
        median = np.nanmedian(y)
        std = np.nanstd(y)
        mad = MAD(y)

        glambda = mean
        g2lambda = std**2 - sigma**2
        g = g2lambda/glambda
        lambda_ = glambda/g

        p = stats.poisson_norm.pdf( y, loc=median, scale=sigma, mu=lambda_, g=g )
        p /= np.nanmax(p)
        print( 'outliers', ' '.join(['{:.4f}'.format(v) for v in [count, sigma, mean, median, std, mad, glambda, g2lambda, lambda_, g]]), np.nanmin(p), pmin, np.nanmax(y) )
        if np.nanmin(p) > pmin:
            break
        y[ np.nonzero( p < pmin ) ] = np.nan
        print( 'removed', np.sum(np.isnan( y )) )
        count += 1
    return y


def crop_std_1d( x, step = 1, abs_tol = .5 ):
    x = np.sort( x.flatten() )
    left, right = 0, len(x)

    std = np.std(x)
    mad = MAD(x)
    mean = np.nanmean(x)
    median = np.nanmedian(x)

    while abs(std - mad) >= abs_tol:
        if mean > median: right -= step
        else: left += step

        mean = np.nanmean( x[left:right] )
        median = np.nanmedian( x[left:right] )
        std = np.std(x[left:right])
        mad = MAD(x[left:right])
        if (right-left) < .7*len(x):
            print('did not converge', mean, median, std, mad, right-left, len(x) )
            break
    print( '\tcrop var', std, mad )
    return std

def crop_mean_1d_partition( x, partitions = 10 ):
    mean = np.nanmean(x)
    median = np.nanmedian(x)

    number_to_delete = 100
    display( [x], 'E', nbins = 100, log=True )
    while abs(mean - median) >= .5:
        bins = np.linspace( x.min(), x.max()+1, partitions )
        bin_indices = np.digitize( x, bins )
        counts, _ = np.histogram( x, bins )
        #print( 'counts before', counts )
        #print( 'counts before', bins )
        indices_to_delete = []
        for i, count in enumerate( counts ):
            indices_in_bin = np.argwhere( bin_indices == (i+1) ).T[0]
            #print( 'indices', i, len(indices_in_bin), counts, number_to_delete )
            #if len(indices_in_bin) == count:
                #indices_to_delete.extend( indices_in_bin )
            if len(indices_in_bin) > 0:
                to_delete = np.random.choice( indices_in_bin, min(len(indices_in_bin), number_to_delete), replace = False )
                #print( len(to_delete ), count )
                indices_to_delete.extend( to_delete )
        x = np.delete(x, indices_to_delete )
        counts, _ = np.histogram( x, bins )

        #mean = np.mean( x )
        #median = np.median( x )
        print( 'stats', np.nanmean( x ), np.nanmedian( x ), np.std( x ), MAD( x ), len(x) )
        display( [x], 'E', nbins = 100, log=True )
        #raw_input()
    return mean

def stats_robust_2(self, gain, lambda_, sigma, tol = 1e-3, partitions = 10):
    x = self.flatten()

    compute_gain = lambda x: (np.std(x)**2 - sigma**2)/np.nanmean(x)
    compute_lambda = lambda x: np.nanmean(x)**2/(np.std(x)**2 - sigma**2)

    mean = [np.nanmean(x)]
    median = [np.nanmedian(x)]
    std = [np.std(x)]
    mad = [MAD(x)]
    gains = [ compute_gain(x) ]
    lambdas = [ compute_lambda(x) ]

    goal_std = np.sqrt( gain**2*lambda_ + sigma**2 )
    goal_mean = gain * lambda_

    while True:
        bins = np.linspace( x.min(), x.max(), 10 )
        bin_indices = np.digitize( x, bins )
        counts, _ = np.histogram( x, bins )
        print( 'counts before', counts )
        number_to_delete = np.min( counts[ counts > 0 ] )
        indices_to_delete = []
        for i, ind in enumerate( bins ):
            indices_in_bin = np.argwhere( bin_indices == i )
            if len(indices_in_bin) > 0:
                to_delete = np.random.choice( indices_in_bin, number_to_delete, replace = False )
                print( len(to_delete ) )
                indices_to_delete.extend( indices_in_bin[ to_delete ] )
        x = np.delete(x, indices_to_delete )
        counts, _ = np.histogram( x, bins )
        print( 'counts after', counts )

        mean.append( np.nanmean( x ) )
        median.append( np.nanmedian( x ) )
        std.append( np.std( x ) )
        mad.append( MAD( x ) )
        gains.append( compute_gain( x ) )
        lambdas.append( compute_lambda( x ) )

        print( 'expected', goal_mean, goal_std )
        print( 'stats', np.nanmean(mean), np.nanmean(median), np.nanmean(std), np.nanmean(mad) )
        print( 'vals', gain, lambda_ )
        print( 'stats2', np.nanmean(gains), np.nanmean(lambdas) )
        print()

    return np.nanmean( mean ), np.nanmean( median ), np.nanmean( std ), np.nanmean( mad ), x


def list_params( fits, param = 'OHDU' ):
    return [ fit.header[ param ] for fit in fits ]

def get_imageHDU( fits, ohdu ):
    index = list_params( fits ).index( ohdu )
    return fits[ index ]

def statistics( data ):
    return (float(np.nanmean(data)),
            float(np.nanmedian(data)),
            float(crop_mean_1d(data, step=10)),
            float(np.std(data)),
            float( MAD(data) ),
            float(crop_std_1d(data, step=10, abs_tol = 1)),
            )

def read_image_from_fits( path, ohdu ):
    paths = glob( expr )[::-1]
    imagesHDU = []
    for path in paths:
        fits_dir = astropy.io.fits.open( path )
        imagesHDU.append( get_imageHDU( fits_dir, ohdu ) )
    return imagesHDU

def label_clusters( condition_image ):
    #s33 = [[1,1,1], [1,1,1], [1,1,1]]
    s33 = [[0,1,0], [1,1,1], [0,1,0]]
    return scipy.ndimage.label( condition_image, structure=s33 )[0]

def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    assert h % nrows == 0, "{} rows is not evenly divisble by {}".format(h, nrows)
    assert w % ncols == 0, "{} cols is not evenly divisble by {}".format(w, ncols)
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def mean_err( x ):
    return np.nanmean(x), np.nanstd(x)/np.sqrt(len(x))

class Part:
    '''
    deals with the image HDU parts
    '''
    ccd_width = 4120
    ccd_height = 4150

    def __init__( self, imageHDU ):
        self.header = dict( imageHDU.header )
        print( 'image.shape', imageHDU.data.shape )
        self.all = Section( imageHDU.data )

        try:
            if 'BIASW' in self.header:
                self.width = imageHDU.data.shape[1]
                self.bias_width = imageHDU.header['BIASW']
                self.bias_height = imageHDU.header['BIASH']
                self.rebin = imageHDU.header['REBIN']
                self.has_right = False
            else:
                self.parse_shape( imageHDU.data.shape )

            # important to remove the first line, but not doing it now
            self.data = Section( imageHDU.data[ None:None, None:(self.width-self.bias_width) ] )
            self.bias = Section( imageHDU.data[ None:None, (self.width-self.bias_width): self.width ] )

            if self.has_right:
                self.dataR = Section( imageHDU.data[:,::-1][ None:None, None:(self.width-self.bias_width) ] )
                self.biasR = Section( imageHDU.data[:,::-1][ None:None, (self.width-self.bias_width):self.width ] )
        except:
            return

    def remove_outliers_from_bias( self, pmin = 1e-5 ):
        self.bias = Section(outliers2nanp( self.bias, axis=1, pmin=pmin ))
        if self.has_right:
            self.biasR = Section(outliers2nanp( self.biasR, axis=1, pmin=pmin ))

    def correct_lines( self, func = np.nanmedian, smoothening_length = 0 ):
        if func is None:
            self.data -= self.bias.get_global_bias( np.nanmedian )
            self.bias -= self.bias.get_global_bias( np.nanmedian )
            if self.has_right:
                self.dataR -= self.biasR.get_global_bias( np.nanmedian )
                self.biasR -= self.biasR.get_global_bias( np.nanmedian )
        else:
            lines_correction = self.bias.get_lines_bias( func, smoothening_length )
            self.data -= lines_correction
            self.bias -= lines_correction
            if self.has_right:
                lines_correctionR = self.biasR.get_lines_bias( func, smoothening_length )
                self.dataR -= lines_correctionR
                self.biasR -= lines_correctionR
        return self

    def parse_shape( self, shape ):
        height, width = shape
        number_of_amplifiers = width//constants.ccd_width
        #print( 'number of amplifiers', number_of_amplifiers )
        side_width = width/number_of_amplifiers
        bias_width = side_width - constants.ccd_width
        #print( 'bias_width', side_width - constants.ccd_width )
        self.parse_height( shape[0] )
        self.parse_width( shape[1] )

    def parse_height( self, height ):
        '''
        Attempts to guess what are the bias height and if the image is rebinned from the image height. If it fails if falls back to 90 of height and no rebinning
        '''
        try:
            self.bias_height, self.rebin = {
                900: [74, 5],
                1055: [90, 1],
                4220: [90, 1],
                }[height]
        except KeyError:
            print( 'image height {} is not recognized. Assuming rebin 1 and bias height of 90'.format(height) )
            self.bias_height, self.rebin = 90, 1
        return

    def parse_width( self, width ):
        '''
        Attempts to guess what are the bias width, image width and if there is a right side image from the width of the input image. If it fails, it falls back to 150 bias width, width and no righ side
        '''
        try:
            self.bias_width, self.has_right, self.width = {
                (Part.ccd_width+450)*2: [450, True, (Part.ccd_width+450)],
                (Part.ccd_width+150)*2: [150, True, (Part.ccd_width+150)],
                (Part.ccd_width+450): [450, False, (Part.ccd_width+450)],
                (Part.ccd_width+150): [150, False, (Part.ccd_width+150)],
                (Part.ccd_width-8+450): [450, False, (Part.ccd_width+450-8)],
                (Part.ccd_width-8+150): [150, False, (Part.ccd_width+150-8)],
                (Part.ccd_width-8+550): [550, False, (Part.ccd_width+550-8)],
                4570: [150, False, 4570],
                }[width]
        except KeyError:
            print( 'image width {} is not recognized. Assuming bias width of 150 and no right side'.format(width) )
            self.bias_width, self.has_right = 150, False
            self.width = width
        return

    def get_median_params_by_line( self ):
        data = self.data.get_half()
        bias = self.bias
        height, width = data.shape

        _, osw = bias.shape
        rebin = self.rebin
        ret = OrderedDict()
        ret['height'] = height
        ret['width'] = width
        ret['osw'] = osw
        ret['rebin'] = rebin

        ret['bias_mu'], ret['bias_mu_err'] = mean_err( np.nanmedian( bias, axis=1 ) )
        ret['bias_sigma'], ret['bias_sigma_err'] = mean_err( MAD( bias, axis=1 ) )

        ret['data_mu'], ret['data_mu_err'] = mean_err( np.nanmedian( data, axis=0 ) )
        ret['data_sigma'], ret['data_sigma_err'] = mean_err( MAD( data, axis=0) )

        ret['g_lambda'], ret['g_lambda_err'] = mean_err( np.median( data, axis=1) - np.median( bias, axis=1) )
        ret['g2_lambda'], ret['g2_lambda_err'] = mean_err( MAD( data, axis=1)**2 - MAD( bias, axis=1)**2 )

        ret['g'] = ret['g2_lambda']/ret['g_lambda']
        ret['lambda'] = ret['g_lambda']/ret['g']
        return ret

class Image:
    def __init__( self, parts ):
        for part in parts:
            try:
                self.all = np.concatenate( [ self.all, part.all ], axis=0 )
            except AttributeError:
                self.all = part.all
                self.header = part.header
            try:
                try:
                    self.data = np.concatenate( [ self.data, part.data ], axis=0 )
                    self.bias = np.concatenate( [ self.bias, part.bias ], axis=0 )
                    if self.has_right:
                        self.dataR = np.concatenate( [ self.dataR, part.dataR ], axis=0 )
                        self.biasR = np.concatenate( [ self.biasR, part.biasR ], axis=0 )
                except AttributeError:
                    self.has_right = part.has_right
                    self.rebin = part.rebin
                    self.header = part.header
                    self.bias_height = part.bias_height
                    self.data = part.data
                    self.bias = part.bias
                    if self.has_right:
                        self.dataR = part.dataR
                        self.biasR = part.biasR
            except:
                pass

        try:
            self.vbias = self.data[ None:self.bias_height, : ].view(Section)
            self.vbias_err = 1

            self.data = self.data[ self.bias_height:None, : ].view(Section)
            self.data_err = 1

            self.dbias = self.bias[ None:self.bias_height, : ].view(Section)
            self.dbias_err = 1

            self.bias = self.bias[ self.bias_height:None, : ].view(Section)
            self.bias_err = 1

            if self.has_right:
                self.vbiasR = self.dataR[ None:self.bias_height, : ].view(Section)
                self.dataR = self.dataR[ self.bias_height:None, : ].view(Section)
                self.dbiasR = self.biasR[ None:self.bias_height, : ].view(Section)
                self.biasR = self.biasR[ self.bias_height:None, : ].view(Section)
        except:
            pass

    def remove_outliers_from_vbias(self, pmin=1e-5):
        self.vbias = Section(outliers2nanp( self.vbias, axis=0, pmin=pmin ))
        if self.has_right:
            self.vbiasR = Section(outliers2nanp( self.vbiasR, axis=0, pmin=pmin ))

    def remove_outliers_from_vbias_second_pass(self, pmin=1e-5):
        self.vbias = Section(outliers2nanp( self.vbias, axis=1, pmin=pmin ))
        if self.has_right:
            self.vbiasR = Section(outliers2nanp( self.vbiasR, axis=1, pmin=pmin ))

    def correct_cols( self, func = np.nanmedian, smoothening_length = 0 ):
        if not func is None:
            data_cols_correction = self.vbias.get_cols_bias( func, smoothening_length )
            bias_cols_correction = self.dbias.get_cols_bias( func, smoothening_length )
            self.data = self.data - data_cols_correction + np.nanmean(data_cols_correction[ data_cols_correction.shape[0]/2: ])
            self.vbias = self.vbias - data_cols_correction + np.nanmean(data_cols_correction[ data_cols_correction.shape[0]/2: ])
            self.bias = self.bias - bias_cols_correction + np.nanmean( bias_cols_correction )
            self.dbias = self.dbias - bias_cols_correction + np.nanmean( bias_cols_correction )
            if self.has_right:
                data_cols_correctionR = self.vbiasR.get_cols_bias( func, smoothening_length )
                self.dataR -= data_cols_correctionR
                self.vbiasR -= data_cols_correctionR
                bias_cols_correctionR = self.dbiasR.get_cols_bias( func, smoothening_length )
                self.biasR -= bias_cols_correctionR
                self.dbiasR -= bias_cols_correctionR
        return self

    def correct_cols2( self, func = np.nanmedian, smoothening_length = 0 ):
        data_cols_correction = self.data.get_cols_bias( func, smoothening_length )

        if not func is None:
            data_cols_correction = self.vbias.get_cols_bias( func, smoothening_length )
            bias_cols_correction = self.dbias.get_cols_bias( func, smoothening_length )
            self.data = self.data - data_cols_correction + np.nanmean(data_cols_correction[ data_cols_correction.shape[0]/2: ])
            self.vbias = self.vbias - data_cols_correction + np.nanmean(data_cols_correction[ data_cols_correction.shape[0]/2: ])
            self.bias = self.bias - bias_cols_correction + np.nanmean( bias_cols_correction )
            self.dbias = self.dbias - bias_cols_correction + np.nanmean( bias_cols_correction )
            if self.has_right:
                data_cols_correctionR = self.vbiasR.get_cols_bias( func, smoothening_length )
                self.dataR -= data_cols_correctionR
                self.vbiasR -= data_cols_correctionR
                bias_cols_correctionR = self.dbiasR.get_cols_bias( func, smoothening_length )
                self.biasR -= bias_cols_correctionR
                self.dbiasR -= bias_cols_correctionR
        return self


    def get_params( self, mode = None, half = False, gain = 1, remove_hits = False, right_side = False, **kwargs ):
        ret = []
        if right_side:
            data = self.dataR
            bias = self.biasR
            vbias = self.vbiasR
            dbias = self.dbias
        else:
            data = self.data
            bias = self.bias
            vbias = self.vbias
            dbias = self.dbias

        if half:
            data = data.get_half()
            vbias = vbias.get_half()

        height, width = data.shape

        _, osw = bias.shape
        osh,_ = vbias.shape
        rebin = self.rebin
        ret.append( ['height', height] )
        ret.append( ['width', width] )
        ret.append( ['osw', osw] )
        ret.append( ['osh', osh] )
        ret.append( ['rebin', rebin] )

        if 'GAIN' in self.header:
            gain = self.header['GAIN']

        if mode == 'fit':
            a, mu, sigma = self.bias.fit_norm_binned()[0]
            print( 'sigma', sigma )
            a, loc, lamb = data.fit_poisson_norm_binned( gain=gain, sigma=sigma )[0]
            return { 'mu':mu, 'sigma':sigma, 'lambda': lamb/factor, 'g': gain }

        if mode == 'mle':
            mu, sigma = self.bias.mle_norm()
            g_lamb = np.nanmean( np.nanmedian( data, axis=1) - np.nanmedian( self.bias, axis=1) )
            loc, lamb = data.mle_poisson_norm( gain=gain, sigma=sigma, p0 = [0, g_lamb/gain] )
            return { 'mu':mu, 'sigma':sigma, 'lambda': lamb/factor, 'g': gain }

        raise Exception( 'mode not recognized {}'.format(mode) )

    def save_pdf( self, output ):
        data = np.concatenate([self.vbias, self.data], axis=0)
        bias = np.concatenate([self.dbias, self.bias], axis=0)
        image = Section( np.concatenate([data, bias], axis=1) )
        image.save_pdf( output )

    def save_sections( self, output ):
        for sec in ['data', 'bias', 'vbias', 'dbias']:
            getattr( self, sec ).save_pdf( '{}_{}'.format(output, sec), title = sec )
            if self.has_right:
                getattr( self, sec+'R' ).save_pdf( '{}_{}'.format(output, sec+'R'), title = sec+'R' )

    def histograms( self ):
        from matplotlib import pylab as plt

        bias = self.bias.get_binned_distribution()
        vbias = self.vbias.get_binned_distribution()
        data = self.data.get_binned_distribution()

        #_, crop_bias = norm_mean_1d( self.bias )
        #hist = crop_bias.get_binned_distribution()
        #print( 'pixels', np.sum(bias[1]), np.sum(hist[1]) )
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.step( bias[0], bias[1]/float(np.sum(bias[1])), where='mid', label='bias' )
        #ax.step( hist[0], hist[1]/float(np.sum(hist[1])), where='mid', label='crop' )


        xlim = ax.get_xlim()
        print('xlim',xlim)
        #ax.step( vbias[0], vbias[1]/float(np.sum(vbias[1])), where='mid', label='vbias' )
        #ax.step( data[0], data[1]/float(np.sum(data[1])), where='mid', label='data' )

        ax.set_xlim((xlim[0],-xlim[0]))
        ax.set_yscale('log')
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        plt.show()

class Section( np.ndarray ):
    def __new__( cls, input_array ):
        if type(input_array) is float: return float(input_array)
        obj = np.asarray(input_array).astype(float).view(cls)
        return obj

    def save_fits(self, output):
        primary = astropy.io.fits.PrimaryHDU()
        fits_image = astropy.io.fits.ImageHDU( self )
        astropy.io.fits.HDUList([primary, fits_image]).writeto( output, overwrite=True )

    def get_half( self ):
        return self[:, self.shape[1]/2:]

    #def mean(self, axis = None):
        #print( self.shape )
        #print( type(self) )
        #return np.mean( self, axis=axis )

    #def std(self, axis = None):
        #return np.std( self, axis=axis )

    def mean_of_lines( self ): return np.nanmedian(self, axis=1)
    def mean_of_halflines( self ): return np.nanmedian(self.get_half(), axis=1)
    def mean_of_cols( self ): return np.nanmedian(self, axis=1)

    def std_of_lines( self ): return np.nanstd(self, axis=1)
    def std_of_halflines( self ): return np.nanstd(self.get_half(), axis=1)
    def std_of_cols( self ): return np.nanstd(self, axis=1)

    def convolve( self, function=lambda x: np.nansum(x), footprint=None, size=(7,7), origin=0 ):
        image = scipy.ndimage.filters.generic_filter(self, function, size=size, footprint=footprint, output=None, mode='constant', cval=np.nan, origin=origin, extra_arguments=(), extra_keywords=None)
        return Section(image)

    def median( self, axis = None ):
        return np.nanmedian( self, axis=axis )

    def MAD( self, axis = None ):
        return MAD( self, axis=axis )

    def get_global_bias( self, func = np.nanmedian ):
        return func( self )

    def get_lines_bias( self, func = np.nanmedian, smoothening_length = 0 ):
        lines_bias = func( self, axis = 1 )
        if smoothening_length > 2:
            lines_bias = scipy.ndimage.filters.uniform_filter1d( lines_bias, size = smoothening_length, mode='nearest' )
        return lines_bias[:,None]

    def get_cols_bias( self, func = np.nanmedian, smoothening_length = 0 ):
        cols_bias = func( self, axis = 0 )
        if smoothening_length > 2:
            cols_bias = scipy.ndimage.filters.uniform_filter1d( cols_bias, size = smoothening_length, mode='nearest' )
        return cols_bias[None,:]

    def get_binned_distribution( self, binsize = 1):
        bins = np.arange( np.nanmin(self), np.nanmax(self), binsize)
        distribution, _ = np.histogram(self, bins)
        return .5*(bins[:-1] + bins[1:]), distribution

    def fit_binned( self, func, p0, bounds, binsize = 1):
        x, y = self.get_binned_distribution( binsize )
        p, pcov = scipy.optimize.curve_fit( func, x, y, p0=p0, bounds=bounds )
        return p, pcov

    def fit_norm_binned( self, binsize = 1):
        func = lambda x, amplitude, mu, sigma: amplitude * scipy.stats.norm.pdf( x, mu, sigma )
        p0 = [ len(self), self.median(), self.MAD() ]
        bounds = [(0, -np.inf, 0), (np.inf, np.inf, np.inf)]
        return self.fit_binned( func, p0, bounds, binsize )

    def fit_poisson_norm_binned( self, sigma, gain, binsize = 1):
        def func( x, amplitude, loc, lamb ):
            print( 'p', amplitude, loc, lamb )
            ret = amplitude * stats.poisson_norm.pdf( x, loc=loc, scale=sigma, mu=lamb, g=gain )
            return ret
        p0 = [ len(self), np.nanmedian(self), 0.1 ]
        bounds = [(0, -np.inf, 1e-3), (np.inf, np.inf, np.inf)]
        return self.fit_binned( func, p0, bounds, binsize )

    def mle( self, negloglikelihood, p0, bounds ):
        params = scipy.optimize.minimize(
            negloglikelihood,
            x0=p0,
            bounds=bounds,
            ).x
        return params

    def mle_norm( self ):
        data = self[~np.isnan(self)].flatten()
        def negloglikelihood( p ):
            return -np.sum( np.log( scipy.stats.norm.pdf(data, p[0], p[1]) ) )
        p0 = [np.nanmedian(data), MAD(data)]
        bounds = [(-np.inf,np.inf), (1e-3, np.inf)]
        return self.mle( negloglikelihood, p0, bounds )

    def mle_poisson_norm(self, sigma, gain, p0):
        data = self[~np.isnan(self)].flatten()
        def negloglikelihood( p ):
            print( 'p', p )
            return -np.sum( np.log( stats.poisson_norm.pdf(data, loc=p[0], scale=sigma, g=gain, mu=p[1], tol=1e-3) ) )
        bounds = [(-np.inf,np.inf), (1e-5, np.inf)]
        return self.mle( negloglikelihood, p0, bounds )

    #def extract_hits( self, mode, **kwargs ):
        #if mode == 'cluster':
            #return self._extract_clusters_( **kwargs )

    def get_clusters( self, threshold, border, verbose=0 ):
        labeled_clusters = label_clusters( self >= threshold )
        if verbose:
            print( 'number of pixels labeled', np.count_nonzero(labeled_clusters) )
            print( 'number of clusters above threshold', labeled_clusters.max() )
        is_cluster = labeled_clusters > 0
        distances_to_cluster = scipy.ndimage.distance_transform_edt( is_cluster == False )
        labeled_clusters = label_clusters( distances_to_cluster <= border*np.sqrt(2) )
        if verbose:
            print( 'number of pixels labeled', np.count_nonzero(labeled_clusters) )
            print( 'number of clusters with border', labeled_clusters.max() )
        return labeled_clusters, distances_to_cluster

    def get_background( self, threshold, border ):
        labeled_clusters, _ = self.get_clusters( threshold, border )
        return self[ labeled_clusters == 0 ]

    def remove_hits( self, threshold, border ):
        labeled_clusters, _ = self.get_clusters( threshold, border )
        data = Section( np.where( labeled_clusters == 0, self, np.nan ) )
        #print('hits removed', len(self.flatten()), len(self[labeled_clusters == 0]) )
        return data

    def get_clusters_noise( self, threshold, border, noise, verbose=0 ):
        p_value_image = .5*erf( self/sqrt(2)/noise )
        z_value_image = sqrt(2)*erfi(1 - 2*p_value_image)

        labeled_clusters = label_clusters( self >= threshold )
        if verbose:
            print( 'number of pixels labeled', np.count_nonzero(labeled_clusters) )
            print( 'number of clusters above threshold', labeled_clusters.max() )
        is_cluster = labeled_clusters > 0
        distances_to_cluster = scipy.ndimage.distance_transform_edt( is_cluster == False )
        labeled_clusters = label_clusters( distances_to_cluster <= border*np.sqrt(2) )
        if verbose:
            print( 'number of pixels labeled', np.count_nonzero(labeled_clusters) )
            print( 'number of clusters with border', labeled_clusters.max() )
        return labeled_clusters, distances_to_cluster

    def extract_hits_noise( self, threshold, border, verbose=0 ):
        labeled_clusters, distances_to_cluster = self.get_clusters( threshold, border, verbose )

        indices = np.unique(labeled_clusters)[1:]
        if verbose:
            print( 'extracted {}'.format(len(indices)) )
        list_of_clusters = scipy.ndimage.labeled_comprehension(
            self,
            labeled_clusters,
            index = indices,
            func = lambda e, p: [e, p],
            out_dtype=list,
            default=-1,
            pass_positions=True
            )
        levels = scipy.ndimage.labeled_comprehension(
            distances_to_cluster,
            labeled_clusters,
            index = indices,
            func = lambda v: v,
            default=-1,
            out_dtype=list
            )
        def process( cluster ):
            ei, level = cluster
            x, y = np.unravel_index(ei[1], self.shape)
            return ei[0], x, y, level.astype(int)
        list_of_clusters = map( process, zip(list_of_clusters, levels) )

        dtype = [
            ('ePix', object),
            ('xPix', object),
            ('yPix', object),
            ('level', object),
            ]
        ret = np.array( list_of_clusters, dtype = dtype ).view(np.recarray)
        return ret

    def extract_hits( self, threshold, border, verbose=0 ):
        labeled_clusters, distances_to_cluster = self.get_clusters( threshold, border, verbose )

        indices = np.unique(labeled_clusters)[1:]
        if verbose:
            print( 'extracted {}'.format(len(indices)) )
        if len(indices) == 0:
            print( 'nothing was extracted' )
            print( 'max value in image', np.amax(self) )
            exit(0)
        list_of_clusters = scipy.ndimage.labeled_comprehension(
            self,
            labeled_clusters,
            index = indices,
            func = lambda e, p: [e, p],
            out_dtype=list,
            default=-1,
            pass_positions=True
            )
        levels = scipy.ndimage.labeled_comprehension(
            distances_to_cluster,
            labeled_clusters,
            index = indices,
            func = lambda v: v,
            default=-1,
            out_dtype=list
            )
        def process( cluster ):
            ei, level = cluster
            x, y = np.unravel_index(ei[1], self.shape)
            return ei[0], x, y, level.astype(int)
        list_of_clusters = map( process, zip(list_of_clusters, levels) )

        dtype = [
            ('ePix', object),
            ('xPix', object),
            ('yPix', object),
            ('level', object),
            ]
        ret = np.array( list_of_clusters, dtype = dtype ).view(np.recarray)
        return ret

    def spectrum( self, output, gain, lambda_, sigma, binsize = 2 ):
        from matplotlib import pylab as plt

        median = np.nanmedian( self )
        mean = self.mean()
        std = self.std()
        #mean_robust, _, std_robust, _, x = self.stats_robust(gain=gain, lambda_=lambda_, sigma=sigma, tol=1e-3, fraction=0.001 )
        #mean_robust, _, std_robust, _, x = self.stats_robust_2(gain=gain, lambda_=lambda_, sigma=sigma, tol=1e-3, partitions = 10 )
        fig = plt.figure()
        ax = fig.add_subplot(111)
        bins = np.arange( self.min(), x.max()*2, binsize )
        ax.hist( self.flatten(), bins = bins, label = 'pix', histtype = 'step' )
        ax.hist( x, bins = bins, label = 'cropped', histtype = 'step' )
        ymin, ymax = ax.get_ylim()
        ax.axvline( median, ymin, ymax, color='r', label = 'median %.3f' % median )
        ax.axvline( mean, ymin, ymax, color='g', label = 'mean %.3f' % mean )
        ax.axvline( std, ymin, ymax, color='m', label = 'std %.3f' % std )
        ax.axvline( mean_robust, ymin, ymax, color='k', label = 'mean_robust %.3f' % mean_robust )
        ax.axvline( std_robust, ymin, ymax, color='b', label = 'std_robust %.3f' % std_robust )
        ax.axvline( std_robust, ymin, ymax, color='b', label = 'gain %.3f' % ((std_robust**2 - sigma**2)/mean_robust) )
        ax.axvline( std_robust, ymin, ymax, color='b', label = 'lambda %.3f' % (mean_robust**2/(std_robust**2 - sigma**2)) )

        ax.set_yscale('log')
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        fig.savefig( output, bbox_extra_artists = (legend,), bbox_inches='tight')

    def save_pdf( self, output, title=None ):
        from matplotlib import pylab as plt

        if not output.endswith('.pdf'): output += '.pdf'
        with Timer( 'saved ' + output ):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if title: ax.set_title(title)
            im = np.log(self - np.nanmin(self) + 1)

            cmap = matplotlib.cm.Blues
            #cmap.set_bad(color='red')
            cmap.set_bad('red',0)
            try:
                ax.imshow( im, cmap=cmap, origin='lower', vmin=np.nanmin(im), vmax=np.nanmax(im) )
            except:
                print( 'im.shape', im.shape, self.shape )
                exit()
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            fig.savefig( output, bbox_inches='tight', pad_inches=0 )
        return

    def save_projection( self, output, title, axis, no_max=False, no_mean=False, no_min=False, do_regress=False ):
        from matplotlib import pylab as plt

        if not output.endswith('.pdf'): output += '.pdf'
        with Timer( 'saved ' + output ):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(title)
            maxs = np.nanmax( self, axis=axis )
            mins = np.nanmin( self, axis=axis )
            stds = np.nanstd( self, axis=axis )
            means = np.nanmean( self, axis=axis )
            medians = np.nanmedian( self, axis=axis )
            mads = MAD( self, axis=axis )

            if not no_max:
                ax.plot( maxs, '.', label='max' )
            ax.plot( medians, '.', label='median' )
            ax.fill_between( range(len(medians)), medians-mads, medians+mads, label='MAD', alpha=.5 )
            if not no_mean:
                ax.plot( means, '.', label='mean' )
                ax.fill_between( range(len(means)), means-stds, means+stds, label='std', alpha=.5 )
            if not no_min:
                ax.plot( mins, '.', label='min' )

            if do_regress:
                x = np.array(range(len(medians)))
                y = medians
                mask = ~np.isnan(y)
                a, b = scipy.stats.linregress(x[mask],y[mask])[:2]
                ax.plot( x[mask], a*x[mask]+b, '-', label='median\na={:.4e}\nb={:.4e}'.format(a, b) )
                x = np.array(range(len(mads)))
                y = mads
                mask = ~np.isnan(y)
                a, b = scipy.stats.linregress(x[mask],y[mask])[:2]
                ax.plot( x[mask], a*x[mask]+b, '-', label='mad\na={:.4e}\nb={:.4e}'.format(a, b) )

            ax.grid()
            ax.legend()
            fig.savefig( output, bbox_inches='tight', pad_inches=0 )

    def save_spectrum( self, output, title, binsize=1 ):
        from matplotlib import pylab as plt

        if not output.endswith('.pdf'): output += '.pdf'
        with Timer( 'saved ' + output ):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(title)

            mean = np.nanmean( self )
            std = np.nanstd( self )
            median = np.nanmedian( self )
            mad = MAD( self )

            w, center = stats.FWHM( self )

            x, y = self.get_binned_distribution( binsize=binsize )
            #halfmax = np.max(y)/2.

            #w = 0
            #w2 = 0
            #f = .01
            #while True:
                #xwidth = x[ abs(y - halfmax) < f*np.max(y) ]
                #if len(xwidth) == 0:
                    #xwidth *= 1.2
                    #continue
                #print_var( 'xwidth', locals() )
                #xwidth_left = xwidth[ xwidth < median ]
                #xwidth_right = xwidth[ xwidth > median ]
                #if len(xwidth_left) == 0 or len(xwidth_right) == 0:
                    #xwidth *= 1.2
                    #continue
                #w = abs( np.mean( xwidth_left ) - np.mean( xwidth_right ) )
                #w2 = abs( np.mean( xwidth_left ) + np.mean( xwidth_right ) )/2.
                #break
            ax.step( x, y, where='mid',
                    label='mean={mean:.4}\nmedian={median:.4}\nMHM={mhm:.4}\nstd={std:.4}\nMAD={mad:.4}\nFWHM={fwhm:.4}\nN={N}'.format(mean=mean, median=median, std=std, mad=mad,N=len(self.flatten()), fwhm = w/2.355, mhm=center)
                    )
            ax.set_yscale('log')
            ax.legend()
            fig.savefig( output, bbox_inches='tight', pad_inches=0 )

    #def display( self, mode, delta = None ):
        #from matplotlib.colors import LogNorm
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #im = -self
        #if mode == 'xy':
            #vmin = None
            #vmax = None
            #if delta:
                #vmin = np.nanmedian(datum)-delta
                #vmax = np.nanmedian(datum)+delta
            #ax.imshow( im, cmap='Blues', vmin=vmin, vmax=vmax)
        #if mode == 'x':
            #ax.plot(self[0,:], '.')
            #ax.plot(self[-1,:], '.')
            #ax.set_ylim( (np.nanmedian(self)-delta,np.nanmedian(self)+delta) )
        #if mode == 'y':
            #ax.plot(self[:, 0], '.')
            #ax.plot(self[:, -1], '.')
            #ax.set_ylim( (np.nanmedian(self)-delta,np.nanmedian(self)+delta) )
        #if mode == 'flat':
            #ax.plot(self, '.')
        #plt.show()


def get_median_params_by_line( data, bias ):
    height, width = data.shape

    _, osw = bias.shape
    ret = OrderedDict()
    ret['height'] = height
    ret['width'] = width
    ret['osw'] = osw

    ret['bias_mu'], ret['bias_mu_err'] = mean_err( np.nanmedian( bias, axis=1 ) )
    ret['bias_sigma'], ret['bias_sigma_err'] = mean_err( MAD( bias, axis=1 ) )

    ret['data_mu'], ret['data_mu_err'] = mean_err( np.nanmedian( data, axis=0 ) )
    ret['data_sigma'], ret['data_sigma_err'] = mean_err( MAD( data, axis=0) )

    ret['g_lambda'], ret['g_lambda_err'] = mean_err( np.nanmedian( data, axis=1) - np.nanmedian( bias, axis=1) )
    ret['g2_lambda'], ret['g2_lambda_err'] = mean_err( MAD( data, axis=1)**2 - MAD( bias, axis=1)**2 )

    ret['g'] = ret['g2_lambda']/ret['g_lambda']
    ret['lambda'] = ret['g_lambda']/ret['g']
    return ret

def get_mean_params_by_line( data, bias ):
    height, width = data.shape

    _, osw = bias.shape
    ret = OrderedDict()
    ret['height'] = height
    ret['width'] = width
    ret['osw'] = osw

    ret['bias_mu'], ret['bias_mu_err'] = mean_err( np.nanmean( bias, axis=1 ) )
    ret['bias_sigma'], ret['bias_sigma_err'] = mean_err( np.nanstd( bias, axis=1 ) )

    ret['g_lambda'], ret['g_lambda_err'] = mean_err( np.nanmedian( data, axis=1 ) - np.nanmean( bias, axis=1 ) )
    ret['g2_lambda'], ret['g2_lambda_err'] = mean_err( MAD( data, axis=1 )**2 - np.nanvar( bias, axis=1 ) )

    ret['g'] = ret['g2_lambda']/ret['g_lambda']
    ret['lambda'] = ret['g_lambda']/ret['g']
    return ret


def add_simulate_options(p):
    Simulation.add_folder_options(p)
    #print( p._actions )
    p.add_argument('--folder', type=str, default = 'test', help = 'name of the simulation folder' )
    #p.add_argument('--overwrite', action='store_true', default=argparse.SUPPRESS, help = 'overwrites folder' )
    p.add_argument('--append', action='store_true', default=argparse.SUPPRESS, help = 'appends to folder' )
    p.set_defaults( func=simulate )

def simulate( args ):
    print( args.folder )
    if os.path.exists( args.folder ):
        #if 'overwrite' in args:
            #print( 'delete folder "{}" already exists. Exiting'.format( args.folder ) )
        print( 'simulation folder "{}" already exists. Exiting'.format( args.folder ) )
        exit(0)
    else:
        os.mkdir( args.folder )

    table = Table()
    table.set_columns(['charge_gain', 'readout_noise', 'dark_current', 'dc*gain', 'type', 'bias_sigma', 'g_lambda', 'g', 'lambda' ])
    with Timer( 'simulate {} images'.format(args.number_of_images) ):

        args.output_fits = None
        args.sim = None
        count = 0
        while count < args.number_of_images:
            count += 1

            print( colored('count %s of %s'%(count, args.number_of_images), 'green') )
            args.charge_gain = np.random.uniform( *args.charge_gain_range )
            args.readout_noise = np.random.uniform( *args.readout_noise_range )
            args.dark_current = np.random.uniform( *args.dark_current_range )
            args.params_mode = 'median'
            print( 'generating values', args.charge_gain, args.readout_noise, args.dark_current )
            with Timer('simulate'):
                sim = Simulation.simulate_events( args )
                HDUlist = sim.generate_image( args )
                args.partlistHDU = HDUlist
            with Timer('get params'):
                params_table = analyse( args )
                for entry in params_table.list():
                    entry['charge_gain'] = args.charge_gain
                    entry['readout_noise'] = args.readout_noise
                    entry['dark_current'] = args.dark_current
                    entry['dc*gain'] = args.charge_gain*args.dark_current
                    table.append(entry)

        print(table)

def monitor( args ):
    import ConnieImage
    print( args.runID, args.ohdu )
    fmt = lambda x: {int:'{:5}', float:'{: .3e}', Section:'{: .3e}'}[type(x)]
    sfmt = lambda x: {int:'{:5.5}', float:'{:10.10}', Section:'{:10.10}'}[type(x)]

    first = True
    for _runID in args.runID:
        for _ohdu in args.ohdu:
            params = []
            params.append([ 'runID', _runID ])
            params.append(['ohdu', _ohdu ])
            params.append(['rn', float(ConnieImage.FullImage( runID=_runID, ohdu=_ohdu, imageType='raw' ).readoutNoiseEstimate()) ])
            dc = float( ConnieImage.FullImage( runID=_runID, ohdu=_ohdu, imageType='raw' ).darkCurrentEstimate() )
            params.append(['dc', dc ])
            params.append(['g_lambda', dc**2 ])
            if first:
                columns = zip(*params)[0]
                max_length = max( map(len, columns) )
                print( ' '.join( [ colored( sfmt(v).format(key), attrs=['bold']) for key, v in params ] ) )
                first = False
            print( ' '.join( [ fmt(v).format(v) for keys, v in params ] ) )

def add_HWFM_options( p ):
    p.set_defaults( func=test_HWFM )
    return

def test_HWFM( args ):
    print('not implemented yet')
    return

def mle_norm( x, axis=None ):
    fit = np.apply_along_axis( scipy.stats.norm.fit, axis, x )
    if axis == 1:
        mu, sigma = fit.T
    else:
        mu, sigma = fit
    return mu, sigma

def mle_poisson_norm( x, axis=None, mu=0, sigma=1, gain=1, fix_mu=False, fix_sigma=False, mode = 1 ):
    func = lambda x: stats.poisson_norm.fit(x, mu=mu, sigma=sigma, gain=gain, fix_mu=fix_mu, fix_sigma=fix_sigma, mode=mode)
    fit = np.apply_along_axis( func, axis, x )
    if axis == 1:
        gain, lamb = fit.T
    else:
        gain, lamb = fit
    return gain, lamb


def group_filenames( input_files, verbose=0 ):
    paths = sorted(glob( input_files ))
    list_of_paths = []
    for path in paths:
        match = re.search(r'(.*?)_p[0-9].fits', path )
        if match is None:
            entry = (path,)
        else:
            base_path = match.groups()[0]
            entry = tuple(sorted(glob( base_path+'_p*.fits*' )))
        if not entry in list_of_paths:
            list_of_paths.append( entry )
            if verbose:
                print( 'added group:' )
                for e in entry:
                    print(e)
    return list_of_paths

def apply_to_files( args ):
    list_of_paths = group_filenames( args.input_files )
    with Timer('apply to files'):
        returnDict = OrderedDict()
        for path_group in list_of_paths:
            partlistHDU = [ astropy.io.fits.open( path ) for path in path_group[::-1] ]
            parts_ohdu = OrderedDict()
            for j, listHDU in enumerate(partlistHDU):
                part_index = j+1
                for i, HDU in enumerate(listHDU):
                    try:
                        ohdu = HDU.header['OHDU']
                        if 'verbose' in args:
                            print( 'opening ohdu', ohdu )
                        if HDU.data is None: continue
                        if 'exclude' in args and args.exclude is not None and HDU.header['OHDU'] in args.exclude: continue
                        if 'ohdu' in args and args.ohdu is not None and HDU.header['OHDU'] not in args.ohdu: continue
                    except:
                        ohdu = i
                        if 'verbose' in args:
                            print( 'opening hdu', i )
                        if HDU.data is None: continue
                        if 'hdu' in args and args.hdu is not None and i not in args.hdu: continue
                    part = Part(HDU)
                    try: parts_ohdu[ohdu].append(part)
                    except KeyError: parts_ohdu[ohdu] = [part]
            returnDict[path_group] = OrderedDict()
            for i, parts in parts_ohdu.items():
                image = Image( parts )
                image.path = path_group
                try:
                    image.ohdu = image.header['OHDU']
                except:
                    image.ohdu = i
                returnDict[image.path][image.ohdu] = args.func( image, args )
                if 'verbose' in args:
                    print( 'applied function to ohdu', image.ohdu )
    return returnDict

def add_display_options( p ):
    p.add_argument('input_file', help = 'fits file input (example: "/share/storage2/connie/data/runs/*/runID_*_03326_*.fits.fz"' )
    p.add_argument(
        '--ohdu',
        type = int,
        default = 2,
        help = 'ohdu to be displayed',
        required = True
    )
    p.add_argument(
        '--plot',
        nargs = '+',
        type = str,
        default = ['image'],
        choices = ['proj', 'spectrum', 'image'],
        required = True,
        help = 'plot type'
    )
    p.add_argument(
        '--png',
        action = 'store_true',
        default = argparse.SUPPRESS,
        help = 'output to png file'
    )

    geom = p.add_argument_group('geometry options')
    geom.add_argument( '--E-range', nargs=2, type=eval, default = [-np.inf, np.inf], help = 'Emin Emax' )
    geom.add_argument( '--E-span', nargs=1, type=eval, default = argparse.SUPPRESS, help = 'mean+-E_span' )
    geom.add_argument( '--x-range', nargs=2, type=eval, default = [None, None], help = 'xmin xmax' )
    geom.add_argument( '--y-range', nargs=2, type=eval, default = [None, None], help = 'ymin ymax' )
    geom.add_argument( '--side', type=str, default = argparse.SUPPRESS, help = 'left or right amplifier' )
    geom.add_argument( '--section', type=str, default = argparse.SUPPRESS, help = 'data or bias' )
    geom.add_argument( '--trim', action='store_true', default = argparse.SUPPRESS, help = 'remove trim' )
    geom.add_argument( '--half', action='store_true', default = argparse.SUPPRESS, help = 'show half' )
    geom.add_argument( '--quarter', action='store_true', default = argparse.SUPPRESS, help = 'show quarter' )
    geom.add_argument( '--vmax', type=float, default = argparse.SUPPRESS, help = 'vmax' )

    corr = p.add_argument_group('correction options')
    corr.add_argument( '--remove', type=float, default=argparse.SUPPRESS, help = 'remove energies above' )
    corr.add_argument( '--global-bias', action='store_true', default=argparse.SUPPRESS, help = 'remove global bias' )
    corr.add_argument( '--correct-side', action='store_true', default=argparse.SUPPRESS, help = 'subtract right side' )
    corr.add_argument( '--smooth-lines', type=int, default=argparse.SUPPRESS, help = 'smooth lines' )
    corr.add_argument( '--sub-half', action='store_true', default=argparse.SUPPRESS, help = 'smooth lines' )
    corr.add_argument( '--overscan-subtraction', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'horizontal overscan subtraction' )
    corr.add_argument( '--vertical-overscan-subtraction', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'vertical overscan subtraction' )

    proj = p.add_argument_group('projection options')
    proj.add_argument( '--axis', type=int, default = 0, help = 'project on axis' )
    proj.add_argument( '--dev', action='store_true', default = argparse.SUPPRESS, help = 'supress the max line at the projection plot' )
    proj.add_argument( '--no-max', action='store_true', default = argparse.SUPPRESS, help = 'supress the max line at the projection plot' )
    proj.add_argument( '--no-min', action='store_true', default = argparse.SUPPRESS, help = 'supress the min line at the projection plot' )
    proj.add_argument( '--no-mean', action='store_true', default = argparse.SUPPRESS, help = 'supress the mean line at the projection plot' )
    proj.add_argument( '--no-center', action='store_true', default = argparse.SUPPRESS, help = 'supress the mean line at the projection plot' )
    proj.add_argument( '--no-median', action='store_true', default = argparse.SUPPRESS, help = 'supress the mean line at the projection plot' )
    proj.add_argument( '--smooth', type=int, default = argparse.SUPPRESS, help = 'smoothening length' )

    spec = p.add_argument_group('spectrum options')
    spec.add_argument( '--binsize', type=float, default=1, help = 'binsize' )

    p.set_defaults( func=display )

def display( args ):
    from matplotlib import pylab as plt

    with Timer('plot'):
        path = glob( args.input_file )
        if len(path) > 1:
            print( 'found {}'.format(path) )
            print( 'using the first path' )
        path = path[0]
        print( colored('path', 'green'), path )
        title = path.split('/')[-1]
        if not os.path.exists(path):
            print( 'file not found:', path )
            exit(0)
        listHDU = astropy.io.fits.open( path )
        imageHDU = None
        for i, HDU in enumerate(listHDU):
            if 'OHDU' in HDU.header:
                if HDU.header['OHDU'] == args.ohdu:
                    imageHDU = HDU
                    print( colored('ohdu', 'green'), args.ohdu )
                    break
            else:
                if i+1 == args.ohdu: imageHDU = HDU
        if imageHDU is None:
            print( 'ohdu {} was not found in {}'.format( args.ohdu, path ) )
            exit(0)

        data = imageHDU.data.astype(float)
        data = data[::-1,:]
        data[data>1e9] = np.nan

        height, width = data.shape
        print( colored('height', 'green'), height )
        print( colored('width', 'green'), width )
        number_of_amplifiers = width // constants.ccd_width
        print( colored('amplifiers', 'green'), number_of_amplifiers )
        if number_of_amplifiers > 0:
            side_width = width//number_of_amplifiers
        else:
            side_width = width
        print( colored('side width', 'green'), side_width )

        bias_width = int( np.ceil( (side_width - constants.ccd_width)/150. ) )*150
        print( colored('bias width', 'green'), bias_width )
        trim = constants.ccd_width + bias_width - side_width
        print( 'trim', trim )
        #bias_width = ( side_width + trim ) % constants.ccd_width

        height_trim = 1 #1

        #data.left = data[None:-height_trim, None:side_width]
        #data.right = data[None:-height_trim, side_width:None][:,::-1]

        if 'overscan_subtraction' in args:
            os_range = args.overscan_subtraction
            print( 'os_range', os_range, os_range[0] )
            correction = np.median( data[ None:None, os_range[0]:os_range[1] ], axis=1 )[:,None]
            data -= correction

        if 'vertical_overscan_subtraction' in args:
            os_range = args.vertical_overscan_subtraction
            print( 'os_range', os_range, os_range[0] )
            correction = np.median( data[ os_range[0]:os_range[1], None:None ], axis=0 )[None,:]
            data -= correction

        if 'side' in args:
            if args.side == 'left' or args.side == '0':
                print( colored('side', 'green'), 'left' )
                if 'correct_side' in args:
                    print( colored('subtract', 'green'), 'right' )
                    data = data[None:-height_trim, None:side_width] - data[None:-height_trim, side_width:None][:,::-1]
                else:
                    data = data[None:-height_trim, None:side_width]
            elif args.side == 'right' or args.side == '1':
                print( colored('side', 'green'), 'right' )
                data = data[None:-height_trim, side_width:None][:,::-1]
                if 'sub_half' in args:
                    print( colored('subtract line means', 'green'), 'true' )
                    correction = np.mean( data[None:None, constants.ccd_width/2:None], axis=1 )
                    if len(args.sub_half) > 0:
                        l = int(args.sub_half[0])
                        #correction =
                    data = data - correction[:,None]
            else:
                print( colored('side', 'green'), 'all' )

        if 'side' in args or number_of_amplifiers == 1:
            if 'half' in args:
                print( colored('half', 'green'), 'true' )
                data = data[None:None, -constants.ccd_width/2:None]
            elif 'quarter' in args:
                print( colored('quarter', 'green'), 'true' )
                data = data[None:None, -constants.ccd_width/4:None]

        if 'smooth_lines' in args:
            print( colored('smooth lines', 'green'), args.smooth_lines )
            data = scipy.ndimage.filters.uniform_filter1d( data, axis=1, size=args.smooth_lines, mode='mirror' )

        if 'trim' in args:
            print( colored('trim', 'green'), '8' )
            data = data[None:None, 8:]

        if 'global_bias' in args:
            print( colored('subtract', 'green'), 'mean(OS)' )
            data = data - np.nanmean( data[None:None, -bias_width:None] )

        if 'section' in args:
            if args.section == 'bias' or args.section == 'os':
                print( colored('section', 'green'), 'overscan' )
                data = data[None:None, -bias_width:None]
            if args.section == 'data' or args.section == 'ac':
                print( colored('section', 'green'), 'active' )
                data = data[None:None, None:-bias_width+5]
                if 'half' in args:
                    print( colored('half', 'green'), 'true' )
                    data = data[None:None, data.shape[1]/2:None]

        if 'remove' in args:
            print( colored('remove above', 'green'), args.remove )
            data[ data > args.remove ] = np.nan


        if 'x_range' in args:
            data = data[:, args.x_range[0]:args.x_range[1]]
            print( colored('x-range', 'green'), args.x_range )
        if 'y_range' in args:
            data = data[args.y_range[0]:args.y_range[1], :]
            print( colored('y-range', 'green'), args.y_range )
        if 'E_span' in args:
            E = np.nanmedian(data)
            mask = np.logical_or( data<E-args.E_span, data>E+args.E_span )
            data[mask] = np.nan
            print( colored('E-span', 'green'), args.E_span, np.nanmin(data), np.nanmax(data) )
        elif 'E_range' in args:
            data[data<=args.E_range[0]] = np.nan
            data[data>args.E_range[1]] = np.nan
            print( colored('E-range', 'green'), args.E_range, np.nanmin(data), np.nanmin(data) )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if args.plot[0] == 'proj':
            print( colored('plot', 'green'), 'projection' )
            axis = args.plot[1]
            if axis in ['lines', '1']:
                axis = 1
            elif axis in ['rows', '0']:
                axis = 0
            else:
                print( 'axis not expected' )
                exit()

            if not 'dev' in args:
                mins = np.nanmin( data, axis=axis )
                if not 'no_mean' in args:
                    means = np.nanmean( data, axis=axis )
                    print( means.shape )
                    ax.plot( means, '.', label='mean $\mu={:.4f}$'.format(np.nanmean(means)) )
                if not 'no_median' in args:
                    medians = np.nanmedian( data, axis=axis )
                    ax.plot( medians, '.', label='median $\mu={:.4f}$'.format(np.nanmean(medians)) )

                if 'smooth' in args:
                    #ax.plot( scipy.ndimage.filters.uniform_filter1d( centers, size=args.smooth, mode='nearest' ), '.', label='center{}'.format(args.smooth) )
                    ax.plot( scipy.ndimage.filters.uniform_filter1d( medians, size=args.smooth, mode='nearest' ), '.', label='median{}'.format(args.smooth) )
                    if not 'no_mean' in args:
                        ax.plot( scipy.ndimage.filters.uniform_filter1d( means, size=args.smooth, mode='nearest' ), '.', label='mean{}'.format(args.smooth) )
                if not 'no_max' in args:
                    maxs = np.nanmax( data, axis=axis )
                    ax.plot( maxs, '.', label='max' )
                if not 'no_min' in args:
                    mins = np.nanmin( data, axis=axis )
                    ax.plot( mins, '.', label='min' )
            else:
                mins = np.nanmin( data, axis=axis )
                if not 'no_mean' in args:
                    stds = np.nanstd( data, axis=axis )
                    ax.plot( stds, '.', label='std $\mu={:.4f}$'.format(np.nanmean(stds)) )
                mads = MAD( data, axis=axis )
                ax.plot( mads, '.', label='MAD $\mu={:.4f}$'.format(np.nanmean(mads)) )
                #fwhm = stats.FWHM( data, axis=axis, f=.001 )[0]
                #ax.plot( fwhm, '.', label='FWHM $\mu={:.4f}$'.format(np.nanmean(fwhm)) )
                sigmas = mle_norm( data, axis=axis )[1]
                ax.plot( sigmas, '.', label='$\sigma$[mle] $\mu={:.4f}$'.format(np.nanmean(sigmas)) )

        elif args.plot[0] == 'spectrum':
            bins = np.arange( np.nanmin(data), np.nanmax(data), args.binsize )
            y, x, dx = stats.make_histogram( data.flatten(), bins )
            ax.step( x, y, where='mid', label = '' )

        elif args.plot[0] == 'image' or args.plot[0] == 'matrix':
            cmap = matplotlib.cm.seismic
            cmap.set_bad(color='black')
            if 'log' in args.plot[1:]:
                im = np.log(data - np.nanmin(data) + 1)
                cmap.set_bad('black',0)
            else:
                im = data
                cmap.set_bad('black',0)
            try:
                vmean = np.nanmedian(im)
                vmax = np.nanmax(np.abs(im - vmean))
                if 'vmax' in args:
                    vmax = args.vmax

                print( 'vmax', vmax )
                if args.plot[0] == 'image':
                    func = ax.imshow
                elif args.plot[0] == 'matrix':
                    func = ax.matshow
                obj = func( im, cmap=cmap, origin='lower', vmin=vmean-vmax, vmax=vmean+vmax, interpolation='none' )
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.1)
                fig.colorbar( obj, cax )
                title += '\n mean {:.4}'.format( np.nanmean( im ) )
                title += '\n std {:.4}'.format( np.nanstd( im ) )
            except:
                print( 'im.shape', im.shape, data.shape )
                exit()
        else:
            print( 'plot option "{}" not implemented'.format( args.plot ) )
            exit(0)

        ax.set_title( title, fontsize=8 )
        ax.grid()
        ax.legend()
    if 'png' in args:
        outfile = args.input_file + '.' + '_'.join(args.plot) + '.png'
        print( 'outfile', outfile )
        fig.savefig( outfile, bbox_inches='tight', pad_inches=0 )
    else:
        plt.show()
    return

def add_analyse_options( p ):
    p.add_argument('name', help = 'fits file input (example: runID3326)' )
    p.add_argument('input_file', help = 'fits file input (example: "/share/storage2/connie/data/runs/*/runID_*_03326_*.fits.fz"' )

    p.add_argument( '--ohdu', nargs='*', type=int, default = None, help = 'ohdus to be analysed' )
    p.add_argument( '--exclude', nargs='*', type=int, default = None, help = 'ohdus not to be analysed' )

    p.add_argument('--params-mode', type=str, default = 'median', help = 'modes for parameter estimation' )
    p.add_argument( '--fix-vbias', action="store_true", default=argparse.SUPPRESS, help = 'mantain the vbais shift' )

    p.add_argument('--remove-hits', nargs=2, type=float, default=argparse.SUPPRESS, help = 'remove hits above ADU with border (example:60 3)' )

    p.add_argument('--find-hits', action="store_true", default=argparse.SUPPRESS, help = 'hits above ADU border' )

    p.add_argument( '--plot-part', action="store_true", default=argparse.SUPPRESS, help = 'plot parts' )
    p.add_argument( '--plot-sections', action="store_true", default=argparse.SUPPRESS, help = 'plot sections' )
    p.add_argument( '--plot-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot spectrum' )
    p.add_argument( '--plot-convolution-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot convolution spectrum' )
    p.add_argument( '--convolution-function', type=str, default=argparse.SUPPRESS, help = 'plot convolution function' )

    p.add_argument( '--plot-block-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot block spectrum' )
    p.add_argument( '--block-function', type=str, default=argparse.SUPPRESS, help = 'function applied to each block' )
    p.add_argument( '--block-function-mean', action="store_true", default=argparse.SUPPRESS, help = 'apply mean to each block' )
    p.add_argument( '--block-function-std', action="store_true", default=argparse.SUPPRESS, help = 'apply std to each block' )

    p.add_argument( '--plot-convolutionfft-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot convolution spectrum' )

    #p.add_argument('--extract', nargs=2, default=argparse.SUPPRESS, help = 'set to False to skip the extraction' )

    p.set_defaults( func=analyse )
    return

class Table:
    def __init__( self ):
        self.entries = []
    def list(self):
        return self.entries
    def append( self, entry ):
        self.entries.append( entry )
    def set_columns( self, cols ):
        self.cols = cols
    def get_entry( self, i ):
        return self.entries[i]
    def __str__( self ):
        ret = ', '.join(self.cols)
        for entry in self.entries:
            ret += '\n' + ', '.join( [
                ( {int:'{}', float:'{:.4f}', str:'{}', Section:'{:.4f}'}[type(entry[col])] ).format( entry[col] )
                for col in self.cols ] )
        return ret

def analyse( args ):
    if 'partlistHDU' in args:
        partlistHDU = [ args.partlistHDU ]
    else:
        paths = glob( args.input_file )
        partlistHDU = [ astropy.io.fits.open( path ) for path in paths[::-1] ]

    if not os.path.exists(args.name):
        os.mkdir(args.name)
    args.name = '{0}/'.format(args.name)

    table = Table()
    table.set_columns( ['ohdu', 'type', 'bias_sigma', 'g_lambda', 'g', 'lambda' ] )
    parts_dict = OrderedDict()
    for j, listHDU in enumerate(partlistHDU):
        part_index = j+1
        for i, HDU in enumerate(listHDU):
            if HDU.data is None: continue
            if 'exclude' in args and args.exclude is not None and HDU.header['OHDU'] in args.exclude: continue
            if 'ohdu' in args and args.ohdu is not None and HDU.header['OHDU'] not in args.ohdu: continue
            if 'plot_part' in args:
                title= 'o{}_p{}'.format(HDU.header['OHDU'], part_index)
                Section(HDU.data).save_pdf( '{0}{1}.pdf'.format(args.name, title) )
            part = Part(HDU)
            params = get_median_params_by_line( part.data.get_half(), part.bias )
            params.update( {'type':'median_by_lines', 'ohdu':'{}p{}'.format( HDU.header['OHDU'], part_index ) } )
            table.append( params )

            part.remove_outliers_from_bias()
            part.correct_lines(np.nanmean)
            try: parts_dict[i].append(part)
            except KeyError: parts_dict[i] = [part]

    first = True
    param_list = OrderedDict()
    for i, parts in parts_dict.items():
        image = Image( parts )
        ohdu = image.header['OHDU']
        args.name += '{}'.format( args.params_mode )
        params = OrderedDict()

        image.remove_outliers_from_vbias()
        params = get_mean_params_by_line( image.data.get_half(), image.bias )
        params.update( {'type': 'mean_by_line', 'ohdu': ohdu } )
        table.append( params )

        image.remove_outliers_from_vbias_second_pass()

        image.correct_cols( np.nanmean )
        #print( 'total charge', np.sum( image.data ), np.nanmean( image.data ) )

        params = get_mean_params_by_line( image.data.get_half().remove_hits( 60, border=3 ), image.bias )
        params.update( {'type': 'mean_by_line_remove60:3', 'ohdu': ohdu } )
        table.append( params )


        if 'remove_hits' in args:
            args.name += '_e{0}b{border}'.format( args.remove_hits[0], border=args.remove_hits[1] )
            image.data = image.data.remove_hits( args.remove_hits[0], border=args.remove_hits[1])

        if 'find_hits' in args:
            bias_std = np.nanstd(image.bias)
            print_var( 'bias_std', locals() )
            outliers2nanPoissonNorm_1d( image.data, sigma = np.nanstd(image.bias), pmin = args.find_hits )

        if 'plot_sections' in args:
            image.save_sections( r'{}_o{}'.format( args.name, image.header['OHDU'] ) )
            image.bias.save_projection( r'{}_o{}_biasProj'.format( args.name, image.header['OHDU'] ), title='bias', axis=1, no_max=True, no_mean = True, no_min = True, do_regress=True )
            image.vbias.save_projection( r'{}_o{}_vbiasProj'.format( args.name, image.header['OHDU'] ), title='vbias', axis=0, no_max=True, no_mean = True, no_min = True, do_regress=True )
            image.data.save_projection( r'{}_o{}_dataProj1'.format( args.name, image.header['OHDU'] ), title='data', axis=1, no_max=True, no_mean = True, no_min = True, do_regress=True )
            image.data.save_projection( r'{}_o{}_dataProj0'.format( args.name, image.header['OHDU'] ), title='data', axis=0, no_max=True, no_mean = True, no_min = True, do_regress=True )

        if 'plot_spectrum' in args:
            for sec in ['dbias', 'vbias', 'bias', 'data']:
                title = r'o{}_{}_spectrum'.format( image.header['OHDU'], sec )
                getattr(image, sec).save_spectrum( '{}_{}'.format(args.name, title ), title=title )

        if 'plot_convolution_spectrum' in args:
            func = np.nansum
            if 'convolution_function' in args:
                func = eval( 'lambda x:' + args.convolution_function )
            for sec in ['dbias', 'vbias', 'bias', 'data']:
                for length in [5,10]:
                    title = r'o{ohdu}_{sec}_convolution{length}'.format( ohdu=image.header['OHDU'], sec=sec, length=length )
                    getattr( image, sec ).convolve( function=func, size=(length,length) )\
                        .save_spectrum( '{}_{}'.format(args.name, title ), title=title )

        if 'plot_block_spectrum' in args:
            func = np.nansum
            if 'block_function' in args:
                func = eval( 'lambda x:' + args.block_function )
            if 'block_function_mean' in args:
                func = np.nanmean
            if 'block_function_std' in args:
                func = lambda x: (lambda y: np.nan if y==0 else y)(np.nanstd(x))

            for sec in ['dbias', 'vbias', 'bias', 'data']:
                for length in [5,10]:
                    title = r'o{ohdu}_{sec}_block{length}'.format( ohdu=image.header['OHDU'], sec=sec, length=length )
                    section = getattr( image, sec )
                    h, w = section.shape
                    section = section[ :-(h%length), :( -(w%length) if (w%length)>0 else None) ]
                    blocks = blockshaped( section, length, length )
                    ret = [ func(block) for block in blocks ]
                    Section(ret).save_spectrum( '{}_{}'.format(args.name, title ), title=title, binsize=.01 )

        if 'plot_dark_current_evolution' in args:
            func = np.nanmedian
            for sec in ['dbias', 'vbias', 'bias', 'data']:
                for length in [5,10]:
                    title = r'o{ohdu}_{sec}_block{length}'.format( ohdu=image.header['OHDU'], sec=sec, length=length )
                    section = getattr( image, sec )
                    medians = func( section, axis = 1 )
                    Section(ret).save_spectrum( '{}_{}'.format(args.name, title ), title=title, binsize=.01 )

        #fmt = lambda x: {int:'{:4}', float:'{: .3e}', Section:'{: .3e}'}[type(x)]
        #sfmt = lambda x: {int:'{:4.4}', float:'{:10.10}', Section:'{:10.10}'}[type(x)]

        #fname = '{}_params.csv'.format( args.name )
        #open( fname, 'w' )
        #if first:
            #columns = zip(*params)[0]
            #max_length = max( map(len, columns) )
            #print( ' '.join( [text('ohdu', mode='B')] + [ text( sfmt(v).format(key), mode='B') for key, v in params ] ) )

            #open( fname, 'a' ).write( '# ' + ', '.join( ['ohdu'] + [ key for key, v in params ] ) + '\n' )
            #first = False
        #column_head = '%4d' % image.header['OHDU']
        #print( ' '.join( [text( column_head, mode='B' )] + [ fmt(v).format(v) for keys, v in params ] ) )
        #open( fname, 'a' ).write( ', '.join( [column_head] + [ fmt(v).format(v) for keys, v in params ] ) + '\n' )
    print(table)
    return table


def add_monitor_options( p, func ):
    p.add_argument( '--runID', nargs='*', type=int, default = None, help = 'runIDs to be analysed' )
    p.add_argument( '--ohdu', nargs='*', type=int, default = [2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14], help = 'ohdus to be analysed' )
    p.add_argument( '--exclude', nargs='*', type=int, default = None, help = 'ohdus not to be analysed' )
    p.set_defaults( func=func )
    return

def superpose( args ):
    print( 'inputs', args.input_files )
    print( 'output', args.output )
    if os.path.exists(args.output):
        print( 'output already exists' )
        exit(0)
    input_files = args.input_files
    subargs = args
    subargs.func = lambda x, y: x
    images = {}
    for input_file in input_files:
        subargs.input_files = input_file
        print( 'input_file', subargs.input_files )
        result = apply_to_files( subargs )
        images[result.keys()[0]] = result.values()[0]

    new_ohdus = None
    for path, ohdus in images.items():
        print( path, ohdus )
        if new_ohdus is None:
            new_ohdus = ohdus
            continue
        for key in new_ohdus.keys():
            new_ohdus[key].all += ohdus[key].all
            print( new_ohdus[key].all, new_ohdus[key].all.shape, type(new_ohdus[key].all) )

    primary = None
    hdu_list = []
    for key, image in new_ohdus.items():
        print( key, image.header )
        header = astropy.io.fits.Header()
        for field, value in image.header.items():
            header[field] = value

        hdu_list += [astropy.io.fits.ImageHDU( image.all, header=header )]

        if primary is None:
            primary_header = header
            del primary_header['OHDU']
            primary = astropy.io.fits.PrimaryHDU( header = header )

    #print( [primary] + hdu_list )
    HDU = astropy.io.fits.HDUList([primary] + hdu_list )
    HDU.writeto( args.output )
    print( 'saved', args.output )
    return

def add_superpose_options( p ):
    p.add_argument( 'input_files', nargs='+' )
    p.add_argument( '--indices', nargs='*', default = None, type=int, help='indexes to be shown' )
    p.add_argument( '--include', nargs='*', default = '', help='fields to be shown' )
    p.add_argument( '--exclude', nargs='*', default = '', help='fields not to be shown' )
    p.add_argument( '--output', default = 'output', help='output name' )
    p.set_defaults( func=superpose )
    return

def add_header_options( p ):
    p.add_argument( 'input_file' )
    p.add_argument( '--indices', nargs='*', default = None, type=int, help='indexes to be shown' )
    p.add_argument( '--include', nargs='*', default = '', help='fields to be shown' )
    p.add_argument( '--exclude', nargs='*', default = '', help='fields not to be shown' )
    p.set_defaults( func=read_header )
    return

def read_header( args ):
    for file in glob( args.input_file ):
        listHDU = astropy.io.fits.open( file )
        for i in range(len(listHDU)) if args.indices is None else args.indices:
            print( colored( 'HDU %s' % i, 'green', attrs=['bold'] ) )
            for key, value in listHDU[i].header.items():
                if key in args.exclude: continue
                if key in args.include or len(args.include) is 0:
                    print( colored( '%s' % key, 'green' ), '%s' % value, end=' ' )
            print()
    return

def postprocess( args ):
    if 'image_mode' in args:
        if args.image_mode is '1':
            args.rebin = [1,1]
            args.horizontal_overscan = 150
            args.vertical_overscan = 90
            args.xyshape = [4150,4120]
        elif args.image_mode is '5':
            args.rebin = [5,1]
            args.horizontal_overscan = 450
            args.vertical_overscan = 70
            args.xyshape = [4150,4120]
    return args

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description = 'image tools',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        )
    subparsers = parser.add_subparsers( help = 'major options' )

    add_simulate_options( subparsers.add_parser('simulate', help='simulate and analyse random') )
    add_analyse_options( subparsers.add_parser('analyse', help='analyse image') )
    add_header_options( subparsers.add_parser('header', help='read image header') )
    add_monitor_options( subparsers.add_parser('monitor', help='monitorViewer functions'), monitor )
    add_HWFM_options( subparsers.add_parser('hwfm', help='test the hwfm algorithm') )
    add_display_options( subparsers.add_parser('display', help='display the image') )
    add_superpose_options( subparsers.add_parser('superpose', help='superpose images') )

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(1)

    args = postprocess( args )
    print_var( vars(args).keys(), vars(args), line_char='\t' )
    with Timer('finished'):
        args.func( args )
