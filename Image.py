# coding: utf-8

from __future__ import print_function
import numpy as np
import astropy.io.fits
import scipy.stats
import scipy.ndimage
from glob import glob
import re
import gi
#gi.require_version('Gtk', '3.0')
#from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib

import matplotlib
matplotlib.use('gtk3agg')
#from matplotlib.backends.backend_gtk3cairo import (FigureCanvasGTK3Cairo as FigureCanvas)

#from SimulateImage import Hits
import Statistics as stats
from Timer import Timer

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

def crop_mean( x, axis = None, mode = None ):
    if mode == None:
        return np.apply_along_axis( crop_mean_1d, axis, x )
    elif mode == 'partition':
        return crop_mean_1d_partition( x )

#def crop_var_1d( x, step = 1 ):
    #x2 = crop_mean_1d(x**2, step)
    #m = crop_mean_1d(x, step)**2
    #return x2 - m**2    

def crop_mean_1d( x, step = 1 ):
    x = np.sort( x.flatten() )
    left, right = 0, len(x)

    mean = np.nanmean(x)
    median = np.nanmedian(x)

    while abs(mean - median) >= .5:
        if mean > median: right -= step
        else: left += step

        mean = np.nanmean( x[left:right] )
        median = np.nanmedian( x[left:right] )
    return mean

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

def parse_image( imagesHDU, correct_lines_by_bias_median = np.nanmedian, smoothening_length = 0, correct_rows_by_bias_median = None, row_correction_function = np.nanmedian ):
    with Timer('read and process') as t:
        header = None
        data = None
        for imageHDU in imagesHDU:
            datum, header = parse_single_image( imageHDU, correct_lines_by_bias_median, smoothening_length )
            if data:
                data = { key: np.concatenate( [ data[key], datum[key] ], axis=0 ) for key in data.keys() }
            else:
                data = datum
        data['bias'] = data['bias'].view(Image)
        #print( 'header', header )
        if 'BIASH' in header:
            vbias_width = header['BIASH']
        else:
            vbias_width = vbias_from_height[ data['data'].shape[0] ]
        data['vbias'] = data['data'][:vbias_width,:].view(Image)
        data['data'] = data['data'][vbias_width:,:].view(Image)
        data['bias'] = data['bias'][vbias_width:,:].view(Image)
        if correct_rows_by_bias_median:
            data['data'], data['vbias'] = correct_rows_by_bias( [ data['data'], data['vbias'] ], data['vbias'], func = correct_rows_by_bias_median, size = smoothening_length )
    return data, dict(header)

def read_image_from_fits( path, ohdu ):
    paths = glob( expr )[::-1]
    imagesHDU = []
    for path in paths:
        fits_dir = astropy.io.fits.open( path )
        imagesHDU.append( get_imageHDU( fits_dir, ohdu ) )
    return imagesHDU
    
def parse_single_image( imageHDU, correct_lines_by_bias_median = np.nanmedian, smoothening_length = 0 ):
    section = {}
    if imageHDU.header['ohdu'] == -1:
        half_width = imageHDU.data.shape[1]
        bias_width = imageHDU.header['biasW']
        bin_width = imageHDU.header['rebin']
    else:
        half_width = imageHDU.data.shape[1]/2
        bias_width = bias_from_width[ imageHDU.data.shape[1] ]
        bin_width = bin_from_width[ imageHDU.data.shape[1] ]

    section['data'] = Section( imageHDU.data[ :-1, :half_width - bias_width ] )#/bin_width
    section['bias'] = Section( imageHDU.data[ :-1, half_width -bias_width: half_width ] )#/bin_width

    if correct_lines_by_bias_median is None:
        section['data'], section['bias'] = correct_global_by_bias( [section['data'], section['bias']], section['bias'] )
    else:
        section['data'], section['bias'] = correct_lines_by_bias( [section['data'], section['bias']], section['bias'], func = correct_lines_by_bias_median, size = smoothening_length )
    return section, dict(imageHDU.header)

def correct_global_by_bias( data, bias ):
    correction = np.nanmedian( bias )
    return [ datum - correction for datum in data ]

def correct_lines_by_bias( data, bias, size = 0, func = np.nanmedian ):
    correction_of_lines = func( bias, axis = 1 )
    if size > 2:
        correction_of_lines = scipy.ndimage.filters.uniform_filter1d( correction_of_lines, size = size, mode='nearest' )
    return [ datum - correction_of_lines[:,None] for datum in data ]

def correct_rows_by_bias( data, bias, size = 0, func = np.nanmedian ):
    correction_of_rows = func( bias, axis = 0 )
    #print( 'correction_of_rows', correction_of_rows.shape )
    if size > 2: 
        correction_of_rows = scipy.ndimage.filters.uniform_filter1d( correction_of_rows, size = size, mode='nearest' )
    return [ datum - correction_of_rows[None,:] for datum in data ]

def label_clusters( condition_image ):
    s33 = [[1,1,1], [1,1,1], [1,1,1]]
    return scipy.ndimage.label( condition_image, structure=s33 )[0]

class Part:
    def __init__( self, imageHDU ):
        self.header = dict( imageHDU.header )
        if imageHDU.header['ohdu'] == -1:
            half_width = imageHDU.data.shape[1]
            bias_width = imageHDU.header['biasW']
            bin_width = imageHDU.header['rebin']
            self.vbias_height = imageHDU.header['biasH']
        else:
            half_width = imageHDU.data.shape[1]/2
            bias_width = Part.get_bias_width_from_width( imageHDU.data.shape[1] )
            bin_width = Part.get_rebin_from_height( imageHDU.data.shape[0] )
            self.vbias_height = Part.get_vbias_height_from_height( imageHDU.data.shape[0] )

        self.data = Section( imageHDU.data[ :-1, :half_width - bias_width ] )#/bin_width
        self.bias = Section( imageHDU.data[ :-1, half_width -bias_width: half_width ] )#/bin_width

    def correct_lines( self, func = np.nanmedian, smoothening_length = 0 ):
        if func is None:
            self.data -= self.bias.get_global_bias( func ) 
            self.bias -= self.bias.get_global_bias( func )
        else:
            lines_correction = self.bias.get_lines_bias( func, smoothening_length )
            self.data -= lines_correction
            self.bias -= lines_correction
        return self
        
    @staticmethod
    def get_rebin_from_height( height ):
        return { 900: 5, 1055: 1 }[height]
    
    @staticmethod
    def get_bias_width_from_width( width ):
        return { 9340: 450, 8540: 150 }[width]

    @staticmethod
    def get_vbias_height_from_height( height ):
        return { 900: 70, 1055: 90 }[height]

class Image:
    def __init__( self, parts ):
        self.data = None
        self.bias = None
        self.header = parts[0].header
        self.vbias_height = parts[0].vbias_height
        for part in parts:
            if self.data is None:
                self.data = part.data
                self.bias = part.bias
            else:
                self.data = np.concatenate( [ self.data, part.data ], axis=0 )
                self.bias = np.concatenate( [ self.bias, part.bias ], axis=0 )

        self.vbias = self.data[:self.vbias_height,:].view(Section)
        self.data = self.data[self.vbias_height:,:].view(Section)
        self.dbias = self.bias[:self.vbias_height,:].view(Section)
        self.bias = self.bias[self.vbias_height:,:].view(Section)
    
    def correct_rows( self, func = np.nanmedian, smoothening_length = 0 ):
        if not func is None:
            data_rows_correction = self.vbias.get_rows_bias( func, smoothening_length )
            self.data -= data_rows_correction
            self.vbias -= data_rows_correction
            bias_rows_correction = self.dbias.get_rows_bias( func, smoothening_length )
            self.bias -= bias_rows_correction
            self.dbias -= bias_rows_correction
        return self
    
    def get_params( self, mode = None, half = False, gain = 1, factor = 1, remove_hits = False, **kwargs ):
        if 'GAIN' in self.header:
            gain = self.header['GAIN']

            #if remove_hits:
                #data['data'] = data['data'].get_background( threshold, border )
        if half: data = self.data.get_half()
        else: data = self.data
        if remove_hits: data = data.remove_hits(60, 3)
        #data /= factor
        
        if mode == 'mean&std':
            mu = np.nanmean( self.bias )
            stds = np.nanstd( self.bias, axis=1 )
            sigma = np.nanmean( stds )
            sigma_err = np.nanstd( stds )/np.sqrt(len(stds))
            
            diff_means = np.nanmean( data, axis=1 ) - np.nanmean( self.bias, axis=1 )
            g_lamb = np.nanmean( diff_means )
            g_lamb_err = np.nanstd( diff_means )/np.sqrt( len( diff_means ) )
            
            diff_vars = np.nanvar( data, axis=1 ) - np.nanvar( self.bias, axis=1 )
            g2_lamb = np.nanmean( diff_vars )
            g2_lamb_err = np.nanstd( diff_vars )/np.sqrt( len( diff_vars ) )
            
            return {'mu': mu, 
                    'sigma': sigma, 
                    'sigma_err':sigma_err, 
                    'g': g2_lamb/g_lamb,
                    'lambda': g_lamb/gain/factor,
                    'lambda_err': g_lamb_err/gain/factor,
                    }

        if mode == 'median&MAD':
            mu = np.nanmedian( np.nanmedian( self.bias, axis=1 ) )
            mads = MAD( self.bias, axis=1)
            sigma = np.nanmean( mads )
            sigma_err = np.nanstd( mads )/np.sqrt(len(mads))
            
            diff_medians = np.nanmedian( data, axis=1) - np.nanmedian( self.bias, axis=1)
            g_lamb = np.nanmean( diff_medians )
            g_lamb_err = np.nanstd( diff_medians )/np.sqrt( len( diff_medians ) )
            
            diff_mads = MAD( data, axis=1) - MAD( self.bias, axis=1)
            g2_lamb = np.nanmean( diff_mads )
            g2_lamb_err = np.nanstd( diff_mads )/np.sqrt( len( diff_mads ) )
            return {'mu':mu, 
                    'sigma':sigma, 
                    'sigma_err':sigma_err, 
                    'g': g2_lamb/g_lamb,
                    'lambda': g_lamb/gain/factor,
                    'lambda_err': g_lamb_err/gain/factor,
                    }
        
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
    
    def save_pdf( self, output ):
        data = np.concatenate([self.vbias, self.data], axis=0)
        bias = np.concatenate([self.dbias, self.bias], axis=0)
        image = Section( np.concatenate([data, bias], axis=1) )
        image.save_pdf( output )
        
    
class Section( np.ndarray ):
    def __new__( cls, input_array ):
        if type(input_array) is float: return input_array
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
    
    def get_rows_bias( self, func = np.nanmedian, smoothening_length = 0 ):
        rows_bias = func( self, axis = 0 )
        if smoothening_length > 2:
            rows_bias = scipy.ndimage.filters.uniform_filter1d( rows_bias, size = smoothening_length, mode='nearest' )
        return rows_bias[None,:]
            
    def get_binned_distribution( self, binsize = 1):
        bins = np.arange( self.min(), self.max(), binsize)
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
        
    def extract_hits( self, mode, **kwargs ):
        if mode == 'cluster':
            return self._extract_clusters_( **kwargs )

    def get_clusters( self, threshold, border ):
        labeled_clusters = label_clusters( self >= threshold )
        print( 'number of clusters above threshold', labeled_clusters.max() )
        is_cluster = labeled_clusters > 0
        distances_to_cluster = scipy.ndimage.distance_transform_edt( is_cluster == False )
        labeled_clusters = label_clusters( distances_to_cluster <= border )
        print( 'number of clusters with border', labeled_clusters.max() )
        return labeled_clusters, distances_to_cluster

    def get_background( self, threshold, border ):
        labeled_clusters, _ = self.get_clusters( threshold, border )
        return self[ labeled_clusters == 0 ]

    def remove_hits( self, threshold, border ):
        labeled_clusters, _ = self.get_clusters( threshold, border )
        data = Section( np.where( labeled_clusters == 0, self, np.nan ) )
        print('hits removed', len(self.flatten()), len(self[labeled_clusters == 0]) )
        return data
        
    def _extract_clusters_( self, threshold, border ):
        labeled_clusters, distances_to_cluster = self.get_clusters( threshold, border )
        
        list_of_clusters = scipy.ndimage.labeled_comprehension(
            self,
            labeled_clusters,
            index = np.unique(labeled_clusters), 
            func = lambda v, p: [v, p], 
            out_dtype=list, 
            default=-1,
            pass_positions=True 
            )
        levels = scipy.ndimage.labeled_comprehension( 
            distances_to_cluster, 
            labeled_clusters, 
            index= np.unique(labeled_clusters),
            func = lambda v: v, 
            default=-1, 
            out_dtype=list 
            )
        def process( cluster ):
            ei, level = cluster
            x, y = np.unravel_index(ei[1], self.shape)
            return ei[0], x, y, level
        list_of_clusters = map( process, zip(list_of_clusters, levels) )
        
        return list_of_clusters, levels, threshold, border

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

    def save_pdf( self, output ):
        output += '.pdf'
        with Timer( 'saved ' + output ):
            from matplotlib import pylab as plt

            fig = plt.figure()
            ax = fig.add_subplot(111)

            im = np.log(self - self.min() + 1)
            ax.imshow( im, cmap='Blues', origin='lower', vmin=im.min(), vmax=im.max() )
            fig.savefig( output )
        return
        
    def display( self, mode, delta = None ):
        from matplotlib import pylab as plt
        from matplotlib.colors import LogNorm
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = -self
        if mode == 'xy':
            vmin = None
            vmax = None
            if delta:
                vmin = np.nanmedian(datum)-delta
                vmax = np.nanmedian(datum)+delta
            ax.imshow( im, cmap='Blues', vmin=vmin, vmax=vmax)
        if mode == 'x':
            ax.plot(self[0,:], '.')
            ax.plot(self[-1,:], '.')
            ax.set_ylim( (np.nanmedian(self)-delta,np.nanmedian(self)+delta) )
        if mode == 'y':
            ax.plot(self[:, 0], '.')
            ax.plot(self[:, -1], '.')
            ax.set_ylim( (np.nanmedian(self)-delta,np.nanmedian(self)+delta) )
        if mode == 'flat':
            ax.plot(self, '.')
        plt.show()

def display( data, mode, delta = None, labels = None, nbins = None, log = False ):
    from matplotlib import pylab as plt
    from matplotlib.colors import LogNorm
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if labels is None:
        labels = map( str, range(len(data)) )
    for i, datum in enumerate(data):
        if mode == 'flat':
            ax.plot( datum, '.', label = labels[i] )
        if mode == 'medianx':
            ax.plot( np.nanmedian(datum, axis=0), '.', label = labels[i] )
        if mode == 'mediany':
            ax.plot( np.nanmedian(datum, axis=1), '.', label = labels[i] )
        if mode == 'xy':
            vmin = datum.min()
            vmax = datum.max()
            if delta:
                vmin = np.nanmedian(datum)-delta
                vmax = np.nanmedian(datum)+delta
            ax.imshow( datum, cmap='Blues', vmin = vmin, vmax = vmax )
            ax.set_title(labels[i])
        if mode == 'E':
            bins = None
            if nbins == None:
                nbins = int(np.sqrt(len(datum)))
                bins = np.linspace( datum.min(), datum.max(), nbins )
            else:
                bins = np.linspace( datum.min(), datum.max(), nbins )
            if delta:
                bins = np.linspace( np.nanmedian(datum)-delta, np.nanmedian(datum)+delta, nbins )
            ax.hist( datum.flatten(), bins = bins, histtype = 'step', label = labels[i] )
            if log:
                ax.set_yscale('log')
    legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper right' )
    plt.show()

def simulate_and_test( correct_lines_by_bias_median, correct_rows_by_bias_median, smoothening_length, mean_function, std_function, median_diff_function, var_diff_function, output_file ):
    pass

def simulate( args ):
    from TerminalColor import text 
    import Simulation
    if 'func' in args: del args.func
    for key, value in vars(args).items():
        print( '\t%s: %s' % ( text(key, mode='B'), value ) ) 
    with Timer( 'simulate %s images' % args.number_of_images ):
        count = 1
        columns = ['n', 'readout_noise', 'dark_current', 'charge_gain', 'mu', 'sigma', 'sigma_err', 'lambda', 'lambda_err', 'g']
        splits = args.output_table.split('.')
        args.output_table = '.'.join( splits[:-1] + [args.params_mode, splits[-1]] )
        open( args.output_table, 'w' ).writelines( '#' + ', '.join(columns) + '\n' )
        args.output_fits = None
        args.sim = None
        while count <= args.number_of_images:
            print( text('count %s of %s'%(count, args.number_of_images), color='green') )
            args.charge_gain = np.random.uniform( *args.charge_gain_range )
            args.readout_noise = np.random.uniform( *args.readout_noise_range )
            args.dark_current = np.random.uniform( *args.dark_current_range )
            print( 'values', args.charge_gain, args.readout_noise, args.dark_current )
            with Timer('simulate'):
                sim = Simulation.simulate_events( args )
                imageHDU = sim.generate_image()
            with Timer('get params'):
                part = Part( imageHDU ).correct_lines()
                image = Image( [part] ).correct_rows()
                params = image.get_params( mode=args.params_mode, factor=args.rebin[0], remove_hits=args.remove_hits )
            if count == 1: image.save_pdf( args.output_table + '.pdf' )
            print( 'params', params['g'], params['sigma'], params['lambda'] )
            for param in params:
                if param is None: continue
            table = vars(args)
            #table['charge_gain*dark_current'] = args.charge_gain*args.dark_current
            #table['charge_gain**2*dark_current'] = args.charge_gain**2*args.dark_current
            table['n'] = count
            entry = [ table[column] if column in table else params[column] for column in columns ]
            open( args.output_table, 'a' ).writelines( ', '.join(map( str, entry)) + '\n' )
            count += 1
    import Plot
    Plot.generate_plot( [args.output_table], xcolumn = 'readout_noise', ycolumns = ['sigma'], yerrcolumns = ['sigma_err'], labels = ['sigma'], output_file = args.output_table + '.readout_noise.pdf', figscale = [1,1], diagonal_line = True, legend_loc = 'lower right' )
    Plot.generate_plot( [args.output_table], xcolumn = 'dark_current', ycolumns = ['lambda'], yerrcolumns = ['lambda_err'], labels = ['lambda'], output_file = args.output_table + '.dark_current.pdf', figscale = [1,1], diagonal_line = True, legend_loc = 'lower right' )
        
        
def analyse( **kargs ):
    pass

def main( args ):

    def median_by_lines( x ):
        medians = np.nanmedian(x, axis=1)
        return np.nanmean( medians )

    def MAD_by_lines2( x ):
        correction = np.nanmedian( x, axis=0 )
        mads_by_line = MAD( x - correction[None,:], axis = 1 )
        return float( np.nanmean( mads_by_line ) )

    def MAD_by_lines( x ):
        return float( np.mean( MAD( x, axis=1 ) ) )
    
    def median_diff_by_half_lines( x, y, g ):
        return float(np.mean( np.median(x[:,x.shape[1]/2:], axis=1) - np.median(y, axis=1) ))

    def MADsqr_diff_by_half_lines( x, y, g ):
        return float(np.mean( MAD(x[:,x.shape[1]/2:], axis=1)**2 - MAD(y, axis=1)**2 ))

    def MAD_by_lines_corrected_rows( x ):
        return float( np.mean( MAD(correct_rows_by_bias( [x], x, func = np.median )[0], axis=1) ) )
    
    def fit_sigma( data ):
        data = data.flatten()
        bins = np.arange(data.min(), data.max(), 1)
        hist, _ = np.histogram(data, bins)
        fit_func = lambda x,A,mu,sig: A*scipy.stats.norm.pdf(x, mu, sig)
        p0 = [np.sum(data),0., 1.]
        try:
            params = scipy.optimize.curve_fit( fit_func, .5*(bins[:-1] + bins[1:]), hist, p0=p0 )[0]
        except TypeError:
            print( 'TypeError' )
            return None
        return float( params[2] )

    def fit_gain_lambda( data, bias, gain ):
        sigma = fit_sigma( bias )
        if sigma is None: return None
        data = data.flatten()
        bins = np.arange(data.min(), data.max(), 1)
        x = .5*(bins[:-1] + bins[1:])
        hist, _ = np.histogram( data, bins )
        
        A = float( np.sum(hist) )
        
        fit_func = lambda x,A,loc,mu: A*stats.poisson_norm.pdf(x, loc, sigma, gain, mu )
        p0 = [A, 0, .001]
        try:
            params = scipy.optimize.curve_fit( fit_func, x, hist, p0=p0, bounds=[(0,-np.inf,0), (np.inf, np.inf, np.inf)] )[0]
        except TypeError:
            print( 'TypeError' )
            return None
        except ValueError:
            print( 'ValueError' )
            return None
        return float( params[2]*gain )
        
    def mle_sigma( data ):
        data = data.flatten()
        fit_func = lambda x, mu, sigma: scipy.stats.norm.pdf(x, mu, sigma)
        negloglikelyhodd = lambda mu, sigma: -np.sum( np.log( fit_func( data, mu, sigma ) ) )
        x0 = [0., 1.]
        try:
            params = scipy.optimize.minimize( 
                fit_func, 
                x0=x0,
                bounds=[(-np.inf, np.inf), (1e-3, np.inf)]
                ).x
        except TypeError:
            print( 'TypeError' )
            return None
        return float( params[1] )
    
    def mle_gain_lambda( data, bias, gain ):
        sigma = mle_sigma( bias )
        if sigma is None: return None
        data = data.flatten()
        
        fit_func = lambda x, loc, lamb: stats.poisson_norm.pdf(x, loc, sigma, gain, lamb )
        negloglikelyhodd = lambda loc, lamb: -np.sum( np.log( fit_func( data, loc, lamb ) ) )
        p0 = [0, .001]
        try:
            params = scipy.optimize.minimize( 
                fit_func, 
                x0=x0, 
                bounds=[(-np.inf,np.inf), (1e-3, np.inf)] 
                ).x
        except TypeError:
            print( 'TypeError' )
            return None
        except ValueError:
            print( 'ValueError' )
            return None
        return float( params[1]*gain )
    
    simulate_and_test( 
        output_file = 'simulate_and_test_pres_dc_mle.csv',
        correct_lines_by_bias_median = np.median, 
        correct_rows_by_bias_median = None,
        smoothening_length = 0,
        mean_function = median_by_lines,
        std_function = mle_sigma,
        #median_diff_function = median_diff_by_half_lines,
        median_diff_function = mle_gain_lambda,
        var_diff_function = MADsqr_diff_by_half_lines,
        )
    
    exit(0)
 
if __name__ == '__main__':
    import argparse
    import sys
    import Simulation

    def tuple_of( type_ ):
        return lambda x: map( type_, eval(x.replace('\"', '')) )

    parser = argparse.ArgumentParser(
        description = 'image tools',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        )
    subparsers = parser.add_subparsers( help = 'major options' )
    parser_simulate = subparsers.add_parser('simulate', help='simulate and analyse random')
    parser_simulate.add_argument('number_of_images', type=int, help = 'charge gain range' )
    parser_simulate.add_argument('--params-mode', type=str, default = 'median&MAD', help = 'modes for parameter estimation' )
    parser_simulate.add_argument('--output-table', type=str, default = 'analysis.csv', help = 'csv output file' )
    parser_simulate.add_argument('--charge-gain-range', type=tuple_of(float), default = '\"[5,9]\"', help = 'charge gain range' )
    parser_simulate.add_argument('--readout-noise-range', type=tuple_of(float), default = '\"[11,14]\"', help = 'readout noise range' )
    parser_simulate.add_argument('--dark-current-range', type=tuple_of(float), default = '\"[0.01,0.3]\"', help = 'dark current range' )
    parser_simulate.add_argument('--rebin', type=tuple_of(int), default = '\"[1,1]\"', help = 'rebin' )
    parser_simulate.add_argument('--image-mode', type=str, default = 'none', help = 'set to "1" to use official 1x1 image geomtry or "5" to 1x5' )
    parser_simulate.add_argument('--remove-hits', action='store_true', help = 'remove hits above 60ADU 3border' )
    parser_simulate.add_argument('--number-of-charges', type=int, default = '4000', help = 'number of charges to be randomly generated' )
    parser_simulate.add_argument('--charge-range', type=tuple_of(int), default = '\"[5,200]\"', help = 'range into which to randomly generate charges' )
    parser_simulate.add_argument('--number-of-Cu-charges',
                        type=int,
                        default = '0',
                        help = 'number of charges to be randomly generated at the Copper fluorescence energy 8.046keV' 
                        )
    parser_simulate.add_argument('--number-of-Cu2-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the secundary Copper fluorescence energy 8.904keV' 
                        )
    parser_simulate.add_argument('--number-of-Si-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the Silicon fluorescence energy 1.740keV' 
                        )
    parser_simulate.add_argument('--horizontal-overscan', type=int, default = '150', help = 'size of the horizontal overscan in pixels' )
    parser_simulate.add_argument('--vertical-overscan', type=int, default = '90', help = 'size of the vertical overscan in pixels' )
    parser_simulate.add_argument('--xyshape', type=tuple_of(int), default = '\"[4000,4000]\"', help = 'shape of the image as 2d pixels' )
    parser_simulate.add_argument('--depth-range', type=tuple_of(int), default = '\"[0,670]\"', help = 'range into which to randomly generate depths' )
    parser_simulate.add_argument('--diffusion-function',
                        type=str, 
                        default = Simulation.default_diffusion_function,
                        help = 'function to map z-depth into sigma' 
                        )
    parser_simulate.add_argument('--charge-efficiency-function',
                        type=str,
                        default = Simulation.default_charge_efficiency_function,
                        help = 'function to map z-depth into sigma' 
                        )
    parser_simulate.add_argument('--vertical-modulation-function',
                        type=str, 
                        default = Simulation.default_vertical_modulation_function,
                        help = 'function to modulate the vertical axis' 
                        )    
    parser_simulate.add_argument('--horizontal-modulation-function',
                        type=str, 
                        default = Simulation.default_horizontal_modulation_function,
                        help = 'function to modulate the horizontal axis' 
                        )
    parser_simulate.set_defaults(func=simulate)
    
    
    
    parser_analyse = subparsers.add_parser('analyse', help='analyse image')
    parser_analyse.add_argument('--input-fits', type=str, default = '/share/storage2/connie/data/runs/029/runID_029_03326_Int-400_Exp-10800_11Mar18_18:06_to_11Mar18_21:10_p*.fits.fz', help = 'fits file input' )
    #parser.add_argument('--input-fits', type=str, default = '/share/storage2/connie/data/runs/047/runID_047_12750_Int-400_Exp-3600_22Mar20_08:05_to_22Mar20_09:09_p1.fits.fz', help = 'fits file input' )
    parser_analyse.add_argument('--ohdu', type=int, default = '3', help = 'ohdu to be read' )
    parser_analyse.add_argument('--extract', type=bool, default = True, help = 'set to False to skip the extraction' )
    parser_analyse.add_argument('--image-energy-spectrum', type=str, default = 'image_energy_spectrum.png', help = 'set to "none" not to plot image energy spectrum' )
    parser_analyse.set_defaults(func=analyse)



    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
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
    
    args.func( args )
    
