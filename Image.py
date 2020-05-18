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

from SimulateImage import Hits
from Timer import Timer

bias_from_width = { 9340: 450, 8540: 150 }
bin_from_width = { 9340: 5, 8540: 1 }
vbias_from_height = { 900: 70, 900-1: 70, 4*1055: 90, 4*1055-4: 90 }

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

def crop_var_1d( x, step = 1 ):
    x2 = crop_mean_1d(x**2, step)
    m = crop_mean_1d(x, step)**2
    return x2 - m**2    

def crop_mean_1d( x, step = 1 ):
    x = np.sort( x.flatten() )
    left, right = 0, len(x)

    mean = np.mean(x)
    median = np.median(x)

    while abs(mean - median) >= .5:
        if mean > median: right -= step
        else: left += step

        mean = np.mean( x[left:right] )
        median = np.median( x[left:right] )
    return mean

def crop_mean_1d_partition( x, partitions = 10 ):
    mean = np.mean(x)
    median = np.median(x)

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
        print( 'stats', np.mean( x ), np.median( x ), np.std( x ), MAD( x ), len(x) )
        display( [x], 'E', nbins = 100, log=True )
        #raw_input()
    return mean

def stats_robust_2(self, gain, lambda_, sigma, tol = 1e-3, partitions = 10):
    x = self.flatten()

    compute_gain = lambda x: (np.std(x)**2 - sigma**2)/np.mean(x)
    compute_lambda = lambda x: np.mean(x)**2/(np.std(x)**2 - sigma**2)

    mean = [np.mean(x)]
    median = [np.median(x)]
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

        mean.append( np.mean( x ) )
        median.append( np.median( x ) )
        std.append( np.std( x ) )
        mad.append( MAD( x ) )
        gains.append( compute_gain( x ) )
        lambdas.append( compute_lambda( x ) )

        print( 'expected', goal_mean, goal_std )
        print( 'stats', np.mean(mean), np.mean(median), np.mean(std), np.mean(mad) )
        print( 'vals', gain, lambda_ )
        print( 'stats2', np.mean(gains), np.mean(lambdas) )
        print()
    
    return np.mean( mean ), np.mean( median ), np.mean( std ), np.mean( mad ), x


def list_params( fits, param = 'OHDU' ):
    return [ fit.header[ param ] for fit in fits ]

def get_imageHDU( fits, ohdu ):
    index = list_params( fits ).index( ohdu )
    return fits[ index ]

def parse_expr( expr, ohdu, correction_mode = 'none', smoothening_length = 0 ):
    with Timer('read and process') as t:
        paths = glob( expr )[::-1]
        header = None
        data = None
        for path in paths:
            datum, header = parse_fits( path, ohdu, correction_mode, smoothening_length )
            if data:
                data = { key: np.concatenate( [ data[key], datum[key] ], axis=0 ) for key in data.keys() }
            else:
                data = datum
        data['bias'] = data['bias'].view(Image)
        vbias_width = vbias_from_height[ data['data'].shape[0] ]
        data['vbias'] = data['data'][:vbias_width,:].view(Image)
        data['data'] = data['data'][vbias_width:,:].view(Image)
        data['data'], data['vbias'] = correct_rows_by_bias( [ data['data'], data['vbias'] ], data['vbias'], mode = correction_mode, size = smoothening_length )
        #print( 'corr', type( data['data'] ) )
    return data, dict(header)

def parse_fits( path, ohdu, correction_mode = 'none', smoothening_length = 0 ):
    fits_dir = astropy.io.fits.open( path )
    imageHDU = get_imageHDU( fits_dir, ohdu )
    #print( 'full data', imageHDU.shape )
    section = {}
    half_width = imageHDU.data.shape[1]/2
    bias_width = bias_from_width[ imageHDU.data.shape[1] ]
    bin_width = bin_from_width[ imageHDU.data.shape[1] ]

    section['data'] = Image( imageHDU.data[ :-1, : half_width - bias_width ] )/bin_width
    section['bias'] = Image( imageHDU.data[ :-1, half_width -bias_width: half_width ] )/bin_width
    
    if correction_mode == 'none':
        section['data'], section['bias'] = correct_global_by_bias( [section['data'], section['bias']], section['bias'] )
    else:
        section['data'], section['bias'] = correct_lines_by_bias( [section['data'], section['bias']], section['bias'], mode = correction_mode, size = smoothening_length )
    return section, dict(imageHDU.header)

def correct_global_by_bias( data, bias ):
    correction = np.median( bias )
    return [ datum - correction for datum in data ]

def correct_lines_by_bias( data, bias, size = None, mode = 'median' ):
    func = np.mean
    if mode == 'median':
        func = np.median
    elif mode == 'crop_mean':
        func = crop_mean
    correction_of_lines = func( bias, axis = 1 )
    if size: 
        correction_of_lines = scipy.ndimage.filters.uniform_filter1d( correction_of_lines, size = size, mode='nearest' )
    return [ datum - correction_of_lines[:,None] for datum in data ]

def correct_rows_by_bias( data, bias, size = None, mode = 'median' ):
    func = np.mean
    if mode == 'median':
        func = np.median
    elif mode == 'crop_mean':
        func = crop_mean
    correction_of_rows = func( bias, axis = 0 )
    if size: 
        correction_of_rows = scipy.ndimage.filters.uniform_filter1d( correction_of_rows, size = size, mode='nearest' )
    return [ datum - correction_of_rows[None,:] for datum in data ]

def MAD( data, axis = None, scale=1.4826 ):
    if axis == 0:
        median = np.median( data, axis = 0 )[None,:]
    elif axis == 1:
        median = np.median( data, axis = 1 )[:,None]
    elif axis is None:
        median = np.median( data, axis = None )
    else:
        raise Exception( 'axis not found', axis )
    return scale * np.median( np.abs( data - median ), axis = axis )

def label_clusters( condition_image ):
    s33 = [[1,1,1], [1,1,1], [1,1,1]]
    return scipy.ndimage.label( condition_image, structure=s33 )[0]

def generate_image_from_fits():
    horizontal_overscan = data[ get_overscan_slice ]

class Image( np.ndarray ):
    def __new__( cls, input_array ):
        obj = np.asarray(input_array).view(cls)
        return obj
        
    def save_fits(self, output):
        primary = astropy.io.fits.PrimaryHDU()
        fits_image = astropy.io.fits.ImageHDU( self )
        astropy.io.fits.HDUList([primary, fits_image]).writeto( output, overwrite=True )

    def sigma( self, mode = 'std' ):
        if mode == 'std':
            return np.std( self.flatten() )
        if mode == 'mad':
            return MAD( self.flatten() )
        if mode == 'crop':
            return np.sqrt( crop_var_1d( self.flatten(), step = 100 ) )
        if mode == 'partition':
            return crop_std( self.flatten(), mode = 'partition' )
    
    def mu( self, mode = 'mean' ):
        if mode == 'mean':
            return np.mean( self.flatten() )
        if mode == 'median':
            return np.median( self.flatten() )
        if mode == 'crop':
            return crop_mean_1d( self.flatten(), step = 1000 )
        if mode == 'partition':
            return crop_mean_1d_partition( self.flatten() )
    
    def median( self ): return np.median(self)
    def mean(self): return np.mean(self)
    def MAD(self): return MAD(self)
    def var(self): return np.var(self)
    
    def get_binned_distribution( self, binsize = 1):
        bins = np.arange( self.min(), self.max(), binsize)
        distribution, _ = np.histogram(data, bins)
        return .5*(bins[:-1] + bins[1:]), distribution
        
    def fit_binned( self, func, p0, bounds, binsize = 1):
        x, y = self.get_binned_distribution( binsize )
        p, pcov = scipy.optimize.curve_fit( func, x, y, p0=p0, bounds=bounds )
        return p, p_err

    def fit_norm_binned( self, binsize = 1):
        func = lambda x, amplitude, mu, sigma: amplitude * scipy.stats.norm.pdf( x, mu, sigma )
        p0 = [ len(self), self.median(), self.MAD() ]
        bounds = [(0, -np.inf, 0), (np.inf, np.inf, np.inf)]
        self.fit_binned( func, p0, bounds, binsize )
        
    
    def extract_hits( self, mode, **kwargs ):
        if mode == 'cluster':
            return self._extract_clusters_( **kwargs )

    #def _label_clusters_above_threshold_( self, threshold ):
        #return label_clusters( self.image >= thr )
    
    #def _compute_distances_to_clusters_( self ):
        #return scipy.ndimage.distance_transform_edt(  )

    def _extract_clusters_( self, threshold, border ):
        labeled_clusters = label_clusters( self >= threshold )
        print( 'number_of_clusters_above_threshold', labeled_clusters.max(), len(np.unique(labeled_clusters)) )
        is_cluster = labeled_clusters > 0
        distances_to_cluster = scipy.ndimage.distance_transform_edt( is_cluster == False )
        labeled_clusters = label_clusters( distances_to_cluster <= border )
        print( 'number_of_clusters_with border', labeled_clusters.max() )
        
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
            func=lambda v: v, 
            default=-1, 
            out_dtype=list 
            )
        def process( cluster ):
            ei, level = cluster
            x, y = np.unravel_index(ei[1], self.shape)
            return ei[0], x, y, level
        list_of_clusters = map( process, zip(list_of_clusters, levels) )
        
        return Hits( list_of_clusters, levels, threshold, border )

    def spectrum( self, output, gain, lambda_, sigma, binsize = 2 ):
        from matplotlib import pylab as plt
        
        median = np.median( self )
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
                vmin = np.median(datum)-delta
                vmax = np.median(datum)+delta
            ax.imshow( im, cmap='Blues', vmin=vmin, vmax=vmax)
        if mode == 'x':
            ax.plot(self[0,:], '.')
            ax.plot(self[-1,:], '.')
            ax.set_ylim( (np.median(self)-delta,np.median(self)+delta) )
        if mode == 'y':
            ax.plot(self[:, 0], '.')
            ax.plot(self[:, -1], '.')
            ax.set_ylim( (np.median(self)-delta,np.median(self)+delta) )
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
            ax.plot( np.median(datum, axis=0), '.', label = labels[i] )
        if mode == 'mediany':
            ax.plot( np.median(datum, axis=1), '.', label = labels[i] )
        if mode == 'xy':
            vmin = None
            vmax = None
            if delta:
                vmin = np.median(datum)-delta
                vmax = np.median(datum)+delta
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
                bins = np.linspace( np.median(datum)-delta, np.median(datum)+delta, nbins )
            ax.hist( datum.flatten(), bins = bins, histtype = 'step', label = labels[i] )
            if log:
                ax.set_yscale('log')
    legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper right' )
    plt.show()

def main( args ):
    #section_none, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'none' )
    #section_mean, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'mean' )
    #section_mean10, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'mean', smoothening_length = 10 )
    section_median, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'median' )
    #section_median10, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'median', smoothening_length = 10 )
    #section_crop, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'crop_mean' )
    section_crop10, header = parse_expr( args['input_fits'], args['ohdu'], correction_mode = 'crop_mean', smoothening_length = 10 )

    ##display( [section_none['data']], 'xy', delta = 5, labels=['none'] )
    ##display( [section_mean['data']], 'xy', delta = 5, labels=['mean'] )
    ###display( [section_mean10['data']], 'xy', delta = 5, labels=['mean10'] )
    #display( [section_median['data']], 'xy', delta = 5, labels=['median'] )
    ###display( [section_median10['data']], 'xy', delta = 5, labels=['median10'] )
    ##display( [section_crop['data']], 'xy', delta = 5, labels=['crop'] )
    #display( [section_crop10['data']], 'xy', delta = 5, labels=['crop10'] )

    print( 'bias mu', section_median['bias'].mu() )
    print( 'data mu', section_median['data'].mu() )
    print( 'bias std', section_median['bias'].sigma() )
    print( 'data std', section_median['data'].sigma() )
    print( 'bias std partition', section_median['bias'].sigma('partition') )
    print( 'data std partition', section_median['data'].sigma('partition') )
    exit(0)
    secs =[
        'bias',
        'data'
        ]
    for sec in secs:
        data = [
                #section_none[ sec ],
                section_median[ sec ],
                #section_median10[ sec ],
                #section_mean[ sec ],
                #section_mean10[ sec ],
                #section_crop[ sec ],
                section_crop10[ sec ],
            ]
        labels = [
            #'none', 
            'median', 
            #'median10',
            #'mean',
            #'mean10',
            #'crop',
            'crop10'
            ]
        #for datum, label in zip(data, labels):
            #print( label, datum.mu() )
            #print( datum.mu( mode = 'median' ) )
            #print( datum.mu( mode = 'crop' ) )
            #print( label, datum.sigma() )
            #print( datum.sigma( mode = 'mad' ) )
            #print( datum.sigma( mode = 'crop' ) )
        display( data, mode = 'E', delta = 30, labels=labels )
        display( data, mode = 'medianx', delta = 10, labels=labels )
        display( data, mode = 'mediany', delta = 10, labels=labels )
    
    exit()
    return
 
if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description = 'read fits image',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        )
    def tuple_of_int( x ):
        return map(int, eval(x))
    #parser.add_argument('--input-fits', type=str, default = '/share/storage2/connie/data/runs/029/runID_029_03326_Int-400_Exp-10800_11Mar18_18:06_to_11Mar18_21:10_p*.fits.fz', help = 'fits file input' )
    parser.add_argument('--input-fits', type=str, default = '/share/storage2/connie/data/runs/047/runID_047_12750_Int-400_Exp-3600_22Mar20_08:05_to_22Mar20_09:09_p1.fits.fz', help = 'fits file input' )
    parser.add_argument('--ohdu', type=int, default = '3', help = 'ohdu to be read' )

    #parser.add_argument('--shape', type=tuple_of_int, default = '[4000,4000,670]', help = 'shape of the image in pixel per pixel per µm' )
    #parser.add_argument('--charge-range', type=tuple_of_int, default = '[5, 200]', help = 'range into which to randomly generate charges' )
    #parser.add_argument('--charge-gain', type=eval, default = '7.25', help = 'factor to convert charges into ADU' )
    #parser.add_argument('--readout-noise', type=eval, default = '0', help = 'sigma of the normal noise distribution in ADU' )
    #parser.add_argument('--dark-current', type=eval, default = '0', help = 'lambda of Poisson distribution dimensionless' )
    #parser.add_argument('--extraction-threshold', type=eval, default = '15*4', help = 'energy threshold for extraction in ADU' )
    #parser.add_argument('--extraction-border', type=int, default = '3', help = 'borders to be added around the axtracted event' )
    #parser.add_argument('--image-fits-output', type=str, default = 'simulation.fits', help = 'set to "none" not to generate a fits output' )
    #parser.add_argument('--image-fits-input', type=str, default = 'none', help = 'set to <path> to load existing fits image' )
    parser.add_argument('--extract', type=bool, default = True, help = 'set to False to skip the extraction' )
    parser.add_argument('--image-energy-spectrum', type=str, default = 'image_energy_spectrum.png', help = 'set to "none" not to plot image energy spectrum' )
    #parser.add_argument('--reconstruction-image', type=str, default = 'reconstruction_image.pdf', help = 'set to "none" not to plot the reconstruction image' )
    #parser.add_argument('--reconstruction-spectra', type=str, default = 'reconstruction_spectra.png', help = 'set to "none" not to plot the reconstruction spectra' )
    #parser.add_argument('--sigma-like-level', type=int, default = '2', help = 'level used to compute the sigma like' )
    #parser.add_argument('--sigma-like-tol', type=float, default = '1e-3', help = 'tolerance used to compute the sigma like' )
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )
    
    main( vars(args) )
