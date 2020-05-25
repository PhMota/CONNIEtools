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
from Timer import Timer

bias_from_width = { 9340: 450, 8540: 150 }
bin_from_width = { 9340: 5, 8540: 1 }
vbias_from_height = { 900: 70, 900-1: 70, 4*1055: 90, 4*1055-4: 90 }

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

    mean = np.mean(x)
    median = np.median(x)

    while abs(mean - median) >= .5:
        if mean > median: right -= step
        else: left += step

        mean = np.mean( x[left:right] )
        median = np.median( x[left:right] )
    return mean

def crop_std_1d( x, step = 1, abs_tol = .5 ):
    x = np.sort( x.flatten() )
    left, right = 0, len(x)

    std = np.std(x)
    mad = MAD(x)
    mean = np.mean(x)
    median = np.median(x)

    while abs(std - mad) >= abs_tol:
        if mean > median: right -= step
        else: left += step

        mean = np.mean( x[left:right] )
        median = np.median( x[left:right] )
        std = np.std(x[left:right])
        mad = MAD(x[left:right])
        if (right-left) < .7*len(x):
            print('did not converge', mean, median, std, mad, right-left, len(x) )
            break
    print( '\tcrop var', std, mad )
    return std

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

def statistics( data ):
    return (float(np.mean(data)), 
            float(np.median(data)), 
            float(crop_mean_1d(data, step=10)), 
            float(np.std(data)), 
            float( MAD(data) ), 
            float(crop_std_1d(data, step=10, abs_tol = 1)),
            )

def parse_image( imagesHDU, correct_lines_by_bias_median = np.median, smoothening_length = 0, correct_rows_by_bias_median = None, row_correction_function = np.median ):
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

def get_params( imagesHDU, correct_lines_by_bias_median = np.median, correct_rows_by_bias_median = None, smoothening_length = 0, remove_hits = False, threshold = 0, border = 0, mean_function = np.median, std_function = MAD, median_diff_function = None, var_diff_function = None ):
    
    data, header = parse_image( imagesHDU, correct_lines_by_bias_median = correct_lines_by_bias_median, smoothening_length = smoothening_length, correct_rows_by_bias_median = correct_rows_by_bias_median )
    
    mu = float( mean_function( data['bias'] ) )
    sigma = float( std_function( data['bias'] ) )
    gain = 1
    if remove_hits:
        data['data'] = data['data'].get_background( threshold, border )
    if 'GAIN' in header:
        gain = header['GAIN']
        #print( 'gain', gain )
    g_lamb = median_diff_function( data['data'], data['bias'] )
    g2_lamb = var_diff_function( data['data'], data['bias'] )
    g = g2_lamb/g_lamb
    return {'mu':mu, 'sigma':sigma, 'g*lambda':g_lamb, 'g**2*lambda':g2_lamb/gain**2, 'g':g}

def read_image_from_fits( path, ohdu ):
    paths = glob( expr )[::-1]
    imagesHDU = []
    for path in paths:
        fits_dir = astropy.io.fits.open( path )
        imagesHDU.append( get_imageHDU( fits_dir, ohdu ) )
    return imagesHDU
    
def parse_single_image( imageHDU, correct_lines_by_bias_median = np.median, smoothening_length = 0 ):
    section = {}
    if imageHDU.header['ohdu'] == -1:
        half_width = imageHDU.data.shape[1]
        bias_width = imageHDU.header['biasW']
        bin_width = imageHDU.header['rebin']
    else:
        half_width = imageHDU.data.shape[1]/2
        bias_width = bias_from_width[ imageHDU.data.shape[1] ]
        bin_width = bin_from_width[ imageHDU.data.shape[1] ]

    section['data'] = Image( imageHDU.data[ :-1, :half_width - bias_width ] )#/bin_width
    section['bias'] = Image( imageHDU.data[ :-1, half_width -bias_width: half_width ] )#/bin_width

    if correct_lines_by_bias_median is None:
        section['data'], section['bias'] = correct_global_by_bias( [section['data'], section['bias']], section['bias'] )
    else:
        section['data'], section['bias'] = correct_lines_by_bias( [section['data'], section['bias']], section['bias'], func = correct_lines_by_bias_median, size = smoothening_length )
    return section, dict(imageHDU.header)

def correct_global_by_bias( data, bias ):
    correction = np.median( bias )
    return [ datum - correction for datum in data ]

def correct_lines_by_bias( data, bias, size = 0, func = np.median ):
    correction_of_lines = func( bias, axis = 1 )
    if size > 2:
        correction_of_lines = scipy.ndimage.filters.uniform_filter1d( correction_of_lines, size = size, mode='nearest' )
    return [ datum - correction_of_lines[:,None] for datum in data ]

def correct_rows_by_bias( data, bias, size = 0, func = np.median ):
    correction_of_rows = func( bias, axis = 0 )
    #print( 'correction_of_rows', correction_of_rows.shape )
    if size > 2: 
        correction_of_rows = scipy.ndimage.filters.uniform_filter1d( correction_of_rows, size = size, mode='nearest' )
    return [ datum - correction_of_rows[None,:] for datum in data ]

def label_clusters( condition_image ):
    s33 = [[1,1,1], [1,1,1], [1,1,1]]
    return scipy.ndimage.label( condition_image, structure=s33 )[0]


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
    
    #def median( self ): return np.median(self)
    #def mean(self): return np.mean(self)
    #def MAD(self): return MAD(self)
    #def var(self): return np.var(self)
    
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
            vmin = datum.min()
            vmax = datum.max()
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

def simulate_and_test( correct_lines_by_bias_median, correct_rows_by_bias_median, smoothening_length, mean_function, std_function, median_diff_function, var_diff_function ):
    import Simulation
    output_file = 'simulate_and_test.csv'
    count = 0
    columns = ['n', 'readout_noise', 'charge_gain*dark_current', 'charge_gain**2*dark_current', 'charge_gain', 'mu', 'sigma', 'g*lambda', 'g**2*lambda', 'g']
    open( output_file, 'w' ).writelines( '#' + ', '.join(columns) + '\n' )
    while count < 100:
        count += 1
        print( 'count', count )
        args = {'simulation_output': 'image_test.csv',
                'image_fits_output': 'image_test.fits',
                'rebin': [1,1],
                'xyshape': [4000,4000],
                'number_of_charges': 0,
                'depth_range': [0,670],
                'charge_range': [1,2000],
                'charge_gain': 10.*np.random.random(),
                'readout_noise': 18.*np.random.random(),
                'dark_current': 0.2*np.random.random(),
                'horizontal_overscan': 150,
                'vertical_overscan': 90,
                'diffusion_function': Simulation.default_diffusion_function,
                'charge_efficiency_function': Simulation.default_charge_efficiency_function,
                'horizontal_modulation_function': Simulation.default_horizontal_modulation_function,
                'vertical_modulation_function': Simulation.default_vertical_modulation_function,
                }
        with Timer('simulate'):
            Simulation.simulate_events( args )
            sim = Simulation.Simulation( args['simulation_output'] )
            imageHDU = sim.generate_image( return_image = True )
        with Timer('get params'):
            params = get_params( 
                [imageHDU],
                correct_lines_by_bias_median = correct_lines_by_bias_median,
                correct_rows_by_bias_median = correct_rows_by_bias_median,
                smoothening_length = smoothening_length,
                mean_function = mean_function,
                std_function = std_function,
                median_diff_function = median_diff_function,
                var_diff_function = var_diff_function,
            )
        args['charge_gain*dark_current'] = args['charge_gain']*args['dark_current']
        args['charge_gain**2*dark_current'] = args['charge_gain']**2*args['dark_current']
        args['n'] = count
        entry = [ args[column] if column in args else params[column] for column in columns ]
        open( output_file, 'a' ).writelines( ', '.join(map( str, entry)) + '\n' )
        #print( 'noise', args['readout_noise'], params['sigma'] )
        #print( 'darkcurrent', args['dark_current'], params['g*lambda/gain'], params['g**2*lambda/gain**2'] )
        #print( 'gain', args['charge_gain'], params['g'] )

def main( args ):

    def median_by_lines( x ):
        medians = np.median(x, axis=1)
        return np.mean( medians )

    def MAD_by_lines2( x ):
        correction = np.median( x, axis=0 )
        #print( 'correction2.shape', correction.shape )
        mads_by_line = MAD( x - correction[None,:], axis = 1 )
        return np.mean( mads_by_line )

    def MAD_by_lines( x ):
        #print( 'correction.shape', correct_rows_by_bias( [x], x, func = np.median)[0].shape )
        return np.mean( MAD(correct_rows_by_bias( [x], x, func = np.median )[0], axis=1) )
    
    def median_diff_by_half_lines( x, y ):
        #display([x], mode = 'medianx')
        #display([x], mode = 'mediany')
        return float(np.mean( np.median(x[:,x.shape[1]/2:], axis=1) - np.median(y, axis=1) ))

    def MADsqr_diff_by_half_lines( x, y ):
        return float(np.mean( MAD(x[:,x.shape[1]/2:], axis=1)**2 - MAD(y, axis=1)**2 ))
    
    simulate_and_test( 
        correct_lines_by_bias_median = np.median, 
        correct_rows_by_bias_median = None,
        smoothening_length = 0,
        mean_function = median_by_lines,
        std_function = MAD_by_lines,
        median_diff_function = median_diff_by_half_lines,
        var_diff_function = MADsqr_diff_by_half_lines,
        )
    
    exit(0)
    
    print( 'official median', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        correct_lines_by_bias_median = np.median, 
        correct_rows_by_bias_median = None,
        smoothening_length = 10,
        mean_function = median_by_lines,
        std_function = MAD_by_lines,
        median_diff_function = median_diff_by_half_lines,
        var_diff_function = MADsqr_diff_by_half_lines,
    ))

    print( 'vertical median', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        correct_lines_by_bias_median = np.median, 
        correct_rows_by_bias_median = np.median,
        smoothening_length = 0,
        mean_function = median_by_lines,
        std_function = MAD_by_lines,
        median_diff_function = median_diff_by_half_lines,
        var_diff_function = MADsqr_diff_by_half_lines,
    ))

    print( 'vertical median', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        correct_lines_by_bias_median = np.median, 
        correct_rows_by_bias_median = np.median,
        smoothening_length = 10,
        mean_function = median_by_lines,
        std_function = MAD_by_lines,
        median_diff_function = median_diff_by_half_lines,
        var_diff_function = MADsqr_diff_by_half_lines,
    ))

    print( 'vertical mean', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        correct_lines_by_bias_median = np.mean, 
        correct_rows_by_bias_median = np.mean,
        smoothening_length = 0,
        mean_function = median_by_lines,
        std_function = MAD_by_lines,
        median_diff_function = median_diff_by_half_lines,
        var_diff_function = MADsqr_diff_by_half_lines,
    ))

    print( 'vertical mean', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        correct_lines_by_bias_median = np.mean,
        correct_rows_by_bias_median = np.mean,
        smoothening_length = 0,
        mean_function = median_by_lines,
        std_function = MAD_by_lines,
        median_diff_function = median_diff_by_half_lines,
        var_diff_function = MADsqr_diff_by_half_lines,
    ))
    
    exit(0)
    
    print( 'correctbias', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        line_correction_function = np.median, 
        correct_vertical = False,
        mean_function = median_by_lines,
        std_function = MAD_by_lines2
    ))
    print( 'correctbias vertical', get_params( 
        args['input_fits'], 
        args['ohdu'], 
        line_correction_function = np.median, 
        correct_vertical = True,
        mean_function = median_by_lines,
        std_function = MAD_by_lines2
    ))
    exit(0)
    print( 'median', get_params( args['input_fits'], args['ohdu'], correction_mode = 'median', correct_vertical = False ) )
    print( 'median100', get_params( args['input_fits'], args['ohdu'], correction_mode = 'median', correct_vertical = True, smoothening_length = 10 ) )
    
    exit(0)
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
    parser.add_argument('--input-fits', type=str, default = '/share/storage2/connie/data/runs/029/runID_029_03326_Int-400_Exp-10800_11Mar18_18:06_to_11Mar18_21:10_p*.fits.fz', help = 'fits file input' )
    #parser.add_argument('--input-fits', type=str, default = '/share/storage2/connie/data/runs/047/runID_047_12750_Int-400_Exp-3600_22Mar20_08:05_to_22Mar20_09:09_p1.fits.fz', help = 'fits file input' )
    parser.add_argument('--ohdu', type=int, default = '3', help = 'ohdu to be read' )

    parser.add_argument('--extract', type=bool, default = True, help = 'set to False to skip the extraction' )
    parser.add_argument('--image-energy-spectrum', type=str, default = 'image_energy_spectrum.png', help = 'set to "none" not to plot image energy spectrum' )
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )
    
    main( vars(args) )
