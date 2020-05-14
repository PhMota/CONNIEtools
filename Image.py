import numpy as np
import astropy.io.fits
import scipy.stats
import scipy.ndimage

from SimulateImage import Hits
from Timer import Timer

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

class Image( np.ndarray ):
    #def __new__( cls, input_array ):
        #obj = np.asarray(input_array).view(cls)
        #return obj

    def __new__( cls, input_array ):
        obj = np.asarray(input_array).view(cls)
        return obj
        
    def save_fits(self, output):
        primary = astropy.io.fits.PrimaryHDU()
        fits_image = astropy.io.fits.ImageHDU( self )
        astropy.io.fits.HDUList([primary, fits_image]).writeto( output, overwrite=True )
    
    def stats_robust(self, gain, lambda_, sigma, tol = 1e-3, fraction = .001):
        x = np.sort(self.flatten())
        left = 0
        right = len(x)
        i = 0
        compute_gain = lambda x: (np.std(x)**2 - sigma**2)/np.mean(x)
        compute_lambda = lambda x: np.mean(x)**2/(np.std(x)**2 - sigma**2)
        mean = [np.mean(x)]
        median = [np.median(x)]
        std = [np.std(x)]
        mad = [ MAD(x) ]
        gains = [ compute_gain(x) ]
        lambdas = [ compute_lambda(x) ]
        goal_std = np.sqrt( gain**2*lambda_ + sigma**2 )
        goal_mean = gain * lambda_
        while True:
            if mean[-1] > median[-1]:
                right -= int( fraction*(right-left) )
            else:
                left += int( fraction*(right-left) )
            x__ = x[left:right]
            mean.append( np.mean( x__ ) )
            median.append( np.median( x__ ) )
            std.append( np.std( x__ ) )
            mad.append( MAD( x__ ) )
            gains.append( compute_gain( x__ ) )
            lambdas.append( compute_lambda( x__ ) )
            print( 'expected', goal_mean, goal_std )
            print( 'lim', left, right )
            print( 'stats', np.mean(mean), np.mean(median), np.mean(std), np.mean(mad) )
            print( 'vals', gain, lambda_ )
            print( 'stats2', np.mean(gains), np.mean(lambdas) )
            print()
            if len(mean) > 3:
                if abs(np.mean(mean) - np.mean(mean[:-1])) < tol:
                    break
            
        return np.mean( mean ), np.mean( median ), np.mean( std ), np.mean( mad ), x[left:right]

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
 
