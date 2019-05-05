import astropy
import astropy.io
import astropy.io.fits

import numpy as np
import scipy
import scipy.stats
from scipy.misc import factorial
import scipy.optimize

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import re
import glob
import datetime
import json
import os
import sys
import StringIO

#import bokeh.plotting
#import statsmodel as sm

def std_w(values, weights):
    """
    Return the weighted standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)

def negLogLikelihood( params, data ):
    return -np.sum( np.log( poisson(data, params[0]) ) )
    
def poisson(k, lambda_, shift = 0):
    #k -= shift
    return np.concatenate( [ np.zeros_like(k[k<0]), np.exp( k[k>=0]*np.log(lambda_) -lambda_ - scipy.special.gammaln(k[k>=0]+1)) ])

def poissonf(k, lamb, shift = 0):
    if k < 0: return 0
    return np.exp( k*np.log(lamb) -lamb - scipy.special.gammaln(k+1))

#def gaussian_poisson( k, lamb = 1, scale = 1 ):
    #return scipy.stats.norm.pdf(k, scale=scale) + poisson( k, lamb=lamb )

def check_new_files():
    pass

def fix_index( index_string ):
    '''
    converts the indexing string of the form [xmin:xmax,ymin:ymax] where 0 < xmin,xmax,ymin,ymax <= MAX (fits indexing) into [xmin-1:xmax,ymin-1:ymax] where 0 <= xmin,ymin < MAX and 0 < xmax,ymax <= MAX (python indexing standard)
    '''
    dims = np.array(map( lambda s: int(s), re.match( r'\[(\d+):(\d+)\,(\d+):(\d+)\]', index_string ).groups() ))
    dims[ np.array([True, False, True, False]) ] -= 1
    return dims
    
def apply_cut( array, mask_string ):
    '''
    access the indexes as strings to numpy indexes and transposes the array to comply with the numpy standard
    '''
    t = tuple(fix_index(mask_string))
    return eval( 'array.T' + '[%i:%i,%i:%i]' % t ), t
    
def manual_conv( x, lamb, loc, scale):
    rvg = scipy.stats.norm(scale = scale)
    rvp = scipy.stats.poisson(lamb)
    y = range(100)
    pmfs = scipy.stats.poisson.pmf(y, lamb)
    print 'pmfs', zip(y,pmfs)
    return

def test_convolution():
    sigma = 10.
    
    x = np.linspace( -5*sigma, 10*sigma, 100 )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rvg = scipy.stats.norm(scale = sigma)
    #ax.plot( x, rvg.pdf(x), lw =2, label='gaussian' )
    
    #conv = lambda x, lamb, loc, scale: np.sum(  )
    
    manual_conv(x, 10, 0, sigma )
    
    for lambda_ in range(1, 10, 1):
    
        
        rvp = scipy.stats.poisson(lambda_)
        print lambda_
        #print rvp.pmf(x).shape, np.sum(rvp.pmf(x))
        
        conv = lambda x: np.sum( rvp.pmf(range(100)) * rvg.pdf( x[np.newaxis,:] - np.array(range(100))[:,np.newaxis] ), axis = 0 )
        conv2 = lambda x: [ scipy.integrate.quad( lambda y: poissonf(y, lambda_)*gaussian(ix-y,scale=sigma), -np.inf, np.inf )[0] for ix in x]
        #ax.plot(x, [ conv(m) for m in x ], label='manual' )
        
        #K = .5
        #ax.plot(x, scipy.stats.exponnorm.pdf(x, K = lambda_/100., scale = sigma ), label='exponnorm' )
        #print conv(x)
        ax.plot(x, conv(x), label='sum' )

        #ax.plot(x, np.convolve( rvp.pmf(range(10)), rvg.pdf(x), mode='same')/(np.sum(rvp.pmf(range(10)))), label = 'conv' )
        
        #frvp = np.fft.fft( rvp.pmf(range(10)), n = len(x) )
        #frvg = np.fft.fft( rvg.pdf(x) )
        
        #ax.plot(x, np.fft.ifft( frvp*frvg )/(np.sum(rvp.pmf(range(10)))), lw = 2, label='fft' )
        
    plt.show()
    return

def compute_values( input_file, database, complete = False ):
    '''
    compute the observables in list_of_observables on each image of an input_file and returns a dict of the form {'input_file':{'observable-name': observable-return, ... }}
    '''
    name = input_file.split('/')[-1].split('_p')[0]
    print 'name', name
    data = {}
    
    if name in database:
        print 'already done'
        data[name] = database[name]
        #return None #remove for updating
    else:
        data[name] = {}
    
    parts, processed = get_all_parts( input_file )
    if len(parts) != 4: 
        print 'somethign wrong with the number of parts'
        return None
    for part in parts:
        print 'parts', part
    print 'processed', processed
    fits_parts = [ astropy.io.fits.open( part ) for part in parts ]
    fits_processed = None
    if not processed is None:
        fits_processed = astropy.io.fits.open( processed )
        
    header = fits_parts[0][1].header
    
    if not 'date' in data[name]:
        data[name]['date'] = datetime.datetime.fromtimestamp( header['EXPSTART'] ).strftime('%d%b%y_%H:%M')
    if not 'expstart' in data[name]:
        data[name]['expstart'] = header['EXPSTART']
        
    if not 'expstop' in data[name]:
        data[name]['expstop'] = header['EXPSTOP']
        
    if not 'exptime' in data[name]:
        data[name]['exptime'] = header['EXPTIME']
        
    if not 'htrmax' in data[name]:
        data[name]['htrmax'] = header['HTRMAX']
        
    if not 'htrmin' in data[name]:
        data[name]['htrmin'] = header['HTRMIN']
        
    if not 'intw' in data[name]:
        data[name]['intw'] = header['INTW']
        
    if not 'presmax' in data[name]:
        data[name]['presmax'] = header['PRESMAX']
        
    if not 'presmin' in data[name]:
        data[name]['presmin'] = header['PRESMIN']
        
    if not 'rdstart' in data[name]:
        data[name]['rdstart'] = header['RDSTART']
        
    if not 'rdstop' in data[name]:
        data[name]['rdstop'] = header['RDSTOP']
        
    if not 'rdtime' in data[name]:
        data[name]['rdtime'] = header['RDTIME']
        
    if not 'runid' in data[name]:
        data[name]['runid'] = header['RUNID']
        
    if not 'tempMax' in data[name]:
        data[name]['tempMax'] = header['TEMPMAX']
        
    if not 'tempMin' in data[name]:
        data[name]['tempMin'] = header['TEMPMIN']
    
    done = True
    for i, hdu in enumerate(fits_parts[0]):
        if hdu.data is None:
            continue
        if not 'OHDU' in hdu.header:
            continue

        if complete: testing( fits_parts[0][5], input_file )
        ccdnum = hdu.header['CCDNUM']
        ohdu = hdu.header['OHDU']
        ohdu_str = 'ohdu%d'%ohdu
        print name, ohdu_str, ccdnum
        
        if not ohdu_str in data[name]:
            data[name][ohdu_str] = {}

        if not 'parts' in data[name][ohdu_str]:
            done = False
            data[name][ohdu_str]['parts'] = {}
            
        if not 'merged' in data[name][ohdu_str]:
            done = False
            data[name][ohdu_str]['merged'] = {}
        
        merge_hdu_header = None
        
        biasA = None
        biasB = None
        dataA = None
        dataB = None
        #if not 'peaks_dataA' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_dataB' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_biasA' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_biasB' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_dist' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_corr' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_sub' in data[name][ohdu_str]['parts'] or\
            #not 'peaks_dataA' in data[name][ohdu_str]['merged'] or\
            #not 'peaks_dataB' in data[name][ohdu_str]['merged'] or\
            #not 'peaks_biasA' in data[name][ohdu_str]['merged'] or\
            #not 'peaks_biasB' in data[name][ohdu_str]['merged'] or\
            #not 'peaks_dist' in data[name][ohdu_str]['merged'] or\
            #not 'peaks_corr' in data[name][ohdu_str]['merged'] or\
            #not 'peaks_sub' in data[name][ohdu_str]['merged']:

        if True:
            hdu_parts = [ part[i] for part in fits_parts ]
            
            #biasA = [ apply_cut( part.data, part.header['BIASSECA'] )[0].astype(float) for part in hdu_parts ]
            #biasB = [ apply_cut( part.data, part.header['BIASSECB'] )[0].astype(float) for part in hdu_parts ]
            #dataA = [ apply_cut( part.data, part.header['DATASECA'] )[0].astype(float) for part in hdu_parts ]
            #dataB = [ apply_cut( part.data, part.header['DATASECB'] )[0].astype(float) for part in hdu_parts ]
            
            
            if not 'peaks_dataA' in data[name][ohdu_str]['parts']:
                dataA = [ apply_cut( part.data, part.header['DATASECA'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_dataA'] = [ computePeakParams( dA ) for dA in dataA ]
                done = False
                print 'dataA', data[name][ohdu_str]['parts']['peaks_dataA'], part.header['DATASECA']

            if not 'peaks_dataB' in data[name][ohdu_str]['parts']:
                dataB = [ apply_cut( part.data, part.header['DATASECB'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_dataB'] = [ computePeakParams( dB ) for dB in dataB ]
                done = False
                print 'dataB', data[name][ohdu_str]['parts']['peaks_dataB'], part.header['DATASECB']
            
            if not 'peaks_biasA' in data[name][ohdu_str]['parts']:
                biasA = [ apply_cut( part.data, part.header['BIASSECA'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_biasA'] = [ computePeakParams( bA ) for bA in biasA ]
                done = False
                print 'biasA', data[name][ohdu_str]['parts']['peaks_biasA'], part.header['BIASSECA']

            if not 'peaks_biasB' in data[name][ohdu_str]['parts']:
                biasB = [ apply_cut( part.data, part.header['BIASSECB'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_biasB'] = [ computePeakParams( bB ) for bB in biasB ]
                done = False
                print 'biasB', data[name][ohdu_str]['parts']['peaks_biasB'], part.header['BIASSECB']

            if not 'peaks_dist' in data[name][ohdu_str]['parts']:
                if dataA is None: dataA = [ apply_cut( part.data, part.header['DATASECA'] )[0].astype(float) for part in hdu_parts ]
                if dataB is None: dataB = [ apply_cut( part.data, part.header['DATASECB'] )[0].astype(float) for part in hdu_parts ]
                if biasA is None: biasA = [ apply_cut( part.data, part.header['BIASSECA'] )[0].astype(float) for part in hdu_parts ]
                if biasB is None: biasB = [ apply_cut( part.data, part.header['BIASSECB'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_dist'] = [ computePeakParams( correctDist( dB, bB, dA, bA ) ) for dB, bB, dA, bA in zip(dataB, biasB, dataA, biasA) ]
                done = False
                print 'dist', data[name][ohdu_str]['parts']['peaks_dist']

            if not 'peaks_corr' in data[name][ohdu_str]['parts']:
                if dataA is None: dataA = [ apply_cut( part.data, part.header['DATASECA'] )[0].astype(float) for part in hdu_parts ]
                if dataB is None: dataB = [ apply_cut( part.data, part.header['DATASECB'] )[0].astype(float) for part in hdu_parts ]
                if biasA is None: biasA = [ apply_cut( part.data, part.header['BIASSECA'] )[0].astype(float) for part in hdu_parts ]
                if biasB is None: biasB = [ apply_cut( part.data, part.header['BIASSECB'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_corr'] = [ computePeakParams( correctDistBiasC( dB, bB, dA, bA ) ) for dB, bB, dA, bA in zip(dataB, biasB, dataA, biasA) ]
                done = False
                print 'distC', data[name][ohdu_str]['parts']['peaks_corr']

            if not 'peaks_sub' in data[name][ohdu_str]['parts']:
                if dataA is None: dataA = [ apply_cut( part.data, part.header['DATASECA'] )[0].astype(float) for part in hdu_parts ]
                if dataB is None: dataB = [ apply_cut( part.data, part.header['DATASECB'] )[0].astype(float) for part in hdu_parts ]
                if biasA is None: biasA = [ apply_cut( part.data, part.header['BIASSECA'] )[0].astype(float) for part in hdu_parts ]
                if biasB is None: biasB = [ apply_cut( part.data, part.header['BIASSECB'] )[0].astype(float) for part in hdu_parts ]
                data[name][ohdu_str]['parts']['peaks_sub'] = [ computePeakParams( correctDist( dB, bB, dA, bA ) ) for dB, bB, dA, bA in zip(dataB, biasB, dataA, biasA) ]
                done = False
                print 'distSub', data[name][ohdu_str]['parts']['peaks_sub']
        
            def make_merge():
                _merge_hdu_ = hdu_parts[0]
                _merge_hdu_.data = np.concatenate( [ part.data for part in hdu_parts ], axis = 0 )
                for tag in [ 'BIASSECA', 'BIASSECB', 'DATASECA', 'DATASECB' ]:
                    _merge_hdu_.header[tag] = _merge_hdu_.header[tag].replace('1055', '4220')
                _merge_hdu_header_ = _merge_hdu_.header['DATASECB']
                return _merge_hdu_, _merge_hdu_header_
            merge_hdu = None
            dataB = None
            dataA = None
            biasB = None
            biasA = None
            #dataB = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
            #dataA = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
            #biasB = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECB'] )[0].astype(float)
            #biasA = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECA'] )[0].astype(float)
            print 'merged'

            if not 'peaks_dataA' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                dataA = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_dataA'] = computePeakParams( dataA )
                done = False
                print 'dataA', data[name][ohdu_str]['merged']['peaks_dataA'], merge_hdu.header['DATASECA']

            if not 'peaks_dataB' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                dataB = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_dataB'] = computePeakParams( dataB )
                done = False
                print 'dataB', data[name][ohdu_str]['merged']['peaks_dataB'], merge_hdu.header['DATASECB']

            if not 'peaks_biasA' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                biasA = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECA'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_biasA'] = computePeakParams( biasA )
                done = False
                print 'biasA', data[name][ohdu_str]['merged']['peaks_biasA'], merge_hdu.header['BIASSECA']

            if not 'peaks_biasB' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                biasB = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECB'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_biasB'] = computePeakParams( biasB )
                done = False
                print 'biasB', data[name][ohdu_str]['merged']['peaks_biasB'], merge_hdu.header['BIASSECB']

            if not 'peaks_dist' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                if dataA is None: dataA = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
                if dataB is None: dataB = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
                if biasA is None: biasA = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECA'] )[0].astype(float)
                if biasB is None: biasB = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECB'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_dist'] = computePeakParams( correctDist( dataB, biasB, dataA, biasA ) )
                done = False
                print 'dist', data[name][ohdu_str]['merged']['peaks_dist']

            if not 'peaks_corr' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                if dataA is None: dataA = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
                if dataB is None: dataB = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
                if biasA is None: biasA = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECA'] )[0].astype(float)
                if biasB is None: biasB = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECB'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_corr'] = computePeakParams( correctDistBiasC( dataB, biasB, dataA, biasA ) )
                done = False
                print 'distC', data[name][ohdu_str]['merged']['peaks_corr']

            if not 'peaks_sub' in data[name][ohdu_str]['merged']:
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                if merge_hdu is None: merge_hdu, merge_hdu_header = make_merge()
                if dataA is None: dataA = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
                if dataB is None: dataB = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
                if biasA is None: biasA = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECA'] )[0].astype(float)
                if biasB is None: biasB = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECB'] )[0].astype(float)
                data[name][ohdu_str]['merged']['peaks_sub'] = computePeakParams( correctDistSub( dataB, biasB, dataA, biasA ) )
                done = False
                print 'distSub', data[name][ohdu_str]['merged']['peaks_sub']
        
        if not fits_processed is None:
            if not 'processed' in data[name][ohdu_str]:
                done = False
                data[name][ohdu_str]['processed'] = {}
            if not 'peaks' in data[name][ohdu_str]['processed']:
                done = False
                proc = [ image for image in fits_processed if image.header['OHDU'] == ohdu ]
                if proc == []:
                    print "no processed image"
                    #print ohdu, [ image.header['OHDU'] for image in fits_processed ]
                    #print ccdnum, [ image.header['CCDNUM'] for image in fits_processed ]
                    continue
                proc = proc[0]
                print 'pre header', merge_hdu_header
                if merge_hdu_header is None:
                    merge_hdu_header = fits_parts[0][1].header['DATASECB'].replace('1055', '4220')
                    print 'header', merge_hdu_header
                cut = merge_hdu_header
                cut = cut.replace('[9:4120', '[1:4112')
                
                dataB = apply_cut( proc.data, cut )[0].astype(float)
                try:
                    vmax = np.max(dataB.flatten())
                    if np.sum( dataB > vmax/2. ) > 0:
                        dataB[ dataB > vmax/2. ] = 0
                    data[name][ohdu_str]['processed']['peaks'] = computePeakParams( dataB )
                    print 'peaks', data[name][ohdu_str]['processed']['peaks']
                    print cut
                except:
                    print 'dataB cut processed', dataB.flatten().shape
                    print cut
        if done: 
            print 'all done for', name
            print 'skipping other hdus'
            data = None
            break
    return data


def gaussian( x, scale = 1, loc = 0 ):
    return scipy.stats.norm.pdf(x, loc = loc, scale = scale )

def Ngaussian( x, scale = 1, loc = 0, N = 1 ):
    return N*scipy.stats.norm.pdf(x, loc = loc, scale = scale )

def poisson_gaussian( x, loc = 0, scale = 1, shift = 0, lambda_ =1, norm = 1, factor = 1 ):
    return norm*gaussian( x, loc = loc, scale = scale ) + factor*poisson( x - shift, lambda_ = lambda_ )

def fit_distribution_poisson_gaussian( dist, tag ):
    ''' fits the distribution using a superposition of a gaussian and a poisson distributions'''
    
    # fix for bad high pixels
    while 1:
        mean = np.mean(dist)
        median = np.median(dist)
        center = .5*(np.min(dist)+np.max(dist))
        std = np.std(dist)
        
        if abs(center - median) > std:
            dist = dist[dist!=np.max(dist)]
        else:
            break
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111)

    mode = float(scipy.stats.mode(dist)[0][0])
    dist_cut_mode = dist[ dist <= mode ] - mode
    dist_mode_symmetric = np.concatenate( [ dist_cut_mode, -dist_cut_mode[ dist_cut_mode != 0 ] ] )
    if len(dist_mode_symmetric) > len(dist):
        print "potential problem!!!"
    gaussian_dist_mode = scipy.stats.norm.fit( dist_mode_symmetric )
    nbins = int(np.sqrt(len(dist))*.1)
    if not nbins > 3: nbins = 3
    hist, bins = np.histogram( dist - mode, bins = nbins )
    dbin = float( bins[1] - bins[0] )
    xbins = (bins[:-1] + bins[1:])*.5
    
    hist_symmetric, tmp = np.histogram( dist_mode_symmetric, bins )
    
    #ax.step(bins[1:], hist, label='hist')
    #ax.step(bins[1:], hist_symmetric, label='histSym')
    #ax.plot( xbins, len(dist_mode_symmetric)*dbin*gaussian(xbins, loc=gaussian_dist_mode[0], scale=gaussian_dist_mode[1] ), label='gaussian' )
    
    hist_poisson2 = hist - hist_symmetric
    hist_poisson2[hist_poisson2<0] = 0
    if np.sum( hist_poisson2 ) <= 0:
        print "negative only hist"
    #ax.step( bins[1:], hist_poisson2, label='histDiff' )
    
    std_poisson = std_w( xbins, hist_poisson2 ) if np.sum(hist_poisson2) > 0 else 1
    
    poisson_f = lambda k,par0,par1,par2: par2*np.sum(hist_poisson2)*poisson( k - par0*std_poisson, lambda_ = par1*std_poisson )
    #print 'dist', dist, max(dist)
    #print 'mode', mode
    #print 'xbins', xbins
    #print 'hist', hist
    #print 'hist_symmetric', hist_symmetric
    #print 'hist_poisson2', hist_poisson2
    try:
        params, cov_matrix = scipy.optimize.curve_fit( poisson_f, xbins, hist_poisson2/dbin )
        print 'params poisson !', params
        print params[0]*std_poisson
        #ax.plot( xbins, poisson_f(xbins, *params)*dbin, label='poisson' )
        
        complete_f = lambda x,par0,par1,par2,par3,par4,par5: poisson_gaussian(x, norm = len(dist_mode_symmetric), loc = par1+gaussian_dist_mode[0], scale = par2*gaussian_dist_mode[1], shift = (par5+params[0])*std_poisson, lambda_ = par4*params[1]*std_poisson, factor = par3*params[2]*np.sum(hist_poisson2) )
        
        params2, cov_matrix = scipy.optimize.curve_fit( complete_f, xbins, hist/dbin )
        #ax.plot( xbins, complete_f(xbins, *params2)*dbin, label='comp' )
        
        #print params2
        #ax.plot( bins, f_complete(bins, *params2), label = 'fA2', lw=2 )
        
    except RuntimeError:
        print 'did not converge'
        params = [0,1,0]
        params2 = [ 0, 1, 1, 1, 1, 0 ]
    except ValueError:
        print 'xbins', xbins
        print 'hist', hist
        print 'hist_symmetric', hist_symmetric
        print 'hist_poisson2', hist_poisson2
        params = [ 0, 1, 0]
        params2 = [ 0, 1, 1, 1, 1, 0 ]
    except TypeError:
        print 'xbins', xbins
        print 'hist', hist
        print 'hist_symmetric', hist_symmetric
        print 'hist_poisson2', hist_poisson2
        params = [ 0, 1, 0]
        params2 = [ 0, 1, 1, 1, 1, 0 ]
    
    #plt.legend()
    #plt.savefig( tag )
    #plt.close(fig)
    return [
        {
            'loc':gaussian_dist_mode[0], 
            'scale': gaussian_dist_mode[1], 
            'shift': params[0]*std_poisson, 
            'lambda_': params[1]*std_poisson, 
            'factor': params[2]*np.sum(hist_poisson2), 
            'norm': len(dist_mode_symmetric) 
         },
        {
            'loc': params2[1]+gaussian_dist_mode[0], 
            'scale': params2[2]*gaussian_dist_mode[1], 
            'shift': (params2[5]+params[0])*std_poisson, 
            'lambda_':params2[4]*params[1]*std_poisson, 
            'factor': params[2]*params2[3]*np.sum(hist_poisson2), 
            'norm': len(dist_mode_symmetric) 
        }
        ]#, complete_f, params2
    
def remove_outliers( dist ):
    
    while 1:
        ordered_indices = np.argsort( dist )
        print dist[ ordered_indices ]
        dist_diff = dist[ ordered_indices[1:] ] - dist[ ordered_indices[:-1] ]
        diff_std = np.std( dist_diff[dist_diff>0] )
        diff_median = np.median( dist_diff[dist_diff>0] )
        diff_mean = np.mean( dist_diff[dist_diff>0] )
        
        print dist_diff
        print diff_std, diff_mean10800
        N = 10
        print dist_diff[ dist_diff - diff_mean > N*diff_std ]
        args_remove = np.argwhere( dist_diff - diff_mean > N*diff_std )
        remove_order = np.argsort( dist_diff[~np.isnan(dist_diff)] )
        print 'order', dist_diff[~np.isnan(dist_diff)][ remove_order ]
        
        if len(args_remove) > 0:
            for remove_arg in args_remove:
                print remove_arg, dist[ ordered_indices[remove_arg] ], dist[ ordered_indices[remove_arg+1] ]
                if remove_arg > .5*len(dist): dist[ ordered_indices[remove_arg+1] ] = np.nan
                else: dist[ ordered_indices[remove_arg] ] = np.nan
                print remove_arg, dist[ ordered_indices[remove_arg] ], dist[ ordered_indices[remove_arg+1] ]
        else: break
    
    print 'done'
    return dist

def subdivideSection( data, n = 1, nMax = 5, E0 = None ):
    halfy = len(data[0,:])/2
    halfx = len(data[:,0])/2
    E = np.sum(data)*4**(n-1)
    if E0 is None: 
        E0 = E
        print 'E0', E0
    if E < E0*.95:
        print 'E', n, E, np.sum(data), E0/(4**(n-1)), data.shape
    if n >= nMax: return
    subdivideSection( data[:halfx,:halfy], n+1, nMax, E0 )
    subdivideSection( data[:halfx,halfy:], n+1, nMax, E0 )
    subdivideSection( data[halfx:,:halfy], n+1, nMax, E0 )
    subdivideSection( data[halfx:,halfy:], n+1, nMax, E0 )

def computePeakParams( dist, n = 1, f = 10., d=1 ):
    flat = dist.flatten()
    hist, bins = np.histogram( flat, bins = d*int(np.sqrt(len(flat))) )
    #if verbose: print hist, bins
    highCounting = np.array(hist) > max(hist)/f
    xbins = (bins[1:] + bins[:-1])*.5
    xmin, xmax = np.min( xbins[highCounting] ), np.max( xbins[highCounting] )
    #if verbose: print xmin, xmax
    subdist = dist[ np.logical_and( dist < n*xmax, dist > n*xmin ) ]
    #params = scipy.stats.norm.fit( subdist, loc = np.mean(subdist), scale = np.std(subdist) )
    try:
        params = scipy.stats.norm.fit( subdist )
    except MemoryError:
        return None
    #if verbose: print params
    return params

def subtractMedians( dist ):
    dist -= np.median( dist.flatten() )
    #dist -= np.median( dist, axis = 1 )[:,np.newaxis]
    #dist -= np.median( dist, axis = 0 )[np.newaxis,:]
    #return dist
    sub = np.zeros_like(dist) + np.median( dist, axis = 1 )[:,np.newaxis] + np.median( dist, axis = 0 )[np.newaxis,:]
    #row = np.zeros_like(dist) + np.median( dist, axis = 1 )[:,np.newaxis]
    #col = np.zeros_like(dist) + np.median( dist, axis = 0 )[np.newaxis,:]
    #dist -= (row + col)*.5
    return dist - sub/2.

def subtractBias( dist, bias ):
    return dist - np.median( bias, axis = 0 )[np.newaxis,:]

def subtractEmpty( dist, empty ):
    dist -= np.median(dist)
    empty -= np.median(empty)
    return dist - np.flipud(empty)

def correctDistMedian(dist, bias1, empty, biasEmpty ):
    return subtractEmpty( 
        subtractBias( subtractMedians(dist), subtractMedians(bias1) ), 
        subtractBias( subtractMedians(empty), subtractMedians(biasEmpty) ) )

def correctDist( dist, bias1, empty, biasEmpty ):
    return subtractEmpty( 
        subtractBias( dist, bias1 ), 
        subtractBias( empty, biasEmpty ) )

def correctDistBiasC(dist, bias1, empty, biasEmpty ):
    M1 = np.median(bias1)
    MEmpty = np.median(biasEmpty)
    bias = subtractEmpty( bias1, biasEmpty )
    
    return subtractEmpty( 
        subtractBias( dist, bias + M1 ), 
        subtractBias( empty, np.flipud(bias + MEmpty) ) )

def extract_hdu( hdu ):
    biasA = apply_cut( hdu.data, hdu.header['BIASSECA'] )[0].astype(float)
    biasB = apply_cut( hdu.data, hdu.header['BIASSECB'] )[0].astype(float)
    dataA = apply_cut( hdu.data, hdu.header['DATASECA'] )[0].astype(float)
    dataB = apply_cut( hdu.data, hdu.header['DATASECB'] )[0].astype(float)
    return ( dataB, biasB, dataA, biasA )

def correctDistSub(dist, bias1, empty, biasEmpty ):
    biasM = np.median(bias1)
    bias = subtractEmpty( bias1, biasEmpty )
    dataM = np.median(dist)
    data = subtractEmpty( dist, empty )
    
    return subtractBias( data + dataM, bias + biasM )

def addmin( dist ): return dist - np.min(dist)*1.1

def testing( hdu, input_file, tag = '' ): 
    print 'testing', hdu.header['OHDU'], hdu.header['CCDNUM']
    
    ohdu = 2
    def gethdu(ohdu_):
        parts_, processed_file_ = get_all_parts(input_file)
        merge_hdu_ = merge_parts(parts_, ohdu_)
        processed_ = astropy.io.fits.open(processed_file_)[ohdu]
        return merge_hdu_, processed_
    
    merge_hdu, processed = gethdu(ohdu)
    
    mdataB = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
    mdataA = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
    mbiasB = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECB'] )[0].astype(float)
    mbiasA = apply_cut( merge_hdu.data, merge_hdu.header['BIASSECA'] )[0].astype(float)

    biasA = apply_cut( hdu.data, hdu.header['BIASSECA'] )[0].astype(float)
    biasB = apply_cut( hdu.data, hdu.header['BIASSECB'] )[0].astype(float)
    dataA = apply_cut( hdu.data, hdu.header['DATASECA'] )[0].astype(float)
    dataB = apply_cut( hdu.data, hdu.header['DATASECB'] )[0].astype(float)

    mdata = correctDist( mdataB, mbiasB, mdataA, mbiasA )
    mdataBiasC = correctDistBiasC( mdataB, mbiasB, mdataA, mbiasA )
    mdataSub = correctDistSub( mdataB, mbiasB, mdataA, mbiasA )

    def plotHists( data, name, bins ):
        print name
        fig = plt.figure()
        ax = fig.add_subplot(221)
        h, b, _ = ax.hist( data, bins = bins, histtype='step' )
        h = np.array(h)
        x = (b[1:]+b[:-1])*.5
        high = h > max(h)/10.
        #print x[high]
        xmin, xmax = x[ high ][0], x[ high ][-1]
        #print xmin, xmax, max(h)
        ax.axhline( max(h)/10. )
        #ax.axvline( xmin, color='r' )
        #ax.axvline( xmax, color='r' )
        par = computePeakParams( data )
        print par
        g = scipy.stats.norm.pdf(x,*par)
        ax.plot( x, g/max(g)*max(h), 'r' )
        ax.set_ylabel(name)
        ax.set_xlabel('adu')
        ax.set_yscale('log')
        ax.set_ylim(.1, ax.get_ylim()[1])
        fig.savefig('%s%s.pdf'%(name,tag), bbox_inches='tight')
    
    plotHists(biasA.flatten(),'biasA', 50)
    plotHists(biasB.flatten(),'biasB', 5000)
    plotHists(dataA.flatten(),'dataA', 1000)
    plotHists(dataB.flatten(),'dataB', 5000)
    exit(1)
    plotHists(mbiasA.flatten(),'mbiasA', 50)
    plotHists(mbiasB.flatten(),'mbiasB', 5000)
    plotHists(mdataA.flatten(),'mdataA', 1000)
    plotHists(mdataB.flatten(),'mdataB', 5000)
    dataproc = apply_cut( processed.data, fix_processed_header( merge_hdu ) )[0].astype(float)
    print np.max(dataproc)
    dataproc[ dataproc == 1.e10 ] = 0
    print np.max(dataproc)
    plotHists(dataproc.flatten(),'processed', 5000)

    plotHists(mdata.flatten(),'mdata', 5000)
    plotHists(mdataBiasC.flatten(),'mdataBiasC', 5000)
    plotHists(mdataSub.flatten(),'mdataSub', 5000)
    #subdivideSection(dataB)
    #exit()
    #fig3 = plt.figure()
    #dataB2 = dataB - np.median(dataB)
    #dB = dataB2.flatten()
    
    #histdataB, binsdataB = np.histogram( dB, bins=10000)
    #dataM = correctDistMedian( dataB, biasB, dataA, biasA )
    #data = correctDist( dataB, biasB, dataA, biasA )
    #dataBiasC = correctDistBiasC( dataB, biasB, dataA, biasA )
    #dataSub = correctDistSub( dataB, biasB, dataA, biasA )
    
    #histdataM, tmp = np.histogram( dataM.flatten(), bins=binsdataB )
    #histdata, tmp = np.histogram( data.flatten(), bins=binsdataB )
    #histdataBiasC, tmp = np.histogram( dataBiasC.flatten(), bins=binsdataB )
    #histdataSub, tmp = np.histogram( dataSub.flatten(), bins=binsdataB )
    
    #dbins = binsdataB[1]-binsdataB[0]
    #xbins = (binsdataB[1:] + binsdataB[:-1])*.5
    #N = float(np.sum(histdataB)*dbins)

    #plots = [
        #lambda ax: ax.step( xbins, histdataB, label = 'dataB' ),
        #lambda ax: ax.step( xbins, histdataM, label = 'dataM', lw=2 ),
        #lambda ax: ax.step( xbins, histdata, label = 'data', lw=2 ),
        #lambda ax: ax.step( xbins, histdataBiasC, label = 'dataBC', lw=2 ),
        #lambda ax: ax.step( xbins, histdataSub, label = 'dataSub', lw=2 ),
        
        #lambda ax: ax.step( xbins, N*scipy.stats.norm.pdf( xbins, *computePeakParams( dataM, f=100, d=100 )), label = 'dataM', lw= 2 ),
        #lambda ax: ax.step( xbins, N*scipy.stats.norm.pdf( xbins, *computePeakParams( data, f=100, d=100 )), label = 'dataP', lw= 2 ),
        #lambda ax: ax.step( xbins, N*scipy.stats.norm.pdf( xbins, *computePeakParams( dataBiasC, f=100, d=100 )), label = 'dataBC', lw= 2 ),
        #lambda ax: ax.step( xbins, N*scipy.stats.norm.pdf( xbins, *computePeakParams( dataSub, f=100, d=100 )), label = 'dataSub', lw= 2 ),
        #]
    
    #ax = fig3.add_subplot(1, 1, 1)
    #ax.set_yscale('log')
    #ax.set_ylim([.1,max(histdataB)*1.5])
    ##ax.set_xlim([min(binsdataB),params[0]+(params[0]-min(binsdataB))])
    #for i, plot in enumerate(plots):
        #plot( ax )
    #plt.legend()
    #plt.show()
    #print 'saving fig3'
    #fig3.savefig('fitting%s.png'%tag)
    #print 'done'
    #exit()
    def plotCCD( data, name, plain=False ):
        print name
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if plain:
            ax.imshow((data-np.median(data)).T)
        else:
            par = computePeakParams( data )
            max_ = par[0] + 2*par[1]
            min_ = par[0] - 2*par[1]
            print max_
            data[ data > max_ ] = max_
            data[ data < min_ ] = min_
            ax.imshow(data.T - par[0])
        ax.set_title(name)
        ax.set_xlabel('horizontal')
        ax.set_ylabel('vertical')
        fig.savefig( '%s%s.pdf'%(name,tag), bbox_inches='tight')
    
    #dataM = correctDistMedian( mdataB, mbiasB, mdataA, mbiasA )
    
    #plotCCD( apply_cut( hdu.data, hdu.header['DATASEC'] )[0].astype(float), 'ccdpartall', plain = True )
    print merge_hdu.data.shape
    print merge_hdu.header['DATASECB']
    #plotCCD( mdataB, 'ccdmerged' )
    #plotCCD( mdataA, 'ccdmergedA' )

    #plotCCD( mdata, 'ccddata' )
    #plotCCD( mdataBiasC, 'ccdBiasC' )
    #plotCCD( mdataSub, 'ccdSub' )
    
    #print fix_processed_header( merge_hdu )
    #print np.sum(processed.data)
    #plotCCD( dataproc, 'ccdprocessed' )
    def plotRows( data, name ):
        print name
        fig = plt.figure()
        ax = fig.add_subplot(221)
        par = computePeakParams( data )
        max_ = par[0] + 2*par[1]
        min_ = par[0] - 2*par[1]
        print max_
        data[ data > max_ ] = par[0]
        data[ data < min_ ] = par[0]
        data -= par[0]
        #x = np.array(range(0, len(data[:,0])) )[:,np.newaxis] + np.zeros_like(data[0,:])[np.newaxis,:]
        x = np.array(range(0, len(data[0,:])) )
        print x.shape, data.shape, np.mean(data, axis=0).shape
        ax.plot( x, np.mean(data, axis=0), '.', alpha=.1 )
        ax.set_title(name)
        ax.set_xlabel('vertical')
        ax.set_ylabel('adu')
        fig.savefig( '%s%s.pdf'%(name,tag), bbox_inches='tight')
    
    plotRows( mdataB, 'rowsdataB' )
    plotRows( mdata, 'rowsdata' )
    plotRows( mdataBiasC, 'rowsBiasC' )
    plotRows( mdataSub, 'rowsSub' )
    plotRows( dataproc, 'rowsProcessed' )
    exit(1)
    
    
    
    print 'minB', np.min(dataB)
    hist, bins = np.histogram( dataB, bins=1000)
    xbins = (bins[:-1] + bins[1:])*.5

    print 'finsished'
    fig2 = plt.figure(figsize=(10,20))
    
        
    
    Nplots = 10
    current_plot = 1
    tmp = dataB2
    fig2.add_subplot(Nplots, 1, current_plot).imshow( np.log(tmp.T) )
    current_plot += 1
    
    print len(tmp[1,:]), len(tmp[:,1]), len(np.mean(tmp,axis=0)), len(np.mean(tmp,axis=1))
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(1,len(tmp[1,:])+1), np.median(dataB2,axis=0)+10, 'b-' )
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(len(tmp[1,:]),len(tmp[1,:])+len(tmp[:,1])), np.median(dataB2,axis=1)+10, 'r-' )
    current_plot += 1    

    tmp = dataM
    fig2.add_subplot(Nplots, 1, current_plot).imshow( np.log(tmp.T) )
    current_plot += 1
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(1,len(tmp[1,:])+1), np.median(dataM,axis=0)+10, 'b-' )
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(len(tmp[1,:]),len(tmp[1,:])+len(tmp[:,1])), np.median(dataM,axis=1)+10, 'r-' )
    current_plot += 1

    tmp = dataBiasC
    fig2.add_subplot(Nplots, 1, current_plot).imshow( np.log(tmp.T) )
    current_plot += 1
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(1,len(tmp[1,:])+1), np.median(dataBiasC,axis=0)+10, 'b-' )
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(len(tmp[1,:]),len(tmp[1,:])+len(tmp[:,1])), np.median(dataBiasC,axis=1)+10, 'r-' )
    current_plot += 1

    tmp = dataSub
    fig2.add_subplot(Nplots, 1, current_plot).imshow( np.log(tmp.T) )
    current_plot += 1
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(1,len(tmp[1,:])+1), np.median(dataSub,axis=0)+10, 'b-' )
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(len(tmp[1,:]),len(tmp[1,:])+len(tmp[:,1])), np.median(dataSub,axis=1)+10, 'r-' )
    current_plot += 1

    tmp = data
    fig2.add_subplot(Nplots, 1, current_plot).imshow( np.log(tmp.T) )
    current_plot += 1
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(1,len(tmp[1,:])+1), np.median(data,axis=0)+10, 'b-' )
    fig2.add_subplot(Nplots, 1, current_plot).semilogy( range(len(tmp[1,:]),len(tmp[1,:])+len(tmp[:,1])), np.median(data,axis=1)+10, 'r-' )
    current_plot += 1
    
    fig2.savefig('ccd%s.pdf'%tag)
    return 


def temperaturePlotter( start_time = None ):
    temp_file2 = "/share/storage1/connie/data/liveMonitoring/logs/logs_MCP9700/log_MCP9700.txt"
    data = np.genfromtxt(temp_file2, delimiter=' ')
    xx = map( datetime.datetime.fromtimestamp, data[:,0])
    yy = data[:,1]
    yy2 = data[:,2]
    #print data
    #exit()
    
    directory = "/share/storage1/connie/data/liveMonitoring/logs/logs_monidae/log_%s*.txt"%(start_time)
    list_of_files = glob.glob(directory)
    x = []
    y = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    y6 = []
    y7 = []
    y8 = []
    y9 = []
    
    show = True
    for input_file in list_of_files:
        tmp = os.path.basename(input_file).split('_')[1].split('.')[0]
        #print tmp
        #file_time = datetime.datetime.strptime( tmp, '%Y%m%d%H%M' )
        #print file_time.strftime('%d%b%y_%H:%M')
        txt = open(input_file).read().replace(' ', '\t')
        data = np.genfromtxt( StringIO.StringIO(txt), dtype=None, delimiter='\t', names=True )#, names = True)
        #print data.dtype
        try:
            file_time = data['Time'][0]
        except:
            continue
        
        if show:
            print data.dtype
            print data[0]
            show = False
        #print datetime.datetime.fromtimestamp( data['Time'][0] ).strftime('%d%b%y_%H:%M')
        #if start_time:
            #if start_time > file_time:
                #print "previous to the start date"
                #print input_file
                #continue
        x += map( datetime.datetime.fromtimestamp, data['Time'].tolist() )
        y += data['Temp'].tolist()
        #y2 += data['relay'].tolist()
        #y3 += data['setTEMP'].tolist()
        y4 += data['vpTempBrng'].tolist()
        y5 += data['vpTempMtr'].tolist()
        y6 += data['vpTempBttm'].tolist()
        y7 += data['vpTempElec'].tolist()
        
        y8 += data['Pressure'].tolist()
        y9 += data['htrPow'].tolist()
        
        #print "processing", input_file
    #datetime.datetime.fromtimestamp( hdu.header['EXPSTART'] ).strftime('%d%b%y_%H:%M')
    fig = plt.figure()
    ax = fig.add_subplot(311)
    #ax.plot( x, y2, label='relay' )
    #ax.plot( x, y3, label='setTEMP' )
    #ax.plot( x, y4, label='vpTempBrng' )
    ax.plot( x, y5, label='vpTempMtr' )
    #ax.plot( x, y6, label='vpTempBttm' )
    #ax.plot( x, y7, label='vpTempElec' )
    
    ax.plot( x, y9, label='htrPow' )
    ax.plot( xx, yy, label='MCP9700[1]' )
    
    plt.legend()
    ax2 = fig.add_subplot(312, sharex = ax)
    ax2.plot( x, y, label='Temp' )
    #ax2.plot( xx, yy2, label='MCP9700[2]' )
    plt.legend()
    
    #ax3 = fig.add_subplot(313, sharex = ax)
    #ax3.semilogy( x, y8, label='Pressure' )
    #plt.legend()
    #vpTempMtr vpTempBttm vpTempElec
    ax.format_xdata = matplotlib.dates.DateFormatter('%m%d')
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%d/%m'))
    #ax.xaxis.set_minor_formatter(matplotlib.dates.DateFormatter('%d'))
    ax.grid(True)
    ax2.grid(True)
    ax2.set_xlabel('t')
    fig.autofmt_xdate()
    plt.savefig( 'plots/temp.png' )
    plt.close(fig)
    return

def get_hdu_list( file_name ):
    actual_name = file_name
    if os.path.exists(file_name):
        pass
    elif os.path.exists(file_name + '.fz'):
        actual_name += '.fz'
    else:
        print 'error', file_name
    return astropy.io.fits.open(actual_name)

def get_all_parts( file_name ):
    actual_name = file_name
    if os.path.exists(file_name):
        pass
    elif os.path.exists(file_name + '.fz'):
        actual_name += '.fz'
    else:
        print 'error', file_name
    
    file_names = re.sub('_p(\d)\.fits', '_p*.fits', actual_name)
    actual_files = [ f for f in glob.glob(file_names) if os.path.exists(f) ]
    processed = glob.glob('/share/storage1/connie/data_analysis/processed_data/runs/*/*/scn/merged/scn_mbs_osi_%s*' % actual_name.split('/')[-1].split('_p')[0])
    
    if processed == []:
        print 'no processed'
        processed = None
    elif not os.path.exists( processed[0] ):
        print processed, 'does not exists'
        processed = None
    else:
        processed = processed[0]
    return actual_files, processed
    
def merge_parts( parts, index ):
    hdu1 = astropy.io.fits.open(parts[0])[index]
    hdu2 = astropy.io.fits.open(parts[1])[index]
    hdu3 = astropy.io.fits.open(parts[2])[index]
    hdu4 = astropy.io.fits.open(parts[3])[index]
    #print hdu1.data.shape
    merge_hdu = hdu1
    merge_hdu.data = np.concatenate( [hdu1.data, hdu2.data, hdu3.data, hdu4.data], axis = 0 )
    for tag in [ 'BIASSECA', 'BIASSECB', 'DATASECA', 'DATASECB', 'DATASEC' ]:
        merge_hdu.header[tag] = merge_hdu.header[tag].replace('1055', '4220')
    print merge_hdu.data.shape
    return merge_hdu

def fix_processed_header( hdu_part ):
    merge_hdu_header = hdu_part.header['DATASECB'].replace('1055', '4220')
    cut = merge_hdu_header
    return cut.replace('[9:4120', '[1:4112')


def plotCommand( output_file, startID = 0, tag='' ):
    line = {}
    #keys = [ 'start', 'tempMax', 'tempMin', 'std_biasA', 'std_gauss_biasA', 'lambda_poisson_biasA', 'std_biasB', 'std_gauss_biasB', 'lambda_poisson_biasB', 'std_dataA', 'std_dataB', ]
    data = json.load( open( output_file, 'r' ) )
    
    #print 
    print len(data.keys())
    
    #f = open('data3.csv', 'w')
    #f.write('# filename, ohdu, start, std_biasA, std_biasD, std_dataA, std_dataB\n')
    #for file_name, entry in data.items():
        #for ohdu, values in entry.items():
            #if not 'ohdu' in ohdu: continue
        
            #if not ohdu in line: line[ohdu] = dict( [ (key,[]) for key in keys ] )
            #print file_name, ohdu, values['std_biasA'], np.average(line[ohdu]['std_biasA']), 5*np.std(line[ohdu]['std_biasA'])
            #hdu_list = None
            
            #actual_name = file_name
            #if os.path.exists(file_name):
                #pass
            #elif os.path.exists(file_name + '.fz'):
                #actual_name += '.fz'
            #else:
                #print 'error', file_name
            
            #parts, processed = get_all_parts( file_name )
            #index = int(ohdu)-1
            #hdu_list = astropy.io.fits.open(actual_name)
            #print hdu_list[index].header['OHDU'], ohdu
            ##ohdu
            #hdu = hdu_list[index]
            
            ##fig = plt.figure()
            ##ax = fig.add_subplot(111)
            
            #dataB = apply_cut( hdu.data, hdu.header['DATASECB'] )[0].astype(float)
            #dataA = apply_cut( hdu.data, hdu.header['DATASECA'] )[0].astype(float)
            ##n, bins, tmp = ax.hist( dataB.flatten(), bins = 5000, histtype = 'step' )
            #print computePeakParams( dataA )
            #print computePeakParams( dataB )
            #print computePeakParams( correctDist( *extract_hdu(hdu)) )
            #print computePeakParams( correctDistBiasC( *extract_hdu(hdu)) )
            #print computePeakParams( correctDistSub( *extract_hdu(hdu)) )
            
            
            
            #hdu1 = astropy.io.fits.open(parts[0])[index]
            #hdu2 = astropy.io.fits.open(parts[1])[index]
            #hdu3 = astropy.io.fits.open(parts[2])[index]
            #hdu4 = astropy.io.fits.open(parts[3])[index]
            ##print hdu1.data.shape
            #merge_hdu = hdu1
            #merge_hdu.data = np.concatenate( [hdu1.data, hdu2.data, hdu3.data, hdu4.data], axis = 0 )
            #for tag in [ 'BIASSECA', 'BIASSECB', 'DATASECA', 'DATASECB' ]:
                #merge_hdu.header[tag] = merge_hdu.header[tag].replace('1055', '4220')
                ##print merge_hdu.header[tag]
            #print merge_hdu.data.shape
            #dataB2 = apply_cut( merge_hdu.data, merge_hdu.header['DATASECB'] )[0].astype(float)
            #dataA2 = apply_cut( merge_hdu.data, merge_hdu.header['DATASECA'] )[0].astype(float)
            ##ax.hist( dataB2.flatten(), bins = bins, histtype = 'step' )
            #print computePeakParams( dataA2 )
            #print computePeakParams( dataB2 )
            #print computePeakParams( correctDist(  *extract_hdu(merge_hdu)) )
            #print computePeakParams( correctDistBiasC(  *extract_hdu(merge_hdu)) )
            #print computePeakParams( correctDistSub(  *extract_hdu(merge_hdu)) )
            
            #proc = astropy.io.fits.open(processed)[index-1]
            #print proc.header['OHDU'], ohdu
            #print proc.data.shape
            #cut = merge_hdu.header['DATASECB']
            #cut = cut.replace('[9:4120', '[1:4112')
            #print cut
            #vmax = np.max(dataB)
            #dataB = apply_cut( proc.data, cut )[0].astype(float)
            #dataB[ dataB > vmax/2. ] = 0
            ##ax.hist( dataB.flatten(), bins = 5000, histtype = 'step' )
            ##plt.show()
            #print computePeakParams( dataB )
            
            #for key in line[ohdu].keys():
                #try:
                    #line[ohdu][key] += [ values[key] ]
                #except:
                    #print values
            #f.write( ' '.join(map(str, [file_name, ohdu] + [values[key] for key in line[ohdu].keys()]) ) + '\n' )
    #f.close()
    
    #print line
    #TOOLS = "hover,crosshair,pan,wheel_zoom,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select"
    #p = bokeh.plotting.figure(tools=TOOLS)

    ohdu_list = map(lambda i: 'ohdu%d'%i, range(2,15))
    parts_peaks_dataA = {}
    parts_peaks_dataB = {}
    parts_peaks_dist = {}
    
    merged_peaks_dataA = {}
    merged_peaks_dataB = {}
    merged_peaks_dist = {}

    proc_peaks_dist = {}
    expstart = []
    runid = []
    tempMax = []
    tempMin = []
    for key, entry in data.items():
        if entry[ 'runid' ] < startID: continue
        skip = False
        for ohdu in ohdu_list:
            if not ohdu in entry: continue
            #print [ key for key in entry.keys() if 'ohdu' in key]
            if not ohdu in parts_peaks_dataA: parts_peaks_dataA[ ohdu ] = []
            if not ohdu in parts_peaks_dataB: parts_peaks_dataB[ ohdu ] = []
            if not ohdu in parts_peaks_dist: parts_peaks_dist[ ohdu ] = []
            if not ohdu in merged_peaks_dataA: merged_peaks_dataA[ ohdu ] = []
            if not ohdu in merged_peaks_dataB: merged_peaks_dataB[ ohdu ] = []
            if not ohdu in merged_peaks_dist: merged_peaks_dist[ ohdu ] = []
            if not ohdu in proc_peaks_dist : proc_peaks_dist[ohdu] = []
            
            #print entry[ ohdu ]['parts']['peaks_dataA']
            #print [ val[1] for val in entry[ ohdu ]['parts']['peaks_dataA'] ]
            vals = [ val[1] for val in entry[ ohdu ]['parts']['peaks_dataA'] ]
            if len(vals) != 4: 
                print 'wrong number of parts'
                print ohdu, key
                skip = True
                break
                #vals += [0]
            
            parts_peaks_dataA[ ohdu ] += [ vals ]
            vals = [ val[1] for val in entry[ ohdu ]['parts']['peaks_dataB'] ]
            parts_peaks_dataB[ ohdu ] += [ vals ]
            
            vals = [ val[1] for val in entry[ ohdu ]['parts']['peaks_dist'] ]
            parts_peaks_dist[ ohdu ] += [ vals ]
            
            if entry[ ohdu ]['merged']['peaks_dataA'] != None:
                merged_peaks_dataA[ ohdu ] += [ entry[ ohdu ]['merged']['peaks_dataA'][1] ]
            else: 
                merged_peaks_dataA[ ohdu ] += [ 0 ]
            if entry[ ohdu ]['merged']['peaks_dataB'] != None:
                merged_peaks_dataB[ ohdu ] += [ entry[ ohdu ]['merged']['peaks_dataB'][1] ]
            else:
                merged_peaks_dataB[ ohdu ] += [ 0 ]
            if entry[ ohdu ]['merged']['peaks_dist'] != None:
                merged_peaks_dist[ ohdu ] += [ entry[ ohdu ]['merged']['peaks_dist'][1] ]
            else:
                merged_peaks_dist[ ohdu ] += [ 0 ]
            
            if 'processed' in entry[ ohdu ]:
                if 'peaks' in entry[ ohdu ]['processed']:
                    proc_peaks_dist[ ohdu ] += [ entry[ ohdu ]['processed']['peaks'][1] ]
                else: 
                    proc_peaks_dist[ ohdu ] += [ 0 ]
            else:
                proc_peaks_dist[ ohdu ] += [ 0 ]
        if skip: continue
        expstart += [ entry[ 'expstart' ] ]
        runid += [ entry[ 'runid' ] ]
        tempMax += [ entry[ 'tempMax' ] ]
        tempMin += [ entry[ 'tempMin' ] ]
    
    
    #print 'parts_peaks_dataA', parts_peaks_dataA['ohdu2']
    #print 'expstart', expstart
    #expstart = np.array(expstart)
        
    #ind = np.array(expstart).argsort()
    
    colors = ['r', 'g', 'b' ]
    axes = {}
    ykeys = ['peaks_dataA', 'peaks_dataB']
    typekeys = ['parts', 'merged', 'processed']
    print ohdu_list
    for i, ohdu in enumerate(ohdu_list):
        fig = plt.figure(figsize=(10,10))
        #M = len(line[ohdu].keys())-1
        M = 7
        #N = int(np.sqrt(M))+1
        #line[ohdu]['start'] = map( datetime.datetime.fromtimestamp, line[ohdu]['start'] )
        
        y = np.array( parts_peaks_dataA[ohdu] )
        s = 5
        print 'x.shape', y.shape, len(expstart), len(runid)
        bbox_to_anchor = (1,.5)
        axes = fig.add_subplot(M,1,1)
        axes.scatter( runid, [ p[0] for p in parts_peaks_dataA[ohdu]], color='r', s=s, label = 'parts dataA' )
        axes.scatter( runid, [ p[1] for p in parts_peaks_dataA[ohdu]], color='g', s=s, label = 'parts dataA' )
        axes.scatter( runid, [ p[2] for p in parts_peaks_dataA[ohdu]], color='b', s=s, label = 'parts dataA' )
        axes.scatter( runid, [ p[3] for p in parts_peaks_dataA[ohdu]], color='k', s=s, label = 'parts dataA' )
        axes.grid(True)
        axes.set_xlabel('runid')
        axes.set_ylabel('adu')
        leg = plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor, numpoints=1)

        axes = fig.add_subplot(M,1,2)
        axes.scatter( runid, [ p[0] for p in parts_peaks_dataB[ohdu]], color='r', s=s, label = 'parts dataB' )
        axes.scatter( runid, [ p[1] for p in parts_peaks_dataB[ohdu]], color='g', s=s, label = 'parts dataB' )
        axes.scatter( runid, [ p[2] for p in parts_peaks_dataB[ohdu]], color='b', s=s, label = 'parts dataB' )
        axes.scatter( runid, [ p[3] for p in parts_peaks_dataB[ohdu]], color='k', s=s, label = 'parts dataB' )
        axes.set_xlabel('runid')
        axes.set_ylabel('adu')
        #plt.gca().set_yscale("log")
        axes.grid(True)
        leg = plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor, numpoints=1)

        axes = fig.add_subplot(M,1,3)
        axes.scatter( runid, [ p[0] for p in parts_peaks_dist[ohdu]], color='r', s=s, label = 'parts dist' )
        axes.scatter( runid, [ p[1] for p in parts_peaks_dist[ohdu]], color='g', s=s, label = 'parts dist' )
        axes.scatter( runid, [ p[2] for p in parts_peaks_dist[ohdu]], color='b', s=s, label = 'parts dist' )
        axes.scatter( runid, [ p[3] for p in parts_peaks_dist[ohdu]], color='k', s=s, label = 'parts dist' )
        axes.set_xlabel('runid')
        axes.set_ylabel('adu')
        #plt.gca().set_yscale("log")
        axes.grid(True)
        #plt.gca().set_ylim([0,20])
        leg = plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor, numpoints=1)
        
        axes = fig.add_subplot(M,1,4)
        axes.scatter( runid, merged_peaks_dataA[ohdu], color='r', s=s, label = 'merged dataA' )
        axes.scatter( runid, merged_peaks_dataB[ohdu], color='g', s=s, label = 'merged dataB' )
        axes.scatter( runid, merged_peaks_dist[ohdu], color='b', s=s, label = 'merged dist' )
        axes.scatter( runid, proc_peaks_dist[ohdu], color='k', s=s, label = 'processed' )
        axes.set_xlabel('runid')
        axes.set_ylabel('adu')
        #plt.gca().set_yscale("log")
        axes.grid(True)
        leg = plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor, numpoints=1)

        axes = fig.add_subplot(M,1,5)
        axes.scatter( runid, tempMax, color='r', s=s, label = 'tempMax' )
        axes.scatter( runid, tempMin, color='g', s=s, label = 'tempMin' )
        axes.set_xlabel('runid')
        axes.set_ylabel('tempMax')
        axes.grid(True)
        leg = plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor, numpoints=1)

        axes = fig.add_subplot(M,1,6)
        #axes.scatter( tempMax, [ p[0] for p in parts_peaks_dataA[ohdu]], color='r', s=s, label = 'parts dataA' )
        #axes.scatter( tempMax, [ p[0] for p in parts_peaks_dataB[ohdu]], color='g', s=s, label = 'parts dataB' )
        #axes.scatter( tempMax, [ p[0] for p in parts_peaks_dist[ohdu]], color='b', s=s, label = 'parts dist' )
        axes.scatter( tempMax, merged_peaks_dataA[ohdu], color='r', s=s, label = 'merged dataA' )
        axes.scatter( tempMax, merged_peaks_dataB[ohdu], color='g', s=s, label = 'merged dataB' )
        axes.scatter( tempMax, merged_peaks_dist[ohdu], color='b', s=s, label = 'merged dist' )
        axes.scatter( tempMax, proc_peaks_dist[ohdu], color='k', s=s, label = 'processed' )
        axes.set_xlabel('tempMax')
        axes.set_ylabel('adu')
        #plt.gca().set_yscale("log")
        axes.grid(True)
        leg = plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor, numpoints=1)
        
        #axes.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%d/%m'))
        #leg = plt.legend(loc='center left', bbox_to_anchor=(1,1), numpoints=1)
        #plt.gca().set_yscale("log")
        #plt.gca().set_ylim([0,50])
        plt.savefig( 'test%s%s.pdf'%(ohdu,tag), bbox_extra_artists=(leg,), bbox_inches='tight' )
        plt.close(fig)
        
    print "save done"

if __name__ == "__main__":
    #input_file = "/share/storage1/connie/data/2016/scienceTest/15Aug2016/00/runID_00_00050_Int-400_Exp-10800_18Aug16_13:48_to_18Aug16_16:52_p1.fits"
    output_file = 'data4.json'
    
    if 'plot' in sys.argv:
        plotCommand( output_file, 470, tag='zoom' )
        plotCommand( output_file, 0, tag='' )
    elif 'temp' in sys.argv:
        temperaturePlotter('201607')
    elif 'conv' in sys.argv:
        test_convolution()
    elif 'file' in sys.argv:
        input_file = sys.argv[-1]
        if not '/share/' in input_file:
            input_file = glob.glob('/share/storage1/connie/data/runs/*/' + input_file + '*')
            if len(input_file) > 0: input_file = input_file[0]
        compute_values( input_file, dict(), complete=True )
    elif 'table' in sys.argv[-1]:
        print 'table option'
        if not os.path.exists(output_file):
            print 'no file', output_file
            print 'exiting...'
            exit(0)
        data = json.load( open( output_file, 'r' ) )
        #print 'keys', data.keys()
        output = ''
        header = '# ' + ', '.join(['date', 'expstart', 'ohdu', 'mean.dataB', 'std.dataB', 'mean.biasB', 'std.biasB'])
        output += header + '\n'
        x = []
        y = []
        for key, entry in data.iteritems():
            if entry['expstart'] < 0: continue
            for skey, sentry in entry.iteritems():
                if 'ohdu' in skey:
                    if 'peaks_biasB' in sentry['merged'] and 'peaks_dataB' in sentry['merged']:
                        #print sentry['merged']['peaks_biasB'][0], sentry['merged']['peaks_biasB'][1]
                        #x += [ entry['expstart'] ]
                        try:
                            y += [ sentry['merged']['peaks_dataB'][0] - sentry['merged']['peaks_biasB'][0] ]
                            x += [ datetime.datetime.fromtimestamp( entry['expstart'] ) ]
                            row = ', '.join([
                                entry['date'],
                                str(entry['expstart']),
                                skey, 
                                str(sentry['merged']['peaks_dataB'][0]),
                                str(sentry['merged']['peaks_dataB'][1]),
                                str(sentry['merged']['peaks_biasB'][0]),
                                str(sentry['merged']['peaks_biasB'][1]),
                                ])
                            output += row + '\n'
                        except:
                            print entry['date']
                            print sentry['merged']
                        #print entry['expstart'], skey, sentry['merged']['peaks_dataB'][0], sentry['merged']['peaks_dataB'][1], sentry['merged']['peaks_biasB'][0], sentry['merged']['peaks_biasB'][1]
        open( 'data4.csv', 'w' ).write( output )
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy( x, y, '.' )
        #ax.format_xdata = matplotlib.dates.DateFormatter('%m%d')
        #ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%d/%m'))
        fig.autofmt_xdate()
        plt.savefig('data4.dc.png')
        print 'expstart', 'ohdu', 'mean.dataB', 'std.dataB', 'mean.biasB', 'std.biasB'
    else:
        if not os.path.exists(output_file):
            open( output_file, 'w' ).write('{}')
        data = json.load( open( output_file, 'r' ) )
        for key in data.keys():
            print key
        print len(data.keys()), 'files done'
        #list_of_files = glob.glob("/share/storage1/connie/data/2016/scienceTest/15Aug2016/*/runID_*Int-400_Exp-10800*_p*.fits*") + glob.glob("/share/storage1/connie/data/runs/*/runID_*Int-400_Exp-10800*_p*.fits*")
        list_of_files = glob.glob("/share/storage1/connie/data/runs/*/runID_*Int-*_Exp-*_p1.fits*") #only take _p1 files since comput_values looks for the other part files and processed files
            #+ glob.glob("/share/storage1/connie/data_analysis/processed_data/runs/*/*/*/*/*.fits")
        #print list_of_files
        #previous_file = None
        for input_file in list_of_files[::-1]:
            #if input_file.rstrip('.fz') in data:
                #print 'already done', input_file
                #continue
            print input_file
            entry = compute_values( input_file, data )
            if entry is None: 
                print 'entry is None'
                print 'not writing'
                continue
            data.update( entry )
            print 'writing...'
            json.dump( data, open( output_file, 'w' ), indent = True, sort_keys = True )
            print '... done'
            #previous_file = input_file.split('_p')
        
        
