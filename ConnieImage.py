# coding: utf-8

import astropy
#import astropy.io
#import astropy.io.fits

import math
import numpy as np
from numpy.lib.recfunctions import append_fields, stack_arrays
import scipy
import scipy.stats
import scipy.signal
import scipy.special
from scipy.misc import factorial
import scipy.optimize

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import patches, colors

import os
import time
import random

import root_numpy

from ConnieDataPath import ConnieDataPath as ConniePaths
import Statistics as stats

class Image:
    __row__ = 0
    __col__ = 1
    
    def __init__( self, matrix ):
        if not len(matrix.shape) == 2:
            raise Exception( 'shape is not 2d' )
        self.matrix = matrix
    
    def get_nrows(self):
        return self.matrix.shape[Image.__row__]

    def get_ncols(self):
        return self.matrix.shape[Image.__col__]
    
    def plot_to_axis( self, ax, **kwargs ):
        ax.hist( self.matrix.flatten(), kwargs )
        return ax

    #def get_mean( self ):
        
class CCDImage:
    def __init__( self ):
        self.data_section = None
        self.overscan_section = None
        self.voverscan_section = None
        self.trim_section = None
    
    def add_section( self, section, _type_ ):
        if not isinstance(section, Image):
            section = Image(section)
        self.__dict__[_type_] = section
        return self
        
    def add_overscan(self, section): return self.add_section(section, 'overscan_section')
    def add_data(self, section): return self.add_section(section, 'data_section')
    def add_voverscan(self, section): return self.add_section(section, 'voverscan_section')
    
def readImage_raw( runID, ohdu ):
    print 'readImage_raw, runID=', runID, 'ohdu=', ohdu, 
    files = ConniePaths.runIDPath( int(runID) )
    print 'path=', files[-1]
    return merge_images( files, ohdu )

def readImage_MB( runID, ohdu ):
    files = ConniePaths.masterBiasPath( int(runID) )
    return merge_images( files, ohdu )

def readImage_processed( runID, ohdu, image='scn', debug=False ):
    #files = ConniePaths.runIDPathProcessed( int(runID), image=image )
    if debug: print '(readImage_processed)runID=%s'%runID
    print 'readImage_processed, runID=', runID, 'ohdu=', ohdu,
    files = ConniePaths.runIDPathProcessed_v3( int(runID), image=image, debug=debug )
    print 'path=', files[-1]
    return merge_images( files, ohdu, debug=debug)

def merge_images( files, ohdu, debug=False ):
    if debug: print '(merge_images)files=%s,ohdu=%s' % (files, ohdu)
    l = []
    for file in files[::-1]:
        fits = astropy.io.fits.open(file)
        entry = None
        for fit in fits:
            if debug: print '(merge_images)HEADER=%s' % (fit.header['OHDU'])
            if int(fit.header['OHDU']) == ohdu: entry = fit.data
        if entry is None: 
            raise Exception('ohdu not found')
        l.append( entry )
    #ohdus = lambda fits: fits[ohdu-1].data
    image = np.concatenate( l, axis=0 )
    #print '(merge_images)image.shape=', image.shape
    return image

class ArrayExtension:
    def __init__(self, array ):
        self.array = array.astype(float)
        self.shape = self.array.shape
        self.size = self.array.size
    
    def median(self, axis=None): return np.median( self.array, axis=axis )

    def MAD(self, axis=None): return stats.MAD( self.array, axis=axis )

    def IQR(self): return stats.IQR( self.array )

    def mean(self, axis=None): 
        return np.mean( self.array, axis=axis )

    def median_half(self, axis=None): return np.median( self.array[self.shape[0]/2:,:], axis = axis )

    def MAD_half(self, axis=None): return stats.MAD( self.array[self.shape[0]/2:,:], axis=axis )

    def std(self, axis=None): return np.std( self.array, axis=axis )

    def flatten(self): return self.array.flatten()

    def mean_medians( self ):
        return .5*( np.mean( self.median(axis=0) ) + np.mean( self.median(axis=1) ) )

    def mean_MADs( self ):
        return .5*( np.mean( self.MAD(axis=0) ) + np.mean( self.MAD(axis=1) ) )

class ImageBase( ArrayExtension ):
    def __init__( self, image, debug=False ):
        if debug: print '(ImageBase.__init__)'
        ArrayExtension.__init__(self, image.astype(float))
        self.image = self.array
        self.shape = self.image.shape
        self.size = self.image.size
        self.T = self.image.T
        self.width = self.shape[1]
        self.height = self.shape[0]

    def __call__( self, mask ): return ImageBase( self.image[mask] )

    def flatten( self ): return self.image.flatten()
    
    def h_medians( self ): return self.median( axis = 1 )
    def h2_medians( self ): return np.median( self.image[:,self.width/2:], axis = 1 )
    def v_medians( self ): return self.median( axis = 0 )

    def h_MADs( self ): return self.MAD( axis = 1 )
    def h2_MADs( self ): return stats.MAD( self.image[:,self.width/2:], axis = 1 )
    
    def biasImage( self ):
        return np.median( self.image, axis=1 )[:,None] + np.median( self.image, axis=0 )[None,:]

    def biasImageCorrected( self ):
        image = ImageBase( self.image - self.median() )
        image.image -= image.biasImage()
        print 'bias shape', image.shape
        return ImageBase( image.image - image.median() )
    
    def biasImage2( self ):
        return np.median( self.image[:,self.width/2:], axis=1 )[:,None] + np.median( self.image, axis=0 )[None,:]

    def biasImage2Corrected( self ):
        image = ImageBase( self.image - self.median() )
        image.image -= image.biasImage2()
        print 'bias2 shape', image.shape
        return ImageBase( image.image - image.median() )
    
    def negatives( self ):
        return hitSumm.Image( -self.image[ self.image<0 ] )
    
    def electron( self ):
        stdNoise = hitSumm.normFit( self.os().image.flatten(), loc=0 )[1]
        stdNoiseNeg = hitSumm.normFit( self.os().negatives().image.flatten(), loc=0 )[1]
        stdDCNoiseNeg0 = hitSumm.normFit( self.negatives().image.flatten(), loc=0 )[1]
        muDCNoiseNeg, stdDCNoiseNeg = hitSumm.normFit( self.negatives().image.flatten() )
        muDCNoise, stdDCNoise = hitSumm.normFit( self.image[ np.abs(self.image) < 3*stdNoise ] .flatten() )
        stdDCNoise0 = hitSumm.normFit( self.image[ self.image < 5*stdNoise ] .flatten(), loc=0 )
        gain = muDCNoise/3/(eV_e/1e3)
        print( '(', muDCNoise, stdDCNoise, stdDCNoise0, stdDCNoiseNeg0, stdNoise, stdNoiseNeg, ')', gain )
        return muDCNoise, stdDCNoise
    
    def extract(self, thr=15):
        return processImage3(self.image, threshold = thr)
    
    def convolve( self, function=lambda x: np.sum(x), footprint=None, size=(7,7), origin=0 ):
        #output = None
        #center = 13 #4
        center = 25
        def _(x):
            #return np.abs( np.sum(x-x[center])/np.sum( np.abs(x) ) )
            #return np.var(x)
            return abs( np.mean(x) - np.median(x) )
        function = _
        print '(convolve)ishape', self.image.shape
        output = scipy.ndimage.filters.generic_filter(self.image, function, size=size, footprint=footprint, output=None, mode='constant', cval=0.0, origin=origin, extra_arguments=(), extra_keywords=None)
        print '(convolve)oshape', output.shape
        return ImageBase( output )
        
    
    
    def add_image_to_axis( self, axes, eMin, eMax ):
        axes.imshow( self.image, vmin = eMin, vmax = eMax, interpolation='none', origin='lower' )
        axes.set_adjustable('box-forced')
        axes.set_xlabel('x[pix]')
        axes.set_ylabel('y[pix]')
        axes.set_frame_on(False)
        #axes.grid( True )

    def add_projection_to_axis( self, axes, bins, axis, align='horizontal' ):
        if axis == 1: 
            image = self.image
            length = self.shape[0]
            other_length = self.shape[1]
        elif axis == 0: 
            image = self.T
            length = self.shape[1]
            other_length = self.shape[0]
        
        axes.grid( True )
        axes.set_frame_on(False)
        x = np.arange(length)
        medians = self.median( axis = axis )
        one_MAD = medians + self.MAD( axis = axis )
        if align == 'horizontal':
            axes.plot( x, medians, color='r' )
            axes.plot( x, one_MAD, color='m' )
            #hist2d, xedges, yedges = np.histogram2d( np.repeat(x, other_length), image.flatten(), bins=[length,bins] )
            #axes.imshow( hist2d.T, interpolation='none', origin='lower' )
            axes.hist2d( np.repeat(x, other_length), (2*np.median(image)-image).flatten(), bins=[length,bins], norm = colors.LogNorm() )
            #axes.hist2d( np.repeat(x, other_length), np.abs(-image).flatten(), bins=[length,bins], norm = colors.LogNorm() )
            axes.set_ylabel('E[adu]')
            axes.set_xlabel('x[pix]')
            axes.xaxis.tick_top()
            axes.xaxis.set_label_position('top')
            
        elif align == 'vertical':
            axes.hist2d( image.flatten(), np.repeat(x, other_length), bins=[bins, length], norm = colors.LogNorm() )
            axes.plot( medians, x[::-1], color='r' )
            axes.plot( one_MAD, x[::-1], color='m' )
            axes.set_xlabel('E[adu]')
            axes.set_ylabel('y[pix]')
            axes.yaxis.tick_right()
            axes.yaxis.set_label_position('right')
        
        return axes
    
    def printFlat( self, fname ):
        fig = plt.figure()
        grid = plt.GridSpec(1, 4, hspace=0, wspace=0)
        main_ax = fig.add_subplot(grid[:,:-2])
        main_ax.grid(True)
        y_hist = fig.add_subplot(grid[:, -1], sharey=main_ax )

        _, bins, _ = y_hist.hist( self.image.flatten(), bins=150, histtype='step', orientation="horizontal" )
        print( 'min', np.min(self.image), np.max(self.image), np.std(self.image) )
        y_hist.set_xscale('log')
        fig.savefig(fname)
        fig.clf()
        plt.close()
        return self
        
    def printImage( self, fname, log=False, func=None, n=10 ):
        fig = plt.figure()
        yticks = [ 0, 1055, 2*1055, 3*1055, 4*1055 ]
        xticks = [ 0, self.image.shape[1]-150, self.image.shape[1] ]
        grid = plt.GridSpec(4, 4, hspace=0, wspace=0)
        main_ax = fig.add_subplot(grid[1:, :-1])
        y_hist = fig.add_subplot(grid[1:, -1], yticks=yticks, yticklabels=[], sharey=main_ax)
        x_hist = fig.add_subplot(grid[0, :-1], xticks=xticks, xticklabels=[], sharex=main_ax)
        zoom_ax = fig.add_subplot(grid[0, -1], yticklabels=[], sharex=y_hist)
        zoom_ax.yaxis.tick_right()
        zoom_ax.xaxis.tick_top()
        
        main_ax.tick_params(direction='in')

        image = self.image
        m, s = normFit( self.os().image.flatten() )
        m0, s0 = self.median(), self.MAD()
        print( '\n(MAD %.2f %.2f %.2f %.2f)'%( m, s, m0, s0 ), )
        if s0 == 0:
            s0 = n*s
        s0_ = n*s0

        mid = lambda _: .5*(_[1:]+_[:-1])
        diff = lambda _: _[1:]-_[:-1]
        
        ebins = np.arange(m0-s0_,m0+s0_,1)
        if log:
            image2 = np.abs(image)
            image2[ image2 == 0] = 1e-2
            logimage = np.log10( image2 )
            main_ax.imshow( logimage, vmin=1e-1, vmax=logimage.max() )
            x_hist.hist2d( np.repeat(range(image.shape[1]), image.shape[0] ), logimage.flatten(), bins=[image.shape[1], ebins], norm = colors.LogNorm() )
            y_hist.hist2d( logimage.flatten(), np.repeat(range(image.shape[0]), image.shape[1] ) , bins=[ebins, image.shape[0]], norm = colors.LogNorm() )
            hist, bins, _ = zoom_ax.hist( logimage.flatten(), bins=ebins )
        else:
            main_ax.imshow( image, vmin=m0-s0, vmax=m0+s0 )
            #print( '(min %.2f %.2f)'%( image.min(), image.max() ), )
            x_hist.hist2d( np.repeat(range(image.shape[1]), image.shape[0] ), image.T.flatten(), bins=[image.shape[1], ebins ], norm = colors.LogNorm() )
            x_hist.plot( range(image.shape[1]), np.median( image, axis = 0 ) )

            y_hist.hist2d( image.flatten(), np.repeat(range(image.shape[0]), image.shape[1] ) , bins=[ebins, image.shape[0]], norm = colors.LogNorm() )
            y_hist.plot( np.median( image, axis = 1 ), range(image.shape[0]) )

            hist, bins, _ = zoom_ax.hist( self.image.flatten(), bins=ebins, histtype='step', color='b' )
            zoom_ax.hist( self.os().image.flatten(), bins=ebins, histtype='step', color='r' )
            zoom_ax.hist( self.vos().image.flatten(), bins=ebins, histtype='step', color='g' )
            
            #zoom_ax2 = zoom_ax.twinx()
            #x = mid(bins)
            #x2 = x[ hist > hist/10 ]
            #x = mid(mid(bins))
            #dhist = diff(hist)
            #zoom_ax2.plot( x, dhist, color='b' )
            #zoom_ax2.plot( mid(mid(mid(bins))), diff(diff(hist)), color='g' )
        zoom_ax.set_yscale('log')
        a, b = zoom_ax.get_ylim()
        zoom_ax.set_ylim((.1, b))
        zoom_ax.grid(True)
        #main_ax.invert_yaxis()
        main_ax.set_xlim( (0,image.shape[1]) )
        main_ax.set_ylim( (0,image.shape[0]) )
        fig.savefig( fname )
        print( 'savefig', fname, )
        fig.clf()
        plt.close()
        return self
    
    def printMinLine( self, fname ):
        fig = plt.figure()
        
        bins = np.arange( -100, 200, 3 )

        std_lines = np.std( self.L().osSub(axis=None).ac().image, axis = 1 )
        minstd_lines = np.argwhere( std_lines == np.min(std_lines) )[0]
        
        std_cols = np.std( self.L().osSub(axis=None).ac().image, axis = 0 )
        minstd_cols = np.argwhere( std_cols == np.min(std_cols) )[0]

        std_colsOS = np.std( self.L().osSub(axis=None).os().image, axis = 0 )
        minstd_colsOS = np.argwhere( std_colsOS == np.min(std_colsOS) )[0]
        
        before_ax = fig.add_subplot(221)
        before_ax.set_title('raw left')
        before_ax.hist( self.L().osSub(axis=None).ac().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='ac lin' )
        before_ax.hist( self.L().osSub(axis=None).ac().image[ :, minstd_cols ], bins=bins, histtype='step', normed=True, label='ac col' )
        before_ax.hist( self.L().osSub(axis=None).os().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='os lin' )
        before_ax.hist( self.L().osSub(axis=None).os().image[ :, minstd_colsOS ], bins=bins, histtype='step', normed=True, label='os col' )
        before_ax.set_yscale('log')
        before_ax.legend( fancybox=True, framealpha=0 )

        before_ax = fig.add_subplot(222)
        before_ax.set_title('raw right')
        before_ax.hist( self.R().osSub(axis=None).ac().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='ac lin' )
        before_ax.hist( self.R().osSub(axis=None).ac().image[ :, minstd_cols ], bins=bins, histtype='step', normed=True, label='ac col' )
        before_ax.hist( self.R().osSub(axis=None).os().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='os lin' )
        before_ax.hist( self.R().osSub(axis=None).os().image[ :, minstd_colsOS ], bins=bins, histtype='step', normed=True, label='os col' )
        before_ax.set_yscale('log')
        before_ax.legend( fancybox=True, framealpha=0 )


        after_ax = fig.add_subplot(223)
        after_ax.set_title('biasSub left')
        after_ax.hist( self.L().biasImageSub().ac().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='ac lin' )
        after_ax.hist( self.L().biasImageSub().ac().image[ :, minstd_cols ], bins=bins, histtype='step', normed=True, label='ac col' )
        after_ax.hist( self.L().biasImageSub().os().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='os lin' )
        after_ax.hist( self.L().biasImageSub().os().image[ :, minstd_colsOS ], bins=bins, histtype='step', normed=True, label='os col' )
        after_ax.set_yscale('log')
        after_ax.legend( fancybox=True, framealpha=0 )
        
        after_ax = fig.add_subplot(224)
        after_ax.set_title('biasSub right')
        after_ax.hist( self.R().biasImageSub().ac().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='ac lin' )
        after_ax.hist( self.R().biasImageSub().ac().image[ :, minstd_cols ], bins=bins, histtype='step', normed=True, label='ac col' )
        after_ax.hist( self.R().biasImageSub().os().image[minstd_lines,:].T, bins=bins, histtype='step', normed=True, label='os lin' )
        after_ax.hist( self.R().biasImageSub().os().image[ :, minstd_colsOS ], bins=bins, histtype='step', normed=True, label='os col' )
        after_ax.set_yscale('log')
        after_ax.legend( fancybox=True, framealpha=0 )

        fig.savefig( fname )
        print( '(fig %s)'%fname, )
        plt.close()
        return self

    
    def printSpectrum( self, fname, gain=2000, N=300 ):
        fig = plt.figure()
        
        bins_ = np.logspace( 1, np.log10(10.*gain), N )

        pixel_ax = fig.add_subplot(311)
        pixel_ax.set_yscale('log')
        pixel_ax.set_xscale('log')
        hist, bins, _ = pixel_ax.hist( self.image.flatten(), bins=bins_ )

        thr = 15
        bins_ = np.linspace( 4*thr, 10*gain, N/6 )
        hits = processImage3( self.image, threshold=thr )
        if hits is not None:
            print( 'E', hits['E0'].min(), hits['E0'].max(), len(hits['E0']), len(hits['E0'][ hits['n0'] > 5 ]), )
            
            spectrum_ax = fig.add_subplot(312)
            spectrum_ax.set_yscale('log')
            hist, bins, _ = spectrum_ax.hist( hits['E0'], bins=bins_, label='all' )
            hist, bins, _ = spectrum_ax.hist( hits['E0'][ hits['n0'] > 10 ], bins=bins_, label='>10' )
            hist, bins, _ = spectrum_ax.hist( hits['E0'][ hits['n0'] > 20 ], bins=bins_, label='>20' )
            hist, bins, _ = spectrum_ax.hist( hits['E0'][ hits['n0'] == 1 ], bins=bins_, label='==1' )
            spectrum_ax.legend( fancybox=True, framealpha=0 )

            spectrum_ax = fig.add_subplot(313)
            spectrum_ax.set_yscale('log')
            hist, bins, _ = spectrum_ax.hist( hits['E3'], bins=bins_, label='all' )
            hist, bins, _ = spectrum_ax.hist( hits['E3'][ hits['n0'] > 10 ], bins=bins_, label='>10' )
            hist, bins, _ = spectrum_ax.hist( hits['E3'][ hits['n0'] > 20 ], bins=bins_, label='>20' )
            hist, bins, _ = spectrum_ax.hist( hits['E3'][ hits['n0'] == 1 ], bins=bins_, label='==1' )
            spectrum_ax.legend( fancybox=True, framealpha=0 )
        
        fig.savefig( fname )
        print( 'savefig', fname )
        plt.close()
        return self

    @staticmethod
    def imageCov( image1, image2 ):
        return np.sum( image1.image.flatten()*image2.image.flatten() )
    
    @staticmethod
    def computeohduMedian( runIDs, ohdus ):
        for ohdu in ohdus:
            images = [ hitSumm.Image(runID=runID, ohdu=ohdu) for runID in runIDs ]
            imagesR = [ image.R().osSub(np.median,axis=1).image for image in images ]
            imagesL = [ image.L().osSub(np.median,axis=1).image for image in images ]
            
            medianR = np.median( imagesR, axis=0 )
            medianL = np.median( imagesL, axis=0 )
            
            writeImageFITS( medianR, header = {'ohdu':ohdu}, fname='medianR_runID%05d-%05d_ohdu%s.fits'%( min(runIDs),max(runIDs), ohdu ) )
            #print( readSingleImageFITS( fname='medianR_runID%05d-%05d_ohdu%s.fits'%( min(runIDs),max(runIDs), ohdu ) )[0].shape )
            writeImageFITS( medianL, header = {'ohdu':ohdu}, fname='medianL_runID%05d-%05d_ohdu%s.fits'%( min(runIDs),max(runIDs), ohdu ) )
            #print( readSingleImageFITS( fname='medianL_runID%05d-%05d_ohdu%s.fits'%( min(runIDs),max(runIDs), ohdu ) )[0].shape )
            
            hitSumm.Image(medianR).printImage('medianR_runID%05d-%05d_ohdu%s.pdf'%( min(runIDs),max(runIDs), ohdu ))
            hitSumm.Image(medianL).printImage('medianL_runID%05d-%05d_ohdu%s.pdf'%( min(runIDs),max(runIDs), ohdu ))
        
    @staticmethod
    def matrixCorrelation( runIDs, ohdus ):
        images = {}
        dark, bias = {}, {}
        thr = 15
        for ohdu in ohdus:
            images[ohdu] = {}
            for runID in runIDs:
                tmp = hitSumm.Image( runID = runID, ohdu = ohdu )
                images[ohdu][runID] = { 'R': tmp.R().osSub(axis=None), 'L': tmp.L().osSub(axis=None) }
            bias[ohdu] = np.median( [ v['R'].image for k,v in images[ohdu].items() ], axis = 0 )
            dark[ohdu] = np.median( [ v['L'].image for k,v in images[ohdu].items() ], axis = 0 )
            
            E = None
            E2 = None
            E3 = None
            
            print( '\nohdu', ohdu )
            hitSumm.Image( dark[ohdu] ).printImage( 'images/darkFrame_runID%05d-%05d_ohdu%s.pdf'%( min(runIDs), max(runIDs), ohdu ) )
            hitSumm.Image( bias[ohdu] ).printImage( 'images/biasFrame_runID%05d-%05d_ohdu%s.pdf'%( min(runIDs), max(runIDs), ohdu ) )
            hitSumm.Image( dark[ohdu] - bias[ohdu] ).printImage( 'images/darkCurrent_runID%05d-%05d_ohdu%s.pdf'%( min(runIDs), max(runIDs), ohdu ) )
            hitSumm.Image( hitSumm.submad( dark[ohdu], bias[ohdu] ) ).printImage( 'images/darkCurrentMAD_runID%05d-%05d_ohdu%s.pdf'%( min(runIDs), max(runIDs), ohdu ) )
            for runID in runIDs:
                CNoise = hitSumm.Image( np.median([ hitSumm.Image( runID = runID, ohdu = ohdu_ ).R().osSub(axis=None).image for ohdu_ in ohdus ], axis=0) ).printImage('images/CN_runID%s.pdf'%(runID) )
                ReadOutCNoise = hitSumm.Image( images[ohdu][runID]['R'].image - CNoise.image ).printImage('images/RCN_runID%s_ohdu%s.pdf'%(runID, ohdu))
                CNoiseSub = hitSumm.Image( images[ohdu][runID]['L'].image - CNoise.image ).printImage('images/CNoiseSub_runID%s_ohdu%s.pdf'%(runID, ohdu))
                ReadOutNoise = hitSumm.Image( images[ohdu][runID]['R'].image - bias[ohdu] ).printImage('images/RN_runID%s_ohdu%s.pdf'%(runID, ohdu))
                BiasSub = hitSumm.Image( images[ohdu][runID]['L'].image - bias[ohdu] ).printImage('images/BiasSub_runID%s_ohdu%s.pdf'%(runID, ohdu))
                DarkSub = hitSumm.Image( images[ohdu][runID]['L'].image - dark[ohdu] ).printImage('images/DarkSub_runID%s_ohdu%s.pdf'%(runID, ohdu))
                BiasSub2 = hitSumm.Image( hitSumm.submad( BiasSub.image, images[ohdu][runID]['R'].image) ).printImage('images/BiasSubCor_runID%s_ohdu%s.pdf'%(runID, ohdu))
                
                hits = BiasSub.extract(thr)
                hits2 = BiasSub2.extract(thr)
                hits3 = DarkSub.extract(thr)
                if E is None:
                    E = hits['E0']
                    E2 = hits2['E0']
                    E3 = hits3['E0']
                else:
                    E = np.concatenate( (E, hits['E0']) )
                    E2 = np.concatenate( (E2, hits2['E0']) )
                    E3 = np.concatenate( (E3, hits3['E0']) )
                hitSumm.printSpectrum( 'images/spectrumBiasSub_ohdu%s.pdf'%ohdu, E, [4*thr, 10*2000] )
                hitSumm.printSpectrum( 'images/spectrumBiasSubCor_ohdu%s.pdf'%ohdu, E2, [4*thr, 10*2000] )
                hitSumm.printSpectrum( 'images/spectrumDarkSub_ohdu%s.pdf'%ohdu, E3, [4*thr, 10*2000] )
        return
    
    @staticmethod
    def _printSpectrum( fname, data, xlim ):
        fig = plt.figure()
        main_ax = fig.add_subplot(211)
        
        n = 100
        
        bins = np.arange( xlim[0], xlim[1], n )
        main_ax.hist( data, bins )
        main_ax.set_yscale('log')
        main_ax.set_xlabel('E[adu]')
        main_ax.set_ylabel('dN/dE[1/adu]')

        main_ax = fig.add_subplot(212)
        
        bins = np.arange( xlim[0], 25*xlim[1], 25*n )
        main_ax.hist( data, bins )
        main_ax.set_yscale('log')
        main_ax.set_xlabel('E[adu]')
        main_ax.set_ylabel('dN/dE[1/adu]')
        
        fig.savefig(fname)
        plt.close()
        return
    
    def spectrumDiff( self, bins ):
        print( 'spectrumDiff' )
        #self.reconstructionEfficiency()
        
        #hits = self.getSCNHits_adu( self.config.runIDmin )
        #print( 'E', hits['E0'].min(), hits['E0'].max(), len(hits['E0']), )
        #fig = plt.figure()
        #ax = fig.add_subplot(211)
        #ax.set_yscale('log')
        #ax.hist( hits['E0'], bins=np.linspace(0, 1000000, 100) )
        #ax = fig.add_subplot(212)
        #ax.set_yscale('log')
        #ax.set_xscale('log')
        #ax.hist( hits['E0'], bins=np.logspace(1, 6, 100) )
        #fig.savefig( 'spm_' + 'scnHits.pdf' )
        
        #hitSumm.computeohduMedian( range(4000,4010), [2,3] )
        Image.matrixCorrelation( range(4000,4010), [2,3,4,5,6,7,8,9,10] )
        exit(0)
        #5825--5828
        histTotal = None
        hist2dTotal = None
        Etotal = None
        ntotal = None
        
        xbins = None
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_yscale('log')
        #for runID in range(6019,6030):
        for runID in range(4000,4101):
            image = Image( self.getRawImage_adu( runID=runID, ohdu=3 ) )
            image.osi_scn()
            hist, xbins, hist2d, dx, dy, E, n = image.L().osSub(func=np.median).printImage( 'image_runID%s.pdf'%runID )
            if histTotal is None:
                xbins = xbins
                histTotal = np.array( hist )
                histTotal2d = np.array( hist2d )
                Etotal = E
                ntotal = n
            else:
                histTotal += np.array( hist )
                histTotal2d += np.array( hist2d )
                Etotal = np.concatenate( ( Etotal, E ) )
                ntotal = np.concatenate( ( ntotal, n ) )
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            mask = Etotal < 2e4
            ax2.hist2d( Etotal[mask], ntotal[mask], bins=[50,50], norm=colors.LogNorm() )
            fig2.savefig( 'spm2d_total031.pdf' )
            fig2.clf()

        ax.step( xbins, histTotal, where='post', label='6019,6030' )
        fig.savefig( 'spm_total031.pdf' )
        exit(0)

        histTotal = None
        for runID in range(5825,5829):
            image = Image( self.getRawImage_adu( runID=runID, ohdu=3 ) )
            hist, xbins, hist2d, dx, dy = image.L().osSub(func=np.median).printImage( 'image_runID%s.pdf'%runID )
            if histTotal is None:
                xbins = xbins
                histTotal = np.array( hist )
            else:
                histTotal += np.array( hist )
        ax.step( xbins, histTotal, where='post', label='5825,5829' )
        ax.legend(fancybox=True, framealpha=0)
        fig.savefig( 'spm_total.pdf' )
        fig.close()

        
        exit(0)
        #image = self.Image(self.getSCNImage_adu( self.config.runIDmin ))
        #print( 'scn.shape', image.image.shape, )
        #image.printImage( 'scn.pdf' )
        for dc in [.1,.2,.4,.6, 2.]:
            self.config.darkCurrent = dc
            dcE = dc*self.exposure
            imageDC = self.getDarkCurrentImage_adu()
            #res0 = self.Image( image=imageDC ).printImage( 'image_%s_%s.pdf'%(dc, 0) )
            for rn in [ 1., 1.8, 2.5]:
                #self.adu_keV = 1000
                print( 'options', dcE, rn, rn*self.adu_keV*self.config['e-[keV]'], self.adu_keV*self.config['e-[keV]'] )
                self.config.readoutNoiseWidth = rn
                
                imageRN = self.getNoiseImage_adu()
                imageCu = self.getCopperImage_adu()
                #res0 = self.Image( image=imageRN ).L().electron()
                #res1 = self.Image( image=imageDC+imageRN ).L().electron()
                #res2 = self.Image( image=imageDC+imageRN ).L().electron()
                #print( 'dc', res0[0]/(3*eV_e/1e3*self.adu_keV), res1[0]/(3*eV_e/1e3*self.adu_keV) )
                self.Image( image=imageDC+imageRN+imageCu ).printImage( 'image_%s_%s.pdf'%(dc, rn) )
        exit(0)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('gain')
        ax.set_ylabel('sigma')
        for runID in range( self.config.runIDmin, self.config.runIDmax+1):
            for ohdu in [2,3,4,5,6,7,8,9]:
                print( runID, ohdu )
                #res = self.Image( runID=runID, ohdu=ohdu ).osi_scn()
                res = self.Image( image=self.getSCNImage_adu( runID, ohdu=ohdu ) ).L().electron()
                gain = self.readGain(runID, ohdu)
                if res > 0:
                    ax.plot( gain, res, 'x', color='C%d'%ohdu )
                    fig.savefig('gain.png')
        exit(0)
        
        x = lambda a: .5*(a[:-1]+a[1:])

        #for key, value in result.items():
            #print( key, scipy.stats.norm.fit(value) )
            
        def std( x, axis=0 ):
            if len(x) == 1:
                return 0
            return np.std(x,axis=axis)
        
        
        meanR_active_line, stdR_active_line = scipy.stats.norm.fit(result['R_active'][3500,1000:-1000].astype(float) )
        meanR_active, stdR_active = scipy.stats.norm.fit(result['R_active'])
        print( 'meanR_active', meanR_active, stdR_active )
        meanL_over, stdL_over = scipy.stats.norm.fit(result['L_over'])
        print( 'meanL_over', meanL_over, stdL_over )
        
        result['L_over_sub'] = result['L_over'] - result['R_over']
        mean_over_sub, std_over_sub = scipy.stats.norm.fit(result['L_over_sub'])
        print( 'L_over_sub', mean_over_sub, std_over_sub )
        print( 'L_over_sub_line', scipy.stats.norm.fit(result['L_over_sub'][1000,10:-10]) )
        
        result['R'] = result['R'].astype(float) - (meanR_active)
        result['R_active_line'] = result['R_active'][3500,1000:-1000].astype(float) - mean_over_sub
        result['R_active'] = result['R_active'].astype(float)
        result['L_active'] = result['L_active'].astype(float) - mean_over_sub
        #result['L_over'] = result['L_over'].astype(float) - (meanL_over)
        
        print( 'L_vOver_line', scipy.stats.norm.fit(result['L_vOver'][100,1000:-1000]) )
        print( 'R_vOver_line', scipy.stats.norm.fit(result['R_vOver'][100,1000:-1000]) )
        print( 'L_vOver0_line', scipy.stats.norm.fit(result['L_vOver0'][0,1000:-1000]) )
        print( 'R_vOver0_line', scipy.stats.norm.fit(result['R_vOver0'][0,1000:-1000]) )
        print( 'L_over_line', scipy.stats.norm.fit(result['L_over'][1000,10:-10]) )
        print( 'R_over_line', scipy.stats.norm.fit(result['R_over'][1000,10:-10]) )
        print( 'R_active_line', scipy.stats.norm.fit( result['R_active_line'] ) )
        print( 'R_active', scipy.stats.norm.fit(result['R_active']) )
        
        result['L_active_sub'] = result['L_active'] - result['R_active']

        fig = plt.figure()
        ax = fig.add_subplot(111)
        data = result['L_active_sub'].astype(float)
        data -= np.min(data)
        data[data<1] = 1
        #data = data[1*1055+1000:2*1055+1110, 0:2000]
        ax.imshow( np.log10(data) )
        fig.savefig('mergedImageCorrected.png')
        

        slices = [ 'R_active', 'L_active_sub', 'L_over' ]
        for key in slices:
            #result[key] = result[key].astype(float)
            #result[key] -= scipy.stats.norm.fit( result[key].flatten() )[0]
            #print( key+'.shape', result[key].shape )
            mean, std = scipy.stats.norm.fit( result[key] )
            print( key, mean, std )
            result[key+'0'] = np.mean( result[key], axis=0 )
            #print( key+'0.shape', result[key+'0'].shape )
            mean, std = scipy.stats.norm.fit( result[key+'0'] )
            print( key+'0', mean, std )
            result[key+'1'] = np.mean( result[key], axis=1 )
            #print( key+'1.shape', result[key+'1'].shape )
            mean, std = scipy.stats.norm.fit( result[key+'1'] )
            print( key+'1', mean, std )
            
        mean, std = scipy.stats.norm.fit( result['L_active_sub'][ result['L_active_sub'] <= 0 ], loc=0 )
        print( 'negs', mean, std )
        negs = result['L_active_sub'].astype(float)
        mean, std = scipy.stats.norm.fit( result['L_active_sub'][ (result['L_active_sub'] - 50) < 20 ] )
        print( 'negs2', mean, std )
        negs -= mean
        #negs = negs[ negs <= -int(std) ]
        negs_ = negs[ negs <= 0 ]
        negs = negs[ negs < 0 ]
        print( 'negs', scipy.stats.norm.fit( negs ) )
        result['L_active_sub_noise'] = np.concatenate( (negs_, -negs) )
        print( 'L_active_sub_noise', scipy.stats.norm.fit( result['L_active_sub_noise'] ) )
        
        slices.append('L_active_sub_noise')
        slices.append('R_active_line')
        
        dE = 1
        for key in slices:
            result[key] = result[key].flatten()

            E0 = dE
            E0 = result[key].min()
            E1 = result[key].max()
            print( 'E', E0, E1 )
            E = np.arange( E0, E1, dE )
            L = int((E1-E0)/dE)
            
            result[key+'.dN_dE'], result[key+'.E'] = np.histogram( result[key], E )
            result[key+'.dN_dE'][ result[key+'.dN_dE'] < 1 ] = 1
            result[key+'.dN_dE'] = np.log10( result[key+'.dN_dE'] )
            
            args = np.argwhere( np.abs(x(result[key+'.E'])) < 100 )
            result[key+'.dN_dE'] = result[key+'.dN_dE'][args]
            result[key+'.E'] = x(result[key+'.E'])[args]
            print( result[key+'.dN_dE'].shape, result[key+'.E'].shape )
        
        max_ = result['L_over_sub'].max()
        print( 'max', max_ )
        
        max_2 = result['L_active_sub'].max()
        
        histogram = lambda y: smooth2( np.histogram( y, bins = np.arange( 0, max_2, dE ) )[0], radius = 0 )
        fourier = lambda y: ( lambda x: [ np.fft.fftfreq( len(x) )/dE, x ] )( np.fft.fft( y ) )

        freqA, dN_dfA = fourier( histogram( result['L_active_sub'][ result['L_active_sub'] > 0 ].flatten() ) )
        indicesA = np.argsort(freqA[freqA>0])
        freqA = freqA[indicesA]
        dN_dfA = dN_dfA[indicesA]
        
        freq, dN_df = fourier( histogram( -result['L_active_sub'][ result['L_active_sub'] < 0 ].flatten() ) )
        indices = np.argsort(freq)
        freq = freq[indices]
        dN_df = dN_df[indices]

        freq0, dN_df0 = fourier( histogram( result['L_over_sub'].flatten() ) )
        indices = np.argsort(freq0[freq0>0])
        freq0 = freq0[indices]
        dN_df0 = dN_df0[indices]
        
        idN_dE = np.fft.ifft( np.abs(dN_df)/np.abs(dN_df0) )
        iE = np.arange( 0, max_, dE )
        print( idN_dE.shape, x(iE).shape )
        
        return [ [ 1./freq, np.log10(dN_df)], [ 1./freq0, np.log10(dN_df0)], [ 1./freq, np.log10(dN_df/dN_df0) ], [ 1./freqA, np.log10(dN_dfA) ] ]
    
        flatten = lambda l: [item for sublist in l for item in sublist]
        print( flatten([ [sl+'0', sl+'1'] for sl in slices ]) )
        return [ [ result[key[0]], result[key[1]] ] for key in [ [sl+'.E', sl+'.dN_dE'] for sl in slices ] ]

    @staticmethod
    def fourier( x, dE ):
        F = np.fft.fft( x )
        return np.fft.fftfreq( len(F), dE ), F
    
    @staticmethod
    def invfourier( F, dE ):
        return np.fft.fftfreq( len(F), 1./dE ), np.fft.ifft(F)
    
    @staticmethod
    def deconvolve( y, y2 ):
        f = np.fft.fftshift( np.fft.fft( y ) )
        f2 = np.fft.fftshift( np.fft.fft( y2 ) )
        print( 'shapes', f.shape, f2.shape )
        return np.fft.fftshift( np.abs(np.fft.ifft( f/f2 )) )

class OverScanImage(ImageBase):
    def __init__(self, image ):
        ImageBase.__init__(self, image)
        self.name = 'overscan'
        if not self.width in [ 150, 550 ]:
            print 'unpredicted overscan shape', self.shape
            exit(1)

class DataImage(ImageBase):
    def __init__( self, image ):
        ImageBase.__init__( self,image )
        self.name = 'data'
        
    def void( self, thr=15 ):
        himage = self.image[:,self.image.shape[1]/2:]
        clusters = scipy.ndimage.label( himage >= 4*thr, structure=[[0,1,0],[1,1,1],[0,1,0]] )[0]
        print 'clusters shape', clusters.shape
        distance_from_thr = scipy.ndimage.distance_transform_edt( clusters==0 )
        #print distance_from_thr.shape
        x, y = np.unravel_index( np.argmax( distance_from_thr ), himage.shape )
        print 'void max dist', distance_from_thr[x,y]
        image = himage[ distance_from_thr >= 2 ]
        if image.size == 0:
            return None
        return ArrayExtension( image )
    
class ActiveImage( ImageBase ):
    def __init__( self, image ):
        ImageBase.__init__(self,image)
        self.width = self.image.shape[1]
        self.name = 'active'
        
        if self.image.shape[0] in [4220, 1055, 3165]:
            self.overscanHeight = 90
        elif self.image.shape[0] == 900:
            self.overscanHeight = 70
        else:
            print 'unpredicted active shape', self.shape
            exit(1)

        self.dataSection = np.s_[ :-self.overscanHeight,:]
        self.halfdataSection = np.s_[ :-self.overscanHeight,self.width/2:]
        self.vos = np.s_[-self.overscanHeight:,:]
        self.hvos = np.s_[-self.overscanHeight:,self.width/2:]
        
    def data( self ):
        return DataImage( self.image[self.dataSection] )

    def halfdata( self ):
        a = DataImage( self.image[self.halfdataSection] )
        a.name = 'hdata'
        return a

    def verticalOverScan( self ):
        a = ImageBase( self.image[self.vos] )
        a.name = 'vertOS'
        return a

    def halfverticalOverScan( self ):
        a = ImageBase( self.image[self.hvos] )
        a.name = 'hvertOS'
        return a

    def median2( self, axis=1 ):
        return np.median( self.image[:,-self.width/2:], axis=axis )
    
class SideImage(ImageBase):
    def __init__( self, image, name=None ):
        ImageBase.__init__(self,image)
        self.width = self.image.shape[1]
        self.trim = 8
        self.name = name
        
        if self.image.shape[1] in [4262, 4350]:
            self.overscanWidth = 150
        elif self.image.shape[1] == 4662:
            self.overscanWidth = 550
        else:
            print( 'side shape not predited', self.image.shape )
            exit(0)
        if self.image.shape[0] in [4220, 1055, 3165, 1800]:
            self.overscanHeight = 70
        elif self.image.shape[0] == 900:
            self.overscanHeight = 70
        else:
            print 'unpredicted active shape', self.shape
            exit(1)

        self.os = np.s_[:,-self.overscanWidth:]
        self.ac = np.s_[:,:-self.overscanWidth]
    
        
    def active( self ):
        return ActiveImage( self.image[self.ac] )
    
    def overscan( self ):
        return OverScanImage( self.image[self.os] )

    def diffMAD( self ):
        image = self.active().image - self.active().v_medians()[None:]
        std = np.mean( stats.MAD( image, axis=1 ) )
        return std**2 - self.readoutNoiseMAD()**2

    def sigmaOS( self ):
        image = self.overscan().image
        sigma = scipy.stats.norm.fit( image.flatten() )[1]
        #print '(SideImage.sigmaOS)sigma=%s'%sigma
        return sigma

    def sigmaOSbin( self ):
        #print '(SideImage.sigmaOSbin)'
        image = self.overscan().image
        #print '(SideImage.sigmaOSbin)image.shape=', image.shape
        data = image.flatten()
        data = data[ data < 1e9 ]
        #print '(SideImage.sigmaOSbin)data.minmax=', data.min(), data.max()
        bins = np.arange(data.min(), data.max(), 1)
        hist, _ = np.histogram(data, bins)
        #print '(SideImage.sigmaOSbin)data.size=', data.size
        fit_func = lambda x,A,mu,sig: A*scipy.stats.norm.pdf(x, mu, sig)
        p0 = [np.sum(data),0., 1.]
        params = scipy.optimize.curve_fit( fit_func, .5*(bins[:-1] + bins[1:]), hist, p0=p0 )[0]
        print '(SideImage.sigmaOSbin)params=', params
        return params[-1]
        
    def sigmaOSMAD( self ):
        image = self.overscan().image
        std = np.mean( stats.MAD( image, axis=1 ) )
        #print '(SideImage.sigmaOSbin)std=%s'%(std)
        return std

    def sigmaOSMAD2( self ):
        image = self.overscan().image
        std = stats.MAD( image.flatten() )
        return std
    
    def lambdaBGbin( self, sigma, gain, thr=15, osi=False ):
        #print '(SideImage.lambdaBGbin)sigma=', sigma, 'gain=', gain, 'thr=', 15
        
        subimage = self.image[:,:-self.overscanWidth ]
        shape = subimage.shape
        subimage = subimage[:,shape[1]/2:]
        if osi: 
            osicorrection = np.median( self.image[:,self.overscanWidth:], axis=1 )
            #print '(lambdaBGbin)shapeOSI=', osicorrection.shape
            subimage -= osicorrection[:,None]
        #print '(lambdaBGbin)shape=', shape, subimage.shape, 

        s33 = [[1,1,1], [1,1,1], [1,1,1]]
        clusters = scipy.ndimage.label( subimage >= 4*thr, structure=s33 )[0]
        #print 'clusters shape', clusters.shape
        distance_from_thr = scipy.ndimage.distance_transform_edt( clusters==0 )
        x, y = np.unravel_index( np.argmax( distance_from_thr ), subimage.shape )
        #print 'void max dist', distance_from_thr[ x, y ]
        data = subimage[ distance_from_thr >= 3 ]
        if data.size == 0:
            return None
        #print '(SideImage.lambdaBGbin)data.size=', data.size, 'subimage.size=', subimage.size, '%.2f' % (float(data.size)/subimage.size)
        
        bins = np.arange(data.min(), data.max(), 1)
        x = .5*(bins[:-1] + bins[1:])
        hist, _ = np.histogram(data, bins)
        maskCount = hist > 10
        
        #fig = plt.figure()
        #ax = fig.add_subplot(211)
        imageOS = self.overscan().image
        if osi: 
            osicorrection = np.median( self.overscan().image, axis=1 )
            #print '(lambdaBGbin)shapeOSI_OS=', osicorrection.shape
            imageOS -= osicorrection[:,None]
            
        dataOS = imageOS.flatten()
        dataOS = dataOS[ dataOS < 1e9 ]
        binsOS = np.arange(dataOS.min(), dataOS.max(), 1)
        histOS, _ = np.histogram(dataOS, binsOS)
        xOS = .5*(binsOS[:-1] + binsOS[1:])
        maskCountOS = histOS > 10
        fit_funcOS = lambda x,A,mu,sig: A*scipy.stats.norm.pdf(x, mu, sig)

        p0OS = [np.sum(histOS[maskCountOS]),0., 1.]
        paramsOS = scipy.optimize.curve_fit( fit_funcOS, xOS[maskCountOS], histOS[maskCountOS], p0=p0OS )[0]
        #print '(SideImage.lambdaBGbin)paramsOS=', paramsOS
        #ax.step( xOS[maskCountOS], histOS[maskCountOS], where='mid' )
        #ax.step( xOS[maskCountOS], paramsOS[0]*scipy.stats.norm.pdf(xOS[maskCountOS],paramsOS[1],paramsOS[2]), where='mid' )

        #ax = fig.add_subplot(212, sharex=ax )
        #ax.step( x[maskCount], hist[maskCount], where='mid' )
        A = float( np.sum(hist[maskCount]) )/np.sum(histOS[maskCountOS]) * paramsOS[0]
        #ax.step( x[maskCount], A * scipy.stats.norm.pdf( x[maskCount], paramsOS[1], paramsOS[2] ), where='mid' )
        #fig.savefig('lambdBGbin.png')
        #print '(SideImage.lambdaBGbin)image saved'
        
        fit_func = lambda x,A,loc,mu: A*stats.poisson_norm.pdf(x, loc, paramsOS[2], gain, mu )
        p0 = [A, paramsOS[1], .001]
        print '(SideImage.lambdaBGbin)fitting'
        try:
            params = scipy.optimize.curve_fit( fit_func, x[maskCount], hist[maskCount], p0=p0, bounds=[(0,-np.inf,0), (np.inf, np.inf, np.inf)] )[0]
            print '(SideImage.lambdaBGbin)params=', params
        except:
            print 'error in lambdaBGbin'
            return None
        
        return params[-1]

    def gainElectronbin( self, thr=15 ):
        subimage = self.image[:,:-self.overscanWidth ]
        shape = subimage.shape
        subimage = subimage[:,shape[1]/2:]
        osicorrection = np.median( self.image[:,self.overscanWidth:], axis=1 )
        subimage -= osicorrection[:,None]

        s33 = [[1,1,1], [1,1,1], [1,1,1]]
        clusters = scipy.ndimage.label( subimage >= 4*thr, structure=s33 )[0]
        #print 'clusters shape', clusters.shape
        distance_from_thr = scipy.ndimage.distance_transform_edt( clusters==0 )
        x, y = np.unravel_index( np.argmax( distance_from_thr ), subimage.shape )
        #print 'void max dist', distance_from_thr[ x, y ]
        data = subimage[ distance_from_thr >= 3 ]
        if data.size == 0:
            return None
        #print '(SideImage.lambdaBGbin)data.size=', data.size, 'subimage.size=', subimage.size, '%.2f' % (float(data.size)/subimage.size)
        
        bins = np.arange(data.min(), data.max(), 1)
        x = .5*(bins[:-1] + bins[1:])
        hist, _ = np.histogram(data, bins)
        maskCount = hist > 10
        
        imageOS = self.overscan().image
        osicorrection = np.median( self.overscan().image, axis=1 )
        imageOS -= osicorrection[:,None]
            
        dataOS = imageOS.flatten()
        dataOS = dataOS[ dataOS < 1e9 ]
        binsOS = np.arange(dataOS.min(), dataOS.max(), 1)
        histOS, _ = np.histogram(dataOS, binsOS)
        xOS = .5*(binsOS[:-1] + binsOS[1:])
        maskCountOS = histOS > 10
        fit_funcOS = lambda x,A,mu,sig: A*scipy.stats.norm.pdf( x, mu, sig )

        p0OS = [ np.sum(histOS[maskCountOS]), 0., 1.]
        paramsOS = scipy.optimize.curve_fit( fit_funcOS, xOS[maskCountOS], histOS[maskCountOS], p0=p0OS )[0]
        print '(SideImage.gainElectronbin)paramsOS=', paramsOS
        A = float( np.sum(hist[maskCount]) )/np.sum(histOS[maskCountOS]) * paramsOS[0]

        fit_func = lambda x,A,mu,sig: A*scipy.stats.norm.pdf( x, mu, sig )
        p0 = [ A, 0, 1. ]
        params = scipy.optimize.curve_fit( fit_func, x[maskCount], hist[maskCount], p0=p0, bounds=[(0,-np.inf,0), (np.inf, np.inf, np.inf)] )[0]
        print '(SideImage.gainElectronbin)params=', params
        
        gain = (params[2]**2 - paramsOS[2]**2)/(params[1] - paramsOS[1])
        print '(SideImage.gainElectronbin)gain=', gain
        return 

    def diffMedianBG( self, thr=15, osi=True ):
        subimage = self.image[:,:-self.overscanWidth ]
        shape = subimage.shape
        subimage = subimage[:,shape[1]/2:]
        if osi: 
            osicorrection = np.median( self.overscan().image, axis=1 )
            subimage -= osicorrection[:,None]

        s33 = [[1,1,1], [1,1,1], [1,1,1]]
        clusters = scipy.ndimage.label( subimage >= 4*thr, structure=s33 )[0]
        #print 'clusters shape', clusters.shape
        distance_from_thr = scipy.ndimage.distance_transform_edt( clusters==0 )
        #x, y = np.unravel_index( np.argmax( distance_from_thr ), subimage.shape )
        #print 'void max dist', distance_from_thr[ x, y ]
        subimage[ distance_from_thr < 3 ] = np.nan
        
        medianAC = np.nanmedian( subimage, axis=1 )
        
        imageOS = self.overscan().image
        if osi: 
            osicorrection = np.median( self.overscan().image, axis=1 )
            imageOS -= osicorrection[:,None]
            #print '(MADBG)shapeOSI_OS=', osicorrection.shape
            
        dataOS = imageOS#.flatten()
        dataOS[ dataOS > 1e9 ] = np.nan
        medianOS = np.nanmedian( dataOS, axis=1 )
        
        diffMedian = np.mean( medianAC - medianOS )
        print '(SideImage.lambdaBGbin)diffMedian=', diffMedian, np.mean(medianAC), np.mean(medianOS)
        if diffMedian < 0: return 0
        return np.sqrt( diffMedian )

    def gainElectron( self, thr=15 ):
        subimage = self.image[:,:-self.overscanWidth ]
        shape = subimage.shape
        subimage = subimage[:,shape[1]/2:]
        osicorrection = np.median( self.overscan().image, axis=1 )
        subimage -= osicorrection[:,None]

        s33 = [[1,1,1], [1,1,1], [1,1,1]]
        clusters = scipy.ndimage.label( subimage >= 4*thr, structure=s33 )[0]
        distance_from_thr = scipy.ndimage.distance_transform_edt( clusters==0 )
        subimage[ distance_from_thr < 3 ] = np.nan
        
        medianAC = np.nanmedian( subimage, axis=1 )
        
        imageOS = self.overscan().image
        osicorrection = np.median( self.overscan().image, axis=1 )
        imageOS -= osicorrection[:,None]
            
        dataOS = imageOS
        dataOS[ dataOS > 1e9 ] = np.nan
        medianOS = np.nanmedian( dataOS, axis=1 )
        
        diffMedian = np.mean( medianAC - medianOS )
        diffVar = np.mean( medianAC**2 - medianOS**2 ) - np.mean( medianAC + medianOS )*diffMedian
        gain = diffVar/diffMedian
        #print '(SideImage.gainElectron)diffMedian=', diffMedian, diffVar
        #print '(SideImage.gainElectron)gain=', gain
        
        diffVarMAD = np.mean( stats.nanMAD( subimage, axis=1 )**2 - stats.nanMAD( imageOS, axis=1 )**2 )
        gain_prime = diffVarMAD/diffMedian
        #print '(SideImage.gainElectron)diffMedian=', diffMedian, diffVarMAD
        #print '(SideImage.gainElectron)gain\'=', gain_prime
        
        #if diffMedian < 0: return 0
        return gain

    def gainElectronMAD( self, thr=15 ):
        subimage = self.image[:,:-self.overscanWidth ]
        shape = subimage.shape
        subimage = subimage[:,shape[1]/2:]
        osicorrection = np.median( self.overscan().image, axis=1 )
        subimage -= osicorrection[:,None]

        s33 = [[1,1,1], [1,1,1], [1,1,1]]
        clusters = scipy.ndimage.label( subimage >= 4*thr, structure=s33 )[0]
        distance_from_thr = scipy.ndimage.distance_transform_edt( clusters==0 )
        subimage[ distance_from_thr < 3 ] = np.nan
        
        medianAC = np.nanmedian( subimage, axis=1 )
        
        imageOS = self.overscan().image
        osicorrection = np.median( self.overscan().image, axis=1 )
        imageOS -= osicorrection[:,None]
            
        dataOS = imageOS
        dataOS[ dataOS > 1e9 ] = np.nan
        medianOS = np.nanmedian( dataOS, axis=1 )
        
        diffMedian = np.mean( medianAC - medianOS )
        diffVarMAD = np.mean( stats.nanMAD( subimage, axis=1 )**2 - stats.nanMAD( imageOS, axis=1 )**2 )
        gain_prime = diffVarMAD/diffMedian

        return gain_prime

    def overscanLineSubtraction(self):
        return SideImage( self.image - self.overscan().h_medians() )

    def readoutNoiseMAD( self ):
        image = self.overscan().image - self.overscan().v_medians()[None:]
        std = np.mean( stats.MAD( image, axis=1 ) )
        return std
    
    def darkCurrentMedian( self ):
        var = np.mean( self.active().h2_medians() - self.overscan().h_medians() )
        if var < 0: return 0
        return np.sqrt( var )

    def biasImageCorrected( self ):
        #print 'side bias', self.active().biasImage2Corrected().image.shape
        return SideImage( np.concatenate( (self.active().biasImage2Corrected().image, self.overscan().biasImageCorrected().image), axis=1 ) )

    def verticalSubtraction(self):
        return SideImage( self.image - np.median( self.image[-self.overscanHeight+1:,:], axis=0 )[None,:] )

    def horizontalSubtraction(self):
        #return SideImage( self.image - np.median( self.image[:,-self.overscanWidth:], axis=1 )[:,None] )
        return SideImage( self.image - np.median( self.image[:,-self.overscanWidth:], axis=1 )[:,None] )

    def horizontalSubtraction(self):
        return SideImage( self.image - np.median( self.image[-self.overscanHeight:,:], axis=0 )[None,:] )

class FullImage(ImageBase):
    def __init__( self, runID, ohdu, imageType, debug=False ):
        self.has_right = True
        self.name = 'full'
        if debug: print '(FullImage.__init__)imageType=%s,runID=%s,ohdu=%s'%(imageType,runID,ohdu)
        if imageType.startswith('raw'):
            ImageBase.__init__(self, readImage_raw( runID=runID, ohdu=ohdu )[::-1,:] )
            if self.width in [ 8540, 9340, 8716, ]:
                self.trim = 8
                self.side =  np.s_[:, self.trim:self.width/2]
                if self.width in [8540, 8716]:
                    self.overscanWidth = 150
                else:
                    self.overscanWidth = 550
            else:
                print 'full shape not predited', self.image.shape, self.imageType
                exit(1)
                
        #elif imageType == 'osi':
            #ImageBase.__init__(self, hitSumm._getOSIImage_adu( runID=runID, ohdu=ohdu )[::-1,:] )
            #if self.width in [ 8524, 9324]:
                #self.trim = 0
                #self.side =  np.s_[:, :self.width/2]
            #else:
                #print( 'full shape not predited', self.image.shape )
                #exit(1)
                
        #elif imageType == 'mbs':
            #ImageBase.__init__(self, hitSumm._getMBSImage_adu( runID=runID, ohdu=ohdu )[::-1,:] )
            #if self.width in [ 8524, 9324 ]:
                #self.trim = 0
                #self.side =  np.s_[:, :self.width/2]
            #else:
                #print( 'full shape not predited', self.image.shape )
                #exit(1)

        elif imageType == 'MB':
            ImageBase.__init__(self, readImage_MB( runID=runID, ohdu=ohdu )[::-1,:] )
            if self.width in [ 8524, 9324 ]:
                self.trim = 0
                self.side =  np.s_[:, :self.width/2]
            else:
                print( 'full shape not predited', self.image.shape )
                exit(1)
        
        elif imageType.startswith('scn'):
            self.has_right = False
            #debug=True
            ImageBase.__init__(self, readImage_processed( runID=runID, ohdu=ohdu, image='scn', debug=debug )[::-1,:] )
            if self.width in [ 4262, 4662 ]:
                self.trim = 0
                self.side =  np.s_[:, :self.width]
            else:
                print( 'full shape not predited', self.image.shape )
                exit(1)
        else:
            print( 'imageType not predited', self.image.shape )
            exit(1)
        
        if imageType in [ 'raw*', 'scn*' ]: 
            self.correct_biasImage()

    def right( self ):
        if not self.has_right: return None
        return SideImage( self.image[:,::-1][self.side], 'right' )
    
    def left( self ):
        return SideImage( self.image[self.side], 'left' )

    def active( self ):
        return self.left().active()
    
    def overscan( self ):
        return self.left().overscan()
    
    def verticalOverScan(self):
        return self.left().active().verticalOverScan()

    def halfverticalOverScan(self):
        return self.left().active().halfverticalOverScan()

    def data( self ):
        return self.left().active().data()
    
    def void( self, thr=15 ):
        return self.left().active().data().void(thr)

    def part(self, i):
        return Image( self.image[i*1055:(i+1)*1055,:] )
    
    def correct_biasImage( self ):
        if not self.has_right: return self.left().biasImageCorrected()
        self.image = np.concatenate( ( self.left().biasImageCorrected().image, self.right().biasImageCorrected().image[:,::-1] ), axis=1 )
        return self
    
    def darkCurrentEstimate( self ):
        return self.left().darkCurrentMedian()

    def readoutNoiseEstimate( self ):
        return self.left().readoutNoiseMAD()
    
    def diffMAD(self):
        return self.left().diffMAD()
    
    def sigmaOS(self):
        return self.left().sigmaOS()

    def sigmaOSbin(self):
        return self.left().sigmaOSbin()

    def lambdaBGbin( self, sigma, gain, osi=False ):
        return self.left().lambdaBGbin( sigma, gain, osi=osi )

    def sigmaOSMAD(self):
        return self.left().sigmaOSMAD()

    def diffMedianBG(self):
        return self.left().diffMedianBG()
    
    def sigmaOSMAD2(self):
        return self.left().sigmaOSMAD2()
    
    def darkCurrentEstimate2( self, medianVoid = None, sigmaVoid=None, sigmaBase = None ):
        #if sigmaVoid is None: sigmaVoid = self.data().mad_half()
        #if medianVoid is None: medianVoid = self.data().median_half()
        #if sigmaBase is None: sigmaBase = self.overscan().mad()
        med = self.data().median_half() - self.overscan().median()
        Var = self.data().mad_half()**2 - self.overscan().MAD()**2
        lamb = med**2/Var #no unit
        gain = Var**2/med**3 #ADU
        sigma = self.overscan().MAD()/gain #no unit
        return gain, lamb, sigma

class FullImage2:
    def __init__(self, ac, os, hos ):
        self.active = ac
        self.overscan = os
        self.hoverscan = hos
    
    def sigma( self ):
        return self.overscan.std()
    
    def gain_lambda( self ):
        return self.active.loc()
    
    def mean_active_minus_overscan( self ):
        return self.active.mean() - self.overscan.mean()
    
    def mean_active_minus_overscan_line_by_line( self ):
        return np.mean( self.active.mean_line_by_line( (-self.overscan.number_of_lines(), None) ) - self.overscan.mean_line_by_line() )
    
def get_overscan_sizes( shape ):
    pass

def read_full_image( image, **kwargs ):
    overscan_lines, vertical_overscan_columns = get_overscan_sizes( image.shapes )
    active_section = image[ : -vertical_overscan_columns, : -overscan_lines ]
    overscan_section = image[ : -vertical_overscan_columns, overscan_lines : ]
    vertical_overscan_section = image[ vertical_overscan_columns : , : ]
    return FullImage2( active_section, overscan_section, vertical_overscan_section )
    
def simulate_full_image( **kwargs ):
    image = FullImage2()
    return image

def read_gainCu( runID, ohdu ):
    path = ConniePaths.get_gainCatalog_path( runID )
    if path is None: return None
    #print '(read_gainCu)runID=%s,ohdu=%s' % (runID, ohdu)
    return read_gainCu_from_gainCatalog( path, runID, ohdu )

def read_gainCu_from_gainCatalog( path, runID, ohdu ):
    #data = root_numpy.root2array( path, treename = 'hitSumm', branches = ['gainCu'], selection = 'runID==%d && ohdu==%d' % (runID, ohdu) )
    data = root_numpy.root2array( path, treename = 'hitSumm', branches = ['gainCu', 'ohdu', 'runID'] )
    #print '(read_gainCu_from_gainCatalog)len(data)', len(data)
    data = data['gainCu'][ np.all( [ data['runID'] == runID, data['ohdu'] == ohdu ], axis=0 ) ]
    #print '(read_gainCu_from_gainCatalog)len(data)', list(data)
    gainCu = list(set(data))
    #print '(read_gainCu_from_gainCatalog)gainCu=', gainCu
    #exit(1)
    if len(gainCu) == 0:
        return None
    return gainCu[0]
