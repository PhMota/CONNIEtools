import astropy
import astropy.io
import astropy.io.fits

import numpy as np
import scipy
import scipy.stats
import scipy.signal
import scipy.special
from scipy.misc import factorial
import scipy.optimize

from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import glob
import os
import sys
import re
import copy
import root_numpy
import datetime
import time

REMOVE = 0
REMOVEOLD = -1e9

class term:
    #HEADER = '\033[95m'
    #OKBLUE = '\033[94m'
    #OKGREEN = '\033[92m'
    #ARNING = '\033[93m'
    #FAIL = '\033[91m'
    cmd = lambda x: '\033[%dm'%x
    Q = cmd(0)
    B = cmd(1)
    I = cmd(3)
    U = cmd(4)
    #r = '\e[31m'
    r = cmd(91)
    g = cmd(32)
    b = cmd(34)
    
    red = staticmethod( lambda x: term.r + x + term.Q )
    bold = staticmethod( lambda x: term.B + x + term.Q )
    italic = staticmethod( lambda x: term.I + x + term.Q )
    
#folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/scn/merged/*'
folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029*/data_*/scn/merged/*'
folder2 = '/share/storage2/connie/data_analysis/processed02_data/runs/009*/data_*/scn/merged/*'

def listFITS( *patterns ):
    '''
    list all files that match the list of patterns given
    '''
    return sortByrunID([ match for pattern in patterns for match in glob.glob(pattern) if not '-with-' in match ])

def listFITS_run( *runs ):
    return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/scn/merged/scn_*.fits'%run for run in runs ] )

def listFITS_runID( *runIDs ):
    return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/scn_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )

def listremovedFITS_runID( *runIDs ):
    if os.access( '/share/storage2/connie/data_analysis/processed02_data/runs/', os.W_OK):
        print 'input path', '/share/storage2/connie/data_analysis/processed02_data/runs'
        return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/*hsub_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    print 'input path', '.'
    return listFITS( *[ 'run*/*hsub_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    #return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/*runID_*_%05d_*'%runID for runID in runIDs ] )

def getrunIDFromPath( path ):
    return re.search(r'runID_[0-9]+_([0-9]+)', os.path.basename( path )).groups()[0]

def getrunFromPath( path ):
    return re.search(r'runs/(.+)/data', path ).groups()[0]

def getrun_runIDFromPath( path ):
    return re.search(r'runID_([0-9]+)_([0-9]+)', os.path.basename( path )).groups()

def sortByrunID( paths ):
    '''
    sort list of files by runID
    '''
    return sorted( paths, key=getrunIDFromPath )

def splitCCD( ccd, noeventsCut = 8, overscanCut = 4112, vtrim=40, voverscanCut = 80, hspacing=12, overscanhspacing = 10 ):
    '''
    takes the full data from the CCD and separates the noevents, active region and overscan applying the proper trimmings
    '''
    noevents = ccd[ voverscanCut+vtrim:-vtrim, 0:noeventsCut ]
    data = ccd[ voverscanCut+vtrim:-vtrim, noeventsCut+hspacing:overscanCut-hspacing ]
    overscan = ccd[ voverscanCut+vtrim:-vtrim, overscanCut+overscanhspacing:-overscanhspacing ]
    #print( ccd.shape, noevents.shape, data.shape, overscan.shape )
    return noevents, data, overscan

def extractFromFITS( path ):
    '''
    extracts the data separated in the different regions from the FITS file
    '''
    ccds = astropy.io.fits.open( path )
    ohdus = { 'OHDU': {ccd.header['OHDU']: splitCCD(ccd.data)  for ccd in ccds}, 'TEMPMAX': ccds[0].header['TEMPMAX'], 'TEMPMIN': ccds[0].header['TEMPMIN'] }
    return ohdus

def getAssociatedCatalog( *paths ):
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/ext/catalog/catalog_data_*.root'%getrunFromPath(path)) if not 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedOSI( *paths ):
    return [ osi for path in paths for osi in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/osi/images/osi_runID_*_%s_*.fits'%(getrunFromPath(path), getrunIDFromPath(path)) ) if not '-with-' in osi ]

def runETA( msg, cmd, eta = None, N = None ):
    sys.stdout.write( msg + '...' )
    sys.stdout.flush()
    if not eta is None and not N is None:
        etastart = time.time()
        eta()
        et = time.time() - etastart
        print 'estimated time', time.strftime("%H:%M:%S", time.gmtime(int(et*N))), '(finishes at %s)'% time.strftime("%H:%M:%S", time.localtime( time.time() + int(et*N)))
    start = time.time()
    res = cmd()
    #print 'elapsed time', time.strftime("%H:%M:%S", time.gmtime(datetime.datetime.time() - start))
    print 'elapsed time', time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
    return res

def readCatalog( ROOTfile, runID = None, verbose = False ):
    min_max = lambda x: (np.min(x), np.max(x))
    if verbose: print 'reading start and stop indexes for runID', runID, 'in ROOT file'
    indexlist = np.argwhere( root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'] == runID ).flatten()
    start, stop = min_max( indexlist )
    if verbose: print 'reading', len(indexlist), 'entries from', start, 'to', stop
    
    readcatalog_once = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu'], selection='runID==%s'%runID, start=start_, stop=stop_+1 )
    if verbose: print 'reading ROOT file'
    hitscatalog = runETA( 
        msg = 'reading ROOTfile',
        cmd = lambda: readcatalog_once(start, stop), #read all entries and return 
        )
    return hitscatalog

def mergeParts( parts ):
    merged = parts[-1][:]
    #print [ (i,ohdu.header['OHDU']) for i, ohdu in enumerate(merged) ]
    for ohduindex in range(len(merged))[::-1]:
        if not merged[ohduindex].data is None:
            merged[ohduindex].data = np.concatenate( [ part[ohduindex].data for part in parts[::-1] ], axis = 0 )
        else:
            merged.pop(ohduindex)
    #print [ (i,ohdu.header['OHDU']) for i, ohdu in enumerate(merged) ]
    return merged

def removeHitsFromrunID( runID, outputfolder = None, ohdu = None, ROOTfile = None, osi = None, verbose = False, image = True ):
    FITSfile = listFITS_runID( runID )[0]
    print 'input FITSfile =', FITSfile
    outputfilenohits = 'hsub_'+os.path.basename(FITSfile)
    outputfilehits = 'hits_'+os.path.basename(FITSfile)
    
    if not outputfolder is None:
        if not os.path.exists(outputfolder):
            print 'creating directory', outputfolder
            os.makedirs( outputfolder )
    else:
        outputfolder = os.path.dirname(FITSfile)
        if os.access( outputfolder, os.W_OK):
            print 'outputfolder =', outputfolder
        else:
            outputfolder = 'run'+getrunFromPath( FITSfile ) + '/'
            print 'outputfolder =', outputfolder
            if not os.path.exists(outputfolder):
                print 'creating directory', outputfolder
                os.makedirs( outputfolder )
    
    removeOHDU = lambda x, ohdu_: map( lambda i: x.pop(i), [ j for j, hdu in enumerate(x) if not int(hdu.header['OHDU']) == int(ohdu) ][::-1] )
    if osi:
        OSIfiles = getAssociatedOSI( FITSfile )
        print 'OSI file ='
        print '\n'.join( OSIfiles )
    
        if verbose: print 'reading and merging OSI part files'
        merged = mergeParts( [ astropy.io.fits.open( OSIfile ) for OSIfile in OSIfiles ] )
        if verbose: print 'shape of merged', merged[-1].data.shape
        if not ohdu is None:
            removeOHDU( merged, ohdu )
    
    if verbose: print 'reading SCN file'
    hdulist = {}

    hdulist['hits'] = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        removeOHDU( hdulist['hits'], ohdu )
        #print 'remaining ohdus', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulist['hits'] ] ) )

    hdulist['nohits'] = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        removeOHDU( hdulist['nohits'], ohdu )
        #print 'remaining ohdus', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulist['nohits'] ] ) )
    
    #if not ohdu is None:
        #removeindex = [ i for i, hdu in enumerate(hdulist['hits']) if not int(hdu.header['OHDU']) == int(ohdu) ]
        #print 'selected ohdu', ohdu, 'removing ohdus', ', '.join(map( str, [ hdulist['nohits'][i].header['OHDU'] for i in removeindex ]))
        #for iohdu in removeindex[::-1]: #this inversion is needed so the removals work
            #hdulist['nohits'].pop(iohdu)
            #hdulist['hits'].pop(iohdu)
            #if osi: merged.pop(iohdu)
    #print 'remaining ohdus', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulist['nohits'] ] ) )
    #print 'remaining ohdus', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulist['hits'] ] ) )
    
    if ROOTfile is None:
        ROOTfile = getAssociatedCatalog( FITSfile )[0]
    print 'input ROOTfile =', ROOTfile
    hitscatalog = readCatalog( ROOTfile, runID=runID)

    if verbose: print 'removing hits'
    for iohdu in range(len(hdulist['hits'])):
        ohdu_ = hdulist['nohits'][iohdu].header['OHDU']
        print 'removing hits from', ohdu_
        hitsohdu = hitscatalog[ hitscatalog['ohdu']==ohdu_ ]
        hdulist['hits'][iohdu].data[:,:] = REMOVE #-1e9
        hits_xy = np.array( [ [iy,ix] for x,y in zip( hitsohdu['xPix'], hitsohdu['yPix']) for ix,iy in zip(x,y) ] )
        hits_e = np.array( [ ie for e in hitsohdu['ePix'] for ie in e ] )
        if len(hits_e) > 0:
            hdulist['nohits'][iohdu].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE #-1e9
            if osi: merged[iohdu].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE
            hdulist['hits'][iohdu].data[hits_xy[:,0],hits_xy[:,1]] = hits_e
        else:
            print term.red('Warning:'), 'empty hits_e on', ohdu_

    if image:
        figs = [ plt.figure() for i in range(2) ]
        map( lambda fig: fig.set_size_inches( np.array(fig.get_size_inches())*[2,1] ), figs)
        map( lambda fig: fig.suptitle('runID %s,ohdu %s'%(runID, hdulist['nohits'][0].header['OHDU'])), figs )
        ax = figs[0].add_subplot(121)
        ax2 = figs[1].add_subplot(111)
        ax.set_title('original')
        _imshow = lambda self, X: self.imshow(np.log(X-np.min(X)+1)[500:700,500:700], origin='lower')
        _hist = lambda self, X, bins, label: self.hist((X[X!=REMOVE]).flatten(), bins, histtype='step', label=label)
        
        X = astropy.io.fits.open(FITSfile)[0].data
        Xmin = X.min()
        h,bins,tmp = _hist(ax2,X[X<1e10],50000,'original')
        X[X==1e10] = 0
        _imshow(ax,X)
        
        ax = figs[0].add_subplot(122)
        ax.set_title('hits subtracted')
        X = hdulist['nohits'][0].data
        h,tmp,tmp2 = _hist(ax2,X[X<1e10],bins, 'subtracted')
        X[X==1e10] = 0
        _imshow(ax,X-Xmin)
        ax2.set_yscale('log')
        #ax2.set_xscale('log')
        #ax2.set_xlim([ bins[h>h.max()*1e-4].min(), bins[h>h.max()*1e-4].max() ])
        ax2.set_xlim([-100,100])
        ax2.legend()
        path = outputfolder + '/ohdu_%s_'%hdulist['nohits'][0].header['OHDU'] + os.path.basename(FITSfile) +'_ccd.png'
        path2 = outputfolder + '/ohdu_%s_'%hdulist['nohits'][0].header['OHDU'] + os.path.basename(FITSfile) +'_hist.png'
        print 'ccd image', path
        figs[0].savefig( path, bbox_inches='tight' )
        figs[1].savefig( path2, bbox_inches='tight' )
    
    hdulisthitsfile = outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilehits
    if os.path.exists(hdulisthitsfile):
        os.remove(hdulisthitsfile)
    if verbose: print 'writing SCN hits'
    hdulist['hits'].writeto( hdulisthitsfile )
    
    hdulistnohitsfile =  outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits 
    if os.path.exists(hdulistnohitsfile):
        os.remove(hdulistnohitsfile)
    if verbose: print 'writing SCN no hits'    
    hdulist['nohits'].writeto( hdulistnohitsfile )
    
    if osi:
        partsfile = outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits.replace('scn_mbs_', '')
        if os.path.exists(partsfile):
            os.remove(partsfile)
        if verbose: print 'writing merged OSI no hits'    
        merged.writeto( partsfile )
    
    print 'created new FITS files'
    print hdulisthitsfile
    print hdulistnohitsfile
    if osi: print partsfile

def fit( x, data, func, sel = None, log=False, p0 = None, labels = None ):
    mask = data > 0
    if sel: mask = np.all( sel+[ data>0 ], axis=0 )
    pp = scipy.optimize.curve_fit( func, x[mask], data[mask] if not log else np.log(data[mask]), p0 )[0]
    hist = func( x, *pp ) if not log else np.exp(func( x, *pp ))
    chisq = scipy.stats.chisquare( hist, f_exp = data )[0]/len(data)
    if labels:
        return {'chisqr': chisq, 'params': OrderedDict([ [label, pp_] for label, pp_ in zip(labels, pp) ] ), 'hist': hist }
    else:
        return pp, chisq, hist

__tab = ' '*4

def printDict( d, indent = 0, skip = None ):
    for key, value in d.items():
        if key in skip: continue
        if type(value) is dict or type(value) is OrderedDict:
            print __tab*indent+key
            printDict( value, indent+1, skip=skip)
        else:
            print __tab*indent+'%s: %s'%(key,value)

def printDictTable( d, indent = 0, skip = None, header = None, line = None ):
    tab = ' '*4
    doheader = True
    for key, value in d.items():
        if key in skip: continue
        if type(value) is dict or type(value) is OrderedDict:
            printDictTable( value, indent+1, skip=skip, header = filter( lambda key: not key in skip, value.keys() ), line = '%s:%s'%(line, key) if line else key )
        else:
            print ' '*len(line), ' '.join(map( lambda s: '%8s'%s, header))
            print line, ' '.join([ '%.2e'%(d[col]) for col in header])
            return
        doheader = False

def cPoisson( x, mu, lamb ):
    #if type(x) is float:
        #return np.abs( np.exp(x*np.log(lamb) -lamb -scipy.special.loggamma(x+1)) ) if x>=0 else 0
    res = np.zeros_like(x)
    x_ = x - mu
    res[x_>=0] = np.abs( np.exp(x_[x_>=0]*np.log(np.abs(lamb)) -np.abs(lamb) -scipy.special.loggamma(x_[x_>=0]+1)) )
    return res

def poisson( x, mu, lamb ):
    x -= mu
    return np.exp(x*np.log(lamb) -lamb -scipy.special.loggamma(x+1))
    
def plotrunID( ohdu, *FITSfiles, **kwargs ):
    fft = False
    if 'fft' in kwargs:
        fft = kwargs['fft']
    plot = True
    if 'plot' in kwargs:
        plot = kwargs['plot']
    runID = None
    runIDflag = False
    if 'runID' in kwargs:
        runID = kwargs['runID']
        runIDflag = True
        
    if runID is None:
        FITSfiles = FITSfiles[0]
    else:
        FITSfiles = listremovedFITS_runID(int(runID))
    
    print 'FITSfiles:\n', '\n'.join(FITSfiles)
    runID = getrunIDFromPath( FITSfiles[0] )
    print 'runID =', runID

    fitfunc = {
        'Gaussian(mu,sigma)*A': lambda x,mu,sigma,A: scipy.stats.norm.pdf(x,mu,sigma)*A,
        'log(Gaussian(mu,sigma)*A)': lambda x,mu,sigma,A: np.log(A/sigma/np.sqrt(2*np.pi)) - .5*(x-mu)**2/sigma**2,
        'Poisson(mu,lambda)*A': lambda x,mu,lamb,A: A*cPoisson(x, mu, lamb),
        }
    
    fitfunc['Gaussian(0,sigma)*A'] = lambda x,b,c: fitfunc['Gaussian(mu,sigma)*A'](x,0,b,c)
    fitfunc['log(Gaussian(0,sigma)*A)'] = lambda x,b,c: fitfunc['log(Gaussian(mu,sigma)*A)'](x,0,b,c)
    fitfunc['Poisson(0,lambda)*A'] = lambda x,lamb,A: fitfunc['Poisson(mu,lambda)*A'](x,0,lamb,A)
    fitfunc['Poisson(-lambda,lambda)*A'] = lambda x,lamb,A: fitfunc['Poisson(mu,lambda)*A'](x,-lamb,lamb,A)

    fitfunc['convolution(Gaussian(0,sigma),Poisson(-lamb,lamb))*A'] = lambda x, sigma, A, lamb: A*np.convolve( fitfunc['Gaussian(0,sigma)*A'](x,sigma,1), fitfunc['Poisson(-lambda,lambda)*A'](x,lamb,1)[::-1], 'same' )
    fitfunc['convolution(Gaussian(mu,sigma),Poisson(-lamb,lamb))*A'] = lambda x, mu, sigma, A, lamb: A*np.convolve( fitfunc['Gaussian(mu,sigma)*A'](x,mu,sigma,1), fitfunc['Poisson(-lambda,lambda)*A'](x,lamb,1)[::-1], 'same' )
    fitfunc['conv(G(0,sigma),P(-lamb,lamb))*A'] = fitfunc['convolution(Gaussian(0,sigma),Poisson(-lamb,lamb))*A']

    result = []
    header = []
    fmt = []
    first = True
    clean = lambda x: x.flatten()[ np.all( [x.flatten()!=REMOVEOLD, x.flatten()!=REMOVE, x.flatten()!=1e10], axis=0) ]

    osi = False
    meanplot = False
    print 'plotting'
    if 1:
        if first: 
            header += ['runID']
            fmt += ['%d']
        line = [ int(runID) ]
        if first: 
            header += ['ohdu']
            fmt += ['%d']
        line += [ int(ohdu) ]
        
        hdulists = [ astropy.io.fits.open( FITSfile ) for FITSfile in FITSfiles ]
        
        if len(hdulists) > 1:
            hdulists = OrderedDict([('scn', hdulists[0]), ('osi', hdulists[1]) ])
            osi = True
        else:
            hdulists = OrderedDict([('scn', hdulists[0])])
        print 'len hdulists', len(hdulists) 
        if first: 
            header += ['tempmin']
            fmt += ['%.2f']
        line += [ float(hdulists['scn'][0].header['TEMPMIN' ]) ]
        if first: 
            header += ['tempmax']
            fmt += ['%.2f']
        line += [ float(hdulists['scn'][0].header['TEMPMAX' ]) ]
        
        print 'ohdus in file', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulists['scn'] ] ))
        
        data = OrderedDict([ (key, hdu.data) for key, hdulist in hdulists.items() for hdu in hdulist if hdu.header['OHDU'] == ohdu ])
        if meanplot:
            print 'generating mean plot'
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.plot( np.min(data['scn'], axis=0) )
            
            if osi: ax2.plot( data['osi'][100,:] ) 
            ax2.set_xlim([1000,5500])
            fig2.savefig('mean.png')
        yccd = 4262
        xccd = 4220
        yos = 4112
        if osi:
            if data['osi'].shape[1] == 2*yccd:
                print 'subtracting idle amplifier'
                data['osisub'] = data['osi'][0:xccd,0:yccd] - np.flip(data['osi'][0:xccd,yccd:2*yccd], 1)
        print 'data len', data.keys()
        dx = 1
        shift = lambda bins: (bins[1:]+bins[:-1])*.5
        print 'bins', min(clean(data['scn'])), max(clean(data['scn']))
        bins = np.r_[min(clean(data['scn'])):max(clean(data['scn'])):dx]

        print 'computing histogram of full ccd'
        vslice = slice(120,4180)
        #print 'vslice', vslice
        hslices = OrderedDict( [('os', slice(4112+10, 4262-10)), ('ac', slice(20,4100))] )
        countCut = 10
        parseHist = lambda h,x: { 'y': h[h>countCut], 'x': (.5*(x[1:]+x[:-1]))[h>countCut] }
        #parseHist = lambda h,x: { 'y': h, 'x': (.5*(x[1:]+x[:-1])) }
        
        hist = OrderedDict( [ (datakey, { 
            #'all': {'hist': np.histogram( clean(datum[0:xccd,0:yccd]), bins = bins )[0] }, 'slice': '0:%s|0:%s'%(xccd,yccd),
            slicekey: {
                'hist': parseHist( *np.histogram( clean(datum[vslice,hslice]), bins = bins ) )
                } for slicekey, hslice in hslices.items()
            })
            for datakey, datum in data.items() 
        ])
        #if first: 
            #header += ['pixelcount']
            #fmt += ['%d']
        #line += [ len(clean(data)) ]

        
        print 'keys', [ (datakey,hist[datakey].keys()) for datakey in data.keys() ]

        for datakey, data in hist.items():
            for slicekey, datum in data.items():
                x = datum['hist']['x']
                y = datum['hist']['y']
                
                datum['fits'] = {}
                datum['fits']['log(Gaussian(mu,sigma)*A)'] = fit( x, y, fitfunc['log(Gaussian(mu,sigma)*A)'], log=True, labels = ['mu','sigma','A'] )
                print datum['fits']['log(Gaussian(mu,sigma)*A)']['params']
                datum['fits']['log(Gaussian(0,sigma)*A)'] = fit( x, y, fitfunc['log(Gaussian(0,sigma)*A)'], sel=[x<0], log=True, labels = ['sigma','A'], p0 = datum['fits']['log(Gaussian(mu,sigma)*A)']['params'].values()[1:] )
                datum['fits']['G($\mu$,$\sigma$)'] = fit( x, y, fitfunc['Gaussian(mu,sigma)*A'], labels = ['mu', 'sigma','A'], p0 = datum['fits']['log(Gaussian(mu,sigma)*A)']['params'].values() ) 
                
                datum['fits'][r'G(0,$\sigma$)'] = fit( x, y, fitfunc['Gaussian(0,sigma)*A'], sel=[x<0], labels = ['sigma','A'], p0 = datum['fits']['log(Gaussian(0,sigma)*A)']['params'].values() ) 
                
                datum['fits'][r'G($\mu$,$\sigma$)*P($\lambda$)'] = fit( x, y, fitfunc['convolution(Gaussian(mu,sigma),Poisson(-lamb,lamb))*A'],  labels = ['mu','sigma','A','lambda'], p0 = datum['fits'][r'G($\mu$,$\sigma$)']['params'].values()+[2.] ) 

                #datum['fits'][r'G(0,$\sigma$)*P($\lambda$)'] = fit( x, y, fitfunc['convolution(Gaussian(0,sigma),Poisson(-lamb,lamb))*A'],  labels = ['sigma','A','lambda'], p0 = datum['fits'][r'G(0,$\sigma$)']['params'].values()+[2.] ) 
                
                for fkey, fval in datum['fits'].items():
                    fval['deconvolve'] = scipy.signal.deconvolve( y, fval['hist'] )[1]
                    
        #printDict( hist, skip=['hist','os', 'osi', 'line', 'hist/%s'%f3 ] )
        #printDictTable( hist, skip=['hist','os', 'osi', 'line', 'hist/%s'%f3 ] )
        
        if plot:
            if osi: del hist['osi']
            print 'initiate figure'
            fig = plt.figure()
            fig.set_size_inches( np.array(fig.get_size_inches())*[len(hist.keys()) ,1] )
            fig.suptitle('runID%s'%runID+' ohdu%s'%ohdu)
            first = True
            if 1:
                gs = GridSpec(5 if fft else 4, len(hist.keys()) )
                ax = { key: fig.add_subplot(gs[:-2 if fft else -1,i]) for i, key in enumerate(hist.keys()) }
                axdiff = { key: fig.add_subplot(gs[-2 if fft else -1,i], sharex=ax[hist.keys()[0]] ) for i, key in enumerate(hist.keys()) }
                #if fft: axfft = [ fig.add_subplot(gs[-1,i]) for i in irange ]
            for ax_ in ax.values(): 
                plt.setp(ax_.get_xticklabels(), visible=False)
        	
            print 'generating plots'
            for key, data in hist.items():
                for slicekey, datum in data.items():
                    if slicekey is 'os': continue
                    
                    ax[key].step( datum['hist']['x'], datum['hist']['y'], where='mid', label = slicekey, lw = 2 )
                    
                    for fitfunc in [
                        #r'G(0,$\sigma$)', 
                        r'G($\mu$,$\sigma$)', 
                        #r'G(0,$\sigma$)*P($\lambda$)', 
                        r'G($\mu$,$\sigma$)*P($\lambda$)']:
                        label = fitfunc
                        label = label.replace(r'$\sigma$', '%.2f'%datum['fits'][fitfunc]['params']['sigma'] )
                        label = label.replace(r'$\lambda$', '' if not 'lambda' in datum['fits'][fitfunc]['params'] else '%.2f'%datum['fits'][fitfunc]['params']['lambda'] ) 
                        label += ' $\chi^2$=%d'%datum['fits'][fitfunc]['chisqr']
                        ax[key].plot( datum['hist']['x'], datum['fits'][fitfunc]['hist'], 
                                     label = label, 
                                     ls = '--' )
                        axdiff[key].step( datum['hist']['x'], datum['hist']['y']/datum['fits'][fitfunc]['hist'], where='mid', label = label )
                        #ax[key].step( datum['hist']['x'], np.abs(datum['fits'][fitfunc]['deconvolve']), where='mid', label = 'dec '+fitfunc )

                    #lamb = 5.
                    #filter_ = cPoisson(datum['hist']['x'], -lamb, lamb)
                    #ax[key].plot( datum['hist']['x'], filter_*np.max(datum['hist']['y']), label = 'poisson', ls = '--' )
        if first: 
            header += ['mu', 'sigma', 'A', 'lambda']
            fmt += ['%.2e', '%.2e', '%.2e', '%.2e']
        line += list(hist['scn']['ac']['fits'][r'G($\mu$,$\sigma$)*P($\lambda$)']['params'].values())
        if first: 
            header += ['chisqr']
            fmt += ['%.2e']
        line += [ hist['scn']['ac']['fits'][r'G($\mu$,$\sigma$)*P($\lambda$)']['chisqr'] ]
        result += [ line ]

        if fft:
            print 'computing Fourier transform'
            for i in irange:
                maskAC = histAClist[i]>countCut
                fftAC = np.fft.fftshift( np.fft.fft( histAClist[i][maskAC] ) )
                fftgAC2neg = np.fft.fftshift( np.fft.fft( gAC2neg[i][maskAC] ) )
                
                axfft[i].step( xbins[maskAC], np.abs(( fftAC/fftgAC2neg )), where='mid', label = 'fftac' )
                #axfft[i].step( xbins[maskAC], np.angle(( fftAC/fftgACneg )), where='mid', label = 'fftac' )
                
                fftpoisson = lambda x,mu,c: np.abs( np.fft.fftshift( np.fft.fft( poisson(x,mu,c))) )
                fftACpp, fftACchi2 = fit( xbins[maskAC], np.abs( fftAC/fftgAC2neg), fftpoisson )
                #print 'fftfit', fftACpp, fftACchi2
                
                #axfft[i].step( xbins[maskAC], np.abs( fftpoisson(xbins[maskAC], 10, 10) ), where='mid', label = 'fftgEx' )
                #axfft[i].step( xbins[maskAC], np.abs( fftpoisson(xbins[maskAC], *fftACpp) ), where='mid', label = 'fftgac' )
                
                maskOS = histOSlist[i]>countCut
                fftOS = np.fft.fftshift( np.fft.fft( histOSlist[i][maskOS] ) )
                fftgOS = np.fft.fftshift( np.fft.fft( gOS[i][maskOS] ) )
                axfft[i].step( xbins[maskOS], np.abs((fftOS/fftgOS)), where='mid', label = 'fftos' )
                #axfft[i].step( xbins[maskOS], np.angle((fftOS/fftgOS)), where='mid', label = 'fftos' )
                maxfreq = xbins[maskAC][ np.abs(fftAC) == np.max(np.abs(fftAC)) ]
                print maxfreq
        first = False
    if plot:
        print 'setting labels and lims'
        for key in hist.keys():
            axdiff[key].set_xlabel('pixel energy [adu]')
            ax[key].set_ylabel('pixel count')
            ax[key].set_yscale('log')
            ax[key].set_title( key )
            if fft: 
                axfft[key].set_yscale('log')
                axfft[key].set_xlim([maxfreq-100,maxfreq+100])
            #ax[key].set_xlim([-120,120])
            ax[key].set_ylim(bottom=countCut)
            ax[key].grid(True)
            axdiff[key].set_yscale('log')
            axdiff[key].grid(True)
            axdiff[key].set_ylim([5e-1,2e0])
            ax[key].legend()
            if fft: axfft[key].legend(loc=2)
            axdiff[key].legend()
        plt.subplots_adjust(hspace=.1, wspace=.25)
        print 'plot generated in PNG file', '%s.png'%(FITSfiles[0])
        fig.savefig('%s.png'%(FITSfiles[0]))
    
    #print 'result'
    #print ', '.join(header)
    #print '\n'.join( map( lambda line:', '.join(map(str,line)), result ))
    #print '\n'.join( map( lambda line:', '.join(map(lambda x: str(type(x)),line)), result ))
    print 'table generated in CSV file', '%s.csv'%(FITSfiles[0])
    np.savetxt('%s.csv'%(FITSfiles[0]), result, header=', '.join(header), fmt=' '.join(fmt), delimiter=', ')
 
if __name__ == "__main__":
    print term.bold('Dark Current Analysis Tool for the CONNIE colaboration')
    print term.bold('by Philipe Mota (philipe.mota@gmail.com)')
    print term.bold('repository https://github.com/PhMota/CONNIEtools/blob/master/fits_reader4.py')
    if len(sys.argv) == 1:
        print 'Usage:'
        print '--runID <runID> --outfolder <path> --removehits'
        print '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in the associated ROOT catalog'
        print
        print '--runID <runID> --outfolder <path> --ROOTfile <path> --removehits'
        print '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in given ROOTfile'
        print 
        print '--ohdu <ohdu> --fft --plot <list of paths>'
        print '\tgenerate plot with FITS files from paths for the given ohdu'
        print
        exit(0)
    runID = None
    outfolder = None
    ohdu = None
    ROOTfile = None
    fft = False
    if '--outfolder' in sys.argv:
        outfolder = sys.argv[sys.argv.index('--outfolder')+1]
        print 'setting outfolder =', outfolder
    if '--ohdu' in sys.argv:
        ohdu = int(sys.argv[sys.argv.index('--ohdu')+1])
        print 'setting ohdu =', ohdu
    if '--runID' in sys.argv:
        runID = int(sys.argv[sys.argv.index('--runID')+1])
        print 'setting runID =', runID
    if '--ROOTfile' in sys.argv:
        ROOTfile = sys.argv[sys.argv.index('--ROOTfile')+1]
        print 'setting ROOTfile =', ROOTfile
    if '--fft' in sys.argv:
        fft = True
        print 'setting fft = True'
        
    if '--removehits' in sys.argv:
        if not runID is None:
            removeHitsFromrunID( runID, outputfolder=outfolder, ohdu = ohdu, ROOTfile=ROOTfile )
            exit(0)
        print term.red('error: runID not set')
    if '--plot' in sys.argv:
        if not ohdu is None:
            if not runID is None:
               plotrunID( ohdu, fft=fft, runID=runID)
               exit(0)
            FITSfiles = sys.argv[sys.argv.index('--plot')+1:]
            print 'input FITSfiles:'
            print '\n'.join(FITSfiles)
            plotrunID( ohdu, FITSfiles, fft=fft )
            exit(0)
        print term.red('error: ohdu not set')
    if '--table' in sys.argv:
        if not ohdu is None:
            FITSfiles = sys.argv[sys.argv.index('--table')+1:]
            print 'input FITSfiles:'
            print '\n'.join(FITSfiles)
            plotrunID( ohdu, FITSfiles, fft=fft, plot=False )
            exit(0)
        print term.red('error: ohdu not set')
