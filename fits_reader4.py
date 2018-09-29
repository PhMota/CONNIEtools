# coding: utf-8

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

from functools import partial

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

def rglob( pattern ):
    print 'trying pattern', pattern
    return rglob( pattern.replace('//','//*/') ) if len(glob.glob( pattern )) == 0 else glob.glob( pattern )

def listFITS( *patterns ):
    '''
    list all files that match the list of patterns given
    '''
    return sortByrunID([ match for pattern in patterns for match in rglob(pattern) if not '-with-' in match ])

def listFITS_run( *runs ):
    return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/scn/merged/scn_*.fits'%run for run in runs ] )

def listrunID_run( *runs ):
    return map( getrunIDFromPath, listFITS_run( runs ) )

def listFITS_runID( *runIDs ):
    return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/scn_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )

def listremovedFITS_runID( *runIDs, **kwargs ):
    #if os.access( '/share/storage2/connie/data_analysis/processed02_data/runs/', os.W_OK):
        #print 'input path', '/share/storage2/connie/data_analysis/processed02_data/runs'
        #return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/*hsub_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    inputpath = '.'
    if 'inputpath' in kwargs:
        inputpath = kwargs['inputpath']
    print 'input path', inputpath
    return listFITS( *[ inputpath.rstrip('/') + '//*hsub_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
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
    #sys.stdout.write( msg + '...' )
    #sys.stdout.flush()
    print 'begin:', msg
    if not eta is None and not N is None:
        etastart = time.time()
        eta()
        et = time.time() - etastart
        print 'estimated time', time.strftime("%H:%M:%S", time.gmtime(int(et*N))), '(finishes at %s)'% time.strftime("%H:%M:%S", time.localtime( time.time() + int(et*N)))
    start = time.time()
    res = cmd()
    #print 'elapsed time', time.strftime("%H:%M:%S", time.gmtime(datetime.datetime.time() - start))
    print 'end:', msg, '[t='+time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start))+']'
    return res

def readCatalog( ROOTfile, runID = None, verbose = False ):
    min_max = lambda x: (np.min(x), np.max(x))
    print root_numpy.list_branches( ROOTfile, treename = 'config' )
    #print [ (branch, root_numpy.root2array( ROOTfile, treename = 'config' )[branch]) for branch in root_numpy.list_branches( ROOTfile, treename = 'config' ) ]
    if verbose: print 'reading start and stop indexes for runID', runID, 'in ROOT file'
    indexlist = np.argwhere( root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'] == runID ).flatten()
    start, stop = min_max( indexlist )
    if verbose: print 'reading', len(indexlist), 'entries from', start, 'to', stop
    
    readcatalog_once = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu'], selection='runID==%s'%runID, start=start_, stop=stop_+1 )
    #if verbose: print 'reading ROOT file'
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

def removeHitsFromrunID( runID, outputfolder = None, ohdu = None, ROOTfile = None, osi = None, verbose = False, image = True, output = True, dohits = True ):
    FITSfile = listFITS_runID( runID )[0]
    print 'input FITSfile =', FITSfile
    print 'output flag =', output
    if output:
        outputfilenohits = 'hsub_'+os.path.basename(FITSfile)
        if dohits: outputfilehits = 'hits_'+os.path.basename(FITSfile)
    
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
    
    if dohits:
        hdulist['hits'] = astropy.io.fits.open(FITSfile)
        if not ohdu is None:
            removeOHDU( hdulist['hits'], ohdu )

    hdulist['nohits'] = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        removeOHDU( hdulist['nohits'], ohdu )
    
    if ROOTfile is None:
        ROOTfile = getAssociatedCatalog( FITSfile )[0]
    print 'input ROOTfile =', ROOTfile
    hitscatalog = readCatalog( ROOTfile, runID=runID)

    if verbose: print 'removing hits'
    
    for iohdu in range(len(hdulist['nohits'])):
        ohdu_ = hdulist['nohits'][iohdu].header['OHDU']
        print 'removing hits from', ohdu_
        hitsohdu = hitscatalog[ hitscatalog['ohdu']==ohdu_ ]
        if dohits: hdulist['hits'][iohdu].data[:,:] = REMOVE #-1e9
        hits_xy = np.array( [ [iy,ix] for x,y in zip( hitsohdu['xPix'], hitsohdu['yPix']) for ix,iy in zip(x,y) ] )
        hits_e = np.array( [ ie for e in hitsohdu['ePix'] for ie in e ] )
        if len(hits_e) > 0:
            hdulist['nohits'][iohdu].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE
            if osi: merged[iohdu].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE
            if dohits: hdulist['hits'][iohdu].data[hits_xy[:,0],hits_xy[:,1]] = hits_e
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

    if output:
        if dohits:
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
    
    return hdulist

def fit( x, data, func, sel = None, log=False, p0 = None, labels = None ):
    mask = data > 0
    if sel: mask = np.all( sel+[ data>0 ], axis=0 )
    pp = scipy.optimize.curve_fit( func, x[mask], data[mask] if not log else np.log(data[mask]), p0, bounds=(0,np.inf) )[0]
    hist = func( x, *pp ) if not log else np.exp(func( x, *pp ))
    chisq = scipy.stats.chisquare( hist, f_exp = data, ddof = len(labels) )[0] #/(len(data) - len(pp))
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
    if type(x) is float:
        x -= mu
        return np.abs( np.exp(x*np.log(lamb) -lamb -scipy.special.loggamma(x+1)) ) if x>=0 else 0
    res = np.zeros_like(x)
    x_ = x - mu
    res[x_>=0] = np.abs( np.exp(x_[x_>=0]*np.log(np.abs(lamb)) -np.abs(lamb) -scipy.special.loggamma(x_[x_>=0]+1)) )
    return res

def poisson( x, mu, lamb ):
    if x < mu: return 0
    x -= mu
    return np.exp(x*np.log(lamb) -lamb -scipy.special.loggamma(x+1))
    
def convolution_GP( x, mu = 0, sigma = 1, A = 1, lamb = 0, Nmax = 500 ):
    lamb = abs(lamb)
    mu = abs(mu)
    #mu += lamb
    #sigma -= lamb
    sp = scipy.stats.norm.pdf(x,mu,sigma)
    if lamb == 0: return A*sp
    sp0 = np.sum( sp )
    for i in range(1,Nmax+1):
        spi = lamb**i*scipy.stats.norm.pdf(x-i,mu,sigma)/factorial(i)
        sp += spi
        if np.sum(spi)/sp0 < 1e-4: break
    if i == Nmax-1: print 'Warning: Nmax reached in convolution'
    return A*np.exp(-lamb)*sp

def convolution_GP2( x, Nmax = 500, **kwargs ):#mu_os, sigma_os, mu_ac, sigma_ac, A, lamb, Nmax = 500 ):
#def convolution_GP2( x, mu_os, sigma_os, mu_ac, sigma_ac, A, lamb, Nmax = 500 ):
    mu_os = kwargs['mu_os']
    mu_ac = kwargs['mu_ac']
    sigma_os = kwargs['sigma_os']
    sigma_ac = kwargs['sigma_ac']
    A = kwargs['A']
    lamb = kwargs['lambda']
    
    gain = kwargs['lambda']
    
    mu = mu_os+mu_ac-lamb
    sigma = np.sqrt( sigma_os**2 + sigma_ac**2 )
    sp = scipy.stats.norm.pdf(x,mu,sigma)
    if lamb == 0.0: return A*sp
    sp0 = np.sum( sp )
    for i in range(1,Nmax+1):
        spi = lamb**i*scipy.stats.norm.pdf(x-float(i)/gain,mu,sigma)/factorial(i)
        sp += spi
        if np.sum(spi)/sp0 < 1e-4: break
    if i == Nmax-1: print 'Warning: Nmax reached in convolution'
    return A*np.exp(-lamb)*sp
    
def fitpartial( x, data, sel = None, log=False, p0 = None, adjust = None, bounds = (0.,np.inf) ):
    fix = filter( lambda k: not k[0] in adjust, p0.items() )
    p_ = [ v for k,v in p0.items() if k in adjust ]
    mask = data > 0
    if sel: mask = np.all( sel+[ data>0 ], axis=0 )
    
    x_ = x[mask]
    if log:
        f = lambda x_, *p: np.log( convolution_GP2(x_, **{ k: ( p[ adjust.index(k) ] if k in adjust else v ) for k,v in p0.items() } ) )
        y_ = np.log( data[mask] )
    else:
        f = lambda x_, *p: convolution_GP2(x_, **{ k: ( p[ adjust.index(k) ] if k in adjust else v ) for k,v in p0.items() } )
        y_ = data[mask]
    pp = scipy.optimize.curve_fit( f, x_, y_, p_, bounds=bounds )[0]
    hist = f( x, *pp ) if not log else np.exp(f( x, *pp ))
    chisq = scipy.stats.chisquare( hist, f_exp = data )[0]/(len(data) - len(pp))
    pp = OrderedDict([ (label, pp[adjust.index(label)] if label in adjust else val) for label, val in p0.items() ] )
    return {'chisqr': chisq, 'params': pp, 'hist': hist }

def computeFits( ohdu, runID, hdulists, verbose=False ):
    result = []
    header = []
    fmt = []
    first = True

    clean = lambda x: x.flatten()[ np.all( [x.flatten()!=REMOVEOLD, x.flatten()!=REMOVE, x.flatten()!=1e10], axis=0) ]

    hdulists = OrderedDict([('scn', hdulists)])
    if verbose: print 'ohdus in file', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulists['scn'] ] ))
    
    data = OrderedDict([ (key, hdu.data) for key, hdulist in hdulists.items() for hdu in hdulist if hdu.header['OHDU'] == ohdu ])
    if verbose: 
        print 'keys', data.keys()
        print 'shapes', data['scn'].shape
        print 'data extremes', data['scn'].min(), data['scn'].max()
    yccd = 4262
    xccd = 4220
    yos = 4112
    #if osi:
        #if data['osi'].shape[1] == 2*yccd:
            #print 'subtracting idle amplifier'
            #data['osisub'] = data['osi'][0:xccd,0:yccd] - np.flip(data['osi'][0:xccd,yccd:2*yccd], 1)
    #if scnsub:
        #print 'subtracting idle overscan'
        #print 'bias correction shape', np.median(data['scn'][0:xccd,yos:], axis=1).shape
        #data['scnsub'] = data['scn'][0:xccd,:] - np.median(data['scn'][0:xccd,yos+10:-10], axis=1)[:,np.newaxis]
    #print 'data len', data.keys()
    dx = 1
    shift = lambda bins: (bins[1:]+bins[:-1])*.5
    #bins = np.r_[min(clean(data['scn'])):max(clean(data['scn'])):dx]
    bins = np.r_[-2*max(clean(data['scn'])):max(clean(data['scn'])):dx]
    if verbose: print 'bins', min(bins), max(bins)

    if verbose: print 'computing histogram of full ccd'
    vslice = slice(120,4180)
    hslices = OrderedDict( [
        ('os', slice(4112+10, 4262-10)), 
        ('ac', slice(20,4100)), 
        #('ac_os', slice(3112+10,3262-10))
        ] )
    countCut = 1e1
    parseHist = lambda h,x: { 'y': h[h/dx>countCut]/dx, 'x': (.5*(x[1:]+x[:-1]))[h/dx>countCut] }
    #parseHist = lambda h,x: { 'y': h, 'x': (.5*(x[1:]+x[:-1])) }
    
    hist = OrderedDict( [ (datakey, { 
        #'all': {'hist': np.histogram( clean(datum[0:xccd,0:yccd]), bins = bins )[0] }, 'slice': '0:%s|0:%s'%(xccd,yccd),
        slicekey: {
            'hist': parseHist( *np.histogram( clean(datum[vslice,hslice]), bins = bins ) ),
            'pixelcount': len(clean(datum[vslice,hslice]))
            } for slicekey, hslice in hslices.items()
        })
        for datakey, datum in data.items() 
    ])
    if verbose: print 'keys', [ (datakey, hist[datakey].keys()) for datakey in data.keys() ]

    def change( d, k, v ):
        d0 = dict(d)
        d0[k] = v
        return d0
    
    first=True
    fitindex = 0
    for datakey, data in hist.items():
        for slicekey in hslices.keys():
            datum = data[slicekey]
            x = datum['hist']['x']
            y = datum['hist']['y']
            N = y.sum()/dx
            #N_ = 2*y[x<0].sum()/dx
            N_ = N
            print 'µ', 'σ', 'λ'
            print datum['pixelcount']
            print y.sum()/dx
            print 2*y[x<0].sum()/dx
            
            datum['fits'] = OrderedDict({})

            p0 = OrderedDict( [ 
                ('mu_os', 0.), 
                ('sigma_os', 1. if slicekey=='os' else 0), 
                ('mu_ac', 0.), 
                ('sigma_ac', 1. if slicekey=='ac' else 0), 
                ('A', N), 
                ('lambda',0.), 
                #('gain', 1.),
                ] )
            print slicekey
            datum['fits'][r'G(µ,σ)'] = fitpartial( x, y, p0=p0, adjust = ['mu_%s'%slicekey, 'sigma_%s'%slicekey] )
            datum['fits'][r'G(µ,σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(µ,σ)']['params'], adjust = ['mu_%s'%slicekey, 'sigma_%s'%slicekey, 'lambda'] )

            datum['fits'][r'G(0,σ)'] = fitpartial( x, y, p0=p0, adjust = ['sigma_%s'%slicekey] )
            datum['fits'][r'G(0,σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)']['params'], adjust = ['sigma_%s'%slicekey, 'lambda'] )
            datum['fits'][r'G(0,σ=σ0)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)']['params'], adjust = ['lambda'] )
            datum['fits'][r'G(0,σ0-)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)']['params'], sel=[x<0], adjust = ['sigma_%s'%slicekey] )
            datum['fits'][r'G(0,σ=σ0-)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ0-)']['params'], adjust = ['lambda'] )
            
            datum['fits'][r'log(G(µ,σ))'] = fitpartial( x, y, p0=datum['fits'][r'G(µ,σ)']['params'], adjust = ['mu_%s'%slicekey, 'sigma_%s'%slicekey], log = True )
            datum['fits'][r'log(G(µ,σ)*P(λ))'] = fitpartial( x, y, p0=datum['fits'][r'G(µ,σ)*P(λ)']['params'], adjust = ['mu_%s'%slicekey, 'sigma_%s'%slicekey, 'lambda'], log = True )
            
            datum['fits'][r'log(G(0,σ))'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)']['params'], adjust = ['sigma_%s'%slicekey], log = True )
            datum['fits'][r'log(G(0,σ)*P(λ))'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)*P(λ)']['params'], adjust = ['sigma_%s'%slicekey, 'lambda'], log = True )
            #datum['fits'][r'logA*G(0,σ)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)']['params'], adjust = ['sigma_%s'%slicekey, 'A'], log = True )
            #datum['fits'][r'logG(0,σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)']['params'], adjust = ['sigma_%s'%slicekey, 'lambda'], log = True )
            #datum['fits'][r'logA*G(0,σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(0,σ)*P(λ)']['params'], adjust = ['sigma_%s'%slicekey, 'A', 'lambda'], log = True )

            #if slicekey == 'ac':
                #p0 = change( data['os']['fits'][r'G(µ,σ)']['params'], 'A', N )
                #p2 = change( p0, 'sigma_ac', 1. )
                #p2a = change(change( p2, 'sigma_os', 0. ), 'mu_os', 0)
                
                #datum['fits'][r'G(µos+µ,σos+σ)'] = fitpartial( x, y, p0=p2, adjust = ['mu_ac', 'sigma_ac' ] )
                #datum['fits'][r'G(µos+µ,σos+σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(µos+µ,σos+σ)']['params'], adjust = ['mu_ac', 'sigma_ac','lambda'] )
                #datum['fits'][r'G(µ=µos+µ,σ=σos+σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(µos+µ,σos+σ)']['params'], adjust = ['lambda'] )
                
                #datum['fits'][r'G(µos+0,σos+σ)'] = fitpartial( x, y, p0=p0, adjust = ['sigma_ac'] )
                #datum['fits'][r'G(µos+0,σos+σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(µos+0,σos+σ)']['params'], adjust = ['sigma_ac', 'lambda'] )
                #datum['fits'][r'G(µ=µos+0,σ=σos+σ)*P(λ)'] = fitpartial( x, y, p0=datum['fits'][r'G(µos+0,σos+σ)']['params'], adjust = ['lambda'] )

            #for sigmakey in ['σ0-','σ0','σ0+']:
                #sel = [x<0] if '-' in sigmakey else ( [x>0] if '+' in sigmakey else [] )
                #fkey_s = r'G(µ=0,%s)'%(sigmakey)
                #if slicekey == 'os':
                    #datum['fits'][r'G(µ=0,%s)'%(sigmakey)] = fitpartial( x, y, sel=sel, p0=p0, adjust = ['sigma_os'] )
                    #datum['fits'][r'G(µ=0,σ=%s)*P(λ)'%(sigmakey)] = fitpartial( x, y, p0=datum['fits'][r'G(µ=0,%s)'%(sigmakey)]['params'], adjust = ['lambda'] )
                #if slicekey == 'ac':
                    #p1 = change(change(data['os']['fits'][r'G(µ=0,%s)'%(sigmakey)]['params'], 'sigma_ac', 1.), 'A', N)

                    #datum['fits'][r'G(0+µ,%s+σ)'%(sigmakey)] = fitpartial( x, y, p0=p1, adjust = ['sigma_ac'] )
                    #datum['fits'][r'G(0+µ,%s+σ)*P(λ)'%(sigmakey)] = fitpartial( x, y, p0=datum['fits'][r'G(0+µ,%s+σ)'%(sigmakey)]['params'], adjust = ['mu_ac', 'sigma_ac', 'lambda'] )
                    #datum['fits'][r'G(0+µ,%s+σ)*P(λ)'%(sigmakey)] = fitpartial( x, y, p0=datum['fits'][r'G(0+µ,%s+σ)'%(sigmakey)]['params'], adjust = ['lambda'] )

                    #datum['fits'][r'G(µ=µos,σos+%s)'%(sigmakey)] = fitpartial( x, y, sel=sel, p0=p2, adjust = ['sigma_ac'] )
                    #datum['fits'][r'G(µ=µos,σos+%s)*P(λ)'%(sigmakey)] = fitpartial( x, y, sel=sel, p0=datum['fits'][r'G(µ=µos,σos+%s)'%(sigmakey)]['params'], adjust = ['lambda'] )
                    
                    #for sigma2key in ['σ0-','σ0','σ0+','σ']:
                        #sel2 = [x<0] if '-' in sigma2key else ( [x>0] if '+' in sigma2key else [] )
                        #datum['fits'][r'G(µ=0,%s+%s)'%(sigmakey, sigma2key)] = fitpartial( x, y, sel=sel2, p0=p1, adjust = ['sigma_ac'] )
                        #datum['fits'][r'G(µ=0,%s+%s)*P(λ)'%(sigmakey, sigma2key)] = fitpartial( x, y, p0=datum['fits'][r'G(µ=0,%s+%s)'%(sigmakey, sigma2key)]['params'], adjust = ['lambda'] )
                        
            tablekeys = datum['fits'].keys()
            fitkeylen = max(map(len,tablekeys))
            print fitkeylen
            for fitkey in tablekeys:
                if not fitkey in datum['fits']: continue
                line = []
                if first: 
                    header += ['fitID']
                    fmt += ['%d']
                line += [ fitindex ]
                if first: 
                    header += ['runID']
                    fmt += ['%d']
                line += [ int(runID) ]
                if first: 
                    header += ['ohdu']
                    fmt += ['%d']
                line += [ int(ohdu) ]
                if first: 
                    header += ['tempmin']
                    fmt += ['%.2f']
                line += [ float(hdulists['scn'][0].header['TEMPMIN' ]) ]
                if first: 
                    header += ['tempmax']
                    fmt += ['%.2f']
                line += [ float(hdulists['scn'][0].header['TEMPMAX' ]) ]
                if first: 
                    header += ['file']
                    fmt += ['%.32s']
                line += [ datakey ]
                if first: 
                    header += ['section']
                    fmt += ['%.32s']
                line += [ slicekey ]
                if first: 
                    header += ['pixelcount']
                    fmt += ['%.2e']
                line += [ datum['pixelcount'] ]
                if first:
                    header += ['fitfunc']
                    fmt += ['%%%ds'%fitkeylen]
                line += [ fitkey ]
                params = ['mu_os', 'sigma_os', 'mu_ac', 'sigma_ac', 'A', 'lambda' ]
                if first:
                    header += params
                    fmt += ['% .4e']*7
                line += [ (datum['fits'][fitkey]['params'][k] if k in datum['fits'][fitkey]['params'] else 0) for k in params ]
                if first: 
                    header += ['chisqr']
                    fmt += ['%.3e']
                line += [ datum['fits'][fitkey]['chisqr'] ]
                result += [ line ]
                first = False
                fitindex += 1
    #sep = ', '
    sep = ' '
    print sep.join( [ str((i,h)) for i, h in enumerate(header)] )
    print sep.join(header)
    print '\n'.join( [ sep.join( [ f%entry for f,entry in zip(fmt,line) ] ) for line in result ] )
    printcsv = lambda outfile: np.savetxt('%s.csv'%(outfile), result, header=sep.join(header), fmt='%s', delimiter=sep)
    return hist, printcsv

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
    inputpath = '.'
    if 'inputpath' in kwargs:
        inputpath = kwargs['inputpath']
        
    if runID is None:
        FITSfiles = FITSfiles[0]
    else:
        FITSfiles = listremovedFITS_runID(int(runID),inputpath=inputpath)
    
    verbose = False
    
    print 'FITSfiles:\n', '\n'.join(FITSfiles)
    runID = getrunIDFromPath( FITSfiles[0] )
        
    hist, printcsv = computeFits( ohdu, runID, astropy.io.fits.open( FITSfiles[0] ), verbose )
    
    if plot:
        if osi: del hist['osi']
        print 'initiate figure'
        fig = plt.figure()
        fig.set_size_inches( np.array(fig.get_size_inches())*[len(hslices.keys()) ,len(hist.keys())] )
        fig.suptitle('runID%s'%runID+' ohdu%s'%ohdu)
        first = True
        if 1:
            span = 5 if fft else 4
            
            gs = GridSpec( span*len(hist.keys()), len(hslices.keys()) )
            ax = { key: { hkey: fig.add_subplot(gs[span*(i-1):span*(i-1)+3,j]) for j, hkey in enumerate(hslices.keys()) } for i, key in enumerate(hist.keys()) }
            axdiff = { key: { hkey: fig.add_subplot(gs[span*(i-1)+3,j], sharex=ax[key][hkey] ) for j, hkey in enumerate(hslices.keys())} for i, key in enumerate(hist.keys()) }
            if fft:
                axfft = { key: { hkey: fig.add_subplot(gs[span*(i-1)+4,j] ) for j, hkey in enumerate(hslices.keys())} for i, key in enumerate(hist.keys()) }
        for ax_ in ax.values():
            for ax__ in ax_.values():
                plt.setp(ax__.get_xticklabels(), visible=False)
        
        if verbose: print 'generating plots'
        for key, data in hist.items():
            for slicekey, datum in data.items():
                #if slicekey is 'ac_os' or slicekey is 'os': continue
                
                ax[key][slicekey].step( datum['hist']['x'], datum['hist']['y'], where='mid', label = slicekey, lw = 2 )
                
                for i, fitfunc_ in enumerate([
                    r'G(0,$\sigma$)', 
                    #r'G($\mu$,$\sigma$)', 
                    r'G(0,$\sigma$)*P($\lambda$)',
                    #r'G(0,$\sigma$\')*P($\lambda$)',
                    r'G(0,$\sigma$os)*P($\lambda$)', 
                    #r'G($\mu$,$\sigma$)*P($\lambda$)',
                    #r'sG($\mu$,$\sigma$)*P($\lambda$)',
                    #r'cG($\mu$,$\sigma$)*P($\lambda$)',
                    ]):
                    label = fitfunc_
                    label = label.replace(r'$\sigma$', '%.2f'%datum['fits'][fitfunc_]['params']['sigma'] )
                    label = label.replace(r'$\lambda$', '' if not 'lambda' in datum['fits'][fitfunc_]['params'] else '%.2g'%abs(datum['fits'][fitfunc_]['params']['lambda']) ) 
                    label += ' $\chi^2$=%.3g'%datum['fits'][fitfunc_]['chisqr']
                    ax[key][slicekey].plot( datum['hist']['x'], datum['fits'][fitfunc_]['hist'], 
                                    label = label, 
                                    ls = '--', color='C%d'%(i+1) )
                    axdiff[key][slicekey].step( datum['hist']['x'], datum['hist']['y']/datum['fits'][fitfunc_]['hist'], where='mid', label = label, color='C%d'%(i+1) )
                    #ax[key].step( datum['hist']['x'], np.abs(datum['fits'][fitfunc]['deconvolve']), where='mid', label = 'dec '+fitfunc )
                
                if fft:
                    fourier = lambda x: np.fft.fftshift( np.fft.fft( x ))
                    y = fourier( datum['hist']['y'] )
                    ly = np.log( y )
                    #axfft[key][slicekey].step( datum['hist']['x'], np.imag(ly), where='mid', label='i' )
                    lamb = np.mean( np.real(ly) )
                    print slicekey, lamb
                    axfft[key][slicekey].step( datum['hist']['x'], np.real(y), where='mid', label='r' )
                    #axfft[key][slicekey].step( datum['hist']['x'], np.abs(ly), where='mid', label='a' )
                    #axfft[key][slicekey].step( datum['hist']['x'], fourier(iy), where='mid', label='fiy' )

    if plot:
        if verbose: print 'setting labels and lims'
        for key in hist.keys():
            for hkey in hslices.keys():
                axdiff[key][hkey].set_xlabel('pixel energy [adu]')
                ax[key][hkey].set_ylabel('pixel count')
                ax[key][hkey].set_yscale('log')
                ax[key][hkey].set_title( key )
                if fft: 
                    #axfft[key][hkey].set_yscale('log')
                    #axfft[key][hkey].set_xlim([maxfreq-100,maxfreq+100])
                    #axfft[key][hkey].set_xlim([-5,5])
                    axfft[key][hkey].legend(loc=2)
                #ax[key].set_xlim([-120,120])
                ax[key][hkey].set_ylim(bottom=countCut)
                ax[key][hkey].grid(True)
                #axdiff[key].set_yscale('log')
                axdiff[key][hkey].grid(True)
                axdiff[key][hkey].set_ylim([5e-1,2e0])
                ax[key][hkey].legend()
                #if fft: 
                #axdiff[key][hkey].legend()
        plt.subplots_adjust(hspace=.1, wspace=.25)
        print 'plot generated in PNG file', '%s.png'%(FITSfiles[0])
        fig.savefig('%s.png'%(FITSfiles[0]), bbox_inches='tight')
    
    print 'table generated in CSV file', '%s.csv'%(FITSfiles[0])
    printcsv(FITSfiles[0])
 
def help_():
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

def getOption( vars_, option, args, f = lambda x:x ):
    if '--%s'%option in args:
        vars_[option] = f(args[args.index('--%s'%option)+1])
        print 'setting %s ='%option, vars_[option]
        return
    vars_[option] = None
    return
    
if __name__ == "__main__":
    print term.bold('Dark Current Analysis Tool for the CONNIE colaboration')
    print term.bold('by Philipe Mota (philipe.mota@gmail.com)')
    print term.bold('repository https://github.com/PhMota/CONNIEtools/blob/master/fits_reader4.py')
    if len(sys.argv) == 1:
        help_()
        exit(0)
    vars_ = {}
    outfolder = None
    ROOTfile = None
    fft = False
    if '--outfolder' in sys.argv:
        outfolder = sys.argv[sys.argv.index('--outfolder')+1]
        print 'setting outfolder =', outfolder
    getOption( vars_, 'ohdu', sys.argv, int )
    getOption( vars_, 'runID', sys.argv, int )
    getOption( vars_, 'run', sys.argv )
    getOption( vars_, 'inputpath', sys.argv )
    if vars_['inputpath'] is None:
        vars_['inputpath'] = '.'
        
    if '--ROOTfile' in sys.argv:
        ROOTfile = sys.argv[sys.argv.index('--ROOTfile')+1]
        print 'setting ROOTfile =', ROOTfile
    if '--fft' in sys.argv:
        fft = True
        print 'setting fft = True'
        
    if '--removehits' in sys.argv:
        if not vars_['runID'] is None:
            runETA( 'remove hits',
                lambda: removeHitsFromrunID( vars_['runID'], outputfolder=outfolder, ohdu = vars_['ohdu'], ROOTfile=ROOTfile )
                )
            exit(0)
        print term.red('error: runID not set')
    if '--plot' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA( 'analyse and plot', lambda: plotrunID( vars_['ohdu'], fft=fft, runID=vars_['runID'], inputpath=vars_['inputpath']) )
                exit(0)
            FITSfiles = sys.argv[sys.argv.index('--plot')+1:]
            print 'input FITSfiles:'
            print '\n'.join(FITSfiles)
            plotrunID( vars_['ohdu'], FITSfiles, fft=fft )
            exit(0)
        print term.red('error: ohdu not set')
    if '--table' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA( 'analyse', lambda: plotrunID( vars_['ohdu'], fft=fft, runID=vars_['runID'], plot=False, inputpath=vars_['inputpath']) )
                exit(0)
            FITSfiles = sys.argv[sys.argv.index('--table')+1:]
            print 'input FITSfiles:'
            print '\n'.join(FITSfiles)
            plotrunID( vars_['ohdu'], FITSfiles, fft=fft, plot=False )
            exit(0)
        print term.red('error: ohdu not set')
    if '--listrunID' in sys.argv:
        print 'list of runID'
        print ', '.join( listrunID_run(vars_['run']) )
        exit(0)
    print 'options not recognized', sys.argv
    help_()
    exit(0)
