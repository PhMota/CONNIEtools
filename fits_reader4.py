# coding: utf-8

import astropy
import astropy.io
import astropy.io.fits

import math
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

class Options:
    def setparam( self, arg, value = None ):
        self.__dict__[arg] = value if value else True
        return
    
    def is_option( self, arg ):
        return ('--' == arg[:2])
    
    def __init__( self, argv ):
        self._cmdline = []
        for i,arg in enumerate(argv):
            if self.is_option(arg): 
                self._cmdline += [ arg ]
                self.setparam(arg, (argv[i+1] if not self.is_option(argv[i+1]) else 'error') if i < len(argv)-1 else None )
        return
    
    def get_cmdline(self):
        for _0 in self._cmdline:
            print _0, self.__dict__[_0]

class Program:
    def __init__(self):
        options = Options( sys.argv )
        options.get_cmdline()
        return

class ROOTcatalog:
    pass

class FITSimage:
    pass
    

REMOVE = 0
REMOVEOLD = -1e9

eIonization = 3.745 #eV/e-

odict = lambda x: OrderedDict(x)
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
    y = cmd(33)
    b = cmd(34)
    
    red = staticmethod( lambda x: term.r + x + term.Q )
    bold = staticmethod( lambda x: term.B + x + term.Q )
    italic = staticmethod( lambda x: term.I + x + term.Q )
    green = staticmethod( lambda x: term.g + x + term.Q )
    yellow = staticmethod( lambda x: term.y + x + term.Q )

def str_with_err(value, error):
    if error == 0.:
        return '{0:e}'.format(value)
    digitsE = -int(math.floor(math.log10(error)))
    digitsV = -int(math.floor(math.log10(abs(value))))
    #print 'inside', value, error, digitsE, digitsV
    if digitsE<digitsV:
        digitsV = digitsE
    return "{0:.{2}f}({1:.0f})e{3:+03d}".format(value*10**digitsV, error*10**(digitsE), digitsE-digitsV, -digitsV )
    
#print str_with_err( 0.0001, 0.0000005)
#print str_with_err( 1, .5)
#print str_with_err( 1, 0)
#print str_with_err( 10, 5)
#print str_with_err( 10000, 500)
#exit()
#folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/scn/merged/*'
folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029*/data_*/scn/merged/*'
folder2 = '/share/storage2/connie/data_analysis/processed02_data/runs/009*/data_*/scn/merged/*'

def rglob( pattern, n=0 ):
    if n>20: return []
    #print 'trying pattern', pattern
    return rglob( pattern.replace('//','//*/'), n+1 ) if len(glob.glob( pattern )) == 0 else glob.glob( pattern )

def listFITS( *patterns ):
    '''
    list all files that match the list of patterns given
    '''
    return sortByrunID([ match for pattern in patterns for match in rglob(pattern) if not '-with-' in match ])

def listFITS_run( *runs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*%s[A-Z]/data_*/scn/merged/scn_*.fits'%run for run in runs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*%s/data_*/scn/merged/scn_*.fits'%run for run in runs ] )
    return l

def listrunID_run( *runs ):
    l = listFITS_run( runs )
    #print '\n'.join(l)
    return map( getrunIDFromPath, l )

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

def getAssociatedCatalogGain( *paths ):
    files = getAssociatedCatalog( *paths )
    #print files[0]
    pattern = re.search( r'(data_[0-9]+_to_[0-9]+)', files[0] ).groups()[0]
    #print pattern
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/nu_processing/scripts/ProcCat/*scn_osi_raw_gain_catalog_%s.root'%pattern) if not 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedCatalog( *paths ):
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/ext/catalog/catalog_data_*.root'%getrunFromPath(path)) if not 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedOSI( *paths ):
    return [ osi for path in paths for osi in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/osi/images/osi_runID_*_%s_*.fits'%(getrunFromPath(path), getrunIDFromPath(path)) ) if not '-with-' in osi ]

def runETA( msg, cmd, eta = None, N = None, loop=None, verbose=None ):
    if verbose: print term.yellow('%s <begin>'%msg)
    if not loop is None:
        start = time.time()
        finishtimes = []
        for n, entry in enumerate(loop):
            print term.yellow('entry %s <begin>'%entry)
            etastart = time.time()
            try:
                cmd( entry )
            except ValueError as error:
                print term.red('Warning: %s'%error)
            et = time.time() - etastart
            print term.green('entry %s <end '%entry + 't=%s>'%time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - etastart)) )
            finishtimes += [ time.time() + int(et*(len(loop) - n - 1)) ]
            print term.yellow('estimated time %s'%time.strftime("%H:%M:%S", time.gmtime( int(np.mean(finishtimes)) - time.time() ) ) + '(finishes at %s)'%time.strftime("%H:%M:%S", time.localtime( int(np.mean(finishtimes)) )) )
        print term.green('%s <end '%msg + 't=%s>'%time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start)) )
        return
    if not eta is None and not N is None:
        etastart = time.time()
        eta()
        et = time.time() - etastart
        print term.yellow('estimated time %s'%time.strftime("%H:%M:%S", time.gmtime(int(et*N))) + '(finishes at %s)'%time.strftime("%H:%M:%S", time.localtime( time.time() + int(et*N))) )
    start = time.time()
    res = cmd()
    if verbose: print term.green('%s <end '%msg + 't=%s>'%time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start)) )
    return res

def readCatalog( ROOTfile, runID = None, verbose = False, readGain=True ):
    min_max = lambda x: (np.min(x), np.max(x))
    if verbose: print 'reading start and stop indexes for runID', runID, 'in ROOT file'
    indexlist = np.argwhere( root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'] == runID ).flatten()
    try:
        start, stop = min_max( indexlist )
    except ValueError:
        print '!!!Error'
        print ROOTfile
        print list(set(root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'] ))
        print 'indexlist', indexlist
        raise
    if verbose: print 'reading', len(indexlist), 'entries from', start, 'to', stop
    
    readcatalog_once = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu','level'] + (['gainCu'] if readGain else []), selection='runID==%s'%runID, start=start_, stop=stop_+1 )
    #if verbose: print 'reading ROOT file'
    hitscatalog = runETA( 
        msg = 'reading ROOTfile',
        cmd = lambda: readcatalog_once(start, stop), #read all entries and return 
        verbose = verbose,
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

def removeHitsFromrunID( runID, outputfolder=None, ohdu=None, ROOTfile=None, osi=None, verbose=False, image=True, output=True, dohits=True, gain=None, crosstalk=False, onlyoccupancy=False ):
    FITSfile = listFITS_runID( runID )[0]
    print 'input FITSfile =', FITSfile
    if verbose: print 'output flag =', output
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
    #print [ (_1,_2) for _1,_2 in hdulist['nohits'][-1].header.items() ]
    ohdus = [ hdulist['nohits'][i].header['OHDU'] for i in range(len(hdulist['nohits'])) ]
    print 'ohdus', ohdus
    badOHDU = [11,12]
    #ohdus = [ ohdu for ohdu in ohdus if not ohdu in badOHDU ]
    #print 'ohdus*', ohdus
    if not ohdu is None:
        removeOHDU( hdulist['nohits'], ohdu )

    if crosstalk:
        hdulist['nohitsall'] = astropy.io.fits.open(FITSfile)
        if not ohdu is None:
            removeOHDU( hdulist['nohitsall'], ohdu )
    
    readGain = True
    if ROOTfile is None:
        ROOTfile = getAssociatedCatalog( FITSfile )[0]
        try:
            ROOTfileGain = getAssociatedCatalogGain( FITSfile )[0]
            print ROOTfileGain
        except:
            readGain = False
            ROOTfileGain = ROOTfile
    
    ROOTfile = ROOTfileGain
    #print root_numpy.list_branches( ROOTfileGain, treename = 'hitSumm' )
    #print root_numpy.root2array( ROOTfileGain, treename = 'hitSumm', branches = ['gainCu'] )
    print 'input ROOTfile =', ROOTfile
    hitscatalog = readCatalog( ROOTfile, runID=runID, readGain=readGain, verbose=verbose)

    if verbose: print 'removing hits'
    
    levelcut = 0
    hitsmask = {}
    gain = {}
    occupancy = {}
    #if crosstalk:
    for iohdu, thisohdu in enumerate( ohdus ):
        print 'mask ohdu', thisohdu, 
        hitsmask[thisohdu] = {}
        occupancy[thisohdu] = {}
        hitsohdu = hitscatalog[ hitscatalog['ohdu'] == thisohdu ]
        if len(hitsohdu['xPix']) > 0:
            x = np.concatenate( hitsohdu['xPix'], axis=0)
            y = np.concatenate( hitsohdu['yPix'], axis=0)
            l = np.concatenate( hitsohdu['level'], axis=0)
            if readGain:
                gain[thisohdu] = np.mean( hitsohdu['gainCu'] )
                print 'gain', gain[thisohdu], 'adu/keV'
            for level in range(4):
                N = 0
                rate = 0
                totalN = (l<=level).sum()
                hitsmask[thisohdu][level] = np.zeros_like( hdulist['hits'][-1].data, dtype=bool )
                hitsmask[thisohdu][level][y[l<=level],x[l<=level]] = True
                N = hitsmask[thisohdu][level].sum()
                rate = float(N)/len(hitsmask[thisohdu][level].flatten())
                #print totalN, N
                occupancy[thisohdu][level] = {'N': N, 'rate': rate, 'gain': gain[thisohdu], 'totalN': totalN }
        else:
            print 'ohdu', thisohdu, 'empty'
    sep=' '
    np.savetxt('%s/runID_%s_occupancy.csv'%(outputfolder,runID), [ [ ohdukey, levelkey, data['N'], data['rate'] ] for ohdukey, thisohdu in occupancy.items() for levelkey, data in thisohdu.iteritems() ], header=sep.join(['ohdu','level','N','rate']), fmt='%s', delimiter=sep)
    
    if onlyoccupancy: return None
    
    for iohdu in range(len(hdulist['nohits'])):
        ohdu_ = hdulist['nohits'][iohdu].header['OHDU']
        print 'crosstalk subtraction'
        for iohdu2, thisohdu in enumerate( [ ohdu for ohdu in ohdus if not ohdu in badOHDU] ): #[:4]
            print ohdu_, thisohdu, '%e'%np.mean(hdulist['nohits'][iohdu].data[ hitsmask[thisohdu][0] ])
    
    for iohdu in range(len(hdulist['nohits'])):
        ohdu_ = hdulist['nohits'][iohdu].header['OHDU']
        print 'from ohdu', ohdu_, 'remove'
        if verbose: print 'removing hits from ohdus', ohdu_
        hitsohdu = hitscatalog[ hitscatalog['ohdu'] == ohdu_ ]
        if dohits: hdulist['hits'][iohdu].data[:,:] = REMOVE #-1e9

        hdulist['nohits'][iohdu].data *= 1e3/gain[ohdu_]/eIonization #e- (electron unit)
        
        hdulist['nohits'][iohdu].data[ hitsmask[ohdu_][3] ] = REMOVE
        if osi: merged[iohdu].data[ hitsmask[ohdu_][3] ] = REMOVE
        
        #if dohits: hdulist['hits'][iohdu].data[ hitsmask[ohdu_] ] = hits_e*1e3/gain/eIonization
        
        #if crosstalk:
            #print 'crosstalk subtraction'
            #for iohdu2, thisohdu in enumerate( [ ohdu for ohdu in ohdus if not ohdu in badOHDU] ): #[:4]
                #print 'from ohdu', ohdu_, 'removing', thisohdu
                #hdulist['nohitsall'][iohdu].data[ hitsmask[thisohdu] ] = REMOVE            

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
    if crosstalk:
        return hdulist['nohitsall']
    else:
        return hdulist['nohits']

def fit( x, data, func, sel = None, log=False, p0 = None, labels = None ):
    mask = data > 0
    if sel: mask = np.all( sel+[ data>0 ], axis=0 )
    pp = scipy.optimize.curve_fit( func, x[mask], data[mask] if not log else np.log(data[mask]), p0, bounds=(0,np.inf) )[0]
    hist = func( x, *pp ) if not log else np.exp(func( x, *pp ))
    chisq = scipy.stats.chisquare( data, f_exp = hist, ddof = len(labels) )[0]#/(len(data) - len(pp))
    if labels:
        return {'chisqr': chisq, 'params': odict([ [label, pp_] for label, pp_ in zip(labels, pp) ] ), 'hist': hist }
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

def convolution_GP2( x, Nmax = 500, eps = 1.e-4, **kwargs ):#mu_os, sigma_os, mu_ac, sigma_ac, A, lamb, Nmax = 500 ):
#def convolution_GP2( x, mu_os, sigma_os, mu_ac, sigma_ac, A, lamb, Nmax = 500 ):
    mu_os = kwargs['mu_os']
    mu_ac = kwargs['mu_ac']
    sigma_os = kwargs['sigma_os']
    sigma_ac = kwargs['sigma_ac']
    A = kwargs['A']
    lamb = kwargs['lamb']
    #gain = kwargs['gain']
    #gain = kwargs['lambda']
    
    #print 'lambda', lamb
    mu = mu_os+mu_ac - lamb
    sigma = np.sqrt( sigma_os**2 + sigma_ac**2 )
    sp = scipy.stats.norm.pdf(x,mu,sigma)
    if lamb == 0.0: return A*sp
    sp0 = np.sum( sp )
    for i_ in range(1,Nmax+1):
        #spi = lamb**i*scipy.stats.norm.pdf(x-float(i), mu, sigma)/factorial(i)
        spi = np.exp(i_*np.log(lamb) - np.log(factorial(i_)))*scipy.stats.norm.pdf(x-float(i_), mu, sigma)
        sp += np.real(spi)
        if np.real(spi.sum())/sp0 < eps: break
    if i_ == Nmax-1: print 'Warning: Nmax reached in convolution'
    return A*np.exp(-lamb)*sp
    
def fitpartial( x, data, sel = None, log=False, p0 = None, adjust = None, bounds = (0.,np.inf), **kwargs ):
    #fix = filter( lambda k: not k[0] in adjust, p0.items() )
    p0_ = odict(p0)
    for key, value in kwargs.items():
        if key in p0_: p0_[key] = value
    
    p_ = [ v for k,v in p0_.items() if k in adjust ]
    #print p_
    mask = data > 0
    if sel: mask = np.all( sel+[ data>0 ], axis=0 )
    
    x_ = x[mask]
    if log:
        f = lambda x__, *p: np.log( convolution_GP2(x__, **{ k: ( p[ adjust.index(k) ] if k in adjust else v ) for k,v in p0_.items() } ) )
        y_ = np.log( data[mask] )
    else:
        f = lambda x__, *p: convolution_GP2(x__, **{ k: ( p[ adjust.index(k) ] if k in adjust else v ) for k,v in p0_.items() } )
        y_ = data[mask]
        
    #pp, pcov = scipy.optimize.curve_fit( f, x_, y_, p_, sigma = y_, bounds=bounds )
    ##pp, pcov = scipy.optimize.curve_fit( f, x_, y_, p_, sigma = np.sqrt(y_), bounds=bounds )
    ##pp, pcov = scipy.optimize.curve_fit( f, x_, y_, p_, sigma = np.ones_like(y_), bounds=bounds )
    #error = []
    #for i in range(len(pp)):
        #try:
          #error.append(np.absolute(pcov[i][i])**0.5)
        #except:
          #error.append( 0.00 )
    #perr = np.array(error)
    #hist2 = f( x, *pp )
    #chisq2 = scipy.stats.chisquare( hist2, data, ddof = len(pp) )[0]/(len(data) - len(pp))
    
    #residual = lambda p, x_, y_: (y_ - f(x_, *p))**2/y_
    #residual = lambda p, x_, y_: (y_ - f(x_, *p))**2
    residual = lambda p, x__, y__: (y__ - f(x__, *p))
    res = scipy.optimize.leastsq( residual, p_, args=(x_,y_), full_output=1 )
    popt, pcov2, infodict, errmsg, success = res
    #hist3 = f( x, *popt )
    #chisq3 = scipy.stats.chisquare( hist3, data, ddof = len(popt) )[0]/(len(data) - len(popt))
    hist3 = f( x_, *popt )
    chisq3 = scipy.stats.chisquare( hist3, y_, ddof = len(popt) )[0]/(len(y_) - len(popt))
    
    #print 'chi2', chisq3, np.sum( (hist3-y_)**2/y_ )/(len(y_) - len(popt))
    #print np.max(hist3), np.max(y_)
    
    if (len(data) > len(popt)) and pcov2 is not None:
        s_sq = ( residual(popt, x_, y_ )**2).sum()/(len(y_)-len(popt))
        pcov2 = pcov2 * s_sq
    else:
        pcov2 = np.inf
        
    error2 = []
    for i in range(len(popt)):
        try:
          error2.append(np.absolute(pcov2[i][i])**0.5)
        except:
          error2.append( 0.00 )
    perr2 = np.array(error2)
    
    #print 'CF', pp, chisq2, perr
    #print 'LS', popt, scipy.stats.chisquare( f(x, *popt), data, ddof = len(popt) )[0]/(len(data) - len(popt)), perr2
    perr = perr2
    pp = popt
    chisq2 = chisq3
    hist2 = hist3
    
    #minresidual = lambda p, x_, y_: np.sum( (y_ - f(x_, *p))**2 )
    #poptmin = scipy.optimize.minimize( minresidual, p_, args=(x,data) ).x
    
    #hist = f( x, *pp ) if not log else np.exp(f( x, *pp ))
    #chisq = scipy.stats.chisquare( data, hist, ddof = len(pp) )[0]/(len(data) - len(pp))
    #pp = odict([ (label, pp[adjust.index(label)] if label in adjust else val) for label, val in p0_.items() ] )
    popt = odict([ (label, pp[adjust.index(label)] if label in adjust else val) for label, val in p0_.items() ] )
    perr = odict([ (label, perr[adjust.index(label)] if label in adjust else 0) for label, val in p0_.items() ] )
    #print pp
    #print popt
    return {'chisqr': chisq2, 'params': popt, 'hist': hist2, 'perr': perr }

def computeFits( ohdu, runID, hdulists, verbose=False, plot=False, gain=None ):
    result = []
    header = []
    fmt = []
    first = True
    
    clean = lambda x: x.flatten()[ np.all( [x.flatten()!=REMOVEOLD, x.flatten()!=REMOVE, x.flatten()!=1e10], axis=0) ]

    hdulists = odict( [('scn', hdulists)] )
    if verbose: print 'ohdus in file', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulists['scn'] ] ))
    
    data = odict([ (key, hdu.data) for key, hdulist in hdulists.items() for hdu in hdulist if hdu.header['OHDU'] == ohdu ])
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
    N = 500
    shift = lambda bins: ( bins[1:] + bins[:-1] )*.5
    #bins = np.r_[min(clean(data['scn'])):max(clean(data['scn'])):dx]
    maxSCN = max(clean(data['scn']))
    bins = np.linspace(-maxSCN, maxSCN, N)
    dx = 2*maxSCN/N
    if verbose: print 'dx', dx
    if verbose: print 'bins', min(bins), max(bins)

    if verbose: print 'computing histogram of full ccd'
    vslice = slice(120,4180)
    hslices = odict( [
        ('os', slice(4112+10, 4262-10)), 
        ('ac', slice(20,4100)), 
        #('ac_os', slice(3112+10,3262-10))
        ] )
    countCut = 1e1/dx
    parseHist = lambda h,x: { 'y': h[h/dx>countCut]/dx, 'x': (.5*(x[1:]+x[:-1]))[h/dx>countCut] }
    #parseHist = lambda h,x: { 'y': h, 'x': (.5*(x[1:]+x[:-1])) }
    
    hist = odict( [ (datakey, { 
        #'all': {'hist': np.histogram( clean(datum[0:xccd,0:yccd]), bins = bins )[0] }, 'slice': '0:%s|0:%s'%(xccd,yccd),
        slicekey: {
            'data': clean(datum[vslice,hslice]),
            'hist': parseHist( *np.histogram( clean(datum[vslice,hslice]), bins = bins ) ),
            'pixelcount': len(clean(datum[vslice,hslice])),
            'totalpixel': len(datum[vslice,hslice].flatten()),
            'mean': np.mean(clean(datum[vslice,hslice])),
            'median': np.median(clean(datum[vslice,hslice])),
            'std': np.std(clean(datum[vslice,hslice])),
            'N': len(clean(datum[vslice,hslice])),
            } for slicekey, hslice in hslices.items()
        })
        for datakey, datum in data.items() 
    ])
    if verbose: print 'keys', [ (datakey, hist[datakey].keys()) for datakey in data.keys() ]
    
    def change( d, k, v ):
        d0 = dict(d)
        d0[k] = v
        return d0
    
    crop = lambda u,w: w[ np.abs(w)<np.max(u) ]
    first=True
    fitindex = 0
    for datakey, thishist in hist.items():
        for slicekey in hslices.keys():
            datum = thishist[slicekey]
            print 'pixelcount', datum['pixelcount'], 'total', datum['totalpixel'], 'occupancy', 1.-float(datum['pixelcount'])/datum['totalpixel']
            x = datum['hist']['x']
            y = datum['hist']['y']
            datum['mean2'] = np.mean(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            datum['median2'] = np.median(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            datum['std2'] = np.std(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            datum['N2'] = len(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            N = y.sum()*dx
            N_ = N
            #print 'µ', 'σ', 'λ'
            #print datum['pixelcount']
            #print y.sum()/dx
            #print 2*y[x<0].sum()/dx

            datum['fits'] = odict({})

            p0 = odict( [ 
                ('mu_os', 0.), 
                ('sigma_os', 0.), 
                ('mu_ac',  0.), 
                ('sigma_ac', 0.), 
                ('A', N), 
                ('lamb', 0.), 
                ('gain', gain),
                ] )

            if verbose: print 'std2', slicekey, datum['std2']
            fits = datum['fits']

            std2os = thishist['os']['std2']
            fits[r'G(0,σ\')'] = fitpartial( x, y, p0=p0, adjust = ['A'], sigma_os = std2os )
            Afit = fits[r'G(0,σ\')']['params']['A']
            fits[r'G(0,σ)A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A'], sigma_os = std2os, A = Afit )
            Afit2 = fits[r'G(0,σ)A']['params']['A']
            fits[r'G(µ,σ)A'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'A'], sigma_os = std2os, A = Afit )
            #print N, Afit
            #if slicekey == 'os':
                #fits[r'G(0,σ)A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A'], sigma_os = std2os, A = Afit )
                #fits[r'log(G(0,σ))'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os'], log=True, sigma_os = std2os, A = Afit )
                #fits[r'log(G(0,σ))A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A'], log=True, sigma_os = std2os, A = Afit )
                #pass
                
            if slicekey == 'ac':
                lamb0 = thishist['ac']['std2']**2 - thishist['os']['std2']**2
                if verbose: print 'lambda', lamb0
                #fits[r'G(µ,σ\')*P(λ)'] = fitpartial( x, y, p0=p0, adjust = ['mu_ac', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=2./3 )
                #fits[r'G(0,σ\')*P(λ\')A'] = fitpartial( x, y, p0=p0, adjust = ['lamb'], sigma_os=std2os, A=Afit2, lamb=lamb0 )
                fits[r'G(0,σ\')*P(λ\')'] = fitpartial( x, y, p0=p0, adjust = ['A'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                fits[r'G(µ,σ)*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                Afit = fits[r'G(0,σ\')*P(λ\')']['params']['A']
                Afiterr = fits[r'G(0,σ\')*P(λ\')']['perr']['A']
                #if verbose: print 'fiterr', Afit, Afiterr, str_with_err(Afit, Afiterr)

                #fits[r'G(0,σ)*P(λ)'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(0,σ)*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(0,σos)*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                fits[r'G(µ,σ)*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                fits[r'G(µ,σ)*P(λ)'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                fits[r'G(µ,σos)*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                fits[r'G(µ,σos)*P(λ)'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(µ,σ)*P(λ)A<'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'A', 'lamb'], sel=[x<0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(µ,σ)*P(λ)<'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'lamb'], sel=[x<0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(µ,σos)*P(λ)A<'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'A', 'lamb'], sel=[x<0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(µ,σos)*P(λ)<'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'lamb'], sel=[x<0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(0,σ\')*P(λ)'] = fitpartial( x, y, p0=p0, adjust = ['lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(0,σ\')*P(λ)p'] = fitpartial( x, y, p0=p0, adjust = ['lamb'], sel=[x>0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(0,σ\')*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['A', 'lamb'], sel=[x>0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'G(0,σ\')*P(λ)g'] = fitpartial( x, y, p0=p0, adjust = ['lamb', 'gain'], sel=[x>0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )

                #fits[r'log(G(0,σ\')*P(λ))'] = fitpartial( x, y, p0=p0, adjust = ['lamb'], log=True, sel=[x>0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'log(G(0,σ\')*P(λ))g'] = fitpartial( x, y, p0=p0, adjust = ['lamb','gain'], log=True, sel=[x>0], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                #fits[r'log(G(0,σl)*P(λ))'] = fitpartial( x, y, p0=p0, adjust = ['lamb'], log=True, sel=[x>0], sigma_os=thishist['os']['fits'][r'log(G(0,σ))A']['params']['sigma_os'], A=Afit, lamb=lamb0, gain=gain )

                #fits[r'log(G(µ,σl)*P(λ))'] = fitpartial( x, y, p0=p0, adjust = ['mu_ac','lamb'], log=True, sigma_os=thishist['os']['fits'][r'log(G(0,σ))A']['params']['sigma_os'], A=Afit, lamb=lamb0, gain=gain )

            tablekeys = datum['fits'].keys()
            fitkeylen = max(map(len,tablekeys))+3
            for fitkey in tablekeys:
                if not fitkey in datum['fits']: continue
                line = []
                if first: 
                    header += ['fitID']
                    fmt += ['%2d']
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
                params = ['mu_os', 'sigma_os', 'mu_ac', 'sigma_ac', 'A', 'lamb', 'gain' ]
                if first:
                    header += params
                    fmt += ['%s']*7
                line += [ (str_with_err(datum['fits'][fitkey]['params'][k], datum['fits'][fitkey]['perr'][k]) if k in datum['fits'][fitkey]['params'] else 0) for k in params ]
                if first: 
                    header += ['chisqr']
                    fmt += ['%.3e']
                line += [ datum['fits'][fitkey]['chisqr'] ]
                result += [ line ]
                first = False
                fitindex += 1
    
    header2 = ['runID', 'ohdu', 'tempmin', 'tempmax', 'pixelcount']
    line2 = [ int(runID), int(ohdu), float(hdulists['scn'][0].header['TEMPMIN' ]), float(hdulists['scn'][0].header['TEMPMAX']), datum['pixelcount'] ]
    header2 += [ 'sigma_os', 'err_sigma_os', 'lambda_os', 'err_lambda_os', 'mu_os', 'err_mu_os', 'A_os', 'err_A_os', 'chisq_os' ]
    datum = hist['scn']['ac']['fits'][r'G(µ,σos)*P(λ)A']
    line2 += [ datum['params']['sigma_os'], datum['perr']['sigma_os'],
             datum['params']['lamb'], datum['perr']['lamb'],
             datum['params']['mu_os'], datum['perr']['mu_os'],
             datum['params']['A'], datum['perr']['A'],
             datum['chisqr']
             ]
    header2 += [ 'sigma_os', 'err_sigma_os', 'lambda_os', 'err_lambda_os', 'mu_os', 'err_mu_os', 'A_os', 'err_A_os', 'chisq_os' ]
    datum = hist['scn']['ac']['fits'][r'G(µ,σ)*P(λ)A']
    line2 += [ datum['params']['sigma_os'], datum['perr']['sigma_os'],
             datum['params']['lamb'], datum['perr']['lamb'],
             datum['params']['mu_os'], datum['perr']['mu_os'],
             datum['params']['A'], datum['perr']['A'],
             datum['chisqr']
             ]
    
    #sep = ', '
    sep = ' '
    if verbose:
        print sep.join( [ str((i,h)) for i, h in enumerate(header)] )
        print sep.join(header)
        print '\n'.join( [ sep.join( [ f%entry for f,entry in zip(fmt,line) ] ) for line in result ] )
    printcsv = lambda outfile: np.savetxt('%s.csv'%(outfile), result, header=sep.join(header), fmt='%s', delimiter=sep)
    printcsv2 = lambda outfile: np.savetxt('%s.csv2'%(outfile), [line2], header=sep.join(header2), fmt='%s', delimiter=sep)
    return hist, printcsv, printcsv2

def makePlot( ohdu, runID, hist, verbose = False ):
    nfiles = len(hist.keys())
    slices = set( [ k for v in hist.values() for k in v.keys() ] )
    print 'slices', slices
    nslices = len(slices)
    fig = plt.figure()
    fig.set_size_inches( np.array(fig.get_size_inches())*[nslices, nfiles] )
    fig.suptitle('runID%s'%runID+' ohdu%s'%ohdu)
    first = True
    if 1:
        span = 5 if fft else 4
        
        gs = GridSpec( span*nfiles, nslices )
        ax = { key: { hkey: fig.add_subplot(gs[span*(i-1):span*(i-1)+3,j]) for j, hkey in enumerate(slices) } for i, key in enumerate(hist.keys()) }
        axdiff = { key: { hkey: fig.add_subplot(gs[span*(i-1)+3,j], sharex=ax[key][hkey] ) for j, hkey in enumerate(slices)} for i, key in enumerate(hist.keys()) }
        if fft:
            axfft = { key: { hkey: fig.add_subplot(gs[span*(i-1)+4,j] ) for j, hkey in enumerate(slices)} for i, key in enumerate(hist.keys()) }
    for ax_ in ax.values():
        for ax__ in ax_.values():
            plt.setp(ax__.get_xticklabels(), visible=False)
    
    if verbose: print 'generating plots'
    for key, data in hist.items():
        for slicekey, datum in data.items():
            x = datum['hist']['x']
            
            ax[key][slicekey].step( datum['hist']['x'], datum['hist']['y'], where='mid', label = slicekey, lw = 2 )
            
            if slicekey == 'os':
                listfits = [
                r'G(0,σ\')', 
                r'G(µ,σ)A',
                #r'G(0,σ)N',
                #r'G(0,σ0-)A',
                #r'G(0,σ0-)N',
                #r'log(G(0,σ))',
                #r'log(G(0,σ))N',
                #r'log(G(0,σ))A',
                ]
            else:
                listfits = [
                #r'G(0,σ\')', 
                #r'G(0,σ)', 
                #r'G(µ,σ\')*P(λ)',
                #r'G(0,σ)*P(λ)',
                #r'G(0,σ\')*P(λ\')',
                #r'G(0,σ\')*P(λ)',
                #r'G(0,σ\')*P(λ)g',
                #r'log(G(0,σ\')*P(λ\'))',
                #r'log(G(0,σl)*P(λ))',
                r'G(µ,σ)A',
                r'G(µ,σ)*P(λ)A',
                r'G(µ,σos)*P(λ)A',
                ]
            n = 0
            for i, fitfunc_ in enumerate(listfits):
                if not fitfunc_ in datum['fits'].keys(): continue
                label = fitfunc_
                label = label.replace(r'µ', '%.2f'%datum['fits'][fitfunc_]['params']['mu_os'] )
                label = label.replace(r'σ', '%.2f'%datum['fits'][fitfunc_]['params']['sigma_os'] )
                label = label.replace(r'λ', '' if not 'lamb' in datum['fits'][fitfunc_]['params'] else '%.2g'%abs(datum['fits'][fitfunc_]['params']['lamb']) ) 
                label += ' $\chi^2$=%.3g'%datum['fits'][fitfunc_]['chisqr']
                
                ax[key][slicekey].plot( datum['hist']['x'], datum['fits'][fitfunc_]['hist'], 
                                label = label, 
                                ls = '--', color='C%d'%(n+1) )
                axdiff[key][slicekey].step( datum['hist']['x'], (datum['hist']['y']-datum['fits'][fitfunc_]['hist'])**2/datum['hist']['y'], where='mid', label = label, color='C%d'%(n+1) )
                n+=1
                #ax[key].step( datum['hist']['x'], np.abs(datum['fits'][fitfunc]['deconvolve']), where='mid', label = 'dec '+fitfunc )
            
            if fft:
                #def fourier( x, y ):
                    #w = np.array(x)
                    #f = np.zeros_like(x)
                    #for w_, f_ in zip(w,f):
                        #f_ = scipy.integrate.quad( lambda x_: y[np.digitize(x_, x)-1]*np.exp(-1.j*w_*x_), -np.inf, np.inf )[0]
                    #return f
                #def fourier( x, y ):
                    #dx = x[1]-x[0]
                    
                #fourier = lambda w: np.fft.fftshift( np.fft.fft( w ))
                y = fourier( datum['hist']['x'], datum['hist']['y'] )
                print y.shape
                ly = np.log( y )
                #axfft[key][slicekey].step( datum['hist']['x'], np.imag(ly), where='mid', label='i' )
                #lamb = np.mean( np.real(ly) )
                print slicekey
                axfft[key][slicekey].step( datum['hist']['x'], np.abs(np.real(y)), where='mid', label='r' )
                axfft[key][slicekey].step( datum['hist']['x'], np.abs(np.imag(y)), where='mid', label='i' )
                
                #y2 = fourier(datum['fits'][r'G(0,σ)']['hist'])
                #axfft[key][slicekey].plot( datum['hist']['x'], np.abs(np.real(y2)), label='rG' )
                #axfft[key][slicekey].plot( datum['hist']['x'], np.abs(np.imag(y2)), label='iG' )

                #axfft[key][slicekey].step( datum['hist']['x'], np.abs(np.real(y/y2)), where='mid', label='r/G' )
                #axfft[key][slicekey].step( datum['hist']['x'], np.abs(np.imag(y/y2)), where='mid', label='i/G' )

                #y2 = fourier(datum['fits'][r'G(0,σ)*P(λ)']['hist'])
                #axfft[key][slicekey].step( datum['hist']['x'], np.abs(np.real(y2)), where='mid', label='rGP' )
                #axfft[key][slicekey].step( datum['hist']['x'], np.abs(np.real(y2)), where='mid', label='iGP' )
                
                #y3 = fourier( cPoisson(x[x>0],-3,3) )
                #axfft[key][slicekey].step( x[x>0], np.abs(np.real(y3)), where='mid', label='rP' )
                #axfft[key][slicekey].step( x[x>0], np.abs(np.real(y3)), where='mid', label='iP' )
                
                #axfft[key][slicekey].step( datum['hist']['x'], np.abs(ly), where='mid', label='a' )
                #axfft[key][slicekey].step( datum['hist']['x'], fourier(iy), where='mid', label='fiy' )

    if verbose: print 'setting labels and lims'
    for key in hist.keys():
        for hkey in slices:
            #axdiff[key][hkey].set_xlabel('pixel energy [adu]')
            axdiff[key][hkey].set_xlabel('pixel electrons [e-]')
            ax[key][hkey].set_ylabel('pixel count')
            ax[key][hkey].set_yscale('log')
            ax[key][hkey].set_title( key )
            if fft: 
                axfft[key][hkey].set_yscale('log')
                #axfft[key][hkey].set_xlim([maxfreq-100,maxfreq+100])
                #axfft[key][hkey].set_xlim([-5,5])
                axfft[key][hkey].legend(loc=2)
            #ax[key].set_xlim([-120,120])
            ax[key][hkey].set_ylim(bottom=min(hist[key][hkey]['hist']['y']))
            ax[key][hkey].grid(True)
            axdiff[key][hkey].set_yscale('log')
            axdiff[key][hkey].grid(True)
            #axdiff[key][hkey].set_ylim(bottom=5e-1,2e0])
            ax[key][hkey].legend()
            #if fft: 
            #axdiff[key][hkey].legend()
    plt.subplots_adjust(hspace=.1, wspace=.05)
    return lambda outfile: fig.savefig('%s.png'%(outfile), bbox_inches='tight')

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
        
    hist, printcsv, printcsv2 = computeFits( ohdu, runID, astropy.io.fits.open( FITSfiles[0] ), verbose, gain )
    print 'tables generated in CSV file', '%s.csv(2)'%(FITSfiles[0])
    printcsv(FITSfiles[0])
    #print 'table generated in CSV file', '%s.csv2'%(FITSfiles[0])
    printcsv(FITSfiles[0])

    if plot:
        saveplot = makePlot( ohdu, runID, hist )
        print 'plot generated in PNG file', '%s.png'%(FITSfiles[0])
        saveplot(FITSfiles[0])

def plot2( outfolder, ohdu, runID, ROOTfile = None, plot=False, gain=None, verbose=True, crosstalk=False, onlyoccupancy=False ):
    if verbose: print 'runID', runID
    outputfile = '%s/runID_%s_ohdu_%s'%(outfolder, runID, ohdu)
    result = removeHitsFromrunID( runID, outputfolder=outfolder, ohdu = ohdu, ROOTfile=ROOTfile, output=False, image=False, gain=gain, verbose=verbose, crosstalk=crosstalk, onlyoccupancy=onlyoccupancy )
    if result is None: return
    hist, printcsv, printcsv2 = computeFits(ohdu, runID, 
                                    result,
                                    plot=plot,
                                    gain=gain,
                                    verbose=verbose,
                                )
    if not os.path.exists( outfolder ):
        print 'creating directory', outfolder
        os.makedirs( outfolder )
    print 'tables generated in CSV file', '%s.csv(2)'%(outputfile)
    printcsv(outputfile)
    #print 'table generated in CSV file', '%s.csv'%(outputfile)
    printcsv2(outputfile)
    
    if plot:
        saveplot = makePlot( ohdu, runID, hist )
        print 'plot generated in PNG file', '%s.png'%(outputfile)
        saveplot(outputfile)
    return

def help_():
    print 'Usage:'
    print sys.argv[0], '--runID <runID> --outfolder <path> --removehits'
    print '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in the associated ROOT catalog'
    print
    print sys.argv[0], '--runID <runID> --outfolder <path> --ROOTfile <path> --removehits'
    print '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in given ROOTfile'
    print 
    print sys.argv[0], '--ohdu <ohdu> --fft --plot <list of paths>'
    print '\tgenerate plot with FITS files from paths for the given ohdu'
    print
    print sys.argv[0], '--ohdu <ohdu> --runID <runID> --outfolder <outfolder> --gain <gain> --table2'
    print '\tgenerate plots and table of function fits for the respective runID and given ohdu'
    print
    print sys.argv[0], '--ohdu <ohdu> --run <run> --outfolder <outfolder> --gain <gain> --table2'
    print '\tgenerate tables of function fits for all the runIDs in the given run and given ohdu'
    print

def getOption( vars_, option, args, f = lambda x:x, flag=False ):
    if '--%s'%option in args:
        vars_[option] = f(args[args.index('--%s'%option)+1]) if not flag else True
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
    p = Program()
    exit(0)
    vars_ = {}
    outfolder = None
    ROOTfile = None
    fft = False
    if '--outfolder' in sys.argv:
        outfolder = sys.argv[sys.argv.index('--outfolder')+1]
        print 'setting outfolder =', outfolder
    getOption( vars_, 'ohdu', sys.argv, int )
    getOption( vars_, 'gain', sys.argv, float )
    getOption( vars_, 'runID', sys.argv, int )
    getOption( vars_, 'run', sys.argv )
    getOption( vars_, 'inputpath', sys.argv )
    getOption( vars_, 'crosstalk', sys.argv, flag=True )
    getOption( vars_, 'onlyoccupancy', sys.argv, flag=True )
    
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
    if '--table2' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA ( 'remove and analyse', lambda: plot2(outfolder, vars_['ohdu'], vars_['runID'], ROOTfile, plot=True, gain=vars_['gain'], crosstalk=vars_['crosstalk'], verbose=False, onlyoccupancy=vars_['onlyoccupancy'] ) )
                exit(0)
            if not vars_['run'] is None:
                l = sorted(list(set( listrunID_run( vars_['run'] ) )))
                print 'from', min(l), 'to', max(l), ' total', len(l)
                runETA( 'remove and analyse full run '+vars_['run'], 
                       cmd = lambda x: plot2(outfolder, vars_['ohdu'], int(x), ROOTfile, gain=vars_['gain'], verbose=False, crosstalk=vars_['crosstalk'], onlyoccupancy=vars_['onlyoccupancy'] ),
                       loop = l
                       )
                exit(0)
            exit(0)
        print term.red('error: ohdu not set')
    if '--table' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA( 'analyse', lambda: plotrunID( vars_['ohdu'], fft=fft, runID=vars_['runID'], plot=False, inputpath=vars_['inputpath'], gain=vars_['gain']) )
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
