# coding: utf-8

import astropy
import astropy.io
import astropy.io.fits

import math
import numpy as np
from numpy.lib.recfunctions import append_fields
import scipy
import scipy.stats
import scipy.signal
import scipy.special
from scipy.misc import factorial
import scipy.optimize

from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('qt4agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import patches

from functools import partial

import glob
import os
import sys
import re
import copy
import root_numpy
import datetime
import time
import random

import json
#np.seterr(all='raise')

class ROOTcatalog:
    pass

class SCNimage:
    def __init__(self, path):
        self.image = astropy.io.fits.open(path)
    
def rglob( pattern, n=0 ):
    if n>20: return []
    return rglob( pattern.replace('//','//*/'), n+1 ) if len(glob.glob( pattern )) == 0 else glob.glob( pattern )

path_connie = '/share/storage2/connie/'
path_processed02data = 'data_analysis/processed02_data/runs/'
path_nuprocessing = 'nu_processing/scripts/ProcCat/'

class Run:
    def initiate(self):
        self.pattern = '*_runID_%s_*'%self.run
        self.path = rglob(path_connie+path_processed02data+self.run+'/data_*/')[0]
        self.range = re.search( r'/data_([0-9]*_to_[0-9]*)/', self.path ).groups()
        print self.path, self.range
        
    def __init__(self, run=None):
        if run:
            if type(run) is int:
                self.run = '%03d*'%run
            else: 
                self.run = run
            self.initiate()
    
    def listAll(self):
        pass
    
    def listrunIDs(self):
        pass

pickfirst = lambda x: None if len(x)==0 else x[0]
class RunID:
    path_scn = path_connie+path_processed02data+'<subrun>/data_<range>/scn/merged/scn_mbs_osi_runID_<run>_<runID>_Int-<integration>_Exp-<exposure>_<start>_to_<stop>.fits'
    path_osi = path_connie+path_processed02data+'<subrun>/data_<range>/osi/images/osi_runID_<run>_<runID>_Int-<integration>_Exp-<exposure>_<start>_to_<stop>_p<part>.fits'
    path_catalog = path_connie+path_processed02data+'<subrun>/data_<range>/ext/catalog/catalog_data_<range>.root'
    path_gaincatalog = path_connie+path_nuprocessing+'hpixP_cut_scn_osi_raw_gain_catalog_data_<range>.root'
    
    def initiate(self):
        self.pattern = '*_runID_*_%05d_*'%self.runID
        self.path_scnmerged = pickfirst( rglob(path_connie+path_processed02data+'*/data_*/scn/merged/'+self.pattern) )
        print 'scn merged:', self.path_scnmerged
        self.subrun, self.range, self.run = re.search( r'runs/([0-9]+.)/data_(.*?)/.*_runID_([0-9]+)_', self.path_scnmerged ).groups()
        print self.subrun, self.range, self.run
        self.path_osiparts = rglob(path_connie+path_processed02data+'%s/data_%s/osi/images/'%(self.subrun,self.range)+self.pattern)
        print 'osi parts:'
        print '\n'.join(self.path_osiparts)
        self.path_catalog = pickfirst( rglob(path_connie+path_processed02data+'%s/data_%s/ext/catalog/catalog_data_*.root'%( self.subrun, self.range ) ) )
        print 'catalog:', self.path_catalog
        self.path_gaincatalog = pickfirst( rglob(path_connie+path_nuprocessing+'*scn_osi_raw_gain_catalog_%s.root'%self.range ) )
        print 'gain catalog:', self.path_gaincatalog
        print '\n'.join( sorted( map( os.path.basename, rglob(path_connie+path_nuprocessing+'*data_*[!skim1].root') ) ) )
        
    def __init__(self, runID=None):
        if runID:
            self.runID = int(runID)
            self.initiate()
    
    @staticmethod
    def _getscn(runID):
        return rglob( re.sub( r'\<.*?\>', r'*', RunID.path_scn.replace('<runID>', '%05d'%runID) ) )[0]
    
    @staticmethod
    def _getosi(runID):
        return rglob( re.sub( r'\<.*?\>', r'*', RunID.path_osi.replace('<runID>', '%05d'%runID) ) )

    @staticmethod
    def _getcatalog(runID):
        subrun, _range, run = RunID._getrun( RunID._getscn( runID ) )
        return rglob( re.sub( r'\<.*?\>', r'*', RunID.path_catalog.replace('<range>', _range) ) )
    
    @staticmethod
    def _getgaincatalog(runID):
        subrun, _range, run = RunID._getrun( RunID._getscn( runID ) )
        return rglob( re.sub( r'\<.*?\>', r'*', RunID.path_gaincatalog.replace('<range>', _range ) ) )
    
    @staticmethod
    def _getrun(s):
        return re.search( r'runs/([0-9]+.)/data_(.*?)/.*_runID_([0-9]+)_', s ).groups()
    
    @staticmethod
    def _getrunID(s):
        return re.search( r'runs/.*/.*/.*_runID_.*_([0-9]+)_', s ).groups()
    
    @staticmethod
    def listall():
        for entry in sorted(rglob( path_connie+path_processed02data+'*/data_*/scn/merged/*' )):
            runID = RunID._getrunID(entry)[0]
            subrun, _range, run = RunID._getrun(entry)
            print runID, subrun, '\t',
        print

#    def listFITS( *patterns ):
#        '''
#        list all files that match the list of patterns given
#        '''
#        return sortByrunID([ match for pattern in patterns for match in rglob(pattern) if not '-with-' in match ])
    
    def getSCN():
        self.image_scn = astropy.io.fits.open( self.path_scnmerged )
        return self.runID

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
                self._cmdline += [ arg[2:] ]
                self.setparam(arg[2:], (argv[i+1] if not self.is_option(argv[i+1]) else True) if i < len(argv)-1 else None )
        return
    
    def get_cmdline(self):
        for _0 in self._cmdline:
            print _0, self.__dict__[_0]
    
    def getlast(self):
        return self._cmdline[-1]

class Program:
    def __init__(self):
        self.options = Options( sys.argv )
        self.options.get_cmdline()
        self._calllast()
        
    def _callbystr(self, s):
        if s in Program.__dict__:
            return Program.__dict__[s](self)
        else:
            return self._help()
    
    def _calllast(self):
        return self._callbystr( self.options.getlast() )
    
    def _help(self):
        print 'usage:'
        for memberfunction in Program.__dict__:
            if not '_' == memberfunction[0]:
                print '--%s:\t'%memberfunction, Program.__dict__[ memberfunction ].__doc__
                print
        
    def test(self):
        """test function"""
        run = Run(35)
        runID = RunID( self.options.runID )
        print RunID._getscn( int(self.options.runID) )
        print RunID._getosi( int(self.options.runID) )
        print RunID._getcatalog( int(self.options.runID) )
        print RunID._getgaincatalog( int(self.options.runID) )
        #RunID.listall()
    
#print 'test'
#p = Program()
#exit(0)
    
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
    #colorscale = staticmethod( lambda x: '\033[48;2;0;0;%dm \033[0m'%int(x) )
    #colorscale = staticmethod( lambda x, min_, max_: '\033[48;5;%dm \033[0m'%int( ( 51 - 16 )*( x - min_)/( max_ - min_ ) + 16 ) )
    colorscaleYR = staticmethod( lambda x, min_, max_: '\033[48;5;%dm \033[0m' % (1+6*int( (( 196/6 - 226/6 )*( x - min_)/( max_ - min_ ) + 226/6 ))) )
    colorscaleGray = staticmethod( lambda x, min_, max_: ('\033[48;5;%dm \033[0m' % (int( (( 255 - 232 )*( x - min_)/( max_ - min_ ) + 232 ))) )*2 )
    
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
    
    readcatalog_once = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu','level','E0','xMin','xMax','yMin','yMax','xBary0','yBary0','flag'] + (['gainCu'] if readGain else []), selection='runID==%s'%runID, start=start_, stop=stop_+1 )
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
    angularplot = True

    if angularplot:
        hitscatalog = root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['flag','E3','E0','xMax','xMin','yMax','yMin','xBary3','yBary3','gainCu','n3','n0','runID','ohdu','xVar3','yVar3','xVar0','yVar0'] )
        #hitscatalog = hitscatalog[ hitscatalog['ohdu'] == 5 ]
        
        get_prop = lambda key: hitscatalog[key][ np.all( [ 
            hitscatalog['flag'] == 0, 
            hitscatalog['E3'] > 0,
            hitscatalog['E0'] > 0,
            hitscatalog['n3'] > 0,
            hitscatalog['n0'] > 0,            
            hitscatalog['xVar3'] > 0,
            hitscatalog['xVar0'] > 0,
            hitscatalog['yVar3'] > 0,
            hitscatalog['yVar0'] > 0,
            ], axis=0) ]
        
        def make_histogram( figure, position, data, xlabel, N_bins, *selections, **kwargs ):
            ax = figure.add_subplot( position )
            ax.set_xlabel(xlabel)
            make_bins = lambda Q, N: np.linspace( Q.min(), Q.max(), N ) if 'xscale' in kwargs else np.geomspace( Q.min(), Q.max(), N ) 
            get_prop_or = lambda prop, alt: kwargs[prop] if prop in kwargs else alt
            ax.set_yscale( get_prop_or('yscale','log') )
            ax.set_xscale( get_prop_or('xscale','log') )
            select_all = data==data
            
            overall_selection = get_prop_or( 'overallselection', select_all )
            #print 'len', len(selections[0])
            selections = [select_all] + list(selections[0])
            #print 'len', len(selections)
            for selection in selections:
                #print selection.shape
                dist, bins = ax.hist( data[selection], 
                                     bins = make_bins( data, N_bins ), 
                                     histtype = get_prop_or('histtype','step'), 
                                     label = r'%.3f %.3f'%(np.mean(data[selection]), np.std(data[selection])), 
                                     normed = get_prop_or('normed',False) 
                                     )[0:2]
                
                if 'fitfunc' in kwargs:
                    x = .5*( bins[:-1] + bins[1:] )
                    x_max = x[dist==dist.max()]
                    pp = scipy.optimize.curve_fit( kwargs['fitfunc'], x[ x > x_max ], dist[ x > x_max ] )[0]
                    ax.plot( x, kwargs['fitfunc']( x, pp ), label='fit' )
            ax.legend(fancybox=True, framealpha=0)
            ax.grid(True)
            return

        
        energy_adu = get_prop('E3')
        energy0_adu = get_prop('E0')
        adu_keV = get_prop('gainCu')
        print 'energy_adu.shape', energy_adu.shape
        print 'adu_keV.shape', adu_keV.shape

        fig = plt.figure(figsize=(20,10))
        #electron_bin_size = 2
        energy_keV = energy_adu/adu_keV

        gcenter_x_pixel = .5*( get_prop('xMax') + get_prop('xMin') )
        gcenter_y_pixel = .5*( get_prop('yMax') + get_prop('yMin'))

        dist_bary_gcenter_pixel = np.sqrt( ( get_prop('xBary3') - gcenter_x_pixel )**2 + ( get_prop('yBary3') - gcenter_y_pixel)**2 )
        dist_bary_gcenter_pixel_selection =  dist_bary_gcenter_pixel<10
        #dist_bary_gcenter_plot = fig.add_subplot(222)
        #dist_bary_gcenter_plot.set_xlabel(r'dist(bary,geo)[px]')
        #dist_bary_gcenter_plot.set_yscale('log')
        #dist_bary_gcenter_bins = dist_bary_gcenter_plot.hist( dist_bary_gcenter_pixel, bins=200, histtype='step' )[1]

        pixel_size_um = 15
        length_x_um = ( get_prop('xMax') - get_prop('xMin') + 1 )*pixel_size_um
        length_y_um = ( get_prop('yMax') - get_prop('yMin') + 1 )*pixel_size_um
        length_xy_um = np.sqrt( length_x_um**2 + length_y_um**2 )
        length_xy_pixel = length_xy_um/pixel_size_um
        print 'min length_xy_pixel', min( length_xy_pixel[ length_xy_pixel > 0 ])
        length_z_um = 675.
        length_um = np.sqrt(length_xy_um**2 + length_z_um**2)
        length_pixel = length_um/pixel_size_um

        dist_bary_gcenter_pixel_shift = dist_bary_gcenter_pixel + 1e-5
        dist_bary_gcenter_length = dist_bary_gcenter_pixel_shift#/length_pixel
        #dist_bary_gcenter_length_log_bins = log_bins( dist_bary_gcenter_length, 200 )
        #dist_bary_gcenter_length_plot = fig.add_subplot(222)
        #dist_bary_gcenter_length_plot.set_xlabel(r'dist(bary,geo)/L[px]')
        #dist_bary_gcenter_length_plot.set_yscale('log')
        #dist_bary_gcenter_length_plot.set_xscale('log')
        #dist_bary_gcenter_length_plot.hist( dist_bary_gcenter_length, bins=dist_bary_gcenter_length_log_bins, histtype='step' )
        
        number_of_pixels = get_prop('n3')
        number0_of_pixels = get_prop('n0')
        var_x_pixel = get_prop('xVar3')
        var_y_pixel = get_prop('yVar3')
        var0_x_pixel = get_prop('xVar0')
        var0_y_pixel = get_prop('yVar0')
        energy_deposition_keV_um = energy_keV/length_um

        rad_degrees = np.pi/180.
        #azimuthal_angle_selection = np.all( [ 
            #energy_deposition_keV_um > .23, 
            #energy_deposition_keV_um < .35, 
            #energy_keV < energy_max_keV,
            #dist_bary_gcenter_pixel > 1., 
            #dist_bary_gcenter_pixel < 2.5, 
            #length_xy_pixel > 300 
            #], axis=0 )
        azimuthal_angle_selection = np.all( [ 
            #energy_keV < energy_max_keV, 
            length_xy_pixel > 300, 
            #dist_bary_gcenter_pixel > .3, 
            #dist_bary_gcenter_pixel < 3, 
            #dist_bary_gcenter_length < 5, 
            #energy_deposition_keV_um < .4,
            #energy_deposition_keV_um > .25,
            #energy_keV/number_of_pixels < 3
            ], axis=0 )
        
        azimuthal_angle_rad = np.arctan2( length_xy_um, length_z_um )
        azimuthal_angle_degrees = azimuthal_angle_rad/rad_degrees
        #print 'min angle', azimuthal_angle_degrees.min()
        min_azimuthal_angle_selection = azimuthal_angle_degrees > azimuthal_angle_degrees.min()
        #print 'min angle', azimuthal_angle_degrees[min_azimuthal_angle_selection].min()
        
        #p = length_pixel > 1
        selection_per_sigma = lambda sigma_times:             ( np.log10( energy_keV/length_um**1.05 ) + 0.38 )**2/(sigma_times*0.12)**2 \
                +( np.log10(energy_keV/length_um**1.05) + 0.45)**2/(sigma_times*0.11)**2 \
                    +( np.log10(energy_keV**1.025/number_of_pixels) - .2)**2/(sigma_times*0.06)**2 \
                            < 1
                        

        Nbins = 500
        basis_sel = length_xy_pixel > 300
        n_std = 4
        
        m = np.mean( np.log10(energy_adu/energy0_adu )[basis_sel] )
        s = n_std*np.std( np.log10(energy_adu/energy0_adu )[basis_sel] )
        
        m1 = np.mean( np.log10(var0_x_pixel/var_x_pixel )[basis_sel] )
        s1 = n_std*np.std( np.log10(var0_x_pixel/var_x_pixel )[basis_sel] )

        m2 = np.mean( np.log10(var0_y_pixel/var_y_pixel )[basis_sel] )
        s2 = n_std*np.std( np.log10(var0_y_pixel/var_y_pixel )[basis_sel] )

        m3 = np.mean( np.log10( energy_keV/length_um )[basis_sel] )
        s3 = n_std*np.std( np.log10( energy_keV/length_um )[basis_sel] )

        m4 = np.mean( np.log10( length_pixel/number_of_pixels )[basis_sel] )
        s4 = n_std*np.std( np.log10( length_pixel/number_of_pixels )[basis_sel] )

        m5 = np.mean( np.log10( energy_keV/number_of_pixels )[basis_sel] )
        s5 = n_std*np.std( np.log10( energy_keV/number_of_pixels )[basis_sel] )

        m6 = np.mean( np.log10( number_of_pixels/number0_of_pixels )[basis_sel] )
        s6 = n_std*np.std( np.log10( number_of_pixels/number0_of_pixels )[basis_sel] )
        
        logvL = np.log10( np.sqrt( var_x_pixel**2 + var_y_pixel**2 ) /length_pixel )
        m7 = np.mean( logvL[basis_sel] )
        s7 = n_std*np.std( logvL[basis_sel] )
        
        
        selections = [ 
            #number_of_pixels > 100,
            #number_of_pixels > 300,
            #number_of_pixels > 500,
            number0_of_pixels > 1000,
            length_xy_pixel > 500, 
            #length_xy_pixel > 400, 
            length_xy_pixel > 300, 
            #length_xy_pixel > 200, 
            length_xy_pixel > 100, 
            #length_xy_pixel > 50,
            #length_xy_pixel < 9,
            #energy_keV < np.max(energy_keV[length_xy_pixel < 9]),
            #dist_bary_gcenter_pixel > 10,
            #selection_per_sigma(1),
            #selection_per_sigma(4)
            #( np.log10(energy_adu/energy0_adu) - m )**2/s**2 < 1,
            np.all( [
                ( np.log10(energy_adu/energy0_adu) - m )**2/s**2 < 1, 
                (np.log10(var0_x_pixel/var_x_pixel ) - m1)**2/s1**2 < 1,
                (np.log10(var0_y_pixel/var_y_pixel ) - m2)**2/s2**2 < 1,
                (np.log10(energy_keV/length_um) - m3)**2/s3**2 < 1,
                (np.log10(length_pixel/number_of_pixels) - m4)**2/s4**2 < 1,
                (np.log10(energy_keV/number_of_pixels) - m5)**2/s5**2 < 1,
                (np.log10(number_of_pixels/number0_of_pixels) - m6)**2/s6**2 < 1,
                (logvL - m7)**2/s7**2 < 1
                ], axis = 0 )
            ]
        
        print 'selection:'
        length_str = 'sqrt( ( xMax - xMin + 1 )^2*15 + ( yMax - yMin + 1 )^2*15 + 675^2 )'
        var_str = 'sqrt( xVar3^2 + yVar3^2 )'
        print 'flag==0 && E3 > 0 && E0 > 0 && n3 > 0 && n0 > 0 && xVar3 > 0 && xVar0 > 0 && yVar3 > 0 && yVar0 > 0 '\
            +'&& ( log10(E3/E0) - %.4f )^2/%.4f^2 < %.4f ' %(m, s, n_std)\
            +'&& ( log10(xVar0/xVar3) %+.4f )^2/%.4f^2 < %.4f ' %(-m1, s1, n_std)\
            +'&& ( log10(yVar0/yVar3) %+.4f )^2/%.4f^2 < %.4f ' %(-m2, s2, n_std)\
            +'&& ( log10(E3/%s) %+.4f )^2/%.4f^2 < %.4f ' %(length_str, -m3, s3, n_std)\
            +'&& ( log10(%s/n3) %+.4f )^2/%.4f^2 < %.4f ' %(length_str, -m4, s4, n_std)\
            +'&& ( log10(E3/n3) %+.4f )^2/%.4f^2 < %.4f ' %( -m5, s5, n_std)\
            +'&& ( log10(n3/n0) %+.4f )^2/%.4f^2 < %.4f ' %( -m6, s6, n_std)\
            +'&& ( log10(%s/%s) %+.4f )^2/%.4f^2 < %.4f ' %(var_str, length_str, -m7, s7, n_std)

        
        
        make_histogram(fig, 331, energy_keV, 'E[keV]', Nbins, selections, xscale='linear' )
        make_histogram(fig, 332, np.log10( energy_keV ), 'E[keV]', Nbins, selections, xscale='linear' )
        make_histogram(fig, 333, np.log10(energy_adu/energy0_adu), r'E3/E0', Nbins, selections, xscale='linear' )

        make_histogram(fig, 334, number_of_pixels/number0_of_pixels, r'n3/n0', Nbins, selections, xscale='linear' )
        make_histogram(fig, 335, np.log10(number_of_pixels/number0_of_pixels), r'n3/n0', Nbins, selections, xscale='linear' )

        #make_histogram(fig, 337, var0_x_pixel/var_x_pixel, r'var0_x/var_x', Nbins, selections, xscale='linear' )
        #make_histogram(fig, 338, np.log( var0_x_pixel/var_x_pixel), r'var0_x/var_x', Nbins, selections, xscale='linear' )

        #make_histogram(fig, 336, var0_y_pixel/var_y_pixel, r'var0_y/var_y', Nbins, selections, xscale='linear' )
        #make_histogram(fig, 339, np.log( var0_y_pixel/ var_y_pixel ), r'var0_y/var_y', Nbins, selections, xscale='linear' )
        
        make_histogram(fig, 336, np.log10(length_pixel/number_of_pixels), r'L/n[px]', Nbins, selections, xscale='linear' )
        make_histogram(fig, 337, np.log10(energy_keV/length_um), r'E/L[keV/$\mu$m]', Nbins, selections, xscale='linear' )
        make_histogram(fig, 338, np.log10(energy_keV/number_of_pixels), r'E/n[keV/px]', Nbins, selections, xscale='linear' )
        
        make_histogram(fig, 339, np.log10( np.sqrt( var_x_pixel**2 + var_y_pixel**2 ) /length_pixel ), r'var/L', Nbins, selections, xscale='linear' )
        #lt5 = azimuthal_angle_degrees>5
        ##make_histogram(fig, 334, azimuthal_angle_degrees[lt5], r'angle[degree]', 30, map( lambda x:x[lt5], selections), xscale='linear', yscale='linear', fitfunc = lambda x,A: 3*A*np.sin(x*rad_degrees)*np.cos(x*rad_degrees)**2 )
        #positive = dist_bary_gcenter_pixel>0
        ##make_histogram(fig, 337, dist_bary_gcenter_pixel[positive]/length_pixel[positive]**1.01, r'dist(bary,geo)/L^{1.01}[px]', Nbins, map( lambda x:x[positive], selections) )
        #make_histogram(fig, 337, np.log10(number_of_pixels/number0_of_pixels), r'n3/n0', Nbins, selections, xscale='linear' )
        #make_histogram(fig, 338, np.abs(np.log10(var_x_pixel/var_y_pixel)), r'var_x/y', Nbins, selections, xscale='linear' )
        #make_histogram(fig, 334, var_x_pixel/(number_of_pixels/length_pixel), r'var_y/(n/L)', Nbins, selections )


        #selections = [
            #selection_per_sigma(1),
            #selection_per_sigma(2),
            #selection_per_sigma(3),
            #selection_per_sigma(4),
            #selection_per_sigma(5)
                #]
        #lt5 = azimuthal_angle_degrees>5
        #make_histogram(fig, 336, azimuthal_angle_degrees[lt5], r'angle[degree]', 30, map( lambda x:x[lt5], selections), xscale='linear', yscale='linear', fitfunc = lambda x,A: 3*A*np.sin(x*rad_degrees)*np.cos(x*rad_degrees)**2, normed=True )


        fig.tight_layout()
        filename = lambda n: 'distribution%s_runID%s.pdf'%(n,runID)
        print 'printing to', filename(0)
        fig.savefig(filename(0))
        dist_bary_gcenter_pixel
        fig2 = plt.figure()
        
        lt5 = azimuthal_angle_degrees>5
        make_histogram(fig2, 111, azimuthal_angle_degrees[lt5], r'angle[degree]', 90, map( lambda x:x[lt5], selections[-1:] ), xscale='linear', yscale='linear', fitfunc = lambda x,A: 3*A*np.sin(x*rad_degrees)*np.cos(x*rad_degrees)**2 )
        
        print 'printing to', filename(1)
        fig2.savefig(filename(1))
        exit(0)

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
            N = 0
            totalN = 0
            if gain[thisohdu] == 0:
                continue
            totalpix = len( hdulist['hits'][-1].data.flatten() )
            Ntiny = np.all( [hitsohdu['E0']/gain[thisohdu] < 0.075], axis = 0).sum()
            Nlow = np.all( [hitsohdu['E0']/gain[thisohdu] > 3, hitsohdu['E0']/gain[thisohdu] < 7], axis = 0).sum()
            Nhigh = np.all( [hitsohdu['E0']/gain[thisohdu] > 250, hitsohdu['E0']/gain[thisohdu] < 350], axis = 0).sum()
            print Nlow, Nhigh
            
                    
            for level in range(4):
                rate = 0
                totalN += (l==level).sum()
                hitsmask[thisohdu][level] = np.zeros_like( hdulist['hits'][-1].data, dtype=bool )
                hitsmask[thisohdu][level][y[l==level],x[l==level]] = True
                N += hitsmask[thisohdu][level].sum()
                rate = float(N)/totalpix
                #print totalN, N
                occupancy[thisohdu][level] = {'S': totalpix, 'N': N, 'rate': rate, 'gain': gain[thisohdu], 'totalN': totalN, 'N3-7': Nlow, 'N250-350': Nhigh, 'N.075': Ntiny }
        else:
            print 'ohdu', thisohdu, 'empty'
    sep=' '
    print '%s/runID_%s_occupancy.csv'%(outputfolder,runID)
    np.savetxt('%s/runID_%s_occupancy.csv'%(outputfolder,runID), [ [ runID, ohdukey, levelkey, data['S'], data['totalN'], data['N'], data['rate'], data['N3-7'], data['N250-350'], data['N.075'] ] for ohdukey, thisohdu in occupancy.items() for levelkey, data in thisohdu.iteritems() ], header=sep.join(['runID', 'ohdu','level', 'S','totalN','N','rate','N3-7','N250-350','N.075']), fmt='%s', delimiter=sep)
    
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

def plotCSV( files, xcol, ycol, xerr = None, yerr = None, title=None ):
    x = []
    y = {}
    for f in sorted(files):
        data = np.genfromtxt(f, names=True)
        runID = int(re.search( 'runID_([0-9]*)',f ).groups()[0] )
        if not ycol in data.dtype.names: continue
        ohdulist = list( set( data['ohdu'] ) )
        x += [ runID ]
        for ohdu in ohdulist:
            if ohdu in y.keys():
                y[ohdu] += [ float( data[data['ohdu'] == ohdu][ycol][-1] ) ]
            else: 
                y[ohdu] = [ float( data[data['ohdu'] == ohdu][ycol][-1] ) ]
    x = np.array(x)
    ymin, ymax = None, None
    for key, y_ in y.items():
        y_ = np.array(y_)
        if ymax is None:
            ymax = np.max(y_)
        else:
            ymax = max( ymax, np.max(y_) )
        if ymin is None:
            ymin = np.min(y_)
        else:
            ymin = min( ymin, np.min(y_) )
    
    print len(x), ymin, ymax
    N = len(y.keys())
    n = int(np.sqrt(N))+1
    m = N/n
    fig = plt.figure()
    fig.suptitle( title )
    
    figON = plt.figure()
    figON.suptitle( title )
    
    axON = figON.add_subplot(111)
    #axON = figON.add_subplot(211)
    #axOFF = figON.add_subplot(212)
    
    maskOFF = np.all( [x>3165, x<3381], axis=0)
    maskON = np.logical_not(maskOFF)
    
    yOFF =  np.array(y[6.])[maskOFF]
    count, bins, _0 = axON.hist( yOFF, histtype='step', label='OFF\nmean=%.2f\nstd=%.2f'%(np.mean(yOFF),np.std(yOFF)), bins = 15, normed=True )
    yON =  np.array(y[6.])[maskON]
    axON.hist( np.array(y[6.])[maskON], histtype='step', label='ON\nmean=%.2f\nstd=%.2f'%(np.mean(yON),np.std(yON)), bins = bins, normed=True )
    
    axON.legend()
    
    figON.savefig('test2.png')

    
    if xerr == 'bin':
        xerr_ = .5*float(x[1]-x[0])
    else:
        xerr_ = 0
    for i, (key, y_) in enumerate(y.items()):
        if yerr == 'Poisson':
            yerr_ = np.sqrt(y_)
        else:
            yerr_ = 0
        ax = fig.add_subplot( n, m, i+1 )
        ax.errorbar( x, y_, xerr=xerr_, yerr=yerr_, fmt='.', label='%d'%key)
        ax.set_ylim([ymin,ymax])
        ax.legend()
        print i%m, i/m, n, m
        ax.tick_params(direction='in')
        ax.grid(True)
        if i%m != 0: 
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel('count')
        if i/m != m: 
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            ax.set_xlabel('runID')
    plt.subplots_adjust(hspace=0, wspace=0)
    fig.savefig('test.png')
    return
    

def plot_track( track ):
    x = track['xPix']
    y = track['yPix']
    z = track['ePix']/track['gainCu'][0]
    
    x -= x.min()
    y -= y.min()
    
    matrix = np.zeros([ x.max()+1, y.max()+1 ])
    print 'matrix.shape', matrix.shape
    matrix[x,y] = z
    
    return matrix

def reduce_matrix( matrix, factor ):
    p = int(2.**factor)
    new_matrix = np.zeros( np.ceil( np.array(matrix.shape)/p ).astype(int) )
    for i in range(new_matrix.shape[0]):
        for j in range(new_matrix.shape[1]):
            new_matrix[i,j] = np.sum( matrix[p*i:p*i+p, p*j:p*j+p] )
    return new_matrix
        
    
def draw_terminal( matrix, log=False, text=[] ):
    min_ = matrix.min()
    max_ = matrix.max()
    for i, line in enumerate(matrix):
        sys.stdout.write( '|' )
        for row in line:
            if log:
                sys.stdout.write( term.colorscaleGray( np.log10(row-min_+1), np.log10(1), np.log10(max_-min_+1) ) )
            else:
                sys.stdout.write( term.colorscaleGray( row, min_, max_ ) )
        sys.stdout.write( '|' )
        if i < len(text): 
            sys.stdout.write( text[i] )
        sys.stdout.write( '\n' )
    return

def array2dict( a ):
    #print 'fields', a.dtype.descr
    return { name: map(float,a[name][0]) if 'O' in type_ else ( str( a[name][0] ) if 'S' in type_ else float( a[name][0] ) ) for name, type_ in a.dtype.descr }

def updatejson( f, entry ):
    data = []
    if os.path.exists(f):
        data = json.load( open(f,'r') )
        if data == '': data = []
    data += [ entry ]
    json.dump( data, open(f,'w'), indent=4, separators=(',', ': ') )

def dict2array( a ):
    keys = [ key for key in a[0].keys() ] # if (not 'Pix' in key) and (not 'level' in key)
    #print 'keys', keys
    m = {float: 'f4', str: 'S256', list: 'O', unicode: 'S256' }
    dtype = [ (str(key), m[ type(a[0][key]) ]) for key in keys] 
    return np.array( [ tuple([ np.array(aa[key]) if type(aa[key]) is list else aa[key] for key in keys ]) for aa in a ], dtype=dtype )

def draw_entry( entry ):
    M = plot_track( entry )
    factor = int(max( M.shape ))/40
    draw_terminal( reduce_matrix(M, factor), log=True, 
                text=[ 
                    'log',
                    'factor=%s'%factor, 
                    'shape=%s,%s'%M.shape, 
                    'flag=%s'%entry['flag'],
                    'E0=%skeV'%(entry['E0']/entry['gainCu'])
                    ] )
    return
    
def interactive_selection( inputcatalog, outputcatalog ):
    number_of_entries = root_numpy.root2array( inputcatalog, treename='hitSumm', branches='runID' ).shape[0]
    print number_of_entries
    get_entry = lambda i: root_numpy.root2array( inputcatalog, treename='hitSumm', start=i, stop=i+1 )[0]
    
    while 1:
        index_ = random.randint(0, number_of_entries )
        entry = get_entry( index_ )

        if not entry['flag'] == 0: continue
        M = plot_track( entry )
        if np.min( M.shape ) <= 7: continue

        entry = append_fields( entry, ('catalog', 'index' ), (inputcatalog, index_), dtypes = ['S256','i4'] )
        factor = int(max( M.shape ))/40
        draw_terminal( reduce_matrix(M, factor), log=True, 
                      text=[ 
                          'log',
                          'factor=%s'%factor, 
                          'shape=%s,%s'%M.shape, 
                          'flag=%s'%entry[0]['flag'],
                          'E0=%skeV'%(entry[0]['E0']/entry[0]['gainCu'])
                          ] )
        char = raw_input('category?([q]uit): ')
        print 'category selected', char
        if char == 'm':
            print 'muon'
            updatejson( 'muon.json', array2dict(entry) )
            continue
        if char == 'n':
            print 'not muon'
            updatejson( 'notmuon.json', array2dict(entry) )
            continue
        if char == 'list':
            entries = json.load( open('muon.json','r') )
            print 'number of muons', len(entries)
            print 'number of notmuons', len(json.load( open('notmuon.json','r') ))
            continue
        if char == 'q':
            print 'quit'
            break
    return

def pca( data_, inv = None, method='invert' ):
    data = data_.copy()
    if not inv is None: 
        datainv = inv.copy()
    #print 'data.shape', data.shape
    for field in data.dtype.names:
        data[field] -= np.mean( data[field] )
        if not inv is None:
            datainv[field] -= np.mean( datainv[field] )
            if method == 'inverteach':
                datainv[field] = 1./datainv[field]
    if not inv is None and method == 'inverteach':
        data = np.append( data, datainv )
    
    term = lambda d1, d2, f1, f2: np.mean( d1[f1]*d2[f2] )/np.sqrt( np.std(d1[f1])*np.std(d2[f2]) )
    invterm = lambda d1, d2, f1, f2: np.mean( 1./(d1[f1]*d2[f2]) )/np.sqrt( np.std(1./d1[f1])*np.std(1./d2[f2]) )
    covmatrix = np.zeros( (len(data.dtype.names),len(data.dtype.names)) )
    for i1, field1 in enumerate(data.dtype.names):
        for i2, field2 in enumerate(data.dtype.names):
            #covmatrix[i1,i2] = np.mean( data[field1]*data[field2] )#/np.sqrt(np.std(data[field1])*np.std(data[field2]))
            covmatrix[i1,i2] = term(data,data,field1,field2)
            if method == 'invert' and not inv is None:
                covmatrix[i1,i2] += 1. / (np.sum( datainv[field1]*datainv[field2])/np.sqrt(np.std(datainv[field1])*np.std(datainv[field2])))
            if method == 'subtract' and not inv is None:
                covmatrix[i1,i2] += invterm( datainv, datainv, field1, field2 )
                #covmatrix[i1,i2] = -term( data, datainv, field1, field2 )

    eigenvalues, eigenvectors = np.linalg.eigh( covmatrix )
    #return eigenvalues, eigenvectors/eigenvalues
    return eigenvalues, eigenvectors

def rotatedata( data, eigenvalues, eigenvectors ):
    print eigenvalues.shape, eigenvectors.shape
    rotateddata = np.array( zip( *[
        np.sum( [ data[ field ] * eigenvectors[ifield,iev] for ifield,field in enumerate( data.dtype.names ) ], axis = 0 ) 
        for iev in range(len(eigenvalues))#[::-1]
        ] ),
        dtype = [ 
            ( ''.join( [ '%+.2e*%s'%(eigenvectors[ifield,iev], field) for ifield, field in enumerate(data.dtype.names)] ) , None) for iev in range(len(eigenvalues)) #[::-1] 
            ]
        )
    return rotateddata

def eccentricity( x, y, N = 2 ):
    e = np.zeros(len(x), dtype=np.complex)
    for i, (ix, iy) in enumerate(zip(x,y)):
        ix -= np.mean(ix)
        iy -= np.mean(iy)
        r2 = ix**2 + iy**2
        theta = np.arctan2( iy, ix )
        e[i] = np.sum( r2 * np.exp( 1.j * N * theta ) )
    return np.abs(e)#, np.angle(e)
    
def muonPixFraction( x, y, e, g ):
    lVar0 = np.zeros( len(x) )
    tVar0 = np.zeros(len(x))
    lLen = np.zeros(len(x))
    tLen = np.zeros(len(x))
    cChiSqr = np.zeros(len(x))
    eChiSqr = np.zeros(len(x))
    c2ChiSqr = np.zeros(len(x))
    meanE = np.zeros(len(x))
    stdE = np.zeros(len(x))
    difflE = np.zeros(len(x))
    difftE = np.zeros(len(x))
    nL = np.zeros(len(x))
    eL = np.zeros(len(x))
    for i, (ix, iy, ie) in enumerate(zip(x,y,e)):
        ix = ix.astype(float)
        iy = iy.astype(float)
        ie /= g[i]
        #print 'len', len(ix),
        ix -= np.mean(ix)
        iy -= np.mean(iy)
        ie -= np.min(ie) - 1
        meanE[i] = np.mean(ie)
        stdE[i] = np.std(ie)
        ixy = np.vstack((ix,iy))
        vals, vecs = np.linalg.eigh( np.cov( ixy, aweights = ie ) )
        it, il = np.dot( vecs, ixy )
        std = lambda x,w: np.sqrt(np.average( np.abs(x - np.average(x, weights = w)), weights = w)**2)
        lVar0[i] = std(il,ie)
        tVar0[i] = std(it,ie)

        tLen[i] = it.max() - it.min() + 1
        lLen[i] = il.max() - il.min() + 1
        nL[i] = len( it[ np.abs(it) < 1 ] )
        eL[i] = np.sum( ie[ np.abs(it) < 1 ] )
        
        if len(it[il<0]) < 2:
            cChiSqr[i] = 1
            eChiSqr[i] = 1
            c2ChiSqr[i] = 1
            difflE[i] = 1
            difftE[i] = 1
            print 'error'
        else:
            t1Var0 = std(it[il<0],ie[il<0])
            l1 = np.mean(il[il<0])
            t2Var0 = std(it[il>0],ie[il>0])
            l2 = np.mean(il[il>0])
            #print 'tVar', tVar0[i], t1Var0, t2Var0, l1, l2, il.min(), il.max()
            difflE[i] = np.abs( np.sum(ie[il<0]) - np.sum(ie[il>0]) )
            difftE[i] = np.abs( np.sum(ie[it<0]) - np.sum(ie[it>0]) )

            tGauss = scipy.stats.norm.pdf(it, 0, tVar0[i])
            lGauss = scipy.stats.norm.pdf(il, 0, lVar0[i])
            tGauss2 = scipy.stats.norm.pdf(it, 0, np.abs( ( t2Var0 - t1Var0 )*(il - l1)/(l2-l1) + t1Var0 ) )#/np.sqrt( np.abs(( t2Var0 - t1Var0 )*(il - l1)/(l2-l1) + t1Var0 ))
            cilinder = tGauss/lLen[i]
            cone = tGauss2/lLen[i]
            circle = tGauss*lGauss
            if tVar0[i] < 0.1 :
                cChiSqr[i] = 1
                eChiSqr[i] = 1
                c2ChiSqr[i] = 1
                continue
            else:
                #cChiSqr[i] = np.sqrt( np.sum( ( ( np.mean(cilinder[cilinder>0])/np.mean(ie[cilinder>0])*ie[cilinder>0] - cilinder[cilinder>0] )**2/( cilinder[cilinder>0] ) ) ) )
                #eChiSqr[i] = np.sqrt( np.sum( ( np.mean(circle[circle>0])/np.mean(ie)*ie[circle>0] - circle[circle>0] )**2/circle[circle>0] ) )
                #c2ChiSqr[i] = np.sqrt( np.sum( ( np.mean(cone[cone>0])/np.mean(ie[cone>0])*ie[cone>0] - cone[cone>0] )**2/cone[cone>0] ) )
                cChiSqr[i] = np.sqrt( np.sum( ( ( np.mean(cilinder[cilinder>0])/np.mean(ie[cilinder>0])*ie[cilinder>0] - cilinder[cilinder>0] )**2 ) ) )
                eChiSqr[i] = np.sqrt( np.sum( ( np.mean(circle[circle>0])/np.mean(ie)*ie[circle>0] - circle[circle>0] )**2 ) )
                c2ChiSqr[i] = np.sqrt( np.sum( ( np.mean(cone[cone>0])/np.mean(ie[cone>0])*ie[cone>0] - cone[cone>0] )**2 ) )
                if cChiSqr[i] > 1e30:
                    cChiSqr[i] = 1
                if c2ChiSqr[i] > 1e30:
                    c2ChiSqr[i] = 1
                if eChiSqr[i] > 1e30:
                    eChiSqr[i] = 1
                
    return np.array( zip( *[ lVar0, tVar0, lLen, tLen, cChiSqr, eChiSqr, c2ChiSqr, meanE, stdE, difflE, difftE, nL, eL ] ), dtype = [('lVar0', 'f4'), ('tVar0', 'f4'), ('lLen', 'f4'), ('tLen', 'f4'), ('cChiSqr', 'f4'), ('eChiSqr', 'f4'), ('c2ChiSqr', 'f4'), ('meanE', 'f4'), ('stdE', 'f4'), ('difflE', 'f4'), ('difftE', 'f4'), ('nL', 'f4'), ('eL', 'f4')] )

def apply_to_structured( data, func ):
    return np.array( zip( *[ func(data[field]) for field in data.dtype.names ]), dtype = data.dtype )

def append_struct( a, b ):
    return append_fields( a, b.dtype.names, [ b[name] for name in b.dtype.names ] )

def findclusters( a ):
    
    pca_ = pca(a)
    

def machinelearning( inputpositive, inputnegative, catalog, outputcatalog ):
    #positivedata = root_numpy.root2array( inputpositive, treename='hitSumm' )
    #negativedata = root_numpy.root2array( inputnegative, treename='hitSumm' )
    load = json.load( open('muon.json','r') )
    positivedata = dict2array( load )
    negativedata = dict2array( json.load( open('notmuon.json','r') ) )
    print len(positivedata)
    print len(negativedata)
    print 'positivedata.descr', positivedata.dtype.descr
    
    extrapos = muonPixFraction( positivedata['xPix'], positivedata['yPix'], positivedata['ePix'], positivedata['gainCu'] )
    #print extrapos['cChiSqr']
    positivedata = append_struct( positivedata, extrapos )
    extraneg = muonPixFraction( negativedata['xPix'], negativedata['yPix'], negativedata['ePix'], positivedata['gainCu'] )
    negativedata = append_struct( negativedata, extraneg )
    fields_ = [ field for field in positivedata.dtype.names if not field in ['xPix', 'yPix', 'ePix', 'level', 'flag', 'hPixFlag', 'osMadSCN', 'xBary0', 'xBary1', 'xBary2', 'xBary3', 'xMax', 'xMin', 'yBary0', 'yBary1', 'yBary2', 'yBary3', 'yMax', 'yMin', 'index', 'catalog', 'osMadRaw', 'dcVar', 'selFlag', 'runID', 'ohdu', 'expoStart', 'resCu', 'gainCu', 'E1', 'n1', 'E2', 'E3', 'n2', 'n3', 'xVar1', 'xVar2', 'xVar3', 'yVar1', 'yVar2', 'yVar3', 'nSat' ] ]
    print 'fields', fields_
    #fields_ = ['E0', 'lLen', 'cChiSqr', 'meanE' ]
    #predim = [ ('log10(%s)'%(key), lambda x: np.log10(x[key]) ) for key in fields_ ]
    #predim += [ ('log10(%s)/log10(%s)'%(key,key2), lambda x: np.log10(x[key])/np.log10(x[key2]) ) for key in fields_ for key2 in fields_ ]
    #print predim, len(predim)
    dimensions1 = dict( ('log10(%s)'%(key), lambda x,key=key: np.log10(x[key]) ) for key in fields_ )
    #dimensions2 = dict( ('%s'%(key), lambda x,key=key: x[key] ) for key in fields_ )
    #dimensions2 = dict( ('log10(%s*%s)'%(key,key2), lambda x,key=key,key2=key2: np.log10(x[key])+np.log(x[key2]) ) for i, key in enumerate(fields_) for j, key2 in enumerate(fields_) if j>i )
    #a
    #print dimensions, len(dimensions)
    dimensionsManual = {
        #'log10(E0)': lambda x: np.log10(x['E0']/x['gainCu']),
        #'log10(E1/E0)': lambda x: np.log10(x['E1']/x['E0']+1),
        #'log10(n1/n0)': lambda x: np.log10(x['n1']/x['n0']+1),
        ##'E2': lambda x: x['E2'],
        ##'E3': lambda x: x['E3'],
        ##'n0': lambda x: x['n0'],
        ##'log10(xVar0)': lambda x: np.log10(x['xVar0']),
        ##'log10(yVar0)': lambda x: np.log10(x['yVar0']),
        ##'Var0': lambda x: np.sqrt( x['yVar0']**2+x['xVar0']**2 ),
        ##'log10(yVar0/xVar0)': lambda x: np.log10(x['yVar0']/x['xVar0']),
        ##'log10(xL)': lambda x: np.log10(x['xMax']-x['xMin']),
        ##'log10(yL)': lambda x: np.log10(x['yMax']-x['yMin']),
        ##'log10(L)': lambda x: np.log10( np.sqrt( ( x['yMax'] - x['yMin'])**2+( x['yMax'] - x['yMin'] )**2 ) ),
        ##'log10(Var0/L)': lambda x: np.log10( np.sqrt( x['yVar0']**2+x['xVar0']**2 )/( np.sqrt( ( x['yMax'] - x['yMin'])**2+( x['yMax'] - x['yMin'] )**2 ) )),
        ##'log10(n0)': lambda x: np.log10(x['n0']),
        #'log10(E0/n0)': lambda x: np.log10( x['E0']/x['gainCu']/x['n0'] ),
        ##'log10(E0/Var0)': lambda x: np.log10(  x['E0']/x['gainCu']/np.sqrt( x['yVar0']**2+x['xVar0']**2) ),
        #'log10(E0/lLen)': lambda x: np.log10( x['E0']/x['gainCu']/x['lLen'] ),
        ##'e2': lambda x: eccentricity( x['xPix'], x['yPix'] ),
        #'log10(eChiSqr)': lambda x: np.log10(x['eChiSqr']),
        #'log10(cChiSqr)': lambda x: np.log10(x['cChiSqr']),
        #'log10(c2ChiSqr)': lambda x: np.log10(x['c2ChiSqr']),
        #'log10(tVar0)': lambda x: np.log10(x['tVar0']),
        #'log10(cChiSqr/eChiSqr)': lambda x: np.log10(x['cChiSqr']/x['eChiSqr']),
        #'log10(cChiSqr/c2ChiSqr)': lambda x: np.log10(x['cChiSqr']/x['c2ChiSqr']),
        #'log10(c2ChiSqr/eChiSqr)': lambda x: np.log10(x['c2ChiSqr']/x['eChiSqr']),
        #'log10(lVar0)': lambda x: np.log10(x['lVar0']),
        #'log10(tVar0/lVar0)': lambda x: np.log10(x['tVar0']/x['lVar0']),
        #'log10(lVar0/lLen)': lambda x: np.log10(x['lVar0']/x['lLen']),
        'log10(cChiSqr/lLen)': lambda x: np.log10(x['cChiSqr']/x['lLen']),
        'log10(cChiSqr/E0)': lambda x: np.log10(x['cChiSqr']/x['E0']),
        'log10(cChiSqr/stdE)': lambda x: np.log10(x['cChiSqr']/x['stdE']),
        'log10(cChiSqr/meanE)': lambda x: np.log10(x['cChiSqr']/x['meanE']),
        
        'log10(eChiSqr/lLen)': lambda x: np.log10(x['eChiSqr']/x['lLen']),
        'log10(eChiSqr/E0)': lambda x: np.log10(x['eChiSqr']/x['E0']),
        'log10(eChiSqr/stdE)': lambda x: np.log10(x['eChiSqr']/x['stdE']),
        'log10(eChiSqr/meanE)': lambda x: np.log10(x['eChiSqr']/x['meanE']),
        #'log10(cChiSqr/diffE)': lambda x: np.log10(x['cChiSqr']/x['diffE']),
        
        #'log10(c2ChiSqr/lLen)': lambda x: np.log10(x['c2ChiSqr']/x['lLen']),
        #'log10(E0/lVar0)': lambda x: np.log10(  x['E0']/x['gainCu']/x['lVar0'] ),
        #'log10(meanE)': lambda x: np.log10(  x['meanE'] ),
        #'log10(meanE/E0)': lambda x: np.log10(  x['meanE']/x['E0'] ),
        #'log10(stdE)': lambda x: np.log10(  x['stdE'] ),
        #'log10(stdE/E0)': lambda x: np.log10(  x['stdE']/x['E0'] ),
        #'log10(stdE/meanE)': lambda x: np.log10(  x['stdE']/x['meanE'] ),
        #'log10(diffE)': lambda x: np.log10(  x['diffE'] ),
        #'log10(diffE/E0)': lambda x: np.log10(  x['diffE']/x['E0'] ),
        #'log10(diffE/stdE)': lambda x: np.log10(  x['diffE']/x['stdE'] ),
        #'log10(E0/L)': lambda x:np.log10(  x['E0']/x['gainCu']/np.sqrt( ( x['yMax'] - x['yMin'])**2+( x['yMax'] - x['yMin'] )**2 ) ),
        #'log10(n0)': lambda x: np.log10(x['n0']),
        }
    #dimensionsManual['log10(eChiSqr/meanE)'] = lambda x: np.log10(x['eChiSqr']/x['meanE'])
    denomfields_ = ['E0', 'meanE', 'stdE', 'lLen', 'n0', 'tLen', 'lVar0', 'tVar0', 'difftE', 'difflE', 'eL', 'nL' ]
    numerfields_ = denomfields_ + [ 'cChiSqr', 'eChiSqr', 'c2ChiSqr' ]
    dimensions3 = {}
    dimensions3 = dict( ('log10(%s/%s)'%(key,key2), lambda x,key=key,key2=key2: np.log10(x[key]/x[key2]) ) for i, key in enumerate(numerfields_) for j, key2 in enumerate(denomfields_) if i>j )
    dimensions = dimensions3
    #dimensions = {}
    dimensions.update( dimensions1 )
    #dimensions.update( dimensions2 )

    #print [ (key) for key in dimensions.keys() ]
    datapos = np.array( zip( *[ dimension(positivedata) for dimension in dimensions.values() ] ), dtype = [ (key, None) for key in dimensions.keys() ] )
    dataneg = np.array( zip( *[ dimension(negativedata) for dimension in dimensions.values() ] ), dtype = [ (key, None) for key in dimensions.keys() ] )
    #print 'isnan', [ (key, np.argwhere( np.isnan(dataneg[key]) ) ) for key in dimensions.keys() ]
    #print 'isinf', [ (key, np.argwhere( np.isinf(dataneg[key]) ) ) for key in dimensions.keys() ]
    #print datapos[146]

    #print datapos.shape
    #datapos = np.delete( datapos, np.argwhere( datapos['log10(E0/Var0)'] > 3 ) )
    #print datapos.shape
    
    meanpos = apply_to_structured( datapos, lambda x: [np.mean(x)] )
    meanneg = apply_to_structured( dataneg, lambda x: [np.mean(x)] )
    #print np.vstack( (meanpos, meanneg))['E0']
    
    #pcacluster = pca( np.vstack( (meanpos, meanneg)) )
    #print 'eigenvalues cluster', list(enumerate(pcacluster[0]))
    
    pcapos = pca(datapos)
    print 'eigenvalues pos', list(enumerate(pcapos[0]))[-3:]

    pcaneg = pca(dataneg)
    print 'eigenvalues neg', list(enumerate(pcaneg[0]))[-3:]

    pcapos_ineg = pca( datapos, inv = dataneg, method = 'subtract')
    print 'eigenvalues pos_ineg', list(enumerate(pcapos_ineg[0]))[-3:]

    pcaneg_ipos = pca( dataneg, inv = datapos, method = 'subtract' )
    print 'eigenvalues neg_ipos', list(enumerate(pcaneg_ipos[0]))[-3:]

    #pca_ = pcaneg_ipos
    pca_ = pcapos_ineg
    #pca_ = pcacluster
    #pca_ = pcaneg
    #print 'eigenvalues', list(enumerate(pca_[0]))
    print 'eigenvectors0', sorted( [ (dimensions.keys()[i], abs(v)) for i, v in enumerate(pca_[1][:,-1]) ], key= lambda x: x[-1] )[-4:]
    print 'eigenvectors1', sorted( [ (dimensions.keys()[i], abs(v)) for i, v in enumerate(pca_[1][:,-2]) ], key= lambda x: x[-1] )[-4:]
    print 'eigenvectors2', sorted( [ (dimensions.keys()[i], abs(v)) for i, v in enumerate(pca_[1][:,-3]) ], key= lambda x: x[-1] )[-4:]
    rpos = rotatedata( datapos, *pca_ )
    rpos_pos = rotatedata( datapos, *pcapos )
    rpos_neg = rotatedata( dataneg, *pcapos )
    rneg = rotatedata( dataneg, *pca_ )
    rneg_pos = rotatedata( datapos, *pcaneg )
    rneg_neg = rotatedata( dataneg, *pcaneg )
    
    #centers, radii = findclusters( rpos )
    #print 'eigenvalues', pca_[0]
    centers = apply_to_structured( rpos, lambda x: [np.mean(x)] )[0]
    stds = apply_to_structured( rpos, lambda x: [3*np.std(x)] )[0]
    #print 'centers', centers
    #print 'stds', stds
    tops = apply_to_structured( rpos, lambda x: [np.max(x)] )[0]
    bottoms = apply_to_structured( rpos, lambda x: [np.min(x)] )[0]
    
    fig = plt.figure()
    
    kwargs = {'fill': False, 'color': 'k' }

    #axis_pairs = [ [-1,-2], [-2,-3], [-3,-4] ]
    axis_pairs = [ [-1,0], [-2,1], [-3,2] ]
    for i, axis_pair in enumerate(axis_pairs):
        ax2 = fig.add_subplot( len(axis_pairs), 1, i+1 )
        axis0 = axis_pair[0]
        axis1 = axis_pair[1]
        #ax2.scatter( rneg[ rneg.dtype.names[axis0] ], rneg[ rneg.dtype.names[axis1] ] )
        #ax2.scatter( rpos[ rpos.dtype.names[axis0] ], rpos[ rpos.dtype.names[axis1] ] )
        
        ax2.scatter( rneg_neg[ rneg_neg.dtype.names[axis0] ], rpos_neg[ rpos_neg.dtype.names[axis1] ] )
        ax2.scatter( rneg_pos[ rneg_pos.dtype.names[axis0] ], rpos_pos[ rpos_pos.dtype.names[axis1] ] )

        ax2.set_xlabel(rneg_neg.dtype.names[axis0])
        ax2.set_ylabel(rpos_pos.dtype.names[axis1])
        #ax2.add_artist( 
            #patches.Rectangle(( bottoms[axis0], bottoms[axis1] ), width=tops[axis0]-bottoms[axis0], height=tops[axis1]-bottoms[axis1], **kwargs )
            #)
        
    #ax2 = fig.add_subplot(211)
    #axis0 = -1
    #axis1 = -2
    #ax2.scatter( rpos[ rpos.dtype.names[axis0] ], rpos[ rpos.dtype.names[axis1] ] )
    #ax2.scatter( rneg[ rneg.dtype.names[axis0] ], rneg[ rneg.dtype.names[axis1] ] )
    #ax2.set_xlabel(rpos.dtype.names[axis0])
    #ax2.set_ylabel(rpos.dtype.names[axis1])
    #ax2.add_artist( 
        #patches.Rectangle(( bottoms[axis0], bottoms[axis1] ), width=tops[axis0]-bottoms[axis0], height=tops[axis1]-bottoms[axis1], **kwargs )
        #)

    #ax = fig.add_subplot(212)
    #axis0_ = -3
    #axis1_ = -2
    #ax.scatter( rpos[ rpos.dtype.names[axis0_] ], rpos[ rpos.dtype.names[axis1_] ] )
    #ax.scatter( rneg[ rneg.dtype.names[axis0_] ], rneg[ rneg.dtype.names[axis1_] ] )
    #ax.set_xlabel(rpos.dtype.names[axis0_])
    #ax.set_ylabel(rpos.dtype.names[axis1_])
    #ax.add_artist( 
        #patches.Rectangle(( bottoms[axis0_], bottoms[axis1_] ), width=tops[axis0_]-bottoms[axis0_], height=tops[axis1_]-bottoms[axis1_], **kwargs )
        #)

    fig.savefig( 'pca.png' )
    
    alldata, rotatedalldata, selecteddata, rotatedselecteddata = selectCluster( catalog, 0, 10000, dimensions, pca_, tops, bottoms )
    
    #for select in selecteddata:
        #draw_entry( select )
    for i, axis_pair in enumerate(axis_pairs):
        axis0 = axis_pair[0]
        axis1 = axis_pair[1]
        fig.axes[i].scatter( rotatedalldata[ rpos.dtype.names[axis0] ], rotatedalldata[ rpos.dtype.names[axis1] ], s=.1 )
        fig.axes[i].scatter( rotatedselecteddata[ rpos.dtype.names[axis0] ], rotatedselecteddata[ rpos.dtype.names[axis1] ], s=.1 )

    fig.savefig( 'pca2.png' )
    
    muonfig = plt.figure()
    angle = lambda x: np.arctan2( x['lLen']*15, 675 )
    x = angle(alldata)
    x = x[ x >.2 ]
    x = sorted(x)
    x1 = sorted(angle(selecteddata))
    muonax = muonfig.add_subplot(111)
    hist, bins, _ = muonax.hist( x, bins=100, histtype='step' )
    hist1, _, _ = muonax.hist( x1, bins=bins, histtype='step' )
    dx = bins[1]-bins[0]
    fitfunc = lambda x: np.sin(x)*np.cos(x)**2
    #Int = scipy.integrate.quad( fitfunc, bins[0], bins[-1] )[0]
    #muonax.plot( x, np.sum(hist)*fitfunc(x)/Int*dx )
    muonax.plot( x, fitfunc(x)/np.mean(fitfunc(x))*np.mean(hist) )
    muonax.plot( x, fitfunc(x)/np.mean(fitfunc(x))*np.mean(hist1) )
    muonfig.savefig( 'muondistributionNew.png' )
    
    return

def selectCluster( catalog, start, stop, dimensions, pca_, tops, bottoms ):
    data = root_numpy.root2array( catalog, treename='hitSumm', start=start, stop=stop, selection='flag==0' )
    extradata = muonPixFraction( data['xPix'], data['yPix'], data['ePix'], data['gainCu'] )
    data = append_struct( data, extradata )

    #print [ (key) for key in dimensions.keys() ]
    datadim = np.array( zip( *[ dimension(data) for dimension in dimensions.values() ] ), dtype = [ (key, None) for key in dimensions.keys() ] )
    rdata = rotatedata( datadim, *pca_ )
    
    mask = np.all( [ rdata[field] > bottoms[i] for i, field in enumerate(rdata.dtype.names) ] + [ rdata[field] < tops[i] for i, field in enumerate(rdata.dtype.names) ], axis = 0 )
    #rdata = rdata[ mask ]
    #data = data[ mask ]
    return data, rdata, data[mask], rdata[mask]
    

def updateCatalog( outputcatalog, inputcatalog ):
    interactive_selection( inputcatalog, outputcatalog )

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
    if '--plotcsv' in sys.argv:
        #plotCSV( rglob('tmp/*_occupancy.csv'), 'runID', 'N250350', xerr='bin', yerr='Poisson', title='count evolution 250--350keV' )
        plotCSV( rglob('tmp/*_occupancy.csv'), 'runID', 'N37', xerr='bin', yerr='Poisson', title='count evolution 3--7keV' )
        exit(0)
        print term.red('error: ohdu not set')
    if '--listrunID' in sys.argv:
        print 'list of runID'
        print ', '.join( listrunID_run(vars_['run']) )
        exit(0)
    if '--muoncatalog' in sys.argv:
        inputcatalog = sys.argv[ sys.argv.index('-i') + 1 ]
        outputcatalog = sys.argv[ sys.argv.index('-o') + 1 ]
        updateCatalog( outputcatalog, inputcatalog )
        exit(0)
    if '--machinelearning' in sys.argv:
        inputpositive = sys.argv[ sys.argv.index('-p') + 1 ]
        inputnegative = sys.argv[ sys.argv.index('-n') + 1 ]
        catalog = sys.argv[ sys.argv.index('-c') + 1 ]
        outputcatalog = sys.argv[ sys.argv.index('-o') + 1 ]
        machinelearning( inputpositive, inputnegative, catalog, outputcatalog )
        exit(0)
    print 'options not recognized', sys.argv
    help_()
    exit(0)
