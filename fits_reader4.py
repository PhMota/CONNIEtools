# coding: utf-8

import astropy
import astropy.io
import astropy.io.fits

import math
import numpy as np
from numpy.lib.recfunctions import append_fields, stack_arrays
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
from matplotlib import patches, colors

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
import timeit
import traceback

import json

import gi
gi.require_version('Gtk', '2.0')
from gi.repository import Gtk, Gdk, GdkPixbuf

#np.seterr(all='raise')

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

class ROOTcatalog:
    pass

class SCNimage:
    def __init__(self, path):
        self.image = astropy.io.fits.open(path)
    
def rglob( pattern, n=0 ):
    if n>20: return []
    return rglob( pattern.replace('//','//*/'), n+1 ) if len(glob.glob( pattern )) == 0 else glob.glob( pattern )

wstd = lambda x, weights: np.sqrt(np.average( np.abs(x - np.average(x, weights = weights)), weights = weights)**2)

path_connie = '/share/storage2/connie/'
path_processed02data = 'data_analysis/processed02_data/runs/'
path_nuprocessing = 'nu_processing/scripts/ProcCat/'

class Run:
    def initiate(self):
        self.pattern = '*_runID_%s_*'%self.run
        self.path = rglob(path_connie+path_processed02data+self.run+'/data_*/')[0]
        self.range = re.search( r'/data_([0-9]*_to_[0-9]*)/', self.path ).groups()
        print( self.path, self.range )
        
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
        print( 'scn merged:', self.path_scnmerged )
        self.subrun, self.range, self.run = re.search( r'runs/([0-9]+.)/data_(.*?)/.*_runID_([0-9]+)_', self.path_scnmerged ).groups()
        print( self.subrun, self.range, self.run )
        self.path_osiparts = rglob(path_connie+path_processed02data+'%s/data_%s/osi/images/'%(self.subrun,self.range)+self.pattern)
        print( 'osi parts:' )
        print( '\n'.join(self.path_osiparts) )
        self.path_catalog = pickfirst( rglob(path_connie+path_processed02data+'%s/data_%s/ext/catalog/catalog_data_*.root'%( self.subrun, self.range ) ) )
        print( 'catalog:', self.path_catalog )
        self.path_gaincatalog = pickfirst( rglob(path_connie+path_nuprocessing+'*scn_osi_raw_gain_catalog_%s.root'%self.range ) )
        print( 'gain catalog:', self.path_gaincatalog )
        print( '\n'.join( sorted( map( os.path.basename, rglob(path_connie+path_nuprocessing+'*data_*[!skim1].root') ) ) ) )
        
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
            print( runID, subrun, '\t', )
        print()

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
            print( _0, self.__dict__[_0])
    
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
        print( 'usage:' )
        for memberfunction in Program.__dict__:
            if not '_' == memberfunction[0]:
                print( '--%s:\t'%memberfunction, Program.__dict__[ memberfunction ].__doc__ )
                print()
        
    def test(self):
        """test function"""
        run = Run(35)
        runID = RunID( self.options.runID )
        print( RunID._getscn( int(self.options.runID) ) )
        print( RunID._getosi( int(self.options.runID) ) )
        print( RunID._getcatalog( int(self.options.runID) ) )
        print( RunID._getgaincatalog( int(self.options.runID) ) )
        #RunID.listall()
    
#print( 'test' )
#p = Program()
#exit(0)
    
REMOVE = 0
REMOVEOLD = -1e9

eIonization = 3.745 #eV/e-
eV_e = eIonization

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
    #print( 'inside', value, error, digitsE, digitsV )
    if digitsE<digitsV:
        digitsV = digitsE
    return "{0:.{2}f}({1:.0f})e{3:+03d}".format(value*10**digitsV, error*10**(digitsE), digitsE-digitsV, -digitsV )
    
#print( str_with_err( 0.0001, 0.0000005) )
#print( str_with_err( 1, .5) )
#print( str_with_err( 1, 0) )
#print( str_with_err( 10, 5) )
#print( str_with_err( 10000, 500) )
#exit()
#folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/scn/merged/*'
folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029*/data_*/scn/merged/*'
folder2 = '/share/storage2/connie/data_analysis/processed02_data/runs/009*/data_*/scn/merged/*'

def listFITS( *patterns ):
    '''
    list all files that match the list of patterns given
    '''
    return sortByrunID([ match for pattern in patterns for match in rglob(pattern) if not '-with-' in match ])

def listAll_runs():
    l = glob.glob( '/share/storage2/connie/data_analysis/processed02_data/runs/*/' )
    return [ re.search(r'runs/(.+)', i ).groups()[0] for i in l ]

def listFITS_run( *runs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*%s[A-Z]/data_*/scn/merged/scn_*.fits'%run for run in runs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*%s/data_*/scn/merged/scn_*.fits'%run for run in runs ] )
    return l

def listFITSosi_run( *runs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*%s[A-Z]/data_*/osi/images/osi_*.fits'%run for run in runs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*%s/data_*/osi/images/osi_*.fits'%run for run in runs ] )
    return l

def listFITSraw_runID( *runIDs ):
    #print( 'runIDs', runIDs )
    l = listFITS( *[ '/share/storage2/connie/data/runs/*/runID_*_%05d_*.fits.fz'%run for run in runIDs ] )
    return l

def listrunID_run( *runs ):
    l = listFITS_run( runs )
    #print( '\n'.join(l) )
    return map( getrunIDFromPath, l )

def listrunID_run( *runs ):
    l = listFITS_run( runs )
    #print( '\n'.join(l) )
    return map( getrunIDFromPath, l )

def listFITSscn_runID( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/scn_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/sim_*/scn/images/scn_mbs_osi_draw_runID_*_%05d_Int*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        return 'file not found'
    return l

def listFITSosi_runID( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/osi/images/osi_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/sim_*/osi/images/osi_draw_runID_*_%05d_Int*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        return 'file not found'
    return l

def listFITSmbs_runID( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/mbs/images/mbs_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/sim_*/mbs/images/mbs_draw_runID_*_%05d_Int*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        return 'file not found'
    return l

def listremovedFITS_runID( *runIDs, **kwargs ):
    #if os.access( '/share/storage2/connie/data_analysis/processed02_data/runs/', os.W_OK):
        #print( 'input path', '/share/storage2/connie/data_analysis/processed02_data/runs' )
        #return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/*hsub_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    inputpath = '.'
    if 'inputpath' in kwargs:
        inputpath = kwargs['inputpath']
    print( 'input path', inputpath )
    return listFITS( *[ inputpath.rstrip('/') + '//*hsub_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    #return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/*runID_*_%05d_*'%runID for runID in runIDs ] )

def getrunIDFromPath( path ):
    return re.search(r'runID_[0-9]+_([0-9]+)', os.path.basename( path )).groups()[0]

def getrunFromPath( path ):
    if type(path) is list: path=path[0]
    #print( '(from %s'%path, )
    run = re.search(r'runs/(.+)/data', path ).groups()[0]
    #print( 'run %s)'%run, )
    return run

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
    #print( files[0] )
    pattern = re.search( r'(data_[0-9]+_to_[0-9]+)', files[0] ).groups()[0]
    #print( pattern )
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/nu_processing/scripts/ProcCat/*scn_osi_raw_gain_catalog_%s.root'%pattern) if not 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedCatalog( *paths ):
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/ext/catalog/catalog_data_*.root'%getrunFromPath(path)) if not 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedCatalogSkim( *paths ):
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/ext/catalog/catalog_data_*.root'%getrunFromPath(path)) if 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedOSI( *paths ):
    return [ osi for path in paths for osi in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/osi/images/osi_runID_*_%s_*.fits'%(getrunFromPath(path), getrunIDFromPath(path)) ) if not '-with-' in osi ]

#class ETA:
    #def __init__( self, function ):
        #self.function = function

    #def loop( self, arrays, args=[], loop_args = [], first = True ):
        #if len(loop_args) >= len(arrays):
            #t0 = time.time()
            #self.function( args+loop_args )
            
            #return
        #this_loop_args = loop_args + [0]
        #for i in arrays[len(loop_args)](*(args+loop_args) ):
            #this_loop_args[-1] = i
            #print( 'args', this_loop_args )
            #self.loop( arrays, args, this_loop_args )
        #return

def runETA( msg, cmd, eta = None, N = None, loop=None, verbose=None ):
    if verbose: print( term.yellow('%s <begin>'%msg) )
    if not loop is None:
        start = time.time()
        finishtimes = []
        for n, entry in enumerate(loop):
            print( term.yellow('entry %s <begin>'%entry) )
            etastart = time.time()
            try:
                cmd( entry )
            except ValueError as error:
                print( term.red('Warning: %s'%error) )
            et = time.time() - etastart
            print( term.green('entry %s <end '%entry + 't=%s>'%time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - etastart)) ) )
            finishtimes += [ time.time() + int(et*(len(loop) - n - 1)) ]
            print( term.yellow('estimated time %s'%time.strftime("%H:%M:%S", time.gmtime( int(np.mean(finishtimes)) - time.time() ) ) + '(finishes at %s)'%time.strftime("%H:%M:%S", time.localtime( int(np.mean(finishtimes)) )) ) )
        print( term.green('%s <end '%msg + 't=%s>'%time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start)) ) )
        return
    if not eta is None and not N is None:
        etastart = time.time()
        eta()
        et = time.time() - etastart
        print( term.yellow('estimated time %s'%time.strftime("%H:%M:%S", time.gmtime(int(et*N))) + '(finishes at %s)'%time.strftime("%H:%M:%S", time.localtime(  time.time() + int(et*N))) ) )
    start = time.time()
    res = cmd()
    if verbose: print( term.green('%s <end '%msg + 't=%s>'%time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start)) ) )
    return res

def readCatalog( ROOTfile, runID = None, verbose = False, readGain=True, selection='' ):
    min_max = lambda x: (np.min(x), np.max(x))
    if verbose: print( 'reading start and stop indexes for runID', runID, 'in ROOT file' )
    start = 0
    stop = 0
    if runID is not None:
        try:
            array = root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )
        except:
            return None
        indexlist = np.argwhere( array['runID'] == runID ).flatten()
        try:
            start, stop = min_max( indexlist )
        except ValueError:
            print( '!!!Error')
            print( ROOTfile )
            print( list(set(root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'] )) )
            print( 'indexlist', indexlist )
            raise
        if verbose: print( 'reading', len(indexlist), 'entries from', start, 'to', stop )
    
        readcatalog_once = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu','level','E0','E1','E2','E3','xMin','xMax','yMin','yMax','xBary0','yBary0','xBary2','yBary2','flag','n0','n2'] + (['gainCu'] if readGain else []), selection='runID==%s'%runID+selection, start=start_, stop=stop_+1 )
    else:
        readcatalog_once = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu','level','E0','E1','E2','E3','xMin','xMax','yMin','yMax','xBary0','yBary0','xBary2','yBary2','flag','n0','n2'] + (['gainCu'] if readGain else []), selection='runID==%s'%runID+selection, )

    hitscatalog = runETA( 
        msg = 'reading ROOTfile',
        cmd = lambda: readcatalog_once(start, stop), #read all entries and return 
        verbose = verbose,
        )
    if hitscatalog is None:
        print( 'failed to read' )
    return hitscatalog

def readCatalogConfig( ROOTfile ):
    print( 'config branches', root_numpy.list_branches( ROOTfile, treename = 'config' ) )
    config = root_numpy.root2array( ROOTfile, treename = 'config', branches = ['startTime','stopTime'] ).flatten()
    return config

def mergeParts( parts ):
    merged = parts[-1][:]
    #print( [ (i,ohdu.header['OHDU']) for i, ohdu in enumerate(merged) ] )
    for ohduindex in range(len(merged))[::-1]:
        if not merged[ohduindex].data is None:
            merged[ohduindex].data = np.concatenate( [ part[ohduindex].data for part in parts[::-1] ], axis = 0 )
        else:
            merged.pop(ohduindex)
    #print( [ (i,ohdu.header['OHDU']) for i, ohdu in enumerate(merged) ] )
    return merged

def removeHitsFromrunID( runID, outputfolder=None, ohdu=None, ROOTfile=None, osi=None, verbose=False, image=True, output=True, dohits=True, gain=None, crosstalk=False, onlyoccupancy=False ):
    FITSfile = listFITSscn_runID( runID )[0]
    print( 'input FITSfile =', FITSfile )
    if verbose: print( 'output flag =', output )
    if output:
        outputfilenohits = 'hsub_'+os.path.basename(FITSfile)
        if dohits: outputfilehits = 'hits_'+os.path.basename(FITSfile)
    
        if not outputfolder is None:
            if not os.path.exists(outputfolder):
                print( 'creating directory', outputfolder )
                os.makedirs( outputfolder )
        else:
            outputfolder = os.path.dirname(FITSfile)
            if os.access( outputfolder, os.W_OK):
                print( 'outputfolder =', outputfolder )
            else:
                outputfolder = 'run'+getrunFromPath( FITSfile ) + '/'
                print( 'outputfolder =', outputfolder )
                if not os.path.exists(outputfolder):
                    print( 'creating directory', outputfolder )
                    os.makedirs( outputfolder )
    
    removeOHDU = lambda x, ohdu_: map( lambda i: x.pop(i), [ j for j, hdu in enumerate(x) if not int(hdu.header['OHDU']) == int(ohdu) ][::-1] )
    if osi:
        OSIfiles = getAssociatedOSI( FITSfile )
        print( 'OSI file =' )
        print( '\n'.join( OSIfiles ) )
    
        if verbose: print( 'reading and merging OSI part files' )
        merged = mergeParts( [ astropy.io.fits.open( OSIfile ) for OSIfile in OSIfiles ] )
        if verbose: print( 'shape of merged', merged[-1].data.shape )
        if not ohdu is None:
            removeOHDU( merged, ohdu )
    
    if verbose: print( 'reading SCN file' )
    hdulist = {}
    
    if dohits:
        hdulist['hits'] = astropy.io.fits.open(FITSfile)
        if not ohdu is None:
            removeOHDU( hdulist['hits'], ohdu )

    hdulist['nohits'] = astropy.io.fits.open(FITSfile)
    #print( [ (_1,_2) for _1,_2 in hdulist['nohits'][-1].header.items() ] )
    ohdus = [ hdulist['nohits'][i].header['OHDU'] for i in range(len(hdulist['nohits'])) ]
    print( 'ohdus', ohdus )
    badOHDU = [11,12]
    #ohdus = [ ohdu for ohdu in ohdus if not ohdu in badOHDU ]
    #print( 'ohdus*', ohdus )
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
            print( ROOTfileGain )
        except:
            readGain = False
            ROOTfileGain = ROOTfile
    
    ROOTfile = ROOTfileGain
    #print( root_numpy.list_branches( ROOTfileGain, treename = 'hitSumm' ) )
    #print( root_numpy.root2array( ROOTfileGain, treename = 'hitSumm', branches = ['gainCu'] ) )
    print( 'input ROOTfile =', ROOTfile )
    hitscatalog = readCatalog( ROOTfile, runID=runID, readGain=readGain, verbose=verbose)

    if verbose: print( 'removing hits' )
    
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
            #print( 'len', len(selections[0]) )
            selections = [select_all] + list(selections[0])
            #print( 'len', len(selections) )
            for selection in selections:
                #print( selection.shape )
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
        print( 'energy_adu.shape', energy_adu.shape )
        print( 'adu_keV.shape', adu_keV.shape )

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
        print( 'min length_xy_pixel', min( length_xy_pixel[ length_xy_pixel > 0 ]) )
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
        #print( 'min angle', azimuthal_angle_degrees.min() )
        min_azimuthal_angle_selection = azimuthal_angle_degrees > azimuthal_angle_degrees.min()
        #print( 'min angle', azimuthal_angle_degrees[min_azimuthal_angle_selection].min() )
        
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
        
        print( 'selection:' )
        length_str = 'sqrt( ( xMax - xMin + 1 )^2*15 + ( yMax - yMin + 1 )^2*15 + 675^2 )'
        var_str = 'sqrt( xVar3^2 + yVar3^2 )'
        print( 'flag==0 && E3 > 0 && E0 > 0 && n3 > 0 && n0 > 0 && xVar3 > 0 && xVar0 > 0 && yVar3 > 0 && yVar0 > 0 '\
            +'&& ( log10(E3/E0) - %.4f )^2/%.4f^2 < %.4f ' %(m, s, n_std)\
            +'&& ( log10(xVar0/xVar3) %+.4f )^2/%.4f^2 < %.4f ' %(-m1, s1, n_std)\
            +'&& ( log10(yVar0/yVar3) %+.4f )^2/%.4f^2 < %.4f ' %(-m2, s2, n_std)\
            +'&& ( log10(E3/%s) %+.4f )^2/%.4f^2 < %.4f ' %(length_str, -m3, s3, n_std)\
            +'&& ( log10(%s/n3) %+.4f )^2/%.4f^2 < %.4f ' %(length_str, -m4, s4, n_std)\
            +'&& ( log10(E3/n3) %+.4f )^2/%.4f^2 < %.4f ' %( -m5, s5, n_std)\
            +'&& ( log10(n3/n0) %+.4f )^2/%.4f^2 < %.4f ' %( -m6, s6, n_std)\
            +'&& ( log10(%s/%s) %+.4f )^2/%.4f^2 < %.4f ' %(var_str, length_str, -m7, s7, n_std) )

        
        
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
        print( 'printing to', filename(0) )
        fig.savefig(filename(0))
        dist_bary_gcenter_pixel
        fig2 = plt.figure()
        
        lt5 = azimuthal_angle_degrees>5
        make_histogram(fig2, 111, azimuthal_angle_degrees[lt5], r'angle[degree]', 90, map( lambda x:x[lt5], selections[-1:] ), xscale='linear', yscale='linear', fitfunc = lambda x,A: 3*A*np.sin(x*rad_degrees)*np.cos(x*rad_degrees)**2 )
        
        print( 'printing to', filename(1) )
        fig2.savefig(filename(1))
        exit(0)

    for iohdu, thisohdu in enumerate( ohdus ):
        print( 'mask ohdu', thisohdu,  )
        hitsmask[thisohdu] = {}
        occupancy[thisohdu] = {}
        hitsohdu = hitscatalog[ hitscatalog['ohdu'] == thisohdu ]
        if len(hitsohdu['xPix']) > 0:
            x = np.concatenate( hitsohdu['xPix'], axis=0)
            y = np.concatenate( hitsohdu['yPix'], axis=0)
            l = np.concatenate( hitsohdu['level'], axis=0)
            if readGain:
                gain[thisohdu] = np.mean( hitsohdu['gainCu'] )
                print( 'gain', gain[thisohdu], 'adu/keV' )
            N = 0
            totalN = 0
            if gain[thisohdu] == 0:
                continue
            totalpix = len( hdulist['hits'][-1].data.flatten() )
            Ntiny = np.all( [hitsohdu['E0']/gain[thisohdu] < 0.075], axis = 0).sum()
            Nlow = np.all( [hitsohdu['E0']/gain[thisohdu] > 3, hitsohdu['E0']/gain[thisohdu] < 7], axis = 0).sum()
            Nhigh = np.all( [hitsohdu['E0']/gain[thisohdu] > 250, hitsohdu['E0']/gain[thisohdu] < 350], axis = 0).sum()
            print( Nlow, Nhigh )
            
                    
            for level in range(4):
                rate = 0
                totalN += (l==level).sum()
                hitsmask[thisohdu][level] = np.zeros_like( hdulist['hits'][-1].data, dtype=bool )
                hitsmask[thisohdu][level][y[l==level],x[l==level]] = True
                N += hitsmask[thisohdu][level].sum()
                rate = float(N)/totalpix
                #print( totalN, N )
                occupancy[thisohdu][level] = {'S': totalpix, 'N': N, 'rate': rate, 'gain': gain[thisohdu], 'totalN': totalN, 'N3-7': Nlow, 'N250-350': Nhigh, 'N.075': Ntiny }
        else:
            print( 'ohdu', thisohdu, 'empty' )
    sep=' '
    print( '%s/runID_%s_occupancy.csv'%(outputfolder,runID) )
    np.savetxt('%s/runID_%s_occupancy.csv'%(outputfolder,runID), [ [ runID, ohdukey, levelkey, data['S'], data['totalN'], data['N'], data['rate'], data['N3-7'], data['N250-350'], data['N.075'] ] for ohdukey, thisohdu in occupancy.items() for levelkey, data in thisohdu.iteritems() ], header=sep.join(['runID', 'ohdu','level', 'S','totalN','N','rate','N3-7','N250-350','N.075']), fmt='%s', delimiter=sep)
    
    if onlyoccupancy: return None
    
    for iohdu in range(len(hdulist['nohits'])):
        ohdu_ = hdulist['nohits'][iohdu].header['OHDU']
        print( 'crosstalk subtraction' )
        for iohdu2, thisohdu in enumerate( [ ohdu for ohdu in ohdus if not ohdu in badOHDU] ): #[:4]
            print( ohdu_, thisohdu, '%e'%np.mean(hdulist['nohits'][iohdu].data[ hitsmask[thisohdu][0] ]) )
    
    for iohdu in range(len(hdulist['nohits'])):
        ohdu_ = hdulist['nohits'][iohdu].header['OHDU']
        print( 'from ohdu', ohdu_, 'remove' )
        if verbose: print( 'removing hits from ohdus', ohdu_ )
        hitsohdu = hitscatalog[ hitscatalog['ohdu'] == ohdu_ ]
        if dohits: hdulist['hits'][iohdu].data[:,:] = REMOVE #-1e9

        hdulist['nohits'][iohdu].data *= 1e3/gain[ohdu_]/eIonization #e- (electron unit)
        
        hdulist['nohits'][iohdu].data[ hitsmask[ohdu_][3] ] = REMOVE
        if osi: merged[iohdu].data[ hitsmask[ohdu_][3] ] = REMOVE
        
        #if dohits: hdulist['hits'][iohdu].data[ hitsmask[ohdu_] ] = hits_e*1e3/gain/eIonization
        
        #if crosstalk:
            #print( 'crosstalk subtraction' )
            #for iohdu2, thisohdu in enumerate( [ ohdu for ohdu in ohdus if not ohdu in badOHDU] ): #[:4]
                #print( 'from ohdu', ohdu_, 'removing', thisohdu )
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
        print( 'ccd image', path )
        figs[0].savefig( path, bbox_inches='tight' )
        figs[1].savefig( path2, bbox_inches='tight' )

    if output:
        if dohits:
            hdulisthitsfile = outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilehits
            if os.path.exists(hdulisthitsfile):
                os.remove(hdulisthitsfile)
            if verbose: print( 'writing SCN hits' )
            hdulist['hits'].writeto( hdulisthitsfile )
        
        hdulistnohitsfile =  outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits 
        if os.path.exists(hdulistnohitsfile):
            os.remove(hdulistnohitsfile)
        if verbose: print( 'writing SCN no hits'     )
        hdulist['nohits'].writeto( hdulistnohitsfile )
        
        if osi:
            partsfile = outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits.replace('scn_mbs_', '')
            if os.path.exists(partsfile):
                os.remove(partsfile)
            if verbose: print( 'writing merged OSI no hits'     )
            merged.writeto( partsfile )
        
        print( 'created new FITS files' )
        print( hdulisthitsfile )
        print( hdulistnohitsfile )
        if osi: print( partsfile )
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
            print( __tab*indent+key )
            printDict( value, indent+1, skip=skip)
        else:
            print( __tab*indent+'%s: %s'%(key,value) )

def printDictTable( d, indent = 0, skip = None, header = None, line = None ):
    tab = ' '*4
    doheader = True
    for key, value in d.items():
        if key in skip: continue
        if type(value) is dict or type(value) is OrderedDict:
            printDictTable( value, indent+1, skip=skip, header = filter( lambda key: not key in skip, value.keys() ), line = '%s:%s'%(line, key) if line else key )
        else:
            print( ' '*len(line), ' '.join(map( lambda s: '%8s'%s, header)) )
            print( line, ' '.join([ '%.2e'%(d[col]) for col in header]) )
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

def poisson_( x, lamb ):
    y = np.zeros_like(x)
    x += lamb
    y[x>0] = np.exp( x[x>0]*np.log(lamb) -lamb -scipy.special.loggamma(x[x>0]+1) )
    return y
    
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
    if i == Nmax-1: print( 'Warning: Nmax reached in convolution' )
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
    
    #print( 'lambda', lamb )
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
    if i_ == Nmax-1: print( 'Warning: Nmax reached in convolution' )
    return A*np.exp(-lamb)*sp
    
def fitpartial( x, data, sel = None, log=False, p0 = None, adjust = None, bounds = (0.,np.inf), **kwargs ):
    #fix = filter( lambda k: not k[0] in adjust, p0.items() )
    p0_ = odict(p0)
    for key, value in kwargs.items():
        if key in p0_: p0_[key] = value
    
    p_ = [ v for k,v in p0_.items() if k in adjust ]
    #print( p_ )
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
    
    #print( 'chi2', chisq3, np.sum( (hist3-y_)**2/y_ )/(len(y_) - len(popt)) )
    #print( np.max(hist3), np.max(y_) )
    
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
    
    #print( 'CF', pp, chisq2, perr )
    #print( 'LS', popt, scipy.stats.chisquare( f(x, *popt), data, ddof = len(popt) )[0]/(len(data) - len(popt)), perr2 )
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
    #print( pp )
    #print( popt )
    return {'chisqr': chisq2, 'params': popt, 'hist': hist2, 'perr': perr }

def curve_fit( func, xdata, ydata, p0 ):
    residual = lambda p, x, y: (y - func(x, *p))
    leastsqReturn = scipy.optimize.leastsq( residual, p0, args=(xdata,ydata), full_output=1 )
    popt, pcov2, infodict, errmsg, success = leastsqReturn
    
    s_sq = ( residual(popt, xdata, ydata )**2).sum()/(len(ydata) - len(popt))
    pcov2 *= s_sq
    perr = [ np.sqrt( np.abs(pcov2[i][i]) ) for i in range(len(popt)) ]
    return popt, perr, pcov2

def computeFits( ohdu, runID, hdulists, verbose=False, plot=False, gain=None ):
    result = []
    header = []
    fmt = []
    first = True
    
    clean = lambda x: x.flatten()[ np.all( [x.flatten()!=REMOVEOLD, x.flatten()!=REMOVE, x.flatten()!=1e10], axis=0) ]

    hdulists = odict( [('scn', hdulists)] )
    if verbose: print( 'ohdus in file', ', '.join( map( str, [ hdu.header['OHDU'] for hdu in hdulists['scn'] ] )) )
    
    data = odict([ (key, hdu.data) for key, hdulist in hdulists.items() for hdu in hdulist if hdu.header['OHDU'] == ohdu ])
    if verbose: 
        print( 'keys', data.keys() )
        print( 'shapes', data['scn'].shape )
        print( 'data extremes', data['scn'].min(), data['scn'].max() )
    yccd = 4262
    xccd = 4220
    yos = 4112
    #if osi:
        #if data['osi'].shape[1] == 2*yccd:
            #print( 'subtracting idle amplifier' )
            #data['osisub'] = data['osi'][0:xccd,0:yccd] - np.flip(data['osi'][0:xccd,yccd:2*yccd], 1)
    #if scnsub:
        #print( 'subtracting idle overscan' )
        #print( 'bias correction shape', np.median(data['scn'][0:xccd,yos:], axis=1).shape )
        #data['scnsub'] = data['scn'][0:xccd,:] - np.median(data['scn'][0:xccd,yos+10:-10], axis=1)[:,np.newaxis]
    #print( 'data len', data.keys() )
    dx = 1
    N = 500
    shift = lambda bins: ( bins[1:] + bins[:-1] )*.5
    #bins = np.r_[min(clean(data['scn'])):max(clean(data['scn'])):dx]
    maxSCN = max(clean(data['scn']))
    bins = np.linspace(-maxSCN, maxSCN, N)
    dx = 2*maxSCN/N
    if verbose: print( 'dx', dx )
    if verbose: print( 'bins', min(bins), max(bins) )

    if verbose: print( 'computing histogram of full ccd' )
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
    if verbose: print( 'keys', [ (datakey, hist[datakey].keys()) for datakey in data.keys() ] )
    
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
            print( 'pixelcount', datum['pixelcount'], 'total', datum['totalpixel'], 'occupancy', 1.-float(datum['pixelcount'])/datum['totalpixel'] )
            x = datum['hist']['x']
            y = datum['hist']['y']
            datum['mean2'] = np.mean(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            datum['median2'] = np.median(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            datum['std2'] = np.std(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            datum['N2'] = len(crop(x,clean(data[datakey][vslice,hslices[slicekey]])))
            N = y.sum()*dx
            N_ = N
            #print( 'µ', 'σ', 'λ' )
            #print( datum['pixelcount'] )
            #print( y.sum()/dx )
            #print( 2*y[x<0].sum()/dx )

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

            if verbose: print( 'std2', slicekey, datum['std2'] )
            fits = datum['fits']

            std2os = thishist['os']['std2']
            fits[r'G(0,σ\')'] = fitpartial( x, y, p0=p0, adjust = ['A'], sigma_os = std2os )
            Afit = fits[r'G(0,σ\')']['params']['A']
            fits[r'G(0,σ)A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A'], sigma_os = std2os, A = Afit )
            Afit2 = fits[r'G(0,σ)A']['params']['A']
            fits[r'G(µ,σ)A'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'A'], sigma_os = std2os, A = Afit )
            #print( N, Afit )
            #if slicekey == 'os':
                #fits[r'G(0,σ)A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A'], sigma_os = std2os, A = Afit )
                #fits[r'log(G(0,σ))'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os'], log=True, sigma_os = std2os, A = Afit )
                #fits[r'log(G(0,σ))A'] = fitpartial( x, y, p0=p0, adjust = ['sigma_os', 'A'], log=True, sigma_os = std2os, A = Afit )
                #pass
                
            if slicekey == 'ac':
                lamb0 = thishist['ac']['std2']**2 - thishist['os']['std2']**2
                if verbose: print( 'lambda', lamb0 )
                #fits[r'G(µ,σ\')*P(λ)'] = fitpartial( x, y, p0=p0, adjust = ['mu_ac', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=2./3 )
                #fits[r'G(0,σ\')*P(λ\')A'] = fitpartial( x, y, p0=p0, adjust = ['lamb'], sigma_os=std2os, A=Afit2, lamb=lamb0 )
                fits[r'G(0,σ\')*P(λ\')'] = fitpartial( x, y, p0=p0, adjust = ['A'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                fits[r'G(µ,σ)*P(λ)A'] = fitpartial( x, y, p0=p0, adjust = ['mu_os', 'sigma_os', 'A', 'lamb'], sigma_os=std2os, A=Afit, lamb=lamb0, gain=gain )
                Afit = fits[r'G(0,σ\')*P(λ\')']['params']['A']
                Afiterr = fits[r'G(0,σ\')*P(λ\')']['perr']['A']
                #if verbose: print( 'fiterr', Afit, Afiterr, str_with_err(Afit, Afiterr) )

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
        print( sep.join( [ str((i,h)) for i, h in enumerate(header)] ) )
        print( sep.join(header) )
        print( '\n'.join( [ sep.join( [ f%entry for f,entry in zip(fmt,line) ] ) for line in result ] ) )
    printcsv = lambda outfile: np.savetxt('%s.csv'%(outfile), result, header=sep.join(header), fmt='%s', delimiter=sep)
    printcsv2 = lambda outfile: np.savetxt('%s.csv2'%(outfile), [line2], header=sep.join(header2), fmt='%s', delimiter=sep)
    return hist, printcsv, printcsv2

def makePlot( ohdu, runID, hist, verbose = False ):
    nfiles = len(hist.keys())
    slices = set( [ k for v in hist.values() for k in v.keys() ] )
    print( 'slices', slices )
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
    
    if verbose: print( 'generating plots' )
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
                print( y.shape )
                ly = np.log( y )
                #axfft[key][slicekey].step( datum['hist']['x'], np.imag(ly), where='mid', label='i' )
                #lamb = np.mean( np.real(ly) )
                print( slicekey )
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

    if verbose: print( 'setting labels and lims' )
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
    
    print( 'FITSfiles:\n', '\n'.join(FITSfiles) )
    runID = getrunIDFromPath( FITSfiles[0] )
        
    hist, printcsv, printcsv2 = computeFits( ohdu, runID, astropy.io.fits.open( FITSfiles[0] ), verbose, gain )
    print( 'tables generated in CSV file', '%s.csv(2)'%(FITSfiles[0]) )
    printcsv(FITSfiles[0])
    #print( 'table generated in CSV file', '%s.csv2'%(FITSfiles[0]) )
    printcsv(FITSfiles[0])

    if plot:
        saveplot = makePlot( ohdu, runID, hist )
        print( 'plot generated in PNG file', '%s.png'%(FITSfiles[0]) )
        saveplot(FITSfiles[0])

def plot2( outfolder, ohdu, runID, ROOTfile = None, plot=False, gain=None, verbose=True, crosstalk=False, onlyoccupancy=False ):
    if verbose: print( 'runID', runID )
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
        print( 'creating directory', outfolder )
        os.makedirs( outfolder )
    print( 'tables generated in CSV file', '%s.csv(2)'%(outputfile) )
    printcsv(outputfile)
    #print( 'table generated in CSV file', '%s.csv'%(outputfile) )
    printcsv2(outputfile)
    
    if plot:
        saveplot = makePlot( ohdu, runID, hist )
        print( 'plot generated in PNG file', '%s.png'%(outputfile) )
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
    
    print( len(x), ymin, ymax )
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
        print( i%m, i/m, n, m )
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
    #z = track['ePix']/track['gainCu'][0]
    z = track['ePix']/track['gainCu']
    
    x -= x.min()
    y -= y.min()
    
    matrix = np.zeros([ x.max()+1, y.max()+1 ])
    print( 'matrix.shape', matrix.shape )
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
    #print( 'fields', a.dtype.descr )
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
    #print( 'keys', keys )
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
    print( number_of_entries )
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
        print( 'category selected', char )
        if char == 'm':
            print( 'muon' )
            updatejson( 'muon.json', array2dict(entry) )
            continue
        if char == 'n':
            print( 'not muon' )
            updatejson( 'notmuon.json', array2dict(entry) )
            continue
        if char == 'list':
            entries = json.load( open('muon.json','r') )
            print( 'number of muons', len(entries) )
            print( 'number of notmuons', len(json.load( open('notmuon.json','r') )) )
            continue
        if char == 'q':
            print( 'quit' )
            break
    return

def pca( data_, inv = None, method='invert' ):
    data = data_.copy()
    if not inv is None: 
        datainv = inv.copy()
    #print( 'data.shape', data.shape )
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
    print( eigenvalues.shape, eigenvectors.shape )
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
        #print( 'len', len(ix), )
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
            print( 'error' )
        else:
            t1Var0 = std(it[il<0],ie[il<0])
            l1 = np.mean(il[il<0])
            t2Var0 = std(it[il>0],ie[il>0])
            l2 = np.mean(il[il>0])
            #print( 'tVar', tVar0[i], t1Var0, t2Var0, l1, l2, il.min(), il.max() )
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
    print( len(positivedata) )
    print( len(negativedata) )
    print( 'positivedata.descr', positivedata.dtype.descr )
    
    extrapos = muonPixFraction( positivedata['xPix'], positivedata['yPix'], positivedata['ePix'], positivedata['gainCu'] )
    #print( extrapos['cChiSqr'] )
    positivedata = append_struct( positivedata, extrapos )
    extraneg = muonPixFraction( negativedata['xPix'], negativedata['yPix'], negativedata['ePix'], positivedata['gainCu'] )
    negativedata = append_struct( negativedata, extraneg )
    fields_ = [ field for field in positivedata.dtype.names if not field in ['xPix', 'yPix', 'ePix', 'level', 'flag', 'hPixFlag', 'osMadSCN', 'xBary0', 'xBary1', 'xBary2', 'xBary3', 'xMax', 'xMin', 'yBary0', 'yBary1', 'yBary2', 'yBary3', 'yMax', 'yMin', 'index', 'catalog', 'osMadRaw', 'dcVar', 'selFlag', 'runID', 'ohdu', 'expoStart', 'resCu', 'gainCu', 'E1', 'n1', 'E2', 'E3', 'n2', 'n3', 'xVar1', 'xVar2', 'xVar3', 'yVar1', 'yVar2', 'yVar3', 'nSat' ] ]
    print( 'fields', fields_ )
    #fields_ = ['E0', 'lLen', 'cChiSqr', 'meanE' ]
    #predim = [ ('log10(%s)'%(key), lambda x: np.log10(x[key]) ) for key in fields_ ]
    #predim += [ ('log10(%s)/log10(%s)'%(key,key2), lambda x: np.log10(x[key])/np.log10(x[key2]) ) for key in fields_ for key2 in fields_ ]
    #print( predim, len(predim) )
    dimensions1 = dict( ('log10(%s)'%(key), lambda x,key=key: np.log10(x[key]) ) for key in fields_ )
    #dimensions2 = dict( ('%s'%(key), lambda x,key=key: x[key] ) for key in fields_ )
    #dimensions2 = dict( ('log10(%s*%s)'%(key,key2), lambda x,key=key,key2=key2: np.log10(x[key])+np.log(x[key2]) ) for i, key in enumerate(fields_) for j, key2 in enumerate(fields_) if j>i )
    #a
    #print( dimensions, len(dimensions) )
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

    #print( [ (key) for key in dimensions.keys() ] )
    datapos = np.array( zip( *[ dimension(positivedata) for dimension in dimensions.values() ] ), dtype = [ (key, None) for key in dimensions.keys() ] )
    dataneg = np.array( zip( *[ dimension(negativedata) for dimension in dimensions.values() ] ), dtype = [ (key, None) for key in dimensions.keys() ] )
    #print( 'isnan', [ (key, np.argwhere( np.isnan(dataneg[key]) ) ) for key in dimensions.keys() ] )
    #print( 'isinf', [ (key, np.argwhere( np.isinf(dataneg[key]) ) ) for key in dimensions.keys() ] )
    #print( datapos[146] )

    #print( datapos.shape )
    #datapos = np.delete( datapos, np.argwhere( datapos['log10(E0/Var0)'] > 3 ) )
    #print( datapos.shape )
    
    meanpos = apply_to_structured( datapos, lambda x: [np.mean(x)] )
    meanneg = apply_to_structured( dataneg, lambda x: [np.mean(x)] )
    #print( np.vstack( (meanpos, meanneg))['E0'] )
    
    #pcacluster = pca( np.vstack( (meanpos, meanneg)) )
    #print( 'eigenvalues cluster', list(enumerate(pcacluster[0])) )
    
    pcapos = pca(datapos)
    print( 'eigenvalues pos', list(enumerate(pcapos[0]))[-3:] )

    pcaneg = pca(dataneg)
    print( 'eigenvalues neg', list(enumerate(pcaneg[0]))[-3:] )

    pcapos_ineg = pca( datapos, inv = dataneg, method = 'subtract')
    print( 'eigenvalues pos_ineg', list(enumerate(pcapos_ineg[0]))[-3:] )

    pcaneg_ipos = pca( dataneg, inv = datapos, method = 'subtract' )
    print( 'eigenvalues neg_ipos', list(enumerate(pcaneg_ipos[0]))[-3:] )

    #pca_ = pcaneg_ipos
    pca_ = pcapos_ineg
    #pca_ = pcacluster
    #pca_ = pcaneg
    #print( 'eigenvalues', list(enumerate(pca_[0])) )
    print( 'eigenvectors0', sorted( [ (dimensions.keys()[i], abs(v)) for i, v in enumerate(pca_[1][:,-1]) ], key= lambda x: x[-1] )[-4:] )
    print( 'eigenvectors1', sorted( [ (dimensions.keys()[i], abs(v)) for i, v in enumerate(pca_[1][:,-2]) ], key= lambda x: x[-1] )[-4:] )
    print( 'eigenvectors2', sorted( [ (dimensions.keys()[i], abs(v)) for i, v in enumerate(pca_[1][:,-3]) ], key= lambda x: x[-1] )[-4:] )
    rpos = rotatedata( datapos, *pca_ )
    rpos_pos = rotatedata( datapos, *pcapos )
    rpos_neg = rotatedata( dataneg, *pcapos )
    rneg = rotatedata( dataneg, *pca_ )
    rneg_pos = rotatedata( datapos, *pcaneg )
    rneg_neg = rotatedata( dataneg, *pcaneg )
    
    #centers, radii = findclusters( rpos )
    #print( 'eigenvalues', pca_[0] )
    centers = apply_to_structured( rpos, lambda x: [np.mean(x)] )[0]
    stds = apply_to_structured( rpos, lambda x: [3*np.std(x)] )[0]
    #print( 'centers', centers )
    #print( 'stds', stds )
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

    #print( [ (key) for key in dimensions.keys() ] )
    datadim = np.array( zip( *[ dimension(data) for dimension in dimensions.values() ] ), dtype = [ (key, None) for key in dimensions.keys() ] )
    rdata = rotatedata( datadim, *pca_ )
    
    mask = np.all( [ rdata[field] > bottoms[i] for i, field in enumerate(rdata.dtype.names) ] + [ rdata[field] < tops[i] for i, field in enumerate(rdata.dtype.names) ], axis = 0 )
    #rdata = rdata[ mask ]
    #data = data[ mask ]
    return data, rdata, data[mask], rdata[mask]
    

def updateCatalog( outputcatalog, inputcatalog ):
    interactive_selection( inputcatalog, outputcatalog )

def help_():
    print( 'Usage:' )
    print( sys.argv[0], '--runID <runID> --outfolder <path> --removehits' )
    print( '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in the associated ROOT catalog' )
    print
    print( sys.argv[0], '--runID <runID> --outfolder <path> --ROOTfile <path> --removehits' )
    print( '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in given ROOTfile' )
    print(  )
    print( sys.argv[0], '--ohdu <ohdu> --fft --plot <list of paths>' )
    print( '\tgenerate plot with FITS files from paths for the given ohdu' )
    print
    print( sys.argv[0], '--ohdu <ohdu> --runID <runID> --outfolder <outfolder> --gain <gain> --table2' )
    print( '\tgenerate plots and table of function fits for the respective runID and given ohdu' )
    print
    print( sys.argv[0], '--ohdu <ohdu> --run <run> --outfolder <outfolder> --gain <gain> --table2' )
    print( '\tgenerate tables of function fits for all the runIDs in the given run and given ohdu' )
    print

def getOption( vars_, option, args, f = lambda x:x, flag=False ):
    if '--%s'%option in args:
        vars_[option] = f(args[args.index('--%s'%option)+1]) if not flag else True
        print( 'setting %s ='%option, vars_[option] )
        return
    vars_[option] = None
    return

def processImage3( image, threshold ):
    image[ image == 1e10 ] = 0

    t0 = time.clock()
    s33 = [[1,1,1],
        [1,1,1],
        [1,1,1]]
    #s33 = [[0,1,0],
        #[1,1,1],
        #[0,1,0]]
    scicluster, num = scipy.ndimage.label( image >= threshold, structure=s33 )
    print( 't/px %.2es'%( float(time.clock()-t0)/len(image.flatten()) ), )
    print( 'N3', len( np.unique(scicluster.flatten()) ), )
    Emax = scipy.ndimage.maximum( image, scicluster, index=np.arange(0,np.max(scicluster)+1) )
    print( 'Emax', np.max(Emax), Emax.shape, )
    if np.max(Emax) < 4*threshold:
        return None
    scicluster[ Emax[ scicluster ] < 4*threshold ] = 0
    labels = np.unique(scicluster.flatten())
    labels = labels[ labels > 0 ]
    print( 'N0', len(labels), )
    
    def xBary( val, pos ):
        print( 'pos', pos, pos.shape )
        return [ np.average( np.unravel_index(pos, image.shape)[0], weights=val, axis=0 ), 1 ]
    def yBary( val, pos ):
        return np.average( np.unravel_index(pos, image.shape)[1], weights=val, axis=0 )
    
    def measures( val, pos ):
        n = val.size
        E = np.sum(val)
        r = np.unravel_index(pos, image.shape)
        xBary = np.average( r[0], weights=val, axis=0 )
        yBary = np.average( r[1], weights=val, axis=0 )
        xVar = np.average( r[0]**2, weights=val, axis=0 ) - xBary**2
        yVar = np.average( r[1]**2, weights=val, axis=0 ) - yBary**2
        return n, E, xBary, yBary, xVar, yVar
    def label_comprehension( input, clusters, func ):
        labels = np.unique(clusters.flatten())
        labels = labels[ labels > 0 ]
        positions = np.arange(input.size).reshape(input.shape)
        return zip( *[ func( input[ clusters == i ], positions[ clusters == i ] ) for i in labels ] )
    def _xVar( val, pos ):
        r = np.unravel_index(pos, image.shape)
        return np.average( r[0]**2, weights=val ) - np.average( r[0], weights=val )**2
    def _yVar( val, pos ):
        r = np.unravel_index(pos, image.shape)
        return np.average( r[1]**2, weights=val ) - np.average( r[1], weights=val )**2
    def _rVar( val, pos ):
        r = np.unravel_index(pos, image.shape)
        r0 = np.average( r[0], weights=val ), np.average( r[1], weights=val )
        return np.average( (r[0]-r0[0])**2 + (r[1]-r0[1])**2, weights=val )
    
    L = 4
    n = {}
    E = {}
    xBary = {}
    yBary = {}
    xVar = {}
    yVar = {}
    rVar = {}
    for i in range( 1, L+1 ):
        image[image<i*threshold] = 0
        j = L-i
        n[j] = scipy.ndimage.labeled_comprehension( image, scicluster, index=labels, func=np.count_nonzero, out_dtype=int, default=0 )
        E[j] = scipy.ndimage.sum( image, scicluster, index=labels )
        xBary[j], yBary[j] = zip( *scipy.ndimage.center_of_mass( image, scicluster, index=labels ) )
        xVar[j] = scipy.ndimage.labeled_comprehension( image, scicluster, index=labels, func=_xVar, out_dtype=int, default=0, pass_positions=True )
        yVar[j] = scipy.ndimage.labeled_comprehension( image, scicluster, index=labels, func=_yVar, out_dtype=int, default=0, pass_positions=True )
        rVar[j] = scipy.ndimage.labeled_comprehension( image, scicluster, index=labels, func=_rVar, out_dtype=int, default=0, pass_positions=True )
    #n3, E3, xBary3, yBary3, xVar3, yVar3 = label_comprehension( image, scicluster, measures )

    print( 'done\n', )
    data = np.array( zip(*[ 
        E[3], E[2], E[1], E[0], 
        n[3], n[2], n[1], n[0], 
        xBary[3], yBary[3], xBary[2], yBary[2], xBary[1], yBary[1], xBary[0], yBary[0], 
        xVar[3], xVar[2], xVar[1], xVar[0], 
        yVar[3], yVar[2], yVar[1], yVar[0],
        rVar[3], rVar[2], rVar[1], rVar[0],
        ]), dtype=[ 
            ('E3', np.float), ('E2', np.float), ('E1', np.float), ('E0', np.float), 
            ('n3', np.int), ('n2', np.int), ('n1', np.int), ('n0', np.int), 
            ('xBary3', np.float), ('yBary3', np.float), ('xBary2', np.float), ('yBary2', np.float), ('xBary1', np.float), ('yBary1', np.float), ('xBary0', np.float), ('yBary0', np.float), 
            ('xVar3', np.float), ('xVar2', np.float), ('xVar1', np.float), ('xVar0', np.float), 
            ('yVar3', np.float), ('yVar2', np.float), ('yVar1', np.float), ('yVar0', np.float),
            ('rVar3', np.float), ('rVar2', np.float), ('rVar1', np.float), ('rVar0', np.float),] )
    #print( 'num', len(E3) )
    return data

def processImage4( image, threshold ):
    image[ image == 1e10 ] = 0

    t0 = time.clock()
    s33 = [[1,1,1],
        [1,1,1],
        [1,1,1]]

    cluster0, num0 = scipy.ndimage.label( image >= 4*threshold, structure=s33 )
    
    print( 't/px %.2es'%( float(time.clock()-t0)/len(image.flatten()) ), )
    print( 'N0', num0, )
    is_bg = cluster0 == 0
    distance_from_track = scipy.ndimage.distance_transform_edt( is_bg )
    
    def _xVar( val, pos ):
        r = np.unravel_index(pos, image.shape)
        return np.average( r[0]**2, weights=val ) - np.average( r[0], weights=val )**2
    def _yVar( val, pos ):
        r = np.unravel_index(pos, image.shape)
        return np.average( r[1]**2, weights=val ) - np.average( r[1], weights=val )**2
    def _rVar( val, pos ):
        r = np.unravel_index(pos, image.shape)
        r0 = np.average( r[0], weights=val ), np.average( r[1], weights=val )
        return np.average( (r[0]-r0[0])**2 + (r[1]-r0[1])**2, weights=val )
    def _flag( val, pos ):
        r = np.unravel_index(pos, image.shape)
        flag = 0
        if 0 in r[0]: flag=1
        if 0 in r[1]: flag=1
        if image.shape[0]-1 in r[0]: flag=1
        if image.shape[1]-1 in r[1]: flag=1
        return flag
    
    def measures( input, clusters, level ):
        labels = np.unique(clusters.flatten())
        labels = labels[ labels > 0 ]
        n = scipy.ndimage.labeled_comprehension( input, clusters, index=labels, func=np.count_nonzero, out_dtype=int, default=0 )
        E = scipy.ndimage.sum( input, clusters, index=labels )
        xBary, yBary = zip( *scipy.ndimage.center_of_mass( input, clusters, index=labels ) )
        xVar = scipy.ndimage.labeled_comprehension( input, clusters, index=labels, func=_xVar, out_dtype=int, default=0, pass_positions=True )
        yVar = scipy.ndimage.labeled_comprehension( input, clusters, index=labels, func=_yVar, out_dtype=int, default=0, pass_positions=True )
        rVar = scipy.ndimage.labeled_comprehension( input, clusters, index=labels, func=_rVar, out_dtype=int, default=0, pass_positions=True )
        flag = scipy.ndimage.labeled_comprehension( input, clusters, index=labels, func=_flag, out_dtype=int, default=0, pass_positions=True )
        return n, E, xBary, yBary, xVar, yVar, rVar, np.ones(len(n)).astype(int)*level, flag
    
    cluster1, num1 = scipy.ndimage.label( distance_from_track < 2, structure=s33 )
    print( num1, )
    cluster2, num2 = scipy.ndimage.label( distance_from_track < 3, structure=s33 )
    print( num2, )
    cluster3, num3 = scipy.ndimage.label( distance_from_track < 4, structure=s33 )
    print( num3, )

    dtype = [ ('n', np.int), ('E', np.float), ('xBary', np.float), ('yBary', np.float), ('xVar', np.float), ('yVar', np.float), ('rVar', np.float), ('level', np.int), ('flag', np.int) ]

    data0 = np.array( zip( *measures(image, cluster0, 0) ), dtype = dtype )
    data1 = np.array( zip( *measures(image, cluster1, 1) ), dtype = dtype )
    data2 = np.array( zip( *measures(image, cluster2, 2) ), dtype = dtype )
    data3 = np.array( zip( *measures(image, cluster3, 3) ), dtype = dtype )
    
    data = stack_arrays((data0, data1, data2, data3))
    print( 'done\n', )
    return data

def processImage2( fname, out=None, **config ):
    i = fname
    if out is None:
        o = os.path.splitext(fname)[0] + '.root'
    else:
        o = out
    if not os.path.exists(o):
        print( 'calling external cmd' )
        #os.system( '/share/storage2/connie/data_analysis/Processor/tools/extract/extract.exe %s -c /home/bonifazi/public/extractConfigFS.xml -o %s'%(i,o) )
        os.system( '/share/storage2/connie/data_analysis/Processor/tools/extract/extract.exe %s -c /home/mota/tests/simulatedFITS/extractConfigFS.xml -o %s'%(i,o) )
        print( 'ended' )
    
    catalogData = readCatalog( o, config['runID'], readGain=False )
    ohduData = catalogData[ catalogData['ohdu'] == int(config['ohdu']) ]
    return ohduData
    
def outputFITS( data, name ):
    print( 'outputFITS' )
    fname = 'simulatedFITS/%s.fits'%name
    if os.path.exists(fname):
        os.remove(fname)
    astropy.io.fits.writeto( fname, data )
    return

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

def writeImageFITS( data, header, outPath='./', tag=None, fname = None ):
    header_ = astropy.io.fits.Header()
    for key, value in header.items():
        #print( key, type(value) )
        if type(value) is list:
            if len(value) > 2: continue
            header_[key] = tuple( header[key] )
        else:
            header_[key] = value
    if fname is None:
        fname = outPath + '%s.fits'%tag
    if os.path.exists(fname):
        os.remove(fname)

    astropy.io.fits.writeto( fname, data, header = header_ )
    return fname
    
def readSingleImageFITS( outPath=None, tag=None, fname = None, i = 0 ):
    if fname is not None:
        fname = fname
    else:
        fname = outPath + '/%s.fits'%tag
    image = astropy.io.fits.open( fname, ignore_missing_end=True )
    return image[i].data, image[i].header

def adu2e( adu, adu_keV ):
    keV = adu/adu_keV
    e = keV*1e3/eV_e
    return e

def e2adu( e, adu_keV ):
    keV = eV_e*e/1e3
    adu = adu_keV*keV
    return adu

z2sigmaModel = {
    '2': lambda x: np.sqrt(-891.074*np.log(1-0.000488203*x))/15,
    #'3': lambda x: np.sqrt(-923.666*np.log(1-0.000441113*x))/15,
    '3': lambda x: np.sqrt(-923.666*np.log1p(-0.000441113*x))/15,
    '3ap': lambda x: np.sqrt(923.666)/15.*np.sqrt(0.000441113*x)*( 1+.25*0.000441113*x ),
    #ext4="sqrt(-820.692*TMath::Log(1-0.00049098*x))/15"
    #ext5="sqrt(-547.614*TMath::Log(1-0.000689114*x))/15"
    #ext6="sqrt(-597.398*TMath::Log(1-0.000633145*x))/15"
    #ext7="sqrt(-62770*TMath::Log(1-6.31885e-06*x))/15"
    #ext8="sqrt(-476.3*TMath::Log(1-0.000759706*x))/15"
    #ext9="sqrt(-565.109*TMath::Log(1-0.000647537*x))/15"
    #ext10="sqrt(-763.933*TMath::Log(1-0.000494423*x))/15"
    #ext11="sqrt(-1016.97*TMath::Log(1-0.000391788*x))/15"
    #ext12="sqrt(-907.968*TMath::Log(1-0.000436115*x))/15"
    #ext13="sqrt(-601.001*TMath::Log(1-0.000625386*x))/15"
    #ext14="sqrt(-556.409*TMath::Log(1-0.00068271*x))/15"
    #ext15="sqrt(-7060.89*TMath::Log(1-6.26018e-05*x))/15"
}

def diffuseElectrons( x, y, z ):
    sigma = z2sigmaModel['3']( z )
    #xoff, yoff = sigma * scipy.stats.norm.rvs( size = len(sigma) ), sigma * scipy.stats.norm.rvs( size = len(sigma) )
    xoff, yoff = sigma * scipy.stats.norm.rvs( size = 2*len(sigma) ).reshape(2,len(sigma))
    return x + xoff, y + yoff

def round2int2( x, y ):
    return np.rint( x ).astype(int), np.rint( y ).astype(int)

def round2int( *x ):
    return [ np.rint( ix ).astype(int) for ix in x ]

def writeArraysTo( labels=None, arrays=None, fmt=None, fname=None ):
    stream = open( fname, 'w' )
    stream.write( '#' + ' '.join(labels) + '\n' )
    for line in zip( *arrays ):
        stream.write( (' '.join(fmt) +'\n' )%line )
    stream.close()

round2next = lambda a,b: a if a%b == 0 else b*int(round(float(a)/b+1))
round2prev = lambda a,b: a if a%b == 0 else b*int(round(float(a)/b))

def parseConfig( configFile ):
    #print( json.load( open( configFile, 'r' ) ) )
    config = json.load( open( configFile, 'r' ) )
    outPath = os.path.splitext( configFile )[0]
    config['outPath'] = outPath + '/'
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    return config
    
class Config:
    def __init__(self, fname ):
        self.config = parseConfig( fname )
        for key, value in self.config.items():
            match = re.match(r"(.*)\[(.*)\]", key)
            if match:
                varname, varunit = re.match(r"(.*)\[(.*)\]", key).groups()
                setattr( self, '%s_unit'%varname, varunit )
            else:
                varname, varunit = key, None
            setattr( self, varname, value )
        return
    
    def __getitem__(self, a):
        return self.config[a]
    
    def __call__(self):
        return self.config

def generatePlots( fname, settings=[], plots=[] ):
    fig = plt.figure()
    n = len(plots)
    m = len(plots[0])
    print( 'plotting', n, m, fname )
    ax = [ [ None for col in line ] for line in plots ]
    for i, line in list(enumerate(plots)):
        for j, column in list(enumerate(line))[::-1]:
            if plots[i][j] is None: continue
            ax[i][j] = plt.subplot2grid( (n, m), (i,j), sharex = ax[0][j], sharey = ax[i][len(line)-1] )
            for plot in plots[i][j]:
                try:
                    plot( ax[i][j] )
                except ValueError as error:
                    traceback.print_exc()
            for setting in settings:
                setting( ax[i][j] )
    print( 'generating', fname )
    plt.subplots_adjust( hspace=0., wspace=0. )
    fig.savefig(fname)
    plt.close()

def peakFinding( x, data, radius, tol=1, num=None ):
    
    print( 'data', data.shape, radius )
    croppedData = data[radius:-radius]
    croppedx = x[radius:-radius]
    print( 'cropped', croppedData.shape )
    movingData = [ data[i:-2*radius+i] for i in range(2*radius) ]
    movingMean = np.mean( movingData, axis=0)
    movingStd = np.std( movingData, axis=0)
    ratio = np.abs(croppedData - movingMean)/movingStd
    if num is None:
        return croppedx[ np.argwhere( croppedData - movingMean > tol*movingStd ) ]
    else:
        return croppedx[ np.argsort( ratio )[::-1][:num] ]

def smooth2( data, radius, func = np.mean ):
    if radius == 0:
        return data
    startData = [ func( data[:radius+i] ) for i in range(radius) ]
    #print( 'len', len(startData) )
    centerData = func( [ data[i:i-2*radius if i < 2*radius else None] for i in range(2*radius+1) ], axis=0 )
    #print( 'len', len(centerData) )
    endData = [ func( data[i-2*radius:] ) for i in range(radius) ]
    #print( 'len', len(endData) )
    movingMean = np.concatenate( (startData, centerData, endData) )
    #print( 'smooth2', len(data), len(movingMean) )
    return movingMean

def find_peaks( data, radius=1, tol=1 ):
    if radius == 0:
        return data, np.zeros_like(data)
    movingMean = smooth2( data, radius, np.mean )
    movingStd = smooth2( data, radius, np.std )
    peaks = np.argwhere( data > tol*movingStd + movingMean )
    return peaks

xlabel = lambda s: ( lambda ax: ax.set_xlabel(s) )
ylabel = lambda s: ( lambda ax: ax.set_ylabel(s) )
ylog = lambda ax: ax.set_yscale('log')
xlog = lambda ax: ax.set_xscale('log')
lineSi = lambda ax: ax.axvline( 1.7, color='c' )
lineCu = lambda ax: ax.axvline( 8.05, color='m' )
legend = lambda ax: ax.legend(fancybox=True, framealpha=0)
grid = lambda ax: [ax.grid(True, which='major', linewidth=1, linestyle=':'), ax.grid(True, which='minor', linewidth=0.2, alpha=.75, linestyle=':')]
minorgrid = lambda ax: ax.minorticks_on()

def efficiencyNeutrinosReconstruction( c ):
    config = Config( c )
    configFileNu = os.path.splitext( c )
    print( configFileNu )
    configNu = Config( configFileNu[0] + 'nu' + configFileNu[1] )
    
    config2 = Config( c )
    config2.readoutNoiseWidth = 0.1

    config3 = Config( c )
    config3.readoutNoiseWidth = 1.

    config4 = Config( c )
    config4.readoutNoiseWidth = 2.5
        
    print( 'rebinningMask', config.rebinningMask )

    adu_keV = 1904.
    
    def generatePlotsExpRebin( fname, settings, preplots = [], plots = [] ):
        res = []
        for exp in config.exposureTimes:
            res += [[]]
            for plot in plots:
                for i, rebin in enumerate(config.rebinningMask):
                    res[-1] += [ lambda ax, exp=exp, rebin=rebin, color='C%s'%i, plot=plot: plot(ax=ax, exp=exp, rebin=rebin, color=color)  ]

        generatePlots( fname,
            settings = settings,
            plots = res
            )
        print( res )
            
    bins = np.arange( min(config.spectrumLims), max(config.spectrumLims), config.spectrumBinsSize )
    
    
    spectrumHist_darkCurrentNoise = lambda config, exposure, rebin, bins: makeHitSumm( config=config, exposure=exposure, darkCurrent=True, noise=True, rebin=rebin ).spectrumHist( bins, normed = True)
    
    spectrumHist_darkCurrentNoiseEmpty = lambda config, exposure, rebin, bins: makeHitSumm( config=config, exposure=exposure, darkCurrent=True, noise=True, rebin=rebin, empty=True ).spectrumHist( bins, normed = True)
    
    absDiffSpectrumHist_darkCurrentNoise = lambda config1, config2, exposure, rebin, bins: np.abs( spectrumHist_darkCurrentNoise( config1, exposure, rebin, bins )
                                - spectrumHist_darkCurrentNoise( config2, exposure, rebin, bins ) )
    if 0:
        generatePlots( 'spectrumDCHitsSub.pdf',
                  settings = [ ylog, xlabel('keV'), ylabel('events/kg/d/keV'), legend, grid ],
                  plots = [
                      [lambda ax, rebin=rebin, bins=bins:
                            ax.step( 
                            bins[:-1], 
                            absDiffSpectrumHist_darkCurrentNoise( config, configNu, 1, rebin, bins ), label='%s'%rebin[0], where='post' )
                            for rebin in config.rebinningMask
                      ],
                      [lambda ax, rebin=rebin, bins=bins: 
                            ax.step( 
                            bins[:-1], 
                            absDiffSpectrumHist_darkCurrentNoise( config, configNu, 2, rebin, bins ), 
                            label='%s'%rebin[0], where='post' )
                            for rebin in config.rebinningMask
                      ],
                      [lambda ax, rebin=rebin, bins=bins: 
                            ax.step( 
                            bins[:-1], 
                            absDiffSpectrumHist_darkCurrentNoise( config, configNu, 3, rebin, bins ),
                            label='%s'%rebin[0], where='post' )
                            for rebin in config.rebinningMask
                      ],
                    ],
                  )
    def signal_noise_rate( rate_bg, rate_nu, T_on, T_off, M, Ebin ):
        signal = rate_nu * T_on * M * Ebin
        noise_on = np.sqrt( (rate_nu + rate_bg) * T_on * M * Ebin )
        noise_off = np.sqrt( rate_bg * T_off * M * Ebin )
        noise = np.sqrt( noise_on**2 + noise_off**2*(T_on/T_off)**2 )
        return signal/noise

    def signal_noise_rate2( rate_off, rate_on, T_on, T_off, M, Ebin ):
        signal = rate_on * T_on * M * Ebin - rate_off1 * T_on * M * Ebin
        noise_on = np.sqrt( rate_on * T_on * M * Ebin )
        noise_off = np.sqrt( rate_off * T_off * M * Ebin )
        noise = np.sqrt( noise_on**2 + noise_off**2*(T_on/T_off)**2 )
        return signal/noise

    def signal_noise_rate3( rate_bg1, rate_bg2, rate_nu, T_on, T_off, M, Ebin ):
        signal = (rate_nu + np.abs(rate_bg1 - rate_bg2) ) * T_on * M * Ebin
        noise_on = np.sqrt( (rate_nu+rate_bg1) * T_on * M * Ebin )
        noise_off = np.sqrt( rate_bg2 * T_off * M * Ebin )
        noise = np.sqrt( noise_on**2 + noise_off**2*(T_on/T_off)**2 )
        return signal/noise
    
    neutrinoRate = readTable('simulatedFITS/neutrinoRate.dat').T
    neutrinoFunction = lambda x: np.interp( x, neutrinoRate[0], neutrinoRate[1] )
    
    spectrumHist_darkCurrentNoiseNeutrino = lambda config, exposure, rebin, bins: makeHitSumm( config=config, exposure=exposure, darkCurrent=True, noise=True, rebin=rebin, neutrinos=True ).spectrumHist( bins, normed = True )

    T_on = 1.
    T_off = 1.
    m = config.CCDmass
    m = 1.
    Ebin = 1.
    make_signal_noise_rate = lambda exposure, rebin, bins: signal_noise_rate( spectrumHist_darkCurrentNoise(config, exposure, rebin, bins), neutrinoFunction(bins[:-1]), T_on, T_off, m, Ebin )
    
    make_signal_noise_rate2 = lambda exposure, rebin, bins:\
        signal_noise_rate2( spectrumHist_darkCurrentNoise(config, exposure, rebin, bins), 
                            spectrumHist_darkCurrentNoiseNeutrino(configNu, exposure, rebin, bins), 
                            T_on, T_off, m, Ebin )
    if 0:
        plots = [
            lambda ax, exp, rebin, color, bins=bins:
                ax.step( bins[:-1], make_signal_noise_rate( exp, rebin, bins ), label='%sh %sx%s'%(exp,rebin[0],rebin[1]), where='post', color=color )
        ]
        generatePlotsExpRebin( 'spectrumDCHitsSNR.pdf',
                  settings = [ ylog, xlabel('keV'), ylabel('SNR'), legend, grid, minorgrid ],
                  plots=plots
                  )

    
    if 0:
        plots = [
                lambda ax, exp, rebin, color, bins=bins: 
                    ax.step( bins[:-1], spectrumHist_darkCurrentNoiseEmpty( config, exp, rebin, bins ), where='post', color=color, lw=1 ),
                lambda ax, exp, rebin, color, bins=bins:
                    ax.step( bins[:-1], spectrumHist_darkCurrentNoise( config, exp, rebin, bins ), label='%sh %sx%s'%(exp, rebin[0], rebin[1]), where='post', color=color )
                ]
        generatePlotsExpRebin( 'spectrumFalse.pdf', settings = [ ylog, xlabel('keV'), ylabel('events/kg/day/keV'), legend, grid, minorgrid ], plots=plots )
        
        
    divideIndexes = lambda a, i, j: a[i]/a[j]
    efficiency_darkCurrentNoise = lambda config, exposure, rebin, bins: divideIndexes( makeHitSumm( config=config, exposure=exposure, darkCurrent=True, noise=True, neutrinos=True, rebin=rebin ).efficiency( bins ), 0, 2 )
    
    make_signal_noise_rate_efficiencyCorrected = lambda config, exposure, rebin, bins: \
        signal_noise_rate( 
            spectrumHist_darkCurrentNoise(config, exposure, rebin, bins), 
            neutrinoFunction(bins[:-1])*efficiency_darkCurrentNoise(config, exposure, rebin, bins), 
            T_on, T_off, m, Ebin )

    if 0:
        plots = [
            lambda ax, exp, rebin, color, bins=bins:
                ax.step( bins[:-1], efficiency_darkCurrentNoise( config, exp, rebin, bins ), label='%sh %sx%s'%(exp,rebin[0],rebin[1]), where='post' )
        ]
    
        generatePlotsExpRebin( 'spectrumDCHitsEfficiency.pdf',
                  settings = [ xlabel('keV'), ylabel('efficiency'), legend, grid ],
                  plots=plots
                  )

    if 0:
        plots = [
            lambda ax, exp, rebin, color, bins=bins:
                ax.step( bins[:-1], make_signal_noise_rate_efficiencyCorrected( config, exp, rebin, bins ), label='%sh %sx%s'%(exp,rebin[0],rebin[1]), where='post' )
        ]
    
        generatePlotsExpRebin( 'spectrumDCHitsSNR_EC.pdf',
                  settings = [ xlabel('keV'), ylabel('SNR'), legend, grid ],
                  plots=plots
                  )

    if 0:
        ref = make_signal_noise_rate_efficiencyCorrected(config, 3, [1,1], bins )
        plots = [
            lambda ax, exp, rebin, color, bins=bins:
                ax.step( bins[:-1], make_signal_noise_rate_efficiencyCorrected( config, exp, rebin, bins )/ref, label='%sh %sx%s'%(exp,rebin[0],rebin[1]), where='post' )
        ]
    
        generatePlotsExpRebin( 'spectrumDCHitsSNR_ECNorm.pdf',
                  settings = [ xlabel('keV'), ylabel('SNR/SNR3'), legend, grid ],
                  plots=plots
                  )

    make_signal_noise_rate_efficiencyCorrected2 = lambda exposure, rebin, bins:\
        signal_noise_rate2( spectrumHist_darkCurrentNoise( config, exposure, rebin, bins),
                            spectrumHist_darkCurrentNoiseNeutrino(exposure, rebin, bins),
                            T_on, T_off, m, Ebin )

    make_signal_noise_rate_efficiencyCorrected2 = lambda exposure, rebin, bins:\
        signal_noise_rate2( spectrumHist_darkCurrentNoise( config, exposure, rebin, bins),
                            spectrumHist_darkCurrentNoiseNeutrino(exposure, rebin, bins),
                            T_on, T_off, m, Ebin )


    if 0:
        plots = [
            lambda ax, exp, rebin, color, bins=bins:
                ax.step( bins[:-1], make_signal_noise_rate_efficiencyCorrected2( exp, rebin, bins ), label='%sh %sx%s'%(exp,rebin[0],rebin[1]), where='post' )
        ]
    
        generatePlotsExpRebin( 'spectrumDCHitsSNR_EC2.pdf',
                  settings = [ ylog, xlabel('keV'), ylabel('SNR'), legend, grid ],
                  plots=plots
                  )

    make_signal_noise_rate_efficiencyCorrected3 = lambda exposure, rebin, bins:\
        signal_noise_rate3( spectrumHist_darkCurrentNoise( config, exposure, rebin, bins),
                            spectrumHist_darkCurrentNoise( configNu, exposure, rebin, bins),
                            neutrinoFunction(bins[:-1])*efficiency_darkCurrentNoise( config, exposure, rebin, bins ),
                            T_on, T_off, m, Ebin )
    if 0:
        plots = [
            lambda ax, exp, rebin, color, bins=bins:
                ax.step( bins[:-1], make_signal_noise_rate_efficiencyCorrected3( exp, rebin, bins ), label='%sh %sx%s'%(exp,rebin[0],rebin[1]), where='post' )
        ]
        generatePlotsExpRebin( 'spectrumDCHitsSNR_EC3.pdf',
                  settings = [ ylog, xlabel('keV'), ylabel('SNR'), legend, grid ],
                  plots=plots
                  )
    return

class SpectrumWindow(Gtk.Window):

    def __del__(self):
        os.remove(self.tempPath)
        Gtk.Window.__del__(self)
        
    def __init__(self):
        self.id = time.time()
        print 'id', self.id
        self.tempPath = 'gui/session%s.png'%self.id
        Gtk.Window.__init__( self, title="Spectrum Viewer [session%s]"%self.id )
        self.set_border_width(3)

        self.paths = []
        self.pathslabel = Gtk.Label()
        self.pathslabel.set_label( 'None' )

        drawspectrumbutton = Gtk.Button(use_underline=True)
        drawspectrumbutton.set_label('_spectrum')
        drawspectrumbutton.connect( 'clicked', self.onspectrumclick )
        
        runsLabel = Gtk.Label()
        runsLabel.set_label('runs')
        
        listRuns = Gtk.VBox()
        for run in listAll_runs()[::-1]:
            runbutton = Gtk.ToggleButton()
            runbutton.set_label( run )
            runbutton.connect( 'toggled', self.ontogglebuttonruns, run )
            listRuns.pack_start( runbutton, False, False, 0 )

        scrolledwindow = Gtk.ScrolledWindow()
        scrolledwindow.add_with_viewport( listRuns )
        
        runIDPanel = Gtk.VBox()
        runIDPanel.pack_start( runsLabel, False, False, 1 )
        runIDPanel.pack_start( scrolledwindow, True, True, 0 )
        
        selectionLabel = Gtk.Label()
        selectionLabel.set_label( 'selection' )
        self.selectionEntry = Gtk.Entry()
        self.selectionEntry.set_text( 'ohdu==2 && E0<50e3' )
        
        expressionLabel = Gtk.Label()
        expressionLabel.set_label( 'expr' )
        self.expressionEntry = Gtk.Entry()
        self.expressionEntry.set_text( 'E0' )
        
        binSizeLabel = Gtk.Label()
        binSizeLabel.set_label( 'binSize' )
        self.binSizeEntry = Gtk.Entry()
        self.binSizeEntry.set_text( 'auto' )
        
        fitFunctionLabel = Gtk.Label()
        fitFunctionLabel.set_label( 'fit function' )
        self.fitFunctionEntry = Gtk.Entry()
        self.fitFunctionEntry.set_text( 'lambda x, a, b: np.log( a/x**3 + b/x + 1 )' )
        
        self.image = Gtk.Image()
        
        firstLine = Gtk.HBox()
        firstLine.pack_start(selectionLabel, False, False, 1)
        firstLine.pack_start(self.selectionEntry, True, True, 1)
        
        secondLine = Gtk.HBox()
        secondLine.pack_start(self.selectionEntry, True, True, 1)
        secondLine.pack_start(expressionLabel, False, False, 1)
        secondLine.pack_start(self.expressionEntry, True, True, 1)
        secondLine.pack_start(binSizeLabel, False, False, 1)
        secondLine.pack_start(self.binSizeEntry, False, False, 1)
        
        thirdLine = Gtk.HBox()
        thirdLine.pack_start( fitFunctionLabel, False, False, 1 )
        thirdLine.pack_start( self.fitFunctionEntry, True, True, 1 )
        
        optionsPanel = Gtk.VBox()
        optionsPanel.pack_start(firstLine, True, True, 1)
        optionsPanel.pack_start(secondLine, True, True, 1)
        optionsPanel.pack_start(thirdLine, True, True, 1)

        imagePanel = Gtk.VBox()
        imagePanel.pack_start( optionsPanel, False, False, 1 )
        imagePanel.pack_start( self.image, True, True, 1 )
        
        body = Gtk.HBox()
        body.pack_start( runIDPanel, False, False, 0 )
        body.pack_end( imagePanel, True, True, 1 )

        footer = Gtk.HBox()
        footer.pack_end( drawspectrumbutton, False, False, 1 )
        footer.pack_start( self.pathslabel, True, True, 1 )
        
        document = Gtk.VBox()
        document.pack_start( body, True, True, 3 )
        document.pack_start( footer, False, False, 3 )
        
        self.add( document )
        self.maximize()
        self.connect( 'destroy', Gtk.main_quit )
        #self.connect( 'key_press_event', self.onkeypressevent )
        self.show_all()

    def ontogglebuttonruns( self, button, run ):
        firstRunID = listFITS_run( run )[0]
        catalog = getAssociatedCatalog( firstRunID )[0]
        if button.get_active() == True:
            self.paths.append( catalog )
        else:
            self.paths.remove( catalog )
        self.pathslabel.set_label( '\n'.join( [ '[%d] %s'%(i,path) for i, path in enumerate(self.paths) ] ) )
        self.show_all()
    
    def onspectrumclick( self, button ):
        fig = plt.figure()
        fig.suptitle( self.selectionEntry.get_text() )
        grid = plt.GridSpec( 4, 1, hspace=0, wspace=0 )
        ax_main = fig.add_subplot( grid[:-1, 0])
        ax_ratio = fig.add_subplot( grid[-1,0])
        
        binMin = None
        binMax = None
        data = []
        nRunIDs = []
        for path in self.paths:
            datum = root_numpy.root2array( path, treename='hitSumm', branches=self.expressionEntry.get_text(), selection=self.selectionEntry.get_text() )
            data.append( datum )
            nRunIDs.append( len( listFITS_run( getrunFromPath(path) ) ) )
            if binMin is None: binMin = np.min(datum)
            else: binMin = min(binMin, np.min( datum ) )
            if binMax is None: binMax = np.max(datum)
            else: binMax = max(binMax, np.max( datum ) )
        if self.binSizeEntry.get_text() == 'auto': 
            bins = np.linspace( binMin, binMax, 100 )
        else: 
            bins = np.arange( binMin, binMax, float(self.binSizeEntry.get_text()) )
            
        for i, (datum, path, nRunID) in enumerate( zip( data, self.paths, nRunIDs ) ):
            hist, _ = np.histogram( datum, bins=bins )
            hist = hist.astype(float)/nRunID
            x = .5*(bins[:-1]+bins[1:])
            median = np.median( hist )
            funcstr = self.fitFunctionEntry.get_text()
            varnames = map( lambda s:s.strip(), funcstr.split(':')[0].split(',')[1:] )
            func = eval( funcstr )
            p = scipy.optimize.curve_fit( func, x[hist>1], np.log( hist[hist>1] ) )[0]
            y = np.exp( func( x, *p ) )
            label = '%s [%d] '%(getrunFromPath(path),nRunID) + ' '.join( [ '%s:%.2g'%(v,pp) for v,pp in zip(varnames, p) ] )
            ax_main.step( x, hist, label=label, where='mid', color='C%d'%i )
            ax_main.plot( x, y, color='C%d'%i, linestyle='--' )
            ax_ratio.step( x, hist/y, where='mid', color='C%d'%i )
        ax_main.set_yscale( 'log' )
        ax_ratio.set_xlabel( self.expressionEntry.get_text() )
        ax_ratio.set_yscale( 'log' )
        ax_main.legend( fancybox=True, framealpha=0 )
        fig.savefig( self.tempPath )
        plt.close()
        self.image.set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.tempPath ) )
        self.show_all()

class ImageWindow(Gtk.Window):

    def __init__(self):
        self.id = time.time()
        print 'id', self.id
        self.tempPath = 'gui/session%s.png'%self.id
        Gtk.Window.__init__( self, title="Image Viewer [session%s]"%self.id )
        self.set_border_width(3)

        self.paths = []
        self.pathslabel = Gtk.Label()
        self.pathslabel.set_label( 'None' )
        self.drawimagebutton = Gtk.Button()
        self.drawimagebutton.set_label('image')
        self.drawimagebutton.connect( 'clicked', self.plotImage )
        
        footer = Gtk.HBox()
        footer.pack_end( self.drawimagebutton, False, False, 1 )
        footer.pack_start( self.pathslabel, True, True, 1 )
        
        self.runID = None
        self.runIDcolumns = {}
        self.runButtons = []
        self.runIDButtons = []
        
        runsLabel = Gtk.Label()
        runsLabel.set_label('runs')
        
        runlabel = Gtk.Label()
        runlabel.set_label( 'runs' )

        runcolumn = Gtk.VBox()
        runcolumn.pack_start( runlabel, False, True, 1 )
        for run in listAll_runs()[::-1]:
            self.runButtons.append( Gtk.ToggleButton() )
            self.runButtons[-1].set_label( run )
            self.runButtons[-1].connect( 'toggled', self.ontogglebuttonruns, run )
            runcolumn.pack_start( self.runButtons[-1], False, False, 0 )

        scrolledwindow = Gtk.ScrolledWindow()
        scrolledwindow.add_with_viewport( runcolumn )
        
        self.selectionLabel = Gtk.Label()
        self.selectionLabel.set_label( 'selection' )
        self.selectionEntry = Gtk.Entry()
        self.selectionEntry.set_text( '' )

        self.radio = {}
        firstLine = Gtk.HBox()
        for key in ['scn','mbs','osi','raw']:
            self.radio[key] = Gtk.ToggleButton()
            self.radio[key].set_label( key )
            self.radio[key].connect( 'toggled', self.onradiotoggled )
            firstLine.pack_end( self.radio[key], False, False, 1 )
        self.radio['scn'].set_active(True)
        
        ohduLabel = Gtk.Label()
        ohduLabel.set_label( 'ohdu' )
        self.ohduEntry = Gtk.Entry()
        self.ohduEntry.set_text('2')
        eWidthLabel = Gtk.Label()
        eWidthLabel.set_label( 'eWidth' )
        self.eWidthEntry = Gtk.Entry()
        self.eWidthEntry.set_text('10')

        xLabel = Gtk.Label()
        xLabel.set_label( 'x' )
        self.xMinEntry = Gtk.Entry()
        self.xMinEntry.set_text('auto')
        self.xMaxEntry = Gtk.Entry()
        self.xMaxEntry.set_text('auto')

        yLabel = Gtk.Label()
        yLabel.set_label( 'x' )
        self.yMinEntry = Gtk.Entry()
        self.yMinEntry.set_text('auto')
        self.yMaxEntry = Gtk.Entry()
        self.yMaxEntry.set_text('auto')
        firstLine.pack_start( ohduLabel, False, False, 1 )
        firstLine.pack_start( self.ohduEntry, False, False, 1 )
        firstLine.pack_start( eWidthLabel, False, False, 1 )
        firstLine.pack_start( self.eWidthEntry, False, False, 1 )
        firstLine.pack_start( xLabel, False, False, 1 )
        firstLine.pack_start( self.xMinEntry, False, False, 1 )
        firstLine.pack_start( self.xMaxEntry, False, False, 1 )
        firstLine.pack_start( yLabel, False, False, 1 )
        firstLine.pack_start( self.yMinEntry, False, False, 1 )
        firstLine.pack_start( self.yMaxEntry, False, False, 1 )
        
        self.imageCanvas = Gtk.Image()
        scrolledwindowImage = Gtk.ScrolledWindow()
        scrolledwindowImage.add_with_viewport( self.imageCanvas )
        
        optionsPanel = Gtk.VBox()
        optionsPanel.pack_start(firstLine, False, True, 1)

        imagePanel = Gtk.VBox()
        imagePanel.pack_start( optionsPanel, False, False, 1 )
        imagePanel.pack_start( scrolledwindowImage, True, True, 1 )

        self.body = Gtk.HBox()
        self.body.pack_start( scrolledwindow, False, False, 1 )
        self.body.pack_end( imagePanel, True, True, 1 )
        
        self.document = Gtk.VBox()
        self.document.pack_start( self.body, True, True, 3 )
        self.document.pack_start( footer, False, False, 3 )
        
        self.add( self.document )
        self.maximize()
        self.connect( 'destroy', Gtk.main_quit )
        self.show_all()
        
    def onradiotoggled( self, radio ):
        if radio.get_active() == True:
            self.imageType = radio.get_label()
            for key, button in self.radio.items():
                if radio != button:
                    button.set_active(False)
            if self.runID is not None:
                self.pathslabel.set_label( '\n'.join( self.getPath( self.runID ) ) )
        self.show_all()
        
    def showListRunIDs( self, run ):
        self.runIDcolumns[ run ] = Gtk.VBox()
        self.body.pack_start( self.runIDcolumns[ run ], False, False, 1 )
        scrolledwindow = Gtk.ScrolledWindow()
        runlabel = Gtk.Label()
        runlabel.set_label( run )
        self.runIDcolumns[ run ].pack_start( runlabel, False, True, 1 )
        self.runIDcolumns[ run ].pack_start( scrolledwindow, True, True, 1 )
        
        subbox = Gtk.VBox()
        scrolledwindow.add_with_viewport( subbox )
        
        self.runIDButtons = []
        for path in listFITS_run( run ):
            runID = getrunIDFromPath(path)
            self.runIDButtons.append( Gtk.ToggleButton() )
            self.runIDButtons[-1].set_label( runID )
            self.runIDButtons[-1].connect( 'toggled', self.ontogglebuttonrunID, runID )
            subbox.pack_start( self.runIDButtons[-1], False, False, 1 )
    
    def getPath( self, runID ):
        if self.imageType == 'raw':
            return listFITSraw_runID( int(runID) )
        elif self.imageType == 'osi':
            return listFITSosi_runID( int(runID) )
        elif self.imageType == 'scn':
            return listFITSscn_runID( int(runID) )
        return ['no file found']
    
    def ontogglebuttonruns( self, button, run ):
        firstRunID = listFITS_run( run )[0]
        catalog = getAssociatedCatalog( firstRunID )[0]
        if button.get_active() == True:
            self.paths.append( catalog )
            self.showListRunIDs( run )
            for runButton in self.runButtons:
                if runButton != button:
                    runButton.set_active(False)
        else:
            self.paths.remove( catalog )
            self.runIDcolumns[run].destroy()
            self.pathslabel.set_label( 'not file selected' )
        self.runID = None
        self.pathslabel.set_label( 'not file selected' )

        self.show_all()
        
    def ontogglebuttonrunID( self, button, runID ):
        if button.get_active() == True:
            self.runID = runID
            for runIDButton in self.runIDButtons:
                if runIDButton != button:
                    runIDButton.set_active(False)
            self.pathslabel.set_label( '\n'.join( self.getPath( self.runID ) ) )
        self.show_all()

    def plotImage( self, button ):
        ohdu = self.ohduEntry.get_text()
        image = Image( runID = self.runID, ohdu = int(ohdu), imageType = self.imageType )
        
        text = lambda s: None if s.get_text() == 'auto' else int(s.get_text())
        mask = np.s_[ text(self.xMinEntry):text(self.xMaxEntry), text(self.yMinEntry):text(self.yMaxEntry) ]
        image = image(mask)
        
        if self.xMaxEntry.get_text() == 'auto':
            self.xMaxEntry.set_text( str( image.shape[0] ) )
        if self.yMaxEntry.get_text() == 'auto':
            self.yMaxEntry.set_text( str( image.shape[1] ) )
        
        w, h = plt.figaspect(1.)
        fig = plt.figure(figsize=(2*w, 2*h))
        fig.suptitle( 'ohdu=%s'%ohdu )
        yticks = [ 0, 1055, 2*1055, 3*1055, 4*1055 ]
        xticks = [ 0, image.shape[1] - image.overscanWidth, image.shape[1] ]
        grid = plt.GridSpec(4, 4, hspace=0, wspace=0)
        
        main_ax = fig.add_subplot(grid[1:, :-1])
        y_hist = fig.add_subplot(grid[1:, -1], yticks=yticks, yticklabels=[], sharey=main_ax)
        x_hist = fig.add_subplot(grid[0, :-1], xticks=xticks, xticklabels=[], sharex=main_ax)
        zoom_ax = fig.add_subplot(grid[0, -1], yticklabels=[], sharex=y_hist)
        zoom_ax.yaxis.tick_right()
        zoom_ax.xaxis.tick_top()
        
        main_ax.tick_params(direction='in')
        
        median = image.median()
        mad = image.mad()

        mid = lambda _: .5*(_[1:]+_[:-1])
        diff = lambda _: _[1:]-_[:-1]
        
        eWidth = float(self.eWidthEntry.get_text())
        eMin = median - eWidth*mad
        eMax = median + eWidth*mad
        ebins = np.arange( eMin, eMax, 1 )

        main_ax.imshow( image.image, vmin = eMin, vmax = eMax )
        x_hist.hist2d( np.repeat(range(image.shape[1]), image.shape[0] ), image.T.flatten(), bins=[image.shape[1], ebins ], norm = colors.LogNorm() )
        x_hist.plot( range(image.shape[1]), image.median( axis = 0 ) )

        y_hist.hist2d( image.flatten(), np.repeat(range(image.shape[0]), image.shape[1] ) , bins=[ebins, image.shape[0]], norm = colors.LogNorm() )
        y_hist.plot( image.median( axis = 1 ), range(image.shape[0]) )

        hist, bins, _ = zoom_ax.hist( image.flatten(), bins=ebins, histtype='step', color='b' )
        if image.overscanWidth != 0 and self.imageType == 'scn':
            zoom_ax.hist( image.os().flatten(), bins=ebins, histtype='step', color='r' )
            
        zoom_ax.set_yscale('log')
        a, b = zoom_ax.get_ylim()
        zoom_ax.set_ylim((.1, b))
        zoom_ax.grid( True )
        main_ax.invert_yaxis()
        main_ax.set_xlim( (0,image.shape[1]) )
        main_ax.set_ylim( (0,image.shape[0]) )

        fig.savefig( self.tempPath )
        plt.close()
        self.imageCanvas.set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.tempPath ) )
        self.show_all()
        
class Callable:
    @staticmethod
    def spectrumViewer():
        win = SpectrumWindow()
        Gtk.main()

    @staticmethod
    def imageViewer():
        win = ImageWindow()
        Gtk.main()

    @staticmethod
    def plotSpectrumSimulated( c ):
        config = Config(c)
        numberOfRunIDs = 10
        numberOfOHDUs = 10
        
        baseLines = { runID: { ohdu: random.randint(1200,2000) for ohdu in range(numberOfOHDUs) } for runID in range(numberOfRunIDs) }
        
        images = { runID: { ohdu: 
                               simulateRawImage( 
                                   baseLineL = random.randint(1200,2000),
                                   baseLineR = random.randint(1200,2000), 
                                   gainL = random.randint(1200,2000),
                                   gainR = random.randint(1200,2000), 
                                   aL = random.randint(100,120),
                                   aR = random.randint(100,120) ) 
                               for ohdu in range(numberOfOHDUs) } for runID in range(numberOfRunIDs) }

    @staticmethod
    def computeMADMatrix():
        runIDs031 = range(3693, 3693+60) #range(3693, 4307)
        runIDs043 = range(6322, 6322+180) #range(6322, 6409)
        runIDs = range(6381, 6385)
        ohdus = [2,3,4,5,6,7,8,9,10,13,14]

        thr = 15
        eMin = thr
        eMax = 2000*10
        binWidth = 250
        
        #processImage4( hitSumm._getSCNImage_adu( 3693, 3 ), 15 )
        #exit(0)
        #images = RunIDImages( runIDs, ohdus )
        #images.generateCatalog( thr=thr, eMin=eMin, eMax=eMax, binWidth=binWidth )
        #exit(0)
        
        images031 = RunIDImages( runIDs031, ohdus, imageType='root' )
        RunIDImages.printSpectrum( images031, eMin=eMin, eMax=eMax, binWidth=binWidth, fname=lambda ohdu: 'spectrum_root_031_%s.pdf'%ohdu )

        images031 = RunIDImages( runIDs031, ohdus, imageType='scn' )
        RunIDImages.printSpectrum( images031, eMin=eMin, eMax=eMax, binWidth=binWidth )
        exit(0)
        
        images043 = RunIDImages( runIDs043, ohdus, imageType='root' )
        RunIDImages.printSpectrum( images043, eMin=eMin, eMax=eMax, binWidth=binWidth )

        RunIDImages.printSpectrum( [images031, images043], eMin=eMin, eMax=eMax, binWidth=binWidth, fname=lambda ohdu: 'spectrum4_031_043_ohdu%s.pdf'%ohdu )

        
    @staticmethod
    def plotBaselineGainNoise():
        ohdus = [2,3,7,10]
        gains = { ohdu:[] for ohdu in ohdus }
        medianL = { ohdu:[] for ohdu in ohdus }
        medianR = { ohdu:[] for ohdu in ohdus }
        meanL = { ohdu:[] for ohdu in ohdus }
        meanR = { ohdu:[] for ohdu in ohdus }
        stdL = { ohdu:[] for ohdu in ohdus }
        stdR = { ohdu:[] for ohdu in ohdus }
        madL = { ohdu:[] for ohdu in ohdus }
        madR = { ohdu:[] for ohdu in ohdus }
        madLa = { ohdu:[] for ohdu in ohdus }
        medianLa = { ohdu:[] for ohdu in ohdus }
        
        x = []
        for runID in range(3615, 4200):
            for ohdu in ohdus:
                gain = hitSumm._readGain(runID,ohdu)
                if gain is None: continue
                gains[ohdu].append( gain )
                res = hitSumm._readBaselineNoise(runID, ohdu)
                meanL[ohdu].append(res[0][0])
                meanR[ohdu].append(res[1][0])
                stdL[ohdu].append(res[0][1])
                stdR[ohdu].append(res[1][1])
                medianL[ohdu].append(res[0][2])
                medianR[ohdu].append(res[1][2])
                madL[ohdu].append(res[0][3])
                madR[ohdu].append(res[1][3])
                
                madLa[ohdu].append(res[2][3])
                medianLa[ohdu].append(res[2][2])
                
            if gain is None: continue
            x.append(runID)
            plots11 = [ lambda ax: ax.scatter(x, np.array(medianL[2]), marker='_', label='medianL'),
                        lambda ax: ax.scatter(x, np.array(medianR[2]), marker='_', label='medianR'),
                        #lambda ax: ax.scatter(x, np.array(meanL[2]), marker='|', label='meanL'),
                        #lambda ax: ax.scatter(x, np.array(meanR[2]), marker='|', label='meanR'), 
                        lambda ax: ax.scatter(x, np.array(medianLa[2]), marker='x', label='medianLa'),
                        lambda ax: ax.scatter(x, np.array(madLa[2]), marker='|', label='madLa'),
                        lambda ax: ax.scatter(x, np.array(gains[2]), marker='o', label='gain'), ]
            
            plots21 = [ lambda ax: ax.scatter(x, np.array(stdL[2]), marker='_', label='stdL'),
                        lambda ax: ax.scatter(x, np.array(stdR[2]), marker='_', label='stdR'),
                        lambda ax: ax.scatter(x, np.array(madL[2]), marker='|', label='madL'),
                        lambda ax: ax.scatter(x, np.array(madR[2]), marker='|', label='madR'),]
                        #lambda ax: ax.scatter(x, np.array(madLa[2]), marker='x', label='madLa'), ]

            plots12 = [ lambda ax: ax.scatter(x, np.array(medianL[7]), marker='_', label='medianL'),
                        lambda ax: ax.scatter(x, np.array(medianR[7]), marker='_', label='medianR'),
                        #lambda ax: ax.scatter(x, np.array(meanL[7]), marker='|', label='meanL'),
                        #lambda ax: ax.scatter(x, np.array(meanR[7]), marker='|', label='meanR'),
                        lambda ax: ax.scatter(x, np.array(medianLa[7]), marker='x', label='medianLa'),
                        lambda ax: ax.scatter(x, np.array(madLa[7]), marker='|', label='madLa'),
                        lambda ax: ax.scatter(x, np.array(gains[7]), marker='o', label='gain'), ]
            
            plots22 = [ lambda ax: ax.scatter(x, np.array(stdL[7]), marker='_', label='stdL'),
                        lambda ax: ax.scatter(x, np.array(stdR[7]), marker='_', label='stdR'),
                        lambda ax: ax.scatter(x, np.array(madL[7]), marker='|', label='madL'),
                        lambda ax: ax.scatter(x, np.array(madR[7]), marker='|', label='madR'),]
                        #lambda ax: ax.scatter(x, np.array(madLa[7]), marker='x', label='madLa'), ]

            plots13 = [ lambda ax: ax.scatter(x, np.array(medianL[10]), marker='_', label='medianL'),
                        lambda ax: ax.scatter(x, np.array(medianR[10]), marker='_', label='medianR'),
                        #lambda ax: ax.scatter(x, np.array(meanL[10]), marker='|', label='meanL'),
                        #lambda ax: ax.scatter(x, np.array(meanR[10]), marker='|', label='meanR'),
                        lambda ax: ax.scatter(x, np.array(medianLa[10]), marker='x', label='medianLa'),
                        lambda ax: ax.scatter(x, np.array(madLa[10]), marker='|', label='madLa'),
                        lambda ax: ax.scatter(x, np.array(gains[10]), marker='o', label='gain'), ]

            
            plots23 = [ lambda ax: ax.scatter(x, np.array(stdL[10]), marker='_', label='stdL'),
                        lambda ax: ax.scatter(x, np.array(stdR[10]), marker='_', label='stdR'),
                        lambda ax: ax.scatter(x, np.array(madL[10]), marker='|', label='madL'),
                        lambda ax: ax.scatter(x, np.array(madR[10]), marker='|', label='madR'),]

            plots14 = [ lambda ax: ax.scatter(x, np.array(medianL[3]), marker='_', label='medianL'),
                        lambda ax: ax.scatter(x, np.array(medianR[3]), marker='_', label='medianR'),
                        #lambda ax: ax.scatter(x, np.array(meanL[3]), marker='|', label='meanL'),
                        #lambda ax: ax.scatter(x, np.array(meanR[3]), marker='|', label='meanR'),
                        lambda ax: ax.scatter(x, np.array(medianLa[3]), marker='x', label='medianLa'),
                        lambda ax: ax.scatter(x, np.array(madLa[3]), marker='|', label='madLa'),
                        lambda ax: ax.scatter(x, np.array(gains[3]), marker='o', label='gain'), ]

            
            plots24 = [ lambda ax: ax.scatter(x, np.array(stdL[3]), marker='_', label='stdL'),
                        lambda ax: ax.scatter(x, np.array(stdR[3]), marker='_', label='stdR'),
                        lambda ax: ax.scatter(x, np.array(madL[3]), marker='|', label='madL'),
                        lambda ax: ax.scatter(x, np.array(madR[3]), marker='|', label='madR'),]

            generatePlots( 'gainCorrelation.pdf',
                settings = [ xlabel('runIDs'), legend, grid ],
                plots=[ [ plots11, plots12, plots13, plots14], 
                        [ plots21, plots22, plots23, plots24] ]
                )

        exit(0)
            
    @staticmethod
    def plotSpectrumRaw( c ):
        config = Config(c)
        #logbins = np.logspace(np.log10(1e-3), np.log10(100), 100)
        linbins = np.linspace(0*7.13, 300*7.13, 50)
        linbins = np.arange( 0*7.13, 300*7.13, .1 )
        linbins = np.arange( 0*7.13, 300*7.13, .1 )
        #ax.set_xticklabels([str(round(float(label), 2)) for label in labels])
        bins = linbins
        #axlim = lambda ax: ax.set_xlim([0,5*7.13])
        axlim = lambda ax: None
        
        spectra = makeHitSumm( config=config, exposure=3, darkCurrent=True, copper=True, neutrinos=True, empty=True, mbs = False ).spectrumDiff(bins)
        
        plots = [
            [   
                lambda ax,i=i,spectrum=spectrum: ax.step(
                    *spectrum, 
                    where='post', color='C%d'%i, label=0 ) for i, spectrum in enumerate(spectra)
                #lambda ax: ax.step( 
                    #*makeHitSumm( config=config2, exposure=3, darkCurrent=True, copper=True, neutrinos=True, noise=True, empty=True, mbs = False ).spectrumDiff(bins)[0], 
                    #where='post', color='g', label=config2.readoutNoiseWidth ),
                #lambda ax: ax.step(
                    #*makeHitSumm( config=config3, exposure=3, darkCurrent=True, copper=True, neutrinos=True, noise=True, empty=True, mbs = False ).spectrumDiff(bins)[0], 
                    #where='post', color='c', label=config3.readoutNoiseWidth ),
                #lambda ax: ax.step( 
                    #*makeHitSumm( config=config, exposure=3, darkCurrent=True, copper=True, neutrinos=True, noise=True, empty=True, mbs = False ).spectrumDiff(bins)[0], 
                    #where='post', color='r', label=config.readoutNoiseWidth ),
                #lambda ax: ax.step( 
                    #*makeHitSumm( config=config4, exposure=3, darkCurrent=True, copper=True, neutrinos=True, noise=True, empty=True, mbs = False ).spectrumDiff(bins)[0], 
                    #where='post', color='m', label=config4.readoutNoiseWidth ),
                ]
            ]
        generatePlots( 'spectrumDiffHitsNu.pdf',
                settings = [ xlabel('diff'), legend, grid, axlim ],
                plots=plots
                )
        
        exit(0)

class CCDImage:
    def getCorrelatedNoise( self, freq ):
        L = image.shape[0]*image.shape[1]
        x = np.arange(0,L,1.)
        return np.cos( x/L/freq ).reshape( image.shape[0], image.shape[1] )
    
    def getReadingDeviation( self, a ):
        x = np.arange( 0, 1000, 1. )
        y = (x/a/1000.)**2[::-1]
        image = np.zeros( image.shape[0], image[1] )
        image[ :, 0:1000 ] -= y
        return image
    
    def composeImage( self ):
        return None

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
def readTable( fname ):
    table = open( fname, 'r' ).readlines()
    data = []
    for line in table:
        if line[0] == '#': continue
        data += [ map( float, line.split() ) ] 
    return np.array( data )

def makeHitSumm( **kwargs ):
    obj = hitSumm( **kwargs )
    return obj

class hitSumm:
    def __init__(self, rebin = None, exposure = None, neutrinos = False, darkCurrent = False, noise = False, empty = False, extract = False, config = None, mbs = True, copper = False ):
        if rebin is None:
            self.rebin = (1,1)
        else:
            self.rebin = rebin
        if self.rebin is not None:
            self.rebin = tuple(self.rebin)
        self.exposure = exposure
        self.neutrinos = neutrinos
        self.darkCurrent = darkCurrent
        self.noise = noise
        self.empty = empty
        self.config = config
        self.extract = extract
        self.mbs = mbs
        self.copper = copper
        
        self.extract_id = 'xthr%s_'%self.config.extractThreshold
        self.scn_id = 'scn_'
        self.raw_id = 'raw_'
        self.mbs_id = 'mbs_'
        self.neutrino_id = lambda n: 'nu%s_%s_%s'%( self.config.neutrinosEnergyRange[0], self.config.neutrinosEnergyRange[1], self.config.numberNeutrinos )
        self.copper_id = lambda n: 'Cu%s_'%n
        self.copper2_id = lambda n: 'Cu2%s_'%n
        self.silicon_id = lambda n: 'Si%s_'%n
        self.darkCurrent_id = 'dc%s_'%(self.config.darkCurrent*self.exposure)
        if self.rebin:
            self.rebin_id = 'rb%dx%d_'%self.rebin
        self.noise_id = 'rn%s_'%self.config.readoutNoiseWidth
        if self.exposure:
            self.exposure_id = 'exp%.1fh_'%self.exposure
        
        self.hits_ext = '.root'
        self.table_ext = '.dat'
        self.image_ext = '.fits'
        
        self.partial_id = config.outPath

        def makeFull_id( runID ):
            full_id = config.outPath
            full_id += self.extract_id
            if self.darkCurrent:
                full_id += self.darkCurrent_id
            if self.rebin is not None:
                full_id += self.rebin_id
            if self.noise:
                full_id += self.noise_id
            if self.neutrinos:
                full_id += self.neutrino_id
            if self.copper:
                full_id += self.copper_id
            if self.mbs:
                full_id += self.mbs_id
            if not self.empty:
                full_id += self.scn_id + str(runID)
            if self.exposure is not None:
                full_id += self.exposure_id
            return full_id.strip('_')
        
        self.full_id = lambda runID: makeFull_id( runID )
        self.readGain( self.config.runID, self.config.ohdu )
        self.nu_dtype = [('xPix', np.float), ('yPix', np.float), ('z_um', np.float), ('ePix_keV', np.float), ('n_e', np.int)]

        self.smoothWidth = 20
        print( 'hitSumm', self.full_id('runID') )
        
        self.adu_keV = self.readGain(self.config.runID, self.config.ohdu)
        if self.adu_keV is None:
            self.getSCNHits_adu(self.config.runIDmin)
        
    def getNeutrinosTable( self, save = False ):
        
        self.partial_id += self.neutrino_id
        path = self.config.outPath + self.neutrino_id( self.config.numberNeutrinos ) + self.table_ext
        if os.path.exists( path ):
            res = readTable( path )
            self.neutrinosTable = res
            return res
        
        xlen_px, ylen_px = self.config.imageShape
        ylen_px -= self.config.overScanCols
        xlen_px -= self.config.overScanLines

        x_px = xlen_px*np.random.rand( self.config.numberNeutrinos )
        y_px = ylen_px*np.random.rand( self.config.numberNeutrinos )
        z_um = self.config.heightCCD*np.random.rand( self.config.numberNeutrinos )
        E_keV = (max( self.config.neutrinosEnergyRange ) - min( self.config.neutrinosEnergyRange ) ) * np.random.rand( self.config.numberNeutrinos ) + min( self.config.neutrinosEnergyRange )
        numberElectrons = np.rint( E_keV/self.config['e-[keV]'] ).astype(int)
        
        table = np.array( zip( x_px, y_px, z_um, E_keV, numberElectrons ) )
        
        if save:
            writeArraysTo(
                arrays = [ x_px, y_px, z_um, E_keV, numberElectrons ], 
                labels = ['x[px]', 'y[px]', 'z[um]', 'E[keV]', 'e-'], 
                fmt = ['%d', '%d', '%.4e', '%.4e', '%d'], 
                fname = path
                )
        return table

    def getCopperImage_e( self, number, save=False, particle='Cu' ):
        print( '+Cu', )
        t = timeit.default_timer()
        path = self.config.outPath + self.copper_id + self.table_ext
        
        xlen_px, ylen_px = self.config.imageShape
        ylen_px -= self.config.overScanCols
        xlen_px -= self.config.overScanLines
        
        if particle == 'Cu':
            #number = self.config.numberCu
            E = self.config.Cu
            self.partial_id += self.copper_id( number )
        elif particle == 'Cu2':
            #number = self.config.numberCu2
            E = self.config.Cu2
            self.partial_id += self.copper2_id( number )
        elif particle == 'Si':
            #number = self.config.numberSi
            E = self.config.Si
            self.partial_id += self.silicon_id( number )
            
        x_px = xlen_px*np.random.rand( number )
        y_px = ylen_px*np.random.rand( number )
        z_um = self.config.heightCCD * np.random.rand( number )
        E_keV = E * np.ones( number )
        nElectrons = np.rint( E_keV/self.config['e-[keV]'] ).astype(int)
        
        ex_px, ey_px, ez_um = np.array( [ [ix_px, iy_px, iz_um] for ix_px, iy_px, iz_um, ne in zip(x_px, y_px, z_um, nElectrons) for _ in range(int(ne)) ] ).T
        iex_px, iey_px = round2int( *diffuseElectrons( ex_px, ey_px, ez_um ) )
        inside_boudary = np.all( [ iex_px >= 0, iex_px < xlen_px, iey_px >= 0 , iey_px < ylen_px ], axis=0 )
        image = np.zeros( self.config.imageShape )
        np.add.at( image, [ iex_px[inside_boudary], iey_px[inside_boudary] ], 1. )
        if save:
            saveImage( image, path )
        print( '%.2fs'%(timeit.default_timer() - t), )
        return image
        
    def getNeutrinosImage_e( self, save = False ):
        self.partial_id += self.neutrino_id( self.config.numberNeutrinos )
        print( '+nu', )
        t = timeit.default_timer()

        path = self.config.outPath + self.neutrino_id + self.image_ext
        if os.path.exists( path ):
            print( '(readnu)', )
            return readImage( path )
        
        xlen_px, ylen_px = self.config.imageShape
        ylen_px -= self.config.overScanCols
        xlen_px -= self.config.overScanLines
        
        res = self.getNeutrinosTable( save=True )
        x_px, y_px, z_um, E_keV, nE = res.T
        ex_px, ey_px, ez_um = np.array( [ [ix_px, iy_px, iz_um] for ix_px, iy_px, iz_um, ne in zip(x_px, y_px, z_um, nE) for _ in range(int(ne)) ] ).T
        iex_px, iey_px = round2int( *diffuseElectrons( ex_px, ey_px, ez_um ) )
        inside_boudary = np.all( [ iex_px >= 0, iex_px < xlen_px, iey_px >= 0 , iey_px < ylen_px ], axis=0 )
        image = np.zeros( self.config.imageShape )
        np.add.at( image, [ iex_px[inside_boudary], iey_px[inside_boudary] ], 1 )
        
        if save:
            saveImage( image, path )
        print( '%.2fs'%(timeit.default_timer() - t), )
        return image

    def makeDarkCurrentImage_e( self, save = False, lambda_e_h = None, exposure_h = 3 ):
        self.partial_id += self.darkCurrent_id
        print( '+dc', )
        t = timeit.default_timer()

        xlen_px, ylen_px = self.config.imageShape
        ylen_px -= self.config.overScanCols
        xlen_px -= self.config.overScanLines
        
        numberElectrons = np.zeros( (xlen_px, ylen_px) )
        #mu = self.config.darkCurrent * self.exposure
        mu = lambda_e_h * exposure_h
        
        numberElectrons0 = scipy.stats.poisson.rvs( mu, size = xlen_px*ylen_px )
        totalElectrons = np.sum( numberElectrons0 )

        electronPositions = np.repeat( np.arange(xlen_px*ylen_px), numberElectrons0 )
        electronPositions = np.unravel_index( electronPositions, (xlen_px, ylen_px) )

        z_um = self.config.heightCCD * np.random.rand( totalElectrons )

        x_px, y_px = electronPositions
        ix_px, iy_px = round2int( *diffuseElectrons( x_px, y_px, z_um ) )

        inside_boudary = np.all( [ ix_px >= 0, ix_px < xlen_px, iy_px >= 0 , iy_px < ylen_px ], axis=0 )
        image = np.zeros( (xlen_px, ylen_px) )

        np.add.at( image, [ ix_px[inside_boudary], iy_px[inside_boudary] ], 1. )
        
        print( '%.2fs'%(timeit.default_timer() - t), )
        return image
    
    def getNoiseImage_e( self ):
        self.partial_id += self.noise_id
        print( '+rn', )
        t = timeit.default_timer()

        nx, ny = self.rebin
        imageShape = self.config.imageShape
        xlen, ylen = round2next( imageShape[0]/nx, nx ), round2next( imageShape[1]/ny, ny )
        
        e = self.config.readoutNoiseWidth * scipy.stats.norm.rvs( size=xlen*ylen ).reshape( xlen, ylen )
        print( '%.2fs'%(timeit.default_timer() - t), )
        return e
    
    def getNoiseImage_adu( self ):
        return self.getNoiseImage_e*self.config['e-[keV]']*self.adu_keV
        
    def saveGain( self, runID, ohdu, gain ):
        fname = self.config.outPath + 'gainCu%s_%s.dat'%(runID, ohdu)
        print( 'saving gain', runID, ohdu, gain )
        data = { 'gainCu': float(gain) }
        json.dump( data, open(fname, 'w') )
        return

    @staticmethod
    def _readBaselineNoise( runID, ohdu ):
        image = Image( runID = runID, ohdu = ohdu )
        R = image.R()
        L = image.L()
        return [ L.os().mean(), L.osSub().os().std(), L.os().median(), L.osSub().os().mad() ], [ R.os().mean(), R.osSub().os().std(), R.os().median(), R.osSub().os().mad() ], [ L.mean(), L.osSub().std(), L.median(), L.osSub().std() ], [ R.mean(), R.osSub().std(), R.median(), R.osSub().std() ]
        
    @staticmethod
    def _readGain( runID, ohdu ):
        scnFile = listFITSscn_runID( int(runID) )
        if len(scnFile) == 0: return None
        scnFile = scnFile[0]
        associatedCatalog = getAssociatedCatalogGain( scnFile )[0]
        catalogData = readCatalog( associatedCatalog, int(runID) )
        hits_adu = catalogData[ catalogData['ohdu'] == int(ohdu) ]
        return hits_adu['gainCu'][0]
    
    def readGain( self, runID, ohdu ):
        fname = self.config.outPath + 'gainCu%s_%s.dat'%(runID, ohdu)
        if not os.path.exists(fname):
            self.getSCNHits_adu(runID, ohdu)
            scnFile = listFITSscn_runID( int(runID) )[0]
            associatedCatalog = getAssociatedCatalogGain( scnFile )[0]
            catalogData = readCatalog( associatedCatalog, int(runID) )
            hits_adu = catalogData[ catalogData['ohdu'] == int(ohdu) ]
            self.adu_keV = hits_adu['gainCu'][0]
            self.saveGain( runID, ohdu, self.adu_keV )
            return self.adu_keV
        data = json.load( open( fname, 'r') )
        self.adu_keV = data['gainCu']
        return self.adu_keV

    @staticmethod
    def _getOSIImage_adu( runID, ohdu ):
        print( 'runID', runID, 'ohdu', ohdu, )
        rawFiles = listFITSosi_runID( int(runID) )
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        image = np.concatenate( [ ohdus(astropy.io.fits.open(rawFile)) for rawFile in rawFiles[::-1] ], axis=0 )
        return image

    @staticmethod
    def _getMBSImage_adu( runID, ohdu ):
        print( 'runID', runID, 'ohdu', ohdu, )
        rawFiles = listFITSmbs_runID( int(runID) )
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        image = np.concatenate( [ ohdus(astropy.io.fits.open(rawFile)) for rawFile in rawFiles[::-1] ], axis=0 )
        return image

    @staticmethod
    def _getRawImage_adu( runID, ohdu ):
        print( 'runID', runID, 'ohdu', ohdu, )
        rawFiles = listFITSraw_runID( int(runID) )
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        image = np.concatenate( [ ohdus(astropy.io.fits.open(rawFile)) for rawFile in rawFiles[::-1] ], axis=0 )
        return image

    def getRawImage_adu( self, runID=None, save = False, ohdu=None ):
        if ohdu is None: ohdu = self.config.ohdu
        if runID is None: runID = self.config.runID
        print( 'runID', runID, 'ohdu', ohdu, )
        path = self.config.outPath + self.raw_id + str(runID) + self.image_ext
        if os.path.exists( path ):
            return readImage( path )
        rawFiles = listFITSraw_runID( int(runID) )
        return hitSumm._getRawImage_adu( runID, ohdu )
    
    @staticmethod
    def _getSCNImage_adu( runID, ohdu ):
        print( 'runID', runID, 'ohdu', ohdu, )
        scnFile = listFITSscn_runID( int(runID) )[0]
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        #print( astropy.io.fits.open(scnFile).header )
        image = ohdus(astropy.io.fits.open(scnFile))
        #print( 'scn.shape', image.shape )
        return image

    def getSCNImage_adu( self, runID=None, ohdu=None ):
        if runID is None: runID = self.config.runID
        if ohdu is None: ohdu = self.config.ohdu
        print( 'runID', runID, 'ohdu', ohdu, )
        path = self.config.outPath + self.scn_id + str(runID) + self.image_ext
        if os.path.exists( path ):
            return readImage( path )
        scnFile = listFITSscn_runID( int(runID) )[0]
        return hitSumm._getSCNImage_adu( runID, ohdu )

    def getSCNHits_adu( self, runID, ohdu=None, save = False ):
        self.partial_id += self.scn_id + str(runID)
        print( '+SCN', )
        t = timeit.default_timer()
        
        if ohdu is None: ohdu = int(self.config.ohdu)
        path = self.config.outPath + self.scn_id + str(runID) + self.hits_ext
        if os.path.exists( path ):
            return readTable( path )
        scnFile = listFITSscn_runID( int(runID) )[0]
        associatedCatalog = getAssociatedCatalogGain( scnFile )[0]
        catalogData = readCatalog( associatedCatalog, int(runID), selection=self.config.selection )
        hits_adu = catalogData[ catalogData['ohdu'] == ohdu ]
        self.adu_keV = hits_adu['gainCu'][0]
        self.saveGain( runID, ohdu, self.adu_keV )
        #print( 'getSCNHits_adu', 'len', len(hits_adu['E0']), 'runID', runID )
        print( '%.2fs'%(timeit.default_timer() - t), )
        return hits_adu

    def getMBSubtracted_adu( self, image ):
        median = np.median( image.flatten() )
        image -= median
        return image
    
    def filterExposure( self, hits ):
        self.partial_id += self.exposure_id
        print( '+exp', )
        if hits is None: 
            return None
        indexlist = range( len( hits ) )
        fraction = float( self.exposure )/self.config.exposureTimeBase
        indexlist = random.sample(range( len( hits )), int( fraction*len(hits) ) )
        hitsFiltered = hits[ indexlist ]
        return hitsFiltered
        
    def getImage( self ):
        path = self.config.outPath + self.full_id + self.image_ext
        if os.path.exists( path ):
            return readSingleImageFITS( '', path )
        
    def convertHits2Image( self, hits, withMask = True ):
        matrix = np.zeros( self.config.imageShape )
        for hit in hits:
            if withMask:
                mask = np.all( [ hit['level'] <= self.config.excludedFromLevelUp, hit['ePix'] < self.config.excludedFromEnergyUp, hit['ePix'] > self.config.excludedFromEnergyDown ], axis=0 )
                try:
                    matrix[ hit['yPix'][mask], hit['xPix'][mask] ] += hit['ePix'][mask]
                except:
                    None
            else:
                matrix[ hit['yPix'], hit['xPix'] ] += hit['ePix']
        return matrix
        
    def convertImage2Hits( self, image, save = None, verbose = False ):
        return processImage3( image, self.config.extractThreshold )
    
        #i = self.partial_id + self.image_ext
        #writeImageHits( image, self.config(), fname = i )
        #if save is not None:
            #o = self.partial_id + self.hits_ext
        #else:
            #o = self.config.outPath + 'temp.root'
        #if not os.path.exists(o):
            #os.system( '/share/storage2/connie/data_analysis/Processor/tools/extract/extract.exe %s -c /home/mota/tests/simulatedFITS/extractConfigFS2.xml -o %s'%(i,o) + ( '>/dev/null' if not verbose else '') )
            #print( 'ended' )
        
        #os.remove( i )
        #ohduData = None
        #try:
            #catalogData = readCatalog( o, None, readGain=False )
            #ohduData = catalogData[ catalogData['ohdu'] == int(self.config.ohdu) ]
        #except:
            #pass
        #if save is None:
            #os.remove( o )
        #return ohduData
    
    def rebinImage( self, image ):
        self.partial_id += self.rebin_id
        print( '+rb', )
        t = timeit.default_timer()

        xlen, ylen = image.shape
        nx, ny = self.rebin
        imageR = np.zeros( ( nx*round2next(image.shape[0]/nx, nx ), ny*round2next(image.shape[1]/ny, ny ) ) )
        imageR[:image.shape[0], :image.shape[1] ] = image
        
        imageR = imageR.reshape([-1, nx, imageR.shape[1]] ).sum( axis = 1 )
        imageR = imageR.reshape([ imageR.shape[0], -1, ny ] ).sum( axis = 2 )
        
        print( '%.2fs'%(timeit.default_timer() - t), )
        return imageR

    def makeSimulatedImage_e( self, lambda_e_h = None, exposure_h = None, numberCu = None, numberCu2 = None, numberSi = None ):
        print( 'making', path )
        
        simulatedImage_e = np.zeros( self.config.imageShape )
        
        if self.darkCurrent:
            simulatedImage_e += self.makeDarkCurrentImage_e( lambda_e_h = lambda_e_h, exposure_h = exposure_h )
            
        if self.neutrinos:
            simulatedImage_e += self.getNeutrinosImage_e()

        if self.copper:
            simulatedImage_e += self.getCopperImage_e('Cu', number = self.config.numberCu )

        if self.copper2:
            simulatedImage_e += self.getCopperImage_e('Cu2')

        if self.silicon:
            simulatedImage_e += self.getCopperImage_e('Si')
            
        return simulatedImage_e
    
    def makeSimulatedImage_adu( self, adu_keV ):
        return self.makeSimulatedImage_e()*self.config['e-[keV]']*adu_keV

    def reconstructSCNImage_adu( self, runID, ohdu ):
        hits_adu = self.getSCNHits_adu( runID, ohdu, save = True )
        if self.exposure:
            hits_adu = self.filterExposure( hits_adu )
        hitsImage_adu = self.convertHits2Image( hits_adu )

        if self.rebin:
            hitsImage_adu = self.rebinImage( hitsImage_adu )

        return hitsImage_adu

    def getHits_adu(self, runID, ohdu, adu_keV, save=False ):
        print( 'image2hit' )
        path = self.outPath + self.partial_id
        simulatedImage_adu = self.getSimulatedImage_adu( adu_keV )
        if not self.empty:
            simulatedImage_adu += self.reconstructSCNImage_adu( runID, ohdu )
        hits_adu = self.convertImage2Hits( simulatedImage_adu, save = True )
        if save:
            root_numpy.array2root(hits_adu, path, treename='hitSumm', mode='recreate')
        return hits_adu

    def spectrum(self):
        print( 'spectrum', )
        if self.empty:
            print( 'empty' )
            hits_adu = self.getHits_adu( None )
            if hits_adu is None:
                return None
            return hits_adu['E2']/self.adu_keV
        E2 = []
        for runID in range( self.config.runIDmin, self.config.runIDmax+1):
            hits_adu = self.getHits_adu( runID )
            if hits_adu is None:
                continue
            E2.extend( hits_adu['E2'] )
        return np.array(E2)/self.adu_keV

    def spectrumDiffPeaks( self, bins ):
        histDiff = self.spectrumDiff( bins )
        peaks = find_peaks( histDiff, radius=2, tol=1 )
        xpeaks = .5*(bins[1:]+bins[:-1])[ peaks ]
        return xpeaks, histDiff[ peaks ]

    def spectrumHist( self, bins, normed = False ):
        # N/day = N/24h = 1/12 N/2h
        if normed:
            #timeFactor = float(self.exposure)/24
            norm = self.config.CCDmass * float(self.exposure) * self.config.spectrumBinsSize
            norm = norm/24
            #norm2 = 1e-3 * 24/float(self.exposure) * self.config.spectrumBinsSize
        else:
            norm = 1
        norm *= (self.config.runIDmax - self.config.runIDmin)+1
        spectrum = self.spectrum()
        if spectrum is None:
            return None
        return np.histogram( spectrum, bins )[0]/norm
    
    def reconstructionEfficiency( self ):
        scn_hits = self.getSCNHits_adu( self.config.runIDmin )
        print( 'scn_hits', len(scn_hits), 'x', np.min(scn_hits['xMin']), np.max(scn_hits['xMax']), 'y', np.min(scn_hits['yMin']), np.max(scn_hits['yMax']), 'E', np.min(scn_hits['E0']), np.min(scn_hits['E1']), np.min(scn_hits['E2']), np.min(scn_hits['E3']) )
        scn_image = hitSumm.Image( self.getSCNImage_adu( self.config.runIDmin ) )
        print( 'scn.shape2', scn_image.image.shape, scn_image.image.shape, )
        rec_hits = self.convertImage2Hits( scn_image.image )
        print( 'rec_hits', len(rec_hits) )
        
        center = lambda hits: ( hits['xBary0'], hits['yBary0'] )
        center2 = lambda hits: ( hits['yBary0'], hits['xBary0'] )
        dist = lambda a,b: np.sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 )
        
        matches = 0.
        print( 'not found:' )
        for scn_hit in scn_hits:
            match = np.argwhere( dist( center(rec_hits), center2(scn_hit) ) < 2 )
            if len(match) == 1:
                matches += 1.
            #else:
                #match = np.argwhere( np.abs(scn_hit['E0'] - rec_hits['E0'])/scn_hit['E0'] < 1e-3 )
                #if len(match) > 0:
                    #print( scn_hit['E0'], rec_hits['E0'][match], center2(scn_hit), [ center(rh) for rh in rec_hits[match] ] )
                ###print( scn_hit['E0'], )
        print( '%.2f'%(matches/len(scn_hits)), 'scn_hits found' )

        matches = 0.
        print( '\nnot found:' )
        for rec_hit in rec_hits:
            match = np.argwhere( dist( center(rec_hit), center2(scn_hits) ) < 2 )
            if len(match) == 1:
                matches += 1.
        print( '%.2f'%(matches/len(rec_hits)), 'rec_hits found' )
        
    def efficiency( self, bins ):
        hits_adu = self.getHits_adu( self.config.runIDmin )
        if hits_adu is None:
            return 0, 0, 0, 0, 0, 0
        
        print( 'all', len(hits_adu) )
        neutrinoHits = readTable( self.config.outPath+self.neutrino_id+self.table_ext )
        print( 'nu', len(neutrinoHits) )
        #neutrinoHits2_adu = self.convertImage2Hits( self.getNeutrinosImage_adu(), threshold_adu=2 )
        #print( 'nuR', len(neutrinoHits2_adu) )
        
        center = lambda hits: ( hits['xBary2'], hits['yBary2'] )
        center2 = lambda hits: ( hits[0]/self.rebin[0], hits[1]/self.rebin[1] )
        dist = lambda a,b: np.sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 )

        reconstructedNeutrinos = np.zeros_like( bins[:-1] )
        reconstructedNeutrinosE2 = np.zeros_like( bins[:-1] )
        totalNeutrinos = np.zeros_like( bins[:-1] )
        #totalNeutrinos2 = np.histogram( neutrinoHits2_adu['E2']/self.adu_keV, bins )[0]
        totalNeutrinos2 = None
        totalEvents = np.histogram( hits_adu['E2']/self.adu_keV, bins )[0]
        
        matches = 0
        extraMatches = 0
        Enu = []
        Erec = []
        for neutrinoHit in neutrinoHits:
            match = np.argwhere( dist( center(hits_adu), center2(neutrinoHit) ) < self.config.matchDistanceTolerance )
            totalNeutrinos[ np.all( [ neutrinoHit[4]*self.config['e-[keV]']>bins[:-1], neutrinoHit[4]*self.config['e-[keV]']<bins[1:] ], axis=0 ) ] += 1
            #totalNeutrinos[ np.all( [ neutrinoHit[3]>bins[:-1], neutrinoHit[3]<bins[1:] ], axis=0 ) ] += 1
            if len( match ) == 1:
                matches += 1
                reconstructedNeutrinosE2[ np.all( [ hits_adu[match[0]]['E2']/self.adu_keV>bins[:-1], hits_adu[match[0]]['E2']/self.adu_keV<bins[1:] ], axis=0 ) ] += 1.
                Enu.append( neutrinoHit[4]*self.config['e-[keV]'] )
                Erec.append( hits_adu[match[0]]['E2']/self.adu_keV )
                reconstructedNeutrinos[ np.all( [ neutrinoHit[4]*self.config['e-[keV]']>bins[:-1], neutrinoHit[4]*self.config['e-[keV]']<bins[1:] ], axis=0 ) ] += 1.
            elif len( match ) > 1:
                extraMatches += 1
        print( 'match percentage', float(matches)/len(neutrinoHits), 'extra', extraMatches, float(extraMatches)/len(neutrinoHits) )
        rate = np.ones_like(reconstructedNeutrinos)
        rate[totalNeutrinos>0] = reconstructedNeutrinos[totalNeutrinos>0]/totalNeutrinos[totalNeutrinos>0]
        
        neutrinoRate = readTable('simulatedFITS/neutrinoRate.dat').T
        neutrinoFunction = lambda x: np.interp( x, neutrinoRate[0], neutrinoRate[1] )

        return reconstructedNeutrinos, reconstructedNeutrinosE2, totalNeutrinos, totalNeutrinos2, Enu, Erec

class RunIDImages:
    def __init__( self, runIDs=None, ohdus=None, imageType='raw' ):
        #images = { runID: { ohdu: Image( runID=runID, ohdu=ohdu ) for ohdu in ohdus } for runID in runIDs }
        #self.lefts = { runID: { ohdu: images[runID][ohdu].L().biasImageSub() for ohdu in ohdus } for runID in runIDs }
        #self.rights = { runID: { ohdu: images[runID][ohdu].R().biasImageSub() for ohdu in ohdus } for runID in runIDs }
        self.imageType = imageType
        self.rights = {}
        self.lefts = {}
        self.ohdus = ohdus
        self.runIDs = runIDs
        self.iM = { runID: None for runID in runIDs }
    
    def loadImages( self, runID ):
        images = { ohdu: Image( runID=runID, ohdu=ohdu, imageType=self.imageType ) for ohdu in self.ohdus }
        if self.imageType == 'raw':
            self.lefts[runID] = { ohdu: images[ohdu].L().biasImageSub() for ohdu in self.ohdus }
            self.rights[runID] = { ohdu: images[ohdu].R().biasImageSub() for ohdu in self.ohdus }
        elif self.imageType == 'scn':
            self.lefts[runID] = { ohdu: images[ohdu] for ohdu in self.ohdus }
        return self

    def unloadImages( self, runID ):
        self.lefts.pop(runID,0)
        self.rights.pop(runID,0)
        return self
        
    def computeDarkImages( self ):
        self.dark = { ohdu: Image( np.median( [ self.lefts[runID][ohdu].image for runID in self.runIDs ], axis=0 ) ).printImage('dark_%s.pdf'%ohdu) for ohdu in self.ohdus }
        self.minsLeft = { ohdu: Image( np.min( [ self.lefts[runID][ohdu].image for runID in self.runIDs ], axis=0 ) ).printImage('minsL_%s.pdf'%ohdu) for ohdu in self.ohdus }
        self.bias = { ohdu: Image( np.median( [ self.rights[runID][ohdu].image for runID in self.runIDs ], axis=0 ) ).printImage('bias_%s.pdf'%ohdu) for ohdu in self.ohdus }
        self.minsRight = { ohdu: Image( np.min( [ self.rights[runID][ohdu].image for runID in self.runIDs ], axis=0 ) ).printImage('minsR_%s.pdf'%ohdu) for ohdu in self.ohdus }
        return
    
    def computeMedianLeft( self ):
        return Image( np.median( [ image.image for image in self.lefts.values() ], axis=0 ) )

    def computeMedianRight( self ):
        return Image( np.median( [ image.image for image in self.rights.values() ], axis=0 ) )
    
    def computeMatrix( self, runID, side='left' ):
        if side == 'left':
            images = self.lefts
        elif side == 'right':
            images = self.rights
        else:
            print( 'unpredicted side', side )
            exit(0)
        print
        M = np.zeros((len(self.ohdus), len(self.ohdus))).astype(float)
        for i, ohdu1 in enumerate(self.ohdus):
            for j, ohdu2 in list(enumerate(self.ohdus))[i:]:
                M[i,j] = madCovSqr( images[runID][ohdu1].image, images[runID][ohdu2].image )
                if i!=j:
                    M[j,i] = M[i,j]
        return M
    
    def invMatrix( self, runID ):
        if self.iM[runID] is None:
            self.iM[runID] = np.linalg.inv( self.computeMatrix(runID, side='right') )
        return self.iM[runID]

    def correctedLeft( self, runID, ohdu ):
        a = self.invMatrix( runID ).dot( [ madCovSqr(self.lefts[runID][ohdu].image, self.rights[runID][ohdu1].image) for ohdu1 in self.ohdus ] )
        return Image( self.lefts[runID][ohdu].image - np.array( [ self.rights[runID][ohdu1].image for ohdu1 in self.ohdus ] ).transpose((1, 2, 0)).dot(a) )#/np.sum(a) )

    def correctedLeftSimple( self, ohdu ):
        return Image( self.lefts[ohdu].image - np.std(self.lefts[ohdu].os().image)/np.std(self.rights[ohdu].os().image)*self.rights[ohdu].image )

    def generateCatalog( self, thr=15, eMin=15, eMax=2000*10, binWidth=100 ):
        first = True
        for runID in self.runIDs:
            self.generateCatalog_runID( runID, thr, eMin, eMax, binWidth, 'recreate' if first else 'update' )
            first = False
        return
    
    def generateCatalog_runID( self, runID, thr=15, eMin=15, eMax=2000*10, binWidth=100, mode='update' ):
        self.loadImages( runID )
        for ohdu in self.ohdus:
            output = self.rootname( ohdu )
            if self.imageType == 'raw': image = self.correctedLeft( runID, ohdu )
            elif self.imageType == 'scn': image = self.lefts[runID][ohdu]
            hits = processImage4( image.image, threshold=thr )
            hits = append_fields( hits, ('runID', 'ohdu' ), (np.ones(len(hits))*runID, np.ones(len(hits))*ohdu), dtypes = ['i4','i4'] )
            root_numpy.array2root( hits, output, treename='hitSumm', mode=mode )
            print( 'root', output, 'mode=%s'%mode )
            RunIDImages.printSpectrum_ohdu( self, ohdu, eMin, eMax, binWidth )
        self.unloadImages( runID )

    def rootname( self, ohdu ):
        return 'processed4_%s_runIDs%s_%s_ohdu%s.root'%(self.imageType, np.min(self.runIDs), np.max(self.runIDs), ohdu )
    
    def figname( self, ohdu ):
        return 'processed4_%s_runIDs%s_%s_ohdu%s.pdf'%(self.imageType, np.min(self.runIDs), np.max(self.runIDs), ohdu )
    
    @staticmethod
    def printSpectrum( images, eMin, eMax, binWidth, fname=None ):
        if type(images) is not list: images = [images]
        for ohdu in images[0].ohdus:
            if images[0].imageType == 'root':
                catalogs = getAssociatedCatalog( listFITSscn_runID( images[0].runIDs[0] ) )
                selection = 'runID>=%d && runID<=%d && ohdu==%d'%( np.min(images[0].runIDs), np.max(images[0].runIDs), ohdu )
                print( 'selection', selection, )
                RunIDImages.printSpectrumCatalog_ohdu( catalogs, ohdu, eMin, eMax, binWidth, fname=fname(ohdu), selection=selection )
            else:
                for image in images:
                    if not os.path.exists( image.rootname(ohdu) ): image.generateCatalog()
                RunIDImages.printSpectrum_ohdu( images, ohdu, eMin, eMax, binWidth, fname=fname )
        return

    @staticmethod
    def printSpectrum_ohdu( images, ohdu, eMin, eMax, binWidth, fname=None ):
        if type(images) is not list: images = [images]
        if fname is None: fname = images[0].figname(ohdu)
        else: fname = fname(ohdu)
        RunIDImages.printSpectrumCatalog_ohdu( [ image.rootname(ohdu) for image in images ], ohdu, eMin, eMax, binWidth, fname=fname )
        
    @staticmethod
    def printSpectrumCatalog_ohdu( catalogs, ohdu, eMin, eMax, binWidth, fname=None, selection=None ):
        if catalogs[0].startswith('processed4_'):
            branches = ['E','level','runID','ohdu']
            E0 = lambda data: data['E'][data['level']==0]
            E1 = lambda data: data['E'][data['level']==1]
            E2 = lambda data: data['E'][data['level']==2]
            E3 = lambda data: data['E'][data['level']==3]
        else:
            branches = ['E0','E1','E2','E3','runID','ohdu']
            E0 = lambda data: data['E0']
            E1 = lambda data: data['E1']
            E2 = lambda data: data['E2']
            E3 = lambda data: data['E3']
        data = [ root_numpy.root2array( catalog, treename='hitSumm', selection=selection, branches=branches ) for catalog in catalogs ]

        for datum in data:
            runIDs = np.unique( datum['runID'] )
            print( 'runIDs', runIDs, len(runIDs) )
            ohdus = np.unique( datum['ohdu'] )
            print( 'ohdus', ohdus, len(ohdus) )
            
        bins = np.arange( eMin, eMax, binWidth )
        fig = plt.figure()
        fig.suptitle('ohdu %s'%ohdu )
        ax = fig.add_subplot(411)
        for datum in data:
            print( 'len', len(E0(datum)) )
            ax.hist( E0(datum), bins=bins, label='E0 %s-%s'%( np.min(datum['runID']), np.max(datum['runID']) ), histtype='step' )
        ax.legend( fancybox=True, framealpha=0 )
        ax.set_yscale('log')
        ax.grid(True)
        
        ax = fig.add_subplot(412)
        for datum in data:
            print( 'len', len(E1(datum)) )
            ax.hist( E1(datum), bins=bins, label='E1 %s-%s'%( np.min(datum['runID']), np.max(datum['runID']) ), histtype='step' )
        ax.legend( fancybox=True, framealpha=0 )
        ax.set_yscale('log')
        ax.grid(True)

        ax = fig.add_subplot(413)
        for datum in data:
            print( 'len', len(E2(datum)) )
            ax.hist( E2(datum), bins=bins, label='E2 %s-%s'%( np.min(datum['runID']), np.max(datum['runID']) ), histtype='step' )
        ax.legend( fancybox=True, framealpha=0 )
        ax.set_yscale('log')
        ax.grid(True)

        ax = fig.add_subplot(414)
        for datum in data:
            print( 'len', len(E3(datum)) )
            ax.hist( E3(datum), bins=bins, label='E3 %s-%s'%( np.min(datum['runID']), np.max(datum['runID']) ), histtype='step' )
        ax.legend( fancybox=True, framealpha=0 )
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_xlabel('ADU')
        
        fig.subplots_adjust(hspace=0, wspace=.05)
        
        fig = plt.savefig( fname )
        print( 'savefig', fname )
        return
            
    
def stdCovMed( data0, data1 ):
    return np.sum( ( data0 - np.median( data0 ))*( data1 - np.median( data1 )) )/len(data0)
    
def madCovSqr( data0, data1 ):
    return np.median( np.abs( (data0 - np.median(data0))*((data1 - np.median(data1))) ) )

def normFit( data, loc=None, scale=None ):
    return scipy.stats.norm.fit( data, floc=loc, fscale=scale )

def mad( data ):
    return np.median( np.abs( data - np.median(data) ) )

class Image:
    def __init__( self, image=None, runID=None, ohdu=None, imageType='raw' ):
        if image is not None:
            self.image = image.astype(float)
        elif imageType == 'raw':
            self.image = hitSumm._getRawImage_adu( runID=runID, ohdu=ohdu )[::-1,:]
            self.image = self.image.astype(float)
        elif imageType == 'osi':
            self.image = hitSumm._getOSIImage_adu( runID=runID, ohdu=ohdu )[::-1,:]
            self.image = self.image.astype(float)
        elif imageType == 'mbs':
            self.image = hitSumm._getMBSImage_adu( runID=runID, ohdu=ohdu )[::-1,:]
            self.image = self.image.astype(float)
        elif imageType == 'scn':
            self.image = hitSumm._getSCNImage_adu( runID=runID, ohdu=ohdu )[::-1,:]
            self.image = self.image.astype(float)
        else:
            print( 'imageType not predicted', self.imageType )
            exit(0)
        
        self.shape = self.image.shape
        self.T = self.image.T
        self.overscanWidth = 0
        if imageType == 'section':
            return
        
        self.overscanWidth = 150
        
        if self.image.shape[1] == 8540 and imageType == 'raw':
            self.right = np.s_[:, 8:4270]
            self.left =  np.s_[:, 8:4270]

        elif self.image.shape[1] == 8524 and (imageType == 'osi' or imageType == 'mbs'):
            self.right = np.s_[:, :4262]
            self.left =  np.s_[:, :4262]
            
        elif self.image.shape[1] == 4262 and imageType == 'scn':
            self.overscan = np.s_[:, -self.overscanWidth: ]
            self.verticalOS = np.s_[ -100:,:]
            self.active = np.s_[:-100, :-self.overscanWidth]
            
        elif self.image.shape[1] == self.overscanWidth:
            pass
        
        elif self.image.shape[1] == 4112:
            self.overscanWidth = 0
        
        elif self.image.shape[1] == 9340:
            self.right = np.s_[:, 8:4670]
            self.left =  np.s_[:, 8:4670]
            
        elif self.image.shape[1] == 4662:
            self.overscanWidth = 550
            self.overscan = np.s_[:, -self.overscanWidth:]
            self.verticalOS = np.s_[-100:,:]
            self.active = np.s_[:-100,:-self.overscanWidth]
            
        elif self.image.shape[1] == 550:
            pass
        else:
            print( 'shape not predited', self.image.shape )
            exit(0)
    
    def __call__( self, mask ):
        return Image( self.image[mask], imageType='section' )
    def flatten( self ):
        return self.image.flatten()
    
    def R( self ):
        return Image( self.image[:,::-1][self.right] )
    
    def L( self ):
        return Image( self.image[self.left] )

    def part(self, i):
        return Image( self.image[i*1055:(i+1)*1055,:] )
    
    def ac( self ):
        return Image( self.image[self.active] )
    
    def os( self ):
        return Image( self.image[self.overscan] )
    
    def vos( self ):
        return Image( self.image[-50:,:] )
    
    def sub( self ):
        return Image( self.L().image - self.R().image )
    
    def subSmooth( self, radius, func=np.mean ):
        return Image( self.image - self.smooth2D( radius, func ).image )
    
    def area( self ):
        return self.image.shape[0]*self.image.shape[1]
    
    def smartSub( self ):
        active_hlines = np.median( self.image[:, 20*(-150):-150], axis=1 )
        overscan_hlines = np.mean( self.image[:, -150:], axis=1 )
        diff_hlines = active_hlines - overscan_hlines
        print( '\n(meanMed %.2f %.2f)'%( np.mean( diff_hlines ), np.std( diff_hlines ) ),  )
        
        image = Image( self.image - overscan_hlines[:,None] )
        
        active_hlines = np.median( image.image[:, 20*(-150):-150], axis=0 )
        overscan_hlines = np.mean( image.image[:, -150:], axis=0 )
        #diff_hlines = active_hlines - overscan_hlines
        print( '(meanMed %.2f %.2f)'%( np.mean( active_hlines ), np.std( active_hlines ) ), )
        print( '(meanMed %.2f %.2f)'%( np.mean( overscan_hlines ), np.std( overscan_hlines ) ), )
        return image
    
    @staticmethod
    def lineFunction( x, x0, x1, a, b, c ):
        y = np.zeros_like(x)
        y[x<x0] = a*( x[x<x0] - x1 )**2 + b
        y[x>x0] = c
        return y

    def biasImageSub( self ):
        image = self.image
        bias = np.median( image, axis=1 )[:,None] + np.median( image, axis=0 )[None,:]
        return Image(image-bias).osSub(axis=None)
    
    def biasImageSub2( self ):
        right = self.R().image
        left = self.L().image
        
        bias = np.median( right, axis=1 )[:,None] + np.median( right, axis=0 )[None,:]
        left = left - bias
        right = right - bias

        return Image(left).osSub(axis=None), Image(right).osSub(axis=None)
        
    def biasImage( self, func = np.median ):
        image = func( self.image, axis=1 )[:,None] + func( self.image, axis=0 )[None,:]
        image -= np.median( image )
        return Image(image)
    
    def vertSub( self ):
        med = np.median( self.image[:, 10*(-150):-150] )
        medOs = self.os().median()
        image = Image( self.image - np.median( self.image, axis=0 )[None,:] )
        image.image[:,:-150] += med
        image.image[:,-150:] += medOs
        return image

    def horiSub( self ):
        med = np.median( self.image[:, 10*(-150):-150] )
        medOs = self.os().median()
        image = Image( self.image - np.median( self.image, axis=1 )[:,None] )
        image.image[:,:-150] += med
        image.image[:,-150:] += medOs
        return image
        
    def vertR( self ):
        a = self.L().os().median()/self.R().os().median()
        return Image( a*np.tile( np.median( self.R().image, axis=0 ), self.image.shape[0] ).reshape(self.R().image.shape[0], self.R().image.shape[1]) )
        
    def shiftCenter( self, func=np.mean ):
        return Image( self.image - func( self.image, axis=1 )[:,None] )
    
    def osSub( self, func=np.median, axis=1 ):
        if axis == 1:
            return Image( self.image - func( self.os().image, axis = 1 )[:,None] )
        else:
            return Image( self.image - func( self.os().image ) )
    
    def hSub( self, func=np.median ):
        #return hitSumm.Image( (self.image - func( self.image, axis = 1 )[:,None])/func( self.os().image, axis = 1 )[:,None] )
        a = self.image - func( self.os().image, axis = 1 )[:,None]
        #r = self.R().image - func( self.R().os().image, axis = 1 )[:,None]
        return hitSumm.Image( a )
    
    def vosSub( self, func=np.median ):
        return Image( self.image - func( self.vos().image, axis = 0 )[None,:] )
    
    def hvosSub( self, func=np.median ):
        return hitSumm.Image( self.image - .5*(func( self.vos().image, axis = 0 )[None,:] + func( self.os().image, axis = 1 )[:,None] ) )
    
    def vos_osSub( self, func=np.median ):
        return self.vosSub( func = func ).osSub( func = func )

    def vos_osSubNorm( self, func=np.median ):
        return hitSumm.Image( self.vosSub( func = func ).osSub( func = func ).image/self.os().std()**2 )
    
    def vSub( self, func=np.median ):
        return hitSumm.Image( self.image - func( self.image, axis = 0 )[None,:] )
            
    def raw( self ):
        return hitSumm.Image( self.L().ac().image - hitSumm.normFit( self.L().os().image )[0] )
    
    def mask( self, m ):
        return hitSumm.Image( self.image[ m(self.image) ] )
    
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
    
    def mean(self, axis=None):
        return np.mean( self.image, axis=axis )
    
    def median(self, axis=None):
        return np.median( self.image, axis=axis )
    
    def std(self, axis=None):
        return np.std( self.image, axis=axis )
    
    def subtract(self, m ):
        return hitSumm.Image( self.image - m )
    
    def submean(self, axis=None ):
        return hitSumm.Image( self.image - self.mean() )
    
    def submedian(self, axis=None ):
        return hitSumm.Image( self.image - self.median() )
    
    def normFit(self):
        return np.array( hitSumm.normFit( self.image ) )
    
    def halfnormFit(self):
        return scipy.stats.halfnorm.fit( self.image )
    
    def subfit(self, axis=None ):
        return hitSumm.Image( self.image - self.normFit()[0] )
    
    def mad( self ):
        return mad( self.image )
    
    def extract(self, thr=15):
        return processImage3(self.image, threshold = thr)
    
    def electron2( self ):
        return hitSumm.normFit( self.negatives().image.flatten(), loc=0 )[1] - hitSumm.normFit( self.os().negatives().image.flatten(), loc=0 )[1]
    
    def osi_scn( self, gain=None ):
        #print( 'LosSub', self.L().osSub(np.median).vosSub().os().printImage('osLosSub.pdf' ).normFit() )
        #print( 'RosSub', self.R().osSub(np.median).vosSub().os().printImage('osRosSub.pdf' ).normFit() )
        a = np.cov( [ self.L().osSub(np.median).vosSub().os().image.flatten(), self.R().osSub().vosSub().os().image.flatten() ] )
        print( 'a', a, a[0,1]/a[1,1] )
        sub = hitSumm.Image( self.L().osSub().vosSub().image - a[0,1]/a[1,1]*self.R().osSub().vosSub().image )
        print( 'os', sub.os().printImage( 'os0.pdf' ) )
        print( 'ac', sub.printImage( 'ac0.pdf' ) )
        
        #print( 'Losmean', self.L().os().subfit().printImage( 'osLsubmean.pdf' ).normFit() )
        #print( 'Rosmean', self.R().os().subfit().printImage( 'osRsubmean.pdf' ).normFit() )
        aFull = np.sum( self.L().image*self.R().image )/np.sum( self.R().image**2 )
        print( 'a*', aFull, aFull*self.R().normFit(), aFull*self.R().os().normFit() )
        aOs = np.sum( self.L().os().image*self.R().os().image )/np.sum( self.R().os().image**2 )
        print( 'a**', aOs, aOs*self.R().normFit(), aOs*self.R().os().normFit() )
        am = self.L().os().mean()/self.R().os().mean()
        print( 'am', am, am*self.R().normFit(), am*self.R().os().normFit() )
        print( 'L', self.L().normFit(), self.L().os().normFit() )
        sub = hitSumm.Image( self.L().image - aFull*self.R().image )
        print( 'os', sub.os().printImage( 'os.pdf' ).normFit() )
        print( 'ac', sub.printImage( 'ac.pdf' ).normFit() )
        print( 'os', sub.osSub(np.median).os().printImage( 'osSub.pdf' ).normFit() )
        print( 'ac', sub.osSub(np.median).printImage( 'acSub.pdf' ).normFit() )
        return L.electron(), L2.electron()

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
        m0, s0 = self.median(), self.mad()
        print( '\n(mad %.2f %.2f %.2f %.2f)'%( m, s, m0, s0 ), )
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
        main_ax.invert_yaxis()
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

    def smooth2D( self, radius, func = np.mean ):
        return hitSumm.Image( hitSumm.smooth( self.image, radius=radius, func=func).reshape( self.image.shape ) )

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
    
    @staticmethod
    def polyFit( data, func=np.median, axis=1 ):
        line = func( data, axis=1 )
        x = range(data.shape[0])
        p = np.poly1d( np.polyfit( range(data.shape[0]), line, 15 ) )
        return p(x)
    
    @staticmethod
    def submad( data0, data1 ):
        return data0 - madCovSqr(data0,data1)/hitSumm.mad(data1)**2*data1
    
    @staticmethod
    def shiftCenter( data, loc=None ):
        if loc is None:
            return data - hitSumm.normFit( data )[0]
        else:
            return data - loc

    @staticmethod
    def smooth( data, radius, func = np.mean ):
        if radius == 0:
            return data
        if len(data.shape) > 1:
            data = data.flatten()
        startData = [ func( data[:radius+i], axis=0 ) for i in range(radius) ]
        centerData = func( [ data[i:i-2*radius if i < 2*radius else None] for i in range(2*radius+1) ], axis=0 )
        endData = [ func( data[i-2*radius:], axis=0 ) for i in range(radius) ]
        movingMean = np.concatenate( (startData, centerData, endData) )
        return movingMean
                
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
        
        #return [ [ x(iE), smooth2(np.abs(idN_dE),1) ], [] ]
        return [ [ 1./freq, np.log10(dN_df)], [ 1./freq0, np.log10(dN_df0)], [ 1./freq, np.log10(dN_df/dN_df0) ], [ 1./freqA, np.log10(dN_dfA) ] ]
    
        flatten = lambda l: [item for sublist in l for item in sublist]
        print( flatten([ [sl+'0', sl+'1'] for sl in slices ]) )
        return [ [ result[key[0]], result[key[1]] ] for key in [ [sl+'.E', sl+'.dN_dE'] for sl in slices ] ]
        #return [[range(len(result[key])), result[key]] for key in flatten([ [sl+'0', sl+'1'] for sl in slices ]) ]
        #return [[x(E), dN_dE], [] ]
        #return ((lambda_, dN_df))#(x(E), dN_dE0), (x(sigmas), dN_dsigma)
        #return (x(E), dN_dE), (lambda_, dN_df)#(x(E), dN_dE0), (x(sigmas), dN_dsigma)

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
    
class Plot:
    def __init__():
        self.lines = {}
    
    def add(self, page, axis, commandFunction, *args ):
        if not page in self.lines.keys():
            self.lines[page] = []
        self.lines[page] += [[ axis, commandFunction, args ]]
    
    def save(self, applytopages = lambda page: None, applytoaxes = lambda axis: None, ext = '.pdf', folder = '.' ):
        for pageKey, page in self.lines.items():
            fig = plt.figure()
            n = np.ceil( np.sqrt( len(page) ) ).astype(int)
            m = np.ceil( len(page)/n ).astype(int)
            for i, [ axis, axisFunc, args ] in enumerate( zip(page) ):
                ax = fig.add_subplot(n,m,i)
                axisFunc( ax, *args )
                applytoaxes(ax)
            applytopages(fig)
            fig.savefig( folder + pageName + ext )
            fig.close()
    
    @staticmethod
    def close():
        plt.close()

if __name__ == "__main__":
    print( term.bold('Dark Current Analysis Tool for the CONNIE colaboration') )
    print( term.bold('by Philipe Mota (philipe.mota@gmail.com)') )
    print( term.bold('repository https://github.com/PhMota/CONNIEtools/blob/master/fits_reader4.py') )
    if len(sys.argv) == 1:
        help_()
        exit(0)
    vars_ = {}
    outfolder = None
    ROOTfile = None
    fft = False
    if '--outfolder' in sys.argv:
        outfolder = sys.argv[sys.argv.index('--outfolder')+1]
        print( 'setting outfolder =', outfolder )
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
        print( 'setting ROOTfile =', ROOTfile )
    if '--fft' in sys.argv:
        fft = True
        print( 'setting fft = True' )
        
    if '--removehits' in sys.argv:
        if not vars_['runID'] is None:
            runETA( 'remove hits',
                lambda: removeHitsFromrunID( vars_['runID'], outputfolder=outfolder, ohdu = vars_['ohdu'], ROOTfile=ROOTfile )
                )
            exit(0)
        print( term.red('error: runID not set') )
    if '--plot' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA( 'analyse and plot', lambda: plotrunID( vars_['ohdu'], fft=fft, runID=vars_['runID'], inputpath=vars_['inputpath']) )
                exit(0)
            FITSfiles = sys.argv[sys.argv.index('--plot')+1:]
            print( 'input FITSfiles:' )
            print( '\n'.join(FITSfiles) )
            plotrunID( vars_['ohdu'], FITSfiles, fft=fft )
            exit(0)
        print( term.red('error: ohdu not set') )
    if '--table2' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA ( 'remove and analyse', lambda: plot2(outfolder, vars_['ohdu'], vars_['runID'], ROOTfile, plot=True, gain=vars_['gain'], crosstalk=vars_['crosstalk'], verbose=False, onlyoccupancy=vars_['onlyoccupancy'] ) )
                exit(0)
            if not vars_['run'] is None:
                l = sorted(list(set( listrunID_run( vars_['run'] ) )))
                print( 'from', min(l), 'to', max(l), ' total', len(l) )
                runETA( 'remove and analyse full run '+vars_['run'], 
                       cmd = lambda x: plot2(outfolder, vars_['ohdu'], int(x), ROOTfile, gain=vars_['gain'], verbose=False, crosstalk=vars_['crosstalk'], onlyoccupancy=vars_['onlyoccupancy'] ),
                       loop = l
                       )
                exit(0)
            exit(0)
        print( term.red('error: ohdu not set') )
    if '--table' in sys.argv:
        if not vars_['ohdu'] is None:
            if not vars_['runID'] is None:
                runETA( 'analyse', lambda: plotrunID( vars_['ohdu'], fft=fft, runID=vars_['runID'], plot=False, inputpath=vars_['inputpath'], gain=vars_['gain']) )
                exit(0)
            FITSfiles = sys.argv[sys.argv.index('--table')+1:]
            print( 'input FITSfiles:' )
            print( '\n'.join(FITSfiles) )
            plotrunID( vars_['ohdu'], FITSfiles, fft=fft, plot=False )
            exit(0)
        print( term.red('error: ohdu not set') )
    if '--plotcsv' in sys.argv:
        #plotCSV( rglob('tmp/*_occupancy.csv'), 'runID', 'N250350', xerr='bin', yerr='Poisson', title='count evolution 250--350keV' )
        plotCSV( rglob('tmp/*_occupancy.csv'), 'runID', 'N37', xerr='bin', yerr='Poisson', title='count evolution 3--7keV' )
        exit(0)
        print( term.red('error: ohdu not set') )
    if '--listrunID' in sys.argv:
        print( 'list of runID' )
        print( ', '.join( listrunID_run(vars_['run']) ) )
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
    if '--simulateneutrinos' in sys.argv:
        config = sys.argv[ sys.argv.index('-c') + 1 ]
        efficiencyNeutrinosReconstruction( config )
        #simulateNeutrinos( config )
        exit(0)
    args = {}
    for iarg, arg in enumerate(sys.argv):
        if arg.startswith('-') and not arg.startswith('--'):
            args[arg.strip('-')] = sys.argv[iarg+1]
    #print( args )
    for arg in sys.argv:
        if arg.startswith('--'):
            #print( Callable.__dict__ )
            getattr(Callable, arg.strip('--'))( **args )
            exit(0)
    print( 'callables', Callable.__dict__ )
    print( 'options not recognized', sys.argv )
    help_()
    exit(0)
