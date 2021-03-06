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
#matplotlib.use('Agg')
#matplotlib.use('qt4agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import patches, colors, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid
#from matplotlib.backends.backend_gtk3agg import (FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.backends.backend_gtk3cairo import (FigureCanvasGTK3Cairo as FigureCanvas)
from matplotlib.backends.backend_gtk3 import (NavigationToolbar2GTK3 as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.ticker as ticker

from functools import partial

import glob
import os
import shutil
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
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib
import threading

from ConnieDataPath import ConnieDataPath as ConniePaths
import MonitorViewer
import ImageViewer
import SimulateImage
import ConnieImage
import Statistics as stats
import Timer
import termcolor
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

#/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/data_*/runs/*/
def list_subruns( run='all' ):
    if run == 'all':
        l = glob.glob( '/share/storage2/connie/data_analysis/processed02_data/runs/*/' )
        #l.extend( glob.glob( '/share/storage2/connie/data_analysis/processed02_data/runs/*/' ) )
    else:
        l = glob.glob( '/share/storage2/connie/data_analysis/processed02_data/runs/%03d*/'%int(run) )
    subruns = [ re.search(r'runs/(.+)/', i ).groups()[0] for i in l ]
    return sorted(subruns)

def list_all_runIDs():
    l = glob.glob( '/share/storage2/connie/data/runs/*/runID_*_*_p1.fits.fz' )
    return l

def list_runIDs( *runs ):
    l = listFITS( *[ '/share/storage2/connie/data/runs/%s/runID_*_p1.fits.fz'%(run) for run in runs ] )
    return l

def list_runIDs_from_subrun( *runs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/scn/merged/scn_*.fits'%(run) for run in runs ] )
    return l

def list_catalogs_skim1( run ):
    catalogs = glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_[0-9]*_to_[0-9]*/ext/catalog/catalog_data_*.skim1.root'%run)
    return catalogs

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

def get_raw_path( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data/runs/*/runID_*_%05d_*.fits.fz'%run for run in runIDs ] )
    return l

def listrunID_run( *runs ):
    l = list_runIDs( runs )
    
    #print( '\n'.join(l) )
    return map( getrunIDFromPath, l )

def getrunIDBounds_from_run( run ):
    runIDs = map( int, listrunID_run( run ) )
    if len(runIDs) == 0: return None
    #print runIDs, run
    return (min(runIDs), max(runIDs))

def get_scn_path( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/scn_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/data_*/scn/images/scn_mbs_osi_runID_*_%05d_Int*.fits'%runID for runID in runIDs ] )
        l = [ il for il in l if 'to_6226/' not in il ]
    if len(l) == 0:
        return ['file not found']
    return l

def get_osi_path( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/osi/images/osi_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/data_*/osi/images/osi_runID_*_%05d_Int*.fits'%runID for runID in runIDs ] )
        l = [ il for il in l if 'to_6226/' not in il ]
    if len(l) == 0:
        return ['file not found']
    return l

def get_mbs_path( *runIDs ):
    l = listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/mbs/images/mbs_*runID_*_%05d_*.fits'%runID for runID in runIDs ] )
    if len(l) == 0:
        l = listFITS( *[ '/share/storage2/connie/nu_processing/temp_guille_paper/connie_proc/connie_*/data_*/mbs/images/mbs_osi_runID_*_%05d_Int*.fits'%runID for runID in runIDs ] )
        l = [ il for il in l if 'to_6226/' not in il ]
    if len(l) == 0:
        return ['file not found']
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

def getrunFromRunID( runID ):
    return getrunFromPath( get_raw_path(runID) )

def getsubrunFromRunID( runID ):
    return getrunFromPath( get_scn_path(runID) )

def getrunFromPath( path ):
    if type(path) is list: path=path[0]
    #print( '(from %s'%path, )
    run = re.search(r'runs/(.+?)/', path ).groups()[0]
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
    print 'input FITSfile =', FITSfile
    if verbose: print( 'output flag =', output )
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

def histogram_fits( data, bins, func, p0, tol = 1, bounds = ()  ):
    class Result: pass
    r = Result()
    r.hist, bins = np.histogram( data, bins )
    mask = r.hist > tol
    x = .5*(bins[:-1]+bins[1:])
    dx = bins[1]-bins[0]
    N = np.sum(r.hist[mask])*dx
    r.popt, pcov = scipy.optimize.curve_fit( lambda x, *a: N*func(x, *a), x[mask], r.hist[mask], p0, sigma=np.sqrt(r.hist[mask]), absolute_sigma=True, bounds=bounds )
    r.perr = np.sqrt(np.diag(pcov))
    r.chisq = stats.chisquare( r.hist[mask], N*func(x[mask], *r.popt), ddof = len(r.popt) )
    return r

#def histogram_fits_lsq( data, bins, func, p0, tol=1, bounds = () ):
    #class Result: pass
    #r = Result()
    #r.hist, bins = np.histogram( data, bins )
    #mask = r.hist > tol
    #N = np.sum(r.hist[mask])
    #x = .5*(bins[:-1]+bins[1:])
    #residual = lambda p: np.sum( (N*func(x[mask], *p) - r.hist[mask])**2/r.hist[mask] )
    
    #r.popt, pcov = scipy.optimize.curve_fit( , x[mask], r.hist[mask], p0, sigma=np.sqrt(r.hist[mask]), absolute_sigma=True, bounds=bounds )
    #r.perr = np.sqrt(np.diag(pcov))
    #r.chisq = stats.chisquare( r.hist[mask], N*func(x[mask], *r.popt), ddof = len(r.popt) )
    #return r

def separate_active_overscan( image ):
    shape = image.val.shape
    vos = 70
    if shape[1] in [4262,4350]:
        os = 150
    elif shape[1] in [4662]:
        os = 550
    elif shape[1] in [9340]:
        os = 550
        width = shape[1]/2
        return image.val[ :-vos, 8 : width - os ], image.val[ :-vos, width - os : width ]
    else:
        print 'shape not predicted', shape
        exit(1)
    overscan = image.val[ :-vos, -os: ]
    active = image.val[ :-vos, :-os ]
    return active, overscan

def divide_image( active, overscan, N, i):
    n = int(active.shape[0]/N)
    active = active[i*n:(i+1)*n,:]
    overscan = overscan[i*n:(i+1)*n,:]
    overscan = overscan[ np.all( [overscan.value != 1e10, overscan.value != REMOVE], axis=0) ]
    active = active[ np.all( [active.value != 1e10, active.value !=REMOVE], axis=0) ]
    return active, overscan

def image_MLE( active, overscan ):
    print 'image_MLE'
    mu_o = stats.ufloat( np.mean(overscan), errSqr=stats.errSqr_mean(overscan) )
    sigma_o = stats.ufloat( np.std(overscan), errSqr=stats.errSqr_std(overscan) )

    mu_a = stats.ufloat( np.mean(active), errSqr=stats.errSqr_mean(active) )
    sigma_a = stats.ufloat( np.std(active), errSqr=stats.errSqr_std(active) )
    
    g = (sigma_a**2 - sigma_o**2)/(mu_a - mu_o)
    lamb = (mu_a-mu_o)/g

    return mu_o, sigma_o, g, lamb

def image_MLE_robust( active, overscan, divide=2 ):
    print 'overscan class', overscan.__class__.__name__
    class Return: pass
    ret = Return()
    ret.g = stats.darray()
    ret.mu = stats.darray()
    ret.sigma = stats.darray()
    ret.lamb = stats.darray()

    for i in range(divide):
        active_, overscan_ = divide_image( active, overscan, divide, i )

        mu = overscan_.median()
        sigma = overscan_.MAD()
        
        mu_a = active_.median()
        sigma_a = active_.MAD()
        
        g = (sigma_a**2 - sigma**2)/(mu_a - mu)
        lamb = (mu_a - mu)/g
        
        print '===divide', i
        print 'mu', mu
        print 'sigma', sigma
        print 'lamb', lamb
        print 'g', g
        ret.mu.append(mu)
        ret.sigma.append(sigma)
        ret.g.append(g)
        ret.lamb.append(lamb)
    
    print '===done'
    print 'g.mean', ret.g.mean()
    #ret.g = stats.mean(ret.g)
    #print 'gain', ret.g
    stats.dfloat.set_conversion('e-', '%s'%ret.mu.unit_str(), ret.g.value )
    print 'MLE robust'
    for key, val in ret.__dict__.items():
        print val.__class__.__name__
        print key, val, val.asunit('e-'), val/ret.g
    
    return ret

def image_MLE_fit( active, overscan, dE=1. ):
    mle = image_MLE( active, overscan )
    
    func = lambda x, loc, scale, gain, mu: stats.poisson_norm.pdf( x, loc, scale, gain, mu, fix_loc=False )
    p0 = ( mle.mu.value, mle.sigma.value, mle.g.value, mle.lamb.value )
    
    loc_bounds = (mle[0][0]-mle[1][0], mle[0][0]+mle[1][0])
    scale_bounds = (mle[0][1]-mle[1][1], mle[0][1]+mle[1][1])
    gain_bounds = (mle[0][2]-mle[1][2], mle[0][2]+mle[1][2])
    lambda_bounds = (mle[0][3]-mle[1][3], mle[0][3]+mle[1][3])
    
    bounds = loc_bounds, scale_bounds, gain_bounds, lambda_bounds
    active, _ = get_regions(image)
    bins = np.arange( active.min(), active.max(), dE )
    active_fit = histogram_fits( active, bins, func, p0 = p0, bounds = zip( *bounds ), tol=10 )
    loc = stats.uncertanty( active_fit.popt[0], active_fit.perr[0]*active_fit.chisq)
    sigma = stats.uncertanty( active_fit.popt[1], active_fit.perr[1]*active_fit.chisq)
    gain = stats.uncertanty( active_fit.popt[2], active_fit.perr[2]*active_fit.chisq)
    lamb = stats.uncertanty( active_fit.popt[3], active_fit.perr[3]*active_fit.chisq)
    
    #print 'chi_sq', active_fit.chisq
    #print 'mu', loc
    #print 'sigma', sigma
    #print 'gain*lambda', gain*lamb
    #print 'gain', gain
    #print 'lambda', lamb
    return (loc.val, sigma.val, gain.val, lamb.val), (loc.err(), sigma.err(), gain.err(), lamb.err())

def image_MLE_fit_robust( active, overscan, dE = 1., divide=1 ):
    loc = []
    sigma = []
    gain = []
    lamb = []
    for i in range(divide):
        active_, overscan_ = divide_image( active, overscan, divide, i )
        mle = image_MLE_robust( active_, overscan_ )
        print 'mle', mle
        func = lambda x, loc_, scale_, gain_, mu_: stats.poisson_norm.pdf( x, loc_, scale_, gain_, mu_, fix_loc=False )
        p0 = ( mle[0].value, mle[1].value, mle[2].value, mle[0].value )
        
        loc_bounds = (mle[0][0]-mle[1][0], mle[0][0]+mle[1][0])
        scale_bounds = (mle[0][1]-mle[1][1], mle[0][1]+mle[1][1])
        gain_bounds = (mle[0][2]-mle[1][2], mle[0][2]+mle[1][2])
        lambda_bounds = (mle[0][3]-mle[1][3], mle[0][3]+mle[1][3])
        bounds = loc_bounds, scale_bounds, gain_bounds, lambda_bounds

        bins = np.arange( active_.min(), active_.max(), dE )
        active_fit = histogram_fits( active_, bins, func, p0 = p0, bounds = zip( *bounds ), tol=10 )
        loc.append( stats.uncertanty( active_fit.popt[0], active_fit.perr[0]) )
        sigma.append( stats.uncertanty( active_fit.popt[1], active_fit.perr[1]) )
        gain.append( stats.uncertanty( active_fit.popt[2], active_fit.perr[2]) )
        lamb.append( stats.uncertanty( active_fit.popt[3], active_fit.perr[3]) )
        print (loc[-1].val, sigma[-1].val, gain[-1].val, lamb[-1].val), (loc[-1].err(), sigma[-1].err(), gain[-1].err(), lamb[-1].err())
    print 'mean', np.mean( loc ), np.sum(loc)/len(loc), np.std( loc )
    print 'sigma', np.mean( sigma ), np.sum(sigma)/len(sigma), np.std( sigma )
    print 'gain', np.mean( gain ), np.std( gain )
    print 'lambda', np.mean( lamb ), np.std( lamb )
    
    exit(0)
    return (loc.val/divide, sigma.val/divide, gain.val/divide, lamb.val/divide), (loc.err()/divide, sigma.err()/divide, gain.err()/divide, lamb.err()/divide)

def compute_params( image ):
    active, overscan = separate_active_overscan( image )
    
    #select half of the active area
    active = active[:, active.shape[1]/2:]

    ##print stats
    #print type(active)
    #print 'shapes', active.shape, overscan.shape
    #print 'mean', np.mean( active ), np.mean( overscan )
    #print 'median', np.median( active ), np.median( overscan )
    #print 'std', np.std( active ), np.std( overscan )
    #print 'MAD', stats.MAD( active ), stats.MAD( overscan )
    
    #subtract overscan lines
    overscan_lines_median = np.median( overscan, axis = 1 )
    #print 'shape', overscan_lines_median.shape
    active = active - overscan_lines_median[:,None]
    overscan = overscan - overscan_lines_median[:,None]

    ##print stats
    #print type(active)
    #print 'shapes', active.shape, overscan.shape
    #print 'mean', np.mean( active ), np.mean( overscan )
    #print 'median', np.mean(np.median( active, axis=1 )), np.mean(np.median( overscan, axis=1 ))
    #print 'std', np.std( active ), np.std( overscan )
    #print 'MAD', stats.MAD( active ), stats.MAD( overscan )
    
    omu = np.mean( np.median( overscan, axis=1 ) )
    osigma = 1.4826 * np.mean(np.median( np.abs(overscan - np.median(overscan, axis=1)[:,None]), axis=1 ))

    #omu_errSqr = osigma/len()
    #osigma_errSqr = stats.errSqr_MAD( overscan )
    
    #print 'omu', omu
    #print 'osig', osigma

    amu = np.mean(np.median(active, axis=1))
    asigma = 1.4826 * np.mean(np.median( np.abs(active - np.median(active, axis=1)[:,None]), axis=1 ))
    
    #print 'amu', amu
    #print 'asig', asigma
    
    g = (asigma**2 - osigma**2)/(amu - omu)
    lamb = (amu - omu)/g

    #print 'g', g
    #print 'lamb', lamb

    func = lambda x, loc_, scale_, gain_, mu_: stats.poisson_norm.pdf( x, loc_, scale_, gain_, mu_, fix_loc=False )
    p0 = ( omu, osigma, g, lamb )

    dE = 1.
    bins = np.arange( active.min(), active.max(), dE )

    loc_bounds = (omu - osigma, omu + osigma)
    scale_bounds = ( (1-.1)*osigma, (1+.1)*osigma )
    gain_bounds = ( (1-.1)*g, (1+.1)*g )
    lambda_bounds = ( (1-.1)*lamb, (1+.1)*lamb )
    bounds = loc_bounds, scale_bounds, gain_bounds, lambda_bounds
    #print 'bounds', bounds
    
    active_fit = histogram_fits( active, bins, func, bounds=zip(*bounds), p0 = p0, tol=10 )
    loc = stats.ufloat( active_fit.popt[0], active_fit.perr[0]*active_fit.chisq )
    sigma = stats.ufloat( active_fit.popt[1], active_fit.perr[1]*active_fit.chisq )
    gain = stats.ufloat( active_fit.popt[2], active_fit.perr[2]*active_fit.chisq )
    lamb = stats.ufloat( active_fit.popt[3], active_fit.perr[3]*active_fit.chisq )
    #print 'loc', loc, active_fit.popt[0], active_fit.perr[0]
    #print 'sigma', sigma, active_fit.popt[1], active_fit.perr[1]
    #print 'g', active_fit.popt[2], active_fit.perr[2]
    #print 'lamb', active_fit.popt[3], active_fit.perr[3]

    return ((active_fit.popt[0], active_fit.perr[0]),
        (active_fit.popt[1], active_fit.perr[1]),
        (active_fit.popt[2], active_fit.perr[2]),
        (active_fit.popt[3], active_fit.perr[3])
        )
        
def table_of_fits( image ):
    print image_MLE( image )
    exit(0)
    
    active, overscan = get_regions( image )
    print 'fits shape', image.shape

    s = ' +- '
    
    dE = 1.
    bins = np.arange( active.min(), active.max(), dE )

    overscan_mean = np.mean(overscan)
    overscan_centered = overscan - overscan_mean
    overscan_std = np.std(overscan)
    overscan_var = np.var(overscan)
    print 'MLE mu', overscan_mean, overscan_std/np.sqrt(len(overscan)), 'sigma', overscan_std, np.sqrt(np.mean( overscan_centered**4 ))/( overscan_std*np.sqrt(len(overscan)) )
    overscan_fit = histogram_fits( overscan, bins, lambda x, loc, scale: stats.norm.pdf( x, loc, scale ), p0 = (overscan_mean, overscan_std), bounds = zip( (-np.inf,np.inf), (0,np.inf) ) )
    print 'fit', 'mu[ADU]', print_value_error( overscan_fit.popt[0], overscan_fit.perr[0], s=s )
    print 'fit', 'sigma[ADU]', print_value_error( overscan_fit.popt[1], overscan_fit.perr[1], s=s )
    print 'chisq', overscan_fit.chisq
    
    active_mean = np.mean(active)
    active_var = np.var(active)
    active_gain = ( active_var - overscan_var )/( active_mean - overscan_mean )
    active_lambda = ( active_mean - overscan_mean )**2/( active_var - overscan_var )
    print 'MLE g', active_gain, 'lambda', active_lambda

    func = lambda x, loc, scale, gain, mu: stats.poisson_norm.pdf( x, loc, scale, gain, mu, fix_loc=False )
    p0 = ( overscan_mean, overscan_std, active_gain, active_lambda )
    
    delta = .1
    amplitude_bounds = (0,np.inf)
    loc_bounds = (-np.inf,np.inf)

    scale_bounds = (0, np.inf)
    gain_bounds = (0, np.inf)
    lambda_bounds = (0,np.inf)
    
    bounds = loc_bounds, scale_bounds, gain_bounds, lambda_bounds
    active_fit = histogram_fits( active, bins, func, p0 = p0, bounds = zip( *bounds ), tol=10 )
    print 'fit', 'mu[ADU]', print_value_error( active_fit.popt[0], active_fit.perr[0], s=s )
    print 'fit', 'sigma[ADU]', print_value_error( active_fit.popt[1], active_fit.perr[1], s=s )
    print 'fit', 'sigma[e]', print_value_error( active_fit.popt[1]/active_fit.popt[2], active_fit.perr[1]/active_fit.popt[2], s=s )
    print 'fit', 'gain[ADU/e]', print_value_error( active_fit.popt[2], active_fit.perr[2], s=s )
    print 'fit', 'gain[ADU/keV]', print_value_error( active_fit.popt[2]/eIonization*1e3, active_fit.perr[2]/eIonization*1e3, s=s )
    print 'fit', 'lambda', print_value_error( active_fit.popt[3], active_fit.perr[3], s=s )
    print 'chisq', active_fit.chisq

    func = lambda x, loc, gain, mu: stats.poisson_norm.pdf( x, loc+gain*mu, overscan_std, gain, mu, fix_loc=False )
    p0 = ( overscan_mean, active_gain, active_lambda )

    bounds = loc_bounds, gain_bounds, lambda_bounds
    active_fit = histogram_fits( active, bins, func, p0 = p0, bounds = zip( *bounds ), tol=10 )
    print 'fit', 'mu[ADU]', print_value_error( active_fit.popt[0], active_fit.perr[0], s=s )
    print 'fit', 'gain[ADU/e]', print_value_error( active_fit.popt[1], active_fit.perr[1], s=s )
    print 'fit', 'gain[ADU/keV]', print_value_error( active_fit.popt[1]/eIonization*1e3, active_fit.perr[1]/eIonization*1e3, s=s )
    print 'fit', 'lambda', print_value_error( active_fit.popt[2], active_fit.perr[2], s=s )
    print 'chisq', active_fit.chisq
    
    exit(0)
    
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

    yccd = data['scn'].shape[1]
    xccd = data['scn'].shape[0]
    if yccd in [4262,4350]:
        yos = yccd-150
    elif yccd == 4662:
        yos = yccd-550
    else:
        print 'shape not predicted, edit this line of the code to include it'
        print 'from the FITS header, BIASSEC =', hdulists[0][0].header['BIASSEC'] 
        exit(0)
    print 'overscan width', yccd-yos
    print 'shape', yccd, xccd

    shift = lambda bins: ( bins[1:] + bins[:-1] )*.5

    maxSCN = max(clean(data['scn']))
    bins = np.linspace(-maxSCN, maxSCN, N)
    dx = 2*maxSCN/N
    if verbose: print 'dx', dx
    if verbose: print 'bins', min(bins), max(bins)

    if verbose: print 'computing histogram of full ccd'
    #vslice = slice(120,4180)
    vslice = slice(120,xccd-70-50)
    hslices = odict( [
        #('os', slice(4112+10, yccd-10)), 
        #('ac', slice(20,4100)), 
        ('os', slice(yos+10, yccd-10)), 
        ('ac', slice(20,yos-20)), 
        #('ac_os', slice(3112+10,3262-10))
        ] )
    countCut = 1e1/dx
    parseHist = lambda h,x: { 'y': h[h/dx>countCut]/dx, 'x': (.5*(x[1:]+x[:-1]))[h/dx>countCut] }
    #parseHist = lambda h,x: { 'y': h, 'x': (.5*(x[1:]+x[:-1])) }
    
    hist = odict( [ (datakey, { 
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
            print 'pixelcount', datum['pixelcount'], 'total', datum['totalpixel'], 'occupancy', 1.-float(datum['pixelcount'])/datum['totalpixel']
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
                line += [ float(fits_image.header['TEMPMIN' ]) ]
                if first: 
                    header += ['tempmax']
                    fmt += ['%.2f']
                line += [ float(fits_image.header['TEMPMAX' ]) ]
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
    line2 = [ int(runID), int(ohdu), float(fits_image.header['TEMPMIN' ]), float(fits_image.header['TEMPMAX']), datum['pixelcount'] ]
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

def processImage3( image, threshold, structure='cross' ):
    image[ image == 1e10 ] = 0

    t0 = time.clock()
    if structure == 'cross':
        s33 = [[0,1,0], [1,1,1], [0,1,0]]
    elif structure == 'square':
        s33 = [[1,1,1], [1,1,1], [1,1,1]]
    else:
        print 'unpredicted structure'
        exit(0)
    
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

def processImage4( image, threshold, structure='cross' ):
    image[ image == 1e10 ] = 0

    t0 = time.clock()
    if structure == 'cross':
        s33 = [[0,1,0], [1,1,1], [0,1,0]]
    elif structure == 'square':
        s33 = [[1,1,1], [1,1,1], [1,1,1]]
    else:
        print 'unpredicted structure'
        exit(0)

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

#class utils:
    #@staticmethod
    #def try_convertion( s, t ):
        #result = None
        #try: result = t(s)
        #except: pass
        #return result
    
class SpectrumWindow(Gtk.Window):

    #def __del__(self):
        #os.remove(self.tempPath)
        #Gtk.Window.__del__(self)
        
    def __init__(self):
        self.id = time.time()
        print 'id', self.id
        self.tempPath = '/home/mota/public/gui/session_spectrum_%s.png'%self.id
        Gtk.Window.__init__( self, title="Spectrum Viewer [session%s]"%self.id )
        self.set_border_width(3)

        self.paths = []
        self.pathslabel = Gtk.Label()
        self.pathslabel.set_label( 'None' )

        drawspectrumbutton = Gtk.Button(use_underline=True)
        drawspectrumbutton.set_label('_spectrum')
        drawspectrumbutton.connect( 'clicked', self.on_spectrum_click )
        
        runsLabel = Gtk.Label()
        runsLabel.set_label('runs')
        
        listRuns = Gtk.VBox()
        for run in list_subruns()[::-1]:
            runbutton = Gtk.ToggleButton()
            runbutton.set_label( run )
            runbutton.connect( 'toggled', self.ontogglebuttonruns, run )
            listRuns.pack_start( runbutton, False, False, 0 )

        scrolledwindow = Gtk.ScrolledWindow()
        scrolledwindow.add_with_viewport( listRuns )
        scrolledwindow.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        
        runIDPanel = Gtk.VBox()
        runIDPanel.pack_start( runsLabel, False, False, 1 )
        runIDPanel.pack_start( scrolledwindow, True, True, 0 )
        
        selectionLabel = Gtk.Label()
        selectionLabel.set_label( 'selection' )
        self.selectionEntry = Gtk.Entry()
        self.selectionEntry.set_text( 'ohdu==3 && E0<22e3 && E0>700' )
        
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
        self.fitFunctionEntry.set_text( 'lambda x, a, b, c: np.log( a/(x**c) + b )' )
        
        self.imageBox = Gtk.VBox()
        #self.image = Gtk.Image()
        #self.imageBox.pack_start(self.image, True,True,0)
        canvas, toolbar = self.build_figure()
        self.imageBox.pack_start(canvas, True, True, 0)
        #self.imageBox.pack_start(toolbar, False, False, 0)
        
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
        
        self.outputPanel = Gtk.Label()
        self.outputPanel.set_label('outputPanel')
        self.outputPanel.set_justify(Gtk.Justification.LEFT)
        box = Gtk.HBox()
        box.pack_start(self.outputPanel, False, False, 0)
 
        imagePanel = Gtk.VBox()
        scroll = Gtk.ScrolledWindow()
        scroll.add_with_viewport(self.imageBox)
        imagePanel.pack_start( optionsPanel, False, False, 0 )
        imagePanel.pack_start( box, False, True, 0 )
        imagePanel.pack_start( scroll, True, True, 0 )
        
        body = Gtk.HBox()
        body.pack_start( runIDPanel, False, False, 0 )
        body.pack_end( imagePanel, True, True, 1 )

        footer = Gtk.HBox()
        footer.pack_end( drawspectrumbutton, False, False, 1 )
        footer.pack_start( self.pathslabel, True, True, 1 )
        
        document = Gtk.VBox()
        document.pack_start( body, True, True, 3 )
        document.pack_start( footer, False, False, 3 )
        
        #self.append_outputPanel('update')
        self.add( document )
        self.maximize()
        self.connect( 'destroy', Gtk.main_quit )
        #self.connect( 'key_press_event', self.onkeypressevent )
        self.show_all()

    #def update_outputPanel(self):
        #text = []
        #text.append( 'catalogs:\n' + '\n'.join( self.paths ) )
        #self.outputPanel.set_label( '\n'.join(text) )
        #self.show_all()
        
    def append_outputPanel(self, text):
        self.outputPanel.set_text( '%s\n%s'%(self.outputPanel.get_label(), text) )
        while Gtk.events_pending():
            Gtk.main_iteration()
        
    def ontogglebuttonruns( self, button, run ):
        catalogs = list_catalogs_skim1( run )
        if catalogs == []:
            print 'no runIDs for run', run
            button.set_active(False)
            return
        if button.get_active() == True:
            self.paths.append( catalogs[0] )
        else:
            self.paths.remove( catalogs[0] )
        print 'paths', self.paths
        self.pathslabel.set_label( '\n'.join( [ '[%d] %s'%(i,path) for i, path in enumerate(self.paths) ] ) )
        self.show_all()
    
    def build_figure(self):
        self.fig = Figure()
        #self.fig.suptitle( self.selectionEntry.get_text() )
        grid = plt.GridSpec( 4, 1, hspace=0, wspace=0 )
        self.ax_main = self.fig.add_subplot( grid[:-1, 0])
        self.ax_ratio = self.fig.add_subplot( grid[-1,0], sharex=self.ax_main)
        canvas = FigureCanvas(self.fig)
        
        #self.imageBox.pack_start( image, True, True, 0 )
        #toolbar = NavigationToolbar(canvas, self)
        toolbar = None
        return canvas, toolbar 
        #self.imageBox.pack_start(toolbar, False, True, 0)

        
    def on_spectrum_click( self, button ):
        self.append_outputPanel('generating spectrum...')
        #fig = plt.figure()
        #fig = Figure()
        if self.selectionEntry.get_text() == '':
            self.append_outputPanel('Error: selection not specified')
            return
        self.fig.suptitle( self.selectionEntry.get_text() )
        self.ax_main.cla()
        self.ax_ratio.cla()
        #grid = plt.GridSpec( 4, 1, hspace=0, wspace=0 )
        #ax_main = fig.add_subplot( grid[:-1, 0])
        #ax_ratio = fig.add_subplot( grid[-1,0], sharex=ax_main)
        
        binMin = None
        binMax = None
        bins = None
        #data = []
        #nRunIDs = []
        for i, path in enumerate(self.paths):
            self.append_outputPanel( 'reading catalog %s'%path )
            datum = root_numpy.root2array( path, treename='hitSumm', branches=self.expressionEntry.get_text(), selection=self.selectionEntry.get_text() )
            self.append_outputPanel( 'read %d events'%len(datum) )
            #data.append( datum )
            #nRunIDs.append( len( list_runIDs_from_subrun( getrunFromPath(path) ) ) )
            nRunID = len( list_runIDs_from_subrun( getrunFromPath(path) ) )
            if binMin is None: binMin = np.min(datum)
            else: binMin = min(binMin, np.min( datum ) )
            if binMax is None: binMax = np.max(datum)
            else: binMax = max(binMax, np.max( datum ) )
            if bins is None:
                if self.binSizeEntry.get_text() == 'auto': 
                    bins = np.linspace( binMin, binMax, 100 )
                else: 
                    bins = np.arange( binMin, binMax, float(self.binSizeEntry.get_text()) )
                self.append_outputPanel( 'bins %f %f %f %d'%(binMin, binMax, bins[1]-bins[0], len(bins)) )
        #for i, (datum, path, nRunID) in enumerate( zip( data, self.paths, nRunIDs ) ):
            hist = np.histogram( datum, bins=bins )[0]
            hist = hist.astype(float)/nRunID
            one_count = 1./nRunID
            x = .5*(bins[:-1]+bins[1:])
            median = np.median( hist )
            funcstr = self.fitFunctionEntry.get_text()
            varnames = map( lambda s:s.strip(), funcstr.split(':')[0].split(',')[1:] )
            try:
                func = eval( funcstr )
                p = scipy.optimize.curve_fit( func, x[hist>one_count], np.log( hist[hist>one_count] ) )[0]
                params_str = ' '.join( [ '%s:%.2g'%(v,pp) for v,pp in zip(varnames, p) ] )
                self.append_outputPanel( 'fit successful with params \n%s'%params_str )
                y = np.exp( func( x, *p ) )
                label = '%s [%d] '%(getrunFromPath(path),nRunID) + params_str
                self.ax_main.plot( x, y, color='C%d'%i, linestyle='--' )
                dy = hist-y
                dy[dy<0] = 0
                self.ax_ratio.step( x, dy, where='mid', color='C%d'%i )
                #try:
                    #peakFunc = lambda x, a, b, c, d, e: b*scipy.stats.norm.pdf(x,loc=abs(a)*1.74+c,scale=15*d) + b*scipy.stats.norm.pdf(x,loc=abs(a)*8.048+c,scale=15*d)+ b*scipy.stats.norm.pdf(x,loc=abs(a)*8.905+c,scale=15*d) + e
                    #pg0 = [1700.,10,10,10]
                    #pg = scipy.optimize.curve_fit( peakFunc, x[hist>one_count], hist[hist>one_count]/y, p0=pg0 )[0]
                    #print 'pg', pg
                    #y = peakFunc(x,*pg)
                    #ax_ratio.plot( x[y>0], y[y>0], color='C%d'%i, linestyle=':' )
                #except:
                    #pass
            except:
                label = '%s [%d] '%(getrunFromPath(path),nRunID) + 'fit failed'
                self.append_outputPanel( 'fit failed' )
            self.ax_main.step( x, hist, label=label, where='mid', color='C%d'%i )
        self.ax_main.set_yscale( 'log' )
        self.ax_ratio.set_xlabel( self.expressionEntry.get_text() )
        #ax_ratio.set_yscale( 'log' )
        self.ax_main.legend( fancybox=True, framealpha=0 )

        #image = FigureCanvas(fig)
        #image.set_size_request(400,400)
        #children = self.imageBox.get_children()
        #children[-1].destroy()
        #children[-2].destroy()
        
        #self.imageBox.pack_start( image, True, True, 0 )
        #toolbar = NavigationToolbar(image, self)
        #self.imageBox.pack_start(toolbar, False, True, 0)
        self.fig.canvas.draw()
        #fig.savefig( self.tempPath )
        #plt.close()
        self.outputPanel.set_text('plot generated successfully')
        #self.image.set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.tempPath ) )
        self.show_all()

def computeEqs( X, labels ):
    cov = np.cov( X )
    eigenW, eigenV = np.linalg.eigh( cov )
    observations_in_newAxes = np.dot( eigenV, X )
    print 'eigenV', eigenV.shape
    print 'newAxes', observations_in_newAxes.shape
    
    #print len(newAxes), len(newAxes[0])
    axesOrdered = np.argsort(np.abs(eigenW))
    #print axesOrdered
    #print eigenW[0], eigenW[axesOrdered[0]]
    #print eigenW[-1], eigenW[axesOrdered[-1]]
    #print np.abs(eigenW).min(), np.abs(eigenW).max()
    #newAxes = []
    #newObservations = []
    #for n in range(axesOrdered):
        #axis = axesOrdered[n]
        #newAxes = [ eigenV[:,axis] ]
        #observationsNewAxes = [ obserservation_in_newAxes[] ]
 
    coefs = []
    variables = []
    stds = []
    eqlabels = []
    reordered_observations = []
    for n in range(len(X)):
        axis = axesOrdered[n]
        reordered_observations.append( observations_in_newAxes[axis] )
        eqlabels.append(labels[axis])
        #c = eigenV[:,axis]/np.abs(eigenV[:,axis]).max()
        args = np.argsort( eigenV[:,axis] )[::-1]
        eigenV[:,axis] /= eigenV[args[-1],axis]
        print 'args', args
        coefs.append( [ eigenV[arg,axis] for arg in args ] ) # if abs(eigenV[arg,axis]) > 1e-3
        variables.append( [ labels[arg] for arg in args ] ) # if abs(eigenV[arg,axis]) > 1e-3
        stds.append( eigenW[axis] )
    print 'coefs', coefs
    print variables
    print eqlabels
    print stds

    fig = plt.figure()
    w, h = fig.get_size_inches()
    fig.set_size_inches(4*w,4*h)
    grid = plt.GridSpec(cov.shape[0], cov.shape[1])
    for i in range(cov.shape[0]):
        for j in range(cov.shape[1]):
            if i>j:
                ax = fig.add_subplot( grid[i,j] )
                ax.scatter( X[i], X[j], s=1. )
                ax.set_xlabel( labels[i] )
                ax.set_ylabel( labels[j] )
            elif i<j:
                ax = fig.add_subplot( grid[i,j] )
                ax.scatter( reordered_observations[i], reordered_observations[j], s=1. )
                ax.set_xlabel( ''.join( [ '%+.2f*%s'%(coefs[i][k], variables[i][k]) for k in range(len(coefs[i])) ] ) )
                ax.set_ylabel( ''.join( [ '%+.2f*%s'%(coefs[j][k], variables[j][k]) for k in range(len(coefs[j])) ] ) )
    fig.savefig('computeEqs.png')

    return (coefs, variables, eqlabels, stds)#, observations_in_newAxes, newAxes, axisDispersions

def print_value_error( x, dx, s='\pm' ):
    power = np.log10(abs(dx))
    power_int = int(np.floor(abs(power)))+1
    #print 'power', power, power_int, x, dx
    if power > 0:
        power_int -= 1
        if power_int == 0:
            pre_str = '%d\pm%d'%(round(x), round(dx))
        else:
            pre_str = '(%d%s%d)10^{%d}'%(round(x/10**(power_int)), s, round(dx/10**(power_int)), power_int)
    else:
        pre_str = '%%.%df%s%%.%df'%(power_int,s,power_int)%(x,dx)
    #print pre_str
    return pre_str

def computeParamsFromFitsFile( filename ):
    fitsfile = astropy.io.fits.open( filename )
    header = ['ohdu', 'mu', 'muE', 'sigma', 'sigmaE', 'gain', 'gainE', 'lambda', 'lambdaE']
    table = []
    for f in fitsfile:
        ohdu = int(f.header['OHDU'])
        if ohdu in [2,3,4,5,6,7,8,9,10,13,14]:
            fitsimage = f.data
            res = compute_params( fitsimage )
            entry = [ohdu] + np.array(res).flatten().tolist()
            table.append(entry)
    return header, table

def plotParamsFromCSVList( filenames, outname ):
    #ax = [fig.add_subplot(111)]
    fulltable = []
    for filename in filenames:
        table = np.genfromtxt(filename, delimiter=',')
        fulltable = table
    #fulltable = np.array( fulltable )
    print fulltable.shape
    for ohdus in [ [2,3,4,5], [6,7,8,9], [10,13,14] ]:
        fig = plt.figure()
        sigma_ax = fig.add_subplot(3,1,1)
        gain_ax = fig.add_subplot(3,1,2)
        lambda_ax = fig.add_subplot(3,1,3)
        sigma_ax.set_ylabel('sigma[ADU]')
        gain_ax.set_ylabel('gain[ADU/e-]')
        lambda_ax.set_ylabel('lambda')
        for ohdu in ohdus:
            mask = np.all( [fulltable[:,1] == ohdu], axis=0 )
            sigma_ax.errorbar( fulltable[:,0][mask], fulltable[:,4][mask], yerr=fulltable[:,5][mask], fmt='.', label='%s'%ohdu )
            gain_ax.errorbar( fulltable[:,0][mask], fulltable[:,6][mask], yerr=fulltable[:,7][mask], fmt='.' )
            lambda_ax.errorbar( fulltable[:,0][mask], fulltable[:,8][mask], yerr=fulltable[:,9][mask], fmt='.' )
        sigma_ax.legend( fancybox=True, framealpha=0 )
        name = outname.split('.')
        fig.savefig( name[0] + '_%02d.' % ohdus[0] + name[1] )
        plt.close()
    return

class NormRVS:
    def __init__( self, mu, sigma, N=1000 ):
        self.mu = mu
        self.sigma = sigma
        self.N = N
        self.x = scipy.stats.norm.rvs( self.mu[0], self.sigma, self.N )
        self.y = scipy.stats.norm.rvs( self.mu[1], self.sigma, self.N )
        self.pos = np.array( zip(self.x, self.y) )
        self.mean = np.mean( self.pos, axis=0 )
        self.var = np.var( self.pos, axis=0 )
        self.std = np.sqrt( self.var )
        #self.fitted = NormFitted( self.x, self.y, self.mean, np.sqrt(self.var) )
        self.pixelated = Pixelated( self.x, self.y )
        self.pixelated5 = Pixelated( self.x, self.y, rebin=(1,5) )

    def plot( self, ax, **kwargs ):
        ax.grid(True)
        ax.set_aspect(aspect=1)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.scatter( self.x, self.y, s=.1, c='b' )
        ax.errorbar( self.mean[0], self.mean[1], xerr=self.std[0], yerr=self.std[1], c='r' )
        ax.set_title(r'$\mu=(%s,%s)\quad \sigma=%s$' %( self.mu[0], self.mu[1], self.sigma ) )
        
        return ax
        
class NormFitted:
    def __init__( self, x, y, mean, std ):
        self.x = x
        self.y = y
        self.mu = [0,0]
        self.sigma = np.sqrt( .5*(std[0]**2 +std[1]**2) )

        self.mu[0], self.mu[1], self.sigma = scipy.optimize.minimize( 
            fun = self.negloglikelihood,
            x0 = [ mean[0], mean[1], std[0] ]
            ).x
        
    def pdf( self, x, y, p ):
        return scipy.stats.norm.pdf( x, p[0], p[2] )*scipy.stats.norm.pdf( y, p[1], p[2] )

    def negloglikelihood( self, p ):
        pdf = self.pdf( self.x, self.y, p )
        return -np.nansum( np.log( pdf ) )
    
    def plot(self, ax, show_original=False ):
        ax.grid(True)
        ax.set_aspect(aspect=1)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        X, Y = np.meshgrid( np.linspace( self.x.min(), self.x.max(), 100), np.linspace( self.y.min(), self.y.max(), 100) )
        z = self.pdf( X, Y, (self.mu[0], self.mu[1], self.sigma) )
        ax.contourf(X, Y, z, cmap=cm.PuBu)
        ax.errorbar( self.mu[0], self.mu[1], xerr=self.sigma, yerr=self.sigma, c='r' )
        ax.set_title(r'$\mu=(%.3f,%.3f)\quad\sigma=%.3f$' %( self.mu[0], self.mu[1], self.sigma ) )
        if show_original:
            ax.scatter(self.x, self.y, s=.1, c='k' )
        return ax
            
        
class Pixelated:
    def __init__( self, x, y, rebin=(1,1) ):
        self.x, self.y = x/rebin[0], y/rebin[1]
        self.xpix = np.floor(self.x).astype(int)
        self.ypix = np.floor(self.y).astype(int)
        
        self.pix, self.pix_count = np.unique( zip(self.xpix, self.ypix), axis=0, return_counts=True )
        self.pixs_x, self.pixs_y = self.pix.T
        self.n_y = self.pixs_y.max() - self.pixs_y.min() + 1
        xbins = np.arange( self.xpix.min(), self.xpix.max()+1.1, 1)
        ybins = np.arange( self.ypix.min(), self.ypix.max()+1.1, 1)
        self.hist, self.xedges, self.yedges = np.histogram2d( self.xpix, self.ypix, bins=[xbins,ybins] )

        self.meanxpix, self.meanypix = np.average( self.pix+.5, weights=self.pix_count, axis=0 )
        self.stdx, self.stdy = np.sqrt( self.var_weighted( self.pix+.5, self.pix_count ) )
        self.avestd = np.mean([self.stdx, self.stdy])
        
        self.fitted = PixelatedNormFitted( self.pixs_x, self.pixs_y, self.pix_count, (self.meanxpix, self.meanypix), (self.avestd, self.avestd), len(self.xpix), fixed_aspect=True )
        self.fittedx = PixelatedNormFitted( self.pixs_x, self.pixs_y, self.pix_count, (self.meanxpix, self.meanypix), (self.stdx, self.stdy), len(self.xpix) )

    def var_weighted( self, x, w ):
        return np.average( x**2, weights=w, axis=0 ) - np.average( x, weights=w, axis=0 )**2

    def plot(self, ax, show_original = False, **kwargs):
        ax.grid(True)
        ax.set_aspect(aspect=1)
        ax.set_xlabel(r'$x[{\rm pix}]$')
        ax.set_ylabel(r'$y[{\rm pix}]$')
        z = self.pix_count.astype(float)/self.pix_count.max()

        for xi, yi, zi in zip( self.pixs_x, self.pixs_y, z):
            rectangle = plt.Rectangle( (xi,yi), width=1, height=1, color=(0,0,1,np.sqrt(zi) ), linewidth=0 )
            ax.add_artist(rectangle)
        
        ax.errorbar( self.meanxpix, self.meanypix, xerr=self.std, yerr=self.std, c='r' )
        ax.set_xlim((self.pixs_x.min(), self.pixs_x.max()+1))
        ax.set_ylim((self.pixs_y.min(), self.pixs_y.max()+1))
        
        if show_original:
            ax.scatter(self.x, self.y, s=.1, c='k' )
        ax.set_title(r'${\rm mean}=(%.3f,%.3f)\quad{\rm std}=(%.3f,%.3f)$' %( self.meanxpix, self.meanypix, self.stdx, self,stdy ) )
        return ax

    def plot_hist(self, ax, show_original = False, **kwargs):
        ax.grid(True)
        ax.set_aspect(aspect=1)
        ax.set_xlabel(r'$x[{\rm pix}]$')
        ax.set_ylabel(r'$y[{\rm pix}]$')
        H = self.hist.T/float(self.hist.max())
        ax.imshow( 
            np.sqrt(H), 
            aspect='equal', 
            alpha=.5, 
            origin='lower', 
            extent=(self.xedges.min(), self.xedges.max(), self.yedges.min(), self.yedges.max()),
            cmap='Blues'
            )
        ax.set_xticks(self.xedges)
        ax.set_yticks(self.yedges)
        if show_original: ax.scatter(self.x, self.y, s=.1, c='k' )
        ax.set_title(r'${\rm mean}=(%.3f,%.3f)\quad{\rm std}=(%.3f,%.3f)$' %( self.meanxpix, self.meanypix, self.stdx, self.stdy ) )
        return #ax
        
class PixelatedNormFitted:
    def __init__( self, pixs_x, pixs_y, pix_count, mean, std, N, fixed_aspect = False ):
        self.fixed_aspect = fixed_aspect
        self.pixs_x = pixs_x
        self.pixs_y = pixs_y
        self.pix_count = pix_count
        self.mean = mean
        self.std = std
        self.N = N
        self.mu = [0,0]
        self.sigma = [0,0]
        self.mu[0], self.mu[1], self.sigma[0], self.sigma[1] = scipy.optimize.minimize(
            fun = self.negloglikelihood,
            x0 = [ self.mean[0], self.mean[1], self.std[0], self.std[1] ],
            bounds = [(-np.inf,np.inf), (-np.inf, np.inf), (0, np.inf), (0, np.inf)] ).x
        
    def integral_norm_pixel( self, xpix, mu, sigma ):
        sqrt2 = np.sqrt(2.)
        return -.5*( scipy.special.erf( -( xpix+1 - mu)/( sqrt2*sigma )) - scipy.special.erf( -( xpix - mu )/( sqrt2*sigma ) ) )
    
    def probability_pixel( self, x, y, p ):
        if self.fixed_aspect: p[3] = p[2]
        return self.integral_norm_pixel( x, mu = p[0], sigma = p[2] ) * self.integral_norm_pixel( y, mu = p[1], sigma = p[3] )
    
    def pdf( self, p ):
        return scipy.stats.binom.pmf( self.pix_count, self.N, self.probability_pixel( self.pixs_x, self.pixs_y, p ))
    
    def negloglikelihood( self, p ):
        pdf = self.pdf( p )
        return -np.nansum( np.log( pdf ) )

    def plot( self, ax, show_original = False, **kwargs ):
        ax.grid(True)
        ax.set_aspect(aspect=1)
        ax.set_xlabel(r'$x[{\rm pix}]$')
        ax.set_ylabel(r'$y[{\rm pix}]$')
        
        X, Y = np.meshgrid(np.arange( self.pixs_x.min(), self.pixs_x.max()+1, 1), np.arange( self.pixs_y.min(), self.pixs_y.max()+1, 1) )
        X = X.flatten()
        Y = Y.flatten()

        Z = self.probability_pixel( X, Y, (self.mu[0], self.mu[1], self.sigma) )
        z = Z.astype(float)/Z.max()
        
        for xi, yi, zi in zip( X, Y, z ):
            rectangle = plt.Rectangle( (xi,yi), width=1, height=1, color=(0,0,1,np.sqrt(zi) ), linewidth=0 )
            ax.add_artist(rectangle)
        ax.set_xlim((self.pixs_x.min(), self.pixs_x.max()+1))
        ax.set_ylim((self.pixs_y.min(), self.pixs_y.max()+1))
        
        ax.errorbar( self.mu[0], self.mu[1], xerr=self.sigma, yerr=self.sigma, c='r' )

        ax.set_title(r'$\mu=(%.3f,%.3f)\quad\sigma=%.3f$' %( self.mu[0], self.mu[1], self.sigma ) )
        return ax

#class NoiseGenerator:
    #def __init__( self, sigma ):
        
    #def generateNoise( self, pixs_x, pixs_y ):
        #X, Y = np.meshgrid(np.arange( pixs_x.min(), pixs_x.max()+1, 1), np.arange( pixs_y.min(), pixs_y.max()+1, 1) )
        

class RandomSamplingFit:
    
    def __init__(self, N_min = 5, N_max = 2000, sigma_max = 1.2):
        self.N_min = N_min
        self.N_max = N_max
        self.sigma_max = sigma_max
        
    def generate_events(self, number_events):
        self.values = {}
        self.values[r'$N_{\rm e}$'] = np.random.randint( self.N_min, self.N_max, size=number_events)
        self.values[r'$\mu$'] = np.random.random( size=(number_events, 2) )
        self.values[r'$\sigma$'] = np.random.random( size=number_events ) * self.sigma_max
        
        self.norm_rvs = [ NormRVS( mu, sigma, N ) for mu, sigma, N in zip(self.values[r'$\mu$'], self.values[r'$\sigma$'], self.values[r'$N_{\rm e}$']) ]
        
        self.values[r'$N_{y\rm pix}$'] = [ event.pixelated.n_y for event in self.norm_rvs ]
        self.values[r'$\sigma_{\rm pix}$'] = [ event.pixelated.fitted.sigma[0] for event in self.norm_rvs ]
        self.values[r'$\sigma_{x\rm pix}$'] = [ event.pixelated.fittedx.sigma[0] for event in self.norm_rvs ]
        self.values[r'$\sigma_{y\rm pix}$'] = [ event.pixelated.fittedx.sigma[1] for event in self.norm_rvs ]
        self.values[r'$\sigma_{\rm pix5}$'] = [ event.pixelated5.fitted.sigma[0] for event in self.norm_rvs ]
        self.values[r'$\sigma_{x\rm pix5}$'] = [ event.pixelated5.fittedx.sigma[0] for event in self.norm_rvs ]
        self.values[r'$\sigma_{y\rm pix5}$'] = [ event.pixelated5.fittedx.sigma[1] for event in self.norm_rvs ]
        self.values[r'$N_{y\rm pix5}$'] = [ event.pixelated5.n_y for event in self.norm_rvs ]
    
    reconstruction_file = 'reconstruction.csv'
    def make_table_reconstruction( self, labels ):
        table = np.array([ self.values[label] for label in labels ]).T
        print 'data.shape', table.shape
        np.savetxt( self.reconstruction_file, table, header = ', '.join(labels), delimiter=', ' )
        
    @classmethod
    def plot_reconstruction( cls, ax, bins = None, show_data = True, ylabels = None, xlabel = None ):
        data = np.genfromtxt( cls.reconstruction_file, names = True, delimiter=',', deletechars='', replace_space=' ' )
        labels = data.dtype.names
        print 'found labels', labels
        if ylabels:
            print 'ylabels', ylabels
            _labels = []
            if not xlabel is None:
                _labels += [ xlabel ]
            else:
                _labels += [ labels[0] ] 
            _labels += ylabels
            labels = _labels
        print 'labels', labels
        if len(labels) == 1:
            raise Exception('no labels to plot')
        ax.set_xlabel( labels[0] )
        ax.set_aspect( aspect = 1 )

        inds = None
        center_bins = None
        width_bins = None
        if not bins is None:
            center_bins = np.mean( [bins[:-1], bins[1:]], axis=0 )
            width_bins = (bins[:-1] - bins[1:])/2
            inds = np.digitize( data[labels[0]], bins )
        if show_data:
            for label in labels[1:]:
                y = data[label]
                if 'y' in label and '5' in label: y *= 5
                ax.plot( data[labels[0]], y, '.', ms=1, label = label )
        
        if not inds is None:
            for label in labels[1:]:
                y = data[label]
                if 'y' in label and '5' in label: y *= 5
                mean = [ np.mean( np.array( y )[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                std = [ np.std( np.array( y )[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                ax.errorbar( center_bins, mean, xerr = width_bins, yerr = std, fmt='.', label = label )

        line = matplotlib.lines.Line2D( ax.get_xlim(), ax.get_xlim(), color='black', linestyle=':')
        ax.add_line( line )                
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        return (legend, )

    @classmethod
    def plot_reconstruction_ny( cls, ax, bins = None, show_data = True, ylabels = None ):
        data = np.genfromtxt( cls.reconstruction_file, names = True, delimiter=',', deletechars='', replace_space=' ' )
        labels = data.dtype.names
        print data.dtype
        selections = [ 
            (r' $N_{y\rm pix5}>1$', data[r'$N_{y\rm pix5}$'] > 1),
            (r' $N_{y\rm pix5}=1$', data[r'$N_{y\rm pix5}$'] == 1),
            ]

        if ylabels:
            print 'ylabels', ylabels
            labels = [labels[0]] + [ label for label in labels[1:] if label in ylabels] 
        print 'labels', labels, data.dtype.names
        if len(labels) == 1:
            raise Exception('no labels to plot')
        ax.set_xlabel( labels[0] )
        ax.set_aspect( aspect = 1 )

        inds = None
        center_bins = None
        width_bins = None
        if show_data:
            for label in labels[1:]:
                for selection in selections:
                    mask = selection[1]
                    ax.plot( data[labels[0]][mask], data[label][mask], '.', ms=1, label = label + selection[0] )
        
        for selection in selections:
            mask = selection[1]
            if not bins is None:
                center_bins = np.mean( [bins[:-1], bins[1:]], axis=0 )
                width_bins = (bins[:-1] - bins[1:])/2
                inds = np.digitize( data[labels[0]][mask], bins )
                for label in labels[1:]:
                    mean = [ np.mean( np.array( data[label][mask] )[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                    std = [ np.std( np.array( data[label][mask] )[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                    ax.errorbar( center_bins, mean, xerr = width_bins, yerr = std, fmt='.', label = label + selection[0] )

        line = matplotlib.lines.Line2D( ax.get_xlim(), ax.get_xlim(), color='black', linestyle=':')
        ax.add_line( line )                
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        return (legend, )
        
    def plot( self, ax, xlabel, ylabel, y2label = False, equal = False, relative_error = False, hist = False, bins = None, efficiency = False, **kwargs ):
        ax.set_xlabel( xlabel )
        x = np.array( self.values[xlabel] )
        if not isinstance(ylabel, list):
            ylabel = [ylabel]

        if relative_error:
            y = self.values[ylabel]
            y2 = self.values[y2label]
            relerr = lambda a, b: 2.*(np.array(a)-np.array(b))/(np.array(a)+np.array(b))
            z = relerr(y,y2)
            
            ax.set_ylabel( r'errRel(%s, %s)'%(ylabel, y2label) )
            if hist:
                ax.hist2d( x, z, bins=100 )
            else:
                dx = .5
                x = np.floor(x/dx).astype(int)
            
            line = matplotlib.lines.Line2D( ax.get_xlim(), (0,0), color='black', linestyle=':')
            ax.add_line( line )
        elif efficiency:
            center_bins = np.mean( [bins[:-1], bins[1:]], axis=0 )
            width_bins = (bins[:-1] - bins[1:])/2
            inds = np.digitize( x, bins )
            for ylabel_ in ylabel:
                mean = [ np.mean( np.array(self.values[ylabel_])[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                std = [ np.std( np.array(self.values[ylabel_])[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                prob = [ scipy.stats.norm.pdf( 0, loc = mean[i]-center_bins[i], scale=std[i] ) for i,_ in enumerate(center_bins)  ]
                ax.errorbar( center_bins, prob, xerr = width_bins, fmt='.', label = ylabel_ )
        else:
            ax.set_aspect(aspect=1)
            inds = None
            center_bins = None
            width_bins = None
            if not bins is None:
                center_bins = np.mean( [bins[:-1], bins[1:]], axis=0 )
                width_bins = (bins[:-1] - bins[1:])/2
                inds = np.digitize( x, bins )
            for ylabel_ in ylabel:
                ax.plot( x, self.values[ylabel_], '.', ms=1, label = ylabel_ )
            if not inds is None:
                for ylabel_ in ylabel:
                    mean = [ np.mean( np.array(self.values[ylabel_])[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                    std = [ np.std( np.array(self.values[ylabel_])[ inds-1 == i ] ) for i,_ in enumerate(center_bins) ]
                    ax.errorbar( center_bins, mean, xerr = width_bins, yerr = std, fmt='.', label = ylabel_ )
            if equal:
                line = matplotlib.lines.Line2D( ax.get_xlim(), ax.get_xlim(), color='black', linestyle=':')
                ax.add_line( line )
                
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        return (legend, )

def plot2file( file, plotfun, **kwargs ):
    fig = plt.figure()
    extra_artists = None
    if not isinstance( plotfun, list ):
        extra_artists = plotfun( fig.add_subplot(111), **kwargs )
    fig.savefig(file, bbox_extra_artists = extra_artists, bbox_inches='tight')
    print 'saved plot', file
    return

class Callable:
    @staticmethod
    def size_like( **kwargs ):
        print 'size_like function'
        generate_events = False
        if generate_events:
            number_events = 10000
            randomSamplingFit = RandomSamplingFit()
            randomSamplingFit.generate_events( number_events )
            randomSamplingFit.make_table_reconstruction(
                labels = [r'$\sigma$', r'$N_{\rm e}$', r'$N_{y\rm pix}$', r'$\sigma_{\rm pix}$', r'$\sigma_{x\rm pix}$', r'$\sigma_{y\rm pix}$', r'$N_{y\rm pix5}$', r'$\sigma_{\rm pix5}$', r'$\sigma_{x\rm pix5}$', r'$\sigma_{y\rm pix5}$']
                )
        else:
            #event = NormRVS( (0,0), 1, 2000)
            #plot2file(
                #'pixelatedB.png',
                #event.pixelated.plot_hist,
                #show_original = True
                #)
            #plot2file(
                #'pixelatedB5.png',
                #event.pixelated5.plot_hist,
                #show_original = True
                #)
            #event = NormRVS( (0,0), .5, 500)
            #plot2file(
                #'pixelatedA.png',
                #event.pixelated.plot_hist,
                #show_original = True
                #)
            #plot2file(
                #'pixelatedA5.png',
                #event.pixelated5.plot_hist,
                #show_original = True
                #)
            #exit(0)
            sets = (
                #(r'$\sigma_{\rm pix}$', r'$\sigma_{x\rm pix}$'),
                #(r'$\sigma_{\rm pix5}$', r'$\sigma_{x\rm pix5}$'),
                #(r'$\sigma_{x\rm pix5}$'),
                (r'$\sigma_{x\rm pix}$', r'$\sigma_{y\rm pix5}$', r'$\sigma_{y\rm pix}$'),
                (r'$\sigma_{x\rm pix}$', r'$\sigma_{y\rm pix}$'),
                (r'$\sigma_{x\rm pix5}$', r'$\sigma_{y\rm pix5}$', r'$\sigma_{x\rm pix}$'),
                )
            for i, ylabels in enumerate(sets):
                print 'ylabels', ylabels
                #plot2file(
                    #'reconstructionNy%d.png'%i, 
                    #RandomSamplingFit.plot_reconstruction_ny,
                    #ylabels = ylabels
                    #)
                #plot2file(
                    #'reconstructionBinNy%d.png'%i, 
                    #RandomSamplingFit.plot_reconstruction_ny,
                    #ylabels = ylabels,
                    #show_data = False,
                    #bins = np.arange( 0, 1.21, .1 )
                    #)
                plot2file(
                    'reconstructionXY%d.png'%i, 
                    RandomSamplingFit.plot_reconstruction,
                    ylabels = ylabels[1:],
                    xlabel = ylabels[0],
                    )
                plot2file(
                    'reconstructionXYBin%d.png'%i, 
                    RandomSamplingFit.plot_reconstruction,
                    ylabels = ylabels[1:],
                    xlabel = ylabels[0],
                    show_data = False,
                    bins = np.arange( 0, 1.21, .1 )
                    )
        
        exit(0)
        
        return
    #@staticmethod
    #def reconstructionViewer( **kwargs ):
        #win = ReconstructionViewer( **kwargs )
        #Gtk.main()
    
    @staticmethod
    def spectrumViewer( **kwargs ):
        win = SpectrumWindow( **kwargs )
        Gtk.main()

    @staticmethod
    def imageViewer( **kwargs ):
        win = ImageViewer.ImageViewer( **kwargs )
        Gtk.main()

    @staticmethod
    def MonitorViewer( **kwargs ):
        Callable.monitorViewer( **kwargs )

    @staticmethod
    def monitorviewer( **kwargs ):
        Callable.monitorViewer( **kwargs )

    @staticmethod
    def monitorViewer( **kwargs ):
        win = MonitorViewer.MonitorViewer( **kwargs )
        try:
            Gtk.main()
        except KeyboardInterrupt:
            print 'KeyboardInterrupt'
            win.__del__()
            
    @staticmethod
    def generateImages( **kwargs ):
        lambda_ = float(kwargs['lambda'])
        sigma = float(kwargs['sigma'])
        v_os = int(kwargs['v_os'])
        os = int(kwargs['os'])
        ohdu = int(kwargs['ohdu'])
        rebin = int(kwargs['rebin'])
        gain = float(kwargs['gain'])
        shape = ( int(  kwargs['ny'])*rebin, int(kwargs['nx']) )
        fname = kwargs['o']
        image = SimulateImage.make_image( shape, v_os, os, ohdu, lambda_, sigma, gain, rebin )
        print 'highest energy in image', image.max()
        if not SimulateImage.writeFITS( image, fname ):
            print 'exiting without extraction'
            exit(0)
        if 'extract' in kwargs:
            SimulateImage.standardExtraction(fname, fname+'.root', verbose=False)
        exit(0)
        
    @staticmethod
    def dcObservations( **kwargs ):
        '''
        find the relation between the parameters of a convoluted poisson-normal distribution and the unbinned estimates by producing random distributions and then computing the principal component analysis to look for correlated variables
        '''
    
        # elections in eV
        eV_e = 3.745#eV/e

        #generate observations
        samples = 6
        min_gain = 500
        max_gain = 2500
        gain_adu_e_list = (max_gain - min_gain) * np.random.random_sample( samples ) + min_gain
        min_lambda = .01
        max_lambda = 2.
        lambda_list = (max_lambda - min_lambda) * np.random.random_sample( samples ) + min_lambda
        min_sigma = .1
        max_sigma = 3.
        sigma_e_list = (max_sigma - min_sigma) * np.random.random_sample( samples ) + min_sigma
        
        #define the desired observables
        observations = {
            'mean': lambda X: np.mean(X),
            'm2': lambda X: np.mean(X**2),
            'm3': lambda X: np.mean(X**3),
            'median': lambda X: np.median(X),
            'mad': lambda X: mad(X,scale=1),
            }
        ptiles = 4
        #loop in the parameters
        header = ['gain_adu_e', 'lambda', 'sigma_e']
        for key in observations.keys():
            header.append('%s_0'%key)
            header.append(key)
        for n in range(2**ptiles-1):
            header.append('q%d_0'%n)
            header.append('q%d'%n)
        print 'header', header, len(header)
        table = []
        for gain_adu_e in gain_adu_e_list:
            for lambda_ in lambda_list:
                for sigma_e in sigma_e_list:
                    #generate the rvs
                    N = int(1e6)
                    sigma_adu = sigma_e*gain_adu_e
                    frozen_pdf = lambda E, s=sigma_adu, mu=lambda_, g=gain_adu_e: stats.poisson_norm.pdf(E, scale=s, mu=mu, g=g )
                    norm_rvs = scipy.stats.norm.rvs(scale=sigma_adu, size=N)
                    Emin_adu, Emax_adu = 1.5*norm_rvs.min(), 3*norm_rvs.max()
                    poisson_norm_rvs = stats.poisson_norm.rvs( frozen_pdf, Emin_adu, Emax_adu, N )
                    
                    #sample multiple times the same rvs
                    for i in range(20):
                        sub_norm_rvs = np.random.choice(norm_rvs, int(1e4))
                        sub_poisson_norm_rvs = np.random.choice(poisson_norm_rvs, int(1e4))
                        
                        #add entries to table
                        entry = [ gain_adu_e, lambda_, sigma_e ]
                        for key, obs in observations.items():
                            entry.append( obs(sub_norm_rvs) )
                            entry.append( obs(sub_poisson_norm_rvs) )
                        sub_norm_tiles = tiles( sub_norm_rvs, ptiles )
                        sub_poisson_norm_tiles = tiles( sub_poisson_norm_rvs, ptiles )
                        #print 'len octatiles', len(sub_norm_tiles)
                        for n in range(2**ptiles-1):
                            entry.append( sub_norm_tiles[n] )
                            entry.append( sub_poisson_norm_tiles[n] )
                        #print 'len entry', len(entry)
                        table.append(entry)
                    print 'len', len(table)

        #append observations in file
        if os.path.exists('darkCurrentTraining.csv'):
            oldtable = np.genfromtxt('darkCurrentTraining.csv', names=True)
            print oldtable.dtype.names
            print len(oldtable), len(oldtable[0])
            oldtable = np.array( [ oldtable[name] for name in oldtable.dtype.names ] ).T
            extended = np.concatenate((oldtable,table), axis=0)
            print len(extended), len(extended[0])
            np.savetxt( 'darkCurrentTraining.csv', extended, header=' '.join(header) )
        else:
            np.savetxt( 'darkCurrentTraining.csv', table, header=' '.join(header) )
        return
    
    @staticmethod
    def dcAnalysis():
        
        table = []
        a = np.random.random_sample(1000)*10
        x = np.random.random_sample(1000)*10
        for a_,x_ in zip(a,x):
            eps = 2*(np.random.random_sample(1)[0]-.5)/10.
            y = a_*x_ + eps
            table.append( [a_, x_, y, y-x_*a_] )
            #for x_, y_ in zip(x,y):
                #table.append([a_, x_, y_,a_*x_])
        table = np.array(table).T
        print table
        print table.shape
        eqs = computeEqs( table, ['a','x','y','(y-a*x)'] )
        for coefs, variables, eqlabel, std in zip(*eqs):
            #print ''.join( [ '%+.4f*%s'%(coefs[i], variables[i]) for i in range(len(coefs)) ] ), std
            print ''.join( [ '%+.4f*%s'%(coefs[i], variables[i]) for i in range(len(coefs)) ] ), '= %+.4f*%s'%(std,eqlabel)
        return
        
        
        
        print 'read table'
        table = np.genfromtxt('darkCurrentTraining.csv', names=True )

        print len(table)
        print table.dtype.names
        print table['lambda']
        
        newVars = {
            'log(gain_adu_e)': lambda x: np.log(x['gain_adu_e']),
            'log(lambda)': lambda x: np.log(x['lambda']),
            #'var': lambda x: x['m2'] - x['mean']**2,
            #'var_0': lambda x: x['m2_0'] - x['mean_0']**2,
            'log(Dvar)': lambda x: np.log( np.abs( (x['m2'] - x['mean']**2) - (x['m2_0'] - x['mean_0']**2) )),
            'log(Dmean)': lambda x: np.log( np.abs( x['mean'] - x['mean_0'] )),
            }
        newtable = np.array( [ func(table) for key, func in newVars.items() ] )
        
        eqs = computeEqs( newtable, newVars.keys() )
        for coefs, variables, std in zip(*eqs):
            print ''.join( [ '%+.4f*%s'%(coefs[i], variables[i]) for i in range(len(coefs)) ] ), std
        return
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        getColumn = lambda s: newtable[ newVars.keys().index(s) ]
        #ax.scatter( 2.18e-4*getColumn('Dmean') + 2e-6*getColumn('gain_adu_e'), getColumn('Dvar') )
        #ax.scatter( getColumn('log(Dvar)') + .69*getColumn('log(Dmean)'), .34*getColumn('log(lambda)') )#.515*getColumn('log(gain_adu_e)') )
        ax.scatter( 2*getColumn('log(Dvar)') - 3*getColumn('log(Dmean)'), getColumn('log(gain_adu_e)') )#.515*getColumn('log(gain_adu_e)') )
        fig.savefig('dc_analysis.pdf')
        return
    
    
        for entry in table:
            #first estimate with mu=<E> overscan
            DeltaMean = entry['mean'] - entry['mean_0'] #ADU
            DeltaVar = entry['m2'] - entry['mean']**2 - (entry['m2_0'] - entry['mean_0']**2) #ADU**2
            l = DeltaMean**2/DeltaVar
            g = DeltaVar**2/DeltaMean**3 # ADU**4/ADU**3 = ADU
            print '1st estimate', entry['lambda'], entry['gain_adu_e'], l, g #ADU/ADU
        
    @staticmethod
    def paperTable():
        count = {}
        #reactor_subruns = {'on': ['029G','029H','029I','029J','031A','031B'], 'off': ['009','029C', '029D', '029E', '029F']}
        reactor_subruns = {'on': [
            #'018A', '018B', '018C', '018D', '018E', '018F',
            #'025A', '025B', '025C', '025D', '025E', '025F',
            #'026A', '026B',
            #'028A', '028B', '028C', '028D', '028E',
            '029A', '029B', '029C', '029D', '029E', '029F', '029G','029H','029I','029J',
            '031A', '031B', '031C', '031D', '031E', '031F', '031G','031H','031I',#'031J',
            #'032A', '032B', '032C', '032D',
            ]}
        selections = {#'3-7keV': 'flag==0 && E2/gainCu>3 && E2/gainCu<7', 
                      '250-350keV': 'flag==0 && E0/gainCu>250 && E0/gainCu<350'}
        for reactor, subruns in reactor_subruns.items():
            for range_, selection in selections.items():
                fullcount = []
                for subrun in subruns:
                    count = []
                    print subrun
                    path = ConniePaths.catalogPath( subrun, gain=True, skim=True )
                    print path
                    data = root_numpy.root2array( path, treename = 'hitSumm', branches = ['runID','E2','gainCu','ohdu', 'expoStart'], selection=selection )
                    runIDs = np.unique( data['runID'] )
                    expoStarts = np.unique( data['expoStart'] )
                    print data['runID']
                    ohdus = np.unique( data['ohdu'] )
                    print ohdus
                    for runID, expoStart in zip(runIDs, expoStarts):
                        for ohdu in ohdus:
                            mask = np.all( [ data['runID'] == runID, data['ohdu'] == ohdu ], axis=0 )
                            count.append( [runID, expoStart, ohdu, len(data[mask])] )
                            fullcount.append( [runID, expoStart, ohdu, len(data[mask])] )
                    np.savetxt( 'count_%s_%s_%s.csv'%(reactor, range_, subrun), count, header='runID expoStart ohdu count', fmt='%s', delimiter=' ' )
                    print 'done', reactor, range_, subrun
                np.savetxt( 'count_%s_%s.csv'%(reactor, range_), count, header='runID expoStart ohdu count', fmt='%s', delimiter=' ' )
                print 'done', reactor, range_

    @staticmethod
    def plotPaper():
        reactors = ['on',#'off'
                    ]
        ranges = [#'3-7keV', 
            '250-350keV']
        for range_ in ranges:
            fig = plt.figure()
            w, h = fig.get_size_inches()
            fig.set_size_inches((w,h*.75))
            ax = fig.add_subplot(111)
            X_cum = []
            XX_cum = []
            for reactor in reactors:
                data = np.genfromtxt('count_%s_%s.csv'%(reactor, range_), names=True )
                ax.plot( data['expoStart'], data['count'], 'b.' ) 

                #if range_ == '250-350keV':
                    ##X = X[ X>600 ]
                    #pass
                #hist, bins = np.histogram( X, bins=np.arange(X.min(), X.max(), 2 if range_ == '3-7keV' else 15) )
                ##hist, bins = np.histogram( X, bins=np.arange(X.min(), X.max(), 1 if range_ == '3-7keV' else 5) )
                #bins = .5*(bins[1:]+bins[:-1])
                #color = 'blue' if reactor == 'on' else 'red'
                #ax.step(bins, hist.astype(float)/np.sum(hist), where='mid', color=color)
                
                ##fit unbinned
                #mean, std = scipy.stats.norm.fit(X)
                #N = len(X)
                #mean_err = std/np.sqrt(N)
                #std_err = np.sqrt( np.mean( (X-mean)**4 )/N )/std

                ##fit binned
                #popt, pcov = scipy.optimize.curve_fit( lambda x, loc, scale, A: A*scipy.stats.norm.pdf(x, loc, scale), bins, hist.astype(float)/np.sum(hist)/(bins[1]-bins[0]), p0=(mean,std,1) )
                #perr = np.sqrt(np.diag(pcov))
                #print 'popt, perr', zip(popt, perr)
                #func = lambda x, *p: (bins[1]-bins[0])*p[2]*scipy.stats.norm.pdf(x, loc=p[0], scale=p[1])
                #chisq = scipy.stats.chisquare( hist, f_exp = np.sum(hist)*func(bins, *popt), ddof = len(popt) )[0]/(len(hist)-len(popt))
                #print 'chisq.fit', chisq
                #chisq = scipy.stats.chisquare( hist, f_exp = np.sum(hist)*func(bins, mean, std, 1), ddof = 2 )[0]/(len(hist)-2)
                #print 'chisq,MLE', chisq
                
                ##label = 'Reactor %s\n$\mu=%s$\n$\sigma=%s$'%(reactor.upper(), print_value_error(mean,mean_err), print_value_error(std, std_err))
                #label = '$\mu=%s$\n$\sigma=%s$'%(print_value_error(mean,mean_err), print_value_error(std, std_err))
                #l = np.arange( bins.min(), bins.max(), .1 )
                #print len(l), len(func(l,*popt))
                #ax.plot( l, (bins[1]-bins[0])*scipy.stats.norm.pdf(l, loc=mean, scale=std), color=color, label= label)
                ##ax.plot( l, func(l, *popt), color=color, ls=':')
                #np.savetxt( 'hist_%s_%s.csv'%(reactor,range_), zip(bins,hist), header='bin hist', fmt='%s', delimiter=' ' )
                #print 'save', 'hist_%s_%s.csv'%(reactor,range_)
            #print X_cum
            #print XX_cum
            ax.minorticks_on()
            ax.set_xlabel('Rate (events/image)', x=1, ha='right')
            ax.set_ylabel('Relative frequency', y=1, ha='right')
            fig.set_tight_layout(True)
            ax.legend( fancybox=True, framealpha=0 )
            fig.savefig( 'timeSeries_%s.pdf'%(range_) )
            print 'save', 'timeSeries_%s.pdf'%(range_)
        print 'done'
            
    
    @staticmethod
    def simulateImageAndFit( **kwargs ):
        simulatedImage = Image.generate( **kwargs )
        simulatedImage.plotSpectrum( **kwargs )
    
    @staticmethod
    def histogramFits( **kwargs ):
        from sklearn import mixture
        import itertools

        def printAverageEstimates( a ):
            E = a.flatten()
            print 'ave', np.mean( E )
            print 'ave', np.sum( E )/len(E)
            print 'aveInv', np.sum( 1./(E**2+1)*E )/np.sum( 1./(E**2+1) )
            print 'aveInv1', np.average( E, weights=1./(E**2+1) )
            print 'aveInv10', np.average( E, weights=1./(E**2+10) )
            print 'aveInv100', np.average( E, weights=1./(E**2+100) )
            print 'aveInv1e3', np.average( E, weights=1./(E**2+1e3) )
            print 'aveInv1e4', np.average( E, weights=1./(E**2+1e4) )

        def printStdEstimates( a ):
            E = a.flatten()
            print 'ave', np.var( E )
            print 'ave', np.sum( E**2 )/len(E) - (np.sum( E )/len(E))**2
            print 'aveInv1', np.average( E**2, weights=1./(E**2+1) ) - np.average( E, weights=1./(E**2+1) )**2
            print 'aveInv10', np.average( E**2, weights=1./(E**2+10) ) - np.average( E, weights=1./(E**2+10) )**2
            print 'aveInv100', np.average( E**2, weights=1./(E**2+100) ) - np.average( E, weights=1./(E**2+100) )**2

        def compute2DParams( E, N=3 ):
            print 'computing 2D params with radius', N
            E0 = np.array(E)
            E1 = np.array(E)
            x = np.arange(0,E1.shape[0],1)[:,None]
            y = np.arange(0,E1.shape[1],1)[None,:]
            E1 -= np.min(E1)
            xE = np.zeros( np.array(E1.shape) - [N,N] )
            yE = np.zeros( np.array(E1.shape) - [N,N] )
            x2E = np.zeros( np.array(E1.shape) - [N,N] )
            y2E = np.zeros( np.array(E1.shape) - [N,N] )
            mE = np.zeros( np.array(E1.shape) - [N,N] )
            Emax = np.zeros( np.array(E1.shape) - [N,N] )
            print xE.shape
            for i in range(0,N):
                for j in range(0,N):
                    xE += x[i:-N+i,:] * E1[i:-N+i,j:-N+j]
                    yE += y[:,j:-N+j] * E1[i:-N+i,j:-N+j]
                    x2E += x[i:-N+i,:]**2 * E1[i:-N+i,j:-N+j]
                    y2E += y[:,j:-N+j]**2 * E1[i:-N+i,j:-N+j]
                    mE += E1[i:-N+i,j:-N+j]
                    Emax = np.amax( [Emax, E0[i:-N+i,j:-N+j]], axis=0 )
            xE = xE/mE
            yE = yE/mE
            print xE.shape, yE.shape, Emax.shape
            sE2 = (x2E + y2E)/mE - (xE**2 + yE**2)
            #print sE2, sE2.shape
            return sE2 
            #return Emax

        def plot( hists, lines = [], xlabel = None, ylabel = None ):
            fig = plt.figure(figsize=(3,3))
            ax = fig.add_subplot(111)
            for ia in hists:
                ax.hist( ia.flatten(), bins=200, histtype='step', normed=True )
            for ia in lines:
                ax.plot( ia[0], ia[1]/factor )
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid(True)
            #ax.get_yaxis().get_major_formatter().set_useOffset(False)
            fig.subplots_adjust()
            fig.tight_layout()
            fig.savefig('sE2.png')
        
        def plotMultiple( m, merr, s, serr, l, lerr, g, gerr):
            fig = plt.figure(figsize=(3,3))
            ax = fig.add_subplot(111)
            ax.errorbar( )
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid(True)
            #ax.get_yaxis().get_major_formatter().set_useOffset(False)
            fig.subplots_adjust()
            fig.tight_layout()
            fig.savefig('stats.png')
        
        def plot_results_1d(X, Y_, means, covariances, weights, index, title):
            color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold', 'darkorange'])

            splot = plt.subplot(2, 1, 1 + index)
            splot = plt.subplot(1, 1, 1)
            _, bins, _ = plt.hist( X, bins=500, histtype='step', normed=True )
            for i, (mean, covar, weight, color) in enumerate(zip(means, covariances, weights, color_iter)):
                print mean, covar
                if not np.any(Y_ == i):
                    continue
                x = .5*(bins[1:]+bins[:-1])
                plt.plot( x, weight*stats.norm.pdf(x,mean,covar[0]) )
            plt.title(title)
            #plt.yscale('log')
            #plt.ylim((1,200))
            plt.savefig('test.png')

        #None
        #print 'matrix'
        #print stats.ufloat.show_conversion_matrix()
        if 'stats' in kwargs:
            v_os = 70
            os_ = 150
            shape = (4150, 4112)
            ohdu = 3
            if 'divide' in kwargs: divide = int( kwargs['divide'] )
            else: divide = 1
            
            if 'rebin' in kwargs: rebin = int( kwargs['rebin'] )
            else: rebin = 5
            
            if 'dE' in kwargs: dE = float( kwargs['dE'] )
            else: dE = 1.
            
            N = int(kwargs['stats'])
            
            print 'N=', N
            if 'robust' in kwargs:
                if 'fit' in kwargs:
                    fbasename = 'stats_fit_robust_%s_dE%s_div%s'%(rebin,dE, divide)
                else:
                    fbasename = 'stats_robust_%s'%rebin
            else:
                if 'fit' in kwargs:
                    fbasename = 'stats_fit_%s_dE%s'%(rebin,dE)
                else:
                    fbasename = 'stats_%s'%rebin
            timer = Timer.LoopTimer(n=N, name='fits')
            for i in range(N):
                unit = 'ADU%s'%ohdu
                timer.eta(i)
                timer.end(i)
                #loc = stats.ufloat( (np.random.random_sample(1)[0]-.5), unit='e-' )
                loc = stats.ufloat( 0, unit='e-' )
                #sigma = stats.ufloat( 3*np.random.random_sample(1)[0], unit='e-' )
                sigma = stats.ufloat( 1., unit='e-' )
                #lambda_ = 1.5*np.random.random_sample(1)[0]
                lambda_ = .5
                #gain = 1000*np.random.random_sample(1)[0] + 1000
                gain = 1500
                mbs = False
                if 'mbs' in kwargs:
                    mbs = True
                g_e = gain*eIonization*1e-3
                stats.ufloat.add_conversion( 'eV', unit, gain*1e-3 )
                
                print 'matrix'
                print stats.ufloat.show_conversion_matrix()
                #fitsimage = SimulateImage.make_image( shape, v_os, os_, ohdu, lambda_, sigma, gain, rebin, mbs )
                fitsimage = ConnieImage.readImage_raw( 7000, 4 )
                #fitsimage += loc
                #fitsimage = fitsimage.asunit(unit)
                #realimage = stats.ufloat( fitsimage.val, unit = 'ADU%s_'%ohdu )
                realimage = stats.ufloat( fitsimage, unit = 'ADU%s_'%ohdu )
                active, overscan = separate_active_overscan( realimage )
                print active.shape, active

                # Fit a Dirichlet process Gaussian mixture using five components
                X_ = active.flatten()
                X_ = X_[X_< 1e9]
                X_ = X_[:int(1e6)]
                print X_.shape
                X = zip(X_)
                dpgmm = mixture.GaussianMixture(n_components=3, covariance_type='full', tol=1e-5, max_iter=200).fit(X)
                print 'done', dpgmm.means_
                print dpgmm.covariances_
                plot_results_1d(X_, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, dpgmm.weights_, 1, 'Bayesian Gaussian Mixture with a Dirichlet process prior')
                
                exit(0)
                
                if 'robust' in kwargs:
                    if 'fit' in kwargs:
                        res = image_MLE_fit_robust( active, overscan, dE, divide )
                    else:
                        res = image_MLE_robust( active, overscan )
                else:
                    if 'fit' in kwargs:
                        res = image_MLE_fit( active, overscan, dE, divide )
                    else:
                        res = image_MLE( active, overscan )
                
                print 'mu', res[0], loc.asunit('ADU3')
                print 'sigma', res[1], sigma.asunit('ADU3')
                print 'lambda_', res[3], lambda_
                print 'gain', res[2], g_e
                
                fig = plt.figure(figsize=(3,3))
                ax = fig.add_subplot(111)
                _, bins, _ = ax.hist( active.flatten(), bins=200, histtype='step', normed=True )
                x = .5*(bins[1:]+bins[:-1])
                y = stats.poisson_norm.pdf(x, res[0].val, res[1].val, res[2].val, res[3].val)
                ax.plot( x, y )
                ax.set_xlabel('E[ADU]')
                ax.set_ylabel('freq')
                ax.grid(True)
                fig.subplots_adjust()
                fig.tight_layout()
                fig.savefig('sE2.png')
                
                exit(0)
                
                entry = '%s  %s %s %s' % (loc, sigma, g_e, lambda_ )
                mle = '%s  %s %s %s' % (res.mu, res.sigma, res.g, res.lamb)
                mle_err = '%s  %s %s %s' % (err[0]/g_e, err[1]/g_e, err[2], err[3]/rebin)
                exit(0)
                if not os.path.exists(fbasename+'.csv'):
                    with open('stats.csv','a+') as f:
                        f.write('#loc scale g lambda mu sigma gain lambda muErr sigmaErr gainErr lambdaErr\n')
                with open(fbasename+'.csv','a+') as f:
                    line = '%s %s %s\n' % ( entry, mle, mle_err )
                    print line
                    f.write(line)
            fig = plt.figure()
            ax_mu = fig.add_subplot(221)
            ax_mu.set_xlabel('mu')
            ax_mu.set_ylabel('err')
            ax_mu.grid(True)
            ax_sigma = fig.add_subplot(222)
            ax_sigma.set_xlabel('sigma')
            ax_sigma.set_ylabel('err')
            ax_sigma.grid(True)
            ax_gain = fig.add_subplot(223)
            ax_gain.set_xlabel('gain')
            ax_gain.set_ylabel('err')
            ax_gain.grid(True)
            ax_lambda = fig.add_subplot(224)
            ax_lambda.set_xlabel('lambda')
            ax_lambda.set_ylabel('err')
            ax_lambda.grid(True)
            data = np.loadtxt(fbasename+'.csv')
            data = data.T

            ax_mu.errorbar(data[0], data[4]-data[0], data[8], fmt='_', elinewidth=.5)
            ax_sigma.errorbar(data[1], data[5]-data[1], data[9], fmt='_', elinewidth=.5)
            ax_gain.errorbar(data[2], data[6]-data[2], data[10], fmt='_', elinewidth=.5)
            ax_lambda.errorbar(data[3], data[7]-data[3], data[11], fmt='_', elinewidth=.5)
            ax_mu.set_ylim((-1,1))
            ax_sigma.set_ylim((-1,1))
            ax_gain.set_ylim((-1,1))
            ax_lambda.set_ylim((-1,1))
            fig.tight_layout()
            fig.savefig(fbasename+'.pdf')
            exit(0)
        else:
            if 'p' in kwargs:
                plotParamsFromCSVList(glob.glob(kwargs['p']), kwargs['o'])
                exit(0)
            elif 'f' in kwargs and 'o' in kwargs:
                filenamelist = sorted( glob.glob(kwargs['f']) )
                print 'input files', filenamelist
                fulltable = []
                timer = Timer.LoopTimer( len(filenamelist), 'loop' )
                for i, filename in enumerate(filenamelist):
                    print '%s(%s)' % ( i+1, len(filenamelist) )
                    runID = int( re.search( r'runID_[0-9]+?_([0-9]+?)_Int', filename ).groups()[0] )
                    header, table = computeParamsFromFitsFile( filename )
                    fullheader = ['runID'] + header
                    fulltable += [ [runID] + row for row in table ]
                    with open( kwargs['o'], 'wb') as output:
                        np.savetxt( output, fulltable, header = ', '.join( fullheader ), delimiter = ', ', fmt='%s' )
                    timer.eta(i)
                    timer.end(i)
                exit(0)
            else:
                v_os = 70
                #os = 550
                os_ = 150
                shape = (4150, 4112)
                ohdu = 3
                rebin = int(kwargs['rebin'])
                sigma = float(kwargs['sigma'])
                lambda_ = float(kwargs['lambda'])
                gain = float(kwargs['gain'])
                mbs = False
                if 'mbs' in kwargs:
                    mbs = True
                fitsimage = SimulateImage.make_image( shape, v_os, os_, ohdu, lambda_, sigma, gain, rebin, mbs )

            #print 'res', res
        exit(0)
        return
    
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
    print( colored('Tools for the CONNIE colaboration', attrs=['bold']) )
    print( colored('by Philipe Mota (philipe.mota@gmail.com)', attrs=['bold']) )
    print( colored('repository https://github.com/PhMota/CONNIEtools', attrs=['bold']) )
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
    option = ' '.join(sys.argv).split('--')
    if len(option) > 1:
        option = option[-1]
        method = option.split(' ')[0]
        subargs = option.split(' -')[1:]
        kwargs = {}
        for arg in subargs:
            pair = arg.split(' ')
            if len(pair) > 1:
                kwargs[pair[0]] = pair[1]
            else:
                kwargs[pair[0]] = None
        print 'called with', method, kwargs
        getattr(Callable, method)( **kwargs )
        exit(0)
    print( 'callables', Callable.__dict__ )
    print( 'options not recognized', sys.argv )
    help_()
    exit(0)
