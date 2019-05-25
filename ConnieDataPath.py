# coding: utf-8

import glob
import os
#import shutil
#import sys
import re
#import copy
#import root_numpy
#import datetime
#import time
#import random
#import timeit
#import traceback

#import json

#import gi
#gi.require_version('Gtk', '3.0')
#from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib
#import threading

def isint(a):
    try:
        int(a)
        return True
    except ValueError:
        return False

class ConnieDataPath:
    connie_folder = '/share/storage2/connie/'
    run_pattern = r'/runs/([0-9]+?)/'
    subrun_pattern = r'/runs/([0-9]+?[A-Z])/'
    runID_pattern = r'runID_([0-9]+?)_([0-9]+?)_.*_p([0-9]).fits.fz'
    range_pattern = r'/data_([0-9]+?_to_[0-9]+)'
    
    @classmethod
    def parse_run(cls, run):
        if run is None: return '*'
        if type(run) is int: return '%03d'%run
        if type(run) is str: return run
        return run
    
    @classmethod
    def run_folder(cls, run):
        return cls.connie_folder+'data/runs/%s/'%( cls.parse_run(run) )
    
    @classmethod
    def run_folder_processed(cls, run):
        return cls.connie_folder+'data_analysis/processed02_data/runs/%s/'%( cls.parse_run(run) )

    @classmethod
    def parse_runID(cls, runID):
        if runID is None: return '*'
        if type(runID) is int: return '%05d'%runID
        return runID
    
    @classmethod
    def raw_pattern( cls, run, runID, part ): 
        return cls.run_folder(run)+'runID_%s_%s_*_p%s.fits.fz'%( cls.parse_run(run), cls.parse_runID(runID), part )

    @classmethod
    def processed_pattern( cls, run, runID, part, image ): 
        return cls.run_folder_processed(run)+'data_[0-9]*_to_[0-9]*/%s/images/%s_*_runID_%s_%s_*_p%s.fits'%( image, image, cls.parse_run(run), cls.parse_runID(runID), part )

    @classmethod
    def masterBias_pattern( cls, run, part ): 
        return cls.run_folder_processed(run)+'data_[0-9]*_to_[0-9]*/mbs/masterBias_*_p%s.fits'%( part )
    
    @classmethod
    def catalog_pattern( cls, subrun, skim=False ):
        return cls.run_folder_processed(subrun)+'data_[0-9]*_to_[0-9]*/ext/catalog/catalog_data_[0-9]*_to_*[0-9][0-9]%s.root'%( '.skim1' if skim else '')
    
    @classmethod
    def runPath(cls, run=None, runID=None):
        if run is None: return sorted( glob.glob( cls.run_folder('*') ) )
        if type(run) is int: return glob.glob( cls.run_folder(run) )
        if type(run) is list: return map( cls.runPaths, run )
        raise Exception( 'type not supported: %s'%( ' '.join([ '%s:%s'%(var, type(val)) for var, val in locals().items() ]) ) )
    
    @classmethod
    def run(cls, path=None, runID=None):
        if path is None and runID is None: return cls.run(cls.runPath())
        if type(path) is str: return int( re.search( cls.run_pattern, path ).groups()[0] )
        if type(path) is list: return map( cls.run, path )
        if type(runID) is str: return int( re.search( cls.runID_pattern, runID ).groups()[0] )
        if type(runID) is int: return cls.run( path=cls.runIDPath(runID=runID)[0] )
        if type(runID) is list: return map( cls.run, cls.runIDPath(runID) )
        if type(runID) is tuple: return sorted(list(set(cls.run(runID=list(runID)))))
        raise Exception('type not supported path=%s runID=%s'%(type(path), type(runID)) )

    @classmethod
    def runIDPath(cls, runID=None, run=None, part='*'):
        if runID is None and run is None:  return sorted( glob.glob( cls.raw_pattern('*','*',part) ) )
        if type(runID) is int: return sorted( glob.glob( cls.raw_pattern('*',runID,part) ) )
        if type(runID) is list: return map( cls.runIDPath, runID )
        if type(run) is int: return sorted( glob.glob( cls.raw_pattern(run,'*',part) ) )
        if type(run) is list: return map( lambda run: cls.runIDPath(run=run), run )
        raise Exception( 'type not supported: %s'%( ' '.join([ '%s:%s'%(var, type(val)) for var, val in locals().items() ]) ) )
    
    @classmethod
    def runID(cls, path=None, run=None, part='*'):
        if path is None and run is None: return cls.runID(cls.runIDPath(part='1'))
        if type(path) is str: return int( re.search( cls.runID_pattern, path ).groups()[1] )
        if type(path) is list: return map( cls.runID, path )
        if type(run) is int: return cls.runID(path=cls.runIDPath(run=run, part='1'))
        if type(run) is list: return map( lambda run: cls.runID(run=run), run )
        if type(run) is tuple: return [ item for line in cls.runID(run=list(run)) for item in line ]
        raise Exception( 'type not supported: %s'%( ' '.join([ '%s:%s'%(var, type(val)) for var, val in locals().items() ]) ) )
    
    @classmethod
    def runIDPathProcessed( cls, runID=None, run=None, part='*', image='*' ):
        if runID is None and run is None:  return sorted( glob.glob( cls.processed_pattern( '*', '*', part, image) ) )
        if type(runID) is int: return sorted( glob.glob( cls.processed_pattern( '*', runID, part, image ) ) )
        if type(runID) is list: return map( lambda r: cls.runIDPathProcessed( r, run, part, image ), runID )
        if type(run) is str: return sorted( glob.glob( cls.processed_pattern( run, '*', part, image ) ) )
        if type(run) is list: return map( lambda r: cls.runIDPathProcessed( run=r, part=part, image=image ), run )
        raise Exception( 'type not supported: %s'%( ' '.join([ '%s:%s'%(var, type(val)) for var, val in locals().items() ]) ) )

    @classmethod
    def runIDProcessed(cls, path=None, subrun=None, part='*', image='*'):
        if path is None and run is None and subrun is None: return cls.runID(cls.runIDPathProcessed( part='1', image=image))
        if type(path) is str: return int( re.search( cls.runID_pattern, path ).groups()[1] )
        if type(path) is list: return map( lambda p: cls.runIDProcessed( p, subrun, part, image ), path )
        if type(subrun) is int: return cls.runIDProcessed( path=cls.runIDPathProcessed( run=run, subrun=None, part='1', image=image ))
        if type(subrun) is list: return map( lambda r: cls.runIDProcessed( run=r, subrun=None, part=part, image=image ), run )
        if type(subrun) is tuple: return [ item for line in cls.runIDProcessed( run=list(run), subrun=None, part=part, image=image ) for item in line ]
        raise Exception( 'type not supported: %s'%( ' '.join([ '%s:%s'%(var, type(val)) for var, val in locals().items() ]) ) )

    @classmethod
    def subrunPath(cls, subrun=None, runID=None):
        if subrun is None: return sorted( glob.glob( cls.run_folder_processed('*') ) )
        if type(subrun) is str: return glob.glob( cls.run_folder_processed(run) )
        if type(subrun) is list: return map( cls.subrunPaths, subrun )
        raise Exception( 'type not supported: %s'%( ' '.join([ '%s:%s'%(var, type(val)) for var, val in locals().items() ]) ) )

    @classmethod
    def subrun(cls, path=None, runID=None):
        if path is None and runID is None: return cls.subrun(cls.subrunPath())
        if type(path) is str: return re.search( cls.subrun_pattern, path ).groups()[0]
        if type(path) is list: return map( cls.subrun, path )
        if type(runID) is str: return int( re.search( cls.subrunID_pattern, runID ).groups()[0] )
        if type(runID) is int: return cls.subrun( path=cls.runIDPathProcessed(runID=runID)[0] )
        if type(runID) is list: return map( cls.run, cls.runIDPathProcessed(runID) )
        if type(runID) is tuple: return sorted(list(set(cls.subrun(runID=list(runID)))))
        raise Exception('type not supported path=%s runID=%s'%(type(path), type(runID)) )
    
    @classmethod
    def masterBiasPath( cls, runID=None, subrun=None, part='*' ):
        if type(subrun) is str: return glob.glob( cls.masterBias_pattern(subrun,part) )
        if type(runID) is int: return cls.masterBiasPath( subrun=cls.subrun(runID=runID), part=part )
        raise Exception('type not supported path=%s runID=%s'%(type(path), type(runID)) )
        
    @classmethod
    def range_( cls, path=None, gain=False ):
        if type(path) is list: #map( lambda path, gain=gain: cls.range_(path,gain), path )
            return cls.range_( path[0], gain )
        if type(path) is str: 
            return re.search( cls.range_pattern, path ).groups()[0]
    
    @classmethod
    def catalogPath( cls, subrun=None, skim=False, gain=False ):
        if gain:
            if type(subrun) is str: 
                return glob.glob( '/share/storage2/connie/nu_processing/scripts/ProcCat/*cut_scn_osi_raw_gain_catalog_data_%s.root'%(cls.range_(cls.catalogPath(subrun,skim,False))) )[0]
        else:
            if type(subrun) is str: 
                return glob.glob( cls.catalog_pattern( subrun, skim ) )
        raise Exception('type not supported path=%s runID=%s'%(type(path), type(runID)) )
    
    @staticmethod
    def test():
        print DataPath.parse_run(None)
        print DataPath.run()
        paths = DataPath.runPath()
        print paths
        #print DataPath.run( DataPath.runPath( DataPath.run(paths)[-1] ) )
        print DataPath.runID(run=[42,43])
        print DataPath.runID(run=(42,43))
 
if __name__ == '__main__':
    #print '\n'.join( ConnieDataPath.runIDPathProcessed(5000) )
    #print ConnieDataPath.subrun( runID=5000 )
    print ConnieDataPath.range_(ConnieDataPath.catalogPath('035A'))
    print '028A', ConnieDataPath.catalogPath('028A')
    print 'range_ 029A', ConnieDataPath.range_(ConnieDataPath.catalogPath('029A'))
    for on in ['029G','029H','029I','029J','031A','031B']:
        print on, ConnieDataPath.catalogPath(on, gain=True, skim=True)
    print
    for off in ['009','029C', '029D', '029E', '029F']:
        print off, ConnieDataPath.catalogPath(off, gain=True, skim=True)
    
