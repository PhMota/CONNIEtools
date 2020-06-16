#!/usr/bin/python

from __future__ import print_function
import os, sys, argparse

from Timer import Timer
import warnings
warnings.filterwarnings("ignore")

try:
    import root_numpy
except ImportError:
    print('missing module, please run')
    print('module load softwares/root/5.34-gnu-5.3')
    exit(0)

import glob
import matplotlib
matplotlib.use('gtk3agg')

import matplotlib.pylab as plt

from numpy.lib.recfunctions import append_fields, stack_arrays
import numpy as np

from ROOT import TFile, TTree, AddressOf

#import rootpy
#import rootpy.io
#from array import array
#from rootpy.tree import Tree, TreeModel, IntCol, IntArrayCol, FloatCol, FloatArrayCol, CharArrayCol
#from rootpy.io import root_open

import Image
from termcolor import colored
from PrintVar import print_var


class HitSumm(np.recarray):
    
    def __new__( cls, list_of_entries, levels, threshold, border ):
        print( 'new hits' )
        dtype = [
            ('ePix', object),
            ('xPix', object),
            ('yPix', object),
            ('level', object),
            ]
        obj = np.array( list_of_entries[1:], dtype = dtype ).view(np.recarray).view(cls)
        obj = rfn.append_fields( obj, 
                                 ('nSavedPix'), 
                                 np.array(map(len, obj.ePix)), 
                                 dtypes=(int), 
                                 asrecarray=True 
                                 ).view(cls)

        obj.border = border
        obj.threshold = threshold
        obj.background = list_of_entries[0]
        return obj
    
    def add_fields( self, names, arrays, types ):
        print( 'added', names )
        return rfn.append_fields( self, names, arrays, dtypes=types, asrecarray=True ).view(Hits)
        
    def compute_number_of_pixel_levels( self, lvl ):
        return np.array(map( lambda entry: len(entry.ePix[ entry.level <= lvl ]), self ))
        
    def add_number_of_pixel_levels( self, lvls ):
        added_fields = []
        for lvl in lvls:
            new_field = 'n%d' % lvl
            self = self.add_fields( 'n%d' % lvl, self.compute_number_of_pixel_levels(lvl), int )
            added_fields.append( new_field )
        return self
        
    def compute_energy_levels( self, lvl ):
        return np.array(map( lambda entry: np.sum(entry.ePix[ entry.level <= lvl ]) , self ))
    
    def add_energy_levels( self, lvls ):
        added_fields = []
        for lvl in lvls:
            new_field = 'E%d' % lvl
            self = self.add_fields( new_field, self.compute_energy_levels(lvl), float )
            added_fields.append( new_field )
        return self
        
    def __average__( self, xy, weights, level, lvl ):
        below_lvl = level <= lvl
        try:
            ret = np.average( xy[:,below_lvl], weights = weights[below_lvl], axis=-1 )
        except:
            ret = np.mean( xy[:,below_lvl], axis=-1 )
        return np.abs(ret)

    def compute_barycenter_levels(self, lvl):
        barycenter = lambda entry: self.__average__( np.array((entry.xPix, entry.yPix)), entry.ePix, entry.level, lvl )
        return np.array(map( barycenter, self ))
        
    def add_barycenter_levels( self, lvls ):
        added_fields = []
        for lvl in lvls:
            new_fields = ['xBary%d' % lvl,'yBary%d' % lvl]
            self = self.add_fields( new_fields, self.compute_barycenter_levels(lvl).T, (float, float) )
            added_fields.extend( new_fields )
        return self
    
    def compute_variance_levels( self, lvl ):
        variance = lambda entry: np.abs(
            self.__average__( np.array((entry.xPix, entry.yPix))**2, entry.ePix, entry.level, lvl ) - self.__average__( np.array((entry.xPix, entry.yPix)), entry.ePix, entry.level, lvl )**2 )
        
        return np.array(map( variance, self ))
        
    def add_variance_levels( self, lvls ):
        added_fields = []
        for lvl in lvls:
            new_fields = ['xVar%d' % lvl,'yVar%d' % lvl]
            self = self.add_fields( new_fields, self.compute_variance_levels( lvl ).T, (float,float) )
            added_fields.extend( new_fields )
        return self
    
    def match_simulated_events( self, events, lvl, length ):
        import itertools
        linked_list = {}
        radius = 1
        shifts = range( -radius, radius+1 )
        listofshifts = [ shifts for i in range(2) ]
        nshifts = tuple(itertools.product(*listofshifts) )

        array_of_projections_events = events.get_xy_rebin()
        array_of_projections_indices_events = get_indices( array_of_projections_events, length=length )
        events_indices_set = set( map(tuple, array_of_projections_indices_events ) )
        
        for i, ni in enumerate(array_of_projections_indices_events):
            for ns in nshifts:
                try:
                    linked_list[ tuple(ni + ns) ].append(i)
                except:
                    linked_list[ tuple(ni + ns) ] = [i]

        array_of_projections_hits = np.array( (self['xBary%d' % lvl], self['yBary%d' % lvl]) ).T
        array_of_projections_indices_hits = get_indices( array_of_projections_hits, length=length )
        
        print( 'max evts', array_of_projections_events[0].min(), array_of_projections_events[0].max(), array_of_projections_events[1].min(), array_of_projections_events[1].max() )
        print( 'max hits', array_of_projections_hits[0].min(), array_of_projections_hits[0].max(), array_of_projections_hits[1].min(), array_of_projections_hits[1].max() )
        
        matched_indices = []
        matched_id = []
        matched_Q = []
        matched_z = []
        matched_distances = []
        dist = lambda a, b: np.sqrt( np.sum( (a-b)**2, axis=-1 ) )
        for hit_ind, hit_ni in enumerate(array_of_projections_indices_hits):
            try:
                event_ind = np.array(linked_list[tuple( hit_ni )])
            except:
                matched_id.append( 0 )
                matched_distances.append( 0 )
                matched_Q.append( 0 )
                matched_z.append( -1 )
                continue
            pos_hit = array_of_projections_hits[ hit_ind ]
            std_max_hit = np.sqrt( np.max( ( self[ hit_ind ]['xVar%d' % lvl], self[ hit_ind ]['yVar%d' % lvl] ) ) )
            pos_event = array_of_projections_events[event_ind]
            ds = dist( pos_hit, pos_event )
            is_matched = ds <= 2*std_max_hit
            event_ind_matched = event_ind[is_matched]
            if len( event_ind_matched ) == 0:
                matched_id.append( 0 )
                matched_distances.append( 0 )
                matched_Q.append( 0 )
                matched_z.append( -1 )
                continue
            
            if len( event_ind_matched ) > 1:
                matched_id.append( 2 )
                matched_distances.append( ds.min() )
                matched_Q.append( np.sum( events.q[event_ind_matched] ) )
                matched_z.append( -1 )
                
            if len( event_ind_matched ) == 1:
                matched_id.append( events.get_id_code()[event_ind_matched[0]] )
                matched_distances.append( ds[0] )
                matched_Q.append( events.q[event_ind_matched[0]] )
                matched_z.append( events.z[event_ind_matched[0]] )
                
            #for ei, d in zip(event_ind_matched, ds[is_matched]):
                #matched_indices.append( (hit_ind, ei, d, is_mult) )

        #print( 'matched id', matched_id )
        #matched = np.array(matched)
        new_fields = ['id', 'Q', 'z', 'dist' ]
        self = self.add_fields( new_fields, (matched_id, matched_Q, matched_z, matched_distances), (int,int,float,float) )
        return self
        
        
    def compute_sizelike_each( self, entry, lvl, gain, sigma_noise, fast = True, tol = 1e-1 ):
        _Q, _mu, _sigma = size_like( entry.ePix/gain, entry.xPix, entry.yPix, 
                                E0 = entry['E%d'%lvl]/gain, 
                                mu0 = [ entry['xBary%d' % lvl], entry['yBary%d' % lvl] ], 
                                sigma0 = np.sqrt( [ entry['xVar%d'%lvl], entry['yVar%d'%lvl] ] ), 
                                sigma_noise0 = sigma_noise,
                                single_sigma = False,
                                fast = fast,
                                tol = tol)
        return _Q*gain, _mu[0], _mu[1], _sigma[0], _sigma[1]

    def add_sizelike( self, lvl, gain, sigma_noise, fast = True, tol = 1e-1 ):
        new_fields = np.array( map( lambda entry: self.compute_sizelike_each( entry, lvl, gain, sigma_noise, fast, tol ), self ) )
        self = self.add_fields( ['EL', 'xMu', 'yMu', 'xSigma', 'ySigma'], new_fields.T, (float, float, float, float, float) )
        return self
        
    #def get_sizeLike( self, lvl, gain, sigma_noise, fast=True, tol = None ):
        #Q = []
        #mu, sigma = [], []
        #for entry in self:
            #_Q, _mu, _sigma = size_like( entry.ePix/gain, entry.xPix, entry.yPix, 
                                  #E0 = entry['E%d'%lvl]/gain, 
                                  #mu0 = [ entry['xBary%d' % lvl], entry['yBary%d' % lvl] ], 
                                  #sigma0 = np.sqrt( [ entry['xVar%d'%lvl], entry['yVar%d'%lvl] ] ), 
                                  #sigma_noise0 = sigma_noise,
                                  #single_sigma=False,
                                  #fast = fast,
                                  #tol = tol)
            #Q.append( _Q )
            #mu.append( _mu )
            #sigma.append( _sigma )
        #return np.array(Q)*gain, np.array(mu), np.array(sigma)
    
    def get_Bary(self):
        return np.array( (self.xBary, self.yBary) ).T

    def get_Var(self):
        return np.array( (self.xVar, self.yVar) ).T

def make_hitSumm( list_of_clusters, levels ):
    dtype = [
        ('ePix', object),
        ('xPix', object),
        ('yPix', object),
        ('level', object),
        ]
    obj = np.array( list_of_clusters, dtype = dtype ).view(np.recarray).view(cls)
    obj = rfn.append_fields( obj, 
                                ('nSavedPix'), 
                                np.array(map(len, obj.ePix)), 
                                dtypes=(int), 
                                asrecarray=True 
                                ).view(HitSumm)
    
    return obj
    
#def list_branches( file, treename = 'hitSumm' ):
    #return root_numpy.list_branches( file, treename = treename )

#def read_catalog( file, treename = 'hitSumm', branches = None, selection = '', indices = None ):
    #if branches is None: branches = list_branches( file, treename )
    #if indices is None:
        #start, stop = None, None
    #else:
        #start, stop = np.min(indices), np.max(indices)+1
    #return root_numpy.root2array( file, treename = treename, branches = branches, selection = selection, start = start, stop = stop )

#def list_runIDs( file ):
    #return np.unique( read_catalog( file, branches = ['runID'] )['runID'] )

#def get_number_of_entries( file ):
    #return len( read_catalog( file, branches = ['runID'] )['runID'] )

#def save_catalog( file, data ):
    #root_numpy.array2root( data, file, treename='hitSumm', mode='recreate' )

#def parse_selection( selection ):
    #pass
    
#def test( *args ):
    #print( args )
    #file = args[1]
    #selection = args[2]
    #print( list_branches(file) )
    #runID_list = list_runIDs( file )
    #print( runID_list, len(runID_list) )
    #print( root_numpy.root2tree(file, treename = 'hitSumm').Scan() )
    
    #with Timer('number of entries') as t:
        #number_of_entries = get_number_of_entries( file )
        #print( number_of_entries )
    #with Timer('read catalog on tenth') as t:
        #data = read_catalog( file, indices = [0,number_of_entries/100] )
    #print( data[0] )
    #print( data.dtype )
    #print( data['xPix'].dtype )
    
    ##data['xPix'].dtype = np.dtype('<i4')
    ##data['xPix'] = data['xPix'].astype( np.dtype('<i4') )
    #print( data['xPix'].dtype )
    #print( data[10]['xPix'] )
    #lens = [ len(data['ePix'][i]) for i in range(len(data)) ]
    #print( min(lens), max(lens) )

    #output = 'test.root'
    #save_catalog( output, data )
    #print( get_number_of_entries(output) )
    #data = read_catalog( output )
    #print( list_branches(output) )
    #print( data[0] )

#def build_new_root_from_existing( input, output, condition = None, selection = '' ):
    #with Timer('open file') as t:
        #input_file = TFile( input )

    #with Timer('number_of_selected_hits') as t:
        #number_of_selected_hits = input_file.hitSumm.GetEntries( selection )
        #print( number_of_selected_hits )
    #if number_of_selected_hits == 0:
        #print( 'selection found no hits, exiting' )
        #return
    
    #output_file = TFile.Open( output, 'recreate' )
    #with Timer('clone config and write') as t:
        #output_config = input_file.config.CloneTree()
        #output_config.Write()

    #with Timer('clone tree') as t:
        #output_hitSumm = input_file.hitSumm.CloneTree(0)
        
    #with Timer('copy tree selection') as t:
        #cut_hitSumm = input_file.hitSumm.CopyTree( selection )

    #with Timer('select') as t:
        #print( cut_hitSumm.GetEntries() )
    
    #if condition:
        #def func_event( event ):
            #if condition(event):
                #output_hitSumm.Fill()            
        #with Timer('copy loop') as t:
            #map( func_event, cut_hitSumm )
    #else:
        #output_hitSumm = cut_hitSumm.CloneTree()

    #print( 'output_file hitSumm entries', output_hitSumm.GetEntries() )
    #with Timer('write and close') as t:    
        #output_hitSumm.Write()
        #output_file.Close()
        #input_file.Close()
    #return

#def build_new_catalog_from_recarray( input, output, model ):
    #N = len(input)
    #print( 'build_new_catalog_from_recarray', N )
    #class Event( TreeModel ):
        #flag = IntCol()
        #nSavedPix = IntCol()
        #runID = IntCol()
        #ohdu = IntCol()
        #nSat = IntCol()
        
        #ePix = FloatArrayCol(10000, length_name = 'nSavedPix')
        #xPix = IntArrayCol(10000, length_name = 'nSavedPix')
        #yPix = IntArrayCol(10000, length_name = 'nSavedPix')
        #level = IntArrayCol(10000, length_name = 'nSavedPix')
    
    #extra_fields = []
    #for name in input.dtype.names:
        #if not name in vars(Event):
            #type_ = getattr( input, name ).dtype
            #eventtype = None
            #if type_ == 'float64':
                #eventtype = FloatCol()
            #elif type_ == 'int64':
                #eventtype = FloatCol()
            #else:
                #raise Exception( 'type not recognized', type_ )
            #setattr( Event, name, eventtype )
            #extra_fields.append( name )
    
    #print( 'fields to be written in catalog', ['flag', 'runID', 'ohdu', 'nSat', 'nSavedPix', 'ePix', 'xPix', 'yPix', 'level'] + extra_fields )

    #rfile = root_open( output, 'w' )
    #tree = Tree( name = 'hitSumm', model = Event )

    #for i in range(N):
        #tree.flag = 0
        #tree.nSat = 0
        #tree.runID = 0
        #tree.ohdu = 0
        #tree.nSavedPix = input.nSavedPix[i]
        
        #for j in range( tree.nSavedPix ):
            #tree.ePix[j] = input.ePix[i][j]
            #tree.xPix[j] = input.xPix[i][j]
            #tree.yPix[j] = input.yPix[i][j]
            #tree.level[j] = int( input.level[i][j] )
        #for name in extra_fields:
            #setattr( tree, name, getattr( input, name)[i] )
        #tree.fill()

    #tree.write()
    #rfile.close()

#def build_new_root_from_existing_rootpy( input, output, condition = None, selection = '' ):
    #with Timer('open file') as t:
        #input_file = rootpy.io.root_open( input )

    #output_file = rootpy.io.root_open( output, 'recreate' )
    #with Timer('clone config and write') as t:
        #output_config = input_file.config.CopyTree()
        #output_config.Write()

    #with Timer('get entries') as t:
        #print( input_file.hitSumm.GetEntries( selection ) )
        
    #with Timer('clone tree') as t:
        #output_hitSumm = input_file.hitSumm.Clone()
        
    #with Timer('copy tree selection') as t:
        #cut_hitSumm = input_file.hitSumm.CopyTree( selection )

    #with Timer('select') as t:
        #print( cut_hitSumm.GetEntries() )
    
    #if condition:
        #def func_event( event ):
            #if condition(event):
                #output_hitSumm.Fill()            
        #with Timer('copy loop') as t:
            #map( func_event, cut_hitSumm )
    #else:
        #output_hitSumm = cut_hitSumm.CloneTree()

    #print( 'output_file hitSumm entries', output_hitSumm.GetEntries() )
    #output_hitSumm.Write()
    #output_file.Close()
    #input_file.Close()

#def status( input ):
    #ifile = TFile( input )
    #print( 'trees', [ key.GetName() for key in ifile.GetListOfKeys() ] )
    #print( 'branches', map( lambda _: _.GetName(), ifile.hitSumm.GetListOfBranches() ) )
    #print( 'number of entries', ifile.hitSumm.GetEntries() )
    #print( ifile.hitSumm.Scan() )
    #exit(0)

#def open_catalog( filename ):
    ##tfile = TFile.Open( filename )
    ##print( type(tfile) )
    #catalog = Catalog( filename )
    #print( 'type', type(catalog) )
    #return catalog

##class Catalog(TFile):
    ##def __init__( self, filename ):
        ##super(Catalog, self).__init__(filename)
        ###self = TFile( filename )
        ##self.filename = filename
        ##self.list_of_treenames = [ tree.GetName() for tree in self.GetListOfKeys() ]
        ###print( '__init__ trees', self.list_of_treenames )
        ##for treename in self.list_of_treenames:
            ##tree = Tree( getattr(self, treename ) )
            ##setattr( self, treename, tree )
            ###self.__dict__[treename] = tree
        ##print( type( self.config ) )
    ##def get_tree( self, treename ):
        ##return Tree( self.Get(treename) )
    
    ##def __enter__( self ):
        ##return self
    
    ##def __exit__(self, type, value, traceback):
        ##print( 'exit called' )
        ###self.Close()

    #def status(self):
        #print( 'status', type(self) )
        #for treename in self.list_of_treenames:
            #print( 'in tree "%s"' % treename )
            #print( 'branches: ', ', '.join( getattr( self, treename ).list_of_branches ) )
            #print( 'number of entries', getattr( self, treename ).number_of_entries )
                
    #def cut(self, selection, treename = None ):
        #if treename in self.treename_list:
            #return getattr(self, treename).CopyTree( selection )
        #raise Exception('tree not selected')


    ##file1 = "/share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root"
    ##file5 = "/share/storage2/connie/DAna/Catalogs/hpix_scn_osi_raw_gain_catalog_data_9096_to_9295_v3.0.root"
    ##build_new_root_from_existing( file1, 'test.root', selection = '&&'.join(conds1_x) )

def extract( args ):
    args.input_files = args.fits_files
    args.ohdu = [3]
    def image_extract( image, args ):
        ret = image.data.extract_hits( args.threshold, args.border )
        print( image.ohdu, len(ret) )
        return ret
    args.func = image_extract
    image = Image.apply_to_files( args )
    print( 'returns' )
    for file_group, entry in image.items():
        for ohdu, value in entry.items():
            print( ohdu, value )
    #clusters = image.extract_hits( args )
    
def add_extract_options(p):
    p.add_argument('fits_files', type=str, help = 'fits file (example: )' )
    p.add_argument('-t', '--threshold', type=float, default=60, help = 'energy threshold in ADU' )
    p.add_argument('-b', '--border', type=int, default=0, help = 'pixel border' )
    #p.add_argument('branches', nargs=2, type=str, default= ['E0','xVar0'], help = 'branches used for x- and y-axis' )
    #p.add_argument('--selection', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'output to file' )
    p.set_defaults(func=extract)


def get_selection( file, branches, selection, start=None, stop=None ):
    return root_numpy.root2array( file, treename = 'hitSumm', branches = list(set(branches)), selection = selection, start=start, stop=stop )

def get_start_stop_from_runID_range( file, runID_range ):
    runIDs = get_selection( file, branches=['runID'], selection=None )['runID']
    start = np.amin(np.argwhere(runIDs==runID_range[0]))
    stop = np.amax(np.argwhere(runIDs==runID_range[1]))+1
    return start, stop

def get_selections( file, branches, selections, global_selection=None, start=None, stop=None, extra_branches = [] ):
    data_selection = {}
    lims = None
    for selection in selections:
        selection_full = selection
        if global_selection:
            selection_full = global_selection + ' and ' + selection
        selection_string = selection_full.replace('and', '&&')
        data_selection[selection] = get_selection( file, branches=branches+extra_branches, selection = 'flag==0 && ' + selection_string, start=start, stop=stop )
        if lims is None:
            lims = [ ( np.min(data_selection[selection][branch]), np.max(data_selection[selection][branch]) ) for branch in branches ]
        print( selection, data_selection[selection].shape )
    return data_selection, lims
    
def scatter( args ):
    with Timer('scatter'):
        file = glob.glob(args.root_file)
        start, stop = None, None
        if 'runID_range' in args:
            start, stop = get_start_stop_from_runID_range(file[0], args.runID_range)
        data_selection, lims = get_selections( file[0], args.branches, args.selections, args.global_selection, start=start, stop=stop, extra_branches=['ePix','xPix','yPix','xVar1'] )
                
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if 'global_selection' in args:
            ax.set_title( args.global_selection )
        if 'selections' in args:
            markers = ['.', '+', 'x', '^']
            for selection, marker in zip(args.selections, markers):
                datum = data_selection[selection]
                above_2sigma = datum['ePix']>30
                bary30 = np.average( (datum['xPix'][above_2sigma], datum['yPix'][above_2sigma]), weights = datum['ePix'][above_2sigma] )
                m2 = np.average( (datum['xPix'][above_2sigma]**2, datum['yPix'][above_2sigma]**2), weights = datum['ePix'][above_2sigma] )
                Var30 = m2 - bary30**2
                ax.scatter( datum[args.branches[0]], datum[args.branches[1]], marker=marker, alpha=.2, label=selection.replace('&&', 'and') )
                ax.scatter( datum[args.branches[0]], var30, marker=marker, alpha=.2, label=selection.replace('&&', 'and') )
        ax.legend()
        ax.set_xlim( lims[0] )
        ax.set_ylim( lims[1] )
        ax.set_xlabel(args.branches[0])
        ax.set_ylabel(args.branches[1])
    if 'ylog' in args:
        ax.set_yscale('log')
    if 'output' in args:
        fig.savefig(args.output+'.pdf')
    else:
        plt.show()
    return


def add_scatter_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('--branches', nargs=2, type=str, default= ['E0','xVar0'], help = 'branches used for x- and y-axis' )
    p.add_argument('--global-selection', type=str, default='0', help = 'global selection' )
    p.add_argument('--selections', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'output to file' )
    p.set_defaults(func=scatter)


def histogram( args ):
    with Timer('histogram'):
        file = glob.glob(args.root_file)
        data = root_numpy.root2array( file[0], treename = 'hitSumm', branches = args.branch, selection = 'flag==0' )
        if 'binsize' in args:
            bins = np.arange( np.min(data), np.max(data), args.binsize )
        elif 'nbins' in args:
            bins = np.linspace( np.min(data), np.max(data), args.nbins )
        else:
            bins = np.linspace( np.min(data), np.max(data), int(np.sqrt(len(data))) )
        
        data_selection = {}
        if 'selection' in args:
            bins = None
            for selection in args.selection:
                selection_string = selection.replace('and', '&&')
                data_selection[selection] = root_numpy.root2array( file[0], treename = 'hitSumm', branches = args.branch, selection = 'flag==0 && ' + selection_string )
                if bins is None:
                    if 'binsize' in args:
                        bins = np.arange( np.min(data_selection[selection]), np.max(data_selection[selection]), args.binsize )
                    elif 'nbins' in args:
                        bins = np.linspace( np.min(data_selection[selection]), np.max(data_selection[selection]), args.nbins )
                    else:    
                        bins = np.linspace( np.min(data_selection[selection]), np.max(data_selection[selection]), int(np.sqrt(len(data_selection[selection]))) )

        #print( 'bins', bins )
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(args.branch)
        ax.hist(data, bins=bins, histtype='step', label='all')
        if 'selection' in args:
            for selection in args.selection:
                ax.hist( data_selection[selection], bins=bins, histtype='step', label=selection.replace('&&', 'and') )
        ax.legend()
        ax.set_xlabel(args.branch)
        ax.set_yscale('log')
    if 'output' in args:
        fig.savefig(args.output+'.pdf')
    else:
        plt.show()
    return

def add_histogram_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('branch', type=str, default= 'E0', help = 'branch used for x-axis' )
    p.add_argument('--selection', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--define', type=str, default=argparse.SUPPRESS, help = 'definitions (ex.: a=E0; b=E1)' )
    p.add_argument('--binsize', type=float, default=argparse.SUPPRESS, help = 'binsize' )
    p.add_argument('--nbins', type=int, default=argparse.SUPPRESS, help = 'number of bins' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'selection' )
    
    p.set_defaults(func=histogram)

def status( args ):
    if not os.path.exists( args.root_file ):
        print( 'file {} does not exist'.format( arg.rootfile ) )
    tfile = TFile( args.root_file )
    for key in tfile.GetListOfKeys():
        tree_name = key.GetName()
        print( colored( tree_name, 'green' ) )
        tree = getattr(tfile, tree_name)
        branches = ', '.join( map( lambda _: _.GetName(), tree.GetListOfBranches() ) )
        print( branches )
        print( tree.GetEntries() )
    return

def add_status_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.set_defaults( func = status )

if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'catalog tools', formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    subparsers = parser.add_subparsers( help = 'major options' )

    add_status_options( subparsers.add_parser('status', help='status', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_histogram_options( subparsers.add_parser('histogram', help='histogram', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_scatter_options( subparsers.add_parser('scatter', help='scatter', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_extract_options( subparsers.add_parser('extract', help='extract', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    args = parser.parse_args()

    func = args.func
    del args.func
    print( colored('using parameters:', 'green', attrs=['bold'] ) )
    print_var( vars(args).keys(), vars(args), line_char='\t' )
    
    with Timer('finished'):
        func(args)
    
