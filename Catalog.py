#!/usr/bin/python

from __future__ import print_function
import os, sys, argparse

from Timer import Timer
#import warnings
#warnings.filterwarnings("ignore")

try:
    import root_numpy
except ImportError:
    print('missing module, please run')
    print('module load softwares/root/5.34-gnu-5.3')
    exit(0)

import matplotlib
matplotlib.use('gtk3agg')

import matplotlib.pylab as plt

from numpy.lib.recfunctions import append_fields, stack_arrays
import numpy as np

#from ROOT import TFile, TTree, AddressOf

#import rootpy
#import rootpy.io
#from array import array
#from rootpy.tree import Tree, TreeModel, IntCol, IntArrayCol, FloatCol, FloatArrayCol, CharArrayCol
#from rootpy.io import root_open


from termcolor import colored
from PrintVar import print_var

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

def histogram( args ):
    data = root_numpy.root2array( args.root_file, treename = 'hitSumm', branches = args.branch, selection = 'flag==0' )
    if 'selection' in args:
        data_selection = root_numpy.root2array( args.root_file, treename = 'hitSumm', branches = args.branch, selection = 'flag==0 && ' + args.selection )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(args.branch)
    bins = None
    if 'selection' in args:
        bins = int(np.sqrt(len(data_selection)))
        hist, bins, _ = ax.hist(data_selection, bins=bins, histtype='step', label=args.selection)
    else:
        bins = int(np.sqrt(len(data_selection)))
    hist, bins, _ = ax.hist(data, bins=bins, histtype='step', label='all')
    ax.legend()
    ax.set_yscale('log')
    plt.show()
    return

def add_histogram_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('branch', type=str, default= 'E0', help = 'branch used for x-axis' )
    p.add_argument('--selection', type=str, default = argparse.SUPPRESS, help = 'selection' )
    
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

    args = parser.parse_args()

    func = args.func
    del args.func
    print( colored('using parameters:', 'green', attrs=['bold'] ) )
    print_var( vars(args).keys(), vars(args), line_char='\t' )
    
    with Timer('finished'):
        func(args)
    
