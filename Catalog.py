#!/usr/bin/python

from __future__ import print_function
import os, sys, argparse, re, glob

from Timer import Timer
import warnings
warnings.filterwarnings("ignore")

try:
    import root_numpy
except ImportError:
    print('missing module, please run')
    print('module load softwares/root/5.34-gnu-5.3')
    exit(0)

import weave
from array import array
from numpy.lib import recfunctions as rfn
import matplotlib
matplotlib.use('gtk3agg')

import matplotlib.pylab as plt

from numpy.lib.recfunctions import append_fields, stack_arrays
import numpy as np
import scipy.special as special
import Statistics as stats

with Timer('import'):
    from ROOT import TFile, TTree, AddressOf

    #import rootpy
    #import rootpy.io
    #from array import array
    #from rootpy.tree import Tree, TreeModel, IntCol, IntArrayCol, FloatCol, FloatArrayCol, CharArrayCol
    #from rootpy.io import root_open

import Image
import Statistics as stats
from termcolor import colored
from PrintVar import print_var

class HitSummary(np.recarray):
    def __new__( cls, recarray ):
        print( '__new__', recarray.xPix[1] )
        return recarray.view(cls)
            
    #def __new__( cls, list_of_entries, levels, threshold, border ):
        #print( 'new hits' )
        #dtype = [
            #('ePix', object),
            #('xPix', object),
            #('yPix', object),
            #('level', object),
            #]
        #obj = np.array( list_of_entries[1:], dtype = dtype ).view(np.recarray).view(cls)
        #obj = rfn.append_fields( obj, 
                                 #('nSavedPix'), 
                                 #np.array(map(len, obj.ePix)), 
                                 #dtypes=(int), 
                                 #asrecarray=True 
                                 #).view(cls)

        #obj.border = border
        #obj.threshold = threshold
        #obj.background = list_of_entries[0]
        #return obj
    
    def add_fields( self, names, types=None, arrays=None, func=None, verbose=False ):
        if verbose: print( 'added', names )
        if type(names) is str:
            if names == 'flag' or names == 'nSat':
                types = np.int32
                arrays = [0]*len(self)
            if names == 'nSavedPix':
                types = np.int32
                arrays = map(len, self.ePix)
            if re.match( r'^n[0-3]$', names ):
                lvl = int(re.search( r'^n([0-3])$', names ).groups()[0])
                types = np.int32
                arrays = self.compute_number_of_pixel_levels(lvl)
            if re.match( r'^E[0-3]$', names ):
                lvl = int(re.search( r'^E([0-3])$', names ).groups()[0])
                types = np.float32
                arrays = self.compute_energy_levels(lvl)
        return rfn.append_fields( self, names, arrays, dtypes=types, asrecarray=True ).view(HitSummary)
    
    def compute_number_of_pixel_levels( self, lvl ):
        return np.array(map( lambda entry: len(entry.ePix[ entry.level <= lvl ]), self ))
        
    def add_number_of_pixel_levels( self, lvls ):
        added_fields = []
        for lvl in lvls:
            new_field = 'n%d' % lvl
            self = self.add_fields( 'n%d' % lvl, int, self.compute_number_of_pixel_levels(lvl) )
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

def build_new_catalog_from_recarray_old( input, output ):
    N = len(input)
    print( 'build_new_catalog_from_recarray', N )
    Nmax = np.max( input.nSavedPix )
    print( 'Nmax', Nmax )
    class Event( TreeModel ):
        flag = IntCol()
        nSavedPix = IntCol()
        runID = IntCol()
        ohdu = IntCol()
        nSat = IntCol()
        
        ePix = FloatArrayCol(Nmax, length_name = 'nSavedPix')
        xPix = IntArrayCol(Nmax, length_name = 'nSavedPix')
        yPix = IntArrayCol(Nmax, length_name = 'nSavedPix')
        level = IntArrayCol(Nmax, length_name = 'nSavedPix')
    
    extra_fields = []
    for name in input.dtype.names:
        if not name in vars(Event):
            type_ = getattr( input, name ).dtype
            eventtype = None
            if type_ == 'float64':
                eventtype = FloatCol()
            elif type_ == 'int64':
                eventtype = FloatCol()
            else:
                raise Exception( 'type not recognized', type_ )
            setattr( Event, name, eventtype )
            extra_fields.append( name )
    
    print( 'fields to be written in catalog', ['flag', 'runID', 'ohdu', 'nSat', 'nSavedPix', 'ePix', 'xPix', 'yPix', 'level'] + extra_fields )

    rfile = root_open( output + '.root', 'w' )
    tree = Tree( name = 'hitSumm', model = Event )

    with Timer('fill') as t:
        N = len(input)
        def fill(entry):
            t.check( 30, N )
            tree.flag = entry.flag
            tree.nSat = entry.nSat
            tree.runID = entry.runID
            tree.ohdu = entry.ohdu
            tree.nSavedPix = entry.nSavedPix

            tree.ePix = entry.ePix
            tree.xPix = entry.xPix
            tree.yPix = entry.yPix
            tree.level = entry.level.astype(int)
            [ setattr( tree, name, getattr( entry, name) ) for name in extra_fields ]

            tree.fill()
        map( fill, input )

    tree.write()
    rfile.close()

def build_new_catalog_from_recarray2( input, output ):
    N = len(input)
    print( 'build_new_catalog_from_recarray2', N )
    Nmax = np.max( input.nSavedPix )
    print( 'Nmax', Nmax )
    #print( 'fields to be written in catalog', ['flag', 'runID', 'ohdu', 'nSat', 'nSavedPix', 'ePix', 'xPix', 'yPix', 'level'] + extra_fields )

    root_file = TFile( output + '.root', 'recreate' )
    tree = TTree( 'hitSumm', 'hitSumm' )
    
    class Model:
        pass
    model = Model()
    
    name = 'nSavedPix'
    setattr( model, name, array( 'i', [0] ) )
    expr = '{}/I'.format(name)
    tree.Branch( name, getattr(model, name), expr )
    for name, type in input.dtype.fields.items():
        print( name, type, type[0] )
        if name == 'nSavedPix':
            continue
        elif name == 'level':
            print( 'level' )
            model.level = array( 'i', [0]*Nmax )
            expr = 'level[nSavedPix]/I'
        elif type[0] == np.float64:
            print( 'float', name )
            setattr( model, name, array( 'f', [0.] ) )
            expr = '{}/F'.format(name)
        elif type[0] == np.int64:
            setattr( model, name, array( 'i', [0] ) )
            expr = '{}/I'.format(name)
        elif type[0] == object:
            print( 'object', name )
            setattr( model, name, array( 'f', [0.]*Nmax ) )
            expr = '{}[nSavedPix]/F'.format(name)
        tree.Branch( name, getattr(model, name), expr )
            
    with Timer('fill') as t:
        N = len(input)
        def fill_tree(entry):
            t.check( 30, N )
            for field in input.dtype.names:
                setattr( model, field, getattr(entry, field) )
            tree.Fill()
        map( fill_tree, input )

    tree.Write()
    root_file.Close()

def build_new_catalog_from_recarray__( input, output ):
    N = len(input)
    print( 'build_new_catalog_from_recarray2', N )
    Nmax = np.max( input.nSavedPix )
    print( 'Nmax', Nmax )
    with Timer('fill') as t:
        root_numpy.array2root( input, output + '.root', treename='hitSumm', mode='recreate' )

    #for fieldname in ['ePix', 'xPix', 'yPix', 'level']:
        #field = np.zeros( (len(input)), dtype=[(fieldname, np.float32, (Nmax,))] ).view(np.recarray)
        #for entry, f in zip(input,field):
            #getattr(f,fieldname)[:len(getattr(entry,fieldname))] = getattr(entry,fieldname)[:len(getattr(entry,fieldname))]
    
        #root_numpy.array2root( field, output+'.root', treename='hitSumm', mode='update' )

def build_new_catalog_from_recarray( input, output ):
    N = len(input)
    print( input.xPix[1] )
    print( 'build_new_catalog_from_recarray2', N )
    Nmax = np.max( input.nSavedPix )
    print( 'Nmax', Nmax )
    
    declarations = 'const Int_t knSavedPix = {};\n'.format(Nmax)
    assignments = ''
    assignmentsArray = ''
    branch = ''
    for name in input.dtype.names[::-1]:
        type_ = getattr( input, name ).dtype
        if name.startswith('n') and name != 'nSavedPix' and name != 'nSat':
            field = name
            var = '&c_{}'.format(name)
            expr = '{}/F'.format(name)
            declarations += 'Float_t c_{};\n'.format(name)
            assignments += '\tc_{} = {}[n];\n'.format( name, name )
        elif type_ == np.int32:
            field = name
            var = '&c_{}'.format(name)
            expr = '{}/I'.format(name)
            declarations += 'Int_t c_{};\n'.format(name)
            assignments += '\tc_{} = {}[n];\n'.format( name, name )
        elif type_ == np.float32:
            field = name
            var = '&c_{}'.format(name)
            expr = '{}/F'.format(name)
            declarations += 'Float_t c_{};\n'.format(name)
            assignments += '\tc_{} = {}[n];\n'.format( name, name )
        elif name in ['level','xPix','yPix']:
            field = '{}[c_nSavedPix]'.format(name)
            field = name
            var = 'c_{}'.format(name)
            expr = '{}[nSavedPix]/I'.format(name)
            declarations += 'Int_t c_{}[knSavedPix];\n'.format(name)
            assignmentsArray += '\t\tc_{}[i] = ((Long_t*) PyArray_DATA(PyList_GetItem({},n))) [i];\n'.format( name, name )
        elif name == 'ePix':
            field = '{}[c_nSavedPix]'.format(name)
            field = name
            var = 'c_{}'.format(name)
            expr = '{}[nSavedPix]/F'.format(name)
            declarations += 'Float_t c_{}[knSavedPix];\n'.format(name)
            assignmentsArray += '\t\tc_{}[i] = ((Double_t*) PyArray_DATA(PyList_GetItem({},n))) [i];\n'.format( name, name )
        branch += 'tree->Branch( "{}", {}, "{}" );\n'.format( field, var, expr )
    
    code = declarations + '\n';
    #code += 'PyArrayObject *ePix_array = convert_to_numpy(PyList_GetItem(ePix,10), "ePix");\n'
    #code += 'conversion_numpy_check_type(ePix_array, PyArray_DOUBLE, "ePix");'
    #code += 'for(int i=0; i<nSavedPix[10]; ++i){\n' 
    #code += '\tstd::cout<<" "<< ((Double_t *) ePix_array->data)[i] <<" ";\n'
    #code += '}\n'
    #code += 'std::cout<<std::endl;\n'    
    #code += 'for(int i=0; i<nSavedPix[10]; ++i){\n' 
    #code += '\tstd::cout<<" "<< ((Long_t *) PyArray_DATA(PyList_GetItem(level,10)))[i] <<" ";\n'
    #code += '}\n'
    #code += 'std::cout<<std::endl;\n'    
    #code += 'for(int i=0; i<nSavedPix[10]; ++i){\n' 
    #code += '\tstd::cout<<" "<< ((Long_t *) PyArray_DATA(PyList_GetItem(xPix,10)))[i] <<" ";\n'
    #code += '}\n'
    #code += 'std::cout<<std::endl;\n'    
    #code += 'for(int i=0; i<nSavedPix[10]; ++i){\n' 
    #code += '\tstd::cout<<" "<< ((Long_t *) PyArray_DATA(PyList_GetItem(yPix,10)))[i] <<" ";\n'
    #code += '}\n'
    #code += 'std::cout<<std::endl;\n'    
    #code += 'std::cout<<"this was ePix[10]"<<std::endl;\n'
    code += 'TFile f( "{}.root", "recreate" );\n'.format( output )
    code += 'TTree *tree = new TTree( "hitSumm", "hitSumm" );\n'
    code += '\n' + branch + '\n';
    #code += 'for(int n=0; n<{}; ++n){{ std::cout<<n<<" "<<(Float_t) nSavedPix[n]<<std::endl; }}\n'.format(len(input))
    code += 'for(int n=0; n<{}; ++n){{\n'.format( len(input) )
    #code += '\tstd::cout<<"n "<<n<<" "<<nSavedPix[n]<<std::endl;\n'
    code += assignments
    code += '\tfor(int i=0; i<nSavedPix[n]; ++i){\n'
    #code += '\t\tstd::cout<<"innerloop"<<nSavedPix[n]<<" "<<i<<std::endl;\n'
    code += assignmentsArray
    code += '\t}\n'
    code += '\ttree->Fill();\n'
    code += '}\n'
    code += '\ntree->Write();\n'

    #print( code )
    
    for name in input.dtype.names:
        exec '{} = np.array(input.{})'.format(name, name)
    #print( E1, type(E1), E1.dtype, E1.shape )
    #print( nSavedPix )

    #fields = [ '{}'.format(name) for name in input.dtype.names]
    #print( fields )

    
    ePix = [ np.array(e) for e in input.ePix ]
    xPix = [ np.array(x) for x in input.xPix ]
    yPix = [ np.array(y) for y in input.yPix ]
    level = [ np.array(l).astype(int) for l in input.level ]
    #for fieldname in ['ePix', 'xPix', 'yPix', 'level']:
    
    #print( nSavedPix )
    weave.inline( code, input.dtype.names,
                        #type_converters=weave.converters.blitz,
                        headers=['"TFile.h"', '"TTree.h"'],
                        libraries=['Core'],
                        include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                        library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                        #extra_compile_args=['-O3'],
                        verbose=1,
                        )    
        
    
    
    
    
    
    
    
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
    if os.path.exists( args.name +'.root' ):
        print( 'catalog {}.root already exists'.format(args.name) )
        exit(1)
    args.input_files = args.fits_files
    #args.ohdu = [3,4]
    def image_extract( image, args ):
        hitSummary = image.data.extract_hits( args.threshold, args.border ).view(HitSummary)
        for field, dtype in [('ohdu', np.int32), ('runID', np.int32), ('gain', np.float32), ('rebin',np.int32), ('dc',np.float32), ('noise',np.float32)]:
            if field.upper() in image.header:
                hitSummary = hitSummary.add_fields(field, dtype, [image.header[field.upper()]]*len(hitSummary) )
        print( 'runID {} ohdu {} extracted hits {}'.format( int(image.header['RUNID']), image.header['OHDU'], len(hitSummary) ) )
        hitSummary = hitSummary.add_fields( 'thr', np.float32, [args.threshold]*len(hitSummary) )
        for field in ['flag', 'nSavedPix', 'nSat' ] + [ fmt.format(lvl) for lvl in range(args.border+1) for fmt in ['n{}', 'E{}'] ]:
            hitSummary = hitSummary.add_fields( field )
        return hitSummary
    args.func = image_extract
    image = Image.apply_to_files( args )
    
    hitSummaries = None
    for file_group, entry in image.items():
        for ohdu, hitSummary in entry.items():
            if hitSummaries is None:
                hitSummaries = hitSummary
            else:
                hitSummaries = rfn.stack_arrays( (hitSummaries, hitSummary), asrecarray=True )
    #print( hitSummaries.dtype, hitSummaries.ohdu, hitSummaries.runID, len(hitSummaries), len(hitSummaries.ePix) )
    build_new_catalog_from_recarray( hitSummaries, args.name )
    #clusters = image.extract_hits( args )
    
def add_extract_options(p):
    p.add_argument('name', type=str, help = 'basename of the output catalog' )
    p.add_argument('fits_files', type=str, help = 'fits file (example: )' )
    p.add_argument('-t', '--threshold', type=float, default=60, help = 'energy threshold in ADU' )
    p.add_argument('-b', '--border', type=int, default=0, help = 'pixel border' )
    #p.add_argument('branches', nargs=2, type=str, default= ['E0','xVar0'], help = 'branches used for x- and y-axis' )
    p.add_argument('--ohdu', nargs='+', type=int, default=range(2,16), help = 'ohdus to be extracted' )
    p.add_argument('--exclude', nargs='+', type=int, default=[11,12,14], help = 'ohdus NOT to be extracted' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'output to file' )
    p.set_defaults(func=extract)


def get_selection( file, branches, selection, start=None, stop=None ):
    return root_numpy.root2array( file, treename = 'hitSumm', branches = list(set(branches)), selection = selection, start=start, stop=stop ).view(np.recarray)

def get_start_stop_from_runID_range( file, runID_range ):
    runIDs = get_selection( file, branches=['runID'], selection=None ).runID
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

def mean_norm( mu, sigma, a, b ):
    primitive = lambda x: .5*mu*special.erf((x - mu)/(np.sqrt(2)*sigma)) - sigma*np.exp(-(x - mu)**2/(2*sigma**2))/np.sqrt(2*np.pi)
    if a == -np.inf and b == np.inf:
        return mu
    elif a == -np.inf:
        return primitive(b) + .5*mu
    elif b == np.inf:
        return .5*mu - primitive(a)
    return primitive(b) - primitive(a)

def var_norm( mu, sigma, a, b ):
    primitive = lambda x: .5*(mu**2 + sigma**2) *special.erf((x - mu)/(np.sqrt(2)*sigma)) \
        - sigma*(mu+x)*np.exp(-(x - mu)**2/(2 *sigma**2))/np.sqrt(2*np.pi)
    if a == -np.inf and b == np.inf:
        return mu
    elif a == -np.inf:
        return primitive(b) + .5*(mu**2 + sigma**2)
    elif b == np.inf:
        return .5*(mu**2 + sigma**2) - primitive(a)
    return primitive(b) - primitive(a)

def sum_norm( mu, sigma, a, b ):
    try:
        print('shapes', a.shape, b.shape, mu.shape, sigma.shape )
    except:
        pass
    primitive = lambda x: .5* special.erf((x - mu)/(np.sqrt(2)*sigma))
    return primitive(b) - primitive(a)

def test_std_correction():
    print( 'mean_norm', mean_norm( 1, 2, -np.inf, np.inf ) )
    print( 'mean_norm', mean_norm( 1, 2, 0, np.inf ) )
    print( 'mean_norm', mean_norm( 1, 2, -np.inf, 0 ) )
    print( 'mean_norm', mean_norm( 1, 2, 0, 10 ) + mean_norm( 1, 2, -5, 0 ) )

    print( 'var_norm', 1**2/var_norm( 0, 1, -3, 3 ) )
    sigma = 3.11
    rvs = stats.norm.rvs(0, sigma, int(1e3))
    p = stats.norm.pdf(rvs)
    f=0
    std = np.std(rvs)
    print( 'stats', std, sigma, np.sum(rvs) )
    rvs = rvs[p>np.max(p)*f]
    print( 'len', len(rvs) )
    std = np.std(rvs)
    mu = np.mean(rvs)
    rvsmin, rvsmax = np.min(rvs), np.max(rvs) 
    print( 'stats factor', mu, std, sigma, std*np.sqrt(std**2/var_norm( mu, std, rvsmin, rvsmax)) )
    print( 'stats add', mu, std, sigma, np.sqrt( (std**2 + var_norm( mu, std, -np.inf, rvsmin) + var_norm( mu, std, rvsmax, np.inf) ))*std**2/var_norm(mu, std, rvsmin, rvsmax) )
    print( 'stats add', mu, std, sigma, np.sqrt( (std**2 + var_norm( mu, std, -np.inf, rvsmin) + var_norm( mu, std, rvsmax, np.inf) )) )
    print( 'stats', mu, std, sigma, std*std**2/var_norm( mu, std, rvsmin, rvsmax) )
    #print( 'stats fwhm', rvsmin, rvsmax, rvsmax-rvsmin, (rvsmax-rvsmin)/( 2*np.sqrt(2*np.log(1./f)) ), stats.FWHM1d(rvs, f=.1) )
    print( np.sum(rvs), sum_norm( mu, sigma, rvsmin, rvsmax ), np.sum(rvs)/sum_norm( mu, sigma, rvsmin, rvsmax ) )
    print( stats.norm.fit(rvs) )
    print( stats.norm.fit2(rvs) )
    y, x, dx = stats.make_histogram( rvs, bins=len(rvs)/20 )
    print( stats.norm.fit2( x, mu=mu, sigma=std, weights=y ) )
    print( 'sum', sum(y) )
    print( stats.norm.fit_curve( y, x, A=sum(y), mu=mu, sigma=std ), dx )
    return

#test_std_correction()

def scatter( args ):
    with Timer('scatter'):
        file = glob.glob(args.root_file)
        start, stop = None, None
        if 'runID_range' in args:
            start, stop = get_start_stop_from_runID_range(file[0], args.runID_range)
        data_selection, lims = get_selections( file[0], args.branches, args.selections, args.global_selection, start=start, stop=stop, extra_branches=['ePix','xPix','yPix','xVar1','E0', 'E1', 'E2', 'n0', 'n1', 'n2'] )
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if 'global_selection' in args:
            ax.set_title( args.global_selection )
        if 'selections' in args:
            markers = ['.', '+', 'x', '^']
            for selection, marker in zip(args.selections, markers):
                datum = data_selection[selection]
                above_2sigma = [ ePix>30 for ePix in datum.ePix ]
                moment = lambda n: np.array( [ np.average( (xPix[mask]**n, yPix[mask]**n), weights = ePix[mask], axis=1 ) for mask, xPix, yPix, ePix in zip(above_2sigma, datum.xPix, datum.yPix, datum.ePix) ] )
                m1 = moment(1)
                m2 = moment(2)
                var30 = (m2 - m1**2).T
                E30, n30 = np.array( [ (np.sum( ePix[mask] ), len(ePix[mask]) ) for mask, ePix in zip(above_2sigma, datum.ePix) ] ).T
                
                #Min = np.array( [ np.min( ( xPix[mask], yPix[mask] ), axis=1) for mask, xPix, yPix in zip(above_2sigma, datum.xPix, datum.yPix) ] ).T
                #Max = np.array( [ np.max( ( xPix[mask], yPix[mask] ), axis=1) for mask, xPix, yPix in zip(above_2sigma, datum.xPix, datum.yPix) ] ).T
                #print( 'min, max', Max[0]-Min[0] )
                #print(Min[0], var30)
                #E30c = E30/( sum_norm(mu=m1[:,0], sigma=var30[0], a=Min[0], b=Max[0]) * sum_norm(mu=m1[:,1], sigma=var30[0], a=Min[1], b=Max[1]) )
                #print( zip(E30c, E30) )
                Eerr = 15
                if False:
                    ax.scatter( datum[args.branches[0]], datum[args.branches[1]], marker=marker, alpha=.2, label=selection.replace('&&', 'and') )
                ax.errorbar( datum['E0'], datum[args.branches[1]], xerr=Eerr*datum['n0'], label='E0', fmt=markers[0] )
                ax.errorbar( datum['E1'], datum[args.branches[1]], xerr=Eerr*datum['n1'], label='E1', fmt=markers[1] )
                #ax.errorbar( datum['E2'], datum[args.branches[1]], xerr=Eerr*datum['n2'], label='E3', fmt=markers[2] )
                ax.errorbar( E30, var30[0], xerr=Eerr*n30, label='E30', fmt=markers[3] )
                #ax.scatter( E30c, var30[0], marker=marker, alpha=.2, label=selection.replace('&&', 'and') )
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
    p.add_argument('--define', nargs='+', type=str, default=argparse.SUPPRESS, help = 'definitions' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'output to file' )
    p.set_defaults(func=scatter)

def energy_threshold( datum, threshold ):
    above_2sigma = [ ePix>threshold for ePix in datum.ePix ]
    moment = lambda n: np.array( [ np.average( (xPix[mask]**n, yPix[mask]**n), weights = ePix[mask], axis=1 ) for mask, xPix, yPix, ePix in zip(above_2sigma, datum.xPix, datum.yPix, datum.ePix) ] )
    m1 = moment(1)
    m2 = moment(2)
    var = (m2 - m1**2).T
    E, n = np.array( [ (np.sum( ePix[mask] ), len(ePix[mask]) ) for mask, ePix in zip(above_2sigma, datum.ePix) ] ).T
    dtype = [ ('{}{}'.format(key,threshold), float) for key in ('E', 'n', 'xVar', 'yVar') ]
    #print( dtype )
    return np.array( zip(E, n, var[0], var[1]), dtype=dtype ).view(np.recarray)

def pdf( x, E, n, sigma ):
    return np.sum( [stats.norm.pdf(x, iE, np.sqrt(in_)*sigma) for iE, in_ in zip(E,n)], axis=0 )

def histogram( args ):
    with Timer('histogram'):
        file = glob.glob(args.root_file)
        data = root_numpy.root2array( file[0], treename = 'hitSumm', branches = args.branch, selection = 'flag==0' )
        #if 'binsize' in args:
            #bins = np.arange( np.min(data), np.max(data), args.binsize )
        #elif 'nbins' in args:
            #bins = np.linspace( np.min(data), np.max(data), args.nbins )
        #else:
            #bins = np.linspace( np.min(data), np.max(data), int(np.sqrt(len(data))) )

        start, stop = None, None
        if 'runID_range' in args:
            start, stop = get_start_stop_from_runID_range(file[0], args.runID_range)
        data_selection, lims = get_selections( file[0], [args.branch], args.selections, args.global_selection, start=start, stop=stop,
                                              extra_branches=[
                                                  #'ePix','xPix','yPix',
                                                  'E0', 'E1', 'n0', 'n1'] )

        first_data_selection = data_selection.values()[0][args.branch]
        if 'binsize' in args:
            bins = np.arange( np.min(first_data_selection), np.max(first_data_selection), args.binsize )
        elif 'nbins' in args:
            bins = np.linspace( np.min(first_data_selection), np.max(first_data_selection), args.nbins )
        else:
            bins = np.linspace( np.min(first_data_selection), np.max(first_data_selection), int(np.sqrt(len(first_data_selection))) )
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(args.branch)
        #ax.hist(data, bins=bins, histtype='step', label='all')
        if 'selections' in args:
            for selection in args.selections:
                print( 'selection', data_selection[selection][args.branch].shape, len(bins) )
                datum = data_selection[selection]

                sigma = 12.5
                #ax.step( bins, pdf( bins, datum.E0, datum.n0, sigma ), where='mid', label='E0')
                ax.step( bins, pdf( bins, datum.E1, datum.n1, sigma ), where='mid', label='E1')

                #ax.step( bins, pdf( bins, datum.E2, datum.n2, sigma ), where='mid', label='E2')

                #dat45 = energy_threshold( datum, 45)
                #ax.step( bins, pdf( bins, dat45.E45, dat45.n45, sigma ), where='mid', label='E45')
                #dat30 = energy_threshold( datum, 30)
                #ax.step( bins, pdf( bins, dat30.E30, dat30.n30, sigma ), where='mid', label='E30')

                #dat15 = energy_threshold( datum, 15)
                #ax.step( bins, pdf( bins, dat15.E15, dat15.n15, sigma ), where='mid', label='E15')
                
                #pdf30 = lambda x: np.sum( [stats.norm.pdf(x, E, n*Eerr) for E, n in zip(E30,n30)], axis=0 )
                #ax.step( bins, pdf30(bins), where='mid', label='E30')
                
                #E30 = E30[ np.logical_or(var30[0]<.3, var30[0]>.8) ]
                #n30 = n30[ np.logical_or(var30[0]<.3, var30[0]>.8) ]
                #pdf30 = lambda x: np.sum( [stats.norm.pdf(x, E, n*Eerr) for E, n in zip(E30,n30)], axis=0 )
                #ax.step( bins, pdf30(bins), where='mid', label='E30*')
        ax.legend()
        ax.set_xlabel(args.branch)
        #ax.set_yscale('log')
    if 'output' in args:
        fig.savefig(args.output+'.pdf')
    else:
        plt.show()
    return

def add_histogram_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('branch', type=str, default= 'E0', help = 'branch used for x-axis' )
    p.add_argument('--selections', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selections' )
    p.add_argument('--global-selection', type=str, default='0', help = 'global selection' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    p.add_argument('--energy-threshold', nargs='+', type=float, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--define', type=str, default=argparse.SUPPRESS, help = 'definitions (ex.: a=E0; b=E1)' )
    p.add_argument('--binsize', type=float, default=argparse.SUPPRESS, help = 'binsize' )
    p.add_argument('--nbins', type=int, default=argparse.SUPPRESS, help = 'number of bins' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'selection' )
    
    p.set_defaults(func=histogram)

def status( args ):
    if not os.path.exists( args.root_file ):
        print( 'file {} does not exist'.format( args.rootfile ) )
    tfile = TFile( args.root_file )
    for key in tfile.GetListOfKeys():
        tree_name = key.GetName()
        print( colored( tree_name, 'green' ) )
        tree = getattr(tfile, tree_name)
        branches = ', '.join( map( lambda _: _.GetName(), tree.GetListOfBranches() ) )
        print( branches )
        print( tree.GetEntries() )
    tfile.Close()
    if 'tree' in args:
        if 'branches' in args:
            data = root_numpy.root2array( args.root_file, treename = args.tree, branches = args.branches ).view(np.recarray)
            for branch in args.branches:
                print( 'branch', getattr( data, branch ) )
    return

def add_status_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('-t', '--tree', type=str, default=argparse.SUPPRESS, help = 'tree to print' )
    p.add_argument('--branches', nargs='+', type=str, default=argparse.SUPPRESS, help = 'branch used for x-axis' )
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
    
