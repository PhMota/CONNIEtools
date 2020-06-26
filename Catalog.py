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
#from array import array
from numpy.lib import recfunctions as rfn
import matplotlib
matplotlib.use('gtk3agg')

import matplotlib.pylab as plt

from numpy.lib.recfunctions import append_fields, stack_arrays
from numpy import *
import Statistics as stats

with Timer('import'):
    from ROOT import TFile, TTree, AddressOf

    #import rootpy
    #import rootpy.io
    #from array import array
    #from rootpy.tree import Tree, TreeModel, IntCol, IntArrayCol, FloatCol, FloatArrayCol, CharArrayCol
    #from rootpy.io import root_open

import Image
import Simulation
import Statistics as stats
from termcolor import colored
from PrintVar import print_var
from scipy.sparse import csr_matrix


def binned_statistic(x, values, func, nbins, range):
    '''The usage is nearly the same as scipy.stats.binned_statistic''' 

    N = len(values)
    r0, r1 = range

    digitized = (float(nbins)/(r1 - r0)*(x - r0)).astype(int)
    S = csr_matrix((values, [digitized, np.arange(N)]), shape=(nbins, N))

    return [func(group) for group in np.split(S.data, S.indptr[1:-1])]

class NeighborIndexHistogram:
    def __init__( self, data, length, shift=1 ):
        import itertools
        self.__data = data
        self.__d = len(self.__data.shape)
        self.shiftsdd = tuple( itertools.product( *[range( -shift, shift+1 )]*self.__d ) )
        self.digitized = rint(self.__data/length).astype(int)
        self.occupied_bins, inv_bins = unique( self.digitized, axis=0, return_inverse = True )
        
        print( 'dig', self.digitized, len(self.digitized) )
        print( 'bins', self.occupied_bins, len(self.occupied_bins) )
        print( 'inv_bins',inv_bins, len(inv_bins) )
        print( 'data', data, len(data) )
        print( 'occupied', self.occupied_bins )
        
        for inv_bin, occupied in enumerate(self.occupied_bins):
            elements = nonzero( inv_bins == inv_bin )[0]
            print( 'elements', elements )
        

class HitSummary:
    def __init__( self, recarray, transpose=False, verbose=False ):
        self.verbose = verbose
        if transpose:
            recarray.xPix, recarray.yPix = recarray.yPix[:], recarray.xPix[:]
        self.__updateattrs__(recarray)

    def __updateattrs__(self, recarray ):
        self.__recarray__ = recarray
        for name in self.__recarray__.dtype.names:
            setattr(self, name, getattr(self.__recarray__, name) )
        
    def __getitem__(self, *args):
        return self.__recarray__[args]

    def __len__(self):
        return len(self.__recarray__)
    
    #def __new__( cls, recarray ):
        #print( '__new__', recarray.xPix[1] )
        #obj = recarray.view(cls)
        #return obj
    
    def extend(self, a):
        self.__updateattrs__(rfn.stack_arrays( (self.__recarray__, a), asrecarray=True ))
    
    def add_fields( self, names, types=None, arrays=None ):
        if self.verbose: print( 'adding', names )
        self.__updateattrs__( rfn.append_fields( self.__recarray__, names, arrays, dtypes=types, asrecarray=True ) )

    def __apply_to_entries( self, func, mask=None ):
        if mask is None:
            return array( map( func, self.__recarray__ ) )
        return array( map( lambda p: func(p[0], p[1]), zip(self.__recarray__, mask) ) )
    
    def __energy_mask(self, eThr):
        return map( lambda e: e >= eThr, self.ePix )

    def __level_mask(self, lvl):
        return map( lambda l: l <= lvl, self.level )

    def __flag( self, entry ): return 0
    
    def __nSat( self, entry ): return 0

    def __nSavedPix( self, entry ): return len(entry.ePix)
    
    def __n( self, entry, mask ): return len(entry.ePix[mask])
    
    def __E( self, entry, mask ): return sum(entry.ePix[mask])

    def __Bary(self, entry, mask, axis ):
        return average( getattr(entry, '{}Pix'.format(axis))[mask], weights=entry.ePix[mask], axis=0 )

    def __Var(self, entry, mask, axis ):
        return average( getattr(entry, '{}Pix'.format(axis))[mask]**2, weights=entry.ePix[mask], axis=0 ) \
                - average( getattr(entry, '{}Pix'.format(axis))[mask], weights=entry.ePix[mask], axis=0 )**2

    def __SqrBary(self, entry, mask, axis ):
        return average( getattr(entry, '{}Pix'.format(axis))[mask]**2, weights=entry.ePix[mask], axis=0 )
        
    def __thr_vars( self, entry, mask ):
        n = self.__n(entry, mask)
        E = self.__E(entry, mask)
        xBary = self.__Bary(entry, mask, 'x')
        yBary = self.__Bary(entry, mask, 'y')
        xVar = self.__SqrBary(entry, mask, 'x') - xBary**2
        yVar = self.__SqrBary(entry, mask, 'y') - yBary**2
        return n, E, xBary, yBary, xVar, yVar
    
    def __basic_vars( self, entry ):
        flag = 0
        nSat = 0
        nSavedPix = self.__nSavedPix(entry)
        return flag, nSat, nSavedPix
        
    def add_basic_fields( self ):
        self.add_fields( ['flag','nSat','nSavedPix'],
                        types=(int,int,int),
                        arrays=self.__apply_to_entries( self.__basic_vars ).T
                        )        
        
    def add_energy_threshold_fields( self, eThr ):
        self.add_fields( [ '{}{}'.format(name,int(eThr)) for name in ['n','E','xBary','yBary','xVar','yVar']],
                        types=(int,float,float,float,float,float),
                        arrays=self.__apply_to_entries( self.__thr_vars, mask=self.__energy_mask(eThr) ).T
                        )
        
    def add_level_fields( self, lvl ):
        self.add_fields( [ '{}{}'.format(name,int(lvl)) for name in ['n','E','xBary','yBary','xVar','yVar']],
                        types=(int,float,float,float,float,float),
                        arrays=self.__apply_to_entries( self.__thr_vars, mask=self.__level_mask(lvl) ).T
                        )
        
    def match_bruteforce( self, events, verbose, basename ):
        
        N = len(events)
        indices = range(N)
        
        if verbose:
            print( 'events.xy', events.xy )
        N_hits = len(self.__recarray__)
        x = zeros(N_hits)
        y = zeros(N_hits)
        z = zeros(N_hits)
        q = zeros(N_hits)
        E = zeros(N_hits)
        id_code = zeros(N_hits)
        sigma = zeros(N_hits)
        for i, hit in enumerate(self.__recarray__):
            xPix0 = hit.xPix[ hit.level==0 ]
            yPix0 = hit.yPix[ hit.level==0 ]
            if verbose:
                print( 'hitSummary.xy', zip(xPix0, yPix0) )
            for j, index in enumerate(indices):
                event = events[index]
                corr_y = event.y - events.vertical_overscan
                if verbose:
                    print( 'event.xy', event.x, event.y )
                if int(event.x) in xPix0 and int(corr_y) in yPix0:
                    if verbose:
                        print( 'match', event.x, event.y )
                    x[i] = event.x
                    y[i] = corr_y
                    z[i] = event.z
                    sigma[i] = event.sigma
                    q[i] = event.q_eff
                    E[i] = event.E_eff
                    id_code[i] = event.id_code
                    del indices[j]
                    break
        
        print( 'matched', count_nonzero(id_code) )
        code = '''
        TFile f( fname.c_str(), "update" );
        TTree *t = (TTree*) f.Get("hitSumm");
        Long64_t nentries = t->GetEntries();
        
        Float_t c_x;
        Float_t c_y;
        Float_t c_z;
        Int_t c_q;
        Float_t c_sigma;
        Float_t c_E;
        Int_t c_id_code;

        TBranch *x_branch = t->Branch("x", &c_x, "x/F");
        TBranch *y_branch = t->Branch("y", &c_y, "y/F");
        TBranch *z_branch = t->Branch("z", &c_z, "z/F");
        TBranch *q_branch = t->Branch("q", &c_q, "q/I");
        TBranch *sigma_branch = t->Branch("sigma", &c_sigma, "sigma/F");
        TBranch *E_branch = t->Branch("E", &c_E, "E/F");
        TBranch *id_code_branch = t->Branch("id_code", &c_id_code, "id_code/I");

        for( int n=0; n<nentries; ++n ){
            c_x = x[n];
            c_y = y[n];
            c_z = z[n];
            c_q = q[n];
            c_sigma = sigma[n];
            c_E = E[n];
            c_id_code = id_code[n];
            
            x_branch->Fill();
            y_branch->Fill();
            z_branch->Fill();
            q_branch->Fill();
            sigma_branch->Fill();
            E_branch->Fill();
            id_code_branch->Fill();
        }
        t->Write("", TObject::kOverwrite);
        '''
        fname = basename + '.root'
        varnames = ['x','y','z','q','sigma','E','id_code','fname']
        weave.inline( code, varnames,
                            headers=['"TFile.h"', '"TTree.h"', '"TObject.h"'],
                            libraries=['Core'],
                            include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                            library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                            extra_compile_args=['-O3', '-Wunused-variable'],
                            compiler='gcc',
                            #verbose=1,
                            )
        
        print( 'added columns {}'.format(varnames) )
        return
        
    def match( self, events, verbose, basename ):
        
        N = len(events)

        if verbose:
            print( 'events.xy', events.xy )
        
        xy2event_ind = {}
        for i, event in enumerate(events):
            key = ( int(event.x), int(event.y) )
            xy2event_ind[key] = i
        if verbose:
            print( xy2event_ind.keys() )
        
        N_hits = len(self.__recarray__)
        x = zeros(N_hits)
        y = zeros(N_hits)
        z = zeros(N_hits)
        q = zeros(N_hits)
        E = zeros(N_hits)
        id_code = zeros(N_hits)
        sigma = zeros(N_hits)

        for i, hit in enumerate(self.__recarray__):
            xPix0 = hit.xPix[ hit.level==0 ]
            yPix0 = hit.yPix[ hit.level==0 ]
            if verbose:
                print( 'hitSummary.xy', zip(xPix0, yPix0) )
            for xPix, yPix in zip(xPix0, yPix0):
                key = ( xPix, yPix )
                try:
                    event_ind = xy2event_ind[key]
                    event = events[event_ind]
                    if verbose:
                        print( 'match', event.x, event.y )
                    if id_code[i] == 0:
                        x[i] = event.x
                        y[i] = event.y
                        z[i] = event.z
                        sigma[i] = event.sigma
                        q[i] = event.q_eff
                        E[i] = event.E_eff
                        id_code[i] = event.id_code
                    else:
                        id_code[i] = 2
                        E[i] += event.E_eff
                        q[i] += event.q_eff
                        x[i] = -1
                        y[i] = -1
                        z[i] = -1
                        sigma[i] = -1
                except KeyError:
                    pass
        
        print( 'matched', count_nonzero(id_code) )
        code = '''
        TFile f( fname.c_str(), "update" );
        TTree *t = (TTree*) f.Get("hitSumm");
        Long64_t nentries = t->GetEntries();
        
        Float_t c_x;
        Float_t c_y;
        Float_t c_z;
        Int_t c_q;
        Float_t c_sigma;
        Float_t c_E;
        Int_t c_id_code;

        TBranch *x_branch = t->Branch("x", &c_x, "x/F");
        TBranch *y_branch = t->Branch("y", &c_y, "y/F");
        TBranch *z_branch = t->Branch("z", &c_z, "z/F");
        TBranch *q_branch = t->Branch("q", &c_q, "q/I");
        TBranch *sigma_branch = t->Branch("sigma", &c_sigma, "sigma/F");
        TBranch *E_branch = t->Branch("E", &c_E, "E/F");
        TBranch *id_code_branch = t->Branch("id_code", &c_id_code, "id_code/I");

        for( int n=0; n<nentries; ++n ){
            c_x = x[n];
            c_y = y[n];
            c_z = z[n];
            c_q = q[n];
            c_sigma = sigma[n];
            c_E = E[n];
            c_id_code = id_code[n];
            
            x_branch->Fill();
            y_branch->Fill();
            z_branch->Fill();
            q_branch->Fill();
            sigma_branch->Fill();
            E_branch->Fill();
            id_code_branch->Fill();
        }
        t->Write("", TObject::kOverwrite);
        '''
        fname = basename + '.root'
        varnames = ['x','y','z','q','sigma','E','id_code','fname']
        weave.inline( code, varnames,
                            headers=['"TFile.h"', '"TTree.h"', '"TObject.h"'],
                            libraries=['Core'],
                            include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                            library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                            extra_compile_args=['-O3', '-Wunused-variable'],
                            compiler='gcc',
                            #verbose=1,
                            )
        
        print( 'added columns {}'.format(varnames) )
        return
        
    def compute_sizelike_each( self, entry, lvl, gain, sigma_noise, fast = True, tol = 1e-1 ):
        _Q, _mu, _sigma = size_like( entry.ePix/gain, entry.xPix, entry.yPix, 
                                E0 = entry['E%d'%lvl]/gain, 
                                mu0 = [ entry['xBary%d' % lvl], entry['yBary%d' % lvl] ], 
                                sigma0 = sqrt( [ entry['xVar%d'%lvl], entry['yVar%d'%lvl] ] ), 
                                sigma_noise0 = sigma_noise,
                                single_sigma = False,
                                fast = fast,
                                tol = tol)
        return _Q*gain, _mu[0], _mu[1], _sigma[0], _sigma[1]

    def add_sizelike( self, lvl, gain, sigma_noise, fast = True, tol = 1e-1 ):
        new_fields = array( map( lambda entry: self.compute_sizelike_each( entry, lvl, gain, sigma_noise, fast, tol ), self ) )
        self = self.add_fields( ['EL', 'xMu', 'yMu', 'xSigma', 'ySigma'], new_fields.T, (float, float, float, float, float) )
        return self
        
    #def get_sizeLike( self, lvl, gain, sigma_noise, fast=True, tol = None ):
        #Q = []
        #mu, sigma = [], []
        #for entry in self:
            #_Q, _mu, _sigma = size_like( entry.ePix/gain, entry.xPix, entry.yPix, 
                                  #E0 = entry['E%d'%lvl]/gain, 
                                  #mu0 = [ entry['xBary%d' % lvl], entry['yBary%d' % lvl] ], 
                                  #sigma0 = sqrt( [ entry['xVar%d'%lvl], entry['yVar%d'%lvl] ] ), 
                                  #sigma_noise0 = sigma_noise,
                                  #single_sigma=False,
                                  #fast = fast,
                                  #tol = tol)
            #Q.append( _Q )
            #mu.append( _mu )
            #sigma.append( _sigma )
        #return array(Q)*gain, array(mu), array(sigma)
    
    def get_Bary(self):
        return array( (self.xBary, self.yBary) ).T

    def get_Var(self):
        return array( (self.xVar, self.yVar) ).T

    def save_catalog( self, basename ):
        N = len(self)
        Nmax = max( self.nSavedPix )
        if self.verbose:
            print( 'number of hits', N, 'nSavedPixmax', Nmax )
        declarations = ''
        assignments = ''
        assignmentsArray = ''
        branch = ''
        setbranch = ''
        for name in self.__recarray__.dtype.names[::-1]:
            type_ = getattr( self.__recarray__, name ).dtype
            if name.startswith('n') and name != 'nSavedPix' and name != 'nSat':
                var = '&c_{}'.format(name)
                expr = '{}/F'.format(name)
                declarations += 'Float_t c_{};\n'.format(name)
                assignments += '\tc_{} = {}[n];\n'.format( name, name )
            elif type_ == int32 or type_ == int:
                var = '&c_{}'.format(name)
                expr = '{}/I'.format(name)
                declarations += 'Int_t c_{};\n'.format(name)
                assignments += '\tc_{} = {}[n];\n'.format( name, name )
            elif type_ == float32 or type_ == float:
                var = '&c_{}'.format(name)
                expr = '{}/F'.format(name)
                declarations += 'Float_t c_{};\n'.format(name)
                assignments += '\tc_{} = {}[n];\n'.format( name, name )
            elif name in ['level','xPix','yPix']:
                var = 'c_{}'.format(name)
                expr = '{}[nSavedPix]/I'.format(name)
                declarations += 'Int_t c_{}[knSavedPix];\n'.format(name)
                assignmentsArray += '\t\tc_{}[i] = ((Long_t*) PyArray_DATA(PyList_GetItem({},n))) [i];\n'.format( name, name )
            elif name == 'ePix':
                var = 'c_{}'.format(name)
                expr = '{}[nSavedPix]/F'.format(name)
                declarations += 'Float_t c_{}[knSavedPix];\n'.format(name)
                assignmentsArray += '\t\tc_{}[i] = ((Double_t*) PyArray_DATA(PyList_GetItem({},n))) [i];\n'.format( name, name )
            else:
                print( 'not implemented', name, type_ )
                exit(0)
            branch += 'tree->Branch( "{}", {}, "{}" );\n'.format( name, var, expr )
            setbranch += 'tree->SetBranchAddress("{}", {});\n'.format(name, var)
        
        code = '''
        const Int_t knSavedPix = Nmax;
        
        TFile f( fname.c_str(), mode.c_str() );
        
        char c_mode[32];
        strcpy(c_mode, mode.c_str());
        TTree *tree;
        
        {declarations}
        
        if( strcmp(c_mode, "recreate") == 0 ){{
            std::cout<<"new TTree"<<std::endl;
            tree = new TTree( "hitSumm", "hitSumm" );
            {branch}
        }}else{{
            tree = (TTree*)f.Get("hitSumm");
            std::cout<<"nEntries "<<tree->GetEntries()<<std::endl;
            {setbranch}
        }}
        for(int n=0; n<N; ++n){{
            {assignments}
            for(int i=0; i<nSavedPix[n]; ++i){{
                {assignmentsArray}
            }}
            tree->Fill();
        }}
        tree->Write(0,TObject::kWriteDelete,0);
        '''.format(declarations=declarations, branch=branch, setbranch=setbranch, assignments=assignments, assignmentsArray=assignmentsArray)
        
        for name in self.__recarray__.dtype.names:
            exec '{} = array(self.__recarray__.{})'.format(name, name)

        ePix = [ array(e) for e in self.__recarray__.ePix ]
        xPix = [ array(x) for x in self.__recarray__.xPix ]
        yPix = [ array(y) for y in self.__recarray__.yPix ]
        level = [ array(l).astype(int) for l in self.__recarray__.level ]
        
        fname = basename+'.root'
        mode = 'recreate'
        if os.path.exists( fname ):
            mode = 'update'
        
        varnames = ['fname', 'mode', 'N', 'Nmax'] + list(self.__recarray__.dtype.names)
        weave.inline( code, varnames,
                            headers=['"TFile.h"', '"TTree.h"', '"TObject.h"'],
                            libraries=['Core'],
                            include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                            library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                            extra_compile_args=['-O3', '-Wunused-variable'],
                            compiler='gcc',
                            #verbose=1,
                            )
        number_of_entries = len(root_numpy.root2array( fname, treename='hitSumm', branches=['E0'] ))
        if mode == 'recreate':
            print( 'saved catalog {} with {}'.format(fname, number_of_entries) )
        elif mode == 'update':
            print( 'updated catalog {} with {}'.format(fname, number_of_entries) )
        return fname
    
def open_HitSummary( file, branches=None, selection=None, start=None, stop=None, runID_range=None ):
    if runID_range:
        runIDs = open_HitSummary( file, branches=['runID'] ).runID
        start = amin(argwhere(runIDs==runID_range[0]))
        stop = amax(argwhere(runIDs==runID_range[1]))+1        
    return HitSummary( root_numpy.root2array( file, treename='hitSumm', branches=branches, selection=selection, start=start, stop=stop).view(recarray) )


#def list_branches( file, treename = 'hitSumm' ):
    #return root_numpy.list_branches( file, treename = treename )

#def read_catalog( file, treename = 'hitSumm', branches = None, selection = '', indices = None ):
    #if branches is None: branches = list_branches( file, treename )
    #if indices is None:
        #start, stop = None, None
    #else:
        #start, stop = min(indices), max(indices)+1
    #return root_numpy.root2array( file, treename = treename, branches = branches, selection = selection, start = start, stop = stop )

#def list_runIDs( file ):
    #return unique( read_catalog( file, branches = ['runID'] )['runID'] )

#def get_number_of_entries( file ):
    #return len( read_catalog( file, branches = ['runID'] )['runID'] )

#def save_catalog( file, data ):
    #root_numpy.array2root( data, file, treename='hitSumm', mode='recreate' )


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


def build_new_catalog_from_recarray__( input, output ):
    N = len(input)
    print( 'build_new_catalog_from_recarray2', N )
    Nmax = max( input.nSavedPix )
    print( 'Nmax', Nmax )
    with Timer('fill') as t:
        root_numpy.array2root( input, output + '.root', treename='hitSumm', mode='recreate' )


#def build_new_catalog_from_recarray( input, output ):
    #N = len(input)
    #print( input.xPix[1] )
    #print( 'build_new_catalog_from_recarray2', N )
    #Nmax = max( input.nSavedPix )
    #print( 'Nmax', Nmax )
    
    #declarations = 'const Int_t knSavedPix = {};\n'.format(Nmax)
    #assignments = ''
    #assignmentsArray = ''
    #branch = ''
    #for name in input.dtype.names[::-1]:
        #type_ = getattr( input, name ).dtype
        #if name.startswith('n') and name != 'nSavedPix' and name != 'nSat':
            #field = name
            #var = '&c_{}'.format(name)
            #expr = '{}/F'.format(name)
            #declarations += 'Float_t c_{};\n'.format(name)
            #assignments += '\tc_{} = {}[n];\n'.format( name, name )
        #elif type_ == int32:
            #field = name
            #var = '&c_{}'.format(name)
            #expr = '{}/I'.format(name)
            #declarations += 'Int_t c_{};\n'.format(name)
            #assignments += '\tc_{} = {}[n];\n'.format( name, name )
        #elif type_ == float32:
            #field = name
            #var = '&c_{}'.format(name)
            #expr = '{}/F'.format(name)
            #declarations += 'Float_t c_{};\n'.format(name)
            #assignments += '\tc_{} = {}[n];\n'.format( name, name )
        #elif name in ['level','xPix','yPix']:
            #field = '{}[c_nSavedPix]'.format(name)
            #field = name
            #var = 'c_{}'.format(name)
            #expr = '{}[nSavedPix]/I'.format(name)
            #declarations += 'Int_t c_{}[knSavedPix];\n'.format(name)
            #assignmentsArray += '\t\tc_{}[i] = ((Long_t*) PyArray_DATA(PyList_GetItem({},n))) [i];\n'.format( name, name )
        #elif name == 'ePix':
            #field = '{}[c_nSavedPix]'.format(name)
            #field = name
            #var = 'c_{}'.format(name)
            #expr = '{}[nSavedPix]/F'.format(name)
            #declarations += 'Float_t c_{}[knSavedPix];\n'.format(name)
            #assignmentsArray += '\t\tc_{}[i] = ((Double_t*) PyArray_DATA(PyList_GetItem({},n))) [i];\n'.format( name, name )
        #branch += 'tree->Branch( "{}", {}, "{}" );\n'.format( field, var, expr )
    
    #code = declarations + '\n';
    #code += 'TFile f( "{}.root", "recreate" );\n'.format( output )
    #code += 'TTree *tree = new TTree( "hitSumm", "hitSumm" );\n'
    #code += '\n' + branch + '\n';
    #code += 'for(int n=0; n<{}; ++n){{\n'.format( len(input) )
    #code += assignments
    #code += '\tfor(int i=0; i<nSavedPix[n]; ++i){\n'
    #code += assignmentsArray
    #code += '\t}\n'
    #code += '\ttree->Fill();\n'
    #code += '}\n'
    #code += '\ntree->Write();\n'
    
    #for name in input.dtype.names:
        #exec '{} = array(input.{})'.format(name, name)
    #ePix = [ array(e) for e in input.ePix ]
    #xPix = [ array(x) for x in input.xPix ]
    #yPix = [ array(y) for y in input.yPix ]
    #level = [ array(l).astype(int) for l in input.level ]
    #weave.inline( code, input.dtype.names,
                        #headers=['"TFile.h"', '"TTree.h"'],
                        #libraries=['Core'],
                        #include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                        #library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                        ##extra_compile_args=['-O3'],
                        #verbose=1,
                        #)    
        
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

def match( args ):
    simulation = Simulation.simulation_from_file( args.basename )
    hitSummary = open_HitSummary( args.catalog, branches = ['ePix', 'xPix', 'yPix', 'level', 'xBary1', 'yBary1', 'xVar1', 'yVar1'] )
    verbose = 0
    if 'verbose' in args:
        verbose = 1
    hitSummary.match( simulation, verbose, args.basename )

def add_match_options(p):
    p.add_argument('basename', type=str, help = 'basename of the simulation csv file' )
    p.add_argument('catalog', type=str, help = 'basename of the simulation catalog' )
    p.add_argument('-v', '--verbose', action="store_true", default=argparse.SUPPRESS, help = 'verbose' )    
    p.set_defaults( _func=match )

def extract( args ):
    verbose = 0
    if 'verbose' in args:
        verbose = 1
    if os.path.exists( args.basename +'.root' ):
        print( 'catalog {}.root already exists'.format(args.basename) )
        exit(1)
    args.input_files = args.fits_files
    def image_extract( image, args ):
        hitSummary = HitSummary( image.data.extract_hits( args.threshold, args.border, verbose ), transpose=False, verbose=verbose )
        for field, dtype in [('ohdu', int32), ('runID', int32), ('gain', float32), ('rebin',int32), ('dc',float32), ('noise',float32)]:
            if field.upper() in image.header:
                hitSummary.add_fields(field, dtype, [image.header[field.upper()]]*len(hitSummary) )
        print( 'runID {} ohdu {} extracted hits {}'.format( int(image.header['RUNID']), image.header['OHDU'], len(hitSummary) ) )
        hitSummary.add_fields( 'thr', float32, [args.threshold]*len(hitSummary) )
        hitSummary.add_basic_fields()
        hitSummary.add_level_fields(0)
        hitSummary.add_level_fields(1)
        hitSummary.save_catalog( args.basename )
        return
    args.func = image_extract
    image = Image.apply_to_files( args )
    
    
def add_extract_options(p):
    p.add_argument('basename', type=str, help = 'basename of the output catalog' )
    p.add_argument('fits_files', type=str, help = 'fits file (example: )' )
    p.add_argument('-t', '--threshold', type=float, default=60, help = 'energy threshold in ADU' )
    p.add_argument('-b', '--border', type=int, default=0, help = 'pixel border' )
    #p.add_argument('branches', nargs=2, type=str, default= ['E0','xVar0'], help = 'branches used for x- and y-axis' )
    p.add_argument('--ohdu', nargs='+', type=int, default=range(2,16), help = 'ohdus to be extracted' )
    p.add_argument('--exclude', nargs='+', type=int, default=[11,12,14], help = 'ohdus NOT to be extracted' )
    p.add_argument('-v', '--verbose', action="store_true", default=argparse.SUPPRESS, help = 'verbose' )
    p.set_defaults(_func=extract)

def get_selections( file, branches, selections, global_selection=None, runID_range=None, extra_branches = [] ):
    data_selection = {}
    lims = None
    for selection in selections:
        selection_full = selection
        if global_selection:
            selection_full = global_selection + ' and ' + selection
        selection_string = selection_full.replace('and', '&&')
        data_selection[selection] = open_HitSummary( file, branches=list(set(branches+extra_branches)), selection='flag==0 && '+selection_string, runID_range=runID_range )
        if lims is None:
            lims = [ ( min(data_selection[selection][branch]), max(data_selection[selection][branch]) ) for branch in branches ]
        print( selection, data_selection[selection].shape )
    return data_selection, lims


def test_std_correction():
    print( 'mean_norm', mean_norm( 1, 2, -inf, inf ) )
    print( 'mean_norm', mean_norm( 1, 2, 0, inf ) )
    print( 'mean_norm', mean_norm( 1, 2, -inf, 0 ) )
    print( 'mean_norm', mean_norm( 1, 2, 0, 10 ) + mean_norm( 1, 2, -5, 0 ) )

    print( 'var_norm', 1**2/var_norm( 0, 1, -3, 3 ) )
    sigma = 3.11
    rvs = stats.norm.rvs(0, sigma, int(1e3))
    p = stats.norm.pdf(rvs)
    f=0
    std = std(rvs)
    print( 'stats', std, sigma, sum(rvs) )
    rvs = rvs[p>max(p)*f]
    print( 'len', len(rvs) )
    std_ = std(rvs)
    mu = mean(rvs)
    rvsmin, rvsmax = min(rvs), max(rvs) 
    print( 'stats factor', mu, std_, sigma, std*sqrt(std**2/var_norm( mu, std_, rvsmin, rvsmax)) )
    print( 'stats add', mu, std_, sigma, sqrt( (std**2 + var_norm( mu, std_, -inf, rvsmin) + var_norm( mu, std, rvsmax, inf) ))*std**2/var_norm(mu, std, rvsmin, rvsmax) )
    print( 'stats add', mu, std, sigma, sqrt( (std**2 + var_norm( mu, std, -inf, rvsmin) + var_norm( mu, std, rvsmax, inf) )) )
    print( 'stats', mu, std, sigma, std*std**2/var_norm( mu, std, rvsmin, rvsmax) )
    #print( 'stats fwhm', rvsmin, rvsmax, rvsmax-rvsmin, (rvsmax-rvsmin)/( 2*sqrt(2*log(1./f)) ), stats.FWHM1d(rvs, f=.1) )
    print( sum(rvs), sum_norm( mu, sigma, rvsmin, rvsmax ), sum(rvs)/sum_norm( mu, sigma, rvsmin, rvsmax ) )
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
        if not 'runID_range' in args:
            args.runID_range = None
        data_selection, lims = get_selections( file[0], args.branches, args.selections, args.global_selection, runID_range=args.runID_range, extra_branches=[ 'ePix', 'xPix', 'yPix', 'xVar1', 'E0', 'E1', 'n0', 'n1' ] )
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if 'global_selection' in args:
            ax.set_title( args.global_selection )
        if 'selections' in args:
            markers = ['.', '+', 'x', '^']
            for selection, marker in zip(args.selections, markers):
                datum = data_selection[selection]
                datum.add_energy_threshold(30)

                Eerr = 15
                if False:
                    ax.scatter( datum[args.branches[0]], datum[args.branches[1]], marker=marker, alpha=.2, label=selection.replace('&&', 'and') )
                ax.errorbar( datum.E0, datum[args.branches[1]], xerr=Eerr*datum.n0, label='E0', fmt=markers[0] )
                ax.errorbar( datum.E1, datum[args.branches[1]], xerr=Eerr*datum.n1, label='E1', fmt=markers[1] )
                #ax.errorbar( datum['E2'], datum[args.branches[1]], xerr=Eerr*datum['n2'], label='E3', fmt=markers[2] )
                ax.errorbar( datum.E30, var30[0], xerr=Eerr*datum.n30, label='E30', fmt=markers[3] )
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
    p.set_defaults(_func=scatter)

def energy_threshold( datum, threshold ):
    above_2sigma = [ ePix>threshold for ePix in datum.ePix ]
    moment = lambda n: array( [ average( (xPix[mask]**n, yPix[mask]**n), weights = ePix[mask], axis=1 ) for mask, xPix, yPix, ePix in zip(above_2sigma, datum.xPix, datum.yPix, datum.ePix) ] )
    m1 = moment(1)
    m2 = moment(2)
    var = (m2 - m1**2).T
    E, n = array( [ (sum( ePix[mask] ), len(ePix[mask]) ) for mask, ePix in zip(above_2sigma, datum.ePix) ] ).T
    dtype = [ ('{}{}'.format(key,threshold), float) for key in ('E', 'n', 'xVar', 'yVar') ]
    #print( dtype )
    return array( zip(E, n, var[0], var[1]), dtype=dtype ).view(recarray)

def pdf( x, E, n, sigma ):
    return sum( [stats.norm.pdf(x, iE, sqrt(in_)*sigma) for iE, in_ in zip(E,n)], axis=0 )

def histogram( args ):
    with Timer('histogram'):
        file = glob.glob(args.root_file)
        data = open_HitSummary(file[0], args.branch, selection='flag==0' )

        start, stop = None, None
        if 'runID_range' in args:
            start, stop = get_start_stop_from_runID_range(file[0], args.runID_range)
        data_selection, lims = get_selections( file[0], [args.branch], args.selections, args.global_selection, start=start, stop=stop,
                                              extra_branches=['ePix','xPix','yPix', 'E0', 'E1', 'n0', 'n1'] )

        first_data_selection = data_selection.values()[0][args.branch]
        if 'binsize' in args:
            bins = arange( min(first_data_selection), max(first_data_selection), args.binsize )
        elif 'nbins' in args:
            bins = linspace( min(first_data_selection), max(first_data_selection), args.nbins )
        else:
            bins = linspace( min(first_data_selection), max(first_data_selection), int(sqrt(len(first_data_selection))) )
        
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
                
                #pdf30 = lambda x: sum( [stats.norm.pdf(x, E, n*Eerr) for E, n in zip(E30,n30)], axis=0 )
                #ax.step( bins, pdf30(bins), where='mid', label='E30')
                
                #E30 = E30[ logical_or(var30[0]<.3, var30[0]>.8) ]
                #n30 = n30[ logical_or(var30[0]<.3, var30[0]>.8) ]
                #pdf30 = lambda x: sum( [stats.norm.pdf(x, E, n*Eerr) for E, n in zip(E30,n30)], axis=0 )
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
    
    p.set_defaults(_func=histogram)

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
            #data = root_numpy.root2array( args.root_file, treename = args.tree, branches = args.branches ).view(recarray)
            data = open_HitSummary( args.root_file, branches = args.branches )
            for branch in args.branches:
                print( 'branch', getattr( data, branch ), len(getattr( data, branch )) )
    return

def add_status_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('-t', '--tree', type=str, default=argparse.SUPPRESS, help = 'tree to print' )
    p.add_argument('--branches', nargs='+', type=str, default=argparse.SUPPRESS, help = 'branch used for x-axis' )
    p.set_defaults( _func = status )

if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'catalog tools', formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    subparsers = parser.add_subparsers( help = 'major options' )

    add_status_options( subparsers.add_parser('status', help='get status from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_histogram_options( subparsers.add_parser('histogram', help='make histogram from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_scatter_options( subparsers.add_parser('scatter', help='make scatter from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_extract_options( subparsers.add_parser('extract', help='extract events from image', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_match_options( subparsers.add_parser('match', help='match catalog with simulation', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    args = parser.parse_args()
    
    _func = args._func
    del args._func
    print( colored('using parameters:', 'green', attrs=['bold'] ) )
    print_var( vars(args).keys(), vars(args), line_char='\t' )
    
    with Timer('finished'):
        _func(args)
    
