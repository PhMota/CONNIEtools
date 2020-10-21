#!/usr/bin/python

from __future__ import print_function
import os, sys, argparse, re, glob, shutil
from argparse import Namespace

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
from numpy.lib import recfunctions as rfn
import matplotlib
matplotlib.use('gtk3agg')



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
from scipy.stats import binned_statistic
from scipy.special import erf, erfinv


def binned_statistic_fast(x, values, func, nbins, range):
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
        self.shape = recarray.shape
        self.size = recarray.size
        self.names = recarray.dtype.names

    def __getitem__(self, *args):
        if type(args[0]) == str:
            return self.__recarray__[args[0]]
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
        weights = entry.ePix[mask]
        pos = weights > 0
        return average( getattr(entry, '{}Pix'.format(axis))[mask][pos], weights=weights[pos], axis=0 )

    def __Var(self, entry, mask, axis ):
        weights = entry.ePix[mask]
        pos = weights > 0
        return average( getattr(entry, '{}Pix'.format(axis))[mask][pos]**2, weights=weights[pos], axis=0 ) \
                - average( getattr(entry, '{}Pix'.format(axis))[mask][pos], weights=weights[pos], axis=0 )**2

    def __SqrBary(self, entry, mask, axis ):
        weights = entry.ePix[mask]
        pos = weights > 0
        return average( getattr(entry, '{}Pix'.format(axis))[mask][pos]**2, weights=weights[pos], axis=0 )

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

    def add_fit_fields( self ):
        self.t_fit = Timer('done')
        self.t_fit.wait_secs = 30
        self.t_fit.total_loop_number = self.size
        self.t_fit.__enter__()
        arrays = self.__apply_to_entries( self.__fit_vars ).T
        del self.t_fit

        print( 'fit fields false', arrays.shape[-1], sum( arrays[-1] == -1 ), float(sum( arrays[-1] == -1 ))/len(arrays[-1]) )
        self.add_fields( ['Efit', 'xMufit', 'yMufit', 'xSigmafit', 'ySigmafit'],
                        types=(float,float,float,float,float),
                        arrays=arrays
                        )

    def __fit_vars( self, entry ):
        try:
            self.t_fit.check()
        except AttributeError:
            pass

        rawNoise = 10
        #if 'rawNoise' in self.names:
            #rawNoise = self.rawNoise

        #with Timer('likevars'):
        try:
            xmu = entry.xBary1
            ymu = entry.yBary1
            xsigma = sqrt(entry.xVar1)
            ysigma = sqrt(entry.yVar1)
        except:
            xmu = entry.xBary
            ymu = entry.yBary
            xsigma = sqrt(entry.xVar)
            ysigma = sqrt(entry.yVar)
            #print( xmu, ymu, xsigma, ysigma )

        xMufit, yMufit, xSigmafit, ySigmafit, Efit = stats.norm2d_norm.fit(
            entry.xPix,
            entry.yPix,
            entry.ePix,
            xmu,
            ymu,
            xsigma,
            ysigma,
            sum(entry.ePix),
            ftol=1e-8,
            mode='fit',
            sigma_e=rawNoise,
            g=7.25,
            lamb=0
        )
        #print( 'fit', Efit, xMufit, yMufit, xSigmafit, ySigmafit )
        return Efit, xMufit, yMufit, xSigmafit, ySigmafit

    def add_like_fields( self ):
        self.t_like = Timer('done')
        self.t_like.wait_secs = 30
        self.t_like.total_loop_number = self.size
        self.t_like.__enter__()
        arrays = self.__apply_to_entries( self.__like_vars ).T
        del self.t_like

        print( 'like fields false', arrays.shape[-1], sum( arrays[-1] == -1 ), float(sum( arrays[-1] == -1 ))/len(arrays[-1]) )

        self.add_fields( ['ELike', 'xMuLike', 'yMuLike', 'xSigmaLike', 'ySigmaLike'],
                        types=(float,float,float,float,float),
                        arrays=arrays
                        )

    def __like_vars( self, entry ):
        try:
            self.t_like.check()
        except AttributeError:
            pass

        rawNoise = 10
        try:
            xmu = entry.xBary1
            ymu = entry.yBary1
            xsigma = sqrt(entry.xVar1)
            ysigma = sqrt(entry.yVar1)
        except:
            xmu = entry.xBary
            ymu = entry.yBary
            xsigma = sqrt(entry.xVar)
            ysigma = sqrt(entry.yVar)

        xMu, yMu, xSigma, ySigma, E = stats.norm2d_binom_norm.mle(
            entry.xPix,
            entry.yPix,
            entry.ePix,
            xmu,
            ymu,
            xsigma,
            ysigma,
            sum(entry.ePix),
            sigma_e = rawNoise,
            g=7.25,
            lamb=0
        )
        return E, xMu, yMu, xSigma, ySigma

    def add_noise_fields( self ):
        self.t_noise = Timer('done')
        self.t_noise.wait_secs = 30
        self.t_noise.total_loop_number = self.size
        self.t_noise.__enter__()
        arrays = self.__apply_to_entries( self.__noise_vars ).T
        del self.t_noise

        print( 'like fields false', arrays.shape[-1], sum( arrays[-1] == -1 ), float(sum( arrays[-1] == -1 ))/len(arrays[-1]) )

        self.add_fields( ['noise', 'noise2'],
                        types=(float, float),
                        arrays=arrays
                        )

    def __noise_vars( self, entry ):
        try:
            self.t_noise.check()
        except AttributeError:
            pass

        e = entry.ePix
        x = entry.xPix
        y = entry.yPix
        E = e[:,None]*e[None,:]
        d = (x[:,None] - x[None,:])**2 + (y[:,None] - y[None,:])**2
        #v = sum(E[d>0]/d[d>0])/sum(e**2)
        v = E[d>0]/d[d>0]
        #v = sum(E[d>0]/d[d>0])/sum(abs(E[d>0])/d[d>0])
        return sum(v), sum(v)/sum(abs(v))


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
            xPix0 = hit.xPix[ hit.level<=1 ]
            yPix0 = hit.yPix[ hit.level<=1 ]
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

    def efficiency( self, events, args ):
        print( 'fields', self.__recarray__.dtype.names )
        E_bins = arange( 0, max(events.E), args.E_binsize )
        simE, x, dx = stats.make_histogram( events.E, E_bins )
        recE, x, dx = stats.make_histogram( self.E[ self.id_code==1 ], E_bins )
        E0, x, dx = stats.make_histogram( self.E0[ self.id_code==1 ], E_bins )
        E1, x, dx = stats.make_histogram( self.E1[ self.id_code==1 ], E_bins )
        print( simE, recE.astype(float)/simE, E0.astype(float)/simE, E1.astype(float)/simE )


        z_bins = arange( 0, max(events.z), args.z_binsize )
        simz, x, dx = stats.make_histogram( events.z, z_bins )
        recz, x, dx = stats.make_histogram( self.z[ self.id_code==1 ], z_bins )
        print( simz, recz.astype(float)/simz )

        sigma_bins = arange( 0, max(events.z), args.sigma_binsize )
        simsigma, x, dx = stats.make_histogram( events.sigma, sigma_bins )
        recsigma, x, dx = stats.make_histogram( self.sigma[ self.id_code==1 ], sigma_bins )
        yVar0, x, dx = stats.make_histogram( self.yVar0[ self.id_code==1 ], sigma_bins )
        xVar0, x, dx = stats.make_histogram( self.xVar0[ self.id_code==1 ], sigma_bins )
        yVar1, x, dx = stats.make_histogram( self.yVar1[ self.id_code==1 ], sigma_bins )
        xVar1, x, dx = stats.make_histogram( self.xVar1[ self.id_code==1 ], sigma_bins )
        print( simsigma, recsigma.astype(float)/simsigma )
        exit()

    def match_slow( self, events, args ):
        if 'sigma' in self.__recarray__.dtype.names:
            print( 'catalog already matched' )
            return
        N = len(events)

        if 'verbose' in args:
            print( 'events.xy', events.xy )

        xy2event_ind = {}
        for i, event in enumerate(events):
            key = ( int(event.x), int(event.y) )
            xy2event_ind[key] = i
        if 'verbose' in args:
            print( xy2event_ind.keys() )

        N_hits = len(self.__recarray__)
        x = zeros(N_hits)
        y = zeros(N_hits)
        z = zeros(N_hits)
        q = zeros(N_hits)
        E = zeros(N_hits)
        #E_eff = zeros(N_hits)
        id_code = zeros(N_hits)
        sigma = zeros(N_hits)

        matched = zeros(N)

        for i, hit in enumerate(self.__recarray__):
            xPix0 = hit.xPix[ hit.level<=1 ]
            yPix0 = hit.yPix[ hit.level<=1 ]
            if 'verbose' in args:
                print( 'hitSummary.xy', zip(xPix0, yPix0) )
            for xPix, yPix in zip(xPix0, yPix0):
                key = ( xPix, yPix )
                try:
                    event_ind = xy2event_ind[key]
                    event = events[event_ind]
                    matched[event_ind] = 1
                    if 'verbose' in args:
                        print( 'match', event.x, event.y )
                    if id_code[i] == 0:
                        x[i] = event.x
                        y[i] = event.y
                        z[i] = event.z
                        sigma[i] = event.sigma
                        q[i] = event.q
                        E[i] = event.E
                        #E_eff[i] = event.E_eff
                        id_code[i] = event.id_code
                    else:
                        id_code[i] = 2
                        E[i] += event.E
                        #E_eff[i] += event.E_eff
                        q[i] += event.q
                        x[i] = -1
                        y[i] = -1
                        z[i] = -1
                        sigma[i] = -1
                except KeyError:
                    pass

        print( 'matched', count_nonzero(id_code) )
        print( 'multiple matched', count_nonzero(id_code==2) )
        print( 'not matched', count_nonzero(matched==0) )
        if 'verbose' in args:
            print( 'not matched events' )
            print( events.xy[matched==0] )
            print( 'not matched hits' )
            print( zip(self.__recarray__.xBary1[ id_code==0 ], self.__recarray__.yBary1[ id_code==0 ]) )

        fname = args.basename + '.root'
        varnames = ['x','y','z','q','sigma','E','id_code']
        self.add_fields( varnames, (float, float, float, int, float, float, int), [x, y, z, q, sigma, E, id_code] )

        if 'no_catalog' in args:
            return

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

        varnames += ['fname']
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

    def match( self, events, args ):
        if 'sigma' in self.__recarray__.dtype.names:
            print( 'catalog already matched' )
            return

        if 'verbose' in args:
            print( 'events.xy', events.xy )

        with Timer('build occupied bins'):
            xy2hit_ind = { ( xPix, yPix ): i for i, hit in enumerate(self.__recarray__) for xPix, yPix, level in zip(hit.xPix, hit.yPix, hit.level)  }

        if 'verbose' in args:
            print( 'all keys', xy2hit_ind.keys() )

        N_hits = len(self.__recarray__)
        xSim = zeros(N_hits)
        ySim = zeros(N_hits)
        zSim = zeros(N_hits)
        qSim = zeros(N_hits)
        ESim = zeros(N_hits)
        idSim = zeros(N_hits)
        sigmaSim = zeros(N_hits)

        matched = zeros( len(events) )

        with Timer('matching'):
            for event_ind, event in enumerate(events):
                #if 'verbose' in args:
                    #print( 'hitSummary.xy', zip(xPix, yPix) )
                key = ( int(event.x/events.rebin[0]), int(event.y/events.rebin[1]) )
                try:
                    #print( 'key', key )
                    hit_ind = xy2hit_ind[key]
                    matched[event_ind] = 1
                    #print( 'matched', key, hit_ind )
                    #print( 'Bary0', self.__recarray__[hit_ind].xBary0, self.__recarray__[hit_ind].yBary0 )
                    #print( 'Bary1', self.__recarray__[hit_ind].xBary1, self.__recarray__[hit_ind].yBary1 )
                    #print( 'xy', event.x, event.y )
                    #if 'verbose' in args:
                        #print( 'match', event.x, event.y )
                    if idSim[hit_ind] == 0:
                        xSim[hit_ind] = float(event.x/events.rebin[0])
                        ySim[hit_ind] = float(event.y/events.rebin[1])
                        zSim[hit_ind] = float(event.z)
                        sigmaSim[hit_ind] = float(event.sigma)
                        qSim[hit_ind] = int(event.q)
                        ESim[hit_ind] = float(event.E)
                        idSim[hit_ind] = int(event.id_code)
                    else:
                        idSim[hit_ind] = 2
                        ESim[hit_ind] += event.E
                        qSim[hit_ind] += event.q
                        xSim[hit_ind] = -1
                        ySim[hit_ind] = -1
                        zSim[hit_ind] = -1
                        sigmaSim[hit_ind] = -1
                except KeyError:
                    #print( 'NO match', event.x, event.y )
                    pass

        #print( 'x', xSim[xSim!=0] )
        #print( 'y', ySim[ySim!=0] )
        print( 'matched', count_nonzero(idSim) )
        print( 'multiple matched', count_nonzero(idSim==2) )
        print( 'not matched', count_nonzero(matched==0) )
        if 'verbose' in args:
            print( 'not matched events' )
            print( events.xy[matched==0] )
            print( 'not matched hits' )
            #print( zip(self.__recarray__.xBary1[ idSim==0 ], self.__recarray__.yBary1[ idSim==0 ]) )

        fname = args.basename + '.root'
        varnames = ['xSim','ySim','zSim','qSim','sigmaSim','ESim','idSim']

        self.add_fields( varnames, (float, float, float, int, float, float, int), [xSim, ySim, zSim, qSim, sigmaSim, ESim, idSim] )

        print( 'added columns {}'.format(varnames) )


        #if 'no_catalog' in args:
            #return

        #code = '''
        #TFile f( fname.c_str(), "update" );
        #TTree *t = (TTree*) f.Get("hitSumm");
        #Long64_t nentries = t->GetEntries();

        #Float_t c_xSim;
        #Float_t c_ySim;
        #Float_t c_zSim;
        #Int_t c_qSim;
        #Float_t c_sigmaSim;
        #Float_t c_ESim;
        #Int_t c_idSim;

        #TBranch *x_branch = t->Branch("xSim", &c_xSim, "xSim/F");
        #TBranch *y_branch = t->Branch("ySim", &c_ySim, "ySim/F");
        #TBranch *z_branch = t->Branch("zSim", &c_zSim, "zSim/F");
        #TBranch *q_branch = t->Branch("qSim", &c_qSim, "qSim/I");
        #TBranch *sigma_branch = t->Branch("sigmaSim", &c_sigmaSim, "sigmaSim/F");
        #TBranch *E_branch = t->Branch("ESim", &c_ESim, "ESim/F");
        #TBranch *id_code_branch = t->Branch("idSim", &c_idSim, "idSim/I");

        #for( int n=0; n<nentries; ++n ){
            #c_xSim = xSim[n];
            #c_ySim = ySim[n];
            #c_zSim = zSim[n];
            #c_qSim = qSim[n];
            #c_sigmaSim = sigmaSim[n];
            #c_ESim = ESim[n];
            #c_idSim = idSim[n];

            #x_branch->Fill();
            #y_branch->Fill();
            #z_branch->Fill();
            #q_branch->Fill();
            #sigma_branch->Fill();
            #E_branch->Fill();
            #id_code_branch->Fill();
        #}
        #t->Write("", TObject::kOverwrite);
        #'''

        #varnames += ['fname']
        #weave.inline( code, varnames,
                            #headers=['"TFile.h"', '"TTree.h"', '"TObject.h"'],
                            #libraries=['Core'],
                            #include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                            #library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                            #extra_compile_args=['-O3', '-Wunused-variable'],
                            #compiler='gcc',
                            ##verbose=1,
                            #)

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

        names = list(self.__recarray__.dtype.names)
        names.remove('nSavedPix')
        names = ['nSavedPix'] + names
        print( names )

        for name in names:
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
            exec( '{} = array(self.__recarray__.{})'.format(name, name) )

        ePix = [ array(e_) for e_ in self.__recarray__.ePix ]
        xPix = [ array(x_) for x_ in self.__recarray__.xPix ]
        yPix = [ array(y_) for y_ in self.__recarray__.yPix ]
        level = [ array(l_).astype(int) for l_ in self.__recarray__.level ]

        fname = basename
        mode = 'recreate'
        if os.path.exists( fname ):
            mode = 'update'

        #print( 'code', code )

        varnames = ['fname', 'mode', 'N', 'Nmax'] + list(self.__recarray__.dtype.names)
        weave.inline( code, varnames,
                            headers=['"TFile.h"', '"TTree.h"', '"TObject.h"'],
                            libraries=['Core'],
                            include_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/include/root/'],
                            library_dirs=['/opt/versatushpc/softwares/root/5.34-gnu-5.3/lib/'],
                            extra_compile_args=['-O3', '-Wunused-variable'],
                            compiler='gcc',
                            verbose=0,
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
    #if not branches is None:
        #branches = list( set(branches) & set(root_numpy.list_branches( file, treename='hitSumm' )) )
    print( file, root_numpy.list_branches( file, treename='hitSumm' ) )
    return HitSummary( root_numpy.root2array( file, treename='hitSumm', branches=branches, selection=selection, start=start, stop=stop).view(recarray) )

def update_catalog( fname, output, hitSummary ):
    existing_branches = root_numpy.list_branches( fname, treename='hitSumm' )
    a = list( set(hitSummary.names) - set(existing_branches) )
    print( 'adding branches', a, 'to catalog', output )
    #print( hitSummary.__recarray__[a] )
    shutil.copyfile(fname, output)
    root_numpy.array2root( hitSummary.__recarray__[a], output, treename='hitSumm', mode='update' )

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

def addbranches( **args ):
    with Timer('add branches'):
        args = Namespace(**args)

        fname = glob.glob(args.basename)[0]
        fout = args.output
        if not os.path.exists(fname):
            raise FileNotFoundError

        if os.path.exists(fout):
            print( 'output already exists', fout )
            raise FileExistsError

        hitSummary = open_HitSummary(fname)
        print( hitSummary.names )
        branches_set = set(args.branches) - set(hitSummary.names)

        fit_set = set(['Efit', 'xSigmafit', 'xMufit', 'ySigmafit', 'yMufit', 'fit'])
        like_set = set( ['ELike', 'xSigmaLike', 'xMuLike', 'ySigmaLike', 'yMuLike', 'like'] )
        noise_set = set( ['noise', 'noise2'] )

        if not branches_set.isdisjoint(fit_set):
            branches_set -= fit_set
            print( 'adding fit fields' )
            hitSummary.add_fit_fields()

        if not branches_set.isdisjoint(like_set):
            branches_set -= like_set
            print( 'adding like fields' )
            hitSummary.add_like_fields()

        if not branches_set.isdisjoint(noise_set):
            branches_set -= noise_set
            print( 'adding noise fields' )
            hitSummary.add_noise_fields()

        if len(branches_set) > 0:
            print( 'branches not recognized', list(branches_set) )

        update_catalog( fname, fout, hitSummary )
    return

def add_addbranches_options(p):
    p.add_argument('basename', type=str, help = 'basename of the simulation csv file' )
    p.add_argument('branches', nargs='+', type=str, default=None, help = 'branches to be added' )
    p.add_argument('-v', '--verbose', action="store_true", default=argparse.SUPPRESS, help = 'verbose' )
    p.add_argument('-o', '--output', type=str, default='added', help = 'output name' )
    p.set_defaults( _func=addbranches )

def simulate( args ):
    if os.path.exists(args.basename):
        print( 'folder already exists', args.basename )
        exit()
    os.mkdir( args.basename )
    print( 'created folder', args.basename )
    basename = str(args.basename)
    args.basename = '{0}/{0}'.format(args.basename)
    args.no_fits = True
    args.png = True
    args.csv = True
    images = {}
    simulations = {}
    with Timer('simulate all') as t:
        for image_mode in args.image_modes:
            images[image_mode] = []
            simulations[image_mode] = []
            for i in range(args.number_of_images):
                args.runID = i
                args.image_mode = image_mode
                args.basename = '{0}/{0}bin{1}_{2}'.format( basename, image_mode, i )
                simulation = Simulation.simulate_events( args )
                HDUlist = simulation.generate_image( args )
                print( 'saved {}.csv'.format( args.basename ) )
                imageHDU = HDUlist[-1]
                part = Image.Part(imageHDU)
                images[image_mode].append( Image.Image([part]) )
                simulations[image_mode].append( simulation )
                t.check( 10, args.number_of_images*len(args.image_modes) )

    thresholds = list(args.threshold)
    with Timer('extract&match all') as t2:
        for image_mode in args.image_modes:
            for threshold in thresholds:
                args.threshold = threshold
                args.basename = '{0}/{0}_bin{1}t{2}'.format( basename, image_mode, threshold )
                args.no_catalog = True
                for image, simulation in zip(images[image_mode],simulations[image_mode]):
                    with Timer('extract'):
                        hitSummary = image_extract( image, args )
                    with Timer('match'):
                        hitSummary.match( simulation, args )
                    with Timer('save catalog'):
                        hitSummary.save_catalog( args.basename + '.root' )
                    t2.check( 60, len(thresholds)*len(images)*len(args.image_modes) )


def add_simulate_options(p):
    p.add_argument('basename', type=str, help = 'basename of the simulation csv file' )
    p.add_argument('-v', '--verbose', action="store_true", default=argparse.SUPPRESS, help = 'verbose' )
    p.add_argument('-t', '--threshold', nargs='+', type=float, default=[60,50,40,30], help = 'energy threshold in ADU' )
    p.add_argument('-b', '--border', type=int, default=1, help = 'pixel border' )
    p.add_argument('--image-modes', nargs='+', type=int, default=[1,5], help = 'image modes' )
    #p.add_argument('--E-binsize', type=float, default=1, help = 'E binsize' )
    #p.add_argument('--z-binsize', type=float, default=10, help = 'z binsize' )
    #p.add_argument('--sigma-binsize', type=float, default=.1, help = 'sigma binsize' )
    #p.add_argument('--efficiency', action='store_true', default=argparse.SUPPRESS, help = 'sigma binsize' )
    p.add_argument('--number-of-images', type=int, default=1, help = 'number of images to each threshold' )
    Simulation.add_params_options(p)
    Simulation.add_geometry_options(p)
    Simulation.add_depth_options(p)
    Simulation.add_charges_options(p)
    Simulation.add_modulation_options(p)

    p.set_defaults( _func=simulate )

def match( args ):
    if os.path.exists(args.output):
        print( 'matched already exists' )
        exit()
    simulation = Simulation.simulation_from_file( args.basename, invert=args.invert, shift=args.x_shift )
    hitSummary = open_HitSummary( args.catalog, branches = ['ePix', 'xPix', 'yPix', 'level'] )
    verbose = 0
    if 'verbose' in args:
        verbose = 1
    hitSummary.match( simulation, args )
    update_catalog( args.catalog, args.output, hitSummary )


def add_match_options(p):
    p.add_argument('basename', type=str, help = 'basename of the simulation csv file' )
    p.add_argument('catalog', type=str, help = 'basename of the simulation catalog' )
    p.add_argument('-v', '--verbose', action="store_true", default=argparse.SUPPRESS, help = 'verbose' )
    p.add_argument('-o', '--output', type=str, default='matched.root', help = 'output' )
    p.add_argument('--invert', type=bool, default=False, help = 'invert x and y' )
    p.add_argument('--x-shift', type=int, default=0, help = 'shift x axis' )
    p.set_defaults( _func=match )

def image_extract( image, args ):
    verbose = 0
    if 'verbose' in args:
        verbose = 1

    if 'overscan_subtraction' in args:
        os_range = args.overscan_subtraction
        print( 'os_range', os_range, os_range[0] )
        correction = median( image.all[ None:None, os_range[0]:os_range[1] ], axis=1 )[:,None]
        image.all -= correction

    if 'vertical_overscan_subtraction' in args:
        os_range = args.vertical_overscan_subtraction
        print( 'os_range', os_range, os_range[0] )
        correction = median( image.all[ os_range[0]:os_range[1], None:None ], axis=0 )[None,:]
        image.all -= correction

    hitSummary = HitSummary( image.all.extract_hits( args.threshold, args.border, verbose ), transpose=False, verbose=verbose )
    for field, dtype in [('ohdu', int32), ('runID', int32), ('gain', float32), ('rebin',int32), ('dc',float32), ('noise',float32)]:
        if field.upper() in image.header:
            hitSummary.add_fields(field, dtype, [image.header[field.upper()]]*len(hitSummary) )
    try:
        print( 'runID {} ohdu {} extracted hits {}'.format( int(image.header['RUNID']), image.header['OHDU'], len(hitSummary) ) )
    except:
        print( 'extracted hits {}'.format( len(hitSummary) ) )
    hitSummary.add_fields( 'thr', float32, [args.threshold]*len(hitSummary) )
    hitSummary.add_basic_fields()
    for lvl in range( args.border+1 ):
        hitSummary.add_level_fields(lvl)
    if not 'no_catalog' in args:
        hitSummary.save_catalog( args.basename + '.root' )
    return hitSummary

def extract( args ):
    if os.path.exists( args.basename +'.root' ):
        print( 'catalog {}.root already exists'.format(args.basename) )
        exit(1)
    args.input_files = args.fits_files
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
    p.add_argument('--hdu', nargs='+', type=int, default=range(1,15), help = 'hdus (indices) to be extracted' )
    p.add_argument('-v', '--verbose', action="store_true", default=argparse.SUPPRESS, help = 'verbose' )
    p.add_argument( '--overscan-subtraction', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'horizontal overscan subtraction' )
    p.add_argument( '--vertical-overscan-subtraction', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'vertical overscan subtraction' )

    p.set_defaults(_func=extract)

def get_selections( file, branches, selections, global_selection=None, runID_range=None, extra_branches = [] ):
    data_selection = {}
    for selection in selections:
        selection_full = selection
        if global_selection:
            selection_full = global_selection + ' and ' + selection
        selection_string = selection_full.replace('and', '&&')
        list_of_branches = list(set(branches+extra_branches))
        # if type(file) is list:
        #     for f in file:
        #     read_data = open_HitSummary( file, branches=list_of_branches, selection='flag==0 && '+selection_string, runID_range=runID_range )
        #     data_selection[selection] = rfn.append_fields(
        #                                     data_selection[selection],
        #                                     list_of_branches,
        #                                     read_data,
        #                                     asrecarray = True
        #                                     )
        # else:
        data_selection[selection] = open_HitSummary( file, branches=list_of_branches, selection='flag==0 && '+selection_string, runID_range=runID_range )
        print( 'read fileds', data_selection[selection].names )
    return data_selection


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

def scatter( **args ):
    args = Namespace(**args)
    import matplotlib.pylab as plt
    with Timer('scatter'):

        file = glob.glob(args.root_file)

        args.xbranches = []
        args.ybranches = []
        args.cbranches = []
        args.selections = []
        has_color = False

        for branch_selection in args.branch_selections:
            args.xbranches.append( branch_selection[0] )
            args.ybranches.append( branch_selection[1] )
            if len(branch_selection) == 4:
                args.cbranches.append( branch_selection[2] )
                has_color = True
            args.selections.append( branch_selection[-1] )
        if not has_color:
            args.cbranches = args.xbranches

        if not 'runID_range' in args:
            args.runID_range = None
        data_selection = get_selections( file[0],
                                        args.xbranches + args.ybranches + args.cbranches,
                                        args.selections,
                                        args.global_selection,
                                        runID_range=args.runID_range,
                                        #extra_branches=[ 'ePix', 'xPix', 'yPix', 'level' ]
                                        )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if 'global_selection' in args:
            ax.set_title( args.global_selection )
        scatter_obj = []
        if 'selections' in args:
            markers = ['.', '+', 'x', '^']
            for xbranch, ybranch, cbranches, selection, marker in zip(args.xbranches, args.ybranches, args.cbranches, args.selections, markers):
                print( 'plot', xbranch, ybranch, cbranches )
                datum = data_selection[selection]

                if has_color:
                    colors = datum[cbranches]
                    cmap = matplotlib.cm.plasma
                    alpha = 1
                else:
                    colors = None
                    cmap = None
                    alpha = (len(datum))**(-.1) if len(datum) > 0 else 1

                if xbranch in datum.names and ybranch in datum.names or True:

                    x = datum[xbranch]
                    y = datum[ybranch]
                    #x = x[ logical_and(args.x_range[0] < x, x < args.x_range[1], axis=0) ]
                    #y = y[ logical_and(args.y_range[0] < x, x < args.y_range[1], axis=0) ]

                    scatter_obj.append( ax.scatter( x, y, label = '{} vs. {} ({})\n{}'.format( ybranch, xbranch, x.size, selection ), marker=marker, c=colors, alpha=alpha, cmap=cmap ) )

                    if 'errorbar' in args:
                        bins = arange( args.x_range[0], args.x_range[1], 20 )
                        bin_means, bin_edges, binnumber = binned_statistic( x, y, statistic='mean', bins=bins )
                        xbins = .5*(bin_edges[1:] + bin_edges[:-1])
                        dx = bin_edges[1] - bin_edges[0]
                        yerr = [0]*len(bin_means)
                        bin_std, bin_edges, binnumber = binned_statistic( x, y, statistic='std', bins=bin_edges )
                        yerr = bin_std
                        ax.errorbar( xbins, bin_means, xerr=dx/2, yerr=yerr, fmt='.' )


        ax.legend()
        ax.grid()
        ax.set_xlim( *args.x_range )
        ax.set_ylim( *args.y_range )
        ax.set_xlabel( args.xbranches[0] )
        ax.set_ylabel( args.ybranches[0] )
        if has_color:
            fig.colorbar( scatter_obj[0] )
    if 'ylog' in args:
        ax.set_yscale('log')

    if 'output' in args:
        if args.output == '':
            args.output = args.root_file
    else:
        args.output = args.root_file

    if 'pdf' in args:
        extra = ''
        fname = args.output+'.scatter.{}.vs.{}{}.pdf'.format(args.ybranches[0].replace('/','_'), args.xbranches[0].replace('/','_'), extra)
        fig.savefig( fname )
        print( 'saved', fname )
    elif 'png' in args:
        #extra = ''
        #if 'fit' in args:
            #extra += '_fit{}'.format(''.join(args.fit))
        #fname = args.output+'.scatter.{}.vs.{}{}.png'.format(args.ybranches[0].replace('/','_'), args.xbranches[0].replace('/','_'), extra)
        fname = args.output+'.png'
        fig.savefig( fname )
        print( 'saved', fname )
    else:
        plt.show()
    return

def add_scatter_options(p):
    p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('-s', '--branch-selections', action='append', nargs='+', type=str, default=argparse.SUPPRESS, help = 'branches used for x- and y-axis' )
    p.add_argument('--global-selection', type=str, default='1', help = 'global selection' )
    #p.add_argument('--selections', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--define', nargs='+', type=str, default=argparse.SUPPRESS, help = 'definitions' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    p.add_argument('--x-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of x' )
    p.add_argument('--y-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of y' )
    p.add_argument('--c-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of color' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'output to file' )
    p.add_argument('--pdf', action='store_true', default=argparse.SUPPRESS, help = 'output to pdf' )
    p.add_argument('--png', action='store_true', default=argparse.SUPPRESS, help = 'output to png' )
    p.add_argument('--fit', nargs='+', type=str, default=argparse.SUPPRESS, help = 'include fit column' )
    p.add_argument('--noise', nargs='+', type=str, default=argparse.SUPPRESS, help = 'include fit column' )
    p.add_argument('--errorbar', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )
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

def histogram( **args ):
    args = Namespace(**args)

    import matplotlib.pylab as plt
    with Timer('histogram'):
        print( 'len', len(args.root_file) )
        if len(args.root_file) == 1:
            args.root_file = args.root_file[0]

            file = glob.glob(args.root_file)
            #data = open_HitSummary(file[0], args.branch, selection='flag==0' )

            args.branches = []
            args.selections = []
            print( args.branch_selections )
            for branch_selection in args.branch_selections:
                args.branches.append( branch_selection[0] )
                args.selections.append( branch_selection[1] )

            if not 'runID_range' in args:
                args.runID_range = None
            data_selection = get_selections( file[0], args.branches, args.selections, args.global_selection, runID_range=args.runID_range )
                                                #extra_branches=['ePix','xPix','yPix', 'E0', 'E1', 'n0', 'n1'] )
        else:
            nargs = len(args.root_file)
            files = args.root_file
            print( args.root_file )
            args.branches = []
            args.selections = []
            data_selection = {}
            for i in range(nargs/3):
                branch = args.root_file[3*i+0]
                args.branches.append(branch)

                file = args.root_file[3*i+1]

                selection = args.root_file[3*i+2]
                args.selections.append('{}:{}:{}'.format(branch,file,selection))
                print( 'file', file )
                print( 'br', branch )
                print( 'sel', selection )
                for f in glob.glob(file):
                    print( 'file', f )
                    data_entry = get_selections( f, [branch], [selection], args.global_selection )

                    print( 'type', data_entry.values()[0], data_entry.values()[0].size )

                    key = '{}:{}:{}'.format(branch,file, data_entry.keys()[0])
                    try:
                        data_selection[key].append( data_entry.values()[0] )
                    except KeyError:
                        data_selection[key] = [ data_entry.values()[0] ]
                    # data_selection.update( { '{}:{}:{}'.format(branch,file, data_entry.keys()[0]): data_entry.values()[0] } )
                print( type(data_selection) )
            print( 'selections', data_selection.keys() )
            print( 'branches', args.branches )
            #exit(0)


        if 'x_range' in args:
            bins = arange( args.x_range[0], args.x_range[1], args.binsize )
        else:
            first_data_selection = data_selection.values()[0][args.branches[0]]
            if 'binsize' in args:
                bins = arange( min(first_data_selection), max(first_data_selection), args.binsize )
            elif 'nbins' in args:
                bins = linspace( min(first_data_selection), max(first_data_selection), args.nbins )
            else:
                bins = linspace( min(first_data_selection), max(first_data_selection), int(sqrt(len(first_data_selection))) )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(', '.join(args.branches))
        #ax.hist(data, bins=bins, histtype='step', label='all')
        if 'selections' in args:
            for i, (branch, selection) in enumerate(zip(args.branches, args.selections)):
                print( 'selection', selection, data_selection[selection].names )
                print( 'selection', data_selection[selection][branch].shape, len(bins) )
                datum = data_selection[selection]

                factor = 1
                if 'factor' in args:
                    factor = args.factor
                hist, x, dx = stats.make_histogram( datum[branch], bins )
                ax.errorbar( x+i*dx/2./len(args.branches), hist*factor, xerr=dx/2., yerr=sqrt(hist)*factor, label='{} ({})'.format(selection,datum[branch].size), fmt='.' )
        ax.legend()
        ax.grid()
        ax.set_xlabel(args.branches[0])
        ax.set_ylabel( r'$\frac{{dN}}{{d {}}}$'.format(args.branches[0]) )
        if 'log' in args:
            ax.set_yscale('log')

    if 'output' in args:
        if args.output == '':
            args.output = args.root_file
    else:
        args.output = args.root_file

    if 'pdf' in args:
        fig.savefig(args.output+'.pdf')
        print( 'saved', args.output+'.pdf' )
    elif 'png' in args:
        fig.savefig(args.output+'.png')
        print( 'saved', args.output+'.png' )
    else:
        plt.show()
    return

def add_histogram_options(p):
    p.add_argument('root_file', nargs='+', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('--branch-selections', action='append', nargs=2, type=str, default=argparse.SUPPRESS, help = 'selections' )
    p.add_argument('--global-selection', type=str, default='1', help = 'global selection' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--energy-threshold', nargs='+', type=float, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--define', type=str, default=argparse.SUPPRESS, help = 'definitions (ex.: a=E0; b=E1)' )
    p.add_argument('--x-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of the x-axis' )
    p.add_argument('--binsize', type=eval, default=argparse.SUPPRESS, help = 'binsize' )
    p.add_argument('--factor', type=eval, default=argparse.SUPPRESS, help = 'factor' )
    p.add_argument('--nbins', type=int, default=argparse.SUPPRESS, help = 'number of bins' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--log', action='store_true', default=argparse.SUPPRESS, help = 'log' )
    p.add_argument('--pdf', action='store_true', default=argparse.SUPPRESS, help = 'output to pdf' )
    p.add_argument('--png', action='store_true', default=argparse.SUPPRESS, help = 'output to png' )

    p.set_defaults(_func=histogram)

def add_ratio_options(p):
    p.add_argument('root_file', nargs='+', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('--branch-selections', action='append', nargs=3, type=str, default=argparse.SUPPRESS, help = 'selections as branch numerator and denominator' )
    p.add_argument('--global-selection', type=str, default='1', help = 'global selection' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--energy-threshold', nargs='+', type=float, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--define', type=str, default=argparse.SUPPRESS, help = 'definitions (ex.: a=E0; b=E1)' )
    p.add_argument('--x-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of the x-axis' )
    p.add_argument('--binsize', type=eval, default=argparse.SUPPRESS, help = 'binsize' )
    p.add_argument('--factor', type=eval, default=argparse.SUPPRESS, help = 'factor' )
    p.add_argument('--nbins', type=int, default=argparse.SUPPRESS, help = 'number of bins' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--pdf', action='store_true', default=argparse.SUPPRESS, help = 'output to pdf' )
    p.add_argument('--png', action='store_true', default=argparse.SUPPRESS, help = 'output to png' )

    p.set_defaults(_func=ratio)

def ratio( **args ):
    args = Namespace(**args)

    import matplotlib.pylab as plt
    with Timer('ratio'):
        print( 'len', len(args.root_file) )
        if len(args.root_file) == 1:
            args.root_file = args.root_file[0]

            file = glob.glob(args.root_file)
            #data = open_HitSummary(file[0], args.branch, selection='flag==0' )

            args.branches = []
            args.numselections = []
            args.denselections = []
            print( args.branch_selections )
            for branch_selection in args.branch_selections:
                args.branches.append( branch_selection[0] )
                args.numselections.append( branch_selection[1] )
                args.denselections.append( branch_selection[2] )

            if not 'runID_range' in args:
                args.runID_range = None
            data_numselection = get_selections( file[0], args.branches, args.numselections, args.global_selection, runID_range=args.runID_range )
            data_denselection = get_selections( file[0], args.branches, args.denselections, args.global_selection, runID_range=args.runID_range )
                                                #extra_branches=['ePix','xPix','yPix', 'E0', 'E1', 'n0', 'n1'] )

        else:
            nargs = len(args.root_file)
            files = args.root_file
            print( args.root_file )
            args.branches = []
            args.numselections = []
            args.denselections = []
            data_numselection = {}
            data_denselection = {}
            for i in range(nargs/4):
                branch = args.root_file[5*i+0]
                args.branches.append(branch)

                file0 = args.root_file[5*i+1]
                selection0 = args.root_file[5*i+2]
                args.numselections.append('{}:{}'.format(file0, selection0))

                file1 = args.root_file[5*i+3]
                selection1 = args.root_file[5*i+4]
                args.denselections.append('{}:{}'.format(file1, selection1))

                #print( 'file', file )
                #print( 'br', branch )
                #print( 'sel', selection )
                data_entry0 = get_selections( file0, [branch], [selection0], args.global_selection )
                data_numselection.update( { '{}:{}'.format(file0, data_entry0.keys()[0]): data_entry0.values()[0] } )

                data_entry1 = get_selections( file1, [branch], [selection1], args.global_selection )
                data_denselection.update( { '{}:{}'.format(file1, data_entry1.keys()[0]): data_entry1.values()[0] } )

            #print( data_selection.keys() )
            #exit(0)



        if 'x_range' in args:
            bins = arange( args.x_range[0], args.x_range[1], args.binsize )
        else:
            first_data_selection = data_numselection.values()[0][args.branches[0]]
            if 'binsize' in args:
                bins = arange( min(first_data_selection), max(first_data_selection), args.binsize )
            elif 'nbins' in args:
                bins = linspace( min(first_data_selection), max(first_data_selection), args.nbins )
            else:
                bins = linspace( min(first_data_selection), max(first_data_selection), int(sqrt(len(first_data_selection))) )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(', '.join(args.branches))
        #ax.hist(data, bins=bins, histtype='step', label='all')
        for branch, numselection, denselection in zip(args.branches, args.numselections, args.denselections):
            #print( 'selection', data_selection[selection][branch].shape, len(bins) )
            numerator = data_numselection[numselection]
            denominator = data_denselection[denselection]

            numhist, x, dx = stats.make_histogram( numerator[branch], bins )
            denhist, x, dx = stats.make_histogram( denominator[branch], bins )
            hist = numhist.astype(float)/denhist
            dy = sqrt( (sqrt(numhist)/denhist)**2 + (numhist*sqrt(denhist)/denhist**2)**2 )
            ax.errorbar( x, hist, xerr=dx/2, yerr=dy, label='{}:{}\n{}'.format(branch,numselection, denselection), fmt=' ' )
        ax.legend()
        ax.set_xlabel(args.branches[0])
        #ax.set_ylabel( r'$\frac{{dN}}{{d {}}}$'.format(args.branches[0]) )

    if 'output' in args:
        if args.output == '':
            args.output = args.root_file
    else:
        args.output = args.root_file

    if 'pdf' in args:
        fig.savefig(args.output+'.pdf')
        print( 'saved', args.output+'.pdf' )
    elif 'png' in args:
        fig.savefig(args.output+'.png')
        print( 'saved', args.output+'.png' )
    else:
        plt.show()
    return


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
    add_simulate_options( subparsers.add_parser('simulation', help='simulate, extract and match', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_addbranches_options( subparsers.add_parser('addbranches', help='add branches to catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_ratio_options( subparsers.add_parser('ratio', help='make ratio histograms from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )

    args = parser.parse_args()

    _func = args._func
    del args._func
    print( colored('using parameters:', 'green', attrs=['bold'] ) )
    print_var( vars(args).keys(), vars(args), line_char='\t' )

    with Timer('finished'):
        try:
            _func(args)
        except TypeError:
            _func(**vars(args))
