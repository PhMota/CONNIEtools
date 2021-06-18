#!/usr/bin/python

from __future__ import print_function
import os, sys, argparse, re, glob, shutil
from argparse import Namespace

from Timer import Timer
# import warnings
# warnings.filterwarnings("ignore")

# try:
#     import root_numpy
# except ImportError:
#     print('missing module, please run')
#     print('module load softwares/root/5.34-gnu-5.3')
#     exit(0)

# import weave
# from numpy.lib import recfunctions as rfn

import matplotlib
    
# matplotlib.use('gtk3agg')
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Helvetica"]})
# matplotlib.rc('text', usetex=True)


from numpy import *
# from numpy.lib.recfunctions import append_fields, stack_arrays
import Statistics as stats

# with Timer('import'):
# from ROOT import TFile, TTree, AddressOf

import Image
import Simulation
import Statistics as stats
from termcolor import colored
from PrintVar import print_var
from scipy.sparse import csr_matrix
from scipy.interpolate import interp1d
# from scipy.stats import binned_statistic
from scipy.special import erf, erfinv
from collections import *
from fstring import F
from progress import progressbar

def binned_statistic_fast(x, values, func, nbins, range):
    '''The usage is nearly the same as scipy.stats.binned_statistic'''

    r0, r1 = range

    digitized = (float(nbins)/(r1 - r0)*(x - r0)).astype(int)
    values = values[digitized < nbins].astype(float)
    digitized = digitized[digitized < nbins]

    N = len(values)
    S = csr_matrix((values, [digitized, arange(N)]), shape=(nbins, N))

    return [func(group) for group in split(S.data, S.indptr[1:-1])]

# def binned_statistic_fast(x, values, func, nbins, range):
#     '''The usage is nearly the same as scipy.stats.binned_statistic'''
#
#     N = len(values)
#     r0, r1 = range
#
#     digitized = (float(nbins)/(r1 - r0)*(x - r0)).astype(int)
#     S = csr_matrix((values, [digitized, arange(N)]), shape=(nbins, N))
#
#     return [func(group) for group in split(S.data, S.indptr[1:-1])]

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
    # print( ', '.join( root_numpy.list_branches( file, treename='hitSumm' )) )
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
        # print( 'read fileds', data_selection[selection].names )
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


def roundRel(v, rel):
    decimals = -log10(v) - rel
    return round( v, int(ceil(decimals)) )
def wmean(v, w=1):
    return average(v, weights=w)
def wstd(v, w=1):
    return sqrt( wmean(v**2, w) - wmean(v, w)**2 )
def unique_unsorted(v):
    _, idx = unique(v, return_index=True)
    x = bins[1:] + bins[:-1]
    dx = bins[1:] - bins[:-1]
    y = v[ sort(idx) ]
    return x/2, y, dx/2, sqrt(y)
class HistogramExpr:
    def __init__(self, x, y, dx, dy):
        self.x, self.y, self.dx, self.dy = x, y, dx, dy
    def astuple(self):
        return self.x, self.y, self.dx, self.dy
    def __tuple__(self):
        return self.astuple()
    def __add__(self, a):
        if isinstance( a, HistogramExpr ):
            self.y += a.y
            self.dy = sqrt( self.dy**2 + a.dy**2 )
            return self
        self.y += a
        return self
    def __mul__(self, a):
        if isinstance( a, HistogramExpr ):
            self.dy = self.y*a.y * sqrt( (self.dy/self.y)**2 + (a.dy/a.y)**2 )
            self.y *= a.y
            return self
        self.y *= a
        self.dy *= a
        return self
    def __rmul__(self, a):
        return self*a
    def __sub__(self, a):
        return self + (-1)*a
    def __div__(self, a):
        if isinstance( a, HistogramExpr ):
            self.dy = self.y/a.y * sqrt( (self.dy/self.y)**2 + (a.dy/a.y)**2 )
            self.y /= a.y
            return self
        self.y /= a
        self.dy /= a
        return self

class Table:
    def __init__(self):
        self.table = []
        self.header = []
    def append( self, header, column ):
        self.table.append(column)
        self.header.append(header)
    def lines(self):
        return [self.header] + zip(*self.table)
    def csv(self):
        return '\n'.join(map(lambda line: ', '.join(map(str,line)),  self.lines() ))
    def save(self, name):
        return savetxt( name, zip(*self.table), header = ', '.join(self.header), delimiter = ', ' )

# table = Table()
# table.append('test', [1,2,3,4])
# table.append('test2', [2,3,4,5])
# table.append('test3', [3,4,5,10])
# print( table.lines() )
# print( table.csv() )
# table.save('test.csv')

def hist(v, bins_, norm=False, edge=False):
    y, x, dx = stats.make_histogram( v, bins_ )
    y = y.astype(float)
    if norm:
        yerr = sqrt(y)/sum(y)
        y /= sum(y)
    else:
        yerr = sqrt(y)
    if edge:
        x = bins_[:-1]
    return HistogramExpr(x, y, dx, yerr)

def hist_pos(v, bins_, norm=False, edge=False):
    x, y, dx, dy = hist(v, bins_, norm, edge=edge).astuple()
    mask = y>0
    return HistogramExpr(x[mask], y[mask], dx, dy[mask])

def waverage( hists ):
    normalization = sum( [ 1./h.dy**2 for h in hists], axis=0 )
    summation = sum( [ h/h.dy**2 for h in hists], axis=0 )
    return summation/normalization

def unique_tuples( x, y, dx, dy ):
    if not hasattr( dx, '__iter__' ):
        dx = [dx]*len(x)
    if not hasattr( dy, '__iter__' ):
        dy = [dy]*len(y)
    return unique( (x, y, dx, dy), axis=-1)

def nu_spectrum(E, file='vSpectrum.csv'):
    '''
        vSpectrum = genfromtxt(file)
        x0, y0 = vSpectrum.T
        y = interp1d(x0,y0)(E)
        return y
    '''
    return interp1d( *genfromtxt(file).T )(E)

# print( nu_spectrum([.4,.401]) )
# exit()
# def nu_spectrum(*a, isotope=None):
#     x = [7.813e-3, 1.563e-2, 3.12e-2, 6.25e-2, .125, .25, .50, .75, 1.0, 1.5, 2.0 ]
#     if isotope == 'U235':
#         a = [1.0461, .87, -.16, -.091]
#         y = [.024, .092, .35, .61, 1.98, 2.16, 2.66, 2.66, 2.41, 1.69, 1.26 ]
#     elif isotope == 'Pu239':
#         a = [1.0527, .896, -.239, -.0981]
#         y = [.14, .56, 2.13, .64, 1.99, 2.08, 2.63, 2.58, 2.32, 1.48, 1.08 ]
#     elif isotope == 'U238':
#         a = [1.0719, .976, -.162, -.079]
#         y = [.089, .35, 1.32, .65, 2.02, 2.18, 2.91, 2.96, 2.75, 1.97, 1.50 ]
#     elif isotope == 'Pu241':
#         a = [1.0818, .793, -0.08, -.1085]
#         y = [.2, .79, 3.00, .59, 1.85, 2.14, 2.82, 2.9, 2.63, 1.75, 1.32 ]
#     def _(E):
#         ret = zeros_like(E)
#         ret[E<2] = scipy.interpolate.interp1d(x, y)(E[E<2])
#         ret[E>=2] = ( lambda _E: a[0]*exp( a[1] + a[2]*_E + a[3]*_E**2 ) )( E[E>2] )
#         return ret
#     return _
#
# def quenching_factor(E):
#     p = [56, 1096, 382, 168, 155 ]
#     return (p[3]*E + p[4]*E**2 + E**3)/(p[0] + p[1]*E + p[2]*E**2)
#
# def combined_nu_spectrum(E):
#     ret = zeros_like(E)
#     for
#
# def differential_rate_recoil(CEvNS_crossSection, v_flux, numberNuclei):
#     Emin = E_recoil + sqrt(E_recoil**2 + 2*M*E_recoil)/2
#     return numberNuclei * quad( v_flux*CEvNS_crossSection, a=Emin, b=inf )




from utils.doc2argparse import doc2argparse
from Catalog.status import status
from Catalog.histogram import histogram
from Catalog.scatter import scatter
from Catalog.ratio import ratio

if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'catalog tools', formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    subparsers = parser.add_subparsers( help = 'major options' )

    doc2argparse( subparsers, status )
    doc2argparse( subparsers, histogram )
    doc2argparse( subparsers, scatter )
    doc2argparse( subparsers, ratio )
#     add_histogram_options( subparsers.add_parser('histogram', help='make histogram from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
#     add_scatter_options( subparsers.add_parser('scatter', help='make scatter from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_extract_options( subparsers.add_parser('extract', help='extract events from image', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_match_options( subparsers.add_parser('match', help='match catalog with simulation', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_simulate_options( subparsers.add_parser('simulation', help='simulate, extract and match', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
    add_addbranches_options( subparsers.add_parser('addbranches', help='add branches to catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )
#     add_ratio_options( subparsers.add_parser('ratio', help='make ratio histograms from catalog', formatter_class = argparse.ArgumentDefaultsHelpFormatter) )

    args = parser.parse_args()
    args.__call__( **vars(args) )
#     del args._func
#     print( colored('using parameters:', 'green', attrs=['bold'] ) )
#     print_var( vars(args).keys(), vars(args), line_char='\t' )

#     with Timer('finished'):
#         try:
#             _func(**vars(args))
#         except TypeError:
#             _func(args)
            
