#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
import numpy as np
import time
import astropy.io.fits
import scipy.stats
import scipy.ndimage
import os
import root_numpy
import Statistics as stats
from numpy.lib import recfunctions as rfn

electron_in_eV = 3.745
electron_in_keV = electron_in_eV*1e-3
Cu_energy_eV = 8046
Cu2_energy_eV = 8904
si_energy_eV = 1740

def z_to_sigma( z, alpha, beta ):
    '''
    [beta] = 1
    [alpha] = px
    '''
    return alpha*np.sqrt(-np.log1p(-beta*z))

alpha = {
    2: np.sqrt(891.074)/15, #1.99005750894
    3: np.sqrt(923.666)/15,
    4: np.sqrt(820.692)/15,
    5: np.sqrt(547.614)/15,
    6: np.sqrt(597.398)/15,
    7: np.sqrt(62770)/15,
    8: np.sqrt(476.3)/15,
    9: np.sqrt(565.109)/15,
    10: np.sqrt(763.933)/15,
    11: np.sqrt(1016.970)/15,
    12: np.sqrt(907.9680)/15,
    13: np.sqrt(601.001)/15,
    14: np.sqrt(556.409)/15,
    15: np.sqrt(7060.89)/15,
        }

beta = {
    2: 0.000488203*675, #0.329537025
    3: 0.000441113*675,
    4: 0.00049098*675,
    5: 0.000689114*675,
    6: 0.000633145*675,
    7: 6.31885e-06*675,
    8: 0.000759706*675,
    9: 0.000647537*675,
    10: 0.000494423*675,
    11: 0.000391788*675,
    12: 0.000436115*675,
    13: 0.000625386*675,
    14: 0.00068271*675,
    15: 6.26018e-05*675,
    }

#beta = { k: stats.dfloat(v,'µm') for k,v in beta_.items() }

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

'''
sigma = np.sqrt(-891.074*np.log(1-0.000488203*z))/15
(sigma*15)**2 = -891.074*np.log(1-0.000488203*z)
(sigma*15)**2/(-891) = np.log(1-0.000488203*z)
np.exp( (sigma*15)**2/(-891) ) = 1-0.000488203*z
1-np.exp( (sigma*15)**2/(-891) ) = 0.000488203*z
(1-np.exp( (sigma*15)**2/(-891) ))/0.000488203 = z
'''

sigma2zModel = {
    '3': lambda s: (1-np.exp( (s*15)**2/(-923.666) ) )/0.000441113,
    }

def diffuse_electrons( image, alpha, beta ):
    '''
    diffuse all electrons
    '''
    t = time.time()
    image_in_e = image.val

    N = int(np.sum( image_in_e ))
    z = np.random.rand(N)
    sigma = z_to_sigma(z, alpha, beta)
    indices = np.repeat( np.arange(image_in_e.size), image_in_e.flatten().astype(int) )
    positions = np.unravel_index( indices, image_in_e.shape )
    offsets = np.rint( sigma * scipy.stats.norm.rvs( size = 2*N ).reshape( 2, N ) ).astype(int)
    x, y = positions + offsets
    inside_boudary = np.all( [ x >= 0, x < image_in_e.shape[0], y >= 0 , y < image_in_e.shape[1] ], axis=0 )
    diffused_image = np.zeros(image_in_e.shape)
    np.add.at( diffused_image, [ x[inside_boudary], y[inside_boudary] ], 1. )
    return stats.ufloat( diffused_image, unit='e-' )

def add_dark_current( image, lambda_ ):
    '''
    adds electrons generated pixel by pixel by a poisson distribution and difused
    '''
    print( 'size', image.val.size )
    size = image.val.size
    shape = image.val.shape
    return image + stats.ufloat( scipy.stats.poisson.rvs( lambda_, size = size ).reshape(shape ), unit='e-' )

def add_readout_noise( image, sigma ):
    '''
    add reaout noise noise
    '''
    size = image.val.size
    shape = image.val.shape
    print( 'sigma', sigma )
    return image + stats.ufloat( sigma.val*scipy.stats.norm.rvs( size=size ).reshape( shape ), errSqr = sigma.val**2, unit ='e-')

def rebinImage( image, n ):
    shape = image.val.shape
    rebinned = np.zeros( ( n * int(np.ceil(shape[0]/n)), shape[1] ) )
    rebinned[ :shape[0], :shape[1] ] = image.val
    rebinned.reshape([-1, n, rebinned.shape[1]] ).sum( axis = 1 )
    return stats.ufloat( rebinned, unit=image.unit )

def add_vOverscan( image, n ):
    shape = image.val.shape
    result = np.zeros(( shape[0]+n, shape[1]))
    result[ :shape[0], :shape[1]] = image.val
    return stats.ufloat( result, unit=image.unit )

def add_overscan( image, n ):
    shape = image.val.shape
    result = np.zeros(( shape[0], shape[1]+n))
    result[ :shape[0], :shape[1]] = image.val
    return stats.ufloat( result, unit=image.unit )

def make_image( shape, v_os=70, os=550, ohdu=2, lambda_=0.1, sigma=2, gain=2000., rebin=1, mbs=False ):
    e_gain = gain*electron_in_keV
    print( 'gain', gain )
    print( 'lambda', lambda_ )
    print( 'sigma', sigma, sigma.asunit('ADU%s'%ohdu) )
    print( 'shape', shape )
    image = stats.ufloat( np.zeros(shape), unit='e-' )
    print( 'zeros shape', image.val.shape, image.unit )
    image = add_dark_current(image, lambda_)
    print( 'dc shape', image.val.shape, image.unit )
    image = diffuse_electrons(image, alpha[ohdu], beta[ohdu] )
    
    image = rebinImage( image, rebin )
    if mbs:
        image -= np.median( image )
    image = add_vOverscan( image, v_os )
    image = add_overscan( image, os )
    image = add_readout_noise( image, sigma )
    print( 'shape after rebin', image.val.shape )
    print( image.__class__.__name__ )
    return image
    #return (image*e_gain).astype(int)

def extractParameters( E, dE=1. ):
    bins = np.arange( E.min(), E.max(), dE )
    histogram = np.histogram( E, bins )
    pass

def writeFITS( image, fname ):
    if os.path.exists(fname):
        os.remove(fname)
    astropy.io.fits.writeto( fname, image, header = astropy.io.fits.Header() )
    return fname

def convertImage2Hits( self, image, save = None, verbose = False ):
    return processImage3( image, self.config.extractThreshold )

def standardExtraction( fname, out, verbose = True ):
    os.system( '/share/storage2/connie/data_analysis/Processor/tools/extract/extract.exe %s -c /home/mota/tests/simulatedFITS/extractConfigFS2.xml -o %s'%(fname,out) + ( '>/dev/null' if not verbose else '') )
    catalogData = root_numpy.root2array( out, treename='hitSumm' )
    print( 'numberExtracted', len(catalogData['ohdu']), 'to file', out )
    #ohduData = catalogData[ catalogData['ohdu'] == int(self.config.ohdu) ]

def makeHitSumm( **kwargs ):
    obj = hitSumm( **kwargs )
    return obj

class SimulateImage:
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
        print( '+Cu' )
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
        print( '%.2fs'%(timeit.default_timer() - t) )
        return image
        
    def getNeutrinosImage_e( self, save = False ):
        self.partial_id += self.neutrino_id( self.config.numberNeutrinos )
        print( '+nu' )
        t = timeit.default_timer()

        path = self.config.outPath + self.neutrino_id + self.image_ext
        if os.path.exists( path ):
            print( '(readnu)' )
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
        print( '%.2fs'%(timeit.default_timer() - t) )
        return image

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
        scnFile = get_scn_path( int(runID) )
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
        #print 'runID', runID, 'ohdu', ohdu,
        rawFiles = get_osi_path( int(runID) )
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        image = np.concatenate( [ ohdus(astropy.io.fits.open(rawFile)) for rawFile in rawFiles[::-1] ], axis=0 )
        return image

    @staticmethod
    def _getMBSImage_adu( runID, ohdu ):
        #print 'runID', runID, 'ohdu', ohdu,
        rawFiles = get_mbs_path( int(runID) )
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        image = np.concatenate( [ ohdus(astropy.io.fits.open(rawFile)) for rawFile in rawFiles[::-1] ], axis=0 )
        return image

    @staticmethod
    def list_ohdus( runID ):
        rawFiles = get_raw_path( int(runID) )
        ohdus = [ int(fit.header['OHDU']) for fit in astropy.io.fits.open(rawFiles[0]) if fit.data is not None ]
        return ohdus

    def getRawImage_adu( self, runID=None, save = False, ohdu=None ):
        if ohdu is None: ohdu = self.config.ohdu
        if runID is None: runID = self.config.runID
        #print 'runID', runID, 'ohdu', ohdu,
        path = self.config.outPath + self.raw_id + str(runID) + self.image_ext
        if os.path.exists( path ):
            return readImage( path )
        rawFiles = get_raw_path( int(runID) )
        return hitSumm._getRawImage_adu( runID, ohdu )
    
    @staticmethod
    def _getSCNImage_adu( runID, ohdu ):
        #print 'runID', runID, 'ohdu', ohdu,
        scnFile = get_scn_path( int(runID) )[0]
        ohdus = lambda fits: [ fit.data for fit in fits if int(fit.header['OHDU']) == ohdu ][0]
        image = ohdus(astropy.io.fits.open(scnFile))
        return image

    def getSCNImage_adu( self, runID=None, ohdu=None ):
        if runID is None: runID = self.config.runID
        if ohdu is None: ohdu = self.config.ohdu
        #print 'runID', runID, 'ohdu', ohdu,
        path = self.config.outPath + self.scn_id + str(runID) + self.image_ext
        if os.path.exists( path ):
            return readImage( path )
        scnFile = get_scn_path( int(runID) )[0]
        return hitSumm._getSCNImage_adu( runID, ohdu )

    def getSCNHits_adu( self, runID, ohdu=None, save = False ):
        self.partial_id += self.scn_id + str(runID)
        print( '+SCN' )
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
        print( '%.2fs'%(timeit.default_timer() - t) )
        return hits_adu

    def getMBSubtracted_adu( self, image ):
        median = np.median( image.flatten() )
        image -= median
        return image
    
    def filterExposure( self, hits ):
        self.partial_id += self.exposure_id
        print( '+exp' )
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


        
def simulate_events( args ):
    shape = args['shape']
    depth_range = args['depth_range']
    charge_range = args['charge_range']
    number_of_charges = args['number_of_charges']
    array_of_positions = np.random.random( number_of_charges*2 ).reshape(-1,2)*(np.array(shape)-1)
    array_of_depths = np.random.random( number_of_charges )*(depth_range[1] - depth_range[0]) + depth_range[0]
    array_of_charges = np.random.random( number_of_charges )*(charge_range[1] - charge_range[0]) + charge_range[0]
    array_of_identities = ['random']*number_of_charges
    if args['number_of_Cu_charges'] > 0:
        number_of_Cu_charges = args['number_of_Cu_charges']
        array_of_positions = np.concatenate( (array_of_positions, np.random.random( number_of_Cu_charges*2 ).reshape(-1,2)*(np.array(shape)-1) ), axis=0 )
        array_of_depths = np.append( array_of_depths, np.random.random( number_of_Cu_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] )
        array_of_charges = np.append( array_of_charges, [Cu_energy_eV/electron_in_eV]*number_of_Cu_charges )
        array_of_identities.extend( ['Cu']*number_of_Cu_charges )
        print( 'Cu', array_of_charges )
    print( 'total charges created', len( array_of_identities ) )
    return Events( array_of_positions, array_of_depths, array_of_charges, array_of_identities, args )

class Events( np.recarray ):
    def __new__(cls, array_of_positions, array_of_depths, array_of_charges, array_of_identities, args ):
        dtype = [
            ('x', float),
            ('y', float),
            #('z', float),
            ]
        obj = np.array( [ tuple(pos) for pos in array_of_positions], dtype = dtype ).view(np.recarray).view(cls)
        obj = rfn.append_fields( obj, 
                                 'z',
                                 array_of_depths, 
                                 dtypes=(float), 
                                 asrecarray=True 
                                 ).view(cls)
        obj = rfn.append_fields( obj, 
                                 'q',
                                 array_of_charges, 
                                 dtypes=(int), 
                                 asrecarray=True 
                                 ).view(cls)
        obj = rfn.append_fields( obj, 
                                 'id',
                                 array_of_identities,
                                 dtypes=None,
                                 asrecarray=True 
                                 ).view(cls)
        obj.xyshape = np.array(args['shape'])
        obj.rebin = np.array(args['rebin'])
        obj.xyrebinshape = obj.xyshape/obj.rebin
        obj.charge_gain = args['charge_gain']
        obj.number_of_Cu_charges = args['number_of_Cu_charges']
        obj.number_of_Cu2_charges = args['number_of_Cu2_charges']
        obj.number_of_Si_charges = args['number_of_Si_charges']
        obj.diffusion_function = args['diffusion_function']
        obj.charge_efficiency_function = args['charge_efficiency_function']
        obj.lambda_ = args['dark_current']
        obj.sigma = args['readout_noise']
        obj.count = len( obj.q )
        print( 'ids', obj.id )
        return obj
    
    def get_charge(self):
        return ( self.charge_efficiency_function(self.z) * self.q ).astype(int)
    
    def get_xy( self ):
        return np.array( (self.x, self.y) ).T

    def get_xy_rebin( self ):
        return (np.array( (self.x, self.y) )/self.rebin[:,None]).T

    def get_E( self ):
        return self.q*self.charge_gain
    
    def get_id_code( self ):
        map_id = {'random': 1, 'Cu': 11, 'Cu2': 12, 'Si':10}
        return map( lambda i: map_id[i], self.id )

    def get_sigma( self ):
        return self.diffusion_function( self.z )
        

def generate_image( events ):
    from Image import Image
    
    if events.count > 0:
        sigma_per_event = events.get_sigma()
        charge_per_event = events.get_charge()
        Q = np.sum(charge_per_event)
        xy_norm = scipy.stats.norm.rvs( size = 2*Q ).reshape( -1, 2 )
        sigma_per_charge = np.repeat( sigma_per_event, charge_per_event )
        xy_per_charge = np.repeat( events.get_xy(), charge_per_event, axis=0 )
        xy = xy_per_charge + sigma_per_charge[:,None] * xy_norm
        bins = [
            np.arange(events.xyshape[0]+1),
            np.arange(events.xyshape[1]+1)
            ]
        image = np.histogramdd( xy, bins = bins )[0] * events.charge_gain
    else:
        image = np.zeros( events.xyshape )
    print( 'shape', image.shape )
        
    if events.lambda_ > 0:
        dark_current = scipy.stats.poisson.rvs( events.lambda_, size = events.xyshape[0]*events.xyshape[1] ).reshape(events.xyshape)
        image += dark_current * events.charge_gain
    
    image = image.reshape( (events.xyrebinshape[0], events.xyrebinshape[1], -1) ).sum(axis = -1)
    print( 'shape rebin', image.shape, events.xyrebinshape )
    
    if events.sigma > 0:
        noise = scipy.stats.norm.rvs( size=events.xyrebinshape[0]*events.xyrebinshape[1] ).reshape(events.xyrebinshape)*events.sigma
        image += noise

    return Image( image )

def integral_norm_pixel( pix_edge, mu, sigma ):
    sqrt2 = np.sqrt(2.)
    return -.5*( scipy.special.erf( -( pix_edge+1 - mu)/( sqrt2*sigma )) - scipy.special.erf( -( pix_edge - mu )/( sqrt2*sigma ) ) )

def probability_pixel( xPix, yPix, mu, sigma, single_sigma = False ):
    if single_sigma: sigma[1] = sigma[0]
    return integral_norm_pixel( xPix, mu = mu[0], sigma = sigma[0] ) * integral_norm_pixel( yPix, mu = mu[1], sigma = sigma[1] )

def pdf( ePix, xPix, yPix, E, mu, sigma, sigma_noise, single_sigma = False):
    prob = probability_pixel( xPix, yPix, mu, sigma, single_sigma )
    prob_i = prob[:,None]
    ePix_i = ePix[:,None]
    q_j = np.arange(E)[None,:]
    j = -1
    return np.sum( scipy.stats.binom.pmf( q_j, E, prob_i ) * scipy.stats.norm.pdf( ePix_i - q_j, scale = sigma_noise ), axis = j )

def pdf_sigma0( ePix, xPix, yPix, E, mu, sigma, single_sigma = False):
    prob = probability_pixel( xPix, yPix, mu, sigma, single_sigma )
    return scipy.stats.binom.pmf( ePix, E, prob )

def negloglikelihood( ePix, xPix, yPix, E, mu, sigma, sigma_noise, single_sigma = False ):
    val = pdf( ePix, xPix, yPix, E, mu, sigma, sigma_noise, single_sigma = single_sigma )
    val = val[ val > 0 ]
    return -np.nansum( np.log( val ) )

binom_pmf = scipy.stats.binom.pmf
norm_pdf = scipy.stats.norm.pdf
np_sum = np.sum
sqrt2 = np.sqrt(2.)
erf = scipy.special.erf
def negloglikelihood_fast( ePix, xPix, yPix, E, mu, sigma, sigma_noise, single_sigma = False ):
    integral_norm_pixel_x = -.5*( erf( -( xPix+1 - mu[0])/( sqrt2*sigma[0] )) - erf( -( xPix - mu[0] )/( sqrt2*sigma[0] ) ) )
    integral_norm_pixel_y = -.5*( erf( -( yPix+1 - mu[1])/( sqrt2*sigma[1] )) - erf( -( yPix - mu[1] )/( sqrt2*sigma[1] ) ) )
    
    prob = integral_norm_pixel_x * integral_norm_pixel_y
    prob_i = prob[:,None]
    ePix_i = ePix[:,None]
    q_j = np.arange(E).astype(int)[None,:]
    j = -1
    val = np_sum( binom_pmf( q_j, E, prob_i ) * norm_pdf( ePix_i - q_j, scale = sigma_noise ), axis = j )
    val = val[ val > 0 ]
    return -np.nansum( np.log( val ) )

def size_like( ePix, xPix, yPix, E0, mu0, sigma0, sigma_noise0, single_sigma = False, fast = False, tol = None ):
    #print( 'size like call', len(ePix) )
    mu = [0,0]
    sigma = [0,0]
    Q = 0
    e = ePix.astype(int)
    x = xPix.astype(int)
    y = yPix.astype(int)
    positive_e = e > 0
    e = e[ positive_e ]
    x = x[ positive_e ]
    y = y[ positive_e ]
    if fast:
        fun = lambda p: negloglikelihood_fast(e, x, y, int(p[0]), [p[1], p[2]], [p[3], p[4]], sigma_noise0, single_sigma)
    else:
        fun = lambda p: negloglikelihood(e, x, y, int(p[0]), [p[1], p[2]], [p[3], p[4]], sigma_noise0, single_sigma)
    
    Q, mu[0], mu[1], sigma[0], sigma[1] = scipy.optimize.minimize(
        fun = fun,
        x0 = [ int(E0), mu0[0]+.5, mu0[1]+.5, sigma0[0], sigma0[1] ],
        tol = tol,
        bounds = [(1, np.inf), (-np.inf,np.inf), (-np.inf, np.inf), (0.001, np.inf), (0.001, np.inf)] ).x
    
    print( 'size like call', (Q, E0), (mu[0], mu[1], mu0), (sigma[0], sigma[1], sigma0) )
    return Q, mu, sigma



class Hits(np.recarray):
    
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
    
def get_indices( array, length ):
    return np.floor(array/length).astype(int)

class Matches( np.recarray ):
    def __new__( cls, hits, events, lvl = 3, length = 3 ):
        matches = match( hits, events, lvl,  length )
        obj = np.core.records.fromrecords( [ tuple(m) for m in matches], dtype=[('rec',int),('sim',int),('dist',float),('mult',bool)] ).view(cls)
        obj.n_events = len(events)
        obj.n_hits = len(hits)
        return obj
    
    def get_unmatched( self, t ):
        all_ = None
        if t == 'sim':
            all_ = set(range(self.n_events))
        elif t == 'rec':
            all_ = set(range(self.n_hits))
        return np.array( list( all_ - set( tuple( self[t] ) ) ) ).astype(int)
    
    def get_matched_indices(self, t):
        return np.unique( self[t] ).astype(int)

    def get_single_matched(self, t):
        return self[t][ self.mult == False ]

    def get_mult_matched(self, t):
        return self[t][ self.mult == True ]
    
def match( hits, events, lvl = 3, length = 4 ):
    import itertools
    linked_list = {}
    radius = 1
    shifts = range( -radius, radius+1 )
    listofshifts = [ shifts for i in range(2) ]
    nshifts = tuple(itertools.product(*listofshifts) )

    array_of_projections_events = events.get_xy()
    array_of_projections_indices_events = get_indices( array_of_projections_events, length=length )
    events_indices_set = set( map(tuple, array_of_projections_indices_events ) )
    
    for i, ni in enumerate(array_of_projections_indices_events):
        for ns in nshifts:
            try:
                linked_list[ tuple(ni + ns) ].append(i)
            except:
                linked_list[ tuple(ni + ns) ] = [i]

    array_of_projections_hits = np.array( (hits['xBary%d' % lvl], hits['yBary%d' % lvl]) ).T
    array_of_projections_indices_hits = get_indices( array_of_projections_hits, length=length )
    
    matched = []
    dist = lambda a, b: np.sqrt( np.sum( (a-b)**2, axis=-1 ) )
    for hit_ind, hit_ni in enumerate(array_of_projections_indices_hits):
        try:
            event_ind = np.array(linked_list[tuple( hit_ni )])
        except:
            continue
        pos_hit = array_of_projections_hits[ hit_ind ]
        std = np.sqrt( ( hits[ hit_ind ]['xVar%d' % lvl], hits[ hit_ind ]['yVar%d' % lvl] ) )
        pos_event = array_of_projections_events[event_ind]

        ds = dist( pos_hit/std, pos_event/std )
        is_matched = ds <= 1
        event_ind_matched = event_ind[is_matched]
        if len( event_ind_matched ) == 0:
            continue
        is_mult = False
        if len( event_ind_matched ) > 1:
            is_mult = True
        for ei, d in zip(event_ind_matched, ds[is_matched]):
            matched.append( (hit_ind, ei, d, is_mult) )

    matched = np.array(matched)
    return matched

def show_reconstruction( output, events, image, hits, level = 1, zoom = None ):
    from matplotlib import pylab as plt
    from matplotlib.patches import Ellipse
    from matplotlib.colors import LogNorm
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = image.T

    ax.imshow( im, origin='lower', cmap='Blues', norm=LogNorm(vmin=0.01, vmax= im.max() ) )
    x, y = events.x-.5, events.y-.5
    ax.scatter( x, y, c='r', s=1 )
    
    first = True
    for hit in hits:
        width = np.sqrt(hit['xVar%d' % level])*5
        height = np.sqrt(hit['yVar%d' % level])*5
        ax.add_artist( Ellipse((hit['xBary%d' % level], hit['yBary%d' % level]), width, height, fill=False ) )
    if zoom:
        ax.set_xlim(zoom[0])
        ax.set_ylim(zoom[1])
    fig.savefig(output)

def spectrum_reconstruction( output, events, image, hits, lvl, tol, binsize ):
    from matplotlib import pylab as plt
    from Timer import Timer
    
    sim_E = events.get_E()
    rec_E = hits[ 'E%d' % lvl ]

    with Timer('size like fast') as t:
        minimization = hits.get_sizeLike( lvl, events.gain, events.sigma, fast = True, tol = tol )

    #with Timer('size like') as t:
        #minimization_fast = hits.get_sizeLike(1, events.gain, events.sigma, fast = False )
    
    size_like_E = minimization[0]
    size_like = minimization[2]
    rec_sigma = size_like[:,0]
    matches = Matches( hits, events, lvl=3, length = hits.border+2 )

    min_E = min( np.min(sim_E), np.min(rec_E) )
    max_E = max( np.max(sim_E), np.max(rec_E) )
    bins = np.arange( min_E, max_E+10, binsize )
    with Timer('histogram E') as t:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.hist( sim_E, bins, histtype='step', label = 'sim' )
        ax.hist( sim_E[ matches.get_unmatched('sim') ], bins, histtype='step', label = 'sim not rec' )
        ax.hist( rec_E, bins, histtype='step', label = 'all rec' )
        ax.hist( rec_E[ matches.get_single_matched('rec')], bins, histtype='step', label = 'rec single' )
        ax.hist( rec_E[ matches.get_mult_matched('rec')], bins, histtype='step', label = 'rec mult' )
        ax.hist( rec_E[ matches.get_unmatched('rec')], bins, histtype='step', label = 'rec fake' )
        
        ax.hist( size_like_E, bins, histtype='step', label = 'size like' )
        ax.set_xlabel(r'$E$[ADU]')
        ax.axvline( hits.threshold, ymin=0, ymax=ax.get_ylim()[1], c='g' )
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        fig.savefig( output, bbox_extra_artists = (legend,), bbox_inches='tight')

    with Timer('scatter') as t:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.scatter( sim_E[ matches.get_single_matched('sim') ], rec_E[ matches.get_single_matched('rec') ], label = r'$E_{\rm rec}$s', s=1 )
        ax.scatter( sim_E[ matches.get_mult_matched('sim') ], rec_E[ matches.get_mult_matched('rec') ], label = r'$E_{\rm rec}$m', s=1 )
        ax.scatter( sim_E[ matches.sim ], size_like_E[ matches.rec ], label = r'$E_{\rm sl}$', s=1 )
        
        ax.plot( (ax.get_xlim()[0], ax.get_xlim()[1]), (ax.get_xlim()[0], ax.get_xlim()[1]), 'k' )
        ax.set_xlabel(r'$E_{\rm sim}$')
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        output2 = (lambda _: _[0]+'sc.'+_[1])(output.split('.'))
        fig.savefig( output2, bbox_extra_artists = (legend,), bbox_inches='tight')
        
    sim_sigma = events.get_sigma()
    rec_std = np.sqrt( hits['xVar%d' % lvl] )
    
    with Timer('scatter') as t:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.scatter( sim_E, sim_sigma, label = r'$\sigma_{\rm sim}$', s=1 )
        ax.scatter( rec_E, rec_std, label = r'${\rm std}_{\rm rec}$', s=1 )
        ax.scatter( rec_E, rec_sigma, label = r'$\sigma_{\rm rec}$', s=1 )

        ax.set_xlabel(r'$E$')
        ax.set_xlim( (sim_E.min(), sim_E.max()*1.2) )
        ax.set_ylim( (sim_sigma.min(), sim_sigma.max()*1.2) )
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        output2 = (lambda _: _[0]+'sigma.'+_[1])(output.split('.'))
        fig.savefig( output2, bbox_extra_artists = (legend,), bbox_inches='tight')

    binsize = .02
    min_sigma = np.min(sim_sigma)
    max_sigma = np.max(sim_sigma)*1.5
    bins = np.arange( min_sigma, max_sigma+1, binsize )

    with Timer('scatter') as t:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist( sim_sigma, bins=bins, label = r'$\sigma_{\rm sim}$', histtype = 'step' )
        ax.hist( rec_std, bins=bins, label = r'${\rm std}_{\rm rec}$', histtype = 'step' )
        ax.hist( rec_sigma, bins=bins, label = r'$\sigma_{\rm rec}$', histtype = 'step' )
        ax.set_xlabel(r'$\sigma$')
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        output2 = (lambda _: _[0]+'sigma2.'+_[1])(output.split('.'))
        fig.savefig( output2, bbox_extra_artists = (legend,), bbox_inches='tight')
    
    with Timer('scatter') as t:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter( sim_sigma[ matches.sim ], rec_std[ matches.rec ], label = r'${\rm std}_{\rm rec}$', s=1 )
        ax.scatter( sim_sigma[ matches.sim ], rec_sigma[ matches.rec ], label = r'$\sigma_{\rm rec}$', s=1 )

        ax.plot( (ax.get_xlim()[0], ax.get_xlim()[1]), (ax.get_xlim()[0], ax.get_xlim()[1]), 'k' )
        ax.set_xlabel(r'$\sigma_{\rm sim}$')
        ax.set_ylim( (sim_sigma.min(), sim_sigma.max()*1.2) )
        legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper left' )
        output2 = (lambda _: _[0] + 'scsigmaPair.' + _[1])(output.split('.'))
        fig.savefig( output2, bbox_extra_artists = (legend,), bbox_inches='tight')

def analysis( args ):
    from Timer import Timer
    import Catalog

    with Timer('analysis') as t:
        if args['image_fits_input'] == 'none':
            with Timer('events generated') as t:
                events = simulate_events( args )
            
            with Timer('image generated') as t:
                image = generate_image( events )

            #if args['image_energy_spectrum'] != 'none':
                #with Timer('plot saved at ' + args['image_energy_spectrum'] ) as t:
                    #image.spectrum( args['image_energy_spectrum'], events.gain, args['dark_current'], args['readout_noise'] )

            if args['image_fits_output'] != 'none':
                with Timer('fits image saved at ' + args['image_fits_output'] ) as t:
                    image.save_fits( args['image_fits_output'] )
        #else:
            #image = read_image_from_file( args['image_fits_input'] )
            
        with Timer('hits extracted') as t:
            hits = image.extract_hits( mode = 'cluster', threshold = args['extraction_threshold'], border = args['extraction_border'] )
            n = args['extraction_border']+1
            hits = hits.add_number_of_pixel_levels( range(n) )
            hits = hits.add_energy_levels( range(n) )
            hits = hits.add_barycenter_levels( range(n) )
            hits = hits.add_variance_levels( range(n) )
            print( 'hits', hits.dtype.names )

        if args['image_fits_input'] == 'none':
            hits = hits.match_simulated_events( events, lvl=1, length=3 )

        if args['catalog_root_output'] != 'none':
            with Timer('root catalog at ' + args['catalog_root_output'] ) as t:
                Catalog.build_new_catalog_from_recarray( hits, args['catalog_root_output'], None )
        
        #with Timer('size like') as t:
            #hits.add_sizelike( 3, args['charge_gain'], args['readout_noise'], fast = False, tol = 1e-3 )
            #t( len(hits) )

        #with Timer('write catalog sizelike') as t:
            #Catalog.build_new_catalog_from_recarray( hits, 'catalog_test.root', None )
            
        if args['reconstruction_image'] != 'none':
            with Timer('plot reconstruction') as t:
                show_reconstruction( args['reconstruction_image'], events, image, hits, zoom=[[0,200],[0,200]] )
            
        #if args['reconstruction_spectra'] != 'none':
            #with Timer('plot spectrum reconstruction') as t:
                #spectrum_reconstruction( args['reconstruction_spectra'], events, image, hits, args['sigma_like_level'], args['sigma_like_tol'], binsize = 50 )


if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description = 'simulate events, generates image and extract hits',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        )
    def tuple_of_int( x ):
        return map( int, eval(x.replace('\"', '')) )
    
    parser.add_argument('number_of_charges', type=int, default = '4000', help = 'number of charges to be randomly generated' )
    parser.add_argument('--number-of-Cu-charges', type=int, default = '0', help = 'number of charges to be randomly generated at the Copper fluorescence energy 8.046keV' )
    parser.add_argument('--number-of-Cu2-charges', type=int, default = '0', help = 'number of charges to be randomly generated at the secundary Copper fluorescence energy 8.904keV' )
    parser.add_argument('--number-of-Si-charges', type=int, default = '0', help = 'number of charges to be randomly generated at the Silicon fluorescence energy 1.740keV' )
    parser.add_argument('--shape', type=tuple_of_int, default = '\"[4000,4000]\"', help = 'shape of the image in pixel per pixel per µm' )
    parser.add_argument('--rebin', type=tuple_of_int, default = '\"[1,1]\"', help = 'shape of the image in pixel per pixel per µm' )
    parser.add_argument('--charge-range', type=tuple_of_int, default = '\"[5,200]\"', help = 'range into which to randomly generate charges' )
    parser.add_argument('--depth-range', type=tuple_of_int, default = '\"[0,670]\"', help = 'range into which to randomly generate depths' )
    parser.add_argument('--charge-gain', type=eval, default = '7.25', help = 'factor to convert charges into ADU' )
    parser.add_argument('--readout-noise', type=eval, default = '0', help = 'sigma of the normal noise distribution in ADU' )
    parser.add_argument('--dark-current', type=eval, default = '0', help = 'lambda of Poisson distribution dimensionless' )
    parser.add_argument('--extraction-threshold', type=eval, default = '15*4', help = 'energy threshold for extraction in ADU' )
    parser.add_argument('--extraction-border', type=int, default = '3', help = 'borders to be added around the axtracted event' )
    parser.add_argument('--image-fits-output', type=str, default = 'simulation.fits', help = 'set to "none" not to generate a fits output' )
    parser.add_argument('--image-fits-input', type=str, default = 'none', help = 'set to <path> to load existing fits image' )
    parser.add_argument('--extract', type=bool, default = True, help = 'set to False to skip the extraction' )
    parser.add_argument('--image-energy-spectrum', type=str, default = 'image_energy_spectrum.png', help = 'set to "none" not to plot image energy spectrum' )
    parser.add_argument('--catalog-root-output', type=str, default = 'catalog_test.root', help = 'set to "none" not to generate root catalog' )
    parser.add_argument('--reconstruction-image', type=str, default = 'reconstruction_image.pdf', help = 'set to "none" not to plot the reconstruction image' )
    parser.add_argument('--reconstruction-spectra', type=str, default = 'reconstruction_spectra.png', help = 'set to "none" not to plot the reconstruction spectra' )
    parser.add_argument('--sigma-like-level', type=int, default = '2', help = 'level used to compute the sigma like' )
    parser.add_argument('--sigma-like-tol', type=float, default = '1e-3', help = 'tolerance used to compute the sigma like' )
    parser.add_argument('--diffusion-function',
                        type=str, 
                        default = 'sqrt(-923.666*log1p(-0.000441113*z))/15 if z < 670 else 0',
                        help = 'function to map z-depth into sigma' 
                        )
    parser.add_argument('--charge-efficiency-function',
                        type=str, 
                        default = '1. if z < 670 else .9',
                        help = 'function to map z-depth into sigma' 
                        )
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )
    from numpy import *
    args.diffusion_function = eval( 'vectorize(lambda z: %s)' % args.diffusion_function )
    args.charge_efficiency_function = eval( 'vectorize(lambda z: %s)' % args.charge_efficiency_function )

    analysis( vars(args) )
