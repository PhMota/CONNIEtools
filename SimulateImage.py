# coding: utf-8

import numpy as np
import time
import astropy.io.fits
import scipy.stats
import os
import root_numpy
import Statistics as stats

electron_in_eV = 3.745
electron_in_keV = electron_in_eV*1e-3

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

#def round2int2( x, y ):
    #return np.rint( x ).astype(int), np.rint( y ).astype(int)

#def round2int( *x ):
    #return [ np.rint( ix ).astype(int) for ix in x ]

#def diffuseElectrons( x, y, sigma ):
    #sigma = z2sigmaModel['3']( z )
    #xoff, yoff = sigma * scipy.stats.norm.rvs( size = 2*len(sigma) ).reshape( 2, len(sigma) )
    #return x + xoff, y + yoff
    #return np.rint(sigma * scipy.stats.norm.rvs( size = 2*len(sigma) ).reshape( 2, len(sigma) ) ).astype(int)

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
    print 'size', image.val.size
    size = image.val.size
    shape = image.val.shape
    return image + stats.ufloat( scipy.stats.poisson.rvs( lambda_, size = size ).reshape(shape ), unit='e-' )

def add_readout_noise( image, sigma ):
    '''
    add reaout noise noise
    '''
    size = image.val.size
    shape = image.val.shape
    print 'sigma', sigma
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
    print 'gain', gain,
    print 'lambda', lambda_
    print 'sigma', sigma, sigma.asunit('ADU%s'%ohdu)
    print 'shape', shape,
    image = stats.ufloat( np.zeros(shape), unit='e-' )
    print 'zeros shape', image.val.shape, image.unit
    image = add_dark_current(image, lambda_)
    print 'dc shape', image.val.shape, image.unit
    image = diffuse_electrons(image, alpha[ohdu], beta[ohdu] )
    
    image = rebinImage( image, rebin )
    if mbs:
        image -= np.median( image )
    image = add_vOverscan( image, v_os )
    image = add_overscan( image, os )
    image = add_readout_noise( image, sigma )
    print 'shape after rebin', image.val.shape
    print image.__class__.__name__
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
    print 'numberExtracted', len(catalogData['ohdu']), 'to file', out
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


def set_above_threshold_to_nan( image, thr, radius ):
    imagecopy = image.copy()
    s33 = [[1,1,1], [1,1,1], [1,1,1]]
    clusters = scipy.ndimage.label( imagecopy >= thr, structure=s33 )[0]
    distance_from_thr = scipy.ndimage.distance_transform_edt( clusters == 0 )
    imagecopy[ distance_from_thr < radius ] = np.nan
    return imagecopy

def label_clusters( condition_image ):
    s33 = [[1,1,1], [1,1,1], [1,1,1]]
    return scipy.ndimage.label( condition_image, structure=s33 )[0]

def compute_distances( image ):
    return scipy.ndimage.distance_transform_edt( image )

def _extract_cluster_( val, pos ):
    return [val, pos]

class Image:
    def _label_clusters_above_threshold_( self, threshold ):
        return label_clusters( self.image >= thr )
    
    def _compute_distances_to_clusters_( self ):
        return scipy.ndimage.distance_transform_edt(  )

    def _extract_clusters_( self, threshold, border ):
        border *= np.sqrt(2)
        labeled_clusters = label_clusters( self.image > thershold )
        is_cluster = labeled_clusters > 0
        distances_to_cluster = scipy.ndimage.distance_transform_edt( is_cluster == False )
        labeled_clusters = label_clusters( distances_to_cluster <= border )
        
        list_of_clusters = scipy.ndimage.labeled_comprehension( 
            self.image, 
            labeled_clusters,
            index=None, 
            func=lambda v, p: [v, p], 
            out_dtype=list, 
            pass_positions=True 
            )
        levels = scipy.ndimage.labeled_comprehension( self.image, labeled_clusters, index=None, func=lambda v: v, out_dtype=list )
        
    def extract_hits( self ):
        pass
    
class Hits:
    def compute_spectrum():
        return Spectrum( self.E )
        
class Spectrum:
    def plot( self ):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.histogram( self.E, self.bins )
        fig.savefig(self.file)

def simulate_events( shape, ranges, number_of_events ):
    pos = np.random.random( number_of_events*2 )
    
def analysis( shape, ranges, numver_of_events ):
    events = simulate_events( shape, ranges, number_of_events )
    image = events.generate_image( shape )
    image.plot()
    image.fits()
    hits = image.extract_hits( mode = 'cluster', threshold = 15*4, border = 3 )
    hits.catalog()
    hits.plot()
    reconstruction_efficiency = hits.compute_efficiency( events )
    spectrum = hits.compute_spectrum()
    spectrum.plot()
    
