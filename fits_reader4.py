import astropy
import astropy.io
import astropy.io.fits

import numpy as np
import scipy
import scipy.stats
from scipy.misc import factorial
import scipy.optimize

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import glob
import os
import sys
import re
import copy
import root_numpy
import datetime
import time

REMOVE = 0
REMOVEOLD = -1e9

class term:
    #HEADER = '\033[95m'
    #OKBLUE = '\033[94m'
    #OKGREEN = '\033[92m'
    #ARNING = '\033[93m'
    #FAIL = '\033[91m'
    cmd = lambda x: '\033[%dm'%x
    Q = cmd(0)
    B = cmd(1)
    I = cmd(3)
    U = cmd(4)
    #r = '\e[31m'
    r = cmd(91)
    g = cmd(32)
    b = cmd(34)
    
    red = staticmethod( lambda x: term.r + x + term.Q )
    bold = staticmethod( lambda x: term.B + x + term.Q )
    italic = staticmethod( lambda x: term.I + x + term.Q )
    
#folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/scn/merged/*'
folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029*/data_*/scn/merged/*'
folder2 = '/share/storage2/connie/data_analysis/processed02_data/runs/009*/data_*/scn/merged/*'

def listFITS( *patterns ):
    '''
    list all files that match the list of patterns given
    '''
    return sortByrunID([ match for pattern in patterns for match in glob.glob(pattern) if not '-with-' in match ])

def listFITS_run( *runs ):
    return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/scn/merged/*'%run for run in runs ] )

def listFITS_runID( *runIDs ):
    return listFITS( *[ '/share/storage2/connie/data_analysis/processed02_data/runs/*/data_*/scn/merged/*runID_*_%05d_*'%runID for runID in runIDs ] )

def getrunIDFromPath( path ):
    return re.search(r'runID_[0-9]+_([0-9]+)', os.path.basename( path )).groups()[0]

def getrunFromPath( path ):
    return re.search(r'runs/(.+)/data', path ).groups()[0]

def getrun_runIDFromPath( path ):
    return re.search(r'runID_([0-9]+)_([0-9]+)', os.path.basename( path )).groups()

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

def getAssociatedCatalog( *paths ):
    return [ catalog for path in paths for catalog in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/ext/catalog/catalog_data_*.root'%getrunFromPath(path)) if not 'skim' in catalog and not '-with-' in catalog ]

def getAssociatedOSI( *paths ):
    return [ osi for path in paths for osi in glob.glob('/share/storage2/connie/data_analysis/processed02_data/runs/%s/data_*/osi/images/osi_runID_*_%s_*.fits'%(getrunFromPath(path), getrunIDFromPath(path)) ) if not '-with-' in osi ]

def runETA( msg, cmd, eta = None, N = None ):
    print msg, '...'
    if not eta is None and not N is None:
        etastart = time.time()
        eta()
        et = time.time() - etastart
        print 'estimated time', time.strftime("%H:%M:%S", time.gmtime(int(et*N))), '(finishes at %s)'% time.strftime("%H:%M:%S", time.localtime( time.time() + int(et*N)))
    start = time.time()
    res = cmd()
    #print 'elapsed time', time.strftime("%H:%M:%S", time.gmtime(datetime.datetime.time() - start))
    print 'elapsed time', time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
    return res

def removeHitsFromrunID( runID, outputfolder = None, ohdu = None, ROOTfile = None ):
    FITSfile = listFITS_runID( runID )[0]
    outputfilenohits = 'hsub_'+os.path.basename(FITSfile)
    outputfilehits = 'hits_'+os.path.basename(FITSfile)
    
    if not outputfolder is None:
        if not os.path.exists(outputfolder):
            print 'creating directory', outputfolder
            os.makedirs( outputfolder )
    else:
        outputfolder = os.path.dirname(FITSfile)
        print 'outputfolder =', outputfolder
    
    OSIfiles = getAssociatedOSI( FITSfile )
    print 'input FITSfile =', FITSfile
    print 'OSI file ='
    print '\n'.join( OSIfiles )
    
    print 'reading OSI parts files'
    parts = [ astropy.io.fits.open( OSIfile ) for OSIfile in OSIfiles ]
    
    print 'merging parts'
    print 'shape of each part', parts[-1][-1].data.shape
    for ohduindex in range(len(parts[-1])):
        if not parts[-1][ohduindex].data is None:
            parts[-1][ohduindex].data = np.concatenate( [ part[ohduindex].data for part in parts[::-1] ], axis = 0 )
    print 'shape of merged', parts[-1][-1].data.shape
    
    print 'reading SCN file for hits'
    hdulisthits = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        print 'removing ohdus'
        for index in range(len(hdulisthits))[::-1]: #this inversion is needed so the removals work
            ohdu_ = hdulisthits[index].header['OHDU']
            if ohdu_!=ohdu:
                hdulisthits.pop(index)

    
    print 'reading SCN file for no hits'
    hdulistnohits = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        print 'removing ohdus'
        for index in range(len(hdulistnohits))[::-1]: #this inversion is needed so the removals work
            ohdu_ = hdulistnohits[index].header['OHDU']
            if ohdu_!=ohdu:
                hdulistnohits.pop(index)
                parts[-1].pop(index+1)
                    

    if ROOTfile is None:
        ROOTfile = getAssociatedCatalog( FITSfile )[0]
    print 'input ROOTfile =', ROOTfile
    min_max = lambda x: (np.min(x), np.max(x))
    print 'reading start and stop indexes for runID', runID, 'in ROOT file'
    start, stop = min_max( np.argwhere( root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'] == runID ) )
    print 'read from', start, 'to', stop
    
    readcatalog = lambda start_, stop_: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu'], selection='runID==%s'%runID, start=start_, stop=stop_+1 )
    print 'reading ROOT file'
    hitscatalog = runETA( 
        msg = 'reading ROOTfile',
        cmd = lambda: readcatalog(start, stop), #read all entries and return 
        #eta = lambda: readcatalog(start, start+1000),
        #N = (stop-start)/1000
        )

    print 'removing hits'
    for index in range(len(hdulisthits)):
        ohdu_ = hdulistnohits[index].header['OHDU']
        partohdu_ = parts[-1][index+1].header['OHDU']
        hitsohdu = hitscatalog[ hitscatalog['ohdu']==ohdu_ ]
        hdulisthits[index].data[:,:] = REMOVE #-1e9
        hits_xy = np.array( [ [iy,ix] for x,y in zip( hitsohdu['xPix'], hitsohdu['yPix']) for ix,iy in zip(x,y) ] )
        hits_e = np.array( [ ie for e in hitsohdu['ePix'] for ie in e ] )
        if len(hits_e) > 0:
            hdulistnohits[index].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE #-1e9
            parts[-1][index+1].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE
            hdulisthits[index].data[hits_xy[:,0],hits_xy[:,1]] = hits_e
        else:
            print term.red('Warning:'), 'empty hits_e on', ohdu_
    
    hdulisthitsfile = outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilehits
    if os.path.exists(hdulisthitsfile):
        os.remove(hdulisthitsfile)
    print 'writing SCN hits'
    hdulisthits.writeto( hdulisthitsfile )
    
    hdulistnohitsfile =  outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits 
    if os.path.exists(hdulistnohitsfile):
        os.remove(hdulistnohitsfile)
    print 'writing SCN no hits'    
    hdulistnohits.writeto( hdulistnohitsfile )
    
    partsfile = outputfolder + '/merged_' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits 
    if os.path.exists(partsfile):
        os.remove(partsfile)
    print 'writing merged OSI no hits'    
    parts[-1].writeto( partsfile )
    
    print 'created 3 new FITS files'
    print hdulisthitsfile
    print hdulistnohitsfile
    print partsfile

def fit( x, data, func, sel = None, log=False ):
    mask = None
    if sel is None:
        mask = data>0
    else:
        mask = np.all( sel+[ data>0 ], axis=0 )
    pp = scipy.optimize.curve_fit( func, x[mask], data[mask] if not log else np.log(data[mask]) )[0]
    chisq = scipy.stats.chisquare( func( x[mask], *pp ), f_exp = data[mask] if not log else np.log( data[mask] ) )[0]/len(data[mask])
    return pp, chisq

def plotrunID( ohdu, *FITSfiles, **kwargs ):
    fft = False
    if 'fft' in kwargs:
        fft = kwargs['fft']
    plot = True
    if 'plot' in kwargs:
        plot = kwargs['plot']
        
    FITSfiles = FITSfiles[0]
    print 'FITSfile =', FITSfiles
    runID = getrunIDFromPath( FITSfiles[0] )
    print 'runID =', runID
    
    if plot:
        print 'initiate figure'
        fig = plt.figure()
        fig.suptitle('runID%s'%runID+' ohdu%s'%ohdu)
        if fft:
            gs = GridSpec(5,1)
            ax = fig.add_subplot(gs[:-2,0])
            axdiff = fig.add_subplot(gs[-2,0], sharex=ax)
            axfft = fig.add_subplot(gs[-1,0])
        plt.setp(ax.get_xticklabels(), visible=False)
    
    gaussian = lambda x,b,c: c*scipy.stats.norm.pdf(x,0,b)
    gaussianx = lambda x,a,b,c: c*scipy.stats.norm.pdf(x,a,b)
    logGauss = lambda x,b,c: np.log(c/b/np.sqrt(2*np.pi)) - .5*x**2/b**2
    logGaussx = lambda x,a,b,c: np.log(c/b/np.sqrt(2*np.pi)) - .5*(x-a)**2/b**2
    poisson = lambda x,b,c: c*scipy.stats.norm.pdf(x,b**2,b)
    exponential = lambda x,b,c,d,e: c*np.exp(b*(x-d))-e

    result = []
    header = []
    fmt = []
    first = True
    clean = lambda x: x.flatten()[ np.all( [x.flatten()!=REMOVEOLD, x.flatten()!=REMOVE, x.flatten()!=1e10], axis=0) ]

    print 'plotting'
    for FITSfile in FITSfiles:
        if first: 
            header += ['runID']
            fmt += ['%d']
        line = [ int(getrunIDFromPath( FITSfile )) ]
        if first: 
            header += ['ohdu']
            fmt += ['%d']
        line += [ ohdu ]
        
        print 'reading FITS file', FITSfile
        hdulist = astropy.io.fits.open( FITSfile )
        if first: 
            header += ['tempmin']
            fmt += ['%.2f']
        line += [ hdulist[0].header['TEMPMIN' ] ]
        if first: 
            header += ['tempmax']
            fmt += ['%.2f']
        line += [ hdulist[0].header['TEMPMAX' ] ]
        
        print [ hdu.header['OHDU'] for hdu in hdulist ]
        
        data = [ hdu.data for hdu in hdulist if hdu.header['OHDU'] == ohdu ][0]
        print 'data shape', data.shape
        yccd = 4262
        xccd = 4220
        if data.shape[1] == 2*yccd:
            print 'subtracting idle amplifier'
            data = data[0:xccd,0:yccd] - np.flip(data[0:xccd,yccd:2*yccd], 1)
        
        dx = 1
        bins = np.r_[min(clean(data)):max(clean(data)):dx]
        xbins = (bins[1:]+bins[:-1])*.5
        print 'computing histogram of full ccd'
        hist, tmp = np.histogram( clean(data[0:4220,0:4262]), bins = bins )
        #hist2, tmp = np.histogram( clean(data2[0:4220,0:4262]), bins = bins )
        #histR, tmp = np.histogram( clean(data3), bins = bins )

        if first: 
            header += ['pixelcount']
            fmt += ['%d']
        line += [ len(clean(data)) ]

        print 'computing histogram of AC and OS'
        histOS, tmp = np.histogram( clean(data[120:4180,4112+10:4262-10]), bins = bins )
        histAC, tmp = np.histogram( clean(data[120:4180,20:4100]), bins = bins )
        
        print 'medians', np.median(histOS), np.median(histAC)
        
        print 'fitting log(histOS>10) with quadratic function'
        elx, elxchi2 = fit( xbins, histOS, logGaussx, sel=[ histOS>10 ], log=True )
        print 'computing histogram of OSx (corrected by the mean)'
        histOSx, tmp = np.histogram( clean(data[120:4180,4112+10:4262-10])-elx[0], bins = bins )
        print 'fitting log(histOSx>10) with quadratic function'
        el, elchi2 = fit( xbins, histOSx, logGauss, sel=[ histOSx>10 ], log=True )
        if first: 
            header += ['os_std', 'os_chisq']
            fmt += ['%.6g','%.6g']
        line += [el[0], elchi2]
        
        countCut = 10
        
        print 'computing histogram of ACx (corrected by the mean of OS)'
        histACx, tmp = np.histogram( clean(data[120:4180,20:4100])-elx[0], bins = bins )
        print 'fitting log(histAC>10,xbins<0) with quadratic function'
        elneg, elnegchi2 = fit( xbins, histAC, logGauss, sel=[ xbins<0, histAC>10 ], log=True )
        if first: 
            header += ['acNeg_std', 'acNeg_chisq']
            fmt += ['%.6g','%.6g']
        line += [elneg[0], elnegchi2]
        
        print 'fitting log(histAC>10,xbins>0) with quadratic function'
        elpos, elposchi2 = fit( xbins, histAC, logGauss, sel=[ xbins>0, histAC>10 ], log=True )
        if first: 
            header += ['acPos_std', 'acPos_chisq']
            fmt += ['%.6g','%.6g']
        line += [elpos[0], elposchi2]
        if plot:
            print 'generating curves'
            ax.step( xbins, hist, where='mid', label = 'all' )
            #ax.step( xbins, hist2, where='mid', label = 'all2' )
            #ax.step( xbins, histR, where='mid', label = 'allR' )
            ax.step( xbins, histOS, where='mid', label = 'os' )
            ax.step( xbins, histAC, where='mid', label = 'ac' )
            gOS = gaussian(xbins,*el)
            axdiff.step( xbins[gOS>countCut], np.abs(histOS[gOS>countCut] - gOS[gOS>countCut])/gOS[gOS>countCut], where='mid', label = 'os(%.2f)'%(el[0]), color = 'C1' )
            gACneg = gaussian(xbins,*elneg)
            axdiff.step( xbins[gACneg>countCut], np.abs(histAC[gACneg>countCut] - gACneg[gACneg>countCut])/gACneg[gACneg>countCut], where='mid', label = 'ac<(%.2f)'%(elneg[0]), color = 'C2' )
            gACpos = gaussian(xbins,*elpos)
            axdiff.step( xbins[gACpos>countCut], np.abs(histAC[gACpos>countCut] - gACpos[gACpos>countCut])/gACpos[gACpos>countCut], where='mid', label = 'ac>(%.2f)'%(elpos[0]), color = 'C3' )
            #exp_pp, exp_chi2 = fit( xbins, (histAC - gACneg), exponential, sel=[xbins>0, gACneg>countCut] )
            #print 'fit', exp_pp, exp_chi2
            #axdiff.plot( xbins[ np.all([xbins>0, gACneg>countCut], axis=0) ], exponential(xbins[ np.all([xbins>0, gACneg>countCut], axis=0) ], *exp_pp), label = 'fit', color = 'C4' )
            
        #print el[0], el_cr[0], noise[0], np.sqrt(noise[0]**2 - el_cr[0]**2)

        result += [ line ]

        if fft:
            print 'computing Fourier transform'
            maskAC = histAC>100
            fftAC = np.fft.fft( histAC[maskAC] )
            axfft.step( xbins[maskAC], np.fft.fftshift( np.abs(np.real( fftAC )) ), where='mid', label = 'fftac' )
            axfft.step( xbins[maskAC], np.fft.fftshift( np.abs(np.imag( fftAC )) ), where='mid', label = 'fftac' )
            maskOS = histOS>100
            fftOS = np.fft.fft( histOS[maskOS] ) 
            axfft.step( xbins[maskOS], np.fft.fftshift( np.abs(np.real(fftOS)) ), where='mid', label = 'fftos' )
            axfft.step( xbins[maskOS], np.fft.fftshift( np.abs(np.imag(fftOS)) ), where='mid', label = 'fftos' )
            maxfreq = xbins[maskAC][ np.fft.fftshift(np.abs(fftAC)) == np.max(np.abs(fftAC)) ]
            #print 'max freq AC:', np.max(fftAC), xbins[ np.fft.fftshift(fftAC) == np.max(fftAC) ]
            #print 'max freq OS:', np.max(fftOS), xbins[ np.fft.fftshift(fftOS) == np.max(fftOS) ]
            #fftdiff = fftAC - fftOS*np.sum(fftAC)/np.sum(fftOS)
            #histdiff = np.abs( np.fft.ifft( fftdiff ) )
            #histdiff[histdiff<1] = 0
            #ax.step( xbins, histdiff, where='mid', label = 'histdiff' )
        first = False
    if plot:
        print 'setting labels and lims'
        axdiff.set_xlabel('pixel energy [adu]')
        ax.set_ylabel('pixel count')
        ax.set_yscale('log')
        if fft: 
            #axfft.set_yscale('log')
            axfft.set_xlim([maxfreq-100,maxfreq+100])
        ax.set_xlim([-120,120])
        ax.grid(True)
        axdiff.set_yscale('log')
        axdiff.grid(True)
        ax.legend()
        if fft: axfft.legend(loc=2)
        axdiff.legend()
        plt.subplots_adjust(hspace=.1, wspace=.1)
        print 'plot generated in PNG file', 'ccd%s_%s.png'%(runID,ohdu)
        fig.savefig('ccd%s_%s.png'%(runID,ohdu))
    
    print 'table generated in CSV file', 'ccd%s_%s.csv'%(runID,ohdu)
    np.savetxt('ccd%s_%s.csv'%(runID,ohdu), result, header=', '.join(header), fmt=' '.join(fmt), delimiter=', ')

def analysis( run, cflag=True ):
    plot=True
    result = []
    for FITSfile in listFITS_run(run):
        run = getrunFromPath( FITSfile )
        try:
            run2 = int(run)
        except:
            run2 = int(run[:-1])
        runID = getrunIDFromPath( FITSfile )
        print('runID', run, runID )
        
        data = astropy.io.fits.open( FITSfile )
        ohdu2 = None
        #print(data[0].header)
        tempmax = data[0].header['TEMPMAX']
        tempmin = data[0].header['TEMPMIN']
        for datum in data:
            if datum.header['OHDU'] == 2: 
                ohdu2 = datum
                break
        #print(f, os.path.basename(f))
        #BIASSECA= '[4271:4420,1:1055]
        #BIASSECB= '[4121:4270,1:1055]
        #POSTSECB= '[9:4120,1:0]'
        ohdu2.data[ohdu2.data==1e10] = 0
        
        ovscut = 4112
        downmargin=120
        upmargin=4180
        leftmargin=20
        rightmargin=4100
        margin=10
        #print('shape', ohdu2.data.shape)
        data = {'data': {
                #'ccd': ohdu2.data[:,:ovscut],
                #'ccdtrim': ohdu2.data[ downmargin:upmargin, leftmargin:rightmargin],
                #'ovs': ohdu2.data[:,ovscut:],
                'prescan': ohdu2.data[ downmargin:upmargin, 0:leftmargin],
                'ovstrim': ohdu2.data[ downmargin:upmargin, ovscut+margin:-margin],
                'vovs': ohdu2.data[ 10:80, leftmargin:rightmargin ],
                'ccdstrip': ohdu2.data[ downmargin:upmargin, leftmargin:leftmargin+130],
                }
            }

        if cflag:
            catalog = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/ext/catalog/catalog_data_3321_to_3380.root'
            readcatalog = root_numpy.root2array( catalog, treename = 'hitSumm', branches = ['xPix','yPix','ePix'], selection='ohdu==2&&runID==%s'%runID) #, start=3500000, stop=3690000 )
            print('number of events', len(readcatalog['xPix']))
            xmax, ymax, xmin, ymin = 0, 0, 0, 0
            for x,y in zip(readcatalog['xPix'], readcatalog['yPix']):
                xmax = max(xmax, max(x))
                xmin = min(xmin, min(x))
                ymax = max(ymax, max(y))
                ymin = min(ymin, min(y))
            print(xmin, xmax, ymin, ymax)
            
            data['data']['hitstrim'] = np.array( [ ie for x,y,e in zip(readcatalog['xPix'], readcatalog['yPix'], readcatalog['ePix']) for ie in list(e[np.all( [x>leftmargin,x<rightmargin,y>downmargin,y<upmargin], axis=0)]) ] )
            data['data']['hitsstrip'] = np.array( [ ie for x,y,e in zip(readcatalog['xPix'], readcatalog['yPix'], readcatalog['ePix']) for ie in list(e[np.all( [x>leftmargin,x<leftmargin+130,y>downmargin,y<upmargin], axis=0)]) ] )
            
        dx = 1
        bins = np.r_[min(ohdu2.data.flatten()):max(ohdu2.data.flatten()):dx]
        xbins = (bins[1:]+bins[:-1])*.5
        fabove = 1e-2
        fabove2 = 1e-4
        
        data['hist'] = {}
        for key in data['data'].keys():
            data['hist'][key], tmp = np.histogram( data['data'][key].flatten(), bins = bins )
        if cflag:
            data['hist']['nsetrim'] = data['hist']['ccdtrim'] - data['hist']['hitstrim']
            data['hist']['nsestrip'] = data['hist']['ccdstrip'] - data['hist']['hitsstrip']
            data['hist']['ovsrsc'] = data['hist']['ovstrim'] * np.sum(data['hist']['nsestrip'])/np.sum(data['hist']['ovstrim'])
        
        if plot:
            print('plotting...')
            fig = plt.figure()
            fig.suptitle(
                'ohdu2 runID'+runID+'ccd[%s:%s,%s:%s]'%(downmargin,upmargin,leftmargin,rightmargin) + ' ovs[%s:%s,%s:%s]'%(downmargin,upmargin,ovscut+margin,ohdu2.data.shape[1]-margin)+'\ndx=%sadu, fit>%smax'%(dx,fabove) )
            ax = fig.add_subplot(211)
            axzoom = fig.add_subplot(212)
            for axis in [ax,axzoom]:
                for n,key in enumerate(['ovstrim', 'ccdstrip', 'prescan', 'vovs']): #data['hist'].keys()):
                    axis.step( xbins, data['hist'][key], 'C%d-'%n, where='mid', label = key )
                axis.set_yscale('log')
            ax.legend()
            top = max(data['hist']['ovstrim'])
            #ax.set_ylim([1e-8,1e-1])
            ax.set_xlim([-100,100])
            axzoom.set_ylim([ fabove2*top, 1.1*top ])
            
            #axzoom.set_xlim([ min(xbins[left]), 0])
            axzoom.set_xlim([-50,50])
            #axzoom.set_yscale('linear')
            
            fig.savefig('ccd%s.png'%runID)
            plot = False

        data['pp'] = {}
        data['chisquare'] = {}
        data['pplinear'] = {}
        data['chisquarelinear'] = {}
        data['ppleft'] = {}
        data['chisquareleft'] = {}
        data['pplinearleft'] = {}
        data['chisquarelinearleft'] = {}
        
        fitfunc = lambda x,b,c: c*scipy.stats.norm.pdf(x,0,b)
        linfitfunc = lambda x,b,c: c - .5*x**2/b**2
                
        for key in data['hist'].keys():
            if 1:# key == 'ccd' or key == 'ovs':
                above = data['hist'][key] > fabove*max(data['hist'][key])
                
                data['pplinear'][key] = scipy.optimize.curve_fit( linfitfunc, xbins[above], np.log(data['hist'][key][above]) )[0]
                data['chisquarelinear'][key] = scipy.stats.chisquare( linfitfunc( xbins[above], *data['pplinear'][key]), f_exp = np.log(data['hist'][key][above]) )[0]/len(data['hist'][key][above])

                data['pp'][key] = scipy.optimize.curve_fit( fitfunc, xbins[above], data['hist'][key][above], p0=data['pplinear'][key] )[0]
                data['chisquare'][key] = scipy.stats.chisquare( fitfunc( xbins[above], *data['pp'][key]), f_exp = data['hist'][key][above] )[0]/len(data['hist'][key][above])

                left = np.all( [xbins < 0, data['hist'][key] > fabove2*max(data['hist'][key]) ], axis=0 )  #data['pp'][key][0] #std

                data['pplinearleft'][key] = scipy.optimize.curve_fit( linfitfunc, xbins[left], np.log(data['hist'][key][left]) )[0]
                data['chisquarelinearleft'][key] = scipy.stats.chisquare( linfitfunc( xbins[left], *data['pplinear'][key]), f_exp = np.log(data['hist'][key][left]) )[0]/len(data['hist'][key][left])
                
                data['ppleft'][key] = scipy.optimize.curve_fit( fitfunc, xbins[left], data['hist'][key][left], p0=data['pp'][key] )[0]
                data['chisquareleft'][key] = scipy.stats.chisquare( fitfunc( xbins[left], *data['ppleft'][key]), f_exp = data['hist'][key][left] )[0]/len(data['hist'][key][left] )
                print('fit', key, data['pp'][key], data['chisquare'][key], data['ppleft'][key], data['chisquareleft'][key], data['pplinear'][key], data['chisquarelinear'][key], data['pplinearleft'][key], data['chisquarelinearleft'][key])
                
        entry = [ run, runID, 
                tempmin, tempmax,
                data['ppleft']['ovstrim'][0], data['pplinearleft']['ovstrim'][0], data['chisquareleft']['ovstrim'], 
                data['ppleft']['ccdstrip'][0], data['pplinearleft']['ccdstrip'][0], data['chisquareleft']['ccdstrip'],
                ]
        entry += ([ data['ppleft']['nsetrim'][0], data['ppleft']['nsetrim'][1], data['chisquareleft']['nsetrim'] ] if cflag else [])
        print(entry)
        result += [entry]
        np.savetxt('ccd%s.csv'%run, result, header='run, runID, tempmin, tempmax, ovs_std, ovs_stda, ovs_chi2, ccd_std, ccd_stda, ccd_chi2' + ', nse_std, nse_amp, nse_chi2' if cflag else '', fmt='%s',delimiter=', ')

        #break
            
    #print result


def darkcurrent( run ):
    data = np.genfromtxt('ccd%s.csv'%run, delimiter=',', dtype=None)
    cols = zip( *data )
    print cols[1]
    #data = data[data[:,1].argsort()]
    fig = plt.figure()
    gs = GridSpec(1,5)
    ax = fig.add_subplot(gs[0,1:])
    axh = fig.add_subplot(gs[0,0], sharey=ax)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax2 = ax.twinx()
    ax2.step( cols[1], cols[2], 'y--', label='tempmin', where='mid' )
    ax.step( cols[1], cols[4], label='ovs', where = 'mid')
    axh.hist( cols[4], histtype='step', orientation="horizontal" )
    ax.step( cols[1], cols[7], label='ccd', where='mid')
    axh.hist( cols[7], histtype='step', orientation="horizontal" )
    dc = np.sqrt( np.array(cols[7])**2 - np.array(cols[4])**2 )
    mask = np.isfinite(dc)
    ax.step( np.array(cols[1])[mask], dc[mask], label='dc', where='mid')
    axh.hist( dc[mask], histtype='step', orientation="horizontal" )

    dc = np.sqrt( np.array(cols[8])**2 - np.array(cols[5])**2 )
    mask = np.isfinite(dc)
    ax.step( np.array(cols[1])[mask], dc[mask], label='dca', where='mid')
    axh.hist( dc[mask], histtype='step', orientation="horizontal" )
    
    ax.grid(True)
    axh.grid(True)
    ax.legend()
    ax2.legend()
    plt.subplots_adjust(hspace=.1, wspace=.1)
    fig.savefig('dc%s.png'%run)
    
if __name__ == "__main__":
    print term.bold('Dark Current Analysis Tool for the CONNIE colaboration')
    print term.bold('by Philipe Mota (philipe.mota@gmail.com)')
    if len(sys.argv) == 1:
        print 'Usage:'
        print '--runID <runID> --outfolder <path> --removehits'
        print '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in the associated ROOT catalog'
        print
        print '--runID <runID> --outfolder <path> --ROOTfile <path> --removehits'
        print '\tcreate new FITS files from the corresponding runID FITS file subtracting the events in given ROOTfile'
        print 
        print '--ohdu <ohdu> --fft --plot <list of paths>'
        print '\tgenerate plot with FITS files from paths for the given ohdu'
        print
        exit(0)
    runID = None
    outfolder = None
    ohdu = None
    ROOTfile = None
    fft = False
    if '--outfolder' in sys.argv:
        outfolder = sys.argv[sys.argv.index('--outfolder')+1]
        print 'setting outfolder =', outfolder
    if '--ohdu' in sys.argv:
        ohdu = int(sys.argv[sys.argv.index('--ohdu')+1])
        print 'setting ohdu =', ohdu
    if '--runID' in sys.argv:
        runID = int(sys.argv[sys.argv.index('--runID')+1])
        print 'setting runID =', runID
    if '--ROOTfile' in sys.argv:
        ROOTfile = sys.argv[sys.argv.index('--ROOTfile')+1]
        print 'setting ROOTfile =', ROOTfile
    if '--fft' in sys.argv:
        fft = True
        print 'setting fft = True'
        
    if '--removehits' in sys.argv:
        if not runID is None and not outfolder is None:
            removeHitsFromrunID( runID, outputfolder=outfolder, ohdu = ohdu, ROOTfile=ROOTfile )
            exit(0)
        print term.red('error: runID or outputfolder not set')
    if '--plot' in sys.argv:
        if not ohdu is None:
            FITSfiles = sys.argv[sys.argv.index('--plot')+1:]
            print 'input FITSfiles:'
            print '\n'.join(FITSfiles)
            plotrunID( ohdu, FITSfiles, fft=fft )
            exit(0)
        print term.red('error: ohdu not set')
    if '--table' in sys.argv:
        if not ohdu is None:
            FITSfiles = sys.argv[sys.argv.index('--table')+1:]
            print 'input FITSfiles:'
            print '\n'.join(FITSfiles)
            plotrunID( ohdu, FITSfiles, fft=fft, plot=False )
            exit(0)
        print term.red('error: ohdu not set')
