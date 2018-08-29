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
            os.makedirs( outputfolder )
    else:
        outputfolder = os.path.dirname(FITSfile)
    
    hdulisthits = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        for index in range(len(hdulisthits))[::-1]: #this inversion is needed so the removals work
            ohdu_ = hdulisthits[index].header['OHDU']
            if ohdu_!=ohdu:
                #print 'ohdu', ohdu_, 'removed'
                hdulisthits.pop(index)

    hdulistnohits = astropy.io.fits.open(FITSfile)
    if not ohdu is None:
        for index in range(len(hdulistnohits))[::-1]: #this inversion is needed so the removals work
            ohdu_ = hdulistnohits[index].header['OHDU']
            if ohdu_!=ohdu:
                #print 'ohdu', ohdu_, 'removed'
                hdulistnohits.pop(index)
                    
    print 'input FITSfile =', FITSfile
    if ROOTfile is None:
        ROOTfile = getAssociatedCatalog( FITSfile )[0]
    print 'input ROOTfile =', ROOTfile
    #runIDs = root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID']
    
    readcatalog = lambda stop: root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu'], selection='runID==%s'%runID, start=0, stop=stop )
    hitscatalog = runETA( 
        msg = 'reading ROOTfile',
        cmd = lambda: readcatalog(None), #read all entries and return 
        eta = lambda: readcatalog(10000),
        N = len(root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['runID'] )['runID'])/10000
        )

    for index in range(len(hdulisthits)):
        ohdu_ = hdulistnohits[index].header['OHDU']
        hitsohdu = hitscatalog[ hitscatalog['ohdu']==ohdu_ ]
        hdulisthits[index].data[:,:] = REMOVE #-1e9
        hits_xy = np.array( [ [iy,ix] for x,y in zip( hitsohdu['xPix'], hitsohdu['yPix']) for ix,iy in zip(x,y) ] )
        hits_e = np.array( [ ie for e in hitsohdu['ePix'] for ie in e ] )
        hdulistnohits[index].data[hits_xy[:,0],hits_xy[:,1]] = REMOVE #-1e9
        hdulisthits[index].data[hits_xy[:,0],hits_xy[:,1]] = hits_e
    
    hdulisthitsfile = outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilehits
    if os.path.exists(hdulisthitsfile):
        os.remove(hdulisthitsfile)
    hdulisthits.writeto( hdulisthitsfile )
    
    hdulistnohitsfile =  outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits 
    if os.path.exists(hdulistnohitsfile):
        os.remove(hdulistnohitsfile)
    
    hdulistnohits.writeto( hdulistnohitsfile )
    print 'created 2 new FITS files'
    print outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilehits
    print outputfolder + '/' + ( '' if ohdu is None else 'ohdu%d_'%ohdu ) + outputfilenohits

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
    runID = getrunIDFromPath( FITSfiles[0] )
    if plot:
        fig = plt.figure()
        fig.suptitle('runID%s'%runID+' ohdu%s'%ohdu)
        gs = GridSpec(4,1)
        ax = fig.add_subplot(gs[:-1,0])
        if fft: axfft = ax.twinx()
        axdiff = fig.add_subplot(gs[-1,0], sharex=ax)
        plt.setp(ax.get_xticklabels(), visible=False)
    
    gaussian = lambda x,b,c: c*scipy.stats.norm.pdf(x,0,b)
    logGauss = lambda x,b,c: np.log(c/b/np.sqrt(2*np.pi)) - .5*x**2/b**2
    poisson = lambda x,b,c: c*scipy.stats.norm.pdf(x,b**2,b)

    result = []
    header = []
    first = True
    clean = lambda x: x.flatten()[ np.all( [x.flatten()!=REMOVEOLD, x.flatten()!=REMOVE, x.flatten()!=1e10], axis=0) ]

    for FITSfile in FITSfiles:
        if first: header += ['runID']
        line = [ int(getrunIDFromPath( FITSfile )) ]
        
        hdulist = astropy.io.fits.open( FITSfile )
        data = [ hdu.data for hdu in hdulist if hdu.header['OHDU'] == ohdu ][0]
        dx = 1
        bins = np.r_[min(clean(data)):max(clean(data)):dx]
        xbins = (bins[1:]+bins[:-1])*.5
        hist, tmp = np.histogram( clean(data), bins = bins )
        
        histOS, tmp = np.histogram( clean(data[120:4180,4112+10:-10]), bins = bins )
        histAC, tmp = np.histogram( clean(data[120:4180,20:4100]), bins = bins )
        
        el, elchi2 = fit( xbins, histOS, logGauss, sel=[ histOS>10 ], log=True )
        if first: header += ['os_std', 'os_chisq']
        line += [el[0], elchi2]
        
        countCut = 10
        
        el_cr, el_crchi2 = fit( xbins, histAC, logGauss, sel=[ xbins<0, histAC>10 ], log=True )
        if first: header += ['acNeg_std', 'acNeg_chisq']
        line += [el_cr[0], el_crchi2]
        
        
        noise, noisechi2 = fit( xbins, histAC, logGauss, sel=[ xbins>0, histAC>10 ], log=True )
        if first: header += ['acPos_std', 'acPos_chisq']
        line += [noise[0], noisechi2]
        if plot:
            ax.step( xbins, hist, where='mid', label = 'all' )
            ax.step( xbins, histOS, where='mid', label = 'os' )
            ax.step( xbins, histAC, where='mid', label = 'ac' )
            gOS = gaussian(xbins,*el)
            axdiff.step( xbins[gOS>countCut], (histOS[gOS>countCut] - gOS[gOS>countCut])/gOS[gOS>countCut], where='mid', label = 'os(%.2f)'%(el[0]), color = 'C1' )
            gAC = gaussian(xbins,*el_cr)
            axdiff.step( xbins[gAC>countCut], (histAC[gAC>countCut] - gAC[gAC>countCut])/gAC[gAC>countCut], where='mid', label = 'ac<(%.2f)'%(el_cr[0]), color = 'C2' )
            gAC = gaussian(xbins,*noise)
            axdiff.step( xbins[gAC>countCut], (histAC[gAC>countCut] - gAC[gAC>countCut])/gAC[gAC>countCut], where='mid', label = 'ac>(%.2f)'%(noise[0]), color = 'C3' )
        #print el[0], el_cr[0], noise[0], np.sqrt(noise[0]**2 - el_cr[0]**2)

        result += [ line ]

        if fft:
            fftAC = np.fft.fftshift( np.abs( np.fft.fft( histAC ) ) )/np.sqrt( len(histAC) )
            axfft.step( xbins, fftAC, where='mid', label = 'fftac' )
            fftOS = np.fft.fftshift( np.abs( np.fft.fft( histOS ) ) )/np.sqrt( len(histOS) )
            axfft.step( xbins, fftOS, where='mid', label = 'fftos' )
        first = False
    if plot:
        axdiff.set_xlabel('pixel energy [adu]')
        ax.set_ylabel('pixel count')
        ax.set_yscale('log')
        if fft: axfft.set_yscale('log')
        ax.set_xlim([-100,100])
        ax.grid(True)
        axdiff.grid(True)
        ax.legend()
        if fft: axfft.legend(loc=2)
        axdiff.legend()
        plt.subplots_adjust(hspace=.1, wspace=.1)
        print 'plot generated in PNG file', 'ccd%s.png'%runID
        fig.savefig('ccd%s.png'%runID)
    
    print 'table generated in CSV file', 'ccd%s.csv'%runID
    np.savetxt('ccd%s.csv'%runID, result, header=', '.join(header), fmt='%s',delimiter=', ')

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
