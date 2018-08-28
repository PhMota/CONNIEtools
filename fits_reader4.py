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

def removeHitsFromrunID( runID, outputfolder = None ):
    FITSfile = listFITS_runID( runID )[0]
    outputfilenohits = 'hsub_'+os.path.basename(FITSfile)
    outputfilehits = 'hits_'+os.path.basename(FITSfile)
    
    if not outputfolder is None:
        if not os.path.exists(outputfolder):
            os.makedirs( outputfolder )
    else:
        outputfolder = os.path.dirname(FITSfile)
    
    hdulisthits = astropy.io.fits.open(FITSfile)
    hdulistnohits = astropy.io.fits.open(FITSfile)
    print FITSfile
    ROOTfile = getAssociatedCatalog( FITSfile )[0]
    print ROOTfile
    hitscatalog = root_numpy.root2array( ROOTfile, treename = 'hitSumm', branches = ['xPix','yPix','ePix','ohdu'], selection='runID==%s'%runID )
    #print hitscatalog.shape
    for index in range(len(hdulisthits)):
        ohdu = hdulisthits[index].header['OHDU']
        #print 'hdulist.data.shape', hdulisthits[index].data.shape
        hitsohdu = hitscatalog[ hitscatalog['ohdu']==ohdu ]
        #print ohdu
        #print 'hitsohdu.shape', hitsohdu.shape
        hdulisthits[index].data[:,:] = -1e9
        #print 'hits', hdulisthits[index].data[ hdulisthits[index].data != -1e9 ].shape
        hits_xy = np.array( [ [iy,ix] for x,y in zip( hitsohdu['xPix'], hitsohdu['yPix']) for ix,iy in zip(x,y) ] )
        #print 'max_xy', np.max( hits_xy[:,0] ), np.max( hits_xy[:,1] )
        #print 'hits_xy.shape', hits_xy.shape
        hits_e = np.array( [ ie for e in hitsohdu['ePix'] for ie in e ] )
        #print 'hits_e.shape', hits_e.shape
        hdulistnohits[index].data[hits_xy[:,0],hits_xy[:,1]] = -1e9
        #print 'nohits', hdulistnohits[index].data[ hdulistnohits[index].data != -1e9 ].shape
        hdulisthits[index].data[hits_xy[:,0],hits_xy[:,1]] = hits_e
        #print 'hits', hdulisthits[index].data[ hdulisthits[index].data != -1e9 ].shape
    
    hdulisthits.writeto( outputfolder + '/' + outputfilehits )
    hdulistnohits.writeto( outputfolder + '/' + outputfilenohits )

def plotrunID( ohdu, *FITSfiles ):
    FITSfiles = FITSfiles[0]
    runID = getrunIDFromPath( FITSfiles[0] )
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for FITSfile in FITSfiles:
        hdulist = astropy.io.fits.open( FITSfile )
        data = [ hdu.data for hdu in hdulist if hdu.header['OHDU'] == ohdu ][0]
        data = data[data!=-1e9]
        data = data[data!=1e10]
        dx = 1
        bins = np.r_[min(data.flatten()):max(data.flatten()):dx]
        xbins = (bins[1:]+bins[:-1])*.5
        hist, tmp = np.histogram( data.flatten(), bins = bins )
        ax.step( xbins, hist, where='mid' )
    ax.set_yscale('log')
    ax.set_xlim([-100,100])
    print 'ccd%s.png'%runID
    fig.savefig('ccd%s.png'%runID)
    

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
    runID = None
    outfolder = None
    ohdu = None
    if '--outfolder' in sys.argv:
        outfolder = sys.argv[sys.argv.index('--outfolder')+1]
        print 'outfolder', outfolder
    if '--ohdu' in sys.argv:
        ohdu = sys.argv[sys.argv.index('--ohdu')+1]
        print 'ohdu', ohdu
    if '--runID' in sys.argv:
        runID = sys.argv[sys.argv.index('--runID')+1]
        print 'runID', runID
    if '--removehits' in sys.argv:
        if not runID is None and not outfolder is None:
            removeHitsFromrunID( int(runID), outputfolder=outfolder )
    if '--plot' in sys.argv:
        if not ohdu is None:
            plotrunID( int(ohdu), sys.argv[sys.argv.index('--plot')+1:] )
