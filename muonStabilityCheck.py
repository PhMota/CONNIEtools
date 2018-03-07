import numpy as np
import sys
import os
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange, epoch2num
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm, chisquare

from scipy.optimize import minimize
from scipy.misc import factorial

import root_numpy

#plt.rc('text', usetex=True)

CuPeak = 8.048
CuPeak2 = 8.905
SiPeak = 1.7
#logEloss = log( E0/sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2 + 675**2  ) )

'atan2(yMax-yMin, xMax-xMin)'
'atan2(yBary0-yMin, xBary0-xMin)'
'atan2( (xBary0-xMin)+(yBary0-yMin)**2,  )'
'abs( atan2(yMax-yMin, xMax-xMin)-atan2(yBary0-yMin, xBary0-xMin) ) < .05'

"""
129M    /share/storage1/connie/data_analysis/processed_data/2016/scienceTest/15Aug2016/00/data_29_to_50/ext/catalog/catalog_29_to_50.skim1.root
1011M /share/storage1/connie/data_analysis/processed_data/runs/003/data_92_to_156/ext/catalog/catalog_92_to_156.skim1.root
392M    /share/storage1/connie/data_analysis/processed_data/runs/009/data_238_to_302/ext/catalog/catalog_238_to_302.skim1.root
237M    /share/storage1/connie/data_analysis/processed_data/runs/010/data_303_to_304/ext/catalog/catalog_303_to_304.skim1.root
5.3G    /share/storage1/connie/data_analysis/processed_data/runs/011/data_305_to_324/ext/catalog/catalog_305_to_324.skim1.root
2.4G    /share/storage1/connie/data_analysis/processed_data/runs/012/data_325_to_332/ext/catalog/catalog_325_to_332.skim1.root
2.7G    /share/storage1/connie/data_analysis/processed_data/runs/013/data_333_to_368/ext/catalog/catalog_333_to_368.skim1.root
196M    /share/storage1/connie/data_analysis/processed_data/runs/014/data_369_to_370/ext/catalog/catalog_369_to_370.skim1.root
396M    /share/storage1/connie/data_analysis/processed_data/runs/015/data_371_to_421/ext/catalog/catalog_371_to_421.skim1.root
271M    /share/storage1/connie/data_analysis/processed_data/runs/016/data_422_to_453/ext/catalog/catalog_422_to_453.skim1.root
1.7G    /share/storage1/connie/data_analysis/processed_data/runs/017/data_454_to_461/ext/catalog/catalog_454_to_461.skim1.root
394M    /share/storage1/connie/data_analysis/processed_data/runs/018/data_462_to_505/ext/catalog/catalog_462_to_505.skim1.root
610M    /share/storage1/connie/data_analysis/processed_data/runs/018/data_506_to_564/ext/catalog/catalog_506_to_564.skim1.root
543M    /share/storage1/connie/data_analysis/processed_data/runs/018/data_565_to_624/ext/catalog/catalog_565_to_624.skim1.root
512M    /share/storage1/connie/data_analysis/processed_data/runs/018/data_625_to_684/ext/catalog/catalog_625_to_684.skim1.root
1.3G    /share/storage1/connie/data_analysis/processed_data/runs/018/data_685_to_744/ext/catalog/catalog_685_to_744.skim1.root
507M    /share/storage1/connie/data_analysis/processed_data/runs/018/data_745_to_804/ext/catalog/catalog_745_to_804.skim1.root
632M    /share/storage1/connie/data_analysis/processed_data/runs/018/data_805_to_860/ext/catalog/catalog_805_to_860.skim1.root
1.5G    /share/storage1/connie/data_analysis/processed_data/runs/019/data_861_to_877/ext/catalog/catalog_861_to_877.skim1.root
929M    /share/storage1/connie/data_analysis/processed_data/runs/021/data_898_to_903/ext/catalog/catalog_898_to_903.skim1.root
2.6G    /share/storage1/connie/data_analysis/processed_data/runs/022/data_2001_to_2016/ext/catalog/catalog_2001_to_2016.skim1.root
903M    /share/storage1/connie/data_analysis/processed_data/runs/023/data_2017_to_2031/ext/catalog/catalog_2017_to_2031.skim1.root
733M    /share/storage1/connie/data_analysis/processed_data/runs/024/data_2032_to_2036/ext/catalog/catalog_2032_to_2036.skim1.root
3.1G    /share/storage1/connie/data_analysis/processed_data/runs/025/data_2037_to_2097/ext/catalog/catalog_2037_to_2097.skim1.root
3.0G    /share/storage1/connie/data_analysis/processed_data/runs/025/data_2098_to_2157/ext/catalog/catalog_2098_to_2157.skim1.root
3.3G    /share/storage1/connie/data_analysis/processed_data/runs/025/data_2158_to_2217/ext/catalog/catalog_2158_to_2217.skim1.root
3.3G    /share/storage1/connie/data_analysis/processed_data/runs/025/data_2218_to_2277/ext/catalog/catalog_2218_to_2277.skim1.root
3.3G    /share/storage1/connie/data_analysis/processed_data/runs/025/data_2278_to_2337/ext/catalog/catalog_2278_to_2337.skim1.root
2.9G    /share/storage1/connie/data_analysis/processed_data/runs/025/data_2338_to_2388/ext/catalog/catalog_2338_to_2388.skim1.root
891M    /share/storage1/connie/data_analysis/processed_data/runs/026/data_2389_to_2453/ext/catalog/catalog_2389_to_2453.skim1.root
674M    /share/storage1/connie/data_analysis/processed_data/runs/026/data_2454_to_2518/ext/catalog/catalog_2454_to_2518.skim1.root
"""

"""
/share/storage2/connie/data_analysis/processed02_data/runs/
/share/storage2/connie/data_analysis/processed02_data/runs/*/*/ext/catalog/*.skim1.root
"""

"""
5506708 3,1G -rw-r--r--  1 mota mota 3,1G Out 23 18:39 stat_gain_catalog_2037_to_2097.skim1.root
5506161 3,0G -rw-r--r--  1 mota mota 3,0G Out 23 18:47 stat_gain_catalog_2098_to_2157.skim1.root
5506190 3,3G -rw-r--r--  1 mota mota 3,3G Out 23 18:52 stat_gain_catalog_2158_to_2217.skim1.root
5506211 3,3G -rw-r--r--  1 mota mota 3,3G Out 23 18:57 stat_gain_catalog_2218_to_2277.skim1.root
5506064 3,3G -rw-r--r--  1 mota mota 3,3G Out 23 19:03 stat_gain_catalog_2278_to_2337.skim1.root
5506661 2,9G -rw-r--r--  1 mota mota 2,9G Out 23 19:08 stat_gain_catalog_2338_to_2388.skim1.root
5506577 889M -rw-r--r--  1 mota mota 889M Out 23 19:09 stat_gain_catalog_2389_to_2453.skim1.root
5506684 391M -rw-r--r--  1 mota mota 391M Out 23 19:10 stat_gain_catalog_238_to_302.skim1.root
5506693 672M -rw-r--r--  1 mota mota 672M Out 23 19:11 stat_gain_catalog_2454_to_2518.skim1.root
5506703 129M -rw-r--r--  1 mota mota 129M Out 23 19:11 stat_gain_catalog_29_to_50.skim1.root
5506662 396M -rw-r--r--  1 mota mota 396M Out 23 19:12 stat_gain_catalog_371_to_421.skim1.root
5506665 395M -rw-r--r--  1 mota mota 395M Out 23 19:12 stat_gain_catalog_462_to_505.skim1.root
5507270 612M -rw-r--r--  1 mota mota 612M Out 23 19:13 stat_gain_catalog_506_to_564.skim1.root
5507066 542M -rw-r--r--  1 mota mota 542M Out 23 19:14 stat_gain_catalog_565_to_624.skim1.root
5507610 514M -rw-r--r--  1 mota mota 514M Out 23 19:15 stat_gain_catalog_625_to_684.skim1.root
5507632 1,3G -rw-r--r--  1 mota mota 1,3G Out 23 19:17 stat_gain_catalog_685_to_744.skim1.root
5507634 507M -rw-r--r--  1 mota mota 507M Out 23 19:18 stat_gain_catalog_745_to_804.skim1.root
5507669 632M -rw-r--r--  1 mota mota 632M Out 23 19:19 stat_gain_catalog_805_to_860.skim1.root
"""
cats = [ [2037,2097], [2098,2157], [2158,2217], [2218,2277], [2278,2337], [2338,2388], [2389,2453], [238,302], [2454,2518], [29,50], [371,421], [462,505], [506,564], [565,624], [625,684], [685,744], [745,804], [805,860] ]

cats2 = [ [0,     29,      50],
[9,     246,     302],
[15,     371,     388],
[15,     392,     421],
[18,     476,     542],
[18,     550,     572],
[18,     584,     660],
[18,     668,     694],
[18,     704,     718],
[18,     730,     860],
[25,     2060,    2106],
[25,     2118,    2171],
[25,     2174,    2235],
[25,     2240,    2301],
[25,     2305,    2388],
[26,     2391,    2518] ]

cats3 = [ 
['000', 1,    29,      50],
['009', 0,    238,     302],
['015', 0,    371,     421],
['018a', 0,    462,     505],
['018b', 1,    506,     564],
['018c', 1,    565,     624],
['018d', 1,    625,     684],
['018e', 1,    685,     744],
['018f', 1,    745,     804],
['018g', 1,    805,     860],
['025a', 1,    2037,    2097],
['025b', 1,    2098,    2157],
['025c', 1,    2158,    2217],
['025d', 1,    2218,    2277],
['025e', 1,    2278,    2337],
['025f', 1,    2338,    2388],
['026a', 1,    2389,    2453],
['026b', 1,    2454,    2518],
]


def summary( var_str ):
    print [ '%s: '%k + str(v) for k,v in globals().iteritems() if k == var_str ][0]

selectionPhilipe = r'flag==0 && abs(log10( E0/sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2 + 675**2  ) )-2.7055039797010396)<2*sqrt(7.326161363760477-2.7055039797010396**2)'

selectionJuan = r'E1>0 && flag==0 && sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))>20 && abs(n0/sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))-4.4)<1'

selectionJuanAngle = r'E1>0 && flag==0 && sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))>20 && abs(n0/sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))-4.4)<1 '

selectionPhilipeAngle = r'flag==0 && abs(log10( E0/sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2 + 675**2  ) )-2.7055039797010396)<2*sqrt(7.326161363760477-2.7055039797010396**2) && atan2( sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2), 675 ) > 0.785398'

#selectionDC = r'&& dcFlag==1'
selectionDC = r''

selection45 = r'&& atan2( sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2), 675 ) > 0.785398'

selections = {
    'Juan': selectionJuan + selectionDC,
    'Philipe': selectionPhilipe + selectionDC,
    'Juan45': selectionJuan + selectionDC + selection45,
    'Philipe45': selectionPhilipe + selectionDC + selection45,
    }

def np_and(*criteria):
    return np.all( [ criterium for criterium in criteria ], axis=0)

def maskOutliers( data, radius = 1 ):
    return np.abs( data - np.median(data) ) < radius*np.std(data)

def removeoutliers( data, radius = 1 ):
    return data[ maskOutliers(data, radius=radius) ]

def getValuesFromCatalog( catalog, treename = 'hitSumm', branches = [], selection = '' ):
    if branches is []: return None
    data = root_numpy.root2array( catalog, treename = treename, branches = branches, selection = selection ).view(np.recarray)
    return np.array([ data[branch] for branch in branches ])

def applySelection( files, selection, branches, output = None, fmt = '%.18e' ):
    variables = {branch: [] for branch in branches}
    for file_ in files:
        base = os.path.basename(file_)
        outfile = 'plots/%s.%s.dat'%(base,output)
        print 'file', file_
        print 'extracting branches', branches
        print "branches available", root_numpy.list_branches(file_, 'hitSumm')
        print "missing branches", [ branch for branch in branches if not branch in root_numpy.list_branches(file_, 'hitSumm') ]
        print "branches available", root_numpy.list_branches(file_, 'config')
        print
        data = getValuesFromCatalog( file_, branches = branches, selection = selection )
        print 'read', len(data[0]), 'entries'
        header = " ".join( branches )
        #data_ = [ variables[branch] for branch in branches ]
        print 'writing into', outfile
        print 'length', len(branches), len(fmt.split())
        np.savetxt( outfile, zip(*data), header=header, fmt=fmt )        
        #for branch, datum in zip(branches, data):
            #variables[branch] = np.append( variables[branch], datum )
    #header = " ".join( branches )
    #data_ = [ variables[branch] for branch in branches ]
    #print 'writing into', 'plots/%s.dat'%output 
    #print 'length', len(branches), len(fmt.split())
    #np.savetxt( 'plots/%s.dat'%output, zip(*data_), header=header, fmt=fmt )
    return

def plotTime( data = None, fnames = None, branch = None, output = None, ohdu = None, labels = None, func = None, ylabel = None, perRunID = True, perCatalog = False ):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel( 'date' )
    ax.set_ylabel( branch )

    runIDnames = [ 'plots/'+os.path.splitext(os.path.basename( fname ))[0]+'.runID.dat' for fname in fnames ]

    if data is None:
        data = []
        for fname in fnames:
            print 'reading file', fname
            data += [ np.genfromtxt( fname, names = True ) ]
    
    for index, datum in enumerate(data):
        mask = np.ones_like( datum[branch] )
        if not ohdu is None:
            mask = np_and( mask, datum['ohdu'] == ohdu )

        if perRunID:
            runIDs = np.array(list(set( datum['runID'] )))
            x = np.array( [ datum[ 'expoStart' ][ np_and(mask, datum['runID'] == runID) ][0] for runID in runIDs ] )
            y = np.array( [ func( datum[ branch ][ np_and(mask, datum['runID'] == runID) ] ) for runID in runIDs ] )
            ax.set_ylabel( ylabel )
            
            header = ' '.join(['runID', 'expoStart', ylabel] )
            data_ = [ runIDs, x, y ]
            print os.splitext(os.basename( fnames[index] ))[0]+'.runID.dat'
            np.savetxt( runIDnames[index] )
            
        else:
            x = epoch2num( datum[ 'expoStart' ][mask] )
            y = data[ branch ][mask]
            if not func is None:
                y = func(y)

        ax.plot_date( epoch2num(x), y, '.', label = labels[index] )
        
        if perCatalog:
            runIDs_ = [ runIDs[ np_and( runIDs>=start, runIDs<end ) ] for n,start,end in cats2 ]
            x_ = [ np.mean( x[ np_and( runIDs>=start, runIDs<end ) ] ) for n,start,end in cats2 ]
            y_ = [ y[ np_and( runIDs>=start, runIDs<end ) ] for n,start,end in cats2 ]
            median = np.median( y_, axis=1 )
            std = np.std( y_, axis=1 )
            print len(median), len(cats2)
    
    ax.xaxis.set_major_formatter( DateFormatter('%d %b %y') )
    fig.autofmt_xdate()
    fig.savefig( output, format='pdf')
    
    return

def plotDistribution( data = None, fnames = None, branch = None, output = None, bins = None ):
    if data is None:
        data = np.genfromtxt( fname, names = True )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel( branch )
    ax.set_ylabel( 'dN/d%s'%branch )
    
    y = data[ branch ]
    
    ax.hist( y, bins = bins, histtype = 'step' )
    
    fig.savefig( output, format='pdf')
    return
    
    
def applyToAll( catalogs, selection = '' ):
    nJ = {'all':[]}
    nJa = {'all':[]}
    nPh = {'all':[]}
    nPha = {'all':[]}
    runID_ = []
    expoStart_ = []
    print 'catalogs', catalogs
    #ohdus = [2,3,4,5,6,8,9]
    ohdus = [2,3,4,5,6,7,8,9,10]
    for catalog in catalogs:
        print catalog
        runIDs, ohdus_, expoStarts = map( set, getValuesFromCatalog( catalog, branches = ['runID','ohdu','expoStart'] ))
        #runIDs = sorted(runIDs)
        ohdus_ = sorted(ohdus_)
        print catalog, os.path.splitext(catalog)
        print runIDs
        print ohdus_
        #print expoStarts
        print 'reading'
        E0J, ohduJ, expoStartJ, runIDJ = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionJuan + selectionDC )
        E0Ja, ohduJa, expoStartJa, runIDJa = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionJuanAngle + selectionDC )
        E0Ph, ohduPh, expoStartPh, runIDPh = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionPhilipe + selectionDC )
        E0Pha, ohduPha, expoStartPha, runIDPha = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionPhilipeAngle + selectionDC )
        print 'done'
        #print len(E0J)
        
        for runID in runIDs:
            print 'expoStart', expoStartJ[runIDJ==runID]
            if len(expoStartJ[runIDJ==runID]) == 0:
                continue
            runID_ += [runID]
            #E0J,ohduJ,expoStart = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart'], selection = 'runID == %d &&' %(runID) + selectionJuan )
            #E0Ph,ohduPh = getValuesFromCatalog( catalog, branches = ['E0','ohdu'], selection = 'runID == %d &&' %(runID) + selectionPhilipe )
            #E0Pha,ohduPha = getValuesFromCatalog( catalog, branches = ['E0','ohdu'], selection = 'runID == %d &&' %(runID) + selectionPhilipeAngle )
            
            print runID, expoStartJ[runIDJ==runID][0] #, expoStartJa[runIDJa==runID][0], expoStartPh[runIDPh==runID][0], expoStartPha[runIDPha==runID][0]
            #print runID, E0J[runIDJ==runID].shape[0], E0Ja[runIDJa==runID].shape, E0Ph[runIDPh==runID].shape, E0Pha[runIDPha==runID].shape, E0J[ np.all([runIDJ==runID, ohduJ==2], axis = 0) ].shape
            print ohduJ[ runIDJ == runID ]
            
            expoStart_ += [ expoStartJ[runIDJ==runID][0] ]
            
            nJ['all'] += [ E0J[runIDJ==runID].shape[0] ]
            nJa['all'] += [ E0Ja[runIDJa==runID].shape[0] ]
            nPh['all'] += [ E0Ph[runIDPh==runID].shape[0] ]
            nPha['all'] += [ E0Pha[runIDPha==runID].shape[0] ]
            
            for ohdu in ohdus:
                mask = np.all( [runIDJ==runID, ohduJ==ohdu], axis=0 )
                try: nJ[ohdu] += [ ohduJ[mask].shape[0] ]
                except: nJ[ohdu] = [ ohduJ[mask].shape[0] ]
                
                mask = np.all( [runIDJa==runID, ohduJa==ohdu], axis=0 )
                try: nJa[ohdu] += [ ohduJa[mask].shape[0] ]
                except: nJa[ohdu] = [ ohduJa[mask].shape[0] ]
                
                mask = np.all( [runIDPh==runID, ohduPh==ohdu], axis=0 )
                try: nPh[ohdu] += [ ohduPh[mask].shape[0] ]
                except: nPh[ohdu] = [ ohduPh[mask].shape[0] ]
                
                mask = np.all( [runIDPha==runID, ohduPha==ohdu], axis=0 )
                try: nPha[ohdu] += [ ohduPha[mask].shape[0] ]
                except: nPha[ohdu] = [ ohduPha[mask].shape[0] ]
                
            data = [runID_] + [expoStart_] + sum([ [nJ[ohdu], nJa[ohdu], nPh[ohdu], nPha[ohdu]] for ohdu in ohdus + ['all'] ], [])
            header = " ".join(["runID"] + ["expoStart"]  + sum([ [ "nJ[%s]"%ohdu, "nJa[%s]"%ohdu, "nPh[%s]"%ohdu, "nPha[%s]"%ohdu] for ohdu in ohdus + ['all'] ], []))
            np.savetxt( 'plots/muonCount.dat', zip(*data), header = header, fmt='%d' )
            

def makePlot( fname ):
    data = np.genfromtxt( fname, names = True )
    catData = [ data[ np.all( [ cat[-2]<data['runID'], data['runID']<=cat[-1] ], axis=0) ] for cat in cats2 ]
    
    pp = PdfPages('plots/muonCount.all.pdf')
    
    labels = {'nJa':r'Juan(theta>45)', 'nPha':r'Philipe2(theta>45)'}
    count = 0

    chi = {}
    mean = {}
    std = {}
    
    ohdus = [2,3,4,5,6,7,8,9,10]
    for ohdu in ohdus:
        fig = plt.figure()
        fig2 = plt.figure()
        
        size = fig.get_size_inches()
        size = [ size[0], size[1]*(1+(len(data.dtype.names)-1)/3) ]
        fig.set_size_inches(size)
        print 'ohdu', ohdu
        
        ax = fig.add_subplot( 1+(len(data.dtype.names)-1)/3, 1 , 1+count)
        ax.set_ylabel(r"Nmu")
        ax.set_xlabel(r"time")
        ax.set_title(r"ohdu=%s"%ohdu)
        ax.grid(True)
        
        #bx2 = bx.twinx()

        count += 1

        nbins = 50
        y_all = {}
        xx2_all = {}
        y_mean = {}
        usedLabels = []
        for index,edges in enumerate(cats2):
            print 'catalog', edges
            for j, key in enumerate( ['nJa','nPha'] ):
                if not key in y_all: 
                    y_all[key] = []
                    xx2_all[key] = []
                    chi[key] = {'all':[], 18:[], 25:[]}
                    y_mean[key] = []

                y = catData[index]["%s%s"%(key,ohdu)]
                xx = catData[index]['expoStart']
                xx2 = catData[index]['runID']
                x = epoch2num(xx)
                y_all[key] = np.append( y_all[key], y )
                xx2_all[key] = np.append( xx2_all[key], xx2 )

                ax.plot_date( x, y, '.', alpha=.1 )

                std_ = np.std(y)
                median_ = np.median(y)
                
                y_mean[key] = np.append( y_mean[key], median_ )

                removeoutliers = abs(y - median_) < 3*std_
                meanR, stdR = np.mean(y[removeoutliers]), np.std(y[removeoutliers])
                removedRunIDs = xx2[ abs( y - meanR) >=3 * stdR ]
                print 'removedRunIDs', removedRunIDs, len(removedRunIDs)*1./len(xx2)
                color = ['k','r'][j]
                #ax.errorbar( epoch2num( .5*(min(xx)+max(xx)) ), np.median(y), yerr = std_ , color=color, fmt='o', alpha=.2 )
                ax.errorbar( epoch2num( .5*(min(xx)+max(xx)) ), meanR, yerr = stdR, color=color, fmt='o', label = labels[key] if not labels[key] in usedLabels else '' )
                usedLabels += [labels[key]]

        ax.legend( fancybox=True, framealpha=0.1)
        dateFormatter = DateFormatter('%d %b %y')
        ax.xaxis.set_major_formatter( dateFormatter )
        fig.autofmt_xdate()
        fig.savefig(pp, bbox_inches='tight', format='pdf')
        
        bx1 = fig2.add_subplot(121)
        bx2 = fig2.add_subplot(122)
        bx1.set_xlabel(r'Nmu')
        bx1.set_ylabel(r'count')
        bx1.set_title(r'distribution ohdu=%d'%ohdu)

        bx2.set_xlabel(r'Nmu')
        bx2.set_ylabel(r'count')
        bx2.set_title(r'distribution ohdu=%d'%ohdu)
        
        for index, key in enumerate(y_all.keys()):
            
            bx = [ bx1, bx2 ][index]
            mean, std = norm.fit( y_all[key] )
            mask = maskOutliers( y_all[key], 3 )
            meanR, stdR = norm.fit( y_all[key][mask] )
            
            count1, edges1, _ = bx.hist( y_all[key], bins = np.arange(min(y_all[key]), max(y_all[key]), 5)-.5 , histtype='step', label = labels[key] )
            h = count1
            
            x = .5*(edges1[1:] + edges1[:-1])
            dx = edges1[1]-edges1[0]
            #print 'dx', dx, len(x)
            Norm_ = np.sum(h)*dx
            print 'raw mean', mean, std, 'removedoutliers', meanR, stdR
            removedRunIDs = xx2_all[key][ abs( y_all[key] - meanR ) >= 3*stdR ]
            print 'removedRunIDs', removedRunIDs, len(removedRunIDs)*1./len(xx2_all[key])
            
            maskbin = np.all( [ x > meanR-3*stdR, x < meanR+3*stdR ], axis=0 )
            
            
            x = x[maskbin]
            h = h[maskbin]
            chi[key]['all'] += [ 1./len(x)*chisquare(h, f_exp=Norm_*norm.pdf(x, meanR, stdR) )[0] ]
            mean[key]['all'] += [ meanR ]
            std[key]['all'] += [ stdR ]
            
            bx.plot( x, Norm_*norm.pdf(x, meanR, stdR), '-', label = 'gaussian chi2r=%.2f'%( 1./len(x)*chisquare(h, f_exp=Norm_*norm.pdf(x, meanR, stdR) )[0] ) )
            bx.plot( x, Norm_*norm.pdf(x, meanR, np.sqrt(meanR)), '-', label = 'gpoisson chi2r=%.2f'%( 1./len(x)*chisquare(h, f_exp=Norm_*norm.pdf(x, meanR, np.sqrt(meanR)) )[0] ) )
            
            bx.legend(  fancybox=True, framealpha=0.1, title=r'mean=%.1f, '%meanR + r'std=%.1f'%stdR )
            
        fig2.savefig(pp, bbox_inches='tight', format='pdf')
        
        #for index, key in enumerate(y_all.keys()):
            #for run in [18,25]:
                #runIDs__ = list(set( sum( [ list(catData[index]['runID']) for index, _ in enumerate(cats2) if cats2[index][0] == run ], [] ) ))
                #count__ = [  for runid in runIDs__ ]
                #h, edges, _ = ax.hist( count__, bins = np.arange(0, np.max(count__), 1)-.5, histtype='step', label = '3--7keV count' )
                #x = .5*(edges[1:] + edges[:-1])
                #dx = edges1[1]-edges1[0]
                #Norm_ = np.sum(h)*dx
                
                #chi[run] += [ 1./len(x)*chisquare(h, f_exp=np.sum(h)*poisson( x, mean ) )[0] ]
        

    pp.close()


def poisson(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    #return (lamb**k/factorial(k)) * np.exp(-lamb)
    return np.exp( k*np.log(lamb) - lamb - np.log(factorial(k)) )

def poisson2(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)
    #return np.exp( k*np.log(lamb) - lamb - np.log(factorial(k)) )


def negLogLikelihood(params, data):
    """ the negative log-Likelohood-Function"""
    lnl = - np.sum(np.log(poisson( data, params[0] )))
    return lnl

def compute_chisquare(f1, f2):
    result = np.sum( ( f1-f2)**2/f2 )
    return result, result/len(f1)


def plotTime2(output, ohdu = 2, ylabel = '', data = [], labels = [] ):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('ohdu=%d' % ohdu)
    ax.set_xlabel('runID')
    ax.set_ylabel( ylabel )
    usedlabel = []
    for index, datum in enumerate(data):
        thisohdu = datum['ohdu'] == ohdu
        catalogs = list(set( datum['catalog'] ) )
        for catalog in catalogs:
            thiscat_ohdu = np.all( [ datum['catalog'] == catalog, thisohdu ], axis=0 )
            xx2 = datum['expoStart'][ thiscat_ohdu ]
            x2 = epoch2num(xx2)
            y2 = datum['N'][thiscat_ohdu]
            if len(x2) > 1:
                ax.plot_date( x2, y2, '.', alpha = .1 )
                ax.errorbar(
                    epoch2num( .5*(min(xx2)+max(xx2)) ), 
                    np.median(y2), 
                    yerr = np.std(y2),
                    color=['k','b'][index],
                    fmt='o',
                    label = labels[index] if not labels[index] in usedlabel else '' )
                usedlabel += [ labels[index] ]
        #ax.errorbar( epoch2num( .5*(min(xx)+max(xx)) ), np.median(y[removeoutliers]), yerr = np.std(y[removeoutliers]) , color=color, fmt='o', label = labels[key] if not labels[key] in usedLabels else '' )
        
    ax.legend( fancybox=True, framealpha=0.1)
    dateFormatter = DateFormatter('%d %b %y')
    ax.xaxis.set_major_formatter( dateFormatter )
    fig.autofmt_xdate()
    fig.savefig(output, format='pdf')
    return

def plotDistributionOutliers( output, ohdu = 2, ylabel = '', data = None, outlierradius = 3, masks = [], labels = [], binsize = 1, table = {} ):
    print 'ohdu', ohdu
    fig = plt.figure()
    thisohdu = data['ohdu'] == ohdu
    
    for index, mask_ in enumerate(masks):
        
        thisrun_ohdu = np.all( [ mask_, thisohdu ], axis = 0)
        print 'number', sum(thisrun_ohdu)
        
        ax = fig.add_subplot(1, len(masks), 1+index )
        ax.set_title('ohdu=%d'%ohdu)
        ax.set_xlabel( ylabel )
        ax.set_ylabel('count')
        
        count = data['N'][thisrun_ohdu]
        
        mean, std = norm.fit( count )
        mask = maskOutliers( count, outlierradius )
        meanR, stdR = norm.fit( count[mask] )
        
        print 'mean, std', mean, std
        h, edges, _ = ax.hist( count, bins = np.arange(np.min(count), np.max(count), binsize)-.5, histtype='step', label = labels[index] )

        x = .5*(edges[1:] + edges[:-1])
        dx = edges[1]-edges[0]
        Norm_ = np.sum(h)*dx
        print 'raw mean', mean, std, 'removedoutliers', meanR, stdR
        runIDs = data['runID'][ thisrun_ohdu ]
        removedRunIDs = np.array(runIDs)[ abs( count - meanR ) >= outlierradius * stdR ]
        print 'removedRunIDs', removedRunIDs, len(removedRunIDs)*1./len(runIDs)
        
        maskbin = np.all( [ x > meanR-3*stdR, x < meanR+3*stdR ], axis=0 )
        
        x = x[maskbin]
        h = h[maskbin]
        
        print 'poisson for meanR', poisson( meanR, meanR ), norm.pdf(meanR, meanR, np.sqrt(meanR))
        chi2rP = 0
        
        #if np.isnan( poisson( meanR, meanR ) ):
        if meanR > 100:
            chi2rP = ( 1./len(x)*chisquare(h, f_exp=Norm_*norm.pdf(x, meanR, np.sqrt(meanR) ) )[0] )
            ax.plot( x, Norm_*norm.pdf(x, meanR, np.sqrt(meanR)), '-', label = 'gpoisson chi2r=%.2f'%chi2rP )
        else:
            chi2rP = ( 1./len(x)*chisquare(h, f_exp=Norm_*poisson( x, meanR ) )[0] )
            ax.plot( x, Norm_*poisson( x, meanR ), '-', label = 'poisson chi2r=%.2f'%chi2rP )
        chi2rG = ( 1./len(x)*chisquare(h, f_exp=Norm_*norm.pdf(x, meanR, stdR) )[0] )
        ax.plot( x, Norm_*norm.pdf(x, meanR, stdR), '-', label = 'gaussian chi2r=%.2f'%chi2rG )
        
        ax.legend( fancybox=True, framealpha=0.1, title=r'N=%d, '%(len(count)) + 'mean=%.1f, '%(meanR) + r'std=%.1f'%stdR )
        
        label = labels[index]
        if not label in table: table[label] = {}
        if not ohdu in table[label]: table[label][ohdu] = {}
        table[label][ohdu]['N' ] = len(count)
        table[label][ohdu]['mean' ] = meanR
        table[label][ohdu]['std' ] = stdR
        table[label][ohdu]['chi2rP' ] = chi2rP
        table[label][ohdu]['chi2rG' ] = chi2rG
        
    fig.savefig( output, format='pdf')
    return table

def computeCounts( data, output = '', ohdus = range(15), mask = None ):
    
    columns = [ 'catalog', 'onflag', 'runID', 'expoStart', 'ohdu', 'N' ]
    #fmt = [ '%s', '%d', '%d', '%d', '%d' ]
    table = { column: [] for column in columns }
    for thiscat in cats3:
        print thiscat
        for thisrunID in range( thiscat[-2], thiscat[-1]+1 ):
            thisdata = data[ data['runID'] == thisrunID ]
            if len(thisdata) == 0:
                continue
            for thisohdu in ohdus:
                table['catalog'] += [ thiscat[0] ]
                table['onflag'] += [ thiscat[1] ]
                table['runID'] += [ int(thisrunID) ]
                table['expoStart'] += [ long(thisdata['expoStart'][0]) ]
                table['ohdu'] += [ int(thisohdu) ]
                table['N'] += [ len(thisdata[ 'E1' ][ np.all( [ mask(thisdata), thisdata['ohdu'] == thisohdu ], axis=0 ) ] ) ]
                    #np.all( [ thisdata['E1'] < thisdata['gainCu']*Emax, thisdata['E1'] >= thisdata['gainCu']*Emin, thisdata['ohdu'] == thisohdu ], axis=0 ) 
                    
    print 'saving table'
    np.savetxt( 'plots/%s.dat'%output, zip(*[ table[column] for column in columns] ), header=' '.join(columns), fmt = '%s' )
    return 

def savetable( table, subs = [], output = '', ohdus = [] ):
    print 'table', table.keys()
    columns = ['N', 'mean', 'std', 'chi2rP', 'chi2rG']
    fmt = ['%d', '%.2f', '%.2f', '%.2e', '%.2e']
    fullfmt = ['%d']
    #subs = ['all', 'run025', 'run018', 'ON', 'OFF']
    data = []
    header = ['ohdu']
    for ohdu_ in ohdus:
        subdata = [ohdu_]
        for index, column in enumerate(columns):
            if ohdu_ == ohdus[0]:
                header += [column]
                fullfmt += [ fmt[index] ]
            subdata += [ table[ohdu_][column] ]
        data += [subdata]
    data = np.array(data)
    print 'shape', data.shape
    np.savetxt( output, data, header=' & '.join(header) + r'\\', fmt = ' & '.join(fullfmt) + r'\\' )
    
def printTable( tables, ohdus = [], labels = [], ylabel = '', output = '' ):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ax.set_xlabel('ohdu')
    ax.set_ylabel(ylabel)
    ax2.set_ylabel('chi2r')
    ax2.set_yscale('log')
    ax.grid(True)
    usedlabel = []
    dx = .5/len(labels)
    colors=['k', 'r', 'b', 'g']
    for ohdu in ohdus:
        for index, table in enumerate(tables):
            ax.errorbar(ohdu + index*dx, table[ohdu]['mean'], yerr = table[ohdu]['std'], fmt = 'o', color=colors[index], label = labels[index] if not labels[index] in usedlabel else '' )
            ax2.plot( ohdu + index*dx, table[ohdu]['chi2rP'], color=colors[index], marker='>', label = 'chir2P' if ohdu == ohdus else '' )
            ax2.plot( ohdu + index*dx, table[ohdu]['chi2rG'], color=colors[index], marker='<', label = 'chir2G' if ohdu == ohdus else '' )
            usedlabel += [ labels[index] ]
    ax.legend( fancybox=True, framealpha=0.1 )
    ax2.legend( fancybox=True, framealpha=0.1 )
    fig.savefig( output, format='pdf')


def plotMuon( outfile, labels = [] ):
    print 'muon'

    #outfile2 = 'plots/muonCount.pdf'
    pp2 = PdfPages( outfile + '.tmp')
    
    ohdus = [2,3,4,5,6,7,8,9,10]

    countData = {}
    table = {}
    keys = ['all', 'run018', 'run025', 'OFF', 'ON']
    
    for label in labels:
        if not os.path.exists('plots/muonCount.%s.dat'%label ):
            data = np.genfromtxt( 'plots/muon.%s.dat'%label, names = True )
            computeCounts( data, output = 'muonCount.%s'%label, ohdus = ohdus, mask = lambda x: np.ones_like(x, dtype = bool) )
        countData[label] = np.genfromtxt( 'plots/muonCount.%s.dat'%label, dtype=None, names = True )
        
        for key in keys:
            table[label][key] = {}
    binsize = 5
    
    for ohdu_ in ohdus:
        
        thisohdu = countData[label]['ohdu'] == ohdu_
        plotTime2( pp2, ylabel = 'muon selection count', data = [countData[label] for label in labels], ohdu = ohdu_, labels = labels )
        
        for label in labels:
            table[label] = plotDistributionOutliers( pp2, ylabel = 'muon %s selection'%label, data = countData[label], ohdu = ohdu_, masks = [np.ones_like(countData[label], dtype = bool)], labels = ['all'], binsize = binsize, table = table[label] )
            
            thisrun = []
            for run in ['18','25']:
                catalogs = [ cat[0] for cat in cats3 if run in cat[0] ]
                thisrun += [ np.any( [ countData[label]['catalog'] == catalog for catalog in catalogs ], axis = 0) ]

            table[label] = plotDistributionOutliers( pp2, ylabel = 'muon %s selection'%label, data = countData[label], ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['run018', 'run025'], 
                                    binsize = binsize,
                                    table = table[label]
                                    )

            thisrun = []
            for onoff in [0,1]:
                catalogs = [ cat[0] for cat in cats3 if cat[1] == onoff ]
                thisrun += [ np.any( [ countData[label]['catalog'] == catalog for catalog in catalogs ], axis = 0) ]
            table[label] = plotDistributionOutliers( pp2, ylabel = 'muon %s selection'%label, data = countData[label], ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['OFF', 'ON'],
                                    binsize = binsize,
                                    table = table[label]
                                    )
        

    pp2.close()
    os.rename( outfile2 + '.tmp', outfile2 )
    for key in keys:
        print 'saving', key
        for label in labels:
            savetable( table[label][key], output = 'plots/summary.muon.%s.%s.dat'%(label,key), ohdus = ohdus )
    #printTable(  )
    
def makePlotCu( fname ):
    if 0: plotMuon('plots/muonCount.pdf', labels = ['Juan45', 'Philipe45'] )
    if 0:
        print 'muon'

        outfile2 = 'plots/muonCount.pdf'
        pp2 = PdfPages( outfile2 + '.tmp')
        
        ohdus = [2,3,4,5,6,7,8,9,10]
        
        if not os.path.exists('plots/muonCount.Juan45.dat'):
            data = np.genfromtxt( 'plots/muon.Juan45.dat', names = True )
            computeCounts( data, output = 'muonCount.Juan45', ohdus = ohdus, mask = lambda x: np.ones_like(x, dtype = bool) )
        countData = np.genfromtxt( 'plots/muonCount.Juan45.dat', dtype=None, names = True )
        print countData['catalog']
        
        if not os.path.exists('plots/muonCount.Philipe45.dat'):
            data2 = np.genfromtxt( 'plots/muon.Philipe45.dat', names = True )
            computeCounts( data2, output = 'muonCount.Philipe45', ohdus = ohdus, mask = lambda x: np.ones_like(x, dtype = bool) )
        countData2 = np.genfromtxt( 'plots/muonCount.Philipe45.dat', dtype=None, names = True )
        print countData2['catalog']
        
        table = {}
        table2 = {}
        keys = ['all', 'run018', 'run025', 'OFF', 'ON']
        for key in keys:
            table[key] = {}
            table2[key] = {}
        binsize = 5
        
        for ohdu_ in ohdus:
            
            thisohdu = countData['ohdu'] == ohdu_
            plotTime2( pp2, ylabel = 'muon selection count', data = [countData, countData2], ohdu = ohdu_, labels = ['Juan45', 'Philipe45'] )
            table = plotDistributionOutliers( pp2, ylabel = 'muon Juan45 selection', data = countData, ohdu = ohdu_, masks = [np.ones_like(countData, dtype = bool)], labels = ['all'], binsize = binsize, table = table )

            table2 = plotDistributionOutliers( pp2, ylabel = 'muon Philipe45 count', data = countData2, ohdu = ohdu_, masks = [np.ones_like(countData2, dtype = bool)], labels = ['all'], binsize = binsize, table = table2 )
                
            thisrun = []
            for run in ['18','25']:
                catalogs = [ cat[0] for cat in cats3 if run in cat[0] ]
                thisrun += [ np.any( [ countData['catalog'] == catalog for catalog in catalogs ], axis = 0) ]

            table = plotDistributionOutliers( pp2, ylabel = 'muon Juan45 selection', data = countData, ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['run018', 'run025'], 
                                    binsize = binsize,
                                    table = table
                                    )
            
            thisrun = []
            for run in ['18','25']:
                catalogs = [ cat[0] for cat in cats3 if run in cat[0] ]
                thisrun += [ np.any( [ countData2['catalog'] == catalog for catalog in catalogs ], axis = 0) ]

            table2 = plotDistributionOutliers( pp2, ylabel = 'muon Philipe45 selection', data = countData2, ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['run018', 'run025'], 
                                    binsize = binsize,
                                    table = table2
                                    )
            

            thisrun = []
            for onoff in [0,1]:
                catalogs = [ cat[0] for cat in cats3 if cat[1] == onoff ]
                thisrun += [ np.any( [ countData['catalog'] == catalog for catalog in catalogs ], axis = 0) ]
            table = plotDistributionOutliers( pp2, ylabel = 'muon Juan45 selection', data = countData, ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['OFF', 'ON'],
                                    binsize = binsize,
                                    table = table
                                    )
            
            thisrun = []
            for onoff in [0,1]:
                catalogs = [ cat[0] for cat in cats3 if cat[1] == onoff ]
                thisrun += [ np.any( [ countData2['catalog'] == catalog for catalog in catalogs ], axis = 0) ]
            table2 = plotDistributionOutliers( pp2, ylabel = 'muon Philipe45 selection', data = countData2, ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['OFF', 'ON'],
                                    binsize = binsize,
                                    table = table2
                                    )

        pp2.close()
        os.rename( outfile2 + '.tmp', outfile2 )
        for key in keys:
            print 'saving', key
            savetable( table[key], output = 'plots/summary.muon.Juan45.%s.dat'%key, ohdus = ohdus )
            savetable( table2[key], output = 'plots/summary.muon.Philipe45.%s.dat'%key, ohdus = ohdus )
        #printTable(  )
        
    if 0:
        print 'high'
        data = np.genfromtxt( 'plots/peakCount.selectionHigh.dat', names = True )

        outfile2 = 'plots/peakAnalysis.stabilityHigh.pdf'
        pp2 = PdfPages( outfile2 + '.tmp')
        
        ohdus = [2,3,4,5,6,7,8,9,10]

        Emin = 250
        Emax = 350
        lowEnergymask = lambda thisdata: np.all( [ thisdata['E1'] < thisdata['gainCu']*Emax, thisdata['E1'] >= thisdata['gainCu']*Emin ], axis=0 )
        computeCounts( data, output = 'peakAnalysis.high', ohdus = ohdus, mask = lowEnergymask )
        countData = np.genfromtxt( 'plots/peakAnalysis.high.dat', dtype=None, names = True )
        print countData['catalog']

        table = {}
        keys = ['all', 'run018', 'run025', 'OFF', 'ON']
        for key in keys:
            table[key] = {}

        binsize = 5
        for ohdu_ in ohdus:
            thisohdu = countData['ohdu'] == ohdu_
            plotTime2( pp2, ylabel = 'count 250--350keV', data = [countData], ohdu = ohdu_, labels = [''] )
            table = plotDistributionOutliers( pp2, ylabel = 'count 250--350keV', data = countData, ohdu = ohdu_, masks = [np.ones_like(countData, dtype = bool)], labels = ['all'], binsize = binsize, table = table )
                
            thisrun = []
            for run in ['18','25']:
                catalogs = [ cat[0] for cat in cats3 if run in cat[0] ]
                thisrun += [ np.any( [ countData['catalog'] == catalog for catalog in catalogs ], axis = 0) ]

            table = plotDistributionOutliers( pp2, ylabel = 'count 250--350keV', data = countData, ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['run018', 'run025'], 
                                    binsize = binsize,
                                    table = table
                                    )
            thisrun = []
            for onoff in [0,1]:
                catalogs = [ cat[0] for cat in cats3 if cat[1] == onoff ]
                thisrun += [ np.any( [ countData['catalog'] == catalog for catalog in catalogs ], axis = 0) ]
            table = plotDistributionOutliers( pp2, ylabel = 'count 250--350keV', data = countData, ohdu = ohdu_, 
                                    masks = thisrun, 
                                    labels = ['OFF', 'ON'],
                                    binsize = binsize,
                                    table = table
                                    )

        pp2.close()
        os.rename( outfile2 + '.tmp', outfile2 )
        for key in keys:
            print 'saving', key
            savetable( table[key], subs = key.split('.'), output = 'plots/summary.peak.high.%s.dat'%key, ohdus = ohdus )
        

    
    data = np.genfromtxt( fname, names = True )

    outfile = 'plots/peakAnalysis.all.pdf'
    pp = PdfPages( outfile + '.tmp')
    outfile2 = 'plots/peakAnalysis.stability.pdf'
    pp2 = PdfPages( outfile2 + '.tmp')
    
    ohdus = [2,3,4,5,6,7,8,9,10]

    Emin = 3
    Emax = 7
    lowEnergymask = lambda thisdata: np.all( [ thisdata['E1'] < thisdata['gainCu']*Emax, thisdata['E1'] >= thisdata['gainCu']*Emin ], axis=0 )
    computeCounts( data, output = 'peakAnalysis.low', ohdus = ohdus, mask = lowEnergymask )
    countData = np.genfromtxt( 'plots/peakAnalysis.low.dat', dtype=None, names = True )
    print countData['catalog']


    table = {}
    keys = ['all', 'run018', 'run025', 'OFF', 'ON']
    for key in keys:
        table[key] = {}
        
    for ohdu_ in ohdus:
        thisohdu = countData['ohdu'] == ohdu_
        print 'low'
        if 0:
            fig_ = plt.figure()
            ax = fig_.add_subplot(111)
            ax.set_title('ohdu=%d'%ohdu_)
            ax.set_xlabel('runID')
            ax.set_ylabel('count 3--7keV')
            ax.plot( countData['runID'][thisohdu], countData['N'][thisohdu], '.' )
            fig_.savefig(pp2, format='pdf')


        if 1:
            plotTime2( pp2, ylabel = 'count 3--7keV', data = [countData], ohdu = ohdu_, labels = [''] )
            
        if 1:
            table = plotDistributionOutliers( pp2, ylabel = 'count 3--7keV', data = countData, ohdu = ohdu_, masks = [np.ones_like(countData, dtype = bool)], labels = ['all'], table = table )
            
            thisrun = []
            for run in ['18','25']:
                catalogs = [ cat[0] for cat in cats3 if run in cat[0] ]
                thisrun += [ np.any( [ countData['catalog'] == catalog for catalog in catalogs ], axis = 0) ]

            table = plotDistributionOutliers( pp2, ylabel = 'count 3--7keV', data = countData, ohdu = ohdu_, 
                                     masks = thisrun, 
                                     labels = ['run018', 'run025'],
                                     table = table
                                     )

            thisrun = []
            for onoff in [0,1]:
                catalogs = [ cat[0] for cat in cats3 if cat[1] == onoff ]
                thisrun += [ np.any( [ countData['catalog'] == catalog for catalog in catalogs ], axis = 0) ]
                
            table = plotDistributionOutliers( pp2, ylabel = 'count 3--7keV', data = countData, ohdu = ohdu_, 
                                     masks = thisrun, 
                                     labels = ['OFF', 'ON'],
                                     table = table
                                     )
            print 'table', table.keys()
            
        if 0:
            print 'ohdu', ohdu_
            fig2 = plt.figure()
            fig6 = plt.figure()

            bins = 200
            N=11
            N=3
            bins2 = [bins, len(set(runID)) ]
            #ax = fig.add_subplot(N,1,1)
            #ax.set_yscale('log')
            #ax.set_title('ohdu=%d'%ohdu_)
            
            Emin = 1e3
            Emax = 20e3
            #print Emin, Emax
            nMin = 4
            print 'nMin=', nMin
            
            #Ecount, Eedges, _ = ax.hist( E1_, bins = bins, histtype = 'step', label = 'all' )
            
            fig2ax = fig2.add_subplot(111)
            fig2ax.set_title(r"ohdu=%d"%ohdu_)
            fig2ax.set_xlabel(r"E[adu]")
            fig2ax.set_ylabel(r"dN/dE [1/adu]")
            fig2ax.grid(True)
            fig2ax.hist( E1_, bins = bins, histtype = 'step', label = 'all', normed=True )
            fig2ax.hist( E1_[ n0_ <= nMin ], bins = bins, histtype = 'step', label = 'n0 <= %d'%nMin, normed=True )
            #fig2ax.hist( catData[1]['E1'][ catData[1]['n0'] <= nMin ], bins = bins, histtype = 'step', label = 'n0 $<= %d$ cat'%nMin, normed=True )
            #fig2ax.hist( E1_[ np.all( [ n0_ <= nMin, runID_ == runID_[0] ], axis=0) ], bins = bins, histtype = 'step', label = 'n0 $<= %d$ runID'%nMin, normed=True )
            fig2ax.legend( fancybox=True, framealpha=0)
            
            
            #ax.axvline( np.mean(gainCu_)*CuPeak, color='g' )
            #ax.axvline( np.mean(gainCu_)*CuPeak2, color='g' )
            #ax.axvline( np.mean(gainCu_)*SiPeak, color='r' )
            #ax.set_xlim([Emin,Emax])
            minVar1 = .15
            ratio = .4
            
            maskC = np.all( [ xVar1_>.01, np.abs(yVar1_/xVar1_ - 1.) > ratio, np.abs(yVar1_/xVar1_ - 1.) < 1, xVar1_ < minVar1, yVar1_ < minVar1 ], axis=0 )
            maskR = np.all( [ xVar1_>.01, np.abs(yVar1_/xVar1_ - 1.) > ratio, np.abs(yVar1_/xVar1_ - 1.) < 1 ], axis=0 )
            maskM = np.all( [ xVar1_ < minVar1, yVar1_ < minVar1 ], axis=0 )
        
            ax5 = fig2ax.twinx()
            ax5.set_ylabel('runID')
        
            meanrunID = []
            meanCu, stdCu = [], []
            meanCu2, stdCu2 = [], []
            meanSi, stdSi = [], []
            meangainCu, meanresCu = [], []

            meanCuA, stdCuA = [], []
            meanCu2A, stdCu2A = [], []
            meanSiA, stdSiA = [], []

            meanCuA, stdCuA = [], []
            meanCu2A, stdCu2A = [], []
            meanSiA, stdSiA = [], []
            runID2 = []
            
            print 'peaks for cats'
            for catData_ in catData:
                E1__ = catData_['E1']
                n0__ = catData_['n0']
                gainCu__ = catData_['gainCu']
                resCu__ = catData_['resCu']
                runID__ = catData_['runID']
                
                meanrunID += [ np.mean(runID__) ]

                meangainCu += [ np.mean(gainCu__) ]
                meanresCu += [ np.mean(resCu__) ]
                mean, std = norm.fit( E1__[ n0__<=nMin ] )
                
                meanSi_, stdSi_ = norm.fit( E1__[ np.all([ n0__<=nMin, E1__< mean ], axis=0) ] )
                meanCu_, stdCu_ = norm.fit( E1__[ np.all([ n0__<=nMin, E1__> mean ], axis=0) ] )
                
                fraction = .8
                while 1:
                    meanCu_2, stdCu_2 = norm.fit( E1__[ np.all([ n0__<=nMin, abs(E1__ - meanCu_)< 1.5*stdCu_ ], axis=0) ] )
                    if stdCu_2 > fraction*stdCu_: break
                    meanCu_, stdCu_ = meanCu_2, stdCu_2
                meanCu += [ meanCu_2 ]
                stdCu += [ stdCu_2 ]
                
                meanSi_ = meanCu_*SiPeak/CuPeak
                stdSi_ = .5e3
                while 1:
                    meanSi_2, stdSi_2 = norm.fit( E1__[ np.all([ n0__<=nMin, abs(E1__ - meanSi_)< 1.5*stdSi_ ], axis=0) ] )
                    if stdSi_2 > fraction*stdSi_: break
                    meanSi_, stdSi_ = meanSi_2, stdSi_2
                    
                meanSi += [ meanSi_2 ]
                stdSi += [ stdSi_2 ]
                
                meanCu2_ = meanCu_*CuPeak2/CuPeak
                stdCu2_ = .5e3

                #meanCu2_, stdCu2_ = norm.fit( E1__[ np.all([ n0__<=nMin, E1__ > meanCu_2 + stdCu_2 ], axis=0) ] )
                while 1:
                    #meanCu2_2, stdCu2_2 = norm.fit( E1__[ np.all([ n0__<=nMin, E1__ > 1.05*meanCu_2, abs(E1__ - meanCu2_)< 1.5*stdCu2_ ], axis=0) ] )
                    meanCu2_2, stdCu2_2 = norm.fit( E1__[ np.all([ n0__<=nMin, abs(E1__ - meanCu2_)< 1.5*stdCu2_ ], axis=0) ] )
                    if stdCu2_2 > fraction*stdCu2_: break
                    meanCu2_, stdCu2_ = meanCu2_2, stdCu2_2
                
                meanCu2 += [ meanCu2_2 ]
                stdCu2 += [ stdCu2_2 ]

                for runid in runID__:
                    runID2 += [runid]
                    r = norm.fit( E1__[ np.all([ n0__<=nMin, runID__==runid, abs(E1__ - meanSi_)< .5e3 ], axis=0) ] )
                    meanSiA += [r[0]]
                    stdSiA += [r[1]]
                    r = norm.fit( E1__[ np.all([ n0__<=nMin, runID__==runid, abs(E1__ - meanCu_)< .5e3 ], axis=0) ] )
                    meanCuA += [r[0]]
                    stdCuA += [r[1]]
                    r = norm.fit( E1__[ np.all([ n0__<=nMin, runID__==runid, abs(E1__ - meanCu2_)< .5e3 ], axis=0) ] )
                    meanCu2A += [r[0]]
                    stdCu2A += [r[1]]
                
            print 'plotting'
            #print meanSiA, stdSiA
            
            ax5.errorbar( np.array(meangainCu)*CuPeak, meanrunID, xerr = np.array(meanresCu)*CuPeak, fmt='o', label = '*Cu=8.04keV' )
            ax5.plot( np.array(meangainCu)*CuPeak2, meanrunID, 'o', label = '*Cu2=8.9keV' )
            ax5.plot( np.array(meangainCu)*SiPeak, meanrunID, 'o', label = '*Si=1.7keV' )
            
            ax5.errorbar( meanCu, meanrunID, xerr=stdCu, fmt='.', label = 'Cu' )
            ax5.errorbar( meanCu2, meanrunID, xerr=stdCu2, fmt='.', label = 'Cu2' )
            ax5.errorbar( meanSi, meanrunID, xerr=stdSi, fmt='.', label = 'Si' )
            ax5.legend( fancybox=True, framealpha=0, loc=4, title="estimates per catalog" )
            
            ax6 = fig6.add_subplot(111)
            ax6.set_title('compare estimates per runID for ohdu=%d'%ohdu_)
            ax6.set_xlabel('runID')
            ax6.set_ylabel('E[adu]')
            ax6.grid(True)

            
            #ax6.errorbar( runID2, meanSiA, yerr=stdSiA, fmt='.', label='Si' )
            #ax6.errorbar( runID2, meanCuA, yerr=stdCuA, fmt='.', label='Cu' )
            #ax6.errorbar( runID2, meanCu2A, yerr=stdCu2A, fmt='.', label='Cu2' )

            ms = 1
            ax6.plot( runID2, meanSiA, '.', ms = ms, label='Si' )
            ax6.plot( runID2, meanCuA, '.', ms = ms, label='Cu' )
            ax6.plot( runID2, meanCu2A, '.', ms = ms, label='Cu2' )

            ax6.plot( runID_, gainCu_*SiPeak, '.', ms = ms, label = '*Si' )
            #ax6.errorbar( runID_, gainCu_*CuPeak, yerr = resCu_*CuPeak, fmt='.', label = '*Cu' )
            ax6.plot( runID_, gainCu_*CuPeak,'.', ms = ms, label = '*Cu' )
            ax6.plot( runID_, gainCu_*CuPeak2, '.', ms = ms, label = '*Cu2' )

            ax6.legend( fancybox=True, framealpha=0 )

            
            if 0:
                fig7 = plt.figure()
                ax7 = fig7.add_subplot(111)
                ax7.set_title('relative difference of the estimates for ohdu=%d'%ohdu_)
                ax7.set_xlabel('runID')
                ax7.set_ylabel(r'(N-N*)/N*')
                ax7.grid(True)
                ax7.plot( runID2, (meanSiA - gainCu_[runID_==runID2]*SiPeak)/(gainCu_[runID_==runID2]*SiPeak), '.', ms = ms, label = 'Si' )
                ax7.plot( runID2, (meanCuA - gainCu_[runID_==runID2]*CuPeak)/(gainCu_[runID_==runID2]*CuPeak),'.', ms = ms, label = 'Cu' )
                ax7.plot( runID2, (meanCu2A - gainCu_[runID_==runID2]*CuPeak2)/(gainCu_[runID_==runID2]*CuPeak2), '.', ms = ms, label = 'Cu2' )
                ax7.legend( fancybox=True, framealpha=0 )
                fig7.savefig(pp,  format='pdf')
            
            fig2.savefig(pp, bbox_inches='tight',  format='pdf')
            fig6.savefig(pp, bbox_inches='tight',  format='pdf')
            
    pp.close()
    pp2.close()
    os.rename( outfile + '.tmp', outfile )
    os.rename( outfile2 + '.tmp', outfile2 )
    for key in keys:
        print 'saving', key
        savetable( table[key], subs = key.split('.'), output = 'plots/summary.peak.low.%s.dat'%key, ohdus = ohdus )
    
    outfile = 'plots/summary.peak.low.pdf'
    pp3 = PdfPages( outfile + '.tmp')
    #printTable( [table['all']], labels = ['all'], ohdus = ohdus, ylabel = '3--7keV', output = pp3 )
    printTable( [table['all'], table['ON'], table['OFF']], labels = ['all', 'ON', 'OFF'], ohdus = ohdus, ylabel = '3--7keV', output = pp3 )
    printTable( [table['all'], table['run018'], table['run025']], labels = ['all', 'run018', 'run025'], ohdus = ohdus, ylabel = '3--7keV', output = pp3 )
    pp3.close()
    os.rename( outfile + '.tmp', outfile )
    #np.savetxt( 'plots/peakAnalysis.stability.dat', zip(*[ table[key] for key in tableorder ]), header=' '.join( tableorder ), fmt='%s')
               #fmt=' '.join([ '%d' if key == 'ohdu' else ( '%s' if key == 'runID' else '%.4f' ) for key in table.keys() ]) )
    return

if __name__ == "__main__":
    catalogs = "/share/storage2/connie/nu_processing/scripts/ProcCat/cut_scn_osi_raw_gain_catalog_data_*.skim1.root"
    if len(sys.argv) > 1:
        if sys.argv[1] == '--h':
            print "usage:"
            print "--h\t\tprints this message"
            print "--muon\t\tprocesses the muon counts"
            print "--muonPlot\tgenetares the muon count plots"
            print "--peak\t\tprocess the peaks info"
            print "--peakPlot\tgenerates the peaks plots"
            exit(0)
        if sys.argv[1] == "--muon":
            print "Processes the muon count"
            catalogs = sys.argv[-1]
            catalogs = glob.glob(catalogs)
            if os.path.splitext( catalogs[0] )[1] == ".root":
                print "trees", root_numpy.list_trees(catalogs[0])
                print "branches", root_numpy.list_branches(catalogs[0], 'hitSumm')
                print "branches", root_numpy.list_branches(catalogs[0], 'config')
                applyToAll( catalogs )
            exit(0)
        if sys.argv[1] == "--muon2":
            print "Processes the muon2 count"
            #catalogs = sys.argv[-1]
            #catalogs = "/share/storage2/connie/data_analysis/processed02_data/runs/*/*/ext/catalog/*.skim1.root"
            catalogs = glob.glob(catalogs)
            branches = ['runID','ohdu','expoStart','E1','gainCu', 'resCu', 'selFlag']
            #branches = ['runID','ohdu','expoStart','E1']
            fmt = ' '.join( ['%d','%d','%d','%.4e','%.4e', '%.4e', '%d'] )
            #fmt = ' '.join( ['%d','%d','%d','%.4e'] )
            if os.path.splitext( catalogs[0] )[1] == ".root":
                for key in ['Juan45', 'Philipe45']:
                    applySelection( catalogs, selection=selections[key], branches=branches, fmt=fmt, output='muon.%s'%key )
            exit(0)
        if sys.argv[1] == "--muon2Plot":
            print "plots the muon count"
            data = sys.argv[2]
            if os.path.splitext( data )[1] == ".dat":
                name = 'plots/new.muonCount.all.pdf'
                pp = PdfPages(name+'.tmp')
                plotTime( fnames = [data], output = pp, branch='E1', ohdu=2, labels=['Philipe'], func = len, ylabel = 'count' )
                pp.close()
                os.rename( name+'.tmp', name )
            exit(0)

        if sys.argv[1] == "--muonPlot":
            print "Processes the muon count"
            data = sys.argv[2]
            if os.path.splitext( data )[1] == ".dat":
                makePlot( data )
            exit(0)
        if sys.argv[1] == "--peak":
            print "Processes the peak count"
            catalogs = sys.argv[-1]
            catalogs = glob.glob(catalogs)
            if os.path.splitext( catalogs[0] )[1] == ".root":
                print "trees", root_numpy.list_trees(catalogs[0])
                print "branches", root_numpy.list_branches(catalogs[0], 'hitSumm')
                print "branches", root_numpy.list_branches(catalogs[0], 'config')
                branches = ['E0', 'E1', 'ohdu','runID', 'selFlag', 'expoStart','xVar1','yVar1', 'gainCu', 'resCu','n0','n1', 'n2','n3']
                fmt = ['%.8e', '%.8e', '%d', '%d', '%d', '%d', '%.8e', '%.8e', '%.8e', '%.8e', '%d', '%d', '%d', '%d']
                print len(branches), len(fmt)
                fmt = ' '.join(fmt)
                Emin = 1e3#adu
                Emax = 20e3#adu
                selection = 'flag==0 && dcFlag==1 && E1> %f && E1 < %f'%(Emin,Emax)
                #applySelection( catalogs, selection = selection, branches = branches, output='peakCount.selection', fmt = fmt )
                
                Emin = 200#keV
                Emax = 400#keV
                selection = 'flag==0 && dcFlag==1 && E1> gainCu*%f && E1 < gainCu*%f'%(Emin,Emax)
                applySelection( catalogs, selection = selection, branches = branches, output='peakCount.selectionHigh', fmt = fmt )
            exit(0)
        if sys.argv[1] == "--peakPlot":
            if len(sys.argv) < 3: 
                print "no file given"
                exit(0)
            print "Processes the peak plots"
            data = sys.argv[2]
            if os.path.splitext( data )[1] == ".dat":
                makePlotCu( data )
            exit(0)
        if sys.argv[1] == "--generic":
            if len(sys.argv) < 6: 
                print "no file given"
                exit(0)
            catalogs = glob.glob(sys.argv[2])
            print 'catalogs:'
            print catalogs
            selection = sys.argv[3]
            print 'selection:'
            print selection
            branches = sys.argv[4].split(',')
            print 'branches:'
            print branches
            output = sys.argv[5]
            print 'output:'
            print output
            applySelection( catalogs, selection, branches, output )
            exit(0)
        print sys.argv[1]," is not a known option"
        exit(0)
    else:
        print "no args given"
