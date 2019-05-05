import root_numpy
import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange, epoch2num
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm, chisquare

from scipy.optimize import minimize
from scipy.misc import factorial

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

def summary( var_str ):
    print [ '%s: '%k + str(v) for k,v in globals().iteritems() if k == var_str ][0]

selectionPhilipe = r'flag==0 && abs(log10( E0/sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2 + 675**2  ) )-2.7055039797010396)<2*sqrt(7.326161363760477-2.7055039797010396**2)'

selectionJuan = r'E1>0 && flag==0 && sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))>20 && abs(n0/sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))-4.4)<1'

selectionJuanAngle = r'E1>0 && flag==0 && sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))>20 && abs(n0/sqrt(pow(xMax-xMin,2)+pow(yMax-yMin,2))-4.4)<1 && atan2( sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2), 675 ) > 0.785398'

selectionPhilipeAngle = r'flag==0 && abs(log10( E0/sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2 + 675**2  ) )-2.7055039797010396)<2*sqrt(7.326161363760477-2.7055039797010396**2) && atan2( sqrt( (15*(xMax-xMin))**2 + (15*(yMax-yMin))**2), 675 ) > 0.785398'

def getValuesFromCatalog( catalog, treename = 'hitSumm', branches = [], selection = '' ):
    if branches is []: return None
    data = root_numpy.root2array( catalog, treename = treename, branches = branches, selection = selection ).view(np.recarray)
    #return data
    return np.array([ data[branch] for branch in branches ])

def applyToAll( catalogs, selection = '' ):
    nJ = {'all':[]}
    nJa = {'all':[]}
    nPh = {'all':[]}
    nPha = {'all':[]}
    runID_ = []
    expoStart_ = []
    print 'catalogs', catalogs
    ohdus = [2,3,4,5,6,8,9]
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
        E0J, ohduJ, expoStartJ, runIDJ = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionJuan )
        E0Ja, ohduJa, expoStartJa, runIDJa = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionJuanAngle )
        E0Ph, ohduPh, expoStartPh, runIDPh = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionPhilipe )
        E0Pha, ohduPha, expoStartPha, runIDPha = getValuesFromCatalog( catalog, branches = ['E0','ohdu','expoStart','runID'], selection = selectionPhilipeAngle )
        print 'done'
        
        for runID in runIDs:
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
    
    labels = {'nJa':r'Juan($\theta>45$)', 'nPha':r'Philipe2($\theta>45$)'}
    count = 0
    #for ohdu in ['all'] + range(20):
    for ohdu in [2,3,5,9]:
        fig = plt.figure()
        fig2 = plt.figure()
        
        size = fig.get_size_inches()
        size = [ size[0], size[1]*(1+(len(data.dtype.names)-1)/3) ]
        fig.set_size_inches(size)
        print ohdu
        
        ax = fig.add_subplot( 1+(len(data.dtype.names)-1)/3, 1 , 1+count)
        ax.set_ylabel(r"$N_\mu$")
        ax.set_xlabel(r"time")
        ax.set_title(r"ohdu=%s"%ohdu)
        ax.grid(True)
        
        #bx2 = bx.twinx()

        count += 1

        nbins = 50
        y_all = {}
        y_mean = {}
        usedLabels = []
        for index,edges in enumerate(cats2):
            for j, key in enumerate( ['nJa','nPha'] ):
                if not key in y_all: y_all[key] = []
                if not key in y_mean: y_mean[key] = []

                y = catData[index]["%s%s"%(key,ohdu)]
                xx = catData[index]['expoStart']
                x = epoch2num(xx)
                y_all[key] = np.append( y_all[key], y )

                ax.plot_date( x, y, '.', alpha=.1 )

                std_ = np.std(y)
                median_ = np.median(y)
                
                y_mean[key] = np.append( y_mean[key], median_ )

                removeoutliers = abs(y - median_) < 2*std_
                color = ['k','r'][j]
                ax.errorbar( epoch2num( .5*(min(xx)+max(xx)) ), np.median(y), yerr = std_ , color=color, fmt='o', alpha=.2 )
                ax.errorbar( epoch2num( .5*(min(xx)+max(xx)) ), np.median(y[removeoutliers]), yerr = np.std(y[removeoutliers]) , color=color, fmt='o', label = labels[key] if not labels[key] in usedLabels else '' )
                usedLabels += [labels[key]]

        ax.legend( fancybox=True, framealpha=0.1)
        dateFormatter = DateFormatter('%d %b %y')
        ax.xaxis.set_major_formatter( dateFormatter )
        fig.autofmt_xdate()
        fig.savefig(pp, bbox_inches='tight', format='pdf')
        
        bx1 = fig2.add_subplot(121)
        bx2 = fig2.add_subplot(122)
        bx1.set_xlabel(r'$N_\mu$')
        bx1.set_ylabel(r'count')
        bx1.set_title(r'distribution of $N_\mu$ ohdu=%d'%ohdu)

        bx2.set_xlabel(r'$N_\mu$')
        bx2.set_ylabel(r'count')
        bx2.set_title(r'distribution of $N_\mu$ per catalog ohdu=%d'%ohdu)
        
        for index, key in enumerate(y_all.keys()):
            bx = [ bx1, bx2 ][index]
            mean, std = norm.fit( y_all[key] )
            count1, edges1, _ = bx.hist( y_all[key], bins = nbins, histtype='step', label = labels[key] )
            #bx.plot( [],[], ' ',label = r'$\mu=%.1f,$'%mean + r'$\sigma=%.1f$'%std )
            bx.legend(  fancybox=True, framealpha=0.1, title=r'$\mu=%.1f,$'%mean + r'$\sigma=%.1f$'%std )
            
        fig2.savefig(pp, bbox_inches='tight', format='pdf')

    pp.close()


def applyToAllCu( catalogs, selection = '' ):
    E1,ohdu,runID,xVar1,yVar1, gainCu, resCu, n0, n1, n2, n3 = [], [], [], [], [], [], [], [], [], [], []
    
    Emin = 1e3
    Emax = 20e3
    count = 0
    for catalog in catalogs:
        print catalog
        E1_, ohdu_, runID_, xVar1_, yVar1_, gainCu_, resCu_, n0_, n1_, n2_, n3_ = getValuesFromCatalog( catalog, branches = ['E0','ohdu','runID','xVar1','yVar1','gainCu','resCu','n0','n1', 'n2','n3'], selection = 'flag == 0 && E1> %f && E1 < %f'%(Emin,Emax) )
        
        E1 = np.append( E1, E1_)
        gainCu = np.append( gainCu, gainCu_)
        resCu = np.append( resCu, resCu_)
        ohdu = np.append( ohdu, ohdu_ )
        runID = np.append( runID, runID_ )
        xVar1 = np.append( xVar1, xVar1_ )
        yVar1 = np.append( yVar1, yVar1_ )
        n0 = np.append( n0, n0_ )
        n1 = np.append( n1, n1_ )
        n2 = np.append( n2, n2_ )
        n3 = np.append( n3, n3_ )
        print E1.shape

    data = [ohdu, runID, E1, xVar1, yVar1, gainCu, resCu, n0, n1, n2, n3 ]
    header = "ohdu runID E1 xVar1 yVar1 gainCu resCu n0 n1 n2 n3"
    np.savetxt( 'plots/peakCount.dat', zip(*data), header = header, fmt='%f' )


def poisson(k, lamb, N=1):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)


def negLogLikelihood(params, data):
    """ the negative log-Likelohood-Function"""
    lnl = - np.sum(np.log(poisson( data, params[0] )))
    return lnl

def compute_chisquare(f1, f2):
    result = np.sum( ( f1-f2)**2/f2 )
    return result, result/len(f1)

def makePlotCu( fname ):
    data = np.genfromtxt( fname, names = True )
    print 'gain', len(data['gainCu'])
    print 'res', len(data['resCu'])
    print 'E1', len(data['E1'])
    outfile = 'plots/peakAnalysis.all.pdf'
    pp = PdfPages( outfile + '.tmp')
    outfile2 = 'plots/peakAnalysis.stability.pdf'
    pp2 = PdfPages( outfile2 + '.tmp')
    
    
    E1 = data['E1']
    ohdu = data['ohdu']
    ohdus = [2,3,5,9]
    runID = data['runID'] 
    gainCu = data['gainCu']
    resCu = data['resCu']
    n0 = data['n0']
    n1 = data['n1']
    n2 = data['n2']
    n3 = data['n3']
    xVar1 = data['xVar1']
    yVar1 = data['yVar1']
    
    print E1.shape
    
    for ohdu_ in ohdus:
        catData = [ data[ np.all( [ cat[-2]<data['runID'], data['runID']<=cat[-1], data['ohdu'] == ohdu_ ], axis=0) ] for cat in cats2 ]
        E1_ = E1[ ohdu == ohdu_ ]
        runID_ = runID[ ohdu == ohdu_ ]
        gainCu_ = gainCu[ ohdu == ohdu_ ]
        resCu_ = resCu[ ohdu == ohdu_ ]
        xVar1_ = xVar1[ ohdu == ohdu_ ]
        yVar1_ = yVar1[ ohdu == ohdu_ ]
        n0_ = n0[ ohdu == ohdu_ ]
        n1_ = n1[ ohdu == ohdu_ ]
        n2_ = n2[ ohdu == ohdu_ ]
        n3_ = n3[ ohdu == ohdu_ ]
        
        #count = [ E1_[ np.all( [ runID_ == runid, E1_ < gainCu_*7, E1_ > gainCu_*3 ], axis = 0) ].shape[0] for runid in runID_ ]
        #print len(E1), len(E1_)
        runIDs = list(set(runID_))
        count = [ len(E1_[ np.all( [ runID_ == runid, E1_ < gainCu_*7, E1_ > gainCu_*3 ], axis = 0) ]) for runid in runIDs ]
        #print len(E1_[ np.all( [ E1_ < gainCu_*7, E1_ > gainCu_*3 ], axis = 0) ])
        #print count
        #print len(count)
        #print set(runID_)
        print len(runIDs), len(count)
        if 1:
            fig_ = plt.figure()
            ax = fig_.add_subplot(111)
            ax.set_title('ohdu=%d'%ohdu_)
            ax.set_xlabel('runID')
            ax.set_ylabel('count (3keV, 7keV)')
            ax.plot( runIDs, count, '.' )
            fig_.savefig(pp2, format='pdf')
            
        if 1:
            fig_ = plt.figure()
            ax = fig_.add_subplot(111)
            ax.set_title('ohdu=%d'%ohdu_)
            ax.set_xlabel('count (3keV, 7keV)')
            ax.set_ylabel('count')
            mean, std = norm.fit( count )
            print 'mean, std', mean, std
            h, edges, _ = ax.hist( count, bins = np.arange(0, np.max(count), 1)-.5, histtype='step', label = '3--7keV count' )
            #result = minimize(negLogLikelihood,  # function to minimize
                  #x0=np.ones(1),     # start value
                  #args=(count,),      # additional arguments for function
                  #method='Powell',   # minimization method, see docs
                  #)
            #print 'result of the minimization', result
            x = .5*(edges[1:] + edges[:-1])
            print 'poisson', compute_chisquare( h, np.sum(h)*poisson( x, mean )  ), compute_chisquare( np.sum(h)*poisson( x, mean ), h )
            print 'gauss', compute_chisquare( h, np.sum(h)*norm.pdf( x, mean, std ) ), compute_chisquare( np.sum(h)*norm.pdf( x, mean, std ), h )
            print chisquare(h, f_exp=np.sum(h)*poisson( x, mean ))
            print chisquare(h, f_exp=np.sum(h)*norm.pdf(x, mean, std))
            print len(x), len(count)
            
            ax.plot( x, np.sum(h)*poisson( x, mean ), '-', label = 'poisson chi2r=%.2f'%( 1./len(x)*chisquare(h, f_exp=np.sum(h)*poisson( x, mean ) )[0] ) )
            ax.plot( x, np.sum(h)*norm.pdf(x, mean, std), '-', label = 'gaussian chi2r=%.2f'%( 1./len(x)*chisquare(h, f_exp=np.sum(h)*norm.pdf(x, mean, std) )[0] ) )
            
            ax.legend( fancybox=True, framealpha=0.1, title=r'mean=%.1f,'%(mean) + r'std=%.1f'%std )
            fig_.savefig(pp2, format='pdf')
        
        
        print 'ohdu', ohdu_
        #fig = plt.figure()
        fig2 = plt.figure()
        fig6 = plt.figure()
        #fig3 = plt.figure()
        #fig4 = plt.figure()
        
        #sizex, sizey = fig.get_size_inches()
        #fig.set_size_inches( sizex, sizey*3 )

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
        fig2ax.hist( E1_[ n0_ <= nMin ], bins = bins, histtype = 'step', label = 'n0 $<= %d$'%nMin, normed=True )
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
        
        #for minVar in [ minVar1 ]:
            #ax.hist( E1_[maskM], bins = bins, histtype = 'step', label = 'Var1$<$%.2f'%minVar )
            #ax.hist( E1_[maskR], bins = bins, histtype = 'step', label = '$.4<yVar1/xVar1<1$' )
            #ax.hist( E1_[maskC], bins = bins, histtype = 'step', label = 'combined', lw = 4, color = 'k' )
        #ax.legend( fancybox=True, framealpha=0)
        
        ###ax = fig.add_subplot(N,1,2)
        ###ax.set_ylabel('runID')
        ###ax.set_title('all')
        ###ax.hist2d( E1_, runID_, bins = bins2, label = '', norm=LogNorm(), cmap = plt.cm.rainbow )
        ###ax.plot( gainCu_*CuPeak, runID_, 'g.' )
        ###ax.plot( gainCu_*CuPeak2, runID_, 'g.' )
        ###ax.plot( gainCu_*SiPeak, runID_, 'r.' )
        ###ax.set_xlim([Emin,Emax])

        ###ax = fig.add_subplot(N,1,3)
        ###ax.set_title('n0$<=%d$'%nMin)
        ###mask = n0_ <= nMin
        ###ax.hist2d( E1_[mask], runID_[mask], bins = bins2, label = '', norm=LogNorm(), cmap = plt.cm.rainbow )
        ####ax.plot( gainCu_*CuPeak, runID_, 'g.' )
        ####ax.plot( gainCu_*CuPeak2, runID_, 'g.' )
        ####ax.plot( gainCu_*SiPeak, runID_, 'r.' )
        ###ax.set_xlim([Emin,Emax])

        #ax2 = fig.add_subplot(N,1,3)
        #ax2.set_ylabel('runID')
        #ax2.set_title(r'xVar1, yVar1 $< %.2f$'%minVar)
        #ax2.hist2d( E1_[ maskM ], runID_[ maskM ], bins = bins2, label = '', norm=LogNorm(), cmap = plt.cm.rainbow )
        #ax2.plot( gainCu_[maskM]*CuPeak, runID_[maskM], 'g.' )
        #ax2.plot( gainCu_[maskM]*CuPeak2, runID_[maskM], 'g.' )
        #ax2.plot( gainCu_[maskM]*SiPeak, runID_[maskM], 'r.' )
        #ax2.set_xlim([Emin,Emax])

        ##print 'map3'
        #ax2 = fig.add_subplot(N,1,4)
        #ax2.set_ylabel('runID')
        ##ax2.set_xlabel('adu')
        #ax2.set_title(r'$.4 < abs(yVar1/xVar1)<1$')
        #ax2.hist2d( E1_[ maskR ], runID_[ maskR ], bins = bins2, label = '', norm=LogNorm(), cmap = plt.cm.rainbow )
        #ax2.plot( gainCu_[maskR]*CuPeak, runID_[maskR], 'g.' )
        #ax2.plot( gainCu_[maskR]*CuPeak2, runID_[maskR], 'g.' )
        #ax2.plot( gainCu_[maskR]*SiPeak, runID_[maskR], 'r.' )
        #ax2.set_xlim([Emin,Emax])

        #ax2 = fig.add_subplot(N,1,5)
        #ax2.set_ylabel('runID')
        #ax2.set_xlabel('adu')
        #ax2.set_title('combined')
        #ax2.hist2d( E1_[ maskC ], runID_[ maskC ], bins = bins2, label = '', norm=LogNorm(), cmap = plt.cm.rainbow )
        #ax2.plot( gainCu_[maskC]*CuPeak, runID_[maskC], 'g.' )
        #ax2.plot( gainCu_[maskC]*CuPeak2, runID_[maskC], 'g.' )
        #ax2.plot( gainCu_[maskC]*SiPeak, runID_[maskC], 'r.' )
        #ax2.set_xlim([Emin,Emax])
        
        #####logyVar1_ = np.log10(yVar1_)
        #####mask = np.isfinite(logyVar1_)
        #####ax3 = fig3.add_subplot(3,1,1)
        #####ax3.set_ylabel('log(yVar1)')
        #####ax3.set_xlabel('E1')
        #####ax3.hist2d( E1_[mask], logyVar1_[mask], bins = [bins,30], norm=LogNorm(), cmap = plt.cm.rainbow )
        #####ax3.set_xlim([Emin,Emax])

        #####logyVar1_ = np.log10(yVar1_)
        #####mask = np.all( [ np.isfinite(logyVar1_), n0_ <= nMin ], axis = 0 )
        #####ax3 = fig3.add_subplot(3,1,2)
        #####ax3.set_ylabel('log(yVar1)')
        #####ax3.set_xlabel('E1')
        #####ax3.hist2d( E1_[mask], logyVar1_[mask], bins = [bins,30], norm=LogNorm(), cmap = plt.cm.rainbow )
        #####ax3.set_xlim([Emin,Emax])
        
        #mask = np.isfinite(logyVar1_)
        #Ex = [ .5*(Eedges[i]+Eedges[i+1]) for i in range(len(Eedges)-1) ]
        #std_ = [ np.std( .5*(xVar1_+yVar1_)[ np.all( [ E1_>Eedges[i], E1_<Eedges[i+1]], axis = 0 ) ] ) for i in range(len(Eedges)-1) ]
        #ax3 = fig.add_subplot(N,1,7)
        #ax3.set_ylabel('std((yVar1+xVar1)/2)')
        #ax3.set_xlabel('E1')
        #ax3.plot( Ex, std_, '.' )
        #ax3.set_xlim([Emin,Emax])

        #mask = np.isfinite(logyVar1_)
        #Ex = [ .5*(Eedges[i]+Eedges[i+1]) for i in range(len(Eedges)-1) ]
        #stdlog_ = [ np.std( np.log10(yVar1_)[ np.all( [ E1_>Eedges[i], E1_<Eedges[i+1], np.isfinite(np.log10(xVar1_)) ], axis = 0 ) ] ) for i in range(len(Eedges)-1) ]
        #ax3 = fig.add_subplot(N,1,8)
        #ax3.set_ylabel('std(log(yVar1))')
        #ax3.set_xlabel('E1')
        #ax3.plot( Ex, stdlog_, '.' )
        #ax3.set_xlim([Emin,Emax])

        ####mask = np.isfinite(np.log10(yVar1_/xVar1_))
        ####ax3 = fig3.add_subplot(3,1,3)
        ####ax3.set_ylabel('log(yVar1/xVar1)')
        ####ax3.set_xlabel('E1')
        ####ax3.hist2d( E1_[mask], np.log10(yVar1_/xVar1_)[mask] , bins = [bins,30], norm=LogNorm(), cmap = plt.cm.rainbow )
        ####ax3.set_xlim([Emin,Emax])

        ####for index, entry in enumerate( [ ['n0', n0_], ['n1', n1_], ['n2', n2_], ['n3', n3_] ] ):
            ####label, value = entry
            ####ax3 = fig4.add_subplot(4,1,1+index)
            ####ax3.set_ylabel(label)
            ####ax3.set_xlabel('E1')
            ####ax3.hist2d( E1_, np.log10( value ) , bins = [bins,30], norm=LogNorm(), cmap = plt.cm.rainbow )
            ####ax3.set_xlim([Emin,Emax])

        #fig5 = plt.figure()
        #ax5 = fig5.add_subplot( 111 )
        ax5 = fig2ax.twinx()
        ax5.set_ylabel('runID')
        #ax5.set_xlabel('runID')

        #fig.savefig(pp, bbox_inches='tight', format='pdf')
        #fig3.savefig(pp, bbox_inches='tight', format='pdf')
        #fig4.savefig(pp, bbox_inches='tight', format='pdf')
        
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

        fig7 = plt.figure()
        ax7 = fig7.add_subplot(111)
        ax7.set_title('relative difference of the estimates for ohdu=%d'%ohdu_)
        ax7.set_xlabel('runID')
        ax7.set_ylabel(r'(N-N*)/N*')
        ax7.grid(True)

        #ax7.plot( runID2, (meanSiA - gainCu_[runID_==runID2]*SiPeak)/(gainCu_[runID_==runID2]*SiPeak), '.', ms = ms, label = 'Si' )
        #ax7.plot( runID2, (meanCuA - gainCu_[runID_==runID2]*CuPeak)/(gainCu_[runID_==runID2]*CuPeak),'.', ms = ms, label = 'Cu' )
        #ax7.plot( runID2, (meanCu2A - gainCu_[runID_==runID2]*CuPeak2)/(gainCu_[runID_==runID2]*CuPeak2), '.', ms = ms, label = 'Cu2' )
        #ax7.legend( fancybox=True, framealpha=0 )
        
        fig2.savefig(pp,  format='pdf')
        fig6.savefig(pp,  format='pdf')
        #fig7.savefig(pp, bbox_inches='tight',  format='pdf')

        #fig5.savefig(pp, bbox_inches='tight', format='pdf')
            
    
    pp.close()
    pp2.close()
    os.rename( outfile + '.tmp', outfile )
    os.rename( outfile2 + '.tmp', outfile2 )
    return

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == '--h':
            print "usage:"
            print "--h\tprints this message"
            print "--muon\tprocesses the muon counts"
            print "--muonPlot\tgenetares the muon count plots"
            print "--peak\tprocess the peaks info"
            print "--peakPlot\tgenerates the peaks plots"
            exit(0)
        if sys.argv[1] == "--muon":
            print "Processes the muon count"
            catalogs = sys.argv[-1]
            catalogs = glob.glob(catalogs)
            exit(0)
        if sys.argv[1] == "--muonPlot":
            print "Processes the muon count"
            exit(0)
        if sys.argv[1] == "--peak":
            print "Processes the muon count"
            catalogs = sys.argv[-1]
            catalogs = glob.glob(catalogs)
            if os.path.splitext( catalogs[0] )[1] == ".root":
                print "trees", root_numpy.list_trees(catalogs[0])
                print "branches", root_numpy.list_branches(catalogs[0], 'hitSumm')
                print "branches", root_numpy.list_branches(catalogs[0], 'config')
                applyToAllCu( catalogs )
            exit(0)
        if sys.argv[1] == "--peakPlot":
            if len(sys.argv) < 3: 
                print "no file given"
                exit(0)
            print "Processes the muon count"
            data = sys.argv[2]
            if os.path.splitext( data )[1] == ".dat":
                makePlotCu( data )
            exit(0)
        
        print sys.argtv[1]," is not a known option"
        exit(0)
    else:
        print "no args given"
