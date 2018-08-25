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

import glob
import os
import re
import copy
import root_numpy

folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/scn/merged/*'

plot=True
result = []
for f in glob.glob(folder):
    #print(f)
    print(os.path.basename(f))
    run, runID = re.search(r'runID_([0-9]+)_([0-9]+)', os.path.basename(f)).groups()
    print('runID', run, runID )
    data = astropy.io.fits.open( f )
    ohdu2 = None
    for datum in data:
        if datum.header['OHDU'] == 2: 
            ohdu2 = datum
            break
    print(f, os.path.basename(f))
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
    print('shape', ohdu2.data.shape)
    data = {'data': {
            #'ccd': ohdu2.data[:,:ovscut],
            'ccdtrim': ohdu2.data[ downmargin:upmargin, leftmargin:rightmargin],
            #'ovs': ohdu2.data[:,ovscut:],
            'ovstrim': ohdu2.data[ downmargin:upmargin, ovscut+margin:-margin]
            }
        }

    activeccd = np.array(ohdu2.data[ downmargin:upmargin, leftmargin:rightmargin])
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
    
    #for x,y,e in zip(readcatalog['xPix'], readcatalog['yPix'], readcatalog['ePix']):
        #mask = np.all( [x>leftmargin,x<rightmargin,y>downmargin,y<upmargin], axis=0)
        #if len(y[mask])>0:
            #activeccd[ y[mask]-downmargin, x[mask]-leftmargin ] = 0
            
    #data['data']['hits'] = np.array( [ ie for x,y,e in zip(readcatalog['xPix'], readcatalog['yPix'], readcatalog['ePix']) for ie in list(e) ] )
    data['data']['hitstrim'] = np.array( [ ie for x,y,e in zip(readcatalog['xPix'], readcatalog['yPix'], readcatalog['ePix']) for ie in list(e[np.all( [x>leftmargin,x<rightmargin,y>downmargin,y<upmargin], axis=0)]) ] )
        
    #data['data']['hits4'] = activeccd.flatten()
    
    dx = 1
    bins = np.r_[min(ohdu2.data.flatten()):max(ohdu2.data.flatten()):dx]
    xbins = (bins[1:]+bins[:-1])*.5
    #print('dx',xbins[1]-xbins[0])
    fabove = 1e-5
    
    fitfunc = lambda x,b,c: c*scipy.stats.norm.pdf(x,0,b)
    
    data['hist'] = {}
    data['pp'] = {}
    data['chisquare'] = {}
    data['ppleft'] = {}
    data['chisquareleft'] = {}
    for key in data['data'].keys():
        data['hist'][key], tmp = np.histogram( data['data'][key].flatten(), bins = bins )
    data['hist']['ccdtrim-hitstrim'] = data['hist']['ccdtrim'] - data['hist']['hitstrim']

    for key in data['hist'].keys():
        if 1:# key == 'ccd' or key == 'ovs':
            above = data['hist'][key] > fabove*max(data['hist'][key])
            data['pp'][key] = scipy.optimize.curve_fit( fitfunc, xbins[above], data['hist'][key][above] )[0]
            
            data['chisquare'][key] = scipy.stats.chisquare( fitfunc( xbins[above], *data['pp'][key]), f_exp = data['hist'][key][above] )[0]/len(data['hist'][key][above])
            
            left = xbins[above] < data['pp'][key][0] #std
            data['ppleft'][key] = scipy.optimize.curve_fit( fitfunc, xbins[above][left], data['hist'][key][above][left] )[0]
            data['chisquareleft'][key] = scipy.stats.chisquare( fitfunc( xbins[above][left], *data['ppleft'][key]), f_exp = data['hist'][key][above][left] )[0]/len(data['hist'][key][above][left] )
    
    entry = [ run, runID, 
             data['ppleft']['ovstrim'][0], data['ppleft']['ovstrim'][1], data['chisquareleft']['ovstrim'], 
             data['ppleft']['ccdtrim'][0], data['ppleft']['ccdtrim'][1], data['chisquareleft']['ccdtrim'],
             data['ppleft']['ccdtrim-hitstrim'][0], data['ppleft']['ccdtrim-hitstrim'][1], data['chisquareleft']['ccdtrim-hitstrim'],
             ]
    print(entry)
    result += [entry]
    np.savetxt('ccd.csv', result, header='run, runID, ovs_std, ovs_amp, ovs_chi2, ccd_std, ccd_amp, ccd_chi2, dc_std, dc_amp, dc_chi2',fmt='%s',delimiter=', ')

    if plot:
        fig = plt.figure()
        fig.suptitle(
            'ohdu2 runID'+runID+'ccd[%s:%s,%s:%s]'%(downmargin,upmargin,leftmargin,rightmargin) + ' ovs[%s:%s,%s:%s]'%(downmargin,upmargin,ovscut+margin,ohdu2.data.shape[1]-margin)+'\ndx=%sadu, fit>%smax'%(dx,fabove) )
        ax = fig.add_subplot(211)
        axzoom = fig.add_subplot(212)
        for axis in [ax,axzoom]:
            for n,key in enumerate(data['hist'].keys()):
                axis.step( xbins, data['hist'][key], 'C%d-'%n, where='mid', label = key )

            #axis.step( xbins, data['hist']['ovs'], 'm-', where='mid', label = 'ovs' )
            #axis.step( xbins, data['hist']['hits4'], 'g-', where='mid', label = 'hits4' )
            ##axis.step( xbins, data['hist']['hits3'], 'y-', where='mid', label = 'hits3' )  
            #axis.step( xbins, data['hist']['hits'], 'r-', where='mid', label = 'hits' )
            #axis.step( xbins, data['hist']['hitstrim'], 'b-', where='mid', label = 'hitstrim' )
            #axis.step( xbins, data['hist']['ccd'], 'c-', where='mid', label = 'ccd' )
            #axis.plot( xbins, fitfunc(xbins,*data['ppleft']['ovs']), 'r--')
            #axis.plot( xbins, fitfunc(xbins,*data['ppleft']['ccd']), 'b--')
            axis.set_yscale('log')
        ax.legend()
        #axzoom.text(0.95, 0.95, 'ovs(%.2f,%.2f)%.2e,%.2e\ndata(%.2f,%.2f)%.2e,%.2e'%(
                #data['pp']['ovs'][0],
                #data['ppleft']['ovs'][0],
                #data['chisquare']['ovs'],
                #data['chisquareleft']['ovs'],
                #data['pp']['ccd'][0],
                #data['ppleft']['ccd'][0],
                #data['chisquare']['ccd'],
                #data['chisquareleft']['ccd'],
                #), transform=ax.transAxes, verticalalignment='top', horizontalalignment='right')
        top = max(data['hist']['ccdtrim'])
        #ax.set_ylim([1e-8,1e-1])
        ax.set_xlim([-100,100])
        axzoom.set_ylim([ fabove*top, 1.1*top ])
        #axzoom.set_xlim([min(xbins[above]),-min(xbins[above])])
        axzoom.set_xlim([-100,100])
        
        fig.savefig('ccd.png')
        plot = False
    #break
        
print result

    
