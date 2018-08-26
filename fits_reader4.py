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

#folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029F/data_3321_to_3380/scn/merged/*'
folder = '/share/storage2/connie/data_analysis/processed02_data/runs/029*/data_*/scn/merged/*'

def analysis( cflag=True ):
    plot=True
    result = []
    for f in sorted(glob.glob(folder)):
        #print(f)
        #print(os.path.basename(f))
        run, runID = re.search(r'runID_([0-9]+)_([0-9]+)', os.path.basename(f)).groups()
        print('runID', run, runID )
        data = astropy.io.fits.open( f )
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
                'ovstrim': ohdu2.data[ downmargin:upmargin, ovscut+margin:-margin],
                'ccdstrip': ohdu2.data[ downmargin:upmargin, leftmargin:leftmargin+130],
                }
            }

        #print('shapes', data['data']['ovstrim'].shape, data['data']['ccdstrip'].shape )
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
                for n,key in enumerate(['ovstrim', 'ccdstrip']): #data['hist'].keys()):
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
            
            fig.savefig('ccd.png')
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
        np.savetxt('ccd.csv', result, header='run, runID, tempmin, tempmax, ovs_std, ovs_stda, ovs_chi2, ccd_std, ccd_stda, ccd_chi2' + ', nse_std, nse_amp, nse_chi2' if cflag else '', fmt='%s',delimiter=', ')

        #break
            
    #print result


def darkcurrent():
    data = np.loadtxt('ccd.csv', delimiter=',')
    data = data[data[:,1].argsort()]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ax2.step( data[:,1], data[:,3], 'y--', label='temp', where='mid' )
    ax.step(data[:,1], data[:,4], label='ovs', where = 'mid')
    ax.step(data[:,1], data[:,7], label='ccd', where='mid')
    ax.step(data[:,1], np.sqrt(data[:,7]**2 - data[:,4]**2), label='dc', where='mid')
    ax.step(data[:,1], np.sqrt(data[:,8]**2 - data[:,5]**2), label='dca', where='mid')
    #ax2.errorbar( data[:,1], .5*(data[:,2]+data[:,3]), xerr=.5, yerr = .5*(data[:,3]-data[:,2]), label='temp', linestyle='None' )
    ax.grid(True)
    ax.legend()
    ax2.legend()
    fig.savefig('dc.png')
    
if __name__ == "__main__":
    #analysis(cflag=False)
    print('dc')
    darkcurrent()
