# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
from Timer import Timer
from glob import glob
from termcolor import colored
import astropy.io.fits
from Statistics import MAD
import Statistics as stats
from scipy.optimize import curve_fit
import argparse

from utils.utils import Namespace
from paths.paths import *
debug = False
    
def norm_p0(data):
    return [np.sum(data), np.median(data), MAD(data) ]

def norm_pdf(x, A, mu, sigma):
    return A*norm.pdf(x, mu, sigma)

def norm_bounds(data):
    bounds = ( 
        [ data.size/10, -np.inf, 0 ], 
        [ data.size, np.inf, MAD(data) ]
        )
    return bounds

def poisson_norm_p0(data):
    return [data.size, np.median(data), MAD(data), 1000, .5 ]

def poisson_norm_bounds(data):
    bounds = ( 
        [ 0, -np.inf, 0, 1000-1, 0.], 
        [data.size, np.inf, MAD(data), 1000+1, np.inf]
        )
    return bounds

def poisson_norm_header():
    return 'A, µ, σ, g, λ'

def poisson_norm_pdf(x, A, mu, sigma, gain, lamb):
    return A*stats.poisson_norm.pdf(x, mu, sigma, gain, lamb)

def hist_fit( data, func, p0, lims=[None,None], binsize=1, bounds=None, header=lambda: '' ):
    xmin = data.min() if lims[0] is None else lims[0]
    xmax = data.max() if lims[1] is None else lims[1]
    bins = np.arange(xmin, xmax, binsize)
    y = np.histogram(data, bins=bins)[0]
    x = .5*(bins[1:] + bins[:-1])
    if not bounds is None:
        popt = curve_fit( func, x, y, p0=p0, bounds=bounds )[0]
    else:
        popt = curve_fit( func, x, y, p0=p0 )[0]
    return popt, header()

def compute_each( HDU, args, print_header ):
    data = HDU.data.astype(float)
    ohdu = HDU.header['OHDU']
    data = data.T
    data[data>1e9] = 0
#     print( 'min,max', np.min(data), np.max(data) )
    if 'side' in args:
        side = int( args.side )
        data = data[side*data.shape[0]:(side+1)*data.shape[0], :]
        
    if 'trim' in args:
        if args.trim == 'trim8':
            data = data[8:,:]
    
    if 'Erange' in args:
        if len(args.Erange) == 1:
            args.Erange.append(None)
        data = data.flatten()
        data = data[ np.all( data > args.Erange[0], data < args.Erange[1]) ]
    
    if 'xrange' in args:
        if len(args.xrange) == 1:
            args.xrange.append(None)
        data = data[args.xrange[0]:args.xrange[1], :]

    if 'yrange' in args:
        if len(args.yrange) == 1:
            args.yrange.append(None)
        data = data[:, args.yrange[0]:args.yrange[1]]
    
    if args.axis == '*':
        data = data.flatten()
    elif args.axis == '1':
        data = data.T
        
    for ax, func in enumerate(args.funcs):
        data, header = eval( func )
    entry = '{:3d}'.format(ohdu) + ', ' + ', '.join(map(lambda _: '{:+.4e}'.format(_), data))
    if print_header:
        print( '#{ohdu}, {header}'.format(ohdu='ohdu', header=header) )
        print( '-'*len(entry) )
    print( entry )
    return
    
def compute( **args ):
    """
    Computes arbitrary functions on the images,
    so that this line is still used
    
    Arguments:
        input_file <input_file>    fits file to be computed
    
    Options:
        --ohdus <ohdus>    ohdus to be loaded
        --axis <a0>, <a1>, ...    order of the axes
        --funcs <f0>, <f1>, ...    functions to be applied on each axis
    """
    args = Namespace( **args )
    
    with Timer('compute'):
        if '*' in args.ohdus:
            ohdus = []
        else:
            ohdus = map( int, args.ohdus )
            
        if 'rawfits' in args.input_file or 'scnfits' in args.input_file:
            paths = eval( args.input_file )
        else:
            paths = glob( args.input_file )
        
        print_header = True
        for path in paths:
            print( colored('path', 'green'), path )
            title = path.split('/')[-1]
            
            listHDU = astropy.io.fits.open( path )
            imageHDU = None
            for i, HDU in enumerate(listHDU):
                if HDU.data is None:
                    continue
                if 'OHDU' in HDU.header:
                    if HDU.header['OHDU'] == '1':
                        continue
                    if HDU.header['OHDU'] in ohdus or ohdus == []:
                        compute_each( HDU, args, print_header )
                        print_header = False
    return

def make_argparse( subparser, func ):
    
    docstr = func.__doc__
    old_docstr = None
    while docstr != old_docstr:
        old_docstr = docstr
        docstr = docstr.replace('\n ', '\n')
    
    help_str = docstr.split('\n\n')[0].replace('\n', '')
    parser = subparser.add_parser( func.__name__, help = help_str )
    
    args_str = docstr.split('Arguments:\n')[1].split('\n\n')[0].split('\n')
    for arg in args_str:
#         key, descr = 
        print( arg.split( ' '*4 ) )
#         if len(key.split(' ')) == 1:
            
    opts_str = docstr.split('Options:\n')[1].split('\n\n')[0].split('\n')
    for opt in opts_str:
        print( opt )
    
    parser.print_help()
    

def add_compute_options( p ):
    p.add_argument(
        'input_file', 
        help='fits file input (example: "/share/storage2/connie/data/runs/*/runID_*_03326_*.fits.fz"' 
    )
    p.add_argument(
        '--ohdus',
        type=str,
        default='*',
        help='ohdus to be used',
        required = True
    )
    
    p.add_argument(
        '--out',
        type=str,
        default = argparse.SUPPRESS,
        help = 'output filename'
    )

    geom = p.add_argument_group('geometry options')
    geom.add_argument( '--xrange', nargs='+', type=eval, default = argparse.SUPPRESS, help = 'xmin xmax' )
    geom.add_argument( '--yrange', nargs='+', type=eval, default = argparse.SUPPRESS, help = 'ymin ymax' )
    geom.add_argument( '--Erange', nargs='+', type=eval, default = argparse.SUPPRESS, help = 'Emin Emax' )
    geom.add_argument( '--side', type=str, default = argparse.SUPPRESS, help = 'left or right amplifier' )
    geom.add_argument( '--section', type=str, default = argparse.SUPPRESS, help = 'data or bias' )
    geom.add_argument( '--trim', action='store_true', default = argparse.SUPPRESS, help = 'remove trim' )
    geom.add_argument( '--half', action='store_true', default = argparse.SUPPRESS, help = 'show half' )
    geom.add_argument( '--quarter', action='store_true', default = argparse.SUPPRESS, help = 'show quarter' )
    geom.add_argument( '--vmax', type=float, default = argparse.SUPPRESS, help = 'vmax' )

    corr = p.add_argument_group('correction options')
    corr.add_argument( '--remove', type=float, default=argparse.SUPPRESS, help = 'remove energies above' )
    corr.add_argument( '--global-bias', action='store_true', default=argparse.SUPPRESS, help = 'remove global bias' )
    corr.add_argument( '--correct-side', action='store_true', default=argparse.SUPPRESS, help = 'subtract right side' )
    corr.add_argument( '--smooth-lines', type=int, default=argparse.SUPPRESS, help = 'smooth lines' )
    corr.add_argument( '--sub-half', action='store_true', default=argparse.SUPPRESS, help = 'smooth lines' )
    corr.add_argument( '--overscan-subtraction', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'horizontal overscan subtraction' )
    corr.add_argument( '--vertical-overscan-subtraction', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'vertical overscan subtraction' )

    quant = p.add_argument_group('quantities options')
    quant.add_argument( '--axis', type=str, default=0, help = 'project on axis' )
    quant.add_argument( '--funcs', nargs='+', type=str, default=argparse.SUPPRESS, help='supress the max line at the projection plot' )

    p.set_defaults( func=compute )