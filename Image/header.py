# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import astropy.io.fits

from termcolor import colored

from utils.utils import Namespace
from paths.paths import *
debug = False

def add_header_options( p ):
    p.add_argument( 'input_file' )
    p.add_argument( '--indices', nargs='*', default = None, type=int, help='indexes to be shown' )
    p.add_argument( '--include', nargs='*', default = '', help='fields to be shown' )
    p.add_argument( '--exclude', nargs='*', default = '', help='fields not to be shown' )
    p.set_defaults( func=header )
    return

def header( **args ):
    args = Namespace( **args )
    
    if 'rawfits' in args.input_file or 'scnfits' in args.input_file:
        args.input_file = eval( args.input_file )
    else:
        args.input_file = glob( args.input_file )
    
    for file in args.input_file:
        listHDU = astropy.io.fits.open( file )
        for i in range(len(listHDU)) if args.indices is None else args.indices:
            print( colored( 'HDU %s' % i, 'green', attrs=['bold'] ) )
            for key, value in listHDU[i].header.items():
                if key in args.exclude: continue
                if key in args.include or len(args.include) is 0:
                    print( colored( '%s' % key, 'green' ), '%s' % value, end=' ' )
            print()
    return
