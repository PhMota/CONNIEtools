# coding: utf-8
from __future__ import print_function
import numpy as np
import astropy.io.fits
import scipy.stats
import scipy.ndimage
from glob import glob
import re
import gi

import matplotlib
matplotlib.use('gtk3agg')
from Timer import Timer

def generate_plot( args ):
    data = None
    if '.csv' in args['input_file']:
        skip_lines = -1
        for file in open( args['input_file'] ):
            line = file.readline()
            if '#' in line:
                skip_lines += 1
            else:
                break
        data = np.genfromtxt( args['input_file'], delimiter = ',', names = True, skip_header = skip_lines )
        data = [ data[column] for column in args['columns'] ]
        #skip_header = 1, dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] )
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter( data[0], data[1] )
    
    if 'output_file' in args:
        fig.savefig(args['output_file')
    else:
        fig.show()
    return

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description = 'plot file',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        )
    def tuple_of_int( x ):
        return map( int, eval(x.replace('\"', '')) )
    
    parser.add_argument('input_file', type=str, help = 'file to be plotted' )
    parser.add_argument('--output-file', type=int, default = 'none', help = 'set to none for X11 display' )
    parser.add_argument('--xlabel', type=int, default = r'$x$', help = 'xlabel' )
    parser.add_argument('--ylabel', type=int, default = r'$y$', help = 'ylabel' )
    parser.add_argument('--xcolumn', type=int, default = 0, help = 'column of the file associated with the x axis' )
    parser.add_argument('--ycolumns', type=int, default = [1], help = 'list of columns of the file associated with the y axis' )

    parser.add_argument('--horizontal-overscan', type=int, default = '150', help = 'size of the horizontal overscan in pixels' )
    parser.add_argument('--vertical-overscan', type=int, default = '90', help = 'size of the vertical overscan in pixels' )
    parser.add_argument('--xyshape', type=tuple_of_int, default = '\"[4000,4000]\"', help = 'shape of the image as 2d pixels' )
    parser.add_argument('--rebin', type=tuple_of_int, default = '\"[1,1]\"', help = '2d rebinning strides' )
    parser.add_argument('--charge-range', type=tuple_of_int, default = '\"[5,200]\"', help = 'range into which to randomly generate charges' )
    parser.add_argument('--depth-range', type=tuple_of_int, default = '\"[0,670]\"', help = 'range into which to randomly generate depths' )
    parser.add_argument('--charge-gain', type=eval, default = '7.25', help = 'factor to convert charges into ADU' )
    parser.add_argument('--readout-noise', type=eval, default = '0', help = 'sigma of the normal noise distribution in ADU' )
    parser.add_argument('--dark-current', type=eval, default = '0', help = 'lambda of Poisson distribution dimensionless' )
    parser.add_argument('--simulation-output', type=str, default = 'simulation.csv', help = 'csv file with the events generation data' )
    parser.add_argument('--image-fits-output', type=str, default = 'simulation.fits', help = 'set to "none" not to generate a fits output' )
    parser.add_argument('--image-pdf-output', type=str, default = 'simulation.pdf', help = 'set to "none" not to generate a pdf output' )
    parser.add_argument('--image-energy-spectrum', type=str, default = 'image_energy_spectrum.png', help = 'set to "none" not to plot image energy spectrum' )
    parser.add_argument('--diffusion-function',
                        type=str, 
                        default = default_diffusion_function,
                        help = 'function to map z-depth into sigma' 
                        )
    parser.add_argument('--charge-efficiency-function',
                        type=str,
                        default = default_charge_efficiency_function,
                        help = 'function to map z-depth into sigma' 
                        )
    parser.add_argument('--vertical-modulation-function',
                        type=str, 
                        default = default_vertical_modulation_function,
                        help = 'function to modulate the vertical axis' 
                        )    
    parser.add_argument('--horizontal-modulation-function',
                        type=str, 
                        default = default_horizontal_modulation_function,
                        help = 'function to modulate the horizontal axis' 
                        )
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )

    generate_plot( vars(args) )
