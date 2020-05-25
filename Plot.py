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
import matplotlib.pylab as plt
from Timer import Timer

def generate_plot( input_file, diagonal_line = False, xcolumn = 0, ycolumns = [1], xlabel = r'$x$', ylabel = r'$y$', marker = '.', output_file = 'none', fig_scale = [1,1], **args ):
    data = None
    if '.csv' in input_file:
        skip_lines = -1
        for line in open( input_file ):
            if '#' in line:
                skip_lines += 1
            else:
                break
        data = np.genfromtxt( input_file, delimiter = ',', names = True, skip_header = skip_lines, deletechars= '#' )
        print( 'available columns in file', data.dtype.names )
        xdatum = data[ xcolumn ]
        ydata = [ data[column] for column in ycolumns ]
    
    fig = plt.figure()
    size = fig.get_size_inches()
    if fig_scale:
        fig.set_size_inches( fig_scale[0]*size[0], fig_scale[1]*size[1] )
    ax = fig.add_subplot(111)
    ax.set_xlabel( xlabel )
    ax.set_ylabel( ylabel )
    
        
    for column in ycolumns:
        ax.scatter( xdatum, data[column], label = column, marker = marker )

    if diagonal_line:
        lims = ax.get_xlim()
        print( 'adding diagonal line', lims )
        ax.plot( [lims[0], lims[1]], [lims[0],lims[1]], 'k-', zorder = 0, alpha = .5 )
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        lenght = lims[1]-lims[0]
        kwargs = {'edgecolor':'b', 'head_width': .02*lenght, 'length_includes_head': True}
        for column in ycolumns:
            y = data[column]
            xarrows = xdatum[ y < lims[0] ]
            offset = .1*lenght
            for xarrow in xarrows: ax.arrow( xarrow, lims[0] + offset, 0, -offset+.1, **kwargs )
            xarrows = xdatum[ y > lims[1] ]
            for xarrow in xarrows: ax.arrow( xarrow, lims[1] - offset, 0, +offset-.1, **kwargs )
            
    
    legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc='upper right' )
    
    if output_file != 'none':
        print( 'save to file ' + output_file )
        fig.savefig( output_file, bbox_extra_artists = (legend,), bbox_inches='tight' )
    else:
        plt.show()
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
    parser.add_argument('--output-file', type=str, default = 'none', help = 'set to none for X11 display' )
    parser.add_argument('--xlabel', type=str, default = r'$x$', help = 'xlabel' )
    parser.add_argument('--ylabel', type=str, default = r'$y$', help = 'ylabel' )
    parser.add_argument('--xcolumn', type=str, default = 0, help = 'column of the file associated with the x axis' )
    parser.add_argument('--ycolumns', type=eval, default = "[1]", help = 'list of columns of the file associated with the y axis' )
    parser.add_argument('--diagonal-line', type=bool, default = False, help = 'add diagonal line to the plot' )
    parser.add_argument('--marker', type=str, default = '.', help = 'marker type' )
    parser.add_argument('--fig-scale', type=eval, default = '[1,1]', help = 'scale the figure size' )

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )

    generate_plot( **vars(args) )
    exit(0)
