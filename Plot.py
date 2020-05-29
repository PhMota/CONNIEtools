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

def generate_plot( input_files, diagonal_line = False, xcolumn = 0, ycolumns = [1], xlabel = None, ylabel = None, marker = '.', output_file = 'none', fig_scale = [1,1], labels = [''], legend_loc = 'upper right', yerrcolumns = None, **args ):
    if type(input_files) != list:
        input_files = [input_files]
    xdata = []
    ydata = []
    yerrdata = []
    if yerrcolumns is None: yerrcolumns = [None]*len(ycolumns)
    
    for input_file in input_files:
        print( 'input_file', input_file )
        if '.csv' in input_file:
            skip_lines = -1
            for line in open( input_file ):
                if '#' in line:
                    skip_lines += 1
                else:
                    break
            data = np.genfromtxt( input_file, delimiter = ',', names = True, skip_header = skip_lines, deletechars= '#' )
            print( 'available columns in file', data.dtype.names )
            xdata.append( data[ xcolumn ] )
            ydata.append( [ data[column] for column in ycolumns ] )
            yerrdata.append( [ None if column is None else data[column] for column in yerrcolumns ] )
            
    fig = plt.figure()
    size = fig.get_size_inches()
    if fig_scale:
        fig.set_size_inches( fig_scale[0]*size[0], fig_scale[1]*size[1] )
    ax = fig.add_subplot(111)
    if xlabel is None: xlabel = xcolumn
    ax.set_xlabel( xlabel )
    if ylabel is None: ylabel = ''
    ax.set_ylabel( ylabel )
    
    colors = ['b', 'r', 'g', 'm', 'c', 'y']
    count = 0
    for i, (xdatum, ydatum, yerrdatum, label) in enumerate(zip(xdata, ydata, yerrdata, labels)):
        for j, (y, yerr, column) in enumerate(zip(ydatum, yerrdatum, ycolumns)):
            print( 'yerr', yerr )
            if yerr is None: yerr = 0
            ax.errorbar( xdatum, y, yerr = yerr, label = label, marker = marker, color = colors[count], linestyle = '' )

        #if diagonal_line:
            #lims = ax.get_xlim()
            #lenght = lims[1]-lims[0]
            #kwargs = {'edgecolor':colors[count], 'facecolor':colors[count], 'head_width': .02*lenght, 'length_includes_head': True}
            #offset = .1*lenght
            #for y in ydatum:
                #xarrows = xdatum[ y < lims[0] ]
                #for xarrow in xarrows: ax.arrow( xarrow, lims[0] + offset, 0, -offset+.1, **kwargs )
                #xarrows = xdatum[ y > lims[1] ]
                #for xarrow in xarrows: ax.arrow( xarrow, lims[1] - offset, 0, +offset-.1, **kwargs )
        count += 1
    if diagonal_line:
        lims = ax.get_xlim()
        print( 'adding diagonal line', lims )
        ax.plot( [lims[0], lims[1]], [lims[0],lims[1]], 'k-', zorder = 0, alpha = .5 )
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    legend = ax.legend( fancybox=True, framealpha=0, bbox_to_anchor=(1.,1.), loc=legend_loc )
    
    if output_file != 'none':
        print( 'save to file ' + output_file )
        fig.savefig( output_file, bbox_extra_artists = (legend,), bbox_inches='tight', pad_inches=0 )
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
    
    parser.add_argument('input_files', type=eval, help = 'list of files to be plotted' )
    parser.add_argument('--output-file', type=str, default = 'none', help = 'set to none for X11 display' )
    parser.add_argument('--xlabel', type=str, default = None, help = 'xlabel' )
    parser.add_argument('--ylabel', type=str, default = None, help = 'ylabel' )
    parser.add_argument('--xcolumn', type=str, default = 0, help = 'column of the file associated with the x axis' )
    parser.add_argument('--ycolumns', type=eval, default = "[1]", help = 'list of columns of the file associated with the y axis' )
    parser.add_argument('--yerrcolumns', type=eval, default = None, help = 'list of columns of the file associated with the yerr axis' )
    parser.add_argument('--diagonal-line', type=bool, default = False, help = 'add diagonal line to the plot' )
    parser.add_argument('--marker', type=str, default = '.', help = 'marker type' )
    parser.add_argument('--fig-scale', type=eval, default = '[1,1]', help = 'scale the figure size' )
    parser.add_argument('--labels', type=eval, default = "['label']", help = 'list of labels' )

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )

    generate_plot( **vars(args) )
    exit(0)
