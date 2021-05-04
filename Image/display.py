# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import astropy.io.fits
from numpy import *

import matplotlib
matplotlib.use('gtk3agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable

# from matplotlib import pylab as plt


def add_display_options( p ):
    p.add_argument('input_file', help = 'fits file input (example: "/share/storage2/connie/data/runs/*/runID_*_03326_*.fits.fz"' )
    p.add_argument(
        '--ohdu',
        type = int,
        default = 2,
        help = 'ohdu to be displayed',
        required = True
    )
    p.add_argument(
        '--plot',
        nargs = '+',
        type = str,
        default = ['image'],
        # choices = ['proj', 'spectrum', 'image'],
        required = True,
        help = 'plot type'
    )
    p.add_argument(
        '--png',
        action = 'store_true',
        default = argparse.SUPPRESS,
        help = 'output to png file'
    )

    geom = p.add_argument_group('geometry options')
    geom.add_argument( '--E-range', nargs=2, type=eval, default = [-inf, inf], help = 'Emin Emax' )
    geom.add_argument( '--E-span', nargs=1, type=eval, default = argparse.SUPPRESS, help = 'mean+-E_span' )
    geom.add_argument( '--x-range', nargs=2, type=eval, default = [None, None], help = 'xmin xmax' )
    geom.add_argument( '--y-range', nargs=2, type=eval, default = [None, None], help = 'ymin ymax' )
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

    proj = p.add_argument_group('projection options')
    proj.add_argument( '--axis', type=int, default = 0, help = 'project on axis' )
    proj.add_argument( '--dev', action='store_true', default = argparse.SUPPRESS, help = 'supress the max line at the projection plot' )
    proj.add_argument( '--no-max', action='store_true', default = argparse.SUPPRESS, help = 'supress the max line at the projection plot' )
    proj.add_argument( '--no-min', action='store_true', default = argparse.SUPPRESS, help = 'supress the min line at the projection plot' )
    proj.add_argument( '--no-mean', action='store_true', default = argparse.SUPPRESS, help = 'supress the mean line at the projection plot' )
    proj.add_argument( '--no-center', action='store_true', default = argparse.SUPPRESS, help = 'supress the mean line at the projection plot' )
    proj.add_argument( '--no-median', action='store_true', default = argparse.SUPPRESS, help = 'supress the mean line at the projection plot' )
    proj.add_argument( '--smooth', type=int, default = argparse.SUPPRESS, help = 'smoothening length' )

    proj.add_argument( '--fit', nargs='append', type=str, default = argparse.SUPPRESS, help = 'fit command' )

    spec = p.add_argument_group('spectrum options')
    spec.add_argument( '--binsize', type=float, default=1, help = 'binsize' )

    p.set_defaults( func=display )




def display( args ):
    args = Namespace(args)
    with Timer('plot'):
        path = glob( args.input_file )
        if len(path) > 1:
            print( 'found {}'.format(path) )
            print( 'using the first path' )
        path = path[0]
        print( colored('path', 'green'), path )
        title = path.split('/')[-1]
        if not os.path.exists(path):
            print( 'file not found:', path )
            exit(0)
        listHDU = astropy.io.fits.open( path )
        imageHDU = None
        for i, HDU in enumerate(listHDU):
            if 'OHDU' in HDU.header:
                if HDU.header['OHDU'] == args.ohdu:
                    imageHDU = HDU
                    print( colored('ohdu', 'green'), args.ohdu )
                    break
            else:
                if i+1 == args.ohdu: imageHDU = HDU
        if imageHDU is None:
            print( 'ohdu {} was not found in {}'.format( args.ohdu, path ) )
            exit(0)

        data = imageHDU.data.astype(float)
        data = data[::-1,:]
        data[data>1e9] = nan

        height, width = data.shape
        print( colored('height', 'green'), height )
        print( colored('width', 'green'), width )
        number_of_amplifiers = width // constants.ccd_width
        print( colored('amplifiers', 'green'), number_of_amplifiers )
        if number_of_amplifiers > 0:
            side_width = width//number_of_amplifiers
        else:
            side_width = width
        print( colored('side width', 'green'), side_width )

        bias_width = int( ceil( (side_width - constants.ccd_width)/150. ) )*150
        print( colored('bias width', 'green'), bias_width )
        trim = constants.ccd_width + bias_width - side_width
        print( 'trim', trim )
        #bias_width = ( side_width + trim ) % constants.ccd_width

        height_trim = 1 #1

        #data.left = data[None:-height_trim, None:side_width]
        #data.right = data[None:-height_trim, side_width:None][:,::-1]

        if 'overscan_subtraction' in args:
            os_range = args.overscan_subtraction
            print( 'os_range', os_range, os_range[0] )
            correction = median( data[ None:None, os_range[0]:os_range[1] ], axis=1 )[:,None]
            data -= correction

        if 'vertical_overscan_subtraction' in args:
            os_range = args.vertical_overscan_subtraction
            print( 'os_range', os_range, os_range[0] )
            correction = median( data[ os_range[0]:os_range[1], None:None ], axis=0 )[None,:]
            data -= correction

        if 'side' in args:
            if args.side == 'left' or args.side == '0':
                print( colored('side', 'green'), 'left' )
                if 'correct_side' in args:
                    print( colored('subtract', 'green'), 'right' )
                    data = data[None:-height_trim, None:side_width] - data[None:-height_trim, side_width:None][:,::-1]
                else:
                    data = data[None:-height_trim, None:side_width]
            elif args.side == 'right' or args.side == '1':
                print( colored('side', 'green'), 'right' )
                data = data[None:-height_trim, side_width:None][:,::-1]
                if 'sub_half' in args:
                    print( colored('subtract line means', 'green'), 'true' )
                    correction = mean( data[None:None, constants.ccd_width/2:None], axis=1 )
                    if len(args.sub_half) > 0:
                        l = int(args.sub_half[0])
                        #correction =
                    data = data - correction[:,None]
            else:
                print( colored('side', 'green'), 'all' )

        if 'side' in args or number_of_amplifiers == 1:
            if 'half' in args:
                print( colored('half', 'green'), 'true' )
                data = data[None:None, -constants.ccd_width/2:None]
            elif 'quarter' in args:
                print( colored('quarter', 'green'), 'true' )
                data = data[None:None, -constants.ccd_width/4:None]

        if 'smooth_lines' in args:
            print( colored('smooth lines', 'green'), args.smooth_lines )
            data = scipy.ndimage.filters.uniform_filter1d( data, axis=1, size=args.smooth_lines, mode='mirror' )

        if 'trim' in args:
            print( colored('trim', 'green'), '8' )
            data = data[None:None, 8:]

        if 'global_bias' in args:
            print( colored('subtract', 'green'), 'mean(OS)' )
            data = data - nanmean( data[None:None, -bias_width:None] )

        if 'section' in args:
            if args.section == 'bias' or args.section == 'os':
                print( colored('section', 'green'), 'overscan' )
                data = data[None:None, -bias_width:None]
            if args.section == 'data' or args.section == 'ac':
                print( colored('section', 'green'), 'active' )
                data = data[None:None, None:-bias_width+5]
                if 'half' in args:
                    print( colored('half', 'green'), 'true' )
                    data = data[None:None, data.shape[1]/2:None]

        if 'remove' in args:
            print( colored('remove above', 'green'), args.remove )
            data[ data > args.remove ] = nan


        if 'x_range' in args:
            data = data[:, args.x_range[0]:args.x_range[1]]
            print( colored('x-range', 'green'), args.x_range )
        if 'y_range' in args:
            data = data[args.y_range[0]:args.y_range[1], :]
            print( colored('y-range', 'green'), args.y_range )
        if 'E_span' in args:
            E = nanmedian(data)
            mask = logical_or( data<E-args.E_span, data>E+args.E_span )
            data[mask] = nan
            print( colored('E-span', 'green'), args.E_span, nanmin(data), nanmax(data) )
        elif 'E_range' in args:
            data[data<=args.E_range[0]] = nan
            data[data>args.E_range[1]] = nan
            print( colored('E-range', 'green'), args.E_range, nanmin(data), nanmax(data) )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if args.plot[0] == 'proj':
            print( colored('plot', 'green'), 'projection' )
            axis = args.plot[1]
            if axis in ['lines', '1', 'y']:
                axis = 1
            elif axis in ['rows', '0', 'x']:
                axis = 0
            else:
                print( 'axis not expected' )
                exit()
            x = arange( 0, data.shape[axis-1] ).astype(float)
            if not 'dev' in args:
                mins = nanmin( data, axis=axis )
                if not 'no_mean' in args:
                    means = nanmean( data, axis=axis )
                    print( means.shape )
                    ax.plot( means, '.', label='mean $\mu={:.4f}$'.format(nanmean(means)) )
                if not 'no_median' in args:
                    medians = nanmedian( data, axis=axis )
                    ax.plot( medians, '.', label='median $\mu={:.4f}$'.format(nanmean(medians)) )

                if 'smooth' in args:
                    #ax.plot( scipy.ndimage.filters.uniform_filter1d( centers, size=args.smooth, mode='nearest' ), '.', label='center{}'.format(args.smooth) )
                    ax.plot( scipy.ndimage.filters.uniform_filter1d( medians, size=args.smooth, mode='nearest' ), '.', label='median{}'.format(args.smooth) )
                    if not 'no_mean' in args:
                        ax.plot( scipy.ndimage.filters.uniform_filter1d( means, size=args.smooth, mode='nearest' ), '.', label='mean{}'.format(args.smooth) )
                if not 'no_max' in args:
                    maxs = nanmax( data, axis=axis )
                    ax.plot( maxs, '.', label='max' )
                if not 'no_min' in args:
                    mins = nanmin( data, axis=axis )
                    ax.plot( mins, '.', label='min' )
#                 if 'fit' in args:
#                     for fit in args.fit:
#                         exec F(fit).str()
            else:
                mins = nanmin( data, axis=axis )
                if not 'no_mean' in args:
                    stds = nanstd( data, axis=axis )
                    ax.plot( stds, '.', label='std $\mu={:.4f}$'.format(nanmean(stds)) )
                mads = MAD( data, axis=axis )
                ax.plot( mads, '.', label='MAD $\mu={:.4f}$'.format(nanmean(mads)) )
                #fwhm = stats.FWHM( data, axis=axis, f=.001 )[0]
                #ax.plot( fwhm, '.', label='FWHM $\mu={:.4f}$'.format(nanmean(fwhm)) )
                sigmas = mle_norm( data, axis=axis )[1]
                ax.plot( sigmas, '.', label='$\sigma$[mle] $\mu={:.4f}$'.format(nanmean(sigmas)) )

        elif args.plot[0] == 'spectrum':
            bins = arange( nanmin(data), nanmax(data), args.binsize )
            y, x, dx = stats.make_histogram( data.flatten(), bins )
            ax.step( x, y, where='mid', label = '' )

        elif args.plot[0] == 'image' or args.plot[0] == 'matrix':
            cmap = matplotlib.cm.seismic
            cmap.set_bad(color='black')
            if 'log' in args.plot[1:]:
                im = log(data - nanmin(data) + 1)
                cmap.set_bad('black',0)
            else:
                im = data
                cmap.set_bad('black',0)
            try:
                vmean = nanmedian(im)
                vmax = nanmax(abs(im - vmean))
                if 'vmax' in args:
                    vmax = args.vmax

                print( 'vmax', vmax )
                if args.plot[0] == 'image':
                    func = ax.imshow
                elif args.plot[0] == 'matrix':
                    func = ax.matshow
                obj = func( im, cmap=cmap, origin='lower', vmin=vmean-vmax, vmax=vmean+vmax, interpolation='none' )
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.1)
                fig.colorbar( obj, cax )
                title += '\n mean {:.4}'.format( nanmean( im ) )
                title += '\n std {:.4}'.format( nanstd( im ) )
            except:
                print( 'im.shape', im.shape, data.shape )
                exit()
        else:
            print( 'plot option "{}" not implemented'.format( args.plot ) )
            exit(0)

        ax.set_title( title, fontsize=8 )
        ax.grid()
        ax.legend()
    if 'png' in args:
        outfile = args.input_file + '.' + '_'.join(args.plot) + '.png'
        print( 'outfile', outfile )
        fig.savefig( outfile, bbox_inches='tight', pad_inches=0 )
    else:
        plt.show()
    return