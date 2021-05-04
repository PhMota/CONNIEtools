# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import astropy.io.fits


def add_analyse_options( p ):
    p.add_argument('name', help = 'fits file input (example: runID3326)' )
    p.add_argument('input_file', help = 'fits file input (example: "/share/storage2/connie/data/runs/*/runID_*_03326_*.fits.fz"' )

    p.add_argument( '--ohdu', nargs='*', type=int, default = None, help = 'ohdus to be analysed' )
    p.add_argument( '--exclude', nargs='*', type=int, default = None, help = 'ohdus not to be analysed' )

    p.add_argument('--params-mode', type=str, default = 'median', help = 'modes for parameter estimation' )
    p.add_argument( '--fix-vbias', action="store_true", default=argparse.SUPPRESS, help = 'mantain the vbais shift' )

    p.add_argument('--remove-hits', nargs=2, type=float, default=argparse.SUPPRESS, help = 'remove hits above ADU with border (example:60 3)' )

    p.add_argument('--find-hits', action="store_true", default=argparse.SUPPRESS, help = 'hits above ADU border' )

    p.add_argument( '--plot-part', action="store_true", default=argparse.SUPPRESS, help = 'plot parts' )
    p.add_argument( '--plot-sections', action="store_true", default=argparse.SUPPRESS, help = 'plot sections' )
    p.add_argument( '--plot-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot spectrum' )
    p.add_argument( '--plot-convolution-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot convolution spectrum' )
    p.add_argument( '--convolution-function', type=str, default=argparse.SUPPRESS, help = 'plot convolution function' )

    p.add_argument( '--plot-block-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot block spectrum' )
    p.add_argument( '--block-function', type=str, default=argparse.SUPPRESS, help = 'function applied to each block' )
    p.add_argument( '--block-function-mean', action="store_true", default=argparse.SUPPRESS, help = 'apply mean to each block' )
    p.add_argument( '--block-function-std', action="store_true", default=argparse.SUPPRESS, help = 'apply std to each block' )

    p.add_argument( '--plot-convolutionfft-spectrum', action="store_true", default=argparse.SUPPRESS, help = 'plot convolution spectrum' )

    #p.add_argument('--extract', nargs=2, default=argparse.SUPPRESS, help = 'set to False to skip the extraction' )

    p.set_defaults( func=analyse )
    return

class Table:
    def __init__( self ):
        self.entries = []
    def list(self):
        return self.entries
    def append( self, entry ):
        self.entries.append( entry )
    def set_columns( self, cols ):
        self.cols = cols
    def get_entry( self, i ):
        return self.entries[i]
    def __str__( self ):
        ret = ', '.join(self.cols)
        for entry in self.entries:
            ret += '\n' + ', '.join( [
                ( {int:'{}', float:'{:.4f}', str:'{}', Section:'{:.4f}'}[type(entry[col])] ).format( entry[col] )
                for col in self.cols ] )
        return ret

def analyse( args ):
    if 'partlistHDU' in args:
        partlistHDU = [ args.partlistHDU ]
    else:
        paths = glob( args.input_file )
        partlistHDU = [ astropy.io.fits.open( path ) for path in paths[::-1] ]

    if not os.path.exists(args.name):
        os.mkdir(args.name)
    args.name = '{0}/'.format(args.name)

    table = Table()
    table.set_columns( ['ohdu', 'type', 'bias_sigma', 'g_lambda', 'g', 'lambda' ] )
    parts_dict = OrderedDict()
    for j, listHDU in enumerate(partlistHDU):
        part_index = j+1
        for i, HDU in enumerate(listHDU):
            if HDU.data is None: continue
            if 'exclude' in args and args.exclude is not None and HDU.header['OHDU'] in args.exclude: continue
            if 'ohdu' in args and args.ohdu is not None and HDU.header['OHDU'] not in args.ohdu: continue
            if 'plot_part' in args:
                title= 'o{}_p{}'.format(HDU.header['OHDU'], part_index)
                Section(HDU.data).save_pdf( '{0}{1}.pdf'.format(args.name, title) )
            part = Part(HDU)
            params = get_median_params_by_line( part.data.get_half(), part.bias )
            params.update( {'type':'median_by_lines', 'ohdu':'{}p{}'.format( HDU.header['OHDU'], part_index ) } )
            table.append( params )

            part.remove_outliers_from_bias()
            part.correct_lines(np.nanmean)
            try: parts_dict[i].append(part)
            except KeyError: parts_dict[i] = [part]

    first = True
    param_list = OrderedDict()
    for i, parts in parts_dict.items():
        image = Image( parts )
        ohdu = image.header['OHDU']
        args.name += '{}'.format( args.params_mode )
        params = OrderedDict()

        image.remove_outliers_from_vbias()
        params = get_mean_params_by_line( image.data.get_half(), image.bias )
        params.update( {'type': 'mean_by_line', 'ohdu': ohdu } )
        table.append( params )

        image.remove_outliers_from_vbias_second_pass()

        image.correct_cols( np.nanmean )
        #print( 'total charge', np.sum( image.data ), np.nanmean( image.data ) )

        params = get_mean_params_by_line( image.data.get_half().remove_hits( 60, border=3 ), image.bias )
        params.update( {'type': 'mean_by_line_remove60:3', 'ohdu': ohdu } )
        table.append( params )


        if 'remove_hits' in args:
            args.name += '_e{0}b{border}'.format( args.remove_hits[0], border=args.remove_hits[1] )
            image.data = image.data.remove_hits( args.remove_hits[0], border=args.remove_hits[1])

        if 'find_hits' in args:
            bias_std = np.nanstd(image.bias)
            print_var( 'bias_std', locals() )
            outliers2nanPoissonNorm_1d( image.data, sigma = np.nanstd(image.bias), pmin = args.find_hits )

        if 'plot_sections' in args:
            image.save_sections( r'{}_o{}'.format( args.name, image.header['OHDU'] ) )
            image.bias.save_projection( r'{}_o{}_biasProj'.format( args.name, image.header['OHDU'] ), title='bias', axis=1, no_max=True, no_mean = True, no_min = True, do_regress=True )
            image.vbias.save_projection( r'{}_o{}_vbiasProj'.format( args.name, image.header['OHDU'] ), title='vbias', axis=0, no_max=True, no_mean = True, no_min = True, do_regress=True )
            image.data.save_projection( r'{}_o{}_dataProj1'.format( args.name, image.header['OHDU'] ), title='data', axis=1, no_max=True, no_mean = True, no_min = True, do_regress=True )
            image.data.save_projection( r'{}_o{}_dataProj0'.format( args.name, image.header['OHDU'] ), title='data', axis=0, no_max=True, no_mean = True, no_min = True, do_regress=True )

        if 'plot_spectrum' in args:
            for sec in ['dbias', 'vbias', 'bias', 'data']:
                title = r'o{}_{}_spectrum'.format( image.header['OHDU'], sec )
                getattr(image, sec).save_spectrum( '{}_{}'.format(args.name, title ), title=title )

        if 'plot_convolution_spectrum' in args:
            func = np.nansum
            if 'convolution_function' in args:
                func = eval( 'lambda x:' + args.convolution_function )
            for sec in ['dbias', 'vbias', 'bias', 'data']:
                for length in [5,10]:
                    title = r'o{ohdu}_{sec}_convolution{length}'.format( ohdu=image.header['OHDU'], sec=sec, length=length )
                    getattr( image, sec ).convolve( function=func, size=(length,length) )\
                        .save_spectrum( '{}_{}'.format(args.name, title ), title=title )

        if 'plot_block_spectrum' in args:
            func = np.nansum
            if 'block_function' in args:
                func = eval( 'lambda x:' + args.block_function )
            if 'block_function_mean' in args:
                func = np.nanmean
            if 'block_function_std' in args:
                func = lambda x: (lambda y: np.nan if y==0 else y)(np.nanstd(x))

            for sec in ['dbias', 'vbias', 'bias', 'data']:
                for length in [5,10]:
                    title = r'o{ohdu}_{sec}_block{length}'.format( ohdu=image.header['OHDU'], sec=sec, length=length )
                    section = getattr( image, sec )
                    h, w = section.shape
                    section = section[ :-(h%length), :( -(w%length) if (w%length)>0 else None) ]
                    blocks = blockshaped( section, length, length )
                    ret = [ func(block) for block in blocks ]
                    Section(ret).save_spectrum( '{}_{}'.format(args.name, title ), title=title, binsize=.01 )

        if 'plot_dark_current_evolution' in args:
            func = np.nanmedian
            for sec in ['dbias', 'vbias', 'bias', 'data']:
                for length in [5,10]:
                    title = r'o{ohdu}_{sec}_block{length}'.format( ohdu=image.header['OHDU'], sec=sec, length=length )
                    section = getattr( image, sec )
                    medians = func( section, axis = 1 )
                    Section(ret).save_spectrum( '{}_{}'.format(args.name, title ), title=title, binsize=.01 )

        #fmt = lambda x: {int:'{:4}', float:'{: .3e}', Section:'{: .3e}'}[type(x)]
        #sfmt = lambda x: {int:'{:4.4}', float:'{:10.10}', Section:'{:10.10}'}[type(x)]

        #fname = '{}_params.csv'.format( args.name )
        #open( fname, 'w' )
        #if first:
            #columns = zip(*params)[0]
            #max_length = max( map(len, columns) )
            #print( ' '.join( [text('ohdu', mode='B')] + [ text( sfmt(v).format(key), mode='B') for key, v in params ] ) )

            #open( fname, 'a' ).write( '# ' + ', '.join( ['ohdu'] + [ key for key, v in params ] ) + '\n' )
            #first = False
        #column_head = '%4d' % image.header['OHDU']
        #print( ' '.join( [text( column_head, mode='B' )] + [ fmt(v).format(v) for keys, v in params ] ) )
        #open( fname, 'a' ).write( ', '.join( [column_head] + [ fmt(v).format(v) for keys, v in params ] ) + '\n' )
    print(table)
    return table
