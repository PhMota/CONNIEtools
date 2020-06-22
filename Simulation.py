#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
try:
    from numpy import *
except ImportError:
    print('missing module, please run')
    print('module load softwares/python/2.7-gnu-5.3')
    exit(0)
    
import scipy.stats
from numpy.lib import recfunctions as rfn
import json
import os
import astropy.io.fits
from Timer import Timer
from termcolor import colored, cprint
from PrintVar import print_var
import argparse
import sys
import constants

electron_in_eV = 3.745
electron_in_keV = electron_in_eV*1e-3
Cu_energy_eV = 8046
Cu2_energy_eV = 8904
Si_energy_eV = 1740


def simulation_from_file( input_file ):
    args = json.loads( open( input_file ).readline().strip(' #').replace('\'','"') )
    data = genfromtxt( input_file, delimiter= ', ', names = True, skip_header = 1, dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] ).view(Simulation)
    return Simulation( data, args )
    
class Simulation( recarray ):
    
    def __new__(cls, data, args ):

        obj = data.view(Simulation)
        
        for key, arg in vars(args).items():
            setattr( obj, key, arg )
        if 'image_mode' in args:
            if args.image_mode == 1:
                obj.ccd_shape = constants.ccd_shape
                obj.horizontal_overscan = 150
                obj.vertical_overscan = 90
                obj.rebin = [1,1]
                obj.expose_hours = 3
            elif args.image_mode == 5:
                obj.ccd_shape = constants.ccd_shape
                obj.horizontal_overscan = 450
                obj.vertical_overscan = 74
                obj.rebin = [5,1]
                obj.expose_hours = 1
            else:
                print( 'image_mode {} not recognized. Ignoring.'.format(obj.image_mode) )
        
        obj.ccd_shape = array(obj.ccd_shape)
        obj.rebinned_ccd_shape = obj.ccd_shape/obj.rebin
        obj.image_shape = obj.rebinned_ccd_shape + [obj.vertical_overscan, obj.horizontal_overscan]
        obj.diffusion_function = eval( 'vectorize(lambda z: {})'.format( obj.diffusion_function ) )
        obj.charge_efficiency_function = eval( 'vectorize(lambda z: {})'.format( obj.charge_efficiency_function ) )
        obj.vertical_modulation_function = eval( 'vectorize(lambda y: {})'.format( obj.vertical_modulation_function ) )
        obj.horizontal_modulation_function = eval( 'vectorize(lambda x: {})'.format( obj.horizontal_modulation_function ) )
        
        try:
            obj.count = len( obj.q )
        except TypeError:
            obj.count = 1
        return obj
    
    def get_charge(self):
        return ( self.charge_efficiency_function(self.z) * self.q ).astype(int)
    
    def get_xy( self ):
        if self.count == 1:
            return array( ([self.x], [self.y]) ).T
        return array( (self.x, self.y) ).T

    #def get_xy_rebin( self ):
        #if self.count == 1:
            #return (array( ([self.x], [self.y]) )/self.rebin[:,None]).T
        #return (array( (self.x, self.y) )/self.rebin[:,None]).T

    def get_E( self ):
        return self.q*self.charge_gain
    
    def get_id_code( self, entry = None ):
        map_id = {'random': 1, 'Cu': 11, 'Cu2': 12, 'Si':10}
        if entry is None:
            return map( lambda i: map_id[i], self.id )
        return map_id[entry.id]

    def get_sigma( self, entry = None ):
        if entry is None:
            return self.diffusion_function( self.z )
        return self.diffusion_function( entry.z )
        
    def generate_image( self, args ):
        if self.count > 1:
            sigma_per_event = self.get_sigma()
            charge_per_event = self.get_charge()
            Q = sum(charge_per_event)
            xy_norm = scipy.stats.norm.rvs( size = 2*Q ).reshape( -1, 2 )
            sigma_per_charge = repeat( sigma_per_event, charge_per_event )
            xy_per_charge = repeat( self.get_xy(), charge_per_event, axis=0 )
            xy = xy_per_charge + sigma_per_charge[:,None] * xy_norm
            bins = [
                arange(self.ccd_shape[0]+1),
                arange(self.ccd_shape[1]+1)
                ]
            image = histogramdd( xy, bins = bins )[0] * self.charge_gain
        elif self.count == 1:
            sigma_per_event = self.get_sigma()
            charge_per_event = self.get_charge()
            Q = sum(charge_per_event)
            xy_norm = scipy.stats.norm.rvs( size = 2*Q ).reshape( -1, 2 )
            sigma_per_charge = repeat( sigma_per_event, charge_per_event )
            xy_per_charge = repeat( self.get_xy(), charge_per_event, axis=0 )
            xy = xy_per_charge + sigma_per_charge[:,None] * xy_norm
            bins = [
                arange(self.ccd_shape[0]+1),
                arange(self.ccd_shape[1]+1)
                ]
            image = histogramdd( xy, bins = bins )[0] * self.charge_gain
        else:
            image = zeros( self.ccd_shape )
        #print('total hit charge', sum(image) )
        #total_hit_pixels = sum(image>0)
        #print('total hit pixels', total_hit_pixels, float(total_hit_pixels)/image.size )

        if 'verbose' in args:
            max_ADU = max( image.flatten() )
            max_charge = max_ADU/self.charge_gain
            print_var( ['max_charge', 'max_ADU'], locals() )
        
        if self.dark_current > 0:
            dark_current_image = scipy.stats.poisson.rvs( self.dark_current, size = self.ccd_shape[0]*self.ccd_shape[1] ).reshape(self.ccd_shape)
            #max_charge = max( dark_current_image.flatten() )
            #total_dark_charge = sum( dark_current_image.flatten() )
            #print( 'total dark charges', total_dark_charge, max_charge )
            #total_pixels = sum(dark_current_image>0)
            #print( 'total dark pixels', total_pixels, float(total_pixels)/dark_current_image.size )
            #print( 'total dark pixels', sum( dark_current_image[:100,:100] ), float(sum( dark_current_image[:100,:100]>0 ))/(100*100) )
            image += dark_current_image * self.charge_gain
            if 'verbose' in args:
                max_ADU = max_charge * self.charge_gain
                print_var( ['max_charge', 'max_ADU', 'total_dark_charge'], locals() )
        
        image = image.reshape( (self.rebinned_ccd_shape[0], -1, self.rebinned_ccd_shape[1]) ).sum(axis = 1)

        new_image = zeros( self.image_shape )
        new_image[self.vertical_overscan:, :-self.horizontal_overscan] = image
        image = new_image

        if 'verbose' in args:
            print( 'rebin' )
            print_var( ['ccd_shape', 'rebinned_ccd_shape', 'image_shape'], vars(self) )
            max_ADU = max( image.flatten() )
            max_charge = max_ADU/self.charge_gain
            print_var( ['max_charge', 'max_ADU'], locals() )
        
        if self.readout_noise > 0:
            readout_noise_image = scipy.stats.norm.rvs( size=self.image_shape[0]*self.image_shape[1] ).reshape(self.image_shape)*self.readout_noise
            image += readout_noise_image

        if 'verbose' in args:
            max_ADU = max( image.flatten() )
            max_charge = max_ADU/self.charge_gain
            print( 'noise' )
            print_var( ['max_charge', 'max_ADU'], locals() )

        vertical_modulation = array([ self.vertical_modulation_function(y) for y in range(self.image_shape[0]) ])[:,None]
        image += vertical_modulation

        horizontal_modulation = array([ self.horizontal_modulation_function(x) for x in range(self.image_shape[1]) ])[None,:]
        image += horizontal_modulation

        if 'verbose' in args:
            print( 'modulations' )
            max_ADU = max( image.flatten() )
            max_charge = max_ADU/self.charge_gain
            print_var( ['max_charge', 'max_ADU'], locals() )

        if 'pdf' in args:
            self.save_pdf( image, args )
        
        hdu_list = self.make_fits( image, args )
        if not 'no_fits' in args: self.save_fits( hdu_list, args )
        return hdu_list

    def make_fits( self, image, args ):
        with Timer('generated fits'):
            header = astropy.io.fits.Header()
            header['biasW'] = self.horizontal_overscan
            header['biasH'] = self.vertical_overscan
            header['rebin'] = max(self.rebin)
            header['gain'] = self.charge_gain
            header['noise'] = self.readout_noise
            header['dc'] = self.dark_current
            header['OHDU'] = -1

            fits_image = astropy.io.fits.ImageHDU( image.astype(args.image_type), header=header )
            header['OHDU'] = 0
            primary = astropy.io.fits.PrimaryHDU( header=header )
            hdu_list = astropy.io.fits.HDUList([primary, fits_image])
            #return fits_image
            return hdu_list
    
    def save_fits( self, hdu_list, args ):
        if 'compress' in args:
            fname = args.name + '.fits.' + args.compress
        else:
            fname = args.name + '.fits'
        with Timer('saved ' + fname ):
            hdu_list.writeto( fname, overwrite=True )
            return
    
    def save_pdf(self, output,  ):
        fname = args.output + '.pdf'
        with Timer( 'saved ' + fname ):
            from matplotlib import pylab as plt
            from matplotlib.patches import Ellipse
            from matplotlib.colors import LogNorm

            fig = plt.figure()
            ax = fig.add_subplot(111)

            #x, y = self.x-.5, self.y-.5
            #ax.scatter( x, y, c='r', s=1 )
            #if self.count == 1:
                #sigma = sqrt( self.get_sigma() )*5
                #ax.add_artist( Ellipse((self.y-.5, self.x-.5), sigma, sigma, fill=False ) )
            #else:
                #for event in self:
                    #sigma = sqrt( self.get_sigma(event) )*5
                    #ax.add_artist( Ellipse((event.y-.5, event.x-.5), sigma, sigma, fill=False ) )
            im = log(image - nanmin(image) + 1)
            #im = image
            ax.imshow( im, cmap='Blues', origin='lower', vmin=nanmin(im), vmax=nanmax(im) )
            fig.savefig( fname )
        return

def add_params_options(p):
    g = p.add_argument_group('parameter options')
    g.add_argument('-g', '--charge-gain', help = 'factor to convert charges into ADU',
                   type=float, default = 7.25 )    
    g.add_argument('-rn', '--readout-noise', help = 'sigma of the normal noise distribution in ADU', 
                   type=float, default = 0 )
    g.add_argument('-dc', '--dark-current', help = 'lambda of Poisson distribution in 1/(e-·h)',
                   type=float, default = 0 )
    g.add_argument('-exp', '--expose-hours', help = 'exposed hours',
                   type=float, default = 1 )
    return

def add_geometry_options(p):
    g = p.add_argument_group('geometry options')
    g.add_argument('-os', '--horizontal-overscan', type=int, default = 150, help = 'size of the horizontal overscan in pixels' )
    g.add_argument('-vos', '--vertical-overscan', type=int, default = 90, help = 'size of the vertical overscan in pixels' )
    #p.add_argument('-hpt', '--horizontal-pretrim', type=int, default = 8, help = 'size of the vertical overscan in pixels' )
    #p.add_argument('-vpt', '--vertical-pretrim', type=int, default = 1, help = 'size of the vertical overscan in pixels' )
    g.add_argument('--ccd-shape', nargs=2, type=int, default = constants.ccd_shape, help = 'shape of the image as 2d pixels' )
    g.add_argument('--rebin', nargs=2, type=int, default = [1,1], help = '2d rebinning strides' )
    g.add_argument('--image-type', type=eval, default='int', help = 'image type' )
    g.add_argument('--image-mode', help = 'set to "1" to use official 1x1 image geomtry or "5" to 1x5', 
                    type=int, default = argparse.SUPPRESS )
    return

def add_depth_options(p):
    g = p.add_argument_group('depth options')
    g.add_argument('--depth-range', nargs=2, type=float, default = [0,670], help = 'range into which to randomly generate depths' )
    g.add_argument('--diffusion-function',
                        type=str, 
                        default = '{}'.format(constants.diffusion_function),
                        help = 'function to map z-depth into xy-sigma' 
                        )
    g.add_argument('--charge-efficiency-function',
                        type=str,
                        default = constants.charge_efficiency_function,
                        help = 'function for charge efficiency dependent of z-depth' 
                        )
    return

def add_charges_options(p):
    g = p.add_argument_group('charge options')
    g.add_argument('-N', '--number-of-charges', type=int, default = 0, help = 'number of charges to be randomly generated' )
    g.add_argument('--charge-range', nargs=2, type=int, default = [5, 200], help = 'range into which to randomly generate charges' )
    g.add_argument('--number-of-Cu-charges',
                        type=int,
                        default = 0,
                        help = 'number of charges to be randomly generated at the Copper fluorescence energy 8.046keV' 
                        )
    g.add_argument('--number-of-Cu2-charges', type=int, default = 0, 
                        help = 'number of charges to be randomly generated at the secundary Copper fluorescence energy 8.904keV' )
    g.add_argument('--number-of-Si-charges', type=int, default = 0, help = 'number of charges to be randomly generated at the Silicon fluorescence energy 1.740keV' )
    return

def add_modulation_options(p):
    g = p.add_argument_group('modulation options')
    g.add_argument('--vertical-modulation-function', help = 'function to modulate the vertical axis', 
                    type=str, default = "0" )
    g.add_argument('--horizontal-modulation-function', help = 'function to modulate the horizontal axis',
                        type=str, default = "0" )
    g.add_argument('--default-vertical-modulation', help = 'set vertical modulation to "{}"'.format(constants.vertical_modulation_function),
                    type=str, default=argparse.SUPPRESS )
    g.add_argument('--default-horizontal-modulation', help = 'set horizontal modulation to "{}"'.format(constants.horizontal_modulation_function), 
                    type=str, default=argparse.SUPPRESS )
    g.add_argument('--default-modulation', help = 'set modulations to "{}" and "{}"'.format(constants.horizontal_modulation_function, constants.vertical_modulation_function), 
                    type=str, default=argparse.SUPPRESS )
    return

def add_output_options(p):
    g = p.add_argument_group('output options')
    g.add_argument('--no-fits', help = 'suppress fits output',
                   action='store_true', default=argparse.SUPPRESS )
    g.add_argument('--pdf', help = 'generate pdf output',
                   action='store_true', default=argparse.SUPPRESS )
    g.add_argument('--spectrum', help = 'generate energy spectrum',
                   action='store_true', default=argparse.SUPPRESS )
    g.add_argument('--verbose', help = 'verbose level', 
                   type=int, default=argparse.SUPPRESS )
    g.add_argument('--csv', help = 'generate csv output',
                   action='store_true', default=argparse.SUPPRESS )
    return

def add_image_options( p ):
    p.add_argument('name', help = 'factor to convert charges into ADU',
                   type=str, default = 'simulated_image' )
    add_params_options(p)
    add_geometry_options(p)
    add_depth_options(p)
    add_charges_options(p)
    add_modulation_options(p)
    add_output_options(p)
    p.set_defaults( func=image )

def image( args ):
    sim = simulate_events( args )
    return sim.generate_image( args )

def simulate_events( args ):
    ccd_shape = args.ccd_shape
    depth_range = args.depth_range
    charge_range = args.charge_range
    number_of_charges = args.number_of_charges
    array_of_positions = random.random( number_of_charges*2 ).reshape(-1,2)*(array(ccd_shape)-1)
    array_of_depths = random.random( number_of_charges )*(depth_range[1] - depth_range[0]) + depth_range[0]
    array_of_charges = random.random( number_of_charges )*(charge_range[1] - charge_range[0]) + charge_range[0]
    array_of_identities = ['random']*number_of_charges
    
    list_of_X_charges = [['number_of_Cu_charges', Cu_energy_eV, 'Cu'], 
                         ['number_of_Cu2_charges', Cu2_energy_eV, 'Cu2'], 
                         ['number_of_Si_charges', Si_energy_eV, 'Si'] ]
    for key, energy_eV, charge_type in list_of_X_charges:
        if key in args:
            number_of_X_charges = getattr(args, key)
            array_of_positions = concatenate( 
                (array_of_positions, 
                 random.random( number_of_X_charges*2 ).reshape(-1,2)*(array(ccd_shape)-1) ), 
                axis=0 )
            array_of_depths = append( 
                array_of_depths, 
                random.random( number_of_X_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] 
                )
            array_of_charges = append( array_of_charges, [energy_eV/electron_in_eV]*number_of_X_charges )
            array_of_identities.extend( [charge_type]*number_of_X_charges )
    
    if 'verbose' in args and args.verbose == 0: 
        print( 'total charges created', len( array_of_identities ) )
    
    cols = ['x', 'y', 'z', 'q', 'id']
    types = [float, float, float, float, 'S16']
    fmt = ['%.4f', '%.4f', '%.4f', '%.4f', '%s']
    dtype = zip( cols, types )
    data = empty( (0,len(dtype)), dtype=dtype )
    if len(array_of_depths) > 0:
        data = array( zip(array_of_positions[:,0], array_of_positions[:,1], array_of_depths, array_of_charges, array_of_identities), dtype=dtype )
        
    if 'csv' in args:
        output = args.name + '.csv'
        header = str(vars(args))
        #print( 'header', header )
        header += '\n' + ', '.join(cols)
        if len(array_of_depths) > 0:
            savetxt( output, data, delimiter=', ', header=header, fmt=fmt )
        else:
            open( output, 'w' ).writelines( header )
    return Simulation( data, args )

def add_folder_options( p ):
    p.add_argument('folder_name', help = 'factor to convert charges into ADU',
                   type=str, default = 'simulated_image' )
    p.add_argument('number_of_images', type=int, default = None, help = 'number of images to be generated' )
    g = p.add_argument_group('parameter range options')
    g.add_argument('--readout-noise-range', nargs=2, type=float, default = [0, 0], help = 'readout_noise range' )
    g.add_argument('--dark-current-range', nargs=2, type=float, default = [0, 0], help = 'dark current range' )
    g.add_argument('--charge-gain-range', nargs=2, type=float, default = [7.25, 7.25], help = 'charge gain range' )
    add_geometry_options(p)
    add_depth_options(p)
    add_charges_options(p)
    add_modulation_options(p)
    add_output_options(p)
    p.set_defaults( func=generate_folder )

def generate_folder( args ):
    if os.path.exists( args.folder_name ):
        print( 'output already exists. Exiting' )
        exit(0)
    os.mkdir( args.folder_name )
    if not 'number_of_images' in args:
        print( 'missing number of images' )
        exit(0)
    for count in range( args.number_of_images ):
        print( 'generating image', count+1, 'out of', args.number_of_images )
        args.readout_noise = random.uniform(*args.readout_noise_range)
        args.dark_current = random.uniform(*args.dark_current_range)
        args.charge_gain = random.unfiform(*args.charge_gain_range)
        args.count = count
        args.name = '{folder_name}/{folder_name}{count:04d}RN{readout_noise:.2f}DC{dark_current:.4f}CG{charge_gain:.1f}'.format(**vars(args)) 
        sim = simulate_events( args )
        sim.generate_image( args )
    return
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'simulation tools', formatter_class = argparse.ArgumentDefaultsHelpFormatter 
                                     )
    subparsers = parser.add_subparsers( help = 'a' )

    add_image_options( subparsers.add_parser('image', help='generate simulated image', formatter_class=argparse.ArgumentDefaultsHelpFormatter) )
    add_folder_options( subparsers.add_parser('folder', help='generate folder with simulated images', formatter_class=argparse.ArgumentDefaultsHelpFormatter) )

    args = parser.parse_args()
    
    with Timer('finished'):
        func = args.func
        del args.func
        print( colored('using parameters:', 'green', attrs=['bold'] ) )
        print_var( vars(args).keys(), vars(args), line_char='\t' )
        func(args)
