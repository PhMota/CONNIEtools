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
import time
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

class to_object:
    def __init__(self, d):
        self.__dict__ = d
    def __contains__(self, a):
        return a in self.__dict__.keys()

def simulation_from_file( basename ):
    fname = basename+'.csv'
    json_str = open( fname ).readline().strip(' #').replace('\'','"')
    args_dict = json.loads( json_str )
    args = to_object(args_dict)
    if 'verbose' in args:
        del args.verbose
    data = genfromtxt( fname, delimiter= ', ', names = True, skip_header = 1, dtype = [('n', int), ('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] ).view(recarray)
    print( 'read events', data.size, data.n )
    if data.size == 1:
        data = array( [(data.n, data.x, data.y, data.z, data.q, data.id)], dtype = [('n', int), ('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] ).view(recarray)
        print( 'read events', data.size, data.n )
        
    return Simulation( data, args )
    
class Simulation:
    def __updateattrs__(self, data, verbose=0 ):
        self.__recarray__ = data
        for name in self.__recarray__.dtype.names:
            if verbose:
                print( name )
            setattr(self, name, getattr(self.__recarray__, name) )
    
    def __init__(self, data, args ):
        if not isinstance(data, recarray):
            print( 'type(data)', type(data) )
        verbose = 0
        if 'verbose' in args:
            verbose = 1
        self.__updateattrs__( data, verbose )
        
        for key, arg in vars(args).items():
            if 'verbose' in args:
                print( key, arg )
            setattr( self, key, arg )
        if 'image_mode' in args:
            if args.image_mode == 1:
                self.ccd_shape = constants.ccd_shape
                self.horizontal_overscan = 150
                self.vertical_overscan = 90
                self.rebin = [1,1]
                self.expose_hours = 3
            elif args.image_mode == 5:
                self.ccd_shape = constants.ccd_shape
                self.horizontal_overscan = 450
                self.vertical_overscan = 74
                self.rebin = [5,1]
                self.expose_hours = 1
            else:
                print( 'image_mode {} not recognized. Ignoring.'.format(self.image_mode) )
        
        self.rebin = array(self.rebin)
        self.ccd_shape = array(self.ccd_shape)
        self.rebinned_ccd_shape = self.ccd_shape/self.rebin
        self.image_shape = self.rebinned_ccd_shape + [self.vertical_overscan, self.horizontal_overscan]
        self.diffusion_function = eval( 'vectorize(lambda z: {})'.format( self.diffusion_function ) )
        self.charge_efficiency_function = eval( 'vectorize(lambda z: {})'.format( self.charge_efficiency_function ) )
        self.vertical_modulation_function = eval( 'vectorize(lambda y: {})'.format( self.vertical_modulation_function ) )
        self.horizontal_modulation_function = eval( 'vectorize(lambda x: {})'.format( self.horizontal_modulation_function ) )
        
        if len( self.q ) > 1:
            q_eff = ( self.charge_efficiency_function(self.z) * self.q ).astype(int)
            E = self.q * self.charge_gain
            E_eff = q_eff * self.charge_gain
            sigma = self.diffusion_function( self.z )
            map_id = {'random': 1, 'Si':10, 'Cu': 11, 'Cu2': 12}
            id_code = map( lambda i: map_id[i], self.__recarray__.id )
        else:
            q_eff = []
            E = []
            E_eff = []
            sigma = []
            id_code = []
        
        try:
            self.count = len( self.q )
        except TypeError:
            print( 'charges', self.q )
            self.count = 1

        self.__updateattrs__( rfn.append_fields( self.__recarray__, ['q_eff', 'E', 'E_eff', 'sigma', 'id_code'], [q_eff, E, E_eff, sigma, id_code], dtypes=(int, float, float, float, int), asrecarray=True ) )
        
        return
    
    def __len__(self):
        return len(self.__recarray__)
    
    def __getitem__(self, i):
        return self.__recarray__[i]
    
    @property
    def xy( self ):
        if not '__xy' in self.__dict__:
            self.__xy = array( (self.x, self.y) ).T
        return self.__xy
        
    def generate_image( self, args ):
        if self.count > 0:
            Q = sum( self.q_eff )
            xy_norm = scipy.stats.norm.rvs( size = 2*Q ).reshape( -1, 2 )
            sigma_per_charge = repeat( self.sigma, self.q_eff )
            xy_per_charge = repeat( self.xy, self.q_eff, axis=0 )
            xy = xy_per_charge + sigma_per_charge[:,None] * xy_norm
            bins = [
                arange(self.ccd_shape[0]+1),
                arange(self.ccd_shape[1]+1)
                ]
            image = histogramdd( xy, bins = bins )[0] * self.charge_gain
            if 'verbose' in args:
                print( 'bin len', len( bins[0]), min(bins[0]), max(bins[0]), bins[0][1] - bins[0][0] )
                print( 'image.shape', image.shape )
                print( 'xy', self.xy )
                print( argwhere( image > 60 ) )
        else:
            image = zeros( self.ccd_shape )

        if 'verbose' in args:
            max_ADU = max( image.flatten() )
            max_charge = max_ADU/self.charge_gain
            print_var( ['max_charge', 'max_ADU'], locals() )
        
        if self.dark_current > 0:
            dark_current_image = scipy.stats.poisson.rvs( self.dark_current, size = self.ccd_shape[0]*self.ccd_shape[1] ).reshape(self.ccd_shape)
            image += dark_current_image * self.charge_gain
            if 'verbose' in args:
                max_ADU = max_charge * self.charge_gain
                #print_var( ['max_charge', 'max_ADU', 'total_dark_charge'], locals() )
        
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
            self.save_image( image, args.basename+'.pdf' )
        if 'png' in args:
            self.save_image( image, args.basename+'.png' )
        
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
            header['RUNID'] = -1
            now = time.time()
            print( 'now', now )
            header['expstart'] = now
            header['expstop'] = now + args.expose_hours*60*60

            fits_image = astropy.io.fits.ImageHDU( image.astype(eval(args.image_type)), header=header )
            header['OHDU'] = 0
            primary = astropy.io.fits.PrimaryHDU( header=header )
            hdu_list = astropy.io.fits.HDUList([primary, fits_image])
            #return fits_image
            return hdu_list
    
    def save_fits( self, hdu_list, args ):
        if 'compress' in args:
            fname = args.basename + '.fits.' + args.compress
        else:
            fname = args.basename + '.fits'
        with Timer('saved ' + fname ):
            hdu_list.writeto( fname, overwrite=True )
            return
    
    def save_image(self, image, fname ):
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
    g.add_argument('--ccd-shape', nargs=2, type=int, default = constants.ccd_shape.tolist(), help = 'shape of the image as 2d pixels' )
    g.add_argument('--rebin', nargs=2, type=int, default = [1,1], help = '2d rebinning strides' )
    g.add_argument('--image-type', type=str, default='int', help = 'image type' )
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
    g.add_argument('--png', help = 'generate png output',
                   action='store_true', default=argparse.SUPPRESS )
    g.add_argument('--spectrum', help = 'generate energy spectrum',
                   action='store_true', default=argparse.SUPPRESS )
    g.add_argument('--verbose', help = 'verbose level', 
                   action='store_true', default=argparse.SUPPRESS )
    g.add_argument('--csv', help = 'generate csv output',
                   action='store_true', default=argparse.SUPPRESS )
    return

def add_image_options( p ):
    p.add_argument('basename', help = 'basename for simulation output',
                   type=str, default = 'simulated_image' )
    add_params_options(p)
    add_geometry_options(p)
    add_depth_options(p)
    add_charges_options(p)
    add_modulation_options(p)
    add_output_options(p)
    p.set_defaults( func=image )

def image( args ):
    print_options(args)
    sim = simulate_events( args )
    return sim.generate_image( args )

def simulate_events( args ):
    ccd_shape = args.ccd_shape
    depth_range = args.depth_range
    charge_range = args.charge_range
    number_of_charges = args.number_of_charges
    array_of_positions = random.random( number_of_charges*2 ).reshape(-1,2)*array(ccd_shape)
    array_of_depths = random.uniform( *(depth_range+[number_of_charges]) )
    array_of_charges = random.uniform( *(charge_range+[number_of_charges]) )
    array_of_identities = ['random']*number_of_charges
    
    list_of_X_charges = [['number_of_Cu_charges', Cu_energy_eV, 'Cu'], 
                         ['number_of_Cu2_charges', Cu2_energy_eV, 'Cu2'], 
                         ['number_of_Si_charges', Si_energy_eV, 'Si'] ]
    for key, energy_eV, charge_type in list_of_X_charges:
        if key in args:
            number_of_X_charges = getattr(args, key)
            array_of_positions = concatenate( 
                ( array_of_positions, random.random( number_of_X_charges*2 ).reshape(-1,2)*array(ccd_shape) ), 
                axis=0 )
            array_of_depths = append( 
                array_of_depths, random.uniform( *(depth_range+[number_of_X_charges]) )
                )
            array_of_charges = append( array_of_charges, [energy_eV/electron_in_eV]*number_of_X_charges )
            array_of_identities.extend( [charge_type]*number_of_X_charges )
    
    if 'verbose' in args and args.verbose == 0: 
        print( 'total charges created', len( array_of_identities ) )
    
    cols = ['n', 'x', 'y', 'z', 'q', 'id']
    types = [int, float, float, float, float, 'S16']
    fmt = ['%d', '%.4f', '%.4f', '%.4f', '%.4f', '%s']
    dtype = zip( cols, types )
    data = empty( (0,len(dtype)), dtype=dtype ).view(recarray)
    if len(array_of_depths) > 0:
        data = array( zip(arange(len(array_of_depths)), array_of_positions[:,0], array_of_positions[:,1], array_of_depths, array_of_charges, array_of_identities), dtype=dtype ).view(recarray)
        
    if 'csv' in args:
        output = args.basename + '.csv'
        header = json.dumps(vars(args))
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
        args.basename = '{folder_name}/{folder_name}{count:04d}RN{readout_noise:.2f}DC{dark_current:.4f}CG{charge_gain:.1f}'.format(**vars(args)) 
        sim = simulate_events( args )
        sim.generate_image( args )
    return
    
def print_options( args ):
    print( colored('using parameters:', 'green', attrs=['bold'] ) )
    print_var( vars(args).keys(), vars(args), line_char='\t' )
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'simulation tools', formatter_class = argparse.ArgumentDefaultsHelpFormatter 
                                     )
    subparsers = parser.add_subparsers( help = 'a' )

    add_image_options( subparsers.add_parser('image', help='generate simulated image', formatter_class=argparse.ArgumentDefaultsHelpFormatter) )
    add_folder_options( subparsers.add_parser('folder', help='generate folder with simulated images', formatter_class=argparse.ArgumentDefaultsHelpFormatter) )

    args = parser.parse_args()
    func = args.func
    del args.func
    with Timer('finished'):
        func(args)
