# coding: utf-8
from __future__ import print_function
from numpy import *
import scipy.stats
from numpy.lib import recfunctions as rfn
import json
import os
import astropy.io.fits
from Timer import Timer
from TerminalColor import text

electron_in_eV = 3.745
electron_in_keV = electron_in_eV*1e-3
Cu_energy_eV = 8046
Cu2_energy_eV = 8904
Si_energy_eV = 1740

def generate_folder( number_of_images, readout_noise_range, dark_current_range, output_pattern, **kwargs ):
    count = 1
    for count in range( number_of_images ):
        print( 'generating image', count, 'out of', number_of_images )
        kwargs['readout_noise'] = (readout_noise_range[1] - readout_noise_range[0])*random.random() + readout_noise_range[0]
        kwargs['dark_current'] = (dark_current_range[1] - dark_current_range[0])*random.random() + dark_current_range[0]
        kwargs['output_file'] = output_pattern.replace('*', 'RN{readout_noise:.3}DC{dark_current:.3}'.format(**kwargs) )
        print( kwargs['image_fits_output'] )
        sim = simulate_events( **kwargs )
        sim.generate_image( output = kwargs['output_file'] )
        count += 1
    return

def simulate_events( args ):
    xyshape = args.xyshape
    depth_range = args.depth_range
    charge_range = args.charge_range
    number_of_charges = args.number_of_charges
    array_of_positions = random.random( number_of_charges*2 ).reshape(-1,2)*(array(xyshape)-1)
    array_of_depths = random.random( number_of_charges )*(depth_range[1] - depth_range[0]) + depth_range[0]
    array_of_charges = random.random( number_of_charges )*(charge_range[1] - charge_range[0]) + charge_range[0]
    array_of_identities = ['random']*number_of_charges
    
    if 'number_of_Cu_charges' in args.__dict__:
        number_of_Cu_charges = args.number_of_Cu_charges
        array_of_positions = concatenate( (array_of_positions, random.random( number_of_Cu_charges*2 ).reshape(-1,2)*(array(xyshape)-1) ), axis=0 )
        array_of_depths = append( array_of_depths, random.random( number_of_Cu_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] )
        array_of_charges = append( array_of_charges, [Cu_energy_eV/electron_in_eV]*number_of_Cu_charges )
        array_of_identities.extend( ['Cu']*number_of_Cu_charges )
        
    if 'number_of_Cu2_charges' in args.__dict__:
        number_of_Cu2_charges = args.number_of_Cu2_charges
        array_of_positions = concatenate( (array_of_positions, random.random( number_of_Cu2_charges*2 ).reshape(-1,2)*(array(xyshape)-1) ), axis=0 )
        array_of_depths = append( array_of_depths, random.random( number_of_Cu2_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] )
        array_of_charges = append( array_of_charges, [Cu2_energy_eV/electron_in_eV]*number_of_Cu2_charges )
        array_of_identities.extend( ['Cu2']*number_of_Cu2_charges )
        
    if 'number_of_Si_charges' in args.__dict__:
        number_of_Si_charges = args.number_of_Si_charges
        array_of_positions = concatenate( (array_of_positions, random.random( number_of_Si_charges*2 ).reshape(-1,2)*(array(xyshape)-1) ), axis=0 )
        array_of_depths = append( array_of_depths, random.random( number_of_Si_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] )
        array_of_charges = append( array_of_charges, [Si_energy_eV/electron_in_eV]*number_of_Si_charges )
        array_of_identities.extend( ['Si']*number_of_Si_charges )

    print( 'total charges created', len( array_of_identities ) )
    
    data = empty( (0,5), dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] )
    if len(array_of_depths) > 0:
        data = array( zip(array_of_positions[:,0], array_of_positions[:,1], array_of_depths, array_of_charges, array_of_identities), dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] )
        
    if args.sim:
        output = args.output_fits + '.csv'
        header = str(vars(args))
        print( 'header', header )
        header += '\nx, y, z, q, id'
        if len(array_of_depths) > 0:
            savetxt( output, data, delimiter = ', ', header = header, fmt = ['%s']*4 + ['%s'] )
        else:
            open( output, 'w' ).writelines( header )
    return Simulation( data, args )

def simulation_from_file( input_file ):
    args = json.loads( open( input_file ).readline().strip(' #').replace('\'','"') )
    data = genfromtxt( input_file, delimiter= ', ', names = True, skip_header = 1, dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] ).view(Simulation)
    return Simulation( data, args )
    
class Simulation( recarray ):
    
    def __new__(cls, data, args ):

        obj = data.view(Simulation)
        
        for key, arg in vars(args).items():
            setattr( obj, key, arg )
            
        obj.xyrebinshape = array(obj.xyshape)/obj.rebin
        obj.fullshape = ( obj.xyrebinshape[0] + obj.vertical_overscan, obj.xyrebinshape[1] + obj.horizontal_overscan )
        obj.diffusion_function = eval( 'vectorize(lambda z: %s)' % obj.diffusion_function )
        obj.charge_efficiency_function = eval( 'vectorize(lambda z: %s)' % obj.charge_efficiency_function )
        obj.vertical_modulation_function = eval( 'vectorize(lambda y: %s)' % obj.vertical_modulation_function )
        obj.horizontal_modulation_function = eval( 'vectorize(lambda x: %s)' % obj.horizontal_modulation_function )
        
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
        
    def generate_image( self, output = None, pdf_output = None ):
        if self.count > 1:
            sigma_per_event = self.get_sigma()
            charge_per_event = self.get_charge()
            Q = sum(charge_per_event)
            xy_norm = scipy.stats.norm.rvs( size = 2*Q ).reshape( -1, 2 )
            sigma_per_charge = repeat( sigma_per_event, charge_per_event )
            xy_per_charge = repeat( self.get_xy(), charge_per_event, axis=0 )
            xy = xy_per_charge + sigma_per_charge[:,None] * xy_norm
            bins = [
                arange(self.xyshape[0]+1),
                arange(self.xyshape[1]+1)
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
                arange(self.xyshape[0]+1),
                arange(self.xyshape[1]+1)
                ]
            image = histogramdd( xy, bins = bins )[0] * self.charge_gain
        else:
            image = zeros( self.xyshape )
            
        if self.dark_current > 0:
            dark_current_image = scipy.stats.poisson.rvs( self.dark_current, size = self.xyshape[0]*self.xyshape[1] ).reshape(self.xyshape)
            image += dark_current_image * self.charge_gain
        
        image = image.reshape( (self.xyrebinshape[0], self.xyrebinshape[1], -1) ).sum(axis = -1)
        
        new_image = zeros( self.fullshape )
        new_image[self.fullshape[0]-self.xyrebinshape[0]:, :self.xyrebinshape[1]] = image
        image = new_image
        
        if self.readout_noise > 0:
            readout_noise_image = scipy.stats.norm.rvs( size=self.fullshape[0]*self.fullshape[1] ).reshape(self.fullshape)*self.readout_noise
            image += readout_noise_image
        
        print( 'image shapes', self.xyshape, self.xyrebinshape, self.fullshape )
        vertical_modulation = array([ self.vertical_modulation_function(y) for y in range(self.fullshape[0]) ])[:,None]
        image += vertical_modulation

        horizontal_modulation = array([ self.horizontal_modulation_function(x) for x in range(self.fullshape[1]) ])[None,:]
        image += horizontal_modulation

        if pdf_output:
            self.save_pdf( output, image )
        
        fits = self.make_fits( image )
        if output: self.save_fits( fits, output )
        return fits

    def make_fits( self, image ):
        with Timer('generated fits'):
            header = astropy.io.fits.Header()
            header['OHDU'] = 0
            header['biasW'] = self.horizontal_overscan
            header['biasH'] = self.vertical_overscan
            header['rebin'] = max(self.rebin)
            header['gain'] = self.charge_gain
            header['noise'] = self.readout_noise
            header['dc'] = self.dark_current
            header['OHDU'] = -1

            fits_image = astropy.io.fits.ImageHDU( image.astype(int), header=header )
            return fits_image
    
    def save_fits( self, fits_image, output ):
        with Timer('saved ' + output ):
            if os.path.exists( output ): os.remove(output)
            header = fits_image.header
            header['OHDU'] = 0
            primary = astropy.io.fits.PrimaryHDU( header=header )
            hdu_list = astropy.io.fits.HDUList([primary, fits_image])
            hdu_list.writeto( output, overwrite=True )
            return
    
    def save_pdf(self, output, image ):
        output += '.pdf'
        with Timer( 'saved ' + output ):
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
            im = log(image - image.min() + 1)
            #im = image
            ax.imshow( im, cmap='Blues', origin='lower', vmin=im.min(), vmax=im.max() )
            fig.savefig( output )
        return

def image( args ):
    del args.func
    for option, value in vars(args).items():
        print( '\t%s: %s' % ( text(option, mode='B'), value) )

    sim = simulate_events( args )
    return sim.generate_image( output=args.output_fits, pdf_output=args.pdf )
    

default_diffusion_function = u'sqrt(-923.666*log1p(-0.000441113*z))/15 if z < 670 else 0'
default_charge_efficiency_function = u'1. if z < 670 else .9'
default_vertical_modulation_function = u'50*cos(y/2/pi/20)'
default_horizontal_modulation_function = u'-1e3*(x/1000. - 1)**2 if x < 1000. else 0'

def tuple_of( type_ ):
    return lambda x: map( type_, eval(x.replace('\"', '')) )

def add_image_options( p, func ):
    p.add_argument('--charge-gain', type=eval, default = '7.25', help = 'factor to convert charges into ADU' )
    p.add_argument('--readout-noise', type=eval, default = '0', help = 'sigma of the normal noise distribution in ADU' )
    p.add_argument('--dark-current', type=eval, default = '0', help = 'lambda of Poisson distribution dimensionless' )
    p.add_argument('--expose-hours', type=float, default = '1', help = 'number of images to be generated' )
    p.add_argument('--output-fits', type=str, default = None, help = 'set to generate a fits output' )
    p.add_argument('--sim', action='store_true', help = 'generate csv output' )
    p.add_argument('--pdf', action='store_true', help = 'generate pdf output' )
    p.add_argument('--spectrum', action='store_true', help = 'generate energy spectrum' )
    add_general_options( p )
    p.set_defaults( func=func )

def add_folder_options( p, func ):
    p.add_argument('--number-of-images', type=int, default = None, help = 'number of images to be generated' )
    p.add_argument('--output-pattern', type=str, default = 'simulation_*.fits', help = 'pattern for output file names' )
    p.add_argument('--readout-noise-range', type=tuple_of(float), default = '\"[11,14]\"', help = 'readout_noise range' )
    p.add_argument('--dark-current-range', type=tuple_of(float), default = '\"[0.01,0.3]\"', help = 'dark current range' )
    add_general_options( p )
    p.set_defaults( func=generate_folder )

def add_general_options( p ):
    p.add_argument('--image-mode', type=str, default = 'none', help = 'set to "1" to use official 1x1 image geomtry or "5" to 1x5' )
    p.add_argument('--number-of-charges', type=int, default = '4000', help = 'number of charges to be randomly generated' )
    p.add_argument('--charge-range', type=tuple_of(int), default = '\"[5,200]\"', help = 'range into which to randomly generate charges' )
    p.add_argument('--number-of-Cu-charges',
                        type=int,
                        default = '0',
                        help = 'number of charges to be randomly generated at the Copper fluorescence energy 8.046keV' 
                        )
    p.add_argument('--number-of-Cu2-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the secundary Copper fluorescence energy 8.904keV' 
                        )
    p.add_argument('--number-of-Si-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the Silicon fluorescence energy 1.740keV' 
                        )
    p.add_argument('--horizontal-overscan', type=int, default = '150', help = 'size of the horizontal overscan in pixels' )
    p.add_argument('--vertical-overscan', type=int, default = '90', help = 'size of the vertical overscan in pixels' )
    p.add_argument('--xyshape', type=tuple_of(int), default = '\"[4000,4000]\"', help = 'shape of the image as 2d pixels' )
    p.add_argument('--rebin', type=tuple_of(int), default = '\"[1,1]\"', help = '2d rebinning strides' )
    p.add_argument('--depth-range', type=tuple_of(int), default = '\"[0,670]\"', help = 'range into which to randomly generate depths' )
    p.add_argument('--diffusion-function',
                        type=str, 
                        default = default_diffusion_function,
                        help = 'function to map z-depth into sigma' 
                        )
    p.add_argument('--charge-efficiency-function',
                        type=str,
                        default = default_charge_efficiency_function,
                        help = 'function to map z-depth into sigma' 
                        )
    p.add_argument('--vertical-modulation-function',
                        type=str, 
                        default = default_vertical_modulation_function,
                        help = 'function to modulate the vertical axis' 
                        )    
    p.add_argument('--horizontal-modulation-function',
                        type=str, 
                        default = default_horizontal_modulation_function,
                        help = 'function to modulate the horizontal axis' 
                        )
    p.add_argument('--no-vertical-modulation', action='store_true', help = 'set vertical modulation to "0"' )
    p.add_argument('--no-horizontal-modulation', action='store_true', help = 'set horizontal modulation to "0"' )
    p.add_argument('--no-modulation', action='store_true', help = 'set modulations to "0"' )

def postprocess( args ):
    if args.image_mode is '1':
        args.rebin = [1,1]
        args.horizontal_overscan = 150
        args.vertical_overscan = 90
        args.xyshape = [4150,4120]
    elif args.image_mode is '5':
        args.rebin = [5,1]
        args.horizontal_overscan = 450
        args.vertical_overscan = 70
        args.xyshape = [4150,4120]
    
    del args.image_mode
    
    if args.no_vertical_modulation:
        args.vertical_modulation_function = '0'
    if args.no_horizontal_modulation: 
        args.horizontal_modulation_function = '0'
    if args.no_modulation:
        args.vertical_modulation_function = '0'
        args.horizontal_modulation_function = '0'

    del args.no_vertical_modulation
    del args.no_horizontal_modulation
    del args.no_modulation
    return args
    
if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser( description = 'simulation tools', formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    subparsers = parser.add_subparsers( help = 'major options' )

    add_image_options( subparsers.add_parser('image', help='generate image'), image )
    add_folder_options( subparsers.add_parser('folder', help='generate folder'), generate_folder )

    args = parser.parse_args()
    args = postprocess( args )
    args.func(args)
    
