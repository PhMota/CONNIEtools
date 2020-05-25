# coding: utf-8
from __future__ import print_function
from numpy import *
import scipy.stats
from numpy.lib import recfunctions as rfn
import json
import astropy.io.fits
from Timer import Timer

electron_in_eV = 3.745
electron_in_keV = electron_in_eV*1e-3
Cu_energy_eV = 8046
Cu2_energy_eV = 8904
si_energy_eV = 1740


def simulate_events( args ):
    xyshape = args['xyshape']
    depth_range = args['depth_range']
    charge_range = args['charge_range']
    number_of_charges = args['number_of_charges']
    array_of_positions = random.random( number_of_charges*2 ).reshape(-1,2)*(array(xyshape)-1)
    array_of_depths = random.random( number_of_charges )*(depth_range[1] - depth_range[0]) + depth_range[0]
    array_of_charges = random.random( number_of_charges )*(charge_range[1] - charge_range[0]) + charge_range[0]
    array_of_identities = ['random']*number_of_charges
    if 'number_of_Cu_charges' in args:
        if args['number_of_Cu_charges'] > 0:
            number_of_Cu_charges = args['number_of_Cu_charges']
            array_of_positions = concatenate( (array_of_positions, random.random( number_of_Cu_charges*2 ).reshape(-1,2)*(array(xyshape)-1) ), axis=0 )
            array_of_depths = append( array_of_depths, random.random( number_of_Cu_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] )
            array_of_charges = append( array_of_charges, [Cu_energy_eV/electron_in_eV]*number_of_Cu_charges )
            array_of_identities.extend( ['Cu']*number_of_Cu_charges )
            #print( 'Cu', array_of_charges )
    print( 'total charges created', len( array_of_identities ) )
    
    if return_simulation:
        data = []
        if len(array_of_depths) > 0:
            data = zip(array_of_positions[:,0], array_of_positions[:,1], array_of_depths, array_of_charges, array_of_identities)
        return Simulation( data, args )
        
    header = str(args)
    header += '\nx, y, z, q, id'
    if len(array_of_depths) > 0:
        data = zip(array_of_positions[:,0], array_of_positions[:,1], array_of_depths, array_of_charges, array_of_identities)
        savetxt( args['simulation_output'], data, delimiter = ', ', header = header, fmt = ['%s']*4 + ['%s'] )
    else:
        open( args['simulation_output'], 'w' ).writelines( header )
    return args['simulation_output']
    #return Events( array_of_positions, array_of_depths, array_of_charges, array_of_identities, args )

def simulation_from_file( input_file ):
    args = json.loads( open( input_file ).readline().strip(' #').replace('\'','"') )
    data = genfromtxt( input_file, delimiter= ', ', names = True, skip_header = 1, dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] ).view(Simulation)
    return Simulation( data, args )
    
class Simulation( recarray ):
    
    def __new__(cls, data, args ):
        #args = json.loads( open( input_file ).readline().strip(' #').replace('\'','"') )
        
        #obj = genfromtxt( input_file, delimiter= ', ', names = True, skip_header = 1, dtype = [('x', float), ('y', float), ('z', float), ('q', float), ('id', 'S16')] ).view(Simulation)

        obj = data.view(Simulation)
        
        for key, arg in args.items():
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

    def get_xy_rebin( self ):
        if self.count == 1:
            return (array( ([self.x], [self.y]) )/self.rebin[:,None]).T
        return (array( (self.x, self.y) )/self.rebin[:,None]).T

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
        
    def generate_image( self, output = None, pdf_output = None, return_image = False ):
        #print( 'count', self.count )
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
        #print( 'shape', image.shape )
            
        if self.dark_current > 0:
            dark_current_image = scipy.stats.poisson.rvs( self.dark_current, size = self.xyshape[0]*self.xyshape[1] ).reshape(self.xyshape)
            image += dark_current_image * self.charge_gain
        
        image = image.reshape( (self.xyrebinshape[0], self.xyrebinshape[1], -1) ).sum(axis = -1)
        #print( 'shape rebin', image.shape )
        
        new_image = zeros( self.fullshape )
        new_image[self.fullshape[0]-self.xyrebinshape[0]:, :self.xyrebinshape[1]] = image
        image = new_image
        
        if self.readout_noise > 0:
            readout_noise_image = scipy.stats.norm.rvs( size=self.fullshape[0]*self.fullshape[1] ).reshape(self.fullshape)*self.readout_noise
            image += readout_noise_image
        
        print( 'image shapes', self.xyshape, self.xyrebinshape, self.fullshape )
        vertical_modulation = array([ self.vertical_modulation_function(y) for y in range(self.fullshape[0]) ])[:,None]
        #print( 'vertical_modulation', vertical_modulation.shape )
        image += vertical_modulation
        horizontal_modulation = array([ self.horizontal_modulation_function(x) for x in range(self.fullshape[1]) ])[None,:]
        #print( 'horizontal_modulation', horizontal_modulation.shape )
        image += horizontal_modulation
        
        if return_image:
            return self.save_fits( output, image, return_image = True )
        else:
            self.save_fits( output, image, return_image = False )
        if 'image_pdf_output' in self.__dict__:
            if self.image_pdf_output != 'none':
                self.save_pdf( pdf_output, image )
        return output

    def save_fits( self, output, image, return_image ):
        header = astropy.io.fits.Header()
        header['OHDU'] = 0
        header['biasW'] = self.horizontal_overscan
        header['biasH'] = self.vertical_overscan
        header['rebin'] = max(self.rebin)
        header['gain'] = self.charge_gain
        primary = astropy.io.fits.PrimaryHDU(header=header)
        header['OHDU'] = -1
        fits_image = astropy.io.fits.ImageHDU( image.astype(int), header=header )
        if return_image:
            return fits_image
        with Timer('writing ' + output ):
            hdu_list = astropy.io.fits.HDUList([primary, fits_image])
            hdu_list.writeto( output, overwrite=True )
            return
    
    def save_pdf(self, output, image ):
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
        
default_diffusion_function = 'sqrt(-923.666*log1p(-0.000441113*z))/15 if z < 670 else 0'
default_charge_efficiency_function = '1. if z < 670 else .9'
default_vertical_modulation_function = '50*cos(y/2/pi/20)'
default_horizontal_modulation_function = '-1e3*(x/1000. - 1)**2 if x < 1000. else 0'

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description = 'simulate events, generates image and extract hits',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        )
    def tuple_of_int( x ):
        return map( int, eval(x.replace('\"', '')) )
    
    parser.add_argument('number_of_charges', type=int, default = '4000', help = 'number of charges to be randomly generated' )
    parser.add_argument('--number-of-Cu-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the Copper fluorescence energy 8.046keV' 
                        )
    parser.add_argument('--number-of-Cu2-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the secundary Copper fluorescence energy 8.904keV' 
                        )
    parser.add_argument('--number-of-Si-charges', 
                        type=int, 
                        default = '0', 
                        help = 'number of charges to be randomly generated at the Silicon fluorescence energy 1.740keV' 
                        )
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

    simulate_events( vars(args) )
    simulation = Simulation( vars(args)['simulation_output'] )
    simulation.generate_image( vars(args)['image_fits_output'], vars(args)['image_pdf_output'] )
