
def simulate_events( args ):
    shape = args['shape']
    depth_range = args['depth_range']
    charge_range = args['charge_range']
    number_of_charges = args['number_of_charges']
    array_of_positions = np.random.random( number_of_charges*2 ).reshape(-1,2)*(np.array(shape)-1)
    array_of_depths = np.random.random( number_of_charges )*(depth_range[1] - depth_range[0]) + depth_range[0]
    array_of_charges = np.random.random( number_of_charges )*(charge_range[1] - charge_range[0]) + charge_range[0]
    array_of_identities = ['random']*number_of_charges
    if args['number_of_Cu_charges'] > 0:
        number_of_Cu_charges = args['number_of_Cu_charges']
        array_of_positions = np.concatenate( (array_of_positions, np.random.random( number_of_Cu_charges*2 ).reshape(-1,2)*(np.array(shape)-1) ), axis=0 )
        array_of_depths = np.append( array_of_depths, np.random.random( number_of_Cu_charges )*(depth_range[1] - depth_range[0]) + depth_range[0] )
        array_of_charges = np.append( array_of_charges, [Cu_energy_eV/electron_in_eV]*number_of_Cu_charges )
        array_of_identities.extend( ['Cu']*number_of_Cu_charges )
        print( 'Cu', array_of_charges )
    print( 'total charges created', len( array_of_identities ) )
    return Events( array_of_positions, array_of_depths, array_of_charges, array_of_identities, args )

class Events( np.recarray ):
    def __new__(cls, array_of_positions, array_of_depths, array_of_charges, array_of_identities, args ):
        dtype = [
            ('x', float),
            ('y', float),
            #('z', float),
            ]
        obj = np.array( [ tuple(pos) for pos in array_of_positions], dtype = dtype ).view(np.recarray).view(cls)
        obj = rfn.append_fields( obj, 
                                 'z',
                                 array_of_depths, 
                                 dtypes=(float), 
                                 asrecarray=True 
                                 ).view(cls)
        obj = rfn.append_fields( obj, 
                                 'q',
                                 array_of_charges, 
                                 dtypes=(int), 
                                 asrecarray=True 
                                 ).view(cls)
        obj = rfn.append_fields( obj, 
                                 'id',
                                 array_of_identities,
                                 dtypes=None,
                                 asrecarray=True 
                                 ).view(cls)
        obj.xyshape = np.array(args['shape'])
        obj.rebin = np.array(args['rebin'])
        obj.xyrebinshape = obj.xyshape/obj.rebin
        obj.charge_gain = args['charge_gain']
        obj.number_of_Cu_charges = args['number_of_Cu_charges']
        obj.number_of_Cu2_charges = args['number_of_Cu2_charges']
        obj.number_of_Si_charges = args['number_of_Si_charges']
        obj.diffusion_function = args['diffusion_function']
        obj.charge_efficiency_function = args['charge_efficiency_function']
        obj.lambda_ = args['dark_current']
        obj.sigma = args['readout_noise']
        obj.count = len( obj.q )
        print( 'ids', obj.id )
        return obj
    
    def get_charge(self):
        return ( self.charge_efficiency_function(self.z) * self.q ).astype(int)
    
    def get_xy( self ):
        return np.array( (self.x, self.y) ).T

    def get_xy_rebin( self ):
        return (np.array( (self.x, self.y) )/self.rebin[:,None]).T

    def get_E( self ):
        return self.q*self.charge_gain
    
    def get_id_code( self ):
        map_id = {'random': 1, 'Cu': 11, 'Cu2': 12, 'Si':10}
        return map( lambda i: map_id[i], self.id )

    def get_sigma( self ):
        return self.diffusion_function( self.z )
        
    def generate_image( self, output ):
        from Image import Image
        
        if events.count > 0:
            sigma_per_event = events.get_sigma()
            charge_per_event = events.get_charge()
            Q = np.sum(charge_per_event)
            xy_norm = scipy.stats.norm.rvs( size = 2*Q ).reshape( -1, 2 )
            sigma_per_charge = np.repeat( sigma_per_event, charge_per_event )
            xy_per_charge = np.repeat( events.get_xy(), charge_per_event, axis=0 )
            xy = xy_per_charge + sigma_per_charge[:,None] * xy_norm
            bins = [
                np.arange(events.xyshape[0]+1),
                np.arange(events.xyshape[1]+1)
                ]
            image = np.histogramdd( xy, bins = bins )[0] * events.charge_gain
        else:
            image = np.zeros( events.xyshape )
        print( 'shape', image.shape )
            
        if events.lambda_ > 0:
            dark_current = scipy.stats.poisson.rvs( events.lambda_, size = events.xyshape[0]*events.xyshape[1] ).reshape(events.xyshape)
            image += dark_current * events.charge_gain
        
        image = image.reshape( (events.xyrebinshape[0], events.xyrebinshape[1], -1) ).sum(axis = -1)
        print( 'shape rebin', image.shape, events.xyrebinshape )
        
        if events.sigma > 0:
            noise = scipy.stats.norm.rvs( size=events.xyrebinshape[0]*events.xyrebinshape[1] ).reshape(events.xyrebinshape)*events.sigma
            image += noise

        return Image( image )

 
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
    parser.add_argument('--number-of-Cu-charges', type=int, default = '0', help = 'number of charges to be randomly generated at the Copper fluorescence energy 8.046keV' )
    parser.add_argument('--number-of-Cu2-charges', type=int, default = '0', help = 'number of charges to be randomly generated at the secundary Copper fluorescence energy 8.904keV' )
    parser.add_argument('--number-of-Si-charges', type=int, default = '0', help = 'number of charges to be randomly generated at the Silicon fluorescence energy 1.740keV' )
    parser.add_argument('--shape', type=tuple_of_int, default = '\"[4000,4000]\"', help = 'shape of the image in pixel per pixel per µm' )
    parser.add_argument('--rebin', type=tuple_of_int, default = '\"[1,1]\"', help = 'shape of the image in pixel per pixel per µm' )
    parser.add_argument('--charge-range', type=tuple_of_int, default = '\"[5,200]\"', help = 'range into which to randomly generate charges' )
    parser.add_argument('--depth-range', type=tuple_of_int, default = '\"[0,670]\"', help = 'range into which to randomly generate depths' )
    parser.add_argument('--charge-gain', type=eval, default = '7.25', help = 'factor to convert charges into ADU' )
    parser.add_argument('--readout-noise', type=eval, default = '0', help = 'sigma of the normal noise distribution in ADU' )
    parser.add_argument('--dark-current', type=eval, default = '0', help = 'lambda of Poisson distribution dimensionless' )
    parser.add_argument('--image-fits-output', type=str, default = 'simulation.fits', help = 'set to "none" not to generate a fits output' )
    parser.add_argument('--image-energy-spectrum', type=str, default = 'image_energy_spectrum.png', help = 'set to "none" not to plot image energy spectrum' )
    parser.add_argument('--diffusion-function',
                        type=str, 
                        default = 'sqrt(-923.666*log1p(-0.000441113*z))/15 if z < 670 else 0',
                        help = 'function to map z-depth into sigma' 
                        )
    parser.add_argument('--charge-efficiency-function',
                        type=str, 
                        default = '1. if z < 670 else .9',
                        help = 'function to map z-depth into sigma' 
                        )
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    
    print( vars(args) )
    from numpy import *
    args.diffusion_function = eval( 'vectorize(lambda z: %s)' % args.diffusion_function )
    args.charge_efficiency_function = eval( 'vectorize(lambda z: %s)' % args.charge_efficiency_function )

    simulate_events( vars(args) )
