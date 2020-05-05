import sys
import astropy.io.fits
import numpy as np
import scipy.ndimage

def read_fits(fits_file):
    return astropy.io.fits.open(fits_file)

def extract_ohdu( fits, ohdu ):
    for image in fits:
        if image.header['OHDU'] == int(ohdu):
            return image.data
    raise Exception('ohdu %s not found' % ohdu)

def label_clusters( condition, edge=1 ):
    s33 = [[edge, 1, edge], [1,1,1], [edge,1,edge]]
    return scipy.ndimage.label( condition, structure=s33 )[0]

def add_border_to_clusters( clusters, border ):
    distance_to_cluster = scipy.ndimage.distance_transform_edt( clusters == 0 )
    return label_clusters( distance_to_cluster <= border )

def index_to_xy( index, shape ):
    return np.unravel_index( index, shape )

class Properties:
    def __init__( self, shape ):
        self.shape = shape
    
    def n( self, values, indices ):
        return len(values)

    def E( self, values, indices ):
        return np.sum(values)
    
    def bary( self, values, indices ):
        r = index_to_xy( indices, self.shape )
        return np.average( r, weights=values, axis=1 )
    
    def xbary( self, values, indices ):
        return self.bary(values, indices)[0]

    def ybary( self, values, indices ):
        return self.bary(values, indices)[1]
    
    def compute( self, props, val, index ):
        return [ Properties.__dict__[prop](self, val, index) for prop in props ]
    
def apply_to_clusters( image, labeled_clusters, props ):
    properties = Properties(image.shape)
    result_per_cluster = scipy.ndimage.labeled_comprehension(
        image,
        labeled_clusters,
        index = np.unique(labeled_clusters),
        func = lambda v, p: properties.compute( props, v, p ),
        out_dtype = list,
        default = 0,
        pass_positions=True,
        )
    result_per_property = zip(*result_per_cluster)
    return { p: r for p, r in zip(props, result_per_property) }

def write_to_file( file, data, keys ):
    np.savetxt( file, zip(*[ data[key] for key in keys ]), delimiter = ', ', header = ', '.join(keys) )

def check_arg( index, type, err_msg ):
    try: 
        return type(sys.argv[index])
    except: 
        print(err_msg)
        exit(1)
    
if __name__ == '__main__':
    if len(sys.argv) < 7:
        print( 'usage: \n\t> %s <fits_file> <ohdu> <threshold> <border> <quantities> <output_file>' % sys.argv[0] )
        print( 'example: \n\t> %s /share/storage2/connie/data_analysis/processed02_data/runs/043C/data_6622_to_6774/scn/images/scn_mbs_osi_runID_043_06661_Int-400_Exp-3600_12May19_07\:45_to_12May19_08\:49_p1.fits 6 60 3 "E, xbary, ybary" output.csv' % sys.argv[0] )
        exit(1)
    fits_file = sys.argv[1]
    ohdu = check_arg(2, int, err_msg = '2nd parameter <ohdu> mut be integer')
    threshold = check_arg(3, float, err_msg = '3rd parameter <threshold> must be float' )
    border = check_arg(4, int, err_msg = '4th parameter <border> must be int' )
    quantities = [ entry.strip() for entry in sys.argv[5].split(',') ]
    output_file = sys.argv[6]
    
    image = extract_ohdu( read_fits( fits_file ), ohdu )
    print( 'extracted image with shape (%d, %d)' % image.shape )
    labeled_clusters = add_border_to_clusters( label_clusters( image > threshold ), border )
    print( '%d clusters found, cluster 0 is the background' % len(np.unique(labeled_clusters)) )
    data = apply_to_clusters( image, labeled_clusters, quantities )
    print( 'computed quantites %s for each cluster' % quantities )
    write_to_file( output_file, data, quantities )
    print( 'results written to %s' % output_file )
