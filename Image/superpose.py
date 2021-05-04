
def superpose( args ):
    print( 'inputs', args.input_files )
    print( 'output', args.output )
    if os.path.exists(args.output):
        print( 'output already exists' )
        exit(0)
    input_files = args.input_files
    subargs = args
    subargs.func = lambda x, y: x
    images = {}
    for input_file in input_files:
        subargs.input_files = input_file
        print( 'input_file', subargs.input_files )
        result = apply_to_files( subargs )
        images[result.keys()[0]] = result.values()[0]

    new_ohdus = None
    for path, ohdus in images.items():
        print( path, ohdus )
        if new_ohdus is None:
            new_ohdus = ohdus
            continue
        for key in new_ohdus.keys():
            new_ohdus[key].all += ohdus[key].all
            print( new_ohdus[key].all, new_ohdus[key].all.shape, type(new_ohdus[key].all) )

    primary = None
    hdu_list = []
    for key, image in new_ohdus.items():
        print( key, image.header )
        header = astropy.io.fits.Header()
        for field, value in image.header.items():
            header[field] = value

        hdu_list += [astropy.io.fits.ImageHDU( image.all, header=header )]

        if primary is None:
            primary_header = header
            del primary_header['OHDU']
            primary = astropy.io.fits.PrimaryHDU( header = header )

    #print( [primary] + hdu_list )
    HDU = astropy.io.fits.HDUList([primary] + hdu_list )
    HDU.writeto( args.output )
    print( 'saved', args.output )
    return

def add_superpose_options( p ):
    p.add_argument( 'input_files', nargs='+' )
    p.add_argument( '--indices', nargs='*', default=None, type=int, help='indexes to be shown' )
    p.add_argument( '--include', nargs='*', default='', help='fields to be shown' )
    p.add_argument( '--exclude', nargs='*', default='', help='fields not to be shown' )
    p.add_argument( '--output', default='output', help='output name' )
    p.set_defaults( func=superpose )
    return