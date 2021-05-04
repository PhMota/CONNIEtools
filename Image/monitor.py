

def add_monitor_options( p, func ):
    p.add_argument( '--runID', nargs='*', type=int, default = None, help = 'runIDs to be analysed' )
    p.add_argument( '--ohdu', nargs='*', type=int, default = [2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14], help = 'ohdus to be analysed' )
    p.add_argument( '--exclude', nargs='*', type=int, default = None, help = 'ohdus not to be analysed' )
    p.set_defaults( func=func )
    return


def monitor( args ):
    import ConnieImage
    print( args.runID, args.ohdu )
    fmt = lambda x: {int:'{:5}', float:'{: .3e}', Section:'{: .3e}'}[type(x)]
    sfmt = lambda x: {int:'{:5.5}', float:'{:10.10}', Section:'{:10.10}'}[type(x)]

    first = True
    for _runID in args.runID:
        for _ohdu in args.ohdu:
            params = []
            params.append([ 'runID', _runID ])
            params.append(['ohdu', _ohdu ])
            params.append(['rn', float(ConnieImage.FullImage( runID=_runID, ohdu=_ohdu, imageType='raw' ).readoutNoiseEstimate()) ])
            dc = float( ConnieImage.FullImage( runID=_runID, ohdu=_ohdu, imageType='raw' ).darkCurrentEstimate() )
            params.append(['dc', dc ])
            params.append(['g_lambda', dc**2 ])
            if first:
                columns = zip(*params)[0]
                max_length = max( map(len, columns) )
                print( ' '.join( [ colored( sfmt(v).format(key), attrs=['bold']) for key, v in params ] ) )
                first = False
            print( ' '.join( [ fmt(v).format(v) for keys, v in params ] ) )

