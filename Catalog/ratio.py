
def add_ratio_options(p):
    p.add_argument('root_file', nargs='+', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('--branch-selections', action='append', nargs=3, type=str, default=argparse.SUPPRESS, help = 'selections as branch numerator and denominator' )
    p.add_argument('--global-selection', type=str, default='1', help = 'global selection' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--energy-threshold', nargs='+', type=float, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--define', type=str, default=argparse.SUPPRESS, help = 'definitions (ex.: a=E0; b=E1)' )
    p.add_argument('--x-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of the x-axis' )
    p.add_argument('--binsize', type=eval, default=argparse.SUPPRESS, help = 'binsize' )
    p.add_argument('--factor', type=eval, default=argparse.SUPPRESS, help = 'factor' )
    p.add_argument('--nbins', type=int, default=argparse.SUPPRESS, help = 'number of bins' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--pdf', action='store_true', default=argparse.SUPPRESS, help = 'output to pdf' )
    p.add_argument('--png', action='store_true', default=argparse.SUPPRESS, help = 'output to png' )

    p.set_defaults(_func=ratio)

def ratio( **args ):
    args = Namespace(**args)

    import matplotlib.pylab as plt
    with Timer('ratio'):
        print( 'len', len(args.root_file) )
        if len(args.root_file) == 1:
            args.root_file = args.root_file[0]

            file = glob.glob(args.root_file)
            #data = open_HitSummary(file[0], args.branch, selection='flag==0' )

            args.branches = []
            args.numselections = []
            args.denselections = []
            print( args.branch_selections )
            for branch_selection in args.branch_selections:
                args.branches.append( branch_selection[0] )
                args.numselections.append( branch_selection[1] )
                args.denselections.append( branch_selection[2] )

            if not 'runID_range' in args:
                args.runID_range = None
            data_numselection = get_selections( file[0], args.branches, args.numselections, args.global_selection, runID_range=args.runID_range )
            data_denselection = get_selections( file[0], args.branches, args.denselections, args.global_selection, runID_range=args.runID_range )
                                                #extra_branches=['ePix','xPix','yPix', 'E0', 'E1', 'n0', 'n1'] )

        else:
            nargs = len(args.root_file)
            files = args.root_file
            print( args.root_file )
            args.branches = []
            args.numselections = []
            args.denselections = []
            data_numselection = {}
            data_denselection = {}
            for i in range(nargs/4):
                branch = args.root_file[5*i+0]
                args.branches.append(branch)

                file0 = args.root_file[5*i+1]
                selection0 = args.root_file[5*i+2]
                args.numselections.append('{}:{}'.format(file0, selection0))

                file1 = args.root_file[5*i+3]
                selection1 = args.root_file[5*i+4]
                args.denselections.append('{}:{}'.format(file1, selection1))

                #print( 'file', file )
                #print( 'br', branch )
                #print( 'sel', selection )
                data_entry0 = get_selections( file0, [branch], [selection0], args.global_selection )
                data_numselection.update( { '{}:{}'.format(file0, data_entry0.keys()[0]): data_entry0.values()[0] } )

                data_entry1 = get_selections( file1, [branch], [selection1], args.global_selection )
                data_denselection.update( { '{}:{}'.format(file1, data_entry1.keys()[0]): data_entry1.values()[0] } )

            #print( data_selection.keys() )
            #exit(0)



        if 'x_range' in args:
            bins = arange( args.x_range[0], args.x_range[1], args.binsize )
        else:
            first_data_selection = data_numselection.values()[0][args.branches[0]]
            if 'binsize' in args:
                bins = arange( min(first_data_selection), max(first_data_selection), args.binsize )
            elif 'nbins' in args:
                bins = linspace( min(first_data_selection), max(first_data_selection), args.nbins )
            else:
                bins = linspace( min(first_data_selection), max(first_data_selection), int(sqrt(len(first_data_selection))) )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(', '.join(args.branches))
        #ax.hist(data, bins=bins, histtype='step', label='all')
        for branch, numselection, denselection in zip(args.branches, args.numselections, args.denselections):
            #print( 'selection', data_selection[selection][branch].shape, len(bins) )
            numerator = data_numselection[numselection]
            denominator = data_denselection[denselection]

            numhist, x, dx = stats.make_histogram( numerator[branch], bins )
            denhist, x, dx = stats.make_histogram( denominator[branch], bins )
            hist = numhist.astype(float)/denhist
            dy = sqrt( (sqrt(numhist)/denhist)**2 + (numhist*sqrt(denhist)/denhist**2)**2 )
            ax.errorbar( x, hist, xerr=dx/2, yerr=dy, label='{}:{}\n{}'.format(branch,numselection, denselection), fmt=' ' )
        ax.legend()
        ax.set_xlabel(args.branches[0])
        #ax.set_ylabel( r'$\frac{{dN}}{{d {}}}$'.format(args.branches[0]) )

    if 'output' in args:
        if args.output == '':
            args.output = args.root_file
    else:
        args.output = args.root_file

    if 'pdf' in args:
        fig.savefig(args.output+'.pdf')
        print( 'saved', args.output+'.pdf' )
    elif 'png' in args:
        fig.savefig(args.output+'.png')
        print( 'saved', args.output+'.png' )
    else:
        plt.show()
    return
