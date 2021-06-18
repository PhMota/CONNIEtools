
def scatter( **args ):
    args = Namespace(**args)
    import matplotlib.pylab as plt
    with Timer('scatter'):
        print( 'len', len(args.root_file) )

        if len(args.root_file) == 1:
            args.root_file = args.root_file[0]

            file = glob.glob(args.root_file)

            args.xbranches = []
            args.ybranches = []
            args.cbranches = []
            args.selections = []
            has_color = False

            print( args.branch_selections )
            for branch_selection in args.branch_selections:
                args.xbranches.append( branch_selection[0] )
                args.ybranches.append( branch_selection[1] )
                if len(branch_selection) == 4:
                    args.cbranches.append( branch_selection[2] )
                    has_color = True
                args.selections.append( branch_selection[-1] )
            if not has_color:
                args.cbranches = args.xbranches

            data_selection = get_selections(
                                file[0],
                                args.branches,
                                args.selections,
                                args.global_selection, runID_range = args.runID_range
            )

            data_selection = get_selections( file[0],
                                        args.xbranches + args.ybranches + args.cbranches,
                                        args.selections,
                                        args.global_selection,
                                        runID_range=args.runID_range,
            )

        else:
            nargs = len(args.root_file)
            print( args.root_file )

            args.xbranches = []
            args.ybranches = []
            args.cbranches = []
            args.selections = []
            args.sel_expr = []
            args.files = []
            has_color = False

            data_selection = {}
            for i in range(nargs/4):
                xbranch = args.root_file[4*i+0]
                args.xbranches.append(xbranch)
                ybranch = args.root_file[4*i+1]
                args.ybranches.append(ybranch)

                file = args.root_file[4*i+2]
                args.files.append(file)
                selection = args.root_file[4*i+3]
                args.sel_expr.append(selection)

                key = 'xbranch: {}, ybranch: {}\nfiles: {}\nselections: {}'.format(xbranch, ybranch, file, selection)

                args.selections.append(key)
                print( 'file', file )
                print( 'br', xbranch, ybranch )
                print( 'sel', selection )

                for f in glob.glob(file):
                    print( 'file', f )
                    data_entry = get_selections( f, [xbranch,ybranch], [selection], args.global_selection )

                    try:
                        data_selection[key].append( data_entry.values()[0] )
                    except KeyError:
                        data_selection[key] = [ data_entry.values()[0] ]
            if not has_color:
                args.cbranches = args.xbranches

            print( 'selections', data_selection.keys() )
            print( 'branches', args.xbranches, args.ybranches )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        title = '{} vs {}'.format(args.xbranches[0], args.ybranches[0])
        if 'global_selection' in args:
            subtitle = args.global_selection
            if len(subtitle) > 50:
                new_subtitle = subtitle[:50]
                for i in range(1, len(subtitle)/50+1):
                    new_subtitle += '\n' + subtitle[i*50:(i+1)*50]
            subtitle = new_subtitle
            title += '\n' + subtitle
        if not 'no_title' in args:
            ax.set_title(title)

        scatter_obj = []

        if 'selections' in args:
            markers = ['.', '+', 'x', '^']
            for index, (xbranch, ybranch, cbranches, files, sel_expr, selection) in enumerate( zip(args.xbranches, args.ybranches, args.cbranches, args.files, args.sel_expr, args.selections) ):
                print( 'plot', xbranch, ybranch, cbranches )
                x_data = None
                for datum in data_selection[selection]:
                    if x_data is None:
                        x_data = datum[xbranch]
                        y_data = datum[ybranch]
                    else:
                        x_data = concatenate((x_data, datum[xbranch]))
                        y_data = concatenate((y_data, datum[ybranch]))

                if has_color:
                    colors = datum[cbranches]
                    cmap = matplotlib.cm.plasma
                    alpha = 1
                else:
                    colors = None
                    cmap = None
                    alpha = (len(datum))**(-.1) if len(datum) > 0 else 1

                label = ''
                if not 'no_label_branch' in args:
                    label += 'x: {}, y: {}\n'.format(xbranch, ybranch)
                if not 'no_label_file' in args:
                    label += 'files: {}\n'.format(files)
                if not 'no_label_selection' in args:
                    label += 'sel: {}'.format(sel_expr)
                label += ' ({})'.format(x_data.size)

                scatter_obj.append( ax.scatter( x_data, y_data, label = label, marker=markers[index%4], c=colors, alpha=alpha, cmap=cmap ) )

                if 'errorbar' in args:
                    bins = arange( args.x_range[0], args.x_range[1], 20 )
                    bin_means, bin_edges, binnumber = binned_statistic_fast( x, y, statistic='mean', bins=bins )
                    xbins = .5*(bin_edges[1:] + bin_edges[:-1])
                    dx = bin_edges[1] - bin_edges[0]
                    yerr = [0]*len(bin_means)
                    bin_std, bin_edges, binnumber = binned_statistic_fast( x, y, statistic='std', bins=bin_edges )
                    yerr = bin_std
                    ax.errorbar( xbins, bin_means, xerr=dx/2, yerr=yerr, fmt='.' )

        ax.legend(frameon=False)
        ax.grid()
        ax.set_xlim( *args.x_range )
        ax.set_ylim( *args.y_range )
        ax.set_xlabel( args.xbranches[0] )
        ax.set_ylabel( args.ybranches[0] )
        if has_color:
            fig.colorbar( scatter_obj[0] )
    if 'ylog' in args:
        ax.set_yscale('log')

    if 'output' in args:
        if args.output == '':
            args.output = args.root_file
        fig.savefig( args.output, bbox_inches='tight' )
    else:
        # args.output = args.root_file
        plt.show()

    # if 'pdf' in args:
    #     # extra = ''
    #     # fname = args.output+'.scatter.{}.vs.{}{}.pdf'.format(args.ybranches[0].replace('/','_'), args.xbranches[0].replace('/','_'), extra)
    #     fname = args.output+'.png'
    #     fig.savefig( fname )
    #     print( 'saved', fname )
    # elif 'png' in args:
    #     fname = args.output+'.png'
    #     fig.savefig( fname, bbox_inches='tight' )
    #     print( 'saved', fname )
    # else:
    #     plt.show()
    return

def add_scatter_options(p):
    p.add_argument('root_file', nargs='+', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('-s', '--branch-selections', action='append', nargs='+', type=str, default=argparse.SUPPRESS, help = 'branches used for x- and y-axis' )
    p.add_argument('--global-selection', type=str, default='1', help = 'global selection' )
    #p.add_argument('--selections', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--define', nargs='+', type=str, default=argparse.SUPPRESS, help = 'definitions' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    p.add_argument('--x-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of x' )
    p.add_argument('--y-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of y' )
    p.add_argument('--c-range', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'range of color' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'output to file' )
    p.add_argument('--pdf', action='store_true', default=argparse.SUPPRESS, help = 'output to pdf' )
    p.add_argument('--png', action='store_true', default=argparse.SUPPRESS, help = 'output to png' )
    p.add_argument('--fit', nargs='+', type=str, default=argparse.SUPPRESS, help = 'include fit column' )
    p.add_argument('--noise', nargs='+', type=str, default=argparse.SUPPRESS, help = 'include fit column' )
    p.add_argument('--errorbar', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )
    p.add_argument('--no-title', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )
    p.add_argument('--no-label-branch', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )
    p.add_argument('--no-label-file', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )
    p.add_argument('--no-label-selection', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )
    p.set_defaults(_func=scatter)
