from utils.utils import Namespace
import uproot4

import matplotlib.pylab as plt

def histogram( **args ):
    """
    make histograms based on .root files
    
    Arguments:
        file <file>    input .root file
        tree <tree>    tree selection
        --branches <b0, ...>    branch selections
        
    Options:
        --cuts <c0, ...>    cut selections
    """
    args = Namespace( **args )
    if not 'cuts' in args: arg.cuts = None
    data = uproot4.concatenate( f"{args.file}:{args.tree}", expressions = args.branches, cuts = arg.cuts )
    print( data )
    
    
def histogram0( **args ):
    args = Namespace(**args)

    with Timer('histogram'):
        print( 'number of input args', len(args.root_file) )
        args.branches = []
        args.files = []
        args.selections = []
        args.labels = []
        args.sel_expr = []
        data = OrderedDict()
        if len(args.root_file) == 0:
            for selection in progressbar(args.selection, msg='selections'):
                branch = selection[0]
                file = selection[1]
                sel = selection[2]
                label = selection[3]
                args.branches.append(branch)
                args.files.append(file)
                args.sel_expr.append(sel)
                args.labels.append(label)

                key = label
                args.selections.append(key)
                for f in progressbar( glob.glob(file), msg='files' ):
                    data_entry = get_selections( f, [branch], [sel], args.global_selection )
                    if len(data_entry.values()[0]) == 0:
                        continue
                    # print( 'len(data_entry)', len(data_entry.values()[0]) )
                    # print( 'type', data_entry.values()[0].__class__.__name__ )
                    try:
                        data[key] = concatenate( (data[key], data_entry.values()[0][branch]) )
                    except KeyError:
                        data[key] = data_entry.values()[0][branch]
                        # print( 'data', data[key] )

        elif len(args.root_file) == 1:
            args.root_file = args.root_file[0]

            file = glob.glob(args.root_file)
            #data = open_HitSummary(file[0], args.branch, selection='flag==0' )

            args.branches = []
            args.selections = []
            print( args.branch_selections )
            for branch_selection in args.branch_selections:
                args.branches.append( branch_selection[0] )
                args.selections.append( branch_selection[1] )

            if not 'runID_range' in args:
                args.runID_range = None
            data_selection = get_selections( file[0], args.branches, args.selections, args.global_selection, runID_range=args.runID_range )
                                                #extra_branches=['ePix','xPix','yPix', 'E0', 'E1', 'n0', 'n1'] )

        if 'labels' not in args:
            args.labels = ['']*len(args.selections)

        # if 'x_range' in args:
        #     bins = arange( float(F(args.x_range[0]).str()), float(F(args.x_range[1]).str()), args.binsize )

        fig = plt.figure()
        ax = fig.add_subplot(111)

        table = Table()
        for key in data.keys():
            for function in args.function:
                function[0] = function[0].replace(key, F('data["{{key}}"]').str() )
                function[1] = function[1].replace(key, F('data["{{key}}"]').str() )
            if 'x_range' in args:
                args.x_range[0] = args.x_range[0].replace(key, F('data["{{key}}"]').str() )
                args.x_range[1] = args.x_range[1].replace(key, F('data["{{key}}"]').str() )
                args.binsize = args.binsize.replace(key, F('data["{{key}}"]').str() )

        if 'function' in args:
            for i, function in enumerate(progressbar(args.function, msg='functions')):
                print( 'args.function', function )
                exec_string = function[0]

                if 'x_range' in args:
                    if type(args.x_range[0]) is str:
                        args.x_range[0] = float(F(args.x_range[0]).str())
                    if type(args.x_range[1]) is str:
                        args.x_range[1] = float(F(args.x_range[1]).str())
                    if type(args.binsize) is str:
                        bins = arange( args.x_range[0], args.x_range[1], float(F(args.binsize).str()) )

                print( 'exec_string', exec_string )
                print( 'label', function[1] )
                try:
                    x, y, xerr, yerr = (eval( exec_string )).astuple()
                except AttributeError:
                    x, y, xerr, yerr = eval( exec_string )

                label = F(function[1]).str()
                kwargs = {
                    'fmt': 'o'
                    , 'ms': 3
                }
                if len(function) > 2:
                    print( function[2] )
                    print( eval(function[2]) )
                    kwargs.update( eval(function[2]) )
                table.append( 'x', x )
                table.append( F('{{label.replace(" ", "_")}}').str(), y )
                table.append( 'xerr', xerr if hasattr(xerr, '__iter__') else [xerr]*len(x) )
                table.append( F('{{label.replace(" ", "_")}}_err').str(), yerr if hasattr(yerr, '__iter__') else [yerr]*len(x) )

                markers, caps, bars = ax.errorbar(
                    x+i*xerr/5./len( args.function )
                    , y
                    , xerr = xerr/2.
                    , yerr = yerr
                    , label = label
                    , **kwargs
                )
                [bar.set_alpha(.2) for bar in bars]


        try:
            # legend_artist = ax.legend(frameon=False)
            if 'legend_title' in args:
                legend_artist = ax.legend(title=F(args.legend_title).str())
            else:
                legend_artist = ax.legend()
        except IndexError:
            legend_artist = None
            print( 'index error in legend')
        # legend_artist = ax.legend(loc='upper center', bbox_to_anchor=(0.5,0))
        ax.grid(alpha=.5)
        if 'xlabel' in args:
            ax.set_xlabel(args.xlabel)
        else:
            ax.set_xlabel(args.branches[0])
        if 'ylabel' in args:
            ax.set_ylabel(args.ylabel)
        else:
            ax.set_ylabel( r'$\frac{{dN}}{{d {}}}$'.format(args.branches[0]) )
        if 'log' in args:
            ax.set_yscale('log')
    if 'x_range' in args:
        ax.set_xlim(*map(float, args.x_range))
    if 'y_range' in args:
        ax.set_ylim(*map(float, args.y_range))
    if 'output' in args:
        if args.output == '':
            args.output = args.root_file
        table.save( args.output + '.csv' )
        if legend_artist != None:
            fig.savefig( args.output, bbox_extra_artists=(legend_artist,), bbox_inches='tight' )
        else:
            fig.savefig( args.output, bbox_inches='tight' )
        print( 'saved', args.output )
    else:
        args.output = args.root_file
        fig.subplots_adjust()
        plt.show()
    return

def add_histogram_options(p):
    p.add_argument('root_file', nargs='*', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
    p.add_argument('--branch-selections', action='append', nargs=2, type=str, default=argparse.SUPPRESS, help = 'selections' )

    p.add_argument('-s','--selection', action='append', nargs='+', type=str, default=argparse.SUPPRESS, help = 'selection' )

    p.add_argument('-f','--function', nargs='+', action='append', type=str, default=argparse.SUPPRESS, help = 'function' )

    p.add_argument('--global-selection', type=str, default='1', help = 'global selection' )
    p.add_argument('--runID-range', nargs=2, type=int, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--energy-threshold', nargs='+', type=float, default=argparse.SUPPRESS, help = 'range of runIDs' )
    #p.add_argument('--define', type=str, default=argparse.SUPPRESS, help = 'definitions (ex.: a=E0; b=E1)' )
    p.add_argument('--xlabel', type=str, default=argparse.SUPPRESS, help = 'xlabel' )

    p.add_argument('--ylabel', type=str, default=argparse.SUPPRESS, help = 'ylabel' )

    p.add_argument('--average', nargs=2, type=eval, default=argparse.SUPPRESS, help = 'average from min to max' )

    p.add_argument('--x-range', nargs=2, type=str, default=argparse.SUPPRESS, help = 'range of the x-axis' )
    p.add_argument('--y-range', nargs=2, type=str, default=argparse.SUPPRESS, help = 'range of the y-axis' )

    p.add_argument('--binsize', type=str, default=argparse.SUPPRESS, help = 'binsize' )
    # p.add_argument('--extra-binsize', type=eval, default=argparse.SUPPRESS, help = 'binsize2' )
    p.add_argument('--factor', type=eval, default=argparse.SUPPRESS, help = 'factor' )
    p.add_argument('--nbins', type=int, default=argparse.SUPPRESS, help = 'number of bins' )
    p.add_argument('-o', '--output', type=str, default=argparse.SUPPRESS, help = 'selection' )
    p.add_argument('--log', action='store_true', default=argparse.SUPPRESS, help = 'log' )
    # p.add_argument('--pdf', action='store_true', default=argparse.SUPPRESS, help = 'output to pdf' )
    # p.add_argument('--png', action='store_true', default=argparse.SUPPRESS, help = 'output to png' )

    p.add_argument('--no-title', action='store_true', default=argparse.SUPPRESS, help = 'add errorbar' )

    p.add_argument('--legend-title', type=str, default=argparse.SUPPRESS, help = 'legend title' )

    p.add_argument('--no-label-branch', action='store_true', default=argparse.SUPPRESS, help = 'hide branchat the label' )
    p.add_argument('--no-label-file', action='store_true', default=argparse.SUPPRESS, help = 'hide the file at the label' )
    p.add_argument('--no-label-selection', action='store_true', default=argparse.SUPPRESS, help = 'hide the selection at the label' )
    p.add_argument('--hide-zeros', action='store_true', default=argparse.SUPPRESS, help = 'hide zero values' )
    p.add_argument('--count-histogram', action='store_true', default=argparse.SUPPRESS, help = 'creat a count histogram' )


    p.set_defaults(_func=histogram)
