from __future__ import print_function
import os

import argparse
# from ROOT import TFile, TTree, AddressOf
import uproot4

from utils.utils import Namespace

def status( **args ):
    """
    show structure of trees and branches of root file
    
    Arguments:
        file <file>    input .root file
        
    Options:
        -t, --tree <name>    select tree name
        -b, --branches <branch1, ...>    select branches
    """
    args = Namespace( **args )
    with uproot4.open( args.file ) as h:
        if "tree" in args:
            h[args.tree].show()
        else:
            for tree in h.keys():
                print(tree)
                h[tree].show()
                print()
    return
    
def status0( args ):
    if not os.path.exists( args.root_file ):
        print( 'file {} does not exist'.format( args.rootfile ) )
    tfile = TFile( args.root_file )
    for key in tfile.GetListOfKeys():
        tree_name = key.GetName()
        print( colored( tree_name, 'green' ) )
        tree = getattr(tfile, tree_name)
        branches = ', '.join( map( lambda _: _.GetName(), tree.GetListOfBranches() ) )
        print( branches )
        print( tree.GetEntries() )
    tfile.Close()
    if 'tree' in args:
        if 'branches' in args:
            #data = root_numpy.root2array( args.root_file, treename = args.tree, branches = args.branches ).view(recarray)
            data = open_HitSummary( args.root_file, branches = args.branches )
            for branch in args.branches:
                print( 'branch', getattr( data, branch ), len(getattr( data, branch )) )
    return

# def add_status_options( subparser ):
#     p = subparser.add_parser(
#         'status', 
#         help='get status from catalog', 
#         formatter_class = argparse.ArgumentDefaultsHelpFormatter
#     )
#     p.add_argument('root_file', type=str, help = 'root file (example: /share/storage2/connie/DAna/Catalogs/hpixP_cut_scn_osi_raw_gain_catalog_data_3165_to_3200.root)' )
#     p.add_argument('-t', '--tree', type=str, default=argparse.SUPPRESS, help = 'tree to print' )
#     p.add_argument('--branches', nargs='+', type=str, default=argparse.SUPPRESS, help = 'branch used for x-axis' )
#     p.set_defaults( _func = status )

