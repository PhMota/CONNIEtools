import argparse
from utils.doc2argparse import doc2argparse
from paths.paths import raw, scn, catalog

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'get path',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser = parser.add_subparsers( help = 'major options' )
    doc2argparse( subparser, raw )
    doc2argparse( subparser, scn )
    doc2argparse( subparser, catalog )
    args = parser.parse_args()
    ret = args.__call__( **vars(args) )
    print( ret[0] )
    
    