import argparse
try:
    from utils.utils import Namespace
except ImportError:
    from utils import Namespace

def doc2argparse( parser, func ):
    docstr = func.__doc__
    old_docstr = None
    while docstr != old_docstr:
        old_docstr = docstr
        docstr = docstr.replace('\n ', '\n')
    
    help_str = docstr.split('\n\n')[0].replace('\n', '')
    sp = parser.add_parser( func.__name__, help = help_str )
    
    args_str = docstr.split('Arguments:\n')[1].split('\n\n')[0].split('\n')
    for arg in args_str:
        arg = arg.strip()
        if arg == "": continue
        kwargs = {}
        def_str, descr = arg.split( ' '*4 )
        kwargs["help"] = descr
        def_list = def_str.split('<')
        names = def_list[0]
        if len(def_list) == 1:
            kwargs["action"] = "store_true"
        elif len(def_list) == 2:
            args_list = def_list[1].split(',')
            if '...' in args_list[-1]:
                nargs = '+'
            else:
                nargs = len(args_list)
            if nargs != 1:
                kwargs["nargs"] = nargs
        names_split = list(map( lambda s: s.strip(), names.split(',')))
        if names_split[0].startswith('-'):
            kwargs["required"] = True
#         print( "name:", names, kwargs )
        sp.add_argument( *names_split, **kwargs )
    
    if 'Options:' in docstr:
        opts_str = docstr.split('Options:\n')[1].split('\n\n')[0].split('\n')
        for opt in opts_str:
            opt = opt.strip()
            if opt == "": continue
            kwargs = {}
            kwargs["default"] = argparse.SUPPRESS
            def_str, descr = opt.split( ' '*4 )
            kwargs["help"] = descr
            def_list = def_str.split('<')
            names = def_list[0]
            if len(def_list) == 1:
                kwargs["action"] = "store_true"
            elif len(def_list) == 2:
                args_list = def_list[1].split(',')
                if '...' in args_list[-1]:
                    nargs = '+'
                else:
                    nargs = len(args_list)
                if nargs != 1:
                    kwargs["nargs"] = nargs
            names_split = list(map( lambda s: s.strip(), names.split(',')))
    #         print( "option:", names_split, kwargs )
            sp.add_argument( *names_split, **kwargs )
    
    sp.set_defaults( __call__=func )
    
if __name__ == "__main__":
    def test( **args ):
        """
        function docstring
        
        Arguments:
            first <first>    1st description
            second <second>    2nd description
            --third    3rd descr
            
        Options:
            --flag    flag description
            -r, --recursive    more description
            --extra <extra>    extra arg
            --args <arg1, arg2>    list args
            --list <...>    list
            --elist <l0, ...>    elist
        """
        args = Namespace(**args)
        print( vars(args) )
        return
    
    parser = argparse.ArgumentParser(
        description = 'parser description',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser = parser.add_subparsers( help = 'major options' )
    doc2argparse( subparser, test )
    args = parser.parse_args()
    args.__call__( **vars(args) )
