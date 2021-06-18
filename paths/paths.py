from glob import glob
from utils.utils import Namespace

debug = False

def raw(**kwargs):
    """
    get raw image paths
    
    Arguments:
        runID <runID>    runID number
    
    Options:
        --part <part>    image part (default:1)
    """
    args = Namespace(**kwargs)
    runID = int(args.runID)
    if not 'part' in args:
        part = 1
    else:
        part = int(args.part)
#     pattern = f'/share/storage2/connie/data/runs/*/runID_*_{runID:05d}_*p{part:1d}.fits.fz'
    pattern = '/share/storage2/connie/data/runs/*/runID_*_{runID:05d}_*p{part:1d}.fits*'.format(runID=runID, part=part)
    if debug:
        print( pattern )
    files = glob(pattern)
    if len(files) == 0:
        raise Exception( 'no raw fits file found for runID {runID}'.format(runID=runID) )
    return files

def scn(**kwargs):
    """
    get scn image paths
    
    Arguments:
        runID <runID>    runID number
    
    Options:
        --part <part>    image part (default:1)
    """
    args = Namespace(**kwargs)
    runID = int(args.runID)
    if not 'part' in args:
        part = 1
    else:
        part = int(args.part)
    
    pattern = '/share/storage2/connie/data_analysis/processed02_data/\
runs/*/data_*/scn/images/scn_mbs_osi_runID_*_{runID:05d}_Int*_p{part:1d}.fits*'.format(runID=runID, part=part)
    files = glob(pattern)
    if len(files) == 0:
        raise Exception( 'no scn fits file found for runID {runID}'.format(runID=runID) )
    return files

def catalog(**kwargs):
    """
    get catalog paths
    
    Arguments:
        runID <runID>    runID number
    
    Options:
        --part <part>    image part (default:1)
    """
    scn_path = scn(**kwargs)
    base_path = scn_path[0].split('/scn/')[0]
    args = Namespace(**kwargs)
    prefix = ''
    if 'prefix' in args:
        prefix = args.prefix
    pattern = f'{base_path}/ext/catalog/{prefix}catalog_*.root'
    files = glob(pattern)
    if len(files) == 0:
        raise Exception( f'no catalog file found for runID {arg.runID}' )
    return files
