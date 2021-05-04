from glob import glob

debug = False

def rawfits(runID, part=1):
#     pattern = f'/share/storage2/connie/data/runs/*/runID_*_{runID:05d}_*p{part:1d}.fits.fz'
    pattern = '/share/storage2/connie/data/runs/*/runID_*_{runID:05d}_*p{part:1d}.fits.fz'.format(runID=runID, part=part)
    if debug:
        print( pattern )
    files = glob(pattern)
    if len(files) == 0:
        raise Exception( 'no raw fits file found for runID {runID}'.format(runID=runID) )
    return files

def scnfits(runID, part=1):
    pattern = '/share/storage2/connie/data_analysis/processed02_data/\
runs/*/data_*/scn/images/scn_mbs_osi_runID_*_{runID:05d}_Int*_p{part:1d}.fits'.format(runID=runID, part=part)
    if debug: 
        print( pattern )
    files = glob(pattern)
    if len(files) == 0:
        raise Exception( 'no scn fits file found for runID {runID}'.format(runID=runID) )
    return files
