from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *

folder = '../analysis_data_split'

if not os.path.exists(folder):
    os.mkdir(folder)

def makefigure( code, file, doc, height=1. ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(folder+'/'+file):
        subprocess.call( ['echo'], shell=True )
        subprocess.call( ['echo %s'%code], shell=True )
        subprocess.call( [code], shell=True )
    b = doc.figure( file, height=height, center=True )
    return ''.join([a,b])

def requiredFile( code, file, doc ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(folder+'/'+file):
        subprocess.call( ['echo'], shell=True )
        subprocess.call( ['echo %s'%code], shell=True )
        subprocess.call( [code], shell=True )
    #b = doc.figure( file, height=height, center=True )
    return a

with Timer('presentation'):
    with openBeamer( folder, 'temporal checks\\\\https://github.com/PhMota/CONNIEtools/wiki' ) as doc:

        doc.par('\n\part{''}\n\\frame{\partpage}\n')

        doc.frame('sanity check',
            'methodology',
            doc.itemize(
                r'divide on19 and on20 into 6 groups',
                r'./catalog histogram E1/gain3Peaks "/share/storage2/connie/DAna/nuCatalogs/shape_fit_lowE_hpix_ndc_qcalib_catalog_data_6*_to_*_v4.0.root" 1 E1/gain3Peaks "/share/storage2/connie/DAna/nuCatalogs/shape_fit_lowE_hpix_ndc_qcalib_catalog_data_14[6,7,9]*_to_*_v4.0.root" 1 --x-range 0 3 --binsize .1 --global-selection "ohdu==3 and E0/gain3Peaks>.05 and xMin>140 and xMax<3960 and yMin>75 and yMax<898 and sizell<.95 and hPixFlag==0"',
                r"https://docs.google.com/document/d/1tiVyRtv_Dn8BBYy1wyq2Co941F3s86kslq-ZsAMm9l0/edit",
'               r"https://connie-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=1937&filename=SpectrumOFF_EnergyBalanceCut_ohdu3.pdf&version=1"
            )
        )
        ohdu_selection = '(ohdu<6 or ohdu==8 or ohdu==9 or ohdu==13 or ohdu==14)'
        geometric_selection = '(xMin>140 and xMax<3960 and yMin>75 and yMax<898 and hPixFlag==0)'
        energy_selection = 'E0/gain3Peaks>.045'
        size_selection = 'sizell<.95'

        global_selection = ' and '.join([ohdu_selection, geometric_selection, energy_selection, size_selection])

        doc.frame('on19',
            '',
            makefigure(
                './catalog histogram \
                 E1/gain3Peaks "/share/storage2/connie/DAna/nuCatalogs/shape_fit_lowE_hpix_ndc_qcalib_catalog_data_6*_to_*_v4.0.root" 1 \
                 E1/gain3Peaks "/share/storage2/connie/DAna/nuCatalogs/shape_fit_lowE_hpix_ndc_qcalib_catalog_data_14[6,7,9]*_to_*_v4.0.root" 1 \
                 --x-range 0.05 3.05 --binsize .2 \
                 --global-selection "{}" \
                 --output {}/hit_test.png'.format(global_selection, folder),
                'hist_test.png',
                doc,
                height = .5
                )
        )
