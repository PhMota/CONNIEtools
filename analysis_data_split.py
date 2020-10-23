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
    if not os.path.exists(folder + '/' + file):
        print()
        print(code)
        print()
        subprocess.call( [code], shell=True )
    if not os.path.exists(folder + '/' + file):
        print()
        print( '!!!not found', folder + '/' + file )
        exit(0)
    b = doc.figure( file, height=height, center=True )
    return ''.join([a,b])

def requiredFile( code, file, doc ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(folder+'/'+file):
        print()
        print(code)
        print()
        subprocess.call( [code], shell=True )
    #b = doc.figure( file, height=height, center=True )
    return a

with Timer('presentation'):
    with openBeamer( folder, 'temporal checks\\\\https://github.com/PhMota/CONNIEtools/wiki' ) as doc:

        doc.par('\n\part{''}\n\\frame{\partpage}\n')

        doc.frame('sanity check',
            'methodology',
            doc.itemize(
                r'divide off19 and off20 into 6 groups',
                r"\url{https://docs.google.com/document/d/1tiVyRtv_Dn8BBYy1wyq2Co941F3s86kslq-ZsAMm9l0/edit}",
                r"\url{https://connie-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=1937&filename=SpectrumOFF_EnergyBalanceCut_ohdu3.pdf&version=1}"
            )
        )
        # ohdu_selection = 'ohdu==3'
        ohdu_selection = '(ohdu<6 || ohdu==8 || ohdu==9 || ohdu==13 || ohdu==14)'
        geometric_selection = '(xMin>140 and xMax<3960 and yMin>75 and yMax<898 and hPixFlag==0)'
        runID_selection = '(runID>6226 and runID<6975) || (runID>14712 and runID<15178)'
        energy_selection = 'E0/gain3Peaks>.045 and E1/gain3Peaks>0.05 and E1/gain3Peaks<3.05'
        size_selection = 'sizell<.95'

        global_selection = r' and '.join([ohdu_selection, geometric_selection, energy_selection, size_selection])

        doc.frame('noise off19 and off20',
            '',
            makefigure(
"""
./catalog scatter \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
    --global-selection "{}" \\
    --no-title --x-range 0 7.5 --y-range 0 0.1 --output {}/scatter_noise_dc --png
""".format(global_selection, folder),
                'scatter_noise_dc.png',
                doc,
                height = .6
                )
        )

        doc.frame('noise off19 zoom',
            '',
            makefigure(
"""\
./catalog scatter \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==2" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==3" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==4" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==5" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==8" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==9" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==13" \\
    scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" "ohdu==14" \\
    --global-selection "{}" \\
    --no-title --no-label-branch --no-label-file --x-range 1.45 2.1 --y-range 0.02 0.06 --output {}/scatter_zoom_noise_dc --png
""".format(global_selection, folder),
                'scatter_zoom_noise_dc.png',
                doc,
                height = .4
                )
        )

        doc.frame('noise off19 and off20',
            '',
            makefigure(
"""
./catalog histogram \\
    scnNoise "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
    scnNoise "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
    --global-selection "{}" \\
    --x-range 0 7.5 --binsize .2 --output {}/hist_noise --png
""".format(global_selection, folder),
                'hist_noise.png',
                doc,
                height = .6
                )
        )

        doc.frame('dc off19 and off20 fancyselection',
            '',
            makefigure(
"""
./catalog histogram \\
    scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
    scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
    --global-selection "({}) and Entry\$%2==0 " \\
    --x-range 0 .1 --binsize .01 --output {}/hist_dc_fancysel --png
""".format(global_selection, folder),
                'hist_dc_fancysel.png',
                doc,
                height = .6
                )
        )

    #E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_14[6,7,9]*_*.root" 1 \\
        doc.frame('off19 and off20',
            '',
            makefigure(
"""
./catalog histogram \\
    E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
    E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
    --global-selection "{}" \\
    --x-range 0.05 3.05 --binsize .2 --output {}/hist_test --png
""".format(global_selection, folder),
                'hist_test.png',
                doc,
                height = .6
                )
        )

        for ohdu in [2,3,4,5,6,8,9,13,14]:
            this_ohdu_selection = 'ohdu==%s' % ohdu
            this_global_selection = r' and '.join([this_ohdu_selection, geometric_selection, energy_selection, size_selection])

            doc.frame('off19 and off20 ohdu=%s' % ohdu,
                '',
                makefigure(
"""
./catalog histogram \\
    E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
    E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
    --global-selection "{}" \\
    --x-range 0.05 3.05 --binsize .2 --output {}/hist_test_ohdu{} --png
""".format(this_global_selection, folder, ohdu),
                    'hist_test_ohdu{}.png'.format(ohdu),
                    doc,
                    height = .6
                    )
            )
