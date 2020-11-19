from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *
from fstring import F

folder = '../analysis_data_split'

print( F("folder") )
print( F("{{folder}}") )

if not os.path.exists(folder):
    os.mkdir(folder)

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
        ohdu_selection20 = '(ohdu<5 || ohdu==9 || ohdu==13 || ohdu==14)'
        geometric_selection = '(xMin>140 and xMax<3960 and yMin>75 and yMax<898 and hPixFlag==0)'
        runID_selection = '(runID>6226 and runID<6975) || (runID>14712 and runID<15178)'

        runID_excluded = '(runID != 6415 and runID != 6475 and runID != 6499 and runID != 6926 and runID != 6927 and runID != 6031 and runID != 6093 and runID != 6096 and runID != 6139 and runID != 7074 and runID != 7222 and runID != 7226 and runID != 7374)'
        runID_excluded20 = '(runID != 14729 and runID != 14816 and runID != 13125 and runID != 13126 and runID != 13315 and runID != 13628 and runID != 13911 and runID != 14069 and runID != 14212 and runID != 14661)'

        energy_selection = 'E0/gain3Peaks>.045 and E1/gain3Peaks>0.05 and E1/gain3Peaks<1.05'
        energy_selection_function = lambda minE, maxE: 'E0/gain3Peaks>.045 and E1/gain3Peaks>{} and E1/gain3Peaks<{}'.format(minE, maxE)
        size_selection = 'sizell<.95 and sizell>0.'

        off19 = "(runID>6226 and runID<6975)"
        on19 = "((runID>6030 and runID<6227) || (runID>6974 and runID<7522))"

        off20ex = "(runID>14712 and runID<15178)"
        off20 = "(runID>14712 and runID<15128)"
        on20 = " ((runID>13013 and runID<13762) || (runID>13910 and runID<14311) || (runID>14510 and runID<14713))"

        global_selection = r' and '.join([ohdu_selection, geometric_selection, energy_selection, size_selection])

        global_selection_hpix = r' and '.join([ohdu_selection, geometric_selection, size_selection])

        global_selection_hpix20 = r' and '.join([ohdu_selection20, geometric_selection, energy_selection])

        nuCatalogs = '/share/storage2/connie/DAna/nuCatalogs/shape*data*.root'
        simCatalogs = '/share/storage2/connie/DAna/nuCatalogs/match*sim*.root'

        ohdus19 = [2,3,4,5,8,9,13,14]

        # doc.frame('noise off19 zoom',
        #     '',
        #     makefigure(
        #         F("""\
        #         ./catalog scatter \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==2" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==3" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==4" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==5" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==8" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==9" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==13" \\
        #             scnNoise scnDC "{{nuCatalogs}}" "ohdu==14" \\
        #             --global-selection "{{global_selection}}" \\
        #             --no-title --no-label-branch --no-label-file --x-range 1.45 2.1 --y-range 0.02 0.06  \\
        #         """),
        #         'scatter_zoom_noise_dc_a.pdf', doc, height = .4 )
        # )

        for ohdu in ['all']+ohdus19:
            ohdu_sel = F('ohdu=={{ohdu}}') if not ohdu is 'all' else 'ohdu==ohdu'
            doc.frame(F('time evolution on19 and off19 0.05--1.05keV ohdu {{ohdu}}'),
                makefigure(
                    F('./catalog histogram ')
                    + F(' -s runID "{{nuCatalogs}}" "{{off19}}" "OFF19"')
                    + F(' -s runID "{{nuCatalogs}}" "{{on19}}" "ON19"')
                    + F(' --global-selection "{{ohdu_sel}} and {{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}}"')
                    + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range 6000 7601 --binsize 1')
                    + F(' --ylabel "counts"')
                    + F(' --xlabel "runID"')
                    + (' -f "hist_pos(OFF19, bins)" "OFF 2019 ({{ len(y) }})\n$\mu={{ roundRel(mean(y),-2) }}$, $\sigma={{ roundRel(std(y),-2) }}$"')
                    + (' -f "hist_pos(ON19, bins)" "ON 2019 ({{ len(y) }})\n$\mu={{ roundRel(mean(y),-2) }}$, $\sigma={{ roundRel(std(y),-2) }}$"')
                    , 'event_time_evo.pdf', doc, height = .6, folder=folder ),
            )

            doc.frame(F('histogram on19 and off19 0.05--1.05keV ohdu {{ohdu}}'),
                makefigure(
                    F('./catalog histogram')
                    + F(' -s runID "{{nuCatalogs}}" "{{off19}}" "OFF19"')
                    + F(' -s runID "{{nuCatalogs}}" "{{on19}}" "ON19"')
                    + F(' --global-selection "{{ohdu_sel}} and {{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}}"')
                    + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range 6000 7601 --binsize 1 --extra-binsize 1')
                    + F(' --ylabel "frequency"')
                    + F(' --xlabel "number of events"')
                    + (' -f "hist_pos( hist_pos(OFF19, bins)[1], arange(0, 25, 1), norm=True, edge=True )" "OFF 2019 ({{ int(sum(hist( hist_pos(OFF19, bins)[1], arange(0, 25, 1) )[1])) }})\n$\mu={{ roundRel(wmean(x, y),-2) }}$, $\sigma={{ roundRel(wstd(x, y), -2) }}$"')
                    + (' -f "hist_pos( hist_pos(ON19, bins)[1], arange(0, 25, 1), norm=True, edge=True )" "ON 2019 ({{ int(sum(hist( hist_pos(ON19, bins)[1], arange(0, 25, 1) )[1])) }})\n$\mu={{ roundRel(wmean(x, y),-2) }}$, $\sigma={{ roundRel(wstd(x, y), -2) }}$"')
                    , 'event_time_evo_count.pdf', doc, height = .6, folder=folder),
            )
