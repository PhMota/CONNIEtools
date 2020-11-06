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

def makefigure( code, filename, doc, height=1. ):
    fullpath = F('{{folder}}/{{filename}}').str()
    print( fullpath )

    code += F(' --output "{{fullpath}}"').str()
    a = doc.code( code, language='Bash' )
    if not os.path.exists(fullpath):
        print()
        print( colored(code, 'green' ))
        print()
        subprocess.call( ["%s" % code], shell=True )
    if not os.path.exists(fullpath):
        print()
        print( '!!!not found', fullpath )
        exit(0)
    b = doc.figure( filename, height=height, center=True )
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
        #             --no-title --no-label-branch --no-label-file --x-range 1.45 2.1 --y-range 0.02 0.06 --png \\
        #         """),
        #         'scatter_zoom_noise_dc_a.png', doc, height = .4 )
        # )

        doc.frame('time evolution on19 and off19 0.05--1.05keV',
            makefigure(
                F('./catalog histogram \
                -s runID "{{nuCatalogs}}" "{{off19}}" "OFF19"\
                -s runID "{{nuCatalogs}}" "{{on19}}" "ON19"\
                --global-selection "{{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}}" \
                --hide-zeros --no-title --no-label-file --no-label-branch --x-range 6000 7600 --binsize 1 --png'),
                'event_time_evo_2.png', doc, height = .6 ),
        )

        doc.frame('histogram on19 and off19 0.05--1.05keV',
            makefigure(
                F('./catalog histogram \
                -s runID "{{nuCatalogs}}" "{{off19}}" "OFF19"\
                -s runID "{{nuCatalogs}}" "{{on19}}" "ON19"\
                --global-selection "{{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}}" \
                --hide-zeros --no-title --no-label-file --no-label-branch --x-range 6000 7600 --binsize 1 --png --count-histogram --extra-binsize 1'),
                'event_time_evo_count4.png', doc, height = .6 ),
        )


        for ohdu in ohdus19:
            this_ohdu_selection = 'ohdu==%s' % ohdu
            this_global_selection = r' and '.join([this_ohdu_selection, geometric_selection, energy_selection, size_selection])

            doc.frame(F('event time evolution on19 and off19 0.05--1.05keV ohdu{{ohdu}}'),
                makefigure(
                    F('./catalog histogram \
                        -s runID "{{nuCatalogs}}" "{{off19}}" "OFF19"\
                        -s runID "{{nuCatalogs}}" "{{on19}}" "ON19"\
                        --global-selection "{{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}} and ohdu=={{ohdu}}" \
                        --hide-zeros --no-title --no-label-file --x-range 6000 7600 --binsize 1 --png'),
                    F('event_time_evo_ohdu{{ohdu}}_2.png'),
                    doc,
                    height = .6
                    )
            )

            doc.frame(F('event time evolution on19 and off19 0.05--1.05keV ohdu{{ohdu}}'),
                makefigure(
                    F('./catalog histogram \
                        -s runID "{{nuCatalogs}}" "{{off19}}" "OFF19"\
                        -s runID "{{nuCatalogs}}" "{{on19}}" "ON19"\
                        --global-selection "{{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}} and ohdu=={{ohdu}}" \
                        --hide-zeros --no-title --no-label-file --x-range 6000 7600 --binsize 1 --count-histogram --extra-binsize 1'),
                    F('event_time_evo_count_ohdu{{ohdu}}_2.png'),
                    doc,
                    height = .6
                    )
            )

#             doc.frame('off19 and off20 ohdu=%s' % ohdu,
#                 '',
#                 makefigure(
# """
# ./catalog histogram \\
#     E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
#     E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
#     --global-selection "{}" \\
#     --x-range 0.05 3.05 --binsize .2 --output {}/hist_test_ohdu{} --png
# """.format(this_global_selection, folder, ohdu),
#                     'hist_test_ohdu{}.png'.format(ohdu),
#                     doc,
#                     height = .6
#                     )
#             )

#         doc.frame('noise off19 and off20',
#             '',
#             makefigure(
# """
# ./catalog scatter \\
#     scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
#     scnNoise scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
#     --global-selection "{}" \\
#     --no-title --x-range 0 7.5 --y-range 0 0.1 --output {}/scatter_noise_dc --png
# """.format(global_selection, folder),
#                 'scatter_noise_dc.png',
#                 doc,
#                 height = .6
#                 )
#         )

    # scnNoise "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\

#         doc.frame('event time evolution off19',
#             '',
#             makefigure(
# """\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
#     --global-selection "{}" \\
#     --no-title --no-label-file --x-range 6300 6950 --binsize 1 --output {}/event_time_evo --png
# """.format(global_selection, folder),
#                 'event_time_evo.png',
#                 doc,
#                 height = .6
#                 )
#         )

#         doc.frame('event time evolution off19',
#             '',
#             makefigure(
# """\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*_d*_6*_*.root" 1 \\
#     --global-selection "{} and {} and {}" \\
#     --no-title --no-label-file --x-range 6200 7000 --binsize 1 --output {}/event_time_evo_hpix_a --png
# """.format(global_selection_hpix, off19, runID_excluded, folder),
#                 'event_time_evo_hpix_a.png',
#                 doc,
#                 height = .6
#                 )
#         )

#         doc.frame('event time evolution on19 and off19 50--1050eV',
#             '',
#             makefigure(
# F("""\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*_d*_[6,7]*_*.root" "{{off19}}" \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*_d*_[6,7]*_*.root" "{{on19}}" \\
#     --global-selection "{{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.05,1.05)}}" \\
#     --hide-zeros --no-title --no-label-file --x-range 6000 7600 --binsize 1 --output {{folder}}/event_time_evo_hpix_on_excl_no0_a --png
# """),#.format(off19, on19, global_selection_hpix, runID_excluded, energy_selection_function(0.05, 1.05), folder),
#                 'event_time_evo_hpix_on_excl_no0_a.png',
#                 doc,
#                 height = .6
#                 )
#         )

#         doc.frame('event time evolution on19 and off19 0.5--1keV',
#             '',
#             makefigure(
# F("""\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*_d*_[6,7]*_*.root" "{{off19}}" \\
#     --global-selection "{{global_selection_hpix}} and {{runID_excluded}} and {{energy_selection_function(0.5,1.)}}" \\
#     --hide-zeros --no-title --no-label-file --x-range 6000 7600 --binsize 1 --output {{folder}}/event_time_evo_hpix_on_excl_a --png
# """),#.format(off19, on19, global_selection_hpix, runID_excluded, energy_selection_function(0.05, 1.05), folder),
#                 'event_time_evo_hpix_on_excl_a.png',
#                 doc,
#                 height = .6
#                 )
#         )

#         doc.frame('event time evolution on20 and off20',
#             '',
#             makefigure(
# """\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*qcalib*_d*_[13,14,15]*_*.root" "{}" \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*qcalib*_d*_[13,14,15]*_*.root" "{}" \\
#     --global-selection "{} and {}" \\
#     --hide-zeros --no-title --no-label-file --x-range 13000 15200 --binsize 1 --output {}/event_time_evo_hpix_on_off20_no0excl5a --png
# """.format(off20, on20, global_selection_hpix20, runID_excluded20, folder),
#                 'event_time_evo_hpix_on_off20_no0excl5a.png',
#                 doc,
#                 height = .6
#                 )
#         )

#         doc.frame('event time evolution off20',
#             '',
#             makefigure(
# """\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*qcalib*_d*_[13,14,15]*_*.root" "{}" \\
#     --global-selection "{} and {}" \\
#     --hide-zeros --no-title --no-label-file --x-range 14700 15200 --binsize 1 --output {}/event_time_evo_hpix_off20_no0excl5a --png
# """.format(off20, global_selection_hpix20, runID_excluded20, folder),
#                 'event_time_evo_hpix_off20_no0excl5a.png',
#                 doc,
#                 height = .6
#                 )
#         )
#
#         doc.frame('event time evolution off20',
#             '',
#             makefigure(
# """\
# ./catalog histogram \\
#     runID "/sh*/s*2/co*/D*/Cat*/hpix*qcalib*_d*_[13,14,15]*_*.root" "{}" \\
#     --global-selection "{} and {}" \\
#     --hide-zeros --no-title --no-label-file --x-range 14700 15200 --binsize 1 --output {}/event_time_evo_hpix_off20_no0excl5b --png
# """.format(off20ex, global_selection_hpix20, runID_excluded20, folder),
#                 'event_time_evo_hpix_off20_no0excl5b.png',
#                 doc,
#                 height = .6
#                 )
#         )
#
#         doc.frame('noise off19 and off20',
#             '',
#             makefigure(
# """
# ./catalog histogram \\
#     scnNoise "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
#     scnNoise "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
#     --global-selection "{}" \\
#     --x-range 0 7.5 --binsize .2 --output {}/hist_noise --png
# """.format(global_selection, folder),
#                 'hist_noise.png',
#                 doc,
#                 height = .6
#                 )
#         )
#
#         doc.frame('dc off19 and off20 fancyselection',
#             '',
#             makefigure(
# """
# ./catalog histogram \\
#     scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
#     scnDC "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
#     --global-selection "({}) and Entry\$%2==0 " \\
#     --x-range 0 .1 --binsize .01 --output {}/hist_dc_fancysel --png
# """.format(global_selection, folder),
#                 'hist_dc_fancysel.png',
#                 doc,
#                 height = .6
#                 )
#         )

#     #E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_14[6,7,9]*_*.root" 1 \\
#         doc.frame('off19 and off20',
#             '',
#             makefigure(
# """
# ./catalog histogram \\
#     E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_6*_*.root" 1 \\
#     E1/gain3Peaks "/sh*/s*2/co*/D*/nu*/sh*_d*_14[7,9]*_*.root" 1 \\
#     --global-selection "{}" \\
#     --x-range 0.05 3.05 --binsize .2 --output {}/hist_test --png
# """.format(global_selection, folder),
#                 'hist_test.png',
#                 doc,
#                 height = .6
#                 )
#         )
