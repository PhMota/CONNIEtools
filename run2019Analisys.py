from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *
from fstring import F

folder = '../run2019Analisys'

if not os.path.exists(folder):
    os.mkdir(folder)

ohdu_selection = '(ohdu<6 || ohdu==8 || ohdu==9 || ohdu==13 || ohdu==14)'
ohdu_selection20 = '(ohdu<5 || ohdu==9 || ohdu==13 || ohdu==14)'

geometric_selection = '(xMin>140 and xMax<3960 and yMin>75 and yMax<898)'
column_selection = 'hPixFlag==0'
runID_selection = '(runID>6226 and runID<6975) || (runID>14712 and runID<15178)'

runID_excluded = '(runID != 6415 and runID != 6475 and runID != 6499 and runID != 6926 and runID != 6927 and runID != 6031 and runID != 6093 and runID != 6096 and runID != 6139 and runID != 7074 and runID != 7222 and runID != 7226 and runID != 7374)'
runID_excluded20 = '(runID != 14729 and runID != 14816 and runID != 13125 and runID != 13126 and runID != 13315 and runID != 13628 and runID != 13911 and runID != 14069 and runID != 14212 and runID != 14661)'

runID_excluded_low2 = '(runID != 6134 and runID != 6981 and runID != 6492)'
runID_excluded_low3 = '(runID != 7140 and runID != 6767)'
runID_excluded_low4 = '(runID != 6135 and runID != 7006 and runID != 7278 and runID != 6503 and runID != 6884 and runID != 6916 and runID != 6953)'
runID_excluded_low5 = '(runID != 7000 and runID != 7135 and runID != 7284 and runID != 7404 and runID != 6781 and runID != 6896)'
runID_excluded_low8 = '(runID != 6362 and runID != 6406 and runID != 6857 and runID != 6145 and runID != 7520)'
runID_excluded_low9 = '(runID != 6473 and runID != 7088)'
runID_excluded_low13 = '(runID != 6327 and runID != 6390 and runID != 6484 and runID != 7412)'
runID_excluded_low14 = '(runID != 6032 and runID != 6119 and runID != 7389)'

runID_excluded_low = ' and '.join([
    runID_excluded_low2,
    runID_excluded_low3,
    runID_excluded_low4,
    runID_excluded_low5,
    runID_excluded_low8,
    runID_excluded_low9,
    runID_excluded_low13,
    runID_excluded_low14,
])

energy_selection = 'E0/gain3Peaks>.045 and E1/gain3Peaks>0.05 and E1/gain3Peaks<2.05'
energy_selection_function = lambda minE, maxE: 'E0/gain3Peaks>.045 and E1/gain3Peaks>{} and E1/gain3Peaks<{}'.format(minE, maxE)
size_selection = 'sizell<.95 and sizell>0.'

unbias_energy_selection_function = lambda minE, maxE: 'E1/gain3Peaks>{} and E1/gain3Peaks<{}'.format(minE, maxE)

off19 = "(runID>6226 and runID<6975)"
on19 = "((runID>6030 and runID<6227) || (runID>6974 and runID<7522))"

off20ex = "(runID>14712 and runID<15178)"
off20 = "(runID>14712 and runID<15128)"
on20 = " ((runID>13013 and runID<13762) || (runID>13910 and runID<14311) || (runID>14510 and runID<14713))"

noff19 = 652.
non19 = 735.

noff_noise19 = 637.
non_noise19 = 718.


nuCatalogs = '/share/storage2/connie/DAna/nuCatalogs/shape*data*.root'
simCatalogs = '/share/storage2/connie/DAna/nuCatalogs/match*sim*.root'

ohdus19 = [2,3,4,5,8,9,13,14]


M = Expr('M')
m = Expr('m')
N = Expr('N')
calN = Expr('\mathcal{N}')
calS = Expr('\mathcal{S}')
calU = Expr('\mathcal{U}')
E = Expr('E')
R = Expr(r'{\rm R}')
calR = Expr('\mathcal{R}')
T = Expr('T')
t = Expr('t')
p = Expr('p')
n = Expr('n')
i = Expr('i')
a = Expr('a')
bg = Expr(r'{\rm bg}')
density_units = Expr(r'\,{\rm g}/{\rm cm}^3')
micrometer = Expr(r'\,\mu{\rm m}')
seconds = Expr(r'\,{\rm s}')
eVolts = Expr(r'\,e{\rm V}')
d_dE = lambda x: frac( d(x), d(E) )
N_S = delim(calN[calS])[i,a]
N_U = delim(calN[calU])[i,a]
N_nu = delim( N[nu] )[i,a]
N_bg = delim( N[bg] )[i,a]
N_nubg = delim( N[nu+bg] )[i,a]
w = Expr('w')


binsizes = [.130, .02]
Emaxs = [2.05, .55]
loop = zip(binsizes, Emaxs)

histogram_command = './catalog histogram'

def load( name, branch = 'E1/gain3Peaks', catalogs = nuCatalogs, selections = [], ohdu_list = ohdus19 ):
    selections_str = ' and '.join(['ohdu=={{ohdu_}}']+selections)
    expr = ''.join( [F(' -s {{branch}} "{{catalogs}}" "{{F(selections_str)}}" "{{name}}_{{ohdu_}}"').str() for ohdu_ in ohdu_list ])
    return expr

exposure_factor = 1./(2.329/1e3*.0015**2*.0675*4120**2/24. )

def rate_func(name, label, binsize, ohdu_list, extra='', N_on = non19, N_off = noff19):
    expr = ', '.join( [ F('hist({{name}}_ON_{{ohdu_}}, bins)*{{exposure_factor}}/{{N_on}} - hist({{name}}_OFF_{{ohdu_}}, bins)*{{exposure_factor}}/{{N_off}}').str() for ohdu_ in ohdu_list ] )
    return F(r' -f "waverage( [{{expr}}] )/{{binsize}}" "{{label}}" {{extra}}')

def rate_defs(Emax, binsize, name, global_selection = runID_excluded):
    return histogram_command + F(' --global-selection "{{global_selection}}" --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "{{name}}"')

def weighted( name, label, binsize, Emax, ohdus, selections = [geometric_selection, size_selection, 'E0/gain3Peaks>.045', column_selection], extra='', N_on = non19, N_off = noff19 ):
    cmd = load(
            F('{{name}}_ON')
            , selections = [on19] + selections
            , ohdu_list = ohdus
        ) + load(
            F('{{name}}_OFF')
            , selections = [off19] + selections
            , ohdu_list = ohdus
        ) + rate_func( name, label, binsize, ohdus, extra=extra, N_on=N_on, N_off=N_off )
    return cmd

height = .45

with Timer('presentation'):
    # with openBeamer( folder, r'run2019 data analysis' ) as doc:
    #
    #     # doc += part('measured event rate')
    #
    #     doc += frame( F('weighted average of the differential rate'),
    #         'the weighted average of the differential event rate is',
    #         Math((
    #             ave(delim(calR[nu])[i,a]) == (
    #                 frac( Sum[a]* w[i,a] * delim(calR[nu])[i,a], Sum[a]* w[i,a] )
    #             ),
    #             w[i,a] == frac(1, delim(delta*delim(calR[nu])[i,a])**2)
    #             )
    #         ),
    #         column( *[ makefigure(
    #             rate_defs( Emax, binsize, 'weighted average' )
    #             + weighted( 'data', 'reactor ON--OFF', binsize, Emax, ohdus19)
    #             , F(r'differentialRateOHDUwAverageBinsize{{binsize}}.pdf')
    #             , height = height
    #             , folder = folder
    #             , nocode = True
    #             ) for binsize, Emax in loop ], widths=[.5,.5]
    #         )
    #     )
    #
    #     doc += frame( F('weighted average of the differential rate'),
    #         'the weighted average of the differential event rate is',
    #         Math((
    #             ave(delim(calR[nu])[i,a]) == (
    #                 frac( Sum[a]* w[i,a] * delim(calR[nu])[i,a], Sum[a]* w[i,a] )
    #             ),
    #             w[i,a] == frac(1, delim(delta*delim(calR[nu])[i,a])**2)
    #             )
    #         ),
    #         column( *[ makefigure(
    #             rate_defs( Emax, binsize, 'weighted average\nreactor ON--OFF', global_selection = '' )
    #             + weighted( 'data', 'low-noise included', binsize, Emax, ohdus19, selections = [geometric_selection, size_selection, 'E0/gain3Peaks>.045', column_selection, runID_excluded])
    #             + weighted( 'data2', 'low-noise excluded', binsize, Emax, ohdus19, selections = [geometric_selection, size_selection, 'E0/gain3Peaks>.045', column_selection, runID_excluded, runID_excluded_low], N_on = non_noise19, N_off = noff_noise19)
    #             , F(r'differentialRateOHDUwAverageLowNoiseXBinsize{{binsize}}.pdf')
    #             , height = height
    #             , folder = folder
    #             , nocode = True
    #             ) for binsize, Emax in loop ], widths=[.5,.5]
    #         )
    #     )
    #
    #     itemize_entry = itemize(
    #                     textcolor( 'orange', 'energy: ' ) + inlinecode( 'E0/gain3Peaks>.045' )
    #                     , textcolor( 'teal', 'geometric: ' ) + inlinecode( geometric_selection )
    #                     , textcolor( 'red', 'columns: ' ) + inlinecode( column_selection )
    #                     , textcolor( 'violet', 'depth: ' ) + inlinecode( size_selection )
    #                 )
    #     doc += frame( F('comparison of the selections'),
    #         'adding selections to the unbiased data',
    #         itemize_entry,
    #         column( *[ makefigure(
    #             rate_defs( Emax, binsize, 'weighted average\nreactor ON--OFF' )
    #             + weighted( 'dataU', 'unbiased', binsize, Emax, ohdus19, selections = [])
    #             + weighted( 'dataE', 'energy selection', binsize, Emax, ohdus19, selections = ['E0/gain3Peaks>.045'])
    #             + weighted( 'dataG', 'geometric selection', binsize, Emax, ohdus19, selections = [geometric_selection])
    #             + weighted( 'dataC', 'column selection', binsize, Emax, ohdus19, selections = [column_selection])
    #             + weighted( 'dataS', 'depth selection', binsize, Emax, ohdus19, selections = [size_selection])
    #             + weighted( 'data', 'combined', binsize, Emax, ohdus19, extra = '"{\'ms\':5}"')
    #             , F(r'differentialRateOHDUwAverageBinsize{{binsize}}.pdf')
    #             , height = height
    #             , folder = folder
    #             , nocode = True
    #             ) for binsize, Emax in loop ], widths=[.5,.5]
    #         )
    #     )
    #
    #     doc += frame( F('comparison of the selections'),
    #         'subtracting selections from the combined data',
    #         itemize_entry,
    #         column( *[ makefigure(
    #             rate_defs( Emax, binsize, 'weighted average\nreactor ON--OFF' )
    #             + weighted( 'dataU', 'unbiased', binsize, Emax, ohdus19, selections = [])
    #             + weighted( 'dataE', 'energy selection', binsize, Emax, ohdus19, selections = [geometric_selection, column_selection, size_selection])
    #             + weighted( 'dataG', 'geometric selection', binsize, Emax, ohdus19, selections = ['E0/gain3Peaks>.045', column_selection, size_selection])
    #             + weighted( 'dataC', 'column selection', binsize, Emax, ohdus19, selections = ['E0/gain3Peaks>.045', geometric_selection, size_selection])
    #             + weighted( 'dataS', 'depth selection', binsize, Emax, ohdus19, selections = ['E0/gain3Peaks>.045', geometric_selection, column_selection])
    #             + weighted( 'data', 'combined', binsize, Emax, ohdus19, selections = ['E0/gain3Peaks>.045', geometric_selection, column_selection, size_selection], extra = '"{\'ms\':5}"')
    #             , F(r'differentialRateOHDUwAverageBinsize{{binsize}}.pdf')
    #             , height = height
    #             , folder = folder
    #             , nocode = True
    #             ) for binsize, Emax in loop ], widths=[.5,.5]
    #         )
    #     )
    #
    #     for ohdu in ohdus19:
    #         doc += frame(F('differential rate ohdu {{ohdu}}'),
    #             'adding selections to the unbiased data',
    #             itemize_entry,
    #             column( *[ makefigure(
    #                 rate_defs( Emax, binsize, F('OHDU{{ohdu}}\nreactor ON--OFF') )
    #                 + weighted( 'dataU', 'unbiased', binsize, Emax, [ohdu], selections = [])
    #                 + weighted( 'dataE', 'energy selection', binsize, Emax, [ohdu], selections = ['E0/gain3Peaks>.045'])
    #                 + weighted( 'dataG', 'geometric selection', binsize, Emax, [ohdu], selections = [geometric_selection])
    #                 + weighted( 'dataC', 'column selection', binsize, Emax, [ohdu], selections = [column_selection])
    #                 + weighted( 'dataS', 'depth selection', binsize, Emax, [ohdu], selections = [size_selection])
    #                 + weighted( 'data', 'combined', binsize, Emax, [ohdu], extra = '"{\'ms\':5}"')
    #                 , F(r'differentialRateOHDU{{ohdu}}binsize{{binsize}}.pdf')
    #                 , height = height
    #                 , folder = folder
    #                 , nocode = True
    #                 ) for binsize, Emax in loop ],
    #             widths=[.5,.5]
    #             )
    #         )
    #
    #     for ohdu in ohdus19:
    #         doc += frame(F('differential rate ohdu {{ohdu}}'),
    #             'subtracting selections from the combined data',
    #             itemize_entry,
    #             column( *[ makefigure(
    #                 rate_defs( Emax, binsize, F('OHDU{{ohdu}}\nreactor ON--OFF') )
    #                 + weighted( 'dataU', 'unbiased', binsize, Emax, [ohdu], selections = [])
    #                 + weighted( 'dataE', 'energy selection', binsize, Emax, [ohdu], selections = [geometric_selection, column_selection, size_selection])
    #                 + weighted( 'dataG', 'geometric selection', binsize, Emax, [ohdu], selections = ['E0/gain3Peaks>.045', column_selection, size_selection])
    #                 + weighted( 'dataC', 'column selection', binsize, Emax, [ohdu], selections = ['E0/gain3Peaks>.045', geometric_selection, size_selection])
    #                 + weighted( 'dataS', 'depth selection', binsize, Emax, [ohdu], selections = ['E0/gain3Peaks>.045', geometric_selection, column_selection])
    #                 + weighted( 'data', 'combined', binsize, Emax, [ohdu], selections = ['E0/gain3Peaks>.045', geometric_selection, column_selection, size_selection], extra = '"{\'ms\':5}"')
    #                 , F(r'differentialRateOHDU{{ohdu}}binsize{{binsize}}.pdf')
    #                 , height = height
    #                 , folder = folder
    #                 , nocode = True
    #                 ) for binsize, Emax in loop ],
    #             widths=[.5,.5]
    #             )
    #         )

    with openBeamer( folder, r'on-chip sources' ) as doc:

        doc += frame('excluded runIDs',
            'excluded runIDs',
            itemize(
                'remove high noise (official): ' + inlinecode(runID_excluded, language='bash'),
                'remove zero-like noise: ' + inlinecode(runID_excluded_low, language='bash')
            )
        )

        doc += part('noise evolution')

        for ohdu in ohdus19:
            doc += frame(F('noise evolution OHDU{{ohdu}}'),
                column(
                makefigure(
                    './catalog histogram'
                    + F(' -s "runID" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "runID_ON_19"')
                    + F(' -s "scnNoise" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "noise_ON_19"')
                    + F(' -s "runID" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "runID_OFF_19"')
                    + F(' -s "scnNoise" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "noise_OFF_19"')
                    + F(' -f "unique_tuples(runID_ON_19, noise_ON_19, 0, 0)"') + (' "reactor ON {{len(unique(runID_ON_19))}}"')
                    + F(' -f "unique_tuples(runID_OFF_19, noise_OFF_19, 0, 0)"') + (' "reactor OFF {{len(unique(runID_OFF_19))}}"')
                    + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --xlabel "runID" --ylabel "noise [e-]" --legend-title "noise"')
                    , F(r'evoNoiseOHDU{{ohdu}}.pdf')
                    , height = height
                    , folder = folder
                    , nocode = True
                    ),
                makefigure(
                    './catalog histogram'
                    + F(' -s "runID" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "runID_ON_19"')
                    + F(' -s "scnNoise" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "noise_ON_19"')
                    + F(' -s "runID" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "runID_OFF_19"')
                    + F(' -s "scnNoise" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "noise_OFF_19"')
                    + F(' -f "hist(unique_tuples(runID_ON_19, noise_ON_19, 0, 0)[1], bins,norm=True)"') + ' "reactor ON {{len(unique(runID_ON_19))}}"'
                    + F(' -f "hist(unique_tuples(runID_OFF_19, noise_OFF_19, 0, 0)[1],bins,norm=True)"') + ' "reactor OFF {{len(unique(runID_OFF_19))}}"'
                    + ' --x-range "{{ mean(noise_OFF_19) - 5*std(noise_OFF_19) }}" "{{ mean(noise_OFF_19) + 5*std(noise_OFF_19) }}" --binsize "{{std(noise_OFF_19)/2}}"'
                    + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --ylabel "frequency" --xlabel "noise [e-]" --legend-title "noise"')
                    , F(r'histNoiseOHDU{{ohdu}}.pdf')
                    , height = height
                    , folder = folder
                    , nocode = True
                    ),
                widths = [.5,.5]
                )
                )

        doc += part('dark current evolution')

        for ohdu in ohdus19:
            doc += frame(F('dark current OHDU{{ohdu}}'),
                column(
                makefigure(
                    './catalog histogram'
                    + F(' -s "runID" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "runID_ON_19"')
                    + F(' -s "scnDC" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "dc_ON_19"')
                    + F(' -s "runID" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "runID_OFF_19"')
                    + F(' -s "scnDC" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "dc_OFF_19"')
                    + F(' -f "unique_tuples(runID_ON_19, dc_ON_19, 0, 0)" "reactor ON"')
                    + F(' -f "unique_tuples(runID_OFF_19, dc_OFF_19, 0, 0)" "reactor OFF"')
                    + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --xlabel "runID" --ylabel "dark current [e-]" --legend-title "dark current"')
                    , F(r'evoDCOHDU{{ohdu}}.pdf')
                    , height = height
                    , folder = folder
                    , nocode = True
                    ),
                makefigure(
                    './catalog histogram'
                    + F(' -s "runID" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "runID_ON_19"')
                    + F(' -s "scnDC" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "dc_ON_19"')
                    + F(' -s "runID" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "runID_OFF_19"')
                    + F(' -s "scnDC" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "dc_OFF_19"')
                    + F(' -f "hist(unique_tuples(runID_ON_19, dc_ON_19, 0, 0)[1], bins,norm=True)"') + ' "reactor ON {{len(unique(runID_ON_19))}}"'
                    + F(' -f "hist(unique_tuples(runID_OFF_19, dc_OFF_19, 0, 0)[1],bins,norm=True)"') + ' "reactor OFF {{len(unique(runID_OFF_19))}}"'
                    + ' --x-range "{{ mean(dc_OFF_19) - 5*std(dc_OFF_19) }}" "{{ mean(dc_OFF_19) + 5*std(dc_OFF_19) }}" --binsize "{{std(dc_OFF_19)/2}}"'
                    + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --ylabel "frequency" --xlabel "dark current [e-]" --legend-title "dark current"')
                    , F(r'histDCOHDU{{ohdu}}.pdf')
                    , height = height
                    , folder = folder
                    , nocode = True
                    ),
                widths = [.5,.5]
                )
                )


        # doc += frame(F('noise vs dark current OHDU{{ohdu}}'),
        #     makefigure(
        #         './catalog histogram'
        #         # + F(' -s "runID" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "runID_ON_19"')
        #         + F(' -s "scnNoise" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "noise_ON_19"')
        #         + F(' -s "scnDC" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu}}" "dc_ON_19"')
        #         # + F(' -s "runID" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "runID_OFF_19"')
        #         # + F(' -s "scnDC" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu}}" "dc_OFF_19"')
        #         + F(' -f "unique_tuples(noise_ON_19, dc_ON_19, 0, 0)" "OHDU{{ohdu}}"')
        #         # + F(' -f "unique_tuples(runID_OFF_19, dc_OFF_19, 0, 0)" "reactor OFF"')
        #         + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --xlabel "noise [e\${}^-\$]" --ylabel "dark current [e\${}^-\$]" --legend-title ""')
        #         , F(r'noiseDCOHDU{{ohdu}}.pdf')
        #         , height = height
        #         , folder = folder
        #         , nocode = True
        #         )
        #     )

        doc += part('combined')

        doc += frame(F('noise vs dark current'),
            column(
            makefigure(
                './catalog histogram'
                + ' '.join( [ F(' -s "scnNoise" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu_}}" "noise_ON_{{ohdu_}}"').str() for ohdu_ in ohdus19 ] )
                + ' '.join( [ F(' -s "scnDC" "{{nuCatalogs}}" "{{on19}} and ohdu=={{ohdu_}}" "dc_ON_{{ohdu_}}"').str() for ohdu_ in ohdus19] )
                + ' '.join( [ F(' -f "unique_tuples(noise_ON_{{ohdu_}}, dc_ON_{{ohdu_}}, 0, 0)" "OHDU {{ohdu_}}"').str() for ohdu_ in ohdus19 ] )
                + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --xlabel "noise [e-]" --ylabel "dark current [e-]" --legend-title "reactor ON" --x-range 1.5 2.1 --binsize 1 --y-range 0.02 0.06')
                , F(r'noiseDCOHDUallON.pdf')
                , height = height
                , folder = folder
                , nocode = True
                ),
            makefigure(
                './catalog histogram'
                + ' '.join( [ F(' -s "scnNoise" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu_}}" "noise_ON_{{ohdu_}}"').str() for ohdu_ in ohdus19 ] )
                + ' '.join( [ F(' -s "scnDC" "{{nuCatalogs}}" "{{off19}} and ohdu=={{ohdu_}}" "dc_ON_{{ohdu_}}"').str() for ohdu_ in ohdus19] )
                + ' '.join( [ F(' -f "unique_tuples(noise_ON_{{ohdu_}}, dc_ON_{{ohdu_}}, 0, 0)" "OHDU {{ohdu_}}"').str() for ohdu_ in ohdus19 ] )
                + F(' --global-selection "{{runID_excluded}} and {{runID_excluded_low}}" --xlabel "noise [e-]" --ylabel "dark current [e-]" --legend-title "reactor OFF" --x-range 1.5 2.1 --binsize 1 --y-range 0.02 0.06')
                , F(r'noiseDCOHDUallOFF.pdf')
                , height = height
                , folder = folder
                , nocode = True
                )
            , widths = [.5,.5]
            )
            )


# ./image display "/share/storage2/connie/data/runs/043/runID_043_07000_Int-400_Exp-3600_28May19_04:12_to_28May19_05:16_p1.fits.fz" --plot proj 0 --ohdu 2 --E-span 500 --side 1 --fit "ax.plot( fitted_curve(lambda x, *p: p[0] + p[1]*(np.exp(-p[2]*x) - p[3]*np.exp(-p[4]*x)), [4300,-1000,1, 100, .1, 1], x, means), '.', label = 'fit')" --no-max --x-range 10 3000 --global-bias
