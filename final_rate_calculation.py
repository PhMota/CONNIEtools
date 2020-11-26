from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *
from fstring import F

folder = '../final_rate_calculation'

if not os.path.exists(folder):
    os.mkdir(folder)

ohdu_selection = '(ohdu<6 || ohdu==8 || ohdu==9 || ohdu==13 || ohdu==14)'
ohdu_selection20 = '(ohdu<5 || ohdu==9 || ohdu==13 || ohdu==14)'
geometric_selection = '(xMin>140 and xMax<3960 and yMin>75 and yMax<898 and hPixFlag==0)'
runID_selection = '(runID>6226 and runID<6975) || (runID>14712 and runID<15178)'

runID_excluded = '(runID != 6415 and runID != 6475 and runID != 6499 and runID != 6926 and runID != 6927 and runID != 6031 and runID != 6093 and runID != 6096 and runID != 6139 and runID != 7074 and runID != 7222 and runID != 7226 and runID != 7374)'
runID_excluded20 = '(runID != 14729 and runID != 14816 and runID != 13125 and runID != 13126 and runID != 13315 and runID != 13628 and runID != 13911 and runID != 14069 and runID != 14212 and runID != 14661)'

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

global_selection = r' and '.join([ohdu_selection, geometric_selection, energy_selection, size_selection])

global_selection_hpix = r' and '.join([ohdu_selection, geometric_selection, size_selection])

global_selection_hpix20 = r' and '.join([ohdu_selection20, geometric_selection, energy_selection])

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


binsizes = [.130, .02]
Emaxs = [2.05, .55]
loop = zip(binsizes, Emaxs)

histogram_command = './catalog histogram'

with Timer('presentation'):
    with openBeamer( folder, r'data analysis run2019' ) as doc:

        doc.par('\n\part{number of measured events}\n\\frame{\partpage}\n')

        for ohdu in [3]:
            # for binsize, Emax in loop:
                doc.frame( F('number of events'),
                    'the data is given as the number of detected events per recontructed energy bin $i$ per ohdu $a$ and its Poisson error',
                    Math(( N[i,a], delta*N[i,a] == sqrt(N[i,a]) )),
                    doc.column(
                    *[ makefigure(
                        histogram_command
                        # + F(' -s runID "{{nuCatalogs}}" "1" "runID"')
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "1" "OFF19"')
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFFSEL19"')
                        + F(' --global-selection "{{off19}} and ohdu=={{ohdu}} and {{runID_excluded}}"')
                        + F(' -f "hist(OFF19, bins)" "all recontructed"')
                        + F(' -f "hist(OFFSEL19, bins)" "selections"')
                        + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "counts" --legend-title "OFF19 OHDU{{ohdu}}" ')
                        , F(r'numberOfEventsOHDU{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True
                    ) for binsize, Emax in loop ],
                    widths = [.5,.5]
                    )
                )

        for ohdu in [3]:
            # for binsize, Emax in loop:
                doc.frame('number of neutrinos',
                    'the number of detected neutrinos (where $n$ is the number of images)',
                    Math(
                        delim(N[i,a])[nu] ==
                            delim(N[i,a])[nu+bg]
                            - frac( n[nu+bg], n[bg] ) * delim(N[i,a])[bg]
                    ),
                    Math(
                        delim(frac(delta * delim(N[i,a])[nu], delim(N[i,a])[nu]))**2  ==
                            frac(1, delim(N[i,a])[nu+bg])
                            + delim( frac( n[nu+bg], n[bg] ) )**2 * frac(1, delim(N[i,a])[bg]),
                    ),
                    doc.column(
                        *[makefigure(
                            histogram_command
                            # + F(' -s runID "{{nuCatalogs}}" "{{off19}}" "runIDOFF"')
                            # + F(' -s runID "{{nuCatalogs}}" "{{on19}}" "runIDON"')
                            + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                            + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                            + F(' --global-selection "ohdu=={{ohdu}} and {{runID_excluded}}"')
                            + F(' -f "hist(ON19, bins) - hist(OFF19, bins)*{{non19}}/{{noff19}}" "reconstructed neutrinos"')
                            + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "number of neutrinos" --legend-title "OHDU{{ohdu}}" ')
                            , F(r'numberOfNeutrinosOHDU{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True
                            ) for binsize, Emax in loop],
                            widths=[.5,.5]
                        ),
                )

        for ohdu in [3]:
                doc.frame('neutrino to background ratio',
                    'the neutrino to background ratio',
                    Math(
                        frac( delim(N[i,a])[nu], delim(N[i,a])[bg] ) + 1 ==
                            frac( n[bg], n[nu+bg] ) * frac(delim(N[i,a])[nu+bg], delim(N[i,a])[bg])
                    ),
                    doc.column(
                        *[ makefigure(
                            histogram_command
                            + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                            + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                            + F(' --global-selection "ohdu=={{ohdu}} and {{runID_excluded}}"')
                            + F(r' -f "hist(ON19, bins)/hist(OFF19, bins) * {{noff19}}/{{non19}}" "nu/bg+1"')
                            + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "neutrinos to background ratio + 1" --legend-title "OHDU{{ohdu}}" ')
                            , F(r'ratioOfNeutrinosOHDU{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True ) for binsize, Emax in loop],
                            widths = [.5,.5]
                        )
                )
        doc.frame('summary',
            doc.itemize(
                'the number of neutrinos measured can be caclulated by subtracting the rescaled OFF spectrum from the ON spectrum',
                'the ratio of the numbers between On and OFF gives us the neutrino to background ratio',
                )
        )

        doc.par('\n\part{measured event rate}\n\\frame{\partpage}\n')

        doc.frame('measured event rate',
            'the measured event rate is',
            Math( calR[i,a] == (
                frac( d(N['',a]), d(M)*d(T)*d(E) )(E[i]) == frac( N[i,a], m[0]*Delta*T*Delta*E )
            )),
            'where',
            doc.itemize(
                F(r'the total mass of the CCD $m_0 = 2.329{{density_units}}\times 125{{micrometer}}^2 \times 675{{micrometer}} \times 4120^2$,'),
                F(r'the total exposition time $\Delta T = n \Delta t = n \times 3600{{seconds}}$ which is the number of images, $n$, times the readout time , $\Delta t$,'),
                F(r'the energy bin $\Delta E = 130{{eVolts}}$.'),
            ),
        )
        factor = 1./(2.329/1e3*.0015**2*.0675*4120**2/24. )
        print('factor', factor, 2.329/1e3*.0015**2*.0675*4120**2)
        for ohdu in [3]:
                doc.frame('measured event rate',
                    'the measured event rate',
                    Math( calR[i,a] == (
                        frac( d(N['',a]), d(M)*d(T)*d(E) )(E[i]) == frac( N[i,a], m[0]*Delta*T*Delta*E )
                    )),
                    doc.column(
                    *[ makefigure(
                        histogram_command
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                        + F(' --global-selection "ohdu=={{ohdu}} and {{runID_excluded}}"')
                        + F(r' -f "hist(OFF19, bins)*{{factor}}/{{binsize}}/{{noff19}}" "reactor OFF"')
                        + F(r' -f "hist(ON19, bins)*{{factor}}/{{binsize}}/{{non19}}" "reactor ON"')
                        + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "OHDU{{ohdu}}" ')
                        , F(r'ratioOfNeutrinosOHDU{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True ) for binsize, Emax in loop ],
                        widths=[.5,.5]
                    )
                )

        doc.frame('differential rate',
            'the differential event rate is',
            Math(
                delim(calR[nu])[i,a] == (
                        delim(calR[nu+bg])[i,a] - delim(calR[bg])[i,a],
                        # Delta*frac( d(N['',a]), d(M)*d(T)*d(E) )(E[i]),
                        delim( frac( d(N['',a]), d(M)*d(T)*d(E) )(E[i]) )[nu+bg]
                        - delim( frac( d(N['',a]), d(M)*d(T)*d(E) )(E[i]) )[bg],
                        frac(1, m[0]*Delta*E*Delta*t)*delim(
                            delim( frac( N[i,a], n ) )[nu+bg]
                            - delim( frac( N[i,a], n ) )[bg]
                        ),
                        frac(1, m[0]*Delta*E*n[nu+bg]*Delta*t)*delim(
                            N_nubg
                            - n[nu+bg] * delim( frac( N[i,a], n ) )[bg]
                        ),
                        frac( N_nu, m[0]*Delta*E*n[nu+bg]*Delta*t ),
                )
            ),
        )

        for ohdu in [3]:
                doc.frame('differential rate',
                    'the differential event rate is',
                    Math(
                        delim(calR[nu])[i,a] == (
                                delim(calR[nu+bg])[i,a] - delim(calR[bg])[i,a],
                        )
                    ),
                    doc.column(
                    *[ makefigure(
                        histogram_command
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                        + F(' --global-selection "ohdu=={{ohdu}} and {{runID_excluded}}"')
                        + F(r' -f "hist(ON19, bins)*{{factor}}/{{binsize}}/{{non19}} - hist(OFF19, bins)*{{factor}}/{{binsize}}/{{noff19}}" "reactor ON--OFF(1)"')
                        + F(r' -f "( hist(ON19, bins) - hist(OFF19, bins)*{{non19}}/{{noff19}} )*{{factor}}/{{binsize}}/{{non19}}" "reactor ON--OFF(2)"')
                        + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "OHDU{{ohdu}}" ')
                        , F(r'ratioOfNeutrinosOHDU{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True
                    ) for binsize, Emax in loop ],
                    widths=[.5,.5]
                    )

                )

        for ohdu in [3]:
                doc.frame('ratio of rates',
                    'the differential event rate is the neutrino to background ratio',
                    Math(
                        frac( delim(calR[nu+bg])[i,a], delim(calR[bg])[i,a] ) == (
                            frac( n[bg], n[nu+bg] ) * frac(delim(N[i,a])[nu+bg], delim(N[i,a])[bg]),
                            frac( delim(N[i,a])[nu], delim(N[i,a])[bg] ) + 1
                        )
                    ),
                    doc.column(
                    *[ makefigure(
                        histogram_command
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                        + F(' --global-selection "ohdu=={{ohdu}} and {{runID_excluded}}"')
                        + F(r' -f "( hist(ON19, bins)*{{factor}}/{{binsize}}/{{non19}} )/( hist(OFF19, bins)*{{factor}}/{{binsize}}/{{noff19}} )" "reactor ON/OFF(1)"')
                        + F(r' -f "( hist(ON19, bins)/{{non19}} )/( hist(OFF19, bins)/{{noff19}} )" "reactor ON/OFF(2)"')
                        + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "OHDU{{ohdu}}" ')
                        , F(r'ratioOfRatesOHDU{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True
                    ) for binsize, Emax in loop ],
                    widths=[.5,.5]
                    )

                )

        doc.frame('summary',
            doc.itemize(
                'the measured -- not corrected for efficiency -- rate was computed',
                'the differential event rate is the neutrino rate',
                'the rate ratio is the same as the measured event ratio, no need to correct for exposure',
                )
        )

        doc.par('\n\part{estimated(!) event rate}\n\\frame{\partpage}\n')

        doc.frame('number of events',
            'the number of estimated events is the number of measured events corrected for the efficiency',
            Math((
                N[i,'*'*a] == frac( N_U, N_S )*N[i,a],
                # delim(frac( (delta*N[i,a]), N[i,a] ))**2 ==
                #     frac( 1, N[i,'*'*a] )
                #     + frac( 1, N_S )
                #     + frac( 1, N_U )
            )),
            'where',
            doc.itemize(
                F(r'${{ N_U }}$ is the {\bf total} number of simulated events in all of the CCD'),
                F(r'${{ N_S }}$ is the number of reconstructed events in the {\bf selection}'),
            ),
            'the efficiency is needed to compare the measured rates with the theoretical models'
        )

        for ohdu in [3]:
            doc.frame( F('reconstructed from simulation'),
                Math((
                    N_U, N_S
                )),
                doc.column(
                *[makefigure(
                    histogram_command
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "1" "SIM"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,2.05)}}" "REC"')
                    + F(' --global-selection "{{off19}} and ohdu=={{ohdu}} and {{runID_excluded}}"')
                    + F(' -f "hist(SIM, bins)" "all recontructed"')
                    + F(' -f "hist(REC, bins)" "reconstructed with selections"')
                    + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "counts" --legend-title "OHDU{{ohdu}}"')
                    , F(r'spectrum_match_19_ohdu{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True
                ) for binsize, Emax in loop ]
                , widths=[.5,.5]
                ),
            )

        for ohdu in [3]:
            doc.frame(F('efficiency as the ratio'),
                Math((
                    N[i,'*'*a] == ( frac( N_U, N_S )*N[i,a], frac( N[i,a], Expr(r'\epsilon' ) ) )
                )),
                doc.column(
                *[ makefigure(
                    histogram_command
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{off19}}" "SIMOFF"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,2.05)}}" "RECOFF"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{on19}}" "SIMON"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,2.05)}}" "RECON"')
                    + F(' --global-selection " ohdu=={{ohdu}} and {{runID_excluded}}"')
                    + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "efficiency" --legend-title "OHDU{{ohdu}}"')
                    + F(' -f "hist(RECOFF, bins)/hist(SIMOFF, bins)" "reactor OFF"')
                    + F(' -f "hist(RECON, bins)/hist(SIMON, bins)" "reactor ON"')
                    , F(r'efficiency_19{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True )
                    for binsize, Emax in loop ]
                    , widths =[.5,.5]
                    )
            )

        for ohdu in [3]:
            doc.frame(F('estimated neutrino rate'),
                'the estimated event rate',
                Math( calR[i,'*'*a] == (
                    frac( N[i,'*'*a], m[0]*Delta*T*Delta*E )
                )),
                doc.column(
                *[ makefigure(
                    histogram_command
                    + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                    + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{off19}}" "SIMOFF"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,2.05)}}" "RECOFF"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{on19}}" "SIMON"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,2.05)}}" "RECON"')
                    + F(' --global-selection " ohdu=={{ohdu}} and {{runID_excluded}}"')
                    + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "efficiency" --legend-title "OHDU{{ohdu}}"')
                    + F(' -f "hist(ON19, bins)*hist(SIMON, bins)/hist(RECON, bins)*{{factor}}/{{binsize}}/{{non19}} - hist(OFF19, bins)*hist(SIMOFF, bins)/hist(RECOFF, bins)*{{factor}}/{{binsize}}/{{noff19}}" "reactor OFF corrected"')
                    + F(' -f "hist(ON19, bins)*{{factor}}/{{binsize}}/{{non19}} - hist(OFF19, bins)*{{factor}}/{{binsize}}/{{noff19}}" "reactor OFF uncorrected"')
                    , F(r'efficiency_19{{ohdu}}.pdf'), doc, height = .5, folder=folder, nocode=True )
                    for binsize, Emax in loop ]
                    , widths =[.5,.5]
                    )
            )

        # doc.par('\n\part{backup}\n\\frame{\partpage}\n')

        doc.frame(F('efficiency denominator revisited'),
            Math((
                N_U, N_S,
            )),
            doc.figure( 'efficiency_drawn.png', height=.6 ),
        )

        doc.frame(F('efficiency as the normalized numerator'),
            # Math((
            #     Expr(r'\varpesilon ') == frac( N_S(1), N_S ),
            # )),
            'the efficiency should be computed with the flat simulated distribution and can be adjusted by a smooth function ',
            doc.figure( 'carla_efficiency_fit.png', height=.6 ),
        )

        doc.frame(F('comparing to theoretical models'),
            'instead of correcting our rate with the binned efficiency, the models are corrected to our experimental conditions using the adjusted efficiency -- as it was done in the previous papers',
            doc.figure( 'paper_model.png', height=.6 ),
        )

        doc.frame(F('simulated and reconstructed energy'),
            'since the models use the simulated energy and our rates are computed with the reconstructed energy, the models must be convolved with our energy uncertainty',
            doc.figure( 'sim_rec_E.png', height=.6 ),
        )

        doc.frame(F('summary'),
            doc.itemize(
                'there was a LOT of confusion during the last month or two about the efficiency meaning and calculation',
                'in the talk I summarize part of this discussion',
                'I compare the differential rate with the neutrino to backgound ratio',
                'there are two ways to compare our rates with the theoretical models'
                + doc.itemize(
                    'either correct the measured rate and compare with the model',
                    'or correct the models and compare to the measured rate',
                ),
                'the first option has the caveat of trying to correct the first bin with a non-constant efficiency',
                'on the second option, one integrates the models after properly correcting for the efficiency',
                'the models use the simulated energy which must be converted to the reconstructed energy'
            )
        )

        doc.frame(F('conclusion'),
            'my proposal',
            doc.itemize(
                'we should go back to not correcting the rates, but compare them with the correcte models',
                'however, the efficiency cannot be computed as the ratio of reconstructed spectra, it needs to use the flat distribution that is fed into the simulation algorithm',
                'by adjusting the efficiency, we can actually normalize the saturated region as 100\%, since we are correcting and not the data',
                'on this note, Guille\'s plots of uncorrected rates are fine! but our efficiency was wrong in the papers'
            )
        )

        doc.frame(F('correcting the model'),
            'the expected neutrino rate for nuclear recoil energy at CONNIE is',
            doc.figure( 'expr_neutrino_rate.png', height=.1 ),
            'converting to ionization energy using the quenching factor',
            doc.figure( 'expr_ion.png', height=.1 ),
            'correcting for the reconstructed (measured) energy',
            doc.figure( 'expr_sim_rec.png', height=.3 ),
        )


        for ohdu in ohdus19:
                doc.frame(F('differential rate ohdu {{ohdu}}'),
                    'the differential event rate is',
                    Math(
                        delim(calR[nu])[i,a] == (
                                delim(calR[nu+bg])[i,a] - delim(calR[bg])[i,a],
                        )
                    ),
                    doc.column(
                    *[ makefigure(
                        histogram_command
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19"')
                        + F(' -s E1/gain3Peaks "{{nuCatalogs}}" "{{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19"')
                        + F(' --global-selection "ohdu=={{ohdu}} and {{runID_excluded}}"')
                        + F(r' -f "hist(ON19, bins)*{{factor}}/{{binsize}}/{{non19}} - hist(OFF19, bins)*{{factor}}/{{binsize}}/{{noff19}}" "reactor ON--OFF"')
                        + F(' --hide-zeros  --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "OHDU{{ohdu}}" ')
                        , F(r'differentialRateOHDU{{ohdu}}binsize{{binsize}}.pdf'), doc, height = .5, folder=folder, nocode=True
                    ) for binsize, Emax in loop ],
                    widths=[.5,.5]
                    )

                )

        w = Expr('w')
        expr = ', '.join( [ F('hist(ON19_{{ohdu}}, bins)*{{factor}}/{{non19}} - hist(OFF19_{{ohdu}}, bins)*{{factor}}/{{noff19}}').str() for ohdu in ohdus19 ] )

        print( 'expr', expr )
        load_off = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{off19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "OFF19_{{ohdu}}"').str() for ohdu in ohdus19 ])
        print( 'expr', load_off )
        load_on = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{on19}} and {{geometric_selection}} and {{size_selection}} and {{energy_selection_function(0.05,2.05)}}" "ON19_{{ohdu}}"').str() for ohdu in ohdus19 ])

        doc.frame(F('weighted average of the differential rate'),
            'the weighted average of the differential event rate is',
            Math((
                ave(delim(calR[nu])[i,a]) == (
                    frac( Sum[a]* w[i,a] * delim(calR[nu])[i,a], Sum[a]* w[i,a] )
                ),
                w[i,a] == frac(1, delim(delta*delim(calR[nu])[i,a])**2)
                )
            ),
            doc.column(
            *[ makefigure(
                histogram_command
                + load_off
                + load_on
                + F(' --global-selection "{{runID_excluded}}"')
                + F(r' -f "waverage( [{{expr}}] )/{{binsize}}" "reactor ON--OFF"')
                + F(' --hide-zeros  --no-title  --no-label-file  --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "weighted average" ')
                , F(r'differentialRateOHDUwAverageBinsize{{binsize}}.pdf'), doc, height = .5, folder=folder, nocode=True
            ) for binsize, Emax in loop ]
            , widths=[.5,.5]
            )
        )

        unbias_load_off = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{off19}} and {{unbias_energy_selection_function(0.05,2.05)}}" "OFF19u_{{ohdu}}"').str() for ohdu in ohdus19 ])
        print( 'expr', load_off )
        unbias_load_on = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{on19}} and {{unbias_energy_selection_function(0.05,2.05)}}" "ON19u_{{ohdu}}"').str() for ohdu in ohdus19 ])
        unbias_expr = ', '.join( [ F('hist(ON19u_{{ohdu}}, bins)*{{factor}}/{{non19}} - hist(OFF19u_{{ohdu}}, bins)*{{factor}}/{{noff19}}').str() for ohdu in ohdus19 ] )

        unbiasE_load_off = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{off19}} and {{unbias_energy_selection_function(0.05,2.05)}} and E0/gain3Peaks>.045" "OFF19uE_{{ohdu}}"').str() for ohdu in ohdus19 ])
        print( 'expr', load_off )
        unbiasE_load_on = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{on19}} and {{unbias_energy_selection_function(0.05,2.05)}} and E0/gain3Peaks>.045" "ON19uE_{{ohdu}}"').str() for ohdu in ohdus19 ])
        unbiasE_expr = ', '.join( [ F('hist(ON19uE_{{ohdu}}, bins)*{{factor}}/{{non19}} - hist(OFF19uE_{{ohdu}}, bins)*{{factor}}/{{noff19}}').str() for ohdu in ohdus19 ] )

        doc.frame(F('weighted average of the differential rate'),
            'the weighted average of the differential event rate is',
            Math((
                ave(delim(calR[nu])[i,a]) == (
                    frac( Sum[a]* w[i,a] * delim(calR[nu])[i,a], Sum[a]* w[i,a] )
                ),
                w[i,a] == frac(1, delim(delta*delim(calR[nu])[i,a])**2)
                )
            ),
            doc.column(
            *[ makefigure(
                histogram_command
                + load_off
                + load_on
                + F(r' -f "waverage( [{{expr}}] )/{{binsize}}" "reactor ON--OFF"')
                + unbias_load_off
                + unbias_load_on
                + F(r' -f "waverage( [{{unbias_expr}}] )/{{binsize}}" "unbiased reactor ON--OFF"')
                + unbiasE_load_off
                + unbiasE_load_on
                + F(r' -f "waverage( [{{unbiasE_expr}}] )/{{binsize}}" "energy selection reactor ON--OFF"')
                + F(' --global-selection "{{runID_excluded}}"')
                + F(' --hide-zeros  --no-title  --no-label-file  --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "weighted average" ')
                , F(r'unbiasedDifferentialRateOHDUwAverageBinsize{{binsize}}.pdf'), doc, height = .5, folder=folder, nocode=True
            ) for binsize, Emax in loop ]
            , widths=[.5,.5]
            )
        )

        unbiasES_load_off = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{off19}} and {{unbias_energy_selection_function(0.05,2.05)}} and E0/gain3Peaks>.045 and {{size_selection}}" "OFF19uES_{{ohdu}}"').str() for ohdu in ohdus19 ])
        print( 'expr', load_off )
        unbiasES_load_on = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{on19}} and {{unbias_energy_selection_function(0.05,2.05)}} and E0/gain3Peaks>.045 and {{size_selection}}" "ON19uES_{{ohdu}}"').str() for ohdu in ohdus19 ])
        unbiasES_expr = ', '.join( [ F('hist(ON19uES_{{ohdu}}, bins)*{{factor}}/{{non19}} - hist(OFF19uES_{{ohdu}}, bins)*{{factor}}/{{noff19}}').str() for ohdu in ohdus19 ] )

        doc.frame(F('weighted average of the differential rate'),
            'the weighted average of the differential event rate is',
            Math((
                ave(delim(calR[nu])[i,a]) == (
                    frac( Sum[a]* w[i,a] * delim(calR[nu])[i,a], Sum[a]* w[i,a] )
                ),
                w[i,a] == frac(1, delim(delta*delim(calR[nu])[i,a])**2)
                )
            ),
            doc.column(
            *[ makefigure(
                histogram_command
                + load_off
                + load_on
                + F(r' -f "waverage( [{{expr}}] )/{{binsize}}" "reactor ON--OFF"')
                + unbiasE_load_off
                + unbiasE_load_on
                + F(r' -f "waverage( [{{unbiasE_expr}}] )/{{binsize}}" "energy selection reactor ON--OFF"')
                + unbiasES_load_off
                + unbiasES_load_on
                + F(r' -f "waverage( [{{unbiasES_expr}}] )/{{binsize}}" "energy \& size selections reactor ON--OFF"')
                + F(' --global-selection "{{runID_excluded}}"')
                + F(' --hide-zeros  --no-title  --no-label-file  --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "weighted average" ')
                , F(r'unbiasedDifferentialRateOHDUwAverageBinsize{{binsize}}.pdf'), doc, height = .5, folder=folder, nocode=True
            ) for binsize, Emax in loop ]
            , widths=[.5,.5]
            )
        )

        unbiasE2_load_off = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{off19}} and {{unbias_energy_selection_function(0.05,2.05)}} and E0/gain3Peaks>.05" "OFF19uE2_{{ohdu}}"').str() for ohdu in ohdus19 ])
        print( 'expr', load_off )
        unbiasE2_load_on = ''.join( [F(' -s E1/gain3Peaks "{{nuCatalogs}}" "ohdu=={{ohdu}} and {{on19}} and {{unbias_energy_selection_function(0.05,2.05)}} and E0/gain3Peaks>.05" "ON19uE2_{{ohdu}}"').str() for ohdu in ohdus19 ])
        unbiasE2_expr = ', '.join( [ F('hist(ON19uE2_{{ohdu}}, bins)*{{factor}}/{{non19}} - hist(OFF19uE2_{{ohdu}}, bins)*{{factor}}/{{noff19}}').str() for ohdu in ohdus19 ] )

        doc.frame(F('weighted average of the differential rate'),
            'the weighted average of the differential event rate is',
            Math((
                ave(delim(calR[nu])[i,a]) == (
                    frac( Sum[a]* w[i,a] * delim(calR[nu])[i,a], Sum[a]* w[i,a] )
                ),
                w[i,a] == frac(1, delim(delta*delim(calR[nu])[i,a])**2)
                )
            ),
            doc.column(
            *[ makefigure(
                histogram_command
                + load_off
                + load_on
                + F(r' -f "waverage( [{{expr}}] )/{{binsize}}" "reactor ON--OFF"')
                + unbiasE_load_off
                + unbiasE_load_on
                + F(r' -f "waverage( [{{unbiasE_expr}}] )/{{binsize}}" "energy selection(1) reactor ON--OFF"')
                + unbiasE2_load_off
                + unbiasE2_load_on
                + F(r' -f "waverage( [{{unbiasE2_expr}}] )/{{binsize}}" "energy selection(2) reactor ON--OFF"')
                + F(' --global-selection "{{runID_excluded}}"')
                + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 {{Emax}} --binsize {{binsize}} --xlabel "energy [keV]" --ylabel "rate [events/kg/day/keV]" --legend-title "weighted average" ')
                , F(r'unbiasedDifferentialRateOHDUwAverageBinsize{{binsize}}.pdf'), doc, height = .5, folder=folder, nocode=True
            ) for binsize, Emax in loop ]
            , widths=[.5,.5]
            )
        )
