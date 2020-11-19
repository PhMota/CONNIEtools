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

with Timer('presentation'):
    with openBeamer( folder, r'final rate calculation\\ https://github.com/PhMota/CONNIEtools/final_rate_calculation.py' ) as doc:

        doc.frame('number of events',
            'the data is given as the number of events per recontructed energy bin $i$ per ohdu $a$ and its Poisson error',
            Math(( N[i,'*'*a], delta*N[i,'*'*a] == sqrt(N[i,'*'*a]) )),
            'the number of events corrected for efficiency and its propagated error are',
            Math((
                N[i,a] == frac( N_S, N_U )*N[i,'*'*a],
                delim(frac( (delta*N[i,a]), N[i,a] ))**2 ==
                    frac( 1, N[i,'*'*a] )
                    + frac( 1, N_S )
                    + frac( 1, N_U )
            )),
            'where',
            doc.itemize(
                F(r'${{ N_S }}$ is the total number of simulated events in all of the CCD for each reconstructed energy bin $i$ and ohdu $a$'),
                F(r'${{ N_U }}$ is the total number of reconstructed events in the active selection of the CCD for each reconstructed energy bin $i$ and ohdu $a$'),
            ),
        )

        doc.frame('number of neutrinos',
            'the number of neutrinos corrected for efficiency',
            Math(
                delim(N[i,a])[nu] ==
                    delim(N[i,a])[nu+bg]
                    - delim( frac( n[nu+bg], n[bg] ) ) * delim(N[i,a])[bg]
            ),
            Math(
                delta * delim(N[i,a])[nu] == sqrt(
                    delim(N[i,a])[nu+bg]
                    - delim( frac( n[nu+bg], n[bg] ) ) * delim(N[i,a])[bg],
                )
            ),
            F("where $n$ is the number of images and it was assumed ${{N_U}}\gg {{N[i,'*'*a]}}$"),
            'the ratio is',
            Math(
                delim(N[i,a])[nu] ==
                    delim(N[i,a])[nu+bg]
                    - delim( frac( n[nu+bg], n[bg] ) ) * delim(N[i,a])[bg]
            ),

        )

        doc.frame('rate from number of events',
            'the number of events corrected for efficiency and its propagated error',
            Math((
                N[i,a],
                delta*N[i,a] == frac( N[i,a], sqrt(N[i,'*'*a]) )
                    + frac( N[i,a], sqrt(N_S) )
            )),
            'the rate of events correctd for efficiency',
            Math( calR[i,a] == (
                frac( d(N['',a]), d(M)*d(T)*d(E) )(E[i]) == frac( N[i,a], m[0]*Delta*T*Delta*E )
            )),
            'where',
            doc.itemize(
                F(r'the total mass of the CCD $m_0 = 2.42{{density_units}}\times 125{{micrometer}}^2 \times 675{{micrometer}}$,'),
                F(r'the total exposition time $\Delta T = n \Delta t = n \times 3600{{seconds}}$ which is the number of images, $n$, times the readout time , $\Delta t$,'),
                F(r'the energy bin $\Delta E = 130{{eVolts}}$.'),
            ),
        )

        doc.frame('differential rate',
            'the differential rate of events is',
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

        for ohdu in [2]:
            doc.frame( F('efficiency off19 0.05--1.05keV ohdu {{ohdu}}'),
                makefigure(
                    F('./catalog histogram')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "1" "SIM"')
                    + F(' -s E1/gain3Peaks "{{simCatalogs}}" "{{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,1.05)}}" "REC"')
                    + F(' --global-selection "{{off19}} and ohdu=={{ohdu}} and {{runID_excluded}}"')
                    + F(' -f "hist(SIM)" "recontructed"')
                    + F(' -f "hist(REC)" "reconstructed with selections"')
                    + F(' --hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 1.05 --binsize .01 --xlabel "energy [keV]" --ylabel "counts" --average 0.3 1')
                    , F(r'spectrum_match_19_ohdu{{ohdu}}.pdf'), doc, height = .6, folder=folder ),
            )

            doc.frame(F('efficiency on19 and off19 0.05--1.05keV ohdu {{ohdu}}'),
                makefigure(
                    F('./catalog histogram')
                    + F('-s E1/gain3Peaks "{{simCatalogs}}" "1" "SIM"')
                    + F('-s E1/gain3Peaks "{{simCatalogs}}" "{{geometric_selection}} and {{size_selection}} and distSim<1.5 and {{energy_selection_function(0.05,1.05)}}" "REC"')
                    + F('--global-selection "{{off19}} and ohdu=={{ohdu}} and {{runID_excluded}}"')
                    + F('--hide-zeros --no-title --no-label-file --no-label-branch --x-range .05 1.05 --binsize .01')
                    + F('-f "divide(hist(REC), hist(SIM)" "selections/all"')'
                    , F(r'efficiency_19{{ohdu}}.pdf'), doc, height = .6, folder=folder ),
            )
