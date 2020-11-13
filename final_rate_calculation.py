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
    with openBeamer( folder, r'final rate calculation\\ https://github.com/PhMota/CONNIEtools/final_rate_calculation.py' ) as doc:

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
