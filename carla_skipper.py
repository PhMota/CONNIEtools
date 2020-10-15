from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *

folder = 'carla_skipper'

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
    with openBeamer( folder, 'skipper study with neutrino rate\\\\https://github.com/PhMota/CONNIEtools/wiki' ) as doc:
        
        doc.par('\n\part{''}\n\\frame{\partpage}\n')

        doc.frame('study',
            'methodology',
            doc.itemize(
                r'generate $\nu$ -- 4k$\times$4k pixels, 675$\mu$m',
                r'2 e$e^-$ noise and dark current of 0.05/h',
                r'using the neutrino rate in counts/day/kg/keV',
            )
        )
        
        geom = '--ccd-shape 4000 4000 --horizontal-overscan 150 --vertical-overscan 50 --charge-gain 7.5'
        
        doc.frame('simulate neutrino rate populations',
            'neutrino images, 20k each',
            *[ requiredFile('./simulation image {0}/nu{2} -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 {1}'.format(folder, geom, n),
                         'nu{0}.fits'.format(n), doc ) for n in range(5)]
        )
        
        doc.frame('simulate noise and dc images',
            '0.3 dc and 2 e$^-$ noise',
            *[ requiredFile('./simulation image {0}/noisedc{2} -N 0 --dark-current 0.15 --readout-noise 15 {1}'.format(folder, geom, n),
                         'noisedc{0}.fits'.format(n), doc ) for n in range(5) ]
        )

        doc.frame('simulate noise and dc images',
            'combine images',
            *[ requiredFile('./image superpose {0}/nu{1}.fits {0}/noisedc{1}.fits --output {0}/sim{1}{1}.fits'.format(folder, n),
                         'sim{0}{0}.fits'.format(n), doc ) for n in range(5)]
        )

        doc.frame('extract',
            'extract threshold 6e$^-$ = 45ADU and 8e$^-$ = 60ADU',
            *[ requiredFile('./catalog extract \
                {0}/ext{1}{1}_{2} {0}/sim{1}{1}.fits \
                -t {2} -b 1 --ohdu -1 -v\
                '.format(folder, n, m),
                'ext{0}{0}_{1}.root'.format(n, m), doc ) for n in range(5) for m in [45, 60] ]
        )
        
        doc.frame('matching',
            *[ requiredFile('./catalog match \
                {0}/nu{1} \
                {0}/ext{1}{1}_{2}.root \
                -o {0}/ext{1}{1}_{2}.match.root\
                '.format(folder, n, m),
                'ext{0}{0}_{1}.match.root'.format(n, m), doc ) for n in range(5) for m in [45, 60] ]
        )

        doc.frame('spectra',
            doc.column( *[ makefigure('./catalog histogram\
                "E0" {0}/ext00_{1}.match.root 1\
                "E1" {0}/ext00_{1}.match.root 1\
                --log --png --binsize 20 \
                --output {0}/ext00_{1}.E20b\
                '.format(folder, n),
                'ext00_{0}.E20b.png'.format(n), doc, height=.4 ) for n in [45,60] ],
                widths=[.4,.4]
            )
        )
        
        doc.frame('spectra 0--100eV',
            doc.column( *[ makefigure('./catalog histogram\
                "E0" {0}/ext00_{1}.match.root 1\
                "E1" {0}/ext00_{1}.match.root 1\
                --log --png --binsize 10 \
                --output {0}/ext00_{1}.E10 --x-range 0 200\
                '.format(folder, n),
                'ext00_{0}.E10.png'.format(n), doc, height=.4 ) for n in [45,60] ],
                widths=[.4,.4]
            )
        )
        
        doc.frame('matched spectrum',
            doc.column( *[ makefigure('./catalog histogram \
                "E0" {0}/ext{1}{1}_{2}.match.root 1\
                "E0" {0}/ext{1}{1}_{2}.match.root "idSim==0"\
                "E0" {0}/ext{1}{1}_{2}.match.root "idSim>0"\
                --log --png --binsize 10 --output {0}/ext{1}{1}_{2}.match.E10d{2} --x-range 0 200'.format(folder, 0, n),
                'ext{0}{0}_{1}.match.E10d{1}.png'.format(0,n), doc, height=.4 ) for n in [45, 60] ],
                widths=[.4,.4]
            )
        )

        doc.frame('ratio',
            doc.column( *[ makefigure('./catalog ratio \
                "E0" {0}/ext{1}{1}_{2}.match.root "idSim==0" {0}/ext{1}{1}_{2}.match.root 1\
                "E0" {0}/ext{1}{1}_{2}.match.root "idSim>0" {0}/ext{1}{1}_{2}.match.root 1\
                --png --binsize 10 --output {0}/ext{1}{1}_{2}.ratio.match.E10{2} --x-range 0 200\
                '.format(folder, 0, n),
                'ext{0}{0}_{1}.ratio.match.E10{1}.png'.format(0, n), doc, height=.4 ) for n in [45, 60] ],
                widths=[.4,.4]
            )
        )

        doc.frame('differential spectrum',
            doc.column( *[ makefigure('./catalog histogram \
                "(E1-E0)/(n1-n0)" {0}/ext{1}{1}_{2}.match.root 1\
                "(E1-E0)/(n1-n0)" {0}/ext{1}{1}_{2}.match.root "idSim==0"\
                "(E1-E0)/(n1-n0)" {0}/ext{1}{1}_{2}.match.root "idSim>0"\
                --log --png --binsize 5 --output {0}/ext{1}{1}_{2}.match.EEnn5_thr{2}\
                --x-range -60 60\
                '.format(folder, 0, n),
                'ext{0}{0}_{1}.match.EEnn5_thr{1}.png'.format(0, n), doc, height=.4 ) for n in [45, 60] ],
                widths=[.4,.4]
            )
        )

        #doc.frame('ratio 60ADU, E0$\gt$50eV',
            #doc.column( *[ makefigure('./catalog ratio \
                #E{2} {0}/ext{1}{1}_60.match.root "idSim==0 and E0>80" {0}/ext{1}{1}_60.match.root "idSim==0"\
                #E{2} {0}/ext{1}{1}_60.match.root "idSim>0 and E0>80" {0}/ext{1}{1}_60.match.root "idSim>0"\
                #--png --binsize 10 --output {0}/ext{1}{1}_60.ratio.match.E10E80a_{2} --x-range 0 200\
                #'.format(folder, 0, n),
                #'ext{0}{0}_60.ratio.match.E10E80a_{1}.png'.format(0, n), doc, height=.4 ) for n in range(2) ],
                #widths=[.4,.4]
            #)
        #)
            
            
            
        ####################### part II #######################
        doc.par('\n\part{''}\n\\frame{\partpage}\n')

        doc.frame('study',
            'methodology',
            doc.itemize(
                r'generate $\nu$ skipper 6k$\times$1k pixels, 675$\mu$m',
                r'0.2 e$e^-$ noise and dark current of 0.05/h',
                r'using the neutrino rate in counts/day/kg/keV',
            )
        )
        geom = '--ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'
        doc.frame('simulate neutrino rate populations',
            'neutrino images, 20k each',
            *[ requiredFile('./simulation image {0}/skipper_nu{2} -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 {1}'.format(folder, geom, n),
                         'skipper_nu{0}.fits'.format(n), doc ) for n in range(5) ]
        )

        doc.frame('simulate noise and dc images',
            '0.05 dc and 0.2 e$^-$ noise',
            *[ requiredFile('./simulation image {0}/skipper_noisedc{2} -N 0 --dark-current 0.05 --readout-noise 1.5 {1}'.format(folder, geom, n),
                         'skipper_noisedc{0}.fits'.format(n), doc ) for n in range(5) ]
        )

        doc.frame('simulate noise and dc images',
            'combine images',
            *[ requiredFile('./image superpose {0}/skipper_nu{1}.fits {0}/skipper_noisedc{1}.fits --output {0}/skipper_sim{1}{1}.fits'.format(folder, n),
                         'skipper_sim{0}{0}.fits'.format(n), doc ) for n in range(5) ]
        )

        doc.frame('simulate noise and dc images',
            'extract threshold 2e$^-$ = 15ADU',
            *[ requiredFile('./catalog extract {0}/skipper_ext{1}{1} {0}/skipper_sim{1}{1}.fits -t 15 -b 1 --ohdu -1 -v'.format(folder, n),
                         'skipper_ext{0}{0}.root'.format(n), doc ) for n in range(5) ]
        )

        doc.frame('spectrum',
            'spectrum 15ADU',
            doc.column( *[ makefigure('./catalog histogram {0}/skipper_ext{1}{1}.root --log --branch-selections E0 1 --branch-selections E1 1 --png --binsize 20 --output {0}/skipper_ext{1}{1}.E20 --x-range 0 400'.format(folder, n),
                         'skipper_ext{0}{0}.E20.png'.format(n), doc, height=.5 ) for n in range(2) ],
                widths=[.4,.4]
            )
        )
        
        doc.frame('spectrum',
            'spectrum 15ADU',
            doc.column( *[ makefigure('./catalog histogram {0}/skipper_ext{1}{1}.root --log --branch-selections E0 1 --branch-selections E1 1 --png --binsize 10 --output {0}/skipper_ext{1}{1}.E10 --x-range 0 200'.format(folder, n),
                         'skipper_ext{0}{0}.E10.png'.format(n), doc, height=.5 ) for n in range(2) ],
                widths=[.4,.4]
            )
        )
