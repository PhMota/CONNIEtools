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

def makefigure( code, file, doc, height=1 ):
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
            *[ requiredFile('./simulation image {0}/nu{0} -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 {1}'.format(folder,geom,n),
                         'nu{0}.fits'.format(n), doc ) for n in range(5)],
        )
        
        exit(0)
        doc.frame('simulate noise and dc images',
            '0.3 dc and 2 e$^-$ noise',
            requiredFile('./simulation image {0}/noisedc0 -N 0 --dark-current 0.15 --readout-noise 15 --ccd-shape 4000 4000 --horizontal-overscan 150 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc0.fits', doc ),
            requiredFile('./simulation image {0}/noisedc1 -N 0 --dark-current 0.15 --readout-noise 15 --ccd-shape 4000 4000 --horizontal-overscan 150 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc1.fits', doc ),
            requiredFile('./simulation image {0}/noisedc2 -N 0 --dark-current 0.15 --readout-noise 15 --ccd-shape 4000 4000 --horizontal-overscan 150 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc2.fits', doc ),
            requiredFile('./simulation image {0}/noisedc3 -N 0 --dark-current 0.15 --readout-noise 15 --ccd-shape 4000 4000 --horizontal-overscan 150 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc3.fits', doc ),
            requiredFile('./simulation image {0}/noisedc4 -N 0 --dark-current 0.15 --readout-noise 15 --ccd-shape 4000 4000 --horizontal-overscan 150 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc4.fits', doc ),
        )

        doc.frame('simulate noise and dc images',
            'combine images',
            requiredFile('./image superpose {0}/sim0.fits {0}/noisedc0.fits --output {0}/sim00.fits'.format(folder),
                         'sim00.fits', doc ),
            requiredFile('./image superpose {0}/sim1.fits {0}/noisedc1.fits --output {0}/sim11.fits'.format(folder),
                         'sim11.fits', doc ),
            requiredFile('./image superpose {0}/sim2.fits {0}/noisedc2.fits --output {0}/sim22.fits'.format(folder),
                         'sim22.fits', doc ),
            requiredFile('./image superpose {0}/sim3.fits {0}/noisedc3.fits --output {0}/sim33.fits'.format(folder),
                         'sim33.fits', doc ),
            requiredFile('./image superpose {0}/sim4.fits {0}/noisedc4.fits --output {0}/sim44.fits'.format(folder),
                         'sim44.fits', doc ),
        )

        doc.frame('simulate noise and dc images',
            'extract threshold 2e$^-$ = 15ADU',
            requiredFile('./catalog extract {0}/ext00 {0}/sim00.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext00.root', doc ),
            requiredFile('./catalog extract {0}/ext11 {0}/sim11.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext11.root', doc ),
            requiredFile('./catalog extract {0}/ext22 {0}/sim22.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext22.root', doc ),
            requiredFile('./catalog extract {0}/ext33 {0}/sim33.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext33.root', doc ),
            requiredFile('./catalog extract {0}/ext44 {0}/sim44.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext44.root', doc ),
        )

        doc.par('\n\part{''}\n\\frame{\partpage}\n')

        doc.frame('study',
            'methodology',
            doc.itemize(
                r'generate $\nu$ skipper 6k$\times$1k pixels, 675$\mu$m',
                r'0.2 e$e^-$ noise and dark current of 0.05/h',
                r'using the neutrino rate in counts/day/kg/keV',
            )
        )
        
        doc.frame('simulate neutrino rate populations',
            'neutrino images, 20k each',
            requiredFile('./simulation image {0}/skipper_nu0 -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'skipper_nu0.fits', doc ),
            requiredFile('./simulation image {0}/skipper_nu1 -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'skipper_nu1.fits', doc ),
            requiredFile('./simulation image {0}/skipper_nu2 -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'skipper_nu2.fits', doc ),
            requiredFile('./simulation image {0}/skipper_nu3 -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'skipper_nu3.fits', doc ),
            requiredFile('./simulation image {0}/skipper_nu4 -N 20000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'skipper_nu4.fits', doc ),
        )

        doc.frame('simulate noise and dc images',
            '0.05 dc and 0.2 e$^-$ noise',
            requiredFile('./simulation image {0}/skipper_noisedc0 -N 0 --dark-current 0.05 --readout-noise 1.5 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc0.fits', doc ),
            requiredFile('./simulation image {0}/skipper_noisedc1 -N 0 --dark-current 0.05 --readout-noise 1.5 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc1.fits', doc ),
            requiredFile('./simulation image {0}/skipper_noisedc2 -N 0 --dark-current 0.05 --readout-noise 1.5 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc2.fits', doc ),
            requiredFile('./simulation image {0}/skipper_noisedc3 -N 0 --dark-current 0.05 --readout-noise 1.5 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc3.fits', doc ),
            requiredFile('./simulation image {0}/skipper_noisedc4 -N 0 --dark-current 0.05 --readout-noise 1.5 --ccd-shape 6000 1000 --horizontal-overscan 100 --vertical-overscan 50 --charge-gain 7.5'.format(folder),
                         'noisedc4.fits', doc ),
        )

        doc.frame('simulate noise and dc images',
            'combine images',
            requiredFile('./image superpose {0}/sim0.fits {0}/noisedc0.fits --output {0}/sim00.fits'.format(folder),
                         'sim00.fits', doc ),
            requiredFile('./image superpose {0}/sim1.fits {0}/noisedc1.fits --output {0}/sim11.fits'.format(folder),
                         'sim11.fits', doc ),
            requiredFile('./image superpose {0}/sim2.fits {0}/noisedc2.fits --output {0}/sim22.fits'.format(folder),
                         'sim22.fits', doc ),
            requiredFile('./image superpose {0}/sim3.fits {0}/noisedc3.fits --output {0}/sim33.fits'.format(folder),
                         'sim33.fits', doc ),
            requiredFile('./image superpose {0}/sim4.fits {0}/noisedc4.fits --output {0}/sim44.fits'.format(folder),
                         'sim44.fits', doc ),
        )

        doc.frame('simulate noise and dc images',
            'extract threshold 2e$^-$ = 15ADU',
            requiredFile('./catalog extract {0}/ext00 {0}/sim00.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext00.root', doc ),
            requiredFile('./catalog extract {0}/ext11 {0}/sim11.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext11.root', doc ),
            requiredFile('./catalog extract {0}/ext22 {0}/sim22.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext22.root', doc ),
            requiredFile('./catalog extract {0}/ext33 {0}/sim33.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext33.root', doc ),
            requiredFile('./catalog extract {0}/ext44 {0}/sim44.fits -t 15 -b 1 --ohdu -1 -v'.format(folder),
                         'ext44.root', doc ),
        )
