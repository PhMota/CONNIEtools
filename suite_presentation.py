from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *

folder = 'suite_presentation'

if not os.path.exists(folder):
    os.mkdir(folder)

def makefigure( code, file, doc, height=1 ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(file):
        subprocess.call( [code], shell=True )
    b = doc.figure( file, height=height, center=True )
    return [a,b]

with Timer('presentation'):
    with openBeamer( folder, 'analysis suite' ) as doc:
    
        doc.frame('simulation',
            'generates fits images',
            doc.itemize(
                'readout noise',
                'dark current',
                'gain',
                'modulation',
                'single hits'+
                doc.itemize(
                    'diffusion',
                    'charge efficiency',
                    'energy PDF',
                )
            )
        )
                
        doc.frame('image',
            'fits image manipulation',
            doc.itemize(
                'header',
                'display'+
                doc.itemize(
                    'overscan subtraction',
                    'vertical overscan subtraction',
                ),
                'extract'+
                doc.itemize(
                    'threshold',
                    'extra border pixels skirts',
                    'neighbor window shape'
                )
            )
        )
                
        doc.frame('catalog',
            'catalog manipulation',
            doc.itemize(
                'status',
                'addbranches'+
                doc.itemize(
                    'size likelihood',
                    'sigma and E fit'
                ),
                'match'+
                doc.itemize(
                    'multiple matching flag',
                    'original charge, depth, width',
                ),
                'histogram',
                'scatter',
                'ratio'
            )
        )


        doc.frame('simulating a blank image',
            'blank image with noise 12 ADU and DC 1e$^-$',
            *makefigure('./simulation image simulation1b_presentation --dark-current 1 --readout-noise 12 --ccd-shape 500 500 --horizontal-overscan 150 --vertical-overscan 150 --charge-gain 7.25 --pdf', 'simulation1b_presentation.pdf', doc, height=.6 )
        )
        
        doc.frame('display projection of the blank image',
            'projection along the horizontal axis shows active and overscan area',
            *makefigure('./image display simulation1b_presentation.fits --ohdu -1 --plot proj 0 --png', 'simulation1b_presentation.fits.proj_0.png', doc, height=.6 )
        )

        doc.frame('display projection of active area',
            'removing the overscan areas, the mean reflects the dc of 1e$^-$ with gain of 7.25 ADU/e$^-$',
            *makefigure('./image display simulation1b_presentation.fits --ohdu -1 --x-range 0 500 --y-range 150 650 --plot proj 0 a --no-max --no-min --png', 'simulation1b_presentation.fits.proj_0_a.png', doc, height=.6 )
        )

        doc.frame('extract hits from blank',
            'extraction with threshold of 60 ADU, hits appear from due to very high dark current plus noise',
            *makefigure('./catalog extract presentation_blank simulation1b_presentation.fits -t 60 -b 3 --ohdu -1 -v; ./catalog histogram presentation_blank.root --branch-selections E1 1 --branch-selections E0 1 --png --binsize 10 --factor 10 --output presentation_blankE1_10', 'presentation_blankE1_10.png', doc, height=.6 )
        )
        


        doc.frame('simulating 100 Cu hits with modulation',
            'simulate 100 Copper events with noise 12 ADU and DC 0.1e$^-$, user defined diffusion function, charge efficiency and modulations',
            *makefigure('./simulation image simulationCu2_presentation --dark-current 0.1 --readout-noise 12 --ccd-shape 500 500 --horizontal-overscan 150 --vertical-overscan 150 --charge-gain 7.25 --number-of-Cu-charges 100 --diffusion-function "sqrt(-258.817238*log1p(-0.000982*z))/15 if z < 670 else 0" --vertical-modulation-function "600*cos(y/pi/10)" --horizontal-modulation-function "-500*(x/100.-1)**2 if x<100 else 0" --depth-range 0 675 --charge-efficiency-function "1. if z < 600 else .9" --pdf --csv', 'simulationCu2_presentation.pdf', doc, height=.6 )
        )

        doc.frame('display projection of modulated image',
            'projection along the horizontal axis',
            *makefigure('./image display simulationCu2_presentation.fits --ohdu -1 --plot proj 0 --png', 'simulationCu2_presentation.fits.proj_0.png', doc, height=.6 )
        )
        
        doc.frame('display projection of modulated image',
            'projection along the vertical axis',
            *makefigure('./image display simulationCu2_presentation.fits --ohdu -1 --plot proj 1 --png', 'simulationCu2_presentation.fits.proj_1.png', doc, height=.6 )
        )        
    
        doc.frame('display projection of subtracted image',
            'horizontal and vertical overscan median subtraction',
            *makefigure('./image display simulationCu2_presentation.fits --ohdu -1 --plot proj 0 b --png --no-max --no-min --overscan-subtraction 500 650 --vertical-overscan-subtraction 0 150', 'simulationCu2_presentation.fits.proj_0_b.png', doc, height=.6 )
        )        
        
        doc.frame('display projection of subtracted image',
            'horizontal and vertical overscan median subtraction',
            *makefigure('./image display simulationCu2_presentation.fits --ohdu -1 --plot proj 1 c --png --no-max --no-min --overscan-subtraction 500 650 --vertical-overscan-subtraction 0 150 --x-range 0 500 --y-range 150 650', 'simulationCu2_presentation.fits.proj_1_c.png', doc, height=.6 )
        )
        
        doc.frame('extract hits from modulated Cu',
            'extracted events from subtracted image clearly shows Cu peak',
            *makefigure('./catalog extract presentationCu2 simulationCu2_presentation.fits -t 60 -b 1 --ohdu -1 -v --overscan-subtraction 500 650 --vertical-overscan-subtraction 0 150; ./catalog histogram presentationCu2.root --branch-selections E1 1 --branch-selections E0 1 --png --binsize 200 --factor 200 --output presentationCu2E1_g50', 'presentationCu2E1_g50.png', doc, height=.6 )
        )
        
        doc.frame('compute sigma fits',
            'scatter with the computed fit for the event width',
            *makefigure('./catalog addbranches presentationCu2.root fit -o presentationCu2_fit.root; ./catalog scatter presentationCu2_fit.root --x-range 5000 20000 --y-range 0 2 --branch-selections E0 xSigmafit 1 --branch-selections E1 xSigmafit 1 --png', 'presentationCu2_fit.root.scatter.xSigmafit.vs.E0.png', doc, height=.6 )
        )
        
        doc.frame('match with simulation',
            'comparison between reconstructed and simulated depths',
            *makefigure('./catalog match simulationCu2_presentation presentationCu2_fit.root -o presentationCu2_fit_matched.root -v --x-shift 150; ./catalog scatter presentationCu2_fit_matched.root --x-range 0 2 --y-range 0 2 --branch-selections sigmaSim sigmaSim 1 --branch-selections xSigmafit sigmaSim 1 --png', 'presentationCu2_fit_matched.root.scatter.sigmaSim.vs.sigmaSim.png', doc, height=.6 )
        )
        
        doc.frame('match with simulation',
            'comparison between reconstructed and simulated depths',
            *makefigure('./catalog scatter presentationCu2_fit_matched.root --x-range 0 2 --y-range -20000 20000 --branch-selections sigmaSim E0 1 --branch-selections sigmaSim ESim 1  --branch-selections sigmaSim \'E0-ESim\' 1 --png', 'presentationCu2_fit_matched.root.scatter.E0.vs.sigmaSim.png', doc, height=.6 )
        )
        
