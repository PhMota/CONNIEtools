from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *

folder = 'dc_study_nu_rate'

if not os.path.exists(folder):
    os.mkdir(folder)

def makefigure( code, file, doc, height=1 ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(folder+'/'+file):
        subprocess.call( [code], shell=True )
    b = doc.figure( file, height=height, center=True )
    return [a,b]

with Timer('presentation'):
    with openBeamer( folder, 'Dark Current study with neutrino rate\\\\https://github.com/PhMota/CONNIEtools/wiki' ) as doc:
        
        doc.frame('study',
            'hypothesis',
            doc.itemize(
                'removing n0$\lt=$2 will take care of all fake events and not affect neutrino spectrum',
            ),
            'methodology',
            doc.itemize(
                r'generate separate $\nu$, noise and dc images',
                r'dc+noise: confirm that all fake events are in n0$\lt=$2',
                r'$\nu$+noise: analyse the effect of this cut in the pure neutrino spectrum',
                r'$\nu$+dc+noise: analyse the effect of this cut in the pure neutrino + fake spectrum',
            )
        )
        
        doc.frame('simulate neutrino rate populations',
            r'25k $\nu$, dc 0 and noise 0 ADU',
            *makefigure('./simulation image {0}/dc_study_nu_rate_test3 -N 25000 --charge-pdf-table rate_neutrino.dat --dark-current 0 --readout-noise 0 --ccd-shape 4130 4120 --horizontal-overscan 450 --vertical-overscan 70 --charge-gain 7.25 --pdf'.format(folder), 'dc_study_nu_rate_test3.pdf', doc, height=.6 )
        )
        
        doc.frame('simulate neutrino rate populations',
            '$x$- and $y$- projections',
            doc.column(
                '\n'.join( makefigure('./image display {0}/dc_study_nu_rate_test3.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_nu_rate_test3.fits.proj_0.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./image display {0}/dc_study_nu_rate_test3.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_nu_rate_test3.fits.proj_1.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('extract neutrino events',
            'extract with threshold 1, full input spectrum',
            doc.column(
                '\n'.join( makefigure('./catalog extract {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_rate_test3.fits -t 1 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_nu_rate_test3.root --log --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2"  --branch-selections E0 "n0==1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_rate_test3.E0_g2'.format(folder), 'dc_study_nu_rate_test3.E0_g2.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog ratio {0}/dc_study_nu_rate_test3.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1  --branch-selections E0 "n0==1" 1 --png --binsize 50 --factor 50 --output {0}/dc_study_nu_rate_test3.ratio.E0_g2'.format(folder), 'dc_study_nu_rate_test3.ratio.E0_g2.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('neutrino spectra',
            'compute xSigmafit and match with simulations',
            doc.column(
                '\n'.join( makefigure('./catalog addbranches {0}/dc_study_nu_rate_test3.root fit -o {0}/dc_study_nu_rate_test3.fit.root; ./catalog scatter {0}/dc_study_nu_rate_test3.fit.root --branch-selections E0 xSigmafit 1 --branch-selections E0 xSigmafit "n0>2" --branch-selections E0 xSigmafit "n0==2" --branch-selections E0 xSigmafit "n0==1" --png --output {0}/dc_study_nu_rate_test3.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_rate_test3.scatter3.E0_g.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_rate_test3.fit.root -o {0}/dc_study_nu_rate_test3.match.fit.root; ./catalog scatter {0}/dc_study_nu_rate_test3.match.fit.root --branch-selections ESim sigmaSim 1 --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0==2" --branch-selections ESim sigmaSim "n0==1" --png --output {0}/dc_study_nu_rate_test3.match.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_rate_test3.match.scatter3.E0_g.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('conclusion',
            doc.itemize(
                'neutrino population generated using a table with $E$ and $dR/dE$ as a command-line input',
                'the population of single and double pixel neutrinos is very small and confined to the surface, as expected',
                'however, these are pure events extracted with no threshold'
            )
        )

        #doc.frame('introducing noise', '')
            
        doc.frame('only noise',
            'image with only noise',
            *makefigure('./simulation image {0}/dc_study_only_noise -N 0 --dark-current 0 --readout-noise 15 --ccd-shape 4130 4120 --horizontal-overscan 450 --vertical-overscan 70 --charge-gain 7.25 --pdf'.format(folder), 'dc_study_only_noise.pdf', doc, height=.6 )                
        )

        doc.frame('only noise',
            'projections',
            doc.column(
                '\n'.join( makefigure('./image display {0}/dc_study_only_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_only_noise.fits.proj_0.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./image display {0}/dc_study_only_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_only_noise.fits.proj_1.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )


        doc.frame('conclusion',
            doc.itemize(
                'a noise of 15 ADU alone is not enough to go over the 60ADU threshold',
            )
        )

            
        doc.frame('only dark current',
            'image with only DC',
            *makefigure('./simulation image {0}/dc_study_only_dc -N 0 --dark-current 0.1 --readout-noise 0 --ccd-shape 4130 4120 --horizontal-overscan 450 --vertical-overscan 70 --charge-gain 7.25 --pdf'.format(folder), 'dc_study_only_dc.pdf', doc, height=.6 )                
        )

        doc.frame('only dark current',
            'projections',
            doc.column(
                '\n'.join( makefigure('./image display {0}/dc_study_only_dc.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_only_dc.fits.proj_0.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./image display {0}/dc_study_only_dc.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_only_dc.fits.proj_1.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('conclusion',
            doc.itemize(
                'two, three and four electron pixels are present, but DC alone cannot pass the 60 ADU threshold either',
            )
        )

        #doc.frame('dark current and noise', '')
                
        doc.frame('dark current and noise',
            'superpose noise 12 ADU and DC 0.1',
            doc.column(
                '\n'.join( makefigure('./image superpose {0}/dc_study_only_noise.fits {0}/dc_study_only_dc.fits --output {0}/dc_study_dc_noise.fits; ./image display {0}/dc_study_dc_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_dc_noise.fits.proj_0.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./image display {0}/dc_study_dc_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_dc_noise.fits.proj_1.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
            )
        )

        doc.frame('dark current and noise',
            'fake events with threshold 60 ADU',
            doc.column(
                '\n'.join( makefigure('./catalog extract {0}/dc_study_dc_noise {0}/dc_study_dc_noise.fits -t 60 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_dc_noise.root --log  --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2" --branch-selections E0 "n0==1" --png --binsize 5 --factor 5 --output {0}/dc_study_dc_noise.E0_g2'.format(folder), 'dc_study_dc_noise.E0_g2.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog ratio {0}/dc_study_dc_noise.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1 --branch-selections E0 "n0==1" 1 --png --binsize 5 --factor 5 --output {0}/dc_study_dc_noise.ratio.E0_g2'.format(folder), 'dc_study_dc_noise.ratio.E0_g2.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )
        
        doc.frame('conclusion',
            doc.itemize(
                'by adding DC 12 ADU and noise 0.1, fake events start showing up (around 30 events) and they all fall in the n0==1 and n0==2 range',
            )
        )


                
        doc.frame('neutrinos and noise',
            'superpose neutrinos and noise',
            doc.column(
                '\n'.join( makefigure('./image superpose {0}/dc_study_only_noise.fits {0}/dc_study_nu_rate_test3.fits --output {0}/dc_study_nu_noise.fits; ./image display {0}/dc_study_nu_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_nu_noise.fits.proj_0.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./image display {0}/dc_study_nu_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_nu_noise.fits.proj_1.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
            )
        )

        doc.frame('neutrinos and noise',
            'spectrum of neutrinos with noise, threshold 60',
            doc.column(
                '\n'.join( makefigure('./catalog extract {0}/dc_study_nu_noise {0}/dc_study_nu_noise.fits -t 60 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_nu_noise.root --log  --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2" --branch-selections E0 "n0==1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_noise.E0_g2'.format(folder), 'dc_study_nu_noise.E0_g2.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog ratio {0}/dc_study_nu_noise.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1 --branch-selections E0 "n0==1" 1 --png --binsize 50 --factor 50 --output {0}/dc_study_nu_noise.ratio.E0_g2'.format(folder), 'dc_study_nu_noise.ratio.E0_g2.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('neutrinos and noise',
            'compute xSigmafit and match the catalog',
            doc.column(
                '\n'.join( makefigure('./catalog addbranches {0}/dc_study_nu_noise.root fit -o {0}/dc_study_nu_noise.fit.root; ./catalog scatter {0}/dc_study_nu_noise.fit.root --branch-selections E0 xSigmafit 1 --branch-selections E0 xSigmafit "n0>2" --branch-selections E0 xSigmafit "n0==2" --branch-selections E0 xSigmafit "n0==1" --png --output {0}/dc_study_nu_noise.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_noise.scatter3.E0_g.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_noise.fit.root -o {0}/dc_study_nu_noise.match.fit.root; ./catalog scatter {0}/dc_study_nu_noise.match.fit.root  --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0==2" --branch-selections ESim sigmaSim "n0==1" --png --output {0}/dc_study_nu_noise.match.scatter3.E0_ga --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_noise.match.scatter3.E0_ga.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )
 
        doc.frame('conclusion',
            doc.itemize(
                'by adding the neutrino spectrum to the noise of 12 ADU and applying a threshold of 60 ADU, the population of n0==1 and n0==2 neutrinos is huge',
                'however, it concentrates in the surface and on very low energy events, as expected'
            )
        )

        #doc.frame('neutrinos dark current and noise', '')
                
        doc.frame('neutrinos, dark current and noise',
            'superpose neutrinos, noise and darkcurrent',
            doc.column(
                '\n'.join( makefigure('./image superpose {0}/dc_study_dc_noise.fits {0}/dc_study_nu_rate_test3.fits --output {0}/dc_study_nu_dc_noise.fits; ./image display {0}/dc_study_nu_dc_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_nu_dc_noise.fits.proj_0.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./image display {0}/dc_study_nu_dc_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_nu_dc_noise.fits.proj_1.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
            )
        )

        doc.frame('neutrinos, dark current and noise',
            'spectra of neutrinos with dark current and noise',
            doc.column(
                '\n'.join( makefigure('./catalog extract {0}/dc_study_nu_dc_noise {0}/dc_study_nu_dc_noise.fits -t 60 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_nu_dc_noise.root --log --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2" --branch-selections E0 "n0==1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.E0_g2'.format(folder), 'dc_study_nu_dc_noise.E0_g2.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog ratio {0}/dc_study_nu_dc_noise.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1 --branch-selections E0 "n0==1" 1 --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.ratio.E0_g2'.format(folder), 'dc_study_nu_dc_noise.ratio.E0_g2.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )
            
        doc.frame('neutrinos, dark current and noise',
            'compute xSigmafit and match the catalog',
            doc.column(
                '\n'.join( makefigure('./catalog addbranches {0}/dc_study_nu_dc_noise.root fit -o {0}/dc_study_nu_dc_noise.fit.root; ./catalog scatter {0}/dc_study_nu_dc_noise.fit.root --branch-selections E0 xSigmafit 1 --branch-selections E0 xSigmafit "n0>2" --branch-selections E0 xSigmafit "n0==2" --branch-selections E0 xSigmafit "n0==1" --png --output {0}/dc_study_nu_dc_noise.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.scatter3.E0_g.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_dc_noise.fit.root -o {0}/dc_study_nu_dc_noise.match.fit.root; ./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root  --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0==2" --branch-selections ESim sigmaSim "n0==1" --png --output {0}/dc_study_nu_dc_noise.match.scatter3.E0_ga --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.match.scatter3.E0_ga.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )
            
        doc.frame('neutrinos, dark current and noise',
            'compute xSigmafit and match the catalog',
            doc.column(
                '\n'.join( makefigure('./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root --branch-selections E0 xSigmafit "1" --branch-selections E0 xSigmafit "idSim==5 and n0>2" --branch-selections E0 xSigmafit "idSim==5 and n0==2" --branch-selections E0 xSigmafit "idSim==5 and n0==1" --png --output {0}/dc_study_nu_dc_noise.match.scatter3.E0_id --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.match.scatter3.E0_id.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root --branch-selections E0 xSigmafit "idSim==0" --branch-selections E0 xSigmafit "idSim==0 and n0>2" --branch-selections E0 xSigmafit "idSim==0 and n0==2" --branch-selections E0 xSigmafit "idSim==0 and n0==1" --png --output {0}/dc_study_nu_dc_noise.match.scatter3.E0_id0b --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.match.scatter3.E0_id0b.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )
            
        doc.frame('neutrinos, dark current and noise',
            'spectra of neutrinos with dark current and noise',
            doc.column(
                '\n'.join( makefigure('./catalog histogram {0}/dc_study_nu_dc_noise.fit.root --log --branch-selections E0 1 --branch-selections E0 "xSigmafit>.3 and n0>2" --branch-selections E0 "xSigmafit>.3 and n0==2" --branch-selections E0 "xSigmafit>.3 and n0==1" --branch-selections E0 "xSigmafit>.3 and n0>1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.E0_g2sigma2'.format(folder), 'dc_study_nu_dc_noise.E0_g2sigma2.png', doc, height=.4 ) )
                ,
                '\n'.join( makefigure('./catalog ratio {0}/dc_study_nu_dc_noise.fit.root --branch-selections E0 "xSigmafit>.3" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0>2" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0==2" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0==1" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0>1" "xSigmafit>.3" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.ratio.E0_g2sigma2'.format(folder), 'dc_study_nu_dc_noise.ratio.E0_g2sigma2.png', doc, height=.4 ) )
                ,
                widths = [.4,.4]
                )
        )
        
        doc.frame('conclusion',
            doc.itemize(
                'by adding the neutrino spectrum to the DC 0.1 and noise 12 ADU and applying a threshold of 60 ADU, the population of n0==1 and n0==2 neutrinos is huge',
                'removing n0==1 and n0==2 succesfully remove all fake events, but it is an overkill and removes MUCH MORE neutrinos',
                'on the other hand, removing surface neutrinos can be dealt with since we will apply a depth cut in the analysis',
                'disregarding the surface events, we end up with a much higher efficiency from 200 ADU up removing n0==1 and from 500 ADU up removing n0==2 aswell',
            )
        )
