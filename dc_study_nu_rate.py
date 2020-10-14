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
    with openBeamer( folder, 'Dark Current study with neutrino rate\\\\https://github.com/PhMota/CONNIEtools/wiki' ) as doc:
        
        doc.par('\n\part{''}\n\\frame{\partpage}\n')

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
                makefigure('./image display {0}/dc_study_nu_rate_test3.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_nu_rate_test3.fits.proj_0.png', doc, height=.4 )
                ,
                makefigure('./image display {0}/dc_study_nu_rate_test3.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_nu_rate_test3.fits.proj_1.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('extract neutrino events',
            'extract with threshold 1, full input spectrum',
            doc.column(
                makefigure('./catalog extract {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_rate_test3.fits -t 1 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_nu_rate_test3.root --log --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2"  --branch-selections E0 "n0==1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_rate_test3.E0_g2'.format(folder), 'dc_study_nu_rate_test3.E0_g2.png', doc, height=.4 )
                ,
                makefigure('./catalog ratio {0}/dc_study_nu_rate_test3.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1  --branch-selections E0 "n0==1" 1 --png --binsize 50 --factor 50 --output {0}/dc_study_nu_rate_test3.ratio.E0_g2'.format(folder), 'dc_study_nu_rate_test3.ratio.E0_g2.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('neutrino spectra',
            'compute xSigmafit and match with simulations',
            doc.column(
                makefigure('./catalog addbranches {0}/dc_study_nu_rate_test3.root fit -o {0}/dc_study_nu_rate_test3.fit.root; ./catalog scatter {0}/dc_study_nu_rate_test3.fit.root --branch-selections E0 xSigmafit 1 --branch-selections E0 xSigmafit "n0>2" --branch-selections E0 xSigmafit "n0==2" --branch-selections E0 xSigmafit "n0==1" --png --output {0}/dc_study_nu_rate_test3.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_rate_test3.scatter3.E0_g.png', doc, height=.4 )
                ,
                makefigure('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_rate_test3.fit.root -o {0}/dc_study_nu_rate_test3.match.fit.root; ./catalog scatter {0}/dc_study_nu_rate_test3.match.fit.root --branch-selections ESim sigmaSim 1 --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0==2" --branch-selections ESim sigmaSim "n0==1" --png --output {0}/dc_study_nu_rate_test3.match.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_rate_test3.match.scatter3.E0_g.png', doc, height=.4 )
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
            makefigure('./simulation image {0}/dc_study_only_noise -N 0 --dark-current 0 --readout-noise 15 --ccd-shape 4130 4120 --horizontal-overscan 450 --vertical-overscan 70 --charge-gain 7.25 --pdf'.format(folder), 'dc_study_only_noise.pdf', doc, height=.6 )                
        )

        doc.frame('only noise',
            'projections',
            doc.column(
                makefigure('./image display {0}/dc_study_only_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_only_noise.fits.proj_0.png', doc, height=.4 )
                ,
                makefigure('./image display {0}/dc_study_only_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_only_noise.fits.proj_1.png', doc, height=.4 )
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
                makefigure('./image display {0}/dc_study_only_dc.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_only_dc.fits.proj_0.png', doc, height=.4 )
                ,
                makefigure('./image display {0}/dc_study_only_dc.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_only_dc.fits.proj_1.png', doc, height=.4 )
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
            'superpose noise 15 ADU and DC 0.1',
            doc.column(
                makefigure('./image superpose {0}/dc_study_only_noise.fits {0}/dc_study_only_dc.fits --output {0}/dc_study_dc_noise.fits; ./image display {0}/dc_study_dc_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_dc_noise.fits.proj_0.png', doc, height=.4 )
                ,
                makefigure('./image display {0}/dc_study_dc_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_dc_noise.fits.proj_1.png', doc, height=.4 )
                ,
                widths = [.4,.4]
            )
        )

        doc.frame('dark current and noise',
            'fake events with threshold 60 ADU',
            doc.column(
                makefigure('./catalog extract {0}/dc_study_dc_noise {0}/dc_study_dc_noise.fits -t 60 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_dc_noise.root --log  --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2" --branch-selections E0 "n0==1" --png --binsize 5 --factor 5 --output {0}/dc_study_dc_noise.E0_g2'.format(folder), 'dc_study_dc_noise.E0_g2.png', doc, height=.4 )
                ,
                makefigure('./catalog ratio {0}/dc_study_dc_noise.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1 --branch-selections E0 "n0==1" 1 --png --binsize 5 --factor 5 --output {0}/dc_study_dc_noise.ratio.E0_g2'.format(folder), 'dc_study_dc_noise.ratio.E0_g2.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )
        
        doc.frame('conclusion',
            doc.itemize(
                'by adding DC 0.1 and noise 15 ADU, fake events start showing up (around 30 events) and they all fall in the n0==1 and n0==2 range',
            )
        )


                
        doc.frame('neutrinos and noise',
            'superpose neutrinos and noise',
            doc.column(
                makefigure('./image superpose {0}/dc_study_only_noise.fits {0}/dc_study_nu_rate_test3.fits --output {0}/dc_study_nu_noise.fits; ./image display {0}/dc_study_nu_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_nu_noise.fits.proj_0.png', doc, height=.4 )
                ,
                makefigure('./image display {0}/dc_study_nu_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_nu_noise.fits.proj_1.png', doc, height=.4 )
                ,
                widths = [.4,.4]
            )
        )

        doc.frame('neutrinos and noise',
            'spectrum of neutrinos with noise, threshold 60',
            doc.column(
                makefigure('./catalog extract {0}/dc_study_nu_noise {0}/dc_study_nu_noise.fits -t 60 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_nu_noise.root --log  --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2" --branch-selections E0 "n0==1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_noise.E0_g2'.format(folder), 'dc_study_nu_noise.E0_g2.png', doc, height=.4 )
                ,
                makefigure('./catalog ratio {0}/dc_study_nu_noise.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1 --branch-selections E0 "n0==1" 1 --png --binsize 50 --factor 50 --output {0}/dc_study_nu_noise.ratio.E0_g2'.format(folder), 'dc_study_nu_noise.ratio.E0_g2.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('neutrinos and noise',
            'compute xSigmafit and match the catalog',
            doc.column(
                makefigure('./catalog addbranches {0}/dc_study_nu_noise.root fit -o {0}/dc_study_nu_noise.fit.root; ./catalog scatter {0}/dc_study_nu_noise.fit.root --branch-selections E0 xSigmafit 1 --branch-selections E0 xSigmafit "n0>2" --branch-selections E0 xSigmafit "n0==2" --branch-selections E0 xSigmafit "n0==1" --png --output {0}/dc_study_nu_noise.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_noise.scatter3.E0_g.png', doc, height=.4 )
                ,
                makefigure('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_noise.fit.root -o {0}/dc_study_nu_noise.match.fit.root; ./catalog scatter {0}/dc_study_nu_noise.match.fit.root  --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0==2" --branch-selections ESim sigmaSim "n0==1" --png --output {0}/dc_study_nu_noise.match.scatter3.E0_ga --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_noise.match.scatter3.E0_ga.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )
 
        doc.frame('conclusion',
            doc.itemize(
                'by adding the neutrino spectrum to the noise of 15 ADU and applying a threshold of 60 ADU, the population of n0==1 and n0==2 neutrinos is huge',
                'however, it concentrates in the surface and on very low energy events, as expected'
            )
        )

        #doc.frame('neutrinos dark current and noise', '')
                
        doc.frame('neutrinos, dark current and noise',
            'superpose neutrinos, noise and darkcurrent',
            doc.column(
                makefigure('./image superpose {0}/dc_study_dc_noise.fits {0}/dc_study_nu_rate_test3.fits --output {0}/dc_study_nu_dc_noise.fits; ./image display {0}/dc_study_nu_dc_noise.fits --ohdu -1 --plot proj 0 --png'.format(folder), 'dc_study_nu_dc_noise.fits.proj_0.png', doc, height=.4 )
                ,
                makefigure('./image display {0}/dc_study_nu_dc_noise.fits --ohdu -1 --plot proj 1 --png'.format(folder), 'dc_study_nu_dc_noise.fits.proj_1.png', doc, height=.4 )
                ,
                widths = [.4,.4]
            )
        )

        doc.frame('neutrinos, dark current and noise',
            'spectra of neutrinos with dark current and noise',
            doc.column(
                makefigure('./catalog extract {0}/dc_study_nu_dc_noise {0}/dc_study_nu_dc_noise.fits -t 60 -b 1 --ohdu -1 -v; ./catalog histogram {0}/dc_study_nu_dc_noise.root --log --branch-selections E0 1 --branch-selections E0 "n0>2" --branch-selections E0 "n0==2" --branch-selections E0 "n0==1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.E0_g2'.format(folder), 'dc_study_nu_dc_noise.E0_g2.png', doc, height=.4 )
                ,
                makefigure('./catalog ratio {0}/dc_study_nu_dc_noise.root --branch-selections E0 1 1 --branch-selections E0 "n0>2" 1 --branch-selections E0 "n0==2" 1 --branch-selections E0 "n0==1" 1 --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.ratio.E0_g2'.format(folder), 'dc_study_nu_dc_noise.ratio.E0_g2.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )
            
        doc.frame('neutrinos, dark current and noise',
            'compute xSigmafit and match the catalog',
            doc.column(
                makefigure('./catalog addbranches {0}/dc_study_nu_dc_noise.root fit -o {0}/dc_study_nu_dc_noise.fit.root; ./catalog scatter {0}/dc_study_nu_dc_noise.fit.root --branch-selections E0 xSigmafit 1 --branch-selections E0 xSigmafit "n0>2" --branch-selections E0 xSigmafit "n0==2" --branch-selections E0 xSigmafit "n0==1" --png --output {0}/dc_study_nu_dc_noise.scatter3.E0_g --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.scatter3.E0_g.png', doc, height=.4 )
                ,
                makefigure('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_dc_noise.fit.root -o {0}/dc_study_nu_dc_noise.match.fit.root; ./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root  --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0>2" --branch-selections ESim sigmaSim "n0==2" --branch-selections ESim sigmaSim "n0==1" --png --output {0}/dc_study_nu_dc_noise.match.scatter3.E0_ga --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.match.scatter3.E0_ga.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )
            
        doc.frame('neutrinos, dark current and noise',
            'compute xSigmafit and match the catalog',
            doc.column(
                makefigure('./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root --branch-selections E0 xSigmafit "1" --branch-selections E0 xSigmafit "idSim==5 and n0>2" --branch-selections E0 xSigmafit "idSim==5 and n0==2" --branch-selections E0 xSigmafit "idSim==5 and n0==1" --png --output {0}/dc_study_nu_dc_noise.match.scatter3.E0_id --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.match.scatter3.E0_id.png', doc, height=.4 )
                ,
                makefigure('./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root --branch-selections E0 xSigmafit "idSim==0" --branch-selections E0 xSigmafit "idSim==0 and n0>2" --branch-selections E0 xSigmafit "idSim==0 and n0==2" --branch-selections E0 xSigmafit "idSim==0 and n0==1" --png --output {0}/dc_study_nu_dc_noise.match.scatter3.E0_id0b --x-range 0 1500 --y-range 0 2'.format(folder), 'dc_study_nu_dc_noise.match.scatter3.E0_id0b.png', doc, height=.4 ) 
                ,
                widths = [.4,.4]
                )
        )
            
        doc.frame('neutrinos, dark current and noise',
            'spectra of neutrinos with dark current and noise',
            doc.column(
                makefigure('./catalog histogram {0}/dc_study_nu_dc_noise.fit.root --log --branch-selections E0 1 --branch-selections E0 "xSigmafit>.3 and n0>2" --branch-selections E0 "xSigmafit>.3 and n0==2" --branch-selections E0 "xSigmafit>.3 and n0==1" --branch-selections E0 "xSigmafit>.3 and n0>1" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.E0_g2sigma2'.format(folder), 'dc_study_nu_dc_noise.E0_g2sigma2.png', doc, height=.4 )
                ,
                makefigure('./catalog ratio {0}/dc_study_nu_dc_noise.fit.root --branch-selections E0 "xSigmafit>.3" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0>2" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0==2" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0==1" "xSigmafit>.3" --branch-selections E0 "xSigmafit>.3 and n0>1" "xSigmafit>.3" --png --binsize 50 --factor 50 --output {0}/dc_study_nu_dc_noise.ratio.E0_g2sigma2'.format(folder), 'dc_study_nu_dc_noise.ratio.E0_g2sigma2.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )
        
        doc.frame('conclusion',
            doc.itemize(
                'by adding the neutrino spectrum to the DC 0.1 and noise 15 ADU and applying a threshold of 60 ADU, the population of n0==1 and n0==2 neutrinos is huge',
                'removing n0==1 and n0==2 succesfully remove all fake events, but it is an overkill and removes MUCH MORE neutrinos',
                'on the other hand, removing surface neutrinos can be dealt with since we will apply a depth cut in the analysis',
                'disregarding the surface events, we end up with a much higher efficiency from 200 ADU up removing n0==1 and from 500 ADU up removing n0==2 aswell',
            )
        )
            
        doc.par('\n\part{''}\n\\frame{\partpage}\n')
        
        doc.frame('fake events',
            doc.itemize(
                '500 images',
                'noise 13.4 ADU',
                'DC 0.043',
                'gain 7 ADU/e$^-$ or 1.87 ADU/eV'
            )
        )
        
        doc.frame('fake events',
            "Javier's extraction algorithm",
            doc.column(
                makefigure( 
                    './catalog histogram /home/bonifazi/public/sims_13.4_0.043_7.0.root --log --branch-selections E0 1 --output {0}/sims_13.4_0.043_7.0.carla.histo --png'.format( folder ),
                    'sims_13.4_0.043_7.0.carla.histo.png', doc, height=.4
                ),
                makefigure( 
                    './catalog scatter /home/bonifazi/public/sims_13.4_0.043_7.0.root --branch-selections E0 n0 1 --x-range 55 135 --y-range 0.5 2.5 --output {0}/sims_13.4_0.043_7.0.carla.scatter --png'.format( folder ),
                    'sims_13.4_0.043_7.0.carla.scatter.png', doc, height=.4
                ),
                widths = [.4,.4]
            )
        )

        doc.frame('fake events',
            requiredFile( 
                './catalog extract {}/sims_13.4_0.043_7.0_b0 "/share/storage2/connie/users/bonifazi/SimsOHDUs/blanks/sims_13.4_0.043_7.0/*.fits" -t 60 -b 0 --ohdu -1 -v'.format(folder), 
                'sims_13.4_0.043_7.0_b0.root'.format(folder), doc 
            ),
            'my extraction algorithm using no skirts',
            doc.column(
                makefigure( 
                    './catalog histogram {0}/sims_13.4_0.043_7.0_b0.root --log --branch-selections E0 1 --output {0}/sims_13.4_0.043_7.0_b0.histo2 --png'.format( folder ),
                    'sims_13.4_0.043_7.0_b0.histo2.png', doc, height=.4
                ),
                makefigure( 
                    './catalog scatter {0}/sims_13.4_0.043_7.0_b0.root --branch-selections E0 n0 1 --x-range 55 135 --y-range 0.5 2.5 --output {0}/sims_13.4_0.043_7.0_b0.scatter2 --png'.format( folder ),
                    'sims_13.4_0.043_7.0_b0.scatter2.png', doc, height=.4
                ),
                widths = [.4,.4]
            )
        )

        doc.frame('fake events',
            requiredFile( 
                './catalog extract {}/sims_13.4_0.043_7.0 "/share/storage2/connie/users/bonifazi/SimsOHDUs/blanks/sims_13.4_0.043_7.0/*.fits" -t 60 -b 0 --ohdu -1 -v'.format(folder), 
                'sims_13.4_0.043_7.0.root'.format(folder), doc 
            ),
            'my extraction algorithm using level 1 skirts',
            doc.column(
                makefigure( 
                    './catalog histogram {0}/sims_13.4_0.043_7.0.root --log --branch-selections E0 1 --output {0}/sims_13.4_0.043_7.0.histo3 --png'.format( folder ),
                    'sims_13.4_0.043_7.0.histo3.png', doc, height=.4
                ),
                makefigure( 
                    './catalog scatter {0}/sims_13.4_0.043_7.0.root --branch-selections E0 n0 1 --x-range 55 135 --y-range 0.5 2.5 --output {0}/sims_13.4_0.043_7.0.scatter3 --png'.format( folder ),
                    'sims_13.4_0.043_7.0.scatter3.png', doc, height=.4
                ),
                widths = [.4,.4]
            )
        )
        
        doc.frame('fake events',
            'all spectra together -- count in $y$-axis',
            makefigure( 
                './catalog histogram /home/bonifazi/public/sims_13.4_0.043_7.0.root E0 1 {0}/sims_13.4_0.043_7.0_b0.root E0 1 {0}/sims_13.4_0.043_7.0.root E0 1 --log --binsize 3 --output {0}/sims_13.4_0.043_7.0.all.histo.range1d --png'.format( folder ),
                'sims_13.4_0.043_7.0.all.histo.range1d.png', doc, height=.6
            )
        )
        
                
        doc.frame('conclusion',
            doc.itemize(
                "the events extracted using my algorithm differs from Javier's both in the minimum and maximum value of the double-pixel population",
                "after close analysis on Javier's code, it is clear that the 1st pixel -- the trigger pixel -- receives a differente treatment (controlled by the seed parameter $=4$) compared to the subsequently ones added to the cluster (addThr$=3.4$)",
                "this explains -- but not justifies -- the fact that Javier's code has a minimum energy for the double-pixel population of 111 ADU and not the expected $120 = 2\\times60$",
                "in my algorithm, the overlapping of the skirts is responsible for the merging of close pixels and produces higher maximal energy for the double-pixel population for level 1 skirt extraction",
                #'proposed cut (n0==1 && E0<0.45) || (n0==2 && E0<70eV) || n0>2'
            )
        )
        
        
        doc.frame('neutrinos thr60',
            'proposed cut (n0==1 and E0$\lt$=100) or (n0==2 and E0$\lt$=130) or n0$\gt$2',
            doc.column(
                makefigure('./catalog scatter {0}/dc_study_nu_dc_noise.fit.root --branch-selections E0 xSigmafit 1 --x-range 0 200 --y-range 0 2 --png --output {0}/dc_study_nu_dc_noise.scatter4aa.E0'.format(folder), 
                'dc_study_nu_dc_noise.scatter4aa.E0.png', doc, height=.4 )
                ,
                makefigure('./catalog scatter {0}/dc_study_nu_dc_noise.match.fit.root --branch-selections E0 xSigmafit "idSim==5 && ( (n0==1 && E0<=100) || (n0==2 && E0<=130) )" --branch-selections E0 xSigmafit "idSim==0 && ( (n0==1 && E0<=100) || (n0==2 && E0<=130) )" --branch-selections E0 xSigmafit "idSim==5 && ( (n0==1 && E0>100) || (n0==2 && E0>130) )" --x-range 0 200 --y-range 0 2 --png --output {0}/dc_study_nu_dc_noise.scatter4bc.E0'.format(folder), 
                           'dc_study_nu_dc_noise.scatter4bc.E0.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )

        doc.frame('extraction with 40 ADU',
            requiredFile('./catalog extract {0}/dc_study_nu_dc_noiset40 {0}/dc_study_nu_dc_noise.fits -t 40 -b 1 --ohdu -1 -v'.format(folder),
                         'dc_study_nu_dc_noiset40.root', doc )
            ,
            requiredFile('./catalog addbranches {0}/dc_study_nu_dc_noiset40.root fit -o {0}/dc_study_nu_dc_noiset40.fit.root'.format(folder),
                         'dc_study_nu_dc_noiset40.fit.root', doc )
            ,
            requiredFile('./catalog match {0}/dc_study_nu_rate_test3 {0}/dc_study_nu_dc_noiset40.fit.root -o {0}/dc_study_nu_dc_noiset40.match.fit.root'.format(folder),
                         'dc_study_nu_dc_noiset40.match.fit.root', doc )
        )
        doc.frame('neutrinos thr40',
            'proposed cut (n0==1 and E0$\lt$=100) or (n0==2 and E0$\lt$=130) or n0$\gt$2',
            doc.column(
                makefigure('./catalog scatter {0}/dc_study_nu_dc_noiset40.match.fit.root --branch-selections E0 xSigmafit 1 --x-range 0 200 --y-range 0 2 --png --output {0}/dc_study_nu_dc_noiset40.match.scatter.E0'.format(folder), 
                'dc_study_nu_dc_noiset40.match.scatter.E0.png', doc, height=.4 )
                ,
                makefigure('./catalog scatter {0}/dc_study_nu_dc_noiset40.match.fit.root --branch-selections E0 xSigmafit "idSim==5 && ( (n0==1 && E0<=100) || (n0==2 && E0<=130) )" --branch-selections E0 xSigmafit "idSim==0 && ( (n0==1 && E0<=100) || (n0==2 && E0<=130) )" --branch-selections E0 xSigmafit "idSim==5 && ( (n0==1 && E0>100) || (n0==2 && E0>130) )" --x-range 0 200 --y-range 0 2 --png --output {0}/dc_study_nu_dc_noiset40.match.scatter3.E0'.format(folder), 'dc_study_nu_dc_noiset40.match.scatter3.E0.png', doc, height=.4 )
                ,
                widths = [.4,.4]
                )
        )
