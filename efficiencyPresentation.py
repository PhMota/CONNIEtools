from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
from numpy import *

folder = 'efficiencyPresentation'

if not os.path.exists(folder):
    os.mkdir(folder)

#def projections_slide( name, cmd ):
    #scale = .3
    #contents = ['projections of the raw image', doc.code( cmd, 'Bash') ]
    #contents += [ doc.figure( name+'/median_o3_{}.pdf'.format(key), scale=scale ) for key in ['dataProj1', 'dataProj0', 'biasProj', 'vbiasProj' ] ]
    #doc.frame( *contents )

def run(cmd):
    subprocess.call( [cmd], shell=True )

def require( fname, cmd, redo = False ):
    if not os.path.exists(fname) or redo:
        print( 'running' )
        print( '>', cmd )
        run(cmd)
    if not os.path.exists(fname):
        print( 'the command' )
        print( '>', cmd )
        print( 'did not generate the file', fname )
        exit(0)
    print( colored('found file', 'green'), fname )
    return True

class Scatter:
    def __init__( self, catalog, x, y, z=None, selection=None, _xrange=None, _yrange=None, fit=None ):
        if z == None:
            if fit is None:
                self.path = '{}.scatter.{}.vs.{}.png'.format( catalog, y.replace('/','_'), x )
                self.cmd = '''
                        ./catalog scatter {} --branch-selections '{}' '{}' '{}' --x-range {} --y-range {} --errorbar --png
                        '''.format(catalog, x, y, selection, _xrange, _yrange)
            else:
                self.path = '{}.scatter.{}.vs.{}_fit{}.png'.format( catalog, y.replace('/','_'), x, ''.join(fit) )
                self.cmd = '''
                        ./catalog scatter {} --branch-selections '{}' '{}' '{}' --x-range {} --y-range {} --errorbar --png --fit {}
                        '''.format(catalog, x, y, selection, _xrange, _yrange, ' '.join(fit) )
        else:
            if fit is None:
                self.path = '{}.scatter.{}.vs.{}.png'.format( catalog, y.replace('/','_'), x )
                self.cmd = '''
                        ./catalog scatter {} --branch-selections '{}' '{}' '{}' '{}' --x-range {} --y-range {} --errorbar --png
                        '''.format(catalog, x, y, z, selection, _xrange, _yrange)
            else:
                self.path = '{}.scatter.{}.vs.{}_fit{}.png'.format( catalog, y.replace('/','_'), x, ''.join(fit) )
                self.cmd = '''
                        ./catalog scatter {} --branch-selections '{}' '{}' '{}' '{}' --x-range {} --y-range {} --errorbar --png --fit {}
                        '''.format(catalog, x, y, z, selection, _xrange, _yrange, ' '.join(fit) )
            
            
    def get_path(self, redo=False):
        require( self.path, self.cmd, redo=redo )
        return self.path

class Catalog:
    def __init__( self, thresholds, b, N, number_of_images, charge_range, image_modes, rn, dc ):
        self.basename = 'newcatalogs_{}N{}max{}rn{}dc{}'.format( number_of_images, N, charge_range[1], rn, dc )
        self.cmd = '''
              ./catalog simulation {} -t {} -b {} -N {} --number-of-images {} --charge-range {} --image-modes {} -rn {} -dc {}
              '''.format( self.basename, ' '.join(thresholds), b, N, number_of_images, ' '.join(charge_range), ' '.join(image_modes), rn, dc )
        
    def scatter(self, image_mode, thr, x, y, z=None, selection=None, _xrange=None, _yrange=None, fit=None, redo=False):
        self.catalogpath = '{0}/{0}_bin{1}t{2}.root'.format( self.basename, image_mode, thr )
        require( self.catalogpath, self.cmd )
        return Scatter( self.catalogpath, x, y, z, selection, _xrange, _yrange, fit=fit ).get_path(redo=redo)

with Timer('presentation'):
    with openBeamer( folder, '' ) as doc:

        doc.par(r'\def\vr{\tt\color{gray}}')
        doc.frame('matching simulation and extraction', 
            doc.itemize(
                r'get x- and y-positions of the simulated charges and match with {\vr xPix} and {\vr yPix}',
                r'include columns to catalog\\ {\vr xSim, ySim, zSim, qSim, ESim, sigmaSim, idSim}',
                r'{\vr idSim==1}: matched hits',
                r'{\vr idSim==2}: mulitple matched hits; {\vr qSim} and {\vr ESim} are the sum and the rest is set to {\vr -1}',
            )
        )
        scale = .35
        doc.frame('matching simulation and extraction', 
            doc.figure('damicviewerbin1.png', scale=scale, center=True)
        )
        doc.frame('matching simulation and extraction', 
            doc.figure('damicviewerEff.png', scale=scale, center=True)
        )
        doc.frame('matching simulation and extraction', 
            doc.figure('damicviewer_noise.png', scale=scale, center=True)
        )
        doc.frame('matching simulation and extraction', 
            doc.figure('damicviewer_double.png', scale=scale, center=True)
        )
        
        doc.frame('usign the code',
            doc.itemize(
                'the code is available at: https://github.com/PhMota/CONNIEtools',
                'tentative usage documentation: https://github.com/PhMota/CONNIEtools/wiki',
                'relevant commands for today\n'+doc.code(
                    './catalog simulation newcatalog -t 60 30 -b 3 -N 1000 --number-of-images 5 --charge-range 7.25 --image-modes 1 5 -rn 12 -dc 0.1', 'sh'
                )+B+doc.code(
                    "./catalog scatter newcatalog/newcatalog_bin1t60.root --branch-selections 'Esim' '(E1-ESim)/ESim' 'zSim' 'idSim==1' --x-range 0 1000 --y-range -.25 .25 --errorbar --png --fit 'fit' 'mle'", 'sh'
                )
            )
        )
        
        figscale = .37
        cat = Catalog(['60.','30.'], 3, 1000, 5, ['1','200'], ['1','5'], rn=0., dc=0. )
        selection = 'idSim==1 and runID==0'
        _xrange = '0 1000'
        _yrange = '-.25 .25'
        redo = False
        figure1x1 = doc.figure(
                        cat.scatter( 1, 60., 'ESim', '(E1-ESim)/ESim', z = 'zSim', selection=selection, _xrange=_xrange, _yrange=_yrange, fit=None, redo=redo ),
                        scale=figscale
                    )
        figure1x5 = doc.figure(
                        cat.scatter( 5, 60., 'ESim', '(E1-ESim)/ESim', z = 'zSim', selection=selection, _xrange=_xrange, _yrange=_yrange, fit=None, redo=redo ),
                        scale=figscale
                    )
        table = [[ '1x1', figure1x1 ],
                [ '1x5', figure1x5 ],
                ]
        table = array(table).T.tolist()
        item = r'energy reconstruction'
        doc.frame('energy reconstruction for noise {} and DC {}'.format(0, 0),
            doc.itemize(item),
            doc.tabular( table, align='c' )
        )

        redo = False
        figure1x1 = doc.figure(
                        cat.scatter( 1, 30., 'ESim', '(E1-ESim)/ESim', z = 'zSim', selection=selection, _xrange=_xrange, _yrange=_yrange, fit=None, redo=redo ),
                        scale=figscale
                    )
        figure1x5 = doc.figure(
                        cat.scatter( 5, 30., 'ESim', '(E1-ESim)/ESim', z = 'zSim', selection=selection, _xrange=_xrange, _yrange=_yrange, fit=None, redo=redo ),
                        scale=figscale
                    )
        table = [[ '1x1', figure1x1 ],
                [ '1x5', figure1x5 ],
                ]
        table = array(table).T.tolist()
        item = r'energy reconstruction'
        doc.frame('energy reconstruction for noise {} and DC {}'.format(0, 0),
            doc.itemize(item),
            doc.tabular( table, align='c' )
        )

        ############################################################################
        doc.frame('',r'{\center\Large fixing energy reconstruction}')
        
        doc.frame('fit',
            doc.itemize(
                r'''
                {\vr Fit}: minimize directly the scaled pixel-Gaussian
                $$
                f(E, \{i\}) = E {\mathcal I}_i
                $$
                $$
                {\mathcal I}_i = \frac{1}{4}[ {\rm erf}( \frac{ \xi - \mu_x }{ \sqrt{2}\sigma_x } ) ]_x^{x+1} 
                [ {\rm erf}( \frac{ \xi - \mu_y }{ \sqrt{2}\sigma_y } ) ]_y^{y+1}
                $$                
                by the cost function
                $$
                {\rm min} \sum_i ( \frac{e_i - f(E,\{i\}) }{\sigma_E} )^2
                $$
                ''',
                r'''
                {\vr Fit2}: improvement by correcting for the missing probability
                $$
                f(E, \{i\}) = E {\mathcal I}_i \sum_i {\mathcal I}_i
                $$
                '''
            )
        )

        xsigma_xlabel = 'sigmaSim'
        xsigma_ylabel = '(sqrt(xVar1)-sigmaSim)/sigmaSim'

        ysigma_xlabel = 'sigmaSim'
        ysigma_ylabel = '(sqrt(yVar1)-sigmaSim)/sigmaSim'

        E_xlabel = 'ESim'
        E_ylabel = '(E1-ESim)/ESim'

        E2_xlabel = 'zSim'
        E2_ylabel = '(E1-ESim)/ESim'
        
        combinations = [[0.,0.],
                        [12., 0.],
                        [12., 0.1],
                        ]
        for index, (xlabel, ylabel) in enumerate([[E_xlabel, E_ylabel], [E2_xlabel, E2_ylabel], 
                                                  #[ysigma_xlabel, ysigma_ylabel], [xsigma_xlabel, xsigma_ylabel] 
                                                  ]):
            for rn, dc in combinations:
                cat = Catalog(['60.','30.'], 3, 1000, 5, ['1','200'], ['1','5'], rn=rn, dc=dc )
                if xlabel == E_xlabel: 
                    item = r'simulated and reconstructed energy relative error'
                    item1 = r'{\color{blue} \tt E1}, {\color{orange} \tt Efit}, {\color{green} \tt Efit2}'
                    _xrange = '0 1000'
                    _yrange = '-.5 .5'
                    redo = False
                elif xlabel == E2_xlabel: 
                    item = r'simulated and reconstructed energy relative error'
                    item1 = r'{\color{blue} \tt E1}, {\color{orange} \tt Efit}, {\color{green} \tt Efit2}'
                    _xrange = '0 700'
                    _yrange = '-.5 .5'
                    redo = False
                elif 'xVar' in ylabel:
                    item = r'simulated and reconstructed x-width relative error {\color{red} swapped axes}'
                    item1 = r'{\color{blue} \tt xVar1}, {\color{orange} \tt xSigmafit}, {\color{green} \tt xSigmafit2}'
                    _xrange = '0 1.2'
                    _yrange = '-1 1'
                    redo = False
                elif 'yVar' in ylabel:
                    item = r'simulated and reconstructed y-width relative error {\color{red} swapped axes}'
                    item1 = r'{\color{blue} \tt yVar1}, {\color{orange} \tt ySigmafit}, {\color{green} \tt ySigmafit2}'
                    _xrange = '0 1.2'
                    _yrange = '-1 1'
                    redo = False
                else:
                    print('combination not predicted', xlabel, ylabel)

                figure1x1 = doc.figure(
                                cat.scatter( 1, 60., xlabel, ylabel, selection=selection, _xrange=_xrange, _yrange=_yrange, fit=['fit', 'fitEsigma'], redo=redo ),
                                scale=figscale
                            )
                figure1x5 = doc.figure(
                                cat.scatter( 5, 60., xlabel, ylabel, selection=selection, _xrange=_xrange, _yrange=_yrange, fit=['fit', 'fitEsigma'], redo=redo ),
                                scale=figscale
                            )
                
                table = [[ '1x1', figure1x1 ],
                        [ '1x5', figure1x5 ],
                        ]
                table = array(table).T.tolist()

                doc.frame('energy reconstruction for noise {} and DC {}'.format(rn, dc),
                    doc.itemize(item,item1),
                    doc.tabular( table, align='c' )
                )

        #cat = Catalog(['60.','30.'], 3, 1000, 5, ['1','200'], ['1','5'], rn=12., dc=0.1 )
        
        #xlabel = 'n0'
        #ylabel = '(E1-ESim)/ESim'
        
        #item = r'simulated and reconstructed energy relative error'
        #item1 = r'{\color{blue} \tt E1}, {\color{orange} \tt Efit}, {\color{green} \tt Efit2}'
        #_xrange = '0 20'
        #_yrange = '-.5 .5'
        #redo = False

        #figure1x1 = doc.figure(
                        #cat.scatter( 1, 60., xlabel, ylabel, selection=selection, _xrange=_xrange, _yrange=_yrange, fit=['fit', 'fitEsigma'], redo=redo ),
                        #scale=figscale
                    #)
        #figure1x5 = doc.figure(
                        #cat.scatter( 5, 60., xlabel, ylabel, selection=selection, _xrange=_xrange, _yrange=_yrange, fit=['fit', 'fitEsigma'], redo=redo ),
                        #scale=figscale
                    #)
        

        ############################################################################
        doc.frame('',r'{\center\Large size like comparion}')
        doc.frame('compare to size like', 
            doc.itemize(
                r'''
                size like approach (Miguel's thesis)\\
                the probability of a charge falling in a given pixel $i$
                $$
                {\mathcal B}(Q, q_i, \{i\}) = \binom{Q}{q_i} {\mathcal I}_i^q( 1 - {\mathcal I}_i )^{Q-q_i}
                $$
                ${\mathcal I}_i$ is integral of the 2d-Gaussian distribution in the pixel
                $$
                {\mathcal I}_i = \frac{1}{4}[ {\rm erf}( \frac{ \xi - \mu_x }{ \sqrt{2}\sigma_x } ) ]_x^{x+1} 
                [ {\rm erf}( \frac{ \xi - \mu_y }{ \sqrt{2}\sigma_y } ) ]_y^{y+1}
                $$
                '''
            )
        )

        doc.frame('note',
            doc.itemize(
                r'''
                further convolution with the energy gaussian distribution
                $$
                {\rm max} \sum_i {\rm log}[ \sum_{\xi=0}^Q {\mathcal B}(Q, \xi, i) {\mathcal G}_e( e_i-\color{red}{g} \xi ) ]
                $$
                high computational cost
                ''',
                'requires charge gain information $q_i=e_i/g$',
                r'takes 2m/1k hits compared to 10s/1k for the fit {\color{red} $\sim \times 12$}'
            )
        )
        
        xsigma_xlabel = 'sigmaSim'
        xsigma_ylabel = '(sqrt(xVar1)-sigmaSim)/sigmaSim'

        ysigma_xlabel = 'sigmaSim'
        ysigma_ylabel = '(sqrt(yVar1)-sigmaSim)/sigmaSim'

        E_xlabel = 'ESim'
        E_ylabel = '(E1-ESim)/ESim'

        combinations = [[0.,0.],
                        [12., 0.],
                        [12., 0.1],
                        ]
        for rn, dc in combinations:
            for index, (xlabel, ylabel) in enumerate([[E_xlabel, E_ylabel], [ysigma_xlabel, ysigma_ylabel], 
                                                      #[xsigma_xlabel, xsigma_ylabel] 
                                                      ]):
                cat = Catalog(['60.','30.'], 3, 1000, 5, ['1','200'], ['1','5'], rn=rn, dc=dc )
                if 'ESim' in xlabel: 
                    item = r'simulated and reconstructed energy relative error'
                    item1 = r'{\color{blue} \tt E1}, {\color{orange} \tt Efit}, {\color{green} \tt Emle}'
                    _xrange = '0 1000'
                    _yrange = '-.5 .5'
                    redo = False
                elif 'xVar' in ylabel:
                    item = r'simulated and reconstructed x-width relative error {\color{red} swapped axes}'
                    item1 = r'{\color{blue} \tt xVar1}, {\color{orange} \tt xSigmafit}, {\color{green} \tt xSigmamle}'
                    _xrange = '0 1.2'
                    _yrange = '-1 1'
                    redo = False
                elif 'yVar' in ylabel:
                    item = r'simulated and reconstructed y-width relative error {\color{red} swapped axes}'
                    item1 = r'{\color{blue} \tt yVar1}, {\color{orange} \tt ySigmafit}, {\color{green} \tt ySigmamle}'
                    _xrange = '0 1.2'
                    _yrange = '-1 1'
                    redo = False
                else:
                    print('combination not predicted', xlabel, ylabel)

                figure1x1 = doc.figure(
                                cat.scatter( 1, 60., xlabel, ylabel, selection=selection, _xrange=_xrange, _yrange=_yrange, fit=['fit', 'mle'], redo=redo ),
                                scale=figscale
                            )
                figure1x5 = doc.figure(
                                cat.scatter( 5, 60., xlabel, ylabel, selection=selection, _xrange=_xrange, _yrange=_yrange, fit=['fit', 'mle'], redo=redo ),
                                scale=figscale
                            )
                
                table = [[ '1x1', figure1x1 ],
                        [ '1x5', figure1x5 ],
                        ]
                table = array(table).T.tolist()

                doc.frame('sigma efficiency for noise {} and DC {}'.format(rn, dc),
                    doc.itemize(item,item1),
                    doc.tabular( table, align='c' )
                )

        doc.frame('summary',
            doc.itemize(
                'new tool for matching catalogs with simulation',
                'new tool for energy and sigma estimation',
                'good agreement with noise-like for sigma, better performance for energy',
                'needs to remove DC contribution from the energy (MB subtraction as explained by Juan)',
                ),
        )

