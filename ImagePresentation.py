from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys

folder = 'ImagePresentation'

if not os.path.exists(folder):
    os.mkdir(folder)

def projections_slide( name, cmd ):
    scale = .3
    contents = ['projections of the raw image', doc.code( cmd, 'Bash') ]
    contents += [ doc.figure( name+'/median_o3_{}.pdf'.format(key), scale=scale ) for key in ['dataProj1', 'dataProj0', 'biasProj', 'vbiasProj' ] ]
    doc.frame( *contents )

with Timer('presentation'):
    with openBeamer( folder ) as doc:

        name = 'median'
        fname = r'/share/storage2/connie/data/runs/*/runID_*_03326_*_p*.fits.fz'
        cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --plot-sections --plot-spectrum'\
            .format( name=name, fname=fname, folder=folder)
        print( cmd )
        func = lambda: subprocess.call( [cmd], shell=True )
        
        doc.set_func(func)
        doc.frame('tool for image analysis', 
                'sample call',
                doc.code( cmd, 'Bash'),
                'other functionalities',
                doc.itemize( 'read header', 'simulate and get params', 'extract hits (next week: comparison with offical extraction)' ),
                'run locally',
                doc.code( '\n'.join(( 'git init', 'git pull https://github.com/PhMota/CONNIEtools')), 'Bash' )
                )
            
        scale = .6
        optsSection = [['data',1], ['bias',1], ['vbias',.7], ['dbias',.1] ]
        figs = [ doc.figure( name+'/median_o3_{}.pdf'.format(key), scale=x*scale) for key, x in optsSection ]
        doc.frame('1x1 raw sections', 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figs ),
                    'vertical modulation is clearly visible',
                    widths = [.7,-1]
                )
                )
        
        scale = .25
        optsProj = ['dataProj1', 'dataProj0', 'biasProj', 'vbiasProj']
        figs = [ doc.figure( name+'/median_o3_{}.pdf'.format(key), scale=scale ) for key in optsProj ]
        doc.frame('projections of the raw image', 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figs ),
                    'vertical modulation, horizontal modulation, hot columns, (global median subtraction)\n\n\n\
                    {\small $\sigma = $mean(MADs)\n\
                    $g\lambda = $mean(med-med)}',
                    widths = [.6,-1]
                ))

        scale = .3
        optsSpectrum = ['data', 'bias', 'vbias', 'dbias']
        figs = [ doc.figure( name+'/median_o3_{}_spectrum.pdf'.format(key), scale=scale ) for key in optsSpectrum ]
        doc.frame('spectra of the raw image', 
                doc.code( cmd, 'Bash'),
                'distributions are crowded with outliers',
                doc.center( *figs ))
        
        name = 'mean'
        cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --plot-sections --plot-spectrum'\
            .format(fname=fname, folder=folder, name=name)
        print( cmd )
        func = lambda: subprocess.call( [cmd], shell=True )

        doc.set_func(func)
        
        doc.frame('Mean outliers analysis with line and col corrections', 
                'bash command',
                doc.code( cmd, 'Bash'),
                'parameters',
                doc.center( doc.table( file=name+'/mean_params.csv', fontsize=5, spacing=6, divide=2 ) ), 
                )

        scale = .6
        figs = [ doc.figure( name+'/mean_o3_{}.pdf'.format(key), scale=x*scale) for key, x in optsSection ]
        doc.frame('sides with line and col corrections', 
                doc.code( cmd, 'Bash'),
                doc.center( *figs ))
                #doc.figure(name+'/mean_o3_data.pdf', scale=scale ), 
                #doc.figure(name+'/mean_o3_bias.pdf', scale=scale ), 
                #doc.figure(name+'/mean_o3_vbias.pdf', scale=.7*scale ), 
                #doc.figure(name+'/mean_o3_dbias.pdf', scale=.1*scale ),
                #))
        
        scale = .3
        figs = [ doc.figure( name+'/mean_o3_{}.pdf'.format(key), scale=scale ) for key in optsProj ]        
        doc.frame('projections with line and col corrections', 
                doc.code( cmd, 'Bash'),
                doc.center( *figs ))
                #doc.figure(name+'/mean_o3_biasProj.pdf', scale=scale ), 
                #doc.figure(name+'/mean_o3_vbiasProj.pdf', scale=scale ),
                #))
        
        scale = .3
        figs = [ doc.figure( name+'/mean_o3_{}_spectrum.pdf'.format(key), scale=scale ) for key in optsSpectrum ]
        doc.frame('spectra with line and col corrections', 
                doc.code( cmd, 'Bash'),
                doc.center( *figs ) )
                #doc.figure(name+'/mean_o3_data_spectrum.pdf', scale=scale ), 
                #doc.figure(name+'/mean_o3_bias_spectrum.pdf', scale=scale ), 
                #doc.figure(name+'/mean_o3_vbias_spectrum.pdf', scale=scale ), 
                #doc.figure(name+'/mean_o3_dbias_spectrum.pdf', scale=scale ), 
                #))


        name = 'mean'
        cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --remove-hits 40 3 --params-mode mean --plot-sections --plot-spectrum'\
            .format(fname=fname, folder=folder, name=name)
        print( cmd )
        func = lambda: subprocess.call( [cmd], shell=True )
        doc.set_func(func)
        
        pre = 'mean_e40.0b3.0'
        
        doc.frame('removed above 40ADU with border 3 params', 
                'bash command',
                doc.code( cmd, 'Bash'),
                'parameters',
                doc.center( doc.table( file=name+'/{}_params.csv'.format(pre), fontsize=5, spacing=6, divide=2 ) ), 
                )
            
        scale = .7
        doc.frame('removed above 40ADU with border 3', doc.center(
                doc.code( cmd, 'Bash'),
                doc.figure(name+'/{}_o3_data.pdf'.format(pre), scale=scale ), 
                doc.figure(name+'/{}_o3_bias.pdf'.format(pre), scale=scale ), 
                doc.figure(name+'/{}_o3_vbias.pdf'.format(pre), scale=.7*scale ), 
                doc.figure(name+'/{}_o3_dbias.pdf'.format(pre), scale=.1*scale ),
                ))
        
        scale = .3
        doc.frame('removed above 40ADU with border 3 spectra', doc.center(
                doc.code( cmd, 'Bash'),
                doc.figure(name+'/{}_o3_data_spectrum.pdf'.format(pre), scale=scale ), 
                doc.figure(name+'/{}_o3_bias_spectrum.pdf'.format(pre), scale=scale ), 
                doc.figure(name+'/{}_o3_vbias_spectrum.pdf'.format(pre), scale=scale ), 
                doc.figure(name+'/{}_o3_dbias_spectrum.pdf'.format(pre), scale=scale ), 
                ),
                    'maybe too agressive',
        )
        
        for threshold, border in [[60.0,3.0], [80.0,2.0], [150.0,0.0], [5000.0,0.0]]:
            name = 'mean'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-spectrum'\
                .format( fname=fname, folder=folder, name=name, threshold=threshold, border=border )
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            doc.set_func(func)
            pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
            
            doc.frame('removed above {}ADU with border {} spectra'.format(threshold, border),
                    'bash command',
                    doc.code( cmd, 'Bash'),
                    'parameters',
                    doc.center( doc.table( file=name+'/{}_params.csv'.format(pre), fontsize=5, spacing=6, divide=2 ) ), 
                    )
                
            scale = .3
            doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure(name+'/{}_o3_data_spectrum.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_bias_spectrum.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_vbias_spectrum.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_dbias_spectrum.pdf'.format(pre), scale=scale), 
                    ),
            )
        
        with Timer('convolution'):
            threshold = 100.0
            border = 3.0
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-convolution-spectrum'\
                .format( fname=fname, folder=folder, name=name, threshold=threshold, border=border )
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            doc.set_func(func)
            
            pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
            scale = .3
            doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure(name+'/{}_o3_data_convolution5.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_bias_convolution5.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_vbias_convolution5.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_dbias_convolution5.pdf'.format(pre), scale=scale), 
                    ),
            )
            doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure(name+'/{}_o3_data_convolution10.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_bias_convolution10.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_vbias_convolution10.pdf'.format(pre), scale=scale), 
                    doc.figure(name+'/{}_o3_dbias_convolution10.pdf'.format(pre), scale=scale), 
                    ),
            )

        with Timer('blocks'):
            threshold = 100.0
            border = 3.0
            for _name_, _func_ in [['blockmean', 'np.nanmean(x)'], ['blockstd', '(lambda y: np.nan if y==0 else y)(np.nanstd(x))']]:
                cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-block-spectrum --block-function "{_func_}"'\
                    .format( fname=fname, folder=folder, name=_name_, threshold=threshold, border=border, _func_=_func_ )
                print( cmd )
                func = lambda: subprocess.call( [cmd], shell=True )
                doc.set_func(func)
                
                pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
                scale = .3
                doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                        doc.center(
                        doc.figure(_name_ + '/{}_o3_data_block5.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_bias_block5.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_vbias_block5.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_dbias_block5.pdf'.format(pre), scale=scale), 
                        ),
                )
                doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                        doc.center(
                        doc.figure(_name_ + '/{}_o3_data_block10.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_bias_block10.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_vbias_block10.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_dbias_block10.pdf'.format(pre), scale=scale), 
                        ),
                )

        with Timer('evolution'):
            threshold = 100.0
            border = 3.0
            _name_ = 'dcevo'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-sections'\
                .format( fname=fname, folder=folder, name=_name_, threshold=threshold, border=border )
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            doc.set_func(func)
            
            pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
            scale = .3
            doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure(_name_ + '/{}_o3_dataProj1.pdf'.format(pre), scale=scale), 
                    doc.figure(_name_ + '/{}_o3_dataProj0.pdf'.format(pre), scale=scale), 
                    doc.figure(_name_ + '/{}_o3_biasProj.pdf'.format(pre), scale=scale), 
                    doc.figure(_name_ + '/{}_o3_vbiasProj.pdf'.format(pre), scale=scale), 
                    ),
            )
                    


        fname = r'/share/storage2/connie/data/runs/*/runID_*_12000_*_p*.fits.fz'
        with Timer('evolution'):
            threshold = 100.0
            border = 3.0
            _name_ = 'dcevo1x5'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-sections --plot-spectrum'\
                .format( fname=fname, folder=folder, name=_name_, threshold=threshold, border=border )
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            doc.set_func(func)
            
            pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
            scale = .3
            doc.frame('removed above {}ADU with border {} spectra 1x5'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure(_name_ + '/{}_o3_dataProj1.pdf'.format(pre), scale=scale), 
                    doc.figure(_name_ + '/{}_o3_dataProj0.pdf'.format(pre), scale=scale), 
                    doc.figure(_name_ + '/{}_o3_biasProj.pdf'.format(pre), scale=scale), 
                    doc.figure(_name_ + '/{}_o3_vbiasProj.pdf'.format(pre), scale=scale), 
                    ),
            )
            scale = .3
            doc.frame('spectra with line and col corrections 1x5', 
                doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure(_name_+'/{}_o3_data_spectrum.pdf'.format(pre), scale=scale ), 
                    doc.figure(_name_+'/{}_o3_bias_spectrum.pdf'.format(pre), scale=scale ), 
                    doc.figure(_name_+'/{}_o3_vbias_spectrum.pdf'.format(pre), scale=scale ), 
                    doc.figure(_name_+'/{}_o3_dbias_spectrum.pdf'.format(pre), scale=scale ), 
                    ))

        with Timer('blocks'):
            threshold = 100.0
            border = 3.0
            for _name_, _func_ in [['blockmean1x5', 'np.nanmean(x)'], ['blockstd1x5', '(lambda y: np.nan if y==0 else y)(np.nanstd(x))']]:
                cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-block-spectrum --block-function "{_func_}"'\
                    .format( fname=fname, folder=folder, name=_name_, threshold=threshold, border=border, _func_=_func_ )
                print( cmd )
                func = lambda: subprocess.call( [cmd], shell=True )
                doc.set_func(func)
                
                pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
                scale = .3
                doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                        doc.center(
                        doc.figure(_name_ + '/{}_o3_data_block5.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_bias_block5.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_vbias_block5.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_dbias_block5.pdf'.format(pre), scale=scale), 
                        ),
                )
                doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                        doc.center(
                        doc.figure(_name_ + '/{}_o3_data_block10.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_bias_block10.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_vbias_block10.pdf'.format(pre), scale=scale), 
                        doc.figure(_name_ + '/{}_o3_dbias_block10.pdf'.format(pre), scale=scale), 
                        ),
                )

        with Timer():
            name = 'median1x5'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --plot-sections --plot-spectrum'\
                .format( name=name, fname=fname, folder=folder)
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            
            doc.set_func(func)
            doc.frame('Monitor Viewer Calculation', 
                    'bash command',
                    doc.code( cmd, 'Bash'),
                    'parameters calculated directly on the raw file',
                    doc.center( doc.table( file= name+'/median_params.csv', fontsize=5, spacing=6, divide=2 ) ),
                    )
                
            scale = .7
            doc.frame('sections of the raw image', 
                    doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure( name+'/median_o3_data.pdf', scale=scale), 
                    doc.figure( name+'/median_o3_bias.pdf', scale=scale ), 
                    doc.figure( name+'/median_o3_vbias.pdf', scale=.7*scale ), 
                    doc.figure( name+'/median_o3_dbias.pdf', scale=.1*scale ),
                    ))
            
            scale = .3
            doc.frame('projections of the raw image', 
                    doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure( name+'/median_o3_dataProj1.pdf', scale=scale ), 
                    doc.figure( name+'/median_o3_dataProj0.pdf', scale=scale ),
                    doc.figure( name+'/median_o3_biasProj.pdf', scale=scale ), 
                    doc.figure( name+'/median_o3_vbiasProj.pdf', scale=scale ),
                    ))

            scale = .3
            doc.frame('spectra of the raw image', 
                    'distributions are crowded with outliers',
                    doc.code( cmd, 'Bash'),
                    doc.center(
                    doc.figure( name+'/median_o3_data_spectrum.pdf', scale=scale ), 
                    doc.figure( name+'/median_o3_bias_spectrum.pdf', scale=scale ), 
                    doc.figure( name+'/median_o3_vbias_spectrum.pdf', scale=scale ), 
                    doc.figure( name+'/median_o3_dbias_spectrum.pdf', scale=scale ), 
                    ))
