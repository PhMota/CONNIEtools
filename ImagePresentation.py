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
            
        scaleSection = .6
        optsSection = [['data',1], ['bias',1], ['vbias',.7], ['dbias',.1] ]
        figsSection = lambda name,pre,optsSection=optsSection: [ doc.figure( name+'/{}_o3_{}.pdf'.format(pre,key), scale=x*scaleSection) for key, x in optsSection ]
        doc.frame('1x1 raw sections', 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figsSection('median', 'median') ),
                    'vertical modulation is clearly visible',
                    widths = [.7,-1]
                )
                )
        
        scaleProj = .25
        optsProj = ['dataProj1', 'dataProj0', 'biasProj', 'vbiasProj']
        figsProj = lambda name, pre: [ doc.figure( name+'/{}_o3_{}.pdf'.format(pre,key), scale=scaleProj ) for key in optsProj ]
        doc.frame('projections of the raw image', 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figsProj('median', 'median') ),
                    'vertical modulation, horizontal modulation, hot columns, (global median subtraction)\n\n\n\
                    {\small $\sigma = $mean(MADs)\n\
                    $g\lambda = $mean(med-med)}',
                    widths = [.6,-1]
                ))

        scaleSpectrum = .3
        optsSpectrum = ['data', 'bias', 'vbias', 'dbias']
        figsSpectrum = lambda name, pre, kind='_spectrum': [ doc.figure( name+'/{}_o3_{}{}.pdf'.format(pre,key,kind), scale=scaleSpectrum ) for key in optsSpectrum ]
        doc.frame('spectra of the raw image', 
                doc.code( cmd, 'Bash'),
                'distributions are crowded with outliers',
                doc.center( *figsSpectrum('median', 'median') ))

        with Timer():
            name = 'mean'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --plot-sections --plot-spectrum'\
                .format(fname=fname, folder=folder, name=name)
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )

            doc.set_func(func)
            
            #doc.frame('Mean outliers analysis with line and col corrections', 
                    #'bash command',
                    #doc.code( cmd, 'Bash'),
                    #'parameters',
                    #doc.center( doc.table( file=name+'/mean_params.csv', fontsize=5, spacing=6, divide=2 ) ), 
                    #)

            doc.frame('sides with line and col corrections', 
                    doc.code( cmd, 'Bash'),
                    doc.center( *figsSection('mean', 'mean') ))
            
            doc.frame('projections with line and col corrections', 
                    doc.code( cmd, 'Bash'),
                    doc.center( *figsProj('mean','mean') ))
            
            doc.frame('spectra with line and col corrections', 
                    doc.code( cmd, 'Bash'),
                    doc.center( *figsSpectrum('mean','mean') ) )


        with Timer():
            name = 'mean'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --remove-hits 40 3 --params-mode mean --plot-sections --plot-spectrum'\
                .format(fname=fname, folder=folder, name=name)
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            doc.set_func(func)
            
            pre = 'mean_e40.0b3.0'
            
            figs = [ doc.figure( name+'/{}_o3_{}.pdf'.format(pre,key), scale=x*scaleSection) for key, x in optsSection ]
            doc.frame('removed above 40ADU with border 3',
                    doc.code( cmd, 'Bash'),
                    doc.center( *figsSection('mean', pre) )
            )
            
            scale = .3
            doc.frame('removed above 40ADU with border 3 spectra',
                    doc.code( cmd, 'Bash'),
                    doc.center( *figsSpectrum('mean', pre) ), 
                    'maybe too agressive',
            )
                    #thr  border noise dc
        thislist = [[60.0,   3.0, 14.4, 2.41, 15.39**2 - 14.5**2 ], 
                    [80.0,   2.0, 14.4, 2.42, 15.39**2 - 14.5**2 ],
                    [150.0,  0.0, 14.4, 2.72, 15.6**2 - 14.5**2 ], 
                    [5000.0, 0.0, 14.4, 3.42, 16.3**2 - 14.5**2 ]]
        for threshold, border, noise, dc, dc2 in thislist:
            name = 'mean'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-spectrum'\
                .format( fname=fname, folder=folder, name=name, threshold=threshold, border=border )
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            doc.set_func(func)
            pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
                            
            scale = .3
            doc.frame(
                '$E<{:.0f}$ADU (+{:.0f} border) spectra'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figsSpectrum('mean', pre) ),
                    'estimations\n\n\
                    $$\sigma={noise}$$\n\
                    $$g\lambda={dc}$$\n\
                    $$g^2\lambda={dc2}$$\n\
                    $$g={g:.2f}$$\n\
                    $$\lambda={lamb:.2f}$$'.format(noise=noise, dc=dc, dc2=dc2, g=dc2/dc, lamb=dc**2/dc2),
                    widths = [.7,-1]
                )
            )
        
        #with Timer('convolution'):
            #threshold = 100.0
            #border = 3.0
            #cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-convolution-spectrum'\
                #.format( fname=fname, folder=folder, name=name, threshold=threshold, border=border )
            #print( cmd )
            #func = lambda: subprocess.call( [cmd], shell=True )
            #doc.set_func(func)
            
            #pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
            #thislist = [[5, 64.6, 52, 77.1],
                        #[10, 127.6, 223.3, 165.7]]
            #for n, noise, dc, dc2 in thislist:
                #dc2 = (dc2**2 - noise**2)
                #doc.frame('convolution sum {D}x{D} $E<{:.0f}$(+{:.0f})'.format(threshold,border,D=n), 
                    #doc.code( cmd, 'Bash'),
                    #doc.column(
                        #doc.center( *figsSpectrum('mean', pre, kind='_convolution%d'%n) ),
                    #'estimations\n\n\
                    #$$\sigma={noise}$$\n\
                    #$$g\lambda={dc}$$\n\
                    #$$g^2\lambda={dc2}$$\n\
                    #$$g={g:.2f}$$\n\
                    #$$\lambda={lamb:.2f}$$'.format(noise=noise, dc=dc, dc2=dc2, g=dc2/dc, lamb=dc**2/dc2),
                        #widths=[.7,-1] )
                #)

        with Timer('blocks'):
            threshold = 100.0
            border = 3.0

            gl5 = 2.54
            sig5 = 14
            g2l5 = 14.8**2-sig5**2
            
            gl10 = 2.55
            sig10 = 14.3
            g2l10 = 15.2**2-sig10**2

            blocks = [
                ['blockmean', 'np.nanmean(x)', 'mean', 5, '$$g\lambda={}$$'.format(gl5)], 
                ['blockstd', '(lambda y: np.nan if y==0 else y)(np.nanstd(x))', 'std', 5, 
                 '$$\sigma={}$$\n$$g^2\lambda={:.0f}$$\n$$g={:.2f}$$\n$$\lambda={:.2f}$$'.format(sig5, g2l5, g2l5/gl5, gl5**2/g2l5 ) ],
                ['blockmean', 'np.nanmean(x)', 'mean', 10,  '$$g\lambda={}$$'.format(gl10)],
                ['blockstd', '(lambda y: np.nan if y==0 else y)(np.nanstd(x))', 'std', 10, 
                 '$$\sigma={}$$\n$$g^2\lambda={:.0f}$$\n$$g={:.2f}$$\n$$\lambda={:.2f}$$'.format(sig10, g2l10, g2l10/gl10, gl10**2/g2l10 )],
                    ]
            for _name_, _func_, funclabel, D, text in blocks:
                cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-block-spectrum --block-function "{_func_}"'\
                    .format( fname=fname, folder=folder, name=_name_, threshold=threshold, border=border, _func_=_func_ )
                print( cmd )
                func = lambda: subprocess.call( [cmd], shell=True )
                doc.set_func(func)
                
                pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
                doc.frame('{funclabel} block {D}x{D} $E<{:.0f}$ADU (+{:.0f})'.format(threshold,border,funclabel=funclabel,D=D), 
                    doc.code( cmd, 'Bash'),
                    doc.column(
                        doc.center( *figsSpectrum(_name_, pre, kind='_block%d'%D) ),
                        text,
                        widths=[.7,-1] ))

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
            doc.frame('evolution of DC through time $E<{}$ADU (+{})'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figsProj('dcevo', pre ) ),
                    '',
                    widths = [.7,1]))
                    

        with Timer():
            name = 'median1x5'
            cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --plot-sections --plot-spectrum'\
                .format( name=name, fname=fname, folder=folder)
            print( cmd )
            func = lambda: subprocess.call( [cmd], shell=True )
            
            scale = 1.2
            optsSection = [['data',1*scale], ['bias',.25*scale], ['vbias',1*scale], ['dbias',.2*scale] ]
            doc.frame('1x5 ', 
                    doc.code( cmd, 'Bash'),
                    doc.center( *figsSection(name, 'median', optsSection=optsSection) )
                    )
            
            text = ''
            doc.frame('projections of the raw image', 
                    doc.code( cmd, 'Bash'),
                    doc.column(
                        doc.center( *figsProj(name, 'median') ),
                        text,
                        widths=[.7,-1] ))

            text = ''
            doc.frame('spectra of the raw image', 
                    'distributions are crowded with outliers',
                    doc.code( cmd, 'Bash'),
                    doc.column(
                        doc.center( *figsSpectrum(name, 'median') ),
                        text,
                        widths=[.7,-1] ))

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

            doc.frame('1x5 spectra with line and col corrections 1x5', 
                doc.code( cmd, 'Bash'),
                doc.column(
                        doc.center( *figsSpectrum('dcevo1x5', pre, kind='_spectrum') ),
                        text,
                        widths=[.7,-1]
                    ))
            
            pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
            doc.frame('1x5 evolution of DC $E<{}$ADU (+{})'.format(threshold,border), 
                doc.code( cmd, 'Bash'),
                doc.column(
                    doc.center( *figsProj('dcevo1x5', pre ) ),
                    '',
                    widths = [.7,1]))


        with Timer('blocks'):
            threshold = 100.0
            border = 3.0
            thislist = [
                ['blockmean1x5', 'np.nanmean(x)',5], 
                ['blockstd1x5', '(lambda y: np.nan if y==0 else y)(np.nanstd(x))',5],
                ['blockmean1x5', 'np.nanmean(x)',10], 
                ['blockstd1x5', '(lambda y: np.nan if y==0 else y)(np.nanstd(x))',10],
                ]
            for _name_, _func_, D in thislist:
                cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --remove-hits {threshold} {border} --plot-block-spectrum --block-function "{_func_}"'\
                    .format( fname=fname, folder=folder, name=_name_, threshold=threshold, border=border, _func_=_func_ )
                print( cmd )
                func = lambda: subprocess.call( [cmd], shell=True )
                doc.set_func(func)
                
                pre = 'mean_e{threshold}b{border}'.format(threshold=threshold, border=border)
                scale = .3
                doc.frame('1x5 block {D}x{D} $E<{}$ADU (+{})'.format(threshold,border,D=D), 
                doc.code( cmd, 'Bash'),
                    doc.column(
                        doc.center( *figsSpectrum(_name_, pre, kind='_block%d'%D) ),
                        text,
                        widths=[.7,-1] ))
                


