from __future__ import print_function
from Beamer import *
import os
import subprocess
import sys

folder = 'ImagePresentation'

if not os.path.exists(folder):
    os.mkdir(folder)

with openBeamer( folder ) as doc:
    
    name = 'median'
    fname = r'/share/storage2/connie/data/runs/*/runID_*_03326_*_p*.fits.fz'
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
    
    #scale = .6
    #doc.frame('raw parts',
              #'raw parts from the fits file',
              #doc.center( doc.figure('parto3p1.pdf', scale=scale, func=func) ),
              #doc.center( doc.figure('parto3p2.pdf', scale=scale, func=func) ),
              #doc.center( doc.figure('parto3p3.pdf', scale=scale, func=func) ),
              #doc.center( doc.figure('parto3p4.pdf', scale=scale, func=func) ),
              #)
    
    scale = .7
    doc.frame('sections of the raw image', 
              doc.center(
              doc.figure( name+'/median_o3_data.pdf', scale=scale), 
              doc.figure( name+'/median_o3_bias.pdf', scale=scale ), 
              doc.figure( name+'/median_o3_vbias.pdf', scale=.7*scale ), 
              doc.figure( name+'/median_o3_dbias.pdf', scale=.1*scale ),
              ))
    
    scale = .3
    doc.frame('projections of the raw image', 
              doc.center(
              doc.figure( name+'/median_o3_biasProj.pdf', scale=scale ), 
              doc.figure( name+'/median_o3_vbiasProj.pdf', scale=scale ),
              ))

    scale = .3
    doc.frame('spectra of the raw image', 
              'distributions are crowded with outliers',
              doc.center(
              doc.figure( name+'/median_o3_data_spectrum.pdf', scale=scale ), 
              doc.figure( name+'/median_o3_bias_spectrum.pdf', scale=scale ), 
              doc.figure( name+'/median_o3_vbias_spectrum.pdf', scale=scale ), 
              doc.figure( name+'/median_o3_dbias_spectrum.pdf', scale=scale ), 
              ))
    
    name = 'mean'
    cmd = 'python Image.py analyse {folder}/{name} "{fname}" --ohdu 3 --params-mode mean --plot-sections --plot-spectrum'\
        .format(fname=fname, folder=folder, name=name)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )

    doc.set_func(func)
    
    doc.frame('Mean outliers analysis with line and row corrections', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file=name+'/mean_params.csv', fontsize=5, spacing=6 ) ), 
              )
        
    scale = .7
    doc.frame('sides with line and row corrections', 
              doc.center(
              doc.figure(name+'/mean_o3_data.pdf', scale=scale ), 
              doc.figure(name+'/mean_o3_bias.pdf', scale=scale ), 
              doc.figure(name+'/mean_o3_vbias.pdf', scale=.7*scale ), 
              doc.figure(name+'/mean_o3_dbias.pdf', scale=.1*scale ),
              ))
    
    scale = .3
    doc.frame('projections with line and row corrections', 
              doc.center(
              doc.figure(name+'/mean_o3_biasProj.pdf', scale=scale ), 
              doc.figure(name+'/mean_o3_vbiasProj.pdf', scale=scale ),
              ))
    
    scale = .3
    doc.frame('spectra with line and row corrections', 
              doc.center(
              doc.figure(name+'/mean_o3_data_spectrum.pdf', scale=scale ), 
              doc.figure(name+'/mean_o3_bias_spectrum.pdf', scale=scale ), 
              doc.figure(name+'/mean_o3_vbias_spectrum.pdf', scale=scale ), 
              doc.figure(name+'/mean_o3_dbias_spectrum.pdf', scale=scale ), 
              ))


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
              doc.center( doc.table( file=name+'/{}_params.csv'.format(pre), fontsize=5, spacing=6, func=func ) ), 
              )
        
    scale = .7
    doc.frame('removed above 40ADU with border 3', doc.center(
              doc.figure(name+'/{}_o3_data.pdf'.format(pre), scale=scale, func=func ), 
              doc.figure(name+'/{}_o3_bias.pdf'.format(pre), scale=scale, func=func ), 
              doc.figure(name+'/{}_o3_vbias.pdf'.format(pre), scale=.7*scale, func=func ), 
              doc.figure(name+'/{}_o3_dbias.pdf'.format(pre), scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('removed above 40ADU with border 3 spectra', doc.center(
              doc.figure(name+'/{}_o3_data_spectrum.pdf'.format(pre), scale=scale, func=func), 
              doc.figure(name+'/{}_o3_bias_spectrum.pdf'.format(pre), scale=scale, func=func), 
              doc.figure(name+'/{}_o3_vbias_spectrum.pdf'.format(pre), scale=scale, func=func), 
              doc.figure(name+'/{}_o3_dbias_spectrum.pdf'.format(pre), scale=scale, func=func), 
              ),
                'maybe too agressive',
    )
    
    for threshold, border in [[60.0,3.0], [80.0,2.0], [150.0,0.0]]:
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
                doc.center( doc.table( file=name+'/{}_params.csv'.format(pre), fontsize=5, spacing=6, func=func ) ), 
                )
            
        scale = .3
        doc.frame('removed above {}ADU with border {} spectra'.format(threshold,border), 
                doc.center(
                doc.figure(name+'/{}_o3_dbias_spectrum.pdf'.format(pre), scale=scale, func=func), 
                doc.figure(name+'/{}_o3_bias_spectrum.pdf'.format(pre), scale=scale, func=func), 
                doc.figure(name+'/{}_o3_vbias_spectrum.pdf'.format(pre), scale=scale, func=func), 
                doc.figure(name+'/{}_o3_data_spectrum.pdf'.format(pre), scale=scale, func=func), 
                ),
        )
    
