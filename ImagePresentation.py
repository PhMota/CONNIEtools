from __future__ import print_function
from Beamer import *
import os
import subprocess
import sys

folder = 'ImagePresentation'

if not os.path.exists(folder):
    os.mkdir(folder)

with openBeamer( folder ) as doc:
        
    fname = r'/share/storage2/connie/data/runs/*/runID_*_03326_*_p*.fits.fz'
    cmd = '''\
python Image.py analyse "{fname}" \
--ohdu 3 \
--output-file "{folder}/params.csv" \
--plot-part "{folder}/part*.pdf" \
--plot-sides "{folder}/side*.pdf" \
--plot-spectrum "{folder}/spectrum*.pdf" \
'''.format(fname=fname, folder=folder)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )
    
    doc.frame('Monitor Viewer Calculation', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file='params.csv', fontsize=5, spacing=6, func=func ) ),               
              )
    
    scale = .6
    doc.frame('raw parts',
              'raw parts from the fits file',
              doc.center( doc.figure('parto3p1.pdf', scale=scale, func=func) ),
              doc.center( doc.figure('parto3p2.pdf', scale=scale, func=func) ),
              doc.center( doc.figure('parto3p3.pdf', scale=scale, func=func) ),
              doc.center( doc.figure('parto3p4.pdf', scale=scale, func=func) ),
              )
    
    scale = .7
    doc.frame('sides', 
              'sections of the image',
              doc.center(
              doc.figure('sideo3data.pdf', scale=scale, func=func ), doc.figure('sideo3bias.pdf', scale=scale, func=func ), 
              doc.figure('sideo3vbias.pdf', scale=.7*scale, func=func ), doc.figure('sideo3dbias.pdf', scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('projections', 
              'vertical and horizontal overscan projections along their longest dimension ',
              doc.center(
              doc.figure('sideo3proj1.pdf', scale=scale, func=func), 
              doc.figure('sideo3proj0.pdf', scale=scale, func=func ),
              ))

    scale = .3
    doc.frame('spectra', 
              'distributions are crowded with outliers',
              doc.center(
              doc.figure('spectrumo3data.pdf', scale=scale, func=func), 
              doc.figure('spectrumo3bias.pdf', scale=scale, func=func), 
              doc.figure('spectrumo3vbias.pdf', scale=scale, func=func), 
              doc.figure('spectrumo3dbias.pdf', scale=scale, func=func), 
              ))
    



    cmd = '''\
python Image.py analyse "{fname}" \
--ohdu 3 \
--params-mode mean \
--output-file "{folder}/meanparams.csv" \
--plot-part "{folder}/meanpart*.pdf" \
--plot-sides "{folder}/meanside*.pdf" \
--plot-spectrum "{folder}/meanspectrum*.pdf" \
'''.format(fname=fname, folder=folder)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )
    
    doc.frame('Mean outliers analysis', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file='meanparams.csv', fontsize=5, spacing=6, func=func ) ), 
              )
        
    scale = .7
    doc.frame('sides', doc.center(
              doc.figure('meansideo3data.pdf', scale=scale, func=func ), doc.figure('meansideo3bias.pdf', scale=scale, func=func ), 
              doc.figure('meansideo3vbias.pdf', scale=.7*scale, func=func ), doc.figure('meansideo3dbias.pdf', scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('projections', doc.center(
              doc.figure('meansideo3proj1.pdf', scale=scale, func=func), 
              doc.figure('meansideo3proj0.pdf', scale=scale, func=func ),
              ))
    
    scale = .3
    doc.frame('spectra', doc.center(
              doc.figure('meanspectrumo3dbias.pdf', scale=scale, func=func), 
              doc.figure('meanspectrumo3bias.pdf', scale=scale, func=func), 
              doc.figure('meanspectrumo3vbias.pdf', scale=scale, func=func), 
              doc.figure('meanspectrumo3data.pdf', scale=scale, func=func), 
              ))



    pre = 'mean_remove40b3'
    cmd = '''\
python Image.py analyse "{fname}" \
--ohdu 3 \
--params-mode mean \
--remove-hits 40 3 \
--output-file "{folder}/{pre}_params.csv" \
--plot-sides "{folder}/{pre}_side*.pdf" \
--plot-spectrum "{folder}/{pre}_spectrum*.pdf" \
'''.format(fname=fname, folder=folder, pre=pre)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )
    
    doc.frame('removed above 40ADU with border 3 params', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file='mean_remove80b2_params.csv', fontsize=5, spacing=6, func=func ) ), 
              )
        
    scale = .7
    doc.frame('removed above 40ADU with border 3', doc.center(
              doc.figure('{pre}_sideo3data.pdf'.format(pre=pre), scale=scale, func=func ), 
              doc.figure('{pre}_sideo3bias.pdf'.format(pre=pre), scale=scale, func=func ), 
              doc.figure('{pre}_sideo3vbias.pdf'.format(pre=pre), scale=.7*scale, func=func ), 
              doc.figure('{pre}_sideo3dbias.pdf'.format(pre=pre), scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('removed above 40ADU with border 3 spectra', doc.center(
              doc.figure('{pre}_spectrumo3dbias.pdf'.format(pre=pre), scale=scale, func=func), 
              doc.figure('{pre}_spectrumo3bias.pdf'.format(pre=pre), scale=scale, func=func), 
              doc.figure('{pre}_spectrumo3vbias.pdf'.format(pre=pre), scale=scale, func=func), 
              doc.figure('{pre}_spectrumo3data.pdf'.format(pre=pre), scale=scale, func=func), 
              ))
    


    cmd = '''\
python Image.py analyse "{fname}" \
--ohdu 3 \
--params-mode mean \
--remove-hits 60 3 \
--output-file "{folder}/mean_remove60b3_params.csv" \
--plot-sides "{folder}/mean_remove60b3_side*.pdf" \
--plot-spectrum "{folder}/mean_remove60b3_spectrum*.pdf" \
'''.format(fname=fname, folder=folder)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )
    
    doc.frame('', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file='mean_remove60b3_params.csv', fontsize=5, spacing=6, func=func ) ), 
              )
        
    scale = .7
    doc.frame('removed above 60ADU with border 3', doc.center(
              doc.figure('mean_remove60b3_sideo3data.pdf', scale=scale, func=func ), doc.figure('mean_remove60b3_sideo3bias.pdf', scale=scale, func=func ), 
              doc.figure('mean_remove60b3_sideo3vbias.pdf', scale=.7*scale, func=func ), doc.figure('mean_remove60b3_sideo3dbias.pdf', scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('removed above 60ADU with border 3 spectra', doc.center(
              doc.figure('mean_remove60b3_spectrumo3dbias.pdf', scale=scale, func=func), 
              doc.figure('mean_remove60b3_spectrumo3bias.pdf', scale=scale, func=func), 
              doc.figure('mean_remove60b3_spectrumo3vbias.pdf', scale=scale, func=func), 
              doc.figure('mean_remove60b3_spectrumo3data.pdf', scale=scale, func=func), 
              ))
    

    cmd = '''\
python Image.py analyse "{fname}" \
--ohdu 3 \
--params-mode mean \
--remove-hits 80 2 \
--output-file "{folder}/mean_remove80b2_params.csv" \
--plot-sides "{folder}/mean_remove80b2_side*.pdf" \
--plot-spectrum "{folder}/mean_remove80b2_spectrum*.pdf" \
'''.format(fname=fname, folder=folder)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )
    
    doc.frame('image tool', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file='mean_remove80b2_params.csv', fontsize=5, spacing=6, func=func ) ), 
              )
        
    scale = .7
    doc.frame('removed above 80ADU with border 2', doc.center(
              doc.figure('mean_remove80b2_sideo3data.pdf', scale=scale, func=func ), doc.figure('mean_remove80b2_sideo3bias.pdf', scale=scale, func=func ), 
              doc.figure('mean_remove80b2_sideo3vbias.pdf', scale=.7*scale, func=func ), doc.figure('mean_remove80b2_sideo3dbias.pdf', scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('removed above 80ADU with border 2 spectra', doc.center(
              doc.figure('mean_remove80b2_spectrumo3dbias.pdf', scale=scale, func=func), 
              doc.figure('mean_remove80b2_spectrumo3bias.pdf', scale=scale, func=func), 
              doc.figure('mean_remove80b2_spectrumo3vbias.pdf', scale=scale, func=func), 
              doc.figure('mean_remove80b2_spectrumo3data.pdf', scale=scale, func=func), 
              ))
    


    pre = 'mean_remove150b1'
    cmd = '''\
python Image.py analyse "{fname}" \
--ohdu 3 \
--params-mode mean \
--remove-hits 150 1 \
--output-file "{folder}/{pre}_params.csv" \
--plot-sides "{folder}/{pre}_side*.pdf" \
--plot-spectrum "{folder}/{pre}_spectrum*.pdf" \
'''.format(fname=fname, folder=folder, pre=pre)
    print( cmd )
    func = lambda: subprocess.call( [cmd], shell=True )
    
    doc.frame('removed above 1500ADU with border 1 params', 
              'bash command',
              doc.code( cmd, 'Bash'),
              'parameters',
              doc.center( doc.table( file='mean_remove80b2_params.csv', fontsize=5, spacing=6, func=func ) ), 
              )
        
    scale = .7
    doc.frame('removed above 150ADU with border 1', doc.center(
              doc.figure('{pre}_sideo3data.pdf'.format(pre=pre), scale=scale, func=func ), 
              doc.figure('{pre}_sideo3bias.pdf'.format(pre=pre), scale=scale, func=func ), 
              doc.figure('{pre}_sideo3vbias.pdf'.format(pre=pre), scale=.7*scale, func=func ), 
              doc.figure('{pre}_sideo3dbias.pdf'.format(pre=pre), scale=.1*scale, func=func ),
              ))
    
    scale = .3
    doc.frame('removed above 150ADU with border 1 spectra', doc.center(
              doc.figure('{pre}_spectrumo3dbias.pdf'.format(pre=pre), scale=scale, func=func), 
              doc.figure('{pre}_spectrumo3bias.pdf'.format(pre=pre), scale=scale, func=func), 
              doc.figure('{pre}_spectrumo3vbias.pdf'.format(pre=pre), scale=scale, func=func), 
              doc.figure('{pre}_spectrumo3data.pdf'.format(pre=pre), scale=scale, func=func), 
              ))
    
