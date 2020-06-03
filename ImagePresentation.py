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
--plot-part "{folder}/part*.pdf" \
--plot-sides "{folder}/side*.pdf" \
'''.format(fname=fname, folder=folder)
    print( cmd )
    subprocess.call( [cmd], shell=True )
    
    doc.frame('this frame', 
              'this text',
              doc.code( cmd, 'Bash'),
              )
    
    doc.frame('raw parts',
              doc.figure('parto3p1.pdf', height=0.2),
              doc.figure('parto3p2.pdf', height=0.2),
              doc.figure('parto3p3.pdf', height=0.2),
              doc.figure('parto3p4.pdf', height=0.2),
              )
    
    scale = .35
    doc.frame('sides',
              doc.column( doc.figure('sideo3data.pdf', scale=scale ), doc.figure('sideo3bias.pdf', scale=scale ),),
              doc.column( doc.figure('sideo3vbias.pdf', scale=scale ), doc.figure('sideo3dbias.pdf', scale=scale ),),
              )
    
    doc.frame('projections',
              doc.figure('sideo3proj1.pdf', scale=scale ), 
              doc.figure('sideo3proj0.pdf', scale=scale ),
              )
    
