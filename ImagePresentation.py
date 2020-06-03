from Beamer import *
import os
import subprocess
import sys
folder = 'ImagePresentation'
if not os.path.exists(folder):
    os.mkdir(folder)

with openBeamer( folder+'/'+folder ) as doc:
        
    fname = r'/share/storage2/connie/data/runs/029/\ runID_029_03326_Int-400_Exp-10800_11Mar18_18:06_to_11Mar18_21:10_p*.fits.fz'
    cmd = 'python Image.py analyse {fname} --ohdu 3 --plot-part "{folder}/part*.pdf"'.format(fname=fname, folder=folder)
    print( cmd )
    subprocess.call( [cmd], shell=True )
    
    doc.frame('this frame', 
              'this text',
              doc.code( cmd, 'Bash'),
              )
    
    doc.frame('raw parts',
              'part 1',
              doc.figure('parto3p1.pdf', height=1./4),
              'part 2',
              doc.figure('parto3p2.pdf', height=1./4),
              'part 3',
              doc.figure('parto3p3.pdf', height=1./4),
              'part 4',
              doc.figure('parto4p3.pdf', height=1./4),
              )
    
    
    
