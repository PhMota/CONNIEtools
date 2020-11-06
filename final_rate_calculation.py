from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *
from fstring import F

folder = '../final_rate_calculation'

print( F("folder") )
print( F("{{folder}}") )

if not os.path.exists(folder):
    os.mkdir(folder)

def makefigure( code, filename, doc, height=1. ):
    fullpath = F('{{folder}}/{{filename}}').str()
    print( fullpath )

    code += F(' --output "{{fullpath}}"').str()
    a = doc.code( code, language='Bash' )
    if not os.path.exists(fullpath):
        print()
        print( colored(code, 'green' ))
        print()
        subprocess.call( ["%s" % code], shell=True )
    if not os.path.exists(fullpath):
        print()
        print( '!!!not found', fullpath )
        exit(0)
    b = doc.figure( filename, height=height, center=True )
    return ''.join([a,b])

def requiredFile( code, file, doc ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(folder+'/'+file):
        print()
        print(code)
        print()
        subprocess.call( [code], shell=True )
    #b = doc.figure( file, height=height, center=True )
    return a

with Timer('presentation'):
    with openBeamer( folder, 'final rate calculation\\\\https://github.com/PhMota/CONNIEtools/final_rate_calculation.py' ) as doc:

        M = Expr('M')

        doc.frame('rate from number of events',
            Math( M == 2 ),
            # doc.itemize(
            #     r'divide off19 and off20 into 6 groups',
            #     r"\url{https://docs.google.com/document/d/1tiVyRtv_Dn8BBYy1wyq2Co941F3s86kslq-ZsAMm9l0/edit}",
            #     r"\url{https://connie-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=1937&filename=SpectrumOFF_EnergyBalanceCut_ohdu3.pdf&version=1}"
            # )
        )
