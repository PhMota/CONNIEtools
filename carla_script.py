import sys, os, glob, pyfits, numpy

for fname in glob.glob( '/share/storage2/connie/data/runs/047/*.fits*' ):
    for hdu in pyfits.open(fname):
        section = hdu.data[ 855:885, 3950:4050 ]
        print hdu.header['RUNID'], hdu.header['OHDU'], numpy.count_nonzero( section > 60 + numpy.median(section))
