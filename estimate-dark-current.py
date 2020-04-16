import matplotlib.pyplot as plt
import ConnieImage
import argparse
import numpy as np
import scipy.stats
import Statistics as stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Estimate dark current value.')
    parser.add_argument('--runid', type=int, help='an integer for the runID')
    parser.add_argument('--ohdu', type=int, help='an integer for the OHDU')
    #parser.add_argument('--sum', dest='accumulate', action='store_const',
                    #const=sum, default=max,
                    #help='sum the integers (default: find the max)')

    kwargs = vars(parser.parse_args())
    print kwargs
    image = ConnieImage.read_full_image( ConnieImage.readImage_raw( kwargs['runid'], kwargs['ohdu'] ) )
    print image.shapes()
    print 'mean', image.mean_active_minus_overscan_line_by_line()
    print '<A-B>', image.median_active_minus_overscan_line_by_line()
    print 'median_sqr', image.median_sqr_active_minus_overscan_line_by_line()
    print 'sqr_median', image.sqr_median_active_minus_overscan_line_by_line()
    
    print '<A**2-B**2>', image.median_active_sqr_minus_overscan_sqr_line_by_line()
    print 'var', np.var(image.active) - np.var(image.overscan)
    #params = image.estimate_params()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.arange(-50,50,1)
    image.plot_spectrum( ax, bins, normed=True )
    sigma = np.mean(image.overscan_MAD_lines())
    ax.plot( bins, scipy.stats.norm.pdf(bins, loc=0, scale=sigma ))
    
    g = 10
    ax.plot( bins, stats.poisson_norm.pdf(bins, loc=0., scale=sigma, g=g, mu=1., fix_loc=True ))
    #ax.plot( bins, stats.poisson_norm.pdf(bins, loc=0., scale=np.mean(image.overscan_MAD_lines()), g = 2, mu=3 ))

    #ax.set_yscale('log')
    fig.savefig('spectrum.png')
    
