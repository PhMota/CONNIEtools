import itertools

import numpy as np
from scipy import linalg, stats
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture


color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold', 'darkorange'])
#color_iter1 = itertools.cycle(['black'])


def plot_results(X, Y_, means, covariances, index, title):
    splot = plt.subplot(2, 1, 1 + index)
    for i, (mean, covar, color, color1) in enumerate(zip(
            means, covariances, color_iter, color_iter1)):
        print mean, covar
        v, w = linalg.eigh(covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color1)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

    plt.xlim(-9., 5.)
    plt.ylim(-3., 6.)
    plt.xticks(())
    plt.yticks(())
    plt.title(title)

def plot_results_1d(X, Y_, means, covariances, index, title):
    splot = plt.subplot(2, 1, 1 + index)
    _, bins, _ = plt.hist(X, bins=100 )
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        print mean, covar
        #v, w = linalg.eigh(covar)
        #v = 2. * np.sqrt(2.) * np.sqrt(v)
        #u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        x = .5*(bins[1:]+bins[:-1])
        plt.plot( x, 300*stats.norm.pdf(x,mean,covar[0]) )
        # Plot an ellipse to show the Gaussian component
        #angle = np.arctan(u[1] / u[0])
        #angle = 180. * angle / np.pi  # convert to degrees
        #plt.plot( ,  )
        #ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color1)
        #ell.set_clip_box(splot.bbox)
        #ell.set_alpha(0.5)
        #splot.add_artist(ell)

    #plt.xlim(-9., 5.)
    #plt.ylim(-3., 6.)
    #plt.xticks(())
    #plt.yticks(())
    plt.title(title)


# Number of samples per component
n_samples = 5000

# Generate random sample, two components
np.random.seed(0)
C = np.array([[0., -0.1], [1.7, .4]])
X = np.r_[np.dot(np.random.randn(n_samples, 2), C),
          .7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]

X = np.r_[np.random.randn(n_samples, 1),
          .7 * np.random.randn(n_samples, 1) - 6]


# Fit a Gaussian mixture with EM using five components
gmm = mixture.GaussianMixture(n_components=10, covariance_type='full').fit(X)
plot_results_1d(X, gmm.predict(X), gmm.means_, gmm.covariances_, 0,
             'Gaussian Mixture')

# Fit a Dirichlet process Gaussian mixture using five components
dpgmm = mixture.BayesianGaussianMixture(n_components=10, covariance_type='full').fit(X)
plot_results_1d(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, 1,
             'Bayesian Gaussian Mixture with a Dirichlet process prior')

print 'here'
plt.savefig('test.png')
#plt.show()
