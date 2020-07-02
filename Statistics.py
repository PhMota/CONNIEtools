# coding: utf-8

from __future__ import print_function
from numpy import *
import scipy.stats
from scipy.special import erf
from Timer import Timer

def stdCovMed( x, y ):
    X = x - median(x)
    Y = y - median(y)
    return mean(X*Y)
    
def madCovSqr( data0, data1 ):
    return median( abs( (data0 - median(data0))*((data1 - median(data1))) ) )

def normFit( data, loc=None, scale=None ):
    return scipy.stats.norm.fit( data, floc=loc, fscale=scale )

def tiles( data, p ):
    m = median(data)
    if p - 1 == 0: return [m]
    return tiles(data[data<m], p-1 ) + [m] + tiles(data[data>m], p-1 )

def quartiles(x):
    return tiles(x,2)

def IQR( data, scale=1.34897 ):
    return scipy.stats.iqr( data )/scale

def MAD( data, axis = None, scale=1.4826 ):
    if axis == 0:
        median = median( data, axis = 0 )[None,:]
    elif axis == 1:
        median = median( data, axis = 1 )[:,None]
    else:
        median = median( data )
    return scale * median( abs( data - median ), axis = axis )

def nanMAD( data, axis = None, scale=1.4826 ):
    if axis == 0:
        median = nanmedian( data, axis = 0 )[None,:]
    elif axis == 1:
        median = nanmedian( data, axis = 1 )[:,None]
    else:
        median = nanmedian( data )
    return scale * nanmedian( abs( data - median ), axis = axis )

def FWHM( x, axis = None, f = 0.5 ):
    fwhm = apply_along_axis( lambda x: FWHM1d(x,f=f), axis, x )
    if axis == 1:
        sigma, center = fwhm.T
    else:
        sigma, center = fwhm
    return sigma, center

def FWHM1d( x, f=.5 ):
    y = array(x)
    if y.size == 0: return nan, nan
    yleft = None
    yright = None
    Max = None
    iteractions = 0
    while True:
        nbins = int(sqrt(y.size))
        yleft_prev, yright_prev = yleft, yright
        bins = linspace( nanmin(y), nanmax(y), nbins )
        hist, edges = histogram( y, bins )
        binwidth = bins[1] - bins[0]

        Max = max(hist)
        above = hist >= f*Max
        
        edge_left = nanmin( edges[:-1][above] )
        edge_right = nanmax( edges[1:][above] )
        yleft = mean( y[ abs(y - edge_left) < binwidth ] )
        yright = mean( y[ abs(y - edge_right) < binwidth ] )
        
        middle = logical_and( y>yleft, y<yright)

        if above[0] == True and above[-1] == True: 
            break
        y = y[ middle ]
        
        if not yleft_prev is None:
            if yleft_prev == yleft and yright_prev == yright: 
                break
        if y.size == 0: 
            print( 'FWHM: ValueWarning: did not converge' )
            return nan, nan
        iteractions += 1

    sigma = (abs(yleft - yright))/( 2*sqrt(2*log(1./f)) )
    center = mean([yleft,yright]) 
    return sigma, center 

def mean_outliers( data, axis=None, n=2 ):
    med = median(data,axis)
    MAD_ = MAD(data,axis)
    return mean( data[ abs(data-med)<n*MAD_ ] )

norm = scipy.stats.norm

def _norm_fit_( X, mu=0, sigma=1., fix_mu=False, fix_sigma=False, weights=1 ):
    if fix_mu:
        pdf = lambda X, sigma: norm.pdf(X, mu, sigma)
        p0 = (sigma,)
        bounds = ( (0, None) )
    elif fix_sigma:
        pdf = lambda X, mu: norm.pdf(X, mu, sigma)
        p0 = (mu,)
        bounds = ( (None, None) )
    else:
        pdf = lambda X, mu, sigma: norm.pdf(X, mu, sigma)
        p0 = (mu, sigma)
        bounds = ( (None,None), (0,None) )
    ret = maxloglikelihood( pdf, X, p0=p0, bounds=bounds, weights=weights )['x']
    return ret

norm.fit2 = _norm_fit_
norm.__primitive0__ = lambda x, mu, sigma: .5*special.erf((x - mu)/(sqrt(2)*sigma))
norm.__primitive1__ = lambda x, mu, sigma: mu*norm.__primitive0__(x, mu, sigma) - sigma**2*norm.pdf(x, mu, sigma)
norm.__primitive2__ = lambda x, mu, sigma: (mu**2+sigma**2)*norm.__primitive0__(x, mu, sigma) - sigma**2*(mu+x)*norm.pdf(x, mu, sigma)
    
def mean_norm( mu, sigma, a, b ):
    primitive = lambda x: .5*mu*special.erf((x - mu)/(sqrt(2)*sigma)) - sigma*exp(-(x - mu)**2/(2*sigma**2))/sqrt(2*pi)
    if a == -inf and b == inf:
        return mu
    elif a == -inf:
        return primitive(b) + .5*mu
    elif b == inf:
        return .5*mu - primitive(a)
    return primitive(b) - primitive(a)

def var_norm( mu, sigma, a, b ):
    primitive = lambda x: .5*(mu**2 + sigma**2) *special.erf((x - mu)/(sqrt(2)*sigma)) \
        - sigma*(mu+x)*exp(-(x - mu)**2/(2 *sigma**2))/sqrt(2*pi)
    if a == -inf and b == inf:
        return (mu**2 + sigma**2)
    elif a == -inf:
        return primitive(b) + .5*(mu**2 + sigma**2)
    elif b == inf:
        return .5*(mu**2 + sigma**2) - primitive(a)
    return primitive(b) - primitive(a)

def sum_norm( mu, sigma, a, b ):
    try:
        print('shapes', a.shape, b.shape, mu.shape, sigma.shape )
    except:
        pass
    primitive = lambda x: .5* special.erf((x - mu)/(sqrt(2)*sigma))
    return primitive(b) - primitive(a)


def _norm_fit_curve_( y, x, A=1, mu=0, sigma=1., fix_A = False, fix_mu=False, fix_sigma=False ):
    if fix_mu:
        pdf = lambda X, A, sigma: A*norm.pdf(X, mu, sigma)
        p0 = (A, sigma)
        bounds = ( (0, None) )
    elif fix_sigma:
        pdf = lambda X, A, mu: A*norm.pdf(X, mu, sigma)
        p0 = (A, mu)
        bounds = ( (None, None) )
    else:
        pdf = lambda X, A, mu, sigma: A*norm.pdf(X, mu, sigma)
        p0 = (A, mu, sigma)
        bounds = ( (0,-inf,0), (inf,inf,inf) )
    ret, pcorr = scipy.optimize.curve_fit( pdf, x, y, p0=p0, bounds=bounds )
    return ret
norm.fit_curve = _norm_fit_curve_

poisson = scipy.stats.poisson

def make_histogram( X, bins=None, binwidth=None ):
    if binwidth:
        bins = arange(min(X), max(X), binwidth)
    hist, edges = histogram( X, bins )
    x = .5*(edges[1:]+edges[:-1])
    return hist, x, edges[1]-edges[0]

binom = scipy.stats.binom
class binom_norm:
    @classmethod
    def pdf( cls, q, Qt, p, sigma_q ):
        Qt = int(Qt)
        if sigma_q == 0:
            return binom.pmf( q.astype(int), Qt, p )
        Q = np.arange(Qt).astype(int)[:,None]
        return sum( binom.pmf( Q, Qt, p[None,:] ) * norm.pdf( q[None,:] - Q, scale = sigma_q ), axis = 0 )

    @classmethod
    def pdf2( cls, q, Qt, p, sigma_q ):
        Qt = int(Qt)
        if sigma_q == 0:
            return binom.pmf( q.astype(int), Qt, p )
        Q = np.arange(Qt).astype(int)[:,None]
        __b = binom.pmf( Q, Qt, p[None,:] )
        __n = norm.pdf( q[None,:] - Q, scale = sigma_q )
        return sum( __b*__n, axis = 0 )
    
class norm2d_binom_norm:
    @classmethod
    def pdf( cls, q, Qt, sigma_q, x, y, mux, muy, sigmax, sigmay ):
        return binom_norm.pdf( q, Qt, self.__p(x, mux, sigmax)*self.__p(y, muy, sigmay), sigma_q )
    
    @classmethod
    def __p( cls, x, mu, sigma ):
        return -.5*( erf( -( x+1 - mu)/( sqrt2*sigma )) - erf( -( x - mu )/( sqrt2*sigma ) ) )
    
    @classmethod
    def fit( cls ):
        pass

def integral_norm_pixel( x, mu, sigma ):
    f = 1./( sqrt(2)*sigma )
    #return -.5*( erf( -f*( x+1 - mu) ) - erf( -f*( x - mu ) ) )
    return .5*( erf( f*( x+1 - mu) ) - erf( f*( x - mu ) ) )

def integral_x_norm_pixel( x, mu, sigma ):
    f = 1./( sqrt(2)*sigma )
    a = .5*mu * erf( f*(x - mu) ) - sigma/sqrt(2*pi) * exp( -f**2*(x - mu)**2 )
    b = .5*mu * erf( f*(x+1 - mu) ) - sigma/sqrt(2*pi) * exp( -f**2*(x+1 - mu)**2 )
    return b - a

def integral_xSqr_norm_pixel( x, mu, sigma ):
    f = 1./( sqrt(2)*sigma )
    a = .5*(mu**2+sigma**2) * erf( f*(x - mu) ) - sigma*(mu+sigma)/sqrt(2*pi) * exp( -f**2*(x - mu)**2 )
    b = .5*(mu**2+sigma**2) * erf( f*(x+1 - mu) ) - sigma*(mu+sigma)/sqrt(2*pi) * exp( -f**2*(x+1 - mu)**2 )
    return b - a

    
class norm2d_norm:
    mu_bounds = (None,None)
    sigma_bounds = ( 0., None)
    @classmethod
    def pdf( cls, x, y, e, x_mu, y_mu, x_sigma, y_sigma, e_sigma, E ):
        x_p = integral_norm_pixel( x, x_mu, x_sigma )
        y_p = integral_norm_pixel( y, y_mu, y_sigma )
        P = sum( x_p*y_p )
        E = sum( e )
        if e_sigma == 0:
            return x_p*y_p*e
        e_p = E*x_p*y_p * norm.pdf( e/E, x_p*y_p, e_sigma/E )
        return e_p

    @classmethod
    def E( cls, x, y, e, x_mu, y_mu, x_sigma, y_sigma ):
        x_p = integral_norm_pixel( x, x_mu, x_sigma )
        y_p = integral_norm_pixel( y, y_mu, y_sigma )
        return sum( e )/sum( x_p*y_p )
        
    @classmethod
    def mle( cls, x, y, e, x_mu, y_mu, x_sigma, y_sigma, e_sigma, E, fix_e_sigma=False, fix_mu=False, fix_sigma=False, fix_E=False ):

        if fix_e_sigma:
            func = lambda p: \
                -sum( log( cls.pdf( x, y, e, p[0], p[1], p[2], p[3], p[4], e_sigma ) ) )
            bounds = ( (x_mu-1, x_mu+1), (y_mu-1, y_mu+1), (x_sigma-.2, x_sigma+.2), (y_sigma-.2, y_sigma+.2), (E, None) )
            p0 = ( x_mu, y_mu, x_sigma, y_sigma, E )

        if fix_e_sigma and fix_mu:
            func = lambda p: \
                -sum( log( cls.pdf( x, y, e, x_mu, y_mu, p[0], p[1], p[2], e_sigma ) ) )
            bounds = ( cls.sigma_bounds, cls.sigma_bounds, (E, None) )
            p0 = ( x_sigma, y_sigma, E )

        if fix_e_sigma and fix_sigma and fix_E:
            func = lambda p: \
                -sum( log( cls.pdf( x, y, e, p[0], p[1], x_sigma, y_sigma, E, e_sigma ) ) )
            bounds = ( (x_mu-1, x_mu+1), (y_mu-1, y_mu+1) )
            p0 = ( x_mu, y_mu )

        if fix_e_sigma and fix_mu and fix_E:
            func = lambda p: \
                -sum( log( cls.pdf( x, y, e, x_mu, y_mu, p[0], p[1], E, e_sigma ) ) )
            bounds = ( cls.sigma_bounds, cls.sigma_bounds )
            p0 = ( x_sigma, y_sigma )
        
        #print( 'bounds', bounds )
        #print( 'p0', p0 )
        #f = lambda *params: -sum( log( func( x, y, e, *params )) )
        ret = scipy.optimize.minimize( func, p0, bounds=bounds )
        #print( ret['fun'] )
        return ret['x']

    @classmethod
    def __fit_func( cls, xy, xmu, ymu, xsigma, ysigma, E ):
        p = integral_norm_pixel( xy[0], xmu, xsigma )*integral_norm_pixel( xy[1], ymu, ysigma )
        return E*p

    @classmethod
    def __fit_func_E( cls, xy, xmu, ymu, xsigma, ysigma, E ):
        p = integral_norm_pixel( xy[0], xmu, xsigma )*integral_norm_pixel( xy[1], ymu, ysigma )
        return E*sum(p)*p

    @classmethod
    def __fit_func_E_sigma( cls, xy, xmu, ymu, xsigma, ysigma, E ):
        p = integral_norm_pixel( xy[0], xmu, xsigma )*integral_norm_pixel( xy[1], ymu, ysigma )
        p = integral_norm_pixel( xy[0], xmu, xsigma*sum(p) )*integral_norm_pixel( xy[1], ymu, ysigma*sum(p) )
        return E*sum(p)*p

    @classmethod
    def __fit_func_E_sigma2( cls, xy, xmu, ymu, xsigma, ysigma, E ):
        pxysigma = sum( integral_xSqr_norm_pixel( xy[0] - xmu, 0, xsigma )*integral_xSqr_norm_pixel( xy[1] - ymu, 0, ysigma ) )
        A = xsigma**2 * ysigma**2
        pxsigma = pxysigma/ysigma**2
        pysigma = pxysigma/xsigma**2
        #print( 'xsigma', _xsigma_**2, _pxsigma_, _xsigma_ - sqrt(_pxsigma_) )
        #print( 'ysigma', _ysigma_**2, _pysigma_, _ysigma_ - sqrt(_pysigma_) )
        deltaxsigma = 0 #_xsigma_ - sqrt(_pxsigma_)
        deltaysigma = 0 #_ysigma_ - sqrt(_pysigma_)
        #a = _xsigma_*( 1 - _xsigma_
        p = integral_norm_pixel( xy[0], xmu, xsigma )*integral_norm_pixel( xy[1], ymu, ysigma )
        p = integral_norm_pixel( xy[0], xmu, xsigma*sum(p) )*integral_norm_pixel( xy[1], ymu, ysigma_*sum(p) )
        return E*sum(p)*p

    @classmethod
    def fit( cls, x, y, e, xmu, ymu, xsigma, ysigma, E, sigma_e, ftol=1e-8, mode=None, loss='linear' ):
        if mode == None or mode == '':
            func = cls.__fit_func
        elif mode == 'E':
            func = cls.__fit_func_E
        elif mode == 'Esigma':
            func = cls.__fit_func_E_sigma
        mins = (xmu-2, ymu-2, 0.01, 0.01, E )
        maxs = (xmu+2, ymu+2, xsigma+.5, ysigma+.5, inf )
        try:
            p, pcov = scipy.optimize.curve_fit( func, [x, y], e, bounds=( mins, maxs ), sigma=sigma_e, ftol=ftol, loss=loss )
        except RuntimeError:
            p = [xmu, ymu, xsigma, ysigma, E]
        return p


def negloglikelihood_fast( ePix, xPix, yPix, E, mu, sigma, sigma_noise, single_sigma = False ):
    integral_norm_pixel_x = -.5*( erf( -( xPix+1 - mu[0])/( sqrt2*sigma[0] )) - erf( -( xPix - mu[0] )/( sqrt2*sigma[0] ) ) )
    integral_norm_pixel_y = -.5*( erf( -( yPix+1 - mu[1])/( sqrt2*sigma[1] )) - erf( -( yPix - mu[1] )/( sqrt2*sigma[1] ) ) )
    
    prob = integral_norm_pixel_x * integral_norm_pixel_y
    prob_i = prob[:,None]
    ePix_i = ePix[:,None]
    q_j = np.arange(E).astype(int)[None,:]
    j = -1
    val = np_sum( binom_pmf( q_j, E, prob_i ) * norm_pdf( ePix_i - q_j, scale = sigma_noise ), axis = j )
    val = val[ val > 0 ]
    return -np.nansum( np.log( val ) )

class poisson_norm:
    verbose = False
    mu_bounds = (None,None)
    sigma_bounds = (0,None)
    lamb_bounds = (1e-4,1)
    gain_bounds = (1,None)
    
    @classmethod
    def pdf3(cls, x, mu=0, sigma=1, gain=1, lamb=1, tol=1e-8, debug=False ):
        result = zeros_like(x).astype(float)
        
        k = 0
        poisson_sum = 0
        while True:
            poisson_weight = poisson.pmf( k, lamb )
            if poisson_weight != poisson_weight: 
                break
            poisson_sum += poisson_weight
            result += poisson_weight * norm.pdf( x, mu+gain*k, sigma )
            if poisson_weight/poisson_sum < tol: break
            k += 1
        return result

    @classmethod
    def pdf2(cls, x, mu=0, sigma=1, gain=1, lamb=1, tol=1e-8 ):
        result = zeros_like(x).astype(float)

        k = 0
        poisson_sum = 0
        poisson_weights = []
        while True:
            poisson_weights.append( poisson.pmf( k, mu=lamb ) )
            poisson_sum += poisson_weights[-1]
            if poisson_weights[-1]/poisson_sum < tol: break
        res = sum([ w*norm.pdf( x, loc = mu + gain*k, scale = sigma ) for k, w in enumerate(poisson_weights) ])
        return result


    @classmethod
    def pdf(cls, x, mu=0, sigma=1, gain=1, lamb=1, tol=1e-6):
        kmax = int( (log(tol)+lamb)/log(lamb) ) + 1
        ks = arange(0, kmax, 1).astype(int)
        poisson_weights = poisson.pmf( ks, mu=lamb )
        
        return sum( [ poisson_weights[k] * norm.pdf( x, mu+gain*k, sigma )[None,:] for k in ks ], axis=0 )

    @classmethod
    def rvs( cls, loc=0, scale=1, g=1, mu=1, size=1 ):
        a = loc - 4*scale + g*mu**2
        b = loc + 4*scale + g*mu**2
        return monteCarlo_rvs( lambda x: cls.pdf(x, loc, scale, g, mu), a, b, size, normed=False, verbose=cls.verbose )

    @classmethod
    def rvs2( cls, loc=0, scale=1, g=1, mu=1, size=1, tol=1 ):
        k = 0
        #poisson_sum = 0
        res = []
        while True:
            n = size*poisson.pmf( k, mu=mu )
            if n < tol: return res
            N = poisson.rvs( n )
            res.extend( norm.rvs( loc=k*g*mu, scale=scale, size=N ) )
            k += 1

    @classmethod
    def fit_binned( cls, X, mu=0, sigma=1., gain=1, lamb=1e-3, fix_mu=False, fix_sigma=False, fix_g=False, fix_lamb=False, binsize=1 ):
        bins = arange(nanmin(X), nanmax(X), binsize)
        hist, edges = histogram( X, bins )
        x = (edges[:-1]+edges[1:])/2.
        
        params = scipy.optimize.curve_fit()
        
    @classmethod
    def fit( cls, X, mu=0, sigma=1., gain=1, lamb=1e-3, fix_mu=False, fix_sigma=False, fix_g=False, fix_lamb=False ):
        func = cls.pdf
        pdf = None
        if fix_mu and fix_sigma:
            pdf = lambda X, gain, lamb: func(X, mu+gain*lamb, sigma, gain, lamb)
            p0 = (gain, lamb)
            bounds = (cls.gain_bounds, cls.lamb_bounds)
        elif fix_sigma:
            pdf = lambda X, mu, gain, lamb: func(X, mu+gain*lamb, sigma, gain, lamb)
            p0 = (mu, gain, lamb)
            bounds = (cls.mu_bounds, cls.gain_bounds, cls.lamb_bounds)
        elif not fix_mu and not fix_sigma and not fix_g and not fix_lamb:
            pdf = lambda X, mu, sigma, gain, lamb: func(X, mu+g*lamb, sigma, gain, lamb)
            p0 = (mu, sigma, gain, lamb)
            bounds = (cls.mu_bounds, cls.sigma_bounds, cls.gain_bounds, cls.lamb_bounds)
        if pdf is None:
            raise Exception('fix combination not implemented')
        ret = maxloglikelihood( pdf, X, p0=p0, bounds=bounds )['x']
        return ret

sqrt2 = sqrt(2.)
class binom_norm2d:
    def integral_norm_pixel( pix_edge, mu, sigma ):
        return -.5*( scipy.special.erf( -( pix_edge+1 - mu)/( sqrt2*sigma )) - scipy.special.erf( -( pix_edge - mu )/( sqrt2*sigma ) ) )

    def probability_pixel( xPix, yPix, mu, sigma, single_sigma = False ):
        if single_sigma: sigma[1] = sigma[0]
        return integral_norm_pixel( xPix, mu = mu[0], sigma = sigma[0] ) * integral_norm_pixel( yPix, mu = mu[1], sigma = sigma[1] )

    def pdf( ePix, xPix, yPix, E, mu_xy, sigma_xy, sigma_E, single_sigma = False):
        prob = probability_pixel( xPix, yPix, mu_xy, sigma_xy, single_sigma )
        if sigma_E == 0:
            return scipy.stats.binom.pmf( ePix, E, prob )
        prob_i = prob[None,:]
        ePix_i = ePix[None,:]
        q_j = arange(E)[:,None]
        j = 0
        return sum( scipy.stats.binom.pmf( q_j, E, prob_i ) * scipy.stats.norm.pdf( ePix_i - q_j, scale = sigma_noise ), axis = j )

    def negloglikelihood( ePix, xPix, yPix, E, mu_xy, sigma_xy, sigma_E, single_sigma = False ):
        val = pdf( ePix, xPix, yPix, E, mu_xy, sigma_xy, sigma_E, single_sigma = single_sigma )
        #val = val[ val > 0 ]
        return -nansum( log( val ) )

    def fit( ePix, xPix, yPix, E0, mu_xy0, sigma_xy0, sigma_E0=None, single_sigma = False ):
        Q, mu[0], mu[1], sigma[0], sigma[1] = scipy.optimize.minimize(
        fun = fun,
        x0 = [ int(E0), mu0[0]+.5, mu0[1]+.5, sigma0[0], sigma0[1] ],
        tol = tol,
        bounds = [(1, inf), (-inf,inf), (-inf, inf), (0.001, inf), (0.001, inf)] ).x


def monteCarlo_rvs( pdf, a, b, N, normed=False, verbose=False ):
    if not normed:
        pdf_max = -scipy.optimize.fmin( lambda x: -pdf(x), 0, full_output=True )[1]
        if verbose: print( 'pdf_max', pdf_max )
        normedpdf = lambda x: pdf(x)/pdf_max
    else:
        normedpdf = pdf
    results = []
    if verbose: print( '(random', end=' ' )
    while 1:
        x = (b - a) * random.random_sample( N ) + a
        tests = random.random_sample( N )
        y = normedpdf(x)
        mask = tests < y
        results = concatenate( (results, x[mask]) )
        if len(results) >= N: break
        if verbose: print( '%.2f'%( float(len(results))/N ), end=' ' )
    if verbose: print( ')' )
    return results[:N]

def negloglikelihood( pdf, params, X, weights ):
    res = -sum( weights*log( pdf( X, *params )) )
    return res

def maxloglikelihood( pdf, X, p0, weights=1, jac=False, bounds=None ):
    return scipy.optimize.minimize( lambda p: negloglikelihood( pdf, p, X, weights ), p0, jac=jac, bounds=bounds, tol=1e-8 )

def chisquare( y, f, yerr=None, ddof=1 ):
    if yerr is None:
        yerr = sqrt(y)
    return sum( ((y-f)/yerr)**2 )/(len(y)-ddof)

def errSqr_mean( x ):
    return var(x)/len(x)

def errSqr_std(x):
    return mean( (x-mean(x))**4 )/var(x)/len(x)

def errSqr_median(x):
    return MAD(x)**2/len(x)

def errSqr_MAD(x):
    return MAD(x)**2/len(x)
    
def sign(x):
    return -1 if x<0 else ( 1 if x>0 else 0 )

def pprint( msg, f ):
    print( msg )
    return f

    def __pow__( self, a ):
        unit = {}
        for key, p in self.unit.items():
            if int(a*p) != a*p:
                raise Exception( 'unit power %s not allowed'%(a*p)) 
            unit[key] = a*p
        return dfloat( self.value**a, unit={ unit: a*p for unit,p in self.unit.items() } )
    
    
    def __div__( self, a ): return self * a**(-1)
    __truediv__ = __div__
    
    def __rdiv__(self,a): return self**(-1) * a
    __rtruediv__ = __rdiv__
    
    def __abs__(self): return dfloat(abs(self.value),self.unit)

    def __float__(self):
        if self.unit != {}:
            raise Exception('cannot convert %s to float'%self.unit )
        return self.value
    
    def __trunc__(self): 
        return int(self.value)
    
    def mean(self):
        if isinstance(self.value[0], ufloat):
            val = array([ entry.val for entry in self.value ])
            w = array([ 1./entry.errSqr for entry in self.value ])
            s = sum(w)
            return dfloat( ufloat( sum(w*val)/s, errSqr=1./s ), unit=self.unit )
        
        return dfloat( ufloat( mean(self.value), errSqr=errSqr_mean(self.value) ), unit=self.unit )

    def median(self):
        return dfloat( ufloat( median(self.value), errSqr=errSqr_median(self.value) ), unit=self.unit )

    def std(self):
        if isinstance(self.value[0], ufloat):
            val = array([ entry.val for entry in self.value ])
            w = array([ 1./entry.errSqr for entry in self.value ])
            s = sum(w)
            return dfloat( ufloat( sqrt( (sum(w*val**2) - sum(w*val)**2)/s ), errSqr=1./s ), unit=self.unit )
        
        return dfloat( ufloat( std(self.value), errSqr=errSqr_std(self.value) ), unit=self.unit )

    def MAD(self):
        return dfloat( ufloat( MAD(self.value), errSqr=errSqr_MAD(self.value) ), unit=self.unit )
    
#for func_name in [ 'sum', '__setitem__' ]:
    #print 'set func', func_name
    #attr = lambda self, func_name=func_name, *args, **kwargs: dfloat( getattr(type(self.value), func_name)(self.value, *args, **kwargs), self.unit )
    #setattr( dfloat, func_name, attr )

def mantissa_exp(x):
    m, e = map( float, ('%e' % (x)).split('e') )
    return m, int(e)

def significant_decimal(x):
    return mantissa_exp(x)[1]

def round2significant( x, sig ):
    return around( x, -sig )
    
class ufloat:
    '''
    class for uncertainties and dimension
    '''
    __name__ = 'ufloat'
    sep = '+-'
    exp = 'e'
    
    __units__ = []
    __conversion_matrix__ = []
    __multipliers__ = { 'f': -15, 'p': -12, 'n': -9, 'µ': -6, '': 0, 'm': -3, 'k': 3, 'M': 6, 'G': 9, 'T': 12, 'P': 15 }
    __powers__ = { v: k for k, v in __multipliers__.items()}
    
    
    def __init__( self, val, err=0, errSqr=None, errRel=None, unit={} ):
        if isinstance( val, ufloat ):
            self.val = val.val
            self.errSqr = val.errSqr
        else:
            self.val = val
            if not errSqr is None:
                self.errSqr = float(errSqr)
            elif not errRel is None:
                self.errSqr = float(val*errRel)**2
            else:
                self.errSqr = float(err)**2
            
        if isinstance( unit, dict ):
            self.unit = unit
        elif isinstance( unit, str ):
            if unit[0] in ufloat.__multipliers__.keys():
                self.val *= 10**(ufloat.__multipliers__[unit[0]])
                self.errSqr *= 10**(2*ufloat.__multipliers__[unit[0]])
                self.unit = {unit[1:]: 1}
            else:
                self.unit = {unit: 1}
        else:
            raise Exception('unit type not recognized: %s %s' % ( type(unit), unit ) )

    def tuple(self):
        if isinstance(self.value, ndarray):
            return map( lambda x: dfloat(x,self.unit).tuple(), self.value )
        if self.unit == {}:
            return (self.value, '')
        power_pair = ('',0)
        if float(self.value) != 0:
            log10 = math.log10( abs(float(self.value)) )
            power = 3*int((log10/self.unit[self.unit.keys()[0]]+1)/3.)
            if power != 0:
                #print 'power', power
                power_pair = ( dfloat.powers[power], power*self.unit[self.unit.keys()[0]] )
        unit_str = power_pair[0]+'·'.join([ '%s^%s'%(unit,p) if p!= 1 else '%s'%unit for unit,p in self.unit.items()if p!=0] )
        return (self.value/10**power_pair[1], unit_str)
    
    def unit_str(self):
        return '·'.join([ '%s^%s'%(unit,p) if p!= 1 else '%s'%unit for unit,p in self.unit.items()if p!=0] )
        
    def __repr__(self):
        if isinstance(self.value, ndarray):
            return '%s%s'%(self.value, self.unit_str())
        return '%s%s' % self.tuple()

    @classmethod
    def add_conversion(cls, unit1, unit2, value ):
        print( 'add', unit1, unit2, value )
        
        new_units = list( cls.__units__ )
        N = len(cls.__units__)
        if unit1 in cls.__units__ and unit2 in cls.__units__:
            raise Exception( 'conversion already set' )
        if not unit1 in new_units:
            new_units.append( unit1 )
        if not unit2 in new_units:
            new_units.append( unit2 )
            
        new_conversion_matrix = zeros(( len(new_units), len(new_units) )).tolist()
        for i1 in range(N):
            for i2 in range(N):
                new_conversion_matrix[i1][i2] = cls.__conversion_matrix__[i1][i2]
        
        i1 = new_units.index( unit1 )
        i2 = new_units.index( unit2 )
        
        new_conversion_matrix[i1][i1] = 1.
        new_conversion_matrix[i2][i2] = 1.
        new_conversion_matrix[i1][i2] = value
        new_conversion_matrix[i2][i1] = 1./value
        
        for i, u in enumerate(new_units):
            if u != unit1 and u != unit2:
                if new_conversion_matrix[ i2 ][ i ] == 0:
                    new_conversion_matrix[ i2 ][ i ] = new_conversion_matrix[ i2 ][ i1 ] * new_conversion_matrix[ i1 ][ i ]
                    new_conversion_matrix[ i ][ i2 ] = 1./new_conversion_matrix[ i2 ][ i ]
                if new_conversion_matrix[ i1 ][ i ] == 0:
                    new_conversion_matrix[ i1 ][ i ] = new_conversion_matrix[ i1 ][ i2 ] * new_conversion_matrix[ i2 ][ i ]
                    new_conversion_matrix[ i ][ i1 ] = 1./new_conversion_matrix[ i1 ][ i ]
        
        cls.__conversion_matrix__ = new_conversion_matrix
        cls.__units__ = new_units
    
    @classmethod
    def show_conversion_matrix(cls):
        print( ' '*5 + ' '.join([ '%5s'%u for u in cls.__units__ ]) )
        for i1, u1 in enumerate(cls.__units__):
            print( '%4s '%u1 + ' '.join([ '%.3f' % cls.__conversion_matrix__[i1][i2] for i2,u2 in enumerate(cls.__units__ )]) )
        return
    
    def asunit( self, unit ):
        new_unit = {}
        if unit not in ufloat.__units__:
            raise Exception('%s is not in conversion_map:\n%s' % ( unit, ufloat.show_conversion_matrix() ) )

        for key in self.unit.keys():
            if key in ufloat.__units__:
                new_unit[unit] = self.unit[key]
                value = self * ufloat.__conversion_matrix__[ ufloat.__units__.index(key) ][ ufloat.__units__.index(unit) ]**(self.unit[key])
            else:
                raise Exception('cannot convert %s into %s: \n%s' % (self.unit, unit, ufloat.show_conversion_matrix()) )
        
        return ufloat( value, unit = new_unit )

    def __unit_str__(self):
        return '·'.join([ '%s^%s' % ( unit, power ) if power!= 1 else '%s' % unit for unit, power in self.unit.items() if power!=0] )

    def pprint( self ):
        '''
        returns a human friendly string of the form (using the sep string that defaults for '+-'. The values are rounded using the nearest even rule.
        ufloat(5,.1) --> (5.0+-0.1)
        ufloat(50,100) --> (0.0+-0.1)k
        if the unit is given:
        ufloat(3.1416,.01,unit='rad') --> (3.14+-0.01)rad
        ufloat(100,50,unit='km') --> (0.10+-0.05)Mm
        '''
        err_decimal = significant_decimal( self.err() )
        val_decimal = significant_decimal( self.val )
        multiplier = int( (val_decimal+1)/3 )
        err_decimal -= multiplier*3
        val_decimal -= multiplier*3
        if err_decimal > 0:
            multiplier = int( (val_decimal+3)/3 )
            err_decimal -= multiplier*3
            val_decimal -= multiplier*3
        if self.errSqr == 0:
            return '%s%s%s'%( self.val/10**(3*multiplier), ufloat.__powers__[multiplier*3], self.__unit_str__() )
        pattern = '(%%.%df%s%%.%df)%s%s'%( -err_decimal, ufloat.sep, -err_decimal, ufloat.__powers__[multiplier*3], self.__unit_str__() )
        return pattern % ( round2significant( self.val/10**(3*multiplier), err_decimal ), round2significant( self.err()/10**(3*multiplier), err_decimal ) )

    def __repr__( self ):
        return self.pprint()
    
    def __str__( self ):
        return self.__repr__()
    
    def set_err(self, err):
        '''
        set squared err as err**2
        '''
        self.errSqr = err**2
    
    def set_errSqr( self, errSqr ):
        '''
        set squared err as errSqr
        '''
        self.errSqr = errSqr
    
    def err(self):
        '''
        get error as sqrt(errSqr)
        '''
        return math.sqrt(self.errSqr)
    
    def errRel(self):
        if self.val == 0: return 0
        return math.sqrt(self.errSqr)/self.val
    
    def __neg__( self ):
        return ufloat( -self.val, errSqr=self.errSqr, unit=self.unit )
    
    def __add__( self, u ):
        if u is self:
            raise Exception( 'cannot add same variable due to correlations' )
        if isinstance( u, ufloat ):
            if self.unit == u.unit:
                return ufloat( self.val + u.val, errSqr = self.errSqr + u.errSqr, unit=self.unit )
            raise Exception( 'cannot add units %s and %s' % (self.unit, a.unit) )
        raise Exception( 'cannot add types %s and %s' % (self.unit, type(a)) )
    
    def __sub__( self, u ):
        return self + (-u)

    def __mul__( self, u ):
        if isinstance( u, ufloat ):
            unit = {}
            for key in set( self.unit.keys() ) & set( u.unit.keys() ):
                p = self.unit[ key ] + u.unit[ key ]
                if p != 0: unit[ key ] = p
            for key in set(self.unit.keys()) - set( u.unit.keys() ):
                unit[ key ] = self.unit[ key ]
            for key in set( u.unit.keys() ) - set( self.unit.keys() ):
                unit[ key ] = u.unit[ key ]
            
            return ufloat( self.val * u.val, errSqr = u.val**2*self.errSqr + self.val**2*u.errSqr, unit = unit )
        return ufloat( u*self.val, errSqr = u**2*self.errSqr, unit = self.unit )
    
    __rmul__ = __mul__
    
    def __div__(self, u ):
        return self * u**(-1)
    
    __truediv__ = __div__
    
    def __rdiv__(self, u ):
        return self**(-1) * u
        
    def __pow__( self, u ):
        '''
        delta[a^u]^2 = delta[exp(u*ln(a))]^2 = ( u*a^(u-1) )**2 * delta[a]^2 + (ln(a)*a^u)**2 * delta[u]^2
        '''
        unit = { key: u*value for key, value in self.unit.items() }
        if isinstance( u, ufloat ):
            if not u.unit is {}:
                raise Exception( 'cannot power to unit %s' % u.unit )
            return ufloat( self.val**u.val, errSqr = (u.val*self.val**(u.val-1))**2*self.errSqr + ( math.log(self.val)*self.val**u.val )**2*u.errSqr, unit = unit )
        return ufloat( self.val**u, errSqr = (u*self.val**(u-1))**2*self.errSqr, unit = unit )

    #def __rpow__( self, u ):
        #if not isinstance( u, ufloat):
            #u = ufloat(u)
        #return ufloat( u.val**self.val, errSqr = (self.val*u.val**(self.val-1))**2*u.errSqr + ( math.log(u.val)*u.val**self.val )**2*self.errSqr )
    
    def __abs__(self):
        return ufloat( abs(self.val), errSqr=self.errSqr, unit = self.unit )
    
    def sqrt(self):
        return self**.5
    
    def floor(self,n=1):
        return self.val - n*self.err()

    def ceil(self,n=1):
        return self.val + n*self.err()
    
    def range(self,n=1):
        return (self.floor(n),self.ceil(n))
    
    def __eq__(self, a):
        c = self-a
        return abs( c.val ) < c.err()
    
    #def __float__(self):
        #return float(self.ceil())

#ufloat.add_conversion( 'e-', 'eV', 3.745 )

class uarray:
    '''
    array of uncertainties
    '''
    def __init__( self, value = None, err = None, errSqr = None, errMode = None, unit = None ):
        '''
        initiate with list of uncertainties
        '''
        if isinstance( value, ndarray ):
            self.unit = unit
            self.val = value
            if errMode == 'poisson':
                self.errSqr = array(self.val)
            elif errMode == 'stats':
                self.errSqr = array(var(self.val))/len(self.val)
            return
        
        self.val = []
        self.errSqr = []
        self.unit = {}
        if not list_ is None:
            for i, entry in enumerate(list_):
                if i == 0:
                    if isinstance( entry, ufloat ):
                        self.unit = entry.unit 
                    else:
                        self.unit = ufloat(entry,unit=unit).unit
                if isinstance( entry, ufloat ):
                    self.val.append( entry.val )
                    self.errSqr.append( entry.errSqr )
                    if not self.unit is entry.unit:
                        raise Exception( 'unmatching units: %s and %s' % ( self.unit, entry.unit ) )
                else:
                    u = ufloat( entry, unit=unit )
                    self.val.append( u.val )
                    self.errSqr.append( u.errSqr )
        
        self.val = array(self.val)
        self.errSqr = array(self.errSqr)
        self.summary()
    
    def set_errSqr( self, errSqr ):
        self.errSqr = zeros_like(self.val)*errSqr
        
    def std_errSqr(self):
        #print 'len', var(self.val), len(self.val), (var(self.val)/len(self.val))**.5
        return var(self.val)/len(self.val)

    def weights(self):
        errSqr = array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        return 1./errSqr/self.weightSum()
    
    def weightSum(self):
        errSqr = array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        return sum( 1./errSqr )
    
    def mean(self):
        errSqr = array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        s = 1./self.weightSum()
        return ufloat( sum(self.val/errSqr)*s, errSqr=s, unit=self.unit )
    
    def mean_list(self):
        sum_ = ufloat(0,0,unit=self.unit)
        for entry in self.tolist():
            sum_ = sum_ + entry
        return sum_/len(self.val)
    
    def mean_stats(self):
        return ufloat( mean(self.val), errSqr = var(self.val)/len(self.val), unit = self.unit )
    
    def std(self):
        errSqr = array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        s = 1./self.weightSum()
        return ufloat( (sum(self.val**2/errSqr)*s - (sum(self.val/errSqr)*s)**2)**.5, errSqr=s, unit=self.unit )

    def std_list(self):
        mean = self.mean_list()
        l = self.tolist()
        sum_diffSqr = (l[0] - mean)**2
        for entry in l[1:]:
            sum_diffSqr = sum_diffSqr + (entry-mean)**2
        return (sum_diffSqr/len(self.val))**.5

    def std_list2(self):
        mean = self.mean_list()
        l = self.tolist()
        sum_Sqr = l[0]**2
        for entry in l[1:]:
            sum_Sqr = sum_Sqr + entry**2
        return (sum_Sqr/len(self.val) - mean**2)**.5
    
    def std_stats(self):
        return ufloat( std(self.val), errSqr = mean((self.val-mean(self.val))**4)/len(self.val)/var(self.val), unit = self.unit )
    
    def summary(self):
        print( 'uarray', self )
        print( 'mean', self.mean(), self.mean_stats(), self.mean_list() )
        print( 'std', self.std(), self.std_stats(), self.std_list(), self.std_list2() )
        
    def append(self, entry):
        self.val = self.val.tolist()
        self.errSqr = self.errSqr.tolist()
        print( 'append', entry )
        if isinstance( entry, ufloat ):
            self.val.append( entry.val )
            self.errSqr.append( entry.errSqr )
            if not self.unit == entry.unit:
                raise Exception('unmatching units: %s and %s' % ( self.unit, entry.unit ) )
            self.val = array(self.val)
            self.errSqr = array(self.errSqr)
            self.summary()
            return
        raise Exception('cannot add type: %s' % ( type(u) ) )

    def tolist(self):
        errSqr = array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        return [ ufloat( val, errSqr = errSqr_, unit=self.unit ) for val, errSqr_ in zip(self.val, errSqr) ]
    
    def __repr__(self):
        return '['+', '.join([ u.__repr__() for u in self.tolist() ])+']'
    
    #def __mul__( self, u ):
        #try:
            #return [

#def sqrt(x):
    #return x**.5

if __name__=='__main__':
    ufloat.add_conversion('ADU', 'eV', ufloat(7,err=.1) )
    a = ufloat(5.123456, err = .055, unit='keV')
    a1 = ufloat(5.123456, err = 50, unit='eV')
    b = ufloat(2, unit = 'eV')
    print( 'a', a )
    print( 'b', b )
    print( 'abs(a)', abs(a) )
    print( 'a+b', a+b )
    print( 'a*2', a*2 )
    print( 'a/1000', a/1e3 )
    print( 'a+a', a+2*a )
    print( 'a*a', a*a )
    print( 'a**2', a**2 )
    print( 'a/(2a)', a/(2*a) )
    print( 'a/b', a/b )
    print( 'sqrt(a*a)', sqrt( a*a ) )
    print( 'sqrt(a**2)', sqrt( a**2 ) )
    #print( 'sqrt(a)', sqrt(a) )
    print( 'a-a', a-a )
    b = ufloat(10, err=.02, unit='ADU')
    print( 'b', b )
    print( '1/b', 1./b )
    print( 'b^{-1}', b**(-1) )
    print( 'a*b', a*b )
    #print( 2**a )
    print( 'a/2',a/2 )
    print( '1/a',1/a )
    print( 'a**2', a**2, (a**2).asunit('ADU') )
    print( '1ADU', ufloat(1,unit='ADU').asunit('eV') )
    print( '1ADU', ufloat(1,unit='ADU').asunit('e-').asunit('eV') )
    #print 'cos', cos( ufloat(1,err=.1) )
    print( ufloat(5,.01) )
    c = uarray( [1, 2, 3, 4, 5.5, 5], unit='keV')
    print( 'c', c )
    print( 'mean', c.mean() )
    print( 'std', c.std() )
    c.append( ufloat(4,err=.1,unit='keV') )
    c.append( ufloat(3.5,err=.1,unit='keV') )
    c.append( ufloat(3.25,err=.1,unit='keV') )
    print( 'c', c )
    print( 'mean', c.mean() )
    print( 'std', c.std() )
    print( '2*c', 2*c )
    print( 'c**2', c**2 )
    exit(0)
    
    
    a = uncertanty(2.5,errRel=.2)
    print( 'a', a )
    print( 'a+a', a+a )
    print( '2*a', 2*a )
    print( '-a',-a )
    print( '9.4*a', 9.5*a )
    print( 10*a )
    print( 100*a )
    print( a**2 )
    print( 2**a )
    print( a**.5 )
    print( sqrt(a) )
    print( 'a/2', a/2 )
    print( '1/a', 1/a )
    print( a/1000 )
    
    exit(0)
    X = norm.rvs(loc=2.,scale=4.,size=int(1e5))
    poisson_norm.verbose = True
    Y = poisson_norm.rvs(loc=3,scale=2.,g=1904, mu=.4, size=1000)
    print( 'pn.rvs', len(X) )
    print( 'norm.fit', norm.fit(X) )
    print( 'MLE', maxloglikelihood( lambda x, loc, scale: norm2(x, loc, scale), X, p0=(0.,1.), jac=None, bounds=((None,None),(0,None)) )['x'] )
    print( 'pn.fit2', poisson_norm.fit2(X) )
    #print 'pn.fit', poisson_norm.fit(X)
    
    print( 'pn.fit', poisson_norm.fit2(Y) )
    exit(0)
