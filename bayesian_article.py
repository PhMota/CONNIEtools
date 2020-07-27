from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
#from numpy import *

folder = 'bayesian_article'

if not os.path.exists(folder):
    os.mkdir(folder)

with Timer('presentation'):
    with openBeamer( folder, '' ) as doc:

    
        mux = r'\mu_x'
        muy = r'\mu_y'
        sigmax = r'\sigma_x'
        sigmay = r'\sigma_y'
        X = r'{\mathbf X}'
        Y = r'{\mathbf Y}'
        xbar = r'{\bar x}'
        ybar = r'{\bar y}'
        sqrt2 = r'\sqrt{2}'
        I = r'{\mathcal I}'
        fN = r'f_{\mathcal N}'
        CN = r'F_{\mathcal N}'
        
        doc.par(r'\def\vr{\tt\color{gray}}')
        doc.frame('intro',
            'following the computation performed at "A Bayesian perspective on estimating mean, variance, and standard-deviation from data"\\ by Travis E. Oliphant',
            r'''\\https://scholarsarchive.byu.edu/cgi/\\
                viewcontent.cgi?article=1277\&context=facpub'''
        )
        
        ################################################################## normal distribution
        f = E|'f'
        X = E|'X'
        x = E|'x'
        Y = E|'Y'
        y = E|'y'
        i = E|'i'
        n = E|'n'
        C = E|'C'
        doc.frame(
            'joint PDF',
            r'''
            Given data $\{{x_1, x_2, \dots, x_n \}}$, the task is to find the mean $\mu$ and standard-deviations $\sigma$ of these data. First, assume that data has a common mean and a common variance. The principle of maximum entropy can then be applied under these constraints (using a flat "ignorance" prior) to choose the distribution''',
            Math(
                f( bf(X)|mu, sigma ) == (
                    Prod[i] * frac(1, sqrt(2*pi)*sigma ) * exp( -frac( (x[i] - mu)**2, 2*sigma**2) ),
                    frac(1, (2*pi)**(n/2) * sigma**(n) ) * exp( -frac( 1, 2*sigma**2) * Sum[i]*(x[i] - mu)**2 ),
                    frac(1, (2*pi)**(n/2) * sigma**(n) ) * exp( -n*frac( (mu-bar(x))**2 + C, 2*sigma**2) ),
                    )
            ),
            'where',
            Math(
                ( bar(x) == 1/n* Sum[i]*x[i], C == 1/n* Sum[i] * (x[i] - bar(x))**2 )
            ),
        )
        print( (x[i]**2 - bar(x)**2), (x[i]**2 - bar(x)**2).__class__ )
        print( Sum[i] * (x[i]**2 - bar(x)**2), (Sum[i] * (x[i]**2 - bar(x)**2) ).__class__ )
        a = E|'a'
        b = E|'b'
        D = E|'D'
        doc.frame(
            'prior',
            r'the posterior distribution is given by the Bayes rule',
            Math( 
                f(mu,sigma|bf(X)) == (
                    frac( f(bf(X)|mu,sigma)*f(mu,sigma), f(bf(X)) ),
                    D[n] * f(bf(X)|mu,sigma) * f(mu, sigma),
                    )
            ), 
            'where $D_n$ is the normalization to keep the integral of the joint pdf to 1.',
            'the prior needs to transform as',
            Math(
                f(mu,sigma) == a*f( mu + b, a *sigma )
            ),
            'so the prior needs to have the form',
            Math(
                f( mu, sigma ) == frac( r'{\rm constant}', sigma )
            )
        )
        doc.frame('a posteriori PDF',
            'the probability of the parameters given the data',
            Math(
                f(mu,sigma|bf(X)) == frac( D[n], sigma**(n+1) ) * exp( -n*frac( (mu-bar(x))**2 + C, 2*sigma**2) )
            )
        )

        doc.frame('normalization',
            'The normalization is needed to ensure the probability denstity interpretation',
            Math(
                D[n] == ( ( Int[0,oo] * d(sigma) * frac( 1, sigma**(n+1) ) * Int[-oo,oo] * d(alpha)
                    * exp( -n*frac( alpha**2 + C, 2 *sigma**2 ) ) )**(-1),
                    sqrt( frac( n**n * C**(n-1), pi * 2**(n-2) ) ) * Gamma(frac(n-1,2))**(-1)
                    )
            )
        )

        doc.frame('joint maximum a posteriori estimator',
            'The parameter estimators given the data are found by maximizing the posteriori distribution',
            Math( 
                hat(mu) *','* hat(sigma) == (
                    (E|r'{\rm arg\ max}')* f(mu,sigma|bf(X)),
                    (E|r'{\rm arg\ min}')( -log( f(mu,sigma|bf(X)) ) ),
                    (E|r'{\rm arg\ min}')( (n+1)*log(sigma) + n*frac( (mu - bar(x))**2 + C**2, 2*sigma**2 ) ),
                    )
            )
        )

        doc.frame('joint maximum a posteriori estimator',
            r'''
            $$
            \frac{{ \hat{{\mu}} - \bar{{x}} }}{{ \hat{{\sigma}}^2/n }} = 0
            \\
            \frac{{n+1}}{{ \hat{{\sigma}} }} + \frac{{ (\hat{{\mu}} - \bar{{x}} )^2 + C }}{{ 2\hat{{\sigma}}^3/n }} = 0
            $$
            solving simulataneously
            $$
            \hat{{\mu}} = \bar{{x}}
            \\
            \hat{{\sigma}}^2 = \frac{{n}}{{n+1}} C
            $$
            '''.format(**locals())
        )
        
        
        ################################################################## normal distribution 2d
        doc.frame('two dimensions',
            'this result can be straightforwardly extended to two dimensions since the probability is the product for each dimension',
            Math(
                f(bf(X), bf(Y)|mu[x], mu[y], sigma[x], sigma[y]) == f(bf(X)|mu[x], sigma[x]) * f(bf(Y)|mu[y], sigma[y])
            )
        )
                    
        ################################################################# normal
        p = E|'p'
        q = E|'q'
        Q = E|'Q'
        _E = E|'E'
        j = E|'j'
        m = E|'m'
        doc.frame('pixel integrated PDF',
            'In the case of pixelated dimensions, one can group repeated results as',
            Math(
                f(bf(X)|mu,sigma) == (
                    Prod[i,n] * p[i],
                    Prod[j,m] * p[j]**q[j]
                    )
            ),
            'where $j$ is the number of pixels and $q_j$ is the count of that specific value in the sample',
            Math(
                f(bf(Q),bf(X)|mu,sigma) == Prod[j]*p[j]**q[i]
            )
        )
            
        doc.frame('pixel integrated PDF',
            'The normalization comes from',
            Math(
                f(bf(Q), bf(X)) ==  Int[-oo,oo]*d(sigma) * Int[0,oo]*d(mu) * Prod[j,m] * p[j]**q[j]
            ),
            'if the probability is a normal distribution',
            Math(
                f(bf(Q), bf(X)|mu,sigma) == (
                    frac(1, (sqrt(2*pi) * sigma)**Q ) * exp( -frac( 1, 2*sigma**2 ) * Sum[j]* q[j] *( x[j]-mu )**2 ),
                    frac(1, (sqrt(2*pi) * sigma)**Q ) * exp( -Q*frac( ( mu-ave(x)[q] )**2 + C[q], 2*sigma**2 ) )
                    )
            ),
            Math(
                ( Q == Sum[j]*q[j], 
                ave(x)[q] == frac( 1, Q )* Sum[j]*( q[j]*x[j] ),
                C[q] == frac( 1, Q )* Sum[j]*( q[j]*x[j]**2 ) - ave(x)[q]**2
                )
            ),    
        )
        doc.frame('pixel integrated PDF',
            Math(
                f(q[j], x[j]) == (
                    frac(1, (2*pi)**(q[j]/2) ) * Int[0,oo] * d(sigma) * frac(1, sigma**q[j]) * Int[-oo,oo]*d(mu) * exp( -q[j]*frac( (x[j] - mu)**2, 2*sigma**2 ) ),
                    frac(1, (2*pi)**(q[j]/2) ) * Int[0,oo] * d(sigma) * frac(1, sigma**q[j]) 
                        * Int[-oo,oo]*d(mu) * exp( -q[j]*frac( x[j]**2 - 2*mu*x[i] + mu**2, 2*sigma**2 ) ),
                    frac(1, (2*pi)**(q[j]/2) ) * Int[0,oo] * d(sigma) * frac(1, sigma**q[j]) 
                        * e**( -q[j]*frac( x[j]**2, 2*sigma**2 ) ) * Int[-oo,oo]*d(mu) * exp( -q[j]*frac( mu**2 - 2*mu*x[i], 2*sigma**2 ) ),
                    frac(sqrt(2*pi*q[j]), (2*pi)**(q[j]/2) ) * Int[0,oo] * d(sigma) * frac(1, sigma**(q[j]-1)) 
                        * e**( -q[j]*frac( x[j]**2, 2*sigma**2 ) ) * exp( q[j]*frac( x[j]**2, 2*sigma**2 ) ),
                    )
            ),
            'the integral does not converge, which means we cannot normalize the single pdf'
        )
        doc.frame('pixel integrated PDF',
            "the joint pdf's normalization",
            Math(
                f(bf(Q), bf(X)) == (
                    frac(1, (2*pi)**(Q/2) ) * Int[0,oo] * d(sigma) * frac(1, sigma**Q) * Int[-oo,oo]*d(mu) * exp( -Q*frac( ( mu-ave(x)[q] )**2 + C[q], 2*sigma**2 ) ),
                    frac(1, (2*pi)**(Q/2) ) * Int[0,oo] * d(sigma) * frac( exp( -Q*frac( C[q], 2*sigma**2 ) ), sigma**Q) * Int[-oo,oo]*d(mu) * exp( -Q*frac( mu**2, 2*sigma**2 ) ),
                    frac( 1, (2*pi)**(Q/2) * sqrt(Q) ) * Int[0,oo] * d(sigma) * frac( exp( -Q*frac( C[q], 2*sigma**2 ) ), sigma**(Q-1) ),
                    frac( 1, (2*pi)**(Q/2) * sqrt(Q) ) * Int[0,oo] * frac(1, xi**2)*d(xi) * xi**(m-1) * exp( -frac( Q*C[q], 2 )* xi**2 ),
                    -frac( 1, (2*pi)**(Q/2) * sqrt(Q*pi) ) * (Q*C[q])**(-(Q-3)/2) * 2**((Q-3)/2) * Gamma((Q-2)/2),
                    )
            )
        )
        
        doc.frame('pixel integrated PDF',
            r'The maximization of the log probability',
            #Math(
                #log*f( bf(Q), bf(X) | mu, sigma ) == log * Gamma( 1 + Sum[i]*q[i] ) - Sum[i]*log*Gamma(q[i]+1) - Sum[i]*q[i]*log*p[i]
            #),
            'where',
            Math( 
                log*p[i] == - frac(1,2)*log(2*pi) - log*sigma - frac( (x[i]-mu)**2, 2*sigma**2 )
            ),
            'the log of a posteriori distribution',
            Math(
                log( f( mu,sigma|bf(Q),bf(X) ) ) == log(D[n]) + log( f(bf(Q),bf(X)|mu,sigma) ) - log( f(mu,sigma) )
            )
        )

        doc.frame('pixel integrated PDF',
            'all variable the terms gathered',
            Math(
                log *f(mu,sigma|bf(Q), bf(X) ) == Sum[j] * q[j] * log(p[j]) - log * f(mu,sigma)
            ),
            'using the normal distribution',
            Math(
                p[j] == frac(1, sqrt(2*pi) * sigma ) * exp( -frac( (x[j] - mu)**2, 2*sigma**2 ) )
            ),
            Math(
                log*p[j] == -frac(1,2)*log(2*pi) - log*sigma - frac( (x[j] - mu)**2, 2*sigma**2 ) 
            ),
            'so',
            Math(
                log * f == (
                    - log( f(mu,sigma) ) - frac(Q,2)*log(2*pi) - Q*log*sigma - Sum[j]* q[j]*frac( (x[j] - mu)**2, 2*sigma**2 ),
                    - log( f(mu,sigma) ) - frac(Q,2)*log(2*pi) - Q*log*sigma - Q*frac( (ave(x)[q] - mu)**2 + C[q], 2*sigma**2 )
                    )
            ),
        )

        doc.frame('joint maximum a posteriori estimators',
            r'the variation in $\mu$',
            Math(
                hat(mu) == ave(x)[q]
            ),
            r'the variation in $\sigma$',
            Math(
                0 == - frac( Q, hat(sigma) ) + Q*frac( (ave(x)[q] - mu)**2 + C[q], hat(sigma)**3 )
            ),
            Math(
                hat(sigma)**2 == (ave(x)[q] - mu)**2 + C[q]
            ),
            r'after using the information of $\hat{{\mu}}$',
            Math(
                hat(sigma)**2 == C[q]
            ),
            r'which is the weighted average of the pixel distribution'
        )

        ################################################################# normal
        doc.frame('pixel integrated PDF',
            'all variable the terms gathered',
            Math(
                log *f(mu,sigma|bf(Q), bf(X) ) == Sum[j] * q[j] * log(p[j]) - log * f(mu,sigma)
            ),
            'using the normal distribution',
            Math(
                p[j] == (
                    Int[x[j],x[j]+1] * d(xi) * frac(1, sqrt(2*pi) * sigma ) * exp( -frac( (xi - mu)**2, 2*sigma**2 ) ),
                    Int[x[j],x[j]+1] * d(xi) * frac(1, sqrt(2*pi) * sigma ) * exp( -frac( (xi - mu)**2, 2*sigma**2 ) ),
                    )
            ),
            Math(
                log*p[j] == -frac(1,2)*log(2*pi) - log*sigma - frac( (x[j] - mu)**2, 2*sigma**2 ) 
            ),
            'so',
            Math(
                log * f == (
                    - log( f(mu,sigma) ) - frac(Q,2)*log(2*pi) - Q*log*sigma - Sum[j]* q[j]*frac( (x[j] - mu)**2, 2*sigma**2 ),
                    - log( f(mu,sigma) ) - frac(Q,2)*log(2*pi) - Q*log*sigma - Q*frac( (ave(x)[q] - mu)**2 + C[q], 2*sigma**2 )
                    )
            ),
        )

        ################################################################# noise
        doc.frame('useful identities',
            'useful identities',
            Math(
                Int[-oo, oo]*d(xi) * exp( -frac(xi**2, 2) ) == sqrt(2*pi)
            ),
            Math(
                Int[-oo, oo]*d(xi) * exp( -frac((xi - mu)**2, 2*sigma**2) ) == sqrt(2*pi)*sigma
            ),
        )
        ################################################################# noise
        Pr = E|r'{\rm Pr}'
        e = E|'e'
        z = E|'z'
        doc.frame('noise probability',
            'the white noise',
            Math(
                f(e[i]|mu, sigma) == frac( 1, sqrt(2*pi)*sigma ) * exp( -frac( e[i]**2, 2*sigma**2 ) )
            ),
            r'the probability is for $e_i>\mu$',
            Math(
                Pr(e[i]) == (
                    Int[e[i],oo] * d(xi) * f(xi|mu, sigma),
                    Int[e[i],oo] * d(xi) * frac( 1, sqrt(2*pi)*sigma ) * exp( -frac( e[i]**2, 2*sigma**2 ) ),
                    frac(1,2)*(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ),
                    )
            ),
        )
        doc.frame('noise probability',
            r'z-score for $e_i>\mu$',
            Math(
                Int[z(e[i]),oo]* d(xi) * f(xi|0, 1) == Pr(e[i])
            ),            
            Math(
                frac(1,2)*(1 - erf( frac( z(e[i]), sqrt(2) ) ) ) == Pr(e[i])
            ),            
            Math(
                erf( frac( z(e[i]), sqrt(2) ) ) == 1 - 2*Pr(e[i])
            ),            
            Math(
                z(e[i]) == sqrt(2) * (erf**(-1))( 1 - 2*Pr(e[i]) )
            ),
            'applied to the single probability',
            Math(
                z(e[i]) == (
                    #sqrt(2) * frac( e[i] - mu, sqrt(2)*sigma ),
                    frac( e[i] - mu, sigma ),
                    )
            ),
            'which measures the deviation in units of standard deviation'
        )
        
        doc.frame('noise probability',
            'the joint probability is obtained by multiplying the single probabilities',
            Math(
                Pr(bf(_E)) == (
                    Prod[i] * frac(1,2)*(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ),
                    frac(1,2**n) * Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ),
                    )
            ),
            'and cannot be obtained directly from the joint distribution is',
            Math(
                f(bf(_E)|mu, sigma) == (
                    Prod[i] * frac( 1, sqrt(2*pi)*sigma ) * exp( -frac( (e[i]-mu)**2, 2*sigma**2 ) ),
                    frac( 1, (2*pi)**(n/2)*sigma**n ) * exp( -frac( 1, 2*sigma**2 )* Sum[i]*(e[i]-mu)**2 ),
                    frac( 1, (2*pi)**(n/2)*sigma**n ) * exp( -frac( (bar(e)-mu)**2 + C[e], 2*sigma**2 ) ),
                    )
            ),
        )

        doc.frame('noise probability',
            'the joint z-score can be found',
            Math(
                frac(1,2**n) * (1 - erf( frac( z(bf(_E)), sqrt(2) ) ) )**n == frac(1,2**n) * Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) )
            ),
            Math(
                z(bf(_E)) == sqrt(2)*(erf**(-1))( 1 - ( Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ) )**(1/n) )
            ),
            'for all values equal',
            Math(
                z(bf(_E)) == (
                    sqrt(2)*(erf**(-1))( 1 - (1 - erf( frac( e - mu, sqrt(2)*sigma ) ) ) ),
                    frac( e - mu, sigma )
                    )
            ),            
        )
            
        
        doc.frame('noise probability',
            'on the other hand, the z-score can come directly from the joint distribution',
            Math(
                Pr( z, n | 0, 1 ) == (
                    Int[z, oo] * d(xi) * frac( 1, (2*pi)**(n/2) ) * exp( -n*frac( xi**2, 2 ) ),
                    frac( 1, 2**((n+1)/2)*pi**((n-1)/2) )*(1 - erf( frac( sqrt(n)*z, sqrt(2) ) ) )
                    )
            ),
            'therefore, we can calculate the joint z-score',
            Math(
                frac( 1, 2**((n+1)/2)*pi**((n-1)/2) )*(1 - erf( frac( sqrt(n)*z, sqrt(2) ) ) ) == frac(1,2**n) * Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) )
            ),
            Math(
                z == sqrt(frac(2,n)) * (erf**(-1))( 1 - frac(2**((n+1)/2)*pi**((n-1)/2),2**n) * Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ) )
            ),
        )
                
        doc.frame('noise probability',
            'for $n=1$ this results reduces to the previous single z-score',
            Math(
                z == ( 
                    sqrt(frac(2,n)) * (erf**(-1))( 1 - frac(2**((n+1)/2)*pi**((n-1)/2),2**n) * Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ) ),
                    sqrt(2) * (erf**(-1))( 1 - (1 - erf( frac( e - mu, sqrt(2)*sigma ) ) ) ),
                    frac( e - mu, sigma ),
                 )
            ),
            'for $n$ equal values',
            Math(
                z == ( 
                    sqrt(frac(2,n)) * (erf**(-1))( 1 - frac(2**((n+1)/2)*pi**((n-1)/2),2**n) * Prod[i] *(1 - erf( frac( e[i] - mu, sqrt(2)*sigma ) ) ) ),
                    sqrt(frac(2,n)) * (erf**(-1))( 1 - frac(2**((n+1)/2)*pi**((n-1)/2),2**n) * (1 - erf( frac( e - mu, sqrt(2)*sigma ) ) )**n ),
                 )
            ),
        )
        
