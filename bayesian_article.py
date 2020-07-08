from __future__ import print_function
from Beamer import *
import os
import subprocess
from Timer import Timer
import sys
from termcolor import colored
from numpy import *

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
        E = r'{\mathbf E}'
        Q = r'{\mathbf Q}'
        xbar = r'{\bar x}'
        ybar = r'{\bar y}'
        sqrt2 = r'\sqrt{2}'
        log = r'{\rm log}'
        exp = r'{\rm exp}'
        erf = r'{\rm erf}'
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
        #doc.frame(
            #'joint PDF',
            #r'''
            #Given data $\{{x_1, x_2, \dots, x_n \}}$, the task is to find the mean $\mu$ and standard-deviations $\sigma$ of these data. First, assume that data has a common mean and a common variance. The principle of maximum entropy can then be applied under these constraints (using a flat "ignorance" prior) to choose the distribution
            #$$
            #f({X}|\mu,\sigma)
            #= \prod_i^n \frac{{1}}{{ \sqrt{{2\pi}}\sigma }}{exp}[ -\frac{{(x_i-\mu)^2 }}{{2\sigma^2}} ]
            #\\
            #= \frac{{1}}{{ (\sqrt{{2\pi}} \sigma )^n }}
            #{exp}[ -\frac{{1}}{{2\sigma^2}} \sum_i( x_i^2-2x_i+\mu^2 ) ]
            #$$
            #'''.format(**locals())
        #)
 
        doc.frame(
            'prior',
            r'''
            the posterior distribution is given by the Bayes rule
            $$
            f(\mu,\sigma|{X})
            = \frac{{ f( {X}|\mu,\sigma ) f(\mu,\sigma) }}{{ f({X}) }}
            \\
            = D_n f({X}|\mu,\sigma) f(\mu,\sigma)
            $$
            where $D_n$ is the normalization to keep the integral of the joint pdf to 1.\\
            the prior needs to transform as
            $$
            f(\mu,\sigma) = a f( \mu + b, a \sigma )
            $$
            so the prior needs to have the form
            $$
            f( \mu, \sigma ) =  \frac{{ \rm constant }}{{ \sigma }}
            $$
            '''.format(**locals())
        )
        doc.frame('a posteriori PDF',
            r'''
            the probability of the parameters given the data
            $$
            f(\mu,\sigma|{Y}) 
            = \frac{{ D_n }}{{ \sigma ^{{n+1}} }}
            {exp}[ -\frac{{1}}{{2\sigma^2}} \sum_i ( x_i - \mu )^2 ]
            \\
            = \frac{{ D_n }}{{ \sigma^{{n+1}} }}
            {exp}[ -\frac{{ (\mu - {xbar})^2 + C }}{{2 \sigma^2/n}} ]
            $$
            where
            $$
            {xbar} = \frac{{1}}{{n}}\sum_i x_i
            \\
            C_x = \frac{{1}}{{n}}\sum_i (x_i - {xbar})^2
            $$
            '''.format(**locals())
        )

        doc.frame('normalization',
            r'''
            The normalization is needed to ensure the probability denstity interpretation
            $$
            D_n
            = [ \int d\alpha d\sigma \frac{{1}}{{ \sigma^{{n+1}} }}
               {exp}( -\frac{{ \alpha^2 + C }}{{2 \sigma^2/n}} ) ]^{{-1}}
            \\
            = \frac{{ n^n \sqrt{{ C_x^{{n-1}}C_y^{{n-1}} }} }}{{ \pi 2^{{n-2}} }} \Gamma(\frac{{n-1}}{{2}})^{{-2}}
            $$
            '''.format(**locals())
        )

        doc.frame('joint maximum a posteriori estimator',
            r'''
            $$
            \hat{{\mu}}_x, \hat{{\sigma}}
            = {{\rm arg\ max}} f( \mu, \sigma|{X} )
            \\
            = {{\rm arg\ min}}[ -{{\rm log}} f( \mu, \sigma|{X} ) ]
            \\
            = {{\rm arg\ min}}[ (n+1){{\rm log}}\sigma
                + \frac{{ (\mu - \bar{{x}})^2 + C }}{{2 \sigma^2/n}} ]
            $$
            '''.format(**locals())
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
        #doc.frame(
            #'joint PDF',
            #r'''
            #Given data $\{{x_1, x_2, \dots, x_n \}}$ and $\{{y_1, y_2, \dots, y_n \}}$, the task is to find the means ${mux}$ and ${muy}$, variances $v_x={sigmax}^2$ and standard-deviations ${sigmax}$ and ${sigmay}$ of these data. First, assume that data has a common mean and a common variance. The principle of maximum entropy can then be applied under these constraints (using a flat "ignorance" prior) to choose the distribution
            #$$
            #f({X},{Y}|{mux},{muy},{sigmax},{sigmay})
            #= \prod_i^n \frac{{1}}{{ \sqrt{{2\pi}} {sigmax}{sigmay} }}{exp}[ -\frac{{x_i-{mux}}}{{2{sigmax}}} - \frac{{y_i-{muy} }}{{2{sigmay}}} ]
            #\\
            #= \frac{{1}}{{ (\sqrt{{2\pi}} {sigmax}{sigmay})^n }}{exp}[ -\frac{{1}}{{2{sigmax}^2}} \sum_i x_i-{mux} -\frac{{1}}{{2{sigmay}^2}} \sum_i y_i-{muy} ]
            #$$
            #'''.format(**locals())
        #)
 
        #doc.frame(
            #'prior',
            #r'''
            #the posterior distribution is given by the Bayes rule
            #$$
            #f({mux},{muy},{sigmax},{sigmay}|{X},{Y})
            #= \frac{{ f({X},{Y}|{mux},{muy},{sigmax},{sigmay}) f({mux},{muy},{sigmax},{sigmay}) }}{{ f({X},{Y}) }}
            #\\
            #= D_n f({X},{Y}|{mux},{muy},{sigmax},{sigmay}) f({mux},{muy},{sigmax},{sigmay})
            #$$
            #the prior needs to transform as
            #$$
            #f({mux},{muy},{sigmax},{sigmay}) = a_xa_y f({mux}+b_x,{muy}+b_y,a_x{sigmax},a_y{sigmay})
            #$$
            #so the prior needs to have the form
            #$$
            #f({mux},{muy},{sigmax},{sigmay}) =  \frac{{ \rm constant }}{{ {sigmax}{sigmay} }}
            #$$
            #'''.format(**locals())
        #)
        #doc.frame('a posteriori PDF',
            #r'''
            #$$
            #f({mux},{muy},{sigmax},{sigmay}|{X},{Y}) 
            #= \frac{{ D_n }}{{ ({sigmax}{sigmay})^{{n+1}} }}
            #\\ \quad\times {exp}[ -\frac{{1}}{{2{sigmax}^2}} \sum_i (x_i-{mux})^2 -\frac{{1}}{{2{sigmay}^2}} \sum_i( y_i-{muy})^2 ]
            #\\
            #= \frac{{ D_n }}{{ ({sigmax}{sigmay})^{{n+1}} }}{exp}[ -\frac{{ ({mux} - {xbar})^2 + C_x }}{{2{sigmax}^2/n}} -\frac{{ ({muy} - {ybar})^2 + C_y }}{{2{sigmay}^2/n}} ]
            #$$
            #where
            #$$
            #{xbar} = \frac{{1}}{{n}}\sum_i x_i
            #\\
            #C_x = \frac{{1}}{{n}}\sum_i (x_i - {xbar})^2
            #$$
            #'''.format(**locals())
        #)
        #doc.frame('normalization',
            #r'''
            #$$
            #D_n
            #= [ \int d\alpha_x d\sigma_x d\alpha_y d\sigma_y \frac{{1}}{{({sigmax}{sigmay})^{{n+1}} }}
               #{exp}( -\frac{{ \alpha_x^2 + C_x }}{{2{sigmax}^2/n}} -\frac{{ \alpha_y^2 + C_y }}{{2{sigmay}^2/n}} ) ]^{{-1}}
            #\\
            #= [ \int d\alpha_x d\sigma_x \frac{{1}}{{ {sigmax}^{{n+1}} }}
                #{exp}( -\frac{{ \alpha_x^2 + C_x }}{{2{sigmax}^2/n}} ) ]^{{-1}} [x y]^{{-1}}
            #\\
            #= \frac{{ n^n \sqrt{{ C_x^{{n-1}}C_y^{{n-1}} }} }}{{ \pi 2^{{n-2}} }} \Gamma(\frac{{n-1}}{{2}})^{{-2}}
            #$$
            #'''.format(**locals())
        #)
        #doc.frame('joint maximum a posteriori estimator',
            #r'''
            #$$
            #\hat{{\mu}}_x, \hat{{\mu}}_y, \hat{{\sigma}}_x, \hat{{\sigma}}_y 
            #= {{\rm arg\ max}} f({mux},{muy},{sigmax},{sigmay}|{X},{Y})
            #\\
            #= {{\rm arg\ min}}[ -{{\rm log}} f({mux},{muy},{sigmax},{sigmay}|{X},{Y}) ]
            #\\
            #= {{\rm arg\ min}}[ (n+1){{\rm log}}\sigma_x + (n+1){{\rm log}}\sigma_x 
                #\\ \quad+ \frac{{ (\mu_x - \bar{{x}})^2 + C_x }}{{2{sigmax}^2/n}} + \frac{{ (\mu_y - \bar{{y}})^2 + C_y }}{{2{sigmay}^2/n}} ]
            #$$
            #'''.format(**locals())
        #)

        #doc.frame('joint maximum a posteriori estimator',
            #r'''
            #$$
            #\frac{{ (\hat{{\mu}}_x - \bar{{x}})}}{{\hat{{\sigma}}_x^2/n}} = 0
            #\\
            #\frac{{n+1}}{{\hat{{\sigma}}_x }} + \frac{{ (\hat{{\mu}}_x - \bar{{x}})^2 + C_x }}{{2\hat{{\sigma}}_x^3/n}} = 0
            #$$
            #solving simulataneously
            #$$
            #\hat{{\mu}}_x = \bar{{x}}
            #\\
            #\hat{{\sigma}}_x^2 = \frac{{n}}{{n+1}}C_x
            #$$
            #'''.format(**locals())
        #)
        
        ################################################################# multinomial
        doc.frame('pixel integrated PDF',
            r'''
            The joint distribution function from the histogram is given by the multinomial distribution
            $$
            f(q_i,x_i|\mu,\sigma) = \frac{{ p_i^{{q_i}} }}{{q_i!}}
            $$
            $$
            f({Q},{X}|\mu,\sigma) = \frac{{ \Gamma(1+\sum_i q_i) }}{{ \prod_i \Gamma(q_i+1) }} \prod_i p_i^{{q_i}}
            $$
            where
            $$
            p_i = \frac{{1 }}{{ \sqrt{{2\pi}}\sigma }}{exp}[ - \frac{{ (x_i - \mu)^2 }}{{ 2\sigma^2 }}]
            $$
            '''.format(**locals())
        )
        doc.frame('pixel integrated PDF',
            r'''
            The maximization of the log probability
            $$
            {log} f({Q},{X}|\mu,\sigma) = {log}\Gamma(1+\sum_i q_i) - \sum_i {log} \Gamma(q_i+1) + \sum_i q_i {log} p_i
            $$
            where
            $$
            {log} p_i = - \frac{{1}}{{2}}{log}(2\pi) - {log}\sigma - \frac{{ (x_i - \mu)^2 }}{{ 2\sigma^2 }}
            $$
            the a posteriori
            $$
            {log} f(\mu,\sigma|{Q},{X}) = {log}D_n + {log} f({Q},{X}|\mu,\sigma) - {log}\sigma
            $$
            where the last term comes from the a priori distribution
            '''.format(**locals())
        )

        doc.frame('pixel integrated PDF',
            r'''
            all variable the terms gathered,
            $$
            {log} f(\mu,\sigma|{Q},{X}) = \sum_i q_i (- {log}\sigma - \frac{{ (x_i - \mu)^2 }}{{ 2\sigma^2 }}) - {log}\sigma
            $$
            
            $$
            {log} f(\mu,\sigma|{Q},{X}) = -\sum_i q_i \frac{{ (x_i - \mu)^2 }}{{ 2\sigma^2 }} - (Q+1){log}\sigma
            $$
            the $\mu$ variation
            $$
            \sum_i q_i \frac{{ (x_i - \mu) }}{{ \sigma^2 }} = 0
            $$
            so 
            $$
            \mu =  \frac{{ \sum_i q_i x_i }}{{ \sum_i q_i }}
            $$
            '''.format(**locals())
        )

        doc.frame('pixel integrated PDF',
            r'''
            the $\sigma$ variation
            $$
            \sum_i q_i \frac{{ (x_i - \mu)^2 }}{{ \sigma^3 }} - (Q+1)\frac{{1 }}{{\sigma }} = 0
            $$
            so 
            $$
            \sum_i q_i(x_i - \mu)^2 = (n+1) \sigma^2
            $$
            $$
            \sigma^2 = \frac{{\sum_i q_i x_i^2 - \mu^2\sum_i q_i }}{{1 + \sum_i q_i}}
            $$
            '''.format(**locals())
        )

        ################################################################# multinomial normal
        doc.frame('pixel integrated PDF',
            r'''
            The joint distribution function from the histogram is given by the multinomial distribution
            $$
            f(q_i,x_i|\mu,\sigma) = \frac{{ p_i^{{q_j}} }}{{q_j!}}
            $$
            $$
            f({Q},{X}|\mu,\sigma) = \Gamma(1+\sum_i q_i) \prod_i \frac{{ p_i^{{q_i}} }}{{ \Gamma(q_i+1) }}
            $$
            where
            $$
            p_i = \frac{{1 }}{{ \sqrt{{2\pi}}\sigma }}\int_0^1d\xi {exp}[ - \frac{{ (x_i + \xi - \mu)^2 }}{{ 2\sigma^2 }}]
            $$
            '''.format(**locals())
        )
        doc.frame('pixel integrated PDF',
            r'''
            The maximization of the log probability
            $$
            {log} f({Q},{X}|\mu,\sigma) = {log}\Gamma(1+\sum_i q_i) - \sum_i {log} \Gamma(q_i+1) + \sum_i q_i {log} p_i
            $$
            where
            $$
            {log} p_i = - \frac{{1}}{{2}}{log}(2\pi) - {log}\sigma + {log} \int_0^1d\xi {exp}[ - \frac{{ (x_i + \xi - \mu)^2 }}{{ 2\sigma^2 }}]
            $$
            the a posteriori
            $$
            {log} f(\mu,\sigma|{Q},{X}) = {log}D_n + {log} f({Q},{X}|\mu,\sigma) - {log}\sigma
            $$
            where the last term comes from the a priori distribution
            '''.format(**locals())
        )

        doc.frame('pixel integrated PDF',
            r'''
            the $\mu$ variation
            $$
            \sum_i q_i \frac{{\delta  p_i}}{{p_i}} = 0
            $$
            and
            $$
            \sum_i\frac{{ q_i }}{{ \sqrt{{2\pi}}\sigma^3p_i }}\int_0^1d\xi (x_i + \xi - \mu) {exp}[ - \frac{{ (x_i + \xi - \mu)^2 }}{{ 2\sigma^2 }}] = 0
            $$

            $$
            \sum_i\frac{{ q_i }}{{ p_i }}\{{ (x_i - \mu)p_i + \int_0^1d\xi \xi {exp}[ - \frac{{ (x_i + \xi - \mu)^2 }}{{ 2\sigma^2 }}] \}}= 0
            $$
            
            $$
            \sum_i q_i(x_i - \mu) = - \sum_i\frac{{ q_i }}{{ p_i }} m_1(x_i-\mu)
            $$
            '''.format(**locals())
        )

        doc.frame('pixel integrated PDF',
            r'''
            the $\sigma$ variation
            $$
            -\frac{{1}}{{\sigma}} + \sum_i q_i \frac{{\delta  p_i}}{{p_i}} = 0
            $$
            and
            $$
            \delta p_i = -\frac{{p_i}}{{\sigma}} + \frac{{1}}{{\sigma}}\int_0^1d\xi \frac{{ (x_i + \xi - \mu)^2 }}{{ \sigma^3 }} {exp}[ - \frac{{ (x_i + \xi - \mu)^2 }}{{ 2\sigma^2 }}]
            $$

            $$
            \delta p_i = -\frac{{p_i}}{{\sigma}} + \frac{{1}}{{\sigma^4}}\int_0^1d\xi ((x_i-\mu)^2 + \xi^2 - 2(x_i-\mu)\xi) {exp}[ - \frac{{ (x_i + \xi - \mu)^2 }}{{ 2\sigma^2 }}]
            $$

            $$
            \delta p_i = -\frac{{p_i}}{{\sigma}} + \frac{{1}}{{\sigma^4}}[(x_i-\mu)^2\sigma p_i - 2(x_i-\mu) m_1(x_i-\mu) + m_2(x_i-\mu) ]
            $$
            
            '''.format(**locals())
        )
        doc.frame('pixel integrated PDF',
            r'''
            the $\sigma$ variation
            $$
            0 = -\frac{{1}}{{\sigma}} -\frac{{1}}{{\sigma}} \sum_i q_i + \sum_i \frac{{q_i}}{{\sigma^4p_i}}[(x_i-\mu)^2\sigma p_i 
            \\ \quad - 2(x_i-\mu) m_1(x_i-\mu) + m_2(x_i-\mu) ]
            $$

            $$
            (1+Q)\sigma^2 = \sum_i q_i(x_i-\mu)^2 - \sum_i \frac{{q_i}}{{\sigma p_i}}2(x_i-\mu) m_1(x_i-\mu) 
            \\ \quad+ \sum_i \frac{{q_i}}{{\sigma p_i}}m_2(x_i-\mu)
            $$
            
            '''.format(**locals())
        )
        ################################################################# multinomial normal
        doc.frame('pixel integrated PDF',
            r'''
            $$
            f({X}|\mu,\sigma)
            \\= \prod_i^n \int_{{0}}^{{1}}d\xi \frac{{ 1 }}{{ \sqrt{{2\pi}} \sigma }}
            {exp}[ -\frac{{(x_i+\xi-\mu)^2}}{{2\sigma^2}} ]
            \\
            = \frac{{ 1 }}{{ (2\pi)^{{n/2}} \sigma^n }}
            \int_{{0}}^{{1}}d\xi {exp}[ -\frac{{1}}{{2\sigma^2}}\sum_i^n (x_i+\xi-\mu)^2 ]
            \\
            = \frac{{ 1 }}{{ (2\pi)^{{n/2}} \sigma^n }}
            \int_{{0}}^{{1}}d\xi {exp}[ -\frac{{(\bar{{x}}+\xi-\mu)^2 + C}}{{2\sigma^2/n}} ]
            $$
            '''.format(**locals())
        )


        doc.frame('a posteriori PDF',
            r'''
            $$
            \sigma^2 - 2\frac{{n}}{{ n+1 }}{exp}[ -\frac{{C}}{{2\sigma^2/n}} ]\int_{{ 0 }}^{{ \frac{{1}}{{2}} }}d\xi [\xi^2+C] {exp}[ -\frac{{\xi^2}}{{2\sigma^2/n}} ] = 0
            $$
            $$
            \sigma^2 - 2\sqrt{{2}}\frac{{n \sigma}}{{ n+1\sqrt{{n}} }} {exp}[ -\frac{{C}}{{2\sigma^2/n}} ]\\ \quad
            \times \int_0^{{ \frac{{ \sqrt{{n/2}} }}{{ 2\sigma }} }}d\phi [ 2\sigma^2\phi^2/n + C] {exp}( -\phi^2 ) = 0
            $$
            the continuos limit show be recovered $\sigma^2/n \rightarrow 0$ the integral extends to infinity,
            $$
            \sigma^2 - 2\sqrt{{2}}\frac{{n \sigma}}{{ n+1\sqrt{{n}} }} {exp}[ -\frac{{C}}{{2\sigma^2/n}} ]\\ \quad
            \int_0^\infty d\phi [ 2\sigma^2\phi^2/n + C] {exp}( -\phi^2 ) = 0
            $$            
            '''.format(**locals())
        )

        doc.frame('the variation for $\sigma$',
            r'''
            $$
            \sigma^2 - \frac{{n}}{{ n+1 }} {exp}[ -\frac{{C}}{{2\sigma^2/n}} ] [\frac{{\sigma^2}}{{n}}g_2( \frac{{\sigma^2}}{{n}} ) + Cg_0(\frac{{\sigma^2}}{{n}})] = 0
            $$
            $$
            g_0(x^2) = \int_{{ -\frac{{1}}{{2}} }}^{{ \frac{{1}}{{2}} }}d\xi{exp}[ -\frac{{\xi^2}}{{2x^2}} ] < 1
            $$
            $$
            x g_2(x) = \int_{{ -\frac{{1}}{{2}} }}^{{ \frac{{1}}{{2}} }}d\xi \xi^2 {exp}[ -\frac{{\xi^2}}{{2x}} ], \quad g_2(x) < 1
            $$
            '''.format(**locals())
        )

        doc.frame('joint maximum a posteriori estimator',
            r'''
            wich provides the two sets of equations to be solve simultaneously
            $$
            \sum_i\frac{{ {fN}(x_i - \mu, \sigma) - {fN}(x_i+1 - \mu, \sigma) }}{{ [{CN}(\xi - \mu, \sigma) ]_{{x_i}}^{{x_i+1}} }} = 0
            $$
            and
            $$
            + \sum_i\frac{{ x_i{fN}(x_i - \mu,\sigma) }}{{ [{CN}(\xi - \mu, \sigma) ]_{{x_i}}^{{x_i+1}} }} 
            \\ \quad - \sum_i\frac{{ x_i{fN}(x_i+1 - \mu, \sigma) }}{{ [{CN}(\xi - \mu, \sigma) ]_{{x_i}}^{{x_i+1}} }} 
            \\ \quad - \sum_i\frac{{ {fN}(x_i+1 - \mu, \sigma) }}{{ [{CN}(\xi - \mu, \sigma) ]_{{x_i}}^{{x_i+1}} }} = \frac{{1}}{{ 2 }} 
            $$
            '''.format(**locals())
        )

        ##############################################################
        doc.frame('2d pixel integrated PDF',
            r'''
            $$
            f({X},{Y},{E}|{mux},{muy},{sigmax},{sigmay})
            \\= \prod_i^n \int_{{x_i}}^{{x_i+1}}dx \int_{{y_i}}^{{y_i+1}}dy e_i {fN}(x-{mux}, {sigmax}) {fN}(y-{muy}, {sigmay})
            %\\= \prod_i^n \int_{{x_i}}^{{x_i+1}}dx \int_{{y_i}}^{{y_i+1}}dy \frac{{ e_i }}{{ \sqrt{{2\pi}} {sigmax}{sigmay} }}
            %{exp}[ -\frac{{x-{mux}}}{{2{sigmax}^2}} - \frac{{y-{muy} }}{{2{sigmay}^2}} ]
            %\\
            %= \prod_i^n \frac{{e_i}}{{ 4 }} 
            %[ {erf}(\frac{{ x_i+1-{mux} }}{{ {sigmax} {sqrt2} }}) - {erf}(\frac{{ x_i-{mux} }}{{ {sigmax} {sqrt2} }}) ]
            %[ x \rightarrow y ]
            \\
            = \prod_i^n e_i [{CN}(\xi - {mux}, {sigmax})]_{{x_i}}^{{x_i+1}} [{CN}(\xi - {muy}, {sigmay})]_{{y_i}}^{{y_i+1}}
            $$
            '''.format(**locals())
        )
        doc.frame('a posteriori PDF',
            r'''
            following the same steps as done for the joint distribution
            $$
            f({mux},{muy},{sigmax},{sigmay}|{X},{Y},{E}) 
            = \frac{{ D_n }}{{ {sigmax}{sigmay} }} \prod_i^n e_i 
            [{CN}(\xi - {mux}, {sigmax})]_{{x_i}}^{{x_i+1}}
            [{CN}(\xi - {muy}, {sigmay})]_{{y_i}}^{{y_i+1}}
            $$
            where
            $$
            D_n = [\int d\alpha_xd\sigma_xd\alpha_yd\sigma_y \frac{{ 1 }}{{ {sigmax}{sigmay} }} \prod_i^n e_i {I}(\alpha_x, {sigmax}) {I}(\alpha_y, {sigmay})]^{{-1}}
            \\
            = [\prod_i^n e_i \int d\alpha_xd\sigma_xd\alpha_yd\sigma_y \frac{{ 1 }}{{ {sigmax}{sigmay} }} {I}(\alpha_x, {sigmax}) {I}(\alpha_y, {sigmay})]^{{-1}}
            $$
            '''.format(**locals())
        )

        doc.frame('a posteriori PDF',
            r'''
            following the same steps as done for the joint distribution
            $$
            {log} f({mux},{muy},{sigmax},{sigmay}|{X},{Y},{E}) 
            = {log} D_n - {log}{sigmax} - {log}{sigmay} 
            \\ \quad + \sum_i^n {log} e_i + \sum_i^n {log}[{CN}(\xi - {mux}, {sigmax})]_{{x_i}}^{{x_i+1}} 
            + \sum_i^n {log} [{CN}(\xi - {muy}, {sigmay})]_{{y_i}}^{{y_i+1}}
            $$
            '''.format(**locals())
        )
        doc.frame('joint maximum a posteriori estimator',
            r'''
            wich provides the two sets of equations to be solve simultaneously
            $$
            \sum_i\frac{{ {fN}(x_i - {mux}, {sigmax}) - {fN}(x_i+1 - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }} = 0
            $$
            and
            $$
            -\frac{{1}}{{{sigmax}}} 
            + \frac{{ 2 }}{{ {sigmax} }} \sum_i\frac{{ (x_i-{mux}){fN}(x_i - {mux},{sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }} 
            \\ \quad - \frac{{ 2 }}{{ {sigmax} }} \sum_i\frac{{ (x_i+1-{mux}){fN}(x_i+1 - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }} = 0
            $$
            '''.format(**locals())
        )

        doc.frame('joint maximum a posteriori estimator',
            r'''
            the constraint of ${mux}$ applied to ${sigmax}$ yields
            $$
            \sum_i\frac{{ {fN}(x_i+1 - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }}
            = \sum_i\frac{{ {fN}(x_i - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }} 
            $$
            and
            $$
            \sum_i\frac{{ x_i {fN}(x_i - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }} 
            \\ \quad - \sum_i\frac{{ x_i{fN}(x_i+1 - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }}
            \\ - \sum_i\frac{{ {fN}(x_i - {mux}, {sigmax}) }}{{ [{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}} }} 
            = \frac{{1}}{{2}}
            $$
            '''.format(**locals())
        )

        doc.frame('joint maximum a posteriori estimator',
            r'''
            if each term is equal
            $$
            {fN}(x_i+1 - {mux}, {sigmax})
            = {fN}(x_i - {mux}, {sigmax})
            $$
            and
            $$
            x_i {fN}(x_i - {mux}, {sigmax}) - x_i{fN}(x_i+1 - {mux}, {sigmax}) - {fN}(x_i - {mux}, {sigmax})
            \\ \quad = \frac{{[{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}}}}{{2N}}
            $$
            
            $$
            {fN}(x_i - {mux}, {sigmax}) = - \frac{{[{CN}(\xi - {mux}, {sigmax}) ]_{{x_i}}^{{x_i+1}}}}{{2N}}
            $$
            '''.format(**locals())
        )
