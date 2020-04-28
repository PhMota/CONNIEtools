import subprocess
from sympy import *
from sympy.core.singleton import S

init_printing()
class PDF:
    def writelines( self, s, mode = 'a' ):
        open( '%s.tex' % self.fname, mode ).writelines( s )
    
    def __init__(self, fname='calculations'):
        self.fname = fname
        preamble = [
            r'\documentclass{article}',
            r'\usepackage{amsmath}',
            r'\usepackage{breqn}',
            r'\delimitershortfall=-1pt',
            r'\def\operatorname{}',
            r'\begin{document}'
            ]
        self.writelines( preamble, mode = 'w' )
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.writelines( [r'\end{document}'] )
        self.pdflatex()
    
    def pdflatex( self ):
        #os.system( 'pdflatex %s.tex' % self.fname )
        subprocess.check_output(['pdflatex', '%s.tex' % self.fname])

    def section( self, text ):
        self.writelines( [r'\section{%s}' % text] )
    
    def math( self, m, label = '' ):
        print 'eq', label
        self.writelines( [r'\begin{equation}', latex(m), r'\label{%s}'%label, r'\end{equation}'] )

    def par( self, s ):
        self.writelines( [s] )

def openPDF( fname ): return PDF( fname )

class Distribution(Function):
    @classmethod
    def pdf(self, x, *args):
        raise Exception('generic distribution not implemented')
    
    @classmethod
    def average( cls, weight, x, *args ):
        return integrate( cls.pdf(x, *args) * weight, (x, -oo, oo) )

    @classmethod
    def moment( cls, n, *args ):
        x = Symbol('x', real=True)
        return cls.average( x**n, x, *args ).doit(deep=True)

    @classmethod
    def mean( cls, *args ):
        return cls.moment( 1, *args )

    @classmethod
    def var( cls, *args ):
        x = Symbol('x', real=True)
        return cls.average( (x - cls.mean(*args) )**2, x, *args )
    
class Norm(Distribution):
    @classmethod
    def pdf( cls, x, mu, sigma ):
        return 1/(sqrt(2*pi))/sigma * exp( -S.Half*(x-mu)**2/sigma**2 )

class Poisson(Function):
    @classmethod
    def pmf( cls, k, lamb ):
        return lamb**k*exp(-lamb)/factorial(k)
    
    @classmethod
    def average( cls, weight, k, lamb ):
        return Sum( cls.pmf(k, lamb) * weight, (k, 0, oo) )

class NormPoisson(Distribution):
    @classmethod
    def pdf( cls, x, lamb, g, mu, sigma ):
        k = Symbol('k', real=True)
        return Sum( Poisson.pmf( k, lamb )*Norm.pdf( x - g*k, mu, sigma ), (k, 0, oo ) )

    @classmethod
    def average( cls, weight, x, *args ):
        k = Symbol('k', real=True)
        return Sum( integrate( cls.pdf(x, *args).args[0] * weight, (x, -oo, oo) ), (k, 0, oo ) )

def equal( a, b ):
    return latex(a) + '=' + latex(b)

def section_continuous( document ):
    document.section('Convolution')
    document.par( 'Let us explore the properties fo teh Norm-Poisson convolution probability distribution function,' )
    x, mu, g = symbols( 'x, mu, g', real=True )
    sigma = Symbol('sigma', real=True, positive=True)
    lamb = Symbol('lambda', real=True, positive=True)
    k = Symbol('k', real=True)

    normpoisson = Sum( Poisson.pmf(k, lamb)*Norm.pdf( x, g*k, sigma ), (k, 0, oo) )
    document.math( normpoisson, label='normpoisson' )

    document.par("the normal distribution's first moment (mean) is")
    document.math( Norm.moment( 1, mu, sigma ), label='norm 1st moment' )
    document.par("the normal distribution's 2nd moment is")
    document.math( Norm.moment( 2, mu, sigma ), label='norm 2nd moment' )
    document.par("the normal distribution's 3rd moment is")
    document.math( Norm.moment( 3, mu, sigma ), label='norm 3rd moment' )
    document.par("the distribution's first moment (mean) is")
    m1_expr = NormPoisson.moment( 1, lamb, g, mu, sigma ).expand().doit(deep=True)
    m1 = Symbol('m_1')
    document.math( equal(m1, m1_expr), label='np 1st moment' )
    #document.par("the distribution's first moment (mean) is")
    #document.math( NormPoisson.mean( lamb, g, mu, sigma ), label='np 1st moment' )
    document.par("the distribution's 2nd moment is")
    m2_expr = NormPoisson.moment( 2, lamb, g, mu, sigma ).expand().doit(deep=True)
    m2 = Symbol('m_2')
    document.math( equal(m2, m2_expr), label='np 2nd moment' )
    document.par("the distribution's variance is")
    var = Symbol(r'{\rm var}')
    var_expr = simplify(m2_expr - m1_expr**2)
    document.math( equal( var, var_expr ), label='np var' )

    document.par(r"using the eqs \ref{np 1st moment} and \ref{np var}, one can solve the system to")
    solved = solve( [m1 - m1_expr, var - var_expr], g, lamb, dict=True )[0]
    document.math( equal(g, solved[g]), label='np solved g' )
    document.math( equal(lamb, solved[lamb]), label='np solved lamb' )

    #document.par("the distribution's 3nd moment (mean) is")
    #m3_expr = NormPoisson.moment( 3, lamb, g, mu, sigma ).expand().doit(deep=True)
    #m3 = Symbol('m_3')
    #document.math( equal(m3, m3_expr) )
    #document.par("the 3rd moment minus all the lower order combinations is")
    #document.math( ((x-mu)**3).expand() )
    #m3_bar = Symbol(r'\bar{m}_3')
    #document.math( equal( m3_bar, (m3 - 3*m1*m2 + 3*m1**3 - m1**3).expand() ) )
    #m3_bar_expr = (m3_expr - 3*m1_expr*m2_expr + 2*m1_expr**3 ).expand().doit(deep=True)
    #document.math( equal( m3_bar, m3_bar_expr ) )
    
PDF.section_continuous = section_continuous

def section_median( document ):
    x, mu, g = symbols( 'x, mu, g', real=True )
    sigma = Symbol('sigma', real=True, positive=True)
    lamb = Symbol('lambda', real=True, positive=True)
    k = Symbol('k', real=True)
    
    document.section('Median')
    document.par("the median is calculated with")
    M = Symbol('M')
    lhs = simplify( integrate( Norm.pdf( x, mu, sigma ), (x, -oo, M) ) ) 
    rhs = simplify( integrate( Norm.pdf( x, mu, sigma ), (x, M, oo) ) )
    median_eq = lhs - rhs
    document.math( equal( median_eq, 0 ) )
    document.math( equal( M, solve( median_eq, M )[0] ) )

    #M = Symbol('M')
    #lhs = simplify( integrate( NormPoisson.pdf( x, mu, sigma, lamb, g ), (x, -oo, M) ) ) 
    #rhs = simplify( integrate( NormPoisson.pdf( x, mu, sigma, lamb, g ), (x, M, oo) ) )
    #median_eq = simplify( lhs - rhs ).doit()
    #document.math( equal( median_eq, 0 ) )
    #document.math( equal( M, solve( median_eq, M ) ) )
    
PDF.section_median = section_median

class Sample:
    def __init__( self, N = None ):
        if N:
            self.N = N
        else: 
            self.N = Symbol('N', real = True )
    
    def average( self, x, i ):
        return Sum( x, (i, 0, self.N -1) )/self.N
    
    def mean( self, x, i ):
        return self.average( x[i], i )

    def var( self, x, i ):
        return self.average( x[i]**2, i ) - self.mean(x, i)**2

    def moment( self, x, n, i ):
        return self.average( x[i]**n, i )
        
    def err( self, expr, x, dx, j ):
        core = ( expr.diff(x[j]) * dx )**2
        s = Sum( core, ( j, 0, self.N-1 ) ).doit(deep=True)
        return refine(s, Q.is_true( j <= self.N-1 ) )

def subs_mean( expr, kernel, symb ):
    #print expr, type(expr), expr.func, expr.args
    if len(expr.args) == 0: return expr
    if not hasattr(expr, 'func'): 
        return expr
    args = []
    for arg in expr.args:
        if arg.func == Sum:
            #print 'found sum', arg, arg.args[0]
            #print 'kernel', arg.args[0], kernel
            if arg.args[0] == kernel:
                #print 'found kernel', arg.args[0]
                args.append( symb )
            else:
                args.append( subs_mean(arg, kernel, symb) )
        else:
            args.append( subs_mean(arg, kernel, symb) )
    return expr.func( *args )

def section_sample( document ):
    document.section('Computing the sample parameters')
    X = IndexedBase('X')
    i = Symbol('i', positive = True, real=True)
    j = Symbol('j', positive = True, real=True)
    
    sample = Sample()
    sigma = Symbol('sigma', real=True)
    mu = Symbol('mu', real=True)
    mu_err = Symbol(r'\mu_{\rm err}', real=True)
    mean = sample.mean(X, i)
    document.math( equal(mu, mean), label = 'mean' )
    document.math( equal(mu_err, simplify( sample.err( mean, X, sigma, j ) ).doit()), label = 'mean.diff' )

    m2_expr = sample.moment(X, 2, i)
    m2 = Symbol('m_2', real=True)
    m2_err = Symbol(r'{m_{2}}_{\rm err}', real=True)
    document.math( equal(m2, m2_expr), label = 'm2' )
    document.math( equal(m2_err, simplify( sample.err( m2_expr, X, sigma, j ) ).doit()), label = 'm2.diff' )

    var_expr = sample.var( X, i )
    var = Symbol(r'{\rm var}')
    document.math( equal(var, var_expr), label = 'var' )
    dvar = refine(var_expr.diff(X[j]).doit(), Q.is_true( j <= sample.N-1) )
    document.math( equal( Derivative(var, X[j]), dvar), label = 'var diff' )
    
    var_err = Symbol( var.name + r'_{\rm err}' )
    var_err_expr = Sum( (dvar*sigma)**2, (j, 0, sample.N-1) ).doit()
    #print dir(var_err_expr)
    #print var_err_expr.variables
    #print var_err_expr.factor()
    document.math( equal( var_err, var_err_expr), label = 'var_expr' )
    document.math( equal( var_err, var_err_expr.factor()), label = 'var factor' )
    
    document.math( equal( var_err, subs_mean(var_err_expr.expand(), X[i], mu).factor() ), label = 'var subs_mean' )

    var_err_expr = 4*sigma**2/sample.N * var
    document.math( equal( var_err, var_err_expr.expand()), label = 'var err expr' )

PDF.section_sample = section_sample

if __name__ == '__main__':
    with openPDF( 'calculations' ) as document:
        document.section('SymPy')
        document.par('All the calculations in this document were automatically performed using the SymPy package.')
        
        document.section_continuous()
        document.section_median()
        document.section_sample()

        
