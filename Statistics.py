# coding: utf-8

import numpy as np
import scipy.stats
import math

def stdCovMed( x, y ):
    X = x - np.median(x)
    Y = y - np.median(y)
    return np.mean(X*Y)
    
def madCovSqr( data0, data1 ):
    return np.median( np.abs( (data0 - np.median(data0))*((data1 - np.median(data1))) ) )

def normFit( data, loc=None, scale=None ):
    return scipy.stats.norm.fit( data, floc=loc, fscale=scale )

def tiles( data, p ):
    m = np.median(data)
    if p - 1 == 0: return [m]
    return tiles(data[data<m], p-1 ) + [m] + tiles(data[data>m], p-1 )

def quartiles(x):
    return tiles(x,2)

def IQR( data, scale=1.34897 ):
    return scipy.stats.iqr( data )/scale

def MAD( data, axis = None, scale=1.4826 ):
    if axis == 0:
        median = np.median( data, axis = 0 )[None,:]
    elif axis == 1:
        median = np.median( data, axis = 1 )[:,None]
    else:
        median = np.median( data )
    return scale * float(np.median( np.abs( data - median ), axis = axis ))

def mean_outliers( data, axis=None, n=2 ):
    med = np.median(data,axis)
    MAD_ = MAD(data,axis)
    return np.mean( data[ np.abs(data-med)<n*MAD_ ] )

norm = scipy.stats.norm
poisson = scipy.stats.poisson

class poisson_norm:
    verbose = False
    loc_bounds = (None,None)
    scale_bounds = (0,None)
    mu_bounds = (0,None)
    
    @classmethod
    def pdf2(cls, x, loc = 0, scale = 1, g = 1, mu = 1, tol=1e-8 ):
        result = np.zeros_like(x)
        mu = abs(mu)
        scale = abs(scale)
        g = abs(g)
        #loc = loc + g*mu**2
        
        k = 0
        poisson_sum = 0
        if cls.verbose: print 'poisson_gaussian_pdf call with', loc, scale, mu,
        poisson_weights = []
        while True:
            poisson_weights.append( poisson.pmf( k, mu=mu ) )
            poisson_sum += poisson_weights[-1]
            if poisson_weights[-1]/poisson_sum < tol: break
        res = np.sum([ w*norm.pdf( x, loc = loc + g*k*mu, scale = scale ) for k, w in enumerate(poisson_weights) ])
        #while True:
            #poisson_weight = poisson.pmf( k, mu=mu )
            #poisson_sum += poisson_weight
            #result += poisson_weight * norm.pdf( x, loc = loc + g*k*mu, scale = scale )
            #if poisson_weight/poisson_sum < tol: break
            #k += 1
        if cls.verbose: print 'took', len(poisson_weights), 'iteractions'
        return result

    @classmethod
    def pdf(cls, x, loc = 0, scale = 1, g = 1, mu = 1, tol=1e-8, fix_loc=False ):
        result = np.zeros_like(x)
        if fix_loc:
            loc = loc + g*mu
        
        k = 0
        poisson_sum = 0
        if cls.verbose: print 'poisson_gaussian_pdf call with', loc, scale, mu,
        while True:
            poisson_weight = poisson.pmf( k, mu=mu )
            poisson_sum += poisson_weight
            result += poisson_weight * norm.pdf( x, loc = loc + g*k, scale = scale )
            if poisson_weight/poisson_sum < tol: break
            k += 1
        if cls.verbose: print 'took', k, 'iteractions'
        return result

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
    def fit( cls, X, loc=0, scale=1., g=1, mu=1., floc=None, fscale=None, fg=None, fmu=None ):
        return tuple(maxloglikelihood( cls.pdf, X, p0=(loc,scale,g,mu), jac=False, bounds=(cls.loc_bounds,cls.scale_bounds,(1,None),cls.mu_bounds) )['x'])

    @classmethod
    def fit2( cls, X, loc=0, scale=1., g=1, gmu=1., floc=None, fscale=None, fg=None, fmu=None ):
        return tuple(maxloglikelihood( cls.pdf2, X, p0=(loc,scale,g,gmu), jac=False, bounds=(cls.loc_bounds,cls.scale_bounds,(1,None),(0,None)) )['x'])

def monteCarlo_rvs( pdf, a, b, N, normed=False, verbose=False ):
    if not normed:
        pdf_max = -scipy.optimize.fmin( lambda x: -pdf(x), 0, full_output=True )[1]
        if verbose: print 'pdf_max', pdf_max
        normedpdf = lambda x: pdf(x)/pdf_max
    else:
        normedpdf = pdf
    results = []
    if verbose: print '(random',
    while 1:
        x = (b - a) * np.random.random_sample( N ) + a
        tests = np.random.random_sample( N )
        y = normedpdf(x)
        mask = tests < y
        results = np.concatenate( (results, x[mask]) )
        if len(results) >= N: break
        if verbose: print '%.2f'%( float(len(results))/N ),
    if verbose: print ')'
    return results[:N]

def negloglikelihood( pdf, params, X ):
    res = -np.sum( np.log( pdf( X, *params )) )
    #print 'params', params, res
    return res

def maxloglikelihood( pdf, X, p0, jac=False, bounds=None ):
    return scipy.optimize.minimize( lambda p: negloglikelihood( pdf, p, X ), p0, jac=jac, bounds=bounds )

def norm2( x, loc, scale ):
    res = norm.pdf( x, loc, scale )
    #print 'loc,scale', loc, scale, res
    return res

def chisquare( y, f, yerr=None, ddof=1 ):
    if yerr is None:
        yerr = np.sqrt(y)
    return np.sum( ((y-f)/yerr)**2 )/(len(y)-ddof)

def errSqr_mean( x ):
    return np.var(x)/len(x)

def errSqr_std(x):
    return np.mean( (x-np.mean(x))**4 )/np.var(x)/len(x)

def errSqr_median(x):
    return MAD(x)**2/len(x)

def errSqr_MAD(x):
    return MAD(x)**2/len(x)
    
def sign(x):
    return -1 if x<0 else ( 1 if x>0 else 0 )

def pprint( msg, f ):
    print msg
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
            val = np.array([ entry.val for entry in self.value ])
            w = np.array([ 1./entry.errSqr for entry in self.value ])
            s = np.sum(w)
            return dfloat( ufloat( np.sum(w*val)/s, errSqr=1./s ), unit=self.unit )
        
        return dfloat( ufloat( np.mean(self.value), errSqr=errSqr_mean(self.value) ), unit=self.unit )

    def median(self):
        return dfloat( ufloat( np.median(self.value), errSqr=errSqr_median(self.value) ), unit=self.unit )

    def std(self):
        if isinstance(self.value[0], ufloat):
            val = np.array([ entry.val for entry in self.value ])
            w = np.array([ 1./entry.errSqr for entry in self.value ])
            s = np.sum(w)
            return dfloat( ufloat( np.sqrt( (np.sum(w*val**2) - np.sum(w*val)**2)/s ), errSqr=1./s ), unit=self.unit )
        
        return dfloat( ufloat( np.std(self.value), errSqr=errSqr_std(self.value) ), unit=self.unit )

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
    return np.around( x, -sig )
    
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
            raise Exception('unit type not recognized: %s' % (type(unit)) )

    def tuple(self):
        if isinstance(self.value, np.ndarray):
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
        if isinstance(self.value, np.ndarray):
            return '%s%s'%(self.value, self.unit_str())
        return '%s%s' % self.tuple()

    @classmethod
    def add_conversion(cls, unit1, unit2, value ):
        print 'add', unit1, unit2, value
        
        N = len(cls.__units__)
        if not unit1 in cls.__units__:
            cls.__units__.append( unit1 )
        if not unit2 in cls.__units__:
            cls.__units__.append( unit2 )
            
        __conversion_matrix__ = np.zeros(( len(cls.__units__), len(cls.__units__) )).tolist()
        for i1 in range(N):
            for i2 in range(N):
                __conversion_matrix__[i1][i2] = cls.__conversion_matrix__[i1][i2]
        
        i1 = cls.__units__.index(unit1)
        i2 = cls.__units__.index(unit2)
        __conversion_matrix__[i1][i1] = 1.
        __conversion_matrix__[i2][i2] = 1.
        __conversion_matrix__[i1][i2] = value
        __conversion_matrix__[i2][i1] = 1./value
        for i, u in enumerate(cls.__units__):
            if u != unit1 and u != unit2:
                i1 = cls.__units__.index(unit1)
                i2 = cls.__units__.index(unit2)
                __conversion_matrix__[i1][i] = __conversion_matrix__[i2][i] * value
                __conversion_matrix__[i2][i] = __conversion_matrix__[i1][i] / value
                __conversion_matrix__[i][i1] = __conversion_matrix__[i][i2] / value
                __conversion_matrix__[i][i2] = __conversion_matrix__[i][i1] * value
        cls.__conversion_matrix__ = __conversion_matrix__
    
    @classmethod
    def show_conversion_matrix(cls):
        print ' '*5 + ' '.join([ '%5s'%u for u in cls.__units__ ])
        for i1, u1 in enumerate(cls.__units__):
            print '%4s '%u1 + ' '.join([ '%.3f' % cls.__conversion_matrix__[i1][i2] for i2,u2 in enumerate(cls.__units__ )])
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

ufloat.add_conversion( 'e-', 'eV', 3.745 )

class uarray:
    '''
    array of uncertainties
    '''
    def __init__( self, value = None, err = None, errSqr = None, errMode = None, unit = None ):
        '''
        initiate with list of uncertainties
        '''
        if isinstace( value, np.ndarray ):
            self.unit = unit
            self.val = value
            if errMode == 'poisson':
                self.errSqr = np.array(self.val)
            elif errMode == 'stats':
                self.errSqr = np.array(np.var(self.val))/len(self.val)
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
                        raise Exception('unmatching units: %s and %s' % ( self.unit, entry.unit ) )
                else:
                    u = ufloat( entry, unit=unit )
                    self.val.append( u.val )
                    self.errSqr.append( u.errSqr )
        
        self.val = np.array(self.val)
        self.errSqr = np.array(self.errSqr)
        self.summary()
    
    def set_errSqr( self, errSqr ):
        self.errSqr = np.zeros_like(self.val)*errSqr
        
    def std_errSqr(self):
        #print 'len', np.var(self.val), len(self.val), (np.var(self.val)/len(self.val))**.5
        return np.var(self.val)/len(self.val)

    def weights(self):
        errSqr = np.array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        return 1./errSqr/self.weightSum()
    
    def weightSum(self):
        errSqr = np.array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        return np.sum( 1./errSqr )
    
    def mean(self):
        errSqr = np.array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        s = 1./self.weightSum()
        return ufloat( np.sum(self.val/errSqr)*s, errSqr=s, unit=self.unit )
    
    def mean_list(self):
        sum_ = ufloat(0,0,unit=self.unit)
        for entry in self.tolist():
            sum_ = sum_ + entry
        return sum_/len(self.val)
    
    def mean_stats(self):
        return ufloat( np.mean(self.val), errSqr = np.var(self.val)/len(self.val), unit = self.unit )
    
    def std(self):
        errSqr = np.array(self.errSqr)
        errSqr[ errSqr==0 ] = self.std_errSqr()
        s = 1./self.weightSum()
        return ufloat( (np.sum(self.val**2/errSqr)*s - (np.sum(self.val/errSqr)*s)**2)**.5, errSqr=s, unit=self.unit )

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
        return ufloat( np.std(self.val), errSqr = np.mean((self.val-np.mean(self.val))**4)/len(self.val)/np.var(self.val), unit = self.unit )
    
    def summary(self):
        print 'uarray', self
        print 'mean', self.mean(), self.mean_stats(), self.mean_list()
        print 'std', self.std(), self.std_stats(), self.std_list(), self.std_list2()
        
    def append(self, entry):
        self.val = self.val.tolist()
        self.errSqr = self.errSqr.tolist()
        print 'append', entry
        if isinstance( entry, ufloat ):
            self.val.append( entry.val )
            self.errSqr.append( entry.errSqr )
            if not self.unit == entry.unit:
                raise Exception('unmatching units: %s and %s' % ( self.unit, entry.unit ) )
            self.val = np.array(self.val)
            self.errSqr = np.array(self.errSqr)
            self.summary()
            return
        raise Exception('cannot add type: %s' % ( type(u) ) )

    def tolist(self):
        errSqr = np.array(self.errSqr)
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
    print 'a', a
    print 'b', b
    print 'abs(a)', abs(a)
    print 'a+b', a+b
    print 'a*2', a*2
    print 'a/1000', a/1e3
    print 'a+a', a+2*a
    print 'a*a', a*a
    print 'a**2', a**2
    print 'a/(2a)', a/(2*a)
    print 'a/b', a/b
    print 'sqrt(a*a)', np.sqrt( a*a )
    print 'sqrt(a**2)', np.sqrt( a**2 )
    #print 'sqrt(a)', np.sqrt(a)
    print 'a-a', a-a
    b = ufloat(10, err=.02, unit='ADU')
    print 'b', b
    print '1/b', 1./b
    print 'b^{-1}', b**(-1)
    print 'a*b', a*b
    #print 2**a
    print 'a/2',a/2
    print '1/a',1/a
    print 'a**2', a**2, (a**2).asunit('ADU')
    print '1ADU', ufloat(1,unit='ADU').asunit('eV')
    print '1ADU', ufloat(1,unit='ADU').asunit('e-').asunit('eV')
    #print 'cos', np.cos( ufloat(1,err=.1) )
    print ufloat(5,.01)
    c = uarray( [1, 2, 3, 4, 5.5, 5], unit='keV')
    print 'c', c
    print 'mean', c.mean()
    print 'std', c.std()
    c.append( ufloat(4,err=.1,unit='keV') )
    c.append( ufloat(3.5,err=.1,unit='keV') )
    c.append( ufloat(3.25,err=.1,unit='keV') )
    print 'c', c
    print 'mean', c.mean()
    print 'std', c.std()
    print '2*c', 2*c
    print 'c**2', c**2
    exit(0)
    
    
    a = uncertanty(2.5,errRel=.2)
    print 'a', a
    print 'a+a', a+a
    print '2*a', 2*a
    print '-a',-a
    print '9.4*a', 9.5*a
    print 10*a
    print 100*a
    print a**2
    print 2**a
    print a**.5
    print sqrt(a)
    print 'a/2', a/2
    print '1/a', 1/a
    print a/1000
    
    exit(0)
    X = norm.rvs(loc=2.,scale=4.,size=int(1e5))
    poisson_norm.verbose = True
    Y = poisson_norm.rvs(loc=3,scale=2.,g=1904, mu=.4, size=1000)
    print 'pn.rvs', len(X)
    print 'norm.fit', norm.fit(X)
    print 'MLE', maxloglikelihood( lambda x, loc, scale: norm2(x, loc, scale), X, p0=(0.,1.), jac=None, bounds=((None,None),(0,None)) )['x']
    print 'pn.fit2', poisson_norm.fit2(X)
    #print 'pn.fit', poisson_norm.fit(X)
    
    print 'pn.fit', poisson_norm.fit2(Y)
    exit(0)
