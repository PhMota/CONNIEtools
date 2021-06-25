import numpy as np
from numpy import sqrt, pi, arange, unique, all, any
import xarray as xr
from xarray import DataArray
import pint
import pint_xarray
from pint_xarray import unit_registry as ureg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic

from uncertainties import unumpy as un

import warnings
warnings.filterwarnings("error")

plt.rcParams.update({
    "image.origin": 'lower',
    "font.family": "serif",
    "font.size": 15,
    "grid.alpha": .5,
    "text.usetex": True,
})

def iprint(*args):
    print(*args)
    return args[-1]

ureg.default_format = '~P'
ureg.auto_reduce_dimensions = True

amu = 9.31e8*ureg.eV
eV = ureg.eV
keV = ureg.keV
MeV = ureg.MeV
GeV = ureg.GeV

from .xArray_ext import *
###############################################################################
def DataArray_curve_fit(self, x, z, func, p0, **kwargs):
    popt_list = []
    if type(func) == str:
        _func = lambda x, *p: eval(func)
    else:
        _func = func
    for sel in self[z]:
        y = self.sel( {z: sel} ).pint.magnitude
        sigma = un.std_devs(y)
        popt, _ = curve_fit( 
            _func, 
            self[x].pint.magnitude, 
            un.nominal_values(y), 
            sigma = np.where( sigma>0, sigma, np.ones_like(sigma) ),
            p0=p0,
            **kwargs
        )
        popt_list.append(popt)
    popt = xr.DataArray( 
        popt_list, 
        dims=[z, "params"],
        coords = {
            z: self[z]
        },
        attrs = {
            "func": func,
            "fit": popt
        }
    )
    self.attrs["fit"] = popt
    return self
xr.DataArray.curve_fit = DataArray_curve_fit
###############################################################################
class Selection:
    def __init__(cls, obj):
        if isinstace(obj, cls):
            return obj
        return super().__new__(cls)

    def __init__(self, s=None):
        self.selection = s
    
    def __repr__(self):
        return self.selection

    def str(self):
        return self.selection

    @classmethod
    def distribute(cls, key, values, op="=="):
        return [ f"{key}{op}{value}" for value in values ]

    
    @classmethod
    def __op__(cls, *args, op=" & "):
        return Selection(
            op.join( [ f"({arg})" for arg in args ] )
        )

    @classmethod
    def prod(cls, *args, key=None, op="=="):
        if key:
            return cls.prod( *cls.distribute(key, args[0], op) )
        return cls.__op__(*args)

    @classmethod
    def sum(cls, *args, key=None, op="=="):
        if key:
            return cls.sum( *cls.distribute(key, args[0], op) )
        return cls.__op__(*args, op=" | ")

    def __mul__(self, other):
        return self.__op__(self, other)

    def __add__(self, other):
        return self.__op__(self, other, op=" | ")
    

###############################################################################
def sel(a, cond):
    return a[cond]

def sels(a, *conds, cmp=any):
    return sel(a, cmp( conds, axis=0 ) )

def cmp_list( a, values, cmp=np.ndarray.__eq__ ):
    return [ cmp(a, value) for value in values ]

scatterByGroups = lambda a, x, y, selections: (
    plt.scatter( 
        groupby(a[x], *[ a[s] for s in selections ] )[0],
        groupby(a[y], *[ a[s] for s in selections ] )[0]
    ), 
    plt.xlabel(x),
    plt.ylabel(y)
)

scatterExcluded = lambda a, x, y, s, excluded: (
    plt.scatter(
        sels(a[x], *cmp_list(a[s], excluded) ),
        sels(a[y], *cmp_list(a[s], excluded) ),
        marker = "x"
    ),
    plt.xlabel(x),
    plt.ylabel(y)
)

scatterByGroupsExcluded = lambda a, x, y, excluded: (
    scatterByGroups( a, x, y, ["runID", "ohdu"] ),
    scatterExcluded( a, x, y, "runID", excluded )
)