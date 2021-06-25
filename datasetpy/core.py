from .defs import *
import numpy as np
from uncertainties import unumpy as un
from uncertainties.core import AffineScalarFunc

par = " "*4
significant_digits = 1
def fmt(a):
    if not isinstance(a, AffineScalarFunc):
        return f".6"
    return f".{significant_digits}uP"

np.set_printoptions(
    threshold = 10,
    formatter = {
        "all": lambda x: f"{x:{fmt(x)}}"
    }
)

def pprint(arr):
    u = f"[{arr.units:~P}]\n"
    return u + str( arr.magnitude )
#     if arr.size > 8:
#         return u + ", ".join( f"{a:{fmt(a)}}" for a in arr[:4].magnitude ) + " ... " + ", ".join( f"{a:{fmt(a)}}" for a in arr[-4:].magnitude )
#     return u + ", ".join( f"{a:{fmt(a)}}" for a in arr.magnitude )

def is_str(arg):
    return isinstance(arg, str)

def is_list_or_tuple(arg):
    return isinstance(arg, list) or isinstance(arg, tuple)

def aslist(arg):
    if not type(arg) in [list, tuple]:
        return [arg]
    return arg

def has_std_devs(arg):
    if not arg.dtype.name == "object":
        return False
    if not hasattr( arg[0], "std_dev" ):
        return False
    if not np.any( un.std_devs(arg) != 0 ):
        return False
    return True

def has_units(arg):
    return hasattr(arg, "units")

def print_dict(d):
    return ", ".join( f"{key}: {value}" for key, value in d.items() )

def print_values(X):
    return f"{par}{X.name}\t({', '.join( X.dims.keys() )}) {pprint( X._data )}"




class Coord:
    def __init__(self, name, data, dims, err = 0, units = "dimensionless", **kwargs):
        self.name = name
        self._data = data
        if not hasattr(self._data[0], "std_dev"):
            self._data = un.uarray( self._data, err )
        if not np.any( self.std_devs != 0 ):
            self._data = self.nominal_values
        if not hasattr(self._data, "units"):
            self._data *= ureg(units)
                
        self.dims = { dim: size for dim, size in zip(aslist(dims), self.shape) }
        self.attrs = kwargs
        
    @property
    def nominal_values(self): return un.nominal_values( self._data )

    @property
    def std_devs(self): return un.std_devs( self._data )

    @property
    def units(self): return self._data.units

    @property
    def n(self): return self.nominal_values

    @property
    def s(self): return self.std_devs

    @property
    def u(self): return self.units

    @property
    def shape(self): return self._data.shape
    
    @property
    def label(self): return f"{self.name}" + ( "" if self.units == "dimensionless" else f" [{self.units:~P}]" )
    
    def print_header(self):
        return [f"<{ self.__class__.__name__ }>\t({ print_dict(self.dims) })"]

    def print_data(self):
        return [ f"{par}{self.name}\t({', '.join(self.dims.keys())}) {pprint(self._data)}" ]
    
    def print_attrs(self):
        if not self.attrs == {}:
            return []
        return ["\n".join([ 
            "Attributes:",
            *[ f"{par}{key}:\t{value}" for key, value in self.attrs.items() ]
        ])]
    
    def __repr__(self):
        return "\n".join(
            self.print_header()
            + self.print_data()
            + self.print_attrs()
        ) + "\n"

    def __eq__(self, other):
        if isinstance(other, Coord):
            return self._data == other._data
        return self._data == other
    
    def __add__(self, other):
        print( isinstance(other, Coord) )
        if isinstance(other, Coord):
            return self._data.__add__( other._data )
        return self._data.__add__( other )

    

class Array:
    """
    describes 
    F(x, y) -> z
    """
    def __init__(self, name, data, dims = None, coords = None, err = 0, units = None, **kwargs):
        self.name = name
        self._data = un.uarray( data, err ) * ureg(units)
        self.coords = { coord.name: coord for coord in aslist(coords) }
        self.attrs = kwargs
        if dims is None:
            self.dims = { dim_name: dim_size for coord in aslist(coords) for dim_name, dim_size in coords.dims.items() }
        
    @property
    def nominal_values(self): return un.nominal_values( self._data )

    @property
    def std_devs(self): return un.std_devs( self._data )

    @property
    def units(self): return self._data.units

    @property
    def n(self): return self.nominal_values

    @property
    def s(self): return self.std_devs

    @property
    def u(self): return self.units

    @property
    def label(self): return f"{self.name}" + ( "" if self.units in ["dimensionless", ""] else f" [{self.units:~P}]" )

    def print_header(self):
        return [f"<{self.__class__.__name__}>\t()"]
    
    def print_data(self):
        return [f"{par}{self.name}\t({', '.join(self.dims.keys())}) {pprint(self._data)}"]
    
    def print_attrs(self):
        if not self.attrs == {}:
            return []
        return ["\n".join([ 
            "Attributes:",
            *[ f"{par}{key}:\t{value}" for key, value in self.attrs.items() ]
        ])]

    def print_coords(self):
        if self.coords is None or len(self.coords) == 0:
            return []
        return ["\n".join([ 
            "Coordinates:",
            *[ coord.print_data()[0] for coord in self.coords.values() ]
        ])]

    def __repr__(self):
        return "\n".join(
            self.print_header()
            + self.print_data()
            + self.print_coords()
            + self.print_attrs()
        ) + "\n"
        
        
    def plot(self):
        ax = plt.figure().add_subplot(111)
        for coord in self.coords.values():
            ax.errorbar( 
                x = coord.n,
                xerr = coord.s,
                y = self.n,
                yerr = self.s,
                fmt = "."
            )
            ax.set_xlabel( coord.label )
            ax.set_ylabel( self.label )
        
    
class Data:
    """
    mimics xarray dataset interface
    """
    def __init__(self, **kwargs):
        self.name = name
        
        self.data = {}
        self.coords = {}
        self.dims = {}
        self.attrs = {}
    
    def __repr__(self):
        s = [f"<datasetpy.Dataset>"]
        s += [ "Dimensions:\t\t" + repr(self.dims) ]
        c = [ "Coordinates:" ]
        return 
    
    def __setitem__(self, attr, value):
        pass
    
    def __getitem__(self, attr):
        if isinstance(attr, str):
            return self._getitem_str(attr)
        
    def plot(self, x=None, y=None, **kwargs):
        """
        interface for lineplot
        """
        x = self._getx(x)
        y = self._gety(y)
        
        