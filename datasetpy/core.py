from .defs import ureg, plt
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

def dict2str(d):
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
                
        self._dims = aslist(dims)
        self.attrs = kwargs
    
    @property
    def dims(self):
        return { dim: size for dim, size in zip(self._dims, self._data.shape) }
    
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
    def label(self): 
        _label = self.name
        if hasattr(self, "latex"):
            _label = self.latex
        return f"{_label}" + ( "" if self.units == "dimensionless" else f" [{self.units:~P}]" )
    
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

    
def parse_array( data, err, units ):
    ret = un.uarray( data, err ) if np.any( err != 0 ) else data
    return ret * ( ureg("dimensionless") if units is None else ureg(units) )

class Array:
    """
    describes 
    F(x, y) -> z
    """
    def __init__(self, name, data, dims = None, coords = None, err = 0, units = "dimensionless", **kwargs):
        self.name = name
        self._data = parse_array(data, err, units)
        _dims = {}
        if coords:
            self._coords = aslist(coords)
            self.__dict__.update( { coord.name: coord for coord in self._coords} )
            for coord in self._coords:
                for dim_name, dim_size in coord.dims.items():
                    if dim_name in _dims:
                        if dim_size != _dims[dim_name]:
                            raise Exception(f"mismatching dimension sizes for {dim_name}: {_dims[dim_name]} and {dim_size}")
                    else:
                        _dims[dim_name] = dim_size
        
        
        self.attrs = kwargs
        if dims is None:
            self._dims = list({ dim for coord in self.coords for dim in coords._dims })

    @property
    def dims(self):
        return { dim: size for dim, size in zip(self._dims, self._data.shape) }

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
    def label(self): 
        _label = self.name
        if hasattr(self, "latex"):
            _label = self.latex
        return f"{_label}" + ( "" if self.units in ["dimensionless", ""] else f" [{self.units:~P}]" )

    def print_header(self):
        return [f"<{self.__class__.__name__}>\t({print_dict(self.dims)})"]
    
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
            *[ coord.print_data()[0] for coord in self.coords ]
        ])]

    def __repr__(self):
        return "\n".join(
            self.print_header()
            + self.print_data()
            + self.print_coords()
            + self.print_attrs()
        ) + "\n"
        
        
        
class Data:
    """
    mimics xarray dataset interface
    """
    def __init__(self, name, arrays, coords = None, **kwargs):
        self.name = name
        arrays = aslist(arrays)
        _dims = {}
        for array in arrays:
            for dim_name, dim_size in array.dims.items():
                if dim_name in _dims:
                    if dim_size != _dims[dim_name]:
                        raise Exception(f"mismatching dimension sizes for {dim_name}: {_dims[dim_name]} and {dim_size}")
                else:
                    _dims[dim_name] = dim_size
        print( _dims )

        self.coords = {}
        self.dims = {}
        self.attrs = {}
    
    @property
    def dims(self):
        return
    
    @property
    def header(self):
        return f"<{self.__class__.__name__}>\t({dict2str(self.dims)})"
    
    def __repr__(self):
        return "\n".join([
            self.header,
            
        ]) + "\n"
    
    def __setitem__(self, attr, value):
        pass
    
    def __getitem__(self, attr):
        if isinstance(attr, str):
            return self._getitem_str(attr)
        
