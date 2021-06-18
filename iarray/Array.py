import xarray as xr
import numpy as np
import csv
import re
from pint import UnitRegistry
from matplotlib import pyplot as plt

from .Coord import Coord

def genfromcsv( 
    fpath, 
    x=None,
    y=None,
    ureg=UnitRegistry() 
):
    with open(fpath, 'r') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
    fieldnames = [ fieldname.strip("# ") for fieldname in fieldnames ]
    fieldprops = [ re.match(r"(.*)\[(.*)\]", fieldname).groups() for fieldname in fieldnames ]
    
    data = np.genfromtxt( fpath, delimiter=", ")
    data = data.T
    
    xcol = x.pop("col",0)
    coord = Coord( 
        x.pop("name", fieldprops[xcol][0]), 
        x.pop("dim", f"dim_{fieldprops[xcol][0]}"), 
        data[xcol] * x.pop("unit", ureg(fieldprops[xcol][1]) )
    )
    ycol = y.pop("col", xcol+1)
    array = Array( 
        y.pop("name", fieldprops[ycol][0]), 
        [coord], 
        data[ ycol ] * y.pop("unit", ureg(fieldprops[xcol][1]) )
    )
    return array

class Array:
    def __init__( self, name, coords, data ):
        """
        constructure with signature
        name: str
        coords: list of icoord
        data: ndarray-like
        """
        assert hasattr( coords, "__iter__" ), \
            f"iarray must be initiated with a sequence of coordinates and not type(coords) = {type(coords)}"
        self._da = xr.DataArray(
            attrs = dict(
                name = name,
                long_name = name,
            ),
            dims = [ coord.dim for coord in coords ],
            data = data,
        )
        for coord in coords:
            self._da[coord.name] = coord.asDataArray()

    def __repr__(self):
        """
        return xarray representation
        """        
        return self._da.__repr__() + "\n"
    
    @property
    def name(self):
        return self._da.attrs["name"]
    
    @property
    def coords(self):
        return self._da.coords
    
    @property
    def coord_names(self):
        return list(self.coords.keys())
    
    def plot(self, yunits=None, xunits=None, **kwargs):
        da = self._da.assign_attrs( long_name = kwargs.pop("ylabel", self.name) )
        if yunits:
            da.pint.to(yunits)
        coord = kwargs.pop( "x", self.coord_names[0] )
        da.pint.dequantify().plot( x=coord, label = self.name, **kwargs )    
        if da.attrs["long_name"] != self.name:
            plt.legend()
        