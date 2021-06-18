import xarray as xr
import pint_xarray 

def crange(name, dim, xmin, xmax, dx, unit ):
    return Coord( name, dim, np.arange(xmin, xmax, dx) * unit )

def asCoord( a ):
    return Coord( None, a.dims[0], a.data )
    
class Coord:
    def __init__( self, name, dim, data ):
        """
        constructor with signature
        name: str
        dim: str
        data: ndarray_like(1d)
        """
        assert isinstance(dim, str), f"icoord must be initiated with a str dim, but type(dim) == {type(dim)}"
        self._da = xr.DataArray(
            attrs = dict(
                name = name,
                math_name = name
            ),
            dims = dim,
            data = data,
        )
        
    def asDataArray(self):
        """
        returns the DataArray
        """
        return self._da
    def __repr__( self ):
        """
        DataArray representation
        """
        return self._da.__repr__()
    
    @property
    def data(self):
        """
        get data
        """        
        return self._da.data
    @property
    def dim(self): 
        """
        get (single) dimension
        """        
        return self._da.dims[0]
    @property
    def name(self):
        """
        get name
        """        
        return self._da.attrs["name"]
    
    def dequantify(self):
        return self._da.pint.dequantify()
    
    def quantify(self):
        return self._da.pint.quantify()

    @property
    def unit(self):
        return self.data.unit if hasattr(self.data, "unit") else 1
    
    def __pow__(self, a):
        a = float(a)
        return asCoord( self._da**a )
    
    def __mul__(self, a):
        da = (self.dequantify()*a.dequantify()).pint.quantify()
        return Coord( None, da.dims, da.data )
    