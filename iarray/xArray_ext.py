import numpy as np

import xarray as xr
from pint_xarray import unit_registry as ureg
from matplotlib import pyplot as plt
import uproot4 as uproot
from uncertainties import unumpy as un
from . import iprint

def genfromcsv( 
    fpath,
    x={},
    y={},
    ureg=ureg,
#     **kwargs
):
    import csv
    import re
    with open(fpath, 'r') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
    fieldnames = [ fieldname.strip("# ") for fieldname in fieldnames ]
    fieldprops = [ re.match(r"(.*)\[(.*)\]", fieldname).groups() for fieldname in fieldnames ]
    
    data = np.genfromtxt( fpath, delimiter=",")
    data = data.T
    
    xcol = x.pop("col",0)
#     coord = Coord( 
#         x.pop("name", fieldprops[xcol][0]), 
#         x.pop("dim", f"dim_{fieldprops[xcol][0]}"), 
#         data[xcol] * x.pop( "unit", ureg(fieldprops[xcol][1]) )
#     )
    ycol = y.pop("col", xcol+1)
#     array = Array( 
#         y.pop("name", fieldprops[ycol][0]), 
#         [coord], 
#         data[ycol] * y.pop( "unit", ureg(fieldprops[ycol][1]) )
#     )
#     print( data[xcol] )
    if True:
        coord = xr.DataArray(
#             name = ,
            dims = [ x.pop("dim", f"dim_{fieldprops[xcol][0]}") ],
            data = data[xcol] * x.pop( "unit", ureg(fieldprops[xcol][1]) )
        )
        array = xr.DataArray(
            name = y.pop("name", fieldprops[ycol][0]),
            dims = coord.dims,
            coords = {x.pop("name", fieldprops[xcol][0]): coord},
            data = data[ycol] * y.pop( "unit", ureg(fieldprops[ycol][1]) )
        )
    if False:
        array = xr.DataArray(
            name = y.pop("name", fieldprops[ycol][0]),
            dims = [ x.pop("dim", f"dim_{fieldprops[xcol][0]}") ],
            coords = {
                x.pop("name", fieldprops[xcol][0]): data[xcol] * x.pop( "unit", ureg(fieldprops[xcol][1]) )
            },
            data = data[ycol] * y.pop( "unit", ureg(fieldprops[ycol][1]) )
        )
    
    return array

def _genfromroot_base( fpath, treename, branches, cut ):
    """
    genfromroot( path, "hitSumm", {"E": ("oEnergy", ureg.eV)} )
    """
    aliases = None
    units = None
    if type(branches) == dict:
        units = { name: None if type(alias) == str else alias[1] for name, alias in branches.items() }
        aliases = { f"br_{idx}": alias if type(alias) == str else alias[0] for idx, alias in enumerate(branches.values()) }
        names = { f"br_{idx}": name for idx, name in enumerate(branches.keys()) }
    elif type(branches) in [tuple, list]:
        pass
    data = uproot.concatenate( 
        f"{fpath}:{treename}", 
        expressions = list(aliases.keys()),
        aliases = aliases,
        cut = cut,
        library="np"
    )
    data = xr.Dataset(
        data_vars = { names[alias]: xr.DataArray(name=names[alias], data=subdata) for alias, subdata in data.items() }
    )
    if units:
        for name, unit in units.items():
            if unit:
                data[name] = data[name] * ureg(unit)
    return data


def genfromroot( fpath, treename, branches, groupby=None, cut=None, verbose=True, tabs=0 ):
    """
    reads root files and constructs a dataset
    """
    if verbose: print( "  "*tabs + f"call genfromroot({branches})")
    _branches = branches
    if groupby:
        _branches.update( { name: name for name in groupby } )
    data = _genfromroot_base( fpath, treename, _branches, cut )
    if verbose: print( "  "*tabs + f"return genfromroot({branches})")
    return data

def groupbyArray( arr, by, func=min ):
    print("by", len(by))
    if type(by) == str:
        by = [by]
    if len(by) == 0:
        return func(arr)
    by, *rest = by
    return [
        groupbyArray( arr[mask], [ r[mask] for r in rest], func=func )
        for value in np.unique( by ) for mask in [by==iprint(value)]
    ]

# def groupby( dataset, by, func=min ):
#     if type(by) == str:
#         by = [by]
#     if by == []:
#         return func(data)
#     firstby, *restby = by
#     unique_values = np.unique( dataset[firstby] )
#     print( firstby, unique_values )
    
#     new_dataset = xr.Dataset(
#         data_vars = None,
#         coords = {firstby: unique_values},
#         attrs = {
#             value: xr.Dataset(
#                 data_vars = { name: sub[ dataset[firstby] == value ] for name, sub in dataset.items() if name != firstby },
#             )
#             for value in unique_values
#         }
#     )
#     return new_dataset

def groupby( dataset, by_list, func=min ):
    unique_lists = [ np.unique( dataset[by] ) for by in by_list ]
    print( "dims", [ len(unique_list) for unique_list in unique_lists ] )
    unique_combinations = np.array( np.meshgrid(*unique_lists) ).reshape((2, -1)).T
    print( unique_combinations )
    print( "combinations", unique_combinations.shape )
    make_selection = lambda unique_comb: np.all( 
        [ dataset[by] == unique for by, unique in zip(by_list, unique_comb ) ], 
        axis=0 
    )
    make_data = lambda name, dataarray: np.array( [ 
        func( dataarray[ make_selection( unique_comb ) ] )
        for unique_comb in unique_combinations
    ] ).reshape([len(u) for u in unique_lists[::-1]])
    data_vars = {
        iprint(name): xr.DataArray( dims = by_list[::-1], data = make_data( name, dataarray ) )
        for name, dataarray in dataset.data_vars.items() if not name in by_list
    }
    coords = { by: unique for by, unique in zip(by_list, unique_lists) }
    return xr.Dataset( data_vars = data_vars, coords = coords )

    
def Coord( name, dims, data ):
    a = xr.DataArray( name = name, dims = dims, data = data )
    if not isinstance(data, xr.DataArray):
        a[name] = a
    return a

def Array( name, coords=None, data=None ):
    if coords is None:
        data.name = name
        return data
    if isinstance(coords, xr.DataArray):
        raise Exception("trying to create new DataArray from existing")
    return xr.DataArray(
        name = name,
        dims = [ dim for coord in coords for dim in coord.dims ],
        coords = { coord.name: coord for coord in coords },
        data = data.variable if hasattr(data, "variable") else data,
    )

xr.DataArray.xcoords = property(
    lambda self: [ self[coord] for coord in self.coords ]
)

def xplot(self, *args, **kwargs):
    dq = self.pint.dequantify()
    unit = kwargs.pop("unit", None)
    if unit:
        dq = self.pint.to(unit).pint.dequantify()
    x = kwargs.pop("x", None)
    y = kwargs.pop("y", None)
    for idx, coord in zip(["x","y"], dq.coords):
        if not idx in kwargs:
            kwargs[idx] = coord

    ylabel = kwargs.pop("ylabel", None)
    if ylabel:
        dq.attrs["long_name"] = ylabel
        kwargs["label"] = dq.name
    if not "x" in kwargs:
        x = dq.coords.values()[0]
    try:
        dq.plot(*args, **kwargs)
    except ValueError:
        y = kwargs.pop("y", None)
        dq.plot(*args, **kwargs)
    if "label" in kwargs:
        plt.legend()
    return
xr.DataArray.xplot = xplot



def xdiff(self, coord=None):
    if coord is None:
        coord = self.xcoords[0].name
#     self[coord].name = ""
#     print( self[coord] )
    dequantified = self.pint.dequantify()
    diff_dequantify = dequantified.differentiate( coord )
    return diff_dequantify.pint.quantify()/self[coord].pint.units * ( self.pint.units or 1 )
xr.DataArray.xdiff = xdiff



def ds_stats(self, name):
    mean = self.mean(name)
    std = self.std(name)
    return xr.Dataset( 
        { 
            da_name: (
                mean[da_name].dims,
                un.uarray( mean[da_name].variable, std[da_name].variable ) 
            ) for da_name in mean.data_vars.keys()       
        },
        coords = mean.coords
    )
xr.Dataset.stats = ds_stats



def ds_errorbar(self, cols_wrap=2, **kwargs):
#     print( "ds_errorbar" )
    z = kwargs.pop("z", None)
    if z:
        size = self[z].size
        nrows, ncols = int(np.ceil(size/cols_wrap)), cols_wrap
        width, height = plt.rcParams["figure.figsize"]
        fig = plt.figure( figsize = (width*ncols, height*nrows) )        
        for i, zi in enumerate(self[z]):
            ax = fig.add_subplot(nrows, ncols, i+1)
            ax.grid(True)
            self.sel({z:zi}).xerrorbar(ax=ax, label=f"{z} {zi.data}", **kwargs)
        fig.tight_layout()
        return
    y = kwargs.pop("y", None)
    x = kwargs.pop("x", None)
    ylabel = kwargs.pop("ylabel", None)
    yunits = kwargs.pop("yunits", None)
    yscale = kwargs.pop("yscale", None)
    ax = kwargs.pop("ax", None)
    if not ax:
        ax = plt.figure().subplots()
    ax.set_xlabel( f"{x}" + (f"[{self[x].pint.units}]" if self[x].pint.units != '' else '') )
    
    if ylabel:
        ax.set_ylabel( f"{ylabel}" + (f"[{yunits}]" if yunits else '') )
    else:
        ax.set_ylabel( f"{y}" + (f"[{self[y].pint.units}]" if self[y].pint.units != '' else '') )

    ds = self.pint.dequantify()
    label = kwargs.pop("label", None)

#     if yunits:
#         ds = self.pint.to(yunits).pint.magnitude
#         ds = ds.pint.dequantify()
#         print( "yunits" )
    for iy in ([y] if type(y) == str else y):
        _x = un.nominal_values(ds[x].data)
        _dx = (_x[1]-_x[0])/2
        args = dict(
            x = _x,
            xerr = _dx,
            y = un.nominal_values(ds[iy].data),
            yerr = un.std_devs(ds[iy].data),
            label = iy
        )
        ax.errorbar(
#             label = ds[iy].name,
            **args,
            **kwargs
        )
    if yscale:
        ax.set_yscale( yscale )
    plt.title(label)
    if label:
        plt.legend()
xr.Dataset.xerrorbar = ds_errorbar



def ds_histogram( self, by, bins=None, centers=None, per_bin=True ):
    if not bins:
        bins = (centers[1:] + centers[:-1])/2
        centers = (bins[1:] + bins[:-1])/2
        binsize = (bins[1] - bins[0]).data
    unique = np.unique(self[by])
    unique_C = Coord("hdu", "hdu", unique)
    def make_hist( name, value ):
        return np.histogram( 
            self[name][ self[by] == value ].pint.to(bins.pint.units).pint.magnitude, 
            bins.pint.magnitude 
        )[0].astype(float)
    
    coords = {by: unique, centers.name: centers}
    data_vars = {}
    for name in self.data_vars.keys():
        if name == by: continue
        variable = []
        for value in unique:
            val = make_hist(name, value)
            variable.append( val )
        variable = np.array( variable )
        da = xr.DataArray( 
            data = un.uarray(variable, np.sqrt(variable) )* ureg("counts")/( binsize if per_bin else 1),
            coords = coords,
            dims = [ by, *centers.dims ]
        )
        data_vars[ fr"$dN/d{centers.name.replace('$', '')}$" ] = da

    return xr.Dataset( data_vars = data_vars, coords = coords )
xr.Dataset.xhistogram = ds_histogram



def DataArray_errorbar(self, x, z, cols_wrap=1, fmt=" ", **kwargs):
    print( "da_errorbar" )
    size = self[z].size
    nrows, ncols = int(np.ceil(size/cols_wrap)), cols_wrap
    width, height = plt.rcParams["figure.figsize"]
    fig = plt.figure( figsize = (width*ncols, height*nrows) )
    yunits = kwargs.pop("yunits", self.pint.units)
    xlim = kwargs.pop("xlim", None)

    for i, sel in enumerate(self[z]):
        ax = fig.add_subplot(nrows, ncols, i+1)
        try:
            y = self.sel( {z: sel} ).pint.to(yunits).pint.magnitude
        except Exception as e:
            print( f"trying to select \n{self}\n\n in coord \n{z}\n\n and selection \n{sel}" )
            raise e
        try:
            ax.errorbar(
                self[x].pint.magnitude, 
                un.nominal_values( y ), 
                xerr = (self[x].pint.magnitude[1] - self[x].pint.magnitude[0])/2,
                yerr = un.std_devs( y ), 
                fmt=fmt, 
                label=f"hdu {sel.data}",
                **kwargs,
            )
            ax.grid(True)
        except Exception as e:
            print( f"error in\n{ self[x].pint.magnitude }\n\n{y} ")
            raise e
        if "fit" in self.attrs:
            popt = self.attrs["fit"].sel({z:sel}).pint.magnitude.data
            func = self.attrs["fit"].attrs["func"]
            if type(func) == str:
                _func = lambda x, *p: eval(func)
            else:
                _func = func
            ax.plot( 
                self[x].pint.magnitude,
                _func( self[x].pint.magnitude, *popt ),
                label = "\n".join([ f"$p_{i}$={p:.2e}" for i, p in enumerate(popt) ])
            )
        ax.set_xlim(xlim)
        ax.set_xlabel( f"{self[x].name} [{self[x].pint.units}]" )
        if str(yunits) == "":
            ax.set_ylabel( f"{self.name}" )
        else:
            ax.set_ylabel( f"{self.name} [{yunits}]" )
        ax.legend()
    plt.tight_layout()
    return self
xr.DataArray.errorbar = DataArray_errorbar

def print_table(da, x, y=None, out=None):
    X = da[x]
    if y:
        Y = da[y]
    else:
        Y = da
        y = da.name
    
    s = [f"# {x}[{X.pint.units}], {y}[{Y.pint.units}]"] + [
        f"{ix:.2f}, {iy:.2f}" for ix, iy in zip( X.data.magnitude, Y.data.magnitude )
    ]
    if out:
        with open(out, "w") as f:
            f.write("\n".join(s))
    else:
        print( "\n".join(s) )
    return da