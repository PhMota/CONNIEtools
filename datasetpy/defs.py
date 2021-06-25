from pint import UnitRegistry

ureg = UnitRegistry()
ureg.default_format = '~P'
ureg.auto_reduce_dimensions = True

import matplotlib
matplotlib.rcParams.update({
    "image.origin": "lower",
    "font.family": "serif",
    "font.size": 15,
    "grid.alpha": .5,
    "axes.grid": True,
    "text.usetex": True,
    "xaxis.labellocation": "right",
    "yaxis.labellocation": "top",
})

import matplotlib.pyplot as plt