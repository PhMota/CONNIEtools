from pint import UnitRegistry, UnitStrippedWarning

ureg = UnitRegistry()
ureg.default_format = '~P'
ureg.auto_reduce_dimensions = True

import warnings
warnings.filterwarnings("error", UnitStrippedWarning)

import matplotlib.pyplot as plt
plt.rcParams.update({
    "image.origin": "lower",
    "font.family": "serif",
    "font.size": 15,
    "grid.alpha": .5,
    "axes.grid": True,
    "text.usetex": True,
    "xaxis.labellocation": "right",
    "yaxis.labellocation": "top",
})
warnings.filterwarnings("ignore", DeprecationWarning)