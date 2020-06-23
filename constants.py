
import numpy as np

energies = {
    'e-': "3.745eV",
    'Si': "1.740keV",
    'Cu': "8.046keV",
    'Cu2': "8.904keV"
    }

ccd_shape = np.array([4130,4120])

diffusion_function = 'sqrt(-258.817238*log1p(-0.000982*z))/15 if z < 670 else 0'
charge_efficiency_function = '1. if z < 670 else .9'
vertical_modulation_function = '50*cos(y/2/pi/20)'
horizontal_modulation_function = '-1e3*(x/1000. - 1)**2 if x < 1000. else 0'

gains = {
    2: 5.71617,
    3: 7.16026,
    4: 6.13546,
    5: 6.38599,
    6: 6.12512,
    7: 5.59342,
    8: 7.00279,
    9: 7.31909,
    10: 6.74754,
    13: 8.14274,
    14: 7.23029,
    15: 8.1933,
    }
