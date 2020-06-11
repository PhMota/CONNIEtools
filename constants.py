

energies = {
    'e-': "3.745eV",
    'Si': "1.740keV",
    'Cu': "8.046keV",
    'Cu2': "8.904keV"
    }

ccd_shape = [4130, 4120]

diffusion_function = 'sqrt(-258.817238*log1p(-0.000982*z))/15 if z < 670 else 0'
charge_efficiency_function = '1. if z < 670 else .9'
vertical_modulation_function = '50*cos(y/2/pi/20)'
horizontal_modulation_function = '-1e3*(x/1000. - 1)**2 if x < 1000. else 0'
